!_______________________________________________________________________________________
!_________________________________________________________________________mohid.Life.1.0



!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Surface Water
! MODULE        : Life
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Marcos Mateus - v1.0
! DESCRIPTION   : It all started in the water...
!                 Module containing a multi-element biogeochemical water-column model
!                 with variable stoichiometry 
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

Module ModuleLife

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleFunctions, only: PhytoLightLimitationFactor, Chunk_I

    implicit none

    private 


!_____________________________________________________________________________________
!   -§-     Defaul value not assigned. Definition required
!   N.D.    Non dimensional
!
!  ## ATTENTION: BOTH VALUES AND UNITS ARE SPECIFIC FOR THE FOLLOWING CONDITIONS:  ##
!  ##                C & Chl : mg / m3                                             ##
!  ##               N, P, Si : mmol / m3                                           ##
!  ##                O2, CO2 : mg / L                                              ##
!_____________________________________________________________________________________

!____________________________________________________________________
!__________________________________________________________GLOBAL____

!---keyword-----------------units----------------value--------------description
!
!   REDFIELD_NC         :   mmol N / mg C         [0.01261]          !Redfield N:C ratio
!   MIN_NC_RATIO        :   mmol N / mg C         [ -§- ]            !Minimal N:C ratio
!   MAX_NC_RATIO        :   mmol N / mg C         [ -§- ]            !Maximal N:C ratio
!   REDFIELD_PC         :   mmol P / mg C         [0.000786]         !Redfield P:C ratio
!   MIN_PC_RATIO        :   mmol P / mg C         [ -§- ]            !Minimal P:C ratio
!   MAX_PC_RATIO        :   mmol P / mg C         [ -§- ]            !Maximal P:C ratio
!   REDFIELD_SiC        :   mmol Si / mg C        [0.03]             !Standard Si:C ratio
!   Q10_VALUE           :   N.D.                  [ -§- ]            !Q10 value for temperature limitation
!   BIO_SI_DISS         :   /day                  [0.01]             !Biogenic silica dissolution rate
!   O2_CARB_CONVERS     :   mg O2 L-1 / mg C m-3  [0.0026666]        !Oxygen to carbon conversion factor
!   CO2_CARB_CONVERS    :   mg CO2 L-1 / mg C m-3 [0.0036641]        !Carbon Dioxide to carbon conversion factor
!   POM_BAC_KS          :   mg C m-3              [200]              !Bacteria mediated POM Hydrolysis MM cosntant
!   POM_BAC_VMAX        :   / day                 [2]                !Vmax for POM Hydrolysis
!   DOMSL_BAC_KS        :   mg C m-3              [400]              !Bacteria mediated DOMsl Hydrolysis MM cosntant
!   DOMSL_BAC_VMAX      :   / day                 [2]                !Vmax for DOMsl Hydrolysis
!   LIGHT_LIM_METHOD    :   ----                  [1]                !1:Iopt variable; 2:Waterquality; 3:Schematic
!   TEMP_LIM_METHOD     :   ----                  [2]                !1:Q10; 2:Moore et al.(2002)
!   MASS_XEK            :   BOLEAN                [0]                !Mass conservation teste activation
!   REF_TEMP            :   C                     [30]               !Reference temperature
!   REF_TEMP_Q10        :   ºC                    [20]               !Reference temperature for Q10 method
!   NITRIFRATE          :   d-1                   [0.04]             !Nitrification rate
!   NITRIFRADLIM        :   W/m2                  [4.0]              !Light radiation bellow which nitrification occurs
!   NIT_IN_COEF         :   L mg-1                [0.6]              !Nitrification inhibition coeficient
!   NIT_O_N_CONV        :   mg L-1 / mmol m-3     [0.0588]           !Nitrification O:N consumption ratio

!______________________________________________________________________
!_________________________________________________________PRODUCERS____
!
!---keyword-----------------units----------------value--------------description
!
!   MAX_ASSIMIL         :   /day                 [ -§- ]            !Maximal assimilation rate
!   EXU_NUT_STRESS      :   N.D.                 [ -§- ]            !Exudation under nutrient stress
!   RESP_BASAL          :   /day                 [ -§- ]            !Basal respiration rate
!   RESP_FRAC_PROD      :   N.D.                 [ -§- ]            !Respired fraction of production
!   MIN_LYSIS           :   /day                 [0.05]             !Minimal lysis rate
!   SED_NUT_STRESS      :   m/day                [ -§- ]            !Nutrient stress sedimentation rate
!   SED_MIN             :   m/day                [0]                !Minimal sedimentation rate
!   NUT_STRESS_THRESH   :   N.D.                 [ -§- ]            !Nutrient stress threshold (sedimentation)
!   MAX_STORE_FILL      :   /day                 [1]                !Maximal rate of storage filling
!   AFFINITY_PO4        :   (mg C)-1 m-3 d-1     [0.0025]           !Affinity for PO4 uptake
!   AFFINITY_NH4        :   (mg C)-1 m-3 d-1     [0.0025]           !Affinity for NH4 uptake 
!   AFFINITY_NO3        :   (mg C)-1 m-3 d-1     [0.0025]           !Affinity for NO3 uptake
!   REL_EXCESS_SI       :   /day                 [1]                !Release rate of excess silicate
!   SI_UPTAKE_KS        :   mmol Si m-3          [0.3]              !Silicate uptake Michaelis constant
!   SILICA_USE          :   BOLEAN               [ -§- ]            !Set Silica use by the producer
!   EXC_DOM_SL_FRAC     :   N.D.                 [0.1]              !DOM diverted to semi-labile pool
!   PHOTOINHIBITION     :   W/m2                 [100]              !Photoinhibition 
!   MAX_CHLN_RATIO      :   mg Chl / mmol N      [3.0]              !Maximal Chl:N ratio
!   ALPHA_CHL           :   mg C m2 / mg Chl W d [0.25*12,01]       !Chl specific initial slop of P vs I curve 
!   CHL_DEGRAD_RATE     :   d-1                  [0.0]              !Chl degradation rate constant
!   MIXOTROPHY          :   BOOLEAN              [0]                !Hability to performe mixotrophy

!______________________________________________________________________
!_________________________________________________________CONSUMERS____
!
!---keyword-----------------units----------------value--------------description
!
!   REST_RESP_@10C      :   / day                [0.02]             !rest respiration @ 10ºC
!   ASSIMIL_EFFIC       :   N.D.                 [ -§- ]            !assimilation efficiency
!   EXCRE_UP_FRAC       :   N.D.                 [0.5]              !excreted fraction of uptake
!   MORT_POM_FRAC       :   N.D.                 [0.5]              !fraction of mortality to POM
!   MORT_O2_DEP         :   /day                 [0.25]             !oxygen-dependent mortality rate
!   MORT_RATE           :   /day                 [0.05]             !temperature-independent mortality rate
!   O2_KS               :   mg O2 / L            [0.25]             !oxygen half saturation constant 
!   MAX_SPEC_UP_@10C    :   / day                [ -§- ]            !maximum specific uptake @ 10ºC
!   GRAZ_UP_KS          :   mg C / m3            [ -§- ]            !half saturation value for uptake    
!   GRAZ_AVAIL          :   N.D.                 [ -§- ]            !availability of Prey X
!   SEDIM_MIN           :   m /day               [ -§- ]            !manimal sedimentation rate
!   SEDIM_NUT_STRESS    :   m /day               [ -§- ]            !nutrient stress sedimentation rate!
!________________________________________________________________________
!_________________________________________________________DECOMPOSERS____
!
!---keyword-----------------units----------------value--------------description
!
!   LYS_REF_CON         :   mg C / m3            [50]               ! ###Lysis_Ref_Con    
!   MORT_DOM_SL_FRAC    :   N.D.                 [0.2]              ! ###DOC_SL_Frac    
!   MORT_POM_FRAC       :   N.D.                 [0.4]              !fraction of mortality to POM
!   O2_KS               :   mg O2 / L            [0.01]             !oxygen half saturation constant
!   O2_LOW_ASS_EFIC     :   mg O2 / L            [1.6]              ![oxygen] bollow which ass efic is low                   
!   MORT_RATE           :   /day                 [0]                !temperature-independent mortality rate
!   DENS_DEP_MORT       :   /day                 [0.5]              !density-dependence mortality rate
!   REST_RESP_@10C      :   / day                [0.01]             !rest respiration @ 10ºC
!   EXCRE_UP_FRAC       :   N.D.                 [0.1]              !excreted fraction of uptake
!   ASS_EFFIC           :   N.D.                 [0.3]              !Asimilation efficiency
!   ASS_EFFIC_LOW_O2    :   N.D.                 [0.3]              !Asimilation efficiency @ low O2
!   MAX_SPEC_UP_@10C    :   / day                [8.28]             !maximum specific uptake @ 10ºC
!   PO4_Ks              :   (mg C)-1 m-3 d-1     [0.0025]           !PO4 uptake affinity     
!   NH4_Ks              :   (mg C)-1 m-3 d-1     [0.0025]           !PO4 uptake affinity     
!   NO3_Ks              :   (mg C)-1 m-3 d-1     [0.0025]           !PO4 uptake affinity
!   DOM_UP_KS           :   mg C / m3            [ -§- ]            !half saturation value for DOM uptake     



!________________________________________________________________________
!___________________________________________________AUTHOR'S COMMENTS____
!
! Mixotrophy: needs update with CO2 dynamics
!





    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructLife
    private ::      AllocateInstance
    private ::      ReadData
    private ::          ConstructGlobalVariables
    private ::          ConstructModelOptions
    private ::          ConstructProducers
    private ::              AddProducer
    private ::              ConstructProducerParameters
    private ::                  ConstructProducerLimitation
    private ::                      ConstructLightLimitation
    private ::                  ConstructProducerUptake
    private ::                  ConstructProducerExudation
    private ::          ConstructConsumers
    private ::              AddConsumer
    private ::              ConstructConsumerParameters
    private ::          ConstructDecomposers
    private ::              AddDecomposer
    private ::              ConstructDecomposerParameters
    private ::                  ConstructDecomposerConsume
    private ::          PropertyIndexNumber
    private ::          ConstructPropertyList

    private ::              ConstructExcretion
    private ::              ConstructRespiration
    private ::              ConstructMortality
    private ::              ConstructRatio
    private ::              ConstructGrazing
    private ::                  AddPrey
    private ::                  ConstructPrey
    private ::              ConstructTemperatureLimitation

    
    !Selector
    public  :: GetDTLife
    public  :: GetLifePropertyList
    public  :: GetLifeSize
    public  :: GetLifePropIndex
    public  :: UnGetLife
                     
    
    !Modifier
    public  :: ModifyLife
    private ::      Producers
    private ::      Consumers
    private ::      Decomposers
    private ::      BioChemicalProcesses




    !Destructor
    public  :: KillLife                                                     
    

    !Management
    private ::      Ready
    private ::          LocateObjLife 
    
    !Interfaces----------------------------------------------------------------


    !Types---------------------------------------------------------------------


!*****_____________________________@ Start @____________________________ ******
    
    private :: T_Modoptions
        type       T_Modoptions
            logical                     :: Massxek          = .false.
        end type   T_Modoptions
    

    private :: T_ID
        type       T_ID
            integer                         :: ID
            character(len=StringLength)     :: Name
            character(len=StringLength)     :: Description
        end type   T_ID
    

    private :: T_Temperature
        type       T_Temperature
            real        :: T0       = FillValueReal             !Reference temperature for Q10 method
            real        :: Q10      = FillValueReal             !Q10 value for temperature limitation
            real        :: Reftemp  = FillValueReal             !Reference temperature
            real        :: Limit    = FillValueReal 
        end type   T_Temperature


    private :: T_Light
        type       T_Light
            real    :: Factor           = FillValueReal
            real    :: PhotoInhibition  = FillValueReal
            real    :: OptMin           = FillValueReal         !Minimum for optimal light
            real    :: OptMax           = FillValueReal         !Maximum for optimal light
            real    :: Adapt_rate       = FillValueReal         !Light adaptation rate
            real    :: MaxPhoto         = FillValueReal         !Maximum photosynthesys value
            real    :: AlphaChl         = FillValueReal         !Chl specific initial slop of P vs I curve       
            real    :: Chl_Synthesis    = FillValueReal   
            real    :: PhotoAclim       = FillValueReal
            real    :: Chl_Degrad_r     = FillValueReal         !Chl degradation rate constant       
        end type   T_Light


    private :: T_Nutrients
        type       T_Nutrients
            real    :: N_Status         =   FillValueReal
            real    :: P_Status         =   FillValueReal
            real    :: Si_Status        =   FillValueReal
            real    :: Chl_Status       =   FillValueReal
            real    :: Int_Nut_Limit    =   FillValueReal   
        end type   T_Nutrients


    private :: T_Limitation
        type       T_Limitation
            real                    :: Factor       = FillValueReal
            real                    :: Oxygen       = FillValueReal
            Type(T_Light      )     :: Light
            type(T_Temperature)     :: Temperature
            type(T_Nutrients  )     :: Nutrients
        end type   T_Limitation


    private :: T_Ratios
        type   T_Ratios
            real    :: NC_max       = FillValueReal             !Maximal N:C ratio
            real    :: NC_min       = FillValueReal             !Minimal N:C ratio
            real    :: NC_Actual    = FillValueReal
            real    :: PC_max       = FillValueReal             !Maximal N:C ratio
            real    :: PC_min       = FillValueReal             !Minimal N:C ratio
            real    :: PC_Actual    = FillValueReal
            real    :: ChlC_Actual  = FillValueReal
            real    :: ChlN_Max     = FillValueReal             !Maximal Chl:N ratio
            real    :: SiC_Actual   = FillValueReal
        end type T_Ratios


    private :: T_Pool
        type       T_Pool
            real    :: C    = FillValueReal
            real    :: N    = FillValueReal
            real    :: P    = FillValueReal
            real    :: Si   = FillValueReal
            real    :: Chl  = FillValueReal
        end type   T_Pool


    private :: T_Respiration
        type     T_Respiration
            real    :: Basal        = FillValueReal              !Basal respiration rate
            real    :: Frac_Prod    = FillValueReal              !Respired fraction of production
            real    :: Rate         = FillValueReal
            real    :: Carbon       = FillValueReal
            real    :: CO2          = FillValueReal
            real    :: Activity     = FillValueReal
            real    :: StandStock   = FillValueReal
            real    :: at10C        = FillValueReal              !rest respiration @ 10ºC
        end type T_Respiration


    private :: T_Exudation
        type     T_Exudation
            real            :: Nut_Stress   = FillValueReal     !Exudation under nutrient stress
            real            :: DOC_SL_Frac  = FillValueReal     !DOM diverted to semi-labile pool
            real            :: Rate         = FillValueReal
            real            :: DOC_SL       = FillValueReal             
            real            :: DOC_L        = FillValueReal
            type(T_Pool)    :: Rel_Exc
        end type T_Exudation


    private :: T_Excretion
        type T_Excretion
            real            :: Up_Frac      = FillValueReal      !excreted fraction of uptake
            real            :: Rate         = FillValueReal
            type(T_Pool)    :: Balance          
            type(T_Pool)    :: Rel_Exc  
            type(T_Pool)    :: DOM
        end type T_Excretion


    private :: T_Mortality
        type     T_Mortality
            real            :: Min_Lysis        = FillValueReal          !Minimal lysis rate
            real            :: Lysis_Ref_Con    = FillValueReal
            real            :: Ref_rate         = FillValueReal
            real            :: DOC_SL_Frac      = FillValueReal
            real            :: Density_Dep      = FillValueReal          !density-dependence mortality rate
            real            :: POM_Frac         = FillValueReal          !fraction of mortality to POM
            real            :: Rate             = FillValueReal          !temperature-independent mortality rate
            real            :: Lysis            = FillValueReal
            real            :: O2_Dep           = FillValueReal          !oxygen-dependent mortality rate
            real            :: O2_Dep_r         = FillValueReal
            real            :: O2_Ks            = FillValueReal          !oxygen half saturation constant 
            type(T_Pool)    :: POM          
            type(T_Pool)    :: DOM_Tot      
            type(T_Pool)    :: DOM_SL       
            type(T_Pool)    :: DOM_L        
        end type T_Mortality


    private :: T_Uptake
        type     T_Uptake
            real            :: Exc_Si           = FillValueReal         !Release rate of excess silicate
            real            :: Max_Ass_Rate     = FillValueReal         !Maximal assimilation rate
            real            :: Assimil          = FillValueReal
            real            :: Ass_Inc          = FillValueReal
            real            :: Ass_Net          = FillValueReal
            real            :: Max_Nut_Stor     = FillValueReal         !Maximal rate of storage filling
            real            :: Up_NH4_r         = FillValueReal         !Affinity for NH4 uptake 
            real            :: Up_NO3_r         = FillValueReal         !Affinity for NO3 uptake
            real            :: Up_PO4_r         = FillValueReal         !Affinity for PO4 uptake
            real            :: Carbon           = FillValueReal
            real            :: Up_NH4           = FillValueReal         
            real            :: Up_NO3           = FillValueReal         
            real            :: Up_PO4           = FillValueReal
            real            :: Up_Si_Ks         = FillValueReal         !Silicate uptake Michaelis constant
            real            :: Up_SiO           = FillValueReal
            real            :: P_Ext            = FillValueReal
            real            :: N_Ext            = FillValueReal
            real            :: NO3_Ext          = FillValueReal
            real            :: NH4_Ext          = FillValueReal
            Type(T_Pool)    :: Balance
            Type(T_Pool)    :: Q_Int        
        end type T_Uptake


    private :: T_Prey
        type    T_Prey
            type(T_ID)                  :: ID
            real                        :: Avail            = FillValueReal         !availability of Prey X
            real                        :: Pot              = FillValueReal
            logical                     :: Use_Silica       = .false.
            logical                     :: Use_Chl          = .false.
            type(T_Pool)                :: Frac
            type(T_Ratios)              :: Ratio
            type(T_Prey), pointer       :: Next
        end type T_Prey


    private :: T_Grazing
        type       T_Grazing
            real                        :: Vmax             = FillValueReal         !maximum specific uptake @ 10ºC
            real                        :: Ks               = FillValueReal         !half saturation value for uptake
            real                        :: Ass_Efic         = FillValueReal         !assimilation efficiency
            real                        :: Spec_Up          = FillValueReal
            real                        :: NC_Bal           = FillValueReal            
            real                        :: PC_Bal           = FillValueReal            
            type(T_Pool)                :: Total
            type(T_Pool)                :: Recyc
            type(T_Pool)                :: Assim
            type(T_Prey), pointer       :: FirstPrey
        end type T_Grazing


    private :: T_Movement
        type       T_Movement
            real                        :: Sed_min          = FillValueReal         !manimal sedimentation rate
            real                        :: Sed_nut_stress   = FillValueReal         !nutrient stress sedimentation rate
            real                        :: Sed_nut_thresh   = FillValueReal         !nutrient stress threshold
            type(T_Pool)                :: Sedimentation
        end type T_Movement



!___look @ here________________look @ here________________look @ here_____________
!
!   Given the unique set of parameters involved in bacterial uptake kinetics,
!   a diferent type (T_Consume) is defined only to decomposers. 
!
!___look @ here________________look @ here________________look @ here_____________

    private :: T_Consume
        type       T_Consume
            real            :: Vmax             = FillValueReal         !maximum specific uptake @ 10ºC
            real            :: Ks               = FillValueReal         !half saturation value for DOM uptake
            real            :: Ass_Efic_Norm    = FillValueReal         !Asimilation efficiency
            real            :: Ass_Efic_LowO2   = FillValueReal         !Asimilation efficiency @ low O2
            real            :: Ass_Efic         = FillValueReal
            real            :: O2_Low_Ass       = FillValueReal         ![oxygen] bollow which ass efic is low
            type(T_Pool)    :: Total
            real            :: DOM_Avail        = FillValueReal         !availability of DOM
            real            :: DOM_Up_Pot       = FillValueReal
            type(T_Pool)    :: DOM_Up_Frac
            real            :: PO4_Ks           = FillValueReal         !PO4 uptake Michaelis constant
            real            :: NH4_Ks           = FillValueReal         !NH4 uptake Michaelis constant
            real            :: NO3_Ks           = FillValueReal         !NO3 uptake Michaelis constant
            real            :: Up_PO4           = FillValueReal
            real            :: Up_NO3           = FillValueReal
            real            :: Up_NH4           = FillValueReal     
        end type T_Consume 
       
    private :: T_PoolIndex
        type       T_PoolIndex
            integer                                 :: Carbon      = null_int
            integer                                 :: Nitrogen    = null_int         
            integer                                 :: Phosphorus  = null_int         
            integer                                 :: Silica      = null_int
            integer                                 :: Chlorophyll = null_int        
        end type T_PoolIndex

    private :: T_Producer
        type       T_Producer      
            type(T_ID)                       :: ID
            type(T_PoolIndex)                :: PoolIndex
            logical                          :: Use_Silica
            logical                          :: Mixotrophy
            type(T_Limitation  )             :: Limitation
            type(T_Uptake      )             :: Uptake
            type(T_Respiration )             :: Respiration
            type(T_Excretion   )             :: Excretion
            type(T_Exudation   )             :: Exudation
            type(T_Mortality   )             :: Mortality
            type(T_Grazing     )             :: Grazing
            type(T_Ratios      )             :: Ratio
            type(T_Movement    )             :: Movement
            type(T_Producer    ), pointer    :: Next     
        end type    T_Producer


    private :: T_Consumer
        type       T_Consumer     
            type(T_ID)                       :: ID
            type(T_PoolIndex)                :: PoolIndex
            type(T_Limitation )              :: Limitation
            type(T_Respiration )             :: Respiration
            type(T_Excretion   )             :: Excretion
            type(T_Mortality   )             :: Mortality
            type(T_Grazing     )             :: Grazing
            type(T_Ratios      )             :: Ratio   
            type(T_Consumer    ), pointer    :: Next
        end type   T_Consumer


    private :: T_Decomposer
        type       T_Decomposer      
            type(T_ID)                       :: ID
            type(T_PoolIndex)                :: PoolIndex
            type(T_Limitation )              :: Limitation
            type(T_Consume     )             :: Consume
            type(T_Respiration )             :: Respiration
            type(T_Excretion   )             :: Excretion
            type(T_Mortality   )             :: Mortality
            type(T_Ratios      )             :: Ratio   
            type(T_Decomposer  ), pointer    :: Next
        end type    T_Decomposer



    private :: T_PropIndex
    type       T_PropIndex
         
        integer                                 :: Ammonia      = null_int        
        integer                                 :: Nitrate      = null_int        
        integer                                 :: Phosphate    = null_int
        integer                                 :: Silicate     = null_int
        integer                                 :: BioSilica    = null_int
        
        integer                                 :: DOC          = null_int         
        integer                                 :: DON          = null_int         
        integer                                 :: DOP          = null_int 
        integer                                 :: DOC_SL       = null_int         
        integer                                 :: DON_SL       = null_int         
        integer                                 :: DOP_SL       = null_int 

        integer                                 :: POC          = null_int         
        integer                                 :: PON          = null_int         
        integer                                 :: POP          = null_int
      
        integer                                 :: CarbonDioxide = null_int
        integer                                 :: Oxygen        = null_int   

    end type T_PropIndex



    private :: T_External
    type       T_External
        real, pointer, dimension(:  )       :: Salinity
        real, pointer, dimension(:  )       :: Temperature
        real, pointer, dimension(:  )       :: ShortWaveTop
        real, pointer, dimension(:  )       :: ShortWaveAverage
        real, pointer, dimension(:  )       :: LightExtCoefField
        real, pointer, dimension(:  )       :: Thickness
        real, pointer, dimension(:,:)       :: Mass
        integer, pointer, dimension(:  )    :: OpenPoints     
    end type T_External


    private :: T_BioChemParam
    type       T_BioChemParam   
        real                :: Redfield_NC      = null_real          !Redfield N:C ratio
        real                :: Redfield_PC      = null_real          !Redfield P:C ratio
        real                :: Redfield_SiC     = null_real          !Standard Si:C ratio
        real                :: BioSi_Diss       = null_real          !Biogenic silica dissolution rate
        real                :: O2C_Conversion   = null_real          !Oxygen to carbon conversion factor
        real                :: CO2C_Conversion  = null_real          !Carbon Dioxide to carbon conversion factor
        logical             :: CO2_on           = .false.            !Activates CO2 cycle
        real                :: POM_bac_ks       = null_real          !Bacteria mediated POM Hydrolysis MM cosntant
        real                :: POM_bac_vmax     = null_real          !Vmax for POM Hydrolysis
        real                :: DOMsl_bac_ks     = null_real          !Bacteria mediated DOMsl Hydrolysis MM cosntant
        real                :: DOMsl_bac_vmax   = null_real          !Vmax for DOMsl Hydrolysis
        real                :: Bac_Temp_Lim     = null_real
        real                :: Bac_conc         = null_real
        integer             :: L_lim_method     = null_int           !Light Limitatiom calculation method
        integer             :: T_lim_method     = null_int           !Temperature Limitatiom calculation method
        real                :: Water_abs_k      = null_real          !Absorption coefficient for water
        real                :: Chla_abs_k       = null_real          !Absorption coefficient for chlorophyll
        real                :: Nitrif           = null_real          
        real                :: Nitrifradiation  = null_real          !Light radiation bellow which nitrification occurs
        real                :: Nitrifrate       = null_real          !Nitrification rate
        real                :: Nit_Inib_coef    = null_real          !Nitrification inhibition coeficient
        real                :: Nit_ON_Conv      = null_real          !Nitrification O:N consumption ratio
        real                :: Nitrif_lim       = null_real
        type(T_Pool)        :: POM_bac_Hyd
        type(T_Pool)        :: DOMsl_bac_Hyd
    end type T_BioChemParam


    private :: T_Life
    type       T_Life
        integer                                      :: InstanceID
        real                                         :: DT
        real                                         :: DT_day
        integer                                      :: JulianDay
        integer, dimension(:), pointer               :: PropertyList
        type(T_Size1D       )                        :: Array
        type(T_Size1D       )                        :: Prop
        type(T_PropIndex    )                        :: PropIndex
        type(T_Producer     ), pointer               :: FirstProducer
        type(T_Consumer     ), pointer               :: FirstConsumer
        type(T_Decomposer   ), pointer               :: FirstDecomposer
        type(T_External     )                        :: ExternalVar
        type(T_BioChemParam )                        :: BioChemPar
        type(T_Modoptions   )                        :: Moptions
        integer                                      :: ObjEnterData = 0
        type(T_Life), pointer                        :: Next
    end type  T_Life


    !Global Module Variables
    type (T_Life), pointer                          :: FirstObjLife
    type (T_Life), pointer                          :: Me

    !--------------------------------------------------------------------------
    
     contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructLife(ObjLifeID, FileName, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjLifeID  
        character(len=StringLength)                     :: FileName
        integer, optional, intent(OUT)                  :: STAT

        !External----------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mLife_)) then
            nullify (FirstObjLife)
            call RegisterModule (mLife_) 
        endif
        
        call Ready(ObjLifeID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructLife - ModuleLife - ERROR #1'

            call ReadData

            call PropertyIndexNumber
        
            call ConstructPropertyList

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructLife - ModuleLife - ERROR #2'

            !Returns ID
            ObjLifeID          = Me%InstanceID

            STAT_ = SUCCESS_
        else 
            
            stop 'ModuleLife - ConstructLife - ERROR #1' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructLife

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Life), pointer                         :: NewObjLife
        type (T_Life), pointer                         :: PreviousObjLife


        !Allocates new instance
        allocate (NewObjLife)
        nullify  (NewObjLife%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjLife)) then
            FirstObjLife         => NewObjLife
            Me                   => NewObjLife
        else
            PreviousObjLife      => FirstObjLife
            Me                   => FirstObjLife%Next
            do while (associated(Me))
                PreviousObjLife  => Me
                Me               => Me%Next
            enddo
            Me                   => NewObjLife
            PreviousObjLife%Next => NewObjLife
        endif

        Me%InstanceID = RegisterNewInstance (mLIFE_)

    end subroutine AllocateInstance


!____________________________________________________________________________________________________
!___________________________________________________________________________________@_KEYWORD READING
    !--------------------------------------------------------------------------
    
    subroutine ReadData

        !Arguments-------------------------------------------------------------
                                                           
        !Local-----------------------------------------------------------------
        

        call ConstructGlobalVariables

        call ConstructModelOptions

        call ConstructProducers

        call ConstructConsumers
        
        call ConstructDecomposers


    end subroutine ReadData
 
    !--------------------------------------------------------------------------
 
    subroutine PropertyIndexNumber

        !Arguments-------------------------------------------------------------


        !Local-----------------------------------------------------------------
        type(T_Producer),      pointer             :: Producer
        type(T_Consumer),      pointer             :: Consumer
        type(T_Decomposer),    pointer             :: Decomposer
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
                Index                               = Index + 1
                Producer%PoolIndex%Nitrogen         = Index
                Index                               = Index + 1
                Producer%PoolIndex%Phosphorus       = Index
                Index                               = Index + 1
                Producer%PoolIndex%Chlorophyll      = Index
               
                if(Producer%Use_Silica)then
                    Index                           = Index + 1
                    Producer%PoolIndex%Silica       = Index
                    Me%Prop%IUB                     = Me%Prop%IUB + 5
                else
                    Me%Prop%IUB                     = Me%Prop%IUB + 4
                end if

                Producer => Producer%Next
            end do
        
        !Consumer index number
            Consumer => Me%FirstConsumer
            do while(associated(Consumer))
                
                Index                               = Index + 1
                Consumer%PoolIndex%Carbon           = Index
                Index                               = Index + 1
                Consumer%PoolIndex%Nitrogen         = Index
                Index                               = Index + 1
                Consumer%PoolIndex%Phosphorus       = Index

                Me%Prop%IUB        = Me%Prop%IUB + 3

                Consumer => Consumer%Next
            end do
   
        !Decomposer index number
            Decomposer => Me%FirstDecomposer
            do while(associated(Decomposer))
                
                Index                               = Index + 1
                Decomposer%PoolIndex%Carbon         = Index
                Index                               = Index + 1
                Decomposer%PoolIndex%Nitrogen       = Index
                Index                               = Index + 1
                Decomposer%PoolIndex%Phosphorus     = Index

                Me%Prop%IUB     = Me%Prop%IUB + 3

                Decomposer => Decomposer%Next
            end do


        !Nitrogen index number
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Ammonia        = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Nitrate        = Me%Prop%IUB


        !Phosphorus index number
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Phosphate      = Me%Prop%IUB

        !Silica index number
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Silicate       = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%BioSilica      = Me%Prop%IUB


        !OrganicMatter index number
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%DOC            = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%DON            = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%DOP            = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%DOC_SL         = Me%Prop%IUB
 
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%DON_SL         = Me%Prop%IUB
      
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%DOP_SL         = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%POC            = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%PON            = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%POP            = Me%Prop%IUB


     If (Me%BioChemPar%CO2_on) then 

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%CarbonDioxide  = Me%Prop%IUB

     end if

            
        !Oxygen index number
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Oxygen         = Me%Prop%IUB

        !----------------------------------------------------------------------

    end subroutine PropertyIndexNumber

    !--------------------------------------------------------------------------

    
    subroutine GetDTLife(Life_ID, DTDay, DT, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: Life_ID
        real,    optional, intent(OUT)      :: DTDay
        real,    optional, intent(OUT)      :: DT
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Life_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(DTDay   )) DTDay    = Me%DT_Day
            if (present(DT)) DT = Me%DT

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDTLife


   !--------------------------------------------------------------------------


    subroutine ConstructPropertyList

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        type(T_Producer),      pointer             :: Producer
        type(T_Consumer),      pointer             :: Consumer
        type(T_Decomposer),    pointer             :: Decomposer
        integer                                    :: Index
        !Local-----------------------------------------------------------------
        
        allocate(Me%PropertyList(Me%Prop%ILB: Me%Prop%IUB))
        
        Index = 0

        !Producer index number      
            Producer => Me%FirstProducer
            do while(associated(Producer))
                
                Me%PropertyList(Producer%PoolIndex%Carbon)     = &
                            GetPropertyIDNumber(trim(Producer%ID%name)//" carbon")
                Me%PropertyList(Producer%PoolIndex%Nitrogen)   = &
                            GetPropertyIDNumber(trim(Producer%ID%name)//" nitrogen")
                Me%PropertyList(Producer%PoolIndex%Phosphorus) = &
                            GetPropertyIDNumber(trim(Producer%ID%name)//" phosphorus")
                Me%PropertyList(Producer%PoolIndex%Chlorophyll)    = &
                            GetPropertyIDNumber(trim(Producer%ID%name)//" chlorophyll")

                if(Producer%Use_Silica)then
                    Me%PropertyList(Producer%PoolIndex%Silica)     = &
                            GetPropertyIDNumber(trim(Producer%ID%name)//" silica")
                end if


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
   
        !Decomposer index number
            Decomposer => Me%FirstDecomposer
            do while(associated(Decomposer))
                
                Me%PropertyList(Decomposer%PoolIndex%Carbon)     = &
                            GetPropertyIDNumber(trim(Decomposer%ID%name)//" carbon")
                Me%PropertyList(Decomposer%PoolIndex%Nitrogen)   = &
                            GetPropertyIDNumber(trim(Decomposer%ID%name)//" nitrogen")
                Me%PropertyList(Decomposer%PoolIndex%Phosphorus) = &
                            GetPropertyIDNumber(trim(Decomposer%ID%name)//" phosphorus")

                Decomposer => Decomposer%Next
            end do

        !Nitrogen index number
            Me%PropertyList(Me%PropIndex%Ammonia)     = Ammonia_
            Me%PropertyList(Me%PropIndex%Nitrate)     = Nitrate_

        !Phosphorus index number
            Me%PropertyList(Me%PropIndex%Phosphate)   = Inorganic_Phosphorus_

        !Silica index number
            Me%PropertyList(Me%PropIndex%Silicate)    = Silicate_
            Me%PropertyList(Me%PropIndex%BioSilica)   = BioSilica_
   
        !OrganicMatter index number
            Me%PropertyList(Me%PropIndex%DOC)         = DOC_
            Me%PropertyList(Me%PropIndex%DON)         = DON_
            Me%PropertyList(Me%PropIndex%DOP)         = DOP_

            Me%PropertyList(Me%PropIndex%DOC_SL)      = DOCsl_
            Me%PropertyList(Me%PropIndex%DON_SL)      = DONsl_
            Me%PropertyList(Me%PropIndex%DOP_SL)      = DOPsl_

            Me%PropertyList(Me%PropIndex%POC)         = POC_
            Me%PropertyList(Me%PropIndex%PON)         = PON_
            Me%PropertyList(Me%PropIndex%POP)         = POP_

     If (Me%BioChemPar%CO2_on) then 

            Me%PropertyList(Me%PropIndex%CarbonDioxide)         = CarbonDioxide_
     
     end if
            
        !Oxygen index number
            Me%PropertyList(Me%PropIndex%Oxygen)      = Oxygen_


        !----------------------------------------------------------------------

    end subroutine ConstructPropertyList


    !----------------------------------------------------------------------

    subroutine ConstructGlobalVariables

        !Arguments-------------------------------------------------------------
        
        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromFile 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromFile = FromFile)
        
        call GetData(Me%BioChemPar%Redfield_NC,             &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'REDFIELD_NC',          &
                     Default      = 0.012568,               &
                     ClientModule = 'ModuleLife',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleLife - ERROR #1'

        call GetData(Me%BioChemPar%Redfield_PC,             &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'REDFIELD_PC',          &
                     Default      = 0.000786,               &
                     ClientModule = 'ModuleLife',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleLife - ERROR #2'

        call GetData(Me%BioChemPar%Redfield_SiC,            &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'REDFIELD_SiC',         &
                     Default      = 0.03,                   &
                     ClientModule = 'ModuleLife',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleLife - ERROR #3'
        
        call GetData(Me%BioChemPar%BioSi_Diss,              &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'BIO_SI_DISS',          &
                     Default      = 0.01,                   &
                     ClientModule = 'ModuleLife',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleLife - ERROR #4'


        call GetData(Me%BioChemPar%O2C_Conversion,          &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'O2_CARB_CONVERS',      &
                     Default      = 2.664,                  &    !setup for biomass in mg C L-1
                     ClientModule = 'ModuleLife',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleLife - ERROR #5a'

        call GetData(Me%BioChemPar%CO2C_Conversion,         &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'CO2_CARB_CONVERS',     &
                     Default      = 3.6640579,              &   !setup for biomass in mg C L-1
                     ClientModule = 'ModuleLife',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleLife - ERROR #5b'

        call GetData(Me%BioChemPar%CO2_on,                  &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'CO2_ON',               &
                     Default      = .false.,                & 
                     ClientModule = 'ModuleLife',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleLife - ERROR #5c'


        call GetData(Me%DT,                                 &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'DT',                   &
                     Default      = 3600.,                  &
                     ClientModule = 'ModuleLife',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleLife - ERROR #6'

        call GetData(Me%BioChemPar%POM_bac_ks,              &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'POM_BAC_KS',           &
                     Default      = 200.,                   &
                     ClientModule = 'ModuleLife',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleLife - ERROR #7'

        call GetData(Me%BioChemPar%POM_bac_vmax,            &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'POM_BAC_VMAX',         &
                     Default      = 2.,                     &
                     ClientModule = 'ModuleLife',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleLife - ERROR #8'

        call GetData(Me%BioChemPar%DOMsl_bac_ks,            &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'DOMSL_BAC_KS',         &
                     Default      = 400.,                   &
                     ClientModule = 'ModuleLife',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleLife - ERROR #9'

        call GetData(Me%BioChemPar%DOMsl_bac_vmax,          &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'DOMSL_BAC_VMAX',       &
                     Default      = 2.,                     &
                     ClientModule = 'ModuleLife',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleLife - ERROR #10'
   
        call GetData(ME%BioChemPar%L_lim_method,                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'LIGHT_LIM_METHOD',                 &
                     Default      = 3,                                  &
                     ClientModule = 'ModuleLife',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables  - ModuleLife - ERROR #11'
        
        call GetData(ME%BioChemPar%T_lim_method,                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'TEMP_LIM_METHOD',                  &
                     Default      = 2,                                  &
                     ClientModule = 'ModuleLife',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables  - ModuleLife - ERROR #12' 

        call GetData(ME%BioChemPar%Nitrifradiation,                     &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'NITRIFRADLIM',                     &
                     Default      = 4.0,                                &
                     ClientModule = 'ModuleLife',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables  - ModuleLife - ERROR #13'
        
        call GetData(ME%BioChemPar%Nitrifrate,                          &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'NITRIFRATE',                       &
                     Default      = 0.04,                               &
                     ClientModule = 'ModuleLife',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables  - ModuleLife - ERROR #14' 

        call GetData(ME%BioChemPar%Nit_Inib_coef,                       &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'NIT_IN_COEF',                      &
                     Default      = 0.6,                                &
                     ClientModule = 'ModuleLife',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables  - ModuleLife - ERROR #15' 

        call GetData(ME%BioChemPar%Nit_ON_Conv,                         &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'NIT_O_N_CONV',                     &
                     Default      = 0.0588,                             &
                     ClientModule = 'ModuleLife',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables  - ModuleLife - ERROR #16' 


!***** ? *****

       !_________converts DT [sec] to DT [day] for model use
        Me%DT_day      = Me%DT / (3600.*24.)


    end subroutine ConstructGlobalVariables



!----------------------------------------------------------------------

    subroutine ConstructModelOptions

        !Arguments-------------------------------------------------------------
        
        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromFile 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromFile = FromFile)
        
        call GetData(Me%Moptions%massxek,                   &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'MASS_XEK',             &
                     Default      = .false.,                &
                     ClientModule = 'ModuleLife',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructModelOptions - ModuleLife - ERROR #1'


    end subroutine ConstructModelOptions




!_______________________________________________________________________________________________________
!_______________________________________________________________________________@_CONSTRUCT_PRODUCERS___ 

    
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

                    call ConstructProducerParameters    (NewProducer, ClientNumber)

                    nullify(NewProducer)

                else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop       'ConstructProducers - ModuleLife - ERROR #1'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                    stop       'ConstructProducers - ModuleLife - ERROR #2'
            else cd1
                    stop       'ConstructProducers - ModuleLife - ERROR #3'
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

    subroutine ConstructProducerParameters (NewProducer, ClientNumber)

        !Arguments-------------------------------------------------------------
        type (T_Producer),      pointer           :: NewProducer
        integer                                   :: ClientNumber

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
                     ClientModule = MohidModules(mLIFE_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleLife - ERROR #1'

        call GetData(NewProducer%ID%Description,                &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'DESCRIPTION',              &
                     ClientModule = MohidModules(mLIFE_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleLife - ERROR #2'

        call GetData(NewProducer%Use_Silica,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'SILICA_USE',               &
                     ClientModule = MohidModules(mLIFE_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleLife - ERROR #3'

        call GetData(NewProducer%Mixotrophy,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'MIXOTROPHY',               &
                     ClientModule = MohidModules(mLIFE_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleLife - ERROR #4'



        if(.not. CheckPropertyName(trim(NewProducer%ID%Name)//" carbon"))&
            stop 'ConstructProducerParameters - ModuleLife - ERROR #4'
          


        call ConstructProducerLimitation    (NewProducer                       )
        call ConstructProducerUptake        (NewProducer                       )
        call ConstructProducerExudation     (NewProducer                       )
        call ConstructRespiration           (NewProducer%Respiration           )
        call ConstructMortality             (NewProducer%Mortality             )
        call ConstructRatio                 (NewProducer%Ratio                 )
        call ConstructMovement              (NewProducer%Movement              )

        If (NewProducer%Mixotrophy) then
                call ConstructGrazing               (NewProducer%Grazing,  ClientNumber)
                call ConstructExcretion             (NewProducer%Excretion             )
        endif


    end subroutine ConstructProducerParameters


!_______________________________________________________________________________________________________
!________________________________________________________________________________@_CONSTRUCT_CONSUMERS__ 

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
                        stop       'ConstructConsumer - ModuleLife - ERROR #1'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                    stop       'ConstructConsumer - ModuleLife - ERROR #2'
            else cd1
                    stop       'ConstructConsumer - ModuleLife - ERROR #3'
            end if cd1
        end do do1

    end subroutine ConstructConsumers


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
                     ClientModule = MohidModules(mLIFE_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleLife - ERROR #1'

        if(.not. CheckPropertyName(trim(NewConsumer%ID%Name)//" carbon"))&
        stop 'ConstructConsumerParameters - ModuleLife - ERROR #2'

        call ConstructExcretion             (NewConsumer%Excretion             )
        call ConstructRespiration           (NewConsumer%Respiration           )
        call ConstructMortality             (NewConsumer%Mortality             )
        call ConstructTemperatureLimitation (NewConsumer%Limitation%Temperature)
        call ConstructGrazing               (NewConsumer%Grazing , ClientNumber)
        call ConstructRatio                 (NewConsumer%Ratio                 )


    end subroutine ConstructConsumerParameters

    !--------------------------------------------------------------------------




!_______________________________________________________________________________________________________
!________________________________________________________________________________@_CONSTRUCT DECOMPOSERS 
    !--------------------------------------------------------------------------

    subroutine ConstructDecomposers

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        type (T_Decomposer),      pointer         :: NewDecomposer
        integer                                   :: ClientNumber, STAT_CALL
        logical                                   :: BlockFound

        !Begin-----------------------------------------------------------------


do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<begin_decomposer>', &
                                        block_end       = '<end_decomposer>',   &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                    
                    call AddDecomposer                  (NewDecomposer)

                    call ConstructDecomposerParameters  (NewDecomposer)

                    nullify(NewDecomposer)

                else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop       'ConstructDecomposers - ModuleLife - ERROR #1'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                    stop       'ConstructDecomposers - ModuleLife - ERROR #2'
            else cd1
                    stop       'ConstructDecomposers - ModuleLife - ERROR #3'
            end if cd1
        end do do1

    end subroutine ConstructDecomposers


    !--------------------------------------------------------------------------


    subroutine AddDecomposer (ObjDecomposer)

        !Arguments-------------------------------------------------------------
        type (T_Decomposer),     pointer           :: ObjDecomposer
        !Local-----------------------------------------------------------------
        type (T_Decomposer),     pointer           :: PreviousDecomposer
        type (T_Decomposer),     pointer           :: NewDecomposer
        integer, save                              :: NextDecomposerID = 1

        !Allocates new Consumer
        allocate (NewDecomposer)
        nullify  (NewDecomposer%Next)

        !Insert new Consumer into list and makes current ?? point to it
        if (.not. associated(Me%FirstDecomposer)) then
            Me%FirstDecomposer            => NewDecomposer
            ObjDecomposer                 => NewDecomposer
        else
            PreviousDecomposer            => Me%FirstDecomposer
            ObjDecomposer                 => Me%FirstDecomposer%Next

            do while (associated(ObjDecomposer))
                PreviousDecomposer        => ObjDecomposer
                ObjDecomposer             => ObjDecomposer%Next
            enddo
            ObjDecomposer                 => NewDecomposer
            PreviousDecomposer%Next       => NewDecomposer
        endif

        !Attributes ID
        ObjDecomposer%ID%ID               = NextDecomposerID

        NextDecomposerID                  = NextDecomposerID + 1


    end subroutine AddDecomposer
    
    !--------------------------------------------------------------------------

    subroutine ConstructDecomposerParameters (NewDecomposer)

        !Arguments-------------------------------------------------------------
        type (T_Decomposer),    pointer           :: NewDecomposer

        !External--------------------------------------------------------------
        integer                                   :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                   :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(NewDecomposer%ID%Name,                     &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'NAME',                     &
                     ClientModule = MohidModules(mLIFE_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructDecomposerParameters - ModuleLife - ERROR #1'

        
        if(.not. CheckPropertyName(trim(NewDecomposer%ID%Name)//" carbon"))&
        stop 'ConstructDecomposerParameters - ModuleLife - ERROR #2'
                

        call ConstructDecomposerConsume     (NewDecomposer                       )
        call ConstructExcretion             (NewDecomposer%Excretion             )
        call ConstructRespiration           (NewDecomposer%Respiration           )
        call ConstructMortality             (NewDecomposer%Mortality             )
        call ConstructTemperatureLimitation (NewDecomposer%Limitation%Temperature)
        call ConstructRatio                 (NewDecomposer%Ratio                 )


    end subroutine ConstructDecomposerParameters

    !--------------------------------------------------------------------------




!_______________________________________________________________________________________________________
!____________________________________________________________________________@_START_READING_KEYWORDS___ 

    subroutine ConstructProducerLimitation (NewProducer)

        !Arguments-------------------------------------------------------------
        type (T_Producer),      pointer           :: NewProducer
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        call ConstructLightLimitation       (NewProducer)
        
        call ConstructTemperatureLimitation (NewProducer%Limitation%Temperature)
        

    end subroutine ConstructProducerLimitation
   
    !--------------------------------------------------------------------------
   
    subroutine ConstructLightLimitation (NewProducer)

        !Arguments-------------------------------------------------------------
        type (T_Producer),      pointer           :: NewProducer

        !External--------------------------------------------------------------
        integer                                   :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                   :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(NewProducer%Limitation%Light%PhotoInhibition,      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'PHOTOINHIBITION',                  &
                     ClientModule = MohidModules(mLIFE_)%Name,          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructLightLimitation - ModuleLife - ERROR #1'


        call GetData(NewProducer%Limitation%Light%AlphaChl,             &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'ALPHA_CHL',                        &
                     ClientModule = MohidModules(mLIFE_)%Name,          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructLightLimitation - ModuleLife - ERROR #2'


        call GetData(NewProducer%Limitation%Light%Chl_Degrad_r,         &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'CHL_DEGRAD_RATE',                  &
                     ClientModule = MohidModules(mLIFE_)%Name,          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructLightLimitation - ModuleLife - ERROR #3'
        
                                             
                                             
                                                                           

    end subroutine ConstructLightLimitation
 
    
   !--------------------------------------------------------------------------
 
   
    subroutine ConstructTemperatureLimitation (Temperature)

        !Arguments-------------------------------------------------------------
        type (T_Temperature)                        :: Temperature

        !External--------------------------------------------------------------
        integer                                     :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                     :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(Temperature%Q10,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'Q10',                              &
                     ClientModule = MohidModules(mLIFE_)%Name,          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTemperatureLimitation - ModuleLife - ERROR #1'

        call GetData(Temperature%Reftemp,                               &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'REF_TEMP',                         &
                     ClientModule = MohidModules(mLIFE_)%Name,          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTemperatureLimitation - ModuleLife - ERROR #2'

       call GetData(Temperature%T0,                                     &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'REF_TEMP_Q10',                     &
                     ClientModule = MohidModules(mLIFE_)%Name,          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTemperatureLimitation - ModuleLife - ERROR #3'



    end subroutine ConstructTemperatureLimitation

    !------------------------------------------------------------------

    subroutine ConstructProducerUptake (NewProducer)

        !Arguments-------------------------------------------------------------
        type (T_Producer)                                :: NewProducer
        
        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(NewProducer%Uptake%Exc_Si,                     &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'REL_EXCESS_SI',                &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerUptake - ModuleLife - ERROR #1'

        call GetData(NewProducer%Uptake%Max_Ass_Rate,               &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MAX_ASSIMIL',                  &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerUptake - ModuleLife - ERROR #2'

        call GetData(NewProducer%Uptake%Max_Nut_Stor,               &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MAX_STORE_FILL',               &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerUptake - ModuleLife - ERROR #3'

        call GetData(NewProducer%Uptake%Up_Si_Ks,                   &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'SI_UPTAKE_KS',                 &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerUptake - ModuleLife - ERROR #4'

        call GetData(NewProducer%Uptake%Up_NH4_r,                   &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'AFFINITY_NH4',                 &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerUptake - ModuleLife - ERROR #5'

        call GetData(NewProducer%Uptake%Up_NO3_r,                   &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'AFFINITY_NO3',                 &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerUptake - ModuleLife - ERROR #6'

        call GetData(NewProducer%Uptake%Up_PO4_r,                   &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'AFFINITY_PO4',                 &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerUptake - ModuleLife - ERROR #7'

        
  
    end subroutine ConstructProducerUptake

    !------------------------------------------------------------------

    subroutine ConstructProducerExudation (NewProducer)

        !Arguments-------------------------------------------------------------
        type (T_Producer)                                :: NewProducer
        
        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(NewProducer%Exudation%Nut_Stress,              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'EXU_NUT_STRESS',               &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerExudation - ModuleLife - ERROR #1'

        call GetData(NewProducer%Exudation%DOC_SL_Frac,             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'EXC_DOM_SL_FRAC',              &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerExudation - ModuleLife - ERROR #2'

               
    end subroutine ConstructProducerExudation
    
    !------------------------------------------------------------------

    subroutine ConstructRatio (Ratio)

        !Arguments-------------------------------------------------------------
        type (T_Ratios)                                 :: Ratio
        
        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(Ratio%NC_max,                                  &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MAX_NC_RATIO',                 &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructRatio - ModuleLife - ERROR #1'

        call GetData(Ratio%NC_min,                                  &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MIN_NC_RATIO',                 &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructRatio - ModuleLife - ERROR #2'

        call GetData(Ratio%PC_max,                                  &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MAX_PC_RATIO',                 &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructRatio - ModuleLife - ERROR #3'

        call GetData(Ratio%PC_min,                                  &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MIN_PC_RATIO',                 &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructRatio - ModuleLife - ERROR #4'

        call GetData(Ratio%ChlN_max,                                &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MAX_CHLN_RATIO',               &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructRatio - ModuleLife - ERROR #5'


    end subroutine ConstructRatio

    !--------------------------------------------------------------------------

    subroutine ConstructExcretion (Excretion)

        !Arguments-------------------------------------------------------------
        type (T_Excretion)                              :: Excretion
        
        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(Excretion%Up_Frac,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'EXCRE_UP_FRAC',                &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructExcretion - ModuleLife - ERROR #1'

        
    end subroutine ConstructExcretion

    !--------------------------------------------------------------------------

    subroutine ConstructRespiration (Respiration)

        !Arguments-------------------------------------------------------------
        type (T_Respiration)                            :: Respiration
        
        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(Respiration%Basal,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'RESP_BASAL',                   &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructRespiration - ModuleLife - ERROR #1'

        call GetData(Respiration%Frac_Prod,                         &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'RESP_FRAC_PROD',               &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructRespiration - ModuleLife - ERROR #2'
   
        call GetData(Respiration%at10C,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'REST_RESP_@10C',               &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructRespiration - ModuleLife - ERROR #3'

  
    end subroutine ConstructRespiration
    
    !--------------------------------------------------------------------------

    subroutine ConstructMortality (Mortality)

        !Arguments-------------------------------------------------------------
        type (T_Mortality)                              :: Mortality
        
        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(Mortality%Min_Lysis,                           &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MIN_LYSIS',                    &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMortality - ModuleLife - ERROR #1'

        call GetData(Mortality%POM_Frac,                            &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MORT_POM_FRAC',                &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMortality - ModuleLife - ERROR #2'
   
        call GetData(Mortality%Rate,                                &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MORT_RATE',                    &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMortality - ModuleLife - ERROR #3'

        call GetData(Mortality%O2_Dep,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MORT_O2_DEP',                  &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMortality - ModuleLife - ERROR #4'

        call GetData(Mortality%O2_Ks,                               &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'O2_KS',                        &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMortality - ModuleLife - ERROR #5'

        call GetData(Mortality%DOC_SL_Frac,                         &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MORT_DOM_SL_FRAC',             &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMortality - ModuleLife - ERROR #5'

        call GetData(Mortality%Lysis_Ref_Con,                       &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'LYS_REF_CON',                  &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMortality - ModuleLife - ERROR #5'

        call GetData(Mortality%Density_Dep,                         &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'DENS_DEP_MORT',                &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMortality - ModuleLife - ERROR #5'
  
    end subroutine ConstructMortality
    

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
        type (T_Prey), pointer                          :: NewPrey

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(Grazing%Vmax,                                  &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MAX_SPEC_UP_@10C',             &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleLife - ERROR #1'


        call GetData(Grazing%Ks,                                    &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'GRAZ_UP_KS',                   &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleLife - ERROR #2'


        call GetData(Grazing%Ass_Efic,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'ASSIMIL_EFFIC',                &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleLife - ERROR #3'

        do
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,         &
                                       '<begin_prey>', '<end_prey>',          &
                                       BlockInBlockFound,                     &
                                       FirstLine = FirstLine,                 &
                                       LastLine  = LastLine,                  &
                                       STAT      = STAT_CALL)
            if      (STAT_CALL .EQ. SUCCESS_) then    
            
                if (BlockInBlockFound) then
            
                    if     (((LastLine + 1) - (FirstLine - 1)) .GE. 1) then

                        call AddPrey        (Grazing, NewPrey)

                        call ConstructPrey  (NewPrey)

                        nullify(NewPrey)

                    else
                        write(*,*)  
                        write(*,*) 'Error counting preys. '
                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleLife - ERROR #4'
                    end if
                else

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleLife - ERROR #5'
                    
                    exit 

                end if 

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBlock. '
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleLife - ERROR #6'
            end if

        end do
        
    end subroutine ConstructGrazing
    
    
    !--------------------------------------------------------------------------


    subroutine ConstructPrey (NewPrey)

        !Arguments-------------------------------------------------------------
        type (T_Prey), pointer                          :: NewPrey

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL, iflag
        integer                                         :: FromBlockInBlock      

        !Local-----------------------------------------------------------------
        type(T_Producer), pointer                       :: Producer

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlockInBlock = FromBlockInBlock)

        call GetData(NewPrey%ID%Name,                               &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlockInBlock,               &
                     keyword      = 'NAME',                         &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPrey - ModuleLife - ERROR #1'
        
        
        if(.not. CheckPropertyName(trim(NewPrey%ID%Name)//" carbon"))&
        stop 'ConstructPrey - ModuleLife - ERROR #2'
                
        call GetData(NewPrey%Avail,                                 &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlockInBlock,               &
                     keyword      = 'GRAZ_AVAIL',                   &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPrey - ModuleLife - ERROR #2'


        Producer => Me%FirstProducer
        do while(associated(Producer))

            if(NewPrey%ID%Name == trim(Producer%ID%Name))then

                NewPrey%Use_Chl = .true.

                if(Producer%Use_Silica)NewPrey%Use_Silica = .true.

            end if

         

            Producer => Producer%Next
        end do

    end subroutine ConstructPrey
  
    !--------------------------------------------------------------------------


    subroutine AddPrey (Grazing, Prey)

        !Arguments-------------------------------------------------------------
        type (T_Grazing)                          :: Grazing
        type (T_Prey),          pointer           :: Prey

        !Local-----------------------------------------------------------------
        type (T_Prey),          pointer           :: PreviousPrey
        type (T_Prey),          pointer           :: NewPrey
        integer, save                             :: NextPreyID = 1

        !Allocates new Producer
        allocate (NewPrey)
        nullify  (NewPrey%Next)

        !Insert new Prey into list and makes current algae point to it
        if (.not. associated(Grazing%FirstPrey)) then
            NextPreyID              = 1
            Grazing%FirstPrey       => NewPrey
            Prey                    => NewPrey
        else
            PreviousPrey            => Grazing%FirstPrey
            Prey                    => Grazing%FirstPrey%Next

            do while (associated(Prey))
                PreviousPrey        => Prey
                Prey                => Prey%Next
            enddo
            Prey                    => NewPrey
            PreviousPrey%Next       => NewPrey
        endif

        !Attributes ID
        Prey%ID%ID               = NextPreyID

        NextPreyID               = NextPreyID + 1


    end subroutine AddPrey

    !------------------------------------------------------------------

    subroutine ConstructDecomposerConsume (NewDecomposer)

        !Arguments-------------------------------------------------------------
        type (T_Decomposer)                              :: NewDecomposer
        
        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(NewDecomposer%Consume%Ass_Efic_Norm,           &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'ASS_EFFIC',                    &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructDecomposerConsume - ModuleLife - ERROR #1'

        call GetData(NewDecomposer%Consume%Ass_Efic_LowO2,          &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'ASS_EFFIC_LOW_O2',             &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructDecomposerConsume - ModuleLife - ERROR #2'

        call GetData(NewDecomposer%Consume%Vmax,                    &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MAX_SPEC_UP_@10C',             &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructDecomposerConsume - ModuleLife - ERROR #3'

        call GetData(NewDecomposer%Consume%Ks,                      &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'DOM_UP_KS',                    &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructDecomposerConsume - ModuleLife - ERROR #4'

        call GetData(NewDecomposer%Consume%NH4_Ks,                  &
                     Me%ObjEnterData, iflag,                        &    
                     SearchType   = FromBlock,                      &
                     keyword      = 'NH4_Ks',                       &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructDecomposerConsume - ModuleLife - ERROR #5'

        call GetData(NewDecomposer%Consume%NO3_Ks,                  &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'NO3_Ks',                       &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructDecomposerConsume - ModuleLife - ERROR #6'

        call GetData(NewDecomposer%Consume%PO4_Ks,                  &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'PO4_Ks',                       &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructDecomposerConsume - ModuleLife - ERROR #7'

        call GetData(NewDecomposer%Consume%DOM_Avail,               &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'DOM_AVAIL',                    &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructDecomposerConsume - ModuleLife - ERROR #8'

        call GetData(NewDecomposer%Consume%O2_Low_Ass,              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'O2_LOW_ASS_EFIC',              &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructDecomposerConsume - ModuleLife - ERROR #9'

  
    end subroutine ConstructDecomposerConsume

    !------------------------------------------------------------------

    subroutine ConstructMovement (Movement)

        !Arguments-------------------------------------------------------------
        type (T_Movement)                               :: Movement
        
        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(Movement%Sed_min,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'SEDIM_MIN',                    &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMovement - ModuleLife - ERROR #1'

        call GetData(Movement%Sed_nut_stress,                       &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'SEDIM_NUT_STRESS',             &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMovement - ModuleLife - ERROR #2'

        call GetData(Movement%Sed_nut_thresh,                       &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'NUT_STRESS_THRESH',            &
                     ClientModule = MohidModules(mLIFE_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMovement - ModuleLife - ERROR #3'
   
    
    end subroutine ConstructMovement
    
    !--------------------------------------------------------------------------




    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    subroutine GetLifeSize(Life_ID, PropLB, PropUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: Life_ID
        integer, optional, intent(OUT)      :: PropLB,PropUB
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Life_ID, ready_)    
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PropLB   )) PropLB    = Me%Prop%ILB
            if (present(PropUB   )) PropUB    = Me%Prop%IUB

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetLifeSize

    !--------------------------------------------------------------------------

    subroutine GetLifePropertyList(Life_ID, PropertyList, STAT)

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

            call Read_Lock(mLIFE_, Me%InstanceID)

            PropertyList =>  Me%PropertyList

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetLifePropertyList

    !--------------------------------------------------------------------------
    
    subroutine GetLifePropIndex (Life_ID, PropIDNumber, PropertyIndex, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: Life_ID
        integer,           intent(IN )      :: PropIDNumber
        integer,           intent(OUT)      :: PropertyIndex
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_, CurrentIndex
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Life_ID, ready_)    
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            do CurrentIndex = Me%Prop%ILB, Me%Prop%IUB

                if (PropIDNumber == Me%PropertyList(CurrentIndex))then
                    PropertyIndex = CurrentIndex
                    STAT_ = SUCCESS_
                    exit
                else
                    STAT_ = NOT_FOUND_ERR_
                end if

            end do

        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))STAT = STAT_

    end subroutine GetLifePropIndex

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



    subroutine UnGetLife(ObjLifeID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjLifeID
        integer, dimension(:), pointer                  :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjLifeID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mLIFE_, Me%InstanceID, "UnGetLife3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetLife

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyLife(ObjLifeID,                &
                            Salinity,               & 
                            Temperature,            &
                            ShortWaveAverage,       &
                            ShortWaveTop,           &
                            LightExtCoefField,      &
                            Thickness,              &
                            Mass,                   &
                            OpenPoints,             &
                            ArraySize,              &
                            JulianDay,              &
                            STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjLifeID
        real,                 pointer, dimension(:  )   :: Salinity
        real,                 pointer, dimension(:  )   :: Temperature
        real,                 pointer, dimension(:  )   :: ShortWaveAverage
        real,                 pointer, dimension(:  )   :: ShortWaveTop
        real,                 pointer, dimension(:  )   :: LightExtCoefField
        real,                 pointer, dimension(:  )   :: Thickness
        real,                 pointer, dimension(:,:)   :: Mass
        integer, optional,    pointer, dimension(:  )   :: OpenPoints
        type(T_Size1D)                                  :: ArraySize
        integer,intent(IN)                              :: JulianDay
        integer, intent(OUT), optional                  :: STAT


        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        logical                                     :: CalcPoint, checkmass
        integer                                     :: index
        integer                                     :: Chunk
        !----------------------------------------------------------------------

        CHUNK = CHUNK_I(Me%Array%ILB, Me%Array%IUB)

        STAT_ = UNKNOWN_

        call Ready(ObjLifeID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            checkmass = .false.
            if(ME%JulianDay /= JulianDay)checkmass = .true.

            ME%JulianDay = JulianDay

            Me%ExternalVar%Salinity                   => Salinity
            if (.NOT. associated(Me%ExternalVar%Salinity))         &
                stop 'Subroutine ModifyLife; Module ModuleLife. ERR10' 

            Me%ExternalVar%Temperature                => Temperature
            if (.NOT. associated(Me%ExternalVar%Temperature))        &
                stop 'Subroutine ModifyLife; Module ModuleLife. ERR20'

            Me%ExternalVar%ShortWaveAverage           => ShortWaveAverage
            if (.NOT. associated(Me%ExternalVar%ShortWaveAverage)) &
                stop 'Subroutine ModifyLife; Module ModuleLife. ERR30'                 

            Me%ExternalVar%ShortWaveTop              => ShortWaveTop
            if (.NOT. associated(Me%ExternalVar%ShortWaveTop)) &
                stop 'Subroutine ModifyLife; Module ModuleLife. ERR40' 
            
            Me%ExternalVar%LightExtCoefField          => LightExtCoefField
            if (.NOT. associated(Me%ExternalVar%LightExtCoefField))  &
                stop 'Subroutine ModifyLife; Module ModuleLife. ERR50' 
            
            Me%ExternalVar%Thickness                  => Thickness
            if (.NOT. associated(Me%ExternalVar%Thickness))          &
                stop 'Subroutine ModifyLife; Module ModuleLife. ERR60'  

            Me%ExternalVar%Mass                       => Mass
            if (.NOT. associated(Me%ExternalVar%Mass))               &
                stop 'Subroutine ModifyLife; Module ModuleLife. ERR70'

            Me%ExternalVar%OpenPoints                 => OpenPoints
            if (.NOT. associated(Me%ExternalVar%OpenPoints))         &
               stop 'Subroutine ModifyLife; Module ModuleLife. ERR80'

            Me%Array%ILB = ArraySize%ILB
            Me%Array%IUB = ArraySize%IUB

            !$OMP PARALLEL PRIVATE(index)

            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
d1:         do index = Me%Array%ILB, Me%Array%IUB
            
                if (present(OpenPoints)) then
                    if (OpenPoints(index) == OpenPoint) then
                        CalcPoint = .true.
                    else
                        CalcPoint = .false.
                    endif
                else
                    CalcPoint = .true.
                endif


i1:             if (CalcPoint) then
                    !____________________________________________
                    !                           START ACTION.....ex nihilo
                    ! _____________________call organisms by name

                    ! Option: to check mass conservation during calc
                    ! file on the same foulder as .exe
                    if (Me%Moptions%massxek .and. checkmass)                            &
                        call MassConservationCheck (index)

                    call Producers            (index)
                        
                    call Decomposers          (index)
                       
                    call Consumers            (index)

                    call BioChemicalProcesses (index)
                               
                end if i1

            end do d1
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyLife



    !--------------------------------------------------------------------------


!_________________________________________
!_________________mass conservation test__

    subroutine MassConservationCheck (index)

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        type(T_Producer),      pointer             :: Producer
        type(T_Consumer),      pointer             :: Consumer
        type(T_Decomposer),    pointer             :: Decomposer

        integer :: PHY_N,    PHY_P,  PHY_Si 
        integer :: ZOO_N,    ZOO_P
        integer :: BAC_N,    BAC_P
        integer :: AM,       NA,      PO,     Si   
        integer :: DON,      DONsl,   PON 
        integer :: DOP,      DOPsl,   POP,    POSi 
        real    :: N_Count,  P_Count, Si_Count 
    !------------------------------------------------------------------------
       
        AM      = Me%PropIndex%Ammonia
        NA      = Me%PropIndex%Nitrate
        PO      = Me%PropIndex%Phosphate
        Si      = Me%PropIndex%Silicate
        DON     = Me%PropIndex%DON
        DOP     = Me%PropIndex%DOP
        DONsl   = Me%PropIndex%DON_SL
        DOPsl   = Me%PropIndex%DOP_SL
        PON     = Me%PropIndex%PON
        POP     = Me%PropIndex%POP
        POSi    = Me%PropIndex%BioSilica


        !__________________reset counter
        N_Count     = 0.
        P_Count     = 0.
        Si_Count    = 0. 


        !____________________________________________________producers count
        Producer => Me%FirstProducer
        do while(associated(Producer))

            PHY_N   = Producer%PoolIndex%Nitrogen
            PHY_P   = Producer%PoolIndex%Phosphorus

            N_Count = N_Count + Me%ExternalVar%Mass (PHY_N, Index)
            P_Count = P_Count + Me%ExternalVar%Mass (PHY_P, Index)

                if(Producer%Use_Silica)then
                    PHY_Si  = Producer%PoolIndex%Silica
                    Si_Count = Si_Count + Me%ExternalVar%Mass (PHY_Si, Index)
                end if

            Producer => Producer%Next
        end do


        !____________________________________________________consumers count
        Consumer => Me%FirstConsumer
        do while(associated(Consumer))

            ZOO_N   = Consumer%PoolIndex%Nitrogen
            ZOO_P   = Consumer%PoolIndex%Phosphorus

            N_Count = N_Count + Me%ExternalVar%Mass (ZOO_N, Index)
            P_Count = P_Count + Me%ExternalVar%Mass (ZOO_P, Index)

            Consumer => Consumer%Next
        end do


        !____________________________________________________decomposers count
        Decomposer => Me%FirstDecomposer
        do while(associated(Decomposer))

            BAC_N   = Decomposer%PoolIndex%Nitrogen
            BAC_P   = Decomposer%PoolIndex%Phosphorus

            N_Count = N_Count + Me%ExternalVar%Mass (BAC_N, Index)
            P_Count = P_Count + Me%ExternalVar%Mass (BAC_P, Index)

            Decomposer => Decomposer%Next
        end do


        !________________________________________________chemical species count


        N_Count = N_Count + Me%ExternalVar%Mass (AM, Index) +       &
                  Me%ExternalVar%Mass (NA, Index)           +       &
                  Me%ExternalVar%Mass (DON, Index)          +       &
                  Me%ExternalVar%Mass (DONsl, Index)        +       &
                  Me%ExternalVar%Mass (PON, Index) 

        P_Count = P_Count + Me%ExternalVar%Mass (PO, Index) +       &
                  Me%ExternalVar%Mass (DOP, Index)          +       &
                  Me%ExternalVar%Mass (DOPsl, Index)        +       &
                  Me%ExternalVar%Mass (POP, Index) 
    
        Si_Count = Si_Count                                 +       &
                   Me%ExternalVar%Mass (Si, Index)          +       &
                   Me%ExternalVar%Mass (POSi, Index) 


!_________________________________________
!___________________print test to screen__

!        IF (ME%JulianDay .eq. 1 ) then
!                    write (*,*) '  '
!                    write (*,*) 'Total N..............', N_Count, '[mmol N m-3]'
!                    write (*,*) 'Total P..............', P_Count, '[mmol P m-3]'
!                    write (*,*) 'Total Si.............', Si_Count,'[mmol Si m-3]' 
!        end if
!
        open(99,FILE='mass_check.dat', STATUS='NEW')
        write (99,69) ME%Julianday, N_Count, P_Count, Si_Count
        69 Format (1x, i3, 1x, f10.7, 1x, f10.7, 1x, f10.7)
!

    end subroutine MassConservationCheck
!--------------------------------------------------------------------------


!_________________________________________
!_________________________Producers calc__

    subroutine Producers (index)
    
    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        type(T_Producer),      pointer             :: Producer
        type(T_Prey),          pointer             :: Prey

        integer :: PHY_C,   PHY_N,  PHY_P,  PHY_Si, PHY_Chl
        integer :: PreyIndexC,  PreyIndexN, PreyIndexP, PreyIndexSi, PreyIndexChl
        integer :: AM,      NA,     PO,     Si   
        integer :: DOC,     DOCsl,  POC 
        integer :: DON,     DONsl,  PON 
        integer :: DOP,     DOPsl,  POP,    POSi 
        integer :: O2,      CO2
        integer :: aux_SiC, aux_NC, aux_PC
        real    :: Avail, Mortality
        real    :: Mix_Ex_DOC, Mix_Ex_DON, Mix_Ex_DOP
        real    :: Mix_Mort_POC, Mix_Mort_PON, Mix_Mort_POP
        real    :: Mix_Resp_Rate  
        real    :: AverageRadiation
    !------------------------------------------------------------------------
       
        AM      = Me%PropIndex%Ammonia
        NA      = Me%PropIndex%Nitrate
        PO      = Me%PropIndex%Phosphate
        Si      = Me%PropIndex%Silicate
        DOC     = Me%PropIndex%DOC
        DON     = Me%PropIndex%DON
        DOP     = Me%PropIndex%DOP
        DOCsl   = Me%PropIndex%DOC_SL
        DONsl   = Me%PropIndex%DON_SL
        DOPsl   = Me%PropIndex%DOP_SL
        POC     = Me%PropIndex%POC
        PON     = Me%PropIndex%PON
        POP     = Me%PropIndex%POP
        POSi    = Me%PropIndex%BioSilica
        O2      = Me%PropIndex%Oxygen

        If (Me%BioChemPar%CO2_on) then 

            CO2     = Me%PropIndex%CarbonDioxide

        end if


        AverageRadiation = Me%ExternalVar%ShortWaveAverage(index) 

        Producer => Me%FirstProducer

d1:     do while(associated(Producer))
            
            PHY_C   = Producer%PoolIndex%Carbon
            PHY_N   = Producer%PoolIndex%Nitrogen
            PHY_P   = Producer%PoolIndex%Phosphorus
            PHY_Chl = Producer%PoolIndex%Chlorophyll



            !_________________________________________
            !_______ckeck for Si use & Si limit calc__
              
i1:         if(Producer%Use_Silica)then

               PHY_Si  = Producer%PoolIndex%Silica

                                                                           ![units] >>  mmol si / mg C
               Producer%Ratio%SiC_Actual = Me%ExternalVar%Mass (PHY_Si, Index)/ &  
                                           Me%ExternalVar%Mass (PHY_C, Index)
               
               Producer%Limitation%Nutrients%Si_Status = Me%ExternalVar%Mass (Si, Index)/    &
                                                         (Me%ExternalVar%Mass (Si, Index)+   &
                                                         Producer%Uptake%Up_Si_Ks)
            end if i1


            !_________________________________________
            !__________________actual element ratios__

            Producer%Ratio%NC_Actual = Me%ExternalVar%Mass (PHY_N, Index)/      &
                                       Me%ExternalVar%Mass (PHY_C, Index)              ![units] >>  mmol N / mg C

            Producer%Ratio%PC_Actual = Me%ExternalVar%Mass (PHY_P, Index)/      &
                                       Me%ExternalVar%Mass (PHY_C, Index)              ![units] >>  mmol P / mg C
    
            Producer%Ratio%ChlC_Actual = Me%ExternalVar%Mass (PHY_Chl, Index)/      &
                                         Me%ExternalVar%Mass (PHY_C, Index)            ![units] >>  mg Chl / mg C


            !_________________________________________
            !___________________temp limitation calc__


s1:         Select case (ME%BioChemPar%T_lim_method)

            case (1)

                Producer%Limitation%Temperature%Limit = Temp_Lim_Q10 (Producer%Limitation%Temperature%Q10,  &
                                                        Producer%Limitation%Temperature%T0, index)

            case (2)
         
                Producer%Limitation%Temperature%Limit = exp(-4000.0 * ((1.0 / (Me%ExternalVar%Temperature (index)            &
                                                        + 273.15)) - (1.0 / (Producer%Limitation%Temperature%Reftemp + 273.15))))

            end select s1


            !_________________________________________
            !______________nutriente limitation calc__

            Producer%Limitation%Nutrients%N_Status = Nutrient_Lim (Producer%Ratio%NC_Actual,            & 
                                                     Producer%Ratio%NC_Min, Me%BiochemPar%Redfield_NC)
    
            Producer%Limitation%Nutrients%P_Status = Nutrient_Lim (Producer%Ratio%PC_Actual,            & 
                                                     Producer%Ratio%PC_Min, Me%BiochemPar%Redfield_PC)


            Producer%Limitation%Nutrients%Int_Nut_Limit = MIN (Producer%Limitation%Nutrients%N_Status,  &
                                                          Producer%Limitation%Nutrients%P_Status)

            !_________________________________________
            !__________________light limitation calc__

s2:         Select case (ME%BioChemPar%L_lim_method)

            case (1) s2

                !_________________________________________
                !___________________assimilation calc #1__                                      [units] >>  d-1

                if(Producer%Use_Silica)then

                        Producer%Limitation%Light%MaxPhoto= Producer%Uptake%Max_Ass_Rate *              &
                                                            Producer%Limitation%Temperature%Limit *     &
                                                            Producer%Limitation%Nutrients%Si_Status
                                                                   
                    else
                        Producer%Limitation%Light%MaxPhoto= Producer%Uptake%Max_Ass_Rate *              &
                                                            Producer%Limitation%Temperature%Limit
                                                    
                     
                end if


                if (Producer%Limitation%Light%MaxPhoto .eq. 0.0) then

                    Producer%Uptake%Assimil = 0.0

                else

                                                                                    !        [units] >>  d-1
                    Producer%Uptake%Assimil = Producer%Limitation%Light%MaxPhoto *                         &
                                              (1.0 - (exp((-1.0 * Producer%Limitation%Light%AlphaChl *     &  
                                              Producer%Ratio%ChlC_Actual *                                 &
                                              AverageRadiation)/                   &
                                              Producer%Limitation%Light%MaxPhoto)))                    
                endif   


            case (2) s2
                Producer%Limitation%Light%Factor =                                                      &
                PhytoLightLimitationFactor(Thickness       = Me%ExternalVar%Thickness(index),           &
                                           TopRadiation    = Me%ExternalVar%ShortWaveTop(index),        &
                                           PExt            = Me%ExternalVar%LightExtCoefField(index),   &
                                           Photoinhibition = Producer%Limitation%Light%Photoinhibition)


                !_________________________________________
                !___________________assimilation calc #1__                                      [units] >>  d-1

                if(Producer%Use_Silica)then

                    Producer%Uptake%Assimil =  Producer%Limitation%Temperature%Limit *      &
                                               Producer%Uptake%Max_Ass_Rate          *      &
                                               Producer%Limitation%Light%Factor      *      &
                                               Producer%Limitation%Nutrients%Si_Status
                                        
                else
                    Producer%Uptake%Assimil =  Producer%Limitation%Temperature%Limit *      &
                                               Producer%Limitation%Light%Factor      *      &
                                               Producer%Uptake%Max_Ass_Rate
                                       
                end if


            case (3) s2
                Producer%Limitation%Light%Factor = (0.4 + SIN((ME%JulianDay - 100.)*(PI/180.)) * 0.3)

                    !_________________________________________
                    !___________________assimilation calc #1__                                      [units] >>  d-1

                if(Producer%Use_Silica)then

                    Producer%Uptake%Assimil =  Producer%Limitation%Temperature%Limit *      &
                                               Producer%Uptake%Max_Ass_Rate          *      &
                                               Producer%Limitation%Light%Factor      *      &
                                               Producer%Limitation%Nutrients%Si_Status
                                        
                else
                    Producer%Uptake%Assimil =  Producer%Limitation%Temperature%Limit *      &
                                               Producer%Limitation%Light%Factor      *      &
                                               Producer%Uptake%Max_Ass_Rate

                                      
                end if


            case default s2
                Producer%Limitation%Light%Factor = (0.4 + SIN((ME%JulianDay - 100.)*(PI/180.)) * 0.3)

            end select s2




            !_________________________________________
            !_________________________DIC limitation__


!            If (Me%BioChemPar%CO2_on) then 
     
!            Producer%Uptake%Assimil_CO2 = Producer%Uptake%Assimil * (

!            Endif




            !_________________________________________
            !_________________mortality related calc__

    
            if (Producer%Ratio%NC_Actual .gt. 0.0  .AND. Producer%Ratio%PC_Actual .gt. 0.0) then

                Producer%Mortality%POM_Frac = MIN(1.0, MIN (Producer%Ratio%NC_Min / Producer%Ratio%NC_Actual,    &
                                              Producer%Ratio%PC_Min / Producer%Ratio%PC_Actual))

                else
                    Producer%Mortality%POM_Frac = 0.0
            endif



            Producer%Mortality%Rate = Producer%Mortality%Min_Lysis *    &
                                      (1. /(Producer%Limitation%Nutrients%Int_Nut_Limit + 0.1))    ![units] >>  d-1


    
            !________ .carbon fraction [part. & dissol(sl & l)]                            [units] >>  mg C m-3 d-1
            Producer%Mortality%POM%C = Producer%Mortality%Rate * Producer%Mortality%POM_Frac *  &
                                       Me%ExternalVar%Mass (PHY_C, Index)   


            Producer%Mortality%DOM_Tot%C = Producer%Mortality%Rate *(1. - Producer%Mortality%POM_Frac)*    &
                                           Me%ExternalVar%Mass (PHY_C, Index)


            Producer%Mortality%DOM_SL%C = Producer%Mortality%DOM_Tot%C * Producer%Mortality%DOC_SL_Frac

            Producer%Mortality%DOM_L%C = Producer%Mortality%DOM_Tot%C * (1. - Producer%Mortality%DOC_SL_Frac)


    
            !_______.nitrogen fraction [part. & dissol(sl & l)]                          [units] >>  mmol N m-3 d-1
            Producer%Mortality%POM%N = Producer%Mortality%Rate * Producer%Mortality%POM_Frac *  &
                                       Me%ExternalVar%Mass (PHY_N, Index)

            Producer%Mortality%DOM_Tot%N = Producer%Mortality%Rate *(1. - Producer%Mortality%POM_Frac)*    &
                                           Me%ExternalVar%Mass (PHY_N, Index)

            Producer%Mortality%DOM_SL%N = Producer%Mortality%DOM_Tot%N * Producer%Mortality%DOC_SL_Frac

            Producer%Mortality%DOM_L%N = Producer%Mortality%DOM_Tot%N * (1. - Producer%Mortality%DOC_SL_Frac)



            !_____.phosphorus fraction [part. & dissol(sl & l)]                          [units] >>  mmol P m-3 d-1
            Producer%Mortality%POM%P = Producer%Mortality%Rate * Producer%Mortality%POM_Frac *  &
                                       Me%ExternalVar%Mass (PHY_P, Index)

            Producer%Mortality%DOM_Tot%P = Producer%Mortality%Rate *(1. - Producer%Mortality%POM_Frac)*    &
                                           Me%ExternalVar%Mass (PHY_P, Index)

            Producer%Mortality%DOM_SL%P = Producer%Mortality%DOM_Tot%P * Producer%Mortality%DOC_SL_Frac

            Producer%Mortality%DOM_L%P = Producer%Mortality%DOM_Tot%P * (1. - Producer%Mortality%DOC_SL_Frac)



            !_________.silica fraction [part.]                                          [units] >>  mmol Si m-3 d-1
            if(Producer%Use_Silica)then

                Producer%Mortality%POM%Si = Producer%Mortality%Rate * Me%ExternalVar%Mass (PHY_Si, Index)

            end if



            !_________________________________________
            !_________________________exudation calc__                             
                                                                                                ! [units] >>  d-1
            Producer%Exudation%Rate = (Producer%Uptake%Assimil * (Producer%Exudation%Nut_Stress +       &
                                      (1. - Producer%Exudation%Nut_Stress)*                             &
                                      (1. - Producer%Limitation%Nutrients%Int_Nut_Limit)))          
   
                                                                                    ! [units] >>  mg C m-3 d-1
            Producer%Exudation%DOC_SL = Producer%Exudation%Rate * Producer%Exudation%DOC_SL_FRAC    &
                                        * Me%ExternalVar%Mass (PHY_C, Index)

            Producer%Exudation%DOC_L = Producer%Exudation%Rate * (1. - Producer%Exudation%DOC_SL_FRAC) &
                                   * Me%ExternalVar%Mass (PHY_C, Index)


            !_________________________________________
            !_______________________respiration calc__ 
                                                                                                ! [units] >>  d-1
            Producer%Respiration%Rate = (Producer%Respiration%Basal*Producer%Limitation%Temperature%Limit)  &
                                        +(Producer%Respiration%Frac_Prod*(Producer%Uptake%Assimil           &
                                        - Producer%Exudation%Rate))                             

                                                                                       ! [units] >>  mg C m-3 d-1
            Producer%Respiration%Carbon = Producer%Respiration%Rate * Me%ExternalVar%Mass (PHY_C, Index)


            !_________________________________________
            !___________________assimilation calc #2__                    

                                                                                            ! [units] >>  mg C m-3 d-1
            Producer%Uptake%Carbon = Producer%Uptake%Assimil * Me%ExternalVar%Mass (PHY_C, Index)

                                                                                        ! [units] >>  d-1
            Producer%Uptake%Ass_Inc = Producer%Uptake%Assimil - Producer%Exudation%Rate

            Producer%Uptake%Ass_Net = Producer%Uptake%Ass_Inc - Producer%Respiration%Rate



            !_________________________________________
            !_______________nutrient uptake kinetics__    


                !_____.internal quota [n & p]                               [units] >>  mmol (N or P) mg C-1 d-1
            !    Producer%Uptake%Q_Int%N = ((Producer%Uptake%Ass_Net*Producer%Ratio%NC_Max) +     &
            !                              (Producer%Ratio%NC_Max-Producer%Ratio%NC_Actual) *     &
            !                              Producer%Uptake%Max_Nut_Stor) * Me%ExternalVar%Mass (PHY_C, Index)


            if (Producer%Ratio%PC_Max .gt. Producer%Ratio%PC_Actual) then

                Producer%Uptake%Q_Int%P = ((Producer%Uptake%Ass_Net * Producer%Ratio%PC_Max) +     &
                                          (Producer%Ratio%PC_Max - Producer%Ratio%PC_Actual) *     &
                                          Producer%Uptake%Max_Nut_Stor) * Me%ExternalVar%Mass (PHY_C, Index)

            else
                Producer%Uptake%Q_Int%P = 0.0
            endif



            if (Producer%Ratio%NC_Max .gt. Producer%Ratio%NC_Actual) then

                Producer%Uptake%Q_Int%N = ((Producer%Uptake%Ass_Net * Producer%Ratio%NC_Max) +     &
                                          (Producer%Ratio%NC_Max - Producer%Ratio%NC_Actual) *     &
                                          Producer%Uptake%Max_Nut_Stor)
            else
                Producer%Uptake%Q_Int%N = 0.0
            endif



 

            !_____.external availability [n & p]                       [units] >>  mmol (N or P) mg C-1 d-1
            !_______.nitrogen
            Producer%Uptake%NH4_Ext = Producer%Uptake%Up_NH4_r * Me%ExternalVar%Mass (AM, Index)

            Producer%Uptake%NO3_Ext = Producer%Uptake%Up_NO3_r * Me%ExternalVar%Mass (NA, Index)



            Producer%Uptake%N_Ext = Producer%Uptake%NH4_Ext + Producer%Uptake%NO3_Ext
    
            !_______.phosphorus
            Producer%Uptake%P_Ext = Producer%Uptake%Up_PO4_r * Me%ExternalVar%Mass (PO, Index) * Me%ExternalVar%Mass (PHY_C, Index)

            !_____.uptake 

            !_______.nitrogen                                                       [units] >>  mmol N m-3 d-1
            Producer%Uptake%Balance%N = MIN (Producer%Uptake%Q_Int%N, Producer%Uptake%N_Ext)        
                                    

            !______.to avoid division by zero when N_Ext = 0
i2:         if (Producer%Uptake%N_Ext .le. 0.) then

                    Producer%Uptake%Up_NH4 = 0.
                    Producer%Uptake%Up_NO3 = 0.
    
            else i2

                if (((Producer%Uptake%NH4_Ext / Producer%Uptake%N_Ext) *        &
                   Producer%Uptake%Balance%N) .gt. 0.) then

                  Producer%Uptake%Up_NH4 = Min(((Producer%Uptake%NH4_Ext / Producer%Uptake%N_Ext) *           &
                                            Producer%Uptake%Balance%N) * Me%ExternalVar%Mass (PHY_C, Index),  &
                                            Me%ExternalVar%Mass (AM, Index)/ Me%DT_day)

                    else
                        Producer%Uptake%Up_NH4 = 0.
                end if


                if (((Producer%Uptake%NO3_Ext / Producer%Uptake%N_Ext) *        &
                   Producer%Uptake%Balance%N) .gt. 0.) then

                   Producer%Uptake%Up_NO3 = Min(((Producer%Uptake%NO3_Ext / Producer%Uptake%N_Ext) *          &
                                            Producer%Uptake%Balance%N) * Me%ExternalVar%Mass (PHY_C, Index),  &
                                            Me%ExternalVar%Mass (NA, Index)/ Me%DT_day)

                    else
                        Producer%Uptake%Up_NO3 = 0.
                end if

            endif i2


            !_______.phosphorus                                                     [units] >>  mmol P m-3 d-1
            Producer%Uptake%Balance%P = MIN (Producer%Uptake%Q_Int%P, Producer%Uptake%P_Ext)        
                                   

            if (Producer%Uptake%Balance%P .gt. 0.) then
                Producer%Uptake%Up_PO4 = Min(Producer%Uptake%Balance%P, Me%ExternalVar%Mass (PO, Index)/ Me%DT_day)
            else
                Producer%Uptake%Up_PO4 = 0.
            end if


            !_______.silica                                                       [units] >>  mmol Si m-3 d-1
i3:         if(Producer%Use_Silica)then

                if (Producer%Uptake%Ass_Net .gt. 0.) then

                  Producer%Uptake%Up_SiO = (MAX (0., Producer%Uptake%Ass_Net * Me%BiochemPar%Redfield_SiC)     &
                                           - (MAX (0., Producer%Ratio%SiC_Actual - Me%BiochemPar%Redfield_SiC) &
                                           * Producer%Uptake%Exc_Si)) * Me%ExternalVar%Mass (PHY_C, Index)

                else
                    Producer%Uptake%Up_SiO = 0.
                end if
            end if i3



            !_________________________________________
            !_________chlorophyll synthesis kinetics__   


i4:         if (ME%BioChemPar%L_lim_method .eq. 1) then

                if (AverageRadiation .gt. 0.0) then
                                                                                              ![units] >>  mg Chl / mmol N
                    Producer%Limitation%Light%Chl_Synthesis = Producer%Ratio%ChlN_Max * (Producer%Uptake%Assimil /      &
                                                              (Producer%Limitation%Light%AlphaChl *                     &  
                                                              Producer%Ratio%ChlC_Actual *                              &
                                                              AverageRadiation))
                else

                    Producer%Limitation%Light%Chl_Synthesis = 0.0

                end if
  

                                                                                              ![units] >>  mg Chl m-3 d-1
                                                                          !already considering Chl mass , real units: d-1
                Producer%Limitation%Light%PhotoAclim = (((Producer%Limitation%Light%Chl_Synthesis *                 &
                                                       ((Producer%Uptake%Up_NH4 + Producer%Uptake%Up_NO3) /         &
                                                       Me%ExternalVar%Mass (PHY_C, Index))) /                       &
                                                       Producer%Ratio%ChlC_Actual) -                                &
                                                       Producer%Limitation%Light%Chl_Degrad_r) *                    &
                                                       Me%ExternalVar%Mass (PHY_Chl, Index)

            endif i4


            !___________________________________________________________
            !_____________________________________myxotrophy processes__  
            !code adaptation (from consumers)

            !   If (Producer%Mixotrophy .AND. ((Producer%Ratio%NC_Actual .lt. Me%BioChemPar%Redfield_NC) .OR.       &
            !        (Producer%Ratio%PC_Actual .lt. Me%BioChemPar%Redfield_PC))) then

            
i5:         if (Producer%Mixotrophy) then
            
 
                !###TEST
                !   open(98,FILE='myxo_check.dat')
                !   write (98,68) ME%Julianday, Producer%Ratio%NC_Actual,    &
                !      Me%BioChemPar%Redfield_NC, Producer%Ratio%PC_Actual, Me%BioChemPar%Redfield_PC
                !   68 Format (1x, i3, 1x, f10.7, 1x, f10.7, 1x, f10.7, 1x, f10.7)



                !_________.set total predated fraction counter to zero         
                Producer%Grazing%Total%C = 0.
                Producer%Grazing%Total%N = 0.
                Producer%Grazing%Total%P = 0.

                Producer%Grazing%Spec_Up = 0. 



                !_________________________________________
                !__oxygen limitation calc.________________
                Producer%Limitation%Oxygen = Oxygen_Lim (Producer%Mortality%O2_Ks, Me%ExternalVar%Mass (O2, Index))


                !_________________________________________
                !________________________predation calc.__
                          
                Prey => Producer%Grazing%FirstPrey
d3:             do while(associated(Prey))

                    PreyIndexC  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" carbon"))
                    PreyIndexN  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" nitrogen"))
                    PreyIndexP  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" phosphorus"))

                    if(Prey%Use_Silica)then
                    PreyIndexSi = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" silica"))
                    end if

                    if(Prey%Use_Chl)then
                    PreyIndexChl  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" chlorophyll"))
                    end if


                    !_________.prey actual element ratio                          [units] >>  mmol (N & P) / mg C
                    Prey%Ratio%NC_Actual = Me%ExternalVar%Mass (PreyIndexN, Index) /      &
                                           Me%ExternalVar%Mass (PreyIndexC, Index)          

                    Prey%Ratio%PC_Actual = Me%ExternalVar%Mass (PreyIndexP, Index) /      &
                                           Me%ExternalVar%Mass (PreyIndexC, Index) 


                    !_________.prey potential grazing
                    Avail = Prey%Avail * Me%ExternalVar%Mass (PreyIndexC, Index)            ![units] >>  mg C m-3

                    Prey%Pot = Producer%Grazing%Vmax * (Avail/(Avail+Producer%Grazing%Ks)) *    &
                               Producer%Limitation%Temperature%Limit                             ![units] >>  d-1


                                                                                ![units] >>  mmol (N & P) m-3 d-1
                    Prey%Frac%N = Prey%Pot * Me%ExternalVar%Mass (PHY_C, Index) * Prey%Ratio%NC_Actual

                    Prey%Frac%P = Prey%Pot * Me%ExternalVar%Mass (PHY_C, Index) * Prey%Ratio%PC_Actual

                    Prey%Frac%C = Prey%Pot* Me%ExternalVar%Mass (PHY_C, Index)          ![units] >>  mg C m-3 d-1


                    Producer%Grazing%Total%C = Producer%Grazing%Total%C + Prey%Frac%C

                    Producer%Grazing%Total%N = Producer%Grazing%Total%N + Prey%Frac%N

                    Producer%Grazing%Total%P = Producer%Grazing%Total%P + Prey%Frac%P

                    Producer%Grazing%Spec_Up = Producer%Grazing%Spec_Up + Prey%Pot                ![units] >>  d-1


                    !___________________________________________________________________________                                                       
                    !_________.prey mass balance equations                                            [ matrix update ]

                    Me%ExternalVar%Mass (PreyIndexC, Index) = Me%ExternalVar%Mass (PreyIndexC, Index)           &
                                                              - (Prey%Frac%C) * Me%DT_day                   ! <- sinks

                    Me%ExternalVar%Mass (PreyIndexN, Index) = Me%ExternalVar%Mass (PreyIndexN, Index)           &
                                                              - (Prey%Frac%N) * Me%DT_day                   ! <- sinks

                    Me%ExternalVar%Mass (PreyIndexP, Index) = Me%ExternalVar%Mass (PreyIndexP, Index)           &
                                                              - (Prey%Frac%P) * Me%DT_day                   ! <- sinks              
    


                    if(ME%BioChemPar%L_lim_method .eq. 1 .AND. Prey%Use_Chl)then
                                                                                           ![units] >>  mg Chl / mg C
                        Prey%Ratio%ChlC_Actual = Me%ExternalVar%Mass (PreyIndexChl, Index) /      &
                                                 Me%ExternalVar%Mass (PreyIndexC, Index) 

                                                                                           ![units] >>  mmol (N & P) m-3 d-1
                        Prey%Frac%Chl = Prey%Pot * Me%ExternalVar%Mass (PHY_C, Index) * Prey%Ratio%ChlC_Actual


                        Me%ExternalVar%Mass (PreyIndexChl, Index) = Me%ExternalVar%Mass (PreyIndexChl, Index)       &
                                                                    - (Prey%Frac%Chl) * Me%DT_day                !<- sinks   

                    end if



                    if(Prey%Use_Silica)then

                        Prey%Ratio%SiC_Actual = Me%ExternalVar%Mass (PreyIndexSi, Index) /      &
                                                Me%ExternalVar%Mass (PreyIndexC, Index) 

                                                                                      ![units] >>  mmol Si m-3 d-1
                        Prey%Frac%Si = Prey%Pot * Me%ExternalVar%Mass (PHY_C, Index) * Prey%Ratio%SiC_Actual

                        Me%ExternalVar%Mass (PreyIndexSi, Index) = Me%ExternalVar%Mass (PreyIndexSi, Index)        &
                                                               - Prey%Frac%Si * Me%DT_day             ! <- sinks

                        Me%ExternalVar%Mass (POSi, Index) = Me%ExternalVar%Mass (POSi, Index)        &
                                                        + Prey%Frac%Si * Me%DT_day                    ! <- sources

                    end if

                    Prey => Prey%Next
                end do d3


                !_________________________________________
                !___________assimilated & recycled frac.__ 
                                            
                Producer%Grazing%NC_Bal = MIN (MAX ((Producer%Ratio%NC_Max - Producer%Ratio%NC_Actual) /     &
                                          (Producer%Ratio%NC_Max - Producer%Ratio%NC_Min), 0.), 1.)

                Producer%Grazing%PC_Bal = MIN (MAX ((Producer%Ratio%PC_Max - Producer%Ratio%PC_Actual) /     &
                                          (Producer%Ratio%PC_Max - Producer%Ratio%PC_Min), 0.), 1.)


                !____________assimilated                                                ![units] >>  mmol N m-3 d-1
                Producer%Grazing%Assim%N = Producer%Grazing%Total%N * Producer%Grazing%NC_Bal

                Producer%Grazing%Assim%P = Producer%Grazing%Total%P * Producer%Grazing%PC_Bal


                !____________recycled (to inorganic fraction)                           ![units] >>  mmol P m-3 d-1
                Producer%Grazing%Recyc%N = Producer%Grazing%Total%N * (1. - Producer%Grazing%NC_Bal)

                Producer%Grazing%Recyc%P = Producer%Grazing%Total%P * (1. - Producer%Grazing%PC_Bal)


                !_________________________________________
                !________________________mortality calc.__ 
                                                                                                           ![units] >>  d-1

                Producer%Mortality%O2_Dep_r = (1.  - Producer%Limitation%Oxygen) * Producer%Mortality%O2_Dep 
                                                                                           

                Mortality = Producer%Mortality%Rate + Producer%Mortality%O2_Dep_r


                !_________________________________________
                !________________respiration & excretion__
                                                                                                          ![units] >>  d-1
                Producer%Respiration%StandStock = Producer%Limitation%Temperature%Limit * Producer%Respiration%at10C

                Producer%Excretion%Rate = Producer%Grazing%Spec_Up *                                 & 
                                          (1.  - Producer%Grazing%Ass_Efic) * Producer%Excretion%Up_Frac

                Producer%Respiration%Activity = Producer%Grazing%Spec_Up *                                       & 
                                                (1.  - Producer%Grazing%Ass_Efic) * (1.  - Producer%Excretion%Up_Frac)

                                                                               
                Mix_Resp_Rate = (Producer%Respiration%StandStock + Producer%Respiration%Activity)   &
                                * Me%ExternalVar%Mass (PHY_C, Index)            ![units] >>  mg C m-3 d-1



                !________ .balance (mortality) to particulated  
                Mix_Mort_POC = Mortality * Producer%Mortality%POM_Frac              &
                                           * Me%ExternalVar%Mass (PHY_C, Index)         ![units] >>  mg C m-3 d-1

                Mix_Mort_PON = Mortality * Producer%Mortality%POM_Frac              &
                                           * Me%ExternalVar%Mass (PHY_N, Index)       ![units] >>  mmol N m-3 d-1
        
                Mix_Mort_POP = Mortality * Producer%Mortality%POM_Frac              &
                                           * Me%ExternalVar%Mass (PHY_P, Index)       ![units] >>  mmol P m-3 d-1


                !________ .balance (mortality + excretion) to dissolved 
                Mix_Ex_DOC = (Mortality * (1.  - Producer%Mortality%POM_Frac)         & 
                             + Producer%Excretion%Rate)                               &
                             * Me%ExternalVar%Mass (PHY_C, Index)         ![units] >>  mg C m-3 d-1

                Mix_Ex_DON = (Mortality * (1.  - Producer%Mortality%POM_Frac)          & 
                             + Producer%Excretion%Rate)                                &
                             * Me%ExternalVar%Mass (PHY_N, Index)       ![units] >>  mmol N m-3 d-1

                Mix_Ex_DOP = (Mortality * (1.  - Producer%Mortality%POM_Frac)          & 
                             + Producer%Excretion%Rate)                                &
                             * Me%ExternalVar%Mass (PHY_P, Index)       ![units] >>  mmol P m-3 d-1


                !___________________________________________________________
                !___________________________________mass balance equations__                            [ matrix update ]

                !_________.producers [c n p]
                Me%ExternalVar%Mass (PHY_C, Index) = Me%ExternalVar%Mass (PHY_C, Index)     &
                                                     + (Producer%Uptake%Carbon              &   ! <- sources
                                                     + Producer%Grazing%Total%C             &
                                                     - (Producer%Mortality%POM%C +          &   ! <- sinks
                                                     Producer%Mortality%DOM_SL%C +          &
                                                     Producer%Mortality%DOM_L%C  +          &
                                                     Producer%Exudation%DOC_SL   +          &
                                                     Producer%Exudation%DOC_L    +          &
                                                     Producer%Respiration%Carbon +          &
                                                     Mix_Mort_POC                +          &
                                                     Mix_Ex_DOC                  +          &
                                                     Mix_Resp_Rate)) * Me%DT_day 

                Me%ExternalVar%Mass (PHY_N, Index) = Me%ExternalVar%Mass (PHY_N, Index)     &
                                                     + (Producer%Uptake%Up_NH4 +            &   ! <- sources
                                                     Producer%Uptake%Up_NO3    +            &
                                                     Producer%Grazing%Assim%N               &
                                                     - (Producer%Mortality%POM%N +          &   ! <- sinks
                                                     Producer%Mortality%DOM_SL%N +          &
                                                     Producer%Mortality%DOM_L%N  +          &
                                                     Mix_Mort_PON                +          &
                                                     Mix_Ex_DON))                           &
                                                     * Me%DT_day 
                                             
                Me%ExternalVar%Mass (PHY_P, Index) = Me%ExternalVar%Mass (PHY_P, Index)     &
                                                     + (Producer%Uptake%Up_PO4 +            &   ! <- sources
                                                     Producer%Grazing%Assim%P               &
                                                     - (Producer%Mortality%POM%P +          &   ! <- sinks
                                                     Producer%Mortality%DOM_SL%P +          &
                                                     Producer%Mortality%DOM_L%P  +          & 
                                                     Mix_Mort_POP                +          &
                                                     Mix_Ex_DOP))                           &
                                                     * Me%DT_day
                                    
                if (ME%BioChemPar%L_lim_method .eq. 1) then

                    Me%ExternalVar%Mass (PHY_Chl, Index) = Me%ExternalVar%Mass (PHY_Chl, Index)     &
                                                         + (Producer%Limitation%Light%PhotoAclim    &   ! <- sources
                                                         - ((Producer%Mortality%POM%C  +            &   ! <- sinks
                                                         Producer%Mortality%DOM_SL%C   +            &
                                                         Producer%Mortality%DOM_L%C    +            & 
                                                         Mix_Mort_POC                  +            &
                                                         Mix_Ex_DOC)                                &
                                                         * Producer%Ratio%ChlC_Actual))             &
                                                         * Me%DT_day
                endif



                !_________.inorganic nutrients [n p]

                Me%ExternalVar%Mass (AM, Index) = Me%ExternalVar%Mass (AM, Index)           &
                                                  + (Producer%Grazing%Recyc%N               &
                                                  - Producer%Uptake%Up_NH4) * Me%DT_day         ! <- sinks
                                        
                Me%ExternalVar%Mass (PO, Index) = Me%ExternalVar%Mass (PO, Index)           &
                                                  + (Producer%Grazing%Recyc%P               &
                                                  - Producer%Uptake%Up_PO4) * Me%DT_day         ! <- sinks
                                        

                !_________.diss organic matter [c n p]
                Me%ExternalVar%Mass (DOC, Index) = Me%ExternalVar%Mass (DOC, Index)   &
                                                   + (Producer%Mortality%DOM_L%C      &
                                                   + Mix_Ex_DOC                       &
                                                   + Producer%Exudation%DOC_L) * Me%DT_day

                Me%ExternalVar%Mass (DON, Index) = Me%ExternalVar%Mass (DON, Index)   &
                                                   + (Producer%Mortality%DOM_L%N      &
                                                   + Mix_Ex_DON) * Me%DT_day

                Me%ExternalVar%Mass (DOP, Index) = Me%ExternalVar%Mass (DOP, Index)   &
                                                   + (Producer%Mortality%DOM_L%P      &
                                                   + Mix_Ex_DOP) * Me%DT_day


                !_________.part organic matter [c n p]
                Me%ExternalVar%Mass (POC, Index) = Me%ExternalVar%Mass (POC, Index)   &
                                                   + (Producer%Mortality%POM%C        &
                                                   + Mix_Mort_POC) * Me%DT_day

                Me%ExternalVar%Mass (PON, Index) = Me%ExternalVar%Mass (PON, Index)   &
                                                   + (Producer%Mortality%POM%N        &
                                                   + Mix_Mort_PON) * Me%DT_day

                Me%ExternalVar%Mass (POP, Index) = Me%ExternalVar%Mass (POP, Index)   &
                                                   + (Producer%Mortality%POM%P        &
                                                   + Mix_Mort_POP) * Me%DT_day


                !_________.oxygen [O]                                       [units] >>  mg O2 L-1 d-1

                Me%ExternalVar%Mass (O2, Index) = Me%ExternalVar%Mass (O2, Index)                &
                                                 + (Producer%Uptake%Carbon                       &
                                                 - Producer%Respiration%Carbon                   &
                                                 - Mix_Resp_Rate) * Me%BioChemPar%O2C_Conversion * Me%DT_day
                    

            !from initial myxotrophy condition....  
            else i5

                !___________________________________________________________
                !___________________________________mass balance equations__                            [ matrix update ]

                !_________.producers [c n p]
                Me%ExternalVar%Mass (PHY_C, Index) = Me%ExternalVar%Mass (PHY_C, Index)     &
                                                     + (Producer%Uptake%Carbon              &   ! <- sources
                                                     - (Producer%Mortality%POM%C +          &   ! <- sinks
                                                     Producer%Mortality%DOM_SL%C +          &
                                                     Producer%Mortality%DOM_L%C  +          &
                                                     Producer%Exudation%DOC_SL   +          &
                                                     Producer%Exudation%DOC_L    +          &
                                                     Producer%Respiration%Carbon)) * Me%DT_day 

                Me%ExternalVar%Mass (PHY_N, Index) = Me%ExternalVar%Mass (PHY_N, Index)     &
                                                     + (Producer%Uptake%Up_NH4 +            &   ! <- sources
                                                     Producer%Uptake%Up_NO3                 &
                                                     - (Producer%Mortality%POM%N +          &   ! <- sinks
                                                     Producer%Mortality%DOM_SL%N +          &
                                                     Producer%Mortality%DOM_L%N))           &
                                                     * Me%DT_day

                Me%ExternalVar%Mass (PHY_P, Index) = Me%ExternalVar%Mass (PHY_P, Index)     &
                                                     + (Producer%Uptake%Up_PO4              &   ! <- sources
                                                     - (Producer%Mortality%POM%P +          &   ! <- sinks
                                                     Producer%Mortality%DOM_SL%P +          &
                                                     Producer%Mortality%DOM_L%P))           &
                                                     * Me%DT_day


                if (ME%BioChemPar%L_lim_method .eq. 1) then

                    Me%ExternalVar%Mass (PHY_Chl, Index) = Me%ExternalVar%Mass (PHY_Chl, Index)     &
                                                         + (Producer%Limitation%Light%PhotoAclim    &   ! <- sources
                                                         - ((Producer%Mortality%POM%C +             &   ! <- sinks
                                                         Producer%Mortality%DOM_SL%C +              &
                                                         Producer%Mortality%DOM_L%C) *              &
                                                         Producer%Ratio%ChlC_Actual))               &
                                                         * Me%DT_day
                endif



                !_________.inorganic nutrients [n p]

                Me%ExternalVar%Mass (AM, Index) = Me%ExternalVar%Mass (AM, Index)           &
                                                  - Producer%Uptake%Up_NH4 * Me%DT_day          ! <- sinks
                                        
                Me%ExternalVar%Mass (PO, Index) = Me%ExternalVar%Mass (PO, Index)           &
                                                  - Producer%Uptake%Up_PO4 * Me%DT_day          ! <- sinks
                                        

                !_________.diss organic matter [c n p]
                Me%ExternalVar%Mass (DOC, Index) = Me%ExternalVar%Mass (DOC, Index)   &
                                                   + (Producer%Mortality%DOM_L%C +    &
                                                   Producer%Exudation%DOC_L) * Me%DT_day

                Me%ExternalVar%Mass (DON, Index) = Me%ExternalVar%Mass (DON, Index)   &
                                                   + Producer%Mortality%DOM_L%N * Me%DT_day

                Me%ExternalVar%Mass (DOP, Index) = Me%ExternalVar%Mass (DOP, Index)   &
                                                   + Producer%Mortality%DOM_L%P * Me%DT_day


                !_________.part organic matter [c n p]
                Me%ExternalVar%Mass (POC, Index) = Me%ExternalVar%Mass (POC, Index)   &
                                                   + Producer%Mortality%POM%C * Me%DT_day

                Me%ExternalVar%Mass (PON, Index) = Me%ExternalVar%Mass (PON, Index)   &
                                                   + Producer%Mortality%POM%N * Me%DT_day

                Me%ExternalVar%Mass (POP, Index) = Me%ExternalVar%Mass (POP, Index)   &
                                                   + Producer%Mortality%POM%P * Me%DT_day


                !_________.oxygen [O]                                       [units] >>  mg O2 L-1 d-1

                Me%ExternalVar%Mass (O2, Index) = Me%ExternalVar%Mass (O2, Index)                &
                                                 + (Producer%Uptake%Carbon                       &
                                                 - Producer%Respiration%Carbon) * Me%BioChemPar%O2C_Conversion * Me%DT_day



                !_________.carbon dioxide [CO2]                            [units] >>  mg CO2 L-1 d-1

                 If (Me%BioChemPar%CO2_on) then 

                    Me%ExternalVar%Mass (CO2, Index) = Me%ExternalVar%Mass (CO2, Index)           &
                                                      + (Producer%Respiration%Carbon              &
                                                      - Producer%Uptake%Carbon) * Me%BioChemPar%CO2C_Conversion * Me%DT_day
                 Endif





            !from initial myxotrophy condition....
            endif i5

                
            !___________________________________________________________
            !___________________________mass balance equations (cont.)__                            [ matrix update ]


            !_________.inorganic nutrients [n p]
            Me%ExternalVar%Mass (NA, Index) = Me%ExternalVar%Mass (NA, Index)           &
                                              - Producer%Uptake%Up_NO3 * Me%DT_day          ! <- sinks


            !_________.silica [si]
i6:         if(Producer%Use_Silica)then

                !_________.producers
                Me%ExternalVar%Mass (PHY_Si, Index) = Me%ExternalVar%Mass (PHY_Si, Index)   &
                                                      + (Producer%Uptake%Up_SiO             &   ! <- sources
                                                      - Producer%Mortality%POM%Si)          &   ! <- sinks
                                                      * Me%DT_day
                
                !_________.biogenic silica
                Me%ExternalVar%Mass (POSi, Index) = Me%ExternalVar%Mass (POSi, Index)       &
                                                    + Producer%Mortality%POM%Si * Me%DT_day     ! <- sources

                !_________.silicate acid
                Me%ExternalVar%Mass (Si, Index) = Me%ExternalVar%Mass (Si, Index)           &
                                                  - Producer%Uptake%Up_SiO * Me%DT_day          ! <- sinks
            end if i6


            !_________.diss organic matter [c n p]
            Me%ExternalVar%Mass (DOCsl, Index) = Me%ExternalVar%Mass (DOCsl, Index)     &
                                                 + (Producer%Mortality%DOM_SL%C +       &
                                                 Producer%Exudation%DOC_SL) * Me%DT_day

            Me%ExternalVar%Mass (DONsl, Index) = Me%ExternalVar%Mass (DONsl, Index)   &
                                                 + Producer%Mortality%DOM_SL%N * Me%DT_day

            Me%ExternalVar%Mass (DOPsl, Index) = Me%ExternalVar%Mass (DOPsl, Index)   &
                                                 + Producer%Mortality%DOM_SL%P * Me%DT_day


            !_________________________________________
            !______________nutrient excess exudation__  


            !_________________________________________
            !__________________updated element ratios__

            aux_NC = Me%ExternalVar%Mass (PHY_N, Index) / Me%ExternalVar%Mass (PHY_C, Index)     ![units] >>  mmol N / mg C

            aux_PC = Me%ExternalVar%Mass (PHY_P, Index) / Me%ExternalVar%Mass (PHY_C, Index)     ![units] >>  mmol P / mg C
    

            !_______.nitrogen                                                           [units] >>  mmol N m-3
            If (aux_NC .gt. Producer%Ratio%NC_Max) then

                Producer%Exudation%Rel_Exc%N = (aux_NC - Producer%Ratio%NC_Max) *    &
                                                Me%ExternalVar%Mass (PHY_C, Index) 
            
            else
                Producer%Exudation%Rel_Exc%N = 0.
            end if


            !_______.phosphorus                                                         [units] >>  mmol P m-3
            if (aux_PC .gt. Producer%Ratio%PC_Max) then

                Producer%Exudation%Rel_Exc%P = (aux_PC - Producer%Ratio%PC_Max) *    &
                                               Me%ExternalVar%Mass (PHY_C, Index) 
            
            else
                Producer%Exudation%Rel_Exc%P = 0.
            end if 
       

            !_____.silica                                                             [units] >>  mmol Si m-3
i7:         if(Producer%Use_Silica)then

                aux_SiC = Me%ExternalVar%Mass (PHY_Si, Index) / Me%ExternalVar%Mass (PHY_C, Index)
                   

                if (aux_SiC .gt. Me%BiochemPar%Redfield_SiC) then

                    Producer%Exudation%Rel_Exc%Si = (Producer%Ratio%SiC_Actual - Me%BiochemPar%Redfield_SiC)  &
                                                    * Me%ExternalVar%Mass (PHY_C, Index)

                else
                    Producer%Exudation%Rel_Exc%Si = 0.
                end if

            end if i7



            !___________________________________________________________
            !___________________________update mass balance equations __                            [ matrix update ]


            Me%ExternalVar%Mass (PHY_N, Index) = Me%ExternalVar%Mass (PHY_N, Index)     &
                                                 - Producer%Exudation%Rel_Exc%N


            Me%ExternalVar%Mass (PHY_P, Index) = Me%ExternalVar%Mass (PHY_P, Index)     &
                                                 - Producer%Exudation%Rel_Exc%P



            !_________.inorganic nutrients [n p]

            Me%ExternalVar%Mass (AM, Index) = Me%ExternalVar%Mass (AM, Index)           &
                                              + Producer%Exudation%Rel_Exc%N            ! <- sources  
                        

            Me%ExternalVar%Mass (PO, Index) = Me%ExternalVar%Mass (PO, Index)           &
                                              + Producer%Exudation%Rel_Exc%P            ! <- sources  
                        

            !_________.silica [si]
            if(Producer%Use_Silica)then

                !_________.producers
                Me%ExternalVar%Mass (PHY_Si, Index) = Me%ExternalVar%Mass (PHY_Si, Index)   &
                                                      - Producer%Exudation%Rel_Exc%Si
    
                !_________.silicate acid
                Me%ExternalVar%Mass (Si, Index) = Me%ExternalVar%Mass (Si, Index)           &
                                                  + Producer%Exudation%Rel_Exc%Si

            end if

            


            !__________________________
            !____TO BE IMPLEMENTED_____

            !_____________________________________________________________________
            !___________________Sedimentation rate                          [units] >>  m d-1

            !    If (Producer%Use_Silica)then

            !        Producer%Movement%Sedimentation%C = Producer%Movement%Sed_nut_stress                   &
            !                                            * max(0., Producer%Movement%Sed_nut_thresh -       &
            !                                            (min(Producer%Limitation%Nutrients%Int_Nut_Limit,  &
            !                                            Producer%Limitation%Nutrients%Si_Status)))         &
            !                                            + Producer%Movement%Sed_min
        
            !        Producer%Movement%Sedimentation%Si = Producer%Movement%Sedimentation%C *         &
            !                                             Producer%Ratio%SiC_Actual 


            !        else

            !            Producer%Movement%Sedimentation%C = Producer%Movement%Sed_nut_stress               &
            !                                                * max(0., Producer%Movement%Sed_nut_thresh     &
            !                                                - Producer%Limitation%Nutrients%Int_Nut_Limit) &
            !                                                + Producer%Movement%Sed_min

            !    endif


            !    Producer%Movement%Sedimentation%N = Producer%Movement%Sedimentation%C *         &
            !                                        Producer%Ratio%NC_Actual 

            !    Producer%Movement%Sedimentation%P = Producer%Movement%Sedimentation%C *         &
            !                                        Producer%Ratio%PC_Actual 

            !    If (ME%BioChemPar%L_lim_method .eq. 1) then

            !        Producer%Movement%Sedimentation%Chl = Producer%Movement%Sedimentation%C *         &
            !                                              Producer%Ratio%ChlC_Actual
    
            !    Endif


            Producer => Producer%Next
        end do d1

    end subroutine Producers
 !--------------------------------------------------------------------------




!_________________________________________
!_________________________Consumers calc__

    subroutine Consumers (index)
    
    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        type(T_Consumer),      pointer             :: Consumer
        type(T_Prey),          pointer             :: Prey

        integer :: ZOO_C,       ZOO_N,      ZOO_P
        integer :: PreyIndexC,  PreyIndexN, PreyIndexP, PreyIndexSi, PreyIndexChl
        integer :: AM,          PO   
        integer :: DOC,         POC 
        integer :: DON,         PON 
        integer :: DOP,         POP,        POSi 
        integer :: O2,          CO2
        real    :: Avail,       Mortality,  aux_NC, aux_PC

    !------------------------------------------------------------------------
       
        AM      = Me%PropIndex%Ammonia
        PO      = Me%PropIndex%Phosphate
        DOC     = Me%PropIndex%DOC
        DON     = Me%PropIndex%DON
        DOP     = Me%PropIndex%DOP
        POC     = Me%PropIndex%POC
        PON     = Me%PropIndex%PON
        POP     = Me%PropIndex%POP
        POSi    = Me%PropIndex%BioSilica
        O2      = Me%PropIndex%Oxygen

        If (Me%BioChemPar%CO2_on) then 

            CO2     = Me%PropIndex%CarbonDioxide

        end if
  

 
        Consumer => Me%FirstConsumer
d1:     do while(associated(Consumer))
            
            ZOO_C   = Consumer%PoolIndex%Carbon
            ZOO_N   = Consumer%PoolIndex%Nitrogen
            ZOO_P   = Consumer%PoolIndex%Phosphorus


            !_________.set total predated fraction counter to zero         
            Consumer%Grazing%Total%C = 0.
            Consumer%Grazing%Total%N = 0.
            Consumer%Grazing%Total%P = 0.

            Consumer%Grazing%Spec_Up = 0. 

            !_________________________________________
            !__________________actual element ratios__


            If (Me%ExternalVar%Mass (ZOO_C, Index) .GT. 0.0) then

            Consumer%Ratio%NC_Actual = Me%ExternalVar%Mass (ZOO_N, Index) /      &
                                       Me%ExternalVar%Mass (ZOO_C, Index)           ![units] >>  mmol N / mg C

            Consumer%Ratio%PC_Actual = Me%ExternalVar%Mass (ZOO_P, Index) /      &
                                       Me%ExternalVar%Mass (ZOO_C, Index)           ![units] >>  mmol P / mg C


            else

                Consumer%Ratio%NC_Actual = 0.0

                Consumer%Ratio%PC_Actual = 0.0

            end if

            !_________________________________________
            !__temperature & oxygen limitation calc.__


s1:         Select case (ME%BioChemPar%T_lim_method)

            case (1) s1
          
                Consumer%Limitation%Temperature%Limit = Temp_Lim_Q10 (Consumer%Limitation%Temperature%Q10,  &
                                                        Consumer%Limitation%Temperature%T0, index)

            case (2) s1
         
                Consumer%Limitation%Temperature%Limit = exp(-4000.0 * ((1.0 / (Me%ExternalVar%Temperature (index)            &
                                                        + 273.15)) - (1.0 / (Consumer%Limitation%Temperature%Reftemp + 273.15))))

            end select s1


            Consumer%Limitation%Oxygen = Oxygen_Lim (Consumer%Mortality%O2_Ks, Me%ExternalVar%Mass (O2, Index))

            !_________________________________________
            !________________________predation calc.__                    

            Prey => Consumer%Grazing%FirstPrey
d2:         do while(associated(Prey))

                PreyIndexC  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" carbon"))
                PreyIndexN  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" nitrogen"))
                PreyIndexP  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" phosphorus"))

                if(Prey%Use_Silica)then
                PreyIndexSi = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" silica"))
                end if

                if(Prey%Use_Chl)then
                PreyIndexChl  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" chlorophyll"))
                end if


                !_________.prey actual element ratio                          [units] >>  mmol (N & P) / mg C
                
                If (Me%ExternalVar%Mass (PreyIndexC, Index) .gt. 0.0) then
                
                        Prey%Ratio%NC_Actual = Me%ExternalVar%Mass (PreyIndexN, Index) /      &
                                               Me%ExternalVar%Mass (PreyIndexC, Index)          

                        Prey%Ratio%PC_Actual = Me%ExternalVar%Mass (PreyIndexP, Index) /      &
                                               Me%ExternalVar%Mass (PreyIndexC, Index) 
                
                    else
                    
                        Prey%Ratio%NC_Actual = 0.0
                        Prey%Ratio%PC_Actual = 0.0
                        
                end if
 

                !_________.prey potential grazing
                Avail = Prey%Avail * Me%ExternalVar%Mass (PreyIndexC, Index)            ![units] >>  mg C m-3

                Prey%Pot = Consumer%Grazing%Vmax * (Avail/(Avail+Consumer%Grazing%Ks)) *    &
                           Consumer%Limitation%Temperature%Limit                             ![units] >>  d-1
    

                                                                            ![units] >>  mmol (N & P) m-3 d-1
                Prey%Frac%N = Prey%Pot * Me%ExternalVar%Mass (ZOO_C, Index) * Prey%Ratio%NC_Actual

                Prey%Frac%P = Prey%Pot * Me%ExternalVar%Mass (ZOO_C, Index) * Prey%Ratio%PC_Actual


                Prey%Frac%C = Prey%Pot* Me%ExternalVar%Mass (ZOO_C, Index)          ![units] >>  mg C m-3 d-1


                Consumer%Grazing%Total%C = Consumer%Grazing%Total%C + Prey%Frac%C

                Consumer%Grazing%Total%N = Consumer%Grazing%Total%N + Prey%Frac%N

                Consumer%Grazing%Total%P = Consumer%Grazing%Total%P + Prey%Frac%P


                Consumer%Grazing%Spec_Up = Consumer%Grazing%Spec_Up + Prey%Pot                ![units] >>  d-1



                !___________________________________________________________________________                                                       
                !_________.prey mass balance equations                                            [ matrix update ]

                Me%ExternalVar%Mass (PreyIndexC, Index) = Me%ExternalVar%Mass (PreyIndexC, Index)           &
                                                          - (Prey%Frac%C) * Me%DT_day                   ! <- sinks

                Me%ExternalVar%Mass (PreyIndexN, Index) = Me%ExternalVar%Mass (PreyIndexN, Index)           &
                                                          - (Prey%Frac%N) * Me%DT_day                   ! <- sinks

                Me%ExternalVar%Mass (PreyIndexP, Index) = Me%ExternalVar%Mass (PreyIndexP, Index)           &
                                                          - (Prey%Frac%P) * Me%DT_day                   ! <- sinks              
        


i1:             if(ME%BioChemPar%L_lim_method .eq. 1 .AND. Prey%Use_Chl)then
                                                                                       ![units] >>  mg Chl / mg C
                    Prey%Ratio%ChlC_Actual = Me%ExternalVar%Mass (PreyIndexChl, Index) /      &
                                             Me%ExternalVar%Mass (PreyIndexC, Index) 

                                                                                       ![units] >>  mmol (N & P) m-3 d-1
                    Prey%Frac%Chl = Prey%Pot * Me%ExternalVar%Mass (ZOO_C, Index) * Prey%Ratio%ChlC_Actual


                    Me%ExternalVar%Mass (PreyIndexChl, Index) = Me%ExternalVar%Mass (PreyIndexChl, Index)       &
                                                                - (Prey%Frac%Chl) * Me%DT_day                !<- sinks   

                end if i1

                
                
i2:             if(Prey%Use_Silica)then
                    
                    Prey%Ratio%SiC_Actual = Me%ExternalVar%Mass (PreyIndexSi, Index) /      &
                                            Me%ExternalVar%Mass (PreyIndexC, Index) 

                                                                                  ![units] >>  mmol Si m-3 d-1
                    Prey%Frac%Si = Prey%Pot * Me%ExternalVar%Mass (ZOO_C, Index) * Prey%Ratio%SiC_Actual
                    
                    Me%ExternalVar%Mass (PreyIndexSi, Index) = Me%ExternalVar%Mass (PreyIndexSi, Index)        &
                                                               - Prey%Frac%Si * Me%DT_day             ! <- sinks

                    Me%ExternalVar%Mass (POSi, Index) = Me%ExternalVar%Mass (POSi, Index)        &
                                                        + Prey%Frac%Si * Me%DT_day                    ! <- sources

                end if i2


                Prey => Prey%Next
            end do d2



            !_________________________________________
            !___________assimilated & recycled frac.__ 
                                            
            Consumer%Grazing%NC_Bal = MIN (MAX ((Consumer%Ratio%NC_Max - Consumer%Ratio%NC_Actual) /     &
                                      (Consumer%Ratio%NC_Max - Consumer%Ratio%NC_Min), 0.), 1.)

            Consumer%Grazing%PC_Bal = MIN (MAX ((Consumer%Ratio%PC_Max - Consumer%Ratio%PC_Actual) /     &
                                      (Consumer%Ratio%PC_Max - Consumer%Ratio%PC_Min), 0.), 1.)


            !____________assimilated                                                ![units] >>  mmol N m-3 d-1
            Consumer%Grazing%Assim%N = Consumer%Grazing%Total%N * Consumer%Grazing%NC_Bal

            Consumer%Grazing%Assim%P = Consumer%Grazing%Total%P * Consumer%Grazing%PC_Bal


            !____________recycled (to inorganic fraction)                           ![units] >>  mmol P m-3 d-1
            Consumer%Grazing%Recyc%N = Consumer%Grazing%Total%N * (1. - Consumer%Grazing%NC_Bal)

            Consumer%Grazing%Recyc%P = Consumer%Grazing%Total%P * (1. - Consumer%Grazing%PC_Bal)



            !_________________________________________
            !________________________mortality calc.__ 
                                                                                                       ![units] >>  d-1

            Consumer%Mortality%O2_Dep_r = (1.  - Consumer%Limitation%Oxygen) * Consumer%Mortality%O2_Dep 
                                                                                           

            Mortality = Consumer%Mortality%Rate + Consumer%Mortality%O2_Dep_r


            !_________________________________________
            !________________respiration & excretion__
                                                                                                  ![units] >>  d-1
            Consumer%Respiration%StandStock = Consumer%Limitation%Temperature%Limit * Consumer%Respiration%at10C

            Consumer%Excretion%Rate = Consumer%Grazing%Spec_Up *                                 & 
                                      (1.  - Consumer%Grazing%Ass_Efic) * Consumer%Excretion%Up_Frac

            Consumer%Respiration%Activity = Consumer%Grazing%Spec_Up *                                       & 
                                            (1.  - Consumer%Grazing%Ass_Efic) * (1.  - Consumer%Excretion%Up_Frac)

                                                                               
            Consumer%Respiration%Rate = (Consumer%Respiration%StandStock + Consumer%Respiration%Activity)   &
                                        * Me%ExternalVar%Mass (ZOO_C, Index)            ![units] >>  mg C m-3 d-1



            !________ .balance (mortality) to particulated  
    !        Consumer%Mortality%POM%C = Mortality * Consumer%Mortality%POM_Frac              &
    !                                   * Me%ExternalVar%Mass (ZOO_C, Index)         ![units] >>  mg C m-3 d-1

    !        Consumer%Mortality%POM%N = Mortality * Consumer%Mortality%POM_Frac              &
    !                                   * Me%ExternalVar%Mass (ZOO_N, Index)       ![units] >>  mmol N m-3 d-1
        
    !        Consumer%Mortality%POM%P = Mortality * Consumer%Mortality%POM_Frac              &
    !                                   * Me%ExternalVar%Mass (ZOO_P, Index)       ![units] >>  mmol P m-3 d-1


            !________ .balance (mortality + excretion) to dissolved 
    !        Consumer%Excretion%DOM%C = (Mortality * (1.  - Consumer%Mortality%POM_Frac)         & 
    !                                   + Consumer%Excretion%Rate)                               &
    !                                   * Me%ExternalVar%Mass (ZOO_C, Index)         ![units] >>  mg C m-3 d-1

    !        Consumer%Excretion%DOM%N = (Mortality * (1.  - Consumer%Mortality%POM_Frac)          & 
    !                                   + Consumer%Excretion%Rate)                                &
    !                                   * Me%ExternalVar%Mass (ZOO_N, Index)       ![units] >>  mmol N m-3 d-1

    !        Consumer%Excretion%DOM%P = (Mortality * (1.  - Consumer%Mortality%POM_Frac)          & 
    !                                   + Consumer%Excretion%Rate)                                &
    !                                   * Me%ExternalVar%Mass (ZOO_P, Index)       ![units] >>  mmol P m-3 d-1




    !__________________________________________________________***** Quadratic mortality option *****_________

     !________ .balance (mortality) to particulated  
            Consumer%Mortality%POM%C = ((Mortality * Consumer%Mortality%POM_Frac)           &
                                       + (Mortality * Consumer%Mortality%POM_Frac           &
                                       * (Me%ExternalVar%Mass (ZOO_C, Index)                &
                                       * Me%ExternalVar%Mass (ZOO_C, Index))))              &
                                       * Me%ExternalVar%Mass (ZOO_C, Index)                ![units] >>  mg C m-3 d-1

            Consumer%Mortality%POM%N = ((Mortality * Consumer%Mortality%POM_Frac)           &
                                       + (Mortality * Consumer%Mortality%POM_Frac           &
                                       * (Me%ExternalVar%Mass (ZOO_N, Index)                &
                                       * Me%ExternalVar%Mass (ZOO_N, Index))))              &
                                       * Me%ExternalVar%Mass (ZOO_N, Index)                ![units] >>  mmol N m-3 d-1
        
            Consumer%Mortality%POM%P = ((Mortality * Consumer%Mortality%POM_Frac)           &
                                       + (Mortality * Consumer%Mortality%POM_Frac           &
                                       * (Me%ExternalVar%Mass (ZOO_P, Index)                &
                                       * Me%ExternalVar%Mass (ZOO_P, Index))))              &
                                       * Me%ExternalVar%Mass (ZOO_P, Index)                ![units] >>  mmol P m-3 d-1


            !________ .balance (mortality + excretion) to dissolved 
            Consumer%Excretion%DOM%C = ((Mortality * (1.  - Consumer%Mortality%POM_Frac))       &
                                       + ((Mortality * (1.  - Consumer%Mortality%POM_Frac)      &
                                       * Me%ExternalVar%Mass (ZOO_C, Index)                     &
                                       * Me%ExternalVar%Mass (ZOO_C, Index))))                  &
                                       * Me%ExternalVar%Mass (ZOO_C, Index)                     &
                                       + (Consumer%Excretion%Rate                               &
                                       * Me%ExternalVar%Mass (ZOO_C, Index))        ![units] >>  mg C m-3 d-1

            Consumer%Excretion%DOM%N = ((Mortality * (1.  - Consumer%Mortality%POM_Frac))       &
                                       + ((Mortality * (1.  - Consumer%Mortality%POM_Frac)      &
                                       * Me%ExternalVar%Mass (ZOO_N, Index)                     &
                                       * Me%ExternalVar%Mass (ZOO_N, Index))))                  &
                                       * Me%ExternalVar%Mass (ZOO_N, Index)                     &
                                       + (Consumer%Excretion%Rate                               &
                                       * Me%ExternalVar%Mass (ZOO_N, Index))        ![units] >>  mmol N m-3 d-1

            Consumer%Excretion%DOM%P = ((Mortality * (1.  - Consumer%Mortality%POM_Frac))       &
                                       + ((Mortality * (1.  - Consumer%Mortality%POM_Frac)      &
                                       * Me%ExternalVar%Mass (ZOO_P, Index)                     &
                                       * Me%ExternalVar%Mass (ZOO_P, Index))))                  &
                                       * Me%ExternalVar%Mass (ZOO_P, Index)                     &
                                       + (Consumer%Excretion%Rate                               &
                                       * Me%ExternalVar%Mass (ZOO_P, Index))        ![units] >>  mmol P m-3 d-1


    !___________________________________________________________________________________________________________



    !___________________________________________________________
    !___________________________________mass balance equations__                            [ matrix update ]

            !_________.consumers [c n p]
            Me%ExternalVar%Mass (ZOO_C, Index) = Me%ExternalVar%Mass (ZOO_C, Index)         &
                                                 + (Consumer%Grazing%Total%C                &   ! <- sources
                                                 - (Consumer%Mortality%POM%C                &   ! <- sinks
                                                 + Consumer%Excretion%DOM%C                 &
                                                 + Consumer%Respiration%Rate)) * Me%DT_day       

            Me%ExternalVar%Mass (ZOO_N, Index) = Me%ExternalVar%Mass (ZOO_N, Index)         &
                                                 + (Consumer%Grazing%Assim%N                &   ! <- sources
                                                 - (Consumer%Mortality%POM%N                &   ! <- sinks
                                                 + Consumer%Excretion%DOM%N)) * Me%DT_day
 
            Me%ExternalVar%Mass (ZOO_P, Index) = Me%ExternalVar%Mass (ZOO_P, Index)         &
                                                 + (Consumer%Grazing%Assim%P                &   ! <- sources
                                                 - (Consumer%Mortality%POM%P                &   ! <- sinks
                                                 + Consumer%Excretion%DOM%P)) * Me%DT_day


            !_________.inorganic nutrients [n p]

            Me%ExternalVar%Mass (AM, Index) = Me%ExternalVar%Mass (AM, Index)           &
                                              + Consumer%Grazing%Recyc%N * Me%DT_day

            Me%ExternalVar%Mass (PO, Index) = Me%ExternalVar%Mass (PO, Index)           &
                                              + Consumer%Grazing%Recyc%P * Me%DT_day


            !_________.diss organic matter [c n p]
            Me%ExternalVar%Mass (DOC, Index) = Me%ExternalVar%Mass (DOC, Index)   &
                                               + Consumer%Excretion%DOM%C * Me%DT_day       ! <- sources

            Me%ExternalVar%Mass (DON, Index) = Me%ExternalVar%Mass (DON, Index)   &
                                               + Consumer%Excretion%DOM%N * Me%DT_day       ! <- sources

            Me%ExternalVar%Mass (DOP, Index) = Me%ExternalVar%Mass (DOP, Index)   &
                                               + Consumer%Excretion%DOM%P * Me%DT_day       ! <- sources


 
            !_________.part organic matter [c n p]
            Me%ExternalVar%Mass (POC, Index) = Me%ExternalVar%Mass (POC, Index)   &
                                               + Consumer%Mortality%POM%C * Me%DT_day       ! <- sources

            Me%ExternalVar%Mass (PON, Index) = Me%ExternalVar%Mass (PON, Index)   &
                                               + Consumer%Mortality%POM%N * Me%DT_day       ! <- sources

            Me%ExternalVar%Mass (POP, Index) = Me%ExternalVar%Mass (POP, Index)   &
                                               + Consumer%Mortality%POM%P * Me%DT_day       ! <- sources



        !_________.oxygen [O]                                                  [units] >>  mg O2 L-1 d-1

            Me%ExternalVar%Mass (O2, Index) = Me%ExternalVar%Mass (O2, Index)             &
                                              - Consumer%Respiration%Rate                 &
                                              * Me%BioChemPar%O2C_Conversion * Me%DT_day


        !_________.carbon dioxide [CO2]                                        [units] >>  mg CO2 L-1 d-1

         If (Me%BioChemPar%CO2_on) then 

            Me%ExternalVar%Mass (CO2, Index) = Me%ExternalVar%Mass (CO2, Index)           &
                                              + Consumer%Respiration%Rate                 &
                                              * Me%BioChemPar%CO2C_Conversion * Me%DT_day
         Endif


            !_________________________________________
            !______________nutrient excess excretion  


            !_________________________________________
            !__________________updated element ratios__

            aux_NC = Me%ExternalVar%Mass (ZOO_N, Index) / Me%ExternalVar%Mass (ZOO_C, Index)     ![units] >>  mmol N / mg C

            aux_PC = Me%ExternalVar%Mass (ZOO_P, Index) / Me%ExternalVar%Mass (ZOO_C, Index)     ![units] >>  mmol P / mg C
    

            !_______.nitrogen                                                           [units] >>  mmol N m-3
            if (aux_NC .gt. Consumer%Ratio%NC_Max) then

                Consumer%Excretion%Rel_Exc%N = (aux_NC - Consumer%Ratio%NC_Max) *    &
                                                Me%ExternalVar%Mass (ZOO_C, Index) 
            
            else
                Consumer%Excretion%Rel_Exc%N = 0.
            end if


            !_______.phosphorus                                                         [units] >>  mmol P m-3
            if (aux_PC .gt. Consumer%Ratio%PC_Max) then

                Consumer%Excretion%Rel_Exc%P = (aux_PC - Consumer%Ratio%PC_Max) *    &
                                                Me%ExternalVar%Mass (ZOO_C, Index) 
            
            else
                Consumer%Excretion%Rel_Exc%P = 0.
            end if

            !___________________________________________________________
            !___________________________update mass balance equations __                            [ matrix update ]


            Me%ExternalVar%Mass (ZOO_N, Index) = Me%ExternalVar%Mass (ZOO_N, Index)     &
                                                 - Consumer%Excretion%Rel_Exc%N


            Me%ExternalVar%Mass (ZOO_P, Index) = Me%ExternalVar%Mass (ZOO_P, Index)     &
                                                 - Consumer%Excretion%Rel_Exc%P



            !_________.inorganic nutrients [n p]

            Me%ExternalVar%Mass (AM, Index) = Me%ExternalVar%Mass (AM, Index)           &
                                              + Consumer%Excretion%Rel_Exc%N            ! <- sources  
                        

            Me%ExternalVar%Mass (PO, Index) = Me%ExternalVar%Mass (PO, Index)           &
                                              + Consumer%Excretion%Rel_Exc%P            ! <- sources  
    


            Consumer => Consumer%Next
        
        end do d1

    end subroutine Consumers
 !--------------------------------------------------------------------------



!_________________________________________
!_______________________Decomposers calc__

    subroutine Decomposers (index)
    
    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        type(T_Decomposer),    pointer             :: Decomposer

        integer :: BAC_C,           BAC_N,          BAC_P
        integer :: AM,              NA,             PO   
        integer :: DOC,             DOCsl,          POC 
        integer :: DON,             PON 
        integer :: DOP,             POP
        integer :: O2,              CO2
        real    :: Mortality,       POM_Frac
        real    :: DOM_NC_RATIO,    DOM_PC_RATIO
        real    :: InN_Lim,         InP_Lim
        real    :: aux_NC,          aux_PC
    !------------------------------------------------------------------------
       
        AM      = Me%PropIndex%Ammonia
        NA      = Me%PropIndex%Nitrate
        PO      = Me%PropIndex%Phosphate
        DOC     = Me%PropIndex%DOC
        DON     = Me%PropIndex%DON
        DOP     = Me%PropIndex%DOP
        DOCsl   = Me%PropIndex%DOC_SL
        POC     = Me%PropIndex%POC
        PON     = Me%PropIndex%PON
        POP     = Me%PropIndex%POP
        O2      = Me%PropIndex%Oxygen

        If (Me%BioChemPar%CO2_on) then 

            CO2     = Me%PropIndex%CarbonDioxide

        end if

 
        Decomposer => Me%FirstDecomposer
d1:     do while(associated(Decomposer))
            
            BAC_C   = Decomposer%PoolIndex%Carbon
            BAC_N   = Decomposer%PoolIndex%Nitrogen
            BAC_P   = Decomposer%PoolIndex%Phosphorus


            !_________________________________________
            !__________________actual element ratios__

            Decomposer%Ratio%NC_Actual = Me%ExternalVar%Mass (BAC_N, Index) /      &
                                       Me%ExternalVar%Mass (BAC_C, Index)           ![units] >>  mmol N / mg C

            Decomposer%Ratio%PC_Actual = Me%ExternalVar%Mass (BAC_P, Index) /      &
                                       Me%ExternalVar%Mass (BAC_C, Index)           ![units] >>  mmol P / mg C


            DOM_NC_RATIO  =   Me%ExternalVar%Mass (DON, Index) /      &
                            Me%ExternalVar%Mass (DOC, Index)                        ![units] >>  mmol N / mg C

            DOM_PC_RATIO  =   Me%ExternalVar%Mass (DOP, Index) /      &
                            Me%ExternalVar%Mass (DOC, Index)                        ![units] >>  mmol P / mg C

            !_________________________________________
            !__temperature & oxygen limitation calc.__


s1:         select case (ME%BioChemPar%T_lim_method)

            case (1) s1

                Decomposer%Limitation%Temperature%Limit = Temp_Lim_Q10 (Decomposer%Limitation%Temperature%Q10,        &
                                                          Decomposer%Limitation%Temperature%T0, index)

            case (2) s1
     
                Decomposer%Limitation%Temperature%Limit = exp(-4000.0 * ((1.0 / (Me%ExternalVar%Temperature (index)            &
                                         + 273.15)) - (1.0 / (Decomposer%Limitation%Temperature%Reftemp + 273.15))))

            end select s1


            Decomposer%Limitation%Oxygen = Oxygen_Lim (Decomposer%Mortality%O2_Ks, Me%ExternalVar%Mass (O2, Index))


            !__________For Biochemical calc use
            ME%BioChemPar%Bac_Temp_Lim = Decomposer%Limitation%Temperature%Limit



            !_________________________________________
            !______________nutriente limitation calc__

            Decomposer%Limitation%Nutrients%N_Status = Nutrient_Lim (Decomposer%Ratio%NC_Actual,            & 
                                                       Decomposer%Ratio%NC_Min, Decomposer%Ratio%NC_Max)
    
            Decomposer%Limitation%Nutrients%P_Status = Nutrient_Lim (Decomposer%Ratio%PC_Actual,            & 
                                                       Decomposer%Ratio%PC_Min, Decomposer%Ratio%PC_Max)


            Decomposer%Limitation%Nutrients%Int_Nut_Limit = MIN (Decomposer%Limitation%Nutrients%N_Status,  &
                                                            Decomposer%Limitation%Nutrients%P_Status)
            !_________________________________________
            !_____________________________DOM uptake__
                                                                                         ![units] >>  d-1
            Decomposer%Consume%DOM_Up_Pot = Decomposer%Consume%Vmax *                   & 
                                            Decomposer%Limitation%Temperature%Limit *   &
                                            (Me%ExternalVar%Mass (DOC, Index) /         &
                                            (Me%ExternalVar%Mass (DOC, Index) + Decomposer%Consume%Ks))

                                                                                    ![units] >>  mg C m-3 d-1
            Decomposer%Consume%DOM_Up_Frac%C = MAX( MIN((Decomposer%Limitation%Nutrients%Int_Nut_Limit *    &
                                               Decomposer%Consume%DOM_Up_Pot *                              &
                                               Me%ExternalVar%Mass (BAC_C, Index)),                         &
                                               Me%ExternalVar%Mass (DOC, Index)/ Me%DT_day), 0.0)

                                                                          ![units] >>  mmol (N or P)  m-3 d-1
            Decomposer%Consume%DOM_Up_Frac%N = MIN ((Decomposer%Consume%DOM_Up_Frac%C * DOM_NC_RATIO),   &
                                               Me%ExternalVar%Mass (DON, Index)/ Me%DT_day) 

            Decomposer%Consume%DOM_Up_Frac%P = MIN ((Decomposer%Consume%DOM_Up_Frac%C * DOM_PC_RATIO),   &
                                               Me%ExternalVar%Mass (DOP, Index)/ Me%DT_day) 

            !_________________________________________
            !________________inorg. nutrients uptake__
                                                                                  ![units] >>  mmol (N or P) m-3 d-1

            InN_Lim = MAX((Decomposer%Ratio%NC_max - Decomposer%Ratio%NC_Actual) /       &
                      (Decomposer%Ratio%NC_max - Decomposer%Ratio%NC_min), 0.0)


            Decomposer%Consume%Up_NH4 = Min(Decomposer%Consume%NH4_Ks * Me%ExternalVar%Mass (AM, Index)    &
                                        * Me%ExternalVar%Mass (BAC_C, Index) * InN_Lim,                    &
                                        Me%ExternalVar%Mass (AM, Index)/ Me%DT_day)


            Decomposer%Consume%Up_NO3 = Min(Decomposer%Consume%NO3_Ks * Me%ExternalVar%Mass (NA, Index)  &
                                        * Me%ExternalVar%Mass (BAC_C, Index) * InN_Lim,                  &
                                        Me%ExternalVar%Mass (NA, Index)/ Me%DT_day)



            InP_Lim = MAX((Decomposer%Ratio%PC_max - Decomposer%Ratio%PC_Actual) /       &
                      (Decomposer%Ratio%PC_max - Decomposer%Ratio%PC_min), 0.0)


            Decomposer%Consume%Up_PO4 = Min(Decomposer%Consume%PO4_Ks * Me%ExternalVar%Mass (PO, Index)  &
                                        * Me%ExternalVar%Mass (BAC_C, Index) * InP_Lim,                  &
                                        Me%ExternalVar%Mass (PO, Index)/ Me%DT_day)



            !_________________________________________
            !______________________respiration calc.__ 

            if (Me%ExternalVar%Mass (O2, Index) .gt. Decomposer%Consume%O2_Low_Ass) then
         
                Decomposer%Consume%Ass_Efic = Decomposer%Consume%Ass_Efic_Norm

            else
                Decomposer%Consume%Ass_Efic = Decomposer%Consume%Ass_Efic_LowO2
            end if


            Decomposer%Respiration%Activity = (1.  - Decomposer%Consume%Ass_Efic) *         &
                                              Decomposer%Consume%DOM_Up_Frac%C *            &
                                              (1.  - Decomposer%Limitation%Oxygen)    ![units] >>  mg C m-3 d-1

            Decomposer%Respiration%StandStock = Decomposer%Respiration%at10C *              &
                                                Decomposer%Limitation%Temperature%Limit     &
                                                * Me%ExternalVar%Mass (BAC_C, Index)  ![units] >>  mg C m-3 d-1

            Decomposer%Respiration%Rate = Decomposer%Respiration%Activity + Decomposer%Respiration%StandStock
                                                                                      ![units] >>  mg C m-3 d-1
            !_________________________________________
            !________________________mortality calc.__ 


            Decomposer%Mortality%Lysis = Decomposer%Mortality%Density_Dep   *                   &
                                         (Me%ExternalVar%Mass (BAC_C, Index) /                  &
                                         Decomposer%Mortality%Lysis_Ref_Con)

            Mortality = Decomposer%Mortality%Rate + Decomposer%Mortality%Lysis       ![units] >>  d-1



            !POM_Frac calc.
            if (Mortality .gt. 0.) then

                POM_Frac = Min (1.  - (Decomposer%Mortality%Lysis / Mortality), Decomposer%Mortality%POM_Frac)

            else
                POM_Frac = Decomposer%Mortality%POM_Frac
            end if


            !________ .balance (mortality) to particulated  
            Decomposer%Mortality%POM%C = Mortality * POM_Frac              &
                                       * Me%ExternalVar%Mass (BAC_C, Index)        ![units] >>  mg C m-3 d-1

            Decomposer%Mortality%POM%N = Mortality * POM_Frac              &
                                       * Me%ExternalVar%Mass (BAC_N, Index)        ![units] >>  mmol N m-3 d-1
        
            Decomposer%Mortality%POM%P = Mortality * POM_Frac              &
                                       * Me%ExternalVar%Mass (BAC_P, Index)        ![units] >>  mmol P m-3 d-1


            !________ .balance (mortality + excretion) to dissolved 
        
            Decomposer%Mortality%DOM_Tot%C = Mortality * (1.  - POM_Frac)               & 
                                             * Me%ExternalVar%Mass (BAC_C, Index)     ![units] >>  mg C m-3 d-1
        
            Decomposer%Mortality%DOM_L%C = Decomposer%Mortality%DOM_Tot%C           &
                                         * (1.  - Decomposer%Mortality%DOC_SL_Frac)    

            Decomposer%Mortality%DOM_SL%C = Decomposer%Mortality%DOM_Tot%C          &
                                            * Decomposer%Mortality%DOC_SL_Frac   


            Decomposer%Mortality%DOM_L%N = Mortality * (1.  - POM_Frac)       & 
                                       * Me%ExternalVar%Mass (BAC_N, Index)        ![units] >>  mmol N m-3 d-1

            Decomposer%Mortality%DOM_L%P = Mortality * (1.  - POM_Frac)       & 
                                       * Me%ExternalVar%Mass (BAC_P, Index)        ![units] >>  mmol P m-3 d-1

    !___________________________________________________________
    !___________________________________mass balance equations__                            [ matrix update ]

            !_________.decomposers [c n p]
            Me%ExternalVar%Mass (BAC_C, Index) = Me%ExternalVar%Mass (BAC_C, Index)         &
                                                 + (Decomposer%Consume%DOM_Up_Frac%C        &   ! <- sources    
                                                 - (Decomposer%Mortality%DOM_Tot%C          &   ! <- sinks
                                                 + Decomposer%Respiration%Rate              &
                                                 + Decomposer%Mortality%POM%C)) * Me%DT_day

            Me%ExternalVar%Mass (BAC_N, Index) = Me%ExternalVar%Mass (BAC_N, Index)         &
                                                 + (Decomposer%Consume%DOM_Up_Frac%N        &   ! <- sources
                                                 + Decomposer%Consume%Up_NH4                &
                                                 + Decomposer%Consume%Up_NO3                &
                                                 - (Decomposer%Mortality%DOM_L%N +          &   ! <- sinks
                                                 Decomposer%Mortality%POM%N)) * Me%DT_day 

            Me%ExternalVar%Mass (BAC_P, Index) = Me%ExternalVar%Mass (BAC_P, Index)         &
                                                 + (Decomposer%Consume%DOM_Up_Frac%P        &   ! <- sources
                                                 + Decomposer%Consume%Up_PO4                &
                                                 - (Decomposer%Mortality%DOM_L%P +          &   ! <- sinks
                                                 Decomposer%Mortality%POM%P)) * Me%DT_day                                          



            !_________.inorganic nutrients [n p]

            Me%ExternalVar%Mass (AM, Index) = Me%ExternalVar%Mass (AM, Index)           &
                                              - Decomposer%Consume%Up_NH4 * Me%DT_day       ! <- sinks
                                        
            Me%ExternalVar%Mass (NA, Index) = Me%ExternalVar%Mass (NA, Index)           &
                                              - Decomposer%Consume%Up_NO3 * Me%DT_day       ! <- sinks

            Me%ExternalVar%Mass (PO, Index) = Me%ExternalVar%Mass (PO, Index)           &
                                              - Decomposer%Consume%Up_PO4 * Me%DT_day       ! <- sinks
                                        

  
            !_________.diss organic matter [c n p]
            Me%ExternalVar%Mass (DOC, Index) = Me%ExternalVar%Mass (DOC, Index)     &
                                               + (Decomposer%Mortality%DOM_L%C      &   ! <- sources
                                               - Decomposer%Consume%DOM_Up_Frac%C)  &   ! <- sinks
                                               * Me%DT_day
                                            

            Me%ExternalVar%Mass (DON, Index) = Me%ExternalVar%Mass (DON, Index)     &
                                               + (Decomposer%Mortality%DOM_L%N      &   ! <- sources
                                               - Decomposer%Consume%DOM_Up_Frac%N)  &   ! <- sinks
                                               * Me%DT_day
                                            

            Me%ExternalVar%Mass (DOP, Index) = Me%ExternalVar%Mass (DOP, Index)     &
                                               + (Decomposer%Mortality%DOM_L%P      &   ! <- sources
                                               - Decomposer%Consume%DOM_Up_Frac%P)  &   ! <- sinks
                                               * Me%DT_day


            Me%ExternalVar%Mass (DOCsl, Index) = Me%ExternalVar%Mass (DOCsl, Index)     &
                                                 + Decomposer%Mortality%DOM_SL%C * Me%DT_day



            !_________.part organic matter [c n p]
            Me%ExternalVar%Mass (POC, Index) = Me%ExternalVar%Mass (POC, Index)   &
                                               + Decomposer%Mortality%POM%C * Me%DT_day

            Me%ExternalVar%Mass (PON, Index) = Me%ExternalVar%Mass (PON, Index)   &
                                               + Decomposer%Mortality%POM%N * Me%DT_day

            Me%ExternalVar%Mass (POP, Index) = Me%ExternalVar%Mass (POP, Index)   &
                                               + Decomposer%Mortality%POM%P * Me%DT_day



            !_________.oxygen [O]                                       [units] >>  mg O2 L-1 d-1

            Me%ExternalVar%Mass (O2, Index) = Me%ExternalVar%Mass (O2, Index)               &
                                              - Decomposer%Respiration%Rate                 &
                                              * Me%BioChemPar%O2C_Conversion * Me%DT_day



            !_________.carbon dioxide [CO2]                             [units] >>  mg CO2 L-1 d-1

             If (Me%BioChemPar%CO2_on) then 

                Me%ExternalVar%Mass (CO2, Index) = Me%ExternalVar%Mass (CO2, Index)           &
                                                  + Decomposer%Respiration%Rate               &
                                                  * Me%BioChemPar%CO2C_Conversion * Me%DT_day
             Endif



            !__________For Biochemical calc use
            ME%BioChemPar%Bac_conc = Me%ExternalVar%Mass (BAC_C, Index)

            !_________________________________________
            !______________nutrient excess exudation__  


            !_________________________________________
            !__________________updated element ratios__

            aux_NC = Me%ExternalVar%Mass (BAC_N, Index) / Me%ExternalVar%Mass (BAC_C, Index)     ![units] >>  mmol N / mg C

            aux_PC = Me%ExternalVar%Mass (BAC_P, Index) / Me%ExternalVar%Mass (BAC_C, Index)     ![units] >>  mmol P / mg C
    

            !_________________________________________
            !______________nutrient excess excretion__  

            !_______.nitrogen                                                           [units] >>  mmol N m-3
            if (aux_NC .gt. Decomposer%Ratio%NC_Max) then

                Decomposer%Excretion%Rel_Exc%N = (aux_NC - Decomposer%Ratio%NC_Max) *    &
                                                 Me%ExternalVar%Mass (BAC_C, Index) 
            
            else
                Decomposer%Excretion%Rel_Exc%N = 0.
            end if


            !_______.phosphorus                                                         [units] >>  mmol P m-3
            if (aux_PC .gt. Decomposer%Ratio%PC_Max) then

                Decomposer%Excretion%Rel_Exc%P = (aux_PC - Decomposer%Ratio%PC_Max) *    &
                                                  Me%ExternalVar%Mass (BAC_C, Index) 
        
            else
                Decomposer%Excretion%Rel_Exc%P = 0.
            end if
     

            !___________________________________________________________
            !___________________________update mass balance equations __                            [ matrix update ]


            Me%ExternalVar%Mass (BAC_N, Index) = Me%ExternalVar%Mass (BAC_N, Index)     &
                                                 - Decomposer%Excretion%Rel_Exc%N


            Me%ExternalVar%Mass (BAC_P, Index) = Me%ExternalVar%Mass (BAC_P, Index)     &
                                                 - Decomposer%Excretion%Rel_Exc%P



            !_________.inorganic nutrients [n p]

            Me%ExternalVar%Mass (AM, Index) = Me%ExternalVar%Mass (AM, Index)           &
                                              + Decomposer%Excretion%Rel_Exc%N            ! <- sources  
                        

            Me%ExternalVar%Mass (PO, Index) = Me%ExternalVar%Mass (PO, Index)           &
                                              + Decomposer%Excretion%Rel_Exc%P            ! <- sources  
                        


            Decomposer => Decomposer%Next
        end do d1

    end subroutine Decomposers
 !--------------------------------------------------------------------------



!_________________________________________
!_______________________BioChemical calc__

    subroutine BioChemicalProcesses (index)
    
    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
     
        integer :: Si,              POSi,       AM,       NA
        integer :: DOCsl,           DONsl,      DOPsl,    O2
        integer :: POC,             PON,        POP
        integer :: DOC,             DON,        DOP
        real    :: POM_NC_RATIO,    POM_PC_RATIO
        real    :: DOMsl_NC_RATIO,  DOMsl_PC_RATIO
        real    :: BioSi_Dissol
        real    :: AverageRadiation

    !------------------------------------------------------------------------
       
        Si      = Me%PropIndex%Silicate
        POSi    = Me%PropIndex%BioSilica
        DOCsl   = Me%PropIndex%DOC_SL
        DONsl   = Me%PropIndex%DON_SL
        DOPsl   = Me%PropIndex%DOP_SL
        POC     = Me%PropIndex%POC
        PON     = Me%PropIndex%PON
        POP     = Me%PropIndex%POP
        DOC     = Me%PropIndex%DOC
        DON     = Me%PropIndex%DON
        DOP     = Me%PropIndex%DOP
        AM      = Me%PropIndex%Ammonia
        NA      = Me%PropIndex%Nitrate
        O2      = Me%PropIndex%Oxygen


        AverageRadiation = Me%ExternalVar%ShortWaveAverage(index) 
!_________________________________________
!__________________BioSilica dissolution__

                                                                           ![units] >>  mmol Si m-3 d-1
        BioSi_Dissol = Me%ExternalVar%Mass (POSi, Index) * Me%BioChemPar%BioSi_Diss
        

!_________________________________________
!________bacteria mediated OM hydrolysis__

       
        !__________organic matter element ratios__   
        DOMsl_NC_RATIO  =   Me%ExternalVar%Mass (DONsl, Index) /      &
                          Me%ExternalVar%Mass (DOCsl, Index)                  ![units] >>  mmol N / mg C

        DOMsl_PC_RATIO  =   Me%ExternalVar%Mass (DOPsl, Index) /      &
                          Me%ExternalVar%Mass (DOCsl, Index)                  ![units] >>  mmol P / mg C


        POM_NC_RATIO  =   Me%ExternalVar%Mass (PON, Index) /      &
                        Me%ExternalVar%Mass (POC, Index)                     ![units] >>  mmol N / mg C

        POM_PC_RATIO  =   Me%ExternalVar%Mass (POP, Index) /      &
                        Me%ExternalVar%Mass (POC, Index)                     ![units] >>  mmol P / mg C


        !__________POM to DOM calc__                                              
                                                                               ![units] >>  mg C m-3 d-1
        Me%BioChemPar%POM_bac_Hyd%C = MAX (MIN (Me%BioChemPar%POM_bac_vmax *          &
                                    (Me%ExternalVar%Mass (POC, Index) /             &
                                    (Me%ExternalVar%Mass (POC, Index) +             &
                                    Me%BioChemPar%POM_bac_ks)) *                    &
                                    ME%BioChemPar%Bac_Temp_Lim *                    &
                                    ME%BioChemPar%Bac_conc,                         &
                                    Me%ExternalVar%Mass (POC, Index) /              &
                                    Me%DT_day), 0.0)

                                                                              ![units] >>  mmol N m-3 d-1
        Me%BioChemPar%POM_bac_Hyd%N = Me%BioChemPar%POM_bac_Hyd%C * POM_NC_RATIO  
                    
                                                                              ![units] >>  mmol P m-3 d-1
        Me%BioChemPar%POM_bac_Hyd%P = Me%BioChemPar%POM_bac_Hyd%C * POM_PC_RATIO



        !__________DOMsl to DOM calc__  
                                                                              ![units] >>  mg C m-3 d-1
        Me%BioChemPar%DOMsl_bac_Hyd%C = MAX (MIN (Me%BioChemPar%DOMsl_bac_vmax *           &
                                     (Me%ExternalVar%Mass (DOCsl, Index) /              &
                                     (Me%ExternalVar%Mass (DOCsl, Index) +              &
                                     Me%BioChemPar%DOMsl_bac_ks)) *                     &
                                     ME%BioChemPar%Bac_Temp_Lim *                       &
                                     ME%BioChemPar%Bac_conc,                            &
                                     Me%ExternalVar%Mass (DOCsl, Index)/                &
                                     Me%DT_day), 0.0)

                                                                              ![units] >>  mmol N m-3 d-1
        Me%BioChemPar%DOMsl_bac_Hyd%N = Me%BioChemPar%DOMsl_bac_Hyd%C * DOMsl_NC_RATIO 

                                                                              ![units] >>  mmol P m-3 d-1
        Me%BioChemPar%DOMsl_bac_Hyd%P = Me%BioChemPar%DOMsl_bac_Hyd%C * DOMsl_PC_RATIO 



        !__________nitrification
                                                                    
        if (AverageRadiation .le. Me%BioChemPar%Nitrifradiation) then 

            Me%BioChemPar%Nitrif_lim = (1.0 - EXP(-ME%BioChemPar%Nit_Inib_coef * Me%ExternalVar%Mass (O2, Index)))

                                                                                                ![units] >>  mmol N m-3 d-1
            Me%BioChemPar%Nitrif = Me%BioChemPar%Nitrifrate * Me%ExternalVar%Mass (AM, Index) * Me%BioChemPar%Nitrif_lim
    
        else

            Me%BioChemPar%Nitrif = 0.0
        endif



        !___________________________________________________________
        !___________________________________mass balance equations__                            [ matrix update ]


        !_________.nitrogen species [NA, AM]
        Me%ExternalVar%Mass (NA, Index) = Me%ExternalVar%Mass (NA, Index)              &
                                          + Me%BioChemPar%Nitrif * Me%DT_day                     ! <- sources

        Me%ExternalVar%Mass (AM, Index) = Me%ExternalVar%Mass (AM, Index)              &
                                          - Me%BioChemPar%Nitrif * Me%DT_day                     ! <- sinks

        !_________.oxygen
        Me%ExternalVar%Mass (O2, Index) = Me%ExternalVar%Mass (O2, Index)              &
                                          - (Me%BioChemPar%Nitrif                      &         ! <- sinks
                                          * ME%BioChemPar%Nit_ON_Conv) * Me%DT_day


                                          
            !_________.carbon dioxide [CO2]                             [units] >>  mg CO2 L-1 d-1

!             If (Me%BioChemPar%CO2_on) then 

!                Me%ExternalVar%Mass (CO2, Index) = Me%ExternalVar%Mass (CO2, Index)           &
!                                                  + Decomposer%Respiration%Rate               &
!                                                  * Me%BioChemPar%CO2C_Conversion * Me%DT_day
!             Endif




        !_________.silica species [Si, BSi]
        Me%ExternalVar%Mass (POSi, Index) = Me%ExternalVar%Mass (POSi, Index)          &
                                             - BioSi_Dissol * Me%DT_day                     ! <- sinks

        Me%ExternalVar%Mass (Si, Index) = Me%ExternalVar%Mass (Si, Index)              &
                                          + BioSi_Dissol * Me%DT_day                        ! <- sources


        !_________.semi-labile diss organic matter [c n p]
        Me%ExternalVar%Mass (DOCsl, Index) = Me%ExternalVar%Mass (DOCsl, Index)     &
                                             - Me%BioChemPar%DOMsl_bac_Hyd%C * Me%DT_day

        Me%ExternalVar%Mass (DONsl, Index) = Me%ExternalVar%Mass (DONsl, Index)     &
                                             - Me%BioChemPar%DOMsl_bac_Hyd%N * Me%DT_day

        Me%ExternalVar%Mass (DOPsl, Index) = Me%ExternalVar%Mass (DOPsl, Index)     &
                                             - Me%BioChemPar%DOMsl_bac_Hyd%P * Me%DT_day


        !_________.part organic matter [c n p]
        Me%ExternalVar%Mass (POC, Index) = Me%ExternalVar%Mass (POC, Index)   &
                                           - Me%BioChemPar%POM_bac_Hyd%C * Me%DT_day

        Me%ExternalVar%Mass (PON, Index) = Me%ExternalVar%Mass (PON, Index)   &
                                           - Me%BioChemPar%POM_bac_Hyd%N * Me%DT_day

        Me%ExternalVar%Mass (POP, Index) = Me%ExternalVar%Mass (POP, Index)   &
                                           - Me%BioChemPar%POM_bac_Hyd%P * Me%DT_day


        !_________.diss organic matter [c n p]
        Me%ExternalVar%Mass (DOC, Index) = Me%ExternalVar%Mass (DOC, Index)     &
                                           + (Me%BioChemPar%DOMsl_bac_Hyd%C     &   ! <- sources
                                           + Me%BioChemPar%POM_bac_Hyd%C)       &   
                                           * Me%DT_day
                                    

        Me%ExternalVar%Mass (DON, Index) = Me%ExternalVar%Mass (DON, Index)     &
                                           + (Me%BioChemPar%DOMsl_bac_Hyd%N     &   ! <- sources
                                           + Me%BioChemPar%POM_bac_Hyd%N)       &   
                                           * Me%DT_day
                                    

        Me%ExternalVar%Mass (DOP, Index) = Me%ExternalVar%Mass (DOP, Index)     &
                                           + (Me%BioChemPar%DOMsl_bac_Hyd%P     &   ! <- sources
                                           + Me%BioChemPar%POM_bac_Hyd%P)       &   
                                           * Me%DT_day


    end subroutine BioChemicalProcesses
    !--------------------------------------------------------------------------






!_________________________________________
!__________temperature limitation factor__

    function Temp_Lim_Q10 (Q10, T0, index)

                real                :: Temp_Lim_Q10
    
    !Arguments---------------------------------------------------------------

        real                :: Q10
        real                :: T0
        integer, intent(IN) :: index
    !------------------------------------------------------------------------
    
        Temp_Lim_Q10 = Q10**((Me%ExternalVar%Temperature (index) - T0) / T0)
 
    end function Temp_Lim_Q10
   
    !--------------------------------------------------------------------------


!_________________________________________
!_______________oxygen limitation factor__

    function Oxygen_Lim (Ks, O2_Conc)

                real                :: Oxygen_Lim
    
    !Arguments---------------------------------------------------------------

        real                :: Ks
        real                :: O2_Conc
    !------------------------------------------------------------------------
    
        Oxygen_Lim = O2_Conc / (O2_Conc + Ks)
 
    end function Oxygen_Lim
   
    !--------------------------------------------------------------------------


!_________________________________________
!_____________nutrient limitation factor__

    function Nutrient_Lim (Ratio_Actual, Ratio_Ref_1, Ratio_Ref_2)

                real                :: Nutrient_Lim
    
    !Arguments---------------------------------------------------------------

        real                :: Ratio_Actual
        real                :: Ratio_Ref_1                     !Minimum ratio
        real                :: Ratio_Ref_2                     !Redfield or maximum ratio
    !------------------------------------------------------------------------
    
        Nutrient_Lim = MIN(1., MAX((Ratio_Actual - Ratio_Ref_1)/(Ratio_Ref_2 - Ratio_Ref_1), 0.))
 
    end function Nutrient_Lim
   
    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillLife(ObjLifeID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjLifeID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_, nUsers             

        !Local-------------------------------------------------------------------
        integer                             :: STAT_            

        !------------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(ObjLifeID, ready_)    

if1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mLIFE_,  Me%InstanceID)
  
            if (nUsers == 0) then

                call DeallocateInstance
                ObjLifeID = 0

                STAT_ = SUCCESS_

            endif

        else 

            STAT_ = ready_

        end if if1

        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine KillLife
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Life), pointer          :: AuxObjLife
        type (T_Life), pointer          :: PreviousObjLife

        !Updates pointers
        if (Me%InstanceID == FirstObjLife%InstanceID) then
            FirstObjLife => FirstObjLife%Next
        else
            PreviousObjLife => FirstObjLife
            AuxObjLife      => FirstObjLife%Next
            do while (AuxObjLife%InstanceID /= Me%InstanceID)
                PreviousObjLife => AuxObjLife
                AuxObjLife      => AuxObjLife%Next
            enddo

            !Now update linked list
            PreviousObjLife%Next => AuxObjLife%Next

            !Deallocates instance
            deallocate (Me)
            nullify    (Me) 

        endif
            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjLifeID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjLifeID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjLifeID > 0) then
            call LocateObjLife (ObjLifeID)
            ready_ = VerifyReadLock (mLIFE_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1
     
        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjLife (ObjLifeID)

        !Arguments-------------------------------------------------------------
        integer                                    :: ObjLifeID

        !Local-----------------------------------------------------------------

        Me => FirstObjLife
        do while (associated (Me))
            if (Me%InstanceID == ObjLifeID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                          &
            stop 'ModuleLife - LocateObjLife - ERROR #1'

    end subroutine LocateObjLife

    !--------------------------------------------------------------------------

end module ModuleLife

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------



!____________________________________________________________________________________________________________
!______________________________________________________________________________________________SPARE PARTS___




!_____________________________________________________________________________### kata.Error.2.0 ###
!If (Me%ExternalVar%Mass (Me%PropIndex%DOC, Index) .lt. 0.0) stop 'DOC boom!  :ModifyLife_pre'


!_____________________________________________________________________________### TEMP TEST: kata.Error.1.0 ###
!If (Me%ExternalVar%Mass (Me%PropIndex%Nitrate, Index) .lt. 0.0) write (*,*) 'Nitrate boom!   :ModifyLife_pre'
!If (Me%ExternalVar%Mass (Me%PropIndex%Ammonia, Index) .lt. 0.0) write (*,*) 'Ammonia boom!   :ModifyLife_pre'
!If (Me%ExternalVar%Mass (Me%PropIndex%DOC, Index) .lt. 0.0) write (*,*) 'DOC boom!   :ModifyLife_pre'
!If (Me%ExternalVar%Mass (Me%PropIndex%DON, Index) .lt. 0.0) write (*,*) 'DON boom!   :ModifyLife_pre'
!If (Me%ExternalVar%Mass (Me%PropIndex%DOP, Index) .lt. 0.0) write (*,*) 'DOP boom!   :ModifyLife_pre'
!If (Me%ExternalVar%Mass (Me%PropIndex%POC, Index) .lt. 0.0) write (*,*) 'POC boom!   :ModifyLife_pre'
!If (Me%ExternalVar%Mass (Me%PropIndex%PON, Index) .lt. 0.0) write (*,*) 'PON boom!   :ModifyLife_pre'
!If (Me%ExternalVar%Mass (Me%PropIndex%POP, Index) .lt. 0.0) write (*,*) 'POP boom!   :ModifyLife_pre'
!If (Me%ExternalVar%Mass (Me%PropIndex%DOC_SL, Index) .lt. 0.0) write (*,*) 'DOC_SL boom!   :ModifyLife_pre'
!If (Me%ExternalVar%Mass (Me%PropIndex%DON_SL, Index) .lt. 0.0) write (*,*) 'DON_SL boom!   :ModifyLife_pre'
!If (Me%ExternalVar%Mass (Me%PropIndex%DOP_SL, Index) .lt. 0.0) write (*,*) 'DOP_SL boom!   :ModifyLife_pre'
!_______________________________________________________________________________________________________________
