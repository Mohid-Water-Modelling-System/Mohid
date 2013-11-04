!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : WaterQuality
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Zero-dimensional model for primary production, nitrogen and phosphorus cycle
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

Module ModuleWaterQuality
    
    use ModuleGlobalData
    use ModuleLUD
    use ModuleEnterData
    use ModuleFunctions, only: OxygenSaturation, PhytoLightLimitationFactor

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartWaterQuality
    private ::      AllocateInstance
    private ::      WQReadData
    private ::          WaterQualityOptions
    private ::          WQPropertyIndexNumber
    private ::          WQReadCalcOptions
    private ::          WQOptionsConsistencyVerif
    private ::          WQConfiguration
    private ::          WQReadFileConstants
    private ::              WQReadSilicaFileConstants   !aqui
    private ::              WQReadDiatomsFileConstants  !aqui    
    private ::      AllocateVariables
    private ::          Add_PropRateFlux
    private ::          Add_EquaRateFlux

    public  ::          Construct_WQRateFlux
   
     !Selector
    public  :: GetDTWQM
    public  :: GetWQOptions
    public  :: GetWaterQualitySize   
    public  :: GetWQPropIndex
    public  :: GetWQPropRateFlux   
    public  :: UnGetWQPropRateFlux              
                  
    public  :: GetNCRatio            
                  
    !Modifier
    public  :: WaterQuality            
    private ::      StartWaterQualityIteration
    private ::      WQCoeficientsCalculation
    private ::          WQOxygen
    private ::              WQOxygenSaturation
    private ::              WQOxygenCalculation
    private ::          WQBOD
    private ::          WQLarvae  
    private ::          WQAge     
    private ::          WQZooplankton
    private ::          WQPhytoplankton
    private ::          WQBacteria
    private ::          WQCiliate
    private ::          WQDiatoms
    private ::          WQSilica
    private ::          WQBiogenicSilica
    private ::          WQDissolvedSilica
    private ::          WQNitrogen
    private ::              WQAmmonia
    private ::              WQNitrite
    private ::              WQNitrate
    private ::              WQOrganicNitrogen
    private ::                  WQParticulateOrganicNitrogen
    private ::                  WQDONRefractory
    private ::                  WQDONNonRefractory
    private ::          WQPhosphorus
    private ::              WQOrganicPhosphorus
    private ::              WQInorganicPhosphorus
    private ::          WQPOMpools
    private ::      WQSystemResolution
   
    private ::      WQRatesCalculation

    !Destructor
    public  ::  KillWaterQuality
    private ::      DeallocateInstance
    

    !Management
    private ::      Ready
    private ::          LocateObjWaterQuality

    !Parameter-----------------------------------------------------------------

    !WQConfigurations definition
    integer, parameter                              :: ZooPhy                  = 1
    integer, parameter                              :: ZooPhyDia               = 2
    integer, parameter                              :: ZooPhyDiaCil            = 3
    integer, parameter                              :: ZooPhyDiaCilBac         = 4
    integer, parameter                              :: ZooDia                  = 5
    integer, parameter                              :: ZooDiaCilBac            = 6
    integer, parameter                              :: ZooPhyCil               = 7
    integer, parameter                              :: ZooPhyCilBac            = 8


    !Types---------------------------------------------------------------------

    type     T_PropIndex
        integer :: Zoo                              = null_int
        integer :: Larvae                           = null_int  
        integer :: Age                              = null_int  
        integer :: Phyto                            = null_int
        integer :: Bacteria                         = null_int
        integer :: Ciliate                          = null_int
        integer :: Diatoms                          = null_int !aqui
        integer :: BiogenicSilica                   = null_int !aqui
        integer :: DissolvedSilica                  = null_int !aqui
        integer :: Ammonia                          = null_int
        integer :: Nitrate                          = null_int
        integer :: Nitrite                          = null_int
        integer :: DissOrganicNitrogenRefractory    = null_int
        integer :: DONNonRefractory                 = null_int
        integer :: PartOrganicNitrogen              = null_int
        integer :: PONitrogen1                      = null_int
        integer :: PONitrogen2                      = null_int
        integer :: PONitrogen3                      = null_int
        integer :: PONitrogen4                      = null_int
        integer :: PONitrogen5                      = null_int
        integer :: POPhosphorus1                    = null_int
        integer :: POPhosphorus2                    = null_int
        integer :: POPhosphorus3                    = null_int
        integer :: POPhosphorus4                    = null_int
        integer :: POPhosphorus5                    = null_int
        integer :: PartOrganicNitrogenRefractory    = null_int
        integer :: Oxygen                           = null_int
        integer :: BOD                              = null_int
        integer :: DissOrganicPhosphorusRefractory  = null_int
        integer :: DOPNonRefractory                 = null_int
        integer :: PartOrganicPhosphorus            = null_int
        integer :: InorganicPhosphorus              = null_int
        
    end type T_PropIndex
   
    type     T_PropRateFlux
        integer                              :: ID    = null_int !initialization: jauch - or should be set to 0 (zero)?
        real, pointer, dimension(:)          :: Field => null()
        type(T_PropRateFlux), pointer        :: Next  => null()
        type(T_PropRateFlux), pointer        :: Prev  => null()
    end type T_PropRateFlux

    type     T_EquaRateFlux
        integer                              :: ID        = null_int !initialization: jauch - or should be set to 0 (zero)?
        real                                 :: scalar    = null_real
        logical                              :: TimeSerie = .false. !initialization: jauch
        type(T_EquaRateFlux), pointer        :: next              => null()
        type(T_EquaRateFlux), pointer        :: prev              => null()
        type(T_PropRateFlux), pointer        :: FirstPropRateFlux => null()
        type(T_PropRateFlux), pointer        :: LastPropRateFlux  => null()

    end type T_EquaRateFlux

    type     T_ExtraRate
        integer                              :: ID    = null_int !initialization: jauch - or should be set to 0 (zero)?
        real, pointer, dimension(:)          :: Field => null()
    end type T_ExtraRate
  

    type     T_PropCalc
        logical :: Zoo        = OFF
        logical :: Larvae     = OFF  
        logical :: Age        = OFF  
        logical :: Phyto      = OFF
        logical :: Bacteria   = OFF
        logical :: Ciliate    = OFF
        logical :: Diatoms    = OFF !aqui
        logical :: Silica     = OFF !aqui
        logical :: Nitrogen   = OFF
        logical :: Phosphorus = OFF
        logical :: Oxygen     = OFF
        logical :: Salinity   = OFF  
        logical :: BOD        = OFF
        logical :: Pompools   = OFF
    end type T_PropCalc


    type       T_CalcMethod
        logical :: ExplicitMethod = OFF
        logical :: ImplicitMethod = OFF
        logical :: SemiImpMethod  = OFF             !ppina
    end type T_CalcMethod

    type       T_External
        real, pointer, dimension(:  ) :: Salinity           => null()
        real, pointer, dimension(:  ) :: Temperature        => null()
        real, pointer, dimension(:  ) :: ShortWaveRadiation => null()
        real, pointer, dimension(:  ) :: LightExtCoefField  => null()
        real, pointer, dimension(:  ) :: Thickness          => null()
        real, pointer, dimension(:,:) :: Mass               => null()
        real, pointer, dimension(:  ) :: FishFood           => null()
    end type T_External

    type      T_SilicaCycle       !aqui
        real :: KSiBiogenicDissRate                         = null_real
        real :: BiogenicDissTCoef                           = null_real
    end type T_SilicaCycle        !aqui
    
    type      T_Diatoms   !aqui 
        real :: DiaGrowMaxRate                              = null_real
        real :: DiaGrossGrowRate                            = null_real
        real :: DiaEndogRepConst                            = null_real
        real :: DiaPhotorespFactor                          = null_real
        real :: DiaExcretionConstant                        = null_real
        real :: DiaMortMaxRate                              = null_real 
        real :: DiaMortSatConst                             = null_real 
        real :: DiaE                                        = null_real  
        real :: DiaNSatConst                                = null_real  
        real :: DiaPSatConst                                = null_real  
        real :: DiaSiSatConst                               = null_real  
        real :: DiaPhotoinhibition                          = null_real
        real :: DiaTOptMin                                  = null_real
        real :: DiaTOptMax                                  = null_real
        real :: DiaTMin                                     = null_real
        real :: DiaTMax                                     = null_real
        real :: DiaK1                                       = null_real
        real :: DiaK2                                       = null_real
        real :: DiaK3                                       = null_real
        real :: DiaK4                                       = null_real
        real :: DiaAlfaNC                                   = null_real
        real :: DiaAlfaPC                                   = null_real
        real :: DiaAlfaSiC                                  = null_real
        real :: DiaSolublInorgExcreFraction                 = null_real
        real :: DiaExcreDissOrgFraction                     = null_real
        real :: GrazDiaMin                                  = null_real
        real :: DiaRatioIngestionZoo                        = null_real
        real :: DiaZooAssimilationRate                      = null_real
        real :: ZooEfficiencyCaptureDiatoms                 = null_real
        real :: DiaLightLimitationFactor                    = null_real
        real :: DiaNutrientsLimitationFactor                = null_real
        real :: DiaNLimitationFactor                        = null_real
        real :: DiaPLimitationFactor                        = null_real
        real :: DiaSiLimitationFactor                       = null_real        
        real :: DiaTempLimitationFactor                     = null_real
        
        type(T_ExtraRate   ), pointer       :: DiaGrossProduction => null()
        type(T_ExtraRate   ), pointer       :: DiaTempLimitation  => null()
        type(T_ExtraRate   ), pointer       :: DiaNutLimitation   => null()
        type(T_ExtraRate   ), pointer       :: DiaNLimitation     => null()
        type(T_ExtraRate   ), pointer       :: DiaSiLimitation    => null()
        type(T_ExtraRate   ), pointer       :: DiaPLimitation     => null()
        type(T_ExtraRate   ), pointer       :: DiaLightLimitation => null()

    end type T_Diatoms   !aqui 

    type      T_WaterQuality      
        private
        integer                             :: InstanceID
        type(T_Size1D      )                :: Prop          
        type(T_PropIndex   )                :: PropIndex
        type(T_PropCalc    )                :: PropCalc
        type(T_CalcMethod  )                :: CalcMethod

        type(T_External    )                :: ExternalVar

        type(T_EquaRateFlux), pointer       :: FirstEquaRateFlux => null()
        type(T_EquaRateFlux), pointer       :: LastEquaRateFlux  => null()

        type(T_ExtraRate   ), pointer       :: GrossProduction   => null()
        type(T_ExtraRate   ), pointer       :: TempLimitation    => null()
        type(T_ExtraRate   ), pointer       :: NutLimitation     => null()
        type(T_ExtraRate   ), pointer       :: NLimitation       => null()
        type(T_ExtraRate   ), pointer       :: PLimitation       => null()
        type(T_ExtraRate   ), pointer       :: LightLimitation   => null()

        !aqui
        type(T_SilicaCycle )                :: SilicaCycle

        !aqui
        type(T_Diatoms     )                :: Diatoms  
        
        real, pointer, dimension(:)         :: NewProd         => null()
        real, pointer, dimension(:)         :: RegeneratedProd => null()
        real, pointer, dimension(:)         :: NitrateUptake   => null()
        real, pointer, dimension(:)         :: AmmoniaUptake   => null()

        real :: DTDay                                       = null_real
        real :: DTSecond                                    = null_real

        real :: NSatConst                                   = null_real
        real :: PhytoNutRegenerationSatConst                = null_real

        real :: AlfaPhytoNC                                 = null_real
        real :: AlfaZooNC                                   = null_real
        real :: AlfaCilNC                                   = null_real
        real :: AlfaBacteriaNC                              = null_real
        real :: OMAlfaNC                                    = null_real
        real :: OMAlfaPC                                    = null_real
        real :: BactAlfaOC                                  = null_real
        real :: PlanktonOxygenCarbonRatio                   = null_real
                        
        real :: KRefrAmmoniaMinRate                         = null_real     !KRefrAmmoniaMinRate
        real :: KNonRefrAmmoniaMinRate                      = null_real     !KNonRefrAmmoniaMinRate
        real :: NonRefrAmmoniaMinRate 
        real :: KDenitrificationRate                        = null_real
        real :: KNitrificationRateK1                        = null_real     !FIRST NITRIFICATION STEP (NH3-»NO2)
        real :: KNitrificationRateK2                        = null_real     !SECOND NITRIFICATION STEP (NO2-»NO3)
        real :: KPartDecompRate                             = null_real
        real :: PhytoAvaibleDecomp                          = null_real

        real :: TRefrAmmoniaMin                             = null_real     !TRefractoryAmmoniaMineralization
        real :: TNonRefrAmmoniaMin                          = null_real     !TNonRefractoryAmmoniaMineralization
        real :: TDenitrification                            = null_real
        real :: TNitrification                              = null_real
        real :: TPartDecomposition                          = null_real

        real :: NitrificationSatConst                       = null_real
        real :: DenitrificationSatConst                     = null_real

        real :: KDOPnrMinRate                               = null_real
        real :: TDOPnrMin                                   = null_real
        real :: KDOPrMinRate                                = null_real
        real :: TDOPrMin                                    = null_real
       
        real :: AlfaPhytoPC                                 = null_real
        real :: AlfaZooPC                                   = null_real
        real :: PSatConst                                   = null_real
        real :: AlfaCilPC                                   = null_real
        real :: AlfaBacteriaPC                              = null_real
        real :: AlfaSubstratPC                              = null_real
        
        real :: PhytoLightLimitationFactor                  = null_real
        real :: PhytoNutrientsLimitationFactor              = null_real
        real :: TPhytoLimitationFactor                      = null_real
        real :: PhytoPLimitationFactor                      = null_real
        real :: PhytoNLimitationFactor                      = null_real


        real :: GrowMaxPhytoRate                            = null_real
        real :: PhytoMortMaxRate                            = null_real
        real :: PhytoGrossGrowRate                          = null_real
        real :: E                                           = null_real
        real :: TOptPhytoMin                                = null_real
        real :: TOptPhytoMax                                = null_real
        real :: TPhytoMin                                   = null_real
        real :: TPhytoMax                                   = null_real
        real :: FK1                                         = null_real
        real :: FK2                                         = null_real
        real :: FK3                                         = null_real
        real :: FK4                                         = null_real
        real :: Photoinhibition                             = null_real
        real :: FMortSatConst                               = null_real
        real :: PhotorespFactor                             = null_real
        real :: PhytoExcretionConstant                      = null_real
        real :: PhytoEndogRepConst                          = null_real     !PhytoEndogenousRepirationConstant
        real :: PhotosynthesisOxygenCarbonRatio             = null_real 
        real :: ONMineralizationRatio                       = null_real     !Ratio between organic nitrogen and 
        real :: MinOxygen                                   = null_real     !minimum oxygen possible to allow aerobic 
                                                                            !reactions
        real :: POPDecompRate                               = null_real
        real :: TPOPDecompRate                              = null_real
        real :: GrowMaxZooRate                              = null_real
        real :: RatioOxygenCarbonZooRespiration             = null_real
        real :: RatioOxygenCarbonCilRespiration             = null_real
        real :: TOptZooMin                                  = null_real
        real :: TOptZooMax                                  = null_real
        real :: TZooMin                                     = null_real
        real :: TZooMax                                     = null_real
        real :: ZK1                                         = null_real
        real :: ZK2                                         = null_real
        real :: ZK3                                         = null_real
        real :: ZK4                                         = null_real
        real :: GrazPhytoMin                                = null_real
        real :: GrazCiliateMin                              = null_real
        real :: GrazPreyMin                                 = null_real
        real :: IvlevGrazConst                              = null_real
        real :: ZPredMortalityRate                          = null_real
        real :: TZooLimitationFactor                        = null_real

        
        real :: ZooExcretionFactor                          = null_real
        real :: ZooExcretionConst                           = null_real 
        real :: ZooIngestionConst                           = null_real
        real :: ZooEfficiencyCapturePhyto                   = null_real 
        real :: ZooEfficiencyCaptureCiliate                 = null_real 
        real :: ZooIngestionMax                             = null_real 
        real :: ZooAssimilationPhytoRate                    = null_real 
        real :: ZooAssimilationCiliateRate                  = null_real 
        real :: PhytoRatioIngestionZoo                      = null_real
        real :: CiliatesRatioIngestionZoo                   = null_real
        real :: BactRatioIngestionCiliates                  = null_real
        real :: PhytoRatioIngestionCiliates                 = null_real  
  
    
        real :: BacteriaNonGrazingMortalityRate             = null_real
        real :: BacteriaExcretionRate                       = null_real
        real :: NitrogenSaturationConstBacteria             = null_real
        real :: BacteriaMaxUptake                           = null_real
        real :: BacteriaMinSubstrate                        = null_real                    
        real :: TOptbacteriaMin                             = null_real
        real :: TOptBacteriaMax                             = null_real
        real :: TBacteriaMin                                = null_real
        real :: TBacteriaMax                                = null_real
        real :: BK1                                         = null_real
        real :: BK2                                         = null_real
        real :: BK3                                         = null_real
        real :: BK4                                         = null_real

        real :: ZooNaturalMortalityRate                     = null_real
        real :: ZooMortalityCoef                            = null_real
        real :: ZooMinMortalityRate                         = null_real
        real :: ZooMaxMortalityRate                         = null_real
        real :: PhytoExcreDissOrgFraction                   = null_real
        real :: ZooExcreDissOrgFraction                     = null_real
        real :: PhytoSolublInorgExcreFraction               = null_real
        real :: ZooSolublInorgExcreFraction                 = null_real


        real :: CiliateNaturalMortalityRate                 = null_real
        real :: CiliateMortalityCoef                        = null_real
        real :: CiliateMinMortalityRate                     = null_real
        real :: CiliateMaxMortalityRate                     = null_real
        real :: CiliateExcretionConst                       = null_real
        real :: CiliateExcretionFactor                      = null_real
        real :: CiliateIngestionConst                       = null_real 
        real :: CiliateEfficiencyCaptureBact                = null_real 
        real :: CiliateEfficiencyCapturePhyto               = null_real 
        real :: CiliateIngestionMax                         = null_real 
        real :: CiliateAssimilationBacteriaRate             = null_real 
        real :: CiliateAssimilationPhytoRate                = null_real 
        real :: CiliateGrazBactMin                          = null_real
        real :: CiliateGrazPhytoMin                         = null_real
        real :: CiliateGrazPreyMin                          = null_real

        real :: CiliateReferenceRespirationRate             = null_real

        real :: BODOxidationCoefficient                     = null_real
        real :: BODOxidationReferenceRate                   = null_real
        real :: BODOxygenSSatConstant                       = null_real
        real :: NConsOxyNitRatio                            = null_real     !NitrateConsumptionOxygenNitrateRatio
        real :: NitrificationK1_ON_ratio                    = null_real     ! NH4 -> NO2 O:N Ratio   
        real :: NitrificationK2_ON_ratio                    = null_real     ! NO2 -> NO3 O:N Ratio         
        real :: PConsOxyPhosphorusRatio                     = null_real     !PhosphateConsumptionOxygenPhosphorusRatio (Rosa)
        real :: PhytoTotalRespirationLossesRate             = null_real
        real :: NitrificationRateK1                         = null_real
        real :: NitrificationRateK2                         = null_real
        real :: BODOxidationRate                            = null_real


        real :: OxyCarbonRatio                              = null_real       !O:C Ratio in Co2

        real :: DenitrificationRate                         = null_real
        real :: PhytoNonGrazingMortalityRate                = null_real


        real :: ZooReferenceRespirationRate                 = null_real
        real :: PhytoExcretionsRate                         = null_real

        integer :: NPhases                                  = null_int
        real :: Awg                                         = null_real
        real :: Bwg                                         = null_real
        real :: Awz                                         = null_real
        real :: Bwz                                         = null_real
        real :: Atg                                         = null_real
        real :: Btg                                         = null_real
        real :: Atz                                         = null_real
        real :: Btz                                         = null_real
        real :: Ldensity                                    = null_real
        real :: Lshape                                      = null_real
        real :: Init_age                                    = null_real
        real :: Inter_age                                   = null_real
        real :: Final_age                                   = null_real
        real :: Init_length                                 = null_real
        real :: Inter_length                                = null_real
        real :: Final_length                                = null_real
        real :: FishFood_ref                                = null_real
        real :: Temperature_ref                             = null_real
        real :: Afg                                         = null_real
        
        real :: POMIngestVmax                               = null_real
        real :: PONIngestKs                                 = null_real
        real :: POPIngestKs                                 = null_real
        real :: PON_CNratio                                 = null_real
        real :: PON_CPratio                                 = null_real
        


        !aqui_10
        integer                                     :: WQConfiguration =0

        double precision, pointer, dimension(:,:)   :: Matrix  => null()
        real,             pointer, dimension(:  )   :: IndTerm => null()
        real, pointer, dimension(:  )               :: NewMass => null()  !Used with Explicit method

        !Instance of Module_EnterData
        integer                                     :: ObjEnterData = 0

        !Instance of ModuleLUD
        integer                                     :: ObjLUD = 0

        type(T_WaterQuality ), pointer              :: Next => null()

    end type T_WaterQuality

    !Global Module Variables
    type (T_WaterQuality), pointer                  :: FirstObjWaterQuality => null()
    type (T_WaterQuality), pointer                  :: Me                   => null()

    !--------------------------------------------------------------------------
    
    contains



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine StartWaterQuality(WaterQualityID, FileName, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: WaterQualityID
        character(LEN = *)              :: FileName    
        integer, optional, intent(OUT)  :: STAT     

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL
        integer                         :: ready_         

        !Local-----------------------------------------------------------------
        integer                         :: STAT_
 
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mWaterQuality_)) then
            nullify (FirstObjWaterQuality)
            call RegisterModule (mWaterQuality_) 
        endif
        
        call Ready(WaterQualityID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then
            
            call AllocateInstance           

            call Nullify_all_Sub_Type_Pointers

            call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartWaterQuality - ModuleWaterQuality - ERR01'
              

            call WQReadData
            call AllocateVariables

cd2:        if (.NOT. Me%CalcMethod%ExplicitMethod) then
                call StartLUD(Me%ObjLUD,                                                            &
                              Me%Prop%ILB,                                                          &
                              Me%Prop%IUB,                                                          &
                              Me%Prop%ILB,                                                          &
                              Me%Prop%IUB,                                                          &
                              STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'StartWaterQuality - ModuleWaterQuality - ERR02'
            end if cd2

                                                                                                        
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'StartWaterQuality - ModuleWaterQuality - ERR03'

            !Returns ID
            WaterQualityID              = Me%InstanceID

            STAT_ = SUCCESS_
        else 
            
            stop 'ModuleWaterQuality - StartWaterQuality - ERR99' 

        end if cd0


        if (present(STAT))                                                                          &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartWaterQuality

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Local-----------------------------------------------------------------
        type (T_WaterQuality), pointer           :: NewObjWaterQuality
        type (T_WaterQuality), pointer           :: PreviousObjWaterQuality


        !Allocates new instance
        allocate (NewObjWaterQuality)
        nullify  (NewObjWaterQuality%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjWaterQuality)) then
            FirstObjWaterQuality        => NewObjWaterQuality
            Me                          => NewObjWaterQuality
        else
            PreviousObjWaterQuality     => FirstObjWaterQuality
            Me                          => FirstObjWaterQuality%Next
            do while (associated(Me))
                PreviousObjWaterQuality => Me
                Me                      => Me%Next
            enddo
            Me                          => NewObjWaterQuality
            PreviousObjWaterQuality%Next=> NewObjWaterQuality
        endif

        Me%InstanceID = RegisterNewInstance (mWATERQUALITY_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    Subroutine Nullify_all_Sub_Type_Pointers

        nullify(Me%ExternalVar%Salinity    )
        nullify(Me%ExternalVar%Temperature )
        nullify(Me%ExternalVar%Mass        )
        nullify(Me%ExternalVar%FishFood    )
        nullify(Me%NewProd                 )
        nullify(Me%RegeneratedProd         )
        nullify(Me%NitrateUptake           )
        nullify(Me%AmmoniaUptake           )
        nullify(Me%Matrix                  )
        nullify(Me%IndTerm                 )
        nullify(Me%NewMass                 )

    end Subroutine Nullify_all_Sub_Type_Pointers

    !--------------------------------------------------------------------------

    subroutine Construct_WQRateFlux(WaterQualityID, WQArrayLB, WQArrayUB, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: WaterQualityID
        type(T_EquaRateFlux)    , pointer   :: EquaRateFluxX
        integer, optional                   :: STAT

        !External----------------------------------------------------------------
        integer                             :: WQArrayLB,WQArrayUB

        !Local-------------------------------------------------------------------
        type (T_PropRateFlux), pointer      :: NewPropRateFlux
        type (T_EquaRateFlux), pointer      :: NewEquaRateFlux
        
        !Begin-------------------------------------------------------------------

        integer                             :: STAT_CALL
        integer                             :: PropLB, PropUB
                                            
        integer                             :: NumAM
        integer                             :: NumPON
        integer                             :: NumPONr
        integer                             :: NumZoo
        integer                             :: NumPhyto
        integer                             :: NumDONnr
        integer                             :: NumNI
        integer                             :: NumNA
        integer                             :: NumDONr
        integer                             :: NumO
        integer                             :: NumBOD
        integer                             :: NumIP
        integer                             :: NumDOPr
        integer                             :: NumDOPnr
        integer                             :: NumPOP
        integer                             :: NumPON1
        integer                             :: NumPON2
        integer                             :: NumPON3
        integer                             :: NumPON4
        integer                             :: NumPON5
        integer                             :: NumPOP1
        integer                             :: NumPOP2
        integer                             :: NumPOP3
        integer                             :: NumPOP4
        integer                             :: NumPOP5
        integer                             :: NumAge
        integer                             :: numLarvae
        integer                             :: NumCiliate
        integer                             :: NumBacteria
        integer                             :: NumDiatom
        integer                             :: NumSiBio !aqui
        integer                             :: NumSiDiss !aqui

        logical, pointer, dimension (:)     :: LogicalEqua

        integer                             :: equa,countequa
        integer                             :: STAT_, ready_
        logical                             :: CheckName
        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(WaterQualityID, ready_)

cd0 :   if (ready_ .EQ. IDLE_ERR_) then

            numPhyto   = Me%PropIndex%Phyto
            numAM      = Me%PropIndex%Ammonia
            numNI      = Me%PropIndex%Nitrite
            numNA      = Me%PropIndex%Nitrate
            numPON     = Me%PropIndex%PartOrganicNitrogen
            numPONr    = Me%PropIndex%PartOrganicNitrogenRefractory
            numDONnr   = Me%PropIndex%DONNonRefractory
            numDONr    = Me%PropIndex%DissOrganicNitrogenRefractory
            numZoo     = Me%PropIndex%Zoo
            numBOD     = Me%PropIndex%BOD
            numO       = Me%PropIndex%Oxygen      
            numDOPr    = Me%PropIndex%DissOrganicPhosphorusRefractory          
            numDOPnr   = Me%PropIndex%DOPNonRefractory 
            numPOP     = Me%PropIndex%PartOrganicPhosphorus 
            numIP      = Me%PropIndex%InorganicPhosphorus 
            numPON1    = Me%PropIndex%PONitrogen1
            numPON2    = Me%PropIndex%PONitrogen2
            numPON3    = Me%PropIndex%PONitrogen3
            numPON4    = Me%PropIndex%PONitrogen4
            numPON5    = Me%PropIndex%PONitrogen5
            numPOP1    = Me%PropIndex%POPhosphorus1
            numPOP2    = Me%PropIndex%POPhosphorus2
            numPOP3    = Me%PropIndex%POPhosphorus3
            numPOP4    = Me%PropIndex%POPhosphorus4
            numPOP5    = Me%PropIndex%POPhosphorus5 
            numAge     = Me%PropIndex%Age
            numLarvae  = Me%PropIndex%Larvae
            NumCiliate = Me%PropIndex%Ciliate
            NumBacteria= Me%PropIndex%Bacteria
            NumDiatom  = Me%PropIndex%Diatoms         !aqui
            NumSiBio   = Me%PropIndex%BiogenicSilica  !aqui
            NumSiDiss  = Me%PropIndex%DissolvedSilica !aqui


            PropUB = Me%Prop%IUB
            PropLB = Me%Prop%ILB

            allocate(LogicalEqua(PropLB:PropUB))
            LogicalEqua =.false.
    
            countequa=0
       
            !oxigen is always computed
            Logicalequa(NumO) =.true.
            countequa = countequa + 1

            if (Me%PropCalc%Nitrogen) then
                Logicalequa(NumAM    )=.true.
                Logicalequa(NumNI    )=.true.
                Logicalequa(NumPON   )=.true.
                Logicalequa(NumDONnr )=.true.
                Logicalequa(NumDONr  )=.true.
                Logicalequa(NumNA    )=.true.
                countequa = countequa + 6
            endif

            if (Me%PropCalc%Age) then
                Logicalequa(NumAge   )=.true.
                countequa = countequa + 1
            endif

            if (Me%PropCalc%Bacteria) then
                Logicalequa(NumBacteria)=.true.
                countequa = countequa + 1
                if(Me%PropCalc%Nitrogen)then               
                    Logicalequa(NumPONr)=.true.
                    countequa = countequa + 1
                end if
            endif

            if (Me%PropCalc%Ciliate) then
                Logicalequa(NumCiliate   )=.true.
                countequa = countequa + 1
            endif

            if (Me%PropCalc%Larvae) then
                Logicalequa(NumLarvae   )=.true.
                countequa = countequa + 1
            endif
       
            if (Me%PropCalc%Phosphorus) then
                Logicalequa(numDOPr )=.true.
                Logicalequa(numDOPnr)=.true.
                Logicalequa(NumPOP  )=.true.
                Logicalequa(NumIP   )=.true.
                countequa = countequa + 4
            endif
            
            if(Me%Propcalc%Pompools) then
                if (Me%PropCalc%Nitrogen) then
                    Logicalequa(NumPON1  )=.true.
                    Logicalequa(NumPON2  )=.true.
                    Logicalequa(NumPON3  )=.true.
                    Logicalequa(NumPON4  )=.true.
                    Logicalequa(NumPON5  )=.true.
                    countequa = countequa + 5
                endif
            
                if (Me%PropCalc%Phosphorus) then
                    Logicalequa(NumPOP1 )=.true.
                    Logicalequa(NumPOP2 )=.true.
                    Logicalequa(NumPOP3 )=.true.
                    Logicalequa(NumPOP4 )=.true.
                    Logicalequa(NumPOP5 )=.true.
                    countequa = countequa + 5
                endif    
            endif

            if (Me%PropCalc%Silica) then
                Logicalequa(NumSiBio )=.true.
                Logicalequa(NumSiDiss)=.true.
                countequa = countequa + 2
            endif


            if (Me%PropCalc%Phyto) then
                Logicalequa(NumPhyto  )=.true.
                countequa = countequa + 1
            endif

            !aqui
            if (Me%PropCalc%Diatoms) then
                Logicalequa(NumDiatom  )=.true.
                countequa = countequa + 1
            endif


            if (Me%PropCalc%Zoo) then
                Logicalequa(NumZOO     )=.true.
                countequa = countequa + 1 
            endif

            if (Me%PropCalc%BOD) then
                Logicalequa(NumBOD     )=.true.
                countequa = countequa + 1 
            endif
       
            if (Me%PropCalc%Phyto) then 
            !if phyto then gross production and growth limitations parameter output are available  
            
                !GrossProduction
                allocate (Me%GrossProduction, STAT = STAT_CALL)            
                
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR01.'

                CheckName = CheckPropertyName('grossprod', number = Me%GrossProduction%ID)

                nullify(Me%GrossProduction%Field)
                allocate(Me%GrossProduction%Field   (WQArrayLB:WQArrayUB))
                !These "rates" need to be initialized or in no openpoints they will be filled with
                !negative values. This may cause errors if in between WaterQuality computations a closed
                !cell turns into open point and uses the negative "rate" values. David
                Me%GrossProduction%Field = 0.0
                
                !TempLimitation
                allocate (Me%TempLimitation, STAT = STAT_CALL)            
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR02.'

                CheckName = CheckPropertyName('temperaturelim', number = Me%TempLimitation%ID)

                nullify(Me%TempLimitation%Field)
                allocate(Me%TempLimitation%Field   (WQArrayLB:WQArrayUB))
                !Limiting factors initialized as zero are consistent with zero gross production and avoid
                !filling with other values (e.g. 1 * DT) that in no water points will not be transformed back to
                ![0-1] in recieving modules (water properties, drainage network). David
                Me%TempLimitation%Field = 0.0
                
                !NutLimitation
                allocate (Me%NutLimitation, STAT = STAT_CALL)            
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR03.'

                CheckName = CheckPropertyName('nutrientlim', number = Me%NutLimitation%ID)  

                nullify(Me%NutLimitation%field)
                allocate(Me%NutLimitation%field   (WQArrayLB:WQArrayUB))
                Me%NutLimitation%field = 0.0
                
                !NLimitation
                allocate (Me%NLimitation, STAT = STAT_CALL)            
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR03.'

                CheckName = CheckPropertyName('nitrogenlim', number = Me%NLimitation%ID)  

                nullify(Me%NLimitation%field)
                allocate(Me%NLimitation%field   (WQArrayLB:WQArrayUB))
                Me%NLimitation%field = 0.0

                !PLimitation
                allocate (Me%PLimitation, STAT = STAT_CALL)            
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR03.'

                CheckName = CheckPropertyName('phosphoruslim', number = Me%PLimitation%ID)  

                nullify(Me%PLimitation%field)
                allocate(Me%PLimitation%field   (WQArrayLB:WQArrayUB))
                Me%PLimitation%field = 0.0

                !LightLimitation
                allocate (Me%LightLimitation, STAT = STAT_CALL)            
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR04.'

                CheckName = CheckPropertyName('lightlim', number = Me%LightLimitation%ID)

                nullify(Me%LightLimitation%field)
                allocate(Me%LightLimitation%field   (WQArrayLB:WQArrayUB))            
                Me%LightLimitation%field = 0.0
                
            endif
       
           !aqui_1

            if (Me%PropCalc%Diatoms) then 
            !if diatoms then gross production and growth limitations parameter output are available  
            
            !DiaGrossProduction
                allocate (Me%Diatoms%DiaGrossProduction, STAT = STAT_CALL)            
                
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR11.'

                CheckName = CheckPropertyName('diagrossprod', number = Me%Diatoms%DiaGrossProduction%ID)

                nullify(Me%Diatoms%DiaGrossProduction%Field)
                allocate(Me%Diatoms%DiaGrossProduction%Field   (WQArrayLB:WQArrayUB))
                Me%Diatoms%DiaGrossProduction%Field = 0.0
                
            !TempLimitation
                allocate (Me%Diatoms%DiaTempLimitation, STAT = STAT_CALL)            
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR12.'

                CheckName = CheckPropertyName('diatemperaturelim', number = Me%Diatoms%DiaTempLimitation%ID)

                nullify(Me%Diatoms%DiaTempLimitation%Field)
                allocate(Me%Diatoms%DiaTempLimitation%Field   (WQArrayLB:WQArrayUB))
                Me%Diatoms%DiaTempLimitation%Field = 0.0
                
            !NutLimitation
                allocate (Me%Diatoms%DiaNutLimitation, STAT = STAT_CALL)            
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR13.'

                CheckName = CheckPropertyName('dianutrientlim', number = Me%Diatoms%DiaNutLimitation%ID)  

                nullify(Me%Diatoms%DiaNutLimitation%field)
                allocate(Me%Diatoms%DiaNutLimitation%field   (WQArrayLB:WQArrayUB))
                Me%Diatoms%DiaNutLimitation%field = 0.0
                
            !DiaNLimitation
                allocate (Me%Diatoms%DiaNLimitation, STAT = STAT_CALL)            
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR13.'

                CheckName = CheckPropertyName('dianitrogenlim', number = Me%Diatoms%DiaNLimitation%ID)  

                nullify(Me%Diatoms%DiaNLimitation%field)
                allocate(Me%Diatoms%DiaNLimitation%field   (WQArrayLB:WQArrayUB))
                Me%Diatoms%DiaNLimitation%field = 0.0
                
            !DiaSiLimitation
                allocate (Me%Diatoms%DiaSiLimitation, STAT = STAT_CALL)            
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR13.'

                CheckName = CheckPropertyName('diasilicalim', number = Me%Diatoms%DiaSiLimitation%ID)  

                nullify(Me%Diatoms%DiaSiLimitation%field)
                allocate(Me%Diatoms%DiaSiLimitation%field   (WQArrayLB:WQArrayUB))
                Me%Diatoms%DiaSiLimitation%field = 0.0

            !DiaPLimitation
                allocate (Me%Diatoms%DiaPLimitation, STAT = STAT_CALL)            
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR13.'

                CheckName = CheckPropertyName('diaphosphoruslim', number = Me%Diatoms%DiaPLimitation%ID)  

                nullify(Me%Diatoms%DiaPLimitation%field)
                allocate(Me%Diatoms%DiaPLimitation%field   (WQArrayLB:WQArrayUB))
                Me%Diatoms%DiaPLimitation%field = 0.0

            !LightLimitation
                allocate (Me%Diatoms%DiaLightLimitation, STAT = STAT_CALL)            
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR14.'

                CheckName = CheckPropertyName('dialightlim', number = Me%Diatoms%DiaLightLimitation%ID)

                nullify(Me%Diatoms%DiaLightLimitation%field)
                allocate(Me%Diatoms%DiaLightLimitation%field   (WQArrayLB:WQArrayUB))            
                Me%Diatoms%DiaLightLimitation%field = 0.0
                
            endif

            if (countequa.ne.PropUB)                                                                & 
                stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR20.'

            do equa = PropLB,PropUB

                if (Logicalequa(equa)) then

                    allocate (NewEquaRateFlux, STAT = STAT_CALL)            
                    if (STAT_CALL .NE. SUCCESS_)                                                    &
                        stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR30.'

                    nullify(NewEquaRateFlux%Prev,NewEquaRateFlux%Next)
                    nullify(NewEquaRateFlux%FirstPropRateFlux,NewEquaRateFlux%LastPropRateFlux)

                    ! Add new Property to the WaterProperties List 
                    call Add_EquaRateFlux(NewEquaRateFlux)         

                    NewEquaRateFlux%ID = equa

                    if (STAT_CALL .NE. SUCCESS_)                                                    &
                        stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR40.' 

                endif
            enddo


            EquaRateFluxX => Me%FirstEquaRateFlux  


do1:        do while (associated(EquaRateFluxX))

                do equa = PropLB,PropUB

                    if (LogicalEqua(equa)) then

                        allocate (NewPropRateFlux, STAT = STAT_CALL)            
                        if (STAT_CALL .NE. SUCCESS_)                                                 &
                            stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR50.'

                        nullify(NewPropRateFlux%field)
                        nullify(NewPropRateFlux%Prev,NewPropRateFlux%Next)
                        allocate(NewPropRateFlux%Field   (WQArrayLB:WQArrayUB))
                        
                        ! Add new Prop  
                        call Add_PropRateFlux(EquaRateFluxX, NewPropRateFlux)         

                        NewPropRateFlux%ID         = equa
                        NewPropRateFlux%Field      = 0.

                        if (STAT_CALL .NE. SUCCESS_)                                                &
                            stop 'Subroutine Construct_WQRateFlux - ModuleWaterQuality. ERR60.' 
                    endif
                enddo
        
                EquaRateFluxX => EquaRateFluxX%Next

            enddo do1
        
            nullify(EquaRateFluxX) 

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd0

        if (present(STAT))                                                                          &
            STAT = STAT_          
          
    end subroutine Construct_WQRateFlux

    !----------------------------------------------------------------------------
    !This subroutine adds a new property rateflux to the Rate Fluxes List
    subroutine Add_EquaRateFlux(NewEquaRateFlux)

        !Arguments-------------------------------------------------------------
        type(T_EquaRateFlux),pointer    :: NewEquaRateFlux

        !----------------------------------------------------------------------

        if (.not.associated(Me%FirstEquaRateFlux)) then

            Me%FirstEquaRateFlux        => NewEquaRateFlux
            Me%LastEquaRateFlux         => NewEquaRateFlux
        else
            NewEquaRateFlux%Prev        => Me%LastEquaRateFlux
            Me%LastEquaRateFlux%Next    => NewEquaRateFlux
            Me%LastEquaRateFlux         => NewEquaRateFlux
        
        end if 

    end subroutine Add_EquaRateFlux 

    !--------------------------------------------------------------------------
  
    subroutine Add_PropRateFlux(EquaRateFluxX, NewPropRateFlux)

        !Arguments-------------------------------------------------------------
        type(T_EquaRateFlux),pointer              :: EquaRateFluxX
        type(T_PropRateFlux),pointer              :: NewPropRateFlux

        !----------------------------------------------------------------------

        if (.not.associated(EquaRateFluxX%FirstPropRateFlux)) then
        
            EquaRateFluxX%FirstPropRateFlux        => NewPropRateFlux
            EquaRateFluxX%LastPropRateFlux         => NewPropRateFlux
        
        else
            
            NewPropRateFlux%Prev                   => EquaRateFluxX%LastPropRateFlux
            EquaRateFluxX%LastPropRateFlux%Next    => NewPropRateFlux
            EquaRateFluxX%LastPropRateFlux         => NewPropRateFlux
        
        end if 

    end subroutine Add_PropRateFlux 

    !--------------------------------------------------------------------------
    
    subroutine WaterQualityOptions  

        !External--------------------------------------------------------------
        integer                         :: flag
        integer                         :: FromFile
        integer                         :: STAT_CALL

        !----------------------------------------------------------------------

        call GetExtractType(FromFile = FromFile)

        call GetData(Me%PropCalc%Nitrogen,                                                          &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='NITROGEN',                                                            &
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WaterQualityOptions - ModuleWaterQuality. ERR01.' 


        call GetData(Me%PropCalc%Phosphorus,                                                        &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='PHOSPHOR',                                                            &
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WaterQualityOptions - ModuleWaterQuality. ERR02.' 


        call GetData(Me%PropCalc%Phyto,                                                             &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='PHYTO',                                                               &
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)                        
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WaterQualityOptions - ModuleWaterQuality. ERR03.' 

        call GetData(Me%PropCalc%Zoo,                                                               &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='ZOO',                                                                 &
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WaterQualityOptions - ModuleWaterQuality. ERR04.' 


        call GetData(Me%PropCalc%Larvae,                                                            &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='LARVAE',                                                              &
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WaterQualityOptions - ModuleWaterQuality. ERR07.' 


        call GetData(Me%PropCalc%Age,                                                               &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='AGE',                                                                 &
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WaterQualityOptions - ModuleWaterQuality. ERR08.'

        call GetData(Me%PropCalc%Oxygen,                                                            &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='OXYGEN',                                                              &
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)                        
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WaterQualityOptions - ModuleWaterQuality. ERR09.' 


        if (Me%PropCalc%Oxygen) then 
            Me%PropCalc%Salinity = .FALSE.
        else
            Me%PropCalc%Salinity = .TRUE.
        endif


        call GetData(Me%PropCalc%BOD,                                                               &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='BOD',                                                                 &
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)                        
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WaterQualityOptions - ModuleWaterQuality. ERR10.' 

        call GetData(Me%PropCalc%Bacteria,                                                          &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='BACTERIA',                                                            &
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)                        
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WaterQualityOptions - ModuleWaterQuality. ERR11.' 


        call GetData(Me%PropCalc%Ciliate,                                                           &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='CILIATE',                                                             &
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)                        
        if (STAT_CALL .NE. SUCCESS_)                                                                & 
            stop 'Subroutine WaterQualityOptions - ModuleWaterQuality. ERR12.'
            
        
        call GetData(Me%PropCalc%Diatoms,                                                           &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='DIATOMS',                                                             &
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)                        
        if (STAT_CALL .NE. SUCCESS_)                                                                & 
            stop 'Subroutine WaterQualityOptions - ModuleWaterQuality. ERR13.' 
        

        
        call GetData(Me%PropCalc%Silica,                                                            &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='SILICA',                                                              &
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)                        
        if (STAT_CALL .NE. SUCCESS_)                                                                & 
            stop 'Subroutine WaterQualityOptions - ModuleWaterQuality. ERR14.' 
            
            
        call GetData(Me%PropCalc%Pompools,                                                          &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='POMPOOLS',                                                            &
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)                        
        if (STAT_CALL .NE. SUCCESS_)                                                                & 
            stop 'Subroutine WaterQualityOptions - ModuleWaterQuality. ERR15.'   
        

    end subroutine WaterQualityOptions         
    
    !--------------------------------------------------------------------------
    
    !A subroutine WQPropertyIndexNumber serve para atribuir indices as propriedades a ser 
    !calculadas no modulo de Qualidade da Agua
    !
    !Ricardo C. Miranda, 1996
    subroutine WQPropertyIndexNumber

        !----------------------------------------------------------------------

        Me%Prop%ILB = 1
        Me%Prop%IUB = 0

        !Nitrogen index number
        if (Me%PropCalc%Nitrogen) then
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%Ammonia                          = Me%Prop%IUB

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%Nitrate                          = Me%Prop%IUB

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%Nitrite                          = Me%Prop%IUB

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%PartOrganicNitrogen              = Me%Prop%IUB

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%DissOrganicNitrogenRefractory    = Me%Prop%IUB

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%DONNonRefractory                 = Me%Prop%IUB
            
         endif   !Nitrogen  
    
        !Phosphorus index number
        if (Me%PropCalc%Phosphorus) then   

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%PartOrganicPhosphorus            = Me%Prop%IUB

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%DissOrganicPhosphorusRefractory  = Me%Prop%IUB

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%DOPNonRefractory                 = Me%Prop%IUB
            
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%InorganicPhosphorus              = Me%Prop%IUB
        endif   !Phosphorus


        if (Me%PropCalc%Pompools) then   
            
            if (Me%PropCalc%Nitrogen) then   
            
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%PONitrogen1                      = Me%Prop%IUB            
            
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%PONitrogen2                      = Me%Prop%IUB 

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%PONitrogen3                      = Me%Prop%IUB 
            
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%PONitrogen4                      = Me%Prop%IUB 

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%PONitrogen5                      = Me%Prop%IUB                                     
            
            endif
            
            if (Me%PropCalc%Phosphorus) then   

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%POPhosphorus1                    = Me%Prop%IUB            
            
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%POPhosphorus2                    = Me%Prop%IUB 

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%POPhosphorus3                    = Me%Prop%IUB 
            
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%POPhosphorus4                    = Me%Prop%IUB 

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%POPhosphorus5                    = Me%Prop%IUB
        
            endif
            
        endif   !POMpools


        !Silica index number
        if (Me%PropCalc%Silica) then   
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%BiogenicSilica                   = Me%Prop%IUB

            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%DissolvedSilica                  = Me%Prop%IUB

        endif   !Silica
        

        !Phytoplankton index number
        if (Me%PropCalc%Phyto) then
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%Phyto                            = Me%Prop%IUB        
        endif   !Phyto



        !Zooplankton index number
        if (Me%PropCalc%Zoo) then 
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%Zoo                              = Me%Prop%IUB
        endif   !Zoo

        
        !Larvae index number
        if (Me%PropCalc%Larvae) then
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%Larvae                           = Me%Prop%IUB
        endif   !Larvae 

        
        !Age index number
        if (Me%PropCalc%Age) then
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%Age                              = Me%Prop%IUB
        endif   !Age 

        !Bacteria index number
        if (Me%PropCalc%Bacteria) then
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%Bacteria                         = Me%Prop%IUB
            
            if (Me%PropCalc%Nitrogen) then
                Me%Prop%IUB                               = Me%Prop%IUB + 1
                Me%PropIndex%PartOrganicNitrogenRefractory= Me%Prop%IUB
            endif

        endif   !Bacteria


        !Ciliate index number
        if (Me%PropCalc%Ciliate) then
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%Ciliate                          = Me%Prop%IUB
        endif   !Ciliate


        !Oxygen index number -> The oxygen is always calculated.
        Me%Prop%IUB                                       = Me%Prop%IUB + 1
        Me%PropIndex%Oxygen                               = Me%Prop%IUB



        !BOD index number
        if (Me%PropCalc%BOD) then   
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%BOD                              = Me%Prop%IUB
        endif   !BOD


        !Diatoms index number
        if (Me%PropCalc%Diatoms) then   
            Me%Prop%IUB                                   = Me%Prop%IUB + 1
            Me%PropIndex%Diatoms                          = Me%Prop%IUB
        endif   !Diatoms
    

    end subroutine WQPropertyIndexNumber

    !--------------------------------------------------------------------------

    subroutine WQReadData

        !Local-----------------------------------------------------------------
        logical :: Consistent

        !----------------------------------------------------------------------
        

        call WaterQualityOptions
        call WQPropertyIndexNumber
        call WQReadCalcOptions

        Consistent = WQOptionsConsistencyVerif ()
    
        if (Consistent) then
            call WQConfiguration
            call WQReadFileConstants
        else
            write(*,*) 
            write(*,*) 'The Water Quality Options were not consistent, verify file data.'
            stop       'SUBROUTINE WQReadData - ModuleWaterQuality. ERR01'
        endif   !Consistent

        !----------------------------------------------------------------------

    end subroutine WQReadData

    !--------------------------------------------------------------------------

    subroutine WQReadCalcOptions

        !External--------------------------------------------------------------
        integer                         :: FromFile
        integer                         :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                         :: flag
        
        !----------------------------------------------------------------------
        call GetExtractType(FromFile = FromFile)

        !Verifica se se pretende calcular usando um metodo EXPLICITO
        call GetData(            Me%CalcMethod%ExplicitMethod,                                      &
                                 Me%ObjEnterData, flag,                                             &
                                 SearchType = FromFile,                                             &
                                 keyword='EXPLICIT',                                                &
                                 ClientModule = 'ModuleWaterQuality',                               &
                                 STAT       = STAT_CALL)                                            
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WQReadCalcOptions - ModuleWaterQuality. ERR01.' 


        !Verifica se se pretende calcular usando um metodo IMPLICITO
        call GetData(Me%CalcMethod%ImplicitMethod,                                                  &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='IMPLICIT',                                                            &    
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)                        
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WQReadCalcOptions - ModuleWaterQuality. ERR02.' 

            !Verifica se se pretende calcular usando um metodo IMPLICITO/EXPLICITO        
        call GetData(Me%CalcMethod%SemiimpMethod,                                                   & 
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='SEMIIMP',                                                             &    
                     ClientModule = 'ModuleWaterQuality',                                           &
                     STAT       = STAT_CALL)                         
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WQReadCalcOptions - ModuleWaterQuality. ERR03.'         

    end subroutine WQReadCalcOptions

    !--------------------------------------------------------------------------

    logical function WQOptionsConsistencyVerif ()

        !Local-----------------------------------------------------------------

        integer :: aux

        !----------------------------------------------------------------------

cd3 :   if (Me%PropCalc%Nitrogen .OR. Me%PropCalc%Phosphorus) then

            !aqui_2 
cd4 :       if (Me%PropCalc%Phyto .OR. Me%PropCalc%Diatoms ) then 
                WQOptionsConsistencyVerif = .TRUE.

            else
cd5 :           if (Me%PropCalc%Zoo) then
                    write(*,*) 
                    write(*,*) 'It is not possible to simulate the Water Quality with zooplankton '
                    write(*,*) 'and without phytoplankton or diatoms.                                        '
                    write(*,*) 'FUNCTION WQOptionsConsistencyVerif - ModuleWaterQuality. WARN01.'
                    write(*,*) 
                    WQOptionsConsistencyVerif = .FALSE.

                else
                    WQOptionsConsistencyVerif = .TRUE.
                end if cd5 
            end if cd4  
           
        else if (Me%PropCalc%Age .or. Me%PropCalc%Larvae .or. (Me%PropCalc%BOD .AND. Me%PropCalc%Oxygen)) then

            WQOptionsConsistencyVerif = .TRUE.

        else
        
        !aqui_22    
            write(*,*) 
            write(*,*) 'It is just possible to simulate the Water Quality with nutrients.'
            write(*,*) 'Or simple AGE or LARVAE'
            write(*,*) 'FUNCTION WQOptionsConsistencyVerif - ModuleWaterQuality. WARN02.'
            write(*,*) 
            WQOptionsConsistencyVerif = .FALSE.
        end if cd3       

        !! ------- M&M -------  
        !Para simular as bacterias é obrigatório simular os ciliados (e vice-versa)
        !A condição para simular estes dois grupos é a existência de Fito e Zoo na simulção
        !e de se estar a simular o ciclo do azoto

cd8 :   if (WQOptionsConsistencyVerif.AND.Me%PropCalc%Zoo.AND.Me%PropCalc%Nitrogen) then

cd9 :       if (Me%PropCalc%Bacteria) then

cd10 :          if (Me%PropCalc%Ciliate) then
                        WQOptionsConsistencyVerif = .TRUE.
                else
                    write(*,*) 
                    write(*,*) 'It is not possible to simulate Bacteria without Ciliate. '
                    write(*,*) 'FUNCTION WQOptionsConsistencyVerif - ModuleWaterQuality. WARN04.'
                    write(*,*) 
                    WQOptionsConsistencyVerif = .FALSE.
                end if cd10
            
            else 

cd11 :          if (Me%PropCalc%Ciliate) then
                    
                    if (Me%PropCalc%Phyto) then
                        WQOptionsConsistencyVerif = .TRUE.
                    else

                    write(*,*) 
                    write(*,*) 'It is not possible to simulate Ciliate without Bacteria or Phytoplankton. '
                    write(*,*) 'FUNCTION WQOptionsConsistencyVerif - ModuleWaterQuality. WARN05.'
                    write(*,*) 
                    WQOptionsConsistencyVerif = .FALSE.
                    endif
                else 
                WQOptionsConsistencyVerif = .TRUE.
                end if cd11

            end if cd9

        else 

cd12 :      if (Me%PropCalc%Bacteria.OR.Me%PropCalc%Ciliate) then
                    write(*,*) 
                    write(*,*) 'It is not possible to simulate Ciliate and Bacteria without Phytoplankton and Zooplankton. '
                    write(*,*) 'FUNCTION WQOptionsConsistencyVerif - ModuleWaterQuality. WARN05.'
                    write(*,*) 
                    WQOptionsConsistencyVerif = .FALSE.
            else
            WQOptionsConsistencyVerif = .TRUE.
            end if cd12

        end if cd8
! ------- M&M -------


cd7 :   if (WQOptionsConsistencyVerif) then
cd2 :   if (Me%PropCalc%BOD .AND. (.NOT. Me%PropCalc%Oxygen)) then
            write(*,*) 
            write(*,*) 'It is not possible to simulate the Water Quality with BOD and without Oxygen.'
            write(*,*) 'FUNCTION WQOptionsConsistencyVerif - ModuleWaterQuality. WARN03.'
            write(*,*) 
            WQOptionsConsistencyVerif = .FALSE.
        end if cd2
        end if cd7

        !aqui_222
        if (WQOptionsConsistencyVerif) then
        
            if (Me%PropCalc%Diatoms .AND. (.NOT. Me%PropCalc%Silica)) then
                write(*,*) 
                write(*,*) 'It is not possible to simulate the Water Quality with Diatoms and without Silica.'
                write(*,*) 'FUNCTION WQOptionsConsistencyVerif - ModuleWaterQuality. WARN04.'
                write(*,*) 
                WQOptionsConsistencyVerif = .FALSE.
            end if 
        
        end if 

cd6 :   if (WQOptionsConsistencyVerif) then
            aux = 0
            if (Me%CalcMethod%ExplicitMethod) aux = aux + 1
            if (Me%CalcMethod%ImplicitMethod) aux = aux + 1
            if (Me%CalcMethod%SemiImpMethod ) aux = aux + 1


cd1 :       if (aux .EQ. 1) then
                WQOptionsConsistencyVerif = .TRUE.

            else 
                WQOptionsConsistencyVerif = .FALSE.
            end if cd1
        end if cd6
           
        
cd90 :  if (WQOptionsConsistencyVerif) then
        
cd91 :       if (Me%PropCalc%Pompools) then
            
                if ((.NOT. Me%PropCalc%Nitrogen) .AND. (.NOT. Me%PropCalc%Phosphorus)) then
                    write(*,*) 
                    write(*,*) 'Impossible to simulate the Water Quality with POM pools without'
                    write(*,*) 'at least one nutrient cycle (N or P).'
                    write(*,*) 'FUNCTION WQOptionsConsistencyVerif - ModuleWaterQuality. WARN06a.'
                    write(*,*) 
                    WQOptionsConsistencyVerif = .FALSE.
                end if
                 
                if (WQOptionsConsistencyVerif .AND. (.NOT. Me%PropCalc%Oxygen)) then
                    write(*,*) 
                    write(*,*) 'Impossible to simulate the Water Quality with POM pools without oxygen.'
                    write(*,*) 'FUNCTION WQOptionsConsistencyVerif - ModuleWaterQuality. WARN06b.'
                    write(*,*) 
                    WQOptionsConsistencyVerif = .FALSE.
                end if 
                          
             end if cd91
                     
        end if cd90
                


        !----------------------------------------------------------------------

    end function WQOptionsConsistencyVerif

    !--------------------------------------------------------------------------

    subroutine WQConfiguration

        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        
        !----------------------------------------------------------------------

cd110 :  if (Me%PropCalc%Zoo) then
        
cd21 :      if (Me%PropCalc%Phyto) then

cd31 :          if (Me%PropCalc%Diatoms) then

cd41 :              if (Me%PropCalc%Ciliate) then

cd51 :                  if (Me%PropCalc%Bacteria) then
                            Me%WQConfiguration = 4                        
                        else cd51
                            Me%WQConfiguration = 3
                        endif cd51
                    
                    else cd41
                        
                        Me%WQConfiguration = 2
                    
                    endif cd41

                else cd31

cd61 :               if (Me%PropCalc%Ciliate) then

cd71 :                  if (Me%PropCalc%Bacteria) then
                            Me%WQConfiguration = 8                        
                        else cd71
                            Me%WQConfiguration = 7
                        endif cd71
                    
                    else cd61
                       
                        Me%WQConfiguration = 1
                    
                    endif cd61
                
                endif cd31

            else cd21

cd81 :          if (Me%PropCalc%Diatoms) then

cd91 :              if (Me%PropCalc%Ciliate) then

cd100 :                 if (Me%PropCalc%Bacteria) then
                            Me%WQConfiguration = 6                        
                        endif cd100
    
                    else cd91
        
                        Me%WQConfiguration = 5
                    
                    endif cd91

                endif cd81
            
            endif cd21
        
        else cd110
        
            Me%WQConfiguration = 0
        
        end if cd110 

    end subroutine WQConfiguration

    !--------------------------------------------------------------------------

    
    subroutine WQReadFileConstants

        !External--------------------------------------------------------------
        integer                         :: FromFile
        integer                         :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                         :: flag

        !----------------------------------------------------------------------

        call GetExtractType(FromFile = FromFile)

cd1 :   if (Me%DTSecond .LE. 0.0) then
            !DTSecond, time step, in seconds, between 2 WaterQuality calls 
            call GetData(           Me%DTSecond,                                    &
                                    Me%ObjEnterData, flag,                          &
                                    SearchType = FromFile, keyword='DTSECONDS',     & 
                                    default    = 60.0 * 60.0,                       & !1 hour
                                    ClientModule = 'ModuleWaterQuality',            &
                                    STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR00.' 

cd22 :      if (flag .EQ. 0) then
                write(*,*) 
                write(*,*) 'Keyword DTSECONDS not found in Water quality data file.'
                write(*,*) 'Subroutine WQReadFileConstants - ModuleWaterQuality. WRN01.'
                write(*,*) 'Assumed ', Me%DTSecond, &
                            'seconds (',  Me%DTSecond / 3600.0, 'hour).'
                write(*,*) 
            end if cd22
        end if cd1

        
        !For compatibility with the rest of the program, !DTSeconds converted to day!
        Me%DTDay = Me%DTSecond / 24.0 / 60.0 / 60.0

        !Reads non specific rates & constants--------------------------------------

        !NSatConst, nitrogen half-saturation constant, mgN/l
        call GetData(           Me%NSatConst,                                       &
                                Me%ObjEnterData, flag,                              &
                                SearchType = FromFile,  keyword ='NSATCONS',        & 
                                default    = 0.014,                                 &
                                ClientModule = 'ModuleWaterQuality',                &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR01.' 


        !Photoinhibition, phytoinhibition, W/m2                                                                                   
        call GetData(           Me%Photoinhibition,                                 &
                                Me%ObjEnterData, flag,                              &
                                SearchType   = FromFile, keyword ='PHOTOIN',        & 
                                default      = 121.0,                               &
                                ClientModule = 'ModuleWaterQuality',                &
                                STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR02.' 

        !PhytoNutRegenerationSatConst, phytoplankton nutrient regeneration half saturation rate, 
        !mgC/l
        call GetData(            Me%PhytoNutRegenerationSatConst,                   &
                                 Me%ObjEnterData, flag,                             &
                                 SearchType = FromFile, keyword='FREGSATC',         &       
                                 default    = 1.0,                                  &
                                 ClientModule = 'ModuleWaterQuality',               &
                                 STAT       = STAT_CALL)
        if (               STAT_CALL .NE. SUCCESS_)                                 &
            stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR03.' 


        !AlfaPhytoNC, phytoplankton ratio between Nitrogen and Carbon, mgN/mgC
        call GetData(            Me%AlfaPhytoNC,                                    &
                                 Me%ObjEnterData, flag,                             &
                                 SearchType = FromFile, keyword='FRATIONC',         & 
                                 default    = 0.18,                                 & !Redfield ratio
                                 ClientModule = 'ModuleWaterQuality',               &
                                 STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR04.' 

        !Reads rates & constants to the Nitrogen simulation------------------------


        !OMAlfaNC, organic matter ratio between Nitrogen and Carbon, mgN/mgC
        call GetData(            Me%OMAlfaNC,                                       &
                                 Me%ObjEnterData, flag,                             &
                                 SearchType = FromFile, keyword='OMRATIONC',        & 
                                 default    = 0.18,                                 & !Redfield ratio
                                 ClientModule = 'ModuleWaterQuality',               &
                                 STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR04a.' 

        !BactAlfaOC, organic matter ratio between Oxygen and Carbon, mgO/mgC
        call GetData(            Me%BactAlfaOC,                                     &
                                 Me%ObjEnterData, flag,                             &
                                 SearchType = FromFile, keyword='BACTRATIOOC',      & 
                                 default    = 1.4,                                  & !Redfield ratio
                                 ClientModule = 'ModuleWaterQuality',               &
                                 STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR04a.' 

        
        
        !OMAlfaPC, organic matter ratio between Nitrogen and Carbon, mgN/mgC
        call GetData(            Me%OMAlfaPC,                                       &
                                 Me%ObjEnterData, flag,                             &
                                 SearchType = FromFile, keyword='OMRATIOPC',        & 
                                 default    = 0.024,                                 & !Redfield ratio
                                 ClientModule = 'ModuleWaterQuality',               &
                                 STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR04b.' 



            !AlfaZooNC, zooplankton ratio between Nitrogen and Carbon, mgN/mgC
            call GetData(            Me%AlfaZooNC,                                  &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword='ZRATIONC',     & 
                                     default    = 0.15,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR07.' 

            !KRefrAmmoniaMinRate, reference ammonia mineralization rate of the
            !   refractory DON, 1/T
            call GetData(            Me%KRefrAmmoniaMinRate,                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword='NMINR',       &
                                     default    = 0.01,                             & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR08.' 

            !KPartDecompRate, reference particulate organic Nitrogen decomposition rate, 1/T
            call GetData(          Me%KPartDecompRate,                              &
                                   Me%ObjEnterData, flag,                           &
                                   SearchType = FromFile,  keyword='NOPREF',        &
                                   default    = 0.1,                                & !1/day
                                   ClientModule = 'ModuleWaterQuality',             &
                                   STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR10.' 

            !PhytoAvaibleDecomp: Fraction of PON available for mineralization
            call GetData(            Me%PhytoAvaibleDecomp,                         &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword='PHDECOMP',    &
                                     default    = 0.7,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR11.' 


            !KDenitrificationRate, reference denitirfication rate, 1/T
            call GetData(            Me%KDenitrificationRate,                       &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword='DENITREF',    &
                                     default    = 0.125,                             & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR12.' 

            !KNitrificationRateK1, reference nitrification rate, 1/T
            call GetData(            Me%KNitrificationRateK1,                       &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword='NITRIREFK1',  &
                                     default    = 0.02,                             & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR13.' 
                
            !KNitrificationRateK2, reference nitrification rate, 1/T
            call GetData(            Me%KNitrificationRateK2,                       &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword='NITRIREFK2',  &
                                     default    = 0.25,                             & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR14.' 



            !TRefrAmmoniaMin, DONren mineralization temperature coefficient 
            !   of the refractory DON
            call GetData(            Me%TRefrAmmoniaMin,                            &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword= 'TMINR',       &
                                     default    = 1.02,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR15.' 


            !TPartDecomposition, particulate organic Nitrogen decomposition temperature 
            !coefficient
            call GetData(            Me%TPartDecomposition,                         &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,   keyword='NOPCOEF',    & 
                                     default    = 1.02,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR16.' 



            !TDenitrification, denitirfication temperature coefficient
            call GetData(            Me%TDenitrification,                           &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword ='TDENCOEF',   & 
                                     default    = 1.045,                            &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR17.' 

            !TNitrification, nitrification temperature coefficient
            call GetData(            Me%TNitrification,                             &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword ='TNITCOEF',   & 
                                     default    = 1.047,                            &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR18.' 

            !NitrificationSatConst, nitrification semi-saturation constant, mgO2/l
            call GetData(            Me%NitrificationSatConst,                      &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword ='NITSATCO',   & 
                                     default    = 2.0,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR19.' 

            !DenitrificationSatConst, denitrification semi-saturation constant, mgO2/l
            call GetData(            Me%DenitrificationSatConst,                    &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,   keyword ='DENSATCO',  & 
                                     default    = 0.1,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR20.' 

            !PhytoSolublInorgExcreFraction, soluble inorganic fraction of the plankton excretions
            call GetData(            Me%PhytoSolublInorgExcreFraction,              &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword='FSOLEXCR',    &  
                                     default    = 0.4,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR21.' 


            !ExcreDissOrgFraction, dissolved organic fraction of the plankton excretions
            call GetData(            Me%PhytoExcreDissOrgFraction,                  &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword='FDISSDON',    &   
                                     default    = 0.5,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR22.' 

            ! ------- M&M -------
            !   Pergunta para determinar se o programa vai ler as taxas e constantes referentes à opção
            !   de simulação com ou sem Bact. & Cil.   

                !PlanktonOxygenCarbonRatio
                call GetData(            Me%PlanktonOxygenCarbonRatio,                  &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword='PLANK_OC_RAT', & 
                                         default    = 32./12.,                            &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR23.' 

                !ZooSolublInorgExcreFraction, soluble inorganic fraction of the zooplankton excretions
                call GetData(            Me%ZooSolublInorgExcreFraction,                &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword='ZSOLEXCR',     &    
                                         default    = 0.4,                             &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR24.' 


                !ZooExcreDissOrgFraction, dissolved organic fraction of the zooplankton excretions
                call GetData(            Me%ZooExcreDissOrgFraction,                    &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword='ZDISSDON',     &
                                         default    = 0.5,                             &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR25.' 


                !KNonRefrAmmoniaMinRate, reference ammonia mineralization rate of the
                !   non refractory DON, 1/T
                call GetData(            Me%KNonRefrAmmoniaMinRate,                     &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword='NMINENR',      &
                                         default    = 0.1,                              & !1/day
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)

                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR26.' 

                !TNonRefrAmmoniaMin, Nitrogen mineralization temperature coefficient 
                !   of the non refractory DON
                call GetData(          Me%TNonRefrAmmoniaMin,                           &
                                       Me%ObjEnterData, flag,                           &
                                       SearchType = FromFile,keyword='TMINNR',          & 
                                       default    = 1.02,                               &
                                       ClientModule = 'ModuleWaterQuality',             &
                                       STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR27.' 





        !Reads rates & constants to the Phosphorus simulation------------------------


            
            !KDOPnrMinRate, reference DOPnr mineralization rate, 1/T
            call GetData(            Me%KDOPnrMinRate,                              &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,   keyword ='PMINNR',    & 
                                     default    = 0.1,                             & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR28.' 

            !TDOPnrMin, Phosphorus mineralization temperature coefficient
            call GetData(            Me%TDOPnrMin,                                  &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,   keyword ='PMINNRCOEF',&
                                     default    = 1.064,                            &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR29.' 


            
            !KDOPrMinRate, reference DOPrefractary mineralization rate, 1/T
            call GetData(            Me%KDOPrMinRate,                               &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,   keyword ='PMINR',     & 
                                     default    = 0.03,                              & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR30.' 



            !TDOPnrMin, Phosphorus mineralization temperature coefficient
            call GetData(            Me%TDOPrMin,                                   &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,   keyword ='PMINRCOEF', &
                                     default    = 1.064,                            &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR31.' 


   
            !POPDecompRate, reference POP mineralization rate, 1/T
            call GetData(            Me%POPDecompRate,                              &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,   keyword ='PPARTMIN',  & 
                                     default    = 0.2,                             & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR32.' 



            !TPOPMin, Phosphorus mineralization temperature coefficient
            call GetData(            Me%TPOPDecompRate,                             & 
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,keyword ='TPPARTMINCOEF',&
                                     default    = 1.08,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR33.' 


            !AlfaPhytoPC, phytoplankton ratio between Phosphorus and Carbon, mgP/mgC
            call GetData(            Me%AlfaPhytoPC,                                &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword ='FRATIOPC',   &
                                     default    = 0.024,                            & !Redfield ratio
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR34.' 


            !AlfaZooPC, zooplankton ratio between Phosphorus and Carbon, mgP/mgC
            call GetData(            Me%AlfaZooPC,                                  &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword ='ZRATIOPC',    &
                                     default    = 0.024,                            &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR35.' 




        !Reads the rates & constants to the Phytoplankton simulation---------------
 

            !PSatConst, Phosphorus half-saturation constant, phosphorus, M/L^3
            call GetData(            Me%PSatConst,                                  &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword ='PSATCONS',    &
                                     default    = 0.001,                            & !mgP/l
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR36.' 



            !GrowMaxPhytoRate, maximum phytoplankton growth rate, 1/T
            call GetData(            Me%GrowMaxPhytoRate,                           &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword ='GROWMAXF',   &
                                     default    = 2.,                              & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR37.' 


            !PhytoMortMaxRate, phytoplankton maximum mortality, carbon, M/(L^3.T)
            call GetData(            Me%PhytoMortMaxRate,                           &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,   keyword ='FMORTMAX',  &  
                                     default    = 0.02,                             & !mgC/l/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR38.' 




            !TOptPhytoMin, minimum temperature of the optimal interval for the phytoplankton 
            !growth, oC
            call GetData(            Me%TOptPhytoMin,                               &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword ='TOPTFMIN',    &
                                     default    = 25.0,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR39.' 


            !TOptPhytoMax, maximum temperature of the optimal interval for the phytoplankton 
            !growth, oC
            call GetData(            Me%TOptPhytoMax,                               &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword ='TOPTFMAX',    &  
                                     default    = 26.5,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR40.' 


            !TPhytoMin, minimum tolerable temperature of the  interval for the phytoplankton 
            !growth, oC
            call GetData(            Me%TPhytoMin,                                  &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword ='TFMIN',      &
                                     default    = 4.0,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR41.' 


            !TPhytoMax, maximum tolerable temperature of the  interval for the phytoplankton 
            !growth, oC
            call GetData(            Me%TPhytoMax,                                  &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword ='TFMAX',       &
                                     default    = 37.0,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR42.' 



            !FK1, constant to control temperature response curve shape
            call GetData(            Me%FK1,                                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword ='TFCONST1',   &
                                     default    = 0.05,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR43.' 



            !FK2, constant to control temperature response curve shape
            call GetData(            Me%FK2,                                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword ='TFCONST2',    &
                                     default    = 0.98,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR44.' 


            !FK3, constant to control temperature response curve shape
            call GetData(            Me%FK3,                                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword ='TFCONST3',    &
                                     default    = 0.98,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR45.' 


            !FK4, constant to control temperature response curve shape
            call GetData(            Me%FK4,                                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword ='TFCONST4',    & 
                                     default    = 0.02,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR46.' 



            !FMortSatConst, mortality half saturation rate, M/(L^3.T)
            call GetData(            Me%FMortSatConst,                              &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword ='FMORTCON',    & 
                                     default    = 0.3,                              & !mgC/l/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR47.' 



            !PhotorespFactor, fraction of actual photosynthesis which is oxidised by
            !photorespiration
            call GetData(            Me%PhotorespFactor,                            &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword ='PHOTORES',    & 
                                     default    = 0.018,                            &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR48.' 




            !PhytoEndogRepConst, Phytoplankton endogenous respiration constant, 
            !1 / T 
            call GetData(            Me%PhytoEndogRepConst,                         &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword ='FENDREPC',    &
                                     default    = 0.0175,                           & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR49.' 

        
        
            !PhytoExcretionConstant, excretion constant
            call GetData(            Me%PhytoExcretionConstant,                     &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'EXCRCONS',   &
                                     default    = 0.08,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR50.' 


            ! ------- M&M -------
            !   Pergunta para determinar se o programa vai ler as taxas e constantes referentes à opção
            !   de simulação com ou sem Bact. & Cil.   

            !E, assimilation efficiency of the phytoplankton by the zooplankton
                call GetData(            Me%E,                                          &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,   keyword ='ASS_EFIC',  &
                                         default    = 0.8,                              & !60%
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR51.' 

        !Reads the rates & constants to the Zooplankton simulation-----------------

            !TOptZooMin, minimum temperature of the optimal interval for the zooplankton growth, 
            !oC
            call GetData(            Me%TOptZooMin,                                 &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TOPTZMIN',   &
                                     default    = 24.8,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR52.' 


            !TOptZooMax, maximum temperature of the optimal interval for the zooplankton growth, 
            !oC
            call GetData(            Me%TOptZooMax,                                 &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TOPTZMAX',   &
                                     default    = 25.1,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR53.' 

            !TZooMin, minimum tolerable temperature of the  interval for the zooplankton growth, 
            !oC
            call GetData(            Me%TZooMin,                                    &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TZMIN',      &
                                     default    = 5.0,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR54.' 


            !TZooMax, maximum tolerable temperature of the  interval for the zooplankton growth, 
            !oC
            call GetData(            Me%TZooMax,                                    &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TZMAX',      &
                                     default    = 35.0,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR55.' 


            !ZK1, constant to control temperature response curve shape
            call GetData(            Me%ZK1,                                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TZCONST1',   &
                                     default    = 0.05,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR56.' 

            !ZK2, constant to control temperature response curve shape
            call GetData(            Me%ZK2,                                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TZCONST2',   &
                                     default    = 0.98,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                             &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR57.' 




            !ZK3, constant to control temperature response curve shape
            call GetData(            Me%ZK3,                                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TZCONST3',   &
                                     default    = 0.98,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR58.' 


            !ZK4, constant to control temperature response curve shape
            call GetData(            Me%ZK4,                                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TZCONST4',   &
                                     default    = 0.02,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR59.' 

            !ZooReferenceRespirationRate, rate of consumption of Carbon by respiration and 
            !non-predatory mortality at the reference temperature, 1/T
            call GetData(            Me%ZooReferenceRespirationRate,                &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'ZREFRESP',   &
                                     default    = 0.02,                            & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR60.' 

            !GrazPhytoMin, minimum phytoplankton concentration for the existence of grazing, mgC/l
            call GetData(            Me%GrazPhytoMin,                               &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'GRAZFITOMIN',&
                                     default    = 0.0045,                            &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR50.'
           
            !GrazPreyMin, minimum phytoplankton concentration for the existence of grazing, mgC/l
            call GetData(            Me%GrazPreyMin,                               &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'ZOOPREYMIN',&
                                     default    = 0.045,                            &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR50.' 
 
        !Reads the rates & constants to the Ciliates and Bacteria simulation-----------------

                !GrazCiliateMin, minimum phytoplankton concentration for the existence of grazing, mgC/l
                call GetData(      Me%GrazCiliateMin,                                   &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword = 'GRAZCILMIN', &
                                         default    = 0.0045,                            &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR61.' 

        

                !ZooExcretionFactor, Zooplankton Excretion Rate
                call GetData(           Me%ZooExcretionFactor,                          &
                                        Me%ObjEnterData, flag,                          &
                                        SearchType = FromFile, keyword =  'ZEXCFAC',    &
                                        default    = 0.02,                            &
                                        ClientModule = 'ModuleWaterQuality',            &
                                        STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR62.' 

                !ZooExcretionConstant, excretion constant
                call GetData(            Me%ZooExcretionConst,                          &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword = 'ZEXCCONS',   &
                                         default    = 1.0305,                           &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR63.' 


                !ZooMortalityCoef  
                call GetData(           Me%ZooMortalityCoef,                            &
                                        Me%ObjEnterData, flag,                          &
                                        SearchType = FromFile, keyword = 'MORTZCOEF',   &
                                        default    = 0.0,                               &
                                        ClientModule = 'ModuleWaterQuality',            &
                                        STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR64.' 



                !ZooMinMortalityRate  
                call GetData(            Me%ZooMinMortalityRate,                        &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword = 'MINMORTZ',   &
                                         default    = 0.0,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR65.' 


                !ZooMaxMortalityRate
                call GetData(            Me%ZooMaxMortalityRate,                        &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword = 'MAXMORTZ',   &
                                         default    = 0.040,                            & !1/day
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR66.' 
    
            
                !ZooIngestionConst, Half-Saturation Constant for Grazing
                call GetData(Me%ZooIngestionConst,                                      &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,keyword = 'INGCONSZ',    &
                                         default    = 0.85,                             &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR67.' 

                                            
                !ZooEfficiencyCapturePhyto  
                call GetData(            Me%ZooEfficiencyCapturePhyto,                  &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword = 'ZOOEFFCAPHY',   &
                                         default    = 0.8,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)

                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR68.' 

                        
                !ZooEfficiencyCaptureCiliate
                 call GetData(           Me%ZooEfficiencyCaptureCiliate,                &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword = 'ZOOEFFCAPCIL' , &
                                         default    = 0.2,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR69.' 

                                
                !ZooIngestionMax    
                call GetData(            Me%ZooIngestionMax,                            &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword =  'ZINGMAX',   &
                                         default    = 2.,                             &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR70.' 
                                    
                !ZooAssimilationPhytoRate       
                call GetData(            Me%ZooAssimilationPhytoRate,                   &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,keyword = 'ZOPHYASS',    &
                                         default    = 0.8,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)

                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR71.' 



                !ZooAssimilationCiliateRate 
                call GetData(            Me%ZooAssimilationCiliateRate,                 &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,keyword = 'ZOCILASS',    &
                                         default    = 0.8,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR72.' 

                !PhytoRatioIngestionZoo 
                call GetData(            Me%PhytoRatioIngestionZoo,                     &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword = 'PHYRATING',  &
                                         default    = 0.3,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR673.' 
        
                !PhytoRatioIngestionZoo 
                call GetData(            Me%CiliatesRatioIngestionZoo,                  &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword = 'CILRATINGZOO',&
                                         default    = 0.3,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR673.' 



    
                !GrowMaxZooRate, maximum zooplankton growth rate, 1/T
                call GetData(            Me%GrowMaxZooRate,                             &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword = 'GROWMAXZ',   &
                                         default    = 0.3,                             & !1/day
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR74.' 

            
                !IvlevGrazConst, Ivlev grazing constant
                call GetData(            Me%IvlevGrazConst,                             &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword = 'IVLEVCON',   &
                                         default    = 1.6,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR75.' 


            
            !ZPredMortalityRate, predatory mortality rate, 1/T
            call GetData(            Me%ZPredMortalityRate,                         &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'ZPREDMOR',   &
                                     default    = 0.0077,                           & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR76.' 




        !Reads the rates & constants to the Oxygen simulation----------------------


            !PhotosynthesisOxygenCarbonRatio, Photosynthesis Oxygen:Carbon ratio, (M/L^3)/(M/L^3)
            call GetData( Me%PhotosynthesisOxygenCarbonRatio,                       &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,keyword = 'PHOTOSOC',    &
                                     default    = 32.0 / 12.0,                      & !mgO2 / mgC
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR77.' 


            !RatioOxygenCarbonZooRespiration, Zooplankton respiration Oxygen:Carbon ratio, 
            !mgO2 / mgC
            call GetData(            Me%RatioOxygenCarbonZooRespiration,            &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,keyword = 'ZOCRATIO',    &
                                     default    = 32.0 / 12.0,                      &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR78.' 


            !NConsOxyNitRatio, secondary Oxygen production due to Nitrate 
            !consumption, (M/L^3)/(M/L^3)
            call GetData(            Me%NConsOxyNitRatio,                           &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,keyword = 'NITONRAT',    &
                                     default    = 48.0 / 14.0,                      & !mgO2 / mgN
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR79a.' 
                
            
            call GetData(            Me%NitrificationK1_ON_ratio,                   &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,keyword = 'NITK1_ONRAT', &
                                     default    = (1.5*32) / 14.0,                      & !mgO2 / mgN
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR79b.'
                
            call GetData(            Me%NitrificationK2_ON_ratio,                   &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,keyword = 'NITK2_ONRAT', &
                                     default    = (0.5*32) / 14.0,                      & !mgO2 / mgN
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR79c.'    
                
                
                
                

            !     if (flag .EQ. 0)                                                   &
            !         Me%NConsOxyNitRatio = 48.0 / 14.0          !mgO2 / mgN
            !        NO3 -> N +O3
            !         Me%NConsOxyNitRatio = 32.0 / 14.0          !mgO2 / mgN

            !Rosa
            !PConsOxyPhosphorusRatio, secondary Oxygen production due to phosphate consumption
            call GetData(            Me%PConsOxyPhosphorusRatio,                    &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,keyword = 'PHOSOPRAT',   &
                                     default    = 64.0 / 31.0,                      & !mgO2 / mgN
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR80.' 
    
            call GetData(            Me%MinOxygen,                                  &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,keyword = 'MINOXYGEN',   &
                                     default    = 10E-5,                            &  !mgO2 /L
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR83.' 


            !O:C Ratio in CO2
            call GetData(            Me%OxyCarbonRatio,                             &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'OCRATIO',    &
                                     default    = 32.0 / 12.0,                      & !mgO2/mgC
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR84a.' 



        !Reads the rates & constants to the BOD simulation-------------------------




            !BODOxidationCoefficient, BOD oxidation coefficient
            call GetData(            Me%BODOxidationCoefficient,                    &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'BODCOEF',    &
                                     default    = 1.047,                            &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR85.' 


            !BODOxidationReferenceRate, reference BOD oxidation, 1/T
            call GetData(            Me%BODOxidationReferenceRate,                  &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'BODREF',     &
                                     default    = 0.18,                             & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR86.' 


            !BODOxygenSSatConstant, Oxygen limitation half-saturation constant, 1/T
            call GetData(            Me%BODOxygenSSatConstant,                      &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword =  'BODOSSAT',  &
                                     default    = 0.5,                              & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR87.' 


            !AlfaCilNC, ciliate ratio between Nitrogen and Carbon, mgN/mgC
            call GetData(            Me%AlfaCilNC,                                  &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'CRATIONC',   &
                                     default    = 0.16,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR88a.' 


            !AlfaCilPC, ciliate ratio between Phosp and Carbon, mgP/mgC
            call GetData(            Me%AlfaCilPC,                                  &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'CRATIOPC',   &
                                     default    = 0.024,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR88b.' 


!Bacteria-----------------------------------------------------------------------------------------
    
            !AlfaBacteriaNC, bacteria ratio between Nitrogen and Carbon, mgN/mgC
            call GetData(            Me%AlfaBacteriaNC,                             &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'BRATIONC',   &
                                     default    = 0.2,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR89.' 


           !Reads the rates & constants to the Bacteria simulation---------------------- 
            !BacteriaNonGrazingMortalityRate
            call GetData(            Me%BacteriaNonGrazingMortalityRate,            &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword ='NATMORB',     &
                                     default    = 0.1,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR91.' 


            !BacteriaExcretionRate   day-1
            call GetData(  Me%BacteriaExcretionRate    ,                            &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'BARESPCO',   &
                                     default    = 2.5,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR92.' 

            !BacteriaMaxUptake
            call GetData(            Me%BacteriaMaxUptake,                          &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,  keyword = 'BMAXUPTA',  &
                                     default    = 0.25,                             & ! 6.6 day-1  /24 
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR93.' 


            !BacteriaMinSubstrate
            call GetData(            Me%BacteriaMinSubstrate,                       &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'BACMINSUB',  &
                                     default    = 0.010,                            &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR95.' 


            !NitrogenSaturationConstBacteria
            call GetData(            Me%NitrogenSaturationConstBacteria,            &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'BACNCONS',   &
                                     default    = 0.0008,                           & ! mgN/l
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR96.' 


            !TOptBacteriaMin, minimum temperature of the optimal interval for the Bacteria growth, 
            !oC
            call GetData(            Me%TOptBacteriaMin,                            &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TOPTBMIN',   &
                                     default    = 24.8,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR97.' 

            !TOptBacteriaMax, maximum temperature of the optimal interval for the Bacteria growth, 
            !oC
            call GetData(            Me%TOptBacteriaMax,                            &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TOPTBMAX',   &
                                     default    = 25.1,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR98.' 


            !TBacteriaMin, minimum tolerable temperature of the  interval for the Bacteria growth, 
            !oC
            call GetData(            Me%TBacteriaMin,                               &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TBMIN',      &
                                     default    = 5.0,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR99.' 


            !TBacteriaMax, maximum tolerable temperature of the  interval for the Bacteria growth, 
            !oC
            call GetData(            Me%TBacteriaMax,                               &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TBMAX',      &
                                     default    = 35.0,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR100.' 


            !BK1, constant to control temperature response curve shape
            call GetData(            Me%BK1,                                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TBCONST1',   &
                                     default    = 0.05,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR101.' 




            !BK2, constant to control temperature response curve shape
            call GetData(            Me%BK2,                                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TBCONST2',   &
                                     default    = 0.98,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR102.' 


            !BK3, constant to control temperature response curve shape
            call GetData(            Me%BK3,                                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TBCONST3',   &
                                     default    = 0.98,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR103.' 


            !BK4, constant to control temperature response curve shape
            call GetData(            Me%BK4,                                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'TBCONST4',   &
                                     default    = 0.02,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR104.' 


            !Reads the rates & constants to the Ciliate simulation-----------------

            !CiliateGrazBactMin, minimum ciliates concentration for the existence of grazing, mgC/l
            call GetData(            Me%CiliateGrazBactMin,                         &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'GRAZBACMIN', &
                                     default    = 0.0045,                            &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR105.' 

            !CiliateGrazPhytoMin, minimum phytoplankton concentration for the existence of grazing, mgC/l
            call GetData(            Me%CiliateGrazPhytoMin,                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'CILGRAZPHYMIN', &
                                     default    = 0.0045,                            &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR105.' 

            !CiliateGrazPreyMin, minimum phytoplankton concentration for the existence of grazing, mgC/l
            call GetData(            Me%CiliateGrazPreyMin,                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'CILPREYMIN', &
                                     default    = 0.0045,                            &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR105.' 


            !CiliateReferenceRespirationRate, rate of consumption of Carbon by respiration and 
            call GetData(            Me%CiliateReferenceRespirationRate,            &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'CREFRESP',   &
                                     default    = 0.02,                            & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR106.' 


            !CiliateExcretionFactor, excretion constant
            call GetData(            Me%CiliateExcretionFactor,                     &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'CEXCFAC',    &
                                     default    = 0.02,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR107.' 


            !CiliateExcretionConstant, excretion constant
            call GetData(           Me%CiliateExcretionConst,                       &
                                    Me%ObjEnterData, flag,                          &
                                    SearchType = FromFile, keyword = 'CEXCCONS',    &
                                    default    = 1.0305,                            &
                                    ClientModule = 'ModuleWaterQuality',            &
                                    STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR108.' 


            !CiliateMortalityCoef    
            call GetData(            Me%CiliateMortalityCoef,                       &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'MORTCICOEF', &
                                     default    = 0.0,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR109.' 


            !CiliateMinMortalityRate  
            call GetData(            Me%CiliateMinMortalityRate,                    &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'MINMORTCI',  &
                                     default    = 0.0,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR110.' 


            !CiliateMaxMortalityRate
            call GetData(            Me%CiliateMaxMortalityRate,                    &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'MAXMORTCI',  &
                                     default    = 0.044,                            & !1/day
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR111.' 
            
            !CiliateIngestionConst, 1/2 sat
           call GetData(             Me%CiliateIngestionConst,                      &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'INGCONSC',   &
                                     default    = 0.85,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR112.' 


            !CiliateEfficiencyCaptureBacteria   
             call GetData( Me%CiliateEfficiencyCaptureBact,                         &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'CILEFFCAPBA',   &
                                     default    = 0.5,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR113.' 
                                
            
            !CiliateEfficiencyCapturePhyto  
             call GetData( Me%CiliateEfficiencyCapturePhyto,                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'CILEFFCAPPHY',  &
                                     default    = 0.5,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR113.'

            !CiliateIngestionMax    
            call GetData(            Me%CiliateIngestionMax,                        &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'CINGMAX',    &
                                     default    = 2.,                             &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR114.' 
                                    
            
            !CiliateAssimilationBacteriaRate    
            call GetData(            Me%CiliateAssimilationBacteriaRate,            &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile, keyword = 'CILBACASS',  &
                                     default    = 0.5,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR115.' 

            !CiliateAssimilationPhytoRate    
            call GetData(            Me%CiliateAssimilationPhytoRate,               &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,                         &
                                     keyword = 'CILPHYASS',                         &
                                     default    = 0.5,                              &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
           
           
            !RatioOxygenCarbonCilRespiration, Ciliates respiration Oxygen:Carbon ratio, 
            !mgO2 / mgC
            call GetData(            Me%RatioOxygenCarbonCilRespiration,            &
                                     Me%ObjEnterData, flag,                         &
                                     SearchType = FromFile,keyword = 'CILOCRATIO',    &
                                     default    = 32.0 / 12.0,                      &
                                     ClientModule = 'ModuleWaterQuality',           &
                                     STAT       = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR78.' 


            !
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR115.' 

                !PhytoRatioIngestionZoo 
                call GetData(            Me%BactRatioIngestionCiliates,                  &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword = 'BACINGCIL',&
                                         default    = 0.5,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR673.' 


                !PhytoRatioIngestionZoo 
                call GetData(            Me%PhytoRatioIngestionCiliates,                &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile, keyword = 'PHYINGCIL',&
                                         default    = 0.5,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR673.' 

        !Reads the rates & constants to the Larvae simulation-------------------------

                ! Awg Growth coefficient dependent of larvae weight 
                call GetData(            Me%Awg,                                        &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'AWG',                               &
                                         default    = 0.1344,                           &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR674.' 

                ! Bwg Growth coefficient dependent of larvae weight 
                call GetData(            Me%Bwg,                                        &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'BWG',                               &
                                         default    = 0.061,                            &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR674.' 

                ! Awz Death coefficient dependent of larvae weight 
                call GetData(            Me%Awz,                                        &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'AWZ',                               &
                                         default    = 1.073,                            &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR675.' 

                ! Bwz Death coefficient dependent of larvae weight 
                call GetData(            Me%Bwz,                                        &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'BWZ',                               &
                                         default    = 0.353,                            &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR675.' 

                ! Atg Growth coefficient dependent of temperature 
                call GetData(            Me%Atg,                                        &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'ATG',                               &
                                         default    = 0.0511,                           &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR676.' 

                ! Btg Growth coefficient dependent of temperature
                call GetData(            Me%Btg,                                        &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'BTG',                               &
                                         default    = 0.0052,                           &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR676.' 

                ! Atz Death coefficient dependent of temperature 
                call GetData(            Me%Atz,                                        &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'ATZ',                               &
                                         default    = 0.0149,                           &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR677.' 

                ! Btz Death coefficient dependent of temperature
                call GetData(            Me%Btz,                                        &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'BTZ',                               &
                                         default    = 0.0129,                           &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR677.' 

                ! Larvae shape factor 
                call GetData(            Me%Lshape,                                     &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'LSHAPE',                            &
                                         default    = 3.78,                             &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR678.' 

                ! Larvae density factor
                call GetData(            Me%Ldensity,                                   &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'LDENSITY',                          &
                                         default    = 2.95,                             &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR678.' 

                ! Number of larvae phases (valid values are 1 and 2) 
                call GetData(            Me%NPhases,                                    &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'NPHASES',                           &
                                         default    = 2,                                &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR679.' 

                ! Larvae Inital Age (days) 
                call GetData(            Me%Init_age,                                   &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'INIT_AGE',                          &
                                         default    = 5.,                               &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR679.' 

                ! Larvae Intermediate Age (days) 
                call GetData(            Me%Inter_age,                                  &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'INTER_AGE',                         &
                                         default    = 11.,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR679.' 

                ! Larvae Final Age (days) 
                call GetData(            Me%Final_age,                                  &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'FINAL_AGE',                         &
                                         default    = 46.,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR679.' 

                ! Larvae Inital Length (mm)
                call GetData(            Me%Init_length,                                &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'INIT_LENGTH',                       &
                                         default    = 5.,                               &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR680.' 

                ! Larvae Intermediate Length (mm)
                call GetData(            Me%Inter_length,                               &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'INTER_LENGTH',                      &
                                         default    = 10.,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR680.' 

                ! Larvae Final Length (mm)
                call GetData(            Me%Final_length,                               &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'FINAL_LENGTH',                      &
                                         default    = 35.,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR680.' 

                ! Reference food availability (mg/m3)
                call GetData(            Me%FishFood_ref,                               &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'FISHFOOD_REF',                      &
                                         default    = 0.5,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR681.' 

                ! Reference temperature (ºC)
                call GetData(            Me%Temperature_ref,                            &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'TEMPERATURE_REF',                   &
                                         default    = 15.,                              &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR682.' 

                ! Afg Growth coefficient dependent of fishfood availability (mg/m3) HalfSaturationConstant
                call GetData(            Me%Afg,                                        &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'AFG',                               &
                                         default    = 0.05,                             &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR682.' 
                    
        
        
        !Reads the rates & constants to the POM Pools  simulation---------------
        
                call GetData(            Me%POMIngestVmax,                              &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'POMINGESTVMAX',                     &
                                         default    = 0.05,                             &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR700.'
                
                call GetData(            Me%PONIngestKs,                                &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'PONINGESTKS',                       &
                                         default    = 0.05,                             &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR701.'
                
                call GetData(            Me%POPIngestKs,                                &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'POPINGESTKS',                       &
                                         default    = 0.05,                             &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR702.'
                
                call GetData(            Me%PON_CNratio,                                &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'PONCNRATIO',                        &
                                         default    = 0.05,                             &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR703.'    
        
                call GetData(            Me%PON_CPratio,                                &
                                         Me%ObjEnterData, flag,                         &
                                         SearchType = FromFile,                         &
                                         keyword = 'PONCPRATIO',                        &
                                         default    = 0.05,                             &
                                         ClientModule = 'ModuleWaterQuality',           &
                                         STAT       = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine WQReadFileConstants - ModuleWaterQuality. ERR704.'
        

        !aqui
        !Silica
        if (Me%PropCalc%Silica) then

            call  WQReadSilicaFileConstants                                         
 
        endif       !Silica

        !Diatoms
        if (Me%PropCalc%Diatoms) then
            call  WQReadDiatomsFileConstants                                        
        endif     !Diatoms

        !aqui

        !----------------------------------------------------------------------

    end subroutine WQReadFileConstants

    !-------------------------------------------------------------------------------

    subroutine WQReadSilicaFileConstants

        !External--------------------------------------------------------------
        integer                         :: FromFile
        integer                         :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                         :: flag

        !Begin----------------------------------------------------------------------

        call GetExtractType    (FromFile = FromFile)

        !Biogenic Silica Dissolution Rate in the water column, d-1
        call GetData(Me%SilicaCycle%KSiBiogenicDissRate,                            &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'SIKDISS',                                      &
                     default      = 0.01,                                           &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadSilicaFileConstants - Module ModuleWaterQuality - ERRO1'

        !Biogenic Silica Dissolution temperature cefficient, adim
        call GetData(Me%SilicaCycle%BiogenicDissTCoef,                              &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'SIDISSTCOEF',                                  &
                     default      = 0.01,                                           &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadSilicaFileConstants - Module ModuleWaterQuality - ERRO1'


     end subroutine WQReadSilicaFileConstants

    !-------------------------------------------------------------------------------

    subroutine WQReadDiatomsFileConstants

        !External--------------------------------------------------------------
        integer                         :: FromFile
        integer                         :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                         :: flag

        !Begin----------------------------------------------------------------------

        call GetExtractType    (FromFile = FromFile)    
        
        !Diatoms Maximum Gross Growth Rate, d-1 , EPA 1985
        call GetData(Me%Diatoms%DiaGrowMaxRate,                                     &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIGROWMAX',                                    &
                     default      = 3.0     ,                                       &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERRO1'

        !Diatoms Endogenous Respiration Constant, d-1
        call GetData(Me%Diatoms%DiaEndogRepConst,                                   &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIFENDREPC',                                   &
                     default      = 0.0175,                                         &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERRO2'

        !Fraction of Diatoms photosynthesis which is oxidized by photorespiration, d-1
        call GetData(Me%Diatoms%DiaPhotorespFactor,                                 &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIPHOTORES',                                   &
                     default      = 0.0175,                                        &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERRO3'

        !Diatoms Excretion Constant, adim
        call GetData(Me%Diatoms%DiaExcretionConstant,                               &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIEXCRCONS',                                   &
                     default      = 0.07,                                           &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERRO4'

        !Diatoms Maximum Mortality Rate, d-1
        call GetData(Me%Diatoms%DiaMortMaxRate,                                     &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIMORTMAX',                                    &
                     default      = 0.03,                                           &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERRO5'


        !Diatoms half-saturation mortality rate, d-1
        call GetData(Me%Diatoms%DiaMortSatConst,                                    &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIMORTCON',                                    &
                     default      = 0.3,                                            &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERRO6'

        !Zooplankton assimilation efficiency for diatoms, adim
        call GetData(Me%Diatoms%DiaE,                                               &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIASS_EFIC',                                   &
                     default      = 0.8,                                            &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERRO7'

        !Nitrogen half-saturation constant for Diatoms, mgN/l
        call GetData(Me%Diatoms%DiaNSatConst,                                       &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DINSATCONS',                                   &
                     default      = 0.015,                                          &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERRO8'

        !Phosphorus half-saturation constant for Diatoms, mgP/l
        call GetData(Me%Diatoms%DiaPSatConst,                                       &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIPSATCONS',                                   &
                     default      = 0.002,                                          &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERRO9'

        !Silica half-saturation constant for Diatoms, mgSi/l
        call GetData(Me%Diatoms%DiaSiSatConst,                                      &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DISISATCONS',                                  &
                     default      = 0.08,                                          &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR10'

        !Optimum light intensity for Diatoms photosynthesis, W/m2
        
        call GetData(Me%Diatoms%DiaPhotoinhibition,                                 &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIPHOTOIN',                                    &
                     default      = 121.,                                           &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR11'

        
        !Minimum temperature of optimal interval for Diatoms photosynthesis, ºC
        call GetData(Me%Diatoms%DiaTOptMin,                                         &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DITOPTMIN',                                    &
                     default      = 25.,                                            &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR12'

        !Maximum temperature of optimal interval for Diatoms photosynthesis, ºC
        call GetData(Me%Diatoms%DiaTOptMax,                                         &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DITOPTMAX',                                    &
                     default      = 26.5,                                           &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR13'

        !Minimum tolerable temperature for Diatoms Growth, ºC
        call GetData(Me%Diatoms%DiaTMin,                                            &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DITMIN',                                       &
                     default      = 4.,                                             &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR14'

        !Maximum tolerable temperature for Diatoms Growth, ºC
        call GetData(Me%Diatoms%DiaTMax,                                            &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DITMAX',                                       &
                     default      = 37.,                                            &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR15'

        !Constant to control Diatoms temperature response curve shape, adim
        call GetData(Me%Diatoms%DiaK1,                                              &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DITCONST1',                                    &
                     default      = 0.05,                                           &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR16'

        !Constant to control Diatoms temperature response curve shape, adim
        call GetData(Me%Diatoms%DiaK2,                                              &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DITCONST2',                                    &
                     default      = 0.98,                                           &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR17'

    
        !Constant to control Diatoms temperature response curve shape, adim
        call GetData(Me%Diatoms%DiaK3,                                              &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DITCONST3',                                    &
                     default      = 0.98,                                           &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR18'

        !Constant to control Diatoms temperature response curve shape, adim
        call GetData(Me%Diatoms%DiaK4,                                              &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DITCONST4',                                    &
                     default      = 0.02,                                           &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR19'

        !Diatoms Nitrogen/Carbon Ratio, mg N/mgC
        call GetData(Me%Diatoms%DiaAlfaNC,                                          &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIRATIONC',                                    &
                     default      = 0.18,                                           &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR20'

        !Diatoms Phosphorus/Carbon Ratio, mg P/mgC
        call GetData(Me%Diatoms%DiaAlfaPC,                                          &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIRATIOPC',                                    &
                     default      = 0.024,                                          &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR21'

        !Diatoms Silica/Carbon Ratio, mg Si/mgC
        call GetData(Me%Diatoms%DiaAlfaSiC,                                         &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIRATIOSIC',                                   &
                     default      = 0.6,                                            &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR22'


        !Fraction of soluble inorganic material excreted by Diatoms, adim
        call GetData(Me%Diatoms%DiaSolublInorgExcreFraction,                        &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DISOLEXCR',                                    &
                     default      = 0.4,                                           &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR23'

        !Fraction of dissolved organic material excreted by Diatoms, adim
        call GetData(Me%Diatoms%DiaExcreDissOrgFraction,                            &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIDISSDON',                                    &
                     default      = 0.5,                                           &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR24'

        !Diatoms minimum concentration for predation, mgC/l
        call GetData(Me%Diatoms%GrazDiaMin,                                         &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIGRAZMIN' ,                                   &
                     default      = 0.0045,                                          &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR25'

        !Proportion of diatoms in mesozooplankton ingestion, adim
        call GetData(Me%Diatoms%DiaRatioIngestionZoo,                               &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIRATINGZOO',                                  &
                     default      = 0.3,                                            &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR26'

        !Assimilation Coefficient o Diatoms by mesozooplankton
        call GetData(Me%Diatoms%DiaZooAssimilationRate,                             &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIZOASS',                                      &
                     default      = 0.8,                                            &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR27'

        !Zooplankton Efficiency Capture of Diatoms
        call GetData(Me%Diatoms%ZooEfficiencyCaptureDiatoms,                        &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DIZOOEFFCAP',                                  &
                     default      = 0.8,                                            &
                     ClientModule = 'ModuleWaterQuality',                           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WQReadDiatomsFileConstants - Module ModuleWaterQuality - ERR28'

        
     end subroutine WQReadDiatomsFileConstants
    
    !----------------------------------------------------------------------------
    
    subroutine AllocateVariables

        !External----------------------------------------------------------------
        integer :: STAT_CALL

        !Local-------------------------------------------------------------------
        integer :: PropLB, PropUB

        !------------------------------------------------------------------------

        PropLB    = Me%Prop%ILB
        PropUB    = Me%Prop%IUB

        allocate(Me%Matrix (PropLB:PropUB, PropLB:PropUB))
        allocate(Me%IndTerm(PropLB:PropUB               ))


cd1 :   if (Me%CalcMethod%ExplicitMethod) then
            allocate(Me%NewMass(PropLB:PropUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'Subroutine AllocateVariables - ModuleWaterQuality. ERR01.'
        end if cd1

        !------------------------------------------------------------------------

    end subroutine AllocateVariables   

    !----------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetNCRatio(WaterQualityID, Property, Ratio, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: WaterQualityID
        real                                :: Ratio
        integer                             :: Property
        integer, optional, intent(OUT)      :: STAT    

        !External--------------------------------------------------------------
        integer                         :: ready_

        !Local-----------------------------------------------------------------
        integer                         :: STAT_
        real                                :: AlphaNTS

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WaterQualityID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !Ratio of gN in g dry substance according to Redfiel: 16mol*14gN/mol / 2749 gTS
            AlphaNTS = 16.*14./2749.
            Ratio = 1.0/AlphaNTS  !for all non-living VSS as PON and PONr
            if (Property == Me%PropIndex%Phyto                  ) Ratio = Me%AlfaPhytoNC   /AlphaNTS
            if (Property == Me%PropIndex%Diatoms                ) Ratio = Me%OMAlfaNC      /AlphaNTS
            if (Property == Me%PropIndex%Zoo                    ) Ratio = Me%AlfaZooNC     /AlphaNTS
            if (Property == Me%PropIndex%Ciliate                ) Ratio = Me%AlfaCilNC     /AlphaNTS
            if (Property == Me%PropIndex%Bacteria               ) Ratio = Me%AlfaBacteriaNC/AlphaNTS

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1
        
        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetNCRatio

    !--------------------------------------------------------------------------


    subroutine GetWaterQualitySize(WaterQualityID, PropLB, PropUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: WaterQualityID
        integer, optional, intent(OUT)  :: PropLB,    PropUB
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------

        integer                         :: ready_              

        !Local-----------------------------------------------------------------
        integer                         :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WaterQualityID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PropLB   )) PropLB    = Me%Prop%ILB
            if (present(PropUB   )) PropUB    = Me%Prop%IUB

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetWaterQualitySize

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine GetWQOptions(WaterQualityID, Zoo,                              &
                                           Larvae,                            &
                                           Age,                               & ! Aires
                                           Phyto,                             &
                                           Nitrogen,                          & 
                                           Phosphorus,                        &
                                           Oxygen,                            &
                                           Salinity,                          &  
                                           BOD,                               &
                                           Bacteria,                          &
                                           Ciliate,                           &
                                           Diatoms,                           &  !aqui
                                           Silica,                            &  !aqui
                                           PomPools,                          &
                                           ExplicitMethod,                    &
                                           ImplicitMethod,                    &
                                           SemiImpMethod, STAT) 

        !Arguments-------------------------------------------------------------
        integer                         :: WaterQualityID
        integer, optional, intent(OUT)  :: STAT
        logical, optional, intent(OUT)  :: Zoo,  Age, Larvae, Phyto, Nitrogen, Phosphorus, Oxygen,   &
                                           Salinity, BOD, Bacteria, Ciliate, Diatoms, Silica, Pompools 
        logical, optional, intent(OUT)  :: ExplicitMethod, ImplicitMethod, SemiImpMethod  

        !External--------------------------------------------------------------
        integer                         :: ready_              

        !Local-----------------------------------------------------------------
        integer                         :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WaterQualityID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(Zoo           )) Zoo            = Me%PropCalc%Zoo
            if (present(Larvae        )) Larvae         = Me%PropCalc%Larvae 
            if (present(Age           )) Age            = Me%PropCalc%Age    
            if (present(Phyto         )) Phyto          = Me%PropCalc%Phyto
            if (present(Nitrogen      )) Nitrogen       = Me%PropCalc%Nitrogen
            if (present(Phosphorus    )) Phosphorus     = Me%PropCalc%Phosphorus
            if (present(Oxygen        )) Oxygen         = Me%PropCalc%Oxygen       
            if (present(Salinity      )) Salinity       = Me%PropCalc%Salinity                            
            if (present(BOD           )) BOD            = Me%PropCalc%BOD 
            if (present(Bacteria      )) Bacteria       = Me%PropCalc%Bacteria         ! so
            if (present(Ciliate       )) Ciliate        = Me%PropCalc%Ciliate  
            if (present(Diatoms       )) Diatoms        = Me%PropCalc%Diatoms        !aqui
            if (present(Silica        )) Silica         = Me%PropCalc%Silica         !aqui
            if (present(Pompools      )) Pompools       = Me%PropCalc%Pompools
            if (present(ExplicitMethod)) ExplicitMethod = Me%CalcMethod%ExplicitMethod
            if (present(ImplicitMethod)) ImplicitMethod = Me%CalcMethod%ImplicitMethod
            if (present(SemiImpMethod )) SemiImpMethod  = Me%CalcMethod%SemiImpMethod    

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetWQOptions

    !--------------------------------------------------------------------------

    subroutine GetWQPropIndex(WaterQualityID, Zoo,                             &
                                             Larvae,                           &  
                                             Age,                              &
                                             Phyto,                            &
                                             Ammonia,                          &
                                             Nitrate,                          &
                                             Nitrite,                          &
                                             DissOrganicNitrogenRefractory,    &
                                             DONNonRefractory,                 &
                                             PartOrganicNitrogen,              &
                                             PartOrganicNitrogenRefractory,    &
                                             PONitrogen1,                      &
                                             PONitrogen2,                      &
                                             PONitrogen3,                      &
                                             PONitrogen4,                      &
                                             PONitrogen5,                      &
                                             Oxygen,                           &
                                             BOD,                              &
                                             Bacteria,                         &
                                             Ciliate,                          &
                                             Diatoms,                          &  !aqui
                                             BiogenicSilica,                   &  !aqui
                                             DissolvedSilica,                  &  !aqui
                                             DissOrganicPhosphorusRefractory,  &
                                             DOPNonRefractory,                 &
                                             PartOrganicPhosphorus,            &
                                             POPhosphorus1,                    &
                                             POPhosphorus2,                    &
                                             POPhosphorus3,                    &
                                             POPhosphorus4,                    &
                                             POPhosphorus5,                    &
                                             InorganicPhosphorus,              &
                                             STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: WaterQualityID
                                        
        integer, optional, intent(OUT)  :: STAT
                                        
        integer, optional, intent(OUT)  :: Zoo                              
        integer, optional, intent(OUT)  :: Larvae                               
        integer, optional, intent(OUT)  :: Age                             
        integer, optional, intent(OUT)  :: Phyto                            
        integer, optional, intent(OUT)  :: Ammonia                          
        integer, optional, intent(OUT)  :: Nitrate                          
        integer, optional, intent(OUT)  :: Nitrite                          
        integer, optional, intent(OUT)  :: DissOrganicNitrogenRefractory    
        integer, optional, intent(OUT)  :: DONNonRefractory 
        integer, optional, intent(OUT)  :: PartOrganicNitrogen
        integer, optional, intent(OUT)  :: PONitrogen1
        integer, optional, intent(OUT)  :: PONitrogen2
        integer, optional, intent(OUT)  :: PONitrogen3
        integer, optional, intent(OUT)  :: PONitrogen4
        integer, optional, intent(OUT)  :: PONitrogen5
        integer, optional, intent(OUT)  :: POPhosphorus1
        integer, optional, intent(OUT)  :: POPhosphorus2 
        integer, optional, intent(OUT)  :: POPhosphorus3 
        integer, optional, intent(OUT)  :: POPhosphorus4 
        integer, optional, intent(OUT)  :: POPhosphorus5               
        integer, optional, intent(OUT)  :: PartOrganicNitrogenRefractory              
        integer, optional, intent(OUT)  :: Oxygen                           
        integer, optional, intent(OUT)  :: BOD    
        integer, optional, intent(OUT)  :: Bacteria
        integer, optional, intent(OUT)  :: Ciliate
        integer, optional, intent(OUT)  :: Diatoms          !aqui
        integer, optional, intent(OUT)  :: BiogenicSilica   !aqui
        integer, optional, intent(OUT)  :: DissolvedSilica  !aqui                         
        integer, optional, intent(OUT)  :: DissOrganicPhosphorusRefractory    
        integer, optional, intent(OUT)  :: DOPNonRefractory 
        integer, optional, intent(OUT)  :: PartOrganicPhosphorus              
        integer, optional, intent(OUT)  :: InorganicPhosphorus              

        !External--------------------------------------------------------------
        integer :: ready_              

        !Local-----------------------------------------------------------------
        integer                         :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WaterQualityID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            if (present(Zoo                  )) Zoo                   = Me%PropIndex%Zoo
            if (present(Larvae               )) Larvae                = Me%PropIndex%Larvae 
            if (present(Age                  )) Age                   = Me%PropIndex%Age     
            if (present(Phyto                )) Phyto                 = Me%PropIndex%Phyto
            if (present(Ammonia              )) Ammonia               = Me%PropIndex%Ammonia
            if (present(Nitrate              )) Nitrate               = Me%PropIndex%Nitrate
            if (present(Nitrite              )) Nitrite               = Me%PropIndex%Nitrite
            if (present(PartOrganicNitrogen  )) PartOrganicNitrogen   = Me%PropIndex%PartOrganicNitrogen
            if (present(PartOrganicNitrogenRefractory  ))   &
                PartOrganicNitrogenRefractory   = Me%PropIndex%PartOrganicNitrogenRefractory             
            if (present(PONitrogen1          )) PONitrogen1           = Me%PropIndex%PONitrogen1    
            if (present(PONitrogen2          )) PONitrogen2           = Me%PropIndex%PONitrogen2
            if (present(PONitrogen3          )) PONitrogen3           = Me%PropIndex%PONitrogen3
            if (present(PONitrogen4          )) PONitrogen4           = Me%PropIndex%PONitrogen4
            if (present(PONitrogen5          )) PONitrogen5           = Me%PropIndex%PONitrogen5
            if (present(POPhosphorus1        )) POPhosphorus1         = Me%PropIndex%POPhosphorus1    
            if (present(POPhosphorus2        )) POPhosphorus2         = Me%PropIndex%POPhosphorus2
            if (present(POPhosphorus3        )) POPhosphorus3         = Me%PropIndex%POPhosphorus3
            if (present(POPhosphorus4        )) POPhosphorus4         = Me%PropIndex%POPhosphorus4
            if (present(POPhosphorus5        )) POPhosphorus5         = Me%PropIndex%POPhosphorus5       
            if (present(Oxygen               )) Oxygen                = Me%PropIndex%Oxygen
            if (present(BOD                  )) BOD                   = Me%PropIndex%BOD 
            if (present(Bacteria             )) Bacteria              = Me%PropIndex%Bacteria      
            if (present(Ciliate              )) Ciliate               = Me%PropIndex%Ciliate
            if (present(Diatoms              )) Diatoms               = Me%PropIndex%Diatoms          !aqui
            if (present(BiogenicSilica       )) BiogenicSilica        = Me%PropIndex%BiogenicSilica   !aqui
            if (present(DissolvedSilica      )) DissolvedSilica       = Me%PropIndex%DissolvedSilica  !aqui      
            if (present(PartOrganicPhosphorus)) PartOrganicPhosphorus = Me%PropIndex%PartOrganicPhosphorus
            if (present(InorganicPhosphorus  )) InorganicPhosphorus   = Me%PropIndex%InorganicPhosphorus

            if (present(DissOrganicNitrogenRefractory   ))  &
                DissOrganicNitrogenRefractory    = Me%PropIndex%DissOrganicNitrogenRefractory
            if (present(DONNonRefractory))                  &
                DONNonRefractory = Me%PropIndex%DONNonRefractory


            if (present(DissOrganicPhosphorusRefractory )) &
                DissOrganicPhosphorusRefractory    = Me%PropIndex%DissOrganicPhosphorusRefractory
            if (present(DOPNonRefractory))                  &
                DOPNonRefractory = Me%PropIndex%DOPNonRefractory


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetWQPropIndex

    !--------------------------------------------------------------------------

    subroutine GetDTWQM(WaterQualityID, DTDay, DTSecond, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: WaterQualityID
        real,    optional, intent(OUT)  :: DTDay
        real,    optional, intent(OUT)  :: DTSecond
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                         :: ready_              

        !Local-----------------------------------------------------------------
        integer                         :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WaterQualityID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(DTDay   )) DTDay    = Me%DTDay
            if (present(DTSecond)) DTSecond = Me%DTSecond

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDTWQM

    !--------------------------------------------------------------------------
 
    subroutine GetWQPropRateFlux(WaterQualityID, Firstprop, Secondprop, PropRateFlux, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: WaterQualityID
        real, dimension(:),     pointer     :: PropRateFlux
        integer                             :: Firstprop,Secondprop
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer         :: ready_        

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
        integer                             :: prop,equa
        logical                             :: found       
        type(T_EquaRateFlux),       pointer :: EquaRateFluxX
        type(T_PropRateFlux),   pointer     :: PropRateFluxX

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WaterQualityID, ready_)    
        
        found=.FALSE.

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            call Read_Lock(mWATERQUALITY_, Me%InstanceID)


                EquaRateFluxX => Me%FirstEquaRateFlux

do1 :           do while(associated(EquaRateFluxX))  
                    
                    PropRateFluxX => EquaRateFluxX%FirstPropRateFlux
                    
                    do while(associated(PropRateFluxX))  
                    
                        equa = EquaRateFluxX%ID
                        prop = PropRateFluxX%ID
                   
                        if(Prop.eq.Firstprop.and.Equa.eq.SecondProp) then

                            PropRateFlux => PropRateFluxX%Field (:)

                            found=.true.
                            exit do1

                        endif
                     
                        PropRateFluxX => PropRateFluxX%Next

                    end do
                 
                    EquaRateFluxX => EquaRateFluxX%Next
          
                end do do1


               !aqui_3
                if (.NOT.found) then

                    if (FirstProp.eq. Me%GrossProduction%ID) then

                        PropRateFlux => Me%GrossProduction%Field
                        found=.TRUE.

                    elseif (FirstProp.eq. Me%TempLimitation%ID) then

                        PropRateFlux => Me%TempLimitation%Field
                        found=.TRUE.

                    elseif (FirstProp.eq. Me%NutLimitation%ID) then

                        PropRateFlux => Me%NutLimitation%Field
                        found=.TRUE.

                    elseif (FirstProp.eq. Me%NLimitation%ID) then

                        PropRateFlux => Me%NLimitation%Field
                        found=.TRUE.

                    elseif (FirstProp.eq. Me%PLimitation%ID) then

                        PropRateFlux => Me%PLimitation%Field
                        found=.TRUE.

                    elseif (FirstProp.eq. Me%LightLimitation%ID) then

                        PropRateFlux => Me%LightLimitation%Field
                        found=.TRUE.
                    
                    elseif (FirstProp.eq. Me%Diatoms%DiaGrossProduction%ID) then

                        PropRateFlux => Me%Diatoms%DiaGrossProduction%Field
                        found=.TRUE.

                    elseif (FirstProp.eq. Me%Diatoms%DiaTempLimitation%ID) then

                        PropRateFlux => Me%Diatoms%DiaTempLimitation%Field
                        found=.TRUE.

                    elseif (FirstProp.eq. Me%Diatoms%DiaNutLimitation%ID) then

                        PropRateFlux => Me%Diatoms%DiaNutLimitation%Field
                        found=.TRUE.

                    elseif (FirstProp.eq. Me%Diatoms%DiaNLimitation%ID) then

                        PropRateFlux => Me%Diatoms%DiaNLimitation%Field
                        found=.TRUE.

                    elseif (FirstProp.eq. Me%Diatoms%DiaPLimitation%ID) then

                        PropRateFlux => Me%Diatoms%DiaPLimitation%Field
                        found=.TRUE.

                    elseif (FirstProp.eq. Me%Diatoms%DiaSiLimitation%ID) then

                        PropRateFlux => Me%Diatoms%DiaSiLimitation%Field
                        found=.TRUE.

                    elseif (FirstProp.eq. Me%Diatoms%DiaLightLimitation%ID) then

                        PropRateFlux => Me%Diatoms%DiaLightLimitation%Field
                        found=.TRUE.

                    endif

                endif

                nullify  (PropRateFluxX,EquaRateFluxX)
           
            if (found) then              
                STAT_ = SUCCESS_
            else
                STAT_ = NOT_FOUND_ERR_
            endif
        else cd1
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetWQPropRateFlux

    !--------------------------------------------------------------------------

    
    
    !--------------------------------------------------------------------------
    subroutine UnGetWQPropRateFlux(WaterQualityID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: WaterQualityID
        integer, optional, intent(OUT)  :: STAT
        real, pointer, dimension(:)     :: Array

        !External--------------------------------------------------------------
        integer                         :: ready_   

        !Local-----------------------------------------------------------------
        integer                         :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WaterQualityID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_Unlock(mWATERQUALITY_, Me%InstanceID, "UnGetWQPropRateFlux")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UnGetWQPropRateFlux

    !--------------------------------------------------------------------------




    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    
    subroutine WaterQuality(WaterQualityID,                                     &
                            Salinity,                                           &
                            Temperature,                                        &
                            ShortWaveRadiation,                                 &
                            LightExtCoefField,                                  &
                            Thickness,                                          &
                            Mass,                                               &
                            WQArrayLB, WQArrayUB,                               &
                            OpenPoints,                                         &               
                            FishFood,                                           &
                            STAT)  

        !Arguments---------------------------------------------------------------
        integer                                       :: WaterQualityID
        real,                 pointer, dimension(:  ) :: Salinity
        real,                 pointer, dimension(:  ) :: Temperature
        real,                 pointer, dimension(:  ) :: ShortWaveRadiation
        real,                 pointer, dimension(:  ) :: LightExtCoefField
        real,                 pointer, dimension(:  ) :: Thickness
        real,                 pointer, dimension(:,:) :: Mass
        integer, optional,    pointer, dimension(:  ) :: OpenPoints
        real,    optional,    pointer, dimension(:  ) :: FishFood
        integer,              intent(IN )             :: WQArrayLB, WQArrayUB  
        integer, optional,    intent(OUT)             :: STAT
         
        !External----------------------------------------------------------------
        integer                                       :: index
        integer                                       :: ready_   
               
        !Local-------------------------------------------------------------------
        integer                                       :: STAT_          
        logical                                       :: CalcPoint

        !real                                          :: totalN
        !real                                          :: totalP
        !real                                          :: totalSi
        !------------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(WaterQualityID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            Me%ExternalVar%Salinity                   => Salinity
            if (.NOT. associated(Me%ExternalVar%Salinity))         &
                stop 'Subroutine WaterQuality - ModuleWaterQuality. ERR01' 


            Me%ExternalVar%Temperature                => Temperature
            if (.NOT. associated(Me%ExternalVar%Temperature))        &
                stop 'Subroutine WaterQuality - ModuleWaterQuality. ERR02'
                
            Me%ExternalVar%ShortWaveRadiation         => ShortWaveRadiation
            if (.NOT. associated(Me%ExternalVar%ShortWaveRadiation)) &
                stop 'Subroutine WaterQuality - ModuleWaterQuality. ERR02' 
            
            Me%ExternalVar%LightExtCoefField          => LightExtCoefField
            if (.NOT. associated(Me%ExternalVar%LightExtCoefField))  &
                stop 'Subroutine WaterQuality - ModuleWaterQuality. ERR02' 
            
            Me%ExternalVar%Thickness                  => Thickness
            if (.NOT. associated(Me%ExternalVar%Thickness))          &
                stop 'Subroutine WaterQuality - ModuleWaterQuality. ERR02'  
             

            Me%ExternalVar%Mass                       => Mass
            if (.NOT. associated(Me%ExternalVar%Mass))               &
                stop 'Subroutine WaterQuality - ModuleWaterQuality. ERR03.'

cd6 :       if (present(FishFood)) then
                Me%ExternalVar%FishFood               => FishFood
                if (.NOT. associated(Me%ExternalVar%FishFood))       &
                    stop 'Subroutine WaterQuality - ModuleWaterQuality. ERR4.'
            else cd6
                nullify(Me%ExternalVar%FishFood)
            end if cd6


            call StartWaterQualityIteration

do1 :       do index = WQArrayLB, WQArrayUB
            
            !If this module is called from the WaterQuality3D module, OpenPoint is present
            !and the WQ module runs for all Openpoints
            !If this module is called from the Lagrangian module, OpenPoint is not present
            !and the WQ module runs for all volumes
            if (present(OpenPoints)) then
                if (OpenPoints(index) == OpenPoint) then
                    CalcPoint = .true.
                else
                    CalcPoint = .false.
                endif
            else
                CalcPoint = .true.
            endif


            if (CalcPoint) then
                call WQCoeficientsCalculation   (index)

                !The rates can just be calculated if the rate flux is associated
                !In the case that this module is used by the lagrangian module
                !the rate fluxes are not calculated
                !Rates must be computed before call to WQSystemResolution to use 
                !old concentrations 
                if (associated(Me%FirstEquaRateFlux)) then
                    call WQRatesCalculation      (index)
                end if

                call WQSystemResolution         (index)
                
        !Para fazer o balanço das propriedades        
        
        !If (Me%ExternalVar%Mass(Me%PropIndex%Nitrate, index) .lt.0) then

        !    totalN = Me%ExternalVar%Mass(Me%PropIndex%Nitrate, index)                              
        !    write (*,*) 'Nitrate_afterWQ', totalN, index
        !endif
        !
        !If (Me%ExternalVar%Mass(Me%PropIndex%DONNonRefractory, index) .lt.0) then                            
        !   write (*,*) 'DONnr_afterWQ', Me%ExternalVar%Mass(Me%PropIndex%DONNonRefractory, index), index
        !endif
        !
        !If (Me%ExternalVar%Mass(Me%PropIndex%Ammonia, index) .lt.0) then
        !
        !    totalN = Me%ExternalVar%Mass(Me%PropIndex%Ammonia, index)                              
        !   write (*,*) 'Ammonia_afterWQ', totalN, index
        !endif
        !
        !        
        !If (Me%PropCalc%Diatoms .and. (.not. Me%PropCalc%Phyto)) then
        !
        ! totalN = Me%ExternalVar%Mass(Me%PropIndex%Zoo, index)                              &
        !      + Me%ExternalVar%Mass(Me%PropIndex%Diatoms, index)                          
        !
        !
        ! write (*,*) 'Total C', totalN, index
        !end if
        !
        !
        !If (Me%PropCalc%Phyto .and. (.not. Me%PropCalc%Diatoms)) then
        !
        ! totalN = Me%ExternalVar%Mass(Me%PropIndex%Zoo, index)                              &
        !       + Me%ExternalVar%Mass(Me%PropIndex%Phyto, index)                          
        !
        !
        ! write (*,*) 'Total C', totalN, index
        !end if 
        !
        !
        !If (Me%PropCalc%Phyto .and. Me%PropCalc%Diatoms) then
        !
        ! totalN = Me%ExternalVar%Mass(Me%PropIndex%Zoo, index)                              &
        !       + Me%ExternalVar%Mass(Me%PropIndex%Phyto, index)                             &
        !       + Me%ExternalVar%Mass(Me%PropIndex%Diatoms, index)                          
        !
        !
        ! write (*,*) 'Total C', totalN, index
        !end if 
        !
        !
        !If (Me%PropCalc%Nitrogen) then 
        !
        !  totalN = Me%ExternalVar%Mass(Me%PropIndex%Ammonia, index)                              &
        !           + Me%ExternalVar%Mass(Me%PropIndex%Nitrate, index)                             &
        !           + Me%ExternalVar%Mass(Me%PropIndex%Nitrite, index)                             &
        !           + Me%ExternalVar%Mass(Me%PropIndex%PartOrganicNitrogen, index)                 &
        !           + Me%ExternalVar%Mass(Me%PropIndex%DONNonRefractory, index)                    &
        !           + Me%ExternalVar%Mass(Me%PropIndex%DissOrganicNitrogenRefractory, index)       &
        !            + Me%ExternalVar%Mass(Me%PropIndex%Zoo, index) * Me%AlfaZooNC     
        !           
        !     If (Me%PropCalc%Diatoms) then
        !       totalN = totalN &
        !              + Me%ExternalVar%Mass(Me%PropIndex%Diatoms, index) * Me%Diatoms%DiaAlfaNC
        !    endif
        !
        !     If (Me%PropCalc%Phyto) then
        !     
        !          totalN = totalN &
        !                 + Me%ExternalVar%Mass(Me%PropIndex%Phyto, index) * Me%AlfaPhytoNC
        !     endif
        !
        !
        !    write (*,*) 's/des...Depois..Total N', totalN, index
        !
        !end if 
        !
        !If (Me%PropCalc%Phosphorus) then     
        !    totalP = Me%ExternalVar%Mass(Me%PropIndex%PartOrganicPhosphorus, index)               &
        !           + Me%ExternalVar%Mass(Me%PropIndex%DissOrganicPhosphorusRefractory, index)     &
        !           + Me%ExternalVar%Mass(Me%PropIndex%DOPNonRefractory, index)                    &
        !           + Me%ExternalVar%Mass(Me%PropIndex%InorganicPhosphorus, index)                 &
        !           + Me%ExternalVar%Mass(Me%PropIndex%Zoo, index) * Me%AlfaZooPC                  
        !
        !   
        !    If (Me%PropCalc%Diatoms) then
        !
        !        totalP = totalP &
        !                 + Me%ExternalVar%Mass(Me%PropIndex%Diatoms, index) * Me%Diatoms%DiaAlfaPC
        !
        !    endif
        !
        !    If (Me%PropCalc%Phyto) then
        !
        !        totalP = totalP &
        !                 + Me%ExternalVar%Mass(Me%PropIndex%Phyto, index) * Me%AlfaPhytoPC
        !
        !    endif
        !
        !
        !    write (*,*) 's/des...Depois..Total P.         ', totalP, index 
        !
        !endif
        !
        !
        !If (Me%PropCalc%Silica) then
        !   
        !    totalSi = Me%ExternalVar%Mass(Me%PropIndex%BiogenicSilica , index)                    &
        !            + Me%ExternalVar%Mass(Me%PropIndex%DissolvedSilica, index)                    &
        !            + Me%ExternalVar%Mass(Me%PropIndex%Diatoms, index) * Me%Diatoms%DiaAlfaSiC                 
        !
        !    write (*,*) 's/des...Depois..TotalSi                    ',totalSi, index
        !
        !endif

        !!----------------------------------------------------------------------


            end if
            end do do1


            nullify(Me%ExternalVar%Salinity   )
            nullify(Me%ExternalVar%Temperature)
            nullify(Me%ExternalVar%Mass       )
            nullify(Me%ExternalVar%FishFood   )

           
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                      &
            STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine WaterQuality

    !----------------------------------------------------------------------------






!********************
!********* ver Sophie
!********************

    !----------------------------------------------------------------------------

    subroutine StartWaterQualityIteration

        !Local-------------------------------------------------------------------
        integer :: PropLB, PropUB
        integer :: i, j

        !------------------------------------------------------------------------

        propLB = Me%Prop%ILB 
        propUB = Me%Prop%IUB 


do1 :   do i = PropLB, PropUB
            Me%IndTerm(i) = 0.0

do2 :       do j = PropLB, PropUB
                Me%Matrix(i,j) = 0.0
            end do do2
        end do do1

        !------------------------------------------------------------------------

    end subroutine StartWaterQualityIteration

    !----------------------------------------------------------------------------




    !--------------------------------------------------------------------------

    subroutine WQCoeficientsCalculation(index)

        !Arguments-------------------------------------------------------------

        integer, intent(IN) :: index

        !real                                          :: totalN
        !real                                          :: totalP
        !real                                          :: totalSi

        
        
        !!Para fazer o balanço das propriedades        
        
        !If (Me%ExternalVar%Mass(Me%PropIndex%Nitrate, index) .lt.0) then
        !
        !    totalN = Me%ExternalVar%Mass(Me%PropIndex%Nitrate, index)                              
        !    write (*,*) 'Nitrate_beforeWQ', totalN, index
        !endif
        !
        !If (Me%ExternalVar%Mass(Me%PropIndex%Ammonia, index) .lt.0) then
        !
        !    totalN = Me%ExternalVar%Mass(Me%PropIndex%Ammonia, index)                              
        !    write (*,*) 'Ammonia_beforeWQ', totalN, index
        !endif
        !
        !
        !If (Me%PropCalc%Diatoms .and. (.not. Me%PropCalc%Phyto)) then
        
        ! totalN = Me%ExternalVar%Mass(Me%PropIndex%Zoo, index)                              &
        !      + Me%ExternalVar%Mass(Me%PropIndex%Diatoms, index)                          
        !
        !
        ! write (*,*) 'Total C', totalN, index
        !end if
        !
        !
        !If (Me%PropCalc%Phyto .and. (.not. Me%PropCalc%Diatoms)) then
        !
        ! totalN = Me%ExternalVar%Mass(Me%PropIndex%Zoo, index)                              &
        !       + Me%ExternalVar%Mass(Me%PropIndex%Phyto, index)                          
        !
        !
        ! write (*,*) 'Total C', totalN, index
        !end if 
        !
        !
        !If (Me%PropCalc%Phyto .and. Me%PropCalc%Diatoms) then
        !
        ! totalN = Me%ExternalVar%Mass(Me%PropIndex%Zoo, index)                              &
        !       + Me%ExternalVar%Mass(Me%PropIndex%Phyto, index)                             &
        !       + Me%ExternalVar%Mass(Me%PropIndex%Diatoms, index)                          
        !
        !
        ! write (*,*) 'Total C', totalN, index
        !end if 
        !
        !
        !If (Me%PropCalc%Nitrogen) then 
        !
        !   totalN = Me%ExternalVar%Mass(Me%PropIndex%Ammonia, index)                              &
        !           + Me%ExternalVar%Mass(Me%PropIndex%Nitrate, index)                             &
        !           + Me%ExternalVar%Mass(Me%PropIndex%Nitrite, index)                             &
        !           + Me%ExternalVar%Mass(Me%PropIndex%PartOrganicNitrogen, index)                 &
        !           + Me%ExternalVar%Mass(Me%PropIndex%DONNonRefractory, index)                    &
        !           + Me%ExternalVar%Mass(Me%PropIndex%DissOrganicNitrogenRefractory, index)       &
        !            + Me%ExternalVar%Mass(Me%PropIndex%Zoo, index) * Me%AlfaZooNC     
        !           
        !     If (Me%PropCalc%Diatoms) then
        !       totalN = totalN &
        !              + Me%ExternalVar%Mass(Me%PropIndex%Diatoms, index) * Me%Diatoms%DiaAlfaNC
        !    endif
        !
        !     If (Me%PropCalc%Phyto) then
        !     
        !          totalN = totalN &
        !                 + Me%ExternalVar%Mass(Me%PropIndex%Phyto, index) * Me%AlfaPhytoNC
        !     endif
        !
        !
        !    write (*,*) 's/des...Depois..Total N', totalN, index
        !
        !end if 
        !
        !If (Me%PropCalc%Phosphorus) then     
        !    totalP = Me%ExternalVar%Mass(Me%PropIndex%PartOrganicPhosphorus, index)               &
        !           + Me%ExternalVar%Mass(Me%PropIndex%DissOrganicPhosphorusRefractory, index)     &
        !           + Me%ExternalVar%Mass(Me%PropIndex%DOPNonRefractory, index)                    &
        !           + Me%ExternalVar%Mass(Me%PropIndex%InorganicPhosphorus, index)                 &
        !           + Me%ExternalVar%Mass(Me%PropIndex%Zoo, index) * Me%AlfaZooPC                  
        !
        !   
        !    If (Me%PropCalc%Diatoms) then
        !
        !        totalP = totalP &
        !                 + Me%ExternalVar%Mass(Me%PropIndex%Diatoms, index) * Me%Diatoms%DiaAlfaPC
        !
        !    endif
        !
        !    If (Me%PropCalc%Phyto) then
        !
        !        totalP = totalP &
        !                 + Me%ExternalVar%Mass(Me%PropIndex%Phyto, index) * Me%AlfaPhytoPC
        !
        !    endif
        !
        !
        !    write (*,*) 's/des...Depois..Total P.         ', totalP, index 
        !
        !endif
        !
        !
        !If (Me%PropCalc%Silica) then
        !   
        !    totalSi = Me%ExternalVar%Mass(Me%PropIndex%BiogenicSilica , index)                    &
        !            + Me%ExternalVar%Mass(Me%PropIndex%DissolvedSilica, index)                    &
        !            + Me%ExternalVar%Mass(Me%PropIndex%Diatoms, index) * Me%Diatoms%DiaAlfaSiC                 
        !
        !    write (*,*) 's/des...Depois..TotalSi                    ',totalSi, index
        !
        !endif
!
        !----------------------------------------------------------------------

                                    call WQOxygen       (index)
        if (Me%PropCalc%BOD       ) call WQBOD          (index)
        if (Me%PropCalc%Phosphorus) call WQPhosphorus   (index)
        if (Me%PropCalc%Nitrogen  ) call WQNitrogen     (index)
        if (Me%PropCalc%Phyto     ) call WQPhytoplankton(index)
        if (Me%PropCalc%Zoo       ) call WQZooplankton  (index)
        if (Me%PropCalc%Bacteria  ) call WQBacteria     (index)
        if (Me%PropCalc%Ciliate   ) call WQCiliate      (index)
        if (Me%PropCalc%Silica    ) call WQSilica       (index)  !aqui
        if (Me%PropCalc%Diatoms   ) call WQDiatoms      (index)  !aqui
        if (Me%PropCalc%Age       ) call WQAge          (index)
        if (Me%PropCalc%Larvae    ) call WQLarvae       (index)  ! Aires
        if (Me%PropCalc%Pompools  ) call WQPompools     (index)

!


        !----------------------------------------------------------------------

    end subroutine WQCoeficientsCalculation

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine WQSystemResolution(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)             :: index

        !Extertnal-------------------------------------------------------------
        integer                         :: STAT_CALL
        real, pointer, dimension(:)     ::  x

        !Local-----------------------------------------------------------------
        integer                         :: PropLB, PropUB
        integer                         :: prop
        integer                         :: equa

        !----------------------------------------------------------------------

        propLB = Me%Prop%ILB 
        propUB = Me%Prop%IUB 

        !Resolution using an explicit method
cd1 :   if (Me%CalcMethod%ExplicitMethod) then
do1 :       do equa = PropLB, PropUB           !Percorre as equacoes
do2 :       do prop = PropLB, PropUB           !Percorre as propriedades
cd2 :       if (Me%Matrix(equa, prop) .NE. 0.)  then
cd3 :           if (equa .EQ. prop) then
                    Me%IndTerm(equa) =  Me%IndTerm         (equa)                      &
                                     -((Me%Matrix          (equa, prop) - 1.0)         &
                                      * Me%ExternalVar%Mass(prop, index))
                else cd3
                    Me%IndTerm(equa) =  Me%IndTerm         (equa      )                &
                                     -  Me%Matrix          (equa, prop)                &
                                     *  Me%ExternalVar%Mass(prop, index)
                end if cd3

                Me%NewMass(equa)     = Me%IndTerm(equa)
            end if cd2
            end do do2
            end do do1


do4 :       do equa = PropLB, PropUB           !Percorre as equacoes
                Me%ExternalVar%Mass(equa, index) = Me%NewMass(equa)
            end do do4

        else if (Me%CalcMethod%ImplicitMethod) then

            !Resolution using an implicit method
            nullify   (x)

            call LUD(Me%ObjLUD,                                               &
                     Me%Matrix,                                               &
                     Me%IndTerm,                                              &
                     x,                                                       &
                     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                      &
                stop 'Subroutine WQSystemResolution - ModuleWaterQuality. ERR03.'




do3 :       do prop = PropLB, PropUB
                Me%ExternalVar%Mass(prop, index) = x(prop)
            end do do3

            nullify   (x)
        
    else if (Me%CalcMethod%SemiimpMethod) then

do31 :       do equa = PropLB, PropUB           !Percorre as equacoes
do32 :       do prop = PropLB, PropUB           !Percorre as propriedades
cd32 :           if (Me%Matrix(equa, prop) .GT. 0.)  then
cd33 :               if (equa .EQ. prop) then
                        Me%IndTerm(equa) =  Me%IndTerm         (equa      )        &
                                         -((Me%Matrix          (equa, prop) - 1.0) &
                                          * Me%ExternalVar%Mass(prop, index))
                                 
                        Me%Matrix(equa, prop) = 1.0
                    else
                        Me%IndTerm(equa) =  Me%IndTerm         (equa      )        &
                                         -  Me%Matrix          (equa, prop)        &
                                         *  Me%ExternalVar%Mass(prop, index)

                        Me%Matrix(equa, prop) = 0.0
                    end if cd33
                end if cd32
            end do do32
            end do do31
       


            !Resolution using an implicit method
            nullify   (x)
            
            call LUD(Me%ObjLUD,                                               &
                     Me%Matrix,                                               &
                     Me%IndTerm,                                              &
                     x,                                                       &
                     STAT = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                      &
                stop 'Subroutine WQSystemResolution - ModuleWaterQuality. ERR04.'


do33 :      do prop = PropLB, PropUB
                Me%ExternalVar%Mass(prop, index) = x(prop)
            end do do33
            
       
            nullify   (x)
        end if cd1

        !----------------------------------------------------------------------

    end subroutine WQSystemResolution

    !--------------------------------------------------------------------------

    subroutine WQRatesCalculation(index)
    
    
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                 :: index

        !Local-----------------------------------------------------------------
        integer                             :: Phyto
        integer                             :: Diatoms    !aqui
        real                                :: DTSec
        integer                             :: PropUB,PropLB,equa,prop
        type(T_PropRateFlux),       pointer :: PropRateFluxX
        type(T_EquaRateFlux),       pointer :: EquaRateFluxX 

        !----------------------------------------------------------------------

        Phyto   = Me%PropIndex%Phyto
        Diatoms = Me%PropIndex%Diatoms          !aqui

        propLB  = Me%Prop%ILB 
        propUB  = Me%Prop%IUB
        DTSec   = Me%DTSecond

        
cd0:    if (Me%PropCalc%Phyto) then
           
           Me%GrossProduction%field(index)= Me%GrowMaxPhytoRate                                  &
                                          * Me%TPhytoLimitationFactor                            &
                                          * Me%PhytoLightLimitationFactor                        &
                                          * Me%PhytoNutrientsLimitationFactor                    &
                                          * Me%ExternalVar%Mass(Phyto,index)                     &
                                          * DTSec
                                       
                         
           Me%NutLimitation%field(index)  = Me%PhytoNutrientsLimitationFactor * DTSec
           Me%NLimitation%field(index)    = Me%PhytoNLimitationFactor         * DTSec
           Me%PLimitation%field(index)    = Me%PhytoPLimitationFactor         * DTSec
           Me%LightLimitation%field(index)= Me%PhytoLightLimitationFactor     * DTSec
           Me%TempLimitation%field(index) = Me%TPhytoLimitationFactor         * DTSec
 
        endif    cd0                    
             
        !aqui_5
        if (Me%PropCalc%Diatoms) then
           
           Me%Diatoms%DiaGrossProduction%field(index)= Me%Diatoms%DiaGrowMaxRate                 &
                                                     * Me%Diatoms%DiaTempLimitationFactor        &
                                                     * Me%Diatoms%DiaLightLimitationFactor       &
                                                     * Me%Diatoms%DiaNutrientsLimitationFactor   &
                                                     * Me%ExternalVar%Mass(Diatoms,index)        &
                                                     * DTSec
                             
           Me%Diatoms%DiaNutLimitation%field(index)  = Me%Diatoms%DiaNutrientsLimitationFactor * DTSec
           Me%Diatoms%DiaNLimitation%field(index)    = Me%Diatoms%DiaNLimitationFactor         * DTSec
           Me%Diatoms%DiaPLimitation%field(index)    = Me%Diatoms%DiaPLimitationFactor         * DTSec
           Me%Diatoms%DiaSiLimitation%field(index)   = Me%Diatoms%DiaSiLimitationFactor        * DTSec
           Me%Diatoms%DiaLightLimitation%field(index)= Me%Diatoms%DiaLightLimitationFactor * DTSec
           Me%Diatoms%DiaTempLimitation%field(index) = Me%Diatoms%DiaTempLimitationFactor  * DTSec
 
        endif
        !aqui_5
        
           
          
        EquaRateFluxX => Me%FirstEquaRateFlux

do1:    do while(associated(EquaRateFluxX))  
                    
            PropRateFluxX => EquaRateFluxX%FirstPropRateFlux
                    
            do while(associated(PropRateFluxX))  
                         
                equa = EquaRateFluxX%ID
                prop = PropRateFluxX%ID

                if (equa.ne.prop) then

                     PropRateFluxX%Field(index)=  -Me%Matrix(equa, prop)                     & 
                                                  * Me%ExternalVar%Mass(prop,index)
                else
                      !the rates were inconsistent with system resolution and results
!                     PropRateFluxX%Field(index)=  -(1-Me%Matrix(equa, prop))                 & 
!                                                  * Me%ExternalVar%Mass(prop,index)
                     PropRateFluxX%Field(index)=  -(Me%Matrix(equa, prop) - 1.)              & 
                                                  * Me%ExternalVar%Mass(prop,index)

                endif

                PropRateFluxX => PropRateFluxX%Next

            enddo

            EquaRateFluxX => EquaRateFluxX%Next

        enddo do1

        !unidades

        !RateFlux = Matrix * Mass
        !mg/l        1      mg/l 

        !Matrix = dt * K
        !   1   = DtDay * 1/Days
        
        nullify(PropRateFluxX)
        nullify(EquaRateFluxX)

    end subroutine WQRatesCalculation

     !----------------------------------------------------------------------------

    subroutine WQBacteria(index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN) :: index

        !Local-----------------------------------------------------------------
        integer                             :: DONnr
        integer                             :: AM
        integer                             :: BAC
        integer                             :: PON
        integer                             :: PONr
        integer                             :: O


        real                                :: s1, s2, xa, xb, ya, yb 
        real                                :: DTDay

        real                                :: TBacteriaLimitationFactor    = null_real
        real                                :: BacteriaDONUptake            = null_real
        real                                :: BacteriaAmmoniaUptake        = null_real
        real                                :: BacteriaPONUptake            = null_real
        real                                :: BacteriaTotalUptake          = null_real
        real                                :: BacteriaExcretionNitrogen    = null_real
        real                                :: DeadBacteria                 = null_real

        !----------------------------------------------------------------------


        ! Bacteria Common Coeficients ------------
        DONnr   = Me%PropIndex%DONNonRefractory
        PON     = Me%PropIndex%PartOrganicNitrogen
        PONr    = Me%PropIndex%PartOrganicNitrogenRefractory
        AM      = Me%PropIndex%Ammonia
        BAC     = Me%PropIndex%Bacteria
        O       = Me%PropIndex%Oxygen

        DTDay   = Me%DTDay

   


         !FoodBacteriaLimitationFactor, food availability limitation factor
!cd1 :       if ((Me%ExternalVar%Mass(numAmmonia, index)   &
!               + Me%ExternalVar%Mass(numDONnr, index)     &
!               - Me%NutMinBact) .LE. 0.0) then
!                 FoodBactLimitationFactor = 0.0
!            else
!                 BactExponent = -Me%IvlevAssimBactConst                   &
!                             * (Me%ExternalVar%Mass(numAmmonia, index)    &
!                             + Me%ExternalVar%Mass(numDONnr, index)       & 
!                             - Me%NutMinBact)

!                 FoodBactLimitationFactor = 1.0 - exp(BactExponent)
!            end if cd1




       !TBacteriaLimitationFactor : limitation by temperature

        s1 = (1. / (Me%TOptBacteriaMin - Me%TBacteriaMin)) &
        * log((Me%BK2 * (1.0 - Me%BK1))                 &
           / (Me%BK1 * (1.0 - Me%BK2)))

        s2 = (1. / (Me%TBacteriaMax - Me%TOptBacteriaMax)) &
        * log((Me%BK3 * (1.0 - Me%BK4))                 &
           / (Me%BK4 * (1.0 - Me%BK3)))

        ya = exp(s1 * (Me%ExternalVar%Temperature(index) - Me%TBacteriaMin))
        yb = exp(s2 * (Me%TBacteriaMax - Me%ExternalVar%Temperature(index)))

        xa = (Me%BK1 * ya) / (1.0 + Me%BK1 * (ya - 1.0))
        xb = (Me%BK4 * yb) / (1.0 + Me%BK4 * (yb - 1.0))

        TBacteriaLimitationFactor = xa * xb

        !Needs C and AM to grow. C is availiable through PON and DONnr
        if (((Me%ExternalVar%Mass(PON,index  ) .GE. Me%BacteriaMinSubstrate)  .OR.       &
            ( Me%ExternalVar%Mass(DONnr,index) .GE. Me%BacteriaMinSubstrate)) .AND.      &
            ( Me%ExternalVar%Mass(AM,index   ) .GE. Me%BacteriaMinSubstrate)) then

            !BacteriaAmmoniaUptake (mgN/L)
            BacteriaAmmoniaUptake = TBacteriaLimitationFactor                            &
                           * Me%BacteriaMaxUptake                                        &
                           * (Me%ExternalVar%Mass(AM, index)                             &
                           /(Me%NitrogenSaturationConstBacteria                          &
                           +  Me%ExternalVar%Mass(AM, index)))    

            !PON exists
            if (Me%ExternalVar%Mass(PON,index  ) .GE. Me%BacteriaMinSubstrate) then

                !BacteriaPONUptake (mgN/L)
                BacteriaPONUptake  = TBacteriaLimitationFactor                              &
                           * Me%BacteriaMaxUptake                                           &
                           * (Me%ExternalVar%Mass(PON, index)                               &
                           / (Me%NitrogenSaturationConstBacteria                            &
                           + Me%ExternalVar%Mass(PON, index))) 

            else

                BacteriaPONUptake  = 0.0

            endif

            !DONnr exists
            if (Me%ExternalVar%Mass(DONnr,index) .GE. Me%BacteriaMinSubstrate) then

                !BacteriaDONUptake, (mgN/L)
                BacteriaDONUptake =  TBacteriaLimitationFactor                              &
                          * Me%BacteriaMaxUptake                                            &
                          * (Me%ExternalVar%Mass(DONnr, index)                              &
                          / (Me%NitrogenSaturationConstBacteria                             &
                          + Me%ExternalVar%Mass(DONnr, index)))

            else

                BacteriaDONUptake = 0.0

            endif

        else

            BacteriaAmmoniaUptake = 0.0
            BacteriaPONUptake     = 0.0
            BacteriaDONUptake     = 0.0

        endif
        
        !BacteriaTotalUptake, uptake in mgC/L
        BacteriaTotalUptake = (BacteriaAmmoniaUptake                                    &
                     / Me%AlfaBacteriaNC                                                &
                     + (BacteriaDONUptake                                               &
                     + BacteriaPONUptake)                                               &
                     / Me%OMAlfaNC )                       
                 


    !BacteriaExcretionNitrogen, excretion of soluble Nitrogen compounds by the bacteria  (mgN/L)

     BacteriaExcretionNitrogen = Me%AlfaBacteriaNC * Me%BacteriaExcretionRate
       
    
            
            
    !DeadBacteria Organic nitrogen from dead bacteria  (mgN/L)
            DeadBacteria = Me%BacteriaNonGrazingMortalityRate * Me%AlfaBacteriaNC 
       
      
                                    
    !Calculation of system coeficients-----------------------------------------
        Me%Matrix(BAC, BAC)   = 1.0 + DTDay * (Me%BacteriaExcretionRate                    &
                                            + Me%BacteriaNonGrazingMortalityRate    &
                                            - BacteriaTotalUptake )

        Me%Matrix(PON, BAC)   = DTDay * BacteriaPONUptake

        Me%Matrix(PONr, BAC)  = - DTDay * DeadBacteria


        Me%Matrix(AM, BAC)    = (BacteriaAmmoniaUptake - BacteriaExcretionNitrogen) &
                              * DTDay     


        Me%Matrix(DONnr, Bac) = + DTDay * BacteriaDONUptake
        
        
        !Oxygen---------------------------------------
   
        if (Me%PropCalc%Oxygen) then                                                                     
            
            Me%Matrix(O,  Bac)     =   DTDay * ((BacteriaPONUptake + BacteriaDONUptake) / Me%OMAlfaNC)    &
                                             *  Me%BactAlfaOC      
   
        endif        
        


    !Independent term
       Me%IndTerm(BAC)        = Me%ExternalVar%Mass(BAC, index) 

    !----------------------------------------------------------------------
 
end subroutine WQBacteria






subroutine WQCiliate(index)

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-----------------------------------------------------------------

        integer :: BAC
        integer :: CIL
        integer :: O
        integer :: PON
        integer :: AM
        integer :: DONnr
        integer :: Phyto
        integer :: IP
        integer :: POP
        integer :: DOPnr

        
        real :: OxygenCiliateRespRate               = null_real
        real :: CiliateGrazBactLimitationFactor     = null_real
        real :: CiliateGrazPhyLimitationFactor      = null_real
        real :: CiliateIngestionBacteria            = null_real
        real :: CiliateIngestionPhyto               = null_real
        real :: CiliateNaturalMortality             = null_real
        real :: DeadCiliateNitrogen                 = null_real  
        real :: DeadCiliatePhosphorus               = null_real
        real :: LostChainNitrogen                   = null_real
        real :: LostChainPhosphorus                 = null_real
        reaL :: LostGrazPhosphorus                  = null_real
        reaL :: LostGrazNitrogen                    = null_real

  
        real :: CiliateRespirationRate              = null_real
        real :: CiliateGrossGrowRate                = null_real
        real :: CiliateExcretionRate                = null_real  
        real :: CiliateExcretionNitrogen            = null_real 
        real :: CiliateExcretionPhosphorus          = null_real 
        real :: CiliateGrazPrey                     = null_real
        real :: x, y, xa, ya
        real :: DTDay
        

    !----------------------------------------------------------------------


        BAC   = Me%PropIndex%Bacteria
        CIL   = Me%PropIndex%Ciliate
        O     = Me%PropIndex%Oxygen
        PON   = Me%PropIndex%PartOrganicNitrogen
        AM    = Me%PropIndex%Ammonia
        DONnr = Me%PropIndex%DONNonRefractory
        Phyto = Me%PropIndex%Phyto
        POP   = Me%PropIndex%PartOrganicPhosphorus
        DOPnr = Me%PropIndex%DOPNonRefractory
        IP    = Me%PropIndex%InorganicPhosphorus


        x       = null_real
        y       = null_real
        xa      = null_real
        ya      = null_real

        
        DTDay = Me%DTDay


    !Ciliates Grossgrowrate/NaturalMortality/ZooGrazPrey/LostChain/LostGraz-----------------
    select case (Me%WQConfiguration)

        case (ZooPhyDiaCilBac,ZooPhyCilBac)

            !Ciliates Grossgrowrate---------------------------------------------------------

                !limitation by phyto (1º) 
                x = Me%CiliateEfficiencyCapturePhyto * Me%ExternalVar%Mass(Phyto, index)          &
                  - Me%CiliateGrazPhytoMin                                                             
                                                                                    
                y = Me%CiliateIngestionConst + (Me%CiliateEfficiencyCapturePhyto                    &
                  * Me%ExternalVar%Mass(Phyto, index) - Me%CiliateGrazPhytoMin)                         

                if  (x .LE. 0.0) then                                                                   
                     CiliateGrazPhyLimitationFactor = 0.                                                
                else
                     CiliateGrazPhyLimitationFactor = x/ y
                end if                                                                                  

                !CiliatesIngestion of phyto
                CiliateIngestionPhyto = CiliateGrazPhyLimitationFactor                              &
                                        * Me%PhytoRatioIngestionCiliates                            &
                                        * Me%CiliateIngestionMax                                    &
                                        * Me%TZooLimitationFactor                                         
                                                    
 
                !limitation by Bacteria (2º)                                                                 
                xa = (Me%CiliateEfficiencyCaptureBact * Me%ExternalVar%Mass(BAC, index))            &
                  - Me%CiliateGrazBactMin                                                               
                                                                                        
                ya = Me%CiliateIngestionConst + (Me%CiliateEfficiencyCaptureBact                    &
                  * Me%ExternalVar%Mass(BAC, index) - Me%CiliateGrazBactMin)                            
                                                                                        
                if  (xa .LE. 0.0) then                                                               
                     CiliateGrazBactLimitationFactor = 0.                                               
                else                                                                                    
                     CiliateGrazBactLimitationFactor = xa/ ya                                             
                end if                                                                               
                                                                                        
                !CiliatesIngestion of Bacteria                
                CiliateIngestionBacteria = CiliateGrazBactLimitationFactor                          &
                                           * Me%BactRatioIngestionCiliates                          &
                                           * (Me%CiliateIngestionMax-CiliateIngestionPhyto)         &
                                           * Me%TZooLimitationFactor 
            
                !Ciliate gross growth rate
                CiliateGrossGrowRate = Me%CiliateAssimilationPhytoRate * CiliateIngestionPhyto      & 
                                     + Me%CiliateAssimilationBacteriaRate * CiliateIngestionBacteria
                         
                if (CiliateGrossGrowRate .LT. 0.0) then                                                     
                    stop 'Subroutine WQCiliate-Module ModuleWaterQuality. ERR10.'
                !No sense on this. David
!                elseif (CiliateGrossGrowRate .EQ. 0.0) then
!                    CiliateGrossGrowRate =-1.0 / null_real   !Avoid division by zero below
                end if 

            !CiliateGrazPrey for NaturalMortality Calculation--------------------------------------
                            
                CiliateGrazPrey  = Me%ExternalVar%Mass(Phyto, index)                                &
                                  + Me%ExternalVar%Mass(BAC, index)                                    
            
            !LostChain, term due to different Phyto, Bacteria and Ciliates ratios----------- 
            !and LostGraz, term due to not 100% assimilation----------------- 
                            
                !Nitrogen
                if (Me%PropCalc%Nitrogen) then

                    LostChainNitrogen = (Me%AlfaPhytoNC-Me%AlfaCilNC)                               &
                                      * (Me%CiliateAssimilationPhytoRate * CiliateIngestionPhyto)   &
                                      +(Me%AlfaBacteriaNC-Me%AlfaCilNC)                             &
                                      * (Me%CiliateAssimilationBacteriaRate * CiliateIngestionBacteria)    
                
                    LostGrazNitrogen = (1- Me%CiliateAssimilationPhytoRate)                         &
                                        * CiliateIngestionPhyto * Me%AlfaPhytoNC                    &
                                      +(1- Me%CiliateAssimilationBacteriaRate)                      &
                                        * CiliateIngestionBacteria * Me%AlfaBacteriaNC   

                    if (LostChainNitrogen .LT. 0.0)                                                 &
                        stop 'Subroutine WQCiliate-Module ModuleWaterQuality. ERR20.'
             
                    if (LostGrazNitrogen .LT. 0.0)                                                  &
                        stop 'Subroutine WQCiliate-Module ModuleWaterQuality. ERR30.'

                endif
                
                
                !Phosphorus
                if (Me%PropCalc%Phosphorus) then

                    LostChainPhosphorus = (Me%AlfaPhytoPC-Me%AlfaCilPC)                             &
                                      * (Me%CiliateAssimilationPhytoRate * CiliateIngestionPhyto)        
                
                    LostGrazPhosphorus = (1- Me%CiliateAssimilationPhytoRate)                       &
                                        * Me%CiliateAssimilationPhytoRate * Me%AlfaPhytoPC                        

                    if (LostChainPhosphorus .LT. 0.0)                                               &
                        stop 'Subroutine WQCiliate-Module ModuleWaterQuality. ERR40.'
             
                    if (LostGrazPhosphorus .LT. 0.0)                                                &
                        stop 'Subroutine WQCiliate-Module ModuleWaterQuality. ERR50.'

                endif
 
        !endcase (ZooPhyDiaBact,ZooPhyCilBac)
        case (ZooPhyDiaCil,ZooPhyCil)

            !Ciliates Grossgrowrate---------------------------------------------------------

                !limitation by phyto 
                x = (Me%CiliateEfficiencyCapturePhyto * Me%ExternalVar%Mass(Phyto, index)          &
                  - Me%CiliateGrazPhytoMin)                                                             
                                                                                    
                y = Me%CiliateIngestionConst + (Me%CiliateEfficiencyCapturePhyto                    &
                  * Me%ExternalVar%Mass(Phyto, index) - Me%CiliateGrazPhytoMin)                         

                if  (x .LE. 0.0) then                                                                   
                     CiliateGrazPhyLimitationFactor = 0.                                                
                else
                     CiliateGrazPhyLimitationFactor = x/ y
                end if                                                                                  

                !CiliatesIngestion of phyto
                CiliateIngestionPhyto = CiliateGrazPhyLimitationFactor                              &
                                        * Me%PhytoRatioIngestionCiliates                            &
                                        * Me%CiliateIngestionMax                                    &
                                        * Me%TZooLimitationFactor                                         
                                                               
                !Ciliate gross growth rate
                CiliateGrossGrowRate = Me%CiliateAssimilationPhytoRate * CiliateIngestionPhyto       
                         
                if (CiliateGrossGrowRate .LT. 0.0) then                                                     
                    stop 'Subroutine WQCiliate-Module ModuleWaterQuality. ERR60.'
                !No sense on this. David
!                elseif (CiliateGrossGrowRate .EQ. 0.0) then
!                    CiliateGrossGrowRate =-1.0 / null_real   !Avoid division by zero below
                end if 

            !CiliateGrazPrey for NaturalMortality Calculation--------------------------------------
                            
                CiliateGrazPrey  = Me%ExternalVar%Mass(Phyto, index)                                
            

            !LostChain, term due to different Phyto and Ciliates ratios----------- 
            !and LostGraz, term due to not 100% assimilation----------------- 
                            
                !Nitrogen
                if (Me%PropCalc%Nitrogen) then

                    LostChainNitrogen = (Me%AlfaPhytoNC-Me%AlfaCilNC)                               &
                                      * (Me%CiliateAssimilationPhytoRate * CiliateIngestionPhyto)      

                    LostGrazNitrogen = (1- Me%CiliateAssimilationPhytoRate)                         &
                                        * CiliateIngestionPhyto * Me%AlfaPhytoNC                    

                    if (LostChainNitrogen .LT. 0.0)                                                 &
                        stop 'Subroutine WQCiliate-Module ModuleWaterQuality. ERR70.'
             
                    if (LostGrazNitrogen .LT. 0.0)                                                  &
                        stop 'Subroutine WQCiliate-Module ModuleWaterQuality. ERR80.'

                endif
                
                
                !Phosphorus
                if (Me%PropCalc%Phosphorus) then

                    LostChainPhosphorus = (Me%AlfaPhytoPC-Me%AlfaCilPC)                             &
                                      * (Me%CiliateAssimilationPhytoRate * CiliateIngestionPhyto)        
                
                    LostGrazPhosphorus = (1- Me%CiliateAssimilationPhytoRate)                       &
                                        * Me%CiliateAssimilationPhytoRate * Me%AlfaPhytoPC                        

                    if (LostChainPhosphorus .LT. 0.0)                                               &
                        stop 'Subroutine WQCiliate-Module ModuleWaterQuality. ERR90.'
             
                    if (LostGrazPhosphorus .LT. 0.0)                                                &
                        stop 'Subroutine WQCiliate-Module ModuleWaterQuality. ERR100.'

                endif

        !end case (ZooPhyDiaCil,ZooPhyCil)
        case (ZooDiaCilBac)

            !Ciliates Grossgrowrate---------------------------------------------------------
                                                    
                !limitation by Bacteria (1º)                                                                 
                xa = (Me%CiliateEfficiencyCaptureBact * Me%ExternalVar%Mass(BAC, index))            &
                  - Me%CiliateGrazBactMin                                                               
                                                                                        
                ya = Me%CiliateIngestionConst + (Me%CiliateEfficiencyCaptureBact                    &
                  * Me%ExternalVar%Mass(BAC, index) - Me%CiliateGrazBactMin)                            
                                                                                        
                if  (xa .LE. 0.0) then                                                               
                     CiliateGrazBactLimitationFactor = 0.                                               
                else                                                                                    
                     CiliateGrazBactLimitationFactor = xa/ ya                                             
                end if                                                                               
                                                                                        
                !CiliatesIngestion of Bacteria                
                CiliateIngestionBacteria = CiliateGrazBactLimitationFactor                          &
                                           * Me%BactRatioIngestionCiliates                          &
                                           * Me%CiliateIngestionMax                                 &
                                           * Me%TZooLimitationFactor 
            
                !Ciliate gross growth rate
                CiliateGrossGrowRate = Me%CiliateAssimilationBacteriaRate * CiliateIngestionBacteria
                         
                if (CiliateGrossGrowRate .LT. 0.0) then                                                     
                    stop 'Subroutine WQCiliate-Module ModuleWaterQuality. ERR110.'
                !No sense on this. David
!                elseif (CiliateGrossGrowRate .EQ. 0.0) then
!                    CiliateGrossGrowRate =-1.0 / null_real   !Avoid division by zero below
                end if 

            !CiliateGrazPrey for NaturalMortality Calculation--------------------------------------
                            
                CiliateGrazPrey  = Me%ExternalVar%Mass(BAC, index)                                    
            

            !LostChain, term due to different  Bacteria and Ciliates ratios----------- 
            !and LostGraz, term due to not 100% assimilation----------------- 
                            
                !Nitrogen
                if (Me%PropCalc%Nitrogen) then

                    LostChainNitrogen = (Me%AlfaBacteriaNC-Me%AlfaCilNC)                             &
                                      * (Me%CiliateAssimilationBacteriaRate * CiliateIngestionBacteria)    
                
                    LostGrazNitrogen = (1- Me%CiliateAssimilationBacteriaRate)                      &
                                        * CiliateIngestionBacteria * Me%AlfaBacteriaNC   

                    if (LostChainNitrogen .LT. 0.0)                                                 &
                        stop 'Subroutine WQCiliate-Module ModuleWaterQuality. ERR120.'
             
                    if (LostGrazNitrogen .LT. 0.0)                                                  &
                        stop 'Subroutine WQCiliate-Module ModuleWaterQuality. ERR130.'

                endif
                
                
                !Phosphorus
                if (Me%PropCalc%Phosphorus) then

                    !Bacteria dont have phosphurus
                    LostChainPhosphorus = 0.0
                    
                    LostGrazPhosphorus  = 0.0  
                endif

        !endcase (ZooDiaCilBac)
        end select

        !CiliatesNaturalMortality----------------------------------------------------------

            if (CiliateGrazPrey .LE. Me%CiliateGrazPreyMin ) then
                CiliateNaturalMortality = Me%CiliateMaxMortalityRate
            else                                                                                    
                CiliateNaturalMortality = (Me%CiliateMortalityCoef / CiliateGrazPrey)               & 
                                          + Me%CiliateMinMortalityRate
            end if  
            
            
            if (Me%PropCalc%Nitrogen) then

                !DeadCiliateNitrogen, Organic Nitrogen from dead Ciliate (mgN/L)
                DeadCiliateNitrogen   = CiliateNaturalMortality * Me%AlfaCilNC
            
            endif

            if (Me%PropCalc%Phosphorus) then

                !DeadZooPhosphorus, Organic Phosphorus from dead Ciliate (mgP/L)
                DeadCiliatePhosphorus = CiliateNaturalMortality * Me%AlfaCilPC 
            
            endif

        !Ciliates RespirationRate -------------------------------------------------------
        
            CiliateRespirationRate = Me%CiliateReferenceRespirationRate * Me%TZooLimitationFactor
       
            !Oxygen loss due to Ciliate respiration
            if (Me%PropCalc%Oxygen)  then
                
                OxygenCiliateRespRate = Me%RatioOxygenCarbonCilRespiration * CiliateRespirationRate
            
                if (MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen).EQ.Me%MinOxygen) then
                    CiliateGrossGrowRate  = -1.0 / null_real
                    OxygenCiliateRespRate = -1.0 / null_real
                endif
            
            endif 

        !Ciliate ExcretionRate -------------------------------------------------------

            CiliateExcretionRate = Me%CiliateExcretionFactor  * Me%CiliateExcretionConst            &
                               ** Me%ExternalVar%Temperature(index)
                
            if (Me%PropCalc%Nitrogen) then

                !CiliateExcretionNitrogen (mg N/L)
                CiliateExcretionNitrogen = CiliateExcretionRate * Me%AlfaCilNC
            
            endif

            if (Me%PropCalc%Phosphorus) then

                !ZooExcretionPhosphorus (mg P/L)
                CiliateExcretionPhosphorus = CiliateExcretionRate * Me%AlfaCilPC 
            
            endif
    
    !Calculation of system coeficients-----------------------------------------             
    
    !Organisms---------------------------------------

        Me%Matrix(CIL, CIL)  = 1.0 + DTDay * (CiliateExcretionRate                                  &
                                                 + CiliateRespirationRate                           &
                                                 + CiliateNaturalMortality                          &
                                                 - CiliateGrossGrowRate)

        if (Me%PropCalc%Bacteria) then

            Me%Matrix(BAC, CIL)=  DTDay * CiliateIngestionBacteria
        
        endif

        if (Me%PropCalc%Phyto) then

            Me%Matrix(Phyto, CIL)=  DTDay * CiliateIngestionPhyto
        
        endif

    !Oxygen---------------------------------------
   
        if (Me%PropCalc%Oxygen) then                                                                     
            Me%Matrix(O,  CIL)     =   DTDay * OxygenCiliateRespRate   
        endif

    !Nitrogen---------------------------------------
        if (Me%PropCalc%Nitrogen) then
        
            Me%Matrix(AM, CIL   ) = - DTDay * (CiliateExcretionNitrogen * Me%ZooSolublInorgExcreFraction &
                                              +CiliateRespirationRate   * Me%AlfaCilNC)    

            Me%Matrix(DONnr, CIL) = - DTDay * (CiliateExcretionNitrogen                             & 
                                               * Me%ZooExcreDissOrgFraction                         & 
                                               * (1.0 - Me%ZooSolublInorgExcreFraction))           
            
            Me%Matrix(PON, CIL  ) = - DTDay * (CiliateExcretionNitrogen                             & 
                                               * (1.0 - Me%ZooExcreDissOrgFraction)                 & 
                                               * (1.0 - Me%ZooSolublInorgExcreFraction)             &
                                              + DeadCiliateNitrogen                                 &
                                              + LostChainNitrogen                                   &
                                              + LostGrazNitrogen)
        endif                                    
 
   !Phosphorus---------------------------------------
        if (Me%PropCalc%Phosphorus) then
        
            Me%Matrix(IP, CIL   ) = -DTDay * (CiliateExcretionPhosphorus * Me%ZooSolublInorgExcreFraction &
                                            +CiliateRespirationRate * Me%AlfaCilPC)    

            
            Me%Matrix(DOPnr, CIL) = -DTDay * (CiliateExcretionPhosphorus                            &
                                              * Me%ZooExcreDissOrgFraction                          & 
                                              * (1.0 - Me%ZooSolublInorgExcreFraction))           

            Me%Matrix(POP, CIL  ) = -DTDay * (CiliateExcretionPhosphorus                            &
                                               * (1.0 - Me%ZooExcreDissOrgFraction)                 & 
                                               * (1.0 - Me%ZooSolublInorgExcreFraction)             &
                                              + DeadCiliatePhosphorus                               &
                                              + LostChainPhosphorus                                 &
                                              + LostGrazPhosphorus)
        endif

    !Independent term
        Me%IndTerm(CIL) = Me%ExternalVar%Mass(CIL, index) 
  
end subroutine WQCiliate

subroutine WQOxygen(index)

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !------------------------------------------------------------------------

cd1 :   if (Me%PropCalc%Oxygen) then
            call WQOxygenCalculation(index)
        else cd1
            call WQOxygenSaturation (index)
        end if cd1

    !------------------------------------------------------------------------

end subroutine WQOxygen

    !----------------------------------------------------------------------------





    !----------------------------------------------------------------------------

subroutine WQOxygenSaturation(index)

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !------------------------------------------------------------------------

        Me%ExternalVar%Mass(Me%PropIndex%Oxygen, index) =                           &
                                OxygenSaturation(Me%ExternalVar%Temperature(index), &
                                                 Me%ExternalVar%Salinity (index))


        !Calculation of system coeficients 
        Me%Matrix(Me%PropIndex%Oxygen, Me%PropIndex%Oxygen) = 1.0 


        !Independent term
        Me%IndTerm(Me%PropIndex%Oxygen) =                                           &
                                Me%ExternalVar%Mass(Me%PropIndex%Oxygen, index)

    !------------------------------------------------------------------------

end subroutine WQOxygenSaturation

    !----------------------------------------------------------------------------


  



    !----------------------------------------------------------------------------
    !OXYGEN
    !
    !SOURCES: - Phytoplankton Oxygen production due to photosynthesis;
    !         - Oxygen production related with the Phytoplankton Nitrate consumption.
    !
    !SINKS:   - Phytoplankton and Zooplankton respiration;
    !         - BOD oxidation (Carencia bioquimica de oxigenio).

subroutine WQOxygenCalculation(index)

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------

        integer :: O

        
    !------------------------------------------------------------------------

        O   = Me%PropIndex%Oxygen



    !Calculation of system coeficients---------------------------------------
        Me%Matrix(O, O) = 1.0 


    !Independent term
        Me%IndTerm(O) = Me%ExternalVar%Mass(O, index)

    !------------------------------------------------------------------------

end subroutine WQOxygenCalculation

    !----------------------------------------------------------------------------


  




    !----------------------------------------------------------------------------

    subroutine WQBOD(index)

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------

        integer :: BOD
        integer :: O

        real    :: DTDay
        real    :: BODOxidationRate    = null_real 

    !------------------------------------------------------------------------

        BOD     = Me%PropIndex%BOD 
        O       = Me%PropIndex%Oxygen 

        DTDay   = Me%DTDay




    !BODOxidationRate, BOD oxidation rate, 1/T
cd7 :   if (Me%PropCalc%Oxygen) then
            
            BODOxidationRate = Me%BODOxidationReferenceRate                         &
                             * Me%BODOxidationCoefficient                           &
                             ** (Me%ExternalVar%Temperature(index) - 20.0)          &
                             * (MAX(Me%ExternalVar%Mass(O, index), Me%MinOxygen)    &
                             /(Me%BODOxygenSSatConstant                             &
                             + MAX(Me%ExternalVar%Mass(O, index), Me%MinOxygen)))
        end if cd7




    !Calculation of system coeficients---------------------------------------
         Me%Matrix(BOD, BOD) = DTDay * BODOxidationRate + 1.0


        if (Me%PropCalc%Oxygen) Me%Matrix(O, BOD) = DTDay * BODOxidationRate            


    !Independent term
        Me%IndTerm(BOD) = Me%ExternalVar%Mass(BOD, index)

    !------------------------------------------------------------------------

    end subroutine WQBOD

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine WQZooplankton(index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN) :: index

        !Local-------------------------------------------------------------------

        integer :: IP
        integer :: POP
        integer :: DOPnr
        integer :: AM
        integer :: PON
        integer :: Zoo
        integer :: Phyto
        integer :: Diatoms
        integer :: BioSi
        integer :: DissSi
        integer :: DONnr
        integer :: O
        integer :: CIL

        real    :: DTDay

        real    :: ZooExcretionPhosphorus           = null_real
        real    :: ZooExcretionNitrogen             = null_real
        reaL    :: LostGrazPhosphorus               = null_real
        reaL    :: LostGrazNitrogen                 = null_real
        real    :: exponent                         = null_real
        real    :: FoodZooLimitationFactor          = null_real
        real    :: ZooRespirationRate               = null_real
        real    :: ZooGrossGrowRate                 = null_real
        real    :: LostChainNitrogen                = null_real
        real    :: LostChainPhosphorus              = null_real
        real    :: OxygenZooRespRate                = null_real
        real    :: ZooGrazPhytoLimitationFactor     = null_real
        real    :: ZooGrazDiatomsLimitationFactor   = null_real
        real    :: ZooGrazCiliateLimitationFactor   = null_real
        real    :: ZooGrazPrey                      = null_real !aqui
        real    :: ZooIngestionPhyto                = null_real 
        real    :: ZooIngestionCiliate              = null_real 
        real    :: ZooIngestionDiatoms              = null_real 
        real    :: ZooNaturalMortality              = null_real
        real    :: DeadZooPhosphorus                = null_real  
        real    :: DeadZooNitrogen                  = null_real  
        real    :: PredatedZooPhosphorus            = null_real
        real    :: PredatedZooNitrogen              = null_real
        real    :: ZooExcretionRate                 = null_real

        real    :: s1, s2, xa, xb, ya, yb, x, y, yac, ybc, xaox, xbox

    !------------------------------------------------------------------------

        POP     = Me%PropIndex%PartOrganicPhosphorus
        DOPnr   = Me%PropIndex%DOPNonRefractory
        IP      = Me%PropIndex%InorganicPhosphorus
        AM      = Me%PropIndex%Ammonia
        PON     = Me%PropIndex%PartOrganicNitrogen
        DONnr   = Me%PropIndex%DONNonRefractory
        Phyto   = Me%PropIndex%Phyto
        Diatoms = Me%PropIndex%Diatoms
        BioSi   = Me%PropIndex%BiogenicSilica   !aqui
        DissSi  = Me%PropIndex%DissolvedSilica  !aqui
        Zoo     = Me%PropIndex%Zoo
        O       = Me%PropIndex%Oxygen
        CIL     = Me%PropIndex%Ciliate
              
        DTDay   = Me%DTDay

        s1      = null_real
        s2      = null_real
        xa      = null_real
        xb      = null_real
        xaox    = null_real
        xbox    = null_real
        ya      = null_real
        yb      = null_real
        x       = null_real
        y       = null_real
        yac     = null_real
        ybc     = null_real

        
    !TZooLimitationFactor, temperature effect on zooplankton assimilation rate
        s1 = (1. / (Me%TOptZooMin - Me%TZooMin)) * log((Me%ZK2 * (1.0 - Me%ZK1))                    &
                                                     / (Me%ZK1 * (1.0 - Me%ZK2)))

        s2 = (1. / (Me%TZooMax - Me%TOptZooMax)) * log((Me%ZK3 * (1.0 - Me%ZK4))                    &
                                                     / (Me%ZK4 * (1.0 - Me%ZK3)))

        ya = exp(s1 * (Me%ExternalVar%Temperature(index) - Me%TZooMin))
        yb = exp(s2 * (Me%TZooMax - Me%ExternalVar%Temperature(index)))

        xa = (Me%ZK1 * ya) / (1.0 + Me%ZK1 * (ya - 1.0))
        xb = (Me%ZK4 * yb) / (1.0 + Me%ZK4 * (yb - 1.0))

        Me%TZooLimitationFactor = xa * xb

    !Zooplankton Grossgrowrate/NaturalMortality_--------------------------------------
    select case (Me%WQConfiguration)

        case (ZooPhyDiaCil, ZooPhyDiaCilBac)

            !Zooplankton Grossgrowrate------------------------------------------------

                !limitation by diatoms (1º) 
                yac = (Me%Diatoms%ZooEfficiencyCaptureDiatoms * Me%ExternalVar%Mass(Diatoms, index) &
                     - Me%Diatoms%GrazDiaMin)                                                               
                                                                                
                ybc = (Me%ZooIngestionConst + Me%Diatoms%ZooEfficiencyCaptureDiatoms                &
                    * Me%ExternalVar%Mass(Diatoms, index) - Me%Diatoms%GrazDiaMin)

                if  (yac .LE. 0.0) then
                    ZooGrazDiatomsLimitationFactor = 0.0
                else           
                    ZooGrazDiatomsLimitationFactor = yac / ybc
                end if 

                !ZooIngestion of diatoms
                ZooIngestionDiatoms = Me%Diatoms%DiaRatioIngestionZoo                               &
                                    * Me%ZooIngestionMax                                            &
                                    * ZooGrazDiatomsLimitationFactor                                &
                                    * Me%TZooLimitationFactor                       
                                                    
                !limitation by phyto (2º)
                x = Me%ZooEfficiencyCapturePhyto * Me%ExternalVar%Mass(Phyto, index)                & 
                - Me%GrazPhytoMin                                                                     
                                                                        
                y = Me%ZooIngestionConst + Me%ZooEfficiencyCapturePhyto                             &
                * Me%ExternalVar%Mass(Phyto, index) - Me%GrazPhytoMin                                 

                if  (x .LE. 0.0) then
                    ZooGrazPhytoLimitationFactor = 0.0
                else 
                    ZooGrazPhytoLimitationFactor = x / y
                end if

                !ZooIngestion of Phyto 
                ZooIngestionPhyto = Me%PhytoRatioIngestionZoo                                       &
                                  * (Me%ZooIngestionMax -  ZooIngestionDiatoms)                     &                         
                                  * ZooGrazPhytoLimitationFactor                                    &
                                  * Me%TZooLimitationFactor                                         

                !limitation by Ciliate 
                yac = (Me%ZooEfficiencyCaptureCiliate * Me%ExternalVar%Mass(CIL, index)             &
                     - Me%GrazCiliateMin)                                                           
                                                                                                     
                ybc = (Me%ZooIngestionConst + Me%ZooEfficiencyCaptureCiliate                        &
                    * Me%ExternalVar%Mass(CIL, index) - Me%GrazCiliateMin)

                if  (yac .LE. 0.0) then
                    ZooGrazCiliateLimitationFactor = 0.0
                else            
                    ZooGrazCiliateLimitationFactor = yac / ybc
                end if 

                !ZooIngestion of Ciliate
                ZooIngestionCiliate = Me%CiliatesRatioIngestionZoo                                  &
                                    * (Me%ZooIngestionMax - ZooIngestionDiatoms-ZooIngestionPhyto)  &
                                    * ZooGrazCiliateLimitationFactor                                &
                                    * Me%TZooLimitationFactor                                           
                                                                                
                !ZooGrossGrowRate, zooplankton gross growth rate                                        
                ZooGrossGrowRate =  (Me%Diatoms%DiaZooAssimilationRate * ZooIngestionDiatoms)       &
                                   +(Me%ZooAssimilationPhytoRate       * ZooIngestionPhyto  )       &
                                   +(Me%ZooAssimilationCiliateRate     * ZooIngestionCiliate)                 
                         
                if (ZooGrossGrowRate .LT. 0.0) then                                                     
                    stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR01.'
                !No sense on this. David
!                elseif (ZooGrossGrowRate .EQ. 0.0) then
!                    ZooGrossGrowRate =-1.0 / null_real   !Avoid division by zero below
                end if 
            
            !ZooGrazPrey for NaturalMortality Calculation--------------------------------------
            
                ZooGrazPrey = Me%ExternalVar%Mass(Diatoms, index)                                   &
                            + Me%ExternalVar%Mass(Phyto, index)                                     &
                            + Me%ExternalVar%Mass(CIL, index)                                    

            
            !LostChain, term due to different Phyto, Diatom and Bacteria and Zoo ratios----------- 
            !and LostGraz, term due to not 100% assimilation----------------- 
                            
                !Nitrogen
                if (Me%PropCalc%Nitrogen) then

                    LostChainNitrogen = (Me%Diatoms%DiaAlfaNC-Me%AlfaZooNC)                         &
                                        * (Me%Diatoms%DiaZooAssimilationRate * ZooIngestionDiatoms) &
                                       +(Me%AlfaPhytoNC-Me%AlfaZooNC)                               &
                                        * (Me%ZooAssimilationPhytoRate       * ZooIngestionPhyto)   & 
                                       +(Me%AlfaCilNC-Me%AlfaZooNC)                                 &
                                        * (Me%ZooAssimilationCiliateRate     * ZooIngestionCiliate)      
                
                    LostGrazNitrogen = (1- Me%Diatoms%DiaZooAssimilationRate)                       &
                                        * ZooIngestionDiatoms * Me%Diatoms%DiaAlfaNC                &
                                      +(1- Me%ZooAssimilationPhytoRate)                             &
                                        * ZooIngestionPhyto * Me%AlfaPhytoNC                        &
                                      +(1- Me%ZooAssimilationCiliateRate)                           &
                                        * ZooIngestionCiliate * Me%AlfaCilNC   

                    if (LostChainNitrogen .LT. 0.0)                                                 &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'
             
                    if (LostGrazNitrogen .LT. 0.0)                                                  &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'

                endif
                
                
                !Phosphorus
                if (Me%PropCalc%Phosphorus) then

                    LostChainPhosphorus = (Me%Diatoms%DiaAlfaPC-Me%AlfaZooPC)                       &
                                        * (Me%Diatoms%DiaZooAssimilationRate * ZooIngestionDiatoms) &
                                       +(Me%AlfaPhytoPC-Me%AlfaZooPC)                               &
                                        * (Me%ZooAssimilationPhytoRate       * ZooIngestionPhyto)   & 
                                       +(Me%AlfaCilPC-Me%AlfaZooPC)                                 &
                                        * (Me%ZooAssimilationCiliateRate     * ZooIngestionCiliate)      

                
                    LostGrazPhosphorus = (1- Me%Diatoms%DiaZooAssimilationRate)                     &
                                        * ZooIngestionDiatoms * Me%Diatoms%DiaAlfaPC                &
                                      +(1- Me%ZooAssimilationPhytoRate)                             &
                                        * ZooIngestionPhyto * Me%AlfaPhytoPC                        &
                                      +(1- Me%ZooAssimilationCiliateRate)                           &
                                        * ZooIngestionCiliate * Me%AlfaCilPC   

                    if (LostChainPhosphorus .LT. 0.0)                                               &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'
             
                    if (LostGrazPhosphorus .LT. 0.0)                                                &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'

                endif
                                       
        !endcase (ZooPhyDiaCil, ZooPhyDiaCilBac)
        case (ZooPhyCil, ZooPhyCilBac)

            !Zooplankton Grossgrowrate------------------------------------------------
                                                    
                !limitation by phyto (1º)
                x = Me%ZooEfficiencyCapturePhyto * Me%ExternalVar%Mass(Phyto, index)                & 
                - Me%GrazPhytoMin                                                                     
                                                                        
                y = Me%ZooIngestionConst + Me%ZooEfficiencyCapturePhyto                             &
                * Me%ExternalVar%Mass(Phyto, index) - Me%GrazPhytoMin                                 

                if  (x .LE. 0.0) then
                    ZooGrazPhytoLimitationFactor = 0.0
                else 
                    ZooGrazPhytoLimitationFactor = x / y
                end if

                !ZooIngestion of Phyto 
                ZooIngestionPhyto = Me%PhytoRatioIngestionZoo                                       &
                                  * (Me%ZooIngestionMax)                                            &                         
                                  * ZooGrazPhytoLimitationFactor                                    &
                                  * Me%TZooLimitationFactor                                         

                !limitation by Ciliate (2º)
                yac = (Me%ZooEfficiencyCaptureCiliate * Me%ExternalVar%Mass(CIL, index)             &
                     - Me%GrazCiliateMin)                                                               
                                                                                
                ybc = (Me%ZooIngestionConst + Me%ZooEfficiencyCaptureCiliate                        &
                    * Me%ExternalVar%Mass(CIL, index) - Me%GrazCiliateMin)

                if  (yac .LE. 0.0) then
                    ZooGrazCiliateLimitationFactor = 0.0
                else            
                    ZooGrazCiliateLimitationFactor = yac / ybc
                end if 

                !ZooIngestion of Ciliate
                ZooIngestionCiliate = Me%CiliatesRatioIngestionZoo                                  &
                                    * (Me%ZooIngestionMax - ZooIngestionPhyto)                      &
                                    * ZooGrazCiliateLimitationFactor                                &
                                    * Me%TZooLimitationFactor                                           
                                                                                
                !ZooGrossGrowRate, zooplankton gross growth rate                                        
                ZooGrossGrowRate =  (Me%ZooAssimilationPhytoRate   * ZooIngestionPhyto)             &
                                   +(Me%ZooAssimilationCiliateRate * ZooIngestionCiliate)                 
                         
                if (ZooGrossGrowRate .LT. 0.0) then                                                     
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR04.'
                !No sense on this. David
!                elseif (ZooGrossGrowRate .EQ. 0.0) then
!                    ZooGrossGrowRate =-1.0 / null_real   !Avoid division by zero below
                end if 
            
            !ZooGrazPrey for NaturalMortality Calculation--------------------------------------
            
                ZooGrazPrey = Me%ExternalVar%Mass(Phyto, index)                                     &
                            + Me%ExternalVar%Mass(CIL, index)                                       

            !LostChain, term due to different Phyto, Ciliates and Zoo ratios----------- 
            !and LostGraz, term due to not 100% assimilation----------------- 
                            
                !Nitrogen
                if (Me%PropCalc%Nitrogen) then

                    LostChainNitrogen = (Me%AlfaPhytoNC-Me%AlfaZooNC)                               &
                                        * (Me%ZooAssimilationPhytoRate       * ZooIngestionPhyto)   & 
                                       +(Me%AlfaCilNC-Me%AlfaZooNC)                                 &
                                        * (Me%ZooAssimilationCiliateRate     * ZooIngestionCiliate)      
                
                    LostGrazNitrogen = (1- Me%ZooAssimilationPhytoRate)                             &
                                        * ZooIngestionPhyto * Me%AlfaPhytoNC                        &
                                      +(1- Me%ZooAssimilationCiliateRate)                           &
                                        * ZooIngestionCiliate * Me%AlfaCilNC   

                    if (LostChainNitrogen .LT. 0.0)                                                 &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'
             
                    if (LostGrazNitrogen .LT. 0.0)                                                  &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'

                endif
                
                
                !Phosphorus
                if (Me%PropCalc%Phosphorus) then

                    LostChainPhosphorus = (Me%AlfaPhytoPC-Me%AlfaZooPC)                              &
                                         * (Me%ZooAssimilationPhytoRate       * ZooIngestionPhyto)   & 
                                        +(Me%AlfaCilPC-Me%AlfaZooPC)                                 &
                                         * (Me%ZooAssimilationCiliateRate     * ZooIngestionCiliate)      
                
                    LostGrazPhosphorus = (1- Me%ZooAssimilationPhytoRate)                           &
                                        * ZooIngestionPhyto * Me%AlfaPhytoPC                        &
                                      +(1- Me%ZooAssimilationCiliateRate)                           &
                                        * ZooIngestionCiliate * Me%AlfaCilPC   
             

                    if (LostChainPhosphorus .LT. 0.0)                                               &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'
             
                    if (LostGrazPhosphorus .LT. 0.0)                                                &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'

                endif

        !endcase (ZooPhyCil, ZooPhyCilBac)
        case (ZooDiaCilBac)

            !Zooplankton Grossgrowrate------------------------------------------------

                !limitation by diatoms (1º) 
                yac = (Me%Diatoms%ZooEfficiencyCaptureDiatoms * Me%ExternalVar%Mass(Diatoms, index) &
                     - Me%Diatoms%GrazDiaMin)                                                               
                                                                                
                ybc = (Me%ZooIngestionConst + Me%Diatoms%ZooEfficiencyCaptureDiatoms                &
                    * Me%ExternalVar%Mass(Diatoms, index) - Me%Diatoms%GrazDiaMin)

                if  (yac .LE. 0.0) then
                    ZooGrazDiatomsLimitationFactor = 0.0
                else           
                    ZooGrazDiatomsLimitationFactor = yac / ybc
                end if 

                !ZooIngestion of diatoms
                ZooIngestionDiatoms = Me%Diatoms%DiaRatioIngestionZoo                               &
                                    * Me%ZooIngestionMax                                            &
                                    * ZooGrazDiatomsLimitationFactor                                &
                                    * Me%TZooLimitationFactor                       
                                                    
                !limitation by Ciliate (2º) 
                yac = (Me%ZooEfficiencyCaptureCiliate * Me%ExternalVar%Mass(CIL, index)             &
                     - Me%GrazCiliateMin)                                                               
                                                                                
                ybc = (Me%ZooIngestionConst + Me%ZooEfficiencyCaptureCiliate                        &
                    * Me%ExternalVar%Mass(CIL, index) - Me%GrazCiliateMin)

                if  (yac .LE. 0.0) then
                    ZooGrazCiliateLimitationFactor = 0.0
                else           
                    ZooGrazCiliateLimitationFactor = yac / ybc
                end if 

                !ZooIngestion of Ciliate
                ZooIngestionCiliate = Me%CiliatesRatioIngestionZoo                                  &
                                    * (Me%ZooIngestionMax - ZooIngestionDiatoms)                    &
                                    * ZooGrazCiliateLimitationFactor                                &
                                    * Me%TZooLimitationFactor                                           
                                                                                
                !ZooGrossGrowRate, zooplankton gross growth rate                                        
                ZooGrossGrowRate =  (Me%Diatoms%DiaZooAssimilationRate * ZooIngestionDiatoms)       &
                                   +(Me%ZooAssimilationCiliateRate     * ZooIngestionCiliate)                 
                         
                if (ZooGrossGrowRate .LT. 0.0) then                                                     
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR07.'
                 !No sense on this. David
!                elseif (ZooGrossGrowRate .EQ. 0.0) then
!                    ZooGrossGrowRate =-1.0 / null_real   !Avoid division by zero below
                end if 
            
            !ZooGrazPrey for NaturalMortality Calculation--------------------------------------
            
                ZooGrazPrey = Me%ExternalVar%Mass(Diatoms, index)                                   &
                            + Me%ExternalVar%Mass(CIL, index)                                    

            !LostChain, term due to different  Diatoms, Ciliates and Zoo ratios----------- 
            !and LostGraz, term due to not 100% assimilation----------------- 
                            
                !Nitrogen
                if (Me%PropCalc%Nitrogen) then

                    LostChainNitrogen = (Me%Diatoms%DiaAlfaNC-Me%AlfaZooNC)                         &
                                        * (Me%Diatoms%DiaZooAssimilationRate * ZooIngestionDiatoms) &
                                       +(Me%AlfaCilNC-Me%AlfaZooNC)                                 &
                                        * (Me%ZooAssimilationCiliateRate     * ZooIngestionCiliate)      
                
                    LostGrazNitrogen = (1- Me%ZooAssimilationPhytoRate)                             &
                                        * ZooIngestionPhyto * Me%AlfaPhytoNC                        &
                                      +(1- Me%ZooAssimilationCiliateRate)                           &
                                        * ZooIngestionCiliate * Me%AlfaCilNC   

                    if (LostChainNitrogen .LT. 0.0)                                                 &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'
             
                    if (LostGrazNitrogen .LT. 0.0)                                                  &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'

                endif
                
                
                !Phosphorus
                if (Me%PropCalc%Phosphorus) then

                    LostChainPhosphorus = (Me%Diatoms%DiaAlfaPC-Me%AlfaZooPC)                       &
                                         * (Me%Diatoms%DiaZooAssimilationRate * ZooIngestionDiatoms)&
                                        +(Me%AlfaCilPC-Me%AlfaZooPC)                                &
                                         * (Me%ZooAssimilationCiliateRate     * ZooIngestionCiliate)      
                
                    LostGrazPhosphorus = (1- Me%Diatoms%DiaZooAssimilationRate)                     &
                                         * ZooIngestionDiatoms * Me%Diatoms%DiaAlfaPC               &
                                        +(1- Me%ZooAssimilationCiliateRate)                         &
                                         * ZooIngestionCiliate * Me%AlfaCilPC   
             

                    if (LostChainPhosphorus .LT. 0.0)                                               &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'
             
                    if (LostGrazPhosphorus .LT. 0.0)                                                &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'

                endif

        !endcase (ZooDiaCilBac)
        case (ZooPhyDia)

            !Zooplankton Grossgrowrate------------------------------------------------

                !limitation by diatoms (1º) 
                yac = (Me%Diatoms%ZooEfficiencyCaptureDiatoms * Me%ExternalVar%Mass(Diatoms, index) &
                     - Me%Diatoms%GrazDiaMin)                                                               
                                                                                
                ybc = (Me%ZooIngestionConst + Me%Diatoms%ZooEfficiencyCaptureDiatoms                &
                    * Me%ExternalVar%Mass(Diatoms, index) - Me%Diatoms%GrazDiaMin)

                if  (yac .LE. 0.0) then
                    ZooGrazDiatomsLimitationFactor = 0.0
                else           
                    ZooGrazDiatomsLimitationFactor = yac / ybc
                end if 

                !ZooIngestion of diatoms
                ZooIngestionDiatoms = Me%Diatoms%DiaRatioIngestionZoo                               &
                                    * Me%ZooIngestionMax                                            &
                                    * ZooGrazDiatomsLimitationFactor                                &
                                    * Me%TZooLimitationFactor                       
                                                    
                !limitation by phyto (2º)
                x = Me%ZooEfficiencyCapturePhyto * Me%ExternalVar%Mass(Phyto, index)                & 
                - Me%GrazPhytoMin                                                                     
                                                                        
                y = Me%ZooIngestionConst + Me%ZooEfficiencyCapturePhyto                             &
                * Me%ExternalVar%Mass(Phyto, index) - Me%GrazPhytoMin                                 

                if  (x .LE. 0.0) then
                    ZooGrazPhytoLimitationFactor = 0.0
                else 
                    ZooGrazPhytoLimitationFactor = x / y
                end if

                !ZooIngestion of Phyto 
                ZooIngestionPhyto = Me%PhytoRatioIngestionZoo                                       &
                                  * (Me%ZooIngestionMax -  ZooIngestionDiatoms)                     &                         
                                  * ZooGrazPhytoLimitationFactor                                    &
                                  * Me%TZooLimitationFactor                                                                        
                                                                                
                !ZooGrossGrowRate, zooplankton gross growth rate                                        
                ZooGrossGrowRate =  (Me%Diatoms%DiaZooAssimilationRate * ZooIngestionDiatoms)       &
                                   +(Me%ZooAssimilationPhytoRate       * ZooIngestionPhyto  )       
                         
                if (ZooGrossGrowRate .LT. 0.0) then 
                    stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR10.'
                !No sense on this. David
!                elseif (ZooGrossGrowRate .EQ. 0.0) then
!                    ZooGrossGrowRate =-1.0 / null_real   !Avoid division by zero below
                end if 
            
            
            !ZooGrazPrey for NaturalMortality Calculation--------------------------------------
            
                ZooGrazPrey = Me%ExternalVar%Mass(Diatoms, index)                                   &
                            + Me%ExternalVar%Mass(Phyto, index)                                     


            !LostChain, term due to different Phyto, Diatoms and Zoo ratios----------- 
            !and LostGraz, term due to not 100% assimilation----------------- 
                            
                !Nitrogen
                if (Me%PropCalc%Nitrogen) then

                    LostChainNitrogen = (Me%Diatoms%DiaAlfaNC-Me%AlfaZooNC)                         &
                                        * (Me%Diatoms%DiaZooAssimilationRate * ZooIngestionDiatoms) &
                                       +(Me%AlfaPhytoNC-Me%AlfaZooNC)                               &
                                        * (Me%ZooAssimilationPhytoRate       * ZooIngestionPhyto)    
                
                    LostGrazNitrogen = (1- Me%Diatoms%DiaZooAssimilationRate)                       &
                                        * ZooIngestionDiatoms * Me%Diatoms%DiaAlfaNC                &
                                      +(1- Me%ZooAssimilationPhytoRate)                             &
                                        * ZooIngestionPhyto * Me%AlfaPhytoNC                        

                    if (LostChainNitrogen .LT. 0.0)                                                 &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'
             
                    if (LostGrazNitrogen .LT. 0.0)                                                  &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'

                endif
                
                
                !Phosphorus
                if (Me%PropCalc%Phosphorus) then

                    LostChainPhosphorus = (Me%Diatoms%DiaAlfaPC-Me%AlfaZooPC)                       &
                                        * (Me%Diatoms%DiaZooAssimilationRate * ZooIngestionDiatoms) &
                                       +(Me%AlfaPhytoPC-Me%AlfaZooPC)                               &
                                        * (Me%ZooAssimilationPhytoRate       * ZooIngestionPhyto)    
                
                    LostGrazPhosphorus = (1- Me%Diatoms%DiaZooAssimilationRate)                     &
                                        * ZooIngestionDiatoms * Me%Diatoms%DiaAlfaPC                &
                                      +(1- Me%ZooAssimilationPhytoRate)                             &
                                        * ZooIngestionPhyto * Me%AlfaPhytoPC                        
             

                    if (LostChainPhosphorus .LT. 0.0)                                               &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'
             
                    if (LostGrazPhosphorus .LT. 0.0)                                                &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'

                endif

        !endcase (ZooPhyDia)
        case (ZooPhy)

            !Zooplankton Grossgrowrate------------------------------------------------
                
                if ((Me%ExternalVar%Mass(Phyto, index) - Me%GrazPhytoMin) .LE. 0.0) then
                    
                    FoodZooLimitationFactor = 0.0
                
                else 

                    exponent = -Me%IvlevGrazConst                                                   &
                             * (Me%ExternalVar%Mass(Phyto, index) - Me%GrazPhytoMin)
                    
                    FoodZooLimitationFactor = 1.0 - exp(exponent)

                end if 

                ZooGrossGrowRate = Me%GrowMaxZooRate                                                &
                                 * Me%TZooLimitationFactor * FoodZooLimitationFactor
                         
                if (ZooGrossGrowRate .LT. 0.0) then
                
                    write(*,*)'ZooGrossGrowRate set to zero     :', ZooGrossGrowRate
                    write(*,*)'Me%GrowMaxZooRate                :', Me%GrowMaxZooRate 
                    write(*,*)'Me%TZooLimitationFactor          :', Me%TZooLimitationFactor
                    write(*,*)'FoodZooLimitationFactor          :', FoodZooLimitationFactor
                    write(*,*)'Me%ExternalVar%Mass(Phyto, index):', Me%ExternalVar%Mass(Phyto, index)
                    write(*,*)'exponent                         :', exponent
                    write(*,*)'index                            :', index
                    
                    ZooGrossGrowRate = 0.0 
                    
!                    ZooGrossGrowRate =-1.0 / null_real   !Avoid division by zero below
                                                                     
                    !stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR13.'
                
                !No sense on this. David
!                elseif (ZooGrossGrowRate .EQ. 0.0) then
!                    ZooGrossGrowRate =-1.0 / null_real   !Avoid division by zero below
                end if 
            
           
            !ZooGrazPrey for NaturalMortality Calculation--------------------------------------
            
                ZooGrazPrey = Me%ExternalVar%Mass(Phyto, index)                                   

 
            !LostChain, term due to different Phyto, Diatoms and Zoo ratios----------- 
            !and LostGraz, term due to not 100% assimilation----------------- 
                            
                !Nitrogen
                if (Me%PropCalc%Nitrogen) then

                    LostChainNitrogen = (Me%AlfaPhytoNC-Me%AlfaZooNC) * ZooGrossGrowRate      
                

                    LostGrazNitrogen = ZooGrossGrowRate                                             &
                                      * Me%AlfaPhytoNC * ((1.0 - Me%E) / Me%E)   

                    if (LostChainNitrogen .LT. 0.0)                                                 &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'
             
                    if (LostGrazNitrogen .LT. 0.0)                                                  &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'

                endif
                
                
                !Phosphorus
                if (Me%PropCalc%Phosphorus) then

                    LostChainPhosphorus = (Me%AlfaPhytoPC-Me%AlfaZooPC) * ZooGrossGrowRate          
                
                    LostGrazPhosphorus = ZooGrossGrowRate                                           &
                                        * Me%AlfaPhytoPC * ((1.0 - Me%E) / Me%E)                                 

                    if (LostChainPhosphorus .LT. 0.0)                                               &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'
             
                    if (LostGrazPhosphorus .LT. 0.0)                                                &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'

                endif

        !endcase (ZooPhy)
        case (ZooDia)

            !Zooplankton Grossgrowrate------------------------------------------------
                
                if ((Me%ExternalVar%Mass(Diatoms, index) - Me%Diatoms%GrazDiaMin) .LE. 0.0) then
                    
                    FoodZooLimitationFactor = 0.0
                
                else 

                    exponent = -Me%IvlevGrazConst                                                   &
                             * (Me%ExternalVar%Mass(Diatoms, index) - Me%Diatoms%GrazDiaMin)
                    
                    FoodZooLimitationFactor = 1.0 - exp(exponent)

                end if 

                ZooGrossGrowRate = Me%GrowMaxZooRate                                                &
                                 * Me%TZooLimitationFactor * FoodZooLimitationFactor
                         
                if (ZooGrossGrowRate .LT. 0.0) then                                                     
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR16.'
                    
                 !No sense on this. David
!                elseif (ZooGrossGrowRate .EQ. 0.0) then
!                    ZooGrossGrowRate =-1.0 / null_real   !Avoid division by zero below
                end if 
            
           
            !ZooGrazPrey for NaturalMortality Calculation--------------------------------------
            
                ZooGrazPrey = Me%ExternalVar%Mass(Diatoms, index)                                   

            !LostChain, term due to different Phyto, Diatoms, Ciliates and Zoo ratios----------- 
            !and LostGraz, term due to not 100% assimilation----------------- 

                !Nitrogen
                if (Me%PropCalc%Nitrogen) then

                    LostChainNitrogen = (Me%Diatoms%DiaAlfaNC-Me%AlfaZooNC) * ZooGrossGrowRate
                
                    LostGrazNitrogen = ZooGrossGrowRate                                             &
                                      * Me%Diatoms%DiaAlfaNC * ((1.0 - Me%Diatoms%DiaE) / Me%Diatoms%DiaE)                

                    if (LostChainNitrogen .LT. 0.0)                                                 &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'
             
                    if (LostGrazNitrogen .LT. 0.0)                                                  &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'

                endif
                
                
                !Phosphorus
                if (Me%PropCalc%Phosphorus) then

                    LostChainPhosphorus = (Me%Diatoms%DiaAlfaPC-Me%AlfaZooPC) * ZooGrossGrowRate
                
                    LostGrazPhosphorus = ZooGrossGrowRate                                           &
                                        * Me%Diatoms%DiaAlfaPC * ((1.0 - Me%Diatoms%DiaE) / Me%Diatoms%DiaE)              

                    if (LostChainPhosphorus .LT. 0.0)                                               &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'
             
                    if (LostGrazPhosphorus .LT. 0.0)                                                &
                        stop 'Subroutine WQZooplankton-Module ModuleWaterQuality. ERR02.'

                endif

        !endcase (ZooDia)
        end select

        !ZooNaturalMortality/ExcretionRate/Zooplankton predation----------------------

            !ZooNaturalMortalityRate
            if (ZooGrazPrey .LE. Me%GrazPreyMin ) then
                ZooNaturalMortality = Me%ZooMaxMortalityRate
            else                                                                                    
                ZooNaturalMortality = (Me%ZooMortalityCoef / ZooGrazPrey)                       & 
                                    + Me%ZooMinMortalityRate
            end if  
            
            !Zooplankton ExcretionRate -------------------------------------------------------
            ZooExcretionRate = Me%ZooExcretionFactor * Me%ZooExcretionConst                         &
                               ** Me%ExternalVar%Temperature(index)

            !Nitrogen
            if (Me%PropCalc%Nitrogen) then
            
                !DeadZooNitrogen, Organic Nitrogen from dead Zooplankton (mgN/L)
                DeadZooNitrogen   = ZooNaturalMortality * Me%AlfaZooNC

                !Zooplankton predated by higher trophic levels---------------------------------
                PredatedZooNitrogen   = Me%ZPredMortalityRate * Me%AlfaZooNC

                !ZooExcretionNitrogen (mg N/L)
                ZooExcretionNitrogen = ZooExcretionRate * Me%AlfaZooNC 

            endif
            
            !Phosphorus
            if (Me%PropCalc%Phosphorus) then
            
                !DeadZooPhosphorus, Organic Phosphorus from dead Zooplankton (mgP/L)
                DeadZooPhosphorus = ZooNaturalMortality * Me%AlfaZooPC 

                !Zooplankton predated by higher trophic levels---------------------------------
                PredatedZooPhosphorus = Me%ZPredMortalityRate * Me%AlfaZooPC

                !ZooExcretionPhosphorus (mg P/L)
                ZooExcretionPhosphorus = ZooExcretionRate * Me%AlfaZooPC 

            endif

        !Zooplankton RespirationRate -------------------------------------------------------
        
            ZooRespirationRate = Me%ZooReferenceRespirationRate * Me%TZooLimitationFactor
       
            !Oxygen loss due to Zooplankton respiration
            if (Me%PropCalc%Oxygen)  then
                
                OxygenZooRespRate = Me%RatioOxygenCarbonZooRespiration * ZooRespirationRate
            
                if (MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen).EQ.Me%MinOxygen) then
                    ZooGrossGrowRate  = -1.0 / null_real
                    OxygenZooRespRate = -1.0 / null_real
                endif
            
            endif 

           
    !Calculation of system coeficients---------------------------------------

    !Organisms---------------------------------------
    
    Me%Matrix(Zoo, Zoo)     = 1.0 + DTDay * ( ZooExcretionRate                                      & 
                                                     + ZooRespirationRate                           & 
                                                     + ZooNaturalMortality                          &                     
                                                     + Me%ZPredMortalityRate                        &
                                                     - ZooGrossGrowRate)
   
    select case (Me%WQConfiguration)

        case (ZooPhy)

            Me%Matrix(Phyto, Zoo)   = (DTDay * ZooGrossGrowRate) / Me%E

        case (ZooDia)

            Me%Matrix(Diatoms, Zoo)=  (DTDay * ZooGrossGrowRate) / Me%Diatoms%DiaE

            Me%Matrix(BioSi, Zoo)  = -(DTDay * ZooGrossGrowRate) / Me%Diatoms%DiaE * Me%Diatoms%DiaAlfaSiC

        case default  
       
            if (Me%PropCalc%Phyto)then
               
                Me%Matrix(Phyto, Zoo)   =   DTDay * ZooIngestionPhyto
            
            endif


            if (Me%PropCalc%Diatoms)then
               
                Me%Matrix(Diatoms, Zoo) =   DTDay * ZooIngestionDiatoms

                Me%Matrix(BioSi, Zoo)   = - DTDay * ZooIngestionDiatoms * Me%Diatoms%DiaAlfaSiC
            
            endif

            if (Me%PropCalc%Ciliate)then
               
                Me%Matrix(CIL, Zoo)     =   DTDay * ZooIngestionCiliate  
            
            endif

    end select
 
    !Oxygen---------------------------------------
   
        if (Me%PropCalc%Oxygen)                                                                     &
            Me%Matrix(O,  Zoo)     =   DTDay * OxygenZooRespRate   

    !Nitrogen---------------------------------------
        if (Me%PropCalc%Nitrogen) then
        
            Me%Matrix(AM, Zoo   ) = - DTDay * (ZooExcretionNitrogen * Me%ZooSolublInorgExcreFraction &
                                              +ZooRespirationRate   * Me%AlfaZooNC)    

            Me%Matrix(DONnr, Zoo) = - DTDay * (ZooExcretionNitrogen                                 & 
                                               * Me%ZooExcreDissOrgFraction                         & 
                                               * (1.0 - Me%ZooSolublInorgExcreFraction))           
            
            Me%Matrix(PON, Zoo  ) = - DTDay * (ZooExcretionNitrogen                                 & 
                                               * (1.0 - Me%ZooExcreDissOrgFraction)                 & 
                                               * (1.0 - Me%ZooSolublInorgExcreFraction)             &
                                              + DeadZooNitrogen                                     &
                                              + PredatedZooNitrogen                                 &
                                              + LostChainNitrogen                                   &
                                              + LostGrazNitrogen)
        endif

    !Phosphorus---------------------------------------
        if (Me%PropCalc%Phosphorus) then
        
            Me%Matrix(IP, Zoo   ) = - DTDay * (ZooExcretionPhosphorus * Me%ZooSolublInorgExcreFraction &
                                              + ZooRespirationRate   * Me%AlfaZooPC) 
            
            
            Me%Matrix(DOPnr, Zoo) = -DTDay * (ZooExcretionPhosphorus                                &
                                              * Me%ZooExcreDissOrgFraction                          & 
                                              * (1.0 - Me%ZooSolublInorgExcreFraction))           

            Me%Matrix(POP, Zoo  ) = -DTDay * (ZooExcretionPhosphorus                                &
                                               * (1.0 - Me%ZooExcreDissOrgFraction)                 & 
                                               * (1.0 - Me%ZooSolublInorgExcreFraction)             &
                                              + DeadZooPhosphorus                                   &
                                              + PredatedZooPhosphorus                               &
                                              + LostChainPhosphorus                                 &
                                              + LostGrazPhosphorus)
        endif

    !Independent term---------------------------------------
        Me%IndTerm(Zoo) = Me%ExternalVar%Mass(Zoo, index) 

    end subroutine WQZooplankton

    
    !----------------------------------------------------------------------------

   
    subroutine WQLarvae(index) ! Aires Santos and João Nogueira

    
    !   This routine computes growth and mortality rates for sardine Larvae. This
    !   rates include temperature, fishfood and weight dependencie parameteres.
    !
        
        
        
        !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index


        !Local-------------------------------------------------------------------

        integer :: Larvae
        integer :: Age
        integer :: NPhases

        real    :: Age_
        real    :: Temperature
        real    :: DTDay
        real    :: CLarvae
        real    :: WLarvae
        real    :: ZLarvae
        real    :: ZLarvaeWeight
        real    :: GLarvaeWeight
        real    :: GLarvae
        real    :: Gfishfood
        real    :: FishFood

        real    :: KGtemp
        real    :: KZtemp
        real    :: KGfishfood
        real    :: Gfishfood_ref
        real    :: GLarvaeTemp_ref
        real    :: ZLarvaeTemp_ref
        real    :: Afg
        real    :: Awg
        real    :: Bwg
        real    :: Awz
        real    :: Bwz
        real    :: Atg
        real    :: Btg
        real    :: Atz
        real    :: Btz
        real    :: Ldensity
        real    :: Lshape
        real    :: Init_age
        real    :: Inter_age
        real    :: Final_age
        real    :: Init_length
        real    :: Inter_length
        real    :: Final_length
        real    :: FishFood_ref
        real    :: Temperature_ref

        !------------------------------------------------------------------------

        Larvae      = Me%PropIndex%Larvae
        Age         = Me%PropIndex%Age
        Age_        = Me%ExternalVar%Mass(Age,index)
        DTDay       = Me%DTDay
        Temperature = Me%ExternalVar%Temperature(index)
        FishFood    = Me%ExternalVar%FishFood(index)

        !Calculation of the temperature coefficients-----------------------------

        Temperature_ref = Me%Temperature_ref                                                   ! ºC
        Atg = Me%Atg
        Btg = Me%Btg
        Atz = Me%Atz
        Btz = Me%Btz

        GLarvaeTemp_ref = Atg + Btg * Temperature_ref                     ! Houde e Zastrow
        ZLarvaeTemp_ref = Atz + Btz * Temperature_ref                     ! Houde e Zastrow

        KGtemp = 1.+(Btg/GLarvaeTemp_ref)*(Temperature-Temperature_ref)
        KZtemp = 1.+(Btz/ZLarvaeTemp_ref)*(Temperature-Temperature_ref)

        !Calculation of the nutrient coefficients (Michaelis-Menten equation)-------------

        Afg = Me%Afg                                          
        FishFood_ref = Me%FishFood_ref                                                      ! mg/m3
        
        Gfishfood_ref = FishFood_ref/(Afg + FishFood_ref)  
        
        Gfishfood = FishFood/(Afg + FishFood)  
        KGfishfood = Gfishfood/Gfishfood_ref


        !Calculation of the mass growth rate------------------------------------- 

        Init_age = Me%Init_age
        Inter_age = Me%Inter_age
        Final_age = Me%Final_age
        Init_length = Me%Init_length
        Inter_length = Me%Inter_length
        Final_length = Me%Final_length
        NPhases = Me%NPhases

cd1 :       if (NPhases .EQ. 1) then

cd2 :           if (Age_ .le. Final_age) then
                    CLarvae = Init_length + ((Final_length - Init_length) / Final_age) * Age_
                else
                    write(*,*)'Larvae too OLD....'
                    stop      'WQLarvae - ModuleWaterQuality - ERR01'
                endif cd2
                
            elseif (NPhases .EQ. 2) then

cd3 :           if (Age_ .le. Inter_age) then
                    CLarvae = Init_length + ((Inter_length - Init_length) / Inter_age) * Age_
                elseif (Age_ .gt. Inter_age .and. Age_ .le. Final_age) then
                    CLarvae = Init_length + ((Final_length - Inter_length) / (Final_age - Inter_age)) * (Age_ - Inter_age)
                else
                    write(*,*)'Larvae too OLD....'
                    stop      'WQLarvae - ModuleWaterQuality - ERR01'
                endif cd3

            else 
                write(*,*)'Valid values for NPHASES are only 1 or 2'
                stop      'WQLarvae - ModuleWaterQuality - ERR02'
            endif cd1

        Ldensity = Me%Ldensity
        Lshape = Me%Lshape
        Awg = Me%Awg   
        Bwg = Me%Bwg 
         
        WLarvae         = exp ((-1 * Ldensity) + Lshape * log(CLarvae))                ! Girão
        GLarvaeWeight   = Awg / (WLarvae**(Bwg))                                ! Houde
        GLarvae         = KGfishfood * KGtemp * GLarvaeWeight 


        !Calculation of the mortality rate--------------------------------------- 

        Awz = Me%Awz   
        Bwz = Me%Bwz         
        
        ZLarvaeWeight  = Awz * (WLarvae**(-1 * Bwz))                             ! Houde
        ZLarvae        = KZtemp * ZLarvaeWeight


         
        !Calculation of system coeficients---------------------------------------
        Me%Matrix(Larvae, Larvae) = 1.0 + DTDay * ZLarvae 


        !Independent term
        Me%IndTerm(Larvae) = Me%ExternalVar%Mass(Larvae, index)       &
                                          * (1+DTDay*GLarvae)


   !------------------------------------------------------------------------

    end subroutine WQLarvae

    !----------------------------------------------------------------------------


    subroutine WQAge(index) ! Aires !index é as particulas


        !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

        !Local-------------------------------------------------------------------

        integer :: Age

        real    :: DTDay

        !------------------------------------------------------------------------

        Age     = Me%PropIndex%Age 

        DTDay   = Me%DTDay

 
        !Calculation of system coeficients---------------------------------------
        Me%Matrix(Age, Age) = 1.0  


        !Independent term
        Me%IndTerm(Age) = Me%ExternalVar%Mass(Age, index) + DTDay

    !------------------------------------------------------------------------

    end subroutine WQAge

    !----------------------------------------------------------------------------


    subroutine WQPhytoplankton(index)

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        
        integer :: Phyto
        integer :: IP
        integer :: POP, DOPnr
        integer :: AM
        integer :: PON
        integer :: O
        integer :: NA
        integer :: DONnr
        integer :: BOD


        real    :: DTDay

        real    :: PhytoExcretionsRate              = null_real
        real    :: PhytoExcretionPhosphorus         = null_real
        real    :: PhytoExcretionNitrogen           = null_real
        real    :: DeadPhytoNitrogen                = null_real
        real    :: DeadPhytoPhosphorus
        real    :: PhytoTotalRespirationLossesRate  = null_real  
        real    :: NitrateOxygen                    = null_real
        real    :: IPOxygen                         = null_real !Rosa
        
        real    :: PhytoAmmoniaPreferenceFactor     = null_real
        real    :: PhytoEndogenousRepiration        = null_real
        real    :: PhotorespirationRate             = null_real
        real    :: PhotoOxygen                      = null_real
        real    :: OxygenPhytoRespRate              = null_real
        real    :: PhytoNonGrazingMortalityRate     = null_real

        real    :: s1, s2, xa, xb, ya, yb
        real    :: x1, x2, x3 ,x4

    !------------------------------------------------------------------------


        POP     = Me%PropIndex%PartOrganicPhosphorus
        DOPnr   = Me%PropIndex%DOPNonRefractory
        Phyto   = Me%PropIndex%Phyto
        IP      = Me%PropIndex%InorganicPhosphorus
        AM      = Me%PropIndex%Ammonia
        NA      = Me%PropIndex%Nitrate
        PON     = Me%PropIndex%PartOrganicNitrogen
        O       = Me%PropIndex%Oxygen
        DONnr   = Me%PropIndex%DONNonRefractory
        BOD     = Me%PropIndex%BOD
        DTDay   = Me%DTDay


        s1                              = null_real
        s2                              = null_real
        xa                              = null_real
        xb                              = null_real
        ya                              = null_real
        yb                              = null_real

        x1                              = null_real
        x2                              = null_real
        x3                              = null_real
        x4                              = null_real


    !------------------------------------------------------------------------

    !PhytoGrossGrowRate, phytoplankton gross grow rate-----------------------
    !PhytoNutrientsLimitationFactor, the nutrients limitation factor is the minimum between 
    !Nitrogen and Phosphorus
cd2 :   if (Me%PropCalc%Nitrogen) then
            Me%PhytoNLimitationFactor   = (Me%ExternalVar%Mass(AM, index)                      &
                                                + Me%ExternalVar%Mass(NA, index))                     &
                                                / (Me%NSatConst                                       &
                                                + Me%ExternalVar%Mass(AM, index)                      &
                                                + Me%ExternalVar%Mass(NA, index))
        else
            Me%PhytoNLimitationFactor   = 1.0   !Nitrogen is not a limiting nutrient
        end if cd2




    !Phosphorus limitation factor
cd3 :  if (Me%PropCalc%Phosphorus.AND.(.NOT.Me%PropCalc%Bacteria)) then
            Me%PhytoPLimitationFactor =  Me%ExternalVar%Mass(IP, index)                      &
                                               / (Me%PSatConst                                        &
                                               + Me%ExternalVar%Mass(IP, index))
        else
            Me%PhytoPLimitationFactor = 1.0   !Phosphorus is not a limiting nutrient
       end if cd3
        

        Me%PhytoNutrientsLimitationFactor = &
                              min(Me%PhytoNLimitationFactor, Me%PhytoPLimitationFactor)


    !TPhytoLimitationFactor, temperature effect on Phytoplankton assimilation rate
        s1 = (1.0 / (Me%TOptPhytoMin - Me%TPhytoMin)) * log((Me%FK2 * (1.0 - Me%FK1))              &
                                                          / (Me%FK1 * (1.0 - Me%FK2)))

        s2 = (1.0 / (Me%TPhytoMax - Me%TOptPhytoMax)) * log((Me%FK3 * (1.0 - Me%FK4))              &
                                                          / (Me%FK4 * (1.0 - Me%FK3)))

        ya = exp(s1 * (Me%ExternalVar%Temperature(index) - Me%TPhytoMin))
        yb = exp(s2 * (Me%TPhytoMax - Me%ExternalVar%Temperature(index)))

        xa = (Me%FK1 * ya) / (1.0 + Me%FK1 * (ya - 1.0))
        xb = (Me%FK4 * yb) / (1.0 + Me%FK4 * (yb - 1.0))

        Me%TPhytoLimitationFactor = xa * xb


    !PhytoLightLimitationFactor, light effect on Phytoplankton assimilation rate

        Me%PhytoLightLimitationFactor =                                                            &
          PhytoLightLimitationFactor(Thickness       = Me%ExternalVar%Thickness(index),            &
                                     TopRadiation    = Me%ExternalVar%ShortWaveRadiation(index),   &
                                     PExt            = Me%ExternalVar%LightExtCoefField(index),    &
                                     Photoinhibition = Me%Photoinhibition)

    !PhytoGrossGrowRate
        Me%PhytoGrossGrowRate = Me%GrowMaxPhytoRate                                                &
                                * Me%TPhytoLimitationFactor                                        &
                                * Me%PhytoLightLimitationFactor                                    &
                                * Me%PhytoNutrientsLimitationFactor


cd21 :  if     (Me%PhytoGrossGrowRate .LT. 0.0) then
            write(*,*)'Me%PhytoGrossGrowRate set to zero:', Me%PhytoGrossGrowRate
            write(*,*)'Me%TPhytoLimitationFactor        :', Me%TPhytoLimitationFactor
            write(*,*)'Me%PhytoLightLimitationFactor    :', Me%PhytoLightLimitationFactor
            write(*,*)'Me%PhytoNutrientsLimitationFactor:', Me%PhytoNutrientsLimitationFactor
            write(*,*)'Me%ExternalVar%Mass(AM, index)   :', Me%ExternalVar%Mass(AM, index)
            write(*,*)'Me%ExternalVar%Mass(NA, index)   :', Me%ExternalVar%Mass(NA, index)
            write(*,*)'Me%ExternalVar%Mass(IP, index)   :', Me%ExternalVar%Mass(IP, index)
            
            !this does not makes sense, crashing the rates. David
!            Me%PhytoGrossGrowRate =-1.0 / null_real   !Avoid division by zero below
!            !stop 'Subroutine WQPhytoplankton - ModuleWaterQuality. ERR01.'
            Me%PhytoGrossGrowRate = 0.0
            
        !this does not makes sense, crashing the rates if one limiting factor is zero
        !(e.g. at night where light limiting factor is zero)
        !Division by zero is avoided by putting the death rate in maximum. David
!        else if (Me%PhytoGrossGrowRate .EQ. 0.0) then cd21
!            Me%PhytoGrossGrowRate =-1.0 / null_real   !Avoid division by zero below
        end if cd21

        if (MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen).eq.Me%MinOxygen) then
            Me%PhytoGrossGrowRate = -1.0 / null_real
        endif


    !PhytoAmmoniaPreferenceFactor
cd22 :  if (Me%PropCalc%Nitrogen) then
          
            x1 = Me%ExternalVar%Mass(AM, index) * Me%ExternalVar%Mass(NA, index)

            x2 = (Me%NSatConst + Me%ExternalVar%Mass(AM, index))                                   &
               * (Me%NSatConst + Me%ExternalVar%Mass(NA, index)) 

            x3 = Me%NSatConst * Me%ExternalVar%Mass(AM, index)

            x4 = (Me%ExternalVar%Mass(AM, index) + Me%ExternalVar%Mass(NA, index))                 &
               * (Me%NSatConst + Me%ExternalVar%Mass(NA, index))

cd45 :      if ((x1 .EQ. 0.0) .AND. (x3 .EQ. 0.0)) then
                PhytoAmmoniaPreferenceFactor = 0.0                 
            else cd45
                PhytoAmmoniaPreferenceFactor = (x1 / x2) + (x3 / x4)
            end if cd45
        end if cd22


    !Respiration losses----------------------------------------
    
        !PhytoEndogenousRepiration, endogenous respiration
            PhytoEndogenousRepiration = Me%PhytoEndogRepConst *                                   &
                                        exp(0.069 * Me%ExternalVar%Temperature(index))
      
            if (MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen).eq.Me%MinOxygen) then
                PhytoEndogenousRepiration = -1.0 / null_real
            endif
        
        !PhotorespirationRate
            PhotorespirationRate = Me%PhotorespFactor * Me%PhytoGrossGrowRate


        !PhytoTotalRespirationLosses, total respiration losses
            PhytoTotalRespirationLossesRate = PhytoEndogenousRepiration + PhotorespirationRate


    !DeadPhyto, phyto mortality from non grazing
        !previous version with absurd values if gross rate was zero (if one limiting factor zero)
        !makes no sense, crashing rates. 
        !The below equation, as gross rate decreases, tends to max mortality. Only at very 
        !low phyto conc (<1E-10) and very low gross rates, the below equation would generate mortality 
        !rates lower than max mortality so it makes sense to define max mortality at zero growth rate
        !avoiding division by zero in a more elegant way and not crashing rates. David
        if (Me%PhytoGrossGrowRate .EQ. 0.0) then
            PhytoNonGrazingMortalityRate = Me%PhytoMortMaxRate
        
        else
            PhytoNonGrazingMortalityRate = Me%PhytoMortMaxRate *                                       &
                                          ((Me%ExternalVar%Mass(Phyto, index)                          &
                                            / Me%PhytoGrossGrowRate)                                   &
                                           / (Me%FMortSatConst                                         &
                                              + (Me%ExternalVar%Mass(Phyto, index)                     &  
                                              / Me%PhytoGrossGrowRate)))

        endif

    !Excretion and Dead losses----------------------------------------

        PhytoExcretionsRate = Me%PhytoExcretionConstant * Me%PhytoGrossGrowRate                    &
                            * (1.0 - Me%PhytoLightLimitationFactor)

cd12 :  if (Me%PropCalc%Phosphorus) then


            PhytoExcretionPhosphorus = Me%AlfaPhytoPC                                              &
                                     * (PhytoTotalRespirationLossesRate                            &
                                     + PhytoExcretionsRate)
       
            !aqui
            !DeadPhyto, Organic Phosphorus from dead Phytoplankton
            DeadPhytoPhosphorus = PhytoNonGrazingMortalityRate * Me%AlfaPhytoPC
     
        end if cd12


cd13 :  if (Me%PropCalc%Nitrogen) then
            
            PhytoExcretionNitrogen =  Me%AlfaPhytoNC                                               &
                                   * (PhytoTotalRespirationLossesRate                              &
                                   + PhytoExcretionsRate)

            !DeadPhyto, Organic Nitrogen from dead Phytoplankton
            DeadPhytoNitrogen = PhytoNonGrazingMortalityRate * Me%AlfaPhytoNC 
        
        end if cd13

    !Oxygen Cycle----------------------------------------

cd6 :   if (Me%PropCalc%Oxygen) then
          
            if (MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen).eq.Me%MinOxygen) then
                NitrateOxygen = -1.0 / null_real
                IPOxygen = -1.0/null_real
                PhotoOxygen = -1.0 / null_real
                OxygenPhytoRespRate = -1.0 / null_real
            else 

               !Phytoplankton Oxygen production, NitrateOxygen + IPOxygen + PhotoOxygen
cd66 :          if (Me%PropCalc%Nitrogen) then
      
                    !Secondary Oxygen production due to Nitrate consumption
                    NitrateOxygen = Me%PhytoGrossGrowRate * Me%NConsOxyNitRatio                    &
                                  * Me%AlfaPhytoNC * (1 - PhytoAmmoniaPreferenceFactor)      

                else
                   
                    NitrateOxygen = 0.0

                end if cd66               

!Rosa
cd666:          if (Me%PropCalc%Phosphorus) then 

                    !Secondary Oxygen production due to Inorganic phosphorus (phosphate) consumption
                    IPOxygen = Me%PhytoGrossGrowRate * Me%PConsOxyPhosphorusRatio * Me%AlfaPhytoPC     
                else
            
                    IPOxygen = 0.0
            
                end if cd666               

                !Oxygen photosynthetic production
                PhotoOxygen = Me%PhytoGrossGrowRate * Me%PhotosynthesisOxygenCarbonRatio


                !Oxygen loss by respiration     
                OxygenPhytoRespRate = Me%PlanktonOxygenCarbonRatio * PhytoTotalRespirationLossesRate
            end if
        
        end if cd6


        !Calculation of system coeficients---------------------------------------
             Me%Matrix(Phyto, Phyto) = 1.0 + DTDay * (PhytoTotalRespirationLossesRate              &
                                                     + PhytoExcretionsRate                         &
                                                     + PhytoNonGrazingMortalityRate                &
                                                     - Me%PhytoGrossGrowRate)

cd4 :   if (Me%PropCalc%Phosphorus) then


            Me%Matrix(IP, Phyto   ) =  DTDay * (Me%AlfaPhytoPC * Me%PhytoGrossGrowRate              &
                                            - (PhytoExcretionPhosphorus * Me%PhytoSolublInorgExcreFraction))


            Me%Matrix(DOPnr, Phyto) = -DTDay * PhytoExcretionPhosphorus                            &
                                             * (1.0 - Me%PhytoSolublInorgExcreFraction)            &
                                             * Me%PhytoExcreDissOrgFraction

            Me%Matrix(POP, Phyto  ) = -DTDay * ((PhytoExcretionPhosphorus                          & 
                                               * (1.0 - Me%PhytoSolublInorgExcreFraction)          &
                                               * (1.0 - Me%PhytoExcreDissOrgFraction))             &
                                               + DeadPhytoPhosphorus)                          
                    
        end if cd4


cd5 :   if (Me%PropCalc%Nitrogen) then

            Me%Matrix(AM, Phyto   ) =  DTDay * ((PhytoAmmoniaPreferenceFactor                      &
                                               * Me%AlfaPhytoNC                                    &
                                               * Me%PhytoGrossGrowRate)-                           &
                                               (PhytoExcretionNitrogen                             &
                                               * Me%PhytoSolublInorgExcreFraction))                   

            Me%Matrix(NA, Phyto   ) =  DTDay * (1. - PhytoAmmoniaPreferenceFactor)                 &
                                             * Me%AlfaPhytoNC                                      &
                                             * Me%PhytoGrossGrowRate 


            Me%Matrix(DONnr, Phyto) = -DTDay * PhytoExcretionNitrogen                              &
                                             * (1.0 - Me%PhytoSolublInorgExcreFraction)            &
                                             * Me%PhytoExcreDissOrgFraction


            Me%Matrix(PON, Phyto  ) = -DTDay * ((PhytoExcretionNitrogen                            & 
                                             * (1.0 - Me%PhytoSolublInorgExcreFraction)            &
                                             * (1.0 - Me%PhytoExcreDissOrgFraction))               &
                                             + DeadPhytoNitrogen)

            
        end if cd5
     

        if (Me%PropCalc%Oxygen) then
        
            Me%Matrix(O, Phyto    ) =  DTDay * (OxygenPhytoRespRate                                &
                                                - PhotoOxygen                                      &
                                                - NitrateOxygen                                    &
                                                - IPOxygen)
        end if            


        !Independent term
         Me%IndTerm(Phyto) = Me%ExternalVar%Mass(Phyto,index)

    !------------------------------------------------------------------------

    end subroutine WQPhytoplankton

    !----------------------------------------------------------------------------
 
    subroutine WQDiatoms(index) 

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        
        integer :: Diatoms
        integer :: BioSi
        integer :: DissSi
        integer :: IP
        integer :: POP, DOPnr
        integer :: AM
        integer :: PON
        integer :: O
        integer :: NA
        integer :: DONnr
        integer :: BOD

        real    :: DTDay

        real    :: DiaExcretionsRate              = null_real
        real    :: DiaExcretionPhosphorus         = null_real
        real    :: DiaExcretionNitrogen           = null_real
        real    :: DeadDiaNitrogen                = null_real
        real    :: DeadDiaPhosphorus              = null_real
        real    :: DiaExcretionSilica             = null_real
        real    :: DeadDiaSilica                  = null_real
        real    :: DiaTotalRespirationLossesRate  = null_real  
        real    :: DiaNitrateOxygen               = null_real
        real    :: DiaIPOxygen                    = null_real 
        
        real    :: DiaAmmoniaPreferenceFactor     = null_real
        real    :: DiaEndogenousRepiration        = null_real
        real    :: DiaPhotorespirationRate        = null_real
        real    :: DiaPhotoOxygen                 = null_real
        real    :: OxygenDiaRespRate              = null_real
        real    :: DiaNonGrazingMortalityRate     = null_real


        real    :: s1, s2, xa, xb, ya, yb
        real    :: x1, x2, x3 ,x4

    !------------------------------------------------------------------------

        POP     = Me%PropIndex%PartOrganicPhosphorus
        DOPnr   = Me%PropIndex%DOPNonRefractory
        IP      = Me%PropIndex%InorganicPhosphorus
        AM      = Me%PropIndex%Ammonia
        NA      = Me%PropIndex%Nitrate
        PON     = Me%PropIndex%PartOrganicNitrogen
        O       = Me%PropIndex%Oxygen
        DONnr   = Me%PropIndex%DONNonRefractory
        BOD     = Me%PropIndex%BOD
        Diatoms = Me%PropIndex%Diatoms
        BioSi   = Me%PropIndex%BiogenicSilica
        DissSi  = Me%PropIndex%DissolvedSilica
        DTDay   = Me%DTDay


        s1       = null_real
        s2       = null_real
        xa       = null_real
        xb       = null_real
        ya       = null_real
        yb       = null_real

        x1       = null_real
        x2       = null_real
        x3       = null_real
        x4       = null_real

    !------------------------------------------------------------------------

    !DiaGrossGrowRate, Diatoms gross grow rate
    
        !DiaNutrientsLimitationFactor, the nutrients limitation factor is the minimum between 
        !Nitrogen, Phosphorus and Silica
cd1 :   if (Me%PropCalc%Nitrogen) then
            
            Me%Diatoms%DiaNLimitationFactor   = (Me%ExternalVar%Mass(AM, index)                        &
                                             + Me%ExternalVar%Mass(NA, index))                     &
                                             / (Me%Diatoms%DiaNSatConst                            &
                                             + Me%ExternalVar%Mass(AM, index)                      &
                                             + Me%ExternalVar%Mass(NA, index))
        else cd1
            Me%Diatoms%DiaNLimitationFactor   = 1.0   !Nitrogen is not a limiting nutrient        
        end if cd1

        !Phosphorus limitation factor
cd2 :   if (Me%PropCalc%Phosphorus.AND.(.NOT.Me%PropCalc%Bacteria)) then
            Me%Diatoms%DiaPLimitationFactor =  Me%ExternalVar%Mass(IP, index)             &
                                            / (Me%Diatoms%DiaPSatConst                             &
                                            + Me%ExternalVar%Mass(IP, index))
        else cd2
            Me%Diatoms%DiaPLimitationFactor = 1.0   !Phosphorus is not a limiting nutrient
        end if cd2

        !Silica limitation factor
cd3 :   if (Me%PropCalc%Silica) then
            Me%Diatoms%DiaSiLimitationFactor =  Me%ExternalVar%Mass(DissSi, index)             &
                                        / (Me%Diatoms%DiaSiSatConst                                &
                                        + Me%ExternalVar%Mass(DissSi, index))
        else cd3
            Me%Diatoms%DiaSiLimitationFactor = 1.0   !Silica is not a limiting nutrient
        end if cd3

       
        Me%Diatoms%DiaNutrientsLimitationFactor = min(Me%Diatoms%DiaNLimitationFactor,     &
                                                      Me%Diatoms%DiaPLimitationFactor,   &        
                                                      Me%Diatoms%DiaSiLimitationFactor     )


    !DiaTempLimitationFactor, temperature effect on Diatom assimilation rate
        s1 = (1.0 / (Me%Diatoms%DiaTOptMin - Me%Diatoms%DiaTMin))                                  &
           * log((Me%Diatoms%DiaK2 * (1.0 - Me%Diatoms%DiaK1))                                     &
           / (Me%Diatoms%DiaK1 * (1.0 - Me%Diatoms%DiaK2)))

        s2 = (1.0 / (Me%Diatoms%DiaTMax - Me%Diatoms%DiaTOptMax))                                  &
           * log((Me%Diatoms%DiaK3 * (1.0 - Me%Diatoms%DiaK4))                                     &
           / (Me%Diatoms%DiaK4 * (1.0 - Me%Diatoms%DiaK3)))

        ya = exp(s1 * (Me%ExternalVar%Temperature(index) - Me%Diatoms%DiaTMin))
        yb = exp(s2 * (Me%Diatoms%DiaTMax - Me%ExternalVar%Temperature(index)))

        xa = (Me%Diatoms%DiaK1 * ya) / (1.0 + Me%Diatoms%DiaK1 * (ya - 1.0))
        xb = (Me%Diatoms%DiaK4 * yb) / (1.0 + Me%Diatoms%DiaK4 * (yb - 1.0))

        Me%Diatoms%DiaTempLimitationFactor = xa * xb


    !DiaLightLimitationFactor, temperature effect on Diatom assimilation rate
        Me%Diatoms%DiaLightLimitationFactor =                                                       &
          PhytoLightLimitationFactor(Thickness      = Me%ExternalVar%Thickness(index),              &
                                     TopRadiation   = Me%ExternalVar%ShortWaveRadiation(index),     &
                                     PExt           = Me%ExternalVar%LightExtCoefField(index),      &
                                     Photoinhibition= Me%Diatoms%DiaPhotoinhibition)


    !DiaGrossGrowRate
        Me%Diatoms%DiaGrossGrowRate = Me%Diatoms%DiaGrowMaxRate                                     &   
                                           * Me%Diatoms%DiaTempLimitationFactor                     &
                                           * Me%Diatoms%DiaLightLimitationFactor                    &
                                           * Me%Diatoms%DiaNutrientsLimitationFactor



cd4 :   if( Me%Diatoms%DiaGrossGrowRate .LT. 0.0) then
            stop 'Subroutine WQDiatoms - Module ModuleWaterQuality - ERR01.'
         
         !No sense on this. David
!        else if (Me%Diatoms%DiaGrossGrowRate .EQ. 0.0) then cd4
!
!            Me%Diatoms%DiaGrossGrowRate = -1.0 / null_real   !Avoid division by zero below
       
        end if cd4

cd5 :   if (MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen).eq.Me%MinOxygen) then
             
            Me%Diatoms%DiaGrossGrowRate = -1.0 / null_real
            
        endif cd5


        !DiatomAmmoniaPreferenceFactor
cd6 :   if (Me%PropCalc%Nitrogen) then
          
            x1 = Me%ExternalVar%Mass(AM, index) * Me%ExternalVar%Mass(NA, index)

            x2 = (Me%Diatoms%DiaNSatConst + Me%ExternalVar%Mass(AM, index))                        &
               * (Me%Diatoms%DiaNSatConst + Me%ExternalVar%Mass(NA, index))

            x3 = Me%Diatoms%DiaNSatConst* Me%ExternalVar%Mass(AM, index)

            x4 = (Me%ExternalVar%Mass(AM, index) + Me%ExternalVar%Mass(NA, index))                 &
               * (Me%Diatoms%DiaNSatConst + Me%ExternalVar%Mass(NA, index))

cd61 :      if ((x1 .EQ. 0.0) .AND. (x3 .EQ. 0.0)) then
                DiaAmmoniaPreferenceFactor = 0.0                 
            else cd61
                DiaAmmoniaPreferenceFactor = (x1 / x2) + (x3 / x4)
            end if cd61
        end if cd6


    !Respiration losses----------------------------------------
       
        !DiaEndogenousRepiration, endogenous respiration
        DiaEndogenousRepiration = Me%Diatoms%DiaEndogRepConst *                                     &
                                  exp(0.069 * Me%ExternalVar%Temperature(index))
      
        if (MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen).eq.Me%MinOxygen) then
            DiaEndogenousRepiration = -1.0 / null_real
        endif
        
        !DiaPhotorespirationRate
        DiaPhotorespirationRate = Me%Diatoms%DiaPhotorespFactor * Me%Diatoms%DiaGrossGrowRate


        !DiaTotalRespirationLosses, total respiration losses
        DiaTotalRespirationLossesRate = DiaEndogenousRepiration + DiaPhotorespirationRate

    !Dead losses----------------------------------------
        !DiaNonGrazingMortalityRate, Dia mortality from non grazing------------
        if (Me%Diatoms%DiaGrossGrowRate .EQ. 0.0) then
            DiaNonGrazingMortalityRate =  Me%Diatoms%DiaMortMaxRate
        
        else
            DiaNonGrazingMortalityRate =  Me%Diatoms%DiaMortMaxRate *                                   &
                                        (( Me%ExternalVar%Mass(Diatoms, index)                          &
                                         /  Me%Diatoms%DiaGrossGrowRate)                                &
                                         / (Me%Diatoms%DiaMortSatConst                                  &
                                         + (Me%ExternalVar%Mass(Diatoms, index)                         &  
                                         /  Me%Diatoms%DiaGrossGrowRate)))
        endif
        
    !Excretion and Dead losses----------------------------------------
        
        !Dias Excretions
        DiaExcretionsRate = Me%Diatoms%DiaExcretionConstant                                         &  
                          * Me%Diatoms%DiaGrossGrowRate                                             &
                          * (1.0 - Me%Diatoms%DiaLightLimitationFactor)
       
        !Silica
            DiaExcretionSilica =  Me%Diatoms%DiaAlfaSiC                                             &
                               * (DiaTotalRespirationLossesRate                                     &
                               + DiaExcretionsRate)

            !DiaNonGrazingMortalityRate, Diatoms Natural Mortality
            DeadDiaSilica = DiaNonGrazingMortalityRate * Me%Diatoms%DiaAlfaSiC

        !Phosphorus        
cd9 :   if (Me%PropCalc%Phosphorus) then

            DiaExcretionPhosphorus = Me%Diatoms%DiaAlfaPC                                           &
                                   * (DiaTotalRespirationLossesRate                                 &
                                   + DiaExcretionsRate)
      
     
            !DeadDiatom, Organic Phosphorus from dead Diatoms
            DeadDiaPhosphorus = DiaNonGrazingMortalityRate * Me%Diatoms%DiaAlfaPC
        
        end if cd9


        !Nitrogen 
cd100 : if (Me%PropCalc%Nitrogen) then

            DiaExcretionNitrogen =  Me%Diatoms%DiaAlfaNC                                            &
                                 * (DiaTotalRespirationLossesRate                                   &
                                 + DiaExcretionsRate)

            !DeadDiatom, Organic Nitrogen from dead Diatom
            DeadDiaNitrogen = DiaNonGrazingMortalityRate * Me%Diatoms%DiaAlfaNC 
        
        end if cd100

    !Oxygen Cycle----------------------------------------

cd7 :   if (Me%PropCalc%Oxygen) then

cd71 :      if (MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen).eq.Me%MinOxygen) then 
               DiaNitrateOxygen  = -1.0 / null_real
               DiaIPOxygen       = -1.0 / null_real
               DiaPhotoOxygen    = -1.0 / null_real
               OxygenDiaRespRate = -1.0 / null_real
           
            else cd71

cd711 :         if (Me%PropCalc%Nitrogen) then
      
                    !Secondary Oxygen production due to Nitrate consumption
                    DiaNitrateOxygen = Me%Diatoms%DiaGrossGrowRate                                  &
                                     * Me%NConsOxyNitRatio                                          &
                                     * Me%Diatoms%DiaAlfaNC * (1 - DiaAmmoniaPreferenceFactor)      
                else cd711
                    DiaNitrateOxygen = 0.0
                end if cd711               

cd712:          if (Me%PropCalc%Phosphorus) then 

                    !Secondary Oxygen production due to Inorganic phosphorus (phosphate) consumption
                    DiaIPOxygen = Me%Diatoms%DiaGrossGrowRate                                       &
                                * Me%PConsOxyPhosphorusRatio                                        &
                                * Me%Diatoms%DiaAlfaPC      
                else cd712
                    DiaIPOxygen = 0.0           
                end if cd712               


                !Oxygen photosynthetic production
                DiaPhotoOxygen = Me%Diatoms%DiaGrossGrowRate * Me%PhotosynthesisOxygenCarbonRatio
    

                !Oxygen loss by respiration 
                OxygenDiaRespRate = Me%PlanktonOxygenCarbonRatio * DiaTotalRespirationLossesRate
            
            end if cd71      
        end if cd7


    !Calculation of system coeficients---------------------------------------
     
            Me%Matrix(Diatoms, Diatoms) = 1.0 + DTDay * ( DiaTotalRespirationLossesRate             &
                                                        + DiaExcretionsRate                         &
                                                        + DiaNonGrazingMortalityRate                &
                                                        - Me%Diatoms%DiaGrossGrowRate)
          
cd101 : if (Me%PropCalc%Phosphorus) then
            
            Me%Matrix(POP, Diatoms)   = - DTDay * ((DiaExcretionPhosphorus                          & 
                                                * (1.0 - Me%Diatoms%DiaSolublInorgExcreFraction)    &
                                                * (1.0 - Me%Diatoms%DiaExcreDissOrgFraction))       &
                                                + DeadDiaPhosphorus)

            Me%Matrix(DOPnr, Diatoms) = - DTDay * DiaExcretionPhosphorus                            &
                                                * (1.0 - Me%Diatoms%DiaSolublInorgExcreFraction)    &
                                                * Me%Diatoms%DiaExcreDissOrgFraction
                          
            Me%Matrix(IP, Diatoms)    =   DTDay * (Me%Diatoms%DiaAlfaPC * Me%Diatoms%DiaGrossGrowRate  &
                                                - (DiaExcretionPhosphorus * Me%Diatoms%DiaSolublInorgExcreFraction))
                    
        end if cd101

cd102 : if (Me%PropCalc%Nitrogen) then

            Me%Matrix(AM, Diatoms)    =  DTDay * ((DiaAmmoniaPreferenceFactor                       &
                                               * Me%Diatoms%DiaAlfaNC                               &
                                               * Me%Diatoms%DiaGrossGrowRate)                       &
                                               - (DiaExcretionNitrogen                              &    
                                               * Me%Diatoms%DiaSolublInorgExcreFraction))             

            Me%Matrix(PON, Diatoms)  = - DTDay * ((DiaExcretionNitrogen                             & 
                                               * (1.0 - Me%Diatoms%DiaSolublInorgExcreFraction)     &
                                               * (1.0 - Me%Diatoms%DiaExcreDissOrgFraction))        &
                                               + DeadDiaNitrogen)

            Me%Matrix(DONnr, Diatoms)= - DTDay * DiaExcretionNitrogen                               &
                                               * (1.0 - Me%Diatoms%DiaSolublInorgExcreFraction)     &
                                               * Me%Diatoms%DiaExcreDissOrgFraction

            Me%Matrix(NA, Diatoms)   =   DTDay * (1. - DiaAmmoniaPreferenceFactor)                  &
                                               * Me%Diatoms%DiaAlfaNC * Me%Diatoms%DiaGrossGrowRate  

        end if cd102     

        if (Me%PropCalc%Oxygen) then
            
            Me%Matrix(O, Diatoms)    =   DTDay * (OxygenDiaRespRate                                 &
                                                  - DiaPhotoOxygen                                  &
                                                  - DiaNitrateOxygen                                &
                                                  - DiaIPOxygen)            
        end if

cd103 : if (Me%PropCalc%Silica) then

            Me%Matrix(DissSi, Diatoms) =   DTDay * Me%Diatoms%DiaAlfaSiC * Me%Diatoms%DiaGrossGrowRate                             
                      
            Me%Matrix(BioSi, Diatoms)  = - DTDay * (DiaExcretionSilica + DeadDiaSilica)                                     

        end if cd103

        !Independent term
         Me%IndTerm(Diatoms) = Me%ExternalVar%Mass(Diatoms,index)

    end subroutine WQDiatoms

    !----------------------------------------------------------------------------

    subroutine WQSilica(index)

        !Arguments-------------------------------------------------------------------
        integer, intent(IN) :: index

        !Local-------------------------------------------------------------------
        real :: SiBiogenicDissRate       = null_real              
    !------------------------------------------------------------------------

        call WQDissolvedSilica  (index,SiBiogenicDissRate)
                                                        
        call WQBiogenicSilica   (index,SiBiogenicDissRate)
    
    end subroutine WQSilica
    
    !------------------------------------------------------------------------
    
    subroutine WQDissolvedSilica(index,SiBiogenicDissRate)

        !Arguments---------------------------------------------------------------
        integer, intent(IN ):: index
        real, intent(OUT)   :: SiBiogenicDissRate

        !Local-------------------------------------------------------------------
        integer             :: DissSi
        integer             :: BioSi
        integer             :: Diatoms
        real                :: DTDay

        !------------------------------------------------------------------------
        DissSi    = Me%PropIndex%DissolvedSilica
        BioSi     = Me%PropIndex%BiogenicSilica
        Diatoms   = Me%PropIndex%Diatoms
        DTDay     = Me%DTDay

       !SiBiogenicDissRate, Biogenic Silica dissolution (1/T)
        SiBiogenicDissRate = Me%SilicaCycle%KSiBiogenicDissRate                                     &
                           * Me%SilicaCycle%BiogenicDissTCoef                                       &
                           **(Me%ExternalVar%Temperature(index) - 20.0) 
                                                       
        !Calculation of system coeficients---------------------------------------
        Me%Matrix(DissSi, BioSi  ) = - SiBiogenicDissRate * DTDay
        Me%Matrix(DissSi, DissSi ) = 1.0

        !Independent term
        Me%IndTerm(DissSi) = Me%ExternalVar%Mass(DissSi, index)     
   
    end subroutine WQDissolvedSilica
   
    !------------------------------------------------------------------------
 
     subroutine WQBiogenicSilica(index,SiBiogenicDissRate)

        !Arguments---------------------------------------------------------------
        integer, intent(IN ):: index
        real, intent(IN)    :: SiBiogenicDissRate


        !Local-------------------------------------------------------------------
        integer             :: BioSi
        real                :: DTDay

        !------------------------------------------------------------------------
        BioSi               = Me%PropIndex%BiogenicSilica
        DTDay               = Me%DTDay

                                                        
        !Calculation of system coeficients---------------------------------------
        Me%Matrix(BioSi, BioSi) = 1.0 + SiBiogenicDissRate * DTDay

        !Independent term
        Me%IndTerm(BioSi) = Me%ExternalVar%Mass(BioSi, index)

    end subroutine WQBiogenicSilica

    !------------------------------------------------------------------------

    subroutine WQNitrogen(index)

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        
        real :: RefrAmmoniaMinRate      = null_real    !RefractoryAmmoniaMineralizationRate
!       real :: NonRefrAmmoniaMinRate   = null_real    !NonRefractoryAmmoniaMineralizationRate
        real :: PartDecompRate          = null_real
        real :: NitrificationRateK1     = null_real 
        real :: NitrificationRateK2     = null_real             

    !------------------------------------------------------------------------


        call WQAmmonia        (index, RefrAmmoniaMinRate,                                           &
!                                     NonRefrAmmoniaMinRate,                                        &
                                      PartDecompRate,                                               &
                                      NitrificationRateK1,                                          &
                                      NitrificationRateK2)

        call WQNitrite        (index, NitrificationRateK1, NitrificationRateK2)
                                                                                                            
        call WQNitrate        (index, NitrificationRateK2)

        call WQOrganicNitrogen(index, RefrAmmoniaMinRate,                                           &
!                                     NonRefrAmmoniaMinRate,                                        &
                                      PartDecompRate)

    !------------------------------------------------------------------------
    end subroutine WQNitrogen


    !----------------------------------------------------------------------------
    !AMMONIA 
    !
    !SOURCES: - excretion by phytoplankton, diatoms and zooplankton;
    !         - decomposition of particulate and dissolved organic matter.
    !
    !SINKS:   - uptake by phytoplankton and diatoms;
    !         - nitrification.

subroutine WQAmmonia(index, RefrAmmoniaMinRate,                                     &
!                           NonRefrAmmoniaMinRate,                                  &
                            PartDecompRate,                                         &
                            NitrificationRateK1,                                    &
                            NitrificationRateK2)
!The mineralization processes assume Phytoplankton N/C Ratio; 


    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

        real, intent(OUT) :: RefrAmmoniaMinRate
!       real, intent(OUT) :: NonRefrAmmoniaMinRate
        real, intent(OUT) :: PartDecompRate                         
        real, intent(OUT) :: NitrificationRateK1                         
        real, intent(OUT) :: NitrificationRateK2

    !Local-------------------------------------------------------------------

        integer :: AM
        integer :: DONnr
        integer :: DONr
        integer :: PON
        integer :: Phyto
        integer :: Diatoms
        integer :: O

        real    :: DTDay

        real    ::  x5                      = null_real
        real    ::  x6                      = null_real
        real    ::  x7                      = null_real 

        real    :: OxygenSinkNitrificationRate     = null_real

    !------------------------------------------------------------------------

        Phyto    = Me%PropIndex%Phyto
        Diatoms  = Me%PropIndex%Diatoms
        AM       = Me%PropIndex%Ammonia
        DONnr    = Me%PropIndex%DONNonRefractory
        DONr     = Me%PropIndex%DissOrganicNitrogenRefractory
        PON      = Me%PropIndex%PartOrganicNitrogen
        O        = Me%PropIndex%Oxygen

        DTDay    = Me%DTDay


       
    !NitrificationRate-------------------------------------------------------
        
        x5 = MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)                                        &
           /(Me%NitrificationSatConst + Me%ExternalVar%Mass(O, index))

        NitrificationRateK1 = Me%KNitrificationRateK1 * Me%TNitrification                               &
                         **(Me%ExternalVar%Temperature(index) - 20.0) * x5

    
        NitrificationRateK2 = Me%KNitrificationRateK2 * Me%TNitrification                           &
                         **(Me%ExternalVar%Temperature(index) - 20.0) * x5
                         
                                            
    !MineralizationRate-------------------------------------------------------

cd2 :   if (Me%PropCalc%Phyto.and.(.not.Me%PropCalc%Diatoms)) then 

cd3 :       if (.NOT.Me%PropCalc%Bacteria) then
            
                !DONnr MineralizationRate:NonRefrAmmoniaMinRate
                x6 = Me%ExternalVar%Mass(Phyto, index)                                              & 
                    / (Me%PhytoNutRegenerationSatConst + Me%ExternalVar%Mass(Phyto, index))

                Me%NonRefrAmmoniaMinRate = Me%KNonRefrAmmoniaMinRate                                &
                                         * Me%TNonRefrAmmoniaMin                                    &
                                         **(Me%ExternalVar%Temperature(index) - 20.0) * x6
            else cd3
               Me%NonRefrAmmoniaMinRate = 0.0
            end if cd3
            
            !DONre MineralizationRate: RefrAmmoniaMinRate
            x7 = Me%ExternalVar%Mass(Phyto, index) / (Me%PhytoNutRegenerationSatConst               &
               + Me%ExternalVar%Mass(Phyto, index))
            
            RefrAmmoniaMinRate = Me%KRefrAmmoniaMinRate * Me%TRefrAmmoniaMin                        &
                               **(Me%ExternalVar%Temperature(index) - 20.0) * x7

        elseif (Me%PropCalc%Phyto.and.Me%PropCalc%Diatoms) then cd2

            if (.NOT.Me%PropCalc%Bacteria) then

                !DONnr MineralizationRate:NonRefrAmmoniaMinRate
                x6 = (Me%ExternalVar%Mass(Phyto, index)+ Me%ExternalVar%Mass(Diatoms, index))       & 
                / (Me%PhytoNutRegenerationSatConst                                                  &
                + (Me%ExternalVar%Mass(Phyto, index)+Me%ExternalVar%Mass(Diatoms, index)))

                Me%NonRefrAmmoniaMinRate = Me%KNonRefrAmmoniaMinRate                                &
                                      * Me%TNonRefrAmmoniaMin                                       &
                                      **(Me%ExternalVar%Temperature(index) - 20.0) * x6
            else 
                Me%NonRefrAmmoniaMinRate = 0.0
            end if 

            !DONre MineralizationRate: RefrAmmoniaMinRate
            x7 = (Me%ExternalVar%Mass(Phyto, index)+ Me%ExternalVar%Mass(Diatoms, index))           &
                 / (Me%PhytoNutRegenerationSatConst                                                 &
                 + (Me%ExternalVar%Mass(Phyto, index)+Me%ExternalVar%Mass(Diatoms, index)))

            RefrAmmoniaMinRate = Me%KRefrAmmoniaMinRate                                             &
                             * Me%TRefrAmmoniaMin                                                   &
                             **(Me%ExternalVar%Temperature(index) - 20.0) * x7

        elseif ((.NOT.Me%PropCalc%Phyto).and.(Me%PropCalc%Diatoms)) then  cd2

            if (.NOT.Me%PropCalc%Bacteria) then

            !DONnr MineralizationRate:NonRefrAmmoniaMinRate
                x6 = Me%ExternalVar%Mass(Diatoms, index)                                            & 
                    / (Me%PhytoNutRegenerationSatConst + Me%ExternalVar%Mass(Diatoms, index))

                Me%NonRefrAmmoniaMinRate = Me%KNonRefrAmmoniaMinRate                                &
                                         * Me%TNonRefrAmmoniaMin                                     &
                                         **(Me%ExternalVar%Temperature(index) - 20.0) * x6
            else 
                Me%NonRefrAmmoniaMinRate = 0.0
            end if 

            !DONre MineralizationRate: RefrAmmoniaMinRate
            x7 = Me%ExternalVar%Mass(Diatoms, index) / (Me%PhytoNutRegenerationSatConst             &
                 + Me%ExternalVar%Mass(Diatoms, index))

            RefrAmmoniaMinRate = Me%KRefrAmmoniaMinRate * Me%TRefrAmmoniaMin                        &
                             **(Me%ExternalVar%Temperature(index) - 20.0) * x7  
        
        else cd2

          Me%NonRefrAmmoniaMinRate = 0.0
          RefrAmmoniaMinRate    = 0.0
 
        end if cd2
       
       !PON MineralizationRate: PartDecompRate (1/T)
       !aqui_5 iste termo não deveria ser só quando não existem bacterias?sim
       PartDecompRate  = Me%KPartDecompRate * Me%TPartDecomposition                                  &
                      **(Me%ExternalVar%Temperature(index) - 20.0) 


    !Oxygen Cycle
    !Oxygen consumption due to nitrification: Ammonia-> Nitrate
    ! NH4 + 1.5O2 -> NO2- +H2O + 2H+
    ! NO2 + 0.5O2 -> NO3
    ! 64/14 O/N

    if (Me%PropCalc%Oxygen)  then  
        
        OxygenSinkNitrificationRate =  NitrificationRateK1 * Me%NitrificationK1_ON_ratio                   
        !MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)
    end if

   !Calculation of system coeficients---------------------------------------

       Me%Matrix(AM, AM   ) = 1.0 + NitrificationRateK1 * DTDay

       if (.NOT.Me%PropCalc%Bacteria) then
       
           Me%Matrix(AM, DONnr) = - Me%NonRefrAmmoniaMinRate * DTDay
 
           Me%Matrix(AM, PON  ) = - PartDecompRate * Me%PhytoAvaibleDecomp * DTDay
       
       end if
        
       Me%Matrix(AM, DONr ) = - RefrAmmoniaMinRate    * DTDay
       

       if (Me%PropCalc%Oxygen) Me%Matrix(O, AM) = DTDay * OxygenSinkNitrificationRate 


    !Independent term
        Me%IndTerm(AM) = Me%ExternalVar%Mass(AM, index)

    !------------------------------------------------------------------------

end subroutine WQAmmonia

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    !NITRITE 
    !
    !SOURCES: - nitrification step 1 (N-NH4-->N-NO2).
    !
    !SINKS:   - nitrification step 2 (N-NO2-->N-NO3).

subroutine WQNitrite(index, NitrificationRateK1, NitrificationRateK2)

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

        real, intent(IN) :: NitrificationRateK1, NitrificationRateK2
                                 

    !Local-------------------------------------------------------------------

        integer :: AM
        integer :: NI
        integer :: O
        
        real :: DTDay
        
        real    :: OxygenSinkNitrificationRateK2     = null_real

    !------------------------------------------------------------------------

        AM = Me%PropIndex%Ammonia
        NI = Me%PropIndex%Nitrite
        O  = Me%PropIndex%Oxygen
        
        DTDay   = Me%DTDay
        
    if (Me%PropCalc%Oxygen)  then  
        
        OxygenSinkNitrificationRateK2 =  NitrificationRateK2 * Me%NitrificationK2_ON_ratio                   
        !MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)
    end if
    

    !Calculation of system coeficients---------------------------------------
        Me%Matrix(NI, NI) =   DTDay * (NitrificationRateK2) + 1.0 
        
        Me%Matrix(NI, AM) =   DTDay * (-NitrificationRateK1)
        
        if (Me%PropCalc%Oxygen) Me%Matrix(O, NI) = DTDay * OxygenSinkNitrificationRateK2 

    !Independent term
        Me%IndTerm(NI) = Me%ExternalVar%Mass(NI, index) 


    !------------------------------------------------------------------------

end subroutine WQNitrite

    !----------------------------------------------------------------------------







    !----------------------------------------------------------------------------
    !NITRATE 
    !
    !SOURCES: - nitrification step 2 (N-NO2-->N-NO3).
    !
    !SINKS:   - denitrification (N-NO3-->N2) under anaerobic conditions;
    !         - phytoplankton uptake.

subroutine WQNitrate(index, NitrificationRateK2)

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

        real, intent(IN) :: NitrificationRateK2                       

    !Local-------------------------------------------------------------------

        integer :: O
        integer :: NA
        integer :: NI
        integer :: BOD

        real    :: ODsourceDenitrificationRate   = null_real
        real    :: DenitrificationRate           = null_real
        real    :: DTDay
        real    :: x1                            = null_real

    !------------------------------------------------------------------------

        O       = Me%PropIndex%Oxygen
        NA      = Me%PropIndex%Nitrate
        NI      = Me%PropIndex%Nitrite
        BOD     = Me%PropIndex%BOD

        DTDay     = Me%DTDay

    !------------------------------------------------------------------------

                         
    !DenitrificationRate, denitrification rate
        x1 = Me%DenitrificationSatConst / (Me%DenitrificationSatConst                               &
                                         + MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen))

        DenitrificationRate = Me%KDenitrificationRate * Me%TDenitrification                         &
                            ** (Me%ExternalVar%Temperature(index) - 20.) * x1


    !BOD not consumed due to anaerobic organic matter decomposition during Denitrification
        !aqui_6 5/4 não percebi...acho que não existe libertação de O2 com a desnitrificação 
        if (Me%PropCalc%Oxygen)                                                                     &
            ODsourceDenitrificationRate = 5.0 / 4.0 * Me%NConsOxyNitRatio                           &
                                                    * DenitrificationRate

    !Calculation of system coeficients---------------------------------------
        Me%Matrix(NA, NA) = 1. + DTDay * (DenitrificationRate)


        Me%Matrix(NA, NI) = -DTDay * NitrificationRateK2

        if (Me%PropCalc%Oxygen)                                                                     &
            Me%Matrix(O, NA) = -DTDay * ODsourceDenitrificationRate   


    !Independent term
        Me%IndTerm(NA) = Me%ExternalVar%Mass(NA, index) 

    !------------------------------------------------------------------------

end subroutine WQNitrate

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

subroutine WQOrganicNitrogen(index, RefrAmmoniaMinRate,                             &
!                                   NonRefrAmmoniaMinRate,                          &
                                    PartDecompRate)

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

        real,    intent(IN) :: RefrAmmoniaMinRate 
!       real,    intent(IN) :: NonRefrAmmoniaMinRate 
        real,    intent(IN) :: PartDecompRate 

    !------------------------------------------------------------------------

        call WQParticulateOrganicNitrogen (index, PartDecompRate)

        call WQDONRefractory   (index, RefrAmmoniaMinRate, PartDecompRate)

        call WQDONNonRefractory(index)
        
        

    !------------------------------------------------------------------------

end subroutine WQOrganicNitrogen

    !----------------------------------------------------------------------------







    !----------------------------------------------------------------------------

subroutine WQParticulateOrganicNitrogen(index, PartDecompRate)
!The mineralization processes assume Phytoplankton N/C Ratio if phyto and diatoms are considered; 

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

        real,    intent(IN) :: PartDecompRate

    !Local-------------------------------------------------------------------

        integer :: PON, PONr, O

        real    :: DTDay

    !------------------------------------------------------------------------

        PON      = Me%PropIndex%PartOrganicNitrogen
        PONr     = Me%PropIndex%PartOrganicNitrogenRefractory
        O        = Me%PropIndex%Oxygen
        DTDay    = Me%DTDay
   

    !Calculation of system coeficients---------------------------------------

        if (Me%PropCalc%Bacteria) then

            Me%Matrix(PON, PON) = 1.0  
            Me%Matrix(PONr, PONr) = 1.0  

        else

            Me%Matrix(PON, PON) = DTDay * PartDecompRate + 1.0 

         !aqui_7
         !CH2O+O2 -> CO2 + H2O
         if (Me%PropCalc%Oxygen) then
         
            if(.NOT. Me%PropCalc%BOD) then
            
                Me%Matrix(O, PON) = DTDay * PartDecompRate * 1/Me%OMAlfaNC * Me%OxyCarbonRatio * &
                MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)                         &
                /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) +0.5)
            endif
         endif
        endif
         !aqui_7

    !Independent term
        Me%IndTerm(PON) = Me%ExternalVar%Mass(PON, index) 


        if (Me%PropCalc%Bacteria) then

            Me%IndTerm(PONr) = Me%ExternalVar%Mass(PONr, index) 

       end if
    !------------------------------------------------------------------------

end subroutine WQParticulateOrganicNitrogen

    !----------------------------------------------------------------------------






    !----------------------------------------------------------------------------

subroutine WQDONNonRefractory(index)
 
    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !External----------------------------------------------------------------

!        real, intent(IN) :: NonRefrAmmoniaMinRate

    !Local-------------------------------------------------------------------

        integer :: DONnr, O

        real    :: DTDay

    !------------------------------------------------------------------------

        DONnr = Me%PropIndex%DONNonRefractory
        O     = Me%PropIndex%Oxygen
        DTDay = Me%DTDay



    !Calculation of system coeficients---------------------------------------

cd1 :   if (Me%PropCalc%Bacteria) then
        
            Me%Matrix(DONnr, DONnr) = 1.0  
        
        else
           
            Me%Matrix(DONnr, DONnr) = 1.0  + DTDay * Me%NonRefrAmmoniaMinRate
        
            if (Me%PropCalc%Oxygen)  then
            
                if(.NOT. Me%PropCalc%BOD) then                                                               
              
                    Me%Matrix(O, DONnr) = DTDay * Me%NonRefrAmmoniaMinRate                          &
                                        * 1/Me%OMAlfaNC * Me%OxyCarbonRatio                         &
                                        * MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)           &
                                        / (MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) +0.5)                
                endif
            endif
        endif cd1


    !Independent term
        Me%IndTerm(DONnr) = Me%ExternalVar%Mass(DONnr, index)

    !------------------------------------------------------------------------

end subroutine WQDONNonRefractory

    !----------------------------------------------------------------------------






    !----------------------------------------------------------------------------

subroutine WQDONRefractory(index, RefrAmmoniaMinRate, PartDecompRate)

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

        real, intent(IN) :: RefrAmmoniaMinRate, PartDecompRate

    !Local-------------------------------------------------------------------

        integer :: DONre
        integer :: PON
        integer :: O
        real    :: DTDay

    !------------------------------------------------------------------------

        DONre     = Me%PropIndex%DissOrganicNitrogenRefractory
        PON     = Me%PropIndex%PartOrganicNitrogen
          O     = Me%PropIndex%Oxygen
        DTDay   = Me%DTDay




    !Calculation of system coeficients---------------------------------------
        Me%Matrix(DONre, DONre) = 1.0 + DTDay * RefrAmmoniaMinRate

        Me%Matrix(DONre, PON) =-DTDay * PartDecompRate * (1.0 - Me%PhytoAvaibleDecomp)

        if (Me%PropCalc%Oxygen) then
        
            if(.NOT. Me%PropCalc%BOD) then
            
                Me%Matrix(O, DONre) = DTDay * RefrAmmoniaMinRate                                    &
                                  * 1/Me%OMAlfaNC * Me%OxyCarbonRatio                               &
                                  * MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)                 &
                                  / (MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) +0.5)
            endif
        endif

    !Independent term
        Me%IndTerm(DONre) = Me%ExternalVar%Mass(DONre, index) 

    !------------------------------------------------------------------------

end subroutine WQDONRefractory

    !----------------------------------------------------------------------------



    
    
    
    !--------------------------------------------------------------------------

subroutine WQPhosphorus(index)

    !Arguments-------------------------------------------------------------

        integer, intent(IN) :: index

    !External--------------------------------------------------------------
        real                :: POPDecompRate
        real                :: PhosphorusDOPrMinRate
        real                :: PhosphorusDOPnrMinRate

    !----------------------------------------------------------------------

        call WQInorganicPhosphorus(index, POPDecompRate,                            &
                                   PhosphorusDOPrMinRate, PhosphorusDOPnrMinRate)
        call WQOrganicPhosphorus  (index, POPDecompRate,                            &
                                   PhosphorusDOPrMinRate, PhosphorusDOPnrMinRate)

    !----------------------------------------------------------------------

end subroutine WQPhosphorus

    !--------------------------------------------------------------------------





    !----------------------------------------------------------------------------
    !INORGANIC PHOSPHORUS 
    !
    !SOURCES: - excretion by phytoplankton and zooplankton;
    !         - decomposition of particulate and dissolved organic matter.
    !
    !SINKS:   - uptake by phytoplankton.

    subroutine WQInorganicPhosphorus(index, POPDecompRate,                          &
                                     PhosphorusDOPrMinRate, PhosphorusDOPnrMinRate)

    !Arguments---------------------------------------------------------------
        real, intent(OUT)                           :: POPDecompRate
        real, intent(OUT)                           :: PhosphorusDOPrMinRate
        real, intent(OUT)                           :: PhosphorusDOPnrMinRate
        integer, intent(IN)                         :: index

        

    !Local-------------------------------------------------------------------

        integer :: Zoo
        integer :: Phyto
        integer :: Diatoms
        integer :: IP
        integer :: POP
        integer :: DOPnr
        integer :: DOPr

        real    :: DTDay
        real    :: x1

    !------------------------------------------------------------------------



        Zoo       = Me%PropIndex%Zoo  
        Phyto     = Me%PropIndex%Phyto
        Diatoms   = Me%PropIndex%Diatoms
        IP        = Me%PropIndex%InorganicPhosphorus
        POP       = Me%PropIndex%PartOrganicPhosphorus
        DOPnr     = Me%PropIndex%DOPNonRefractory
        DOPr      = Me%PropIndex%DissOrganicPhosphorusRefractory

        DTDay     = Me%DTDay
        
        !this relation was lacking and unabling implicit LUD computation (division by zero on this diagonal)
        Me%Matrix(IP, IP) = 1.0

    !PhosphorusMineralizationRate, rate of Phosphorus mineralizationRate
  
cd4:    if (Me%PropCalc%Phyto) then
            
            if (Me%PropCalc%Diatoms) then
        
                x1 = (Me%ExternalVar%Mass(Phyto, index)+ Me%ExternalVar%Mass(Diatoms, index))       &
                    /(Me%PhytoNutRegenerationSatConst + (Me%ExternalVar%Mass(Phyto, index)          &
                    +Me%ExternalVar%Mass(Diatoms, index)))
       
            else
        
                x1 = Me%ExternalVar%Mass(Phyto, index)                                              &
                   / (Me%PhytoNutRegenerationSatConst + Me%ExternalVar%Mass(Phyto, index))

            endif

        else 
        
            if (Me%PropCalc%Diatoms) then 
 
                x1 = Me%ExternalVar%Mass(Diatoms, index)                                            &
                   / (Me%PhytoNutRegenerationSatConst + Me%ExternalVar%Mass(Diatoms, index))
            endif
        endif cd4
            
        !Non refractary phosphorus
        PhosphorusDOPnrMinRate       = Me%KDOPnrMinRate                                             &
                                     * Me%TDOPnrMin                                                 &
                                     **(Me%ExternalVar%Temperature(index) - 20.0) * x1
        !Refractary phosphorus !aqui_8        
        PhosphorusDOPrMinRate        = Me%KDOPrMinRate                                              &
                                     * Me%TDOPrMin                                                  &
                                     **(Me%ExternalVar%Temperature(index) - 20.0) * x1

        POPDecompRate = Me%POPDecompRate * Me%TPOPDecompRate                                        &
                    **(Me%ExternalVar%Temperature(index) - 20.0)

    
    !Calculation of system coeficients---------------------------------------

        

        Me%Matrix(IP, POP)   = - DTDay * POPDecompRate * Me%PhytoAvaibleDecomp
        
        Me%Matrix(IP, DOPnr) = - DTDay * PhosphorusDOPnrMinRate

        Me%Matrix(IP, DOPr)  = - DTDay * PhosphorusDOPrMinRate


        !Independent term
        Me%IndTerm(IP) = Me%ExternalVar%Mass(IP, index)
 
 

    !------------------------------------------------------------------------

end subroutine WQInorganicPhosphorus

    !----------------------------------------------------------------------------






    !----------------------------------------------------------------------------
    !ORGANIC PHOSPHORUS 
    !
    !SOURCES: - excretion by phytoplankto, diatoms and zooplankton;
    !         - decomposition of particulate and dissolved organic matter.
    !
    !SINKS:   - mineralization.

    subroutine WQOrganicPhosphorus(index, POPDecompRate,                                            &
                                   PhosphorusDOPrMinRate, PhosphorusDOPnrMinRate)

        !Arguments---------------------------------------------------------------
        integer, intent(IN)                         :: index
        real, intent(IN)                            :: POPDecompRate
        real, intent(IN)                            :: PhosphorusDOPrMinRate
        real, intent(IN)                            :: PhosphorusDOPnrMinRate

        !Local-------------------------------------------------------------------
        integer :: POP
        integer :: DOPnr
        integer :: DOPr
        integer :: IP
        integer :: O
        real    :: DTDay

        !------------------------------------------------------------------------

        IP      = Me%PropIndex%InorganicPhosphorus
        O       = Me%PropIndex%Oxygen
        DTDay   = Me%DTDay
        POP     = Me%PropIndex%PartOrganicPhosphorus
        DOPnr   = Me%PropIndex%DOPNonRefractory
        DOPr    = Me%PropIndex%DissOrganicPhosphorusRefractory


        !Calculation of system coeficients for POP---------------------------------------
        Me%Matrix(POP, POP) = DTDay * POPDecompRate + 1.0 

         if (Me%PropCalc%Oxygen)     then

            if(.NOT. Me%PropCalc%BOD) then
            
                if(.NOT. Me%PropCalc%Nitrogen) then
                
                Me%Matrix(O, POP) = DTDay * POPDecompRate * 1/Me%OMAlfaPC * Me%OxyCarbonRatio *     &
                                   MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)                  &
                                   /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) +0.5)
                endif   
            endif
        endif

        !Calculation of system coeficients for DOPnr---------------------------------------       
        Me%Matrix(DOPnr, DOPnr) = 1.0  + DTDay * PhosphorusDOPnrMinRate

        if (Me%PropCalc%Oxygen)     then
        
            if(.NOT. Me%PropCalc%BOD) then
            
                if(.NOT. Me%PropCalc%Nitrogen) then

                Me%Matrix(O, DOPnr) = DTDay * PhosphorusDOPnrMinRate * 1/Me%OMAlfaPC * Me%OxyCarbonRatio * &
                                   MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)                 &
                                   /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) +0.5)
                endif
            endif
        endif

        !Calculation of system coeficients for DOPr---------------------------------------
        Me%Matrix(DOPr, DOPr) = 1.0 + DTDay * PhosphorusDOPrMinRate

        Me%Matrix(DOPr, POP) =-DTDay * POPDecompRate * (1.0 - Me%PhytoAvaibleDecomp)
       
        if (Me%PropCalc%Oxygen)     then
        
            if(.NOT. Me%PropCalc%BOD) then
            
                if(.NOT. Me%PropCalc%Nitrogen) then

                Me%Matrix(O, DOPr) = DTDay * PhosphorusDOPrMinRate * 1/Me%OMAlfaPC * Me%OxyCarbonRatio *     &
                                   MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)                 &
                                   /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) +0.5)
                endif
            endif
        endif

    !Independent term
        Me%IndTerm(DOPr) = Me%ExternalVar%Mass(DOPr, index) 

    !Independent term
        Me%IndTerm(DOPnr) = Me%ExternalVar%Mass(DOPnr, index)
    
    !Independent term
        Me%IndTerm(POP) = Me%ExternalVar%Mass(POP, index) 

    !------------------------------------------------------------------------

    end subroutine WQOrganicPhosphorus
    
    
!------------------------------------------------------------------------    
    
    
    subroutine WQPompools(index)                             ! by M Mateus @ MARETEC oct_2009

        !Arguments---------------------------------------------------------------
        integer, intent(IN) :: index

        !Local-----------------------------------------------------------------
        integer                             :: AM
        integer                             :: IP
        integer                             :: PON1, PON2, PON3, PON4, PON5
        integer                             :: POP1, POP2, POP3, POP4, POP5
        integer                             :: DOPr, DONre
        integer                             :: Zoo
        integer                             :: O

        real                                :: DTDay
        real                                :: PONIngestPot, PONIngest
        real                                :: POPIngestPot, POPIngest
        

        !External--------------------------------------------------------------
        real                :: POPDecompRate
        real                :: PartDecompRate

        !------------------------------------------------------------------------

        IP      = Me%PropIndex%InorganicPhosphorus
        AM      = Me%PropIndex%Ammonia
        Zoo     = Me%PropIndex%Zoo
        O       = Me%PropIndex%Oxygen
        DOPr    = Me%PropIndex%DissOrganicPhosphorusRefractory
        DONre   = Me%PropIndex%DissOrganicNitrogenRefractory 
        
        DTDay   = Me%DTDay
        
        PONIngestPot       = null_real
        PONIngest          = null_real
        POPIngestPot       = null_real
        POPIngest          = null_real
        
        
        if (Me%PropCalc%Nitrogen) then
        
            PON1    = Me%PropIndex%PONitrogen1
            PON2    = Me%PropIndex%PONitrogen2
            PON3    = Me%PropIndex%PONitrogen3
            PON4    = Me%PropIndex%PONitrogen4
            PON5    = Me%PropIndex%PONitrogen5
        
        end if
        
        if (Me%PropCalc%Phosphorus) then
        
            POP1    = Me%PropIndex%POPhosphorus1
            POP2    = Me%PropIndex%POPhosphorus2
            POP3    = Me%PropIndex%POPhosphorus3
            POP4    = Me%PropIndex%POPhosphorus4
            POP5    = Me%PropIndex%POPhosphorus5
        
        end if
                
        !----------------------------------------------------------------------
    
    !Calculation of system coeficients for POM pools------------------------------
        
        if (Me%PropCalc%Nitrogen) then
        
        PartDecompRate  = Me%KPartDecompRate * Me%TPartDecomposition                                  &
                      **(Me%ExternalVar%Temperature(index) - 20.0)
                   
        PONIngestPot = Me%TZooLimitationFactor * Me%POMIngestVmax *                                   &
                      (Me%ExternalVar%Mass(PON1, index) / (Me%ExternalVar%Mass(PON1, index) + Me%PONIngestKs))
                     
        PONIngest    = MIN(PONIngestPot * Me%ExternalVar%Mass(Zoo, index),                            &
                       Me%ExternalVar%Mass(PON1, index) / DTDay)  
       
       
            Me%Matrix(PON1, PON1) =  DTDay * (PartDecompRate + PONIngest) + 1.0
            Me%Matrix(PON2, PON2) =  DTDay * PartDecompRate + 1.0
            Me%Matrix(PON3, PON3) =  DTDay * PartDecompRate + 1.0
            Me%Matrix(PON4, PON4) =  DTDay * PartDecompRate + 1.0
            Me%Matrix(PON5, PON5) =  DTDay * PartDecompRate + 1.0              
                                    
            
    !Calculation of system coeficients---------------------------------------
            
            Me%Matrix(O, PON1) = DTDay * PartDecompRate * 1/Me%OMAlfaNC * Me%OxyCarbonRatio *     &
                                 MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)                  &
                                 /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) +0.5)
            
            Me%Matrix(O, PON2) = DTDay * PartDecompRate * 1/Me%OMAlfaNC * Me%OxyCarbonRatio *     &
                                 MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)                  &
                                 /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) +0.5)
            
            Me%Matrix(O, PON3) = DTDay * PartDecompRate * 1/Me%OMAlfaNC * Me%OxyCarbonRatio *     &
                                 MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)                  &
                                 /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) +0.5)
            
            Me%Matrix(O, PON4) = DTDay * PartDecompRate * 1/Me%OMAlfaNC * Me%OxyCarbonRatio *     &
                                 MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)                  &
                                 /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) +0.5)
            
            Me%Matrix(O, PON5) = DTDay * PartDecompRate * 1/Me%OMAlfaNC * Me%OxyCarbonRatio *     &
                                 MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen)                  &
                                 /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) +0.5)
            
            Me%Matrix(AM, PON1  ) = - PartDecompRate * Me%PhytoAvaibleDecomp * DTDay
            
            Me%Matrix(AM, PON2  ) = - PartDecompRate * Me%PhytoAvaibleDecomp * DTDay
            
            Me%Matrix(AM, PON3  ) = - PartDecompRate * Me%PhytoAvaibleDecomp * DTDay
            
            Me%Matrix(AM, PON4  ) = - PartDecompRate * Me%PhytoAvaibleDecomp * DTDay
           
            Me%Matrix(AM, PON5  ) = - PartDecompRate * Me%PhytoAvaibleDecomp * DTDay
            
            
            Me%Matrix(DONre, PON1) = -DTDay * PartDecompRate * (1.0 - Me%PhytoAvaibleDecomp)
            
            Me%Matrix(DONre, PON2) = -DTDay * PartDecompRate * (1.0 - Me%PhytoAvaibleDecomp)
            
            Me%Matrix(DONre, PON3) = -DTDay * PartDecompRate * (1.0 - Me%PhytoAvaibleDecomp)
            
            Me%Matrix(DONre, PON4) = -DTDay * PartDecompRate * (1.0 - Me%PhytoAvaibleDecomp)
            
            Me%Matrix(DONre, PON5) = -DTDay * PartDecompRate * (1.0 - Me%PhytoAvaibleDecomp)
            
            
            Me%Matrix(Zoo, PON1) = -DTDay * (PONIngest * Me%PON_CNratio)
            
            
            !Independent term
            Me%IndTerm(PON1) = Me%ExternalVar%Mass(PON1, index)
            Me%IndTerm(PON2) = Me%ExternalVar%Mass(PON2, index)
            Me%IndTerm(PON3) = Me%ExternalVar%Mass(PON3, index)
            Me%IndTerm(PON4) = Me%ExternalVar%Mass(PON4, index)
            Me%IndTerm(PON5) = Me%ExternalVar%Mass(PON5, index)
        
        end if
        
           
        
        
        if (Me%PropCalc%Phosphorus) then
                        
            !Calculation of system coeficients for POM pools------------------------------
                       
            POPDecompRate = Me%POPDecompRate * Me%TPOPDecompRate                                        &
                    **(Me%ExternalVar%Temperature(index) - 20.0)
                    
            POPIngestPot = Me%TZooLimitationFactor * Me%POMIngestVmax *                                   &
                      (Me%ExternalVar%Mass(POP1, index) / (Me%ExternalVar%Mass(POP1, index) + Me%POPIngestKs))
                     
            POPIngest    = MIN(POPIngestPot * Me%ExternalVar%Mass(Zoo, index),                            &
                           Me%ExternalVar%Mass(POP1, index) / DTDay)  
       
         
         
            Me%Matrix(POP1, POP1) =  DTDay * (POPDecompRate + POPIngest) + 1.0
            Me%Matrix(POP2, POP2) =  DTDay * POPDecompRate + 1.0
            Me%Matrix(POP3, POP3) =  DTDay * POPDecompRate + 1.0
            Me%Matrix(POP4, POP4) =  DTDay * POPDecompRate + 1.0
            Me%Matrix(POP5, POP5) =  DTDay * POPDecompRate + 1.0
            
    
    !Calculation of system coeficients---------------------------------------

        Me%Matrix(O, POP1) = DTDay * POPDecompRate * 1/Me%OMAlfaPC * Me%OxyCarbonRatio *     &
                                 MAX(Me%ExternalVar%Mass(O, index), Me%MinOxygen)            &
                                 /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) + 0.5) 
        
        Me%Matrix(O, POP2) = DTDay * POPDecompRate * 1/Me%OMAlfaPC * Me%OxyCarbonRatio *     &
                                 MAX(Me%ExternalVar%Mass(O, index), Me%MinOxygen)            &
                                 /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) + 0.5)
        
        Me%Matrix(O, POP3) = DTDay * POPDecompRate * 1/Me%OMAlfaPC * Me%OxyCarbonRatio *     &
                                 MAX(Me%ExternalVar%Mass(O, index), Me%MinOxygen)            &
                                 /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) + 0.5)
        
        Me%Matrix(O, POP4) = DTDay * POPDecompRate * 1/Me%OMAlfaPC * Me%OxyCarbonRatio *     &
                                 MAX(Me%ExternalVar%Mass(O, index), Me%MinOxygen)            &
                                 /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) + 0.5)
        
        Me%Matrix(O, POP5) = DTDay * POPDecompRate * 1/Me%OMAlfaPC * Me%OxyCarbonRatio *     &
                                 MAX(Me%ExternalVar%Mass(O, index), Me%MinOxygen)            &
                                 /(MAX(Me%ExternalVar%Mass(O, index),Me%MinOxygen) + 0.5)
        
               
        Me%Matrix(IP, POP1)   = - DTDay * POPDecompRate * Me%PhytoAvaibleDecomp
        
        Me%Matrix(IP, POP2)   = - DTDay * POPDecompRate * Me%PhytoAvaibleDecomp
        
        Me%Matrix(IP, POP3)   = - DTDay * POPDecompRate * Me%PhytoAvaibleDecomp
        
        Me%Matrix(IP, POP4)   = - DTDay * POPDecompRate * Me%PhytoAvaibleDecomp
        
        Me%Matrix(IP, POP5)   = - DTDay * POPDecompRate * Me%PhytoAvaibleDecomp
        
        
        Me%Matrix(DOPr, POP1) = - DTDay * POPDecompRate * (1.0 - Me%PhytoAvaibleDecomp)
        
        Me%Matrix(DOPr, POP2) = - DTDay * POPDecompRate * (1.0 - Me%PhytoAvaibleDecomp)
        
        Me%Matrix(DOPr, POP3) = - DTDay * POPDecompRate * (1.0 - Me%PhytoAvaibleDecomp)
        
        Me%Matrix(DOPr, POP4) = - DTDay * POPDecompRate * (1.0 - Me%PhytoAvaibleDecomp)
        
        Me%Matrix(DOPr, POP5) = - DTDay * POPDecompRate * (1.0 - Me%PhytoAvaibleDecomp)
        
        
        Me%Matrix(Zoo, POP1) = -DTDay * (POPIngest * Me%PON_CPratio)
            
            
            !Independent term             
            Me%IndTerm(POP1) = Me%ExternalVar%Mass(POP1, index)
            Me%IndTerm(POP2) = Me%ExternalVar%Mass(POP2, index)
            Me%IndTerm(POP3) = Me%ExternalVar%Mass(POP3, index)
            Me%IndTerm(POP4) = Me%ExternalVar%Mass(POP4, index)
            Me%IndTerm(POP5) = Me%ExternalVar%Mass(POP5, index)
        
        end if
        
    end subroutine WQPompools   
    
    
    
    
    
    
    
    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillWaterQuality(WaterQualityID, STAT)

        !Arguments---------------------------------------------------------------
        integer                        :: WaterQualityID
        integer, optional, intent(OUT) :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_CALL
        integer                                     :: STAT_     
        integer                                     :: nUsers
        !------------------------------------------------------------------------                      

        STAT_ = UNKNOWN_

        call Ready(WaterQualityID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mWATERQUALITY_,  Me%InstanceID)
  
            if (nUsers == 0) then

                if(Me%ObjLUD /= 0) then
                    call KillLUD(Me%ObjLUD, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Subroutine KillWaterQuality - ModuleWaterQuality. ERR01.'
                end if

                deallocate(Me%IndTerm, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'Subroutine Kill_WaterQuality - ModuleWaterQuality. ERR02.'
                nullify(Me%IndTerm)


                deallocate(Me%Matrix, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'Subroutine Kill_WaterQuality - ModuleWaterQuality. ERR03.'
                nullify(Me%Matrix)


cd4 :           if (associated(Me%NewMass)) then
                    deallocate(Me%NewMass, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Subroutine Kill_WaterQuality - ModuleWaterQuality. ERR04.'
                    nullify(Me%NewMass)
                end if cd4


                call DeallocateInstance

                WaterQualityID = 0


                STAT_ = SUCCESS_

            endif

        else 

            STAT_ = ready_


        end if cd1


        if (present(STAT)) STAT = STAT_

    !------------------------------------------------------------------------

    end subroutine KillWaterQuality

    !------------------------------------------------------------------------
    
    
    !------------------------------------------------------------------------

    subroutine DeallocateInstance

        !Local-----------------------------------------------------------------
        type (T_WaterQuality), pointer           :: AuxObjWaterQuality
        type (T_WaterQuality), pointer           :: PreviousObjWaterQuality

        !Updates pointers
        if (Me%InstanceID == FirstObjWaterQuality%InstanceID) then
            FirstObjWaterQuality => FirstObjWaterQuality%Next
        else
            PreviousObjWaterQuality => FirstObjWaterQuality
            AuxObjWaterQuality      => FirstObjWaterQuality%Next
            do while (AuxObjWaterQuality%InstanceID /= Me%InstanceID)
                PreviousObjWaterQuality => AuxObjWaterQuality
                AuxObjWaterQuality      => AuxObjWaterQuality%Next
            enddo

            !Now update linked list
            PreviousObjWaterQuality%Next => AuxObjWaterQuality%Next

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

    subroutine Ready (WaterQualityID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: WaterQualityID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (WaterQualityID > 0) then
            call LocateObjWaterQuality (WaterQualityID)
            ready_ = VerifyReadLock (mWATERQUALITY_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1
     
        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjWaterQuality (WaterQualityID)

        !Arguments-------------------------------------------------------------
        integer                                     :: WaterQualityID

        !Local-----------------------------------------------------------------

        Me => FirstObjWaterQuality
        do while (associated (Me))
            if (Me%InstanceID == WaterQualityID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleWaterQuality - LocateObjWaterQuality - ERR01'

    end subroutine LocateObjWaterQuality

    !--------------------------------------------------------------------------

end module ModuleWaterQuality

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
