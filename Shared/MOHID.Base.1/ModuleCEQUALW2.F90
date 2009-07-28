!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : CEQUALW2
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Pedro Pina, Luis Fernandes - v4.0
! DESCRIPTION   : U.S. Army Corps of Engineers zero-dimensional model for primary production 
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
Module ModuleCEQUALW2
    
    use ModuleGlobalData
    use ModuleFunctions, only: OxygenSaturationCEQUALW2
    use ModuleEnterData

    implicit none

    private 

    !subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartCEQUALW2
    private ::      AllocateInstance
    private ::      ReadWaterColumnData
    private ::      ReadBenthicData
    private ::          ReadBenthicGlobalVariables
    private ::          ReadSODData    
    private ::          ReadDetritusParameters        
    private ::          ReadAmmoniaParameters      
    private ::          ReadPhosphorusParameters 
    private ::          ReadDsilicaParameters          
    private ::          ReadICarbonParameters
    private ::          PropertyIndexNumber
    private ::          BenthicPropertyIndexNumber
    private ::          RateIndexNumber
    private ::          ReadGlobalVariables
    private ::          ConstructAlgaeClasses
    private ::              AddAlgae
    private ::              ReadAlgaeParameters
    private ::          ConstructEpiphytonClasses
    private ::              AddEpiphyton
    private ::              ReadEpiphytonParameters 
    private ::          ReadOMParameters
    private ::          ReadSilicaParameters
    private ::          ReadOxygenParameters
    private ::          ReadNitrogenParameters
    private ::          ConstructPropertyList
    private ::          ConstructBenthicPropertyList
   
    !Selector
    public  :: GetDTCEQUALW2
    public  :: GetCEQUALW2Options
    public  :: GetCEQUALW2Size   
    public  :: GetCEQUALW2PropIndex
    public  :: GetCEQUALW2PropertyList
    public  :: GetCEQUALW2RateFlux
    public  :: UnGetCEQUALW2RateFlux
    public  :: UngetCEQUALW2
    
    
    !Modifier
    
               ! Water Column Processes
    public  :: CEQUALW2            
               
    private ::      TemperatureRateMultipliers    
    private ::      ComputeKineticRates            
    private ::      ComputePhosphorus         
    private ::      ComputeAmmonia            
    private ::      ComputeNitrate            
    private ::      ComputeDissolvedSilica    
    private ::      ComputeParticulateSilica  
    private ::      ComputeRefDOM
    private ::      ComputeLabDOM
    private ::      ComputeLabPOM
    private ::      ComputeRefPOM
    private ::      ComputeAlgae
    private ::      ComputeBOD 
    private ::      ComputeDissolvedOxygen
    private ::      ComputeICarbon
    private ::      ComputeDetritus
    private ::      ComputeEpiphyton
    private ::      ComputepH_CO2
    
               ! Benthic Processes

    public  ::  CEQUALW2Benthic 
       
    private ::      ComputeBenthicKineticRates    
    private ::      ComputeBenthicPhosphorus         
    private ::      ComputeBenthicAmmonia            
    private ::      ComputeBenthicDissolvedSilica    
    private ::      ComputeBenthicDissolvedOxygen
    private ::      ComputeBenthicICarbon
    private ::      ComputeBenthicDetritus

               ! Auxiliar Functions
    private ::      Rising
    private ::      Falling    

    !Destructor
    public  ::  KillCEQUALW2
    private ::      DeallocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjCEQUALW2
    
    !Interfaces
    private :: UngetCEQUALW2_1D_Int
    interface  UngetCEQUALW2
        module procedure UngetCEQUALW2_1D_Int
    end interface UngetCEQUALW2

    !Types---------------------------------------------------------------------
    
    private :: T_ID
    type       T_ID
        character(len=StringLength)                 :: Name
        integer                                     :: IDNumber
    end type   T_ID                                 
                                                    
    private :: T_Size                               
    type       T_Size                               
        integer                                     :: PropLB       = null_int
        integer                                     :: PropUB       = null_int
        integer                                     :: ArrayLB      = null_int
        integer                                     :: ArrayUB      = null_int
    end type   T_Size

    private :: T_RateIndex
    type       T_RateIndex
        integer                                     :: ANLim      
        integer                                     :: APLim              
        integer                                     :: ASLim      
        integer                                     :: ALightLim     
        integer                                     :: AOverallLim
        integer                                     :: ENLim      
        integer                                     :: EPLim      
        integer                                     :: ESLim      
        integer                                     :: ELightLim  
        integer                                     :: EOverallLim          
        integer                                     :: NH4D      
        integer                                     :: NO3D      
        integer                                     :: LDOMD     
        integer                                     :: RDOMD     
        integer                                     :: LPOMD     
        integer                                     :: RPOMD     
        integer                                     :: LRDOMD    
        integer                                     :: LRPOMD    
        integer                                     :: CBODD                                         
        integer                                     :: PO4ER     
        integer                                     :: PO4EG     
        integer                                     :: PO4AR     
        integer                                     :: PO4AG     
        integer                                     :: PO4OM     
        integer                                     :: PO4BOD    
        integer                                     :: NH4ER     
        integer                                     :: NH4EG     
        integer                                     :: NH4AR     
        integer                                     :: NH4AG     
        integer                                     :: NH4OM     
        integer                                     :: NH4BOD    
        integer                                     :: NO3AG     
        integer                                     :: NO3EG     
        integer                                     :: DSIAG     
        integer                                     :: DSIEG     
        integer                                     :: DSID      
        integer                                     :: PSIAM     
        integer                                     :: PSID      
        integer                                     :: LDOMAP    
        integer                                     :: LDOMEP    
        integer                                     :: LPOMAP    
        integer                                     :: DOAP      
        integer                                     :: DOEP      
        integer                                     :: DOAR      
        integer                                     :: DOER      
        integer                                     :: DOOM      
        integer                                     :: DONIT     
        integer                                     :: ICarbonAP 
        integer                                     :: ICarbonEP 
        integer                                     :: ICarbonBOD

    end type T_RateIndex

    private :: T_Rate                               
    type       T_Rate                               
        integer                                     :: LB       = null_int
        integer                                     :: UB       = null_int
        integer, pointer, dimension(:  )            :: Match
        real, pointer, dimension   (:,:)            :: Value
        logical                                     :: compute  =.false.
        type(T_RateIndex )                          :: MohidIndex
        type(T_RateIndex )                          :: CequalIndex    
    end type   T_Rate


    private :: T_PropIndex
    type       T_PropIndex
        integer                                     :: ICarbon      = null_int 
        integer                                     :: pomref       = null_int !Refractory Particulate Organic Matter        
        integer                                     :: pomlab       = null_int !Labile POM        
        integer                                     :: domlab       = null_int !Labile Dissolved OM        
        integer                                     :: domref       = null_int !Refractory DOM        
        integer                                     :: Ammonia      = null_int        
        integer                                     :: Nitrate      = null_int        
        integer                                     :: Phosphorus   = null_int     
        integer                                     :: Oxygen       = null_int         
        integer                                     :: BOD          = null_int !Biochemical Oxygen Demand               
        integer                                     :: sipart       = null_int !Particulate Silica        
        integer                                     :: sidiss       = null_int !Dissolved Silica        
        integer                                     :: pH           = null_int
        integer                                     :: CO2          = null_int
        integer                                     :: CO3          = null_int
        integer                                     :: HCO3         = null_int
        integer                                     :: Detritus     = null_int
    end type T_PropIndex
   
    

    private :: T_Compute
    type       T_Compute
        logical                                     ::  Algae           = OFF 
        logical                                     ::  Epiphyton       = OFF  
        logical                                     ::  Nitrogen        = OFF  
        logical                                     ::  Phosphorus      = OFF
        logical                                     ::  OrganicMatter   = OFF
        logical                                     ::  ICarbon         = OFF
        logical                                     ::  BOD             = OFF
        logical                                     ::  Oxygen          = ON
        logical                                     ::  Silica          = OFF
        logical                                     ::  Detritus        = OFF
    end type   T_Compute 
    
    private :: T_BenthicCompute
    type       T_BenthicCompute
        logical                                     ::  Detritus        = OFF 
        logical                                     ::  Ammonia         = OFF  
        logical                                     ::  Phosphorus      = OFF  
        logical                                     ::  Dsilica         = OFF
        logical                                     ::  Oxygen          = OFF
        logical                                     ::  ICarbon         = OFF
    end type   T_BenthicCompute      
    
    private :: T_SOD
    type       T_SOD
        logical                                     :: UseSOD
        real, pointer, dimension(:  )               :: Rate
        real                                        :: T1
        real                                        :: T2
        real                                        :: K1
        real                                        :: K2
        real                                        :: PO4R
        real                                        :: NH4R
        real                                        :: SiR
        real                                        :: CO2R 
        real                                        :: O2Sink
        real                                        :: TRM
        logical                                     :: DefaultO2                
    end type T_SOD
                 
    private :: T_External
    type       T_External
        real, pointer, dimension(:  )               :: Salinity
        real, pointer, dimension(:  )               :: Temperature
        real, pointer, dimension(:  )               :: Oxygen
        real, pointer, dimension(:  )               :: Alkalinity
        real, pointer, dimension(:  )               :: ShortWaveRadiation
        real, pointer, dimension(:  )               :: LightExtCoefField
        real, pointer, dimension(:  )               :: Thickness
        real, pointer, dimension(:,:)               :: Mass      
    end type T_External

    private :: T_Algae
    type       T_Algae
        type(T_ID)                                  :: ID
        integer                                     :: PropIndex
        
        ! Algal Maximum Biological Rates [day^-1]                               
        real                                        :: AG !Growth
        real                                        :: AR !Respiration     
        real                                        :: AE !Excretion
        real                                        :: AM !Mortality             
        
        !Algal saturating light intensity at maximum phtosynthetic rate [W m^-2]
        real                                        :: ASAT 

        !Algal half-saturation coefficients [g m^-3]
        real                                        :: AHSP  !for Phosphurus   
        real                                        :: AHSN  !for Nitrate + Ammonium   
        real                                        :: AHSSI !for Silica   
            
                                                
        !Algal Temperature Rate Coefficients    
        real                                        :: AT1 !Lower temperature for algal growth (ºC)
        real                                        :: AT2 !Lower temperature for maximum algal growth (ºC)  
        real                                        :: AT3 !Upper temperature for maximum algal growth     
        real                                        :: AT4 !Upper temperature for algal growth   
        real                                        :: AK1 !Fraction of algal growth rate at AT1    
        real                                        :: AK2 !Fraction of maximum algal growth rate at AT2     
        real                                        :: AK3 !Fraction of maximum algal growth rate at AT3     
        real                                        :: AK4 !Fraction of algal growth rate at AT4     
                                                    
        !Algal Stoichiometry                        
        real                                        :: AP   !Algal stoichiometric coefficient for phosphorus     
        real                                        :: AN   !Algal stoichiometric coefficient for nitrogen     
        real                                        :: AC   !Algal stoichiometric coefficient for carbon     
        real                                        :: ASI  !Algal stoichiometric coefficient for silica       
        real                                        :: APOM !Algal stoichiometric coefficient for POM 
        real                                        :: O2AR 
        real                                        :: O2AG  
        
        !Algal Ammonia Preference
        integer                                     :: ANEQN !Equation for preference factor (either 1 or 2)    
        real                                        :: ANPR  !half-saturation copreference constant  
                                                    
        !Algal Temperature Rate Multipliers [1]              
        real                                        :: ATRM  !Algae Temperature Rate Multiplier
        real                                        :: ATRMR !ATMR for rising limb of curve
        real                                        :: ATRMF !ATMR for falling limb of curve
                                                    
        !Algal light extintion [m^-1]                
        real                                        :: EXA
                                                
        !Algal Growth Limitations                     
        real,pointer, dimension(:)                  :: NLim !Nitrogen
        real,pointer, dimension(:)                  :: PLim !Phosphorus
        real,pointer, dimension(:)                  :: SLim !Silica
        real,pointer, dimension(:)                  :: LightLim 
        real,pointer, dimension(:)                  :: OverallLim
                                                    
        !Algal Biological Rates [day^-1]                                     
        real                                        :: AGR !Growth
        real                                        :: ARR !Respiration
        real                                        :: AER !Excretion
        real                                        :: AMR !Mortality

        !Algal Source/sink Term 
        real                                        :: ASS 
                                                    
        !Collection of algae                        
        type(T_Algae), pointer                      :: Next
    end type T_Algae                                
        
        
                                                    
    private :: T_Epiphyton                          
    type       T_Epiphyton                          
        type(T_ID)                                  :: ID
        integer                                     :: PropIndex

        ! Epiphyte Maximum Biological Rates [day^-1]                       
        real                                        :: EG !Growth   
        real                                        :: ER !Respiration   
        real                                        :: EE !Excretion   
        real                                        :: EM !Mortality   
        
        
        !Epiphyte half-saturation coefficients [g m^-3]
        real                                        :: EHSP  !Half-saturation coefficient for phosphorus
        real                                        :: EHSN  !for nitrates
        real                                        :: EHSSI !for silica
                                                    
                      
        !Epiphyte saturating light intensity at maximum phtosynthetic rate [W m^-2]
        real                                        :: ESAT     
        
        !Epiphyte Amonia preference
        real                                        :: ENPR  !half-saturation coefficient [gm^-3]   
        integer                                     :: ENEQN !Equation for preference factor (1 or 2)
                                                    
        ! Epiphyte Temperature Rate Coefficients    
        real                                        :: ET1 
        real                                        :: ET2 
        real                                        :: ET3 
        real                                        :: ET4 
        real                                        :: EK1 
        real                                        :: EK2 
        real                                        :: EK3 
        real                                        :: EK4 
                                                
        !Epiphyte Stoichiometry                 
        real                                        :: EP  !Epiphyte stoichiom. coef. for phosphorus    
        real                                        :: EN  !Epiphyte stoichiom. coef. for nitrogen   
        real                                        :: EC  !Epiphyte stoichiom. coef. for carbon   
        real                                        :: ESI !Epiphyte stoichiom. coef. for silica      
        real                                        :: EPOM !Epiphyte stoichiom. coef. for POM
        real                                        :: O2ER
        real                                        :: O2EG     
                                                    
        !Epiphyte Temperature Rate Multipliers               
        real                                        :: ETRM  !Temperature Rate Multiplier  
        real                                        :: ETRMR !TRM for rising limb of curve
        real                                        :: ETRMF !TRM for falling lim of curve
                                                    
        !Epiphyte Growth Limitations                         
        real,pointer, dimension(:)                  :: NLim  !Nitrogen
        real,pointer, dimension(:)                  :: PLim  !Phosphorus
        real,pointer, dimension(:)                  :: SLim  !Silica
        real,pointer, dimension(:)                  :: LightLim
        real,pointer, dimension(:)                  :: OverallLim
                                                
        !Epiphyte Biological Rates [s^-1]                                
        real                                        :: EGR  !Growth
        real                                        :: ERR  !Respiration
        real                                        :: EER  !Excretion
        real                                        :: EMR  !Mortality

        !Epiphyte Source/sink Term [g m^-3 s^-1]
        real                                        :: ESS  
        
        !Collection of epiphyton
        type(T_Epiphyton), pointer                  :: Next
    end type T_Epiphyton



    private :: T_CEQUALW2
    type      T_CEQUALW2
        private
        integer                                     :: InstanceID
        integer, dimension(:), pointer              :: PropertyList
        type(T_Size      )                          :: Size
        type(T_PropIndex )                          :: PropIndex
        type(T_Rate      )                          :: Rate
        type(T_Compute   )                          :: Compute
        type(T_BenthicCompute)                      :: BenthicCompute
        type(T_SOD       )                          :: SOD
        type(T_External  )                          :: ExternalVar
        type(T_Algae     ),pointer                  :: FirstAlgae
        type(T_Epiphyton ),pointer                  :: FirstEpiphyton
        real                                        :: DTDay     = null_real
        real                                        :: DTSecond  = null_real
        real, dimension(:), pointer                 :: SinksSources 
                                                            
        !Extinction Coefficients                            
        !real                                        ::  EXH2O    = null_real !Extinction for pure water, m-1
        !real                                        ::  EXSS     = null_real !Extinction due to inorganic suspended solids, m-1
        !real                                        ::  EXOM     = null_real !Extinction due to organic suspended solids, m-1
        !Logical                                     ::  EXC      = .false.   !Read extinction coefficients, ON or OFF
        !Logical                                     ::  EXIC     = .false.   !Interpolate extinction coefficients, ON or OFF
                                                            
        !Dissolved Organic Matter                           
        real                                        ::  LDOMDK   = null_real !LDOM decay rate [day^-1]
        real                                        ::  RDOMDK   = null_real !RDOM decay rate [day^-1]
        real                                        ::  LRDDK    = null_real !Labile to refractory DOM dekay rate [day^-1]
                                                     
        !Particulate Organic Matter                  
        real                                        ::  LPOMDK   = null_real !LPOM decay rate [day^-1]
        real                                        ::  RPOMDK   = null_real !RPOM decay rate [day^-1]
        real                                        ::  LRPDK    = null_real
                                                            
        !Organic Matter Stoichiometry                       
        real                                        ::  ORGP     = null_real !Stoichiometric coef. for phosphorus
        real                                        ::  ORGN     = null_real !Stoichiometric coef. for nitrogen
        real                                        ::  ORGC     = null_real !Stoichiometric coef. for carbon
        real                                        ::  ORGSI    = null_real !Stoichiometric coef. for silica

        !Organic Matter Temperature Rate Multipliers
        real                                        ::  OMT1     = null_real 
        real                                        ::  OMT2     = null_real 
        real                                        ::  OMK1     = null_real 
        real                                        ::  OMK2     = null_real

        
        !Detritus Temperature Rate Multipliers
        real                                        ::  DETTRM    = null_real
        real                                        ::  DETT1     = null_real 
        real                                        ::  DETT2     = null_real 
        real                                        ::  DETK1     = null_real 
        real                                        ::  DETK2     = null_real
        real                                        ::  SDK       = null_real 
 
                                                    
        !Carbonaceous Biochemical Oxygen Demand     
        real                                        ::  KBOD     = null_real !CBOD Decay Rate [day^-1]
        real                                        ::  TBOD     = null_real !BOD Temperature Rate Multiplier
        real                                        ::  RBOD     = null_real
                                                    
        !CBOD Stoichiometry                         
        real                                        ::  BODP     = null_real !Phosphorus/CBOD stochiometric ratio
        real                                        ::  BODN     = null_real !Nitrogen/CBOD stochiometric ratio
        real                                        ::  BODC     = null_real !Carbon/CBOD stochiometric ratio
                                                            
                                                            
        !Ammonium                                   
        real                                        ::  NH4DK    = null_real !Ammonium Decay Rate [day^-1]
                                                    
        !Ammonium Temperature Rate Multipliers      
        real                                        ::  NH4T1    = null_real !Minimum Temperature (T1)
        real                                        ::  NH4T2    = null_real !Optimal Temperature (T2)
        real                                        ::  NH4K1    = null_real !Multiplier factor for T1
        real                                        ::  NH4K2    = null_real !Multiplier factor for T2
                                                    
        !Nitrate                                    
        real                                        ::  NO3DK    = null_real !Nitrate Decay Rate [day^-1]
                                                                    
        !Nitrate Temperature Rate Multipliers               
        real                                        ::  NO3T1    = null_real !Minimum temperature (T1)
        real                                        ::  NO3T2    = null_real !Optimal temperature (T2)
        real                                        ::  NO3K1    = null_real !Multiplier factor for T1
        real                                        ::  NO3K2    = null_real !Multiplier factor for T2
                                                    
        !Silica                                     
        real                                        ::  PSIDK    = null_real
                                                            
                                                            
                                                           
        !Oxygen Stoichiometry 1                     
        real                                        ::  O2NH4    = null_real
        real                                        ::  O2OM     = null_real
                                                    
       
                                                    
        !Oxygen Limit                               
        real                                        ::  O2LIM    = null_real
                                                    
        !Oxygen Kinetic Flux                        
        !real                                        ::  DOAE     = null_real
                                                            
                       
        !More Temperature Rate Multipliers          
        real                                        ::  NH4TRM   = null_real
        real                                        ::  NO3TRM   = null_real
        real                                        ::  OMTRM    = null_real
                                                    
        !Auxiliar variables and arrays              
        real                                        ::  DO1      = null_real
        real                                        ::  DO2      = null_real
        real                                        ::  DO3      = null_real
                                                    
        !Decay Rates                                
        real                                        ::  NH4D     = null_real
        real                                        ::  NO3D     = null_real
        real                                        ::  LDOMD    = null_real
        real                                        ::  RDOMD    = null_real
        real                                        ::  LPOMD    = null_real
        real                                        ::  RPOMD    = null_real
        real                                        ::  LRDOMD   = null_real
        real                                        ::  LRPOMD   = null_real
        real                                        ::  CBODD    = null_real
        !real                                        ::  CBODDK   = null_real
        real                                        ::  DETD     = null_real
        real                                        ::  SODD     = null_real
        !Instance of Module_EnterData
        integer                                     :: ObjEnterData = 0
        
        character(StringLength)                     :: Model
        type(T_CEQUALW2), pointer                   :: Next

    end type T_CEQUALW2
    
    !Global Module Variables
    type (T_CEQUALW2), pointer                      :: FirstObjCEQUALW2
    type (T_CEQUALW2), pointer                      :: Me

    !--------------------------------------------------------------------------
    
    contains



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartCEQUALW2(CEQUALW2_ID, ArrayLB, ArrayUB, FileName, Model, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: CEQUALW2_ID
        integer, optional, intent(IN )      :: ArrayLB, ArrayUB
        character(len=*)                    :: FileName
        character(StringLength),intent(IN ) :: Model
        integer, optional, intent(OUT)      :: STAT     

        !External--------------------------------------------------------------
        integer                             :: ready_         

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
         
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mCEQUALW2_)) then
            nullify (FirstObjCEQUALW2)
            call RegisterModule (mCEQUALW2_) 
        endif

        
        
        call Ready(CEQUALW2_ID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%Size%ArrayLB = ArrayLB
            Me%Size%ArrayUB = ArrayUB
            Me%Model        = Model
          select case (Model)

          case (CEQUALW2Model)
            
            call ReadWaterColumnData (FileName)
          
          case (BenthicCEQUALW2Model)

            call ReadBenthicData (FileName)

          case default
                write(*,*) 
                write(*,*) 'Defined sinks and sources model was not recognised.'
                stop 'StartCEQUALW2 - ModuleCeQualW2 - ERR01' 
          end select

            !Returns ID
            CEQUALW2_ID = Me%InstanceID

            STAT_ = SUCCESS_
        else 
            
            stop 'StartCEQUALW2 - ModuleCEQUALW2 - ERR00' 

        end if cd0


        if (present(STAT))STAT = STAT_

    end subroutine StartCEQUALW2

    !--------------------------------------------------------------------------

    subroutine AllocateInstance 

       
        !Local-----------------------------------------------------------------
        type (T_CEQUALW2), pointer           :: NewObjCEQUALW2
        type (T_CEQUALW2), pointer           :: PreviousObjCEQUALW2


        !Allocates new instance
        allocate (NewObjCEQUALW2)
        nullify  (NewObjCEQUALW2%Next)

        nullify (NewObjCEQUALW2%FirstAlgae)
        nullify (NewObjCEQUALW2%FirstEpiphyton)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjCEQUALW2)) then
            FirstObjCEQUALW2        => NewObjCEQUALW2
            Me                      => NewObjCEQUALW2
        else
            PreviousObjCEQUALW2     => FirstObjCEQUALW2
            Me                      => FirstObjCEQUALW2%Next
            do while (associated(Me))
                PreviousObjCEQUALW2 => Me
                Me                  => Me%Next
            enddo
            Me                      => NewObjCEQUALW2
            PreviousObjCEQUALW2%Next=> NewObjCEQUALW2
        endif

        Me%InstanceID = RegisterNewInstance (mCEQUALW2_)

    end subroutine AllocateInstance

    
    !--------------------------------------------------------------------------

  
    subroutine PropertyIndexNumber

       
        !Local-----------------------------------------------------------------
        type(T_Algae),      pointer             :: Algae
        type(T_Epiphyton),  pointer             :: Epiphyton
        integer                                 :: Index
        !Local-----------------------------------------------------------------
        
        Me%Size%PropLB = 1
        Me%Size%PropUB = 0

        Index                   = 0

        !Algae index number
        if (Me%Compute%Algae) then
            
            Algae => Me%FirstAlgae
            do while(associated(Algae))
                
                Index                   = Index + 1
                Algae%PropIndex         = Index
                Me%Size%PropUB = Me%Size%PropUB + 1

                Algae => Algae%Next
            end do
        endif   
        
        !Epiphyton index number
        if (Me%Compute%Epiphyton) then

            Epiphyton => Me%FirstEpiphyton
            do while(associated(Epiphyton))
                
                Index                   = Index + 1
                Epiphyton%PropIndex     = Index
                Me%Size%PropUB = Me%Size%PropUB + 1

                Epiphyton => Epiphyton%Next
            end do
        endif   


        !Nitrogen index number
        if (Me%Compute%Nitrogen) then
            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%Ammonia       = Me%Size%PropUB

            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%Nitrate       = Me%Size%PropUB
        endif

        !Phosphorus index number
        if (Me%Compute%Phosphorus) then   
            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%Phosphorus    = Me%Size%PropUB
        endif   


        !OrganicMatter index number
        if (Me%Compute%OrganicMatter) then
            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%pomref        = Me%Size%PropUB

            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%pomlab        = Me%Size%PropUB

            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%domlab        = Me%Size%PropUB

            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%domref        = Me%Size%PropUB
        endif

        if (Me%Compute%Silica) then
            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%sipart        = Me%Size%PropUB
      
            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%sidiss        = Me%Size%PropUB
        endif

        
        if (Me%Compute%ICarbon) then

            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%ICarbon       = Me%Size%PropUB

            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%CO2           = Me%Size%PropUB

            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%pH            = Me%Size%PropUB

            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%HCO3          = Me%Size%PropUB

            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%CO3           = Me%Size%PropUB

        endif

        !Oxygen index number -> The oxygen is always calculated.
        Me%Size%PropUB                 = Me%Size%PropUB + 1
        Me%PropIndex%Oxygen            = Me%Size%PropUB

        !BOD index number
        if (Me%Compute%BOD) then   
            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%BOD           = Me%Size%PropUB
        endif

        !Detritus index number
        if (Me%Compute%Detritus) then   
            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%Detritus           = Me%Size%PropUB
        endif




        !----------------------------------------------------------------------

    end subroutine PropertyIndexNumber


    !--------------------------------------------------------------------------

  
    subroutine RateIndexNumber

       
        !Local-----------------------------------------------------------------
        integer                                 :: countrate, i
        integer                                 :: STAT_CALL
        integer, dimension(100)                 :: LocalMatch
        
        countrate = 0
        
        
        if (CheckPropertyName('ANLIM',        number = Me%rate%MohidIndex%ANLim      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%ANLim      
          Me%Rate%CeQualIndex%ANLim      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR1'
        endif


        if (CheckPropertyName('APLIM',        number = Me%rate%MohidIndex%APLim      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%APLim      
          Me%Rate%CeQualIndex%APLim      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR2'
        endif


        if (CheckPropertyName('ASLIM',        number = Me%rate%MohidIndex%ASLim      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%ASLim      
          Me%Rate%CeQualIndex%ASLim      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR3'
        endif


        if (CheckPropertyName('ALIGHTLIM',        number = Me%rate%MohidIndex%ALightLim  )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%ALightLim  
          Me%Rate%CeQualIndex%ALightLim  = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR4'
        endif


        if (CheckPropertyName('AOVERALLLIM',       number = Me%rate%MohidIndex%AOverallLim)) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%AOverallLim
          Me%Rate%CeQualIndex%AOverallLim= countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR6'
        endif


        if (CheckPropertyName( 'ENLIM',        number = Me%rate%MohidIndex%ENLim      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%ENLim      
          Me%Rate%CeQualIndex%ENLim      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR7'
        endif


        if (CheckPropertyName( 'EPLIM',        number = Me%rate%MohidIndex%EPLim      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%EPLim      
          Me%Rate%CeQualIndex%EPLim      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR8'
        endif


        if (CheckPropertyName( 'ESLIM',        number = Me%rate%MohidIndex%ESLim      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%ESLim      
          Me%Rate%CeQualIndex%ESLim      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR9'
        endif


        if (CheckPropertyName( 'ELIGHTLIM',        number = Me%rate%MohidIndex%ELightLim  )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%ELightLim  
          Me%Rate%CeQualIndex%ELightLim  = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR10'
        endif


        if (CheckPropertyName( 'EOVERALLLIM',        number = Me%rate%MohidIndex%EOverallLim)) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%EOverallLim
          Me%Rate%CeQualIndex%EOverallLim= countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR12'
        endif


        if (CheckPropertyName('NH4D',        number = Me%rate%MohidIndex%NH4D      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%NH4D      
          Me%Rate%CeQualIndex%NH4D      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR13'
        endif


        if (CheckPropertyName('NO3D',        number = Me%rate%MohidIndex%NO3D      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%NO3D      
          Me%Rate%CeQualIndex%NO3D      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR14'
        endif


        if (CheckPropertyName('LDOMD',        number = Me%rate%MohidIndex%LDOMD     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%LDOMD     
          Me%Rate%CeQualIndex%LDOMD     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR15'
        endif


        if (CheckPropertyName('RDOMD',        number = Me%rate%MohidIndex%RDOMD     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%RDOMD     
          Me%Rate%CeQualIndex%RDOMD     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR16'
        endif


        if (CheckPropertyName('LPOMD',        number = Me%rate%MohidIndex%LPOMD     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%LPOMD     
          Me%Rate%CeQualIndex%LPOMD     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR17'
        endif


        if (CheckPropertyName('RPOMD',        number = Me%rate%MohidIndex%RPOMD     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%RPOMD     
          Me%Rate%CeQualIndex%RPOMD     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR18'
        endif


        if (CheckPropertyName('LRDOMD',        number = Me%rate%MohidIndex%LRDOMD    )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%LRDOMD    
          Me%Rate%CeQualIndex%LRDOMD    = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR19'
        endif


        if (CheckPropertyName('LRPOMD',        number = Me%rate%MohidIndex%LRPOMD    )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%LRPOMD    
          Me%Rate%CeQualIndex%LRPOMD    = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR20'
        endif


        if (CheckPropertyName('CBODD',        number = Me%rate%MohidIndex%CBODD      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%CBODD      
          Me%Rate%CeQualIndex%CBODD      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR21'
        endif


        if (CheckPropertyName('PO4ER',        number = Me%rate%MohidIndex%PO4ER     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%PO4ER     
          Me%Rate%CeQualIndex%PO4ER     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR22'
        endif


        if (CheckPropertyName('PO4EG',        number = Me%rate%MohidIndex%PO4EG     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%PO4EG     
          Me%Rate%CeQualIndex%PO4EG     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR23'
        endif


        if (CheckPropertyName('PO4AR',        number = Me%rate%MohidIndex%PO4AR     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%PO4AR     
          Me%Rate%CeQualIndex%PO4AR     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR24'
        endif


        if (CheckPropertyName('PO4AG',        number = Me%rate%MohidIndex%PO4AG     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%PO4AG     
          Me%Rate%CeQualIndex%PO4AG     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR25'
        endif


        if (CheckPropertyName('PO4OM',        number = Me%rate%MohidIndex%PO4OM     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%PO4OM     
          Me%Rate%CeQualIndex%PO4OM     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR26'
        endif


        if (CheckPropertyName('PO4BOD',        number = Me%rate%MohidIndex%PO4BOD    )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%PO4BOD    
          Me%Rate%CeQualIndex%PO4BOD    = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR27'
        endif


        if (CheckPropertyName('NH4ER',        number = Me%rate%MohidIndex%NH4ER     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%NH4ER     
          Me%Rate%CeQualIndex%NH4ER     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR28'
        endif


        if (CheckPropertyName('NH4EG',        number = Me%rate%MohidIndex%NH4EG     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%NH4EG     
          Me%Rate%CeQualIndex%NH4EG     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR29'
        endif


        if (CheckPropertyName('NH4AR',        number = Me%rate%MohidIndex%NH4AR     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%NH4AR     
          Me%Rate%CeQualIndex%NH4AR     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR30'
        endif


        if (CheckPropertyName('NH4AG',        number = Me%rate%MohidIndex%NH4AG     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%NH4AG     
          Me%Rate%CeQualIndex%NH4AG     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR31'
        endif


        if (CheckPropertyName('NH4OM',        number = Me%rate%MohidIndex%NH4OM     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%NH4OM     
          Me%Rate%CeQualIndex%NH4OM     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR32'
        endif


        if (CheckPropertyName('NH4BOD',        number = Me%rate%MohidIndex%NH4BOD    )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%NH4BOD    
          Me%Rate%CeQualIndex%NH4BOD    = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR33'
        endif


        if (CheckPropertyName('NO3AG',        number = Me%rate%MohidIndex%NO3AG     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%NO3AG     
          Me%Rate%CeQualIndex%NO3AG     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR34'
        endif


        if (CheckPropertyName('NO3EG',        number = Me%rate%MohidIndex%NO3EG     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%NO3EG     
          Me%Rate%CeQualIndex%NO3EG     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR35'
        endif


        if (CheckPropertyName('DSIAG',        number = Me%rate%MohidIndex%DSIAG     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%DSIAG     
          Me%Rate%CeQualIndex%DSIAG     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR36'
        endif


        if (CheckPropertyName('DSIEG',        number = Me%rate%MohidIndex%DSIEG     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%DSIEG     
          Me%Rate%CeQualIndex%DSIEG     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR37'
        endif


        if (CheckPropertyName('DSID ',        number = Me%rate%MohidIndex%DSID      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%DSID      
          Me%Rate%CeQualIndex%DSID      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR38'
        endif


        if (CheckPropertyName('PSIAM',        number = Me%rate%MohidIndex%PSIAM     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%PSIAM     
          Me%Rate%CeQualIndex%PSIAM     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR39'
        endif


        if (CheckPropertyName('PSID ',        number = Me%rate%MohidIndex%PSID      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%PSID      
          Me%Rate%CeQualIndex%PSID      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR40'
        endif


        if (CheckPropertyName('LDOMAP',        number = Me%rate%MohidIndex%LDOMAP    )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%LDOMAP    
          Me%Rate%CeQualIndex%LDOMAP    = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR41'
        endif


        if (CheckPropertyName('LDOMEP',        number = Me%rate%MohidIndex%LDOMEP    )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%LDOMEP    
          Me%Rate%CeQualIndex%LDOMEP    = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR42'
        endif


        if (CheckPropertyName('LPOMAP',        number = Me%rate%MohidIndex%LPOMAP    )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%LPOMAP    
          Me%Rate%CeQualIndex%LPOMAP    = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR43'
        endif


        if (CheckPropertyName('DOAP',        number = Me%rate%MohidIndex%DOAP      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%DOAP      
          Me%Rate%CeQualIndex%DOAP      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR44'
        endif


        if (CheckPropertyName('DOEP',        number = Me%rate%MohidIndex%DOEP      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%DOEP      
          Me%Rate%CeQualIndex%DOEP      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR45'
        endif


        if (CheckPropertyName('DOAR',        number = Me%rate%MohidIndex%DOAR      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%DOAR      
          Me%Rate%CeQualIndex%DOAR      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR46'
        endif


        if (CheckPropertyName('DOER',        number = Me%rate%MohidIndex%DOER      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%DOER      
          Me%Rate%CeQualIndex%DOER      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR47'
        endif


        if (CheckPropertyName('DOOM',        number = Me%rate%MohidIndex%DOOM      )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%DOOM      
          Me%Rate%CeQualIndex%DOOM      = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR48'
        endif


        if (CheckPropertyName('DONIT',        number = Me%rate%MohidIndex%DONIT     )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%DONIT     
          Me%Rate%CeQualIndex%DONIT     = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR49'
        endif


        if (CheckPropertyName('ICARBONAP',        number = Me%rate%MohidIndex%ICarbonAP )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%ICarbonAP 
          Me%Rate%CeQualIndex%ICarbonAP = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR50'
        endif


        if (CheckPropertyName('ICARBONEP',        number = Me%rate%MohidIndex%ICarbonEP )) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%ICarbonEP 
          Me%Rate%CeQualIndex%ICarbonEP = countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR51'
        endif


        if (CheckPropertyName('ICARBONBOD',        number = Me%rate%MohidIndex%ICarbonBOD)) then
          countrate = countrate +1
          LocalMatch(countrate) = Me%Rate%MohidIndex%ICarbonBOD
          Me%Rate%CeQualIndex%ICarbonBOD= countrate
        else
          stop 'RateIndexNumber - ModuleCEQUALW2 - ERR52'
        endif

        Me%Rate%LB =1
        Me%Rate%UB = countrate
        
        allocate(Me%Rate%Value(Me%Rate%LB:Me%Rate%UB,Me%Size%ArrayLB:Me%Size%ArrayUB), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "PropertyIndexNumber - ModuleCEQUALW2 - ERR53"
    
        allocate (ME%Rate%Match(countrate))

        do i=1,countrate
          ME%Rate%Match(i) = LocalMatch(i)
        enddo

        if (countrate.ne.0) then
           Me%Rate%Compute=.true.
        endif
        
    end subroutine RateIndexNumber

        !--------------------------------------------------------------------------       
    
    subroutine ReadBenthicData(FileName)

        !Arguments-------------------------------------------------------------
        character(len=*)                     :: FileName

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadBenthicData - ModuleCEQUALW2 - ERR01'

        call ReadBenthicGlobalVariables  
        
        call ReadSODData  

        call ReadOMParameters

        call ReadOxygenParameters
        
        call ReadDetritusParameters        

        call ReadAmmoniaParameters      

        call ReadPhosphorusParameters 
        
        call ReadDsilicaParameters          

        call ReadICarbonParameters
        
        call BenthicPropertyIndexNumber
        
        !call RateIndexNumber   
        
        call ConstructBenthicPropertyList 

        
        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadBenthicData - ModuleCEQUALW2 - ERR03'

    end subroutine ReadBenthicData

    !--------------------------------------------------------------------------
    

    
    subroutine ReadBenthicGlobalVariables

        
        !External--------------------------------------------------------------
        integer                                     :: FromFile
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: flag

        !--------------------------------------------------------------------------

        call GetExtractType(FromFile = FromFile)

        if (Me%DTSecond .LE. 0.0) then

            !DTSecond, time step, in seconds, between two CEQUALW2 calls 
            call GetData(Me%DTSecond,                                       &
                         Me%ObjEnterData, flag,                             &
                         SearchType     = FromFile,                         &
                         keyword        ='DTSECONDS',                       & 
                         Default        = 60.0 * 60.0,                      &
                         ClientModule   = MohidModules(mCEQUALW2_)%Name,    &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'ReadBenthicGlobalVariables - ModuleCEQUALW2 - ERR01' 

            if (flag .EQ. 0) then
                write(*,*) 
                write(*,*) 'Keyword DTSECONDS not found in CEQUALW2 data file.'
                write(*,*) 'ReadBenthicGlobalVariables - ModuleCEQUALW2 - WRN01'
                write(*,*) 'Assumed ', Me%DTSecond, 'seconds (',  Me%DTSecond / 3600.0, 'hour).'
                write(*,*) 
            end if
        end if


        !For compatibility with the rest of the program
        Me%DTDay = Me%DTSecond / 24.0 / 60.0 / 60.0
     
       
   
        

    end subroutine ReadBenthicGlobalVariables

    ! -------------------------------------------------------------------------

    subroutine ReadSODData
        
        !External--------------------------------------------------------------
        integer                                     :: FromBlock, ClientNumber
        logical                                     :: BlockFound
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: flag

        !--------------------------------------------------------------------------

        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,              &
                                        '<begin_SOD>', '<end_SOD>', BlockFound, &
                                        STAT = STAT_CALL)
        if(STAT_CALL .EQ. SUCCESS_)then    
            if (BlockFound) then
                
                Me%SOD%UseSOD = .true.
                
                call GetExtractType(FromBlock = FromBlock)
                
                !Temp 1 for temperature rate multiplier
                call GetData(   Me%SOD%T1,                                          &
                                Me%ObjEnterData, flag,                              &
                                SearchType    =  FromBlock,                         &
                                keyword       = 'SODT1',                            &
                                default       =  4.   ,                             &
                                ClientModule  =  MohidModules(mCEQUALW2_)%Name,     &
                                STAT          =  STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop "ReadSODData - ModuleCEQUALW2 - ERR01"
                
                !Temp 2 for temperature rate multiplier
                call GetData(   Me%SOD%T2,                                          &
                                Me%ObjEnterData, flag,                              &
                                SearchType    =  FromBlock,                         &
                                keyword       = 'SODT2',                            &
                                default       =  35.   ,                            &
                                ClientModule  =  MohidModules(mCEQUALW2_)%Name,     &
                                STAT          =  STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop "ReadSODData - ModuleCEQUALW2 - ERR02"
                
                !K1 for temperature rate multiplier
                call GetData(   Me%SOD%K1,                                          &
                                Me%ObjEnterData, flag,                              &
                                SearchType    =  FromBlock,                         &
                                keyword       = 'SODK1',                            &
                                default       =  0.1   ,                            &
                                ClientModule  =  MohidModules(mCEQUALW2_)%Name,     &
                                STAT          =  STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop "ReadSODData - ModuleCEQUALW2 - ERR03"
                
                !K1 for temperature rate multiplier
                call GetData(   Me%SOD%K2,                                          &
                                Me%ObjEnterData, flag,                              &
                                SearchType    =  FromBlock,                         &
                                keyword       = 'SODK2',                            &
                                default       =  0.99  ,                            &
                                ClientModule  =  MohidModules(mCEQUALW2_)%Name,     &
                                STAT          =  STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop "ReadSODData - ModuleCEQUALW2 - ERR04"
                
                !Release of PO4 by SOD rate
                call GetData(   Me%SOD%PO4R,                                        &
                                Me%ObjEnterData, flag,                              &
                                SearchType    =  FromBlock,                         &
                                keyword       = 'PO4R',                             &
                                default       =  0.001,                             &
                                ClientModule  =  MohidModules(mCEQUALW2_)%Name,     &
                                STAT          =  STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop "ReadSODData - ModuleCEQUALW2 - ERR05"
                
                !Release of NH4 by SOD rate
                call GetData(   Me%SOD%NH4R,                                        &
                                Me%ObjEnterData, flag,                              &
                                SearchType    =  FromBlock,                         &
                                keyword       = 'NH4R',                             &
                                default       =  0.001,                             &
                                ClientModule  =  MohidModules(mCEQUALW2_)%Name,     &
                                STAT          =  STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop "ReadSODData - ModuleCEQUALW2 - ERR06"
                
                !Release of Silica by SOD rate
                call GetData(   Me%SOD%SiR,                                         &
                                Me%ObjEnterData, flag,                              &
                                SearchType    =  FromBlock,                         &
                                keyword       = 'SiR',                              &
                                default       =  0.1  ,                             &
                                ClientModule  =  MohidModules(mCEQUALW2_)%Name,     &
                                STAT          =  STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop "ReadSODData - ModuleCEQUALW2 - ERR07"

                !Release of CO2 by SOD rate
                call GetData(   Me%SOD%CO2R,                                        &
                                Me%ObjEnterData, flag,                              &
                                SearchType    =  FromBlock,                         &
                                keyword       = 'CO2R',                             &
                                default       =  0.1  ,                             &
                                ClientModule  =  MohidModules(mCEQUALW2_)%Name,     &
                                STAT          =  STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop "ReadSODData - ModuleCEQUALW2 - ERR08"


                !Sink of O2 by SOD rate
                call GetData(   Me%SOD%O2Sink,                                      &
                                Me%ObjEnterData, flag,                              &
                                SearchType    =  FromBlock,                         &
                                keyword       = 'O2Consumption',                    &
                                default       =  1.  ,                              &
                                ClientModule  =  MohidModules(mCEQUALW2_)%Name,     &
                                STAT          =  STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop "ReadSODData - ModuleCEQUALW2 - ERR09"                
                
                !In CEQUAL SOD uptakes O2 even above O2 min if this is turned off
                !SOD will inly use O2 below O2 min
                call GetData(   Me%SOD%DefaultO2,                                   &
                                Me%ObjEnterData, flag,                              &
                                SearchType    =  FromBlock,                         &
                                keyword       = 'DefaultO2',                        &
                                default       =  .true.   ,                         &
                                ClientModule  =  MohidModules(mCEQUALW2_)%Name,     &
                                STAT          =  STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop "ReadSODData - ModuleCEQUALW2 - ERR10"
              
                
            else                
                 Me%SOD%UseSOD = .false.                          
                                
            endif
            
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if(STAT_CALL .ne. SUCCESS_) stop "ReadSODData - ModuleCEQUALW2 - ERR011"        
        endif
                              
    end subroutine ReadSODData

    ! -------------------------------------------------------------------------
       
    subroutine ReadDetritusParameters

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: FromBlock

        !Local-----------------------------------------------------------------
        logical                                 :: BlockFound
        integer                                 :: ClientNumber
        integer                                 :: flag

        !----------------------------------------------------------------------

do1 :   do

            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,              &
                                        '<begin_det>', '<end_det>', BlockFound,     &
                                        STAT = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then    
cd2 :           if (BlockFound) then
                    
                    Me%Compute%Detritus =        .True.
                    Me%BenthicCompute%Detritus = .True. 
                    
                    call GetExtractType(FromBlock = FromBlock)
 

                    !SDK Sediment decay rate,  [day^-1]
                    call GetData(Me%SDK,                                                &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'DET_DECAY',                           &
                                 default       =  0.1 ,                                 &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadDetritusParameters - ModuleCEQUALW2 - ERR01"

                     
                    call GetData(Me%DETT1,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'DET_T1',                              &
                                 default       =  4.   ,                                 &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadDetritusParameters - ModuleCEQUALW2 - ERR02"

                    call GetData(Me%DETT2,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'DET_T2',                              &
                                 default       =  30.   ,                                 &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadDetritusParameters - ModuleCEQUALW2 - ERR03"

                    call GetData(Me%DETK1,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'DET_K1',                              &
                                 default       =  0.1   ,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadDetritusParameters - ModuleCEQUALW2 - ERR04"

                    call GetData(Me%DETK2,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'DET_K2',                              &
                                 default       =  0.99   ,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadDetritusParameters - ModuleCEQUALW2 - ERR05"

                    else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadDetritusParameters - ModuleCEQUALW2 - ERR006"

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                stop "ReadDetritusParameters - ModuleCEQUALW2 - ERR007"
       
            end if cd1
       
        end do do1


    end subroutine ReadDetritusParameters

       
   ! -------------------------------------------------------------------------
       
    subroutine ReadAmmoniaParameters

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: FromBlock

        !Local-----------------------------------------------------------------
        logical                                 :: BlockFound
        integer                                 :: ClientNumber
        integer                                 :: flag

        !----------------------------------------------------------------------

do1 :   do

            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,              &
                                        '<begin_ammonia>', '<end_ammonia>', BlockFound,     &
                                        STAT = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then    
cd2 :           if (BlockFound) then
                    
                    Me%BenthicCompute%Ammonia = .true.
                    
                    call GetExtractType(FromBlock = FromBlock)
 

                    !ORGN - Stoichiometric equivalent between organic matter and nitrogen  
                    call GetData(Me%ORGN,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'OM_STOICHIOMETRY_N',                  &
                                 default       =  0.0800,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadAmmoniaParameters - ModuleCEQUALW2 - ERR01"
                 
                    else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadAmmoniaParameters - ModuleCEQUALW2 - ERR002"

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                stop "ReadAmmoniaParameters - ModuleCEQUALW2 - ERR007"
       
            end if cd1
       
        end do do1


    end subroutine ReadAmmoniaParameters
    
    !--------------------------------------------------------------------------       
    
        subroutine ReadPhosphorusParameters

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: FromBlock

        !Local-----------------------------------------------------------------
        logical                                 :: BlockFound
        integer                                 :: ClientNumber
        integer                                 :: flag

        !----------------------------------------------------------------------

do1 :   do

            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,              &
                                        '<begin_phos>', '<end_phos>', BlockFound,     &
                                        STAT = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then    
cd2 :           if (BlockFound) then
                    
                    Me%BenthicCompute%Phosphorus = .True.
                    Me%Compute%Phosphorus = .True.
                    
                    call GetExtractType(FromBlock = FromBlock)
 

                select case (Me%Model)

                  case (BenthicCEQUALW2Model)
                    !ORGN - Stoichiometric equivalent between organic matter and phosphorus  
                    call GetData(Me%ORGP,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'OM_STOICHIOMETRY_P',                  &
                                 default       =  0.0050,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadPhosphorusParameters - ModuleCEQUALW2 - ERR01"
                   end select
                    else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadPhosphorusParameters - ModuleCEQUALW2 - ERR002"

                    exit do1    !No more blocks
                
                  
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                stop "ReadPhosphorusParameters - ModuleCEQUALW2 - ERR007"
       
            end if cd1
       
        end do do1


    end subroutine ReadPhosphorusParameters
    
    !-------------------------------------------------------------------------- 

    
        
        subroutine ReadDsilicaParameters

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: FromBlock

        !Local-----------------------------------------------------------------
        logical                                 :: BlockFound
        integer                                 :: ClientNumber
        integer                                 :: flag

        !----------------------------------------------------------------------

do1 :   do

            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,              &
                                        '<begin_dsi>', '<end_dsi>', BlockFound,     &
                                        STAT = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then    
cd2 :           if (BlockFound) then
                    
                    Me%BenthicCompute%DSilica=.true.
                    
                    call GetExtractType(FromBlock = FromBlock)
 

                    !ORGN - Stoichiometric equivalent between organic matter and phosphorus  
                    call GetData(Me%ORGSI,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'OM_STOICHIOMETRY_SI',                  &
                                 default       =  0.1800,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadDsilicaParameters - ModuleCEQUALW2 - ERR01"
                 
                    else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadDsilicaParameters - ModuleCEQUALW2 - ERR002"

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                stop "ReadDsilicaParameters - ModuleCEQUALW2 - ERR007"
       
            end if cd1
       
        end do do1


    end subroutine ReadDsilicaParameters
    
    !-------------------------------------------------------------------------- 

            
     subroutine ReadICarbonParameters

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: FromBlock

        !Local-----------------------------------------------------------------
        logical                                 :: BlockFound
        integer                                 :: ClientNumber
        integer                                 :: flag

        !----------------------------------------------------------------------

do1 :   do

            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,              &
                                        '<begin_ic>', '<end_ic>', BlockFound,     &
                                        STAT = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then    
cd2 :           if (BlockFound) then
                    
                    Me%BenthicCompute%ICarbon =.true.
                    Me%Compute%ICarbon        =.true.
                    
                    call GetExtractType(FromBlock = FromBlock)

                      select case (Me%Model)

                         case (BenthicCEQUALW2Model)
 

                    !ORGN - Stoichiometric equivalent between organic matter and phosphorus  
                    call GetData(Me%ORGC,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'OM_STOICHIOMETRY_C',                  &
                                 default       =  0.4500,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadICarbonParameters - ModuleCEQUALW2 - ERR01"
                  
                       end select

                    else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadICarbonParameters - ModuleCEQUALW2 - ERR002"

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                stop "ReadICarbonParameters - ModuleCEQUALW2 - ERR007"
       
            end if cd1
       
        end do do1


    end subroutine ReadICarbonParameters
    
    !-------------------------------------------------------------------------- 

  
    subroutine BenthicPropertyIndexNumber

               
        Me%Size%PropLB = 1
        Me%Size%PropUB = 0
                
  
  
        if (Me%BenthicCompute%Detritus) then
            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%Detritus      = Me%Size%PropUB
        endif

     
        if (Me%BenthicCompute%Ammonia) then
            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%Ammonia       = Me%Size%PropUB
        endif

        !Phosphorus index number
        if (Me%BenthicCompute%Phosphorus) then   
            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%Phosphorus    = Me%Size%PropUB
        endif   


        if (Me%BenthicCompute%DSilica) then      
            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%sidiss        = Me%Size%PropUB
        endif

        
        if (Me%BenthicCompute%ICarbon) then
            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%ICarbon       = Me%Size%PropUB
        endif

        if (Me%BenthicCompute%Oxygen) then
            Me%Size%PropUB             = Me%Size%PropUB + 1
            Me%PropIndex%Oxygen        = Me%Size%PropUB
        endif

        !----------------------------------------------------------------------

    end subroutine BenthicPropertyIndexNumber


    !--------------------------------------------------------------------------

        subroutine ConstructBenthicPropertyList

    ! -------------------------------------------------------------------------       
        
        allocate(Me%PropertyList(Me%Size%PropLB: Me%Size%PropUB))
        
        if (Me%BenthicCompute%Detritus) then
            Me%PropertyList(Me%PropIndex%Detritus)     = Detritus_
        endif

        
        if (Me%BenthicCompute%Ammonia) then
            Me%PropertyList(Me%PropIndex%Ammonia)     = Ammonia_
        endif

        !Phosphorus
        if (Me%BenthicCompute%Phosphorus) then   
            Me%PropertyList(Me%PropIndex%Phosphorus)  = Inorganic_Phosphorus_
        endif   

        if (Me%BenthicCompute%DSilica) then
            Me%PropertyList(Me%PropIndex%sidiss)      = DSilica_
        endif

        if (Me%BenthicCompute%ICarbon) then
            Me%PropertyList(Me%PropIndex%ICarbon)     = ICarbon_
        endif

        if (Me%BenthicCompute%Oxygen) then
            Me%PropertyList(Me%PropIndex%Oxygen)      = Oxygen_
        endif
  
        !----------------------------------------------------------------------

    end subroutine ConstructBenthicPropertyList


    
    subroutine ReadWaterColumnData(FileName)

        !Arguments-------------------------------------------------------------
        character(len=*)                     :: FileName

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWaterColumnData - ModuleCEQUALW2 - ERR01'


        call ReadGlobalVariables    

        call ConstructAlgaeClasses       

        call ConstructEpiphytonClasses   

        call ReadOMParameters            

        call ReadNitrogenParameters  
        
        call ReadPhosphorusParameters    

        call ReadICarbonParameters

        call ReadBODParameters 

        call ReadDetritusParameters
        
        call ReadOxygenParameters        

        call ReadSilicaParameters
        
        call PropertyIndexNumber
        
        call RateIndexNumber   

        call ConstructPropertyList 

        
        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWaterColumnData - ModuleCEQUALW2 - ERR03'

    end subroutine ReadWaterColumnData

    !--------------------------------------------------------------------------

    
    subroutine ReadGlobalVariables

        
        !External--------------------------------------------------------------
        integer                                     :: FromFile
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: flag

        !--------------------------------------------------------------------------

        call GetExtractType(FromFile = FromFile)

        if (Me%DTSecond .LE. 0.0) then

            !DTSecond, time step, in seconds, between two CEQUALW2 calls 
            call GetData(Me%DTSecond,                                       &
                         Me%ObjEnterData, flag,                             &
                         SearchType     = FromFile,                         &
                         keyword        ='DTSECONDS',                       & 
                         Default        = 60.0 * 60.0,                      &
                         ClientModule   = MohidModules(mCEQUALW2_)%Name,    &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'ReadGlobalVariables - ModuleCEQUALW2 - ERR01' 

            if (flag .EQ. 0) then
                write(*,*) 
                write(*,*) 'Keyword DTSECONDS not found in CEQUALW2 data file.'
                write(*,*) 'ReadGlobalVariables - ModuleCEQUALW2 - WRN01'
                write(*,*) 'Assumed ', Me%DTSecond, 'seconds (',  Me%DTSecond / 3600.0, 'hour).'
                write(*,*) 
            end if
        end if


        !For compatibility with the rest of the program
        Me%DTDay = Me%DTSecond / 24.0 / 60.0 / 60.0

        call GetData(Me%Compute%ICarbon,                                    &
                         Me%ObjEnterData, flag,                             &
                         SearchType     = FromFile,                         &
                         keyword        ='COMPUTE_ICARBON',                 & 
                         Default        = OFF,                              &
                         ClientModule   = MohidModules(mCEQUALW2_)%Name,    &
                         STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ReadGlobalVariables - ModuleCEQUALW2 - ERR04' 
        

    end subroutine ReadGlobalVariables

    
    !----------------------------------------------------------------------
    
    
    subroutine ConstructAlgaeClasses 

                
        !Local-----------------------------------------------------------------
        type (T_Algae),         pointer           :: NewAlgae
        integer                                   :: ClientNumber, STAT_CALL
        logical                                   :: BlockFound

        !Begin-----------------------------------------------------------------


do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<begin_algae>',      &
                                        block_end       = '<end_algae>',        &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then        

                    Me%Compute%Algae = .true.                                          
                    
                    call AddAlgae               (NewAlgae)

                    call ReadAlgaeParameters    (NewAlgae)

                    nullify(NewAlgae)

                else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructAlgaeClasses - ModuleCEQUALW2 - ERR01'
                        
                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'ConstructAlgaeClasses - ModuleCEQUALW2 - ERR02'
            else cd1
                stop       'ConstructAlgaeClasses - ModuleCEQUALW2 - ERR03'
            end if cd1
        end do do1

    end subroutine ConstructAlgaeClasses

    !--------------------------------------------------------------------------

    subroutine AddAlgae (ObjAlgae)

        !Arguments-------------------------------------------------------------
        type (T_Algae),      pointer           :: ObjAlgae

        !Local-----------------------------------------------------------------
        type (T_Algae),      pointer           :: PreviousAlgae
        type (T_Algae),      pointer           :: NewAlgae

        !Allocates new Algae
        allocate (NewAlgae)
        nullify  (NewAlgae%Next)

        !Insert new Algae into list and makes current algae point to it
        if (.not. associated(Me%FirstAlgae)) then
            Me%FirstAlgae       => NewAlgae
            ObjAlgae                 => NewAlgae
        else
            PreviousAlgae            => Me%FirstAlgae
            ObjAlgae                 => Me%FirstAlgae%Next

            do while (associated(ObjAlgae))
                PreviousAlgae        => ObjAlgae
                ObjAlgae             => ObjAlgae%Next
            enddo
            ObjAlgae                 => NewAlgae
            PreviousAlgae%Next       => NewAlgae
        endif

    end subroutine AddAlgae

    !--------------------------------------------------------------------------

    subroutine ReadAlgaeParameters(NewAlgae)

        !Arguments--------------------------------------------------------------
        type (T_Algae),   pointer               :: NewAlgae

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: FromBlock
        
        !Local-----------------------------------------------------------------
        integer                                 :: flag
        integer                                 :: ArrayLB, ArrayUB

        !----------------------------------------------------------------------
        
        ArrayLB = Me%Size%ArrayLB
        ArrayUB = Me%Size%ArrayUB

        call GetExtractType(FromBlock = FromBlock)

        call GetData(NewAlgae%ID%Name,                                  &
                     Me%ObjEnterData, flag,                             &
                     SearchType     = FromBlock,                        &
                     keyword        ='NAME',                            &
                     ClientModule   = MohidModules(mCEQUALW2_)%Name,    &
                     STAT           = STAT_CALL)
        if (STAT_CALL .ne. SUCCESS_) stop 'ReadAlgaeParameters - ModuleCEQUALW2 - ERR01'

        if (flag==0) then
            write (*,*)'Property without name'
            stop 'ReadAlgaeParameters - ModuleCEQUALW2 - ERR02'
        endif

        if (.not. CheckPropertyName (NewAlgae%ID%Name, NewAlgae%ID%IDnumber))then      
            write (*,*)'The following property isnt recognized by the model :'
            write (*,*)trim(NewAlgae%ID%Name)
            stop 'ReadAlgaeParameters - ModuleCEQUALW2- ERR03'     
        endif

        !AG - Algal maximum growth rate [day^-1]
        call GetData(NewAlgae%AG,                                                       &
                     Me%ObjEnterData, flag,                                             &
                     SearchType = FromBlock,                                            &
                     keyword        = 'A_GROWTH',                                       &
                     default        = 2.0000,                                           &
                     ClientModule   = MohidModules(mCEQUALW2_)%Name,                    &
                     STAT           = STAT_CALL)
        if (STAT_CALL .ne. SUCCESS_) stop 'ReadAlgaeParameters - ModuleCEQUALW2 - ERR04'

        !AR - Algal maximum respiration rate [day^-1]
        call GetData(NewAlgae%AR,                                                       &
                     Me%ObjEnterData, flag,                                             &
                     SearchType = FromBlock,                                            &
                     keyword        = 'A_RESPIRATION',                                  &
                     default        = 0.0400,                                           &
                     ClientModule   = MohidModules(mCEQUALW2_)%Name,                    &
                     STAT           = STAT_CALL)
        if (STAT_CALL .ne. SUCCESS_) stop 'ReadAlgaeParameters - ModuleCEQUALW2 - ERR05'

        !AE - Algal maximum excretion rate [day^-1]
        call GetData(NewAlgae%AE,                                                       &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_EXCRETION',                                     &
                     default       =  0.0400,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR06"

        !AM - Algal maximum mortality rate [day^-1]
        call GetData(NewAlgae%AM,                                                       &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_MORTALITY',                                     &
                     default       =  0.1000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR07"


        !AHSP - Algal half-saturation for phosphorus limited growth [g m^-3]
        call GetData(NewAlgae%AHSP,                                                     &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_HALFSAT_P',                                     &
                     default       =  0.0030,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR08"



        !AHSN - Algal half-saturation for nitrogen limited growth [g m^-3]
        call GetData(NewAlgae%AHSN,                                                     &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_HALFSAT_N',                                     &
                     default       =  0.0140,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR09"



        !AHSSI - Algal half-saturation for silica limited growth [g m^-3]
        call GetData(NewAlgae%AHSSI,                                                    &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_HALFSAT_SI',                                    &
                     default       =  0.0000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR010"


        !ASAT - Algal light saturation intensity at maximum phtosynthetic rate [W m^-2]
        call GetData(NewAlgae%ASAT,                                                     &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_LIGHT_SAT',                                     &
                     default       =  75.0000,                                          &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR011"

        
        
        !Algal temperature rate multipliers

        !AT1 - Lower temperature for algal growth (ºC)
        call GetData(NewAlgae%AT1,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_T1',                                            &
                     default       =  5.0000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR012"

        !AT2 - Lower temperature for maximum algal growth (ºC)
        call GetData(NewAlgae%AT2,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_T2',                                            &
                     default       =  25.0000,                                          &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR013"

        !AT3 - Upper temperature for maximum algal growth (ºC)
        call GetData(NewAlgae%AT3,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_T3',                                            &
                     default       =  35.0000,                                          &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR014"

        !AT4 - Upper temperature for algal growth (ºC)
        call GetData(NewAlgae%AT4,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_T4',                                            &
                     default       =  40.0000,                                          &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR015"

        !AK1 - Fraction of algal growth rate at AT1
        call GetData(NewAlgae%AK1,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_K1',                                            &
                     default       =  0.1000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR016"

        !AK2 - Fraction of maximum algal growth rate at AT2
        call GetData(NewAlgae%AK2,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_K2',                                            &
                     default       =  0.9900,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR017"

        
        !AK3 - Fraction of maximum algal growth rate at AT3
        call GetData(NewAlgae%AK3,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_K3',                                            &
                     default       =  0.9900,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR018"

        !AK4 - Fraction of algal growth rate at AT4
        call GetData(NewAlgae%AK4,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_K4',                                            &
                     default       =  0.1000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR019"

        
        !AP - Stoichiometric equivalent between algal biomass and phosphorus
        call GetData(NewAlgae%AP,                                                       &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_STOICHIOMETRY_P',                               &
                     default       =  0.0050,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR020"


        !AN - toichiometric equivalent between algal biomass and nitrogen  
        call GetData(NewAlgae%AN,                                                       &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_STOICHIOMETRY_N',                               &
                     default       =  0.0800,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR021"


        !AC - Stoichiometric equivalent between algal biomass and carbon
        call GetData(NewAlgae%AC,                                                       &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_STOICHIOMETRY_C',                               &
                     default       =  0.4500,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR022"


        !ASI - Stoichiometric equivalent between algal biomass and silica
        call GetData(NewAlgae%ASI,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_STOICHIOMETRY_Si',                              &
                     default       =  0.1800,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR023"


        !APOM - Fraction of algal biomass that is not converted to POM when algae die
        call GetData(NewAlgae%APOM,                                                     &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_POM',                                           &
                     default       =  0.8000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR025"

        !ANEQN - Equation number for algal ammonium preference (either 1 or 2)
        call GetData(NewAlgae%ANEQN,                                                    &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_NEQUATIONNUMBER',                                           &
                     default       =  2,                                                &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR026"
        if ((NewAlgae%ANEQN.NE.1).and.(NewAlgae%ANEQN.NE.2)) then
            write (*,*) "Possible values for equation number: 1 or 2!"
            stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR026A" 
        end if
         
            
        !ANPR - Algal half saturation preference constant for ammonium
        call GetData(NewAlgae%ANPR,                                                     &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_AMMONIUM_PREF',                                 &
                     default       =  0.0010,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR027"


        !O2AR - Oxygen stoichiometry for algal respiration
        call GetData(NewAlgae%O2AR,                                                     &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'O2_A_RESPIRATION',                                &
                     default       =  1.1000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR028"


        !02AG - Oxygen stoichiometry for algal primary prodution
        call GetData(NewAlgae%O2AG,                                                     &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'O2_A_GROWTH',                                     &
                     default       =  1.4000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR029"


        !EXA - Algal light extinction [m^-1]
        call GetData(NewAlgae%EXA,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'A_LIGHT_EXTINTION',                               &
                     default       =  0.2,                                              &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR030"



        allocate(NewAlgae%NLim(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR031"

        allocate(NewAlgae%PLim(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR032"

        allocate(NewAlgae%SLim(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR033"

        allocate(NewAlgae%LightLim(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR034"


        allocate(NewAlgae%OverallLim(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadAlgaeParameters - ModuleCEQUALW2 - ERR036"

    end subroutine ReadAlgaeParameters

    !----------------------------------------------------------------------------
 
    subroutine ConstructEpiphytonClasses 

        
        !Local-----------------------------------------------------------------
        type (T_Epiphyton),     pointer           :: NewEpiphyton
        integer                                   :: ClientNumber, STAT_CALL
        logical                                   :: BlockFound

        !Begin-----------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<begin_epiphyton>',  &
                                        block_end       = '<end_epiphyton>',    &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                    
                    Me%Compute%Epiphyton=.true.

                    call AddEpiphyton               (NewEpiphyton)

                    call ReadEpiphytonParameters    (NewEpiphyton)

                    nullify(NewEpiphyton)

                else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEpiphytonClasses - ModuleCEQUALW2 - ERR01'
                        
                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'ConstructEpiphytonClasses - ModuleCEQUALW2 - ERR02'
            else cd1
                stop       'ConstructEpiphytonClasses - ModuleCEQUALW2 - ERR03'
            end if cd1
        end do do1

    end subroutine ConstructEpiphytonClasses

    !----------------------------------------------------------------------------

    subroutine AddEpiphyton (ObjEpiphyton)

        !Arguments-------------------------------------------------------------
        type (T_Epiphyton),      pointer           :: ObjEpiphyton

        !Local-----------------------------------------------------------------
        type (T_Epiphyton),      pointer           :: PreviousEpiphyton
        type (T_Epiphyton),      pointer           :: NewEpiphyton

        !Allocates new Epiphyton
        allocate (NewEpiphyton)
        nullify  (NewEpiphyton%Next)

        !Insert new Epiphyton into list and makes current algae point to it
        if (.not. associated(Me%FirstEpiphyton)) then
            Me%FirstEpiphyton   => NewEpiphyton
            ObjEpiphyton                 => NewEpiphyton
        else
            PreviousEpiphyton            => Me%FirstEpiphyton
            ObjEpiphyton                 => Me%FirstEpiphyton%Next

            do while (associated(ObjEpiphyton))
                PreviousEpiphyton        => ObjEpiphyton
                ObjEpiphyton             => ObjEpiphyton%Next
            enddo
            ObjEpiphyton                 => NewEpiphyton
            PreviousEpiphyton%Next       => NewEpiphyton
        endif

    end subroutine AddEpiphyton

    
    !--------------------------------------------------------------------------
    
    
    subroutine ReadEpiphytonParameters (NewEpiphyton)

        !Arguments--------------------------------------------------------------
        type (T_Epiphyton),     pointer             :: NewEpiphyton

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: FromBlock
        !Local-----------------------------------------------------------------
        integer                                     :: flag
        integer                                     :: ArrayLB, ArrayUB

        !----------------------------------------------------------------------
        
        ArrayLB = Me%Size%ArrayLB
        ArrayUB = Me%Size%ArrayUB

        call GetExtractType(FromBlock = FromBlock)

        call GetData(NewEpiphyton%ID%Name,                                              &
                     Me%ObjEnterData, flag,                                             &
                     FromBlock,                                                         &
                     keyword        ='NAME',                                            &
                     ClientModule   = MohidModules(mCEQUALW2_)%Name,                    &
                     STAT           = STAT_CALL)
        if (STAT_CALL.ne.SUCCESS_) stop  'ReadEpiphytonParameters - ModuleCEQUALW2 - ERR01'
        if (flag==0) then
            write (*,*)'Property without name'
            stop  'ReadEpiphytonParameters - ModuleCEQUALW2 - ERR02'
        endif


        if (.not. CheckPropertyName (NewEpiphyton%ID%Name,NewEpiphyton%ID%IDNumber)) then      
            write (*,*)'The following property isnt recognized by the model :'
            write (*,*)trim(NewEpiphyton%ID%Name)
            stop  'ReadEpiphytonParameters - ModuleCEQUALW2 - ERR03'     
        endif

        !Default data from W2Con.npt

        !EG - Maximum epiphyton growth rate [day^-1]
        call GetData(NewEpiphyton%EG,                                                   &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_GROWTH',                                        &
                     default       =  2.0000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR04"


        !ER - Maximum epiphyton respiration rate [day^-1]
        call GetData(NewEpiphyton%ER,                                                   &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_RESPIRATION',                                   &
                     default       =  0.0400,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR05"

        
        !EE - Maximum epiphyton excretion rate [day^-1]        
        call GetData(NewEpiphyton%EE,                                                   &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_EXCRETION',                                     &
                     default       =  0.0400,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR06"

        
        !EM - Maximum epiphyton mortality rate [day^-1]        
        call GetData(NewEpiphyton%EM,                                                   &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_MORTALITY',                                     &
                     default       =  0.1000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR07"


        !EHSP - Epiphyton half-saturation for phosphorus limited growth [g m^-3]
        call GetData(NewEpiphyton%EHSP,                                                 &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_HALFSAT_P',                                     &
                     default       =  0.0030,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR08"


        !EHSN - Epiphyton half-saturation for nitrogen limited growth [g m^-3]
        call GetData(NewEpiphyton%EHSN,                                                 &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_HALFSAT_N',                                     &
                     default       =  0.0140,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR09"


        !EHSSI - Epiphyton half-saturation for silica limited growth [g m^-3]
        call GetData(NewEpiphyton%EHSSI,                                                &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_HALFSAT_SI',                                    &
                     default       =  0.0000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR010"


        !ESAT - Light saturation intensity at maximum photosynthetic rate  [W m^-2]
        call GetData(NewEpiphyton%ESAT,                                                 &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_LIGHT_SAT',                                     &
                     default       =  75.0000,                                          &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR011"

        
        !ENEQN - Ammonia preference factor equation for epiphyton (either 1 or 2)
        call GetData(NewEpiphyton%ENEQN,                                                &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_NEQUATIONNUMBER',                               &
                     default       =  2,                                                &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR013"
        if ((NewEpiphyton%ENEQN.NE.1).and.(NewEpiphyton%ENEQN.NE.2)) then
            write (*,*) "Possible values for equation number: 1 or 2!"
            stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR012A" 
        end if

        !ENPR - N preference half-saturation constant (only used if EPEQN = 2) [mg/l]
        call GetData(NewEpiphyton%ENPR,                                                 &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_AMMONIUM_PREF',                                 &
                     default       =  0.0010,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR014"

       
       !ET1 - Lower temperature for epiphyton growth [ºC]
        call GetData(NewEpiphyton%ET1,                                                  &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_T1',                                            &
                     default       =  5.0000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR015"


       !ET2 - Lower temperature for maximum epiphyton growth [ºC]
        call GetData(NewEpiphyton%ET2,                                                  &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_T2',                                            &
                     default       =  25.0000,                                          &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR016"


       !ET3 - Upper temperature for maximum epiphyton growth [ºC]
        call GetData(NewEpiphyton%ET3,                                                  &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_T3',                                            &
                     default       =  35.0000,                                          &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR017"

        
       !ET4 - Upper temperature for epiphyton growth [ºC]        
        call GetData(NewEpiphyton%ET4,                                                  &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_T4',                                            &
                     default       =  40.0000,                                          &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR018"

       
       !EK1 - Fraction of epiphyton growth rate at ET1
        call GetData(NewEpiphyton%EK1,                                                  &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_K1',                                            &
                     default       =  0.1000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR019"


       !EK2 - Fraction of maximum epiphyton growth rate at ET2
        call GetData(NewEpiphyton%EK2,                                                  &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_K2',                                            &
                     default       =  0.9900,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR020"

       !EK3 - Fraction of maximum epiphyton growth rate at ET3
        call GetData(NewEpiphyton%EK3,                                                  &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_K3',                                            &
                     default       =  0.9900,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR021"


       !EK4 - Fraction of epiphyton growth rate at ET4
        call GetData(NewEpiphyton%EK4,                                                  &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_K4',                                            &
                     default       =  0.1000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR022"

        !EP - Stoichiometric equivalent between epiphyton biomass and phosphorus
        call GetData(NewEpiphyton%EP,                                                   &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_STOICHIOMETRY_P',                               &
                     default       =  0.0050,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR023"

        !EN - Stoichiometric equivalent between epiphyton biomass and nitrogen
        call GetData(NewEpiphyton%EN,                                                   &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_STOICHIOMETRY_N',                               &
                     default       =  0.0800,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR024"


        !EC - Stoichiometric equivalent between epiphyton biomass and carbon
        call GetData(NewEpiphyton%EC,                                                   &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_STOICHIOMETRY_C',                               &
                     default       =  0.4500,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR025"


        !ESI - Stoichiometric equivalent between epiphyton biomass and silica
        call GetData(NewEpiphyton%ESI,                                                  &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_STOICHIOMETRY_SI',                              &
                     default       =  0.1800,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR026"


        !Fraction of epiphyton biomass that is not converted to POM when epiphyton die
        call GetData(NewEpiphyton%EPOM,                                                 &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'E_POM',                                           &
                     default       =  0.8000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR028"

        !O2ER - Oxygen stoichiometry for epiphyton respiration 
        call GetData(NewEpiphyton%O2ER,                                                 &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'O2_E_RESPIRATION',                                &
                     default       =  1.1000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR029"

        !O2EG - Oxygen stoichiometry for epiphyton primary prodution
        call GetData(NewEpiphyton%O2EG,                                                 &
                     Me%ObjEnterData, flag,                                             &
                     SearchType    =  FromBlock,                                        &
                     keyword       = 'O2_E_GROWTH',                                     &
                     default       =  1.4000,                                           &
                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,                    &
                     STAT          =  STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR030"


        allocate(NewEpiphyton%NLim(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR031"

        allocate(NewEpiphyton%PLim(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR032"

        allocate(NewEpiphyton%SLim(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR033"

        allocate(NewEpiphyton%LightLim(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR034"


        allocate(NewEpiphyton%OverallLim(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop "ReadEpiphytonParameters - ModuleCEQUALW2 - ERR036"


    end subroutine ReadEpiphytonParameters

   !----------------------------------------------------------------------------
    
    subroutine ReadOMParameters 

       
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        logical                                 :: BlockFound
        integer                                 :: ClientNumber
        integer                                 :: flag
        integer                                 :: FromBlock

        !----------------------------------------------------------------------

do1 :   do

            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,            &
                                        '<begin_om>', '<end_om>', BlockFound,     &
                                        STAT = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_) then    
cd2 :           if (BlockFound) then
                    
                    Me%Compute%OrganicMatter=.true.

                    call GetExtractType(FromBlock = FromBlock)

   
                    !DOM Parameters---------------------------------------------  
                    
                    
                    !LDOMDK - Labile DOM decay rate [day^-1]                      
                    call GetData(Me%LDOMDK,                                               &
                                 Me%ObjEnterData, flag,                                   &
                                 SearchType    =  FromBlock,                              &
                                 keyword       = 'LDOM_DECAY',                            &
                                 default       =  0.1000,                                 &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,          &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR01"


                    !RDOMDK - Refractory DOM decay rate [day^-1] 
                    call GetData(Me%RDOMDK,                                             &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'RDOM_DECAY',                          &
                                 default       =  0.0010,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR02"


                    !LRDDK - Labile to refractory DOM decay rate [day^-1]
                    call GetData(Me%LRDDK,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'LRDOM_DECAY',                         &
                                 default       =  0.0100,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR03"

                    
                    !LPOMDK - Labile POM decay rate [day^-1]
                    call GetData(Me%LPOMDK,                                             &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'LPOM_DECAY',                          &
                                 default       =  0.0800,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR04"

                    
                    !RPOMDK - Refractory POM decay rate [day^-1]
                    call GetData(Me%RPOMDK,                                             &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'RPOM_DECAY',                          &
                                 default       =  0.0010,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR05"

                    
                    !LRPDK - Labile to refractory POM decay rate [day^-1]
                    call GetData(Me%LRPDK,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'LRPOM_DECAY',                         &
                                 default       =  0.0100,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR06"

                    
                    !ORGP - Stoichiometric equivalent between organic matter and phosphorus
                    call GetData(Me%ORGP,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'OM_STOICHIOMETRY_P',                  &
                                 default       =  0.0050,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR07"

  
                    !ORGN - Stoichiometric equivalent between organic matter and nitrogen  
                    call GetData(Me%ORGN,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'OM_STOICHIOMETRY_N',                  &
                                 default       =  0.0800,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR08"


                    !ORGC - Stoichiometric equivalent between organic matter and carbon
                    call GetData(Me%ORGC,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'OM_STOICHIOMETRY_C',                  &
                                 default       =  0.4500,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR09"

 
                    !ORGSI - Stoichiometric equivalent between organic matter and silica
                    call GetData(Me%ORGSI,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'OM_STOICHIOMETRY_SI',                 &
                                 default       =  0.1800,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR010"


                    !OMT1 - Lower temperature for organic matter decay [ºC]
                    call GetData(Me%OMT1,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'OM_T1',                               &
                                 default       =  4.0000,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR011"

                   !OMT2 - Upper temperature for organic matter decay [ºC]
                    call GetData(Me%OMT2,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'OM_T2',                               &
                                 default       =  25.0000,                              &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR012"

  
                   !OMK1 - Fraction of organic matter decay rate at OMT1
                    call GetData(Me%OMK1,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'OM_K1',                               &
                                 default       =  0.1000,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR013"

                    
                    !OMK2 - Fraction of organic matter decay rate at OMT2
                    call GetData(Me%OMK2,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'OM_K2',                               &
                                 default       =  0.9900,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR014"

                                     
                   
                    !O2OM - Oxygen stoichiometry for organic matter decay
                    call GetData(Me%O2OM,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'O2_OM',                               &
                                 default       =  1.4000,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR015"



                else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOMParameters - ModuleCEQUALW2 - ERR016"

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                stop "ReadOMParameters - ModuleCEQUALW2 - ERR017"
       
            end if cd1
       
        end do do1

    end subroutine ReadOMParameters   
   
    
    !----------------------------------------------------------------------------
    
    subroutine ReadBODParameters 

        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: FromBlock

        !Local-----------------------------------------------------------------
        logical                                 :: BlockFound
        integer                                 :: ClientNumber
        integer                                 :: flag

        !----------------------------------------------------------------------

do1 :   do

            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,              &
                                        '<begin_BOD>', '<end_BOD>', BlockFound,     &
                                        STAT = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then    
cd2 :           if (BlockFound) then
                    
                    Me%Compute%BOD=.true.
                    call GetExtractType(FromBlock = FromBlock)
 

                    !KBOD - 5-day decay rate at 20ºC [day^-1]
                    call GetData(Me%KBOD,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'BOD_DECAY',                           &
                                 default       =  0.25,                                 &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadBODParameters - ModuleCEQUALW2 - ERR01"


                    !TBOD - Temperature coefficient
                    call GetData(Me%TBOD,                                                &
                                 Me%ObjEnterData, flag,                                  &
                                 SearchType    =  FromBlock,                             &
                                 keyword       = 'BOD_T_COEF',                           &
                                 default       =  1.0147,                                &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,         &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadBODParameters - ModuleCEQUALW2 - ERR02"

                    
                    !RBOD - Ratio of CBOD5 to ultimate CBOD
                    call GetData(Me%RBOD,                                                 &
                                 Me%ObjEnterData, flag,                                   &
                                 SearchType    =  FromBlock,                              &
                                 keyword       = 'BOD_RATIO',                             &
                                 default       =  1.8500,                                 &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,          &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadBODParameters - ModuleCEQUALW2 - ERR03"

                    
                    !BODP - P stoichiometry for CBOD decay
                    call GetData(Me%BODP,                                                 &
                                 Me%ObjEnterData, flag,                                   &
                                 SearchType    =  FromBlock,                              &
                                 keyword       = 'BOD_STOICHIOMETRY_P',                   &
                                 default       =  0.0040,                                 &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,          &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadBODParameters - ModuleCEQUALW2 - ERR04"

                    
                    !BODN - N stoichiometry for CBOD decay
                    call GetData(Me%BODN,                                                 &
                                 Me%ObjEnterData, flag,                                   &
                                 SearchType    =  FromBlock,                              &
                                 keyword       = 'BOD_STOICHIOMETRY_N',                   &
                                 default       =  0.0600,                                 &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,          &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadBODParameters - ModuleCEQUALW2 - ERR05"


                    !BODC - C stoichiometry for CBOD decay
                    call GetData(Me%BODC,                                                 &
                                 Me%ObjEnterData, flag,                                   &
                                 SearchType    =  FromBlock,                              &
                                 keyword       = 'BOD_STOICHIOMETRY_C',                   &
                                 default       =  0.3200,                                 &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,          &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadBODParameters - ModuleCEQUALW2 - ERR06"

                          
                else cd2
                    
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadBODParameters - ModuleCEQUALW2 - ERR07"

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                stop "ReadBODParameters - ModuleCEQUALW2 - ERR08"
       
            end if cd1
       
        end do do1

    end subroutine ReadBODParameters   
   
    
    !----------------------------------------------------------------------------    
   
    
    subroutine ReadSilicaParameters 
        
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: FromBlock

        !Local-----------------------------------------------------------------
        
        character(len=StringLength), parameter :: block_begin      = '<beginsilica>'
        character(len=StringLength), parameter :: block_end        = '<endsilica>'
 
        logical :: BlockFound
        integer :: ClientNumber
        integer :: flag

        !----------------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                  &
                                        '<begin_silica>', '<end_silica>', BlockFound,   &
                                        STAT = STAT_CALL)

cd1 :       if (STAT_CALL .EQ. SUCCESS_)then    
cd2 :           if (BlockFound) then  
                    
                    Me%Compute%Silica=.true.
                    call GetExtractType(FromBlock = FromBlock)


                   !Particulate biogenic silica decay rate [day^-1]
                    call GetData(Me%PSIDK,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'PARTSI_DECAY',                        &
                                 default       =  0.3000,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadSilicaParameters - ModuleCEQUALW2 - ERR03"

                                       
                else cd2
                    
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadSilicaParameters - ModuleCEQUALW2 - ERR05"

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                stop "ReadSilicaParameters - ModuleCEQUALW2 - ERR06"
       
            end if cd1
       
        end do do1

    end subroutine ReadSilicaParameters   
   
    
    !----------------------------------------------------------------------------

    
    
    

    subroutine ReadOxygenParameters 
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: FromBlock

        !Local-----------------------------------------------------------------
        logical                                 :: BlockFound
        integer                                 :: ClientNumber
        integer                                 :: flag

        !----------------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,     &
                                        '<begin_oxygen>', '<end_oxygen>', BlockFound, &
                                        STAT = STAT_CALL)

cd1 :       if(STAT_CALL .EQ. SUCCESS_) then    
cd2 :           if (BlockFound) then
                        
                        call GetExtractType(FromBlock = FromBlock)
                       
                        Me%BenthicCompute%Oxygen =.True.

                        !O2LIM - Dissolved oxygen concentration at which anaerobic processes begin [g m^-3]
                        call GetData(Me%O2LIM,                                                  &
                                     Me%ObjEnterData, flag,                                     &
                                     SearchType    =  FromBlock,                                &
                                     keyword       = 'O2LIM',                                   &
                                     default       =  0.1000,                                   &
                                     ClientModule  =  MohidModules(mCEQUALW2_)%Name,            &
                                     STAT          =  STAT_CALL)
                        if(STAT_CALL .ne. SUCCESS_) stop "ReadOxygenParameters - ModuleCEQUALW2 - ERR01"
                         
                else cd2
                    
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadOxygenParameters - ModuleCEQUALW2 - ERR02"

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                stop "ReadOxygenParameters - ModuleCEQUALW2 - ERR03"
       
            end if cd1
       
        end do do1

        !-----------------------------------------------------------------------

    end subroutine ReadOxygenParameters   

   


    
    !----------------------------------------------------------------------------
   
    
    subroutine ReadNitrogenParameters 
        
       
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: FromBlock

        !Local-----------------------------------------------------------------
        logical                                 :: BlockFound
        integer                                 :: ClientNumber
        integer                                 :: flag

        !----------------------------------------------------------------------

do1 :   do

            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                  &
                                        '<begin_nitro>', '<end_nitro>', BlockFound,     &
                                        STAT = STAT_CALL)

cd1 :       if(STAT_CALL .EQ. SUCCESS_)then    
cd2 :           if (BlockFound) then      
                    
                    Me%Compute%Nitrogen=.true.

                    call GetExtractType(FromBlock = FromBlock)

                    call GetData(Me%NH4DK,                                                &
                                 Me%ObjEnterData, flag,                                   &
                                 SearchType    =  FromBlock,                              &
                                 keyword       = 'NH4_DECAY',                             &
                                 default       =  0.1200,                                 &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,          &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadNitrogenParameters - ModuleCEQUALW2 - ERR02"

                    
                    !NH4T1 - Lower temperature for ammonia decay [ºC]
                    call GetData(Me%NH4T1,                                                &
                                 Me%ObjEnterData, flag,                                   &
                                 SearchType    =  FromBlock,                              &
                                 keyword       = 'NH4_T1',                                &
                                 default       =  5.0000,                                 &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,          &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadNitrogenParameters - ModuleCEQUALW2 - ERR03"

  
                    !NH4T2 - Upper temperature for ammonia decay [ºC]  
                    call GetData(Me%NH4T2,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'NH4_T2',                              &
                                 default       =  25.0000,                              &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadNitrogenParameters - ModuleCEQUALW2 - ERR04"

                    
                    
                    !NH4K1 - Fraction of nitrification rate at NH4T1
                    call GetData(Me%NH4K1,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'NH4_K1',                              &
                                 default       =  0.1000,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadNitrogenParameters - ModuleCEQUALW2 - ERR05"

 
                    !NH4K2 - Fraction of nitrification rate at NH4T2 
                    call GetData(Me%NH4K2,                                                &
                                 Me%ObjEnterData, flag,                                   &
                                 SearchType    =  FromBlock,                              &
                                 keyword       = 'NH4_K2',                                &
                                 default       =  0.9900,                                 &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,          &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadNitrogenParameters - ModuleCEQUALW2 - ERR06"

 
                    !NO3DK - Nitrate decay rate [day^-1]
                    call GetData(Me%NO3DK,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'NO3_DECAY',                           &
                                 default       =  0.0300,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadNitrogenParameters - ModuleCEQUALW2 - ERR07"


                    !NO3T1 - Lower temperature for nitrate decay [ºC] 
                    call GetData(Me%NO3T1,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'NO3_T1',                              &
                                 default       =  5.0000,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadNitrogenParameters - ModuleCEQUALW2 - ERR09"


                    !NO3T2 - Upper temperature for nitrate decay [ºC]  
                    call GetData(Me%NO3T2,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'NO3_T2',                              &
                                 default       =  25.0000,                              &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadNitrogenParameters - ModuleCEQUALW2 - ERR010"

 
                    !NO3K1 - Fraction of denitrification rate at NO3T1 
                    call GetData(Me%NO3K1,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'NO3_K1',                              &
                                 default       =  0.1000,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadNitrogenParameters - ModuleCEQUALW2 - ERR011"


                    !NO3K2 - Fraction of denitrification rate at NO3T2
                    call GetData(Me%NO3K2,                                              &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType    =  FromBlock,                            &
                                 keyword       = 'NO3_K2',                              &
                                 default       =  0.9900,                               &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,        &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadNitrogenParameters - ModuleCEQUALW2 - ERR012"

 
                    !O2NH4 - Oxygen stoichiometry for nitrification
                    call GetData(Me%O2NH4,                                                &
                                 Me%ObjEnterData, flag,                                   &
                                 SearchType    =  FromBlock,                              &
                                 keyword       = 'O2_NH4',                                &
                                 default       =  4.5700,                                 &
                                 ClientModule  =  MohidModules(mCEQUALW2_)%Name,          &
                                 STAT          =  STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadNitrogenParameters - ModuleCEQUALW2 - ERR013"


                else cd2

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop "ReadNitrogenParameters - ModuleCEQUALW2 - ERR14"

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                stop "ReadNitrogenParameters - ModuleCEQUALW2 - ERR15"
       
            end if cd1
       
        end do do1

    end subroutine ReadNitrogenParameters   
   
    !----------------------------------------------------------------------------

    subroutine ConstructPropertyList

       
        !Local-----------------------------------------------------------------
        type(T_Algae    ),      pointer             :: Algae
        type(T_Epiphyton),      pointer             :: Epiphyton
        integer                                     :: Index
        !Local-----------------------------------------------------------------
        
        allocate(Me%PropertyList(Me%Size%PropLB: Me%Size%PropUB))
        
        Index = 0

        !Producer index number      
        Algae => Me%FirstAlgae
        do while(associated(Algae))
            
            Me%PropertyList(Algae%PropIndex)     = Algae%ID%IDNumber

            Algae => Algae%Next
        end do
        
        Epiphyton => Me%FirstEpiphyton
        do while(associated(Epiphyton))
            
            Me%PropertyList(Epiphyton%PropIndex) = Epiphyton%ID%IDNumber

            Epiphyton => Epiphyton%Next
        end do


        !Nitrogen 
        if (Me%Compute%Nitrogen) then
            Me%PropertyList(Me%PropIndex%Ammonia)     = Ammonia_
            Me%PropertyList(Me%PropIndex%Nitrate)     = Nitrate_
        endif

        !Phosphorus
        if (Me%Compute%Phosphorus) then   
            Me%PropertyList(Me%PropIndex%Phosphorus)  = Inorganic_Phosphorus_
        endif   

        !OrganicMatter
        if (Me%Compute%OrganicMatter) then
            Me%PropertyList(Me%PropIndex%pomref)      = RPOM_
            Me%PropertyList(Me%PropIndex%pomlab)      = LPOM_
            Me%PropertyList(Me%PropIndex%domlab)      = LDOM_
            Me%PropertyList(Me%PropIndex%domref)      = RDOM_
        endif

        if (Me%Compute%Silica) then
            Me%PropertyList(Me%PropIndex%sipart)      = PSilica_
            Me%PropertyList(Me%PropIndex%sidiss)      = DSilica_
        endif

     

        if (Me%Compute%ICarbon) then
            Me%PropertyList(Me%PropIndex%ICarbon)     = ICarbon_
            Me%PropertyList(Me%PropIndex%CO2)         = CarbonDioxide_
            Me%PropertyList(Me%PropIndex%pH)          = pH_
            Me%PropertyList(Me%PropIndex%HCO3)        = HCO3_
            Me%PropertyList(Me%PropIndex%CO3)         = CO3_
        endif

        !Oxygen
        Me%PropertyList(Me%PropIndex%Oxygen)          = Oxygen_

        !BOD index number
        if (Me%Compute%BOD) then   
            Me%PropertyList(Me%PropIndex%BOD)         = BOD_
        endif

        if (Me%Compute%Detritus) then   
            Me%PropertyList(Me%PropIndex%Detritus)    = Detritus_
        endif

        !----------------------------------------------------------------------

    end subroutine ConstructPropertyList


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetCEQUALW2Size(CEQUALW2_ID, PropLB, PropUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: CEQUALW2_ID
        integer, optional, intent(OUT)      :: PropLB,PropUB
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
       
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(CEQUALW2_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PropLB   )) PropLB    = Me%Size%PropLB
            if (present(PropUB   )) PropUB    = Me%Size%PropUB

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetCEQUALW2Size

    
    !--------------------------------------------------------------------------

    
    subroutine GetCEQUALW2Options(CEQUALW2_ID,                         &
                                   Algae,                              &
                                   Epiphyton,                          & 
                                   Nitrogen,                           &
                                   Phosphorus,                         & 
                                   OrganicMatter,                      &
                                   ICarbon,                            &
                                   BOD,                                &
                                   Oxygen,                             &
                                   Silica,                             &
                                   Detritus,                           &
                                   STAT) 

        !Arguments-------------------------------------------------------------
        integer                             :: CEQUALW2_ID
        integer, optional, intent(OUT)      :: STAT
        logical, optional, intent(OUT)      :: Algae
        logical, optional, intent(OUT)      :: Nitrogen
        logical, optional, intent(OUT)      :: Phosphorus
        logical, optional, intent(OUT)      :: Oxygen
        logical, optional, intent(OUT)      :: BOD
        logical, optional, intent(OUT)      :: OrganicMatter
        logical, optional, intent(OUT)      :: Epiphyton
        logical, optional, intent(OUT)      :: ICarbon
        logical, optional, intent(OUT)      :: Silica 
        logical, optional, intent(OUT)      :: Detritus 

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_    
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(CEQUALW2_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            if (present(Algae          )) Algae          = Me%Compute%Algae    
            if (present(Epiphyton      )) Epiphyton      = Me%Compute%Epiphyton    
            if (present(Nitrogen       )) Nitrogen       = Me%Compute%Nitrogen    
            if (present(Phosphorus     )) Phosphorus     = Me%Compute%Phosphorus
            if (present(OrganicMatter  )) OrganicMatter  = Me%Compute%OrganicMatter
            if (present(ICarbon        )) Icarbon        = Me%Compute%ICarbon
            if (present(BOD            )) BOD            = Me%Compute%BOD
            if (present(Oxygen         )) Oxygen         = Me%Compute%Oxygen 
            if (present(Silica         )) Silica         = Me%Compute%Silica
            if (present(Detritus       )) Detritus       = Me%Compute%Detritus
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetCEQUALW2Options

    
    !--------------------------------------------------------------------------

    
    subroutine GetCEQUALW2PropIndex (CEQUALW2_ID, PropertyIDNumber, PropertyIndex, STAT)

                                     

        !Arguments-------------------------------------------------------------
        integer                             :: CEQUALW2_ID
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

        call Ready(CEQUALW2_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            found = .false.
            do CurrentIndex = Me%Size%PropLB,Me%Size%PropUB

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

    end subroutine GetCEQUALW2PropIndex


    !--------------------------------------------------------------------------

    
    subroutine GetDTCEQUALW2(CEQUALW2_ID, DTDay, DTSecond, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: CEQUALW2_ID
        real,    optional, intent(OUT)      :: DTDay
        real,    optional, intent(OUT)      :: DTSecond
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(CEQUALW2_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(DTDay   )) DTDay    = Me%DTDay
            if (present(DTSecond)) DTSecond = Me%DTSecond

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDTCEQUALW2
   
    !----------------------------------------------------------------------

    subroutine GetCEQUALW2PropertyList(CEQUALW2_ID, PropertyList, STAT)

        !Arguments-------------------------------------------------------------
        integer                                                 :: CEQUALW2_ID
        integer, dimension(:), pointer                          :: PropertyList
        integer, optional, intent(OUT)                          :: STAT

        !External--------------------------------------------------------------
        integer                                                 :: ready_              

        !Local-----------------------------------------------------------------
        integer                                                 :: STAT_
       
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(CEQUALW2_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mCEQUALW2_, Me%InstanceID)

            PropertyList =>  Me%PropertyList

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetCEQUALW2PropertyList

    
    !----------------------------------------------------------------------
    
    
    subroutine UnGetCEQUALW2_1D_Int(CEQUALW2_ID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                             :: CEQUALW2_ID
        integer, dimension(:), pointer:: Array
        integer, intent(OUT), optional                      :: STAT

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_, ready_
       
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(CEQUALW2_ID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mCEQUALW2_, Me%InstanceID, "UnGetCEQUALW2_1D_Int")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetCEQUALW2_1D_Int


    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine GetCEQUALW2RateFlux(CEQUALW2_ID, CequalRateIndex , CequalRateFlux, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: CEQUALW2_ID
        real, dimension(:),     pointer     :: CequalRateFlux
        integer                             :: CequalRateIndex
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer         :: ready_        

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
        logical                             :: found  
        integer                             :: nrate     

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(CEQUALW2_ID, ready_)    
        
        found=.false.

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


          call Read_Lock(mCEQUALW2_, Me%InstanceID)

           do nrate = Me%Rate%LB, Me%Rate%UB
             if (CequalRateIndex.eq.me%rate%match(nrate)) then
              found=.true.
               CequalRateFlux => Me%Rate%Value (nrate, :)
             endif
           enddo
           
           if (.not.found) &
             stop 'GetCEQUALW2RateFlux - ModuleCEQUALW2 - ERR01'
           

 
           STAT_ = SUCCESS_
           
        else cd1
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetCEQUALW2RateFlux

    !--------------------------------------------------------------------------

    
    
    !--------------------------------------------------------------------------
    subroutine UnGetCEQUALW2RateFlux(CEQUALW2_ID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: CEQUALW2_ID
        integer, optional, intent(OUT)  :: STAT
        real, pointer, dimension(:)     :: Array

        !External--------------------------------------------------------------
        integer                         :: ready_   

        !Local-----------------------------------------------------------------
        integer                         :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(CEQUALW2_ID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_Unlock(mCEQUALW2_, Me%InstanceID, "CEQUALW2RateFlux")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UnGetCEQUALW2RateFlux

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    
    subroutine CEQUALW2(CEQUALW2_ID,                                        &
                        Salinity,                                           &
                        Temperature,                                        &
                        Alkalinity,                                         &
                        ShortWaveRadiation,                                 &
                        LightExtCoefField,                                  &
                        Thickness,                                          &
                        Mass,                                               &
                        OpenPoints,                                         &
                        STAT)  

        !Arguments---------------------------------------------------------------
        integer                                       :: CEQUALW2_ID
        real,                 pointer, dimension(:  ) :: Salinity
        real,                 pointer, dimension(:  ) :: Temperature
        real,                 pointer, dimension(:  ) :: Alkalinity
        real,                 pointer, dimension(:  ) :: ShortWaveRadiation
        real,                 pointer, dimension(:  ) :: LightExtCoefField
        real,                 pointer, dimension(:  ) :: Thickness
        real,                 pointer, dimension(:,:) :: Mass
        integer, optional,    pointer, dimension(:  ) :: OpenPoints
        integer, optional,    intent(OUT)             :: STAT
         
        !External----------------------------------------------------------------
        integer                                       :: index
        integer                                       :: ready_   
               
        !Local-------------------------------------------------------------------
        integer                                       :: STAT_, propi          
        logical                                       :: CalcPoint
        !------------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(CEQUALW2_ID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            Me%ExternalVar%Salinity                   => Salinity
            if (.NOT. associated(Me%ExternalVar%Salinity))             &
                stop 'CEQUALW2 - ModuleCEQUALW2 - ERR01'


            Me%ExternalVar%Temperature                => Temperature
            if (.NOT. associated(Me%ExternalVar%Temperature))          &
                stop 'CEQUALW2 - ModuleCEQUALW2 - ERR02'


            Me%ExternalVar%Alkalinity                => Alkalinity
            if (.NOT. associated(Me%ExternalVar%Alkalinity))           &
                stop 'CEQUALW2 - ModuleCEQUALW2 - ERR03'

                
            Me%ExternalVar%ShortWaveRadiation         => ShortWaveRadiation
            if (.NOT. associated(Me%ExternalVar%ShortWaveRadiation))   &
                stop 'CEQUALW2 - ModuleCEQUALW2 - ERR04'
            
            Me%ExternalVar%LightExtCoefField          => LightExtCoefField
            if (.NOT. associated(Me%ExternalVar%LightExtCoefField))    &
                stop 'CEQUALW2 - ModuleCEQUALW2 - ERR05'
            
            Me%ExternalVar%Thickness                  => Thickness
            if (.NOT. associated(Me%ExternalVar%Thickness))            &
                stop 'CEQUALW2 - ModuleCEQUALW2 - ERR06'

            Me%ExternalVar%Mass                       => Mass
            if (.NOT. associated(Me%ExternalVar%Mass))                 &
                stop 'CEQUALW2 - ModuleCEQUALW2 - ERR07'
             
            allocate (Me%SinksSources(Me%Size%PropLB:Me%Size%PropUB)) 
                          
do1 :       do index = Me%Size%ArrayLB, Me%Size%ArrayUB
            
                !If this module is called from the CEQUALW23D module, OpenPoint is present
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


cd0:            if (CalcPoint) then

                    call TemperatureRateMultipliers     (index)

                    call ComputeKineticRates            (index)
                    
                    
                    

                   

                    !Phosphorus
                    if (Me%Compute%Phosphorus  )then
                        call ComputePhosphorus          (index)
                    end if       
                    
                    
                    !Nitrogen
                    if (Me%Compute%Nitrogen    )then
                        call ComputeAmmonia             (index)
                        call ComputeNitrate             (index)   
                    endif

                    !Silica
                    if (Me%Compute%Silica      )then
                        call ComputeDissolvedSilica     (index)
                        call ComputeParticulateSilica   (index)
                    endif
                    
                    !OrganicMatter
                    if (Me%Compute%OrganicMatter)then  
                        call ComputeLabDOM              (index)
                        call ComputeRefDOM              
                        call ComputeLabPOM              (index)
                        call ComputeRefPOM              
                    endif
                    
                                        
                    !Oxygen
                    if (Me%Compute%Oxygen      )then
                        call ComputeDissolvedOxygen     (index)
                    end if
                    

                     !Algae
                    if (Me%Compute%Algae       )then
                        call ComputeAlgae               (index)
                    end if
                    
                    !BOD
                    if (Me%Compute%BOD         )then
                        call ComputeBOD                 (index)
                    end if

                    if (Me%Compute%Detritus    )then
                        call ComputeDetritus
                    end if
                   
                    
                    !Icarbon
                    if (Me%Compute%Icarbon     )then 
                        call ComputeICarbon             (index)
                        call ComputepH_CO2              (index)
                    endif


                    !Epiphyton
                    if (Me%Compute%Epiphyton   )then
                        call ComputeEpiphyton           (index)
                    end if
                    
                    do propi=Me%Size%PropLB, Me%Size%PropUB
                        Me%ExternalVar%Mass(propi,index) = Me%ExternalVar%Mass(propi,index) + Me%SinksSources(propi) * Me%DtDay
                    end do
                
                end if cd0
            end do do1

            nullify(Me%ExternalVar%Salinity    )
            nullify(Me%ExternalVar%Temperature )
            nullify(Me%ExternalVar%Alkalinity  )
            nullify(Me%ExternalVar%Mass        )
           
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine CEQUALW2

    
    !----------------------------------------------------------------------------
        
    subroutine TemperatureRateMultipliers(index)
    !Calculates the temperature dependence for biological and chemical rates

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index

        !Local-----------------------------------------------------------------
        real                                    :: lam1, lam2  ! Temporary auxiliar variables
        real                                    :: Temperature
!        real                                    :: NH4T1, NH4T2, NH4K1, NH4K2, NO3T1, NO3T2, NO3K1, NO3K2 
!        real                                    :: OMT1, OMT2, OMK1, OMK2
        type(T_Algae),      pointer             :: Algae
        type(T_Epiphyton),  pointer             :: Epiphyton
                                                
        ! ---------------------------------------------------------------------

        Temperature             = Me%ExternalVar%Temperature (index)

       ! NH4T1                   = Me%NH4T1 
       ! NH4T2                   = Me%NH4T2 
       ! NH4K1                   = Me%NH4K1 
       ! NH4K2                   = Me%NH4K2
                                
      !  NO3T1                   = Me%NO3T1
       ! NO3T2                   = Me%NO3T2
        !NO3K1                   = Me%NO3K1
        !NO3K2                   = Me%NO3K2
                                
       ! OMT1                    = Me%OMT1
       ! OMT2                    = Me%OMT2
       ! OMK1                    = Me%OMK1
       ! OMK2                    = Me%OMK2

        lam1                    = Rising(Temperature,Me%NH4T1,Me%NH4T2,Me%NH4K1,Me%NH4K2)
        Me%NH4TRM               = lam1/(1.0+lam1-Me%NH4K1)

        lam1                    = Rising(Temperature,Me%NO3T1,Me%NO3T2,Me%NO3K1,Me%NO3K2)
        Me%NO3TRM               = lam1/(1.0+lam1-Me%NO3K1)

        lam1                    = Rising(Temperature,Me%OMT1,Me%OMT2,Me%OMK1,Me%OMK2)
        Me%OMTRM                = lam1/(1.0+lam1-Me%OMK1)
   
        lam1          = Rising(Temperature,Me%DETT1,Me%DETT2,Me%DETK1,Me%DETK2)
        Me%DETTRM   = lam1/(1.0+lam1-Me%DETK1)

      
        Algae => Me%FirstAlgae
        !Calculation of TRM for each algae
        do while(associated(Algae))

            lam1                = Rising(Temperature, Algae%AT1, Algae%AT2, Algae%AK1, Algae%AK2)
            lam2                = Falling(Temperature, Algae%AT3, Algae%AT4, Algae%AK3, Algae%AK4)

            Algae%ATRMR         = lam1/(1.0+lam1-Algae%AK1)
            Algae%ATRMF         = lam2/(1.0+lam2-Algae%AK4)
            Algae%ATRM          = Algae%ATRMR * Algae%ATRMF

            Algae => Algae%Next
        end do



        Epiphyton => Me%FirstEpiphyton
        !Calculation of TRM for each epiphyton
        do while(associated(Epiphyton))
                                
            lam1                = Rising(Temperature, Epiphyton%ET1, Epiphyton%ET2, Epiphyton%EK1, Epiphyton%EK2)
            lam2                = Falling(Temperature, Epiphyton%ET3, Epiphyton%ET4, Epiphyton%EK3, Epiphyton%EK4)

            Epiphyton%ETRMR     = lam1/(1.0+lam1-Epiphyton%EK1)
            Epiphyton%ETRMF     = lam2/(1.0+lam2-Epiphyton%EK4)
            Epiphyton%ETRM      = Epiphyton%ETRMR * Epiphyton%ETRMF
                                      
            Epiphyton => Epiphyton%Next
        end do
      

    end subroutine TemperatureRateMultipliers

    
    !--------------------------------------------------------------------------


    subroutine ComputeKineticRates(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index

        !Local-----------------------------------------------------------------
        real                                    :: lam1, lam2  ! Temporary auxiliar variables
        real                                    :: Temperature
        integer                                 :: O, POMREF, POMLAB, DOMLAB, DOMREF, AMM
        integer                                 :: NIT, PHOSP, SIPART, SIDISS,DET
        integer                                 :: AlgaeIndex, EpiphytonIndex
        real                                    :: Shade, TopRadiation, Thickness, Gamma
        real                                    :: ASat, ESat, LTCoef 
        real                                    :: Lightlim, ATRM, ATRMR, ATRMF,AG,AP,AN,AR,AM,AE
        real                                    :: ETRM, ETRMR, ETRMF,EG,EP,EN,ER,EM,EE
        real                                    :: Nitrate, Ammonia, Algae_, Epiphyton_, Phosphorus     
        real                                    :: aux1, aux2, aux3
        real                                    :: NONZERO = 1.0E-20
        type(T_Algae),      pointer             :: Algae
        type(T_Epiphyton),  pointer             :: Epiphyton


        !-----------------------------------------------------------------------

        O           = Me%PropIndex%Oxygen
        POMREF      = Me%PropIndex%POMREF    
        POMLAB      = Me%PropIndex%POMLAB    
        DOMLAB      = Me%PropIndex%DOMLAB    
        DOMREF      = Me%PropIndex%DOMREF    
        AMM         = Me%PropIndex%Ammonia   
        NIT         = Me%PropIndex%Nitrate
        PHOSP       = Me%PropIndex%Phosphorus
        SIPART      = Me%PropIndex%SIPART
        SIDISS      = Me%PropIndex%SIDISS
        DET         = Me%PropIndex%Detritus

        Temperature = Me%ExternalVar%Temperature (index)

        !Auxiliary variables (related with oxygen limitation)for the calculation of decay rates                
        
        Me%DO1    = (1.0+SIGN(1.0, Me%ExternalVar%Mass(O,index) - Me%O2LIM))  * 0.5
        Me%DO2    = (1.0+SIGN(1.0, Me%O2LIM  - Me%ExternalVar%Mass(O,index))) * 0.5
        Me%DO3    = (1.0+SIGN(1.0, Me%ExternalVar%Mass(O,index)-1.E-10)) * 0.5

        
        
        !Decay Rates

        Me%NH4D   =  Me%NH4TRM * Me%NH4DK * &
                              Me%ExternalVar%Mass(AMM,index) * Me%DO1

        Me%NO3D   =  Me%NO3TRM * Me%NO3DK * &
                            Me%ExternalVar%Mass(NIT,index) *Me%DO2

        Me%LDOMD  =  Me%OMTRM  * Me%LDOMDK* &
                            Me%ExternalVar%Mass(DOMLAB,index)*Me%DO3

        Me%RDOMD  =  Me%OMTRM  * Me%RDOMDK* &
                            Me%ExternalVar%Mass(DOMREF,index)*Me%DO3

        Me%LPOMD  =  Me%OMTRM  * Me%LPOMDK* &
                            Me%ExternalVar%Mass(POMLAB,index)*Me%DO3

        Me%RPOMD  =  Me%OMTRM  * Me%RPOMDK* &
                            Me%ExternalVar%Mass(POMREF,index)*Me%DO3

        Me%LRDOMD =  Me%OMTRM  * Me%LRDDK * &
                            Me%ExternalVar%Mass(DOMLAB,index)*Me%DO3

        Me%LRPOMD =  Me%OMTRM  * Me%LRPDK * &
                            Me%ExternalVar%Mass(POMLAB,index)*Me%DO3

        Me%CBODD  =  Me%KBOD *Me%TBOD**(Temperature-20.0)*Me%DO3

        Me%DETD   =  Me%DETTRM * Me%SDK *Me%ExternalVar%Mass(DET,index) * Me%DO3


        !Me%CBODDK  = Me%CBODDK + Me%CBODD

        Me%Rate%Value(Me%Rate%CequalIndex%NH4D,  index) = Me%NH4D    / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%NO3D,  index) = Me%NO3D    / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%LDOMD, index) = Me%LDOMD   / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%RDOMD, index) = Me%RDOMD   / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%LPOMD, index) = Me%LPOMD   / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%RPOMD, index) = Me%RPOMD   / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%LRDOMD,index) = Me%LRDOMD  / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%LRPOMD,index) = Me%LRPOMD  / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%CBODD, index) = Me%CBODD   / 24.0 / 60.0 / 60.0

       

! Algal rates

        Shade           = 1. !Shade effect can be computed in CEQUALW2 we have chosen to consider no shade(Shade=1)
        TopRadiation    = Me%ExternalVar%ShortWaveRadiation(index)
        Thickness       = Me%ExternalVar%Thickness(index)
        Gamma           = Me%ExternalVar%LightExtCoefField(index)

  
        if (Me%Compute%Algae) then

            Algae => Me%FirstAlgae

            do while(associated(Algae))

                ASat   = Algae%Asat
                
                LTCoef = TopRadiation*Shade/ASat
                                                          

                !**** Growth Limitations


                Lam1           = LTCoef
                Lam2           = LTCoef*exp(-Gamma*Thickness)
                

                Algae%Lightlim(index) = 2.718282*(exp(-Lam2)-exp(-Lam1))/(Gamma*Thickness)
                Algae%Plim(index)     = 1.0                                               
                Algae%Nlim(index)     = 1.0                                               
                Algae%Slim(index)     = 1.0                                               

                !Phosphorus Limitation
                if (Algae%AHSP  /= 0.0) &
                Algae%Plim(index) = Me%ExternalVar%Mass(PHOSP,index)/      &
                (Me%ExternalVar%Mass(PHOSP,index)+Algae%AHSP)                                

                ! Nitrogen Limitation
                if (Algae%AHSN  /= 0.0) &
                Algae%Nlim(index) = (Me%ExternalVar%Mass(AMM,index)+Me%ExternalVar%Mass(NIT,index))/ &
                (Me%ExternalVar%Mass(AMM,index)+Me%ExternalVar%Mass(NIT,index)+Algae%AHSN)                      

                !Silica Limitattion
                if (Algae%AHSSI /= 0.0) &
                Algae%Slim(index) = Me%ExternalVar%Mass(SIDISS,index)/     &
                (Me%ExternalVar%Mass(SIDISS,index)+Algae%AHSSI)                                          


                Algae%OverallLim (index)    = min(Algae%Plim(index), Algae%Nlim(index),Algae%Slim(index),Algae%Lightlim(index))

                AlgaeIndex = Algae%PropIndex   

                ATRM    =  Algae%ATRM
                ATRMR   =  Algae%ATRMR
                ATRMF   =  Algae%ATRMF
                AG      =  Algae%AG
                AP      =  Algae%AP
                AN      =  Algae%AN
                AR      =  Algae%AR
                AM      =  Algae%AM
                AE      =  Algae%AE

                Nitrate =  Me%ExternalVar%Mass(NIT,index)
                Ammonia =  Me%ExternalVar%Mass(AMM,index)
                Algae_ =  Me%ExternalVar%Mass(AlgaeIndex,index)
                Phosphorus =  Me%ExternalVar%Mass(PHOSP,index)

                Lightlim =  Algae%Lightlim(index)


                Aux1 = ATRM * AG * Algae%OverallLim (index)
                Aux2 = Phosphorus/(AP * Me%DTDay * Algae_ + NONZERO) 
                Aux3 = (Ammonia + Nitrate)/(AN * Me%DTDay * Algae_ + NONZERO)  

                Algae%AGR =  min (Aux1, Aux2, Aux3) 

                Algae%ARR =  ATRM * AR * Me%DO3

                Algae%AMR = (ATRMR + 1.0-ATRMF) * AM
                
                Algae%AER =  MIN((1.0 - Lightlim)*AE * ATRM , Algae%AGR)                                       

                !corrigir mais tarde (so da rates da ultima alga)!!!pina
                Me%Rate%Value(Me%Rate%CequalIndex%ANLim,        index)  = Algae%NLim(index)       
                Me%Rate%Value(Me%Rate%CequalIndex%APLim,        index)  = Algae%PLim(index)       
                Me%Rate%Value(Me%Rate%CequalIndex%ASLim,        index)  = Algae%SLim(index)       
                Me%Rate%Value(Me%Rate%CequalIndex%ALightLim,    index)  = Algae%LightLim(index)   
                Me%Rate%Value(Me%Rate%CequalIndex%AOverallLim,  index)  = Algae%OverallLim(index) 

                
                Algae => Algae%Next
            end do
        endif



!Epiphyton rates
    if (Me%Compute%Epiphyton) then

            
        Epiphyton => Me%FirstEpiphyton
        
        do while(associated(Epiphyton))

            EpiphytonIndex = Epiphyton%PropIndex 
            Epiphyton_     = Me%ExternalVar%Mass(EpiphytonIndex,index)

            ESat   = Epiphyton%Esat
            LTCoef = TopRadiation*Shade/ESat

            Lam1           = LTCoef
            Lam2           = LTCoef*exp(-Gamma*Thickness)
            

            Epiphyton%lightlim(index) = 2.718282*(exp(-Lam2)-exp(-Lam1))/(Gamma*Thickness)
            Epiphyton%Plim(index)     = 1.0
            Epiphyton%Nlim(index)     = 1.0
            Epiphyton%Slim(index)     = 1.0

            !Phosphorus Limitation
            if (Epiphyton%EHSP  /= 0.0) &
            Epiphyton%Plim(index) = Me%ExternalVar%Mass(PHOSP,index)/      &
            (Me%ExternalVar%Mass(PHOSP,index)+Epiphyton%EHSP)                                
      
            ! Nitrogen Limitation
            if (Epiphyton%EHSN  /= 0.0) &
            Epiphyton%Nlim(index) = (Me%ExternalVar%Mass(AMM,index)+Me%ExternalVar%Mass(NIT,index))/ &
            (Me%ExternalVar%Mass(AMM,index)+Me%ExternalVar%Mass(NIT,index)+Epiphyton%EHSN)
      
            !Silica Limitattion
            if (Epiphyton%EHSSI /= 0.0) &
            Epiphyton%Slim(index) = Me%ExternalVar%Mass(SIDISS,index)/     &
            (Me%ExternalVar%Mass(SIDISS,index)+Epiphyton%EHSSI)                                          
         
      
            Epiphyton%OverallLim (index)   = min(Epiphyton%Plim(index),      &
                                                 Epiphyton%Nlim(index),      &
                                                 Epiphyton%Slim(index),      &
                                                 Epiphyton%Lightlim(index))
            
            ETRM =  Epiphyton%ETRM
            ETRMR =  Epiphyton%ETRMR
            ETRMF =  Epiphyton%ETRMF
            EG =  Epiphyton%EG
            EP =  Epiphyton%EP
            EN =  Epiphyton%EN
            ER =  Epiphyton%ER
            EM =  Epiphyton%EM
            EE =  Epiphyton%EE

            Phosphorus =  Me%ExternalVar%Mass(PHOSP,index)     
            Nitrate =  Me%ExternalVar%Mass(NIT,index)
            Ammonia =  Me%ExternalVar%Mass(AMM,index)
            Lightlim =  Epiphyton%Lightlim(index)
      
            Epiphyton%EGR = min( ETRM * EG * Epiphyton%OverallLim(index),   &
                                Phosphorus/(EP * Me%DTDay * Epiphyton_ + NONZERO),  &
                                (Ammonia + Nitrate)/(EN * Me%DTDay * Epiphyton_ + NONZERO))



            Epiphyton%ERR =  ETRM * ER * Me%DO3

            Epiphyton%EMR = (ETRMR + 1.0-ETRMF) * EM

            Epiphyton%EER =  MIN((1.0 - lightlim)*EE * ETRM , Epiphyton%EGR)

            Me%Rate%Value(Me%Rate%CequalIndex%ENLim,        index)  = Epiphyton%NLim(index) 
            Me%Rate%Value(Me%Rate%CequalIndex%EPLim,        index)  = Epiphyton%PLim(index)  
            Me%Rate%Value(Me%Rate%CequalIndex%ESLim,        index)  = Epiphyton%SLim(index)  
            Me%Rate%Value(Me%Rate%CequalIndex%ELightLim,    index)  = Epiphyton%LightLim(index) 
            Me%Rate%Value(Me%Rate%CequalIndex%EOverallLim,  index)  = Epiphyton%OverallLiM(index) 

            Epiphyton => Epiphyton%Next
        end do

         

    
    end if


        !----------------------------------------------------------------------

    end subroutine ComputeKineticRates

    
    !--------------------------------------------------------------------------

    
    subroutine ComputePhosphorus(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index

        !Local-----------------------------------------------------------------
        integer                                 :: O, POMREF, POMLAB, DOMLAB, DOMREF, AM, NIT
        integer                                 :: PHOSP, SIPART, SIDISS, BOD
        integer                                 :: AlgaeIndex, EpiphytonIndex
        real                                    :: PO4BOD, PO4AR, PO4AG, PO4ER, PO4EG 
        real                                    :: PO4EP,PO4AP,PO4POM,PO4DOM,PO4OM, PO4SS
        type(T_Algae),      pointer             :: Algae
        type(T_Epiphyton),  pointer             :: Epiphyton

        !-----------------------------------------------------------------------

        O       = Me%PropIndex%Oxygen
        POMREF  = Me%PropIndex%pomref    
        POMLAB  = Me%PropIndex%pomlab    
        DOMLAB  = Me%PropIndex%domlab    
        DOMREF  = Me%PropIndex%domref    
        AM      = Me%PropIndex%Ammonia   
        NIT     = Me%PropIndex%Nitrate
        PHOSP   = Me%PropIndex%Phosphorus
        SIPART  = Me%PropIndex%sipart    
        SIDISS  = Me%PropIndex%sidiss    
        BOD     = Me%PropIndex%BOD

        if (Me%Compute%BOD) then
           PO4BOD = Me%CBODD  * Me%ExternalVar%Mass(BOD,index) * Me%BODP
        else
           PO4BOD = 0.0
        endif

        PO4AR = 0.0
        PO4AG = 0.0
        PO4ER = 0.0
        PO4EG = 0.0
        
        Algae => Me%FirstAlgae

        do while(associated(Algae))

            AlgaeIndex = Algae%PropIndex   

            PO4AG = PO4AG + Algae%AGR * &
                             Me%ExternalVar%Mass(AlgaeIndex,index) * Algae%AP
                 
            PO4AR = PO4AR + Algae%ARR * &
                             Me%ExternalVar%Mass(AlgaeIndex,index) * Algae%AP

            Algae => Algae%Next
        enddo


        Epiphyton => Me%FirstEpiphyton

        do while(associated(Epiphyton))

            EpiphytonIndex = Epiphyton%PropIndex  
            
            PO4EG = PO4EG + Epiphyton%EGR * &
                                Me%ExternalVar%Mass(EpiphytonIndex,index) * Epiphyton%EP

            PO4ER = PO4ER +Epiphyton%ERR * &
                                Me%ExternalVar%Mass(EpiphytonIndex,index) * Epiphyton%EP

            Epiphyton => Epiphyton%Next
        enddo
      
        PO4EP  = PO4ER  - PO4EG
        PO4AP  = PO4AR  - PO4AG
        PO4POM = Me%ORGP * (Me%LPOMD + Me%RPOMD)
        PO4DOM = Me%ORGP * (Me%LDOMD + Me%RDOMD)

        PO4OM  = PO4POM + PO4DOM

        
        PO4SS  = PO4AP  + PO4EP + PO4OM + PO4BOD
        Me%SinksSources(PHOSP) = PO4SS
                           
       
        Me%Rate%Value(Me%Rate%CequalIndex%PO4ER, index)         = PO4ER / 24.0 / 60.0 / 60.0 
        Me%Rate%Value(Me%Rate%CequalIndex%PO4EG, index)         = PO4EG / 24.0 / 60.0 / 60.0 
        Me%Rate%Value(Me%Rate%CequalIndex%PO4AR, index)         = PO4AR / 24.0 / 60.0 / 60.0 
        Me%Rate%Value(Me%Rate%CequalIndex%PO4AG, index)         = PO4AG / 24.0 / 60.0 / 60.0 
        Me%Rate%Value(Me%Rate%CequalIndex%PO4OM, index)         = PO4OM / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%PO4BOD,index)         = PO4BOD / 24.0 / 60.0 / 60.0
    
    
    end subroutine ComputePhosphorus

    
    !--------------------------------------------------------------------------
    
    
    subroutine ComputeAmmonia(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index

        !Local-----------------------------------------------------------------
        integer                                 :: BOD, AM, NIT, POMREF, POMLAB, DOMLAB, DOMREF
        integer                                 :: AlgaeIndex, EpiphytonIndex
        real                                    :: Ammonia, Nitrate, Algae_, Epiphyton_
        real                                    :: NONZERO = 1.0E-20
        real                                    :: NH4PR,   NH4AG, NH4AR, NH4ER, NH4EG, NH4BOD
        real                                    :: ANPR, AHSN 
        real                                    :: ENPR, EHSN 
        integer                                 :: ANEQN
        integer                                 :: ENEQN
        real                                    :: NH4EP, NH4AP, NH4DOM, NH4POM, NH4OM, NH4SS
        type(T_Algae),      pointer             :: Algae
        type(T_Epiphyton),  pointer             :: Epiphyton

        !-----------------------------------------------------------------------

        BOD     = Me%PropIndex%BOD
        AM      = Me%PropIndex%Ammonia
        NIT     = Me%PropIndex%Nitrate
        POMREF  = Me%PropIndex%pomref    
        POMLAB  = Me%PropIndex%pomlab    
        DOMLAB  = Me%PropIndex%domlab    
        DOMREF  = Me%PropIndex%domref    


        if (Me%Compute%BOD) then 
          NH4BOD  = Me%CBODD * &
                  Me%ExternalVar%Mass(BOD,index) * Me%BODN 
        else
          NH4BOD  = 0.0
        
        endif                                                        
    
        Ammonia = Me%ExternalVar%Mass(AM,index)
        Nitrate = Me%ExternalVar%Mass(NIT,index)

        NH4AR = 0.0
        NH4AG = 0.0
        NH4ER = 0.0
        NH4EG = 0.0

        Algae   => Me%FirstAlgae

        do while(associated(Algae))

            ANPR       = Algae%ANPR   !Algal ammonium preference factor
            ANEQN      = Algae%ANEQN  !Equation number for algal ammonium preference (either 1 or 2)
            AHSN       = Algae%AHSN   !Algal half-saturation for nitrogen limited growth

            AlgaeIndex = Algae%PropIndex   
            Algae_     = Me%ExternalVar%Mass(AlgaeIndex,index)

            if (ANEQN.EQ.1) then

                NH4PR = Ammonia / ( Ammonia + Nitrate + NONZERO)

            else if (ANEQN.EQ.2) then
                NH4PR = Ammonia * Nitrate /((ANPR + Ammonia)*(ANPR + Nitrate)) + &
                        Ammonia * ANPR /((Nitrate + Ammonia + NONZERO)*(ANPR+Nitrate))
            end if

            if (AHSN > 0.0) then
                NH4AG = NH4AG + Algae%AGR * Algae_ * Algae%AN * NH4PR
                NH4AR = NH4AR + Algae%ARR * Algae_ * Algae%AN
            endif

            Algae => Algae%Next
        enddo



        Epiphyton => Me%FirstEpiphyton

        do while(associated(Epiphyton))

            ENPR       = Epiphyton%ENPR
            ENEQN      = Epiphyton%ENEQN  !Equation number for algal ammonium preference (either 1 or 2)
            EHSN       = Epiphyton%EHSN   !Epiphyton half-saturation for nitrogen limited growth

            EpiphytonIndex = Epiphyton%PropIndex    
            Epiphyton_     = Me%ExternalVar%Mass(EpiphytonIndex,index)

            if (ENEQN == 1) then
                NH4PR = Ammonia/(Ammonia + Nitrate + NONZERO)
            else if (ENEQN == 2) then
                NH4PR = Ammonia * Nitrate/((ENPR + Ammonia)*(ENPR + Nitrate))+Ammonia * ENPR /(( Nitrate +Ammonia + NONZERO)*(ENPR &
                + Nitrate))
            end if                                                                                                         

            NH4EG = NH4EG + Epiphyton%EGR * Epiphyton_ * Epiphyton%EN * NH4PR                                                      
            NH4ER = NH4ER + Epiphyton%ERR * Epiphyton_ * Epiphyton%EN

            Epiphyton => Epiphyton%Next

        enddo

      
        NH4EP  =  NH4ER - NH4EG
        NH4AP  =  NH4AR - NH4AG
        NH4POM =  Me%ORGN * (Me%LPOMD + Me%RPOMD)
        NH4DOM =  Me%ORGN * (Me%LDOMD + Me%RDOMD)
        NH4OM  =  NH4DOM + NH4POM

        NH4SS   =  NH4AP + NH4EP + NH4OM  + NH4BOD - Me%NH4D   
        
        Me%SinksSources(AM) = NH4SS
        
        
        Me%Rate%Value(Me%Rate%CequalIndex%NH4ER, index)         = NH4ER  / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%NH4EG, index)         = NH4EG  / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%NH4AR, index)         = NH4AR  / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%NH4AG, index)         = NH4AG  / 24.0 / 60.0 / 60.0 
        Me%Rate%Value(Me%Rate%CequalIndex%NH4OM, index)         = NH4OM  / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%NH4BOD,index)         = NH4BOD / 24.0 / 60.0 / 60.0

    end subroutine ComputeAmmonia

    
    !---------------------------------------------------------------------------
    
    
    subroutine ComputeNitrate(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: index

        !Local-----------------------------------------------------------------
        integer                                     :: BOD, AM, NIT, POMREF, POMLAB, DOMLAB, DOMREF
        integer                                     :: AlgaeIndex, EpiphytonIndex
        real                                        :: Ammonia, Algae_, Nitrate, Epiphyton_
        real                                        :: NONZERO = 1.0E-20
        real                                        :: NO3PR, NO3AG, NO3EG, NO3SS
        real                                        :: ANPR, AHSN  
        real                                        :: ENPR, EHSN 
        integer                                     :: ANEQN
        integer                                     :: ENEQN
        type(T_Algae),      pointer                 :: Algae
        type(T_Epiphyton),  pointer                 :: Epiphyton
        !-----------------------------------------------------------------------

        BOD     = Me%PropIndex%BOD
        AM      = Me%PropIndex%Ammonia
        NIT     = Me%PropIndex%Nitrate
        POMREF  = Me%PropIndex%pomref    
        POMLAB  = Me%PropIndex%pomlab    
        DOMLAB  = Me%PropIndex%domlab    
        DOMREF  = Me%PropIndex%domref    
  
        
   
        Ammonia = Me%ExternalVar%Mass(AM,index)
        Nitrate = Me%ExternalVar%Mass(NIT,index)


        NO3AG = 0.0
        NO3EG = 0.0
        Algae => Me%FirstAlgae

        do while(associated(Algae))

            ANPR       = Algae%ANPR
            ANEQN      = Algae%ANEQN
            AHSN       = Algae%AHSN

            AlgaeIndex = Algae%PropIndex   
            Algae_     = Me%ExternalVar%Mass(AlgaeIndex,index)


            if (ANEQN.EQ.1) then                                                                                      

                NO3PR = 1.0 - Ammonia/(Ammonia + Nitrate + NONZERO)                                                     

            else if (ANEQN.EQ.2) then                                                                             

                NO3PR = 1.0-(Ammonia*Nitrate/((ANPR+Ammonia)*(ANPR+Nitrate))+&
                            Ammonia*ANPR/((Nitrate + Ammonia + NONZERO)  &
                            *(ANPR+Nitrate)))                                                                               
            end if                                                                                                              

            if(AHSN.GT.0.0) then
                NO3AG = NO3AG + Algae%AGR * Algae_ * NO3PR * Algae%AN                              
            endif

            Algae => Algae%Next
        enddo


        Epiphyton => Me%FirstEpiphyton

        do while(associated(Epiphyton))

            ENPR       = Epiphyton%ENPR
            ENEQN      = Epiphyton%ENEQN  !Equation number for EPIPHYTON ammonium preference (either 1 or 2)
            EHSN       = Epiphyton%EHSN   !Epiphyton half-saturation for nitrogen limited growth

            EpiphytonIndex = Epiphyton%PropIndex 
            Epiphyton_     = Me%ExternalVar%Mass(EpiphytonIndex,index)


            if (ENEQN == 1) then 

                NO3PR = 1.0-Ammonia/(Ammonia+Nitrate+NONZERO)  

            else if (ENEQN.EQ.2) then                                                                                   

                NO3PR = 1.0-(Ammonia*Nitrate/((ENPR+Ammonia)*(ENPR+Nitrate)) &
                            +Ammonia*ENPR/((Nitrate+Ammonia+NONZERO)     &
                        *(ENPR+Nitrate)))                                                                                
            end if

            NO3EG = NO3EG + Epiphyton%EGR * Epiphyton_ * NO3PR * Epiphyton%EN

            Epiphyton => Epiphyton%Next

        enddo
   
        NO3SS = Me%NH4D  - Me%NO3D  - NO3AG - NO3EG 

        Me%SinksSources(NIT) = NO3SS

        
        Me%Rate%Value(Me%Rate%CequalIndex%NO3AG, index)         = NO3AG / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%NO3EG, index)         = NO3EG / 24.0 / 60.0 / 60.0

    end subroutine ComputeNitrate

    
    !--------------------------------------------------------------------------
    
    
    subroutine ComputeDissolvedSilica(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: index

        !Local-----------------------------------------------------------------
        integer                                     :: PSI, DSI
        integer                                     :: AlgaeIndex, EpiphytonIndex
        real                                        :: Algae_, Epiphyton_, SilicaPart
        real                                        :: DSIAG, DSIEG, DSID, DSISS
        type(T_Algae),      pointer                 :: Algae
        type(T_Epiphyton),  pointer                 :: Epiphyton
        !-----------------------------------------------------------------------
        
        PSI        = Me%PropIndex%sipart
        DSI        = Me%PropIndex%sidiss 
        SilicaPart = Me%ExternalVar%Mass(PSI,index)
        
        DSIAG = 0.0
        DSIEG = 0.0

        Algae => Me%FirstAlgae

        do while(associated(Algae))
             
            AlgaeIndex = Algae%PropIndex   
            Algae_     = Me%ExternalVar%Mass(AlgaeIndex,index)

            DSIAG = DSIAG + Algae%AGR * Algae_ * Algae%ASI

            Algae => Algae%Next
        enddo


        Epiphyton => Me%FirstEpiphyton

        do while(associated(Epiphyton))

            EpiphytonIndex = Epiphyton%PropIndex 
            Epiphyton_     = Me%ExternalVar%Mass(EpiphytonIndex,index)

            DSIEG = DSIEG + Epiphyton%EGR * Epiphyton_ * Epiphyton%ESI

            Epiphyton => Epiphyton%Next

        enddo

        DSID  =  Me%PSIDK * SilicaPart
        
        DSISS  =  DSID  - DSIAG - DSIEG 

        Me%SinksSources(DSI) = DSISS

       
        Me%Rate%Value(Me%Rate%CequalIndex%DSIAG, index)         = DSIAG  / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%DSIEG, index)         = DSIEG  / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%DSID,  index)         = DSID   / 24.0 / 60.0 / 60.0

    
    end subroutine ComputeDissolvedSilica 

    
    !--------------------------------------------------------------------------
    

    subroutine ComputeParticulateSilica(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index

        !Local-----------------------------------------------------------------
        integer                                 :: PSI
        real                                    :: SilicaPart
        real                                    :: PSIAM, PSID, PSISS
        type(T_Algae),      pointer             :: Algae
        
        !-----------------------------------------------------------------------

        PSI        = Me%PropIndex%SIPART
        SilicaPart = Me%ExternalVar%Mass(PSI,index)

        PSIAM = 0.0

        Algae => Me%FirstAlgae

        do while(associated(Algae))
             
            PSIAM = PSIAM + Algae%AMR * SilicaPart * Algae%ASI 

            Algae => Algae%Next
        enddo

        PSID   =  Me%PSIDK * SilicaPart
        
        PSISS =  PSIAM - PSID 

        Me%SinksSources(PSI) = PSISS
        
        Me%Rate%Value(Me%Rate%CequalIndex%PSIAM, index)         = PSIAM / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%PSID,  index)         = PSID  / 24.0 / 60.0 / 60.0
       
    
    end subroutine ComputeParticulateSilica
    !--------------------------------------------------------------------------
    
    !IRON - Only vertical physsical procesess
    
    !--------------------------------------------------------------------------

    subroutine ComputeLabDOM(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: index

        !Local-----------------------------------------------------------------
        integer                                     :: EpiphytonIndex, DOMLAB 
        real                                        :: Algae_, Epiphyton_
        real                                        :: LDOMAP, LDOMEP, LDOMSS
        type(T_Algae),      pointer                 :: Algae
        type(T_Epiphyton),  pointer                 :: Epiphyton

        !-----------------------------------------------------------------------

        DOMLAB  = Me%PropIndex%DomLab

        LDOMAP = 0.0
        LDOMEP = 0.0

        Algae => Me%FirstAlgae

        do while(associated(Algae))
                
            Algae_     = Me%ExternalVar%Mass(Algae%PropIndex, index)

            LDOMAP = LDOMAP +(Algae%AER +(1.0 - Algae%APOM) * Algae%AMR) * Algae_

            Algae => Algae%Next
        enddo


        Epiphyton => Me%FirstEpiphyton

        do while(associated(Epiphyton))

            EpiphytonIndex = Epiphyton%PropIndex 
            Epiphyton_     = Me%ExternalVar%Mass(EpiphytonIndex,index)

            LDOMEP = LDOMEP + (Epiphyton%EER+(1.0- Epiphyton%EPOM) * Epiphyton%EMR)*Epiphyton_
            
            Epiphyton => Epiphyton%Next

        enddo

        LDOMSS = LDOMAP + LDOMEP - Me%LDOMD - Me%LRDOMD

        Me%SinksSources(DOMLAB) = LDOMSS        
        
        Me%Rate%Value(Me%Rate%CequalIndex%LDOMAP, index)         = LDOMAP / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%LDOMEP, index)         = LDOMEP / 24.0 / 60.0 / 60.0
   
    end subroutine ComputeLabDOM

    !--------------------------------------------------------------------------

    subroutine ComputeRefDOM

        !Local-----------------------------------------------------------------
        integer                                 :: DOMREF
        real                                    :: RDOMSS
        !----------------------------------------------------------------------
 
        DOMREF  = Me%PropIndex%DomRef

        RDOMSS = Me%LRDOMD - Me%RDOMD

        Me%SinksSources(DOMREF) = RDOMSS
       
    end subroutine ComputeRefDOM

    
    !--------------------------------------------------------------------------
      
    
    subroutine ComputeLabPOM(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index

        !Local-----------------------------------------------------------------
        integer                                 :: AlgaeIndex, POMLAB
        real                                    :: Algae_
        real                                    :: LPOMAP, LPOMSS
        type(T_Algae),      pointer             :: Algae

        !-----------------------------------------------------------------------
        
        POMLAB  = Me%PropIndex%pomlab


        LPOMAP = 0.0
        Algae => Me%FirstAlgae

        do while(associated(Algae))
            AlgaeIndex = Algae%PropIndex   
            Algae_     = Me%ExternalVar%Mass(AlgaeIndex,index)

            LPOMAP = LPOMAP +  Algae%APOM * Algae%AMR * Algae_
            Algae => Algae%Next
        enddo

        
        LPOMSS  =  LPOMAP - Me%LPOMD 

       
        Me%SinksSources(POMLAB) = LPOMSS
        
        
        Me%Rate%Value(Me%Rate%CequalIndex%LPOMAP, index)         = LPOMAP / 24.0 / 60.0 / 60.0

    end subroutine ComputeLabPOM
      
    !---------------------------------------------------------------------------
      
    subroutine ComputeRefPOM

                
        !Local-----------------------------------------------------------------
        integer                                 :: POMREF
        real                                    :: RPOMSS

        !-----------------------------------------------------------------------

        POMREF = Me%PropIndex%pomref

       
        RPOMSS =  Me%LRPOMD - Me%RPOMD 


        Me%SinksSources(POMREF) = RPOMSS
       

    end subroutine ComputeRefPOM
    
    !--------------------------------------------------------------------------
      
    subroutine ComputeAlgae(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index

        !Local-----------------------------------------------------------------
        integer                                 :: AlgaeIndex
        real                                    :: Algae_, ASS
        type(T_Algae),      pointer             :: Algae
        !-----------------------------------------------------------------------
        
        Algae => Me%FirstAlgae

        do while(associated(Algae))
            AlgaeIndex = Algae%PropIndex   
            Algae_     = Me%ExternalVar%Mass(AlgaeIndex,index)

            ASS = (Algae%AGR - Algae%AER - Algae%AMR - Algae%ARR) * Algae_
            
            
            Me%SinksSources(AlgaeIndex) = ASS                             
            
            Algae => Algae%Next
        enddo

    end subroutine ComputeAlgae

    !--------------------------------------------------------------------------

    subroutine ComputeBOD(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index

        !Local-----------------------------------------------------------------
        real                                    :: CBODSS
        integer                                 :: BOD 
        !-----------------------------------------------------------------------

        BOD    = Me%PropIndex%BOD
        CBODSS = - Me%CBODD * Me%ExternalVar%Mass(BOD,index)


        Me%SinksSources(BOD) = CBODSS
        

    end subroutine ComputeBOD
    
    
    !--------------------------------------------------------------------------


    subroutine ComputeDissolvedOxygen(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: index

        !Local-----------------------------------------------------------------
        integer                                     :: BOD 
        integer                                     :: AlgaeIndex, EpiphytonIndex, O
        real                                        :: DOAP, DOAR, DOEP, DOER, DOBOD
        real                                        :: Algae_, Epiphyton_, Temperature, Salinity 
        real                                        :: DOPOM, DODOM, DOOM, DONIT, DOSS !,DOSAT
        logical                                     :: SALTWATER
        type(T_Algae),      pointer                 :: Algae
        type(T_Epiphyton),  pointer                 :: Epiphyton


        !-----------------------------------------------------------------------

        BOD                = Me%PropIndex%BOD
        O                  = Me%PropIndex%Oxygen
        
        if (Me%Compute%BOD) then
            DOBOD =  Me%RBOD * Me%CBODD * Me%ExternalVar%Mass(BOD,index)
        else
            DOBOD = 0.0
        endif


        DOAP = 0.0
        DOAR = 0.0
        DOER = 0.0
        DOEP = 0.0  

        Algae => Me%FirstAlgae

        do while(associated(Algae))
            AlgaeIndex = Algae%PropIndex   
            Algae_     = Me%ExternalVar%Mass(AlgaeIndex,index)

            DOAP = DOAP + Algae%AGR * Algae_ * Algae%O2AG
            DOAR = DOAR + Algae%ARR * Algae_ * Algae%O2AR
            
            Algae => Algae%Next
        enddo

        Epiphyton => Me%FirstEpiphyton

        do while(associated(Epiphyton))

            EpiphytonIndex = Epiphyton%PropIndex 
            Epiphyton_     = Me%ExternalVar%Mass(EpiphytonIndex,index)

            DOEP  = DOEP + Epiphyton%EGR * Epiphyton_ * Epiphyton%O2EG
            DOER  = DOER + Epiphyton%ERR * Epiphyton_ * Epiphyton%O2ER
            
            Epiphyton => Epiphyton%Next

        enddo

        DOPOM      = (Me%LPOMD + Me%RPOMD) * Me%O2OM
        DODOM      = (Me%LDOMD + Me%RDOMD) * Me%O2OM
        DOOM       =  DOPOM + DODOM + DOBOD
        DONIT      =  Me%NH4D * Me%O2NH4
        
        
        DOSS  =  DOAP + DOEP -DOAR -DOER -DOOM -DONIT 
        
        Me%SinksSources(O) = DOSS
        
        Temperature = Me%ExternalVar%Temperature(index)
        Salinity    = Me%ExternalVar%Salinity(index)
        if (Salinity .le. 0.5) SALTWATER =.false.

                                                                      
        !DOSAT = OxygenSaturationCEQUALW2(Temperature, Salinity, Me%PALT, SALTWATER)

       
        Me%Rate%Value(Me%Rate%CequalIndex%DOAP, index)         = DOAP / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%DOEP, index)         = DOEP / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%DOAR, index)         = DOAR / 24.0 / 60.0 / 60.0 
        Me%Rate%Value(Me%Rate%CequalIndex%DOER, index)         = DOER / 24.0 / 60.0 / 60.0 
        Me%Rate%Value(Me%Rate%CequalIndex%DOOM, index)         = DOOM / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%DONIT,index)         = DONIT / 24.0 / 60.0 / 60.0
 
    end subroutine ComputeDissolvedOxygen   

    !--------------------------------------------------------------------------

    subroutine ComputeICarbon(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: index

        !Local-----------------------------------------------------------------
        integer                                     :: BOD 
        integer                                     :: AlgaeIndex, EpiphytonIndex, ICarbonIndex
        real                                        :: Epiphyton_, Algae_
        real                                        :: ICarbonAP, ICarbonEP, ICarbonBOD, ICarbonSS
        type(T_Algae),      pointer                 :: Algae
        type(T_Epiphyton),  pointer                 :: Epiphyton

        !-----------------------------------------------------------------------

        BOD          = Me%PropIndex%BOD
        ICarbonIndex = Me%PropIndex%ICarbon

        if (Me%Compute%BOD) then
            ICarbonBOD   = Me%CBODD * Me%ExternalVar%Mass(BOD,index) * Me%BODC
        else 
            ICarbonBOD =0.0
        endif

        ICarbonAP = 0.0
        ICarbonEP = 0.0

        Algae => Me%FirstAlgae

        do while(associated(Algae))
            AlgaeIndex = Algae%PropIndex   
            Algae_     = Me%ExternalVar%Mass(AlgaeIndex,index)
            
            ICarbonAP = ICarbonAP + Algae%AC * (Algae%ARR - Algae%AGR) * Algae_
           
            Algae => Algae%Next

        enddo



        Epiphyton => Me%FirstEpiphyton

        do while(associated(Epiphyton))

            EpiphytonIndex = Epiphyton%PropIndex 
            Epiphyton_     = Me%ExternalVar%Mass(EpiphytonIndex,index)

            ICarbonEP  = ICarbonEP + Epiphyton%EC * (Epiphyton%ERR - Epiphyton%EGR) * Epiphyton_

            Epiphyton => Epiphyton%Next

        enddo


        ICarbonSS = ICarbonAP + ICarbonEP + ICarbonBOD + &
                       Me%ORGC * (Me%LPOMD + Me%RPOMD&
                       + Me%LDOMD + Me%RDOMD )                    
                                                                             

        Me%SinksSources(ICarbonindex) = ICarbonSS
        
        Me%Rate%Value(Me%Rate%CequalIndex%ICarbonAP, index)         = ICarbonAP / 24.0 / 60.0 / 60.0
        Me%Rate%Value(Me%Rate%CequalIndex%ICarbonEP, index)         = ICarbonEP / 24.0 / 60.0 / 60.0 
        Me%Rate%Value(Me%Rate%CequalIndex%ICarbonBOD,index)         = ICarbonBOD / 24.0 / 60.0 / 60.0

    end subroutine ComputeICarbon
   
    !--------------------------------------------------------------------------

    subroutine ComputeDetritus

        !Local-----------------------------------------------------------------

        integer                                 :: DetIndex

        !-----------------------------------------------------------------------

        DetIndex = Me%PropIndex%Detritus


        Me%SinksSources(DetIndex) = -Me%DETD
        

    end subroutine ComputeDetritus
   
    !--------------------------------------------------------------------------

    
    subroutine ComputeEpiphyton(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: index

        !Local-----------------------------------------------------------------
        integer                                     :: EpiphytonIndex
        real                                        :: Epiphyton_ 
        real                                        :: ESS 
        type(T_Epiphyton),  pointer                 :: Epiphyton

        !-----------------------------------------------------------------------
       
        

        Epiphyton => Me%FirstEpiphyton
        
        do while(associated(Epiphyton))

            EpiphytonIndex = Epiphyton%PropIndex 
            Epiphyton_     = Me%ExternalVar%Mass(EpiphytonIndex,index)


            
            ESS =   (Epiphyton%EGR - Epiphyton%ERR - Epiphyton%EMR - Epiphyton%EER) * Epiphyton_


            Me%SinksSources(EpiphytonIndex) = ESS
                                                       
            Epiphyton => Epiphyton%Next

        enddo


   end subroutine ComputeEpiphyton
  
   !--------------------------------------------------------------------------

    subroutine ComputepH_CO2(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index

        !Local-----------------------------------------------------------------
        integer                                 :: pHindex, ICarbonindex, CO2index, CO3index, HCO3index
        integer                                 :: iter, n
        real                                    :: ICarbon, Alkalinity, pH, Temperature, Salinity, CO2, HCO3, CO3
        real                                    :: S2, CART, ALKT, T1K, SQRS2 
        real                                    :: DH1, DH2, H2CO3T, CO3T, OH, HCO3T
        real                                    :: KW, K1, K2, PHT, INCR, HION, BICART, F
        !-----------------------------------------------------------------------
        pHindex        = Me%PropIndex%pH
        ICarbonindex   = Me%PropIndex%ICarbon
        CO2index       = Me%PropIndex%CO2
        CO3index       = Me%PropIndex%CO3
        HCO3index      = Me%PropIndex%HCO3

        pH             = Me%ExternalVar%Mass(pHindex,index)
        ICarbon        = Me%ExternalVar%Mass(ICarbonindex,index)

        Alkalinity     = Me%ExternalVar%Alkalinity(index)
        Temperature    = Me%ExternalVar%Temperature(index)
        Salinity       = Me%ExternalVar%Salinity(index)


        !pH and carbonate species
        CART = ICarbon/12000.0
        ALKT = Alkalinity/5.0E+04
        T1K  = Temperature + 273.15

        !**** Ionic strength
        if (Salinity < 0.5) then
            S2 = 2.5E-05*Salinity
        else 
            S2 = 1.47E-3+1.9885E-2 * Salinity+3.8E-4 * Salinity**2
        endif

        !**** Debye-Huckel terms and activity coefficients

        SQRS2  =  sqrt(S2)
        DH1    = -0.5085*SQRS2/(1.0+1.3124*SQRS2)+4.745694E-03+4.160762E-02*S2-9.284843E-03*S2*S2
        DH2    = -2.0340*SQRS2/(1.0+1.4765*SQRS2)+1.205665E-02+9.715745E-02*S2-2.067746E-02*S2*S2
        H2CO3T =  10.0**(0.0755*S2)
        HCO3T  =  10.0**DH1
        CO3T   =  10.0**DH2
        OH     =  HCO3T

        !**** Temperature adjustment

        KW = 10.0**(-5242.39/T1K+35.3944-8.350E-3*T1K-11.8261*LOG10(T1K))/OH
        K1 = 10.0**(-3404.71/T1K+14.8435-0.032786*T1K)*H2CO3T/HCO3T
        K2 = 10.0**(-2902.39/T1K+ 6.4980-0.023790*T1K)*HCO3T/CO3T

        !**** pH evaluation
        PHT = -pH - 2.1
        if (pH <= 0.0) PHT = -14.0
        INCR = 10.0

        do n=1,3
            F    = 1.0
            INCR = INCR/10.0
            ITER = 0
            
            do while (F > 0.0 .AND. ITER < 12)
                PHT    = PHT+INCR
                HION   = 10.0**PHT
                BICART = CART*K1*HION/(K1*HION+K1*K2+HION*HION)
                F      = BICART*(HION+2.0*K2)/HION+KW/HION-ALKT-HION/OH
                ITER   = ITER+1
            end do
            
            PHT = PHT-INCR
        enddo

        !**** pH, carbon dioxide, bicarbonate, and carbonate concentrations
        HION      =  10.0**PHT
        pH        = -PHT
        CO2       =  ICarbon /(1.0+K1/HION+K1*K2/(HION*HION))
        HCO3      =  ICarbon /(1.0+HION/K1+K2/HION)
        CO3       =  ICarbon /((HION*HION)/(K1*K2)+HION/K2+1.0)


        Me%ExternalVar%Mass(pHindex,index)   = pH
        Me%ExternalVar%Mass(CO2index,index)  = CO2
        Me%ExternalVar%Mass(CO3index,index)  = CO3
        Me%ExternalVar%Mass(HCO3index,index) = HCO3
        Me%SinksSources(pHindex)   = 0.0
        Me%SinksSources(CO2index)  = 0.0
        Me%SinksSources(CO3index)  = 0.0
        Me%SinksSources(HCO3index) = 0.0

    end subroutine ComputepH_CO2

   !-----------------------------------------------------------------------------

  
    subroutine CEQUALW2Benthic(CEQUALW2_ID,                                    &
                           Temperature,                                        &
                           Oxygen,                                             &
                           Mass,                                               &                           
                           OpenPoints,                                         &
                           SODRate,                                            &
                           STAT)  

        !Arguments---------------------------------------------------------------
        integer                                       :: CEQUALW2_ID
        real,                 pointer, dimension(:,:) :: Mass
        real   , optional,    pointer, dimension(:  ) :: SODRate
        real,                 pointer, dimension(:  ) :: Temperature
        real,                 pointer, dimension(:  ) :: Oxygen
        integer, optional,    pointer, dimension(:  ) :: OpenPoints
        integer, optional,    intent(OUT)             :: STAT
         
        !External----------------------------------------------------------------
        integer                                       :: index
        integer                                       :: ready_   
               
        !Local-------------------------------------------------------------------
        integer                                       :: STAT_         
        logical                                       :: CalcPoint
        !------------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(CEQUALW2_ID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then
                            
            Me%ExternalVar%Mass                         => Mass
            if (.NOT. associated(Me%ExternalVar%Mass))                  &
                stop 'CEQUALW2Benthic - ModuleCEQUALW2 - ERR01'
          
            Me%ExternalVar%Temperature                  => Temperature
            if (.NOT. associated(Me%ExternalVar%Temperature))           &
                stop 'CEQUALW2Benthic - ModuleCEQUALW2 - ERR02'
            
            Me%ExternalVar%Oxygen                       => Oxygen
            if (.NOT. associated(Me%ExternalVar%Oxygen))                &
                stop 'CEQUALW2Benthic - ModuleCEQUALW2 - ERR03'
                
            if (Me%SOD%UseSod .and. .NOT. present(SODRate))         then
                write(*,*) 'SOD defined in BentichCEQUALW2 but not in Interface Sediment Water.'
                stop 'CEQUALW2Benthic - ModuleCEQUALW2 - ERR04'
            
            else if (.NOT. Me%SOD%UseSod .and. present(SODRate))    then
                    write(*,*) 'SOD defined in InterfaceSedimentWater but not in BentichCEQUALW2 .'
                stop 'CEQUALW2Benthic - ModuleCEQUALW2 - ERR05'
            
            else if (Me%SOD%UseSod .and. present(SODRate))          then
                    Me%SOD%Rate => SODRate
            endif
             
            allocate (Me%SinksSources(Me%Size%PropLB:Me%Size%PropUB)) 
                          
do2 :       do index = Me%Size%ArrayLB, Me%Size%ArrayUB
            
                !If this module is called from the CEQUALW23D module, OpenPoint is present
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


cd3:            if (CalcPoint) then

                    call ComputeBenthicKineticRates           (index)
                    

                    !Phosphorus
                    if (Me%BenthicCompute%Phosphorus)then
                        call ComputeBenthicPhosphorus         (index)          
                    end if       
                    
                    
                    !Nitrogen
                    if (Me%BenthicCompute%Ammonia)then
                        call ComputeBenthicAmmonia            (index)  
                    endif

                    !Silica
                    if (Me%BenthicCompute%Dsilica)then
                        call ComputeBenthicDissolvedSilica    (index)   
                    endif
                                        
                                        
                    !Oxygen
                    if (Me%BenthicCompute%Oxygen)then
                        call ComputeBenthicDissolvedOxygen    (index)
                    endif

                    !Icarbon
                    if (Me%BenthicCompute%ICarbon)then 
                        call ComputeBenthicICarbon            (index)   
                    endif
   
                
                    if (Me%BenthicCompute%Detritus)then
                       call ComputeBenthicDetritus            (index)
                    endif
                end if cd3
            end do do2

            nullify(Me%ExternalVar%Salinity    )
            nullify(Me%ExternalVar%Temperature )
            nullify(Me%ExternalVar%Oxygen      )
            nullify(Me%ExternalVar%Alkalinity  )
            nullify(Me%ExternalVar%Mass        )
           
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine CEQUALW2Benthic

    
    !----------------------------------------------------------------------------
        

    subroutine ComputeBenthicKineticRates(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index

        !Local-----------------------------------------------------------------
        integer                                 :: DET
        real                                    :: temperature, Oxygen, Lam1


        !-----------------------------------------------------------------------
        
        DET         = Me%PropIndex%Detritus   
        
        ! temperature rate multiplier

        Temperature = Me%ExternalVar%Temperature (index)
        Oxygen      = Me%ExternalVar%Oxygen      (index)
   
        lam1        = Rising(Temperature,Me%DETT1,Me%DETT2,Me%DETK1,Me%DETK2)
        Me%DETTRM   = lam1/(1.0+lam1-Me%DETK1)
        
        !Auxiliary variables (related with oxygen limitation)for the calculation of decay rates                       
        Me%DO3      = (1.0+SIGN(1.0, Oxygen -1.E-10)) * 0.5
        
        Me%DETD     =  Me%DETTRM  * Me%SDK      *Me%ExternalVar%Mass(DET,index) * Me%DO3                         
        
        !Calculate SOD TRM if needed
        if (Me%SOD%UseSOD) then
            
            lam1        = Rising(Temperature,Me%SOD%T1,Me%SOD%T2,Me%SOD%K1,Me%SOD%K2)
            
            Me%SOD%TRM  = lam1/(1.0+lam1-Me%SOD%K1)
            
            Me%DO2      = (1.0+SIGN(1.0, Me%O2LIM  - Oxygen)) * 0.5
            
            Me%SODD     =  Me%SOD%TRM * Me%SOD%Rate(index) * Me%DO2
    
        endif                

                        
        !----------------------------------------------------------------------

    end subroutine ComputeBenthicKineticRates

    
    !--------------------------------------------------------------------------

    
    subroutine ComputeBenthicPhosphorus(index)
    
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index



        !Local-----------------------------------------------------------------
        integer                                 :: PHOSP
        real                                    :: PO4DET 

        !-----------------------------------------------------------------------

        PHOSP       = Me%PropIndex%Phosphorus
        
        PO4DET  = Me%DETD * Me%ORGP

        Me%ExternalVar%Mass(PHOSP,index) = Me%ExternalVar%Mass(PHOSP,index) + PO4DET * Me%DtDay
        
        
        if (Me%SOD%UseSOD )then
        
            Me%ExternalVar%Mass(PHOSP,index) = Me%ExternalVar%Mass(PHOSP,index) + Me%SODD * Me%SOD%PO4R
        
        endif
        
                            
    
    end subroutine ComputeBenthicPhosphorus

    
    !--------------------------------------------------------------------------
    
    
    subroutine ComputeBenthicAmmonia(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index


        !Local-----------------------------------------------------------------
        integer                                 :: AM
        real                                    :: NH4DET

        !-----------------------------------------------------------------------

        AM      = Me%PropIndex%Ammonia
      
        NH4DET  =  Me%DETD * Me%ORGN
        
        Me%ExternalVar%Mass(AM,index) = Me%ExternalVar%Mass(AM,index) + NH4DET * Me%DtDay
        
        if (Me%SOD%UseSOD )then
        
            Me%ExternalVar%Mass(AM,index) = Me%ExternalVar%Mass(AM,index) + Me%SODD * Me%SOD%NH4R
        
        endif

        
    end subroutine ComputeBenthicAmmonia

    
    !---------------------------------------------------------------------------
           
    
    subroutine ComputeBenthicDissolvedSilica(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index


        !Local-----------------------------------------------------------------
        integer                                     ::  DSI
        real                                        ::  DSIDDET
        !-----------------------------------------------------------------------
        
        DSI        = Me%PropIndex%sidiss 
        

        DSIDDET  =  Me%DETD * Me%ORGSI
        

        Me%ExternalVar%Mass(DSI,index) = Me%ExternalVar%Mass(DSI,index) + DSIDDET * Me%DtDay
        
        if (Me%SOD%UseSOD )then
        
            Me%ExternalVar%Mass(DSI,index) = Me%ExternalVar%Mass(DSI,index) + Me%SODD * Me%SOD%SiR
        
        endif
    
    end subroutine ComputeBenthicDissolvedSilica 

    
    !--------------------------------------------------------------------------
    

    subroutine ComputeBenthicDissolvedOxygen(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index


        !Local-----------------------------------------------------------------
        integer                                     :: O 
        real                                        :: DODET
        real                                        :: SODDL

        !-----------------------------------------------------------------------

        O                  = Me%PropIndex%Oxygen
        
        DODET  =  Me%DETD  * Me%O2OM
        
        Me%ExternalVar%Mass(O,index) = Me%ExternalVar%Mass(O,index) - DODET * Me%DtDay
   
        if (Me%SOD%UseSOD )then
            
            if (Me%SOD%DefaultO2)then
                !remove Relplace DO2 with DO3 in SODD (this is how it's done in the original cequal code)
                !Doesn't make mutch cense because SOD uses O2 even if it still doesn't release any nutrients
                SODDL = Me%SODD * Me%DO3 / Me%DO2                
            else
                SODDL = Me%SODD
            endif
            
            Me%ExternalVar%Mass(O,index) = Me%ExternalVar%Mass(O,index) - SODDL * Me%SOD%O2Sink
        
        endif                                                                       
 
    end subroutine ComputeBenthicDissolvedOxygen   

    !--------------------------------------------------------------------------

    subroutine ComputeBenthicICarbon(index)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index


        !Local-----------------------------------------------------------------
        real                                        :: ICarbonDET
        integer                                     :: ICarbonIndex

        !-----------------------------------------------------------------------

        ICarbonIndex = Me%PropIndex%ICarbon

        ICarbonDET =  Me%ORGC * Me%DETD
                                                              
        Me%ExternalVar%Mass(ICarbonIndex,index) = Me%ExternalVar%Mass(ICarbonIndex,index) - ICarbonDET * Me%DtDay        
        
        if (Me%SOD%UseSOD )then
        
            Me%ExternalVar%Mass(ICarbonIndex,index) = Me%ExternalVar%Mass(ICarbonIndex,index) + Me%SODD * Me%SOD%CO2R
        
        endif       
                

    end subroutine ComputeBenthicICarbon
   

    !------------------------------------------------------------------------


    subroutine ComputeBenthicDetritus(index)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: index


        !Local-----------------------------------------------------------------

        integer                                     :: DetIndex

        !-----------------------------------------------------------------------

        DetIndex = Me%PropIndex%Detritus
                                                              
        Me%ExternalVar%Mass(DetIndex,index) = Me%ExternalVar%Mass(DetIndex,index) - Me%DETD * Me%DtDay
        

    end subroutine ComputeBenthicDetritus
   

    !------------------------------------------------------------------------
    !This function calculates the Rising limb temperature rate
    
    real function Rising (TT,TT1,TT2,SK1,SK2)
        
        !Arguments--------------------------------------------------------- 
        real, intent(IN) ::  TT,TT1,TT2,SK1,SK2
        
        !Begin-------------------------------------------------------------
        
        Rising = SK1*EXP(LOG(SK2*(1.0-SK1)/(SK1*(1.0-SK2)))/(TT2-TT1)*(TT-TT1))        !TC 2/4/01
        

    end function Rising             !TC 2/4/01
  
    
    !------------------------------------------------------------------------
    !This function calculates the Falling limb temperature rate
    
    real function Falling (TT,TT3,TT4,SK3,SK4)
        
        !Arguments--------------------------------------------------------- 
        real, intent(IN) ::  TT,TT3,TT4,SK3,SK4
        
        !Begin-------------------------------------------------------------

        Falling = SK4*EXP(LOG(SK3*(1.0-SK4)/(SK4*(1.0-SK3)))/(TT4-TT3)*(TT4-TT))
    
    end function Falling   

    !------------------------------------------------------------------------
 
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillCEQUALW2(CEQUALW2_ID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: CEQUALW2_ID
        integer, optional, intent(OUT)      :: STAT


        !External----------------------------------------------------------------
        integer                             :: ready_              
        integer                             :: STAT_CALL, nUsers

        !Local-------------------------------------------------------------------
        integer                             :: STAT_ 
        type(T_Algae),      pointer         :: Algae
        type(T_Epiphyton),  pointer         :: Epiphyton

        !Begin-------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(CEQUALW2_ID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mCEQUALW2_,  Me%InstanceID)
  
            if (nUsers == 0) then

                Algae => Me%FirstAlgae

                do while(associated(Algae))

                    deallocate(Algae%NLim, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'subroutine Kill_CEQUALW2; module ModuleCEQUALW2. ERR05.'
                    nullify(Algae%NLim)

                    deallocate(Algae%PLim, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'subroutine Kill_CEQUALW2; module ModuleCEQUALW2. ERR05.'
                    nullify(Algae%PLim)

                    deallocate(Algae%SLim, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'subroutine Kill_CEQUALW2; module ModuleCEQUALW2. ERR05.'
                    nullify(Algae%SLim)


                    deallocate(Algae%LightLim, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'subroutine Kill_CEQUALW2; module ModuleCEQUALW2. ERR06.'
                    nullify(Algae%LightLim)


                    deallocate(Algae%OverallLim, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'subroutine Kill_CEQUALW2; module ModuleCEQUALW2. ERR08.'
                    nullify(Algae%OverallLim)

                    Algae => Algae%Next

                end do


                Epiphyton => Me%FirstEpiphyton

                do while(associated(Epiphyton))

                    deallocate(Epiphyton%NLim, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'subroutine Kill_CEQUALW2; module ModuleCEQUALW2. ERR012.'
                    nullify(Epiphyton%NLim)

                    deallocate(Epiphyton%LightLim, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'subroutine Kill_CEQUALW2; module ModuleCEQUALW2. ERR013.'
                    nullify(Epiphyton%LightLim)


                    deallocate(Epiphyton%SLim, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'subroutine Kill_CEQUALW2; module ModuleCEQUALW2. ERR014a.'
                    nullify(Epiphyton%SLim)


                    deallocate(Epiphyton%OverallLim, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'subroutine Kill_CEQUALW2; module ModuleCEQUALW2. ERR015.'
                    nullify(Epiphyton%OverallLim)

                    Epiphyton => Epiphyton%Next

                end do

                if (ME%Rate%Compute) then

                    deallocate(Me%Rate%Value, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'subroutine Kill_CEQUALW2; module ModuleCEQUALW2. ERR016.'
                    nullify(Me%Rate%Value)

                    deallocate(Me%Rate%Match, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'subroutine Kill_CEQUALW2; module ModuleCEQUALW2. ERR017.'
                    nullify(Me%Rate%Match)

                 endif

                call DeallocateInstance

                CEQUALW2_ID = 0

                STAT_ = SUCCESS_

            endif

        else 

            STAT_ = ready_


        end if cd1


        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine KillCEQUALW2




    !------------------------------------------------------------------------

    subroutine DeallocateInstance 

      
        !Local-----------------------------------------------------------------
        type (T_CEQUALW2), pointer           :: AuxObjCEQUALW2
        type (T_CEQUALW2), pointer           :: PreviousObjCEQUALW2

        !Updates pointers
        if (Me%InstanceID == FirstObjCEQUALW2%InstanceID) then
            FirstObjCEQUALW2 => FirstObjCEQUALW2%Next
        else
            PreviousObjCEQUALW2 => FirstObjCEQUALW2
            AuxObjCEQUALW2      => FirstObjCEQUALW2%Next
            do while (AuxObjCEQUALW2%InstanceID /= Me%InstanceID)
                PreviousObjCEQUALW2 => AuxObjCEQUALW2
                AuxObjCEQUALW2      => AuxObjCEQUALW2%Next
            enddo

            !Now update linked list
            PreviousObjCEQUALW2%Next => AuxObjCEQUALW2%Next

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

    subroutine Ready (CEQUALW2_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: CEQUALW2_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (CEQUALW2_ID > 0) then
            call LocateObjCEQUALW2 (CEQUALW2_ID)
            ready_ = VerifyReadLock (mCEQUALW2_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1
     
        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjCEQUALW2 (CEQUALW2_ID)

        !Arguments-------------------------------------------------------------
        integer                                     :: CEQUALW2_ID

        !Local-----------------------------------------------------------------

        Me => FirstObjCEQUALW2
        do while (associated (Me))
            if (Me%InstanceID == CEQUALW2_ID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                                        &
            stop 'ModuleCEQUALW2 - LocateObjCEQUALW2 - ERR01'

    end subroutine LocateObjCEQUALW2

    !--------------------------------------------------------------------------


end module ModuleCEQUALW2

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------


