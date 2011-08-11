!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Sediment Quality
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Pedro Galvão - v4.0
! DESCRIPTION   : Zero-dimensional model for primary production, nitrogen and carbon cycle
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
!Properties Units:
!dissolved                 - as from called - in MOHID Land mg/L
!adsorbed/particulated     - as from called - in MOHID Land mg/kgsoil
!gases:
!CO2 and CH4               - same as adsorbed/particulated properties
!N2                        - same as dissolved properties
!exceptions:
!oxygen                    - mol/L
!hydrogen                  - mol/L
!microorganisms population - #org/kgsoil
!wind                      - km/day
!
!OXYGEN   	                :  [0/1]       0      !Connect/Disconnect Oxygen computation 
!SOL_BACTERIA 	            :  [0/1]       0	  !Connect/Disconnect Solubilizing bacteria computation 
!
!DTSECONDS                  :  [s]       86400.   !dt to evaluate
!NO3_LIMIT                  :  [mg/L]      0.     !Minimum value for denitrification or maximum value for 
!                                                   !methane production in organic matter decay
!NEW_RATES                  : [0/1]        0      !Connect/Disconnect new rates formulation using maximum * factors
!IMOBILIZATION              : [0/1]        1      !Connect/Disconnect immobilization
!
!EXPLICIT                   : [0/1]        1      !Explicit computation for time dicretization


Module ModuleSedimentQuality
    use ModuleFunctions            
    use ModuleLUD                 
    use ModuleEnterData
    use ModuleGlobalData             

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartSedimentQuality
    private ::      AllocateInstance
    private ::          NullifyAllSubTypePointers
    private ::      SQReadData
    private ::          SedimentQualityOptions
    private ::          SQPropertyIndexNumber
    private ::          SQReadCalcOptions
    private ::          SQOptionsConsistencyVerif
    private ::          SQReadFileConstants
    private ::      AllocateVariables
    private ::          Add_PropRateFlux
    private ::          Add_EquaRateFlux

    public  ::      Construct_SQRateFlux

    !Selector
    public  :: GetDTSedimentQuality
    public  :: GetSQOptions
    public  :: GetSedimentQualitySize   
    public  :: GetPropIndex
    public  :: GetPropRateFlux   
    public  :: UngetPropRateFlux
                  
    
    !Modifier
    public  :: SedimentQuality 
    private ::      AnaerobioseCalculation 
    private ::      OxygenCalculation
    private ::      HydrogenCalculation
    private ::      PotentialRatesCalc
    private ::          CalcTterm       !function
    private ::          CalcTtermDeath
    private ::          CalcPopulation          
    private ::      LogicalImobilization 
    private ::      LogicalImobilization_P   !!!Lúcia   
    private ::      LogicalLimitation
    private ::      StartSedimentQualityIteration
    private ::      SedimentQualityCoefCalculation    
    private ::          SQOxygen
    private ::          SQCarbon
    private ::              SQLabilOrganicCarbon
    private ::              SQRefractOrganicCarbon   
    private ::              SQHeterotrophicC          
    private ::              SQAutotrophicC        
    private ::              SQAnaerobicC
    private ::          SQNitrogen
    private ::              SQAmmonia
    private ::              SQNitrate
    private ::              SQLabilOrganicNitrogen   
    private ::              SQRefractOrganicNitrogen 
    private ::              SQHeterotrophicN          
    private ::              SQAutotrophicN        
    private ::              SQAnaerobicN             
    private ::              SQNgas
    private ::      SystemResolution
   
    private ::      SedimentQualityRatesCalculation

    private ::     CalcOxygenTerm
    private ::     CalcpHTerm


    !Destructor
    public  ::  KillSedimentQuality                                                     
    private ::      DeallocateInstance


    !Management
    private ::      Ready
    private ::          LocateObjSedimentQuality
    

    !Types---------------------------------------------------------------------
    type       T_PropIndex
        integer :: HeterotrophicN                   = null_int
        integer :: HeterotrophicC                   = null_int
        integer :: AutotrophicN                     = null_int
        integer :: AutotrophicC                     = null_int
        integer :: AnaerobicN                       = null_int
        integer :: AnaerobicC                       = null_int
        integer :: Labil_OM_C                       = null_int
        integer :: Labil_OM_N                       = null_int
        integer :: RefractOM_C                      = null_int
        integer :: RefractOM_N                      = null_int
        integer :: Ammonia                          = null_int
        integer :: Nitrate                          = null_int
        integer :: Ngas                             = null_int
        integer :: Oxygen                           = null_int
        integer :: HeterotrophicP                   = null_int
        integer :: HeterotrophicPop                 = null_int  !#org/kgsoil
        integer :: AutotrophicP                     = null_int
        integer :: AutotrophicPop                   = null_int  !#org/kgsoil     
        integer :: AnaerobicP                       = null_int
        integer :: AnaerobicPop                     = null_int  !#org/kgsoil
        integer :: Labil_OM_P                       = null_int
        integer :: RefractOM_P                      = null_int
        integer :: Inorganic_P_soluble              = null_int
        integer :: Inorganic_P_fix                  = null_int  
        integer :: SolC                             = null_int 
        integer :: SolN                             = null_int
        integer :: SolP                             = null_int
        integer :: SolPop                           = null_int
        integer :: CO2                              = null_int
        integer :: Urea                             = null_int
        integer :: AmmoniaGas                       = null_int
        integer :: methane                          = null_int
        
    
    end type T_PropIndex
    type T_Files
        integer                                 :: AsciiUnit
    end type T_Files    

    type           T_PropRateFlux
        integer                              :: ID
        real, pointer, dimension(:)          :: Field
        type(T_PropRateFlux), pointer        :: next,prev
    end type T_PropRateFlux

   type           T_EquaRateFlux
        integer                              :: ID
        real                                 :: scalar = null_real
        logical                              :: TimeSerie
        type(T_EquaRateFlux), pointer        :: next,prev
        type(T_PropRateFlux), pointer        :: FirstPropRateFlux
        type(T_PropRateFlux), pointer        :: LastPropRateFlux
    end type T_EquaRateFlux

    type           T_ExtraRate
        integer                              :: ID
        real, pointer, dimension(:)          :: Field
    end type T_ExtraRate
   
    type       T_PropCalc
        logical :: Carbon       = OFF
        logical :: Nitrogen     = OFF
        logical :: Oxygen       = OFF
        logical :: Phosphorus   = OFF           !!! lucia
        logical :: Sol_Bacteria = OFF           !!!Lúcia
    end type T_PropCalc

    type       T_CalcMethod
        logical :: ExplicitMethod = OFF
        logical :: ImplicitMethod = OFF
        logical :: SemiImpMethod  = OFF
    end type T_CalcMethod

 
    type       T_External
        real, pointer, dimension(:  )       :: Temperature
        real, pointer, dimension(:  )       :: ThetaF         !water content
        real, pointer, dimension(:,:)       :: Mass
        real, pointer, dimension(:  )       :: DissolvedToParticulate 
!        real                                :: ParticleDensity
        real, pointer, dimension(:  )       :: SoilDryDensity
!        real                                :: DrySoilVolume
        real, pointer, dimension(:  )       :: Salinity
        real, pointer, dimension(:  )       :: pH
        real, pointer, dimension(:  )       :: Ionic
        real , pointer, dimension(:  )      :: Pai
        real , pointer, dimension(:  )      :: Wind
        real , pointer, dimension(:  )      :: Oxygen
    end type T_External



    type        T_Coeficients
        real                                :: ActivationE          = null_real     ![kcal/mol]Initial Activation energy
        real                                :: Acoef                = null_real     ![s.day-1.pop-1] Specific coeficient
        real                                :: Kp                   = null_real     ![L/mol] Salinity Coefficient
        real                                :: OptimumTemperature   = null_real     !Optimum temperature
        real                                :: Value                = null_real     !specific rate value     
        integer                             :: RateIndex            = null_int      !Rate Index number
        
        !New computation
        real                                :: ConcMinCarbon        = null_real     ![mg/kgsoil]Minimum conc. below, rate is max
        real                                :: ConcMinNitrate       = null_real
        real                                :: ConcMinAmmonia       = null_real
        real                                :: ConcMinPF            = null_real
        real                                :: OptimumpH            = null_real     ! optimum pH
        real                                :: ConcOptO2            = null_real     ! [mg/L] Optimum conc. above, rate is maximum
        real                                :: ConcOptCarbon        = null_real     ! [mg/L] Optimum conc. above, rate is maximum
    end type    T_Coeficients


    !all the specific rates that are explicitly calculated
    type       T_Rates
        type(T_Coeficients)                :: Labil_OM_C
        type(T_Coeficients)                :: RefractOM_C
        type(T_Coeficients)                :: AmmoniaToNitrate
        type(T_Coeficients)                :: AmmoniaImobilization         
        type(T_Coeficients)                :: NitrateToNgas
        type(T_Coeficients)                :: NitrateImobilization
        type(T_Coeficients)                :: PhosphorusImobilization       !!!!Lúcia
        type(T_Coeficients)                :: Solubilizing       !!!!!!
        type(T_Coeficients)                :: Heterotrophs
        type(T_Coeficients)                :: Autotrophs
        type(T_Coeficients)                :: Anaerobic
        type(T_Coeficients)                :: MethaneProduction
        type(T_Coeficients)                :: UreaHydrolysis        
        type(T_Coeficients)                :: Sol      !!!!Lúcia              
    end type T_Rates

    type        T_Constants
        real                                :: CNRatio          = null_real     !Microorganisms C/N ratio
        real                                :: CPRatio          = null_real     !Microorganisms C/P ratio       !!!Lúcia
        real                                :: CPopRatio        = null_real     !Carbon to # microorganisms ratio
        real                                :: Population       = null_real     !Population
        integer                             :: MicroIndex       = null_int      !Propertie number
        real                                :: EficiencyC       = null_real     !Microorganisms eficiency 
                                                                                !(ex: remainings goes off as CO2, NO3)    
        real                                :: EficiencyN       = null_real     !Only used for the anaerobic 
                                                                                !Microorganisms      
        real                                :: MinimumPop       = null_real     !Minimum possible population     
        logical                             :: LogicalMinumumPOP                !True or false
    end type    T_Constants    
   

    type        T_Microorganisms
        type(T_Constants        )          :: Heterotrophs
        type(T_Constants        )          :: Autotrophs
        type(T_Constants        )          :: Anaerobic 
        type(T_Constants        )          :: Sols                     !!!Lúcia    
    end type    T_Microorganisms

    
    type      T_SedimentQuality
        private

        integer                                             :: InstanceID
        type(T_Size1D    )                                  :: Prop   
        type(T_PropIndex )                                  :: PropIndex
        type(T_PropCalc  )                                  :: PropCalc
        type(T_CalcMethod)                                  :: CalcMethod

        type(T_External     )                               :: ExternalVar
        type(T_EquaRateFlux ),   pointer                    :: FirstEquaRateFlux
        type(T_EquaRateFlux ),   pointer                    :: LastEquaRateFlux

        type(T_Rates            )                           :: SpecificRates
        type(T_Microorganisms   )                           :: Microorganisms
        type(T_Files            )                           :: Files   
        logical                                             :: Imobilization        = OFF       !!! o nome tem de ser mudado
        logical                                             :: Imobilization_P      = OFF       !!!Lúcia
        logical                                             :: NLimitation          = OFF   
        logical                                             :: PLimitation          = OFF       !!!Lúcia
        integer                                             :: Select               = null_int  !!!Lúcia
        real                                                :: Partition            = null_real
        real                                                :: AnaerobicPartition   = null_real
        real                                                :: Solpartition        = null_real  !!!
        double precision,       pointer, dimension(:,:)     :: Matrix
        real,                   pointer, dimension(:  )     :: IndTerm
        real,                   pointer, dimension(:  )     :: NewMass          !Used with Explicit method
        double precision,       pointer, dimension(:,:)     :: OxygenTab
        
        real                                                :: DTDay                = null_real
        real                                                :: DTSecond             = null_real
        
        real                                                :: Aerobiose            = null_real
        real                                                :: Anaerobiose          = null_real
        real                                                :: LabiOM_CN_Ratio      = null_real
        real                                                :: LabilOM_CP_ratio     = null_real !!!Lúcia
        real                                                :: RefractOM_CN_ratio   = null_real
        real                                                :: RefractOM_CP_ratio   = null_real !!!Lúcia
        real                                                :: Oxygen               = null_real
        real                                                :: Hydrogen             = null_real
        real                                                :: Khn
        real                                                :: NO3limit
        logical                                             :: ChangeRates
        logical                                             :: ComputeImobilization
        logical                                             :: NewRates
        logical                                             :: OxygenForcing
        real                                                :: ADJ
        !Instance of Module_EnterData
        integer                                             :: ObjEnterData = 0

        !Instance of ModuleLUD
        integer                                             :: ObjLUD = 0

        !Collection of instances
        type(T_SedimentQuality), pointer                    :: Next
 
    end type T_SedimentQuality


    !Global Module Variables
    type (T_SedimentQuality), pointer                       :: FirstObjSedimentQuality
    type (T_SedimentQuality), pointer                       :: Me


    !--------------------------------------------------------------------------
    !Constants
    real,    parameter :: Boltzman                 = 1.383E-23      ![J.ºK-1]
    real,    parameter :: Planck                   = 6.63E-34       ![J.s]
    real,    parameter :: UnivGC                   = 1.99E-3        ![Kcal.mole-1.ºK-1]


    contains



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartSedimentQuality(SedimentQualityID, Filename, STAT)
        
        !Arguments-------------------------------------------------------------
        integer                             :: SedimentQualityID
        character(LEN = *)                  :: FileName 
        integer, optional, intent(OUT)      :: STAT     

        !External--------------------------------------------------------------
        integer :: STAT_CALL
        integer :: ready_         

        !Local-----------------------------------------------------------------
        integer :: STAT_

        !----------------------------------------------------------------------
        
        STAT_ = UNKNOWN_
        

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mSedimentQuality_)) then
            nullify (FirstObjSedimentQuality)
            call RegisterModule (mSedimentQuality_) 
        endif
        
        call Ready(SedimentQualityID, ready_)

cd0 :   if (ready_ .EQ. OFF_ERR_) then
            
            call AllocateInstance 
            call NullifyAllSubTypePointers


            !Associate EnterData
            call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine StartSedimentQuality; module ModuleSedimentQuality. ERR02.'
            
            call ConstructAsciiOutPut

            call SQReadData
            call AllocateVariables

             
cd1 :       if (.NOT. Me%CalcMethod%ExplicitMethod) then
                call StartLUD(Me%ObjLUD,                            &
                              Me%Prop%ILB,                          &
                              Me%Prop%IUB,                          &
                              Me%Prop%ILB,                          &
                              Me%Prop%IUB,                          &
                              STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'Subroutine StartSedimentQuality; ModuleSedimentQuality. ERR01.'
            end if cd1 


            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'Subroutine StartSedimentQuality; ModuleSedimentQuality. ERR02.'

            !Returns ID
            SedimentQualityID    = Me%InstanceID

            STAT_ = SUCCESS_
        else  cd0
            
            stop 'ModuleSedimentQuality - StartSedimentQuality - ERR03'

        end if cd0


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartSedimentQuality

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
                 FILE   = '..\res\SQ_Situation_'//trim(adjustl(Number))//'.log', &
                 STATUS = "REPLACE",                                      &
                 IOSTAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                exit do1
            else
                Counter = Counter + 1
            end if
        enddo do1

        write (Me%Files%AsciiUnit, FMT=*) ' Situation Aerobiose Anaerobiose '
    
    end subroutine ConstructAsciiOutPut
    
    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------
    subroutine AllocateInstance

        !Local-----------------------------------------------------------------
        type (T_SedimentQuality), pointer           :: NewObjSedimentQuality
        type (T_SedimentQuality), pointer           :: PreviousObjSedimentQuality


        !Allocates new instance
        allocate (NewObjSedimentQuality)
        nullify  (NewObjSedimentQuality%Next)


        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjSedimentQuality)) then
            FirstObjSedimentQuality    => NewObjSedimentQuality
            Me                         => NewObjSedimentQuality
        else
            PreviousObjSedimentQuality => FirstObjSedimentQuality
            Me                         => FirstObjSedimentQuality%Next
            do while (associated(Me))
                PreviousObjSedimentQuality  => Me
                Me                          => Me%Next
            enddo
            Me                              => NewObjSedimentQuality
            PreviousObjSedimentQuality%Next => NewObjSedimentQuality
        endif


        Me%InstanceID = RegisterNewInstance (mSEDIMENTQUALITY_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------

    Subroutine NullifyAllSubTypePointers

        !----------------------------------------------------------------------

        nullify(Me%ExternalVar%Temperature)
        nullify(Me%ExternalVar%Mass       )
        nullify(Me%ExternalVar%ThetaF     )

        nullify(Me%Matrix                 )
        nullify(Me%IndTerm                )
        nullify(Me%NewMass                )

        !----------------------------------------------------------------------

    end subroutine NullifyAllSubTypePointers

    !--------------------------------------------------------------------------   
   
   
    !----------------------------------------------------------------------------
    subroutine Construct_SQRateFlux(SedimentQualityID, ArrayLB, ArrayUB, STAT)


        !Arguments---------------------------------------------------------------
        integer                             :: SedimentQualityID
        type(T_EquaRateFlux     ), pointer  :: EquaRateFluxX
        integer, optional                   :: STAT

        !External----------------------------------------------------------------
        integer                             :: ArrayLB,ArrayUB

        !Local-------------------------------------------------------------------
        type (T_PropRateFlux), pointer      :: NewPropRateFlux
        type (T_EquaRateFlux), pointer      :: NewEquaRateFlux
        type (T_PropIndex   ), pointer      :: PIndex
        integer                             :: STAT_CALL, STAT_, ready_
        integer                             :: PropLB, PropUB
        logical, dimension(:), pointer      :: LogicalEqua
        integer                             :: equa,countequa

        !----------------------------------------------------------------------

        
        STAT_ = UNKNOWN_

        call Ready(SedimentQualityID, ready_)

cd0 :   if (ready_ .EQ. IDLE_ERR_) then
            
                  
        
            PIndex => Me%PropIndex
   
            PropUB = Me%Prop%IUB
            PropLB = Me%Prop%ILB


            allocate(LogicalEqua(PropLB:PropUB))
            LogicalEqua =.false.
        
            countequa=0
       
            !Oxygen is always computed 
            Logicalequa(PIndex%Oxygen) =.true.
            countequa = countequa + 1.

            !Microorganisms are always computed
            Logicalequa(PIndex%HeterotrophicPop) =.true.
            countequa = countequa + 1.

            Logicalequa(PIndex%AutotrophicPop) =.true.
            countequa = countequa + 1.
            
            Logicalequa(PIndex%AnaerobicPop) =.true.
            countequa = countequa + 1.            

            if (Me%PropCalc%Sol_Bacteria) then
                Logicalequa(PIndex%SolPop) =.true.
                countequa = countequa + 1.    
            endif
      
            if (Me%PropCalc%Nitrogen) then

                Logicalequa(PIndex%HeterotrophicN     )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%AutotrophicN       )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%AnaerobicN         )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%Labil_OM_N         )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%RefractOM_N        )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%Ammonia            )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%Nitrate            )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%Ngas               )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%Urea               )=.true.
                countequa = countequa + 1.

!                Logicalequa(PIndex%AmmoniaGas         )=.true.
!                countequa = countequa + 1.
      
                !!!!                
                if (Me%PropCalc%Sol_Bacteria) then

                    Logicalequa(PIndex%SolN           )=.true.
                    countequa = countequa + 1.
                
                endif
  
            endif

       
            if (Me%PropCalc%Carbon) then

                Logicalequa(PIndex%HeterotrophicC     )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%AutotrophicC       )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%AnaerobicC         )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%Labil_OM_C         )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%RefractOM_C        )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%CO2                )=.true.
                countequa = countequa + 1.

!                Logicalequa(PIndex%methane            )=.true.
!                countequa = countequa + 1.



!!!!!!
                if (Me%PropCalc%Sol_Bacteria) then

                    Logicalequa(PIndex%SolC               )=.true.
                    countequa = countequa + 1.
                
                endif
    
            endif

          !!! o mesmo para o fósforo
          
            if (Me%PropCalc%Phosphorus) then

                Logicalequa(PIndex%HeterotrophicP     )=.true.  !!!Lúcia
                countequa = countequa + 1.

                Logicalequa(PIndex%AutotrophicP       )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%AnaerobicP         )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%Labil_OM_P         )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%RefractOM_P        )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%Inorganic_P_soluble)=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%Inorganic_P_fix    )=.true.  !!!Lúcia
                countequa = countequa + 1.

!!!
                if (Me%PropCalc%Sol_Bacteria) then

                    Logicalequa(PIndex%SolP           )=.true.   !!!Lúcia
                    countequa = countequa + 1.
                
                endif
            endif

            if (countequa.ne.PropUB) stop 'SubRoutine Construct_SQRateFlux ModuleSedimentQuality - ERR01'
                       
                  do equa = PropLB,PropUB
              
                      if (Logicalequa(equa)) then
                         allocate (NewEquaRateFlux, STAT = STAT_CALL)            
                         if (STAT_CALL .NE. SUCCESS_)                                                       &
                         stop 'Subroutine Construct_SQRateFlux; Module ModuleSedimentQuality. ERR02.'
                    
                         nullify(NewEquaRateFlux%Prev,NewEquaRateFlux%Next)
                         nullify(NewEquaRateFlux%FirstPropRateFlux,NewEquaRateFlux%LastPropRateFlux)
                     
                         ! Add new Property to the SedimentQuality List 
                         Call Add_EquaRateFlux(NewEquaRateFlux)         
                     
                         NewEquaRateFlux%ID         = equa
                     
                         if (STAT_CALL .NE. SUCCESS_)                                                       &
                         stop 'Subroutine Construct_SQRateFlux; Module ModuleSedimentQuality. ERR03.' 

                                      
                      endif
                  enddo

                !to_change foram alocados rate fluxes suf
              EquaRateFluxX => Me%FirstEquaRateFlux  


do1:         do while (associated(EquaRateFluxX))
           
                 do equa = PropLB,PropUB
              
                      if (LogicalEqua(equa)) then
                         allocate (NewPropRateFlux, STAT = STAT_CALL)            
                         if (STAT_CALL .NE. SUCCESS_)       &
                         stop 'Subroutine Construct_SQRateFlux; Module ModuleSedimentQuality. ERR04.'
                    
                         nullify(NewPropRateFlux%field)
                         nullify(NewPropRateFlux%Prev,NewPropRateFlux%Next)
                         allocate(NewPropRateFlux%Field   (ArrayLB:ArrayUB))
                     
                         ! Add new Prop  
                         Call Add_PropRateFlux(EquaRateFluxX, NewPropRateFlux)         
                     
                         NewPropRateFlux%ID         = equa
                         NewPropRateFlux%Field      = 0.
                     
                         if (STAT_CALL .NE. SUCCESS_)       &
                         stop 'Subroutine Construct_SQRateFlux; Module ModuleSedimentQuality. ERR05.' 

                                      
                      endif
                  enddo
        
              EquaRateFluxX => EquaRateFluxX%Next
          
              enddo  do1
        
              nullify(EquaRateFluxX)
            
           
            STAT_ = SUCCESS_
        
        else cd0               
        
            STAT_ = ready_
        
        end if cd0


        if (present(STAT))  STAT = STAT_          

        !------------------------------------------------------------------------

    end subroutine Construct_SQRateFlux

    !----------------------------------------------------------------------------

    !--------------------------------------------------------------------------

      ! This subroutine adds a new property rateflux to the Rate Fluxes List  

    subroutine Add_EquaRateFlux(NewEquaRateFlux)

        !Arguments-------------------------------------------------------------
        type(T_EquaRateFlux   ), pointer        :: NewEquaRateFlux

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstEquaRateFlux)) then
        
            Me%FirstEquaRateFlux      => NewEquaRateFlux
            Me%LastEquaRateFlux       => NewEquaRateFlux
        
        else
            
            NewEquaRateFlux%Prev      => Me%LastEquaRateFlux
            Me%LastEquaRateFlux%Next  => NewEquaRateFlux
            Me%LastEquaRateFlux       => NewEquaRateFlux
        
        end if 

        !----------------------------------------------------------------------

    end subroutine Add_EquaRateFlux 

    !--------------------------------------------------------------------------

    
    !--------------------------------------------------------------------------

      ! This subroutine adds a new property rateflux to the Rate Fluxes List  

    subroutine Add_PropRateFlux(EquaRateFluxX, NewPropRateFlux)

        !Arguments-------------------------------------------------------------
        type(T_EquaRateFlux),     pointer :: EquaRateFluxX
        type(T_PropRateFlux    ), pointer :: NewPropRateFlux

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(EquaRateFluxX%FirstPropRateFlux)) then
        
            EquaRateFluxX%FirstPropRateFlux        => NewPropRateFlux
            EquaRateFluxX%LastPropRateFlux         => NewPropRateFlux
        
        else
            
            NewPropRateFlux%Prev                   => EquaRateFluxX%LastPropRateFlux
            EquaRateFluxX%LastPropRateFlux%Next    => NewPropRateFlux
            EquaRateFluxX%LastPropRateFlux         => NewPropRateFlux
        
        end if 

        !----------------------------------------------------------------------

    end subroutine Add_PropRateFlux 

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine SedimentQualityOptions

        !External--------------------------------------------------------------
        integer                     :: FromFile
        integer                     :: STAT_CALL
        integer                     :: flag
        !Begin-----------------------------------------------------------------
        
        call GetExtractType (FromFile = FromFile)


        call GetData(Me%PropCalc%Nitrogen                   ,           &
                     Me%ObjEnterData, flag                  ,           &
                     SearchType   = FromFile                ,           &
                     keyword      = 'NITROGEN'              ,           &
                     default      = .true.                  ,           & 
                     ClientModule = 'ModuleSedimentQuality' ,           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'Subroutine SedimentQualityOptions; Module ModuleSedimentQuality. ERR01.' 


        call GetData(Me%PropCalc%Carbon                     ,           &
                     Me%ObjEnterData, flag                  ,           &
                     SearchType   = FromFile                ,           &
                     keyword      = 'CARBON'                ,           &
                     default      = .true.                  ,           & 
                     ClientModule = 'ModuleSedimentQuality' ,           &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'Subroutine SedimentQualityOptions; Module ModuleSedimentQuality. ERR02.'

!!!!Opção válida para o fosforo tb!!!
        call GetData(Me%PropCalc%Phosphorus                 ,           &    !!!Lúcia
                     Me%ObjEnterData, flag                  ,           &    !!!Lúcia
                     SearchType   = FromFile                ,           &    !!!Lúcia
                     keyword      = 'PHOSPHORUS'            ,           &    !!!Lúcia
                     default      = .true.                  ,           &    !!!Lúcia
                     ClientModule = 'ModuleSedimentQuality' ,           &    !!!Lúcia
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'Subroutine SedimentQualityOptions; Module ModuleSedimentQuality. ERR03.'

!!!!Opção válida para as bacterias solubilizadoras tb!!!
        call GetData(Me%PropCalc%Sol_Bacteria               ,           &    !!!Lúcia
                     Me%ObjEnterData, flag                  ,           &    !!!Lúcia
                     SearchType   = FromFile                ,           &    !!!Lúcia
                     keyword      = 'SOL_BACTERIA'          ,           &    !!!Lúcia
                     default      = .false.                  ,          &    !!!Lúcia
                     ClientModule = 'ModuleSedimentQuality' ,           &    !!!Lúcia
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'Subroutine SedimentQualityOptions; Module ModuleSedimentQuality. ERR04.'


        call GetData(Me%NO3limit                            ,           &    
                     Me%ObjEnterData, flag                  ,           &    
                     SearchType   = FromFile                ,           &    
                     keyword      = 'NO3_LIMIT'             ,           &    
                     default      = 0.                      ,           &    
                     ClientModule = 'ModuleSedimentQuality' ,           &    
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'Subroutine SedimentQualityOptions; Module ModuleSedimentQuality. ERR05.'


        call GetData(Me%ChangeRates                         ,           &    
                     Me%ObjEnterData, flag                  ,           &    
                     SearchType   = FromFile                ,           &    
                     keyword      = 'CHANGE_RATES'          ,           &    
                     default      = .false.                 ,           &    
                     ClientModule = 'ModuleSedimentQuality' ,           &    
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'Subroutine SedimentQualityOptions; Module ModuleSedimentQuality. ERR06.'

        call GetData(Me%NewRates                            ,           &    
                     Me%ObjEnterData, flag                  ,           &    
                     SearchType   = FromFile                ,           &    
                     keyword      = 'NEW_RATES'             ,           &    
                     default      = .false.                 ,           &    
                     ClientModule = 'ModuleSedimentQuality' ,           &    
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'Subroutine SedimentQualityOptions; Module ModuleSedimentQuality. ERR07.'

        call GetData(Me%ComputeImobilization                ,           &    
                     Me%ObjEnterData, flag                  ,           &    
                     SearchType   = FromFile                ,           &    
                     keyword      = 'IMOBILIZATION'         ,           &    
                     default      = .true.                  ,           &    
                     ClientModule = 'ModuleSedimentQuality' ,           &    
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'Subroutine SedimentQualityOptions; Module ModuleSedimentQuality. ERR08.'

        !----------------------------------------------------------------------

    end subroutine SedimentQualityOptions         
    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------
    subroutine SQPropertyIndexNumber

        !----------------------------------------------------------------------

        Me%Prop%ILB = 1.
        Me%Prop%IUB = 0

        !Carbon index number
        if (Me%PropCalc%Carbon) then
            Me%Prop%IUB                                 = Me%Prop%IUB + 1
            Me%PropIndex%HeterotrophicC                 = Me%Prop%IUB
            Me%Microorganisms%Heterotrophs%MicroIndex   = Me%Prop%IUB

            Me%Prop%IUB                                 = Me%Prop%IUB + 1
            Me%PropIndex%AutotrophicC                   = Me%Prop%IUB
            Me%Microorganisms%Autotrophs%MicroIndex     = Me%Prop%IUB

            Me%Prop%IUB                                 = Me%Prop%IUB + 1
            Me%PropIndex%AnaerobicC                     = Me%Prop%IUB
            Me%Microorganisms%Anaerobic%MicroIndex      = Me%Prop%IUB    

            Me%Prop%IUB                                 = Me%Prop%IUB + 1
            Me%PropIndex%Labil_OM_C                     = Me%Prop%IUB

            Me%Prop%IUB                                 = Me%Prop%IUB + 1
            Me%PropIndex%RefractOM_C                    = Me%Prop%IUB
  !!!!          
            Me%Prop%IUB                                 = Me%Prop%IUB + 1
            Me%PropIndex%CO2                            = Me%Prop%IUB

!            Me%Prop%IUB                                 = Me%Prop%IUB + 1
!            Me%PropIndex%methane                        = Me%Prop%IUB

            if (Me%PropCalc%Sol_Bacteria) then
                Me%Prop%IUB                                 = Me%Prop%IUB + 1
                Me%PropIndex%SolC                           = Me%Prop%IUB
                Me%Microorganisms%Sols%MicroIndex           = Me%Prop%IUB
            endif

          
        endif  

        !Nitrogen index number
        if (Me%PropCalc%Nitrogen) then
            Me%Prop%IUB                     = Me%Prop%IUB + 1
            Me%PropIndex%HeterotrophicN     = Me%Prop%IUB

            Me%Prop%IUB                     = Me%Prop%IUB + 1
            Me%PropIndex%AutotrophicN       = Me%Prop%IUB

            Me%Prop%IUB                     = Me%Prop%IUB + 1
            Me%PropIndex%AnaerobicN         = Me%Prop%IUB

            Me%Prop%IUB                     = Me%Prop%IUB + 1
            Me%PropIndex%Labil_OM_N         = Me%Prop%IUB

            Me%Prop%IUB                     = Me%Prop%IUB + 1
            Me%PropIndex%RefractOM_N        = Me%Prop%IUB

            Me%Prop%IUB                     = Me%Prop%IUB + 1
            Me%PropIndex%Ammonia            = Me%Prop%IUB
        
            Me%Prop%IUB                     = Me%Prop%IUB + 1
            Me%PropIndex%Nitrate            = Me%Prop%IUB

            Me%Prop%IUB                     = Me%Prop%IUB + 1
            Me%PropIndex%Ngas               = Me%Prop%IUB
        
            Me%Prop%IUB                     = Me%Prop%IUB + 1
            Me%PropIndex%Urea               = Me%Prop%IUB
 
!            Me%Prop%IUB                     = Me%Prop%IUB + 1
!            Me%PropIndex%AmmoniaGas         = Me%Prop%IUB

 !!!!
            if (Me%PropCalc%Sol_Bacteria) then
 
                Me%Prop%IUB                     = Me%Prop%IUB + 1
                Me%PropIndex%SolN               = Me%Prop%IUB

            endif
               
        endif    


        !Oxygen index number -> The oxygen is always calculated 
        Me%Prop%IUB         = Me%Prop%IUB + 1
        Me%PropIndex%Oxygen = Me%Prop%IUB
        
        !Micro Populations are always computed
        Me%Prop%IUB                   = Me%Prop%IUB + 1
        Me%PropIndex%HeterotrophicPop = Me%Prop%IUB
        
        Me%Prop%IUB                   = Me%Prop%IUB + 1
        Me%PropIndex%AutotrophicPop   = Me%Prop%IUB
       
        Me%Prop%IUB                   = Me%Prop%IUB + 1
        Me%PropIndex%AnaerobicPop     = Me%Prop%IUB

        if (Me%PropCalc%Sol_Bacteria) then
            Me%Prop%IUB               = Me%Prop%IUB + 1
            Me%PropIndex%SolPop       = Me%Prop%IUB
       endif

        ! o mesmo para o fósforo

        if (Me%PropCalc%Phosphorus) then
            Me%Prop%IUB                      = Me%Prop%IUB + 1  !!!Lúcia
            Me%PropIndex%HeterotrophicP      = Me%Prop%IUB

            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%AutotrophicP        = Me%Prop%IUB

            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%AnaerobicP          = Me%Prop%IUB

            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%Labil_OM_P          = Me%Prop%IUB

            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%RefractOM_P         = Me%Prop%IUB

            Me%Prop%IUB                      = Me%Prop%IUB + 1
            Me%PropIndex%Inorganic_P_soluble = Me%Prop%IUB

            Me%Prop%IUB                      = Me%Prop%IUB + 1  !!!Lúcia
            Me%PropIndex%Inorganic_P_fix     = Me%Prop%IUB

    !!!!!!
            if (Me%PropCalc%Sol_Bacteria) then

                Me%Prop%IUB                      = Me%Prop%IUB + 1
                Me%PropIndex%SolP                = Me%Prop%IUB
            
            endif

        endif

        !----------------------------------------------------------------------
    end subroutine SQPropertyIndexNumber

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine SQReadData

       
        !Local-----------------------------------------------------------------
        logical :: Consistent

        !----------------------------------------------------------------------


        call SedimentQualityOptions
        call SQPropertyIndexNumber
        call SQReadCalcOptions 


        Consistent = SQOptionsConsistencyVerif ()


        if (Consistent) then
            call SQReadFileConstants
        else
            write(*,*) 
            write(*,*) 'The SedimentQuality Options were not consistent, verify file data.'
            stop 'Subroutine SQReadData; ModuleSedimentQuality. ERR01.'
        endif   !Consistent


        !----------------------------------------------------------------------

    end subroutine SQReadData

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine SQReadCalcOptions

        !External--------------------------------------------------------------
        integer                     :: FromFile
        integer                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                     :: flag
        
        !----------------------------------------------------------------------
        
        call GetExtractType (FromFile = FromFile)

        !Verifica se se pretende calcular usando um metodo EXPLICITO
        call GetData(Me%CalcMethod%ExplicitMethod                   ,   &
                     Me%ObjEnterData, flag                          ,   &
                     SearchType   = FromFile                        ,   &
                     keyword      = 'EXPLICIT'                      ,   &
                     ClientModule = 'ModuleSedimentQuality'         ,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'Subroutine SQReadCalcOptions; Module ModuleSedimentQuality. ERR01.'
        
        !Verifica se se pretende calcular usando um metodo IMPLICITO/EXPLICITO        
        call GetData(Me%CalcMethod%SemiImpMethod                    ,   &
                     Me%ObjEnterData, flag                          ,   &
                     SearchType   = FromFile                        ,   &
                     keyword      = 'SEMIIMP'                       ,   &
                     ClientModule = 'ModuleSedimentQuality'         ,   &
                     STAT         = STAT_CALL)  
        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'Subroutine SQReadCalcOptions; Module ModuleSedimentQuality. ERR02.' 
        
        !It has to be inspected why implicit and semi implicit methods are giving
        !different results from explicit. Check the hipotesis that Matrix in Water quality has the 
        !signals different and that has to be accounted when calling ModuleLUD
        if (Me%CalcMethod%SemiImpMethod) then
            write (*,*) 'For now only explicit computation is possible'
            write (*,*) 'Use EXPLICIT : 1'
            stop 'Subroutine SQReadCalcOptions; Module ModuleSedimentQuality. ERR020.' 
        endif
        
        !Verifica se se pretende calcular usando um metodo IMPLICITO
        call GetData(Me%CalcMethod%ImplicitMethod                   ,   &
                     Me%ObjEnterData, flag                          ,   &
                     SearchType   = FromFile                        ,   &
                     keyword      = 'IMPLICIT'                      ,   &
                     ClientModule = 'ModuleSedimentQuality'         ,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'Subroutine SQReadCalcOptions; Module ModuleSedimentQuality. ERR03.' 

        !It has to be inspected why implicit and semi implicit methods are giving
        !different results from explicit. Check the hipotesis that Matrix in Water quality has the 
        !signals different and that has to be accounted when calling ModuleLUD
        if (Me%CalcMethod%ImplicitMethod) then
            write (*,*) 'For now only explicit computation is possible'
            write (*,*) 'Use EXPLICIT : 1'
            stop 'Subroutine SQReadCalcOptions; Module ModuleSedimentQuality. ERR030.' 
        endif

        !----------------------------------------------------------------------

    end subroutine SQReadCalcOptions

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    logical function SQOptionsConsistencyVerif ()

        !Local-----------------------------------------------------------------

        integer :: aux

        !----------------------------------------------------------------------

            aux = 0

            if (Me%CalcMethod%ExplicitMethod) aux = aux + 1.
            if (Me%CalcMethod%ImplicitMethod) aux = aux + 1.
            if (Me%CalcMethod%SemiImpMethod ) aux = aux + 1.


cd1 :       if (aux .EQ. 1.) then
                SQOptionsConsistencyVerif = .TRUE.

            else 
                SQOptionsConsistencyVerif = .FALSE.
            end if cd1

        !----------------------------------------------------------------------

    end function SQOptionsConsistencyVerif

    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------
    subroutine SQReadFileConstants

        !External--------------------------------------------------------------
        
        integer :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                         :: FromFile
        
        integer                         :: Flag
        integer                         :: blockpointer

        type(T_Coeficients) , pointer   :: RateinProgress
        type(T_Constants)   , pointer   :: MicroorginProgress

        character(LEN = StringLength)   :: block_begin 
        character(LEN = StringLength)   :: block_end

        integer                                     :: ClientNumber
        integer                                     :: status
        logical                                     :: BlockFound

        !--------------------------------------------------------------------------
        
        call GetExtractType (FromFile = FromFile)
        
        call GetData(Me%DTSecond                             ,   &
                     Me%ObjEnterData, flag                   ,   &
                     SearchType   = FromFile                 ,   &
                     keyword      ='DTSECONDS'               ,   & 
                     default      = 3600.                    ,   & 
                     ClientModule = 'ModuleSedimentQuality'  ,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                             &
            stop 'Subroutine SQReadFileConstants; Module ModuleSedimentQuality. ERR01.' 

cd1:   if (flag .EQ. 0) then
            write(*,*) 
            write(*,*) 'Keyword DTSECONDS not found in Water quality data file.'
            write(*,*) 'Subroutine SQReadFileConstants; Module ModuleSedimentQuality. WRN01.'
            write(*,*) 'Assumed ', Me%DTSecond, 'seconds (',  Me%DTSecond / 60.0, 'hour).'
            write(*,*) 
        end if cd1
        
        !For compatibility with the rest of the program,  
        Me%DTDay = Me%DTSecond / 24.0 / 60.0 / 60.0


        !Reads non specific rates & constants--------------------------------------
        blockpointer = 0

do1 :   do
            blockpointer = blockpointer + 1.

            Call SelectRateBlock(RateinProgress, block_begin, block_end, blockpointer)

            if (blockpointer .EQ. null_int) exit
            
            !be able to start from beggining of file so that rate order is not important
            call RewindBuffer (Me%ObjEnterData, STAT = status)
            if (status /= SUCCESS_) &
            call SetError(FATAL_, INTERNAL_, 'SQReadFileConstants - ModuleSedimentQuality - ERR010')
                
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,      &
                                        block_begin, block_end, BlockFound,                 &
                                        STAT = status)

cd2 :       if      (status .EQ. SUCCESS_      ) then    
cd3 :           if (BlockFound) then                                                  
                    ! Construct a New Rate 
                    Call ConstructRate(RateinProgress) 
                
                end if cd3

            else if (status .EQ. BLOCK_END_ERR_) then cd2
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                call SetError(FATAL_, INTERNAL_, 'SQReadFileConstants - ModuleSedimentQuality - ERR020') 
            end if cd2
        end do do1

    
    !Read Microorganisms constants--------------------------------------
        blockpointer = 0

do2 :   do
            blockpointer = blockpointer + 1.

            Call SelectMicroBlock(MicroorginProgress, block_begin, block_end, blockpointer)
            if (blockpointer .EQ. null_int) exit
                
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,              &
                                        block_begin, block_end, BlockFound,         &
                                        STAT = status)

cd4 :       if      (status .EQ. SUCCESS_      ) then    
cd5 :           if (BlockFound) then                                                  
                    ! Construct a Micororganisms 
                    Call ConstructMicroorg(MicroorginProgress)

                end if cd5

            else if (status .EQ. BLOCK_END_ERR_) then cd4
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                call SetError(FATAL_, INTERNAL_, 'SQReadFileConstants - ModuleSedimentQuality - ERR030') 
            end if cd4
        end do do2

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = status) 
        if (status /= SUCCESS_) &
            call SetError(FATAL_, INTERNAL_, 'SQReadFileConstants - ModuleSedimentQuality - ERR040')

    end subroutine SQReadFileConstants
    !--------------------------------------------------------------------------


    !----------------------------------------------------------------------------
    subroutine SelectRateBlock (RateinProgress, block_begin, block_end, blockpointer)

        !Arguments---------------------------------------------------------------
        
        type(T_Coeficients), pointer    :: RateinProgress
        
        character(LEN = StringLength)   :: block_begin 
        character(LEN = StringLength)   :: block_end
        
        integer                         :: blockpointer
        !------------------------------------------------------------------------
        

        if (Me%PropCalc%Carbon) then

            Select case (blockpointer)
        
            Case (1)
                RateinProgress  =>  Me%SpecificRates%Labil_OM_C
                block_begin     =   '<begin_Labil_OM_C_Rate>'
                block_end       =   '<end_Labil_OM_C_Rate>'
                return 

            Case (2) 
                RateinProgress  =>  Me%SpecificRates%RefractOM_C
                block_begin     =   '<begin_RefractOM_C_Rate>'
                block_end       =   '<end_RefractOM_C_Rate>'
                return

            Case (3)
                RateinProgress  =>  Me%SpecificRates%Heterotrophs
                block_begin     =   '<begin_Heterotrophs_Rate>'
                block_end       =   '<end_Heterotrophs_Rate>'
                return

            Case (4)
                RateinProgress  =>  Me%SpecificRates%Autotrophs
                block_begin     =   '<begin_Autotrophs_Rate>'
                block_end       =   '<end_Autotrophs_Rate>'
                return

            Case (5)
                RateinProgress  =>  Me%SpecificRates%Anaerobic 
                block_begin     =   '<begin_Anaerobic_Rate>'
                block_end       =   '<end_Anaerobic_Rate>'
                return

            Case (6)
                RateinProgress  =>  Me%SpecificRates%MethaneProduction 
                block_begin     =   '<begin_methane_production>'
                block_end       =   '<end_methane_production>'
                return



!!!!!!!!!!!!                    
            end select
    end if
            
        
        if (Me%PropCalc%Nitrogen) then

            Select case (blockpointer)

            Case (7)
                RateinProgress  =>  Me%SpecificRates%AmmoniaToNitrate
                block_begin     =   '<begin_AmmoniaToNitrate_Rate>'
                block_end       =   '<end_AmmoniaToNitrate_Rate>'
                return

            Case (8)
                RateinProgress  =>  Me%SpecificRates%AmmoniaImobilization 
                block_begin     =   '<begin_AmmoniaImobilization_Rate>'
                block_end       =   '<end_AmmoniaImobilization_Rate>'
                return

            Case (9)
                RateinProgress  =>  Me%SpecificRates%NitrateToNgas
                block_begin     =   '<begin_NitrateToNgas_Rate>'
                block_end       =   '<end_NitrateToNgas_Rate>'
                return

            Case (10)        
                RateinProgress  =>  Me%SpecificRates%NitrateImobilization               
                block_begin     =   '<begin_NitrateImobilization_Rate>'
                block_end       =   '<end_NitrateImobilization_Rate>'
                return

            Case (11)        
                RateinProgress  =>  Me%SpecificRates%UreaHydrolysis               
                block_begin     =   '<begin_Urea_Hydrolysis>'
                block_end       =   '<end_Urea_Hydrolysis>'
                return
            end select
         
        end if


        !!!!!!O mesmo para o fosforo... e esta em construcçção!!!
        
        if (Me%PropCalc%Phosphorus) then                    !!!Lúcia

            Select case (blockpointer)                      !!!Lúcia

            Case (12)        
                RateinProgress  =>  Me%SpecificRates%PhosphorusImobilization !!!Lúcia              
                block_begin     =   '<begin_PhosphorusImobilization_Rate>'   !!!Lúcia
                block_end       =   '<end_PhosphorusImobilization_Rate>'
                return

            end Select  
         
        end if
                    
        if (Me%PropCalc%Sol_Bacteria) then

            Select case(blockpointer)
            
            Case (13)
                RateinProgress  =>  Me%SpecificRates%Sol 
                block_begin     =   '<begin_Sol_Rate>'
                block_end       =   '<end_Sol_Rate>'
                return
               
            Case (14)        
                RateinProgress  =>  Me%SpecificRates%Solubilizing 
                block_begin     =   '<begin_Solubilizing_Rate>'
                block_end       =   '<end_Solubilizing_Rate>'
                return

            end select
        
        end if  
   
        blockpointer = null_int 
        !------------------------------------------------------------------------

    end subroutine SelectRateBlock   
    !----------------------------------------------------------------------------  


    !----------------------------------------------------------------------------
    subroutine SelectMicroBlock(MicroorginProgress, block_begin, block_end, blockpointer)

        !Arguments---------------------------------------------------------------
        
        type(T_Constants)   , pointer   :: MicroorginProgress
        
        character(LEN = StringLength)   :: block_begin 
        character(LEN = StringLength)   :: block_end
        
        integer                         :: blockpointer
        !------------------------------------------------------------------------
        
        if (Me%PropCalc%Carbon) then

            Select case (blockpointer)
        
            Case (1)
                MicroorginProgress  =>  Me%Microorganisms%Heterotrophs
                block_begin         =   '<begin_Heterotrophs>'
                block_end           =   '<end_Heterotrophs>'
                return 

            Case (2) 
                MicroorginProgress  => Me%Microorganisms%Autotrophs
                block_begin         =   '<begin_Autotrophs>'
                block_end           =   '<end_Autotrophs>'
                return

            Case (3)
                MicroorginProgress  =>  Me%Microorganisms%Anaerobic
                block_begin         =   '<begin_Anaerobic>'
                block_end           =   '<end_Anaerobic>'
                return          

            end select
        end if    
        
        if (Me%PropCalc%Sol_Bacteria) then

            Select case (blockpointer)
    
            Case (4)
                MicroorginProgress  =>  Me%Microorganisms%Sols
                block_begin         =   '<begin_Sols>'
                block_end           =   '<end_Sols>'
                return
            
          end select

        end If
        
        blockpointer = null_int 
        !------------------------------------------------------------------------

    end subroutine SelectMicroBlock   
    !---------------------------------------------------------------------------- 


    !----------------------------------------------------------------------------
    subroutine ConstructRate(RateinProgress)

        !Arguments---------------------------------------------------------------
        
        type(T_Coeficients) , pointer   :: RateinProgress

        !External----------------------------------------------------------------
        integer                                     :: status

        !Local-------------------------------------------------------------------
        integer :: FromBlock
        real    :: value
        integer :: iflag
        !------------------------------------------------------------------------
        
        call GetExtractType(FromBlock = FromBlock)
        
        call GetData(   value                                    ,   & 
                        Me%ObjEnterData, iflag                   ,   &
                        SearchType      = FromBlock              ,   &
                        Keyword         = 'AE'                   ,   &
                        ClientModule    = 'ModuleSedimentQuality',   & 
                        default         = 0.                     ,   &
                        STAT            = status)
       
        if (status /= SUCCESS_)                                      &
           call SetError(FATAL_, INTERNAL_, 'ConstructRate - ModuleSedimentQuality - ERR01') 
        
        RateinProgress%ActivationE = value
       
        call GetData(   value                                    ,   &
                        Me%ObjEnterData, iflag                   ,   &
                        SearchType      = FromBlock              ,   &
                        keyword         ='Acoef'                 ,   &
                        ClientModule    = 'ModuleSedimentQuality',   &
                        default         = 0.                     ,   &
                        STAT           = status)
       
        if (status /= SUCCESS_)                                      &
           call SetError(FATAL_, INTERNAL_, 'ConstructRate - ModuleSedimentQuality - ERR02') 
        
        RateinProgress%Acoef = value

        call GetData(   value                                    ,   &
                        Me%ObjEnterData, iflag                   ,   &
                        SearchType      = FromBlock              ,   &
                        keyword         ='kp'                 ,   &
                        ClientModule    = 'ModuleSedimentQuality',   &
                        default         = 0.                     ,   &
                        STAT           = status)
       
        if (status /= SUCCESS_)                                      &
           call SetError(FATAL_, INTERNAL_, 'ConstructRate - ModuleSedimentQuality - ERR02') 
        
        RateinProgress%kp = value


        call GetData(   value                                    ,   & 
                        Me%ObjEnterData, iflag                   ,   &
                        SearchType      = FromBlock              ,   &
                        keyword         ='Temperature'           ,   &
                        ClientModule    = 'ModuleSedimentQuality',   &
                        default         = 0.                     ,   &
                        STAT            = status)
       
        if (status /= SUCCESS_)                                      &
           call SetError(FATAL_, INTERNAL_, 'ConstructRate - ModuleSedimentQuality - ERR03') 
        
        RateinProgress%OptimumTemperature = value


        call GetData(   value                                    ,   & 
                        Me%ObjEnterData, iflag                   ,   &
                        SearchType      = FromBlock              ,   &
                        keyword         ='OptimumpH'             ,   &
                        ClientModule    = 'ModuleSedimentQuality',   &
                        default         = 7.                     ,   &
                        STAT            = status)
       
        if (status /= SUCCESS_)                                      &
           call SetError(FATAL_, INTERNAL_, 'ConstructRate - ModuleSedimentQuality - ERR04') 
        
        RateinProgress%OptimumpH = value
        
        !For death rates
        call GetData(   value                                    ,   & 
                        Me%ObjEnterData, iflag                   ,   &
                        SearchType      = FromBlock              ,   &
                        keyword         ='ConcMinCarbon'         ,   &
                        ClientModule    = 'ModuleSedimentQuality',   &
                        default         = 0.                     ,   &
                        STAT            = status)
       
        if (status /= SUCCESS_)                                      &
           call SetError(FATAL_, INTERNAL_, 'ConstructRate - ModuleSedimentQuality - ERR05') 
        
        RateinProgress%ConcMinCarbon = value

        call GetData(   value                                    ,   & 
                        Me%ObjEnterData, iflag                   ,   &
                        SearchType      = FromBlock              ,   &
                        keyword         ='ConcMinNitrate'       ,   &
                        ClientModule    = 'ModuleSedimentQuality',   &
                        default         = 0.                     ,   &
                        STAT            = status)
       
        if (status /= SUCCESS_)                                      &
           call SetError(FATAL_, INTERNAL_, 'ConstructRate - ModuleSedimentQuality - ERR06') 
        
        RateinProgress%ConcMinNitrate = value

        call GetData(   value                                    ,   & 
                        Me%ObjEnterData, iflag                   ,   &
                        SearchType      = FromBlock              ,   &
                        keyword         ='ConcMinAmmonia'        ,   &
                        ClientModule    = 'ModuleSedimentQuality',   &
                        default         = 0.                     ,   &
                        STAT            = status)
       
        if (status /= SUCCESS_)                                      &
           call SetError(FATAL_, INTERNAL_, 'ConstructRate - ModuleSedimentQuality - ERR07') 
        
        RateinProgress%ConcMinAmmonia = value

        call GetData(   value                                    ,   & 
                        Me%ObjEnterData, iflag                   ,   &
                        SearchType      = FromBlock              ,   &
                        keyword         ='ConcMinPF'             ,   &
                        ClientModule    = 'ModuleSedimentQuality',   &
                        default         = 0.                     ,   &
                        STAT            = status)
       
        if (status /= SUCCESS_)                                      &
           call SetError(FATAL_, INTERNAL_, 'ConstructRate - ModuleSedimentQuality - ERR08') 
        
        RateinProgress%ConcMinPF = value

        !For denitrification
        call GetData(   value                                    ,   & 
                        Me%ObjEnterData, iflag                   ,   &
                        SearchType      = FromBlock              ,   &
                        keyword         ='ConcOptCarbon'         ,   &
                        ClientModule    = 'ModuleSedimentQuality',   &
                        default         = 0.                     ,   &
                        STAT            = status)
       
        if (status /= SUCCESS_)                                      &
           call SetError(FATAL_, INTERNAL_, 'ConstructRate - ModuleSedimentQuality - ERR09') 
        
        RateinProgress%ConcOptCarbon = value

        call GetData(   value                                    ,   & 
                        Me%ObjEnterData, iflag                   ,   &
                        SearchType      = FromBlock              ,   &
                        keyword         ='ConcOptO2'             ,   &
                        ClientModule    = 'ModuleSedimentQuality',   &
                        default         = 0.                     ,   &
                        STAT            = status)
       
        if (status /= SUCCESS_)                                      &
           call SetError(FATAL_, INTERNAL_, 'ConstructRate - ModuleSedimentQuality - ERR010') 
        
        RateinProgress%ConcOptO2 = value

        !------------------------------------------------------------------------

    end subroutine ConstructRate   
    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------
    subroutine ConstructMicroorg(MicroorginProgress   )

        !Arguments---------------------------------------------------------------
        
        type(T_Constants)       , pointer   :: MicroorginProgress

        !External----------------------------------------------------------------
        
        integer                                     :: status
        
        !Local-------------------------------------------------------------------
        
        integer :: FromBlock
        integer :: iflag
        !------------------------------------------------------------------------
        
        call GetExtractType(FromBlock = FromBlock)
        
        call GetData(   MicroorginProgress%CNRatio                  ,   & 
                        Me%ObjEnterData, iflag                      ,   &
                        SearchType      = FromBlock                 ,   &
                        Keyword         = 'CN_RATIO'                ,   &
                        ClientModule    = 'ModuleSedimentQuality'   ,   & 
                        default         = 0.                        ,   &
                        STAT            = status)
        if (status /= SUCCESS_)                                         &
           call SetError(FATAL_, INTERNAL_, 'ConstructMicroorg - ModuleSedimentQuality - ERR01') 
        

        !!! inserir a razão CP dos microorganismos                          
        call GetData(   MicroorginProgress%CPRatio                  ,   & !!!Lúcia
                        Me%ObjEnterData, iflag                      ,   &
                        SearchType      = FromBlock                 ,   &
                        Keyword         = 'CP_RATIO'                ,   &
                        ClientModule    = 'ModuleSedimentQuality'   ,   & 
                        default         = 0.                        ,   & !!!Lúcia
                        STAT            = status)
        if (status /= SUCCESS_)                                         &
           call SetError(FATAL_, INTERNAL_, 'ConstructMicroorg - ModuleSedimentQuality - ERR01') 



       call GetData(   MicroorginProgress%CPopRatio                 ,   & 
                        Me%ObjEnterData, iflag                      ,   &
                        SearchType      = FromBlock                 ,   &
                        Keyword         = 'POPULATION_CARBON_RATIO' ,   &
                        ClientModule    = 'ModuleSedimentQuality'   ,   & 
                        default         = 0.                        ,   &
                        STAT            = status)       
        if (status /= SUCCESS_)                                         &
           call SetError(FATAL_, INTERNAL_, 'ConstructMicroorg - ModuleSedimentQuality - ERR02') 
       
       
         call GetData(   MicroorginProgress%EficiencyC              ,   & 
                        Me%ObjEnterData, iflag                      ,   &
                        SearchType      = FromBlock                 ,   &
                        Keyword         = 'CARBON_EFICIENCY'        ,   &
                        ClientModule    = 'ModuleSedimentQuality'   ,   & 
                        default         = 0.                        ,   &
                        STAT            = status)
        if (status /= SUCCESS_)                                         &
           call SetError(FATAL_, INTERNAL_, 'ConstructMicroorg - ModuleSedimentQuality - ERR03')
           
        
        call GetData(   MicroorginProgress%EficiencyN               ,   & 
                        Me%ObjEnterData, iflag                      ,   &
                        SearchType      = FromBlock                 ,   &
                        Keyword         = 'NITROGEN_EFICIENCY'      ,   &
                        ClientModule    = 'ModuleSedimentQuality'   ,   & 
                        default         = 0.                        ,   &
                        STAT            = status)
        if (status /= SUCCESS_)                                         &
           call SetError(FATAL_, INTERNAL_, 'ConstructMicroorg - ModuleSedimentQuality - ERR04')  
           
       
        call GetData(   MicroorginProgress%MinimumPop               ,   & 
                        Me%ObjEnterData, iflag                      ,   &
                        SearchType      = FromBlock                 ,   &
                        Keyword         = 'MINIMUM_POPULATION'      ,   &
                        ClientModule    = 'ModuleSedimentQuality'   ,   & 
                        default         = 0.                        ,   &
                        STAT            = status)
        if (status /= SUCCESS_)                                         &
           call SetError(FATAL_, INTERNAL_, 'ConstructMicroorg - ModuleSedimentQuality - ERR05')  
                 
       
        !------------------------------------------------------------------------

    end subroutine ConstructMicroorg   

    !----------------------------------------------------------------------------


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
                stop 'Subroutine AllocateVariables; module ModuleSedimentQuality. ERR01.'
        end if cd1

        !------------------------------------------------------------------------

    end subroutine AllocateVariables   
    !----------------------------------------------------------------------------



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------
    subroutine GetSedimentQualitySize(SedimentQualityID, PropLB, PropUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: SedimentQualityID
        integer, optional, intent(OUT)      :: PropLB,    PropUB
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentQualityID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PropLB   )) PropLB    = Me%Prop%ILB
            if (present(PropUB   )) PropUB    = Me%Prop%IUB

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetSedimentQualitySize
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine GetSQOptions(SedimentQualityID,                 &
                            Nitrogen,                          &
                            Oxygen,                            &
                            Carbon,                            &
                            Phosphorus,                        &  !!!
                            Sol_Bacteria,                      &  !!!Lúcia
                            ExplicitMethod,                    &
                            ImplicitMethod,                    &
                            SemiImpMethod, STAT) 

        !Arguments-------------------------------------------------------------
        integer                        :: SedimentQualityID
        integer, optional, intent(OUT) :: STAT
        logical, optional, intent(OUT) :: Nitrogen, Oxygen, Carbon , Phosphorus,Sol_Bacteria        !!!Lúcia
        logical, optional, intent(OUT) :: ExplicitMethod, ImplicitMethod, SemiImpMethod  


        !External--------------------------------------------------------------

        integer :: ready_              

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentQualityID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.     &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(Nitrogen      )) Nitrogen       = Me%PropCalc%Nitrogen
            if (present(Oxygen        )) Oxygen         = Me%PropCalc%Oxygen
            if (present(Carbon        )) Carbon         = Me%PropCalc%Carbon    
            if (present(Phosphorus    )) Phosphorus     = Me%PropCalc%Phosphorus    !!!Lúcia                
            if (present(Sol_Bacteria  )) Sol_Bacteria   = Me%PropCalc%Sol_Bacteria  !!!
            if (present(ExplicitMethod)) ExplicitMethod = Me%CalcMethod%ExplicitMethod
            if (present(ImplicitMethod)) ImplicitMethod = Me%CalcMethod%ImplicitMethod
            if (present(SemiImpMethod )) SemiImpMethod  = Me%CalcMethod%SemiImpMethod    

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetSQOptions
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine GetPropIndex(SedimentQualityID,  HeterotrophicN,              &
                                                HeterotrophicC,              &
                                                AutotrophicN,            &
                                                AutotrophicC,            &
                                                AnaerobicN,                 &
                                                AnaerobicC,                 &
                                                Labil_OM_C,                 &
                                                Labil_OM_N,                 &
                                                RefractOM_C,                &
                                                RefractOM_N,                &
                                                Ammonia,                    &
                                                Nitrate,                    &
                                                Ngas,                       &
                                                Oxygen,                     &
                                                HeterotrophicP,              &  !!!Lúcia
                                                AutotrophicP,               &
                                                AnaerobicP,                 &
                                                Labil_OM_P,                 &
                                                RefractOM_P,                &
                                                Inorganic_P_soluble,        &
                                                Inorganic_P_fix,            &
                                                SolC,                       &
                                                SolN,                       &
                                                SolP,                       &
                                                CO2,                        &
                                                Urea,                       &
!                                                Ammoniagas,                 &
!                                                Methane,                    &
                                                HeterotrophicPop,           &
                                                AutotrophicPop,             &
                                                AnaerobicPop,               &
                                                SolPop,                     &
                                                STAT)


        !Arguments-------------------------------------------------------------
        integer                             :: SedimentQualityID

        integer, optional, intent(OUT)      :: HeterotrophicN
        integer, optional, intent(OUT)      :: HeterotrophicC
        integer, optional, intent(OUT)      :: HeterotrophicP    !!!Lúcia
        integer, optional, intent(OUT)      :: HeterotrophicPop

        integer, optional, intent(OUT)      :: AutotrophicN
        integer, optional, intent(OUT)      :: AutotrophicC
        integer, optional, intent(OUT)      :: AutotrophicP !!!Lúcia
        integer, optional, intent(OUT)      :: AutotrophicPop

    
        integer, optional, intent(OUT)      :: AnaerobicN
        integer, optional, intent(OUT)      :: AnaerobicC
        integer, optional, intent(OUT)      :: AnaerobicP   !!!Lúcia
        integer, optional, intent(OUT)      :: AnaerobicPop
        integer, optional, intent(OUT)      :: SolPop
        
        integer, optional, intent(OUT)      :: Labil_OM_C
        integer, optional, intent(OUT)      :: Labil_OM_N
        integer, optional, intent(OUT)      :: Labil_OM_P   !!!Lúcia



        integer, optional, intent(OUT)      :: RefractOM_C
        integer, optional, intent(OUT)      :: RefractOM_N
        integer, optional, intent(OUT)      :: RefractOM_P  !!!Lúcia
    
    
        
        integer, optional, intent(OUT)      :: Ammonia
        integer, optional, intent(OUT)      :: Nitrate
        integer, optional, intent(OUT)      :: Ngas
        integer, optional, intent(OUT)      :: Urea
!        integer, optional, intent(OUT)      :: AmmoniaGas
!        integer, optional, intent(OUT)      :: Methane
        integer, optional, intent(OUT)      :: Inorganic_P_soluble  !!!Lúcia
        integer, optional, intent(OUT)      :: CO2
        integer, optional, intent(OUT)      :: Inorganic_P_fix      !!!Lúcia

        integer, optional, intent(OUT)      :: SolC      !!!Lúcia
        integer, optional, intent(OUT)      :: SolN      !!!Lúcia
        integer, optional, intent(OUT)      :: SolP      !!!Lúcia
        

        integer, optional, intent(OUT) :: Oxygen
        
        integer, optional, intent(OUT) :: STAT
 
        !External--------------------------------------------------------------
        integer :: ready_              

        !Local-----------------------------------------------------------------
        integer :: STAT_              !Auxiliar local variable
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentQualityID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                     &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            if (present(HeterotrophicN       )) HeterotrophicN        = Me%PropIndex%HeterotrophicN
            if (present(HeterotrophicC       )) HeterotrophicC        = Me%PropIndex%HeterotrophicC
            if (present(HeterotrophicP       )) HeterotrophicP        = Me%PropIndex%HeterotrophicP    !!!Lúcia
            if (present(HeterotrophicPop     )) HeterotrophicPop      = Me%PropIndex%HeterotrophicPop    

            if (present(AutotrophicN        )) AutotrophicN         = Me%PropIndex%AutotrophicN
            if (present(AutotrophicC        )) AutotrophicC         = Me%PropIndex%AutotrophicC
            if (present(AutotrophicP        )) AutotrophicP         = Me%PropIndex%AutotrophicP !!!Lúcia
            if (present(AutotrophicPop      )) AutotrophicPop       = Me%PropIndex%AutotrophicPop

            if (present(AnaerobicN          )) AnaerobicN           = Me%PropIndex%AnaerobicN
            if (present(AnaerobicC          )) AnaerobicC           = Me%PropIndex%AnaerobicC
            if (present(AnaerobicP          )) AnaerobicP           = Me%PropIndex%AnaerobicP   !!!Lúcia
            if (present(AnaerobicPop        )) AnaerobicPop         = Me%PropIndex%AnaerobicPop
            
            if (present(Labil_OM_C          )) Labil_OM_C           = Me%PropIndex%Labil_OM_C
            if (present(Labil_OM_N          )) Labil_OM_N           = Me%PropIndex%Labil_OM_N
            if (present(Labil_OM_P          )) Labil_OM_P           = Me%PropIndex%Labil_OM_P   !!!Lúcia
            
            
            if (present(Labil_OM_C          )) RefractOM_C          = Me%PropIndex%RefractOM_C
            if (present(Labil_OM_N          )) RefractOM_N          = Me%PropIndex%RefractOM_N
            if (present(Labil_OM_P          )) RefractOM_P          = Me%PropIndex%RefractOM_P  !!!Lúcia
            

            if (present(Ammonia             )) Ammonia              = Me%PropIndex%Ammonia
            if (present(Nitrate             )) Nitrate              = Me%PropIndex%Nitrate
            if (present(Ngas                )) Ngas                 = Me%PropIndex%Ngas
            if (present(CO2                 )) CO2                  = Me%PropIndex%CO2
            if (present(Urea                )) Urea                 = Me%PropIndex%Urea
!            if (present(AmmoniaGas          )) AmmoniaGas            = Me%PropIndex%AmmoniaGas
!            if (present(methane             )) Methane              = Me%PropIndex%Methane
            
            if (present(Inorganic_P_soluble )) Inorganic_P_soluble  = Me%PropIndex%Inorganic_P_soluble  !!!Lúcia
            if (present(Inorganic_P_fix     )) Inorganic_P_fix      = Me%PropIndex%Inorganic_P_fix  !!!Lúcia
   
            if (present(SolC                )) SolC                 = Me%PropIndex%SolC !!!Lúcia
            if (present(SolN                )) SolN                 = Me%PropIndex%SolN !!!Lúcia
            if (present(SolP                )) SolP                 = Me%PropIndex%SolP !!!Lúcia
            if (present(SolPop              )) SolPop               = Me%PropIndex%SolPop

            if (present(Oxygen              )) Oxygen               = Me%PropIndex%Oxygen    

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
        !----------------------------------------------------------------------

    end subroutine GetPropIndex
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine GetDTSedimentQuality(SedimentQualityID, DTDay, DTSecond, STAT)

        !Arguments-------------------------------------------------------------
        integer                        :: SedimentQualityID
        real,    optional, intent(OUT) :: DTDay
        real,    optional, intent(OUT) :: DTSecond
        integer, optional, intent(OUT) :: STAT

        !External--------------------------------------------------------------
        integer :: ready_              

        !Local-----------------------------------------------------------------
        integer                          :: STAT_              !Auxiliar local variable
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentQualityID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.       &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(DTDay   )) DTDay    = Me%DTDay
            if (present(DTSecond)) DTSecond = Me%DTSecond

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDTSedimentQuality
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine GetPropRateFlux(SedimentQualityID, Firstprop, Secondprop, PropRateFlux, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: SedimentQualityID

        real, dimension(:), pointer         :: PropRateFlux

        integer                             :: Firstprop,Secondprop
        integer, optional, intent(OUT)      :: STAT
                
        type(T_PropRateFlux    ), pointer   :: PropRateFluxX
        type(T_EquaRateFlux    ), pointer   :: EquaRateFluxX

        !External--------------------------------------------------------------
        integer :: ready_        

        !Local-----------------------------------------------------------------
        integer                            :: STAT_              !Auxiliar local variable
        integer                            :: prop,equa
        logical                            :: found       

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentQualityID, ready_)
        
        found=.FALSE.
       

cd1 :   if ( (ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_) ) then
                
                call Read_Lock(mSEDIMENTQUALITY_, Me%InstanceID)


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

                nullify  (PropRateFluxX,EquaRateFluxX)
           
            if (found) then              
              STAT_ = SUCCESS_
            else
              STAT_ = NOT_FOUND_ERR_
            endif
        
        else  cd1
        
            STAT_ = ready_
        
        end if cd1


        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetPropRateFlux
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine UngetPropRateFlux(SedimentQualityID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: SedimentQualityID
        integer, optional, intent(OUT)      :: STAT
        real, pointer, dimension(:)         :: Array

        !External--------------------------------------------------------------
        integer :: ready_   

        !Local-----------------------------------------------------------------
        integer :: STAT_            
        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(SedimentQualityID, ready_)

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            
            call Read_UnLock(mSEDIMENTQUALITY_, Me%InstanceID, "UngetPropRateFlux")
            
            nullify(Array)

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetPropRateFlux
    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    !--------------------------------------------------------------------------
    subroutine SedimentQuality (    SedimentQualityID       ,                   &
                                    Temperature             ,                   &
                                    Mass                    ,                   &
                                    ThetaF                  ,                   &
                                    DissolvedToParticulate  ,                   &
                                    ArrayLB, ArrayUB        ,                   &
!                                    particledensity         ,                   &
                                    SoilDryDensity          ,                   &
!                                    DrySoilVolume           ,                   &
                                    Salinity                ,                   &
                                    pH                      ,                   &
                                    Ionic                   ,                   &
                                    Pai                     ,                   &
                                    Wind                    ,                   &
                                    Oxygen                  ,                   &
                                    OpenPoints              ,                   &
                                    STAT)  

        !Arguments---------------------------------------------------------------
        integer                                     :: SedimentQualityID
        real,               pointer, dimension(:  ) :: Temperature
        real,               pointer, dimension(:  ) :: ThetaF
        real,               pointer, dimension(:  ) :: DissolvedToParticulate
        real,               pointer, dimension(:,:) :: Mass
        integer, optional,  pointer, dimension(:  ) :: OpenPoints
        integer,            intent(IN )             :: ArrayLB, ArrayUB
!        real,               intent (IN)             :: ParticleDensity
!        real,               intent (IN)             :: DrySoilVolume
        real,   pointer, dimension(:  ), optional   :: SoilDryDensity
        real,               pointer, dimension(:  ) :: Salinity
        real,               pointer, dimension(:  ) :: pH
        real,               pointer, dimension(:  ) :: Ionic
        real,   pointer, dimension(:  ), optional   :: Pai
        real,   pointer, dimension(:  ), optional   :: Wind
        real,   pointer, dimension(:  ), optional   :: Oxygen               
        integer, optional,  intent(OUT)             :: STAT
        !External----------------------------------------------------------------
        integer                                     :: index
        integer                                     :: ready_   
!        integer                                     :: conta      
        !Local-------------------------------------------------------------------
        integer                                     :: STAT_          
        logical                                     :: CalcPoint

        !------------------------------------------------------------------------                         
            
        STAT_ = UNKNOWN_


        call Ready(SedimentQualityID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            Me%ExternalVar%Temperature => Temperature
            if (.NOT. associated(Me%ExternalVar%Temperature) )                  &
                stop 'Subroutine SedimentQuality; Module ModuleSedimentQuality. ERR01' 

            Me%ExternalVar%Mass => Mass
            if (.NOT. associated(Me%ExternalVar%Mass) )                         &
                stop 'Subroutine SedimentQuality; Module ModuleSedimentQuality. ERR02.'

            Me%ExternalVar%ThetaF => ThetaF !to_change
            if (.NOT. associated(Me%ExternalVar%ThetaF) )                       &
                stop 'Subroutine SedimentQuality; Module ModuleSedimentQuality. ERR03.'

            Me%ExternalVar%DissolvedToParticulate => DissolvedToParticulate 
            if (.NOT. associated(Me%ExternalVar%DissolvedToParticulate) )       &
                stop 'Subroutine SedimentQuality; Module ModuleSedimentQuality. ERR04.'
            
            if (present(SoilDryDensity)) Me%ExternalVar%SoilDryDensity => SoilDryDensity
            if (.NOT. associated(Me%ExternalVar%SoilDryDensity) )                       &
                stop 'Subroutine SedimentQuality; Module ModuleSedimentQuality. ERR05.'
            
            if (Me%PropCalc%Phosphorus) then

                if (present(Pai))            Me%ExternalVar%Pai            => Pai
                if (.NOT. associated(Me%ExternalVar%SoilDryDensity) )                       &
                    stop 'Subroutine SedimentQuality; Module ModuleSedimentQuality. ERR06.'
            endif


            if (Me%PropCalc%Nitrogen) then
                if (present(Wind))           Me%ExternalVar%Wind           => Wind
                if (.NOT. associated(Me%ExternalVar%SoilDryDensity) )                       &
                    stop 'Subroutine SedimentQuality; Module ModuleSedimentQuality. ERR07.'
            endif

            !Salinity fo oxygen computation
            Me%ExternalVar%Salinity       => Salinity
            if (.NOT. associated(Me%ExternalVar%SoilDryDensity) )                       &
                stop 'Subroutine SedimentQuality; Module ModuleSedimentQuality. ERR08.'
            
            Me%ExternalVar%pH             => pH
            if (.NOT. associated(Me%ExternalVar%SoilDryDensity) )                       &
                stop 'Subroutine SedimentQuality; Module ModuleSedimentQuality. ERR09.'
            
            Me%ExternalVar%Ionic          => Ionic
            if (.NOT. associated(Me%ExternalVar%SoilDryDensity) )                       &
                stop 'Subroutine SedimentQuality; Module ModuleSedimentQuality. ERR010.'

            
            if (present(Oxygen)) then
                Me%OxygenForcing = .true.
                Me%ExternalVar%Oxygen     => Oxygen
            endif


            call StartSedimentQualityIteration !zeros the matrix and independent term



do1 :       do index = ArrayLB, ArrayUB
            
                !If this module is called from the Interface module, OpenPoint is present
                !and the  module runs for all Openpoints
                !If this module is called from the Lagrangian module, OpenPoint is not present
                !and the  module runs for all volumes
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
                
                
                    call OxygenCalculation              (index           )
                
                    call AnaerobioseCalculation         (ThetaF, index   )
                
                    call HydrogenCalculation            (index           )        
                    
                    if (.not. Me%NewRates) then
                        call PotentialRatesCalc         (index           )
                    else
                        call PotentialRatesCalc_New     (index           )
                    endif
                    
                    call SelectSituation                (index           )   !!!Lúcia

                    call SedimentQualityCoefCalculation (index           )
                    
                    !The rates can just be calculated if the rate flux is associated
                    !In the case that this module is used by the lagrangian module
                    !the rate fluxes are not calculated
                    !Rates must be computed before call to WQSystemResolution to use 
                    !old concentrations 
                    
                    if (associated(Me%FirstEquaRateFlux)) then
                        call SedimentQualityRatesCalculation   (index)
                    endif                
                    
                    call SystemResolution               (index           )
                
                    call CorrectOM                      (index           )

                   
                    call LogSituation
                    

                endif  
              
!                conta                   = conta +1
                Me%Matrix    (:, :)     = 0.
                Me%Imobilization        = OFF
                Me%Imobilization_P      = OFF !!!Lúcia
                Me%NLimitation          = OFF
                ME%PLimitation          = OFF !!!Lúcia
                Me%Select               = 0.  !!!Lúcia
                Me%AnaerobicPartition   = 0.
                Me%Partition            = 0. 
                   

            end do do1


            nullify(Me%ExternalVar%Temperature   )
            nullify(Me%ExternalVar%Mass          )
            nullify(Me%ExternalVar%ThetaF        )
                     
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1

!        Write(*,*) 'Factor Aerobiose' , Me%Aerobiose
!        Write(*,*) 'Factor Anaerobiose' , Me%Anaerobiose



        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine SedimentQuality
    !----------------------------------------------------------------------------


    subroutine LogSituation

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------
                
       
        write (Me%Files%AsciiUnit, fmt=1000) Me%Select, Me%Aerobiose, Me%Anaerobiose

        1000 format(i3, f12.5, f12.5)

    end subroutine LogSituation

    !--------------------------------------------------------------------------


    ! Calculation of the aerobic and anaerobic functions for Microbiological action.
    ! When a soil atmosphere algoritm is implemented this will no longer be necessary
    !----------------------------------------------------------------------------
    subroutine AnaerobioseCalculation (ThetaF, index)

        !Arguments---------------------------------------------------------------

        integer             , intent(IN)                    :: index  
        real                , dimension(:  ), intent(IN)    :: ThetaF
        !Local-------------------------------------------------------------------

        real                                                :: PWFPC
!        real                                                :: ParticleDensity
!        real                                                :: SoilDensity
        !real                                                :: conversion
      !------------------------------------------------------------------------
        
!        ParticleDensity = Me%ExternalVar%ParticleDensity

!        SoilDensity = Me%ExternalVar%SoilDryDensity 

        PWFPC = ThetaF(index)*100
                        
                
! coeficiente de anaerobiose

        if (PWFPC <60) then 

            Me%Anaerobiose = 0.00001            ! porque o SedimteQuality nao esta preparado para receber FANA =0

        else

            Me%Anaerobiose = 0.000304*EXP(0.0815*PWFPC)

        end if


! coeficiente de aerobiose


        if (PWFPC<=20) then

            Me%Aerobiose = 0.75

        else if (PWFPC>20 .AND. PWFPC<59) then

            Me%Aerobiose = -0.253+0.0203*PWFPC


        else

            Me%Aerobiose = 41.1*EXP(-0.0625*PWFPC)


        end If




! vou por aqui o NO3 limit, so para não fazer confusao
    
!        Conversion = Me%ExternalVar%DissolvedToParticulate (index)
!
!        Me%NO3limit = 0.00000000000000005*1000/conversion
       


!       Me%aerobiose = 0.7164           !versao tste comparação RZWQM

!       Me%anaerobiose = 5.97E-2

        !------------------------------------------------------------------------

    end subroutine AnaerobioseCalculation
    !----------------------------------------------------------------------------


    !Reads or Computes oxygen concentration (mol/L) depending of the temperature and salinity 
    !Metcalf and Eddy (1978).
    !----------------------------------------------------------------------------
    subroutine OxygenCalculation (index)

        !Arguments---------------------------------------------------------------

        integer                 , intent(IN)    :: index  

        !Local-------------------------------------------------------------------
        integer :: IOXI
        real    :: Temp
        real    :: oxy
        real    :: sal

        !------------------------------------------------------------------------
       
       ! values in Mol/L
       
        !If forcing then read the passed values, else, compute from temperature and salinity
        if (Me%OxygenForcing) then
        
            oxy = Me%ExternalVar%Oxygen(index)
        else
         
            Temp                        = Me%ExternalVar%Temperature(index)
            IOXI                        = Me%PropIndex%OXYGEN   
            sal                         = Me%ExternalVar%Salinity(index)               

            ! values in Mol/L

            if (sal == 0) Then

                oxy =  -2E-09*temp**3 + 2E-07*temp**2 - 1E-05*temp + 0.0005

            else if (sal == 5 .OR. sal == 10 .OR. sal == 15 .OR. sal == 20 &
                        .OR. sal == 25 .OR. sal == 30 .OR. sal == 35) then


                oxy = -1E-09*temp**3 + 2E-07*temp**2 - 1E-05*temp + 0.0004

            else if  (sal == 40) then

                oxy =  -1E-09*temp**3 + 2E-07*temp**2 - 9E-06*temp + 0.0003

            else

                oxy = -1E-09*temp**3 + 1E-07*temp**2 - 8E-06*temp + 0.0003

            end if
        endif
        
        Me%Oxygen= oxy
        Me%ExternalVar%Mass(Me%PropIndex%Oxygen,index) = oxy

   
    end subroutine OxygenCalculation 
    !----------------------------------------------------------------------------


    ! Calculation of the Hydrogen activity (mol/L) for specific rates calculations
    !----------------------------------------------------------------------------
    subroutine HydrogenCalculation(index)

        !Argument----------------------------------------------------------------
        integer                 , intent(IN)    :: index  

        !Local-------------------------------------------------------------------
        real    :: Hyd,Khn
        real    :: pH
        !------------------------------------------------------------------------
       
        pH    = Me%ExternalVar%pH(index)
    
        Hyd   =  1. * 10. ** (-pH )
                            
        Me%Hydrogen = Hyd

           
        
        if (pH <= 7) then
        
            khn = 0.167

            Me%Adj  = 1
            
        else
        
            khn = -0.333            
        
            Me%Adj = 3159.7

        end if


        Me%Khn = Khn
    
    
    end subroutine HydrogenCalculation
    !----------------------------------------------------------------------------
    

    ! Check if the Heterotrophs will have to immobilizate Mineral N
    !----------------------------------------------------------------------------
    subroutine LogicalImobilization (index)

        !Arguments---------------------------------------------------------------
        integer                 , intent(IN)    :: index  

        !Local-------------------------------------------------------------------
        integer :: NI, AMI              !N property index
                
        real    :: LC                   !Labil OM C content
        real    :: LCN                  !Labil CN ratio
        real    :: Pot_LOM_C_Decay      !Potential OM C Decay           [massC / Time]

        real    :: RC                   !Refractary OM C content
        real    :: RCN                  !Refractary CN ratio
        real    :: Pot_ROM_C_Decay      !Potential Refract OM C Decay   [massC / Time]
        
        real    :: N, AM, LN,RN         !inorganic, Labil or Refractary N content
        integer :: LCI,RNI,LNI,RCI          
        real    :: obtained             !Obtained N during OM decay    
        real    :: needed               !Needed N for the obtained C in OM decay
        real    :: PotAMImobil      !Potential Ammonia immobilizaion rate
        real    :: PotNIImobil      !Potential Nitrate immobilizaion rate
        !------------------------------------------------------------------------
        
        !calculate the Labil OM CN ratio
        LCI  = Me%PropIndex%Labil_OM_C
        LNI  = Me%PropIndex%Labil_OM_N           
        AMI  = Me%PropIndex%AMMONIA
        NI   = Me%PropIndex%Nitrate
        LC   = Me%ExternalVar%Mass(LCI, index)
        LN   = Me%ExternalVar%Mass(LNI, index)    !
        AM   = Me%ExternalVar%Mass(AMI, index)
        N    = Me%ExternalVar%Mass(NI, index)

        Me%LabiOM_CN_Ratio   = LC/LN 
        LCN                  = Me%LabiOM_CN_Ratio
        
        !calculate the Refractary OM CN ratio
        RCI  = Me%PropIndex%RefractOM_C
        RNI  = Me%PropIndex%RefractOM_N

        RC  = Me%ExternalVar%Mass(RCI, index)
        RN   = Me%ExternalVar%Mass(RNI, index)

        Me%RefractOM_CN_Ratio = RC/RN
        RCN                   = Me%RefractOM_CN_Ratio
   

     
        !Calculate the OM_C decay rates
        Pot_LOM_C_Decay     = LC * Me%SpecificRates%Labil_OM_C%Value
        Pot_ROM_C_Decay     = RC * Me%SpecificRates%RefractOM_C%Value


        !Potential Mineral N imobilization
        PotNIImobil         = N* Me%SpecificRates%NitrateImobilization%Value
        PotAMImobil         = AM * Me%SpecificRates%AmmoniaImobilization%Value


        !Obtained N during OM decay
        obtained            = Pot_LOM_C_Decay / LCN + Pot_ROM_C_Decay / RCN
        
        !Needed N for de C obtained in the OM decay              
        needed              =   (Pot_LOM_C_Decay + Pot_ROM_C_Decay) /                  & 
                                Me%Microorganisms%Heterotrophs%CNRatio  
        !Immobilization test
        if (obtained < needed ) then
            Me%Imobilization = ON    !the Heterotrophs will have to incorporate Mineral N
                            
        end if
        
        !Set Anaeribic Partition

        Me%AnaerobicPartition = (Pot_ROM_C_Decay)/(Pot_LOM_C_Decay )

        Me%Partition = (PotNIImobil) / (PotAMImobil )

        !------------------------------------------------------------------------

    end subroutine LogicalImobilization 
    !----------------------------------------------------------------------------



    ! Check If heterotrophs have to immobilize P
    !----------------------------------------------------------------------------
        
    subroutine LogicalImobilization_P(index)

        !Arguments---------------------------------------------------------------
        integer                 , intent(IN)    :: index  

        !Local-------------------------------------------------------------------
                
        real    :: LC                   !Labil OM C content
        real    :: LCP                  !Labil CP ratio
        real    :: Pot_LOM_C_Decay      !Potential OM C Decay           [massC / Time]

        real    :: RC                   !Refractary OM C content
        real    :: RCP                  !Refractary CP ratio
        real    :: Pot_ROM_C_Decay      !Potential Refract OM C Decay   [massC / Time]
        integer :: RPI,LPI,LCI,RCI       
        real    :: LP,RP    
          
        real    :: obtained             !Obtained P during OM decay    
        real    :: needed               !Needed P for the obtained C in OM decay
        !------------------------------------------------------------------------

        !calculate the Labil OM CP ratio

        LCI  = Me%PropIndex%Labil_OM_C
        LPI   = Me%PropIndex%Labil_OM_P

        LC   = Me%ExternalVar%Mass(LCI,index)
        LP   = Me%ExternalVar%Mass(LPI,index)

        Me%LabilOM_CP_Ratio = LC/LP
        LCP   = Me%LabilOM_CP_Ratio

        !calculate the Refractory OM CP ratio

        RCI  = Me%PropIndex%RefractOM_C
        RPI = Me%PropIndex%RefractOM_P

        RC  = Me%ExternalVar%Mass(RCI,index)
        RP  = Me%ExternalVar%Mass(RPI,index)

        Me%RefractOM_CP_Ratio = RC/RP
        RCP = Me%RefractOM_CP_Ratio 

        !Calculate the OM_C decay rates    assumindo o potencial!!!
        Pot_LOM_C_Decay     = LC * Me%SpecificRates%Labil_OM_C%Value
        Pot_ROM_C_Decay     = RC * Me%SpecificRates%RefractOM_C%Value


        !Obtained P during OM decay 

        obtained = (Pot_LOM_C_Decay/ LCP) + (Pot_ROM_C_Decay/RCP)
    
    
        !Needed P for de C obtained in the OM decay    
        
        needed   = (Pot_LOM_C_Decay+Pot_ROM_C_Decay)/Me%Microorganisms%Heterotrophs%CPRatio 
        
        if( obtained < needed) then

            Me%Imobilization_P = ON    !the Heterotrophs will have to incorporate Mineral P
                            
        end if
        
        !!! não sei se é importante ou não, mas vou por de novo a equaçao do partition      

!        Me%AnaerobicPartition = (Pot_ROM_C_Decay)/(Pot_LOM_C_Decay )



    end subroutine LogicalImobilization_P

    !------------------------------------------------------------------------


    ! Check if the Heterotrphs are limited by C or N
    !----------------------------------------------------------------------------
    subroutine LogicalLimitation (index)

        !Arguments---------------------------------------------------------------

        integer                 , intent(IN)    :: index  

        !Local-------------------------------------------------------------------
        
        integer :: RCI              !Refractary OM C index
        integer :: LCI              !Labil OM C index
        integer :: AMI              !Ammonia index
        integer :: NII              !Nitrate index
        real    :: RC               !Refractary OM C concentration
        real    :: RCN              !Refractary OM CN ratio
        real    :: Pot_ROM_C_Decay  !Potential Refractary OM C decay rate

        real    :: LC               !Labil OM C concentration
        real    :: LCN              !Labil OM CN ratio
        real    :: Pot_LOM_C_Decay  !Potential Labil OM C decay rate

        real    :: AM               !Ammonia concentration
        real    :: PotAMImobil      !Potential Ammonia immobilizaion rate

        real    :: NI               !Nitrate concentration
        real    :: PotNIImobil      !Potential Nitrate immobilizaion rate
                       
        real    :: MCN              !Heterotrophs OM CN ratio
                       
        real    :: ImobilizationS   !N immobilization rate
        real    :: CarbonS          !Carbon immobilization rate
        real    :: Conversion       !From dissolved to particulate
        !------------------------------------------------------------------------
        
        !get the OM values
        RCI = Me%PropIndex%RefractOM_C
        LCI = Me%PropIndex%Labil_OM_C

        RC  = Me%ExternalVar%Mass(RCI, index)
        LC  = Me%ExternalVar%Mass(LCI, index)
       

        !get de mineral forms of N values
        AMI = Me%PropIndex%Ammonia
        NII = Me%PropIndex%Nitrate

        AM  = Me%ExternalVar%Mass(AMI, index)
        NI  = Me%ExternalVar%Mass(NII, index)


        !get de mineral forms of N values
        LCN = Me%LabiOM_CN_Ratio
        RCN = Me%RefractOM_CN_Ratio
        MCN = Me%Microorganisms%Heterotrophs%CNRatio

                
        !Potential C uptake
        Pot_LOM_C_Decay     = LC * Me%SpecificRates%Labil_OM_C%Value
        Pot_ROM_C_Decay     = RC * Me%SpecificRates%RefractOM_C%Value


        !Potential Mineral N imobilization
        PotNIImobil         = NI * Me%SpecificRates%NitrateImobilization%Value
        PotAMImobil         = AM * Me%SpecificRates%AmmoniaImobilization%Value

        
        !Conversion from dissolved to Particulate
        Conversion      = Me%ExternalVar%DissolvedToParticulate (index) 

       
        !Rate that the Heterotrophs can uptake mineral N in mg kg-1 soil   
        ImobilizationS  = ( PotNIImobil + PotAMImobil ) * Conversion 


        !Rate that the Heterotrophs need extra Nitrogen (not in the OM) assuming the potential rates mg kg-1 soil
        CarbonS         = (Pot_LOM_C_Decay + Pot_ROM_C_Decay) / MCN -       & 
                           Pot_LOM_C_Decay / LCN - Pot_ROM_C_Decay / RCN    

        !Imobilization test        
        if (ImobilizationS < CarbonS ) then
            !The organic mather decay rates are function of de N immobilization rates
            Me%NLimitation = ON    
            !The immobilization rates will be a function of Potential OM carbon decay        
            !The OM carbon decay rates will be a function of Potential immobilization rates 
        
        end if
           

    end subroutine LogicalLimitation
    !----------------------------------------------------------------------------

  

     ! Check if the Heterotrphs are limited by C or P 
      !----------------------------------------------------------------------------
    subroutine LogicalLimitation_P (index)

        !Arguments---------------------------------------------------------------

        integer                 , intent(IN)    :: index  

        !Local-------------------------------------------------------------------
        
        integer :: RCI              !Refractary OM C index
        integer :: LCI              !Labil OM C index
        integer :: PI               !Phosphorus Soluble index
        real    :: RC               !Refractary OM C concentration
        real    :: RCP              !Refractary OM CP ratio
        real    :: Pot_ROM_C_Decay  !Potential Refractary OM C decay rate

        real    :: LC               !Labil OM C concentration
        real    :: LCP              !Labil OM CP ratio
        real    :: Pot_LOM_C_Decay  !Potential Labil OM C decay rate

        real    :: P                !Phosphorus soluble concentratio
        real    :: PotPImobil       !Potential Soluble Phosphorus immobilizaion rate

        real    :: MCP              !Heterotrophs OM CP ratio
                       
        real    :: ImobilizationS   !P immobilization rate
        real    :: CarbonS          !Carbon immobilization rate
        real    :: Conversion       !From dissolved to particulate
        !------------------------------------------------------------------------
 
         !get the OM values
 
        RCI = Me%PropIndex%RefractOM_C
        LCI = Me%PropIndex%Labil_OM_C
        PI  = Me%PropIndex%Inorganic_P_soluble
        
        !get the mineral form of P values
        
        RC  = Me%ExternalVar%Mass(RCI,index)
        LC  = Me%ExternalVar%Mass(LCI,index)
        P   = Me%ExternalVar%Mass(PI,index)

        !get the ratios values

        LCP = Me%LabilOM_CP_Ratio
        RCP = Me%RefractOM_CP_Ratio
        MCP = Me%Microorganisms%Heterotrophs%CPRatio
          
        !Potential C uptake

        Pot_LOM_C_Decay     = LC * Me%SpecificRates%Labil_OM_C%Value
        Pot_ROM_C_Decay     = RC * Me%SpecificRates%RefractOM_C%Value

        !Potential Phosphorus Imobilization
            
        PotPImobil          = P * Me%SpecificRates%PhosphorusImobilization%Value

        !Conversion from dissolved to Particulate
        Conversion      = Me%ExternalVar%DissolvedToParticulate (index) 


        !Rate that the Heterotrophs can uptake mineral P in ug kg-1 soil   
        ImobilizationS  = (PotPImobil  ) * Conversion 



        !Rate that the Heterotrophs need extra Phosphorus (not in the OM) assuming the potential rates ug kg-1 soil
        CarbonS         = (Pot_LOM_C_Decay + Pot_ROM_C_Decay) / MCP -       & 
                           Pot_LOM_C_Decay / LCP - Pot_ROM_C_Decay / RCP    



        !Imobilization test        
        if (ImobilizationS < CarbonS ) then
            !The organic mather decay rates are function of de P immobilization rates
            Me%PLimitation = ON    
              
        end if

    end subroutine LogicalLimitation_P 
  
    !----------------------------------------------------------------------------


    ! Check for special case when both immobilizations occur
    !----------------------------------------------------------------------------
    subroutine BothImobilization ()


        !Arguments---------------------------------------------------------------

!        integer                 , intent(IN)    :: index  

        !----------------------------------------------------------------------------

        if (Me%NLimitation .AND. Me%PLimitation) then
                
            Me%select  = 1              ! OM decay        : special case!
                                        ! N imobilization : potential
                                        ! P imobilization : potential   

        end if

        if (Me%NLimitation .AND. .NOT. Me%PLimitation ) then 
                
            Me%select = 2               ! OM decay        : Real N
                                        ! N imobilization : Potential
                                        ! P imobilization : real imobilization


        end if

        if (.NOT.Me%NLimitation .AND. Me%PLimitation )  then
                    

            Me%select = 4               ! OM decay        : real P
                                        ! N imobilization : real imobilization
                                        ! P imobilization : potential


        end if

    
        if (.NOT.Me%NLimitation .AND. .NOT.Me%PLimitation ) then
        
            Me%select = 5               ! OM decay        : potential
                                        ! N imobilization : real imobilization
                                        ! P imobilization : real imobilization

        end if

    end subroutine BothImobilization 

    !----------------------------------------------------------------------------



    ! Select the situation where organic matter decay and hete growth is defined
    !----------------------------------------------------------------------------
    subroutine SelectSituation (index)
    
        !Arguments---------------------------------------------------------------

        integer                 , intent(IN)    :: index  

         !------------------------------------------------------------------------


        if (Me%ComputeImobilization) then
        
            call LogicalImobilization (index)

            call LogicalImobilization_P (index)

            if (Me%Imobilization) then 
                
                call LogicalLimitation (index)
            end if

            if (Me%Imobilization_P) then

                call LogicalLimitation_P (index)
            end if
        
        endif
        
        !!! começamos aqui a selecionar os casos


        if (.NOT.Me%Imobilization .AND. .NOT.Me%Imobilization_P) then

            Me%select=9                 ! OM decay        : potential
                                        ! N imobilization : 0
                                        ! P imobilization : 0
        end if


        if (.NOT.Me%Imobilization .AND. Me%Imobilization_P) then
        
            if (Me%PLimitation) then 
        
                Me%select =7            ! OM decay        : Real P
                                        ! N imobilization : 0
                                        ! P imobilization : Potencial
            
            else

                Me%select = 8           ! OM decay        : potential
                                        ! N imobilization : 0
                                        ! P imobilization : Real Imobilization
        
            end if

        end if

        if (Me%Imobilization .AND. .NOT.Me%Imobilization_P) then

            if (Me%NLimitation) then
           
                Me%select =3            ! OM decay        : Real N
                                        ! N imobilization : Potential
                                        ! P imobilization : 0
            
            else

                Me%select =6            ! OM decay        : potential
                                        ! N imobilization : real imobilization
                                        ! P imobilization : 0
            end if
        end if
    
    
        if (Me%Imobilization .AND. Me%Imobilization_P) then 
    
            call BothImobilization ()

        end if 


    end subroutine SelectSituation 

    !----------------------------------------------------------------------------

    
    ! Calculates the Specific rates - day-1
    !----------------------------------------------------------------------------
    subroutine PotentialRatesCalc (index )

        !Arguments---------------------------------------------------------------
        integer, intent(IN)                 :: index

        !Local-------------------------------------------------------------------
        integer                             :: CLI
        integer                             :: CRI
        integer                             :: AMI,PFI,NI
        real(8)                             :: Tterm            ![day-1.pop-1] Temperature Rate
        real                                :: Aerobiose        ![0-1] Factor for considering aerobiose conditions
        real                                :: Anaerobiose      ![0-1] Factor for considering anaerobiose conditions
        real                                :: Oxygen           ![mol/L] Oxygen Concentration
        real                                :: Hydrogen         ![mol/L] Hydrogen Concentration
        real                                :: Temp             ![ºC] Temperature
        real                                :: Ionic            ![mol/L] Inonic Strenght
        real                                :: PF
        real                                :: CSUBST

        real                                :: Ammonia, Nitrate, Limit
        real                                :: Khn              !Exponent for Hydrogen ion dependent on pH (to balance units?)
        real                                :: Adj              !Coefficient dependent on pH (to balance units?)

        type(T_Coeficients)     , pointer   :: Rate
        type(T_Microorganisms)  , pointer   :: Micro
        !------------------------------------------------------------------------
        
        Aerobiose   = Me%Aerobiose
        Anaerobiose = Me%Anaerobiose

        Oxygen      = Me%Oxygen 
        Hydrogen    = Me%Hydrogen
        
        Temp        = Me%ExternalVar%Temperature(index)
        Ionic       = Me%ExternalVar%Ionic(index)
        Khn         = Me%Khn
        Adj         = Me%Adj
        Micro       => Me%Microorganisms        
        
        CLI         = Me%PropIndex%Labil_OM_C
        CRI         = Me%PropIndex%RefractOM_C
        CSUBST      = Me%ExternalVar%Mass(CLI, index) +         &
                      Me%ExternalVar%Mass(CRI, index) 
        
        PFI         = Me%PropIndex%Inorganic_P_fix
   
        PF          = Me%ExternalVar%Mass(PFI, index) 
        
        AMI         = Me%PropIndex%Ammonia
        NI          = Me%PropIndex%Nitrate
        Ammonia     = Me%ExternalVar%Mass(AMI, index)
        Nitrate     = Me%ExternalVar%Mass(NI, index) 
        


        ! Limit for specific rate : ammonia, nitrate and Csubs 
        Limit = 1. !1400000.000000                      


        call CalcPopulation (index)

        if (Me%PropCalc%Carbon) then
             
        !Calculates the OM C specific decay Rate                 
            
            Rate        => Me%SpecificRates%Labil_OM_C
           
            ![day-1.pop-1] where [pop] = [#org . kgsoil-1]
            Tterm       =   CalcTterm (Rate, Temp,Ionic)                                                 
            
            !units are only consistent if Khn and Adj balance O2 and H units terms
            ![day-1] = [-] * [day-1.pop-1] * [molO2/L] / [molH/L]^[-] * [pop] * [-]
            Rate%Value  =   Aerobiose * Tterm * Oxygen / Hydrogen**Khn       & 
                            * Micro%Heterotrophs%Population * Adj
            
        !Calculates the Refractary OM C specific decay Rate
            Rate        =>  Me%SpecificRates%RefractOM_C
            Tterm       =   CalcTterm (Rate, Temp,Ionic)                   
            Rate%Value  =   Aerobiose * Tterm * Oxygen / Hydrogen**Khn       &
                            * Micro%Heterotrophs%Population * Adj
        
        !Calculates the Heterotrophs C specific decay (death) Rate
        !ATTENTION! the added of CSUBST changes rate units!!
            if (.not. Me%ChangeRates) then
                Rate        =>  Me%SpecificRates%Heterotrophs
                Tterm       =   CalcTtermDeath (Rate, Temp,Ionic)
                Rate%Value  =   1./ Aerobiose * Tterm * Hydrogen**Khn /                   & 
                                ( Oxygen * Max(CSUBST,Limit)) * Micro%Heterotrophs%Population / Adj 
            else
                !in order to maintain units consistency CSUBST was removed
                Rate        =>  Me%SpecificRates%Heterotrophs
                Tterm       =   CalcTtermDeath (Rate, Temp,Ionic)
                Rate%Value  =   1./ Aerobiose * Tterm * Hydrogen**Khn /                   & 
                                 Oxygen * Micro%Heterotrophs%Population / Adj 
            endif
                        
        !Calculates the Autotrophs C specific decay (death) Rate
        !ATTENTION! the added of Ammonia changes rate units!!
            if (.not. Me%ChangeRates) then
                Rate        =>  Me%SpecificRates%Autotrophs
                Tterm       =   CalcTtermDeath (Rate, Temp,Ionic)
                Rate%Value  =   1./ Aerobiose * Tterm * Hydrogen**Khn /                   & 
                                ( Oxygen * Max(Ammonia,Limit)) * Micro%Autotrophs%Population / Adj 
            else
                !in order to maintain units consistency Ammonia was removed    
                Rate        =>  Me%SpecificRates%Autotrophs
                Tterm       =   CalcTtermDeath (Rate, Temp,Ionic)
                Rate%Value  =   1./ Aerobiose * Tterm * Hydrogen**Khn /                   & 
                                Oxygen * Micro%Autotrophs%Population / Adj 
            endif        
        
        !Calculates the Anaerobic C specific decay (death) Rate
        !ATTENTION! the added of CSUBST and nitrate changes rate units!!
            if (.not. Me%ChangeRates) then
                Rate        =>  Me%SpecificRates%Anaerobic
                Tterm       =   CalcTtermDeath (Rate, Temp,Ionic)
                Rate%Value  =   (1./ Anaerobiose * Tterm * Hydrogen**Khn /                 & 
                                ( Max(Nitrate,Limit)* Max(CSUBST,Limit) )) * Micro%anaerobic%Population / Adj  
            else
                !in order to maintain units consistency CSUBST and nitrate was removed
                !ATTENTION! without oxygen the units also are unbalanced so it was also removed hidrogen and Adj
                Rate        =>  Me%SpecificRates%Anaerobic
                Tterm       =   CalcTtermDeath (Rate, Temp,Ionic)
                Rate%Value  =   (1./ Anaerobiose * Tterm )                 & 
                                 * Micro%Anaerobic%Population   
            endif
            
        !Calculates the production of methane rate
            !ATTENTION! without oxygen the units also are unbalanced 
            if (.not. Me%ChangeRates) then
                Rate        =>  Me%SpecificRates%MethaneProduction
                Tterm       =   CalcTterm (Rate, Temp,Ionic)
                Rate%Value  =   (Anaerobiose * Tterm * Micro%Anaerobic%Population   &
                                /Hydrogen**Khn) * Adj
            else
                !in order to maintain units consistency Hidrogen and Adj were removed
                Rate        =>  Me%SpecificRates%MethaneProduction
                Tterm       =   CalcTterm (Rate, Temp,Ionic)
                Rate%Value  =   Anaerobiose * Tterm * Micro%Anaerobic%Population 
                                                                        
            endif
              
        end if      
                                                
        
        if (Me%PropCalc%Sol_Bacteria) then

        !Calculates the Solubilizing bacteria C specific decay (death) Rate
         !ATTENTION! the added of PF and lack of khn and sols population changes rate units!!
            if (.not. Me%ChangeRates) then
                Rate        =>  Me%SpecificRates%Sol
                Tterm       =   CalcTtermDeath (Rate, Temp,Ionic)
                Rate%Value  =   1./ Aerobiose * Tterm * Hydrogen /                   & 
                                ( Oxygen * (0.5+PF) ) /Adj 
            else
                Rate        =>  Me%SpecificRates%Sol
                Tterm       =   CalcTtermDeath (Rate, Temp,Ionic)
                Rate%Value  =   1./ Aerobiose * Tterm * Hydrogen**Khn /                   & 
                                ( Oxygen * Micro%Sols%Population ) /Adj 
            endif
    
        end if


        if (Me%PropCalc%Nitrogen) then
        
        !Calculates the AmmoniaToNitrate (nitrification) specific Rate  
            Rate        =>  Me%SpecificRates%AmmoniaToNitrate                   
            Tterm       =   CalcTterm (Rate, Temp,Ionic)
            Rate%Value  =   Aerobiose * Tterm * Oxygen**0.5 / Hydrogen**Khn        &
                            * Micro%Autotrophs%Population * Adj
                                                
        !Calculates the AmmoniaImobilization specific Rate
            Rate        =>  Me%SpecificRates%AmmoniaImobilization 
            Rate%Value  =   0.45
                                                        
        !Calculates the NitrateToNgas specific Rate
            !ATTENTION! the added of CSUBST changes rate units!!
            if (.not. Me%ChangeRates) then
                Rate        =>  Me%SpecificRates%NitrateToNgas       
                Tterm       =   CalcTterm (Rate, Temp,Ionic)
                Rate%Value  =   Anaerobiose * Tterm * CSUBST / Hydrogen **Khn             &
                                * Micro%Anaerobic%Population * Adj
           else 
                !in order to maintain units consistency CSUBST and nitrate was removed
                !ATTENTION! without oxygen the units also are unbalanced so it was also removed hidrogen and Adj
                Rate        =>  Me%SpecificRates%NitrateToNgas       
                Tterm       =   CalcTterm (Rate, Temp,Ionic)
                Rate%Value  =   Anaerobiose * Tterm * Micro%Anaerobic%Population
            endif
            
        !Calculates the NitrateImobilization specific Rate
            Rate        =>  Me%SpecificRates%NitrateImobilization 
            Rate%Value  =   0.45

        !Calculates the Urea hydrolysis rate
            !ATTENTION! day-1.pop-1. rate? Need population however this is done by an enzyme not 
            !a microorganism simulated?
            Rate        =>  Me%SpecificRates%UreaHydrolysis 
            Tterm       =   CalcTterm (Rate, Temp,Ionic)
            Rate%Value  =   Aerobiose * Tterm

                                                        
        end if  


        if (Me%PropCalc%Phosphorus) then                              !!!Lúcia

        !Calculates the Phosphorus Immobilizationspecific Rate  
                                                                        !!!Lúcia
            Rate        =>  Me%SpecificRates%PhosphorusImobilization    !!!Lúcia
            Tterm       =   CalcTterm (Rate, Temp,Ionic)
            Rate%Value  =   Aerobiose * Tterm * Oxygen / Hydrogen ** Khn               &
                            * Micro%Heterotrophs%Population * Adj

            if (Me%PropCalc%Sol_Bacteria) then
                    
                Rate        =>  Me%SpecificRates%Solubilizing                  
                Tterm       =   CalcTterm (Rate, Temp,Ionic)
                Rate%Value  =   Aerobiose * Tterm * Oxygen**0.5 / Hydrogen**Khn         &
                                * Micro%Sols%Population * Adj

             end if

          end if  


        !------------------------------------------------------------------------

    end subroutine PotentialRatesCalc
    !----------------------------------------------------------------------------


    ! Calculates the Specific rates - day-1
    !----------------------------------------------------------------------------
    subroutine PotentialRatesCalc_New (index )

        !Arguments---------------------------------------------------------------
        integer, intent(IN)                 :: index

        !Local-------------------------------------------------------------------
        integer                             :: CLI
        integer                             :: CRI
        integer                             :: AMI,PFI,NI
        real(8)                             :: MaximumRate      ![day-1] with Tterm [day-1.pop-1] Temperature Rate
        real                                :: AerobioseTerm    ![0-1] Factor for considering aerobiose conditions
        real                                :: AnaerobioseTerm  ![0-1] Factor for considering anaerobiose conditions
        real                                :: Oxygen           ![mol/L] Oxygen Concentration
        real                                :: Hydrogen         ![mol/L] Hydrogen Concentration
        real                                :: Temp             ![ºC] Temperature
        real                                :: Ionic            ![mol/L] Inonic Strenght
        real                                :: PF
        real                                :: CSUBST
        real                                :: Tterm, OxygenTerm, pHterm ![-]

        real                                :: Ammonia, Nitrate, pH
        real                                :: SubstrateTerm, SubstrateTerm1, SubstrateTerm2

        type(T_Coeficients)     , pointer   :: Rate
        type(T_Microorganisms)  , pointer   :: Micro
        !------------------------------------------------------------------------
        
        AerobioseTerm   = Me%Aerobiose
        AnaerobioseTerm = Me%Anaerobiose

        Oxygen      = Me%Oxygen 
        Hydrogen    = Me%Hydrogen
        pH          = Me%ExternalVar%pH(index)
        
        Temp        = Me%ExternalVar%Temperature(index)
        Ionic       = Me%ExternalVar%Ionic(index)
        Micro       => Me%Microorganisms        
        
        CLI         = Me%PropIndex%Labil_OM_C
        CRI         = Me%PropIndex%RefractOM_C
        CSUBST      = Me%ExternalVar%Mass(CLI, index) +         &
                      Me%ExternalVar%Mass(CRI, index) 
        
        PFI         = Me%PropIndex%Inorganic_P_fix
   
        PF          = Me%ExternalVar%Mass(PFI, index) 
        
        AMI         = Me%PropIndex%Ammonia
        NI          = Me%PropIndex%Nitrate
        Ammonia     = Me%ExternalVar%Mass(AMI, index)
        Nitrate     = Me%ExternalVar%Mass(NI, index) 
        

        call CalcPopulation (index)


        if (Me%PropCalc%Carbon) then
             
        !Calculates the OM C specific decay Rate                 
            
            Rate        => Me%SpecificRates%Labil_OM_C
           
            !maximum rate [day-1] = [day-1.pop-1] . [pop],  where [pop] = [#org . kgsoil-1]
            Tterm       =   CalcTterm (Rate, Temp,Ionic)
            MaximumRate =   Tterm * Micro%Heterotrophs%Population 
            
            ![-] factors between 0 and 1 affecting maximum rate
            OxygenTerm  =   CalcOxygenTerm (Rate, Oxygen)
            pHTerm      =   CalcpHTerm (Rate, pH)
            
            ![day-1] = [day-1] * [-] * [-] * [-]
            Rate%Value  =   MaximumRate * AerobioseTerm * OxygenTerm * pHTerm 
            
        !Calculates the Refractary OM C specific decay Rate
            Rate        =>  Me%SpecificRates%RefractOM_C
            Tterm       =   CalcTterm (Rate, Temp,Ionic)
            MaximumRate =   Tterm * Micro%Heterotrophs%Population               
            OxygenTerm  =   CalcOxygenTerm (Rate, Oxygen)
            pHTerm      =   CalcpHTerm (Rate, pH)          
            Rate%Value  =   MaximumRate * AerobioseTerm * OxygenTerm * pHTerm       
        
        !CSUBST
        !Calculates the Heterotrophs C specific decay (death) Rate
            Rate        =>  Me%SpecificRates%Heterotrophs
            Tterm       =   CalcTtermDeath (Rate, Temp,Ionic)
            MaximumRate =   Tterm * Micro%Heterotrophs%Population
            OxygenTerm  =   CalcOxygenTerm (Rate, Oxygen)
            pHTerm      =   CalcpHTerm (Rate, pH)
            if (CSUBST .gt. 0.) then
                SubstrateTerm = min(Rate%ConcMinCarbon/CSUBST, 1.) 
            else !maximum death rate
                SubstrateTerm = 1.
            endif
            Rate%Value  =   MaximumRate * (1./ AerobioseTerm) * (1. - pHTerm) * (1. -  OxygenTerm) * SubstrateTerm
        
        !Ammonia                
        !Calculates the Autotrophs C specific decay (death) Rate
            Rate        =>  Me%SpecificRates%Autotrophs
            Tterm       =   CalcTtermDeath (Rate, Temp,Ionic)
            MaximumRate =   Tterm * Micro%Autotrophs%Population
            OxygenTerm  =   CalcOxygenTerm (Rate, Oxygen)
            pHTerm      =   CalcpHTerm (Rate, pH)
            if (Ammonia .gt. 0.) then
                SubstrateTerm = min(Rate%ConcMinAmmonia/Ammonia, 1.)
            else
                SubstrateTerm = 1.
            endif
            Rate%Value  =   MaximumRate * (1./ AerobioseTerm) * (1. - pHTerm) * (1. -  OxygenTerm) * SubstrateTerm 
                            
        !Nitrate CSUBST
        !Calculates the Anaerobic C specific decay (death) Rate
            Rate        =>  Me%SpecificRates%Anaerobic
            Tterm       =   CalcTtermDeath (Rate, Temp,Ionic)
            MaximumRate =   Tterm * Micro%Anaerobic%Population
            pHTerm      =   CalcpHTerm (Rate, pH)            
            if (CSUBST .gt. 0.) then
                SubstrateTerm1 = min(Rate%ConcMinCarbon/CSUBST, 1.)
            else
                SubstrateTerm1 = 1.
            endif
            if (Nitrate .gt. 0.) then
                SubstrateTerm2 = min(Rate%ConcMinNitrate/Nitrate, 1.)
            else
                SubstrateTerm2 = 1.
            endif
            Rate%Value  =   MaximumRate * (1./ AnaerobioseTerm)  * (1. - pHTerm) * SubstrateTerm1 * SubstrateTerm2
                            
            
        !Calculates the production of methane rate
            Rate        =>  Me%SpecificRates%MethaneProduction
            Tterm       =   CalcTterm (Rate, Temp,Ionic)
            MaximumRate =   Tterm * Micro%Anaerobic%Population 
            pHTerm      =   CalcpHTerm (Rate, pH)                  
            Rate%Value  =   MaximumRate * AnaerobioseTerm  * pHTerm

              
        end if      
                                                
        
        if (Me%PropCalc%Sol_Bacteria) then

        !Calculates the Solubilizing bacteria C specific decay (death) Rate
            Rate        =>  Me%SpecificRates%Sol
            Tterm       =   CalcTtermDeath (Rate, Temp,Ionic)
            MaximumRate =   Tterm * Micro%Sols%Population
            OxygenTerm  =   CalcOxygenTerm (Rate, Oxygen)
            pHTerm      =   CalcpHTerm (Rate, pH)
            if (PF .gt. 0.) then
                SubstrateTerm = min(Rate%ConcMinPF/PF, 1.)
            else
                SubstrateTerm = 1.
            endif            
            Rate%Value  =   MaximumRate * (1./ AerobioseTerm) * (1. - pHTerm) * (1. -  OxygenTerm) * SubstrateTerm    
    
        end if


        if (Me%PropCalc%Nitrogen) then
        
        !Calculates the AmmoniaToNitrate (nitrification) specific Rate  
            Rate        =>  Me%SpecificRates%AmmoniaToNitrate   
            Tterm       =   CalcTterm (Rate, Temp,Ionic)                
            MaximumRate =   Tterm * Micro%Autotrophs%Population 
            OxygenTerm  =   CalcOxygenTerm (Rate, Oxygen)
            pHTerm      =   CalcpHTerm (Rate, pH)                              
            Rate%Value  =   MaximumRate * AerobioseTerm  * OxygenTerm * pHTerm  
                            
                                                
        !Calculates the AmmoniaImobilization specific Rate
            Rate        =>  Me%SpecificRates%AmmoniaImobilization 
            Rate%Value  =   0.45
        
        !Substrate CSUBST                                                        
        !Calculates the NitrateToNgas specific Rate
            Rate        =>  Me%SpecificRates%NitrateToNgas      
            pHTerm      =   CalcpHTerm (Rate, pH)            
            Tterm       =   CalcTterm (Rate, Temp,Ionic)                               
            MaximumRate =   Tterm * Micro%Anaerobic%Population
            !with increasing substrate, substrate term is closer to 1.
            if (CSUBST .gt. 0.) then
                SubstrateTerm = min(CSUBST/Rate%ConcOptCarbon, 1.)
            else
                SubstrateTerm = 0.
            endif            
            Rate%Value  =   MaximumRate * AnaerobioseTerm * pHTerm * SubstrateTerm
                            
        !Calculates the NitrateImobilization specific Rate
            Rate        =>  Me%SpecificRates%NitrateImobilization 
            Rate%Value  =   0.45

        !Calculates the Urea hydrolysis rate
            !ATTENTION! day-1.pop-1. rate? Need population however this is done by an enzyme not 
            !a microorganism simulated?
            Rate        =>  Me%SpecificRates%UreaHydrolysis 
            Tterm       =   CalcTterm (Rate, Temp,Ionic)
            MaximumRate =   Tterm
            Rate%Value  =   AerobioseTerm * MaximumRate

                                                        
        end if  


        if (Me%PropCalc%Phosphorus) then                              

        !Calculates the Phosphorus Immobilizationspecific Rate  
                                                                        
            Rate        =>  Me%SpecificRates%PhosphorusImobilization
            Tterm       =   CalcTterm (Rate, Temp,Ionic)    
            MaximumRate =   Tterm * Micro%Heterotrophs%Population
            OxygenTerm  =   CalcOxygenTerm (Rate, Oxygen)
            pHTerm      =   CalcpHTerm (Rate, pH)                                                       
            Rate%Value  =   MaximumRate * AerobioseTerm * OxygenTerm * pHTerm
                            

            if (Me%PropCalc%Sol_Bacteria) then
                    
                Rate        =>  Me%SpecificRates%Solubilizing  
                Tterm       =   CalcTterm (Rate, Temp,Ionic)                
                MaximumRate =   Tterm * Micro%Sols%Population
                OxygenTerm  =   CalcOxygenTerm (Rate, Oxygen)
                pHTerm      =   CalcpHTerm (Rate, pH)                                                       
                Rate%Value  =   MaximumRate * AerobioseTerm * OxygenTerm * pHTerm 
                                

             end if

          end if  


        !------------------------------------------------------------------------

    end subroutine PotentialRatesCalc_New
    !----------------------------------------------------------------------------


    ! Calculates the temperature term of the Specific rates - [day-1 . pop-1] where [pop] = [#org . kgsoil-1]
    !----------------------------------------------------------------------------
    function   CalcTterm (Coeficient, Temperature,Ionic)
    real(8)     :: CalcTterm

        !Arguments---------------------------------------------------------------
        type(T_Coeficients) , pointer                       :: Coeficient
        real                                                :: Temperature
        real                                                :: Ionic        
        !Local-------------------------------------------------------------------
        real    :: OptimumTemp                               ! [ºC] Optimum Temperature 
        real    :: tcalc                                     ! [ºK] Rate Temperature
        real    :: AE                                        ! [kcal.mole-1] Initial Activation Energy -> constant
        real    :: A                                         ! [s.day-1.pop-1] Rate Coefficient
        real    :: kp                                        ! [L/mol] Salinity Coefficient
        real    :: I                                         ! [mol/L] Ionic Strenght
        real    :: Eactivation                               ! [kcal.mole-1] Activation Energy
        !------------------------------------------------------------------------
        
            OptimumTemp     = Coeficient%OptimumTemperature  !rate optimum temperature
            tcalc           = Temperature                    !first equal to soil temperature
            if (tcalc > OptimumTemp) tcalc = 2 * OptimumTemp - tcalc !now updated

            AE          = Coeficient%ActivationE        
            kp          = Coeficient%kp                 
            A           = Coeficient%Acoef              
            tcalc       = tcalc + 273.15                
            I           = Ionic

            Eactivation = AE + kp * I           ! True activation energy based on Ionic strenght
 
            
            CalcTterm = (Boltzman * tcalc / Planck) * A * &
                        exp (- Eactivation / (UnivGC * tcalc) ) 
        !------------------------------------------------------------------------

    end function CalcTterm
    !----------------------------------------------------------------------------

    ! Calculates the temperature term of the Specific rates - [day-1 . pop-1] where [pop] = [#org . kgsoil-1]
    !----------------------------------------------------------------------------
    function   CalcTtermDeath (Coeficient, Temperature,Ionic)
    real(8)     :: CalcTtermDeath

        !Arguments---------------------------------------------------------------
        type(T_Coeficients) , pointer                       :: Coeficient
        real                                                :: Temperature
        real                                                :: Ionic
        
        !Local-------------------------------------------------------------------
        real    :: OptimumTemp
        real    :: tcalc
        real    :: AE
        real    :: A
        real    :: kp
        real    :: I
        real    :: Eactivation
        !------------------------------------------------------------------------
        
            OptimumTemp        = Coeficient%OptimumTemperature 
            tcalc              = Temperature
            if (tcalc > OptimumTemp) tcalc = 2 * OptimumTemp - tcalc

            AE          = Coeficient%ActivationE
            kp          = Coeficient%kp
            A           = Coeficient%Acoef
            tcalc       = tcalc + 273.15
            I           = Ionic

            Eactivation = AE + kp * I
 
 
            CalcTtermDeath = (Boltzman * tcalc / Planck) * A / &
                        exp (- Eactivation / (UnivGC * tcalc) ) 
        !------------------------------------------------------------------------

    end function CalcTtermDeath
    !----------------------------------------------------------------------------


    ! Calculates the oxygen term of the Specific rates - [-]
    !if conc O2 lower than optimum than value lower than 1. if higher, value is 1
    !----------------------------------------------------------------------------
    function   CalcOxygenTerm (Coeficient, Oxygen)
    real(8)     :: CalcOxygenTerm

        !Arguments---------------------------------------------------------------
        type(T_Coeficients) , pointer                       :: Coeficient
        real                                                :: Oxygen, ConcOpt
        !------------------------------------------------------------------------
        
        !mol/L = mg/L * 1E-3 g/mg / 32g/mol
        ConcOpt = Coeficient%ConcOptO2 * 1E-3 / 32.
        
        CalcOxygenTerm = min (Oxygen / ConcOpt, 1.)
        !------------------------------------------------------------------------

    end function CalcOxygenTerm
    !----------------------------------------------------------------------------

    ! Calculates the pH term of the Specific rates - [-]
    !if pH lower or higher than optimum than value lower than 1
    !----------------------------------------------------------------------------
    function   CalcpHTerm (Coeficient, pH)
    real(8)     :: CalcpHTerm

        !Arguments---------------------------------------------------------------
        type(T_Coeficients) , pointer                       :: Coeficient
        real                                                :: pH, optimumpH
        !------------------------------------------------------------------------
        real                                                :: pHread
        
        pHread = pH
        OptimumpH = Coeficient%OptimumpH
        
        if (pHread > OptimumpH) pHread = 2*OptimumpH - pHread  !adaptation from temperature to make triangle with optimum
        
        !from here all the pH were transformed in pH lower than optimum (mirror around optimum)
        !optimum gives 1 and lower give lower accordingly
        CalcpHTerm = pHread / (2 * OptimumpH - pHread)
       
       
        !------------------------------------------------------------------------

    end function CalcpHTerm
    !----------------------------------------------------------------------------

    ! Calculates the temperature term of the Specific rates
    !----------------------------------------------------------------------------
!    function   CalcActivationEnergy (coeficient,I,Ea)
!    real(8)     :: CalcActivationEnergy

        !Arguments---------------------------------------------------------------
!        type(T_Coeficients) , pointer                       :: Coeficient
!       real                                                :: I        
!       real                                                :: Ea
        !Local-------------------------------------------------------------------

!       real    ::  E0i
!       real    ::  kp
        !----------------------------------------------------------------------------




!       E0i = Coeficient%ActivationE
!       kp  = Coeficient%kp
!       I   = Me%ExternalVar%Ionic


!       Ea = E0i +kp*I


    !----------------------------------------------------------------------------
!       end function CalcActivationEnergy 
    !----------------------------------------------------------------------------



    ! Calculates the microorganisms population - #org.kgsoil-1
    !----------------------------------------------------------------------------
    subroutine CalcPopulation (index)

        !Arguments---------------------------------------------------------------
        integer                 , intent(IN)                    :: index

        !Local-------------------------------------------------------------------
        real                           :: Carbon
        real                           :: CpopRatio
        integer                        :: propindex
        integer                        :: poppropindex  
        
        type(T_Constants)  , pointer   :: Microorganisms 
        !------------------------------------------------------------------------
        
        Microorganisms    => Me%Microorganisms%Heterotrophs

            propindex   = Microorganisms%MicroIndex
            Carbon      = Me%ExternalVar%Mass(propindex, index) 
            CpopRatio   = Microorganisms%CPopRatio

            Microorganisms%Population = Carbon * CpopRatio
            
            !To go to external module
            poppropindex                             = Me%PropIndex%HeterotrophicPop
            Me%ExternalVar%Mass(poppropindex, index) = Microorganisms%Population
            Me%IndTerm(poppropindex)                 = Me%ExternalVar%Mass(poppropindex, index)
            Me%Matrix(poppropindex, poppropindex)    = 1.0
            
            if (Microorganisms%Population .LE.  Microorganisms%MinimumPop   ) then
                Microorganisms%LogicalMinumumPOP = ON
            else
                Microorganisms%LogicalMinumumPOP = OFF
            end if             
        
        
        Microorganisms    => Me%Microorganisms%Autotrophs        

            propindex   = Microorganisms%MicroIndex
            Carbon      = Me%ExternalVar%Mass(propindex, index) 
            CpopRatio   = Microorganisms%CPopRatio

            Microorganisms%Population = Carbon * CpopRatio

            !To go to external module
            poppropindex                             = Me%PropIndex%AutotrophicPop
            Me%ExternalVar%Mass(poppropindex, index) = Microorganisms%Population
            Me%IndTerm(poppropindex)                 = Me%ExternalVar%Mass(poppropindex, index)
            Me%Matrix(poppropindex, poppropindex)    = 1.0
            
            if (Microorganisms%Population .LE.  Microorganisms%MinimumPop   ) then
                Microorganisms%LogicalMinumumPOP = ON 
            else
                Microorganisms%LogicalMinumumPOP = OFF
            end if                       


        Microorganisms    => Me%Microorganisms%Anaerobic

            propindex   = Microorganisms%MicroIndex
            Carbon      = Me%ExternalVar%Mass(propindex, index) 
            CpopRatio   = Microorganisms%CPopRatio

            Microorganisms%Population = Carbon * CpopRatio

            !To go to external module
            poppropindex                             = Me%PropIndex%AnaerobicPop
            Me%ExternalVar%Mass(poppropindex, index) = Microorganisms%Population
            Me%IndTerm(poppropindex)                 = Me%ExternalVar%Mass(poppropindex, index)
            Me%Matrix(poppropindex, poppropindex)    = 1.0
            
            if (Microorganisms%Population .LE.  Microorganisms%MinimumPop   ) then
                Microorganisms%LogicalMinumumPOP = ON
            else
                Microorganisms%LogicalMinumumPOP = OFF
            end if             
 



        if ( Me%PropCalc%Sol_Bacteria) then

            Microorganisms    => Me%Microorganisms%Sols

                propindex   = Microorganisms%MicroIndex
                Carbon      = Me%ExternalVar%Mass(propindex, index) 
                CpopRatio   = Microorganisms%CPopRatio

                Microorganisms%Population = Carbon * CpopRatio

                !To go to external module
                poppropindex                             = Me%PropIndex%SolPop
                Me%ExternalVar%Mass(poppropindex, index) = Microorganisms%Population
                Me%IndTerm(poppropindex)                 = Me%ExternalVar%Mass(poppropindex, index)
                Me%Matrix(poppropindex, poppropindex)    = 1.0
                
                if (Microorganisms%Population .LE.  Microorganisms%MinimumPop   ) then
                    Microorganisms%LogicalMinumumPOP = ON
                else
                    Microorganisms%LogicalMinumumPOP = OFF
                end if             

        end if
        
        nullify(Microorganisms        )                     
        !------------------------------------------------------------------------

    end subroutine CalcPopulation
    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------
    subroutine StartSedimentQualityIteration


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

    end subroutine StartSedimentQualityIteration
    !----------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine SedimentQualityCoefCalculation(index)        !!!Lúcia

        !Arguments-------------------------------------------------------------
        integer, intent(IN)         :: index
        !----------------------------------------------------------------------

                                      call SQOxygen       (index)
        if (Me%PropCalc%Carbon    )   call SQCarbon       (index)
        if (Me%PropCalc%Nitrogen  )   call SQNitrogen     (index)
        if (Me%PropCalc%Phosphorus)   call SQPhosphorus   (index)         !!!Lúcia
        if (Me%PropCalc%Sol_Bacteria) call SQSol_Bacteria (index)         !!!Lúcia
        !----------------------------------------------------------------------

    end subroutine 
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine SystemResolution(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)              :: index

        !Extertnal-------------------------------------------------------------
        integer                     :: STAT_CALL
        real, pointer, dimension(:) ::  x

        !Local-----------------------------------------------------------------
        integer :: PropLB, PropUB
        integer :: prop
        integer :: equa
        !----------------------------------------------------------------------

        propLB = Me%Prop%ILB 
        propUB = Me%Prop%IUB 

        !Resolution using an explicit method
cd1 :   if (Me%CalcMethod%ExplicitMethod) then
do1 :       do equa = PropLB, PropUB           !Percorre as equacoes
do2 :       do prop = PropLB, PropUB           !Percorre as propriedades

cd2 :           if (Me%Matrix(equa, prop) .NE. 0.)  then

cd3 :               if (equa .EQ. prop) then
                    
                        Me%IndTerm(equa) =  Me%IndTerm (equa )                          &
                                         + ((Me%Matrix (equa, prop ) - 1. )             &
                                         * Me%ExternalVar%Mass (prop, index) )
                    else cd3
                    
                        Me%IndTerm(equa) =  Me%IndTerm (equa)                           &
                                         +  Me%Matrix    (equa, prop)                   &
                                         *  Me%ExternalVar%Mass(prop, index)
                    end if cd3

                    Me%NewMass(equa) = Me%IndTerm(equa)
            
                end if cd2
            
            end do do2
            end do do1


do4 :       do equa = PropLB, PropUB           !Percorre as equacoes
                Me%ExternalVar%Mass(equa, index) = Me%NewMass(equa)
            end do do4

        else if (Me%CalcMethod%ImplicitMethod) then

            !Resolution using an implicit method
            nullify   (x)
            
            call LUD(Me%ObjLUD, Me%Matrix, Me%IndTerm, x, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine SystemResolution; module ModuleSedimentQuality. ERR03.'

do3 :       do prop = PropLB, PropUB
                Me%ExternalVar%Mass(prop, index) = x(prop)
            end do do3

            nullify   (x)
        
        else if (Me%CalcMethod%SemiimpMethod) then

do31 :       do equa = PropLB, PropUB           !Percorre as equacoes
do32 :       do prop = PropLB, PropUB           !Percorre as propriedades

cd32 :           if (Me%Matrix(equa, prop) .GT. 0.)  then

cd33 :               if (equa .EQ. prop) then
                        Me%IndTerm(equa) = Me%IndTerm (equa)                        &
                                         + ((Me%Matrix (equa, prop) - 1.0)          &
                                         * Me%ExternalVar%Mass(prop, index))
                                 
                        Me%Matrix(equa, prop) = 1.0
                    else cd33
                        Me%IndTerm(equa) = Me%IndTerm (equa)                        &
                                         + Me%Matrix (equa, prop)                   &
                                         * Me%ExternalVar%Mass(prop, index)

                        Me%Matrix(equa, prop) = 0.0
                    end if cd33

                end if cd32

            end do do32
            end do do31
       

            !Resolution using an implicit method
            nullify   (x)
            
            call LUD(Me%ObjLUD, Me%Matrix, Me%IndTerm, x, STAT = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine SystemResolution; module ModuleSedimentQuality. ERR04.'


do33 :      do prop = PropLB, PropUB
                Me%ExternalVar%Mass(prop, index) = x(prop)
            end do do33
            
       
            nullify   (x)

        end if cd1
        !----------------------------------------------------------------------

    end subroutine SystemResolution
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine SedimentQualityRatesCalculation(index)
        
        !Arguments-------------------------------------------------------------
        integer, intent(IN) :: index

        type(T_PropRateFlux         ), pointer :: PropRateFluxX
        type(T_EquaRateFlux         ), pointer :: EquaRateFluxX 

        !Local-----------------------------------------------------------------
        integer :: PropUB,PropLB,equa,prop
        !----------------------------------------------------------------------

        propLB = Me%Prop%ILB 
        propUB = Me%Prop%IUB

          
        EquaRateFluxX => Me%FirstEquaRateFlux

do1 :   do while(associated(EquaRateFluxX))  
            
            PropRateFluxX => EquaRateFluxX%FirstPropRateFlux
            
            do while(associated(PropRateFluxX))  
                 
                equa = EquaRateFluxX%ID
                prop = PropRateFluxX%ID
                                                            

                if (equa.ne.prop) then
            
                    PropRateFluxX%Field(index)  = Me%Matrix(equa, prop)         & 
                                                * Me%ExternalVar%Mass(prop,index)
                else
                    !this can not be as in waterquality - to matrix has to be removed 1.
                    !just as system resolution
                    PropRateFluxX%Field(index)  = (Me%Matrix(equa, prop) - 1.)   & 
                                                * Me%ExternalVar%Mass(prop,index)                       
                endif

                PropRateFluxX => PropRateFluxX%Next
            
            end do
         
            EquaRateFluxX => EquaRateFluxX%Next

         end do do1

          !unidades
          
          !RateFlux = Matrix * Mass
          !mg/l        1.      mg/l 

          !Matrix = dt * K
          !   1.   = DtDay * 1./Days          
        
        nullify(PropRateFluxX)
        nullify(EquaRateFluxX)
        ! ---------------------------------------------------------------------------

    end subroutine SedimentQualityRatesCalculation
    !----------------------------------------------------------------------------

    !OXYGEN - for now it is constant and sources and sinks are not accounted
    ! mol/L
    !SOURCES: - 
    !SINKS:   - 
    !----------------------------------------------------------------------------
    subroutine SQOxygen (index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN)         :: index
        
        
        !Local-------------------------------------------------------------------
        integer                     :: O
        !------------------------------------------------------------------------

        O   = Me%PropIndex%Oxygen

        !Calculation of system coeficients---------------------------------------
        Me%Matrix(O, O) = 1.0 

        !Independent term
        Me%IndTerm(O) = Me%ExternalVar%Mass(O, index)


    end subroutine SQOxygen
    !----------------------------------------------------------------------------
    

    !----------------------------------------------------------------------------
    subroutine SQCarbon(index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN) :: index
        !------------------------------------------------------------------------
                                                        
        call SQLabilOrganicCarbon     (index)

        call SQRefractOrganicCarbon   (index)

        call SQHeterotrophicC       (index)

        call SQAutotrophicC         (index)

        call SQAnaerobicC           (index)
        
        call SQCO2                  (index)     

!       call SQMethane              (index)       
    !------------------------------------------------------------------------

    end subroutine SQCarbon
    !---------------------------------------------------------------------------


    !Labil Organic Carbon
    !
    !SOURCES: - Microorganisms death 
    !SINKS:   - Heterotrophs Aerobic and Anaerobic uptake
    !----------------------------------------------------------------------------    
    subroutine SQLabilOrganicCarbon (index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN)         :: index

        !Local-------------------------------------------------------------------
        integer :: AMI, NII,LOM_CI,ROM_CI,PI            !!!Lúcia
        
        real    :: ImobilizationRateNH4
        real    :: DenitrificationRate
        real    :: ImobilizationRateNO3
        real    :: ImobilizationRateP
        real    :: solubilizingRate
        real    :: LOM_C_Srate

        real    :: LCN
        real    :: LCP                          !!!Lúcia
        real    :: RCN
        real    :: RCP                          !!!Lúcia

        real    :: HeteroDeathRate, HetCN , HetCP    !!!Lúcia
        integer :: HetI                
        
        real    :: AnaDeathRate, AnaCN, AnaCP     !!!Lúcia   
        integer :: AnaI                  

        real    :: AutoDeathRate, AutoCN , AutoCP !!!Lúcia
        integer :: AutoI               

        real    :: Partition, AnaPartition , seletion,solpartition   !!!Lúcia
        real    :: Conversion, DTDay
        integer :: PFI

        real    :: AM,N,P,Real_N,Real_P         !!!Lúcia
        real    :: NO3,NO3limit,MethaneProduction
        !------------------------------------------------------------------------
        LOM_C_Srate         = Me%SpecificRates%Labil_OM_C%Value
        LCN                 = Me%LabiOM_CN_Ratio
        LOM_CI              = Me%PropIndex%Labil_OM_C
        ROM_CI              = Me%PropIndex%RefractOM_C
        LCP                 = Me%LabilOM_CP_Ratio       !!!Lúcia  
        
        PFI                 = Me%PropIndex%Inorganic_P_fix  
    
            
        RCN                 = Me%RefractOM_CN_Ratio   
        RCP                 = Me%RefractOM_CP_Ratio         !!!Lúcia
                
        DenitrificationRate     = Me%SpecificRates%NitrateToNgas%Value
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value
        ImobilizationRateP      = Me%SpecificRates%PhosphorusImobilization%Value  !!!Lúcia
        solubilizingrate        = Me%SpecificRates%Solubilizing%Value           
            
        NII                 = Me%PropIndex%Nitrate
        AMI                 = Me%PropIndex%Ammonia
        PI                  = Me%PropIndex%Inorganic_P_soluble   !!!Lúcia

        Partition           = Me%Partition
        AnaPartition        = Me%AnaerobicPartition
        solpartition        = Me%Solpartition
        seletion            = Me%Select                          !!!Lúcia

        HeteroDeathRate     = Me%SpecificRates%Heterotrophs%Value
        HetI                = Me%PropIndex%HeterotrophicC
        HetCN               = Me%Microorganisms%Heterotrophs%CNRatio
        HetCP               = Me%Microorganisms%Heterotrophs%CPRatio        !!!Lúcia

        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaI                = Me%PropIndex%AnaerobicC
        AnaCN               = Me%Microorganisms%Anaerobic%CNRatio
        AnaCP               = Me%Microorganisms%Anaerobic%CPRatio

        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutoI               = Me%PropIndex%AutotrophicC
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio
        AutoCP              = Me%Microorganisms%Autotrophs%CPRatio

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)

        DTDay               = Me%DTDay

        AM   = Me%ExternalVar%Mass(AMI,index)   !!!
        N    = Me%ExternalVar%Mass(NII,index)   !!!
        P    = Me%ExternalVar%Mass(PI,index)    !!!

        NO3limit            = Me%NO3limit
        methaneproduction   = Me%SpecificRates%MethaneProduction%Value
        NO3                 = Me%ExternalVar%Mass(NII,index)



        Me%Matrix(LOM_CI, LOM_CI)  = 1.


        !Sinks  : Heterotrophs uptake:  will depende on the type of biomass growth
        
            
        if (seletion==5.OR.seletion==6 .OR. seletion==8 .OR. seletion==9 )  then 

            ! biomass potential growth
        
            Me%Matrix(LOM_CI, LOM_CI)  = Me%Matrix(LOM_CI, LOM_CI) - DTDay * LOM_C_Srate  

        end if

        if (seletion==2.OR. seletion==3) then 

            ! biomass real growth limited by nitrogen 

            Me%Matrix(LOM_CI, AMI)  = Me%Matrix(LOM_CI, AMI) - DTDay                & 
                                    * ImobilizationRateNH4  * Conversion            &
                                    * (1./( (1./HetCN - 1./LCN ) + AnaPartition     &
                                    * (1./HetCN - 1./RCN))) 
                
            Me%Matrix(LOM_CI, NII)  = Me%Matrix(LOM_CI, NII) - DTDay                &
                                    * ImobilizationRateNO3  * Conversion            & 
                                    * (1./( (1./HetCN - 1./LCN ) + AnaPartition     &
                                    * (1./HetCN - 1./RCN)))
         end if

        if (seletion==4.OR. seletion==7) then

            ! biomass real growth limited by phosphorus
                    
            Me%Matrix(LOM_CI, PI)  = Me%Matrix(LOM_CI, PI) - DTDay                  &
                                    * ImobilizationRateP  * Conversion              & 
                                    * (1./( (1./HetCP - 1./LCP ) + AnaPartition     &
                                    * (1./HetCP - 1./RCP)))

        end if

        
        if(seletion==1) then
    
            ! special case

            Real_N = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)               &
                      *Conversion)*(1./( (1./HetCN - 1./LCN )                        &
                      + AnaPartition * (1./HetCN - 1./RCN)))

            Real_P = ((ImobilizationRateP*P)*Conversion)                             &
                    * (1./( (1./HetCP - 1./LCP ) +                                   &
                    AnaPartition * (1./HetCP - 1./RCP)))

    
            if( Real_N<Real_P) then

                  Me%Matrix(LOM_CI, AMI)  = Me%Matrix(LOM_CI, AMI) - DTDay           & 
                                    * ImobilizationRateNH4  * Conversion             &
                                    * (1./( (1./HetCN - 1./LCN ) +                   &
                                    AnaPartition * (1./HetCN - 1./RCN))) 
                
                  Me%Matrix(LOM_CI, NII)  = Me%Matrix(LOM_CI, NII) - DTDay           &
                                    * ImobilizationRateNO3  * Conversion             & 
                                    * (1./( (1./HetCN - 1./LCN )                     &
                                    + AnaPartition * (1./HetCN - 1./RCN)))  
                            
                                else

                  Me%Matrix(LOM_CI, PI)  = Me%Matrix(LOM_CI, PI) - DTDay             &
                                    * ImobilizationRateP  * Conversion               & 
                                    * (1./( (1./HetCP - 1./LCP )                     &
                                    + AnaPartition * (1./HetCP - 1./RCP)))

            end if
        end if
                
        
        !Sinks    : Anaerobic uptake, vai depender se respiram NO3 ou CO2 
        if (Me%ComputeImobilization) then

            if (NO3>NO3limit) then    ! respiram NO3

                Me%Matrix(LOM_CI, NII)     = Me%Matrix(LOM_CI, NII) - DTDay          & 
                                         * DenitrificationRate * Conversion        &
                                         * ( 0.1/( AnaPartition +1.) )
      
            else                   ! respiram CO2
                 
                Me%Matrix(LOM_CI, LOM_CI)     = Me%Matrix(LOM_CI, LOM_CI) - DTDay    & 
                                         *  methaneproduction         &
                                         * ( 0.1/( AnaPartition +1.) )

                Me%Matrix(LOM_CI, ROM_CI)     = Me%Matrix(LOM_CI, ROM_CI) - DTDay    & 
                                         *  methaneproduction                     &
                                         * ( 0.1/( AnaPartition +1.) )
                 
            end if              
        else
            
            if (NO3>NO3limit) then    ! respiram NO3

                Me%Matrix(LOM_CI, NII)     = Me%Matrix(LOM_CI, NII) - DTDay          & 
                                         * DenitrificationRate * Conversion        &
                                         * ( 0.1 )
      
            else                   ! respiram CO2
                 
                Me%Matrix(LOM_CI, LOM_CI)     = Me%Matrix(LOM_CI, LOM_CI) - DTDay    & 
                                         *  methaneproduction         &
                                         * ( 0.1 )

                Me%Matrix(LOM_CI, ROM_CI)     = Me%Matrix(LOM_CI, ROM_CI) - DTDay    & 
                                         *  methaneproduction                     &
                                         * ( 0.1 )
                 
            end if                      
        
        endif
        
        !Sources  : Anaerobic death

        if (.NOT. Me%Microorganisms%Anaerobic%LogicalMinumumPOP)                    &
            Me%Matrix(LOM_CI, AnaI) = Me%Matrix(LOM_CI, AnaI) + DTDay * AnaDeathRate
                                              
        !Sources  : Heterotrophic death

        if (.NOT. Me%Microorganisms%Heterotrophs%LogicalMinumumPOP )                &
            Me%Matrix(LOM_CI, HetI) = Me%Matrix(LOM_CI, HetI) + DTDay * HeteroDeathRate
                                              
        !Sources  : Autotrophic death

        if (.NOT. Me%Microorganisms%Autotrophs%LogicalMinumumPOP)                   &
            Me%Matrix(LOM_CI, AutoI)= Me%Matrix(LOM_CI, AutoI) + DTDay * AutoDeathRate


        if (Me%PropCalc%Sol_Bacteria) then
       
            !Sink: Sol uptake
            Me%Matrix(LOM_CI, PFI) = Me%Matrix(LOM_CI, PFI) - DTDay                       &
                                  * solubilizingRate * Conversion                       &
                                  * (0.1 /( (SolPartition + 1. ))) 

        end if

       
        !Independent term
        Me%IndTerm(LOM_CI) = Me%ExternalVar%Mass(LOM_CI, index)                                    
    !------------------------------------------------------------------------

    end subroutine SQLabilOrganicCarbon 
    !----------------------------------------------------------------------------

    !Refractary Organic Carbon
    !
    !SOURCES: - Nitrification eficiency 
    !SINKS:   - Death
    !----------------------------------------------------------------------------    
    subroutine SQRefractOrganicCarbon (index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN)         :: index

       !Local-------------------------------------------------------------------
        integer :: AMI, NII,ROM_CI,PI,PFI,LOM_CI            !!!Lúcia
       
        real    :: ImobilizationRateNH4
        real    :: DenitrificationRate
        real    :: ImobilizationRateNO3
        real    :: ImobilizationRateP         !!!Lúcia
        real    :: solubilizingrate
        real    :: ROM_C_Srate, RCN,RCP

        real    :: LCN
        real    :: LCP                          !!!Lúcia
        real    :: HeteroDeathRate, HetCN , HetCP    !!!Lúcia
        
        integer :: HetI                
        
        real    :: AnaDeathRate, AnaCN     !!!Lúcia  
        integer :: AnaI                

        real    :: AutoDeathRate, AutoCN !!!Lúcia
        integer :: AutoI               

        real    :: Partition, AnaPartition,seletion,solpartition    !!!Lúcia
        real    :: Conversion, DTDay

        real    :: AM,N,P,Real_N,Real_P         !!!Lúcia
        real    :: NO3,NO3limit,methaneproduction
        !------------------------------------------------------------------------
        
        ROM_C_Srate         = Me%SpecificRates%RefractOM_C%Value
        RCN                 = Me%RefractOM_CN_Ratio   
        ROM_CI              = Me%PropIndex%RefractOM_C
        LOM_CI              = Me%PropIndex%Labil_OM_C

        LCP                 = Me%LabilOM_CP_Ratio       !!!Lúcia  
        PFI                 = Me%PropIndex%Inorganic_P_fix  
            
        LCN                 = Me%LabiOM_CN_Ratio  
        RCP                 = Me%RefractOM_CP_Ratio         !!!Lúcia 
                
        DenitrificationRate     = Me%SpecificRates%NitrateToNgas%Value
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value
        ImobilizationRateP      = Me%SpecificRates%PhosphorusImobilization%Value  !!!Lúcia
        solubilizingrate        = Me%SpecificRates%Solubilizing%Value           
                      
        NII                 = Me%PropIndex%Nitrate
        AMI                 = Me%PropIndex%Ammonia
        PI                  = Me%PropIndex%Inorganic_P_soluble   !!!Lúcia

        Partition           = Me%Partition
        AnaPartition        = Me%AnaerobicPartition
        solpartition        = Me%Solpartition

        HeteroDeathRate     = Me%SpecificRates%Heterotrophs%Value
        HetI                = Me%PropIndex%HeterotrophicC
        HetCN               = Me%Microorganisms%Heterotrophs%CNRatio
        HetCP               = Me%Microorganisms%Heterotrophs%CPRatio        !!!Lúcia

        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaI                = Me%PropIndex%AnaerobicC
        AnaCN               = Me%Microorganisms%Anaerobic%CNRatio

        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutoI               = Me%PropIndex%AutotrophicC
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)
        seletion            = Me%Select                          !!!Lúcia       
        DTDay               = Me%DTDay


        AM   = Me%ExternalVar%Mass(AMI,index)   !!!
        N    = Me%ExternalVar%Mass(NII,index)   !!!
        P    = Me%ExternalVar%Mass(PI,index)    !!!
       
        NO3limit            = Me%NO3limit
        methaneproduction   = Me%SpecificRates%MethaneProduction%Value
        NO3                 = Me%ExternalVar%Mass(NII,index)


        
        Me%Matrix(ROM_CI, ROM_CI)   = 1.

        !Sinks: Heterotrophs uptake           
        

        if (seletion==5.OR.seletion==6 .OR. seletion==8 .OR. seletion==9 )  then 
        
             Me%Matrix(ROM_CI, ROM_CI)  = Me%Matrix(ROM_CI, ROM_CI) - DTDay * ROM_C_Srate  

        end if

        
        if (seletion==2.OR. seletion==3) then 


            Me%Matrix(ROM_CI, AMI)  = Me%Matrix(ROM_CI, AMI) - DTDay                & 
                                    * ImobilizationRateNH4  * Conversion            &
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)   &
                                    +  (1./HetCN - 1./RCN))) 
                
            Me%Matrix(ROM_CI, NII)  = Me%Matrix(ROM_CI, NII) - DTDay                &
                                    * ImobilizationRateNO3  * Conversion            & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)   &
                                    +  (1./HetCN - 1./RCN)))
         end if

        if (seletion==4.OR. seletion==7) then

            
            Me%Matrix(ROM_CI, PI)  = Me%Matrix(ROM_CI, PI) - DTDay                  &
                                    * ImobilizationRateP  * Conversion              & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)   &
                                    + (1./HetCP - 1./RCP)))

        end if

        
        if(seletion==1) then
    
            Real_N = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)*Conversion)  &
                        *(1./( (1./HetCN - 1./LCN )*(1/AnaPartition) +              & 
                        (1./HetCN - 1./RCN)))

            Real_P = ((ImobilizationRateP*P)*Conversion)                            &
                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)                   &
                    +  (1./HetCP - 1./RCP)))

    
            if( Real_N<Real_P) then

                  Me%Matrix(ROM_CI, AMI)  = Me%Matrix(ROM_CI, AMI) - DTDay          & 
                                    * ImobilizationRateNH4  * Conversion            &
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)   &
                                    +  (1./HetCN - 1./RCN))) 
                
                  Me%Matrix(ROM_CI, NII)  = Me%Matrix(ROM_CI, NII) - DTDay          &
                                    * ImobilizationRateNO3  * Conversion            & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)   &
                                    +  (1./HetCN - 1./RCN)))    
                            
                                else

                  Me%Matrix(ROM_CI, PI)  = Me%Matrix(ROM_CI, PI) - DTDay            &
                                    * ImobilizationRateP  * Conversion              & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)   &
                                    + (1./HetCP - 1./RCP)))

            end if
        end if


        !Sinks    : Anaerobic uptake, vai depender se respiram NO3 ou CO2 
        if (Me%ComputeImobilization) then
            if (NO3>NO3limit) then    ! respiram NO3


                Me%Matrix(ROM_CI, NII) = Me%Matrix(ROM_CI, NII) - DTDay                     &
                                   * DenitrificationRate * Conversion                   &
                                   * (0.1/( (1./AnaPartition) +1.) )     

            else                   ! respiram CO2

                Me%Matrix(ROM_CI, LOM_CI)     = Me%Matrix(ROM_CI, LOM_CI) - DTDay    & 
                                         *  methaneproduction                       &
                                         * ( 0.1/( (1/AnaPartition) +1.) )

                Me%Matrix(ROM_CI, ROM_CI)     = Me%Matrix(ROM_CI, ROM_CI) - DTDay    & 
                                         *  methaneproduction                        &
                                         * ( 0.1/( (1/AnaPartition) +1.) )
                 
            end if              
        else
        
            if (NO3>NO3limit) then    ! respiram NO3


                Me%Matrix(ROM_CI, NII) = Me%Matrix(ROM_CI, NII) - DTDay                     &
                                   * DenitrificationRate * Conversion                   &
                                   * (0.1 )     

            else                   ! respiram CO2

                Me%Matrix(ROM_CI, LOM_CI)     = Me%Matrix(ROM_CI, LOM_CI) - DTDay    & 
                                         *  methaneproduction                       &
                                         * ( 0.1 )

                Me%Matrix(ROM_CI, ROM_CI)     = Me%Matrix(ROM_CI, ROM_CI) - DTDay    & 
                                         *  methaneproduction                        &
                                         * ( 0.1 )
                 
            end if                      
        
        endif

                                             
        !Independent term
        Me%IndTerm(ROM_CI)     = Me%ExternalVar%Mass(ROM_CI, index)

        if (Me%PropCalc%Sol_Bacteria) then
        
            !Sink: Sol uptake
            Me%Matrix(ROM_CI, PFI) = Me%Matrix(ROM_CI, PFI) -DTDay                      &
                                  * solubilizingRate * Conversion                       &
                                  * (0.1 /( ((1/SolPartition) + 1. ))) 
        end if      
                                                 
    !----------------------------------------------------------------------------


    end subroutine SQRefractOrganicCarbon 
    !----------------------------------------------------------------------------


    !Heterotrophic Carbon
    !
    !SOURCES: - Aerobic labile and refractory decay,  
    !SINKS:   - Death, Respiration
    !----------------------------------------------------------------------------    
    subroutine SQHeterotrophicC (index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN)         :: index

        !Local-------------------------------------------------------------------
        real    :: HeteroDeathRate, HetCN, HetEF,HetCP          !!!Lúcia
        integer :: HetCI
        
        real    :: ROM_C_Srate, RCN,RCP      !!!Lúcia          
        integer :: ROM_CI                          

        real    :: LOM_C_Srate, LCN, LCP    !!!Lúcia            
        integer :: LOM_CI             
       
        real    :: ImobilizationRateNO3
        integer :: NII, AMI


        real    :: ImobilizationRateP       !!!Lúcia
        integer :: PI                   !!!Lúcia
        
        real    :: ImobilizationRateNH4

        real    :: Real_N_ROM , Real_P_ROM,Real_N_LOM, Real_P_LOM, AM,N, P   !!!Lúcia
        
        real    :: AnaPartition,partition, Conversion, DTDay,seletion          !!!Lúcia
        !------------------------------------------------------------------------

        HeteroDeathRate      = Me%SpecificRates%Heterotrophs%Value
        HetCI                = Me%PropIndex%HeterotrophicC
        HetCN                = Me%Microorganisms%Heterotrophs%CNRatio
        HetCP                = Me%Microorganisms%Heterotrophs%CPRatio !!!Lúcia
        HetEF                = Me%Microorganisms%Heterotrophs%EficiencyC

        ROM_C_Srate          = Me%SpecificRates%RefractOM_C%Value
        RCN                  = Me%RefractOM_CN_Ratio 
        RCP                  = Me%RefractOM_CP_Ratio  !!!Lúcia
        ROM_CI               = Me%PropIndex%RefractOM_C
            
        LOM_C_Srate          = Me%SpecificRates%Labil_OM_C%Value
        LCN                  = Me%LabiOM_CN_Ratio
        LCP                  = Me%LabilOM_CP_Ratio    !!!Lúcia
        LOM_CI               = Me%PropIndex%Labil_OM_C
        
        NII                  = Me%PropIndex%Nitrate
        ImobilizationRateNO3 = Me%SpecificRates%NitrateImobilization%Value
        N                    = Me%ExternalVar%Mass(NII, index)       !!!Lúcia

        AMI                  = Me%PropIndex%Ammonia
        ImobilizationRateNH4 = Me%SpecificRates%AmmoniaImobilization%Value
        AM                   = Me%ExternalVar%Mass(AMI, index)    !!!Lúcia


        PI                   = Me%PropIndex%Inorganic_P_soluble
        ImobilizationRateP   = Me%SpecificRates%PhosphorusImobilization%Value
        P                    = Me%ExternalVar%Mass(PI, index)     !!!Lúcia

        AnaPartition         = Me%AnaerobicPartition
        Partition            = Me%Partition
        Conversion           = Me%ExternalVar%DissolvedToParticulate (index)

        DTDay                = Me%DTDay
        seletion             = Me%Select    !!!Lúcia

        
        Me%Matrix(HetCI, HetCI) = 1.
        
        !Sink   : Heterotrophic death

        if (.NOT. Me%Microorganisms%Heterotrophs%LogicalMinumumPOP)                 &
           Me%Matrix(HetCI, HetCI) = Me%Matrix(HetCI, HetCI) - DTDay * HeteroDeathRate


        !Souce  : Aerobic labile and refractory decay
        !Sink   : Breath


        if (seletion==5.OR.seletion==6 .OR. seletion==8 .OR. seletion==9 )  then 
        
            Me%Matrix(HetCI, LOM_CI)  = Me%Matrix(HetCI, LOM_CI) +DTDay                &
                                         * LOM_C_Srate  

            Me%Matrix(HetCI, ROM_CI)  = Me%Matrix(HetCI, ROM_CI) +DTDay                &
                                         * ROM_C_Srate 
        
             !Breath
            Me%Matrix(HetCI, LOM_CI)  = Me%Matrix(HetCI, LOM_CI) -DTDay *              &
                                         LOM_C_Srate*(HetEF)

            Me%Matrix(HetCI, ROM_CI)  = Me%Matrix(HetCI, ROM_CI) -DTDay *              &
                                         ROM_C_Srate*(HetEF)
               
        end if

        if (seletion==2.OR. seletion==3) then 

            Me%Matrix(HetCI, AMI)  = Me%Matrix(HetCI, AMI) + DTDay                      & 
                                   * ImobilizationRateNH4  * Conversion                &
                                   * (1/((1/HetCN-1/LCN)*(1/Anapartition)              &
                                   +(1/HetCN-1/RCN)))

            Me%Matrix(HetCI, AMI)  = Me%Matrix(HetCI, AMI) + DTDay                      & 
                                   * ImobilizationRateNH4  * Conversion                &
                                   * (1/((1/HetCN-1/LCN)+                              &
                                   Anapartition*(1/HetCN-1/RCN)))


            Me%Matrix(HetCI, NII)  = Me%Matrix(HetCI, NII) + DTDay                      &
                                   * ImobilizationRateNO3  * Conversion                &
                                   * (1/((1/HetCN-1/LCN)+                              & 
                                   Anapartition*(1/HetCN-1/RCN)))                           


            Me%Matrix(HetCI, NII)  = Me%Matrix(HetCI, NII) + DTDay                      &
                                   * ImobilizationRateNO3  * Conversion                &
                                   * (1/((1/HetCN-1/LCN)*(1/Anapartition)              &
                                   +(1/HetCN-1/RCN)))                           


            !Breath

            Me%Matrix(HetCI, AMI)  = Me%Matrix(HetCI, AMI) - DTDay                     &
                                   * ImobilizationRateNH4 * Conversion * ( HetEF)      &
                                   * (1/((1/HetCN-1/LCN)*(1/Anapartition)              &
                                   +(1/HetCN-1/RCN))) 

            Me%Matrix(HetCI, AMI)  = Me%Matrix(HetCI, AMI) - DTDay                     &
                                   * ImobilizationRateNH4 * Conversion * ( HetEF)      &
                                   * (1/((1/HetCN-1/LCN)             &
                                   +Anapartition*(1/HetCN-1/RCN))) 
                
            Me%Matrix(HetCI, NII)  = Me%Matrix(HetCI, NII) - DTDay                     & 
                                   * ImobilizationRateNO3 * Conversion * ( HetEF)      & 
                                   * (1/((1/HetCN-1/LCN)              &
                                   +Anapartition*(1/HetCN-1/RCN)))

            Me%Matrix(HetCI, NII)  = Me%Matrix(HetCI, NII) - DTDay                     & 
                                   * ImobilizationRateNO3 * Conversion * ( HetEF)      & 
                                   * (1/((1/HetCN-1/LCN)*(1/Anapartition)              &
                                   +(1/HetCN-1/RCN)))
        
         end if


        if (seletion==4.OR. seletion==7) then

            
            Me%Matrix(HetCI, PI)  = Me%Matrix(HetCI, PI) + DTDay                      &
                                    * ImobilizationRateP  * Conversion                & 
                                    * (1./((1./HetCP - 1./LCP )*(1/Anapartition)      & 
                                    + (1./HetCP - 1./RCP))) 
                                    
            Me%Matrix(HetCI, PI)  = Me%Matrix(HetCI, PI) + DTDay                      &
                                    * ImobilizationRateP  * Conversion                & 
                                    * (1./((1./HetCP - 1./LCP )                       & 
                                    + Anapartition*(1./HetCP - 1./RCP)))    
                                                                                    
                                    
            !Breath

            Me%Matrix(HetCI, PI)  = Me%Matrix(HetCI, PI) - DTDay                      &
                                    * ImobilizationRateP  * Conversion*(HetEF)        & 
                                    * (1./((1./HetCP - 1./LCP )*(1/Anapartition)      & 
                                    + (1./HetCP - 1./RCP))) 
                        
            Me%Matrix(HetCI, PI)  = Me%Matrix(HetCI, PI) - DTDay                      &
                                    * ImobilizationRateP  * Conversion*(HetEF)        & 
                                    * (1./((1./HetCP - 1./LCP )                       & 
                                    + Anapartition*(1./HetCP - 1./RCP)))    
            
                                                                
        end if

        

        if(seletion==1) then
    
            Real_N_ROM = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)            &
                         *Conversion)*(1./( (1./HetCN - 1./LCN )                      &
                         *(1/AnaPartition) +  (1./HetCN - 1./RCN)))

            Real_P_ROM = ((ImobilizationRateP*P)*Conversion)                          &
                         * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)                &
                         +  (1./HetCP - 1./RCP)))



            Real_N_LOM = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)            &
                         *Conversion) *(1./( (1./HetCN - 1./LCN )                     &
                         +  AnaPartition*(1./HetCN - 1./RCN)))


            Real_P_LOM = ((ImobilizationRateP*P)*Conversion)                          &
                         * (1./( (1./HetCP - 1./LCP )                                 &
                         +  AnaPartition*(1./HetCP - 1./RCP)))


    
            if( Real_N_ROM<Real_P_ROM) then

                Me%Matrix(HetCI, AMI)  = Me%Matrix(HetCI, AMI) + DTDay               & 
                                    * ImobilizationRateNH4  * Conversion               &
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)      &
                                    +  (1./HetCN - 1./RCN))) 
                
                Me%Matrix(HetCI, NII)  = Me%Matrix(HetCI, NII) + DTDay               &
                                    * ImobilizationRateNO3  * Conversion               & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)      & 
                                    +  (1./HetCN - 1./RCN)))    
                !Breath 
                
                        
                Me%Matrix(HetCI, AMI)  = Me%Matrix(HetCI, AMI) - DTDay               & 
                                    * ImobilizationRateNH4  * Conversion*(HetEF)       &
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)      &
                                    +  (1./HetCN - 1./RCN))) 
                
                Me%Matrix(HetCI, NII)  = Me%Matrix(HetCI, NII) - DTDay               &
                                    * ImobilizationRateNO3  * Conversion *(HetEF)      & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)      & 
                                    +  (1./HetCN - 1./RCN)))    


            else

                Me%Matrix(HetCI, PI)  = Me%Matrix(HetCI, PI) + DTDay                   &
                                    * ImobilizationRateP  * Conversion                 & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)      &
                                    + (1./HetCP - 1./RCP)))


                !Breath

                Me%Matrix(HetCI, PI)  = Me%Matrix(HetCI, PI) - DTDay                  &
                                    * ImobilizationRateP  * Conversion*(HetEF)         & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)      &
                                    + (1./HetCP - 1./RCP)))

            end if

            if( Real_N_LOM<Real_P_LOM) then

                Me%Matrix(HetCI, AMI)  = Me%Matrix(HetCI, AMI) + DTDay                & 
                                    * ImobilizationRateNH4  * Conversion                &
                                    * (1./( (1./HetCN - 1./LCN )                        &
                                    + AnaPartition* (1./HetCN - 1./RCN))) 
                
                Me%Matrix(HetCI, NII)  = Me%Matrix(HetCI, NII) + DTDay                &
                                    * ImobilizationRateNO3  * Conversion                & 
                                    * (1./( (1./HetCN - 1./LCN )                        & 
                                    +  AnaPartition*(1./HetCN - 1./RCN)))   
                !Breath 
                
                        
                Me%Matrix(HetCI, AMI)  = Me%Matrix(HetCI, AMI) - DTDay                & 
                                    * ImobilizationRateNH4  * Conversion*(HetEF)        &
                                    * (1./( (1./HetCN - 1./LCN )                        &
                                    +  AnaPartition*(1./HetCN - 1./RCN))) 
                
                Me%Matrix(HetCI, NII)  = Me%Matrix(HetCI, NII) - DTDay                &
                                    * ImobilizationRateNO3  * Conversion *(HetEF)       & 
                                    * (1./( (1./HetCN - 1./LCN )                        & 
                                    + AnaPartition*(1./HetCN - 1./RCN)))    


            else


                Me%Matrix(HetCI, PI) = Me%Matrix(HetCI, PI) + DTDay                   &
                                    * ImobilizationRateP  * Conversion                  & 
                                    * (1./( (1./HetCP - 1./LCP )                        &
                                    + AnaPartition*(1./HetCP - 1./RCP)))

                 !Breath 


                Me%Matrix(HetCI, PI)  = Me%Matrix(HetCI, PI) - DTDay                  &
                                    * ImobilizationRateP  * Conversion*(HetEF)          & 
                                    * (1./( (1./HetCP - 1./LCP )                        &
                                    + AnaPartition*(1./HetCP - 1./RCP)))

            end if
    
        end if
            
            !Independent term
        Me%IndTerm(HetCI)       = Me%ExternalVar%Mass(HetCI, index) 
    !------------------------------------------------------------------------

    end subroutine SQHeterotrophicC 
    !----------------------------------------------------------------------------


    !Autotrophic Carbon
    !
    !SOURCES: - Nitrification eficiency 
    !SINKS:   - Death
    !----------------------------------------------------------------------------    
    subroutine SQAutotrophicC (index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN)         :: index

        !Local-------------------------------------------------------------------
        real    :: DTDay
          
        real    :: AutoCN, AutoEf, AutoDeathRate
        integer :: AutoCI
      
        real    :: NitificationRate, Conversion
        
        integer :: AMI    
        !------------------------------------------------------------------------
        
        DTDay               = Me%DTDay

        AutoCI              = Me%PropIndex%AutotrophicC
        AMI                 = Me%PropIndex%Ammonia
        
        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio
        
        NitificationRate    = Me%SpecificRates%AmmoniaToNitrate%Value
        AutoEf              = Me%Microorganisms%Autotrophs%EficiencyN

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)

        Me%Matrix(AutoCI, AutoCI)   =   1.
        
        !Sink: Death
        if (.NOT.Me%Microorganisms%Autotrophs%LogicalMinumumPOP)                    &
            Me%Matrix(AutoCI, AutoCI) = Me%Matrix(AutoCI, AutoCI) - DTDay * AutoDeathRate
        
        !Sources: Nitrification (CO2 uptake)
        Me%Matrix(AutoCI, AMI) = Me%Matrix(AutoCI, AMI) + DTDay                    &
                               * NitificationRate * Conversion * AutoEf    !*AutoCN
        
        !Independent term
        Me%IndTerm(AutoCI)     = Me%ExternalVar%Mass(AutoCI, index) 
    !------------------------------------------------------------------------

    end subroutine SQAutotrophicC 
    !----------------------------------------------------------------------------


    !Anaerobic Carbon
    !
    !SOURCES: - Mineralization, CO2 assimilation
    !SINKS:   - Death, respiration
    !----------------------------------------------------------------------------    
    subroutine SQAnaerobicC (index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN)         :: index
        !------------------------------------------------------------------------

        real    :: AnaCN, AnaCEf, AnaDeathRate
        integer :: AnaCI,LOM_CI,ROM_CI
        real    :: DenitrificationRate, Conversion ,AnaPartition,methaneproduction
        integer :: NII 
        real    :: DTDay   
        real    :: NO3limit
        real    :: Eficiencia
        real    :: NO3 
        !------------------------------------------------------------------------
        
        DTDay               = Me%DTDay

        AnaCI               = Me%PropIndex%AnaerobicC
        NII                 = Me%PropIndex%Nitrate
        LOM_CI              = Me%PropIndex%Labil_OM_C
        ROM_CI              = Me%PropIndex%RefractOM_C
        NO3                 = Me%ExternalVar%Mass(NII,index)
        
        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaCN               = Me%Microorganisms%Anaerobic%CNRatio
        
        DenitrificationRate = Me%SpecificRates%NitrateToNgas%Value
        AnaCEf              = Me%Microorganisms%Anaerobic%EficiencyC


        AnaPartition        = Me%AnaerobicPartition
        Conversion          = Me%ExternalVar%DissolvedToParticulate (index) 
        
        NO3limit            = Me%NO3limit
        methaneproduction   = Me%SpecificRates%MethaneProduction%Value
       
        Eficiencia = 0.1        

        Me%Matrix(AnaCI, AnaCI) = 1.

        !Sink   : Death
        if (.NOT.Me%Microorganisms%Anaerobic%LogicalMinumumPOP )                    &
            Me%Matrix(AnaCI, AnaCI) = Me%Matrix(AnaCI, AnaCI) - DTDay * AnaDeathRate
        

! Agora vamos ter dois casos: se existe nitrato para ser respirado ou se respiram CO2


        if (NO3>NO3limit)   then

            !Sources: degração de matéria orgãnica labile e refractaria
            Me%Matrix(AnaCI, NII) = Me%Matrix(AnaCI, NII) + DTDay * DenitrificationRate &
                                  * Conversion *0.1 


            ! Sink:  libertaçao de CO2, segundo a respiração de NO3

            Me%Matrix(AnaCI, NII) = Me%Matrix(AnaCI, NII) - DTDay * DenitrificationRate &
                                  * Conversion  *0.1* AnaCEf

        else

            !Source :  degradação de matéria orgânica e respiração de CO2

                !( vem da labile)    
            Me%Matrix(AnaCI, LOM_CI)     = Me%Matrix(AnaCI, LOM_CI) + DTDay    & 
                                     *  methaneproduction       &
                                     *  0.1

                ! (vem da refractory)

            Me%Matrix(AnaCI, ROM_CI)     = Me%Matrix(AnaCI, ROM_CI) + DTDay    & 
                                     *  methaneproduction         &
                                     *  0.1


            ! Source : assimilaçao directa de CO2

            Me%Matrix(AnaCI, LOM_CI)     = Me%Matrix(AnaCI, LOM_CI) + DTDay    & 
                                     *  methaneproduction       &
                                     *  0.1

            Me%Matrix(AnaCI, ROM_CI)     = Me%Matrix(AnaCI, ROM_CI) + DTDay    & 
                                     *  methaneproduction         &
                                     *  0.1



            ! sink: excreção de CO2 e CH4, após a respiração de CO2

            Me%Matrix(AnaCI, LOM_CI)     = Me%Matrix(AnaCI, LOM_CI) - DTDay    & 
                                        *  methaneproduction      &
                                        *0.1*AnaCEf

            Me%Matrix(AnaCI, ROM_CI)     = Me%Matrix(AnaCI, ROM_CI) - DTDay    & 
                                        *  methaneproduction     &
                                        *0.1*AnaCEf

        end if
                                    
                                   
        !Independent term
        Me%IndTerm(AnaCI)     = Me%ExternalVar%Mass(AnaCI, index) 


    !------------------------------------------------------------------------

    end subroutine SQAnaerobicC 
    !----------------------------------------------------------------------------

    !CO2
    !
    !SOURCES: - Aerobic and anaerobic Excretation
    !
    !----------------------------------------------------------------------------    
    subroutine SQCO2 (index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN)         :: index
        !------------------------------------------------------------------------

        integer     :: ICO2
        real        :: HetCN, HetEF,HetCP           !!!Lúcia

        real        :: ROM_C_Srate, RCN,RCP      !!!Lúcia          
        integer     :: ROM_CI                          

        real        :: LOM_C_Srate, LCN, LCP    !!!Lúcia            
        integer     :: LOM_CI             
       
        real        :: ImobilizationRateNO3
        integer     :: NII, AMI
        real        :: DenitrificationRate
        real        :: AnaCEf
        real        :: ImobilizationRateP       !!!Lúcia
        integer     :: PI                   !!!Lúcia
        
        real        :: ImobilizationRateNH4

        real        :: Real_N_ROM , Real_P_ROM,Real_N_LOM, Real_P_LOM, AM,N, P   !!!Lúcia


        real        :: AnaPartition,partition, Conversion, DTDay,seletion          !!!Lúcia
        real        :: NO3,NO3limit,methaneproduction
        !------------------------------------------------------------------------
        
        HetCN                = Me%Microorganisms%Heterotrophs%CNRatio
        HetCP                = Me%Microorganisms%Heterotrophs%CPRatio !!!Lúcia
        HetEF                = Me%Microorganisms%Heterotrophs%EficiencyC
        ROM_C_Srate          = Me%SpecificRates%RefractOM_C%Value
        RCN                  = Me%RefractOM_CN_Ratio 
        RCP                  = Me%RefractOM_CP_Ratio  !!!Lúcia
        ROM_CI               = Me%PropIndex%RefractOM_C
        LOM_C_Srate          = Me%SpecificRates%Labil_OM_C%Value
        LCN                  = Me%LabiOM_CN_Ratio
        LCP                  = Me%LabilOM_CP_Ratio    !!!Lúcia
        LOM_CI               = Me%PropIndex%Labil_OM_C
        
        NII                  = Me%PropIndex%Nitrate
        ImobilizationRateNO3 = Me%SpecificRates%NitrateImobilization%Value
        N                    = Me%ExternalVar%Mass(NII, index)       !!!Lúcia

        AMI                  = Me%PropIndex%Ammonia
        ImobilizationRateNH4 = Me%SpecificRates%AmmoniaImobilization%Value
        AM                   = Me%ExternalVar%Mass(AMI, index)    !!!Lúcia


        PI                   = Me%PropIndex%Inorganic_P_soluble
        ImobilizationRateP   = Me%SpecificRates%PhosphorusImobilization%Value
        P                    = Me%ExternalVar%Mass(PI, index)     !!!Lúcia

        DenitrificationRate  = Me%SpecificRates%NitrateToNgas%Value
        AnaCEf               = Me%Microorganisms%Anaerobic%EficiencyC


        AnaPartition         = Me%AnaerobicPartition
        Partition            = Me%Partition
        Conversion           = Me%ExternalVar%DissolvedToParticulate (index)

        DTDay                = Me%DTDay
        seletion             = Me%Select    !!!Lúcia

        ICO2        = Me%PropIndex%CO2


        NO3limit            = Me%NO3limit
        methaneproduction   = Me%SpecificRates%MethaneProduction%Value
        NO3                 = Me%ExternalVar%Mass(NII,index)

       
         Me%Matrix(ICO2, ICO2) =  1.


        ! Source: aerobic Excretation
        if (seletion==5.OR.seletion==6 .OR. seletion==8 .OR. seletion==9 )  then 
        
        
            Me%Matrix(ICO2, LOM_CI)  = Me%Matrix(ICO2, LOM_CI) + DTDay *            &
                                         LOM_C_Srate*(HetEF)

            Me%Matrix(ICO2, ROM_CI)  = Me%Matrix(ICO2, ROM_CI) + DTDay *            &
                                         ROM_C_Srate*(HetEF)
               
                  
        end if

        if (seletion==2.OR. seletion==3) then 

            Me%Matrix(ICO2, AMI)  = Me%Matrix(ICO2, AMI) + DTDay                     &
                                   * ImobilizationRateNH4 * Conversion * ( HetEF)      &
                                   * (1/((1/HetCN-1/LCN)*(1/Anapartition)              &
                                   +(1/HetCN-1/RCN))) 

            Me%Matrix(ICO2, AMI)  = Me%Matrix(ICO2, AMI) + DTDay                     &
                                   * ImobilizationRateNH4 * Conversion * ( HetEF)      &
                                   * (1/((1/HetCN-1/LCN)             &
                                   +Anapartition*(1/HetCN-1/RCN))) 
                
            Me%Matrix(ICO2, NII)  = Me%Matrix(ICO2, NII) + DTDay                     & 
                                   * ImobilizationRateNO3 * Conversion * ( HetEF)      & 
                                   * (1/((1/HetCN-1/LCN)              &
                                   +Anapartition*(1/HetCN-1/RCN)))

            Me%Matrix(ICO2, NII)  = Me%Matrix(ICO2, NII) + DTDay                     & 
                                   * ImobilizationRateNO3 * Conversion * ( HetEF)      & 
                                   * (1/((1/HetCN-1/LCN)*(1/Anapartition)              &
                                   +(1/HetCN-1/RCN)))
        
         end if


        if (seletion==4.OR. seletion==7) then

            

            Me%Matrix(ICO2, PI)  = Me%Matrix(ICO2, PI) + DTDay                      &
                                    * ImobilizationRateP  * Conversion*(HetEF)        & 
                                    * (1./((1./HetCP - 1./LCP )*(1/Anapartition)      & 
                                    + (1./HetCP - 1./RCP))) 
                        
            Me%Matrix(ICO2, PI)  = Me%Matrix(ICO2, PI) + DTDay                      &
                                    * ImobilizationRateP  * Conversion*(HetEF)        & 
                                    * (1./((1./HetCP - 1./LCP )                       & 
                                    + Anapartition*(1./HetCP - 1./RCP)))    
            
                                                                
        end if

        

        if(seletion==1) then
    
            Real_N_ROM = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)            &
                         *Conversion)*(1./( (1./HetCN - 1./LCN )                      &
                         *(1/AnaPartition) +  (1./HetCN - 1./RCN)))

            Real_P_ROM = ((ImobilizationRateP*P)*Conversion)                          &
                         * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)                &
                         +  (1./HetCP - 1./RCP)))



            Real_N_LOM = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)            &
                         *Conversion) *(1./( (1./HetCN - 1./LCN )                     &
                         +  AnaPartition*(1./HetCN - 1./RCN)))


            Real_P_LOM = ((ImobilizationRateP*P)*Conversion)                          &
                         * (1./( (1./HetCP - 1./LCP )                                 &
                         +  AnaPartition*(1./HetCP - 1./RCP)))


    
            if( Real_N_ROM<Real_P_ROM) then

                Me%Matrix(ICO2, AMI)  = Me%Matrix(ICO2, AMI) + DTDay               & 
                                    * ImobilizationRateNH4  * Conversion*(HetEF)       &
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)      &
                                    +  (1./HetCN - 1./RCN))) 
                
                Me%Matrix(ICO2, NII)  = Me%Matrix(ICO2, NII) + DTDay               &
                                    * ImobilizationRateNO3  * Conversion *(HetEF)      & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)      & 
                                    +  (1./HetCN - 1./RCN)))    


                                else

                Me%Matrix(ICO2, PI)  = Me%Matrix(ICO2, PI) + DTDay                  &
                                    * ImobilizationRateP  * Conversion*(HetEF)         & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)      &
                                    + (1./HetCP - 1./RCP)))

            end if

            if( Real_N_LOM<Real_P_LOM) then

                
                        
                Me%Matrix(ICO2, AMI)  = Me%Matrix(ICO2, AMI) + DTDay                & 
                                    * ImobilizationRateNH4  * Conversion*(HetEF)        &
                                    * (1./( (1./HetCN - 1./LCN )                        &
                                    +  AnaPartition*(1./HetCN - 1./RCN))) 
                
                Me%Matrix(ICO2, NII)  = Me%Matrix(ICO2, NII) + DTDay                &
                                    * ImobilizationRateNO3  * Conversion *(HetEF)       & 
                                    * (1./( (1./HetCN - 1./LCN )                        & 
                                    + AnaPartition*(1./HetCN - 1./RCN)))    


            else


                Me%Matrix(ICO2, PI)  = Me%Matrix(ICO2, PI) + DTDay                  &
                                    * ImobilizationRateP  * Conversion*(HetEF)          & 
                                    * (1./( (1./HetCP - 1./LCP )                        &
                                    + AnaPartition*(1./HetCP - 1./RCP)))

            end if
    
        end if

        ! Source : anaerobic Excretation

        if (NO3> NO3limit) then 

            Me%Matrix(ICO2, NII) = Me%Matrix(ICO2, NII) + DTDay * DenitrificationRate &
                              * Conversion  *0.1* AnaCEf

        else

            Me%Matrix(ICO2, LOM_CI)     = Me%Matrix(ICO2, LOM_CI) + DTDay    & 
                                        *  methaneproduction     &
                                        *0.1*AnaCEf

            Me%Matrix(ICO2, ROM_CI)     = Me%Matrix(ICO2, ROM_CI) + DTDay    & 
                                        *  methaneproduction     &
                                        *0.1*AnaCEf

        end if


        !Independent term
        Me%IndTerm(ICO2)     = Me%ExternalVar%Mass(ICO2, index) 


    end subroutine SQCO2
    !----------------------------------------------------------------------------


    subroutine SQMethane(index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN) :: index
        !------------------------------------------------------------------------

        integer  ::    ICH4,LOM_CI,ROM_CI,NII
        real     ::    methaneproduction
        real     ::    DTDay
        real     ::    NO3,NO3limit
        real     ::    conversion
        
        !------------------------------------------------------------------------


        ICH4                = Me%PropIndex%methane
        LOM_CI              = Me%PropIndex%Labil_OM_C
        ROM_CI              = Me%PropIndex%RefractOM_C
        methaneproduction   = Me%SpecificRates%MethaneProduction%Value
        DTDay               = Me%DTDay
        NII                 = Me%PropIndex%Nitrate
        NO3                 = Me%ExternalVar%Mass(NII,index)
        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)



        NO3limit            = Me%NO3limit

            
        Me%Matrix(ICH4, ICH4) =  1.

        if (NO3<=NO3limit) then 


            ! source: anaerobic respiration when nitrate is not available
            Me%Matrix(ICH4, LOM_CI)     = Me%Matrix(ICH4, LOM_CI) + DTDay    & 
                                        *  methaneproduction     &
                                        *(1-0.1)

            Me%Matrix(ICH4, ROM_CI)     = Me%Matrix(ICH4, ROM_CI) + DTDay    & 
                                        *  methaneproduction     &
                                        *(1-0.1)

        end if
                                
        !Independent term
        Me%IndTerm(ICH4)        = Me%ExternalVar%Mass(ICH4, index) 

    end subroutine SQMethane

    !----------------------------------------------------------------------------
 
 
    subroutine SQNitrogen(index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN) :: index
        !------------------------------------------------------------------------

        call SQAmmonia                (index)

        call SQNitrate                (index)
                                                        
        call SQLabilOrganicNitrogen   (index)

        call SQRefractOrganicNitrogen (index)

        call SQHeterotrophicN         (index)

        call SQAutotrophicN           (index)

        call SQAnaerobicN             (index)
        
        call SQNgas                   (index)

        call SQUrea                   (index)

!       call SQAmmoniaGas              (index)
    !------------------------------------------------------------------------

    end subroutine SQNitrogen
    !---------------------------------------------------------------------------


    !AMMONIA 
    !
    !SOURCES    : Heterophs
    !SINKS      : Heterophs, NO3 (Autotrophs), N gas (Anaerobic), Plants (to_do)                        
    !----------------------------------------------------------------------------       
    subroutine SQAmmonia(index)
    
        !Arguments---------------------------------------------------------------
        integer, intent(IN)         :: index            
        
        !Local-------------------------------------------------------------------
        integer :: AMI, NII ,PI,PFI,IUrea         !Lúcia
        real    :: AM, N , P              !!!Lúcia 
        integer :: LOM_CI
        real    :: LOM_C_Srate, LCN, LCP, CurrentLabilC  !!!Lúcia

        integer :: ROM_CI 
        real    :: ROM_C_Srate, RCN,RCP, CurrentRefractC
        
        real    :: solCN,SolCEf
        real    :: HetCN, HetEF, AnaCN, AnaNEF, AnaCEf, AutoEf
        real    :: HetCP, AnaCP    !!!Lúcia
        
        real    :: NitrificationRate, DenitrificationRate
        real    :: ImobilizationRateNH4, ImobilizationRateNO3, ImobilizationRateP       
        real    :: Partition, AnaPartition
        real    :: Real_N_LOM, Real_P_LOM, Real_N_ROM, Real_P_ROM   !!!Lúcia
        real    :: Conversion, DTDay, seletion     !!!Lúcia
        real    :: solubilizingrate,SolPartition
        real    :: NO3,NO3limit,MethaneProduction
        
        real    :: wind, XKG,TF,EK,PNH3,PANH3,TK,XK1,XK
        real    :: Kvol,temp,H,ureahydrolysis      
                        
        !------------------------------------------------------------------------
        
        AMI                   = Me%PropIndex%Ammonia
        NII                   = Me%PropIndex%Nitrate
        PI                    = Me%PropIndex%Inorganic_P_soluble             !!!Lúcia
        PFI                   = Me%PropIndex%Inorganic_P_fix
        IUrea                 = Me%PropIndex%Urea                   
        HetCN                 = Me%Microorganisms%Heterotrophs%CNRatio
        HetCP                 = Me%Microorganisms%Heterotrophs%CPRatio     !!!Lúcia

        HetEF                 = Me%Microorganisms%Heterotrophs%EficiencyC
        AnaCN                 = Me%Microorganisms%Anaerobic%CNRatio
        AnaCP                 = Me%Microorganisms%Anaerobic%CPRatio         !!!Lúcia
        AnaCEf                = Me%Microorganisms%Anaerobic%EficiencyC
        solCEf                = Me%Microorganisms%Heterotrophs%EficiencyC
        SolCN                 = Me%Microorganisms%Sols%CNRatio

        AnaNEF                = Me%Microorganisms%Anaerobic%EficiencyN       
        AutoEf               = Me%Microorganisms%Autotrophs%EficiencyN
        
        LOM_CI                = Me%PropIndex%Labil_OM_C
        LOM_C_Srate           = Me%SpecificRates%Labil_OM_C%Value
        LCN                   = Me%LabiOM_CN_Ratio
        LCP                   = Me%LabilOM_CP_Ratio                          !!!Lúcia
        CurrentLabilC         = Me%ExternalVar%Mass(LOM_CI, index)

        ROM_CI                = Me%PropIndex%RefractOM_C        
        ROM_C_Srate           = Me%SpecificRates%RefractOM_C%Value
        RCN                   = Me%RefractOM_CN_Ratio
        RCP                   = Me%RefractOM_CP_Ratio
        CurrentRefractC       = Me%ExternalVar%Mass(ROM_CI, index)
       
        AM                    = Me%ExternalVar%Mass(AMI, index)    !!!Lúcia
        N                     = Me%ExternalVar%Mass(NII, index)       !!!Lúcia
        ImobilizationRateP    = Me%SpecificRates%PhosphorusImobilization%Value
        P                     = Me%ExternalVar%Mass(PI, index)     !!!Lúcia

        ImobilizationRateNH4  = Me%SpecificRates%AmmoniaImobilization%Value
        ImobilizationRateNO3  = Me%SpecificRates%NitrateImobilization%Value
        solubilizingrate      = Me%SpecificRates%Solubilizing%Value   
        NitrificationRate   = Me%SpecificRates%AmmoniaToNitrate%Value 
        DenitrificationRate = Me%SpecificRates%NitrateToNgas%Value   

        seletion              = Me%Select    !!!Lúcia       
        Partition             = Me%Partition
        AnaPartition          = Me%AnaerobicPartition
        Conversion            = Me%ExternalVar%DissolvedToParticulate (index)        
        solpartition          = Me%Solpartition

        DTDay                 = Me%DTDay

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)
        methaneproduction   = Me%SpecificRates%MethaneProduction%Value
        NO3                 = Me%ExternalVar%Mass(NII,index)
        NO3limit            = Me%NO3limit       
        ureahydrolysis      = Me%SpecificRates%ureahydrolysis%Value
    
        wind    =   Me%ExternalVar%Wind(index)
        Temp    =   Me%ExternalVar%Temperature(index)
        H       =   Me%Hydrogen




        !------------------------------------------------------------------------
        
      
        Me%Matrix(AMI, AMI) = 1.   
          

        ! Sink : Nitrification

        !Me%Matrix(AMI, AMI) = Me%Matrix(AMI, AMI) - DTDay * NitrificationRate 
        !The above is not consistent with nitrate routine
        Me%Matrix(AMI, AMI) = Me%Matrix(AMI, AMI) - DTDay * NitrificationRate* (1. - AutoEf)
    
       ! Sources : Anaerobic excretion
       ! Agora vamos ter dois casos: se existe nitrato para ser respirado ou se respiram CO2


        if (NO3>NO3limit)   then

            Me%Matrix(AMI, NII) = Me%Matrix(AMI, NII) + DTDay                &
                              * DenitrificationRate*0.1                  & 
                              * AnaCEf/ AnaCN                                    

        else 


            Me%Matrix(AMI, LOM_CI)     = Me%Matrix(AMI, LOM_CI) + DTDay    & 
                                        *  (methaneproduction/conversion)      &
                                        *0.1*AnaCEf/AnaCN

            Me%Matrix(AMI, ROM_CI)     = Me%Matrix(AMI, ROM_CI) + DTDay       & 
                                        *  (methaneproduction/conversion)    &
                                        *0.1*AnaCEf/AnaCN

        end If 
 



        !Sources : Heterotrophs excretion

        if (seletion==5.OR.seletion==6 .OR. seletion==8 .OR. seletion==9 )  then 


            Me%Matrix(AMI, LOM_CI)  = Me%Matrix(AMI, LOM_CI) +DTDay                  &
                                       * LOM_C_Srate* (HetEF/HetCN)*(1/conversion)

            Me%Matrix(AMI, ROM_CI)  = Me%Matrix(AMI, ROM_CI) +DTDay                  &
                                       * ROM_C_Srate* (HetEF/HetCN)*(1/conversion)
                      
        end if

        if (seletion==2.OR. seletion==3) then 


            Me%Matrix(AMI, AMI)  = Me%Matrix(AMI, AMI) + DTDay                         &   
                                   * ImobilizationRateNH4  *(HetEF/HetCN)             &
                                   * 1./((1./HetCN - 1./LCN)*(1/Anapartition)         &
                                   + (1./HetCN - 1./RCN))               
                                   
            Me%Matrix(AMI, AMI)  = Me%Matrix(AMI, AMI) + DTDay                         &   
                                   * ImobilizationRateNH4  *(HetEF/HetCN)             &
                                   * 1./((1./HetCN - 1./LCN)        &
                                   + anapartition*(1./HetCN - 1./RCN))               


            Me%Matrix(AMI, NII)  = Me%Matrix(AMI, NII) + DTDay                         &
                                   * ImobilizationRateNO3 *(HetEF/HetCN)              &
                                   * 1./((1./HetCN - 1./LCN)*(1/anapartition)         &
                                   + (1./HetCN - 1./RCN))      
                                   
            Me%Matrix(AMI, NII)  = Me%Matrix(AMI, NII) + DTDay                         &
                                   * ImobilizationRateNO3 *(HetEF/HetCN)              &
                                   * 1./((1./HetCN - 1./LCN)                          &
                                   + Anapartition*(1./HetCN - 1./RCN))      
                                   
                                   
                                   
                                            
                                                                                 
        
        end if

        if (seletion==4.OR. seletion==7) then

            Me%Matrix(AMI, PI)  = Me%Matrix(AMI, PI) + DTDay                         &
                                    * ImobilizationRateP *( HetEF/HetCN)             & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/anapartition)    &                     
                                    + (1./HetCP - 1./RCP)))                         

            Me%Matrix(AMI, PI)  = Me%Matrix(AMI, PI) + DTDay                         &
                                    * ImobilizationRateP *( HetEF/HetCN)             & 
                                    * (1./( (1./HetCP - 1./LCP )                     &                     
                                    + Anapartition*(1./HetCP - 1./RCP)))                            


        end if

        if(seletion==1) then
    
            Real_N_ROM = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)           &
                          *Conversion) *(1./( (1./HetCN - 1./LCN )*(1/AnaPartition)  &
                          +  (1./HetCN - 1./RCN)))

            Real_P_ROM = ((ImobilizationRateP*P)*Conversion)                         &
                          * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)              &
                          +  (1./HetCP - 1./RCP)))



            Real_N_LOM = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)           &
                         *Conversion) *(1./( (1./HetCN - 1./LCN )                    &
                         +  AnaPartition*(1./HetCN - 1./RCN)))


            Real_P_LOM = ((ImobilizationRateP*P)*Conversion)                         &
                         * (1./( (1./HetCP - 1./LCP ) +                              &
                         AnaPartition*(1./HetCP - 1./RCP)))


    
            if( Real_N_ROM<Real_P_ROM) then

                Me%Matrix(AMI, AMI)  = Me%Matrix(AMI, AMI) + DTDay                 & 
                                    * ImobilizationRateNH4  * (HetEF/HetCN)          &
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)    &
                                    +  (1./HetCN - 1./RCN))) 
                
                Me%Matrix(AMI, NII)  = Me%Matrix(AMI, NII) + DTDay                 &
                                    * ImobilizationRateNO3  * (HetEF/HetCN)          & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)    & 
                                    +  (1./HetCN - 1./RCN)))    


            else

                Me%Matrix(AMI, PI)  = Me%Matrix(AMI, PI) + DTDay                   &
                                    * ImobilizationRateP  * (HetEF/ HetCN)           & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)    &
                                    + (1./HetCP - 1./RCP)))

            end if

            if( Real_N_LOM<Real_P_LOM) then

                Me%Matrix(AMI, AMI)  = Me%Matrix(AMI, AMI) + DTDay                 & 
                                    * ImobilizationRateNH4  * (HetEF/HetCN)          &
                                    * (1./( (1./HetCN - 1./LCN )                     &
                                    +  AnaPartition*(1./HetCN - 1./RCN))) 
                
                Me%Matrix(AMI, NII)  = Me%Matrix(AMI, NII) + DTDay                 &
                                    * ImobilizationRateNO3  * (HetEF/HetCN)          & 
                                    * (1./( (1./HetCN - 1./LCN )                     & 
                                    +  AnaPartition*(1./HetCN - 1./RCN)))   


            else

                Me%Matrix(AMI, PI)  = Me%Matrix(AMI, PI) + DTDay                   &
                                    * ImobilizationRateP  * (HetEF/HetCN)            & 
                                    * (1./( (1./HetCP - 1./LCP )                     &
                                    + AnaPartition*(1./HetCP - 1./RCP)))

            end if
    
        end if


        !Sinks : NH4 imobilization


        if (seletion == 4 .OR. seletion ==5  .OR. seletion == 6) then


            Me%Matrix(AMI, LOM_CI) = Me%Matrix(AMI, LOM_CI) - DTDay             &
                                         * (LOM_C_Srate * (1./HetCN - 1./LCN))   &     
                                         *(1/(Partition+1.))*(1/Conversion)
            
            Me%Matrix(AMI, ROM_CI) = Me%Matrix(AMI, ROM_CI) - DTDay             &
                                         * (ROM_C_Srate * (1./HetCN - 1./RCN))   &
                                         *(1/(Partition+1.))*(1/Conversion) 

        
        end if

        if (seletion == 1 .OR. seletion == 2 .OR. seletion == 3) then

            !NH4 immobilization       
            Me%Matrix(AMI, AMI) = Me%Matrix(AMI, AMI) - DTDay                  &
                                      * ImobilizationRateNH4                     

        end if
   
        if (Me%PropCalc%Sol_Bacteria) then

            ! Source: haverá a excreção destas bacterias!

            !Sink: Excretion
            Me%Matrix(AMI, PFI) =   Me%Matrix(AMI, PFI) - DTDay                       &
                                  * solubilizingRate*solCEf*(1/LCN)       &
                                  * (0.1 /( (SolPartition + 1. )))  
                          
            Me%Matrix(AMI, PFI) = Me%Matrix(AMI, PFI) - DTDay                       &
                                  * SolubilizingRate *(1/RCN)              &
                                  *solCEf* (0.1 /( (1./SolPartition + 1. ) ))                                   

        end if


        ! Source: Urea Hydrolysis


        Me%Matrix(AMI,Iurea) = Me%Matrix(AMI,Iurea) + DTDay *ureahydrolysis


            ! Sink : volatilization

        wind    =   Me%ExternalVar%Wind(index)
        if (wind .gt. 0.0) then
            
            Temp    =   Me%ExternalVar%Temperature(index)
            AMI     =   Me%Propindex%Ammonia
            AM      =   Me%ExternalVar%Mass(AMI,index)        !User units here assumed mg/L
            H       =   Me%Hydrogen                           !Mol/L
            
            !Convert assumed ammonia mg/L in mol/L
            !Molecular weight NH4 - 18.0385 g/mol
            ![mol/L] = [mg/L] * 1g/1000mg / [g/mol]
            AM      = AM / (1000. * 18.0385)
            
            TK      =  29447
            PANH3   =  2.45E-8
            EK      =  8.79E-12
            XK1     =  1000
            XK      =  -0.25
            TF      = TK * EXP(-6/(1.99E-3 * (temp + 273.15)))

            PNH3    = EK * AM * 7.14286E-11 / H       
        
            XKG     = XK1 * log(Wind) * exp(XK)

            Kvol    = XKG * TF * (PNH3 - PANH3)



    !        Me%Matrix(AMI,AMI)= Me%Matrix(AMI,AMI) + DTDay * Kvol
        
        endif


        !Independent term
        Me%IndTerm(AMI) = Me%ExternalVar%Mass(AMI, index)
                          

    end subroutine SQAmmonia
    !----------------------------------------------------------------------------


    !NITRATE 
    !
    !SOURCES: - Nitrification 
    !SINKS:   - Immobilization, denitrification
    !----------------------------------------------------------------------------
    subroutine SQNitrate(index)

        !Arguments---------------------------------------------------------------
        integer                 , intent(IN) :: index

        !Local-------------------------------------------------------------------
        integer :: AMI, NII
        
        real    :: NitrificationRate, DenitrificationRate, ImobilizationRateNO3

        real    :: LOM_C_Srate, LCN
        integer :: LOM_CI

        real    :: ROM_C_Srate, RCN
        integer :: ROM_CI

        real    :: HetCN, AutoEf, AnaNEf

        real    :: Partition, Conversion, DTDay, seletion    !!!Lúcia
        real    :: NO3,NO3limit
        !------------------------------------------------------------------------

        NII     = Me%PropIndex%Nitrate
        AMI     = Me%PropIndex%Ammonia
        DTDay   = Me%DTDay

        NitrificationRate    = Me%SpecificRates%AmmoniaToNitrate%Value
        DenitrificationRate  = Me%SpecificRates%NitrateToNgas%Value
        ImobilizationRateNO3 = Me%SpecificRates%NitrateImobilization%Value
        
        AutoEf               = Me%Microorganisms%Autotrophs%EficiencyN
        Conversion           = Me%ExternalVar%DissolvedToParticulate (index)
        seletion             = Me%Select          !!!Lúcia
        
        Partition   = Me%Partition
                
        LOM_C_Srate = Me%SpecificRates%Labil_OM_C%Value
        LCN         = Me%LabiOM_CN_Ratio          
        LOM_CI      = Me%PropIndex%Labil_OM_C

        ROM_C_Srate = Me%SpecificRates%RefractOM_C%Value
        RCN         = Me%RefractOM_CN_Ratio   
        ROM_CI      = Me%PropIndex%RefractOM_C

        HetCN       = Me%Microorganisms%Heterotrophs%CNRatio
        
        NO3                 = Me%ExternalVar%Mass(NII,index)
        NO3limit            = Me%NO3limit
        
        AnaNEf              = Me%Microorganisms%Anaerobic%EficiencyN


        Me%Matrix(NII, NII) = 1.
        
            
        !Sources :Nitrification

        Me%Matrix(NII, AMI) = Me%Matrix(NII, AMI) + DTDay                           &
                            * NitrificationRate * (1. - AutoEf) 
        
        if (NO3> NO3limit) then  
            !Sinks:Denitrification
            
            !This equation is inconsistent with the one in SQNGas
            !Me%Matrix(NII, NII) = Me%Matrix(NII, NII) - DTDay * DenitrificationRate 
            !so it was changed to match (Anaerobic organisms work)
            Me%Matrix(NII, NII) = Me%Matrix(NII, NII) - DTDay * DenitrificationRate * (1. - AnaNEf)
            
        end If

        !Sink:Immobilization NO3

        if (seletion == 4 .OR. seletion == 5 .OR. seletion == 6) then


            Me%Matrix(NII, LOM_CI) = Me%Matrix(NII, LOM_CI) - DTDay*(1/conversion)    &
                                     * (LOM_C_Srate * (1./HetCN - 1./LCN))            &
                                     *(1/((1/partition)+1))
        
            Me%Matrix(NII, ROM_CI) = Me%Matrix(NII, ROM_CI) - DTDay *(1/conversion)   &    !Lúcia
                                     * (ROM_C_Srate * (1./HetCN - 1./RCN))            &
                                     *(1/((1/partition)+1))
        end if

        if (seletion == 1 .OR. seletion == 2 .OR. seletion == 3) then

              !NO3 immobilization
            Me%Matrix(NII, NII) = Me%Matrix(NII, NII) - DTDay                   &
                                   * ImobilizationRateNO3 

        end if


        !Independent term
        Me%IndTerm(NII) = Me%ExternalVar%Mass(NII, index) 


    end subroutine SQNitrate
    !----------------------------------------------------------------------------

  
  
    !Labil Organic Nitrogen 
    !
    !SOURCES: - Microorganisms death
    !SINKS:   - Heterotrphs uptake, including anaerobic uptake 
    !----------------------------------------------------------------------------
    subroutine SQLabilOrganicNitrogen (index)

        !Arguments---------------------------------------------------------------
        integer                , intent(IN)         :: index

        !Local-------------------------------------------------------------------
        integer :: AMI, NII,LOM_CI,LOM_NI,PI,PFI,ROM_CI             !!!Lúcia
        
        real    :: ImobilizationRateNH4, DenitrificationRate, ImobilizationRateNO3 !!!Lúcia
        real    :: solubilizingrate
        real    :: ImobilizationRateP 
        real    :: LOM_C_Srate, LCN,LCP             !!!Lúcia

        real    :: RCN
        real    :: RCP                              !!!Lúcia
        real    :: HeteroDeathRate, HetCN , HetCP   !!!Lúcia   
        integer :: HetI                
        
        real    :: AnaDeathRate, AnaCN       !!!Lúcia   
        integer :: AnaI                
        
        real    :: AutoDeathRate, AutoCN    !!!Lúcia
        integer :: AutI               

        real    :: Partition, AnaPartition, Conversion, DTDay, seletion
        real    :: solpartition                       !!!Lúcia

        real    :: AM,N,P,Real_N,Real_P         !!!Lúcia
        real    :: NO3,NO3limit,MethaneProduction

        !------------------------------------------------------------------------
        
        LOM_C_Srate         = Me%SpecificRates%Labil_OM_C%Value
        LCN                 = Me%LabiOM_CN_Ratio
        LCP                 = Me%LabilOM_CP_Ratio       !!!Lúcia 
        LOM_CI              = Me%PropIndex%Labil_OM_C
        ROM_CI              = Me%PropIndex%RefractOM_C

        LOM_NI              = Me%PropIndex%Labil_OM_N
            
        RCN                 = Me%RefractOM_CN_Ratio   
        RCP                 = Me%RefractOM_CP_Ratio         !!!Lúcia
               
        DenitrificationRate     = Me%SpecificRates%NitrateToNgas%Value
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value
        ImobilizationRateP      = Me%SpecificRates%PhosphorusImobilization%Value  !!!Lúcia
        solubilizingrate        = Me%SpecificRates%Solubilizing%Value           
                      
        NII                 = Me%PropIndex%Nitrate
        AMI                 = Me%PropIndex%Ammonia
        PI                  = Me%PropIndex%Inorganic_P_soluble   !!!Lúcia
        PFI                 = Me%PropIndex%Inorganic_P_fix  

        Partition           = Me%Partition
        AnaPartition        = Me%AnaerobicPartition
        seletion            = Me%Select                          !!!Lúcia
        solpartition        = Me%Solpartition

        HeteroDeathRate     = Me%SpecificRates%Heterotrophs%Value
        HetI                = Me%PropIndex%HeterotrophicC
        HetCN               = Me%Microorganisms%Heterotrophs%CNRatio
        HetCP               = Me%Microorganisms%Heterotrophs%CPRatio    !!!Lúcia

        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaI                = Me%PropIndex%AnaerobicC
        AnaCN               = Me%Microorganisms%Anaerobic%CNRatio

        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutI                = Me%PropIndex%AutotrophicC
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio

        DTDay               = Me%DTDay

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)

        AM   = Me%ExternalVar%Mass(AMI,index)   !!!
        N    = Me%ExternalVar%Mass(NII,index)   !!!
        P    = Me%ExternalVar%Mass(PI,index)    !!!
       
        methaneproduction      = Me%SpecificRates%MethaneProduction%Value
        NO3                    = Me%ExternalVar%Mass(NII,index)

        NO3limit               = Me%NO3limit



        Me%Matrix(LOM_NI, LOM_NI)   = 1.
        
      
       !Sinks: Heterotrophs uptake 
        

        if (seletion==5.OR.seletion==6 .OR. seletion==8 .OR. seletion==9 )  then 
        
            Me%Matrix(LOM_NI, LOM_CI)  = Me%Matrix(LOM_NI, LOM_CI) - DTDay          &
                                          * LOM_C_Srate * 1./ LCN  

        end if

        if (seletion==2.OR. seletion==3) then 


            Me%Matrix(LOM_NI, AMI)  = Me%Matrix(LOM_NI, AMI) - DTDay                 & 
                                    * ImobilizationRateNH4  * Conversion* (1./ LCN)  &
                                    * (1./( (1./HetCN - 1./LCN ) +AnaPartition       &
                                    * (1./HetCN - 1./RCN))) 
                
            Me%Matrix(LOM_NI, NII)  = Me%Matrix(LOM_NI, NII) - DTDay                 &
                                    * ImobilizationRateNO3  * Conversion*(1./ LCN)   & 
                                    * (1./( (1./HetCN - 1./LCN ) + AnaPartition      &
                                    * (1./HetCN - 1./RCN)))
         end if

        if (seletion==4.OR. seletion==7) then

            
            Me%Matrix(LOM_NI, PI)  = Me%Matrix(LOM_NI, PI) - DTDay                   &
                                    * ImobilizationRateP  * Conversion*( 1./ LCN )   & 
                                    * (1./( (1./HetCP - 1./LCP ) + AnaPartition      &
                                    * (1./HetCP - 1./RCP)))

        end if

        
        if(seletion==1) then
    
            Real_N = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)               &
                     *Conversion)*(1./( (1./HetCN - 1./LCN )                         &
                     + AnaPartition * (1./HetCN - 1./RCN)))

            Real_P = ((ImobilizationRateP*P)*Conversion)                             &
                     * (1./( (1./HetCP - 1./LCP )                                    &
                     + AnaPartition * (1./HetCP - 1./RCP)))

    
            if( Real_N<Real_P) then

                    Me%Matrix(LOM_NI, AMI)  = Me%Matrix(LOM_NI, AMI) - DTDay           & 
                                    * ImobilizationRateNH4  * Conversion* (1./ LCN ) &
                                    * (1./( (1./HetCN - 1./LCN ) + AnaPartition      &
                                    * (1./HetCN - 1./RCN))) 
                
                    Me%Matrix(LOM_NI, NII)  = Me%Matrix(LOM_NI, NII) - DTDay           &
                                    * ImobilizationRateNO3  * Conversion * (1./ LCN) & 
                                    * (1./( (1./HetCN - 1./LCN ) + AnaPartition      &
                                    * (1./HetCN - 1./RCN))) 
                            
                                else

                    Me%Matrix(LOM_NI, PI)  = Me%Matrix(LOM_NI, PI) - DTDay             &
                                    * ImobilizationRateP  * Conversion * (1./ LCN )  & 
                                    * (1./( (1./HetCP - 1./LCP ) + AnaPartition      &
                                    * (1./HetCP - 1./RCP)))

            end if
        end if


        !Sinks    : Anaerobic uptake, vai depender se respiram NO3 ou CO2 
        if (Me%ComputeImobilization) then

            if (NO3>NO3limit) then    ! respiram NO3


                Me%Matrix(LOM_NI, NII) = Me%Matrix(LOM_NI, NII) - DTDay                      & 
                               * DenitrificationRate * Conversion * (1./ LCN )       &
                               * ( 0.1/( AnaPartition +1.))  
                               
            else                   ! respiram CO2
                               
                Me%Matrix(LOM_NI, LOM_CI)     = Me%Matrix(LOM_NI, LOM_CI) - DTDay     & 
                                     *  methaneproduction *(1/LCN)  &
                                     * ( 0.1/( AnaPartition +1.) )

                Me%Matrix(LOM_NI, ROM_CI)     = Me%Matrix(LOM_NI, ROM_CI) - DTDay     & 
                                     *  methaneproduction *(1/LCN)   &
                                     * ( 0.1/( AnaPartition +1.) )
             
            end if              
        else
        
            if (NO3>NO3limit) then    ! respiram NO3


                Me%Matrix(LOM_NI, NII) = Me%Matrix(LOM_NI, NII) - DTDay                      & 
                               * DenitrificationRate * Conversion * (1./ LCN )       &
                               * ( 0.1)  
                               
            else                   ! respiram CO2
                               
                Me%Matrix(LOM_NI, LOM_CI)     = Me%Matrix(LOM_NI, LOM_CI) - DTDay     & 
                                     *  methaneproduction *(1/LCN)  &
                                     * ( 0.1 )

                Me%Matrix(LOM_NI, ROM_CI)     = Me%Matrix(LOM_NI, ROM_CI) - DTDay     & 
                                     *  methaneproduction *(1/LCN)   &
                                     * ( 0.1 )
             
            end if                      
        endif                       
                                                                                                                      
        !Sources: Anaerobic death 
        if (.NOT. Me%Microorganisms%Anaerobic%LogicalMinumumPOP)                     &        
            Me%Matrix(LOM_NI, AnaI) = Me%Matrix(LOM_NI, AnaI) + DTDay                &
                                    * AnaDeathRate *( 1. / AnaCN )
                                              
        !Sources: Heterotrophic death
        if ( .NOT. Me%Microorganisms%Heterotrophs%LogicalMinumumPOP)                 &        
            Me%Matrix(LOM_NI, HetI) = Me%Matrix(LOM_NI, HetI) + DTDay                &
                                    * HeteroDeathRate * (1. / HetCN)
                                              
        !Sources: Autotrophic death
        if ( .NOT. Me%Microorganisms%Autotrophs%LogicalMinumumPOP)                   &        
            Me%Matrix(LOM_NI, AutI) = Me%Matrix(LOM_NI, AutI) + DTDay                &
                                    * AutoDeathRate * (1. / AutoCN)                                       


        if (Me%PropCalc%Sol_Bacteria) then

            !Sinks: Labil N
            Me%Matrix(LOM_NI, PFI) = Me%Matrix(LOM_NI, PFI) - DTDay                     &
                                  * solubilizingRate * Conversion *(1/LCN)              &
                                  * (0.1 /( (SolPartition + 1. ))) 
        end if

       
        !Independent term
        Me%IndTerm(LOM_NI) = Me%ExternalVar%Mass(LOM_NI, index)  
        
                               
    end subroutine SQLabilOrganicNitrogen 
    !----------------------------------------------------------------------------



    !Refractary Organic Nitrogen
    !
    !SOURCES: - INTERPOOL TRANSFORMATION
    !SINKS:   - Heterotrphs uptake, including anaerobic uptake 
    !----------------------------------------------------------------------------
    subroutine SQRefractOrganicNitrogen (index)

        !Arguments---------------------------------------------------------------
        integer                , intent(IN) :: index

        !Local-------------------------------------------------------------------
        integer :: AMI, NII,PI,PFI              !!!Lucia
        
        real    :: ImobilizationRateNH4, DenitrificationRate, ImobilizationRateNO3  !!!Lúcia
        real    :: ImobilizationRateP
        real    :: solubilizingrate
        real    :: ROM_C_Srate, RCN , RCP    !!!Lúcia
        integer :: ROM_CI, ROM_NI,LOM_CI

        real    :: LCN
        real    :: LCP                          !!!Lúcia

        real    :: HeteroDeathRate, HetCN , HetCP   !!!Lúcia   
        integer :: HetI                
        
        real    :: AnaDeathRate, AnaCN     !!!Lúcia        
        integer :: AnaI                

        real    :: AutoDeathRate, AutoCN    !!!Lúcia
        integer :: AutoI               

        real    :: Partition, AnaPartition, Conversion, DTDay ,seletion
        real    :: solpartition !!!Lúcia
        real    :: AM,N,P,Real_N,Real_P         !!!Lúcia
        real    :: NO3,NO3Limit,MethaneProduction

        !------------------------------------------------------------------------
        ROM_C_Srate         = Me%SpecificRates%RefractOM_C%Value
        RCN                 = Me%RefractOM_CN_Ratio 
        RCP                 = Me%RefractOM_CP_Ratio         !!!Lúcia 
 
        ROM_CI              = Me%PropIndex%RefractOM_C
        LOM_CI              = Me%PropIndex%Labil_OM_C
        ROM_NI              = Me%PropIndex%RefractOM_N
        PFI                 = Me%PropIndex%Inorganic_P_fix  
           
        LCN                 = Me%LabiOM_CN_Ratio   
        LCP                 = Me%LabilOM_CP_Ratio       !!!Lúcia  
               
        DenitrificationRate     = Me%SpecificRates%NitrateToNgas%Value
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value
        ImobilizationRateP      = Me%SpecificRates%PhosphorusImobilization%Value  !!!Lúcia
        solubilizingrate        = Me%SpecificRates%Solubilizing%Value           
             
        NII                 = Me%PropIndex%Nitrate
        AMI                 = Me%PropIndex%Ammonia
        PI                  = Me%PropIndex%Inorganic_P_soluble   !!!Lúcia

        Partition           = Me%Partition
        AnaPartition        = Me%AnaerobicPartition
        solpartition        = Me%Solpartition

        HeteroDeathRate     = Me%SpecificRates%Heterotrophs%Value
        HetI                = Me%PropIndex%HeterotrophicC
        HetCN               = Me%Microorganisms%Heterotrophs%CNRatio
        HetCP               = Me%Microorganisms%Heterotrophs%CPRatio        !!!Lúcia

        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaI                = Me%PropIndex%AnaerobicC
        AnaCN               = Me%Microorganisms%Anaerobic%CNRatio

        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutoI               = Me%PropIndex%AutotrophicC
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio

        DTDay               = Me%DTDay
        seletion            = Me%Select                          !!!Lúcia       

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)
        
        AM   = Me%ExternalVar%Mass(AMI,index)   !!!
        N    = Me%ExternalVar%Mass(NII,index)   !!!
        P    = Me%ExternalVar%Mass(PI,index)    !!!

        NO3limit            = Me%NO3limit
        methaneproduction   = Me%SpecificRates%MethaneProduction%Value
        NO3                 = Me%ExternalVar%Mass(NII,index)




        Me%Matrix(ROM_NI, ROM_NI)  = 1.
        
        !Sinks Heterotrophs uptake           
       
        if (seletion==5.OR.seletion==6 .OR. seletion==8 .OR. seletion==9 )  then 
        
             Me%Matrix(ROM_NI, ROM_CI)  = Me%Matrix(ROM_NI, ROM_CI) - DTDay *         &
                                          ROM_C_Srate *(1./ RCN)  

        end if

        if (seletion==2.OR. seletion==3) then 


            Me%Matrix(ROM_NI, AMI)  = Me%Matrix(ROM_NI, AMI) - DTDay                  & 
                                    * ImobilizationRateNH4  * Conversion*(1./ RCN)    &
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)     &
                                    +  (1./HetCN - 1./RCN))) 
                
            Me%Matrix(ROM_NI, NII)  = Me%Matrix(ROM_NI, NII) - DTDay                  &
                                    * ImobilizationRateNO3  * Conversion *(1./ RCN)   & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)     &
                                    +  (1./HetCN - 1./RCN)))
         end if

        if (seletion==4.OR. seletion==7) then

            
            Me%Matrix(ROM_NI, PI)  = Me%Matrix(ROM_NI, PI) - DTDay                    &
                                    * ImobilizationRateP  * Conversion * (1./ RCN )   & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)     &
                                    + (1./HetCP - 1./RCP)))

        end if

        
        if(seletion==1) then
    
            Real_N = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)                &
                     *Conversion)*(1./( (1./HetCN - 1./LCN )*(1/AnaPartition)         &
                     +  (1./HetCN - 1./RCN)))

            Real_P = ((ImobilizationRateP*P)*Conversion)                              &
                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)                     &
                    +  (1./HetCP - 1./RCP)))

    
            if( Real_N<Real_P) then

                Me%Matrix(ROM_NI, AMI)  = Me%Matrix(ROM_NI, AMI) - DTDay            & 
                                    * ImobilizationRateNH4  * Conversion *(1./ RCN )  &
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)     &
                                    +  (1./HetCN - 1./RCN))) 
                
                Me%Matrix(ROM_NI, NII)  = Me%Matrix(ROM_NI, NII) - DTDay            &
                                    * ImobilizationRateNO3  * Conversion *(1./ RCN )  & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)     &
                                    +  (1./HetCN - 1./RCN)))    
                            
            else

                Me%Matrix(ROM_NI, PI)  = Me%Matrix(ROM_NI, PI) - DTDay              &
                                    * ImobilizationRateP  * Conversion *(1./ RCN)     & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)     &
                                    + (1./HetCP - 1./RCP)))

            end if
        end if


        !Sinks    : Anaerobic uptake, vai depender se respiram NO3 ou CO2 
        if (Me%ComputeImobilization) then
            if (NO3>NO3limit) then    ! respiram NO3

                Me%Matrix(ROM_NI, NII) = Me%Matrix(ROM_NI, NII) - DTDay                       &
                                   * DenitrificationRate * Conversion * (1./RCN)          &
                                   * (0.1/( (1./AnaPartition) +1.) )
                                   
            else                   ! respiram CO2

                Me%Matrix(ROM_NI, LOM_CI)     = Me%Matrix(ROM_NI, LOM_CI) - DTDay    & 
                                         *  methaneproduction *(1./RCN)          &
                                         * ( 0.1/( (1/AnaPartition) +1.) )

                Me%Matrix(ROM_NI, ROM_CI)     = Me%Matrix(ROM_NI, ROM_CI) - DTDay    & 
                                         *  methaneproduction *(1./RCN)          &
                                         * ( 0.1/( (1/AnaPartition) +1.) )
                 
            end if              
        else
            if (NO3>NO3limit) then    ! respiram NO3

                Me%Matrix(ROM_NI, NII) = Me%Matrix(ROM_NI, NII) - DTDay                       &
                                   * DenitrificationRate * Conversion * (1./RCN)          &
                                   * (0.1 )
                                   
            else                   ! respiram CO2

                Me%Matrix(ROM_NI, LOM_CI)     = Me%Matrix(ROM_NI, LOM_CI) - DTDay    & 
                                         *  methaneproduction *(1./RCN)          &
                                         * ( 0.1 )

                Me%Matrix(ROM_NI, ROM_CI)     = Me%Matrix(ROM_NI, ROM_CI) - DTDay    & 
                                         *  methaneproduction *(1./RCN)          &
                                         * ( 0.1 )
                 
            end if                      
        endif
                                    
        if (Me%PropCalc%Sol_Bacteria) then
    
            !Sinks: Refract N
            Me%Matrix(ROM_NI, PFI) = Me%Matrix(ROM_NI, PFI) - DTDay                       &
                                  * SolubilizingRate * Conversion *(1/RCN)              &
                                  * (0.1 /( ((1./SolPartition) + 1. ) )) 
        end if

                                                 
        !Independent term
        Me%IndTerm(ROM_NI)     = Me%ExternalVar%Mass(ROM_NI, index)                                     


    end subroutine SQRefractOrganicNitrogen 
    !----------------------------------------------------------------------------

    
    !Heterotrophic N
    !
    !SOURCES: - Organic matter N decay
    !SINKS:   - Heterotrophic Death, excretion
    !----------------------------------------------------------------------------
    subroutine SQHeterotrophicN (index)

        !Arguments---------------------------------------------------------------
        integer                 , intent(IN) :: index

        !Local-------------------------------------------------------------------
        real    :: HeteroDeathRate, HetCN, HetEF , HetCP   !!!Lúcia
        integer :: HetCI, HetNI           
        
        real    :: ROM_C_Srate, RCN, RCP     !!!Lúcia
        integer :: ROM_CI, ROM_NI              

        real    :: LOM_C_Srate, LCN, LCP    !!!Lúcia          
        integer :: LOM_CI, LOM_NI

        
        real    :: ImobilizationRateNO3, ImobilizationRateNH4 ,ImobilizationRateP !!!Lúcia

        integer :: NII, AMI ,PI       !!!Lúcia
        
        real    :: Real_N_ROM , Real_P_ROM,Real_N_LOM, Real_P_LOM, AM,N, P   !!!Lúcia

        real    :: AnaPartition,Partition, Conversion, DTDay ,seletion          !!!Lúcia
        !------------------------------------------------------------------------

        HeteroDeathRate         = Me%SpecificRates%Heterotrophs%Value
        HetCI                   = Me%PropIndex%HeterotrophicC
        HetNI                   = Me%PropIndex%HeterotrophicN
        HetCN                   = Me%Microorganisms%Heterotrophs%CNRatio
        HetCP                   = Me%Microorganisms%Heterotrophs%CPRatio !!!Lúcia
        HetEF                   = Me%Microorganisms%Heterotrophs%EficiencyC

        ROM_C_Srate             = Me%SpecificRates%RefractOM_C%Value
        RCN                     = Me%RefractOM_CN_Ratio 
        RCP                     = Me%RefractOM_CP_Ratio  !!!Lúcia         
        ROM_CI                  = Me%PropIndex%RefractOM_C
        ROM_NI                  = Me%PropIndex%RefractOM_N
            
        LOM_C_Srate             = Me%SpecificRates%Labil_OM_C%Value
        LCN                     = Me%LabiOM_CN_Ratio
        LCP                     = Me%LabilOM_CP_Ratio    !!!Lúcia
        LOM_CI                  = Me%PropIndex%Labil_OM_C
        LOM_NI                  = Me%PropIndex%Labil_OM_N
        
        NII                     = Me%PropIndex%Nitrate
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value
        N                       = Me%ExternalVar%Mass(NII, index)       !!!Lúcia


        AMI                     = Me%PropIndex%Ammonia
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value
        AM                      = Me%ExternalVar%Mass(AMI, index)    !!!Lúcia

        PI                      = Me%PropIndex%Inorganic_P_soluble
        ImobilizationRateP      = Me%SpecificRates%PhosphorusImobilization%Value
        P                       = Me%ExternalVar%Mass(PI, index)     !!!Lúcia


        Partition               = Me%Partition
        AnaPartition            = Me%AnaerobicPartition
        
        DTDay                   = Me%DTDay

        Conversion              = Me%ExternalVar%DissolvedToParticulate (index)
        seletion                = Me%Select    !!!Lúcia


        
        Me%Matrix(HetNI, HetNI) = 1.
        
        !Sink: Heterotrophic death
        if ( .not. Me%Microorganisms%Heterotrophs%LogicalMinumumPOP)                        & 
                    Me%Matrix(HetNI, HetCI) = Me%Matrix(HetNI, HetCI) -                     &
                                              DTDay * HeteroDeathRate *( 1. / HetCN)

        !Sink : Breath
        !Source : Labile and Refractory uptake 

        if (seletion==5.OR.seletion==6 .OR. seletion==8 .OR. seletion==9 )  then 
        
            Me%Matrix(HetNI, LOM_CI)  = Me%Matrix(HetNI, LOM_CI) +DTDay *          &
                                         LOM_C_Srate * (1/LCN)

            Me%Matrix(HetNI, ROM_CI)  = Me%Matrix(HetNI, ROM_CI) +DTDay *          &
                                         ROM_C_Srate * (1/RCN)
        
                            
            !Breath

            Me%Matrix(HetNI, LOM_CI)  = Me%Matrix(HetNI, LOM_CI) -DTDay *          &
                                     LOM_C_Srate* (HetEF/HetCN)

            Me%Matrix(HetNI, ROM_CI)  = Me%Matrix(HetNI, ROM_CI) -DTDay *          &
                                     ROM_C_Srate* (HetEF/HetCN)
                      
        end if

        if (seletion==2.OR. seletion==3) then 

            Me%Matrix(HetNI, AMI)  = Me%Matrix(HetNI, AMI) + DTDay                   & 
                                   * ImobilizationRateNH4  * Conversion *(1/LCN )   &
                                   * 1./((1./HetCN - 1./LCN)                        &
                                   + AnaPartition * (1./HetCN - 1./RCN))               
                       

            Me%Matrix(HetNI, AMI)  = Me%Matrix(HetNI, AMI) + DTDay                   & 
                                   * ImobilizationRateNH4  * Conversion *(1/RCN )   &
                                   * 1./((1./HetCN - 1./LCN ) * (1./AnaPartition)   &
                                   + (1./HetCN - 1./RCN))


            Me%Matrix(HetNI, NII)  = Me%Matrix(HetNI, NII) + DTDay                   &
                                   * ImobilizationRateNO3  * Conversion *(1/LCN)    &
                                   * 1./((1./HetCN - 1./LCN)                        &
                                   + AnaPartition * (1./HetCN - 1./RCN))               
                                                


            Me%Matrix(HetNI, NII)  = Me%Matrix(HetNI, NII) + DTDay                   &
                                   * ImobilizationRateNO3  * Conversion *(1/RCN )   &
                                   * 1./((1./HetCN - 1./LCN ) * (1./ AnaPartition ) &
                                   + (1./HetCN - 1./RCN))                    


            !Breath

            Me%Matrix(HetNI, AMI)  = Me%Matrix(HetNI, AMI) - DTDay                   & 
                                   * ImobilizationRateNH4  * Conversion             &
                                   *(HetEF/HetCN )                                  &
                                   * 1./((1./HetCN - 1./LCN)                        &
                                   + AnaPartition * (1./HetCN - 1./RCN))               
                       

            Me%Matrix(HetNI, AMI)  = Me%Matrix(HetNI, AMI) - DTDay                   & 
                                   * ImobilizationRateNH4  * Conversion             &
                                   *(HetEF/HetCN)                                   &
                                   * 1./((1./HetCN - 1./LCN ) * (1./AnaPartition)   &
                                   + (1./HetCN - 1./RCN))


            Me%Matrix(HetNI, NII)  = Me%Matrix(HetNI, NII) - DTDay                   &
                                   * ImobilizationRateNO3  * Conversion             &
                                   *(HetEF/HetCN )                                  &
                                   * 1./((1./HetCN - 1./LCN)                        &
                                   + AnaPartition * (1./HetCN - 1./RCN))               
                                                


            Me%Matrix(HetNI, NII)  = Me%Matrix(HetNI, NII) - DTDay                   &
                                   * ImobilizationRateNO3  * Conversion             &
                                   *(HetEF/HetCN )                                  &
                                   * 1./((1./HetCN - 1./LCN ) * (1./ AnaPartition)  &
                                   + (1./HetCN - 1./RCN))                    

        end if


        if (seletion==4.OR. seletion==7) then

            
            Me%Matrix(HetNI, PI)  = Me%Matrix(HetNI, PI) + DTDay                     &
                                    * ImobilizationRateP  * Conversion * (1/RCN )    & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)    & 
                                    + (1./HetCP - 1./RCP)))                         


            Me%Matrix(HetNI, PI)  = Me%Matrix(HetNI, PI) + DTDay                     &
                                    * ImobilizationRateP  * Conversion * (1/LCN )    & 
                                    * (1./( (1./HetCP - 1./LCP )                     &
                                    + AnaPartition*(1./HetCP - 1./RCP)))


            !Breath

            Me%Matrix(HetNI, PI)  = Me%Matrix(HetNI, PI) - DTDay                     &
                                    * ImobilizationRateP  * Conversion               &
                                    * (HetEF/HetCN )                                 & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/Anapartition)    & 
                                    + (1./HetCP - 1./RCP)))                         

            Me%Matrix(HetNI, PI)  = Me%Matrix(HetNI, PI) - DTDay                     &
                                    * ImobilizationRateP  * Conversion               &
                                    * (HetEF/HetCN )                                 & 
                                    * (1./( (1./HetCP - 1./LCP )                     & 
                                    + Anapartition*(1./HetCP - 1./RCP)))                            



        end if

        

        if(seletion==1) then
    
            Real_N_ROM = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)           &
                        *Conversion) *(1./( (1./HetCN - 1./LCN )*(1/AnaPartition)    &
                         +  (1./HetCN - 1./RCN)))

            Real_P_ROM = ((ImobilizationRateP*P)*Conversion)                         &
                        * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)                &
                        +  (1./HetCP - 1./RCP)))



            Real_N_LOM = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)           &
                         *Conversion) *(1./( (1./HetCN - 1./LCN )                    &
                         +  AnaPartition*(1./HetCN - 1./RCN)))


            Real_P_LOM = ((ImobilizationRateP*P)*Conversion)                         &
                         * (1./( (1./HetCP - 1./LCP ) +                              &
                         AnaPartition*(1./HetCP - 1./RCP)))


    
            if( Real_N_ROM<Real_P_ROM) then

                  Me%Matrix(HetNI, AMI)  = Me%Matrix(HetNI, AMI) + DTDay                 & 
                                           * ImobilizationRateNH4  *                     &
                                           Conversion *(1/RCN)                           &
                                           * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition) &
                                            +  (1./HetCN - 1./RCN))) 
                

                  Me%Matrix(HetNI, NII)  = Me%Matrix(HetNI, NII) + DTDay                 &
                                          * ImobilizationRateNO3  *                      &
                                          Conversion * (1/RCN)                           & 
                                          * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)  & 
                                          +  (1./HetCN - 1./RCN)))
                                                
                    
                !Breath 
                
                        
                  Me%Matrix(HetNI, AMI)  = Me%Matrix(HetNI, AMI) - DTDay                 & 
                                    * ImobilizationRateNH4  * Conversion*(HetEF/HetCN)   &
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)        &
                                    +  (1./HetCN - 1./RCN))) 
                
                  Me%Matrix(HetNI, NII)  = Me%Matrix(HetNI, NII) - DTDay                 &
                                    * ImobilizationRateNO3  * Conversion *(HetEF/HetCN)  & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)        & 
                                    +  (1./HetCN - 1./RCN)))    


            else


                 Me%Matrix(HetNI, PI)  = Me%Matrix(HetNI, PI) + DTDay                    &
                                    * ImobilizationRateP  * Conversion*(1/RCN)           & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)        &
                                    + (1./HetCP - 1./RCP)))


                  !Breath

                  Me%Matrix(HetNI, PI)  = Me%Matrix(HetNI, PI) - DTDay                   &
                                    * ImobilizationRateP  * Conversion*(HetEF/ HetCN)    & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)        &
                                    + (1./HetCP - 1./RCP)))

            end if



            if( Real_N_LOM<Real_P_LOM) then

                Me%Matrix(HetNI, AMI)  = Me%Matrix(HetNI, AMI) + DTDay                 & 
                                    * ImobilizationRateNH4  * Conversion *(1/LCN )       &
                                    * (1./( (1./HetCN - 1./LCN )                         &
                                    + AnaPartition* (1./HetCN - 1./RCN))) 
                
                Me%Matrix(HetNI, NII)  = Me%Matrix(HetNI, NII) + DTDay                 &
                                    * ImobilizationRateNO3  * Conversion * (1/LCN )      & 
                                    * (1./( (1./HetCN - 1./LCN )                         & 
                                    +  AnaPartition*(1./HetCN - 1./RCN)))   
                                            

                !Breath 
                        
                Me%Matrix(HetNI, AMI)  = Me%Matrix(HetNI, AMI) - DTDay                 & 
                                * ImobilizationRateNH4  * Conversion*(HetEF/HetCN)   &
                                * (1./( (1./HetCN - 1./LCN )                         &
                                +  AnaPartition*(1./HetCN - 1./RCN))) 

                Me%Matrix(HetNI, NII)  = Me%Matrix(HetNI, NII) - DTDay                 &
                                * ImobilizationRateNO3  * Conversion *(HetEF/HetCN)  & 
                                * (1./( (1./HetCN - 1./LCN )                         & 
                                +  AnaPartition*(1./HetCN - 1./RCN)))   


                            else

                Me%Matrix(HetNI, PI)  = Me%Matrix(HetNI, PI) + DTDay                  &
                                * ImobilizationRateP  * Conversion*(1/LCN)           & 
                                * (1./( (1./HetCP - 1./LCP )                         &
                                + AnaPartition*(1./HetCP - 1./RCP)))


                !Breath

                Me%Matrix(HetNI, PI)  = Me%Matrix(HetNI, PI) - DTDay                   &
                                * ImobilizationRateP  * Conversion*(HetEF/HetCN)     & 
                                * (1./( (1./HetCP - 1./LCP )                         &
                                + AnaPartition*(1./HetCP - 1./RCP)))

            end if
    
        end if



        !Source  : imobilization of NO3 and NH4

        if (seletion == 4 .OR. seletion == 5 .OR. seletion == 6) then
        
             !NO3 immobilization
            Me%Matrix(HetNI, LOM_CI) = Me%Matrix(HetNI, LOM_CI) + DTDay          &
                                         * (LOM_C_Srate * (1./HetCN - 1./LCN))     &
                                         *(1/((1/partition)+1))
            
            Me%Matrix(HetNI, ROM_CI) = Me%Matrix(HetNI, ROM_CI) + DTDay          &
                                         * (ROM_C_Srate * (1./HetCN - 1./RCN))     &
                                         *(1/((1/partition)+1))
              !NH4 immobilization
            Me%Matrix(HetNI, LOM_CI) = Me%Matrix(HetNI, LOM_CI) + DTDay         &
                                         * (LOM_C_Srate * (1./HetCN - 1./LCN))     &
                                         *(1/(Partition+1.))
            
            Me%Matrix(HetNI, ROM_CI) = Me%Matrix(HetNI, ROM_CI) + DTDay         &
                                         * (ROM_C_Srate * (1./HetCN - 1./RCN))     &
                                         *(1/(Partition+1.)) 

        end if


        if (seletion == 1 .OR. seletion == 2 .OR. seletion == 3) then

            !NH4 immobilization
            Me%Matrix(HetNI, AMI) = Me%Matrix(HetNI, AMI) + DTDay                &
                          * ImobilizationRateNH4 * Conversion 

            !NO3 immobilization
            Me%Matrix(HetNI, NII) = Me%Matrix(HetNI, NII) + DTDay                &
                          * ImobilizationRateNO3 * Conversion 

        end if
        
        !Independent term
        Me%IndTerm(HetNI) = Me%ExternalVar%Mass(HetNI,index) 
           

    end subroutine SQHeterotrophicN 
    !----------------------------------------------------------------------------

    
    !Anaerobic N
    !
    !SOURCES: - Denitrification eficiency
    !SINKS:   - Death
    !----------------------------------------------------------------------------    
    subroutine SQAnaerobicN (index)

        !Arguments---------------------------------------------------------------
        integer                , intent(IN) :: index

        !Local-------------------------------------------------------------------
        real    :: AnaCN, AnaNEf, AnaCEf, AnaDeathRate
        integer :: AnaNI, AnaCI
        real    :: DenitrificationRate 
        real    :: LCN, RCN
        real    :: AnaPartition, Conversion, DTDay       
        integer :: NII
        real    :: NO3,NO3limit,methaneproduction
        integer :: LOM_CI,ROM_CI
            
        !------------------------------------------------------------------------
        
        DTDay               = Me%DTDay

        AnaNI               = Me%PropIndex%AnaerobicN
        AnaCI               = Me%PropIndex%AnaerobicC
        NII                 = Me%PropIndex%Nitrate
        ROM_CI              = Me%PropIndex%RefractOM_C
        LOM_CI              = Me%PropIndex%Labil_OM_C

        
        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaCN               = Me%Microorganisms%Anaerobic%CNRatio
        
        DenitrificationRate  = Me%SpecificRates%NitrateToNgas%Value
        AnaNEf              = Me%Microorganisms%Anaerobic%EficiencyN
        AnaCEf              = Me%Microorganisms%Anaerobic%EficiencyC

        LCN                 = Me%LabiOM_CN_Ratio
        RCN                 = Me%RefractOM_CN_Ratio
    
        AnaPartition        = Me%AnaerobicPartition

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)

        NO3                 = Me%ExternalVar%Mass(NII,index)
        NO3limit            = Me%NO3limit
        methaneproduction   = Me%SpecificRates%MethaneProduction%Value


        Me%Matrix(AnaNI, AnaNI) = 1.
        !Sink: Death
        if ( .not. Me%Microorganisms%Anaerobic%LogicalMinumumPOP )                  &        
            Me%Matrix(AnaNI, AnaCI) = Me%Matrix(AnaNI, AnaCI) - DTDay               &
                                     * AnaDeathRate / AnaCN  
        

        if (NO3>NO3limit) then
       
            !Sources: DeNitrification (direct assimilation)
            Me%Matrix(AnaNI, NII) = Me%Matrix(AnaNI, NII) + DTDay                       &
                                  * DenitrificationRate * Conversion * AnaNEf
            if (Me%ComputeImobilization) then
                !Sources: Labil N
                Me%Matrix(AnaNI, NII) = Me%Matrix(AnaNI, NII) + DTDay                       &
                                      * DenitrificationRate * Conversion *(1/LCN)           &
                                      * (0.1 /( (AnaPartition + 1. ))) 
                !Sources: Refract N
                Me%Matrix(AnaNI, NII) = Me%Matrix(AnaNI, NII) + DTDay                       &
                                      * DenitrificationRate * Conversion *(1/RCN)           &
                                      * (0.1 /( (1./AnaPartition + 1. ) ))
            else
            
                !Sources: Labil N
                Me%Matrix(AnaNI, NII) = Me%Matrix(AnaNI, NII) + DTDay                       &
                                      * DenitrificationRate * Conversion *(1/LCN)           &
                                      * (0.1 ) 
                !Sources: Refract N
                Me%Matrix(AnaNI, NII) = Me%Matrix(AnaNI, NII) + DTDay                       &
                                      * DenitrificationRate * Conversion *(1/RCN)           &
                                      * (0.1 )            
            endif                                                   
            !Sink: Excretion
            Me%Matrix(AnaNI, NII) = Me%Matrix(AnaNI, NII) - DTDay                       &
                                  * DenitrificationRate * Conversion                    &
                                  * AnaCEf*0.1/ AnaCN                       
 
        else 

            !Source :  degradação de matéria orgânica e respiração de CO2
        
                !( vem da labile)    
            if (Me%ComputeImobilization) then                
                Me%Matrix(AnaNI, LOM_CI)     = Me%Matrix(AnaNI, LOM_CI) + DTDay    & 
                                         *  methaneproduction *(1/LCN)         &
                                         * ( 0.1/( AnaPartition +1.) )

                Me%Matrix(AnaNI, ROM_CI)     = Me%Matrix(AnaNI, ROM_CI) + DTDay    & 
                                         *  methaneproduction *(1/LCN)         &
                                         * ( 0.1/( AnaPartition +1.) )

                    ! (vem da refractory)

                Me%Matrix(AnaNI, LOM_CI)     = Me%Matrix(AnaNI, LOM_CI) + DTDay    & 
                                         *  methaneproduction *(1/RCN)         &
                                         * ( 0.1/( (1/AnaPartition) +1.) )

                Me%Matrix(AnaNI, ROM_CI)     = Me%Matrix(AnaNI, ROM_CI) + DTDay    & 
                                         *  methaneproduction *(1/RCN)         &
                                         * ( 0.1/( (1/AnaPartition) +1.) )
            else

                Me%Matrix(AnaNI, LOM_CI)     = Me%Matrix(AnaNI, LOM_CI) + DTDay    & 
                                         *  methaneproduction *(1/LCN)         &
                                         * ( 0.1 )

                Me%Matrix(AnaNI, ROM_CI)     = Me%Matrix(AnaNI, ROM_CI) + DTDay    & 
                                         *  methaneproduction *(1/LCN)         &
                                         * ( 0.1 )

                    ! (vem da refractory)

                Me%Matrix(AnaNI, LOM_CI)     = Me%Matrix(AnaNI, LOM_CI) + DTDay    & 
                                         *  methaneproduction *(1/RCN)         &
                                         * ( 0.1 )

                Me%Matrix(AnaNI, ROM_CI)     = Me%Matrix(AnaNI, ROM_CI) + DTDay    & 
                                         *  methaneproduction *(1/RCN)         &
                                         * ( 0.1 )
            
            endif
            ! sink: excreção de NH4+

            Me%Matrix(AnaNI, LOM_CI)     = Me%Matrix(AnaNI, LOM_CI) - DTDay    & 
                                        *  methaneproduction      &
                                        *0.1*AnaCEf/AnaCN

            Me%Matrix(AnaNI, ROM_CI)     = Me%Matrix(AnaNI, ROM_CI) - DTDay    & 
                                        *  methaneproduction     &
                                        *0.1*AnaCEf/AnaCN
    

        end if 
 

       
        !Independent term
        Me%IndTerm(AnaNI)    = Me%ExternalVar%Mass(AnaNI, index) 


    end subroutine SQAnaerobicN 
    !----------------------------------------------------------------------------

    
    !Autotrophic N
    !
    !SOURCES: - Nitrification eficiency 
    !SINKS:   - Death
    !----------------------------------------------------------------------------    
    subroutine SQAutotrophicN (index)

        !Arguments---------------------------------------------------------------
        integer                 , intent(IN)    :: index

        !Local-------------------------------------------------------------------
        real    :: DTDay
        
        real    :: AutoCN, AutoEf, AutoDeathRate
        integer :: AutoNI, AutoCI

        real    :: NitificationRate, Conversion 

        integer :: AMI    
        !------------------------------------------------------------------------
        
        DTDay               = Me%DTDay

        AutoNI              = Me%PropIndex%AutotrophicN
        AutoCI              = Me%PropIndex%AutotrophicC
        AMI                 = Me%PropIndex%Ammonia
        
        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio
        
        NitificationRate    = Me%SpecificRates%AmmoniaToNitrate%Value
        AutoEf              = Me%Microorganisms%Autotrophs%EficiencyN
        
        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)

        Me%Matrix(AutoNI, AutoNI) = 1.

        !Sink: Death
        if ( .not. Me%Microorganisms%Autotrophs%LogicalMinumumPOP)                  &        
            Me%Matrix(AutoNI, AutoCI) = Me%Matrix(AutoNI, AutoCI) - DTDay           &
                                        * AutoDeathRate/AutoCN
        
        !Sources: Nitrification
        Me%Matrix(AutoNI, AMI) = Me%Matrix(AutoNI, AMI) + DTDay                     &
                                 * NitificationRate * Conversion * AutoEf/AutoCN
        
        !Independent term
        Me%IndTerm(AutoNI) = Me%ExternalVar%Mass(AutoNI, index) 


    end subroutine SQAutotrophicN 
    !----------------------------------------------------------------------------


      
    !N gas (N20 N2) - same units as dissolved properties
    !
    !SOURCES: - Denitrification
    !SINKS:   - none for the moment
    !----------------------------------------------------------------------------
    subroutine SQNgas (index)

        !Arguments---------------------------------------------------------------
        integer                , intent(IN) :: index
    
        !Local-------------------------------------------------------------------
        real    :: DTDay
        
        integer :: NGI, NII
        
        real    :: DenitrificationRate, AnaNEf
        real    :: NO3,NO3limit
        real    :: conversion
        
        !------------------------------------------------------------------------
        
        DTDay       = Me%DTDay

        NGI         = Me%PropIndex%Ngas
        NII         = Me%PropIndex%Nitrate

        DenitrificationRate     = Me%SpecificRates%NitrateToNgas%Value
        AnaNEf                  = Me%Microorganisms%Anaerobic%EficiencyN
        NO3                     = Me%ExternalVar%Mass(NII,index)
        Conversion              = Me%ExternalVar%DissolvedToParticulate (index)
        
  
        NO3limit            = Me%NO3limit
  
  
        Me%Matrix(NGI, NGI) =  1.

        if (NO3> NO3limit) then 
        
            !Sources: Denitrification

            Me%Matrix(NGI, NII) = Me%Matrix(NGI, NII) + DTDay * DenitrificationRate * (1-AnaNEf)

        end if 
        
        !Independent term
        Me%IndTerm(NGI)     = Me%ExternalVar%Mass(NGI, index) 


    end subroutine SQNgas 

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    subroutine SQUrea (index)

        !Arguments---------------------------------------------------------------
        integer                , intent(IN) :: index
    
        !Local-------------------------------------------------------------------
        real    :: DTDay
        
        integer :: Iurea
        
        real    :: ureahydrolysis
        
        !------------------------------------------------------------------------
  
        DTDay           = Me%DTDay
        Iurea           = Me%Propindex%urea
        ureahydrolysis  = Me%SpecificRates%UreaHydrolysis%Value

        
        Me%Matrix(Iurea, Iurea) =  1.


        ! sink: Urea Hydrolysis 

        Me%Matrix(Iurea,Iurea) = Me%Matrix(Iurea,Iurea) - DTDay *ureahydrolysis

        !Independent term
        Me%IndTerm(Iurea)     = Me%ExternalVar%Mass(Iurea, index) 


        end subroutine SQUrea
    !----------------------------------------------------------------------------

    subroutine SQAmmoniaGas (index)

        !Arguments---------------------------------------------------------------
        integer                , intent(IN) :: index
    
        !Local-------------------------------------------------------------------
        real    :: DTDay
        integer :: IAmmoniaGas,AMI
        real    :: wind, XKG,TF,EK,PNH3,PANH3,TK,XK1,XK
        real    :: Kvol,temp,AM,H      
        !------------------------------------------------------------------------

        DTDay                = Me%DTDay
        IAmmoniaGas           = Me%Propindex%AmmoniaGas

        wind    =   Me%ExternalVar%Wind(index)
    
        if (wind .gt. 0.0) then
            
            Temp    =   Me%ExternalVar%Temperature(index)
            AMI     =   Me%Propindex%Ammonia
            AM      =   Me%ExternalVar%Mass(AMI,index)        !User units here assumed mg/L
            H       =   Me%Hydrogen                           !Mol/L
            
            !Convert assumed ammonia mg/L in mol/L
            !Molecular weight NH4 - 18.0385 g/mol
            ![mol/L] = [mg/L] * 1g/1000mg / [g/mol]
            AM      = AM / (1000. * 18.0385)
            
            TK      =  29447
            PANH3   =  2.45E-8
            EK      =  8.79E-12
            XK1     =  1000
            XK      =  -0.25
            TF      = TK * EXP(-6/(1.99E-3 * (temp + 273.15)))

            PNH3    = EK * AM * 7.14286E-11 / H       
        
            XKG     = XK1 * log(Wind) * exp(XK)

            Kvol    = XKG * TF * (PNH3 - PANH3)

            ! source: volatilization
            Me%Matrix(IAmmoniaGas,AMI)= Me%Matrix(IAmmoniaGas,AMI)    &
                                            + DTDay * Kvol
        
        endif
        
        Me%Matrix(IAmmoniaGas,IAmmoniaGas)= 1


        !Independent term
        Me%IndTerm(IAmmoniaGas)    = Me%ExternalVar%Mass(IAmmoniaGas, index) 

    end subroutine SQAmmoniaGas


    !----------------------------------------------------------------------------
    subroutine SQPhosphorus(index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN) :: index
        !------------------------------------------------------------------------

        call SQLabilOrganicPhosphorus     (index)
  
        call SQRefractOrganicPhosphorus   (index)

        call SQHeterotrophicP              (index)

        call SQAutotrophicP               (index)

        call SQAnaerobicP                 (index)
        
        call SQInorganicPhosphorusSoluble (index)

        call SQInorganicPhosphorusFix     (index)        
    !------------------------------------------------------------------------

    end subroutine SQPhosphorus
    !---------------------------------------------------------------------------




    !Labil Organic Phosphorus 
    !
    !SOURCES: - Microorganisms death
    !SINKS:   - Heterotrphs uptake, including anaerobic uptake 
    !----------------------------------------------------------------------------
    subroutine SQLabilOrganicPhosphorus (index)

        !Arguments---------------------------------------------------------------
        integer                , intent(IN)         :: index


        !Local-------------------------------------------------------------------
        integer :: AMI, NII,LOM_CI,LOM_PI,PI,PFI,ROM_CI             !!!Lúcia
        
        real    :: ImobilizationRateNH4, DenitrificationRate, ImobilizationRateNO3 !!!Lúcia
        real    :: ImobilizationRateP 
        real    :: LOM_C_Srate, LCN,LCP
        real    :: solubilizingrate             !!!Lúcia

        real    :: RCN
        real    :: RCP                              !!!Lúcia
        real    :: HeteroDeathRate, HetCN , HetCP   !!!Lúcia   
        integer :: HetI                
        
        real    :: AnaDeathRate, AnaCP      !!!Lúcia   
        integer :: AnaI                
        
        real    :: AutoDeathRate, AutoCP   !!!Lúcia
        integer :: AutI               

        real    :: Partition, AnaPartition, Conversion, DTDay, seletion
        real    :: solpartition    !!!Lúcia

        real    :: AM,N,P,Real_N,Real_P         !!!Lúcia
        real    :: NO3,NO3limit, methaneProduction

        !------------------------------------------------------------------------
        
        LOM_C_Srate         = Me%SpecificRates%Labil_OM_C%Value
        LCN                 = Me%LabiOM_CN_Ratio
        LCP                 = Me%LabilOM_CP_Ratio       !!!Lúcia 
        LOM_CI              = Me%PropIndex%Labil_OM_C
        ROM_CI              = Me%PropIndex%RefractOM_C
        LOM_PI              = Me%PropIndex%Labil_OM_P
        PFI                 = Me%PropIndex%Inorganic_P_fix
 
        RCN                 = Me%RefractOM_CN_Ratio   
        RCP                 = Me%RefractOM_CP_Ratio         !!!Lúcia
               

        DenitrificationRate     = Me%SpecificRates%NitrateToNgas%Value
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value
        ImobilizationRateP      = Me%SpecificRates%PhosphorusImobilization%Value  !!!Lúcia
        solubilizingrate        = Me%SpecificRates%Solubilizing%Value   
                      
        NII                 = Me%PropIndex%Nitrate
        AMI                 = Me%PropIndex%Ammonia
        PI                  = Me%PropIndex%Inorganic_P_soluble   !!!Lúcia

        Partition           = Me%Partition
        AnaPartition        = Me%AnaerobicPartition
        seletion            = Me%Select                          !!!Lúcia
        solpartition        = Me%Solpartition

        HeteroDeathRate     = Me%SpecificRates%Heterotrophs%Value
        HetI                = Me%PropIndex%HeterotrophicC
        HetCN               = Me%Microorganisms%Heterotrophs%CNRatio
        HetCP               = Me%Microorganisms%Heterotrophs%CPRatio    !!!Lúcia


        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaI                = Me%PropIndex%AnaerobicC
        AnaCP               = Me%Microorganisms%Anaerobic%CPRatio

        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutI                = Me%PropIndex%AutotrophicC
        AutoCP              = Me%Microorganisms%Autotrophs%CPRatio



        DTDay               = Me%DTDay

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)

        AM   = Me%ExternalVar%Mass(AMI,index)   !!!
        N    = Me%ExternalVar%Mass(NII,index)   !!!
        P    = Me%ExternalVar%Mass(PI,index)    !!!


        methaneproduction      = Me%SpecificRates%MethaneProduction%Value
        NO3                    = Me%ExternalVar%Mass(NII,index)

        NO3limit               = Me%NO3limit



        Me%Matrix(LOM_PI, LOM_PI)   = 1.
        
           
       !Sinks Heterotrophs uptake 
        

        if (seletion==5.OR.seletion==6 .OR. seletion==8 .OR. seletion==9 )  then 
        
             Me%Matrix(LOM_PI, LOM_CI)  = Me%Matrix(LOM_PI, LOM_CI) -               &
                                          DTDay * LOM_C_Srate * 1./ LCP  


        end if


        if (seletion==2.OR. seletion==3) then 


            Me%Matrix(LOM_PI, AMI)  = Me%Matrix(LOM_PI, AMI) - DTDay                & 
                                    * ImobilizationRateNH4  * Conversion* 1./ LCP   &
                                    * (1./( (1./HetCN - 1./LCN ) + AnaPartition     &
                                    * (1./HetCN - 1./RCN))) 
                
            Me%Matrix(LOM_PI, NII)  = Me%Matrix(LOM_PI, NII) - DTDay                &
                                    * ImobilizationRateNO3  * Conversion* 1./ LCP   & 
                                    * (1./( (1./HetCN - 1./LCN ) + AnaPartition     &
                                    * (1./HetCN - 1./RCN)))
         end if


        if (seletion==4.OR. seletion==7) then

            
            Me%Matrix(LOM_PI, PI)  = Me%Matrix(LOM_PI, PI) - DTDay                  &
                                    * ImobilizationRateP  * Conversion*(1./ LCP )   & 
                                    * (1./( (1./HetCP - 1./LCP ) + AnaPartition     &
                                    * (1./HetCP - 1./RCP)))

        end if


        if(seletion==1) then
    
            Real_N = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)              &
                      *Conversion)*(1./( (1./HetCN - 1./LCN )                       &
                      + AnaPartition * (1./HetCN - 1./RCN)))

            Real_P = ((ImobilizationRateP*P)*Conversion)                            &
                    * (1./( (1./HetCP - 1./LCP ) +                                  &
                    AnaPartition * (1./HetCP - 1./RCP)))

    
            if( Real_N<Real_P) then

                  Me%Matrix(LOM_PI, AMI)  = Me%Matrix(LOM_PI, AMI) - DTDay          & 
                                    * ImobilizationRateNH4  * Conversion* 1./ LCP   &
                                    * (1./( (1./HetCN - 1./LCN ) + AnaPartition     &
                                    * (1./HetCN - 1./RCN))) 
                
                  Me%Matrix(LOM_PI, NII)  = Me%Matrix(LOM_PI, NII) - DTDay          &
                                    * ImobilizationRateNO3  * Conversion * 1./ LCP  & 
                                    * (1./( (1./HetCN - 1./LCN ) + AnaPartition     &
                                    * (1./HetCN - 1./RCN))) 
                            
            else

                  Me%Matrix(LOM_PI, PI)  = Me%Matrix(LOM_PI, PI) - DTDay             &
                                    * ImobilizationRateP  * Conversion * 1./ LCP     & 
                                    * (1./( (1./HetCP - 1./LCP ) + AnaPartition      &
                                    * (1./HetCP - 1./RCP)))

            end if
        end if


        !Sinks    : Anaerobic uptake, vai depender se respiram NO3 ou CO2 
        if (Me%ComputeImobilization) then
        
            if (NO3>NO3limit) then    ! respiram NO3


                Me%Matrix(LOM_PI, NII) = Me%Matrix(LOM_PI, NII) - DTDay                     & 
                                       * DenitrificationRate * Conversion * (1./ LCP)         &
                                       * ( 0.1/( AnaPartition +1.))

            else                   ! respiram CO2

                                   
                Me%Matrix(LOM_PI, LOM_CI)     = Me%Matrix(LOM_PI, LOM_CI) - DTDay     & 
                                         *  methaneproduction *(1/LCP)  &
                                         * ( 0.1/( AnaPartition +1.) )

                Me%Matrix(LOM_PI, ROM_CI)     = Me%Matrix(LOM_PI, ROM_CI) - DTDay     & 
                                         *  methaneproduction *(1/LCP)   &
                                         * ( 0.1/( AnaPartition +1.) )
                 
            end if              
        else

            if (NO3>NO3limit) then    ! respiram NO3


                Me%Matrix(LOM_PI, NII) = Me%Matrix(LOM_PI, NII) - DTDay                     & 
                                       * DenitrificationRate * Conversion * (1./ LCP)         &
                                       * ( 0.1)

            else                   ! respiram CO2

                                   
                Me%Matrix(LOM_PI, LOM_CI)     = Me%Matrix(LOM_PI, LOM_CI) - DTDay     & 
                                         *  methaneproduction *(1/LCP)  &
                                         * ( 0.1 )

                Me%Matrix(LOM_PI, ROM_CI)     = Me%Matrix(LOM_PI, ROM_CI) - DTDay     & 
                                         *  methaneproduction *(1/LCP)   &
                                         * ( 0.1 )
                 
            end if              
        
        endif
                                    
        !Sources: Anaerobic death
        if (.NOT. Me%Microorganisms%Anaerobic%LogicalMinumumPOP)                    &        
            Me%Matrix(LOM_PI, AnaI) = Me%Matrix(LOM_PI, AnaI) + DTDay               &
                                    * AnaDeathRate * 1. / AnaCP 


        !Sources: Heterotrophic death
        if ( .NOT. Me%Microorganisms%Heterotrophs%LogicalMinumumPOP)                &        
            Me%Matrix(LOM_PI, HetI) = Me%Matrix(LOM_PI, HetI) + DTDay               &
                                    * HeteroDeathRate * 1. / HetCP
                                              
        !Sources: Autotrophic death
        if ( .NOT. Me%Microorganisms%Autotrophs%LogicalMinumumPOP)                  &        
            Me%Matrix(LOM_PI, AutI) = Me%Matrix(LOM_PI, AutI) + DTDay               &
                                    * AutoDeathRate * 1. / AutoCP   
                                    
                                                                        
        if (Me%PropCalc%Sol_Bacteria) then

            !Sinks: Labil P
            Me%Matrix(LOM_PI, PFI) = Me%Matrix(LOM_PI, PFI) - DTDay                       &
                                  * solubilizingRate * Conversion *(1/LCP)              &
                                  * (0.1 /( (SolPartition + 1. ))) 
        end if


        !Independent term
        Me%IndTerm(LOM_PI) = Me%ExternalVar%Mass(LOM_PI, index)                     
    !----------------------------------------------------------------------------

    end subroutine SQLabilOrganicPhosphorus 
    !----------------------------------------------------------------------------


    !Refractary Organic Phosphorus
    !
    !SOURCES: - INTERPOOL TRANSFORMATION
    !SINKS:   - Heterotrphs uptake, including anaerobic uptake 
    !----------------------------------------------------------------------------
    subroutine SQRefractOrganicPhosphorus (index)

        !Arguments---------------------------------------------------------------
        integer                , intent(IN) :: index
       !Local-------------------------------------------------------------------
        integer :: AMI, NII,PI,PFI              !!!Lucia
        
        real    :: ImobilizationRateNH4, DenitrificationRate, ImobilizationRateNO3 !!!Lúcia
        real    :: ImobilizationRateP 
        real    :: solubilizingrate
        real    :: ROM_C_Srate, RCN , RCP    !!!Lúcia
        integer :: ROM_CI, ROM_PI,LOM_CI

        real    :: LCN
        real    :: LCP                          !!!Lúcia

        real    :: HeteroDeathRate, HetCN , HetCP   !!!Lúcia   
        integer :: HetI                
        
        real    :: AnaDeathRate, AnaCN     !!!Lúcia        
        integer :: AnaI                

        real    :: AutoDeathRate, AutoCN   !!!Lúcia
        integer :: AutoI               

        real    :: Partition, AnaPartition, Conversion, DTDay ,seletion,solpartition !!!Lúcia
        real    :: AM,N,P,Real_N,Real_P         !!!Lúcia
        real    :: NO3,NO3limit,methaneproduction

        !------------------------------------------------------------------------

        !------------------------------------------------------------------------
        ROM_C_Srate         = Me%SpecificRates%RefractOM_C%Value
        RCN                 = Me%RefractOM_CN_Ratio 
        RCP                 = Me%RefractOM_CP_Ratio         !!!Lúcia 
 
        ROM_CI              = Me%PropIndex%RefractOM_C
        LOM_CI              = Me%PropIndex%Labil_OM_C
        ROM_PI              = Me%PropIndex%RefractOM_P
        PFI                 = Me%PropIndex%Inorganic_P_fix



        LCN                 = Me%LabiOM_CN_Ratio   
        LCP                 = Me%LabilOM_CP_Ratio       !!!Lúcia  
               
        DenitrificationRate     = Me%SpecificRates%NitrateToNgas%Value
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value
        ImobilizationRateP      = Me%SpecificRates%PhosphorusImobilization%Value  !!!Lúcia
        solubilizingrate        = Me%SpecificRates%Solubilizing%Value   
             
        NII                 = Me%PropIndex%Nitrate
        AMI                 = Me%PropIndex%Ammonia
        PI                  = Me%PropIndex%Inorganic_P_soluble   !!!Lúcia

        Partition           = Me%Partition
        AnaPartition        = Me%AnaerobicPartition
        solpartition        = Me%Solpartition

        HeteroDeathRate     = Me%SpecificRates%Heterotrophs%Value
        HetI                = Me%PropIndex%HeterotrophicC
        HetCN               = Me%Microorganisms%Heterotrophs%CNRatio
        HetCP               = Me%Microorganisms%Heterotrophs%CPRatio        !!!Lúcia

        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaI                = Me%PropIndex%AnaerobicC
        AnaCN               = Me%Microorganisms%Anaerobic%CNRatio

        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutoI               = Me%PropIndex%AutotrophicC
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio

        DTDay               = Me%DTDay
        seletion            = Me%Select                          !!!Lúcia       

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)
        
        AM   = Me%ExternalVar%Mass(AMI,index)   !!!
        N    = Me%ExternalVar%Mass(NII,index)   !!!
        P    = Me%ExternalVar%Mass(PI,index)    !!!


        NO3limit            = Me%NO3limit
        methaneproduction   = Me%SpecificRates%MethaneProduction%Value
        NO3                 = Me%ExternalVar%Mass(NII,index)





        Me%Matrix(ROM_PI, ROM_PI)  = 1.


        
        !Sinks Heterotrophs uptake           
       
        if (seletion==5.OR.seletion==6 .OR. seletion==8 .OR. seletion==9 )  then 
        
            Me%Matrix(ROM_PI, ROM_CI)  = Me%Matrix(ROM_PI, ROM_CI) -               &
                                          DTDay * ROM_C_Srate *(1./ RCP)  

        end if

        if (seletion==2.OR. seletion==3) then 


            Me%Matrix(ROM_PI, AMI)  = Me%Matrix(ROM_PI, AMI) - DTDay                & 
                                    * ImobilizationRateNH4  * Conversion*(1./ RCP)  &
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)   &
                                    +  (1./HetCN - 1./RCN))) 
                
            Me%Matrix(ROM_PI, NII)  = Me%Matrix(ROM_PI, NII) - DTDay                &
                                    * ImobilizationRateNO3  * Conversion *(1./ RCP) & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)   &
                                    +  (1./HetCN - 1./RCN)))
         end if

        if (seletion==4.OR. seletion==7) then

            
            Me%Matrix(ROM_PI, PI)  = Me%Matrix(ROM_PI, PI) - DTDay                  &
                                    * ImobilizationRateP  * Conversion *(1./ RCP)   & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)   &
                                    + (1./HetCP - 1./RCP)))

        end if

        
        if(seletion==1) then
    
            Real_N = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)              &
                      *Conversion)*(1./( (1./HetCN - 1./LCN )*(1/AnaPartition)      &
                      +  (1./HetCN - 1./RCN)))

            Real_P = ((ImobilizationRateP*P)*Conversion)                            &
                      * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)                 &
                      +  (1./HetCP - 1./RCP)))

    
            if( Real_N<Real_P) then

                Me%Matrix(ROM_PI, AMI)  = Me%Matrix(ROM_PI, AMI) - DTDay          & 
                                    * ImobilizationRateNH4  * Conversion *1./ RCP   &
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)   &
                                    +  (1./HetCN - 1./RCN))) 
                
                Me%Matrix(ROM_PI, NII)  = Me%Matrix(ROM_PI, NII) - DTDay          &
                                    * ImobilizationRateNO3  * Conversion *1./ RCP   & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)   &
                                    +  (1./HetCN - 1./RCN)))    
                            
            else

                Me%Matrix(ROM_PI, PI)  = Me%Matrix(ROM_PI, PI) - DTDay            &
                                    * ImobilizationRateP  * Conversion *1./ RCP     & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)   &
                                    + (1./HetCP - 1./RCP)))

            end if
        end if


        !Sinks    : Anaerobic uptake, vai depender se respiram NO3 ou CO2 
        if (Me%ComputeImobilization) then
            if (NO3>NO3limit) then    ! respiram NO3

                Me%Matrix(ROM_PI, NII) = Me%Matrix(ROM_PI, NII) - DTDay                     &
                                   * DenitrificationRate * Conversion * (1./RCP)          &
                                   * (0.1/( 1./AnaPartition +1.) )     


             else                   ! respiram CO2

                Me%Matrix(ROM_PI, LOM_CI)     = Me%Matrix(ROM_PI, LOM_CI) - DTDay    & 
                                         *  methaneproduction *(1./RCP)          &
                                         * ( 0.1/( (1/AnaPartition) +1.) )

                Me%Matrix(ROM_PI, ROM_CI)     = Me%Matrix(ROM_PI, ROM_CI) - DTDay    & 
                                         *  methaneproduction *(1./RCP)          &
                                         * ( 0.1/( (1/AnaPartition) +1.) )
                 
            end if              
        else

            if (NO3>NO3limit) then    ! respiram NO3

                Me%Matrix(ROM_PI, NII) = Me%Matrix(ROM_PI, NII) - DTDay                     &
                                   * DenitrificationRate * Conversion * (1./RCP)          &
                                   * (0.1 )     


             else                   ! respiram CO2

                Me%Matrix(ROM_PI, LOM_CI)     = Me%Matrix(ROM_PI, LOM_CI) - DTDay    & 
                                         *  methaneproduction *(1./RCP)          &
                                         * ( 0.1 )

                Me%Matrix(ROM_PI, ROM_CI)     = Me%Matrix(ROM_PI, ROM_CI) - DTDay    & 
                                         *  methaneproduction *(1./RCP)          &
                                         * ( 0.1 )
                 
            end if              
        
        endif


        if (Me%PropCalc%Sol_Bacteria) then

            !Sinks: Refract P

            Me%Matrix(ROM_PI, PFI) = Me%Matrix(ROM_PI, PFI) - DTDay                       &
                                  * SolubilizingRate * Conversion *(1/RCP)              &
                                  * (0.1 /( (1./SolPartition + 1. ) )) 

        end if
                                             
                                  
        !Independent term
        Me%IndTerm(ROM_PI)     = Me%ExternalVar%Mass(ROM_PI, index)                                     
    !----------------------------------------------------------------------------

    end subroutine SQRefractOrganicPhosphorus 
    !----------------------------------------------------------------------------



    
    !Heterotrophic P
    !
    !SOURCES: - Organic matter P decay
    !SINKS:   - Heterotrophic Death, excretion
    !----------------------------------------------------------------------------

    subroutine SQHeterotrophicP (index)

        !Arguments---------------------------------------------------------------
        integer                 , intent(IN) :: index

        !Local-------------------------------------------------------------------
        real    :: HeteroDeathRate, HetCN, HetEF , HetCP   !!!Lúcia
        integer :: HetCI, HetPI           
        
        real    :: ROM_C_Srate, RCN, RCP     !!!Lúcia
        integer :: ROM_CI, ROM_PI              

        real    :: LOM_C_Srate, LCN, LCP    !!!Lúcia          
        integer :: LOM_CI, LOM_PI

        
        real    :: ImobilizationRateNO3, ImobilizationRateNH4 ,ImobilizationRateP !!!Lúcia

        integer :: NII, AMI ,PI       !!!Lúcia
        
        real    :: Real_N_ROM , Real_P_ROM,Real_N_LOM, Real_P_LOM, AM,N, P   !!!Lúcia

        real    :: Partition,AnaPartition, Conversion, DTDay ,seletion          !!!Lúcia

        !------------------------------------------------------------------------


        HeteroDeathRate         = Me%SpecificRates%Heterotrophs%Value
        HetCI                   = Me%PropIndex%HeterotrophicC
        HetPI                   = Me%PropIndex%HeterotrophicP
        HetCN                   = Me%Microorganisms%Heterotrophs%CNRatio
        HetCP                   = Me%Microorganisms%Heterotrophs%CPRatio !!!Lúcia
        HetEF                   = Me%Microorganisms%Heterotrophs%EficiencyC


        ROM_C_Srate             = Me%SpecificRates%RefractOM_C%Value
        RCN                     = Me%RefractOM_CN_Ratio 
        RCP                     = Me%RefractOM_CP_Ratio  !!!Lúcia         
        ROM_CI                  = Me%PropIndex%RefractOM_C
        ROM_PI                  = Me%PropIndex%RefractOM_P

        LOM_C_Srate             = Me%SpecificRates%Labil_OM_C%Value
        LCN                     = Me%LabiOM_CN_Ratio
        LCP                     = Me%LabilOM_CP_Ratio    !!!Lúcia
        LOM_CI                  = Me%PropIndex%Labil_OM_C
        LOM_PI                  = Me%PropIndex%Labil_OM_P

        
        NII                     = Me%PropIndex%Nitrate
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value
        N                       = Me%ExternalVar%Mass(NII, index)       !!!Lúcia



        AMI                     = Me%PropIndex%Ammonia
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value
        AM                      = Me%ExternalVar%Mass(AMI, index)    !!!Lúcia

        PI                      = Me%PropIndex%Inorganic_P_soluble
        ImobilizationRateP      = Me%SpecificRates%PhosphorusImobilization%Value
        P                       = Me%ExternalVar%Mass(PI, index)     !!!Lúcia


        Partition               = Me%Partition
        AnaPartition            = Me%AnaerobicPartition
        
        DTDay                   = Me%DTDay

        Conversion              = Me%ExternalVar%DissolvedToParticulate (index)
        seletion                = Me%Select    !!!Lúcia


        
        Me%Matrix(HetPI, HetPI) = 1.
        

        !Sink: Heterotrophic death
        if ( .not. Me%Microorganisms%Heterotrophs%LogicalMinumumPOP)                & 
            Me%Matrix(HetPI, HetCI) = Me%Matrix(HetPI, HetCI) - DTDay               &
                                     * HeteroDeathRate * (1. / HetCP)

        !Sink : Breath
        !Source : Labile and Refractory uptake 

        if (seletion==5.OR.seletion==6 .OR. seletion==8 .OR. seletion==9 )  then 
         
            Me%Matrix(HetPI, LOM_CI)  = Me%Matrix(HetPI, LOM_CI) +DTDay *          &
                                         LOM_C_Srate * (1/LCP)

            Me%Matrix(HetPI, ROM_CI)  = Me%Matrix(HetPI, ROM_CI) +DTDay *          &
                                         ROM_C_Srate * (1/RCP)
                                
            !Breath

            Me%Matrix(HetPI, LOM_CI)  = Me%Matrix(HetPI, LOM_CI) -DTDay            &
                                        * LOM_C_Srate* (HetEF/HetCP)

            Me%Matrix(HetPI, ROM_CI)  = Me%Matrix(HetPI, ROM_CI) -DTDay            &
                                        * ROM_C_Srate* (HetEF/HetCP)
                      
        end if

        if (seletion==2.OR. seletion==3) then 

            Me%Matrix(HetPI, AMI)  = Me%Matrix(HetPI, AMI) + DTDay                   & 
                                   * ImobilizationRateNH4  * Conversion *(1/LCP)    &
                                   * 1./((1./HetCN - 1./LCN)                        &
                                   + AnaPartition * (1./HetCN - 1./RCN))               
                                   

            Me%Matrix(HetPI, AMI)  = Me%Matrix(HetPI, AMI) + DTDay                   & 
                                   * ImobilizationRateNH4  * Conversion *(1/RCP)    &
                                   * 1./((1./HetCN - 1./LCN ) * (1./AnaPartition)   &
                                   + (1./HetCN - 1./RCN))


            Me%Matrix(HetPI, NII)  = Me%Matrix(HetPI, NII) + DTDay                   &
                                   * ImobilizationRateNO3  * Conversion *(1/LCP)    &
                                   * 1./((1./HetCN - 1./LCN)                        &
                                   + AnaPartition * (1./HetCN - 1./RCN))               
                                                            


            Me%Matrix(HetPI, NII)  = Me%Matrix(HetPI, NII) + DTDay                   &
                                   * ImobilizationRateNO3  * Conversion *(1/RCP)    &
                                   * 1./((1./HetCN - 1./LCN ) *(1./ AnaPartition)   &
                                   + (1./HetCN - 1./RCN))                    


            !Breath

            Me%Matrix(HetPI, AMI)  = Me%Matrix(HetPI, AMI) - DTDay                   & 
                                   * ImobilizationRateNH4  *                        &
                                   Conversion *(HetEF/HetCP)                        &
                                   * 1./((1./HetCN - 1./LCN)                        &
                                   + AnaPartition * (1./HetCN - 1./RCN))               
                                   

            Me%Matrix(HetPI, AMI)  = Me%Matrix(HetPI, AMI) - DTDay                   & 
                                   * ImobilizationRateNH4  *                        &
                                   Conversion *(HetEF/HetCP)                        &
                                   * 1./((1./HetCN - 1./LCN ) * (1./AnaPartition)   &
                                   + (1./HetCN - 1./RCN))


            Me%Matrix(HetPI, NII)  = Me%Matrix(HetPI, NII) - DTDay                   &
                                   * ImobilizationRateNO3  *                        &
                                   Conversion *(HetEF/HetCP)                        &
                                   * 1./((1./HetCN - 1./LCN)                        &
                                   + AnaPartition * (1./HetCN - 1./RCN))               
                                                            


            Me%Matrix(HetPI, NII)  = Me%Matrix(HetPI, NII) - DTDay                   &
                                   * ImobilizationRateNO3  * Conversion             &
                                   *(HetEF/HetCP)                                   &
                                   * 1./((1./HetCN - 1./LCN ) *(1./ AnaPartition)   &
                                   + (1./HetCN - 1./RCN))                    
        
        end if


        if (seletion==4.OR. seletion==7) then

            
            Me%Matrix(HetPI, PI)  = Me%Matrix(HetPI, PI) + DTDay                     &
                                    * ImobilizationRateP  * Conversion * (1/RCP)     & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)    & 
                                    + (1./HetCP - 1./RCP)))                         


            Me%Matrix(HetPI, PI)  = Me%Matrix(HetPI, PI) + DTDay                     &
                                    * ImobilizationRateP  * Conversion * (1/LCP)     & 
                                    * (1./( (1./HetCP - 1./LCP )                     &
                                    + AnaPartition*(1./HetCP - 1./RCP)))


            !Breath

            Me%Matrix(HetPI, PI)  = Me%Matrix(HetPI, PI) - DTDay                     &
                                    * ImobilizationRateP  * Conversion               &
                                    * (HetEF/HetCP)                                  & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/Anapartition)    & 
                                    + (1./HetCP - 1./RCP)))                         


            Me%Matrix(HetPI, PI)  = Me%Matrix(HetPI, PI) - DTDay                     &
                                    * ImobilizationRateP  * Conversion               &
                                    * (HetEF/HetCP)                                  & 
                                    * (1./( (1./HetCP - 1./LCP )                     & 
                                    + Anapartition*(1./HetCP - 1./RCP)))                            


        end if

        

        if(seletion==1) then
    
            Real_N_ROM = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)           &
                         *Conversion) *(1./( (1./HetCN - 1./LCN )*(1/AnaPartition)   &
                         +  (1./HetCN - 1./RCN)))

            Real_P_ROM = ((ImobilizationRateP*P)*Conversion)                         &          
                        * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)                & 
                        +  (1./HetCP - 1./RCP)))



            Real_N_LOM = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)           &
                         *Conversion) *(1./( (1./HetCN - 1./LCN )                    &
                         +  AnaPartition*(1./HetCN - 1./RCN)))


            Real_P_LOM = ((ImobilizationRateP*P)*Conversion)                         &
                    * (1./( (1./HetCP - 1./LCP )                                     &
                    +  AnaPartition*(1./HetCP - 1./RCP)))


    
            if( Real_N_ROM<Real_P_ROM) then

                Me%Matrix(HetPI, AMI)  = Me%Matrix(HetPI,AMI) + DTDay              & 
                                           * ImobilizationRateNH4                    &
                                           * Conversion *1/RCP                       &
                                           * (1./( (1./HetCN - 1./LCN )              &
                                           *(1/AnaPartition)                         &
                                            +  (1./HetCN - 1./RCN))) 
                

                Me%Matrix(HetPI, NII)  = Me%Matrix(HetPI, NII) + DTDay             &
                                    * ImobilizationRateNO3  * Conversion * 1/RCP     & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)    & 
                                    +  (1./HetCN - 1./RCN)))
                                    
                                
                !Breath 
                
                        
                Me%Matrix(HetPI, AMI)  = Me%Matrix(HetPI, AMI) - DTDay              & 
                                    * ImobilizationRateNH4  * Conversion              &
                                    *(HetEF/HetCP)                                    &
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)     &
                                    +  (1./HetCN - 1./RCN))) 
                
                Me%Matrix(HetPI, NII)  = Me%Matrix(HetPI, NII) - DTDay              &
                                    * ImobilizationRateNO3  * Conversion              &
                                    *(HetEF/HetCP)                                    & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)     & 
                                    +  (1./HetCN - 1./RCN)))    


            else


                Me%Matrix(HetPI, PI)  = Me%Matrix(HetPI, PI) + DTDay                 &
                                    * ImobilizationRateP  * Conversion*(1/RCP)        & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)     &
                                    + (1./HetCP - 1./RCP)))


                  !Breath

                Me%Matrix(HetPI, PI)  = Me%Matrix(HetPI, PI) - DTDay                 &
                                    * ImobilizationRateP  * Conversion*(HetEF/ HetCP)  & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)      &
                                    + (1./HetCP - 1./RCP)))

            end if



            if( Real_N_LOM<Real_P_LOM) then

                Me%Matrix(HetPI, AMI)  = Me%Matrix(HetPI, AMI) + DTDay               & 
                                    * ImobilizationRateNH4  * Conversion * 1/LCP       &
                                    * (1./( (1./HetCN - 1./LCN )                       &
                                    + AnaPartition* (1./HetCN - 1./RCN))) 
                
                Me%Matrix(HetPI, NII)  = Me%Matrix(HetPI, NII) + DTDay               &
                                    * ImobilizationRateNO3  * Conversion * 1/LCP       & 
                                    * (1./( (1./HetCN - 1./LCN )                       & 
                                    +  AnaPartition*(1./HetCN - 1./RCN)))                                               

                !Breath 
                
                        
                Me%Matrix(HetPI, AMI)  = Me%Matrix(HetPI, AMI) - DTDay               & 
                                    * ImobilizationRateNH4  * Conversion               &
                                    *(HetEF/HetCP)                                     &
                                    * (1./( (1./HetCN - 1./LCN )                       &
                                    +  AnaPartition*(1./HetCN - 1./RCN))) 
                
                Me%Matrix(HetPI, NII)  = Me%Matrix(HetPI, NII) - DTDay               &
                                    * ImobilizationRateNO3  * Conversion               &
                                    *(HetEF/HetCP)                                     & 
                                    * (1./( (1./HetCN - 1./LCN )                       & 
                                    +  AnaPartition*(1./HetCN - 1./RCN)))   


            else

                Me%Matrix(HetPI, PI)  = Me%Matrix(HetPI, PI) + DTDay                &
                                    * ImobilizationRateP  * Conversion*(1/LCP)         & 
                                    * (1./( (1./HetCP - 1./LCP )                       &
                                    + AnaPartition*(1./HetCP - 1./RCP)))


                !Breath

                Me%Matrix(HetPI, PI)  = Me%Matrix(HetPI, PI) - DTDay                 &
                                    * ImobilizationRateP  * Conversion*(HetEF/HetCP)   & 
                                    * (1./( (1./HetCP - 1./LCP )                       &
                                    + AnaPartition*(1./HetCP - 1./RCP)))

            end if
    
        end if


            !Source  : imobilization of Soluble Phosphorus



        if (seletion == 2 .OR. seletion == 5 .OR. seletion == 8) then
        
             ! RealP imobilization

            Me%Matrix(HetPI, LOM_CI) = Me%Matrix(HetPI, LOM_CI) + DTDay         &
                                         * (LOM_C_Srate * (1./HetCP - 1./LCP))
            
            Me%Matrix(HetPI, ROM_CI) = Me%Matrix(HetPI, ROM_CI) + DTDay         &
                                         * (ROM_C_Srate * (1./HetCP - 1./RCP))
        
        end if

             ! Potential P immobilization

        if (seletion == 1 .OR. seletion == 4 .OR. seletion == 7) then


            Me%Matrix(HetPI, PI) = Me%Matrix(HetPI, PI) + DTDay                &            
                              * ImobilizationRateP * Conversion 


        end if


        !Independent term
        Me%IndTerm(HetPI) = Me%ExternalVar%Mass(HetPI,index) 
            


    !----------------------------------------------------------------------------

    end subroutine SQHeterotrophicP 
    !----------------------------------------------------------------------------


    !Anaerobic P
    !
    !SOURCES: - Denitrification eficiency - soluble phosphorus
    !SINKS:   - Death
    !----------------------------------------------------------------------------    
    subroutine SQAnaerobicP (index)

        !Arguments---------------------------------------------------------------
        integer                , intent(IN) :: index


        !Local-------------------------------------------------------------------
        real    :: AnaCP, AnaNEf, AnaCEf, AnaDeathRate
        integer :: AnaPI, AnaCI,LOM_CI,ROM_CI

        real    :: DenitrificationRate 

        real    :: LCP, RCP

        real    :: AnaPartition, Conversion, DTDay
        
        integer :: NII    
        real    :: NO3,NO3limit,methaneproduction
        !------------------------------------------------------------------------

        LOM_CI              = Me%PropIndex%Labil_OM_C
        ROM_CI              = Me%PropIndex%RefractOM_C
        DTDay               = Me%DTDay

        AnaPI               = Me%PropIndex%AnaerobicP
        AnaCI               = Me%PropIndex%AnaerobicC
        NII                 = Me%PropIndex%Nitrate
        
        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaCP               = Me%Microorganisms%Anaerobic%CPRatio
        
        DenitrificationRate  = Me%SpecificRates%NitrateToNgas%Value
        AnaNEf              = Me%Microorganisms%Anaerobic%EficiencyN
        AnaCEf              = Me%Microorganisms%Anaerobic%EficiencyC

        LCP                 = Me%LabilOM_CP_Ratio
        RCP                 = Me%RefractOM_CP_Ratio
    
        AnaPartition        = Me%AnaerobicPartition

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)


        NO3                 = Me%ExternalVar%Mass(NII,index)
        NO3limit            = Me%NO3limit
        methaneproduction   = Me%SpecificRates%MethaneProduction%Value



        Me%Matrix(AnaPI, AnaPI) = 1.


        !Sink: Death
        if ( .not. Me%Microorganisms%Anaerobic%LogicalMinumumPOP )                    &        
            Me%Matrix(AnaPI, AnaCI) = Me%Matrix(AnaPI, AnaCI) - DTDay * AnaDeathRate / AnaCP  


        if (NO3>NO3limit) then

            if (Me%ComputeImobilization) then

                !Sources: Labil P
                Me%Matrix(AnaPI, NII) = Me%Matrix(AnaPI, NII) + DTDay                         &
                                      * DenitrificationRate * Conversion *(1/LCP)             &
                                      * (0.1 /( (AnaPartition + 1. ))) 

               !Sources: Refract P
                Me%Matrix(AnaPI, NII) = Me%Matrix(AnaPI, NII) + DTDay                         &
                                      * DenitrificationRate * Conversion *(1/RCP)             &
                                      * (0.1 /( (1./AnaPartition + 1. ))) 
            
            else

                !Sources: Labil P
                Me%Matrix(AnaPI, NII) = Me%Matrix(AnaPI, NII) + DTDay                         &
                                      * DenitrificationRate * Conversion *(1/LCP)             &
                                      * (0.1 ) 

               !Sources: Refract P
                Me%Matrix(AnaPI, NII) = Me%Matrix(AnaPI, NII) + DTDay                         &
                                      * DenitrificationRate * Conversion *(1/RCP)             &
                                      * (0.1 ) 
            
            endif
            
            !Sink: Excretion
            Me%Matrix(AnaPI, NII) = Me%Matrix(AnaPI, NII) - DTDay                         &
                                  * DenitrificationRate * Conversion                      &
                                  *AnaCEf*0.1/ AnaCP       

        else 

            !Source :  degradação de matéria orgânica e respiração de CO2
        
                !( vem da labile)    
            if (Me%ComputeImobilization) then    
                Me%Matrix(AnaPI, LOM_CI)     = Me%Matrix(AnaPI, LOM_CI) + DTDay    & 
                                         *  methaneproduction *(1/LCP)         &
                                         * ( 0.1/( AnaPartition +1.) )

                Me%Matrix(AnaPI, ROM_CI)     = Me%Matrix(AnaPI, ROM_CI) + DTDay    & 
                                         *  methaneproduction *(1/LCP)         &
                                         * ( 0.1/( AnaPartition +1.) )
        
                    ! (vem da refractory)

                Me%Matrix(AnaPI, LOM_CI)     = Me%Matrix(AnaPI, LOM_CI) + DTDay    & 
                                         *  methaneproduction *(1/RCP)         &
                                         * ( 0.1/( (1/AnaPartition) +1.) )

                Me%Matrix(AnaPI, ROM_CI)     = Me%Matrix(AnaPI, ROM_CI) + DTDay    & 
                                         *  methaneproduction *(1/RCP)         &
                                         * ( 0.1/( (1/AnaPartition) +1.) )
            else

                Me%Matrix(AnaPI, LOM_CI)     = Me%Matrix(AnaPI, LOM_CI) + DTDay    & 
                                         *  methaneproduction *(1/LCP)         &
                                         * ( 0.1 )

                Me%Matrix(AnaPI, ROM_CI)     = Me%Matrix(AnaPI, ROM_CI) + DTDay    & 
                                         *  methaneproduction *(1/LCP)         &
                                         * ( 0.1 )
        
                    ! (vem da refractory)

                Me%Matrix(AnaPI, LOM_CI)     = Me%Matrix(AnaPI, LOM_CI) + DTDay    & 
                                         *  methaneproduction *(1/RCP)         &
                                         * ( 0.1 )

                Me%Matrix(AnaPI, ROM_CI)     = Me%Matrix(AnaPI, ROM_CI) + DTDay    & 
                                         *  methaneproduction *(1/RCP)         &
                                         * ( 0.1 )
            
            endif
    
            ! sink: excreção de fósforo soluvel

            Me%Matrix(AnaPI, LOM_CI)     = Me%Matrix(AnaPI, LOM_CI) - DTDay    & 
                                        *  methaneproduction      &
                                        *0.1*AnaCEf/AnaCP

            Me%Matrix(AnaPI, ROM_CI)     = Me%Matrix(AnaPI, ROM_CI) - DTDay    & 
                                        *  methaneproduction      &
                                        *0.1*AnaCEf/AnaCP
    
        end if 



       !Independent term
        Me%IndTerm(AnaPI)    = Me%ExternalVar%Mass(AnaPI, index) 
    !----------------------------------------------------------------------------

    end subroutine SQAnaerobicP 
    !----------------------------------------------------------------------------


    !Autotrophic P
    !
    !SOURCES: - uptake phosphorus
    !SINKS:   - Death
    !----------------------------------------------------------------------------    
    subroutine SQAutotrophicP (index)

        !Arguments---------------------------------------------------------------
        integer                 , intent(IN)    :: index


        !Local-------------------------------------------------------------------
        real    :: DTDay
        
        real    :: AutoCP, AutoEf, AutoDeathRate,AutoCN
        integer :: AutoPI, AutoCI

        real    :: NitificationRate, Conversion 

        integer :: AMI,PI
        
        real    :: AutoNP     !!!nova variavel  
        !------------------------------------------------------------------------

        

        DTDay               = Me%DTDay
        AutoPI              = Me%PropIndex%AutotrophicP
        AutoCI              = Me%PropIndex%AutotrophicC
        AMI                 = Me%PropIndex%Ammonia
        PI                  = Me%PropIndex%Inorganic_P_soluble
        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutoCP              = Me%Microorganisms%Autotrophs%CPRatio
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio
        AutoNP              = AutoCP/AutoCN             !!!!! atenção a esta formula!!!

        NitificationRate    = Me%SpecificRates%AmmoniaToNitrate%Value
        AutoEf              = Me%Microorganisms%Autotrophs%EficiencyN
        
        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)

        Me%Matrix(AutoPI, AutoPI) = 1.

        !Sink: Death
        if ( .not. Me%Microorganisms%Autotrophs%LogicalMinumumPOP)                  &        
            Me%Matrix(AutoPI, AutoCI) = Me%Matrix(AutoPI, AutoCI) - DTDay *         &
                                        AutoDeathRate / AutoCP  


        !Sources: uptake phosphorus
        Me%Matrix(AutoPI, AMI) = Me%Matrix(AutoPI, AMI) + DTDay *                   &
                                 NitificationRate * Conversion * AutoEf/AutoCP


        
        !Independent term
        Me%IndTerm(AutoPI) = Me%ExternalVar%Mass(AutoPI, index) 
    !----------------------------------------------------------------------------


    end subroutine SQAutotrophicP 



    subroutine SQInorganicPhosphorusSoluble (index)


       !Arguments---------------------------------------------------------------
        integer                 , intent(IN)    :: index


        !Local-------------------------------------------------------------------

        real    :: DTDay
        integer :: PI, NII, AMI, LOM_CI, ROM_CI, PFI
        real    :: AM, N, P ,PF
        real    :: DenitrificationRate, NitrificationRate
        real    :: ImobilizationRateNH4
        real    :: ImobilizationRateP,ROM_C_Srate,LOM_C_Srate,ImobilizationRateNO3
        real    :: AnaNEf,AnaCEf, AutoEf 
        real    :: AnaCP, AutoNP, AutoCN,AutoCP
        real    :: LCP,RCP, LCN, RCN
        real    :: HetCN,HetCP, HetEf
        real    :: Real_N_LOM,Real_N_ROM,Real_P_LOM,Real_P_ROM
        real    :: pai
        real    :: SolCEf,SolCP
        real    :: solubilizingrate 
        real    :: conversion, anaPartition, seletion,partition,solpartition
        real    :: NO3,NO3limit,methaneproduction
        real    :: SoilDens
        !------------------------------------------------------------------------

        DTDay               = Me%DTDay
        
        AMI                 = Me%PropIndex%Ammonia
        PI                  = Me%PropIndex%Inorganic_P_soluble
        NII                 = Me%PropIndex%Nitrate
        PFI                 = Me%PropIndex%Inorganic_P_fix

        AM                      = Me%ExternalVar%Mass(AMI, index)       !!!Lúcia
        N                       = Me%ExternalVar%Mass(NII, index)       !!!Lúcia
        P                       = Me%ExternalVar%Mass(PI, index)       !!!Lúcia
        PF                      = Me%ExternalVar%Mass(PFI, index)


        LOM_C_Srate             = Me%SpecificRates%Labil_OM_C%Value
        LCN                     = Me%LabiOM_CN_Ratio
        LCP                     = Me%LabilOM_CP_Ratio    !!!Lúcia
        LOM_CI                  = Me%PropIndex%Labil_OM_C


        ROM_C_Srate             = Me%SpecificRates%RefractOM_C%Value
        RCN                     = Me%RefractOM_CN_Ratio 
        RCP                     = Me%RefractOM_CP_Ratio  !!!Lúcia         
        ROM_CI                  = Me%PropIndex%RefractOM_C

        HetCP                   = Me%Microorganisms%Heterotrophs%CPRatio
        HetCN                   = Me%Microorganisms%Heterotrophs%CNRatio
        HetEF                   = Me%Microorganisms%Heterotrophs%EficiencyC

        AnaNEf              = Me%Microorganisms%Anaerobic%EficiencyN
        AnaCEf              = Me%Microorganisms%Anaerobic%EficiencyC
        AnaCP               = Me%Microorganisms%Anaerobic%CPRatio

        AutoEf              = Me%Microorganisms%Autotrophs%EficiencyN
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio
        AutoCP              = Me%Microorganisms%Autotrophs%CPRatio
        AutoNP              = AutoCP/AutoCN             !!!!! atenção a esta formula!!!
        solCEf              = Me%Microorganisms%Heterotrophs%EficiencyC
        SolCP               = Me%Microorganisms%Sols%CPRatio

        DenitrificationRate     = Me%SpecificRates%NitrateToNgas%Value
        ImobilizationRateP      = Me%SpecificRates%PhosphorusImobilization%Value
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value
        NitrificationRate       = Me%SpecificRates%AmmoniaToNitrate%Value
        solubilizingrate        = Me%SpecificRates%Solubilizing%Value   

        AnaPartition        = Me%AnaerobicPartition
        solpartition        = Me%Solpartition

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)

        seletion            = Me%Select    

        Partition               = Me%Partition

        SoilDens             = Me%ExternalVar%SoilDryDensity(index)

        NO3                 = Me%ExternalVar%Mass(NII,index)
        NO3limit            = Me%NO3limit
        methaneproduction   = Me%SpecificRates%MethaneProduction%Value



        PAI= Me%ExternalVar%PAI(index)

        Me%Matrix(PI, PI) = 1.   

        !Source: Excretion of Soluble Phosphorus vt anaerobics

        if (NO3>NO3limit) then

 
            Me%Matrix(PI, NII) =  Me%Matrix(PI, NII) + DTDay                       &
                              * DenitrificationRate                          &
                              * AnaCEf*0.1/ AnaCP
        else 

              Me%Matrix(PI, LOM_CI)     = Me%Matrix(PI, LOM_CI) + DTDay    & 
                                            *  (methaneproduction/conversion)      &
                                            *0.1*AnaCEf/AnaCP

              Me%Matrix(PI, ROM_CI)     = Me%Matrix(PI, ROM_CI) + DTDay    & 
                                            *  (methaneproduction/conversion)      &
                                            *0.1*AnaCEf/AnaCP
        end if

        
        !Sink: uptake phosphorus autotrophics

        Me%Matrix(PI, AMI) = Me%Matrix(PI, AMI) - DTDay * NitrificationRate *AutoEf/AutoCP


        !Source: Heterotrophic excretion 

        if (seletion==5.OR.seletion==6 .OR. seletion==8 .OR. seletion==9 )  then 


            Me%Matrix(PI, LOM_CI)  = Me%Matrix(PI, LOM_CI) +DTDay *                &
                                 LOM_C_Srate* (HetEF/HetCP)*(1/conversion)

            Me%Matrix(PI, ROM_CI)  = Me%Matrix(PI, ROM_CI) +DTDay *                &
                                 ROM_C_Srate* (HetEF/HetCP)*(1/conversion)
          
        end if


        if (seletion==2.OR. seletion==3) then 

            Me%Matrix(PI, AMI)  = Me%Matrix(PI, AMI) + DTDay                        & 
                                   * ImobilizationRateNH4 *(HetEF/HetCP)           &
                                   * (1./((1./HetCN - 1./LCN)                      &
                                   + AnaPartition * (1./HetCN - 1./RCN)) )              
                       

            Me%Matrix(PI, AMI)  = Me%Matrix(PI, AMI) + DTDay                        & 
                                   * ImobilizationRateNH4 *(HetEF/HetCP)           &
                                   * 1./((1./HetCN - 1./LCN ) * (1./AnaPartition)  &
                                   + (1./HetCN - 1./RCN))


            Me%Matrix(PI, NII)  = Me%Matrix(PI, NII) + DTDay                        &
                                   * ImobilizationRateNO3 *(HetEF/HetCP)           &
                                   *( 1./((1./HetCN - 1./LCN)                      &
                                   + AnaPartition * (1./HetCN - 1./RCN)) )              
                                                


            Me%Matrix(PI, NII)  = Me%Matrix(PI, NII) + DTDay                        &
                                   * ImobilizationRateNO3  *(HetEF/HetCP )         &
                                   *( 1./((1./HetCN - 1./LCN )*(1./ AnaPartition)  &
                                   + (1./HetCN - 1./RCN)))                       

        end if


        if (seletion==4.OR. seletion==7) then

            Me%Matrix(PI, PI)  =    Me%Matrix(PI, PI) + DTDay                      &
                                    * ImobilizationRateP* (HetEF/HetCP )           & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/Anapartition)  & 
                                    + (1./HetCP - 1./RCP))) 
                                    
            Me%Matrix(PI, PI)  =    Me%Matrix(PI, PI) + DTDay                      &
                                    * ImobilizationRateP* (HetEF/HetCP )           & 
                                    * (1./( (1./HetCP - 1./LCP )                   & 
                                    + Anapartition*(1./HetCP - 1./RCP)))    
                                    
                                                            
                    
        end if

        if(seletion==1) then
    
            Real_N_ROM = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)         &
                        *Conversion) *(1./( (1./HetCN - 1./LCN )*(1/AnaPartition)  &
                        +  (1./HetCN - 1./RCN)))

            Real_P_ROM = ((ImobilizationRateP*P)*Conversion)                       &
                        * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)              &
                         +  (1./HetCP - 1./RCP)))



            Real_N_LOM = ((ImobilizationRateNH4*AM+ImobilizationRateNO3*N)         &
                        *Conversion) *(1./( (1./HetCN - 1./LCN )                   &
                        +  AnaPartition*(1./HetCN - 1./RCN)))


            Real_P_LOM = ((ImobilizationRateP*P)*Conversion)                       &
                        * (1./( (1./HetCP - 1./LCP ) +                             &
                        AnaPartition*(1./HetCP - 1./RCP)))



            if( Real_N_ROM<Real_P_ROM) then

                Me%Matrix(PI, AMI)  = Me%Matrix(PI, AMI) + DTDay                 & 
                                    * ImobilizationRateNH4*(HetEF/HetCP)           &
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)  &
                                    +  (1./HetCN - 1./RCN))) 
                
                Me%Matrix(PI, NII)  = Me%Matrix(PI, NII) + DTDay                 &
                                    * ImobilizationRateNO3*(HetEF/HetCP)           & 
                                    * (1./( (1./HetCN - 1./LCN )*(1/AnaPartition)  & 
                                    +  (1./HetCN - 1./RCN)))    

            else

                Me%Matrix(PI, PI)  = Me%Matrix(PI, PI) + DTDay                   &
                                    * ImobilizationRateP *(HetEF/ HetCP)           & 
                                    * (1./( (1./HetCP - 1./LCP )*(1/AnaPartition)  &
                                    + (1./HetCP - 1./RCP)))

            end if

            if( Real_N_LOM<Real_P_LOM) then


                        
                Me%Matrix(PI, AMI)  = Me%Matrix(PI, AMI) + DTDay                 & 
                                    * ImobilizationRateNH4*(HetEF/HetCP)           &
                                    * (1./( (1./HetCN - 1./LCN )                   &
                                    +  AnaPartition*(1./HetCN - 1./RCN))) 
                
                Me%Matrix(PI, NII)  = Me%Matrix(PI, NII) + DTDay                 &
                                    * ImobilizationRateNO3  *(HetEF/HetCP)         & 
                                    * (1./( (1./HetCN - 1./LCN )                   & 
                                    +  AnaPartition*(1./HetCN - 1./RCN)))   

            else


                Me%Matrix(PI, PI)  = Me%Matrix(PI, PI) + DTDay                   &
                                    * ImobilizationRateP*(HetEF/HetCP)             & 
                                    * (1./( (1./HetCP - 1./LCP )                   &
                                    + AnaPartition*(1./HetCP - 1./RCP)))

            end if
    
        end if


        !Sink  : imobilization of Soluble Phosphorus by heterotrophics


        if (seletion == 2 .OR. seletion == 5 .OR. seletion == 8) then

             ! RealP imobilization

            Me%Matrix(PI, LOM_CI) = Me%Matrix(PI, LOM_CI) - DTDay               &
                                      * (LOM_C_Srate * (1./HetCP - 1./LCP))    &
                                      *(1/conversion)
 
            Me%Matrix(PI, ROM_CI) = Me%Matrix(PI, ROM_CI) - DTDay               &
                                      * (ROM_C_Srate * (1./HetCP - 1./RCP))    &
                                      *(1/conversion)

        end if


        if (seletion == 1 .OR. seletion == 4 .OR. seletion == 7) then

             ! Potential P immobilization  
       
            Me%Matrix(PI, PI) = Me%Matrix(PI, PI) - DTDay                    &
                           * ImobilizationRateP 

        end if

        !!! passagem de fosforo soluvel para inorganico, ou seja fixação
        if( P> PF*(pai/(1-pai)) )  then


            Me%Matrix(PI, PI)  = Me%Matrix(PI, PI)-DTDay*100/(SoilDens*1)

            Me%Matrix(PI, PFI) = Me%Matrix(PI, PFI)+DTDay*(pai/(1-pai))*100/(SoilDens*1)

        end if


        !!! escolha entre a passagem meramente quimica ou com a intervenção dos microorganismos

        if (Me%PropCalc%Sol_Bacteria) then 


            Me%Matrix(PI, PFI) = Me%Matrix(PI, PFI) + DTDay                       &
                           * SolubilizingRate*(1-0.1)                          
                                    

        else

            if( P< PF*(pai/(1-pai)) )  then


                Me%Matrix(PI, PI)  = Me%Matrix(PI, PI)-0.1*DTDay*100/(SoilDens*1)

                Me%Matrix(PI, PFI) = Me%Matrix(PI, PFI)+0.1*DTDay*(pai/(1-pai))*100/(SoilDens*1)

            end if

        end if

 
         !Independent term
         Me%IndTerm(PI) = Me%ExternalVar%Mass(PI,index) 
       


    end subroutine SQInorganicPhosphorusSoluble  

  !----------------------------------------------------------------------------


       subroutine SQInorganicPhosphorusFix (index)

       !Arguments---------------------------------------------------------------
        integer                   :: index


       !Local-------------------------------------------------------------------

        real    :: DTDay
        integer :: PI,PFI
        real    :: pai, P,PF
        real    :: solubilizingrate
        real    :: solCEf,SolCP,solpartition
        real    :: conversion
        real    :: SoilDens

       !------------------------------------------------------------------------

        DTDay               = Me%DTDay!
        PFI                 = Me%PropIndex%Inorganic_P_fix
        PI                  = Me%PropIndex%Inorganic_P_soluble
        P                   = Me%ExternalVar%Mass(PI, index)
        PF                  = Me%ExternalVar%Mass(PFI, index)
        solubilizingrate    = Me%SpecificRates%Solubilizing%Value   
        solpartition        = Me%Solpartition
        solCEf              = Me%Microorganisms%Heterotrophs%EficiencyC
        SolCP               = Me%Microorganisms%Sols%CPRatio
        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)
        SoilDens            = Me%ExternalVar%SoilDryDensity(index)
    
        PAI  = Me%EXternalVar%PAI(index)


        Me%Matrix(PFI, PFI) = 1.   


        if( P> PF*(pai/(1-pai)) )  then


            Me%Matrix(PFI, PI)  = Me%Matrix(PFI, PI)+DTDay*100/SoilDens

            Me%Matrix(PFI, PFI) = Me%Matrix(PFI, PFI)-DTDay*(pai/(1-pai))*100/(SoilDens*1)

        end if

        !!! escolha entre a passagem meramente quimica ou com a intervenção dos microorganismos

        if (Me%PropCalc%Sol_Bacteria) then 


            Me%Matrix(PFI, PFI) = Me%Matrix(PFI, PFI) - DTDay                     &
                          * SolubilizingRate                                   &
                          * (1-0.1)/ SolCP                      

        else

            if( P< PF*(pai/(1-pai)) )  then


                Me%Matrix(PI, PI)  = Me%Matrix(PI, PI)+0.1*DTDay*100/(SoilDens*1)

                Me%Matrix(PI, PFI) = Me%Matrix(PI, PFI)-0.1*DTDay*(pai/(1-pai))*100/(SoilDens*1)

            end if

        end if

  
  
         !Independent term
        Me%IndTerm(PFI) = Me%ExternalVar%Mass(PFI,index) 
       


    end subroutine SQInorganicPhosphorusFix 

   !----------------------------------------------------------------------------


    subroutine SQSol_Bacteria(index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN) :: index
        !------------------------------------------------------------------------
            
            call SQsol_Carbon      (index)

            call SQsol_Nitrogen   (index)

            call SQsol_Phosphorus  (index)


        !------------------------------------------------------------------------

    end subroutine SQSol_Bacteria
    !---------------------------------------------------------------------------



    subroutine SQsol_Carbon (index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN) :: index
        !------------------------------------------------------------------------
       !Local-------------------------------------------------------------------

        integer     :: solCI,PFI
        real        :: SolDeathRate, solubilizingrate
        real        :: DTDay
        real        :: conversion
        real        :: solCEf

       !----------------------------------------------------------------------------

        SolDeathRate      = Me%SpecificRates%sol%Value
        SolubilizingRate  = Me%SpecificRates%Solubilizing%Value   
       
        SolCI             = Me%PropIndex%SolC
        PFI               = Me%PropIndex%Inorganic_P_Fix
        solCEf            = Me%Microorganisms%Heterotrophs%EficiencyC             
       
        Conversion        = Me%ExternalVar%DissolvedToParticulate (index)
        DTDay             = Me%DTDay


        Me%Matrix(SolCI, SolCI)   =1

         
        !Sink   : Heterotrophic death

        if (.NOT. Me%Microorganisms%Sols%LogicalMinumumPOP)                 &
           Me%Matrix(SolCI, SolCI) = Me%Matrix(SolCI, SolCI) - DTDay * SolDeathRate


        ! source : Labil 

        Me%Matrix(SolCI, PFI) = Me%Matrix(SolCI, PFI) + DTDay                       &
                              * solubilizingRate * Conversion                       

                              
                               
        ! sink   : excrecion of CO2

        Me%Matrix(SolCI, PFI) = Me%Matrix(SolCI, PFI) - DTDay                       &
                              * solubilizingRate * Conversion*SolCEf                


                              
        !Independent term
        Me%IndTerm(solCI)     = Me%ExternalVar%Mass(solCI, index) 
                              
                              
                               
        !------------------------------------------------------------------------

     end subroutine SQsol_Carbon

        !------------------------------------------------------------------------


        subroutine SQsol_Nitrogen (index)
        
        !Arguments---------------------------------------------------------------
        integer, intent(IN) :: index
        !------------------------------------------------------------------------
       !Local-------------------------------------------------------------------

        integer     :: solNI,PFI
        real        :: SolDeathRate, solubilizingrate
        real        :: DTDay
        real        :: conversion
        real        :: solCN,anapartition
        real        :: LCN,RCN
        real        :: solCEf
       !----------------------------------------------------------------------------

        PFI                 = Me%PropIndex%Inorganic_P_fix
            
        SolNI               = Me%PropIndex%SolN
        SolDeathRate        = Me%SpecificRates%sol%Value
        solubilizingrate    = Me%SpecificRates%Solubilizing%Value   
        solCEf              = Me%Microorganisms%Heterotrophs%EficiencyC

        SolCN               = Me%Microorganisms%Sols%CNRatio
        LCN                 = Me%LabiOM_CN_Ratio
        RCN                 = Me%RefractOM_CN_Ratio
        AnaPartition        = Me%AnaerobicPartition
        Conversion        = Me%ExternalVar%DissolvedToParticulate (index)
        DTDay             = Me%DTDay


        Me%Matrix(SolNI, SolNI)    = 1
         
        !Sink   : Heterotrophic death

        if (.NOT. Me%Microorganisms%Sols%LogicalMinumumPOP)                         &
           Me%Matrix(SolNI, SolNI) = Me%Matrix(SolNI, SolNI) - DTDay * SolDeathRate/SolCN

        !Sources: Labil N
        Me%Matrix(SolNI, PFI) = Me%Matrix(SolNI, PFI) + DTDay                       &
                              * solubilizingRate * Conversion *(1/LCN)              &
                              * (0.1 /( (anapartition + 1. ))) 
        !Sources: Refract N
        Me%Matrix(SolNI, PFI) = Me%Matrix(SolNI, PFI) + DTDay                       &
                              * SolubilizingRate * Conversion *(1/RCN)              &
                              * (0.1 /( (1./anapartition + 1. ) )) 

        !Sink: Excretion
        Me%Matrix(SolNI, PFI) = Me%Matrix(SolNI, PFI) + DTDay                       &
                              * solubilizingRate * Conversion *solCEf*(1/LCN)       &
                              * (0.1 /( (anapartition + 1. )))  
                              
        Me%Matrix(SolNI, PFI) = Me%Matrix(SolNI, PFI) + DTDay                       &
                              * SolubilizingRate * Conversion *(1/RCN)              &
                              *solCEf* (0.1 /( (1./anapartition + 1. ) ))                                   

                              
        !Independent term
        Me%IndTerm(solNI)     = Me%ExternalVar%Mass(solNI, index) 



        !------------------------------------------------------------------------

      end subroutine SQsol_Nitrogen

        !------------------------------------------------------------------------



        subroutine SQsol_Phosphorus (index)
        
        !Arguments---------------------------------------------------------------
        integer, intent(IN) :: index
        !------------------------------------------------------------------------
       !Local-------------------------------------------------------------------

        integer     :: solPI,PFI
        real        :: SolDeathRate, solubilizingrate
        real        :: DTDay
        real        :: conversion
        real        :: solCP
        real        :: LCP,RCP
        real        :: Anapartition
       !----------------------------------------------------------------------------

        PFI                 = Me%PropIndex%Inorganic_P_fix
            
        SolPI               = Me%PropIndex%SolP
        SolDeathRate        = Me%SpecificRates%sol%Value
        solubilizingrate    = Me%SpecificRates%Solubilizing%Value   


        SolCP               = Me%Microorganisms%Sols%CPRatio
        LCP                 = Me%LabilOM_CP_Ratio
        RCP                 = Me%RefractOM_CP_Ratio
        AnaPartition        = Me%AnaerobicPartition
        Conversion        = Me%ExternalVar%DissolvedToParticulate (index)
        DTDay             = Me%DTDay



        Me%Matrix(SolPI, SolPI)    = 1


        !Sink   : Heterotrophic death

        if (.NOT. Me%Microorganisms%Sols%LogicalMinumumPOP)                         &
           Me%Matrix(SolPI, SolPI) = Me%Matrix(SolPI, SolPI) - DTDay * SolDeathRate/SolCP


        !Sources: Labil P
        Me%Matrix(SolPI, PFI) = Me%Matrix(SolPI, PFI) + DTDay                       &
                              * solubilizingRate * Conversion *(1/LCP)              &
                              * (0.1 /( (anapartition + 1. ))) 
        !Sources: Refract P
        Me%Matrix(SolPI, PFI) = Me%Matrix(SolPI, PFI) + DTDay                       &
                              * SolubilizingRate * Conversion *(1/RCP)              &
                              * (0.1 /( (1./anapartition + 1. ) )) 

        !Independent term
        Me%IndTerm(solPI)     = Me%ExternalVar%Mass(solPI, index) 

       !----------------------------------------------------------------------------

       end subroutine SQsol_Phosphorus

       !----------------------------------------------------------------------------



    subroutine CorrectOM (index)

       !Local-------------------------------------------------------------------

        real    :: conversion
        integer :: index
        real    :: AnaC,AnaN,HetC,HetN,AnaP,HetP,solC,solN,solP
        real    :: HetCN,HetCP,AnaCN,AnaCP,solCN,solCP
        real    :: N,P
        real    :: excessoHN,excessoHP,excessoAN,ExcessoAP,Excesso_solP,Excesso_solN 
        integer :: AnaCI,AnaNI,AnaPI,HetCI,HetNI,HetPI,PI,NI,solCI,solNI,solPI

       !----------------------------------------------------------------------------

!!! indices das propriedades


        AnaCI = Me%PropIndex%anaerobicC
        AnaNI = Me%PropIndex%anaerobicN
        AnaPI = Me%PropIndex%anaerobicP

        HetCI = Me%PropIndex%heterotrophicC
        HetNI = Me%PropIndex%heterotrophicN
        HetPI = Me%PropIndex%heterotrophicP

        NI = Me%PropIndex%ammonia 
        PI = Me%PropIndex%Inorganic_P_soluble

        solCI = Me%PropIndex%solC 
        solNI = Me%PropIndex%solN
        solPI = Me%PropIndex%solP           

            
        conversion = Me%ExternalVar%DissolvedToParticulate (index) 

        if (Me%PropCalc%Carbon) then

            AnaC = Me%ExternalVar%Mass(anaCI, index)
            HetC = Me%ExternalVar%Mass(HetCI, index)
            
            if (Me%PropCalc%Sol_Bacteria) then
                solC = Me%ExternalVar%Mass(solCI, index) 
            endif

        endif

        if (Me%PropCalc%Nitrogen) then

            AnaN = Me%ExternalVar%Mass(AnaNI, index)
            HetN = Me%ExternalVar%Mass(HetNI, index)
            N = Me%ExternalVar%Mass(NI, index)
 
            if (Me%PropCalc%Sol_Bacteria) then

                solN = Me%ExternalVar%Mass(solNI, index)

            endif

        endif

        if (Me%PropCalc%Phosphorus) then

            AnaP = Me%ExternalVar%Mass(AnaPI, index)
            HetP = Me%ExternalVar%Mass(HetPI, index)
            P = Me%ExternalVar%Mass(PI, index)

            if (Me%PropCalc%Sol_Bacteria) then

                solP = Me%ExternalVar%Mass(solPI, index)

            endif
        
        endif

        HetCN = Me%Microorganisms%Heterotrophs%CNratio
        HetCP = Me%Microorganisms%Heterotrophs%CPratio

        AnaCN = Me%Microorganisms%Anaerobic%CNRatio
        AnaCP = Me%Microorganisms%Anaerobic%CPRatio

        solCN = Me%Microorganisms%sols%CNRatio
        solCP = Me%Microorganisms%sols%CPRatio



!! Excess organic matter in heterotrophicP

        if(HetC/ HetCP< HetP)then

            excessoHP = (HetP-(HetC/ HetCP))
            HetP = HetP  - excessoHP

            P = P+(1/conversion)*excessoHP


        end if

!!! Excess organic matter in heterotrophicN


        if(HetC/ HetCN< HetN)then 

            excessoHN = (HetN-(HetC/ HetCN))

            HetN = HetN  - excessoHN

            N = N + (1/conversion)*excessoHN


        end if


!!! excess organic matter in anaerobic N



        excessoAN = (AnaN-(AnaC/ AnaCN))

        AnaN = AnaN  - excessoAN

        N = N+(1/conversion)*excessoAN





!!! excess organic matter in anaerobicP


        ExcessoAP  = (AnaP-(AnaC/ AnaCP))
        AnaP = AnaP  -  ExcessoAP

        P = P+(1/conversion)*ExcessoAP




!!! excess organic matter in sol N

        if(solC/solCN< solN)then

            excesso_solN = (solN-(solC/ solCN))

            solN = solN  - excesso_solN

            N = N+(1/conversion)*excesso_solN


        end if

!!! excess organic matter in solP

        if(solC/solCP< solP)then

            Excesso_solP  = (solP-(solC/ solCP))
            solP = solP  -  Excesso_solP

            P = P+(1/conversion)*Excesso_solP


        end if

        if (Me%PropCalc%Nitrogen) then

            Me%ExternalVar%Mass(AnaNI, index) = AnaN 
            Me%ExternalVar%Mass(HetNI, index) = HetN
            Me%ExternalVar%Mass(NI, index)    = N
            
            if (Me%PropCalc%Sol_Bacteria) then
                
                Me%ExternalVar%Mass(solNI, index) = solN 
            
            endif
        endif

        if (Me%PropCalc%Phosphorus) then

            Me%ExternalVar%Mass(AnaPI, index) = AnaP 
            Me%ExternalVar%Mass(HetPI, index) = HetP 
            Me%ExternalVar%Mass(PI, index)    = P 

            if (Me%PropCalc%Sol_Bacteria) then

                Me%ExternalVar%Mass(solPI, index) = solP
                
            endif 

        endif




    end subroutine CorrectOM

    !----------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !----------------------------------------------------------------------------
    subroutine KillSedimentQuality(SedimentQualityID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: SedimentQualityID
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer :: ready_              
        integer :: STAT_CALL

        !Local-------------------------------------------------------------------
        integer :: STAT_
        integer                                     :: nUsers 

        !------------------------------------------------------------------------                      

        STAT_ = UNKNOWN_

        call Ready(SedimentQualityID, ready_)

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mSEDIMENTQUALITY_,  Me%InstanceID)

cd2 :       if (nUsers == 0) then


cd3 :           if(Me%ObjLUD /= 0) then
                    call KillLUD(Me%ObjLUD, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Subroutine KillSedimentQuality; module ModuleSedimentQuality. ERR01.'
                end if cd3


                deallocate(Me%IndTerm, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'Subroutine Kill_SedimentQuality; module ModuleSedimentQuality. ERR02.'
                nullify(Me%IndTerm)


                deallocate(Me%Matrix, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'Subroutine Kill_SedimentQuality; module ModuleSedimentQuality. ERR03.'
                nullify(Me%Matrix)


cd4 :           if (associated(Me%NewMass)) then
                    deallocate(Me%NewMass, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Subroutine Kill_SedimentQuality; module ModuleSedimentQuality. ERR04.'
                    nullify(Me%NewMass)
                end if cd4


                call DeallocateInstance 

                SedimentQualityID = 0


                STAT_ = SUCCESS_

            end if cd2
        
        else cd1
        
            STAT_ = ready_
        
        end if cd1


        if (present(STAT)) STAT = STAT_
    !------------------------------------------------------------------------

    end subroutine KillSedimentQuality
    !------------------------------------------------------------------------


    !------------------------------------------------------------------------
    subroutine DeallocateInstance 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_SedimentQuality), pointer           :: AuxObjSedimentQuality
        type (T_SedimentQuality), pointer           :: PreviousObjSedimentQuality

        !Updates pointers
        if (Me%InstanceID == FirstObjSedimentQuality%InstanceID) then
            FirstObjSedimentQuality        => FirstObjSedimentQuality%Next
        else
            PreviousObjSedimentQuality => FirstObjSedimentQuality
            AuxObjSedimentQuality      => FirstObjSedimentQuality%Next
            do while (AuxObjSedimentQuality%InstanceID /= Me%InstanceID)
                PreviousObjSedimentQuality => AuxObjSedimentQuality
                AuxObjSedimentQuality      => AuxObjSedimentQuality%Next
            enddo

            !Now update linked list
            PreviousObjSedimentQuality%Next => AuxObjSedimentQuality%Next

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


    !----------------------------------------------------------------------------
    subroutine Ready (SedimentQualityID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: SedimentQualityID
        integer                                     :: ready_
        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (SedimentQualityID > 0) then
            
            call LocateObjSedimentQuality (SedimentQualityID)
            ready_ = VerifyReadLock (mSEDIMENTQUALITY_, Me%InstanceID)

        else

            ready_ = OFF_ERR_

        end if cd1
        !----------------------------------------------------------------------

    end subroutine Ready
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine LocateObjSedimentQuality (SedimentQualityID)

        !Arguments-------------------------------------------------------------
        integer                                     :: SedimentQualityID
    !--------------------------------------------------------------------------

        Me => FirstObjSedimentQuality
        do while (associated (Me))
            if (Me%InstanceID == SedimentQualityID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))   &
            stop 'ModuleSedimentQuality - LocateObjSedimentQuality - ERR01'
    !--------------------------------------------------------------------------
    
    end subroutine LocateObjSedimentQuality
    !--------------------------------------------------------------------------


end module ModuleSedimentQuality

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
