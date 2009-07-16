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

Module ModuleSedimentQuality
!to_change - coisas a mudar
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
    private ::          CalcPopulation          
    private ::      LogicalImobilization    
    private ::      LogicalLimitation
    private ::      StartSedimentQualityIteration
    private ::      SedimentQualityCoefCalculation    
    private ::          Oxygen
    private ::          Carbon
    private ::              LabilOrganicCarbon
    private ::              RefractOrganicCarbon   
    private ::              HetrotrophicC          
    private ::              AutotrotrophicC        
    private ::              AnaerobicC
    private ::          Nitrogen
    private ::              Ammonia
    private ::              Nitrate
    private ::              LabilOrganicNitrogen   
    private ::              RefractOrganicNitrogen 
    private ::              HetrotrophicN          
    private ::              AutotrotrophicN        
    private ::              AnaerobicN             
    private ::              Ngas
    private ::      SystemResolution
   
    private ::      SedimentQualityRatesCalculation



    !Destructor
    public  ::  KillSedimentQuality                                                     
    private ::      DeallocateInstance


    !Management
    private ::      Ready
    private ::          LocateObjSedimentQuality
    

    !Types---------------------------------------------------------------------
    type       T_PropIndex
        integer :: HetrotrophicN                    = null_int
        integer :: HetrotrophicC                    = null_int
        integer :: AutotrotrophicN                  = null_int
        integer :: AutotrotrophicC                  = null_int
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
    end type T_PropIndex


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
        logical :: Carbon     = OFF
        logical :: Nitrogen   = OFF
        logical :: Oxygen     = OFF
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
    end type T_External

    type        T_Coeficients
        real                                :: ActivationE          = null_real     !Activation energy
        real                                :: Acoef                = null_real     !Specific coeficient
        real                                :: optimumTemperature   = null_real     !Optimum temperature
        real                                :: Value                = null_real     !specific rate value     
        integer                             :: RateIndex            = null_int      !Rate Index number
    end type    T_Coeficients


    !all the specific rates that are explicitly calculated
    type       T_Rates
        type(T_Coeficients)                :: Labil_OM_C
        type(T_Coeficients)                :: RefractOM_C
        type(T_Coeficients)                :: AmmoniaToNitrate
        type(T_Coeficients)                :: AmmoniaImobilization         
        type(T_Coeficients)                :: NitrateToNgas
        type(T_Coeficients)                :: NitrateImobilization
        type(T_Coeficients)                :: Heterotrophs
        type(T_Coeficients)                :: Autotrophs
        type(T_Coeficients)                :: Anaerobic      
    end type T_Rates

    type        T_Constants
        real                                :: CNRatio          = null_real     !Microorganisms C/N ratio
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
        logical                                             :: Imobilization        = OFF
        logical                                             :: NLimitation          = OFF 
        real                                                :: Partition            = null_real
        real                                                :: AnaerobicPartition   = null_real
        double precision,       pointer, dimension(:,:)     :: Matrix
        real,                   pointer, dimension(:  )     :: IndTerm
        real,                   pointer, dimension(:  )     :: NewMass          !Used with Explicit method
        
        real                                                :: DTDay                = null_real
        real                                                :: DTSecond             = null_real
        
        real                                                :: Aerobiose            = null_real
        real                                                :: Anaerobiose          = null_real
        real                                                :: LabiOM_CN_Ratio      = null_real
        real                                                :: RefractOM_CN_Ratio   = null_real
        real                                                :: Oxygen               = null_real
        real                                                :: Hydrogen             = null_real

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
    real,    parameter :: Boltzman                 = 1.383E-23      ![J ºK-1.]
    real,    parameter :: Planck                   = 6.63E-34       ![J s]
    real,    parameter :: UnivGC                   = 1.99E-3        ![Kca mole-1. ºK-1.]


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
      
            if (Me%PropCalc%Nitrogen) then

                Logicalequa(PIndex%HetrotrophicN      )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%AutotrotrophicN    )=.true.
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
       
            endif

       
            if (Me%PropCalc%Carbon) then

                Logicalequa(PIndex%HetrotrophicC      )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%AutotrotrophicC    )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%AnaerobicC         )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%Labil_OM_C         )=.true.
                countequa = countequa + 1.

                Logicalequa(PIndex%RefractOM_C        )=.true.
                countequa = countequa + 1.
       
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
            
            NewEquaRateFlux%Prev                      => Me%LastEquaRateFlux
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
            Me%PropIndex%HetrotrophicC                  = Me%Prop%IUB
            Me%Microorganisms%Heterotrophs%MicroIndex   = Me%Prop%IUB

            Me%Prop%IUB                                 = Me%Prop%IUB + 1
            Me%PropIndex%AutotrotrophicC                = Me%Prop%IUB
            Me%Microorganisms%Autotrophs%MicroIndex     = Me%Prop%IUB

            Me%Prop%IUB                                 = Me%Prop%IUB + 1
            Me%PropIndex%AnaerobicC                     = Me%Prop%IUB
            Me%Microorganisms%Anaerobic%MicroIndex      = Me%Prop%IUB    

            Me%Prop%IUB                                 = Me%Prop%IUB + 1
            Me%PropIndex%Labil_OM_C                     = Me%Prop%IUB

            Me%Prop%IUB                                 = Me%Prop%IUB + 1
            Me%PropIndex%RefractOM_C                    = Me%Prop%IUB
            
         
        endif  

        !Nitrogen index number
        if (Me%PropCalc%Nitrogen) then
            Me%Prop%IUB                     = Me%Prop%IUB + 1
            Me%PropIndex%HetrotrophicN      = Me%Prop%IUB

            Me%Prop%IUB                     = Me%Prop%IUB + 1
            Me%PropIndex%AutotrotrophicN    = Me%Prop%IUB

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
        
            
        endif    


        !Oxygen index number -> The oxygen is always calculated.
        Me%Prop%IUB         = Me%Prop%IUB + 1
        Me%PropIndex%Oxygen = Me%Prop%IUB


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

        !Verifica se se pretende calcular usando um metodo IMPLICITO
        call GetData(Me%CalcMethod%ImplicitMethod                   ,   &
                     Me%ObjEnterData, flag                          ,   &
                     SearchType   = FromFile                        ,   &
                     keyword      = 'IMPLICIT'                      ,   &
                     ClientModule = 'ModuleSedimentQuality'         ,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'Subroutine SQReadCalcOptions; Module ModuleSedimentQuality. ERR03.' 

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
                call SetError(FATAL_, INTERNAL_, 'SQReadFileConstants - ModuleSedimentQuality - ERR01') 
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
                call SetError(FATAL_, INTERNAL_, 'SQReadFileConstants - ModuleSedimentQuality - ERR01') 
            end if cd4
        end do do2

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = status) 
        if (status /= SUCCESS_) &
            call SetError(FATAL_, INTERNAL_, 'ConstructPropertyList - ModuleSurface - ERR01')

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

            end select
        
        end if             
       
        if (Me%PropCalc%Nitrogen) then

            Select case (blockpointer)

            Case (6)
                RateinProgress  =>  Me%SpecificRates%AmmoniaToNitrate
                block_begin     =   '<begin_AmmoniaToNitrate_Rate>'
                block_end       =   '<end_AmmoniaToNitrate_Rate>'
                return

            Case (7)
                RateinProgress  =>  Me%SpecificRates%AmmoniaImobilization 
                block_begin     =   '<begin_AmmoniaImobilization_Rate>'
                block_end       =   '<end_AmmoniaImobilization_Rate>'
                return

            Case (8)
                RateinProgress  =>  Me%SpecificRates%NitrateToNgas
                block_begin     =   '<begin_NitrateToNgas_Rate>'
                block_end       =   '<end_NitrateToNgas_Rate>'
                return

            Case (9)        
                RateinProgress  =>  Me%SpecificRates%NitrateImobilization               
                block_begin     =   '<begin_NitrateImobilization_Rate>'
                block_end       =   '<end_NitrateImobilization_Rate>'
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
                        keyword         ='Temperature'           ,   &
                        ClientModule    = 'ModuleSedimentQuality',   &
                        default         = 0.                     ,   &
                        STAT            = status)
       
        if (status /= SUCCESS_)                                      &
           call SetError(FATAL_, INTERNAL_, 'ConstructRate - ModuleSedimentQuality - ERR03') 
        
        RateinProgress%optimumTemperature = value

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
                            ExplicitMethod,                    &
                            ImplicitMethod,                    &
                            SemiImpMethod, STAT) 

        !Arguments-------------------------------------------------------------
        integer                        :: SedimentQualityID
        integer, optional, intent(OUT) :: STAT
        logical, optional, intent(OUT) :: Nitrogen, Oxygen, Carbon
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
    subroutine GetPropIndex(SedimentQualityID,  HetrotrophicN,              &
                                                HetrotrophicC,              &
                                                AutotrotrophicN,            &
                                                AutotrotrophicC,            &
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
                                                STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: SedimentQualityID

        integer, optional, intent(OUT)      :: HetrotrophicN
        integer, optional, intent(OUT)      :: HetrotrophicC
        integer, optional, intent(OUT)      :: AutotrotrophicN
        integer, optional, intent(OUT)      :: AutotrotrophicC
        integer, optional, intent(OUT)      :: AnaerobicN
        integer, optional, intent(OUT)      :: AnaerobicC
        
        integer, optional, intent(OUT)      :: Labil_OM_C
        integer, optional, intent(OUT)      :: Labil_OM_N

        integer, optional, intent(OUT)      :: RefractOM_C
        integer, optional, intent(OUT)      :: RefractOM_N
        
        integer, optional, intent(OUT)      :: Ammonia
        integer, optional, intent(OUT)      :: Nitrate
        integer, optional, intent(OUT)      :: Ngas
        
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
            
            if (present(HetrotrophicN       )) HetrotrophicN        = Me%PropIndex%HetrotrophicN
            if (present(HetrotrophicC       )) HetrotrophicC        = Me%PropIndex%HetrotrophicC
            if (present(AutotrotrophicN     )) AutotrotrophicN      = Me%PropIndex%AutotrotrophicN
            if (present(AutotrotrophicC     )) AutotrotrophicC      = Me%PropIndex%AutotrotrophicC
            if (present(AnaerobicN          )) AnaerobicN           = Me%PropIndex%AnaerobicN
            if (present(AnaerobicC          )) AnaerobicC           = Me%PropIndex%AnaerobicC
            
            if (present(Labil_OM_C          )) Labil_OM_C           = Me%PropIndex%Labil_OM_C
            if (present(Labil_OM_N          )) Labil_OM_N           = Me%PropIndex%Labil_OM_N
            
            if (present(Labil_OM_C          )) RefractOM_C          = Me%PropIndex%RefractOM_C
            if (present(Labil_OM_N          )) RefractOM_N          = Me%PropIndex%RefractOM_N

            if (present(Ammonia             )) Ammonia              = Me%PropIndex%Ammonia
            if (present(Nitrate             )) Nitrate              = Me%PropIndex%Nitrate
            if (present(Ngas                )) Ngas                 = Me%PropIndex%Ngas
            
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
        integer, optional,  intent(OUT)             :: STAT
         
        !External----------------------------------------------------------------
        integer                                     :: index
        integer                                     :: ready_   
               
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
                call HydrogenCalculation           
                
                call PotentialRatesCalc             (index)
                call LogicalImobilization           (index           )

                if (Me%Imobilization) then
                    call LogicalLimitation          (index           )
                end if
                
                call SedimentQualityCoefCalculation (index           )
                
                call SystemResolution               (index           )
                
                !The rates can just be calculated if the rate flux is associated
                !In the case that this module is used by the lagrangian module
                !the rate fluxes are not calculated
                if (associated(Me%FirstEquaRateFlux)) then
                    call SedimentQualityRatesCalculation   (index)
                endif
            endif    
            
            Me%Matrix    (:, :)     = 0.
            Me%Imobilization        = OFF
            Me%NLimitation          = OFF
            Me%Partition            = 0.
            Me%AnaerobicPartition   = 0.

            end do do1


            nullify(Me%ExternalVar%Temperature   )
            nullify(Me%ExternalVar%Mass          )
            nullify(Me%ExternalVar%ThetaF        )
                     
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine SedimentQuality
    !----------------------------------------------------------------------------


    ! Calculation of the aerobic and anaerobic functions for Microbiological action.
    ! When a soil atmosphere algoritm is implemented this will no longer be necessary
    !----------------------------------------------------------------------------
    subroutine AnaerobioseCalculation (ThetaF, index)

        !Arguments---------------------------------------------------------------

        integer             , intent(IN)                    :: index  
        real                , dimension(:  ), intent(IN)    :: ThetaF

        !Local-------------------------------------------------------------------

        real :: thetaef
        real :: m
        real :: b
        !------------------------------------------------------------------------
        
        thetaef = ThetaF(index) * 100.

        if      (thetaef .LE. 60                     )   then
                    Me%Aerobiose     = 1.
                    Me%Anaerobiose   = 0.01

        elseif  (thetaef .GE. 60 .AND. thetaef .LE. 70)   then
                    
                    m                = 0.6 / (60.-70.)
                    b                = 1. - m * 60.
                    Me%Aerobiose     = m * thetaef + b
                    
                    m                = 0.13 / (80.-60.)
                    b                = 0
                    Me%Anaerobiose   = m * thetaef + b
        
        elseif  (thetaef .GE. 70 .AND. thetaef .LE. 80)   then
        
                    m                = 0.3 / (70.- 80.)
                    b                = 0.4 - m * 70.
                    Me%Aerobiose     = m * thetaef + b
                    
                    m                = 0.13 / (80.-60.)
                    b                = 0
                    Me%Anaerobiose   = m * thetaef + b
        else                
                    Me%Aerobiose     = 0.1
                    Me%Anaerobiose   = 1.
        end if  
        
        !------------------------------------------------------------------------

    end subroutine AnaerobioseCalculation
    !----------------------------------------------------------------------------


    ! Calculation of the saturation [O2] in water (saturation atm pressure is used)
    ! Anaerobiose calculatons will be made using the anaerobiose function
    !----------------------------------------------------------------------------
    subroutine OxygenCalculation (index)

        !Arguments---------------------------------------------------------------

        integer                 , intent(IN)    :: index  

        !Local-------------------------------------------------------------------
        real    :: Temp
        !------------------------------------------------------------------------
       
        Temp                        = Me%ExternalVar%Temperature(index)
        Me%Oxygen                   = OxygenSaturationHenry (Temp)
        
        !------------------------------------------------------------------------

    end subroutine OxygenCalculation 
    !----------------------------------------------------------------------------


    ! Calculation of the Hydrogen activity for specific rates calculations

    !----------------------------------------------------------------------------
    subroutine HydrogenCalculation 

        Me%Hydrogen = 1. * 10. ** (-7.0)
           
    end subroutine HydrogenCalculation
    !----------------------------------------------------------------------------
    

    ! Check if the Heterotrophs will have to immobilizate Mineral N
    !----------------------------------------------------------------------------
    subroutine LogicalImobilization (index)

        !Arguments---------------------------------------------------------------
        integer                 , intent(IN)    :: index  

        !Local-------------------------------------------------------------------
        integer :: CI                   !C property index
        integer :: NI                   !N property index
                
        real    :: LC                   !Labil OM C content
        real    :: LCN                  !Labil CN ratio
        real    :: Pot_LOM_C_Decay      !Potential OM C Decay           [massC / Time]

        real    :: RC                   !Refractary OM C content
        real    :: RCN                  !Refractary CN ratio
        real    :: Pot_ROM_C_Decay      !Potential Refract OM C Decay   [massC / Time]
        
        real    :: N                    !Labil or Refractary N content
              
        real    :: obtained             !Obtained N during OM decay    
        real    :: needed               !Needed N for the obtained C in OM decay
        !------------------------------------------------------------------------
        
        !calculate the Labil OM CN ratio
        CI  = Me%PropIndex%Labil_OM_C
        NI  = Me%PropIndex%Labil_OM_N

        LC  = Me%ExternalVar%Mass(CI, index)
        N   = Me%ExternalVar%Mass(NI, index)
        
        Me%LabiOM_CN_Ratio   = LC/N 
        LCN                  = Me%LabiOM_CN_Ratio
        
        !calculate the Refractary OM CN ratio
        CI  = Me%PropIndex%RefractOM_C
        NI  = Me%PropIndex%RefractOM_N

        RC  = Me%ExternalVar%Mass(CI, index)
        N   = Me%ExternalVar%Mass(NI, index)

        Me%RefractOM_CN_Ratio = RC/N
        RCN                   = Me%RefractOM_CN_Ratio
    
        !Calculate the OM_C decay rates
        Pot_LOM_C_Decay     = LC * Me%SpecificRates%Labil_OM_C%Value
        Pot_ROM_C_Decay     = RC * Me%SpecificRates%RefractOM_C%Value

        !Obtained N during OM decay
        obtained            = Pot_LOM_C_Decay / LCN + Pot_ROM_C_Decay / RCN
        
        !Needed N for de C obtained in the OM decay              
        needed              =   (Pot_LOM_C_Decay + Pot_ROM_C_Decay) /                  & 
                                Me%Microorganisms%Heterotrophs%CNRatio  
        !Immobilization test
        if (obtained .LT. needed ) then
            Me%Imobilization = ON    !the Heterotrophs will have to incorporate Mineral N
                            
        end if
        
        !Set Anaeribic Partition
        Me%AnaerobicPartition = (Pot_ROM_C_Decay)/(Pot_LOM_C_Decay )

        !------------------------------------------------------------------------

    end subroutine LogicalImobilization 
    !----------------------------------------------------------------------------


    ! Check if the Heterotrphs are limited by C or N
    !----------------------------------------------------------------------------
    subroutine LogicalLimitation (index)

        !Arguments---------------------------------------------------------------

        integer                 , intent(IN)    :: index  

        !Local-------------------------------------------------------------------
        
        integer :: RCI              !Refractary OM C index
        real    :: RC               !Refractary OM C concentration
        real    :: RCN              !Refractary OM CN ratio
        real    :: Pot_ROM_C_Decay  !Potential Refractary OM C decay rate

        integer :: LCI              !Labil OM C index
        real    :: LC               !Labil OM C concentration
        real    :: LCN              !Labil OM CN ratio
        real    :: Pot_LOM_C_Decay  !Potential Labil OM C decay rate

        integer :: AMI              !Ammonia index
        real    :: AM               !Ammonia concentration
        real    :: PotAMImobil      !Potential Ammonia immobilizaion rate

        integer :: NII              !Nitrate index
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


        !get the CN ratios
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

       
        !Rate that the Heterotrophs can uptake mineral N in ug kg-1 soil   
        ImobilizationS  = ( PotNIImobil + PotAMImobil ) * Conversion 


        !Rate that the Heterotrophs need extra Nitrogen (not in the OM) assuming the potential rates ug kg-1 soil
        CarbonS         = (Pot_LOM_C_Decay + Pot_ROM_C_Decay) / MCN -       & 
                           Pot_LOM_C_Decay / LCN - Pot_ROM_C_Decay / RCN    

        
        !Imobilization test        
        if (ImobilizationS .LE. CarbonS ) then
            !The organic mather decay rates are function of de N immobilization rates
            Me%NLimitation = ON    
            !The immobilization rates will be a function of Potential OM carbon decay
            Me%Partition = (Pot_ROM_C_Decay)/(Pot_LOM_C_Decay )
        
        else
            !The OM carbon decay rates will be a function of Potential immobilization rates 
            Me%Partition = (PotNIImobil) / (PotAMImobil )
        
        end if
           
        !------------------------------------------------------------------------

    end subroutine LogicalLimitation
    !----------------------------------------------------------------------------

    
    ! Calculates the Specific rates
    !----------------------------------------------------------------------------
    subroutine PotentialRatesCalc (index )

        !Arguments---------------------------------------------------------------
        integer, intent(IN)                 :: index

        !Local-------------------------------------------------------------------
        real(8)                             :: Tterm
        real                                :: Aerobiose
        real                                :: Anaerobiose
        real                                :: Oxygen
        real                                :: Hydrogen
        real                                :: Temp

        integer                             :: CLI
        integer                             :: CRI
        real                                :: CSUBST

        integer                             :: AMI
        real                                :: Ammonia

        type(T_Coeficients)     , pointer   :: Rate
        type(T_Microorganisms)  , pointer   :: Micro
        !------------------------------------------------------------------------
        
        Aerobiose   = Me%Aerobiose
        Anaerobiose = Me%Anaerobiose

        Oxygen      = Me%Oxygen 
        Hydrogen    = Me%Hydrogen
        
        Temp        = Me%ExternalVar%Temperature(index)

        CLI         = Me%PropIndex%Labil_OM_C
        CRI         = Me%PropIndex%RefractOM_C
        CSUBST      = Me%ExternalVar%Mass(CLI, index) +         &
                      Me%ExternalVar%Mass(CRI, index) 
        
        AMI         = Me%PropIndex%Ammonia
        Ammonia     = Me%ExternalVar%Mass(AMI, index) 

        Rate        => Me%SpecificRates%Labil_OM_C
        Micro       => Me%Microorganisms

        call CalcPopulation (index)

        if (Me%PropCalc%Carbon) then
             
        !Calculates the OM C specific decay Rate                 
            Tterm       =   CalcTterm (Rate, Temp)                                                 
            Rate%Value  =   Aerobiose * Tterm * Oxygen / Hydrogen**0.167        & 
                            * Micro%Heterotrophs%Population
            
        !Calculates the Refractary OM C specific decay Rate
            Rate        =>  Me%SpecificRates%RefractOM_C
            Tterm       =   CalcTterm (Rate, Temp)                    
            Rate%Value  =   Aerobiose * Tterm * Oxygen / Hydrogen**0.167        &
                            * Micro%Heterotrophs%Population
        
        !Calculates the Heterotrophs C specific decay (death) Rate
            Rate        =>  Me%SpecificRates%Heterotrophs
            Tterm       =   CalcTterm (Rate, Temp)
            Rate%Value  =   1./ Aerobiose * Tterm * Hydrogen /                   & 
                            ( Oxygen * CSUBST )  
                        
        !Calculates the Autotrophs C specific decay (death) Rate
            Rate        =>  Me%SpecificRates%Autotrophs
            Tterm       =   CalcTterm (Rate, Temp)
            Rate%Value  =   1./ Aerobiose * Tterm * Hydrogen /                   & 
                            ( Oxygen * (Ammonia + 0.5) )  
                        
        !Calculates the Anaerobic C specific decay (death) Rate
            Rate        =>  Me%SpecificRates%Anaerobic
            Tterm       =   CalcTterm (Rate, Temp)
            Rate%Value  =   1./ Anaerobiose * Tterm * Hydrogen /                 & 
                            ( Oxygen * CSUBST )        
        
        end if


        if (Me%PropCalc%Nitrogen) then
        
        !Calculates the AmmoniaToNitrate (nitrification) specific Rate  
            Rate        =>  Me%SpecificRates%AmmoniaToNitrate                   
            Tterm       =   CalcTterm (Rate, Temp)
            Rate%Value  =   Aerobiose * Tterm * Oxygen ** 0.5 / Hydrogen**0.167         &
                            * Micro%Autotrophs%Population
                                                
        !Calculates the AmmoniaImobilization specific Rate
            Rate        =>  Me%SpecificRates%AmmoniaImobilization 
            Tterm       =   CalcTterm (Rate, Temp)
            Rate%Value  =   Aerobiose * Tterm * Oxygen / Hydrogen **0.167               &
                            * Micro%Heterotrophs%Population 
                                                        
        !Calculates the NitrateToNgas specific Rate
            Rate        =>  Me%SpecificRates%NitrateToNgas       
            Tterm       =   CalcTterm (Rate, Temp)
            Rate%Value  =   Anaerobiose * Tterm * CSUBST / Hydrogen **0.167             &
                            * Micro%Anaerobic%Population

        !Calculates the NitrateImobilization specific Rate
            Rate        =>  Me%SpecificRates%NitrateImobilization 
            Tterm       =   CalcTterm (Rate, Temp)
            Rate%Value  =   Aerobiose * Tterm * Oxygen / Hydrogen **0.167               &
                            * Micro%Heterotrophs%Population
                                                        
        end if  
        !------------------------------------------------------------------------

    end subroutine PotentialRatesCalc
    !----------------------------------------------------------------------------


    ! Calculates the temperature term of the Specific rates
    !----------------------------------------------------------------------------
    function   CalcTterm (Coeficient, Temperature)
    real(8)     :: CalcTterm

        !Arguments---------------------------------------------------------------
        type(T_Coeficients) , pointer                       :: Coeficient
        real                                                :: Temperature
        
        !Local-------------------------------------------------------------------
        real    :: optimumt
        real    :: tcalc
        real    :: AE
        real    :: A
        !------------------------------------------------------------------------
        
            optimumt        = Coeficient%optimumTemperature  !rate optimum temperature
            tcalc           = Temperature
            tcalc           = 30
            if (tcalc > optimumt) tcalc = 2 * optimumt - tcalc

            AE          = Coeficient%ActivationE
            A           = Coeficient%Acoef
            tcalc       = tcalc + 273.

            CalcTterm = Boltzman * tcalc / Planck * A * &
                        exp (- AE / (UnivGC) * tcalc  ) 
        !------------------------------------------------------------------------

    end function CalcTterm
    !----------------------------------------------------------------------------


    ! Calculates the microorganisms population
    !----------------------------------------------------------------------------
    subroutine CalcPopulation (index)

        !Arguments---------------------------------------------------------------
        integer                 , intent(IN)                    :: index

        !Local-------------------------------------------------------------------
        real                           :: Carbon
        real                           :: CpopRatio
        integer                        :: propindex 
        
        type(T_Constants)  , pointer   :: Microorganisms 
        !------------------------------------------------------------------------
        
        Microorganisms    => Me%Microorganisms%Heterotrophs

            propindex   = Microorganisms%MicroIndex
            Carbon      = Me%ExternalVar%Mass(propindex, index) 
            CpopRatio   = Microorganisms%CPopRatio

            Microorganisms%Population = Carbon * CpopRatio
            
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

            if (Microorganisms%Population .LE.  Microorganisms%MinimumPop   ) then
                Microorganisms%LogicalMinumumPOP = ON
            else
                Microorganisms%LogicalMinumumPOP = OFF
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
    subroutine SedimentQualityCoefCalculation(index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)         :: index
        !----------------------------------------------------------------------

                                    call Oxygen       (index)
        if (Me%PropCalc%Carbon    ) call Carbon       (index)
        if (Me%PropCalc%Nitrogen  ) call Nitrogen     (index)
        !----------------------------------------------------------------------

    end subroutine SedimentQualityCoefCalculation
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

cd2 :       if (Me%Matrix(equa, prop) .NE. 0.)  then

cd3 :           if (equa .EQ. prop) then
                    
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
                                         - ((Me%Matrix (equa, prop) - 1.0)          &
                                         * Me%ExternalVar%Mass(prop, index))
                                 
                        Me%Matrix(equa, prop) = 1.0
                    else cd33
                        Me%IndTerm(equa) = Me%IndTerm (equa)                        &
                                         - Me%Matrix (equa, prop)                   &
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

do1 :           do while(associated(EquaRateFluxX))  
                    
                    PropRateFluxX => EquaRateFluxX%FirstPropRateFlux
                    
                    do while(associated(PropRateFluxX))  
                         
                         equa = EquaRateFluxX%ID
                         prop = PropRateFluxX%ID
                                                                    

                    if (equa.ne.prop) then
                    
                        PropRateFluxX%Field(index)  = Me%Matrix(equa, prop)         & 
                                                    * Me%ExternalVar%Mass(prop,index)
                    else
                        
                        PropRateFluxX%Field(index)  = -(1.-Me%Matrix(equa, prop))   & 
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

    !OXYGEN
    !
    !SOURCES: - 
    !SINKS:   - 
    !----------------------------------------------------------------------------
    subroutine Oxygen (index)

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
        ! vou correr com 0xigénio constante
        !------------------------------------------------------------------------

    end subroutine Oxygen
    !----------------------------------------------------------------------------
    

    !----------------------------------------------------------------------------
    subroutine Carbon(index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN) :: index
        !------------------------------------------------------------------------
                                                        
        call LabilOrganicCarbon     (index)

        call RefractOrganicCarbon   (index)

        call HetrotrophicC          (index)

        call AutotrotrophicC        (index)

        call AnaerobicC             (index)       
    !------------------------------------------------------------------------

    end subroutine Carbon
    !---------------------------------------------------------------------------


    !Labil Organic Carbon
    !
    !SOURCES: - Microorganisms death 
    !SINKS:   - Heterotrophs Aerobic and Anaerobic uptake
    !----------------------------------------------------------------------------    
    subroutine LabilOrganicCarbon (index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN)         :: index

        !Local-------------------------------------------------------------------
        integer :: AMI, NII
        
        real    :: ImobilizationRateNH4
        real    :: DenitrificationRate
        real    :: ImobilizationRateNO3

        real    :: LOM_C_Srate,LCN
        integer :: LOM_CI

        real    :: RCN

        real    :: HeteroDeathRate, HetCN     
        integer :: HetI                
        
        real    :: AnaDeathRate, AnaCN        
        integer :: AnaI                  

        real    :: AutoDeathRate, AutoCN 
        integer :: AutoI               

        real    :: Partition, AnaPartition
        real    :: Conversion, DTDay
        !------------------------------------------------------------------------
        LOM_C_Srate         = Me%SpecificRates%Labil_OM_C%Value
        LCN                 = Me%LabiOM_CN_Ratio
        LOM_CI              = Me%PropIndex%Labil_OM_C
            
        RCN                 = Me%RefractOM_CN_Ratio   
                
        DenitrificationRate     = Me%SpecificRates%NitrateToNgas%Value
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value
                      
        NII                 = Me%PropIndex%Nitrate
        AMI                 = Me%PropIndex%Ammonia

        Partition           = Me%Partition
        AnaPartition        = Me%AnaerobicPartition

        HeteroDeathRate     = Me%SpecificRates%Heterotrophs%Value
        HetI                = Me%PropIndex%HetrotrophicC
        HetCN               = Me%Microorganisms%Heterotrophs%CNRatio

        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaI                = Me%PropIndex%AnaerobicC
        AnaCN               = Me%Microorganisms%Anaerobic%CNRatio

        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutoI               = Me%PropIndex%AutotrotrophicC
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)

        DTDay               = Me%DTDay

        Me%Matrix(LOM_CI, LOM_CI)  = 1.

        !Sinks Heterotrophs uptake           
        if (Me%NLimitation)  then 
        !The OM C uptake rate is controled by the potential immobilization rate
            
            Me%Matrix(LOM_CI, AMI)  = Me%Matrix(LOM_CI, AMI) - DTDay                & 
                                    * ImobilizationRateNH4  * Conversion            &
                                    * (1./( (1./HetCN - 1./LCN ) + Partition * (1./HetCN - 1./RCN))) 
                
            Me%Matrix(LOM_CI, NII)  = Me%Matrix(LOM_CI, NII) - DTDay                &
                                    * ImobilizationRateNO3  * Conversion            & 
                                    * (1./( (1./HetCN - 1./LCN ) + Partition * (1./HetCN - 1./RCN))) 

        else
        
            Me%Matrix(LOM_CI, LOM_CI)  = Me%Matrix(LOM_CI, LOM_CI) - DTDay * LOM_C_Srate  
        end if
        
        !Sinks: Anaerobic uptake 
        Me%Matrix(LOM_CI, NII)      = Me%Matrix(LOM_CI, NII) - DTDay                & 
                                    * DenitrificationRate * Conversion              &
                                    * ( 0.1/( AnaPartition +1.) )
        
        !Sources: Anaerobic death
        if (.NOT. Me%Microorganisms%Anaerobic%LogicalMinumumPOP)                    &
            Me%Matrix(LOM_CI, AnaI) = Me%Matrix(LOM_CI, AnaI) + DTDay * AnaDeathRate
                                              
        !Sources: Hetrotrophic death
        if (.NOT. Me%Microorganisms%Heterotrophs%LogicalMinumumPOP )                &
            Me%Matrix(LOM_CI, HetI) = Me%Matrix(LOM_CI, HetI) + DTDay * HeteroDeathRate
                                              
        !Sources: Autotrotrophic death
        if (.NOT. Me%Microorganisms%Autotrophs%LogicalMinumumPOP)                   &
            Me%Matrix(LOM_CI, AutoI)= Me%Matrix(LOM_CI, AutoI) + DTDay * AutoDeathRate
        
        !Independent term
        Me%IndTerm(LOM_CI) = Me%ExternalVar%Mass(LOM_CI, index)                                    
    !------------------------------------------------------------------------

    end subroutine LabilOrganicCarbon 
    !----------------------------------------------------------------------------


    !Refractary Organic Carbon
    !
    !SOURCES: - Nitrification eficiency 
    !SINKS:   - Death
    !----------------------------------------------------------------------------    
    subroutine RefractOrganicCarbon (index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN)         :: index

       !Local-------------------------------------------------------------------
        integer :: AMI, NII
        
        real    :: ImobilizationRateNH4
        real    :: DenitrificationRate
        real    :: ImobilizationRateNO3

        real    :: ROM_C_Srate, RCN
        integer :: ROM_CI

        real    :: LCN

        real    :: HeteroDeathRate, HetCN     
        integer :: HetI                
        
        real    :: AnaDeathRate, AnaCN  
        integer :: AnaI                

        real    :: AutoDeathRate, AutoCN 
        integer :: AutoI               

        real    :: Partition, AnaPartition

        real    :: Conversion, DTDay
        !------------------------------------------------------------------------
        
        ROM_C_Srate         = Me%SpecificRates%RefractOM_C%Value
        RCN                 = Me%RefractOM_CN_Ratio   
        ROM_CI              = Me%PropIndex%RefractOM_C
            
        LCN                 = Me%LabiOM_CN_Ratio   
                
        DenitrificationRate     = Me%SpecificRates%NitrateToNgas%Value
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value
                      
        NII                 = Me%PropIndex%Nitrate
        AMI                 = Me%PropIndex%Ammonia

        Partition           = Me%Partition
        AnaPartition        = Me%AnaerobicPartition

        HeteroDeathRate     = Me%SpecificRates%Heterotrophs%Value
        HetI                = Me%PropIndex%HetrotrophicC
        HetCN               = Me%Microorganisms%Heterotrophs%CNRatio

        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaI                = Me%PropIndex%AnaerobicC
        AnaCN               = Me%Microorganisms%Anaerobic%CNRatio

        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutoI               = Me%PropIndex%AutotrotrophicC
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)
        
        DTDay               = Me%DTDay

        
        Me%Matrix(ROM_CI, ROM_CI)   = 1.

        !Sinks Heterotrophs uptake           
        if (Me%NLimitation)  then 
        !The OM N uptake rate is controlrd by the potential immobilization rate
            
            Me%Matrix(ROM_CI, AMI)   = Me%Matrix(ROM_CI, AMI) - DTDay               & 
                                     * ImobilizationRateNH4  * Conversion           &
                                     * (1./( (1./HetCN - 1./LCN )                   &
                                     * 1./ Partition + (1./HetCN - 1./RCN)))       
                
            Me%Matrix(ROM_CI, NII)   = Me%Matrix(ROM_CI, NII) - DTDay               & 
                                     * ImobilizationRateNO3  * Conversion           &
                                     * (1./( (1./HetCN - 1./LCN )                   &
                                     * 1./ Partition +  (1./HetCN - 1./RCN)))     

        else
        
            Me%Matrix(ROM_CI,ROM_CI) = Me%Matrix(ROM_CI, ROM_CI) - DTDay * ROM_C_Srate
        end if
        
        !Sinks: Anaerobic uptake 
        Me%Matrix(ROM_CI, NII)       = Me%Matrix(ROM_CI, NII) - DTDay               &
                                     * DenitrificationRate * Conversion             &
                                     * ( 0.1/( 1./ AnaPartition +1.) )             
                                                                                 
        !Independent term
        Me%IndTerm(ROM_CI)       =  Me%ExternalVar%Mass(ROM_CI, index) 
        !------------------------------------------------------------------------

    end subroutine RefractOrganicCarbon 
    !----------------------------------------------------------------------------


    !Heterotrophic Carbon
    !
    !SOURCES: - Nitrification eficiency 
    !SINKS:   - Death
    !----------------------------------------------------------------------------    
    subroutine HetrotrophicC (index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN)         :: index

        !Local-------------------------------------------------------------------
        real    :: HeteroDeathRate, HetCN, HetEF
        integer :: HetCI
        
        real    :: ROM_C_Srate, RCN                 
        integer :: ROM_CI                          

        real    :: LOM_C_Srate, LCN                 
        integer :: LOM_CI             
       
        real    :: ImobilizationRateNO3
        
        integer :: NII, AMI
        
        real    :: ImobilizationRateNH4
        
        real    :: Partition, Conversion, DTDay
        !------------------------------------------------------------------------

        HeteroDeathRate      = Me%SpecificRates%Heterotrophs%Value
        HetCI                = Me%PropIndex%HetrotrophicC
        HetCN                = Me%Microorganisms%Heterotrophs%CNRatio
        HetEF                = Me%Microorganisms%Heterotrophs%EficiencyC

        ROM_C_Srate          = Me%SpecificRates%RefractOM_C%Value
        RCN                  = Me%RefractOM_CN_Ratio   
        ROM_CI               = Me%PropIndex%RefractOM_C
            
        LOM_C_Srate          = Me%SpecificRates%Labil_OM_C%Value
        LCN                  = Me%LabiOM_CN_Ratio
        LOM_CI               = Me%PropIndex%Labil_OM_C
        
        NII                  = Me%PropIndex%Nitrate
        ImobilizationRateNO3 = Me%SpecificRates%NitrateImobilization%Value

        AMI                  = Me%PropIndex%Ammonia
        ImobilizationRateNH4 = Me%SpecificRates%AmmoniaImobilization%Value

        Partition            = Me%Partition
        
        Conversion           = Me%ExternalVar%DissolvedToParticulate (index)

        DTDay                = Me%DTDay

        
        Me%Matrix(HetCI, HetCI) = 1.
        
        !Sink   : Hetrotrophic death
        if (.NOT. Me%Microorganisms%Heterotrophs%LogicalMinumumPOP)                 &
           Me%Matrix(HetCI, HetCI) = Me%Matrix(HetCI, HetCI) - DTDay * HeteroDeathRate

         if (Me%NLimitation)  then 
         !The OM N uptake rate is controlrd by the potential immobilization rate
            
        !Source : Labil and Refract C uptake
            Me%Matrix(HetCI, AMI)  = Me%Matrix(HetCI, AMI) + DTDay                  & 
                                   * ImobilizationRateNH4  * Conversion             &
                                   * (1./((1./HetCN - 1./LCN)                       &
                                   + Partition * (1./HetCN - 1./RCN))               &
                                   + 1./((1./HetCN - 1./LCN ) * 1./Partition        &
                                   + (1./HetCN - 1./RCN)))
                
            Me%Matrix(HetCI, NII)  = Me%Matrix(HetCI, NII) + DTDay                  &
                                   * ImobilizationRateNO3  * Conversion             &
                                   * (1./((1./HetCN - 1./LCN)                       &
                                   + Partition * (1./HetCN - 1./RCN))               &
                                   + 1./((1./HetCN - 1./LCN ) * 1./ Partition       &
                                   + (1./HetCN - 1./RCN)))               
                                                   
        !Sink   : Labil And refract C breathe
            Me%Matrix(HetCI, AMI)  = Me%Matrix(HetCI, AMI) - DTDay                  &
                                   * ImobilizationRateNH4 * Conversion * (1.- HetEF)&
                                   * (1./( (1./HetCN - 1./LCN)                      &
                                   + Partition * (1./HetCN - 1./RCN))               &
                                   + 1./((1./HetCN - 1./LCN ) * 1./ Partition       &
                                   + (1./HetCN - 1./RCN)))  
                
            Me%Matrix(HetCI, NII)  = Me%Matrix(HetCI, NII) - DTDay                  & 
                                   * ImobilizationRateNO3 * Conversion * (1.- HetEF)& 
                                   * (1./((1./HetCN - 1./LCN )                      &
                                   + Partition * (1./HetCN - 1./RCN))               &
                                   + 1./((1./HetCN - 1./LCN ) * 1./ Partition       &
                                   + (1./HetCN - 1./RCN)))

        else ! OM C uptake equals the potential rate

            !Source : Refract C uptake
            Me%Matrix(HetCI, ROM_CI) = Me%Matrix(HetCI, ROM_CI) + DTDay * ROM_C_Srate
            
            !Sink   : Refract C breathe
            Me%Matrix(HetCI, ROM_CI) = Me%Matrix(HetCI, ROM_CI) - DTDay * ROM_C_Srate * (1.- HetEF)          
            
            !Source : Labil C uptake
            Me%Matrix(HetCI, LOM_CI) = Me%Matrix(HetCI, LOM_CI) + DTDay * LOM_C_Srate
            
            !Sink   : Labil C breathe
            Me%Matrix(HetCI, LOM_CI) = Me%Matrix(HetCI, LOM_CI) - DTDay * LOM_C_Srate * (1.- HetEF)
            
        end if
            
            !Independent term
            Me%IndTerm(HetCI)       = Me%ExternalVar%Mass(HetCI, index) 
    !------------------------------------------------------------------------

    end subroutine HetrotrophicC 
    !----------------------------------------------------------------------------


    !Autotrotrophic Carbon
    !
    !SOURCES: - Nitrification eficiency 
    !SINKS:   - Death
    !----------------------------------------------------------------------------    
    subroutine AutotrotrophicC (index)

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

        AutoCI              = Me%PropIndex%AutotrotrophicC
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
                               * NitificationRate * Conversion * AutoEf * AutoCN
        
        !Independent term
        Me%IndTerm(AutoCI)     = Me%ExternalVar%Mass(AutoCI, index) 
    !------------------------------------------------------------------------

    end subroutine AutotrotrophicC 
    !----------------------------------------------------------------------------


    !Anaerobic Carbon
    !
    !SOURCES: - Nitrification eficiency 
    !SINKS:   - Death, respiration
    !----------------------------------------------------------------------------    
    subroutine AnaerobicC (index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN)         :: index
        !------------------------------------------------------------------------

        real    :: AnaCN, AnaCEf, AnaDeathRate
        integer :: AnaCI
        real    :: DenitrificationRate, Conversion 
        
        integer :: NII, DTDay    
        !------------------------------------------------------------------------
        
        DTDay               = Me%DTDay

        AnaCI               = Me%PropIndex%AnaerobicC
        NII                 = Me%PropIndex%Nitrate
        
        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaCN               = Me%Microorganisms%Anaerobic%CNRatio
        
        DenitrificationRate = Me%SpecificRates%NitrateToNgas%Value
        AnaCEf              = Me%Microorganisms%Anaerobic%EficiencyC

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)        

        Me%Matrix(AnaCI, AnaCI) = 1.

        !Sink   : Death
        if (.NOT.Me%Microorganisms%Anaerobic%LogicalMinumumPOP )                    &
            Me%Matrix(AnaCI, AnaCI) = Me%Matrix(AnaCI, AnaCI) - DTDay * AnaDeathRate
        
        !Sources: Labil and Refract C
        Me%Matrix(AnaCI, NII) = Me%Matrix(AnaCI, NII) + DTDay * DenitrificationRate &
                              * Conversion * 0.1
        !Sink   : Respiration 
        Me%Matrix(AnaCI, NII) = Me%Matrix(AnaCI, NII) - DTDay * DenitrificationRate &
                              * Conversion * 0.1 * AnaCEf
        !Independent term
        Me%IndTerm(AnaCI)     = Me%ExternalVar%Mass(AnaCI, index) 
    !------------------------------------------------------------------------

    end subroutine AnaerobicC 
    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------
    subroutine Nitrogen(index)

        !Arguments---------------------------------------------------------------
        integer, intent(IN) :: index
        !------------------------------------------------------------------------

        call Ammonia                (index)

        call Nitrate                (index)
                                                        
        call LabilOrganicNitrogen   (index)

        call RefractOrganicNitrogen (index)

        call HetrotrophicN          (index)

        call AutotrotrophicN        (index)

        call AnaerobicN             (index)
        
        call Ngas                   (index)
    !------------------------------------------------------------------------

    end subroutine Nitrogen
    !---------------------------------------------------------------------------



    !AMMONIA 
    !
    !SOURCES    : Heterophs
    !SINKS      : Heterophs, NO3 (Autotrophs), N gas (Anaerobic), Plants (to_do)                        
    !----------------------------------------------------------------------------       
    subroutine Ammonia(index)
    
        !Arguments---------------------------------------------------------------
        integer, intent(IN)         :: index            
        
        !Local-------------------------------------------------------------------
        integer :: AMI, NII         !Ammonia and Nitrate

        integer :: LOM_CI
        real    :: LOM_C_Srate, LCN, CurrentLabilC, potLCuptake

        integer :: ROM_CI 
        real    :: ROM_C_Srate, RCN, CurrentRefractC, potRCuptake
        
        real    :: neededNuptake
        
        real    :: HetCN, HetEF, AnaCN, AnaNEF
        
        real    :: NitrificationRate, DenitrificationRate
        real    :: ImobilizationRateNH4, ImobilizationRateNO3
                
        real    :: Partition, AnaPartition

        real    :: Conversion, DTDay, totalNuptake        
        !------------------------------------------------------------------------
        
        AMI             = Me%PropIndex%Ammonia
        NII             = Me%PropIndex%Nitrate
                
        HetCN           = Me%Microorganisms%Heterotrophs%CNRatio
        HetEF           = Me%Microorganisms%Heterotrophs%EficiencyC
        AnaCN           = Me%Microorganisms%Anaerobic%CNRatio
        AnaNEF          = Me%Microorganisms%Anaerobic%EficiencyN       
        
        LOM_CI          = Me%PropIndex%Labil_OM_C
        LOM_C_Srate     = Me%SpecificRates%Labil_OM_C%Value
        LCN             = Me%LabiOM_CN_Ratio
        CurrentLabilC   = Me%ExternalVar%Mass(LOM_CI, index)

        ROM_CI          = Me%PropIndex%RefractOM_C        
        ROM_C_Srate     = Me%SpecificRates%RefractOM_C%Value
        RCN             = Me%RefractOM_CN_Ratio
        CurrentRefractC = Me%ExternalVar%Mass(ROM_CI, index)
                
        Partition       = Me%Partition
        AnaPartition    = Me%AnaerobicPartition
        Conversion      = Me%ExternalVar%DissolvedToParticulate (index)        

        DTDay           = Me%DTDay
        !------------------------------------------------------------------------

        Me%Matrix(AMI, AMI) = 1.   
        !Sinks : NO3 formation
        NitrificationRate   = Me%SpecificRates%AmmoniaToNitrate%Value        
        Me%Matrix(AMI, AMI) = Me%Matrix(AMI, AMI) - DTDay * NitrificationRate 
        
        !Sources : Anaerobic excretion
        DenitrificationRate = Me%SpecificRates%NitrateToNgas%Value   
        
        Me%Matrix(AMI, NII) = Me%Matrix(AMI, NII) + DTDay * DenitrificationRate     &
                            * (AnaNEF + 0.1 / ( (AnaPartition + 1. ) * LCN )        &
                            + 0.1 / (( 1./AnaPartition + 1. ) * RCN) - 0.1 / AnaCN)


        if (.NOT. Me%Imobilization .OR. .NOT. Me%NLimitation) then       
        
        !Sources : Hetrotrophs excretion
            Me%Matrix(AMI, LOM_CI) = Me%Matrix(AMI, LOM_CI) + DTDay                 &
                                   * (LOM_C_Srate * 1 / Conversion * HetEF/HetCN) 

            Me%Matrix(AMI, ROM_CI) = Me%Matrix(AMI, ROM_CI) + DTDay                 &
                                   * (ROM_C_Srate * 1 / Conversion * HetEF/HetCN)
        end if
              
                                                   
        if (Me%Imobilization) then !immobilization

            ImobilizationRateNH4   = Me%SpecificRates%AmmoniaImobilization%Value
            ImobilizationRateNO3   = Me%SpecificRates%NitrateImobilization%Value
                   
            if (Me%NLimitation)  then    !The N immobilization rate is the potential rate
        
        !Sinks : NH4 imobilization

               Me%Matrix(AMI, AMI) = Me%Matrix(AMI, AMI) - DTDay * ImobilizationRateNH4
        
        !Sources : Hetrotrophs labil and refract N excretion
                   
               Me%Matrix(AMI, AMI) = Me%Matrix(AMI, AMI) + DTDay                    &
                                   * ImobilizationRateNH4 * HetEF/HetCN             &
                                   * (1./((1./HetCN - 1./LCN)                       &
                                   + Partition * (1./HetCN - 1./RCN))               &
                                   + 1./((1./HetCN - 1./LCN) * 1. / Partition       &
                                   + (1./HetCN - 1./RCN))) 
                
               Me%Matrix(AMI, NII) = Me%Matrix(AMI, NII) + DTDay                    &
                                   * ImobilizationRateNO3 * HetEF/HetCN             &
                                   * (1./((1./HetCN - 1./LCN )                      &
                                   + Partition * (1./HetCN - 1./RCN))               &
                                   + 1./((1./HetCN - 1./LCN) * 1. / Partition       &
                                   + (1./HetCN - 1./RCN)))
           
            else 
         
                !Sinks : NH4 imobilization                                    
               Me%Matrix(AMI, LOM_CI) = Me%Matrix(AMI, LOM_CI) - DTDay              &
                                      * (LOM_C_Srate * 1. / Conversion              &
                                      * (1./HetCN - 1./LCN) ) / (Partition + 1.)    
                
               Me%Matrix(AMI, ROM_CI) = Me%Matrix(AMI, ROM_CI) - DTDay              &
                                      * (ROM_C_Srate * 1. / Conversion              &
                                      * (1./HetCN - 1./RCN) ) / ( Partition + 1.)
        
            end if
        
        end if
        

        !Independent term
        Me%IndTerm(AMI) = Me%ExternalVar%Mass(AMI, index)

        !Exessive N in the OM
            if (.NOT. Me%Imobilization)      then

                potLCuptake     = LOM_C_Srate * Me%ExternalVar%Mass(LOM_CI, index) 
                potRCuptake     = ROM_C_Srate * Me%ExternalVar%Mass(ROM_CI, index)  
                
                totalNuptake    = potLCuptake /LCN + potRCuptake / RCN 
                neededNuptake   = (potLCuptake + potRCuptake) / HetCN  
                
                if  (totalNuptake > neededNuptake)  then  
                    Me%IndTerm(AMI) = Me%IndTerm(AMI)                               &
                                    + (totalNuptake - neededNuptake)                &
                                    * 1./ Conversion 
                end if 
                           
            end if
    !------------------------------------------------------------------------

    end subroutine Ammonia
    !----------------------------------------------------------------------------


    !NITRATE 
    !
    !SOURCES: - Nitrification 
    !SINKS:   - Immobilization, denitrification
    !----------------------------------------------------------------------------
    subroutine Nitrate(index)

        !Arguments---------------------------------------------------------------
        integer                 , intent(IN) :: index

        !Local-------------------------------------------------------------------
        integer :: AMI, NII
        
        real    :: NitrificationRate, DenitrificationRate, ImobilizationRateNO3

        real    :: LOM_C_Srate, LCN
        integer :: LOM_CI

        real    :: ROM_C_Srate, RCN
        integer :: ROM_CI

        real    :: HetCN, AutoEf

        real    :: Partition, Conversion, DTDay
        !------------------------------------------------------------------------

        NII     = Me%PropIndex%Nitrate
        AMI     = Me%PropIndex%Ammonia
        DTDay   = Me%DTDay

        NitrificationRate   = Me%SpecificRates%AmmoniaToNitrate%Value
        DenitrificationRate = Me%SpecificRates%NitrateToNgas%Value
        
        AutoEf              = Me%Microorganisms%Autotrophs%EficiencyN
        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)
        
        Me%Matrix(NII, NII) = 1.    
        !Sources :
        !Nitrification
        Me%Matrix(NII, AMI) = Me%Matrix(NII, AMI) + DTDay                           &
                            * NitrificationRate * (1. - AutoEf) 
        
        !Sinks:
        !Denitrification
        Me%Matrix(NII, NII) = Me%Matrix(NII, NII) - DTDay * DenitrificationRate 
        
        !Immobilization
        if (Me%Imobilization) then !immobilization

            ImobilizationRateNO3 = Me%SpecificRates%NitrateImobilization%Value
                   
            if (Me%NLimitation)  then    
            !The N immobilization rate is the potential rate
                
                Me%Matrix(NII, NII) = Me%Matrix(NII, NII) - DTDay * ImobilizationRateNO3
           
            else 

                Partition   = Me%Partition
                
                LOM_C_Srate = Me%SpecificRates%Labil_OM_C%Value
                LCN         = Me%LabiOM_CN_Ratio          
                LOM_CI      = Me%PropIndex%Labil_OM_C

                ROM_C_Srate = Me%SpecificRates%RefractOM_C%Value
                RCN         = Me%RefractOM_CN_Ratio   
                ROM_CI      = Me%PropIndex%RefractOM_C

                HetCN       = Me%Microorganisms%Heterotrophs%CNRatio
                                    
                Me%Matrix(NII, LOM_CI) = Me%Matrix(NII, LOM_CI) - DTDay             &
                                       * (LOM_C_Srate * 1. / Conversion             &
                                       * (1./HetCN - 1./LCN)) / (1./Partition + 1.)

                Me%Matrix(NII, ROM_CI) = Me%Matrix(NII, ROM_CI) - DTDay             &
                                       * (ROM_C_Srate * 1. / Conversion             &
                                       * (1./HetCN - 1./RCN)) / (1./Partition + 1.) 
        
            end if
        end if               

        !Independent term
        Me%IndTerm(NII) = Me%ExternalVar%Mass(NII, index) 
    !----------------------------------------------------------------------------

    end subroutine Nitrate
    !----------------------------------------------------------------------------

    
    !Labil Organic Nitrogen 
    !
    !SOURCES: - Microorganisms death
    !SINKS:   - Heterotrphs uptake, including anaerobic uptake 
    !----------------------------------------------------------------------------
    subroutine LabilOrganicNitrogen (index)

        !Arguments---------------------------------------------------------------
        integer                , intent(IN)         :: index

        !Local-------------------------------------------------------------------
        integer :: AMI, NII
        
        real    :: ImobilizationRateNH4, DenitrificationRate, ImobilizationRateNO3

        real    :: LOM_C_Srate, LCN
        integer :: LOM_CI, LOM_NI

        real    :: RCN

        real    :: HeteroDeathRate, HetCN     
        integer :: HetI                
        
        real    :: AnaDeathRate, AnaCN        
        integer :: AnaI                
        
        real    :: AutoDeathRate, AutoCN
        integer :: AutI               

        real    :: Partition, AnaPartition, Conversion, DTDay
        !------------------------------------------------------------------------
        
        LOM_C_Srate         = Me%SpecificRates%Labil_OM_C%Value
        LCN                 = Me%LabiOM_CN_Ratio
        LOM_CI              = Me%PropIndex%Labil_OM_C
        LOM_NI              = Me%PropIndex%Labil_OM_N
            
        RCN                 = Me%RefractOM_CN_Ratio   
                
        DenitrificationRate     = Me%SpecificRates%NitrateToNgas%Value
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value
                      
        NII                 = Me%PropIndex%Nitrate
        AMI                 = Me%PropIndex%Ammonia

        Partition           = Me%Partition
        AnaPartition        = Me%AnaerobicPartition

        HeteroDeathRate     = Me%SpecificRates%Heterotrophs%Value
        HetI                = Me%PropIndex%HetrotrophicC
        HetCN               = Me%Microorganisms%Heterotrophs%CNRatio

        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaI                = Me%PropIndex%AnaerobicC
        AnaCN               = Me%Microorganisms%Anaerobic%CNRatio

        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutI                = Me%PropIndex%AutotrotrophicC
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio

        DTDay               = Me%DTDay

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)

        Me%Matrix(LOM_NI, LOM_NI)   = 1.
        
        !Sinks Heterotrophs uptake           
        
        if (Me%NLimitation)  then !The OM N uptake rate is controlrd by the potential immobilization rate
            
            Me%Matrix(LOM_NI, AMI)    = Me%Matrix(LOM_NI, AMI) - DTDay              &
                                      * ImobilizationRateNH4 * Conversion * 1./ LCN & 
                                      * (1./((1./HetCN - 1./LCN )                   &
                                      + Partition * (1./HetCN - 1./RCN))) 
                
            Me%Matrix(LOM_NI, NII)    = Me%Matrix(LOM_NI, NII) - DTDay              &
                                      * ImobilizationRateNO3 * Conversion * 1./ LCN &
                                      * (1./((1./HetCN - 1./LCN )                   &
                                      + Partition * (1./HetCN - 1./RCN)))

        else
        
            Me%Matrix(LOM_NI, LOM_CI) = Me%Matrix(LOM_NI, LOM_CI) - DTDay * LOM_C_Srate / LCN  
        end if
        
        !Sinks: Anaerobic uptake 
        Me%Matrix(LOM_NI, NII) = Me%Matrix(LOM_NI, NII) - DTDay                     & 
                               * DenitrificationRate * Conversion * 1./ LCN         &
                               * ( 0.1/( AnaPartition +1.))     
        !Sources: Anaerobic death
        if (.NOT. Me%Microorganisms%Anaerobic%LogicalMinumumPOP)                    &        
            Me%Matrix(LOM_NI, AnaI) = Me%Matrix(LOM_NI, AnaI) + DTDay               &
                                    * AnaDeathRate * 1. / AnaCN 
                                              
        !Sources: Hetrotrophic death
        if ( .NOT. Me%Microorganisms%Heterotrophs%LogicalMinumumPOP)                &        
            Me%Matrix(LOM_NI, HetI) = Me%Matrix(LOM_NI, HetI) + DTDay               &
                                    * HeteroDeathRate * 1. / HetCN
                                              
        !Sources: Autotrotrophic death
        if ( .NOT. Me%Microorganisms%Autotrophs%LogicalMinumumPOP)                  &        
            Me%Matrix(LOM_NI, AutI) = Me%Matrix(LOM_NI, AutI) + DTDay               &
                                    * AutoDeathRate * 1. / AutoCN                                       
       
        !Independent term
        Me%IndTerm(LOM_NI) = Me%ExternalVar%Mass(LOM_NI, index)                     
    !----------------------------------------------------------------------------

    end subroutine LabilOrganicNitrogen 
    !----------------------------------------------------------------------------


    !Refractary Organic Nitrogen
    !
    !SOURCES: - INTERPOOL TRANSFORMATION
    !SINKS:   - Heterotrphs uptake, including anaerobic uptake 
    !----------------------------------------------------------------------------
    subroutine RefractOrganicNitrogen (index)

        !Arguments---------------------------------------------------------------
        integer                , intent(IN) :: index

        !Local-------------------------------------------------------------------
        integer :: AMI, NII
        
        real    :: ImobilizationRateNH4, DenitrificationRate, ImobilizationRateNO3

        real    :: ROM_C_Srate, RCN
        integer :: ROM_CI, ROM_NI

        real    :: LCN

        real    :: HeteroDeathRate, HetCN   
        integer :: HetI                
        
        real    :: AnaDeathRate, AnaCN        
        integer :: AnaI                

        real    :: AutoDeathRate, AutoCN
        integer :: AutoI               

        real    :: Partition, AnaPartition, Conversion, DTDay
        !------------------------------------------------------------------------
        ROM_C_Srate         = Me%SpecificRates%RefractOM_C%Value
        RCN                 = Me%RefractOM_CN_Ratio   
        ROM_CI              = Me%PropIndex%RefractOM_C
        ROM_NI              = Me%PropIndex%RefractOM_N
            
        LCN                 = Me%LabiOM_CN_Ratio   
                
        DenitrificationRate     = Me%SpecificRates%NitrateToNgas%Value
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value
                      
        NII                 = Me%PropIndex%Nitrate
        AMI                 = Me%PropIndex%Ammonia

        Partition           = Me%Partition
        AnaPartition        = Me%AnaerobicPartition

        HeteroDeathRate     = Me%SpecificRates%Heterotrophs%Value
        HetI                = Me%PropIndex%HetrotrophicC
        HetCN               = Me%Microorganisms%Heterotrophs%CNRatio

        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaI                = Me%PropIndex%AnaerobicC
        AnaCN               = Me%Microorganisms%Anaerobic%CNRatio

        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutoI               = Me%PropIndex%AutotrotrophicC
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio

        DTDay               = Me%DTDay
        
        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)
        
        Me%Matrix(ROM_NI, ROM_NI)  = 1.
        
        !Sinks Heterotrophs uptake           
        if (Me%NLimitation)  then 
        !The OM N uptake rate is controlrd by the potential immobilization rate
            
            Me%Matrix(ROM_NI, AMI)  = Me%Matrix(ROM_NI, AMI) - DTDay                &
                                    * ImobilizationRateNH4 * Conversion * 1./RCN    &
                                    * (1./((1./HetCN - 1./LCN ) * 1./ Partition     &
                                    + (1./HetCN - 1./RCN))) 
                
            Me%Matrix(ROM_NI, NII)  = Me%Matrix(ROM_NI, NII) - DTDay                &
                                    * ImobilizationRateNO3 * Conversion * 1./RCN    &
                                    * (1./((1./HetCN - 1./LCN )* 1./ Partition      &
                                    +  (1./HetCN - 1./RCN)))

        else
        
            Me%Matrix(ROM_NI, ROM_CI) = Me%Matrix(ROM_NI, ROM_CI) - DTDay * ROM_C_Srate / RCN  
        end if
        
        !Sinks: Anaerobic uptake 
        Me%Matrix(ROM_NI, NII) = Me%Matrix(ROM_NI, NII) - DTDay                     &
                               * DenitrificationRate * Conversion * 1./RCN          &
                               * (0.1/( 1./AnaPartition +1.) )     
                                             
        !Independent term
        Me%IndTerm(ROM_NI)     = Me%ExternalVar%Mass(ROM_NI, index)                                     
    !----------------------------------------------------------------------------

    end subroutine RefractOrganicNitrogen 
    !----------------------------------------------------------------------------

    
    !Heterotrophic N
    !
    !SOURCES: - Organic matter N decay
    !SINKS:   - Hetrotrophic Death, excretion
    !----------------------------------------------------------------------------
    subroutine HetrotrophicN (index)

        !Arguments---------------------------------------------------------------
        integer                 , intent(IN) :: index

        !Local-------------------------------------------------------------------
        real    :: HeteroDeathRate, HetCN, HetEF
        integer :: HetCI, HetNI           
        
        real    :: ROM_C_Srate, RCN, potRCuptake
        integer :: ROM_CI, ROM_NI              

        real    :: LOM_C_Srate, LCN, potLCuptake          
        integer :: LOM_CI, LOM_NI

        real    :: totalNuptake, neededNuptake
        
        real    :: ImobilizationRateNO3, ImobilizationRateNH4
        integer :: NII, AMI
        
        real    :: Partition, Conversion, DTDay
        !------------------------------------------------------------------------

        HeteroDeathRate = Me%SpecificRates%Heterotrophs%Value
        HetCI           = Me%PropIndex%HetrotrophicC
        HetNI           = Me%PropIndex%HetrotrophicN
        HetCN           = Me%Microorganisms%Heterotrophs%CNRatio
        HetEF           = Me%Microorganisms%Heterotrophs%EficiencyC

        ROM_C_Srate         = Me%SpecificRates%RefractOM_C%Value
        RCN                 = Me%RefractOM_CN_Ratio   
        ROM_CI              = Me%PropIndex%RefractOM_C
        ROM_NI              = Me%PropIndex%RefractOM_N
            
        LOM_C_Srate         = Me%SpecificRates%Labil_OM_C%Value
        LCN                 = Me%LabiOM_CN_Ratio
        LOM_CI              = Me%PropIndex%Labil_OM_C
        LOM_NI              = Me%PropIndex%Labil_OM_N
        
        NII                     = Me%PropIndex%Nitrate
        ImobilizationRateNO3    = Me%SpecificRates%NitrateImobilization%Value

        AMI                     = Me%PropIndex%Ammonia
        ImobilizationRateNH4    = Me%SpecificRates%AmmoniaImobilization%Value

        Partition           = Me%Partition
        
        DTDay               = Me%DTDay

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)

        
        Me%Matrix(HetNI, HetNI) = 1.
        
        !Sink: Hetrotrophic death
        if ( .not. Me%Microorganisms%Heterotrophs%LogicalMinumumPOP)                & 
            Me%Matrix(HetNI, HetCI) = Me%Matrix(HetNI, HetCI) - DTDay * HeteroDeathRate * 1. / HetCN

        !Sink Heterotrophs excretion
        if (.NOT. Me%NLimitation) then       
                                                
                Me%Matrix(HetNI, LOM_CI) = Me%Matrix(HetNI, LOM_CI) - DTDay         &
                                         * ( LOM_C_Srate * HetEF/HetCN ) 

                Me%Matrix(HetNI, ROM_CI) = Me%Matrix(HetNI, ROM_CI) - DTDay         &
                                         * ( ROM_C_Srate * HetEF/HetCN )
        else
            
                Me%Matrix(HetNI, AMI)    = Me%Matrix(HetNI, AMI) - DTDay            &
                                         * ImobilizationRateNH4 * Conversion * HetEF/HetCN &
                                         * (1./((1./HetCN - 1./LCN)                 &
                                         + Partition * (1./HetCN - 1./RCN))         &
                                         + 1./((1./HetCN - 1./LCN) * 1. / Partition &
                                         + (1./HetCN - 1./RCN))) 
                
                Me%Matrix(HetNI, NII)    = Me%Matrix(HetNI, NII) - DTDay            &
                                         * ImobilizationRateNO3 * Conversion * HetEF/HetCN &
                                         * (1./((1./HetCN - 1./LCN )                &
                                         + Partition * (1./HetCN - 1./RCN))         &
                                         + 1./((1./HetCN - 1./LCN ) * 1./ Partition & 
                                         + (1./HetCN - 1./RCN))) 
           
        end if
        
        !Source: OM N uptake
         if (Me%NLimitation)  then 
         !The OM N uptake rate is controlrd by the potential immobilization rate
            
            !Labil N uptake
            Me%Matrix(HetNI, AMI) = Me%Matrix(HetNI, AMI) + DTDay                   &
                                  * ImobilizationRateNH4 * Conversion * 1./LCN      &
                                  * (1./((1./HetCN - 1./LCN)                        &
                                  + Partition * (1./HetCN - 1./RCN))) 
                
            Me%Matrix(HetNI, NII) = Me%Matrix(HetNI, NII) + DTDay                   &
                                  * ImobilizationRateNO3 * Conversion * 1./LCN      &
                                  *(1./((1./HetCN - 1./LCN)                         &
                                  + Partition * (1./HetCN - 1./RCN))) 
            !Refract N uptake
            Me%Matrix(HetNI, AMI) = Me%Matrix(HetNI, AMI) + DTDay                   &
                                  * ImobilizationRateNH4 * Conversion * 1./RCN      &
                                  * (1./((1./HetCN - 1./LCN) * 1./ Partition        &
                                  + (1./HetCN - 1./RCN))) 
                
            Me%Matrix(HetNI, NII) = Me%Matrix(HetNI, NII) + DTDay                   &
                                  * ImobilizationRateNO3 * Conversion * 1./RCN      &
                                  * (1./((1./HetCN - 1./LCN ) * 1./ Partition       &
                                  + (1./HetCN - 1./RCN))) 
            !NH4 immobilization
            Me%Matrix(HetNI, AMI) = Me%Matrix(HetNI, AMI) + DTDay                   &
                                  * ImobilizationRateNH4 * Conversion 

            !NO3 immobilization
            Me%Matrix(HetNI, NII) = Me%Matrix(HetNI, NII) + DTDay                   &
                                  * ImobilizationRateNO3 * Conversion 

        else
            !Labil N uptake    
            Me%Matrix(HetNI, ROM_CI) = Me%Matrix(HetNI, ROM_CI) + DTDay * ROM_C_Srate / RCN 
            !Refract N uptake
            Me%Matrix(HetNI, LOM_CI) = Me%Matrix(HetNI, LOM_CI) + DTDay * LOM_C_Srate / LCN
                
                if (Me%Imobilization) then
                    !NO3 immobilization
                    Me%Matrix(HetNI, LOM_CI) = Me%Matrix(HetNI, LOM_CI) + DTDay     &
                                             * (LOM_C_Srate * (1./HetCN - 1./LCN))/(1./Partition+1.)
                
                    Me%Matrix(HetNI, ROM_CI) = Me%Matrix(HetNI, ROM_CI) + DTDay     &
                                             * (ROM_C_Srate * (1./HetCN - 1./RCN))/(1./Partition+1.)


                    !NH4 immobilization
                    Me%Matrix(HetNI, LOM_CI) = Me%Matrix(HetNI, LOM_CI) + DTDay     &
                                             * (LOM_C_Srate * (1./HetCN - 1./LCN))/(Partition+1.)
                
                    Me%Matrix(HetNI, ROM_CI) = Me%Matrix(HetNI, ROM_CI) + DTDay     &
                                             * (ROM_C_Srate * (1./HetCN - 1./RCN))/(Partition+1.)
                end if                                                                         
        end if
            
            !Independent term
            Me%IndTerm(HetNI) = Me%ExternalVar%Mass(HetNI,index) 
            
            !Exessive N in the OM
            if (.NOT. Me%Imobilization)      then

                potLCuptake     = LOM_C_Srate * Me%ExternalVar%Mass(LOM_CI, index) 
                potRCuptake     = ROM_C_Srate * Me%ExternalVar%Mass(ROM_CI, index)  
                
                totalNuptake    = potLCuptake /LCN + potRCuptake / RCN 
                neededNuptake   = (potLCuptake + potRCuptake) / HetCN  
                
                if  (totalNuptake > neededNuptake)  then  
                    Me%IndTerm(HetNI) = Me%IndTerm(HetNI) - (totalNuptake - neededNuptake ) 
                end if 
                           
            end if
    !----------------------------------------------------------------------------

    end subroutine HetrotrophicN 
    !----------------------------------------------------------------------------

    
    !Anaerobic N
    !
    !SOURCES: - Denitrification eficiency
    !SINKS:   - Death
    !----------------------------------------------------------------------------    
    subroutine AnaerobicN (index)

        !Arguments---------------------------------------------------------------
        integer                , intent(IN) :: index

        !Local-------------------------------------------------------------------
        real    :: AnaCN, AnaNEf, AnaCEf, AnaDeathRate
        integer :: AnaNI, AnaCI

        real    :: DenitrificationRate 

        real    :: LCN, RCN

        real    :: AnaPartition, Conversion, DTDay
        
        integer :: NII    
        !------------------------------------------------------------------------
        
        DTDay               = Me%DTDay

        AnaNI               = Me%PropIndex%AnaerobicN
        AnaCI               = Me%PropIndex%AnaerobicC
        NII                 = Me%PropIndex%Nitrate
        
        AnaDeathRate        = Me%SpecificRates%Anaerobic%Value
        AnaCN               = Me%Microorganisms%Anaerobic%CNRatio
        
        DenitrificationRate  = Me%SpecificRates%NitrateToNgas%Value
        AnaNEf              = Me%Microorganisms%Anaerobic%EficiencyN
        AnaCEf              = Me%Microorganisms%Anaerobic%EficiencyC

        LCN                 = Me%LabiOM_CN_Ratio
        RCN                 = Me%RefractOM_CN_Ratio
    
        AnaPartition        = Me%AnaerobicPartition

        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)

        Me%Matrix(AnaNI, AnaNI) = 1.
        !Sink: Death
        if ( .not. Me%Microorganisms%Anaerobic%LogicalMinumumPOP )                  &        
            Me%Matrix(AnaNI, AnaCI) = Me%Matrix(AnaNI, AnaCI) - DTDay * AnaDeathRate / AnaCN  
        
        !Sources: Nitrification (direct assimilation)
        Me%Matrix(AnaNI, NII) = Me%Matrix(AnaNI, NII) + DTDay                       &
                              * DenitrificationRate * Conversion * AnaNEf
        !Sources: Labil N
        Me%Matrix(AnaNI, NII) = Me%Matrix(AnaNI, NII) + DTDay                       &
                              * DenitrificationRate * Conversion                    &
                              * (0.1 /( (AnaPartition + 1. ) * LCN)) 
        !Sources: Refract N
        Me%Matrix(AnaNI, NII) = Me%Matrix(AnaNI, NII) + DTDay                       &
                              * DenitrificationRate * Conversion                    &
                              * (0.1 /( (1./AnaPartition + 1. ) * RCN )) 
        !Sink: Excretion
        Me%Matrix(AnaNI, NII) = Me%Matrix(AnaNI, NII) - DTDay                       &
                              * DenitrificationRate * Conversion                    &
                              * (AnaNEf + 0.1/((AnaPartition + 1.) * LCN )          &
                              + 0.1/(( 1./AnaPartition + 1.)* RCN)                  &
                              - 0.1 * (1 - AnaCEf)/ AnaCN)
        
        !Independent term
        Me%IndTerm(AnaNI)    = Me%ExternalVar%Mass(AnaNI, index) 
    !----------------------------------------------------------------------------

    end subroutine AnaerobicN 
    !----------------------------------------------------------------------------


    
    !Autotrotrophic N
    !
    !SOURCES: - Nitrification eficiency 
    !SINKS:   - Death
    !----------------------------------------------------------------------------    
    subroutine AutotrotrophicN (index)

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

        AutoNI              = Me%PropIndex%AutotrotrophicN
        AutoCI              = Me%PropIndex%AutotrotrophicC
        AMI                 = Me%PropIndex%Ammonia
        
        AutoDeathRate       = Me%SpecificRates%Autotrophs%Value
        AutoCN              = Me%Microorganisms%Autotrophs%CNRatio
        
        NitificationRate    = Me%SpecificRates%AmmoniaToNitrate%Value
        AutoEf              = Me%Microorganisms%Autotrophs%EficiencyN
        
        Conversion          = Me%ExternalVar%DissolvedToParticulate (index)

        Me%Matrix(AutoNI, AutoNI) = 1.

        !Sink: Death
        if ( .not. Me%Microorganisms%Autotrophs%LogicalMinumumPOP)                  &        
            Me%Matrix(AutoNI, AutoCI) = Me%Matrix(AutoNI, AutoCI) - DTDay * AutoDeathRate / AutoCN  
        
        !Sources: Nitrification
        Me%Matrix(AutoNI, AMI) = Me%Matrix(AutoNI, AMI) + DTDay * NitificationRate * Conversion * AutoEf
        
        !Independent term
        Me%IndTerm(AutoNI) = Me%ExternalVar%Mass(AutoNI, index) 
    !----------------------------------------------------------------------------

    end subroutine AutotrotrophicN 
    !----------------------------------------------------------------------------

      
    !N gas (N20 N2)
    !
    !SOURCES: - Denitrification
    !SINKS:   - none for the moment
    !----------------------------------------------------------------------------
    subroutine Ngas (index)

        !Arguments---------------------------------------------------------------
        integer                , intent(IN) :: index
    
        !Local-------------------------------------------------------------------
        real    :: DTDay
        
        integer :: NGI, NII
        
        real    :: DenitrificationRate, AnaNEf
        
        !------------------------------------------------------------------------
        
        DTDay       = Me%DTDay

        NGI         = Me%PropIndex%Ngas
        NII         = Me%PropIndex%Nitrate

        DenitrificationRate     = Me%SpecificRates%NitrateToNgas%Value
        AnaNEf                  = Me%Microorganisms%Anaerobic%EficiencyN
        
        Me%Matrix(NGI, NGI) =  1.

        !Sources: Denitrification
        Me%Matrix(NGI, NII) = Me%Matrix(NGI, NII) + DTDay * DenitrificationRate * (1-AnaNEf)
        
        !Independent term
        Me%IndTerm(NGI)     = Me%ExternalVar%Mass(NGI, index) 
    !----------------------------------------------------------------------------

    end subroutine Ngas 
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
