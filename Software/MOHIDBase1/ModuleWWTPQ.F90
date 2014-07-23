!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : WWTPQ
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : October 2013
! REVISION      : v1.0
! DESCRIPTION   : Zero-dimensional model for Water Quality in activated sludge
! process at WasteWater Treatment Plant.
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

!Module WWTPQ uses ModuleGlobalData, ModuleLUD, ModuleEnterData
!?uses use ModuleFunctions, only: OxygenSaturation, PhytoLightLimitationFactor, what does it mean? what for? 
!?uses all these modules in the end?

Module ModuleWWTPQ

    use ModuleGlobalData
    use ModuleLUD
    use ModuleEnterData

    implicit none

    private 
!?what does it mean?

    !Subroutines---------------------------------------------------------------
!? The commented subroutines are relly not needed? what do they do? If they are not needed, please remove them. 

    !Constructor
    public  :: StartWWTPQ
    private ::      AllocateInstance
    private ::      Nullify_all_Sub_Type_Pointers
    private ::      WWTPQReadData
!    private ::          WWTPQOptions
    private ::          WWTPQPropertyIndexNumber
    private ::          WWTPQReadCalcOptions
!    private ::          WWTPQOptionsConsistencyVerif
!    private ::          WWTPQConfiguration
    private ::          WWTPQReadFileConstants
    private ::      AllocateVariables
!    private ::          Add_PropRateFlux
!    private ::          Add_EquaRateFlux

!    public  ::          Construct_WWTPQRateFlux

     !Selector
    public  :: GetDTWWTPQM
!    public  :: GetWWTPQOptions
    public  :: GetWWTPQSize   
    public  :: GetWWTPQPropIndex
!    public  :: GetWWTPQPropRateFlux   
!    public  :: UnGetWWTPQPropRateFlux              

!    public  :: GetWWTPQNCRatio            

    !Modifier
    public  :: WWTPQ            
    private ::      StartWWTPQIteration
!!    private ::      WWTPQCoeficientsCalculation
!    private ::          WWTPQOxygen
!    private ::              WWTPQOxygenSaturation
!    private ::              WWTPQOxygenCalculation
!    private ::      WWTPQSystemResolution
!
!    private ::      WWTPQRatesCalculation
!
    !Destructor
    public  ::  KillWWTPQ
    private ::      DeallocateInstance
!
!
    !Management
    private ::      Ready
    private ::          LocateObjWWTPQ

    !Types---------------------------------------------------------------------
    
    !Stechiometric Parameters
    type T_StechParameters
        real :: HetBioYield                         = null_real !heterotrophic biomass yield, G COD G-1 COD
        real :: AutoBioYield                        = null_real !autotrophic biomass yield, G COD G-1 N
        real :: FracBioPartProd                     = null_real !fraction of biomass leading to particulate products, G COD G-1 COD
        real :: NCODBioMassRatio                    = null_real !mass N/mass COD in biomass, G N G-1 COD
        real :: NCODPartProdBioMassRatio            = null_real !mass N/mass COD in products from biomass, G N G-1 COD
    end type T_StechParameters

    !Kinetic Parameters
    type T_KineticParameters
    
    !heterotrophic organisms
        real :: MaxSpecGrowthRateHetBio             = null_real !maximum specific growth rate for heterotrophic biomass,
                                                                ! D-1, value at 20 ºC
        real :: HalfSatCoefHetBio                   = null_real !half-saturation coefficient for heterotrophic biomass, G COD M-3
        real :: OxygHalfSatCoefHetBio               = null_real !oxygen half-saturation coefficient for heterotrophic biomass, 
                                                                !G O2 M-3
        real :: NitHalfSatCoefDenHetBio             = null_real !nitrate half-saturation coefficient for denitrifying
                                                                ! heterotrophic biomass,G NO3-N M-3
        real :: DecayCoefHetBio                     = null_real !decay coefficient for heterotrophic biomass, D-1, value at 20 ºC
        real :: CorFacAnoxGrowthHetBio              = null_real !correction factor for anoxic growth of heterotrophs, dimensionless
        real :: AmmonifRate                         = null_real !ammonification rate, M3 COD (G.D)-1, value at 20 ºC
        real :: MaxSpecHydroRate                    = null_real !maximum specific hydrolysis rate, 
                                                                !G slowly biodegradable (COD G cell COD D)-1, value at 20 ºC
        real :: HalfSatCoefHydroSlowlyBioSub        = null_real !half-saturation coefficient for hydrolysis of slowly biodegradable 
                                                                !substrate, G slowly biodegradable COD G cell COD-1, value at 20 ºC
        real :: CorFacAnoxCond                      = null_real !correction factor under anoxic conditions, dimensionless
     
     !autotrophic organisms   
        real :: MaxSpecGrowthRateAutBio             = null_real !maximum specific growth rate for autotrophic biomass,
                                                                ! D-1, value at 20 ºC 
        real :: AmmoniaHalfSatCoefAutBio            = null_real !ammonia half-saturation coefficient for autotrophic biomass, 
                                                                !G NH3-N M-3
        real :: OxyHalfSatCoefAutBio                = null_real !oxygen half-saturation coefficient for autotrophic biomass,
                                                                ! G O2 M-3
        real :: DecayCoefAutBio                     = null_real !decay coefficient for autotrophic biomass, M3 COD (G D)-1
        
    end type T_KineticParameters
    
    !state variables
    type    T_PropIndex
        
        integer :: SolInertOrgMat                   = null_int !SI, soluble inert organic matter, G COD M-3
        integer :: ReadilyBioSub                    = null_int !SS, readily biodegradable substrate, G COD M-3
        integer :: PartInertOrgMar                  = null_int !XI, particualte inert organic matter, G COD M-3 
        integer :: SlowlyBioSub                     = null_int !XS, slowly biodegradable substrate, G COD M-3
        integer :: HetBio                           = null_int !XBH, activate heterotrophic biomass, G COD M-3
        integer :: AutBio                           = null_int !XBA, activate heterotrophic biomass, G COD M-3
        integer :: PartProd                         = null_int !XP, particulate products arising from biomass decay, G COD M-3 
        integer :: Oxygen                           = null_int !SO, oxygen, G (-COD) M-3 
        integer :: Nitrate                          = null_int !SNO, nitrate and nitrite, G NO3-N M-3  
        integer :: Ammonia                          = null_int !SNH, NH4+ + NH3 nitrogen, G NH3-N M-3
        integer :: SolBioOrgNitrogen                = null_int !SND, soluble biodegradable organic nitrogen, G N M-3 
        integer :: PartBioOrgNitrogen               = null_int !XND, particulate biodegradable organic matter, G N M-3 
        integer :: Alkalinity                       = null_int !SALK, alkalinity, mol M-3 
 
         end type    T_PropIndex
   
   !para ter como output as taxas ao longo do tempo
!    type    T_PropRateFlux
!        integer                              :: ID
!        real, pointer, dimension(:)          :: Field
!        type(T_PropRateFlux), pointer        :: Next,Prev !?o que +e isto?
!    end type T_PropRateFlux

   
!    type    T_EquaRateFlux
!        integer                              :: ID   !??? porquê duas vezes declardo o ID?
!        real                                 :: scalar = null_real
!        logical                              :: TimeSerie
!        type(T_EquaRateFlux), pointer        :: next,prev
!        type(T_PropRateFlux), pointer        :: FirstPropRateFlux
!        type(T_PropRateFlux), pointer        :: LastPropRateFlux
!
!    end type T_EquaRateFlux

   !o que é isto?
    type           T_ExtraRate
        integer                              :: ID
        real, pointer, dimension(:)          :: Field
    end type T_ExtraRate

   !o que é isto?
    type       T_PropCalc
        logical :: Nitrogen   = OFF
        logical :: Phosphorus = OFF
        logical :: Oxygen     = OFF
    end type T_PropCalc

   !o que é isto?
    type       T_CalcMethod
        logical :: ExplicitMethod = OFF
        logical :: ImplicitMethod = OFF
        logical :: SemiImpMethod  = OFF            
    end type T_CalcMethod

   !o que é isto?
    type       T_External
        real, pointer, dimension(:  ) :: Salinity
        real, pointer, dimension(:  ) :: Temperature
!        real, pointer, dimension(:  ) :: ShortWaveRadiation
!        real, pointer, dimension(:  ) :: LightExtCoefField
!        real, pointer, dimension(:  ) :: Thickness
        real, pointer, dimension(:,:) :: Mass      
    end type T_External

   !o que é isto?
    type      T_WWTPQ      
!        private
        integer                             :: InstanceID
        type(T_Size1D      )                :: Prop          
        type(T_PropIndex   )                :: PropIndex
        type(T_PropCalc    )                :: PropCalc
        type(T_CalcMethod  )                :: CalcMethod

        type(T_External    )                :: ExternalVar

!        type(T_EquaRateFlux), pointer       :: FirstEquaRateFlux
!        type(T_EquaRateFlux), pointer       :: LastEquaRateFlux

        type(T_ExtraRate   ), pointer       :: GrossProduction
        type(T_ExtraRate   ), pointer       :: TempLimitation
        type(T_ExtraRate   ), pointer       :: NutLimitation
        type(T_ExtraRate   ), pointer       :: NLimitation
        type(T_ExtraRate   ), pointer       :: PLimitation
        type(T_ExtraRate   ), pointer       :: LightLimitation

        real :: DTDay                                       = null_real
        real :: DTSecond                                    = null_real 
   
   !o que é isto?     
         type(T_KineticParameters)                   :: KineticParameters
         type(T_StechParameters)                     :: StechParameters

!
   !o que é isto?
        double precision, pointer, dimension(:,:)   :: Matrix
        real,             pointer, dimension(:  )   :: IndTerm
        real, pointer, dimension(:  )               :: NewMass   !Used with Explicit method
!
!        !Instance of Module_EnterData
        integer                                     :: ObjEnterData = 0
!
!        !Instance of ModuleLUD
        integer                                     :: ObjLUD = 0
!
        type(T_WWTPQ ), pointer              :: Next
!
    end type T_WWTPQ
!
!    !Global Module Variables
    type (T_WWTPQ), pointer                  :: FirstObjWWTPQ
    type (T_WWTPQ), pointer                  :: Me

    !--------------------------------------------------------------------------
    
    contains

!What is the meaning of contains? 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine StartWWTPQ(WWTPQID, FileName, STAT)
    
    !? What does it subroutine do?

        !Arguments-------------------------------------------------------------
        integer                         :: WWTPQID
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
        if (.not. ModuleIsRegistered(mWWTPQ_)) then
            nullify (FirstObjWWTPQ)
            call RegisterModule (mWWTPQ_) 
        endif
        
        call Ready(WWTPQID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then
            
            call AllocateInstance           

            call Nullify_all_Sub_Type_Pointers

            call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartWWTPQ - ModuleWWTPQ - ERR01'
              

            call WWTPQReadData
            call AllocateVariables

cd2:        if (.NOT. Me%CalcMethod%ExplicitMethod) then
                call StartLUD(Me%ObjLUD,                                                            &
                              Me%Prop%ILB,                                                          &
                              Me%Prop%IUB,                                                          &
                              Me%Prop%ILB,                                                          &
                              Me%Prop%IUB,                                                          &
                              STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'StartWWTPQ - ModuleWWTPQ - ERR02'
            end if cd2

                                                                                                        
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'StartWWTPQ - ModuleWWTPQ - ERR03'

            !Returns ID
            WWTPQID              = Me%InstanceID

            STAT_ = SUCCESS_
        else 
            
            stop 'ModuleWWTPQ - StartWWTPQ - ERR99' 

        end if cd0


        if (present(STAT))                                                                          &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartWWTPQ

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

 !? What does it subroutine do?

        !Local-----------------------------------------------------------------
        type (T_WWTPQ), pointer           :: NewObjWWTPQ
        type (T_WWTPQ), pointer           :: PreviousObjWWTPQ


        !Allocates new instance
        allocate (NewObjWWTPQ)
        nullify  (NewObjWWTPQ%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjWWTPQ)) then
            FirstObjWWTPQ        => NewObjWWTPQ
            Me                          => NewObjWWTPQ
        else
            PreviousObjWWTPQ     => FirstObjWWTPQ
            Me                          => FirstObjWWTPQ%Next
            do while (associated(Me))
                PreviousObjWWTPQ => Me
                Me                      => Me%Next
            enddo
            Me                          => NewObjWWTPQ
            PreviousObjWWTPQ%Next=> NewObjWWTPQ
        endif

        Me%InstanceID = RegisterNewInstance (mWWTPQ_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    Subroutine Nullify_all_Sub_Type_Pointers

 !? What does it subroutine do?

        nullify(Me%ExternalVar%Salinity    )
        nullify(Me%ExternalVar%Temperature )
        nullify(Me%ExternalVar%Mass        )
        nullify(Me%Matrix                  )
        nullify(Me%IndTerm                 )
        nullify(Me%NewMass                 )

    end Subroutine Nullify_all_Sub_Type_Pointers

    !----------------------------------------------------------------------------
!    !This subroutine adds a new property rateflux to the Rate Fluxes List
!    subroutine Add_EquaRateFlux(NewEquaRateFlux)
!
!!? What does it subroutine do?
!
!        !Arguments-------------------------------------------------------------
!        type(T_EquaRateFlux),pointer    :: NewEquaRateFlux
!
!        !----------------------------------------------------------------------
!
!        if (.not.associated(Me%FirstEquaRateFlux)) then
!
!            Me%FirstEquaRateFlux        => NewEquaRateFlux
!            Me%LastEquaRateFlux         => NewEquaRateFlux
!        else
!            NewEquaRateFlux%Prev        => Me%LastEquaRateFlux
!            Me%LastEquaRateFlux%Next    => NewEquaRateFlux
!            Me%LastEquaRateFlux         => NewEquaRateFlux
!        
!        end if 
!
!    end subroutine Add_EquaRateFlux 

    !--------------------------------------------------------------------------
  
!    subroutine Add_PropRateFlux(EquaRateFluxX, NewPropRateFlux)
!
!!? What does it subroutine do?
!
!        !Arguments-------------------------------------------------------------
!        type(T_EquaRateFlux),pointer              :: EquaRateFluxX
!        type(T_PropRateFlux),pointer              :: NewPropRateFlux
!
!        !----------------------------------------------------------------------
!
!        if (.not.associated(EquaRateFluxX%FirstPropRateFlux)) then
        
!            EquaRateFluxX%FirstPropRateFlux        => NewPropRateFlux
!            EquaRateFluxX%LastPropRateFlux         => NewPropRateFlux
!        
!        else
!            
!            NewPropRateFlux%Prev                   => EquaRateFluxX%LastPropRateFlux
!            EquaRateFluxX%LastPropRateFlux%Next    => NewPropRateFlux
!            EquaRateFluxX%LastPropRateFlux         => NewPropRateFlux
!        
!        end if 
!
!    end subroutine Add_PropRateFlux 
!
!    !--------------------------------------------------------------------------
    
!    subroutine WWTPQOptions  
!    
!    !? What does it subroutine do? Why nitrogen? phosphorous and oxygen? 
!
!        !External--------------------------------------------------------------
!        integer                         :: flag
!        integer                         :: FromFile
!        integer                         :: STAT_CALL
!
!        !----------------------------------------------------------------------
!
!        call GetExtractType(FromFile = FromFile)
!
!        call GetData(Me%PropCalc%Nitrogen,                                                          &
!                     Me%ObjEnterData, flag,                                                         &
!                     SearchType = FromFile,                                                         &
!                     keyword='NITROGEN',                                                            &
!                     ClientModule = 'ModuleWWTPQ',                                           &
!                     STAT       = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                                                &
!            stop 'Subroutine WWTPQOptions; Module ModuleWWTPQ. ERR01.' 
!
!
!        call GetData(Me%PropCalc%Phosphorus,                                                        &
!                     Me%ObjEnterData, flag,                                                         &
!                     SearchType = FromFile,                                                         &
!                     keyword='PHOSPHOR',                                                            &
!                     ClientModule = 'ModuleWWTPQ',                                           &
!                     STAT       = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                                                &
!            stop 'Subroutine WWTPQOptions; Module ModuleWWTPQ. ERR02.' 
!
!
!        call GetData(Me%PropCalc%Oxygen,                                                            &
!                     Me%ObjEnterData, flag,                                                         &
!                     SearchType = FromFile,                                                         &
!                     keyword='OXYGEN',                                                              &
!                     ClientModule = 'ModuleWWTPQ',                                           &
!                     STAT       = STAT_CALL)                        
!        if (STAT_CALL .NE. SUCCESS_)                                                                &
!            stop 'Subroutine WWTPQOptions; Module ModuleWWTPQ. ERR09.' 
!
!    end subroutine WWTPQOptions         
    
    !--------------------------------------------------------------------------
    
    !A subroutine WWTPQPropertyIndexNumber serve para atribuir indices as propriedades a ser 
    !calculadas no modulo de Qualidade da Agua
   
    subroutine WWTPQPropertyIndexNumber

!? What does it subroutine do? Give IndexNumber to each property

        Me%Prop%ILB = 1
        Me%Prop%IUB = 0
                    
        Me%Prop%IUB                                   = Me%Prop%IUB + 1
        Me%PropIndex%SolInertOrgMat                   = Me%Prop%IUB

        Me%Prop%IUB                                   = Me%Prop%IUB + 1
        Me%PropIndex%ReadilyBioSub                    = Me%Prop%IUB

        Me%Prop%IUB                                   = Me%Prop%IUB + 1
        Me%PropIndex%PartInertOrgMar                  = Me%Prop%IUB

        Me%Prop%IUB                                   = Me%Prop%IUB + 1
        Me%PropIndex%SlowlyBioSub                     = Me%Prop%IUB
                
        Me%Prop%IUB                                   = Me%Prop%IUB + 1
        Me%PropIndex%HetBio                           = Me%Prop%IUB

        Me%Prop%IUB                                   = Me%Prop%IUB + 1
        Me%PropIndex%AutBio                           = Me%Prop%IUB
            
        Me%Prop%IUB                                   = Me%Prop%IUB + 1
        Me%PropIndex%PartProd                         = Me%Prop%IUB

        Me%Prop%IUB                                   = Me%Prop%IUB + 1
        Me%PropIndex%Oxygen                           = Me%Prop%IUB

        Me%Prop%IUB                                   = Me%Prop%IUB + 1
        Me%PropIndex%Nitrate                          = Me%Prop%IUB

        Me%Prop%IUB                                   = Me%Prop%IUB + 1
        Me%PropIndex%Ammonia                          = Me%Prop%IUB
            
        Me%Prop%IUB                                   = Me%Prop%IUB + 1
        Me%PropIndex%SolBioOrgNitrogen                = Me%Prop%IUB

        Me%Prop%IUB                                   = Me%Prop%IUB + 1
        Me%PropIndex%PartBioOrgNitrogen               = Me%Prop%IUB

        Me%Prop%IUB                                   = Me%Prop%IUB + 1
        Me%PropIndex%Alkalinity                       = Me%Prop%IUB
                                      
    end subroutine WWTPQPropertyIndexNumber

    !--------------------------------------------------------------------------

    subroutine WWTPQReadData

!? What does it subroutine do? 

        !Local-----------------------------------------------------------------
        logical :: Consistent

        !----------------------------------------------------------------------
        

 !       call WWTPQOptions
        call WWTPQPropertyIndexNumber
        call WWTPQReadCalcOptions

        !Consistent = WWTPQOptionsConsistencyVerif ()
    
        !if (Consistent) then
            !call WWTPQConfiguration
            call WWTPQReadFileConstants
!        else
!            write(*,*) 
!            write(*,*) 'The Water Quality Options were not consistent, verify file data.'
!            stop       'SUBROUTINE WWTPQReadData; Module ModuleWWTPQ. ERR01'
!        endif   !Consistent

        !----------------------------------------------------------------------

    end subroutine WWTPQReadData

    !--------------------------------------------------------------------------

    subroutine WWTPQReadCalcOptions
    
    !? What does it subroutine do? Verify which method is used to calculate (explicit, implicit, semi-explicit?)

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
                                 ClientModule = 'ModuleWWTPQ',                                      &
                                 STAT       = STAT_CALL)                                            
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WWTPQReadCalcOptions; Module ModuleWWTPQ. ERR01.' 


        !Verifica se se pretende calcular usando um metodo IMPLICITO
        call GetData(Me%CalcMethod%ImplicitMethod,                                                  &
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='IMPLICIT',                                                            &    
                     ClientModule = 'ModuleWWTPQ',                                                  &
                     STAT       = STAT_CALL)                        
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WWTPQReadCalcOptions; Module ModuleWWTPQ. ERR02.' 

            !Verifica se se pretende calcular usando um metodo IMPLICITO/EXPLICITO        
        call GetData(Me%CalcMethod%SemiimpMethod,                                                   & 
                     Me%ObjEnterData, flag,                                                         &
                     SearchType = FromFile,                                                         &
                     keyword='SEMIIMP',                                                             &    
                     ClientModule = 'ModuleWWTPQ',                                                  &
                     STAT       = STAT_CALL)                         
        if (STAT_CALL .NE. SUCCESS_)                                                                &
            stop 'Subroutine WWTPQReadCalcOptions; Module ModuleWWTPQ. ERR03.'         

    end subroutine WWTPQReadCalcOptions

    !--------------------------------------------------------------------------
    
    subroutine WWTPQReadFileConstants

!? What does it subroutine do? Read constants

        !External--------------------------------------------------------------
        integer                         :: FromFile
        integer                         :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                         :: flag

        !----------------------------------------------------------------------

        call GetExtractType(FromFile = FromFile)


!?Why call for dtsecond is different than for other constants? 

cd1 :   if (Me%dtsecond .le. 0.0) then
            !dtsecond, time step, in seconds, between 2 wwtpq calls 
            call GetData(           Me%dtsecond,                                    &
                                    Me%objenterdata, flag,                          &
                                    searchtype = fromfile, keyword='DTSECONDS',     & 
                                    default    = 60.0 * 60.0,                       & !1 hour
                                    clientmodule = 'modulewwtpq',                   &
                                    stat       = stat_call)
            if (stat_call .ne. success_)                                            &
                stop 'subroutine wwtpqreadfileconstants; module modulewwtpq. err00.' 

cd22 :      if (flag .eq. 0) then
                write(*,*) 
                write(*,*) 'keyword dtseconds not found in water quality data file.'
                write(*,*) 'subroutine wwtpqreadfileconstants; module modulewwtpq. wrn01.'
                write(*,*) 'assumed ', me%dtsecond, &
                            'seconds (',  me%dtsecond / 3600.0, 'hour).'
                write(*,*) 
            end if cd22
        end if cd1

        
        !for compatibility with the rest of the program, !dtseconds converted to day!
        me%dtday = me%dtsecond / 24.0 / 60.0 / 60.0

        !Reads non specific rates & constants--------------------------------------

        !HetBioYield, heterotrophic biomass yield, G COD G-1 COD
        call GetData(           Me%StechParameters%HetBioYield,                     &
                                Me%ObjEnterData, flag,                              &
                                SearchType = FromFile,                              &
                                keyword ='HET_BIO_YIELD',                           & 
                                default    = 0.670,                                 &
                                ClientModule = 'ModuleWWTPQ',                       &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.' 

        !AutoBioYield, autotrophic biomass yield, G COD G-1 N
        call GetData(           Me%StechParameters%AutoBioYield,                    &
                                Me%ObjEnterData, flag,                              &
                                SearchType = FromFile,  keyword ='AUT_BIO_YIELD',   & 
                                default    = 0.240,                                 &
                                ClientModule = 'ModuleWWTPQ',                       &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.' 

        !FracBioPartProd, fraction of biomass leading to particulate products, G COD G-1 COD
        call GetData(           Me%StechParameters%FracBioPartProd,                            &
                                Me%ObjEnterData, flag,                                         &
                                SearchType = FromFile,  keyword ='FRAC_BIO_PART_PROD',         & 
                                default    = 0.080,                                            &
                                ClientModule = 'ModuleWWTPQ',                                  &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                           &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.' 

        !NCODBioMassRatio, mass N/mass COD in biomass, G N G-1 COD
        call GetData(           Me%StechParameters%NCODBioMassRatio,                          &
                                Me%ObjEnterData, flag,                                        &
                                SearchType = FromFile,  keyword ='MASSN_MASSCOD_BIO_RATIO',   & 
                                default    = 0.086,                                           &
                                ClientModule = 'ModuleWWTPQ',                                 &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                          &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.' 

        !NCODPartProdBioMassRatio, mass N/mass COD in products from biomass, G N G-1 COD
        call GetData(           Me%StechParameters%NCODPartProdBioMassRatio,                            &
                                Me%ObjEnterData, flag,                                                  &
                                SearchType = FromFile,  keyword ='MASSN_MASSCOD_PART_PROD_BIO_RATIO',   & 
                                default    = 0.060,                                                     &
                                ClientModule = 'ModuleWWTPQ',                                           &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                    &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.' 

        !MaxSpecGrowthRateHetBio, maximum specific growth rate for heterotrophic biomass, D-1, value at 20 ºC
        call GetData(           Me%KineticParameters%MaxSpecGrowthRateHetBio,                                       &
                                Me%ObjEnterData, flag,                                                              &
                                SearchType = FromFile,  keyword ='MAX_SPEC_GROWTH_RATE_HET_BIO',                    & 
                                default    = 6.0,                                                                   &
                                ClientModule = 'ModuleWWTPQ',                                                       &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                                &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.' 
            
        !HalfSatCoefHetBio, half-saturation coefficient for heterotrophic biomass, G COD M-3
        call GetData(           Me%KineticParameters%HalfSatCoefHetBio,                         &
                                Me%ObjEnterData, flag,                                          &
                                SearchType = FromFile,  keyword ='HALF_SAT_COEF_HET_BIO',       & 
                                default    = 20.0,                                              &
                                ClientModule = 'ModuleWWTPQ',                                   &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                            &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.' 
  
        !OxygHalfSatCoefHetBio, oxygen half-saturation coefficient for heterotrophic biomass, G O2 M-3
        call GetData(           Me%KineticParameters%OxygHalfSatCoefHetBio,                               &
                                Me%ObjEnterData, flag,                                                    &
                                SearchType = FromFile,  keyword ='OXYG_HALF_SAT_COEF_HET_BIO',            & 
                                default    = 0.20,                                                        &
                                ClientModule = 'ModuleWWTPQ',                                             &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                      &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.' 
            
        !NitHalfSatCoefDenHetBio, nitrate half-saturation coefficient for denitrifying heterotrophic biomass, G NO3-N M-3
        call GetData(           Me%KineticParameters%NitHalfSatCoefDenHetBio,                                                 &
                                Me%ObjEnterData, flag,                                                                        &
                                SearchType = FromFile,  keyword ='NIT_HALF_SAT_COEF_DEN_HET_BIO',                             & 
                                default    = 0.50,                                                                            &
                                ClientModule = 'ModuleWWTPQ',                                                                 &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                                          &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.' 
            
        !DecayCoefHetBio, decay coefficient for heterotrophic biomass, D-1, value at 20 ºC
        call GetData(           Me%KineticParameters%DecayCoefHetBio,                         &
                                Me%ObjEnterData, flag,                                        &
                                SearchType = FromFile,  keyword ='DECAY_COEF_HET_BIO',        & 
                                default    = 0.62,                                            &
                                ClientModule = 'ModuleWWTPQ',                                 &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                          &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.' 
          
        !CorFacAnoxGrowthHetBio, correction factor for anoxic growth of heterotrophs, dimensionless
        call GetData(           Me%KineticParameters%CorFacAnoxGrowthHetBio,                               &
                                Me%ObjEnterData, flag,                                                     &
                                SearchType = FromFile,  keyword ='COR_FAC_ANO_GROWTH_HET_BIO',             & 
                                default    = 0.8,                                                          &
                                ClientModule = 'ModuleWWTPQ',                                              &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                       &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.' 
 
    
        !AmmonifRate, ammonification rate, M3 COD (G.D)-1, value at 20 ºC
        call GetData(           Me%KineticParameters%AmmonifRate,                   &
                                Me%ObjEnterData, flag,                              &
                                SearchType = FromFile,  keyword ='AMMONIF_RATE',    & 
                                default    = 0.08,                                  &
                                ClientModule = 'ModuleWWTPQ',                       &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.'          
      
       !MaxSpecHydroRate, maximum specific hydrolysis rate, G slowly biodegradable (COD G cell COD D)-1, value at 20 ºC
        call GetData(           Me%KineticParameters%MaxSpecHydroRate,                                                     &
                                Me%ObjEnterData, flag,                                                                     &
                                SearchType = FromFile,  keyword ='MAX_SPEC_HYDRO_RATE',                                    & 
                                default    = 3.0,                                                                          &
                                ClientModule = 'ModuleWWTPQ',                                                              &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                                       &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.' 
       
  
       !HalfSatCoefHydroSlowlyBioSub, half-saturation coefficient for hydrolysis of slowly biodegradable substrate, 
       !G slowly biodegradable COD G cell COD-1, value at 20 ºC
        call GetData(           Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub,                                       &
                                Me%ObjEnterData, flag,                                                                   &
                                SearchType = FromFile,  keyword ='HALF_SAT_COEF_HYDRO_SLOWLY_BIO_SUB',                   & 
                                default    = 0.03,                                                                       &
                                ClientModule = 'ModuleWWTPQ',                                                            &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                                     &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.' 
            
        !CorFacAnoxCond, correction factor under anoxic conditions, dimensionless
        call GetData(           Me%KineticParameters%CorFacAnoxCond,                       &
                                Me%ObjEnterData, flag,                                     &
                                SearchType = FromFile,  keyword ='COR_FAC_ANOX_COND',      & 
                                default    = 0.4,                                          &
                                ClientModule = 'ModuleWWTPQ',                              &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                       &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.'    
            
        !MaxSpecGrowthRateAutBio, maximum specific growth rate for autotrophic biomass, D-1, value at 20 ºC 
        call GetData(           Me%KineticParameters%MaxSpecGrowthRateAutBio,                                       &
                                Me%ObjEnterData, flag,                                                              &
                                SearchType = FromFile,  keyword ='MAX_SPEC_GROWTH_RATE_AUT_BIO',                    & 
                                default    = 0.80,                                                                  &
                                ClientModule = 'ModuleWWTPQ',                                                       &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                                &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.'    
            
        !AmmoniaHalfSatCoefAutBio, ammonia half-saturation coefficient for autotrophic biomass, G NH3-N M-3 
        call GetData(           Me%KineticParameters%AmmoniaHalfSatCoefAutBio,                                     &
                                Me%ObjEnterData, flag,                                                             &
                                SearchType = FromFile,  keyword ='AMMONIA_HALF_SAT_COEF_AUT_BIO',                  & 
                                default    = 1.0,                                                                  &
                                ClientModule = 'ModuleWWTPQ',                                                      &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                               &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.'  
            
        !OxyHalfSatCoefAutBio, oxygen half-saturation coefficient for autotrophic biomass, G O2 M-3 
        call GetData(           Me%KineticParameters%OxyHalfSatCoefAutBio,                             &
                                Me%ObjEnterData, flag,                                                 &
                                SearchType = FromFile,  keyword ='OXY_HALF_SAT_COEF_AUT_BIO',          & 
                                default    = 0.4,                                                      &
                                ClientModule = 'ModuleWWTPQ',                                          &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                                   &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.'  
            
   
        !DecayCoefAutBio, decay coefficient for autotrophic biomass, M3 COD (G D)-1
        call GetData(           Me%KineticParameters%DecayCoefAutBio,                      &
                                Me%ObjEnterData, flag,                                     &
                                SearchType = FromFile,  keyword ='DECAY_COEF_AUT_BIO',     & 
                                default    = 0.05,                                         &
                                ClientModule = 'ModuleWWTPQ',                &
                                STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine WWTPQReadFileConstants; Module ModuleWWTPQ. ERR01.'  
        
    
     end subroutine WWTPQReadFileConstants
   
    !----------------------------------------------------------------------------
    
    subroutine AllocateVariables
    
   !? What does it subroutine do? Create a matrix? 

        !External----------------------------------------------------------------
        integer :: STAT_CALL

        !Local-------------------------------------------------------------------
        integer :: PropLB, PropUB

        !------------------------------------------------------------------------

        PropLB    = Me%Prop%ILB
        PropUB    = Me%Prop%IUB

        allocate(Me%Matrix (PropLB:PropUB, PropLB:PropUB)) !rates matrix
        allocate(Me%IndTerm(PropLB:PropUB               )) !old concentration proprieties colunm


cd1 :   if (Me%CalcMethod%ExplicitMethod) then
            allocate(Me%NewMass(PropLB:PropUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'Subroutine AllocateVariables; module ModuleWWTPQ. ERR01.'
        end if cd1

        !------------------------------------------------------------------------

    end subroutine AllocateVariables   

    !----------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetWWTPQSize(WWTPQID, PropLB, PropUB, STAT)
    
    !? What does it subroutine do? 

        !Arguments-------------------------------------------------------------
        integer                         :: WWTPQID
        integer, optional, intent(OUT)  :: PropLB,    PropUB
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------

        integer                         :: ready_              

        !Local-----------------------------------------------------------------
        integer                         :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WWTPQID, ready_)    
        
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

    end subroutine GetWWTPQSize
!
!    !--------------------------------------------------------------------------

    subroutine GetWWTPQPropIndex(WWTPQID,                                   &
                                      SolInertOrgMat,                       &
                                      ReadilyBioSub,                        &
                                      PartInertOrgMar,                      & 
                                      SlowlyBioSub,                         &  
                                      HetBio,                               &
                                      AutBio,                               &
                                      PartProd,                             &              
                                      Oxygen,                               &
                                      Ammonia,                              &
                                      Nitrate,                              & 
                                      SolBioOrgNitrogen,                    &
                                      PartBioOrgNitrogen,                   &
                                      Alkalinity,                           &
                                      STAT)    

!? What does it subroutine do? Get the PropIndex of each propriety. 
!? Why some subroutine have arguments and others do not? 


        !Arguments-------------------------------------------------------------
        integer                         :: WWTPQID
                                        
        integer, optional, intent(OUT)  :: STAT
!       
        integer, optional, intent(OUT)  :: SolInertOrgMat   
        integer, optional, intent(OUT)  :: ReadilyBioSub 
        integer, optional, intent(OUT)  :: PartInertOrgMar
        integer, optional, intent(OUT)  :: SlowlyBioSub
        integer, optional, intent(OUT)  :: HetBio 
        integer, optional, intent(OUT)  :: AutBio
        integer, optional, intent(OUT)  :: PartProd
        integer, optional, intent(OUT)  :: Oxygen
        integer, optional, intent(OUT)  :: Ammonia
        integer, optional, intent(OUT)  :: Nitrate
        integer, optional, intent(OUT)  :: SolBioOrgNitrogen
        integer, optional, intent(OUT)  :: PartBioOrgNitrogen
        integer, optional, intent(OUT)  :: Alkalinity

        !External--------------------------------------------------------------
        integer :: ready_              

        !Local-----------------------------------------------------------------
        integer                         :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WWTPQID, ready_)    
!        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
           
                   
             if (present(SolInertOrgMat                 )) SolInertOrgMat          = Me%PropIndex%SolInertOrgMat
             if (present(ReadilyBioSub                  )) ReadilyBioSub           = Me%PropIndex%ReadilyBioSub
             if (present(PartInertOrgMar                )) PartInertOrgMar         = Me%PropIndex%PartInertOrgMar
             if (present(SlowlyBioSub                   )) SlowlyBioSub            = Me%PropIndex%SlowlyBioSub
             if (present(HetBio                         )) HetBio                  = Me%PropIndex%HetBio
             if (present(AutBio                         )) AutBio                  = Me%PropIndex%AutBio
             if (present(PartProd                       )) PartProd                = Me%PropIndex%PartProd
             if (present(Oxygen                         )) Oxygen                  = Me%PropIndex%Oxygen
             if (present(Ammonia                        )) Ammonia                 = Me%PropIndex%Ammonia
             if (present(Nitrate                        )) Nitrate                 = Me%PropIndex%Nitrate
             if (present(SolBioOrgNitrogen              )) SolBioOrgNitrogen       = Me%PropIndex%SolBioOrgNitrogen
             if (present(PartBioOrgNitrogen             )) PartBioOrgNitrogen      = Me%PropIndex%PartBioOrgNitrogen
             if (present(Alkalinity                     )) Alkalinity              = Me%PropIndex%Alkalinity
             
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetWWTPQPropIndex
!
!    !--------------------------------------------------------------------------
!
    subroutine GetDTWWTPQM(WWTPQID, DTDay, DTSecond, STAT)
    
  !? What does it subroutine do?  

        !Arguments-------------------------------------------------------------
        integer                         :: WWTPQID
        real,    optional, intent(OUT)  :: DTDay
        real,    optional, intent(OUT)  :: DTSecond
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                         :: ready_              

        !Local-----------------------------------------------------------------
        integer                         :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WWTPQID, ready_)    
        
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

    end subroutine GetDTWWTPQM

!--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine WWTPQ(WWTPQID,                                                   &
                            Salinity,                                           &
                            Temperature,                                        &
                            ShortWaveRadiation,                                 &
                            LightExtCoefField,                                  &
                            Thickness,                                          &
                            Mass,                                               &
                            WWTPQArrayLB, WWTPQArrayUB,                         &
                            OpenPoints,                                         &               
                            FishFood,                                           &
                            STAT)  

!? What does it subroutine do?

        !Arguments---------------------------------------------------------------
        integer                                       :: WWTPQID
        real,                 pointer, dimension(:  ) :: Salinity
        real,                 pointer, dimension(:  ) :: Temperature
        real,                 pointer, dimension(:  ) :: ShortWaveRadiation
        real,                 pointer, dimension(:  ) :: LightExtCoefField
        real,                 pointer, dimension(:  ) :: Thickness
        real,                 pointer, dimension(:,:) :: Mass
        integer, optional,    pointer, dimension(:  ) :: OpenPoints
        real,    optional,    pointer, dimension(:  ) :: FishFood
        integer,              intent(IN )             :: WWTPQArrayLB, WWTPQArrayUB  
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

        call Ready(WWTPQID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then
!
!            Me%ExternalVar%Salinity                   => Salinity
!            if (.NOT. associated(Me%ExternalVar%Salinity))         &
!                stop 'Subroutine WWTPQ; Module ModuleWWTPQ. ERR01' 
!
!
            Me%ExternalVar%Temperature                => temperature
            if (.not. associated(me%externalvar%temperature))        &
                stop 'subroutine wwtpq; module modulewwtpq. err02'

            Me%ExternalVar%Mass                       => Mass
            if (.NOT. associated(Me%ExternalVar%Mass))               &
                stop 'Subroutine WWTPQ; Module ModuleWWTPQ. ERR03.'

            call StartWWTPQIteration

do1 :       do index = WWTPQArrayLB, WWTPQArrayUB
            
            !If this module is called from the WWTPQ3D module, OpenPoint is present
            !and the WWTPQ module runs for all Openpoints
            !If this module is called from the Lagrangian module, OpenPoint is not present
            !and the WWTPQ module runs for all volumes
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
                call WWTPQCoeficientsCalculation   (index)

                !The rates can just be calculated if the rate flux is associated
                !In the case that this module is used by the lagrangian module
                !the rate fluxes are not calculated
                !Rates must be computed before call to WWTPQSystemResolution to use 
                !old concentrations 
!                if (associated(Me%FirstEquaRateFlux)) then
!!                    call WWTPQRatesCalculation      (index)
!               end if

                call WWTPQSystemResolution         (index)
                
            end if
            end do do1
          
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                      &
            STAT = STAT_

!        !------------------------------------------------------------------------
!
    end subroutine WWTPQ
!
    !----------------------------------------------------------------------------

    subroutine StartWWTPQIteration

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

    end subroutine StartWWTPQIteration

    !----------------------------------------------------------------------------

    subroutine WWTPQCoeficientsCalculation(index)
    
    !? What does it subroutine do?

        !Arguments-------------------------------------------------------------

        integer, intent(IN) :: index

        call WWTPQSolInertOrgMat        (index)
        call WWTPQReadilyBioSub         (index)
        call WWTPQPartInertOrgMar       (index)
        call WWTPQSlowlyBioSub          (index)
        call WWTPQHetBio                (index)
        call WWTPQAutBio                (index)
        call WWTPQPartProd              (index)
        call WWTPQOxygen                (index)
        call WWTPQNitrate               (index)  
        call WWTPQAmmonia               (index)  
        call WWTPQSolBioOrgNitrogen     (index)
        call WWTPQPartBioOrgNitrogen    (index)  
        call WWTPQAlkalinity            (index)

      !----------------------------------------------------------------------

    end subroutine WWTPQCoeficientsCalculation

!    !--------------------------------------------------------------------------

    subroutine WWTPQSystemResolution(index)
    
    !? What does it subroutine do? Calculate each new propriety using the method choosen

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
                stop 'Subroutine WWTPQSystemResolution; module ModuleWWTPQ. ERR03.'




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
                stop 'Subroutine WWTPQSystemResolution; module ModuleWWTPQ. ERR04.'


do33 :      do prop = PropLB, PropUB
                Me%ExternalVar%Mass(prop, index) = x(prop)
            end do do33
            
       
            nullify   (x)
        end if cd1

        !----------------------------------------------------------------------

    end subroutine WWTPQSystemResolution


 subroutine WWTPQSolInertOrgMat  (index)
    
    !Subroutine to calculate the propriety SolInertOrgMat
    

    !Arguments---------------------------------------------------------------

!? what is the meaning of arguments?

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
    
    !? what is the meaning of local?
    
    !State variables
        integer :: SolInertOrgMat
        integer :: ReadilyBioSub
        integer :: PartInertOrgMar
        integer :: SlowlyBioSub
        integer :: HetBio
        integer :: AutBio
        integer :: PartProd
        integer :: Oxygen
        integer :: Nitrate
        integer :: Ammonia
        integer :: SolBioOrgNitrogen
        integer :: PartBioOrgNitrogen
        integer :: Alkalinity

        real :: HetBioYield
        real :: AutoBioYield
        real :: FracBioPartProd
        real :: NCODBioMassRatio
        real :: NCODPartProdBioMassRatio

        real :: MaxSpecGrowthRateHetBio
        real :: HalfSatCoefHetBio
        real :: OxygHalfSatCoefHetBio
        real :: NitHalfSatCoefDenHetBio
        real :: DecayCoefHetBio
        real :: CorFacAnoxGrowthHetBio
        real :: AmmonifRate
        real :: MaxSpecHydroRate
        real :: HalfSatCoefHydroSlowlyBioSub
        real :: CorFacAnoxCond
        real :: MaxSpecGrowthRateAutBio
        real :: AmmoniaHalfSatCoefAutBio
        real :: OxyHalfSatCoefAutBio
        real :: DecayCoefAutBio
 
         real :: DTDay
    
!------------------------------------------------------------------------
!State variables
        SolInertOrgMat                 = Me%PropIndex%SolInertOrgMat
        ReadilyBioSub                  = Me%PropIndex%ReadilyBioSub
        PartInertOrgMar                = Me%PropIndex%PartInertOrgMar 
        SlowlyBioSub                   = Me%PropIndex%SlowlyBioSub
        HetBio                         = Me%PropIndex%HetBio
        AutBio                         = Me%PropIndex%AutBio
        PartProd                       = Me%PropIndex%PartProd
        Oxygen                         = Me%PropIndex%Oxygen
        Nitrate                        = Me%PropIndex%Nitrate
        Ammonia                        = Me%PropIndex%Ammonia
        SolBioOrgNitrogen              = Me%PropIndex%SolBioOrgNitrogen
        PartBioOrgNitrogen             = Me%PropIndex%PartBioOrgNitrogen 
        Alkalinity                     = Me%PropIndex%Alkalinity 

!Stecheo Parameters
        HetBioYield                    = Me%StechParameters%HetBioYield
        AutoBioYield                   = Me%StechParameters%AutoBioYield
        FracBioPartProd                = Me%StechParameters%FracBioPartProd
        NCODBioMassRatio               = Me%StechParameters%NCODBioMassRatio
        NCODPartProdBioMassRatio       = Me%StechParameters%NCODPartProdBioMassRatio

!Kinetic Parameters
        MaxSpecGrowthRateHetBio        = Me%KineticParameters%MaxSpecGrowthRateHetBio
        HalfSatCoefHetBio              = Me%KineticParameters%HalfSatCoefHetBio
        OxygHalfSatCoefHetBio          = Me%KineticParameters%OxygHalfSatCoefHetBio
        NitHalfSatCoefDenHetBio        = Me%KineticParameters%NitHalfSatCoefDenHetBio
        DecayCoefHetBio                = Me%KineticParameters%DecayCoefHetBio
        CorFacAnoxGrowthHetBio         = Me%KineticParameters%CorFacAnoxGrowthHetBio
        AmmonifRate                    = Me%KineticParameters%AmmonifRate
        MaxSpecHydroRate               = Me%KineticParameters%MaxSpecHydroRate
        HalfSatCoefHydroSlowlyBioSub   = Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub
        CorFacAnoxCond                 = Me%KineticParameters%CorFacAnoxCond
        MaxSpecGrowthRateAutBio        = Me%KineticParameters%MaxSpecGrowthRateAutBio
        AmmoniaHalfSatCoefAutBio       = Me%KineticParameters%AmmoniaHalfSatCoefAutBio
        OxyHalfSatCoefAutBio           = Me%KineticParameters%OxyHalfSatCoefAutBio
        DecayCoefAutBio                = Me%KineticParameters%DecayCoefAutBio

        DTDay                          = Me%DTDay

                 
 !Calculation of system coeficients---------------------------------------
         Me%Matrix(SolInertOrgMat, SolInertOrgMat) = 1.0
         
 !Independent term
        Me%IndTerm(SolInertOrgMat) = Me%ExternalVar%Mass(SolInertOrgMat, index)
        
 !------------------------------------------------------------------------

    end subroutine WWTPQSolInertOrgMat 

    subroutine WWTPQReadilyBioSub (index)
    
    !Subroutine to calculate the propriety ReadilyBioSub

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        integer :: SolInertOrgMat
        integer :: ReadilyBioSub
        integer :: PartInertOrgMar
        integer :: SlowlyBioSub
        integer :: HetBio
        integer :: AutBio
        integer :: PartProd
        integer :: Oxygen
        integer :: Nitrate
        integer :: Ammonia
        integer :: SolBioOrgNitrogen
        integer :: PartBioOrgNitrogen
        integer :: Alkalinity

        real :: HetBioYield
        real :: AutoBioYield
        real :: FracBioPartProd
        real :: NCODBioMassRatio
        real :: NCODPartProdBioMassRatio

        real :: MaxSpecGrowthRateHetBio
        real :: HalfSatCoefHetBio
        real :: OxygHalfSatCoefHetBio
        real :: NitHalfSatCoefDenHetBio
        real :: DecayCoefHetBio
        real :: CorFacAnoxGrowthHetBio
        real :: AmmonifRate
        real :: MaxSpecHydroRate
        real :: HalfSatCoefHydroSlowlyBioSub
        real :: CorFacAnoxCond
        real :: MaxSpecGrowthRateAutBio
        real :: AmmoniaHalfSatCoefAutBio
        real :: OxyHalfSatCoefAutBio
        real :: DecayCoefAutBio
 
        real :: DTDay
        real :: ReadilyBioSubRateHetBio             = null_real
        
        real :: ReadilyBioSubRateHetBioAerGrowth    = null_real
        real :: ReadilyBioSubRateHetBioAnoGrowth    = null_real
        real :: ReadilyBioSubRateHetBioHydrolysis   = null_real
        
        !?why some = null_real and other not? because after you say that other are equal to something? 

    !------------------------------------------------------------------------
!State variables
        SolInertOrgMat                 = Me%PropIndex%SolInertOrgMat
        ReadilyBioSub                  = Me%PropIndex%ReadilyBioSub
        PartInertOrgMar                = Me%PropIndex%PartInertOrgMar 
        SlowlyBioSub                   = Me%PropIndex%SlowlyBioSub
        HetBio                         = Me%PropIndex%HetBio
        AutBio                         = Me%PropIndex%AutBio
        PartProd                       = Me%PropIndex%PartProd
        Oxygen                         = Me%PropIndex%Oxygen
        Nitrate                        = Me%PropIndex%Nitrate
        Ammonia                        = Me%PropIndex%Ammonia
        SolBioOrgNitrogen              = Me%PropIndex%SolBioOrgNitrogen
        PartBioOrgNitrogen             = Me%PropIndex%PartBioOrgNitrogen 
        Alkalinity                     = Me%PropIndex%Alkalinity 

!Stecheo Parameters
        HetBioYield                    = Me%StechParameters%HetBioYield
        AutoBioYield                   = Me%StechParameters%AutoBioYield
        FracBioPartProd                = Me%StechParameters%FracBioPartProd
        NCODBioMassRatio               = Me%StechParameters%NCODBioMassRatio
        NCODPartProdBioMassRatio       = Me%StechParameters%NCODPartProdBioMassRatio

!Kinetic Parameters
        MaxSpecGrowthRateHetBio        = Me%KineticParameters%MaxSpecGrowthRateHetBio
        HalfSatCoefHetBio              = Me%KineticParameters%HalfSatCoefHetBio
        OxygHalfSatCoefHetBio          = Me%KineticParameters%OxygHalfSatCoefHetBio
        NitHalfSatCoefDenHetBio        = Me%KineticParameters%NitHalfSatCoefDenHetBio
        DecayCoefHetBio                = Me%KineticParameters%DecayCoefHetBio
        CorFacAnoxGrowthHetBio         = Me%KineticParameters%CorFacAnoxGrowthHetBio
        AmmonifRate                    = Me%KineticParameters%AmmonifRate
        MaxSpecHydroRate               = Me%KineticParameters%MaxSpecHydroRate
        HalfSatCoefHydroSlowlyBioSub   = Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub
        CorFacAnoxCond                 = Me%KineticParameters%CorFacAnoxCond
        MaxSpecGrowthRateAutBio        = Me%KineticParameters%MaxSpecGrowthRateAutBio
        AmmoniaHalfSatCoefAutBio       = Me%KineticParameters%AmmoniaHalfSatCoefAutBio
        OxyHalfSatCoefAutBio           = Me%KineticParameters%OxyHalfSatCoefAutBio
        DecayCoefAutBio                = Me%KineticParameters%DecayCoefAutBio

        DTDay                          = Me%DTDay

!ReadilyBioSubRateHetBio, ReadilyBioSub Rate in function of HetBio, D-1
              
          ReadilyBioSubRateHetBioAerGrowth = (-1./HetBioYield)*MaxSpecGrowthRateHetBio*                                   &
           (Me%ExternalVar%Mass(ReadilyBioSub, index)/(HalfSatCoefHetBio+Me%ExternalVar%Mass(ReadilyBioSub, index)))      &
           *(Me%ExternalVar%Mass(Oxygen, index)/(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))      
                                       
          ReadilyBioSubRateHetBioAnoGrowth = (-1./HetBioYield)*MaxSpecGrowthRateHetBio*                                   &
          (Me%ExternalVar%Mass(ReadilyBioSub, index)                                                                      &
           /(HalfSatCoefHetBio+Me%ExternalVar%Mass(ReadilyBioSub, index)))*(OxygHalfSatCoefHetBio                         &
           /(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))                                                   &
           *(Me%ExternalVar%Mass(Nitrate, index)/(NitHalfSatCoefDenHetBio+Me%ExternalVar%Mass(Nitrate, index)))           &
           *CorFacAnoxGrowthHetBio  
                       
           if (Me%ExternalVar%Mass(HetBio, index) == 0) then
           
           ReadilyBioSubRateHetBioHydrolysis = 0
           
           else
           
            ReadilyBioSubRateHetBioHydrolysis = MaxSpecHydroRate*((Me%ExternalVar%Mass(SlowlyBioSub, index)               &
          /Me%ExternalVar%Mass(HetBio, index))                                                                            &
           /(HalfSatCoefHydroSlowlyBioSub+Me%ExternalVar%Mass(SlowlyBioSub, index)/Me%ExternalVar%Mass(HetBio, index)))   &
           *((Me%ExternalVar%Mass(Oxygen, index)                                                                          &
           /(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))+CorFacAnoxCond*(OxygHalfSatCoefHetBio             &
           /(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))*(Me%ExternalVar%Mass(Nitrate, index)              &
           /(NitHalfSatCoefDenHetBio+Me%ExternalVar%Mass(Nitrate, index))))    
           
           endif
            
             
           ReadilyBioSubRateHetBio = ReadilyBioSubRateHetBioAerGrowth + ReadilyBioSubRateHetBioAnoGrowth                  &
           + ReadilyBioSubRateHetBioHydrolysis
                                                  
 !Calculation of system coeficients---------------------------------------
         Me%Matrix(ReadilyBioSub,HetBio ) = DTDay * (-ReadilyBioSubRateHetBio) 
         Me%Matrix(ReadilyBioSub, ReadilyBioSub) = 1.0

 !Independent term
        Me%IndTerm(ReadilyBioSub) = Me%ExternalVar%Mass(ReadilyBioSub, index)

 !------------------------------------------------------------------------

    end subroutine WWTPQReadilyBioSub

    subroutine WWTPQPartInertOrgMar  (index)
    
    !Subroutine to calculate the propriety PartInertOrgMar
    
    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        !State variables
        integer :: SolInertOrgMat
        integer :: ReadilyBioSub
        integer :: PartInertOrgMar
        integer :: SlowlyBioSub
        integer :: HetBio
        integer :: AutBio
        integer :: PartProd
        integer :: Oxygen
        integer :: Nitrate
        integer :: Ammonia
        integer :: SolBioOrgNitrogen
        integer :: PartBioOrgNitrogen
        integer :: Alkalinity

        real :: HetBioYield
        real :: AutoBioYield
        real :: FracBioPartProd
        real :: NCODBioMassRatio
        real :: NCODPartProdBioMassRatio

        real :: MaxSpecGrowthRateHetBio
        real :: HalfSatCoefHetBio
        real :: OxygHalfSatCoefHetBio
        real :: NitHalfSatCoefDenHetBio
        real :: DecayCoefHetBio
        real :: CorFacAnoxGrowthHetBio
        real :: AmmonifRate
        real :: MaxSpecHydroRate
        real :: HalfSatCoefHydroSlowlyBioSub
        real :: CorFacAnoxCond
        real :: MaxSpecGrowthRateAutBio
        real :: AmmoniaHalfSatCoefAutBio
        real :: OxyHalfSatCoefAutBio
        real :: DecayCoefAutBio
        
        real :: DTDay
 
!------------------------------------------------------------------------
!State variables
        SolInertOrgMat                 = Me%PropIndex%SolInertOrgMat
        ReadilyBioSub                  = Me%PropIndex%ReadilyBioSub
        PartInertOrgMar                = Me%PropIndex%PartInertOrgMar 
        SlowlyBioSub                   = Me%PropIndex%SlowlyBioSub
        HetBio                         = Me%PropIndex%HetBio
        AutBio                         = Me%PropIndex%AutBio
        PartProd                       = Me%PropIndex%PartProd
        Oxygen                         = Me%PropIndex%Oxygen
        Nitrate                        = Me%PropIndex%Nitrate
        Ammonia                        = Me%PropIndex%Ammonia
        SolBioOrgNitrogen              = Me%PropIndex%SolBioOrgNitrogen
        PartBioOrgNitrogen             = Me%PropIndex%PartBioOrgNitrogen 
        Alkalinity                     = Me%PropIndex%Alkalinity 

!Stecheo Parameters
        HetBioYield                    = Me%StechParameters%HetBioYield
        AutoBioYield                   = Me%StechParameters%AutoBioYield
        FracBioPartProd                = Me%StechParameters%FracBioPartProd
        NCODBioMassRatio               = Me%StechParameters%NCODBioMassRatio
        NCODPartProdBioMassRatio       = Me%StechParameters%NCODPartProdBioMassRatio

!Kinetic Parameters
        MaxSpecGrowthRateHetBio        = Me%KineticParameters%MaxSpecGrowthRateHetBio
        HalfSatCoefHetBio              = Me%KineticParameters%HalfSatCoefHetBio
        OxygHalfSatCoefHetBio          = Me%KineticParameters%OxygHalfSatCoefHetBio
        NitHalfSatCoefDenHetBio        = Me%KineticParameters%NitHalfSatCoefDenHetBio
        DecayCoefHetBio                = Me%KineticParameters%DecayCoefHetBio
        CorFacAnoxGrowthHetBio         = Me%KineticParameters%CorFacAnoxGrowthHetBio
        AmmonifRate                    = Me%KineticParameters%AmmonifRate
        MaxSpecHydroRate               = Me%KineticParameters%MaxSpecHydroRate
        HalfSatCoefHydroSlowlyBioSub   = Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub
        CorFacAnoxCond                 = Me%KineticParameters%CorFacAnoxCond
        MaxSpecGrowthRateAutBio        = Me%KineticParameters%MaxSpecGrowthRateAutBio
        AmmoniaHalfSatCoefAutBio       = Me%KineticParameters%AmmoniaHalfSatCoefAutBio
        OxyHalfSatCoefAutBio           = Me%KineticParameters%OxyHalfSatCoefAutBio
        DecayCoefAutBio                = Me%KineticParameters%DecayCoefAutBio

        DTDay                          = Me%DTDay
              
 !Calculation of system coeficients---------------------------------------
         Me%Matrix(PartInertOrgMar, PartInertOrgMar) = 1.0
       
 !Independent term
        Me%IndTerm(PartInertOrgMar) = Me%ExternalVar%Mass(PartInertOrgMar, index)

 !------------------------------------------------------------------------

    end subroutine WWTPQPartInertOrgMar 

    subroutine WWTPQSlowlyBioSub (index)
    
    !Subroutine to calculate the propriety SlowlyBioSub

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
    
    !State variables
        integer :: SolInertOrgMat
        integer :: ReadilyBioSub
        integer :: PartInertOrgMar
        integer :: SlowlyBioSub
        integer :: HetBio
        integer :: AutBio
        integer :: PartProd
        integer :: Oxygen
        integer :: Nitrate
        integer :: Ammonia
        integer :: SolBioOrgNitrogen
        integer :: PartBioOrgNitrogen
        integer :: Alkalinity

        real :: HetBioYield
        real :: AutoBioYield
        real :: FracBioPartProd
        real :: NCODBioMassRatio
        real :: NCODPartProdBioMassRatio

        real :: MaxSpecGrowthRateHetBio
        real :: HalfSatCoefHetBio
        real :: OxygHalfSatCoefHetBio
        real :: NitHalfSatCoefDenHetBio
        real :: DecayCoefHetBio
        real :: CorFacAnoxGrowthHetBio
        real :: AmmonifRate
        real :: MaxSpecHydroRate
        real :: HalfSatCoefHydroSlowlyBioSub
        real :: CorFacAnoxCond
        real :: MaxSpecGrowthRateAutBio
        real :: AmmoniaHalfSatCoefAutBio
        real :: OxyHalfSatCoefAutBio
        real :: DecayCoefAutBio
 
        real :: DTDay
        real :: SlowlyBioSubRateHetBio             = null_real
        real :: SlowlyBioSubRateAutBio             = null_real
        
        real :: SlowlyBioSubRateHetBioDecay        = null_real
        real :: SlowlyBioSubRateHetBioHydrolysis   = null_real
        real :: SlowlyBioSubRateAutBioDecay        = null_real

!------------------------------------------------------------------------

!State variables
        SolInertOrgMat                 = Me%PropIndex%SolInertOrgMat
        ReadilyBioSub                  = Me%PropIndex%ReadilyBioSub
        PartInertOrgMar                = Me%PropIndex%PartInertOrgMar 
        SlowlyBioSub                   = Me%PropIndex%SlowlyBioSub
        HetBio                         = Me%PropIndex%HetBio
        AutBio                         = Me%PropIndex%AutBio
        PartProd                       = Me%PropIndex%PartProd
        Oxygen                         = Me%PropIndex%Oxygen
        Nitrate                        = Me%PropIndex%Nitrate
        Ammonia                        = Me%PropIndex%Ammonia
        SolBioOrgNitrogen              = Me%PropIndex%SolBioOrgNitrogen
        PartBioOrgNitrogen             = Me%PropIndex%PartBioOrgNitrogen 
        Alkalinity                     = Me%PropIndex%Alkalinity 

!Stecheo Parameters
        HetBioYield                    = Me%StechParameters%HetBioYield
        AutoBioYield                   = Me%StechParameters%AutoBioYield
        FracBioPartProd                = Me%StechParameters%FracBioPartProd
        NCODBioMassRatio               = Me%StechParameters%NCODBioMassRatio
        NCODPartProdBioMassRatio       = Me%StechParameters%NCODPartProdBioMassRatio

!Kinetic Parameters
        MaxSpecGrowthRateHetBio        = Me%KineticParameters%MaxSpecGrowthRateHetBio
        HalfSatCoefHetBio              = Me%KineticParameters%HalfSatCoefHetBio
        OxygHalfSatCoefHetBio          = Me%KineticParameters%OxygHalfSatCoefHetBio
        NitHalfSatCoefDenHetBio        = Me%KineticParameters%NitHalfSatCoefDenHetBio
        DecayCoefHetBio                = Me%KineticParameters%DecayCoefHetBio
        CorFacAnoxGrowthHetBio         = Me%KineticParameters%CorFacAnoxGrowthHetBio
        AmmonifRate                    = Me%KineticParameters%AmmonifRate
        MaxSpecHydroRate               = Me%KineticParameters%MaxSpecHydroRate
        HalfSatCoefHydroSlowlyBioSub   = Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub
        CorFacAnoxCond                 = Me%KineticParameters%CorFacAnoxCond
        MaxSpecGrowthRateAutBio        = Me%KineticParameters%MaxSpecGrowthRateAutBio
        AmmoniaHalfSatCoefAutBio       = Me%KineticParameters%AmmoniaHalfSatCoefAutBio
        OxyHalfSatCoefAutBio           = Me%KineticParameters%OxyHalfSatCoefAutBio
        DecayCoefAutBio                = Me%KineticParameters%DecayCoefAutBio

        DTDay                          = Me%DTDay


!SlowlyBioSubRateHetBio, SlowlyBioSub Rate in function of HetBio, D-1

        SlowlyBioSubRateHetBioDecay = (1-FracBioPartProd)*DecayCoefHetBio  
           
           if (Me%ExternalVar%Mass(HetBio, index) == 0) then
           
           SlowlyBioSubRateHetBioHydrolysis=0
           
           else
           
        SlowlyBioSubRateHetBioHydrolysis = -MaxSpecHydroRate*((Me%ExternalVar%Mass(SlowlyBioSub, index)                  &
        /Me%ExternalVar%Mass(HetBio, index))/                                                                            &
           (HalfSatCoefHydroSlowlyBioSub+Me%ExternalVar%Mass(SlowlyBioSub, index)/Me%ExternalVar%Mass(HetBio, index)))   &
           *((Me%ExternalVar%Mass(Oxygen, index)/(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))             &
           +CorFacAnoxCond*(OxygHalfSatCoefHetBio                                                                        &
           /(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))                                                  &
           *(Me%ExternalVar%Mass(Nitrate, index)/(NitHalfSatCoefDenHetBio+Me%ExternalVar%Mass(Nitrate, index))))

          end if

           SlowlyBioSubRateHetBio = SlowlyBioSubRateHetBioDecay + SlowlyBioSubRateHetBioHydrolysis
           

!SlowlyBioSubRateAutBio, SlowlyBioSub Rate in function of AutBio, D-1   
        
           SlowlyBioSubRateAutBioDecay=(1-FracBioPartProd)*DecayCoefAutBio
                   
           SlowlyBioSubRateAutBio = SlowlyBioSubRateAutBioDecay
                   
 !Calculation of system coeficients---------------------------------------
         Me%Matrix(SlowlyBioSub, HetBio) = DTDay * (-SlowlyBioSubRateHetBio) 
         Me%Matrix(SlowlyBioSub, AutBio) = DTDay * (-SlowlyBioSubRateAutBio) 
         Me%Matrix(SlowlyBioSub, SlowlyBioSub) = 1.0

 !Independent term
        Me%IndTerm(SlowlyBioSub) = Me%ExternalVar%Mass(SlowlyBioSub, index)
   
 !------------------------------------------------------------------------    
     
    end subroutine WWTPQSlowlyBioSub

    subroutine WWTPQHetBio (index)
          
!Subroutine to calculate HetBio

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------

        integer :: SolInertOrgMat
        integer :: ReadilyBioSub
        integer :: PartInertOrgMar
        integer :: SlowlyBioSub
        integer :: HetBio
        integer :: AutBio
        integer :: PartProd
        integer :: Oxygen
        integer :: Nitrate
        integer :: Ammonia
        integer :: SolBioOrgNitrogen
        integer :: PartBioOrgNitrogen
        integer :: Alkalinity

        real :: HetBioYield
        real :: AutoBioYield
        real :: FracBioPartProd
        real :: NCODBioMassRatio
        real :: NCODPartProdBioMassRatio

        real :: MaxSpecGrowthRateHetBio
        real :: HalfSatCoefHetBio
        real :: OxygHalfSatCoefHetBio
        real :: NitHalfSatCoefDenHetBio
        real :: DecayCoefHetBio
        real :: CorFacAnoxGrowthHetBio
        real :: AmmonifRate
        real :: MaxSpecHydroRate
        real :: HalfSatCoefHydroSlowlyBioSub
        real :: CorFacAnoxCond
        real :: MaxSpecGrowthRateAutBio
        real :: AmmoniaHalfSatCoefAutBio
        real :: OxyHalfSatCoefAutBio
        real :: DecayCoefAutBio
        
        real :: DTDay
        real :: HetBioRateHetBio               = null_real
        
        real :: HetBioRateHetBioAerGrowth      = null_real
        real :: HetBioRateHetBioAnoGrowth      = null_real
        real :: HetBioRateHetBioDecay          = null_real

     !------------------------------------------------------------------------

!State variables
        SolInertOrgMat                 = Me%PropIndex%SolInertOrgMat
        ReadilyBioSub                  = Me%PropIndex%ReadilyBioSub
        PartInertOrgMar                = Me%PropIndex%PartInertOrgMar 
        SlowlyBioSub                   = Me%PropIndex%SlowlyBioSub
        HetBio                         = Me%PropIndex%HetBio
        AutBio                         = Me%PropIndex%AutBio
        PartProd                       = Me%PropIndex%PartProd
        Oxygen                         = Me%PropIndex%Oxygen
        Nitrate                        = Me%PropIndex%Nitrate
        Ammonia                        = Me%PropIndex%Ammonia
        SolBioOrgNitrogen              = Me%PropIndex%SolBioOrgNitrogen
        PartBioOrgNitrogen             = Me%PropIndex%PartBioOrgNitrogen 
        Alkalinity                     = Me%PropIndex%Alkalinity 

!Stecheo Parameters
        HetBioYield                    = Me%StechParameters%HetBioYield
        AutoBioYield                   = Me%StechParameters%AutoBioYield
        FracBioPartProd                = Me%StechParameters%FracBioPartProd
        NCODBioMassRatio               = Me%StechParameters%NCODBioMassRatio
        NCODPartProdBioMassRatio       = Me%StechParameters%NCODPartProdBioMassRatio

!Kinetic Parameters
        MaxSpecGrowthRateHetBio        = Me%KineticParameters%MaxSpecGrowthRateHetBio
        HalfSatCoefHetBio              = Me%KineticParameters%HalfSatCoefHetBio
        OxygHalfSatCoefHetBio          = Me%KineticParameters%OxygHalfSatCoefHetBio
        NitHalfSatCoefDenHetBio        = Me%KineticParameters%NitHalfSatCoefDenHetBio
        DecayCoefHetBio                = Me%KineticParameters%DecayCoefHetBio
        CorFacAnoxGrowthHetBio         = Me%KineticParameters%CorFacAnoxGrowthHetBio
        AmmonifRate                    = Me%KineticParameters%AmmonifRate
        MaxSpecHydroRate               = Me%KineticParameters%MaxSpecHydroRate
        HalfSatCoefHydroSlowlyBioSub   = Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub
        CorFacAnoxCond                 = Me%KineticParameters%CorFacAnoxCond
        MaxSpecGrowthRateAutBio        = Me%KineticParameters%MaxSpecGrowthRateAutBio
        AmmoniaHalfSatCoefAutBio       = Me%KineticParameters%AmmoniaHalfSatCoefAutBio
        OxyHalfSatCoefAutBio           = Me%KineticParameters%OxyHalfSatCoefAutBio
        DecayCoefAutBio                = Me%KineticParameters%DecayCoefAutBio

        DTDay                          = Me%DTDay

!HetBioRateHetBio, HetBio Rate in function of HetBio, D-1

        HetBioRateHetBioAerGrowth=MaxSpecGrowthRateHetBio*(Me%ExternalVar%Mass(ReadilyBioSub, index)    &
           /(HalfSatCoefHetBio+Me%ExternalVar%Mass(ReadilyBioSub, index)))                              &
           *(Me%ExternalVar%Mass(Oxygen, index)                                                         &
           /(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index))) 
   
  
        HetBioRateHetBioAnoGrowth=MaxSpecGrowthRateHetBio*(Me%ExternalVar%Mass(ReadilyBioSub, index)    & 
           /(HalfSatCoefHetBio+                                                                         &
           Me%ExternalVar%Mass(ReadilyBioSub, index)))*(OxygHalfSatCoefHetBio                           &
           /(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))                                 &
           *(Me%ExternalVar%Mass(Nitrate, index)                                                        &
           /(NitHalfSatCoefDenHetBio+Me%ExternalVar%Mass(Nitrate, index)))*CorFacAnoxGrowthHetBio

                   
        HetBioRateHetBioDecay= -DecayCoefHetBio      


        HetBioRateHetBio =HetBioRateHetBioAerGrowth+HetBioRateHetBioAnoGrowth                           &
           +HetBioRateHetBioDecay                                   
                
 !Calculation of system coeficients---------------------------------------
         Me%Matrix(HetBio,HetBio ) = DTDay * (-HetBioRateHetBio) + 1.0

 !Independent term
         Me%IndTerm(HetBio) = Me%ExternalVar%Mass(HetBio, index)
    
 !------------------------------------------------------------------------    

    end subroutine WWTPQHetBio
   
   subroutine WWTPQAutBio (index)
    
    !Subroutine to calculate the propriety AutBio

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------

        integer :: SolInertOrgMat
        integer :: ReadilyBioSub
        integer :: PartInertOrgMar
        integer :: SlowlyBioSub
        integer :: HetBio
        integer :: AutBio
        integer :: PartProd
        integer :: Oxygen
        integer :: Nitrate
        integer :: Ammonia
        integer :: SolBioOrgNitrogen
        integer :: PartBioOrgNitrogen
        integer :: Alkalinity
        
        real :: HetBioYield
        real :: AutoBioYield
        real :: FracBioPartProd
        real :: NCODBioMassRatio
        real :: NCODPartProdBioMassRatio

        real :: MaxSpecGrowthRateHetBio
        real :: HalfSatCoefHetBio
        real :: OxygHalfSatCoefHetBio
        real :: NitHalfSatCoefDenHetBio
        real :: DecayCoefHetBio
        real :: CorFacAnoxGrowthHetBio
        real :: AmmonifRate
        real :: MaxSpecHydroRate
        real :: HalfSatCoefHydroSlowlyBioSub
        real :: CorFacAnoxCond
        real :: MaxSpecGrowthRateAutBio
        real :: AmmoniaHalfSatCoefAutBio
        real :: OxyHalfSatCoefAutBio
        real :: DecayCoefAutBio
 
        real :: DTDay
        real :: AutBioRateAutBio             = null_real
        
        real :: AutBioRateAutBioAerGrowth    = null_real
        real :: AutBioRateAutBioDecay        = null_real
        
  
!------------------------------------------------------------------------
!State variables 
        SolInertOrgMat                 = Me%PropIndex%SolInertOrgMat
        ReadilyBioSub                  = Me%PropIndex%ReadilyBioSub
        PartInertOrgMar                = Me%PropIndex%PartInertOrgMar 
        SlowlyBioSub                   = Me%PropIndex%SlowlyBioSub
        HetBio                         = Me%PropIndex%HetBio
        AutBio                         = Me%PropIndex%AutBio
        PartProd                       = Me%PropIndex%PartProd
        Oxygen                         = Me%PropIndex%Oxygen
        Nitrate                        = Me%PropIndex%Nitrate
        Ammonia                        = Me%PropIndex%Ammonia
        SolBioOrgNitrogen              = Me%PropIndex%SolBioOrgNitrogen
        PartBioOrgNitrogen             = Me%PropIndex%PartBioOrgNitrogen 
        Alkalinity                     = Me%PropIndex%Alkalinity 

!Stecheo Parameters
        HetBioYield                    = Me%StechParameters%HetBioYield
        AutoBioYield                   = Me%StechParameters%AutoBioYield
        FracBioPartProd                = Me%StechParameters%FracBioPartProd
        NCODBioMassRatio               = Me%StechParameters%NCODBioMassRatio
        NCODPartProdBioMassRatio       = Me%StechParameters%NCODPartProdBioMassRatio

!Kinetic Parameters
        MaxSpecGrowthRateHetBio        = Me%KineticParameters%MaxSpecGrowthRateHetBio
        HalfSatCoefHetBio              = Me%KineticParameters%HalfSatCoefHetBio
        OxygHalfSatCoefHetBio          = Me%KineticParameters%OxygHalfSatCoefHetBio
        NitHalfSatCoefDenHetBio        = Me%KineticParameters%NitHalfSatCoefDenHetBio
        DecayCoefHetBio                = Me%KineticParameters%DecayCoefHetBio
        CorFacAnoxGrowthHetBio         = Me%KineticParameters%CorFacAnoxGrowthHetBio
        AmmonifRate                    = Me%KineticParameters%AmmonifRate
        MaxSpecHydroRate               = Me%KineticParameters%MaxSpecHydroRate
        HalfSatCoefHydroSlowlyBioSub   = Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub
        CorFacAnoxCond                 = Me%KineticParameters%CorFacAnoxCond
        MaxSpecGrowthRateAutBio        = Me%KineticParameters%MaxSpecGrowthRateAutBio
        AmmoniaHalfSatCoefAutBio       = Me%KineticParameters%AmmoniaHalfSatCoefAutBio
        OxyHalfSatCoefAutBio           = Me%KineticParameters%OxyHalfSatCoefAutBio
        DecayCoefAutBio                = Me%KineticParameters%DecayCoefAutBio

        DTDay                          = Me%DTDay

!AutBioRateAutBio, AutBioRateAutBio in function of AutBio, D-1

         AutBioRateAutBioAerGrowth= MaxSpecGrowthRateAutBio*(Me%ExternalVar%Mass(Ammonia, index)   &
           /(AmmoniaHalfSatCoefAutBio+Me%ExternalVar%Mass(Ammonia, index)))  &
           *(Me%ExternalVar%Mass(Oxygen, index)/(OxyHalfSatCoefAutBio     &
           +Me%ExternalVar%Mass(Oxygen, index)))
           
        AutBioRateAutBioDecay=-DecayCoefAutBio  

        AutBioRateAutBio = AutBioRateAutBioAerGrowth+AutBioRateAutBioDecay                            
           
 !Calculation of system coeficients---------------------------------------
         Me%Matrix(AutBio,AutBio) = DTDay * (-AutBioRateAutBio) + 1.0

 !Independent term
         Me%IndTerm(AutBio) = Me%ExternalVar%Mass(AutBio, index)

 !------------------------------------------------------------------------ 

    end subroutine WWTPQAutBio

   subroutine WWTPQPartProd (index)
  
     !Subroutine to calculate the propriety PartProd

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
    
    !State variables
        integer :: SolInertOrgMat
        integer :: ReadilyBioSub
        integer :: PartInertOrgMar
        integer :: SlowlyBioSub
        integer :: HetBio
        integer :: AutBio
        integer :: PartProd
        integer :: Oxygen
        integer :: Nitrate
        integer :: Ammonia
        integer :: SolBioOrgNitrogen
        integer :: PartBioOrgNitrogen
        integer :: Alkalinity

        real :: HetBioYield
        real :: AutoBioYield
        real :: FracBioPartProd
        real :: NCODBioMassRatio
        real :: NCODPartProdBioMassRatio

        real :: MaxSpecGrowthRateHetBio
        real :: HalfSatCoefHetBio
        real :: OxygHalfSatCoefHetBio
        real :: NitHalfSatCoefDenHetBio
        real :: DecayCoefHetBio
        real :: CorFacAnoxGrowthHetBio
        real :: AmmonifRate
        real :: MaxSpecHydroRate
        real :: HalfSatCoefHydroSlowlyBioSub
        real :: CorFacAnoxCond
        real :: MaxSpecGrowthRateAutBio
        real :: AmmoniaHalfSatCoefAutBio
        real :: OxyHalfSatCoefAutBio
        real :: DecayCoefAutBio
 
        real    :: DTDay
        real    :: PartProdRateHetBio      = null_real
        real    :: PartProdRateAutBio      = null_real
        
        real    :: PartProdRateHetBioDecay      = null_real
        real    :: PartProdRateAutBioDecay      = null_real
        
!------------------------------------------------------------------------

!State variables
        SolInertOrgMat                 = Me%PropIndex%SolInertOrgMat
        ReadilyBioSub                  = Me%PropIndex%ReadilyBioSub
        PartInertOrgMar                = Me%PropIndex%PartInertOrgMar 
        SlowlyBioSub                   = Me%PropIndex%SlowlyBioSub
        HetBio                         = Me%PropIndex%HetBio
        AutBio                         = Me%PropIndex%AutBio
        PartProd                       = Me%PropIndex%PartProd
        Oxygen                         = Me%PropIndex%Oxygen
        Nitrate                        = Me%PropIndex%Nitrate
        Ammonia                        = Me%PropIndex%Ammonia
        SolBioOrgNitrogen              = Me%PropIndex%SolBioOrgNitrogen
        PartBioOrgNitrogen             = Me%PropIndex%PartBioOrgNitrogen 
        Alkalinity                     = Me%PropIndex%Alkalinity 

!Stecheo Parameters
        HetBioYield                    = Me%StechParameters%HetBioYield
        AutoBioYield                   = Me%StechParameters%AutoBioYield
        FracBioPartProd                = Me%StechParameters%FracBioPartProd
        NCODBioMassRatio               = Me%StechParameters%NCODBioMassRatio
        NCODPartProdBioMassRatio       = Me%StechParameters%NCODPartProdBioMassRatio

!Kinetic Parameters
        MaxSpecGrowthRateHetBio        = Me%KineticParameters%MaxSpecGrowthRateHetBio
        HalfSatCoefHetBio              = Me%KineticParameters%HalfSatCoefHetBio
        OxygHalfSatCoefHetBio          = Me%KineticParameters%OxygHalfSatCoefHetBio
        NitHalfSatCoefDenHetBio        = Me%KineticParameters%NitHalfSatCoefDenHetBio
        DecayCoefHetBio                = Me%KineticParameters%DecayCoefHetBio
        CorFacAnoxGrowthHetBio         = Me%KineticParameters%CorFacAnoxGrowthHetBio
        AmmonifRate                    = Me%KineticParameters%AmmonifRate
        MaxSpecHydroRate               = Me%KineticParameters%MaxSpecHydroRate
        HalfSatCoefHydroSlowlyBioSub   = Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub
        CorFacAnoxCond                 = Me%KineticParameters%CorFacAnoxCond
        MaxSpecGrowthRateAutBio        = Me%KineticParameters%MaxSpecGrowthRateAutBio
        AmmoniaHalfSatCoefAutBio       = Me%KineticParameters%AmmoniaHalfSatCoefAutBio
        OxyHalfSatCoefAutBio           = Me%KineticParameters%OxyHalfSatCoefAutBio
        DecayCoefAutBio                = Me%KineticParameters%DecayCoefAutBio

        DTDay                          = Me%DTDay


!PartProdRateHetBio, PartProd Rate in function of HetBio, D-1

           PartProdRateHetBioDecay=FracBioPartProd*DecayCoefHetBio
   
           PartProdRateHetBio = PartProdRateHetBioDecay

!PartProdRateAutBio, PartProd Rate in function of AutBio, D-1   
                
           PartProdRateAutBioDecay = FracBioPartProd*DecayCoefAutBio
                
           PartProdRateAutBio = PartProdRateAutBioDecay 
           
  !Calculation of system coeficients---------------------------------------
         Me%Matrix(PartProd, HetBio) = DTDay * (-PartProdRateHetBio) 
         Me%Matrix(PartProd, AutBio) = DTDay * (-PartProdRateAutBio) 
         Me%Matrix(PartProd, PartProd) = 1.0

 !Independent term
        Me%IndTerm(PartProd) = Me%ExternalVar%Mass(PartProd, index)

 !------------------------------------------------------------------------ 

    end subroutine WWTPQPartProd

  subroutine WWTPQOxygen (index)
     
    !Subroutine to calculate the propriety Oxygen

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
    
    !State variables
        integer :: SolInertOrgMat
        integer :: ReadilyBioSub
        integer :: PartInertOrgMar
        integer :: SlowlyBioSub
        integer :: HetBio
        integer :: AutBio
        integer :: PartProd
        integer :: Oxygen
        integer :: Nitrate
        integer :: Ammonia
        integer :: SolBioOrgNitrogen
        integer :: PartBioOrgNitrogen
        integer :: Alkalinity

        real :: HetBioYield
        real :: AutoBioYield
        real :: FracBioPartProd
        real :: NCODBioMassRatio
        real :: NCODPartProdBioMassRatio

        real :: MaxSpecGrowthRateHetBio
        real :: HalfSatCoefHetBio
        real :: OxygHalfSatCoefHetBio
        real :: NitHalfSatCoefDenHetBio
        real :: DecayCoefHetBio
        real :: CorFacAnoxGrowthHetBio
        real :: AmmonifRate
        real :: MaxSpecHydroRate
        real :: HalfSatCoefHydroSlowlyBioSub
        real :: CorFacAnoxCond
        real :: MaxSpecGrowthRateAutBio
        real :: AmmoniaHalfSatCoefAutBio
        real :: OxyHalfSatCoefAutBio
        real :: DecayCoefAutBio
 
        real :: DTDay
        real :: OxygenRateHetBio               = null_real
        real :: OxygenRateAutBio               = null_real
        real :: OxygenRateHetBioAerGrowth      = null_real
        real :: OxygenRateAutBioAerGrowth      = null_real

!------------------------------------------------------------------------

!State variables
        SolInertOrgMat                 = Me%PropIndex%SolInertOrgMat
        ReadilyBioSub                  = Me%PropIndex%ReadilyBioSub
        PartInertOrgMar                = Me%PropIndex%PartInertOrgMar 
        SlowlyBioSub                   = Me%PropIndex%SlowlyBioSub
        HetBio                         = Me%PropIndex%HetBio
        AutBio                         = Me%PropIndex%AutBio
        PartProd                       = Me%PropIndex%PartProd
        Oxygen                         = Me%PropIndex%Oxygen
        Nitrate                        = Me%PropIndex%Nitrate
        Ammonia                        = Me%PropIndex%Ammonia
        SolBioOrgNitrogen              = Me%PropIndex%SolBioOrgNitrogen
        PartBioOrgNitrogen             = Me%PropIndex%PartBioOrgNitrogen 
        Alkalinity                     = Me%PropIndex%Alkalinity 

!Stecheo Parameters
        HetBioYield                    = Me%StechParameters%HetBioYield
        AutoBioYield                   = Me%StechParameters%AutoBioYield
        FracBioPartProd                = Me%StechParameters%FracBioPartProd
        NCODBioMassRatio               = Me%StechParameters%NCODBioMassRatio
        NCODPartProdBioMassRatio       = Me%StechParameters%NCODPartProdBioMassRatio

!Kinetic Parameters
        MaxSpecGrowthRateHetBio        = Me%KineticParameters%MaxSpecGrowthRateHetBio
        HalfSatCoefHetBio              = Me%KineticParameters%HalfSatCoefHetBio
        OxygHalfSatCoefHetBio          = Me%KineticParameters%OxygHalfSatCoefHetBio
        NitHalfSatCoefDenHetBio        = Me%KineticParameters%NitHalfSatCoefDenHetBio
        DecayCoefHetBio                = Me%KineticParameters%DecayCoefHetBio
        CorFacAnoxGrowthHetBio         = Me%KineticParameters%CorFacAnoxGrowthHetBio
        AmmonifRate                    = Me%KineticParameters%AmmonifRate
        MaxSpecHydroRate               = Me%KineticParameters%MaxSpecHydroRate
        HalfSatCoefHydroSlowlyBioSub   = Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub
        CorFacAnoxCond                 = Me%KineticParameters%CorFacAnoxCond
        MaxSpecGrowthRateAutBio        = Me%KineticParameters%MaxSpecGrowthRateAutBio
        AmmoniaHalfSatCoefAutBio       = Me%KineticParameters%AmmoniaHalfSatCoefAutBio
        OxyHalfSatCoefAutBio           = Me%KineticParameters%OxyHalfSatCoefAutBio
        DecayCoefAutBio                = Me%KineticParameters%DecayCoefAutBio

        DTDay                          = Me%DTDay


!OxygenRateHetBio, Oxygen Rate in fucntion of HetBio, D-1

        OxygenRateHetBioAerGrowth=-((1-HetBioYield)/HetBioYield)                                                            &
                   *MaxSpecGrowthRateHetBio*                                                                                &
           (Me%ExternalVar%Mass(ReadilyBioSub, index)/(HalfSatCoefHetBio+Me%ExternalVar%Mass(ReadilyBioSub, index)))        &
           *(Me%ExternalVar%Mass(Oxygen, index)/(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index))) 
       
           OxygenRateHetBio = OxygenRateHetBioAerGrowth        
          
!OxygenRateAutBio, Oxygen Rate in function of AutBio, D-1  
           
        OxygenRateAutBioAerGrowth      = -((4.57-AutoBioYield)/AutoBioYield)*MaxSpecGrowthRateAutBio*                       &
           (Me%ExternalVar%Mass(Ammonia, index)/(AmmoniaHalfSatCoefAutBio+Me%ExternalVar%Mass(Ammonia, index)))             &
           *(Me%ExternalVar%Mass(Oxygen, index)/(OxyHalfSatCoefAutBio+Me%ExternalVar%Mass(Oxygen, index)))
           
           OxygenRateAutBio =  OxygenRateAutBioAerGrowth 
           
 !Calculation of system coeficients---------------------------------------
         Me%Matrix(Oxygen, HetBio) = DTDay * (-OxygenRateHetBio) 
         Me%Matrix(Oxygen, AutBio) = DTDay * (-OxygenRateAutBio)
         Me%Matrix(Oxygen, Oxygen) = 1.0 

 !Independent term
        Me%IndTerm(Oxygen) = Me%ExternalVar%Mass(Oxygen, index)

 !------------------------------------------------------------------------ 

    end subroutine WWTPQOxygen

subroutine WWTPQNitrate (index)
      
    !Subroutine to calculate the propriety Nitrate

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
    
    !State variables
        integer :: SolInertOrgMat
        integer :: ReadilyBioSub
        integer :: PartInertOrgMar
        integer :: SlowlyBioSub
        integer :: HetBio
        integer :: AutBio
        integer :: PartProd
        integer :: Oxygen
        integer :: Nitrate
        integer :: Ammonia
        integer :: SolBioOrgNitrogen
        integer :: PartBioOrgNitrogen
        integer :: Alkalinity

        real :: HetBioYield
        real :: AutoBioYield
        real :: FracBioPartProd
        real :: NCODBioMassRatio
        real :: NCODPartProdBioMassRatio

        real :: MaxSpecGrowthRateHetBio
        real :: HalfSatCoefHetBio
        real :: OxygHalfSatCoefHetBio
        real :: NitHalfSatCoefDenHetBio
        real :: DecayCoefHetBio
        real :: CorFacAnoxGrowthHetBio
        real :: AmmonifRate
        real :: MaxSpecHydroRate
        real :: HalfSatCoefHydroSlowlyBioSub
        real :: CorFacAnoxCond
        real :: MaxSpecGrowthRateAutBio
        real :: AmmoniaHalfSatCoefAutBio
        real :: OxyHalfSatCoefAutBio
        real :: DecayCoefAutBio
 
        real    :: DTDay
        real    :: NitrateRateHetBio      = null_real
        real    :: NitrateRateAutBio      = null_real
        
        real    :: NitrateRateHetBioAnoGrowth      = null_real
        real    :: NitrateRateAutBioAerGrowth      = null_real

!------------------------------------------------------------------------

!State variables     
        SolInertOrgMat                 = Me%PropIndex%SolInertOrgMat
        ReadilyBioSub                  = Me%PropIndex%ReadilyBioSub
        PartInertOrgMar                = Me%PropIndex%PartInertOrgMar 
        SlowlyBioSub                   = Me%PropIndex%SlowlyBioSub
        HetBio                         = Me%PropIndex%HetBio
        AutBio                         = Me%PropIndex%AutBio
        PartProd                       = Me%PropIndex%PartProd
        Oxygen                         = Me%PropIndex%Oxygen
        Nitrate                        = Me%PropIndex%Nitrate
        Ammonia                        = Me%PropIndex%Ammonia
        SolBioOrgNitrogen              = Me%PropIndex%SolBioOrgNitrogen
        PartBioOrgNitrogen             = Me%PropIndex%PartBioOrgNitrogen 
        Alkalinity                     = Me%PropIndex%Alkalinity 

!Stecheo Parameters
        HetBioYield                    = Me%StechParameters%HetBioYield
        AutoBioYield                   = Me%StechParameters%AutoBioYield
        FracBioPartProd                = Me%StechParameters%FracBioPartProd
        NCODBioMassRatio               = Me%StechParameters%NCODBioMassRatio
        NCODPartProdBioMassRatio       = Me%StechParameters%NCODPartProdBioMassRatio

!Kinetic Parameters
        MaxSpecGrowthRateHetBio        = Me%KineticParameters%MaxSpecGrowthRateHetBio
        HalfSatCoefHetBio              = Me%KineticParameters%HalfSatCoefHetBio
        OxygHalfSatCoefHetBio          = Me%KineticParameters%OxygHalfSatCoefHetBio
        NitHalfSatCoefDenHetBio        = Me%KineticParameters%NitHalfSatCoefDenHetBio
        DecayCoefHetBio                = Me%KineticParameters%DecayCoefHetBio
        CorFacAnoxGrowthHetBio         = Me%KineticParameters%CorFacAnoxGrowthHetBio
        AmmonifRate                    = Me%KineticParameters%AmmonifRate
        MaxSpecHydroRate               = Me%KineticParameters%MaxSpecHydroRate
        HalfSatCoefHydroSlowlyBioSub   = Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub
        CorFacAnoxCond                 = Me%KineticParameters%CorFacAnoxCond
        MaxSpecGrowthRateAutBio        = Me%KineticParameters%MaxSpecGrowthRateAutBio
        AmmoniaHalfSatCoefAutBio       = Me%KineticParameters%AmmoniaHalfSatCoefAutBio
        OxyHalfSatCoefAutBio           = Me%KineticParameters%OxyHalfSatCoefAutBio
        DecayCoefAutBio                = Me%KineticParameters%DecayCoefAutBio

        DTDay                          = Me%DTDay

!NitrateRateHetBio, Nitrate Rate in realtion of HetBio, units
      
         NitrateRateHetBioAnoGrowth      = -((1.-HetBioYield)/(2.86*HetBioYield))*MaxSpecGrowthRateHetBio*              &
          (Me%ExternalVar%Mass(ReadilyBioSub, index)                                                                    &
           /(HalfSatCoefHetBio+Me%ExternalVar%Mass(ReadilyBioSub, index)))*(OxygHalfSatCoefHetBio                       &
           /(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))                                                 &
           *(Me%ExternalVar%Mass(Nitrate, index)/(NitHalfSatCoefDenHetBio+Me%ExternalVar%Mass(Nitrate, index)))         &
           *CorFacAnoxGrowthHetBio  
           
           NitrateRateHetBio = NitrateRateHetBioAnoGrowth  

!NitrateRateAutBio, Nitrate Rate in realtion of AutBio, units   
           
         NitrateRateAutBioAerGrowth = (1./AutoBioYield)*MaxSpecGrowthRateAutBio*(Me%ExternalVar%Mass(Ammonia, index)   &
           /(AmmoniaHalfSatCoefAutBio+Me%ExternalVar%Mass(Ammonia, index)))                                             &
           *(Me%ExternalVar%Mass(Oxygen, index)/(OxyHalfSatCoefAutBio                                                   &
           +Me%ExternalVar%Mass(Oxygen, index)))
           
           NitrateRateAutBio = NitrateRateAutBioAerGrowth  
 
 !Calculation of system coeficients---------------------------------------
         Me%Matrix(Nitrate, HetBio) = DTDay * (-NitrateRateHetBio) 
         Me%Matrix(Nitrate, AutBio) = DTDay * (-NitrateRateAutBio)
         Me%Matrix(Nitrate, Nitrate) = 1.0 

 !Independent term
        Me%IndTerm(Nitrate) = Me%ExternalVar%Mass(Nitrate, index)

 !------------------------------------------------------------------------
 
    end subroutine WWTPQNitrate

    subroutine WWTPQAmmonia (index)
  
        !Subroutine to calculate the propriety Ammonia

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
    
    !State variables
        integer :: SolInertOrgMat
        integer :: ReadilyBioSub
        integer :: PartInertOrgMar
        integer :: SlowlyBioSub
        integer :: HetBio
        integer :: AutBio
        integer :: PartProd
        integer :: Oxygen
        integer :: Nitrate
        integer :: Ammonia
        integer :: SolBioOrgNitrogen
        integer :: PartBioOrgNitrogen
        integer :: Alkalinity

        real :: HetBioYield
        real :: AutoBioYield
        real :: FracBioPartProd
        real :: NCODBioMassRatio
        real :: NCODPartProdBioMassRatio

        real :: MaxSpecGrowthRateHetBio
        real :: HalfSatCoefHetBio
        real :: OxygHalfSatCoefHetBio
        real :: NitHalfSatCoefDenHetBio
        real :: DecayCoefHetBio
        real :: CorFacAnoxGrowthHetBio
        real :: AmmonifRate
        real :: MaxSpecHydroRate
        real :: HalfSatCoefHydroSlowlyBioSub
        real :: CorFacAnoxCond
        real :: MaxSpecGrowthRateAutBio
        real :: AmmoniaHalfSatCoefAutBio
        real :: OxyHalfSatCoefAutBio
        real :: DecayCoefAutBio
 
        real :: DTDay
        real :: AmmoniaRateHetBio                 = null_real
        real :: AmmoniaRateAutBio                 = null_real
        
        real :: AmmoniaRateHetBioAerGrowth        = null_real
        real :: AmmoniaRateHetBioAnoGrowth        = null_real
        real :: AmmoniaRateHetBioAmmonification   = null_real
        
        real :: AmmoniaRateAutBioAerGrowth        = null_real
        
!------------------------------------------------------------------------

!State variables
       
        SolInertOrgMat                 = Me%PropIndex%SolInertOrgMat
        ReadilyBioSub                  = Me%PropIndex%ReadilyBioSub
        PartInertOrgMar                = Me%PropIndex%PartInertOrgMar 
        SlowlyBioSub                   = Me%PropIndex%SlowlyBioSub
        HetBio                         = Me%PropIndex%HetBio
        AutBio                         = Me%PropIndex%AutBio
        PartProd                       = Me%PropIndex%PartProd
        Oxygen                         = Me%PropIndex%Oxygen
        Nitrate                        = Me%PropIndex%Nitrate
        Ammonia                        = Me%PropIndex%Ammonia
        SolBioOrgNitrogen              = Me%PropIndex%SolBioOrgNitrogen
        PartBioOrgNitrogen             = Me%PropIndex%PartBioOrgNitrogen 
        Alkalinity                     = Me%PropIndex%Alkalinity 

!Stecheo Parameters
        HetBioYield                    = Me%StechParameters%HetBioYield
        AutoBioYield                   = Me%StechParameters%AutoBioYield
        FracBioPartProd                = Me%StechParameters%FracBioPartProd
        NCODBioMassRatio               = Me%StechParameters%NCODBioMassRatio
        NCODPartProdBioMassRatio       = Me%StechParameters%NCODPartProdBioMassRatio

!Kinetic Parameters
        MaxSpecGrowthRateHetBio        = Me%KineticParameters%MaxSpecGrowthRateHetBio
        HalfSatCoefHetBio              = Me%KineticParameters%HalfSatCoefHetBio
        OxygHalfSatCoefHetBio          = Me%KineticParameters%OxygHalfSatCoefHetBio
        NitHalfSatCoefDenHetBio        = Me%KineticParameters%NitHalfSatCoefDenHetBio
        DecayCoefHetBio                = Me%KineticParameters%DecayCoefHetBio
        CorFacAnoxGrowthHetBio         = Me%KineticParameters%CorFacAnoxGrowthHetBio
        AmmonifRate                    = Me%KineticParameters%AmmonifRate
        MaxSpecHydroRate               = Me%KineticParameters%MaxSpecHydroRate
        HalfSatCoefHydroSlowlyBioSub   = Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub
        CorFacAnoxCond                 = Me%KineticParameters%CorFacAnoxCond
        MaxSpecGrowthRateAutBio        = Me%KineticParameters%MaxSpecGrowthRateAutBio
        AmmoniaHalfSatCoefAutBio       = Me%KineticParameters%AmmoniaHalfSatCoefAutBio
        OxyHalfSatCoefAutBio           = Me%KineticParameters%OxyHalfSatCoefAutBio
        DecayCoefAutBio                = Me%KineticParameters%DecayCoefAutBio

        DTDay                          = Me%DTDay

!AmmoniaRateHetBio, Ammonia Rate in function of HetBio, D-1

           AmmoniaRateHetBioAerGrowth =-NCODBioMassRatio*MaxSpecGrowthRateHetBio*                                         &
           (Me%ExternalVar%Mass(ReadilyBioSub, index)/(HalfSatCoefHetBio+Me%ExternalVar%Mass(ReadilyBioSub, index)))      &
           *(Me%ExternalVar%Mass(Oxygen, index)/(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))   
       
        AmmoniaRateHetBioAnoGrowth = -NCODBioMassRatio*MaxSpecGrowthRateHetBio*                                        &
           (Me%ExternalVar%Mass(ReadilyBioSub, index)/(HalfSatCoefHetBio+Me%ExternalVar%Mass(ReadilyBioSub, index)))      &
           *(OxygHalfSatCoefHetBio/(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))               &
           *(Me%ExternalVar%Mass(Nitrate, index)                                                                          &
           /(NitHalfSatCoefDenHetBio+Me%ExternalVar%Mass(Nitrate, index)))*CorFacAnoxGrowthHetBio
           
           AmmoniaRateHetBioAmmonification      = AmmonifRate*Me%ExternalVar%Mass(SolBioOrgNitrogen, index)

           AmmoniaRateHetBio = AmmoniaRateHetBioAerGrowth +AmmoniaRateHetBioAnoGrowth  +AmmoniaRateHetBioAmmonification
            

!AmmoniaRateAutBio, Ammonia Rate in function of AutBio, D-1 
           
           AmmoniaRateAutBioAerGrowth      = -(NCODBioMassRatio+(1./AutoBioYield))*MaxSpecGrowthRateAutBio*               &
           (Me%ExternalVar%Mass(Ammonia, index)/(AmmoniaHalfSatCoefAutBio+Me%ExternalVar%Mass(Ammonia, index)))           &
           *(Me%ExternalVar%Mass(Oxygen, index)/(OxyHalfSatCoefAutBio+Me%ExternalVar%Mass(Oxygen, index)))
           
           
           AmmoniaRateAutBio =  AmmoniaRateAutBioAerGrowth 
         
 
 !Calculation of system coeficients---------------------------------------
         Me%Matrix(Ammonia, HetBio) = DTDay * (-AmmoniaRateHetBio) 
         Me%Matrix(Ammonia, AutBio) = DTDay * (-AmmoniaRateAutBio) 
         Me%Matrix(Ammonia, Ammonia) = 1.0

 !Independent term
        Me%IndTerm(Ammonia) = Me%ExternalVar%Mass(Ammonia, index)

 !------------------------------------------------------------------------
 
    end subroutine WWTPQAmmonia

    subroutine WWTPQSolBioOrgNitrogen (index)
      
    !Subroutine to calculate the propriety SolBioOrgNitrogen

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
    
    !State variables
        integer :: SolInertOrgMat
        integer :: ReadilyBioSub
        integer :: PartInertOrgMar
        integer :: SlowlyBioSub
        integer :: HetBio
        integer :: AutBio
        integer :: PartProd
        integer :: Oxygen
        integer :: Nitrate
        integer :: Ammonia
        integer :: SolBioOrgNitrogen
        integer :: PartBioOrgNitrogen
        integer :: Alkalinity

        real :: HetBioYield
        real :: AutoBioYield
        real :: FracBioPartProd
        real :: NCODBioMassRatio
        real :: NCODPartProdBioMassRatio

        real :: MaxSpecGrowthRateHetBio
        real :: HalfSatCoefHetBio
        real :: OxygHalfSatCoefHetBio
        real :: NitHalfSatCoefDenHetBio
        real :: DecayCoefHetBio
        real :: CorFacAnoxGrowthHetBio
        real :: AmmonifRate
        real :: MaxSpecHydroRate
        real :: HalfSatCoefHydroSlowlyBioSub
        real :: CorFacAnoxCond
        real :: MaxSpecGrowthRateAutBio
        real :: AmmoniaHalfSatCoefAutBio
        real :: OxyHalfSatCoefAutBio
        real :: DecayCoefAutBio
 
        real    :: DTDay
        real    :: SolBioOrgNitrogenRateHetBio                    = null_real
        
        real    :: SolBioOrgNitrogenRateHetBioAmmonification      = null_real
        real    :: SolBioOrgNitrogenRateHetBioHydrolysis          = null_real

!------------------------------------------------------------------------

!State variables
       
        SolInertOrgMat                 = Me%PropIndex%SolInertOrgMat
        ReadilyBioSub                  = Me%PropIndex%ReadilyBioSub
        PartInertOrgMar                = Me%PropIndex%PartInertOrgMar 
        SlowlyBioSub                   = Me%PropIndex%SlowlyBioSub
        HetBio                         = Me%PropIndex%HetBio
        AutBio                         = Me%PropIndex%AutBio
        PartProd                       = Me%PropIndex%PartProd
        Oxygen                         = Me%PropIndex%Oxygen
        Nitrate                        = Me%PropIndex%Nitrate
        Ammonia                        = Me%PropIndex%Ammonia
        SolBioOrgNitrogen              = Me%PropIndex%SolBioOrgNitrogen
        PartBioOrgNitrogen             = Me%PropIndex%PartBioOrgNitrogen 
        Alkalinity                     = Me%PropIndex%Alkalinity 

!Stecheo Parameters
        HetBioYield                    = Me%StechParameters%HetBioYield
        AutoBioYield                   = Me%StechParameters%AutoBioYield
        FracBioPartProd                = Me%StechParameters%FracBioPartProd
        NCODBioMassRatio               = Me%StechParameters%NCODBioMassRatio
        NCODPartProdBioMassRatio       = Me%StechParameters%NCODPartProdBioMassRatio

!Kinetic Parameters
        MaxSpecGrowthRateHetBio        = Me%KineticParameters%MaxSpecGrowthRateHetBio
        HalfSatCoefHetBio              = Me%KineticParameters%HalfSatCoefHetBio
        OxygHalfSatCoefHetBio          = Me%KineticParameters%OxygHalfSatCoefHetBio
        NitHalfSatCoefDenHetBio        = Me%KineticParameters%NitHalfSatCoefDenHetBio
        DecayCoefHetBio                = Me%KineticParameters%DecayCoefHetBio
        CorFacAnoxGrowthHetBio         = Me%KineticParameters%CorFacAnoxGrowthHetBio
        AmmonifRate                    = Me%KineticParameters%AmmonifRate
        MaxSpecHydroRate               = Me%KineticParameters%MaxSpecHydroRate
        HalfSatCoefHydroSlowlyBioSub   = Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub
        CorFacAnoxCond                 = Me%KineticParameters%CorFacAnoxCond
        MaxSpecGrowthRateAutBio        = Me%KineticParameters%MaxSpecGrowthRateAutBio
        AmmoniaHalfSatCoefAutBio       = Me%KineticParameters%AmmoniaHalfSatCoefAutBio
        OxyHalfSatCoefAutBio           = Me%KineticParameters%OxyHalfSatCoefAutBio
        DecayCoefAutBio                = Me%KineticParameters%DecayCoefAutBio

        DTDay                          = Me%DTDay


!SolBioOrgNitrogenRateHetBio, SolBioOrgNitrogen Rate in function of HetBio, D-1

   SolBioOrgNitrogenRateHetBioAmmonification = -AmmonifRate*Me%ExternalVar%Mass(SolBioOrgNitrogen, index)
 
      if (Me%ExternalVar%Mass(HetBio, index) == 0) then
   
   SolBioOrgNitrogenRateHetBioHydrolysis=0
   
   else
     
       SolBioOrgNitrogenRateHetBioHydrolysis= MaxSpecHydroRate*((Me%ExternalVar%Mass(SlowlyBioSub, index)                   &
            /Me%ExternalVar%Mass(HetBio, index))                                                                            &
           /(HalfSatCoefHydroSlowlyBioSub+(Me%ExternalVar%Mass(SlowlyBioSub, index)/Me%ExternalVar%Mass(HetBio, index)))    &
           *(((Me%ExternalVar%Mass(Oxygen, index)/(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))               &
           +CorFacAnoxCond*(OxygHalfSatCoefHetBio/(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index))))              &
           *(Me%ExternalVar%Mass(Nitrate, index)/(NitHalfSatCoefDenHetBio+Me%ExternalVar%Mass(Nitrate, index)))             & 
           *(Me%ExternalVar%Mass(PartBioOrgNitrogen, index)/Me%ExternalVar%Mass(SlowlyBioSub, index))))
endif

           SolBioOrgNitrogenRateHetBio =  SolBioOrgNitrogenRateHetBioAmmonification                                         &
           +SolBioOrgNitrogenRateHetBioHydrolysis                 
        

 !Calculation of system coeficients---------------------------------------
         Me%Matrix(SolBioOrgNitrogen, HetBio) = DTDay * (-SolBioOrgNitrogenRateHetBio) 
         Me%Matrix(SolBioOrgNitrogen, SolBioOrgNitrogen) = 1.0
         
 !Independent term
        Me%IndTerm(SolBioOrgNitrogen) = Me%ExternalVar%Mass(SolBioOrgNitrogen, index)
        
 !------------------------------------------------------------------------

    end subroutine WWTPQSolBioOrgNitrogen

subroutine WWTPQPartBioOrgNitrogen (index)
  
       !Subroutine to calculate the propriety PartBioOrgNitrogen

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
    
    !State variables
        integer :: SolInertOrgMat
        integer :: ReadilyBioSub
        integer :: PartInertOrgMar
        integer :: SlowlyBioSub
        integer :: HetBio
        integer :: AutBio
        integer :: PartProd
        integer :: Oxygen
        integer :: Nitrate
        integer :: Ammonia
        integer :: SolBioOrgNitrogen
        integer :: PartBioOrgNitrogen
        integer :: Alkalinity

        real :: HetBioYield
        real :: AutoBioYield
        real :: FracBioPartProd
        real :: NCODBioMassRatio
        real :: NCODPartProdBioMassRatio

        real :: MaxSpecGrowthRateHetBio
        real :: HalfSatCoefHetBio
        real :: OxygHalfSatCoefHetBio
        real :: NitHalfSatCoefDenHetBio
        real :: DecayCoefHetBio
        real :: CorFacAnoxGrowthHetBio
        real :: AmmonifRate
        real :: MaxSpecHydroRate
        real :: HalfSatCoefHydroSlowlyBioSub
        real :: CorFacAnoxCond
        real :: MaxSpecGrowthRateAutBio
        real :: AmmoniaHalfSatCoefAutBio
        real :: OxyHalfSatCoefAutBio
        real :: DecayCoefAutBio
 
        real :: DTDay
        real :: PartBioOrgNitrogenRateHetBio             = null_real
        real :: PartBioOrgNitrogenRateAutBio             = null_real
        
        real :: PartBioOrgNitrogenRateHetBioDecay        = null_real
        real :: PartBioOrgNitrogenRateHetBioHydrolysis   = null_real
        real :: PartBioOrgNitrogenRateAutBioDecay        = null_real

!------------------------------------------------------------------------

!State variables
       
        SolInertOrgMat                 = Me%PropIndex%SolInertOrgMat
        ReadilyBioSub                  = Me%PropIndex%ReadilyBioSub
        PartInertOrgMar                = Me%PropIndex%PartInertOrgMar 
        SlowlyBioSub                   = Me%PropIndex%SlowlyBioSub
        HetBio                         = Me%PropIndex%HetBio
        AutBio                         = Me%PropIndex%AutBio
        PartProd                       = Me%PropIndex%PartProd
        Oxygen                         = Me%PropIndex%Oxygen
        Nitrate                        = Me%PropIndex%Nitrate
        Ammonia                        = Me%PropIndex%Ammonia
        SolBioOrgNitrogen              = Me%PropIndex%SolBioOrgNitrogen
        PartBioOrgNitrogen             = Me%PropIndex%PartBioOrgNitrogen 
        Alkalinity                     = Me%PropIndex%Alkalinity 

!Stecheo Parameters
        HetBioYield                    = Me%StechParameters%HetBioYield
        AutoBioYield                   = Me%StechParameters%AutoBioYield
        FracBioPartProd                = Me%StechParameters%FracBioPartProd
        NCODBioMassRatio               = Me%StechParameters%NCODBioMassRatio
        NCODPartProdBioMassRatio       = Me%StechParameters%NCODPartProdBioMassRatio

!Kinetic Parameters
        MaxSpecGrowthRateHetBio        = Me%KineticParameters%MaxSpecGrowthRateHetBio
        HalfSatCoefHetBio              = Me%KineticParameters%HalfSatCoefHetBio
        OxygHalfSatCoefHetBio          = Me%KineticParameters%OxygHalfSatCoefHetBio
        NitHalfSatCoefDenHetBio        = Me%KineticParameters%NitHalfSatCoefDenHetBio
        DecayCoefHetBio                = Me%KineticParameters%DecayCoefHetBio
        CorFacAnoxGrowthHetBio         = Me%KineticParameters%CorFacAnoxGrowthHetBio
        AmmonifRate                    = Me%KineticParameters%AmmonifRate
        MaxSpecHydroRate               = Me%KineticParameters%MaxSpecHydroRate
        HalfSatCoefHydroSlowlyBioSub   = Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub
        CorFacAnoxCond                 = Me%KineticParameters%CorFacAnoxCond
        MaxSpecGrowthRateAutBio        = Me%KineticParameters%MaxSpecGrowthRateAutBio
        AmmoniaHalfSatCoefAutBio       = Me%KineticParameters%AmmoniaHalfSatCoefAutBio
        OxyHalfSatCoefAutBio           = Me%KineticParameters%OxyHalfSatCoefAutBio
        DecayCoefAutBio                = Me%KineticParameters%DecayCoefAutBio

        DTDay                          = Me%DTDay


!PartBioOrgNitrogenRateHetBio, PartBioOrgNitrogen Rate in function of HetBio, D-1

          PartBioOrgNitrogenRateHetBioDecay = (NCODBioMassRatio-FracBioPartProd*NCODPartProdBioMassRatio)*DecayCoefHetBio 
 
       if (Me%ExternalVar%Mass(HetBio, index) == 0) then
 
            PartBioOrgNitrogenRateHetBioHydrolysis=0
 
      else
 
           PartBioOrgNitrogenRateHetBioHydrolysis = - MaxSpecHydroRate*((Me%ExternalVar%Mass(SlowlyBioSub, index)              &
       /Me%ExternalVar%Mass(HetBio, index))                                                                                    &
           /(HalfSatCoefHydroSlowlyBioSub+(Me%ExternalVar%Mass(SlowlyBioSub, index)/Me%ExternalVar%Mass(HetBio, index)))       &
           *(((Me%ExternalVar%Mass(Oxygen, index)/(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))                  &
           +CorFacAnoxCond*(OxygHalfSatCoefHetBio/(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index))))                 &
           *(Me%ExternalVar%Mass(Nitrate, index)/(NitHalfSatCoefDenHetBio+Me%ExternalVar%Mass(Nitrate, index)))                & 
           *(Me%ExternalVar%Mass(PartBioOrgNitrogen, index)/Me%ExternalVar%Mass(SlowlyBioSub, index))))
 
     endif
        
           PartBioOrgNitrogenRateHetBio = PartBioOrgNitrogenRateHetBioDecay + PartBioOrgNitrogenRateHetBioHydrolysis

!PartBioOrgNitrogenRateAutBio, PartBioOrgNitrogen Rate in relation of AutBio, units   
           
           PartBioOrgNitrogenRateAutBioDecay = (NCODBioMassRatio-FracBioPartProd*NCODPartProdBioMassRatio)*DecayCoefAutBio 
           
           PartBioOrgNitrogenRateAutBio = PartBioOrgNitrogenRateAutBioDecay
           
 
 !Calculation of system coeficients---------------------------------------
         Me%Matrix(PartBioOrgNitrogen, HetBio) = DTDay * (-PartBioOrgNitrogenRateHetBio) 
         Me%Matrix(PartBioOrgNitrogen, AutBio) = DTDay * (-PartBioOrgNitrogenRateAutBio) 
         Me%Matrix(PartBioOrgNitrogen, PartBioOrgNitrogen) = 1.0

 !Independent term
        Me%IndTerm(PartBioOrgNitrogen) = Me%ExternalVar%Mass(PartBioOrgNitrogen, index)

 !------------------------------------------------------------------------
 
    end subroutine WWTPQPartBioOrgNitrogen

    subroutine WWTPQAlkalinity  (index)
  
       !Subroutine to calculate Alkalinity

    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
    
    !State variables
        integer :: SolInertOrgMat
        integer :: ReadilyBioSub
        integer :: PartInertOrgMar
        integer :: SlowlyBioSub
        integer :: HetBio
        integer :: AutBio
        integer :: PartProd
        integer :: Oxygen
        integer :: Nitrate
        integer :: Ammonia
        integer :: SolBioOrgNitrogen
        integer :: PartBioOrgNitrogen
        integer :: Alkalinity

        real :: HetBioYield
        real :: AutoBioYield
        real :: FracBioPartProd
        real :: NCODBioMassRatio
        real :: NCODPartProdBioMassRatio

        real :: MaxSpecGrowthRateHetBio
        real :: HalfSatCoefHetBio
        real :: OxygHalfSatCoefHetBio
        real :: NitHalfSatCoefDenHetBio
        real :: DecayCoefHetBio
        real :: CorFacAnoxGrowthHetBio
        real :: AmmonifRate
        real :: MaxSpecHydroRate
        real :: HalfSatCoefHydroSlowlyBioSub
        real :: CorFacAnoxCond
        real :: MaxSpecGrowthRateAutBio
        real :: AmmoniaHalfSatCoefAutBio
        real :: OxyHalfSatCoefAutBio
        real :: DecayCoefAutBio
 
        real :: DTDay
        real :: AlkalinityRateHetBio                 = null_real
        real :: AlkalinityRateAutBio                 = null_real
        
        real :: AlkalinityRateHetBioAerGrowth        = null_real
        real :: AlkalinityRateHetBioAnoGrowth        = null_real
        real :: AlkalinityRateHetBioAmmonification   = null_real
        real :: AlkalinityRateAutBioAerGrowth        = null_real

!------------------------------------------------------------------------

!State variables
       
        SolInertOrgMat                 = Me%PropIndex%SolInertOrgMat
        ReadilyBioSub                  = Me%PropIndex%ReadilyBioSub
        PartInertOrgMar                = Me%PropIndex%PartInertOrgMar 
        SlowlyBioSub                   = Me%PropIndex%SlowlyBioSub
        HetBio                         = Me%PropIndex%HetBio
        AutBio                         = Me%PropIndex%AutBio
        PartProd                       = Me%PropIndex%PartProd
        Oxygen                         = Me%PropIndex%Oxygen
        Nitrate                        = Me%PropIndex%Nitrate
        Ammonia                        = Me%PropIndex%Ammonia
        SolBioOrgNitrogen              = Me%PropIndex%SolBioOrgNitrogen
        PartBioOrgNitrogen             = Me%PropIndex%PartBioOrgNitrogen 
        Alkalinity                     = Me%PropIndex%Alkalinity 

!Stecheo Parameters
        HetBioYield                    = Me%StechParameters%HetBioYield
        AutoBioYield                   = Me%StechParameters%AutoBioYield
        FracBioPartProd                = Me%StechParameters%FracBioPartProd
        NCODBioMassRatio               = Me%StechParameters%NCODBioMassRatio
        NCODPartProdBioMassRatio       = Me%StechParameters%NCODPartProdBioMassRatio

!Kinetic Parameters
        MaxSpecGrowthRateHetBio        = Me%KineticParameters%MaxSpecGrowthRateHetBio
        HalfSatCoefHetBio              = Me%KineticParameters%HalfSatCoefHetBio
        OxygHalfSatCoefHetBio          = Me%KineticParameters%OxygHalfSatCoefHetBio
        NitHalfSatCoefDenHetBio        = Me%KineticParameters%NitHalfSatCoefDenHetBio
        DecayCoefHetBio                = Me%KineticParameters%DecayCoefHetBio
        CorFacAnoxGrowthHetBio         = Me%KineticParameters%CorFacAnoxGrowthHetBio
        AmmonifRate                    = Me%KineticParameters%AmmonifRate
        MaxSpecHydroRate               = Me%KineticParameters%MaxSpecHydroRate
        HalfSatCoefHydroSlowlyBioSub   = Me%KineticParameters%HalfSatCoefHydroSlowlyBioSub
        CorFacAnoxCond                 = Me%KineticParameters%CorFacAnoxCond
        MaxSpecGrowthRateAutBio        = Me%KineticParameters%MaxSpecGrowthRateAutBio
        AmmoniaHalfSatCoefAutBio       = Me%KineticParameters%AmmoniaHalfSatCoefAutBio
        OxyHalfSatCoefAutBio           = Me%KineticParameters%OxyHalfSatCoefAutBio
        DecayCoefAutBio                = Me%KineticParameters%DecayCoefAutBio

        DTDay                          = Me%DTDay


!AlkalinityRateHetBio, Alkalinity Rate in function of HetBio, D-1

         AlkalinityRateHetBioAerGrowth      = (-NCODBioMassRatio/14.)*MaxSpecGrowthRateHetBio                             &
           *(Me%ExternalVar%Mass(ReadilyBioSub, index)/(HalfSatCoefHetBio+Me%ExternalVar%Mass(ReadilyBioSub, index)))     &
           *(Me%ExternalVar%Mass(Oxygen, index)                                                                           &
           /(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))
            
         AlkalinityRateHetBioAnoGrowth= ((1.-HetBioYield)/(14.*2.86*HetBioYield)-NCODBioMassRatio/14.)                    &
           *MaxSpecGrowthRateHetBio*                                                                                      &
          (Me%ExternalVar%Mass(ReadilyBioSub, index)                                                                      &
           /(HalfSatCoefHetBio+Me%ExternalVar%Mass(ReadilyBioSub, index)))*(OxygHalfSatCoefHetBio                         &
           /(OxygHalfSatCoefHetBio+Me%ExternalVar%Mass(Oxygen, index)))                                                   &
           *(Me%ExternalVar%Mass(Nitrate, index)/(NitHalfSatCoefDenHetBio+Me%ExternalVar%Mass(Nitrate, index)))           &
           *CorFacAnoxGrowthHetBio  
        
         AlkalinityRateHetBioAmmonification = (1./14.)*AmmonifRate*Me%ExternalVar%Mass(SolBioOrgNitrogen , index)
        
           AlkalinityRateHetBio = AlkalinityRateHetBioAerGrowth +  AlkalinityRateHetBioAnoGrowth                          &
           +   AlkalinityRateHetBioAmmonification                             
             
!AlkalinityRateAutBio, Alkalinity Rate in function of AutBio, D-1 
                    
           AlkalinityRateAutBioAerGrowth      = (-NCODBioMassRatio/14.-1./(7.*AutoBioYield))*MaxSpecGrowthRateAutBio      &
           *(Me%ExternalVar%Mass(Ammonia , index)/(AmmoniaHalfSatCoefAutBio+Me%ExternalVar%Mass(Ammonia , index)))        &
            *(Me%ExternalVar%Mass(Oxygen , index)/(OxyHalfSatCoefAutBio+Me%ExternalVar%Mass(Oxygen , index)))
  
           AlkalinityRateAutBio = AlkalinityRateAutBioAerGrowth  
           
 
 !Calculation of system coeficients---------------------------------------
         Me%Matrix(Alkalinity, HetBio) = DTDay * (-AlkalinityRateHetBio) 
         Me%Matrix(Alkalinity, AutBio) = DTDay * (-AlkalinityRateAutBio) 
         Me%Matrix(Alkalinity, Alkalinity) = 1.0

 !Independent term
        Me%IndTerm(Alkalinity) = Me%ExternalVar%Mass(Alkalinity, index)

 !------------------------------------------------------------------------

    end subroutine WWTPQAlkalinity 


!    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR
!
!    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
!
    subroutine KillWWTPQ(WWTPQID, STAT)
    
    !?What does it do?

        !Arguments---------------------------------------------------------------
        integer                        :: WWTPQID
        integer, optional, intent(OUT) :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_CALL
        integer                                     :: STAT_     
        integer                                     :: nUsers
        !------------------------------------------------------------------------                      

        STAT_ = UNKNOWN_

        call Ready(WWTPQID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mWWTPQ_,  Me%InstanceID)
  
            if (nUsers == 0) then

                if(Me%ObjLUD /= 0) then
                    call KillLUD(Me%ObjLUD, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Subroutine KillWWTPQ; module ModuleWWTPQ. ERR01.'
                end if

                deallocate(Me%IndTerm, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'Subroutine Kill_WWTPQ; module ModuleWWTPQ. ERR02.'
                nullify(Me%IndTerm)

                deallocate(Me%Matrix, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'Subroutine Kill_WWTPQ; module ModuleWWTPQ. ERR03.'
                nullify(Me%Matrix)

cd4 :           if (associated(Me%NewMass)) then
                    deallocate(Me%NewMass, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Subroutine Kill_WWTPQ; module ModuleWWTPQ. ERR04.'
                    nullify(Me%NewMass)
                end if cd4

                call DeallocateInstance

                WWTPQID = 0

                STAT_ = SUCCESS_

            endif

        else 

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    !------------------------------------------------------------------------
!
    end subroutine KillWWTPQ
!
!    !------------------------------------------------------------------------
!    
!    
!    !------------------------------------------------------------------------
!
    subroutine DeallocateInstance

        !Local-----------------------------------------------------------------
        type (T_WWTPQ), pointer           :: AuxObjWWTPQ
        type (T_WWTPQ), pointer           :: PreviousObjWWTPQ

        !Updates pointers
        if (Me%InstanceID == FirstObjWWTPQ%InstanceID) then
            FirstObjWWTPQ => FirstObjWWTPQ%Next
        else
            PreviousObjWWTPQ => FirstObjWWTPQ
            AuxObjWWTPQ      => FirstObjWWTPQ%Next
            do while (AuxObjWWTPQ%InstanceID /= Me%InstanceID)
                PreviousObjWWTPQ => AuxObjWWTPQ
                AuxObjWWTPQ      => AuxObjWWTPQ%Next
            enddo

            !Now update linked list
            PreviousObjWWTPQ%Next => AuxObjWWTPQ%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 
            
    end subroutine DeallocateInstance

!    !--------------------------------------------------------------------------

!    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME
!
!    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
    subroutine Ready (WWTPQID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: WWTPQID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (me)

cd1:    if (wwtpqid > 0) then
            call locateobjwwtpq (wwtpqid)
            ready_ = verifyreadlock (mwwtpq_, me%instanceid)
        else
            ready_ = off_err_
        end if cd1
     
        !----------------------------------------------------------------------

    end subroutine ready

    !--------------------------------------------------------------------------

    subroutine locateobjwwtpq (wwtpqid)

        !arguments-------------------------------------------------------------
        integer                                     :: wwtpqid

        !local-----------------------------------------------------------------

        me => firstobjwwtpq
        do while (associated (me))
            if (me%instanceid == wwtpqid) exit
            me => me%next
        enddo

        if (.not. associated(me)) stop 'modulewwtpq - locateobjwwtpq - err01'

    end subroutine locateobjwwtpq

    !--------------------------------------------------------------------------

end module ModuleWWTPQ

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
