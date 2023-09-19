!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid HNS
! PROJECT       : Mohid Water
! MODULE        : ModuleHNS
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : October 2013
! REVISION      : Rodrigo Fernandes
! DESCRIPTION   : Module responsbile for computing physical processes of HNS spills
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
 
Module ModuleHNS
!BOP
!
! !MODULE: ModuleHNS

!    !DESCRIPTION: 
!     This model is responsible  for determine specific behaviour of chemical lagrangian Tracers
 

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleOil_0D
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartHNS
    public  :: ReadFinalHNS
    private :: FindHNSBehaviourClass

    !Selector
    public :: GetHNSAtmosphericDispersion
    public :: GetHNSSolarRadiationNeeded
    public :: GetHNSSedimentation
    public :: GetHNSEntrainment
    public :: GetHNSSpreading
    public :: GetHNSInitialMass
    public :: GetHNSDensity
    public :: GetInitialArea
    public :: GetHNSSpreadingDiffVel
    public :: GetHNSBehaviourClass
    public  :: GetHNSInitialState
    public :: ExportDropletOptions
    public :: GetHNSAirTurbHorizontalVel
    PUBLIC :: GetHNSAirTurbVerticalVel
    public :: GetHNSPAC
    public :: GetHNSPNEC
    public :: GetHNSLC50

    !Modifier
    public  ::  HNSInternalProcesses    
    private :: FindInitialDropletDiameter

    !Destructor
    public  ::  KillHNS 
    public  ::  WriteFinalHNS

    !Management
    private ::  Ready
    private ::  LocateObjHNS


    ! Parameters
    real, parameter     :: WaterMolecularWeight     = 18.015 !kg/kmol
    real, parameter     :: WaterDiffusionCoef       = 2.3e-9 !m2/s
    real, parameter     :: AirCinematicVisc         = 1.5e-5 !m2/s
    
    !based on observations by Audunson et al. (1984), Sørstrøm and Johansen (1985), 
    !JBF (1976), Sørstrøm (1989), and Reed et al (1990, 1992) - original value 5 x 108 day-1
    real, parameter     :: SpreadingRateConstant    = 5787.037 !s-1

    integer, parameter :: Gas_                        = 1
    integer, parameter :: GasDissolver_               = 2 
    integer, parameter :: Evaporator_                 = 3
    integer, parameter :: EvaporatorDissolver_        = 4
    integer, parameter :: FloaterEvaporator_          = 5
    integer, parameter :: FloaterEvaporatorDissolver_ = 6
    integer, parameter :: Floater_                    = 7
    integer, parameter :: FloaterDissolver_           = 8
    integer, parameter :: DissolverEvaporator_        = 9
    integer, parameter :: Dissolver_                  = 10
    integer, parameter :: SinkerDissolver_            = 11
    integer, parameter :: Sinker_                     = 12
    integer, parameter :: ClassUnknown_               = 13

    integer, parameter :: Class_A_                    = 1
    integer, parameter :: Class_B_                    = 2
    integer, parameter :: Class_C_                    = 3
    integer, parameter :: Class_D_                    = 4
    integer, parameter :: Class_E_                    = 5
    integer, parameter :: Class_F_                    = 6
    
    character(LEN = StringLength), parameter    :: LC50_begin         = '<<<BeginLC50>>>'
    character(LEN = StringLength), parameter    :: LC50_end           = '<<<EndLC50>>>'

    integer, parameter :: LC50Surface_                    = 1
    integer, parameter :: LC50WaterColumn_                = 2
    integer, parameter :: LC50Bottom_                     = 3
    
    !Types---------------------------------------------------------------------

    type       T_State
        logical                                              :: FirstStepIP = ON
    end type T_State

    type       T_Files
        character(LEN = StringLength)                        :: ConstructData  = null_str
    end type T_Files

    type       T_Var
        !Time
        real                                                 :: DTHNSInternalProcesses  = null_real
        real                                                 :: Time                 = null_real
        real                                                 :: Mass                 = null_real
        real                                                 :: Density              = null_real
        real                                                 :: Viscosity            = null_real
        real                                                 :: ViscCin              = null_real
        real                                                 :: MolecularWeight      = null_real
        real                                                 :: VaporPressure        = null_real
        real                                                 :: WaterSolubility      = null_real
        real                                                 :: OWPartitionCoef      = null_real
        real                                                 :: AirDegradationRate   = null_real
        real                                                 :: WaterDegradationRate = null_real
        real                                                 :: SedimentDegradationRate  = null_real
        integer                                              :: HNSBehaviourClass    = null_int
        real                                                 :: InitialMass          = null_int
        integer                                              :: HNSParticleState     = null_int
        
        !Force chemical to behave like a gas (exclusively atmospheric processes)
        logical                                              :: HNSGas              = OFF
        
        !compute PAC (exclusively atmospheric processes)
        logical                                              :: PAC                = OFF
        integer                                              :: PAC1               = null_int
        integer                                              :: PAC2               = null_int
        integer                                              :: PAC3               = null_int

        !compute PNEC 
        real                                                 :: PNEC               = null_real
        
        !compute LC50 
        logical                                              :: LC50               = OFF
!        real                                                 :: LC50Hours          = null_real
        integer                                              :: LC50Number           = null_int
        character(LEN = StringLength), dimension(:), pointer :: listLC50Specie      => null()        
        real, dimension(:), pointer                          :: listLC50Value       => null()        
        real, dimension(:), pointer                          :: listLC50Hours       => null()        
        integer, dimension(:), pointer                       :: listLC50Type        => null()        

        !Spreading
        logical                                              :: HNSSpreading         = OFF
        
        !Evaporation
        logical                                              :: HNSEvaporation       = OFF
        real                                                 :: MEvaporatedDT        = null_real
        real                                                 :: MEvaporated          = null_real
        logical                                              :: HNSAtmosphericDispersion = OFF
        real                                                 :: AirHorizontalDifCoef = null_real
        real                                                 :: AirVerticalDifCoef   = null_real
        integer                                              :: StabilityClass       = null_int
        logical                                              :: SolarRadiationNeeded = OFF
        
        !Volatilization
        logical                                              :: HNSVolatilization    = OFF
        real                                                 :: MVolatilizedDT       = null_real
        real                                                 :: MVolatilized         = null_real

        !Entrainment
        logical                                              :: HNSEntrainment       = OFF
        real                                                 :: DropletDiameter      = null_real
        real                                                 :: MEntrainedDT         = null_real
        real                                                 :: MEntrained           = null_real
       
        !Dissolution
        logical                                              :: HNSDissolution       = OFF
        real                                                 :: MDissolvedDT         = null_real
        real                                                 :: MDissolved           = null_real
        logical                                              :: Organic              = OFF
        
        !Sedimentation
        logical                                              :: HNSSedimentation     = OFF
        real                                                 :: MSedimentedDT        = null_real
        real                                                 :: MSedimented          = null_real
        
        !Degradation
        logical                                              :: HNSDegradation       = OFF
        real                                                 :: MDegradedDT          = null_real
        real                                                 :: MDegraded            = null_real
        
    end type   T_Var
    
    type       T_External
        type(T_Time)                                         :: Now
        type(T_Time)                                         :: BeginTime
        type(T_Time)                                         :: EndTime

        
        real                                                 :: Wind                 = null_real
        real                                                 :: WaterTemperature     = null_real
        real                                                 :: WaterDensity         = null_real
        real                                                 :: SPM                  = null_real
        real                                                 :: Currents             = null_real
        real                                                 :: AirTemperature       = null_real
        real                                                 :: AtmPressure          = null_real
        real                                                 :: WaveHeight           = null_real
        real                                                 :: WavePeriod           = null_real
        real                                                 :: WaterSolubility      = null_real
        real                                                 :: ParticleState        = null_real
        real                                                 :: Depth                = null_real
        real                                                 :: Area                 = null_real
        logical                                              :: BackTracking         = OFF
        real                                                 :: DropletsD50          = null_real
        integer                                              :: MethodBWDropletsDiameter = null_int
    end type   T_External
    
    type      T_HNS
        integer                                              :: InstanceID = 0
        type(T_State    )                                    :: State
        type(T_Files    )                                    :: Files
        type(T_Var      )                                    :: Var
        type(T_External )                                    :: ExternalVar

        type(T_Time)                                         :: NextInternalComputeTime
        type(T_Time)                                         :: NextActiveComputeTime
        type(T_Time)                                         :: Now
        type(T_Time)                                         :: BeginTime

        !Instance of ModuleTime
        integer                                              :: ObjTime = 0

        !Instance of Module_EnterData
        integer                                              :: ObjEnterData = 0

        type (T_HNS), pointer                                :: Next         => null()

    end type T_HNS

    !Global Module Variables
    type (T_HNS), pointer                                    :: FirstObjHNS  => null()
    type (T_HNS), pointer                                    :: Me           => null()



    !----------------------------------------------------------------------------

    contains



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine StartHNS (HNSID,                                                 &
                         TimeID,                                                &
                         EnterDataID,                                           &
                         DT,                                                    &
                         ContCalc,                                              &
                         ClientNumber,                                          &
                         ExtractType,                                           &
                         STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: HNSID
        integer                                     :: TimeID
        integer                                     :: EnterDataID
        integer,           intent(IN )              :: ClientNumber    
        integer, optional, intent(IN )              :: ExtractType    
        integer, optional, intent(OUT)              :: STAT 
        real,              intent(IN )              :: DT
        
        logical, intent(IN)  :: ContCalc 

        !External----------------------------------------------------------------

        integer :: STAT_CALL
        integer :: ready_ 
        integer :: FromFile
        integer :: ExtractType_  
        integer :: nUsers
        !Local-------------------------------------------------------------------

        integer :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mHNS_)) then
            nullify (FirstObjHNS)
            call RegisterModule (mHNS_) 
        endif
        
        call Ready(HNSID, ready_) 

        if (ready_ .EQ. OFF_ERR_) then

            !Allocates Instance
            call AllocateInstance


            !Associates External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjEnterData      = AssociateInstance (mENTERDATA_,      EnterDataID     )

            !Gets Time
            call GetComputeTimeLimits(Me%ObjTime,                                       &
                                      BeginTime = Me%ExternalVar%BeginTime,             &
                                      EndTime   = Me%ExternalVar%EndTime,               &
                                      STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartHNS - ModuleHNS - ERR10'

            ! Actualized the time
            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'StartHNS - ModuleHNS - ERR20'

            ! Check if the simulation goes backward in time or forward in time (default mode)
            call GetBackTracking(Me%ObjTime,                                            &
                                 Me%ExternalVar%BackTracking, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'StartHNS - ModuleHNS - ERR30'

cd2 :       if (present(ExtractType)) then
                ExtractType_ = ExtractType
            else
                call GetExtractType(FromFile = FromFile)
                ExtractType_ = FromFile
            end if cd2

            call HNSOptions       (DT, ContCalc, ClientNumber, ExtractType_)

            !Initialization of integrated values
            
ifContCalc: if (.NOT. ContCalc ) then
!                Me%Var%Time              = 0.0
                Me%Var%MEvaporated              = 0.0
                Me%Var%MEntrained               = 0.0
                Me%Var%MSedimented              = 0.0
                Me%Var%MDissolved               = 0.0
                Me%Var%MVolatilized             = 0.0
                Me%Var%MDegraded                = 0.0

!                Me%Var%Volume            = 0.0
                Me%State%FirstStepIP     = ON

            end if ifContCalc
            

            Me%NextInternalComputeTime = Me%ExternalVar%Now + Me%Var%DTHNSInternalProcesses
            Me%NextActiveComputeTime   = Me%ExternalVar%Now + DT
                  
            nUsers = DeassociateInstance (mENTERDATA_,      Me%ObjEnterData     )
            if (nUsers == 0) stop 'StartHNS - ModuleHNS - ERR98'

            HNSID = Me%InstanceID

            STAT_ = SUCCESS_
        
        else 
            
            stop 'StartOil - ModuleHNS - ERR99' 

        end if 

        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine StartHNS

    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance
                                                    
        !Local-----------------------------------------------------------------
        type (T_HNS), pointer                  :: NewObjHNS
        type (T_HNS), pointer                  :: PreviousObjHNS


        !Allocates new instance
        allocate (NewObjHNS)
        nullify  (NewObjHNS%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjHNS)) then
            FirstObjHNS         => NewObjHNS
            Me                  => NewObjHNS
        else
            PreviousObjHNS      => FirstObjHNS
            Me                  => FirstObjHNS%Next
            do while (associated(Me))
                PreviousObjHNS  => Me
                Me              => Me%Next
            enddo
            Me                  => NewObjHNS
            PreviousObjHNS%Next => NewObjHNS
        endif

        Me%InstanceID = RegisterNewInstance (mHNS_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine HNSOptions(DT, ContCalc, ClientNumber, ExtractType) 

        !Arguments---------------------------------------------------------------
        real   , intent(IN)     :: DT
        logical, intent(IN)     :: ContCalc
        integer, intent(IN)     :: ClientNumber

        !External----------------------------------------------------------------
        integer                 :: flag
        integer                 :: ExtractType
        integer                 :: STAT_CALL


        !Internal------------------------------------------------------------------------

!        character(LEN = StringLength), parameter :: Char_Amount                = trim(adjustl('Amount'))

        !--------------------------------------------------------------------------------


        Me%NextInternalComputeTime = Me%ExternalVar%Now
        Me%NextActiveComputeTime   = Me%ExternalVar%Now
        Me%BeginTime               = Me%ExternalVar%Now

notcc : if (.NOT. ContCalc ) then
            call GetData(Me%Var%DTHNSInternalProcesses,                                 &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType = ExtractType,                                      &    
                         keyword    = 'DT_HNS_INTPROCESSES',                            &         
                         ClientModule ='ModuleHNS',                                  &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OilOptions - ModuleHNS - ERR10'


cd1 :       if (flag .EQ. 0) then

                Me%Var%DTHNSInternalProcesses = DT

            end if cd1                  

            call GetData(Me%Var%HNSGas,                                                 &
                            Me%ObjEnterData,                                            &
                            flag,                                                       &
                            SearchType     = ExtractType,                               & 
                            keyword        = 'HNS_GAS',                                 &
                            Default        = OFF,                                       &
                            ClientModule   ='ModuleHNS',                                &
                            STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR21'

            call GetData(Me%Var%HNSSpreading,                                           &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'HNS_SPREADING',                              &
                         Default        = OFF,                                          &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR15'

            call GetData(Me%Var%HNSEvaporation,                                         &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'HNS_EVAPORATION',                            &
                         Default        = OFF,                                          &
                         ClientModule   ='ModuleHNS',                                &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR20'

            call GetData(Me%Var%HNSEntrainment,                                         &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'HNS_ENTRAINMENT',                            &
                         Default        = OFF,                                          &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'

            call GetData(Me%Var%HNSDissolution,                                          &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'HNS_DISSOLUTION',                            &
                         Default        = OFF,                                          &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'

            if (Me%Var%HNSDissolution) then
                call GetData(Me%Var%HNSVolatilization,                                  &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType     = ExtractType,                              & 
                             keyword        = 'HNS_VOLATILIZATION',                     &
                             Default        = OFF,                                      &
                             ClientModule   ='ModuleHNS',                               &
                             STAT           = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR25'           
                call GetData(Me%Var%Organic,                                          &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             SearchType     = ExtractType,                                  & 
                             keyword        = 'ORGANIC_COMPOUND',                           &
                             Default        = OFF,                                          &
                             ClientModule   ='ModuleHNS',                                   &
                             STAT           = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'
            endif

            if (Me%Var%HNSEvaporation .OR. Me%Var%HNSVolatilization) then
                call GetData(Me%Var%HNSAtmosphericDispersion,                               &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             SearchType     = ExtractType,                                  & 
                             keyword        = 'HNS_ATMOSPHERIC_DISPERSION',                 &
                             Default        = OFF,                                          &
                             ClientModule   ='ModuleHNS',                                   &
                             STAT           = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR21'
            endif

            if (Me%Var%HNSAtmosphericDispersion .OR. Me%Var%HNSGas) then
                call GetData(Me%Var%PAC,                                                   &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             SearchType   = FromBlock,                                      &
                             keyword      ='PAC',                                          &
                             ClientModule ='ModuleHNS',                                     &
                             Default      = .FALSE.,                                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR22'
                if (Me%Var%PAC) then
                    call GetData(Me%Var%PAC1,                                                  &
                                 Me%ObjEnterData,                                               &
                                 flag,                                                          &
                                 SearchType   = FromBlock,                                      &
                                 keyword      ='PAC_1',                                        &
                                 ClientModule ='ModuleHNS',                                     &
                                 Default      = FillValueInt,                                  &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR23'
                    call GetData(Me%Var%PAC2,                                                  &
                                 Me%ObjEnterData,                                               &
                                 flag,                                                          &
                                 SearchType   = FromBlock,                                      &
                                 keyword      ='PAC_2',                                        &
                                 ClientModule ='ModuleHNS',                                     &
                                 Default      = FillValueInt,                                  &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR24'
                    call GetData(Me%Var%PAC3,                                                  &
                                 Me%ObjEnterData,                                               &
                                 flag,                                                          &
                                 SearchType   = FromBlock,                                      &
                                 keyword      ='PAC_3',                                        &
                                 ClientModule ='ModuleHNS',                                     &
                                 Default      = FillValueInt,                                  &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR25'
                endif
            endif

            call GetData(Me%Var%PNEC,                                                   &
                            Me%ObjEnterData,                                            &
                            flag,                                                          &
                            SearchType   = FromBlock,                                      &
                            keyword      ='PNEC',                                          &
                            ClientModule ='ModuleHNS',                                     &
                            Default      = FillValueReal,                                  &
                            STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR25a'

            !call GetData(Me%Var%LC50,                                                   &
            !                Me%ObjEnterData,                                            &
            !                flag,                                                          &
            !                SearchType   = FromBlock,                                      &
            !                keyword      ='LC50',                                          &
            !                ClientModule ='ModuleHNS',                                     &
            !                Default      = FillValueReal,                                  &
            !                STAT         = STAT_CALL)
            !if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR25b'
            !
            !if (Me%Var%LC50 .GT. 0) then
            !    call GetData(Me%Var%LC50Hours,                                                 &
            !                    Me%ObjEnterData,                                               &
            !                    flag,                                                          &
            !                    SearchType   = FromBlock,                                      &
            !                    keyword      ='LC50_HOURS',                                    &
            !                    ClientModule ='ModuleHNS',                                     &
            !                    Default      = 96.,                                             &
            !                    STAT         = STAT_CALL)
            !    if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR25c'
            !endif
            
            call GetData(Me%Var%HNSSedimentation,                                       &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'HNS_SEDIMENTATION',                          &
                         Default        = OFF,                                          &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'

            call GetData(Me%Var%HNSDegradation,                                         &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'HNS_DEGRADATION',                            &
                         Default        = OFF,                                          &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'
            
            ! kg/m3
            call GetData(Me%Var%Density,                                                &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'DENSITY',                                    &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'

ifentr: if  (Me%Var%HNSEntrainment) then
            ! cP at 20ºC
            call GetData(Me%Var%Viscosity,                                              &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'VISCOSITY',                                  &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'

           !units - cSt
            Me%Var%ViscCin    = 1.0e6 * (Me%Var%Viscosity/1000.) / Me%Var%Density 

        end if ifentr
        
ifevap: if  (Me%Var%HNSEvaporation) then
            ! kg/kmol
            call GetData(Me%Var%MolecularWeight,                                        &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'MOLECULARWEIGHT',                            &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'


        end if ifevap

        ! Pa
        call GetData(Me%Var%VaporPressure,                                          &
                     Me%ObjEnterData,                                               &
                     flag,                                                          &
                     SearchType     = ExtractType,                                  & 
                     keyword        = 'VAPORPRESSURE',                              &
                     ClientModule   ='ModuleHNS',                                   &
                     STAT           = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'

!ifdiss: if  (Me%Var%HNSDissolution) then
            ! mg/L
            call GetData(Me%Var%WaterSolubility,                                      &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'WATERSOLUBILITY',                            &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'
!        end if ifdiss
        
ifsed:  if  (Me%Var%HNSSedimentation) then
            call GetData(Me%Var%OWPartitionCoef,                                        &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'OWPARTITIONCOEF',                            &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'
        end if ifsed

ifdegr:  if  (Me%Var%HNSDegradation) then
            call GetData(Me%Var%AirDegradationRate,                                     &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'AIR_DEGRADATIONRATE',                        &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'
            call GetData(Me%Var%WaterDegradationRate,                                   &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'WATER_DEGRADATIONRATE',                      &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'
            call GetData(Me%Var%SedimentDegradationRate,                                &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'SEDIMENT_DEGRADATIONRATE',                   &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'
            if ((Me%Var%AirDegradationRate .LT. 0) .OR.                                 &
                (Me%Var%WaterDegradationRate .LT. 0) .OR.                               &
                (Me%Var%SedimentDegradationRate .LT. 0)) then
                if (Me%Var%AirDegradationRate .LT. 0) then
                    write(*,*) 'AIR_DEGRADATIONRATE has a negative value. A default value of 0 will be assumed.'
                    Me%Var%AirDegradationRate = 0.
                endif
                if (Me%Var%WaterDegradationRate .LT. 0) then
                    write(*,*) 'WATER_DEGRADATIONRATE has a negative value. A defaultvalue of 0 will be assumed.'
                    Me%Var%WaterDegradationRate = 0.
                endif
                if (Me%Var%SedimentDegradationRate .LT. 0) then
                    write(*,*) 'SEDIMENT_DEGRADATIONRATE has a negative value. A defaultvalue of 0 will be assumed.'
                    Me%Var%SedimentDegradationRate = 0.
                endif
                write(*,*) 'HNSOptions - ModuleHNS - WRN009'             
            endif
            
        end if ifdegr

        if (Me%Var%HNSGas) then
            Me%Var%HNSAtmosphericDispersion = .true.
            Me%Var%HNSBehaviourClass        = Gas_
            Me%Var%HNSParticleState         = Air_Evaporated_
            Me%Var%PNEC                     = FillValueReal
        endif

        if (Me%Var%HNSAtmosphericDispersion) then
            call GetData(Me%Var%AirHorizontalDifCoef,                                   &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'HORIZONTALDIFCOEF',                          &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'

            call GetData(Me%Var%AirVerticalDifCoef,                                     &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'VERTICALDIFCOEF',                            &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'

            call GetData(Me%Var%StabilityClass,                                         &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType     = ExtractType,                                  & 
                         keyword        = 'STABILITY_CLASS',                            &
                         ClientModule   ='ModuleHNS',                                   &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'
            if (Me%Var%StabilityClass .LT. 0.) then
                Me%Var%SolarRadiationNeeded = .true.     
            endif
            
        endif

            if (Me%ExternalVar%BackTracking) then
                Me%Var%HNSSpreading             = .false.
                Me%Var%HNSEvaporation           = .false.
                Me%Var%HNSEntrainment           = .false.
                Me%Var%HNSDissolution           = .false.
                Me%Var%HNSVolatilization        = .false.
                Me%Var%HNSSedimentation         = .false.
                Me%Var%HNSDegradation           = .false.
                
                write(*,*) "Backtracking option is ON - all HNS processes were disconnected"
                write(*,*) "Subroutine HNSOptions - ModuleHNS - WRN010"  
            endif
                    
            Call CountLC50(ClientNumber)
            
        end if notcc

    end subroutine HNSOptions
    
    subroutine CountLC50(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer                                         :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, nLC50
        logical                                         :: LC50Found


        !Begin-----------------------------------------------------------------

        Me%Var%LC50 = .false.

        nLC50 = 0

DOLC50: do 

            call ExtractBlockFromBlockFromBlock(Me%ObjEnterData, ClientNumber,                   &
                                       LC50_begin, LC50_end,                            &
                                       LC50Found, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CountLC50 - ModuleHNS - ERR10'
            
i1:         if (LC50Found) then
                Me%Var%LC50 = .true.

                nLC50 = nLC50 + 1
 
            else i1
                call RewindBlockInBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'CountLC50 - ModuleHNS - ERR20'
                exit
            endif i1

        enddo DOLC50
        
        Me%Var%LC50Number = nLC50

        allocate(Me%Var%listLC50Specie(nLC50))
        allocate(Me%Var%listLC50Value(nLC50))
        allocate(Me%Var%listLC50Hours(nLC50))
        allocate(Me%Var%listLC50Type(nLC50))
        

        if (Me%Var%LC50) then
            
            call ReadLC50s(ClientNumber)
        endif

    end subroutine CountLC50
    !--------------------------------------------------------------------------    
    
    subroutine ReadLC50s(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer                                         :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, nLC50, flag
        logical                                         :: LC50Found
        character(len=StringLength)                     :: AuxChar
        

        !Begin-----------------------------------------------------------------

        call RewindBlockInBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'CountLC50 - ModuleHNS - ERR20'

DOLC50: do nLC50 = 1, Me%Var%LC50Number

            call ExtractBlockFromBlockFromBlock(Me%ObjEnterData, ClientNumber,          &
                                       LC50_begin, LC50_end,                            &
                                       LC50Found, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CountLC50s - ModuleHNS - ERR10'
            
i1:         if (LC50Found) then
 
                call GetData(AuxChar,                                                       &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             SearchType     = FromBlockInBlockInBlock,                      & 
                             keyword        = 'LC50_SPECIE',                                &
                             ClientModule   ='ModuleHNS',                                   &
                             STAT           = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR30'
                Me%Var%listLC50Specie(nLC50) = trim(adjustl(AuxChar))
                
                call GetData(Me%Var%listLC50Value(nLC50),                                   &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             SearchType     = FromBlockInBlockInBlock,                      & 
                             keyword        = 'LC50_VALUE',                                 &
                             ClientModule   ='ModuleHNS',                                   &
                             STAT           = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR40'

                call GetData(Me%Var%listLC50Hours(nLC50),                                   &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             SearchType     = FromBlockInBlockInBlock,                      & 
                             keyword        = 'LC50_HOURS',                                 &
                             ClientModule   ='ModuleHNS',                                   &
                             STAT           = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR50'
                call GetData(AuxChar,                                                       &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             SearchType     = FromBlockInBlockInBlock,                      & 
                             keyword        = 'LC50_TYPE',                                  &
                             ClientModule   ='ModuleHNS',                                   &
                             STAT           = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'HNSOptions - ModuleHNS - ERR60'
                select case (AuxChar)
                case("Surface", "surface", "SURFACE")
                    Me%Var%listLC50Type(nLC50) = LC50Surface_
                case("WaterColumn", "watercolumn", "WATERCOLUMN")
                    Me%Var%listLC50Type(nLC50) = LC50WaterColumn_
                case("Bottom", "bottom", "BOTTOM")
                    Me%Var%listLC50Type(nLC50) = LC50Bottom_
                case default
                    Me%Var%listLC50Type(nLC50) = LC50Surface_
                end select
            else i1
            
                stop 'CountLC50s - ModuleHNS - ERR40'

            endif i1

        enddo DOLC50

        call RewindBlockInBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'CountLC50s - ModuleHNS - ERR20'

    end subroutine ReadLC50s

    !--------------------------------------------------------------------------    
    
        subroutine ExportDropletOptions(HNSID, MethodBWDropletsDiameter, DropletsD50, STAT) 

        !Arguments---------------------------------------------------------------
        integer, intent(IN )                        :: HNSID     
        integer, intent(IN )                        :: MethodBWDropletsDiameter     
        real,    intent(IN)                         :: DropletsD50     
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_
        
        !--------------------------------------------------------------------------------     

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_) 

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Me%ExternalVar%MethodBWDropletsDiameter = MethodBWDropletsDiameter
            Me%ExternalVar%DropletsD50 = DropletsD50
            
            If ((MethodBWDropletsDiameter .NE. Computed_Classes_Random_) .AND. Me%Var%HNSEntrainment) then
                write(*,*)'Cannot compute entrainment with only one class of droplet diameters.'
                write(*,*)'Change keyword METHOD_BW_DROPLETS_DIAMETER to 3 (computed classes).'
                stop      'ExportDropletOptions - ModuleHNS - ERR10'
            Endif

            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
        
        end subroutine ExportDropletOptions

        !------------------------------------------------------------------------
    
    subroutine FindHNSBehaviourClass(VaporPressure,WaterSolubility,Density, WaterDensity, HNSBehaviourClass)

        !Arguments---------------------------------------------------------------
        real,           intent(IN )                           :: VaporPressure     
        real,           intent(IN )                           :: WaterSolubility     
        real,           intent(IN )                           :: Density     
        real, optional, intent(IN )                           :: WaterDensity     
        integer,        intent(OUT)                           :: HNSBehaviourClass

        !Local-----------------------------------------------------------------
        real                                                  :: WaterDensity_

        !--------------------------------------------------------------------------------     
        
        if (present(WaterDensity)) then
            WaterDensity_ = WaterDensity
        else
            WaterDensity_ = 1023.0
        end if

        if ((VaporPressure > 101.3E3) .and. (WaterSolubility <= 100000.)) then
                HNSBehaviourClass = Gas_
        elseif ((VaporPressure > 101.3E3) .and. (WaterSolubility > 100000.)) then
                HNSBehaviourClass = GasDissolver_
        elseif ((VaporPressure > 3.E3) .and. (WaterSolubility <= 10000.) .and. (Density <= WaterDensity_)) then
                HNSBehaviourClass = Evaporator_
        elseif ((VaporPressure > 3.E3) .and. ((WaterSolubility > 10000.) .and.       &
                (WaterSolubility <= 50000.)) .and. (Density <= 1025.)) then
                HNSBehaviourClass = EvaporatorDissolver_
        elseif ((VaporPressure > 0.3E3) .and. (VaporPressure <= 3.E3) .and.             &
                (WaterSolubility <= 1000.) .and. (Density <= 1025.)) then
                HNSBehaviourClass = FloaterEvaporator_
        elseif ((VaporPressure > 0.3E3) .and. (VaporPressure <= 3.E3) .and.             &
                (WaterSolubility > 1000.) .and. (WaterSolubility <= 50000.) .and. (Density <= WaterDensity_)) then
                HNSBehaviourClass = FloaterEvaporatorDissolver_
        elseif ((VaporPressure <= 0.3E3) .and. (WaterSolubility <= 1000.) .and. (Density <= WaterDensity_)) then
                HNSBehaviourClass = Floater_
        elseif ((VaporPressure <= 0.3E3) .and. ((WaterSolubility > 10000.) .and.    &
                (WaterSolubility <= 50000.)) .and. (Density <= 1025.)) then
                HNSBehaviourClass = FloaterDissolver_
        elseif ((VaporPressure > 10.E3) .and. (WaterSolubility > 50000.) .and. (Density <= WaterDensity_)) then
                HNSBehaviourClass = DissolverEvaporator_
        elseif ((VaporPressure <= 10.E3) .and. (WaterSolubility > 50000.)) then
                HNSBehaviourClass = Dissolver_
        elseif ((WaterSolubility > 1000.)) then
                HNSBehaviourClass = SinkerDissolver_
        elseif ((WaterSolubility <= 1000.)) then
                HNSBehaviourClass = Sinker_
        else
                HNSBehaviourClass = Unknown_
        end if

    end subroutine FindHNSBehaviourClass
    

    !----------------------------------------------------------------------------
    
    subroutine GetHNSInitialState(HNSID, HNSBehaviourClass, EmissionAtSurface, HNSInitialState, STAT)

        !Arguments---------------------------------------------------------------
        integer, intent (IN)                         :: HNSID
        integer, intent (IN)                         :: HNSBehaviourClass
        logical, intent (IN)                         :: EmissionAtSurface
        integer, intent(OUT)                         :: HNSInitialState
        integer, optional, intent(OUT)               :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_
        
        !--------------------------------------------------------------------------------     

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (HNSBehaviourClass .EQ. Gas_)    then

                If (EmissionAtSurface) then
                    HNSInitialState = Air_Evaporated_
                else
                    HNSInitialState = WaterColumn_Droplet_
                endif

            else
                If (EmissionAtSurface) then
                    HNSInitialState = Surface_
                else
                    HNSInitialState = WaterColumn_Droplet_
                endif
                
!            elseif ((HNSBehaviourClass .EQ. GasDissolver_) .OR.                         &
!                    (HNSBehaviourClass .EQ. Evaporator_) .OR.                           &
!                    (HNSBehaviourClass .EQ. EvaporatorDissolver_) .OR.                  &
!                    (HNSBehaviourClass .EQ. FloaterEvaporator_) .OR.                    &
!                    (HNSBehaviourClass .EQ. FloaterEvaporatorDissolver_) .OR.           &
!                    (HNSBehaviourClass .EQ. Floater_) .OR.                              &
!                    (HNSBehaviourClass .EQ. FloaterDissolver_) .OR.                     &
!                    (HNSBehaviourClass .EQ. DissolverEvaporator_)                       &
!                   )    then
!                  
!                If (EmissionAtSurface) then
!                    HNSInitialState = Surface_
!                else
!                    HNSInitialState = WaterColumn_Droplet_
!                endif
!
!            elseif (HNSBehaviourClass .EQ. Dissolver_) then
!
!                HNSInitialState = WaterColumn_Dissolved_
!
!            elseif (HNSBehaviourClass .EQ. SinkerDissolver_) then
!            
!                HNSInitialState = WaterColumn_Droplet_
!                
!            elseif (HNSBehaviourClass .EQ. Sinker_) then
!            
!                HNSInitialState = Bottom_Deposited_

            end if
        
            if ((HNSInitialState .EQ. WaterColumn_Droplet_) .AND. (Me%ExternalVar%DropletsD50 < 0.)) then
                write(*,*)'Median Particle Droplets (keyword "DROPLETS_D50") need to be defined'
                stop      'GetHNSInitialState - ModuleHNS - ERR10'
            end if
            
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetHNSInitialState
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetHNSAtmosphericDispersion(HNSID, HNSAtmosphericDispersion, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HNSID
        logical,              intent(OUT)           :: HNSAtmosphericDispersion
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            HNSAtmosphericDispersion = Me%Var%HNSAtmosphericDispersion          
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetHNSAtmosphericDispersion

    !--------------------------------------------------------------------------

    subroutine GetHNSSolarRadiationNeeded(HNSID, SolarRadiationNeeded, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HNSID
        logical,              intent(OUT)           :: SolarRadiationNeeded
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            SolarRadiationNeeded = Me%Var%SolarRadiationNeeded          
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetHNSSolarRadiationNeeded

    !----------------------------------------------------------------------

    subroutine GetHNSSedimentation(HNSID, HNSSedimentation, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HNSID
        logical,              intent(OUT)           :: HNSSedimentation
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            HNSSedimentation = Me%Var%HNSSedimentation          
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetHNSSedimentation

    !--------------------------------------------------------------------------

    subroutine GetHNSEntrainment(HNSID, HNSEntrainment, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HNSID
        logical,              intent(OUT)           :: HNSEntrainment
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            HNSEntrainment = Me%Var%HNSEntrainment          
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetHNSEntrainment

    !--------------------------------------------------------------------------

    subroutine GetHNSSpreading(HNSID, HNSSpreading, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HNSID
        logical,              intent(OUT)           :: HNSSpreading
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            HNSSpreading = Me%Var%HNSSpreading          
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetHNSSpreading

    !--------------------------------------------------------------------------

    subroutine GetHNSBehaviourClass(HNSID, HNSBehaviourClass, STAT)

        !Arguments-------------------------------------------------------------
        integer, intent (IN)                        :: HNSID
        integer, intent(OUT)                        :: HNSBehaviourClass
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            HNSBehaviourClass = Me%Var%HNSBehaviourClass          
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetHNSBehaviourClass

    !--------------------------------------------------------------------------

    subroutine GetHNSInitialMass(HNSID, Volume, InitialMass, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HNSID
        real,              intent(IN)               :: Volume
        real,              intent(OUT)              :: InitialMass
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        real                                        :: Density
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_) 

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            Density = Me%Var%Density
            InitialMass = Density * Volume
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetHNSInitialMass

    !--------------------------------------------------------------------------

    subroutine GetHNSDensity(HNSID, Density, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HNSID
        real,              intent(OUT)              :: Density
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_) 

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            Density = Me%Var%Density
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetHNSDensity

    !----------------------------------------------------------------------------

    subroutine GetInitialArea(HNSID, VolInic, AreaTotal, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: HNSID
        real, intent(IN)                    :: VolInic
        real, intent(OUT)                   :: AreaTotal
        integer, optional,intent(OUT)       :: STAT

        !Local-----------------------------------------------------------------
        integer                    :: ready_, STAT_
        !----------------------------------------------------------------------
        
        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !OILPOL METHOD
            !Rabeh, A.H., Lardner, R.W., Gunay, N. 2000           
            AreaTotal = Pi * (2.81 * sqrt(VolInic))**2.         

            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
        
    end subroutine GetInitialArea
    
        !--------------------------------------------------------------------------

    subroutine GetHNSSpreadingDiffVel(HNSID,                                       &
                                       DT,                                         &
                                       AreaTotal,                                  &
                                       AreaParticle,                               &
                                       VolumeParticle,                             &
                                       SpreadingDiffVel,                           &
                                       SpreadingDiffCoef,                          &
                                       STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: HNSID
        real, intent(IN)                   :: DT
        real, intent(IN)                   :: AreaTotal
        real, intent(IN)                   :: AreaParticle
        real, intent(IN)                   :: VolumeParticle
        real, intent(OUT)                  :: SpreadingDiffCoef
        real, intent(OUT)                  :: SpreadingDiffVel
        integer, optional, intent(OUT)     :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_
       
        call Ready(HNSID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then

        ! From Mackay et al., 1980 and Kolluru (1992)
        !Spreading diffusion coefficient is for each particle; main difference to Kolluru 1992
        !is that we are multiplying Mackay formulation with a ratio of areas instead a ratio of radius
        SpreadingDiffCoef = SpreadingRateConstant * AreaParticle**(1./3.) * (VolumeParticle / AreaParticle)**(4./3.)    &
                            * (AreaParticle/AreaTotal)**(2./3.)

        SpreadingDiffVel  = sqrt(6.0 * SpreadingDiffCoef / DT)
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetHNSSpreadingDiffVel
    
    !--------------------------------------------------------------------------

    subroutine GetHNSPAC(HNSID, PAC, PAC1, PAC2, PAC3, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HNSID
        logical,              intent(OUT)           :: PAC
        integer,              intent(OUT)           :: PAC1
        integer,              intent(OUT)           :: PAC2
        integer,              intent(OUT)           :: PAC3
        integer, optional,    intent(OUT)           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            PAC  = Me%Var%PAC          
            PAC1 = Me%Var%PAC1          
            PAC2 = Me%Var%PAC2          
            PAC3 = Me%Var%PAC3          
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !--------------------------------------------------------------------------

    end subroutine GetHNSPAC

    !--------------------------------------------------------------------------

    subroutine GetHNSPNEC(HNSID, PNEC, STAT)

        !Arguments-------------------------------------------------------------
        integer                                       :: HNSID
        real,                   intent(OUT)           :: PNEC
        integer, optional,      intent(OUT)           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            PNEC  = Me%Var%PNEC          
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetHNSPNEC

        !----------------------------------------------------------------------

    subroutine GetHNSLC50(HNSID, LC50, LC50Number, LC50Species, LC50Values, LC50Types, LC50Hours, STAT)

!        Arguments-------------------------------------------------------------
        integer                                             :: HNSID
        logical,              intent(OUT)                    :: LC50
        integer,              intent(OUT)                    :: LC50Number
        character(LEN = StringLength), dimension(:), pointer  :: LC50Species
        real, dimension(:), pointer                          :: LC50Values
        integer, dimension(:), pointer                      :: LC50Types
        real, dimension(:), pointer                          :: LC50Hours
        integer, optional,    intent(OUT)                   :: STAT
 !       integer                                             :: nLC50

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

       ! ----------------------------------------------------------------------

        STAT_ = UNKNOWN_
        
        nullify (LC50Species)
        nullify (LC50Values)
        nullify (LC50Types)
        nullify (LC50Hours)

        call Ready(HNSID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.           &
            (ready_ .EQ. READ_LOCK_ERR_)) then
    
            !allocate(LC50Species(1:LC50Number))
            !allocate(LC50Values(1:LC50Number))
            !allocate(LC50Types(1:LC50Number))
            !allocate(LC50Hours(1:LC50Number))

            LC50               = Me%Var%LC50
            LC50Number         = Me%Var%LC50Number
            LC50Values(1:LC50Number)  => Me%Var%listLC50Value         
            LC50Species(1:LC50Number) => Me%Var%listLC50Specie         
            LC50Types(1:LC50Number)   => Me%Var%listLC50Type
            LC50Hours(1:LC50Number)   => Me%Var%listLC50Hours         

            !nLC50 = 0
            !do nLC50 = 1, Me%Var%LC50Number 
            !    LC50Values(nLC50)  = Me%Var%listLC50Value(nLC50)         
            !    LC50Species(nLC50) = Me%Var%listLC50Specie(nLC50)         
            !    LC50Types(nLC50)   = Me%Var%listLC50Type(nLC50)         
            !    LC50Hours(nLC50)   = Me%Var%listLC50Hours(nLC50)         
            !enddo            
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

     !   ----------------------------------------------------------------------

    end subroutine GetHNSLC50

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine HNSInternalProcesses(HNSID,                                        &
                                    Wind,                                         &
                                    Currents,                                     &
                                    WaterTemperature,                             &
                                    WaterDensity,                                 &
                                    SPM,                                          &
                                    AirTemperature,                               &
                                    AtmPressure,                                  &
                                    WaveHeight,                                   &
                                    WavePeriod,                                   &
                                    Area,                                         &
                                    InitialMass,                                  &
                                    MassIN,                                       &
                                    HNSParticleStateIN,                           &
                                    DropletsDiameterIN,                           &
                                    Depth,                                        &
                                    Density,                                      &
                                    MassOUT,                                      &
                                    VolumeOUT,                                    &
                                    DropletsDiameterOUT,                          &
                                    HNSParticleStateOUT,                          &
                                    STAT)
                               
        !Arguments---------------------------------------------------------------
        integer                                     :: HNSID
        real,              intent(IN )              :: Wind
        real,              intent(IN )              :: Currents
        real,              intent(IN )              :: WaterTemperature
        real,              intent(IN )              :: WaterDensity
        real,              intent(IN )              :: SPM
        real,              intent(IN)               :: AirTemperature
        real,              intent(IN)               :: AtmPressure       
        real,              intent(IN)               :: WaveHeight       
        real,              intent(IN)               :: WavePeriod       
        real,              intent(IN)               :: Area       
        real,              intent(IN)               :: InitialMass
        real,              intent(IN)               :: MassIN
        integer,           intent(IN)               :: HNSParticleStateIN
        real,              intent(IN)               :: DropletsDiameterIN
        real,              intent(IN)               :: Depth
        real,              intent(OUT)              :: Density
        real,              intent(OUT)              :: MassOUT
        real,              intent(OUT)              :: VolumeOUT
        real,              intent(OUT)              :: DropletsDiameterOUT
        integer,           intent(OUT)              :: HNSParticleStateOUT
        integer, optional, intent(OUT)              :: STAT


        !Local-------------------------------------------------------------------
        integer                                     :: ready_ 
        integer                                     :: STAT_ 

        !Parameter---------------------------------------------------------------
        character(LEN = StringLength), parameter :: Char_Gas                          = 'Gas (GD)'
        character(LEN = StringLength), parameter :: Char_GasDissolver                 = 'Gas - Dissolver (GD)'
        character(LEN = StringLength), parameter :: Char_Evaporator                   = 'Evaporator (E)'
        character(LEN = StringLength), parameter :: Char_EvaporatorDissolver          = 'Evaporator - Dissolver (ED)'
        character(LEN = StringLength), parameter :: Char_FloaterEvaporator            = 'Floater - Evaporator (FE)'
        character(LEN = StringLength), parameter :: Char_FloaterEvaporatorDissolver   = 'Floater - Evaporator - Dissolver (FDE)'
        character(LEN = StringLength), parameter :: Char_Floater                      = 'Floater (F)'
        character(LEN = StringLength), parameter :: Char_FloaterDissolver             = 'Floater - Dissolver (FD)'
        character(LEN = StringLength), parameter :: Char_DissolverEvaporator          = 'Dissolver - Evaporator (DE)'
        character(LEN = StringLength), parameter :: Char_Dissolver                    = 'Dissolver (D)'
        character(LEN = StringLength), parameter :: Char_SinkerDissolver              = 'Sinker - Dissolver (SD)'
        character(LEN = StringLength), parameter :: Char_Sinker                       = 'Sinker (S)'
        character(LEN = StringLength), parameter :: Char_Unknown                      = 'Unknown'

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_)
        
cd1 :   if (ready_ .EQ. IDLE_ERR_) then

!cd2 :       if (LagrangianTime .GE. Me%NextInternalComputeTime) then

!                do while (LagrangianTime .GE. Me%NextInternalComputeTime) 
!                    Me%Now                      = Me%NextInternalComputeTime            
!                    Me%Var%Time                 = Me%Var%Time + Me%Var%DTHNSInternalProcesses
        
!                    !Output into an ASCII file 
!    cd11 :          if (Me%State%FirstStepIP .and. Me%State%TimeSerie) then
!                        call TimeSerieOutput(CurrentTime = Me%BeginTime)
!
!                    end if cd11          

                    Me%ExternalVar%Wind                 = Wind
                    Me%ExternalVar%AirTemperature       = AirTemperature
                    Me%ExternalVar%AtmPressure          = AtmPressure
                    Me%ExternalVar%WaveHeight           = WaveHeight
                    Me%ExternalVar%WavePeriod           = WavePeriod
                    Me%ExternalVar%Area                 = Area
                    Me%Var%HNSParticleState             = HNSParticleStateIN
                    Me%ExternalVar%WaterTemperature     = WaterTemperature
                    Me%ExternalVar%WaterDensity         = WaterDensity
                    Me%ExternalVar%Wind                 = Wind
                    Me%ExternalVar%Currents             = Currents
                    Me%ExternalVar%SPM                  = SPM
                    Me%ExternalVar%Depth                = Depth
                    Me%Var%InitialMass                  = InitialMass
                    Me%Var%Mass                         = MassIN
                    Me%Var%DropletDiameter              = DropletsDiameterIN
                    
                    Call InitializeVariables
                    
                    if (Me%State%FirstStepIP) then
                        if (Me%Var%HNSGas) then
                            Me%Var%HNSBehaviourClass = Gas_
                            Me%Var%HNSParticleState = Air_Evaporated_
                        else
                        
                            !Find HNS behaviour class
                            call FindHNSBehaviourClass(VaporPressure     = Me%Var%VaporPressure,        &
                                                      WaterSolubility    = Me%Var%WaterSolubility,      &
                                                      Density            = Me%Var%Density,              &
                                                      WaterDensity       = Me%ExternalVar%WaterDensity, &
                                                      HNSBehaviourClass  = Me%Var%HNSBehaviourClass     )
                        endif                        
                        select case(Me%Var%HNSBehaviourClass)
                        case(Gas_)
                            write(*,*) 'The chemical used is classified as ',  Char_Gas                

                        case(GasDissolver_)
                            write(*,*) 'The chemical used is classified as ',  Char_GasDissolver                        

                        case(Evaporator_)
                            write(*,*) 'The chemical used is classified as ',  Char_Evaporator
                
                        case(EvaporatorDissolver_)
                            write(*,*) 'The chemical used is classified as ',  Char_EvaporatorDissolver  
                
                        case(FloaterEvaporator_)
                            write(*,*) 'The chemical used is classified as ',  Char_FloaterEvaporator
                
                        case(FloaterEvaporatorDissolver_)
                            write(*,*) 'The chemical used is classified as ',  Char_FloaterEvaporatorDissolver            
                
                        case(Floater_)
                            write(*,*) 'The chemical used is classified as ',  Char_Floater   
                
                        case(FloaterDissolver_)
                            write(*,*) 'The chemical used is classified as ',  Char_FloaterDissolver
                
                        case(DissolverEvaporator_)
                            write(*,*) 'The chemical used is classified as ',  Char_DissolverEvaporator          
                
                        case(Dissolver_)
                            write(*,*) 'The chemical used is classified as ',  Char_Dissolver    
                
                        case(SinkerDissolver_)
                            write(*,*) 'The chemical used is classified as ',  Char_SinkerDissolver             
                
                        case(Sinker_)
                            write(*,*) 'The chemical used is classified as ',  Char_Sinker            
                
                        case(ClassUnknown_)
                            write(*,*) 'The chemical used is classified as ',  Char_Unknown
                
                        end select
                    endif
                    
                    if (Me%Var%HNSGas) then
                        Me%Var%HNSBehaviourClass = Gas_
                        Me%Var%HNSParticleState = Air_Evaporated_
                    endif
                    
                    if (Me%Var%HNSParticleState .EQ. Surface_) then
                        if (Me%Var%HNSEvaporation) call Evaporation
                        if (Me%Var%HNSEntrainment) call Entrainment
                        ! no surface dissolution to gases discharged at surface
                        if ((Me%Var%HNSBehaviourClass .NE. Gas_) .AND. &
                            (Me%Var%HNSBehaviourClass .NE. GasDissolver_)) then
                            if (Me%Var%HNSDissolution) then
                                call Dissolution
                                if (Me%Var%HNSSedimentation) call Sedimentation
                            endif
                        endif
                    elseif (Me%Var%HNSParticleState .EQ. WaterColumn_Droplet_) then
                        !check if ParticleDropletDiameter is already initialized
                        if (Me%Var%DropletDiameter < 0.) call FindInitialDropletDiameter
                        !check if particle is dissolved (< 0.2um)
                        if (Me%Var%DropletDiameter .LT. 2.0E-7) then
                            Me%Var%HNSParticleState = WaterColumn_Dissolved_
                        else
                            if (Me%Var%HNSDissolution) then
                                call Dissolution
                                if (Me%Var%HNSSedimentation) call Sedimentation
                            endif
                        endif
                    elseif (Me%Var%HNSParticleState .EQ. WaterColumn_Dissolved_) then
                        if (Me%Var%HNSVolatilization) call Volatilization
                    endif
                    
                    call UpdateHNSParticleState
                    HNSParticleStateOUT = Me%Var%HNSParticleState                   
                    
                    if (Me%Var%HNSDegradation) call Degradation
                    MassOUT   = max(0., MassIN - Me%Var%MDegradedDT)
                    
                    call DegradeDropletDiameter(MassIN, MassOUT)
                    
                    Density   = Me%Var%Density
                    VolumeOUT = MassOUT / Me%Var%Density 
                    DropletsDiameterOUT = Me%Var%DropletDiameter
!                    !Output into an ASCII file 
!    cd3 :           if (Me%State%TimeSerie) then
!                        call TimeSerieOutput
!                    end if cd3
                    
!                    Me%NextInternalComputeTime = Me%NextInternalComputeTime + Me%Var%DTHNSInternalProcesses

!                enddo 

!            else cd2

!            end if cd2

            Me%State%FirstStepIP = OFF

            STAT_ = SUCCESS_
            
        else cd1
         
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                                          &
            STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine HNSInternalProcesses

    !----------------------------------------------------------------------------

    subroutine InitializeVariables
            Me%Var%MEvaporatedDT   = 0.0
            Me%Var%MEntrainedDT    = 0.0
            Me%Var%MDissolvedDT    = 0.0
            Me%Var%MVolatilizedDT  = 0.0
            Me%Var%MSedimentedDT   = 0.0   
            Me%Var%MDegradedDT     = 0.0   
    end subroutine InitializeVariables

    !----------------------------------------------------------------------------

    subroutine UpdateHNSParticleState
    
        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        real :: RAND
        real :: InitialMass
        real :: MEvaporatedDT, MEntrainedDT
        real :: SumFractionEvaporatedDT, SumFractionEntrainedDT
        real :: MDissolvedDT, MVolatilizedDT, MSedimentedDT
        real :: SumFractionDissolvedDT, SumFractionVolatilizedDT, SumFractionSedimentedDT
        integer :: HNSParticleState
        !real:: MinimumProbability
        !integer::SampleSize
        
        !------------------------------------------------------------------------
        MEvaporatedDT = 0.0
        MEntrainedDT  = 0.0
        MDissolvedDT  = 0.0
        MVolatilizedDT = 0.0
        MVolatilizedDT = 0.0
        HNSParticleState = Me%Var%HNSParticleState
        If (Me%Var%HNSEvaporation)      MEvaporatedDT = Me%Var%MEvaporatedDT
        If (Me%Var%HNSEntrainment)      MEntrainedDT = Me%Var%MEntrainedDT
        If (Me%Var%HNSDissolution)      MDissolvedDT = Me%Var%MDissolvedDT
        If (Me%Var%HNSVolatilization)   MVolatilizedDT = Me%Var%MVolatilizedDT
        If (Me%Var%HNSSedimentation)    MSedimentedDT = Me%Var%MSedimentedDT

        InitialMass = Me%Var%InitialMass

        Call RANDOM_NUMBER(RAND)

        select case(HNSParticleState)
        case(Surface_)
        
            SumFractionEvaporatedDT = MEvaporatedDT / InitialMass
            SumFractionEntrainedDT = (MEvaporatedDT + MEntrainedDT) / InitialMass
            SumFractionDissolvedDT = (MEvaporatedDT + MEntrainedDT + MDissolvedDT) /InitialMass
            SumFractionSedimentedDT = (MEvaporatedDT + MEntrainedDT + MDissolvedDT + MSedimentedDT) /InitialMass
            !MinimumProbability = min(SumFractionEvaporatedDT, SumFractionEntrainedDT, SumFractionDissolvedDT, &
            !                         SumFractionSedimentedDT)
            !SampleSize = nint(1./ MinimumProbability)
            
            
            if      (RAND < SumFractionEvaporatedDT) then

                HNSParticleState = Air_Evaporated_
                
            else if (RAND < SumFractionEntrainedDT .and. RAND >= SumFractionEvaporatedDT) then

                HNSParticleState = WaterColumn_Droplet_

            else if (RAND < SumFractionDissolvedDT .and. RAND >= SumFractionEntrainedDT) then
        
                HNSParticleState = WaterColumn_Dissolved_

            else if (RAND < SumFractionSedimentedDT .and. RAND >= SumFractionDissolvedDT) then
        
                HNSParticleState = WaterColumn_Sedimented_
            else
                HNSParticleState = Surface_
            endif
                       
        case (WaterColumn_Droplet_) 

            SumFractionDissolvedDT = MDissolvedDT / InitialMass
            SumFractionSedimentedDT = (MDissolvedDT + MSedimentedDT) / InitialMass

            if      (RAND < SumFractionDissolvedDT) then

                HNSParticleState = WaterColumn_Dissolved_
                
            else if (RAND < SumFractionSedimentedDT .and. RAND >= SumFractionDissolvedDT) then

                HNSParticleState = WaterColumn_Sedimented_
            else
            
                HNSParticleState = WaterColumn_Droplet_
            endif

        case (WaterColumn_Dissolved_) 

            SumFractionVolatilizedDT = MVolatilizedDT / InitialMass

            if      (RAND < SumFractionVolatilizedDT) then

                HNSParticleState = Air_Volatilized_                
            else
                HNSParticleState = WaterColumn_Dissolved_
            endif

        case (WaterColumn_Sedimented_)
                           
        case (Bottom_Deposited_)
            If (Me%ExternalVar%Currents > 0.2)  then
                HNSParticleState = WaterColumn_Sedimented_
            endif

        end select
                    
        Me%Var%HNSParticleState = HNSParticleState
    end subroutine UpdateHNSParticleState

    
        !----------------------------------------------------------------------------

    subroutine DegradeDropletDiameter(MassIN, MassOUT)

        !Arguments---------------------------------------------------------------
        real, intent(IN) :: MassIN
        real, intent(IN) :: MassOUT
        
        !Local-------------------------------------------------------------------
        real            :: MassRatio
        real            :: DropletIndividualMassIN
        real            :: DropletIndividualMassOUT
        real            :: DropletDiameterIN
        real            :: DropletDiameterOUT
        real            :: aux
    !----------------------------------------------------------------------------
        DropletDiameterIN           = Me%Var%DropletDiameter
        MassRatio                   = MassOUT / MassIN
        aux                         = Me%Var%Density * Pi / 6
        
        DropletIndividualMassIN     = aux * DropletDiameterIN**3
        DropletIndividualMassOUT    = MassRatio * DropletIndividualMassIN
        DropletDiameterOUT          = (DropletIndividualMassOUT / aux)**(1./3.)
        Me%Var%DropletDiameter      = DropletDiameterOUT
        
    end subroutine DegradeDropletDiameter
    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    subroutine Evaporation

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------

        real :: KMassTransf
        real :: Wind
        real :: Diameter
        real :: Dm 
        real :: MolecularWeight
        real :: SchmidtNbr
        real :: EvaporationRate
        real :: CorrectedEvaporationRate
        real :: Area     
        real :: VaporPressure
        real :: AirTemperatureInKelvin 
        real :: AtmPressure
        real :: MEvaporated
        real :: MEvaporatedDT
        real :: VolatileCorrection
        !------------------------------------------------------------------------

        Area = Me%ExternalVar%Area 
        Diameter = sqrt( (4.*Area)/ Pi )
        Wind = Me%ExternalVar%Wind
        MolecularWeight = Me%Var%MolecularWeight
        VaporPressure = Me%Var%VaporPressure
        AirTemperatureInKelvin = Me%ExternalVar%AirTemperature + AbsoluteZero
        AtmPressure = Me%ExternalVar%AtmPressure
        MEvaporated = Me%Var%MEvaporated
        !Molecular Diffusion Coef - Graham's Law - Thibodeaux 1979
        Dm = WaterDiffusionCoef * (WaterMolecularWeight / MolecularWeight) **(1./2.)

        !Schmidt number, unitless ratio
        SchmidtNbr = AirCinematicVisc / Dm
        
        !Mass transfer coefficient, Mackay and Matsugu, 1973 - units m/s
        KMassTransf = 0.0048 * Wind **(7./9.) * Diameter **(-1./9.) * SchmidtNbr **(-2./3.)

        ! evaporation rate, Kawamura and Mackay 1985 - units kg/s
        !next R is multiplied by 1.e3 in order to get R in units J/kmol.K
        EvaporationRate = Area * KMassTransf * ((MolecularWeight * VaporPressure) / (R * 1.e3 * AirTemperatureInKelvin))
        
        !Correction for Volatilization (Brighton, 1985, 1990; Reynolds 1992)
        VolatileCorrection = - (AtmPressure / VaporPressure) * log(Max(AlmostZero,(1. - (VaporPressure / AtmPressure ))))
        
        CorrectedEvaporationRate = EvaporationRate * VolatileCorrection
        
        MEvaporatedDT = CorrectedEvaporationRate * Me%Var%DTHNSInternalProcesses
        Me%Var%MEvaporatedDT = MEvaporatedDT
        
    end subroutine Evaporation
    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    ! Surface_ entrainment due to breaking waves
    subroutine Entrainment

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        real            :: Area
        real            :: ViscCin
        real            :: QdTotal
!        real            :: DropletDiameter(5), Qd(5)
        real            :: WaveHeight
        real            :: WaterDensity
        real            :: Wind
        real            :: WavePeriod
        real            :: DropletsD50
        real            :: DropletDiameter
        
        

        WaveHeight = Me%ExternalVar%WaveHeight
        WaterDensity = Me%ExternalVar% WaterDensity
        Wind = Me%ExternalVar%Wind
        WavePeriod = Me%ExternalVar%WavePeriod
        DropletsD50 = Me%ExternalVar%DropletsD50
        Area = Me%ExternalVar%Area 
        ViscCin = Me%Var%ViscCin
        
        call GetDropletDiameterOrQdTotal(&
              MethodBWDropletsDiameter    = Computed_Classes_Random_, &
              D50                         = DropletsD50, &
              ParticleViscCin             = ViscCin, &
              WaveHeight                  = WaveHeight, &
              WaterDensity                = WaterDensity, &
              Wind                        = Wind, &
              WavePeriod                  = WavePeriod, &
              DropletDiameter             = DropletDiameter, &
              QdTotal                     = QdTotal)

        Me%Var%DropletDiameter = DropletDiameter
        
        Me%Var%MEntrainedDT           = QdTotal * Area * Me%Var%DTHNSInternalProcesses

    end subroutine Entrainment
    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    ! Dissolution from entrained chemical at water column
    subroutine Dissolution

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        real :: Diffusivity
        real :: WaterTemperatureInKelvin
        real :: WaterCinVisc_cSt
        real :: WaterDynamicVisc_
        real :: WaterSolubility
        real :: MolarVolume_       
        real :: MolarVolumeLeBas_
        real :: Density
        real :: KTransferMass
        real :: WaterTemperature
        real :: Diameter, SurfaceArea, ReynoldsVelocity, ReynoldsNbr, SherwoodNbr, SchmidtNbr, WaterDensity
!        real  :: R_converted
        integer :: ParticleState
        real :: MDissolvedDT 
        
        ParticleState = Me%Var%HNSParticleState
        
        !Density in g/cm3
        Density = Me%Var%Density / 1000.
        
        !Molar Volume in cm3/mol
        MolarVolumeLeBas_ = MolarVolumeLeBas(Me%Var%Organic, Me%Var%MolecularWeight)
        MolarVolume_      = Me%Var%MolecularWeight / Density
        
        WaterTemperature = Me%ExternalVar%WaterTemperature
        WaterTemperatureInKelvin = Me%ExternalVar%WaterTemperature + AbsoluteZero

        WaterDensity = Me%ExternalVar%WaterDensity
        
        !Dynamic Viscocity in cP (mPa/s)
        WaterDynamicVisc_ = GetWaterDynamicVisc(WaterTemperature)
        !Water Kynematic Visc in cSt (mm2/s)
        WaterCinVisc_cSt = WaterDynamicVisc_ / WaterDensity
        
        WaterSolubility = Me%Var%WaterSolubility

        !universal gas constant R in atm-m3/mol-K
        !R_converted = 8.206e-5

        !Diffusivity (Hayduk & Minhas, 1982) in cm2/s
        !Diffusivity = 1.25e-8 * ((MolarVolume_**(-0.19)) - 0.292) * (WaterTemperatureInKelvin **1.52) * WaterDynamicVisc_

        !Diffusivity (Hayduk & Laudie method) - in (Lyman et al., 1982 apud Hines & Maddox, 1985)  in cm2/s
        !Diffusivity = (R_converted * 1.E-5 * WaterTemperatureInKelvin) / ((WaterDynamicVisc_**1.14) * (MolarVolumeLeBas_**(0.589)))

        !Diffusivity (Hayduk & Laudie method) - in (Hayduk, W. and H. Laudie, 1974)  in cm2/s
        Diffusivity = (13.26 * 1.E-5) / ((WaterDynamicVisc_**1.14) * (MolarVolumeLeBas_**(0.589)))
        
        SchmidtNbr = WaterCinVisc_cSt / Diffusivity

        If (ParticleState .EQ. Surface_) then

            !Surface_ Slick Dissolution = Disc Plate
            
            SurfaceArea = Me%ExternalVar%Area 
            
            Diameter = sqrt( (4.*SurfaceArea)/ Pi )
           
            ReynoldsVelocity = Me%ExternalVar%Wind 
                      
            ReynoldsNbr = ReynoldsVelocity * Diameter / (WaterCinVisc_cSt*1e-6)

            SherwoodNbr = .578 * (SchmidtNbr **(1./3.)) * (ReynoldsNbr ** 0.5)
                 
        Else
        
            !Dissolution at water column = Sphere

            Diameter = Me%Var%DropletDiameter
            SurfaceArea = Pi * (Diameter **2.)
            
            ReynoldsVelocity = Me%ExternalVar%Currents 

            ReynoldsNbr = ReynoldsVelocity * Diameter / (WaterCinVisc_cSt*1.e-6)

            SherwoodNbr = 2. + 0.347 * (ReynoldsNbr**(0.62)) * (SchmidtNbr **(0.31))

        EndIf        

        !next 1.e-4 is to convert diffusivity units from cm2/s to m2/s
        KTransferMass = (SherwoodNbr * Diffusivity * 1.e-4) / Diameter  

        !next 1.e-3 is to convert mass units, because WaterSolubility is in mg/l = g/m3
        !and it should be in kg/m3 to get MDissolvedDT in kg/s
        MDissolvedDT = KTransferMass * WaterSolubility * 1.e-3 * SurfaceArea * Me%Var%DTHNSInternalProcesses
        
        Me%Var%MDissolvedDT = MDissolvedDT
!        Me%Var%MDissolved = max(MDissolved - MVolatilizedDT * Me%Var%DTHNSInternalProcesses, 0.)
        
    end subroutine Dissolution

    !----------------------------------------------------------------------------

    ! Droplet Diameter for particles initially entrained
    subroutine FindInitialDropletDiameter

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        real            :: Area
        real            :: ViscCin
        real            :: WaveHeight
        real            :: WaterDensity
        real            :: Wind
        real            :: WavePeriod
        real            :: DropletsD50
        real            :: DropletDiameter
        integer         :: MethodBWDropletsDiameter
        
        

        WaveHeight = Me%ExternalVar%WaveHeight
        WaterDensity = Me%ExternalVar% WaterDensity
        Wind = Me%ExternalVar%Wind
        WavePeriod = Me%ExternalVar%WavePeriod
        DropletsD50 = Me%ExternalVar%DropletsD50
        Area = Me%ExternalVar%Area 
        ViscCin = Me%Var%ViscCin
        MethodBWDropletsDiameter = Me%ExternalVar%MethodBWDropletsDiameter
        
        call GetDropletDiameterOrQdTotal(&
              MethodBWDropletsDiameter    = MethodBWDropletsDiameter, &
              D50                         = DropletsD50, &
              ParticleViscCin             = ViscCin, &
              WaveHeight                  = WaveHeight, &
              WaterDensity                = WaterDensity, &
              Wind                        = Wind, &
              WavePeriod                  = WavePeriod, &
              DropletDiameter             = DropletDiameter)

        Me%Var%DropletDiameter = DropletDiameter
        

    end subroutine FindInitialDropletDiameter
    !----------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    real function MolarVolumeLeBas(Organic, MolecularWeight)

        !Arguments-------------------------------------------------------------

        logical :: Organic
        real    :: MolecularWeight
        
        !---------------------------------------------------------------------
        if (Organic) then
            MolarVolumeLeBas = 4.9807 * (MolecularWeight ** 0.6963)
        else
            MolarVolumeLeBas = 2.8047 * (MolecularWeight ** 0.651)
        endif
    
    end function MolarVolumeLeBas
    
    !----------------------------------------------------------------------------
    
    ! Water Dynamic Viscosity (cP) -> Van Veltzen et al., Yaws et.al, and Duhne (Reid, 1987)
    real function GetWaterDynamicVisc(WaterTemperature)

        !Arguments-------------------------------------------------------------

        real    :: WaterTemperature
        
        !Local-------------------------------------------------------------
        
        real    :: WaterTemperatureInKelvin
        !---------------------------------------------------------------------

        WaterTemperatureInKelvin = WaterTemperature + AbsoluteZero
        GetWaterDynamicVisc = exp(-24.71 + (4209. /WaterTemperatureInKelvin) + 0.04527 * WaterTemperatureInKelvin &
                                - 3.376e-5 * (WaterTemperatureInKelvin ** 2.)) 
    
    end function GetWaterDynamicVisc
    
    !----------------------------------------------------------------------------
   
    ! volatilization from the water column
    subroutine Volatilization

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------

        real :: H
        real :: H_Dimensionless
        real :: MolecularWeight
        real :: VaporPressure
        real :: AirTemperatureInKelvin 
        real :: WaveHeight
        real :: WaterSolubility
        real :: K5, K6, K7       
        real :: MixDepth
        real :: MVolatilizedDT
        real :: MDissolved
        real :: Depth

        !------------------------------------------------------------------------

        MolecularWeight = Me%Var%MolecularWeight
        VaporPressure = Me%Var%VaporPressure
        AirTemperatureInKelvin = Me%ExternalVar%AirTemperature + AbsoluteZero
        WaterSolubility = Me%Var%WaterSolubility
        MDissolved = Me%Var%Mass
        WaveHeight = Me%ExternalVar%WaveHeight
        Depth = Me%ExternalVar%Depth
        
        !next 1.e-3 is to convert solubility units to kg/m3, in order to get H in 
        ! m3*Pa/kmol or J/kmol
        H = VaporPressure / (WaterSolubility * 1.e-3 / MolecularWeight)
        
if1:    If ( H < 30.3975) then
            MVolatilizedDT = 0
        Else if1

            MixDepth = WaveHeight / 2.

            If (Depth > MixDepth) then
                MVolatilizedDT = 0
            Else
                !next 1.e3 is to put R in J/kmol.K
                H_Dimensionless = H / (R * 1.e3 * AirTemperatureInKelvin)
                
                !liquid-phase exchange coefficient (m/s)
                K5 = 20. * sqrt(44./ MolecularWeight) * (0.01 / 3600.) 

                !gas-phase exchange coefficient (m/s)
                K6 = 3000. * sqrt (18. / MolecularWeight)  * (0.01 / 3600.)
                
                !overall mass transfer coefficient (m/s)
                K7 = (H_Dimensionless * K5 * K6)/(H_Dimensionless * K6 + K5)
                MVolatilizedDT = min(MDissolved, (K7 * MDissolved / MixDepth) * Me%Var%DTHNSInternalProcesses)
            End If
        End If if1

        Me%Var%MVolatilizedDT = MVolatilizedDT
!        Me%Var%MVolatilized = Me%Var%MVolatilized + MVolatilizedDT * Me%Var%DTHNSInternalProcesses

!        Me%Var%MDissolved = max(MDissolved - MVolatilizedDT * Me%Var%DTHNSInternalProcesses, 0.)
        
    end subroutine Volatilization
    
    !----------------------------------------------------------------------------
    
    !----------------------------------------------------------------------------
    ! Adsorption of dissolved particles to suspended sediments at water column
    subroutine Sedimentation()

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        
        real :: SedimentsConcentration
        real :: logKow, Koc, MSedimentedDT, MDissolvedDT
        real :: SedimentedFactor, ActualMassDissolvedDT
        !------------------------------------------------------------------------

        SedimentsConcentration = Me%ExternalVar%SPM
        MDissolvedDT = Me%Var%MDissolvedDT

        !octanol-water partition coefficient        
        logKow = Me%Var%OWPartitionCoef

        !Sorption Coefficient
        Koc = 10**(0.983 * logKow + 0.00028)
        
        !Cads = Koc * Css * Cdiss => We want to know fraction adsorbed (Cads / (Cdiss + Cads))
        
!        Volume = MDissolvedDT / Me%Var%Density
        
!        !Mass of Sediments = SPM available in a particle volume 
!        SedimentsMass = (SedimentsConcentration * 1.E-3) / Volume
!
!        !1 - find sedimentation factor to multiply by dissolved mass
!        SedimentedFactor  = Koc * SedimentsMass 
        

        !1 - find sedimentation factor to multiply by dissolved mass
        SedimentedFactor  = Koc * (SedimentsConcentration * 1.E-3)

        !2 - find actual dissolved mass in equilibrium with sedimented mass
        !(dissolved mass is conditioned by the part sedimented)
        !SedimentedFactor * MDissolvedDT + MDissolvedDT = InitialMass
        ActualMassDissolvedDT = MDissolvedDT / (SedimentedFactor + 1)

        !3 - Update mass dissolved and sedimented, based on actual dissolved mass
        MSedimentedDT  =   SedimentedFactor  * ActualMassDissolvedDT
        MDissolvedDT   =   ActualMassDissolvedDT

        Me%Var%MSedimentedDT = MSedimentedDT
        Me%Var%MDissolvedDT   = MDissolvedDT

        Me%Var%MSedimentedDT = MSedimentedDT

    end subroutine Sedimentation

    !----------------------------------------------------------------------------

    ! degradation of mass for every parcel
    !routine reduces mass of every parcel
    subroutine Degradation

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------

        real :: AirDegradationRate, WaterDegradationRate, SedimentDegradationRate
        real :: DT
        real :: AirDegradationFactor, WaterDegradationFactor, SedimentDegradationFactor
        real :: MDegradedDT, MEntrained, MDissolved, MEvaporated, MVolatilized
        real :: MSedimented
        real :: HNSParticleState
        real :: Mass
        !------------------------------------------------------------------------
        
        Mass = Me%Var%Mass
        HNSParticleState = Me%Var%HNSParticleState
        !Degradation Rate in s-1
        AirDegradationRate =        Me%Var%AirDegradationRate / 86400.
        WaterDegradationRate =      Me%Var%WaterDegradationRate / 86400.
        SedimentDegradationRate =   Me%Var%SedimentDegradationRate / 86400.
        
        DT = Me%Var%DTHNSInternalProcesses
        MEntrained = Me%Var%MEntrained
        MDissolved = Me%Var%MDissolved
        MEvaporated = Me%Var%MEvaporated
        MVolatilized = Me%Var%MVolatilized
        MSedimented = Me%Var%MSedimented
        MDegradedDT = Me%Var%MDegradedDT
        
        AirDegradationFactor = exp(-AirDegradationRate * DT)
        WaterDegradationFactor = exp(-WaterDegradationRate * DT)
        SedimentDegradationFactor = exp(-SedimentDegradationRate * DT)
                     
        If ((HNSParticleState .EQ. Air_Evaporated_) .OR.            &
            (HNSParticleState .EQ. Air_Volatilized_)) then
            MDegradedDT = (1. - AirDegradationFactor) * Mass
        ElseIf (HNSParticleState .EQ. Bottom_Deposited_) then
            MDegradedDT = (1. - SedimentDegradationFactor) * Mass
        ElseIf ((HNSParticleState .EQ. Surface_) .OR.                &
                (HNSParticleState .EQ. WaterColumn_Droplet_) .OR.    &
                (HNSParticleState .EQ. WaterColumn_Dissolved_) .OR.  &
                (HNSParticleState .EQ. WaterColumn_Sedimented_)      &
               )    then
            MDegradedDT = (1. - WaterDegradationFactor) * Mass
        EndIf
                      
        Me%Var%MDegradedDT = MDegradedDT
        
        !the degraded mass must be removed from the particle mass  
              
    end subroutine Degradation   
    
    !----------------------------------------------------------------------------   

    subroutine GetHNSAirTurbHorizontalVel(WindSpeed, Radiation, Latitude, Longitude, CloudCover, DT, HNSAirTurbHorizontalVel)
       !Arguments---------------------------------------------------------------
        real, intent(IN) :: WindSpeed
        real, intent(IN) :: Radiation
        real, intent(IN) :: Latitude
        real, intent(IN) :: Longitude
        real, intent(IN) :: CloudCover
        real, intent(IN) :: DT
        real, intent(OUT) :: HNSAirTurbHorizontalVel

        !Local-------------------------------------------------------------------
        real :: DifCoef
        integer :: StabilityClass
        real :: DifParameter
        logical                                     :: Night
        
        !------------------------------------------------------------------------

        if (Me%Var%AirHorizontalDifCoef .GE. 0.) then
            DifCoef = Me%Var%AirHorizontalDifCoef
        else
            if (Me%Var%StabilityClass .GE. 0.) then
                StabilityClass = Me%Var%StabilityClass
            else
                Night = Night_(Latitude, Longitude)
                !WindSpeed, Night, CloudCover, Radiation
                StabilityClass = StabilityClass_(WindSpeed,Night, CloudCover, Radiation)
            endif
            
            select case(StabilityClass)
            case (Class_A_)
                DifParameter = 0.22 * DT / SQRT(1. + (0.0001 * DT))
            case (Class_B_)
                DifParameter = 0.16 * DT / SQRT(1. + (0.0001 * DT))
            case (Class_C_)
                DifParameter = 0.11 * DT / SQRT(1. + (0.0001 * DT))
            case (Class_D_)
                DifParameter = 0.08 * DT / SQRT(1. + (0.0001 * DT))
            case (Class_E_)
                DifParameter = 0.06 * DT / SQRT(1. + (0.0001 * DT))
            case (Class_F_)
                DifParameter = 0.04 * DT / SQRT(1. + (0.0001 * DT))
            end select
            
            DifCoef = DifParameter*DifParameter/(2*DT)            
        endif

        HNSAirTurbHorizontalVel = SQRT(6.*DifCoef/DT)
        
    end subroutine
    
    !----------------------------------------------------------------------------   

    subroutine GetHNSAirTurbVerticalVel(WindSpeed, Radiation, Latitude, Longitude, CloudCover, DT, HNSAirTurbVerticalVel)
       !Arguments---------------------------------------------------------------
        real, intent(IN) :: WindSpeed
        real, intent(IN) :: Radiation
        real, intent(IN) :: Latitude
        real, intent(IN) :: Longitude
        real, intent(IN) :: CloudCover
        real, intent(IN) :: DT
        real, intent(OUT) :: HNSAirTurbVerticalVel

        !Local-------------------------------------------------------------------
        real :: DifCoef
        integer :: StabilityClass
        real :: DifParameter
        logical                                     :: Night
        
        !------------------------------------------------------------------------
    
        if (Me%Var%AirVerticalDifCoef .GE. 0.) then
            DifCoef = Me%Var%AirVerticalDifCoef
        else
            if (Me%Var%StabilityClass .GE. 0.) then
                StabilityClass = Me%Var%StabilityClass
            else
                Night = Night_(Latitude, Longitude)
                !WindSpeed, Night, CloudCover, Radiation
                StabilityClass = StabilityClass_(WindSpeed,Night, CloudCover, Radiation)
            endif
            
            select case(StabilityClass)
             case (Class_A_)
                DifParameter = 0.2 * DT
            case (Class_B_)
                DifParameter = 0.12 * DT
            case (Class_C_)
                DifParameter = 0.08 * DT / SQRT(1. + (0.0002 * DT))
            case (Class_D_)
                DifParameter = 0.06 * DT / SQRT(1. + (0.0015 * DT))
            case (Class_E_)
                DifParameter = 0.03 * DT / (1. + (0.0003 * DT))
            case (Class_F_)
                DifParameter = 0.016 * DT / (1. + (0.0003 * DT))
            end select
            
            DifCoef = DifParameter*DifParameter/(2*DT)            
        endif

        HNSAirTurbVerticalVel = SQRT(6.*DifCoef/DT)

    end subroutine
    
    !--------------------------------------------------------------------------
    
    logical function Night_(Latitude, Longitude)
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        real                                        :: HourAngle, SunHighAngle
        real                                        :: Declination
        real                                        :: Hour, Minute, Second
        real                                        :: LatitudePI, LongitudePI
        real                                        :: Latitude, Longitude
        integer                                     :: JulDay, STAT_CALL
        real                                        :: GmtReference, RacingWithTheSun

        !GMT Reference
        call GetGmtReference(Me%ObjTime, GmtReference, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeCloudCover - ModuleAtmosphere - ERR00'

        call ExtractDate(Me%ExternalVar%Now, Hour = Hour, Minute = Minute, Second = Second )
        call JulianDay  (Me%ExternalVar%Now, Julday)

        Hour = Hour + Minute/60. + Second/3600.

        !Sun "Racing"
        RacingWithTheSun = RacingWithTheSun_(JulDay)

        !Declination of the Sun
        Declination = SunDeclination_(JulDay)

        HourAngle = HourAngle_(Hour, LongitudePI, GmtReference, RacingWithTheSun)
        LatitudePI  = Latitude * PI / 180.0
        LongitudePI = Longitude

        !use a sunset angle of 10 as safety factor
        SunHighAngle =  Asin( sin(LatitudePI) * sin(Declination)            +   &
                        cos(LatitudePI) * cos(Declination) * cos(HourAngle)) *   &
                        180./PI
    
        if (SunHighAngle .LT. 10.) then
            Night_ = .true.
        else
            Night_ = .false.
        endif

    end function Night_

    !--------------------------------------------------------------------------
    
    integer function StabilityClass_(WindSpeed, Night, CloudCover, Radiation)
    
        !Arguments-------------------------------------------------------------
        real, intent    (IN)                        :: WindSpeed, CloudCover, Radiation
        logical, intent (IN)                        :: Night
        
        !Local-----------------------------------------------------------------
        
        If (CloudCover .GT. 0.95) then
            StabilityClass_ = 4
        elseif (CloudCover .LE. 0.95) then
            if (Night) then
                if(CloudCover .GE. 0.5) then
                    if (WindSpeed .LE. 3) then
                        StabilityClass_ = 5
                    elseif (WindSpeed .GT. 3) then
                        StabilityClass_ = 4
                    endif
                elseif(CloudCover .LT. 0.5) then
                    if (WindSpeed .LE. 3) then
                        StabilityClass_ = 6
                    elseif ((WindSpeed .GT. 3) .AND. (WindSpeed .LE. 5)) then
                        StabilityClass_ = 5
                    elseif (WindSpeed .GT. 5) then
                        StabilityClass_ = 4
                    endif
                endif
            else
               if(Radiation .GT. 851.) then
                    if       (WindSpeed .LE. 2) then
                        StabilityClass_ = 1
                    elseif   ((WindSpeed .GT. 2) .AND. (WindSpeed .LE. 5)) then
                        StabilityClass_ = 2
                    elseif  (WindSpeed .GT. 5) then
                        StabilityClass_ = 3
                    endif
               elseif ((Radiation .GT. 526) .AND. (Radiation .LE. 851)) then
                    if       (WindSpeed .LE. 2) then
                        StabilityClass_ = 1
                    elseif   ((WindSpeed .GT. 2) .AND. (WindSpeed .LE. 3)) then
                        StabilityClass_ = 2
                    elseif  ((WindSpeed .GT. 3) .AND. (WindSpeed .LE. 6)) then
                        StabilityClass_ = 3
                    elseif  (WindSpeed .GT. 6) then
                        StabilityClass_ = 4
                    endif
               elseif ((Radiation .GT. 176) .AND. (Radiation .LE. 526)) then
                    if       (WindSpeed .LE. 2) then
                        StabilityClass_ = 2
                    elseif   ((WindSpeed .GT. 2) .AND. (WindSpeed .LE. 5)) then
                        StabilityClass_ = 3
                    elseif  (WindSpeed .GT. 5) then
                        StabilityClass_ = 4
                    endif
               elseif (Radiation .LE. 176) then
                    StabilityClass_ = 4
               endif
            endif
        endif
        
    end function StabilityClass_
    
    !--------------------------------------------------------------------------
    
    !This function computes the declination of the Sun which is the angle between &
    !the rays of the sun and the plane of the earth's equator. 
    !Possibly taken from Deas, M.L. and Lowney C.L. - "Water Temperature Modelling Review" Sept. 2000.
    real function SunDeclination_(JulDay)

        !Arguments-------------------------------------------------------------
        integer                                     :: JulDay

        !Local-----------------------------------------------------------------

        !Sun Declination
        SunDeclination_ = 23.45 * cos((172.0 - JulDay) * 2.0 * PI / 365.0) * PI / 180.0

    end function SunDeclination_

    !--------------------------------------------------------------------------

    !This function computes radian angles according to the (solar) hour for sun height computation.
    !At 12h, angle is zero (maximum co-sine in sun height computation) and increases towards sunrise & 
    !and sunset (co-sine in sun height computation decreases).
    !Reference? Possibly adapted from Deas, M.L. and Lowney C.L. - "Water Temperature Modelling Review" Sept. 2000.
    real function HourAngle_ (Hour, LongitudePI, GmtReference, RacingWithTheSun)

        !Arguments-------------------------------------------------------------
        real                                        :: LongitudePI, GmtReference
        real                                        :: RacingWithTheSun, Hour

        !Local-----------------------------------------------------------------
        real                                        :: HourC            
        
        !Corrected Hour for longitude, timezone and small disturbances (Racing with the sun)
        ! h   =  h   +    degrees    / deg/h -       h       +       h
        HourC = Hour + (LongitudePI / 15.0) - GmtReference + RacingWithTheSun 
        
        !Hour angle (to_change passar para SCT)
        if (HourC .LT. 12) then
            HourAngle_ = (HourC + 12.0) * PI / 12.0 
        else
            HourAngle_ = (HourC - 12.0) * PI / 12.0 
        end if

    end function HourAngle_

    !--------------------------------------------------------------------------

    !This function accounts for disturbances in the earth’s rotation rate that affect 
    !the time the suns takes to go through the longitude differences. 
    !Taken from "Evapotranspiration Technical Manual" eq. 53 and eq. 54 but reference not found.
    !It adds maximum 10 min or takes maximum 15 min to the hour depending on day of year.
    real function RacingWithTheSun_ (JulDay)

        !Arguments-------------------------------------------------------------
        integer                                     :: JulDay

        !Local-----------------------------------------------------------------
        real                                        :: Aux

        Aux = 2 * PI * (JulDay - 81)/364. 
        RacingWithTheSun_ = 0.1645 * sin(2*Aux) - 0.1255 * cos(Aux) - 0.025 * Aux 

    end function RacingWithTheSun_

    !--------------------------------------------------------------------------

    
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillHNS(HNSID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HNSID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_ 
        integer                                     :: STAT_, nUsers 

        !---------------------------------------------------------------------- 

        STAT_ = UNKNOWN_
        
        call Ready(HNSID, ready_) 

        if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mHNS_,  Me%InstanceID)

            if (nUsers == 0) then
                
                nUsers = DeassociateInstance (mTIME_,           Me%ObjTime          )
                if (nUsers == 0) stop 'KillHNS - ModuleHNS - ERR01'

                !Deallocates Instance
                call DeallocateInstance ()

                HNSID  = 0
                STAT_  = SUCCESS_

            end if

        else 
            STAT_ = ready_
        end if 


        if (present(STAT)) STAT = STAT_


        !------------------------------------------------------------------------

    end subroutine KillHNS

    !--------------------------------------------------------------------------


    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_HNS), pointer          :: AuxObjHNS
        type (T_HNS), pointer          :: PreviousObjHNS

        !Updates pointers
        if (Me%InstanceID == FirstObjHNS%InstanceID) then
            FirstObjHNS    => FirstObjHNS%Next
        else
            PreviousObjHNS => FirstObjHNS
            AuxObjHNS      => FirstObjHNS%Next
            do while (AuxObjHNS%InstanceID /= Me%InstanceID)
                PreviousObjHNS => AuxObjHNS
                AuxObjHNS      => AuxObjHNS%Next
            enddo

            !Now update linked list
            PreviousObjHNS%Next => AuxObjHNS%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------

    subroutine WriteFinalHNS(HNSID, UnitID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HNSID
        integer                                     :: UnitID
        integer, optional                           :: STAT
 
        !Local-------------------------------------------------------------------
        integer                                     :: ready_ 
        integer                                     :: STAT_ 

        !----------------------------------------------------------------------
       
        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            write (UnitID)  Me%Var%Time
            write (UnitID) Me%State%FirstStepIP
            write (UnitID) Me%Var%DTHNSInternalProcesses

            STAT_ = SUCCESS_

        else 
         
            STAT_ = ready_
        end if 


        if (present(STAT))  STAT = STAT_
       

    end subroutine WriteFinalHNS
    
    !--------------------------------------------------------------------------

    subroutine ReadFinalHNS(HNSID, UnitID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HNSID
        integer                                     :: UnitID
        integer, optional                           :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_ 
        integer                                     :: STAT_ 

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HNSID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            read (UnitID) Me%Var%Time
            read (UnitID) Me%State%FirstStepIP
            read (UnitID) Me%Var%DTHNSInternalProcesses

            STAT_ = SUCCESS_

        else 
         
            STAT_ = ready_
        end if 


        if (present(STAT))  STAT = STAT_

    end subroutine ReadFinalHNS

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine Ready (ObjHNS_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjHNS_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjHNS_ID > 0) then
            call LocateObjHNS (ObjHNS_ID)
            ready_ = VerifyReadLock (mHNS_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjHNS (ObjHNSID)

        !Arguments-------------------------------------------------------------
        integer                    :: ObjHNSID

        !Local-----------------------------------------------------------------

        Me => FirstObjHNS
        do while (associated (Me))
            if (Me%InstanceID == ObjHNSID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleHNS - LocateObjHNS - ERR01'

    end subroutine LocateObjHNS

    !--------------------------------------------------------------------------

end module ModuleHNS

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

