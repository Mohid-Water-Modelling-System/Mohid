!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Oil
! PROJECT       : Mohid Water
! MODULE        : ModuleOil
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jun 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module responsbile for computing Oil processe (Spill / dissolution, etc)
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
 
Module ModuleOil

!BOP
!
! !MODULE: ModuleOil

!    !DESCRIPTION: 
!     This model is responsible  for determine spreading velocities of oil Lagrangian Tracers,
!     and for weathering the oil slick
 
! !REVISION HISTORY: 
!    Apr 2001   Ricardo Miranda
!    Sep 2001   Rodrigo Fernandes
!
!
!Input DataFile
!-Lagrangian DataFile:
!       
!    (parameters from the ModuleLagrangian)
!    
!    <BeginOrigin>
!    
!    (parameters from the ModuleLagrangian)
!    
!       <<BeginProperty>>
!    
!       (parameters from the ModuleLagrangian)
!    
!       <<EndProperty>>
!
!       <<BeginOil>>
!       OIL_TIMESERIE         : char                                    [-]         !Name of the Output results file
!       DT_OUTPUT_TIME        : real                                    [-]         !Time between output results
!
!       
!       OIL_SPREADING         : 0/1                                     [1]         !Oil Spreading Process
!       SPREADINGMETHOD       : Fay/ThicknessGradient                   []          !Method for Spreading 
!       USERCOEFVELMANCHA     : real                                    [10.0]      !Empirical Thickness Gradient 
!                                                                                    (typical values 5-30)
!                                                                                    Spreading Vel. Coef. 
!
!       OIL_EVAPORATION       : 0/1                                     [0]         !Oil Evaporation Process
!       EVAPORATIONMETHOD     : EvaporativeExposure/PseudoComponents    
!                               /Fingas                                 []          !Method for Evaporation
!       OIL_DISPERSION        : 0/1                                     [0]         !Oil Dispersion Process
!       DISPERSIONMETHOD      : Delvigne/Mackay                         []          !Method for Dispersion
!
!       OIL_EMULSIFICATION    : 0/1                                     [0]         !Oil Emulsification Process
!       EMULSIFICATIONMETHOD  : Mackay/Rasmussen                        []          !Method for Emulsification
!
!       OIL_DISSOLUTION       : 0/1                                     [0]         !Oil Dissolution Process
!
!       OIL_SEDIMENTATION     : 0/1                                     [0]         !Oil Sedimentation Process
!
!       OILTYPE               : Crude/Refined                           []          !Oil Type
!       API                   : real                                    [-]         !American Petroleum Institute (API) Gravity
!       POURPOINT             : real (ºC)                               [-]         !Pour Point
!       CEMULS                : real (%)                                [0.0]       !Emulsification Constant 
!                                                                                    (% of evaporated oil before emulsification
!                                                                                     begins) 
!       MAXVWATERCONTENT      : real (%)                                [null_real] !Maximum Volume Water Content   
!       ASPHALTENECONTENT     : real (%)                                [-]         !Asphaltene Content
!       WAXCONTENT            : real (%)                                [-]         !Wax Content
!       EMULSPARAMETER        : real                                    [1.6E-6]    !Water Uptake Parameter 
!                                                                                    (typical values 1.0E-6 to 2.0E-6) 
!       TEMPVISCREF           : real (ºC)                               [-]         !Temperature of Reference Viscosity
!       VISCREF               : real (cP)                               [-]         !Reference Dynamic Viscosity
!       VISCCINREF            : real (cSt)                              [-]         !Reference Cinematic Viscosity
!       OWINTERFACIALTENSION  : real (Dyne/cm)                          [-]         !Oil-Water Interfacial Tension
!
!       (the following 3 keywords are only necessary when Evaporation Method = PseudoComponents)
!       NBRDISTCUTS           : int                                     [0]         !Number of Distillation Cuts
!       TDISTEXP              : list(real) (ºC)                         [-]         !Vapour Temperature of Distillate 
!       CPDISTEXP             : list(real) (%)                          [-]         !Cumulative Volume Fraction of Oil Distilled  
!
!       (the following 5 keywords are only necessary when Evaporation Method = Fingas)
!       FINGAS_EVAP_EQTYPE    : Logarithmic / SquareRoot                []          !Evaporation Equation Type
!       FINGAS_EVAP_EMP_DATA  : 0/1                                     [0]         !Knowledge of Empirical Data for Evaporation 
!       FINGAS_EVAP_CONST1    : real                                    [-]         !Empirical Constant 1 
!                                                                                    (Necessary If Fingas_Evap_Emp_Data = 1)
!       FINGAS_EVAP_CONST2    : real                                    [-]         !Empirical Constant 2 
!                                                                                    (Necessary If Fingas_Evap_Emp_Data = 1)
!       PERC_MASSDIST180      : real (%)                                [-]         !%(Wheight) of Oil Evaporated until 180ºC 
!                                                                                    (Necessary If Fingas_Evap_Emp_Data = 0)

!       OIL_CHEM_DISPERSION   : 0/1                                     [0]         !Chemical Dispersants Application
!       P_AREA_SPRAYED        : real (%)                                [-]         !% of Spill Area sprayed whit dispersant
!       EFFICIENCY            : real (%)                                [-]         !% of Area sprayed effectively dispersed
!       START_CHEM_DISPERSION : YYYY MM DD HH MM SS                     [BeginModel]!Starting Time of Dispersant Application
!       END_CHEM_DISPERSION   : YYYY MM DD HH MM SS                     [EndModel]  !Ending Time of Dispersant Application

!       OIL_MEC_CLEANUP       : 0/1                                     [0]         !Mechanical Cleanup Operation
!       START_MEC_CLEANUP     : YYYY MM DD HH MM SS                     [BeginModel]!Starting Time of Mechanical Cleanup Operation
!       End_MEC_CLEANUP       : YYYY MM DD HH MM SS                     [EndModel]  !End Time of Mechanical Cleanup Operation
!       RECOVERY              : real (l/h or l)                         [-]         !rate or volume of Emulsion Recovered
!       RECOVERY_DATAFORM     : Rate / Amount                           []          !DataForm of emulsion recovered
    
!       <<EndOil>>
!    
!    <EndOrigin>
!    
!Properties also needed and present in other Input DataFiles:
!-WaterProperties DataFile
!   
!   Temperature
!   Salinity
!   Cohesive Sediment
!
!-Surface DataFile
!
!   Wind Velocity
!   Atmospheric Pressure
!   Wave Period
!   Wave Height


!Output Results File:
!
!
!    
! to do 
!   -revision of pseudo-components evaporation process (when Nbrdistcuts<5)
!   -shoreline deposition
!   -calculus continuation
!   -better callibration of spreading velocities

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleTimeSerie
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, UngetHorizontalGrid      
    use ModuleGeometry,         only : GetGeometrySize               
    use ModuleMap,              only : GetComputeFaces3D, UngetMap

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartOil
    public  :: ReadFinalOil 
    private :: OilOptions 
    private :: OilOptionsAPI
    private :: EvapPropINI
    private :: F_MaxVWaterContent      !Function 
    private :: F_Cdisp                 !Function
    private :: F_IBP                   !Function
    private :: F_Slope                 !Function
    private :: ReadTimeSerieFile
    private :: AllocateVariables 


    !Selector
    public  :: GetOilDensityOil
    public  :: GetOilAPI
    public  :: GetOilSpreadingVelocity
    public  :: GetOilSpreadingList
    public  :: GetOilSpreading
    public  :: UngetOil


    !Modifier
    public  ::  OilInternalProcesses    !Zero dimensional processes
    private ::      Density
    private ::      Viscosity
    private ::      Evaporation
    private ::      Dispersion
    private ::      Sedimentation
    private ::      Dissolution
    private ::      Emulsification
    private ::      ChemDispersion
    private ::      MecCleanup
    private ::      TimeSerieOutput
    public  ::      OilActiveProcesses      !Influence on movement, growth, etc.
    private ::      OilPropIni
    private ::      AreaTeoric
    private ::      F_CVisc_E               !Function
    private ::      F_ThicknessLimit        !Function
    public  ::      F_FayArea               !Function


    !Destructor
    public  ::  KillOil 
    public  ::      WriteFinalOil

    !Management
    private ::  Ready
    private ::      LocateObjOil

    !Parameter-----------------------------------------------------------------


    !Oilinternalprocesses

    real, parameter :: PETIT               = 1.0E-10
    real, parameter :: Temp15              = 15.0          
        
    !Density

    real, parameter :: CDens_DE            = 0.18       
    real, parameter :: CDens_DT            = 0.0008    
    real, parameter :: FreshWaterDensity15 = 999.0 

    !Viscosity

    real, parameter :: CVisc_V             = 2.5                  
    real, parameter :: CVisc_T             = 5000.0                 
    real, parameter :: CVisc_M             = 0.65                   

    !Evaporation

    real, parameter :: CEvap_A             = 6.3    
    real, parameter :: CEvap_B             = 10.3   
    real, parameter :: CEvap_deltaZ        = 0.97  
 
    !Dispersion

    real, parameter :: CDisp_d0            = 37.5E-6
    real, parameter :: CDisp_Deltad        = 65.0E-6
    real, parameter :: CDisp_Uvi           = 4.0
    real, parameter :: CDisp_b             = 0.032         

    !Sedimentation

    real, parameter :: CSed_d0             = 135.0E-6
    real, parameter :: CSed_Deltad         = 130.0E-6
    real, parameter :: CSed_Sticking       = 1.0E-4     !m3/kg
    !Dissolution

    real, parameter :: SolubilityFreshOil  = 0.030      !kg/m3
    real, parameter :: CDiss_KTransfMass   = 0.01       !m/hr
    real, parameter :: CDiss_decayment     = 0.1           

    !Emulsification

    real, parameter :: CEmuls_1            = 5.0E-7   
    real, parameter :: CEmuls_2            = 1.2E-5   
 

    !OilActiveProcesses
    real, parameter :: CFay_1              = 1.14
    real, parameter :: CFay_2              = 1.45


    !TimeSerie's columns
    
    integer, parameter :: ColMassOil       =  1
    integer, parameter :: ColVolOilBeached =  2
    integer, parameter :: ColVolumeBeached =  3
    integer, parameter :: ColVolumeOil     =  4
    integer, parameter :: ColVolume        =  5
    integer, parameter :: ColArea          =  6
    integer, parameter :: ColAreaTeoric    =  7
    integer, parameter :: ColThickness     =  8
    integer, parameter :: ColMEvaporated   =  9
    integer, parameter :: ColVEvaporated   = 10
    integer, parameter :: ColFMEvaporated  = 11
    integer, parameter :: ColMDispersed    = 12
    integer, parameter :: ColVDispersed    = 13
    integer, parameter :: ColFMDispersed   = 14
    integer, parameter :: ColMSedimented   = 15
    integer, parameter :: ColVSedimented   = 16
    integer, parameter :: ColFMSedimented  = 17
    integer, parameter :: ColMDissolved    = 18
    integer, parameter :: ColVDissolved    = 19
    integer, parameter :: ColFMDissolved   = 20
    integer, parameter :: ColMChemDispersed= 21
    integer, parameter :: ColVChemDispersed= 22
    integer, parameter :: ColFMChemDispersed=23
    integer, parameter :: ColMOilRecovered = 24
    integer, parameter :: ColVOilRecovered = 25
    integer, parameter :: ColFMOilRecovered= 26
    integer, parameter :: ColMWaterContent = 27
    integer, parameter :: ColVWaterContent = 28
    integer, parameter :: ColDensity       = 29
    integer, parameter :: ColViscosity     = 30
    integer, parameter :: ColNbr           = 30


    integer, parameter :: Fay_                = 1
    integer, parameter :: ThicknessGradient_  = 2
    integer, parameter :: Refined             = 1
    integer, parameter :: Crude               = 2
    integer, parameter :: PseudoComponents    = 1
    integer, parameter :: EvaporativeExposure = 2
    integer, parameter :: Fingas              = 3
    integer, parameter :: Logarithmic         = 1
    integer, parameter :: SquareRoot          = 2
    integer, parameter :: Delvigne            = 1
    integer, parameter :: Mackay              = 2
    integer, parameter :: Rasmussen           = 1
    integer, parameter :: Rate                = 1
    integer, parameter :: Amount              = 2


    !Types---------------------------------------------------------------------

    type       T_State
        logical :: TimeSerie  = OFF
        logical :: FirstStepIP = ON
        logical :: FirstStepAP = ON
    end type T_State

    type       T_TimeSerie
        character(LEN = StringLength) :: TimeSerieFile = null_str

        character(LEN = StringLength), dimension(:), pointer :: PropertyList
        real,                          dimension(:), pointer :: DataLine
    end type T_TimeSerie


    type       T_Files
        character(LEN = StringLength) :: ConstructData = null_str
    end type T_Files


    type       T_External
        type(T_Time)                        :: Now
        type(T_Time)                        :: BeginTime
        type(T_Time)                        :: EndTime


        real                                :: Area
        real,    pointer, dimension(:,:  )  :: GridThickness


        !ObjMapp
        integer, pointer, dimension(:,:,:)  :: ComputeFacesU3D
        integer, pointer, dimension(:,:,:)  :: ComputeFacesV3D


        !HorizontalGrid
        real,    pointer, dimension(:,:  )  :: DZX
        real,    pointer, dimension(:,:  )  :: DZY


        !ObjGeometry
        type (T_Size3D)                     :: Size, WorkSize

        !ObjParticProp
        real                                :: Wind
        real                                :: AtmosphericPressure
        real                                :: WaterTemperature
        real                                :: WaterDensity
        real                                :: SPM
        real                                :: WaveHeight
        real                                :: WavePeriod
        real                                :: VolOilBeached        
        real                                :: VolumeBeached        
   end type T_External


    type       T_Var
        
        
        integer                             :: OilType              = null_int
        real                                :: PourPoint            = null_real


        !Time
        real                                :: DTOilInternalProcesses

        
        !Mass
        real                                :: MassOil              = null_real
        real                                :: MassINI              = null_real


        !Volume
        real                                :: Volume               = null_real
        real                                :: VolumeOil            = null_real
        real                                :: VolInic              = null_real


        !Thickness
        real                                :: SlickThickness       = null_real
        real                                :: OilThickness         = null_real


        !Density
        real                                :: Density              = null_real
        real                                :: API                  = null_real
        real                                :: Density15            = null_real
        logical                             :: OilIsFloating        = ON

        !Viscosity
        real                                :: Viscosity            = null_real
        real                                :: CVisc_E              = null_real
        real                                :: ViscRef              = null_real
        real                                :: TempViscRef          = null_real
        real                                :: DensTempViscRef      = null_real
        real                                :: ViscCinRef           = null_real
        real                                :: ViscINI              = null_real
        real                                :: ViscCin              = null_real

        !Spreading
        logical                             :: OilSpreading         = OFF
        integer                             :: SpreadingMethod      = null_int
        real,    pointer, dimension(:,:  )  :: SpreadingVelocityX
        real,    pointer, dimension(:,:  )  :: SpreadingVelocityY
        real                                :: UserCoefVelMancha    = null_real
        real                                :: DiffVelocity         = null_real
        real                                :: ThicknessLimit       = null_real
        real                                :: alfa                 = null_real
        real                                :: alfa_old             = null_real
        real                                :: DiffCoef             = null_real
               
        real                                :: AreaTeoric           = null_real
        real                                :: beta                 = null_real
        real                                :: beta_old             = null_real


        !Evaporation
        logical                             :: OilEvaporation       = OFF
        integer                             :: EvaporationMethod    = null_int
        real                                :: MEvaporatedDT        = null_real
        real                                :: VEvaporatedDT        = null_real
        real                                :: MEvaporated          = null_real
        real                                :: VEvaporated          = null_real
        real                                :: FMEvaporated         = null_real
        real                                :: IBP                  = null_real
        real                                :: Slope                = null_real
        integer                             :: NbrDistCuts          = null_int
        integer                             :: NbrPC                = null_int
        real,pointer,dimension(:)           :: TDistExp          
        real,pointer,dimension(:)           :: CPDistExp         
        real,pointer,dimension(:)           :: VolInicPC          
        real,pointer,dimension(:)           :: VolPC          
        real,pointer,dimension(:)           :: VEvaporatedPCDT          
        real,pointer,dimension(:)           :: PvapPC         
        real,pointer,dimension(:)           :: VmrelPC         
        real,pointer,dimension(:)           :: MWPC         
        logical                             :: Fingas_Evap_Emp_Data = OFF
        integer                             :: Fingas_Evap_EqType   = null_int
        real                                :: Fingas_Evap_Const1   = null_real
        real                                :: Fingas_Evap_Const2   = null_real
        real                                :: Perc_MassDist180     = null_real
        real                                :: Time                 = null_real
        
        !Dispersion
        logical                             :: OilDispersion        = OFF
        integer                             :: DispersionMethod     = null_int
        real                                :: MDispersedDT         = null_real
        real                                :: VDispersedDT         = null_real
        real                                :: MDispersed           = null_real  
        real                                :: VDispersed           = null_real
        real                                :: FMDispersed          = null_real
        real                                :: Cdisp                = null_real
        real                                :: OWInterfacialTension = null_real

        
        !Sedimentation
        logical                             :: OilSedimentation     = OFF
        real                                :: MSedimentedDT        = null_real
        real                                :: VSedimentedDT        = null_real
        real                                :: MSedimented          = null_real  
        real                                :: VSedimented          = null_real
        real                                :: FMSedimented         = null_real

        
        !Emulsification
        logical                             :: OilEmulsification    = OFF
        integer                             :: EmulsificationMethod = null_int
        real                                :: MWaterContentDT      = null_real
        real                                :: VWaterContentDT      = null_real
        real                                :: MWaterContent        = null_real
        real                                :: VWaterContent        = null_real
        real                                :: Cemuls               = null_real
        real                                :: MaxVWaterContent     = null_real
        real                                :: AsphalteneContent    = null_real
        real                                :: EmulsParameter       = null_real
        real                                :: WaxContent           = null_real

        
        !Dissolution
        logical                             :: OilDissolution       = OFF
        real                                :: MDissolvedDT         = null_real
        real                                :: VDissolvedDT         = null_real
        real                                :: MDissolved           = null_real  
        real                                :: VDissolved           = null_real
        real                                :: FMDissolved          = null_real
        real                                :: SolubilityOilInWater = null_real    

        !Removal:

        ! - Chemical Dispersion
        logical                             :: OilChemDispersion    = OFF
        real                                :: MChemDispersedDT     = null_real
        real                                :: VChemDispersedDT     = null_real
        real                                :: MChemDispersed       = null_real  
        real                                :: VChemDispersed       = null_real
        real                                :: FMChemDispersed      = null_real
        real                                :: P_AreaSprayed        = null_real
        real                                :: Efficiency           = null_real
        real                                :: MassOilStartChemDisp = null_real
        type (T_Time)                       :: Start_ChemDispersion
        type (T_Time)                       :: End_ChemDispersion
        
        ! - Mechanical Cleanup
        logical                             :: OilMecCleanup        = OFF
        real                                :: VEmulsionRecoveryRate= null_real
        type (T_Time)                       :: Start_Mec_Cleanup
        type (T_Time)                       :: End_Mec_Cleanup
        real                                :: Recovery             = null_real
        integer                             :: Recovery_DataForm    = null_int
        real                                :: MOilRecoveredDT      = null_real
        real                                :: VOilRecoveredDT      = null_real
        real                                :: MOilRecovered        = null_real
        real                                :: VOilRecovered        = null_real
        real                                :: FMOilRecovered       = null_real

     end type T_Var


    type      T_Oil
        integer                 :: InstanceID
        type(T_State    )       :: State
        type(T_Files    )       :: Files
        type(T_TimeSerie)       :: TimeSerie
        type(T_Var      )       :: Var
        type(T_External )       :: ExternalVar

        type(T_Time)            :: NextInternalComputeTime
        type(T_Time)            :: NextActiveComputeTime

        !Instance of ModuleTime
        integer                 :: ObjTime = 0

        !Instance of ModuleTimeSerie
        integer                 :: ObjTimeSerie = 0
         
        !Instance of ModuleHorizontalGrid
        integer                 :: ObjHorizontalGrid = 0
       
        !Instance of ModuleGeometry
        integer                 :: ObjGeometry = 0

        !Instance of ModuleMap
        integer                 :: ObjMap = 0
 
        !Instance of Module_EnterData
        integer                 :: ObjEnterData = 0

        type (T_Oil), pointer   :: Next

    end type T_Oil

    !Global Module Variables
    type (T_Oil), pointer              :: FirstObjOil
    type (T_Oil), pointer              :: Me


    !----------------------------------------------------------------------------

    contains



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine StartOil (OilID,                                                 &
                         TimeID,                                                &
                         EnterDataID,                                           &
                         HorizontalGridID,                                      &
                         GeometryID,                                            &
                         MapID,                                                 &
                         DT,                                                    &
                         ContCalc,                                              &
                         PropertyListIN,                                        &
                         ExtractType,                                           &
                         STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: OilID
        integer                                     :: TimeID
        integer                                     :: EnterDataID
        integer                                     :: HorizontalGridID
        integer                                     :: MapID
        integer                                     :: GeometryID
        integer, optional, intent(IN )              :: ExtractType    
        integer, optional, intent(OUT)              :: STAT 
        real,              intent(IN )              :: DT
        character(*), optional, dimension(:), pointer :: PropertyListIN
        
        logical, intent(IN)  :: ContCalc 

        !External----------------------------------------------------------------

        integer :: STAT_CALL
        integer :: ready_ 
        integer :: FromFile
        integer :: ExtractType_  
        integer :: nUsers
        !Local-------------------------------------------------------------------

        integer :: STAT_

        real :: MecCleanupTime
        
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mOil_)) then
            nullify (FirstObjOil)
            call RegisterModule (mOil_) 
        endif
        
        call Ready(OilID, ready_) 

        if (ready_ .EQ. OFF_ERR_) then

            !Allocates Instance
            call AllocateInstance

            nullify (Me%TimeSerie%PropertyList  )
            nullify (Me%TimeSerie%DataLine      )


            !Associates External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjEnterData      = AssociateInstance (mENTERDATA_,      EnterDataID     )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjMap            = AssociateInstance (mMAP_,            MapID           )


            !Gets Time
            call GetComputeTimeLimits(Me%ObjTime,                                        &
                                      Me%ExternalVar%EndTime,                            &
                                      Me%ExternalVar%BeginTime,                          &
                                      STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartOil - ModuleOil - ERR01'

            ! Actualized the time
            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'StartOil - ModuleOil - ERR02'


            call GetGeometrySize(Me%ObjGeometry,                            &
                                 Size           = Me%ExternalVar%Size,      &
                                 WorkSize       = Me%ExternalVar%WorkSize,  &
                                 STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartOil - ModuleOil - ERR03'


cd2 :       if (present(ExtractType)) then
                ExtractType_ = ExtractType
            else
                call GetExtractType(FromFile = FromFile)
                ExtractType_ = FromFile
            end if cd2
                

            call OilOptions       (DT, ContCalc, ExtractType_)

cd3 :       if (present(PropertyListIN)) then
                call ReadTimeSerieFile(ExtractType_, PropertyListIN)
            else cd3
                call ReadTimeSerieFile(ExtractType_)
            end if cd3


            call AllocateVariables


            !Initialization of integrated values
            
ifContCalc: if (.NOT. ContCalc ) then
                Me%Var%Time              = 0.0
                Me%Var%MEvaporated       = 0.0
                Me%Var%VEvaporated       = 0.0
                Me%Var%FMEvaporated      = 0.0
                Me%Var%MDispersed        = 0.0
                Me%Var%VDispersed        = 0.0
                Me%Var%FMDispersed       = 0.0
                Me%Var%MSedimented       = 0.0
                Me%Var%VSedimented       = 0.0
                Me%Var%FMSedimented      = 0.0
                Me%Var%MDissolved        = 0.0
                Me%Var%VDissolved        = 0.0
                Me%Var%FMDissolved       = 0.0
                Me%Var%MChemDispersed    = 0.0  
                Me%Var%VChemDispersed    = 0.0
                Me%Var%FMChemDispersed   = 0.0
                Me%Var%MOilRecovered     = 0.0  
                Me%Var%VOilRecovered     = 0.0
                Me%Var%FMOilRecovered    = 0.0
                Me%Var%MWaterContent     = 0.0
                Me%Var%VWaterContent     = 0.0

                Me%State%FirstStepIP     = ON
                Me%State%FirstStepAP     = ON

                !initialization of some properties
                if (Me%Var%OilEmulsification) Me%Var%MaxVWaterContent = F_MaxVWaterContent()
            
                Me%Var%SolubilityOilInWater = SolubilityFreshOil

                Me%Var%ThicknessLimit       = F_ThicknessLimit ()
                    
cd5:            if (Me%Var%OilEvaporation) then
cd6:                if (Me%Var%EvaporationMethod .EQ. PseudoComponents) then
                        
cd7:                    if (Me%Var%OilType .EQ. Refined) then
                            write (*, *) 'Refined Oil + PseudoComponents Method : '
                            write (*, *) 'Rough aproximation in evaporation process'
                            write (*, *) '(Refined products usually don`t demonstrate linear distilattion curves)'
                        end if  cd7
                                            
                   else if (Me%Var%EvaporationMethod .EQ. EvaporativeExposure) then cd6

cd8:                    if (Me%Var%OilType .EQ. Refined) then
                            write (*, *) 'Refined Oil + Evaporative Exposure Method : '
                            write (*, *) 'Rough aproximation in evaporation process'
                            write (*, *) '(Refined products usually don`t demonstrate linear distilattion curves)'
                        end if cd8
                    
                       Me%Var%IBP              = F_IBP  ()
                       Me%Var%Slope            = F_Slope()
                   end if cd6

                end if cd5


cd9:            if (Me%Var%OilMecCleanup) then

cd10:               if (Me%Var%Recovery_DataForm .EQ. Rate) then
                
                        Me%Var%VEmulsionRecoveryRate  = Me%Var%Recovery / (1000.0 * 3600.0)
            
                    else if (Me%Var%Recovery_DataForm .EQ. Amount)  then cd10
                
                        MecCleanupTime                = Me%Var%End_Mec_Cleanup - Me%Var%Start_Mec_Cleanup
                        Me%Var%VEmulsionRecoveryRate  = Me%Var%Recovery / (1000.0 * MecCleanupTime)
            
                    end if cd10

                end if cd9

            end if ifContCalc

            nUsers = DeassociateInstance (mENTERDATA_,      Me%ObjEnterData     )
            if (nUsers == 0) stop 'StartOil - ModuleOil - ERR98'

            OilID = Me%InstanceID

            STAT_ = SUCCESS_
        
        else 
            
            stop 'StartOil - ModuleOil - ERR99' 

        end if 


        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine StartOil

    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance
                                                    
        !Local-----------------------------------------------------------------
        type (T_Oil), pointer                  :: NewObjOil
        type (T_Oil), pointer                  :: PreviousObjOil


        !Allocates new instance
        allocate (NewObjOil)
        nullify  (NewObjOil%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjOil)) then
            FirstObjOil         => NewObjOil
            Me                  => NewObjOil
        else
            PreviousObjOil      => FirstObjOil
            Me                  => FirstObjOil%Next
            do while (associated(Me))
                PreviousObjOil  => Me
                Me              => Me%Next
            enddo
            Me                  => NewObjOil
            PreviousObjOil%Next => NewObjOil
        endif

        Me%InstanceID = RegisterNewInstance (mOIL_)


    end subroutine AllocateInstance
    
    !----------------------------------------------------------------------------

    subroutine OilOptions(DT, ContCalc, ExtractType) 

        !Arguments---------------------------------------------------------------
        real   , intent(IN)     :: DT
        logical, intent(IN)     :: ContCalc

        !External----------------------------------------------------------------
        integer                 :: flag
        integer                 :: ExtractType
        integer                 :: STAT_CALL
        character(LEN = StringLength) :: String


        !Internal------------------------------------------------------------------------

        character(LEN = StringLength), parameter :: Char_Fay                   = trim(adjustl('Fay'))
        character(LEN = StringLength), parameter :: Char_ThicknessGradient     = trim(adjustl('ThicknessGradient'))
        character(LEN = StringLength), parameter :: Char_Refined               = trim(adjustl('Refined'))
        character(LEN = StringLength), parameter :: Char_Crude                 = trim(adjustl('Crude'))
        character(LEN = StringLength), parameter :: Char_PseudoComponents      = trim(adjustl('PseudoComponents'))
        character(LEN = StringLength), parameter :: Char_EvaporativeExposure   = trim(adjustl('EvaporativeExposure'))
        character(LEN = StringLength), parameter :: Char_Fingas                = trim(adjustl('Fingas'))
        character(LEN = StringLength), parameter :: Char_Logarithmic           = trim(adjustl('Logarithmic'))
        character(LEN = StringLength), parameter :: Char_SquareRoot            = trim(adjustl('SquareRoot'))
        character(LEN = StringLength), parameter :: Char_Delvigne              = trim(adjustl('Delvigne'))
        character(LEN = StringLength), parameter :: Char_Mackay                = trim(adjustl('Mackay'))
        character(LEN = StringLength), parameter :: Char_Rasmussen             = trim(adjustl('Rasmussen'))
        character(LEN = StringLength), parameter :: Char_Rate                  = trim(adjustl('Rate'))
        character(LEN = StringLength), parameter :: Char_Amount                = trim(adjustl('Amount'))

        !--------------------------------------------------------------------------------

        Me%NextInternalComputeTime = Me%ExternalVar%Now
        Me%NextActiveComputeTime   = Me%ExternalVar%Now

        if (.NOT. ContCalc ) then

            call GetData(Me%Var%DTOilInternalProcesses,                                  &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType = ExtractType,                                       &    
                         keyword    = 'DT_OIL_INTPROCESSES',                             &         
                         ClientModule ='ModuleOil',                                      &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OilOptions - ModuleOil - ERR01'


cd1 :       if (flag .EQ. 0) then

                Me%Var%DTOilInternalProcesses = DT

            else cd1

cd2 :           if ((mod(Me%Var%DTOilInternalProcesses, DT) .EQ. 0.0      ) .AND.       &
                    (    Me%Var%DTOilInternalProcesses      .GE. DT)) then

                    !OK

                else cd2

                    stop 'OilOptions - ModuleOil - ERR02'

                end if cd2   

            end if cd1                  
                                 

            call GetData(String,                                                         &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = ExtractType,                                     &
                         keyword      ='OILTYPE',                                        &
                         ClientModule ='ModuleOil',                                      &
                         Default      = Char_Crude,                                      &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'OilOptions - ModuleOil - ERR04'

case1 :     select case(trim(adjustl(String)))
            case(Char_Refined)

                Me%Var%OilType = Refined

            case(Char_Crude)

                Me%Var%OilType = Crude

            case default

                if (STAT_CALL /= SUCCESS_) stop 'OilOptions - ModuleOil - ERR05'

            end select case1


            call GetData(Me%Var%OilEvaporation,                                          &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType     = ExtractType,                                   & 
                         keyword        = 'OIL_EVAPORATION',                             &
                         Default        = OFF,                                           &
                         ClientModule   ='ModuleOil',                                    &
                         STAT           = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'OilOptions - ModuleOil - ERR06'


ifevap: if (Me%Var%OilEvaporation) then

            call GetData(String,                                                           &
                         Me%ObjEnterData,                                               &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        &
                         keyword      ='EVAPORATIONMETHOD',                                 &
                         ClientModule ='ModuleOil',                                         &
                         Default      = Char_EvaporativeExposure,                           &
                         STAT         = STAT_CALL)        
            if (STAT_CALL .NE. SUCCESS_)                                                    &
                call SetError(FATAL_, KEYWORD_,                                               &
                             "Subroutine OilOptions; Module ModuleOil. ERR07") 

case2 :         select case(trim(adjustl(String)))
                case(Char_PseudoComponents)

                    Me%Var%EvaporationMethod = PseudoComponents

                case(Char_EvaporativeExposure)

                    Me%Var%EvaporationMethod = EvaporativeExposure

                case(Char_Fingas)

                    Me%Var%EvaporationMethod = Fingas

                case default

                    call SetError(FATAL_, KEYWORD_,                                       &
                                 "Subroutine OilOptions; Module ModuleOil. ERR08") 
            end select case2


cd5:        if (Me%Var%EvaporationMethod  .EQ. PseudoComponents) then

                call GetData(Me%Var%NbrDistCuts,                                        &
                             Me%ObjEnterData,                                           &
                             flag,                                                          &
                             SearchType   = ExtractType,                                    & 
                             keyword      = 'NBRDISTCUTS',                                  &
                             Default      = 0,                                              & 
                             ClientModule ='ModuleOil',                                     &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    call SetError(FATAL_, KEYWORD_,                                           &
                                   "Subroutine OilOptions; Module ModuleOil. ERR26") 



               allocate(Me%Var%TDistExp(1:Me%Var%NbrDistCuts), STAT = STAT_CALL)
               if (STAT_CALL .NE. SUCCESS_)                                                 &
                   call SetError(FATAL_, INTERNAL_,                                           &
                                 "Subroutine OilOptions; Module ModuleOil. ERR27") 


               allocate(Me%Var%CPDistExp(1:Me%Var%NbrDistCuts), STAT = STAT_CALL)
               if (STAT_CALL .NE. SUCCESS_)                                                 &
                   call SetError(FATAL_, INTERNAL_,                                           &
                                 "Subroutine OilOptions; Module ModuleOil. ERR28") 


               call GetData(Me%Var%CPDistExp,                                           &
                            Me%ObjEnterData,                                            &
                            flag         = flag,                                            &
                            SearchType   = ExtractType,                                     &
                            keyword      = 'CPDISTEXP',                                     &
                            ClientModule ='ModuleOil',                                      &
                            STAT         = STAT_CALL)
cd6 :          if   (STAT_CALL .EQ. SIZE_ERR_)  then
                     write(*,*) 
                     write(*,*) 'Error calling GetData.  '
                     write(*,*) 'Number of distillation cuts is incorrect:'
                     write(*,*) '    NbrDistCuts  =', Me%Var%NbrDistCuts
                     write(*,*) '   CPDistExpData =', flag
                     stop       'Subroutine OilOptions; Module ModuleOil. ERR29.'           

               else if ((STAT_CALL .NE. SIZE_ERR_) .AND.                                     &
                        (STAT_CALL .NE. SUCCESS_)) then cd6                                                                    
                    stop 'Subroutine OilOptions; Module ModuleOil. ERR30.'           
               end if cd6           

               call GetData(Me%Var%TDistExp,                                            &
                                    Me%ObjEnterData,                                    &
                                    flag,                                                   &
                                    SearchType   = ExtractType,                             &
                                    keyword      = 'TDISTEXP',                              &
                                    ClientModule ='ModuleOil',                              &
                                    STAT         = STAT_CALL)
cd7 :           if   (STAT_CALL .EQ. SIZE_ERR_)  then
                     write(*,*) 
                     write(*,*) 'Error calling GetData.  '
                     write(*,*) 'Number of distillation cuts is incorrect:'
                     write(*,*) '    NbrDistCuts  =', Me%Var%NbrDistCuts
                     write(*,*) '    TDistExpData =', flag
                     stop       'Subroutine OilOptions; Module ModuleOil. ERR31.'           

                else if ((STAT_CALL .NE. SIZE_ERR_) .AND.                                    &
                                 (STAT_CALL .NE. SUCCESS_)) then cd7            
                    stop 'Subroutine OilOptions; Module ModuleOil. ERR32.'           
                end if cd7           

            else if (Me%Var%EvaporationMethod  .EQ. Fingas) then cd5
                
                call GetData(Me%Var%Fingas_Evap_Emp_Data,                               &
                             Me%ObjEnterData,                                           &
                             flag,                                                          &
                             SearchType   = ExtractType,                                    & 
                             keyword      = 'FINGAS_EVAP_EMP_DATA',                         &
                             Default      = OFF,                                            &
                             ClientModule ='ModuleOil',                                     &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    call SetError(FATAL_, KEYWORD_,                                           &
                                   "Subroutine OilOptions; Module ModuleOil. ERR26") 
              
                
cd76:           if (Me%Var%Fingas_Evap_Emp_Data) then

                    call GetData(Me%Var%Fingas_Evap_Const1,                             &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                      &
                                 SearchType   = ExtractType,                                & 
                                 keyword      = 'FINGAS_EVAP_CONST1',                       &
                                 ClientModule ='ModuleOil',                                 &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                              &
                        call SetError(FATAL_, KEYWORD_,                                       &
                                       "Subroutine OilOptions; Module ModuleOil. ERR26") 

                    call GetData(Me%Var%Fingas_Evap_Const2,                             &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                      &
                                 SearchType   = ExtractType,                                & 
                                 keyword      = 'FINGAS_EVAP_CONST2',                       &
                                 ClientModule ='ModuleOil',                                 &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                              &
                        call SetError(FATAL_, KEYWORD_,                                       &
                                       "Subroutine OilOptions; Module ModuleOil. ERR26") 

                else if (.NOT. Me%Var%Fingas_Evap_Emp_Data) then    cd76

                    call GetData(Me%Var%Perc_MassDist180,                               &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                      &
                                 SearchType   = ExtractType,                                & 
                                 keyword      = 'PERC_MASSDIST180',                         &
                                 ClientModule ='ModuleOil',                                 &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                              &
                        call SetError(FATAL_, KEYWORD_,                                       &
                                       "Subroutine OilOptions; Module ModuleOil. ERR26") 
                end if  cd76


                call GetData(String,                                                        &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             SearchType   = ExtractType,                                    &
                             keyword      ='FINGAS_EVAP_EQTYPE',                            &
                             ClientModule ='ModuleOil',                                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL .NE. SUCCESS_)                                                &
                    call SetError(FATAL_, KEYWORD_,                                         &
                                 "Subroutine OilOptions; Module ModuleOil. ERR07") 

cd77 :          if (flag .EQ. 1) then
case77 :            select case(trim(adjustl(String)))
                        case(Char_Logarithmic)

                            Me%Var%Fingas_Evap_EqType = Logarithmic

                        case(Char_SquareRoot)

                            Me%Var%Fingas_Evap_EqType = SquareRoot

                        case default

                            call SetError(FATAL_, KEYWORD_,                                   &
                                         "Subroutine OilOptions; Module ModuleOil. ERR08") 
                    end select case77

                else cd77

                    call SetError(FATAL_, KEYWORD_,                                           &
                                 "Subroutine OilOptions; Module ModuleOil. ERR09") 
                end if cd77


            end if cd5
                
        end if ifevap

        call GetData(Me%Var%OilDispersion,                                               &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType     = ExtractType,                                       &    
                     keyword        = 'OIL_DISPERSION',                                  &
                     Default        = OFF,                                               &
                     ClientModule   = 'ModuleOil',                                       &
                     STAT           = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                                        &
            call SetError(FATAL_, INTERNAL_, "Subroutine OilOptions; Module ModuleOil. ERR35")

ifdisp: if  (Me%Var%OilDispersion) then

            call GetData(String,                                                            &
                         Me%ObjEnterData,                                                   &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        &
                         keyword      ='DISPERSIONMETHOD',                                  &
                         ClientModule ='ModuleOil',                                         &
                         Default      = Char_Delvigne,                                      &
                         STAT         = STAT_CALL)        
            if (STAT_CALL .NE. SUCCESS_)                                                    &
                call SetError(FATAL_, KEYWORD_,                                               &
                             "Subroutine OilOptions; Module ModuleOil. ERR10") 

case3 :         select case(trim(adjustl(String)))
                case(Char_Delvigne)

                    Me%Var%DispersionMethod = Delvigne

                case(Char_Mackay)

                    Me%Var%DispersionMethod = Mackay

                case default

                    call SetError(FATAL_, KEYWORD_,                                       &
                                 "Subroutine OilOptions; Module ModuleOil. ERR11") 
            end select case3


cd9:        if (Me%Var%DispersionMethod .EQ. Mackay) then

                call GetData(Me%Var%OWInterfacialTension,                               &
                             Me%ObjEnterData,                                           &
                             flag,                                                          &
                             SearchType   = ExtractType,                                    & 
                             keyword      = 'OWINTERFACIALTENSION',                         &
                             ClientModule ='ModuleOil',                                     &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    call SetError(FATAL_, KEYWORD_,                                           &
                                 "Subroutine OilOptions; Module ModuleOil. ERR25") 

            end if cd9

        end if ifdisp






        call OilOptionsAPI(Me%ObjEnterData,                                             &
                           ExtractType  = ExtractType,                                      &
                           API          = Me%Var%API)




        call GetData(Me%Var%PourPoint,                                                  &
                     Me%ObjEnterData,                                                   &
                     flag,                                                                  &
                     SearchType   = ExtractType,                                            & 
                     keyword      = 'POURPOINT',                                            &     
                     ClientModule ='ModuleOil',                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                          &
            call SetError(FATAL_, KEYWORD_,                                                   &
                         "Subroutine OilOptions; Module ModuleOil. ERR17") 

        
        call GetData(Me%Var%OilEmulsification,                                           &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType     = ExtractType,                                          & 
                     keyword        = 'OIL_EMULSIFICATION',                                 &
                     Default        = OFF,                                                  &
                     ClientModule   ='ModuleOil',                                           &
                     STAT           = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                                        &
            call SetError(FATAL_, INTERNAL_, "Subroutine OilOptions; Module ModuleOil. ERR37")

ifemul: if (Me%Var%OilEmulsification) then

            call GetData(String,                                                           &
                         Me%ObjEnterData,                                               &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        &
                         keyword      ='EMULSIFICATIONMETHOD',                              &
                         Default      = Char_Rasmussen,                                     &
                         ClientModule ='ModuleOil',                                         &
                         STAT         = STAT_CALL)        
            if (STAT_CALL .NE. SUCCESS_)                                                    &
                call SetError(FATAL_, KEYWORD_,                                               &
                             "Subroutine OilOptions; Module ModuleOil. ERR13") 

case4 :         select case(trim(adjustl(String)))
                case(Char_Rasmussen)

                    Me%Var%EmulsificationMethod = Rasmussen

                case(Char_Mackay)

                    Me%Var%EmulsificationMethod = Mackay

                case default

                    call SetError(FATAL_, KEYWORD_,                                       &
                                 "Subroutine OilOptions; Module ModuleOil. ERR14") 
            end select case4

        
            call GetData(Me%Var%Cemuls,                                                 &
                         Me%ObjEnterData,                                               &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        & 
                         keyword      = 'CEMULS',                                           &
                         Default      = 0.0,                                                &
                         ClientModule ='ModuleOil',                                         &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                call SetError(FATAL_, KEYWORD_,                                               &
                             "Subroutine OilOptions; Module ModuleOil. ERR18") 


            call GetData(Me%Var%MaxVWaterContent,                                       &
                         Me%ObjEnterData,                                               &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        & 
                         keyword      = 'MAXVWATERCONTENT',                                 &
                         Default      = null_real,                                          &
                         ClientModule ='ModuleOil',                                         &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                call SetError(FATAL_, KEYWORD_,                                               &
                             "Subroutine OilOptions; Module ModuleOil. ERR19") 


cd11:        if (Me%Var%EmulsificationMethod .EQ. Rasmussen) then

                call GetData(Me%Var%AsphalteneContent,                                  &
                             Me%ObjEnterData,                                           &
                             flag,                                                          &
                             SearchType   = ExtractType,                                    & 
                             keyword      = 'ASPHALTENECONTENT',                            &
                             ClientModule ='ModuleOil',                                     &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    call SetError(FATAL_, KEYWORD_,                                           &
                                 "Subroutine OilOptions; Module ModuleOil. ERR20") 


                call GetData(Me%Var%WaxContent,                                         &
                             Me%ObjEnterData,                                           &
                             flag,                                                          &
                             SearchType   = ExtractType,                                    &  
                             keyword      = 'WAXCONTENT',                                   &
                             ClientModule ='ModuleOil',                                     &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    call SetError(FATAL_, KEYWORD_,                                           &
                                 "Subroutine OilOptions; Module ModuleOil. ERR21") 
        
            else if (Me%Var%EmulsificationMethod .EQ. Mackay) then cd11
            
                call GetData(Me%Var%EmulsParameter,                                     &
                             Me%ObjEnterData,                                           &
                             flag,                                                          &
                             SearchType   = ExtractType,                                    &  
                             keyword      = 'EmulsParameter',                               &
                             Default      = 1.6E-6,                                         &
                             ClientModule ='ModuleOil',                                     &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    call SetError(FATAL_, KEYWORD_,                                           &
                                 "Subroutine OilOptions; Module ModuleOil. ERR21") 

            end if cd11

        end if ifemul

        call GetData(Me%Var%TempViscRef,                                                &
                     Me%ObjEnterData,                                                   &
                     flag,                                                                  &
                     SearchType   = ExtractType,                                            &  
                     keyword      = 'TEMPVISCREF',                                          &
                     ClientModule ='ModuleOil',                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                          &
            call SetError(FATAL_, KEYWORD_,                                                   &
                         "Subroutine OilOptions; Module ModuleOil. ERR22") 


        call GetData(Me%Var%ViscRef,                                                    &
                     Me%ObjEnterData,                                                   &
                     flag,                                                                  &
                     SearchType   = ExtractType,                                            & 
                     keyword      = 'VISCREF',                                              &
                     ClientModule ='ModuleOil',                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                          &
            call SetError(FATAL_, KEYWORD_,                                                   &
                         "Subroutine OilOptions; Module ModuleOil. ERR23") 


        call GetData(Me%Var%ViscCinRef,                                                 &
                     Me%ObjEnterData,                                                   &
                     flag,                                                                  &
                     SearchType   = ExtractType,                                            &    
                     keyword      = 'VISCCINREF',                                           &
                     ClientModule ='ModuleOil',                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                          &
            call SetError(FATAL_, KEYWORD_,                                                   &
                         "Subroutine OilOptions; Module ModuleOil. ERR24") 






        call GetData(Me%Var%OilSpreading,                                                &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType     = ExtractType,                                       &    
                     keyword        = 'OIL_SPREADING',                                   &
                     Default        = ON,                                                &
                     ClientModule   ='ModuleOil',                                        &
                     STAT           = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                                        &
            call SetError(FATAL_, INTERNAL_, "Subroutine OilOptions; Module ModuleOil. ERR33")


ifspr:  if (Me%Var%OilSpreading) then

            call GetData(String,                                                            &
                         Me%ObjEnterData,                                                   &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        &
                         keyword      ='SPREADINGMETHOD',                                   &
                         Default      = Char_ThicknessGradient,                             &
                         ClientModule ='ModuleOil',                                         &
                         STAT         = STAT_CALL)        
            if (STAT_CALL .NE. SUCCESS_)                                                    &
                call SetError(FATAL_, KEYWORD_,                                               &
                             "Subroutine OilOptions; Module ModuleOil. ERR39") 

case5 :         select case(trim(adjustl(String)))
                case(Char_Fay)

                    Me%Var%SpreadingMethod = Fay_

                case(Char_ThicknessGradient)

                    Me%Var%SpreadingMethod = ThicknessGradient_

                case default

                    call SetError(FATAL_, KEYWORD_,                                       &
                                 "Subroutine OilOptions; Module ModuleOil. ERR40") 
            end select case5



cd13:       if (Me%Var%SpreadingMethod  .EQ. ThicknessGradient_) then

                call GetData(Me%Var%UserCoefVelMancha,                                  &
                             Me%ObjEnterData,                                           &
                             flag,                                                          &
                             SearchType   = ExtractType,                                    & 
                             keyword      = 'USERCOEFVELMANCHA',                            &
                             Default      = 10.0,                                           &
                             ClientModule = 'ModuleOil',                                    &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    call SetError(FATAL_, KEYWORD_,                                           &
                                   "Subroutine OilOptions; Module ModuleOil. ERR42") 

            end if cd13

        end if ifspr

        call GetData(Me%Var%OilSedimentation,                                               &
                     Me%ObjEnterData,                                                       &
                     flag,                                                                  &
                     SearchType     = ExtractType,                                          & 
                     keyword        = 'OIL_SEDIMENTATION',                                  &
                     Default        = OFF,                                                  &
                     ClientModule   = 'ModuleOil',                                          &
                     STAT           = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                                        &
            call SetError(FATAL_, INTERNAL_, "Subroutine OilOptions; Module ModuleOil. ERR36")


        call GetData(Me%Var%OilDissolution,                                                 &
                     Me%ObjEnterData,                                                       &
                     flag,                                                                  &
                     SearchType     = ExtractType,                                          & 
                     keyword        = 'OIL_DISSOLUTION',                                    &
                     Default        = OFF,                                                  &
                     ClientModule   = 'ModuleOil',                                          &
                     STAT           = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                                        &
            call SetError(FATAL_, INTERNAL_, "Subroutine OilOptions; Module ModuleOil. ERR38")



    !Removal


        call GetData(Me%Var%OilChemDispersion,                                              &
                     Me%ObjEnterData,                                                       &
                     flag,                                                                  &
                     SearchType     = ExtractType,                                          & 
                     keyword        = 'OIL_CHEM_DISPERSION',                                &
                     Default        = OFF,                                                  &
                     ClientModule   = 'ModuleOil',                                          &
                     STAT           = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                                        &
            call SetError(FATAL_, INTERNAL_, "Subroutine OilOptions; Module ModuleOil. ERR43")


ifcdis: if (Me%Var%OilChemDispersion) then

            call GetData(Me%Var%Start_ChemDispersion,                                   &
                         Me%ObjEnterData,                                               &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        & 
                         keyword      = 'START_CHEM_DISPERSION',                            &
                         Default      = Me%ExternalVar%BeginTime,                       &
                         ClientModule = 'ModuleOil',                                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                call SetError(FATAL_, KEYWORD_,                                               &
                               "Subroutine OilOptions; Module ModuleOil. ERR43a") 


            call GetData(Me%Var%End_ChemDispersion,                                     &
                         Me%ObjEnterData,                                               &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        & 
                         keyword      = 'END_CHEM_DISPERSION',                              &
                         Default      = Me%ExternalVar%EndTime,                         &
                         ClientModule = 'ModuleOil',                                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                call SetError(FATAL_, KEYWORD_,                                               &
                               "Subroutine OilOptions; Module ModuleOil. ERR43b") 


            call GetData(Me%Var%P_AreaSprayed,                                          &
                         Me%ObjEnterData,                                               &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        & 
                         keyword      = 'P_AREA_SPRAYED',                                   &
                         ClientModule = 'ModuleOil',                                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                call SetError(FATAL_, KEYWORD_,                                               &
                               "Subroutine OilOptions; Module ModuleOil. ERR44") 

            call GetData(Me%Var%Efficiency,                                             &
                         Me%ObjEnterData,                                               &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        & 
                         keyword      = 'EFFICIENCY',                                       &
                         ClientModule = 'ModuleOil',                                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                call SetError(FATAL_, KEYWORD_,                                               &
                               "Subroutine OilOptions; Module ModuleOil. ERR45") 

        end if ifcdis

        call GetData(Me%Var%OilMecCleanup,                                                  &
                     Me%ObjEnterData,                                                       &
                     flag,                                                                  &
                     SearchType     = ExtractType,                                          &    
                     keyword        = 'OIL_MEC_CLEANUP',                                    &
                     Default        = OFF,                                                  &
                     ClientModule   = 'ModuleOil',                                          &
                     STAT           = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                                        &
            call SetError(FATAL_, INTERNAL_, "Subroutine OilOptions; Module ModuleOil. ERR46")


ifmcle: if (Me%Var%OilMecCleanup) then
        
            call GetData(Me%Var%Start_Mec_Cleanup,                                      &
                         Me%ObjEnterData,                                               &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        & 
                         keyword      = 'START_MEC_CLEANUP',                                &
                         Default      = Me%ExternalVar%BeginTime,                       &
                         ClientModule = 'ModuleOil',                                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                call SetError(FATAL_, KEYWORD_,                                               &
                               "Subroutine OilOptions; Module ModuleOil. ERR46a") 


            call GetData(Me%Var%End_Mec_Cleanup,                                        &
                         Me%ObjEnterData,                                               &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        & 
                         keyword      = 'END_MEC_CLEANUP',                                  &
                         Default      = Me%ExternalVar%EndTime,                         &
                         ClientModule ='ModuleOil',                                         &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                call SetError(FATAL_, KEYWORD_,                                               &
                               "Subroutine OilOptions; Module ModuleOil. ERR46b") 


            call GetData(Me%Var%Recovery,                                               &
                         Me%ObjEnterData,                                               &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        & 
                         keyword      = 'RECOVERY',                                         &
                         ClientModule ='ModuleOil',                                         &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                call SetError(FATAL_, KEYWORD_,                                               &
                               "Subroutine OilOptions; Module ModuleOil. ERR47") 

            call GetData(String,                                                            &
                         Me%ObjEnterData,                                                   &
                         flag,                                                              &
                         SearchType   = ExtractType,                                        &
                         keyword      ='RECOVERY_DATAFORM',                                 &
                         ClientModule ='ModuleOil',                                         &
                         STAT         = STAT_CALL)        
            if (STAT_CALL .NE. SUCCESS_)                                                    &
                call SetError(FATAL_, KEYWORD_,                                               &
                             "Subroutine OilOptions; Module ModuleOil. ERR48") 

cd14 :       if (flag .EQ. 1) then

case6 :         select case(trim(adjustl(String)))
                    case(Char_Rate)

                        Me%Var%Recovery_DataForm = Rate

                    case(Char_Amount)

                        Me%Var%Recovery_DataForm = Amount

                    case default

                        call SetError(FATAL_, KEYWORD_,                                       &
                                     "Subroutine OilOptions; Module ModuleOil. ERR49") 
                end select case6

            else cd14

                call SetError(FATAL_, KEYWORD_,                                               &
                             "Subroutine OilOptions; Module ModuleOil. ERR50") 
            end if cd14

        end if ifmcle

        end if 
    !------------------------------------------------------------------------

    end subroutine OilOptions 

    !----------------------------------------------------------------------------

    subroutine OilOptionsAPI(ObjEnterData, ExtractType, API)

        !Arguments---------------------------------------------------------------
        integer                                     :: ObjEnterData
        integer, intent(IN )                        :: ExtractType     
        real,    intent(OUT)                        :: API     

        !External----------------------------------------------------------------

        integer :: flag, STAT_CALL

        call GetData(API,                                                                   &
                     ObjEnterData,                                                          &
                     flag,                                                                  &
                     SearchType   = ExtractType,                                            & 
                     keyword      = 'API',                                                  &
                     ClientModule ='ModuleOil',                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                          &
            call SetError(FATAL_, KEYWORD_,                                                 &
                         "Subroutine OilOptionsAPI; Module ModuleOil. ERR01") 

        !------------------------------------------------------------------------

        end subroutine OilOptionsAPI

    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------

    subroutine ReadTimeSerieFile(ExtractType, PropertyListIN)

        !Arguments---------------------------------------------------------------
        integer, intent(IN )                        :: ExtractType     
        character(*), optional, dimension(:), pointer :: PropertyListIN

        !External----------------------------------------------------------------
        integer                                     :: flag, STAT_CALL
        character(LEN = StringLength)               :: TimeSerieDataFile = null_str


        !Local-------------------------------------------------------------------
        integer :: aux, Prop

        !------------------------------------------------------------------------

        call GetData(Me%TimeSerie%TimeSerieFile,                                        &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = ExtractType,                                            &
                     keyword      ='OIL_TIMESERIE',                                         &
                     ClientModule ='ModuleOil',                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                          &
            call SetError(FATAL_, KEYWORD_,                                                   &
                         "Subroutine ReadTimeSerieFile; Module ModuleOil. ERR01") 


cd1 :   if (flag == 1) then
            Me%State%TimeSerie = ON

        else cd1
         
            Me%State%TimeSerie = OFF
        end if cd1


        !Start TimeSerie Input
cd2 :   if (Me%State%TimeSerie) then

cd3 :       if (present(PropertyListIN)) then
                aux = ColNbr + size(PropertyListIN)

            else

                aux = ColNbr
            end if cd3


            allocate(Me%TimeSerie%PropertyList(aux))
            
            Me%TimeSerie%PropertyList(ColMassOil       ) = trim(adjustl('MassOil'      ))
            Me%TimeSerie%PropertyList(ColVolOilBeached ) = trim(adjustl('VolOilBeached'))
            Me%TimeSerie%PropertyList(ColVolumeBeached ) = trim(adjustl('VolumeBeached'))
            Me%TimeSerie%PropertyList(ColVolumeOil     ) = trim(adjustl('VolumeOil'    ))
            Me%TimeSerie%PropertyList(ColVolume        ) = trim(adjustl('Volume'       ))
            Me%TimeSerie%PropertyList(ColArea          ) = trim(adjustl('Area'         ))
            Me%TimeSerie%PropertyList(ColAreaTeoric    ) = trim(adjustl('TeoricalArea' ))
            Me%TimeSerie%PropertyList(ColThickness     ) = trim(adjustl('Thickness'    ))
            Me%TimeSerie%PropertyList(ColMEvaporated   ) = trim(adjustl('MEvaporated'  ))
            Me%TimeSerie%PropertyList(ColVEvaporated   ) = trim(adjustl('VEvaporated'  ))
            Me%TimeSerie%PropertyList(ColFMEvaporated  ) = trim(adjustl('FMEvaporated' ))
            Me%TimeSerie%PropertyList(ColMDispersed    ) = trim(adjustl('MDispersed'   ))
            Me%TimeSerie%PropertyList(ColVDispersed    ) = trim(adjustl('VDispersed'   ))
            Me%TimeSerie%PropertyList(ColFMDispersed   ) = trim(adjustl('FMDispersed'  ))
            Me%TimeSerie%PropertyList(ColMSedimented   ) = trim(adjustl('MSedimented'  ))
            Me%TimeSerie%PropertyList(ColVSedimented   ) = trim(adjustl('VSedimented'  ))
            Me%TimeSerie%PropertyList(ColFMSedimented  ) = trim(adjustl('FMSedimented' ))
            Me%TimeSerie%PropertyList(ColMDissolved    ) = trim(adjustl('MDissolved'   ))
            Me%TimeSerie%PropertyList(ColVDissolved    ) = trim(adjustl('VDissolved'   ))
            Me%TimeSerie%PropertyList(ColFMDissolved   ) = trim(adjustl('FMDissolved'  ))
            Me%TimeSerie%PropertyList(ColMChemDispersed) = trim(adjustl('MChemDisp'    ))
            Me%TimeSerie%PropertyList(ColVChemDispersed) = trim(adjustl('VChemDisp'    ))
            Me%TimeSerie%PropertyList(ColFMChemDispersed)= trim(adjustl('FMChemDisp'   ))
            Me%TimeSerie%PropertyList(ColMOilRecovered) = trim(adjustl('MOilRecovered' ))
            Me%TimeSerie%PropertyList(ColVOilRecovered) = trim(adjustl('VOilRecovered' ))
            Me%TimeSerie%PropertyList(ColFMOilRecovered)= trim(adjustl('FMOilRecovered'))
            Me%TimeSerie%PropertyList(ColMWaterContent ) = trim(adjustl('MWaterContent'))
            Me%TimeSerie%PropertyList(ColVWaterContent ) = trim(adjustl('VWaterContent'))
            Me%TimeSerie%PropertyList(ColDensity       ) = trim(adjustl('Density'      ))
            Me%TimeSerie%PropertyList(ColViscosity     ) = trim(adjustl('Viscosity'    ))

cd4 :       if (present(PropertyListIN)) then
do1 :       do Prop = (ColNbr+1), aux
                Me%TimeSerie%PropertyList(Prop) = trim(adjustl(PropertyListIN(aux-ColNbr)))
            end do do1
            end if cd4


            allocate(Me%TimeSerie%DataLine(aux),                                       &
                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                call SetError(FATAL_, KEYWORD_,                                               &
                             "Subroutine ReadTimeSerieFile; Module ModuleOil. ERR03") 
            
            Me%TimeSerie%DataLine(:) = null_real



            call GetFileName(Me%ObjEnterData,                                           &
                             TimeSerieDataFile,                                             &
                             STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                call SetError(FATAL_, KEYWORD_,                                               &
                             "Subroutine ReadTimeSerieFile; Module ModuleOil. ERR04") 

            call StartTimeSerie(Me%ObjTimeSerie,                                        &
                                Me%ObjTime,                                             &
                                TimeSerieDataFile  = TimeSerieDataFile,                 &
                                PropertyList       = Me%TimeSerie%PropertyList,         &
                                Extension          ='sro',                              &
                                ResultFileName     = Me%TimeSerie%TimeSerieFile,        &
                                STAT               = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                call SetError(FATAL_, KEYWORD_,                                         &
                             "Subroutine ReadTimeSerieFile; Module ModuleOil. ERR05") 
        end if cd2

        !------------------------------------------------------------------------

    end subroutine ReadTimeSerieFile 

    !----------------------------------------------------------------------------



    !----------------------------------------------------------------------------

    Subroutine AllocateVariables

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------

        integer :: ILB, IUB 
        integer :: JLB, JUB

        !------------------------------------------------------------------------

        ILB = Me%ExternalVar%Size%ILB
        IUB = Me%ExternalVar%Size%IUB

        JLB = Me%ExternalVar%Size%JLB
        JUB = Me%ExternalVar%Size%JUB

        allocate(Me%Var%SpreadingVelocityX(ILB:IUB, JLB:JUB))
        allocate(Me%Var%SpreadingVelocityY(ILB:IUB, JLB:JUB))

        Me%Var%SpreadingVelocityX(:,:) = null_real
        Me%Var%SpreadingVelocityY(:,:) = null_real

        !------------------------------------------------------------------------

    End Subroutine AllocateVariables




    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetOilDensityOil(OilID, OilDensity, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: OilID
        real,              intent(OUT)              :: OilDensity
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(OilID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            OilDensity = Me%Var%Density
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetOilDensityOil


    !--------------------------------------------------------------------------

    subroutine GetOilAPI(OilID, API, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: OilID
        real,              intent(OUT)              :: API
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(OilID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            API = Me%Var%API
            STAT_ = SUCCESS_
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetOilAPI

    !--------------------------------------------------------------------------

    subroutine GetOilSpreadingVelocity(OilID,                                              &
                                       SpreadingVelocityX,                                 &
                                       SpreadingVelocityY,                                 &
                                       DiffVelocity,                                       &
                                       STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: OilID
        real,    optional, pointer, dimension(:,:)  :: SpreadingVelocityX
        real,    optional, pointer, dimension(:,:)  :: SpreadingVelocityY
        real,    optional                           :: DiffVelocity
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(OilID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then

cd2 :       if      (Me%Var%SpreadingMethod .EQ. ThicknessGradient_) then

                call Read_Lock(mOil_, Me%InstanceID)

                SpreadingVelocityX => Me%Var%SpreadingVelocityX

                call Read_Lock(mOil_, Me%InstanceID)

                SpreadingVelocityY => Me%Var%SpreadingVelocityY

            else if (Me%Var%SpreadingMethod .EQ. Fay_              ) then cd2

                DiffVelocity = Me%Var%DiffVelocity

            end if cd2


            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetOilSpreadingVelocity

    !----------------------------------------------------------------------------

    subroutine GetOilSpreadingList(ThicknessGradient, Fay)

        !Arguments---------------------------------------------------------------
        integer, optional, intent(OUT) :: ThicknessGradient 
        integer, optional, intent(OUT) :: Fay 

        !------------------------------------------------------------------------

        if (present(ThicknessGradient)) ThicknessGradient = ThicknessGradient_
        if (present(Fay              )) Fay               = Fay_

        !------------------------------------------------------------------------

    end subroutine GetOilSpreadingList

    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------

    subroutine GetOilSpreading(OilID, SpreadingMethod,  STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: OilID
        integer,           intent(OUT)              :: SpreadingMethod 
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(OilID, ready_) 

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            SpreadingMethod = Me%Var%SpreadingMethod

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine GetOilSpreading

    !----------------------------------------------------------------------------

    subroutine UngetOil(OilID, Array, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: OilID    
        real, pointer, dimension(:,:)               :: Array
        integer, optional, intent (OUT)             :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(OilID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_UnLock(mOIL_, Me%InstanceID, "UngetOil")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine UngetOil

    !----------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    !Zero dimensional processes
    subroutine OilInternalProcesses(OilID,                                               &
                                    Wind,                                                &
                                    AtmosphericPressure,                                 &
                                    WaterTemperature,                                    &
                                    WaterDensity,                                        &
                                    SPM,                                                 &
                                    VWaterContent,                                       & 
                                    MWaterContent,                                       & 
                                    MDispersed,                                          &
                                    OilDensity,                                          &
                                    OilViscosity,                                        &
                                    FMDispersed,                                         &
                                    FMEvaporated,                                        & 
                                    VolTotOilBeached,                                    &
                                    VolTotBeached,                                       &
                                    VolumeTotalIN,                                       &
                                    VolumeTotalOUT,                                      &
                                    AreaTotal,                                           &
                                    AreaTotalOUT,                                        &
                                    DataLineIN,                                          &
                                    WaveHeight,                                          &
                                    WavePeriod,                                          &    
                                    STAT)
        
        !Arguments---------------------------------------------------------------
        integer                                     :: OilID
        real,              intent(IN )              :: Wind
        real,              intent(IN )              :: AtmosphericPressure
        real,              intent(IN )              :: WaterTemperature
        real,              intent(IN )              :: WaterDensity
        real,              intent(IN )              :: SPM
        real,              intent(IN )              :: WaveHeight
        real,              intent(IN )              :: WavePeriod
        real,              intent(IN )              :: VolTotOilBeached
        real,              intent(IN )              :: VolTotBeached
        real,              intent(IN )              :: VolumeTotalIN
        real,              intent(OUT)              :: VolumeTotalOUT
        real,              intent(OUT)              :: VWaterContent
        real,              intent(OUT)              :: MWaterContent
        real,              intent(OUT)              :: MDispersed
        real,              intent(OUT)              :: OilDensity
        real,              intent(OUT)              :: OilViscosity      
        real,              intent(OUT)              :: FMDispersed      
        real,              intent(OUT)              :: FMEvaporated      
        real,              intent(OUT)              :: AreaTotalOUT
        real,              intent(IN )              :: AreaTotal
        real,    optional, dimension(:), pointer    :: DataLineIN
        integer, optional, intent(OUT)              :: STAT


        !Local-------------------------------------------------------------------
        integer                                     :: ready_ 
        integer                                     :: STAT_CALL 
        integer                                     :: STAT_ 

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(OilID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            ! Actualized the time
            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                                              &
                call SetError(FATAL_, INTERNAL_,                                                      &
                             "Subroutine OilInternalProcesses; Module ModuleOil. ERR01") 

cd2 :       if (Me%ExternalVar%Now .GE. Me%NextInternalComputeTime) then
                
                VolumeTotalOUT                  = null_real
                Me%Var%VolumeOil                = VolumeTotalIN
                

                if (AreaTotal >= 0.) then

                    Me%ExternalVar%Area         = AreaTotal

                else

                    Me%ExternalVar%Area         = Me%Var%AreaTeoric

                endif

                Me%ExternalVar%VolOilBeached        = VolTotOilBeached
                Me%ExternalVar%VolumeBeached        = VolTotBeached
                Me%ExternalVar%Wind                 = Wind
                Me%ExternalVar%AtmosphericPressure  = AtmosphericPressure
                Me%ExternalVar%WaterTemperature     = WaterTemperature
                Me%ExternalVar%WaterDensity         = WaterDensity
                Me%ExternalVar%SPM                  = SPM
                Me%ExternalVar%WaveHeight           = WaveHeight
                Me%ExternalVar%WavePeriod           = WavePeriod


                if ( (Me%State%FirstStepIP) .AND. (Me%Var%OilEvaporation)           &
                     .AND. (Me%Var%EvaporationMethod .EQ. PseudoComponents) )           &
                    call EvapPropINI()
                                            
       

                Me%Var%SlickThickness = Me%Var%Volume    / (max(AllmostZero,Me%ExternalVar%Area))
                Me%Var%OilThickness   = Me%Var%VolumeOil / (max(AllmostZero,Me%ExternalVar%Area))

                if ( ((Me%Var%OilDispersion) .AND. (Me%Var%DispersionMethod .EQ. Delvigne)) &
                    .OR. (Me%Var%OilSedimentation) )                                            &
                    Me%Var%Cdisp      = F_Cdisp ()

                !Integrated Values
                if ((Me%Var%OilEmulsification) .and. (Me%Var%OilIsFloating))                                                   &
                    call Emulsification

                if ((Me%Var%OilEvaporation) .and. (Me%Var%OilIsFloating))                       &
                    call Evaporation

                if ((Me%Var%OilDispersion)  .and. (Me%Var%OilIsFloating))                       &
                    call Dispersion 

                if (Me%Var%OilSedimentation)                                                    &
                    call Sedimentation 

                if (Me%Var%OilDissolution)                                                      &
                    call Dissolution 

                if (Me%Var%OilChemDispersion)                                                   &
                    call ChemDispersion 

                if (Me%Var%OilMecCleanup)                                                       &
                    call MecCleanup                            

                call Density   

                !when oil is denser than surrounding water, oil will sink 
                if (.not. Me%Var%OilIsFloating) then                 
                    if (Me%Var%Density .GT. Me%ExternalVar%WaterDensity)                        &
                        write(*,*)  
                        write(*,*) 'Oil is denser than surrounding water, and will sink'
                        write(*,*) 'Evaporation,dispersion and emulsification processes will be interrupted'
                        write(*,*)  
                        Me%Var%OilIsFloating  = OFF
                end if

                call Viscosity 




                !Output into an ASCII file 
cd3 :           if (Me%State%TimeSerie) then
cd4 :           if (present(DataLineIN  )) then
                    call TimeSerieOutput(DataLineIN)

                else cd4

                    call TimeSerieOutput
                end if cd4
                end if cd3


                Me%NextInternalComputeTime = Me%ExternalVar%Now + Me%Var%DTOilInternalProcesses


                VolumeTotalOUT                 = Me%Var%VolumeOil 
                VWaterContent                  = Me%Var%VWaterContent
                MWaterContent                  = Me%Var%MWaterContent

                !AreaTotalOUT, OilDensity e MDispersed necessaria para o calculo da concentracao no ModuleLagrangian
                AreaTotalOUT                   = Me%ExternalVar%Area
                OilDensity                     = Me%Var%Density
                MDispersed                     = Me%Var%MDispersed
                
                OilViscosity                    = Me%Var%Viscosity        
                FMEvaporated                    = Me%Var%FMEvaporated        
                FMDispersed                     = Me%Var%FMDispersed    
            end if cd2


            Me%State%FirstStepIP = OFF


            call null_time   (Me%ExternalVar%Now)

            STAT_ = SUCCESS_

        else cd1
         
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                                          &
            STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine OilInternalProcesses

    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------

    subroutine Density

        !Arguments---------------------------------------------------------------

        !------------------------------------------------------------------------
        
        Me%Var%Density  = Me%Var%VWaterContent * Me%ExternalVar%WaterDensity            &
                              + (1- Me%Var%VWaterContent) * Me%Var%Density15                &
                              * (1 - CDens_DT*(Me%ExternalVar%WaterTemperature - Temp15))       &  
                              * (1 + CDens_DE*Me%Var%FMEvaporated) 


        !------------------------------------------------------------------------

    end subroutine Density

    !----------------------------------------------------------------------------

    subroutine Viscosity

        !Arguments---------------------------------------------------------------

        !------------------------------------------------------------------------
               
                       
        Me%Var%Viscosity  = 1000.*(Me%Var%ViscRef/1000.) * exp( (Me%Var%CVisc_E         &
                                * Me%Var%FMEvaporated) + (CVisc_V * Me%Var%VWaterContent)   &
                                 / (1. - CVisc_M * Me%Var%VWaterContent)                        &
                                + CVisc_T*(1./(Me%ExternalVar%WaterTemperature+273.15)          &
                                - 1./(Me%Var%TempViscRef+273.15))  )
        !a viscosidade anterior é calculada em cP
                  
        Me%Var%ViscCin    = 1.0e6 * (Me%Var%Viscosity/1000.) / Me%Var%Density !cSt
        !unidades - cSt
                                          

        !------------------------------------------------------------------------

    end subroutine Viscosity

    !----------------------------------------------------------------------------



    !----------------------------------------------------------------------------

    subroutine Evaporation

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer :: n
        real :: KTransfMass
        real :: AverageMW
        real :: InternalConstant
        REAL,DIMENSION(1:Me%Var%NbrPC):: MolesPC
        REAL,DIMENSION(1:Me%Var%NbrPC):: MolarFractionPC
        REAL,DIMENSION(1:Me%Var%NbrPC):: VolumeFractionPC
        !------------------------------------------------------------------------


cd1:    if (Me%Var%EvaporationMethod  .EQ. PseudoComponents) then
        
do1:        do n=1,Me%Var%NbrPC 
                MolesPC(n)      = Me%Var%VolPC(n) / Me%Var%VmrelPC(n)
            end do  do1
        
do2:        do n=1,Me%Var%NbrPC 
cd2:            if (Me%Var%VolPC(n) .EQ. 0) then
                    MolarFractionPC(n)            = 0.0
                    VolumeFractionPC(n)           = 0.0
                            
                else    cd2
                    MolarFractionPC(n)            = MolesPC(n) / Sum(MolesPC)
                    VolumeFractionPC(n)           = Me%Var%VolPC(n) / Sum(Me%Var%VolPC)
                    AverageMW                     = AverageMW +MolarFractionPC(n)*Me%Var%MWPC(n)

                end if  cd2

            end do  do2

            KTransfMass =  0.0048*Me%ExternalVar%Wind **(7./9.)                                 &    
                           * ( (1.3676)*(sqrt(0.018/AverageMW))**(2./3.) )                          &
                           * (sqrt(4.*Me%ExternalVar%Area/Pi))**(-1./9.)

do22:       do n=1,Me%Var%NbrPC 

                Me%Var%VEvaporatedPCDT(n)     =(1-Me%Var%VWaterContent)*KTransfMass         &
                                                   * Me%Var%VolumeOil*Me%Var%PvapPC(n)      &
                                                   * Me%Var%VmrelPC(n)*VolumeFractionPC(n)      &
                                                   /(R*Me%Var%SlickThickness*(Me%ExternalVar%WaterTemperature+273.15))

            end do  do22


    
            Me%Var%VEvaporatedDT          = Sum(Me%Var%VEvaporatedPCDT)
            Me%Var%MEvaporatedDT          = Me%Var%VEvaporatedDT * Me%Var%Density

do3:        do n=1,Me%Var%NbrPC 
cd3:           if (Me%Var%VolPC(n)-(Me%Var%VEvaporatedPCDT(n))                              &
                   * Me%Var%DTOilInternalProcesses .GT. 0) then            

                   Me%Var%VolPC(n)        = Me%Var%VolPC(n)-(Me%Var%VEvaporatedPCDT(n)) &
                                                * Me%Var%DTOilInternalProcesses
               else cd3
                   Me%Var%VolPC(n)        = 0.0
               end if   cd3

            end do  do3


        else if (Me%Var%EvaporationMethod  .EQ. EvaporativeExposure) then  cd1
            
            KTransfMass                       = 1.5e-3 * Me%ExternalVar%Wind **0.78

            Me%Var%VEvaporatedDT          = (1-Me%Var%VWaterContent)                        &
                                                *(KTransfMass*Me%ExternalVar%Area)* exp(CEvap_A &
                                                -(CEvap_B/(Me%ExternalVar%WaterTemperature      &
                                                +273.15))*(Me%Var%IBP + Me%Var%Slope        &
                                                *(Me%Var%VEvaporated/ Me%Var%VolInic)))
            
            Me%Var%MEvaporatedDT          = Me%Var%VEvaporatedDT * Me%Var%Density



        else if (Me%Var%EvaporationMethod  .EQ. Fingas) then  cd1

            if (Me%Var%Fingas_Evap_EqType .EQ.Logarithmic) then
                
                Me%Var%Time = Me%Var%Time + Me%Var%DTOilInternalProcesses

                if (Me%Var%Fingas_Evap_Emp_Data) then
                
                    InternalConstant         = (Me%Var%MassINI/100.0)                           &         
                                               * (Me%Var%Fingas_Evap_Const1                     &
                                               + Me%Var%Fingas_Evap_Const2                      &
                                               * Me%ExternalVar%WaterTemperature)

                elseif (.NOT. Me%Var%Fingas_Evap_Emp_Data) then
                    
                    InternalConstant         = (Me%Var%MassINI/100.0)                           &
                                               * (0.165 * Me%Var%Perc_MassDist180               &
                                               + 0.045 * (Me%ExternalVar%WaterTemperature - 15))
                    
                end if
                    
                if (Me%Var%Time .LE. 60 * exp(1.0))   then
                    
                    Me%Var%MEvaporatedDT =  InternalConstant / (exp(1.0) * 60.0)
                
                else

                    Me%Var%MEvaporatedDT =  InternalConstant                                    &
                                                * exp(-Me%Var%MEvaporated/InternalConstant)/60.0
                
                end if   
                
                Me%Var%VEvaporatedDT     = Me%Var%MEvaporatedDT / Me%Var%Density

            
            else if (Me%Var%Fingas_Evap_EqType .EQ.SquareRoot) then

                if (Me%Var%Fingas_Evap_Emp_Data) then
                
                    InternalConstant         = (Me%Var%MassINI/100.0)                           &
                                                * (Me%Var%Fingas_Evap_Const1                    &
                                                + Me%Var%Fingas_Evap_Const2                     &
                                                * Me%ExternalVar%WaterTemperature)

                elseif (.NOT. Me%Var%Fingas_Evap_Emp_Data) then
                    
                    InternalConstant         = (Me%Var%MassINI/100.0)                           &
                                               * (0.0254 * Me%Var%Perc_MassDist180              &
                                               + 0.01 * (Me%ExternalVar%WaterTemperature - 15))
                    
                end if

                if (Me%State%FirstStepIP)   then
                    
                    Me%Var%MEvaporatedDT = InternalConstant                                     &
                                               * sqrt(Me%Var%DTOilInternalProcesses/60.)        &
                                               / Me%Var%DTOilInternalProcesses
                
                else

                    Me%Var%MEvaporatedDT = InternalConstant * InternalConstant                  &
                                               / (2*Me%Var%MEvaporated*60.)

                end if   
                
                Me%Var%VEvaporatedDT     = Me%Var%MEvaporatedDT / Me%Var%Density


            end if 

        end if cd1

    
cd4:   if  (Me%Var%MassOil - (Me%Var%MEvaporatedDT) * Me%Var%DTOilInternalProcesses      &
           .GT. 0.0)  then  
            
            Me%Var%MassOil            = Me%Var%MassOil - (Me%Var%MEvaporatedDT)          &
                                            * Me%Var%DTOilInternalProcesses
            
            Me%Var%VolumeOil          = max(0.0,Me%Var%VolumeOil - (Me%Var%VEvaporatedDT)        &
                                            * Me%Var%DTOilInternalProcesses)
        
        else    cd4
            Me%Var%MEvaporatedDT      = Me%Var%MassOil / Me%Var%DTOilInternalProcesses
            Me%Var%VEvaporatedDT      = Me%Var%MEvaporatedDT / Me%Var%Density
            Me%Var%MassOil            = 0.0
            Me%Var%VolumeOil          = 0.0
        end if  cd4
 
        Me%Var%Volume         = PETIT + Me%Var%VolumeOil /  (1 - Me%Var%VWaterContent)

        Me%Var%SlickThickness = Me%Var%Volume / (max(AllmostZero,Me%ExternalVar%Area))
        Me%Var%OilThickness   = Me%Var%VolumeOil / (max(AllmostZero,Me%ExternalVar%Area))

        Me%Var%MEvaporated    = Me%Var%MEvaporated + Me%Var%MEvaporatedDT                &
                                    * Me%Var%DTOilInternalProcesses
        Me%Var%VEvaporated    = Me%Var%VEvaporated + Me%Var%VEvaporatedDT                &
                                    *  Me%Var%DTOilInternalProcesses
        Me%Var%FMEvaporated   = Me%Var%MEvaporated/Me%Var%MassINI

        !------------------------------------------------------------------------

    end subroutine Evaporation

    !----------------------------------------------------------------------------



    !----------------------------------------------------------------------------

    subroutine Dispersion

        !Arguments---------------------------------------------------------------


        !Internal----------------------------------------------------------------
        real            :: Hrms
        real            :: Dba
        real            :: Fwc
        integer         :: n

        !------------------------------------------------------------------------




cd1:    if (Me%Var%DispersionMethod .EQ. Delvigne)  then

            Hrms                              = Me%ExternalVar%WaveHeight / sqrt(2.0)
            Dba                               = 0.0034 * Me%ExternalVar%WaterDensity            &
                                                * Gravity * Hrms**2

cd2:        if (CDisp_Uvi .GT. Me%ExternalVar%Wind) then

                Fwc = 0

            else    cd2

                Fwc = CDisp_b*(Me%ExternalVar%Wind - CDisp_Uvi) / Me%ExternalVar%WavePeriod

            end if  cd2

            Me%Var%MDispersedDT           = (1-Me%Var%MWaterContent) * (Me%Var%Cdisp    &
                                                * Dba**0.57 * Fwc * CDisp_d0**0.7 * CDisp_Deltad)   &
                                                * Me%ExternalVar%Area 
        
        else if (Me%Var%DispersionMethod .EQ. Mackay)  then   cd1
   
            Me%Var%MDispersedDT           = (0.11/3600) * Me%Var%MassOil                    &
                                                * ( Me%ExternalVar%Wind + 1 )**2                &
                                                / ( 1 + 50 * sqrt(Me%Var%Viscosity) *           &
                                                (100.0 * Me%Var%SlickThickness) &
                                                * Me%Var%OWInterfacialTension)
       
        end if cd1
        
        Me%Var%VDispersedDT               = Me%Var%MDispersedDT / Me%Var%Density
          

ifPC:   if ((Me%Var%OilEvaporation) .AND. (Me%Var%EvaporationMethod  .EQ. PseudoComponents)) then

do1:        do n=1,Me%Var%NbrPC 
       
cd3:           if  (Me%Var%VolPC(n)-(Me%Var%VDispersedDT/Me%Var%NbrPC)                      &
                                                    * Me%Var%DTOilInternalProcesses .GT. 0) then 

                    Me%Var%VolPC(n)           = Me%Var%VolPC(n)-(Me%Var%VDispersedDT        &
                                                    /Me%Var%NbrPC) * Me%Var%DTOilInternalProcesses

               else cd3

                    Me%Var%VolPC(n)           = 0.0

               end if   cd3

            end do  do1

        end if ifPC

cd4:    if (Me%Var%MassOil - (Me%Var%MDispersedDT) * Me%Var%DTOilInternalProcesses      &
            .GT. 0) then

            Me%Var%MassOil                = Me%Var%MassOil - (Me%Var%MdispersedDT)      &
                                                * Me%Var%DTOilInternalProcesses
            Me%Var%VolumeOil              = max(0.0,Me%Var%VolumeOil - (Me%Var%VDispersedDT)    &
                                                * Me%Var%DTOilInternalProcesses)
        
        else    cd4

            Me%Var%MDispersedDT           = Me%Var%MassOil / Me%Var%DTOilInternalProcesses
            Me%Var%VDispersedDT           = Me%Var%MDispersedDT / Me%Var%Density
            Me%Var%MassOil                = 0.0
            Me%Var%VolumeOil              = 0.0

        end if  cd4

        Me%Var%Volume                     = PETIT + Me%Var%VolumeOil                        &
                                                /  (1 - Me%Var%VWaterContent)

        Me%Var%SlickThickness             = Me%Var%Volume / (max(AllmostZero,Me%ExternalVar%Area))
        Me%Var%OilThickness               = Me%Var%VolumeOil / (max(AllmostZero,Me%ExternalVar%Area))


        Me%Var%MDispersed                 = Me%Var%MDispersed + Me%Var%MDispersedDT     &
                                                * Me%Var%DTOilInternalProcesses
        Me%Var%VDispersed                 = Me%Var%VDispersed  + Me%Var%VDispersedDT    &
                                                * Me%Var%DTOilInternalProcesses
        Me%Var%FMDispersed                = Me%Var%MDispersed/Me%Var%MassINI

        !------------------------------------------------------------------------

    end subroutine Dispersion
        !------------------------------------------------------------------------


        !------------------------------------------------------------------------

    subroutine Sedimentation

        !Arguments---------------------------------------------------------------


        !Internal----------------------------------------------------------------
        real            :: Hrms
        real            :: Dba                      !J/m2
        real            :: Fwc
        real            :: IntrusionDepth
        real            :: DissipationRate          !J/(m3.s)
        real            :: OilConcentrationDT       !Kg/(m3.s)
        real            :: OilConcentration         !Kg/m3
        integer         :: n
        !------------------------------------------------------------------------

        Hrms                              = Me%ExternalVar%WaveHeight / sqrt(2.0)
        Dba                               = 0.0034 * Me%ExternalVar%WaterDensity                &
                                            * Gravity * Hrms**2
        
        
cd1:    if (CDisp_Uvi .GT. Me%ExternalVar%Wind) then
            Fwc                           = 0
        else
            Fwc                           = CDisp_b*(Me%ExternalVar%Wind - CDisp_Uvi)           &
                                            / Me%ExternalVar%WavePeriod
        end if  cd1

        IntrusionDepth                    = 1.5 * Me%ExternalVar%WaveHeight


        OilConcentrationDT                = (1-Me%Var%MWaterContent)* (Me%Var%Cdisp         &
                                            * Dba**0.57 * Fwc * CSed_d0**0.7  * CSed_Deltad)        &
                                            / IntrusionDepth 

        OilConcentration                  = OilConcentrationDT * Me%Var%DTOilInternalProcesses

        DissipationRate                   = Dba / (IntrusionDepth * Me%ExternalVar%WavePeriod)
        
        Me%Var%MSedimentedDT          = 1.3 * sqrt(DissipationRate/WaterDynamicVisc)            &
                                            * CSed_Sticking * OilConcentration                      &
                                            * Me%ExternalVar%SPM * (IntrusionDepth              &
                                            *Me%ExternalVar%Area)

        Me%Var%VSedimentedDT          = Me%Var%MSedimentedDT / Me%Var%Density
        

ifPC:   if ((Me%Var%OilEvaporation) .AND. (Me%Var%EvaporationMethod  .EQ. PseudoComponents)) then

do1:        do n=1,Me%Var%NbrPC 
       
cd2:           if  (Me%Var%VolPC(n)-(Me%Var%VSedimentedDT/Me%Var%NbrPC)                     &
                    * Me%Var%DTOilInternalProcesses .GT. 0) then 

                    Me%Var%VolPC(n)       = Me%Var%VolPC(n)-(Me%Var%VSedimentedDT           &
                                                /Me%Var%NbrPC) * Me%Var%DTOilInternalProcesses

               else cd2

                    Me%Var%VolPC(n) = 0.0

               end if   cd2

            end do do1

        end if ifPC

cd3:    if (Me%Var%MassOil - (Me%Var%MSedimentedDT) * Me%Var%DTOilInternalProcesses     &
            .GT. 0) then
        
            Me%Var%MassOil            = Me%Var%MassOil - (Me%Var%MSedimentedDT)         &
                                            * Me%Var%DTOilInternalProcesses
            Me%Var%VolumeOil          = max(0.0,Me%Var%VolumeOil - (Me%Var%VSedimentedDT)       &
                                            * Me%Var%DTOilInternalProcesses)
        else    cd3
            Me%Var%MSedimentedDT      = Me%Var%MassOil / Me%Var%DTOilInternalProcesses
            Me%Var%VSedimentedDT      = Me%Var%MSedimentedDT / Me%Var%Density
            Me%Var%MassOil            = 0.0
            Me%Var%VolumeOil          = 0.0
        end if  cd3


        Me%Var%Volume                 = PETIT + Me%Var%VolumeOil /  (1 - Me%Var%VWaterContent)

        Me%Var%SlickThickness         = Me%Var%Volume / (max(AllmostZero,Me%ExternalVar%Area))
        Me%Var%OilThickness           = Me%Var%VolumeOil / (max(AllmostZero,Me%ExternalVar%Area))


        Me%Var%MSedimented            = Me%Var%MSedimented + Me%Var%MSedimentedDT       &
                                            * Me%Var%DTOilInternalProcesses
        Me%Var%VSedimented            = Me%Var%VSedimented  + Me%Var%VSedimentedDT      &
                                            * Me%Var%DTOilInternalProcesses
        Me%Var%FMSedimented           = Me%Var%MSedimented/Me%Var%MassINI



    end subroutine Sedimentation

        !---------------------------------------------------------------------------

      
    subroutine Dissolution

        !Arguments---------------------------------------------------------------

        !Internal----------------------------------------------------------------
        integer         ::  n

        !------------------------------------------------------------------------

        Me%Var%SolubilityOilInWater = Me%Var%SolubilityOilInWater * exp(-CDiss_decayment*   &
                                      ((Me%Var%DTOilInternalProcesses)/3600.0))

        Me%Var%MDissolvedDT         = (CDiss_KTransfMass/3600) * (1 - Me%Var%MWaterContent) &
                                           * Me%ExternalVar%Area * Me%Var%SolubilityOilInWater

        Me%Var%VDissolvedDT         = Me%Var%MDissolvedDT / Me%Var%Density


ifPC:   if ((Me%Var%OilEvaporation) .AND. (Me%Var%EvaporationMethod  .EQ. PseudoComponents)) then

do1:        do n=1,Me%Var%NbrPC 
       
cd1:            if  (Me%Var%VolPC(n) - (Me%Var%VDissolvedDT/Me%Var%NbrPC)                   &
                     * Me%Var%DTOilInternalProcesses .GT. 0) then 

                    Me%Var%VolPC(n)     = Me%Var%VolPC(n)-(Me%Var%VDissolvedDT              &
                                              / Me%Var%NbrPC) * Me%Var%DTOilInternalProcesses
                else    cd1
                
                    Me%Var%VolPC(n)     = 0.0
            
                end if  cd1
        
            end do  do1

        end if ifPC

cd2:    if (Me%Var%MassOil - (Me%Var%MdissolvedDT) * Me%Var%DTOilInternalProcesses      &
            .GT. 0) then

            Me%Var%MassOil          = Me%Var%MassOil - (Me%Var%MdissolvedDT)            &
                                          * Me%Var%DTOilInternalProcesses
            Me%Var%VolumeOil        = max(0.0,Me%Var%VolumeOil - (Me%Var%VDissolvedDT)          &
                                          * Me%Var%DTOilInternalProcesses)
        else    cd2 
            Me%Var%MDissolvedDT     = Me%Var%MassOil / Me%Var%DTOilInternalProcesses
            Me%Var%VDissolvedDT     = Me%Var%MDissolvedDT / Me%Var%Density
            Me%Var%MassOil          = 0.0
            Me%Var%VolumeOil        = 0.0
        end if  cd2

        Me%Var%Volume               = PETIT + Me%Var%VolumeOil                              &
                                          /(1 - Me%Var%VWaterContent)

        Me%Var%SlickThickness       = Me%Var%Volume / Me%ExternalVar%Area
        Me%Var%OilThickness         = Me%Var%VolumeOil / Me%ExternalVar%Area


        Me%Var%MDissolved           = Me%Var%MDissolved + Me%Var%MDissolvedDT           &
                                          * Me%Var%DTOilInternalProcesses
        Me%Var%VDissolved           = Me%Var%VDissolved + Me%Var%VDissolvedDT           &
                                          * Me%Var%DTOilInternalProcesses
        Me%Var%FMDissolved          = Me%Var%MDissolved / Me%Var%MassINI

    end subroutine Dissolution                                      

    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------


    subroutine Emulsification

        !Arguments---------------------------------------------------------------
        
     
        !------------------------------------------------------------------------

 
cd1:    if  (Me%Var%FMEvaporated  .GE. (Me%Var%Cemuls/100)) then
        
cd2:        if (Me%Var%EmulsificationMethod .EQ. Rasmussen) then
            
                Me%Var%VWaterContentDT = ((CEmuls_1 / (Me%Var%ViscINI/1000.))               &
                                             * (1 + Me%ExternalVar%Wind)**2                     &
                                             * ((Me%Var%MaxVWaterContent/100.)                  &
                                             - Me%Var%VWaterContent))                           &
                                             - ((Me%Var%VWaterContent * CEmuls_2)               &
                                             /(Me%Var%AsphalteneContent * Me%Var%WaxContent & 
                                             * (Me%Var%ViscINI/1000.)))                          
            
            else if (Me%Var%EmulsificationMethod .EQ. Mackay) then   cd2
            
                Me%Var%VWaterContentDT = Me%Var%EmulsParameter * (1+Me%ExternalVar%Wind)**2        &
                                             * (1 - Me%Var%VWaterContent                        &
                                             / (Me%Var%MaxVWaterContent/100.))
            
            end if cd2
        
        else cd1
            
            Me%Var%VWaterContentDT     = 0.0
        
        end if cd1

        Me%Var%MWaterContentDT = Me%Var%VWaterContentDT / (Me%Var%Density               &
                                     / Me%ExternalVar%WaterDensity) 

        Me%Var%MWaterContent   = Me%Var%MWaterContent + (Me%Var%MWaterContentDT         &
                                     * Me%Var%DTOilInternalProcesses)
        Me%Var%VWaterContent   = Me%Var%VWaterContent + (Me%Var%VWaterContentDT         &
                                     * Me%Var%DTOilInternalProcesses)
        

        
        Me%Var%Volume          = PETIT + Me%Var%VolumeOil /  (1 - Me%Var%VWaterContent)

        Me%Var%SlickThickness  = Me%Var%Volume / (max(AllmostZero,Me%ExternalVar%Area))
        Me%Var%OilThickness    = Me%Var%VolumeOil / (max(AllmostZero,Me%ExternalVar%Area))




        !------------------------------------------------------------------------

    end subroutine Emulsification

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------


    subroutine ChemDispersion

        !Arguments---------------------------------------------------------------


        !Internal----------------------------------------------------------------

        integer :: n
        real    :: ChemDispersionTime
        !------------------------------------------------------------------------
  
cd1:    if ((Me%ExternalVar%Now .GT. Me%Var%Start_ChemDispersion) .AND.                     &
            (Me%ExternalVar%Now .LE. Me%Var%End_ChemDispersion))  then
            
cd2:        if (Me%Var%MassOilStartChemDisp .LT. Me%Var%MassOil)    then
                
                Me%Var%MassOilStartChemDisp  = Me%Var%MassOil
                ChemDispersionTime               = Me%Var%End_ChemDispersion                         &
                                                   - Me%Var%Start_ChemDispersion
            
                Me%Var%MChemDispersedDT      = (Me%Var%P_AreaSprayed/100.0)  &
                                                   * Me%Var%MassOilStartChemDisp & 
                                                   *(Me%Var%Efficiency/100.0) /  ChemDispersionTime
            end if  cd2

            Me%Var%VChemDispersedDT          = Me%Var%MChemDispersedDT / Me%Var%Density

ifPC:       if ((Me%Var%OilEvaporation) .AND. (Me%Var%EvaporationMethod  .EQ. PseudoComponents)) then

do1:            do n=1,Me%Var%NbrPC 
       
cd3:                if  (Me%Var%VolPC(n) - (Me%Var%VChemDispersedDT/Me%Var%NbrPC)           &
                         * Me%Var%DTOilInternalProcesses .GT. 0) then 

                        Me%Var%VolPC(n)      = Me%Var%VolPC(n)-(Me%Var%VChemDispersedDT      &
                                                  / Me%Var%NbrPC) * Me%Var%DTOilInternalProcesses
                    else    cd3
                
                        Me%Var%VolPC(n)      = 0.0
            
                    end if  cd3
        
                end do  do1

            end if ifPC

cd4:        if (Me%Var%MassOil - (Me%Var%MChemDispersedDT)                                  &
                * Me%Var%DTOilInternalProcesses .GT. 0) then

                Me%Var%MassOil          = Me%Var%MassOil - (Me%Var%MChemDispersedDT)    &
                                              * Me%Var%DTOilInternalProcesses
                Me%Var%VolumeOil        = max(0.0,Me%Var%VolumeOil - (Me%Var%VChemDispersedDT)  &
                                              * Me%Var%DTOilInternalProcesses)
            else    cd4 
                Me%Var%MChemDispersedDT = Me%Var%MassOil / Me%Var%DTOilInternalProcesses
                Me%Var%VChemDispersedDT = Me%Var%MChemDispersedDT / Me%Var%Density
                Me%Var%MassOil          = 0.0
                Me%Var%VolumeOil        = 0.0
            end if  cd4

            Me%Var%Volume               = PETIT + Me%Var%VolumeOil                         &
                                              /(1 - Me%Var%VWaterContent)

            Me%Var%SlickThickness       = Me%Var%Volume / Me%ExternalVar%Area
            Me%Var%OilThickness         = Me%Var%VolumeOil / Me%ExternalVar%Area


            Me%Var%MChemDispersed       = Me%Var%MChemDispersed + Me%Var%MChemDispersedDT  &
                                              * Me%Var%DTOilInternalProcesses

            Me%Var%VChemDispersed       = Me%Var%VChemDispersed + Me%Var%VChemDispersedDT  &
                                              * Me%Var%DTOilInternalProcesses

            Me%Var%FMChemDispersed      = Me%Var%MChemDispersed / Me%Var%MassINI
        
        else cd1
            
            Me%Var%MChemDispersedDT     = 0.0 

            Me%Var%VChemDispersedDT     = 0.0
        
        end if  cd1
            
    end subroutine ChemDispersion
  
    !----------------------------------------------------------------------------
    
    !----------------------------------------------------------------------------

    subroutine MecCleanup

        !Arguments---------------------------------------------------------------

        !Internal----------------------------------------------------------------

        integer :: n
        !------------------------------------------------------------------------
  
cd1:    if ((Me%ExternalVar%Now .GT. Me%Var%Start_Mec_Cleanup) .AND. &
            (Me%ExternalVar%Now .LE. Me%Var%End_Mec_Cleanup))  then
            
            Me%Var%VOilRecoveredDT     = Me%Var%VEmulsionRecoveryRate           & 
                                             * (1 - Me%Var%VWaterContent)

            Me%Var%MOilRecoveredDT     = Me%Var%VOilRecoveredDT / Me%Var%Density        

ifPC:       if ((Me%Var%OilEvaporation) .AND. (Me%Var%EvaporationMethod  .EQ. PseudoComponents)) then

do1:            do n=1,Me%Var%NbrPC 
       
cd3:                if  (Me%Var%VolPC(n) - (Me%Var%VOilRecoveredDT/Me%Var%NbrPC)            &
                         * Me%Var%DTOilInternalProcesses .GT. 0) then 

                        Me%Var%VolPC(n)    = Me%Var%VolPC(n)-(Me%Var%VOilRecoveredDT        &
                                                 / Me%Var%NbrPC) * Me%Var%DTOilInternalProcesses
                    else    cd3
                
                        Me%Var%VolPC(n)    = 0.0
            
                    end if  cd3
        
                end do  do1

            end if ifPC

cd4:        if (Me%Var%MassOil - Me%Var%MOilRecoveredDT * Me%Var%DTOilInternalProcesses &
                .GT. 0) then

                Me%Var%MassOil         = Me%Var%MassOil - (Me%Var%MOilRecoveredDT)      &
                                             * Me%Var%DTOilInternalProcesses
              
                Me%Var%VolumeOil       = max(0.0,Me%Var%VolumeOil - (Me%Var%VOilRecoveredDT)    &
                                             * Me%Var%DTOilInternalProcesses)
            else    cd4 
                Me%Var%MOilRecoveredDT = Me%Var%MassOil / Me%Var%DTOilInternalProcesses
                Me%Var%VOilRecoveredDT = Me%Var%MOilRecoveredDT / Me%Var%Density
                Me%Var%MassOil         = 0.0
                Me%Var%VolumeOil       = 0.0
            end if  cd4

            Me%Var%Volume              = PETIT + Me%Var%VolumeOil                           &
                                             /(1 - Me%Var%VWaterContent)
            
            Me%Var%SlickThickness      = Me%Var%Volume / Me%ExternalVar%Area
            Me%Var%OilThickness        = Me%Var%VolumeOil / Me%ExternalVar%Area


            Me%Var%MOilRecovered       = Me%Var%MOilRecovered + Me%Var%MOilRecoveredDT  &
                                             * Me%Var%DTOilInternalProcesses
            
            Me%Var%VOilRecovered       = Me%Var%VOilRecovered + Me%Var%VOilRecoveredDT  &
                                             * Me%Var%DTOilInternalProcesses
            
            Me%Var%FMOilRecovered      = Me%Var%MOilRecovered / Me%Var%MassINI
        
        else cd1
            
            Me%Var%MOilRecoveredDT     = 0.0 

            Me%Var%VOilRecoveredDT     = 0.0
        
        end if  cd1
            
    end subroutine MecCleanup
  
    !----------------------------------------------------------------------------
    
    !----------------------------------------------------------------------------

    real function F_MaxVWaterContent ()

        !Arguments-------------------------------------------------------------
    

cd1 :   if  (Me%Var%MaxVWaterContent .LT. 0.0) then
    
cd2 :       if (Me%Var%OilType .EQ. refined) then

                F_MaxVWaterContent      = 0.0

            else cd2

                F_MaxVWaterContent      = 70.0


            end if cd2
        
        else cd1

            F_MaxVWaterContent             = Me%Var%MaxVWaterContent
        
        end if cd1

    end function F_MaxVWaterContent

    !----------------------------------------------------------------------------
 
    !------------------------------------------------------------------------

    real function F_Cdisp ()

        !Arguments---------------------------------------------------------------

        !------------------------------------------------------------------------

!        F_Cdisp = max(0.0 , -390.07*log(Me%Var%ViscCin/1.0e6) - 2779.4 )

         F_Cdisp = max(0.0 , -312.25*log(Me%Var%ViscCin) + 2509.8 )

        !------------------------------------------------------------------------

    end function F_Cdisp

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    
    real function F_IBP ()

        !Arguments---------------------------------------------------------------

        !------------------------------------------------------------------------

        if (Me%Var%OilType .EQ. Refined) then
           F_IBP = 654.45 - 4.6588 * Me%Var%API
        else
            F_IBP = 532.98 - 3.1295 * Me%Var%API
        end if            

        !como as equações anteriores dao resultados mto baixos, usa-se uma 
        !aproximacao de mackay 1980
        
!        f_IBP = 542.6-30.275*Me%Var%API+1.565*Me%Var%API**2-0.03439*     &
!                Me%Var%API**3 + 0.0002604*Me%Var%API**4    
        !------------------------------------------------------------------------

    end function F_IBP

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    real function F_Slope ()

        !Arguments---------------------------------------------------------------


        !------------------------------------------------------------------------
!a próxima correlação é do ADIOS 1, para ser combinada com ktransfmass=2.5E-3
!esta correlação fornece valores muito altos para API>42
!        if (Me%Var%OilType .EQ. Refined) then
!            F_Slope = 388.19 - 3.8725 * Me%Var%API
!        else
!            F_Slope = 985.62 - 13.597 * Me%Var%API
!        end if

        !na nova versao do ADIOS em vez das correlacoes anteriores,e usada uma nova
        !correlacao para estimar a slope, e e a seguinte(deve ser usada com ktransfmass = 1.5E-3):
        
        F_Slope = 1356.7 - 247.36 * log(Me%Var%API)
        
        !se o petroleo for refinado (oiltype .eq. 'Refined')
        !a aproximacao anterior e grosseira.os produtos refinados muitas vezes 
        !demonstram curvas de destilacao nao-lineares.o melhor sera arranjar 
        !outra formula para estes casos.


        !------------------------------------------------------------------------

    end function F_Slope

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    real function F_CVisc_E ()

         !Arguments---------------------------------------------------------------
  
         !Internal----------------------------------------------------------------

         real            :: ViscCin15 = null_real

         !------------------------------------------------------------------------

         ViscCin15 = Me%Var%ViscCin * exp( CVisc_T*((1/(Temp15+273.15))                         &
                     -(1/(Me%ExternalVar%WaterTemperature+273.15))) )                             


!         if (ViscCin15  .LT.  0.5 )                                                                 &
!             F_CVisc_E =1.0
!         if (ViscCin15  .GT. 1500 )                                                                 &
!             F_CVisc_E =10.0
!         if ((ViscCin15 .GE.  0.5 ) .and. (ViscCin15 .LE. 1500.0 ))                                 &
!             F_CVisc_E =9.969e-1 + 6.002e-3 * ViscCin15
                    

         if (ViscCin15  .GT. 38. )                                                                    &
             F_CVisc_E =10.0
         if (ViscCin15 .LE. 38.  )                                                                 &
             F_CVisc_E =-0.0059*ViscCin15**2 + 0.4461*ViscCin15 + 1.413 





    end function F_CVisc_E


    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    
    subroutine EvapPropINI

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------
        
        integer :: STAT_CALL

        !Local-------------------------------------------------------------------

        REAL,DIMENSION(:),pointer::TDist
        REAL,DIMENSION(:),pointer::CFDist
        REAL,DIMENSION(:),pointer::BPPC
        REAL,DIMENSION(:),pointer::CPC
        REAL,DIMENSION(:),pointer::DeltaSPC
        integer                      ::n
        integer                      ::NbrPoints
        real                         ::TDistInic
        real                         ::dTDistdCFDist


    !------------------------------------------------------------------------




        Me%Var%TDistExp = Me%Var%TDistExp + 273.15
        Me%Var%CPDistExp= Me%Var%CPDistExp/100.0
cd1:    if (Me%Var%NbrDistCuts .LT. 2) then
                !se o petroleo for refinado (oiltype .eq. Refined)
                !a aproximacao e grosseira.os produtos refinados muitas vezes 
                !demonstram curvas de destilacao nao-lineares.o melhor sera arranjar 
                !outra formula para estes casos.

            TDistInic = 457.16 - 3.3447*Me%Var%API
   
            dTDistdCFDist = 1356.7 - 247.36 * log(Me%Var%API)

            Me%Var%NbrPC = 5

            allocate(BPPC(1:Me%Var%NbrPC), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                            &
                call SetError(FATAL_, INTERNAL_,                                                      &
                             "Subroutine EvapPropINI; Module ModuleOil. ERR01") 

            allocate(Me%Var%VolInicPC(1:Me%Var%NbrPC), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                            &
                call SetError(FATAL_, INTERNAL_,                                                      &
                             "Subroutine EvapPropINI; Module ModuleOil. ERR02") 

            do n=1,Me%Var%NbrPC
        
                BPPC(n)                 = TDistInic + dTDistdCFDist * ((n - 1/2)/5.0)
                Me%Var%VolInicPC(n) = Me%Var%VolInic / 5.0
    
            end do

        else  cd1

            NbrPoints = Me%Var%NbrDistCuts +2

            allocate(TDist(1:NbrPoints), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                            &
                call SetError(FATAL_, INTERNAL_,                                                      &
                             "Subroutine EvapPropINI; Module ModuleOil. ERR03") 

            allocate(CFDist(1:NbrPoints), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                            &
                call SetError(FATAL_, INTERNAL_,                                                      &
                             "Subroutine EvapPropINI; Module ModuleOil. ERR04") 


            TDist(1)   = Me%Var%TDistExp(1) - Me%Var%CPDistExp(1)                     &
                         * ((Me%Var%TDistExp(2) - Me%Var%TDistExp(1)) /                     &
                         (Me%Var%CPDistExp(2) - Me%Var%CPDistExp(1)))
            
            CFDist(1)  = 0.0

            do n=2,NbrPoints-1

                TDist(n) = Me%Var%TDistExp(n-1)
                CFDist(n)= Me%Var%CPDistExp(n-1)

            end do

            TDist(NbrPoints) = Me%Var%TDistExp(Me%Var%NbrDistCuts) +                        &
                               (1 - Me%Var%CPDistExp(Me%Var%NbrDistCuts))                   &
                               * ((Me%Var%TDistExp(Me%Var%NbrDistCuts) -                    &
                               Me%Var%TDistExp(1)) /                                            &
                               (Me%Var%CPDistExp(Me%Var%NbrDistCuts) -                      &
                               Me%Var%CPDistExp(1)))
            
            CFDist(NbrPoints)= 1.0

            Me%Var%NbrPC = NbrPoints -1

            allocate(BPPC(1:Me%Var%NbrPC), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                            &
                call SetError(FATAL_, INTERNAL_,                                                      &
                             "Subroutine EvapPropINI; Module ModuleOil. ERR05") 

            allocate(Me%Var%VolInicPC(1:Me%Var%NbrPC), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                            &
                call SetError(FATAL_, INTERNAL_,                                                      &
                             "Subroutine EvapPropINI; Module ModuleOil. ERR06") 

            do n=1,Me%Var%NbrPC

                BPPC(n)                 = TDist(n+1)
                Me%Var%VolInicPC(n) = (CFDist(n+1) - CFDist(n)) * Me%Var%VolInic

            end do

        end if  cd1

       allocate(CPC(1:Me%Var%NbrPC), STAT = STAT_CALL)
       if (STAT_CALL .NE. SUCCESS_)                                                                 &
          call SetError(FATAL_, INTERNAL_,                                                            &
                       "Subroutine EvapPropINI; Module ModuleOil. ERR07") 

       allocate(DeltaSPC(1:Me%Var%NbrPC), STAT = STAT_CALL)
       if (STAT_CALL .NE. SUCCESS_)                                                                 &
           call SetError(FATAL_, INTERNAL_,                                                           &
                        "Subroutine EvapPropINI; Module ModuleOil. ERR08") 

       allocate(Me%Var%VmrelPC(1:Me%Var%NbrPC), STAT = STAT_CALL)
       if (STAT_CALL .NE. SUCCESS_)                                                                 &
           call SetError(FATAL_, INTERNAL_,                                                           &
                             "Subroutine EvapPropINI; Module ModuleOil. ERR09") 

       allocate(Me%Var%MWPC(1:Me%Var%NbrPC), STAT = STAT_CALL)
       if (STAT_CALL .NE. SUCCESS_)                                                                 &
           call SetError(FATAL_, INTERNAL_,                                                           &
                             "Subroutine EvapPropINI; Module ModuleOil. ERR09a") 

       allocate(Me%Var%PvapPC(1:Me%Var%NbrPC), STAT = STAT_CALL)
       if (STAT_CALL .NE. SUCCESS_)                                                                 &    
           call SetError(FATAL_, INTERNAL_,                                                           &
                         "Subroutine EvapPropINI; Module ModuleOil. ERR10") 

       allocate(Me%Var%VolPC(1:Me%Var%NbrPC), STAT = STAT_CALL)
       if (STAT_CALL .NE. SUCCESS_)                                                                 &
           call SetError(FATAL_, INTERNAL_,                                                           &
                         "Subroutine EvapPropINI; Module ModuleOil. ERR11") 

        do n=1,Me%Var%NbrPC
        
            CPC(n)                  = 0.19 * BPPC(n) - 18.0
            DeltaSPC(n)             = 8.75 + 1.987 * log10(BPPC(n))
            Me%Var%VmrelPC(n)   = 7e-5 - (2.102e-7 * BPPC(n)) + (1e-9 * BPPC(n)**2.0)
            Me%Var%MWPC(n)      = 0.04132 - (1.985e-4 * BPPC(n)) + (9.494e-7 * BPPC(n)**2.0)

            !AtmosphericPressure in Pa
            Me%Var%PvapPC(n)    = (Me%ExternalVar%AtmosphericPressure)          &
                                      * exp(((1.0/(BPPC(n)-CPC(n))) - (1.0                          &
                                      /((Me%ExternalVar%WaterTemperature+273.15)-CPC(n))))      &
                                      * DeltaSPC(n)*(BPPC(n)-CPC(n))**2.0/(CEvap_deltaZ*R*BPPC(n)))
            Me%Var%VolPC(n)     = Me%Var%VolInicPC(n) 
        
        end do


       allocate(Me%Var%VEvaporatedPCDT(1:Me%Var%NbrPC), STAT = STAT_CALL)
       if (STAT_CALL .NE. SUCCESS_)                                                                 &
           call SetError(FATAL_, INTERNAL_,                                                           &
                         "Subroutine EvapPropINI; Module ModuleOil. ERR12") 


        deallocate(BPPC, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                                  &
            call SetError(FATAL_, INTERNAL_,                                                          &
                          "Subroutine EvapPropINI; Module ModuleOil. ERR13") 
        nullify   (BPPC)

        deallocate(TDist, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                                  &
            call SetError(FATAL_, INTERNAL_,                                                          &
                          "Subroutine EvapPropINI; Module ModuleOil. ERR14") 
        nullify   (TDist)

        deallocate(CFDist, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                                  &
            call SetError(FATAL_, INTERNAL_,                                                          &
                          "Subroutine EvapPropINI; Module ModuleOil. ERR15") 
        nullify   (CFDist)

        deallocate(CPC, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                                  &
            call SetError(FATAL_, INTERNAL_,                                                          &
                          "Subroutine EvapPropINI; Module ModuleOil. ERR16") 
        nullify   (CPC)

        deallocate(DeltaSPC, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                                  &
            call SetError(FATAL_, INTERNAL_,                                                          &
                          "Subroutine EvapPropINI; Module ModuleOil. ERR17") 
        nullify   (DeltaSPC)

    end subroutine EvapPropINI


    !----------------------------------------------------------------------------


    subroutine TimeSerieOutput(DataLineIN)

        !Arguments---------------------------------------------------------------
        real,    optional, dimension(:), pointer :: DataLineIN

        !External--------------------------------------------------------------
        
        integer :: STAT_CALL

        !Local-------------------------------------------------------------------

        integer :: aux, Prop

        !------------------------------------------------------------------------

        Me%TimeSerie%DataLine(ColMassOil       ) = Me%Var%MassOil
        Me%TimeSerie%DataLine(ColVolOilBeached ) = Me%ExternalVar%VolOilBeached
        Me%TimeSerie%DataLine(ColVolumeBeached ) = Me%ExternalVar%VolumeBeached
        Me%TimeSerie%DataLine(ColVolumeOil     ) = Me%Var%VolumeOil
        Me%TimeSerie%DataLine(ColVolume        ) = Me%Var%Volume
        Me%TimeSerie%DataLine(ColArea          ) = Me%ExternalVar%Area
        Me%TimeSerie%DataLine(ColAreaTeoric    ) = Me%Var%AreaTeoric
        Me%TimeSerie%DataLine(ColThickness     ) = Me%Var%OilThickness * 1000.
        Me%TimeSerie%DataLine(ColMEvaporated   ) = Me%Var%MEvaporated
        Me%TimeSerie%DataLine(ColVEvaporated   ) = Me%Var%VEvaporated
        Me%TimeSerie%DataLine(ColFMEvaporated  ) = Me%Var%FMEvaporated
        Me%TimeSerie%DataLine(ColMDispersed    ) = Me%Var%MDispersed
        Me%TimeSerie%DataLine(ColVDispersed    ) = Me%Var%VDispersed
        Me%TimeSerie%DataLine(ColFMDispersed   ) = Me%Var%FMDispersed
        Me%TimeSerie%DataLine(ColMSedimented   ) = Me%Var%MSedimented 
        Me%TimeSerie%DataLine(ColVSedimented   ) = Me%Var%VSedimented
        Me%TimeSerie%DataLine(ColFMSedimented  ) = Me%Var%FMSedimented
        Me%TimeSerie%DataLine(ColMDissolved    ) = Me%Var%MDissolved
        Me%TimeSerie%DataLine(ColVDissolved    ) = Me%Var%VDissolved
        Me%TimeSerie%DataLine(ColFMDissolved   ) = Me%Var%FMDissolved
        Me%TimeSerie%DataLine(ColMChemDispersed) = Me%Var%MChemDispersed
        Me%TimeSerie%DataLine(ColVChemDispersed) = Me%Var%VChemDispersed
        Me%TimeSerie%DataLine(ColFMChemDispersed)= Me%Var%FMChemDispersed
        Me%TimeSerie%DataLine(ColMOilRecovered)  = Me%Var%MOilRecovered
        Me%TimeSerie%DataLine(ColVOilRecovered)  = Me%Var%VOilRecovered
        Me%TimeSerie%DataLine(ColFMOilRecovered) = Me%Var%FMOilRecovered
        Me%TimeSerie%DataLine(ColMWaterContent ) = Me%Var%MWaterContent
        Me%TimeSerie%DataLine(ColVWaterContent ) = Me%Var%VWaterContent
        Me%TimeSerie%DataLine(ColDensity       ) = Me%Var%Density
        Me%TimeSerie%DataLine(ColViscosity     ) = Me%Var%Viscosity

cd4 :   if (present(DataLineIN)) then
            aux = ColNbr + size(DataLineIN)

do1 :       do Prop = (ColNbr+1), aux
                Me%TimeSerie%DataLine(Prop) = DataLineIN(aux-ColNbr)
            end do do1
        end if cd4


        call WriteTimeSerieLine(Me%ObjTimeSerie,                                &
                                DataLine     = Me%TimeSerie%DataLine,           &
                                STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
                call SetError(FATAL_, INTERNAL_,                                &
                             "Subroutine TimeSerieOutput; Module ModuleOil. ERR01") 

        !------------------------------------------------------------------------

    end subroutine TimeSerieOutput

    !----------------------------------------------------------------------------



    !----------------------------------------------------------------------------
    !Influence on movement, growth, etc.


    subroutine OilActiveProcesses(OilID,                          &
                                  GridThickness,                  &
                                  WaterTemperature,               &
                                  WaterDensity,                   &
                                  VolInic,                        &
                                  DT,                             &
                                  AreaTotal,                      &
                                  STAT)


        !Arguments---------------------------------------------------------------
        integer                         :: OilID
        real, intent(in)                :: WaterTemperature
        real, intent(in)                :: WaterDensity
        real, intent(in)                :: DT

        real, intent(in)                :: VolInic
        real, intent(in), optional      :: AreaTotal 
        integer, optional, intent(OUT)  :: STAT

        real, pointer, dimension(:,:)   :: GridThickness


        !External----------------------------------------------------------------

        integer :: ready_ 
        integer :: STAT_CALL 

        !Local-------------------------------------------------------------------

        integer :: STAT_ 

        integer :: ILB, IUB 
        integer :: JLB, JUB
        integer :: KLB, KUB
        integer :: I,   J
               
        real    :: CoefVelMancha
        real    :: Delta
        real    :: TransitionalTime
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(OilID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then


            ! Actualized the time
            call GetComputeCurrentTime(Me%ObjTime,                          &
                                       Me%ExternalVar%Now,                  &
                                       STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                          &
                call SetError(FATAL_, INTERNAL_,                                  &
                             "Subroutine OilActiveProcesses; Module ModuleOil. ERR01") 




            ILB = Me%ExternalVar%WorkSize%ILB
            IUB = Me%ExternalVar%WorkSize%IUB

            JLB = Me%ExternalVar%WorkSize%JLB
            JUB = Me%ExternalVar%WorkSize%JUB

            KLB = Me%ExternalVar%WorkSize%KLB
            KUB = Me%ExternalVar%WorkSize%KUB

            

cd2 :       if (Me%ExternalVar%Now .GE. Me%NextActiveComputeTime) then
                
                
                Me%ExternalVar%GridThickness => GridThickness
                
                if (.NOT. associated(Me%ExternalVar%GridThickness))             &
                    call SetError(FATAL_, INTERNAL_,                                  &
                                 "Subroutine OilActiveProcesses; Module ModuleOil. ERR02") 



                call GetComputeFaces3D (Me%ObjMap,                                         &
                                        ComputeFacesU3D = Me%ExternalVar%ComputeFacesU3D,  &
                                        ComputeFacesV3D = Me%ExternalVar%ComputeFacesV3D,  &
                                        STAT            = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_)                                              &
                    call SetError(FATAL_, INTERNAL_,                                      &
                                 "Subroutine OilActiveProcesses; Module ModuleOil. ERR03") 
                


                call GetHorizontalGrid(Me%ObjHorizontalGrid,                &
                                      DZX  = Me%ExternalVar%DZX,            &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                      &
                    call SetError(FATAL_, INTERNAL_,                              &
                                 "Subroutine OilActiveProcesses; Module ModuleOil. ERR04") 


                call GetHorizontalGrid(Me%ObjHorizontalGrid,                &
                                      DZY  = Me%ExternalVar%DZY,            &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                      &
                    call SetError(FATAL_, INTERNAL_,                              &
                                 "Subroutine OilActiveProcesses; Module ModuleOil. ERR05") 

                    Me%ExternalVar%WaterTemperature = WaterTemperature


cd3:            if (Me%State%FirstStepAP)  then                                     

cd3a:                if (Me%Var%PourPoint .GT. Me%ExternalVar%WaterTemperature) then
                        write (*, *) 'Water Temperature is colder than Pour Point:'
                        write (*, *) 'Calculations may be unreliable'

                     end if    cd3a

                    Me%Var%VolInic   = VolInic
                    
                    call OilPropIni
                    ! esta rotina serve para inicializar uma série de variáveis 
                    !q nao puderam ser inicializadas no constructor, devido a 
                    !nao ser conhecida a WaterTemperature

                                    
cd4:                if (present(AreaTotal)) then
    
                        Me%ExternalVar%Area          = AreaTotal
                    
                    else    cd4
                    
                        Me%ExternalVar%Area          = F_FayArea(VolInic = Me%Var%VolInic, API = Me%Var%API )                

                    end if  cd4

                    Me%Var%OilThickness              = Me%Var%Volinic / Me%ExternalVar%Area
                    Me%Var%SlickThickness            = Me%Var%Volinic / Me%ExternalVar%Area

                end if  cd3



                Me%Var%SpreadingVelocityX(:,:) = 0.0
                Me%Var%SpreadingVelocityY(:,:) = 0.0

                call AreaTeoric(DT,WaterDensity,Delta)

cd6:            if (Me%Var%OilSpreading ) then

cd7:                if (Me%Var%SpreadingMethod .EQ. ThicknessGradient_) then

cd8:                    if (Me%Var%OilThickness .GE. Me%Var%ThicknessLimit) then

                            CoefVelMancha = Me%Var%UserCoefVelMancha * Gravity * (Delta**(1.6)) *     &
                                            (Me%Var%VolInic* Me%Var%VolInic /          &
                                            SQRT(WaterCinematicVisc))**(1./6.)

                            do j = JLB, JUB + 1
                            do i = ILB, IUB

                                if (Me%ExternalVar%ComputeFacesU3D (i, j, KUB) == Compute) then

                                    Me%Var%SpreadingVelocityX(i, j) =  CoefVelMancha  *  &
                                            (GridThickness(i, j-1) - GridThickness(i, j)) /  &
                                             Me%ExternalVar%DZX(i,  j-1)
                                else

                                    Me%Var%SpreadingVelocityX(i, j) = 0.

                                endif

                            enddo
                            enddo

                            do j = JLB, JUB
                            do i = ILB, IUB + 1

                                if (Me%ExternalVar%ComputeFacesV3D (i, j, KUB) == Compute) then

                                    Me%Var%SpreadingVelocityY(i, j) = CoefVelMancha   *  &
                                            (GridThickness(i-1, j) - GridThickness(i, j)) /  &
                                             Me%ExternalVar%DZY(i,  j-1)
                                else

                                    Me%Var%SpreadingVelocityY(i, j) = 0.

                                endif

                            enddo
                            enddo



!do1 :                       do J = JLB+1,JUB-1
!do2 :                       do I = ILB+1,IUB-1

!cd9 :                       if (Me%ExternalVar%OpenPoints3D(I,J,KUB) .EQ. 1) then 

                            
!                                Delta         = (WaterDensity - Me%Var%Density) / WaterDensity
!                                CoefVelMancha = Me%Var%UserCoefVelMancha * (Delta * VolInic * VolInic / &
!                                                 SQRT(WaterCinematicVisc))**(1./6.)

                                !As instrucoes comentadas referem-se a um calculo da derivada com 5 pontos.
    
                                !UMancha0=2*(Mancha(II,JJ-1,nPetroleo)-Mancha(II,JJ-2,nPetroleo))/(DX2D(II,JJ-2)+DX2D(II,JJ-1))
                                !UMancha3=2*(Mancha(II,JJ+2,nPetroleo)-Mancha(II,JJ+1,nPetroleo))/(DX2D(II,JJ)+DX2D(II,JJ+1))
!                                UMancha1 = 2.0 * (GridThickness(I,  J  ) - GridThickness(I,  J-1)) / &
!                                           (Me%ExternalVar%DUX(I,  J-1) + Me%ExternalVar%DUX(I,  J  ))
!                                UMancha2 = 2.0 * (GridThickness(I,  J+1) - GridThickness(I,  J  )) / &
!                                           (Me%ExternalVar%DUX(I,  J  ) + Me%ExternalVar%DUX(I,  J+1))
                                !UMancha=-(UMancha0+2*UMancha1+2*UMancha2+UMancha3)/6*CoefVelMancha
!                                Me%Var%SpreadingVelocityX(I,J) =-(UMancha1 + UMancha2) / 2.0 * CoefVelMancha
    
                                !VMancha0=2*(Mancha(II-1,JJ,nPetroleo)-Mancha(II-2,JJ,nPetroleo))/(DY2D(II-2,JJ)+DX2D(II-1,JJ))
                                 !VMancha3=2*(Mancha(II+2,JJ,nPetroleo)-Mancha(II+1,JJ,nPetroleo))/(DY2D(II+1,JJ)+DX2D(II+2,JJ))
!                                VMancha1 = 2.0 * (GridThickness(I,  J  ) - GridThickness(I-1,J  )) / &
!                                                 (Me%ExternalVar%DVY(I-1,J  ) + Me%ExternalVar%DVY(I,  J  ))
!                                VMancha2 = 2.0 * (GridThickness(I+1,J  ) - GridThickness(I,  J  )) / &
!                                                 (Me%ExternalVar%DVY(I,  J  ) + Me%ExternalVar%DVY(I+1,J  ))
                                !VMancha=-(VMancha0+2*VMancha1+2*VMancha2+VMancha3)/6*CoefVelMancha
!                                Me%Var%SpreadingVelocityY(I,J)=-(VMancha1 + VMancha2) / 2.0 * CoefVelMancha
    
                                !if (H(II,JJ+2)<-55..or.H(II,JJ-2)<-55.) UMancha=-(UMancha1+UMancha2)/2*CoefVelMancha
                                !if (H(II+2,JJ)<-55..or.H(II-2,JJ)<-55.) VMancha=-(VMancha1+VMancha2)/2*CoefVelMancha
    
!                                if (Me%ExternalVar%OpenPoints3D(I,  J+1,KUB) .EQ. 0)    &
!                                    Me%Var%SpreadingVelocityX(I,J) =-2.0 * (GridThickness(I,  J  ) -   &
!                                                                        GridThickness(I,  J-1)) /          &
!                                                                       (Me%ExternalVar%DUX(I,  J  ) +  &
!                                                                        Me%ExternalVar%DUX(I,  J-1)) * &
!                                                                        CoefVelMancha
    
!                                if (Me%ExternalVar%OpenPoints3D(I,  J-1,KUB) .EQ. 0)    &
!                                    Me%Var%SpreadingVelocityX(I,J) =-2.0 * (GridThickness(I,  J+1) -   &
!                                                                       GridThickness(I,  J  )) /           &
!                                                                       (Me%ExternalVar%DUX(I,  J+1) +  &
!                                                                        Me%ExternalVar%DUX(I,  J  )) * &
!                                                                        CoefVelMancha
    
!                                if (Me%ExternalVar%OpenPoints3D(I+1,J  ,KUB) .EQ. 0)    &
!                                    Me%Var%SpreadingVelocityY(I,J) =-2.0 * (GridThickness(I,  J  ) -   &
!                                                                        GridThickness(I-1,J  )) /          &
!                                                                       (Me%ExternalVar%DVY(I,  J  ) +  &
!                                                                        Me%ExternalVar%DVY(I-1,J  )) * &
!                                                                        CoefVelMancha
    
!                                if (Me%ExternalVar%OpenPoints3D(I-1,J  ,KUB) .EQ. 0)    &
!                                    Me%Var%SpreadingVelocityY(I,J) =-2.0 * (GridThickness(I+1,J  ) -   &
!                                                                         GridThickness(I,  J  )) /         &
!                                                                         (Me%ExternalVar%DVY(I+1,J  ) + &
!                                                                         Me%ExternalVar%DVY(I,  J  )) * &
!                                                                         CoefVelMancha


!                             end if cd9
                    
                        
!                             end do do2
!                            end do do1
                    
                        else    cd8
                    
                            Me%Var%SpreadingVelocityX(:,:) = 0.0
                            Me%Var%SpreadingVelocityY(:,:) = 0.0
                    
                        end if  cd8
                    
                    else if (Me%Var%SpreadingMethod .EQ. Fay_) then cd7
                
cd10:                    if (Me%Var%OilThickness .GE. Me%Var%ThicknessLimit) then

                             Me%Var%alfa         = (CFay_2*CFay_2/32.0) * &
                                                       (Delta*Gravity*Me%Var%VolInic*Me%Var%VolInic/          &
                                                       sqrt(WaterCinematicVisc))**(1.0/3.0)

cd11:                        if (Me%State%FirstStepAP) then

                                TransitionalTime     = (CFay_2/CFay_1)**4 * (Me%Var%VolInic /(Gravity * Delta *                 &
                                                       WaterCinematicVisc ))**(1./3.)

                                Me%Var%DiffCoef      = Me%Var%alfa * (1/sqrt(TransitionalTime)) 

                                Me%Var%alfa_old      = Me%Var%alfa

                            else   cd11

    
                                Me%Var%DiffCoef      = Me%Var%alfa * ((Me%Var%alfa_old / Me%Var%DiffCoef)**2 + DT)**(-1.0/2.0)

                                Me%Var%alfa_old      = Me%Var%alfa
      
                            end if  cd11


                            Me%Var%DiffVelocity = sqrt(2.0 * Me%Var%DiffCoef / DT)

                        else    cd10 
                            
                            Me%Var%DiffVelocity = 0.0
                        
                        end if  cd10
                    
                    end if cd7

                end if cd6




                Me%NextActiveComputeTime = Me%ExternalVar%Now + DT

                Me%State%FirstStepAP = OFF
                




                call UngetHorizontalGrid(Me%ObjHorizontalGrid,              &
                                         Me%ExternalVar%DZX,                &
                                         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                      &
                    call SetError(FATAL_, INTERNAL_,                              &
                                 "Subroutine OilActiveProcesses; Module ModuleOil. ERR06") 


                call UngetHorizontalGrid(Me%ObjHorizontalGrid,              &
                                         Me%ExternalVar%DZY,                &
                                         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                      &
                    call SetError(FATAL_, INTERNAL_,                              &
                                 "Subroutine OilActiveProcesses; Module ModuleOil. ERR07") 


                call UngetMap(Me%ObjMap,                                    &
                              Me%ExternalVar%ComputeFacesU3D,               &
                              STAT = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_)                                      &
                    call SetError(FATAL_, INTERNAL_,                              &
                                 "Subroutine OilActiveProcesses; Module ModuleOil. ERR08a") 
                call UngetMap(Me%ObjMap,                                    &
                              Me%ExternalVar%ComputeFacesV3D,               &
                              STAT = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_)                                      &
                    call SetError(FATAL_, INTERNAL_,                              &
                                 "Subroutine OilActiveProcesses; Module ModuleOil. ERR08b") 


                nullify(Me%ExternalVar%GridThickness)
            end if cd2

            call null_time   (Me%ExternalVar%Now)

            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                      &
            STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine OilActiveProcesses

    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------

    subroutine OilPropIni

        !Arguments---------------------------------------------------------------

        !------------------------------------------------------------------------


            Me%Var%Density15          = FreshWaterDensity15 * 141.5                    &
                                            / (131.5 + Me%Var%API)

            Me%Var%Density            = Me%Var%Density15 * (1 - CDens_DT           &
                                            * (Me%ExternalVar%WaterTemperature         &
                                            - Temp15)) 
            Me%Var%MassINI            = Me%Var%Density * Me%Var%VolInic
            Me%Var%MassOil            = Me%Var%MassINI
            Me%Var%VolumeOil          = Me%Var%VolInic



cd1:        if (Me%Var%ViscRef .LT. 0) then 
                Me%Var%DensTempViscRef = Me%Var%Density15 * (1 - 0.0008             &
                                             *(Me%Var%TempViscRef - Temp15))
                Me%Var%ViscRef         = (Me%Var%ViscCinRef/1.e6)                   &
                                             * Me%Var%DensTempViscRef * 1000     !cP
            
            end if cd1
                  
                       
            Me%Var%Viscosity           = 1000.*(Me%Var%ViscRef/1000.) * exp(CVisc_T *   &
                                             (1./(Me%ExternalVar%WaterTemperature+273.15)   &
                                             - 1./(Me%Var%TempViscRef+273.15))  )
            !a viscosidade anterior é calculada em cP
                  
            Me%Var%ViscCin             = 1.0e6 * (Me%Var%Viscosity/1000.)           &
                                             / Me%Var%Density !cSt
            !unidades - cSt
          
            Me%Var%ViscINI             = Me%Var%Viscosity !cP

            Me%Var%CVisc_E             = F_CVisc_E ()           

    end subroutine OilPropIni

    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------

    subroutine AreaTeoric(DT,WaterDensity,Delta)
   
        !Arguments---------------------------------------------------------------


        real, intent(in)       :: DT
        real, intent(in)       :: WaterDensity
        real, intent(out)      :: Delta

        !Local-------------------------------------------------------------------
        real                   :: aux1, aux2

        !------------------------------------------------------------------------

cd1:    if (Me%Var%OilThickness .GE. Me%Var%ThicknessLimit) then
           
           Delta                      = max(1e-9,(WaterDensity - Me%Var%Density)) / WaterDensity

           Me%Var%beta                       = Pi *(CFay_2**2/4.) * (Delta*Gravity*(Me%Var%VolInic**2)/          &
                                        sqrt(WaterCinematicVisc))**(1.0/3.0)
           
cd2:       if (Me%State%FirstStepAP) then
     
               Me%Var%AreaTeoric = Pi * (CFay_2**4/(CFay_1 * CFay_1)) * (1./4.) *                                  &
                                       (Me%Var%VolInic**5 * Gravity * Delta /                                                  &
                                       (WaterCinematicVisc * WaterCinematicVisc))**(1./6.)
           
               Me%Var%beta_old              = Me%Var%beta

           else cd2
           
               aux1 = Me%Var%beta_old

               aux2 = Me%Var%AreaTeoric/aux1

               Me%Var%AreaTeoric = Me%Var%beta * sqrt((aux2)**2 + DT)

               Me%Var%beta_old              = Me%Var%beta

           end if   cd2  
        
        else    cd1
          
           if (Me%State%FirstStepAP)                                        &
               Me%Var%AreaTeoric = Me%ExternalVar%Area
        
        end if  cd1

    end subroutine AreaTeoric
    
    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------

    real function F_ThicknessLimit ()

        
        !Arguments---------------------------------------------------------------

        !------------------------------------------------------------------------

        if (Me%Var%API .LT. 17.5 )                           &
            F_ThicknessLimit =1.0e-4
        if (Me%Var%API .GT. 45.0 )                           &
            F_ThicknessLimit =1.0e-5
        if ((Me%Var%API .GE. 17.5 ) .and. (Me%Var%API .LE. 45.0 ))    F_ThicknessLimit = 1.5727e-4 - 3.2727e-6 * Me%Var%API
            
        !------------------------------------------------------------------------
                                
    end function F_ThicknessLimit




!    !----------------------------------------------------------------------------
!
!    function F_CoefVelMancha(ObjOil)
!    real ::  F_CoefVelMancha
!
!        !Arguments---------------------------------------------------------------
!
!        type(T_Oil), pointer :: ObjOil  
!
!        !External----------------------------------------------------------------
!        
!        integer :: STAT_CALL
!
!        !------------------------------------------------------------------------
!
!        F_CoefVelMancha = 5E4
!
!        !------------------------------------------------------------------------
!
!    end function F_CoefVelMancha
!
!    !----------------------------------------------------------------------------



    !----------------------------------------------------------------------------

    !The model uses a modification of the well-known Fay spreading formulas to 
    ! compute the initial area of the slick. Fay broke the spreading process into 
    ! three stages. The first stage allows the oil to spread simply due to its 
    ! gravitational potential and occurs rather rapidly for even quite large spills. 
    !The value for t0 is quite short for most spills, usually 5 minutes or less. ADIOS2 
    ! assumes that, during this gravity-inertial spreading, none of the weathering 
    ! processes are taking place. In effect, the model has t0 as its starting time. 
    ! The initial area, A0, at the end of this gravity-inertial phase and the beginning
    ! of the weathering processes is computed to be

    function F_FayArea(VolInic, API, ObjEnterData)
    real ::  F_FayArea

        !Arguments---------------------------------------------------------------
        real, intent (IN)                           :: VolInic
        real, intent (IN), optional                 :: API
        integer, optional                           :: ObjEnterData

        !External----------------------------------------------------------------
        integer :: FromFile
        integer :: ExtractType

        real    :: lAPI

        !Internal----------------------------------------------------------------

        real            :: Delta
        real            :: Density15
        real            :: Density
        !tentar alterar as próximas constantes para variáveis
        real, parameter :: ConstantWaterDensity     = 1027.0
        real, parameter :: ConstantWaterTemperature = 18.0
        !------------------------------------------------------------------------


        if (present(API)) then

            lAPI = API

        else if (present(ObjEnterData)) then
            call GetExtractType(FromFile = FromFile)
            ExtractType = FromFile

            call OilOptionsAPI(ObjEnterData, ExtractType = ExtractType, API = lAPI)
        else
            call SetError(FATAL_, INTERNAL_, "Function  F_FayArea; Module ModuleOil. ERR00") 
        endif    


        Density15 = FreshWaterDensity15 * 141.5 / (131.5 + lAPI)
        
        Density   = Density15 * (1 - CDens_DT * (ConstantWaterTemperature - Temp15)) 

        !when oil is denser than surrounding water, oil will sink (in this case                     & 
        !oil is denser than an aproximation of surrounding water density
        if (Density .GT. ConstantWaterDensity)                                                      &
            call SetError(FATAL_, INTERNAL_,                                                          &
                          "Function  F_FayArea; Module ModuleOil. ERR01 ") 

        Delta     = (ConstantWaterDensity - Density) / ConstantWaterDensity

        
        F_FayArea = Pi * (CFay_2**4/(CFay_1 * CFay_1)) * (1./4.) *                                  &
                   (VolInic**5 * Gravity * Delta /                                                  &
                    (WaterCinematicVisc * WaterCinematicVisc))**(1./6.)

        !------------------------------------------------------------------------


    end function F_FayArea






    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillOil(OilID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: OilID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_ 
        integer                                     :: STAT_CALL
        integer                                     :: STAT_, nUsers 

        !---------------------------------------------------------------------- 

        STAT_ = UNKNOWN_
        
        call Ready(OilID, ready_) 

        if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mOIL_,  Me%InstanceID)

            if (nUsers == 0) then
                
                nUsers = DeassociateInstance (mTIME_,           Me%ObjTime          )
                if (nUsers == 0) stop 'KillOil - ModuleOil - ERR01'

                !nUsers = DeassociateInstance (mENTERDATA_,      Me%ObjEnterData     )
                !if (nUsers == 0) stop 'KillOil - ModuleOil - ERR01'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillOil - ModuleOil - ERR02'

                nUsers = DeassociateInstance (mGEOMETRY_,       Me%ObjGeometry      )
                if (nUsers == 0) stop 'KillOil - ModuleOil - ERR03'

                nUsers = DeassociateInstance (mMAP_,            Me%ObjMap           )
                if (nUsers == 0) stop 'KillOil - ModuleOil - ERR04'



                deallocate(Me%Var%SpreadingVelocityX)
                deallocate(Me%Var%SpreadingVelocityY)

                if (associated(Me%Var%VolInicPC)) then
                    deallocate(Me%Var%VolInicPC)
                endif

                if (associated(Me%Var%VmrelPC)) then
                    deallocate(Me%Var%VmrelPC)
                endif

                if (associated(Me%Var%MWPC)) then
                    deallocate(Me%Var%MWPC)
                endif

                if (associated(Me%Var%PvapPC)) then
                    deallocate(Me%Var%PvapPC)
                endif

                if (associated(Me%Var%VEvaporatedPCDT)) then
                    deallocate(Me%Var%VEvaporatedPCDT)
                endif

                if (associated(Me%Var%VolPC)) then
                    deallocate(Me%Var%VolPC)
                endif


                if (associated(Me%Var%TDistExp)) then
                    deallocate(Me%Var%TDistExp)
                endif

                if (associated(Me%Var%CPDistExp)) then
                    deallocate(Me%Var%CPDistExp)
                endif


cd3 :           if (Me%State%TimeSerie) then
                    deallocate(Me%TimeSerie%PropertyList)
                    deallocate(Me%TimeSerie%DataLine)

                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (nUsers == 0) stop 'KillOil - ModuleOil - ERR05'
                end if cd3

                !Deallocates Instance
                call DeallocateInstance ()

                OilID  = 0
                STAT_  = SUCCESS_

            end if

        else 
            STAT_ = ready_
        end if 


        if (present(STAT)) STAT = STAT_


        !------------------------------------------------------------------------

    end subroutine KillOil

    !--------------------------------------------------------------------------


    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Oil), pointer          :: AuxObjOil
        type (T_Oil), pointer          :: PreviousObjOil

        !Updates pointers
        if (Me%InstanceID == FirstObjOil%InstanceID) then
            FirstObjOil    => FirstObjOil%Next
        else
            PreviousObjOil => FirstObjOil
            AuxObjOil      => FirstObjOil%Next
            do while (AuxObjOil%InstanceID /= Me%InstanceID)
                PreviousObjOil => AuxObjOil
                AuxObjOil      => AuxObjOil%Next
            enddo

            !Now update linked list
            PreviousObjOil%Next => AuxObjOil%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------

    subroutine WriteFinalOil(OilID, UnitID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: OilID
        integer                                     :: UnitID
        integer, optional                           :: STAT
 
        !Local-------------------------------------------------------------------
        integer                                     :: ready_ 
        integer                                     :: STAT_ 

        !----------------------------------------------------------------------
       
        STAT_ = UNKNOWN_

        call Ready(OilID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            write (UnitID) Me%Var%OilType
            write (UnitID) Me%Var%API
            write (UnitID) Me%Var%PourPoint

            write (UnitID) Me%State%FirstStepAP
            write (UnitID) Me%State%FirstStepIP
            write (UnitID) Me%Var%DTOilInternalProcesses

            !Mass
            write (UnitID) Me%Var%MassOil
            write (UnitID) Me%Var%MassINI

            !Volume
            write (UnitID) Me%Var%VolumeOil
            write (UnitID) Me%Var%VolInic

            !Density
            write (UnitID) Me%Var%Density
            write (UnitID) Me%Var%Density15

            !Viscosity
            write (UnitID) Me%Var%ViscRef
            write (UnitID) Me%Var%TempViscRef
            write (UnitID) Me%Var%Viscosity
            write (UnitID) Me%Var%ViscCin
            write (UnitID) Me%Var%ViscINI
            write (UnitID) Me%Var%CVisc_E

            !Spreading
            write (UnitID) Me%Var%beta_old
            write (UnitID) Me%Var%AreaTeoric
            write (UnitID) Me%Var%OilThickness         
            write (UnitID) Me%Var%SlickThickness         
            write (UnitID) Me%Var%ThicknessLimit

            write (UnitID) Me%Var%OilSpreading

ifspre:     if (Me%Var%OilSpreading) then
            
                write (UnitID) Me%Var%SpreadingMethod
            
                if (Me%Var%SpreadingMethod .EQ. Fay_)   then        
                
                    write (UnitID) Me%Var%alfa_old
                    write (UnitID) Me%Var%DiffCoef
            
                else if (Me%Var%SpreadingMethod .EQ. ThicknessGradient_) then
            
                    write (UnitID) Me%Var%UserCoefVelMancha
            
                end if
        
            end if ifspre

            !Evaporation
            write (UnitID) Me%Var%OilEvaporation
            write (UnitID) Me%Var%MEvaporated
            write (UnitID) Me%Var%VEvaporated
            write (UnitID) Me%Var%FMEvaporated

ifevap:     if (Me%Var%OilEvaporation) then
            
                write (UnitID) Me%Var%EvaporationMethod

                if (Me%Var%EvaporationMethod .EQ. PseudoComponents) then
            
                    write (UnitID) Me%Var%NbrPC
                    write (UnitID) Me%Var%VolPC
                    write (UnitID) Me%Var%VmrelPC
                    write (UnitID) Me%Var%MWPC
                    write (UnitID) Me%Var%PvapPC
            
                else if (Me%Var%EvaporationMethod .EQ. EvaporativeExposure) then
                
                    write (UnitID) Me%Var%IBP
                    write (UnitID) Me%Var%Slope

                else if (Me%Var%EvaporationMethod .EQ. Fingas) then
                
                    write (UnitID)  Me%Var%Fingas_Evap_EqType
             
                    if (Me%Var%Fingas_Evap_EqType .EQ.Logarithmic)                      &
                        write (UnitID)  Me%Var%Time
             
                    write (UnitID) Me%Var%Fingas_Evap_Emp_Data
                
                    if (Me%Var%Fingas_Evap_Emp_Data) then
                
                        write (UnitID)  Me%Var%Fingas_Evap_Const1
                        write (UnitID)  Me%Var%Fingas_Evap_Const2
                
                    else
                
                        write (UnitID)   Me%Var%Perc_MassDist180 
                
                    end if

                end if 

            end if ifevap
         
            !Dispersion
            write (UnitID) Me%Var%OilDispersion
            write (UnitID) Me%Var%MDispersed
            write (UnitID) Me%Var%VDispersed
            write (UnitID) Me%Var%FMDispersed

ifdisp:     if (Me%Var%OilDispersion) then

                write (UnitID) Me%Var%DispersionMethod
            
                if (Me%Var%DispersionMethod .EQ. Mackay)                        &
                    write (UnitID) Me%Var%OWInterfacialTension

            end if ifdisp
         
            !Sedimentation
            write (UnitID) Me%Var%OilSedimentation
            write (UnitID) Me%Var%MSedimented
            write (UnitID) Me%Var%VSedimented
            write (UnitID) Me%Var%FMSedimented

            !Emulsification
            write (UnitID) Me%Var%OilEmulsification
            write (UnitID) Me%Var%MWaterContent
            write (UnitID) Me%Var%VWaterContent

ifemuls:    if (Me%Var%OilEmulsification) then

                write (UnitID) Me%Var%EmulsificationMethod
                write (UnitID) Me%Var%CEmuls
                write (UnitID) Me%Var%MaxVWaterContent

                if (Me%Var%EmulsificationMethod .EQ. Rasmussen) then

                    write (UnitID) Me%Var%AsphalteneContent
                    write (UnitID) Me%Var%WaxContent

                else if (Me%Var%EmulsificationMethod .EQ. Mackay) then 

                    write (UnitID) Me%Var%EmulsParameter

                end if

            end if ifemuls           
         
            !Dissolution
            write (UnitID) Me%Var%OilDissolution
            write (UnitID) Me%Var%MDissolved
            write (UnitID) Me%Var%VDissolved
            write (UnitID) Me%Var%FMDissolved

ifdiss:     if (Me%Var%OilDissolution) then
                write (UnitID) Me%Var%SolubilityOilInWater
            end if ifdiss

            !Chemical Dispersion
            write (UnitID) Me%Var%OilChemDispersion
            write (UnitID) Me%Var%MChemDispersed
            write (UnitID) Me%Var%VChemDispersed
            write (UnitID) Me%Var%FMChemDispersed

            !Mechanical Cleanup
            write (UnitID) Me%Var%OilMecCleanup
            write (UnitID) Me%Var%MOilRecovered
            write (UnitID) Me%Var%VOilRecovered
            write (UnitID) Me%Var%FMOilRecovered
            write (UnitID) Me%Var%VEmulsionRecoveryRate

            STAT_ = SUCCESS_

        else 
         
            STAT_ = ready_
        end if 


        if (present(STAT))  STAT = STAT_
       

    end subroutine WriteFinalOil
    
    !--------------------------------------------------------------------------
        
    subroutine ReadFinalOil(OilID, UnitID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: OilID
        integer                                     :: UnitID
        integer, optional                           :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_ 
        integer                                     :: STAT_ 

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(OilID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            read (UnitID) Me%Var%OilType
            read (UnitID) Me%Var%API
            read (UnitID) Me%Var%PourPoint

            read (UnitID) Me%State%FirstStepAP
            read (UnitID) Me%State%FirstStepIP
            read (UnitID) Me%Var%DTOilInternalProcesses

            !Mass
            read (UnitID) Me%Var%MassOil
            read (UnitID) Me%Var%MassINI

            !Volume
            read (UnitID) Me%Var%VolumeOil
            read (UnitID) Me%Var%VolInic

            !Density
            read (UnitID) Me%Var%Density
            read (UnitID) Me%Var%Density15

            !Viscosity
            read (UnitID) Me%Var%ViscRef
            read (UnitID) Me%Var%TempViscRef
            read (UnitID) Me%Var%Viscosity
            read (UnitID) Me%Var%ViscCin
            read (UnitID) Me%Var%ViscINI
            read (UnitID) Me%Var%CVisc_E

            !Spreading
            read (UnitID) Me%Var%beta_old
            read (UnitID) Me%Var%AreaTeoric
            read (UnitID) Me%Var%OilThickness         
            read (UnitID) Me%Var%SlickThickness         
            read (UnitID) Me%Var%ThicknessLimit

            read (UnitID) Me%Var%OilSpreading

ifspre:     if (Me%Var%OilSpreading) then
            
                read (UnitID) Me%Var%SpreadingMethod
                if (Me%Var%SpreadingMethod .EQ. Fay_)   then        
                    read (UnitID) Me%Var%alfa_old
                    read (UnitID) Me%Var%DiffCoef
                else if (Me%Var%SpreadingMethod .EQ. ThicknessGradient_) then
                    read (UnitID) Me%Var%UserCoefVelMancha
                end if
            
            end if ifspre

            !Evaporation
            read (UnitID) Me%Var%OilEvaporation
            read (UnitID) Me%Var%MEvaporated
            read (UnitID) Me%Var%VEvaporated
            read (UnitID) Me%Var%FMEvaporated

ifevap:     if (Me%Var%OilEvaporation) then
            
                read (UnitID) Me%Var%EvaporationMethod

                if (Me%Var%EvaporationMethod .EQ. PseudoComponents) then
            
                    read (UnitID) Me%Var%NbrPC
                    read (UnitID) Me%Var%VolPC
                    read (UnitID) Me%Var%VmrelPC
                    read (UnitID) Me%Var%MWPC
                    read (UnitID) Me%Var%PvapPC
            
                else if (Me%Var%EvaporationMethod .EQ. EvaporativeExposure) then
                
                    read (UnitID) Me%Var%IBP
                    read (UnitID) Me%Var%Slope

                else if (Me%Var%EvaporationMethod .EQ. Fingas) then
                
                    read (UnitID)  Me%Var%Fingas_Evap_EqType
             
                    if (Me%Var%Fingas_Evap_EqType .EQ.Logarithmic) then
                        read (UnitID)  Me%Var%Time
                    endif

                    read (UnitID) Me%Var%Fingas_Evap_Emp_Data
                
                    if (Me%Var%Fingas_Evap_Emp_Data) then
                
                        read (UnitID)  Me%Var%Fingas_Evap_Const1
                        read (UnitID)  Me%Var%Fingas_Evap_Const2
                
                    else
                
                        read (UnitID)   Me%Var%Perc_MassDist180 
                
                    end if
                
                end if 

            end if ifevap
         
            !Dispersion
            read (UnitID) Me%Var%OilDispersion
            read (UnitID) Me%Var%MDispersed
            read (UnitID) Me%Var%VDispersed
            read (UnitID) Me%Var%FMDispersed

ifdisp:     if (Me%Var%OilDispersion) then

                read (UnitID) Me%Var%DispersionMethod
            
                if (Me%Var%DispersionMethod .EQ. Mackay)                        &
                    read (UnitID) Me%Var%OWInterfacialTension

            end if ifdisp
         
            !Sedimentation
            read (UnitID) Me%Var%OilSedimentation
            read (UnitID) Me%Var%MSedimented
            read (UnitID) Me%Var%VSedimented
            read (UnitID) Me%Var%FMSedimented

            !Emulsification
            read (UnitID) Me%Var%OilEmulsification
            read (UnitID) Me%Var%MWaterContent
            read (UnitID) Me%Var%VWaterContent

ifemuls:    if (Me%Var%OilEmulsification) then

                read (UnitID) Me%Var%EmulsificationMethod
                read (UnitID) Me%Var%CEmuls
                read (UnitID) Me%Var%MaxVWaterContent

                if (Me%Var%EmulsificationMethod .EQ. Rasmussen) then

                    read (UnitID) Me%Var%AsphalteneContent
                    read (UnitID) Me%Var%WaxContent

                else if (Me%Var%EmulsificationMethod .EQ. Mackay) then 

                    read (UnitID) Me%Var%EmulsParameter
         
                end if

            end if ifemuls           
         
            !Dissolution
            read (UnitID) Me%Var%OilDissolution
            read (UnitID) Me%Var%MDissolved
            read (UnitID) Me%Var%VDissolved
            read (UnitID) Me%Var%FMDissolved

            if (Me%Var%OilDissolution) then
                read (UnitID) Me%Var%SolubilityOilInWater
            end if 

            !Chemical Dispersion
            read (UnitID) Me%Var%OilChemDispersion
            read (UnitID) Me%Var%MChemDispersed
            read (UnitID) Me%Var%VChemDispersed
            read (UnitID) Me%Var%FMChemDispersed

            !Mechanical Cleanup
            read (UnitID) Me%Var%OilMecCleanup
            read (UnitID) Me%Var%MOilRecovered
            read (UnitID) Me%Var%VOilRecovered
            read (UnitID) Me%Var%FMOilRecovered
            read (UnitID) Me%Var%VEmulsionRecoveryRate


            STAT_ = SUCCESS_

        else 
         
            STAT_ = ready_
        end if 


        if (present(STAT))  STAT = STAT_

    end subroutine ReadFinalOil

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine Ready (ObjOil_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjOil_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjOil_ID > 0) then
            call LocateObjOil (ObjOil_ID)
            ready_ = VerifyReadLock (mOil_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjOil (ObjOilID)

        !Arguments-------------------------------------------------------------
        integer                    :: ObjOilID

        !Local-----------------------------------------------------------------

        Me => FirstObjOil
        do while (associated (Me))
            if (Me%InstanceID == ObjOilID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleOil - LocateObjOil - ERR01'

    end subroutine LocateObjOil

    !--------------------------------------------------------------------------

end module ModuleOil

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------


