!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : Jet
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module which calculates the near field of an submarine outfall
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

Module ModuleJet

!BOP
!
! !MODULE: ModuleLagrangian

!    !DESCRIPTION: 
!     This model is responsible to compute the initial dilution of a submarine jet
 
! !REVISION HISTORY: 
!    
!
! !FILES USED:   
!   Do not use any file. All the files are defined by the client module (ModuleLagrangian)
!   
!              
!
! !SEE ALSO:    
!  http://www.mohid.com
!


!
!DataFile
!
!   PORT_DIAMETER           : real (m)         [0.1]      !Diameter of each port
!   PORT_BOTTOM_DISTANCE    : real (m)         [0]        !Port distance from the bottom
!   PORT_ANGLE_XY           : real (Deg.)      [0]        !Port horizontal angle
!   PORT_ANGLE_HZ           : real (Deg.)      [0]        !Port vertical angle
!   OUTFALL_LENGTH          : real (m)         [100]      !Outfall length
!   OUTFALL_ANGLE           : real (Deg.)      [0.]       !Outfall length
!   PORTS_NUMBER            : integer ()       [1]        !Number of Ports
!   RUN_MAX_PERIOD          : real (s)         [3600]     !Maximum run period
!   RUN_MIN_PERIOD          : real (s)         [60]       !Minimum run period
!   DT_OUTPUT               : real (s)         [5]        !Time interval between outputs
!   MAX_DV                  : real (%)         [1]        !Maximum volume variation between time steps
!   MAX_DT                  : real (s)         [10]       !Maximum time step interval
!   PARAMETERIZATION        : string ()        [CORJET]   !Parametrization used to simulate the 
                                                           !entrainmenet process
!   OUTPUT_TYPE             : string ()        [JET]      !The output can be made given the exact !
                                                           !information in specific output times or
                                                           !a could of particles for each output time.
!   PARTICLES_NUMBER        : integer ()       [40]       !In case of OUTPUT_TYPE = CLOUD this is the number
                                                           !of output tracer per output time interval
!   INITIAL_TRACER_CONCENTRATION : real (M/V)  [1e6]      !Initial concentration of generic tracer 
!   MAX_PLUME_DIAMETER      : real (m)         [1e3]      !Plume diameter from which initial dilution stops
                                                          !This value is used to simulate the jets overlapping
!   DISCONNECT_VERT_LIMIT   : logical          [false]    !With this option on the initial field dilution do not end
                                                          !when the jet arrives to the surface or bottom boundaries. 
                                                          !When the jet arrives to one of these boundaries a similiar
                                                          !methodology used to simulate jets overlapping is used. 

!   LOCAL_TYPE              : string ()        [FIELD3D]  !Methodology to define the ambient variables
                                                           !FIELD3D - 3D field generate by the MOHID system 
                                                           !UNIFORM - The user wants a uniform water column
                                                           !LINEAR  - The user wants a water column where 
                                                           !the density and velocity have a linear profile
!   DEFAULT_SALINITY        : real (psu)       [36]       !ambient salinity when a UNIFORM water column is admitted
!   DEFAULT_TEMPERATURE     : real (ºC)        [16]       !ambient temperature when a UNIFORM water column is admitted
!   DEFAULT_VELU            : real (m/s)       [0.2]      !ambient velocity U when a UNIFORM water column is admitted
!   DEFAULT_VELV            : real (m/s)       [0]        !ambient velocity V when a UNIFORM water column is admitted
!   BOTTOM_SALINITY         : real (psu)       [36]       !ambient bottom salinity when a LINEAR water column is admitted
!   BOTTOM_TEMPERATURE      : real (ºC)        [16]       !ambient bottom temperature when a LINEAR water column is admitted
!   BOTTOM_VELU             : real (m/s)       [0.2]      !ambient bottom velocity U when a LINEAR water column is admitted
!   BOTTOM_VELV             : real (m/s)       [0]        !ambient bottom velocity V when a LINEAR water column is admitted
!   SURFACE_SALINITY        : real (psu)       [36]       !ambient surface salinity when a LINEAR water column is admitted
!   SURFACE_TEMPERATURE     : real (ºC)        [16]       !ambient surface temperature when a LINEAR water column is admitted 
!   SURFACE_VELU            : real (m/s)       [0.2]      !ambient surface velocity U when a LINEAR water column is admitted
!   SURFACE_VELV            : real (m/s)       [0]        !ambient surface velocity V when a LINEAR water column is admitted

    use ModuleGlobalData
    use ModuleEnterData,      only : ConstructEnterData, GetData, KillEnterData, GetExtractType
    use ModuleFunctions,      only : SigmaUNESCO
    use ModuleHorizontalGrid, only : GetCellRotation    

    implicit none 

    private

!   !PUBLIC SUBROUTINES

    !Constructor
        !  Construct_Jet
    !Modifier
        !  ModifyJet
    !Selector 
        !  GetPlumeDilution
        !  GetOutPutMatrix
        !  GetPlumeLocation
        !  GetPlumeDensity
        !  GetPlumeTemperature
        !  GetPlumeSalinity
        !  GetPlumeVelocity
        !  GetPlumeMixingHorLength
        !  UngetJet
    !Destructor
        !  KillJet

    !PUBLIC TYPES: 
        !  T_Jet

 !EOP


    !subroutines---------------------------------------------------------------

    !Constructor
    public  :: Construct_Jet
    private ::      AllocateInstance
    private ::      ReadJetData
    private ::      ComputePortProperties

    !Selector
    public  :: GetPlumeDilution
    public  :: GetOutPutMatrix
    public  :: GetOutPutHeader
    public  :: GetPlumeLocation
    public  :: GetPlumeDensity
    public  :: GetPlumeTemperature
    public  :: GetPlumeSalinity
    public  :: GetPlumeVelocity
    public  :: GetPlumeMixingHorLength
    public  :: GetPlumeThickness
    public  :: GetPortGeometry
    public  :: UngetJet

    private :: UngetJetChar1D
    private :: UngetJetReal2D
    interface  UngetJet
        module procedure UngetJetChar1D
        module procedure UngetJetReal2D
    end interface UngetJet

    !Modifier
    public  :: ModifyJet
    private ::      ResetVariables
    private ::      ComputeInitialJetProperties
    private ::          ComputeInitialDt
    private ::      ComputeStationaryPlume
    private ::          ComputeOneTimeStep
    private ::              ShearEntrainmentCoef   
    private ::              ComputeShearEntrainment
    private ::              ComputeDragEntrainment
    private ::              VolumeVariation       
    private ::              DensityVariation      
    private ::              LocalAmbientProp          
    private ::              VelocityVariation         
    private ::              AuxiliarVectors           
    private ::              ThicknessVariation        
    private ::              Plumetracjectory          
    private ::              ComputeGeometry           
    private ::              TimeControl    
    private ::                  JetOutPut 
    private ::                  CloudOutPut 

    !Destructor
    public  :: KillJet
    private ::      DeallocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjJet

    !Parameter-----------------------------------------------------------------
    integer,                       parameter :: CORJET = 1, JETLAG = 2, MOHIDJET = 3, JET = 1, CLOUD = 2
    integer,                       parameter :: Field3D = 1, Uniform = 2, Linear = 3
    character(LEN = StringLength), parameter :: Char_CORJET = "CORJET", Char_JETLAG = "JETLAG",     &
                                                Char_MOHIDJET = "MOHIDJET", Char_JET = "JET",       &
                                                Char_Cloud = "CLOUD", Char_Field3D = "FIELD3D",     &
                                                Char_Uniform = "UNIFORM", Char_Linear = "LINEAR"
    real,                          parameter :: SigmaD       = 0.5
    real,                          parameter :: Precision    = 1e-4

   
    !Types---------------------------------------------------------------------
    type T_Port
        real                            :: HZAngle                      = null_real
        real                            :: XYAngle                      = null_real
        real                            :: TotalFlow                    = null_real
        real                            :: Flow                         = null_real                  
        real                            :: Diameter                     = null_real
        real                            :: x                            = null_real
        real                            :: y                            = null_real
        real                            :: z                            = null_real
        real                            :: Salinity                     = null_real
        real                            :: Temperature                  = null_real
        real                            :: Density                      = null_real
        real                            :: BottomDistance               = null_real
        real                            :: VelU                         = null_real
        real                            :: VelV                         = null_real
        real                            :: VelW                         = null_real
        real                            :: VelModulus                   = null_real
        real                            :: ex                           = null_real
        real                            :: ey                           = null_real
        real                            :: ez                           = null_real
        real                            :: OutfallLength                = null_real
        real                            :: OutfallAngle                 = null_real
        integer                         :: Number                       = null_int
        real                            :: MinPortVelocity              = null_real 
    end type T_Port

    type T_Ambient
        real, dimension(:,:,:), pointer :: Density                      => null()
        real, dimension(:,:,:), pointer :: Salinity                     => null()
        real, dimension(:,:,:), pointer :: Temperature                  => null()
        real, dimension(:,:,:), pointer :: VelU                         => null()
        real, dimension(:,:,:), pointer :: VelV                         => null()
        real, dimension(:,:,:), pointer :: VelW                         => null()
        real, dimension(:,:,:), pointer :: SZZ                          => null()
        integer                         :: BottomLayer                  = null_int
        integer                         :: SurfaceLayer                 = null_int
        integer                         :: I                            = null_int
        integer                         :: J                            = null_int
        real                            :: LocalDensity                 = null_real
        real                            :: LocalSalinity                = null_real
        real                            :: LocalTemperature             = null_real
        real                            :: LocalVelU                    = null_real
        real                            :: LocalVelV                    = null_real
        real                            :: LocalVelW                    = null_real
        real                            :: LocalVelModulus              = null_real
        integer                         :: LocalType                    = null_int
        real                            :: DefaultSalinity              = null_real 
        real                            :: DefaultTemperature           = null_real       
        real                            :: DefaultVelU                  = null_real
        real                            :: DefaultVelV                  = null_real
        real                            :: SurfaceSalinity              = null_real
        real                            :: SurfaceTemperature           = null_real 
        real                            :: SurfaceVelU                  = null_real 
        real                            :: SurfaceVelV                  = null_real 
        real                            :: BottomSalinity               = null_real 
        real                            :: BottomTemperature            = null_real 
        real                            :: BottomVelU                   = null_real 
        real                            :: BottomVelV                   = null_real 
    end type T_Ambient

    type T_NumericalOptions
        real                            :: MaxSimulationPeriod          = null_real
        real                            :: MinSimulationPeriod          = null_real
        real                            :: MaxDV                        = null_real
        real                            :: MaxDT                        = null_real
        integer                         :: Parametrization              = null_int
    end type T_NumericalOptions

    type T_OutPut
        real                                                :: DT                = null_real
        real                                                :: TimeOut           = null_real
        real, dimension(:,:), pointer                       :: Matrix            => null()
        character(LEN=StringLength), dimension(:), pointer  :: Header            => null()
        integer                                             :: Number            = null_int
        character(LEN=StringLength)                         :: AxisX             = null_str
        character(LEN=StringLength)                         :: AxisY             = null_str
        integer                                             :: FormatType        = null_int
        logical                                             :: OK                = .false.
        integer                                             :: ParticlesNumber   = null_int
        integer                                             :: OutColumns        = null_int
        real                                                :: Concentration     = null_int
        real                                                :: VelModulus_Old    = null_real
        real                                                :: Dilution_Old      = null_real
        real                                                :: Diameter_Old      = null_real
        real                                                :: Salinity_Old      = null_real
        real                                                :: Temperature_Old   = null_real
        real                                                :: Density_Old       = null_real

        real                                                :: x_old             = null_real
        real                                                :: y_old             = null_real
        real                                                :: z_old             = null_real
        real                                                :: ex_old            = null_real
        real                                                :: ey_old            = null_real
        real                                                :: ez_old            = null_real
    end type T_OutPut

    type T_Evolution
        real                            :: HZAngle                      = null_real
        real                            :: XYAngle                      = null_real
        real                            :: Volume                       = null_real
        real                            :: VolumeOld                    = null_real
        real                            :: DV                           = null_real
        real                            :: InitialVolume                = null_real
        real                            :: Dilution                     = null_real
        real                            :: Diameter                     = null_real
        real                            :: Dilution_Old                 = null_real
        real                            :: Diameter_Old                 = null_real
        real                            :: Thickness                    = null_real
        real                            :: ContactArea                  = null_real
        real                            :: VelU                         = null_real
        real                            :: VelUOld                      = null_real
        real                            :: VelV                         = null_real
        real                            :: VelVOld                      = null_real
        real                            :: VelW                         = null_real
        real                            :: VelWOld                      = null_real
        real                            :: VelModulus                   = null_real
        real                            :: DVel                         = null_real
        real                            :: x                            = null_real
        real                            :: y                            = null_real
        real                            :: Xrand                        = null_real
        real                            :: Yrand                        = null_real
        real                            :: z                            = null_real
        real                            :: x_old                        = null_real
        real                            :: y_old                        = null_real
        real                            :: z_old                        = null_real
        real                            :: Time                         = null_real
        real                            :: DT                           = null_real
        real                            :: ex                           = null_real
        real                            :: ey                           = null_real
        real                            :: ez                           = null_real
        real                            :: AmbientUt                    = null_real
        real                            :: AmbientVt                    = null_real
        real                            :: AmbientWt                    = null_real
        real                            :: AmbientModulusT              = null_real
        real                            :: AmbientUn                    = null_real
        real                            :: AmbientVn                    = null_real
        real                            :: AmbientWn                    = null_real
        real                            :: AmbientModulusN              = null_real
        real                            :: Density                      = null_real
        real                            :: DensityOld                   = null_real
        real                            :: Temperature                  = null_real
        real                            :: TemperatureOld               = null_real
        real                            :: Salinity                     = null_real
        real                            :: SalinityOld                  = null_real
        real                            :: dro                          = null_real
        real                            :: Fl2                          = null_real
        real                            :: ShearEntrainment             = null_real
        real                            :: DragEntrainment              = null_real
        real                            :: SigmaS                       = null_real
        logical                         :: EndRun                       = .false.
        logical                         :: VertBoundContact             = .false.
        logical                         :: NotFirstContact              = .false.
        logical                         :: DisconnectVertLimit          = .false.                       
        logical                         :: DragEntrainmentON            = .false.    
        logical                         :: InversionZ                   = .false.
        real                            :: MaxHorizLengthScale          = null_real
        real                            :: VertLengthScale              = null_real
        character(len=StringLength)     :: EndRunType                   = null_str
    end type T_Evolution


    type      T_Jet 
        integer                         :: InstanceID        = null_int
        type(T_Size3D )                 :: Size
        type(T_Size3D )                 :: WorkSize
        type(T_Port)                    :: Port
        type(T_Ambient)                 :: Ambient
        type(T_NumericalOptions)        :: NumericalOptions
        type(T_OutPut)                  :: OutPut
        type(T_Evolution)               :: Evolution
        integer                         :: ObjEnterData      = 0
        integer                         :: ObjHorizontalGrid = 0

        !Collection of instances
        type(T_Jet  ), pointer          :: Next

    end type T_Jet        

    !Global Module Variables
    type (T_Jet), pointer               :: FirstJet
    type (T_Jet), pointer               :: Me


    !--------------------------------------------------------------------------
    
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Construct_Jet(JetID, FileName, HorizontalGridID, PositionX, PositionY,  &
                             Flow, Salinity, Temperature, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: JetID
        character(LEN=*)                            :: FileName
        integer                                     :: HorizontalGridID
        real                                        :: PositionX, PositionY
        real                                        :: Flow, Salinity, Temperature
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_         
                       
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mJet_)) then
            nullify (FirstJet)
            call RegisterModule (mJet_) 
        endif

        call Ready(JetID, ready_)    

if0 :   if (ready_ .EQ. OFF_ERR_) then

            !Allocates a new Instance
            call AllocateInstance
            
            if (HorizontalGridID > 0) then
                Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)            
            endif

            call ConstructEnterData (Me%ObjEnterData, FileName, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'Construct_Jet - ModuleJet - ERR01'

            Me%Port%x           = PositionX
            Me%Port%y           = PositionY
            Me%Port%Salinity    = Salinity
            Me%Port%Temperature = Temperature
            Me%Port%TotalFlow   = Flow

            call ReadJetData

            call ComputePortProperties

            !Returns ID
            JetID    = Me%InstanceID
            STAT_    = SUCCESS_

        else

            stop 'Construct_Jet - ModuleJet - ERR02'

        end if if0

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine Construct_Jet

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Jet), pointer           :: NewJet
        type (T_Jet), pointer           :: PreviousJet


        !Allocates new instance
        allocate (NewJet)
        nullify  (NewJet%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstJet)) then
            FirstJet          => NewJet
            Me                => NewJet
        else
            PreviousJet       => FirstJet
            Me                => FirstJet%Next
            do while (associated(Me))
                PreviousJet   => Me
                Me            => Me%Next
            enddo
            Me                => NewJet
            PreviousJet%Next  => NewJet
        endif

        Me%InstanceID = RegisterNewInstance (mJET_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ReadJetData

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: flag
        integer                                     :: FromFile
        integer                                     :: STAT_CALL
        character(LEN=StringLength)                 :: StringAux

        !Begin-----------------------------------------------------------------

        !Gets extraction type
        call GetExtractType    (FromFile = FromFile)

        !Diameter of each port
        call GetData(Me%Port%Diameter,                                                   &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='PORT_DIAMETER',                                      &
                     ClientModule ='ModuleJet',                                          &
                     Default      = 0.1,                                                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR01'

        !Port distance from the bottom
        call GetData(Me%Port%Bottomdistance,                                             &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='PORT_BOTTOM_DISTANCE',                               &
                     ClientModule ='ModuleJet',                                          &
                     Default      = 0.,                                                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR02'

        !Port horizontal angle
        call GetData(Me%Port%XYAngle,                                                    &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='PORT_ANGLE_XY',                                      &
                     ClientModule ='ModuleJet',                                          &
                     Default      = 0.,                                                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR05'

        Me%Port%XYAngle = Me%Port%XYAngle * Pi / 180.

        !Port vertical angle
        call GetData(Me%Port%HZAngle,                                                    &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='PORT_ANGLE_HZ',                                      &
                     ClientModule ='ModuleJet',                                          &
                     Default      = 0.,                                                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR06'

        Me%Port%HZAngle = Me%Port%HZAngle * Pi / 180.

        !Number of Ports
        call GetData(Me%Port%Number,                                                     &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='PORTS_NUMBER',                                       &
                     ClientModule ='ModuleJet',                                          &
                     Default      = 1,                                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR11'

        !Maximum run period
        call GetData(Me%NumericalOptions%MaxSimulationPeriod,                            &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='RUN_MAX_PERIOD',                                     &
                     ClientModule ='ModuleJet',                                          &
                     Default      = 3600.,                                               &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR12'

        !Minimum run period
        call GetData(Me%NumericalOptions%MinSimulationPeriod,                            &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='RUN_MIN_PERIOD',                                     &
                     ClientModule ='ModuleJet',                                          &
                     Default      = 60.,                                                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR13'


        if (Me%NumericalOptions%MinSimulationPeriod >                                    &
            Me%NumericalOptions%MaxSimulationPeriod) then
            stop 'ReadJetData - ModuleJet - ERR130'
        endif

        !Time interval between outputs
        call GetData(Me%OutPut%DT,                                                       &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='DT_OUTPUT',                                          &
                     ClientModule ='ModuleJet',                                          &
                     Default      = 5.,                                                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR14'

        !Maximum volume variation between time steps [%]
        call GetData(Me%NumericalOptions%MaxDV,                                          &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='MAX_DV',                                             &
                     ClientModule ='ModuleJet',                                          &
                     Default      = 1.,                                                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR15'

        !Parametrization used to simulate the entrainment process 
        !Options (CORJET, JETLAG)
        call GetData(StringAux,                                                          &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='PARAMETERIZATION',                                   &
                     ClientModule ='ModuleJet',                                          &
                     Default      = Char_CORJET,                                         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR16'

        if (StringAux == Char_CORJET) then

            Me%NumericalOptions%Parametrization = CORJET

        else if (StringAux == Char_JETLAG) then

            Me%NumericalOptions%Parametrization = JETLAG

        else
            stop 'ReadJetData - ModuleJet - ERR17'
        endif

        !Maximum time step interval [s]
        call GetData(Me%NumericalOptions%MaxDT,                                          &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='MAX_DT',                                             &
                     ClientModule ='ModuleJet',                                          &
                     Default      = 10.,                                                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR18'

        !The output can be made given the exact information in specific 
        !output times or a could of particle for each output time.
        !The options are JET or CLOUD
        call GetData(StringAux,                                                          &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='OUTPUT_TYPE',                                        &
                     ClientModule ='ModuleJet',                                          &
                     Default      = Char_CLOUD,                                          &
                     STAT         = STAT_CALL)        

        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR19'

        if (StringAux == Char_Jet) then

            Me%OutPut%FormatType = JET

        else if (StringAux == Char_CLOUD) then

            Me%OutPut%FormatType = CLOUD

        else
            stop 'ReadJetData - ModuleJet - ERR30'
        endif

        if ( Me%OutPut%FormatType == CLOUD) then

            !If the type of output is cloud than a number of particle by output must be specified
            call GetData(Me%OutPut%ParticlesNumber,                                      &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='PARTICLES_NUMBER',                               &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 40,                                              &
                         STAT         = STAT_CALL)        

            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR20'

            call GetData(Me%OutPut%Concentration,                                        &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='INITIAL_TRACER_CONCENTRATION',                   &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 1e6,                                             &
                         STAT         = STAT_CALL)        


            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR21'

        endif

        if (Me%OutPut%FormatType == JET) then

            Me%OutPut%OutColumns = 13

            allocate(Me%OutPut%Matrix(int(Me%NumericalOptions%MaxSimulationPeriod/Me%OutPut%DT)+1,Me%OutPut%OutColumns))
            allocate(Me%OutPut%Header(Me%OutPut%OutColumns))

            Me%OutPut%Header(1)  = "Time"

            Me%OutPut%Header(2)  = "X"

            Me%OutPut%Header(3)  = "Y"

            Me%OutPut%Header(4)  = "Z"

            Me%OutPut%Header(5)  = "Dilution"

            Me%OutPut%Header(6)  = "1/Dilution"

            Me%OutPut%Header(7)  = "Length_H"

            Me%OutPut%Header(8)  = "Velocity"

            Me%OutPut%Header(9)  = "Density_Gradient"

            Me%OutPut%Header(10) = "Plume_Density"

            Me%OutPut%Header(11) = "Ambient_Density"

            Me%OutPut%Header(12) = "Thickness"

            Me%OutPut%Header(13) = "Length_V"

        else if (Me%OutPut%FormatType == CLOUD) then

            Me%OutPut%OutColumns = 16

            allocate(Me%OutPut%Matrix(int(Me%NumericalOptions%MaxSimulationPeriod/Me%OutPut%DT)* &
                                          Me%OutPut%ParticlesNumber+1,Me%OutPut%OutColumns))
            allocate(Me%OutPut%Header(Me%OutPut%OutColumns))

            Me%OutPut%Header(1)  = "Time"

            Me%OutPut%Header(2)  = "X"

            Me%OutPut%Header(3)  = "Y"

            Me%OutPut%Header(4)  = "Z"

            Me%OutPut%Header(5)  = "Dilution"

            Me%OutPut%Header(6)  = "1 / Dilution"

            Me%OutPut%Header(7)  = "Radius"

            Me%OutPut%Header(8)  = "Velocity"

            Me%OutPut%Header(9)  = "Salinity"

            Me%OutPut%Header(10) = "Temperature"

            Me%OutPut%Header(11) = "Density"

            Me%OutPut%Header(12) = "AuxRandT"

            Me%OutPut%Header(13) = "AuxRandN1"

            Me%OutPut%Header(14) = "AuxRandN2"

            Me%OutPut%Header(15) = "dr"

            Me%OutPut%Header(16) = "Length V"

        endif


        !The local ambient properties can be defined using 3D fields or 
        !considering a uniform environment or a vertical linear evolution 
        call GetData(StringAux,                                                          &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='LOCAL_TYPE',                                         &
                     ClientModule ='ModuleJet',                                          &
                     Default      = Char_Field3D,                                        &
                     STAT         = STAT_CALL)        

        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR22'

        if (StringAux == Char_Field3D) then

            Me%Ambient%LocalType = Field3D

        else if (StringAux == Char_Uniform) then

            Me%Ambient%LocalType = Uniform

        else if (StringAux == Char_Linear) then

            Me%Ambient%LocalType = Linear

        else
            stop 'ReadJetData - ModuleJet - ERR23'
        endif

   
        if (Me%Ambient%LocalType == Uniform) then

           
            call GetData(Me%Ambient%DefaultSalinity,                                     &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='DEFAULT_SALINITY',                               &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 36.,                                             &
                         STAT         = STAT_CALL)        

            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR24'


            call GetData(Me%Ambient%DefaultTemperature,                                  &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='DEFAULT_TEMPERATURE',                            &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 16.,                                             &
                         STAT         = STAT_CALL)        

            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR25'

            call GetData(Me%Ambient%DefaultVelU,                                         &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='DEFAULT_VELU',                                   &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 0.2,                                             &
                         STAT         = STAT_CALL)        

            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR26'


           call GetData(Me%Ambient%DefaultVelV,                                          &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='DEFAULT_VELV',                                   &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 0.0,                                             &
                         STAT         = STAT_CALL)        

            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR27'


        else if (Me%Ambient%LocalType == Linear) then

            call GetData(Me%Ambient%BottomSalinity,                                      &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='BOTTOM_SALINITY',                                &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 36.,                                             &
                         STAT         = STAT_CALL)        

            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR28'


            call GetData(Me%Ambient%BottomTemperature,                                   &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='BOTTOM_TEMPERATURE',                             &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 16.,                                             &
                         STAT         = STAT_CALL)        

            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR29'

            call GetData(Me%Ambient%BottomVelU,                                          &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='BOTTOM_VELU',                                    &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 0.2,                                             &
                         STAT         = STAT_CALL)        

            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR30'


           call GetData(Me%Ambient%BottomVelV,                                           &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='BOTTOM_VELV',                                    &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 0.0,                                             &
                         STAT         = STAT_CALL)        

            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR31'

            call GetData(Me%Ambient%SurfaceSalinity,                                     &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='SURFACE_SALINITY',                               &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 36.,                                             &
                         STAT         = STAT_CALL)        

            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR32'


            call GetData(Me%Ambient%SurfaceTemperature,                                  &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='SURFACE_TEMPERATURE',                            &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 16.,                                             &
                         STAT         = STAT_CALL)        

            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR33'

            call GetData(Me%Ambient%SurfaceVelU,                                         &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='SURFACE_VELU',                                   &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 0.2,                                             &
                         STAT         = STAT_CALL)        

            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR34'


           call GetData(Me%Ambient%SurfaceVelV,                                          &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromFile,                                        &
                         keyword      ='SURFACE_VELV',                                   &
                         ClientModule ='ModuleJet',                                      &
                         Default      = 0.0,                                             &
                         STAT         = STAT_CALL)        

            if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR35'

        end if 
        
        call GetData(Me%Evolution%MaxHorizLengthScale,                                   &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='MAX_PLUME_DIAMETER',                                 &
                     ClientModule ='ModuleJet',                                          &
                     Default      = 1e3,                                                 &
                     STAT         = STAT_CALL)        

        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR36'

        call GetData(Me%Port%OutfallLength,                                              &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='OUTFALL_LENGTH',                                     &
                     ClientModule ='ModuleJet',                                          &
                     Default      = 100.,                                                &
                     STAT         = STAT_CALL)        

        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR37'

        call GetData(Me%Port%OutfallAngle,                                               &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     SearchType   = FromFile,                                            &
                     keyword      ='OUTFALL_ANGLE',                                      &
                     ClientModule ='ModuleJet',                                          &
                     Default      = 0.,                                                  &
                     STAT         = STAT_CALL)        

        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR38'

        Me%Port%OutfallAngle = Me%Port%OutfallAngle * Pi / 180.


        call GetData(Me%Evolution%DisconnectVertLimit,                                  &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='DISCONNECT_VERT_LIMIT',                             &
                     ClientModule ='ModuleJet',                                         &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)        

        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR48'


        call GetData(Me%Port%MinPortVelocity,                                           &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='MIN_PORT_VELOCITY',                                 &
                     ClientModule ='ModuleJet',                                         &
                     Default      = 0.,                                                 &
                     STAT         = STAT_CALL)        

        if (STAT_CALL /= SUCCESS_) stop 'ReadJetData - ModuleJet - ERR60'
        


    end subroutine ReadJetData


    subroutine ComputePortProperties

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real :: PortArea

        !Begin-----------------------------------------------------------------

        PortArea = Pi * (Me%Port%Diameter/2)**2

        Me%Port%Flow = Me%Port%TotalFlow / Me%Port%Number

        Me%Port%VelModulus = Me%Port%Flow / PortArea

        !Defining a minimum port velocity the user can avoid unrealistic high initial dilutions
        !By default is zero but for very small flows/velocities dilutions values can be very large order of 1e5
        !If a minimum value of 0.5 m/s is considered this high dilution values canbe avoid. 
        Me%Port%VelModulus = max(Me%Port%VelModulus, Me%Port%MinPortVelocity)

        Me%Port%ex = cos(Me%Port%HZAngle) * cos(Me%Port%XYAngle)
        Me%Port%ey = cos(Me%Port%HZAngle) * sin(Me%Port%XYAngle)
        Me%Port%ez = sin(Me%Port%HZAngle)

        Me%Port%VelU = Me%Port%ex * Me%Port%VelModulus
        Me%Port%VelV = Me%Port%ey * Me%Port%VelModulus
        Me%Port%VelW = Me%Port%ez * Me%Port%VelModulus

        Me%Port%Density = SigmaUNESCO (Me%Port%Temperature, Me%Port%Salinity) + SigmaDensityReference


    end subroutine ComputePortProperties

    !--------------------------------------------------------------------------

    
    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    subroutine ModifyJet(JetID, Salinity, Temperature, VelU, VelV, VelW, SZZ,           &
                         I, J, BottomLayer, SurfaceLayer, OutPutOK, JetFlow,            &
                         JetTemperature, JetSalinity, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: JetID
        real, dimension(:,:,:), pointer             :: Salinity, Temperature,  &
                                                       VelU, VelV, VelW, SZZ
        integer                                     :: I, J, BottomLayer, SurfaceLayer
        logical                                     :: OutPutOK
        real                                        :: JetFlow, JetSalinity, JetTemperature
        integer, optional, intent(OUT)              :: STAT
   
        !Local-----------------------------------------------------------------
        integer                                     :: ready_
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

       
        STAT_ = UNKNOWN_

        call Ready(JetID, ready_)

if1 :   if (ready_ .EQ. IDLE_ERR_) then

            Me%Ambient%Salinity    => Salinity
            Me%Ambient%Temperature => Temperature
            Me%Ambient%VelU        => VelU
            Me%Ambient%VelV        => VelV
            Me%Ambient%VelW        => VelW
            Me%Ambient%SZZ         => SZZ


            Me%Ambient%I           =  I
            Me%Ambient%J           =  J
            Me%Ambient%BottomLayer =  BottomLayer
            Me%Ambient%SurfaceLayer=  SurfaceLayer

            Me%OutPut%OK           =  OutPutOK

            !GRiflet : Actualiza caudal, temperatura e salinidade
            Me%Port%Salinity    = JetSalinity
            Me%Port%Temperature = JetTemperature
            Me%Port%TotalFlow   = JetFlow
            call ComputePortProperties
            
            call ResetVariables             

            call ComputeInitialJetProperties

            call ComputeStationaryPlume     

            nullify(Me%Ambient%Salinity   )
            nullify(Me%Ambient%Temperature)
            nullify(Me%Ambient%VelU       )
            nullify(Me%Ambient%VelV       )
            nullify(Me%Ambient%VelW       )
            nullify(Me%Ambient%SZZ        )


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ModifyJet

    !--------------------------------------------------------------------------    

    subroutine ResetVariables

        !Arguments-------------------------------------------------------------
   
        !Local-----------------------------------------------------------------

        Me%Evolution%HZAngle                      = null_real
        Me%Evolution%XYAngle                      = null_real
        Me%Evolution%Volume                       = null_real
        Me%Evolution%VolumeOld                    = null_real
        Me%Evolution%DV                           = null_real
        Me%Evolution%InitialVolume                = null_real
        Me%Evolution%Dilution                     = null_real
        Me%Evolution%Diameter                     = null_real
        Me%Evolution%Dilution_Old                 = null_real
        Me%Evolution%Diameter_Old                 = null_real
        Me%Evolution%Thickness                    = null_real
        Me%Evolution%ContactArea                  = null_real
        Me%Evolution%VelU                         = null_real
        Me%Evolution%VelUOld                      = null_real
        Me%Evolution%VelV                         = null_real
        Me%Evolution%VelVOld                      = null_real
        Me%Evolution%VelW                         = null_real
        Me%Evolution%VelWOld                      = null_real
        Me%Evolution%VelModulus                   = null_real
        Me%Evolution%DVel                         = null_real
        Me%Evolution%x                            = null_real
        Me%Evolution%y                            = null_real
        Me%Evolution%z                            = null_real
        Me%Evolution%x_old                        = null_real
        Me%Evolution%y_old                        = null_real
        Me%Evolution%z_old                        = null_real
        Me%Evolution%Time                         = null_real
        Me%Evolution%DT                           = null_real
        Me%Evolution%ex                           = null_real
        Me%Evolution%ey                           = null_real
        Me%Evolution%ez                           = null_real
        Me%Evolution%AmbientUt                    = null_real
        Me%Evolution%AmbientVt                    = null_real
        Me%Evolution%AmbientWt                    = null_real
        Me%Evolution%AmbientModulusT              = null_real
        Me%Evolution%AmbientUn                    = null_real
        Me%Evolution%AmbientVn                    = null_real
        Me%Evolution%AmbientWn                    = null_real
        Me%Evolution%AmbientModulusN              = null_real
        Me%Evolution%Density                      = null_real
        Me%Evolution%DensityOld                   = null_real
        Me%Evolution%Temperature                  = null_real
        Me%Evolution%TemperatureOld               = null_real
        Me%Evolution%Salinity                     = null_real
        Me%Evolution%SalinityOld                  = null_real
        Me%Evolution%dro                          = null_real
        Me%Evolution%Fl2                          = null_real
        Me%Evolution%ShearEntrainment             = null_real
        Me%Evolution%DragEntrainment              = null_real
        Me%Evolution%SigmaS                       = null_real
        Me%Evolution%EndRun                       = .false.                       
        Me%Evolution%DragEntrainmentON            = .false.    
        Me%Evolution%InversionZ                   = .false.
        Me%Evolution%VertBoundContact             = .false.
        Me%Evolution%NotFirstContact              = .false.

        Me%OutPut%TimeOut                         = 0.
        Me%OutPut%Number                          = 0
        Me%Output%VelModulus_Old                  = null_real
        Me%Output%Dilution_Old                    = null_real
        Me%Output%Diameter_Old                    = null_real
        Me%Output%Salinity_Old                    = null_real
        Me%Output%Temperature_Old                 = null_real
        Me%Output%Density_Old                     = null_real
        Me%Output%x_old                           = null_real
        Me%Output%y_old                           = null_real
        Me%Output%z_old                           = null_real
        Me%Output%ex_old                          = null_real
        Me%Output%ey_old                          = null_real
        Me%Output%ez_old                          = null_real


    end subroutine ResetVariables
    
    !--------------------------------------------------------------------------    

    subroutine ComputeInitialJetProperties

        !Arguments-------------------------------------------------------------
   
        !Local-----------------------------------------------------------------
        real :: PortArea
        !----------------------------------------------------------------------

        Me%Evolution%HZAngle    = Me%Port%HZAngle
        Me%Evolution%XYAngle    = Me%Port%XYAngle

        Me%Evolution%ex         = Me%Port%ex 
        Me%Evolution%ey         = Me%Port%ey 
        Me%Evolution%ez         = Me%Port%ez 

        Me%Evolution%VelModulus = Me%Port%VelModulus 
        Me%Evolution%VelU       = Me%Port%VelU
        Me%Evolution%VelV       = Me%Port%VelV
        Me%Evolution%VelW       = Me%Port%VelW

        Me%Evolution%x          = Me%Port%x
        Me%Evolution%y          = Me%Port%y
        Me%Evolution%z          = Me%Ambient%SZZ(Me%Ambient%I,               &
                                                 Me%Ambient%J,               &
                                                 Me%Ambient%BottomLayer-1) - &
                                                 Me%Port%BottomDistance

        Me%Evolution%Diameter   = Me%Port%Diameter

        Me%Evolution%Salinity   = Me%Port%Salinity
        Me%Evolution%Temperature= Me%Port%Temperature
        Me%Evolution%Density    = Me%Port%Density

        Me%Evolution%InversionZ = .false.

        Call LocalAmbientProp    
        Call AuxiliarVectors     
        Call ShearEntrainmentCoef
        Call ComputeInitialDt    

        

        !Me%Evolution%Volume        = Me%Port%Flow * Me%Evolution%dt
        PortArea                   = Pi * (Me%Port%Diameter/2)**2
        Me%Evolution%Volume        = Me%Port%VelModulus * PortArea * Me%Evolution%dt
        Me%Evolution%Thickness     = Me%Port%VelModulus * Me%Evolution%dt
        Me%Evolution%InitialVolume = Me%Evolution%Volume
        Me%Evolution%ContactArea   = Pi * Me%Port%Diameter * Me%Evolution%Thickness

        Me%Evolution%Dilution      = 1

        Me%Evolution%VolumeOld     = Me%Evolution%Volume        
        Me%Evolution%DV            = 0



        Me%Output%VelModulus_Old   = Me%Evolution%VelModulus
        Me%Output%Dilution_Old     = Me%Evolution%Dilution
        Me%Output%Diameter_Old     = Me%Evolution%Diameter
        Me%Output%Salinity_Old     = Me%Evolution%Salinity
        Me%Output%Temperature_Old  = Me%Evolution%Temperature
        Me%Output%Density_Old      = Me%Evolution%Density
        Me%Output%x_old            = Me%Evolution%x
        Me%Output%y_old            = Me%Evolution%y
        Me%Output%z_old            = Me%Evolution%z
        Me%Output%ex_old           = Me%Evolution%ex
        Me%Output%ey_old           = Me%Evolution%ey
        Me%Output%ez_old           = Me%Evolution%ez


       
        !----------------------------------------------------------------------

    end subroutine ComputeInitialJetProperties

    !--------------------------------------------------------------------------    

    subroutine ComputeStationaryPlume

        !Arguments-------------------------------------------------------------
   
        !Local-----------------------------------------------------------------
        
        Me%Evolution%EndRun = .false.

        Me%Evolution%Time = 0

        do while (.not. Me%Evolution%EndRun .or.                                     &
                   Me%Evolution%Time < Me%NumericalOptions%MinSimulationPeriod)

            call ComputeOneTimeStep
            
        enddo
       
        !----------------------------------------------------------------------

    end subroutine ComputeStationaryPlume

    !--------------------------------------------------------------------------    

    subroutine ComputeOneTimeStep

        !Arguments-------------------------------------------------------------
   
        !Local-----------------------------------------------------------------
        
        !----------------------------------------------------------------------
        Call ShearEntrainmentCoef      
        Call ComputeShearEntrainment   
        If (Me%Evolution%DragEntrainmentON) Then
            Call ComputeDragEntrainment
        Else
            Me%Evolution%DragEntrainment = 0
        end If
        Call VolumeVariation          
        Call DensityVariation         
        Call LocalAmbientProp         
        Call VelocityVariation        
        Call AuxiliarVectors          
        Call ThicknessVariation       
        Call Plumetracjectory         
        Call ComputeGeometry          
        Call TimeControl              

        !----------------------------------------------------------------------

    end subroutine ComputeOneTimeStep

    !--------------------------------------------------------------------------    

    subroutine LocalAmbientProp

        !Arguments-------------------------------------------------------------
   
        !Local-----------------------------------------------------------------
        integer :: I, J, K, Kbottom, Ksurface
        real    :: a, CellRotationX, CellRotationY
        logical :: FoundLayer
        integer :: STAT_CALL
        
        !----------------------------------------------------------------------

        I        = Me%Ambient%I
        J        = Me%Ambient%J
        Kbottom  = Me%Ambient%BottomLayer - 1
        Ksurface = Me%Ambient%SurfaceLayer

        if (Me%Ambient%LocalType == Field3D) then

            do K = Me%Ambient%BottomLayer, Me%Ambient%SurfaceLayer

                if (Me%Evolution%z <= Me%Ambient%Szz(I, J, K-1) .and. Me%Evolution%z >= Me%Ambient%Szz(I, J, K)) then
                    FoundLayer = .true.
                    exit

                endif

            enddo

            if (Me%Evolution%z < Me%Ambient%Szz(I, J, Ksurface)) then
                FoundLayer = .true.
                K          = Ksurface
            endif

            if (Me%Evolution%z > Me%Ambient%Szz(I, J, Kbottom)) then
                stop 'ModuleJet - LocalAmbientProp - ERR10'
            endif

            if (FoundLayer) then

                Me%Ambient%LocalSalinity    = Me%Ambient%Salinity(I, J, K)
                Me%Ambient%LocalTemperature = Me%Ambient%Temperature(I, J, K)
                
                if (Me%ObjHorizontalGrid > 0) then
                
                    call GetCellRotation(Me%ObjHorizontalGrid, I, J,                        &
                                         CellRotationX = CellRotationX,                     &
                                         CellRotationY = CellRotationY,                     &
                                         STAT          = STAT_CALL)
                                     
                    if (STAT_CALL /= SUCCESS_) then
                        stop 'ModuleJet - LocalAmbientProp - ERR20'
                    endif          

                else
                    CellRotationX = 0.
                    CellRotationY = Pi/2.
                endif                    

                Me%Ambient%LocalVelU = Me%Ambient%VelU(I, J, K) * cos(CellRotationX) +  &
                                       Me%Ambient%VelV(I, J, K) * cos(CellRotationY)
                Me%Ambient%LocalVelV = Me%Ambient%VelU(I, J, K) * sin(CellRotationX) +  &
                                       Me%Ambient%VelV(I, J, K) * sin(CellRotationY)
                Me%Ambient%LocalVelW = Me%Ambient%VelW(I, J, K)

           
            else
                
                write(*,*) 'Z Surface =', Me%Ambient%Szz(I, J, Ksurface)    
                write(*,*) 'Z Plume   =', Me%Evolution%z
                write(*,*) 'Z Bottom  =', Me%Ambient%Szz(I, J, Kbottom)
                stop 'LocalAmbientProp - ModuleJet - ERR30'
            endif

        else if (Me%Ambient%LocalType == Uniform) then

            Me%Ambient%LocalSalinity    = Me%Ambient%DefaultSalinity
            Me%Ambient%LocalTemperature = Me%Ambient%DefaultTemperature
                                              

            Me%Ambient%LocalVelU        = Me%Ambient%DefaultVelU
            Me%Ambient%LocalVelV        = Me%Ambient%DefaultVelV
            Me%Ambient%LocalVelW        = 0


        else if (Me%Ambient%LocalType == Linear) then

            a = (Me%Ambient%Szz(I, J, Kbottom) - Me%Evolution%z)/                &
                (Me%Ambient%Szz(I, J, Kbottom) - Me%Ambient%Szz(I, J, Ksurface)) 

            Me%Ambient%LocalSalinity    = a * Me%Ambient%SurfaceSalinity     + (1-a) * Me%Ambient%BottomSalinity
            Me%Ambient%LocalTemperature = a * Me%Ambient%SurfaceTemperature  + (1-a) * Me%Ambient%BottomTemperature
                                              

            Me%Ambient%LocalVelU        = a * Me%Ambient%SurfaceVelU         + (1-a) * Me%Ambient%BottomVelU 
            Me%Ambient%LocalVelV        = a * Me%Ambient%SurfaceVelV         + (1-a) * Me%Ambient%BottomVelV
            Me%Ambient%LocalVelW        = 0


        end if 


        Me%Ambient%LocalDensity     = SigmaUNESCO (Me%Ambient%LocalTemperature, Me%Ambient%LocalSalinity) + SigmaDensityReference
        Me%Evolution%dro            = Me%Ambient%LocalDensity -                  &
                                      Me%Evolution%Density

        Me%Ambient%LocalVelModulus  = sqrt(Me%Ambient%LocalVelU**2.           + &
                                           Me%Ambient%LocalVelV**2.           + &
                                           Me%Ambient%LocalVelW**2.)

        if (Me%Ambient%LocalVelModulus>0) then

            Me%Evolution%DragEntrainmentON = .true.
        
        else

            Me%Evolution%DragEntrainmentON = .false.

        endif

    end subroutine LocalAmbientProp

    !--------------------------------------------------------------------------    

    subroutine AuxiliarVectors

        !Arguments-------------------------------------------------------------
           
        !Local-----------------------------------------------------------------
        real ::  AuxAngle1, AuxAngle2, VelHorz, Aux

        !----------------------------------------------------------------------


        
        Me%Evolution%VelModulus = Sqrt(Me%Evolution%VelU ** 2 +                  &
                                           Me%Evolution%VelV ** 2 +                  &
                                           Me%Evolution%VelW ** 2)
        VelHorz   =  Sqrt(Me%Evolution%VelU ** 2 + Me%Evolution%VelV ** 2)


        Me%Evolution%HZAngle=ComputeAngleXY (Me%Evolution%VelModulus, VelHorz, Me%Evolution%VelW)


        Me%Evolution%XYAngle=ComputeAngleXY (VelHorz, Me%Evolution%VelU, Me%Evolution%VelV)

        Me%Evolution%ex = Cos(Me%Evolution%XYAngle) * Cos(Me%Evolution%HZAngle)
        Me%Evolution%ey = Sin(Me%Evolution%XYAngle) * Cos(Me%Evolution%HZAngle)
        Me%Evolution%ez = Sin(Me%Evolution%HZAngle)

        Me%Evolution%AmbientModulusT = Me%Ambient%LocalVelU * Me%Evolution%ex &
                                         + Me%Ambient%LocalVelV * Me%Evolution%ey &
                                         + Me%Ambient%LocalVelW * Me%Evolution%ez

        Me%Evolution%DVel = Abs(Me%Evolution%VelModulus - Me%Evolution%AmbientModulusT)

        If (Me%Evolution%DragEntrainmentON) Then

            aux = Me%Evolution%AmbientModulusT / Me%Ambient%LocalVelModulus

            if (aux> 1.) aux =  1.
            if (aux<-1.) aux = -1.

            AuxAngle1 =  Acos(Aux)
            AuxAngle2 = Pi / 2 - AuxAngle1

            Me%Evolution%AmbientUn = Me%Ambient%LocalVelU * Cos(AuxAngle2)
            Me%Evolution%AmbientVn = Me%Ambient%LocalVelV * Cos(AuxAngle2)
            Me%Evolution%AmbientWn = Me%Ambient%LocalVelW * Cos(AuxAngle2)

            Me%Evolution%AmbientModulusN =  Sqrt(Me%Evolution%AmbientUn ** 2 +   &
                                                     Me%Evolution%AmbientVn ** 2 +   &
                                                     Me%Evolution%AmbientWn ** 2)
        Else
            Me%Evolution%AmbientUn       = 0
            Me%Evolution%AmbientVn       = 0
            Me%Evolution%AmbientWn       = 0
            Me%Evolution%AmbientModulusN = 0
        end If
    end subroutine

    !--------------------------------------------------------------------------    

    subroutine ShearEntrainmentCoef

        !Arguments-------------------------------------------------------------
        
   
        !Local-----------------------------------------------------------------
        real ::  AuxGravity        
        !----------------------------------------------------------------------


        ![m/s]^2 = [m/s^2*m]
        AuxGravity = Gravity * Me%Evolution%Diameter * Me%Evolution%dro / Me%Evolution%Density

        If (AuxGravity == 0) Then
            Me%Evolution%Fl2 = 100000000.0
        Else
            Me%Evolution%Fl2 = Me%Evolution%VelModulus ** 2 / AuxGravity
        end If
        Me%Evolution%Fl2 = Max(4.66 ** 2, Me%Evolution%Fl2)

        If (Me%NumericalOptions%Parametrization == CORJET) Then
            Me%Evolution%SigmaS = 0.055 + 0.6 * Sin(Me%Evolution%HZAngle) / Me%Evolution%Fl2 + &
                     0.055 * Max(Me%Evolution%AmbientModulusT / Me%Evolution%VelModulus, 0.)
            !Me%Evolution%SigmaS = 0.055 + 0.6 / 4.66 ** 2
        ElseIf (Me%NumericalOptions%Parametrization == JETLAG) Then
            Me%Evolution%SigmaS =  Sqrt(2.) * (0.057 + 0.554 * Sin(Me%Evolution%HZAngle) / Me%Evolution%Fl2) &
                        * (2 * Me%Evolution%VelModulus / (Me%Evolution%VelModulus + Me%Evolution%DVel))
        ElseIf (Me%NumericalOptions%Parametrization == MOHIDJET) Then

        end If

    end subroutine

    !--------------------------------------------------------------------------    

    subroutine ComputeInitialDt

        !Arguments-------------------------------------------------------------
           
        !Local-----------------------------------------------------------------
        
        !----------------------------------------------------------------------


        Me%Evolution%DT = Me%NumericalOptions%MaxDV / 100 * Me%Evolution%Diameter &
                              / 4 / (Me%Evolution%SigmaS * Me%Evolution%DVel + SigmaD * Me%Evolution%AmbientModulusN)


    end subroutine

    !--------------------------------------------------------------------------    

    subroutine ComputeDragEntrainment

        !Arguments-------------------------------------------------------------
           
        !Local-----------------------------------------------------------------
        
        !----------------------------------------------------------------------


        If (Me%Evolution%DragEntrainmentON) Then

            If (Me%NumericalOptions%Parametrization == CORJET) Then
                Me%Evolution%DragEntrainment = SigmaD          *    &
                Me%Evolution%AmbientModulusN * Me%Evolution%ContactArea /  Sqrt(2.)
            ElseIf (Me%NumericalOptions%Parametrization == JETLAG) Then
                Me%Evolution%DragEntrainment = Me%Evolution%AmbientModulusN *    &
                !Me%Evolution%Diameter * Me%Evolution%Thickness =
                Me%Evolution%ContactArea / Pi
            ElseIf (Me%NumericalOptions%Parametrization == MOHIDJET) Then

            end If


        Else

            Me%Evolution%DragEntrainment = 0

        end If

    end subroutine ComputeDragEntrainment

    !--------------------------------------------------------------------------    

    subroutine ComputeShearEntrainment

        !Arguments-------------------------------------------------------------
           
        !Local-----------------------------------------------------------------
        
        !----------------------------------------------------------------------


        If (Me%NumericalOptions%Parametrization == CORJET) Then
            Me%Evolution%ShearEntrainment = Me%Evolution%SigmaS *                &
            Me%Evolution%ContactArea /  Sqrt(2.) * Me%Evolution%DVel
        ElseIf (Me%NumericalOptions%Parametrization == JETLAG) Then
            Me%Evolution%ShearEntrainment = Me%Evolution%SigmaS *                &
            Me%Evolution%ContactArea * Me%Evolution%DVel
        ElseIf (Me%NumericalOptions%Parametrization == MOHIDJET) Then

        end If


    end subroutine

    !--------------------------------------------------------------------------    

    subroutine VolumeVariation

        !Arguments-------------------------------------------------------------
           
        !Local-----------------------------------------------------------------
        real                        :: Dilution
        !----------------------------------------------------------------------


        Me%Evolution%VolumeOld = Me%Evolution%Volume
        Me%Evolution%DV        = Me%Evolution%DT * (Me%Evolution%ShearEntrainment + &
                                                            Me%Evolution%DragEntrainment)
        Me%Evolution%Volume    = Me%Evolution%VolumeOld + Me%Evolution%DV

        Me%Evolution%Dilution_old = Me%Evolution%Dilution
        Dilution                    = Me%Evolution%Volume / Me%Evolution%InitialVolume
        Me%Evolution%Dilution       = Dilution

    end subroutine VolumeVariation

    !--------------------------------------------------------------------------    

    subroutine DensityVariation

        !Arguments-------------------------------------------------------------
        
   
        !Local-----------------------------------------------------------------
        
        !----------------------------------------------------------------------


        Me%Evolution%SalinityOld = Me%Evolution%Salinity
        Me%Evolution%Salinity = (Me%Evolution%SalinityOld * Me%Evolution%VolumeOld + &
                         Me%Ambient%LocalSalinity * Me%Evolution%DV) / Me%Evolution%Volume

        Me%Evolution%TemperatureOld = Me%Evolution%Temperature
        Me%Evolution%Temperature = (Me%Evolution%TemperatureOld * Me%Evolution%VolumeOld + &
                         Me%Ambient%LocalTemperature * Me%Evolution%DV) / Me%Evolution%Volume

        Me%Evolution%DensityOld = Me%Evolution%Density
        Me%Evolution%Density    = SigmaUNESCO (Me%Evolution%Temperature, Me%Evolution%Salinity) + SigmaDensityReference


    end subroutine

    !--------------------------------------------------------------------------    

    subroutine VelocityVariation
    
        !Arguments-------------------------------------------------------------
        
   
        !Local-----------------------------------------------------------------

        real ::  Impulsion, Cd, Aux
        !----------------------------------------------------------------------

        If (Me%NumericalOptions%Parametrization == CORJET) Then
            Cd = 1.3
            !Aux = Cd * Me%Evolution%Diameter * Me%Evolution%Thickness / 2  &
            Aux = Cd * Me%Evolution%ContactArea / Pi / 2.                  &
                     * Me%Ambient%LocalDensity                             &
                     * Me%Evolution%AmbientModulusN  * Me%Evolution%DT
        ElseIf (Me%NumericalOptions%Parametrization == JETLAG) Then
            Cd = 1.3
            !Aux = Cd *  Sqrt(2.) * Me%Evolution%Diameter * Me%Evolution%Thickness / 2 &
            Aux = Cd *  Sqrt(2.) * Me%Evolution%ContactArea / Pi / 2.                 &
                     *  Me%Ambient%LocalDensity *                                     &
                     Me%Evolution%AmbientModulusN  * Me%Evolution%DT
        ElseIf (Me%NumericalOptions%Parametrization == MOHIDJET) Then
            !
        end If


        Me%Evolution%VelUOld = Me%Evolution%VelU
        Me%Evolution%VelU = (Me%Evolution%VelUOld * Me%Evolution%VolumeOld * Me%Evolution%DensityOld &
                    + Me%Ambient%LocalVelU * Me%Evolution%DV * Me%Ambient%LocalDensity                   &
                    + Aux * Me%Evolution%AmbientUn) / (Me%Evolution%Volume * Me%Evolution%Density)


        Me%Evolution%VelVOld = Me%Evolution%VelV
        Me%Evolution%VelV = (Me%Evolution%VelVOld * Me%Evolution%VolumeOld * Me%Evolution%DensityOld &
                    + Me%Ambient%LocalVelV * Me%Evolution%DV * Me%Ambient%LocalDensity &
                    + Aux * Me%Evolution%AmbientVn) / (Me%Evolution%Volume * Me%Evolution%Density)

        if (.not. Me%Evolution%VertBoundContact) then

            Impulsion = Gravity * (Me%Ambient%LocalDensity - Me%Evolution%Density) / Me%Evolution%Density

            Me%Evolution%VelWOld = Me%Evolution%VelW
            Me%Evolution%VelW = (Me%Evolution%VelWOld * Me%Evolution%VolumeOld * Me%Evolution%DensityOld &
                        + Me%Ambient%LocalVelW * Me%Evolution%DV * Me%Ambient%LocalDensity &
                        + Aux * Me%Evolution%AmbientWn) / (Me%Evolution%Volume * Me%Evolution%Density) + Impulsion * Me%Evolution%DT
        endif


    end subroutine

    !--------------------------------------------------------------------------    
    
    subroutine ThicknessVariation


        !Arguments-------------------------------------------------------------
        
   
        !Local-----------------------------------------------------------------
        
        real ::  K 
        !----------------------------------------------------------------------

        K = 0

        If ( Abs(Me%Evolution%VelU) > 0) Then
            K = K + abs(Me%Evolution%ex) * (1 - Me%Evolution%VelUOld / Me%Evolution%VelU)
        end If

        If ( Abs(Me%Evolution%VelV) > 0) Then
            K = K + abs(Me%Evolution%ey) * (1 - Me%Evolution%VelVOld / Me%Evolution%VelV)
        end If

        If ( Abs(Me%Evolution%VelW) > 0) Then
            K = K + abs(Me%Evolution%ez) * (1 - Me%Evolution%VelWOld / Me%Evolution%VelW)
        end If


        If (K > 0) Then
            Me%Evolution%Thickness = Me%Evolution%Thickness * (1 + K)
        Else
            Me%Evolution%Thickness = Me%Evolution%Thickness / (1 - K)
        end If

    end subroutine

    !--------------------------------------------------------------------------    

    subroutine ComputeGeometry

        !Arguments-------------------------------------------------------------
        
   
        !Local-----------------------------------------------------------------
        real    :: RestrainedLengthScale, DensityFrontSpeed
        !----------------------------------------------------------------------


        ![m] = [m^3/m]**.5
        Me%Evolution%Diameter_old = Me%Evolution%Diameter


i1:     if (Me%Evolution%Diameter < Me%Evolution%MaxHorizLengthScale .and.              &
            .not. Me%Evolution%VertBoundContact) then

            !Disk geometry
            Me%Evolution%Diameter     =  Sqrt(Me%Evolution%Volume * 4 / Me%Evolution%Thickness / Pi)

            !plume contact area with the exterior
            ![m^2] = [m] * [m]
            Me%Evolution%ContactArea = Pi * Me%Evolution%Diameter * Me%Evolution%Thickness

        
        
        else i1
            !There is a jet overlapping or the plume arrive to the surface or bottom

i2:        if (Me%Evolution%Diameter > Me%Evolution%MaxHorizLengthScale .and. Me%Evolution%VertBoundContact) then

                ! Plume is not able to grow 
                Me%Evolution%Diameter     =  Me%Evolution%Diameter
                Me%Evolution%ContactArea  =  Me%Evolution%ContactArea

            else i2

i3:             if (Me%Evolution%Diameter > Me%Evolution%MaxHorizLengthScale) then

                    RestrainedLengthScale = Me%Evolution%MaxHorizLengthScale

                else if (Me%Evolution%VertBoundContact) then i3

                    DensityFrontSpeed = sqrt(abs(Me%Evolution%dro) / Me%Ambient%LocalDensity * &
                                              Gravity * Me%Evolution%VertLengthScale)

                    !Density front propagation 
                    Me%Evolution%VertLengthScale = Me%Evolution%VertLengthScale / (2. * DensityFrontSpeed*Me%Evolution%DT&
                                                                                  /Me%Evolution%Diameter  + 1.)

                    RestrainedLengthScale = Me%Evolution%VertLengthScale                    

                endif i3

                !A rectangular geometry is admitted 
                !In this case the variable Diameter represents the length of the not limited side of the rectangular and
                !MaxLengthScale the length of the limited side.
                !Example: 
                !   When there is jets overlapping the largest side is the interface with the other JETS 
                !   the smallest side is the interface with the "exterior"
                Me%Evolution%Diameter     =  Me%Evolution%Volume / RestrainedLengthScale / Me%Evolution%Thickness
                !plume contact area with the exterior
                ![m^2] = [m] * [m]
                Me%Evolution%ContactArea  = 2 * RestrainedLengthScale * Me%Evolution%Thickness

            endif i2

        endif i1



    end subroutine

    !--------------------------------------------------------------------------    

    subroutine Plumetracjectory

        !Arguments-------------------------------------------------------------
        
   
        !Local-----------------------------------------------------------------
        real :: dx, dy, dz, Aux, SurfaceLevel, BottomLevel, AuxRandX, AuxRandY, AuxRand
        integer :: I, J, K, Kbottom
        !----------------------------------------------------------------------

        I = Me%Ambient%I
        J = Me%Ambient%J
        K = Me%Ambient%SurfaceLayer
        KBottom = Me%Ambient%BottomLayer-1

        SurfaceLevel = Me%Ambient%Szz(I, J, K) 
        BottomLevel  = Me%Ambient%Szz(I, J, KBottom)

        dx = Me%Evolution%VelU * Me%Evolution%DT
        dy = Me%Evolution%VelV * Me%Evolution%DT
        dz = Me%Evolution%VelW * Me%Evolution%DT

        Me%Evolution%x_old = Me%Evolution%x
        Me%Evolution%y_old = Me%Evolution%y
        Me%Evolution%z_old = Me%Evolution%z

        Me%Evolution%x     = Me%Evolution%x + dx
        Me%Evolution%y     = Me%Evolution%y + dy
        Me%Evolution%z     = Me%Evolution%z - dz


i1:     if (.not. Me%Evolution%VertBoundContact) then

            If (Me%NumericalOptions%Parametrization == CORJET)  Then
                Aux = Me%Evolution%Diameter * Cos(Me%Evolution%HZangle) / 2
                If (Me%Evolution%z <= SurfaceLevel + Aux) Then
                    Me%Evolution%z = SurfaceLevel + Aux + Precision
                    Me%Evolution%VelW = 0
                    Me%Evolution%VertBoundContact = .true.
                end If


            ElseIf (Me%NumericalOptions%Parametrization == JETLAG) Then
                If (Me%Evolution%z <= SurfaceLevel) Then
                    Me%Evolution%z    = SurfaceLevel + Precision
                    Me%Evolution%VelW = 0
                    Me%Evolution%VertBoundContact = .true.
                end If

            ElseIf (Me%NumericalOptions%Parametrization == MOHIDJET) Then

            end If

            If (Me%Evolution%z >= BottomLevel) Then
                Me%Evolution%z    = BottomLevel - Precision
                Me%Evolution%VelW = 0
                Me%Evolution%VertBoundContact = .true.
            end If

        endif i1

i2:     if (Me%Evolution%VertBoundContact) then
        
            if (.not.Me%Evolution%DisconnectVertLimit) then
                Me%Evolution%EndRun = .true.
                Me%Evolution%EndRunType = "Plume Arrive to surface"
            endif

            if (.not. Me%Evolution%NotFirstContact) then
                Me%Evolution%NotFirstContact    = .true.
            endif

        else 

            Me%Evolution%VertLengthScale = Me%Evolution%Diameter

        endif i2


        if (Me%Evolution%z_old > Me%Evolution%z .and.                            &
            Me%Evolution%Time > Me%NumericalOptions%MinSimulationPeriod)         &
            Me%Evolution%InversionZ = .true.

        if (Me%Evolution%z     > Me%Evolution%z_old .and. Me%Evolution%InversionZ) then
            Me%Evolution%EndRun = .true.
            Me%Evolution%EndRunType = "Plume invert the vertical trajectory"            
        endif

        if (Me%Port%Number > 1) then
            call RANDOM_NUMBER(AuxRand)
            AuxRand = (AuxRand - 0.5) * Me%Port%OutfallLength
            AuxRandX = AuxRand * cos(Me%Port%OutfallAngle)
            AuxRandY = AuxRand * sin(Me%Port%OutfallAngle)
            !In this way the trajectory is compute for one jet located in the centered of the outfall 
            !and the position of the tracer in the far field will locate random along the outfall
            Me%Evolution%Xrand = Me%Evolution%x + AuxRandX
            Me%Evolution%Yrand = Me%Evolution%y + AuxRandY
        else
            Me%Evolution%Xrand = Me%Evolution%x
            Me%Evolution%Yrand = Me%Evolution%y
        endif

    end subroutine

    !--------------------------------------------------------------------------    

    subroutine TimeControl

        !Arguments-------------------------------------------------------------
        
   
        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------



        Me%Evolution%DT = Me%Evolution%Volume * Me%NumericalOptions%MaxDV / 100 / &
                         (Me%Evolution%ShearEntrainment + Me%Evolution%DragEntrainment)

        If (Me%Evolution%DT > Me%NumericalOptions%MaxDT) Then

            Me%Evolution%DT = Me%NumericalOptions%MaxDT

        end If

        Me%Evolution%Time = Me%Evolution%Time + Me%Evolution%DT

        if (Me%Evolution%Time > Me%NumericalOptions%MaxSimulationPeriod) then
            Me%Evolution%EndRun = .true.
            Me%Evolution%EndRunType = "Max run period"
        endif

        If (((Me%Evolution%Time > Me%OutPut%TimeOut .and. Me%Evolution%Time < Me%NumericalOptions%MaxSimulationPeriod) &
            .or. Me%Evolution%EndRun).and. Me%OutPut%Ok) Then

            if (Me%OutPut%FormatType == JET) then

                call JetOutPut

            else if (Me%OutPut%FormatType == CLOUD) then

                call CloudOutPut

            endif

            Me%OutPut%TimeOut = Me%OutPut%TimeOut + Me%OutPut%DT

        end If

    end subroutine

    Subroutine JetOutPut

        !Arguments-------------------------------------------------------------
        
   
        !Local-----------------------------------------------------------------
        integer :: I, J, Kbottom
        !----------------------------------------------------------------------

        I       = Me%Ambient%I
        J       = Me%Ambient%J
        Kbottom = Me%Ambient%BottomLayer

        Me%OutPut%Number = Me%OutPut%Number + 1

        Me%OutPut%Matrix(Me%OutPut%Number, 1) = Me%Evolution%Time

        Me%OutPut%Matrix(Me%OutPut%Number, 2) = Me%Evolution%x

        Me%OutPut%Matrix(Me%OutPut%Number, 3) = Me%Evolution%y

        Me%OutPut%Matrix(Me%OutPut%Number, 4) = Me%Evolution%z

        Me%OutPut%Matrix(Me%OutPut%Number, 5) = Me%Evolution%Dilution

        Me%OutPut%Matrix(Me%OutPut%Number, 6) = 1 / Me%Evolution%Dilution
                
        Me%OutPut%Matrix(Me%OutPut%Number, 7) = Me%Evolution%Diameter

        Me%OutPut%Matrix(Me%OutPut%Number, 8) = Me%Evolution%VelModulus

        Me%OutPut%Matrix(Me%OutPut%Number, 9) = Me%Evolution%dro

        Me%OutPut%Matrix(Me%OutPut%Number,10) = Me%Evolution%Density

        Me%OutPut%Matrix(Me%OutPut%Number,11) = Me%Ambient%LocalDensity

        Me%OutPut%Matrix(Me%OutPut%Number,12) = Me%Evolution%Thickness

        Me%OutPut%Matrix(Me%OutPut%Number, 13)= Me%Evolution%VertLengthScale


    end subroutine JetOutPut


    Subroutine CloudOutPut

        !Arguments-------------------------------------------------------------
        
   
        !Local-----------------------------------------------------------------
        real                :: dx, dy, dz, dxyz
        real                :: dxrandT, dyrandT, dzrandT, dxyzrandT, AuxRandT
        real                :: DilutionCenter, VelModulusCenter, SalinityCenter,         &
                               TemperatureCenter, DensityCenter, ConcentrationCenter
        real                :: AuxRandN1, AuxRandN2
        real                :: dxrandN1, dyrandN1, dzrandN1, dxrandN2, dyrandN2, dzrandN2
        real                :: Diameter, dr, drandN1, drandN2, dxrand, dyrand, dzrand
        real                :: zposition, SurfaceLevel
        real                :: AuxRandLengthX, AuxRandLengthY, AuxRandLength
        real, dimension(3)  :: eT, eN1, eN2
        !Local-----------------------------------------------------------------
        integer             :: I, J, Kbottom, iP
        !----------------------------------------------------------------------

        I       = Me%Ambient%I
        J       = Me%Ambient%J
        Kbottom = Me%Ambient%BottomLayer

        do ip=1, Me%OutPut%ParticlesNumber

            Me%OutPut%Number = Me%OutPut%Number + 1

            dx = Me%Evolution%x - Me%OutPut%x_old
            dy = Me%Evolution%y - Me%OutPut%y_old
            dz = Me%Evolution%z - Me%OutPut%z_old

            dxyz = sqrt(dx**2+dy**2+dz**2)

            call RANDOM_NUMBER(AuxRandT)


            dxyzrandT = dxyz * AuxRandT

            eT(1)    = Me%Evolution%ex * AuxRandT + (1.-AuxRandT) * Me%OutPut%ex_old
            eT(2)    = Me%Evolution%ey * AuxRandT + (1.-AuxRandT) * Me%OutPut%ey_old
            eT(3)    = Me%Evolution%ez * AuxRandT + (1.-AuxRandT) * Me%OutPut%ez_old

            eT(1:3)  = eT(1:3) / sqrt(eT(1)**2+eT(2)**2+eT(3)**2)

            dxrandT = eT(1) * dxyzrandT
            dyrandT = eT(2) * dxyzrandT
            dzrandT = eT(3) * dxyzrandT

            Diameter =  Me%Evolution%Diameter * AuxRandT + (1.- AuxRandT) * Me%Output%Diameter_Old

            eN1(1) = - eT(2)
            eN1(2) =   eT(1)
            eN1(3) =   0.

            eN1(1:3)  = eN1(1:3) / sqrt(eN1(1)**2+eN1(2)**2+eN1(3)**2)

            AuxRandN1 = 1.41*(2 * GaussRand() - 1)

            dxrandN1 = eN1(1) * AuxRandN1 * Diameter/2
            dyrandN1 = eN1(2) * AuxRandN1 * Diameter/2
            dzrandN1 = 0

            drandN1 = sqrt(dxrandN1**2+dyrandN1**2+dzrandN1**2)

            !Computes the cross product of eT with eN1
            call Cross_Product(eT, eN1, eN2)

            eN2(1:3)  = eN2(1:3) / sqrt(eN2(1)**2+eN2(2)**2+eN2(3)**2)

            AuxRandN2 = 1.41*(2 * GaussRand() -1)

            dxrandN2 = eN2(1) * AuxRandN2 * Diameter/2
            dyrandN2 = eN2(2) * AuxRandN2 * Diameter/2
            dzrandN2 = eN2(3) * AuxRandN2 * Diameter/2

            drandN2 = sqrt(dxrandN2**2+dyrandN2**2+dzrandN2**2)

            if (Me%Port%Number > 1) then
                call RANDOM_NUMBER(AuxRandLength)
                AuxRandLength = (AuxRandLength - 0.5) * Me%Port%OutfallLength
                AuxRandLengthX = AuxRandLength * cos(Me%Port%OutfallAngle)
                AuxRandLengthY = AuxRandLength * sin(Me%Port%OutfallAngle)
            else
                AuxRandLengthX = 0.
                AuxRandLengthY = 0.
            endif


            dxrand = dxrandT + dxrandN1 + dxrandN2 + AuxRandLengthX
            dyrand = dyrandT + dyrandN1 + dyrandN2 + AuxRandLengthY
            dzrand = dzrandT + dzrandN1 + dzrandN2

            dr = sqrt(drandN1**2 + drandN2**2)

            !linear evolution along the plume trajectory
            DilutionCenter =  Me%Evolution%Dilution     *AuxRandT + (1.-AuxRandT) * Me%Output%Dilution_Old

            ConcentrationCenter = Me%OutPut%Concentration / DilutionCenter

            VelModulusCenter = Me%Evolution%VelModulus  *AuxRandT + (1.-AuxRandT) * Me%Output%VelModulus_Old

            SalinityCenter = Me%Evolution%Salinity      *AuxRandT + (1.-AuxRandT) * Me%Output%Salinity_Old

            TemperatureCenter = Me%Evolution%Temperature*AuxRandT + (1.-AuxRandT) * Me%Output%Temperature_Old

            DensityCenter = Me%Evolution%Density        *AuxRandT + (1.-AuxRandT) * Me%Output%Density_Old

            zposition = Me%OutPut%z_old + dzrand

            SurfaceLevel = Me%Ambient%SZZ(Me%Ambient%I, Me%Ambient%J,        &
                                              Me%Ambient%SurfaceLayer)

            if (zposition < SurfaceLevel) zposition = SurfaceLevel

            !Gaussian distribution normal to the plume trajectory
            Me%OutPut%Matrix(Me%OutPut%Number, 1) = Me%Evolution%Time

            Me%OutPut%Matrix(Me%OutPut%Number, 2) = Me%OutPut%x_old + dxrand

            Me%OutPut%Matrix(Me%OutPut%Number, 3) = Me%OutPut%y_old + dyrand

            Me%OutPut%Matrix(Me%OutPut%Number, 4) = zposition

            Me%OutPut%Matrix(Me%OutPut%Number, 5) = GaussianEvolution(DilutionCenter, 0., Diameter/2., dr)

            Me%OutPut%Matrix(Me%OutPut%Number, 6) = GaussianEvolution(ConcentrationCenter, 0., Diameter/2., dr)
            
            Me%OutPut%Matrix(Me%OutPut%Number, 7) = Diameter / 2.

            Me%OutPut%Matrix(Me%OutPut%Number, 8) = GaussianEvolution(VelModulusCenter, &
                                                                      Me%Evolution%AmbientModulusT, Diameter/2., dr)

            Me%OutPut%Matrix(Me%OutPut%Number, 9) = GaussianEvolution(SalinityCenter, Me%Ambient%LocalSalinity, Diameter/2., dr)

            Me%OutPut%Matrix(Me%OutPut%Number,10) = GaussianEvolution(TemperatureCenter, &
                                                                      Me%Ambient%LocalTemperature, Diameter/2., dr)

            Me%OutPut%Matrix(Me%OutPut%Number,11) = GaussianEvolution(DensityCenter, Me%Ambient%LocalDensity, Diameter/2., dr)

            Me%OutPut%Matrix(Me%OutPut%Number,12) = AuxRandT
            
            Me%OutPut%Matrix(Me%OutPut%Number,13) = AuxRandN1

            Me%OutPut%Matrix(Me%OutPut%Number,14) = AuxRandN2

            Me%OutPut%Matrix(Me%OutPut%Number,15) = dr
        
        enddo

        Me%OutPut%TimeOut = Me%OutPut%TimeOut + Me%OutPut%DT

        Me%OutPut%VelModulus_Old    = Me%Evolution%VelModulus
        Me%OutPut%Dilution_Old      = Me%Evolution%Dilution
        Me%OutPut%Diameter_Old      = Me%Evolution%Diameter
        Me%OutPut%Salinity_Old      = Me%Evolution%Salinity
        Me%OutPut%Temperature_Old   = Me%Evolution%Temperature
        Me%OutPut%Density_Old       = Me%Evolution%Density
        Me%OutPut%x_old             = Me%Evolution%x
        Me%OutPut%y_old             = Me%Evolution%y
        Me%OutPut%z_old             = Me%Evolution%z
        Me%OutPut%ex_old            = Me%Evolution%ex
        Me%OutPut%ey_old            = Me%Evolution%ey
        Me%OutPut%ez_old            = Me%Evolution%ez




    end subroutine CloudOutPut

    !--------------------------------------------------------------------------    

    real Function ComputeAngleXY(Modulus, X, Y) 

        !Arguments-------------------------------------------------------------
        real                                        :: X, Y, Modulus, Aux

        If (Modulus > 0) Then

            aux = X / Modulus

            if (aux > 1.) aux =  1.
            if (aux <-1.) aux = -1.

            ComputeAngleXY = Acos(Aux)

            If (Y < 0) Then

                ComputeAngleXY = -ComputeAngleXY

            End If

        Else
            ComputeAngleXY = 0
        End If

    End Function

    !--------------------------------------------------------------------------    

    real Function GaussianEvolution(CenterValue, AmbientValue, Radius, r) 

        !Arguments-------------------------------------------------------------
        real                                        :: CenterValue, AmbientValue, Radius, r

        GaussianEvolution = (CenterValue - AmbientValue) * exp(-r**2/Radius**2) + AmbientValue


    End Function

    ! returns a random number between 0 and 1 with a gaussian distribution 
    real Function GaussRand() 

        !Local-------------------------------------------------------------
        real    :: Sum, Aux
        integer :: i

        Sum=0

        do i = 1, 20
            call RANDOM_NUMBER(Aux)
            Sum = Sum + Aux
        enddo

        GaussRand = Sum / 20.


    End Function
!    =================================================================
!   |                                                                |
!   |     SUBROUTINE THAT CALCULATES CROSS_PRODUCT                   |
!   |                                                                |
!   | COMMENTS : ALL LOCAL VARIABLES                                 |
!   |                                                                |
!   |                                                                |
!   |                                                                |
!    ================================================================

!>====================================================================
    
    SUBROUTINE CROSS_PRODUCT(AV,BV,CV)

!   > VARIABLES DECLARATIONS ::

!   > THE TWO VECTORS ::

!   > VECTOR A, B
        
        real, dimension(1:3), INTENT(IN)    :: AV, BV
    
    
!   > THE CROSS PRODUCT SOLUTION ::
        
        real, dimension(1:3), INTENT(OUT)   :: CV


        
!>============================================================================================================
!>============================================================================================================

!   > THE CROSS PRODUCT ROUTINE ::

!       | i j   k  |
!       |          |
!       | AX    AY  AZ | = [AY*BZ - AZ*BY]i - [AX*BZ - AZ*BX]j + [AX*BY - AY*BX]k
!       |          |
!       | BX    BY  BZ |

!   > (AX + AY + AZ ) X (BX + BY + BZ ) ::

    
!   > THE CROSS PRODUCT ::  

        CV(1) =   AV(2)*BV(3) - AV(3)*BV(2)      ! THE i direction
    
        CV(2) = -(AV(1)*BV(3) - AV(3)*BV(1))     ! the j dierction
    
        CV(3) =   AV(1)*BV(2) - AV(2)*BV(1)      ! the k direction
        
            
    END SUBROUTINE CROSS_PRODUCT




    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetPlumeDilution(JetID, PlumeDilution, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: JetID
        real                                        :: PlumeDilution
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(JetID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            PlumeDilution = Me%Evolution%Dilution

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------
     
    end subroutine GetPlumeDilution
    
    !--------------------------------------------------------------------------

    subroutine GetPlumeLocation(JetID, PlumeX, PlumeY, PlumeZ, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: JetID
        real                                        :: PlumeX, PlumeY, PlumeZ
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(JetID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            PlumeX = Me%Evolution%XRand
            PlumeY = Me%Evolution%YRand
            PlumeZ = Me%Evolution%z

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------
     
    end subroutine GetPlumeLocation

    !--------------------------------------------------------------------------

    subroutine GetPlumeDensity(JetID, PlumeDensity, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: JetID
        real                                        :: PlumeDensity
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(JetID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            PlumeDensity = Me%Evolution%Density

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------
     
    end subroutine GetPlumeDensity

    !--------------------------------------------------------------------------

    subroutine GetPlumeSalinity(JetID, PlumeSalinity, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: JetID
        real                                        :: PlumeSalinity
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(JetID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            PlumeSalinity = Me%Evolution%Salinity

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------
     
    end subroutine GetPlumeSalinity

    !--------------------------------------------------------------------------

    subroutine GetPlumeTemperature(JetID, PlumeTemperature, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: JetID
        real                                        :: PlumeTemperature
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(JetID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            PlumeTemperature = Me%Evolution%Temperature

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------
     
    end subroutine GetPlumeTemperature

    !--------------------------------------------------------------------------

    subroutine GetPlumeVelocity(JetID, PlumeVelU, PlumeVelV, PlumeVelW, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: JetID
        real                                        :: PlumeVelU, PlumeVelV, PlumeVelW
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(JetID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            PlumeVelU = Me%Evolution%VelU
            PlumeVelV = Me%Evolution%VelV
            PlumeVelW = Me%Evolution%VelW

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------
     
    end subroutine GetPlumeVelocity

    !----------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine GetPlumeMixingHorLength(JetID, PlumeHorLength, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: JetID
        real                                        :: PlumeHorLength
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(JetID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            PlumeHorLength = sqrt((Me%Evolution%x - Me%Port%x)**2 +           &
                                  (Me%Evolution%y - Me%Port%y)**2)

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------
     
    end subroutine GetPlumeMixingHorLength

    !----------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetPlumeThickness(JetID, PlumeThickness, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: JetID
        real                                        :: PlumeThickness
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(JetID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            PlumeThickness = Me%Evolution%Diameter

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------
     
    end subroutine GetPlumeThickness

    !----------------------------------------------------------------------



    subroutine GetOutPutMatrix(JetID, OutPutMatrix, OutPutLines, OutPutColumns, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: JetID
        real, dimension(:,:), pointer               :: OutPutMatrix
        integer, intent(OUT), optional              :: OutPutLines
        integer, intent(OUT), optional              :: OutPutColumns
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(JetID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            if (present(OutPutLines  )) OutPutLines   = Me%OutPut%Number
            if (present(OutPutColumns)) OutPutColumns = Me%OutPut%OutColumns
            call Read_Lock(mJET_, Me%InstanceID)
            OutPutMatrix => Me%OutPut%Matrix

  
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------
     
    end subroutine GetOutPutMatrix
    !----------------------------------------------------------------------

    !----------------------------------------------------------------------



    subroutine GetOutPutHeader(JetID, OutPutHeader, OutPutColumns, STAT)

        !Arguments-------------------------------------------------------------
        integer                                            :: JetID
        character(LEN=StringLength), dimension(:), pointer :: OutPutHeader
        integer, intent(OUT), optional                     :: OutPutColumns
        integer, optional, intent(OUT)                     :: STAT

        !Local-----------------------------------------------------------------
        integer                                            :: ready_          
        integer                                            :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(JetID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            if (present(OutPutColumns)) OutPutColumns = Me%OutPut%OutColumns
            call Read_Lock(mJET_, Me%InstanceID)
            OutPutHeader => Me%OutPut%Header

  
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------
     
    end subroutine GetOutPutHeader
    !----------------------------------------------------------------------



    subroutine GetPortGeometry(JetID, HZAngle, XYAngle, Diameter, Number, STAT)

        !Arguments-------------------------------------------------------------
        integer                                            :: JetID
        real   , optional, intent(OUT)                     :: HZAngle, XYAngle, Diameter
        integer, optional, intent(OUT)                     :: Number
        integer, optional, intent(OUT)                     :: STAT        
        !Local-----------------------------------------------------------------
        integer                                            :: ready_          
        integer                                            :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(JetID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            HZAngle     = Me%Port%HZAngle
            XYAngle     = Me%Port%XYAngle
            Diameter    = Me%Port%Diameter
            Number      = Me%Port%Number   
  
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT))                                                    &
            STAT = STAT_

    end subroutine GetPortGeometry
    
    !----------------------------------------------------------------------
    subroutine UngetJetReal2D(JetID, Matrix2D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: JetID
        real, dimension(:,:), pointer               :: Matrix2D
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(JetID, ready_) 

if1 :   if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Matrix2D)
            call Read_UnLock(mJET_, Me%InstanceID, "UngetJet")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetJetReal2D


    !----------------------------------------------------------------------
    subroutine UngetJetChar1D(JetID, Matrix1D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                   :: JetID
        Character(LEN=*), dimension(:), pointer   :: Matrix1D
        integer, optional, intent(OUT)            :: STAT

        !Local-----------------------------------------------------------------
        integer                                   :: ready_          
        integer                                   :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(JetID, ready_) 

if1 :   if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Matrix1D)
            call Read_UnLock(mJET_, Me%InstanceID, "UngetJet")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetJetChar1D



    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillJet(JetID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: JetID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_,ready_
        integer                                     :: nUsers

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(JetID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mJET_,  Me%InstanceID)

            if (nUsers == 0) then

                if (Me%ObjHorizontalGrid > 0) then
                    nUsers = DeassociateInstance(mHORIZONTALGRID_,  Me%ObjHorizontalGrid)
                    if (nUsers == 0) stop 'KillJet - ModuleJet -ERR10'
                endif    
                
                write(*,*) trim(Me%Evolution%EndRunType)

                !Deallocates the output matrix
                deallocate (Me%OutPut%Matrix)
                nullify    (Me%OutPut%Matrix)

                deallocate (Me%OutPut%Header)
                nullify    (Me%OutPut%Header)

                !Deallocates Instance
                call DeallocateInstance 


                JetID   = 0
                STAT_   = SUCCESS_

            end if

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine KillJet

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance 

        !Arguments-------------------------------------------------------------
        type (T_Jet), pointer           :: ObjJet

        !Local-----------------------------------------------------------------
        type (T_Jet), pointer           :: AuxJet
        type (T_Jet), pointer           :: PreviousJet

        !Updates pointers
        if (Me%InstanceID == FirstJet%InstanceID) then
            FirstJet => FirstJet%Next
        else
            PreviousJet       => FirstJet
            AuxJet            => FirstJet%Next
            do while (AuxJet%InstanceID /= Me%InstanceID)
                PreviousJet   => AuxJet
                AuxJet        => AuxJet%Next
            enddo

            !Now update linked list
            PreviousJet%Next  => AuxJet%Next

        endif

        !Deallocates instance
        deallocate (ObjJet)
        nullify    (ObjJet) 
            
    end subroutine DeallocateInstance


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine Ready (ObjJet_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjJet_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjJet_ID > 0) then
            call LocateObjJet(ObjJet_ID)
            ready_ = VerifyReadLock (mJet_,  Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjJet (ObjJetID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjJetID

        !Local-----------------------------------------------------------------

        Me => FirstJet
        do while (associated (Me))
            if (Me%InstanceID == ObjJetID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleJet - LocateObjJet - ERR01'

    end subroutine LocateObjJet

    !--------------------------------------------------------------------------


end Module ModuleJet

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
