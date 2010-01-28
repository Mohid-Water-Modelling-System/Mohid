!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : PorousMedia
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - Complete Revision, 
! DESCRIPTION   : Simulates Water Flow in variable saturated soils
!
!------------------------------------------------------------------------------

! Keywords read in the Data File
!
! Keyword                   : Data Type         Default     !Comment
!
! BOTTOM_FILE               : char              -           !Path to Bottom Topography File
! START_WITH_FIELD          : logical           1           !Sets Theta initial Field Capacity
! CONTINUOUS                : logical           0           !Continues from previous run
! STOP_ON_WRONG_DATE        : logical           1           !Stops if previous run end is different from actual
!                                                           !Start
! OUTPUT_TIME               : sec. sec. sec.    -           !Output Time
! SURFACE_OUTPUT_TIME       : sec. sec. sec.    -           !Output Time of surface layer
! TIME_SERIE_LOCATION       : char              -           !Path to File which defines Time Series
! CONTINUOUS_OUTPUT_FILE    : logical           1           !Writes "famous" iter.log
! CONDUTIVITYFACE           : integer           1           !Way to interpolate conducivity face
!                                                           !1 - Average, 2 - Maximum, 3 - Minimum, 4 - Weigthed
! HORIZONTAL_K_FACTOR       : real              1.0         !Factor for Horizontal Conductivity = Kh / Kv
! CUT_OFF_THETA_LOW         : real              1e-6        !Disables calculation when Theta is near ThetaR
! CUT_OFF_THETA_HIGH        : real              1e-15       !Set Theta = ThetaS when Theta > ThetaS - CUT_OFF_THETA_HIGH
! MIN_ITER                  : integer           2           !Number of iterations below which the DT is increased
! MAX_ITER                  : integer           3           !Number of iterations above which the DT is decreased
! LIMIT_ITER                : integer           50          !Number of iterations of a time step (for restart)
! THETA_TOLERANCE           : real              0.001       !Converge Parameter
! INCREASE_DT               : real              1.25        !Increase of DT when iter < MIN_ITER
! DECREASE_DT               : real              0.70        !Decrease of DT when iter > MAX_ITER
!
!
!
!<beginproperty>
! NAME                      : Theta / waterlevel 
!
! see Module FillMatrix for more options
!
!<endproperty>

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


Module ModulePorousMedia

    use ModuleGlobalData
    use ModuleStopWatch
    use ModuleFunctions
    use ModuleTime
    use ModuleHDF5
    use ModuleEnterData
    use ModuleProfile,          only : StartProfile, WriteProfile, KillProfile
    use ModuleGridData,         only : ConstructGridData, GetGridData, UngetGridData,    &
                                       KillGridData           
    use ModuleTimeSerie,        only : StartTimeSerie, WriteTimeSerie, KillTimeSerie         
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, GetGridCellArea,               &
                                       WriteHorizontalGrid, UnGetHorizontalGrid
    use ModuleBasinGeometry,    only : GetBasinPoints, GetRiverPoints, GetCellSlope,     &
                                       UnGetBasin
    use ModuleGeometry,         only : ConstructGeometry, GetGeometrySize,               &
                                       GetGeometryDistances, GetGeometryKFloor,          &
                                       UnGetGeometry, ComputeInitialGeometry,            &
                                       ComputeVerticalGeometry, GetGeometryVolumes,      &
                                       GetGeometryAreas, KillGeometry 
    use ModuleMap,              only : ConstructMap, GetWaterPoints3D, GetOpenPoints3D,  &
                                       GetComputeFaces3D, UnGetMap,                      &
                                       UpdateComputeFaces3D, KillMap         
    use ModuleFillMatrix,       only : ConstructFillMatrix, KillFillMatrix
    use ModuleDrainageNetwork,  only : GetChannelsWaterLevel, GetChannelsBottomLevel,    &
                                       GetChannelsBottomWidth, GetChannelsOpenProcess,   &
                                       GetChannelsNodeLength, UnGetDrainageNetwork


    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  ::  ConstructPorousMedia
    private ::      AllocateInstance
    private ::      ReadDataFile
    private ::      ConstructBottomTopography
    private ::      AllocateVariables
    private ::      ReadSoilTypes
    private ::      InitialFields
    private ::      ReadInitialSoilFile
    private ::      ConstructHDF5Output    
    private ::      ConstructTimeSerie

    !Selector
    public  ::  GetNextPorousMediaDT
    public  ::  GetPotentialInfiltration
    public  ::  GetInfiltration
    public  ::  GetEfectiveEVTP
    public  ::  GetGWFlowToChannels
    public  ::  GetGWLayer
    public  ::  GetTotalStoredVolume
    public  ::  GetFluxU
    public  ::  GetFluxV
    public  ::  GetFluxW
    public  ::  GetUnsatU
    public  ::  GetUnsatV
    public  ::  GetUnsatW
    public  ::  GetWaterColumn
    public  ::  GetWaterContent
    public  ::  GetHead
    public  ::  GetThetaR
    public  ::  GetThetaS
    public  ::  GetOldWaterContent
    public  ::  GetThetaField
    public  ::  GetComputeSoilField
    public  ::  GetLimitThetaLow
    public  ::  GetUnsatK
    public  ::  GetEvaporation
    public  ::  UnGetPorousMedia
    
    !Modifier
    public  ::  ModifyPorousMedia      
    private ::      VariableSaturatedFlow
    private ::          Condutivity_Face
    private ::              CondutivityAverage
    private ::              CondutivityMaximum
    private ::              CondutivityMinimum
    private ::              CondutivityGeometricAverage
    private ::          SoilWaterVelocity
    private ::          SoilParameters
    private ::          VerticalContinuity
    private ::          variation_test
    private ::      ExchangeWithDrainageNetwork
    private ::      PorousMediaOutput
    private ::          OutPutTimeSeries
    private ::      CalculateTotalStoredVolume

    !Destructor
    public  ::  KillPorousMedia                                                     
    private ::      DeAllocateInstance
    private ::      WriteFinalSoilFile

    !Management
    private ::      Ready
    private ::          LocateObjPorousMedia 
    
    !Interfaces----------------------------------------------------------------

    interface  UnGetPorousMedia
        module procedure UnGetPorousMedia_R4
        module procedure UnGetPorousMedia_R8
        module procedure UnGetPorousMedia_R
        module procedure UnGetPorousMedia_R1
        module procedure UnGetPorousMedia_RI
        module procedure UnGetPorousMedia_AI
        module procedure UnGetPorousMedia_AI2D
    end interface UnGetPorousMedia

    !Parameter-------------------------------------------------------------------

    !Conductivity Face
    integer, parameter :: Average           = 1
    integer, parameter :: Maximum           = 2
    integer, parameter :: Minimum           = 3
    integer, parameter :: Weighted          = 4
    integer, parameter :: GeometricAvg      = 5

    !Evapotranspiration Method
    integer, parameter :: SingleEvapotranspiration   = 1
    integer, parameter :: SeparateEvapotranspiration = 2
    !Min Thickness UG Watercolumn
    real, parameter    :: MinUGThickness    = 0.10
    
    !Water Contents
    character(LEN = StringLength), parameter :: char_Theta        = trim(adjustl('Theta'            ))

    !Conductivity
    character(LEN = StringLength), parameter :: char_SoilID       = trim(adjustl('SoilID'))                        

    !Waterlevel
    character(LEN = StringLength), parameter :: char_waterlevel   = trim(adjustl('waterlevel'       ))

    !Types---------------------------------------------------------------------
    type T_OutPut
        type (T_Time), pointer, dimension(:)    :: OutTime
        type (T_Time), dimension(:), pointer    :: RestartOutTime
        type (T_Time), dimension(:), pointer    :: SurfaceOutTime
        integer                                 :: NextOutPut
        logical                                 :: Yes = .false.
        logical                                 :: TimeSerieON
        logical                                 :: ProfileON
        logical                                 :: WriteRestartFile     = .false.
        logical                                 :: SurfaceOutput        = .false.
        logical                                 :: RestartOverwrite     = .false.
        integer                                 :: NextRestartOutput    = 1
        integer                                 :: NextSurfaceOutput    = 1
    end type T_OutPut

    type T_Files
        character(PathLength)                   :: DataFile
        character(PathLength)                   :: InitialFile
        character(PathLength)                   :: FinalFile
        character(PathLength)                   :: TransientHDF
        character(PathLength)                   :: BottomFile
        integer                                 :: AsciiUnit
    end type T_Files    


    type T_ExtVar
        integer, dimension(:,:), pointer        :: BasinPoints
!        real   , dimension(:,:), pointer        :: PermeableFraction
        
        !ObjGeometry
        real   , pointer, dimension(:,:  )      :: DUX, DVY
        real   , pointer, dimension(:,:  )      :: DZX, DZY
        real   , pointer, dimension(:,:,:)      :: DZZ, DWZ
        real   , pointer, dimension(:,:,:)      :: CenterCell
        real   , pointer, dimension(:,:  )      :: Area
        
        real   , pointer, dimension(:,:  )      :: Topography  
        real   , pointer, dimension(:,:  )      :: BottomTopoG               
        
        real   , pointer, dimension(:,:,:)      :: AreaU, AreaV

        integer, dimension(:,:), pointer        :: RiverPoints
        integer, dimension(:,:), pointer        :: KFloor      

        real(8), pointer, dimension(:,:,:)      :: CellVolume
        real ,   pointer, dimension(:,:,:)      :: SZZ
                
        
        !Map 
        integer, pointer, dimension(:,:,:)      :: WaterPoints3D
        integer, pointer, dimension(:,:,:)      :: OpenPoints3D
        integer, pointer, dimension(:,:,:)      :: ComputeFacesU3D
        integer, pointer, dimension(:,:,:)      :: ComputeFacesV3D
        integer, pointer, dimension(:,:,:)      :: ComputeFacesW3D      
        
        !Vegetation Stuff
        real, dimension(:,:), pointer           :: RootDepth
        real, dimension(:,:), pointer           :: RootFeddesH1
        real, dimension(:,:), pointer           :: RootFeddesH2
        real, dimension(:,:), pointer           :: RootFeddesH3
        real, dimension(:,:), pointer           :: RootFeddesH4
        
        real(8), dimension(:,:  ), pointer      :: InfiltrationColumn
        real, dimension(:,:,:), pointer         :: TranspirationFlux        
        real, dimension(:,:  ), pointer         :: PotentialEvaporationFlux         
        logical                                 :: ConstructEvaporation
        logical                                 :: ConstructTranspiration   

        !Time
        type (T_Time)                           :: Now
        real                                    :: DT
    end type T_ExtVar
    
    !Unsaturated Zone Types
    type T_SoilOptions
        logical :: CalcHorizontal
        integer :: CondutivityFace
        logical :: Continuous
        logical :: StopOnWrongDate
        logical :: CheckGlobalMass
        logical :: StartWithFieldCapacity
        logical :: ComputeSoilField
        logical :: RemoveWater              !Lúcia
        real    :: HCondFactor
        logical :: LimitEVAPWaterVelocity
        logical :: LimitEVAPHead
        real    :: HeadLimit
    end type T_SoilOptions

    type T_SoilType
        real                                :: ThetaR           = null_real
        real                                :: ThetaS           = null_real
        real                                :: nfit             = null_real
        real                                :: mfit             = null_real
        real                                :: alfa             = null_real
        real                                :: lfit             = null_real
        real                                :: SatK             = null_real
        real                                :: OverSatSlope     = null_real
    end type T_SoilType

    type T_Retention !Main parameters in the Mualem-van Genuchten retention and conductivity cuves
        real, dimension(:,:,:), pointer     :: ThetaR           => null()     !Minimum water content
        real, dimension(:,:,:), pointer     :: ThetaS           => null()     !Saturated water content
        real, dimension(:,:,:), pointer     :: ThetaF           => null()     !(Theta-ThetaR)/(ThetaS-ThetaR)
    end type T_Retention

    type T_Converge
        integer      :: MinIter             = null_int
        integer      :: MaxIter             = null_int
        integer      :: IterStep            = null_int
        integer      :: LimitIter           = null_int
        real         :: ThetaTolerance      = null_real
        real         :: IncreaseDT          = null_real
        real         :: DecreaseDT          = null_real
        real         :: PredictedDT         = null_real 
        real         :: CurrentDT           = null_real 
        real         :: MinDT               = null_real
        real         :: LimitThetaLo
        real         :: LimitThetaHi
        real         :: ThetaHydroCoef      = null_real

        real,    pointer, dimension(:,:,:) :: ThetaOld
        real,    pointer, dimension(:,:,:) :: ThetaIni
!        real,    pointer, dimension(:,:,:) :: HeadOld
        real,    pointer, dimension(:,:,:) :: HeadIni
    end type T_Converge
    
    type T_Property
        type (T_PropertyID)                     :: ID
        real, dimension(:,:,:), pointer         :: Concentration            => null()
        real, dimension(:,:,:), pointer         :: ConcentrationIni         => null()
        real, dimension(:,:,:), pointer         :: ConcentrationOld         => null()
        real, dimension(:,:  ), pointer         :: UpperConcentration       => null()
        type (T_Property), pointer              :: Next                     => null()
    end type T_Property

    type       T_PorousMedia        
        !Instaces of other Objects
        integer                                 :: InstanceID
        character(len=StringLength)             :: ModelName
        integer                                 :: ObjBasinGeometry         = 0
        integer                                 :: ObjTime                  = 0
        integer                                 :: ObjGeometry              = 0
        integer                                 :: ObjMap                   = 0
        integer                                 :: ObjHorizontalGrid        = 0
        integer                                 :: ObjHorizontalMap         = 0
        integer                                 :: ObjTopography            = 0
        integer                                 :: ObjTimeSerie             = 0
        integer                                 :: ObjHDF5                  = 0
        integer                                 :: ObjDrainageNetwork       = 0
        integer                                 :: ObjBottomTopography      = 0
        integer                                 :: ObjEnterData             = 0
        integer                                 :: ObjProfile               = 0
        real,    pointer, dimension(:,:,:)      :: ThetaField                => null() !!FieldCapacity [m3/m3]                

        type (T_OutPut)                         :: OutPut
        type (T_ExtVar)                         :: ExtVar
        type (T_Files)                          :: Files
        type (T_Time)                           :: BeginTime
        type (T_Time)                           :: EndTime
        real(8),    pointer, dimension(:,:  )   :: WaterColumn              => null()
        real(8),    pointer, dimension(:,:  )   :: Infiltration             => null()
        real(8),    pointer, dimension(:,:  )   :: EfectiveEVTP             => null()

        !Watertable Properties
        real,    dimension(:,:), pointer        :: UGWaterLevel2D           => null()
        real,    dimension(:,:), pointer        :: UGWaterDepth2D           => null()
        integer, dimension(:,:), pointer        :: UGCell                   => null()
        
        !Exchange with channels
        real,    dimension(:,:), pointer        :: lFlowToChannels          => null()
        real,    dimension(:,:), pointer        :: iFlowToChannels          => null()

        !Velocities
        real,    dimension(:,:,:), pointer      :: UnsatVelU                => null()
        real,    dimension(:,:,:), pointer      :: UnsatVelV                => null()
        real,    dimension(:,:,:), pointer      :: UnsatVelW                => null()
        real,    dimension(:,:,:), pointer      :: UnsatVelWOld             => null()

        !infiltration 
        real,   pointer, dimension(:,:)         :: InfiltrationVelocity     => null()
        real,   pointer, dimension(:,:)         :: ImpermeableFraction

        !Fluxes
        real(8), dimension(:,:,:), pointer      :: FluxU                    => null()
        real(8), dimension(:,:,:), pointer      :: FluxV                    => null()
        real(8), dimension(:,:,:), pointer      :: FluxW                    => null()
        real(8), dimension(:,:,:), pointer      :: FluxWFinal               => null()  !Flux Corrected with the routine Vertical continuity
        real,    dimension(:,:  ), pointer      :: EvaporationFlux          => null()
        !Flow Properties
        real,    pointer, dimension(:,:,:)      :: Theta                    => null() !water content on each cell [m3/m3]                
        real,    pointer, dimension(:,:,:)      :: Head                     => null() !Suction Head on each cell 
        real,    pointer, dimension(:,:,:)      :: HydroPressure            => null() !Hydrostatic pressure
        real,    pointer, dimension(:,:,:)      :: FinalHead                => null() !Sum of Suction, Hydrostatic and Topography

        !Common Properties
        real,    pointer, dimension(:,:,:)      :: SatK                     => null()
        integer, pointer, dimension(:,:,:)      :: SoilID                   => null()
        real,    pointer, dimension(:,:,:)      :: UnSatK                   => null()
        real,    pointer, dimension(:,:,:)      :: UnSatK_X                 => null()
        real,    pointer, dimension(:,:,:)      :: UnSatK_Y                 => null()
        real,    pointer, dimension(:,:,:)      :: UnSatK_Z                 => null()

        !Auxiliar SpeedUp Matrixes          
        logical, pointer, dimension(:,:,:)      :: CalculateHead            => null()
    
        logical                                 :: TranspirationExists            = .false.
        logical                                 :: EvaporationExists              = .false.

        type (T_Retention       )               :: RC           !retention curve
        type (T_Converge        )               :: CV           !Converge data 

        !Options
        real                                    :: NextDT
        !Unsaturated Options
        type (T_SoilOptions )                   :: SoilOpt

        real(8)                                 :: TotalStoredVolume    = 0.0
        real(8)                                 :: LossToGround         = 0.0

        !Grid size
        type (T_Size3D)                         :: Size,   WorkSize       
        type (T_Size2D)                         :: Size2D
        
        !Soil Types
        type (T_SoilType), dimension(:), pointer :: SoilTypes       => null()
        
        !Properties
        type (T_Property), pointer              :: FirstProperty    => null()
                
        type(T_PorousMedia), pointer            :: Next             => null()
        
        !Accumulated fluxes for Average flux computation (for Advection diffusion) 
        real(8), dimension(:,:,:), pointer      :: FluxUAcc       => null()
        real(8), dimension(:,:,:), pointer      :: FluxVAcc       => null()
        real(8), dimension(:,:,:), pointer      :: FluxWAcc       => null()
        real(8), dimension(:,:,:), pointer      :: FluxWAccFinal  => null()        
    end type  T_PorousMedia

    !Global Module Variables
    type (T_PorousMedia), pointer               :: FirstObjPorousMedia
    type (T_PorousMedia), pointer               :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructPorousMedia(ModelName,                                  &
                                    ObjPorousMediaID,                           &
                                    ComputeTimeID,                              &
                                    HorizontalGridID,                           &
                                    HorizontalMapID,                            &
                                    TopographyID,                               &
                                    BasinGeometryID,                            &
                                    DrainageNetworkID,                          &
                                    CheckGlobalMass,                            &
                                    ConstructEvaporation,                       &
                                    ConstructTranspiration,                     &
                                    GeometryID,                                 &
                                    MapID,                                      &
                                    STAT)

        !Arguments---------------------------------------------------------------
        character(len=*)                                :: ModelName
        integer                                         :: ObjPorousMediaID 
        integer                                         :: ComputeTimeID
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: TopographyID
        integer                                         :: BasinGeometryID
        integer                                         :: DrainageNetworkID
        logical                                         :: CheckGlobalMass
        logical                                         :: ConstructEvaporation
        logical                                         :: ConstructTranspiration
        integer, intent (OUT)                           :: GeometryID
        integer, intent (OUT)                           :: MapID
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, STAT_CALL

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_


        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mPorousMedia_)) then
            
            nullify (FirstObjPorousMedia)
            call RegisterModule (mPorousMedia_) 
        
        endif

        call Ready(ObjPorousMediaID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%ModelName = ModelName

            !Associate External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           ComputeTimeID   )
            Me%ObjTopography     = AssociateInstance (mGRIDDATA_,       TopographyID    ) 
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjBasinGeometry  = AssociateInstance (mBASINGEOMETRY_,  BasinGeometryID )

            if (DrainageNetworkID /= 0)                                                     &
                Me%ObjDrainageNetwork  = AssociateInstance (MDRAINAGENETWORK_,  DrainageNetworkID )
                
            Me%SoilOpt%CheckGlobalMass       = CheckGlobalMass
            Me%ExtVar%ConstructEvaporation   = ConstructEvaporation
            Me%ExtVar%ConstructTranspiration = ConstructTranspiration

            !Time
            call GetComputeCurrentTime  (Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMedia - ModulePorousMedia - ERR01'
                    
            call GetComputeTimeLimits   (Me%ObjTime, BeginTime = Me%BeginTime,              &
                                         EndTime = Me%EndTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMedia - ModulePorousMedia - ERR01a'            

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMedia - ModulePorousMedia - ERR01b'            

           
            !Read data files
            call ReadDataFile
                        
            call ConstructBottomTopography

            !Build 3D domain
            call VerticalDiscretization
            
            !After 3D build send to ModuleBasin IDs to be associated with ModulePorousMediaProperties and 
            !ModuleVegetation
            GeometryID = Me%ObjGeometry
            MapID      = Me%ObjMap

            call ReadLockExternalVar 
           
            call AllocateVariables
            
            call ReadSoilTypes
            
            !Build Initial fields
            call InitialFields

            if (Me%SoilOpt%Continuous)  call ReadInitialSoilFile
            
            !Calculates initial Theta
            if (.not. Me%SoilOpt%Continuous) then
                call StartWithFieldCapacity
            endif    
            
            !vegetation model growth needs field capacity computation
            if (Me%SoilOpt%ComputeSoilField) then
                call ComputeSoilFieldCapacity
            endif   

            !Set initial time steps
            Me%CV%CurrentDT   = Me%ExtVar%DT
            Me%CV%PredictedDT = Me%ExtVar%DT

            !Calculate initial heads
            call SoilParameters (Me%ExtVar%BasinPoints)
            
            !Calculates Initial GW Cell
            call CalculateUGWaterLevel
            
            if (Me%OutPut%Yes .or. Me%Output%SurfaceOutput) then
                call ConstructHDF5Output
            endif

            call ConstructTimeSerie

            call ConstructProfileOutput
            
            call ConstructASCIIOutput

            call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMedia - ModulePorousMedia - ERR02'

            if (Me%SoilOpt%CheckGlobalMass) then
                call CalculateTotalStoredVolume
            endif

            call ReadUnLockExternalVar

            !Returns ID
            ObjPorousMediaID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModulePorousMedia - ConstructPorousMedia - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructPorousMedia
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance        
                                                    
        !Local-----------------------------------------------------------------
        type (T_PorousMedia), pointer       :: NewObjPorousMedia
        type (T_PorousMedia), pointer       :: PreviousObjPorousMedia


        !Allocates new instance
        allocate (NewObjPorousMedia)
        nullify  (NewObjPorousMedia%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjPorousMedia)) then
           
            FirstObjPorousMedia         => NewObjPorousMedia
            Me                          => NewObjPorousMedia
        
        else
        
            PreviousObjPorousMedia      => FirstObjPorousMedia
            Me                          => FirstObjPorousMedia%Next
            do while (associated(Me))
                PreviousObjPorousMedia  => Me
                Me                      => Me%Next
            enddo
            Me                          => NewObjPorousMedia
            PreviousObjPorousMedia%Next => NewObjPorousMedia
        
        endif

        Me%InstanceID = RegisterNewInstance (mPorousMedia_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------
    
    subroutine ReadDataFile        

        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Reads the name of the data file from nomfich
        call ReadFileName ('POROUS_DATA', Me%Files%DataFile, "PorousMedia Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR01'

        !Reads the name of the transient HDF file from nomfich
        call ReadFileName ('POROUS_HDF', Me%Files%TransientHDF, "PorousMedia HDF File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR01b'
                
        !Reads the name of the file where to store final data
        call ReadFileName ('POROUS_FIN', Me%Files%FinalFile, "PorousMedia Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR01c'

        !Constructs the DataFile
        call ConstructEnterData (Me%ObjEnterData, Me%Files%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR02'

        !Botom file
        call GetData(Me%Files%BottomFile,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BOTTOM_FILE',                                      &
                     ClientModule = 'ModulePorousMedia',                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR02a'

        !General Options        
        call GetData(Me%SoilOpt%StartWithFieldCapacity,                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType = FromFile,                                             &
                     keyword    = 'START_WITH_FIELD',                                   &
                     Default    = .true.,                                               &                                           
                     ClientModule ='ModulePorousMedia',                                 &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR03'

        call GetData(Me%SoilOpt%ComputeSoilField,                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType = FromFile,                                             &
                     keyword    = 'COMPUTE_SOIL_FIELD',                                 &
                     Default    = .false.,                                              &                                           
                     ClientModule ='ModulePorousMedia',                                 &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR04'

        call GetData(Me%SoilOpt%RemoveWater,                                            &   
                     Me%ObjEnterData, iflag,                                            &
                     SearchType = FromFile,                                             &
                     keyword    = 'REMOVE_WATER',                                       &
                     Default    = .true.,                                               &                                           
                     ClientModule ='ModulePorousMedia',                                 &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR05'
   
        call GetData(Me%SoilOpt%LimitEVAPHead,                                          &     
                     Me%ObjEnterData, iflag,                                            &
                     SearchType     = FromFile,                                         &
                     keyword        ='LIMIT_EVAP_HEAD',                                 &
                     Default        = .false.,                                          &
                     ClientModule   ='ModulePorousMedia',                               &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'GetUnSaturatedOptions - ModulePorousMedia - ERR06'

        if (Me%SoilOpt%LimitEVAPHead) then
                
            call GetData(Me%SoilOpt%HeadLimit,                                         &     
                         Me%ObjEnterData, iflag,                                       &
                         SearchType     = FromFile,                                    &
                         keyword        ='HEAD_LIMIT',                                 &
                         Default        = -100.0,                                      &
                         ClientModule   ='ModulePorousMedia',                          &
                         STAT           = STAT_CALL)             
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR6a'
        
        endif

        call GetData(Me%SoilOpt%Continuous,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType = FromFile,                                             &
                     keyword    = 'CONTINUOUS',                                         &
                     Default    = .false.,                                              &                                           
                     ClientModule ='ModulePorousMedia',                                 &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR07'
        
        if (Me%SoilOpt%Continuous) then
            !Reads the name of the file where to read initial data
            call ReadFileName ('POROUS_INI', Me%Files%InitialFile, "PorousMedia Initial File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR008'

            call GetData(Me%SoilOpt%StopOnWrongDate,                                    &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType = FromFile,                                         &
                         keyword    = 'STOP_ON_WRONG_DATE',                             &
                         Default    = .true.,                                           &                                           
                         ClientModule ='ModulePorousMedia',                             &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR09'

        endif


        !Output Options--------------------------------------------------------

        !Gets Output Time 
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime = Me%ExtVar%Now,                                 &
                           EndTime     = Me%EndTime,                                    &
                           keyword     = 'OUTPUT_TIME',                                 &
                           SearchType  = FromFile,                                      &
                           OutPutsTime = Me%OutPut%OutTime,                             &
                           OutPutsOn   = Me%OutPut%Yes,                                 &
                           STAT        = STAT_CALL)
        Me%OutPut%NextOutPut = 1
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR010'

        !Output for restart
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime  = Me%ExtVar%Now,                                &
                           EndTime      = Me%EndTime,                                   &
                           keyword      = 'RESTART_FILE_OUTPUT_TIME',                   &
                           SearchType   = FromFile,                                     &
                           OutPutsTime  = Me%OutPut%RestartOutTime,                     &
                           OutPutsOn    = Me%OutPut%WriteRestartFile,                   &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR11a'

        call GetData(Me%OutPut%RestartOverwrite,                                        &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESTART_FILE_OVERWRITE',                           &
                     Default      = .true.,                                             &
                     ClientModule = 'ModulePorousMedia',                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleBasin - ERR11c'

        !Output for surface output
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime  = Me%ExtVar%Now,                                &
                           EndTime      = Me%EndTime,                                   &
                           keyword      = 'SURFACE_OUTPUT_TIME',                        &
                           SearchType   = FromFile,                                     &
                           OutPutsTime  = Me%OutPut%SurfaceOutTime,                     &
                           OutPutsOn    = Me%OutPut%SurfaceOutput,                      &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR11b'

        !Checks consistency
        if (Me%OutPut%Yes .and. Me%OutPut%SurfaceOutput) then
            write(*,*)'Only normal output or 2D output can be active'
            write(*,*)'OUTPUT_TIME or SURFACE_OUTPUT_TIME'
            stop 'ReadDataFile - ModuleBasin - ERR11d'
        endif

        !Directional Options---------------------------------------------------        

        call GetData(Me%SoilOpt%CalcHorizontal,                                 &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='CALC_HORIZONTAL',                         &
                     Default        =.TRUE.,                                    &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'GetUnSaturatedOptions - ModulePorousMedia - ERR02'

       
        ! 1 - AVERAGE of the conductivity in the cells
        ! 2 - MAXIMUM of the conductivity in the cells
        ! 3 - MINIMUM of the conductivity in the cells
        ! 4 - WEIGTHED of the conductivity in the cells
        ! 5 - GEOMETRIC AVERAGE of the conductivity in the cells

        call GetData(Me%SoilOpt%CondutivityFace,                                &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='CONDUTIVITYFACE',                         &
                     Default        = 1    ,                                    &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'GetUnSaturatedOptions - ModulePorousMedia - ERR13'


        call GetData(Me%SoilOpt%HCondFactor,                                    &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='HORIZONTAL_K_FACTOR',                     &
                     Default        = 1.0,                                      &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'GetUnSaturatedOptions - ModulePorousMedia - ERR16'

        call GetData(Me%SoilOpt%LimitEVAPWaterVelocity,                         &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'LIMIT_EVAP_WATER_VEL',                       &
                     Default    = .false.,                                      &                                           
                     ClientModule ='ModulePorousMedia',                         &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'GetUnSaturatedOptions - ModulePorousMedia - ERR17'

        !Number of iterations below which the DT is reduce (lower optimal iteration range)
        call GetData(Me%CV%MinIter,                                             &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='MIN_ITER',                                &
                     Default        = 2,                                        &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ConvergeOptions; ModulePorousMedia. ERR02") 

        !Number of iterations above which the DT is increased (upper optimal iteration
        !range)
        call GetData(Me%CV%MaxIter,                                             &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='MAX_ITER',                                &
                     Default        = 3,                                        &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ConvergeOptions; ModulePorousMedia. ERR03") 

        !Limit number of iteractions (when it is reached current dt is devided 
        !by 3 and the iteration starts from the begining
        call GetData(Me%CV%IterStep,                                            &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='ITER_STEP',                               &
                     Default        = 1,                                        &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ConvergeOptions; ModulePorousMedia. ERR04") 

        call GetData(Me%CV%LimitIter,                                           &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='LIMIT_ITER',                              &
                     Default        = 50,                                       &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ConvergeOptions; ModulePorousMedia. ERR05") 



        !Maximum diference allowed between THETA of each iteration
        call GetData(Me%CV%ThetaTolerance,                                      &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='THETA_TOLERANCE',                         &
                     Default        = 0.001,                                    &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ConvergeOptions; ModulePorousMedia. ERR06") 

        !Increase factor of currentDT when the number of iterations is smaller
        !then MinIter (lower time step multiplication factor)
        call GetData(Me%CV%IncreaseDT,                                          &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='INCREASE_DT',                             &
                     Default        = 1.20,                                     &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ConvergeOptions; ModulePorousMedia. ERR08") 

        !Increase factor of currentDT when the number of iterations is smaller
        !then MinIter (upper time step multiplication factor)
        call GetData(Me%CV%DecreaseDT,                                          &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='DECREASE_DT',                             &
                     Default        = 0.7,                                      &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ConvergeOptions; ModulePorousMedia. ERR09") 

        ! This value says how near the ThetaR the calculation is disconected. 
        ! Disables calculation when Theta is near ThetaR
        call GetData(Me%CV%LimitThetaLo,                                        &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='CUT_OFF_THETA_LOW',                       &
                     Default        = 0.001,                                    &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ConvergeOptions; ModulePorousMedia. ERR10") 
      
        ! This value says when Theta is converted to ThetaS
        ! Set Theta = ThetaS
        call GetData(Me%CV%LimitThetaHi,                                        &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='CUT_OFF_THETA_HIGH',                      &
                     Default        = 1.0e-15,                                  &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ConvergeOptions; ModulePorousMedia. ERR11") 
            
        ! This value says from which thetaS hydrostatic pressure is to be consider
        ! Set Theta = ThetaS
        call GetData(Me%CV%ThetaHydroCoef,                                      &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='THETA_HYDRO_COEF',                        &
                     Default        = 0.98,                                     &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ConvergeOptions; ModulePorousMedia. ERR12") 
        
    
    end subroutine ReadDataFile
    
    !--------------------------------------------------------------------------

    subroutine ReadInitialSoilFile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                                        :: Year_File, Month_File, Day_File
        real                                        :: Hour_File, Minute_File, Second_File
        integer                                     :: InitialFile
        type (T_Time)                               :: BeginTime, EndTimeFile, EndTime
        real                                        :: DT_error
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------


        call UnitsManager(InitialFile, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialSoilFile - ModulePorousMedia - ERR01'

        open(Unit = InitialFile, File = Me%Files%InitialFile, Form = 'UNFORMATTED', status = 'OLD', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialSoilFile - ModulePorousMedia - ERR02'

        !Reads Date
        read(InitialFile) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
        call SetDate(EndTimeFile, Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File)

        call GetComputeTimeLimits(Me%ObjTime, BeginTime = BeginTime, EndTime = EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialSoilFile - ModulePorousMedia - ERR03'
        
        DT_error = EndTimeFile - BeginTime

        !Avoid rounding erros
        if (abs(DT_error) >= 0.01) then
            
            write(*,*) 'The end time of the previous run is different from the start time of this run'
            write(*,*) 'Date in the file'
            write(*,*) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
            write(*,*) 'DT_error', DT_error
            if (Me%SoilOpt%StopOnWrongDate) stop 'ReadInitialSoilFile - ModulePorousMedia - ERR04'   

        endif

        read(InitialFile)Me%Theta

        call UnitsManager(InitialFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialSoilFile - ModulePorousMedia - ERR05'
        

    end subroutine ReadInitialSoilFile

    !--------------------------------------------------------------------------

    subroutine ConstructBottomTopography        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Constructs GridData
        call ConstructGridData      (Me%ObjBottomTopography, Me%ObjHorizontalGrid,  &
                                     FileName = Me%Files%BottomFile,                &
                                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR05'


    end subroutine ConstructBottomTopography

    !--------------------------------------------------------------------------
    
    subroutine VerticalDiscretization
                
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !get the topography
        call GetGridData      (Me%ObjTopography, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalDiscretization - ModulePorousMedia - ERR00'

         !Constructs GridData
        call ConstructGeometry (GeometryID       = Me%ObjGeometry,                 &
                                GridDataID       = Me%ObjBottomTopography,         &
                                HorizontalGridID = Me%ObjHorizontalGrid,           &
                                HorizontalMapID  = Me%ObjHorizontalMap,            &
                                ActualTime       = Me%ExtVar%Now,                  &
                                SurfaceElevation = Me%ExtVar%Topography,           &
                                BathymTopoFactor = -1.0,                           &
                                STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - VerticalDiscretization - ERR01'
        
        !Map - Soil Column            
        call ConstructMap       (Map_ID           = Me%ObjMap,                      &
                                 GeometryID       = Me%ObjGeometry,                 &
                                 HorizontalMapID  = Me%ObjHorizontalMap,            &
                                 TimeID           = Me%ObjTime,                     &
                                 STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - VerticalDiscretization - ERR02'
        
        !Get water points
        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - VerticalDiscretization - ERR03'


        !Geometry Size
        call GetGeometrySize    (Me%ObjGeometry,             &    
                                 Size     = Me%Size,         &
                                 WorkSize = Me%WorkSize,     &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalDiscretization - ModulePorousMedia - ERR05'

        Me%Size2D%ILB = Me%Size%ILB
        Me%Size2D%IUB = Me%Size%IUB
        Me%Size2D%JLB = Me%Size%JLB
        Me%Size2D%JUB = Me%Size%JUB

        !Initial geometry
        call ComputeInitialGeometry(Me%ObjGeometry,                                 &
                                    WaterPoints3D       = Me%ExtVar%WaterPoints3D,  &
                                    SurfaceElevation    = Me%ExtVar%Topography,     &
                                    ContinuesCompute    = .false.,                  &
                                    ActualTime          = Me%ExtVar%Now,            &
                                    STAT                = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalDiscretization - ModulePorousMedia - ERR06'
        
       
        !Checks Vertical Discretization
        call GetGeometryDistances(Me%ObjGeometry,                                       &
                                  DWZ         = Me%ExtVar%DWZ,                          &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalDiscretization - ModulePorousMedia. ERR07'

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%DWZ, STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'VerticalDiscretization - ModulePorousMedia. ERR09'

        call UnGetMap(Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "VerticalDiscretization; ModulePorousMedia. ERR10") 

        call UngetGridData (Me%ObjTopography, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalDiscretization - ModulePorousMedia - ERR11'

    end subroutine VerticalDiscretization

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
                 FILE   = '..\res\iter.soi_'//trim(adjustl(Number))//'.log', &
                 STATUS = "REPLACE",                                      &
                 IOSTAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                exit do1
            else
                Counter = Counter + 1
            end if
        enddo do1

        write (Me%Files%AsciiUnit, FMT=*) 'YY     MM   DD   HH   MM     SS       Iter  Time_Step '
        write (Me%Files%AsciiUnit, FMT=*) '                                                 s    '
    
    end subroutine ConstructAsciiOutPut
    
    !--------------------------------------------------------------------------

    subroutine AllocateVariables        
        
        !Local-----------------------------------------------------------------        
        integer                                         :: ILB, IUB, JLB,  JUB 
        integer                                         :: KLB, KUB 

        !Bounds
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KUB = Me%Size%KUB
        KLB = Me%Size%KLB
        
               
        !Water Content---------------------------------------------------------
        allocate (Me%Theta          (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%Head           (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%HydroPressure  (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%FinalHead      (ILB:IUB,JLB:JUB,KLB:KUB))
        
        allocate(Me%UGWaterLevel2D  (ILB:IUB,JLB:JUB)) 
        allocate(Me%UGWaterDepth2D  (ILB:IUB,JLB:JUB))
        allocate(Me%UGCell          (ILB:IUB,JLB:JUB)) 
        allocate(Me%WaterColumn     (ILB:IUB,JLB:JUB))
        allocate(Me%Infiltration    (ILB:IUB,JLB:JUB))
        allocate(Me%EfectiveEVTP    (ILB:IUB,JLB:JUB))
        allocate(Me%ImpermeableFraction (ILB:IUB,JLB:JUB))

        if (Me%SoilOpt%ComputeSoilField) then
            allocate (Me%ThetaField          (ILB:IUB,JLB:JUB,KLB:KUB))
        endif
          
        Me%Theta                = null_real
        Me%Head                 = null_real
        Me%HydroPressure        = null_real
        Me%FinalHead            = null_real

        Me%UGWaterLevel2D       = null_real
        Me%UGWaterDepth2D       = null_real
        Me%UGCell               = null_int
        Me%WaterColumn          = null_real
        Me%Infiltration         = null_real
        Me%EfectiveEVTP         = null_real
        Me%ImpermeableFraction  = null_real

        !Conductivities--------------------------------------------------------
        allocate (Me%SatK               (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%UnsatK             (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%UnsatK_X           (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%UnsatK_Y           (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%UnsatK_Z           (ILB:IUB,JLB:JUB,KLB:KUB))

        !SoilID
        allocate (Me%SoilID             (ILB:IUB,JLB:JUB,KLB:KUB))        
        
        !Speed Up Matrtix
        allocate (Me%CalculateHead      (ILB:IUB,JLB:JUB,KLB:KUB))
        
        
        Me%SatK            = null_real
        Me%UnsatK          = null_real
        Me%UnsatK_X        = null_real
        Me%UnsatK_Y        = null_real
        Me%UnsatK_Z        = null_real
        
        Me%SoilID          = null_int
        
        Me%CalculateHead   = .false.

        !Retention Curve-------------------------------------------------------
        allocate (Me%RC%ThetaR (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%RC%ThetaS (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%RC%ThetaF (ILB:IUB,JLB:JUB,KLB:KUB))
        
        Me%RC%ThetaR        = null_real
        Me%RC%ThetaS        = null_real
        Me%RC%ThetaF        = null_real
        
        !Converge method arrays------------------------------------------------
        allocate(Me%CV%HeadIni          (ILB:IUB,JLB:JUB,KLB:KUB))

        allocate(Me%CV%ThetaIni         (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%CV%ThetaOld         (ILB:IUB,JLB:JUB,KLB:KUB))
                
        !Velocities------------------------------------------------------------
        allocate(Me%UnsatVelU           (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%UnsatVelV           (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%UnsatVelW           (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%UnsatVelWOld        (ILB:IUB,JLB:JUB,KLB:KUB))
        
        allocate(Me%FluxU               (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxV               (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxW               (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxWFinal          (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxUAcc            (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxVAcc            (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxWAcc            (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxWAccFinal       (ILB:IUB,JLB:JUB,KLB:KUB))
        
        if (Me%ExtVar%ConstructEvaporation) then
            allocate(Me%EvaporationFlux     (ILB:IUB,JLB:JUB        ))
        endif
        
        allocate(Me%InfiltrationVelocity(ILB:IUB,JLB:JUB        ))

        Me%UnsatVelU            = 0.0
        Me%UnsatVelV            = 0.0
        Me%UnsatVelW            = 0.0
        Me%UnsatVelWOld         = 0.0
        
        Me%FluxU                = 0.0
        Me%FluxV                = 0.0
        Me%FluxW                = 0.0
        Me%FluxWFinal           = 0.0
        if (Me%ExtVar%ConstructEvaporation) then
            Me%EvaporationFlux      = 0.0
        endif

        Me%InfiltrationVelocity = null_real

        !Flow to Channel
        allocate(Me%iFlowToChannels                 (ILB:IUB,JLB:JUB))
        allocate(Me%lFlowToChannels                 (ILB:IUB,JLB:JUB))
        Me%iFlowToChannels     = 0.0
        Me%lFlowToChannels     = 0.0

        Me%FluxUAcc      = 0.0
        Me%FluxVAcc      = 0.0
        Me%FluxWAcc      = 0.0
        Me%FluxWAccFinal = 0.0
    end subroutine AllocateVariables

    !--------------------------------------------------------------------------
    
    subroutine ReadSoilTypes
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, nSoils
        logical                                     :: SoilTypeFound
        integer                                     :: ClientNumber, iflag, SoilID
        real                                        :: thf
        
        !Counts the number of SoilTypes
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        nSoils = 0
doS:    do
            !Gets soils type
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<beginsoiltype>',                    &
                                        block_end       = '<endsoiltype>',                      &
                                        BlockFound      = SoilTypeFound,                        &   
                                        STAT            = STAT_CALL)
SF:         if (STAT_CALL == SUCCESS_ .and. SoilTypeFound) then
                nSoils = nSoils + 1
           else

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR02'
               
                exit doS

            end if SF
        enddo doS            
        
        allocate(Me%SoilTypes(nSoils))
        
        
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)

doH:    do 
            !Gets soils type
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<beginsoiltype>',                    &
                                        block_end       = '<endsoiltype>',                      &
                                        BlockFound      = SoilTypeFound,                        &   
                                        STAT            = STAT_CALL)
HF:         if (STAT_CALL == SUCCESS_ .and. SoilTypeFound) then

                !Reads ID
                call GetData(SoilID, Me%ObjEnterData,  iflag,                                   &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'ID',                                             &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR10'
                if (iflag /= 1) then
                    write(*,*)'Missing ID in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR15'
                endif

                !Reads Theta R
                call GetData(Me%SoilTypes(SoilID)%ThetaR, Me%ObjEnterData,  iflag,              &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'THETA_R',                                        &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR20'
                if (iflag /= 1) then
                    write(*,*)'Missing THETA_R in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR30'
                endif
    
                !Reads Theta S
                call GetData(Me%SoilTypes(SoilID)%ThetaS, Me%ObjEnterData,  iflag,              &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'THETA_S',                                        &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR30'
                if (iflag /= 1) then
                    write(*,*)'Missing THETA_S in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR40'
                endif

                !Reads NFit
                call GetData(Me%SoilTypes(SoilID)%NFit, Me%ObjEnterData,  iflag,                &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'N_FIT',                                          &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR50'
                if (iflag /= 1) then
                    write(*,*)'Missing N_FIT in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR60'
                endif


                !Reads LFit
                call GetData(Me%SoilTypes(SoilID)%LFit, Me%ObjEnterData,  iflag,                &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'L_FIT',                                          &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR70'
                if (iflag /= 1) then
                    write(*,*)'Missing L_FIT in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR80'
                endif
                
                !Reads Alfa
                call GetData(Me%SoilTypes(SoilID)%Alfa, Me%ObjEnterData,  iflag,                &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'ALPHA',                                          &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR90'
                if (iflag /= 1) then
                    write(*,*)'Missing ALPHA in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR100'
                endif

                !Reads SatK
                call GetData(Me%SoilTypes(SoilID)%SatK, Me%ObjEnterData,  iflag,                &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'SAT_K',                                          &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR110'
                if (iflag /= 1) then
                    write(*,*)'Missing SAT_K in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR120'
                endif

                !Calculates MFIT
                Me%SoilTypes(SoilID)%MFit = 1.0 - (1.0 / Me%SoilTypes(SoilID)%NFit)
                
                !Calculates Oversat Head Slope
                thf = ThetaF_ (Me%SoilTypes(SoilID)%ThetaS - Me%CV%LimitThetaHi, SoilID)
                Me%SoilTypes(SoilID)%OverSatSlope = - Head_(thf, SoilID) / Me%CV%LimitThetaHi

           else

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR16'
               
                exit doH

            end if HF
        enddo doH
                
    
    end subroutine ReadSoilTypes

    !--------------------------------------------------------------------------   
    
    subroutine InitialFields

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, WLNumber
        integer                                     :: KLB, KUB, ClientNumber
        logical                                     :: HorizonFound, BlockFound
        logical                                     :: AllOK
        integer                                     :: nProps, iflag
        type (T_PropertyID)                         :: WaterLevelID
        type (T_PropertyID)                         :: ImpermeableFractionID
        integer, allocatable, dimension(:)          :: LayerControl
        integer, dimension(:,:,:), pointer          :: AuxPointsToFill
        real, dimension(:,:,:), pointer             :: AuxSoilID
        integer                                     :: i, j, k
        character(LEN = StringLength)               :: string
        type (T_PropertyID)                         :: ID
        
        
        allocate(LayerControl(Me%WorkSize%KLB: Me%WorkSize%KUB))
        LayerControl = 0

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)

doH:    do 
            !Gets soils horizons
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<beginhorizon>',                     &
                                        block_end       = '<endhorizon>',                       &
                                        BlockFound      = HorizonFound,                         &   
                                        STAT            = STAT_CALL)
HF:         if (STAT_CALL == SUCCESS_ .and. HorizonFound) then

                !Reads lower layer of horizon
                call GetData(KLB, Me%ObjEnterData,  iflag,                                      &
                            SearchType     = FromBlock,                                         &
                            keyword        = 'KLB',                                             &
                            ClientModule   = 'ModuleFillMatrix',                                &
                            STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR08'
                if (iflag /= 1) then
                    write(*,*)'Missing KLB in Horizon definition'
                    stop 'InitialFields - ModulePorousMedia - ERR09'
                endif

                !Reads upper layer of horizon       
                call GetData(KUB, Me%ObjEnterData,  iflag,                                      &
                            SearchType     = FromBlock,                                         &
                            keyword        = 'KUB',                                             &
                            ClientModule   = 'ModuleFillMatrix',                                &
                            STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR10'
                if (iflag /= 1) then
                    write(*,*)'Missing KUB in Horizon definition'
                    stop 'InitialFields - ModulePorousMedia - ERR11'
                endif

                do k = KLB, KUB
                    if (LayerControl(k) /= 0) then
                        write(*,*)'Inconsistent horizon definition. Layer:',k
                    else
                        LayerControl(k) = 1
                    endif
                enddo
                
                !Allocates Aux Matrix
                allocate(AuxPointsToFill (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))

                do i = Me%Size%ILB, Me%Size%IUB
                do j = Me%Size%JLB, Me%Size%JUB
                do k = Me%Size%KLB, Me%Size%KUB
                    if (k >= KLB .and. k <= KUB) then
                        AuxPointsToFill(i, j, k) = Me%ExtVar%WaterPoints3D(i, j, k)
                    else
                        AuxPointsToFill(i, j, k) = 0
                    endif
                enddo
                enddo
                enddo                

                nProps         = 0  !Control Variable
doSP:           do
                    !Gets properties of horizon
                    call ExtractBlockFromBlock(Me%ObjEnterData,                             &
                                                ClientNumber        = ClientNumber,         &
                                                block_begin         = '<beginproperty>',    &
                                                block_end           = '<endproperty>',      &
                                                BlockInBlockFound   = BlockFound,           &   
                                                STAT                = STAT_CALL)
                    if (STAT_CALL .EQ. SUCCESS_ .and. BlockFound) then                                                  

                        !Get Name of Property
                        call GetData(String, Me%ObjEnterData, iflag,        &
                                    SearchType     = FromBlockInBlock,     &
                                    keyword        ='NAME',                &
                                    ClientModule   ='ModulePOrousMedia',   &
                                    STAT           = STAT_CALL)  

                        select case (trim(adjustl(string)))
            
                        !Water contents
                        case (char_Theta        )                 
        
                            call ConstructFillMatrix(PropertyID           = ID,                               &
                                                     EnterDataID          = Me%ObjEnterData,                  &
                                                     TimeID               = Me%ObjTime,                       &
                                                     HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                                     GeometryID           = Me%ObjGeometry,                   &
                                                     ExtractType          = FromBlockInBlock,                 &
                                                     PointsToFill3D       = AuxPointsToFill,                  &
                                                     Matrix3D             = Me%Theta,                         &
                                                     TypeZUV              = TypeZ_,                           &
                                                     STAT                 = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR97'
        
                            call KillFillMatrix (ID%ObjFillMatrix, STAT = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR97a'

                        !SoilID
                        case (char_SoilID)

                            !Construct Fill Matrix cant read integer. So allocate aux matrix here.
                            allocate(AuxSoilID (Me%Size%ILB:Me%Size%IUB,                                      &
                                                Me%Size%JLB:Me%Size%JUB,                                      &
                                                Me%Size%KLB:Me%Size%KUB))
                                                
                            AuxSoilID = null_real

                            call ConstructFillMatrix(PropertyID           = ID,                               &
                                                     EnterDataID          = Me%ObjEnterData,                  &
                                                     TimeID               = Me%ObjTime,                       &
                                                     HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                                     GeometryID           = Me%ObjGeometry,                   &
                                                     ExtractType          = FromBlockInBlock,                 &
                                                     PointsToFill3D       = AuxPointsToFill,                  &
                                                     Matrix3D             = AuxSoilID,                        &
                                                     TypeZUV              = TypeZ_,                           &
                                                     STAT                 = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR98'
        
                            call KillFillMatrix (ID%ObjFillMatrix, STAT = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR98a'

                            do k = Me%Size%KLB, Me%Size%KUB
                            do j = Me%Size%JLB, Me%Size%JUB
                            do i = Me%Size%ILB, Me%Size%IUB
                                if (AuxPointsToFill(i, j, k) == 1) then
                                    Me%SoilID(i, j, k) = NINT(AuxSoilID(i, j, k))
                                endif
                            enddo
                            enddo
                            enddo
                            
                            deallocate (AuxSoilID)

                        !Invalid
                        case default

                            write(*,*)'Invalid Property', trim(adjustl(string))
                            stop 'InitialFields - ModulePorousMedia - ERR99'

                        end select                    
                    
                        nProps = nProps + 1
                    else
                        if (nProps /= 2) then
                            write(*,*)'Soil hydraulic properties incorrected defined'
                            stop 'InitialFields - ModulePorousMedia - ERR15'
                        endif
                        exit doSP
                    endif
                enddo doSP
                
                deallocate (AuxPointsToFill) 
                !Fills Constant Soil Variables
                

            else

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR16'
               
                exit doH

            end if HF
        enddo doH

        AllOK = .true.
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            if (LayerControl(k) /= 1) then
                AllOK = .false.
            endif
        enddo
        
        deallocate(LayerControl)

        if (.not. AllOK) then
            write(*,*)'Inconsistent horizon definition.'
            stop 'InitialFields - ModulePorousMedia - ERR61'
        endif
        
        do k = Me%Size%KLB, Me%Size%KUB
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB
            if (Me%ExtVar%Waterpoints3D(i, j, k) == 1) then
                if (Me%SoilID(i, j, k) < 0) then
                    write(*,*)'Soils not defined for [i,j,k]', i, j, k
                    AllOK = .false.
                endif
            endif
        enddo
        enddo
        enddo
        
        if (.not. AllOK) then
            write(*,*)'Inconsistent soil definition.'
            stop 'InitialFields - ModulePorousMedia - ERR61a'
        endif
        

        !Sets Matrixes of ThetaR and ThetaS        
        do k = Me%Size%KLB, Me%Size%KUB
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB
            if (Me%ExtVar%Waterpoints3D(i, j, k) == 1) then
                Me%RC%ThetaS(i, j, k) = Me%SoilTypes(Me%SoilID(i, j, k))%ThetaS
                Me%RC%ThetaR(i, j, k) = Me%SoilTypes(Me%SoilID(i, j, k))%ThetaR
                Me%SatK     (i, j, k) = Me%SoilTypes(Me%SoilID(i, j, k))%SatK
            endif
        enddo
        enddo
        enddo
        
      
        
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR70'

        !Constructs Water Level
        call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                    ClientNumber    = WLNumber,                             &
                                    block_begin     = '<beginwaterlevel>',                  &
                                    block_end       = '<endwaterlevel>',                    &
                                    BlockFound      = BlockFound,                           &   
                                    STAT            = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            if (.not. BlockFound) then
                write(*,*)'Missing Block <beginwaterlevel> / <endwaterlevel>'
                stop 'InitialFields - ModulePorousMedia - ERR80'
            endif
            
            call ConstructFillMatrix  ( PropertyID           = WaterLevelID,                &
                                        EnterDataID          = Me%ObjEnterData,             &
                                        TimeID               = Me%ObjTime,                  &
                                        HorizontalGridID     = Me%ObjHorizontalGrid,        &
                                        ExtractType          = FromBlock,                   &
                                        PointsToFill2D       = Me%ExtVar%BasinPoints,       &
                                        Matrix2D             = Me%UGWaterLevel2D,           &
                                        TypeZUV              = TypeZ_,                      &
                                        STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR81'
            
            call KillFillMatrix       (WaterLevelID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR82'

        else
            stop 'InitialFields - ModulePorousMedia - ERR90'
        endif
        
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR70'

        !Constructs Impermeable Fraction
        call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                    ClientNumber    = WLNumber,                             &
                                    block_begin     = '<beginimpermeablefraction>',         &
                                    block_end       = '<endimpermeablefraction>',           &
                                    BlockFound      = BlockFound,                           &   
                                    STAT            = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            if (.not. BlockFound) then
                write(*,*)'Missing Block <beginimpermeablefraction> / <endimpermeablefraction>'
                stop 'InitialFields - ModulePorousMedia - ERR100'
            endif
            
            call ConstructFillMatrix  ( PropertyID           = ImpermeableFractionID,       &
                                        EnterDataID          = Me%ObjEnterData,             &
                                        TimeID               = Me%ObjTime,                  &
                                        HorizontalGridID     = Me%ObjHorizontalGrid,        &
                                        ExtractType          = FromBlock,                   &
                                        PointsToFill2D       = Me%ExtVar%BasinPoints,       &
                                        Matrix2D             = Me%ImpermeableFraction,      &
                                        TypeZUV              = TypeZ_,                      &
                                        STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR81'
            
            call KillFillMatrix       (ImpermeableFractionID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR82'

        else
            stop 'InitialFields - ModulePorousMedia - ERR110'
        endif


    end subroutine InitialFields

    !--------------------------------------------------------------------------

    subroutine StartWithFieldCapacity

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i,j,k
        real                                        :: Hinf
        integer                                     :: WTCell
        real                                        :: inf_border 
        
        do j= Me%WorkSize%JLB, Me%WorkSize%JUB
        do i= Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(i,j,Me%WorkSize%KUB) == 1) then
            
                !Searches cell were initial WT is located
                WTCell  = 0
                do k = Me%WorkSize%KUB, Me%ExtVar%KFloor(i, j), -1
                    inf_border = - Me%ExtVar%SZZ(i,j,k-1)
                    if (Me%UGWaterLevel2D(i,j) .gt. inf_border ) then
                        WTCell  = k
                        exit
                    endif                
                enddo        
                
                if (WTCell < Me%ExtVar%KFloor(i, j)) WTCell = Me%ExtVar%KFloor(i, j)
            
                !Saturates cells below WTable
                do k = Me%ExtVar%KFloor(i, j), WTCell-1
                    Me%Theta(i, j, k) = Me%RC%ThetaS(i, j, k)
                enddo


                if (Me%SoilOpt%StartWithFieldCapacity) then
                
                    !Sets cell above WT to Field Capacity
                    do k = WTCell, Me%WorkSize%KUB

                        if (k == WTCell) then
                            !Half of the cell distance
                            Hinf = - Me%ExtVar%DWZ(i,j,k) * 0.5
                        else
                            Hinf = - (Me%ExtVar%DZZ(i,j,k-1) - Hinf)
                        endif
                        
                        Me%Theta(i, j, k) = Theta_(Hinf, Me%SoilID(i, j, k))
                        
                    enddo
                
                endif
                    
            endif
        enddo
        enddo        
                                

    end subroutine StartWithFieldCapacity

    !--------------------------------------------------------------------------

    subroutine ComputeSoilFieldCapacity

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: i,j,k
        real                                                :: Hinf
        !Begin-----------------------------------------------------------------

        do j= Me%WorkSize%JLB, Me%WorkSize%JUB
        do i= Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(i,j,Me%WorkSize%KUB) == 1) then
            
                do k = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB

                    if (k == Me%ExtVar%KFloor(i, j)) then
                        !Half of the cell distance
                        Hinf = - Me%ExtVar%DWZ(i,j,k) * 0.5
                    else
                        Hinf = - (Me%ExtVar%DZZ(i,j,k-1) - Hinf)
                    endif
                    
                    Me%ThetaField(i, j, k) = Theta_(Hinf, Me%SoilID(i, j, k))
                    
                enddo
                
                    
            endif
        enddo
        enddo        

    end subroutine ComputeSoilFieldCapacity

    !--------------------------------------------------------------------------


    subroutine ConstructHDF5Output        

        !Local-----------------------------------------------------------------
        integer                                             :: ILB,IUB,JLB,JUB,KLB,KUB    
        integer                                             :: STAT_CALL
        integer                                             :: HDF5_CREATE
        real, dimension(:, :), pointer                      :: BottomData

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%TransientHDF)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR01'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR02'

        !Sets limits for next write operations
        call HDF5SetLimits      (Me%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR03'

        !Writes the Grid
        call HDF5WriteData      (Me%ObjHDF5, "/Grid", "Topography", "m",                    &
                              Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR04'

        call GetGridData(Me%ObjBottomTopography, BottomData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR05'

        call HDF5WriteData      (Me%ObjHDF5, "/Grid", "Bathymetry", "m",                    &
                              Array2D = BottomData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR06'

        call HDF5WriteData      (Me%ObjHDF5, "/Grid", "BasinPoints", "-",                   &
                              Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR07'

        !Water Points
        if (Me%OutPut%Yes) then
            call HDF5WriteData   ( Me%ObjHDF5,  "/Grid", "WaterPoints3D", "-",                  &
                                   Array3D = Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR08'
        endif              
                
        !Flushes All pending HDF5 commands
        call HDF5FlushMemory    (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR09'

        call UnGetGridData(Me%ObjBottomTopography, BottomData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR10'

    end subroutine ConstructHDF5Output

    !--------------------------------------------------------------------------

    subroutine ConstructTimeSerie

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        integer                                             :: nProperties
        integer                                             :: STAT_CALL
        integer                                             :: iflag, i
        character(len=StringLength)                         :: TimeSerieLocationFile
        !Begin-----------------------------------------------------------------
        
        nProperties = 10
        if (Me%ExtVar%ConstructEvaporation) then
            nProperties = nProperties + 1
        endif
        if (Me%ExtVar%ConstructTranspiration) then
            nProperties = nProperties + 1
        endif


        !Allocates PropertyList
        allocate(PropertyList(nProperties))

        !Fills up PropertyList
        PropertyList(1)  = 'Theta'
        PropertyList(2)  = 'relative water content'
        PropertyList(3)  = 'VelW'
        PropertyList(4)  = 'InF_Vel'
        PropertyList(5)  = 'Head'
        PropertyList(6)  = 'Conductivity'
        PropertyList(7)  = 'level water table'
        PropertyList(8)  = 'water table depth'
        PropertyList(9)  = 'Hydro Pressure'
        PropertyList(10) = 'Final Head'
!        PropertyList(11) = 'PM Water Column m'
        i = 10
        if (Me%ExtVar%ConstructEvaporation) then
            i = i + 1
            PropertyList(i) = 'Surface Evaporation Flux m3/s'
        endif
        if (Me%ExtVar%ConstructTranspiration) then
            i = i + 1
            PropertyList(i) = 'Transpiration Flux m3/s'
        endif

        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModulePorousMedia',                                &
                     Default      = Me%Files%DataFile,                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR01' 

        if (iflag == 1) then
            Me%OutPut%TimeSerieON = .true.
        else
            Me%OutPut%TimeSerieON = .false.
        endif

        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                            TimeSerieLocationFile,                                      &
                            PropertyList, "srp",                                        &
                            WaterPoints3D = Me%ExtVar%WaterPoints3D,                    &
                            ModelName = Me%ModelName,                                   &
                            STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR02' 

        !Deallocates PropertyList
        deallocate(PropertyList)
       
    end subroutine ConstructTimeSerie

    !--------------------------------------------------------------------------

    subroutine ConstructProfileOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL, iflag
        character(len=StringLength)                         :: TimeSerieLocationFile
        character(len=StringLength), dimension(:,:), pointer:: PropertyList
        integer                                             :: nProperties

        nProperties = 3 !Theta, ThetaF, Head

        !Allocates PropertyList
        allocate(PropertyList(nProperties, 2))

        !Fills up PropertyList
        PropertyList(1, 1) = 'PropI'
        PropertyList(2, 1) = 'PropII'

        PropertyList(1, 2) = "m3/m3"
        PropertyList(2, 2) = "m3/m3"


        !----------------------------------------------------------------------

        call GetData(Me%OutPut%ProfileON,                                               &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'OUTPUT_PROFILE',                                   &
                     ClientModule = 'ModulePorousMedia',                                &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileOutput - ModulePorousMedia - ERR01' 


        if (Me%OutPut%ProfileON) then

            call GetData(TimeSerieLocationFile,                                             &
                         Me%ObjEnterData,iflag,                                             &
                         SearchType   = FromFile,                                           &
                         keyword      = 'TIME_SERIE_LOCATION',                              &
                         ClientModule = 'ModulePorousMedia',                                &
                         Default      = Me%Files%DataFile,                                  &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileOutput - ModulePorousMedia - ERR02' 
            
            !Starts Profile for Theta / ThetaF
            call StartProfile  (ProfileID       = Me%ObjProfile,                            &
                                ObjTime         = Me%ObjTime,                               &
                                ProfileDataFile = trim(TimeSerieLocationFile),              &
                                WaterPoints2D   = Me%ExtVar%BasinPoints,                    &
                                nProperties     = 3,                                        &
                                PropertyList    = PropertyList,                             &
                                KUB             = Me%WorkSize%KUB,                          &
                                ClientName      = "PorousMedia",                            &
                                STAT            = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileOutput - ModulePorousMedia - ERR03' 

        endif

        deallocate (PropertyList)
        
    end subroutine ConstructProfileOutput


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !---------------------------------------------------------------------------

    subroutine GetNextPorousMediaDT (ObjPorousMediaID, DT, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real, intent(OUT)                               :: DT
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            DT        = Me%NextDT

            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetNextPorousMediaDT

    !--------------------------------------------------------------------------
    
    subroutine GetPotentialInfiltration (ObjPorousMediaID, PotInfiltration, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:, :), pointer               :: PotInfiltration
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            PotInfiltration => Me%ExtVar%InfiltrationColumn

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetPotentialInfiltration

    !--------------------------------------------------------------------------    

    subroutine GetInfiltration (ObjPorousMediaID, Infiltration, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:, :), pointer               :: Infiltration
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            Infiltration => Me%Infiltration

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetInfiltration

    !--------------------------------------------------------------------------

    subroutine GetEfectiveEVTP (ObjPorousMediaID, EfectiveEVTP,STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:, :), pointer               :: EfectiveEVTP
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            EfectiveEVTP => Me%EfectiveEVTP

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetEfectiveEVTP

    !--------------------------------------------------------------------------

    subroutine GetGWFlowToChannels (ObjPorousMediaID, FlowToChannels, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real, dimension(:, :), pointer                  :: FlowToChannels
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        

        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            FlowToChannels => Me%iFlowToChannels

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetGWFlowToChannels

    !--------------------------------------------------------------------------

    subroutine GetGWLayer (ObjPorousMediaID, GWLayer, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer, dimension(:,:), pointer                :: GWLayer
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        

        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            GWLayer => Me%UGCell

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetGWLayer

    !--------------------------------------------------------------------------

    subroutine GetTotalStoredVolume (ObjPorousMediaID, TotalStoredVolume, LossToGround, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8)                                         :: TotalStoredVolume
        real(8), optional                               :: LossToGround
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            TotalStoredVolume = Me%TotalStoredVolume

            if (present(LossToGround)) LossToGround      = Me%LossToGround

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetTotalStoredVolume

    !-------------------------------------------------------------------------
    subroutine GetGeometryInstance (ObjPorousMediaID, GeometryID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer, intent(OUT), optional                  :: STAT
        integer                                         :: GeometryID

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)  !Changed 

            GeometryID = Me%ObjGeometry

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetGeometryInstance

    !--------------------------------------------------------------------------


    subroutine GetFluxU (ObjPorousMediaID, FlowU, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:,:,:), pointer              :: FlowU
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            FlowU => Me%FluxUAcc

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetFluxU

    !--------------------------------------------------------------------------

    subroutine GetFluxV (ObjPorousMediaID, FlowV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:,:,:), pointer              :: FlowV
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            FlowV => Me%FluxVAcc

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetFluxV

    !--------------------------------------------------------------------------

    subroutine GetFluxW (ObjPorousMediaID, FluxWFinal, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:,:,:), pointer              :: FluxWFinal
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            ! Lucia
            !FlowWOld => Me%FluxWold
            
            !changed by Eduardo Jauch
            FluxWFinal => Me%FluxWAccFinal

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetFluxW

    !--------------------------------------------------------------------------

    subroutine GetUnsatV (ObjPorousMediaID, UnsatV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,dimension(:,:,:), pointer                  :: UnsatV                        
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            !UnsatW => Me%UnsatVelWOld
            UnsatV => Me%UnsatVelV

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetUnsatV
    
    !--------------------------------------------------------------------------

    subroutine GetUnsatU (ObjPorousMediaID, UnsatU, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,dimension(:,:,:), pointer                  :: UnsatU                        
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            !UnsatW => Me%UnsatVelWOld
            UnsatU => Me%UnsatVelU

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetUnsatU
    
    !--------------------------------------------------------------------------

    subroutine GetUnsatW (ObjPorousMediaID, UnsatW, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,dimension(:,:,:), pointer                  :: UnsatW                        
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            UnsatW => Me%UnsatVelW

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetUnsatW

    !---------------------------------------------------------------------------

    subroutine GetWaterColumn (ObjPorousMediaID, WaterColumn, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                          :: ObjPorousMediaID
        real,    pointer, dimension(:,:) :: WaterColumn
        integer, intent(OUT), optional   :: STAT

        !Local-----------------------------------------------------------------
        integer :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            WaterColumn => Me%WaterColumn

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_    
    
    end subroutine GetWaterColumn

    !---------------------------------------------------------------------------

    subroutine GetWaterContent (ObjPorousMediaID, WC, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: WC
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            WC => Me%Theta

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetWaterContent

    !--------------------------------------------------------------------------

    subroutine GetOldWaterContent (ObjPorousMediaID, WCold, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: WCold
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            WCold => Me%CV%ThetaIni !Me%CV%ThetaOld

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetOldWaterContent

    !--------------------------------------------------------------------------

    subroutine GetHead (ObjPorousMediaID, Head, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: Head
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            Head => Me%Head

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine GetHead

    !--------------------------------------------------------------------------

    subroutine GetThetaR (ObjPorousMediaID, ThetaR, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: ThetaR
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            ThetaR => Me%RC%ThetaR

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine GetThetaR

    !--------------------------------------------------------------------------

    subroutine GetThetaS (ObjPorousMediaID, ThetaS, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: ThetaS
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            ThetaS => Me%RC%ThetaS

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine GetThetaS

    !--------------------------------------------------------------------------

    subroutine GetUnsatK (ObjPorousMediaID, UnsatK, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: UnsatK
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            UnsatK => Me%UnsatK

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine GetUnsatK

    !--------------------------------------------------------------------------

    subroutine GetEvaporation (ObjPorousMediaID, Evaporation, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real   ,    pointer, dimension(:,:)             :: Evaporation
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            Evaporation => Me%EvaporationFlux

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine GetEvaporation

    !--------------------------------------------------------------------------

    subroutine GetThetaField (ObjPorousMediaID, ThetaField, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: ThetaField
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            ThetaField => Me%ThetaField

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetThetaField

    !--------------------------------------------------------------------------

    subroutine GetComputeSoilField (ObjPorousMediaID, ComputeSoilField, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        logical                                         :: ComputeSoilField
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            ComputeSoilField = Me%SoilOpt%ComputeSoilField

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetComputeSoilField

    !--------------------------------------------------------------------------

    subroutine GetLimitThetaLow (ObjPorousMediaID, LimitThetaLow, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real                                            :: LimitThetaLow
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            LimitThetaLow = Me%CV%LimitThetaLo

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetLimitThetaLow

    !--------------------------------------------------------------------------

    subroutine UnGetPorousMedia_R4(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(4), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_R4")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMedia_R4

    !--------------------------------------------------------------------------

    subroutine UnGetPorousMedia_R8(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_R8")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMedia_R8
       
    !----------------------------------------------------------------------------
    subroutine UnGetPorousMedia_R(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), pointer, dimension(:,:,:)               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_R8")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMedia_R
    !-----------------------------------------------------------------------------

    subroutine UnGetPorousMedia_R1(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(4), pointer, dimension(:,:,:)              :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_R8")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMedia_R1

!-------------------------------------------------------------------------------

    subroutine UnGetPorousMedia_RI(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer                                         :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            array = 0
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_RI")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        end subroutine UnGetPorousMedia_RI

    ! ---------------------------------------------------------------------!

     subroutine UnGetPorousMedia_AI(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer, pointer, dimension(:,:,:)              :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_AI")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        end subroutine UnGetPorousMedia_AI

    ! ---------------------------------------------------------------------!

     subroutine UnGetPorousMedia_AI2D(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer, pointer, dimension(:,:)                :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_AI2D")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        end subroutine UnGetPorousMedia_AI2D

    ! ---------------------------------------------------------------------!
      
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyPorousMedia(ObjPorousMediaID,                         &
                                 InfiltrationColumn,                       &
                                 PotentialEvaporation,                     &
                                 ActualTranspiration,                      &
                                 STAT)

        !Arguments---------------------------------------------------------------
        integer,                            intent(IN)           :: ObjPorousMediaID       !IN
        real(8), dimension(:,:  ), pointer                       :: InfiltrationColumn     !IN
        real,    dimension(:,:  ), pointer, optional             :: PotentialEvaporation   !IN
        real,    dimension(:,:,:), pointer, optional             :: ActualTranspiration    !IN  
        integer, optional,                  intent(OUT)          :: STAT                   !OUT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_, ready_, STAT_CALL
        integer                                     :: DummyI
        real                                        :: DummyR
        !------------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "ModifyPorousMedia")

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then                        

            !Updates time
            call GetComputeCurrentTime  (Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
            
            if (STAT_CALL /= SUCCESS_) stop 'ModifyPorousMedia - ModulePorousMedia - ERR01'
            
            !Gets Time Step
            call GetComputeTimeStep (Me%ObjTime, DT = Me%ExtVar%DT)                


            !Update unsaturated computefaces
            call UpdateComputeFaces3D(  Map_ID         = Me%ObjMap,                       &
                                        DummyR         = DummyR,                          &
                                        DummyI         = DummyI,                          &
                                        STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyPorousMedia - ModulePorousMedia - ERR02'


            !Sets External Variables
            call ReadLockExternalVar
            
            !Points to Arguments
            Me%TranspirationExists = .false.
            Me%EvaporationExists   = .false.
            if (present(ActualTranspiration)) then
                Me%TranspirationExists               = .true.
                !m3/s
                Me%ExtVar%TranspirationFlux          => ActualTranspiration
            endif
            if (present(PotentialEvaporation)) then
                Me%EvaporationExists                 = .true.
                !m/s
                Me%ExtVar%PotentialEvaporationFlux   => PotentialEvaporation
            endif
            
            Me%ExtVar%InfiltrationColumn => InfiltrationColumn
       
            !Calculate flow in unsaturated part of soil
            call VariableSaturatedFlow (InfiltrationColumn)
            
            call CalculateUGWaterLevel

            !Output
            if (Me%OutPut%Yes .or. Me%OutPut%SurfaceOutput) call PorousMediaOutput                        
            if (Me%OutPut%TimeSerieON)                      call OutPutTimeSeries
            if (Me%OutPut%ProfileON  )                      call ProfileOutput
            
            !Restart Output
            if (Me%Output%WriteRestartFile .and. .not. (Me%ExtVar%Now == Me%EndTime)) then
                if(Me%ExtVar%Now >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then
                    call WriteFinalSoilFile
                    Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1
                endif
            endif

            
            if (Me%SoilOpt%CheckGlobalMass) then
                call CalculateTotalStoredVolume
            endif

            call ReadUnLockExternalVar
            
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "ModifyPorousMedia")

    end subroutine ModifyPorousMedia

    !--------------------------------------------------------------------------
    
    subroutine VariableSaturatedFlow(InfiltrationColumn)
    
        !Arguments-------------------------------------------------------------
        real(8), dimension(:,:), pointer            :: InfiltrationColumn
        
        !Local-----------------------------------------------------------------
        logical                                     :: StrongVariation            
        integer                                     :: Niteration, iteration
        real                                        :: SumDT
        real                                        :: Zero = 0.0
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "VariableSaturatedFlow")

        !Stores initial values
        call SetMatrixValue (Me%CV%ThetaIni,      Me%Size,   Me%Theta,          Me%ExtVar%WaterPoints3D)
        call SetMatrixValue (Me%CV%HeadIni,       Me%Size,   Me%Head,           Me%ExtVar%WaterPoints3D)
        call SetMatrixValue (Me%WaterColumn,      Me%Size2D, InfiltrationColumn,Me%ExtVar%BasinPoints)

        !Time Integrated Values
        call SetMatrixValue (Me%Infiltration,     Me%Size2D, dble(0.0),         Me%ExtVar%BasinPoints)
        call SetMatrixValue (Me%EfectiveEVTP,     Me%Size2D, dble(0.0),         Me%ExtVar%BasinPoints)
!        call SetMatrixValue (Me%EfectiveEVTP2,     Me%Size2D, dble(0.0),         Me%ExtVar%BasinPoints)       
        call SetMatrixValue (Me%iFlowToChannels,  Me%Size2D, Zero,              Me%ExtVar%BasinPoints)
        
!        call SetMatrixValue (Me%FluxWAcc, Me%Size, dble(0.0), Me%ExtVar%WaterPoints3D) 
!        call SetMatrixValue (Me%FluxVAcc, Me%Size, dble(0.0), Me%ExtVar%WaterPoints3D) 
!        call SetMatrixValue (Me%FluxUAcc, Me%Size, dble(0.0), Me%ExtVar%WaterPoints3D) 
!        call SetMatrixValue (Me%FluxWAccFinal, Me%Size, dble(0.0), Me%ExtVar%WaterPoints3D) 
        
        iteration         = 1
        Niteration        = 1
        Me%CV%CurrentDT   = Me%ExtVar%DT / Niteration
        SumDT             = 0.0
        StrongVariation   = .false.

        Me%FluxWAcc = 0.
        Me%FluxVAcc = 0.
        Me%FluxUAcc = 0.
        Me%FluxWAccFinal = 0.
        
dConv:  do while (iteration <= Niteration)
        
            !Convergence Test
            call SetMatrixValue (Me%CV%ThetaOld, Me%Size, Me%Theta, Me%ExtVar%WaterPoints3D)

            !Calculates Face Conductivities
            call Condutivity_Face
            
            !Calculates Water velocity
            call SoilWaterVelocity
           
            !Calculates Water Flux
            call SoilWaterFlux
            
            !Calculates Flux to channels
            if (Me%ObjDrainageNetwork /= 0) then
                call ExchangeWithDrainageNetwork ()
            endif
            
            if (Me%EvaporationExists) then
                call EvaporationFlux    ()
            endif
            
            !Calculates New Theta
            call CalculateNewTheta

            call variation_test     (StrongVariation)

            !Vertical Continuty            
            if (StrongVariation) then
            
                Niteration        = Niteration + Me%CV%IterStep
                Me%CV%CurrentDT   = Me%ExtVar%DT / Niteration
                call WriteDTLog ('ModulePorousMedia', Niteration, Me%CV%CurrentDT)                    
                iteration         = 1
                SumDT             = 0.0
                StrongVariation   = .false.
                
                !Restores Initial Values
                call SetMatrixValue (Me%Theta,          Me%Size,   Me%CV%ThetaIni,      Me%ExtVar%WaterPoints3D)
                call SetMatrixValue (Me%Head,           Me%Size,   Me%CV%HeadIni,       Me%ExtVar%WaterPoints3D)
                call SetMatrixValue (Me%WaterColumn,    Me%Size2D, InfiltrationColumn,  Me%ExtVar%BasinPoints)

                !Resets Time Integrated Values
                call SetMatrixValue (Me%Infiltration,   Me%Size2D, dble(0.0),           Me%ExtVar%BasinPoints)
                call SetMatrixValue (Me%EfectiveEVTP,   Me%Size2D, dble(0.0),           Me%ExtVar%BasinPoints)
                call SetMatrixValue (Me%iFlowToChannels,Me%Size2D, Zero,                Me%ExtVar%BasinPoints)
                
                !Resets Accumulated flows
!                call SetMatrixValue (Me%FluxWAcc, Me%Size, dble(0.0), Me%ExtVar%WaterPoints3D) 
!                call SetMatrixValue (Me%FluxVAcc, Me%Size, dble(0.0), Me%ExtVar%WaterPoints3D) 
!                call SetMatrixValue (Me%FluxUAcc, Me%Size, dble(0.0), Me%ExtVar%WaterPoints3D) 
!                call SetMatrixValue (Me%FluxWAccFinal, Me%Size, dble(0.0), Me%ExtVar%WaterPoints3D)  

                Me%FluxWAcc = 0.
                Me%FluxVAcc = 0.
                Me%FluxUAcc = 0.
                Me%FluxWAccFinal = 0.
                              
            else
                
                call VerticalContinuity

                !Calulates Heads / Conductivities from new Theta values
                call SoilParameters    (Me%ExtVar%BasinPoints)

                !Removes water due to infiltration
                call IntegrateValuesInTime(SumDT)
                
                !Advection of Properties
!                call CalculateAdvection
                
                !Accumulate flows
                call AccumulateFlows
                
                
                SumDT       = SumDT + Me%CV%CurrentDT
                iteration   = iteration + 1
                
            endif
            
            if (Niteration > Me%CV%LimitIter) then
                    call SetError(WARNING_, INTERNAL_, "Strong Variation after LIMIT ITER", OFF)
            endif
            
        enddo dConv
        
        call CalculateMeanFlows (Niteration)
        call InsertInfiltrationOnFluxMatrix
        
        call LogDT (Niteration)
        
        call PredictDT(Niteration)

        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "VariableSaturatedFlow")   

    end subroutine VariableSaturatedFlow
    
    !--------------------------------------------------------------------------
    
    subroutine InsertInfiltrationOnFluxMatrix
    
        !Local-----------------------------------------------------------------
        integer :: i, j, k
        real    :: infiltration_flux

        !----------------------------------------------------------------------
        
        k = Me%WorkSize%KUB + 1
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB 
        
            if (Me%ExtVar%BasinPoints (i, j) == 1) then
        
                infiltration_flux = Me%Infiltration(i, j) * Me%ExtVar%Area(i, j) / Me%ExtVar%DT
                Me%FluxWAccFinal(i, j, k) = Me%FluxWAccFinal(i, j, k) - Infiltration_flux
                
            endif
        
        enddo
        enddo           
        
        !----------------------------------------------------------------------
    
    end subroutine InsertInfiltrationOnFluxMatrix
    
    !--------------------------------------------------------------------------
    
    subroutine CalculateMeanFlows (iterations)
    
        !Arguments-------------------------------------------------------------
        integer :: iterations

        !Local-----------------------------------------------------------------
        integer :: i, j, k

        !----------------------------------------------------------------------
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB            
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        
            if (Me%ExtVar%ComputeFacesU3D(I,J,K) .EQ. Compute) then
                
               Me%FluxUAcc(i, j, k) = Me%FluxUAcc(i, j, k) / iterations                 
            
            endif
            
            if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then
            
                Me%FluxVAcc(i, j, k) = Me%FluxVAcc(i, j, k) / iterations   
            
            endif

            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then
            
                Me%FluxWAcc(i, j, k) = Me%FluxWAcc(i, j, k) / iterations
                Me%FluxWAccFinal(i, j, k) = Me%FluxWAccFinal(i, j, k) / iterations
            
            endif
            
        enddo
        enddo
        enddo
        
        Me%FluxWAcc(i, j, Me%WorkSize%KUB+1) = Me%FluxWAcc(i, j, Me%WorkSize%KUB+1) / iterations
             
        !----------------------------------------------------------------------
        
    
    end subroutine CalculateMeanFlows
    
    subroutine AccumulateFlows

        !Local-----------------------------------------------------------------
        integer :: i, j, k

        !----------------------------------------------------------------------
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB            
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        
            if (Me%ExtVar%ComputeFacesU3D(I,J,K) .EQ. Compute) then
                
               Me%FluxUAcc(i, j, k) = Me%FluxUAcc(i, j, k) + Me%FluxU(i, j, k)                 
            
            endif
            
            if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then
            
                Me%FluxVAcc(i, j, k) = Me%FluxVAcc(i, j, k) + Me%FluxV(i, j, k)    
            
            endif

            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then
            
                Me%FluxWAcc(i, j, k) = Me%FluxWAcc(i, j, k) + Me%FluxW(i, j, k)
                Me%FluxWAccFinal(i, j, k) = Me%FluxWAccFinal(i, j, k) + Me%FluxWFinal(i, j, k)
            
            endif
            
        enddo
        enddo
        enddo
        
        Me%FluxWAcc(i, j, Me%WorkSize%KUB+1) = Me%FluxWAcc(i, j, Me%WorkSize%KUB+1) + Me%FluxW(i, j, Me%WorkSize%KUB+1)
             
        !----------------------------------------------------------------------
        
    end subroutine AccumulateFlows
    
    !--------------------------------------------------------------------------

    subroutine EffectiveVelocity

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
    
    
    end subroutine EffectiveVelocity

    !--------------------------------------------------------------------------

    subroutine PredictDT(Niter)

        !Arguments-------------------------------------------------------------
        integer                                     :: Niter

        !Local-----------------------------------------------------------------

        Me%NextDT = Me%ExtVar%DT
        
        if (Niter < Me%CV%MinIter) then
            Me%NextDT = Me%NextDT * Me%CV%IncreaseDT
        else if (Niter > Me%CV%MaxIter) then
            Me%NextDT = Me%NextDT * Me%CV%DecreaseDT
        endif
   
    end subroutine PredictDT 

    !----------------------------------------------------------------------------

    subroutine SoilWaterVelocity

        !Arguments---------------------------------------------------------------
                  
        !Local-------------------------------------------------------------------
        integer                                     :: I, J, K                
        integer                                     :: Chunk

        !------------------------------------------------------------------------          

        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "SoilWaterVelocity")

        call FinalHead

        !$OMP PARALLEL PRIVATE(I,J,K)

        !Horizontal Velocities
        if (Me%SoilOpt%CalcHorizontal) then
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB            

                if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then
                    Me%UnsatVelV(I,J,K) = BuckinghamDarcyEquation               &
                                       (con      = Me%UnsatK_Y(I,  J,K),        &
                                        hinf     = Me%FinalHead (I-1,J,K),      &
                                        hsup     = Me%FinalHead (I,  J,K),      &
                                        delta    = Me%ExtVar%DZY (I-1,J  ))
                end if

                if (Me%ExtVar%ComputeFacesU3D(I,J,K) .EQ. Compute) then
                    Me%UnsatVelU(I,J,K) = BuckinghamDarcyEquation               &
                                       (con      = Me%UnsatK_X(I,J,K  ),        &
                                        hinf     = Me%FinalHead (I,J-1,K),      &
                                        hsup     = Me%FinalHead (I,J,  K),      &
                                        delta    = Me%ExtVar%DZX (I,J-1  ))
                end if
            end do
            end do
            end do
            !$OMP END DO NOWAIT
        endif
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB            
        
            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then
            
!               Me%UnsatVelWold(I,J,K) = Me%UnsatVelW(I,J,K)
                
                Me%UnsatVelW(I,J,K) = BuckinghamDarcyEquation           &
                                   (con      = Me%UnsatK_Z(I,J,K  ),    &
                                    hinf     = Me%FinalHead (I,J,K-1),  &
                                    hsup     = Me%FinalHead (I,J,K  ),  &
                                    delta    = Me%ExtVar%DZZ (I,J,K-1))
            else

!               Me%UnsatVelWold(I,J,K) = Me%UnsatVelW(I,J,K)

                Me%UnsatVelW(I,J,K) = 0.0

            end if
        end do
        end do
        end do
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

        call InfiltrationVelocity


        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "SoilWaterVelocity")

    end subroutine SoilWaterVelocity

    !--------------------------------------------------------------------------
    
    subroutine FinalHead

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        integer                                     :: CHUNK
        real                                        :: AccumPressure, Coef
        real                                        :: CenterVelocityW
        integer                                     :: KUB

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        KUB   = Me%WorkSize%KUB

        !$OMP PARALLEL PRIVATE(I,J,K,AccumPressure, Coef)
        !$OMP DO SCHEDULE    (DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints (i, j) == 1) then
            
                !Initial pressure = watercolumn - move downwards
                AccumPressure = Me%WaterColumn(i, j)
                do k = Me%WorkSize%KUB, Me%ExtVar%KFloor(i, j), -1
                
                    if (Me%Theta(i, j, k) > Me%RC%ThetaS(i, j, k) * Me%CV%ThetaHydroCoef) then
                    
                        !Velocity at the center of the cell
                        CenterVelocityW  = (Me%UnsatVelW(i, j, k) + Me%UnsatVelW(i, j, k+1)) / 2.0

                        
                        Coef = 0.5 * Me%ExtVar%DWZ(i, j, k) * (1.0 - CenterVelocityW / Me%SatK(i, j, k)) * &
                               LinearInterpolation(Me%RC%ThetaS(i, j, k) * Me%CV%ThetaHydroCoef, 0.0, Me%RC%ThetaS(i, j, k), 1.0, Me%Theta(i, j, k))

                        AccumPressure = AccumPressure + Coef
                        Me%HydroPressure(i,j,k) = AccumPressure
                        AccumPressure = AccumPressure + Coef
                        
                    else

                        Me%HydroPressure(i,j,k) = 0.0
                        AccumPressure = 0.0
                    
                    endif
                
                enddo
                
                !Final Head = Topography + Hydropressure + Soil Suction
                do k = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB
                    Me%FinalHead(i, j, k) = Me%ExtVar%CenterCell(i, j, k) + Me%HydroPressure(i,j,k) + Me%Head(i, j, k)
                enddo

            endif
            
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    end subroutine FinalHead

    !----------------------------------------------------------------------------
    
    subroutine InfiltrationVelocity
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, KUB
        real                                        :: hsup_aux
        
        KUB = Me%WorkSize%KUB
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints (i, j) == 1) then

                Me%UnsatVelWold(I,J,KUB+1) = Me%UnsatVelW(I,J,KUB+1)
                                    
                hsup_aux = Me%ExtVar%Topography(i, j) + Me%WaterColumn(i, j)
            
                Me%UnsatVelW(i, j, KUB+1) = -1.0 * BuckinghamDarcyEquation                      &
                                               (con      = Me%SatK (i, j, KUB),                 &
                                                hinf     = Me%FinalHead(i, j, KUB),             &
                                                hsup     = hsup_aux,                &
                                                delta    = Me%ExtVar%DWZ(i, j, KUB)/2.0)
 

                    
                                                                   
                if (Me%UnsatVelW(i, j, KUB+1) > Me%WaterColumn(i, j) / Me%CV%CurrentDT) then  
                    
                    Me%UnsatVelWold(I,J,KUB+1) = Me%UnsatVelW(I,J,KUB+1)

                    Me%UnsatVelW(i, j, KUB+1)  = Me%WaterColumn(i, j) / Me%CV%CurrentDT
               
                endif

            endif
            

        enddo
        enddo
    
    end subroutine InfiltrationVelocity

    !----------------------------------------------------------------------------

    subroutine SoilWaterFlux

        !Arguments---------------------------------------------------------------
                  
        !Local-------------------------------------------------------------------
        integer                                     :: I, J, K                
        integer                                     :: Chunk

        !------------------------------------------------------------------------          

        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "SoilWaterFlux")

        !$OMP PARALLEL PRIVATE(I,J,K)

        !Horizontal Velocities
        if (Me%SoilOpt%CalcHorizontal) then
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB            

                if (Me%ExtVar%ComputeFacesU3D(I,J,K) .EQ. Compute) then
                    Me%FluxU    (I,J,K) = Me%UnsatVelU(I,J,K) * Me%ExtVar%AreaU(I,J,K)
                end if

                if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then
                    Me%FluxV    (I,J,K) = Me%UnsatVelV(I,J,K) * Me%ExtVar%AreaV(I,J,K)
                end if

            end do
            end do
            end do
            !$OMP END DO NOWAIT
        endif
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB            
            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then
                if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then
                                        
                    Me%FluxW(i,j,k)      = Me%UnsatVelW(i,j,k) * Me%ExtVar%Area(i,j)
                    Me%FluxWFinal(i,j,k) = Me%FluxW(i,j,k)

                end if
            end if
        end do
        end do
        end do
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "SoilWaterFlux")

    end subroutine SoilWaterFlux
    
    !--------------------------------------------------------------------------
    
    subroutine Condutivity_Face                          

        !Local-----------------------------------------------------------------                        

        !----------------------------------------------------------------------        

        select case (Me%SoilOpt%CondutivityFace)

        case (Average )

            call CondutivityAverage

        case (Maximum )

            call CondutivityMaximum            

        case (Minimum )

            call CondutivityMinimum

        case (Weighted)

            write (*,*)'Not Implemented'
            stop 'Condutivity_Face - ModulePorousMedia - ERR01'

        case (GeometricAvg)

            call CondutivityGeometricAverage

        end select

    end subroutine Condutivity_Face

    !--------------------------------------------------------------------------

    subroutine CondutivityAverage

        !Local-----------------------------------------------------------------                        
        integer                                     :: i, j, k
        integer                                     :: Chunk        

        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        !$OMP PARALLEL PRIVATE(I,J,K)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesU3D(i, j, k) .EQ. Compute) then
                
                 Me%UnsatK_X(i, j, k) = (Me%UnSatK(i,j-1,k) * Me%ExtVar%DUX(i,j  )       +  &
                                         Me%UnSatK(i,j  ,k) * Me%ExtVar%DUX(i,j-1))      /  &
                                        (Me%ExtVar%DUX (i,j-1  ) + Me%ExtVar%DUX(i,j  ))      *  &
                                         Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then

                 Me%UnsatK_Y(i, j, k) = (Me%UnSatK(i-1,j,k) * Me%ExtVar%DVY(i  ,j)       +  &
                                         Me%UnSatK(i,j  ,k) * Me%ExtVar%DVY(i-1,j))      /  &
                                        (Me%ExtVar%DVY (i-1,j  ) + Me%ExtVar%DVY(i,j  ))      *  &
                                         Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Z
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then

                Me%UnsatK_Z(i, j, k) = (Me%UnSatK(i,j,k-1) * Me%ExtVar%DWZ(i,j,k  )      +  &
                                        Me%UnSatK(i,j,k)   * Me%ExtVar%DWZ(i,j,k-1))     /  &
                                       (Me%ExtVar%DWZ  (i,j,k-1) + Me%ExtVar%DWZ(i,j  ,k))

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

    end subroutine CondutivityAverage

    !--------------------------------------------------------------------------

    subroutine CondutivityMaximum

        !Local-----------------------------------------------------------------                        
        integer                                     :: i, j, k
        integer                                     :: Chunk        

        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        !X
        !$OMP PARALLEL PRIVATE(I,J,K)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesU3D(i, j, k) .EQ. Compute) then
                
                 Me%UnsatK_X(i, j, k) = max (Me%UnSatK(i,j-1,k), Me%UnSatK(i,j,k)) *  &
                                                Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT
        
        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then

                 Me%UnsatK_Y(i, j, k) = max (Me%UnSatK(i-1,j,k), Me%UnSatK(i,j,k)) *  &
                                                Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Z
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then


                Me%UnsatK_Z(i, j, k) = max (Me%UnSatK(i,j,k-1), Me%UnSatK(i,j,k))

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

    end subroutine CondutivityMaximum

    !--------------------------------------------------------------------------

    subroutine CondutivityMinimum

        !Local-----------------------------------------------------------------                        
        integer                                     :: i, j, k
        integer                                     :: Chunk        

        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        !X
        !$OMP PARALLEL PRIVATE(I,J,K)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesU3D(i, j, k) .EQ. Compute) then
                
                 Me%UnsatK_X(i, j, k) = min (Me%UnSatK(i,j-1,k), Me%UnSatK(i,j,k)) *  &
                                                Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then

                 Me%UnsatK_Y(i, j, k) = min (Me%UnSatK(i-1,j,k), Me%UnSatK(i,j,k)) *  &
                                                Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Z
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then
                Me%UnsatK_Z(i, j, k) = min (Me%UnSatK(i,j,k-1), Me%UnSatK(i,j,k))
            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

    end subroutine CondutivityMinimum

    !--------------------------------------------------------------------------

    subroutine CondutivityGeometricAverage

        !Local-----------------------------------------------------------------                        
        integer                                     :: i, j, k
        integer                                     :: Chunk        

        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        

        !X
        !$OMP PARALLEL PRIVATE(I,J,K)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesU3D(i, j, k) .EQ. Compute) then
                
                 Me%UnsatK_X(i, j, k) = (Me%UnSatK(i,j-1,k) ** Me%ExtVar%DUX(i,j  )     *  &
                                            Me%UnSatK(i,j  ,k) ** Me%ExtVar%DUX(i,j-1))    **  &
                                            (1.0/(Me%ExtVar%DUX (i,j-1) + Me%ExtVar%DUX(i,j)))   *  &
                                            Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then

                 Me%UnsatK_Y(i, j, k) = (Me%UnSatK(i-1,j,k) ** Me%ExtVar%DVY(i  ,j)       *   &
                                            Me%UnSatK(i,j  ,k) ** Me%ExtVar%DVY(i-1,j))      **  &
                                            (1.0/(Me%ExtVar%DVY (i-1,j) + Me%ExtVar%DVY(i,j)))     *   &
                                            Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Z
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then

                Me%UnsatK_Z(i, j, k) = (Me%UnSatK(i,j,k-1) ** Me%ExtVar%DWZ(i,j,k  )      *  &
                                           Me%UnSatK(i,j,k)   ** Me%ExtVar%DWZ(i,j,k-1))     ** &
                                          (1.0/(Me%ExtVar%DWZ  (i,j,k-1) + Me%ExtVar%DWZ(i,j  ,k)))

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

    end subroutine CondutivityGeometricAverage

    !--------------------------------------------------------------------------

    subroutine EvaporationFlux()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------                 
        integer                                     :: i, j, k 
        real                                        :: WaterVolume
        real                                        :: VelocityVolume
        real(8)                                     :: EvapoVolume, SoilVolume
        real                                        ::  NewTheta 
        real                                        :: TotalCol, HeadLimit
        !Begin-----------------------------------------------------------------

        !Set EvapFlux to zero
        Me%EvaporationFlux(:,:) = 0.0

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if (Me%ExtVar%WaterPoints3D(i,j,Me%WorkSize%KUB) == 1) then
        
                if (Me%WaterColumn(i, j) < AllmostZero) then             
                   
                    k = Me%WorkSize%KUB   ! evaporation just at the surface
                    
                    !m = m/s * s
                    TotalCol = Me%ExtVar%PotentialEvaporationFlux(i, j) * Me%CV%CurrentDT    ! available water for evaporation
                    ! m3               
                    WaterVolume = TotalCol * Me%ExtVar%Area (i,j)      

                    !Velocity Volume
                    if (Me%SoilOpt%LimitEVAPWaterVelocity) then
                        VelocityVolume = Me%UnSatK (i,j,k) * Me%CV%CurrentDT * Me%ExtVar%Area (i,j)
                        EvapoVolume    = min(WaterVolume, VelocityVolume)
                    else
                        EvapoVolume    = WaterVolume
                    endif
                            
                    !Avaliable Soil Water volume in layer
                    SoilVolume  = (Me%Theta(i,j,k) - Me%RC%ThetaR(i,j,k)) * Me%ExtVar%CellVolume(i,j,k) 
                    EvapoVolume = min(EvapoVolume, SoilVolume)
                            
                  
                    !Estimates new Theta
                    NewTheta = Me%Theta(i,j,k) - EvapoVolume/ Me%ExtVar%CellVolume(i,j,k)
                            
                    if (Me%SoilOpt%LimitEVAPHead) then     
                        HeadLimit= Me%Soilopt%HeadLimit
                       
                        if (NewTheta > Theta_(HeadLimit, Me%SoilID(i,j,k)))   then        
                                
                            !Evaporation Flux
                            Me%EvaporationFlux(i,j) = EvapoVolume / Me%CV%CurrentDT
                        else 
                    
                            Me%EvaporationFlux(i,j) = 0.0

                        endif
                    else
                        
                        !Just uses new theta if not all dry... stability reasons
                        if (NewTheta > Me%RC%ThetaR(i,j,k) + Me%CV%LimitThetaLo) then
                        
                            !Evaporation Flux
                           Me%EvaporationFlux(i,j) = EvapoVolume / Me%CV%CurrentDT

                        else 
                    
                            Me%EvaporationFlux(i,j) = 0.0

                        endif 
                    endif           
                
                endif
            
            endif
            
        enddo
        enddo            


    end subroutine EvaporationFlux

    !--------------------------------------------------------------------------

    subroutine CalculateNewTheta

        !Arguments-------------------------------------------------------------  

        !Local-----------------------------------------------------------------        
        integer                             :: CHUNK, I, J, K
        
        
        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "CalculateNewTheta")

        !$OMP PARALLEL PRIVATE(I,J,K)

        !Horizontal Fluxes
        if (Me%SoilOpt%CalcHorizontal) then

            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                    Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) +        &
                                        (Me%FluxU(i,j,k)  * Me%ExtVar%ComputeFacesU3D(i,j,k) -      &
                                         Me%FluxU(i,j+1,k)* Me%ExtVar%ComputeFacesU3D(i,j+1,k)) *   &
                                         Me%CV%CurrentDT) / Me%ExtVar%CellVolume(i, j, k)
                endif
            enddo
            enddo
            enddo
            !$OMP END DO
            
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                    Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) +        &
                                        (Me%FluxV(i,j,k)  * Me%ExtVar%ComputeFacesV3D(i,j,k) -      &
                                         Me%FluxV(i+1,j,k)* Me%ExtVar%ComputeFacesV3D(i+1,j,k)) *   &
                                         Me%CV%CurrentDT) / Me%ExtVar%CellVolume(i, j, k)
                endif
            enddo
            enddo
            enddo
            !$OMP END DO

        endif
        
        !Vertical Flux
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) +        &
                                     (Me%FluxW(i,j,k)  * Me%ExtVar%ComputeFacesW3D(i,j,k) -     &
                                      Me%FluxW(i,j,k+1)* Me%ExtVar%ComputeFacesW3D(i,j,k+1)) *  &
                                      Me%CV%CurrentDT) / Me%ExtVar%CellVolume(i, j, k)
            endif
        enddo
        enddo
        enddo
        !$OMP END DO
        
        !Evapotranspiration
        
        if (Me%TranspirationExists) then
            !Transpiration
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                    Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) -        & 
                                         Me%ExtVar%TranspirationFlux(i, j, k) * Me%CV%CurrentDT)  / &
                                         Me%ExtVar%CellVolume(i,j,k)
                endif
            enddo
            enddo
            enddo
            !$OMP END DO
        endif

        if (Me%EvaporationExists) then
            !Evaporation
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%BasinPoints(I,J) == WaterPoint) then
                    k = Me%WorkSize%KUB
                    Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) -       & 
                                         Me%EvaporationFlux(i, j) * Me%CV%CurrentDT)  /            &
                                         Me%ExtVar%CellVolume(i,j,k)
                endif
            enddo
            enddo
            !$OMP END DO
        endif
        
        !Infiltration
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i,j) == WaterPoint) then
                k = Me%WorkSize%KUB
                Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) +                &
                                     Me%UnsatVelW(i, j, k+1) * Me%ExtVar%Area(i, j) *                   &
                                     (1.0 - Me%ImpermeableFraction(i, j)) * Me%CV%CurrentDT) /             &
                                     Me%ExtVar%CellVolume(i, j, k)
            endif
        enddo
        enddo
        !$OMP END DO
        
        
        !Exchange with River
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%RiverPoints(i,j) == WaterPoint) then
                k = Me%UGCell(i,j)
                Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) -                &
                                     Me%lFlowToChannels(i, j) * Me%CV%CurrentDT) /                      &
                                     Me%ExtVar%CellVolume(i, j, k)
            endif
        enddo
        enddo
        !$OMP END DO
        
        !$OMP END PARALLEL
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "CalculateNewTheta")
        

    end subroutine CalculateNewTheta

    !--------------------------------------------------------------------------

    subroutine IntegrateValuesInTime(SumDT)
    
        !Arguments-------------------------------------------------------------  
        real                                :: SumDT

        !Local-----------------------------------------------------------------        

        !Updates water column and infiltration
        call UpdateWaterColumnInfiltration
        
        !Updates Efective EVTP
        call UpdateEfectiveEVTP
        
        !Integrates Flow
        call IntegrateFlow (SumDT)

    end subroutine IntegrateValuesInTime

    !--------------------------------------------------------------------------
    ! esta função é para desaparecer do porous Media porque passou para o PorousMediaProperties
    subroutine CalculateAdvection
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, k, CHUNK
        real                                        :: OldMass, NewMass
        
        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)

        !$OMP PARALLEL PRIVATE(I,J,K, OldMass, NewMass)
        
        CurrProperty => Me%FirstProperty
        do while (associated(CurrProperty))
        
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                
                    !Old Mass
                    OldMass = Me%CV%ThetaOld(i, j, k) * Me%ExtVar%CellVolume(i, j, k) * &
                              CurrProperty%ConcentrationOld(i, j, k)

                    !New Mass                              
                    NewMass = OldMass + Me%FluxU(i, j,   k) * CurrProperty%ConcentrationOld(i, j-1, k) - &
                                        Me%FluxU(i, j+1, k) * CurrProperty%ConcentrationOld(i, j+1, k) + &
                                        Me%FluxV(i, j,   k) * CurrProperty%ConcentrationOld(i-1, j, k) - &
                                        Me%FluxV(i+1, j, k) * CurrProperty%ConcentrationOld(i+1, j, k) + &
                                        Me%FluxW(i, j,   k) * CurrProperty%ConcentrationOld(i, j, k-1) - &
                                        Me%FluxW(i, j, k+1) * CurrProperty%ConcentrationOld(i, j, k+1)
                                        
                    !UpperBoundary
                    NewMass = NewMass + Me%UnsatVelW(i, j, Me%WorkSize%KUB+1) * (1.0 - Me%ImpermeableFraction(i, j)) *    &
                                        Me%ExtVar%Area(i, j) * CurrProperty%UpperConcentration(i, j)
                    
                    !New Concentration
                    CurrProperty%Concentration(i, j, k) = NewMass / (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k))
                    
                    
                
                endif
            enddo
            enddo
            enddo
            !$OMP END DO
           
        
            CurrProperty => CurrProperty%Next
        enddo

        !$OMP END PARALLEL

    end subroutine CalculateAdvection
    
    !--------------------------------------------------------------------------

    subroutine UpdateWaterColumnInfiltration 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: chunk
        real                                        :: dh
     
        !Begin-----------------------------------------------------------------
        
        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J, dh)

        !Update Water Column
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB,     Me%WorkSize%JUB
        do I = Me%WorkSize%ILB,     Me%WorkSize%IUB
                                    
            if (Me%ExtVar%BasinPoints (i,j) == 1) then
                
                !Variation in height
                dh                  = Me%UnsatVelW(i, j, Me%WorkSize%KUB+1) * (1.0 - Me%ImpermeableFraction(i, j)) * &
                                      Me%CV%CurrentDT
                
                !Just reduce water column due to Infiltration (PermeableFraction only)
                Me%WaterColumn(i,j) = Me%WaterColumn(i,j)  - dh
                
                Me%Infiltration(i,j)= Me%Infiltration(i,j) + dh
                 
                if (abs(Me%WaterColumn(i, j)) < AllmostZero) Me%WaterColumn(i, j) = 0.0
                
                if (Me%WaterColumn(i, j) < 0.0) then
                    write(*,*)'Bug Infiltration', i, j, Me%WaterColumn(i, j)
                endif
                
            endif
            
        enddo
        enddo            
        !$OMP END DO
        !$OMP END PARALLEL


    end subroutine UpdateWaterColumnInfiltration

    !--------------------------------------------------------------------------

    subroutine UpdateEfectiveEVTP

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        integer                                     :: chunk
     
        !Begin-----------------------------------------------------------------
        
        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)

        if (Me%TranspirationExists) then
            !$OMP PARALLEL SHARED(CHUNK) PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB,     Me%WorkSize%KUB
            do J = Me%WorkSize%JLB,     Me%WorkSize%JUB
            do I = Me%WorkSize%ILB,     Me%WorkSize%IUB
                                
                if (Me%ExtVar%Waterpoints3D (i,j,k) == 1) then
          
                    Me%EfectiveEVTP(i,j) = Me%EfectiveEVTP(i,j) + Me%ExtVar%TranspirationFlux(i, j, k) * Me%CV%CurrentDT/ &
                                           Me%ExtVar%Area(i, j)
                endif

            enddo
            enddo            
            enddo  
              
            !$OMP END DO
            !$OMP END PARALLEL     
        
        endif
        
        if (Me%EvaporationExists) then
            !!!$OMP PARALLEL SHARED(CHUNK) PRIVATE(I,J,K)
            !!!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do J = Me%WorkSize%JLB,     Me%WorkSize%JUB
            do I = Me%WorkSize%ILB,     Me%WorkSize%IUB
                                
                if (Me%ExtVar%BasinPoints (i,j) == 1) then
          
                    Me%EfectiveEVTP(i,j) = Me%EfectiveEVTP(i,j) + Me%EvaporationFlux(i, j) * Me%CV%CurrentDT / &
                                           Me%ExtVar%Area(i, j)
                endif
        
            enddo
            enddo            
              
            !!!$OMP END DO
            !!!$OMP END PARALLEL     
        
        endif
        

    end subroutine UpdateEfectiveEVTP

    !--------------------------------------------------------------------------
    
    subroutine IntegrateFlow(SumDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: SumDT

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: chunk
     
        CHUNK = CHUNK_J(Me%WorkSize%KLB, Me%WorkSize%JUB)

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%RiverPoints (i,j) == 1) then
                Me%iFlowToChannels(i, j) = (Me%iFlowToChannels(i, j) * SumDT +                  &
                                            Me%lFlowToChannels(i, j) * Me%CV%CurrentDT) /       &
                                           (SumDT + Me%CV%CurrentDT)
            endif
            
        enddo
        enddo

    end subroutine IntegrateFlow

    !--------------------------------------------------------------------------

    subroutine SoilParameters (Mapping)

        !Arguments-------------------------------------------------------------        
        integer, dimension(:,:), pointer            :: Mapping

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        integer                                     :: Chunk

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "SoilParameters")
                
        !$OMP PARALLEL PRIVATE(I,J,K)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

cd1 :   if (Mapping(i, j) == 1) then

            do k = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB

            !Verifies if cell is saturated (or very close). If so set it to
            !saturation
            !Otherwise calculate ThetaF Value                           !0.00001
            if (abs(Me%Theta(i,j,k) - Me%RC%ThetaS(i,j,k)) < Me%CV%LimitThetaHi) then

                Me%Theta      (i, j, k)         = Me%RC%ThetaS(i,j,k)
                Me%Head       (i, j, k)         = 0.0 
                Me%RC%ThetaF  (i, j, k)         = 1.0
                Me%UnSatK     (i, j, k)         = Me%SatK(i, j, k)
                Me%CalculateHead(i, j, k)       = .false.

            !Close or below residual value
            else if (Me%Theta(i,j,k) < Me%RC%ThetaR(i,j,k) + Me%CV%LimitThetaLo) then
            
                !This creates mass...
                Me%RC%ThetaF    (i, j, k) = Me%CV%LimitThetaLo / (Me%SoilTypes(Me%SoilID(I,J,K))%ThetaS - Me%SoilTypes(Me%SoilID(I,J,K))%ThetaR)
                Me%CalculateHead(i, j, k) = .true.
                
                call SetError(WARNING_, INTERNAL_, "Mass Created, SoilParameters", OFF)

            !Normal Case
            else
                Me%RC%ThetaF  (i, j, k)   = ThetaF_ (Me%Theta (i, j, k), Me%SoilID(i, j, k))
                Me%CalculateHead(i, j, k) = .true.
            endif
            
            end do

        end if cd1
        
        end do
        end do
        !$OMP END DO


        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

cd2 :   if (Mapping(i, j) == 1) then

            do k = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB

            if (Me%CalculateHead(i, j, k)) then

                !Over saturation
                if (Me%RC%ThetaF(i, j, k) > 1.0) then

                    Me%Head   (i, j, k) = Me%SoilTypes(Me%SoilID(I,J,K))%OverSatSlope * (Me%Theta (i, j, k) - Me%RC%ThetaS (i, j, k))
                    Me%UnSatK (i, j, k) = Me%SatK(i, j, k)

                !0 < Theta < 1
                else if (Me%RC%ThetaF(i, j, k) > 0.0 .and. Me%RC%ThetaF(i, j, k) < 1.0) then

                    Me%Head  (i, j, k) = Head_   (Me%RC%ThetaF (i, j, k), Me%SoilID(i, j, k))

                    Me%UnSatK(i, j, k) = UnsatK_ (Me%RC%ThetaF (i, j, k), Me%SoilID(i, j, k))



                !Theta <= 0
                else

                    Me%Head(I,J,K) = null_real

                endif
                
            endif

            end do

        end if cd2
        
        end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
            

        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "SoilParameters")

    end subroutine SoilParameters

    !--------------------------------------------------------------------------    
    
    subroutine CheckCompleteStupidSolution(IsStupid, CorrectTheta)

  
        !Arguments-------------------------------------------------------------
        logical, intent(OUT)                        :: IsStupid        
        logical, intent(IN)                         :: CorrectTheta

        !Local-----------------------------------------------------------------        
        integer                                     :: I, J, K
        integer                                     :: nStupid

        !----------------------------------------------------------------------               

        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "CheckCompleteStupidSolution")

        nStupid = 0

        !Tests variations
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(I,J,K) == 1) then
                if (Me%Theta(i, j, k) < Me%RC%ThetaR(i, j, k)) then
                    nStupid = nStupid + 1
                endif
            endif
        enddo
        enddo
        enddo    
        
        if (nStupid > 0) then
            IsStupid = .true.
        else
            IsStupid = .false.
        endif
        
        if (CorrectTheta .and. IsStupid) then
            IsStupid = .false.
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%WaterPoints3D(I,J,K) == 1) then
                    if (Me%Theta(i, j, k) < Me%RC%ThetaR(i, j, k)) then
                        Me%Theta(i, j, k) = Me%RC%ThetaR(i, j, k)
                        call SetError(WARNING_, INTERNAL_, "Mass Created, STUPID SOLUTION", OFF)
                    endif
                endif
            enddo
            enddo
            enddo    
        endif
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "CheckCompleteStupidSolution")
        
    end subroutine CheckCompleteStupidSolution        

    !--------------------------------------------------------------------------    

    subroutine variation_test( StrongVariation)
        
        !Arguments-------------------------------------------------------------
        logical, intent(OUT)                        :: StrongVariation        

        !Local-----------------------------------------------------------------        
        integer                                     :: I, J, K
        integer                                     :: nStrongs

        !----------------------------------------------------------------------               

        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "variation_test")

        nStrongs = 0

        !Tests variations
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(I,J) == 1) then
            
                !In the upper layer the module can't check the convergence criteria because it is imposing
                !a boundary value... 
                !If the water column isn't suficient, DT for infiltration velocity and DT to converge are cutted
                !in the same way.... 
                do K = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB - 1
                    if (abs(Me%CV%ThetaOld(I,J,K)-Me%Theta(I,J,K)) > Me%CV%ThetaTolerance) then
                        nStrongs = nStrongs + 1
                    endif
                enddo
            endif
        enddo
        enddo
        
        if (nStrongs > 0) then
            StrongVariation = .true.
        else
            StrongVariation = .false.
        endif
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "variation_test")

    end subroutine variation_test
    
    !--------------------------------------------------------------------------

    subroutine VerticalContinuity
    
        !Arguments-------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        integer                                     :: chunk
        real                                        :: ExcessVolume, AvaliableVolume, dh
        
        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J,K, ExcessVolume, AvaliableVolume, dh)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(I,J) == 1) then
            
                do k = Me%WorkSize%KUB, (Me%ExtVar%KFloor(i, j)+1), -1               
                    if (Me%Theta(i,j,k) .gt. Me%RC%ThetaS(i,j,k)) then
                        !If cell is oversaturated, set to saturation and put water downwards
                        ExcessVolume         = (Me%Theta(i,j,k) - Me%RC%ThetaS (i,j,k)) * Me%ExtVar%CellVolume(i,j,k)
                        Me%Theta(i,j,k-1)    = ((Me%Theta(i,j,k-1) * Me%ExtVar%CellVolume(i,j,k-1)) + ExcessVolume) / Me%ExtVar%CellVolume(i,j,k-1)
                        Me%FluxWFinal(i,j,k) = Me%FluxWFinal(i,j,k) + (-ExcessVolume / Me%CV%CurrentDT)                       
                        Me%Theta(i,j,k)      = Me%RC%ThetaS (i,j,k)
                    endif
                enddo
                
                !Invert process and put the rest upwards
                do k = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB-1
                    if (Me%Theta(i,j,k) .gt. Me%RC%ThetaS(i,j,k)) then
                        !If cell is oversaturated, set to saturation and put water upwards
                        ExcessVolume           = (Me%Theta(i,j,k) - Me%RC%ThetaS (i,j,k)) * Me%ExtVar%CellVolume(i,j,k)
                        Me%Theta(i,j,k+1)      = ((Me%Theta(i,j,k+1) * Me%ExtVar%CellVolume(i,j,k+1)) + ExcessVolume) / Me%ExtVar%CellVolume(i,j,k+1)
                        Me%FluxWFinal(i,j,k+1) = Me%FluxWFinal(i,j,k+1) + (ExcessVolume / Me%CV%CurrentDT)                       
                        Me%Theta(i,j,k)        = Me%RC%ThetaS (i,j,k)
                    endif
                enddo
                
                
                !Put remaing volume to watercolumn
                k = Me%WorkSize%KUB
                if (Me%Theta(i,j,k) .gt. Me%RC%ThetaS(i,j,k)) then
                    !If cell is oversaturated, set to saturation and put water on water column
                    ExcessVolume           = ((Me%Theta(i,j,k) - Me%RC%ThetaS (i,j,k)) * Me%ExtVar%CellVolume(i,j,k))
                    Me%FluxWFinal(i,j,k+1) = Me%FluxWFinal(i,j,k+1) + (ExcessVolume / Me%CV%CurrentDT)                                           
                    dh                     = ExcessVolume  / Me%ExtVar%Area(i, j)
                    Me%WaterColumn  (i,j)  = Me%WaterColumn(i,j) + dh
                    Me%Infiltration (i,j)  = Me%Infiltration (i,j) - dh
                    Me%Theta(i,j,k)        = Me%RC%ThetaS (i,j,k)
                endif
            endif
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    end subroutine VerticalContinuity

    !--------------------------------------------------------------------------

    subroutine LogDT (iteration)

        !Arguments-------------------------------------------------------------
        integer                                     :: iteration

        !Local-----------------------------------------------------------------
        real                                        :: Year, Month, Day
        real                                        :: Hour, Minute, Second

        call ExtractDate(Me%ExtVar%Now, Year, Month, Day, Hour, Minute, Second)
        
        write (Me%Files%AsciiUnit, fmt=1000) Year, Month, Day, Hour, Minute, Second, &
                                             iteration, Me%CV%CurrentDT

        1000 format(f5.0, f5.0, f5.0, f5.0, f5.0, f12.5, i3, f12.5)

    end subroutine LogDT

    !--------------------------------------------------------------------------

    subroutine CalculateUGWaterLevel
    
        !Arguments-------------------------------------------------------------
            
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, CHUNK
        real                                        :: CellBottomLevel, DZInCell
        real                                        :: FieldHead, FieldTheta

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J,K, CellBottomLevel, DZInCell, FieldHead, FieldTheta)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(I,J) == 1) then

                !Find 1. non saturated cell            
doK:            do K = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB
                    if (Me%Theta(i,j,k) < Me%RC%ThetaS(i, j, k) - Me%CV%LimitThetaHi) then
                        exit doK
                    endif
                enddo doK
                
                k = min(k, Me%WorkSize%KUB)
                
                CellBottomLevel = - Me%ExtVar%SZZ(i,j,k-1)
                FieldHead       = - 0.5 * Me%ExtVar%DWZ(i, j, k)
                FieldTheta      = Theta_(FieldHead, Me%SoilID(i,j,k))
                
                if (Me%Theta(i,j,k) < FieldTheta) then
                    DZInCell        = 0.0
                else
                    DZInCell        = LinearInterpolation(FieldTheta, 0.0, Me%RC%ThetaS(i,j,k), &
                                                          Me%ExtVar%DWZ(i,j,k), Me%Theta(i,j,k))
                endif
                
                Me%UGWaterLevel2D(i, j) = CellBottomLevel + DZInCell
                Me%UGWaterDepth2D(i, j) = Me%ExtVar%Topography(i, j) -Me%UGWaterLevel2D(i, j)
                Me%UGCell        (i, j) = k
                
            endif
            
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
                
    end subroutine CalculateUGWaterLevel

    !--------------------------------------------------------------------------

    real function Head_(thf, SoilID)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: thf
        integer, intent(IN)                         :: SoilID

        !Local-------------------------------------------------------------------

        Head_= - abs( ( ( thf**( -1.0 / Me%SoilTypes(SoilID)%MFit)) -1.0 )**( 1.0 / Me%SoilTypes(SoilID)%NFit ) ) / &
                Me%SoilTypes(SoilID)%Alfa

    end function Head_   

    !--------------------------------------------------------------------------

    real function Head_slope_(th, SoilID)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: th
        integer, intent(IN)                         :: SoilID
        
        !Local-------------------------------------------------------------------
        real                                        :: thf

        !------------------------------------------------------------------------

        if (th < Me%SoilTypes(SoilID)%ThetaS - Me%CV%LimitThetaHi) then

            if (th > Me%SoilTypes(SoilID)%ThetaR) then
                thf = ThetaF_(th, SoilID)
            else 
                thf = Me%CV%LimitThetaLo / (Me%SoilTypes(SoilID)%ThetaS - Me%SoilTypes(SoilID)%ThetaR)
            endif

            Head_slope_= (1./Me%SoilTypes(SoilID)%Alfa/Me%SoilTypes(SoilID)%NFit) * &
                        (thf**(-1./Me%SoilTypes(SoilID)%MFit)-1.) ** ((1.-Me%SoilTypes(SoilID)%NFit) / &
                        Me%SoilTypes(SoilID)%NFit) * (1./Me%SoilTypes(SoilID)%MFit) * &
                        thf ** (-(1.+ Me%SoilTypes(SoilID)%MFit) / Me%SoilTypes(SoilID)%MFit) * 1. /   &
                        ( Me%SoilTypes(SoilID)%ThetaS - Me%SoilTypes(SoilID)%ThetaR)
            return
        else
            
            Head_slope_ = Me%SoilTypes(SoilID)%OverSatSlope
            return

        end if
            

    end function Head_slope_
    
    !--------------------------------------------------------------------------

    real function Theta_ (head, SoilID)

        !Arguments-------------------------------------------------------------------
        real, intent(IN)                            :: head
        integer, intent(IN)                         :: SoilID

        !----------------------------------------------------------------------

        Theta_ = Me%SoilTypes(SoilID)%ThetaR + ( Me%SoilTypes(SoilID)%ThetaS - Me%SoilTypes(SoilID)%ThetaR ) / &
                ((1.0 + (Me%SoilTypes(SoilID)%alfa * (- head) ) ** Me%SoilTypes(SoilID)%Nfit ) **              &
                Me%SoilTypes(SoilID)%MFit )

    end function Theta_

    !----------------------------------------------------------------------------
    
    real function ThetaF_ (th, SoilID)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: th
        integer, intent(IN)                         :: SoilID
        
        ThetaF_ = (th - Me%SoilTypes(SoilID)%ThetaR) / (Me%SoilTypes(SoilID)%ThetaS - Me%SoilTypes(SoilID)%ThetaR)     

    end function ThetaF_ 

    !----------------------------------------------------------------------------

    real function UnsatK_(thf, SoilID)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: thf
        integer, intent(IN)                         :: SoilID

        !Local-----------------------------------------------------------------

        UnsatK_ = Me%SoilTypes(SoilID)%SatK*((thf)**Me%SoilTypes(SoilID)%lfit) * &
                (1.0d0-(1.0d0-thf**(1.0d0/Me%SoilTypes(SoilID)%MFit))**Me%SoilTypes(SoilID)%MFit)**2.0d0
                    
    end function UnsatK_

    !--------------------------------------------------------------------------

    subroutine ExchangeWithDrainageNetwork2

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, STAT_CALL
        real                                        :: Area
        real,   dimension(:, :), pointer            :: ChannelsWaterLevel
        real,   dimension(:, :), pointer            :: ChannelsBottomLevel 
        real,   dimension(:, :), pointer            :: ChannelsBottomWidth
        real,   dimension(:, :), pointer            :: ChannelsNodeLength
        integer,dimension(:, :), pointer            :: ChannelsOpenProcess

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR01'

        call GetChannelsBottomLevel (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR02'

        call GetChannelsBottomWidth (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR03'

        call GetChannelsOpenProcess (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR04'

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR05'

        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%RiverPoints(i, j) == OpenPoint) then
            
                !Channel will loose water
                if (ChannelsWaterLevel(i, j) > Me%UGWaterLevel2D(i, j)) then

                    
                    if (ChannelsOpenProcess(i, j) == 1) then

                        !The Area is always two times the lateral area and
                        !one time the bottom area
                        Area  = (2. * (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j))  & 
                                     + ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)
                                     
                        !Positive Flow to channel when Channel "gains" water
                        ![m3/s]                    [m2]   [m/s]                     
                        Me%lFlowToChannels(i, j) = -1.0 * Area * Me%UnSatK(i, j, Me%UGCell(i,j))

                        if (-1. * Me%lFlowToChannels(i, j) * Me%ExtVar%DT > (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j)) * Area) then
                            Me%lFlowToChannels(i, j) = -0.5 * (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j)) * Area / Me%ExtVar%DT
                            call SetError(WARNING_, INTERNAL_, "Flow to channel corrected. FLOW TO CHANNEL", OFF)
                        endif
                        
                    else
                    
                        Me%lFlowToChannels(i, j) = 0.0
                    
                    endif
                        
                    
                else
            
                    Area = (2. * (Me%UGWaterLevel2D(i, j) - ChannelsBottomLevel(i, j)) + &
                            ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)

                    !Positive Flow to channel when Channel "gains" water
                    ![m3/s]                    [m2]   [m/s]                     
                    Me%lFlowToChannels(i, j) = Area * Me%UnSatK(i, j, Me%UGCell(i,j))
                    
                endif

            else

                Me%lFlowToChannels(i, j) = 0.0

            endif

        enddo
        enddo
        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR07'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR08'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR09'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR10'

    end subroutine ExchangeWithDrainageNetwork2
    
    !--------------------------------------------------------------------------

    subroutine ExchangeWithDrainageNetwork

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, STAT_CALL
        real                                        :: dH, Area
        real                                        :: MinHead, FieldTheta, MaxVolume
        real,   dimension(:, :), pointer            :: ChannelsWaterLevel
        real,   dimension(:, :), pointer            :: ChannelsBottomLevel 
        real,   dimension(:, :), pointer            :: ChannelsBottomWidth
        real,   dimension(:, :), pointer            :: ChannelsNodeLength
        integer,dimension(:, :), pointer            :: ChannelsOpenProcess

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR01'

        call GetChannelsBottomLevel (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR02'

        call GetChannelsBottomWidth (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR03'

        call GetChannelsOpenProcess (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR04'

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR05'

        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%RiverPoints(i, j) == OpenPoint) then
            
            
                !Channel will loose water
                if (ChannelsWaterLevel(i, j) > Me%UGWaterLevel2D(i, j)) then

                    !The Area is always two times the lateral area and
                    !one time the bottom area
                    Area  = (2. * (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j))  & 
                                 + ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j) 

                    if (Me%UGWaterLevel2D(i, j) > ChannelsBottomLevel(i, j)) then

                        !Negative dh -> Flux from channels to porous media
                        dH    = Me%UGWaterLevel2D(i, j) - ChannelsWaterLevel(i, j) 

                    else

                        !Negative dh -> Flux from channels to porous media
                        dH    = ChannelsBottomLevel(i, j) - ChannelsWaterLevel(i, j)

                    endif

                else
            
                    !REVER AQUI
                    if (Me%UGWaterLevel2D(i, j) - Me%ExtVar%BottomTopoG(i, j) > MinUGThickness) then

                        !Postive dh -> Flux from porous media to channels
                        !This logical has to be rechecked - What happens if channels are deeper then soil layer??????
                        dH   = Me%UGWaterLevel2D(i, j) - max(ChannelsWaterLevel(i, j), Me%ExtVar%BottomTopoG(i, j))
                        !2* Lateral Area + Bottom Area 
                        Area = (2. * (Me%UGWaterLevel2D(i, j) - ChannelsBottomLevel(i, j)) + &
                                ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)

                    else

                        dH   = 0.0
                        Area = 0.0

                    endif

                endif

                !Positive Flow to channel when Channel "gains" water
                ![m3/s]                  [m]  [m2]   [m/s]                     [m]
                Me%lFlowToChannels(i, j) = dH * Area * Me%SoilOpt%HCondFactor * Me%UnSatK(i, j, Me%UGCell(i,j)) / &
                                           max(ChannelsBottomWidth(i, j) / 2.0, 1.0)
                !Me%lFlowToChannels(i, j) = Area * Me%SatK(i, j, Me%UGCell(i,j))

                !If the channel looses water (infiltration), then set max flux so that volume in channel does not get 
                !negative
                if (dH < 0) then
                    if (-1. * Me%lFlowToChannels(i, j) * Me%ExtVar%DT > (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j)) * Area) then
                        Me%lFlowToChannels(i, j) = -0.5 * (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j)) * Area / Me%ExtVar%DT
                        write(*,*)'FlowToChannels corrected - ModulePorousMedia'
                    endif
                endif
                
                !If soil looses water set flow so that cell stays at least with field theta
                !THIS DOESN'T WORK
                if (dH > 0) then
                    k           = Me%UGCell(i,j)
                    MinHead     = -0.5 * Me%ExtVar%DWZ(i, j, k)
                    FieldTheta  = Theta_(MinHead, Me%SoilID(i,j,k))
                    MaxVolume   = (Me%RC%ThetaS (i,j,k) - Me%Theta(i,j,k)) * Me%ExtVar%CellVolume(i,j,k)
                    Me%lFlowToChannels(i,j) = min(Me%lFlowToChannels(i,j), MaxVolume / Me%ExtVar%DT)
                endif
                

            else

                Me%lFlowToChannels(i, j) = 0.0

            endif

        enddo
        enddo
        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR07'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR08'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR09'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR10'

    end subroutine ExchangeWithDrainageNetwork

    !--------------------------------------------------------------------------

    function BuckinghamDarcyEquation(con, hinf, hsup, delta)
    real ::  BuckinghamDarcyEquation

        !Local-------------------------------------------------------------------
        real, intent(IN) :: con
        real, intent(IN) :: hinf,  hsup
        real, intent(IN) :: delta

        !------------------------------------------------------------------------

        BuckinghamDarcyEquation = -1.0 * con * ((hsup - hinf) / delta)
                       
        !------------------------------------------------------------------------

    end function BuckinghamDarcyEquation

    !----------------------------------------------------------------------------

    subroutine PorousMediaOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: KLB, KUB, k
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePointer       
        real, dimension(:,:,:), pointer             :: Modulus3D
        real, dimension(:,:,:), pointer             :: CenterU3D, CenterV3D, CenterW3D
        real, dimension(:,:), pointer               :: SurfaceSlice

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        !NOTE: This test below is only valid if either normal output or 2D is active. 
        !That way during startup it is checked if only one of the output types is activated
        !Both can't ba active at the same time beause of the TIME output

        if (Me%OutPut%Yes) then
            
            !3D "Normal" Output
            if (Me%ExtVar%Now >= Me%OutPut%OutTime(Me%OutPut%NextOutPut)) then

                !----------------------------------------------------------------------------
                !Common Output---------------------------------------------------------------
                !----------------------------------------------------------------------------

                !Time
                call ExtractDate   (Me%ExtVar%Now , AuxTime(1), AuxTime(2),         &
                                                    AuxTime(3), AuxTime(4),         &  
                                                    AuxTime(5), AuxTime(6))
                TimePointer => AuxTime

                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)            
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR01'


                call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time",                   &
                                     "YYYY/MM/DD HH:MM:SS",                         &
                                     Array1D      = TimePointer,                    &
                                     OutputNumber = Me%OutPut%NextOutPut,           &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR02'


                !Limits 
                call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, KLB-1, KUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR03'


                !Vertical 
                call HDF5WriteData  ( Me%ObjHDF5,  "/Grid/VerticalZ",               & 
                                     "Vertical",   "m"              ,               & 
                                      Array3D      = Me%ExtVar%SZZ  ,               &
                                      OutputNumber = Me%OutPut%NextOutPut,          &
                                      STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR01'


                !Limits 
                call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR03'

                !Open Points
                call HDF5WriteData   ( Me%ObjHDF5,  "/Grid/OpenPoints" ,            &
                                      "OpenPoints", "-"                 ,           &
                                       Array3D      = Me%ExtVar%OpenPoints3D ,      &
                                       OutputNumber = Me%OutPut%NextOutPut,         &
                                       STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR04'


                !----------------------------------------------------------------------------
                !Saturated Output------------------------------------------------------------
                !----------------------------------------------------------------------------
                
                !UGWaterLevel
                call HDF5WriteData   (Me%ObjHDF5, "/Results/WaterLevel",            &
                                      "WaterLevel", "m",                            &
                                      Array2D      = Me%UGWaterLevel2D,             &
                                      OutputNumber = Me%OutPut%NextOutPut,          &
                                      STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR05'

                !UGWaterDepth
                call HDF5WriteData   (Me%ObjHDF5, "/Results/water table depth",     &
                                      "water table depth", "m",                     &
                                      Array2D      = Me%UGWaterDepth2D,             &
                                      OutputNumber = Me%OutPut%NextOutPut,          &
                                      STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR05'

                !----------------------------------------------------------------------------
                !Unsaturated Output----------------------------------------------------------
                !----------------------------------------------------------------------------

                call HDF5WriteData ( Me%ObjHDF5,    "/Results/water content",       &
                                    "water content", "m3water/m3soil"            ,  &
                                     Array3D      =  Me%Theta             ,         &
                                     OutputNumber =  Me%OutPut%NextOutPut    ,      &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR01'
                
                call HDF5WriteData ( Me%ObjHDF5,    "/Results/relative water content",     &
                                    "relative water content", "m3water/m3water"        ,   &
                                     Array3D      =  Me%RC%ThetaF         ,         &
                                     OutputNumber =  Me%OutPut%NextOutPut    ,      &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR01'

                call HDF5WriteData ( Me%ObjHDF5, "/Results/Head"            ,       &
                                    "Head"     , 'm'                            ,   &
                                     Array3D      = Me%Head               ,         &
                                     OutputNumber = Me%OutPut%NextOutPut        ,   & 
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR01'

                call HDF5WriteData ( Me%ObjHDF5, "/Results/FinalHead"    ,          &
                                    "FinalHead"     , 'm'                       ,   &
                                     Array3D      = Me%FinalHead,                   &
                                     OutputNumber = Me%OutPut%NextOutPut,           & 
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR01'

                call HDF5WriteData ( Me%ObjHDF5, "/Results/HydroPressure"    ,   &
                                    "HydroPressure"     , 'm'                            ,   &
                                     Array3D      = Me%HydroPressure,   &
                                     OutputNumber = Me%OutPut%NextOutPut        ,   & 
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR01'

                !Write unsat velocities
                allocate(CenterU3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB), &
                         CenterV3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB), &
                         CenterW3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB), &
                         Modulus3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
                CenterU3D  = 0.0
                CenterV3D  = 0.0
                CenterW3D  = 0.0
                Modulus3D  = 0.0
                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB
                    if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                        CenterU3D (i, j, k) = (Me%UnsatVelU(i, j, k) + Me%UnsatVelU(i, j+1, k)) / 2.0
                        CenterV3D (i, j, k) = (Me%UnsatVelV(i, j, k) + Me%UnsatVelV(i+1, j, k)) / 2.0
                        CenterW3D (i, j, k) = (Me%UnsatVelW(i, j, k) + Me%UnsatVelW(i, j, k+1)) / 2.0
                        Modulus3D (i, j, k) = sqrt(CenterU3D (i, j, k)**2.0 + CenterV3D (i, j, k)**2.0 + CenterW3D(i,j,k)**2.0)
                    endif
                enddo
                enddo            
                enddo

                call HDF5WriteData ( Me%ObjHDF5, "/Results/Velocity/X",         & 
                                    "Vel. X",    'm/s',                             &
                                     Array3D      = CenterU3D,                      &
                                     OutputNumber = Me%OutPut%NextOutPut,           &  
                                     STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR01'
            
                call HDF5WriteData ( Me%ObjHDF5, "/Results/Velocity/Y",         & 
                                    "Vel. Y",    'm/s',                             &
                                     Array3D      = CenterV3D,                      &
                                     OutputNumber = Me%OutPut%NextOutPut,           &  
                                     STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR01'

                call HDF5WriteData ( Me%ObjHDF5, "/Results/Velocity/Z",         & 
                                    "Vel. Z",    'm/s',                             &
                                     Array3D      = CenterW3D,                      &
                                     OutputNumber = Me%OutPut%NextOutPut,           &  
                                     STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR01'

                call HDF5WriteData ( Me%ObjHDF5, "/Results/Velocity/Modulus",   & 
                                    "Modulus",    'm/s',                            &
                                     Array3D      = Modulus3D,                      &
                                     OutputNumber = Me%OutPut%NextOutPut,           &  
                                     STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR01'


                deallocate (CenterU3D, CenterV3D, CenterW3D, Modulus3D)

                Me%OutPut%NextOutPut = Me%OutPut%NextOutPut + 1

                !----------------------------------------------------------------------------
                !Write everything to disk----------------------------------------------------
                !----------------------------------------------------------------------------

                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR09'


            endif
            
        endif
        
        !2D Surface output
        if (Me%OutPut%SurfaceOutput) then
            
            if (Me%ExtVar%Now >= Me%OutPut%SurfaceOutTime(Me%OutPut%NextSurfaceOutput)) then


                !----------------------------------------------------------------------------
                !Common Output---------------------------------------------------------------
                !----------------------------------------------------------------------------

                !Time
                call ExtractDate   (Me%ExtVar%Now , AuxTime(1), AuxTime(2),         &
                                                    AuxTime(3), AuxTime(4),         &  
                                                    AuxTime(5), AuxTime(6))
                TimePointer => AuxTime

                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)            
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR01'


                call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time",                   &
                                     "YYYY/MM/DD HH:MM:SS",                         &
                                     Array1D      = TimePointer,                    &
                                     OutputNumber = Me%OutPut%NextSurfaceOutput,    &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR02'


                !Limits 
                call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR03'


                !Write unsat velocities
                allocate(SurfaceSlice(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                
                !Theta
                do j = JLB, JUB
                do i = ILB, IUB
                    SurfaceSlice(i, j) =  Me%RC%ThetaF(i, j, KUB)
                enddo            
                enddo                
                
                !Please don't change the name of the data sets, since they are used by 
                !MOHID Land Operational
                call HDF5WriteData (Me%ObjHDF5, "/Results_2D/relative water content",               &
                                    "relative water content",                                       &
                                    "m3water/m3water",                                              &
                                    Array2D      =  SurfaceSlice,                                   &
                                    OutputNumber =  Me%OutPut%NextSurfaceOutput,                    &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR01'
                
                deallocate(SurfaceSlice)
                
                !Please don't change the name of the data sets, since they are used by 
                !MOHID Land Operational
                call HDF5WriteData   (Me%ObjHDF5, "/Results_2D/water table depth",  &
                                      "water table depth", "m",                     &
                                      Array2D      = Me%UGWaterDepth2D,             &
                                      OutputNumber = Me%OutPut%NextSurfaceOutput,   &
                                      STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR05'


                Me%OutPut%NextSurfaceOutput = Me%OutPut%NextSurfaceOutput + 1
                

                !----------------------------------------------------------------------------
                !Write everything to disk----------------------------------------------------
                !----------------------------------------------------------------------------

                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR09'
    
            
            endif


        endif

       
    end subroutine PorousMediaOutput

    !--------------------------------------------------------------------------

    subroutine OutPutTimeSeries

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                 :: i, j, STAT_CALL

        !Calculates real Infiltration Velocity
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j) == 1) then
                Me%InfiltrationVelocity (i, j) = - Me%Infiltration(i, j) / Me%ExtVar%DT
            endif
        enddo
        enddo

        !Theta
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = Me%Theta,                                    &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR01'

        !ThetaF
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = Me%RC%ThetaF,                                &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR02'

        !velocity for now......
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = Me%UnsatVelW ,                                     &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR02'

        !infiltration velocity for now
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%InfiltrationVelocity,                           &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR02'
        

        !Head
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = Me%Head,                                     &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR03'

        !Conductivity
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = Me%UnSatK,                                   &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR04'

        !Level Water Table
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%UGWaterLevel2D,                                 &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR07'

        !Depth Water Table
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%UGWaterDepth2D,                                 &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR07'

        !Hydro Pressure
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = Me%HydroPressure,                                  &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR08'

        !Final Head
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = Me%FinalHead,                                      &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR08'

        !Porous Media Water column - to compare to Basin water column
!        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
!                            Data2D = Me%WaterColumn,                                    &
!                            STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR08'
        
        if (Me%ExtVar%ConstructEvaporation) then
            !Evaporation
            call WriteTimeSerie(Me%ObjTimeSerie,                                        &
                                Data2D = Me%EvaporationFlux,                            &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR09'
        endif

        if (Me%ExtVar%ConstructTranspiration) then
            !Transpiration
            call WriteTimeSerie(Me%ObjTimeSerie,                                        &
                                Data3D = Me%ExtVar%TranspirationFlux,                   &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR10'        
        endif


    end subroutine OutPutTimeSeries

    !--------------------------------------------------------------------------

    subroutine ProfileOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Theta
        call WriteProfile(Me%ObjProfile,                                        &
                          Data3D = Me%Theta,                              &
                          SZZ    = Me%ExtVar%SZZ,                               &
                          STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ProfileOutput - ModulePorousMedia - ERR01'

        !ThetaF
        call WriteProfile(Me%ObjProfile,                                        &
                          Data3D = Me%RC%ThetaF,                          &
                          SZZ    = Me%ExtVar%SZZ,                               &
                          STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ProfileOutput - ModulePorousMedia - ERR02'

        !Head
        call WriteProfile(Me%ObjProfile,                                        &
                          Data3D = Me%Head,                               &
                          SZZ    = Me%ExtVar%SZZ,                               &
                          STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ProfileOutput - ModulePorousMedia - ERR02'

    end subroutine ProfileOutput

    !--------------------------------------------------------------------------

    subroutine CalculateTotalStoredVolume(WriteOut)

        !Arguments-------------------------------------------------------------
        logical, optional                           :: WriteOut !Debug only

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        
        Me%TotalStoredVolume = 0.0

        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
            if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then

                Me%TotalStoredVolume = Me%TotalStoredVolume + Me%Theta(i,j,k)             * &
                                       Me%ExtVar%CellVolume(i,j,k)
            endif

        enddo
        enddo
        enddo

        if (present(WriteOut)) then
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%BasinPoints(i, j) == 1) then
                    Me%TotalStoredVolume = Me%TotalStoredVolume + Me%WaterColumn(i, j)* Me%ExtVar%Area(i, j)
                endif
            enddo
            enddo

            write(*,*)Me%TotalStoredVolume
        endif
      
    end subroutine CalculateTotalStoredVolume

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine KillPorousMedia(ObjPorousMediaID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjPorousMediaID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers, STAT_CALL           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mPorousMedia_,  Me%InstanceID)

            if (nUsers == 0) then
                
                !Write Output for continuous computation
                call WriteFinalSoilFile

                !Kills the TimeSerie
                if (Me%ObjTimeSerie /= 0) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillPorousMedia - Porousmedia - ERR05'
                endif

                !Kills Profile Output
                if (Me%OutPut%ProfileON) then
                    call KillProfile(Me%ObjProfile, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillPorousMedia - Porousmedia - ERR01'
                endif
                                    
                !Deassociates External Instances
                if (Me%ObjDrainageNetwork /= 0) then
                    nUsers = DeassociateInstance (mDRAINAGENETWORK_, Me%ObjDrainageNetwork)
                    if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR06'
                endif                

                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR07'

                nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR08'

                nUsers = DeassociateInstance (mGRIDDATA_, Me%ObjTopography)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR09'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR10'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR11'
                
                call KillMap (Me%ObjMap, STAT = STAT_)
                if (STAT_ /= SUCCESS_) stop 'KillPorousMedia - Porousmedia - ERR02'

                call KillGeometry (Me%Objgeometry, STAT = STAT_)
                if (STAT_ /= SUCCESS_) stop 'KillPorousMedia - Porousmedia - ERR03'

                call KillGridData (Me%ObjBottomTopography, STAT = STAT_)
                if (STAT_ /= SUCCESS_) stop 'KillPorousMedia - Porousmedia - ERR01'                
                
                if (Me%OutPut%Yes .or. Me%Output%SurfaceOutput) then                    
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_)
                    if (STAT_ /= SUCCESS_) stop 'KillPorousMedia - Porousmedia - ERR05'
                endif


                !Deallocates Instance
                call DeallocateInstance ()

                ObjPorousMediaID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
                   
    end subroutine KillPorousMedia        

    !--------------------------------------------------------------------------

    subroutine WriteFinalSoilFile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                                        :: Year_File, Month_File, Day_File
        real                                        :: Hour_File, Minute_File, Second_File
        integer                                     :: FinalFile
        integer                                     :: STAT_CALL
        character(LEN = PathLength)                 :: FileName

        !----------------------------------------------------------------------

        if (Me%ExtVar%Now == Me%EndTime) then
            FileName = Me%Files%FinalFile
        else
            FileName = ChangeSuffix(Me%Files%FinalFile,                                 &
                            "_"//trim(TimeToString(Me%ExtVar%Now))//".fin")
        endif            


        call UnitsManager(FinalFile, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalSoilFile - ModulePorousMedia - ERR01'

        open(Unit = FinalFile, File = FileName, Form = 'UNFORMATTED', status = 'UNKNOWN', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalSoilFile - ModulePorousMedia - ERR02'

        !Writes Date
        call ExtractDate(Me%ExtVar%Now, Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File)
        write(FinalFile) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File

        write(FinalFile)Me%Theta

        call UnitsManager(FinalFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalSoilFile - ModulePorousMedia - ERR03'
        

    end subroutine WriteFinalSoilFile

    !------------------------------------------------------------------------    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_PorousMedia), pointer          :: AuxObjPorousMedia
        type (T_PorousMedia), pointer          :: PreviousObjPorousMedia

        !Updates pointers
        if (Me%InstanceID == FirstObjPorousMedia%InstanceID) then
            
            FirstObjPorousMedia => FirstObjPorousMedia%Next
        
        else

            PreviousObjPorousMedia => FirstObjPorousMedia
            AuxObjPorousMedia      => FirstObjPorousMedia%Next
            
            do while (AuxObjPorousMedia%InstanceID /= Me%InstanceID)                
                PreviousObjPorousMedia => AuxObjPorousMedia
                AuxObjPorousMedia      => AuxObjPorousMedia%Next            
            enddo

            !Now update linked list
            PreviousObjPorousMedia%Next => AuxObjPorousMedia%Next

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

    subroutine Ready (ObjPorousMedia_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPorousMedia_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

        if (ObjPorousMedia_ID > 0) then
            call LocateObjPorousMedia (ObjPorousMedia_ID)
            ready_ = VerifyReadLock (mPorousMedia_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjPorousMedia (ObjPorousMediaID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPorousMediaID        

        Me => FirstObjPorousMedia
        do while (associated (Me))
            if (Me%InstanceID == ObjPorousMediaID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModulePorousMedia - LocateObjPorousMedia - ERR01'

    end subroutine LocateObjPorousMedia

    !--------------------------------------------------------------------------

    subroutine ReadLockExternalVar                

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL


        !Time------------------------------------------------------------------
        
        call GetComputeCurrentTime  (Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR01'

        call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR02'
        


        !Topography------------------------------------------------------------

        call GetGridData  (Me%ObjTopography, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR03'


        !BottomTopography------------------------------------------------------
        
        call GetGridData  (Me%ObjBottomTopography, Me%ExtVar%BottomTopoG, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR04'



        !Basin Geometry--------------------------------------------------------                        
        
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR05'  

        call GetRiverPoints   (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR06'


        !Drainage network
        
                      

        !Horizontal Grid-------------------------------------------------------
        
        call GetHorizontalGrid(Me%ObjHorizontalGrid,                                    &
                               DUX = Me%ExtVar%DUX, DVY = Me%ExtVar%DVY,                &
                               DZX = Me%ExtVar%DZX, DZY = Me%ExtVar%DZY,                &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR07'
        

        call GetGridCellArea (Me%ObjHorizontalGrid,                                     & 
                              GridCellArea = Me%ExtVar%Area,                            & 
                              STAT = STAT_CALL )    
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR08'



        !Geometry--------------------------------------------------------------

        call GetGeometryDistances(Me%ObjGeometry,                                       &
                                  SZZ         = Me%ExtVar%SZZ,                          &
                                  DZZ         = Me%ExtVar%DZZ,                          &
                                  DWZ         = Me%ExtVar%DWZ,                          &
                                  ZCellCenter = Me%ExtVar%CenterCell,                   &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR09")


        call GetGeometryVolumes(Me%ObjGeometry,                                         &
                                VolumeZ    = Me%ExtVar%CellVolume,                      &
                                STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR010")


        call GetGeometryKFloor(Me%ObjGeometry,                                          &
                               Z    = Me%ExtVar%KFloor,                                 &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR11")

        call GetGeometryAreas(Me%ObjGeometry,                                           &
                              Me%ExtVar%AreaU, Me%ExtVar%AreaV,                         &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR12")

        call GetGeometryAreas(Me%ObjGeometry,                                           &
                              Me%ExtVar%AreaV, Me%ExtVar%AreaV,                         &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR13")

        


        !Map-------------------------------------------------------------------
        
        call GetWaterPoints3D(Me%ObjMap,                                                &
                              Me%ExtVar%WaterPoints3D,                                  &
                              STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR14")


        call GetOpenPoints3D(Me%ObjMap,                                                 &
                             Me%ExtVar%OpenPoints3D,                                    &
                             STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR15")


        call GetComputeFaces3D(Me%ObjMap,                                               &
                               ComputeFacesU3D = Me%ExtVar%ComputeFacesU3D,             &
                               ComputeFacesV3D = Me%ExtVar%ComputeFacesV3D,             &
                               ComputeFacesW3D = Me%ExtVar%ComputeFacesW3D,             &
                               STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR16")


        !Uderground HorizontalMap----------------------------------------------


    end subroutine ReadLockExternalVar

    !--------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar                

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        
        !Topography------------------------------------------------------------

        call UngetGridData (Me%ObjTopography, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR01'


        !BottomTopography------------------------------------------------------
        
        call UngetGridData (Me%ObjBottomTopography, Me%ExtVar%BottomTopoG, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR02'

        
        !Basin Geometry--------------------------------------------------------                        
          
        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR03'

        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR04'               


        !Horizontal Grid-------------------------------------------------------
        
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DUX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR05'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DVY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR06'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR09'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR10'
        
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%Area, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR11'               
        


        !Geometry--------------------------------------------------------------

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%SZZ,          STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR12")

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%DZZ,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR14") 

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%DWZ,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR15") 

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%CenterCell,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR15a") 


        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%CellVolume,   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR16")

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%KFloor,       STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR18")

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%AreaU,        STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR19")

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%AreaV,       STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR20")
        

        
        !Map-------------------------------------------------------------------
        
        call UnGetMap(Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR21") 
        
        call UnGetMap(Me%ObjMap, Me%ExtVar%OpenPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR22") 


        call UnGetMap(Me%ObjMap, Me%ExtVar%ComputeFacesU3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR23") 


        call UnGetMap(Me%ObjMap, Me%ExtVar%ComputeFacesV3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR24") 


        call UnGetMap(Me%ObjMap, Me%ExtVar%ComputeFacesW3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR25")         


    end subroutine ReadUnLockExternalVar
    
end module ModulePorousMedia

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 

