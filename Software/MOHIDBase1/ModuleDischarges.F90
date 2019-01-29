!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Discharges
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to Impose Discharges
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

Module ModuleDischarges

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions, only : InterpolateValueInTime, ConstructPropertyID
    use ModuleTimeSerie, only : StartTimeSerieInput, StartTimeSerie, GetTimeSerieValue, &
                                WriteTimeSerieLine, KillTimeSerie, GetTimeSerieDataColumns
    use ModuleDrawing

    implicit none 

    private

    !Subroutines & Functions---------------------------------------------------

    !Constructor
    public  :: Construct_Discharges
    private ::      AllocateInstance
    private ::      Construct_DischargeList
    private ::          Add_Discharge
    private ::          Construct_Discharge
    private ::              Construct_Discharge_ID
    private ::              ConstDischargeLoc
    private ::              Read_DataBaseFile
    private ::              Construct_FlowValues
    private ::              Construct_VelocityValues
    private ::              Construct_PropertyList
    private ::                  Construct_Property
    private ::                      Construct_PropertyValues
    private ::                      Add_Property
    private :: ConstructIntakeDischarges
    private :: ConstructLog


    !Selector
    public  :: GetDischargesNumber                                  !Returns the number of discharge points
    public  :: GetDischargesGridLocalization
    public  :: GetDischargesNodeID
    public  :: GetDischargesReservoirID
    public  :: GetIsReservoirOutflow    
    public  :: GetDischargesIDName
    public  :: GetDischargeWaterFlow
    public  :: SetDischargeWaterFlow
    public  :: SetDischargeInterceptionRatio
    public  :: GetDischargeInterceptionRatio    
    public  :: GetDischargeFlowVelocity
    public  :: GetDischargeParameters
    public  :: GetDischargeConcentration
    public  :: GetByPassON
    public  :: GetByPassConcIncrease    
    public  :: GetDischargeFromIntakeON
    public  :: GetIntakePosition
    private ::    Search_Discharge
    private ::    Search_Discharge_ByName
    private ::    Search_Property
    public  :: GetDischargeSpatialEmission
    public  :: GetDischargeSpatialType
    public  :: GetDischargeFlowDistribuiton
    public  :: GetDischargeON
    public  :: IsUpscaling
    public  :: GetDistributionCoefMass    
    public  :: SetLocationCellsZ
    public  :: SetLayer
    public  :: SetDistributionCoefMass
    public  :: UngetDischarges
    !Modifier
    public  :: CorrectsCellsDischarges
    public  :: CorrectsByPassCellsDischarges
    public  :: TryIgnoreDischarge
    !Destructor
    public  :: Kill_Discharges
    private ::      KillIndividualDischarge
    private ::      DeallocateInstance


    !Management
    private ::      Ready

    private :: UngetDischarges1Dreal
    private :: UngetDischarges1Dinteger
    interface  UngetDischarges
        module procedure UngetDischarges1Dinteger
        module procedure UngetDischarges1Dreal
    end interface UngetDischarges


    !Parameter-----------------------------------------------------------------

    !STAT
    integer, parameter :: NO_ID_          =  7    ! No property ID was specified

    ! Time Series 
    integer, parameter :: TIME_COLUMN     = 1
    integer, parameter :: SECONDS         = 1
    integer, parameter :: MINUTES         = 2
    integer, parameter :: HOURS           = 3
    integer, parameter :: DAYS            = 4
    integer, parameter :: MONTHS          = 5
    integer, parameter :: YEARS           = 6

    !Direction
    integer, parameter :: DirectionX_   = 1
    integer, parameter :: DirectionY_   = 2

    !Discharges type 
    integer, parameter :: Normal        = 1
    integer, parameter :: FlowOver      = 2
    integer, parameter :: Valve         = 3
    integer, parameter :: OpenMILink    = 4
    integer, parameter :: RatingCurve   = 5

    !Valve side
    integer, parameter :: SideA         = 1
    integer, parameter :: SideB         = 2
    
    !Valve - section type
    integer, parameter :: circular_area     = 1
    integer, parameter :: rectangular_area  = 2
    


    character(len=StringLength), parameter :: block_begin               = '<begindischarge>'
    character(len=StringLength), parameter :: block_end                 = '<enddischarge>'
    character(len=StringLength), parameter :: beginproperty             = '<<beginproperty>>'
    character(len=StringLength), parameter :: endproperty               = '<<endproperty>>'

    !Types---------------------------------------------------------------------

    type       T_ID
        integer                                 :: IDnumber    = FillValueInt
        character(LEN = StringLength)           :: name        = Space
        character(LEN = StringLength)           :: description = Space
        character(LEN = StringLength)           :: units       = Space
    end type T_ID

    type       T_Property
        type (T_PropertyID)                     :: ID
        logical                                 :: Variable         = .false.
        integer                                 :: ConcColumn       = null_int 
        real                                    :: scalar           = FillValueReal
        logical                                 :: TimeSerieON      = .false. 
        integer                                 :: TimeSerie        = 0
        integer                                 :: TimeSerieMaxCol  = 0        
        logical                                 :: PropTimeSerie    = .false.
        logical                                 :: FromIntake       = .false.
        real                                    :: IncreaseValue    = FillValueReal 
        type (T_Property), pointer              :: Next => null(), &
                                                   Prev => null()
    end type T_Property

    type       T_WaterFlow
        logical                                 :: Variable      = .false.
        integer                                 :: FlowColumn    = null_int 
        real                                    :: scalar        = FillValueReal
        logical                                 :: Upscaling     = .false.
    end type T_WaterFlow


    type       T_WaterVelocity
        logical                                 :: UVariable    = .false.  
        logical                                 :: VVariable    = .false.  
        logical                                 :: WVariable    = .false.  
        integer                                 :: UColumn      = null_int 
        integer                                 :: VColumn      = null_int 
        integer                                 :: WColumn      = null_int 
        real                                    :: Uscalar = FillValueReal
        real                                    :: Vscalar = FillValueReal
        real                                    :: Wscalar = FillValueReal        
    end type T_WaterVelocity

    type       T_FlowOver
        real                                    :: WeirLength           = null_real 
        real                                    :: DischargeCoeficient  = null_real 
        real                                    :: CrestHeigth          = null_real 
    end  type T_FlowOver

    type       T_RatingCurve
        real, dimension(:), pointer             :: Level
        real, dimension(:), pointer             :: Flow
        integer                                 :: nValues
    end  type T_RatingCurve
    
    type       T_Valve
        integer                                 :: SectionType          = null_int
        real                                    :: Diameter             = null_real 
        real                                    :: Height               = null_real 
        real                                    :: Width                = null_real         
        real                                    :: DischargeCoeficient  = null_real 
        real                                    :: AxisHeigth           = null_real 
        real                                    :: SillHeigth           = null_real
        real                                    :: PipeLength           = null_real 
        real                                    :: PipeManning          = null_real 
        real                                    :: AreaInTime           = null_real 
    end  type T_Valve

    type T_GridCoordinates
        integer                                 :: I    = FillValueInt
        integer                                 :: J    = FillValueInt
        integer                                 :: K    = FillValueInt
        integer                                 :: OldI = FillValueInt
        integer                                 :: OldJ = FillValueInt
    end type T_GridCoordinates

    type       T_Localization 
         type (T_GridCoordinates)               :: GridCoordinates
         logical                                :: AlternativeLocations             = .false.
         real                                   :: MinimumDischargeDepth            = null_real 
         logical                                :: StartFromLastDischargeLocation   = .false.
         logical                                :: Location2D                       = .false. 
         integer                                :: DischVertical                    = FillValueInt
         real                                   :: Kdepth                           = FillValueReal
         integer                                :: NodeID                           = FillValueInt
         integer                                :: ReservoirID                      = FillValueInt
         logical                                :: TrackLocation                    = .false. 
         character(len=StringLength)            :: TrackLocationFile                = null_str 
         integer                                :: TrackLocationFileUnitNumber      = null_int 
         logical                                :: UseDischargePathFile             = .false.
         character(len=StringLength)            :: DischargePathFile                = null_str 
         real                                   :: CoordinateX                      = null_real, & 
                                                   CoordinateY                      = null_real    
         integer                                :: XColumn                          = null_int, & 
                                                   YColumn                          = null_int    
         logical                                :: VariableX                        = .false., & 
                                                   VariableY                        = .false.    
         logical                                :: CoordinatesON                    = .false.
         integer                                :: HorizontalType                   = FillValueInt
         logical                                :: CellCorrect                      = .false.
         integer                                :: SpatialEmission                  = null_int    
         character(len=StringLength)            :: SpatialFile                      = null_str 
         type (T_Polygon), pointer              :: Polygon                          => null()
         type (T_Lines),   pointer              :: Line                             => null()
         type (T_XYZPoints), pointer            :: XYZPoints                        => null()         
         integer                                :: nCells                           = 1
         integer, dimension(:), pointer         :: VectorI                          => null(), &
                                                   VectorJ                          => null(), &
                                                   VectorK                          => null()
         integer                                :: FlowDistribution                 = DischByCell_
         integer                                :: kmin                             = FillValueInt
         integer                                :: kmax                             = FillValueInt

         !Important for the domain decomposition approach
         !is the ratio of the XYZPoints or Line or Polygon that intercepts the model domain
         real                                   :: InterceptionRatio               = FillValueReal
         
         real,    dimension(:    ), pointer     :: DistributionCoefMass             => null()
    end  type T_Localization

    type      T_FromIntake
        character(len=PathLength)               :: IntakeName                       = null_str
        integer                                 :: IntakeID                         = FillValueInt
        logical                                 :: ON                               = .false.
        logical                                 :: AssociateFlow                    = .false.
        real                                    :: FlowFraction                     = 1.0
    end type  T_FromIntake

    type      T_ByPass 
         integer                                :: i        = null_int
         integer                                :: j        = null_int    
         integer                                :: k        = null_int
         real                                   :: X        = null_real
         real                                   :: Y        = null_real
         logical                                :: ON       = .false.
         logical                                :: OneWay   = .false.    
         integer                                :: Side     = null_int 
    end  type T_ByPass


    type       T_IndividualDischarge
         type(T_ID                 )            :: ID
         type(T_Localization       )            :: Localization
         integer                                :: PropertiesNumber = FillValueInt
         character(len=PathLength)              :: DataBaseFile     = null_str  
         character(len=PathLength)              :: OutPutFile       = null_str    
         logical                                :: TimeSerieON      = .false.  
         integer                                :: TimeSerie        = 0
         integer                                :: TimeSerieMaxCol  = 0       
         logical                                :: UseOriginalValues    = .false.  
         type(T_WaterFlow          )            :: WaterFlow   
         type(T_WaterVelocity      )            :: VelocityFlow
         integer                                :: DischargeType    = null_int  
         type(T_Valve   )                       :: Valve
         type(T_FlowOver)                       :: FlowOver
         type(T_RatingCurve)                    :: RatingCurve
         type(T_Property           ), pointer   :: FirstProperty    => null()
         type(T_Property           ), pointer   :: LastProperty     => null()
         type(T_Property           ), pointer   :: CurrProperty     => null()
         type(T_IndividualDischarge), pointer   :: Next             => null()
         type(T_IndividualDischarge), pointer   :: Prev             => null()
         type(T_ByPass             )            :: ByPass
         type(T_FromIntake         )            :: FromIntake       
         logical                                :: IgnoreON         = .false.  
         logical                                :: IsReservoirOutflow = .false.
    end type T_IndividualDischarge    

    type      T_Discharges
         integer                                :: InstanceID       = null_int             
         integer                                :: ObjEnterData     = 0
         integer                                :: ObjTime          = 0
         character(len=Pathlength)              :: DataFile         = null_str 
         integer                                :: DischargesNumber = FillValueInt
         type (T_IndividualDischarge), pointer  :: FirstDischarge   => null()
         type (T_IndividualDischarge), pointer  :: LastDischarge    => null()
         type (T_IndividualDischarge), pointer  :: CurrentDischarge => null()
         type (T_Discharges), pointer           :: Next             => null()
         logical                                :: IgnoreON         = .false.
         integer                                :: ReferentialZ     = FillValueInt          
    end type T_Discharges

    !Global Variables
    type (T_Discharges), pointer                :: FirstDischarges  => null()
    type (T_Discharges), pointer                :: Me               => null()
    
    !--------------------------------------------------------------------------

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Construct_Discharges(DischargesID, ObjTime, DataFile, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: DischargesID
        integer                                     :: ObjTime
        character(len=*), optional                  :: DataFile
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_

        integer                                     :: STAT_CALL, flag
                                 
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mDischarges_)) then
            nullify (FirstDischarges)
            call RegisterModule (mDischarges_) 
        endif

        call Ready(DischargesID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            !Allocates Instance
            call AllocateInstance

            !Associates Time
            Me%ObjTime = AssociateInstance   (mTIME_, ObjTime)


            if (present(DataFile)) then
                Me%DataFile = DataFile
            else
                call ReadFileName('DISCHARG', Me%DataFile, Message = "Discharges Data File", STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    write(*,*    ) 'DISCHARG keyword not found in Nomfich'
                    stop 'Construct_Discharges - ModuleDischarges - ERR10'
                endif
            endif
        
            ! Construct one instance to use the moduleEnterData
            call ConstructEnterData(Me%ObjEnterData, Me%DataFile, STAT = STAT_CALL)

cd1 :       if      ( STAT_CALL .EQ. FILE_NOT_FOUND_ERR_) then

                write(*,*    ) 
                write(*,*    ) 'Fatal error ! Discharges data file not found' 
                write(*,'(A)') 'File : ', trim(adjustl(Me%DataFile))
                write(*,*    ) 'look at DISCHARGES KeyWord at nomfich.dat file  '
                stop 'Construct_Discharges - ModuleDischarges - ERR20'  

            else if ((STAT_CALL .NE. FILE_NOT_FOUND_ERR_) .AND.                &
                     (STAT_CALL .NE. SUCCESS_          )) then cd1
                stop 'Subroutine Construct_Discharges - ModuleDischarges. ERR30.' 
            end if cd1

            call GetData(Me%IgnoreON,                                           &
                         Me%ObjEnterData,                                       &
                         flag,                                                  &
                         SearchType   = FromFile,                               &
                         keyword      ='IGNORE_ON',                             &
                         Default      = .false.,                                &
                         ClientModule ='ModuleDischarges',                      &
                         STAT         = STAT_CALL)        

            if (STAT_CALL  /= SUCCESS_) stop 'Construct_Discharges - ModuleDischarges - ERR40'

            call GetData(Me%ReferentialZ,                                       &
                         Me%ObjEnterData,                                       &
                         flag,                                                  &
                         SearchType   = FromFile,                               &
                         keyword      ='REFERENTIAL_Z',                         &
                         Default      = Hydrographic_,                          &
                         ClientModule ='ModuleDischarges',                      &
                         STAT         = STAT_CALL)        

            if (STAT_CALL  /= SUCCESS_) stop 'Construct_Discharges - ModuleDischarges - ERR42'
            
            if (Me%ReferentialZ /= Hydrographic_ .and. Me%ReferentialZ /=  Topographic_) then
                stop 'Construct_Discharges - ModuleDischarges - ERR44'
            endif


            ! Constructs the discharge list 
            call Construct_DischargeList

            call ConstructIntakeDischarges

            !User Feed-Back
            call ConstructLog

            call KillEnterData  (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'Construct_Discharges - ModuleDischarges - ERR50'

            !Returns ID
            DischargesID    = Me%InstanceID

            STAT_ = SUCCESS_
        else 

            stop 'ModuleDischarges - Construct_Discharges - ERR60' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end Subroutine Construct_Discharges

    !--------------------------------------------------------------------------


    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Discharges), pointer                :: NewDischarges
        type (T_Discharges), pointer                :: PreviousDischarges


        !Allocates new instance
        allocate (NewDischarges)
        nullify  (NewDischarges%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstDischarges)) then
            FirstDischarges         => NewDischarges
            Me                      => NewDischarges
        else
            PreviousDischarges      => FirstDischarges
            Me                      => FirstDischarges%Next
            do while (associated(Me))
                PreviousDischarges  => Me
                Me                  => Me%Next
            enddo
            Me                      => NewDischarges
            PreviousDischarges%Next => NewDischarges
        endif

        Me%InstanceID = RegisterNewInstance (mDISCHARGES_)

        ! Initialize the  discharges number   
        Me%DischargesNumber = 0

        ! Initialize the water discharges list   
        nullify (Me%FirstDischarge)
        nullify (Me%LastDischarge)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine Construct_DischargeList

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
                                  
        integer :: STAT_CALL
        integer :: IDnumber 
        integer :: ClientNumber

        logical :: BlockFound
                                  
        !Local-----------------------------------------------------------------

        type (T_IndividualDischarge), pointer :: NewDischarge

        !----------------------------------------------------------------------

        IDnumber = 0

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber, &
                                        block_begin, block_end, BlockFound,       &
                                        STAT = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_      ) then    
cd2 :           if (BlockFound) then                                  
                
                    ! Increments the number of Discharges
                    IDnumber = IDnumber + 1
                    ! Construct a New Property 
                    Call Construct_Discharge(NewDischarge, IDnumber, ClientNumber)
                                             

                    ! Add new Property to the discharge List 
                    Call Add_Discharge(NewDischarge)

                else
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                               &
                        stop 'Subroutine Construct_DischargeList - ModuleDischarges. ERR01.'

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                stop 'Subroutine Construct_DischargeList - ModuleDischarges. ERR02.'
            end if cd1
        end do do1

        !----------------------------------------------------------------------

    end subroutine Construct_DischargeList

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct a new property.           

    subroutine Construct_Discharge(NewDischarge, IDnumber, ClientNumber)

        !Arguments--------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: NewDischarge
        integer, intent(IN)                         :: ClientNumber
        integer, intent(IN)                         :: IDnumber 
         
        !----------------------------------------------------------------------


        allocate (NewDischarge)
        nullify  (NewDischarge%Next)
        nullify  (NewDischarge%Prev)


        nullify (NewDischarge%FirstProperty )
        nullify (NewDischarge%LastProperty  )

        !Construct Discharge ID
        call Construct_Discharge_ID           (NewDischarge, IDNumber)

        !Construct Discharge Localization
        call ConstDischargeLoc                (NewDischarge)

        !Read_DataBaseFile
        call Read_DataBaseFile                (NewDischarge)

        !Construct Discharge Flow values
        call Construct_FlowValues             (NewDischarge, ClientNumber)

        !Construct Discharge Velocity values
        call Construct_VelocityValues         (NewDischarge)

        !Construct Discharge Variable Location
        call Construct_VariableLocation       (NewDischarge)

        !Construct Property List
        call Construct_PropertyList           (NewDischarge, ClientNumber)

        !----------------------------------------------------------------------

    end subroutine Construct_Discharge

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------
    ! This subroutine adds a new discharge to the discharge List  

    subroutine Add_Discharge(NewDischarge)

        !Arguments--------------------------------------------------------------
                
        type(T_IndividualDischarge), pointer :: NewDischarge

        !----------------------------------------------------------------------

        ! Add to the discharge list a new discharge
        if (.not.associated(Me%FirstDischarge)) then
            Me%DischargesNumber   = 1
            Me%FirstDischarge     => NewDischarge
            Me%LastDischarge      => NewDischarge
        else
            NewDischarge%Prev     => Me%LastDischarge
            Me%LastDischarge%Next => NewDischarge
            Me%LastDischarge      => NewDischarge
            Me%DischargesNumber   = Me%DischargesNumber + 1
        end if 

        !----------------------------------------------------------------------

    end subroutine Add_Discharge 

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct the property ID          

    subroutine Construct_Discharge_ID(NewDischarge,IDNumber)
        
        !Arguments--------------------------------------------------------------
        type(T_IndividualDischarge), pointer    :: NewDischarge
        integer,                     intent(IN) :: IDnumber

        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------

        integer :: flag

        !----------------------------------------------------------------------
     


        ! The property ID number is incremented when a new discharge is created
        NewDischarge%ID%IDnumber = IDnumber 


        ! Property name and is units don't have a default value. The program stops when a
        ! is not sepecified the property name and untis 
        call GetData(NewDischarge%ID%name,                                              &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     FromBlock,                                                         &
                     keyword='NAME',                                                    &
                     ClientModule = 'ModuleDischarges')

        if (flag==0) stop 'Discharges without name - Construct_Discharge'


        ! The property description is a character*132 where the user can
        ! store information about the property, example: 
        ! spring climatologic temperature field for the North Atlantic
        call GetData(NewDischarge%ID%Description,                                       &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     FromBlock,                                                         &
                     keyword='DESCRIPTION',                                             &
                     default='No description was given.',                               &
                     ClientModule = 'ModuleDischarges')


        !----------------------------------------------------------------------

     end subroutine Construct_Discharge_ID

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    Subroutine ConstDischargeLoc (NewDischarge)
        
        !Arguments--------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: NewDischarge

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                     :: flag
        logical                                     :: AuxLog
        character(LEN = StringLength)               :: RootPath, DischargeName, AuxName
        
        !----------------------------------------------------------------------
        !X Location
        call GetData(NewDischarge%Localization%CoordinateX,                         &
                     Me%ObjEnterData,                                               &
                     flag,                                                          &
                     FromBlock,                                                     &
                     keyword      ='COORD_X',                                       &
                     ClientModule = 'ModuleDischarges',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstDischargeLoc - ModuleDischarges - ERR110.'

        !When X_COORDINATE is found, the model assumes that the discharge location is given using coordinates
        if (flag/=0) then
            NewDischarge%Localization%CoordinatesON = .true.
            NewDischarge%Localization%Location2D    = .true.
        else
            NewDischarge%Localization%CoordinatesON = .false.
            NewDischarge%Localization%Location2D    = .false.
        endif

i3:     if (NewDischarge%Localization%CoordinatesON) then

            !Y Location
            call GetData(NewDischarge%Localization%CoordinateY,                     &
                         Me%ObjEnterData,                                           &
                         flag,                                                      &
                         FromBlock,                                                 &
                         keyword      ='COORD_Y',                                   &
                         ClientModule = 'ModuleDischarges',                         &
                         STAT         = STAT_CALL)

            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR120.'
            if (flag==0)                                                            &
                stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR130.'

        endif i3
        
i4:     if (.not. NewDischarge%Localization%CoordinatesON) then
            !I Location
            call GetData(NewDischarge%Localization%GridCoordinates%I,                       &
                         Me%ObjEnterData,                                                   &
                         flag,                                                              &
                         FromBlock,                                                         &
                         keyword      ='I_CELL',                                            &
                         ClientModule = 'ModuleDischarges',                                 &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstDischargeLoc - ModuleDischarges - ERR10.'

            !When I_CELL is found, the model assumes that the discharge is grid based
            if (flag/=0) then
                NewDischarge%Localization%Location2D = .true.
            else
                NewDischarge%Localization%Location2D = .false.
            endif

i8:         if (NewDischarge%Localization%Location2D) then

                !J Location
                call GetData(NewDischarge%Localization%GridCoordinates%J,                   &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             FromBlock,                                                     &
                             keyword      ='J_CELL',                                        &
                             ClientModule = 'ModuleDischarges',                             &
                             STAT         = STAT_CALL)

                if (STAT_CALL .NE. SUCCESS_)                                                &
                    stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR20.'
                if (flag==0)                                                                &
                    stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR30.'
        
            endif i8

        endif i4

i1:     if (NewDischarge%Localization%Location2D) then


            call GetData(NewDischarge%Localization%DischVertical,                       &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='VERTICAL_DISCHARGE',                            &
                         ClientModule = 'ModuleDischarges',                             &
                         default      = DischLayer_,                                    &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR40.'

            call GetData(AuxLog,                                                        &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='DISCHARGE_UNIFORM',                             &
                         ClientModule = 'ModuleDischarges',                             &
                         default      = .false.,                                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR50.'

            if (AuxLog) NewDischarge%Localization%DischVertical = DischUniform_

            !Alway read a default layer 
            call GetData(NewDischarge%Localization%GridCoordinates%K,                   &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='K_CELL',                                        &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR60.'
            
            if (NewDischarge%Localization%DischVertical .eq. DischLayer_ .and. flag .eq. 0)then
                write(*,*)"You must define the K_CELL in the discharge: "//trim(NewDischarge%ID%Name)
                stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR61.'
            endif

            select case (NewDischarge%Localization%DischVertical)

                case (DischLayer_)
                    !do not do nothing                 
                case (DischDepth_)

                    call GetData(NewDischarge%Localization%Kdepth,                      &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                  &
                                 FromBlock,                                             &
                                 keyword      ='K_DEPTH',                               &
                                 ClientModule = 'ModuleDischarges',                     &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                        stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR70.'

                case (DischUniform_)


                    call GetData(NewDischarge%Localization%kmin,                        &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                  &
                                 FromBlock,                                             &
                                 keyword      ='K_MIN',                                 &
                                 default      = FillValueInt,                           &
                                 ClientModule = 'ModuleDischarges',                     &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                        stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR74.'


                    call GetData(NewDischarge%Localization%kmax,                        &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                  &
                                 FromBlock,                                             &
                                 keyword      ='K_MAX',                                 &
                                 default      = FillValueInt,                           &
                                 ClientModule = 'ModuleDischarges',                     &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                        stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR76.'

                case (DischBottom_, DischSurf_)
                    !do not do nothing 
                case default
                    write(*,*) "VERTICAL DISCHARGE option not known ", NewDischarge%Localization%DischVertical

                    write(*,*) "The known options are : "," Bottom=",DischBottom_," Surface=",DischSurf_,&
                                                          " Layer =",DischLayer_, " Depth  =",DischDepth_,&
                                                          " Uniform=",DischUniform_
                    stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR80'

            end select




            !Searches for alternative locations. These are place where the discharge should
            !happen if the disharge location at a given time is not an open point.
            call GetData(NewDischarge%Localization%AlternativeLocations,                &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='ALTERNATIVE_LOCATIONS',                         &
                         ClientModule = 'ModuleDischarges',                             &
                         default      = .false.,                                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR90.'

i2:         if (NewDischarge%Localization%AlternativeLocations) then
                call GetData(NewDischarge%Localization%MinimumDischargeDepth,           &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             FromBlock,                                                 &
                             keyword      ='MINIMUM_DISCHARGE_DEPTH',                   &
                             ClientModule = 'ModuleDischarges',                         &
                             Default      = 1.0,                                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR100'


                 call GetData(NewDischarge%Localization%StartFromLastDischargeLocation, &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             FromBlock,                                                 &
                             keyword      ='START_FROM_LAST_LOCATION',                  &
                             ClientModule = 'ModuleDischarges',                         &
                             Default      = .false.,                                    &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR110'
                                    
                
                call GetData(NewDischarge%Localization%TrackLocation,                   &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             FromBlock,                                                 &
                             keyword      ='TRACKLOCATION',                             &
                             ClientModule ='ModuleDischarges',                          &
                             Default      = .false.,                                    &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR120'
                
                if (NewDischarge%Localization%TrackLocation) then
                    
                    !Gets the root path from the file nomfich.dat
                    call ReadFileName("ROOT_SRT", RootPath, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) then
                        call ReadFileName("ROOT", RootPath, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) then
                            call ReadFileName("RAIZ", RootPath, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                              &
                                stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR130'
                        endif
                    endif

                    !Construct name for the file                    
                    DischargeName = "location"//"_"//trim(NewDischarge%ID%name)

                    NewDischarge%Localization%TrackLocationFile =   trim(adjustl(RootPath      ))// &
                                                                    trim(adjustl(DischargeName ))// &
                                                                    ".txt"  
                    
                    !Opens the file
                    call UnitsManager(NewDischarge%Localization%TrackLocationFileUnitNumber,        &
                                        OPEN_FILE, STAT = STAT_CALL)
                    
                    open(UNIT   = NewDischarge%Localization%TrackLocationFileUnitNumber,            &
                         FILE   = NewDischarge%Localization%TrackLocationFile, STATUS  = "UNKNOWN", &
                         IOSTAT  = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) then
                        write(*,*) 'Error opening diacharge location track file ',                  &
                                trim(   NewDischarge%Localization%TrackLocationFile)
                        stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR140'
                    endif

                endif

            endif i2   

            call GetData(AuxName,                                                       &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     FromBlock,                                                         &
                     keyword      ='SPATIAL_EMISSION',                                  &
                     ClientModule = 'ModuleDischarges',                                 &
                     Default      = "Point",                                            &
                     STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR150'

            select case (AuxName)

                case ("Point")
                    NewDischarge%Localization%SpatialEmission = DischPoint_
                case ("Line")
                    NewDischarge%Localization%SpatialEmission = DischLine_
                case ("Polygon")
                    NewDischarge%Localization%SpatialEmission = DischPolygon_
                case ("XYZPoints")
                    NewDischarge%Localization%SpatialEmission = DischXYZPoints_                    
                case default
                    write(*,*) "SPATIAL EMISSION option not known ",trim(AuxName)," ????"
                    write(*,*) "The known options are : ","Point ", "Line ", "Polygon"
                    stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR160'

            end select

            if (NewDischarge%Localization%SpatialEmission == DischLine_ .or.            &
                NewDischarge%Localization%SpatialEmission == DischPolygon_ .or.         &
                NewDischarge%Localization%SpatialEmission == DischXYZPoints_) then
                call GetData(NewDischarge%Localization%SpatialFile,                     &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='SPATIAL_FILE',                                  &
                         ClientModule = 'ModuleDischarges',                             &
                         Default      = "Point",                                        &
                         STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR170'

                if (flag == 0)                                                          &
                    stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR180'

                if (NewDischarge%Localization%SpatialEmission == DischLine_)            &
                    call New(NewDischarge%Localization%Line, NewDischarge%Localization%SpatialFile)

                if (NewDischarge%Localization%SpatialEmission == DischPolygon_)         &
                    call New(NewDischarge%Localization%Polygon, NewDischarge%Localization%SpatialFile)

                if (NewDischarge%Localization%SpatialEmission == DischXYZPoints_)         &
                    call New(NewDischarge%Localization%XYZPoints, NewDischarge%Localization%SpatialFile)

                call GetData(AuxName,                                                   &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='FLOW_DISTRIBUTION',                             &
                         ClientModule = 'ModuleDischarges',                             &
                         Default      = "by cell",                                      &
                         STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR190'

                select case (AuxName)

                    case ("by cell")
                        NewDischarge%Localization%FlowDistribution = DischByCell_
                    case ("by water column")
                        NewDischarge%Localization%FlowDistribution = DischByWaterColumn_
                    case ("by volume")
                        NewDischarge%Localization%FlowDistribution = DischByVolume_
                    case default
                        write(*,*) "FLOW_DISTRIBUTION option not known ",trim(AuxName)," ????"
                        write(*,*) "The known options are : ","by cell ", "by water column ", "by volume"
                        stop 'Subroutine ConstDischargeLoc - ModuleDischarges. ERR200'

                end select

            endif

        else    i1

            !Darinage Network Dishcarge
            call GetData(NewDischarge%Localization%NodeID,                              &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='NODE_ID',                                       &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstDischargeLoc - ModuleDischarges - ERR210'
            
            !MOHID Land Reservoirs Discharge
            if (flag == 0) then
               call GetData(NewDischarge%Localization%ReservoirID,                          &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             FromBlock,                                                     &
                             keyword      ='RESERVOIR_ID',                                  &
                             ClientModule = 'ModuleDischarges',                             &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstDischargeLoc - ModuleDischarges - ERR220'        
                
                !is an imposed outflow? to distinguish from other discharge types (abstraction and input)
               call GetData(NewDischarge%IsReservoirOutflow,                                &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             FromBlock,                                                     &
                             keyword      ='IS_OUTFLOW',                                    &
                             ClientModule = 'ModuleDischarges',                             &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstDischargeLoc - ModuleDischarges - ERR230'                      
            endif
            
        endif   i1


        !----------------------------------------------------------------------

     End Subroutine ConstDischargeLoc

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine Read_DataBaseFile(NewDischarge)

        !Arguments-------------------------------------------------------------
        type(T_IndividualDischarge), pointer                :: NewDischarge

        !External--------------------------------------------------------------
        integer                                             :: flag, STAT_CALL
        !----------------------------------------------------------------------


        !Looks for a definition of the data base file. If there is one, this module assumes 
        !that the discharge is time variable
        call GetData(NewDischarge%DataBaseFile,                                         &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     FromBlock,                                                         &
                     keyword      = 'DATA_BASE_FILE',                                   &
                     ClientModule = 'ModuleDischarges',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Read_DataBaseFile - ModuleDischarges - ERR10'

        if (flag == 1) then
            NewDischarge%TimeSerieON = .true.
        else 
            NewDischarge%TimeSerieON = .false.
        endif


        !Start TimeSerie Input
        if (NewDischarge%TimeSerieON) then
            call StartTimeSerieInput(NewDischarge%TimeSerie, NewDischarge%DataBaseFile, &
                                     Me%ObjTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Read_DataBaseFile - ModuleDischarges - ERR20'
            
            call GetData(NewDischarge%UseOriginalValues,                                &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromFile,                                                      &
                         keyword      = 'USE_ORIGINAL_VALUES',                          &
                         default      = .false.,                                        &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Read_DataBaseFile - ModuleDischarges - ERR30'
            
            call GetTimeSerieDataColumns (TimeSerieID = NewDischarge%TimeSerie,         &
                                          DataColumns = NewDischarge%TimeSerieMaxCol,   &
                                          STAT        = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Read_DataBaseFile - ModuleDischarges - ERR40'

        end if 
        

        !----------------------------------------------------------------------

    end subroutine Read_DataBaseFile 

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    
    subroutine Construct_FlowValues(NewDischarge, ClientNumber)

        !Arguments--------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: NewDischarge
        integer                                     :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                                     :: flag, STAT_CALL
        character(len = StringLength), parameter    :: beginratingcurve = '<<begin_rating_curve>>'
        character(len = StringLength), parameter    :: endratingcurve   = '<<end_rating_curve>>'
        real, dimension(:), pointer                 :: BufferLine
        logical                                     :: BlockLayersFound
        integer                                     :: FirstLine, LastLine
        integer                                     :: iValue, iLine
        
        
        !----------------------------------------------------------------------


        !Searches for the default flow value, by default the value is zero
        call GetData(NewDischarge%WaterFlow%scalar,                                     &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     FromBlock,                                                         &
                     keyword      ='DEFAULT_FLOW_VALUE',                                &
                     ClientModule = 'ModuleDischarges',                                 &
                     default      = 0.0,                                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR10'



        !Searches for the flow column (if the discharge is read from a data base)
i1:     if (NewDischarge%TimeSerieON) then
            call GetData(NewDischarge%WaterFlow%FlowColumn,                             &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='FLOW_COLUMN',                                   &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR20'

            if (flag == 1) then

                NewDischarge%WaterFlow%Variable = .true.
                
                if (NewDischarge%WaterFlow%FlowColumn > NewDischarge%TimeSerieMaxCol .or.   &
                    NewDischarge%WaterFlow%FlowColumn < 2) then
                    
                    write(*,*)'Look at file', trim(Me%DataFile)
                    write(*,*)'FLOW_COLUMN not valid in discharge ', trim(NewDischarge%ID%Name)
                    stop 'Construct_FlowValues - ModuleDischarges - ERR20'
                    
                endif
                
            else
                write(*,*)
                write(*,*)'Look at file                        : ', trim(Me%DataFile)
                write(*,*)'FLOW_COLUMN not found in discharge  : ', trim(NewDischarge%ID%Name)
                write(*,*)'Discharge flow will be assumed constant'
                write(*,*)
                NewDischarge%WaterFlow%Variable = .false.            
            endif
            
        endif i1


        !Check if the discharge is of the type flow over
        call GetData(NewDischarge%DischargeType,                                        &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     FromBlock,                                                         &
                     keyword      ='DISCHARGE_TYPE',                                    &
                     ClientModule = 'ModuleDischarges',                                 &
                     default      = Normal,                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR30'

        if      (NewDischarge%DischargeType /= Normal       .and.                       &
                 NewDischarge%DischargeType /= FlowOver     .and.                       &
                 NewDischarge%DischargeType /= Valve        .and.                       &
                 NewDischarge%DischargeType /= OpenMILink   .and.                       & 
                 NewDischarge%DischargeType /= RatingCurve)   then
                 stop 'Construct_FlowValues - ModuleDischarges - ERR40'
        endif
        
i2:     if (NewDischarge%DischargeType == FlowOver) then

            call GetData(NewDischarge%FlowOver%WeirLength,                              &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='WEIR_LENGTH',                                   &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR50'
            if (flag /= 1) then
                write(*,*) 'Weir Length Missing'
                stop ' Construct_FlowValues - ModuleDischarges - ERR60'
            endif

            call GetData(NewDischarge%FlowOver%DischargeCoeficient,                     &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='WEIR_COEF',                                     &
                         ClientModule = 'ModuleDischarges',                             &
                         default      = 0.4,                                            &
                        STAT         = STAT_CALL)
            
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR70'

            call GetData(NewDischarge%FlowOver%CrestHeigth,                             &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='CREST_HEIGTH',                                  &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR80'
            
            if (flag /= 1) then
                write(*,*) 'Crest Height Missing'
                stop ' Construct_FlowValues - ModuleDischarges - ERR90'
            endif
            
            if (Me%ReferentialZ == Hydrographic_) then
                NewDischarge%FlowOver%CrestHeigth =  - NewDischarge%FlowOver%CrestHeigth
            endif                
            
        else if (NewDischarge%DischargeType == Valve) then i2
        

            call GetData(NewDischarge%Valve%SectionType,                                &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='VALVE_SECTION_TYPE',                            &
                         default      = circular_area,                                      &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR100'
            
            if      (NewDischarge%Valve%SectionType == circular_area) then

                call GetData(NewDischarge%Valve%Diameter,                                   &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             FromBlock,                                                     &
                             keyword      ='VALVE_DIAMETER',                                &
                             ClientModule = 'ModuleDischarges',                             &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR102'
                
                if (flag /= 1) then
                    write(*,*) 'Valve diameter missing - circular section'
                    stop ' Construct_FlowValues - ModuleDischarges - ERR104'
                endif
                
            elseif  (NewDischarge%Valve%SectionType == rectangular_area) then
            

                call GetData(NewDischarge%Valve%Height,                                     &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             FromBlock,                                                     &
                             keyword      ='VALVE_HEIGHT',                                  &
                             ClientModule = 'ModuleDischarges',                             &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR106'
                
                if (flag /= 1) then
                    write(*,*) 'Valve height missing - rectangular section'
                    stop ' Construct_FlowValues - ModuleDischarges - ERR107'
                endif
                
                call GetData(NewDischarge%Valve%Width,                                  &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             FromBlock,                                                 &
                             keyword      ='VALVE_WIDTH',                               &
                             ClientModule = 'ModuleDischarges',                         &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR108'
                
                if (flag /= 1) then
                    write(*,*) 'Valve width missing - rectangular section'
                    stop ' Construct_FlowValues - ModuleDischarges - ERR109'
                endif
                                
            else
            
                stop ' Construct_FlowValues - ModuleDischarges - ERR110'
            
            endif

            call GetData(NewDischarge%Valve%DischargeCoeficient,                        &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='VALVE_COEF',                                    &
                         ClientModule = 'ModuleDischarges',                             &
                         default      = 1.,                                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR120'

            call GetData(NewDischarge%Valve%SillHeigth,                                 &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='VALVE_SILL_HEIGTH',                             &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR130'

            if (Me%ReferentialZ == Hydrographic_) then
                NewDischarge%Valve%SillHeigth =  - NewDischarge%Valve%SillHeigth
            endif                
            

            if (flag /= 1) then

                call GetData(NewDischarge%Valve%AxisHeigth,                                 &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             FromBlock,                                                     &
                             keyword      ='VALVE_AXIS_HEIGTH',                             &
                             ClientModule = 'ModuleDischarges',                             &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR135'

                if (Me%ReferentialZ == Hydrographic_) then
                    NewDischarge%Valve%AxisHeigth =  - NewDischarge%Valve%AxisHeigth
                endif  
                
                if      (NewDischarge%Valve%SectionType == circular_area    ) then
                
                    NewDischarge%Valve%SillHeigth  = NewDischarge%Valve%AxisHeigth - NewDischarge%Valve%Diameter / 2.
                    
                elseif  (NewDischarge%Valve%SectionType == rectangular_area ) then
                
                    NewDischarge%Valve%SillHeigth  = NewDischarge%Valve%AxisHeigth - NewDischarge%Valve%Height / 2.
                
                endif                    
            
                if (flag /= 1) then
                    write(*,*) 'Valve axis Missing'
                    stop ' Construct_FlowValues - ModuleDischarges - ERR138'
                endif
                                    
            endif
            
            call GetData(NewDischarge%Valve%PipeLength,                                 &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='PIPE_LENGTH',                                   &
                         ClientModule = 'ModuleDischarges',                             &
                         default      = 0.,                                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR136'

            !manning coefficient 
            call GetData(NewDischarge%Valve%PipeManning,                                &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='PIPE_MANNING',                                  &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR137'      
            
            if (NewDischarge%Valve%PipeLength > 0. .and. flag == 0) then
                write(*,*) 'need to define PIPE_MANNING'
                stop 'Construct_FlowValues - ModuleDischarges - ERR138'
            endif     

        else if (NewDischarge%DischargeType == RatingCurve) then i2
            
            !Get Block with rating curves values
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,                   &
                            beginratingcurve, endratingcurve,                           &
                            BlockLayersFound,                                           &
                            FirstLine = FirstLine,                                      &
                            LastLine  = LastLine,                                       &
                            STAT      = STAT_CALL)

cd1 :       if (STAT_CALL .EQ. SUCCESS_  .and. BlockLayersFound) then

    
                    allocate(NewDischarge%RatingCurve%Level(LastLine - FirstLine -1))
                    allocate(NewDischarge%RatingCurve%Flow(LastLine - FirstLine -1))
                    NewDischarge%RatingCurve%nValues = LastLine - FirstLine -1
                    allocate(BufferLine(2))
                    
                    iValue = 1;                   
                    do  iLine = FirstLine+1, LastLine-1

                        call GetData(BufferLine,                                            &
                                     Me%ObjEnterData,                                       &
                                     flag,                                                  &
                                     Buffer_Line = iLine,                                   &
                                     STAT = STAT_CALL) 
                        if (STAT_CALL /= SUCCESS_ .or. flag /= 2)             &
                            stop "Read Rating Curve Values - ModuleDischarges - ERR139"

                        NewDischarge%RatingCurve%Level(iValue) = BufferLine (1)
                        NewDischarge%RatingCurve%Flow(iValue) = BufferLine (2)
                    
                        iValue = iValue + 1
                        
                    enddo

                    deallocate (BufferLine)
                    
            else 

                stop "Read Rating Curve Values - ModuleDischarges - ERR140"

            endif cd1            
            
            
            
        endif i2

        call GetData(NewDischarge%ByPass%ON,                                            &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     FromBlock,                                                         &
                     keyword      ='BYPASS_ON',                                         &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleDischarges',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR140'
        
        if (NewDischarge%DischargeType == Valve .and. .not. NewDischarge%ByPass%ON) then
            write(*,*) 'In the case of a type "valve" discharge the discharge must also be "bypass"'
            stop 'Construct_FlowValues - ModuleDischarges - ERR145'                    
        endif

i3:     if (NewDischarge%ByPass%ON) then

i4:         if (NewDischarge%Localization%CoordinatesON) then

                call GetData(NewDischarge%ByPass%X,                                     &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             FromBlock,                                                 &
                             keyword      ='BYPASS_X',                                  &
                             ClientModule = 'ModuleDischarges',                         &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR150'

                if (flag /= 1) then
                    write(*,*) 'Bypass I cell missing'
                    stop ' Construct_FlowValues - ModuleDischarges - ERR160'
                endif

                call GetData(NewDischarge%ByPass%Y,                                     &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             FromBlock,                                                 &
                             keyword      ='BYPASS_Y',                                  &
                             ClientModule = 'ModuleDischarges',                         &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR160'

                if (flag /= 1) then
                    write(*,*) 'Bypass J cell missing'
                    stop ' Construct_FlowValues - ModuleDischarges - ERR170'
                endif


            else i4
            
                call GetData(NewDischarge%ByPass%i,                                     &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             FromBlock,                                                 &
                             keyword      ='BYPASS_I',                                  &
                             ClientModule = 'ModuleDischarges',                         &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR180'

                if (flag /= 1) then
                    write(*,*) 'Bypass I cell missing'
                    stop ' Construct_FlowValues - ModuleDischarges - ERR190'
                endif

                call GetData(NewDischarge%ByPass%j,                                     &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             FromBlock,                                                 &
                             keyword      ='BYPASS_J',                                  &
                             ClientModule = 'ModuleDischarges',                         &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR160'

                if (flag /= 1) then
                    write(*,*) 'Bypass J cell missing'
                    stop ' Construct_FlowValues - ModuleDischarges - ERR200'
                endif

            endif i4
            
            
            call GetData(NewDischarge%ByPass%k,                                         &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='BYPASS_K',                                      &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR210'

            if (flag /= 1) then
                write(*,*) 'Bypass K cell missing'
                stop ' Construct_FlowValues - ModuleDischarges - ERR220'
            endif

            call GetData(NewDischarge%ByPass%Side,                                      &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='BYPASS_SIDE',                                   &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR230'

            if (flag /= 1) then
                write(*,*) 'Bypass side missing'
                stop ' Construct_FlowValues - ModuleDischarges - ERR240'
            endif

            if      (NewDischarge%ByPass%Side /= SideA   .and.                          &
                     NewDischarge%ByPass%Side /= SideB) then
                     stop 'Construct_FlowValues - ModuleDischarges - ERR250'
            endif

            call GetData(NewDischarge%ByPass%OneWay,                                    &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='BYPASS_ONEWAY',                                 &
                         default      = .false.,                                        &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR260'

        endif i3         

        call GetData(NewDischarge%FromIntake%ON,                                        &
                     Me%ObjEnterData, flag,                                             &
                     FromBlock,                                                         &
                     keyword      = 'FROM_INTAKE',                                      &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleDischarges',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR270'

        if(NewDischarge%FromIntake%ON)then

            call GetData(NewDischarge%FromIntake%IntakeName,                            &
                         Me%ObjEnterData, flag,                                         &
                         FromBlock,                                                     &
                         keyword      = 'INTAKE_NAME',                                  &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR280'

            if(flag == 0)then
                write(*,*)"Must define INTAKE_NAME"
                write(*,*)"in discharge ", trim(adjustl(NewDischarge%ID%Name))
                stop 'Construct_FlowValues - ModuleDischarges - ERR290'
            endif

            call GetData(NewDischarge%FromIntake%AssociateFlow,                         &
                         Me%ObjEnterData, flag,                                         &
                         FromBlock,                                                     &
                         keyword      = 'FLOW_FROM_INTAKE',                             &
                         default      = .false.,                                        &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR300'


            if(NewDischarge%FromIntake%AssociateFlow)then

                call GetData(NewDischarge%FromIntake%FlowFraction,                      &
                             Me%ObjEnterData, flag,                                     &
                             FromBlock,                                                 &
                             keyword      = 'FLOW_FRACTION',                            &
                             default      = 1.0,                                        &
                             ClientModule = 'ModuleDischarges',                         &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR310'

                if(NewDischarge%FromIntake%FlowFraction < 0)then

                    write(*,*)"Discharge flow percentage of the water intake cannot be negative"
                    write(*,*)"Discharge name ", trim(adjustl(NewDischarge%ID%Name))
                
                end if

            end if

        end if
        
        if (NewDischarge%DischargeType == OpenMILink) then
        
            NewDischarge%WaterFlow%Variable = .true.
        
        endif
        !Joao Sobrinho
        call GetData(NewDischarge%WaterFlow%Upscaling,                         &
                Me%ObjEnterData,                                               &
                flag,                                                          &
                FromBlock,                                                     &
                keyword      ='UPSCALING_DISCHARGE',                           &
                ClientModule = 'ModuleDischarges',                             &
                default      = .false.,                                        &
                STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR320'
        


     End Subroutine Construct_FlowValues  

    !--------------------------------------------------------------------------


     Subroutine Construct_VariableLocation(NewDischarge)

        !Arguments--------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: NewDischarge

        !External-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: flag

        !----------------------------------------------------------------------

        if (NewDischarge%TimeSerieON) then
            !Searches for the X location column (if the discharge as associated an input time serie)
            call GetData(NewDischarge%Localization%XColumn,                              &
                         Me%ObjEnterData, flag,                                          &
                         keyword    = 'X_COLUMN',                                        &
                         default    = FillValueInt,                                      &
                         SearchType = FromBlock,                                         &
                         ClientModule = 'ModuleDischarges',                              &
                         STAT       = STAT_CALL)
                         
            if (STAT_CALL /= SUCCESS_) stop 'Construct_VariableLocation - ModuleDischarges - ERR10'

            if (flag == 1) then

                NewDischarge%Localization%VariableX = .true. 
                
                if (NewDischarge%Localization%XColumn > NewDischarge%TimeSerieMaxCol .or.   &
                    NewDischarge%Localization%XColumn < 2) then
                    
                    write(*,*)'Look at file', trim(Me%DataFile)
                    write(*,*)'X_COLUMN not valid in discharge ', trim(NewDischarge%ID%Name)
                    stop 'Construct_VariableLocation - ModuleDischarges - ERR20'
                    
                endif
                           
            else
                NewDischarge%Localization%VariableX = .false.            
            endif
            
            !Searches for the Y location column (if the discharge as associated an input time serie)
            call GetData(NewDischarge%Localization%YColumn,                              &
                         Me%ObjEnterData, flag,                                          &
                         keyword    = 'Y_COLUMN',                                        &
                         default    = FillValueInt,                                      &
                         SearchType = FromBlock,                                         &
                         ClientModule = 'ModuleDischarges',                              &
                         STAT       = STAT_CALL)
                         
            if (STAT_CALL /= SUCCESS_) stop 'Construct_VariableLocation - ModuleDischarges - ERR30'

            if (flag == 1) then

                NewDischarge%Localization%VariableY = .true.            
                
                if (NewDischarge%Localization%YColumn > NewDischarge%TimeSerieMaxCol .or.   &
                    NewDischarge%Localization%YColumn < 2) then
                    
                    write(*,*)'Look at file', trim(Me%DataFile)
                    write(*,*)'Y_COLUMN not valid in discharge ', trim(NewDischarge%ID%Name)
                    stop 'Construct_VariableLocation - ModuleDischarges - ERR40'
                    
                endif
                
            else
                NewDischarge%Localization%VariableY = .false.            
            endif
            
        endif

        !----------------------------------------------------------------------

    end subroutine Construct_VariableLocation


    !--------------------------------------------------------------------------

     Subroutine Construct_VelocityValues(NewDischarge)

        !Arguments--------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: NewDischarge

        !External-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: flag
        real, dimension(:), pointer                 :: Aux 

        !----------------------------------------------------------------------


        !By default the property value in the domain is zero
        allocate(aux(3))
        call GetData(Aux, Me%ObjEnterData, flag,                                         &
                     keyword    = 'DEFAULT_VELOCITY_VALUE',                              &
                     default    = 0.,                                                    &
                     SearchType = FromBlock,                                             &
                     ClientModule = 'ModuleDischarges',                                  &
                     STAT       = STAT_CALL)
                     
        if   (flag /= 0 .and. flag /=2 .and. flag /=3 ) then    
            stop 'Construct_VelocityValues - ModuleDischarges - ERR10'
        endif                 

        NewDischarge%VelocityFlow%UScalar = Aux(1)
        NewDischarge%VelocityFlow%VScalar = Aux(2)
        if     (flag == 2) then
            NewDischarge%VelocityFlow%WScalar = 0.
        elseif (flag == 0 .or. flag == 3) then
            NewDischarge%VelocityFlow%WScalar = Aux(3)
        endif
        deallocate(aux)

        !Searches for the velocity U column (if the discharge is read from a data base)
        if (NewDischarge%TimeSerieON) then
            call GetData(NewDischarge%VelocityFlow%UColumn,                             &
                         Me%ObjEnterData, flag,                                         &
                         keyword        = 'U_COLUMN',                                   &
                         default        = FillValueInt,                                 &
                         SearchType     = FromBlock,                                    &
                         ClientModule   = 'ModuleDischarges',                           &
                         STAT           = STAT_CALL)
                         
            if (STAT_CALL /= SUCCESS_) then
                stop 'Construct_VelocityValues - ModuleDischarges - ERR20'
            endif                

            if (flag == 1) then

                NewDischarge%VelocityFlow%UVariable = .true.            
                
                if (NewDischarge%VelocityFlow%UColumn > NewDischarge%TimeSerieMaxCol .or. &
                    NewDischarge%VelocityFlow%UColumn < 2) then
                    
                    write(*,*)'Look at file', trim(Me%DataFile)
                    write(*,*)'U_COLUMN not valid in discharge ', trim(NewDischarge%ID%Name)
                    stop 'Construct_VelocityValues - ModuleDischarges - ERR30'
                    
                endif                
                
            else
                NewDischarge%VelocityFlow%UVariable = .false.            
            endif

            call GetData(NewDischarge%VelocityFlow%VColumn,                              &
                         Me%ObjEnterData, flag,                                          &
                         keyword    = 'V_COLUMN',                                        &
                         default    = FillValueInt,                                      &
                         SearchType = FromBlock,                                         &
                         ClientModule = 'ModuleDischarges',                              &
                         STAT       = STAT_CALL)
                         
            if (STAT_CALL /= SUCCESS_) then
                stop 'Construct_VelocityValues - ModuleDischarges - ERR40'
            endif                              

            if (flag == 1) then
                NewDischarge%VelocityFlow%VVariable = .true.
                
                if (NewDischarge%VelocityFlow%VColumn > NewDischarge%TimeSerieMaxCol .or. &
                    NewDischarge%VelocityFlow%VColumn < 2) then
                    
                    write(*,*)'Look at file', trim(Me%DataFile)
                    write(*,*)'V_COLUMN not valid in discharge ', trim(NewDischarge%ID%Name)
                    stop 'Construct_VelocityValues - ModuleDischarges - ERR50'
                    
                endif                
                
            else
                NewDischarge%VelocityFlow%VVariable = .false.
            endif

            call GetData(NewDischarge%VelocityFlow%WColumn,                              &
                         Me%ObjEnterData, flag,                                          &
                         keyword    = 'W_COLUMN',                                        &
                         default    = FillValueInt,                                      &
                         SearchType = FromBlock,                                         &
                         ClientModule = 'ModuleDischarges',                              &
                         STAT       = STAT_CALL)
                         
            if (STAT_CALL /= SUCCESS_) then
                stop 'Construct_VelocityValues - ModuleDischarges - ERR60'
            endif                           

            if (flag == 1) then
                NewDischarge%VelocityFlow%WVariable = .true.
                
                if (NewDischarge%VelocityFlow%WColumn > NewDischarge%TimeSerieMaxCol .or. &
                    NewDischarge%VelocityFlow%WColumn < 2) then
                    
                    write(*,*)'Look at file', trim(Me%DataFile)
                    write(*,*)'W_COLUMN not valid in discharge ', trim(NewDischarge%ID%Name)
                    stop 'Construct_VelocityValues - ModuleDischarges - ERR70'
                    
                endif                          
            else
                NewDischarge%VelocityFlow%WVariable = .false.
            endif

        endif
        
        !----------------------------------------------------------------------

    end subroutine Construct_VelocityValues

    !--------------------------------------------------------------------------

    subroutine Construct_PropertyList(NewDischarge, ClientNumber)

        !Arguments--------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: NewDischarge
        integer, intent(IN)                         :: ClientNumber

        !External--------------------------------------------------------------
        logical                                     :: BlockFound
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: NewProperty

        !----------------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,            &
                                       beginproperty, endproperty, BlockFound,   &
                                       STAT = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_      ) then    
cd2 :           if (BlockFound) then                                                  
                    ! Construct a New Property 
                    call Construct_Property(NewDischarge,NewProperty)

                    ! Add new Property to the WaterProperties List 
                    call Add_Property(NewDischarge,NewProperty)
                else
                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                stop 'Subroutine Construct_PropertyList - ModuleDischarges. ERR02.'
            end if cd1
        end do do1

        !----------------------------------------------------------------------

    end subroutine Construct_PropertyList

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct a new property.           

    subroutine Construct_Property(NewDischarge,NewProperty)

        !Arguments--------------------------------------------------------------
        type(T_property), pointer                   :: NewProperty
        type(T_IndividualDischarge), pointer        :: NewDischarge
        !----------------------------------------------------------------------

        allocate(NewProperty)

        nullify (NewProperty%Prev,NewProperty%Next)


        !Construct property ID
        call ConstructPropertyID     (NewProperty%ID, Me%ObjEnterData, FromBlockInBlock)


        !Construct property values
        call Construct_PropertyValues(NewDischarge, NewProperty)

        !----------------------------------------------------------------------

    end subroutine Construct_Property

    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct the property values       
    ! in the domain and in the boundaries            

    subroutine Construct_PropertyValues(NewDischarge,NewProperty)

        !Arguments--------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: NewDischarge
        type(T_property),   pointer                 :: NewProperty

        !External-----------------------------------------------------------------
        character(Len = PathLength)                 :: FileName
        integer                                     :: flag, STAT_CALL

        !----------------------------------------------------------------------
   
        call GetData(NewProperty%ConcColumn,                                            &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     FromBlockinBlock,                                                  &
                     keyword      ='TIME_SERIE_COLUMN',                                 &
                     ClientModule ='ModuleDischarges')

        if (flag == 1) then

            NewProperty%Variable = .true.
            
        else
            NewProperty%Variable = .false.
        endif

ifvar:  if (NewProperty%Variable) then

            call GetData(FileName,                                                       &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         FromBlockinBlock,                                               &
                         keyword      ='FILENAME',                                       &
                         ClientModule ='ModuleDischarges')

            if (flag == 1) then
                NewProperty%PropTimeSerie = .true.
            else
                NewProperty%PropTimeSerie = .false.
            endif

            if (NewProperty%PropTimeSerie) then

                call StartTimeSerieInput(NewProperty%TimeSerie, FileName,  Me%ObjTime, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleDischarges - ERR10'
                
                call GetTimeSerieDataColumns (TimeSerieID = NewProperty%TimeSerie,      &
                                              DataColumns = NewProperty%TimeSerieMaxCol,&
                                              STAT        = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleDischarges - ERR20'
                
                if (NewProperty%ConcColumn > NewProperty%TimeSerieMaxCol .or.               &
                    NewProperty%ConcColumn < 2) then
                    
                    write(*,*)'Look at file', trim(Me%DataFile)
                    write(*,*)'TIME_SERIE_COLUMN not valid in dicharge ', trim(NewDischarge%ID%Name)
                    write(*,*) 'Property name =',trim(NewProperty%ID%Name)
                    stop 'Construct_PropertyValues - ModuleDischarges - ERR30'
                    
                endif                     
                

            else
                if (NewDischarge%TimeSerieON) then
                    NewProperty%TimeSerie = NewDischarge%TimeSerie
                else
                    stop 'Construct_PropertyValues - ModuleDischarges - ERR40'
                endif
            endif

        endif ifvar


        call GetData(NewProperty%FromIntake,                                            &
                     Me%ObjEnterData, flag,                                             &
                     FromBlockInBlock,                                                  &
                     keyword      = 'PROP_FROM_INTAKE',                                 &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleDischarges',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleDischarges - ERR50'


        if(NewProperty%FromIntake)then

            if(.not. NewDischarge%FromIntake%ON)then
                
                write(*,*)"Property concentration cannot be determined from intake"
                write(*,*)"because discharge is not based on intake"
                write(*,*)"Property  : ", trim(adjustl(NewProperty%ID%Name))
                write(*,*)"Discharge : ", trim(adjustl(NewDischarge%ID%Name))
                stop      'Construct_PropertyValues - ModuleDischarges - ERR60'

            end if
            

        endif
        
        if (NewProperty%FromIntake .or. NewDischarge%ByPass%ON) then

            call GetData(NewProperty%IncreaseValue,                                     &
                         Me%ObjEnterData, flag,                                         &
                         FromBlockInBlock,                                              &
                         keyword      = 'INCREASE_VALUE',                               &
                         default      = 0.0,                                            &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if(STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleDischarges - ERR70'
            
            if(flag .eq. 0 .and. NewProperty%FromIntake)then

                write(*,*)"Discharge from intake property. Please define INCREASE_VALUE for :"
                write(*,*)trim(NewProperty%ID%Name)
                write(*,*)"Discharge name : ", trim(NewDischarge%ID%Name)
                stop 'Construct_PropertyValues - ModuleDischarges -  ERR80'

            end if 

        end if

        !By default the property value in the domain is zero
        call GetData(NewProperty%scalar,                                                &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     FromBlockinBlock,                                                  &
                     keyword      ='DEFAULTVALUE',                                      &
                     ClientModule ='ModuleDischarges',                                  &
                     default      = 0.0)
                     
        if((.not. NewProperty%Variable) .and. (.not. NewProperty%FromIntake) .and. (flag .eq. 0))then

            write(*,*)"Please define concentration default value for property :"
            write(*,*)trim(NewProperty%ID%Name)
            write(*,*)"Discharge name : ", trim(NewDischarge%ID%Name)
            stop 'Construct_PropertyValues - ModuleDischarges -  ERR90'

        end if

        !----------------------------------------------------------------------

    end subroutine Construct_PropertyValues

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    ! This subroutine adds a new property to the Water Property List  

    subroutine Add_Property(NewDischarge,NewProperty)

        !Arguments 
          
        type(T_IndividualDischarge), pointer    :: NewDischarge
        type(T_property),            pointer    :: NewProperty


        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(NewDischarge%FirstProperty)) then
            NewDischarge%PropertiesNumber = 1
            NewDischarge%FirstProperty    => NewProperty
            NewDischarge%LastProperty     => NewProperty
        else
            NewProperty%Prev                     => NewDischarge%LastProperty
            NewDischarge%LastProperty%Next => NewProperty
            NewDischarge%LastProperty      => NewProperty
            NewDischarge%PropertiesNumber  = &
                   NewDischarge%PropertiesNumber + 1
        end if 

        !----------------------------------------------------------------------

    end subroutine Add_Property 

    !--------------------------------------------------------------------------

    subroutine ConstructIntakeDischarges

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: CurrentDischarge
        type(T_IndividualDischarge), pointer        :: Intake
        integer                                     :: STAT_CALL

        CurrentDischarge => Me%FirstDischarge
        do while (associated (CurrentDischarge))

            if(CurrentDischarge%FromIntake%ON)then

                call Search_Discharge_ByName(Intake, STAT_CALL, trim(adjustl(CurrentDischarge%FromIntake%IntakeName)))
                if (STAT_CALL/=SUCCESS_) then 
                    write(*,*)'Can not find discharge with name ', trim(adjustl(CurrentDischarge%FromIntake%IntakeName)), '.'
                    stop      'Subroutine ConstructIntakeDischarges - ModuleDischarges. ERR01.'
                else
                    CurrentDischarge%FromIntake%IntakeID = Intake%ID%IDNumber
                    
                    !if intake flow is variable, and discharge depends on intake flow
                    !then it must also be variable
                    if(CurrentDischarge%FromIntake%AssociateFlow)then
                        CurrentDischarge%WaterFlow%Variable = Intake%WaterFlow%Variable
                    endif
                    
                endif

            end if


            CurrentDischarge => CurrentDischarge%Next
        enddo


    end subroutine ConstructIntakeDischarges

    !--------------------------------------------------------------------------

    subroutine ConstructLog

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
#ifndef _OUTPUT_OFF_
        type(T_IndividualDischarge), pointer        :: CurrentDischarge

        write(*, *)"----------------------- DISCHARGES -----------------------"
        write(*, *)
        write(*, *)"Number of Discharges : ", Me%DischargesNumber
        write(*, *)

        CurrentDischarge => Me%FirstDischarge
        do while (associated (CurrentDischarge))
            write(*, *)"----Discharge        : ", trim(adjustl(CurrentDischarge%ID%Name))
            write(*, *)"----Num of Properties: ", max(CurrentDischarge%PropertiesNumber, 0)
            CurrentDischarge => CurrentDischarge%Next
        enddo
#endif

    end subroutine ConstructLog



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                            
    !--------------------------------------------------------------------------

    Subroutine GetDischargesNumber(DischargesID, DischargesNumber, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: DischargesID
        integer,           intent(OUT)              :: DischargesNumber
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_
 
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then
           
            DischargesNumber = Me%DischargesNumber

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end Subroutine GetDischargesNumber

    !--------------------------------------------------------------------------

    subroutine GetDischargesGridLocalization(DischargesID, DischargeIDNumber, Igrid,    &
                                             JGrid, KGrid, IByPass, JByPass, KByPass,   &
                                             DischVertical, KDepth, WaterColumnZ,       &
                                             Bathymetry, OpenPoints3D, CoordinateX,     &
                                             CoordinateY, CoordinatesON, TimeX,         &
                                             XByPass, YByPass, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: DischargesID
        integer,           intent(IN )                  :: DischargeIDNumber
        integer, optional, intent(OUT)                  :: IGrid, JGrid, KGrid
        integer, optional, intent(OUT)                  :: IByPass, JByPass, KByPass
        integer, optional, intent(OUT)                  :: DischVertical
        integer, optional, dimension(:, :, :), pointer  :: OpenPoints3D
        real   , optional, dimension(:, :   ), pointer  :: WaterColumnZ
        real   , optional, dimension(:, :   ), pointer  :: Bathymetry
        
        real,    optional, intent(OUT)                  :: CoordinateX, CoordinateY, Kdepth
        logical, optional, intent(OUT)                  :: CoordinatesON        
        type (T_Time), optional, intent (IN)            :: TimeX
        real,    optional, intent(OUT)                  :: XByPass, YByPass
                
        integer, optional, intent(OUT)                  :: STAT

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL
        type(T_IndividualDischarge), pointer            :: DischargeX

        !Local-----------------------------------------------------------------
        integer                                         :: ready_         
        integer                                         :: STAT_
        integer                                         :: WorkSizeILB, WorkSizeIUB
        integer                                         :: WorkSizeJLB, WorkSizeJUB        

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeIDNumber, '.'
                stop       'Subroutine GetDischargesGridLocalization - ModuleDischarges. ERR01.'
            endif

            !if (.not. DischargeX%Localization%Location2D) then
            !    write(*,*)'Discharge location not given in Grid Coordinates'
            !    write(*,*)trim(adjustl(DischargeX%ID%Name))
            !    stop 'GetDischargesGridLocalization - ModuleDischarges - ERR01a'
            !endif
            
            if (DischargeX%Localization%VariableX .and. present(TimeX)) then
                DischargeX%Localization%CoordinateX  =                                  &
                    TimeSerieValue(DischargeX%TimeSerie, DischargeX%UseOriginalValues,  &
                                   TimeX, DischargeX%Localization%XColumn)            
            endif

            if (DischargeX%Localization%VariableY .and. present(TimeX)) then
                DischargeX%Localization%CoordinateY  =                                  &
                    TimeSerieValue(DischargeX%TimeSerie, DischargeX%UseOriginalValues,  &
                                   TimeX, DischargeX%Localization%YColumn)            
            endif

            if (present(IGrid           )) IGrid              = DischargeX%Localization%GridCoordinates%I
            if (present(JGrid           )) JGrid              = DischargeX%Localization%GridCoordinates%J
            

            if (present(CoordinateX     )) CoordinateX        = DischargeX%Localization%CoordinateX
            if (present(CoordinateY     )) CoordinateY        = DischargeX%Localization%CoordinateY

            if (present(CoordinatesON   )) CoordinatesON      = DischargeX%Localization%CoordinatesON

            if (present(DischVertical   )) DischVertical      = DischargeX%Localization%DischVertical
            if (present(KGrid           )) KGrid              = DischargeX%Localization%GridCoordinates%K
            if (present(Kdepth          )) Kdepth             = DischargeX%Localization%Kdepth

            if (present(IByPass         )) IByPass            = DischargeX%ByPass%I
            if (present(JByPass         )) JByPass            = DischargeX%ByPass%J
            if (present(KByPass         )) KByPass            = DischargeX%ByPass%K
            if (present(XByPass         )) XByPass            = DischargeX%ByPass%X
            if (present(YByPass         )) YByPass            = DischargeX%ByPass%Y


            !If OpenPoints is present an DischargeX have alternative locations, 
            !check out which locations to use.
            if (present(WaterColumnZ)) then

                WorkSizeILB = lbound(WaterColumnZ, dim = 1) + 1
                WorkSizeIUB = ubound(WaterColumnZ, dim = 1) - 1

                WorkSizeJLB = lbound(WaterColumnZ, dim = 2) + 1
                WorkSizeJUB = ubound(WaterColumnZ, dim = 2) - 1

                if (IGrid < WorkSizeILB .or. IGrid > WorkSizeIUB)                       &
                    stop  'GetDischargesGridLocalization - ModuleDischarges. ERR02.'

                if (JGrid < WorkSizeJLB .or. JGrid > WorkSizeJUB)                       &
                    stop  'GetDischargesGridLocalization - ModuleDischarges. ERR03.'
                
                if (WaterColumnZ(IGrid, JGrid) <                                        &
                    DischargeX%Localization%MinimumDischargeDepth .and.                 &
                    DischargeX%Localization%AlternativeLocations) then                    
                                        
                    
                    if (.not. DischargeX%Localization%StartFromLastDischargeLocation .or. &
                        DischargeX%Localization%GridCoordinates%OldI == FillValueInt) then
                        call FindDischargePointFromStart (DischargeX, OpenPoints3D, WaterColumnZ, Bathymetry, IGrid, JGrid,  &
                                                            WorkSizeILB, WorkSizeIUB, WorkSizeJLB, WorkSizeJUB  )
                    else
                        
                        call FindDischargePointFromLast  (DischargeX, OpenPoints3D, WaterColumnZ, Bathymetry, IGrid, JGrid,  &
                                                            WorkSizeILB, WorkSizeIUB, WorkSizeJLB, WorkSizeJUB  )

                    endif
                        
                    DischargeX%Localization%GridCoordinates%OldI = IGrid
                    DischargeX%Localization%GridCoordinates%OldJ = JGrid
                                                
                    if (DischargeX%Localization%TrackLocation) then
                        !escrever discharge location
                        write (DischargeX%Localization%TrackLocationFileUnitNumber,'(A25, I10, I10, F10.5, F10.5, F10.5)')  &
                            DischargeX%ID%Name, IGrid, JGrid,  WaterColumnZ(IGrid, JGrid), WaterColumnZ(IGrid+1, JGrid),    &
                            WaterColumnZ(IGrid+2, JGrid)
                    endif

                endif
            endif
        
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDischargesGridLocalization

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine GetByPassON(DischargesID, DischargeIDNumber, ByPassON, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: DischargesID
        integer,           intent(IN )                  :: DischargeIDNumber
        logical,           intent(OUT)                  :: ByPassON
        integer, optional, intent(OUT)                  :: STAT

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge), pointer            :: DischargeX
        integer                                         :: STAT_CALL, STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeIDNumber, '.'
                stop       'Subroutine GetDischargesGridLocalization - ModuleDischarges. ERR01.'
            endif

            ByPassON = DischargeX%ByPass%ON
                    
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetByPassON


    !--------------------------------------------------------------------------

    subroutine GetDischargeFromIntakeON(DischargesID, DischargeIDNumber, DischargeFromIntakeON, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: DischargesID
        integer,           intent(IN )                  :: DischargeIDNumber
        logical,           intent(OUT)                  :: DischargeFromIntakeON
        integer, optional, intent(OUT)                  :: STAT

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge), pointer            :: DischargeX
        integer                                         :: STAT_CALL, STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeIDNumber, '.'
                stop       'Subroutine GetDischargesGridLocalization - ModuleDischarges. ERR01.'
            endif

            DischargeFromIntakeON = DischargeX%FromIntake%ON
                    
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDischargeFromIntakeON


    !--------------------------------------------------------------------------

    
    !Starts at the initial discharge point and searches for nearby locations that are open points and 
    !have a minimum watercolumn
    subroutine FindDischargePointFromStart (DischargeX, OpenPoints3D, WaterColumnZ, Bathymetry, newI, newJ, &
                                            WorkSizeILB, WorkSizeIUB, WorkSizeJLB, WorkSizeJUB)
        
        !Arguments-------------------------------------------------------------                
        integer, optional, intent(OUT)                  :: newI, newJ        
        integer, optional, dimension(:, :, :), pointer  :: OpenPoints3D
        real   , optional, dimension(:, :   ), pointer  :: WaterColumnZ
        real   , optional, dimension(:, :   ), pointer  :: Bathymetry 
                       
        type(T_IndividualDischarge), pointer            :: DischargeX

        integer                                         :: WorkSizeILB, WorkSizeIUB
        integer                                         :: WorkSizeJLB, WorkSizeJUB

        !Local-----------------------------------------------------------------
        logical, dimension(:,:), allocatable            :: PossibleDischargePoint
        integer                                         :: MoveDirection
        logical                                         :: NewLocationFound
        integer, parameter                              :: NONE_    = 1
        integer, parameter                              :: NORTH_   = 2
        integer, parameter                              :: EAST_    = 3
        integer, parameter                              :: SOUTH_   = 4
        integer, parameter                              :: WEST_    = 5
        real                                            :: Depth, DepthCurrent
        real                                            :: Bottom, BtomCurrent
        integer                                         :: i,j,k
        
        !----------------------------------------------------------------------

        !Allocates temporary Matrix of possible discharge points
        allocate(PossibleDischargePoint(WorkSizeILB:WorkSizeIUB, WorkSizeJLB:WorkSizeJUB))
        PossibleDischargePoint = .true.
                    
        !Search the lowest point in the sourring points        
        NewLocationFound = .false.
        i = DischargeX%Localization%GridCoordinates%I
        j = DischargeX%Localization%GridCoordinates%J
        k = DischargeX%Localization%GridCoordinates%K
                    
        DepthCurrent  = WaterColumnZ(i, j)
        BtomCurrent   = Bathymetry   (i, j)

        PossibleDischargePoint(i, j) = .false.

        do while (.not. NewLocationFound)

            MoveDirection = NONE_                    
                            
cd1 :       if (OpenPoints3D(i,j,k) == OpenPoint) then 
                !neste caso compara-se as colunas de agua 
                                
                !Point in the North
                if (i < WorkSizeIUB) then                            
                    Depth = WaterColumnZ(i+1, j)                               
                    if (Depth >= DepthCurrent .and. PossibleDischargePoint(i+1, j)) then
                        MoveDirection = NORTH_
                        DepthCurrent  = Depth
                    endif

                endif

                !Point in the East
                if (j < WorkSizeJUB) then
                    Depth = WaterColumnZ(i, j+1)
                    if (Depth >= DepthCurrent .and. PossibleDischargePoint(i, j+1)) then
                        MoveDirection = EAST_
                        DepthCurrent  = Depth
                    endif
                endif

               !Point in the South
                if (i > WorkSizeILB) then
                    Depth = WaterColumnZ(i-1, j)
                    if (Depth >= DepthCurrent .and. PossibleDischargePoint(i-1, j)) then
                        MoveDirection = SOUTH_
                        DepthCurrent  = Depth
                    endif
                endif

                !Point in the West
                if (j > WorkSizeJLB) then
                    Depth = WaterColumnZ(i, j-1)
                    if (Depth >= DepthCurrent .and. PossibleDischargePoint(i, j-1)) then
                        MoveDirection = WEST_
                        DepthCurrent  = Depth
                    endif
                endif

            else cd1
            !neste caso comparam-se as cotas (os openpoints n ficam todos bem com a minWC)
                                
                !Point in the North
                if (i < WorkSizeIUB) then                            
                    Bottom = Bathymetry(i+1, j)                               
                    if (Bottom >= BtomCurrent .and. PossibleDischargePoint(i+1, j)) then
                        MoveDirection = NORTH_
                        BtomCurrent   = Bottom
                    endif
                endif

                !Point in the East
                if (j < WorkSizeJUB) then
                    Bottom = Bathymetry(i, j+1)
                    if (Bottom >= BtomCurrent .and. PossibleDischargePoint(i, j+1)) then
                        MoveDirection = EAST_
                        BtomCurrent   = Bottom
                    endif
                endif

                !Point in the South
                if (i > WorkSizeILB) then
                    Bottom = Bathymetry(i-1, j)
                    if (Bottom >= BtomCurrent .and. PossibleDischargePoint(i-1, j)) then
                        MoveDirection = SOUTH_
                        BtomCurrent   = Bottom
                    endif
                endif

                !Point in the West
                if (j > WorkSizeJLB) then
                    Bottom = Bathymetry (i, j-1)
                    if (Bottom >= BtomCurrent .and. PossibleDischargePoint(i, j-1)) then
                        MoveDirection = WEST_
                        BtomCurrent   = Bottom
                    endif
                endif

            endif cd1

                                            
            select case (MoveDirection)

                case (NORTH_)
                    i = i + 1
                    PossibleDischargePoint(i, j) = .false.

                case (EAST_ )
                    j = j + 1
                    PossibleDischargePoint(i, j) = .false.

                case (SOUTH_)
                    i = i - 1
                    PossibleDischargePoint(i, j) = .false.

                case (WEST_ )
                    j = j - 1
                    PossibleDischargePoint(i, j) = .false.

                case (NONE_)
                    !a descarga n pode ir para sitio nenhum fica no mm sitio
                    NewLocationFound = .true.
                    newI            = i
                    newJ            = j

            end select

            !Checks out if the new defined point is an openpoint
            if (WaterColumnZ(i, j) > DischargeX%Localization%MinimumDischargeDepth) then
                NewLocationFound = .true.
                newI            = i
                newJ            = j                         
            endif

        enddo
        
        !Deallocates temporary matrix
        deallocate(PossibleDischargePoint)

    end subroutine FindDischargePointFromStart

    !--------------------------------------------------------------------------
        
    subroutine FindDischargePointFromLast (DischargeX, OpenPoints3D, WaterColumnZ, Bathymetry, newI, newJ, &
                                            WorkSizeILB, WorkSizeIUB, WorkSizeJLB, WorkSizeJUB)
        
        !Arguments-------------------------------------------------------------                
        integer, optional, intent(OUT)                  :: newI, newJ        
        integer, optional, dimension(:, :, :), pointer  :: OpenPoints3D
        real   , optional, dimension(:, :   ), pointer  :: WaterColumnZ
        real   , optional, dimension(:, :   ), pointer  :: Bathymetry 
                       
        type(T_IndividualDischarge), pointer            :: DischargeX

        integer                                         :: WorkSizeILB, WorkSizeIUB
        integer                                         :: WorkSizeJLB, WorkSizeJUB

        !Local-----------------------------------------------------------------
        logical, dimension(:,:), allocatable            :: PossibleDischargePoint
        integer                                         :: MoveDirection
        logical                                         :: NewLocationFound
        integer, parameter                              :: NONE_    = 1
        integer, parameter                              :: NORTH_   = 2
        integer, parameter                              :: EAST_    = 3
        integer, parameter                              :: SOUTH_   = 4
        integer, parameter                              :: WEST_    = 5
        real                                            :: Depth, DepthCurrent
        real                                            :: Bottom, BtomCurrent
        real                                            :: MinDischargeDepth
        integer                                         :: i,j,k
        
        !----------------------------------------------------------------------

    
        !Allocates temporary Matrix of possible discharge points
        allocate(PossibleDischargePoint(WorkSizeILB:WorkSizeIUB, WorkSizeJLB:WorkSizeJUB))
        PossibleDischargePoint = .true.
                    
        !Search the lowest point in the sourring points        
        NewLocationFound = .false.
        i = DischargeX%Localization%GridCoordinates%OldI
        j = DischargeX%Localization%GridCoordinates%OldJ
        k = DischargeX%Localization%GridCoordinates%K
                    
        DepthCurrent        = WaterColumnZ(i, j)
        MinDischargeDepth   = DischargeX%Localization%MinimumDischargeDepth


        PossibleDischargePoint(i, j) = .false.

        do while (.not. NewLocationFound)

            MoveDirection = NONE_                    

cd1 :       if (OpenPoints3D(i,j,k) == OpenPoint) then                                                     
                                
                !Point in the North
                if (i < WorkSizeIUB .and. OpenPoints3D(i+1,j,k) == OpenPoint) then                            
                    Depth = WaterColumnZ(i+1, j)                               
                    if ( MinDischargeDepth <= Depth .and. Depth < DepthCurrent .and. PossibleDischargePoint(i+1, j)) then
                        MoveDirection = NORTH_
                        DepthCurrent  = Depth
                    endif

                endif

                !Point in the East
                if (j < WorkSizeJUB .and. OpenPoints3D(i,j+1,k) == OpenPoint) then
                    Depth = WaterColumnZ(i, j+1)
                    if ( MinDischargeDepth <= Depth .and. Depth< DepthCurrent .and. PossibleDischargePoint(i, j+1)) then
                        MoveDirection = EAST_
                        DepthCurrent  = Depth
                    endif
                endif

               !Point in the South
                if (i > WorkSizeILB .and. OpenPoints3D(i-1,j,k) == OpenPoint) then
                    Depth = WaterColumnZ(i-1, j)
                    if ( MinDischargeDepth <= Depth .and. Depth < DepthCurrent .and. PossibleDischargePoint(i-1, j)) then
                        MoveDirection = SOUTH_
                        DepthCurrent  = Depth
                    endif
                endif

                !Point in the West
                if (j > WorkSizeJLB .and. OpenPoints3D(i,j-1,k) == OpenPoint) then
                    Depth = WaterColumnZ(i, j-1)
                    if ( MinDischargeDepth <= Depth .and. Depth < DepthCurrent .and. PossibleDischargePoint(i, j-1)) then
                        MoveDirection = WEST_
                        DepthCurrent  = Depth
                    endif
                endif
            
            else cd1
            !neste caso comparam-se as cotas (os openpoints n ficam todos bem com a minWC), 
            !o antigo ponto de descarga deixou de ser Open Point.

                !Point in the North
                if (i < WorkSizeIUB) then                            
                    Bottom = Bathymetry(i+1, j)                               
                    if (Bottom >= BtomCurrent .and. PossibleDischargePoint(i+1, j)) then
                        MoveDirection = NORTH_
                        BtomCurrent   = Bottom
                    endif
                endif

                !Point in the East
                if (j < WorkSizeJUB) then
                    Bottom = Bathymetry(i, j+1)
                    if (Bottom >= BtomCurrent .and. PossibleDischargePoint(i, j+1)) then
                        MoveDirection = EAST_
                        BtomCurrent   = Bottom
                    endif
                endif

                !Point in the South
                if (i > WorkSizeILB) then
                    Bottom = Bathymetry(i-1, j)
                    if (Bottom >= BtomCurrent .and. PossibleDischargePoint(i-1, j)) then
                        MoveDirection = SOUTH_
                        BtomCurrent   = Bottom
                    endif
                endif

                !Point in the West
                if (j > WorkSizeJLB) then
                    Bottom = Bathymetry (i, j-1)
                    if (Bottom >= BtomCurrent .and. PossibleDischargePoint(i, j-1)) then
                        MoveDirection = WEST_
                        BtomCurrent   = Bottom
                    endif
                endif
                
            endif cd1                                  
            
            select case (MoveDirection)

                case (NORTH_)
                    i = i + 1
                    PossibleDischargePoint(i, j) = .false.

                case (EAST_ )
                    j = j + 1
                    PossibleDischargePoint(i, j) = .false.

                case (SOUTH_)
                    i = i - 1
                    PossibleDischargePoint(i, j) = .false.

                case (WEST_ )
                    j = j - 1
                    PossibleDischargePoint(i, j) = .false.

                case (NONE_)
                    !a descarga n pode ir para sitio nenhum fica no mm sitio
                    NewLocationFound = .true.
                    newI            = i
                    newJ            = j

            end select

            !Checks out if the new defined point is an openpoint
            if (WaterColumnZ(i, j) > DischargeX%Localization%MinimumDischargeDepth) then
                NewLocationFound = .true.
                newI            = i
                newJ            = j                         
            endif

        enddo
        
        !Deallocates temporary matrix
        deallocate(PossibleDischargePoint)

    end subroutine FindDischargePointFromLast

    !--------------------------------------------------------------------------

    subroutine GetDischargesIDName (DischargesID, DischargeIDNumber, IDName, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: DischargesID
        integer,           intent(IN )                  :: DischargeIDNumber
        character(len=*),  intent(OUT)                  :: IDName
        integer, optional, intent(OUT)                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, STAT_
        type(T_IndividualDischarge), pointer            :: DischargeX


        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeIDNumber
                stop       'GetDischargesIDName - ModuleDischarges - ERR10'
            endif

            IDName = DischargeX%ID%name

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_


    end subroutine GetDischargesIDName

    !--------------------------------------------------------------------------

    subroutine GetDischargesNodeID (DischargesID, DischargeIDNumber, NodeID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: DischargesID
        integer,           intent(IN )                  :: DischargeIDNumber
        integer,           intent(OUT)                  :: NodeID
        integer, optional, intent(OUT)                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, STAT_
        type(T_IndividualDischarge), pointer            :: DischargeX


        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeIDNumber
                stop       'GetDischargesNodeID - ModuleDischarges - ERR01'
            endif

            !Code below does not makes sense for drainage network. 
            !If you want to impose discharge in Drainage Network and RunOff the code will not allow you to run.
            !Can be deleted in 2019 - Frank
            !if (DischargeX%Localization%Location2D) then
            !    write(*,*)'Discharge location not given as Node ID'
            !    write(*,*)trim(adjustl(DischargeX%ID%Name))
            !    stop 'GetDischargesNodeID - ModuleDischarges - ERR01a'
            !endif

            NodeID = DischargeX%Localization%NodeID

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_


    end subroutine GetDischargesNodeID

    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    subroutine GetDischargesReservoirID (DischargesID, DischargeIDNumber, ReservoirID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: DischargesID
        integer,           intent(IN )                  :: DischargeIDNumber
        integer,           intent(OUT)                  :: ReservoirID
        integer, optional, intent(OUT)                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, STAT_
        type(T_IndividualDischarge), pointer            :: DischargeX


        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeIDNumber
                stop       'GetDischargesReservoirID - ModuleDischarges - ERR01'
            endif

            if (DischargeX%Localization%Location2D) then
                write(*,*)'Discharge location not given as Reservoir ID'
                write(*,*)trim(adjustl(DischargeX%ID%Name))
                stop 'GetDischargesReservoirID - ModuleDischarges - ERR01a'
            endif

            ReservoirID = DischargeX%Localization%ReservoirID

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_


    end subroutine GetDischargesReservoirID

    !--------------------------------------------------------------------------    

   subroutine GetIsReservoirOutflow(DischargesID, DischargeIDNumber, IsOutflow, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: DischargesID
        integer,           intent(IN )                  :: DischargeIDNumber
        logical,           intent(OUT)                  :: IsOutflow
        integer, optional, intent(OUT)                  :: STAT

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge), pointer            :: DischargeX
        integer                                         :: STAT_CALL, STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeIDNumber, '.'
                stop       'Subroutine GetIsReservoirOutflow - ModuleDischarges. ERR01.'
            endif

            IsOutflow = DischargeX%IsReservoirOutflow
                    
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetIsReservoirOutflow


    !--------------------------------------------------------------------------

    Subroutine GetDischargeSpatialEmission(DischargesID, DischargeIDNumber, LineX, PolygonX, &
                                           SpatialEmission, XYZPointsX, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer,                        intent(IN ) :: DischargeIDNumber
        type (T_Lines),  pointer                    :: LineX
        type (T_Polygon), pointer                   :: PolygonX
        integer,                        intent(OUT) :: SpatialEmission
        type (T_XYZPoints), pointer, optional       :: XYZPointsX        
        integer, optional,              intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_
        integer                                     :: STAT_
        type(T_IndividualDischarge), pointer        :: DischargeX
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_)) then



            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)

cd3 :       if (STAT_CALL/=SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                stop       'Subroutine GetDischargeSpatialEmission - ModuleDischarges. ERR01.'
            end if cd3

            SpatialEmission = DischargeX%Localization%SpatialEmission

            PolygonX   => DischargeX%Localization%Polygon
            LineX      => DischargeX%Localization%Line
            if (present(XYZPointsX)) XYZPointsX => DischargeX%Localization%XYZPoints

            nullify(DischargeX)


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                              &
            STAT = STAT_

        !----------------------------------------------------------------------

    end Subroutine GetDischargeSpatialEmission

    !--------------------------------------------------------------------------

    Subroutine GetDischargeSpatialType(DischargesID, DischargeIDNumber, SpatialEmission, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer,                        intent(IN ) :: DischargeIDNumber
        integer,                        intent(OUT) :: SpatialEmission
        integer, optional,              intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_
        integer                                     :: STAT_
        type(T_IndividualDischarge), pointer        :: DischargeX
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_)) then



            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)

cd3 :       if (STAT_CALL/=SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                stop       'Subroutine GetDischargeSpatialEmission - ModuleDischarges. ERR01.'
            end if cd3

            SpatialEmission = DischargeX%Localization%SpatialEmission

            nullify(DischargeX)


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                              &
            STAT = STAT_

        !----------------------------------------------------------------------

    end Subroutine GetDischargeSpatialType

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    Subroutine GetDischargeFlowDistribuiton(DischargesID, DischargeIDNumber,            &
                                            nCells, FlowDistribution,                   &
                                            VectorI, VectorJ, VectorK, kmin, kmax, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer,                        intent(IN ) :: DischargeIDNumber
        integer,                        intent(OUT) :: nCells
        integer, optional,              intent(OUT) :: FlowDistribution
                        
        integer, dimension(:), pointer, optional    :: VectorI, VectorJ, VectorK
        integer,                        optional    :: kmin, kmax
        integer, optional,              intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_
        integer                                     :: STAT_
        type(T_IndividualDischarge), pointer        :: DischargeX
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_)) then



            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)

cd3 :       if (STAT_CALL/=SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                    stop       'Subroutine GetDischargeFlowDistribution - ModuleDischarges. ERR01.'
            end if cd3

            nCells = DischargeX%Localization%nCells

            if (present(FlowDistribution)) then

                FlowDistribution = DischargeX%Localization%FlowDistribution

            endif

            if (present(VectorI)) then

                call Read_Lock(mDischarges_, Me%InstanceID)
                VectorI => DischargeX%Localization%VectorI

            endif

            if (present(VectorJ)) then

                call Read_Lock(mDischarges_, Me%InstanceID)
                VectorJ => DischargeX%Localization%VectorJ
            
            endif

            if (present(VectorK)) then

                call Read_Lock(mDischarges_, Me%InstanceID)
                VectorK => DischargeX%Localization%VectorK

            endif
            
            if (present(kmin)) then

                kmin = DischargeX%Localization%kmin

            endif    
            
            if (present(kmax)) then

                kmax = DischargeX%Localization%kmax

            endif                       

            nullify(DischargeX)


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end Subroutine GetDischargeFlowDistribuiton
    
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    Subroutine GetDistributionCoefMass(DischargesID, DischargeIDNumber,                 &
                                       DistributionCoefMass, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer,                        intent(IN ) :: DischargeIDNumber
        real,    dimension(:), pointer              :: DistributionCoefMass
        integer, optional,              intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_
        integer                                     :: STAT_
        type(T_IndividualDischarge), pointer        :: DischargeX
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_)) then



            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)

cd3 :       if (STAT_CALL/=SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                    stop       'Subroutine GetDistributionCoefMass - ModuleDischarges. ERR01.'
            end if cd3


            call Read_Lock(mDischarges_, Me%InstanceID)
            DistributionCoefMass => DischargeX%Localization%DistributionCoefMass

            nullify(DischargeX)


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end Subroutine GetDistributionCoefMass

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine GetDischargeWaterFlow(DischargesID, TimeX, DischargeIDNumber,            &
                                     SurfaceElevation, Flow, SurfaceElevation2,         &
                                     FlowDistribution, FlowArea, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        type(T_Time),                   intent(IN ) :: TimeX
        integer,                        intent(IN ) :: DischargeIDNumber
        real   ,                        intent(IN)  :: SurfaceElevation
        real   ,                        intent(OUT) :: Flow
        real   , optional,              intent(IN)  :: SurfaceElevation2
        real   , optional,              intent(IN)  :: FlowDistribution        
        real   , optional,              intent(OUT) :: FlowArea
        integer, optional,              intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_
        integer                                     :: STAT_
        type(T_IndividualDischarge), pointer        :: DischargeX
        real                                        :: H, C, A
        integer                                     :: STAT_CALL
        logical                                     :: AssociateIntakeFlowON
        real                                        :: FlowFraction
        real                                        :: UpstreamH, DownstreamH, FlowDir
        real                                        :: D, TopValveH, BottomValveH, Haux, Theta
        real                                        :: PipeFriction, P, Rh
        real                                        :: dQ1, dQ2
        integer                                     :: iAux
         !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
cd3 :       if (STAT_CALL/=SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                stop       'SetDischargeWaterFlow - ModuleDischarges - ERR01'
            end if cd3
            
            if (present(FlowArea)) FlowArea = -99. 


            if(DischargeX%FromIntake%ON .and. DischargeX%FromIntake%AssociateFlow)then

                AssociateIntakeFlowON    = ON
                FlowFraction                = DischargeX%FromIntake%FlowFraction

                !DischargeX becomes the intake (flow is multiplied by -1.0 after it is determined
                call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeX%FromIntake%IntakeID)
cd31:           if (STAT_CALL/=SUCCESS_) then 
                    write(*,*) 'Can not find intake discharge number ', DischargeIDNumber
                    stop       'SetDischargeWaterFlow - ModuleDischarges - ERR02'
                end if cd31

            else

                AssociateIntakeFlowON = OFF

            end if


cd2:        if (DischargeX%DischargeType == Normal .and. DischargeX%WaterFlow%Variable) then


                Flow = TimeSerieValue(DischargeX%TimeSerie, DischargeX%UseOriginalValues, &
                                     TimeX, DischargeX%WaterFlow%FlowColumn)

            elseif (DischargeX%DischargeType == FlowOver) then

                !Q = cv * b * sqrt(2*g) * H^(1.5)            
                if (DischargeX%ByPass%ON .and. DischargeX%ByPass%Side == SideB) then
                    !Downstream of the Weir
                    H = SurfaceElevation2 - DischargeX%FlowOver%CrestHeigth

                else
                    !Upstream of the Weir
                    H = SurfaceElevation  - DischargeX%FlowOver%CrestHeigth
                endif            
            
                if (H > 0) then
                    Flow = -sqrt(19.6) * DischargeX%FlowOver%DischargeCoeficient *       &
                                         DischargeX%FlowOver%WeirLength  * H ** 1.5
                else
                    Flow = 0.
                endif
                
                
                if (DischargeX%ByPass%ON .and. DischargeX%ByPass%Side == SideB) then
                    Flow = - Flow
                endif
                
            elseif (DischargeX%DischargeType == RatingCurve) then

                
                if (SurfaceElevation < DischargeX%RatingCurve%Level(1)) then
                    Flow = 0.0
                else if (SurfaceElevation > DischargeX%RatingCurve%Level(DischargeX%RatingCurve%nValues)) then
                    Flow = -1.0 * DischargeX%RatingCurve%Flow(DischargeX%RatingCurve%nValues)
                else

                    iAux = 2
                    do while (SurfaceElevation >= DischargeX%RatingCurve%Level(iAux))
                        iAux = iAux + 1
                    enddo
                        
                        
                    dQ1 = SurfaceElevation - DischargeX%RatingCurve%Level(iAux-1)
                    dQ2 = DischargeX%RatingCurve%Level(iAux) - SurfaceElevation
                    
                    Flow = -1.0 * (dQ1 * DischargeX%RatingCurve%Flow(iAux)             +           &
                            dQ2 * DischargeX%RatingCurve%Flow(iAux-1) )         /           &
                           (DischargeX%RatingCurve%Level(iAux)                  -           &
                            DischargeX%RatingCurve%Level(iAux-1)) 
                    
                endif
                    
                
                

                                    
            elseif (DischargeX%DischargeType == Valve) then


                !Q = C * A * sqrt(2*g) * sqrt(H)
                if (SurfaceElevation > SurfaceElevation2) then
                    UpstreamH    = SurfaceElevation
                    DownstreamH  = SurfaceElevation2
                    FlowDir      = -1
                else
                    UpstreamH    = SurfaceElevation2
                    DownstreamH  = SurfaceElevation
                    FlowDir      = 1
                endif
                
                
                if      (DischargeX%Valve%SectionType == circular_area) then
                
                    D  = DischargeX%Valve%Diameter                    
                    
                elseif  (DischargeX%Valve%SectionType == rectangular_area) then
                
                    D = DischargeX%Valve%Height
                
                endif                 

                TopValveH    =   DischargeX%Valve%SillHeigth + D
                BottomValveH =   DischargeX%Valve%SillHeigth
                
                H            = UpstreamH - max(DownstreamH, BottomValveH)
                
                C            = DischargeX%Valve%DischargeCoeficient                
                

                !if the axis valve minus the radius is above the water level in both sides than there is no flow
                if ( UpstreamH <= BottomValveH) then
                
                    H    =  0.
                    
                    Flow =  0.

                else
                    !pressure conditions    
                    if (UpstreamH >= TopValveH) then
                    
                        if      (DischargeX%Valve%SectionType == circular_area) then
                        
                            A    = Pi * (D/2.)**2.
                            P    = Pi * D
                            
                        elseif (DischargeX%Valve%SectionType == rectangular_area) then
                        
                            A    = D * DischargeX%Valve%Width
                            P    = 2. * (D + DischargeX%Valve%Width)

                        endif                                     

                    !free surface
                    else
                    
                        Haux  = UpstreamH - BottomValveH
                        
                        if      (DischargeX%Valve%SectionType == circular_area) then
                        
                            Theta = 2.* acos(1.-2.*Haux/D)
                            A     = (Theta - sin(Theta))* D**2. / 8.
                            P     = Theta * D
                            
                        elseif  (DischargeX%Valve%SectionType == rectangular_area) then

                            A    = Haux * DischargeX%Valve%Width
                            P    = 2 * Haux + DischargeX%Valve%Width
                        
                        endif                         
                        
                        
                        if (Haux < 0.01) then
                            !no flow in the pipe is assumed
                            C = 0.
                        endif                    
                        
                    
                    endif
                    
                    if (P > 0.) then
                        Rh  = A / P
                    else
                        Rh  = 0. 
                    endif                        
                    
                    
                    if (DischargeX%Valve%PipeLength > 0.) then
                        if (Rh > 0.) then
                        !based in BASIC HYDRAULIC PRINCIPLES OF OPEN-CHANNEL FLOW By Harvey E. Jobson and David C. Froehlich
                        !U.S. GEOLOGICAL SURVEY Open-File Report 88-707

                        ![ ]     =  [m/s2] * [s/m^0.333]^2 * [m] / [m]^1.333] = [m] / [m^0.6667] * [m] / [m]^1.333 =
                        !          [m]^(1 -0.6667 + 1 - 1.333] = [m]^0= [ ]
                            PipeFriction = 2. * Gravity * C**2. * DischargeX%Valve%PipeManning ** 2. * &
                                           DischargeX%Valve%PipeLength / Rh**1.3333
                            H = H / (1. + PipeFriction)
                        endif                            
                    endif    

                    Flow = sqrt(2.* Gravity) * C * A * sqrt(H)
                    
                    if (present(FlowArea)) FlowArea = A

                    Flow = Flow * FlowDir
                
                    if     (DischargeX%ByPass%OneWay) then

                        if ((DischargeX%ByPass%Side == SideA .and. FlowDir > 0.) .or.         &
                            (DischargeX%ByPass%Side == SideB .and. FlowDir < 0.)) Flow = 0.
                    endif
                endif

            else

                Flow = DischargeX%WaterFlow%Scalar

            end if cd2

            if(AssociateIntakeFlowON)then

                if(Flow <= 0)then 
                    !If the discharge comes from an intake then the flow
                    !from the intake must be negative, so that the flow
                    !of the discharge is positive
                    Flow = Flow * (-1.0) * FlowFraction
                else 
                    write(*,*)"Discharge has flow based on intake"
                    write(*,*)"However, intake flow is positive and should be negative"
                    stop 'GetDischargeWaterFlow - ModuleDischarges - ERR03'
                end if 
               
            end if
            
            
            if (present(FlowDistribution)) then
                Flow = Flow * FlowDistribution
            endif


            nullify(DischargeX)


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDischargeWaterFlow

    !--------------------------------------------------------------------------
    
    real function TimeSerieValue(TimeSerieID, UseOriginalValues, TimeX, XColumn)
    
        !Arguments-------------------------------------------------------------
        integer,                        intent(IN ) :: TimeSerieID
        logical,                        intent(IN ) :: UseOriginalValues
        type(T_Time),                   intent(IN ) :: TimeX
        integer, optional,              intent(IN ) :: XColumn

        !Local-----------------------------------------------------------------
        type(T_Time)                                :: Time1, Time2
        real                                        :: Value1, Value2, NewValue
        logical                                     :: TimeCycle
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------        
    
        call GetTimeSerieValue(TimeSerieID, TimeX, XColumn, Time1, Value1,              &
                               Time2, Value2, TimeCycle, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'TimeSerieValue - ModuleDischarges - ERR20'

        if (TimeCycle) then
            NewValue = Value1
        else
            
            if(UseOriginalValues)then
            
                NewValue = Value1
                
            else
        
                !Interpolates Value for current instant
                call InterpolateValueInTime(TimeX, Time1, Value1, Time2, Value2,        &
                                            NewValue)
            end if
            
        endif
        
        TimeSerieValue = NewValue
                
    end function TimeSerieValue

    !--------------------------------------------------------------------------

    subroutine SetDischargeWaterFlow(DischargesID, DischargeIDNumber,                   &
                                     Flow, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer,                        intent(IN ) :: DischargeIDNumber
        real   ,                        intent(IN)  :: Flow
        integer, optional,              intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_
        integer                                     :: STAT_
        type(T_IndividualDischarge), pointer        :: DischargeX
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
cd3 :       if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ',DischargeIDNumber
                stop       'SetDischargeWaterFlow - ModuleDischarges - ERR01'
            end if cd3


            if (DischargeX%DischargeType == OpenMILink) then

                DischargeX%WaterFlow%Scalar = Flow

            else
            
                write(*,*)  'SetDischargeWaterFlow can only be called if discharge type is OpenMILink (4)'
                stop        'SetDischargeWaterFlow - ModuleDischarges - ERR02'
            
            endif


            nullify(DischargeX)


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SetDischargeWaterFlow

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine SetDischargeInterceptionRatio(DischargesID, DischargeIDNumber,           &
                                             InterceptionRatio, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer,                        intent(IN ) :: DischargeIDNumber
        real   ,                        intent(IN)  :: InterceptionRatio
        integer, optional,              intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_
        integer                                     :: STAT_
        type(T_IndividualDischarge), pointer        :: DischargeX
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
cd3 :       if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ',DischargeIDNumber
                stop       'SetDischargeInterceptionRatio - ModuleDischarges - ERR01'
            end if cd3

            DischargeX%Localization%InterceptionRatio = InterceptionRatio

            nullify(DischargeX)


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SetDischargeInterceptionRatio

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetDischargeInterceptionRatio(DischargesID, DischargeIDNumber,           &
                                             InterceptionRatio, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer,                        intent(IN ) :: DischargeIDNumber
        real   ,                        intent(OUT) :: InterceptionRatio
        integer, optional,              intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_
        integer                                     :: STAT_
        type(T_IndividualDischarge), pointer        :: DischargeX
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
cd3 :       if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ',DischargeIDNumber
                stop       'GetDischargeInterceptionRatio - ModuleDischarges - ERR01'
            end if cd3

            InterceptionRatio = DischargeX%Localization%InterceptionRatio

            nullify(DischargeX)


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDischargeInterceptionRatio

    !--------------------------------------------------------------------------
        
    

    subroutine GetDischargeFlowVelocity(DischargesID, TimeX, DischargeIDNumber,          &
                                        VelocityU, VelocityV, VelocityW, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        real,    optional, intent(OUT)              :: VelocityU, VelocityV, VelocityW
        integer, optional, intent(OUT)              :: STAT
        integer,           intent(IN )              :: DischargeIDNumber
        type(T_Time),      intent(IN )              :: TimeX

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_
        type(T_IndividualDischarge), pointer        :: DischargeX
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (.not. present(VelocityU) .and. .not. present(VelocityV) .and. .not. present(VelocityW))   &
                stop 'Subroutine GetDischargeFlowVelocity - ModuleDischarges. ERR01.'

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber = DischargeIDNumber)
                

cd3 :       if (STAT_CALL /= SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                    stop 'Subroutine GetDischargeFlowVelocity - ModuleDischarges. ERR02.'
            end if cd3


cd4 :       if (Present(VelocityU)) then 

cd2:            if (DischargeX%VelocityFlow%UVariable) then

                    VelocityU = TimeSerieValue(DischargeX%TimeSerie,                    &
                                               DischargeX%UseOriginalValues, TimeX,     &
                                               DischargeX%VelocityFlow%UColumn)

                else

                    VelocityU = DischargeX%VelocityFlow%UScalar

                end if cd2

            endif cd4

cd5 :       if (Present(VelocityV)) then 

cd6:            if (DischargeX%VelocityFlow%VVariable) then

                    VelocityV = TimeSerieValue(DischargeX%TimeSerie,                    &
                                               DischargeX%UseOriginalValues, TimeX,     &
                                               DischargeX%VelocityFlow%VColumn)
                else
                    VelocityV = DischargeX%VelocityFlow%VScalar

                end if cd6

            endif cd5


cd7 :       if (Present(VelocityW)) then 

cd8:            if (DischargeX%VelocityFlow%WVariable) then

                    VelocityW = TimeSerieValue(DischargeX%TimeSerie,                    &
                                               DischargeX%UseOriginalValues, TimeX,     &
                                               DischargeX%VelocityFlow%WColumn)
                else
                    VelocityW = DischargeX%VelocityFlow%WScalar

                end if cd8

            endif cd7

            nullify(DischargeX)

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDischargeFlowVelocity

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine GetDischargeON(DischargesID, DischargeIDNumber, IgnoreON, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer,           intent(IN )              :: DischargeIDNumber
        logical,           intent(OUT)              :: IgnoreON
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_
        type(T_IndividualDischarge), pointer        :: DischargeX
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber = DischargeIDNumber)
                

cd3 :       if (STAT_CALL /= SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                    stop 'Subroutine GetDischargeON - ModuleDischarges. ERR10.'
            end if cd3

            IgnoreON = DischargeX%IgnoreON

            nullify(DischargeX)

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDischargeON
    !--------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Checks if a discharge is of the type "upscalling"
    !>@param[in] DischargesID, number, STAT_CALL    
    logical function IsUpscaling(DischargesID, number)
    
    !Arguments-----------------------------------------------------------------
    integer                                     :: DischargesID, number
    !local---------------------------------------------------------------------
    integer                                     :: ready_, STAT_CALL
    type(T_IndividualDischarge), pointer        :: DischargeX
    !--------------------------------------------------------------------------
    
    call Ready(DischargesID, ready_)    
        
    if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

        call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=number)
        if (STAT_CALL/=SUCCESS_) then 
            write(*,*) 'Can not find discharge number ', number, '.'
            stop       'Function GetIsUpscaling - ModuleDischarges. ERR01.'
        endif

        
        IsUpscaling = DischargeX%WaterFlow%Upscaling
    else 
        IsUpscaling = .false.
    end if
           
    return
    
    end function IsUpscaling
    
    !--------------------------------------------------------------------------

    subroutine GetDischargeConcentration(DischargesID, TimeX,DischargeIDNumber,         &
                                         Concentration, PropertyIDNumber,               &
                                         PropertyFromIntake, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        type(T_Time),      intent(IN)               :: TimeX
        integer,           intent(IN)               :: DischargeIDNumber
        real,              intent(OUT)              :: Concentration
        integer,           intent(IN)               :: PropertyIDNumber
        logical, optional, intent(OUT)              :: PropertyFromIntake
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: DischargeX
        type(T_Property),            pointer        :: PropertyX
        integer                                     :: ready_
        integer                                     :: STAT_
        integer                                     :: STAT_CALL
        real                                        :: PropertyIncreaseValue
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                stop  'GetDischargeConcentration - ModuleDischarges - ERR01'
            endif
             
            call Search_Property(DischargeX, PropertyX, STAT_CALL, PropertyXIDNumber=PropertyIDNumber)
            if (STAT_CALL/=SUCCESS_) then 
                !If the proeprty is not found the program don't stop is return a error 
                !not found
                if (STAT_CALL /= NOT_FOUND_ERR_) then 
                    stop  'GetDischargeConcentration - ModuleDischarges - ERR02'
                endif
            endif
            
cd2 :       if (STAT_CALL == SUCCESS_) then

                if(PropertyX%FromIntake)then
                    
                    if (present(PropertyFromIntake)) PropertyFromIntake      = ON
                    PropertyIncreaseValue   = PropertyX%IncreaseValue


                    !DischargeX becomes the      discharge
                    call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeX%FromIntake%IntakeID)
                        if (STAT_CALL/=SUCCESS_) then 
                        write(*,*) 'Can not find intake discharge number ', DischargeIDNumber
                        stop       'GetDischargeConcentration - ModuleDischarges - ERR03'
                    end if

                    call Search_Property(DischargeX, PropertyX, STAT_CALL, PropertyXIDNumber=PropertyIDNumber)
                    if (STAT_CALL/=SUCCESS_) then 
                        !If the proeprty is not found the program don't stop is return a error 
                        !not found
                        if (STAT_CALL /= NOT_FOUND_ERR_) then 
                            stop  'GetDischargeConcentration - ModuleDischarges - ERR04'
                        endif
                    endif
                else
                    
                    if (present(PropertyFromIntake)) PropertyFromIntake = OFF

                end if
                                                                         
                if (PropertyX%Variable) then

                    Concentration = TimeSerieValue(PropertyX%TimeSerie,                 &
                                                   DischargeX%UseOriginalValues, TimeX, &
                                                   PropertyX%ConcColumn)
                else
 
                    Concentration = PropertyX%Scalar

                endif
                
                if (present(PropertyFromIntake)) then
                    if(PropertyFromIntake) then
                        Concentration = Concentration + PropertyIncreaseValue
                    end if
                endif

                nullify(PropertyX)

                nullify(DischargeX)

                STAT_ = SUCCESS_

            else if (STAT_CALL == NOT_FOUND_ERR_) then cd2


                nullify(PropertyX)

                nullify(DischargeX)

                STAT_ = NOT_FOUND_ERR_


            end if cd2
            

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end Subroutine GetDischargeConcentration

    !--------------------------------------------------------------------------
    
    subroutine GetIntakePosition(DischargesID, DischargeIDNumber, IntakeI, IntakeJ, IntakeK, STAT)
                                      

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer,           intent(IN)               :: DischargeIDNumber
        integer,           intent(OUT)              :: IntakeI, IntakeJ, IntakeK
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: DischargeX
        type(T_IndividualDischarge), pointer        :: IntakeX
        integer                                     :: ready_
        integer                                     :: STAT_
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                stop  'GetIntakePosition - ModuleDischarges - ERR01'
            endif

            call Search_Discharge(IntakeX, STAT_CALL, DischargeXIDNumber=DischargeX%FromIntake%IntakeID)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                stop  'GetIntakeConcentration - ModuleDischarges - ERR03'
            endif

            IntakeI = IntakeX%Localization%GridCoordinates%I
            IntakeJ = IntakeX%Localization%GridCoordinates%J
            IntakeK = IntakeX%Localization%GridCoordinates%K
                

            nullify(IntakeX)

            nullify(DischargeX)

            STAT_ = SUCCESS_

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end Subroutine GetIntakePosition

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetByPassConcIncrease(DischargesID, DischargeIDNumber,                   &
                                     PropertyIDNumber, ConcIncrease, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer,           intent(IN)               :: DischargeIDNumber
        integer,           intent(IN)               :: PropertyIDNumber
        real,              intent(OUT)              :: ConcIncrease        
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: DischargeX
        type(T_Property),            pointer        :: PropertyX
        integer                                     :: ready_
        integer                                     :: STAT_
        integer                                     :: STAT_CALL
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                stop  'GetByPassConcIncrease - ModuleDischarges - ERR01'
            endif
             
            call Search_Property(DischargeX, PropertyX, STAT_CALL, PropertyXIDNumber=PropertyIDNumber)
            if (STAT_CALL/=SUCCESS_) then 
                !If the proeprty is not found the program don't stop is return a error 
                !not found
                if (STAT_CALL /= NOT_FOUND_ERR_) then 
                    stop  'GetByPassConcIncrease - ModuleDischarges - ERR02'
                endif
            endif
            
cd2 :       if (STAT_CALL == SUCCESS_) then

                if (DischargeX%ByPass%ON) then
                    ConcIncrease   = PropertyX%IncreaseValue
                else
                    ConcIncrease   = 0.
                endif                    

                nullify(PropertyX)

                nullify(DischargeX)

                STAT_ = SUCCESS_

            else if (STAT_CALL == NOT_FOUND_ERR_) then cd2


                nullify(PropertyX)

                nullify(DischargeX)

                STAT_ = NOT_FOUND_ERR_


            end if cd2
            

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end Subroutine GetByPassConcIncrease

    !--------------------------------------------------------------------------




    !--------------------------------------------------------------------------

    subroutine Search_Discharge(DischargeX, STAT_, DischargeXIDNumber)

        !Arguments--------------------------------------------------------------
        type(T_IndividualDischarge),  pointer         :: DischargeX
        integer,                      intent (OUT)    :: STAT_
        integer,                      intent (IN)     :: DischargeXIDNumber

        !Begin-----------------------------------------------------------------

        STAT_  = NOT_FOUND_ERR_
        
        !Tries a "quick find". Improves performance if the model has many discharges (e.g. RiverNetwork - SWAT)
        DischargeX => Me%CurrentDischarge
        do while (associated(DischargeX)) 
            if (DischargeX%ID%IDNumber==DischargeXIDNumber) then
                STAT_ = SUCCESS_
                Me%CurrentDischarge => DischargeX
                return
            else
                DischargeX => DischargeX%Next
            endif
        enddo
        
        !Nothing found yet, start from beginning
        DischargeX => Me%FirstDischarge

        do while (associated(DischargeX)) 
            if (DischargeX%ID%IDNumber==DischargeXIDNumber) then
                STAT_ = SUCCESS_
                Me%CurrentDischarge => DischargeX
                exit
            else
                DischargeX => DischargeX%Next
            endif
        enddo

    end subroutine Search_Discharge

    !--------------------------------------------------------------------------

    subroutine Search_Discharge_ByName(DischargeX, STAT_, DischargeXIDName)

        !Arguments--------------------------------------------------------------
        type(T_IndividualDischarge),  pointer         :: DischargeX
        integer,                      intent (OUT)    :: STAT_
        character(len=*),             intent (IN)     :: DischargeXIDName

        !Begin-----------------------------------------------------------------

        STAT_  = NOT_FOUND_ERR_

        DischargeX => Me%FirstDischarge

        do while (associated(DischargeX)) 
            if (DischargeX%ID%Name == DischargeXIDName) then
                STAT_ = SUCCESS_
                Me%CurrentDischarge => DischargeX
                exit
            else
                DischargeX => DischargeX%Next
            endif
        enddo



    end subroutine Search_Discharge_ByName

    !--------------------------------------------------------------------------

    subroutine Search_Property(DischargeX, PropertyX, STAT_, PropertyXIDNumber)

        !External--------------------------------------------------------------

        type(T_IndividualDischarge),pointer           :: DischargeX
        type(T_Property),           pointer           :: PropertyX
        integer         ,           intent (OUT)      :: STAT_
        integer         ,           intent (IN)       :: PropertyXIDNumber

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        STAT_  = UNKNOWN_

        !Tries a "quick find". Improves performance if the model has many discharges / properties
        !(e.g. RiverNetwork - SWAT)
        PropertyX => DischargeX%CurrProperty
        do while (associated(PropertyX)) 
            if (PropertyX%ID%IDNumber==PropertyXIDNumber) then
                STAT_ = SUCCESS_
                DischargeX%CurrProperty => PropertyX
                return
            else
                PropertyX => PropertyX%Next
            endif
        enddo
        
        
        !Nothing found yet, start from beginning
        PropertyX => DischargeX%FirstProperty

        do while (associated(PropertyX)) 
           if (PropertyX%ID%IDNumber==PropertyXIDNumber) then
                STAT_ = SUCCESS_
                DischargeX%CurrProperty => PropertyX
                return
           else
               PropertyX => PropertyX%Next                 
           endif
        enddo

        STAT_  = NOT_FOUND_ERR_ ! The PropertyX wasn't found
                

        !----------------------------------------------------------------------
      
    end subroutine Search_Property

    !--------------------------------------------------------------------------
  
    subroutine GetDischargeParameters (DirectionX, DirectionY)

        !Arguments-------------------------------------------------------------

        integer, optional, intent(OUT) :: DirectionX      
        integer, optional, intent(OUT) :: DirectionY       
    
        !----------------------------------------------------------------------

        if (present(DirectionX           )) DirectionX        = DirectionX_
        if (present(DirectionY           )) DirectionY        = DirectionY_
                                                         
        !----------------------------------------------------------------------

    end subroutine GetDischargeParameters

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine SetLocationCellsZ (DischargeID, DischargeIDNumber, nCells, VectorI, VectorJ, VectorK, STAT)

        !Arguments--------------------------------------------------------------
        integer,                      intent (IN)     :: DischargeID, DischargeIDNumber
        integer,                      intent (IN)     :: nCells                         
        integer, dimension(:),        pointer         :: VectorI, VectorJ, VectorK
        integer, optional,            intent (OUT)    :: STAT

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge),pointer           :: DischargeX
        integer                                       :: ready_              
        integer                                       :: STAT_, STAT_CALL           

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargeID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)

            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                stop  'Sub. SetLocationCellsZ - ModuleDischarges - ERR01'
            endif


            DischargeX%Localization%VectorI => VectorI
            DischargeX%Localization%VectorJ => VectorJ
            DischargeX%Localization%VectorK => VectorK

            DischargeX%Localization%nCells  =  nCells

            nullify(DischargeX)
   
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine SetLocationCellsZ

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine SetLayer (DischargeID, DischargeIDNumber, K, STAT)

        !Arguments--------------------------------------------------------------
        integer,                      intent (IN)     :: DischargeID, DischargeIDNumber
        integer,                      intent (IN)     :: K
        integer, optional,            intent (OUT)    :: STAT

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge),pointer           :: DischargeX
        integer                                       :: ready_              
        integer                                       :: STAT_, STAT_CALL           

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargeID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)

            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                stop  'Sub. SetLayer - ModuleDischarges - ERR10'
            endif

            DischargeX%Localization%GridCoordinates%K = K

            nullify(DischargeX)
   
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine SetLayer

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine SetDistributionCoefMass (DischargeID, DischargeIDNumber, DistributionCoefMass, STAT)

        !Arguments--------------------------------------------------------------
        integer,                      intent (IN)     :: DischargeID, DischargeIDNumber                       
        real,    dimension(:),        pointer         :: DistributionCoefMass
        integer, optional,            intent (OUT)    :: STAT

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge),pointer           :: DischargeX
        integer                                       :: ready_              
        integer                                       :: STAT_, STAT_CALL           

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargeID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)

            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                stop  'Sub. SetDistributionCoefMass - ModuleDischarges - ERR01'
            endif

            DischargeX%Localization%DistributionCoefMass => DistributionCoefMass
 
            nullify(DischargeX)
   
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine SetDistributionCoefMass

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------    

    subroutine UngetDischarges1Dinteger(DischargesID, Array, STAT)

        !Arguments-------------------------------------------------------------

        integer,               intent(IN ) :: DischargesID
        integer, pointer,   dimension(:  ) :: Array
        integer, optional,     intent(OUT) :: STAT
   
        

        !External--------------------------------------------------------------

        integer :: ready_   

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mDischarges_, Me%InstanceID, "UngetDischarges1Dinteger")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetDischarges1Dinteger

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    
    subroutine UngetDischarges1Dreal(DischargesID, Array, STAT)

        !Arguments-------------------------------------------------------------

        integer,               intent(IN ) :: DischargesID
        real,    pointer,   dimension(:  ) :: Array
        integer, optional,     intent(OUT) :: STAT
   
        

        !External--------------------------------------------------------------

        integer :: ready_   

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mDischarges_, Me%InstanceID, "UngetDischarges1Dreal")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetDischarges1Dreal

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine CorrectsCellsDischarges(DischargesID, NDischarge, i, j, STAT)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )              :: DischargesID
        integer,           intent(IN )              :: NDischarge, i, j
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge),  pointer       :: DischargeX
        integer                                     :: ready_, STAT_, STAT2_ 
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            call Search_Discharge(DischargeX, STAT2_, NDischarge)

                if (STAT2_ == SUCCESS_) then

                    DischargeX%Localization%GridCoordinates%I = i
                    DischargeX%Localization%GridCoordinates%J = j
          
                    STAT_ = SUCCESS_

                endif
        else              
         
            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine CorrectsCellsDischarges

    !--------------------------------------------------------------------------
    subroutine CorrectsByPassCellsDischarges(DischargesID, NDischarge, i, j, STAT)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )              :: DischargesID
        integer,           intent(IN )              :: NDischarge, i, j
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge),  pointer       :: DischargeX
        integer                                     :: ready_, STAT_, STAT2_ 
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            call Search_Discharge(DischargeX, STAT2_, NDischarge)

                if (STAT2_ == SUCCESS_) then

                    DischargeX%ByPass%I = i
                    DischargeX%ByPass%J = j
          
                    STAT_ = SUCCESS_

                endif
        else              
         
            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine CorrectsByPassCellsDischarges

    !--------------------------------------------------------------------------


    subroutine TryIgnoreDischarge(DischargesID, NDischarge, IgnoreOK, STAT)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )              :: DischargesID
        integer,           intent(IN )              :: NDischarge
        logical,           intent(OUT)              :: IgnoreOK
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge),  pointer       :: DischargeX
        integer                                     :: ready_, STAT_, STAT2_ 
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            call Search_Discharge(DischargeX, STAT2_, NDischarge)

                if (STAT2_ == SUCCESS_) then

                    if (Me%IgnoreON) then
                        DischargeX%IgnoreON = .true. 
                        IgnoreOK            = .true.
                    else
                        IgnoreOK            = .false.
                    endif
                                  
                    STAT_ = SUCCESS_

                endif
        else              
         
            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine TryIgnoreDischarge
    
!----------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Kill_Discharges(DischargesID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_            
        type(T_IndividualDischarge), pointer        :: DischargeX, DischargeToKill
        integer                                     :: nUsers

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mDISCHARGES_,  Me%InstanceID)
  
            if (nUsers == 0) then
                
                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'Kill_Discharges - ModuleDischarges - ERR01'

                ! Deallocates all the Discharges
                DischargeX=> Me%FirstDischarge

                do while(associated(DischargeX))  

                    DischargeToKill => DischargeX
                    DischargeX      => DischargeX%Next
                    
                    if (DischargeToKill%Localization%TrackLocation)  &
                        call UnitsManager (DischargeToKill%Localization%TrackLocationFileUnitNumber, CLOSE_FILE)
                        
                    call KillIndividualDischarge(DischargeToKill)                                        

                end do
                nullify   (Me%FirstDischarge,Me%LastDischarge)

                call DeallocateInstance

                DischargesID = 0
            
                STAT_ = SUCCESS_
            
            end if

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine Kill_Discharges

    !--------------------------------------------------------------------------

    subroutine KillIndividualDischarge(DischargeToDelete)

        !Arguments-------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: DischargeToDelete

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: CurrentDischarge
        type(T_IndividualDischarge), pointer        :: Prev
        type(T_IndividualDischarge), pointer        :: Next
        type(T_Property           ), pointer        :: PropertyX
        integer                                     :: STAT_CALL


        CurrentDischarge => Me%FirstDischarge

        do 
            if (CurrentDischarge%ID%IDnumber == DischargeToDelete%ID%IDnumber) then

                Prev => CurrentDischarge%Prev
                Next => CurrentDischarge%Next

                !Updates foward pointer
                if (associated(CurrentDischarge%Prev)) then
                    Prev%Next  => Next
                endif

                !Updates backward pointer
                if (associated(CurrentDischarge%Next)) then
                    Next%Prev      => Prev
                endif

                !Updates first domain
                if (DischargeToDelete%ID%IDnumber == Me%FirstDischarge%ID%IDnumber) then
                    Me%FirstDischarge => Next
                endif

                !Updates last domain
                if (DischargeToDelete%ID%IDnumber == Me%LastDischarge%ID%IDnumber) then
                    Me%LastDischarge => Prev
                endif


                if (DischargeToDelete%TimeSerie /= 0) then

                    call KillTimeSerie(DischargeToDelete%TimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Subroutine Kill_Discharges; ModuleDischarges. ERR10.' 

                end if


                PropertyX => DischargeToDelete%FirstProperty


                do while(associated(PropertyX))  

                    if (PropertyX%PropTimeSerie) then
                        call KillTimeSerie(PropertyX%TimeSerie, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Subroutine Kill_Discharges; ModuleDischarges. ERR20.' 
                    endif

                    if (associated(PropertyX%Next)) then
                        PropertyX => PropertyX%Next
                        deallocate(PropertyX%Prev, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'Subroutine Kill_Discharges; ModuleDischarges. ERR30.'

                        nullify    (PropertyX%Prev)
                    else 
                        deallocate(PropertyX, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'Subroutine Kill_Discharges; ModuleDischarges. ERR40.' 
                                           
                        nullify    (PropertyX)
                    end if
                end do

                if (associated(DischargeToDelete%Localization%VectorI)) then

                    deallocate(DischargeToDelete%Localization%VectorI)
                    nullify   (DischargeToDelete%Localization%VectorI)

                endif

                if (associated(DischargeToDelete%Localization%VectorJ)) then

                    deallocate(DischargeToDelete%Localization%VectorJ)
                    nullify   (DischargeToDelete%Localization%VectorJ)

                endif


                if (associated(DischargeToDelete%Localization%VectorK)) then

                    deallocate(DischargeToDelete%Localization%VectorK)
                    nullify   (DischargeToDelete%Localization%VectorK)

                endif
                

                if (associated(DischargeToDelete%Localization%DistributionCoefMass)) then

                    deallocate(DischargeToDelete%Localization%DistributionCoefMass)
                    nullify   (DischargeToDelete%Localization%DistributionCoefMass)

                endif                

                nullify(DischargeToDelete%FirstProperty, DischargeToDelete%LastProperty)


                nullify   (DischargeToDelete%Next)
                nullify   (DischargeToDelete%Prev)
                deallocate(DischargeToDelete)
                nullify   (DischargeToDelete)

                return

            endif
        
            CurrentDischarge => CurrentDischarge%Next

        enddo

    end subroutine KillIndividualDischarge

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Discharges), pointer                :: AuxDischarges
        type (T_Discharges), pointer                :: PreviousDischarges

        !Updates pointers
        if (Me%InstanceID == FirstDischarges%InstanceID) then
            FirstDischarges => FirstDischarges%Next
        else
            PreviousDischarges => FirstDischarges
            AuxDischarges      => FirstDischarges%Next
            do while (AuxDischarges%InstanceID /= Me%InstanceID)
                PreviousDischarges => AuxDischarges
                AuxDischarges      => AuxDischarges%Next
            enddo

            !Now update linked list
            PreviousDischarges%Next => AuxDischarges%Next

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

    subroutine Ready (DischargesID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (DischargesID > 0) then
            call LocateObjDischarges (DischargesID)
            ready_ = VerifyReadLock (mDISCHARGES_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjDischarges (DischargesID)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID

        !Local-----------------------------------------------------------------

        Me => FirstDischarges
        do while (associated (Me))
            if (Me%InstanceID == DischargesID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                                          &
            stop 'ModuleDischarges - LocateObjDischarges - ERR01'

    end subroutine LocateObjDischarges

    !--------------------------------------------------------------------------

#ifdef _OPENMI_

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetNumberOfDischarges
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETNUMBEROFDISCHARGES"::GetNumberOfDischarges
    !DEC$ ENDIF
    !Return the number of Error Messages
    integer function GetNumberOfDischarges(DischargesID)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         

        if (.not. associated(FirstDischarges)) then
            GetNumberOfDischarges = - 99

        else            
        
            call Ready(DischargesID, ready_)    
        
            if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
                GetNumberOfDischarges = Me%DischargesNumber
            else 
                GetNumberOfDischarges = - 99
            end if
       
        end if 
        
        return
 
    end function GetNumberOfDischarges

    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetDischargeName
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETDISCHARGENAME"::GetDischargeName
    !DEC$ ENDIF
    logical function GetDischargeName(DischargeID, DischargeNumber, DischargeName)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DischargeID
        integer                                     :: DischargeNumber
        character(len=*)                            :: DischargeName        
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
        type(T_IndividualDischarge), pointer        :: DischargeX

        call Ready(DischargeID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeNumber, '.'
                stop       'Subroutine GetDischargeType - ModuleDischarges. ERR01.'
            endif

            DischargeName =  trim(DischargeX%ID%Name)

            GetDischargeName = .true.
        else 
            GetDischargeName = .false.
        end if
           
        return

    end function GetDischargeName
    
    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetDischargeType
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETDISCHARGETYPE"::GetDischargeType
    !DEC$ ENDIF
    integer function GetDischargeType(DischargeID, DischargeNumber)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DischargeID
        integer                                     :: DischargeNumber
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
        type(T_IndividualDischarge), pointer        :: DischargeX

        call Ready(DischargeID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeNumber, '.'
                stop       'Subroutine GetDischargeType - ModuleDischarges. ERR01.'
            endif

            GetDischargeType = DischargeX%DischargeType
        else 
            GetDischargeType = - 99
        end if
           
        return

    end function GetDischargeType

    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetDischargeXCoordinate
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETDISCHARGEXCOORDINATE"::GetDischargeXCoordinate
    !DEC$ ENDIF
    real(8) function GetDischargeXCoordinate(DischargeID, DischargeNumber)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DischargeID
        integer                                     :: DischargeNumber
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
        type(T_IndividualDischarge), pointer        :: DischargeX

        call Ready(DischargeID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeNumber, '.'
                stop       'Subroutine GetDischargeXCoordinate - ModuleDischarges. ERR01.'
            endif

            GetDischargeXCoordinate = DischargeX%Localization%CoordinateX
        else 
            GetDischargeXCoordinate = - 99.0
        end if
           
        return

    end function GetDischargeXCoordinate
    
    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetDischargeYCoordinate
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETDISCHARGEYCOORDINATE"::GetDischargeYCoordinate
    !DEC$ ENDIF
    real(8) function GetDischargeYCoordinate(DischargeID, DischargeNumber)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DischargeID
        integer                                     :: DischargeNumber
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
        type(T_IndividualDischarge), pointer        :: DischargeX

        call Ready(DischargeID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeNumber, '.'
                stop       'Subroutine GetDischargeYCoordinate - ModuleDischarges. ERR01.'
            endif

            GetDischargeYCoordinate = DischargeX%Localization%CoordinateY
        else 
            GetDischargeYCoordinate = - 99.0
        end if
           
        return

    end function GetDischargeYCoordinate
    
    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::SetDischargeFlow
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_SETDISCHARGEFLOW"::SetDischargeFlow
    !DEC$ ENDIF
    logical function SetDischargeFlow(DischargeID, DischargeNumber, Flow)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DischargeID
        integer                                     :: DischargeNumber
        real(8)                                     :: Flow
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
        type(T_IndividualDischarge), pointer        :: DischargeX

        call Ready(DischargeID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeNumber, '.'
                stop       'Subroutine SetDischargeFlow - ModuleDischarges. ERR01.'
            endif
            
            DischargeX%WaterFlow%Scalar = Flow
           
            SetDischargeFlow = .true.

        else
        
            SetDischargeFlow = .false.

        end if
           
        return

    end function SetDischargeFlow

   
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetNumberOfDischargeProperties
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETNUMBEROFDISCHARGEPROPERTIES"::GetNumberOfDischargeProperties
    !DEC$ ENDIF
    !Return the number of Error Messages
    integer function GetNumberOfDischargeProperties(DischargeID, DischargeNumber)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DischargeID
        integer                                     :: DischargeNumber
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_
        type(T_IndividualDischarge), pointer        :: DischargeX
                 

        call Ready(DischargeID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeNumber, '.'
                stop       'Subroutine GetNumberOfDischargeProperties - ModuleDischarges. ERR01.'
            endif
            GetNumberOfDischargeProperties = DischargeX%PropertiesNumber
        else 
            GetNumberOfDischargeProperties = -99
        end if
           
    end function GetNumberOfDischargeProperties
    
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetDischargePropertyID
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETDISCHARGEPROPERTYID"::GetDischargePropertyID
    !DEC$ ENDIF
    !Return the number of Error Messages
    integer function GetDischargePropertyID(DischargeID, DischargeNumber, idx)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DischargeID
        integer                                     :: DischargeNumber
        integer                                     :: idx
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
        type (T_Property), pointer                  :: Property
        integer                                     :: iProp
        type(T_IndividualDischarge), pointer        :: DischargeX

        call Ready(DischargeID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeNumber, '.'
                stop       'Subroutine GetDischargePropertyID - ModuleDischarges. ERR01.'
            endif

            Property => DischargeX%FirstProperty
            iProp = 1
            do while (associated (Property))
                 
                 if (iProp == idx) then
                 
                    GetDischargePropertyID = Property%ID%IDNumber
                    return
                 
                 endif
                 
                 Property => Property%Next
                 iProp = iProp + 1
            enddo
        
            GetDischargePropertyID = -99
        else 
            GetDischargePropertyID = -99
        end if
           
    end function GetDischargePropertyID
    

    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::SetDischargeConcentration
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_SETDISCHARGECONCENTRATION"::SetDischargeConcentration
    !DEC$ ENDIF
    logical function SetDischargeConcentration(DischargeID, DischargeNumber, PropertyIDNumber, Concentration)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DischargeID
        integer                                     :: DischargeNumber
        integer                                     :: PropertyIDNumber
        real(8)                                     :: Concentration
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
        type(T_IndividualDischarge), pointer        :: DischargeX
        type(T_Property), pointer                   :: PropertyX

        call Ready(DischargeID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeNumber)
            if (STAT_CALL/=SUCCESS_) then 
                write(*,*) 'Can not find discharge number ', DischargeNumber, '.'
                stop       'Subroutine SetDischargeConcentration - ModuleDischarges. ERR01.'
            endif
            
            call Search_Property(DischargeX, PropertyX, STAT_CALL, PropertyXIDNumber=PropertyIDNumber)
            if (STAT_CALL==SUCCESS_) then
            
                PropertyX%Scalar = Concentration
            
                SetDischargeConcentration = .true.
            else
            
                SetDischargeConcentration = .false.
            endif 

        else
        
            SetDischargeConcentration = .false.

        end if
           
        return

    end function SetDischargeConcentration    
    
    !--------------------------------------------------------------------------

#endif


end module ModuleDischarges

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

