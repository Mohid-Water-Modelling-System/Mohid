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
    use ModuleTimeSerie, only : StartTimeSerieInput, GetTimeSerieValue, KillTimeSerie
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


    !Selector
    public  :: GetDischargesNumber                                  !Returns the number of discharge points
    public  :: GetDischargesGridLocalization
    public  :: GetDischargesNodeID
    public  :: GetDischargesIDName
    public  :: GetDischargeWaterFlow
    public  :: GetDischargeFlowVelocity
    public  :: GetDischargeParameters
    public  :: GetDischargeConcentration
    public  :: GetByPassON
    private ::    Search_Discharge
    private ::    Search_Property
    public  :: GetDischargeSpatialEmission
    public  :: GetDischargeFlowDistribuiton
    public  :: GetDischargeON
    public  :: SetLocationCellsZ
    public  :: SetLayer
    public  :: UngetDischarges
    !Modifier
    public  :: CorrectsCellsDischarges
    public  :: TryIgnoreDischarge
    !Destructor
    public  :: Kill_Discharges
    private ::      KillIndividualDischarge
    private ::      DeallocateInstance


    !Management
    private ::      Ready


    private :: UngetDischarges1Dinteger
    interface  UngetDischarges
        module procedure UngetDischarges1Dinteger
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

    !Valve side
    integer, parameter :: SideA         = 1
    integer, parameter :: SideB         = 2


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
        logical                                 :: Variable       = .false.
        integer                                 :: ConcColumn
        real                                    :: scalar         = FillValueReal
        logical                                 :: TimeSerieON
        integer                                 :: TimeSerie      = 0
        logical                                 :: PropTimeSerie  = .false. 
        type (T_Property), pointer              :: Next,Prev
    end type T_Property

    type       T_WaterFlow
        logical                                 :: Variable      = .false.
        integer                                 :: FlowColumn
        real                                    :: scalar        = FillValueReal
    end type T_WaterFlow


    type       T_WaterVelocity
        logical                                 :: UVariable
        logical                                 :: VVariable
        integer                                 :: UColumn               
        integer                                 :: VColumn               
        real                                    :: Uscalar = FillValueReal
        real                                    :: Vscalar = FillValueReal
    end type T_WaterVelocity

    type       T_FlowOver
        real                                    :: WeirLength
        real                                    :: DischargeCoeficient
        real                                    :: CrestHeigth
    end  type T_FlowOver

    type       T_Valve
        real                                    :: Diameter
        real                                    :: DischargeCoeficient
        real                                    :: AxisHeigth
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
         real                                   :: MinimumDischargeDepth
         logical                                :: StartFromLastDischargeLocation   = .false.
         logical                                :: Location2D
         integer                                :: DischVertical                    = FillValueInt
         real                                   :: Kdepth                           = FillValueReal
         integer                                :: NodeID = FillValueInt
         logical                                :: TrackLocation
         character(len=StringLength)            :: TrackLocationFile
         integer                                :: TrackLocationFileUnitNumber       
         logical                                :: UseDischargePathFile             = .false.
         character(len=StringLength)            :: DischargePathFile
         real                                   :: CoordinateX, CoordinateY
         logical                                :: CoordinatesON                    = .false.
         integer                                :: HorizontalType                   = FillValueInt
         logical                                :: CellCorrect                      = .false.
         integer                                :: SpatialEmission
         character(len=StringLength)            :: SpatialFile
         type (T_Polygon), pointer              :: Polygon
         type (T_Lines),   pointer              :: Line
         integer                                :: nCells                           = 1
         integer, dimension(:), pointer         :: VectorI, VectorJ, VectorK
         integer                                :: FlowDistribution                 = DischByCell_

    end  type T_Localization


    type       T_ByPass 
         integer                                :: i, j
         logical                                :: ON, OneWay
         integer                                :: Side
    end  type T_ByPass


    type       T_IndividualDischarge
         type(T_ID                 )            :: ID
         type(T_Localization       )            :: Localization
         integer                                :: PropertiesNumber = FillValueInt
         character(len=PathLength)              :: DataBaseFile
         logical                                :: TimeSerieON
         integer                                :: TimeSerie        = 0
         type(T_WaterFlow          )            :: WaterFlow   
         type(T_WaterVelocity      )            :: VelocityFlow
         integer                                :: DischargeType
         type(T_Valve   )                       :: Valve
         type(T_FlowOver)                       :: FlowOver
         type(T_Property           ), pointer   :: FirstProperty    => null()
         type(T_Property           ), pointer   :: LastProperty     => null()
         type(T_Property           ), pointer   :: CurrProperty     => null()
         type(T_IndividualDischarge), pointer   :: Next             => null()
         type(T_IndividualDischarge), pointer   :: Prev             => null()
         type(T_ByPass             )            :: ByPass           
         logical                                :: IgnoreON
    end type T_IndividualDischarge    

    type      T_Discharges
         integer                                :: InstanceID            
         integer                                :: ObjEnterData     = 0
         integer                                :: ObjTime          = 0
         character(len=Pathlength)              :: DataFile
         integer                                :: DischargesNumber= FillValueInt
         type (T_IndividualDischarge), pointer  :: FirstDischarge   => null()
         type (T_IndividualDischarge), pointer  :: LastDischarge    => null()
         type (T_IndividualDischarge), pointer  :: CurrentDischarge => null()
         type (T_Discharges), pointer           :: Next
         logical                                :: IgnoreON
    end type T_Discharges

    !Global Variables
    type (T_Discharges), pointer                :: FirstDischarges
    type (T_Discharges), pointer                :: Me
    
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
                if (STAT_CALL /= SUCCESS_) stop 'Construct_Discharges - ModuleDischarges - ERR10'
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
                stop 'Subroutine Construct_Discharges; ModuleDischarges. ERR30.' 
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


            ! Constructs the discharge list 
            call Construct_DischargeList


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
                        stop 'Subroutine Construct_DischargeList; ModuleDischarges. ERR01.'

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                stop 'Subroutine Construct_DischargeList; ModuleDischarges. ERR02.'
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
        call Construct_FlowValues             (NewDischarge)

        !Construct Discharge Velocity values
        call Construct_VelocityValues         (NewDischarge)

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
                stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR120.'
            if (flag==0)                                                            &
                stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR130.'

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
                    stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR20.'
                if (flag==0)                                                                &
                    stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR30.'
        
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
                stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR40.'

            call GetData(AuxLog,                                                        &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='DISCHARGE_UNIFORM',                             &
                         ClientModule = 'ModuleDischarges',                             &
                         default      = .false.,                                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR50.'

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
                stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR60.'
            
            if (NewDischarge%Localization%DischVertical .eq. DischLayer_ .and. flag .eq. 0)then
                write(*,*)"You must define the K_CELL in the discharge: "//trim(NewDischarge%ID%Name)
                stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR61.'
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
                        stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR70.'

                case (DischBottom_, DischSurf_, DischUniform_)
                    !do not do nothing 
                case default
                    write(*,*) "VERTICAL DISCHARGE option not known ", NewDischarge%Localization%DischVertical

                    write(*,*) "The known options are : "," Bottom=",DischBottom_," Surface=",DischSurf_,&
                                                          " Layer =",DischLayer_, " Depth  =",DischDepth_,&
                                                          " Uniform=",DischUniform_
                    stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR80'

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
                stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR90.'

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
                    stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR100'


                 call GetData(NewDischarge%Localization%StartFromLastDischargeLocation, &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             FromBlock,                                                 &
                             keyword      ='START_FROM_LAST_LOCATION',                  &
                             ClientModule = 'ModuleDischarges',                         &
                             Default      = .false.,                                    &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR110'
                                    
                
                call GetData(NewDischarge%Localization%TrackLocation,                   &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             FromBlock,                                                 &
                             keyword      ='TRACKLOCATION',                             &
                             ClientModule ='ModuleDischarges',                          &
                             Default      = .false.,                                    &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR120'
                
                if (NewDischarge%Localization%TrackLocation) then
                    
                    !Gets the root path from the file nomfich.dat
                    call ReadFileName("ROOT_SRT", RootPath, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) then
                        call ReadFileName("ROOT", RootPath, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) then
                            call ReadFileName("RAIZ", RootPath, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                              &
                                stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR130'
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
                        stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR140'
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
                stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR150'

            select case (AuxName)

                case ("Point")
                    NewDischarge%Localization%SpatialEmission = DischPoint_
                case ("Line")
                    NewDischarge%Localization%SpatialEmission = DischLine_
                case ("Polygon")
                    NewDischarge%Localization%SpatialEmission = DischPolygon_
                case default
                    write(*,*) "SPATIAL EMISSION option not known ",trim(AuxName)," ????"
                    write(*,*) "The known options are : ","Point ", "Line ", "Polygon"
                    stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR160'

            end select

            if (NewDischarge%Localization%SpatialEmission == DischLine_ .or.            &
                NewDischarge%Localization%SpatialEmission == DischPolygon_) then
                call GetData(NewDischarge%Localization%SpatialFile,                     &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='SPATIAL_FILE',                                  &
                         ClientModule = 'ModuleDischarges',                             &
                         Default      = "Point",                                        &
                         STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR170'

                if (flag == 0)                                                          &
                    stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR180'

                if (NewDischarge%Localization%SpatialEmission == DischLine_)            &
                    call New(NewDischarge%Localization%Line, NewDischarge%Localization%SpatialFile)

                if (NewDischarge%Localization%SpatialEmission == DischPolygon_)         &
                    call New(NewDischarge%Localization%Polygon, NewDischarge%Localization%SpatialFile)

                call GetData(AuxName,                                                   &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='FLOW_DISTRIBUTION',                             &
                         ClientModule = 'ModuleDischarges',                             &
                         Default      = "by cell",                                      &
                         STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR190'

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
                        stop 'Subroutine ConstDischargeLoc; Module ModuleDischarges. ERR200'

                end select

            endif

        else    i1

            !I Location
            call GetData(NewDischarge%Localization%NodeID,                              &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='NODE_ID',                                       &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstDischargeLoc - ModuleDischarges - ERR210'
            
        endif   i1


        !----------------------------------------------------------------------

     End Subroutine ConstDischargeLoc

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine Read_DataBaseFile(NewDischarge)

        !Arguments-------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: NewDischarge

        !External--------------------------------------------------------------
        integer                                     :: flag, STAT_CALL

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
            call StartTimeSerieInput(NewDischarge%TimeSerie, NewDischarge%DataBaseFile,  &
                                     Me%ObjTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Read_DataBaseFile - ModuleDischarges - ERR20'
        end if 

        !----------------------------------------------------------------------

    end subroutine Read_DataBaseFile 

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    
    subroutine Construct_FlowValues(NewDischarge)

        !Arguments--------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: NewDischarge

        !Local-----------------------------------------------------------------
        integer                                     :: flag, STAT_CALL

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
            else
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

        if      (NewDischarge%DischargeType /= Normal   .and.                           &
                 NewDischarge%DischargeType /= FlowOver .and.                           &
                 NewDischarge%DischargeType /= Valve   ) then
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

            
        else if (NewDischarge%DischargeType == Valve) then i2

            call GetData(NewDischarge%Valve%Diameter,                                   &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='VALVE_DIAMETER',                                &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR100'
            
            if (flag /= 1) then
                write(*,*) 'Valve diameter missing'
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

            call GetData(NewDischarge%Valve%AxisHeigth,                                 &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='VALVE_AXIS_HEIGTH',                             &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR130'

            if (flag /= 1) then
                write(*,*) 'Valve axis Missing'
                stop ' Construct_FlowValues - ModuleDischarges - ERR140'
            endif

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

i3:     if (NewDischarge%ByPass%ON) then

            call GetData(NewDischarge%ByPass%i,                                         &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='BYPASS_I',                                      &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR150'

            if (flag /= 1) then
                write(*,*) 'Bypass I cell missing'
                stop ' Construct_FlowValues - ModuleDischarges - ERR160'
            endif

            call GetData(NewDischarge%ByPass%j,                                         &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='BYPASS_J',                                      &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR160'

            if (flag /= 1) then
                write(*,*) 'Bypass J cell missing'
                stop ' Construct_FlowValues - ModuleDischarges - ERR170'
            endif


            call GetData(NewDischarge%ByPass%Side,                                      &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='BYPASS_SIDE',                                   &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR180'

            if (flag /= 1) then
                write(*,*) 'Bypass side missing'
                stop ' Construct_FlowValues - ModuleDischarges - ERR190'
            endif

            if      (NewDischarge%ByPass%Side /= SideA   .and.                          &
                     NewDischarge%ByPass%Side /= SideB) then
                     stop 'Construct_FlowValues - ModuleDischarges - ERR200'
            endif

            call GetData(NewDischarge%ByPass%OneWay,                                    &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         FromBlock,                                                     &
                         keyword      ='BYPASS_ONEWAY',                                 &
                         default      = .false.,                                        &
                         ClientModule = 'ModuleDischarges',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_FlowValues - ModuleDischarges - ERR210'

        endif i3         
        !----------------------------------------------------------------------

        !----------------------------------------------------------------------

     End Subroutine Construct_FlowValues  

    !--------------------------------------------------------------------------

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
        allocate(aux(2))
        call GetData(Aux, Me%ObjEnterData, flag,                                         &
                     keyword    = 'DEFAULT_VELOCITY_VALUE',                              &
                     default    = 0.,                                                    &
                     SearchType = FromBlock,                                             &
                     ClientModule = 'ModuleDischarges',                                  &
                     STAT       = STAT_CALL)

        NewDischarge%VelocityFlow%UScalar = Aux(1)
        NewDischarge%VelocityFlow%VScalar = Aux(2)
        deallocate(aux)

        !Searches for the velocity U column (if the discharge is read from a data base)
        if (NewDischarge%TimeSerieON) then
            call GetData(NewDischarge%VelocityFlow%UColumn,                              &
                         Me%ObjEnterData, flag,                                          &
                         keyword    = 'U_COLUMN',                                        &
                         default    = FillValueInt,                                      &
                         SearchType = FromBlock,                                         &
                         ClientModule = 'ModuleDischarges',                              &
                         STAT       = STAT_CALL)

            if (flag == 1) then
                NewDischarge%VelocityFlow%UVariable = .true.            
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

            if (flag == 1) then
                NewDischarge%VelocityFlow%VVariable = .true.
            else
                NewDischarge%VelocityFlow%VVariable = .false.
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
                stop 'Subroutine Construct_PropertyList; ModuleDischarges. ERR02.'
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

        if (NewProperty%Variable) then

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

            else
                if (NewDischarge%TimeSerieON) then
                    NewProperty%TimeSerie = NewDischarge%TimeSerie
                else
                    stop 'Construct_PropertyValues - ModuleDischarges - ERR20'
                endif
            endif

        endif

        !By default the property value in the domain is zero
        call GetData(NewProperty%scalar,                                                 &
                     Me%ObjEnterData,                                                    &
                     flag,                                                               &
                     FromBlockinBlock,                                                   &
                     keyword      ='DEFAULTVALUE',                                       &
                     ClientModule ='ModuleDischarges',                                   &
                     default      = 0.0)
                     
        if((.not. NewProperty%Variable) .and. (flag .eq. 0))then
            write(*,*)"Please define concentration default value for property :"
            write(*,*)trim(NewProperty%ID%Name)
            write(*,*)"Discharge name : ", trim(NewDischarge%ID%Name)
            stop 'Construct_PropertyValues - ModuleDischarges -  ERR30'
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

    subroutine ConstructLog

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
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

    subroutine GetDischargesGridLocalization(DischargesID, DischargeIDNumber, Igrid,        &
                                             JGrid, KGrid, IByPass, JByPass, DischVertical, &
                                             KDepth, WaterColumnZ, Bathymetry, OpenPoints3D, &
                                             CoordinateX, CoordinateY, CoordinatesON, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: DischargesID
        integer,           intent(IN )                  :: DischargeIDNumber
        integer, optional, intent(OUT)                  :: IGrid, JGrid, KGrid
        integer, optional, intent(OUT)                  :: IByPass, JByPass
        integer, optional, intent(OUT)                  :: DischVertical
        integer, optional, dimension(:, :, :), pointer  :: OpenPoints3D
        real   , optional, dimension(:, :   ), pointer  :: WaterColumnZ
        real   , optional, dimension(:, :   ), pointer  :: Bathymetry
        
        real,    optional, intent(OUT)                  :: CoordinateX, CoordinateY, Kdepth
        logical, optional, intent(OUT)                  :: CoordinatesON        
                
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
                stop       'Subroutine GetDischargesGridLocalization; Module ModuleDischarges. ERR01.'
            endif

            if (.not. DischargeX%Localization%Location2D) then
                write(*,*)'Discharge location not given in Grid Coordinates'
                write(*,*)trim(adjustl(DischargeX%ID%Name))
                stop 'GetDischargesGridLocalization - ModuleDischarges - ERR01a'
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


            !If OpenPoints is present an DischargeX have alternative locations, 
            !check out which locations to use.
            if (present(WaterColumnZ)) then

                WorkSizeILB = lbound(WaterColumnZ, dim = 1) + 1
                WorkSizeIUB = ubound(WaterColumnZ, dim = 1) - 1

                WorkSizeJLB = lbound(WaterColumnZ, dim = 2) + 1
                WorkSizeJUB = ubound(WaterColumnZ, dim = 2) - 1

                if (IGrid < WorkSizeILB .or. IGrid > WorkSizeIUB)                       &
                    stop  'GetDischargesGridLocalization; ModuleDischarges. ERR02.'

                if (JGrid < WorkSizeJLB .or. JGrid > WorkSizeJUB)                       &
                    stop  'GetDischargesGridLocalization; ModuleDischarges. ERR03.'
                
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
                stop       'Subroutine GetDischargesGridLocalization; Module ModuleDischarges. ERR01.'
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
                stop       'GetDischargesGridLocalization - ModuleDischarges - ERR01'
            endif

            if (DischargeX%Localization%Location2D) then
                write(*,*)'Discharge location not given as Node ID'
                write(*,*)trim(adjustl(DischargeX%ID%Name))
                stop 'GetDischargesGridLocalization - ModuleDischarges - ERR01a'
            endif

            NodeID = DischargeX%Localization%NodeID

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_


    end subroutine GetDischargesNodeID

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    Subroutine GetDischargeSpatialEmission(DischargesID, DischargeIDNumber, LineX, PolygonX, &
                                           SpatialEmission, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer,                        intent(IN ) :: DischargeIDNumber
        type (T_Lines),  pointer                    :: LineX
        type (T_Polygon), pointer                   :: PolygonX
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
                stop       'Subroutine GetDischargeSpatialEmission; Module ModuleDischarges. ERR01.'
            end if cd3

            SpatialEmission = DischargeX%Localization%SpatialEmission

            PolygonX => DischargeX%Localization%Polygon
            LineX    => DischargeX%Localization%Line

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


    !--------------------------------------------------------------------------

    Subroutine GetDischargeFlowDistribuiton(DischargesID, DischargeIDNumber,            &
                                            nCells, FlowDistribution,                   &
                                            VectorI, VectorJ, VectorK, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        integer,                        intent(IN ) :: DischargeIDNumber
        integer,                        intent(OUT) :: nCells
        integer, optional,              intent(OUT) :: FlowDistribution
                        
        integer, dimension(:), pointer, optional    :: VectorI, VectorJ, VectorK
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
                    stop       'Subroutine GetDischargeFlowDistribution; Module ModuleDischarges. ERR01.'
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
    Subroutine GetDischargeWaterFlow(DischargesID, TimeX, DischargeIDNumber,            &
                                     SurfaceElevation, Flow, SurfaceElevation2, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        type(T_Time),                   intent(IN ) :: TimeX
        integer,                        intent(IN ) :: DischargeIDNumber
        real   ,                        intent(IN)  :: SurfaceElevation
        real   ,                        intent(OUT) :: Flow
        real   , optional,              intent(IN)  :: SurfaceElevation2
        integer, optional,              intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_
        integer                                     :: STAT_
        type(T_IndividualDischarge), pointer        :: DischargeX
        type(T_Time)                                :: Time1, Time2
        real                                        :: Value1, Value2, NewValue
        real                                        :: H, C, A
        logical                                     :: TimeCycle
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then



            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber=DischargeIDNumber)

cd3 :       if (STAT_CALL/=SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                    stop       'Subroutine GetDischargeWaterFlow; Module ModuleDischarges. ERR01.'
            end if cd3


cd2:        if (DischargeX%DischargeType == Normal .and. DischargeX%WaterFlow%Variable) then

                call GetTimeSerieValue(DischargeX%TimeSerie, TimeX,                      &
                                       DischargeX%WaterFlow%FlowColumn, Time1, Value1,   &
                                       Time2, Value2, TimeCycle, STAT)

                if (TimeCycle) then
                    NewValue = Value1
                else
                    !Interpolates Value for current instant
                    call InterpolateValueInTime(TimeX, Time1, Value1, Time2, Value2,     &
                                                NewValue)
                endif

                Flow = NewValue

            elseif (DischargeX%DischargeType == FlowOver) then

                !Q = cv * b * sqrt(2*g) * H^(1.5)
                H = SurfaceElevation - DischargeX%FlowOver%CrestHeigth
                if (H > 0) then
                    Flow = -sqrt(19.6) * DischargeX%FlowOver%DischargeCoeficient *       &
                                         DischargeX%FlowOver%WeirLength  * H ** 1.5
                else
                    Flow = 0.
                endif

            elseif (DischargeX%DischargeType == Valve) then


                !Q = C * A * sqrt(2*g) * sqrt(H)
                H = SurfaceElevation - SurfaceElevation2

                !if the axis valve is above the water level in both sides than there is no flow
                if (- DischargeX%Valve%AxisHeigth > max(SurfaceElevation, SurfaceElevation2)) then
                
                    H = 0.

                else

                    A    = Pi * (DischargeX%Valve%Diameter/2.)**2
                    C    = DischargeX%Valve%DischargeCoeficient

                    Flow = - sqrt(19.6) * C * A * sqrt(abs(H))

                    if (H<0.) Flow = - Flow
                
                    if     (DischargeX%ByPass%OneWay) then

                        if ((DischargeX%ByPass%Side == SideA .and. H < 0.) .or.         &
                            (DischargeX%ByPass%Side == SideB .and. H > 0.)) Flow = 0.
                    endif
                endif

            else

                Flow = DischargeX%WaterFlow%Scalar

            end if cd2

            nullify(DischargeX)


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end Subroutine GetDischargeWaterFlow

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetDischargeFlowVelocity(DischargesID, TimeX, DischargeIDNumber,          &
                                        VelocityU, VelocityV ,STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        real,    optional, intent(OUT)              :: VelocityU, VelocityV
        integer, optional, intent(OUT)              :: STAT
        integer,           intent(IN )              :: DischargeIDNumber
        type(T_Time),      intent(IN )              :: TimeX

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_
        type(T_IndividualDischarge), pointer        :: DischargeX
        type(T_Time)                                :: Time1, Time2
        real                                        :: Value1, Value2, NewValue
        logical                                     :: TimeCycle
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DischargesID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (.not. present(VelocityU) .and. .not. present(VelocityV))                 &
                stop 'Subroutine GetDischargeFlowVelocity; Module ModuleDischarges. ERR01.'

            call Search_Discharge(DischargeX, STAT_CALL, DischargeXIDNumber = DischargeIDNumber)
                

cd3 :       if (STAT_CALL /= SUCCESS_) then 
                write(*,*) ' can not find discharge number ',DischargeIDNumber
                    stop 'Subroutine GetDischargeFlowVelocity; Module ModuleDischarges. ERR02.'
            end if cd3


cd4 :       if (Present(VelocityU)) then 

cd2:            if (DischargeX%VelocityFlow%UVariable) then

                    call GetTimeSerieValue(DischargeX%TimeSerie, TimeX,                      &
                                           DischargeX%VelocityFlow%UColumn, Time1, Value1,   &
                                           Time2, Value2, TimeCycle, STAT)

                    if (TimeCycle) then
                        NewValue = Value1
                    else
                        !Interpolates Value for current instant
                        call InterpolateValueInTime(TimeX, Time1, Value1, Time2, Value2,     &
                                                    NewValue)
                    endif

                    VelocityU = NewValue

                else

                    VelocityU = DischargeX%VelocityFlow%UScalar

                end if cd2

            endif cd4

cd5 :       if (Present(VelocityV)) then 

cd6:            if (DischargeX%VelocityFlow%VVariable) then

                    call GetTimeSerieValue(DischargeX%TimeSerie, TimeX,                      &
                                           DischargeX%VelocityFlow%VColumn, Time1, Value1,   &
                                           Time2, Value2, TimeCycle, STAT)

                    if (TimeCycle) then
                        NewValue = Value1
                    else
                        !Interpolates Value for current instant
                        call InterpolateValueInTime(TimeX, Time1, Value1, Time2, Value2,     &
                                                    NewValue)
                    endif

                    VelocityV = NewValue

                else
                    VelocityV = DischargeX%VelocityFlow%VScalar

                end if cd6

            endif cd5


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
                    stop 'Subroutine GetDischargeON; Module ModuleDischarges. ERR10.'
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
    subroutine GetDischargeConcentration(DischargesID, TimeX,DischargeIDNumber,          &
                                         Concentration, PropertyIDNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: DischargesID
        type(T_Time),      intent(IN)               :: TimeX
        integer,           intent(IN)               :: DischargeIDNumber
        real,              intent(OUT)              :: Concentration
        integer,           intent(IN)               :: PropertyIDNumber
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        type(T_IndividualDischarge), pointer        :: DischargeX
        type(T_Property),            pointer        :: PropertyX
        type(T_Time)                                :: Time1, Time2
        real                                        :: Value1, Value2, NewValue
        logical                                     :: TimeCycle
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
                stop  'Sub. GetDischargeConcentration - ModuleDischarges - ERR01'
            endif
             

            call Search_Property(DischargeX, PropertyX, STAT_CALL, PropertyXIDNumber=PropertyIDNumber)

            if (STAT_CALL/=SUCCESS_) then 
                !If the proeprty is not found the program don't stop is return a error 
                !not found
                if (STAT_CALL /= NOT_FOUND_ERR_) then 
                    stop  'Sub. GetDischargeConcentration - ModuleDischarges - ERR05'
                endif

            endif


            
cd2 :       if (STAT_CALL == SUCCESS_) then
                                                                         
                if (PropertyX%Variable) then

                    call GetTimeSerieValue(PropertyX%TimeSerie, TimeX,                   &
                                           PropertyX%ConcColumn, Time1, Value1,          &
                                           Time2, Value2, TimeCycle, STAT)

                    if (TimeCycle) then
                        NewValue = Value1
                    else
                        !Interpolates Value for current instant
                        call InterpolateValueInTime(TimeX, Time1, Value1, Time2, Value2, &
                                                    NewValue)
                    endif

                    Concentration = NewValue

                else
 
                    Concentration = PropertyX%Scalar

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

end module ModuleDischarges

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

