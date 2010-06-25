!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Land, Mohid SurfaceWater
! MODULE        : Atmosphere
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : July 2003
! REVISION      : Frank & Pedro Chambel Leitao, Luis Fernandes - v4.0
! DESCRIPTION   : Module used to read and calculate values for atmosphere 
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

Module ModuleAtmosphere

    use ModuleGlobalData
    use ModuleTime
    use ModuleHDF5
    use ModuleFunctions,      only : ConstructPropertyID, CHUNK_J
    use ModuleFillMatrix,     only : ConstructFillMatrix, ModifyFillMatrix, KillFillMatrix,  &
                                     GetIfMatrixRemainsConstant, GetFillMatrixDTPrediction
    use ModuleTimeSerie,      only : StartTimeSerie, WriteTimeSerie, KillTimeSerie,     &
                                     GetTimeSerieLocation, CorrectsCellsTimeSerie,      &
                                     GetNumberOfTimeSeries, TryIgnoreTimeSerie
    use ModuleEnterData,      only : ReadFileName, ConstructEnterData, GetData,              &
                                     ExtractBlockFromBuffer, Block_Unlock, GetOutPutTime,    &
                                     KillEnterData
    use ModuleGridData,       only : GetGridData, UngetGridData         
    use ModuleHorizontalGrid, only : GetHorizontalGrid, GetHorizontalGridSize, GetGridAngle, &
                                     GetGridLatitudeLongitude, WriteHorizontalGrid,          &
                                     UnGetHorizontalGrid, RotateVectorFieldToGrid,           &
                                     GetGridCellArea, GetXYCellZ
    use ModuleStatistic,      only : ConstructStatistic, GetStatisticMethod,                 &
                                     GetStatisticParameters, ModifyStatistic, KillStatistic
    use ModuleStopWatch,      only: StartWatch, StopWatch

    implicit none

    private 

    !Subroutines & Functions---------------------------------------------------

    !Constructor
    public  :: StartAtmosphere
    private ::      AllocateInstance
    private ::      ReadPropertiesFilesName
    private ::      ConstructPropertyList
    private ::          ConstructProperty  
    private ::          Add_Property  
    private ::              Construct_PropertyValues
    private ::              ConstructOutPut
    private ::              ConstructSurfStatistics
    private ::      OpenHDF5OutPutFile

    private ::      ConstructTimeSerie

    private ::      RotateAtmosphereVectorFields


    !Selector
    public  :: GetAtmosphereProperty
    private ::      SearchProperty
    public  :: AtmospherePropertyExists
    public  :: GetAtmospherenProperties
    public  :: GetAtmospherePropertiesIDByIdx
    public  :: GetAtmosphereDTPrediction

    public  :: UngetAtmosphere

    private :: ReadLockExternalVar
    private :: ReadUnlockExternalVar

    !Modifier
    public  ::  ModifyAtmosphere
    private ::      ModifyAirTemperature
    private ::      ModifyWindVelocity
    private ::          ComputeWindVelocity
    private ::      ModifyWindDirection
    private ::      ModifyWindModulus
    private ::      ModifyPrecipitation
    private ::      ModifySolarRadiation
    private ::          ClimatologicSolarRadiation            
    private ::              RandomCloud                   !Function
    private ::              TOARadiation                  !Function                                    
    private ::      ModifyAtmosphericPressure
    private ::      ModifyCo2AtmosphericPressure
    private ::      ModifyO2AtmosphericPressure
    private ::      ModifyPropByIrri
    private ::      ModifyPropByRain
    private ::      ModifySunHours
    private ::      ModifyCloudCover
    private ::      ModifyRelativeHumidity
    private ::      ModifyIrrigation
    private ::      ModifyRandom
    private ::      ModifyOutput
    private ::          OutPutResultsHDF5
    private ::          OutPut_TimeSeries
    private ::          OutPut_Statistics

    !Destructor
    public  ::  KillAtmosphere
    private ::      DeallocateInstance
    private ::      DeallocateVariables  
    

    !Management
    private ::      Ready
    private ::      LocateObjAtmosphere

    !Interfaces----------------------------------------------------------------
    private :: UngetAtmosphere1D
    private :: UngetAtmosphere2D
    interface  UngetAtmosphere
        module procedure UngetAtmosphere1D
        module procedure UngetAtmosphere2D
    end interface UngetAtmosphere

    !Parameter-----------------------------------------------------------------
    
    !Sun constant (W / m**2)
    real,    parameter                              :: KSun                 = 1367.0     

    character(LEN = StringLength), parameter        :: block_begin          = '<beginproperty>'
    character(LEN = StringLength), parameter        :: block_end            = '<endproperty>'

    !Parameter
    integer, parameter                              :: Radiation_MOHID      = 1
    integer, parameter                              :: Radiation_CEQUALW2   = 2

    integer, parameter                              :: CloudFromSunHours    = 1
    integer, parameter                              :: CloudFromRandom      = 2
    integer, parameter                              :: CloudFromRadiation   = 3



    !Types---------------------------------------------------------------------
    type       T_External
        integer, dimension(:,:), pointer            :: MappingPoints2D      => null()
        real,    dimension(:,:), pointer            :: GridCellArea         => null()
        real,    dimension(:,:), pointer            :: Latitude             => null()
        real,    dimension(:,:), pointer            :: Longitude            => null()
    end type   T_External


    type       T_OutPut
         type (T_Time), pointer, dimension(:)       :: OutTime              => null()
         logical                                    :: True
         integer                                    :: NextHDF5
         integer                                    :: Number
    end type T_OutPut


    type       T_Statistics
         integer                                    :: ID, IDx, IDy
         character(LEN = Pathlength)                :: File
         logical                                    :: ON
    end type T_Statistics


    type       T_Property
        type (T_PropertyID)                         :: ID
        real, dimension(:,:), pointer               :: Field                => null()
        real, dimension(:,:), pointer               :: FieldGrid            => null()
        real                                        :: RandomValue
        logical                                     :: HasRandomComponent   = .false.
        logical                                     :: PropAddedByIrri      = .false.
        logical                                     :: PropAddedByRain      = .false.
        logical                                     :: FirstActualization   = .true.
        real                                        :: RandomComponent      = FillValueReal
        logical                                     :: TimeSerie            = .false.
        logical                                     :: BoxTimeSerie         = .false.
        logical                                     :: OutputHDF            = .false.
        logical                                     :: Constant             = .false.
        logical                                     :: NoInterpolateValueInTime = .false.
        type (T_Statistics)                         :: Statistics
        type (T_Property), pointer                  :: Next                 => null()
        type (T_Property), pointer                  :: Prev                 => null()
    end type T_Property

    type       T_Files
         character(len=Pathlength)                  :: ConstructData
         character(len=Pathlength)                  :: Results
    end type T_Files

    type      T_Atmosphere
        integer                                     :: InstanceID
        character(PathLength)                       :: ModelName
        type(T_Size2D)                              :: Size
        type(T_Size2D)                              :: WorkSize
        type(T_External)                            :: ExternalVar
        type(T_Files)                               :: Files
        type(T_Property), pointer                   :: FirstAtmosphereProp       => null()
        type(T_Property), pointer                   :: LastAtmosphereProp        => null()
        type(T_OutPut)                              :: OutPut
        type(T_Time     )                           :: BeginTime
        type(T_Time     )                           :: EndTime
        type(T_Time     )                           :: ActualTime
        type(T_Time     )                           :: NextCompute
        type(T_Time     )                           :: LastOutPutHDF5

        real                                        :: PredictedDT              = -null_real
        real                                        :: DTForNextEvent           = -null_real
        
        integer                                     :: RadiationMethod          = 1
        integer                                     :: CloudCoverMethod
        real,    pointer, dimension(:,:  )          :: LastRadiation             => null()
        integer                                     :: LastCalculateRandomCloud = null_int
        integer                                     :: PropertiesNumber         = FillValueInt
        integer                                     :: CurrentIndex             = 2
        
        logical                                     :: PropsAddedByRain         = .false.
        logical                                     :: PropsAddedByIrri         = .false.
        
        !Instance of Module HDF5
        integer                                     :: ObjHDF5 = 0

        !Instance of Module_EnterData
        integer                                     :: ObjEnterData     = 0    !Data File - ConstructData
       
        !Instance of ModuleGridData
        integer                                     :: ObjGridData     = 0
      
        !Instance of ModuleHorizontalGrid
        integer                                     :: ObjHorizontalGrid = 0

        !Instance of ModuleTime
        integer                                     :: ObjTime           = 0

        !Instance of ModuleBoxDif
        integer                                     :: ObjBoxDif         = 0

        !Instance of ModuleTimeSerie
        integer                                     :: ObjTimeSerie      = 0

        !Collection of instances                
        type(T_Atmosphere), pointer                 :: Next              => null()

    end type T_Atmosphere

    !Global Module Variables
    type (T_Atmosphere), pointer                    :: FirstObjAtmosphere  => null()
    type (T_Atmosphere), pointer                    :: Me                  => null()


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine StartAtmosphere(ModelName,                                       &
                               AtmosphereID,                                    &
                               TimeID,                                          &
                               GridDataID,                                      &
                               HorizontalGridID,                                &
                               MappingPoints,                                   &
                               STAT)

        !Arguments--------------------------------------------------------------
        character(Len=*)                            :: ModelName
        integer                                     :: AtmosphereID
        integer                                     :: TimeID         
        integer                                     :: GridDataID     
        integer                                     :: HorizontalGridID
        integer, dimension(:, :), pointer           :: MappingPoints
        integer, optional, intent(OUT)              :: STAT  

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_
 
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mAtmosphere_)) then
            nullify (FirstObjAtmosphere)
            call RegisterModule (mAtmosphere_) 
        endif


       call Ready(AtmosphereID, ready_) 

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            !Allocates a new Instance
            call AllocateInstance 

            Me%ModelName = ModelName

            !Associates External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjGridData       = AssociateInstance (mGRIDDATA_,       GridDataID      )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)

            Me%ExternalVar%MappingPoints2D => MappingPoints
            call ReadLockExternalVar

            !Read the name file of the Atmosphere module
            call ReadPropertiesFilesName

            !Construct enter data 
            call ConstructEnterData(Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'StartAtmosphere - ModuleAtmosphere - ERR01' 

            call ConstructGlobalVariables

            if (Me%CloudCoverMethod == CloudFromRadiation) then

                allocate( Me%LastRadiation (    Me%WorkSize%ILB:Me%WorkSize%IUB,    &
                                                Me%WorkSize%JLB:Me%WorkSize%JUB)    )  
                !Me%LastRadiation = null_real
                Me%LastRadiation = 0.

            endif

            !Constructs the property list 
            call ConstructPropertyList

            call ConstructTimeSerie

            call ConstructOutput

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartAtmosphere - ModuleAtmosphere - ERR02' 

            ! By default a output file is always open in the construction phase
            if (Me%OutPut%True) call OpenHDF5OutPutFile

            call null_time(Me%ActualTime)

            nullify(Me%ExternalVar%MappingPoints2D)
            call ReadUnlockExternalVar


            !Returns ID
            AtmosphereID = Me%InstanceID

            STAT_ = SUCCESS_

        else  cd0

            stop 'StartAtmosphere - ModuleAtmosphere - ERR03'

        end if cd0


        if (present(STAT)) STAT = STAT_


    end subroutine StartAtmosphere


    !--------------------------------------------------------------------------

    
    subroutine AllocateInstance 

        !Local-----------------------------------------------------------------
        type (T_Atmosphere), pointer           :: NewAtmosphere          => null()
        type (T_Atmosphere), pointer           :: PreviousAtmosphere     => null()


        !Allocates new instance
        allocate (NewAtmosphere)
        nullify  (NewAtmosphere%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjAtmosphere)) then
            FirstObjAtmosphere      => NewAtmosphere
            Me                      => NewAtmosphere
        else
            PreviousAtmosphere      => FirstObjAtmosphere
            Me                      => FirstObjAtmosphere%Next
            do while (associated(Me))
                PreviousAtmosphere  => Me
                Me                  => Me%Next
            enddo
            Me                      => NewAtmosphere
            PreviousAtmosphere%Next => NewAtmosphere
        endif

        Me%InstanceID = RegisterNewInstance (mATMOSPHERE_)

    end subroutine AllocateInstance


    !--------------------------------------------------------------------------

    subroutine ConstructGlobalVariables

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL, iflag

        !Begin------------------------------------------------------------------

        call GetComputeTimeLimits(Me%ObjTime, BeginTime = Me%BeginTime, &
                                  EndTime = Me%EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleAtmosphere - ERR01'
        
        !Actualize the time
        Me%ActualTime  = Me%BeginTime
        Me%NextCompute = Me%ActualTime

        ! Sets the last output equal to zero 
        call SetDate(Me%LastOutPutHDF5, 0, 0, 0, 0, 0, 0)

        call GetHorizontalGridSize(Me%ObjHorizontalGrid, Size = Me%Size, &
                                   WorkSize = Me%WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleAtmosphere - ERR02'


        call GetData(Me%RadiationMethod,                                                &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RADIATION_METHOD',                                 &
                     ClientModule = 'ModuleAtmosphere',                                 &
                     Default      = Radiation_Mohid,                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleAtmosphere - ERR03'


        call GetData(Me%CloudCoverMethod,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'CLOUD_COVER_METHOD',                               &
                     ClientModule = 'ModuleAtmosphere',                                 &
                     Default      = CloudFromRandom,                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleAtmosphere - ERR04'

    end subroutine ConstructGlobalVariables

    !--------------------------------------------------------------------------
    !Read the name of the files need to construct and modify
    ! the Atmosphere properties 

    subroutine ReadPropertiesFilesName

        !External--------------------------------------------------------------
        character(len = StringLength)           :: Message
        integer                                 :: STAT_CALL

        !Begin------------------------------------------------------------------

        !Opens the Atmosphere data file 
        ! ASCII file used to construct new properties
        Message  ='Atmosphere Data Properties.'
        Message  = trim(Message)

        call ReadFileName('SURF_DAT', Me%Files%ConstructData, Message = Message, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadPropertiesFilesName - ModuleAtmosphere - ERR01'

        ! ---> File in HDF format where is written instant fields of Atmosphere properties
        Message   ='Instant fields of Atmosphere properties in HDF format.'
        Message   = trim(Message)

        call ReadFileName('SURF_HDF', Me%Files%Results, Message = Message, &
                           TIME_END = Me%EndTime, Extension = 'sur', STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadPropertiesFilesName - ModuleAtmosphere - ERR02'


    end subroutine ReadPropertiesFilesName

    !--------------------------------------------------------------------------


    subroutine OpenHDF5OutPutFile

        !Local-----------------------------------------------------------------
        real, pointer, dimension(:, :)              :: GridData
        integer                                     :: STAT_CALL
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Bounds
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 


        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%Results)//"5", &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAtmosphere - ERR01'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5,       &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAtmosphere - ERR02'

        !Gets a pointer to GridData
        call GetGridData      (Me%ObjGridData, GridData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAtmosphere - ERR03'

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,        &
                              WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAtmosphere - ERR05'

        !Writes the GridData
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",            &
                              Array2D = GridData,                                      &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAtmosphere - ERR06'

        !Writes the WaterPoints
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "MappingPoints2D", "-",         &
                              Array2D = Me%ExternalVar%MappingPoints2D,            &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAtmosphere - ERR07'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAtmosphere - ERR08'

        !Ungets the GridData
        call UngetGridData (Me%ObjGridData, GridData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleAtmosphere - ERR09'


    end subroutine OpenHDF5OutPutFile

    !--------------------------------------------------------------------------


    subroutine ConstructTimeSerie

        !External--------------------------------------------------------------
        integer                                             :: iflag, STAT_CALL
         
        !Local-----------------------------------------------------------------
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK

        integer                                             :: dn, Id, Jd, TimeSerieNumber
        type(T_Property), pointer                           :: PropertyX
        integer                                             :: nProperties
        character(len=PathLength)                           :: TimeSerieLocationFile
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Begin-----------------------------------------------------------------


        !First checks out how many properties will have time series
        PropertyX   => Me%FirstAtmosphereProp
        nProperties =  0
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) nProperties = nProperties + 1
            PropertyX=>PropertyX%Next
        enddo

        if (nProperties > 0) then

            !Allocates PropertyList
            allocate(PropertyList(nProperties), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAtmosphere - ERR10'

            !Fills up PropertyList
            PropertyX   => Me%FirstAtmosphereProp
            nProperties =  0
            do while (associated(PropertyX))
                if (PropertyX%TimeSerie) then
                    nProperties = nProperties + 1
                    PropertyList(nProperties) = trim(adjustl(PropertyX%ID%name))
                endif
                PropertyX=>PropertyX%Next
            enddo


            call GetData(TimeSerieLocationFile,                                             &
                         Me%ObjEnterData, iflag,                                            &
                         SearchType   = FromFile,                                           &
                         keyword      = 'TIME_SERIE_LOCATION',                              &
                         ClientModule = 'ModuleWaterProperties',                            &
                         Default      = Me%Files%ConstructData,                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                    &
                stop 'Construct_Time_Serie - ModuleAtmosphere - ERR20' 


            !Constructs TimeSerie
            call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                            &
                                trim(TimeSerieLocationFile),                            &
                                PropertyList, "srs",                                    &
                                WaterPoints2D = Me%ExternalVar%MappingPoints2D,         &
                                ModelName     = Me%ModelName,                           &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAtmosphere - ERR30'

            !Deallocates PropertyList
            deallocate(PropertyList, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAtmosphere - ERR40'

            !Corrects if necessary the cell of the time serie based in the time serie coordinates
            call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAtmosphere - ERR50'

            do dn = 1, TimeSerieNumber

                call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                          CoordX   = CoordX,                                &
                                          CoordY   = CoordY,                                & 
                                          CoordON  = CoordON,                               &
                                          STAT     = STAT_CALL)
                if (CoordON) then
                    call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAtmosphere - ERR60'

                    if (Id < 0 .or. Jd < 0) then
                
                        call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAtmosphere - ERR70'

                        if (IgnoreOK) then
                            cycle
                        else
                            stop 'ConstructTimeSerie - ModuleAtmosphere - ERR80'
                        endif

                    endif

                    call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleAtmosphere - ERR90'
                endif

            enddo
        endif

    end subroutine ConstructTimeSerie

    !--------------------------------------------------------------------------

    subroutine ConstructPropertyList

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: NewProperty      => null()
        type (T_Property), pointer                  :: PropertyX        => null()
        type (T_Property), pointer                  :: PropertyY        => null()
        integer                                     :: ClientNumber
        integer                                     :: i, j, STAT_CALL
        logical                                     :: BlockFound

        !----------------------------------------------------------------------

        ! Initialize the Atmosphere properties number   
        Me%PropertiesNumber = 0

        ! Initialize the Atmosphere properties list   
        nullify (Me%FirstAtmosphereProp)
        nullify (Me%LastAtmosphereProp)

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,          &
                                        block_begin, block_end, BlockFound,     &
                                        STAT = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_      ) then    
cd2 :           if (BlockFound) then                                                  
                    ! Construct a New Property 
                    Call ConstructProperty(NewProperty)

                    ! Add new Property to the Atmosphere List 
                    Call Add_Property(NewProperty)
                else
                    
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ConstructPropertyList - ModuleAtmosphere - ERR01'

                    exit do1    !No more blocks

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ConstructPropertyList - ModuleAtmosphere - ERR02'
            end if cd1
        end do do1


        !Verifies Wind conistence
        call SearchProperty(PropertyX, WindVelocityX_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            call SearchProperty(PropertyY, WindVelocityY_, .false., STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then 
                if ((      PropertyX%ID%SolutionFromFile .and. .not. PropertyY%ID%SolutionFromFile) .or. &
                    (.not. PropertyX%ID%SolutionFromFile .and.       PropertyY%ID%SolutionFromFile)) then
                    write (*,*) 'wind velocity X must be given in the same way as wind velocity Y'
                    stop 'ConstructPropertyList - ModuleAtmosphere - ERR99'
                endif
            endif
        endif

        !Rotates Vectores
        call RotateAtmosphereVectorFields(Constructing = .true.)

        !Checks if Relative Humidity is between 0 and 1
        call SearchProperty(PropertyX, RelativeHumidity_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExternalVar%MappingPoints2D(i, j) == 1) then
                    if (PropertyX%Field(i, j) > 1.0) then
                        write(*,*)'Relative Humidity must be given between 0 and 1'
                        stop 'ConstructPropertyList - ModuleAtmosphere - ERR98'
                    endif
                endif    
            enddo
            enddo
        endif


        !If Solar radiation exists, add ATMTransmitivity
        call SearchProperty(PropertyX, SolarRadiation_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then

            allocate (NewProperty, STAT = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyList - ModuleAtmosphere - ERR99'

            !nullify(NewProperty%Field    )
            !nullify(NewProperty%FieldGrid)
            !nullify(NewProperty%Next     )
            !nullify(NewProperty%Prev     )

            allocate(NewProperty%Field(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            NewProperty%Field(:,:) = null_real

            NewProperty%FieldGrid => NewProperty%Field

            NewProperty%ID%IDnumber = AtmTransmitivity_
            NewProperty%ID%Name     = GetPropertyName (AtmTransmitivity_)

            NewProperty%Constant    = .false.

            call Add_Property(NewProperty)

        endif

        call CheckForObsoleteNames

    end subroutine ConstructPropertyList

    !--------------------------------------------------------------------------

    subroutine CheckForObsoleteNames

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX        => null()
        logical                                     :: UsingObseletePropertyName

        !Begin-----------------------------------------------------------------

        UsingObseletePropertyName = .false.

        PropertyX => Me%FirstAtmosphereProp

        do while (associated(PropertyX)) 
            if    (PropertyX%ID%IDNumber==WindAngle_    ) then
                
                UsingObseletePropertyName = .true.
                
                write(*,*)
                write(*,*)"You are trying to use property wind angle in ModuleAtmosphere"
                write(*,*)"This property name is obsolete and is now called wind direction" 
                write(*,*)"Please update your Atmosphere input data file"
                write(*,*)

            elseif(PropertyX%ID%IDNumber==WindModulos_  ) then
                
                UsingObseletePropertyName = .true.

                write(*,*)
                write(*,*)"You are trying to use property wind modulos in ModuleAtmosphere"
                write(*,*)"This property name is obsolete and is now called wind modulus" 
                write(*,*)"Please update your Atmosphere input data file"
                write(*,*)

            end if
            PropertyX => PropertyX%Next
        end do 

        if(UsingObseletePropertyName)then
            write(*,*)'CheckForObsoleteNames - ModuleAtmosphere - ERR01'
            stop
        endif

        nullify(PropertyX)

    end subroutine CheckForObsoleteNames

    !--------------------------------------------------------------------------

    subroutine RotateAtmosphereVectorFields(Constructing)
        !Arguments-------------------------------------------------------------
        logical                                     :: Constructing

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX        => null()
        type (T_Property), pointer                  :: PropertyY        => null()
        integer                                     :: STAT_CALL


        !----------------------------------------------------------------------


        call SearchProperty(PropertyX, WindVelocityX_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            call SearchProperty(PropertyY, WindVelocityY_, .false., STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then 

                if (Constructing) then 
                    nullify (PropertyX%FieldGrid) 
                    allocate(PropertyX%FieldGrid(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

                    PropertyX%FieldGrid(:,:) = null_real

                    nullify (PropertyY%FieldGrid) 
                    allocate(PropertyY%FieldGrid(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

                    PropertyY%FieldGrid(:,:) = null_real

                endif


                call RotateVectorFieldToGrid(HorizontalGridID  = Me%ObjHorizontalGrid, &
                                             VectorInX         = PropertyX%Field,                      &
                                             VectorInY         = PropertyY%Field,                      &
                                             VectorOutX        = PropertyX%FieldGrid,                  &
                                             VectorOutY        = PropertyY%FieldGrid,                  &   
                                             WaterPoints2D     = Me%ExternalVar%MappingPoints2D,       &
                                             RotateX           = .true.,                               &
                                             RotateY           = .true.,                               &
                                             STAT              = STAT_CALL)
            endif
        endif


    end subroutine RotateAtmosphereVectorFields

    
    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct a new property.           
    subroutine ConstructProperty(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer       :: NewProperty

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL

        !----------------------------------------------------------------------
             

        !Allocates new property
        allocate (NewProperty, STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ConstructProperty - ModuleAtmosphere - ERR01'

        nullify(NewProperty%Field    )
        nullify(NewProperty%FieldGrid)
        nullify(NewProperty%Next     )
        nullify(NewProperty%Prev     )

        !Construct property ID
        call ConstructPropertyID        (NewProperty%ID, Me%ObjEnterData, FromBlock)

        !Construct property values
        call Construct_PropertyValues   (NewProperty)

        !Defines the property output
        call Construct_PropertyOutPut   (NewProperty)

        !Constructs Statistics
        call ConstructSurfStatistics    (NewProperty) 

        !----------------------------------------------------------------------

    end subroutine ConstructProperty

    !--------------------------------------------------------------------------
    
    !This subroutine reads all the information needed to construct the property values       
    ! in the domain and in the boundaries            
    subroutine Construct_PropertyValues (NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),   pointer                 :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: iflag
        integer                                     :: SizeILB, SizeIUB
        integer                                     :: SizeJLB, SizeJUB 
        integer                                     :: i, j
        
        !----------------------------------------------------------------------

        SizeILB = Me%Size%ILB
        SizeIUB = Me%Size%IUB
        SizeJLB = Me%Size%JLB
        SizeJUB = Me%Size%JUB

        !Fills Matrix
        allocate (NewProperty%Field (SizeILB:SizeIUB, SizeJLB:SizeJUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleAtmosphere - Construct_PropertyValues - ERR00'

        NewProperty%Field(:,:) = null_real

        NewProperty%FieldGrid => NewProperty%Field

        !Properties added in a specific time interval (normally a short interval)
        !This is usefull for near instantaneous events (ex.:high precipitation in
        !very small time period). In this option instead of interpolation one 
        !calculates the exact amount of property in a time period. This requires
        !variableDT.
        call GetData(NewProperty%NoInterpolateValueInTime,                      &
                     Me%ObjEnterData, iflag,                                    &
                     Default      = .false.,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      ='NO_INTERPOLATION',                          &
                     ClientModule = 'ModuleAtmosphere',                         &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Construct_PropertyValues - ModuleAtmosphere - ERR00'

        call ConstructFillMatrix(PropertyID         = NewProperty%ID,                    &
                                 EnterDataID        = Me%ObjEnterData,                   &
                                 TimeID             = Me%ObjTime,                        &
                                 HorizontalGridID   = Me%ObjHorizontalGrid,              &
                                 ExtractType        = FromBlock,                         &
                                 PointsToFill2D     = Me%ExternalVar%MappingPoints2D,    &
                                 Matrix2D           = NewProperty%Field,                 &
                                 TypeZUV            = TypeZ_,                            &
                                 STAT               = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleAtmosphere - ERR01'

        call GetIfMatrixRemainsConstant(FillMatrixID    = NewProperty%ID%ObjFillMatrix,     &
                                        RemainsConstant = NewProperty%Constant,             &
                                        STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                          &
            stop 'Construct_PropertyValues - ModuleAtmosphere - ERR02'

        if (.not. NewProperty%ID%SolutionFromFile) then
            call KillFillMatrix (NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleAtmosphere - ERR05'
        endif

        !By default property don't have random component
        NewProperty%HasRandomComponent = .false.

        call GetData(NewProperty%RandomComponent,               &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'RANDOM_COMPONENT',         &
                     ClientModule = 'ModuleAtmosphere',         &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Construct_PropertyValues - ModuleAtmosphere - ERR11'
        if (iflag == 1) then
            NewProperty%HasRandomComponent = .true.
            NewProperty%RandomValue        = 0.
        endif

        if(NewProperty%HasRandomComponent .and. NewProperty%Constant)then
            write(*,*)
            write(*,*)'WARNING - Atmosphere property has random component and is defined as constant'
            write(*,*)'Property name : ', trim(NewProperty%ID%Name)
            write(*,*)
        end if


        !By default property aren't added by irrigation
        NewProperty%PropAddedByIrri = .false.
        call GetData(NewProperty%PropAddedByIrri,                               &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      ='IRRIGATION',                                &
                     ClientModule = 'ModuleAtmosphere',                         &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Construct_PropertyValues - ModuleAtmosphere - ERR12'
        if (NewProperty%PropAddedByIrri) then
            Me%PropsAddedByIrri = .true.
        endif

        !By default property aren't added by precipitation
        NewProperty%PropAddedByRain = .false.
        call GetData(NewProperty%PropAddedByRain,                               &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      ='PRECIPITATION',                             &
                     ClientModule = 'ModuleAtmosphere',                         &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Construct_PropertyValues - ModuleAtmosphere - ERR13'
        if (NewProperty%PropAddedByRain) then
            Me%PropsAddedByRain = .true.
        endif
        
        if (NewProperty%ID%IDNumber == WindDirection_) then


            !A rotation of the wind direction is done, so that
            !  0 - Wind from the North
            ! 90 - Wind from the West
            !180 - Wind from the South
            !270 - Wind from the East
            do j = Me%Size%JLB, Me%Size%JUB 
            do i = Me%Size%ILB, Me%Size%IUB 
                if (Me%ExternalVar%MappingPoints2D(i, j) == WaterPoint) then
                    NewProperty%Field(i, j) = 270. - NewProperty%Field(i, j)
                else
                    NewProperty%Field(i, j) = FillValueReal
                endif
            enddo
            enddo


        endif

        !----------------------------------------------------------------------

    end subroutine Construct_PropertyValues

    !--------------------------------------------------------------------------
    
    subroutine ConstructOutPut
        
        !External-----------------------------------------------------------------
        integer :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetOutPutTime(Me%ObjEnterData,                                     &
                           CurrentTime   = Me%ActualTime,                       &
                           EndTime       = Me%EndTime,                          &
                           keyword       = 'OUTPUT_TIME',                       &
                           SearchType    = FromFile,                            &
                           OutPutsTime   = Me%OutPut%OutTime,                   &
                           OutPutsOn     = Me%OutPut%True,                      &
                           OutPutsNumber = Me%OutPut%Number,                    &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOutPut - ModuleAtmosphere - ERR01' 

        if (Me%OutPut%True) Me%OutPut%NextHDF5 = 1

    end subroutine ConstructOutPut

    !--------------------------------------------------------------------------

    subroutine Construct_PropertyOutPut(NewProperty)

        !Arguments------------------------------------------------------------
        type(T_property),           pointer     :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: iflag

        !Begin----------------------------------------------------------------

        call GetData(NewProperty%OutputHDF,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'OUTPUT_HDF',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleAtmosphere',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleAtmosphere - ERR01'
           
        call GetData(NewProperty%TimeSerie,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'TIME_SERIE',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleAtmosphere',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleAtmosphere - ERR02'


        call GetData(NewProperty%BoxTimeSerie,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'BOX_TIME_SERIE',                                 &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleAtmosphere',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleAtmosphere - ERR03'

    end subroutine Construct_PropertyOutPut

    !--------------------------------------------------------------------------

    subroutine ConstructSurfStatistics(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),   pointer     :: NewProperty

        !Local-----------------------------------------------------------------
        integer                         :: STAT_CALL
        integer                         :: iflag
        integer                         :: ILB, IUB, JLB, JUB
        integer                         :: WILB, WIUB, WJLB, WJUB
       
        !Begin-----------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        WILB = Me%WorkSize%ILB 
        WIUB = Me%WorkSize%IUB 
        WJLB = Me%WorkSize%JLB 
        WJUB = Me%WorkSize%JUB 

        !<BeginKeyword>
            !Keyword          : STATISTICS
            !<BeginDescription>       
               ! 
               ! Checks out if the user pretends the statistics of this property
               ! 
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : DISPQUAL
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%Statistics%ON,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'STATISTICS',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleAtmosphere',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSurfStatistics - Atmosphere - ERR01'

        
cd2:    if (NewProperty%Statistics%ON) then

            !<BeginKeyword>
                !Keyword          : STATISTICS_FILE
                !<BeginDescription>       
                   ! 
                   ! The statistics definition file of this property
                   ! 
                !<EndDescription>
                !Type             : Character
                !Default          : Do not have
                !File keyword     : DISPQUAL
                !Multiple Options : Do not have
                !Search Type      : FromBlock
                !Begin Block      : <beginproperty>
                !End Block        : <endproperty>
            !<EndKeyword>
            call GetData(NewProperty%Statistics%File,                                   &
                         Me%ObjEnterData, iflag,                                        &
                         Keyword    = 'STATISTICS_FILE',                                &
                         SearchType = FromBlock,                                        &
                         ClientModule = 'ModuleAtmosphere',                             &
                         STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_ .or. iflag /= 1)                                  &
                stop 'ConstructSurfStatistics - Atmosphere - ERR02'

            call ConstructStatistic (StatisticID      = NewProperty%Statistics%ID,   &
                                     ObjTime          = Me%ObjTime,                  &
                                     ObjHDF5          = Me%ObjHDF5,                  &
                                     Size             = T_Size3D(ILB, IUB, JLB, JUB,0,0),       &
                                     WorkSize         = T_Size3D(WILB, WIUB, WJLB, WJUB,0,0),   &
                                     DataFile         = NewProperty%Statistics%File, &
                                     Name             = NewProperty%ID%name,         &
                                     STAT             = STAT_CALL)                                 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSurfStatistics - Atmosphere - ERR03'

        endif cd2


    end subroutine ConstructSurfStatistics

    !--------------------------------------------------------------------------
    ! This subroutine adds a new property to the Water Property List  

    subroutine Add_Property(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),   pointer     :: NewProperty

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstAtmosphereProp)) then
            Me%PropertiesNumber = 1
            Me%FirstAtmosphereProp    => NewProperty
            Me%LastAtmosphereProp     => NewProperty
        else
            NewProperty%Prev                     => Me%LastAtmosphereProp
            Me%LastAtmosphereProp%Next => NewProperty
            Me%LastAtmosphereProp      => NewProperty
            Me%PropertiesNumber     = Me%PropertiesNumber + 1
        end if 

        !----------------------------------------------------------------------

    end subroutine Add_Property 


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetAtmosphereProperty(AtmosphereID, Scalar, ID, STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: AtmosphereID
        real, dimension(:,:), pointer               :: Scalar
        integer                                     :: ID
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        
        type (T_Property), pointer                  :: PropertyX    
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(AtmosphereID, ready_) 

        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            nullify(PropertyX)
            call SearchProperty(PropertyX, ID , .true., STAT = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'GetAtmosphereProperty - ModuleAtmosphere - ERR01'     

            call Read_Lock(mATMOSPHERE_, Me%InstanceID)

            Scalar  => PropertyX%FieldGrid

          
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetAtmosphereProperty

    !--------------------------------------------------------------------------

    logical function AtmospherePropertyExists(AtmosphereID, PropertyNumber)

        !Arguments--------------------------------------------------------------
        integer                                 :: AtmosphereID
        integer                                 :: PropertyNumber
     
        !External--------------------------------------------------------------
        integer                                 :: ready_          

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, STAT_CALL
        type(T_Property), pointer               :: PropertyX

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(AtmosphereID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
                        
            AtmospherePropertyExists = .false.

            call SearchProperty(PropertyX, PropertyNumber, .false., STAT_CALL)
            if (STAT_CALL == SUCCESS_) AtmospherePropertyExists = .true.

            nullify(PropertyX)

        else 
            stop 'AtmospherePropertyExists - ModuleAtmosphere - ERR01'
        end if cd1


    end function AtmospherePropertyExists

    !--------------------------------------------------------------------------

    subroutine GetAtmospherenProperties (AtmosphereID, nProperties, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: AtmosphereID
        integer                                         :: nProperties
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(AtmosphereID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            nProperties       = Me%PropertiesNumber
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetAtmospherenProperties

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetAtmospherePropertiesIDByIdx (AtmosphereID, Idx, ID,PropRain,PropIrri, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: AtmosphereID
        integer, intent(IN)                             :: Idx
        integer, intent(OUT)                            :: ID
        logical, intent(OUT)                            :: PropRain, PropIrri
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, i
        type (T_Property), pointer                      :: CurrProp

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(AtmosphereID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            CurrProp => Me%FirstAtmosphereProp
            do i = 1, idx - 1
                CurrProp => CurrProp%Next
            enddo

            ID        = CurrProp%ID%IDNumber
            PropRain  = CurrProp%PropAddedByRain
            PropIrri  = CurrProp%PropAddedByIrri
            
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetAtmospherePropertiesIDByIdx

    !---------------------------------------------------------------------------

    subroutine GetAtmosphereDTPrediction (AtmosphereID, PredictedDT, DTForNextEvent, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: AtmosphereID
        real, intent(OUT)                               :: PredictedDT, DTForNextEvent
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AtmosphereID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            PredictedDT     = Me%PredictedDT
            DTForNextEvent  = Me%DTForNextEvent

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetAtmosphereDTPrediction

    !--------------------------------------------------------------------------

    subroutine UngetAtmosphere1D(AtmosphereID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: AtmosphereID
        real, pointer, dimension(:)     :: Array
        integer, optional, intent (OUT) :: STAT
   
        !External--------------------------------------------------------------
        integer                         :: ready_   

        !Local-----------------------------------------------------------------
        integer                         :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AtmosphereID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mATMOSPHERE_, Me%InstanceID, "UngetSurface1D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetAtmosphere1D


    !--------------------------------------------------------------------------

    subroutine UngetAtmosphere2D(AtmosphereID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: AtmosphereID
        real, pointer, dimension(:,:)   :: Array
        integer, optional, intent (OUT) :: STAT
   
        !External--------------------------------------------------------------
        integer                         :: ready_   

        !Local-----------------------------------------------------------------
        integer                         :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AtmosphereID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mATMOSPHERE_, Me%InstanceID, "UngetAtmosphere2D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetAtmosphere2D

    !--------------------------------------------------------------------------

    subroutine SearchProperty(PropertyX, PropertyXIDNumber, PrintWarning, STAT)


        !Arguments-------------------------------------------------------------
        type(T_Property), optional, pointer         :: PropertyX
        integer         , optional, intent (IN)     :: PropertyXIDNumber
        logical,          optional, intent (IN)     :: PrintWarning
        integer         , optional, intent (OUT)    :: STAT

        !Local-----------------------------------------------------------------

        integer                                     :: STAT_ 
        
        !----------------------------------------------------------------------

        STAT_  = UNKNOWN_

        PropertyX => Me%FirstAtmosphereProp

do2 :   do while (associated(PropertyX)) 
if5 :       if (PropertyX%ID%IDNumber==PropertyXIDNumber) then
                exit do2 
            else
                PropertyX => PropertyX%Next                 
            end if if5
        end do do2

       !A PropertyX was found
       if (associated(PropertyX)) then
            STAT_ = SUCCESS_  
        else
            if (present(PrintWarning)) then
                if (PrintWarning) write (*,*)'Property Not Found in Module Atmosphere ', &
                                              trim(GetPropertyName(PropertyXIDNumber))
            endif
            STAT_  = NOT_FOUND_ERR_  
        end if


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SearchProperty

    
    !----------------------------------------------------------------------
    
    subroutine ReadLockExternalVar
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !GridCellArea
        call GetGridCellArea (Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleAtmosphere - ERR01'

        !Lat / Lon
        call GetGridLatitudeLongitude(Me%ObjHorizontalGrid,                     &
                                      GridLatitude  = Me%ExternalVar%Latitude,  &
                                      GridLongitude = Me%ExternalVar%Longitude, &
                                      STAT          = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleAtmosphere - ERR02'



    end subroutine ReadLockExternalVar

    !----------------------------------------------------------------------

    subroutine ReadUnlockExternalVar
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !GridCellArea
        call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleAtmosphere - ERR01'

        !Latitude
        call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%Latitude, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleAtmosphere - ERR02'

        !Longitude
        call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%Longitude, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleAtmosphere - ERR03'


    end subroutine ReadUnlockExternalVar


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine ModifyAtmosphere(AtmosphereID, MappingPoints, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: AtmosphereID
        integer, dimension(:, :), pointer           :: MappingPoints
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                   :: PropertyX    => null()
        type(T_Property), pointer                   :: PropertyY    => null()
        integer                                     :: STAT_, ready_

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin------------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(AtmosphereID, ready_)

cd0:    if (ready_ .EQ. IDLE_ERR_) then

            !Gets Current Time
            call GetComputeCurrentTime(Me%ObjTime, Me%ActualTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyAtmosphere - ModuleAtmosphere - ERR01'

            Me%ExternalVar%MappingPoints2D => MappingPoints
            call ReadLockExternalVar

            call SearchProperty(PropertyX, WindModulus_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_) call ModifyWindModulus (PropertyX)

            call SearchProperty(PropertyX, WindDirection_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_) call ModifyWindDirection (PropertyX)

            call SearchProperty(PropertyX, WindVelocityX_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then
                call SearchProperty(PropertyY, WindVelocityY_, STAT = STAT_CALL) 
                if (STAT_CALL == SUCCESS_) call ModifyWindVelocity (PropertyX, PropertyY)
            end if
            
            call SearchProperty(PropertyX, AirTemperature_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_) call ModifyAirTemperature (PropertyX)

            call SearchProperty(PropertyX, RelativeHumidity_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_) call ModifyRelativeHumidity (PropertyX)

            call SearchProperty(PropertyX, SunHours_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_) call ModifySunHours (PropertyX)

            call SearchProperty(PropertyX, CloudCover_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_) call ModifyCloudCover (PropertyX)

            call SearchProperty(PropertyX, Irrigation_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_) call ModifyIrrigation (PropertyX)

            call SearchProperty(PropertyX, Precipitation_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_) call ModifyPrecipitation (PropertyX)

            call SearchProperty(PropertyX, SolarRadiation_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_) call ModifySolarRadiation (PropertyX)

            call SearchProperty(PropertyX, AtmosphericPressure_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_) call ModifyAtmosphericPressure (PropertyX)

            call SearchProperty(PropertyX, CO2AtmosphericPressure_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_) call ModifyCO2AtmosphericPressure (PropertyX)

            call SearchProperty(PropertyX, O2AtmosphericPressure_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_) call ModifyO2AtmosphericPressure (PropertyX)


            call ModifyRandom
            
            if (Me%PropsAddedByIrri) then
                call ModifyPropByIrri
            endif
            
            if (Me%PropsAddedByRain) then
                call ModifyPropByRain
            endif

            call ModifyOutPut

            call RotateAtmosphereVectorFields(Constructing = .false.)

            nullify (Me%ExternalVar%MappingPoints2D)
            call ReadUnlockExternalVar


            STAT_ = SUCCESS_
        else   cd0
            STAT_ = ready_
        end if cd0

        if (present(STAT)) STAT = STAT_


    end subroutine ModifyAtmosphere

    !----------------------------------------------------------------------

    subroutine ModifyRandom

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropertyX

        !Local-----------------------------------------------------------------
        real                                        :: RandomValue
        integer                                     :: IUB, ILB, JUB, JLB, i, j
        integer                                        :: CHUNK
        
        !Begin------------------------------------------------------------------------

        !Begin - Shorten variables name 
        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB
        !End - Shorten variables name 

        PropertyX => Me%FirstAtmosphereProp

do1 :   do while (associated(PropertyX)) 

            !Add random component
            if (PropertyX%HasRandomComponent)then
        
                call random_number(RandomValue)
        
                RandomValue = (RandomValue - 0.5) * PropertyX%RandomComponent

                if (MonitorPerformance) then
                    call StartWatch ("ModuleAtmosphere", "ModifyRandom")
                endif

                CHUNK = CHUNK_J(JLB, JUB)
                !$OMP PARALLEL PRIVATE(i,j)
                !Substract previous random field    
                !$OMP DO SCHEDULE(STATIC,CHUNK)
                do j = JLB, JUB
                do i = ILB, IUB
                    PropertyX%Field(i, j) = PropertyX%Field(i, j) - PropertyX%RandomValue
                enddo
                enddo
                !$OMP END DO
                
                !Add new random value
                !$OMP DO SCHEDULE(STATIC,CHUNK)
                do j = JLB, JUB
                do i = ILB, IUB
                    PropertyX%Field(i, j) = PropertyX%Field(i, j) + RandomValue
                enddo
                enddo
                !$OMP END DO
                !$OMP END PARALLEL

                if (MonitorPerformance) then
                    call StopWatch ("ModuleAtmosphere", "ModifyRandom")
                endif

                !Stores random value
                PropertyX%RandomValue  = RandomValue

            endif
    
            PropertyX => PropertyX%Next                 
        end do do1


    end subroutine ModifyRandom

    !----------------------------------------------------------------------


    subroutine ModifyOutPut

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropertyX

        !Begin------------------------------------------------------------------------

        !Output HDF
        if (Me%OutPut%True) then
            call OutPutResultsHDF5 
        endif


        !Output TimeSerie / Statistics
        PropertyX => Me%FirstAtmosphereProp

do2 :   do while (associated(PropertyX)) 

   
            !Output TimeSerie
            call OutPut_TimeSeries(PropertyX)

            !OutPut Statistics
            if (PropertyX%Statistics%ON) then

                call OutPut_Statistics(PropertyX%Field, PropertyX%Statistics%ID)

            endif

            PropertyX%FirstActualization = .false.

            PropertyX => PropertyX%Next

        enddo  do2

    end subroutine ModifyOutPut

    !----------------------------------------------------------------------

    subroutine ModifySolarRadiation(PropSolarRadiation)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropSolarRadiation

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: i,j
        integer                                        :: CHUNK

        !Begin-----------------------------------------------------------------
       
        if (PropSolarRadiation%ID%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = PropSolarRadiation%ID%ObjFillMatrix, &
                                   Matrix2D       = PropSolarRadiation%Field,            &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,      &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifySolarRadiation - ModuleAtmosphere - ERR01'

             !call ClimatologicSolarRadiation (PropSolarRadiation%Field)

            !Save the last radiation of the day
            if (Me%CloudCoverMethod == CloudFromRadiation ) then
           
                   if (MonitorPerformance) then
                    call StartWatch ("ModuleAtmosphere", "ModifySolarRadiation")
                endif
           
                CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
                !$OMP PARALLEL PRIVATE(i,j)
                !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExternalVar%MappingPoints2D(i, j) == 1) then
                        if (PropSolarRadiation%Field(i,j) .gt. 0.0) then                        
                            Me%LastRadiation(i,j) = PropSolarRadiation%Field(i,j)
                        endif
                    endif
                enddo
                enddo
                !$OMP END DO
                !$OMP END PARALLEL

                   if (MonitorPerformance) then
                    call StopWatch ("ModuleAtmosphere", "ModifySolarRadiation")
                endif

            endif


        elseif (.not. PropSolarRadiation%Constant) then

            select case(Me%RadiationMethod)
                
                case(Radiation_MOHID)
                        
                    call ClimatologicSolarRadiation (PropSolarRadiation)

                case(Radiation_CEQUALW2)

                    call CEQUALW2SolarRadiation     (PropSolarRadiation)

            end select

        endif

    end subroutine ModifySolarRadiation

    !--------------------------------------------------------------------------
    
    subroutine ModifyAtmosphericPressure(PropAtmosphericPressure)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropAtmosphericPressure

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-----------------------------------------------------------------

        if (PropAtmosphericPressure%ID%SolutionFromFile .and. .not. &
            PropAtmosphericPressure%Constant) then

            call ModifyFillMatrix (FillMatrixID   = PropAtmosphericPressure%ID%ObjFillMatrix,&
                                   Matrix2D       = PropAtmosphericPressure%Field,           &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,          &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyAtmosphericPressure - ModuleAtmosphere - ERR01'

        endif

    end subroutine ModifyAtmosphericPressure

    
    !--------------------------------------------------------------------------


    subroutine ModifyCO2AtmosphericPressure(PropCO2AtmosphericPressure)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropCO2AtmosphericPressure

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-----------------------------------------------------------------

        if (PropCO2AtmosphericPressure%ID%SolutionFromFile .and. .not. &
            PropCO2AtmosphericPressure%Constant) then

            call ModifyFillMatrix (FillMatrixID   = PropCO2AtmosphericPressure%ID%ObjFillMatrix,&
                                   Matrix2D       = PropCO2AtmosphericPressure%Field,           &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,             &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyCO2AtmosphericPressure - ModuleAtmosphere - ERR01'

        endif

    end subroutine ModifyCO2AtmosphericPressure

    
    !--------------------------------------------------------------------------


    subroutine ModifyO2AtmosphericPressure(PropO2AtmosphericPressure)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropO2AtmosphericPressure

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-----------------------------------------------------------------

        if (PropO2AtmosphericPressure%ID%SolutionFromFile .and. .not. &
            PropO2AtmosphericPressure%Constant) then

            call ModifyFillMatrix (FillMatrixID   = PropO2AtmosphericPressure%ID%ObjFillMatrix,&
                                   Matrix2D       = PropO2AtmosphericPressure%Field,           &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,            &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyO2AtmosphericPressure - ModuleAtmosphere - ERR01'

        endif

    end subroutine ModifyO2AtmosphericPressure

    
    !--------------------------------------------------------------------------

    
    subroutine ModifyWindVelocity(PropWindVelocityX, PropWindVelocityY)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropWindVelocityX, PropWindVelocityY

        !Begin-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        if (PropWindVelocityX%ID%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = PropWindVelocityX%ID%ObjFillMatrix,  &
                                   Matrix2D       = PropWindVelocityX%Field,             &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,      &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifySolarRadiation - ModuleAtmosphere - ERR01'

        end if

        if (PropWindVelocityY%ID%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = PropWindVelocityY%ID%ObjFillMatrix,  &
                                   Matrix2D       = PropWindVelocityY%Field,             &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,      &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifySolarRadiation - ModuleAtmosphere - ERR01'

        end if


        !Computes Wind Velocity from Modulus and Direction
        if (.not. PropWindVelocityX%Constant            .and. &
            .not. PropWindVelocityY%Constant            .and. &
            .not. PropWindVelocityX%ID%SolutionFromFile .and. &
            .not. PropWindVelocityY%ID%SolutionFromFile ) then

            call ComputeWindVelocity (PropWindVelocityX, PropWindVelocityY)

        endif

    end subroutine ModifyWindVelocity

    !--------------------------------------------------------------------------

    subroutine ComputeWindVelocity (PropWindVelocityX, PropWindVelocityY)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropWindVelocityX
        type(T_Property), pointer                   :: PropWindVelocityY

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                   :: PropWindModulus
        type(T_Property), pointer                   :: PropWindDirection
        integer                                     :: STAT_CALL

        !Begin-------------------------------------------------------------------
        
        !Updates and gets pointer to PropWindModulus
        call SearchProperty(PropWindModulus, WindModulus_ , .true., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeWindVelocity - ModuleAtmosphere - ERR01'
            

        !Updates and gets pointer to PropWindModulus
        call SearchProperty(PropWindDirection, WindDirection_, .true., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeWindVelocity - ModuleAtmosphere - ERR02'

        PropWindVelocityX%Field = PropWindModulus%Field * cos(PropWindDirection%Field * PI / 180.)
        PropWindVelocityY%Field = PropWindModulus%Field * sin(PropWindDirection%Field * PI / 180.)

    end subroutine ComputeWindVelocity

    !--------------------------------------------------------------------------

    subroutine ModifyPrecipitation(PropPrecipitation)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropPrecipitation

        !Local-----------------------------------------------------------------
        real                                        :: ConversionFactor
        integer                                     :: i, j, STAT_CALL
        integer                                        :: CHUNK

        !Begin-----------------------------------------------------------------

        if (PropPrecipitation%ID%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = PropPrecipitation%ID%ObjFillMatrix,  &
                                   Matrix2D       = PropPrecipitation%Field,             &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,      &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyPrecipitation - ModuleAtmosphere - ERR01'

            call GetFillMatrixDTPrediction (PropPrecipitation%ID%ObjFillMatrix, Me%PredictedDT,    &
                                            Me%DTForNextEvent, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyPrecipitation - ModuleAtmosphere - ERR02'

        endif

        if (PropPrecipitation%FirstActualization .or. .not. PropPrecipitation%Constant) then

            if(PropPrecipitation%NoInterpolateValueInTime) then
                if (PropPrecipitation%ID%Units /= 'mm') then
                    write(*,*)'Invalid Precipitation Units for accumulated rain'
                    write(*,*)'Use mm'
                    stop 'ModifyPrecipitation - ModuleAtmosphere - ERR01b'
                endif
            endif

            if (trim(adjustl(PropPrecipitation%ID%Units)) /= 'm3/s') then

                select case (PropPrecipitation%ID%Units)

                case ('mm/day')

                    ConversionFactor = 1. / 86400.        !In mm/s

                case ('mm/hour')

                    ConversionFactor = 1. / 3600.         !In mm/s
                
                case ('mm')

                    ConversionFactor = 1.                 !In mm

                case default

                    write(*,*)'Invalid Precipitation Units'
                    write(*,*)'Use m3/s, mm/day, mm/hour or mm'
                    stop 'ModifyPrecipitation - ModuleAtmosphere - ERR02'


                end select

                   if (MonitorPerformance) then
                    call StartWatch ("ModuleAtmosphere", "ModifyPrecipitation")
                endif

                CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
                !$OMP PARALLEL PRIVATE(i,j)
                !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExternalVar%MappingPoints2D(i, j) == 1) then
                        PropPrecipitation%Field(i, j) = PropPrecipitation%Field(i, j)      * &
                                                        (Me%ExternalVar%GridCellArea(i, j) * &
                                                         ConversionFactor) / 1000. !In m3/s
                    endif
                enddo
                enddo
                !$OMP END DO
                !$OMP END PARALLEL
        
                if (MonitorPerformance) then
                    call StopWatch ("ModuleAtmosphere", "ModifyPrecipitation")
                endif
        
            endif

        endif

    end subroutine ModifyPrecipitation

    !--------------------------------------------------------------------------

    subroutine ModifyPropByRain

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropertyX

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        PropertyX => Me%FirstAtmosphereProp

do1 :   do while (associated(PropertyX)) 

            if (PropertyX%PropAddedByRain) then

                if (PropertyX%ID%SolutionFromFile) then

                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,     &
                                           Matrix2D       = PropertyX%Field,                &
                                           PointsToFill2D = Me%ExternalVar%MappingPoints2D, &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModifyPropByRain - ModuleAtmosphere - ERR01'

                endif
                
            endif

            PropertyX => PropertyX%Next

        end do do1

    end subroutine ModifyPropByRain

    !--------------------------------------------------------------------------

    subroutine ModifySunHours (PropSunHours)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropSunHours

        !Begin-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        if (PropSunHours%ID%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = PropSunHours%ID%ObjFillMatrix,       &
                                   Matrix2D       = PropSunHours%Field,                  &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,      &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifySunHours - ModuleAtmosphere - ERR01'

        endif

    end subroutine ModifySunHours

    !------------------------------------------------------------------------------
    

    subroutine ModifyAirTemperature (PropAirTemperature)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropAirTemperature

        !Begin-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        if (PropAirTemperature%ID%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = PropAirTemperature%ID%ObjFillMatrix, &
                                   Matrix2D       = PropAirTemperature%Field,            &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,      &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyAirTemperature - ModuleAtmosphere - ERR01'

        endif

    end subroutine ModifyAirTemperature

    !------------------------------------------------------------------------------

    subroutine ModifyRelativeHumidity (PropRelativeHumidity)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropRelativeHumidity

        !Begin-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        if (PropRelativeHumidity%ID%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = PropRelativeHumidity%ID%ObjFillMatrix,&
                                   Matrix2D       = PropRelativeHumidity%Field,           &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,       &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyAirTemperature - ModuleAtmosphere - ERR01'

        endif

    end subroutine ModifyRelativeHumidity

    !------------------------------------------------------------------------------

    subroutine ModifyWindModulus (PropWindModulus)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropWindModulus
        type(T_Property), pointer                   :: PropWindVelX
        type(T_Property), pointer                   :: PropWindVelY


        !Begin-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        if (PropWindModulus%ID%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = PropWindModulus%ID%ObjFillMatrix,    &
                                   Matrix2D       = PropWindModulus%Field,               &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,      &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyWindModulus - ModuleAtmosphere - ERR01'

        elseif (.not. PropWindModulus%Constant) then

            !Searches Wind X
            call SearchProperty(PropWindVelX, WindVelocityX_, .false., STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
            
                !Seaches Wind Y
                call SearchProperty(PropWindVelY, WindVelocityY_, .false., STAT = STAT_CALL)
                if (STAT_CALL == SUCCESS_) then

                    PropWindModulus%Field = sqrt(PropWindVelX%Field ** 2.0 + PropWindVelY%Field ** 2.0)

                endif
            endif

        endif

    end subroutine ModifyWindModulus

    !----------------------------------------------------------------------------

    subroutine ModifyWindDirection (PropWindDirection)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropWindDirection

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------


        if (PropWindDirection%ID%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = PropWindDirection%ID%ObjFillMatrix,     &
                                   Matrix2D       = PropWindDirection%Field,                &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,         &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyWindDirection - ModuleAtmosphere - ERR01'

        endif

    end subroutine ModifyWindDirection

    !------------------------------------------------------------------------------

    subroutine ModifyCloudCover (PropCloudCover)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropCloudCover

        !Begin-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: i,j
        type (T_Property), pointer                  :: PropTransmitivity
        integer                                        :: CHUNK


        !Points to ATMTransmitivity
        call SearchProperty(PropTransmitivity, AtmTransmitivity_, .true., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyCloudCover - ModuleAtmosphere - ERR01'


        if (PropCloudCover%ID%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = PropCloudCover%ID%ObjFillMatrix,     &
                                   Matrix2D       = PropCloudCover%Field,                &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,      &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyCloudCover - ModuleAtmosphere - ERR02'

        else

            if (.not. PropCloudCover%Constant) then
                
                !Calculate Cloud Cover Field
                call ComputeCloudCover (PropCloudCover)
                
            endif

        endif

        if (MonitorPerformance) then
            call StartWatch ("ModuleAtmosphere", "ModifyCloudCover")
        endif

        !Update Transmitivity
        select case (Me%CloudCoverMethod)

            case (CloudFromRadiation)

                CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
                !$OMP PARALLEL PRIVATE(i,j)
                !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%ExternalVar%MappingPoints2D(i, j) == 1) then
                        PropTransmitivity%Field(i, j) = PropCloudCover%Field (i, j )
                    endif
                enddo
                enddo
                !$OMP END DO
                !$OMP END PARALLEL
                    
            case default
               
                CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
                !$OMP PARALLEL PRIVATE (i,j)
                !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%ExternalVar%MappingPoints2D(i, j) == 1) then
                        PropTransmitivity%Field(i, j) = 1.0 - 0.65 * PropCloudCover%Field (i, j ) ** 2.
                    endif
                enddo
                enddo
                !$OMP END DO
                !$OMP END PARALLEL
                
        end select        

        if (MonitorPerformance) then
            call StopWatch ("ModuleAtmosphere", "ModifyCloudCover")
        endif        

    end subroutine ModifyCloudCover

    !----------------------------------------------------------------------------


    subroutine ComputeCloudCover (PropCloudCover)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropCloudCover

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: Julday, STAT_CALL
        real                                        :: Declination
        real                                        :: Hour, Minute, Second
        real                                        :: LatitudePI, LongitudePI
        real                                        :: SunstW, PossibleSunHours
        real                                        :: GmtReference, RacingWithTheSun
        real                                        :: HourAngle, QSO, SunHighAngle
        type(T_Property), pointer                   :: PropSunHours, PropSolarRadiation       
        integer                                        :: CHUNK
        

        call JulianDay  (Me%ActualTime, Julday)
               
        select case (Me%CloudCoverMethod)

        case (CloudFromSunHours)

            !GMT Reference
            call GetGmtReference(Me%ObjTime, GmtReference, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeCloudCover - ModuleAtmosphere - ERR00'

            !Points to Sun hours
            call SearchProperty(PropSunHours, SunHours_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeCloudCover - ModuleAtmosphere - ERR01'

            Declination = SunDeclination_(JulDay)

            !Sun "Racing"
            RacingWithTheSun = RacingWithTheSun_(JulDay)

            call ExtractDate(Me%ActualTime, Hour = Hour, Minute = Minute, Second = Second )        
            Hour = Hour + Minute/60. + Second/3600.

            if (MonitorPerformance) then
                call StartWatch ("ModuleAtmosphere", "ComputeCloudCover")
            endif

            CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
            !$OMP PARALLEL PRIVATE(i,j,LatitudePI,LongitudePI,SunstW,PossibleSunHours)
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExternalVar%MappingPoints2D(i, j) == 1) then
            
                    LatitudePI  = Me%ExternalVar%Latitude(i, j) * PI / 180.0
                    LongitudePI = Me%ExternalVar%Longitude(i, j) 

                    !Angstrong formula [Rs = (a+b * n/N) * QSO]
                    !Sunset angle
                    SunstW              = Acos ( -tan (LatitudePI) * tan (Declination) )

                    !Possible sun hours
                    PossibleSunHours    = 24/PI * SunstW

                    !total atenuationatenuation
                    PropCloudCover%Field(i,j) = 1.0 - ( PropSunHours%Field (i, j ) / PossibleSunHours )
                
                endif

            enddo
            enddo   
            !$OMP END DO
            !$OMP END PARALLEL
        
            if (MonitorPerformance) then
                call StopWatch ("ModuleAtmosphere", "ComputeCloudCover")
            endif
        
        case (CloudFromRadiation)
        
            !GMT Reference
            call GetGmtReference(Me%ObjTime, GmtReference, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeCloudCover - ModuleAtmosphere - ERR02'

            !Points to Solar Radiation File
            call SearchProperty(PropSolarRadiation, SolarRadiation_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeCloudCover - ModuleAtmosphere - ERR03'

            call ExtractDate(Me%ActualTime, Hour = Hour, Minute = Minute, Second = Second )
            call JulianDay  (Me%ActualTime, Julday)

            Hour = Hour + Minute/60. + Second/3600.

            !Sun "Racing"
            RacingWithTheSun = RacingWithTheSun_(JulDay)

            !Declination of the Sun
            Declination = SunDeclination_(JulDay)
       
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
                if (Me%ExternalVar%MappingPoints2D(i, j) == 1) then

                    LatitudePI  = Me%ExternalVar%Latitude (i, j) * PI / 180.0
                    LongitudePI = Me%ExternalVar%Longitude(i, j) 
          
                    !Hour angle 
                    HourAngle = HourAngle_ (Hour, LongitudePI, GmtReference, RacingWithTheSun)
                   
                    !use a sunset angle of 10 as safety factor
                    SunHighAngle =  Asin( sin(LatitudePI) * sin(Declination)            +   &
                                   cos(LatitudePI) * cos(Declination) * cos(HourAngle)) *   &
                                   180./PI

                    if (SunHighAngle <= 10. )then
                        
                        !It's nigth use sunset angle of 10º                        
                        QSO       = TOARadiation (LatitudePI  = LatitudePI ,    &
                                                  Declination = Declination,    &
                                                  Julday      = Julday )     
                                                               
                        PropCloudCover%Field (i,j) = Me%LastRadiation(i,j) / QSO                    
                    else
                                                                                             
                        !Clar day radiation
                        QSO = TOARadiation (LatitudePI, Declination, HourAngle, Julday)
                       
                        if (QSO==0.0) then
                            PropCloudCover%Field (i,j) = 0.595
                        else
                            PropCloudCover%Field (i,j) = PropSolarRadiation%Field (i,j)/QSO
                        endif
                    endif

                    !correct the field
                    if (PropCloudCover%Field (i,j) < 0) then

                        PropCloudCover%Field (i,j) = 0.595

                    elseif (PropCloudCover%Field (i,j) > 1) then

                        PropCloudCover%Field (i,j) = 1.0
                    
                    endif
                endif

            end do
            end do



        case (CloudFromRandom)

            call JulianDay  (Me%ActualTime, Julday)
            if (Julday /= Me%LastCalculateRandomCloud) then
                PropCloudCover%Field = RandomCloud(Julday)
                Me%LastCalculateRandomCloud = Julday
            endif

        end select

    end subroutine ComputeCloudCover

    !--------------------------------------------------------------------------

    subroutine ModifyIrrigation(PropIrrigation)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropIrrigation

        !Local-----------------------------------------------------------------
        real                                        :: ConversionFactor
        integer                                     :: ILB, IUB, JLB, JUB, i, j, STAT_CALL

        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        if (PropIrrigation%ID%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = PropIrrigation%ID%ObjFillMatrix,     &
                                   Matrix2D       = PropIrrigation%Field,                &
                                   PointsToFill2D = Me%ExternalVar%MappingPoints2D,      &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyIrrigation - ModuleAtmosphere - ERR01'

        endif



        if (trim(adjustl(PropIrrigation%ID%Units)) /= 'm3/s') then

            select case (PropIrrigation%ID%Units)

            case ('mm/day')

                ConversionFactor = 1. / 86400.        !In mm/s

            case ('mm/month')

                ConversionFactor = 12. / 86400 / 365. !In mm/s

            case ('mm')

                if(.not.PropIrrigation%NoInterpolateValueInTime) then
                    write(*,*)'Invalid Irrigation Units'
                    write(*,*)'Use m3/s, mm/day or mm/month'
                    stop 'ModifyIrrigation - ModuleAtmosphere - ERR01b'
                endif

                ConversionFactor = 1.                 !In mm

            case default

                write(*,*)'Invalid Irrigation Units'
                write(*,*)'Use m3/s, mm/day or mm/month'
                stop 'ModifyIrrigation - ModuleAtmosphere - ERR02'

            end select

            do j = JLB, JUB
            do i = ILB, IUB
                            
                PropIrrigation%Field(i, j) = PropIrrigation%Field(i, j)         * &
                                            (Me%ExternalVar%GridCellArea(i, j)  * & 
                                             ConversionFactor) / 1000. !In m3/s

            enddo
            enddo
        
        endif


    end subroutine ModifyIrrigation
    
    !------------------------------------------------------------------------------

    subroutine ModifyPropByIrri 

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropertyX

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        PropertyX => Me%FirstAtmosphereProp

do1 :   do while (associated(PropertyX)) 

            if (PropertyX%PropAddedByIrri) then

                if (PropertyX%ID%SolutionFromFile) then

                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,     &
                                           Matrix2D       = PropertyX%Field,                &
                                           PointsToFill2D = Me%ExternalVar%MappingPoints2D, &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModifyPropByIrri - ModuleAtmosphere - ERR01'

                endif
                
            endif

            PropertyX => PropertyX%Next

        end do do1

    end subroutine ModifyPropByIrri

    !--------------------------------------------------------------------------

    subroutine ClimatologicSolarRadiation(PropSolarRadiation)

        !Arguments---------------------------------------------------------------
        type(T_Property), pointer     :: PropSolarRadiation

        !External----------------------------------------------------------------
        type (T_Property), pointer    :: PropCloudCover
        type (T_Property), pointer    :: PropTransmitivity
        real                          :: Hour
        real                          :: Minute
        real                          :: Second
        integer                       :: JulDay
        
        !Local-------------------------------------------------------------------
        real                          :: LatitudePI, Declination, LongitudePI 
        real                          :: HourAngle
        real                          :: SunHigh, QSO
        real                          :: AtmTransf
        real                          :: RacingWithTheSun
        real                          :: GmtReference                  !Time zone: GMT +X
        integer                       :: i, j, STAT_CALL

        !Begin-------------------------------------------------------------------

        !GMT Reference
        call GetGmtReference(Me%ObjTime, GmtReference, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ClimatologicSolarRadiation - ModuleAtmosphere - ERR02' 

        !Points to ATMTransmitivity
        call SearchProperty(PropTransmitivity, AtmTransmitivity_, .true., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ClimatologicSolarRadiation - ModuleAtmosphere - ERR03'

        !Gets CloudCover
        call SearchProperty(PropCloudCover, CloudCover_, .true., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ClimatologicSolarRadiation - ModuleAtmosphere - ERR04'
 

        call ExtractDate(Me%ActualTime, Hour = Hour, Minute = Minute, Second = Second )
        call JulianDay  (Me%ActualTime, Julday)

        Hour = Hour + Minute/60. + Second/3600.

        !Sun "Racing"
        RacingWithTheSun = RacingWithTheSun_(JulDay)

        !Declination of the Sun
        Declination = SunDeclination_(JulDay)
       
do1 :   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2 :   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%ExternalVar%MappingPoints2D(i, j) == 1) then

                LatitudePI  = Me%ExternalVar%Latitude(i, j) * PI / 180.0
                LongitudePI = Me%ExternalVar%Longitude(i, j) 
          
                !Hour angle 
                HourAngle = HourAngle_ (Hour, LongitudePI, GmtReference, RacingWithTheSun)

                !Sun high
                SunHigh = sin(LatitudePI) * sin(Declination) + cos(LatitudePI) * cos(Declination) * cos(HourAngle)
    
                !Top atmosphere sun radiation
                QSO = TOARadiation (LatitudePI, Declination, HourAngle, Julday)

                !Fraction absorbed by Atmosferic 
                AtmTransf = AtmosphericTransmitivity_(SunHigh)

                !Final Radiation which hits surface (Albedo must be considered later)
                PropSolarRadiation%Field(i, j) = QSO * AtmTransf * PropTransmitivity%Field(i, j) 

            endif

        end do do2
        end do do1


    end subroutine ClimatologicSolarRadiation

    
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
    
    !This function computes sun height as a funtion of the local (latitude), &
    !the day (declination) and the (solar) hour.
    !Reference?
    real function SunHigh_(LatitudePI, Declination, HourAngle)

        !Arguments-------------------------------------------------------------
        real                                        :: LatitudePI, Declination
        real                                        :: HourAngle

        !Local-----------------------------------------------------------------

        !Sun high
        SunHigh_ = sin(LatitudePI) * sin(Declination) + cos(LatitudePI) * cos(Declination) * cos(HourAngle)   

    end function SunHigh_

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

    real function AtmosphericTransmitivity_(SunHigh)

        !Arguments-------------------------------------------------------------
        real                                        :: SunHigh

        !Local-----------------------------------------------------------------
        real                                        :: AtmDir, AtmDif

        !Fraction absorbed by Atmosferic 
        !STM method ? Reference for this
        if (SunHigh .lt. 1.0E-3) then !Flavio
          AtmDir    = 0.0
        else
          AtmDir    = 0.74 ** (1.0 / SunHigh)
        endif 
        AtmDif      = (0.91 - AtmDir) / 2.0
        AtmosphericTransmitivity_   = AtmDir + AtmDif

    end function AtmosphericTransmitivity_

    !--------------------------------------------------------------------------

    subroutine CEQUALW2SolarRadiation(PropSolarRadiation)

        !Arguments---------------------------------------------------------------
        type(T_Property), pointer     :: PropSolarRadiation

        !External----------------------------------------------------------------
        type (T_Property), pointer    :: PropCloudCover
        real                          :: RandomCloud_ = FillValueReal
        integer                       :: JulDay

        !Local-------------------------------------------------------------------
        real                          :: Hour, Minute, Second, JDay
        real                          :: Nebulosity
        real                          :: sro, A0, Taud, Declination, Standart
        real                          :: EQTNew, HH
        integer                       :: Iday  
        integer                       :: ILB, IUB, JLB, JUB, i, j, STAT_CALL
        logical                       :: UseRandomCloud

        !Begin-------------------------------------------------------------------

        call ExtractDate(Me%ActualTime, Hour = Hour, Minute = Minute, Second = Second)
        
        call JulianDay  (Me%ActualTime, Julday)

        call GetGmtReference(Me%ObjTime, GmtReference = Standart)

        !Gets CloudCover
        call SearchProperty(PropertyX = PropCloudCover,                         &
                            PropertyXIDNumber = CloudCover_ , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= NOT_FOUND_ERR_)            &
            stop 'CEQUALW2SolarRadiation - ModuleAtmosphere - ERR01' 

        if (STAT_CALL == NOT_FOUND_ERR_) then
            RandomCloud_   = RandomCloud(Julday)
            UseRandomCloud = .true.
        else
            UseRandomCloud = .false.
        endif


        ILB  = Me%WorkSize%ILB
        IUB  = Me%WorkSize%IUB
        JLB  = Me%WorkSize%JLB
        JUB  = Me%WorkSize%JUB

        JDay = Julday + Hour/24. + Minute/24./60. + Second/24./3600.

        Standart      =  15.0*Standart !15 degres per hour (Earths Velocity)

        Hour          =  (JDay-int(JDay))*24.

        Iday          =  JDay-((int(JDay/365))*365) 

        Iday          =  Iday+int(int(JDay/365)/4)

        Taud          =  (2*PI*(Iday-1))/365

        EQTNew        =  0.170*sin(4*PI*(Iday-80)/373)-0.129*sin(2*PI*(Iday-8)/355)

        Declination   =  0.006918 - 0.399912 * cos(Taud) + 0.070257             * &
                      sin(Taud) - 0.006758 * cos(2.0 * Taud) + 0.000907      * &
                      sin(2.0 * Taud) - 0.002697* cos(3.0 * Taud) + 0.001480 * sin(3.0*Taud)

        do j=JLB, JUB
        do i=ILB, IUB
             
             !Original equation
             HH            =  0.261799*(Hour+(Me%ExternalVar%Longitude(i, j) - Standart)*0.0666667+EQTNew-12.0) 
             
             !CEQUAL team change the sign of the longitude contribution for HH 
             !HH            =   0.261799*(Hour-(Me%ExternalVar%Longitude(i, j) - Standart)*0.0666667+EQTNew-12.0)
             !The latter formulation is inconsistent, producing wrong solar radiation computation
             !With the latter, sun rises on west and sets on east. The original formulation, in other end, gives correct
             !results. As so, original formulation was recovered
             
             !Sun High  - 0.0174533 = PI/180
             A0            =  57.2957795 * asin(sin(Me%ExternalVar%Latitude(i,j) * 0.0174533) * sin(Declination) + &
                              cos(Me%ExternalVar%Latitude(i,j) * 0.0174533) * cos(Declination) * cos(HH))
 
             if (A0 > 0.0) then !If <0 =>sun below horizon
                

!                Something made by Pedro Pina that we was not able to justified or explain
!
!                if (UseRandomCloud) then
!                    if (RandomCloud_ .LT. 0.9) then
!                        Nebulosity = (0.650 * RandomCloud_ ** 2.0) * 100
!                    else
!                        Nebulosity = (0.585 * RandomCloud_) * 100 
!                    end if
!                else
!                    Nebulosity = 0.585 * PropCloudCover%Field(i, j)
!                endif 

!               Nebulosity in the CEQUAL W2 is from 0 to 10 
!               Cloud cover is from 0 to 1.

                if (UseRandomCloud) then
                    Nebulosity = RandomCloud_ * 10.
                else
                    Nebulosity = PropCloudCover%Field(i, j) * 10.
                endif
                
                sro  =  2.044*A0+0.1296*A0**2-1.941e-3*A0**3+7.591e-6*A0**4  !Solar contribution
                sro  = (1.0-0.0065*Nebulosity**2)*sro*24.0                   !Cloud cover
                sro  =  sro*0.1314                                           !convert btu/(ft2 day) to w/m2

                PropSolarRadiation%Field(i, j) =  sro
                !PropSolarRadiation%Field(i, j) =  sro*0.94 !Albedo is quantified outside this subroutine
                                                            !In this old version Albedo was considered twice 
              
              else
              
                sro  = 0.0
                PropSolarRadiation%Field(i, j) = 0.0
              
              endif

                             
        end do
        end do


    end subroutine CEQUALW2SolarRadiation


    !--------------------------------------------------------------------------

    function TOARadiation(LatitudePI, Declination, HourAngle, Julday)
        real TOARadiation

        !Arguments-------------------------------------------------------------
        integer , intent(IN)                :: Julday   !Julday = 1 -> 1st of January
        real    , intent(IN)                :: LatitudePI
        real    , intent(IN)                :: Declination
        real    , intent(IN), optional      :: HourAngle
        
        !Local-----------------------------------------------------------------
        real    :: SunHigh
        real    :: ROrbit 

        !--------------------------------------------------------------------------

        !Sun high  
        if (present (HourAngle))  then    
            SunHigh = sin(LatitudePI) * sin(Declination) + cos(LatitudePI) * cos(Declination) * cos(HourAngle)        
        else
            SunHigh = sin (10 * PI / 180.)
        endif
        !Atmosphere radiation
        if (SunHigh .LT. 0.0) then     !night
            
            TOARadiation = 0.0

        else                           !day
            
            ROrbit    = 1.0 + 0.017 * cos((186.0 - JulDay) * 2.0 * PI / 365.0)
            !Top atmosphere sun radiation
            TOARadiation       = KSun * SunHigh / ROrbit ** 2.0
        end if

        !----------------------------------------------------------------------

    end function TOARadiation

    !--------------------------------------------------------------------------

    real function RandomCloud(Julday)

        !Arguments-------------------------------------------------------------
        integer, intent(IN) :: Julday   !Julday = 1 -> 1st of January

        !Local-----------------------------------------------------------------
        real                :: ran
        real                :: x

        !--------------------------------------------------------------------------

        call RANDOM_NUMBER(ran)

        ran = ran - 0.5

        x = Julday * 2 * PI / 365.0

        RandomCloud = ran + ((cos(x) + 1.0) / 2.0)
        RandomCloud = min(RandomCloud, 1.0)
        RandomCloud = max(RandomCloud, 0.0)

    end function RandomCloud

    !--------------------------------------------------------------------------

    subroutine OutPutResultsHDF5

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: STAT_CALL
        type(T_Time)                                :: Actual, LastTime, EndTime
        integer                                     :: OutPutNumber
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
             
        !Begin----------------------------------------------------------------

        Actual   = Me%ActualTime
        EndTime  = Me%EndTime
        LastTime = Me%LastOutPutHDF5
         
        OutPutNumber = Me%OutPut%NextHDF5

TNum:   if (OutPutNumber <= Me%OutPut%Number)            then 
TOut:   if (Actual .GE. Me%OutPut%OutTime(OutPutNumber)) then 
                
            call ExtractDate   (Actual, AuxTime(1), AuxTime(2), AuxTime(3),          &
                                AuxTime(4), AuxTime(5), AuxTime(6))

First:      if (LastTime.LT.Actual) then 
                    
                !Writes Time
                TimePtr => AuxTime
                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutResultsHDF5 - ModuleAtmosphere - ERR00'

                call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                                     Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutResultsHDF5 - ModuleAtmosphere - ERR01'
           
                Me%LastOutPutHDF5 = Actual
       
            endif First

            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5,                                &
                                  Me%WorkSize%ILB,                           &
                                  Me%WorkSize%IUB,                           &
                                  Me%WorkSize%JLB,                           &
                                  Me%WorkSize%JUB,                           &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutResultsHDF5 - ModuleAtmosphere - ERR02'

            PropertyX => Me%FirstAtmosphereProp
            do while (associated(PropertyX))

                if (PropertyX%OutputHDF) then

                    call HDF5WriteData   (Me%ObjHDF5,                            &
                                          "/Results/"//trim(PropertyX%ID%Name),          &
                                          trim(PropertyX%ID%Name),                       &
                                          trim(PropertyX%ID%Units),                      &
                                          Array2D = PropertyX%Field,                     &
                                          OutputNumber = OutPutNumber,                   &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutResultsHDF5 - ModuleAtmosphere - ERR03'

                endif

                PropertyX => PropertyX%Next

            enddo

            Me%OutPut%NextHDF5 = OutPutNumber + 1

            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutResultsHDF5 - ModuleAtmosphere - ERR06'
            
        endif  TOut
        endif  TNum

    end subroutine OutPutResultsHDF5

    !--------------------------------------------------------------------------

    subroutine OutPut_Statistics (Value2D, StatisticsID)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer :: Value2D
        integer                       :: StatisticsID

        !Local-----------------------------------------------------------------
        integer                       :: MethodStatistic, Value2DStat2D
        integer                       :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetStatisticMethod (StatisticsID, MethodStatistic, STAT = STAT_CALL)                                     
                                                                                    
        if (STAT_CALL /= SUCCESS_)                                                      & 
            call SetError (FATAL_, INTERNAL_, 'OutPut_Statistics - Atmosphere - ERR01')
                                                                                    
        call GetStatisticParameters (StatisticsID,                                   &
                                     Value2DStat2D = Value2DStat2D,                  &
                                     STAT          = STAT_CALL)                        
                                                                                    
        if (STAT_CALL /= SUCCESS_)                                                       &
            call SetError (FATAL_, INTERNAL_, 'OutPut_Statistics - Atmosphere - ERR02')
                                                                                    
                                                                                    
        if (MethodStatistic /= Value2DStat2D)                                            &
            call SetError (FATAL_, INTERNAL_, 'OutPut_Statistics - Atmosphere - ERR03')
                                                                                    
                                                                                    
        call ModifyStatistic (StatisticsID,                                              &
                              Value2D       = Value2D,                                   &
                              WaterPoints2D = Me%ExternalVar%MappingPoints2D,            &
                              STAT          = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_)                                                       &
            call SetError (FATAL_, INTERNAL_, 'OutPut_Statistics - Atmosphere - ERR04')

    end subroutine OutPut_Statistics

    !--------------------------------------------------------------------------

    subroutine OutPut_TimeSeries(PropertyX)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        if (PropertyX%TimeSerie) then

            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = PropertyX%Field, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleAtmosphere - ERR01'
                
        endif

    end subroutine OutPut_TimeSeries

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillAtmosphere(AtmosphereID, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: AtmosphereID
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                         :: ready_
        integer                         :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                         :: STAT_, nUsers

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_


        call Ready(AtmosphereID, ready_) 

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mATMOSPHERE_,  Me%InstanceID)

cd2 :       if (nUsers == 0) then

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillAtmosphere - ModuleAtmosphere - ERR01'

                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillAtmosphere - ModuleAtmosphere - ERR02'

                nUsers = DeassociateInstance (mGRIDDATA_, Me%ObjGridData)
                if (nUsers == 0) stop 'KillAtmosphere - ModuleAtmosphere - ERR03'

                if (Me%OutPut%True) then
                    !Kills HDF 5
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillAtmosphere - ModuleAtmosphere - ERR06'
                endif

                !Kills the TimeSerie
                if ((Me%ObjTimeSerie > 0)) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                                 &
                        stop 'KillAtmosphere - ModuleAtmosphere - ERR07'
                endif

                call DeallocateVariables

                !ObjAtmosphere
                call DeallocateInstance 

                STAT_ = SUCCESS_

            end if cd2
        else cd1
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_


    end subroutine KillAtmosphere

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance 

        !Local-----------------------------------------------------------------
        type (T_Atmosphere), pointer           :: AuxAtmosphere
        type (T_Atmosphere), pointer           :: PreviousAtmosphere

        !Updates pointers
        if (Me%InstanceID == FirstObjAtmosphere%InstanceID) then
            FirstObjAtmosphere => FirstObjAtmosphere%Next
        else
            PreviousAtmosphere => FirstObjAtmosphere
            AuxAtmosphere      => FirstObjAtmosphere%Next
            do while (AuxAtmosphere%InstanceID /= Me%InstanceID)
                PreviousAtmosphere => AuxAtmosphere
                AuxAtmosphere      => AuxAtmosphere%Next
            enddo

            !Now update linked list
            PreviousAtmosphere%Next => AuxAtmosphere%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance



    !--------------------------------------------------------------------------

    Subroutine DeAllocateVariables

        !Local------------------------------------------------------------------
        type (T_Property), pointer :: PropertyX
        integer                    :: STAT_CALL

        
        !----------------------------------------------------------------------


        ! Deallocates all the water properties 

        PropertyX => Me%FirstAtmosphereProp

do1 :   do while(associated(PropertyX))  

            if (associated(PropertyX%Field, PropertyX%FieldGrid)) then

                deallocate(PropertyX%Field    )

                nullify   (PropertyX%Field    )
                nullify   (PropertyX%FieldGrid)           
            else
                deallocate(PropertyX%Field    )
                deallocate(PropertyX%FieldGrid)
                nullify   (PropertyX%Field    )
                nullify   (PropertyX%FieldGrid)
            endif

            !Kill Statistics
            if (PropertyX%Statistics%ON) then

                call KillStatistic (PropertyX%Statistics%ID, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    call SetError(FATAL_, INTERNAL_, 'DeAllocateVariables - ModuleAtmosphere - ERR01') 

            endif

            if (PropertyX%ID%ObjFillMatrix /= 0) then
                call KillFillMatrix (PropertyX%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'DeAllocateVariables - ModuleAtmosphere - ERR02'
            endif

           
            PropertyX => PropertyX%Next

        end do do1

        !Sets the number of properties equal to the FillValueInt
        Me%PropertiesNumber = FillValueInt

        Nullify   (Me%FirstAtmosphereProp,Me%LastAtmosphereProp)


        !----------------------------------------------------------------------

    End Subroutine DeAllocateVariables   

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Ready (AtmosphereID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                            :: AtmosphereID
        integer                            :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (AtmosphereID > 0) then
            call LocateObjAtmosphere(AtmosphereID)
            ready_ = VerifyReadLock (mATMOSPHERE_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjAtmosphere (AtmosphereID)

        !Arguments-------------------------------------------------------------
        integer                            :: AtmosphereID

        !Local-----------------------------------------------------------------

        Me => FirstObjAtmosphere
        do while (associated (Me))
            if (Me%InstanceID == AtmosphereID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                                        &
            stop 'ModuleAtmosphere - LocateObjAtmosphere - ERR01'

    end subroutine LocateObjAtmosphere

    !--------------------------------------------------------------------------


end module ModuleAtmosphere

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
