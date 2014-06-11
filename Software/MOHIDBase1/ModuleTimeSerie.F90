!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Time Serie
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to Write/Read Time Series
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

Module ModuleTimeSerie

    use ModuleGlobalData
    use ModuleTime                 
    use ModuleEnterData 
    use ModuleDrawing           

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartTimeSerie
    private ::      AllocateInstance
    private ::      ReadNumberOfTimeSeries
    private ::      ReadTimeSeriesLocation
    private ::      ReadTimeSeriesTimes
    private ::      VerifyTimeSerieLocation 
    private ::      ReadOnlyOneTimeSerieData
    private ::      AllocateTimeSerieBuffer
    private ::      OpenTimeSerieFiles

    public  :: StartTimeSerieInput


    !Modifier
    public  :: WriteTimeSerie
    public  :: WriteTimeSerieLineNow
    public  :: WriteTimeSerieLine
    private ::      WriteBufferToFile
    public  :: WriteSpecificTimeSerieLine
    public  :: CorrectsCellsTimeSerie
    public  :: TryIgnoreTimeSerie
    public  :: ReseatCurrentIndex

    !Selector
    public  :: GetNumberOfTimeSeries
    public  :: GetTimeSerieLocation
    public  :: GetTimeSerieValue
    public  :: GetTimeSerieTimeLimits
    public  :: GetTimeSerieIntegral
    public  :: GetTimeSerieDT
    public  :: GetTimeSerieHeader
    public  :: GetTimeSerieName
    public  :: GetTimeSerieNextOutput
    public  :: GetTimeSerieTimeUnits
    public  :: GetTimeSerieInitialData
    public  :: GetTimeSerieDTForNextEvent
    public  :: GetTimeSerieDataMatrix
    public  :: GetTimeSerieDataColumns
    public  :: GetTimeSerieDataValues
    public  :: GetTimeSerieTimeFrameIndexes
    public  :: GetTimeSerieCycle
    public  :: GetTimeSerieCheckDate
    public  :: GetTimeSerieTimeOfDataset
    public  :: GetTimeSerieTimeOfNextDataset
    public  :: GetTimeSerieValueForIndex

    !Destructor
    public  ::  KillTimeSerie
    private ::      CloseTimeSerieFiles
    private ::      DeallocateInstance

    !Management
    private ::      Ready

    !Interface-----------------------------------------------------------------

    private :: WriteTimeSerie1
    private :: WriteTimeSerie3
    private :: WriteTimeSerie4
    interface  WriteTimeSerie
        module procedure WriteTimeSerie1
        module procedure WriteTimeSerie3
        module procedure WriteTimeSerie4
    end interface WriteTimeSerie

    !Parameter-----------------------------------------------------------------
#ifdef _BIG_MAX_COLUMNS
    integer, parameter                       :: MaxColumns      = 500 !Maximum number of input columns
#else
    integer, parameter                       :: MaxColumns      = 200 !Maximum number of input columns
#endif


    character(len = StringLength), parameter :: block_begin     = '<BeginTimeSerie>'
    character(len = StringLength), parameter :: block_end       = '<EndTimeSerie>'

    character(len = StringLength), parameter :: begin_residual  = '<BeginResidual>'
    character(len = StringLength), parameter :: end_residual    = '<EndResidual>'


    !Types---------------------------------------------------------------------

    type       T_TimeSerie
        real                                        :: DT                 = null_real
        real                                        :: FirstDT            = 0
        integer                                     :: TotalOutPutsNumber = null_int
        integer                                     :: BufferSize         = null_int      !Number of instantes which can kept
        integer                                     :: BufferCount        = null_int      !Current number of instantes
        integer                                     :: UnitNumber         = null_int
        integer                                     :: LocalizationI      = null_int
        integer                                     :: LocalizationJ      = null_int
        integer                                     :: LocalizationK      = null_int
        real                                        :: Latitude           = null_real
        real                                        :: Longitude          = null_real
        real                                        :: CoordX             = null_real
        real                                        :: CoordY             = null_real
        logical                                     :: CoordON            = .false.
        real                                        :: DepthLevel         = null_real
        logical                                     :: DepthON            = .false.
        type (T_Time)                               :: BeginOutPut                  !Limit
        type (T_Time)                               :: EndOutPut                    !Limit
        type (T_Time)                               :: NextOutput
        real, dimension(:, :), pointer              :: TimeSerieData      => null()
        type (T_Time), dimension(:), pointer        :: TimeBuffer
        real, dimension(:), pointer                 :: ResidualValues     => null()
        real                                        :: ResidualTime       = null_real
        type (T_Time)                               :: LastResidual
        character(len=PathLength)                   :: FromBlockFileName  = null_str
        character(len=PathLength)                   :: FileName           = null_str
        logical                                     :: IgnoreON           = .false.
    end type T_TimeSerie

    type      T_TimeSerieInOutPut

        !Instance ID
        integer                                     :: InstanceID         = null_int
        character(PathLength)                       :: ModelName          = null_str
        logical                                     :: ModelNameON        = .false.
        logical                                     :: ReplacePathON      = .false.
        character(PathLength)                       :: ReplacePath        = null_str

        !Time Serie Input
        real,    dimension(:,:), pointer            :: DataMatrix         => null()
        integer, dimension(:  ), pointer            :: ColumnsRead        => null()
        integer, dimension(:  ), pointer            :: FileColumns        => null()
        integer                                     :: DataValues         = null_int
        integer                                     :: DataColumns        = null_int

        !Time Serie Output
        integer                                     :: NumberOfProperties   = null_int
        integer                                     :: InternalPropertyCount= null_int
        integer                                     :: NumberOfTimeSeries   = 0
        integer                                     :: MaxBufferSize        = null_int
        logical                                     :: Points               = IDLE
        logical                                     :: TimeSerie3D          = .false.
        logical                                     :: TimeSerie2D          = .false.
        logical                                     :: TimeSerie1D          = .false.
        logical                                     :: ComputeResidual      = .true.
        !logical                                     :: IgnoreON             = .false.

        !TimeSerieInput 
        logical                                     :: TimeCycle            = .false.
        character(len=StringLength)                 :: CharTimeUnits        = null_str
        character(len=line_length)                  :: Header               = null_str
        type (T_Time)                               :: InitialData
        type (T_Time)                               :: PreviousTime
        type (T_Time)                               :: NextTime
        type (T_Time)                               :: TimeOfNextDataset
        real                                        :: DTForNextEvent       = null_real
        real                                        :: DTForNextDataset     = null_real
        integer                                     :: NextInstant          = null_int
        integer                                     :: PreviousInstant      = null_int
        integer                                     :: InstantOfNextDataset = null_int
        integer                                     :: CurrentIndex         = 2    
        integer                                     :: StartIndex           = 1 
        integer                                     :: EndIndex             = 1

        !TimeSerieOutput
        type(T_TimeSerie), dimension(:), pointer    :: TimeSerie            => null()
        
        type (T_Polygon), pointer                   :: ModelDomain          => null()
        logical                                     :: ModelDomainON        = .false.

        !Instance of Module_EnterData
        integer                                     :: ObjEnterData         = 0

        !Instance of ModuleTime
        integer                                     :: ObjTime              = 0

        type (T_TimeSerieInOutPut), pointer         :: Next                 => null()

        logical                                     :: UseTabulatedData     = .true.

    end type T_TimeSerieInOutPut

    !Global Variables
    type (T_TimeSerieInOutPut), pointer             :: FirstTimeSerie       => null()
    type (T_TimeSerieInOutPut), pointer             :: Me                   => null()

    !--------------------------------------------------------------------------

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartTimeSerie(TimeSerieID, ObjTime,                                      &
                              TimeSerieDataFile, PropertyList, Extension, WaterPoints3D, &
                              WaterPoints2D, WaterPoints1D, ResultFileName, Instance,    &
                              ModelName, CoordX, CoordY, UseTabulatedData,               &
                              HavePath, Comment, ModelDomain, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: TimeSerieID
        integer                                         :: ObjTime
        character(len=*), intent(IN )                   :: TimeSerieDataFile
        character(len=*), dimension(:), pointer         :: PropertyList
        character(len=*), intent(IN )                   :: Extension
        integer, dimension(:,:,:), optional, pointer    :: WaterPoints3D
        integer, dimension(:,:  ), optional, pointer    :: WaterPoints2D
        integer, dimension(:    ), optional, pointer    :: WaterPoints1D
        character(len=*), optional, intent(IN )         :: ResultFileName
        character(len=*), optional, intent(IN )         :: Instance  
        character(len=*), optional, intent(IN )         :: ModelName
        real, optional                                  :: CoordX
        real, optional                                  :: CoordY
        logical, optional, intent(IN )                  :: UseTabulatedData
        logical, optional, intent(IN )                  :: HavePath
        character(len=*), optional, intent(IN )         :: Comment
        type (T_Polygon), pointer, optional, intent(IN ):: ModelDomain
        integer, optional, intent(OUT)                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: ready_ , STAT_
        integer                                         :: FromFile
        integer                                         :: flag, ret
        integer                                         :: iTimeSerie, j
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mTimeSerie_)) then
            nullify (FirstTimeSerie)
            call RegisterModule (mTimeSerie_) 
        endif

        call Ready(TimeSerieID, ready_)    

if0 :   if (ready_ .EQ. OFF_ERR_) then

            !Allocates Instance
            call AllocateInstance

            nullify (Me%DataMatrix)
            nullify (Me%ColumnsRead)
            nullify (Me%FileColumns)
            nullify (Me%TimeSerie)

            !Associates module Time
            Me%ObjTime = AssociateInstance   (mTIME_, ObjTime)

            if (present(ModelName)) then
                Me%ModelName    = ModelName
                Me%ModelNameON  = .true.
            endif

            if (present(UseTabulatedData)) then
                Me%UseTabulatedData = UseTabulatedData
            endif
            
            if (present(ModelDomain)) then
                Me%ModelDomainON = .true.
                Me%ModelDomain   => ModelDomain
            else
                Me%ModelDomainON = .false.
            endif

            !Constructs EnterData
            call ConstructEnterData(Me%ObjEnterData, TimeSerieDataFile, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartTimeSerie - ModuleTimeSerie - ERR10'

            !Get Extract type
            call GetExtractType    (FromFile  = FromFile )

            !The maximum BufferSize is set here to 0.1Mb (for each property)
            !This lets perfrom 25000 outputs to the buffer (considering each output of 4 bytes)
            call GetData(Me%MaxBufferSize,                                      &
                         Me%ObjEnterData,                                       &
                         flag,                                                  &
                         SearchType   = FromFile,                               &
                         keyword      ='MAX_BUFFER_SIZE',                       &
                         Default      = 100000,                                 &
                         ClientModule ='ModuleTimeSerie',                       &
                         STAT         = STAT_CALL)        
            if (STAT_CALL .NE. SUCCESS_)                                        &
                call SetError(FATAL_, KEYWORD_, "Subroutine StartTimeSerie - ModuleTimeSerie. ERR20") 

            call GetData(Me%ComputeResidual,                                    &
                         Me%ObjEnterData,                                       &
                         flag,                                                  &
                         SearchType   = FromFile,                               &
                         keyword      ='COMPUTE_RESIDUAL',                      &
                         Default      = .true.,                                 &
                         ClientModule ='ModuleTimeSerie',                       &
                         STAT         = STAT_CALL)        
            if (STAT_CALL .NE. SUCCESS_)                                        &
                call SetError(FATAL_, KEYWORD_, "Subroutine StartTimeSerie - ModuleTimeSerie. ERR30") 

            call GetData(Me%ReplacePath,                                        &
                         Me%ObjEnterData,                                       &
                         flag,                                                  &
                         SearchType   = FromFile,                               &
                         keyword      ='REPLACE_PATH',                          &
                         Default      = '****',                                 &
                         ClientModule ='ModuleTimeSerie',                       &
                         STAT         = STAT_CALL)        
            if (STAT_CALL .NE. SUCCESS_)                                        &
                call SetError(FATAL_, KEYWORD_, "Subroutine StartTimeSerie - ModuleTimeSerie. ERR40") 

            if (flag > 0) Me%ReplacePathON = .true. 

            !call GetData(Me%IgnoreON,                                           &
            !             Me%ObjEnterData,                                       &
            !             flag,                                                  &
            !             SearchType   = FromFile,                               &
            !             keyword      ='IGNORE_ON',                             &
            !             Default      = .false.,                                &
            !             ClientModule ='ModuleTimeSerie',                       &
            !             STAT         = STAT_CALL)        
            !if (STAT_CALL .NE. SUCCESS_)                                        &
            !    call SetError(FATAL_, KEYWORD_, "Subroutine StartTimeSerie - ModuleTimeSerie. ERR50") 

            !Stores the number of properties
            Me%NumberOfProperties = size(PropertyList)

            if (present(WaterPoints1D)) Me%TimeSerie1D = .true.
            if (present(WaterPoints2D)) Me%TimeSerie2D = .true.
            if (present(WaterPoints3D)) Me%TimeSerie3D = .true.

            if (present(ResultFileName)) then


                !Reads the data of only one time serie 
                if (present(CoordX) .and. present(CoordY)) then
                    call ReadOnlyOneTimeSerieData(ResultFileName, CoordX, CoordY)
                else
                    call ReadOnlyOneTimeSerieData(ResultFileName)
                endif

            
            else


                !Reads the data
                Me%NumberOfTimeSeries = ReadNumberOfTimeSeries ()
                if (Me%NumberOfTimeSeries == -1 ) then
                    stop 'StartTimeSerie - ModuleTimeSerie - ERR60'
                elseif (Me%NumberOfTimeSeries > 0) then

                    !Allocate TimeSerie
                    allocate (Me%TimeSerie(Me%NumberOfTimeSeries), STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)  stop 'StartTimeSerie - ModuleTimeSerie - ERR70' 

                    do iTimeSerie = 1, Me%NumberOfTimeSeries
                        nullify (Me%TimeSerie(iTimeSerie)%TimeSerieData)
                        nullify (Me%TimeSerie(iTimeSerie)%TimeBuffer)
                        nullify (Me%TimeSerie(iTimeSerie)%ResidualValues)
                    enddo
                
                    ret = ReadTimeSeriesLocation ()
                    if (ret == SUCCESS_) then
                        ret = ReadTimeSeriesTimes ()
                        if (ret /= SUCCESS_) stop 'StartTimeSerie - ModuleTimeSerie - ERR80'
                    else
                        stop 'StartTimeSerie - ModuleTimeSerie - ERR90'
                    endif
                endif

                Me%Points = .true.
                
                !Verifies the location of the time series
                if (Me%Points) then
                    if (Me%TimeSerie3D        ) call VerifyTimeSerieLocation(WaterPoints3D = WaterPoints3D)
                    if (Me%TimeSerie2D        ) call VerifyTimeSerieLocation(WaterPoints2D = WaterPoints2D)
                    if (Me%TimeSerie1D        ) call VerifyTimeSerieLocation(WaterPoints1D = WaterPoints1D)
                    if (Me%TimeSerie1D .and. Me%TimeSerie2D) stop 'StartTimeSerie - ModuleTimeSerie - ERR100'
                    if (Me%TimeSerie1D .and. Me%TimeSerie3D) stop 'StartTimeSerie - ModuleTimeSerie - ERR110'
                    if (Me%TimeSerie2D .and. Me%TimeSerie3D) stop 'StartTimeSerie - ModuleTimeSerie - ERR120'
                endif

            endif

            !Check if there is any time serie with equal names
            do iTimeSerie = 1, Me%NumberOfTimeSeries
                do j=1, Me%NumberOfTimeSeries
                    if (iTimeSerie /= j) then
                        if (trim(Me%TimeSerie(iTimeSerie)%FromBlockFileName) == &
                            trim(Me%TimeSerie(         j)%FromBlockFileName)) then
                            write(*,*) 'Time series can not have equal names'
                            stop 
                        endif 
                    endif
                enddo
            enddo



            !Allocates the Buffer
            call AllocateTimeSerieBuffer

            !Constructs the time serie files
            if (present(Comment)) then
                if (present(ResultFileName)) then
                    if (present(HavePath)) then 
                        call OpenTimeSerieFiles(Extension, &
                                                PropertyList, &
                                                ResultFileName = ResultFileName, &
                                                HavePath = HavePath, &
                                                Comment = Comment)
                    else
                        call OpenTimeSerieFiles(Extension, PropertyList, ResultFileName = ResultFileName, Comment = Comment)
                    endif
                else if (present(Instance)) then
                    call OpenTimeSerieFiles(Extension, PropertyList, Instance = Instance, Comment = Comment)
                else 
                    call OpenTimeSerieFiles(Extension, PropertyList, Comment = Comment)
                endif
            else
                if (present(ResultFileName)) then
                    if (present(HavePath)) then 
                        call OpenTimeSerieFiles(Extension, PropertyList, ResultFileName = ResultFileName, HavePath = HavePath)
                    else
                        call OpenTimeSerieFiles(Extension, PropertyList, ResultFileName = ResultFileName)
                    endif
                else if (present(Instance)) then
                    call OpenTimeSerieFiles(Extension, PropertyList, Instance = Instance)
                else 
                    call OpenTimeSerieFiles(Extension, PropertyList)
                endif
            endif
            


            !Inits internal counts
            Me%InternalPropertyCount = 0


            !Kills EnterData
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartTimeSerie - ModuleTimeSerie - ERR130'

            !Returns ID
            TimeSerieID    = Me%InstanceID

            STAT_ = SUCCESS_

        else
         
            stop 'ModuleTimeSerie - StartTimeSerie - ERR140' 

        end if if0


        if (present(STAT)) STAT = STAT_

                    
    end subroutine StartTimeSerie

    !--------------------------------------------------------------------------

    subroutine AllocateInstance 

            !Local-----------------------------------------------------------------
        type (T_TimeSerieInOutPut), pointer         :: NewTimeSerie
        type (T_TimeSerieInOutPut), pointer         :: PreviousTimeSerie


        !Allocates new instance
        allocate (NewTimeSerie)
        nullify  (NewTimeSerie%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstTimeSerie)) then
            FirstTimeSerie    => NewTimeSerie
            Me                => NewTimeSerie
        else
            PreviousTimeSerie => FirstTimeSerie
            Me                => FirstTimeSerie%Next
            do while (associated(Me))
                PreviousTimeSerie  => Me
                Me                 => Me%Next
            enddo
            Me                     => NewTimeSerie
            PreviousTimeSerie%Next => NewTimeSerie
        endif

        Me%InstanceID = RegisterNewInstance (mTIMESERIE_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    integer function ReadNumberOfTimeSeries ()

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, NumberOfTimeSeries
        logical                             :: BlockFound
        integer                             :: ClientNumber

        !Sees how many Time Series exists
        NumberOfTimeSeries = 0
        BlockFound = .true.
do1:    do while (BlockFound)
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,              &
                                        block_begin, block_end, BlockFound,         &
                                        STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) then
                ReadNumberOfTimeSeries = -1
                return
            endif
            if (BlockFound) then
                NumberOfTimeSeries = NumberOfTimeSeries + 1
            else
                call Block_Unlock(Me%ObjEnterData, ClientNumber)
            endif
        enddo do1

        ReadNumberOfTimeSeries = NumberOfTimeSeries

    end function ReadNumberOfTimeSeries

    !--------------------------------------------------------------------------

    integer function ReadTimeSeriesLocation ()

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL
        logical                             :: BlockFound
        integer                             :: ClientNumber
        integer                             :: iTimeSerie, iflag
        integer                             :: FromBlock, FromFile
        type (T_PointF), pointer            :: TimeSeriesXY


        !Gets parameter from the module EnterData
        call GetExtractType(FromBlock = FromBlock, FromFile = FromFile)


d1:     do iTimeSerie = 1, Me%NumberOfTimeSeries

            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,              &
                                        block_begin, block_end, BlockFound,         &
                                        STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) then
                ReadTimeSeriesLocation = -1
                return
            endif


            !Searches for NAME
            call GetData(Me%TimeSerie(iTimeSerie)%FromBlockFileName,                &
                         Me%ObjEnterData,                                           &
                         iflag,                                                     &
                         SearchType   = FromBlock,                                  &
                         keyword      ='NAME',                                      &
                         Default      = null_str,                                   &
                         ClientModule ='ModuleTimeSerie',                           &
                         STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) then
                ReadTimeSeriesLocation = -1
                return
            endif

                        !Searches for the depth
            call GetData(Me%TimeSerie(iTimeSerie)%DepthLevel,                           &
                         Me%ObjEnterData,                                               &
                         iflag,                                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      ='DEPTH_LEVEL',                                   &
                         ClientModule ='ModuleTimeSerie',                               &
                         STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ ) then
                ReadTimeSeriesLocation = -1
                return
            endif

            if (iflag == 1) then
                Me%TimeSerie(iTimeSerie)%DepthON = .true. 
            endif


                        !Searches for the Latitude
            call GetData(Me%TimeSerie(iTimeSerie)%CoordX,                               &
                         Me%ObjEnterData,                                               &
                         iflag,                                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      ='COORD_X',                                       &
                         ClientModule ='ModuleTimeSerie',                               &
                         STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ ) then
                ReadTimeSeriesLocation = -1
                return
            endif

            if (iflag == 1) then

                        !Searches for the Longitude
                call GetData(Me%TimeSerie(iTimeSerie)%CoordY,                           &
                             Me%ObjEnterData,                                           &
                             iflag,                                                     &
                             SearchType   = FromBlock,                                  &
                             keyword      ='COORD_Y',                                   &
                             ClientModule ='ModuleTimeSerie',                           &
                             STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) then
                    ReadTimeSeriesLocation = -1
                    return
                endif

                if (iflag == 1) then
                    Me%TimeSerie(iTimeSerie)%CoordON = .true.
                else
                    Me%TimeSerie(iTimeSerie)%CoordON = .false.
                    ReadTimeSeriesLocation = -1
                    return
                endif
                
            endif
            
i12:        if (Me%TimeSerie(iTimeSerie)%CoordON .and. Me%ModelDomainON) then
            
                nullify (TimeSeriesXY)
                allocate(TimeSeriesXY)
            
                TimeSeriesXY%X = Me%TimeSerie(iTimeSerie)%CoordX
                TimeSeriesXY%Y = Me%TimeSerie(iTimeSerie)%CoordY
            
                if (.not. IsPointInsidePolygon(TimeSeriesXY, Me%ModelDomain)) then
                    Me%TimeSerie(iTimeSerie)%IgnoreON = .true. 
                endif
                
                deallocate(TimeSeriesXY)
                
            endif i12

i8:         if (.not. Me%TimeSerie(iTimeSerie)%CoordON) then

                if(Me%TimeSerie1D.or.Me%TimeSerie2D.or.Me%TimeSerie3D) then
                    !Searches for the Localization I
                    call GetData(Me%TimeSerie(iTimeSerie)%LocalizationI,                &
                                 Me%ObjEnterData,                                       &
                                 iflag,                                                 &
                                 SearchType   = FromBlock,                              &
                                 keyword      ='LOCALIZATION_I',                        &
                                 ClientModule ='ModuleTimeSerie',                       &
                                 STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_ .or. iflag /= 1) then
                        ReadTimeSeriesLocation = -1
                        return
                    endif
                endif

                if(Me%TimeSerie2D.or.Me%TimeSerie3D) then
                    !Searches for the Localization J
                    call GetData(Me%TimeSerie(iTimeSerie)%LocalizationJ,                &
                                 Me%ObjEnterData,                                       &
                                 iflag,                                                 &
                                 SearchType   = FromBlock,                              &
                                 keyword      ='LOCALIZATION_J',                        &
                                 ClientModule ='ModuleTimeSerie',                       &
                                 STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_ .or. iflag /= 1) then
                        ReadTimeSeriesLocation = -1
                        return
                    endif
                endif

            endif i8

i9:         if (.not. Me%TimeSerie(iTimeSerie)%DepthON) then

                if(Me%TimeSerie3D) then
                    !Searches for the Localization K
                    call GetData(Me%TimeSerie(iTimeSerie)%LocalizationK,                    &
                                 Me%ObjEnterData,                                           &
                                 iflag,                                                     &
                                 SearchType   = FromBlock,                                  &
                                 keyword      ='LOCALIZATION_K',                            &
                                 ClientModule ='ModuleTimeSerie',                           &
                                 STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_ .or. iflag /= 1) then
                        ReadTimeSeriesLocation = -1
                        return
                    endif
                endif

            endif i9

                        !Searches for the Latitude
            call GetData(Me%TimeSerie(iTimeSerie)%Latitude,                             &
                         Me%ObjEnterData,                                               &
                         iflag,                                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      ='LATITUDE',                                      &
                         ClientModule ='ModuleTimeSerie',                               &
                         STAT = STAT_CALL)                                              
            if (STAT_CALL .NE. SUCCESS_ ) then                                          
                ReadTimeSeriesLocation = -1                                             
                return                                                                  
            endif                                                                       
                                                                                        
                        !Searches for the Longitude                                     
            call GetData(Me%TimeSerie(iTimeSerie)%Longitude,                            &
                         Me%ObjEnterData,                                               &
                         iflag,                                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      ='LONGITUDE',                                     &
                         ClientModule ='ModuleTimeSerie',                               &
                         STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ ) then
                ReadTimeSeriesLocation = -1
                return
            endif

            !Unlocks the block
            call Block_Unlock(Me%ObjEnterData, ClientNumber)


        enddo d1

        ReadTimeSeriesLocation = SUCCESS_

    end function ReadTimeSeriesLocation

    !--------------------------------------------------------------------------

    integer function ReadTimeSeriesTimes ()

        !Local-----------------------------------------------------------------
        type (T_Time)                       :: AuxTime, DummyTime
        integer                             :: STAT_CALL
        logical                             :: BlockFound
        integer                             :: ClientNumber
        integer                             :: iTimeSerie, iflag
        integer                             :: FromBlock, FromFile        

        !Gets parameter from the module EnterData
        call GetExtractType(FromBlock = FromBlock, FromFile = FromFile)

        do iTimeSerie = 1, Me%NumberOfTimeSeries
        
            if (Me%TimeSerie(iTimeSerie)%IgnoreON) cycle

            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,              &
                                        block_begin, block_end, BlockFound,         &
                                        STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) then
                ReadTimeSeriesTimes = -1
                return
            endif

            !Searches for the first output DT
            call GetData(Me%TimeSerie(iTimeSerie)%FirstDT,                          &
                         Me%ObjEnterData,                                           &
                         iflag,                                                     &
                         SearchType   = FromFile,                                   &
                         keyword      ='FIRST_OUTPUT_DT',                           &
                         default      = 0.,                                         &
                         ClientModule ='ModuleTimeSerie',                           &
                         STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            & 
                stop 'ReadTimeSeriesTimes - ModuleTimeSerie - ERR10'

            !Searches for the first output time
            call GetData(AuxTime,                                                   &
                         Me%ObjEnterData,                                           &
                         iflag,                                                     &
                         SearchType   = FromFile,                                   &
                         keyword      ='FIRST_OUTPUT_TIME',                         &
                         ClientModule ='ModuleTimeSerie',                           &
                         STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            & 
                stop 'ReadTimeSeriesTimes - ModuleTimeSerie - ERR10'

            !If the first output time is not specified assmues the current one
            if (iflag == 6) then
                Me%TimeSerie(iTimeSerie)%BeginOutPut  = AuxTime
            elseif (iflag == 0) then
                call GetComputeCurrentTime(Me%ObjTime,                              &
                                           Me%TimeSerie(iTimeSerie)%BeginOutPut,    &
                                           STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ReadTimeSeriesTimes - ModuleTimeSerie - ERR10a'
            else
                write(*,*)"Error specifing first instant of the a time serie"
                write(*,*)"Assuming current time"
                write(*,*)"ReadTimeSeriesTimes - ModuleTimeSerie - WRN01."
                call GetComputeCurrentTime(Me%ObjTime,                              &
                                           Me%TimeSerie(iTimeSerie)%BeginOutPut,    &
                                           STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ReadTimeSeriesTimes - ModuleTimeSerie - ERR10b'
            endif

            !Searches for the DT
            call GetData(Me%TimeSerie(iTimeSerie)%DT,                               &
                         Me%ObjEnterData,                                           &
                         iflag,                                                     &
                         SearchType   = FromFile,                                   &
                         keyword      ='DT_OUTPUT_TIME',                            &
                         ClientModule ='ModuleTimeSerie',                           &
                         STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'ReadTimeSeriesTimes - ModuleTimeSerie - ERR11'

            if (iflag == 0) then
                call GetComputeTimeStep(Me%ObjTime,                                 &
                                        Me%TimeSerie(iTimeSerie)%DT,                &
                                        STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ReadTimeSeriesTimes - ModuleTimeSerie - ERR11a'
            endif

            !Assumes that the last instant of the output is at the end of the simulation
            call GetComputeTimeLimits(Me%ObjTime,                                   &
                                      BeginTime = DummyTime,                        &
                                      EndTime   = Me%TimeSerie(iTimeSerie)%EndOutPut,&
                                      STAT      = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'ReadTimeSeriesTimes - ModuleTimeSerie - ERR11b'

            !Computes the total number of outputs
            Me%TimeSerie(iTimeSerie)%TotalOutPutsNumber                             &
                        = (Me%TimeSerie(iTimeSerie)%EndOutPut                       &
                        -  Me%TimeSerie(iTimeSerie)%BeginOutPut)                    &
                        /  Me%TimeSerie(iTimeSerie)%DT + 1

            !Inits NextOutput
            Me%TimeSerie(iTimeSerie)%NextOutPut =                                   &
                Me%TimeSerie(iTimeSerie)%BeginOutPut + Me%TimeSerie(iTimeSerie)%FirstDT

            !Unlocks the block
            call Block_Unlock(Me%ObjEnterData, ClientNumber)


        enddo

        ReadTimeSeriesTimes = SUCCESS_

    end function ReadTimeSeriesTimes

    !--------------------------------------------------------------------------

    subroutine ReadOnlyOneTimeSerieData(ResultFileName, CoordX, CoordY)

        !Arguments-------------------------------------------------------------
        character(len=*)                    :: ResultFileName
        real, optional                      :: CoordX, CoordY

        !External--------------------------------------------------------------

        type (T_Time)                       :: AuxTime

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, iflag
        integer                             :: FromFile
        integer                             :: ClientNumber
        type (T_Time)                       :: DummyTime
        integer, parameter                  :: OneTimeSerie = 1

        !----------------------------------------------------------------------

        !Gets parameter from the module EnterData
        call GetExtractType(FromFile = FromFile)

        !Allocate only one TimeSerie
        allocate (Me%TimeSerie(OneTimeSerie), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ReadOnlyOneTimeSerieData - ModuleTimeSerie - ERR03'

        nullify (Me%TimeSerie(OneTimeSerie)%TimeSerieData)
        nullify (Me%TimeSerie(OneTimeSerie)%TimeBuffer)
        nullify (Me%TimeSerie(OneTimeSerie)%ResidualValues)
        
        !Stores only one TimeSerie
        Me%NumberOfTimeSeries = OneTimeSerie
        
        
        !Name of the Time Series
        Me%TimeSerie(OneTimeSerie)%FromBlockFileName = ResultFileName
        
        !Sets CoordON
        if (present(CoordX) .and. present(CoordY)) then
            Me%TimeSerie(OneTimeSerie)%CoordON = .true.
            Me%TimeSerie(OneTimeSerie)%CoordX  = CoordX
            Me%TimeSerie(OneTimeSerie)%CoordY  = CoordY
        endif

        !Searches for the first output time
        call GetData(AuxTime,                                                       &
                     Me%ObjEnterData,                                               &
                     iflag,                                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'FIRST_OUTPUT_TIME',                            &
                     ClientModule ='ModuleTimeSerie',                               &
                     STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ReadOnlyOneTimeSerieData - ModuleTimeSerie - ERR10'

        !If the first output time is not specified assmues the current one
        if (iflag == 1) then
            Me%TimeSerie(OneTimeSerie)%BeginOutPut  = AuxTime
        elseif (iflag == 0) then
            call GetComputeCurrentTime(Me%ObjTime,                                  &
                                       Me%TimeSerie(OneTimeSerie)%BeginOutPut,      &
                                       STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'ReadOnlyOneTimeSerieData - ModuleTimeSerie - ERR10a'
        else
            write(*,*)"Error specifing first instant of the a time serie"
            write(*,*)"Assuming current time"
            write(*,*)"ReadOnlyOneTimeSerieData - ModuleTimeSerie - WRN01."
            call GetComputeCurrentTime(Me%ObjTime,                                  &
                                       Me%TimeSerie(OneTimeSerie)%BeginOutPut,      &
                                       STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'ReadOnlyOneTimeSerieData - ModuleTimeSerie - ERR10b'
        endif

        !Searches for the DT
        call GetData(Me%TimeSerie(OneTimeSerie)%DT,                                 &
                     Me%ObjEnterData,                                               &
                     iflag,                                                         &
                     SearchType   = FromFile,                                       &
                     keyword      = 'DT_OUTPUT_TIME'      ,                         &
                     ClientModule ='ModuleTimeSerie',                               &
                     STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ReadOnlyOneTimeSerieData - ModuleTimeSerie - ERR11'

        if (iflag == 0) then
            call GetComputeTimeStep(Me%ObjTime,                                     &
                                    Me%TimeSerie(OneTimeSerie)%DT,                  &
                                    STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            & 
                stop 'ReadOnlyOneTimeSerieData - ModuleTimeSerie - ERR11a'
        endif

        !Assumes that the last instant of the output is at the end of the simulation
        call GetComputeTimeLimits(Me%ObjTime,                                       &
                                  BeginTime = DummyTime,                            &
                                  EndTime   = Me%TimeSerie(OneTimeSerie)%EndOutPut, &
                                  STAT      = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ReadOnlyOneTimeSerieData - ModuleTimeSerie - ERR11b'

        !Computes the total number of outputs
        Me%TimeSerie(OneTimeSerie)%TotalOutPutsNumber                               &
                    = (Me%TimeSerie(OneTimeSerie)%EndOutPut                         &
                    -  Me%TimeSerie(OneTimeSerie)%BeginOutPut)                      &
                    /  Me%TimeSerie(OneTimeSerie)%DT + 1

        !Inits NextOutput
        Me%TimeSerie(OneTimeSerie)%NextOutPut =                                     &
            Me%TimeSerie(OneTimeSerie)%BeginOutPut

        !Unlocks the block
        call Block_Unlock(Me%ObjEnterData, ClientNumber)

    end subroutine ReadOnlyOneTimeSerieData

    !--------------------------------------------------------------------------

    subroutine VerifyTimeSerieLocation(WaterPoints3D, WaterPoints2D, WaterPoints1D)

        !Arguments-------------------------------------------------------------
        integer, dimension(:, :, :), pointer,optional    :: WaterPoints3D
        integer, dimension(:, :   ), pointer,optional    :: WaterPoints2D
        integer, dimension(:      ), pointer,optional    :: WaterPoints1D

        !Local-----------------------------------------------------------------
        integer                                 :: i, j, k
        integer                                 :: iTimeSerie
        logical                                 :: AllTestsPassed, InBounds

        !----------------------------------------------------------------------

        AllTestsPassed = .true.
cd1:    if(present(WaterPoints3D)) then

            do iTimeSerie = 1, Me%NumberOfTimeSeries
            
                if (Me%TimeSerie(iTimeSerie)%IgnoreON) cycle
            
                if (Me%TimeSerie(iTimeSerie)%CoordON .or.                               &
                    Me%TimeSerie(iTimeSerie)%DepthON) cycle

                i = Me%TimeSerie(iTimeSerie)%LocalizationI
                j = Me%TimeSerie(iTimeSerie)%LocalizationJ
                k = Me%TimeSerie(iTimeSerie)%LocalizationK
                
                InBounds = .true.

                !Checks the bounds of the matrix
                if (i < lbound(WaterPoints3D, dim = 1) .or. &
                    i > ubound(WaterPoints3D, dim = 1) .or. &
                    j < lbound(WaterPoints3D, dim = 2) .or. &
                    j > ubound(WaterPoints3D, dim = 2) .or. &
                    k < lbound(WaterPoints3D, dim = 3) .or. &
                    k > ubound(WaterPoints3D, dim = 3)) then
                    write(*,*)"Out of the DOMAIN [i, j, k]", i, j, k
                    InBounds = .false.
                    AllTestsPassed = .false.
                endif

                !Checks the Waterpoints
                if (InBounds) then
                    if (WaterPoints3D(i, j, k) /= 1) then
                        write(*,*)"No WATER POINT [i, j, k]", i, j, k
                        AllTestsPassed = .false.
                    endif
                endif
            enddo

        else if (present(WaterPoints2D)) then cd1

            do iTimeSerie = 1, Me%NumberOfTimeSeries

                if (Me%TimeSerie(iTimeSerie)%CoordON) cycle

                InBounds = .true.

                i = Me%TimeSerie(iTimeSerie)%LocalizationI
                j = Me%TimeSerie(iTimeSerie)%LocalizationJ

                !Checks the bounds of the matrix
                if (i < lbound(WaterPoints2D, dim = 1) .or. &
                    i > ubound(WaterPoints2D, dim = 1) .or. &
                    j < lbound(WaterPoints2D, dim = 2) .or. &
                    j > ubound(WaterPoints2D, dim = 2)) then
                    write(*,*)"Out of the DOMAIN [i, j]", i, j
                    InBounds = .false.
                     AllTestsPassed = .false.
                endif

                !Checks the Waterpoints
                if (InBounds) then
                    if (WaterPoints2D(i, j) /= 1) then
                        write(*,*)"No WATER POINT [i, j]", i, j
                        AllTestsPassed = .false.
                    endif
                endif
            enddo

        else if (present(WaterPoints1D)) then cd1


            do iTimeSerie = 1, Me%NumberOfTimeSeries

                i = Me%TimeSerie(iTimeSerie)%LocalizationI

                !Checks the bounds of the matrix
                if (i < lbound(WaterPoints1D, dim = 1) .or. &
                    i > ubound(WaterPoints1D, dim = 1)) then
                    write(*,*)"Out of the DOMAIN [i]", i
                    InBounds = .false.
                    AllTestsPassed = .false.
                endif

                !Checks the Waterpoints
                if (InBounds) then
                    if (WaterPoints1D(i) /= 1) then
                        write(*,*)"No WATER POINT [i]", i
                        AllTestsPassed = .false.
                    endif
                endif
            enddo

        endif cd1

        if (.not. AllTestsPassed) then
            write (*,*)'Review the Time Series Location Data File'
            stop 'VerifyTimeSerieLocation - ModuleTimeSerie - ERR01'
        endif

    end subroutine VerifyTimeSerieLocation

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine AllocateTimeSerieBuffer

        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                 :: iTimeSerie
        integer                                 :: nProperties, nOutputs
        integer                                 :: STAT_CALL
        real, parameter                         :: AuxTypeReal = 8

        !----------------------------------------------------------------------

        do iTimeSerie = 1, Me%NumberOfTimeSeries

            if (Me%TimeSerie(iTimeSerie)%IgnoreON) cycle

            nProperties = Me%NumberOfProperties
            nOutputs    = Me%TimeSerie(iTimeSerie)%TotalOutPutsNumber

            !Calculates the size of the buffer
            if (nProperties * nOutputs * AuxTypeReal > Me%MaxBufferSize) then
                Me%TimeSerie(iTimeSerie)%BufferSize  = int(Me%MaxBufferSize / (nProperties * AuxTypeReal))
                Me%TimeSerie(iTimeSerie)%BufferCount = 0
            else
                Me%TimeSerie(iTimeSerie)%BufferSize  = nOutputs
                Me%TimeSerie(iTimeSerie)%BufferCount = 0
            endif

            !Allocates the TimeSerie Data Buffer
            allocate(Me%TimeSerie(iTimeSerie)%TimeSerieData(nProperties,            &
                     Me%TimeSerie(iTimeSerie)%BufferSize), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'AllocateTimeSerieBuffer - ModuleTimeSerie - ERR01'
            Me%TimeSerie(iTimeSerie)%TimeSerieData = null_real

            !Allocates the TimeSerie Time Buffer
            allocate(Me%TimeSerie(iTimeSerie)%TimeBuffer(                           &
                     Me%TimeSerie(iTimeSerie)%BufferSize), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'AllocateTimeSerieBuffer - ModuleTimeSerie - ERR02'

            if(Me%ComputeResidual)then
            
                !Allocates the TimeSerie Residual Values Buffer
                allocate(Me%TimeSerie(iTimeSerie)%ResidualValues(nProperties),          &
                         STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'AllocateTimeSerieBuffer - ModuleTimeSerie - ERR03'
                Me%TimeSerie(iTimeSerie)%ResidualValues = 0.

                Me%TimeSerie(iTimeSerie)%ResidualTime   = 0.

            end if

        enddo

    end subroutine AllocateTimeSerieBuffer

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine OpenTimeSerieFiles(Extension, PropertyList, ResultFileName,              &
                                  Instance, HavePath, Comment)

        !Arguments-------------------------------------------------------------

        character(len=*)                       , intent(IN) :: Extension
        character(len=*), dimension(:), pointer             :: PropertyList
        character(len=*), intent(IN), optional              :: ResultFileName
        character(len=*), intent(IN), optional              :: Instance
        logical         , intent(IN), optional              :: HavePath
        character(len=*), intent(IN), optional              :: Comment
 
        !Local-----------------------------------------------------------------

        integer                                 :: iTimeSerie, PropNumber
        integer                                 :: STAT_CALL, iP, unit, nProp
        character(len=3)                        :: AuxI, AuxJ, AuxK, Aux
        character(len=PathLength)               :: RootPath, FileName
        type (T_Time)                           :: CurrentTime
        integer                                 :: iCh, i
        logical                                 :: HavePath_
        character(len=StringLength)             :: Action_


        !----------------------------------------------------------------------

        PropNumber = size (PropertyList, DIM = 1)

        if (PropNumber > 1000) stop 'OpenTimeSerieFiles - ModuleTimeSerie - ERR10'
        
        if (present(HavePath)) then
            HavePath_ = HavePath
        else
            HavePath_ = .false. 
        endif

        if (.not. HavePath_) then
            !Gets the root path from the file nomfich.dat
            call ReadFileName("ROOT_SRT", RootPath, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                call ReadFileName("ROOT", RootPath, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    call ReadFileName("RAIZ", RootPath, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) RootPath = ' '
                endif
            endif
        else
            RootPath=" "
        endif
        
d1:     do iTimeSerie = 1, Me%NumberOfTimeSeries

            if (Me%TimeSerie(iTimeSerie)%IgnoreON) cycle

cd1 :       if (present(ResultFileName)) then
                !Constructs the name of the file
                Me%TimeSerie(iTimeSerie)%FileName = trim(adjustl(RootPath      ))//      &
                                                    trim(adjustl(ResultFileName))//"."// &
                                                    trim(adjustl(Extension     ))
            
            elseif(Me%TimeSerie(iTimeSerie)%FromBlockFileName .ne. null_str)then

            
                if (Me%ReplacePathON) RootPath = trim(Me%ReplacePath)

                !Constructs the name of the file
                Me%TimeSerie(iTimeSerie)%FileName = trim(adjustl(RootPath      ))//      &
                        trim(adjustl(Me%TimeSerie(iTimeSerie)%FromBlockFileName))

                if (Me%ReplacePathON) then
                    Me%TimeSerie(iTimeSerie)%FileName = trim(Me%TimeSerie(iTimeSerie)%FileName)//&
                                                        "_"//trim(Me%ModelName)//"."//           &
                                                        trim(adjustl(Extension     ))
                else
                    Me%TimeSerie(iTimeSerie)%FileName = trim(Me%TimeSerie(iTimeSerie)%FileName)//&
                                                        "."//trim(adjustl(Extension))               
                endif            
            else
                !Constructs the name of the file
                AuxI = "   "
                AuxJ = "   "
                AuxK = "   "
                write(AuxI,'(i3)')Me%TimeSerie(iTimeSerie)%LocalizationI
                write(AuxJ,'(i3)')Me%TimeSerie(iTimeSerie)%LocalizationJ
                write(AuxK,'(i3)')Me%TimeSerie(iTimeSerie)%LocalizationK
                Me%TimeSerie(iTimeSerie)%FileName     = trim(adjustl(RootPath))
                if(present(Instance)) then
                    Me%TimeSerie(iTimeSerie)%FileName = trim(adjustl(Me%TimeSerie(iTimeSerie)%FileName)) //&
                                                        trim(adjustl(Instance ))//"_"
                endif                                        
                if (Me%TimeSerie1D .or. Me%TimeSerie2D .or. Me%TimeSerie3D) then
                    Me%TimeSerie(iTimeSerie)%FileName = trim(adjustl(Me%TimeSerie(iTimeSerie)%FileName)) // &
                                                        trim(adjustl(AuxI))
                endif
                if (Me%TimeSerie2D .or. Me%TimeSerie3D) then
                    Me%TimeSerie(iTimeSerie)%FileName = trim(adjustl(Me%TimeSerie(iTimeSerie)%FileName)) //"_"// &
                                                            trim(adjustl(AuxJ))
                endif
                if (Me%TimeSerie3D) then
                    Me%TimeSerie(iTimeSerie)%FileName = trim(adjustl(Me%TimeSerie(iTimeSerie)%FileName)) //"_"// &
                                                            trim(adjustl(AuxK))
                endif
                Me%TimeSerie(iTimeSerie)%FileName     = trim(adjustl(Me%TimeSerie(iTimeSerie)%FileName)) //"."// &
                                                        trim(adjustl(Extension))

            end if cd1
            !Opens the file
            call UnitsManager(Me%TimeSerie(iTimeSerie)%UnitNumber, OPEN_FILE,       &
                              STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) then
                write(*,*) 'Error opening time series file ',                       &
                            trim(Me%TimeSerie(iTimeSerie)%FileName)
                stop 'OpenTimeSerieFiles - ModuleTimeSerie - ERR20'
            endif
            
            i = 0
            
            FileName = Me%TimeSerie(iTimeSerie)%FileName
            
            do 

                open (file   = FileName,                                                &
                      unit   = Me%TimeSerie(iTimeSerie)%UnitNumber,                     & 
                      status = "unknown",                                               &
                      form   = "formatted", IOSTAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) then
                    write(*,*) 'Error opening time series file ',                       &
                                trim(Me%TimeSerie(iTimeSerie)%FileName)
                    stop 'OpenTimeSerieFiles - ModuleTimeSerie - ERR50'
                endif            
            
                inquire (file      = FileName,                                          &
                         ACTION    = Action_,                                           & 
                         IOSTAT    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then                     
                    write(*,*) 'Error inquiring if time series file has already opened', &
                                trim(Me%TimeSerie(iTimeSerie)%FileName)
                    stop 'OpenTimeSerieFiles - ModuleTimeSerie - ERR30'
                endif
                
                i = i + 1
                
                if (i >100) then
                    stop 'OpenTimeSerieFiles - ModuleTimeSerie - ERR40'
                endif                    
                
                if (Action_ == 'READ') then
                    
                    if (i == 1) then
                        write(*,*) 'Warnning file already opened'
                        write(*,*) 'FileName =',Trim(Me%TimeSerie(iTimeSerie)%FileName)
                    endif
                    
                    write(Aux,'(I3)') i

                    FileName = Trim(Me%TimeSerie(iTimeSerie)%FileName)//"_copy_"//trim(adjustl(Aux))
                    
                    if (i ==1) then
                        write(*,*) 'NewFileName =',Trim(Me%TimeSerie(iTimeSerie)%FileName)
                        write(*,*) 'OpenTimeSerieFiles - ModuleTimeSerie - WRN10'
                    endif
                    
                    !Close the file
                    call UnitsManager(Me%TimeSerie(iTimeSerie)%UnitNumber, CLOSE_FILE,  &
                                      STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) then
                        write(*,*) 'Error closing time series file ',                   &
                                    trim(Me%TimeSerie(iTimeSerie)%FileName)
                        stop 'OpenTimeSerieFiles - ModuleTimeSerie - ERR50'
                    endif                    
                    
                else
                    exit
                endif                    

            enddo

            
            !Writes the header
            unit  = Me%TimeSerie(iTimeSerie)%UnitNumber
            nProp = Me%NumberOfProperties
            call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT = STAT_CALL)

            if (present(Comment)) then
                call WriteDataLine(unit, "!"//trim(Comment))
            endif

            call WriteDataLine(unit, "Time Serie Results File")

            call WriteDataLine(unit, "NAME", Me%TimeSerie(iTimeSerie)%FromBlockFileName)

            if (Me%TimeSerie1D .or. Me%TimeSerie2D .or. Me%TimeSerie3D)& 
            call WriteDataLine(unit, "LOCALIZATION_I", Me%TimeSerie(iTimeSerie)%LocalizationI)
            if (Me%TimeSerie2D .or. Me%TimeSerie3D)& 
            call WriteDataLine(unit, "LOCALIZATION_J", Me%TimeSerie(iTimeSerie)%LocalizationJ)
            if (Me%TimeSerie3D)& 
            call WriteDataLine(unit, "LOCALIZATION_K", Me%TimeSerie(iTimeSerie)%LocalizationK)
            call WriteDataLine(unit, 'SERIE_INITIAL_DATA', CurrentTime)
            call WriteDataLine(unit, 'TIME_UNITS', 'SECONDS')

            if (Me%ModelNameON) then

                call WriteDataLine(unit, 'MODEL_DOMAIN', trim(Me%ModelName))

            endif


            if (Me%TimeSerie(iTimeSerie)%CoordON) then
                Write(unit,*) "COORD_X    :", Me%TimeSerie(iTimeSerie)%CoordX
                Write(unit,*) "COORD_Y    :", Me%TimeSerie(iTimeSerie)%CoordY
            endif


            if (Me%TimeSerie(iTimeSerie)%DepthON) then
                Write(unit,*) "DEPTH_LEVEL:", Me%TimeSerie(iTimeSerie)%DepthLevel
            endif

            !Changes white spaces to underscores in the property list (Turns Input to Excel more easy)
            do iP = 1, nProp
                do iCh = 1, len_trim(PropertyList(iP))
                    if (PropertyList(iP)(iCh:iCh) == ' ') then
                        PropertyList(iP)(iCh:iCh) =  '_'
                    endif
                enddo
            enddo

            if (Me%UseTabulatedData) then
                write(unit, fmt=1000)(trim(PropertyList(iP)), iP = 1, nProp)
            else
                write(unit, fmt=2000)(trim(PropertyList(iP)), iP = 1, nProp)
            endif
            
            call WriteDataLine(unit, block_begin)

            if(Me%ComputeResidual)then
                !Stores LastResidual Calculation
                Me%TimeSerie(iTimeSerie)%LastResidual = CurrentTime
            end if

       enddo d1

        1000 format(1x,"     Seconds   YY  MM  DD  hh  mm       ss", 2x, 1000(1x, A43))
        2000 format(1x,"     Seconds   YY  MM  DD  hh  mm       ss", 2x, 1000(1x, A))
!        1000 format(1x," Seconds    YY   MM  DD  HH  MM  SS", 6x, 1x, 1000(1x, A12))

    end subroutine OpenTimeSerieFiles

    !--------------------------------------------------------------------------

    subroutine StartTimeSerieInput(TimeSerieID, TimeSerieDataFile, ObjTime,              &
                                   ColumnsToRead, DT_TimeSerie, CheckDates,              &
                                   ReadAll, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        character(len=*)                            :: TimeSerieDataFile
        integer, optional                           :: ObjTime
        integer, optional, pointer, dimension (:)   :: ColumnsToRead
        real   , optional, intent(IN )              :: DT_TimeSerie
        logical, optional, intent(IN )              :: CheckDates
        logical, optional, intent(IN )              :: ReadAll        
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        real   ,  dimension (:,:), allocatable      :: AuxMatrix
        integer                                     :: STAT_CALL
        integer                                     :: ready_ , STAT_        
        integer                                     :: FromFile, FromBlock
        integer                                     :: flag, CurrentDataValue
        integer                                     :: i, iLine
        integer                                     :: StartLine, EndLine, ClientNumber
        logical                                     :: BlockFound, DoCheck
        real(8), dimension(:), allocatable          :: BufferLine, BufferLineNext
        real                                        :: Conversion_Seconds
        type (T_Time)                               :: StartTimeModel, EndTimeModel
        type (T_Time)                               :: StartTimeTimeSerie, EndTimeTimeSerie
        type (T_Time)                               :: CurrentTimeSerie, NextValueTimeSerie
        type (T_Time)                               :: NextTimeSerie, PreviousTimeSerie
        real                                        :: DT
        real                                        :: Year, Month, Day, Hour, Minute, Second
        logical                                     :: ReadAll_        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mTimeSerie_)) then
            nullify (FirstTimeSerie)
            call RegisterModule (mTimeSerie_) 
        endif
   
        call Ready(TimeSerieID, ready_)    

if0 :   if (ready_ .EQ. OFF_ERR_) then

            !Allocates Instance
            call AllocateInstance
  
            nullify (Me%DataMatrix )
            nullify (Me%ColumnsRead)
            nullify (Me%FileColumns)
            nullify (Me%TimeSerie  )
  
            !Associates module Time
            if (present(ObjTime)) then
                Me%ObjTime = AssociateInstance   (mTIME_, ObjTime)
                call GetComputeTimeLimits(Me%ObjTime,                                   &
                                          BeginTime = StartTimeModel,                   &
                                          EndTime   = EndTimeModel, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR01'
            endif
            
            if (present(ReadAll)) then
                ReadAll_ = ReadAll
            else
                ReadAll_ = .false.
            endif


            !Constructs EnterData
            call ConstructEnterData(Me%ObjEnterData, TimeSerieDataFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR02'

            !Get Extract Types
            call GetExtractType    (FromFile  = FromFile )
            call GetExtractType    (FromBlock = FromBlock)

            !Gets the time units
            call GetData(Me%CharTimeUnits,                                          &
                         Me%ObjEnterData,                                           &
                         flag,                                                      &
                         SearchType   = FromFile,                                   &
                         keyword      ='TIME_UNITS',                                &
                         ClientModule ='ModuleTimeSerie',                           &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR03'
            if (flag      /= 1       ) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR04'


            !Gets start time
            call GetData(Me%InitialData,                                            &
                         Me%ObjEnterData,                                           &
                         flag,                                                      &
                         SearchType   = FromFile,                                   &
                         keyword      ='SERIE_INITIAL_DATA',                        &
                         ClientModule ='ModuleTimeSerie',                           &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR05'
            if (flag      /= 1       ) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR06'

            !Time cycle
            call GetData(Me%TimeCycle,                                              &
                         Me%ObjEnterData,                                           &
                         flag,                                                      &
                         SearchType   = FromFile,                                   &
                         default      = .false.,                                    &
                         keyword      ='TIME_CYCLE',                                &
                         ClientModule ='ModuleTimeSerie',                           &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR07'

            
            !Extracts the block from the file
            call ExtractBlockFromBuffer(Me%ObjEnterData,                            &
                                        ClientNumber    = ClientNumber,             &
                                        block_begin     = block_begin,              &
                                        block_end       = block_end,                &
                                        BlockFound      = BlockFound,               &
                                        STAT            = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR08'
            if (.not. BlockFound)      stop 'StartTimeSerieInput - ModuleTimeSerie - ERR09'


            !Gets the number of values in the block
            call GetBlockSize(Me%ObjEnterData, ClientNumber, StartLine, EndLine)

            if(StartLine-1 > 0)then
                call GetFullBufferLine(Me%ObjEnterData, StartLine-1, Me%Header, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR09.1'
            endif

            !Number of data values (Block size contains BeginBlock and EndBlock)
            Me%DataValues = EndLine - StartLine - 1


            !Converts the time units of Me%DataMatrix to seconds
            select case (trim(Me%CharTimeUnits))
                case ('SECONDS')
                    Conversion_Seconds = 1.
                    if (Me%TimeCycle) then
                        write(*,*)'Time Units SECONDS cant be an a cycle'
                        stop 'StartTimeSerieInput - ModuleTimeSerie - ERR10'
                    endif

                case ('MINUTES')
                    Conversion_Seconds = 60.
                    if (Me%TimeCycle) then
                        write(*,*)'Time Units MINUTES cant be an a cycle'
                        stop 'StartTimeSerieInput - ModuleTimeSerie - ERR11'
                    endif

                case ('HOURS')
                    Conversion_Seconds = 3600.
                    if (Me%TimeCycle .and. Me%DataValues /= 24) then
                        write(*,*)'For HOUR cycle you must supply 24 values'
                        stop 'StartTimeSerieInput - ModuleTimeSerie - ERR12'
                    endif

                case ('DAYS')
                    Conversion_Seconds = 86400.
                    if (Me%TimeCycle .and. Me%DataValues /= 366) then
                        write(*,*)'For DAY cycle you must supply 366 values - Julian days'
                        stop 'StartTimeSerieInput - ModuleTimeSerie - ERR13'
                    endif

                case ('MONTHS')
                    Conversion_Seconds = 86400. * 31.
                    if (Me%TimeCycle .and. Me%DataValues /= 12) then
                        write(*,*)'For MONTHS cycle you must supply 12 values'
                        stop 'StartTimeSerieInput - ModuleTimeSerie - ERR14'
                    endif

                case default
                    stop 'StartTimeSerieInput - ModuleTimeSerie - ERR15'
            end select



            !Allocates Buffer Line
            allocate(BufferLine(1:MaxColumns))

            !Get a new line
            call GetData(BufferLine,                                                &
                         Me%ObjEnterData,                                           &
                         flag,                                                      &
                         Buffer_Line = StartLine+1,                                 &
                         STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= SIZE_ERR_) then
                stop 'StartTimeSerieInput - ModuleTimeSerie - ERR16'
            endif
            if (flag       < 1         ) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR17'
            if (flag       > MaxColumns) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR18'

            !Deallocates Buffer Line
            deallocate(BufferLine)
 
            allocate(Me%FileColumns(flag))

            Me%FileColumns(:)  = Null_Int

            if (present(ColumnsToRead)) then
                                         ! Time column + Column values to read       
                Me%DataColumns = 1 + Size(ColumnsToRead)

                allocate(Me%ColumnsRead(Me%DataColumns))

                ! Time
                Me%ColumnsRead(1) = 1

                Me%ColumnsRead(2:Me%DataColumns) = ColumnsToRead(:)
            else

                Me%DataColumns = flag

                allocate(Me%ColumnsRead(Me%DataColumns))

                do i = 1, flag

                  Me%ColumnsRead(i) = i

                enddo

            endif

            do i = 1, Me%DataColumns

               Me%FileColumns(Me%ColumnsRead(i)) = i

            enddo
            
            allocate(BufferLine        (1:flag))
            allocate(BufferLineNext    (1:flag))

           !Allocates the auxiliar matrix to hold the data
            allocate(AuxMatrix(1:Me%DataValues, 1:Me%DataColumns))

            CurrentDataValue = 0

            call null_Time(PreviousTimeSerie)

            NextValueTimeSerie = Me%InitialData

            do iLine = StartLine+1, EndLine - 1
                
                !Get a new line
                call GetData(BufferLine,                                            &
                             Me%ObjEnterData,                                       &
                             flag,                                                  &
                             Buffer_Line = iLine,                                   &
                             STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= SIZE_ERR_) then
                    stop 'StartTimeSerieInput - ModuleTimeSerie - ERR19'
                endif
                if (flag       < 1       ) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR20'

                !Conversion to seconds
                BufferLine(1) = BufferLine(1) * Conversion_Seconds

                CurrentTimeSerie = Me%InitialData + BufferLine(1)
                
                if (.not. Me%TimeCycle) then

                    if (iLine + 1 < EndLine) then
                        call GetData(BufferLineNext,                                &
                                     Me%ObjEnterData,                               &
                                     flag,                                          &
                                     Buffer_Line = iLine + 1,                       &
                                     STAT = STAT_CALL) 
                        if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= SIZE_ERR_) then
                            stop 'StartTimeSerieInput - ModuleTimeSerie - ERR21'
                        endif
                        if (flag       < 1       ) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR22'

                        !Conversion to seconds
                        BufferLineNext(1) = BufferLineNext(1) * Conversion_Seconds

                        NextTimeSerie = Me%InitialData + BufferLineNext(1)
                    else
                        BufferLineNext = BufferLine
                        NextTimeSerie  = CurrentTimeSerie
                    endif

                    if (present(ObjTime) .and. .not. ReadAll_) then
                        !Stops reading the time serie data if the values correspond to a time after 
                        ! the end time of the model
                        if (CurrentTimeSerie  > EndTimeModel .and.                      &
                            PreviousTimeSerie > EndTimeModel) then

                            exit

                        endif

                        !Begin reading the time serie data if the values correspond to a time after
                        ! the begin time of the model
                        if (CurrentTimeSerie < StartTimeModel .and.                     &
                            NextTimeSerie    < StartTimeModel) then

                            cycle

                        endif
                    endif

                    !Reads the time serie data 
                    if (present(DT_TimeSerie))                  then
                    if (CurrentTimeSerie >= NextValueTimeSerie) then

                        NextValueTimeSerie = NextValueTimeSerie + DT_TimeSerie

                    else

                        cycle

                    endif 
                    endif


                endif

                
                CurrentDataValue = CurrentDataValue + 1

                !Stores the data in the matrix
                do i = 1, Me%DataColumns
                    AuxMatrix(CurrentDataValue, i) = BufferLine(Me%ColumnsRead(i))
                enddo
                
                PreviousTimeSerie = CurrentTimeSerie
            enddo

            deallocate(BufferLine)
           
            if (.not. Me%TimeCycle) deallocate(BufferLineNext)

            !Verifies if he found some useful information
            if (CurrentDataValue == 0) then
                write (*, *)'The TimeSerie data contains no data for the simulation period.'
                write (*, *)'File : ',trim(adjustl(TimeSerieDataFile))
                stop 'StartTimeSerieInput - ModuleTimeSerie - ERR23'
            endif
    
            !Allocates the matrix to hold the data
            allocate(Me%DataMatrix(1:CurrentDataValue,                              &
                                             1:Me%DataColumns))

            Me%DataMatrix(1:CurrentDataValue, 1:Me%DataColumns)=                    &
                          AuxMatrix(1:CurrentDataValue, 1:Me%DataColumns)


            Me%DataValues = CurrentDataValue

            deallocate(AuxMatrix)

            !Unlocks Block
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR24'

            !Kills EnterData
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR25'



            !If the time serie isnt a cycle, verify if the time limits are within the model
            !run time limits
            if (.not. Me%TimeCycle .and. present(ObjTime)) then

                if (present(CheckDates)) then
                    DoCheck = CheckDates
                else
                    DoCheck = .true.    
                endif

                if (DoCheck) then

                    call GetComputeTimeStep (Me%ObjTime, DT, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'StartTimeSerieInput - ModuleTimeSerie - ERR26'

                    StartTimeTimeSerie = Me%InitialData + Me%DataMatrix(1, 1)
                    EndTimeTimeSerie   = Me%InitialData + Me%DataMatrix(Me%DataValues, 1)
                

                    if (StartTimeModel + DT .lt. StartTimeTimeSerie .or.  &
                        EndTimeModel        .gt. EndTimeTimeSerie) then

                        write(*,*)'Time serie defined between:'

                        call ExtractDate(StartTimeTimeSerie, Year, Month, Day, Hour, Minute, Second)
                        write(*,fmt=100)Year, Month, Day, Hour, Minute, Second

                        call ExtractDate(EndTimeTimeSerie,   Year, Month, Day, Hour, Minute, Second)
                        write(*,fmt=100)Year, Month, Day, Hour, Minute, Second

                        write(*,*)'Cant interpolate for:'
                        call ExtractDate(StartTimeModel + DT, Year, Month, Day, Hour, Minute, Second)
                        write(*,fmt=100)Year, Month, Day, Hour, Minute, Second

                        call ExtractDate(EndTimeModel, Year, Month, Day, Hour, Minute, Second)
                        write(*,fmt=100)Year, Month, Day, Hour, Minute, Second

                        write(*,*)'File                :', trim(adjustl(TimeSerieDataFile))

                        stop 'StartTimeSerieInput - ModuleTimeSerie - ERR27'
            
        100 format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)
                    endif
                endif    

            endif


            !Returns ID
            TimeSerieID    = Me%InstanceID


            STAT_ = SUCCESS_

        else
         
            stop 'ModuleTimeSerie - StartTimeSerieInput - ERR28' 

        end if if0


        if (present(STAT)) STAT = STAT_


    end subroutine StartTimeSerieInput

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine WriteTimeSerie1(TimeSerieID, Data1D,     Data2D,     Data3D,    &
                                            Data1D_8,   Data2D_8,   Data3D_8,  &
                                            Data1D_Int, Data2D_Int, Data3D_Int,&
                                            factor, STAT)

        !Arguments-------------------------------------------------------------
        integer                                           :: TimeSerieID
        real,    dimension(:),       optional, pointer    :: Data1D
        real,    dimension(:, :),    optional, pointer    :: Data2D
        real,    dimension(:, :, :), optional, pointer    :: Data3D
        real(8), dimension(:),       optional, pointer    :: Data1D_8
        real(8), dimension(:, :),    optional, pointer    :: Data2D_8
        real(8), dimension(:, :, :), optional, pointer    :: Data3D_8
        integer, dimension(:),       optional, pointer    :: Data1D_Int
        integer, dimension(:, :),    optional, pointer    :: Data2D_Int
        integer, dimension(:, :, :), optional, pointer    :: Data3D_Int
        real,    optional, intent(IN)                     :: factor
        integer, optional, intent(OUT)                    :: STAT


        !Local-----------------------------------------------------------------
        integer                                     :: ready_ , STAT_        
        type (T_Time)                               :: CurrentTime
        integer                                     :: IPC, IBC, iTimeSerie
        integer                                     :: i, j, k
        real                                        :: DT_Residual, DataValue, DTaux, factor_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    

        !Shorten variavel name
        IPC = Me%InternalPropertyCount

cd1 :   if (ready_ .EQ. IDLE_ERR_) then
cd2 :       if (Me%Points) then

                !Actualize the Property Count
                IPC = IPC + 1
                if (IPC > Me%NumberOfProperties) IPC = 1
            
                !Stores the actual compute time
                call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT)

                if (present(factor))  then
                    factor_ = factor
                else
                    factor_ = 1.
                endif


                !Calculates the residual values
                do iTimeSerie = 1, Me%NumberOfTimeSeries
                
                    if (Me%TimeSerie(iTimeSerie)%IgnoreON) cycle

                    i = Me%TimeSerie(iTimeSerie)%LocalizationI
                    j = Me%TimeSerie(iTimeSerie)%LocalizationJ
                    k = Me%TimeSerie(iTimeSerie)%LocalizationK
                    
                    if (Me%TimeSerie(iTimeSerie)%IgnoreON) then
                        DataValue = FillValueReal
                    
                    else

                        if (present(Data1D))     DataValue = Data1D     (i)
                        if (present(Data2D))     DataValue = Data2D     (i, j)
                        if (present(Data3D))     DataValue = Data3D     (i, j, k)
                        if (present(Data1D_8))   DataValue = Data1D_8   (i)
                        if (present(Data2D_8))   DataValue = Data2D_8   (i, j)
                        if (present(Data3D_8))   DataValue = Data3D_8   (i, j, k)
                        if (present(Data1D_Int)) DataValue = Data1D_Int (i)
                        if (present(Data2D_Int)) DataValue = Data2D_Int (i, j)
                        if (present(Data3D_Int)) DataValue = Data3D_Int (i, j, k)

                         DataValue = DataValue * factor_


                        if(Me%ComputeResidual)then

                            DT_Residual = CurrentTime - Me%TimeSerie(iTimeSerie)%LastResidual
                            DTaux       = Me%TimeSerie(iTimeSerie)%ResidualTime + DT_Residual

                            if (DTaux > 0) then

                                !Updates the Residual Values
                                Me%TimeSerie(iTimeSerie)%ResidualValues(IPC) =              &
                                    (Me%TimeSerie(iTimeSerie)%ResidualValues(IPC)           &
                                   * Me%TimeSerie(iTimeSerie)%ResidualTime                  &
                                   + DataValue * DT_Residual)/ DTaux
                        
                                !Updates the times
                                if (IPC == Me%NumberOfProperties) then
                                    Me%TimeSerie(iTimeSerie)%ResidualTime = DTaux
                                    Me%TimeSerie(iTimeSerie)%LastResidual = CurrentTime
                                endif

                            end if

                        end if
                    
                    endif

                    !Stores the data in the buffer
                    if (CurrentTime .ge. Me%TimeSerie(iTimeSerie)%NextOutPut) then
                    
                        !Increments the internal buffer count
                        if (IPC == 1) then
                            Me%TimeSerie(iTimeSerie)%BufferCount =                  &
                                Me%TimeSerie(iTimeSerie)%BufferCount + 1 
                        endif

                        !Shorten Variable
                        IBC = Me%TimeSerie(iTimeSerie)%BufferCount

                        !Stores the current time
                        if (IPC == 1)                                         &
                            Me%TimeSerie(iTimeSerie)%TimeBuffer(IBC) = CurrentTime

                        if (Present(Data1D))                                        &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = DataValue

                        if (Present(Data2D))                                        &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = DataValue

                        if (Present(Data3D))                                        &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = DataValue

                        if (Present(Data1D_8))                                      &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = DataValue

                        if (Present(Data2D_8))                                      &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = DataValue

                        if (Present(Data3D_8))                                      &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = DataValue

                        if (Present(Data1D_Int))                                    &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = DataValue

                        if (Present(Data2D_Int))                                    &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = DataValue

                        if (Present(Data3D_Int))                                    &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = DataValue


                        !Sets next output time
                        if (IPC == Me%NumberOfProperties)                           &
                            !Me%TimeSerie(iTimeSerie)%NextOutPut =                   &
                            !    CurrentTime + Me%TimeSerie(iTimeSerie)%DT
                            Me%TimeSerie(iTimeSerie)%NextOutPut =                   &
                                Me%TimeSerie(iTimeSerie)%NextOutPut + Me%TimeSerie(iTimeSerie)%DT                        
                    endif

                    !Verifies if the buffer is full
                    if ((Me%TimeSerie(iTimeSerie)%BufferCount  ==                   &
                         Me%TimeSerie(iTimeSerie)%BufferSize)  .and.                &
                        (IPC == Me%NumberOfProperties))        then
                        call WriteBufferToFile(Me%TimeSerie(iTimeSerie), Me%NumberOfProperties)
                        Me%TimeSerie(iTimeSerie)%BufferCount = 0
                    endif

                enddo

                STAT_ = SUCCESS_
            else
                STAT_ = UNKNOWN_
            end if cd2
        else               

            STAT_ = ready_

        end if cd1


        Me%InternalPropertyCount = IPC

        if (present(STAT))                                                          &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine WriteTimeSerie1

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine WriteTimeSerieLine(TimeSerieID, DataLine, ExternalCurrentTime, LineStored, CheckTime, NextDT, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        real,    dimension(:), pointer              :: DataLine
        type (T_Time), optional, intent(IN)         :: ExternalCurrentTime
        logical, optional, intent(OUT)              :: LineStored
        logical, optional, intent(IN)               :: CheckTime
        real, optional, intent(IN)                  :: NextDT        
        integer, optional, intent(OUT)              :: STAT
    
        !External--------------------------------------------------------------
        integer                                     :: ready_ , STAT_        

        !Local-----------------------------------------------------------------
        type (T_Time)                               :: CurrentTime
        integer                                     :: IPC, IBC, iTimeSerie
        real                                        :: DT_Residual
        real                                        :: NextDT_
        logical                                     :: CheckTime_          

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    

        !Shorten variavel name
        IPC = Me%InternalPropertyCount

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            if (present(ExternalCurrentTime)) then
                CurrentTime = ExternalCurrentTime
            else                
                !Stores the actual compute time
                call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT)
            endif
            
            if (present(CheckTime)) then
                CheckTime_ = CheckTime
            else
                CheckTime_ = .true.
            endif                    
            
do1:        do iTimeSerie = 1, Me%NumberOfTimeSeries

                if (Me%TimeSerie(iTimeSerie)%IgnoreON) cycle

                if(Me%ComputeResidual)then

                    !Calculates the residual values
                    DT_Residual = CurrentTime - Me%TimeSerie(iTimeSerie)%LastResidual

                    if (DT_Residual + Me%TimeSerie(iTimeSerie)%ResidualTime > 0) then

                        !Updates the Residual Values
                        do IPC = 1, Me%NumberOfProperties
                            Me%TimeSerie(iTimeSerie)%ResidualValues(IPC) =                  &
                                (Me%TimeSerie(iTimeSerie)%ResidualValues(IPC)               &
                               * Me%TimeSerie(iTimeSerie)%ResidualTime                      &
                               + DataLine(IPC) * DT_Residual)                               &
                               / (Me%TimeSerie(iTimeSerie)%ResidualTime + DT_Residual)
                        enddo
                    endif

                    Me%TimeSerie(iTimeSerie)%ResidualTime =                             &
                    Me%TimeSerie(iTimeSerie)%ResidualTime + DT_Residual

                    Me%TimeSerie(iTimeSerie)%LastResidual = CurrentTime

                end if

                !Stores the data in the buffer
cd2:            if ((.not. CheckTime_) .or.                                     &
                    (CurrentTime .ge. Me%TimeSerie(iTimeSerie)%NextOutPut) .or. &
                    (CurrentTime .eq. Me%TimeSerie(iTimeSerie)%EndOutPut)) then
                    
                    if (present(LineStored)) LineStored = .true.

do2:                do IPC = 1, Me%NumberOfProperties

                        if (IPC == 1) then
                            !Increments the internal buffer count
                            Me%TimeSerie(iTimeSerie)%BufferCount =                  &
                                Me%TimeSerie(iTimeSerie)%BufferCount + 1 
                            IBC = Me%TimeSerie(iTimeSerie)%BufferCount
                            
                            !Stores the current time
                            Me%TimeSerie(iTimeSerie)%TimeBuffer(IBC) = CurrentTime
                        endif

                        !Shorten Variable
                        IBC = Me%TimeSerie(iTimeSerie)%BufferCount

                        !Locates the place of the time serie to store
                        Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC,IBC) = DataLine(IPC)
                        
                        !Sets next output time
                        if (IPC == Me%NumberOfProperties) then
                        
                            !Me%TimeSerie(iTimeSerie)%NextOutPut =                   &
                            !    CurrentTime + Me%TimeSerie(iTimeSerie)%DT
                            if (present(NextDT)) then
                                NextDT_ = NextDT
                            else
                                NextDT_ = Me%TimeSerie(iTimeSerie)%DT
                            endif    
                            
                            Me%TimeSerie(iTimeSerie)%NextOutPut =                   &
                                Me%TimeSerie(iTimeSerie)%NextOutPut + NextDT_                        

                        endif
                        
                        !Verifies if the buffer is full
                        if ((Me%TimeSerie(iTimeSerie)%BufferCount  ==               &
                             Me%TimeSerie(iTimeSerie)%BufferSize)  .and.            &
                            (IPC == Me%NumberOfProperties))        then
                            call WriteBufferToFile(Me%TimeSerie(iTimeSerie),Me%NumberOfProperties)
                            Me%TimeSerie(iTimeSerie)%BufferCount = 0
                        endif

                    enddo do2

                else
                
                    if (present(LineStored)) LineStored = .false.
                    
                endif cd2

            enddo do1

            STAT_ = SUCCESS_

        else               

            STAT_ = ready_

        end if cd1


        Me%InternalPropertyCount = IPC

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine WriteTimeSerieLine

    !--------------------------------------------------------------------------
    !This subroutine is implemented so that user decides when to write line
    !e.g. for hourly and daily values that do not need to be linked to timeserie dt
    
    subroutine WriteTimeSerieLineNow(TimeSerieID, DataLine, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        real,    dimension(:), pointer              :: DataLine
        integer, optional, intent(OUT)              :: STAT
    
        !External--------------------------------------------------------------
        integer                                     :: ready_ , STAT_        

        !Local-----------------------------------------------------------------
        type (T_Time)                               :: CurrentTime
        integer                                     :: IPC, IBC, iTimeSerie
        real                                        :: DT_Residual

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    

        !Shorten variavel name
        IPC = Me%InternalPropertyCount

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            !Stores the actual compute time
            call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT)


do1:        do iTimeSerie = 1, Me%NumberOfTimeSeries

                if (Me%TimeSerie(iTimeSerie)%IgnoreON) cycle

                if(Me%ComputeResidual)then

                    !Calculates the residual values
                    DT_Residual = CurrentTime - Me%TimeSerie(iTimeSerie)%LastResidual

                    if (DT_Residual + Me%TimeSerie(iTimeSerie)%ResidualTime > 0) then

                        !Updates the Residual Values
                        do IPC = 1, Me%NumberOfProperties
                            Me%TimeSerie(iTimeSerie)%ResidualValues(IPC) =                  &
                                (Me%TimeSerie(iTimeSerie)%ResidualValues(IPC)               &
                               * Me%TimeSerie(iTimeSerie)%ResidualTime                      &
                               + DataLine(IPC) * DT_Residual)                               &
                               / (Me%TimeSerie(iTimeSerie)%ResidualTime + DT_Residual)
                        enddo
                    endif

                    Me%TimeSerie(iTimeSerie)%ResidualTime =                             &
                    Me%TimeSerie(iTimeSerie)%ResidualTime + DT_Residual

                    Me%TimeSerie(iTimeSerie)%LastResidual = CurrentTime

                end if

                !Always Stores the data in the buffer!!!
!cd2:            if (CurrentTime .ge. Me%TimeSerie(iTimeSerie)%NextOutPut .or.           &
!                    CurrentTime .eq. Me%TimeSerie(iTimeSerie)%EndOutPut) then
                    
do2:                do IPC = 1, Me%NumberOfProperties



                        if (IPC == 1) then
                            !Increments the internal buffer count
                            Me%TimeSerie(iTimeSerie)%BufferCount =                  &
                                Me%TimeSerie(iTimeSerie)%BufferCount + 1 
                            IBC = Me%TimeSerie(iTimeSerie)%BufferCount
                            
                            !Stores the current time
                            Me%TimeSerie(iTimeSerie)%TimeBuffer(IBC) = CurrentTime
                        endif


                        !Shorten Variable
                        IBC = Me%TimeSerie(iTimeSerie)%BufferCount


                        !Locates the place of the time serie to store
                        Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC,IBC) = DataLine(IPC)


                        !Sets next output time
                        if (IPC == Me%NumberOfProperties)                           &
                            !Me%TimeSerie(iTimeSerie)%NextOutPut =                   &
                            !    CurrentTime + Me%TimeSerie(iTimeSerie)%DT
                            Me%TimeSerie(iTimeSerie)%NextOutPut =                   &
                                Me%TimeSerie(iTimeSerie)%NextOutPut + Me%TimeSerie(iTimeSerie)%DT                        

                        !Verifies if the buffer is full
                        if ((Me%TimeSerie(iTimeSerie)%BufferCount  ==               &
                             Me%TimeSerie(iTimeSerie)%BufferSize)  .and.            &
                            (IPC == Me%NumberOfProperties))        then
                            call WriteBufferToFile(Me%TimeSerie(iTimeSerie),Me%NumberOfProperties)
                            Me%TimeSerie(iTimeSerie)%BufferCount = 0
                        endif

                    enddo do2

!                endif cd2

            enddo do1

            STAT_ = SUCCESS_

        else               

            STAT_ = ready_

        end if cd1


        Me%InternalPropertyCount = IPC

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine WriteTimeSerieLineNow

    !--------------------------------------------------------------------------    
    
    subroutine WriteSpecificTimeSerieLine(TimeSerieID, iTimeSerie, DataLine, CheckTime, NextDT, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        integer                                     :: iTimeSerie
        real, dimension(:), pointer                 :: DataLine        
        logical, optional, intent(IN)               :: CheckTime
        real, optional, intent(IN)                  :: NextDT
        integer, optional, intent(OUT)              :: STAT
    
        !External--------------------------------------------------------------
        integer                                     :: ready_ , STAT_        

        !Local-----------------------------------------------------------------
        type (T_Time)                               :: CurrentTime
        integer                                     :: IPC, IBC
        real                                        :: DT_Residual
        real                                        :: NextDT_
        logical                                     :: CheckTime_        

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    

        !Shorten variavel name
        IPC = Me%InternalPropertyCount

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            if (present(CheckTime)) then
                CheckTime_ = CheckTime
            else
                CheckTime_ = .true.
            endif
            
            if (present(NextDT)) then
                NextDT_ = NextDT
            else
                NextDT_ = Me%TimeSerie(iTimeSerie)%DT
            endif

            !Stores the actual compute time
            call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT)

            if(Me%ComputeResidual)then

                !Calculates the residual values
                DT_Residual = CurrentTime - Me%TimeSerie(iTimeSerie)%LastResidual

                !Updates the Residual Values
                do IPC = 1, Me%NumberOfProperties
                    Me%TimeSerie(iTimeSerie)%ResidualValues(IPC) =                  &
                        (Me%TimeSerie(iTimeSerie)%ResidualValues(IPC)               &
                       * Me%TimeSerie(iTimeSerie)%ResidualTime                      &
                       + DataLine(IPC) * DT_Residual)                               &
                       / (Me%TimeSerie(iTimeSerie)%ResidualTime + DT_Residual)
                enddo

                Me%TimeSerie(iTimeSerie)%ResidualTime =                             &
                Me%TimeSerie(iTimeSerie)%ResidualTime + DT_Residual

                Me%TimeSerie(iTimeSerie)%LastResidual = CurrentTime

            end if

            !Stores the data in the buffer
cd2:        if ((.not. CheckTime_) .or.                                     &
                (CurrentTime .ge. Me%TimeSerie(iTimeSerie)%NextOutPut) .or. &
                (CurrentTime .eq. Me%TimeSerie(iTimeSerie)%EndOutPut)) then
                
do2:            do IPC = 1, Me%NumberOfProperties



                    if (IPC == 1) then
                        !Increments the internal buffer count
                        Me%TimeSerie(iTimeSerie)%BufferCount =                  &
                            Me%TimeSerie(iTimeSerie)%BufferCount + 1 
                        IBC = Me%TimeSerie(iTimeSerie)%BufferCount
                        
                        !Stores the current time
                        Me%TimeSerie(iTimeSerie)%TimeBuffer(IBC) = CurrentTime
                    endif


                    !Shorten Variable
                    IBC = Me%TimeSerie(iTimeSerie)%BufferCount


                    !Locates the place of the time serie to store
                    Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC,IBC) = DataLine(IPC)


                    !Sets next output time
                    if (IPC == Me%NumberOfProperties)                           &
                        !Me%TimeSerie(iTimeSerie)%NextOutPut =                   &
                        !    CurrentTime + Me%TimeSerie(iTimeSerie)%DT
                        Me%TimeSerie(iTimeSerie)%NextOutPut =                   &
                            Me%TimeSerie(iTimeSerie)%NextOutPut + NextDT_                    

                    !Verifies if the buffer is full
                    if ((Me%TimeSerie(iTimeSerie)%BufferCount  ==               &
                         Me%TimeSerie(iTimeSerie)%BufferSize)  .and.            &
                        (IPC == Me%NumberOfProperties))        then
                        call WriteBufferToFile(Me%TimeSerie(iTimeSerie),Me%NumberOfProperties)
                        Me%TimeSerie(iTimeSerie)%BufferCount = 0
                    endif

                enddo do2

            endif cd2

            STAT_ = SUCCESS_

        else               

            STAT_ = ready_

        end if cd1


        Me%InternalPropertyCount = IPC

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine WriteSpecificTimeSerieLine
    
    !--------------------------------------------------------------------------


    subroutine WriteTimeSerie3(TimeSerieID, Data2D, Data3D, Data2D_8, Data3D_8, &
                               Data2D_Int, Data3D_Int, Icell, Jcell, Kcell, factor, STAT)

        !Arguments-------------------------------------------------------------
        integer                                           :: TimeSerieID
        real,    dimension(:, :),    optional, pointer    :: Data2D
        real,    dimension(:, :, :), optional, pointer    :: Data3D
        real(8), dimension(:, :),    optional, pointer    :: Data2D_8
        real(8), dimension(:, :, :), optional, pointer    :: Data3D_8
        integer, dimension(:, :),    optional, pointer    :: Data2D_Int
        integer, dimension(:, :, :), optional, pointer    :: Data3D_Int
        integer, dimension(:),                 intent(IN) :: Icell, Jcell
        integer, dimension(:),        optional,intent(IN) :: Kcell
        real,    optional, intent(IN )                    :: factor
        integer, optional,                     intent(OUT):: STAT


        !Local-----------------------------------------------------------------
        integer                                     :: ready_ , STAT_        
        type (T_Time)                               :: CurrentTime
        integer                                     :: IPC, IBC, iTimeSerie
        integer                                     :: i, j, k
        real                                        :: DT_Residual, DataValue, factor_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    

        !Shorten variavel name
        IPC = Me%InternalPropertyCount

cd1 :   if (ready_ .EQ. IDLE_ERR_) then
cd2 :       if (Me%Points) then

                !Actualize the Property Count
                IPC = IPC + 1
                if (IPC > Me%NumberOfProperties) IPC = 1
            
                !Stores the actual compute time
                call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT)

                if (present(factor))  then
                    factor_ = factor
                else
                    factor_ = 1.
                endif

                !Calculates the residual values
                do iTimeSerie = 1, Me%NumberOfTimeSeries
                
                    if (Me%TimeSerie(iTimeSerie)%IgnoreON) cycle

                    i = Icell(iTimeSerie)
                    j = Jcell(iTimeSerie)
                    k = Kcell(iTimeSerie)

                    DataValue = -1

                    if (i.ne.-1.and.J.ne.-1) then

                        if (present(Data2D))     DataValue = Data2D     (i, j)
                        if (present(Data3D))     DataValue = Data3D     (i, j, k)
                        if (present(Data2D_8))   DataValue = Data2D_8   (i, j)
                        if (present(Data3D_8))   DataValue = Data3D_8   (i, j, k)
                        if (present(Data2D_Int)) DataValue = Data2D_Int (i, j)
                        if (present(Data3D_Int)) DataValue = Data3D_Int (i, j, k)

                        DataValue = DataValue * factor_

                    endif
                    if(Me%ComputeResidual)then

                        DT_Residual = CurrentTime - Me%TimeSerie(iTimeSerie)%LastResidual


                        !Updates the Residual Values
                        Me%TimeSerie(iTimeSerie)%ResidualValues(IPC) =              &
                            (Me%TimeSerie(iTimeSerie)%ResidualValues(IPC)           &
                           * Me%TimeSerie(iTimeSerie)%ResidualTime                  &
                           + DataValue * DT_Residual)                               &
                           / (Me%TimeSerie(iTimeSerie)%ResidualTime + DT_Residual)
                    
                        !Updates the times
                        if (IPC == Me%NumberOfProperties) then
                            Me%TimeSerie(iTimeSerie)%ResidualTime =                 &
                                Me%TimeSerie(iTimeSerie)%ResidualTime + DT_Residual

                            Me%TimeSerie(iTimeSerie)%LastResidual = CurrentTime
                        endif

                    end if


                    !Stores the data in the buffer
                    if (CurrentTime .ge. Me%TimeSerie(iTimeSerie)%NextOutPut) then
                    
                        !Increments the internal buffer count
                        if (IPC == 1) then
                            Me%TimeSerie(iTimeSerie)%BufferCount =                  &
                                Me%TimeSerie(iTimeSerie)%BufferCount + 1 
                        endif

                        !Shorten Variable
                        IBC = Me%TimeSerie(iTimeSerie)%BufferCount

                        !Stores the current time
                        if (IPC == 1)                                         &
                            Me%TimeSerie(iTimeSerie)%TimeBuffer(IBC) = CurrentTime


                        if (Present(Data2D))                                        &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = Data2D(i, j) * factor_

                        if (Present(Data3D))                                        &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = Data3D(i, j, k) * factor_

                        if (Present(Data2D_8))                                      &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = Data2D_8(i, j) * factor_

                        if (Present(Data3D_8))                                      &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = Data3D_8(i, j, k) * factor_

                        if (Present(Data3D_Int))                                    &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = Data3D_Int(i, j, k) * factor_

                        if (Present(Data2D_Int))                                    &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = Data2D_Int(i, j) * factor_

                        !Sets next output time
                        if (IPC == Me%NumberOfProperties)                           &
                            !Me%TimeSerie(iTimeSerie)%NextOutPut =                   &
                            !    CurrentTime + Me%TimeSerie(iTimeSerie)%DT
                            Me%TimeSerie(iTimeSerie)%NextOutPut =                   &
                                Me%TimeSerie(iTimeSerie)%NextOutPut + Me%TimeSerie(iTimeSerie)%DT                        
                    endif

                    !Verifies if the buffer is full
                    if ((Me%TimeSerie(iTimeSerie)%BufferCount  ==                   &
                         Me%TimeSerie(iTimeSerie)%BufferSize)  .and.                &
                        (IPC == Me%NumberOfProperties))        then
                        call WriteBufferToFile(Me%TimeSerie(iTimeSerie), Me%NumberOfProperties)
                        Me%TimeSerie(iTimeSerie)%BufferCount = 0
                    endif

                enddo

                STAT_ = SUCCESS_
            else
                STAT_ = UNKNOWN_
            end if cd2
        else               

            STAT_ = ready_

        end if cd1


        Me%InternalPropertyCount = IPC

        if (present(STAT))                                                          &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine WriteTimeSerie3

    !--------------------------------------------------------------------------
    
    subroutine WriteTimeSerie4(TimeSerieID, Data2D, Data3D, Data2D_8, Data3D_8, &
                               Data2D_Int, Data3D_Int, I_Cell, J_Cell, K_Cell, factor, STAT)

        !Arguments-------------------------------------------------------------
        integer                                             :: TimeSerieID
        real,    dimension(:, :),    optional, pointer      :: Data2D
        real,    dimension(:, :, :), optional, pointer      :: Data3D
        real(8), dimension(:, :),    optional, pointer      :: Data2D_8
        real(8), dimension(:, :, :), optional, pointer      :: Data3D_8
        integer, dimension(:, :),    optional, pointer      :: Data2D_Int
        integer, dimension(:, :, :), optional, pointer      :: Data3D_Int
        integer,           intent(IN )                      :: I_Cell, J_Cell
        integer, optional, intent(IN )                      :: K_Cell
        real,    optional, intent(IN )                      :: factor
        integer, optional, intent(OUT)                      :: STAT


        !Local-----------------------------------------------------------------
        integer                                             :: ready_ , STAT_        
        type (T_Time)                                       :: CurrentTime
        integer                                             :: IPC, IBC, iTimeSerie
        integer                                             :: i, j, k
        integer                                             :: localKcell
        real                                                :: DT_Residual, DataValue, factor_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    

        !Shorten variavel name
        IPC = Me%InternalPropertyCount

cd1 :   if (ready_ .EQ. IDLE_ERR_) then
cd2 :       if (Me%Points) then

                !Actualize the Property Count
                IPC = IPC + 1
                if (IPC > Me%NumberOfProperties) IPC = 1
            
                if (present(K_cell)) then
                    localKcell = K_cell
                else 
                    localKcell = 0
                endif


                !Stores the actual compute time
                call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT)

                if (present(factor))  then
                    factor_ = factor
                else
                    factor_ = 1.
                endif

                !Calculates the residual values
                do iTimeSerie = 1, Me%NumberOfTimeSeries
                
                    if (Me%TimeSerie(iTimeSerie)%IgnoreON) cycle

                    i = I_cell
                    j = J_cell
                    k = localKcell

                    if (present(Data2D))     DataValue = Data2D     (i, j)
                    if (present(Data3D))     DataValue = Data3D     (i, j, k)
                    if (present(Data2D_8))   DataValue = Data2D_8   (i, j)
                    if (present(Data3D_8))   DataValue = Data3D_8   (i, j, k)
                    if (present(Data2D_Int)) DataValue = Data2D_Int (i, j)
                    if (present(Data3D_Int)) DataValue = Data3D_Int (i, j, k)

                    DataValue =  DataValue * factor_


                    if(Me%ComputeResidual)then

                        DT_Residual = CurrentTime - Me%TimeSerie(iTimeSerie)%LastResidual


                        !Updates the Residual Values
                        Me%TimeSerie(iTimeSerie)%ResidualValues(IPC) =              &
                            (Me%TimeSerie(iTimeSerie)%ResidualValues(IPC)           &
                           * Me%TimeSerie(iTimeSerie)%ResidualTime                  &
                           + DataValue * DT_Residual)                               &
                           / (Me%TimeSerie(iTimeSerie)%ResidualTime + DT_Residual)
                    
                        !Updates the times
                        if (IPC == Me%NumberOfProperties) then
                            Me%TimeSerie(iTimeSerie)%ResidualTime =                 &
                                Me%TimeSerie(iTimeSerie)%ResidualTime + DT_Residual

                            Me%TimeSerie(iTimeSerie)%LastResidual = CurrentTime
                        endif

                    end if


                    !Stores the data in the buffer
                    if (CurrentTime .ge. Me%TimeSerie(iTimeSerie)%NextOutPut) then
                    
                        !Increments the internal buffer count
                        if (IPC == 1) then
                            Me%TimeSerie(iTimeSerie)%BufferCount =                  &
                                Me%TimeSerie(iTimeSerie)%BufferCount + 1 
                        endif
                        
                        !Shorten Variable
                        IBC = Me%TimeSerie(iTimeSerie)%BufferCount

                        !Stores the current time
                        if (IPC == 1)                                         &
                            Me%TimeSerie(iTimeSerie)%TimeBuffer(IBC) = CurrentTime


                        if (Present(Data2D))                                        &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = Data2D(i, j) * factor_

                        if (Present(Data3D))                                        &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = Data3D(i, j, k) * factor_

                        if (Present(Data2D_8))                                      &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = Data2D_8(i, j) * factor_

                        if (Present(Data3D_8))                                      &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = Data3D_8(i, j, k) * factor_

                        if (Present(Data3D_Int))                                    &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = Data3D_Int(i, j, k) * factor_

                        if (Present(Data2D_Int))                                    &
                            Me%TimeSerie(iTimeSerie)%TimeSerieData(IPC, IBC) = Data2D_Int(i, j) * factor_

                        !Sets next output time
                        if (IPC == Me%NumberOfProperties)                           &
                            !Me%TimeSerie(iTimeSerie)%NextOutPut =                   &
                            !    CurrentTime + Me%TimeSerie(iTimeSerie)%DT
                            Me%TimeSerie(iTimeSerie)%NextOutPut =                   &
                                Me%TimeSerie(iTimeSerie)%NextOutPut + Me%TimeSerie(iTimeSerie)%DT                        
                    endif

                    !Verifies if the buffer is full
                    if ((Me%TimeSerie(iTimeSerie)%BufferCount  ==                   &
                         Me%TimeSerie(iTimeSerie)%BufferSize)  .and.                &
                        (IPC == Me%NumberOfProperties))        then
                        call WriteBufferToFile(Me%TimeSerie(iTimeSerie), Me%NumberOfProperties)
                        Me%TimeSerie(iTimeSerie)%BufferCount = 0
                    endif

                enddo

                STAT_ = SUCCESS_
            else
                STAT_ = UNKNOWN_
            end if cd2
        else               

            STAT_ = ready_

        end if cd1


        Me%InternalPropertyCount = IPC

        if (present(STAT))                                                          &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine WriteTimeSerie4


    !--------------------------------------------------------------------------
    
    
    subroutine WriteBufferToFile(TimeSerie, nProperties)

        !Arguments-------------------------------------------------------------
        type(T_TimeSerie)                           :: TimeSerie
        integer                                     :: nProperties

        !Local-----------------------------------------------------------------
        integer                                     :: iP, iB, PropNumber
        integer                                     :: unit
        real                                        :: Year, Month, Day, Hour, Minute, Second
        
        unit = TimeSerie%UnitNumber

        PropNumber = size (TimeSerie%TimeSerieData, DIM = 1)

        if (PropNumber > 1000)                                                      &   
            stop 'WriteBufferToFile - ModuleTimeSerie - ERR01'

        do iB = 1, TimeSerie%BufferCount

            !Writes date in the form of 1999 08 08 23 59 23
            call ExtractDate(TimeSerie%TimeBuffer(iB), Year = Year, Month = Month,   &
                             Day = Day, Hour = Hour, Minute = Minute, Second = Second)


            !Writes time in seconds since the beginning of the output
            !Writes all properties in the buffer in exp format
            if (Me%UseTabulatedData) then
                write(unit, fmt=1000) TimeSerie%TimeBuffer(iB) - TimeSerie%BeginOutPut,      &
                                    int(Year), int(Month), int(Day), int(Hour), int(Minute), &
                                    Second,                                                  &
                                    (TimeSerie%TimeSerieData(iP, iB), iP = 1, nProperties)
            else
                write(unit, fmt=2000) TimeSerie%TimeBuffer(iB) - TimeSerie%BeginOutPut,      &
                                    int(Year), int(Month), int(Day), int(Hour), int(Minute), &
                                    Second,                                                  &
                                    (TimeSerie%TimeSerieData(iP, iB), iP = 1, nProperties)
            endif

        enddo

    1000 format(1x, f13.2, 1x, i4, 2x, i2, 2x, i2, 2x, i2, 2x, i2, 2x, f7.4, 1x, 1000(24x, e20.12e3))
    2000 format(1x, f13.2, 1x, i4, 2x, i2, 2x, i2, 2x, i2, 2x, i2, 2x, f7.4, 1x, 1000(1x, e20.12e3))


    end subroutine WriteBufferToFile

    !---------------------------------------------------------------------------

    subroutine CorrectsCellsTimeSerie(TimeSerieID, it, i, j, k, STAT)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )              :: TimeSerieID
        integer,           intent(IN )              :: it
        integer, optional, intent(IN )              :: i, j, k
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            if (present(i)) Me%TimeSerie(it)%LocalizationI = i
            if (present(j)) Me%TimeSerie(it)%LocalizationJ = j
            if (present(k)) Me%TimeSerie(it)%LocalizationK = k
  
            STAT_ = SUCCESS_

        else              
         
            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine CorrectsCellsTimeSerie

    !--------------------------------------------------------------------------

   !---------------------------------------------------------------------------

    subroutine ReseatCurrentIndex(TimeSerieID, STAT)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )              :: TimeSerieID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, STAT_
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            Me%CurrentIndex = 2
  
            STAT_ = SUCCESS_

        else              
         
            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ReseatCurrentIndex

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine TryIgnoreTimeSerie(TimeSerieID, Number, IgnoreOK, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: TimeSerieID       
        integer, intent(IN)                         :: Number
        logical, intent(OUT)                        :: IgnoreOK
        integer, intent(OUT), optional              :: STAT


                  
        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            IgnoreOK = Me%TimeSerie(Number)%IgnoreON

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine TryIgnoreTimeSerie

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetNumberOfTimeSeries(TimeSerieID, NumberOfTimeSeries, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: TimeSerieID
        integer, intent(OUT)                        :: NumberOfTimeSeries
        integer, intent(OUT), optional              :: STAT
                  
        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 
 
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then
 
            NumberOfTimeSeries = Me%NumberOfTimeSeries

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetNumberOfTimeSeries

    !--------------------------------------------------------------------------

    subroutine GetTimeSerieLocation(TimeSerieID, Number, LocalizationI,                 &
                                    LocalizationJ, LocalizationK,                       &
                                    Latitude, Longitude,                                &
                                    CoordX, CoordY, CoordON,                            &
                                    DepthLevel, DepthON,                                &
                                    STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: TimeSerieID       
        integer, intent(IN)                         :: Number
        integer, intent(OUT), optional              :: LocalizationI
        integer, intent(OUT), optional              :: LocalizationJ
        integer, intent(OUT), optional              :: LocalizationK
        real,    intent(OUT), optional              :: Latitude
        real,    intent(OUT), optional              :: Longitude
        real,    intent(OUT), optional              :: CoordX
        real,    intent(OUT), optional              :: CoordY
        logical, intent(OUT), optional              :: CoordON
        real,    intent(OUT), optional              :: DepthLevel
        logical, intent(OUT), optional              :: DepthON
        integer, intent(OUT), optional              :: STAT


                  
        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then
 
            if (present(LocalizationI)) LocalizationI = Me%TimeSerie(Number)%LocalizationI
            if (present(LocalizationJ)) LocalizationJ = Me%TimeSerie(Number)%LocalizationJ
            if (present(LocalizationK)) LocalizationK = Me%TimeSerie(Number)%LocalizationK

            if (present(Longitude    )) Longitude = Me%TimeSerie(Number)%Longitude
            if (present(Latitude     )) Latitude  = Me%TimeSerie(Number)%Latitude
            if (present(CoordX       )) CoordX    = Me%TimeSerie(Number)%CoordX
            if (present(CoordY       )) CoordY    = Me%TimeSerie(Number)%CoordY
            if (present(CoordON      )) CoordON   = Me%TimeSerie(Number)%CoordON

            if (present(DepthLevel   )) DepthLevel= Me%TimeSerie(Number)%DepthLevel
            if (present(DepthON      )) DepthON   = Me%TimeSerie(Number)%DepthON


            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetTimeSerieLocation

    !--------------------------------------------------------------------------

    subroutine GetTimeSerieDT(TimeSerieID, iTimeSerie, DT, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: TimeSerieID
        integer, intent(IN)                         :: iTimeSerie
        integer, intent(OUT)                        :: DT
        integer, intent(OUT), optional              :: STAT
                  

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then
 
            DT = Me%TimeSerie(iTimeSerie)%DT

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetTimeSerieDT

    !--------------------------------------------------------------------------

    subroutine GetTimeSerieHeader(TimeSerieID, Header, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: TimeSerieID
        character(len=line_length), intent(OUT)     :: Header
        integer, intent(OUT), optional              :: STAT
                  

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then
 
            Header = Me%Header

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

    end subroutine GetTimeSerieHeader

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetTimeSerieName(TimeSerieID, iTimeSerie, Name, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: TimeSerieID
        integer,                    intent(IN)      :: iTimeSerie
        character(len=*),           intent(OUT)     :: Name
        integer, intent(OUT), optional              :: STAT
                  

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then
 
            Name = trim(Me%TimeSerie(iTimeSerie)%FromBlockFileName)
            
            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

    end subroutine GetTimeSerieName

    !--------------------------------------------------------------------------

    subroutine GetTimeSerieNextOutput(TimeSerieID, NextOutput, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: TimeSerieID
        type (T_Time),intent(OUT)                   :: NextOutput
        integer, intent(OUT), optional              :: STAT
                  
        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 
 
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then
 
            NextOutput = Me%TimeSerie(1)%NextOutput

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetTimeSerieNextOutput

    !--------------------------------------------------------------------------

    subroutine GetTimeSerieCycle(TimeSerieID, TimeCycle, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        logical                                     :: TimeCycle
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            TimeCycle = Me%TimeCycle

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_
        
    end subroutine GetTimeSerieCycle

    !--------------------------------------------------------------------------

    subroutine GetTimeSerieValue(TimeSerieID, CurrentTime, DataColumn, Time1, Value1,   &
                                 Time2, Value2, TimeCycle, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        type(T_Time),      intent(IN)               :: CurrentTime
        integer,           intent(IN)               :: DataColumn
        type(T_Time),      intent(OUT)              :: Time1
        real,              intent(OUT)              :: Value1
        type(T_Time),      intent(OUT)              :: Time2
        real,              intent(OUT)              :: Value2
        logical,           intent(OUT)              :: TimeCycle
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 
        real                                        :: Year, Month, Day
        real                                        :: Hour, Minute, Second
        integer                                     :: StoredColumn, JulDay


        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !Gets the date
            call ExtractDate(CurrentTime, Year, Month, Day, Hour, Minute, Second)

            StoredColumn = Me%FileColumns(DataColumn)

            if (Me%TimeCycle) then

                TimeCycle = .true.
                select case (trim(Me%CharTimeUnits))
                    case ('HOURS')
                        !Hours range from 0 - 23
                        Value1 = Me%DataMatrix(int(Hour+1), StoredColumn)
                    case ('DAYS')
                        !Days range from 1 - 366 - Julian Day
                        call JulianDay(CurrentTime, JulDay)
                        Value1 = Me%DataMatrix(JulDay,    StoredColumn)
                    case ('MONTHS')
                        !Month range from 1- 12
                        Value1 = Me%DataMatrix(int(Month),  StoredColumn)
                    case default
                        stop 'GetTimeSerieValue - ModuleTimeSerie - ERR01'
                end select

            else

                TimeCycle = .false.

                !Me%CurrentIndex = 2  <- Decrease execution speed with time. 
                do while (Me%InitialData + Me%DataMatrix(Me%CurrentIndex, 1) .lt. &
                          CurrentTime)
                    Me%CurrentIndex = Me%CurrentIndex + 1
                    
                    !Last Instant
                    if (Me%CurrentIndex >= Me%DataValues) then
                        Me%CurrentIndex = Me%DataValues
                        exit 
                    endif
                enddo
                
                Value1 = Me%DataMatrix(Me%CurrentIndex-1,   StoredColumn)
                Time1  = Me%InitialData + Me%DataMatrix(Me%CurrentIndex-1, 1)

                Value2 = Me%DataMatrix(Me%CurrentIndex, StoredColumn)
                Time2  = Me%InitialData + Me%DataMatrix(Me%CurrentIndex, 1)

            endif

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_
        
    end subroutine GetTimeSerieValue


    !--------------------------------------------------------------------------
    
   subroutine GetTimeSerieTimeFrameIndexes(TimeSerieID, StartTime, EndTime, StartIndex, EndIndex, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        type(T_Time),      intent(IN)               :: StartTime,  EndTime
        integer,           intent(OUT)              :: StartIndex, EndIndex
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 
        !Begin-----------------------------------------------------------------        


        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !Find Start Index
            do while (Me%InitialData + Me%DataMatrix(Me%StartIndex, 1) .lt. &
                      StartTime)
                Me%StartIndex = Me%StartIndex + 1
                
                !Last Instant
                if (Me%StartIndex >= Me%DataValues) then
                    Me%StartIndex = Me%DataValues
                    exit 
                endif
            enddo


            !Find End Index
            do while (Me%InitialData + Me%DataMatrix(Me%EndIndex, 1) .lt. &
                      EndTime)
                Me%EndIndex = Me%EndIndex + 1
                
                !Last Instant
                if (Me%EndIndex >= Me%DataValues) then
                    Me%EndIndex = Me%DataValues
                    exit 
                endif
            enddo                
            
            if (Me%InitialData + Me%DataMatrix(Me%EndIndex, 1) == EndTime) then
                 Me%EndIndex = Me%EndIndex - 1
            endif
            
            StartIndex = Me%StartIndex
            EndIndex   = Me%EndIndex            

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_
        
    end subroutine GetTimeSerieTimeFrameIndexes

    !--------------------------------------------------------------------------
    
    subroutine GetTimeSerieValueForIndex (TimeSerieID, Index, DataColumn, Value, STAT) 

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                 :: TimeSerieID
        integer, intent(IN)                 :: Index 
        integer, intent(IN)                 :: DataColumn
        real, intent(OUT)                   :: Value        
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_        

        !----------------------------------------------------------------------
        
        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if ((Index < 1) .or. (Index > Me%DataValues)) then
            
                STAT_ = OUT_OF_BOUNDS_ERR_
            
            else
                        
                Value = Me%DataMatrix(index, Me%FileColumns(DataColumn))
                STAT_ = SUCCESS_
            
            endif
            
        else 

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_
        
    end subroutine GetTimeSerieValueForIndex

    !--------------------------------------------------------------------------    
    
    subroutine GetTimeSerieTimeOfDataset (TimeSerieID, Index, Time, STAT)
    
        !Arguments------------------------------------------------------------- 
        integer, intent(IN)             :: TimeSerieID       
        integer, intent(IN)             :: Index
        type (T_Time), intent(OUT)      :: Time
        integer, optional, intent(OUT)  :: STAT        
        
        !Local-----------------------------------------------------------------
        integer :: STAT_
        integer :: ready_
        
        !----------------------------------------------------------------------
               
        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.   &
            (ready_ .EQ. READ_LOCK_ERR_)) then               
                                       
            if ((Index < 1) .OR. (Index > Me%DataValues)) then
            
                STAT_ = OUT_OF_BOUNDS_ERR_
            
            else
            
                Time = Me%InitialData + Me%DataMatrix(index, 1)
                STAT_ = SUCCESS_
            
            endif
           
        else 

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_
                    
        !----------------------------------------------------------------------
            
    end subroutine GetTimeSerieTimeOfDataset 
    
    !--------------------------------------------------------------------------
    
    subroutine GetTimeSerieTimeOfNextDataset (TimeSerieID, ActualTime, NextTime, STAT)        
    
        !Arguments------------------------------------------------------------- 
        integer, intent(IN)             :: TimeSerieID       
        Type (T_Time), intent(IN)       :: ActualTime
        type (T_Time), intent(OUT)      :: NextTime
        integer, optional, intent(OUT)  :: STAT        
        
        !Local-----------------------------------------------------------------
        integer :: index
        integer :: STAT_
        integer :: ready_         
        
        !----------------------------------------------------------------------
               
        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.   &
            (ready_ .EQ. READ_LOCK_ERR_)) then               
                                       
            NextTime = Me%InitialData + Me%DataMatrix(Me%CurrentIndex, 1)
            index = Me%CurrentIndex
                                               
do1:        do while (NextTime <= ActualTime)              
                if (index < Me%DataValues) then
                    index = index + 1
                    NextTime = Me%InitialData + Me%DataMatrix(index, 1)
                else
                    exit do1
                endif
            enddo do1
            
            STAT_ = SUCCESS_
            
        else 

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_
                    
        !----------------------------------------------------------------------
    
    end subroutine GetTimeSerieTimeOfNextDataset
        

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------


    subroutine GetTimeSerieTimeLimits(TimeSerieID, StartTime, EndTime, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        type(T_Time),      intent(OUT)              :: StartTime, EndTime
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            StartTime = Me%InitialData + Me%DataMatrix(1,             1)

            EndTime   = Me%InitialData + Me%DataMatrix(Me%DataValues, 1)

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_
        
    end subroutine GetTimeSerieTimeLimits

    !--------------------------------------------------------------------------


    real function GetTimeSerieIntegral(TimeSerieID, StartTime, EndTime, DataColumn, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        type(T_Time),      intent(IN)               :: StartTime, EndTime
        integer,           intent(IN)               :: DataColumn
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 
        type(T_Time)                                :: StartTimeSerie, EndTimeSerie, TSi, TEi, TSi1, TEi1
        integer                                     :: StoredColumn, StartIndex, EndIndex, i
        real                                        :: IntegAux, dt1, dt2, PS, PE, PSi, PSi1, PEi, PEi1


        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            StartTimeSerie = Me%InitialData + Me%DataMatrix(1,             1)

            EndTimeSerie   = Me%InitialData + Me%DataMatrix(Me%DataValues, 1)
            
            StoredColumn = Me%FileColumns(DataColumn)            
            
            if (Me%TimeCycle) then
                
                GetTimeSerieIntegral = TimeSerieCycleIntegral(StartTime, EndTime, StoredColumn) 
            
            else

                IntegAux = 0.

                if (StartTime > EndTime       ) stop 'GetTimeSerieIntegral - TimeSerie - ERR10'
                if (StartTime < StartTimeSerie) stop 'GetTimeSerieIntegral - TimeSerie - ERR20'
                if (EndTime   > EndTimeSerie  ) stop 'GetTimeSerieIntegral - TimeSerie - ERR30'


                do i = 1, Me%DataValues

                    if ((Me%InitialData + Me%DataMatrix(i+1, 1)) >= StartTime) then
                        StartIndex = i
                        TSi        = Me%InitialData + Me%DataMatrix(StartIndex  , 1)
                        TSi1       = Me%InitialData + Me%DataMatrix(StartIndex+1, 1)
                        exit
                    endif
                enddo

                do i = 1, Me%DataValues

                    if ((Me%InitialData + Me%DataMatrix(i+1, 1)) >= EndTime  ) then
                        EndIndex = i
                        TEi        = Me%InitialData + Me%DataMatrix(EndIndex,   1)
                        TEi1       = Me%InitialData + Me%DataMatrix(EndIndex+1, 1)
                        exit
                    endif
                enddo

                dt1 = StartTime - TSi
                dt2 = TSi1      - TSi
                PSi = Me%DataMatrix(StartIndex  , StoredColumn)
                PSi1= Me%DataMatrix(StartIndex+1, StoredColumn)

                PS  = (PSi * (dt2-dt1) + PSi1 * dt1) / dt2

                dt1 = EndTime   - TEi
                dt2 = TEi1      - TEi
                PEi = Me%DataMatrix(EndIndex  ,   StoredColumn)
                PEi1= Me%DataMatrix(EndIndex+1,   StoredColumn)

                PE  = (PEi * (dt2-dt1) + PEi1 * dt1) / dt2

i1:             if (StartIndex == EndIndex) then

                    IntegAux = (PS + PE) / 2. * (EndTime - StartTime)

                else i1

                    IntegAux =            (PS  + PSi1) / 2. * (TSi1    - StartTime)
                    IntegAux = IntegAux + (PEi + PE  ) / 2. * (EndTime - TEi      )

                    do i=StartIndex + 1, EndIndex-1
                        IntegAux = IntegAux + (Me%DataMatrix(i+1,StoredColumn) + Me%DataMatrix(i,StoredColumn)) / 2. * &
                                              (Me%DataMatrix(i+1,1           ) - Me%DataMatrix(i,1           ))
                    enddo

                endif i1
                
                GetTimeSerieIntegral = IntegAux
                
            endif                

            


            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_
        
    end function GetTimeSerieIntegral

    !--------------------------------------------------------------------------

    real function TimeSerieCycleIntegral(StartTime, EndTime, StoredColumn) 

        !Arguments-------------------------------------------------------------
        type(T_Time)                                :: StartTime, EndTime
        integer                                     :: StoredColumn

        !Local-----------------------------------------------------------------
        type (T_Time)                               :: CurrentTime
        real                                        :: DTAux, Year, Month, Day, Value1
        real                                        :: Hour, Minute, Second, IntegAux
        integer                                     :: JulDay
        real, parameter                             :: DTinteg = 3600.

        !Begin-----------------------------------------------------------------

        IntegAux    = 0. 
        CurrentTime = StartTime
        
        if (DTinteg > EndTime - StartTime) then
            DTaux = EndTime - StartTime
        else
            DTaux = DTinteg
        endif            

        do while (CurrentTime < EndTime) 

            CurrentTime = CurrentTime + DTaux 

            if (CurrentTime > EndTime) then
                DTaux       = EndTime - CurrentTime
                CurrentTime = EndTime 
            endif
            
            call ExtractDate(CurrentTime, Year, Month, Day, Hour, Minute, Second)            
            
            select case (trim(Me%CharTimeUnits))
                case ('HOURS')
                    !Hours range from 0 - 23
                    Value1 = Me%DataMatrix(int(Hour+1), StoredColumn)
                case ('DAYS')
                    !Days range from 1 - 366 - Julian Day
                    call JulianDay(CurrentTime, JulDay)
                    Value1 = Me%DataMatrix(JulDay,    StoredColumn)
                case ('MONTHS')
                    !Month range from 1- 12
                    Value1 = Me%DataMatrix(int(Month),  StoredColumn)
                case default
                    stop 'GetTimeSerieValue - ModuleTimeSerie - ERR01'
            end select

            IntegAux = IntegAux + Value1 * DTaux            
            
        enddo                        
                
        TimeSerieCycleIntegral = IntegAux
                
        
    end function TimeSerieCycleIntegral


    !--------------------------------------------------------------------------
    
    
    !--------------------------------------------------------------------------

    subroutine GetTimeSerieTimeUnits(TimeSerieID, TimeUnits, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        character(len=*)                            :: TimeUnits
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            TimeUnits = Me%CharTimeUnits
            
            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

    end subroutine GetTimeSerieTimeUnits

    !--------------------------------------------------------------------------

    logical function GetTimeSerieCheckDate (TimeSerieID, CurrentTime, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        type(T_Time)                                :: CurrentTime
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        type(T_Time)                                :: StartTime, EndTime
        integer                                     :: ready_         
        integer                                     :: STAT_ 

        STAT_ = UNKNOWN_
        
        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            StartTime  = Me%InitialData + Me%DataMatrix(1            , 1)
            EndTime    = Me%InitialData + Me%DataMatrix(Me%DataValues, 1)
            
            if (CurrentTime >= StartTime .and. CurrentTime <= EndTime) then
                GetTimeSerieCheckDate = .true.
            else
                GetTimeSerieCheckDate = .false.
            endif

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

    end function GetTimeSerieCheckDate

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetTimeSerieInitialData (TimeSerieID, InitialData, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        type(T_Time)                                :: InitialData
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 

        STAT_ = UNKNOWN_
        
        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            InitialData = Me%InitialData

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

    end subroutine GetTimeSerieInitialData

    !--------------------------------------------------------------------------

    subroutine GetTimeSerieDTForNextEvent     (TimeSerieID, Value, DataColumn,           &
                                               ActualTime, PredictedDT, DTForNextEvent,  &
                                               STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        real, intent(in)                            :: Value
        integer                                     :: DataColumn
        type (T_TIME)                               :: ActualTime
        real, intent(out)                           :: PredictedDT
        real, intent(out)                           :: DTForNextEvent
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 
        integer                                     :: index2, StoredColumn
        real                                        :: aux1, aux2
        type(T_Time)                                :: Time1, Time2, TimeWet

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            StoredColumn = Me%FileColumns(DataColumn)

            Time1  = Me%InitialData + Me%DataMatrix(Me%CurrentIndex-1, 1)
            Time2  = Me%InitialData + Me%DataMatrix(Me%CurrentIndex, 1)

            !MAXIMUM ALLOWED DT
            if (Value /= 0.0) then

                if (Time2 /= ActualTime) then
                    aux1 = Time2 - ActualTime
                else
                    aux1 = -FillValueReal
                endif

                !to ensure that DT does not go through two intervals of atmosphere
                aux2 = Me%DataMatrix(Me%CurrentIndex, 1) - Me%DataMatrix(Me%CurrentIndex-1, 1)

                PredictedDT     = min(aux1, aux2)
                DTForNextEvent  = 0.0

            else

                if (Me%CurrentIndex < Me%DataValues) then
                    !this else is active for dry periods (allowing higher DT until wetting events)
                    index2 = Me%CurrentIndex + 1
                    !find which is the next time in time serie where there is a wetting event 
                    do while (Me%DataMatrix(index2, StoredColumn) == 0.0)
                        index2 = index2 + 1
                        if (index2 > Me%DataValues) then
                            !No more values /= 0 until the end of the run
                            !Leave Predicted DT unchanged and return
                            DTForNextEvent = -null_real
                            return
                        end if
                    enddo

                    !instant where wetting event starts
                    TimeWet  = Me%InitialData + Me%DataMatrix(index2 - 1, 1)
      
                    if (TimeWet /= ActualTime) then
                        !To ensure that the model stops at the begining of the wetting event
                        aux1            = TimeWet - ActualTime
                        DTForNextEvent  = TimeWet - ActualTime
                    else
                        !To ensure that global DT is not higher then the Wetting event period
                        aux1            = Me%DataMatrix(index2, 1) - Me%DataMatrix(index2 - 1, 1)
                        DTForNextEvent  = 0.0
                    endif

                    PredictedDT = aux1


                else
                
                    DTForNextEvent = -null_real                    
                    
                endif

            endif


            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

    end subroutine GetTimeSerieDTForNextEvent

    !--------------------------------------------------------------------------

    subroutine GetTimeSerieDataMatrix (TimeSerieID, DataMatrix, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        real, dimension(:,:),      pointer          :: DataMatrix
        integer, optional                           :: STAT    

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            DataMatrix => Me%DataMatrix

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

    end subroutine GetTimeSerieDataMatrix

    !--------------------------------------------------------------------------

    subroutine GetTimeSerieDataColumns (TimeSerieID, DataColumns, STAT) 
        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        integer                                     :: DataColumns
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            DataColumns = Me%DataColumns

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

    end subroutine GetTimeSerieDataColumns

    !--------------------------------------------------------------------------

    subroutine GetTimeSerieDataValues (TimeSerieID, DataValues, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        integer                                     :: DataValues
        integer, optional                           :: STAT    

        !Local-----------------------------------------------------------------
        integer                                     :: ready_         
        integer                                     :: STAT_ 

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            DataValues = Me%DataValues

            STAT_ = SUCCESS_
        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

    end subroutine GetTimeSerieDataValues



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillTimeSerie(TimeSerieID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                 :: TimeSerieID
        integer, optional, intent(OUT)          :: STAT

        !Local-----------------------------------------------------------------
        integer                                 :: ready_ , STAT_             
        integer                                 :: STAT_CALL
        integer                                 :: iTimeSerie
        integer                                 :: nUsers

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(TimeSerieID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mTIMESERIE_,  Me%InstanceID)
  
            if (nUsers == 0) then


                !Disposes buffer
                do iTimeSerie = 1, Me%NumberOfTimeSeries
                
                    if (Me%TimeSerie(iTimeSerie)%IgnoreON) cycle
                    
                    call WriteBufferToFile(Me%TimeSerie(iTimeSerie), Me%NumberOfProperties)
                    
                enddo

                !Closes files
                call CloseTimeSerieFiles

                if (Me%ObjTime /= 0) then
                    nUsers = DeassociateInstance(mTIME_,            Me%ObjTime)
                    if (nUsers == 0) stop 'KillTimeSerie - ModuleTimeSerie - ERR00'
                endif    

                !Deallocates Time Serie Buffer
                do iTimeSerie = 1, Me%NumberOfTimeSeries
                
                    if (Me%TimeSerie(iTimeSerie)%IgnoreON) cycle
                
                    deallocate(Me%TimeSerie(iTimeSerie)%TimeSerieData, STAT = STAT_CALL)
                    !Intel Compiler gives here an STAT_CALL = 493
!                    if (STAT_CALL .NE. SUCCESS_) stop 'KillTimeSerie - ModuleTimeSerie - ERR01'
                    nullify(Me%TimeSerie(iTimeSerie)%TimeSerieData)

                    deallocate(Me%TimeSerie(iTimeSerie)%TimeBuffer,                 &
                               STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'KillTimeSerie - ModuleTimeSerie - ERR02'
                    nullify(Me%TimeSerie(iTimeSerie)%TimeBuffer)

                    if(Me%ComputeResidual)then
                        deallocate(Me%TimeSerie(iTimeSerie)%ResidualValues,STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                &
                            stop 'KillTimeSerie - ModuleTimeSerie - ERR02a'
                        nullify(Me%TimeSerie(iTimeSerie)%ResidualValues)
                    end if

                enddo


                !Deallocates the data values read
                if (associated(Me%DataMatrix)) then
                    deallocate(Me%DataMatrix,                                       &
                               STAT = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'KillTimeSerie - ModuleTimeSerie - ERR03'
                    nullify(Me%DataMatrix)


                endif

                !Deallocates the column numbers read
                if (associated(Me%ColumnsRead)) then
                    deallocate(Me%ColumnsRead,                                      &
                               STAT = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'KillTimeSerie - ModuleTimeSerie - ERR04'
                    nullify(Me%ColumnsRead)
                endif


                !Deallocates the column numbers read
                if (associated(Me%FileColumns)) then
                    deallocate(Me%FileColumns,                                      &
                               STAT = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'KillTimeSerie - ModuleTimeSerie - ERR05'
                    nullify(Me%FileColumns)
                endif


                !Deallocates TimeSeries
cd3 :           if (associated(Me%TimeSerie)) then
                    deallocate(Me%TimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'KillTimeSerie - ModuleTimeSerie - ERR06'
                    nullify(Me%TimeSerie)
                end if cd3


                call DeallocateInstance 

                TimeSerieID = 0

                STAT_ = SUCCESS_

            end if

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine KillTimeSerie

    !--------------------------------------------------------------------------

    subroutine CloseTimeSerieFiles

        !Local-----------------------------------------------------------------
        integer                                     :: iTimeSerie, STAT_CALL, unit, iP
        integer                                     :: nProperties
        real                                        :: Year, Month, Day, Hour, Minute, Second

        do iTimeSerie = 1, Me%NumberOfTimeSeries
        
            if (Me%TimeSerie(iTimeSerie)%IgnoreON) cycle

            nProperties = Me%NumberOfProperties
            unit = Me%TimeSerie(iTimeSerie)%UnitNumber
            call WriteDataLine(unit, block_end)

            if(Me%ComputeResidual)then

                call ExtractDate(Me%TimeSerie(iTimeSerie)%LastResidual,                 &
                                 Year = Year, Month = Month,                            &
                                 Day = Day, Hour = Hour, Minute = Minute, Second = Second)


                call WriteDataLine(unit, " ")
                call WriteDataLine(unit, " ")
                call WriteDataLine(unit, begin_residual)

                write(unit, fmt=1000) Me%TimeSerie(iTimeSerie)%ResidualTime,                    &
                                      int(Year), int(Month), int(Day), int(Hour), int(Minute),  &
                                      Second,                                                   &
                                      (Me%TimeSerie(iTimeSerie)%ResidualValues(iP),             &
                                       iP = 1, nProperties)

                call WriteDataLine(unit, end_residual)

            end if


            call UnitsManager(unit, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                write(*,*)'Error closing Time Serie data file',                     &
                           trim(Me%TimeSerie(iTimeSerie)%FileName)
                write(*,*)'CloseTimeSerieFiles - ModuleTimeSeries - WRN01'
            endif

        enddo

   1000 format(1x, f12.2, 1x, i4, 2x, i2, 2x, i2, 2x, i2, 2x, i2, 2x, f7.4, 1x, 1000(1x, e20.12e3))

    end subroutine CloseTimeSerieFiles

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance 

        !Local-----------------------------------------------------------------
        type (T_TimeSerieInOutPut), pointer         :: AuxTimeSerie
        type (T_TimeSerieInOutPut), pointer         :: PreviousTimeSerie

        !Updates pointers
        if (Me%InstanceID == FirstTimeSerie%InstanceID) then
            FirstTimeSerie => FirstTimeSerie%Next
        else
            PreviousTimeSerie => FirstTimeSerie
            AuxTimeSerie      => FirstTimeSerie%Next
            do while (AuxTimeSerie%InstanceID /= Me%InstanceID)
                PreviousTimeSerie => AuxTimeSerie
                AuxTimeSerie      => AuxTimeSerie%Next
            enddo

            !Now update linked list
            PreviousTimeSerie%Next => AuxTimeSerie%Next

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

    subroutine Ready (TimeSerieID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (TimeSerieID > 0) then
            call LocateObjTimeSerie (TimeSerieID)
            ready_ = VerifyReadLock (mTIMESERIE_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjTimeSerie (TimeSerieID)

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeSerieID

        !Local-----------------------------------------------------------------

        Me => FirstTimeSerie
        do while (associated (Me))
            if (Me%InstanceID == TimeSerieID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleTimeSerie - LocateObjTimeSerie - ERR01'

    end subroutine LocateObjTimeSerie

    !--------------------------------------------------------------------------

end module ModuleTimeSerie

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
