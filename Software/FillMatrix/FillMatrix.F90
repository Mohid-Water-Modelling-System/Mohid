!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : FillMatrix
! PROGRAM       : FillMatrix
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : April 2004
! REVISION      : Luis Fernandes
! DESCRIPTION   : Program to create Grid Data file from XYZ or TimeSeries
!
!------------------------------------------------------------------------------


!DataFile
!   PROPERTY_NAME           : char              -           !Name of the property whose values are being interpolated
!   PROPERTY_UNITS          : char              -           !Units of property whose values are being interpolated
!   VARIABLE_IN_TIME        : 0/1              [0]          !States if fields to be created are variable in time
!                                                           !If they are variable in time a HDF5 file is created
!   START                   : YYYY MM DD HH MM SS           !If variable in time, states begin time
!   END                     : YYYY MM DD HH MM SS           !If variable in time, states end time
!   MAX_TIME_SPAN           : real             [-]          !Maximum span for a time serie value to be considered valid
!   OUTPUT_TIME             : sec. sec. sec.    -           !Output Time
!   TIME_INSTANT            : YYYY MM DD HH MM SS           !If NOT variable in time, states time instant in which to 
!                                                           !extract from time serie
!   OUTPUT_FILE             : char              -           !Path to output GridData file
!   GRID_DATA_FILE          : char              -           !Path to Grid Data file(Bathymetry/Topography)
!   XYZ_FILE                : char              -           !Path to XYZ file (Mohid GIS format)
!   FILL_OUTSIDE_POINTS     : 0/1              [0]          !Check to perform extrapolation
!   WRITE_TRIANGLES         : 0/1              [0]          !Writes a polygon file with the created triangles
!   TRIANGLES_FILE          : char              -           !Path to the triangles file to create
!   INTERPOLATION_METHOD    : int              [-]          !0 - No Interpolation / 1 - Triangulation / 2 - IWD
!                                                                !if no interpolation, then spatial distribution is done by
!                                                                !providing grid with Fill ID. Model searches for this ID in
!                                                                !the Fill ID list (stations) to fill with correspondent TS 
!                                                                !or value
!   FILL_ID_FILE            : char              -           !Grid with fill IDs to make correspondence to fill ID's list (stations)
!
!   <begin_station>
!   NAME                    : char              -           !Name of the station
!   X                       : real              -           !XX coordinate of the station (in case of interpolation)
!   Y                       : real              -           !YY coordinate of the station (in case of interpolation)
!   ID                      : integer           -           !ID to search in fill ID grid (in case of no interpolation)
!   VALUE_TYPE              : char        ['SINGLE_VALUE']  !Type of information: 'SINGLE_VALUE', 'TIMESERIE' 
!   VALUE                   : real              -           !Property value of the station
!   FILENAME                : char              -           !Station time series file 
!   DATA_COLUMN             : int               -           !Column where base information is stored
!   <end_station>



program FillMatrix

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime
    use ModuleTimeSerie
    use ModuleDrawing
    use ModuleTriangulation
    use ModuleHorizontalGrid
    use ModuleHorizontalMap
    use ModuleGridData
    use ModuleFunctions
    use ModuleHDF5

    implicit none



    !Parameters---------------------------------------------------------
    integer, parameter                          :: NoInterpolation_= 0
    integer, parameter                          :: Triangulation_  = 1
    integer, parameter                          :: InverseWeight_  = 2


    !Types--------------------------------------------------------------
    type     T_File
        character(len=StringLength)             :: FileName
    end type T_File

    type     T_ExternalVar
        real, pointer, dimension(:,:)           :: XX_IE, YY_IE 
        type(T_Size2D)                          :: WorkSize
        integer, dimension(:,:  ), pointer      :: WaterPoints2D
        real,    dimension(:,:,:), pointer      :: WaterPoints3D
        real,    dimension(:,:  ), pointer      :: GridData
    end type T_ExternalVar


    type     T_FromTimeSerie
        integer                                 :: ObjTimeSerie         = 0
        type(T_File)                            :: File
        integer                                 :: DataColumn           = null_int
    end type T_FromTimeSerie

    type     T_Station
        character(len=StringLength)             :: Name
        type(T_PointF)                          :: Location
        integer                                 :: ID
        real                                    :: Value
        character(len=StringLength)             :: ValueType            = null_str
        logical                                 :: TimeSerieHasData     = .false.
        type(T_FromTimeSerie)                   :: TimeSerie
        type(T_Station), pointer                :: Next
    end type T_Station

    type     T_Nodes
        real, dimension(:), pointer             :: X
        real, dimension(:), pointer             :: Y
        real, dimension(:), pointer             :: Z
    end type T_Nodes

    type     T_TimeOptions
        type(T_Time)                            :: BeginTime
        type(T_Time)                            :: EndTime
        type(T_Time)                            :: OneInstant
        type(T_Time), pointer, dimension(:)     :: OutTime
        integer                                 :: Next                 = 0
        real                                    :: MaxTimeSpan
    end type T_TimeOptions

    type     T_FillMatrix
        type(T_File  )                          :: Input
        type(T_File  )                          :: Output    
        type(T_File  )                          :: Grid
        type(T_File  )                          :: XYZ
        type(T_File  )                          :: Triangles
        integer                                 :: ObjTime              = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjTriangulation     = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjHorizontalMap     = 0
        integer                                 :: ObjGridData          = 0
        integer                                 :: NumberOfStations     = 0
        type(T_Nodes)                           :: Nodes
        type(T_TimeOptions)                     :: Time
        type(T_ExternalVar)                     :: ExtVar
        character(len=StringLength)             :: PropertyName
        character(len=StringLength)             :: PropertyUnits
        character(len=StringLength)             :: FillIDFile
        type(T_PointF),dimension(:,:), pointer  :: GridPoint
        real,          dimension(:,:), pointer  :: GridData
        real, dimension(:,:), pointer           :: FillIDGrid
        logical                                 :: OutputTriangles      = .false.
        logical                                 :: VariableInTime       = .false.
        logical                                 :: HasOneTimeInstant    = .false.
        logical                                 :: FillOutsidePoints    = .false.
        logical                                 :: SkipNullValues       = .false.
        logical                                 :: HasTypes             = .false.
        integer                                 :: InterpolationMethod
        real                                    :: MaxDistance
        real                                    :: IWDn      
        type(T_XYZPoints),             pointer  :: XYZPoints
        type(T_Station), pointer                :: FirstStation
    end type T_FillMatrix

    type(T_FillMatrix)                          :: Me

    !Parameters---------------------------------------------------------
    character(len=StringLength), parameter      :: DataFile             = 'FillMatrix.dat'

    character(len=9 ),           parameter      :: FromTimeSerie        = 'TIMESERIE' 
    character(len=12),           parameter      :: SingleValue          = 'SINGLE_VALUE'
    
    !Begin--------------------------------------------------------------

    call StartUpMohid("FillMatrix")

    call ReadOptions

    if(Me%VariableInTime)then
        
        if (Me%InterpolationMethod /= NoInterpolation_) then

            call ConvertStationsToHDF5

        else

            call ConvertFromGridIDToHDF5
        
        endif

    else
        
        if (Me%InterpolationMethod /= NoInterpolation_) then

            call ConvertStationsToGridData

        else

            call ConvertFromGridIDToGridData
        
        endif

    end if

    !------------------------------------------------------------------

    contains

    !Subroutines-------------------------------------------------------
    
    subroutine ReadOptions

        !Local---------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin---------------------------------------------------------

        call ConstructEnterData(Me%ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - FillMatrix - ERR10'
    
        call ReadGlobalOptions

        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - FillMatrix - ERR20'

    
    end subroutine ReadOptions
    
    !------------------------------------------------------------------
    
    subroutine ConstructStations

        !Local---------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ClientNumber, iflag
        logical                                     :: BlockFound
        type(T_Station), pointer                    :: NewStation

        !Begin---------------------------------------------------------

        write(*,*)'Reading stations file...'
        write(*,*)


        nullify(Me%FirstStation, NewStation)

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,              &
                                        '<begin_station>', '<end_station>',         &
                                        BlockFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructStations - FillMatrix - ERR10'

if1 :       if(STAT_CALL .EQ. SUCCESS_) then 
   
if2 :           if (BlockFound) then

                    call AddStation(Me%FirstStation, NewStation)


                    call GetData(NewStation%Name,                                   &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'NAME',                             &
                                 ClientModule = 'FillMatrix',                       &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructStations - FillMatrix - ERR20'

                    write(*,*)'Station: '//trim(NewStation%Name)

                    if (Me%InterpolationMethod /= NoInterpolation_) then
                        call GetData(NewStation%Location%X,                         &
                                     Me%ObjEnterData, iflag,                        &
                                     SearchType   = FromBlock,                      &
                                     keyword      = 'X',                            &
                                     ClientModule = 'FillMatrix',                   &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructStations - FillMatrix - ERR30'

                        call GetData(NewStation%Location%Y,                         &
                                     Me%ObjEnterData, iflag,                        &
                                     SearchType   = FromBlock,                      &
                                     keyword      = 'Y',                            &
                                     ClientModule = 'FillMatrix',                   &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructStations - FillMatrix - ERR40'
                    
                    else

                        call GetData(NewStation%ID,                                 &
                                     Me%ObjEnterData, iflag,                        &
                                     SearchType   = FromBlock,                      &
                                     keyword      = 'ID',                           &
                                     ClientModule = 'FillMatrix',                   &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructStations - FillMatrix - ERR41'
                        if (iflag == 0)then
                            write(*,*)'Must specify station ID'
                            stop 'ConstructStations - FillMatrix - ERR42'
                        end if
                    endif

                    call GetData(NewStation%ValueType,                              &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'VALUE_TYPE',                       &
                                 Default      = SingleValue,                        &
                                 ClientModule = 'FillMatrix',                       &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructStations - FillMatrix - ERR50'

                    select case(trim(NewStation%ValueType))
                        
                        case(SingleValue)
            
                            call GetData(NewStation%Value,                          &
                                         Me%ObjEnterData, iflag,                    &
                                         SearchType   = FromBlock,                  &
                                         keyword      = 'VALUE',                    &
                                         ClientModule = 'FillMatrix',               &
                                         STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructStations - FillMatrix - ERR60'

                        case(FromTimeSerie)

                            call ReadStationTimeSerie(NewStation)

                    end select


                else

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop 'ConstructStations - FillMatrix - ERR70'
                        
                    exit do1    !No more blocks

                end if if2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                if(STAT_CALL .ne. SUCCESS_)stop 'ConstructStations - FillMatrix - ERR80'
                    
            end if if1
        end do do1


    end subroutine ConstructStations

    !------------------------------------------------------------------


    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag, ObjGD
        logical                                     :: FoundOutputTime
        real                                        :: DT

        !Begin---------------------------------------------------------

        write(*,*)'Reading global options...'
        write(*,*)

        call GetData(Me%PropertyName,                                               &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'PROPERTY_NAME',                                &
                     ClientModule = 'FillMatrix',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR10'

        call GetData(Me%PropertyUnits,                                              &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'PROPERTY_UNITS',                               &
                     ClientModule = 'FillMatrix',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR20'

        call GetData(Me%Output%FileName,                                            &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'OUTPUT_FILE',                                  &
                     ClientModule = 'FillMatrix',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR30'
        if (iflag == 0)then
            write(*,*)'Must specify type of file to convert'
            stop 'ReadGlobalOptions - FillMatrix - ERR40'
        end if
        

        call GetData(Me%VariableInTime,                                             &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'VARIABLE_IN_TIME',                             &
                     Default      = .false.,                                        &
                     ClientModule = 'FillMatrix',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR50'


        if(Me%VariableInTime)then


            call GetData(Me%Time%BeginTime,                                         &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'START',                                    &
                         ClientModule = 'FillMatrix',                               &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR80'

            call GetData(Me%Time%EndTime,                                           &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'END',                                      &
                         ClientModule = 'FillMatrix',                               &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR90'

            call GetData(Me%Time%MaxTimeSpan,                                       &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'MAX_TIME_SPAN',                            &
                         ClientModule = 'FillMatrix',                               &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR91'
            if (iflag == 0) then
                write(*,*)'Time Span not given'
                write(*,*)'Use Keyword MAX_TIME_SPAN'
                stop 'ReadGlobalOptions - FillMatrix - ERR91'
            endif
            
            call GetData(Me%SkipNullValues,                                         &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'SKIP_NULLVALUES',                          &
                         Default      = .false.,                                    &
                         ClientModule = 'FillMatrix',                               &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR91'
            

            nullify(Me%Time%OutTime)

            call GetOutPutTime(Me%ObjEnterData,                                     &
                               CurrentTime = Me%Time%BeginTime,                     &
                               EndTime     = Me%Time%EndTime,                       &
                               keyword     = 'OUTPUT_TIME',                         &
                               SearchType  = FromFile,                              &
                               OutPutsTime = Me%Time%OutTime,                       &
                               OutPutsOn   = FoundOutputTime,                       &
                               STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR60'
            if(.not. FoundOutputTime)  stop 'ReadGlobalOptions - FillMatrix - ERR70'


            DT = Me%Time%OutTime(2)-Me%Time%OutTime(1)

            call StartComputeTime(Me%ObjTime,                                       &
                                  Me%Time%BeginTime,                                &
                                  Me%Time%BeginTime,                                &
                                  Me%Time%EndTime,                                  &
                                  DT,                                               &
                                  .false., STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR92'





        else
        
            call GetData(Me%Time%OneInstant,                                        &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'TIME_INSTANT',                             &
                         ClientModule = 'FillMatrix',                               &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR100'

            if(iflag .ne. 0) Me%HasOneTimeInstant = .true.

        end if


        call GetData(Me%Grid%FileName,                                              &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'GRID_DATA_FILE',                               &
                     ClientModule = 'FillMatrix',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR110'
        if (iflag == 0)then
            write(*,*)'Must specify type of file to convert'
            stop 'ReadGlobalOptions - FillMatrix - ERR120'
        end if

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, trim(Me%Grid%FileName), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR130'

        call ConstructGridData (GridDataID       = Me%ObjGridData,                  &
                                HorizontalGridID = Me%ObjHorizontalGrid,            &
                                FileName         = Me%Grid%FileName,                &
                                STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR140'

        call ConstructHorizontalMap (HorizontalMapID  = Me%ObjHorizontalMap,        &
                                     GridDataID       = Me%ObjGridData,             &
                                     HorizontalGridID = Me%ObjHorizontalGrid,       &
                                     ActualTime       = Me%Time%OneInstant,         & 
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR160'

        

        call GetData (Me%InterpolationMethod, Me%ObjEnterData, iflag,               &
                      SearchType   = FromFile_,                                     &
                      keyword      ='INTERPOLATION_METHOD',                         &
                      ClientModule ='DigitalTerrainCreator',                        &
                      STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR172'
        if (iflag == 0) then
            write(*,*)'Interpolation Method not given'
            write(*,*)'Use Keyword INTERPOLATION_METHOD [0 - No interpolation / '   
            write(*,*)' 1 - Triangulation / 2 - IWD]'
            stop 'ReadGlobalOptions - FillMatrix - ERR173'
        endif


        select case (Me%InterpolationMethod)

        case (Triangulation_)

            call GetData(Me%FillOutsidePoints,                                          &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'FILL_OUTSIDE_POINTS',                          &
                         ClientModule = 'FillMatrix',                                   &
                         default      = .false.,                                        &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR170'


            call GetData (Me%OutputTriangles, Me%ObjEnterData, iflag,                   &
                          SearchType   = FromFile_,                                     &
                          keyword      ='WRITE_TRIANGLES',                              &
                          default      = .false.,                                       &
                          ClientModule ='DigitalTerrainCreator',                        &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR171'

            call GetData (Me%Triangles%FileName, Me%ObjEnterData, iflag,                &
                          SearchType   = FromFile_,                                     &
                          keyword      ='TRIANGLES_FILE',                               &
                          ClientModule ='DigitalTerrainCreator',                        &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR172'

        case (InverseWeight_)

            call GetData(Me%MaxDistance,                                                &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'MAX_DISTANCE',                                 &
                         ClientModule = 'FillMatrix',                                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR173'
            if (iflag == 0) then
                write(*,*)'Max Distance not given'
                write(*,*)'Use Keyword MAX_DISTANCE'
                stop 'ReadGlobalOptions - FillMatrix - ERR174'
            endif
            
            call GetData(Me%IWDn,                                                       &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'IWD_N',                                        &
                         default      = 2.0,                                            &
                         ClientModule = 'FillMatrix',                                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR175'
        
        case (NoInterpolation_)
            
            ! if no interpolation, then spatial distribution is done by providing grid with
            ! Fill ID. Model searches for this ID in the Fill ID list (stations) to fill with 
            ! correspondent TS or value
            call GetData(Me%FillIDFile,                                                 &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'FILL_ID_FILE',                                 &
                         ClientModule = 'FillMatrix',                                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR176'
            if (iflag == 0)then
                write(*,*)'Must specify grid with fill ID'
                stop 'ReadGlobalOptions - FillMatrix - ERR42'
            end if

            !Gets Fill IDs 
            ObjGD = 0
            call ConstructGridData  (GridDataID       = ObjGD,                        &
                                     HorizontalGridID = Me%ObjHorizontalGrid,         &
                                     FileName         = Me%FillIDFile,                &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR177'
        
            call GetGridData (ObjGD, Me%FillIDGrid, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR162'

        case default
            write(*,*)'Invalid Interpolation Method'
            write(*,*)'Use Keyword INTERPOLATION_METHOD [0 - No Interpolation / '
            write(*,*)'1 - Triangulation / 2 - IWD]'
            stop 'ReadGlobalOptions - FillMatrix - ERR174'
        end select

        if (Me%InterpolationMethod /= NoInterpolation_) then
            call GetData(Me%XYZ%FileName,                                               &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'XYZ_FILE',                                     &
                         ClientModule = 'FillMatrix',                                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - FillMatrix - ERR180'
        

            if (iflag .ne. 0) then
                call ConstructFromXYZFile
            else
                call ConstructStations
            end if
        else
            call ConstructStations
        endif

        write(*,*)'Number of stations...', Me%NumberOfStations
        write(*,*)
        


    end subroutine ReadGlobalOptions


    !--------------------------------------------------------------------------
    
    
    subroutine ReadStationTimeSerie(NewStation)   
    
        !Arguments-------------------------------------------------------------
        type (T_Station), pointer                   :: NewStation

        !Local-----------------------------------------------------------------
        integer                                     :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------


        call GetData(NewStation%TimeSerie%File%FileName,                &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'FILENAME',                         &
                     ClientModule = 'FillMatrix',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadStationTimeSerie - FillMatrix - ERR01'


        call GetData(NewStation%TimeSerie%DataColumn,                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'DATA_COLUMN',                      &
                     ClientModule = 'FillMatrix',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadStationTimeSerie - FillMatrix - ERR02'


        write(*,*)'Reading time serie for '//trim(NewStation%Name)

        if(Me%VariableInTime)then

            call StartTimeSerieInput(NewStation%TimeSerie%ObjTimeSerie,     &
                                     NewStation%TimeSerie%File%FileName,    &
                                     Me%ObjTime,                            &
                                     CheckDates = .false.,                  &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadStationTimeSerie - FillMatrix - ERR03'

        else
            call StartTimeSerieInput(NewStation%TimeSerie%ObjTimeSerie,     &
                                     NewStation%TimeSerie%File%FileName,    &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadStationTimeSerie - FillMatrix - ERR04'

        end if

        !Starts Time Serie

        write(*,*)'Done.'
        write(*,*)

    end subroutine ReadStationTimeSerie
    
    !------------------------------------------------------------------


    subroutine ConvertStationsToHDF5


        !Local-----------------------------------------------------------------
        type (T_Station), pointer                   :: Station
        integer                                     :: i, j, iStation, STAT_CALL
        integer                                     :: Instant
        type(T_Time)                                :: CurrentTime, PrevTime
        real, dimension(:,:), pointer               :: Field, PrevField
        character(5)                                :: AuxChar
        real                                        :: Nominator, Denominator, dist
        logical                                     :: WriteThisInstant, WritePrevInstant, PrevInstantWritten
        
        
        !Begin-----------------------------------------------------------------

        call ReadLockExternalVar

        call Open_HDF5_OutPut_File

        call DefineGridPoints

        allocate(Field(Me%ExtVar%WorkSize%ILB-1:Me%ExtVar%WorkSize%IUB+1, &
                       Me%ExtVar%WorkSize%JLB-1:Me%ExtVar%WorkSize%JUB+1))

        if (Me%SkipNullValues) then
            allocate(PrevField(Me%ExtVar%WorkSize%ILB-1:Me%ExtVar%WorkSize%IUB+1, &
                               Me%ExtVar%WorkSize%JLB-1:Me%ExtVar%WorkSize%JUB+1))
        endif


        do Instant = 1, size(Me%Time%OutTime)

            Me%Time%Next = Instant

            CurrentTime = Me%Time%OutTime(Me%Time%Next)
        
            !Updates station values and see if there are valid
            Me%NumberOfStations = 0
            Station             => Me%FirstStation
            do while(associated(Station))

                call UpDateStationValue(Station, CurrentTime)
               
                if (Station%TimeSerieHasData) then
                    Me%NumberOfStations = Me%NumberOfStations + 1
                endif

                Station => Station%Next
            end do

            !Select the interpolation method
            select case (Me%InterpolationMethod)

            case (Triangulation_)

                if (Me%NumberOfStations < 3) then
                    write (*,*)'Insufficient data avaliable'
                    write (*,*)'Increase MAX_TIME_SPAN or get more data ;-)'
                    stop 'ConvertStationsToHDF5 - FillMatrix - ERR01'
                endif
                
                !Allocates Variables for Triangulation
                allocate(Me%Nodes%X(1:Me%NumberOfStations))
                allocate(Me%Nodes%Y(1:Me%NumberOfStations))
                allocate(Me%Nodes%Z(1:Me%NumberOfStations))
            
                !Fills Station Values
                iStation            = 0
                Station             => Me%FirstStation
                do while(associated(Station))
                    if (Station%TimeSerieHasData) then                
                        iStation = iStation + 1
                        Me%Nodes%X(iStation) = Station%Location%X
                        Me%Nodes%Y(iStation) = Station%Location%Y
                        Me%Nodes%Z(iStation) = Station%Value
                    endif
                    Station => Station%Next
                end do

                !Constructs Triangulation
                call ConstructTriangulation (Me%ObjTriangulation,    &
                                             Me%NumberOfStations,    &
                                             Me%Nodes%X,             &
                                             Me%Nodes%Y,             &
                                             Me%OutputTriangles,     &
                                             STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConvertStationsToHDF5 - FillMatrix - ERR02'

                !Sets Height Values
                call SetHeightValues(Me%ObjTriangulation, Me%Nodes%Z, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConvertStationsToHDF5 - FillMatrix - ERR03'

                !Interpolates Values
                do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB
                do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB

                    if(Me%ExtVar%WaterPoints2D(i, j) == WaterPoint) then

                        Field(i, j) = InterPolation(Me%ObjTriangulation,          &
                                                    Me%GridPoint(i,j)%X,          &
                                                    Me%GridPoint(i,j)%Y,          &
                                                    Me%FillOutsidePoints,         &
                                                    Default = null_real,          &
                                                    STAT    = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConvertStationsToHDF5 - FillMatrix - ERR04'

                    else

                        Field(i, j) = null_real

                    end if

                enddo
                enddo

                write(AuxChar, fmt='(i5)')Me%Time%Next
                if(Me%OutputTriangles) call WriteTriangles(trim(adjustl(Me%Triangles%FileName))//"_"//trim(adjustl(AuxChar)))


                !Kills Triangulation
                call KillTriangulation (Me%ObjTriangulation, STAT    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConvertStationsToHDF5 - FillMatrix - ERR05'

                !Deallocates Matrixes 
                deallocate(Me%Nodes%X)
                deallocate(Me%Nodes%Y)
                deallocate(Me%Nodes%Z)
                

            case (InverseWeight_)

                !Interpolates Values
                do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB
                do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB

                    if(Me%ExtVar%WaterPoints2D(i, j) == WaterPoint) then

                        Nominator   = 0.0
                        Denominator = 0.0
                        Station     => Me%FirstStation
DoStations:             do while(associated(Station))
                            if (Station%TimeSerieHasData) then
                                dist = sqrt((Me%GridPoint(i,j)%X - Station%Location%X)**2.0 + &
                                            (Me%GridPoint(i,j)%Y - Station%Location%Y)**2.0)

                                if (dist > 0.0) then    
                                    if (dist < Me%MaxDistance) then 
                                        Nominator   = Nominator   + Station%Value / (dist ** Me%IWDn)
                                        Denominator = Denominator + 1.0 / (dist ** Me%IWDn)
                                    endif
                                else
                                    Nominator   = Station%Value
                                    Denominator = 1.0
                                    exit DoStations
                                endif
                            endif
                            Station => Station%Next
                        enddo DoStations

                        if (Denominator .lt. AllMostZero) then
                            write (*,*)'Insufficient data avaliable'
                            write (*,*)'Increase MAX_TIME_SPAN, MAX_DISTANCE or get more data ;-)'
                            write (*,*)'Point [i, j]', i, j
                            stop 'ConvertStationsToHDF5 - FillMatrix - ERR01'
                        endif
                        Field(i, j) = Nominator / Denominator

                    else

                        Field(i, j) = null_real

                    end if

                enddo
                enddo


            end select

            !Verifies if Field is to be written
            if (Me%SkipNullValues) then

                !The logic behind this algorithm is:
                !The user does not want to output zero values (util for accumulated Rain)
                !For accumulated rain, one must ensure that the zero value imediatly before
                !a rain event is also written, so the values are given correctly in module
                !atmosfere.

                !Always writes first and last instant
                if (Instant == 1 .or. Instant == size(Me%Time%OutTime)) then
                    WriteThisInstant = .true.
                    WritePrevInstant = .false.
                else
                    
                    WriteThisInstant = .false.
                    !Checks for maximum value
doj:                do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB
                    do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB
                        if(Me%ExtVar%WaterPoints2D(i, j) == WaterPoint) then
                            if (Field(i, j) > AllmostZero) then
                                WriteThisInstant = .true.
                                exit doj
                            endif
                        endif
                    enddo
                    enddo doj
                    
                    if (WriteThisInstant .and. .not. PrevInstantWritten) then
                        WritePrevInstant = .true.
                    else
                        WritePrevInstant = .false.
                    endif
                    
                endif
            else
                !Don't skip Zero values -> write all values, never write previous values
                WriteThisInstant = .true.
                WritePrevInstant = .false.
            endif
            
            if (Me%SkipNullValues) then
                if (WritePrevInstant) then
                    call WriteFieldToHDF5(PrevField, PrevTime)
                endif
            endif
            
             !Writes Interpolated Field to HDF File
            if (WriteThisInstant) then
                call WriteFieldToHDF5(Field, CurrentTime)
                PrevInstantWritten  = .true.
            else
                PrevInstantWritten  = .false.
            endif
            
            if (Me%SkipNullValues) then
                PrevField           = Field
                PrevTime            = CurrentTime
            endif

        end do


        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertStationsToHDF5 - FillMatrix - ERR40'

        write(*,*)
        write(*,*)'Finished...'
        write(*,*)

    end subroutine ConvertStationsToHDF5
    
    
    !----------------------------------------------------------------------


    subroutine ConvertFromGridIDToHDF5


        !Local-----------------------------------------------------------------
        type (T_Station), pointer                   :: Station
        integer                                     :: i, j, STAT_CALL
        integer                                     :: Instant
        type(T_Time)                                :: CurrentTime, PrevTime
        real, dimension(:,:), pointer               :: Field, PrevField
!        character(5)                                :: AuxChar
!        real                                        :: Nominator, Denominator, dist
        logical                                     :: WriteThisInstant, WritePrevInstant, PrevInstantWritten
        
        
        !Begin-----------------------------------------------------------------

        call ReadLockExternalVar

        call Open_HDF5_OutPut_File

!        call DefineGridPoints

        allocate(Field(Me%ExtVar%WorkSize%ILB-1:Me%ExtVar%WorkSize%IUB+1, &
                       Me%ExtVar%WorkSize%JLB-1:Me%ExtVar%WorkSize%JUB+1))

        if (Me%SkipNullValues) then
            allocate(PrevField(Me%ExtVar%WorkSize%ILB-1:Me%ExtVar%WorkSize%IUB+1, &
                               Me%ExtVar%WorkSize%JLB-1:Me%ExtVar%WorkSize%JUB+1))
        endif


        do Instant = 1, size(Me%Time%OutTime)

            Me%Time%Next = Instant

            CurrentTime = Me%Time%OutTime(Me%Time%Next)
        
            !Updates station values and see if there are valid
            Station             => Me%FirstStation
            do while(associated(Station))

                call UpDateStationValue(Station, CurrentTime)
               
                Station => Station%Next
            end do

            !No interpolation needed
            
            do j = Me%ExtVar%WorkSize%JLB , Me%ExtVar%WorkSize%JUB
            do i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB
                
                if(Me%ExtVar%WaterPoints2D(i, j) == WaterPoint) then
                   
                    Station     => Me%FirstStation
DoStations:         do while(associated(Station))
                        if (Station%TimeSerieHasData) then
                            
                            if (Station%ID == NINT(Me%FillIDGrid(i,j))) then
                                Field(i,j) = Station%Value        
                            endif
                        
                        endif

                        Station => Station%Next
                    enddo DoStations
                    
                    if (Field(i,j) .lt. 0) then
                        write(*,*) 'Fatal error, Type ID was not found'
                        write(*,*) 'in cell:', i,j
                    endif                
                
                else
                    
                    Field(i, j) = null_real
               
                endif
            
            enddo
            enddo


            !Verifies if Field is to be written
            if (Me%SkipNullValues) then

                !The logic behind this algorithm is:
                !The user does not want to output zero values (util for accumulated Rain)
                !For accumulated rain, one must ensure that the zero value imediatly before
                !a rain event is also written, so the values are given correctly in module
                !atmosfere.

                !Always writes first and last instant
                if (Instant == 1 .or. Instant == size(Me%Time%OutTime)) then
                    WriteThisInstant = .true.
                    WritePrevInstant = .false.
                else
                    
                    WriteThisInstant = .false.
                    !Checks for maximum value
doj:                do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB
                    do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB
                        if(Me%ExtVar%WaterPoints2D(i, j) == WaterPoint) then
                            if (Field(i, j) > AllmostZero) then
                                WriteThisInstant = .true.
                                exit doj
                            endif
                        endif
                    enddo
                    enddo doj
                    
                    if (WriteThisInstant .and. .not. PrevInstantWritten) then
                        WritePrevInstant = .true.
                    else
                        WritePrevInstant = .false.
                    endif
                    
                endif
            else
                !Don't skip Zero values -> write all values, never write previous values
                WriteThisInstant = .true.
                WritePrevInstant = .false.
            endif
            
            if (Me%SkipNullValues) then
                if (WritePrevInstant) then
                    call WriteFieldToHDF5(PrevField, PrevTime)
                endif
            endif
            
             !Writes Interpolated Field to HDF File
            if (WriteThisInstant) then
                call WriteFieldToHDF5(Field, CurrentTime)
                PrevInstantWritten  = .true.
            else
                PrevInstantWritten  = .false.
            endif
            
            if (Me%SkipNullValues) then
                PrevField           = Field
                PrevTime            = CurrentTime
            endif

        end do


        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertStationsToHDF5 - FillMatrix - ERR40'

        write(*,*)
        write(*,*)'Finished...'
        write(*,*)

    end subroutine ConvertFromGridIDToHDF5
    
    
    !----------------------------------------------------------------------
    
    subroutine WriteFieldToHDF5(Field, CurrentTime)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer               :: Field
        type(T_Time)                                :: CurrentTime

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        real,    dimension(6), target               :: AuxTime
        real,    dimension(:), pointer              :: TimePtr
        integer, save                               :: OutNumber = 1

        !----------------------------------------------------------------------

        call ExtractDate   (CurrentTime,                                &
                            AuxTime(1), AuxTime(2), AuxTime(3),         &
                            AuxTime(4), AuxTime(5), AuxTime(6))

        write(*,fmt=10)AuxTime(1), AuxTime(2), AuxTime(3), AuxTime(4), AuxTime(5), AuxTime(6)
        
    10  format('Processing Field ID: ', f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)

        TimePtr => AuxTime

        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteFieldToHDF5 - FillMatrix - ERR10'

        call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                             "Time", "YYYY/MM/DD HH:MM:SS",             &
                             Array1D        = TimePtr,                  &
                             OutputNumber   = OutNumber,                &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteFieldToHDF5 - FillMatrix - ERR20'
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB,    &
                             Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteFieldToHDF5 - FillMatrix - ERR30'

        call HDF5WriteData(Me%ObjHDF5,                                  &
                           "/Results/"//trim(Me%PropertyName),          &
                           trim(Me%PropertyName),                       &
                           trim(Me%PropertyUnits),                      &
                           Array2D      = Field,                        &
                           OutputNumber = OutNumber,                    &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteFieldToHDF5 - FillMatrix - ERR40'

        OutNumber = OutNumber + 1

        !Writes everything to disk
        !call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'WriteFieldToHDF5 - FillMatrix - ERR50'
    
    end subroutine WriteFieldToHDF5


    !----------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%Output%FileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - FillMatrix - ERR10'

        call HDF5SetLimits  (Me%ObjHDF5, Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB,    &
                             Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR20'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",                       &
                              Array2D =  Me%ExtVar%GridData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR30'            

        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR40'

        call HDF5SetLimits  (Me%ObjHDF5, Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB,    &
                             Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR50'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints", "-",                      &
                              Array2D = Me%ExtVar%WaterPoints2D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR60'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR70'

        

    end subroutine Open_HDF5_OutPut_File

    !------------------------------------------------------------------


    subroutine ConvertStationsToGridData

        !Local-----------------------------------------------------------------
        type (T_Station), pointer                   :: Station
        integer                                     :: i,j, STAT_CALL
        real                                        :: Nominator, Denominator, dist
        
        !Begin-----------------------------------------------------------------

        
        call ReadLockExternalVar


        select case (Me%InterpolationMethod)

        case (Triangulation_)

            allocate(Me%Nodes%X(1:Me%NumberOfStations))
            allocate(Me%Nodes%Y(1:Me%NumberOfStations))
            allocate(Me%Nodes%Z(1:Me%NumberOfStations))

            Station => Me%FirstStation
        
            i = 0

            write(*,*)'Filling triangulation nodes...'
            write(*,*)

            do while(associated(Station))

                i = i + 1

                call UpDateStationValue(Station, Me%Time%OneInstant)
                Me%Nodes%X(i) = Station%Location%X
                Me%Nodes%Y(i) = Station%Location%Y
                Me%Nodes%Z(i) = Station%Value
            

                Station => Station%Next

            end do


            write(*,*)'Constructing triangles network...'
            write(*,*)

            !Constructs Triangulation
            call ConstructTriangulation (Me%ObjTriangulation,                           &
                                         Me%NumberOfStations,                           &
                                         Me%Nodes%X,                                    &
                                         Me%Nodes%Y,                                    &
                                         Me%OutputTriangles,                            &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConvertStationsToGridData - FillMatrix - ERR10'


            call SetHeightValues(Me%ObjTriangulation, Me%Nodes%Z, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConvertStationsToGridData - FillMatrix - ERR20'


            allocate(Me%GridData(Me%ExtVar%WorkSize%ILB-1:Me%ExtVar%WorkSize%IUB+1, &
                                 Me%ExtVar%WorkSize%JLB-1:Me%ExtVar%WorkSize%JUB+1))

            call DefineGridPoints

            write(*,*)'Interpolating to grid data...'
            write(*,*)

            do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB
            do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB

                if(Me%ExtVar%WaterPoints2D(i, j) == WaterPoint) then

                    Me%GridData(i, j) = InterPolation(Me%ObjTriangulation,              &
                                                      Me%GridPoint(i,j)%X,              &
                                                      Me%GridPoint(i,j)%Y,              &
                                                      Me%FillOutsidePoints,             &
                                                      Default = null_real,              &
                                                       STAT    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConvertStationsToGridData - FillMatrix - ERR30'

                end if

            enddo
            enddo

            if(Me%OutputTriangles) call WriteTriangles(Me%Triangles%FileName)

        case (InverseWeight_)

            allocate(Me%GridData(Me%ExtVar%WorkSize%ILB-1:Me%ExtVar%WorkSize%IUB+1, &
                                 Me%ExtVar%WorkSize%JLB-1:Me%ExtVar%WorkSize%JUB+1))

            call DefineGridPoints

            write(*,*)'Interpolating to grid data...'
            write(*,*)


            !Interpolates Values
            do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB
            do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB

                if(Me%ExtVar%WaterPoints2D(i, j) == WaterPoint) then

                    Nominator   = 0.0
                    Denominator = 0.0
                    Station     => Me%FirstStation
DoStations:         do while(associated(Station))
                        if (Station%TimeSerieHasData) then
                            dist = sqrt((Me%GridPoint(i,j)%X - Station%Location%X)**2.0 + &
                                        (Me%GridPoint(i,j)%Y - Station%Location%Y)**2.0)

                            if (dist > 0.0) then    
                                if (dist < Me%MaxDistance) then 
                                    Nominator   = Nominator   + Station%Value / (dist ** Me%IWDn)
                                    Denominator = Denominator + 1.0 / (dist ** Me%IWDn)
                                endif
                            else
                                Nominator   = Station%Value
                                Denominator = 1.0
                                exit DoStations
                            endif
                        endif
                        Station => Station%Next
                    enddo DoStations

                    if (Denominator .lt. AllMostZero) then
                        write (*,*)'Insufficient data avaliable'
                        write (*,*)'Increase MAX_TIME_SPAN, MAX_DISTANCE or get more data ;-)'
                        write (*,*)'Point [i, j]', i, j
                        stop 'ConvertStationsToHDF5 - FillMatrix - ERR01'
                    endif
                    Me%GridData(i, j) = Nominator / Denominator

                else

                    Me%GridData(i, j) = null_real

                end if

            enddo
            enddo

        end select


        write(*,*)'Writing grid data...'
        write(*,*)

        call WriteGridData(FileName         = Me%Output%FileName,                   &
                           COMENT1          = null_str,                             &
                           COMENT2          = null_str,                             &
                           HorizontalGridID = Me%ObjHorizontalGrid,                 &
                           FillValue        = -99.,                                 &
                           Overwrite        = .true.,                               &
                           GridData2D_Real  = Me%GridData,                          &
                           STAT             = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConvertStationsToGridData - FillMatrix - ERR40'




        write(*,*)'Finished...'
        write(*,*)


    end subroutine ConvertStationsToGridData
    
    
    !------------------------------------------------------------------


    subroutine ConvertFromGridIDToGridData

        !Local-----------------------------------------------------------------
        type (T_Station), pointer                   :: Station
        integer                                     :: i,j, STAT_CALL
!        real                                        :: Nominator, Denominator, dist
        
        !Begin-----------------------------------------------------------------

        
        call ReadLockExternalVar



        allocate(Me%GridData(Me%ExtVar%WorkSize%ILB-1:Me%ExtVar%WorkSize%IUB+1, &
                             Me%ExtVar%WorkSize%JLB-1:Me%ExtVar%WorkSize%JUB+1))

!            call DefineGridPoints

        write(*,*)'Computing values to grid data...'
        write(*,*)


        Station => Me%FirstStation
    
        do while(associated(Station))

            call UpDateStationValue(Station, Me%Time%OneInstant)
       
            Station => Station%Next

        end do

        !No interpolation
        do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB
        do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB

            if(Me%ExtVar%WaterPoints2D(i, j) == WaterPoint) then

                Station     => Me%FirstStation
DoStations:     do while(associated(Station))
                    if (Station%TimeSerieHasData) then
                        if (Station%ID == NINT(Me%FillIDGrid(i,j))) then
                            Me%GridData(i,j) = Station%Value        
                        endif
                    endif
                   Station => Station%Next
                enddo DoStations
                
                if (Me%GridData(i,j) .lt. 0) then
                    write(*,*) 'Fatal error, Type ID was not found'
                    write(*,*) 'in cell:', i,j
                endif                        

            else

                Me%GridData(i, j) = null_real

            end if

        enddo
        enddo


        write(*,*)'Writing grid data...'
        write(*,*)

        call WriteGridData(FileName         = Me%Output%FileName,                   &
                           COMENT1          = null_str,                             &
                           COMENT2          = null_str,                             &
                           HorizontalGridID = Me%ObjHorizontalGrid,                 &
                           FillValue        = -99.,                                 &
                           Overwrite        = .true.,                               &
                           GridData2D_Real  = Me%GridData,                          &
                           STAT             = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConvertStationsToGridData - FillMatrix - ERR40'




        write(*,*)'Finished...'
        write(*,*)


    end subroutine ConvertFromGridIDToGridData
    
    
    !------------------------------------------------------------------
    
    subroutine DefineGridPoints

        !Local-----------------------------------------------------------------
        integer                             :: i, j
        real                                :: XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW

        !Begin-----------------------------------------------------------------

        write(*,*)"Defining grid points..."


        allocate(Me%GridPoint(Me%ExtVar%WorkSize%ILB-1:Me%ExtVar%WorkSize%IUB+1, &
                              Me%ExtVar%WorkSize%JLB-1:Me%ExtVar%WorkSize%JUB+1))

        do i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB
        do j = Me%ExtVar%WorkSize%JLB , Me%ExtVar%WorkSize%JUB
                       
            XSW = Me%ExtVar%XX_IE(i, j)
            YSW = Me%ExtVar%YY_IE(i, j)
            XSE = Me%ExtVar%XX_IE(i, j + 1)
            YSE = Me%ExtVar%YY_IE(i, j + 1)
            XNE = Me%ExtVar%XX_IE(i + 1, j + 1)
            YNE = Me%ExtVar%YY_IE(i + 1, j + 1)
            XNW = Me%ExtVar%XX_IE(i + 1, j)
            YNW = Me%ExtVar%YY_IE(i + 1, j)

            Me%GridPoint(i,j)%X = (XSW+XNW+XNE+XSE) / 4.
            Me%GridPoint(i,j)%Y = (YSW+YNW+YNE+YSE) / 4.


        end do
        end do

    end subroutine DefineGridPoints
    
    !------------------------------------------------------------------

    subroutine ConstructFromXYZFile

        !Local-----------------------------------------------------------------
        integer                                     :: i
        type (T_Station), pointer                   :: NewStation

        !Begin-----------------------------------------------------------------

        write(*,*)'Reading XYZ file...'
        write(*,*)


        call New(Me%XYZPoints, trim(Me%XYZ%FileName), DepthReferential = 1.)

        do i = 1, Me%XYZPoints%Count

            call AddStation(Me%FirstStation, NewStation)

            NewStation%Location%X       = Me%XYZPoints%X(i)
            NewStation%Location%Y       = Me%XYZPoints%Y(i)
            NewStation%Value            = Me%XYZPoints%Z(i)
            NewStation%ValueType        = SingleValue
            NewStation%TimeSerieHasData = .true.

        end do
    
    end subroutine ConstructFromXYZFile
    
    !------------------------------------------------------------------
    
    subroutine UpDateStationValue(Station, CurrentTime)
        
        !Arguments-------------------------------------------------------------
        type(T_Station), pointer                    :: Station
        type(T_Time)                                :: CurrentTime

        !Local-----------------------------------------------------------------
        logical                                     :: TimeCycle
        type (T_Time)                               :: Time1, Time2, InitialDate
        real                                        :: Value1, Value2
        real                                        :: dt1, dt2
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        
        select case(trim(Station%ValueType))
            
            case(SingleValue)

                !Remains Constant
                Station%TimeSerieHasData = .true.

            case(FromTimeSerie)

                call GetTimeSerieInitialData(Station%TimeSerie%ObjTimeSerie, InitialDate, &
                                             STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR09'

                
                if (CurrentTime >= InitialDate) then

                    !Gets Value for current Time
                    call GetTimeSerieValue (Station%TimeSerie%ObjTimeSerie,                   &
                                            CurrentTime,                                      &
                                            Station%TimeSerie%DataColumn,                     &
                                            Time1, Value1, Time2, Value2, TimeCycle,          &
                                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR10'

                    if (TimeCycle) then
                        Station%Value            = Value1
                        Station%TimeSerieHasData = .true.
                    else

                        dt1 = CurrentTime - Time1
                        dt2 = Time2 - CurrentTime

                        if (dt1 <= Me%Time%MaxTimeSpan .and. dt2 <= Me%Time%MaxTimeSpan) then

                            Station%TimeSerieHasData = .true.

                            !Interpolates Value for current instant
                            call InterpolateValueInTime(CurrentTime, Time1, Value1, Time2, Value2, Station%Value)
                        
                        else

                            Station%TimeSerieHasData = .false.

                        endif
                    
                                                
                    endif

                else

                    Station%TimeSerieHasData = .false.

                endif

            case default

        end select

    end subroutine UpDateStationValue

    !------------------------------------------------------------------

    
    subroutine ReadLockExternalVar
        
        !Local-----------------------------------------------------------------
        integer                                     :: GEOG, UTM, MIL_PORT, SIMPLE_GEOG, NLRD
        integer                                     :: GRID_COORD, CoordType, STAT_CALL
                
        !Begin-----------------------------------------------------------------

        call GetHorizontalGridSize(Me%ObjHorizontalGrid, WorkSize = Me%ExtVar%WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - FillMatrix - ERR10'
        
        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%ExtVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - FillMatrix - ERR20'

        !Gets Coordinates in use
        call GetGridCoordType(Me%ObjHorizontalGrid, CoordType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - FillMatrix - ERR30'
        
        call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,             &
                               SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD,      &
                               NLRD = NLRD)

        if    (CoordType == SIMPLE_GEOG)then
            
            call GetGridLatitudeLongitude(Me%ObjHorizontalGrid,                 &
                                          GridLatitudeConn  = Me%ExtVar%YY_IE,  &
                                          GridLongitudeConn = Me%ExtVar%XX_IE,  &
                                          STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - FillMatrix - ERR40'

        elseif(CoordType == UTM        .or. CoordType == MIL_PORT .or. &
               CoordType == GRID_COORD .or. CoordType == NLRD)then

            call GetHorizontalGrid(Me%ObjHorizontalGrid,                        &
                                   XX_IE = Me%ExtVar%XX_IE,                     &
                                   YY_IE = Me%ExtVar%YY_IE,                     &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - FillMatrix - ERR50'

        else

            write(*,*)'GEOG coordinate type cannot be used in digital terrain generation'
            stop 'Read_Lock_External_Var - FillMatrix - ERR60'

        end if

        call GetGridData(Me%ObjGridData, Me%ExtVar%GridData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR70'




    end subroutine ReadLockExternalVar


    !------------------------------------------------------------------


    subroutine AddStation (FirstStation, ObjStation)

        !Arguments-------------------------------------------------------------
        type (T_Station), pointer                   :: FirstStation
        type (T_Station), pointer                   :: ObjStation

        !Local-----------------------------------------------------------------
        type (T_Station), pointer                   :: NewStation
        type (T_Station), pointer                   :: PreviousStation
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewStation)
        nullify  (NewStation%Next)

        Me%NumberOfStations = Me%NumberOfStations + 1

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstStation)) then
            FirstStation         => NewStation
            ObjStation           => NewStation
        else
            PreviousStation      => FirstStation
            ObjStation           => FirstStation%Next
            do while (associated(ObjStation))
                PreviousStation  => ObjStation
                ObjStation       => ObjStation%Next
            enddo
            ObjStation           => NewStation
            PreviousStation%Next => NewStation
        endif


    end subroutine AddStation
    
    !--------------------------------------------------------------------------

    subroutine WriteTriangles (TrianglesFileName)

        !Arguments-------------------------------------------------------------
        character(len=*)                :: TrianglesFileName
        
        !Local-----------------------------------------------------------------
        integer                         :: STAT_CALL, nNodes
        integer                         :: UnitNumber, iT, nTriangles
        real,    dimension(:), pointer  :: XT, YT, ZT
        integer, dimension(:), pointer  :: V1, V2, V3

        !Begin-----------------------------------------------------------------

        write(*,*)"Writing triangles file..."

        !Get the number of triangles
        call GetNumberOfTriangles   (Me%ObjTriangulation, nTriangles)

        !Allocates space for the Triangle vertices and gets them
        allocate(V1(nTriangles))
        allocate(V2(nTriangles))
        allocate(V3(nTriangles))

        call GetTriangleList (Me%ObjTriangulation, v1, v2, v3, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteTriangles - FillMatrix - ERR10'


        !Gets nodes effictive used and the reordered nodes 
        call GetNumberOfNodes (Me%ObjTriangulation, nNodes, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteTriangles - FillMatrix - ERR20'

        allocate(XT(nNodes))
        allocate(YT(nNodes))
        allocate(ZT(nNodes))

        call GetNodesList   (Me%ObjTriangulation, XT, YT, ZT, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteTriangles - FillMatrix - ERR30'


        call UnitsManager (UnitNumber, OPEN_FILE, STAT = STAT_CALL)
        open (unit=UnitNumber, status = 'unknown', file = TrianglesFileName)
        do iT = 1, nTriangles
            write(UnitNumber,*)'<beginpolygon>'
            write(UnitNumber,*)XT(V1(iT)), YT(V1(iT))
            write(UnitNumber,*)XT(V2(iT)), YT(V2(iT))
            write(UnitNumber,*)XT(V3(iT)), YT(V3(iT))
            write(UnitNumber,*)XT(V1(iT)), YT(V1(iT))
            write(UnitNumber,*)'<endpolygon>'
        enddo
        call UnitsManager (UnitNumber, CLOSE_FILE, STAT = STAT_CALL)

        deallocate(XT, YT, ZT, v1, v2, v3)


    end subroutine WriteTriangles
        

end program FillMatrix





