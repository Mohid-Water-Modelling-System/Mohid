!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Inverted Barometer
! MODULE        : InvertedBarometer
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : October 2005
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Module to create mean water level time serie for gauge 
!                 locations based on inverted barometer hypothesis
!
!------------------------------------------------------------------------------

!DataFile
!
!   <BeginHDF5File>                                                                 
!   NAME                    : char                  [-]         !Name of HDF5 file 
!   <EndHDF5File>                                               !(specify one block for each file)
!                                                                                    
!   START_TIME              : YYYY MM DD HH MM SS   [-]         !Start date of time series
!   END_TIME                : YYYY MM DD HH MM SS   [-]         !End date of time series
!
!   <BeginGaugeFile>                                            
!   NAME                    : char                  [-]         !Path to gauge file
!
!   TIMESERIE_NAME          : char                  [-]         !Name for output time serie file (without extension) 
!                                                               !(if more than one gauge in file is first part of name)
!                                                               !(specify one block per gauge file)                    
!
!   WRITE_MOHID_INPUT       : 0/1                   [0]         !Create a twin of the gauge file with new keywords to
!                                                               !allow MOHID Water to consider gauge reference level
!                                                               !evolution according with created time series
!
!   OUTPUTNAME              : char                  [-]         !Name of the improved twin gauge file to create
!   <EndGaugeFile>                                                                                                           
!
!   MAX_BUFFER_SIZE         : integer               [100000]    !Maximum size for time serie buffer of each gauge location
!
!   GRID_FILENAME           : char                  [-]         !Path to HDF5 grid file

Module ModuleInvertedBarometer

    use ModuleGlobalData
    use ModuleTime              
    use ModuleEnterData,            only : ConstructEnterData, KillEnterData, GetData,      &
                                           ExtractBlockFromBuffer, Block_Unlock,            &
                                           GetBlockSize, GetFullBufferLine, WriteDataLine
    use ModuleHDF5,                 only : GetHDF5FileAccess, ConstructHDF5,                &
                                           GetHDF5GroupNumberOfItems, GetHDF5GroupID,       &
                                           HDF5SetLimits, HDF5ReadData, KillHDF5
    use ModuleHorizontalGrid,       only : UnGetHorizontalGrid, ConstructHorizontalGrid,    &
                                           GetHorizontalGridSize, GetGridCoordType,         & 
                                           GetGridLatitudeLongitude
    use ModuleDrawing,              only : IsPointInsidePolygon, T_Polygon, T_PointF

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartInvertedBarometer
    private ::      ConstructInvertedBarometerLevel
    private ::          ReadKeywords
    private ::              ReadGlobalData
    private ::              ReadGaugesFileName
    private ::                  AddGaugeFile
    private ::                  ConstructGaugeFile
    private ::              ReadHDF5FileName
    private ::                  AddHDF5File
    private ::                  ConstructHDF5File
    private ::          OpenAndDateHDF5Files
    private ::              HDF5Evaluator
    private ::                  SamePeriodHDF5
    private ::              HDF5TimeInstant
    private ::              ObtainInstantsTimes
    private ::              AddTSHDF5File
    private ::                  CreateTSHDF5File
    private ::              AddTSInstantsTimes
    private ::          ConstructHDF5Grid
    private ::              ConstructGridPolygon
    private ::                  SetLimitsPolygon
    private ::          ConstructGaugeList
    private ::              AddGaugeLocation
    private ::              ConstructGaugeLocation
    private ::                  Calc_Decimal_geo_coord
    private ::              CheckInsideGridPolygon
    private ::          ConstructGaugeGridLocation
    private ::              ConstructCellPolygon
    private ::          OpenOutputFiles
    private ::              AllocateTimeSerieBuffer
    private ::              CreateGaugeBlocks   

    !Selector
    
    !Modifier
    public  :: ModifyInvertedBarometerLevel
    private ::      OpenAndReadHDF5File
    private ::          ReadPressureFields
    private ::              AddField
    private ::      ComputeInvertedBarometer
    private ::      WriteIBTimeSerie
    private ::          WriteBufferToFile
    private ::      KillFields
    private ::      KillGaugeTimeSeries
    private ::          CloseTimeSerieFile

    !Destructor
    public  :: KillInvertedBarometerLevel
    private ::      KillIndividualHDF5File                                                     

    !Management
    
    !Interfaces----------------------------------------------------------------

    !Types---------------------------------------------------------------------

    type       T_GaugeLocation
        character(len=PathLength)                   :: Name                 = null_str
        real, dimension(3)                          :: Longitude        = FillValueReal
        real, dimension(3)                          :: Latitude         = FillValueReal
        real                                        :: DecimalLatitude  = FillValueReal
        real                                        :: DecimalLongitude = FillValueReal
        real                                        :: Metric_X         = FillValueReal
        real                                        :: Metric_Y         = FillValueReal
        logical                                     :: DefinedMetric        = .false.
        real                                        :: ReferenceLevel   = FillValueReal
        integer                                     :: ObjTimeSerie         = 0
        integer                                     :: I
        integer                                     :: J
        integer                                     :: BufferCount          = null_int
        logical                                     :: CellFound = .false.
        real, dimension(:),             pointer     :: TimeSerieData
        type (T_Time), dimension(:),    pointer     :: TimeBuffer
        type(T_GaugeLocation),          pointer                 :: Next     => null()
    end type  T_GaugeLocation

    type       T_GaugeFile
        character(len=PathLength)                   :: Name                 = null_str
        character(len=PathLength)                   :: TimeSerieName        = null_str
        character(len=PathLength)                   :: OutputName           = null_str
        integer                                     :: ObjEnterData         = 0
        integer                                     :: ObjOutput            = 0
        integer                                     :: TotalGauges          = 0
        integer                                     :: FoundGauges          = 0
        logical                                     :: WriteMOHIDInput
        type (T_GaugeLocation ),        pointer     :: FirstGaugeLocation
        type(T_GaugeFile),              pointer                 :: Next     => null()
    end type  T_GaugeFile

    ! Definition of type T_TimeSerieTime
    type       T_TimeSerieTime
        type(T_Time)                                :: Time
        type(T_TimeSerieTime),          pointer                 :: Next     => null()
    end type  T_TimeSerieTime

    ! Definition of type T_Field
    type       T_Field
        character(len=StringLength)                 :: Name                 = null_str
        real, dimension(:,:  ),         pointer     :: Values2D
        type(T_Field),                  pointer                 :: Next     => null()
    end type  T_Field

    type       T_HDF5File
        integer                                     :: HDFID                = 0
        character(len=StringLength)                 :: Name                 = null_str
        type(T_Time)                                :: StartTime
        type(T_Time)                                :: EndTime
        type(T_Time)                                :: StartFieldTime
        type(T_Time)                                :: EndFieldTime
        integer                                     :: NumberOfInstants
        integer                                     :: StartInstant, EndInstant
        type(T_TimeSerieTime),          pointer     :: FirstInstantTime
        type(T_Field),                  pointer     :: FirstField
        type(T_Field),                  pointer     :: CurrentField
        type(T_HDF5File),               pointer                 :: Next     => null()
    end type  T_HDF5File
    
    private :: T_InvertedBarometer
    type       T_InvertedBarometer
        integer                                     :: ObjEnterData         = 0
        integer                                     :: ObjTime              = 0
        integer                                     :: ObjHorizontalGrid    = 0
        integer                                     :: TotalGauges          = 0
        integer                                     :: FoundGauges          = 0
        integer                                     :: MaxBufferSize        = null_int
        integer                                     :: BufferSize           = null_int
        integer                                     :: TotalOutPutsNumber   = 0
        integer                                     :: GridCoordType        = null_int
        logical                                     :: VariableDT
        character(PathLength)                       :: DataFile, GridFile   = null_str
        type (T_Size2D)                             :: Size2D, WorkSize2D
        type(T_Time)                                :: StartTSTime, EndTSTime
        type(T_Time)                                :: CurrentTime
        type(T_GaugeFile),              pointer     :: FirstGaugeFile
        type(T_HDF5File),               pointer     :: FirstHDF5File
        type(T_HDF5File),               pointer     :: FirstTSHDF5File
        type(T_HDF5File),               pointer     :: LastTSHDF5File
        type(T_TimeSerieTime),          pointer     :: FirstTimeSerieTime
        type(T_Polygon),                pointer     :: GridPolygon
        real, dimension(:,:),           pointer     :: GridLongitudeConn 
        real, dimension(:,:),           pointer     :: GridLatitudeConn
        type(T_InvertedBarometer),      pointer                 :: Next     => null()
    end type  T_InvertedBarometer

    !Global Module Variables
    type (T_InvertedBarometer),         pointer     :: Me

    integer                                         :: mInvertedBarometer_ = 0 !just to compile

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartInvertedBarometer(ObjInvertedBarometerID, DataFile)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjInvertedBarometerID 
        character(PathLength), intent(IN)               :: DataFile 

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------

        !------------------------------------------------------------------------

        !Assures nullification of the global variable
        allocate(Me)

        !Returns ID
        ObjInvertedBarometerID          = 1

        !Atribute the name of data file            
        Me%DataFile = DataFile

        call ConstructInvertedBarometerLevel

        !----------------------------------------------------------------------

    end subroutine StartInvertedBarometer
 
    !--------------------------------------------------------------------------
   
    subroutine ConstructInvertedBarometerLevel

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        ! Read keyword file and HDF5 file
        call ReadKeywords

        ! Open and obtain key features of the HDF5 files
        call OpenAndDateHDF5Files

        ! Construct grid for communication with gauge locations
        call ConstructHDF5Grid
        ! (it is assumed that all HDF5 files have the same grid!)        

        ! Construct gauge location list
        call ConstructGaugeList

        ! Construct gauge grid location
        call ConstructGaugeGridLocation

        ! Open and write header of Output Files
        call OpenOutputFiles

    end subroutine ConstructInvertedBarometerLevel

    !--------------------------------------------------------------------------
  
    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call ConstructEnterData (Me%ObjEnterData, Me%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'ReadKeywords - ModuleInvertedBarometer - ERR10'

        call ReadGlobalData

        call ReadGaugesFileName

        call ReadHDF5FileName

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadKeywords - ModuleInvertedBarometer - ERR20'
        end if

    end subroutine ReadKeywords

    !--------------------------------------------------------------------------

    subroutine ReadGlobalData

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------

        ! Obtain the start and end times for the Time Serie
        ! Start Time
        call GetData(Me%StartTSTime, Me%ObjEnterData, iflag,                    &
                     keyword      = 'START_TIME',                               &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'InvertedBarometer',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                              &
        stop 'ReadGlobalData - ModuleInvertedBarometer - ERR10'   

        ! End Time 
        call GetData(Me%EndTSTime,   Me%ObjEnterData, iflag,                    &
                     keyword      = 'END_TIME',                                 &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'InvertedBarometer',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                              &
        stop 'ReadGlobalData - ModuleInvertedBarometer - ERR20'   

        ! Verifies Time Variables
        if (Me%EndTSTime .lt. Me%StartTSTime) then
            write (*,*) 'Time Serie End Time is BEFORE Time Serie Start Time'
            write (*,*) 'Module :','InvertedBarometer'
            stop 'ReadGlobalData - ModuleInvertedBarometer - ERR30'
        endif

        if (Me%EndTSTime .eq. Me%StartTSTime) then
            write (*,*) 'Time Serie End Time is EQUAL Time Serie Start Time'
            write (*,*) 'Module :','InvertedBarometer'
            stop 'ReadGlobalData - ModuleInvertedBarometer - ERR40'
        endif

        !Obtain grid file name
        call GetData(Me%GridFile, Me%ObjEnterData, iflag,                       &
                     keyword      = 'GRID_FILENAME',                            &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'InvertedBarometer',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                              &
        stop 'ReadGlobalData - ModuleInvertedBarometer - ERR50'   

        !The maximum BufferSize is set here to 0.1Mb (for each property)
        !This lets perform 25000 outputs to the buffer (considering each output of 4 bytes)
        call GetData(Me%MaxBufferSize,                                          &
                     Me%ObjEnterData,                                           &
                     iflag,                                                     &
                     SearchType   = FromFile,                                   &
                     keyword      ='MAX_BUFFER_SIZE',                           &
                     Default      = 100000,                                     &
                     ClientModule ='InvertedBarometer',                         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'ReadGlobalData - ModuleInvertedBarometer - ERR60'   

    end subroutine ReadGlobalData

    !--------------------------------------------------------------------------

    subroutine ReadGaugesFileName

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, ClientNumber
        type (T_GaugeFile),       pointer           :: NewGaugeFile
        logical                                     :: BlockFound
        logical                                     :: AtLeastOneBlock = .false.

        !Begin-----------------------------------------------------------------

        ! Obtain gauges files for the Time Serie
do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<BeginGaugeFile>',   &
                                        block_end       = '<EndGaugeFile>',     &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  

                    AtLeastOneBlock = .true.
                    
                    call AddGaugeFile                   (NewGaugeFile)

                    call ConstructGaugeFile             (NewGaugeFile)

                    nullify(NewGaugeFile)

                else cd2
                    call Block_Unlock(Me%ObjEnterData,                          & 
                                      ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ReadGaugesFileName - ModuleInvertedBarometer - ERR10'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ReadGaugesFileName - ModuleInvertedBarometer - ERR20'
            else cd1
                stop 'ReadGaugesFileName - ModuleInvertedBarometer - ERR30'
            end if cd1

        end do do1

        if (.not. AtLeastOneBlock) then                                            
            write(*,*)
            write(*,*) 'No gauge file block is indicated in input file. '
            write(*,*) 'At least one block is needed.'
            stop 'ReadGaugesFileName - ModuleInvertedBarometer - ERR40'
        end if

    end subroutine ReadGaugesFileName

    !--------------------------------------------------------------------------

    subroutine AddGaugeFile (ObjGaugeFile)

        !Arguments-------------------------------------------------------------
        type (T_GaugeFile),     pointer             :: ObjGaugeFile

        !Local-----------------------------------------------------------------
        type (T_GaugeFile),     pointer             :: PreviousGaugeFile
        type (T_GaugeFile),     pointer             :: NewGaugeFile

        !Begin-----------------------------------------------------------------

        !Allocates new gauge file
        allocate (NewGaugeFile)
        nullify  (NewGaugeFile%Next)

        !Insert new gauge file into list and makes current point to it
        if (.not. associated(Me%FirstGaugeFile)) then
            Me%FirstGaugeFile         => NewGaugeFile
            ObjGaugeFile              => NewGaugeFile
        else
            PreviousGaugeFile         => Me%FirstGaugeFile
            ObjGaugeFile              => Me%FirstGaugeFile%Next
            do while (associated(ObjGaugeFile))
                PreviousGaugeFile     => ObjGaugeFile
                ObjGaugeFile          => ObjGaugeFile%Next
            enddo
            ObjGaugeFile              => NewGaugeFile
            PreviousGaugeFile%Next    => NewGaugeFile
        end if

    end subroutine AddGaugeFile

    !--------------------------------------------------------------------------

    subroutine ConstructGaugeFile (NewGaugeFile)

        !Arguments-------------------------------------------------------------
        type (T_GaugeFile),      pointer            :: NewGaugeFile

        !External--------------------------------------------------------------
        integer                                     :: iflag, STAT_CALL
        logical                                     :: exist
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        
        ! Obtain gauge file name
        call GetData(NewGaugeFile%Name,                                         &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'NAME',                                     &
                     ClientModule = 'InvertedBarometer',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_ .or. iflag .EQ. 0)                          &
        stop 'ConstructGaugeFile - ModuleInvertedBarometer - ERR10'

        !Verifies if file exists
        inquire(FILE = NewGaugeFile%Name, EXIST = exist)
        if (.not. exist) then
            write(*,*)
            write(*,*)'Gauge file does not exist:'//trim(NewGaugeFile%Name)
            stop 'ConstructGaugeFile - ModuleInvertedBarometer - ERR20'
        endif

        ! Obtain intended time serie name for this gauge file
        call GetData(NewGaugeFile%TimeSerieName,                                &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'TIMESERIE_NAME',                           &
                     ClientModule = 'InvertedBarometer',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_ .or. iflag .EQ. 0)                          &
        stop 'ConstructGaugeFile - ModuleInvertedBarometer - ERR30'

        ! Obtain intended time serie name for this gauge file
        call GetData(NewGaugeFile%WriteMOHIDInput,                              &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      ='WRITE_MOHID_INPUT',                         &
                     Default      = .false.,                                    &
                     ClientModule ='InvertedBarometer',                         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'ConstructGaugeFile - ModuleInvertedBarometer - ERR40'   

        if (NewGaugeFile%WriteMOHIDInput) then

            ! Obtain intended name for output gauge file to create
            call GetData(NewGaugeFile%OutputName,                               &
                         Me%ObjEnterData, iflag,                                &
                         SearchType   = FromBlock,                              &
                         keyword      ='OUTPUTNAME',                            &
                         ClientModule ='InvertedBarometer',                     &
                         STAT         = STAT_CALL)        
            if (STAT_CALL .NE. SUCCESS_ .or. iflag .EQ. 0)                      &
            stop 'ConstructGaugeFile - ModuleInvertedBarometer - ERR50'   

        endif

    end subroutine ConstructGaugeFile 

    !--------------------------------------------------------------------------

    subroutine ReadHDF5FileName

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, ClientNumber
        type (T_HDF5File),       pointer            :: NewHDF5File
        logical                                     :: BlockFound
        logical                                     :: AtLeastOneBlock = .false.

        !Begin-----------------------------------------------------------------

do2 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<BeginHDF5File>',    &
                                        block_end       = '<EndHDF5File>',      &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd3 :       if(STAT_CALL .EQ. SUCCESS_)then
cd4 :           if (BlockFound) then                                                  
                    
                    AtLeastOneBlock = .true.

                    call AddHDF5File                     (NewHDF5File)

                    call ConstructHDF5File               (NewHDF5File)

                    nullify(NewHDF5File)

                else cd4
                    call Block_Unlock(Me%ObjEnterData,                          & 
                                      ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ReadHDF5FileName - ModuleInvertedBarometer - ERR10'

                    exit do2

                end if cd4

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd3
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ReadHDF5FileName - ModuleInvertedBarometer - ERR20'
            else cd3
                stop 'ReadHDF5FileName - ModuleInvertedBarometer - ERR30'
            end if cd3

        end do do2

        if (.not. AtLeastOneBlock) then                                            
            write(*,*)
            write(*,*) 'No HDF5 file block is indicated in input file. '
            write(*,*) 'At least one block is needed. '
            stop 'ReadHDF5FileName - ModuleInvertedBarometer - ERR40'
        end if

    end subroutine ReadHDF5FileName

    !--------------------------------------------------------------------------

    subroutine AddHDF5File(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type (T_HDF5File),     pointer              :: ObjHDF5File

        !Local-----------------------------------------------------------------
        type (T_HDF5File),     pointer              :: PreviousHDF5File
        type (T_HDF5File),     pointer              :: NewHDF5File

        !Begin-----------------------------------------------------------------

        !Allocates new HDF5File
        allocate (NewHDF5File)
        nullify  (NewHDF5File%Next)

        !Insert new file into list and makes current point to it
        if (.not. associated(Me%FirstHDF5File)) then
            Me%FirstHDF5File         => NewHDF5File
            ObjHDF5File              => NewHDF5File
        else
            PreviousHDF5File         => Me%FirstHDF5File
            ObjHDF5File              => Me%FirstHDF5File%Next
            do while (associated(ObjHDF5File))
                PreviousHDF5File     => ObjHDF5File
                ObjHDF5File          => ObjHDF5File%Next
            enddo
            ObjHDF5File              => NewHDF5File
            PreviousHDF5File%Next    => NewHDF5File
        end if

    end subroutine AddHDF5File

    !--------------------------------------------------------------------------

    subroutine ConstructHDF5File(NewHDF5File)

        !Arguments-------------------------------------------------------------
        type (T_HDF5File),      pointer             :: NewHDF5File

        !External--------------------------------------------------------------
        integer                                     :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        
        ! Obtain HDF5 file name
        call GetData(NewHDF5File%Name,                                          &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'NAME',                                     &
                     ClientModule = 'InvertedBarometer',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_ .or. iflag .EQ. 0)                          &
        stop 'ConstructHDF5File - ModuleInvertedBarometer - ERR10'

    end subroutine ConstructHDF5File

    !--------------------------------------------------------------------------

    subroutine OpenAndDateHDF5Files

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, HDF5_READ
        integer                                     :: CurrentInstant, AuxNumberInstants
        integer                                     :: HDF5PressureRank
        logical                                     :: exist, FirstTime, Relevant
        type (T_HDF5File),          pointer         :: HDF5FileX
        type(T_Time), dimension(:), pointer         :: AuxInstantsArray 
        real                                        :: Year, Month, Day, Hour 
        real                                        :: Minute, Second
        real, parameter                             :: AuxTypeReal = 8
        character(len=StringLength)                 :: HDF5PressureUnits    = null_str
        character(len=StringLength)                 :: HDF5PressureName     = null_str
      
        !Begin-----------------------------------------------------------------

        FirstTime = .true.

        HDF5FileX => Me%FirstHDF5File
        
        !In a DO cycle open all HDF5 files provided by the user
        do while (associated(HDF5FileX))

            !Verifies if file exists
            inquire(FILE = HDF5FileX%Name, EXIST = exist)
            if (.not. exist) then
                write(*,*)
                write(*,*)'HDF5 file does not exist:'//trim(HDF5FileX%Name)
                stop 'OpenAndDateHDF5Files - ModuleInvertedBarometer - ERR10'
            endif

            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

            !Open HDF5 file
            call ConstructHDF5 (HDF5FileX%HDFID, trim(HDF5FileX%Name),          &
                                HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'OpenAndDateHDF5Files - ModuleInvertedBarometer - ERR20'

            !Obtain start and end times of HDF5 file
            !(obtain number of instants) 
            call GetHDF5GroupNumberOfItems(HDF5FileX%HDFID, "/Time",            &
                                           HDF5FileX%NumberOfInstants,          & 
                                           STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        & 
            stop 'OpenAndDateHDF5Files - ModuleInvertedBarometer - ERR30'

            !(obtain HDF5 start time)
            HDF5FileX%StartTime = HDF5TimeInstant(1, HDF5FileX)

            !(obtain HDF5 end time)
            HDF5FileX%EndTime = HDF5TimeInstant(HDF5FileX%NumberOfInstants,     & 
                                                HDF5FileX)

            !Check if file contains proper atmospheric pressure data
            call GetHDF5GroupID(HDF5FileX%HDFID,                                & 
                            "/Results/"//GetPropertyName(AtmosphericPressure_), &
                            1, HDF5PressureName, HDF5PressureUnits,             &
                            HDF5PressureRank,                                   &
                            STAT = STAT_CALL)                                

            !check if pressure data exists
            if (STAT_CALL .NE. SUCCESS_) then  
                write(*,*)'HDF5 file do not contain pressure data:'             &
                           //trim(HDF5FileX%Name)
                stop 'OpenAndDateHDF5Files - ModuleInvertedBarometer - ERR40'
            end if

            !check if rank and units are the required ones
            if (trim(HDF5PressureUnits) .NE. 'Pa') then

                write(*,*)'Pressure data units not Pascals as required.'
                stop 'OpenAndDateHDF5Files - ModuleInvertedBarometer - ERR50'

            elseif (HDF5PressureRank .NE. 2) then

                write(*,*)
                write(*,*)'HDF5 file not 2D for pressure as required:'          &
                           //trim(HDF5FileX%Name)
                write(*,*)'Sea level atmospheric pressure is required.'
                stop 'OpenAndDateHDF5Files - ModuleInvertedBarometer - ERR60'

            endif 

            !Check if the HDF5 file is relevant for the time series
            call HDF5Evaluator(HDF5FileX, Relevant)

            !If HDF5 file is relevant then obtain key parameters and instants
            if (Relevant) then

                !Get useful time information from file:
                !Set instant array for this file 
                allocate(AuxInstantsArray(1:HDF5FileX%NumberOfInstants))

                !Fill array with instants
                do CurrentInstant = 1, HDF5FileX%NumberOfInstants

                   AuxInstantsArray(CurrentInstant) =                           &
                        HDF5TimeInstant(CurrentInstant, HDF5FileX)

                end do

                !Get start and end instants for this file
                !select time window begin
                do CurrentInstant = 1, HDF5FileX%NumberOfInstants

                    if (AuxInstantsArray(CurrentInstant)                        &
                        .ge. HDF5FileX%StartFieldTime) then

                        HDF5FileX%StartInstant = CurrentInstant
                        HDF5FileX%StartFieldTime =                              &
                            HDF5TimeInstant(CurrentInstant, HDF5FileX)
                
                        exit

                    end if

                end do
    
                !select time window end
                do CurrentInstant = HDF5FileX%StartInstant,                     &
                                            HDF5FileX%NumberOfInstants

                    if (AuxInstantsArray(CurrentInstant)                        &
                        .eq. HDF5FileX%EndFieldTime) then

                        HDF5FileX%EndInstant = (CurrentInstant)
                        HDF5FileX%EndFieldTime =                                &
                            HDF5TimeInstant(CurrentInstant, HDF5FileX)

                        exit

                    end if

                    if (AuxInstantsArray(CurrentInstant)                        &
                        .gt. HDF5FileX%EndFieldTime) then

                        HDF5FileX%EndInstant = (CurrentInstant-1)
                        HDF5FileX%EndFieldTime =                                &
                            HDF5TimeInstant(CurrentInstant-1, HDF5FileX)

                        exit 

                    end if

                end do

                !Check to see the presence of only one time in the time serie
                if (HDF5FileX%StartFieldTime .eq. HDF5FileX%EndFieldTime) then
                    if (Me%StartTSTime >= HDF5FileX%StartTime .AND.             &
                        Me%EndTSTime <= HDF5FileX%EndTime) then          
                        write(*,*)
                        write(*,*) 'Time series has only one time:' 
100                     format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x,     &
                                f3.0, 1x, f3.0)
                        call ExtractDate(HDF5FileX%StartFieldTime, Year, Month, & 
                             Day, Hour, Minute, Second)
                        write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
                        write(*,*) 'This is not allowed.'
                        stop 'OpenAndDateHDF5Files - ModuleInvertedBarometer - ERR70'
                    end if
                end if

                AuxNumberInstants = HDF5FileX%NumberOfInstants

                HDF5FileX%NumberOfInstants = HDF5FileX%EndInstant -             &
                                             HDF5FileX%StartInstant + 1

                !get instant times 
                call ObtainInstantsTimes(HDF5FileX)               

                !Add file to list of relevant files
                call AddTSHDF5File(HDF5FileX)

                !next run of cycle is not the first one                
                FirstTime = .false.
           
                deallocate(AuxInstantsArray)
                nullify(AuxInstantsArray)

            end if

            call killhdf5(HDF5FileX%HDFID)           

            HDF5FileX => HDF5FileX%Next

        end do

        !assume VariableDT
        Me%VariableDT = .True.

        !For each HDF5 file needed for the time serie
        HDF5FileX => Me%FirstTSHDF5File

        !Check for data lacking at the beginning of period
        do while(associated(HDF5FileX))

            !Get instants' times and put them in list 
            call AddTSInstantsTimes(HDF5FileX)

            HDF5FileX => HDF5FileX%Next            

        end do

        if (FirstTime) then

            !No HDF5 file is provided with suitable data
            call ExtractDate(Me%StartTSTime, Year, Month, Day, Hour,            &  
                             Minute, Second)        
            write(*,*)
            write(*,*) 'Data lacking from'      
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
            write(*,*) 'to'      
            call ExtractDate(Me%EndTSTime, Year, Month, Day, Hour, Minute, Second)
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
            write(*,*) 'There are not requested data in the HDF5 files provided.'      
            stop 'OpenAndDateHDF5Files - ModuleInvertedBarometer - ERR80'                

        else         
            !Check data lacking at the beginning of time serie
            if (Me%FirstTimeSerieTime%Time > Me%StartTSTime) then
                call ExtractDate(Me%StartTSTime, Year, Month, Day, Hour,        &  
                                 Minute, Second)        
                write(*,*)
                write(*,*) 'Data lacking from'      
                write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
                write(*,*) 'to'      
                call ExtractDate(Me%FirstTimeSerieTime%Time, Year, Month, Day,  &
                                 Hour, Minute, Second)
                write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
            end if
        end if

        !Calculates the size of the buffer
        if (Me%TotalOutPutsNumber * AuxTypeReal > Me%MaxBufferSize) then
            Me%BufferSize  = int(Me%MaxBufferSize / (AuxTypeReal))
        else
            Me%BufferSize  = Me%TotalOutPutsNumber
        endif

    end subroutine OpenAndDateHDF5Files

    !--------------------------------------------------------------------------

    subroutine HDF5Evaluator(HDF5FileX, Relevant)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX
        logical                                     :: Relevant

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_HDF5File), pointer                   :: OtherFile      

        !Begin-----------------------------------------------------------------

        Relevant = .FALSE.

        !Check if there is a file with same period
        call SamePeriodHDF5(HDF5FileX, OtherFile, STAT_CALL)

        if (STAT_CALL .EQ. SUCCESS_) then 
            write(*,*)
            write(*,*)'Same period is contained in two HDF5 files:'
            write(*,*) HDF5FileX%Name
            write(*,*) OtherFile%Name
            write(*,*) 'This is not allowed.'
            stop 'HDF5Evaluator - ModuleInvertedBarometer - ERR10'
        end if

        !See if the file is to be considered for time serie 
        !Check start and end time
        if ((HDF5FileX%EndTime >= Me%StartTSTime) .and.                         &
            (HDF5FileX%EndTime <= Me%EndTSTime)) then

            !End time serie time is between start and end times of file
            if (HDF5FileX%StartTime < Me%StartTSTime) then

                !Start time serie time is after start time of file
                HDF5FileX%StartFieldTime = Me%StartTSTime
                HDF5FileX%EndFieldTime = HDF5FileX%EndTime

                Relevant = .TRUE.

            else 

                !Start time serie time is before start time of file
                HDF5FileX%StartFieldTime = HDF5FileX%StartTime
                HDF5FileX%EndFieldTime = HDF5FileX%EndTime

                Relevant = .TRUE.

            end if

        else if ((HDF5FileX%StartTime >= Me%StartTSTime) .and.                  &
                 (HDF5FileX%StartTime <= Me%EndTSTime)) then

            !End time serie time is before end time of file
            HDF5FileX%StartFieldTime = HDF5FileX%StartTime
            HDF5FileX%EndFieldTime = Me%EndTSTime

            Relevant = .TRUE.

        else if ((HDF5FileX%StartTime < Me%StartTSTime) .and.                   &
                 (HDF5FileX%EndTime > Me%EndTSTime)) then

            !Time serie is contained in file
            HDF5FileX%StartFieldTime = Me%StartTSTime
            HDF5FileX%EndFieldTime = Me%EndTSTime

            Relevant = .TRUE.

        end if

    end subroutine HDF5Evaluator

    !--------------------------------------------------------------------------

    subroutine SamePeriodHDF5(HDF5FileX, HDF5FileAux, STAT)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX
        type(T_HDF5File), pointer                   :: HDF5FileAux
        integer, intent(OUT)                        :: STAT

        !Local-----------------------------------------------------------------
      
        !Begin-----------------------------------------------------------------

        STAT = UNKNOWN_

        HDF5FileAux => Me%FirstHDF5File

        do while (associated(HDF5FileAux))

            if (HDF5FileX%Name .NE. HDF5FileAux%Name) then 
            !(not the same HDF5 file)

                !Check if the same period is in more than one file
                if (((HDF5FileX%StartTime >= HDF5FileAux%StartTime)             &
                    .and. (HDF5FileX%StartTime <= HDF5FileAux%EndTime))         &
                    .or. ((HDF5FileX%EndTime >= HDF5FileAux%StartTime)          &
                    .and. (HDF5FileX%EndTime <= HDF5FileAux%EndTime))           &
                    .or. ((HDF5FileX%StartTime <= HDF5FileAux%StartTime)        &
                    .and. (HDF5FileX%EndTime >= HDF5FileAux%EndTime))) then

                    STAT = SUCCESS_

                end if

            end if

            HDF5FileAux => HDF5FileAux%Next

        end do

    end subroutine SamePeriodHDF5

    !--------------------------------------------------------------------------

    type(T_Time) function HDF5TimeInstant(Instant, ObjHDF5File)

        !Arguments-------------------------------------------------------------
        integer                                     :: Instant
        type(T_HDF5File),   pointer                 :: ObjHDF5File
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer                 :: TimeVector

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (ObjHDF5File%HDFID, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = ObjHDF5File%HDFID,                &
                             GroupName      = "/Time",                          &
                             Name           = "Time",                           &
                             Array1D        = TimeVector,                       &
                             OutputNumber   = Instant,                          &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'HDF5TimeInstant - ModuleInvertedBarometer - ERR10'

        call SetDate(HDF5TimeInstant, Year  = TimeVector(1),                    &
                     Month  = TimeVector(2), Day      = TimeVector(3),          &
                     Hour   = TimeVector(4), Minute   = TimeVector(5),          &
                     Second = TimeVector(6))

        deallocate(TimeVector)
        nullify   (TimeVector)

    end function HDF5TimeInstant

    !--------------------------------------------------------------------------

    subroutine ObtainInstantsTimes(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: ObjHDF5File
        
        !Local-----------------------------------------------------------------
        type (T_TimeSerieTime), pointer             :: NewTime, ObjTime
        type (T_TimeSerieTime), pointer             :: PreviousTime
        integer                                     :: CurrentInstant
        
        !Begin-----------------------------------------------------------------
        
        do CurrentInstant = ObjHDF5File%StartInstant, ObjHDF5File%EndInstant

            !Allocates new instance
            allocate (NewTime)
            nullify  (NewTime%Next)

            NewTime%Time = HDF5TimeInstant(CurrentInstant, ObjHDF5File)

            !Insert New Instance into list and makes Current point to it
            if (.not. associated(ObjHDF5File%FirstInstantTime)) then
            !FirstField should be the same for all HDF5 files 
                ObjHDF5File%FirstInstantTime    => NewTime
                ObjTime                         => NewTime
            else
                PreviousTime                    => ObjHDF5File%FirstInstantTime
                ObjTime                         => ObjHDF5File%FirstInstantTime%Next
                do while (associated(ObjTime))
                    PreviousTime                => ObjTime
                    ObjTime                     => ObjTime%Next
                enddo
                ObjTime                         => NewTime
                PreviousTime%Next               => NewTime
            endif

        end do 

    end subroutine ObtainInstantsTimes

    !--------------------------------------------------------------------------

    subroutine AddTSHDF5File(HDF5FileX)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX

        !Local-----------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileAux
        type(T_HDF5File), pointer                   :: PreviousHDF5File, LastHDF5File
      
        !Begin-----------------------------------------------------------------

        if (.not. associated(Me%FirstTSHDF5File)) then

            call CreateTSHDF5File(Me%FirstTSHDF5File, HDF5FileX)
            call CreateTSHDF5File(Me%LastTSHDF5File, HDF5FileX)
            deallocate(Me%FirstTSHDF5File%Next) 
            nullify(Me%FirstTSHDF5File%Next)

        else

            if (HDF5FileX%StartFieldTime < Me%FirstTSHDF5File%StartFieldTime) then
                !current file should be the first file in list

                !save the previous list 
                allocate(HDF5FileAux)
                call CreateTSHDF5File(HDF5FileAux, Me%FirstTSHDF5File)
                HDF5FileAux%Next => Me%FirstTSHDF5File%Next              

                !make the first element in list of relevant files to be current file
                call CreateTSHDF5File(Me%FirstTSHDF5File, HDF5FileX)
                Me%FirstTSHDF5File%Next => HDF5FileAux

            else
                !check next files in list

                !locate previous file in the first file
                allocate(PreviousHDF5File)
                PreviousHDF5File => Me%FirstTSHDF5File                   

                do while(associated(PreviousHDF5File))

                    if (.not. associated(PreviousHDF5File%Next)) then
        
                        !current file is the last file in the list of relevant files
                        call CreateTSHDF5File(PreviousHDF5File%Next, HDF5FileX)
                        allocate(LastHDF5File)
                        LastHDF5File => PreviousHDF5File%Next
                        deallocate(LastHDF5File%Next)
                        nullify(LastHDF5File%Next)
                        call CreateTSHDF5File(Me%LastTSHDF5File, HDF5FileX)

                        !current file was added to list
                        exit

                    else

                        !check if current file should be located before the next file
                        if (HDF5FileX%StartFieldTime <                          & 
                            PreviousHDF5File%Next%StartFieldTime) then
                            !current file should be located before next file

                            !save the previous list beginning in PreviousHDF5File%Next
                            allocate(LastHDF5File)
                            LastHDF5File => PreviousHDF5File%Next
                            allocate(HDF5FileAux)

                            call CreateTSHDF5File(HDF5FileAux, LastHDF5File)
                            HDF5FileAux%Next => LastHDF5File%Next

                            call CreateTSHDF5File(LastHDF5File, HDF5FileX)

                            PreviousHDF5File%Next => LastHDF5File
                            LastHDF5File%Next => HDF5FileAux

                            !current file was added to list
                            exit

                        end if

                    end if

                    !check next file in list
                    PreviousHDF5File => PreviousHDF5File%Next     

                end do

            end if

        end if

    end subroutine AddTSHDF5File

    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Grid

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
      
        !Begin-----------------------------------------------------------------

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid,                      & 
                                     Me%GridFile, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             &
            stop 'ConstructHDF5Grid - ModuleInvertedBarometer - ERR10'

        call GetHorizontalGridSize(Me%ObjHorizontalGrid,                        &
                                   WorkSize = Me%WorkSize2D,                    &
                                   Size     = Me%Size2D,                        &
                                   STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             &
            stop 'ConstructHDF5Grid - ModuleInvertedBarometer - ERR20'

        call GetGridCoordType(Me%ObjHorizontalGrid,                             &
                              CoordType = Me%GridCoordType,                     &
                              STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             &
            stop 'ConstructHDF5Grid - ModuleInvertedBarometer - ERR30'

        ! Construct polygon for the whole grid area
        call ConstructGridPolygon

    end subroutine ConstructHDF5Grid

    !--------------------------------------------------------------------------

    subroutine CreateTSHDF5File(HDF5FileNew, HDF5FileX)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileNew
        type(T_HDF5File), pointer                   :: HDF5FileX

        !Local-----------------------------------------------------------------
      
        !Begin-----------------------------------------------------------------

        !This subroutine atributes the values of fields of HDF5File to another HDF5File
        allocate(HDF5FileNew)

        HDF5FileNew%Name             =  HDF5FileX%Name
        HDF5FileNew%StartTime        =  HDF5FileX%StartTime
        HDF5FileNew%EndTime          =  HDF5FileX%EndTime
        HDF5FileNew%StartFieldTime   =  HDF5FileX%StartFieldTime
        HDF5FileNew%EndFieldTime     =  HDF5FileX%EndFieldTime
        HDF5FileNew%NumberOfInstants =  HDF5FileX%NumberOfInstants
        HDF5FileNew%StartInstant     =  HDF5FileX%StartInstant
        HDF5FileNew%EndInstant       =  HDF5FileX%EndInstant
        HDF5FileNew%FirstInstantTime => HDF5FileX%FirstInstantTime

    end subroutine CreateTSHDF5File

    !--------------------------------------------------------------------------

    subroutine AddTSInstantsTimes(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File),       pointer             :: ObjHDF5File
        
        !Local-----------------------------------------------------------------
        type (T_TimeSerieTime), pointer             :: NewTime, ObjTime
        type (T_TimeSerieTime), pointer             :: PreviousTime, ObjInstantTime
        
        !Begin-----------------------------------------------------------------

        ObjInstantTime => ObjHDF5File%FirstInstantTime
        
        do while (associated(ObjInstantTime))

            !Allocates new instance
            allocate (NewTime)
            nullify  (NewTime%Next)

            NewTime%Time = ObjInstantTime%Time

            !Insert New Instance into list and makes Current point to it
            if (.not. associated(Me%FirstTimeSerieTime)) then
            !FirstField should be the same for all HDF5 files 
                Me%FirstTimeSerieTime   => NewTime
                ObjTime                 => NewTime
            else
                PreviousTime            => Me%FirstTimeSerieTime
                ObjTime                 => Me%FirstTimeSerieTime%Next
                do while (associated(ObjTime))
                    PreviousTime        => ObjTime
                    ObjTime             => ObjTime%Next
                enddo
                ObjTime                 => NewTime
                PreviousTime%Next       => NewTime
            endif

            Me%TotalOutPutsNumber = Me%TotalOutPutsNumber + 1

            ObjInstantTime => ObjInstantTime%Next          

        end do 

    end subroutine AddTSInstantsTimes

    !--------------------------------------------------------------------------

    subroutine ConstructGaugeList

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        type(T_GaugeFile),      pointer             :: GaugeFileX
        logical                                     :: BlockFound
        logical                                     :: AtLeastOneBlock = .false.
        integer                                     :: STAT_CALL, ClientNumber
        integer                                     :: TSCount = 0
        type(T_GaugeLocation),  pointer             :: NewGaugeLocation
        
        !Begin-----------------------------------------------------------------

        GaugeFileX => Me%FirstGaugeFile
      
        !In a DO cycle open all gauge files provided by the user
        do while (associated(GaugeFileX))

            !Open file
            call ConstructEnterData (GaugeFileX%ObjEnterData, GaugeFileX%Name,  & 
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ConstructGaugeList - ModuleInvertedBarometer - ERR10'

            ! Obtain gauge locations for Time Serie
do3 :       do
                call ExtractBlockFromBuffer(GaugeFileX%ObjEnterData,            &
                                            ClientNumber    = ClientNumber,     &
                                            block_begin     = '<begingauge>',   &
                                            block_end       = '<endgauge>',     &
                                            BlockFound      = BlockFound,       &
                                            STAT            = STAT_CALL)
cd5 :           if(STAT_CALL .EQ. SUCCESS_)then
cd6 :               if (BlockFound) then                                                  

                        AtLeastOneBlock = .true.
                    
                        TSCount = TSCount + 1

                        call AddGaugeLocation(GaugeFileX, NewGaugeLocation)

                        call ConstructGaugeLocation(GaugeFileX,                 & 
                                                    NewGaugeLocation, TSCount)

                        ! Check if gauge is inside grid polygon
                        call CheckInsideGridPolygon (GaugeFileX, NewGaugeLocation)

                        nullify(NewGaugeLocation)

                    else cd6
                        call Block_Unlock(GaugeFileX%ObjEnterData,              & 
                                          ClientNumber, STAT = STAT_CALL) 

                        if (STAT_CALL .NE. SUCCESS_)                            &
                        stop 'ConstructGaugeList - ModuleInvertedBarometer - ERR20'

                        exit do3
                    end if cd6

                else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd5
                    write(*,*)  
                    write(*,*)'Error calling ExtractBlockFromBuffer.'
                    stop 'ConstructGaugeList - ModuleInvertedBarometer - ERR30'
                else cd5
                    stop 'ConstructGaugeList - ModuleInvertedBarometer - ERR40'
                end if cd5

            end do do3

            if (.not. AtLeastOneBlock) then                                            
                write(*,*)
                write(*,*)'No gauge block is indicated in gauge file: '         &
                          //trim(GaugeFileX%Name)
                write(*,*)'At least one block is needed.'
                stop 'ConstructGaugeList - ModuleInvertedBarometer - ERR50'
            end if

            Me%TotalGauges = Me%TotalGauges + GaugeFileX%TotalGauges

            TSCount = 0 

            call KillEnterData (GaugeFileX%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructGaugeList - ModuleInvertedBarometer - ERR60'
            end if

            GaugeFileX => GaugeFileX%Next

        enddo

        if (Me%TotalGauges .eq. 0) then                                            
            write(*,*)
            write(*,*)'All gauge locations are outside HDF5 grid limits.'
            write(*,*)'No data is available to calculate gauge mean level.'
            stop 'ConstructGaugeList - ModuleInvertedBarometer - ERR70'
        end if

        !Discard Gridpolygon
        deallocate(Me%GridPolygon%VerticesF, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ConstructGaugeList - ModuleInvertedBarometer - ERR80'  
        nullify(Me%GridPolygon%VerticesF)
        deallocate(Me%GridPolygon)
        nullify(Me%GridPolygon)
    
    end subroutine ConstructGaugeList

    !--------------------------------------------------------------------------

    subroutine AddGaugeLocation(ObjGaugeFile,ObjGaugeLocation)

        !Arguments-------------------------------------------------------------
        type(T_GaugeLocation),     pointer          :: ObjGaugeLocation
        type(T_GaugeFile),         pointer          :: ObjGaugeFile

        !Local-----------------------------------------------------------------
        type (T_GaugeLocation),    pointer          :: PreviousGaugeLocation
        type (T_GaugeLocation),    pointer          :: NewGaugeLocation

        !Begin-----------------------------------------------------------------

        !Allocates new gauge file
        allocate (NewGaugeLocation)
        nullify  (NewGaugeLocation%Next)

        !Insert new gauge file into list and makes current point to it
        if (.not. associated(ObjGaugeFile%FirstGaugeLocation)) then
            ObjGaugeFile%FirstGaugeLocation => NewGaugeLocation
            ObjGaugeLocation                => NewGaugeLocation
        else
            PreviousGaugeLocation           => ObjGaugeFile%FirstGaugeLocation
            ObjGaugeLocation                => ObjGaugeFile%FirstGaugeLocation%Next
            do while (associated(ObjGaugeLocation))
                PreviousGaugeLocation       => ObjGaugeLocation
                ObjGaugeLocation            => ObjGaugeLocation%Next
            enddo
            ObjGaugeLocation                => NewGaugeLocation
            PreviousGaugeLocation%Next      => NewGaugeLocation
        end if

    end subroutine AddGaugeLocation

    !--------------------------------------------------------------------------

    subroutine ConstructGaugeLocation(ObjGaugeFile, NewGaugeLocation, TSCount)

        !Arguments-------------------------------------------------------------
        type(T_GaugeFile),      pointer             :: ObjGaugeFile
        type(T_GaugeLocation),  pointer             :: NewGaugeLocation 
        integer                                     :: TSCount

        !External--------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        character(len=PathLength)                   :: AuxTSFileName        = null_str    
        character(StringLength)                     :: AuxNum, AuxGaugeName
        integer                                     :: control = 0
        integer                                     :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------
        
        ! Obtain gauge name
        call GetData(AuxGaugeName,                                              &
                     ObjGaugeFile%ObjEnterData, iflag,                          &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'NAME',                                     &
                     default      = 'Gauge X',                                  &
                     text         = 'Gauge name missing',                       &
                     ClientModule = 'InvertedBarometer',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
        stop 'ConstructGaugeLocation - ModuleInvertedBarometer - ERR10'

        !Construct gauge time serie name
        write(AuxNum, fmt=*) TSCount
        
        AuxTSFileName = trim(ObjGaugeFile%TimeSerieName)//trim(AuxGaugeName)    &
                        //trim(adjustl(AuxNum))
        NewGaugeLocation%Name = trim(AuxTSFileName)//".bts"

        ! Obtain gauge longitude
        call GetData(NewGaugeLocation%Longitude,                                &
                     ObjGaugeFile%ObjEnterData, iflag,                          &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'LONGITUDE',                                &
                     text         = 'Gauge Longitude missing',                  &
                     ClientModule = 'InvertedBarometer',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
        stop 'ConstructGaugeLocation - ModuleInvertedBarometer - ERR20'

        if ((iflag > 0).and.(iflag < 3)) then
            write (*,*)
            write (*,*) '3 values for LONGITUDE are expected!'
            write (*,*) 'Gauge file: '//trim(ObjGaugeFile%Name)
            write (*,*) 'Gauge location name: '//trim(NewGaugeLocation%Name)
            stop 'ConstructGaugeLocation - ModuleInvertedBarometer - ERR30'
        endif
        if (iflag == 3) control = 1

        ! Obtain gauge latitude
        call GetData(NewGaugeLocation%Latitude,                                 &
                     ObjGaugeFile%ObjEnterData, iflag,                          &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'LATITUDE',                                 &
                     text         = 'Gauge Latitude missing',                   &
                     ClientModule = 'InvertedBarometer',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
        stop 'ConstructGaugeLocation - ModuleInvertedBarometer - ERR40'

if5 :   if (iflag < 3) then
            write (*,*)
            write (*,*) '3 values for LATITUDE were expected!'
            write (*,*) 'Gauge file: '//trim(ObjGaugeFile%Name)
            write (*,*) 'Gauge location name: '//trim(NewGaugeLocation%Name)
            stop 'ConstructGaugeLocation - ModuleInvertedBarometer - ERR50'
        end if if5

        call Calc_Decimal_geo_coord(NewGaugeLocation)

        call GetData(NewGaugeLocation%Metric_X,                                &
                     ObjGaugeFile%ObjEnterData, iflag,                          &
                     keyword    = 'METRIC_X',                                   &
                     default    = FillValueReal,                                &
                     SearchType = FromBlock,                                    &
                     ClientModule ='InvertedBarometer',                         &
                     STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
        stop 'ConstructGaugeLocation - ModuleInvertedBarometer - ERR60'

if1 :   if (iflag == 1) then 

            call GetData(NewGaugeLocation%Metric_Y,                             &
                         ObjGaugeFile%ObjEnterData, iflag,                      &
                         keyword    = 'METRIC_Y',                               &
                         default    = FillValueReal,                            &
                         SearchType = FromBlock,                                &
                         ClientModule ='InvertedBarometer',                     &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ConstructGaugeLocation - ModuleInvertedBarometer - ERR70'

            if (iflag == 1) then

                NewGaugeLocation%DefinedMetric = .true.
                if (control == 0) control = 1

            endif

        endif if1

if2 :   if (control == 0) then
            write (*,*) 'No gauge location provided'
            stop 'ConstructGaugeLocation - ModuleInvertedBarometer - ERR80'
        else 

            if (Me%GridCoordType /= 1 .and. Me%GridCoordType /= 4) then

                !Required metric coordinates
                if (.not. NewGaugeLocation%DefinedMetric) then
                    write (*,*)    
                    write (*,*) 'Gauge location not provided in metric coordinates:'
                    write (*,*) trim(AuxTSFileName)
                    write (*,*) 'In gauge file:'
                    write (*,*) trim(ObjGaugeFile%Name)
                    write (*,*) 'This is required because HDF5 grid in metric coordinates.'
                    stop 'ConstructGaugeLocation - ModuleInvertedBarometer - ERR90'
                end if

            end if

        end if if2

        ! Obtain gauge reference level
        call GetData(NewGaugeLocation%ReferenceLevel,                           &
                     ObjGaugeFile%ObjEnterData, iflag,                          &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'REF_LEVEL',                                &
                     default    = 0.,                                           &
                     text       = 'Warning. Gauge Level missing. Zero assumed', &
                     ClientModule = 'InvertedBarometer',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
        stop 'ConstructGaugeLocation - ModuleInvertedBarometer - ER100'

    end subroutine ConstructGaugeLocation

    !--------------------------------------------------------------------------

    subroutine CheckInsideGridPolygon(ObjGaugeFile, ObjGaugeLocation)

        !Arguments-------------------------------------------------------------
        type(T_GaugeFile),      pointer             :: ObjGaugeFile
        type(T_GaugeLocation),  pointer             :: ObjGaugeLocation 

        !External--------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        type(T_PointF),         pointer             :: AuxPoint

        !Begin-----------------------------------------------------------------

        nullify (AuxPoint)
        allocate(AuxPoint)

        ! Construct point
        if (Me%GridCoordType == 1 .or. Me%GridCoordType == 4) then
            AuxPoint%X = ObjGaugeLocation%DecimalLongitude
            AuxPoint%Y = ObjGaugeLocation%DecimalLatitude
        else
            AuxPoint%X = ObjGaugeLocation%Metric_X
            AuxPoint%Y = ObjGaugeLocation%Metric_Y
        endif

        ! Check if inside grid cell
        if (.not. IsPointInsidePolygon(AuxPoint, Me%GridPolygon)) then
            
            write (*,*)
            write (*,*) 'Gauge location outside HDF5 grid'
            write (*,*) 'Gauge location name: '//trim(ObjGaugeLocation%Name)
            write (*,*) 'X: ', AuxPoint%X
            write (*,*) 'Y: ', AuxPoint%Y

        else

            ObjGaugeFile%TotalGauges = ObjGaugeFile%TotalGauges + 1

        endif

    end subroutine CheckInsideGridPolygon

    !--------------------------------------------------------------------------

    subroutine Calc_Decimal_geo_coord(ObjGaugeLocation)

        !Arguments------------------------------------------------------------
        type (T_GaugeLocation), pointer             :: ObjGaugeLocation 

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        if(ObjGaugeLocation%Latitude (1) .ge. 0) then
            ObjGaugeLocation%DecimalLatitude = ObjGaugeLocation%latitude (1) +          &
                                               ObjGaugeLocation%latitude (2)/60.0 +     &
                                               ObjGaugeLocation%latitude (3)/3600.0
                                               

        else
            ObjGaugeLocation%DecimalLatitude = ObjGaugeLocation%Latitude (1) -          &
                                               ObjGaugeLocation%Latitude (2)/60.0 -     &
                                               ObjGaugeLocation%latitude (3)/3600.0

        endif

        if(ObjGaugeLocation%Longitude (1) .ge. 0)then

            ObjGaugeLocation%DecimalLongitude = ObjGaugeLocation%Longitude(1) +         &
                                                ObjGaugeLocation%Longitude(2)/60.0 +    &
                                                ObjGaugeLocation%longitude(3)/3600.0

        else
            
            ObjGaugeLocation%DecimalLongitude = ObjGaugeLocation%Longitude(1) -         &
                                                ObjGaugeLocation%Longitude(2)/60.0 -    &
                                                ObjGaugeLocation%longitude(3)/3600.0

        end if

    end subroutine Calc_Decimal_geo_coord 

    !------------------------------------------------------------------------

    subroutine ConstructGridPolygon

        !Arguments-------------------------------------------------------------
        
        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: CurrentVertix, i, j
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        ! Get coordinates of vertices of grid cells
        call GetGridLatitudeLongitude(Me%ObjHorizontalGrid,                     &
                                      GridLatitudeConn = Me%GridLatitudeConn,   & 
                                      GridLongitudeConn = Me%GridLongitudeConn, &
                                      STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             & 
            stop 'ConstructGridPolygon -  ModuleInvertedBarometer - ERR10'

        allocate(Me%GridPolygon)
        nullify(Me%GridPolygon%Next)

        !Construct the vertixes of polygon (in outer cells outer vertixes):
        !get polygon number of vertixes
        Me%GridPolygon%Count = 2 * Me%Size2D%IUB +                              &
                               2 * (Me%Size2D%JUB - 2) + 1

        allocate(Me%GridPolygon%VerticesF(1:Me%GridPolygon%Count))

        CurrentVertix = 1
        
        j = Me%Size2D%JLB + 1

        do i = Me%Size2D%ILB + 1, Me%Size2D%IUB
            
            Me%GridPolygon%VerticesF(CurrentVertix)%X = Me%GridLongitudeConn(i,j)
            Me%GridPolygon%VerticesF(CurrentVertix)%Y = Me%GridLatitudeConn(i,j)

            CurrentVertix = CurrentVertix + 1

        end do

        i = Me%Size2D%IUB

        do j = Me%Size2D%JLB + 2, Me%Size2D%JUB-1
            
            Me%GridPolygon%VerticesF(CurrentVertix)%X  = Me%GridLongitudeConn(i,j)
            Me%GridPolygon%VerticesF(CurrentVertix)%Y  = Me%GridLatitudeConn(i,j)

            CurrentVertix = CurrentVertix + 1

        end do

        j = Me%Size2D%JUB

        do i = Me%Size2D%ILB + 1, Me%Size2D%IUB
            
            Me%GridPolygon%VerticesF(CurrentVertix)%X  =                        &
                Me%GridLongitudeConn(Me%Size2D%IUB+Me%Size2D%ILB-i+1,j)
            Me%GridPolygon%VerticesF(CurrentVertix)%Y  =                        &
                Me%GridLatitudeConn(Me%Size2D%IUB+Me%Size2D%ILB-i+1,j)

            CurrentVertix = CurrentVertix + 1

        end do

        i = Me%Size2D%ILB + 1

        do j = Me%Size2D%JLB + 2, Me%Size2D%JUB-1
            
            Me%GridPolygon%VerticesF(CurrentVertix)%X  =                        &
                Me%GridLongitudeConn(i,Me%Size2D%JUB+Me%Size2D%JLB-j+1)
            Me%GridPolygon%VerticesF(CurrentVertix)%Y  =                        &
                Me%GridLatitudeConn(i,Me%Size2D%JUB+Me%Size2D%JLB-j+1)

            CurrentVertix = CurrentVertix + 1

        end do

        !close polygon
        Me%GridPolygon%VerticesF(CurrentVertix)%X  = Me%GridPolygon%VerticesF(1)%X
        Me%GridPolygon%VerticesF(CurrentVertix)%Y  = Me%GridPolygon%VerticesF(1)%Y

        call SetLimitsPolygon(Me%GridPolygon)
       
    end subroutine ConstructGridPolygon

    !--------------------------------------------------------------------------

    subroutine ConstructGaugeGridLocation

        !Arguments------------------------------------------------------------

        !External--------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, STAT_CALL
        type(T_Polygon),        pointer             :: CellPolygon
        type(T_PointF),         pointer             :: AuxPoint
        type(T_GaugeFile),      pointer             :: ObjGaugeFile
        type (T_GaugeLocation), pointer             :: ObjGaugeLocation

        !Begin-----------------------------------------------------------------

        ! Obtain the grid cell where each gauge location is

        ! Get grid cell
do5 :   do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
        do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB

            ! Construct cell polygon
            call ConstructCellPolygon(i,j, CellPolygon)

            ObjGaugeFile => Me%FirstGaugeFile
    
            !In a DO cycle open all gauge files provided by the user
            do while (associated(ObjGaugeFile))      

                ObjGaugeLocation => ObjGaugeFile%FirstGaugeLocation

do4 :           do while (associated(ObjGaugeLocation))

                    if (.not. ObjGaugeLocation%CellFound) then

                        nullify (AuxPoint)
                        allocate(AuxPoint)

                        ! Construct point
                        if (Me%GridCoordType == 1 .or. Me%GridCoordType == 4) then
                            AuxPoint%X = ObjGaugeLocation%DecimalLongitude
                            AuxPoint%Y = ObjGaugeLocation%DecimalLatitude
                        else
                            AuxPoint%X = ObjGaugeLocation%Metric_X
                            AuxPoint%Y = ObjGaugeLocation%Metric_Y
                        endif

                        ! Check if inside grid cell
                        if (IsPointInsidePolygon(AuxPoint, CellPolygon)) then
                            
                            ObjGaugeLocation%I = i
                            ObjGaugeLocation%J = j
                
                            ObjGaugeLocation%CellFound = .true.
                            
                            ObjGaugeFile%FoundGauges = ObjGaugeFile%FoundGauges &
                                                       + 1    

                            Me%FoundGauges = Me%FoundGauges + 1

                        endif

                        deallocate(AuxPoint)

                    endif

                    if (ObjGaugeFile%FoundGauges .eq. ObjGaugeFile%TotalGauges) & 
                        exit do4

                    ObjGaugeLocation => ObjGaugeLocation%Next

                enddo do4

                if (Me%FoundGauges .eq. Me%TotalGauges) exit do5

                ObjGaugeFile => ObjGaugeFile%Next
            
            enddo 

            !Discard CellPolygon
            deallocate(CellPolygon%VerticesF, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'ConstructGaugeGridLocation - ModuleInvertedBarometer - ERR10'  
            nullify(CellPolygon%VerticesF)
            deallocate(CellPolygon)
            nullify(CellPolygon)

        enddo
        enddo do5

        !Discard Me%GridLongitudeConn (made here because not required by Modifier)
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%GridLongitudeConn,    & 
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructGaugeGridLocation - ModuleInvertedBarometer - ERR20'

        !Discard Me%GridLatitudeConn (made here because not required by Modifier)
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%GridLatitudeConn,     & 
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructGaugeGridLocation - ModuleInvertedBarometer - ERR30'

    end subroutine ConstructGaugeGridLocation 

    !--------------------------------------------------------------------------

    subroutine ConstructCellPolygon(i,j, CellPolygon)

        !Arguments-------------------------------------------------------------
        integer                                     :: i, j
        type(T_Polygon),    pointer                 :: CellPolygon
          
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        allocate(CellPolygon)
        nullify(CellPolygon%Next)

        !Construct the vertixes of polygon:
        !assume 5 vertixes (4 real and 1 virtual)
        CellPolygon%Count = 5

        allocate(CellPolygon%VerticesF(1:CellPolygon%Count))
       
        CellPolygon%VerticesF(1)%X  = Me%GridLongitudeConn(i,j)
        CellPolygon%VerticesF(1)%Y  = Me%GridLatitudeConn(i,j)
        CellPolygon%VerticesF(2)%X  = Me%GridLongitudeConn(i+1,j)
        CellPolygon%VerticesF(2)%Y  = Me%GridLatitudeConn(i+1,j)
        CellPolygon%VerticesF(3)%X  = Me%GridLongitudeConn(i+1,j+1)
        CellPolygon%VerticesF(3)%Y  = Me%GridLatitudeConn(i+1,j+1)
        CellPolygon%VerticesF(4)%X  = Me%GridLongitudeConn(i,j+1)
        CellPolygon%VerticesF(4)%Y  = Me%GridLatitudeConn(i,j+1)

        !close polygon
        CellPolygon%VerticesF(5)%X  = CellPolygon%VerticesF(1)%X
        CellPolygon%VerticesF(5)%Y  = CellPolygon%VerticesF(1)%Y

        call SetLimitsPolygon(CellPolygon)

    end subroutine ConstructCellPolygon

    !--------------------------------------------------------------------------
   
    subroutine SetLimitsPolygon(Polygon)

        !Arguments-------------------------------------------------------------
        type (T_Polygon)                            :: Polygon
        
        !Begin-----------------------------------------------------------------

        Polygon%Limits%Left   = minval(Polygon%VerticesF%X)
        Polygon%Limits%Right  = maxval(Polygon%VerticesF%X)
        Polygon%Limits%Bottom = minval(Polygon%VerticesF%Y)
        Polygon%Limits%Top    = maxval(Polygon%VerticesF%Y)

    end subroutine SetLimitsPolygon
    
    !--------------------------------------------------------------------------

    subroutine AllocateTimeSerieBuffer(ObjGaugeLocation)

        !Arguments-------------------------------------------------------------
        type (T_GaugeLocation), pointer             :: ObjGaugeLocation

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        real, parameter                             :: AuxTypeReal = 8

        !----------------------------------------------------------------------

        ObjGaugeLocation%BufferCount = 0

        !Allocates the TimeSerie Data Buffer
        allocate(ObjGaugeLocation%TimeSerieData(Me%BufferSize), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'AllocateTimeSerieBuffer - ModuleInvertedBarometer - ERR10'
        ObjGaugeLocation%TimeSerieData = null_real

        !Allocates the TimeSerie Time Buffer
        allocate(ObjGaugeLocation%TimeBuffer(Me%BufferSize), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'AllocateTimeSerieBuffer - ModuleInvertedBarometer - ERR20'

    end subroutine AllocateTimeSerieBuffer

    !--------------------------------------------------------------------------

    subroutine OpenOutputFiles

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: NoInformationUnit, iCh
        logical                                     :: FirstTime = .true.
        type(T_GaugeFile),      pointer             :: ObjGaugeFile
        type (T_GaugeLocation), pointer             :: ObjGaugeLocation
        character (len=StringLength)                :: ElevationName        = null_str

        !Begin-----------------------------------------------------------------

        !Call StartComputeTime for the whole Time Serie
        call StartComputeTime(Me%ObjTime, Me%FirstTSHDF5File%StartFieldTime,    &
                              Me%LastTSHDF5File%EndFieldTime,                   &
                              DT = 0.0, VariableDT = Me%VariableDT,             &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'OpenOutputFiles - ModuleInvertedBarometer - ERR10'
        !(times are common to all time series to be produced)
        !(DT is assumed null but any other value would be just as adequate because
        !this variable is not used after)

        ObjGaugeFile => Me%FirstGaugeFile
        
        !In a DO cycle open all gauge files provided by the user
        do while (associated(ObjGaugeFile))      

            ObjGaugeLocation => ObjGaugeFile%FirstGaugeLocation

            do while (associated(ObjGaugeLocation))

                if (ObjGaugeLocation%CellFound) then

                    call AllocateTimeSerieBuffer(ObjGaugeLocation)

                    ! Write header for output time serie files
                    call UnitsManager (ObjGaugeLocation%ObjTimeSerie,           &
                                       OPEN_FILE, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                  &
                        stop 'OpenOutputFiles - ModuleInvertedBarometer - ERR20'

                    open(ObjGaugeLocation%ObjTimeSerie,                         &
                         file = ObjGaugeLocation%Name, Status = 'unknown')
                    
                    !write information header
                    call WriteDataLine(ObjGaugeLocation%ObjTimeSerie,           & 
                        "Time Serie created by InvertedBarometer.exe")
                    call WriteDataLine(ObjGaugeLocation%ObjTimeSerie,           &
                        "Based on pressure HDF5 2D data on grid by "//trim(Me%GridFile))
                    
                    !write localization
                    write(ObjGaugeLocation%ObjTimeSerie, *)
                    call WriteDataLine(ObjGaugeLocation%ObjTimeSerie,           &
                                       "GRID_LOCALIZATION_I", ObjGaugeLocation%I)
                    call WriteDataLine(ObjGaugeLocation%ObjTimeSerie,           &
                                       "GRID_LOCALIZATION_J", ObjGaugeLocation%J)

                    call WriteDataLine(ObjGaugeLocation%ObjTimeSerie,           & 
                                       "SERIE_INITIAL_DATA",                    &
                                       Me%FirstTSHDF5File%StartFieldTime)        
                                                 
                    call WriteDataLine(ObjGaugeLocation%ObjTimeSerie,           &
                                       "TIME_UNITS", "SECONDS")

                    write(ObjGaugeLocation%ObjTimeSerie, *)                     &
                          "LONGITUDE               : ", ObjGaugeLocation%Longitude
                    write(ObjGaugeLocation%ObjTimeSerie, *)                     &
                          "LATITUDE                : ", ObjGaugeLocation%Latitude

                    write(ObjGaugeLocation%ObjTimeSerie, *)
                    write(ObjGaugeLocation%ObjTimeSerie, *)

                    ElevationName = GetPropertyName(WaterLevel_)

                    !Changes white spaces to underscores in elevation name 
                    !(turns input to Excel more easy)
                    do iCh = 1, len_trim(ElevationName)
                        if (ElevationName(iCh:iCh) == ' ') then
                            ElevationName(iCh:iCh) =  '_'
                        endif
                    enddo

                    write(ObjGaugeLocation%ObjTimeSerie, fmt=1000) ElevationName

 1000               format(1x,"     Seconds   YY  MM  DD  HH  MM       SS", 3x, &
                           A19, "   open_point")

                    write(ObjGaugeLocation%ObjTimeSerie, *) '<BeginTimeSerie>'

                else
    
                    ! Write gauge location reference to no information file
                    if (FirstTime) then
                        
                        call UnitsManager (NoInformationUnit, OPEN_FILE,        & 
                                       STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                        stop 'OpenOutputFiles - ModuleInvertedBarometer - ERR30'

                        FirstTime = .false.

                    endif

                    open(NoInformationUnit, file = 'Noinformation.dat',         & 
                         Status = 'unknown')

                    if (Me%GridCoordType == 1 .or. Me%GridCoordType == 4) then

                        write(NoInformationUnit, *)                             &
                              "LOCATION :  ",                                   &
                              ObjGaugeLocation%DecimalLongitude,                & 
                              ObjGaugeLocation%DecimalLatitude

                    else

                        write(NoInformationUnit, *)                             &
                              "LOCATION :  ",                                   &
                              ObjGaugeLocation%Metric_X,                        & 
                              ObjGaugeLocation%Metric_Y

                    endif
                          
                    write(NoInformationUnit, *)                                 &                         
                          "FILE :  ",                                           &
                          trim(ObjGaugeFile%Name)
                                                                                   
                endif

                ObjGaugeLocation => ObjGaugeLocation%Next

            enddo

            if (ObjGaugeFile%WriteMOHIDInput) then

                !Get unit new gauge file to be created
                call UnitsManager (ObjGaugeFile%ObjOutput,                      & 
                                   OPEN_FILE, STAT = STAT_CALL)
                                   
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'OpenOutputFiles - ModuleInvertedBarometer - ERR40'

                open(ObjGaugeFile%ObjOutput,                                    &
                     file = ObjGaugeFile%OutputName, Status = 'unknown')

                call ConstructEnterData (ObjGaugeFile%ObjEnterData,             &
                                         ObjGaugeFile%Name,                     & 
                                         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'OpenOutputFiles - ModuleInvertedBarometer - ERR50'

                !Transform data lines from gauge file
                call CreateGaugeBlocks(ObjGaugeFile)
                !(this is to read all lines of block, add new ones and 
                !put in memory)

                call UnitsManager(ObjGaugeFile%ObjOutput, CLOSE_FILE,           &
                          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'OpenOutputFiles - ModuleInvertedBarometer - ERR60'

            endif

            ObjGaugeFile => ObjGaugeFile%Next

        enddo

    end subroutine OpenOutputFiles

    !--------------------------------------------------------------------------

    subroutine CreateGaugeBlocks(ObjGaugeFile)

        !Arguments-------------------------------------------------------------
        type(T_GaugeFile),      pointer             :: ObjGaugeFile

        !Local-----------------------------------------------------------------
        logical                                     :: BlockFound      
        integer                                     :: STAT_CALL, ClientNumber
        integer                                     :: CurrentLineNumber
        integer                                     :: StartLine, EndLine
        character(len=line_length)                  :: FullBufferLine
        type (T_GaugeLocation), pointer             :: ObjGaugeLocation

        !Begin-----------------------------------------------------------------

        !Read the data lines for each location (even cell less locations)
        ObjGaugeLocation => ObjGaugeFile%FirstGaugeLocation

do6 :   do
            call ExtractBlockFromBuffer(ObjGaugeFile%ObjEnterData,              &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<begingauge>',       &
                                        block_end       = '<endgauge>',         &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd7 :       if(STAT_CALL .EQ. SUCCESS_)then
cd8 :           if (BlockFound) then                                              

                    !Get line number of block
                    call GetBlockSize(ObjGaugeFile%ObjEnterData,                &
                                      ClientNumber,                             &
                                      StartLine,                                &
                                      EndLine,                                  &
                                      FromBlock,                                &
                                      STAT = STAT_CALL)

                    !Read data lines of block
                    do CurrentLineNumber = StartLine , EndLine

                        call GetFullBufferLine(ObjGaugeFile%ObjEnterData,       &
                                               CurrentLineNumber,               &
                                               FullBufferLine,                  &
                                               STAT = STAT_CALL)
                        if(STAT_CALL .ne. SUCCESS_)                             &
                            stop 'CreateGaugeBlocks - ModuleInvertedBarometer - ERR10'
            
                        write(ObjGaugeFile%ObjOutput,*) trim(FullBufferLine)

                        if ((CurrentLineNumber == (EndLine - 1)) .and.          &
                            ObjGaugeLocation%CellFound) then

                            !Write the new keywords
                            write(ObjGaugeFile%ObjOutput,*)                     &
                                "REF_EVOLUTION : Time Serie"
                            write(ObjGaugeFile%ObjOutput,*)                     &
                                "TIME_SERIE_FILE : "//trim(ObjGaugeLocation%Name)
                            write(ObjGaugeFile%ObjOutput,*)                     &
                                "REFLEVEL_COLUMN : 8"
                            write(ObjGaugeFile%ObjOutput,*)                     &
                                "COVERED_COLUMN : 9"

                        endif

                    enddo

                    ObjGaugeLocation => ObjGaugeLocation%Next

                else cd8
                    call Block_Unlock(ObjGaugeFile%ObjEnterData,                & 
                                      ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'CreateGaugeBlocks - ModuleInvertedBarometer - ERR20'

                    exit do6
                end if cd8

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd7
                write(*,*)  
                write(*,*)'Error calling ExtractBlockFromBuffer.'
                stop 'CreateGaugeBlocks - ModuleInvertedBarometer - ERR30'
            else cd7
                stop 'CreateGaugeBlocks - ModuleInvertedBarometer - ERR40'
            end if cd7

        end do do6

    end subroutine CreateGaugeBlocks

    !--------------------------------------------------------------------------



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    
    !--------------------------------------------------------------------------



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    !--------------------------------------------------------------------------

    subroutine ModifyInvertedBarometerLevel

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                     :: Running, FirstFile
        integer                                     :: STAT_CALL, Count = 0
        integer                                     :: i, j
        real                                        :: DT, NewIBLevel
        type (T_TimeSerieTime), pointer             :: ObjTimeSerieTime
        type (T_HDF5File),      pointer             :: ObjHDF5File
        type (T_GaugeLocation), pointer             :: ObjGaugeLocation
        type (T_GaugeFile),     pointer             :: ObjGaugeFile

        !Begin-----------------------------------------------------------------

        ObjTimeSerieTime => Me%FirstTimeSerieTime

        !Cycle each relevant HDF5 file
        FirstFile =.true.

        !For each HDF5 file needed for the time serie
        ObjHDF5File => Me%FirstTSHDF5File
        
        do while(associated(ObjHDF5File))

            !Open and read relevant data
            call OpenAndReadHDF5File(FirstFile, ObjHDF5File)
            !(HDF5 data must be read in this cycle for memory to be released
            !after writing in the time serie)

            if (FirstFile) FirstFile = .false.

            Running      = .true.

            Me%CurrentTime  = ObjTimeSerieTime%Time
     
            write(*,*)
            write(*,*)'Calculating inverted barometer level correction'
            write(*,*)'and writing in time serie files...'

            do while (Running)

                Count = Count + 1 
            
                if (Count .EQ. 1) then 

                    ObjHDF5File%CurrentField => ObjHDF5File%FirstField 

                end if

                ObjGaugeFile => Me%FirstGaugeFile
        
                !For each gauge file
                do while (associated(ObjGaugeFile))      

                    ObjGaugeLocation => ObjGaugeFile%FirstGaugeLocation

                    !For each gauge location with data
                    do while (associated(ObjGaugeLocation))           

                        if (ObjGaugeLocation%CellFound) then
                        
                            i = ObjGaugeLocation%I
                            j = ObjGaugeLocation%J

                            !Calculate inverted barometer level
                            call ComputeInvertedBarometer(                      &
                                 ObjHDF5File%CurrentField%Values2D(i,j),        & 
                                 ObjGaugeLocation%ReferenceLevel, NewIBLevel)
        
                            !Write the time serie
                            call WriteIBTimeSerie(NewIBLevel, ObjGaugeLocation)

                        endif

                        ObjGaugeLocation => ObjGaugeLocation%Next

                    enddo

                    ObjGaugeFile => ObjGaugeFile%Next

                enddo

                ObjHDF5File%CurrentField => ObjHDF5File%CurrentField%Next               

                if (associated(ObjTimeSerieTime%Next)) then

                    DT = ObjTimeSerieTime%Next%Time - Me%CurrentTime

                    ObjTimeSerieTime => ObjTimeSerieTime%Next

                end if

                !Actualization of time            
                Me%CurrentTime = Me%CurrentTime + DT            

                call ActualizeCurrentTime(Me%ObjTime, DT, STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)                                     &
                stop 'ModifyInvertedBarometerLevel - ModuleInvertedBarometer - ERR10'

                !Running dependent of the last time of file
                if (Me%CurrentTime <= ObjHDF5File%EndFieldTime) then
                    Running = .true.
                else
                    Running = .false.
                end if

            end do

            Count = 0

            call KillFields(ObjHDF5File)

            ObjHDF5File => ObjHDF5File%Next           

        end do 

        call KillGaugeTimeSeries

    end subroutine ModifyInvertedBarometerLevel

    !--------------------------------------------------------------------------

    subroutine  OpenAndReadHDF5File(FirstFile, ObjHDF5File)

        !Arguments-------------------------------------------------------------
        logical                                     :: FirstFile
        type(T_HDF5File), pointer                   :: ObjHDF5File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, HDF5_READ 
        real                                        :: Year, Month, Day, Hour 
        real                                        :: Minute, Second 

        !Begin-----------------------------------------------------------------

        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Open HDF5 file
        call ConstructHDF5 (ObjHDF5File%HDFID, trim(ObjHDF5File%Name),          & 
                            HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) then
            write(*,*)
            write(*,*) 'HDF5 file cannot be opened'//ObjHDF5File%Name                
            stop 'OpenAndReadHDF5File - ModuleInvertedBarometer - ERR10'
        end if

        !Read fields data
        write(*,*)
        write(*,*)'Reading pressure from HDF5 file:'//trim(ObjHDF5File%Name) 

        if (FirstFile) then
        !it is the first file accessed

            nullify(ObjHDF5File%FirstField)
            !nullify have to be here because has to be 
            !nullified for the first file
        
        end if

        !Allocates first field for this parameter
        allocate (ObjHDF5File%FirstField)
        nullify (ObjHDF5File%FirstField)

        call ReadPressureFields(ObjHDF5File)

        !check if there are data lacking after the last data from last file 
        if (.not. associated(ObjHDF5File%Next)) then
           !(last file)
           if (ObjHDF5File%EndFieldTime < (Me%EndTSTime)) then
100             format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)
                call ExtractDate(ObjHDF5File%EndFieldTime, Year, Month, Day,    &
                                 Hour, Minute, Second)
                write(*,*)
                write(*,*) 'Data lacking from'      
                write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
                write(*,*) 'to'      
                call ExtractDate(Me%EndTSTime, Year, Month, Day, Hour, Minute,  &
                                 Second)
                write(*,fmt=100) Year, Month, Day, Hour, Minute, Second                    
            end if
        end if

        call killhdf5(ObjHDF5File%HDFID)

    end subroutine OpenAndReadHDF5File

    !--------------------------------------------------------------------------

    subroutine ReadPressureFields(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File),   pointer                 :: ObjHDF5File
          
        !Local-----------------------------------------------------------------
        integer                                     :: NewCurrentInstant, Count = 1 
        integer                                     :: CurrentInstant, STAT_CALL
        integer, dimension(7)                       :: Dimensions
        type (T_Field),     pointer                 :: NewField
        
        !Begin-----------------------------------------------------------------

        !read/copy pressure fields
        NewCurrentInstant = 0

        do CurrentInstant = ObjHDF5File%StartInstant, ObjHDF5File%EndInstant 

            NewCurrentInstant = NewCurrentInstant + 1
                
            !construct new fields
            call AddField(ObjHDF5File%FirstField, NewField)
            Count = Count + 1

            !get field ID, Rank and Dimensions
            !(this needs to be done for each instant)
            call GetHDF5GroupID(ObjHDF5File%HDFID,                              & 
                                trim("/Results/atmospheric pressure"),          &
                                CurrentInstant, NewField%Name,                  &
                                Dimensions = Dimensions,                        &
                                STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_)                                        &  
            stop 'ReadPressureFields - ModuleInvertedBarometer - ERR10'
                          
            !(assumed 2D data)
            !calculate Size2D
            Me%Size2D%ILB = 1
            Me%Size2D%IUB = Dimensions(1)
            Me%Size2D%JLB = 1
            Me%Size2D%JUB = Dimensions(2)

            !allocate field
            nullify (NewField%Values2D)
            allocate(NewField%Values2D(Me%Size2D%ILB:Me%Size2D%IUB,             &
                                       Me%Size2D%JLB:Me%Size2D%JUB))
               
            call HDF5SetLimits (ObjHDF5File%HDFID, Me%Size2D%ILB,               &
                                Me%Size2D%IUB, Me%Size2D%JLB,Me%Size2D%JUB,     &
                                STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        & 
            stop 'ReadPressureFields - ModuleInvertedBarometer - ERR20'

            !read field
            call HDF5ReadData(ObjHDF5File%HDFID,                                & 
                              trim("/Results/atmospheric pressure"),            &
                              trim(NewField%Name),                              &
                              Array2D      = NewField%Values2D,                 &
                              STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'ReadPressureFields - ModuleInvertedBarometer - ERR30'

        end do

    end subroutine ReadPressureFields

    !--------------------------------------------------------------------------

    subroutine AddField (FirstField, ObjField)

        !Arguments-------------------------------------------------------------
        type (T_Field), pointer                     :: FirstField
        type (T_Field), pointer                     :: ObjField
        
        !Local-----------------------------------------------------------------
        type (T_Field), pointer                     :: NewField, PreviousField
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewField)
        nullify  (NewField%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstField)) then
        !FirstField should be the same for all HDF5 files 
            FirstField              => NewField
            ObjField                => NewField
        else
            PreviousField           => FirstField
            ObjField                => FirstField%Next
            do while (associated(ObjField))
                PreviousField       => ObjField
                ObjField            => ObjField%Next
            enddo
            ObjField                => NewField
            PreviousField%Next      => NewField
        endif

    end subroutine AddField

    !--------------------------------------------------------------------------

    subroutine ComputeInvertedBarometer(ObjPressure, ObjRefLevel, IBLevel)

        !Arguments-------------------------------------------------------------
        real                                        :: ObjPressure
        real                                        :: ObjRefLevel
        real, intent(OUT)                           :: IBLevel
      
        !Local-----------------------------------------------------------------
        real                                        :: PressureRef = 101330.0
        !(Dorandeu, J., P. Le Traon, 1999,"Effects of Global Mean Atmospheric
        !Pressure Variations on Mean Sea Level Changes from TOPEX/Poseidon",
        !Journal of Atmospheric and Oceanic Technology, Vol. 16, No. 9, pp.
        !1279-1283; units are Pa)  
        real                                        :: LevelDif, PressureDif

        !Begin-----------------------------------------------------------------

        !(deltalevel = -(1/(density*g))*(P-Pref)) (Dorandeu, Le Traon, 1999)
        ! Because water density is not known, the rule of thumb
        ! is used: a 1 mb decrease of P below Pref results in a
        ! 1.01 cm increase in level.
        ! (Kantha, L., K. Whitmer, G. Born, 1994, "The Inverted Barometer Effect
        ! in Altimetry: A Study in the North Pacific".)

        !Get pressure difference
        PressureDif = ObjPressure - PressureRef

        !Get level difference
        LevelDif = - (0.0101 * PressureDif)/100.0

        !Get inverse barometer level
        IBLevel = ObjRefLevel + LevelDif

    end subroutine ComputeInvertedBarometer

    !--------------------------------------------------------------------------

    subroutine WriteIBTimeSerie(IBLevel, ObjGaugeLocation)

        !Arguments-------------------------------------------------------------
        real                                        :: IBLevel
        type (T_GaugeLocation), pointer             :: ObjGaugeLocation

        !Local-----------------------------------------------------------------
        type (T_Time)                               :: CurrentTime
        integer                                     :: STAT_CALL, IBC
        
        !Begin-----------------------------------------------------------------

        !(this is made similarly to time serie subroutines)
        !(the buffer is TimeSerieData in each location)

        !Stores the actual compute time
        call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'WriteIBTimeSerie - ModuleInvertedBarometer - ERR10'

        !Write time and level to buffer
        !Increments the internal buffer count
        ObjGaugeLocation%BufferCount = ObjGaugeLocation%BufferCount + 1 

        !Shorten Variable
        IBC = ObjGaugeLocation%BufferCount

        !Stores the current time
        ObjGaugeLocation%TimeBuffer(IBC) = CurrentTime

        !Stores the level
        ObjGaugeLocation%TimeSerieData(IBC) = IBLevel

        !If buffer is full then write to file
        if (ObjGaugeLocation%BufferCount  == Me%BufferSize) then
           call WriteBufferToFile(ObjGaugeLocation)
           ObjGaugeLocation%BufferCount = 0
        endif

    end subroutine WriteIBTimeSerie

    !--------------------------------------------------------------------------

    subroutine WriteBufferToFile(ObjGaugeLocation)

        !Arguments-------------------------------------------------------------
        type (T_GaugeLocation), pointer             :: ObjGaugeLocation

        !Local-----------------------------------------------------------------
        integer                                     :: iB, unit
        real                                        :: Year, Month, Day, Hour 
        real                                        :: Minute, Second
        
        !Begin-----------------------------------------------------------------

        unit = ObjGaugeLocation%ObjTimeSerie

        do iB = 1, ObjGaugeLocation%BufferCount

            !Writes date in the form of 1999 08 08 23 59 23
            call ExtractDate(ObjGaugeLocation%TimeBuffer(iB), Year = Year,      &
                             Month = Month, Day = Day, Hour = Hour,             &
                             Minute = Minute, Second = Second)

            !Writes time in seconds since the beginning of the output
            !Writes elevation in the buffer in exp format
            !Writes covered (1.0, in last column) always
            write(unit, fmt=1001) ObjGaugeLocation%TimeBuffer(iB) -             &
                                  Me%FirstTimeSerieTime%Time, int(Year),        &
                                  int(Month), int(Day), int(Hour), int(Minute), &
                                  Second, ObjGaugeLocation%TimeSerieData(iB),   &
                                  1.0
                                
        enddo

    1001 format(f13.2, 1x, i4, 2x, i2, 2x, i2, 2x, i2, 2x, i2, 2x, f7.4,        & 
                2x, e20.12e3, 3x, f3.1)

    end subroutine WriteBufferToFile

    !--------------------------------------------------------------------------

    subroutine KillFields(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File),   pointer                 :: ObjHDF5File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_Field),      pointer                 :: FieldToKill, CurrentField
        
        !Begin-----------------------------------------------------------------

        CurrentField => ObjHDF5File%FirstField

        do while(associated(CurrentField))

            FieldToKill => CurrentField
            CurrentField => CurrentField%Next

            deallocate(FieldToKill%Values2D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'KillFields - ModuleInvertedBarometer - ERR10'  
            nullify(FieldToKill%Values2D)

        end do 
        
    end subroutine KillFields

    !--------------------------------------------------------------------------

    subroutine KillGaugeTimeSeries

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type (T_GaugeLocation), pointer             :: ObjGaugeLocation
        type(T_GaugeFile),      pointer             :: ObjGaugeFile

        !Begin-----------------------------------------------------------------

        ObjGaugeFile => Me%FirstGaugeFile

        !For each gauge file
        do while (associated(ObjGaugeFile))      

            ObjGaugeLocation => ObjGaugeFile%FirstGaugeLocation

            !For each gauge location with data
            do while (associated(ObjGaugeLocation))           

                if (ObjGaugeLocation%CellFound) then

                    !Disposes buffer
                    call WriteBufferToFile(ObjGaugeLocation)

                    !Closes files
                    call CloseTimeSerieFile(ObjGaugeLocation)

                    deallocate(ObjGaugeLocation%TimeSerieData, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'KillGaugeTimeSeries - ModuleInvertedBarometer - ERR10'
                    nullify(ObjGaugeLocation%TimeSerieData)

                    deallocate(ObjGaugeLocation%TimeBuffer,                     &
                               STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'KillGaugeTimeSeries - ModuleInvertedBarometer - ERR20'
                    nullify(ObjGaugeLocation%TimeBuffer)

                endif

                ObjGaugeLocation => ObjGaugeLocation%Next

            enddo

            ObjGaugeFile => ObjGaugeFile%Next

        enddo

    end subroutine KillGaugeTimeSeries

    !--------------------------------------------------------------------------

    subroutine CloseTimeSerieFile(ObjGaugeLocation)

        !Arguments-------------------------------------------------------------
        type (T_GaugeLocation), pointer             :: ObjGaugeLocation

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        write(ObjGaugeLocation%ObjTimeSerie, *) '<EndTimeSerie>'

        call UnitsManager(ObjGaugeLocation%ObjTimeSerie, CLOSE_FILE,            &
                          STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)
            write(*,*)'Error closing Time Serie data file',                     &
                       trim(ObjGaugeLocation%Name)
            write(*,*)'CloseTimeSerieFile - ModuleInvertedBarometer - WRN01'
        endif

    end subroutine CloseTimeSerieFile

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine KillInvertedBarometerLevel

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        type(T_HDF5File),       pointer             :: HDF5FileX, HDF5FileToKill
        type(T_GaugeFile),      pointer             :: GaugeFileX, GaugeFileToKill
        type (T_GaugeLocation), pointer             :: GaugeLocationX, GaugeLocationToKill
        type (T_TimeSerieTime), pointer             :: TSTimeX, TSTimeToKill
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Kill all HDF5Files
        HDF5FileX => Me%FirstTSHDF5File

        do while(associated(HDF5FileX))  

            HDF5FileToKill  => HDF5FileX
            HDF5FileX       => HDF5FileX%Next
            call KillIndividualHDF5File(HDF5FileToKill)

        end do
        nullify(Me%FirstTSHDF5File)

        HDF5FileX=> Me%FirstHDF5File

        do while(associated(HDF5FileX))  

            HDF5FileToKill  => HDF5FileX
            HDF5FileX       => HDF5FileX%Next
            deallocate(HDF5FileToKill, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'KillInvertedBarometerLevel - ModuleInvertedBarometer - ERR10'
            nullify(HDF5FileToKill)

        end do
        nullify(Me%FirstHDF5File)

        deallocate(Me%LastTSHDF5File, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'KillInvertedBarometerLevel - ModuleInvertedBarometer - ERR20'
        nullify(Me%LastTSHDF5File)

        !Kill gauge files
        GaugeFileX => Me%FirstGaugeFile

        do while(associated(GaugeFileX))

            !Kill gauge locations
            GaugeLocationX => GaugeFileX%FirstGaugeLocation

            do while(associated(GaugeLocationX))

                GaugeLocationToKill => GaugeLocationX
                GaugeLocationX => GaugeLocationX%Next
                deallocate(GaugeLocationToKill, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                stop 'KillInvertedBarometerLevel - ModuleInvertedBarometer - ERR30'
                nullify(GaugeLocationToKill)

            end do
            nullify(GaugeFileX%FirstGaugeLocation)

            GaugeFileToKill => GaugeFileX
            GaugeFileX => GaugeFileX%Next
            deallocate(GaugeFileToKill, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'KillInvertedBarometerLevel - ModuleInvertedBarometer - ERR40'
            nullify(GaugeFileToKill)

        end do
        nullify(Me%FirstGaugeFile)

        !Kill Time Serie times
        TSTimeX => Me%FirstTimeSerieTime
        
        do while(associated(TSTimeX))

            TSTimeToKill => TSTimeX
            TSTimeX => TSTimeX%Next
            deallocate(TSTimeToKill, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'KillInvertedBarometerLevel - ModuleInvertedBarometer - ERR50'
            nullify(TSTimeToKill)

        end do
        nullify(Me%FirstTimeSerieTime)

        deallocate(Me)
        nullify(Me)

        !------------------------------------------------------------------------

    end subroutine KillInvertedBarometerLevel

    !--------------------------------------------------------------------------

    subroutine KillIndividualHDF5File(HDF5ToDispose)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File),       pointer             :: HDF5ToDispose

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_TimeSerieTime),  pointer             :: InstantTimeX, InstantTimeToKill

        !Begin-----------------------------------------------------------------

        !destroy HDF5ToDispose%FirstInstantTime and associated list
        InstantTimeX => HDF5ToDispose%FirstInstantTime
        
        do while(associated(InstantTimeX))

            InstantTimeToKill => InstantTimeX
            InstantTimeX => InstantTimeX%Next
            deallocate(InstantTimeToKill, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'KillIndividualHDF5File - ModuleInvertedBarometer - ERR10'
            nullify(InstantTimeToKill)

        end do
        nullify(HDF5ToDispose%FirstInstantTime)

        deallocate(HDF5ToDispose, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
        stop 'KillIndividualHDF5File - ModuleInvertedBarometer - ERR20'
        nullify(HDF5ToDispose)

    end subroutine KillIndividualHDF5File

    !------------------------------------------------------------------------


    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

end module ModuleInvertedBarometer









