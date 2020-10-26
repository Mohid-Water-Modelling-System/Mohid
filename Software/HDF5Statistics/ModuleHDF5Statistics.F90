!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : HDF5Statistics
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2005
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Module to create Statistics HDF5 from HDF5 files
!
!------------------------------------------------------------------------------

!DataFile
!
!   <BeginHDF5File>                                                                 
!   NAME                    : char                  [-]         !Name of HDF5 file with data for statistics 
!   START_TIME              : YYYY MM DD HH MM SS   [-]         !Date from which the statiscally analysis will start for this file
!   END_TIME                : YYYY MM DD HH MM SS   [-]         !Date from which the statiscally analysis will end for this file
!   <EndHDF5File>                                               !(specify one block for each HDF5 file to use)

!                                                                                    
!   START_TIME              : YYYY MM DD HH MM SS   [-]         !Start date of time window for statistics' calculation
!   END_TIME                : YYYY MM DD HH MM SS   [-]         !End date of time window for statistics' calculation
!
!   METHOD_STATISTIC        : 1/3                   [-]         !Statistic method for calculation (ModuleStatistics)
!                                                               !1 - Value3DStat3D
!                                                               !3 - Value2DStat2D
!                                                               (option 2 - Value3DStatLayers not yet implemented)
!
!   GLOBAL_STATISTIC        : 0/1                   [0]         !Calculate global statistics for time window
!   DAILY_STATISTIC         : 0/1                   [0]         !Calculate daily statistics for time window
!   MONTHLY_STATISTIC       : 0/1                   [0]         !Calculate monthly statistics for time window
!   SPECIFIC_HOUR_STATISTIC : 0/1                   [0]         !Calculate statistics for a specific hour of day
!   SPECIFIC_HOUR           : int                   [12]        !Specific Hour
!
!   GEOMETRIC_MEAN          : 0/1                   [0]         !Calculate also geometric statistics 
!                                                               !(average and standard deviation)
!   ACCUMULATED             : 0/1                   [0]         !Calculate also accumulated value
! 
!   <BeginParameter>                                            
!   HDF_GROUP               : char                  [-]         !Path of the HDF5 group in HDF5 file for the property
!   PROPERTY                : char                  [-]         !Property name (should be equal to the one specified 
!   <EndParameter>                                              !in ModuleGlobalData)
!                                                               !(specify one block for each property)                                                                                 
!
!   OUTPUTFILENAME          : char                  [-]         !Name of produced HDF5 file
!
!   3D_HDF5                 : 0/1                   [0]         !Input HDF5 files are 3D
!                                                               !0 - 2D ; 1 - 3D
!
!   HDF5_MAP_ITEM           : char                  [-]         !HDF5 mapping item name


Module ModuleHDF5Statistics

    use ModuleGlobalData
    use ModuleTime               
    use ModuleEnterData,         only : ConstructEnterData, KillEnterData,              &
                                        GetData, ExtractBlockFromBuffer, Block_Unlock
    use ModuleFunctions         
    use ModuleHDF5

    use ModuleStatistic,         only : ConstructStatistic, GetStatisticMethod,         &
                                        GetStatisticParameters, ModifyStatistic,        &
                                        GetStatisticLayerDef, AddStatisticLayers,       &
                                        GetStatisticLayersNumber, KillStatistic       

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartHDF5StatsCreator
    private ::      ConstructHDF5Statistics
    private ::          ReadKeywords
    private ::              ReadGlobalData
    private ::              ReadParameters
    private ::                  AddParameter
    private ::                  ConstructParameters
    private ::              ReadHDF5FileName
    private ::                  AddHDF5File
    private ::                  ConstructHDF5File  
    private ::          OpenAndDateHDF5Files
    private ::              HDF5TimeInstant
    private ::              HDF5Evaluator
    private ::                  SamePeriodHDF5
    private ::              ObtainInstantsTimes
    private ::              AddStatHDF5File
    private ::                  CreateStatHDF5File
    private ::              AddStatInstantsTimes
    private ::          OpenOutputFiles
    private ::              ConstructHDF5Grid
    private ::              KillGridFields
    private ::          ConstructStatisticGroup

    !Selector
    
    !Modifier
    public  :: ModifyHDF5Statistics
    private ::      OpenAndReadHDF5File
    private ::          ReadParameterFields
    private ::              AddField
    private ::      CalculateHDF5Statistics2D
    private ::      CalculateHDF5Statistics3D
    private ::      KillIndividualParameterFields

    !Destructor
    public  :: KillHDF5Statistics
    private ::      KillIndividualHDF5File
    private ::      KillIndividualParameter                                                     

    !Management
    
    !Interfaces----------------------------------------------------------------

    !Types---------------------------------------------------------------------
    
    ! Definition of type T_Field
    type       T_Field
        character(len=StringLength)                         :: Name
        character(len=StringLength)                         :: Units
        integer                                             :: IDNumber
        type(T_Time)                                        :: Date
        real, dimension(:    ),     pointer                 :: Values1D        
        real, dimension(:,:  ),     pointer                 :: Values2D
        real, dimension(:,:,: ),    pointer                 :: Values3D 
        type(T_Field),              pointer                 :: Next                     => null()
    end type  T_Field

    type       T_Grid
        character(len=StringLength)                         :: Name
        character(len=StringLength)                         :: AditionalName
        character(len=StringLength)                         :: Units
        real, dimension(:,:  ),     pointer                 :: RealValues2D
        real, dimension(:,:,: ),    pointer                 :: RealValues3D
        integer,    dimension(:,:,:), pointer               :: IntegerValues3D
        integer,    dimension(:,:  ), pointer               :: IntegerValues2D
        integer                                             :: Position        
        logical                                             :: Exist = .true.
    end type  T_Grid

    type       T_Statistics
         integer                                            :: ID
    end type T_Statistics

    ! Definition of type T_TimeSerieTime
    type       T_StatisticsTime
        type(T_Time)                                        :: Time
        type(T_StatisticsTime), pointer                     :: Next                     => null()
    end type  T_StatisticsTime

    ! Definition of type T_Parameter
    type       T_Parameter
        character(len=StringLength)                         :: Name
        character(len=PathLength)                           :: Group
        type(T_Field), pointer                              :: FirstField
        type(T_Field), pointer                              :: CurrentField
        integer                                             :: Rank
        character(len=StringLength)                         :: Units
        integer                                             :: LastRank
        character(len=StringLength)                         :: LastName
        type (T_Statistics)                                 :: Statistics
        type(T_Parameter), pointer                          :: Next                     => null()
    end type  T_Parameter

    ! Definition of type T_HDF5File
    type       T_HDF5File
        integer                                             :: HDFID                    = 0
        character(len=StringLength)                         :: Name
        type(T_Time)                                        :: StartTime
        type(T_Time)                                        :: EndTime
        type(T_Time)                                        :: StartFieldTime
        type(T_Time)                                        :: EndFieldTime
        logical                                             :: StartTimeDefault
        logical                                             :: EndTimeDefault
        integer                                             :: Rank
        integer                                             :: NumberOfInstants
        type(T_Time), dimension(:), pointer                 :: InstantsArray 
        integer                                             :: StartInstant, EndInstant
        type(T_StatisticsTime), pointer                     :: FirstInstantTime
        type(T_HDF5File), pointer                           :: Next                     => null()
    end type  T_HDF5File

    private :: T_HDF5Statistics
    type       T_HDF5Statistics
        type (T_Size3D)                                     :: Size, WorkSize
        character(PathLength)                               :: DataFile
        type(T_Time)                                        :: StartTime
        type(T_Time)                                        :: EndTime
        type(T_HDF5File), pointer                           :: FirstHDF5File
        type (T_Parameter), pointer                         :: FirstParameter
        integer                                             :: ParameterNumber          = 0
        integer                                             :: ObjEnterData             = 0
        type(T_HDF5File), pointer                           :: FirstStatHDF5File
        type(T_HDF5File), pointer                           :: LastStatHDF5File
        integer                                             :: ObjHDF5                  = 0
        integer                                             :: ObjStatHDF5              = 0
        integer                                             :: ObjTime                  = 0
        real                                                :: DT, RegularDT
        logical                                             :: VariableDT
        logical                                             :: CheckRelevance           = .true.
        type(T_StatisticsTime), pointer                     :: FirstStatisticsTime
        character(len=PathLength)                           :: OutputFileName
        type(T_Grid)                                        :: Bathymetry
        type(T_Grid)                                        :: ConnectionX
        type(T_Grid)                                        :: ConnectionY
        type(T_Grid)                                        :: Mapping
        type(T_Grid)                                        :: Latitude
        type(T_Grid)                                        :: Longitude
        logical                                             :: File3D    = .false.
        type(T_Time)                                        :: CurrentTime
        logical                                             :: AditionalMap = .false.
        character(len=StringLength)                         :: StatisticGroupName
    
        real,   dimension(:,:,:), pointer                   :: DZ3D    
        logical                                             :: ExistVerticalZ        
        real,   dimension(:,:,:), pointer                   :: VerticalZ
        
        type(T_HDF5Statistics), pointer                     :: Next                     => null()
    end type  T_HDF5Statistics

    !Global Module Variables
    type (T_HDF5Statistics), pointer                        :: Me                       => null()

    integer                                                 :: mHDF5Statistics_ = 0 !just to compile

    !Other Stuff
    type (T_Time)                                           :: InitialSystemTime, FinalSystemTime
    integer, dimension(8)                                   :: F95Time


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartHDF5StatsCreator(ObjHDF5StatisticsID, DataFile)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjHDF5StatisticsID 
        character(PathLength), intent(IN)               :: DataFile    

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------

        !------------------------------------------------------------------------

        !Assures nullification of the global variable
        allocate(Me)

        !Returns ID
        ObjHDF5StatisticsID          = 1

        !Atribute the name of data file            
        Me%DataFile = DataFile

        call ConstructHDF5Statistics

        !----------------------------------------------------------------------

    end subroutine StartHDF5StatsCreator

    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Statistics

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        type (T_Parameter), pointer                     :: ObjParameter
        type(T_Time)                                    :: RunStatBeginTime
        real                                            :: DT

        !Begin-----------------------------------------------------------------

        call StartupMohid ("HDF5 Statistics")

        !Gets the actual time
        call date_and_time(Values = F95Time)
        call SetDate      (InitialSystemTime, float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)

        ! Read keyword file
        call ReadKeywords

        if (Me%ParameterNumber > 0) then

            ! Open and obtain key features of the HDF5 files
            call OpenAndDateHDF5Files 

            ! Open and write header of Output Files
            call OpenOutputFiles

            ! Construct statistic group name
            call ConstructStatisticGroup

            !Construct time reducing DT from start stat time
            DT = - Me%RegularDT !(because there isn't a minus DT subroutine)

            RunStatBeginTime = Me%FirstStatHDF5File%InstantsArray(1) + DT
            
            call StartComputeTime(Me%ObjTime, InitialSystemTime, RunStatBeginTime,      &
                                  Me%LastStatHDF5File%InstantsArray                     &
                                  (Me%LastStatHDF5File%NumberOfInstants),               &
                                  Me%RegularDT, Me%VariableDT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ConstructHDF5Statistics - ModuleHDF5Statistics - ERR01'

            !Make the following inside a parameter loop
            ObjParameter => Me%FirstParameter

            do while(associated(ObjParameter))

                call ConstructStatistic (StatisticID   = ObjParameter%Statistics%ID,    &
                                         ObjTime       = Me%ObjTime,                    &
                                         ObjHDF5       = Me%ObjStatHDF5,                &
                                         Size          = Me%Size,                       &
                                         WorkSize      = Me%WorkSize,                   &
                                         DataFile      = Me%DataFile,                   &
                                         Name          = ObjParameter%Name,             &
!                                         GroupName     = trim(Me%StatisticGroupName),   &
                                         Rank          = ObjParameter%Rank,             &
                                         STAT          = STAT_CALL)                                 
                if (STAT_CALL /= SUCCESS_)                                              & 
                    stop 'ConstructHDF5Statistics - ModuleHDF5Statistics - ERR01'

                ObjParameter => ObjParameter%Next               

            end do

        end if

    end subroutine ConstructHDF5Statistics

    !--------------------------------------------------------------------------
    
    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call ConstructEnterData (Me%ObjEnterData, Me%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadKeywords - ModuleHDF5Statistics - ERR01'

        call ReadGlobalData

        call ReadParameters

        call ReadHDF5FileName

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadKeywords - ModuleHDF5Statistics - ERR02'

    end subroutine ReadKeywords

    !--------------------------------------------------------------------------

    subroutine ReadGlobalData

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------

        ! Obtain the start and end times for the statistics' calculation
        ! Start Time
        call GetData(Me%StartTime, Me%ObjEnterData, iflag,                              &
                     keyword      = 'START_TIME',                                       &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'HDF5Statistics',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
        stop 'ReadGlobalData - ModuleHDF5Statistics - ERR10'   

        ! End Time 
        call GetData(Me%EndTime,   Me%ObjEnterData, iflag,                              &
                     keyword      = 'END_TIME',                                         &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'HDF5Statistics',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
        stop 'ReadGlobalData - ModuleHDF5Statistics - ERR20'   

        ! Verifies Time Variables
        if (Me%EndTime .lt. Me%StartTime) then
            write (*,*) 'Statistics End Time is BEFORE Statistics Start Time'
            write (*,*) 'Module :','HDF5Statistics'
            stop 'ReadGlobalData - ModuleHDF5Statistics - ERR30'
        endif

        if (Me%EndTime .eq. Me%StartTime) then
            write (*,*) 'Statistics End Time is EQUAL Statistics Start Time'
            write (*,*) 'Module :','HDF5Statistics'
            stop 'ReadGlobalData - ModuleHDF5Statistics - ERR40'
        endif
    end subroutine ReadGlobalData

    !--------------------------------------------------------------------------

    subroutine ReadParameters

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, ClientNumber
        type (T_Parameter),       pointer               :: NewParameter
        logical                                         :: BlockFound
        logical                                         :: AtLeastOneBlock = .false.

        !Begin-----------------------------------------------------------------

        ! Obtain Parameters for the statistics' calculation
do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                        ClientNumber    = ClientNumber,                 &
                                        block_begin     = '<BeginParameter>',           &
                                        block_end       = '<EndParameter>',             &
                                        BlockFound      = BlockFound,                   &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  

                    AtLeastOneBlock = .true.
                    
                    call AddParameter                    (NewParameter)

                    call ConstructParameters    (NewParameter)

                    nullify(NewParameter)

                else cd2
                    call Block_Unlock(Me%ObjEnterData,                                  & 
                                      ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ReadParameters - ModuleHDF5Statistics - ERR01'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ReadParameters - ModuleHDF5Statistics - ERR02'
            else cd1
                stop 'ReadParameters - ModuleHDF5Statistics - ERR03'
            end if cd1

        end do do1

        if (.not. AtLeastOneBlock) then                                            
            write(*,*) 'No property block is indicated in input file. '
            stop 'ReadParameters - ModuleHDF5Statistics - ERR04'
        end if

    end subroutine ReadParameters

    !--------------------------------------------------------------------------

    subroutine ReadHDF5FileName

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, ClientNumber, iflag, i
        type (T_HDF5File),       pointer                :: NewHDF5File
        logical                                         :: BlockFound
        logical                                         :: AtLeastOneBlock = .false.

        !Begin-----------------------------------------------------------------

        i = 0

        !Read input HDF5 files names
do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                        ClientNumber    = ClientNumber,                 &
                                        block_begin     = '<BeginHDF5File>',            &
                                        block_end       = '<EndHDF5File>',              &
                                        BlockFound      = BlockFound,                   &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                    
                    AtLeastOneBlock = .true.

                    call AddHDF5File                     (NewHDF5File)
                    
                    i = i + 1

                    call ConstructHDF5File               (NewHDF5File, i)

                    nullify(NewHDF5File)

                else cd2
                    call Block_Unlock(Me%ObjEnterData,                                  & 
                                      ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ReadHDF5FileName - ModuleHDF5Statistics - ERR01'

                    exit do1

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ReadHDF5FileName - ModuleHDF5Statistics - ERR02'
            else cd1
                stop 'ReadHDF5FileName - ModuleHDF5Statistics - ERR03'
            end if cd1

        end do do1

        if (.not. AtLeastOneBlock) then                                            
            write(*,*) 'No HDF5 file block is indicated in input file. '
            stop 'ReadHDF5FileName - ModuleHDF5Statistics - ERR04'
        end if

        !Read output HDF5 file name
        call GetData(Me%OutputFileName, Me%ObjEnterData, iflag,                         &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'HDF5Statistics',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
        stop 'ReadHDF5FileName - ModuleHDF5Statistics - ERR05'   

    end subroutine ReadHDF5FileName

    !--------------------------------------------------------------------------

    subroutine AddParameter (ObjParameter)

        !Arguments-------------------------------------------------------------
        type (T_Parameter),     pointer           :: ObjParameter

        !Local-----------------------------------------------------------------
        type (T_Parameter),     pointer           :: PreviousParameter
        type (T_Parameter),     pointer           :: NewParameter

        !Begin-----------------------------------------------------------------

        !Allocates new Parameter
        allocate (NewParameter)
        nullify  (NewParameter%Next)

        !Insert new Parameter into list and makes current ?? point to it
        if (.not. associated(Me%FirstParameter)) then
            Me%FirstParameter         => NewParameter
            ObjParameter              => NewParameter
            Me%ParameterNumber = 1
        else
            PreviousParameter         => Me%FirstParameter
            ObjParameter              => Me%FirstParameter%Next
            do while (associated(ObjParameter))
                PreviousParameter     => ObjParameter
                ObjParameter          => ObjParameter%Next
            enddo
            ObjParameter              => NewParameter
            PreviousParameter%Next    => NewParameter
            ! Count number of parameters in list
            Me%ParameterNumber = Me%ParameterNumber + 1
        end if

    end subroutine AddParameter

    !--------------------------------------------------------------------------

    subroutine ConstructParameters (NewParameter)

        !Arguments-------------------------------------------------------------
        type (T_Parameter),      pointer          :: NewParameter

        !External--------------------------------------------------------------
        integer                                   :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        
        ! Obtain parameter name
        call GetData(NewParameter%Name,                         &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'PROPERTY',                 &
                     ClientModule = 'HDF5Statistics',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                            &
        stop 'ConstructParameters - ModuleHDF5Statistics - ERR01'

        if (.not.CheckPropertyName(NewParameter%Name)) then
            write(*,*)
            write(*,*) 'The property name is not recognised by the model:'
            write(*,*) trim(NewParameter%Name)
            !stop 'ConstructParameters - ModuleHDF5Statistics - ERR02'
            !No stop is made because in some HDF5 files parameters are not
            !registred in Module GlobalData
        end if

        ! Obtain parameter group
        call GetData(NewParameter%Group,                        &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'HDF_GROUP',                &
                     default      = "/Results/"//trim(NewParameter%Name),&                     
                     ClientModule = 'HDF5Statistics',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                            &
        stop 'ConstructParameters - ModuleHDF5Statistics - ERR03'

    end subroutine ConstructParameters

    !--------------------------------------------------------------------------

    subroutine AddHDF5File(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type (T_HDF5File),     pointer           :: ObjHDF5File

        !Local-----------------------------------------------------------------
        type (T_HDF5File),     pointer           :: PreviousHDF5File
        type (T_HDF5File),     pointer           :: NewHDF5File

        !Begin-----------------------------------------------------------------

        !Allocates new HDF5File
        allocate (NewHDF5File)
        nullify  (NewHDF5File%Next)

        !Insert new Parameter into list and makes current ?? point to it
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

    subroutine ConstructHDF5File(NewHDF5File, i)

        !Arguments-------------------------------------------------------------
        type (T_HDF5File),      pointer           :: NewHDF5File
        integer                                   :: i

        !External--------------------------------------------------------------
        real,   dimension(6)                      :: Aux6
        integer                                   :: iflag, STAT_CALL, HDF5_READ, ObjHDF5
        logical                                   :: ReadWaterPointsName, Exist
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        !Instant from which the statistcs analysis will start for this hdf5 file
        call GetData(Aux6,                                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'START_TIME_FILE',                                  &
                     ClientModule = 'HDF5Statistics',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'ConstructHDF5File - ModuleHDF5Statistics - ERR10'
        
        if (iflag == 6) then
            call SetDate (NewHDF5File%StartTime, Aux6(1), Aux6(2),                      &
                                                 Aux6(3), Aux6(4),                      &
                                                 Aux6(5), Aux6(6))
            NewHDF5File%StartTimeDefault = .false.
        else
            NewHDF5File%StartTimeDefault = .true.
        endif

        !Instant from which the statistcs analysis will end for this hdf5 file
        call GetData(Aux6,                                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'END_TIME_FILE',                                    &
                     ClientModule = 'HDF5Statistics',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'ConstructHDF5File - ModuleHDF5Statistics - ERR20'
        
        if (iflag == 6) then
            call SetDate (NewHDF5File%EndTime, Aux6(1), Aux6(2),                        &
                                               Aux6(3), Aux6(4),                        &
                                               Aux6(5), Aux6(6))
            NewHDF5File%EndTimeDefault = .false.
        else
            NewHDF5File%EndTimeDefault = .true.
        endif
                
       
        ! Obtain HDF5 file name
        call GetData(NewHDF5File%Name,                                  &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'NAME',                             &
                     ClientModule = 'HDF5Statistics',                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                    &
        stop 'ConstructHDF5File - ModuleHDF5Statistics - ERR30'
        
        
        
        if (i==1) then

            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
            
            ObjHDF5 = 0

            !Open HDF5 file
            call ConstructHDF5 (ObjHDF5, trim(NewHDF5File%Name), HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDF5File - ModuleHDF5Statistics - ERR40'
            
            ReadWaterPointsName = .true.

            call GetHDF5DataSetExist (ObjHDF5, DataSetName ="/Grid/WaterPoints2D",      &
                                      Exist = Exist, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5File - ModuleHDF5Statistics - ERR50'

            if (Exist) then
                Me%Mapping%Name        = "WaterPoints2D"
                ReadWaterPointsName    = .false. 
                Me%File3D             = .false.
            endif      
            

            call GetHDF5DataSetExist (ObjHDF5, DataSetName ="/Grid/WaterPoints",      &
                                      Exist = Exist, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5File - ModuleHDF5Statistics - ERR60'

            if (Exist) then
                Me%Mapping%Name        = "WaterPoints"
                ReadWaterPointsName    = .false. 
                Me%File3D             = .false.
            endif                   
            
            call GetHDF5DataSetExist (ObjHDF5, DataSetName ="/Grid/BasinPoints",        &
                                      Exist = Exist, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5File - ModuleHDF5Statistics - ERR70'

            if (Exist) then
                Me%Mapping%Name        = "BasinPoints"
                ReadWaterPointsName    = .false. 
                Me%File3D             = .false.
            endif                  

            call GetHDF5DataSetExist (ObjHDF5, DataSetName ="/Grid/WaterPoints3D",      &
                                      Exist = Exist, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5File - ModuleHDF5Statistics - ERR80'

            if (Exist) then
                Me%Mapping%Name        = "WaterPoints3D"
                ReadWaterPointsName    = .false. 
                Me%File3D             = .true.                
            endif  
            
            if (ReadWaterPointsName) then

                ! Obtain HDF5 file's Map variable name
                call GetData(Me%Mapping%Name, Me%ObjEnterData, iflag,                   &
                             keyword      = 'HDF5_MAP_ITEM',                            &
                             SearchType   = FromFile,                                   &
                             ClientModule = 'HDF5Statistics',                           &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                              &
                stop 'ConstructHDF5File - ModuleHDF5Statistics - ERR90'            

                ! Obtain HDF5 file's dimension (3D or 2D)
                call GetData(Me%File3D, Me%ObjEnterData, iflag,                         &
                             keyword      = '3D_HDF5',                                  &
                             default      = .false.,                                    &
                             SearchType   = FromFile,                                   &
                             ClientModule = 'HDF5Statistics',                           &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                              &
                stop 'ConstructHDF5File - ModuleHDF5Statistics - ERR110'            

            endif            
    
            !Kill HDF5 file
            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDF5File - ModuleHDF5Statistics - ERR120'
            
        endif
        

    end subroutine ConstructHDF5File

    !--------------------------------------------------------------------------

    subroutine OpenAndDateHDF5Files

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type (T_HDF5File), pointer                  :: HDF5FileX
        character(len=StringLength)                 :: ParameterName
        logical                                     :: exist, FirstTime
        integer                                     :: HDF5_READ
        type (T_Parameter), pointer                 :: ObjParameter
        type(T_Time), dimension(:), pointer         :: AuxInstantsArray 
        integer                                     :: CurrentInstant, AuxNumberInstants
        logical                                     :: Relevant, GroupExist 
        real                                        :: LastDT, AuxDT
        real                                        :: Year, Month, Day, Hour 
        real                                        :: Minute, Second
      
        !Begin-----------------------------------------------------------------

        FirstTime = .true.

        HDF5FileX => Me%FirstHDF5File
        
        !In a DO cycle open all HDF5 files provided by the user
        do while (associated(HDF5FileX))

            !Verifies if file exists
            inquire(FILE = HDF5FileX%Name, EXIST = exist)
            if (.not. exist) then
                write(*,*)'HDF5 file does not exist:'//trim(HDF5FileX%Name)
                stop 'OpenAndDateHDF5Files - ModuleHDF5Statistics - ERR01'
            endif

            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

            !Open HDF5 file
            call ConstructHDF5 (HDF5FileX%HDFID, trim(HDF5FileX%Name),              &
                                HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'OpenAndDateHDF5Files - ModuleHDF5Statistics - ERR02'
            
            !Obtain start and end times of HDF5 file
            !(obtain number of instants) 
            call GetHDF5GroupNumberOfItems(HDF5FileX%HDFID, "/Time",                &
                                           HDF5FileX%NumberOfInstants, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            & 
            stop 'OpenAndDateHDF5Files - ModuleHDF5Statistics - ERR03'
            
            if (HDF5FileX%StartTimeDefault) then
                !(obtain HDF5 start time)
                HDF5FileX%StartTime = HDF5TimeInstant(1, HDF5FileX)
            endif

            if (HDF5FileX%EndTimeDefault) then
                !(obtain HDF5 end time)
                HDF5FileX%EndTime = HDF5TimeInstant(HDF5FileX%NumberOfInstants, HDF5FileX)
            endif

            !Get info about the rank and variables present
            !(only data needed for checking are obtained from each file)
            ObjParameter => Me%FirstParameter

            do while(associated(ObjParameter))

                call GetHDF5GroupExist (HDF5FileX%HDFID, ObjParameter%Group, GroupExist)

                !check if file contains parameter required
                if (.NOT. GroupExist) then  
                    write(*,*)'HDF5 file do not contain parameter required:'            &
                               //trim(HDF5FileX%Name)
                    write(*,*)'Parameter required:'//trim(ObjParameter%Name)
                    stop 'OpenAndDateHDF5Files - ModuleHDF5Statistics - ERR04'
                end if

                !get field ID, Rank
                call GetHDF5GroupID(HDF5FileX%HDFID, ObjParameter%Group,            &
                                1, ParameterName, ObjParameter%Units,               &
                                ObjParameter%Rank,                                  &
                                STAT = STAT_CALL)                                
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndDateHDF5Files - ModuleHDF5Statistics - ERR05'

                if (.not. FirstTime) then
                    !check if name of parameter is consistent with the previous
                    if (ParameterName .NE. ObjParameter%LastName) then
                        write(*,*)'HDF5 file do not contain parameter required:'    &
                                   //trim(HDF5FileX%Name)
                        write(*,*)'Parameter required:'//trim(ObjParameter%Name)
                        stop 'OpenAndDateHDF5Files - ModuleHDF5Statistics - ERR06'
                    end if

                    !check if rank of parameter is consistent with the previous
                    if (ObjParameter%Rank .NE. ObjParameter%LastRank) then
                        write(*,*)'File rank not consistent to previous rank:'      &
                                   //trim(HDF5FileX%Name)
                        write(*,*)'Parameter:'//trim(ObjParameter%Name)
                        stop 'OpenAndDateHDF5Files - ModuleHDF5Statistics - ERR07'
                    end if

                end if 

                ObjParameter%LastRank = ObjParameter%Rank
                ObjParameter%LastName = ParameterName               
                
                ObjParameter => ObjParameter%Next

            end do

            !Check if the HDF5 file is relevant for the statistics' calculation
            call HDF5Evaluator(HDF5FileX, Relevant)

            !If HDF5 file is relevant then obtain key parameters and instants
            if (Relevant) then

                !Get useful time information from file:
                !Set instant array for this file 
                allocate(AuxInstantsArray(1:HDF5FileX%NumberOfInstants))

                !Fill array with instants
                do CurrentInstant = 1, HDF5FileX%NumberOfInstants

                   AuxInstantsArray(CurrentInstant) = HDF5TimeInstant(CurrentInstant, & 
                                                                      HDF5FileX)

                end do

                !Get start and end instants for this file
                !select time window begin
                do CurrentInstant = 1, HDF5FileX%NumberOfInstants

                    if (AuxInstantsArray(CurrentInstant)                            &
                        .ge. HDF5FileX%StartFieldTime) then

                        HDF5FileX%StartInstant = CurrentInstant
                        HDF5FileX%StartFieldTime = HDF5TimeInstant(CurrentInstant,  & 
                                                                   HDF5FileX)
                
                        exit

                    end if

                end do
    
                !select time window end
                do CurrentInstant = HDF5FileX%StartInstant, HDF5FileX%NumberOfInstants

                    if (AuxInstantsArray(CurrentInstant)                            &
                        .eq. HDF5FileX%EndFieldTime) then

                        HDF5FileX%EndInstant = (CurrentInstant)
                        HDF5FileX%EndFieldTime = HDF5TimeInstant(CurrentInstant,    & 
                                                                 HDF5FileX)

                        exit

                    end if

                    if (AuxInstantsArray(CurrentInstant)                            &
                        .ge. HDF5FileX%EndFieldTime) then

                        HDF5FileX%EndInstant = (CurrentInstant-1)
                        HDF5FileX%EndFieldTime = HDF5TimeInstant(CurrentInstant-1,  & 
                                                                 HDF5FileX)

                        exit 

                    end if

                end do

                !Check to see the presence of only one time in the statistic
                if (HDF5FileX%StartFieldTime .eq. HDF5FileX%EndFieldTime) then
                    if (Me%StartTime >= HDF5FileX%StartTime .AND.                   &
                        Me%EndTime <= HDF5FileX%EndTime) then          
                        write(*,*) 'Statistic has only one time:' 
100                     format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)
                        call ExtractDate(HDF5FileX%StartFieldTime, Year, Month,     & 
                             Day, Hour, Minute, Second)
                        write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
                        write(*,*) 'This is not allowed.'
                        stop 'OpenAndDateHDF5Files - ModuleHDF5Statistics - ERR08'
                    end if
                end if

                AuxNumberInstants = HDF5FileX%NumberOfInstants

                HDF5FileX%NumberOfInstants = HDF5FileX%EndInstant -                 &
                                             HDF5FileX%StartInstant + 1

                allocate(HDF5FileX%InstantsArray(1:HDF5FileX%NumberOfInstants))

                HDF5FileX%InstantsArray = AuxInstantsArray(HDF5FileX%StartInstant:  &
                                                           HDF5FileX%EndInstant)

                !get instant times 
                call ObtainInstantsTimes(HDF5FileX)               

                !Calculate DT
                if (HDF5FileX%NumberOfInstants .ge. 2) then
                    Me%DT = HDF5FileX%InstantsArray(2) - HDF5FileX%InstantsArray(1)
                else 
                    write(*,*) 'HDF5 file with less than 2 time instants:'                
                    write(*,*) HDF5FileX%Name
                end if

                !Calculate DT for consistency checking
                if (AuxNumberInstants .ge. 3) then
                    AuxDT = AuxInstantsArray(3) - AuxInstantsArray(2)
                    Me%RegularDT = AuxDT
                else
                    !Check if this DT is equal to the DT of last file: they must equal!
                    if ((AuxDT .NE. LastDT) .and. (.not. FirstTime) .and.           &
                        (AuxNumberInstants .ge. 3)) then
                        write(*,*) 'HDF5 files do not have the same DT'                
                        write(*,*) Me%FirstStatHDF5File%Name
                        write(*,*) HDF5FileX%Name
                        stop 'OpenAndDateHDF5Files - ModuleHDF5Statistics - ERR08'                
                    end if 
                    LastDT = AuxDT
                    Me%RegularDT = Me%DT
                end if

                !Add file to list of relevant files
                call AddStatHDF5File(HDF5FileX)

                !next run of cycle (next parameter) is not the first one                
                FirstTime = .false.
           
                deallocate(AuxInstantsArray)
                nullify(AuxInstantsArray)

            end if

            call killhdf5(HDF5FileX%HDFID)           

            HDF5FileX => HDF5FileX%Next

        end do

        if (FirstTime) then

            !No HDF5 file is provided with suitable data
            call ExtractDate(Me%StartTime, Year, Month, Day, Hour,                  &  
                             Minute, Second)        
            write(*,*) 'Data lacking from'      
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
            write(*,*) 'to'      
            call ExtractDate(Me%EndTime, Year, Month, Day, Hour, Minute, Second)
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
            write(*,*) 'There are not requested data in the HDF5 files provided.'      
            stop 'OpenAndDateHDF5Files - ModuleHDF5Statistics - ERR09'                

        else

            !assume VariableDT
            Me%VariableDT = .True.

            !Check for only two statistic times, from same file (RegularDT must equal DT):
            if (Me%FirstStatHDF5File%Name == Me%LastStatHDF5File%Name) then
                if (Me%FirstStatHDF5File%EndFieldTime .le.                          &
                    (Me%FirstStatHDF5File%StartFieldTime + Me%RegularDT)) then
                    Me%RegularDT = Me%DT
                endif
            endif

            !For each HDF5 file needed for the statistics' calculation
            HDF5FileX => Me%FirstStatHDF5File

            do while(associated(HDF5FileX))

                !Get instants' times and put them in list 
                call AddStatInstantsTimes(HDF5FileX)

                HDF5FileX => HDF5FileX%Next            

            end do
           
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
            write(*,*)'Same period is contained in two HDF5 files:'
            write(*,*) HDF5FileX%Name
            write(*,*) OtherFile%Name
            stop 'HDF5Evaluator - ModuleHDF5Statistics - ERR01'
        end if

        !See if the file is to be considered for statistics 
        !Check start and end time
        if ((HDF5FileX%EndTime >= Me%StartTime) .and.                     &
            (HDF5FileX%EndTime <= Me%EndTime)) then

            !End statistics time is between start and end times of file
            if (HDF5FileX%StartTime < Me%StartTime) then

                !Start statistics time is after start time of file
                HDF5FileX%StartFieldTime = Me%StartTime
                HDF5FileX%EndFieldTime   = HDF5FileX%EndTime

                Relevant = .TRUE.

            else 

                !Start statistics time is before start time of file
                HDF5FileX%StartFieldTime = HDF5FileX%StartTime
                HDF5FileX%EndFieldTime   = HDF5FileX%EndTime

                Relevant = .TRUE.

            end if

        else if ((HDF5FileX%StartTime >= Me%StartTime) .and.              &
                 (HDF5FileX%StartTime <= Me%EndTime)) then

            !End statistics time is before end time of file
            HDF5FileX%StartFieldTime = HDF5FileX%StartTime
            HDF5FileX%EndFieldTime   = Me%EndTime

            Relevant = .TRUE.

        else if ((HDF5FileX%StartTime < Me%StartTime) .and.               &
                 (HDF5FileX%EndTime > Me%EndTime)) then

            !Statistics period is contained in file
            HDF5FileX%StartFieldTime = Me%StartTime
            HDF5FileX%EndFieldTime   = Me%EndTime

            Relevant = .TRUE.

        end if

    end subroutine HDF5Evaluator
 
    !--------------------------------------------------------------------------

    type(T_Time) function HDF5TimeInstant(Instant, ObjHDF5File)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        type(T_HDF5File), pointer               :: ObjHDF5File
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer             :: TimeVector

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (ObjHDF5File%HDFID, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadWindow   (HDF5ID         = ObjHDF5File%HDFID,        &
                             GroupName      = "/Time",                  &
                             Name           = "Time",                   &
                             Array1D        = TimeVector,               &
                             OutputNumber   = Instant,                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                      &
        stop 'HDF5TimeInstant - ModuleHDF5Statistics - ERR01'

        call SetDate(HDF5TimeInstant, Year  = TimeVector(1),            &
                     Month  = TimeVector(2), Day      = TimeVector(3),  &
                     Hour   = TimeVector(4), Minute   = TimeVector(5),  &
                     Second = TimeVector(6))


        deallocate(TimeVector)
        nullify   (TimeVector)

    end function HDF5TimeInstant

    !--------------------------------------------------------------------------

    subroutine OpenOutputFiles

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE
        real,    dimension(:), pointer              :: TimePtr
        real,    dimension(6), target               :: AuxTime

        !Begin-----------------------------------------------------------------

        ! Get grid values for HDF5
        call ConstructHDF5Grid

        !Create HDF5 file
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjStatHDF5, trim(Me%OutputFileName),               &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR01'

        !Set grid
        call HDF5SetLimits(Me%ObjStatHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB+1,          &
                           Me%WorkSize%JLB, Me%WorkSize%JUB+1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR02'
        
        if (Me%ConnectionX%Exist) then        

            call HDF5WriteData   (Me%ObjStatHDF5, "/Grid", Me%ConnectionX%Name,         &
                                  Me%ConnectionX%Units,                                 &
                                  Array2D = Me%ConnectionX%RealValues2D,                &
                                  STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR03'

        endif        
        
        if (Me%ConnectionY%Exist) then        

            call HDF5WriteData   (Me%ObjStatHDF5, "/Grid", Me%ConnectionY%Name,             &
                                  Me%ConnectionY%Units,                                     &
                                  Array2D = Me%ConnectionY%RealValues2D,                    &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR04'
        endif

        call HDF5WriteData   (Me%ObjStatHDF5, "/Grid", Me%Latitude%Name,                &
                              Me%Latitude%Units,                                        &
                              Array2D = Me%Latitude%RealValues2D,                       &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR05'

        call HDF5WriteData   (Me%ObjStatHDF5, "/Grid", Me%Longitude%Name,               &
                              Me%Longitude%Units,                                       &
                              Array2D = Me%Longitude%RealValues2D,                      &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR06'

        !Set bathymetry          
        call HDF5SetLimits(Me%ObjStatHDF5, Me%WorkSize%ILB,Me%WorkSize%IUB,             &
                   Me%WorkSize%JLB,Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR07'

        call HDF5WriteData(Me%ObjStatHDF5,                                              &
                   "/Grid",                                                             &
                   Me%Bathymetry%Name, Me%Bathymetry%Units,                             &
                   Array2D      = Me%Bathymetry%RealValues2D,                           &
                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR08'

        !Set map
        if (Me%File3D) then 
            call HDF5SetLimits(Me%ObjStatHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,        &
                               Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                               Me%WorkSize%KLB, Me%WorkSize%KUB,                        &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &      
            stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR09'

            call HDF5WriteData(Me%ObjStatHDF5, "/Grid", trim(Me%Mapping%Name),          &
                               trim(Me%Mapping%Units),                                  & 
                               Array3D      = Me%Mapping%IntegerValues3D,               &
                               STAT         = STAT_CALL)      
            if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR10'
            
            if (Me%ExistVerticalZ) then
                !added set limits to represent all vertical layers (needed the Me%WorkSize%KLB-1)
                call HDF5SetLimits(Me%ObjStatHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,        &
                                   Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                                   Me%WorkSize%KLB-1, Me%WorkSize%KUB,                      &
                                   STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &      
                stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR11'
                
                call HDF5WriteData(HDF5ID       = Me%ObjStatHDF5,                               &
                                   GroupName    = "/Grid/VerticalZ",                            &
                                   Name         = "Vertical",                                   &
                                   units        = "m",                                          & 
                                   Array3D      = Me%VerticalZ,                                 &
                                   OutputNumber = 1,                                            &
                                   STAT         = STAT_CALL)      
                if (STAT_CALL /= SUCCESS_)                                                      &
                    stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR12'            
                
            endif
            

            if (Me%AditionalMap) then

                call HDF5SetLimits(Me%ObjStatHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,    &
                                   Me%WorkSize%JLB, Me%WorkSize%JUB,                    &
                                   STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &      
                stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR13'

                call HDF5WriteData(Me%ObjStatHDF5, "/Grid",                             &
                                   trim(Me%Mapping%AditionalName),                      &
                                   trim(Me%Mapping%Units),                              & 
                                   Array2D      = Me%Mapping%IntegerValues2D,           &
                                   STAT         = STAT_CALL)      
                if (STAT_CALL /= SUCCESS_)                                              &
                stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR14'

            endif 
            
        else
            call HDF5WriteData(Me%ObjStatHDF5, "/Grid", trim(Me%Mapping%Name),          &
                               trim(Me%Mapping%Units),                                  & 
                               Array2D      = Me%Mapping%IntegerValues2D,               &
                               STAT         = STAT_CALL)    
            if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR15'
        endif

        
        !Write First and Last Statistics Time to know the begin and end of statistics period

        call ExtractDate (Me%FirstStatHDF5File%StartFieldTime,                          &
                          AuxTime(1), AuxTime(2), AuxTime(3),                           &
                          AuxTime(4), AuxTime(5), AuxTime(6))


        TimePtr => AuxTime

        call HDF5SetLimits  (Me%ObjStatHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR16'

        call HDF5WriteData  (Me%ObjStatHDF5, "/Time",                                   &
                 "Time", "YYYY/MM/DD HH:MM:SS",                                         &
                 Array1D = TimePtr,                                                     &
                 OutputNumber = 1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR17'


        call HDF5FlushMemory(Me%ObjStatHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
        stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR18'      

        
        call ExtractDate   (Me%LastStatHDF5File%EndFieldTime,                                    &
                           AuxTime(1), AuxTime(2), AuxTime(3),                         &
                           AuxTime(4), AuxTime(5), AuxTime(6))

 
        TimePtr => AuxTime

        call HDF5SetLimits  (Me%ObjStatHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR19'

        call HDF5WriteData  (Me%ObjStatHDF5, "/Time",                                   &
                 "Time", "YYYY/MM/DD HH:MM:SS",                                         &
                 Array1D = TimePtr,                                                     &
                 OutputNumber = 2, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR20'
        ! ----- \\\\ //// ------

        call HDF5FlushMemory(Me%ObjStatHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
        stop 'OpenOutputFiles - ModuleHDF5Statistics - ERR21'      


        call KillGridFields

    end subroutine OpenOutputFiles

    !--------------------------------------------------------------------------

    subroutine ObtainInstantsTimes(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                 :: ObjHDF5File
        
        !Local-----------------------------------------------------------------
        type (T_StatisticsTime), pointer           :: NewTime
        type (T_StatisticsTime), pointer           :: ObjTime
        type (T_StatisticsTime), pointer           :: PreviousTime
        integer                                   :: CurrentInstant
        
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

    subroutine AddStatHDF5File(HDF5FileX)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX

        !Local-----------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileAux
        type(T_HDF5File), pointer                   :: PreviousHDF5File, LastHDF5File
      
        !Begin-----------------------------------------------------------------

        if (.not. associated(Me%FirstStatHDF5File)) then

            call CreateStatHDF5File(Me%FirstStatHDF5File, HDF5FileX)
            call CreateStatHDF5File(Me%LastStatHDF5File, HDF5FileX)
            deallocate(Me%FirstStatHDF5File%Next) 
            nullify(Me%FirstStatHDF5File%Next)

        else

            if (HDF5FileX%StartFieldTime < Me%FirstStatHDF5File%StartFieldTime) then
                !current file should be the first file in list

                !save the previous list 
                allocate(HDF5FileAux)
                call CreateStatHDF5File(HDF5FileAux, Me%FirstStatHDF5File)
                HDF5FileAux%Next => Me%FirstStatHDF5File%Next              

                !make the first element in the list of relevant files equal to current file
                call CreateStatHDF5File(Me%FirstStatHDF5File, HDF5FileX)
                Me%FirstStatHDF5File%Next => HDF5FileAux

            else
                !check next files in list

                !locate previous file in the first file
                allocate(PreviousHDF5File)
                PreviousHDF5File => Me%FirstStatHDF5File                   

                do while(associated(PreviousHDF5File))

                    if (.not. associated(PreviousHDF5File%Next)) then
        
                        !current file is the last file in the list of relevant files
                        call CreateStatHDF5File(PreviousHDF5File%Next, HDF5FileX)
                        allocate(LastHDF5File)
                        LastHDF5File => PreviousHDF5File%Next
                        deallocate(LastHDF5File%Next)
                        nullify(LastHDF5File%Next)
                        call CreateStatHDF5File(Me%LastStatHDF5File, HDF5FileX)

                        !current file was added to list
                        exit

                    else

                        !check if current file should be located before the next file
                        if (HDF5FileX%StartFieldTime < PreviousHDF5File%Next%StartFieldTime) then
                            !current file should be located before next file

                            !save the previous list begining in PreviousHDF5File%Next
                            allocate(LastHDF5File)
                            LastHDF5File => PreviousHDF5File%Next
                            allocate(HDF5FileAux)

                            call CreateStatHDF5File(HDF5FileAux, LastHDF5File)
                            HDF5FileAux%Next => LastHDF5File%Next

                            call CreateStatHDF5File(LastHDF5File, HDF5FileX)

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

    end subroutine AddStatHDF5File

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine SamePeriodHDF5(HDF5FileX, HDF5FileAux, STAT)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX
        integer, intent(OUT)                        :: STAT
        type(T_HDF5File), pointer                   :: HDF5FileAux

        !Local-----------------------------------------------------------------
      
        !Begin-----------------------------------------------------------------

        STAT = UNKNOWN_

        HDF5FileAux => Me%FirstHDF5File

        do while (associated(HDF5FileAux))

            if (HDF5FileX%Name .NE. HDF5FileAux%Name) then 
            !(not the same HDF5 file)

                !Check if the same period is in more than one file
                if (((HDF5FileX%StartTime >= HDF5FileAux%StartTime)             &
                    .and. (HDF5FileX%StartTime < HDF5FileAux%EndTime))          &
                    .or. ((HDF5FileX%EndTime > HDF5FileAux%StartTime)           &
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

    subroutine CreateStatHDF5File(HDF5FileNew, HDF5FileX)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX
        type(T_HDF5File), pointer                   :: HDF5FileNew

        !Local-----------------------------------------------------------------
      
        !Begin-----------------------------------------------------------------

        !This subroutine atributes the values of the fields of a HDF5File to another HDF5File
        allocate(HDF5FileNew)

        HDF5FileNew%Name             =  HDF5FileX%Name
        HDF5FileNew%StartTime        =  HDF5FileX%StartTime
        HDF5FileNew%EndTime          =  HDF5FileX%EndTime
        HDF5FileNew%StartFieldTime   =  HDF5FileX%StartFieldTime
        HDF5FileNew%EndFieldTime     =  HDF5FileX%EndFieldTime
        HDF5FileNew%Rank             =  HDF5FileX%Rank
        HDF5FileNew%NumberOfInstants =  HDF5FileX%NumberOfInstants
        HDF5FileNew%StartInstant     =  HDF5FileX%StartInstant
        HDF5FileNew%EndInstant       =  HDF5FileX%EndInstant
        HDF5FileNew%InstantsArray    => HDF5FileX%InstantsArray
        HDF5FileNew%FirstInstantTime => HDF5FileX%FirstInstantTime

    end subroutine CreateStatHDF5File

    !--------------------------------------------------------------------------

    subroutine AddStatInstantsTimes(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                 :: ObjHDF5File
        
        !Local-----------------------------------------------------------------
        type (T_StatisticsTime), pointer           :: NewTime
        type (T_StatisticsTime), pointer           :: ObjTime
        type (T_StatisticsTime), pointer           :: PreviousTime
        type (T_StatisticsTime), pointer           :: ObjInstantTime
        
        !Begin-----------------------------------------------------------------

        ObjInstantTime => ObjHDF5File%FirstInstantTime
        
        do while (associated(ObjInstantTime))

            !Allocates new instance
            allocate (NewTime)
            nullify  (NewTime%Next)

            NewTime%Time = ObjInstantTime%Time

            !Insert New Instance into list and makes Current point to it
            if (.not. associated(Me%FirstStatisticsTime)) then
            !FirstField should be the same for all HDF5 files 
                Me%FirstStatisticsTime   => NewTime
                ObjTime                 => NewTime
            else
                PreviousTime            => Me%FirstStatisticsTime
                ObjTime                 => Me%FirstStatisticsTime%Next
                do while (associated(ObjTime))
                    PreviousTime        => ObjTime
                    ObjTime             => ObjTime%Next
                enddo
                ObjTime                 => NewTime
                PreviousTime%Next       => NewTime
            endif

            ObjInstantTime => ObjInstantTime%Next          

        end do 

    end subroutine AddStatInstantsTimes

    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Grid

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_HDF5File), pointer                  :: HDF5FileX
        integer                                     :: Rank
        integer, dimension(7)                       :: Dimensions
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_READ
        integer                                     :: GridVarNumber, n
        character(len=StringLength)                 :: GridVariableName
        character(len=StringLength)                 :: GridVariableUnits
        type (T_Parameter), pointer                 :: ObjParameter
        !logical                                     :: ExistVerticalZ
        !real,   dimension(:,:,:), pointer           :: VerticalZ
        integer                                     :: i, j, k

        !Begin-----------------------------------------------------------------

        !Extract from first relevant file: grid, bathymetry, mapping
        !(this assuming that these are equal in all files)
        HDF5FileX => Me%FirstStatHDF5File

        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Open HDF5 file
        call ConstructHDF5 (HDF5FileX%HDFID, trim(HDF5FileX%Name),                      &
                            HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR10'

        Me%ConnectionX%Name = trim("ConnectionX")
        Me%ConnectionY%Name = trim("ConnectionY")
        Me%Latitude%Name = trim("Latitude")
        Me%Longitude%Name = trim("Longitude")
        Me%Bathymetry%Name = trim("Bathymetry")

        !Get position of grid variables to acess later
        call GetHDF5GroupNumberOfItems(HDF5FileX%HDFID, "/Grid", GridVarNumber,         & 
                                       STAT = STAT_CALL)
        do n = 1, GridVarNumber
            call GetHDF5GroupID(HDF5FileX%HDFID, trim("/Grid"),                         &
                            n, GridVariableName,                                        &
                            GridVariableUnits,                                          &
                            STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_)                                                &  
            stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR20'
            if (GridVariableName == trim(Me%Bathymetry%Name)) then
                Me%Bathymetry%Position = n
                Me%Bathymetry%Units = GridVariableUnits
            elseif (GridVariableName == trim(Me%ConnectionX%Name)) then
                Me%ConnectionX%Position = n
                Me%ConnectionX%Units = GridVariableUnits
            elseif (GridVariableName == trim(Me%ConnectionY%Name)) then
                Me%ConnectionY%Position = n
                Me%ConnectionY%Units = GridVariableUnits
            elseif (GridVariableName == trim(Me%Latitude%Name)) then
                Me%Latitude%Position = n
                Me%Latitude%Units = GridVariableUnits
            elseif (GridVariableName == trim(Me%Longitude%Name)) then
                Me%Longitude%Position = n
                Me%Longitude%Units = GridVariableUnits
            elseif (GridVariableName == trim(Me%Mapping%Name)) then
                Me%Mapping%Position = n
                Me%Mapping%Units = GridVariableUnits
           endif
        end do

        !Get grid:
        !Mapping dimensions and units
        call GetHDF5GroupID(HDF5FileX%HDFID, trim("/Grid"),                             &
                            Me%Mapping%Position, Me%Mapping%Name,                       &
                            Me%Mapping%Units, Rank,                                     &
                            Dimensions,                                                 &
                            STAT = STAT_CALL)                                
        if (STAT_CALL .NE. SUCCESS_)                                                    &  
            stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR30'

        !Get file dimensions
        !if (Me%SubSetanalysisON) then
        !    Me%WorkSize%ILB = ILBout
        !    Me%WorkSize%IUB = IUBout
        !    Me%WorkSize%JLB = JLBout
        !    Me%WorkSize%JUB = JUBout
        !
        !else                    
            Me%WorkSize%ILB = 1  
            Me%WorkSize%IUB = Dimensions(1)
            Me%WorkSize%JLB = 1
            Me%WorkSize%JUB = Dimensions(2)
        !endif            
            
        if (Me%File3D) then           
            Me%WorkSize%KLB = 1
        else 
            Me%WorkSize%KLB = 0
        endif
        Me%WorkSize%KUB = Dimensions(3)
        
        Me%Size%ILB = Me%WorkSize%ILB - 1
        Me%Size%JLB = Me%WorkSize%JLB - 1
        Me%Size%IUB = Me%WorkSize%IUB + 1 
        Me%Size%JUB = Me%WorkSize%JUB + 1
        Me%Size%KLB = Me%WorkSize%KLB - 1
        Me%Size%KUB = Me%WorkSize%KUB + 1

        !Allocate variable data matrixes
        nullify(Me%ConnectionX%RealValues2D, Me%ConnectionY%RealValues2D)
        nullify(Me%Longitude%RealValues2D, Me%Latitude%RealValues2D,                    & 
                Me%Bathymetry%RealValues2D)
        
        allocate(Me%ConnectionX%RealValues2D(Me%Size%ILB:Me%Size%IUB,                   &
                 Me%Size%JLB:Me%Size%JUB))
        allocate(Me%ConnectionY%RealValues2D(Me%Size%ILB:Me%Size%IUB,                   & 
                 Me%Size%JLB:Me%Size%JUB))
        allocate(Me%Longitude%RealValues2D(Me%Size%ILB:Me%Size%IUB,                     & 
                 Me%Size%JLB:Me%Size%JUB))
        allocate(Me%Latitude%RealValues2D(Me%Size%ILB:Me%Size%IUB,                      & 
                 Me%Size%JLB:Me%Size%JUB))
        allocate(Me%Bathymetry%RealValues2D(Me%Size%ILB:Me%Size%IUB,                    &   
                 Me%Size%JLB:Me%Size%JUB))
                 !(allocate always with size, but read and write with adequate!)
        
        if (Me%File3D) then 
            nullify(Me%Mapping%IntegerValues3D) 
            allocate(Me%Mapping%IntegerValues3D(Me%Size%ILB:Me%Size%IUB,                &
                     Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))
        else 
            nullify(Me%Mapping%IntegerValues2D) !, MapPoints2D)
            allocate(Me%Mapping%IntegerValues2D(Me%Size%ILB:Me%Size%IUB,                &
                     Me%Size%JLB:Me%Size%JUB))
                     !(allocate always with size!)
        endif

        !Read grid values
        !connections and coordinates
        call HDF5SetLimits (HDF5FileX%HDFID, Me%WorkSize%ILB,                           &
                            Me%WorkSize%IUB+1, Me%WorkSize%JLB,                         &
                            Me%WorkSize%JUB+1, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    & 
            stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR40'
        
        call GetHDF5DataSetExist (HDF5FileX%HDFID,                                      &
                                  DataSetName ="/Grid/"//trim(Me%ConnectionX%Name),     &
                                  Exist = Me%ConnectionX%Exist, STAT= STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR50'
        endif            
        
        if (Me%ConnectionX%Exist) then

            call HDF5ReadWindow(HDF5FileX%HDFID, "/Grid",                                 &
                              trim(Me%ConnectionX%Name),                                &
                              Array2D      = Me%ConnectionX%RealValues2D,               &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR60'
            
        endif   
        
        call GetHDF5DataSetExist (HDF5FileX%HDFID,                                      &
                                  DataSetName ="/Grid/"//trim(Me%ConnectionY%Name),     &
                                  Exist = Me%ConnectionY%Exist, STAT= STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR70'
        endif            
                
        
        if (Me%ConnectionY%Exist) then                 

            call HDF5ReadWindow(HDF5FileX%HDFID, "/Grid",                                 &
                              trim(Me%ConnectionY%Name),                                &
                              Array2D      = Me%ConnectionY%RealValues2D,               &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR80'

        endif

        call HDF5ReadWindow(HDF5FileX%HDFID, "/Grid",                                     &
                          trim(Me%Latitude%Name),                                       &
                          Array2D      = Me%Latitude%RealValues2D,                      &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR90'

        call HDF5ReadWindow(HDF5FileX%HDFID, "/Grid",                                     &
                          trim(Me%Longitude%Name),                                      &
                          Array2D      = Me%Longitude%RealValues2D,                     &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR100'

        !bathymetry
        call HDF5SetLimits (HDF5FileX%HDFID, Me%WorkSize%ILB,                           &
                            Me%WorkSize%IUB, Me%WorkSize%JLB,Me%WorkSize%JUB,           &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    & 
        stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR110'

        call HDF5ReadWindow(HDF5FileX%HDFID, "/Grid",                                     &
                          trim(Me%Bathymetry%Name),                                     &
                          Array2D      = Me%Bathymetry%RealValues2D,                    &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR120'

        !mapping
        if (Me%File3D) then 
        
            call HDF5SetLimits (HDF5FileX%HDFID, Me%WorkSize%ILB,                       &
                                Me%WorkSize%IUB, Me%WorkSize%JLB,Me%WorkSize%JUB,       &
                                Me%WorkSize%KLB,Me%WorkSize%KUB,                        &
                                STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                & 
            stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR130'
            
            call HDF5ReadWindow(HDF5FileX%HDFID, "/Grid",                                 &
                              trim(Me%Mapping%Name),                                    &
                              Array3D      = Me%Mapping%IntegerValues3D,                &
                              STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR140'        
        
            call GetHDF5DataSetExist (HDF5FileX%HDFID,                                  &
                                  DataSetName ="/Grid/VerticalZ/Vertical_00001",        &
                                  Exist = Me%ExistVerticalZ, STAT= STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)  then
                stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR150'
            endif                                       
                                  
            if (Me%ExistVerticalZ) then
            
                allocate(Me%DZ3D  (Me%Size%ILB:Me%Size%IUB,                             &
                                   Me%Size%JLB:Me%Size%JUB,                             &
                                   Me%Size%KLB:Me%Size%KUB))
                allocate(Me%VerticalZ(Me%Size%ILB:Me%Size%IUB,                          &
                                      Me%Size%JLB:Me%Size%JUB,                          &
                                      Me%Size%KLB:Me%Size%KUB))
                    
                Me%VerticalZ(:,:,:) = FillValueReal 
                Me%DZ3D  (:,:,:) = FillValueReal
                
                call HDF5SetLimits(HDF5FileX%HDFID, Me%WorkSize%ILB, Me%WorkSize%IUB,   &
                                    Me%WorkSize%JLB, Me%WorkSize%JUB,                   &
                                    Me%WorkSize%KLB-1, Me%WorkSize%KUB,                 &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)  then
                    stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR160'
                endif                    

                
                call HDF5ReadWindow(HDF5FileX%HDFID, "/Grid/VerticalZ",                   &
                                  "Vertical_00001",                                     &
                                  Array3D      = Me%VerticalZ,                             &
                                  OffSet3       = 0,                                  &
                                  STAT = STAT_CALL) 

                if (STAT_CALL /= SUCCESS_)  then
                    stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR170'
                endif       
                                  
                do k =  Me%WorkSize%KLB, Me%WorkSize%KUB
                do j =  Me%WorkSize%JLB, Me%WorkSize%JUB
                do i =  Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%Mapping%IntegerValues3D(i, j, k) == WaterPoint) then
                        Me%DZ3D  (i,j,k)  = Me%VerticalZ(i,j,k-1) - Me%VerticalZ(i,j,k)
                    endif
                enddo
                enddo
                enddo
                
                !deallocate(Me%VerticalZ)
                    
            endif
        
        
            call HDF5SetLimits (HDF5FileX%HDFID, Me%WorkSize%ILB,                       &
                                Me%WorkSize%IUB, Me%WorkSize%JLB,Me%WorkSize%JUB,       &
                                Me%WorkSize%KLB,Me%WorkSize%KUB,                        &
                                STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) then
                stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR180'
            endif                
            

            !For each parameter calculate statistic 
            ObjParameter => Me%FirstParameter

            !Check if there are 2D parameters, if so a 2D mapping must be created
do2 :       do while(associated(ObjParameter))

                if (ObjParameter%Rank == 2) then

                    Me%Mapping%AditionalName = trim("MappingPoints2D")

                    nullify(Me%Mapping%IntegerValues2D)
                    allocate(Me%Mapping%IntegerValues2D(Me%Size%ILB:Me%Size%IUB,        &
                             Me%Size%JLB:Me%Size%JUB))

                    Me%Mapping%IntegerValues2D =                                        & 
                                Me%Mapping%IntegerValues3D(Me%WorkSize%ILB:             &
                                Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB,        &
                                Me%Size%KUB)
                    !(assume that relevant mapping is for the upper layer)

                    Me%AditionalMap = .true.

                    write(*,*) 'Aditional grid variable MappingPoints2D created'
                    write(*,*) 'based on upper 3D layer mapping.'

                    exit do2
        
                endif

                ObjParameter => ObjParameter%Next               

            end do do2

        else 

            call HDF5ReadWindow(HDF5FileX%HDFID, "/Grid",                             &
                              trim(Me%Mapping%Name),                                &
                              Array2D      = Me%Mapping%IntegerValues2D,            &
                              STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ConstructHDF5Grid - ModuleHDF5Statistics - ERR190'
            
        endif

        call killhdf5(HDF5FileX%HDFID)           

    end subroutine ConstructHDF5Grid

    !--------------------------------------------------------------------------

    subroutine ConstructStatisticGroup

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        character(StringLength)                     :: AuxSYear, AuxSMonth, AuxSDay
        character(StringLength)                     :: AuxSHour, AuxSMin, AuxSSec
        character(StringLength)                     :: AuxEYear, AuxEMonth, AuxEDay
        character(StringLength)                     :: AuxEHour, AuxEMin, AuxESec
        real                                        :: Year, Month, Day, Hour 
        real                                        :: Minute, Second
        character(StringLength)                     :: AuxName1, AuxName2, AuxName3
        character(StringLength)                     :: AuxName4
        
        !Begin-----------------------------------------------------------------

        !Get start and end statistic dates in pieces
        !(these are the actual start and end dates, not the ones specified by the user)
        call ExtractDate(Me%FirstStatHDF5File%StartFieldTime, Year, Month,          & 
                         Day, Hour, Minute, Second)

        !convert date itens to chars
        write(AuxSYear, fmt=*) int(Year)
        write(AuxSMonth, fmt=*) int(Month)
        write(AuxSDay, fmt=*) int(Day)
        write(AuxSHour, fmt=*) int(Hour)
        write(AuxSMin, fmt=*) int(Minute)
        write(AuxSSec, fmt=*) int(Second)

        !Adjust to put zeros when needed
        if (Month < 10) AuxSMonth = "0"//trim(adjustl(AuxSMonth))
        if (Day <10) AuxSDay = "0"//trim(adjustl(AuxSDay))
        if (Hour <10) AuxSHour = "0"//trim(adjustl(AuxSHour))
        if (Minute <10) AuxSMin = "0"//trim(adjustl(AuxSMin))
        if (Second <10) AuxSSec = "0"//trim(adjustl(AuxSSec))     

        call ExtractDate(Me%LastStatHDF5File%EndFieldTime, Year, Month,             & 
                         Day, Hour, Minute, Second)

        !convert date itens to chars
        write(AuxEYear, fmt=*) int(Year)
        write(AuxEMonth, fmt=*) int(Month)
        write(AuxEDay, fmt=*) int(Day)
        write(AuxEHour, fmt=*) int(Hour)
        write(AuxEMin, fmt=*) int(Minute)
        write(AuxESec, fmt=*) int(Second)

        !Adjust to put zeros when needed
        if (Month < 10) AuxEMonth = "0"//trim(adjustl(AuxEMonth))
        if (Day <10) AuxEDay = "0"//trim(adjustl(AuxEDay))
        if (Hour <10) AuxEHour = "0"//trim(adjustl(AuxEHour))
        if (Minute <10) AuxEMin = "0"//trim(adjustl(AuxEMin))
        if (Second <10) AuxESec = "0"//trim(adjustl(AuxESec))     


        AuxName1 = trim(adjustl(AuxSYear))//trim(adjustl(AuxSMonth))//              &
                   trim(adjustl(AuxSDay))

        AuxName2 = "_"//trim(adjustl(AuxSHour))//trim(adjustl(AuxSMin))//           &
                   trim(adjustl(AuxSSec))

        AuxName3 = "-"//trim(adjustl(AuxEYear))//trim(adjustl(AuxEMonth))//         &
                    trim(adjustl(AuxEDay))

        AuxName4 = "_"//trim(adjustl(AuxEHour))//trim(adjustl(AuxEMin))//           &
                   trim(adjustl(AuxESec))    

        Me%StatisticGroupName = trim(AuxName1)//trim(AuxName2)//trim(AuxName3)      &
                                //trim(AuxName4)//"/"                                                         

    end subroutine ConstructStatisticGroup

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyHDF5Statistics

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_HDF5File), pointer                   :: ObjHDF5File
        type (T_Parameter), pointer                 :: ObjParameter, ParameterToKill
        logical                                     :: FirstFile
        type(T_Time)                                :: LastStartTime
        logical                                     :: Running
        real                                        :: LastDT
        integer                                     :: Count = 0
        integer                                     :: STAT_CALL
        type (T_StatisticsTime), pointer            :: ObjStatisticsTime
        real                                        :: DT

        !----------------------------------------------------------------------

        ObjStatisticsTime => Me%FirstStatisticsTime

        !Cycle each relevant HDF5 file
        FirstFile =.true.

        !For each HDF5 file needed for the time serie
        ObjHDF5File => Me%FirstStatHDF5File
        
        !Initialize variables needed to check lacks in data
        LastStartTime = Me%StartTime
        LastDT = 0

        do while(associated(ObjHDF5File))

            !Open and read relevant data
            call OpenAndReadHDF5File(FirstFile, ObjHDF5File, LastDT, LastStartTime)
            !(HDF5 parameter data must be read in this cycle for memory to be released
            !after writing statistics)

            Me%CurrentTime  = ObjStatisticsTime%Time

            if (FirstFile) then
                FirstFile = .false. !next file is not first
                
                call ActualizeCurrentTime(Me%ObjTime, Me%RegularDT, STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)                                             &
                stop 'ModifyHDF5Statistics - ModuleHDF5Statistics - ERR01'
            endif

            !Cycle each parameter
            Running      = .true.
     
            write(*,*)'Calculating statistics...'

            do while (Running)

                !For each parameter calculate statistic 
                ObjParameter => Me%FirstParameter
            
                Count = Count + 1 
            
                if (Count .EQ. 1) then 

                    ObjParameter%CurrentField => ObjParameter%FirstField 

                end if

                do while(associated(ObjParameter))
            

                    select case (ObjParameter%Rank)

                        !Values in the begining of files in middle of time serie are not 
                        !writen. These values are very close in time to last value of 
                        !previous file.

                        case(2)

                            call CalculateHDF5Statistics2D(ObjParameter%            &
                                                           CurrentField%Values2D,   & 
                                                           ObjParameter%Statistics%ID)
 
                        case(3)

                            call CalculateHDF5Statistics3D(ObjParameter%            &
                                                           CurrentField%Values3D,   & 
                                                           ObjParameter%Statistics%ID)

                    case default 
            
                        write(*,*)'Statistics created only for 2D or 3D HDF5 parameters.'
                        stop 'ModifyHDF5Statistics - ModuleHDF5Statistics - ERR02'
            
                    end select

                    ObjParameter%CurrentField => ObjParameter%CurrentField%Next               

                    ObjParameter => ObjParameter%Next               

                end do

                if (associated(ObjStatisticsTime%Next)) then

                    DT = ObjStatisticsTime%Next%Time - Me%CurrentTime

                    ObjStatisticsTime => ObjStatisticsTime%Next

                end if

                !Actualization of time            
                Me%CurrentTime = Me%CurrentTime + DT            

                call ActualizeCurrentTime(Me%ObjTime, DT, STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)                                         &
                    stop 'ModifyHDF5Statistics - ModuleHDF5Statistics - ERR03'

                !Running dependent of the last time of file
                if ((Me%CurrentTime <= ObjHDF5File%EndFieldTime) .and. (DT .ne. 0)) then
                    Running = .true.
                else
                    Running = .false.
                end if

            end do

            Count = 0

            !Kill field space for each parameter 
            ObjParameter => Me%FirstParameter

            do while(associated(ObjParameter))  

                ParameterToKill => ObjParameter
                call KillIndividualParameterFields(ObjParameter)
          
                ObjParameter    => ObjParameter%Next

            end do

            ObjHDF5File => ObjHDF5File%Next           

        end do 

    end subroutine ModifyHDF5Statistics

    !--------------------------------------------------------------------------

    subroutine  OpenAndReadHDF5File(FirstFile, ObjHDF5File, LastDT, LastStartTime)

        !Arguments-------------------------------------------------------------
        logical                                 :: FirstFile
        type(T_HDF5File), pointer               :: ObjHDF5File
        real                                    :: LastDT
        type(T_Time)                            :: LastStartTime

        !Local-----------------------------------------------------------------
        type(T_Parameter), pointer              :: ObjParameter
        integer                                 :: STAT_CALL
        integer                                 :: HDF5_READ
        real                                    :: Year, Month, Day, Hour 
        real                                    :: Minute, Second 

        !Begin-----------------------------------------------------------------

        !Objective: read relevant data for statistics from file

        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Open HDF5 file
        call ConstructHDF5 (ObjHDF5File%HDFID, trim(ObjHDF5File%Name),              & 
                            HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) then
            write(*,*) 'HDF5 file cannot be opened'//ObjHDF5File%Name                
            stop 'OpenAndReadHDF5File - ModuleHDF5Statistics - ERR01'
        end if

        !Check if there are time periods lacking
        if (ObjHDF5File%StartFieldTime > (LastStartTime + LastDT)) then
            !which DT should be? The new or the last? Suppose the last.
100         format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)
            call ExtractDate(LastStartTime, Year, Month, Day, Hour, Minute, Second)
            write(*,*) 'Data lacking from'      
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
            write(*,*) 'to'      
            call ExtractDate(ObjHDF5File%StartFieldTime, Year, Month, Day, Hour,    & 
                             Minute, Second)
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
        end if

        !Read fields data
        write(*,*)'Reading data from HDF5 file:'//trim(ObjHDF5File%Name) 

        !getting parameter fields data
        ObjParameter => Me%FirstParameter

        do while(associated(ObjParameter))

            if (FirstFile) then
            !it is the first file accessed

                nullify(ObjParameter%FirstField)
                !nullify have to be here because each new parameter has to be 
                !nullified for the first file
            
            end if

            !Allocates first field for this parameter
            allocate (ObjParameter%FirstField)
            nullify (ObjParameter%FirstField)

            call ReadParameterFields(ObjHDF5File, ObjParameter)

            ObjParameter%CurrentField => ObjParameter%FirstField

            ObjParameter => ObjParameter%Next

        end do

        LastStartTime = ObjHDF5File%EndFieldTime 
        LastDT = Me%RegularDT

        !check if there are data lacking after the last data from last file 
        if (.not. associated(ObjHDF5File%Next)) then
            !(last file)
            if (ObjHDF5File%EndFieldTime < (Me%EndTime)) then
                !no DT is considered here
                call ExtractDate(ObjHDF5File%EndFieldTime, Year, Month, Day, Hour,  &  
                                 Minute, Second)
                write(*,*) 'Data lacking from'      
                write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
                write(*,*) 'to'      
                call ExtractDate(Me%EndTime, Year, Month, Day, Hour, Minute, Second)
                write(*,fmt=100) Year, Month, Day, Hour, Minute, Second                    
            end if
        end if

        call killhdf5(ObjHDF5File%HDFID)

    end subroutine OpenAndReadHDF5File

    !--------------------------------------------------------------------------

    subroutine ReadParameterFields(ObjHDF5File, ObjParameter)

        !Arguments-------------------------------------------------------------
        type(T_Parameter), pointer              :: ObjParameter
        type(T_HDF5File), pointer               :: ObjHDF5File
          
        !Local-----------------------------------------------------------------
        integer                                 :: NewCurrentInstant
        integer                                 :: CurrentInstant
        type (T_Field), pointer                 :: NewField
        integer                                 :: Count = 1
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !read/copy parameter fields
        write(*,*)'Reading '//trim(ObjParameter%Name)//' fields'

        NewCurrentInstant = 0

        do CurrentInstant = ObjHDF5File%StartInstant, ObjHDF5File%EndInstant 

            NewCurrentInstant = NewCurrentInstant + 1
                
            !construct new fields
            call AddField(ObjParameter%FirstField, NewField)
            NewField%IDNumber = Count
            Count = Count + 1

            !get field ID, Rank and Dimensions
            !(this needs to be done for each instant)
            call GetHDF5GroupID(ObjHDF5File%HDFID, ObjParameter%Group,                  &
                                CurrentInstant, NewField%Name,                          &
                                NewField%Units, ObjParameter%Rank,                      &
                                STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_)                                                &  
            stop 'ReadParameterFields - ModuleHDF5Statistics - ERR01'
               
            NewField%Units = trim(NewField%Units)
            NewField%Date  = ObjHDF5File%InstantsArray(NewCurrentInstant)

            !get field values
            select case (ObjParameter%Rank)

                case(2)
                ! The HDF5 file contains 2D data

                    !allocate field
                    nullify (NewField%Values2D)
                    allocate(NewField%Values2D(Me%Size%ILB:Me%Size%IUB,                 &
                                               Me%Size%JLB:Me%Size%JUB))
                       
                    call HDF5SetLimits (ObjHDF5File%HDFID, Me%WorkSize%ILB,             &
                                        Me%WorkSize%IUB, Me%WorkSize%JLB,               &
                                        Me%WorkSize%JUB,                                &
                                        STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        & 
                    stop 'ReadParameterFields - ModuleHDF5Statistics - ERR02'
        
                    !read field
                    call HDF5ReadWindow(ObjHDF5File%HDFID, ObjParameter%Group,          &
                                      trim(NewField%Name),                              &
                                      Array2D      = NewField%Values2D,                 &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ReadParameterFields - ModuleHDF5Statistics - ERR03'

                case(3)
                ! The HDF5 file contains 3D data
                     
                    !allocate field
                    nullify (NewField%Values3D)
                    allocate(NewField%Values3D(Me%Size%ILB:Me%Size%IUB,                 &
                                               Me%Size%JLB:Me%Size%JUB,                 &
                                               Me%Size%KLB:Me%Size%KUB))
                        
                    call HDF5SetLimits  (ObjHDF5File%HDFID, Me%WorkSize%ILB,            &
                                         Me%WorkSize%IUB, Me%WorkSize%JLB,              &
                                         Me%WorkSize%JUB,                               &
                                         Me%WorkSize%KLB,Me%WorkSize%KUB, STAT = STAT_CALL)                
                    if (STAT_CALL .NE. SUCCESS_)                                        & 
                    stop 'ReadParameterFields - ModuleHDF5Statistics - ERR04'
        
                    !read field
                    call HDF5ReadWindow(ObjHDF5File%HDFID, ObjParameter%Group,            &
                                      trim(NewField%Name),                              &
                                      Array3D      = NewField%Values3D,                 &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ReadParameterFields - ModuleHDF5Statistics - ERR05'
 
            case default 
                    
                write(*,*)'Statistics created only for 2D or 3D HDF5 files.'
                stop 'ReadParameterFields - ModuleHDF5Statistics - ERR06'
                    
            end select

        end do

    end subroutine ReadParameterFields

    !--------------------------------------------------------------------------

    subroutine AddField (FirstField, ObjField)

        !Arguments-------------------------------------------------------------
        type (T_Field), pointer                   :: ObjField
        type (T_Field), pointer                   :: FirstField
        
        !Local-----------------------------------------------------------------
        type (T_Field), pointer                   :: NewField
        type (T_Field), pointer                   :: PreviousField
        
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

    subroutine CalculateHDF5Statistics2D(Field2D, StatisticsID)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer :: Field2D
        integer                       :: StatisticsID
        !Local-----------------------------------------------------------------
        integer                       :: MethodStatistic, Value2DStat2D
        integer                       :: STAT_CALL

        !----------------------------------------------------------------------

        call GetStatisticMethod (StatisticsID, MethodStatistic, STAT = STAT_CALL)                                     
                                                                                    
        if (STAT_CALL /= SUCCESS_)                                                      & 
            call SetError (FATAL_, INTERNAL_,                                           &
                           'CalculateHDF5Statistics2D - ModuleHDF5Statistics - ERR01')
                                                                                    
        call GetStatisticParameters (StatisticsID,                                      &
                                     Value2DStat2D = Value2DStat2D,                     &
                                     STAT          = STAT_CALL)                        
                                                                                    
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, INTERNAL_,                                           &
                           'CalculateHDF5Statistics2D - ModuleHDF5Statistics - ERR02')
                                                                                    
                                                                                    
        if (MethodStatistic /= Value2DStat2D)                                           &
            call SetError (FATAL_, INTERNAL_,                                           & 
                           'CalculateHDF5Statistics2D - ModuleHDF5Statistics - ERR03')
                                                                                    
                                                                                    
        call ModifyStatistic (StatisticsID,                                             &
                              Value2D       = Field2D,                                  &
                              WaterPoints2D = Me%Mapping%IntegerValues2D,               &
                              STAT          = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, INTERNAL_,                                           & 
                           'CalculateHDF5Statistics2D - ModuleHDF5Statistics - ERR04')

    end subroutine CalculateHDF5Statistics2D

    !--------------------------------------------------------------------------

    subroutine CalculateHDF5Statistics3D(Field3D, StatisticsID)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:,:), pointer :: Field3D
        integer                         :: StatisticsID
        !Local-----------------------------------------------------------------
        integer                                 :: MethodStatistic, Value3DStat3D 
        integer                                 :: Value3DStatLayers
        integer                                 :: STAT_CALL
        integer                                 :: LayerDefinition                      
        integer                                 :: Depth, Layer
        integer                                 :: LayersNumber, ln

        !----------------------------------------------------------------------

        call GetStatisticMethod (StatisticsID, MethodStatistic, STAT = STAT_CALL)                                     
                                                                                    
        if (STAT_CALL /= SUCCESS_)                                                      & 
            call SetError (FATAL_, INTERNAL_,                                           &
                           'CalculateHDF5Statistics3D - ModuleHDF5Statistics - ERR01')
                                                                                    
        call GetStatisticParameters (StatisticsID,                                      &
                                     Value3DStat3D = Value3DStat3D,                     &
                                     STAT          = STAT_CALL)                        
                                                                                    
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, INTERNAL_,                                           &
                           'CalculateHDF5Statistics3D - ModuleHDF5Statistics - ERR02')
                           
        call GetStatisticParameters (StatisticsID,                                      &
                                     Value3DStatLayers = Value3DStatLayers,             &
                                     Depth             = Depth,                         &
                                     Layer             = Layer,                         &
                                     STAT          = STAT_CALL)                        

        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, INTERNAL_,                                           &
                           'CalculateHDF5Statistics3D - ModuleHDF5Statistics - ERR03')
                                                                                    
        if ((MethodStatistic /= Value3DStat3D) .and.                                    & 
            (MethodStatistic /= Value3DStatLayers))                                     &
            call SetError (FATAL_, INTERNAL_,                                           &
                           'CalculateHDF5Statistics3D - ModuleHDF5Statistics - ERR04')
        
        if (MethodStatistic == Value3DStatLayers) then


            !Get layer specification
            !nullify(MLD_Surf)
            !nullify(KFloorZ)
            
            call GetStatisticLayersNumber(StatisticsID, LayersNumber, STAT = STAT_CALL)
                                            
            if (STAT_CALL /= SUCCESS_)                                              &
                    call SetError (FATAL_, KEYWORD_,                                & 
                                   'CalculateHDF5Statistics3D - ModuleHDF5Statistics - ERR05')

            do ln=1, LayersNumber

                call GetStatisticLayerDef(StatisticsID, ln, LayerDefinition, STAT = STAT_CALL)
                                           
                if (STAT_CALL /= SUCCESS_)                                          &
                        call SetError (FATAL_, KEYWORD_,                            & 
                                       'CalculateHDF5Statistics3D - ModuleHDF5Statistics - ERR06')
                                       

            !    !Statistic of properties values along the bottom 
                if (LayerDefinition == Layer) then 
                
                    call AddStatisticLayers (StatisticID    = StatisticsID,         &
                                             Value3D        = Field3D,              &
                                             WaterPoints3D  = Me%Mapping%IntegerValues3D, &
                                             DZ3D           = Me%DZ3D,              &
                                             LayerNumber    = ln,                   &
                                             STAT= STAT_CALL) 

                    if (STAT_CALL /= SUCCESS_)                                      &
                        call SetError (FATAL_, KEYWORD_,                            & 
                                       'OutPut_Statistics - ModuleHDF5Statistics - ERR07')


            !    
                else if (LayerDefinition == Depth) then 

                    call AddStatisticLayers (StatisticID    = StatisticsID,         &
                                             Value3D        = Field3D,              &
                                             WaterPoints3D  = Me%Mapping%IntegerValues3D, &
                                             DZ3D           = Me%DZ3D,              &
                                             LayerNumber    = ln,                   &
                                             STAT= STAT_CALL) 
                
                endif
           
            enddo

        endif
                                                                                    
        call ModifyStatistic (StatisticsID,                                             &
                              Value3D       = Field3D,                                  &
                              WaterPoints3D = Me%Mapping%IntegerValues3D,               &
                              STAT          = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, INTERNAL_,                                           &
                           'CalculateHDF5Statistics3D - ModuleHDF5Statistics - ERR06')

    end subroutine CalculateHDF5Statistics3D

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillHDF5Statistics

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        type(T_HDF5File), pointer               :: HDF5FileX, HDF5FileToKill
        type(T_Parameter), pointer              :: ParameterX, ParameterToKill
        integer                                 :: STAT_CALL
        real                                    :: ElapsedSeconds, TotalCPUTime

        !------------------------------------------------------------------------

        !Kill the mapping field
        if (associated(Me%Mapping%IntegerValues3D)) then
            deallocate(Me%Mapping%IntegerValues3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'KillHDF5Statistics - ModuleHDF5Statistics - ERR01'  
            nullify(Me%Mapping%IntegerValues3D)
        endif

        if (associated(Me%Mapping%IntegerValues2D)) then
            deallocate(Me%Mapping%IntegerValues2D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'KillHDF5Statistics - ModuleHDF5Statistics - ERR02'  
            nullify(Me%Mapping%IntegerValues2D)
        endif

        !Kill all HDF5Files from the list of the ones relevant for Statistics
        HDF5FileX=> Me%FirstStatHDF5File

        do while(associated(HDF5FileX))  

            HDF5FileToKill  => HDF5FileX
            HDF5FileX       => HDF5FileX%Next
            call KillIndividualHDF5File(HDF5FileToKill)

        end do
        nullify(Me%FirstStatHDF5File)

        !Kill all parameters from the list for Statistics
        ParameterX=> Me%FirstParameter

        do while(associated(ParameterX))  

            !copy statistics for HDF5 file
            call KillStatistic (ParameterX%Statistics%ID, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                call SetError(FATAL_, INTERNAL_,                                    &
                'KillHDF5Statistics - ModuleHDF5Statistics - ERR03')

            ParameterToKill => ParameterX
            ParameterX      => ParameterX%Next
            call KillIndividualParameter(ParameterToKill)

        end do
        nullify(Me%FirstParameter)

        call HDF5FlushMemory (Me%ObjStatHDF5, STAT = STAT_CALL)          
        if (STAT_CALL /= SUCCESS_)                                                  &
            call SetError(FATAL_, INTERNAL_,                                        &
            'KillHDF5Statistics - ModuleHDF5Statistics - ERR04')
            
        if (associated(Me%DZ3D)) then
            deallocate(Me%DZ3D)
        endif

        deallocate(Me)
        nullify(Me)

        call date_and_time(Values = F95Time)
        call SetDate      (FinalSystemTime,   float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)
        call cpu_time(TotalCPUTime)
        ElapsedSeconds = FinalSystemTime - InitialSystemTime

        call ShutdownMohid ("HDF5 Statistics", ElapsedSeconds, TotalCPUTime)


        !------------------------------------------------------------------------

    end subroutine KillHDF5Statistics

    !--------------------------------------------------------------------------

    subroutine KillIndividualHDF5File(HDF5ToDispose)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5ToDispose

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Deallocate the InstantsArray
        if (associated(HDF5ToDispose%InstantsArray)) then 
            deallocate(HDF5ToDispose%InstantsArray, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'KillIndividualHDF5File - ModuleHDF5Statistics - ERR01'               
            nullify(HDF5ToDispose%InstantsArray)
        end if

        deallocate(HDF5ToDispose, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
        stop 'KillIndividualHDF5File - ModuleHDF5Statistics - ERR02'
        nullify(HDF5ToDispose)

    end subroutine KillIndividualHDF5File

    !--------------------------------------------------------------------------        

    subroutine KillIndividualParameter(ParameterToDispose)

        !Arguments-------------------------------------------------------------
        type(T_Parameter), pointer                  :: ParameterToDispose

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !The fields of the parameters have already been killed
        deallocate(ParameterToDispose, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
        stop 'KillIndividualParameter - ModuleHDF5Statistics - ERR01'
        nullify(ParameterToDispose)

    end subroutine KillIndividualParameter

    !------------------------------------------------------------------------

    subroutine KillGridFields

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Connection X
        deallocate(Me%ConnectionX%RealValues2D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'KillGridFields - ModuleHDF5Statistics - ERR01'  
        nullify(Me%ConnectionX%RealValues2D)

        !Connection Y
        deallocate(Me%ConnectionY%RealValues2D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'KillGridFields - ModuleHDF5Statistics - ERR02'  
        nullify(Me%ConnectionY%RealValues2D)

        !Latitude 
        deallocate(Me%Latitude%RealValues2D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'KillGridFields - ModuleHDF5Statistics - ERR03'  
        nullify(Me%Latitude%RealValues2D)

        !Longitude 
        deallocate(Me%Longitude%RealValues2D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'KillGridFields - ModuleHDF5Statistics - ERR04'  
        nullify(Me%Longitude%RealValues2D)

        !Bathymetry 
        deallocate(Me%Bathymetry%RealValues2D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'KillGridFields - ModuleHDF5Statistics - ERR05'  
        nullify(Me%Bathymetry%RealValues2D)

    end subroutine KillGridFields

    !--------------------------------------------------------------------------

    subroutine KillIndividualParameterFields(ParameterToDispose)

        !Arguments-------------------------------------------------------------
        type(T_Parameter), pointer                  :: ParameterToDispose

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_Field), pointer                      :: FieldToKill, CurrentField

        !Begin-----------------------------------------------------------------

        CurrentField => ParameterToDispose%FirstField

        do while(associated(CurrentField))

            FieldToKill => CurrentField
            CurrentField => CurrentField%Next

            if (associated(FieldToKill%Values2D)) then
                deallocate(FieldToKill%Values2D, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'KillIndividualParameterFields - ModuleHDF5Statistics - ERR01'  
                nullify(FieldToKill%Values2D)
            end if

            if (associated(FieldToKill%Values3D)) then
                deallocate(FieldToKill%Values3D, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'KillIndividualParameterFields - ModuleHDF5Statistics - ERR02'  
                nullify(FieldToKill%Values3D)
            end if

        end do 

    end subroutine KillIndividualParameterFields

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module ModuleHDF5Statistics









