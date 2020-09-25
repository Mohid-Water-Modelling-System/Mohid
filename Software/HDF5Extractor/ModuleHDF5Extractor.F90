!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : HDF5Extractor
! MODULE        : HDF5Extractor
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : October 2005
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Module to create a HDF5 file from a time window of a supplied 
!                 HDF5 file.
!
!------------------------------------------------------------------------------

!DataFile
!
!   FILENAME                : char                  [-]         !Name of HDF5 file from which
!                                                               !a time window is to be 
!                                                               !extracted
!
!   OUTPUTFILENAME          : char                  [-]         !Name of the output HDF5 file
!
!   START_TIME              : YYYY MM DD HH MM SS   [-]         !Start date of time window
!   END_TIME                : YYYY MM DD HH MM SS   [-]         !End date of time window
!
!
!   XY_WINDOW_OUTPUT        : logical               [false]     !Spatial Window from each you want to extract data
!
!   XY_WINDOW_LIMITS        : 4*integer             [4*FillValueInt] ! ilb, jlb, iub, jub
!
!   LAYERS_OUTPUT           : logical               [false]     !Check if the user wants to extract spefici layers
!
!   LAYERS_MIN_MAX          : 2*integer             [2*FillValueInt] ! klb, kub
!
!   INTERVAL                : logical               [false/0]   !Extraction by time interval
!
!   DT_INTERVAL             : integer               [-]         !Time interval (seconds) for extraction
!
!   CONVERT_V3_TO_V4        : logical               [false]     ! If you are extracting from a V3 file,
                                                                ! this is the option of convert it to V4 format
!   <BeginParameter>
!   HDF_GROUP_V3            : char                  [-]         !Path of the group of property in HDF5 file in V3 format
                                                                ! (Only if CONVERT_V3_TO_V4) = 1 
!   PROPERTY_V3             : char                  [-]         !Property name in HD5 file in V3 format  
                                                                ! (Only if CONVERT_V3_TO_V4) = 1 
!   HDF_GROUP               : char                  [-]         !Path of the group of property in HDF5 file 
!   PROPERTY                : char                  [-]         !Property name (should be equal to one specified 
!   <EndParameter>                                              !in ModuleGlobalData)
!                                                               !(specify one block for each property)      

Module ModuleHDF5Extractor

    use HDF5
    use ModuleGlobalData        

    use ModuleHDF5,              only : GetHDF5FileAccess, ConstructHDF5,       &
                                        GetHDF5GroupNumberOfItems,              &
                                        HDF5SetLimits, HDF5ReadData,            &
                                        GetHDF5GroupID, KillHDF5,               &
                                        HDF5WriteData, HDF5FlushMemory,         &
                                        GetHDF5FileID

    use ModuleEnterData,         only : ConstructEnterData, KillEnterData,      &
                                        GetData, ExtractBlockFromBuffer,        &
                                        Block_Unlock
                                        
    use ModuleTime
    
    use ModuleDrawing,           only : ArrayPolygonWindow

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartHDF5Extractor
    private ::      ConstructExtractHDF5
    private ::          ReadKeywords
    private ::              ReadParameters
    private ::                  AddItem
    private ::                  ConstructParameters
    private ::          OpenAndDateHDF5File
    private ::              HDF5TimeInstant
    private ::              AddIntervalHDF5
    private ::              InquireParameters
    private ::              InquireSubGroup
    private ::                  ReadItemField
    private ::          PrepareOutput
    private ::              OpenOutputFile

    !Selector                  
    
    !Modifier
    public  :: ModifyExtractHDF5
    private ::      OpenAndReadHDF5File
    private ::          OutputInstants
    private ::          ReadTimeDependentFields
    private ::              AddField
    private ::      WriteItemFields
    private ::      KillItemFields
    private ::          KillIndividualField

    !Destructor
    public  :: KillExtractHDF5
    private ::      KillIndividualItem                   

    !Management
    
    !Interfaces----------------------------------------------------------------

    !Types---------------------------------------------------------------------

    ! Definition of IntervalHDF5
    type       T_IntervalHDF5
        character(len=PathLength)                   :: Name
        integer                                     :: ObjEnterData         = 0
        type(T_IntervalHDF5),       pointer         :: Next  => null()
    end type T_IntervalHDF5

    ! Definition of Interval
    type       T_Interval
        logical                                     :: Interval_ON
        integer                                     :: DT
        type(T_IntervalHDF5),       pointer         :: FirstHDF5
        integer                                     :: NumberIntervalHDF5   = 0
        integer                                     :: CountIntervalHDF5    = 0
    end type T_Interval

    ! Definition of SpatialWindow
    type       T_SpatialWindow
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        logical                                     :: XY_ON, Layers_ON, SurfaceLayer
        logical                                     :: XY_COORD
        real                                        :: Xmin, Ymin, Xmax, Ymax
    end type T_SpatialWindow

    ! Definition of type T_Field
    type       T_Field
        character(len=StringLength)                 :: Name
        real, dimension(:      ),   pointer         :: Values1D
        integer, dimension(:   ),   pointer         :: IntValues1D !For integer valued
        real, dimension(:,:  ),     pointer         :: Values2D
        integer, dimension(:,: ),   pointer         :: IntValues2D !For integer valued
        real, dimension(:,:,:),     pointer         :: Values3D
        integer, dimension(:,:,:),  pointer         :: IntValues3D !For integer valued
        type(T_Field),              pointer               :: Next  => null()
    end type  T_Field

    ! Definition of type T_Item
    type       T_Item
        character(len=PathLength)                   :: Group
        character(len=StringLength)                 :: Name
        character(len=PathLength)                   :: GroupV3
        character(len=StringLength)                 :: NameV3
        character(len=PathLength)                   :: GroupV4  = null_str
        character(len=StringLength)                 :: NameV4   = null_str
        integer                                     :: Rank
        character(len=StringLength)                 :: Units
        integer, dimension(7)                       :: Dimensions
        logical                                     :: TimeDependent
        integer(HID_T)                              :: NumType     = 0
        type(T_Field),              pointer         :: FirstField
        type(T_Item),               pointer         :: Next  => null()
    end type T_Item

    ! Definition of type T_ExtractHDF5   
    private :: T_ExtractHDF5
    type       T_ExtractHDF5
        logical                                     :: ConvertV3toV4
        integer                                     :: InstanceID
        type (T_Size3D)                             :: WorkSize
        character(PathLength)                       :: DataFile
        character(PathLength)                       :: HDF5File, OutputFile, TimeGroup
        integer                                     :: ObjEnterData         = 0
        integer                                     :: ObjOutputFile        = 0
        integer                                     :: ObjHDF5File          = 0
        type(T_Time)                                :: StartTime, EndTime
        type(T_SpatialWindow)                       :: SpatialWindow
        type(T_Interval)                            :: Interval
        integer                                     :: HDF5StartInstant
        integer                                     :: HDF5EndInstant
        type(T_Time), dimension(:), pointer         :: InstantsArray
        integer                                     :: InstantsNumber
        type(T_Item),               pointer         :: FirstIndependentItem
        type(T_Item),               pointer         :: FirstDependentItem
        type(T_Item),               pointer         :: FirstParameter
        character(StringLength)                     :: CurrentGroup, LastSubGroup
        integer(HID_T)                              :: HDFNativeReal, HDFNativeInteger
        integer(HID_T)                              :: FileID
        type(T_ExtractHDF5),        pointer         :: Next  => null()
    end type  T_ExtractHDF5

    !Global Module Variables
    type (T_ExtractHDF5),           pointer         :: FirstObjExtractHDF5
    type (T_ExtractHDF5),           pointer         :: Me

    integer                                         :: mExtractHDF5_ = 0 !just to compile

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartHDF5Extractor(ObjHDF5ExtractorID, DataFile)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjHDF5ExtractorID 
        character(PathLength), intent(IN)               :: DataFile

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------

        !------------------------------------------------------------------------

        !Assures nullification of the global variable
        allocate(Me)

        !Returns ID
        ObjHDF5ExtractorID          = 1

        !Atribute the name of data file            
        Me%DataFile = DataFile

        call ConstructExtractHDF5

        !----------------------------------------------------------------------

    end subroutine StartHDF5Extractor
 
    !--------------------------------------------------------------------------

    subroutine ConstructExtractHDF5

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------

        !------------------------------------------------------------------------

        ! Read keywords file
        call ReadKeywords
        
        !Construct input and output hdf5 file
        call ConstructInAndOut

        ! Open supplied HDF5 file
        call OpenAndDateHDF5File
        
        ! Open HDF5 output file
        call PrepareOutput
        ! (and copy the time independent items' values to the new file)

        !----------------------------------------------------------------------

    end subroutine ConstructExtractHDF5

    !--------------------------------------------------------------------------
    
    subroutine ConstructInAndOut
    
        !Arguments---------------------------------------------------------------

    
        !Local-------------------------------------------------------------------
        integer                                     :: HDF5_READ, STAT_CALL
        logical                                     :: exist
                
        !Begin-------------------------------------------------------------------        
    
       !Verifies if file exists
        inquire(FILE = Me%HDF5File, EXIST = exist)
        if (.not. exist) then
            write(*,*)
            write(*,*)'HDF5 file does not exist:'//trim(Me%HDF5File)
            stop 'ConstructInAndOut - ModuleHDF5Extractor - ERR10'
        endif

        !(there has been a kill of this file before)
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Open HDF5 file
        call ConstructHDF5 (Me%ObjHDF5File, trim(Me%HDF5File),                  & 
                            HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) then
            write(*,*) 'HDF5 file cannot be opened'//Me%HDF5File                
            stop 'ConstructInAndOut - ModuleHDF5Extractor - ERR20'
        end if

    
    end subroutine ConstructInAndOut
    
    !--------------------------------------------------------------------------
  
    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag
        integer, dimension(4)                       :: Aux4
        integer, dimension(2)                       :: Aux2
        real,    dimension(4)                       :: AuxR
        !Begin-----------------------------------------------------------------
        
        call ConstructEnterData (Me%ObjEnterData, Me%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'ReadKeywords - ModuleHDF5Extractor - ERR10'

        ! Obtain input HDF5 file name
        call GetData(Me%HDF5File, Me%ObjEnterData, iflag,                       &
                     keyword      = 'FILENAME',                                 &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'HDF5Extractor',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                              &
        stop 'ReadKeywords - ModuleHDF5Extractor - ERR20'   

        ! Obtain output HDF5 file name
        call GetData(Me%OutputFile, Me%ObjEnterData, iflag,                     &
                     keyword      = 'OUTPUTFILENAME',                           &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'HDF5Extractor',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                              &
        stop 'ReadKeywords - ModuleHDF5Extractor - ERR30'   


        Aux4(:) = FillValueInt
        Aux2(:) = FillValueInt
        AuxR(:) = FillValueReal

        ! Check if the user want to extract a XY spatial window
        call GetData(Me%SpatialWindow%XY_ON, Me%ObjEnterData, iflag,            &
                     keyword      = 'XY_WINDOW_OUTPUT',                         &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'HDF5Extractor',                            &
                     default      = .false.,                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'ReadKeywords - ModuleHDF5Extractor - ERR40'  


        if (Me%SpatialWindow%XY_ON) then

            call GetData(Me%SpatialWindow%XY_COORD, Me%ObjEnterData, iflag,     &
                         keyword      = 'XY_COORD',                             &
                         SearchType   = FromFile,                               &
                         ClientModule = 'HDF5Extractor',                        &
                         default      = .false.,                                &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ReadKeywords - ModuleHDF5Extractor - ERR50'  

            
            if (Me%SpatialWindow%XY_COORD) then
            
                call GetData(AuxR, Me%ObjEnterData, iflag,                          &
                             keyword      = 'XY_WINDOW_LIMITS',                     &
                             SearchType   = FromFile,                               &
                             ClientModule = 'HDF5Extractor',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ReadKeywords - ModuleHDF5Extractor - ERR60'  

                if (iflag /= 4)                                                     &
                stop 'ReadKeywords - ModuleHDF5Extractor - ERR70'  


                Me%SpatialWindow%Xmin = AuxR(1)
                Me%SpatialWindow%Ymin = AuxR(2)
                Me%SpatialWindow%Xmax = AuxR(3)
                Me%SpatialWindow%Ymax = AuxR(4)
            
            
            else

                call GetData(Aux4, Me%ObjEnterData, iflag,                          &
                             keyword      = 'XY_WINDOW_LIMITS',                     &
                             SearchType   = FromFile,                               &
                             ClientModule = 'HDF5Extractor',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ReadKeywords - ModuleHDF5Extractor - ERR80'  

                if (iflag /= 4)                                                     &
                stop 'ReadKeywords - ModuleHDF5Extractor - ERR90'  


                Me%SpatialWindow%ILB = Aux4(1)
                Me%SpatialWindow%JLB = Aux4(2)
                Me%SpatialWindow%IUB = Aux4(3)
                Me%SpatialWindow%JUB = Aux4(4)

            endif
        endif

        ! Check if the user want to extract some layers
        call GetData(Me%SpatialWindow%Layers_ON, Me%ObjEnterData, iflag,        &
                     keyword      = 'LAYERS_OUTPUT',                            &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'HDF5Extractor',                            &
                     default      = .false.,                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'ReadKeywords - ModuleHDF5Extractor - ERR100'  


        if (Me%SpatialWindow%Layers_ON) then
        
            call GetData(Me%SpatialWindow%SurfaceLayer, Me%ObjEnterData, iflag,         &
                         keyword      = 'SURFACE_LAYER',                                &
                         SearchType   = FromFile,                                       &
                         ClientModule = 'HDF5Extractor',                                &
                         default      = .false.,                                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ReadKeywords - ModuleHDF5Extractor - ERR110'  
            
            if (.not. Me%SpatialWindow%SurfaceLayer) then

                call GetData(Aux2, Me%ObjEnterData, iflag,                              &
                             keyword      = 'LAYERS_MIN_MAX',                           &
                             SearchType   = FromFile,                                   &
                             ClientModule = 'HDF5Extractor',                            &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                stop 'ReadKeywords - ModuleHDF5Extractor - ERR120'  

                if (iflag /= 2)                                                         &
                stop 'ReadKeywords - ModuleHDF5Extractor - ERR130'  

                Me%SpatialWindow%KLB = Aux2(1)
                Me%SpatialWindow%KUB = Aux2(2)
                
            endif                

        endif

        ! Check if the user want to extract by time interval
        call GetData(Me%Interval%Interval_ON, Me%ObjEnterData, iflag,           &
                     keyword      = 'INTERVAL',                                 &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'HDF5Extractor',                            &
                     default      = .false.,                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'ReadKeywords - ModuleHDF5Extractor - ERR140'  

        if (Me%Interval%Interval_ON) then
            !get DT
            call GetData(Me%Interval%DT, Me%ObjEnterData, iflag,                &
                         keyword      = 'DT_INTERVAL',                          &
                         SearchType   = FromFile,                               &
                         ClientModule = 'HDF5Extractor',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_ .or. iflag == 0) then
                write(*,*) 'A DT_INTERVAL must be provided!'
                stop 'ReadKeywords - ModuleHDF5Extractor - ERR150'  
            endif

        else
            ! Obtain the start and end times for the time window
            ! Start Time
            call GetData(Me%StartTime, Me%ObjEnterData, iflag,                  &
                         keyword      = 'START_TIME',                           &
                         SearchType   = FromFile,                               &
                         ClientModule = 'HDF5Extractor',                        &
                         STAT         = STAT_CALL)
            
            if (iflag == 0) then
                
                call GetData(Me%StartTime, Me%ObjEnterData, iflag,                  &
                             keyword      = 'START',                           &
                             SearchType   = FromFile,                               &
                             ClientModule = 'HDF5Extractor',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_ .or. iflag == 0) then
                stop 'ReadKeywords - ModuleHDF5Extractor - ERR160a'
                endif
            else
                if (STAT_CALL /= SUCCESS_) then    
                stop 'ReadKeywords - ModuleHDF5Extractor - ERR160b'  
            endif
            endif
            
            ! End Time 
            call GetData(Me%EndTime,   Me%ObjEnterData, iflag,                  &
                         keyword      = 'END_TIME',                             &
                         SearchType   = FromFile,                               &
                         ClientModule = 'HDF5Extractor',                        &
                         STAT         = STAT_CALL)
            if (iflag == 0) then
                call GetData(Me%EndTime,   Me%ObjEnterData, iflag,                  &
                             keyword      = 'END',                             &
                             SearchType   = FromFile,                               &
                             ClientModule = 'HDF5Extractor',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_ .or. iflag == 0)  then
                stop 'ReadKeywords - ModuleHDF5Extractor - ERR170a'
                endif
            else
                if (STAT_CALL /= SUCCESS_) then
                stop 'ReadKeywords - ModuleHDF5Extractor - ERR170b'       
            endif
            endif
            
            ! Verifies Time Variables
            if (Me%EndTime .lt. Me%StartTime) then
                write (*,*) 'End Time is BEFORE Start Time'
                write (*,*) 'Module :','HDF5Extractor'
                stop 'ReadKeywords - ModuleHDF5Extractor - ERR180'
            endif

        endif

        ! Checks if it is to convert a v3 file to v4 format
        call GetData(Me%ConvertV3toV4,                                          &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromFile,                                   &
                     keyword      = 'CONVERT_V3_TO_V4',                         &
                     default      = .false.,                                    &
                     ClientModule = 'HDF5Extractor',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
        stop 'ConstructParameters - ModuleHDF5Extractor - ERR190'


        call GetData(Me%TimeGroup,                                              &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromFile,                                   &
                     keyword      = 'TIME_GROUP',                               &
                     default      = "Time",                                    &
                     ClientModule = 'HDF5Extractor',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
        stop 'ConstructParameters - ModuleHDF5Extractor - ERR200'

        

        ! Obtain parameters
        call ReadParameters 

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadKeywords - ModuleHDF5Extractor - ERR210'
        end if

    end subroutine ReadKeywords

    !--------------------------------------------------------------------------

    subroutine ReadParameters

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, ClientNumber
        logical                                     :: BlockFound
        logical                                     :: AtLeastOneBlock = .false.
        type(T_Item), pointer                       :: NewParameter

        !Begin-----------------------------------------------------------------

        ! Obtain time independent itens for new HDF5
do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<BeginParameter>',   &
                                        block_end       = '<EndParameter>',     &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  

                    AtLeastOneBlock = .true.
                    
                    call AddItem(Me%FirstParameter, NewParameter)

                    call ConstructParameters(NewParameter)

                    nullify(NewParameter)

                else cd2
                    call Block_Unlock(Me%ObjEnterData,                          & 
                                      ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ReadParameters - ModuleHDF5Extractor - ERR10'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ReadParameters - ModuleHDF5Extractor - ERR20'
            else cd1
                stop 'ReadParameters - ModuleHDF5Extractor - ERR30'
            end if cd1

        end do do1

        if (.not. AtLeastOneBlock) then                                            
            write(*,*)
            write(*,*) 'No parameter block is indicated in input file.'
            write(*,*) 'Must be indicated at least one parameter block.'
            stop 'ReadParameters - ModuleHDF5Extractor - ERR40'
        end if

    end subroutine ReadParameters

    !--------------------------------------------------------------------------

    subroutine ConstructParameters (NewParameter)

        !Arguments-------------------------------------------------------------
        type(T_Item), pointer                       :: NewParameter

        !External--------------------------------------------------------------
        integer                                     :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        
        ! Obtain parameter name
        call GetData(NewParameter%Name,                                         &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'PROPERTY',                                 &
                     ClientModule = 'HDF5Extractor',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
        stop 'ConstructParameters - ModuleHDF5Extractor - ERR10'

        if (.not.CheckPropertyName(NewParameter%Name)) then
            write(*,*)
            write(*,*) 'The property name is not recognised by the model: '     &
                        //trim(NewParameter%Name)
            write(*,*) 'ConstructParameters - ModuleHDF5Extractor - WARN10' 
        end if

        ! Obtain parameter group
        call GetData(NewParameter%Group,                                        &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'HDF_GROUP',                                &
                     ClientModule = 'HDF5Extractor',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
        stop 'ConstructParameters - ModuleHDF5Extractor - ERR20'

        if (Me%ConvertV3toV4) then

            ! Obtain V3 parameter name
            call GetData(NewParameter%NameV3,                                   &
                         Me%ObjEnterData, iflag,                                &
                         SearchType   = FromBlock,                              &
                         keyword      = 'PROPERTY_V3',                          &
                         default      = NewParameter%Name,                      &
                         ClientModule = 'HDF5Extractor',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'ConstructParameters - ModuleHDF5Extractor - ERR30'

            ! Obtain V3 parameter group
            call GetData(NewParameter%GroupV3,                                  &
                         Me%ObjEnterData, iflag,                                &
                         SearchType   = FromBlock,                              &
                         keyword      = 'HDF_GROUP_V3',                         &
                         default      = NewParameter%Group,                     &
                         ClientModule = 'HDF5Extractor',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'ConstructParameters - ModuleHDF5Extractor - ERR40'

            NewParameter%NameV4  = NewParameter%Name
            NewParameter%GroupV4 = NewParameter%Group

            NewParameter%Name  = NewParameter%NameV3
            NewParameter%Group = NewParameter%GroupV3

        end if

    end subroutine ConstructParameters 

    !--------------------------------------------------------------------------

    subroutine OpenAndDateHDF5File

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
!        logical                                     :: exist
!        integer                                     :: HDF5_READ
        integer                                     :: STAT_CALL
        real                                        :: Year, Month, Day, Hour 
        real                                        :: Minute, Second
        type(T_Time), dimension(:), pointer         :: AuxInstantsArray 
        integer                                     :: CurrentInstant, AuxNumberInstants
        type(T_Time)                                :: HDF5StartFieldTime
        type(T_Time)                                :: HDF5EndFieldTime
        type(T_Time)                                :: HDF5StartTime, HDF5EndTime
        type(T_Time)                                :: StartIntervalTime
        character(PathLength)                       :: AuxOutputName, AuxChar
        type(T_IntervalHDF5),       pointer         :: NewIntervalHDF5

        !Begin-----------------------------------------------------------------

!        !Verifies if file exists
!        inquire(FILE = Me%HDF5File, EXIST = exist)
!        if (.not. exist) then
!            write(*,*)
!            write(*,*)'HDF5 file does not exist:'//trim(Me%HDF5File)
!            stop 'OpenAndDateHDF5File - ModuleHDF5Extractor - ERR10'
!        endif
!
!        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
!
!        !Open HDF5 file
!        call ConstructHDF5 (Me%ObjHDF5File, trim(Me%HDF5File),                  &
!                            HDF5_READ, STAT = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                            &
!        stop 'OpenAndDateHDF5File - ModuleHDF5Extractor - ERR20'

        !Obtain start and end times of HDF5 file
        !(obtain number of instants) 
        call GetHDF5GroupNumberOfItems(Me%ObjHDF5File, "/"//trim(Me%TimeGroup), &
                                       Me%InstantsNumber,                       & 
                                       STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            & 
        stop 'OpenAndDateHDF5File - ModuleHDF5Extractor - ERR30'

        !(obtain HDF5 start time)
        HDF5StartTime = HDF5TimeInstant(1)

        !(obtain HDF5 end time)
        HDF5EndTime = HDF5TimeInstant(Me%InstantsNumber)

        if (.not. Me%Interval%Interval_ON) then

            !See if the file contains required period 
            !Check start and end time
            if ((HDF5EndTime >= Me%StartTime) .and.                             &
                (HDF5EndTime <= Me%EndTime)) then

                !End time is between start and end times of file

                !Start time is after start time of file
                HDF5StartFieldTime = Me%StartTime
                HDF5EndFieldTime = HDF5EndTime


                if (HDF5StartTime >= Me%StartTime) then

                    !No extraction is needed: all data included in requested period
                    write(*,*)
                    write(*,*) 'All HDF5 data is already included in requested period.'
                    write(*,*) 'Time Extraction is not needed only spatial!'      
                    write(*,*) 'OpenAndDateHDF5File - ModuleHDF5Extractor - WARN10'                

                endif

                !Check if data is lacking for the end of requested period
                if (HDF5EndTime < Me%EndTime) then

                    !Data lacking for the end of period
                    write(*,*)
                    call ExtractDate(HDF5EndTime, Year, Month, Day, Hour,       &  
                                     Minute, Second)        
                    write(*,*) 'Data lacking from'      
                    write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
                    write(*,*) 'to'      
                    call ExtractDate(Me%EndTime, Year, Month, Day, Hour,        &
                                     Minute, Second)
                    write(*,fmt=100) Year, Month, Day, Hour, Minute, Second

                endif

            else if ((HDF5StartTime >= Me%StartTime) .and.                      &
                     (HDF5StartTime <= Me%EndTime)) then

                !End time is before end time of file
                HDF5StartFieldTime = HDF5StartTime
                HDF5EndFieldTime = Me%EndTime

                !Check if data is lacking for the beginning of requested period
                if (HDF5StartTime > Me%StartTime) then

                    !Data lacking for the beginning of period
                    write(*,*)
                    call ExtractDate(Me%StartTime, Year, Month, Day, Hour,      &  
                                     Minute, Second)        
                    write(*,*) 'Data lacking from'      
                    write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
                    write(*,*) 'to'      
                    call ExtractDate(HDF5StartTime, Year, Month, Day, Hour,     &
                                     Minute, Second)
                    write(*,fmt=100) Year, Month, Day, Hour, Minute, Second

                endif

            else if ((HDF5StartTime < Me%StartTime) .and.                       &
                     (HDF5EndTime > Me%EndTime)) then

                !Period required is contained in file
                HDF5StartFieldTime = Me%StartTime
                HDF5EndFieldTime = Me%EndTime

            else

                !HDF5 file does not contains required period
                write(*,*)
                call ExtractDate(Me%StartTime, Year, Month, Day, Hour,          &  
                                 Minute, Second)        
                write(*,*) 'Data lacking from'      
                write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
                write(*,*) 'to'      
                call ExtractDate(Me%EndTime, Year, Month, Day, Hour, Minute, Second)
                write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
                write(*,*) 'There is no requested data in the HDF5 file provided.'      
                stop 'OpenAndDateHDF5File - ModuleHDF5Extractor - ERR40'                

            endif

        endif

        !Get useful time information from file:
        !Set instant array for this file 
        allocate(AuxInstantsArray(1:Me%InstantsNumber))

        !Fill array with instants
        do CurrentInstant = 1, Me%InstantsNumber

           AuxInstantsArray(CurrentInstant) =                                   &
                HDF5TimeInstant(CurrentInstant)

        end do

        if (.not. Me%Interval%Interval_ON) then

            !Get start and end instants for this file
            !select time window begin
            do CurrentInstant = 1, Me%InstantsNumber

                if (AuxInstantsArray(CurrentInstant)                            &
                    .ge. HDF5StartFieldTime) then

                    Me%HDF5StartInstant = CurrentInstant
                    HDF5StartFieldTime = HDF5TimeInstant(CurrentInstant)
        
                    exit

                end if

            end do

            !select time window end
            do CurrentInstant = Me%HDF5StartInstant, Me%InstantsNumber

                if (AuxInstantsArray(CurrentInstant)                            &
                    .eq. HDF5EndFieldTime) then

                    Me%HDF5EndInstant = (CurrentInstant)
                    HDF5EndFieldTime = HDF5TimeInstant(CurrentInstant)

                    exit

                end if

                if (AuxInstantsArray(CurrentInstant)                            &
                    .gt. HDF5EndFieldTime) then

                    Me%HDF5EndInstant = (CurrentInstant-1)
                    HDF5EndFieldTime = HDF5TimeInstant(CurrentInstant-1)

                    exit 

                end if

            end do

            !Check to see the presence of only one time in the time window
            if (HDF5StartFieldTime .eq. HDF5EndFieldTime) then
                if (Me%StartTime >= HDF5StartTime .AND.                         &
                    Me%EndTime <= HDF5EndTime) then          

                    write(*,*)
                    write(*,*) 'File has only one relevant time:' 
100                 format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0,   &
                            1x, f3.0)
                    call ExtractDate(HDF5StartFieldTime, Year, Month,           & 
                         Day, Hour, Minute, Second)
                    write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
                end if
            end if


        else

            Me%HDF5StartInstant = 1
            Me%HDF5EndInstant = Me%InstantsNumber

            if (Me%HDF5StartInstant == Me%HDF5EndInstant) then
                write(*,*)
                write(*,*) 'HDF5 file provided has only one time.'
                write(*,*) 'HDF5 must have at least 2 times for interval extraction.'
                stop 'OpenAndDateHDF5File - ModuleHDF5Extractor - ERR50'                
            endif                

            StartIntervalTime = AuxInstantsArray(Me%HDF5StartInstant)

            !Get OutputFile name without extension
            AuxOutputName = Me%OutputFile(:(len_trim(Me%OutputFile)-5))

            !Prepare list of output files for extraction by interval
            do while (StartIntervalTime .le. AuxInstantsArray(Me%HDF5EndInstant))

                Me%Interval%NumberIntervalHDF5 = Me%Interval%NumberIntervalHDF5 &
                                                 + 1 

                call AddIntervalHDF5(NewIntervalHDF5)

                !Construct gauge time serie name
                write(AuxChar, fmt=*) Me%Interval%NumberIntervalHDF5
        
                NewIntervalHDF5%Name = trim(AuxOutputName)//                    &
                                       trim(adjustl(AuxChar))//".hdf5"

                StartIntervalTime = StartIntervalTime + Me%Interval%DT

            end do

            if (Me%Interval%NumberIntervalHDF5 == 1) then
                write(*,*)
                write(*,*) 'Interval is longer than provided HDF5 time period.'
                write(*,*) 'Interval must be shorter for interval extraction.'      
                stop 'OpenAndDateHDF5File - ModuleHDF5Extractor - ERR60'                
            endif

        endif

        AuxNumberInstants = Me%HDF5EndInstant - Me%HDF5StartInstant + 1

        allocate(Me%InstantsArray(1:AuxNumberInstants))

        Me%InstantsArray = AuxInstantsArray(Me%HDF5StartInstant: Me%HDF5EndInstant)

        !Deallocate AuxInstantsArray
        deallocate(AuxInstantsArray)
        nullify(AuxInstantsArray)

        call GetHDF5FileID (Me%ObjHDF5File, Me%FileID, STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            & 
            stop 'OpenAndDateHDF5File - ModuleHDF5Extractor - ERR70'

        !Get the native types of HDF5 for data reading and writing
        Me%HDFNativeReal = H5T_NATIVE_REAL
        Me%HDFNativeInteger = H5T_NATIVE_INTEGER
        !(this is needed because these may change with the HDF5 session)

        ! Get reference of the parameters indicated by user to extract
        call InquireParameters

        write(*,*)
        write(*,*) 'Reading time independent data...'

        ! Get all HDF5 reference of time independent/grid itens
        ! (search all file for itens!)
        call InquireSubGroup(Me%FileID, "/", 1)
        

        if (Me%SpatialWindow%XY_ON) then
            if (Me%SpatialWindow%XY_COORD) then
                call ConvertWindow_XY_2_IJ
            endif                
        endif     

        !call KillHDF5(Me%ObjHDF5File)
        
    end subroutine OpenAndDateHDF5File

    !--------------------------------------------------------------------------
    
    subroutine ConvertWindow_XY_2_IJ

        !Local---------------------------------------------------------------------      
        real,    dimension(:,:), pointer :: LonStag, LatStag
        real,    dimension(2,2)          :: WindowLimitsXY = null_real
        integer, dimension(:,:), pointer :: WindowDomain 
        logical                          :: WindowWithData   
        type(T_Item), pointer            :: ObjItem        

        !Begin---------------------------------------------------------------------  
        
        !Point for items independent of time 
        ObjItem => Me%FirstIndependentItem
       
        do while (associated(ObjItem))
                        
            if (trim(ObjItem%Name) == trim("Latitude")) then

                Me%WorkSize%ILB = 1
                Me%WorkSize%IUB = ObjItem%Dimensions(1) - 1
                Me%WorkSize%JLB = 1
                Me%WorkSize%JUB = ObjItem%Dimensions(2) - 1
                
                LatStag => ObjItem%FirstField%Values2D                
                
            endif                
            

            if( trim(ObjItem%Name) == trim("Longitude")) then

                LonStag => ObjItem%FirstField%Values2D
                
            endif        
        
            ObjItem => ObjItem%Next
        enddo
        
        
        nullify(ObjItem)
        

        WindowLimitsXY(2,1) = Me%SpatialWindow%Ymin
        WindowLimitsXY(2,2) = Me%SpatialWindow%Ymax
        WindowLimitsXY(1,1) = Me%SpatialWindow%Xmin
        WindowLimitsXY(1,2) = Me%SpatialWindow%Xmax

        allocate (WindowDomain(2,2))
        
        WindowDomain(:,:) = FillValueInt
    
        call ArrayPolygonWindow(XX              = LonStag,                              &
                                YY              = LatStag,                              &
                                WIn             = WindowLimitsXY,                       &
                                ILB             = Me%WorkSize%ILB,                      &
                                IUB             = Me%WorkSize%IUB+1,                    &
                                JLB             = Me%WorkSize%JLB,                      &
                                JUB             = Me%WorkSize%JUB+1,                    &
                                WOut            = WindowDomain,                         &
                                WindowWithData  = WindowWithData)
                                
        if (WindowWithData) then
        
            Me%SpatialWindow%ILB = max(WindowDomain(1,1),Me%WorkSize%ILB)
            Me%SpatialWindow%IUB = min(WindowDomain(1,2),Me%WorkSize%IUB)
            Me%SpatialWindow%JLB = max(WindowDomain(2,1),Me%WorkSize%JLB)
            Me%SpatialWindow%JUB = min(WindowDomain(2,2),Me%WorkSize%JUB)
            
        else

            Me%SpatialWindow%ILB =Me%WorkSize%ILB
            Me%SpatialWindow%IUB =Me%WorkSize%IUB
            Me%SpatialWindow%JLB =Me%WorkSize%JLB
            Me%SpatialWindow%JUB =Me%WorkSize%JUB

            write(*,*) 'Input file do not intersect the model domain'
            write(*,*) 'The window dimension assumed equal to the original domain'
            write(*,*) 'ConvertWindow_XY_2_IJ - ModuleHDF5Extractor - WRN10'
        
        endif
        
        
        deallocate (WindowDomain)

        nullify(LonStag, LatStag)
        
    end subroutine ConvertWindow_XY_2_IJ
    
    !--------------------------------------------------------------------------    

    type(T_Time) function HDF5TimeInstant(Instant)

        !Arguments-------------------------------------------------------------
        integer                                     :: Instant
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer                 :: TimeVector

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (Me%ObjHDF5File, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            & 
        stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR10'

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = Me%ObjHDF5File,                   &
                             GroupName      = "/"//trim(Me%TimeGroup),          &
                             Name           = trim(Me%TimeGroup),               &
                             Array1D        = TimeVector,                       &
                             OutputNumber   = Instant,                          &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR20'

        call SetDate(HDF5TimeInstant, Year  = TimeVector(1),                    &
                     Month  = TimeVector(2), Day      = TimeVector(3),          &
                     Hour   = TimeVector(4), Minute   = TimeVector(5),          &
                     Second = TimeVector(6))

    end function HDF5TimeInstant

    !--------------------------------------------------------------------------

    subroutine OpenOutputFile(ObjFileName, ObjEnterData)

        !Arguments-------------------------------------------------------------
        character(PathLength)                       :: ObjFileName
        integer                                     :: ObjEnterData
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: MANDATORY_ITEMS
        integer                                     :: HDF5_CREATE
        type(T_Item), pointer                       :: ObjItem

        !Begin-----------------------------------------------------------------

        if (Me%Interval%Interval_ON) then

            Me%Interval%CountIntervalHDF5 = Me%Interval%CountIntervalHDF5 + 1
            
        endif 

        !Create HDF5 file
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (ObjEnterData, trim(ObjFileName),               &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              & 
        stop 'OpenOutputFile - ModuleHDF5Extractor - ERR10'

        !Copy groups for items independent of time 
        ObjItem => Me%FirstIndependentItem

        !Initializes variable to detect mandatory items in HDF5 files
        MANDATORY_ITEMS = 0
        
        do while (associated(ObjItem))

            Me%WorkSize%ILB = 1
            Me%WorkSize%IUB = ObjItem%Dimensions(1)
            Me%WorkSize%JLB = 1
            Me%WorkSize%JUB = ObjItem%Dimensions(2)

            if (Me%SpatialWindow%XY_ON) then

                Me%WorkSize%ILB = Me%SpatialWindow%ILB 
                Me%WorkSize%JLB = Me%SpatialWindow%JLB 
                Me%WorkSize%IUB = Me%SpatialWindow%IUB 
                Me%WorkSize%JUB = Me%SpatialWindow%JUB 

                if( trim(ObjItem%Name) == trim("Latitude"            ) .or.             &
                    trim(ObjItem%Name) == trim("Longitude"           ) .or.             &
                    trim(ObjItem%Name) == trim("ConnectionX"         ) .or.             &
                    trim(ObjItem%Name) == trim("ConnectionY"         ) .or.             &
                    trim(ObjItem%Name) == trim("googlemaps_x"        ) .or.             &
                    trim(ObjItem%Name) == trim("googlemaps_y"       )) then
                                                        
                    Me%WorkSize%IUB = Me%WorkSize%IUB + 1
                    Me%WorkSize%JUB = Me%WorkSize%JUB + 1
                    MANDATORY_ITEMS = MANDATORY_ITEMS + 1
                    
                endif

            endif

            select case (ObjItem%Rank)
                case (1)

                    call HDF5SetLimits (ObjEnterData, Me%WorkSize%ILB,          &   
                                        Me%WorkSize%IUB, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                & 
                    stop 'OpenOutputFile - ModuleHDF5Extractor - ERR20'

                    if (ObjItem%NumType == Me%HDFNativeReal) then

                        call HDF5WriteData (ObjEnterData, ObjItem%Group,        & 
                                            ObjItem%Name, ObjItem%Units,        &
                                            Array1D =  ObjItem%FirstField%Values1D, &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                        stop 'OpenOutputFile - ModuleHDF5Extractor - ERR30'            

                        if (Me%Interval%Interval_ON .and.                       &
                            (Me%Interval%CountIntervalHDF5 ==                   &
                            Me%Interval%NumberIntervalHDF5)) then

                            deallocate(ObjItem%FirstField%Values2D)
                            nullify   (ObjItem%FirstField%Values2D)

                        endif

                    else if (ObjItem%NumType == Me%HDFNativeInteger) then

                        call HDF5WriteData (ObjEnterData, ObjItem%Group,        & 
                                            ObjItem%Name, ObjItem%Units,        &
                                            Array1D =  ObjItem%FirstField%IntValues1D,  &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                        stop 'OpenOutputFile - ModuleHDF5Extractor - ERR40'            

                        if (Me%Interval%Interval_ON .and.                       &
                            (Me%Interval%CountIntervalHDF5 ==                   &
                            Me%Interval%NumberIntervalHDF5)) then

                            deallocate(ObjItem%FirstField%IntValues1D)
                            nullify   (ObjItem%FirstField%IntValues1D)

                        endif

                    endif

            
            
                case (2)

                    call HDF5SetLimits (ObjEnterData, Me%WorkSize%ILB,          &   
                                        Me%WorkSize%IUB, Me%WorkSize%JLB,       &
                                        Me%WorkSize%JUB, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                & 
                    stop 'OpenOutputFile - ModuleHDF5Extractor - ERR20'

                    if (ObjItem%NumType == Me%HDFNativeReal) then

                        call HDF5WriteData (ObjEnterData, ObjItem%Group,        & 
                                            ObjItem%Name, ObjItem%Units,        &
                                            Array2D =  ObjItem%FirstField%Values2D, &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                        stop 'OpenOutputFile - ModuleHDF5Extractor - ERR30'            

                        if (Me%Interval%Interval_ON .and.                       &
                            (Me%Interval%CountIntervalHDF5 ==                   &
                            Me%Interval%NumberIntervalHDF5)) then

                            deallocate(ObjItem%FirstField%Values2D)
                            nullify   (ObjItem%FirstField%Values2D)

                        endif

                    else if (ObjItem%NumType == Me%HDFNativeInteger) then

                        call HDF5WriteData (ObjEnterData, ObjItem%Group,        & 
                                            ObjItem%Name, ObjItem%Units,        &
                                            Array2D =  ObjItem%FirstField%IntValues2D,  &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                        stop 'OpenOutputFile - ModuleHDF5Extractor - ERR40'            

                        if (Me%Interval%Interval_ON .and.                       &
                            (Me%Interval%CountIntervalHDF5 ==                   &
                            Me%Interval%NumberIntervalHDF5)) then

                            deallocate(ObjItem%FirstField%IntValues2D)
                            nullify   (ObjItem%FirstField%IntValues2D)

                        endif

                    endif

                case (3)

                    !calculate Size3D
                    Me%WorkSize%KLB = 1
                    Me%WorkSize%KUB = ObjItem%Dimensions(3)

                    if (Me%SpatialWindow%Layers_ON) then
                        if (Me%SpatialWindow%SurfaceLayer) then
                            Me%WorkSize%KLB = ObjItem%Dimensions(3)
                            Me%WorkSize%KUB = ObjItem%Dimensions(3)
                        else
                            Me%WorkSize%KLB = Me%SpatialWindow%KLB 
                            Me%WorkSize%KUB = Me%SpatialWindow%KUB 
                        endif
                    endif

                    call HDF5SetLimits (ObjEnterData, Me%WorkSize%ILB,          &
                                        Me%WorkSize%IUB, Me%WorkSize%JLB,       &
                                        Me%WorkSize%JUB, Me%WorkSize%KLB,       &
                                        Me%WorkSize%KUB, STAT = STAT_CALL)        
                    if (STAT_CALL .NE. SUCCESS_)                                & 
                    stop 'OpenOutputFile - ModuleHDF5Extractor - ERR50'

                    if (ObjItem%NumType == Me%HDFNativeReal) then

                        call HDF5WriteData (ObjEnterData, ObjItem%Group,        & 
                                            ObjItem%Name, ObjItem%Units,        &
                                            Array3D =  ObjItem%FirstField%Values3D, &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                        stop 'OpenOutputFile - ModuleHDF5Extractor - ERR60'            

                        if (Me%Interval%Interval_ON .and.                       &
                            (Me%Interval%CountIntervalHDF5 ==                   &
                            Me%Interval%NumberIntervalHDF5)) then

                            deallocate(ObjItem%FirstField%Values3D)
                            nullify   (ObjItem%FirstField%Values3D)

                        endif

                    else if (ObjItem%NumType == Me%HDFNativeInteger) then

                        call HDF5WriteData (ObjEnterData, ObjItem%Group,        & 
                                            ObjItem%Name, ObjItem%Units,        &
                                            Array3D =  ObjItem%FirstField%IntValues3D,  &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                        stop 'OpenOutputFile - ModuleHDF5Extractor - ERR70'            

                        if (Me%Interval%Interval_ON .and.                       &
                            (Me%Interval%CountIntervalHDF5 ==                   &
                            Me%Interval%NumberIntervalHDF5)) then

                            deallocate(ObjItem%FirstField%IntValues3D)
                            nullify   (ObjItem%FirstField%IntValues3D)

                        endif

                    endif

            end select

            ObjItem => ObjItem%Next

        end do
        
        !Mandatory items to be detected should be a total of 4 :
        !ConnectionX, ConnectionY or/and Latitude, Longitude
        if (Me%SpatialWindow%XY_ON .and. MANDATORY_ITEMS < 2)                                              &
        stop 'OpenOutputFile - ModuleHDF5Extractor - ERR80'

        !Writes everything to disk
        call HDF5FlushMemory (ObjEnterData, STAT = STAT_CALL)                         
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'OpenOutputFile - ModuleHDF5Extractor - ERR90'

        !(Groups for items dependent of time will be read and written in Modify!)

    end subroutine OpenOutputFile

    !--------------------------------------------------------------------------

    subroutine AddItem (FirstItem, ObjItem)

        !Arguments-------------------------------------------------------------
        type(T_Item), pointer                       :: ObjItem
        type(T_Item), pointer                       :: FirstItem

        !Local-----------------------------------------------------------------
        type(T_Item), pointer                       :: PreviousItem
        type(T_Item), pointer                       :: NewItem

        !Begin-----------------------------------------------------------------

        !Allocates new item
        allocate (NewItem)
        nullify  (NewItem%Next)

        !Insert new item into list and makes current point to it
        if (.not. associated(FirstItem)) then
            FirstItem            => NewItem
            ObjItem              => NewItem
        else
            PreviousItem         => FirstItem
            ObjItem              => FirstItem%Next
            do while (associated(ObjItem))
                PreviousItem     => ObjItem
                ObjItem          => ObjItem%Next
            enddo
            ObjItem              => NewItem
            PreviousItem%Next    => NewItem
        end if

    end subroutine AddItem

    !--------------------------------------------------------------------------

    subroutine ReadItemField (ObjItem, ObjField)

        !Arguments-------------------------------------------------------------
        type(T_Item),   pointer                     :: ObjItem
        type (T_Field), pointer                     :: ObjField

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        select case (ObjItem%Rank)
            case (1)
                !Get 1D values!

                !calculate Size1D
                Me%WorkSize%ILB = 1
                Me%WorkSize%IUB = ObjItem%Dimensions(1)

                call HDF5SetLimits (Me%ObjHDF5File, Me%WorkSize%ILB,            &
                                    Me%WorkSize%IUB, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    & 
                stop 'ReadItemField - ModuleHDF5Extractor - ERR10'

                if (ObjItem%NumType == Me%HDFNativeReal) then
                    !allocate field
                    nullify (ObjField%Values1D)
                    allocate(ObjField%Values1D(Me%WorkSize%ILB:Me%WorkSize%IUB))

                    !read field
                    call HDF5ReadData(Me%ObjHDF5File, ObjItem%Group,            &
                                      trim(ObjField%Name),                      &
                                      Array1D      = ObjField%Values1D,         &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ReadItemField - ModuleHDF5Extractor - ERR20'
                
                elseif (ObjItem%NumType == Me%HDFNativeInteger) then

                    !allocate field
                    nullify (ObjField%IntValues1D)
                    allocate(ObjField%IntValues1D(Me%WorkSize%ILB:Me%WorkSize%IUB))

                    !read field
                    call HDF5ReadData(Me%ObjHDF5File, ObjItem%Group,            &
                                      trim(ObjField%Name),                      &
                                      Array1D      = ObjField%IntValues1D,      &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ReadItemField - ModuleHDF5Extractor - ERR30'
 
                endif


            case (2)
                !Get 2D values!

                !calculate Size2D
                Me%WorkSize%ILB = 1
                Me%WorkSize%IUB = ObjItem%Dimensions(1)
                Me%WorkSize%JLB = 1
                Me%WorkSize%JUB = ObjItem%Dimensions(2)

                call HDF5SetLimits (Me%ObjHDF5File, Me%WorkSize%ILB,            &
                                    Me%WorkSize%IUB, Me%WorkSize%JLB,           &
                                    Me%WorkSize%JUB, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    & 
                stop 'ReadItemField - ModuleHDF5Extractor - ERR40'

                if (ObjItem%NumType == Me%HDFNativeReal) then
                    !allocate field
                    nullify (ObjField%Values2D)
                    allocate(ObjField%Values2D(Me%WorkSize%ILB:Me%WorkSize%IUB, &
                                               Me%WorkSize%JLB:Me%WorkSize%JUB))

                    !read field
                    call HDF5ReadData(Me%ObjHDF5File, ObjItem%Group,            &
                                      trim(ObjField%Name),                      &
                                      Array2D      = ObjField%Values2D,         &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ReadItemField - ModuleHDF5Extractor - ERR50'
                
                elseif (ObjItem%NumType == Me%HDFNativeInteger) then

                    !allocate field
                    nullify (ObjField%IntValues2D)
                    allocate(ObjField%IntValues2D(Me%WorkSize%ILB:Me%WorkSize%IUB,  &
                                                 Me%WorkSize%JLB:Me%WorkSize%JUB))

                    !read field
                    call HDF5ReadData(Me%ObjHDF5File, ObjItem%Group,            &
                                      trim(ObjField%Name),                      &
                                      Array2D      = ObjField%IntValues2D,      &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ReadItemField - ModuleHDF5Extractor - ERR60'
 
                endif

            case (3)
                !Get 3D values!

                !calculate Size3D
                Me%WorkSize%ILB = 1
                Me%WorkSize%IUB = ObjItem%Dimensions(1)
                Me%WorkSize%JLB = 1
                Me%WorkSize%JUB = ObjItem%Dimensions(2)
                Me%WorkSize%KLB = 1
                Me%WorkSize%KUB = ObjItem%Dimensions(3)

                call HDF5SetLimits  (Me%ObjHDF5File, Me%WorkSize%ILB,           &
                                     Me%WorkSize%IUB, Me%WorkSize%JLB,          &
                                     Me%WorkSize%JUB, Me%WorkSize%KLB,          &
                                     Me%WorkSize%KUB, STAT = STAT_CALL)                
                if (STAT_CALL .NE. SUCCESS_)                                    & 
                stop 'ReadItemField - ModuleHDF5Extractor - ERR70'

                if (ObjItem%NumType == Me%HDFNativeReal) then
           
                    !allocate field
                    nullify (ObjField%Values3D)
                    allocate(ObjField%Values3D(Me%WorkSize%ILB:Me%WorkSize%IUB, &
                                              Me%WorkSize%JLB:Me%WorkSize%JUB,  &
                                              Me%WorkSize%KLB:Me%WorkSize%KUB))
            
                    !read field
                    call HDF5ReadData(Me%ObjHDF5File, ObjItem%Group,            &
                                      trim(ObjField%Name),                      &
                                      Array3D      = ObjField%Values3D,         &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ReadItemField - ModuleHDF5Extractor - ERR80'

                elseif (ObjItem%NumType == Me%HDFNativeInteger) then

                    !allocate field
                    nullify (ObjField%IntValues3D)
                    allocate(ObjField%IntValues3D(Me%WorkSize%ILB:Me%WorkSize%IUB,  &
                                              Me%WorkSize%JLB:Me%WorkSize%JUB,  &
                                              Me%WorkSize%KLB:Me%WorkSize%KUB))

                    !read field
                    call HDF5ReadData(Me%ObjHDF5File, ObjItem%Group,            &
                                      trim(ObjField%Name),                      &
                                      Array3D      = ObjField%IntValues3D,      &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ReadItemField - ModuleHDF5Extractor - ERR90'

                endif

        end select

    end subroutine ReadItemField

    !--------------------------------------------------------------------------

    recursive subroutine InquireSubGroup (ID, GroupName, Level)

        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: ID
        character(len=*)                            :: GroupName
        integer                                     :: Level

        !Local-----------------------------------------------------------------
        character(StringLength)                     :: obj_name
        integer                                     :: obj_type
        integer(HID_T)                              :: gr_id, dset_id
        integer(HID_T)                              :: datatype_id, class_id !, size        
        integer                                     :: STAT_CALL
        character(StringLength)                     :: NewGroupName
        integer                                     :: ItensNumber
        integer                                     :: i
        character(StringLength)                     :: Name
        logical                                     :: TimeIndependent = .false.
        type(T_Item), pointer                       :: NewItem

        !Begin-----------------------------------------------------------------

        !Get the number of members in the Group
        call h5gn_members_f(ID, GroupName, ItensNumber, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) return
  
        do i = 1, ItensNumber

            !Gets information about the group
            call h5gget_obj_info_idx_f(ID, GroupName, i-1, obj_name, obj_type,  & 
                                       STAT_CALL)
            if (STAT_CALL /= SUCCESS_) return

            if (GroupName == "/") then

                Me%CurrentGroup = trim(obj_name)

            endif

            if (Me%CurrentGroup == "Residual") then
                ! All itens are time independent but meaningless for extracting purposes
                !(residual values are assessed for the original file period)

                write(*,*)
                write(*,*) 'Residual group is present in file but not extracted.'

            elseif (Me%CurrentGroup == "Statistics") then
                ! All itens are time independent but meaningless for extracting purposes
                !(statistic values are assessed for the original file period)

                write(*,*)
                write(*,*) 'Statistics group is present in file but not extracted.'

                
            elseif (Me%CurrentGroup == "Generic4D") then
                ! All itens are time independent but meaningless for extracting purposes
                !(Generic4D values are assessed for the original file period)

                write(*,*)
                write(*,*) 'Generic4D group is present in file but not extracted.'

            elseif (Me%CurrentGroup == "Grid"     .or. Me%CurrentGroup == "Nodes" .or. &
                    Me%CurrentGroup == "Reaches"  .or. Me%CurrentGroup == "ID") then

                if (obj_type == H5G_DATASET_F) then

                    ! All itens except OpenPoints and VerticalZ are time independent
                    if ( Me%LastSubGroup /= "OpenPoints"  .and.                 &
                         Me%LastSubGroup /= "VerticalZ"   .and.                 &
                         Me%LastSubGroup /= "ScraperPosition") then 

                        !Add to list of time independent itens
                        nullify(NewItem)
                        call AddItem(Me%FirstIndependentItem, NewItem)
            
                        TimeIndependent = .true.

                        NewItem%Name = obj_name
                        
                        NewItem%TimeDependent = .false.

                    else                           

                        !Add to list of time dependent itens
                        nullify(NewItem)
                        call AddItem(Me%FirstDependentItem, NewItem)

                        NewItem%Name = Me%LastSubGroup
                        
                        NewItem%TimeDependent = .true.

                    endif

                    NewItem%Group = GroupName

                    !Get item specifics
                    call GetHDF5GroupID(Me%ObjHDF5File, NewItem%Group, i,       &
                                        Name, Rank = NewItem%Rank,              &
                                        Dimensions = NewItem%Dimensions,        &
                                        Units = NewItem%Units,                  &   
                                        STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                  &
                    stop 'InquireSubGroup - ModuleHDF5Extractor - ERR10'

                    !Get data type (integer or real)
                    !(for time dependent itens assumed that data type equal for all fields)
                    !Opens data set
                    call h5dopen_f(ID, trim(adjustl(obj_name)), dset_id, STAT_CALL)
                    
                    !Gets datatype
                    call h5dget_type_f (dset_id, datatype_id,   STAT_CALL)
                    !call h5tget_size_f (datatype_id, size,      STAT_CALL)
                    call h5tget_class_f(datatype_id, class_id,  STAT_CALL) 

                    if (class_id == H5T_FLOAT_F) then
                        NewItem%NumType = Me%HDFNativeReal
                    elseif (class_id == H5T_INTEGER_F) then
                        NewItem%NumType = Me%HDFNativeInteger
                    else
                        stop 'InquireSubGroup - ModuleHDF5Extractor - ERR20'
                    end if

                    if (TimeIndependent) then

                        nullify(NewItem%FirstField)

                        !Allocates first and only field for this item
                        allocate(NewItem%FirstField)
                        nullify(NewItem%FirstField%Next)

                        NewItem%FirstField%Name = NewItem%Name

                        call ReadItemField(NewItem, NewItem%FirstField) 

                        TimeIndependent = .false.

                    else

                        !do nothing and return to skip the datasets for each time
                        Me%LastSubGroup = obj_name
                        return

                    endif

                 elseif (obj_type == H5G_GROUP_F) then

                    Me%LastSubGroup = obj_name

                    if (GroupName == "/") then
                        NewGroupName = GroupName//trim(adjustl(obj_name))
                    else
                        NewGroupName = GroupName//"/"//trim(adjustl(obj_name))
                    endif
                    call h5gopen_f (ID, trim(adjustl(NewGroupName)), gr_id,     &    
                                    STAT_CALL)
                    call InquireSubGroup (gr_id, trim(adjustl(NewGroupName)),   &
                                          Level + 1)
                    call h5gclose_f (gr_id, STAT_CALL)

                endif

            else
                ! All other groups are assumed to be results groups: time dependent itens
                ! For economy the path to these itens (except group /Time) has to be 
                ! provided by the user and they are read somewhere else 

                ! Code assuming that all itens in these groups are time dependent is 
                ! presented commented here for information purposes only

                !! Itens assumed time dependent
                !! Time is not added to list of dependent itens
                !if (Me%CurrentGroup /= "Time") then
            
                !    nullify(NewItem)
                !    call AddDependentItem(NewItem)
                
                !    !The name is constructed with LastSubGroup name
                !    NewItem%Name = Me%LastSubGroup

                !else

                !    ! Time values are not read; they are treated afterwards
                !    return 

                !endif

                !return

            endif
            
        enddo

    end subroutine InquireSubGroup

    !--------------------------------------------------------------------------

    subroutine InquireParameters

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_Item), pointer                       :: ObjParameter
        integer                                     :: STAT_CALL
        character(StringLength)                     :: Name
        integer                                     :: ItensNumber
        integer(HID_T)                              :: gr_id
        

        !Begin-----------------------------------------------------------------

        !Get info about the parameters in input HDF5 file
        ObjParameter => Me%FirstParameter

        do while(associated(ObjParameter))

            !check if file contains group required
            call h5gopen_f(Me%FileID, trim(adjustl(ObjParameter%Group)), gr_id, & 
                           STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) then  

                write(*,*)
                write(*,*)'HDF5 file: '//trim(Me%HDF5File)
                write(*,*)'Does not contain group required: '                   &
                           //trim(ObjParameter%Group)
                write(*,*)'For parameter: '//trim(ObjParameter%Name)
                stop 'InquireParameters - ModuleHDF5Extractor - ERR10' 
            
            endif 

            !Get parameter number of items
            call GetHDF5GroupNumberOfItems(Me%ObjHDF5File, ObjParameter%Group,  &
                                           ItensNumber, STAT = STAT_CALL)

            !check if the parameter is time dependent
            !(it is assumed time dependent if number of itens is equal to
            !HDF5 file number of instants, previously obtained)
            if (ItensNumber /= Me%InstantsNumber) then

                !This is not a time dependent parameter: extraction is meaningless
                write(*,*)
                write(*,*)'Parameter to extract not time dependent: '           &
                          //trim(ObjParameter%Name)
                ObjParameter%TimeDependent = .false. 
                !stop 'InquireParameters - ModuleHDF5Extractor - ERR20' 
            else
                ObjParameter%TimeDependent = .true.                 
                
            endif

            !Get item specifics
            call GetHDF5GroupID(Me%ObjHDF5File, ObjParameter%Group, 1,          &
                                Name, Rank = ObjParameter%Rank,                 &
                                Dimensions = ObjParameter%Dimensions,           &
                                Units = ObjParameter%Units,                     &   
                                STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'InquireParameters - ModuleHDF5Extractor - ERR30' 

            !Assume real data type
            ObjParameter%NumType = Me%HDFNativeReal
            
            if (Name(1:15)=="ScraperPosition") then
                ObjParameter%NumType = Me%HDFNativeInteger
            endif

            ObjParameter => ObjParameter%Next

        enddo

    end subroutine InquireParameters

    !--------------------------------------------------------------------------

    subroutine AddIntervalHDF5 (ObjHDF5)

        !Arguments-------------------------------------------------------------
        type(T_IntervalHDF5), pointer               :: ObjHDF5

        !Local-----------------------------------------------------------------
        type(T_IntervalHDF5), pointer               :: PreviousHDF5
        type(T_IntervalHDF5), pointer               :: NewHDF5

        !Begin-----------------------------------------------------------------

        !Allocates new item
        allocate (NewHDF5)
        nullify  (NewHDF5%Next)

        !Insert new item into list and makes current point to it
        if (.not. associated(Me%Interval%FirstHDF5)) then
            Me%Interval%FirstHDF5            => NewHDF5
            ObjHDF5              => NewHDF5
        else
            PreviousHDF5         => Me%Interval%FirstHDF5
            ObjHDF5              => Me%Interval%FirstHDF5%Next
            do while (associated(ObjHDF5))
                PreviousHDF5     => ObjHDF5
                ObjHDF5          => ObjHDF5%Next
            enddo
            ObjHDF5              => NewHDF5
            PreviousHDF5%Next    => NewHDF5
        end if

    end subroutine AddIntervalHDF5

    !--------------------------------------------------------------------------

    subroutine PrepareOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_IntervalHDF5), pointer               :: ObjHDF5

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*) 'Writing time independent data...'

        if (Me%Interval%Interval_ON) then

            !cycle all output files, create them and write independent items
            ObjHDF5 => Me%Interval%FirstHDF5

            do while (associated(ObjHDF5))

                call OpenOutputFile(ObjHDF5%Name, ObjHDF5%ObjEnterData)

                ObjHDF5 => ObjHDF5%Next

            enddo

        else

            !create output file and write independent items
            call OpenOutputFile(Me%OutputFile, Me%ObjOutputFile)

        endif

    end subroutine PrepareOutput

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


    !--------------------------------------------------------------------------

    subroutine ModifyExtractHDF5

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_Item), pointer                       :: ObjItem
        type(T_Item), pointer                       :: ObjParameter

        !----------------------------------------------------------------------

        ! Read all period data for each time dependent item and put in memory
        call OpenAndReadHDF5File

        write(*,*)
        write(*,*) 'Writing time dependent data...'

        ! For each time dependent item write fields in the output HDF5 file 
        ObjItem => Me%FirstDependentItem

        !(select the item)
        do while(associated(ObjItem))
        
            if (ObjItem%TimeDependent) then

                call WriteItemFields(ObjItem)

                call KillItemFields(ObjItem)
            
            endif
            
            ObjItem => ObjItem%Next

        end do

        ! For each parameter write fields in the output HDF5 file 
        ObjParameter => Me%FirstParameter

        !(select the item)
        do while(associated(ObjParameter))
        
            if (ObjParameter%TimeDependent) then

                ! Converts V3 in V4 Names and Groups
                if(Me%ConvertV3toV4) then
                    ObjParameter%Name  = ObjParameter%NameV4
                    ObjParameter%Group = ObjParameter%GroupV4
                end if

                call WriteItemFields(ObjParameter)
                call KillItemFields(ObjParameter)
                
            endif                

            ObjParameter => ObjParameter%Next

        end do

    end subroutine ModifyExtractHDF5

    !--------------------------------------------------------------------------

    subroutine  OpenAndReadHDF5File

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_READ
        type(T_Item), pointer                       :: ObjItem
        type(T_Item), pointer                       :: ObjParameter
        integer                                     :: CurrentInstant
        integer                                     :: OutputNumber = 1
        type(T_Time)                                :: StartIntervalTime
        type(T_Time)                                :: NextTime
        type(T_IntervalHDF5), pointer               :: ObjHDF5
        integer                                     :: ObjEnterData

        !Begin-----------------------------------------------------------------
!
!        !(there has been a kill of this file before)
!        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
!
!        !Open HDF5 file
!        call ConstructHDF5 (Me%ObjHDF5File, trim(Me%HDF5File),                  & 
!                            HDF5_READ, STAT = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_) then
!            write(*,*) 'HDF5 file cannot be opened'//Me%HDF5File                
!            stop 'OpenAndReadHDF5File - ModuleHDF5Extractor - ERR10'
!        end if

        !Read fields data
        write(*,*)
        write(*,*)'Reading time dependent data...' 

        if (Me%Interval%Interval_ON) then

            ObjHDF5 => Me%Interval%FirstHDF5

            ObjEnterData = ObjHDF5%ObjEnterData

            !initialize StartIntervalTime
            StartIntervalTime = Me%InstantsArray(1)

        else
            
            ObjEnterData = Me%ObjOutputFile

        endif

        !Time cycle
        do CurrentInstant = Me%HDF5StartInstant, Me%HDF5EndInstant 

            !1. Write time (done here to use the time cycle)
            if (.not. Me%Interval%Interval_ON) then

                call OutputInstants(OutputNumber, OutputNumber, ObjEnterData)

            else

                !write to the correspondent output HDF5 files
                call OutputInstants(CurrentInstant, OutputNumber, ObjEnterData)
                
                if (CurrentInstant /= Me%HDF5EndInstant) then
                    
                    NextTime = Me%InstantsArray(CurrentInstant + 1)

                    !check if next time is beyond the interval
                    if (NextTime > StartIntervalTime + Me%Interval%DT) then

                        !Writes times to disk
                        call HDF5FlushMemory (ObjEnterData, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                              &
                            stop 'OpenAndReadHDF5File - ModuleHDF5Extractor - ERR20'

                        !Set variables to next output HDF5 file
                        ObjHDF5 => ObjHDF5%Next

                        ObjEnterData = ObjHDF5%ObjEnterData

                        StartIntervalTime = StartIntervalTime + Me%Interval%DT

                        OutputNumber = 1
                
                        !write same time to next output HDF5 file
                        call OutputInstants(CurrentInstant, OutputNumber,       &
                                            ObjEnterData)

                    endif

                endif

            endif

            !2. Read time dependent grid items fields
            ObjItem => Me%FirstDependentItem

            do while(associated(ObjItem))

                call ReadTimeDependentFields(ObjItem, CurrentInstant)

                ObjItem => ObjItem%Next

            end do

            !3. Read parameter fields
            ObjParameter => Me%FirstParameter

            do while(associated(ObjParameter))

                call ReadTimeDependentFields(ObjParameter, CurrentInstant)

                ObjParameter => ObjParameter%Next

            end do
            
            OutputNumber = OutputNumber + 1

        end do

        !Writes times to disk
        call HDF5FlushMemory (ObjEnterData, STAT = STAT_CALL)                         
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OpenAndReadHDF5File - ModuleHDF5Extractor - ERR30'

        !call killhdf5(Me%ObjHDF5File)

    end subroutine OpenAndReadHDF5File

    !--------------------------------------------------------------------------

    subroutine ReadTimeDependentFields (ObjItem, CurrentInstant)

        !Arguments-------------------------------------------------------------
        type(T_Item), pointer                       :: ObjItem
        integer                                     :: CurrentInstant

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type (T_Field), pointer                     :: NewField

        !Begin-----------------------------------------------------------------
        
        if (ObjItem%TimeDependent) then

            !Read item fields
            if (CurrentInstant == Me%HDF5StartInstant) then
            
                !(just to inform the user of which items are being actually read)    
                write(*,*)'Reading '//trim(ObjItem%Name)//' fields'
            
            endif

            !construct new fields
            call AddField(ObjItem%FirstField, NewField)

            !get field name
            call GetHDF5GroupID(Me%ObjHDF5File, ObjItem%Group,                      &
                                CurrentInstant, NewField%Name,                      &
                                STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_)                                            &  
            stop 'ReadTimeDependentFields - ModuleHDF5Extractor - ERR10'

            !get field values
            call ReadItemField(ObjItem, NewField)

        endif
    end subroutine ReadTimeDependentFields

    !--------------------------------------------------------------------------

    subroutine AddField (FirstField, ObjField)

        !Arguments-------------------------------------------------------------
        type (T_Field), pointer                     :: FirstField
        type (T_Field), pointer                     :: ObjField
        
        !Local-----------------------------------------------------------------
        type (T_Field), pointer                     :: NewField
        type (T_Field), pointer                     :: PreviousField
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new field
        allocate (NewField)
        nullify  (NewField%Next)

        !Insert new field into list and makes current point to it
        if (.not. associated(FirstField)) then
            nullify(FirstField)
            allocate (FirstField)
            nullify (FirstField%Next)

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

    subroutine OutputInstants(CurrentInstant, OutputNumber, ObjEnterData)

        !Arguments-------------------------------------------------------------
        integer                                     :: CurrentInstant
        integer                                     :: OutputNumber
        integer                                     :: ObjEnterData

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        real,    dimension(6), target               :: AuxTime
        real,    dimension(:), pointer              :: TimePtr
        type(T_Time)                                :: CurrentDate

        !Begin-----------------------------------------------------------------
        
        CurrentDate = Me%InstantsArray(CurrentInstant)

        call ExtractDate   (CurrentDate,                                        &
                            AuxTime(1), AuxTime(2), AuxTime(3),                 &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime

        call HDF5SetLimits  (ObjEnterData, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputInstants - ModuleHDF5Extractor - ERR10'


        call HDF5WriteData  (ObjEnterData, "/Time",                             &
                             "Time", "YYYY/MM/DD HH:MM:SS",                     &
                             Array1D = TimePtr,                                 &
                             OutputNumber = OutputNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputInstants - ModuleHDF5Extractor - ERR20'


    end subroutine OutputInstants

    !--------------------------------------------------------------------------

    subroutine WriteItemFields (ObjItem)

        !Arguments-------------------------------------------------------------
        type(T_Item), pointer                       :: ObjItem

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: OutputNumber, InstantNumber
        type(T_Field), pointer                      :: CurrentField
        integer                                     :: ObjEnterData
        type(T_IntervalHDF5), pointer               :: ObjOutputFile
        type(T_Time)                                :: StartIntervalTime
        type(T_Time)                                :: NextTime

        !Begin-----------------------------------------------------------------

        !calculate Size2D
        Me%WorkSize%ILB = 1
        Me%WorkSize%IUB = ObjItem%Dimensions(1)
        Me%WorkSize%JLB = 1
        Me%WorkSize%JUB = ObjItem%Dimensions(2)

        if (Me%SpatialWindow%XY_ON) then

            Me%WorkSize%ILB = Me%SpatialWindow%ILB 
            Me%WorkSize%JLB = Me%SpatialWindow%JLB 
            Me%WorkSize%IUB = Me%SpatialWindow%IUB 
            Me%WorkSize%JUB = Me%SpatialWindow%JUB 

            if( trim(ObjItem%Name) == trim("Latitude"            ) .or.                 &
                trim(ObjItem%Name) == trim("Longitude"           ) .or.                 &
                trim(ObjItem%Name) == trim("ConnectionX"         ) .or.                 &
                trim(ObjItem%Name) == trim("ConnectionY"         ) .or.                 &
                trim(ObjItem%Name) == trim("Spherical_Mercator_X") .or.                 &
                trim(ObjItem%Name) == trim("Spherical_Mercator_Y")) then
                                        
                Me%WorkSize%IUB = Me%WorkSize%IUB + 1
                Me%WorkSize%JUB = Me%WorkSize%JUB + 1
                    
            endif

        endif
        
        if (Me%Interval%Interval_ON) then

            ObjOutputFile => Me%Interval%FirstHDF5
            ObjEnterData = ObjOutputFile%ObjEnterData
            StartIntervalTime = Me%InstantsArray(1)

        else

            ObjEnterData = Me%ObjOutputFile

        endif
        
        select case (ObjItem%Rank)
            case (1)

                call HDF5SetLimits (ObjEnterData, Me%WorkSize%ILB,              &
                                    Me%WorkSize%IUB, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    & 
                stop 'WriteItemFields - ModuleHDF5Extractor - ERR10'

                write(*,*)'Writing '//trim(ObjItem%Name)//' fields'

                CurrentField => ObjItem%FirstField

                OutputNumber = 1

                InstantNumber = 1

                !(select the field)
do1:            do while(associated(CurrentField))

                    !(write field according to number type)
                    if (ObjItem%NumType == Me%HDFNativeReal) then

                        call HDF5WriteData(ObjEnterData, ObjItem%Group,         & 
                                           ObjItem%Name,                        & 
                                           ObjItem%Units,                       &
                                           Array1D = CurrentField%Values1D,     &
                                           OutputNumber = OutputNumber,         &
                                           STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                            stop 'OutputField - ModuleHDF5Extractor - ERR20'            

                    else if (ObjItem%NumType == Me%HDFNativeInteger) then

                        call HDF5WriteData(ObjEnterData, ObjItem%Group,         & 
                                           ObjItem%Name,                        &
                                           ObjItem%Units,                       &
                                           Array1D = CurrentField%IntValues1D,  &
                                           OutputNumber = OutputNumber,         &
                                           STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                            stop 'WriteItemFields - ModuleHDF5Extractor - ERR30'            

                    endif

                    if (Me%Interval%Interval_ON) then

                        if (InstantNumber < Me%InstantsNumber) then

                            NextTime = Me%InstantsArray(InstantNumber + 1)

                            if (NextTime > StartIntervalTime + Me%Interval%DT) then

                                !Writes everything to disk
                                call HDF5FlushMemory (ObjEnterData, STAT = STAT_CALL)                         
                                if (STAT_CALL /= SUCCESS_)                      &
                                    stop 'WriteItemFields - ModuleHDF5Extractor - ERR40'

                                !A new HDF5 file must be writen 
                                ObjOutputFile => ObjOutputFile%Next

                                ObjEnterData = ObjOutputFile%ObjEnterData
                        
                                StartIntervalTime = StartIntervalTime + Me%Interval%DT

                                OutputNumber = 1

                                call HDF5SetLimits (ObjEnterData, Me%WorkSize%ILB, &
                                    Me%WorkSize%IUB, STAT = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_)                    & 
                                    stop 'WriteItemFields - ModuleHDF5Extractor - ERR50'

                            else
                    
                                CurrentField => CurrentField%Next

                                OutputNumber = OutputNumber + 1

                                InstantNumber = InstantNumber + 1

                            endif

                        else 

                            exit do1

                        endif
                    
                    else

                        CurrentField => CurrentField%Next

                        OutputNumber = OutputNumber + 1

                    endif


                end do do1


            case (2)

                call HDF5SetLimits (ObjEnterData, Me%WorkSize%ILB,              &
                                    Me%WorkSize%IUB, Me%WorkSize%JLB,           &
                                    Me%WorkSize%JUB, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    & 
                stop 'WriteItemFields - ModuleHDF5Extractor - ERR10'

                write(*,*)'Writing '//trim(ObjItem%Name)//' fields'

                CurrentField => ObjItem%FirstField

                OutputNumber = 1

                InstantNumber = 1

                !(select the field)
do2:            do while(associated(CurrentField))

                    !(write field according to number type)
                    if (ObjItem%NumType == Me%HDFNativeReal) then

                        call HDF5WriteData(ObjEnterData, ObjItem%Group,         & 
                                           ObjItem%Name,                        & 
                                           ObjItem%Units,                       &
                                           Array2D = CurrentField%Values2D,     &
                                           OutputNumber = OutputNumber,         &
                                           STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                            stop 'OutputField - ModuleHDF5Extractor - ERR20'            

                    else if (ObjItem%NumType == Me%HDFNativeInteger) then

                        call HDF5WriteData(ObjEnterData, ObjItem%Group,         & 
                                           ObjItem%Name,                        &
                                           ObjItem%Units,                       &
                                           Array2D = CurrentField%IntValues2D,  &
                                           OutputNumber = OutputNumber,         &
                                           STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                            stop 'WriteItemFields - ModuleHDF5Extractor - ERR30'            

                    endif

                    if (Me%Interval%Interval_ON) then

                        if (InstantNumber < Me%InstantsNumber) then

                            NextTime = Me%InstantsArray(InstantNumber + 1)

                            if (NextTime > StartIntervalTime + Me%Interval%DT) then

                                !Writes everything to disk
                                call HDF5FlushMemory (ObjEnterData, STAT = STAT_CALL)                         
                                if (STAT_CALL /= SUCCESS_)                      &
                                    stop 'WriteItemFields - ModuleHDF5Extractor - ERR40'

                                !A new HDF5 file must be writen 
                                ObjOutputFile => ObjOutputFile%Next

                                ObjEnterData = ObjOutputFile%ObjEnterData
                        
                                StartIntervalTime = StartIntervalTime + Me%Interval%DT

                                OutputNumber = 1

                                call HDF5SetLimits (ObjEnterData, Me%WorkSize%ILB, &
                                    Me%WorkSize%IUB, Me%WorkSize%JLB,           &
                                    Me%WorkSize%JUB, STAT = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_)                    & 
                                    stop 'WriteItemFields - ModuleHDF5Extractor - ERR50'

                            else
                    
                                CurrentField => CurrentField%Next

                                OutputNumber = OutputNumber + 1

                                InstantNumber = InstantNumber + 1

                            endif

                        else 

                            exit do2

                        endif
                    
                    else

                        CurrentField => CurrentField%Next

                        OutputNumber = OutputNumber + 1

                    endif


                end do do2
            
            case (3)

                !calculate Size3D
                Me%WorkSize%KLB = 1
                Me%WorkSize%KUB = ObjItem%Dimensions(3)

                if (Me%SpatialWindow%Layers_ON) then
                
                        if (Me%SpatialWindow%SurfaceLayer) then

                            Me%WorkSize%KUB = ObjItem%Dimensions(3)
                            Me%WorkSize%KLB = Me%WorkSize%KUB
                            
                            if (trim(ObjItem%Name)=="VerticalZ") then
                                Me%WorkSize%KLB = Me%WorkSize%KUB - 1
                            endif  
                            
                        else                
                            Me%WorkSize%KLB = Me%SpatialWindow%KLB 
                            Me%WorkSize%KUB = Me%SpatialWindow%KUB 

                            if (trim(ObjItem%Name)=="VerticalZ") then
                                Me%WorkSize%KUB = Me%WorkSize%KUB + 1
                            endif                          

                        endif
 
                endif

                if (trim(ObjItem%Name)=="VerticalZ") then
                    ObjItem%Name = "Vertical"
                endif

                call HDF5SetLimits (ObjEnterData, Me%WorkSize%ILB,              &
                                    Me%WorkSize%IUB, Me%WorkSize%JLB,           &
                                    Me%WorkSize%JUB, Me%WorkSize%KLB,           &
                                    Me%WorkSize%KUB, STAT = STAT_CALL)        
                if (STAT_CALL .NE. SUCCESS_)                                    & 
                stop 'WriteItemFields - ModuleHDF5Extractor - ERR60'

                write(*,*)'Writing '//trim(ObjItem%Name)//' fields'
                


                CurrentField => ObjItem%FirstField

                OutputNumber = 1

                InstantNumber = 1

                !(select the field)
do3:            do while(associated(CurrentField))

                    !(write field according to number type)
                    if (ObjItem%NumType == Me%HDFNativeReal) then

                        call HDF5WriteData(ObjEnterData, ObjItem%Group,         & 
                                           ObjItem%Name,                        &
                                           ObjItem%Units,                       &
                                           Array3D = CurrentField%Values3D,     &
                                           OutputNumber = OutputNumber,         &
                                           STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                            stop 'WriteItemFields - ModuleHDF5Extractor - ERR70'            

                    else if (ObjItem%NumType == Me%HDFNativeInteger) then

                        call HDF5WriteData(ObjEnterData, ObjItem%Group,         & 
                                           ObjItem%Name,                        &
                                           ObjItem%Units,                       &
                                           Array3D =  CurrentField%IntValues3D, &
                                           OutputNumber = OutputNumber,         &
                                           STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                            stop 'WriteItemFields - ModuleHDF5Extractor - ERR80'            

                    endif

                    if (Me%Interval%Interval_ON) then

                        if (InstantNumber < Me%InstantsNumber) then

                            NextTime = Me%InstantsArray(InstantNumber+1)

                            if (NextTime > StartIntervalTime +                  &
                                Me%Interval%DT) then

                                !Writes everything to disk
                                call HDF5FlushMemory (ObjEnterData, STAT = STAT_CALL)                         
                                if (STAT_CALL /= SUCCESS_)                      &
                                    stop 'WriteItemFields - ModuleHDF5Extractor - ERR90'

                                !A new HDF5 file must be writen 
                                ObjOutputFile => ObjOutputFile%Next

                                ObjEnterData = ObjOutputFile%ObjEnterData
                        
                                StartIntervalTime = StartIntervalTime + Me%Interval%DT

                                OutputNumber = 1

                                call HDF5SetLimits (ObjEnterData,               &
                                                    Me%WorkSize%ILB,            &
                                                    Me%WorkSize%IUB,            &
                                                    Me%WorkSize%JLB,            &
                                                    Me%WorkSize%JUB,            &
                                                    Me%WorkSize%KLB,            &
                                                    Me%WorkSize%KUB, STAT = STAT_CALL)        
                                if (STAT_CALL .NE. SUCCESS_)                    & 
                                stop 'WriteItemFields - ModuleHDF5Extractor - ERR100'

                            else
                    
                                CurrentField => CurrentField%Next

                                OutputNumber = OutputNumber + 1

                                InstantNumber = InstantNumber + 1

                            endif

                        else 

                            exit do3

                        endif
                    
                    else

                        CurrentField => CurrentField%Next

                        OutputNumber = OutputNumber + 1

                    endif

                enddo do3

        end select

        !Writes everything to disk
        call HDF5FlushMemory (ObjEnterData, STAT = STAT_CALL)                         
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'WriteItemFields - ModuleHDF5Extractor - ERR110'

    end subroutine WriteItemFields

    !--------------------------------------------------------------------------



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    !--------------------------------------------------------------------------

    subroutine KillExtractHDF5

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_Item), pointer                       :: ItemX, ItemToKill 
        type(T_IntervalHDF5), pointer               :: HDF5FileX, HDF5FileToKill

        !------------------------------------------------------------------------
        
        call KillHDF5(Me%ObjHDF5File)

        !Deallocate the InstantsArray
        if (associated(Me%InstantsArray)) then 
            deallocate(Me%InstantsArray, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'KillExtractHDF5 - ModuleHDF5Extractor - ERR10'               
            nullify(Me%InstantsArray)
        end if

        if (Me%Interval%Interval_ON) then
            !Deallocate interval HDF5 files
            HDF5FileX => Me%Interval%FirstHDF5

            do while(associated(HDF5FileX))  

                HDF5FileToKill  => HDF5FileX
                HDF5FileX       => HDF5FileX%Next
                call killhdf5(HDF5FileToKill%ObjEnterData)
                deallocate(HDF5FileToKill, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                stop 'KillExtractHDF5 - ModuleInvertedBarometer - ERR20'
                nullify(HDF5FileToKill)

            end do

        else

            call killhdf5(Me%ObjOutputFile)

        endif

        !Kill all items from the several lists
        ItemX=> Me%FirstParameter

        do while(associated(ItemX))  

            ItemToKill => ItemX
            ItemX      => ItemX%Next
            call KillIndividualItem(ItemToKill)

        end do
        nullify(Me%FirstParameter)

        if (associated(Me%FirstDependentItem)) then
            ItemX=> Me%FirstDependentItem

            do while(associated(ItemX))  

                ItemToKill => ItemX
                ItemX      => ItemX%Next
                call KillIndividualItem(ItemToKill)

            end do
            nullify(Me%FirstDependentItem)
        endif

        ItemX=> Me%FirstIndependentItem

        do while(associated(ItemX))  

            ItemToKill => ItemX
            ItemX      => ItemX%Next
            call KillIndividualItem(ItemToKill)

        end do
        nullify(Me%FirstIndependentItem)

        deallocate(Me)
        nullify(Me)

        !------------------------------------------------------------------------

    end subroutine KillExtractHDF5
        
    !--------------------------------------------------------------------------

    subroutine KillItemFields(ItemToDispose)

        !Arguments-------------------------------------------------------------
        type(T_Item), pointer                       :: ItemToDispose

        !Local-----------------------------------------------------------------
        type(T_Field), pointer                      :: FieldToKill, CurrentField

        !Begin-----------------------------------------------------------------

        CurrentField => ItemToDispose%FirstField

        do while(associated(CurrentField))

            FieldToKill => CurrentField
            CurrentField => CurrentField%Next

            call KillIndividualField(FieldToKill)

        end do 

    end subroutine KillItemFields

    !--------------------------------------------------------------------------

    subroutine KillIndividualItem(ItemToDispose)

        !Arguments-------------------------------------------------------------
        type(T_Item), pointer                       :: ItemToDispose

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !The fields of the parameters have already been killed
        deallocate(ItemToDispose, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
        stop 'KillIndividualItem - ModuleHDF5Extractor - ERR10'
        nullify(ItemToDispose)

    end subroutine KillIndividualItem

    !--------------------------------------------------------------------------  

    subroutine KillIndividualField(FieldToDispose)
    
        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                      :: FieldToDispose

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (associated(FieldToDispose%Values2D)) then
            deallocate(FieldToDispose%Values2D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'KillIndividualField - ModuleHDF5Extractor - ERR10'  
            nullify(FieldToDispose%Values2D)
        end if

        if (associated(FieldToDispose%Values3D)) then
            deallocate(FieldToDispose%Values3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'KillIndividualField - ModuleHDF5Extractor - ERR20'  
            nullify(FieldToDispose%Values3D)
        end if

        if (associated(FieldToDispose%IntValues2D)) then
            deallocate(FieldToDispose%IntValues2D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'KillIndividualField - ModuleHDF5Extractor - ERR30'  
            nullify(FieldToDispose%IntValues2D)
        end if
        
        if (associated(FieldToDispose%IntValues3D)) then
            deallocate(FieldToDispose%IntValues3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'KillIndividualField - ModuleHDF5Extractor - ERR40'  
            nullify(FieldToDispose%IntValues3D)
        end if    

    end subroutine KillIndividualField

    !--------------------------------------------------------------------------  


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

end module ModuleHDF5Extractor









