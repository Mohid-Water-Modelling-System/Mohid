!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid ConvertToHDF5
! MODULE        : UpscaleHDF5
! URL           : http://www.mohid.com
! AFFILIATION   : ColabAtlantic
! DATE          : April 2022
! REVISION      : Joao Sobrinho - v1.0
! DESCRIPTION   : Module to fill grid a parent domain with several HDF5 files with specific grids/areas through upscaling
!
!------------------------------------------------------------------------------
!DataFile
!
!   <begin_file>
!
!   ACTION                      : char                  [-]     !'UPSCALE HDF5 FILES' to use this module
!   UPSCALE_METHOD              : int                   [-]     !Type of upscaling to use. 1 = Average
!                                                               !(only implemented 3 - Triangulation)
!
!   OUTPUTFILENAME              : char                  [-]     !Path to HDF5 file to be created
!
!   <<begin_child>>                                            
!   CHILD_FILENAME              : char                  [-]     !Path to father HDF5 file
!   CHILD_GRID_FILENAME         : char                  [-]     !Path to father grid file
!   <<end_child>> 
!   <<BeginFields>>
!   temperature                 : char                  [-]     ! name of the HDF5 field to be upscaled
!   <<EndFields>>
!   <end_file>                                             

Module ModuleUpscaleHDF5

    use ModuleGlobalData
    use ModuleTime
    use ModuleDrawing
    use ModuleHDF5
    use ModuleEnterData,        only : GetData, ExtractBlockFromBlock, Block_Unlock,    &
                                       GetBlockSize, RewindBlock
    use ModuleHorizontalGrid  
    use ModuleGridData,         only : ConstructGridData, GetGridData, UngetGridData,   &
                                       KillGridData 
    use ModuleHorizontalMap,    only : ConstructHorizontalMap, GetWaterPoints2D,        &
                                       UngetHorizontalMap, KillHorizontalMap
    use ModuleTriangulation,    only : ConstructTriangulation, SetHeightValues,         & 
                                       InterPolation, KillTriangulation
    use ModuleGeometry
    use ModuleMap

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartUpscaleHDF5
    private ::      ReadOptions
    private ::          ReadChildKeywords
    private ::              AddChildFile
    private ::                  ConstructChildFile
    private ::                  CopyChildFile
    private ::      ConstructChildGrid
    private ::      ConstructGridPolygon
    private ::          SetLimitsPolygon
    private ::      ConstructNewGrid
    private ::      set_module_time
    !private ::      ConstructRelevantChild
    private ::      FatherSonCommunication
    private ::      Open_HDF5_OutPut_File
    private ::      OpenAndReadChildHDF5Files
    private ::          get_timeinstants
    private ::          get_nproperties
    private ::          HDF5TimeInstant
    private ::          OutputFields
    private ::          OutputFields3D
    private ::      KillUpscaleHDF5
    private ::          KillChildGrid

    !Selector                    
    
    !Modifier

    !Destructor

    !Management
    
    !Interfaces----------------------------------------------------------------

    !Parameters----------------------------------------------------------------

    !Types---------------------------------------------------------------------
    type       T_Field
        character(len=StringLength)                         :: Name
        character(len=StringLength)                         :: Units
        integer                                             :: IDNumber
        type(T_Time)                                        :: Date
        real, dimension(:,:  ),     pointer                 :: Values2D
        real, dimension(:,:,:),     pointer                 :: Values3D
        type(T_Field),              pointer                 :: Next
    end type  T_Field

    type       T_Grid
        integer                                             :: ObjHorizontalGrid    = 0
        integer                                             :: ObjHDF5              = 0
        integer                                             :: ObjHorizontalMap     = 0
        integer                                             :: ObjBathymetry        = 0
        integer                                             :: ObjGeometry          = 0
        integer                                             :: ObjMap               = 0
        character(len=StringLength)                         :: Name
        character(len=StringLength)                         :: Units
        character(len=PathLength)                           :: OutputFileName
        character(len=PathLength)                           :: InputFileName, InputFileGrid
        character(len=PathLength)                           :: GridFileName
        character(len=PathLength)                           :: GeometryFileName
        logical                                             :: Batim_from_file = .false.
        logical                                             :: GridFromHDF5 = .false.
        type(T_Time), dimension(:), pointer                 :: InstantsArray
        integer                                             :: NumberOfInstants
        integer                                             :: NumberOfProperties
        type(T_Size2D)                                      :: Size2D
        type(T_Size3D)                                      :: Size3D
        type(T_Size2D)                                      :: WorkSize2D
        type(T_Size3D)                                      :: WorkSize3D
        integer, dimension(:,:),    pointer                 :: WaterPoints2D
        integer, dimension(:,:, :), pointer                 :: WaterPoints3D
        real,    dimension(:,:),    pointer                 :: ConnectionX, ConnectionY
        integer                                             :: Level
        type(T_Polygon), pointer                            :: Polygon
        integer                                             :: NRemoveFrame
        integer, dimension(:), allocatable                  :: FatherPoint_ID_I, FatherPoint_ID_J
        real,    dimension(:  ), allocatable                :: NodeX, NodeY, NodeZ
        !integer, dimension(:), allocatable                  :: State_open
        integer, dimension(:,:), allocatable                :: nPointsInside
        real, dimension(:,:), pointer                       :: Bathym
    end type  T_Grid

    type       T_Child
        type(T_Grid)                                        :: Info
        integer                                             :: NumberOfNodes
        type(T_Field), pointer                              :: NewField
        type(T_Child), pointer                              :: Next
    end type  T_Child
    
    private :: T_UpscaleHDF5
    type       T_UpscaleHDF5
        integer                                             :: ObjEnterData         = 0
        integer                                             :: ObjTime              = 0
        integer                                             :: upscale_method = 1
        logical                                             :: Upscale3D = .false.
        type(T_Grid )                                       :: New, Father
        type(T_Time)                                        :: BeginTime, EndTime
        logical                                             :: TimeWindow = .true.
        real,    dimension(:,:), pointer                    :: SonCenterX, SonCenterY
        type(T_Child), pointer                              :: FirstChildFile
        character(len=StringLength), dimension(:), pointer  :: FieldsToUpscale
        integer                                             :: nFieldsToUpscale
        logical                                             :: ConvertAllFields
    end type  T_UpscaleHDF5

    !Global Module Variables
    type (T_UpscaleHDF5), pointer                        :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartUpscaleHDF5(ObjUpscaleHDF5ID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: ObjUpscaleHDF5ID
        integer,           intent(IN )                  :: ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !External----------------------------------------------------------------
        
        
        !Local-------------------------------------------------------------------
        type (T_Child),     pointer                    :: ObjChildFile

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, ObjUpscaleHDF5ID)

        !Read instructions from user
        call ReadOptions(ClientNumber)

        !Cycle all the childs to construct grid
        ObjChildFile => Me%FirstChildFile

        do while (associated(ObjChildFile))
            call ConstructChildGrid(ObjChildFile)

            call ConstructGridPolygon(ObjChildFile%Info)

            ObjChildFile => ObjChildFile%Next
        enddo

        call ConstructNewGrid
 
        !!Construct array of relevant child files (one for each cell)
        !call ConstructRelevantChild 
        !
        !Construct relation between the fathers and the new grid
        call FatherSonCommunication
        
        !Open Output File (construct header data)
        if(Me%Upscale3D) then
            call Open_HDF5_OutPut_File3D
        else 
            call Open_HDF5_OutPut_File
        endif

        !Read data from father files
        call OpenAndReadChildHDF5Files

        call KillUpscaleHDF5

        STAT = SUCCESS_

        !----------------------------------------------------------------------

    end subroutine StartUpscaleHDF5

    !------------------------------------------------------------------------

    subroutine ReadOptions(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer, intent(in)                         :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------
        
        write(*,*)'Reading instructions...'

        call GetData(Me%New%OutputFileName,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUTFILENAME',                   &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)'Error reading OUTPUTFILENAME keyword'
            stop 'ReadOptions - ModuleUpscaleHDF5 - ERR10'
        end if
        
        if (iflag == 0)then
            write(*,*)'Must specify name of output file'
            stop 'ReadOptions - ModuleUpscaleHDF5 - ERR20'
        end if
        call GetData(Me%Father%InputFileName,                           &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'INPUTFILENAME',                    &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)'Error reading INPUTFILENAME keyword'
            stop 'ReadOptions - ModuleUpscaleHDF5 - ERR30'
        end if
        
        if (iflag == 0)then
            write(*,*)'Must specify name of output file'
            stop 'ReadOptions - ModuleUpscaleHDF5 - ERR40'
        end if
        
        call GetData(Me%upscale_method,                                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'UPSCALE_METHOD',                   &  
                     Default      = 1,                                  &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)'Error reading UPSCALE_METHOD keyword. Check if you are using an integer'
            stop 'ReadOptions - ModuleUpscaleHDF5 - ERR50'
        end if
        
        if (iflag == 0)then
            write(*,*)'Must specify the upscaling methodology to use'
            stop 'ReadOptions - ModuleUpscaleHDF5 - ERR60'
        end if
        
        call GetData(Me%Upscale3D,                                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'UPSCALE_3D',                       &
                     ClientModule = 'ConvertToHDF5',                    &
                     Default      = .false.,                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleUpscaleHDF5 - ERR70'
        
        call GetData(Me%Father%GeometryFileName,                           &
                        Me%ObjEnterData, iflag,                            &
                        SearchType   = FromBlock,                          &
                        keyword      = 'FATHER_GEOMETRY',                  &
                        ClientModule = 'ConvertToHDF5', STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)'Error reading FATHER_GEOMETRY keyword'
            stop 'ReadOptions - ModuleUpscaleHDF5 - ERR80'
        end if
        
        if (Me%Upscale3D) then
            if (iflag == 0)then
                write(*,*)'Must specify name of FATHER_GEOMETRY file'
                stop 'ReadOptions - ModuleUpscaleHDF5 - ERR90'
            end if
        end if
        
        call GetData(Me%Father%Batim_from_file,                            &
                        Me%ObjEnterData, iflag,                            &
                        SearchType   = FromBlock,                          &
                        keyword      = 'FATHERBATIM_FROM_FILE',            &
                        ClientModule = 'ConvertToHDF5', STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)'Error reading FATHERBATIM_FROM_FILE keyword'
            stop 'ReadOptions - ModuleUpscaleHDF5 - ERR100'
        end if
        
        Me%Father%GridFromHDF5 = .true.
        
        if (.not. Me%Father%Batim_from_file) then
            call GetData(Me%Father%InputFileGrid,                           &
                         Me%ObjEnterData, iflag,                            &
                         SearchType   = FromBlock,                          &
                         keyword      = 'INPUTFILEGRID',                    &
                         ClientModule = 'ConvertToHDF5',                    &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) then
                write(*,*)'Error reading INPUTFILEGRID keyword'
                stop 'ReadOptions - ModuleUpscaleHDF5 - ERR110'
            end if
            if (iflag == 0)then
                write(*,*)'Must specify name of INPUTFILEGRID file'
                stop 'ReadOptions - ModuleUpscaleHDF5 - ERR120'
            end if
            Me%Father%GridFromHDF5 = .false.
        end if
        
        call ReadChildKeywords (ClientNumber) 

        call ReadFieldsToConvert(ClientNumber)

    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    subroutine ReadFieldsToConvert(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer, intent(in)                 :: ClientNumber
        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL, StartLine, EndLine, Count, CurrentLineNumber
        logical                             :: BlockFound
        character(len=StringLength)         :: PropertyName

        !Begin-----------------------------------------------------------------
        Me%ConvertAllFields = .true. 
        call RewindBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModuleUpscaleHDF5 - ERR00'

        call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,             &
                                    '<<BeginFields>>', '<<EndFields>>',       &
                                    BlockFound, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)'Problem reading <<BeginFields>> / <<EndFields>> block'
            stop 'ReadFieldsToConvert - ModuleUpscaleHDF5 - ERR10'
        end if
        
        if(BlockFound)then
            Me%ConvertAllFields = .false.
    
            call GetBlockSize(Me%ObjEnterData, ClientNumber, StartLine,EndLine, &
                              FromBlockInBlock, STAT = STAT_CALL)

            Me%nFieldsToUpscale = EndLine - StartLine - 1

            allocate(Me%FieldsToUpscale(1:Me%nFieldsToUpscale))
        
            Count = 1

            do CurrentLineNumber = StartLine + 1 , EndLine - 1
                call GetData(PropertyName, Me%ObjEnterData, flag, SearchType  = FromBlock_, &
                             Buffer_Line = CurrentLineNumber, STAT = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) then
                    write(*,*)'Problem reading property name given in block <<BeginFields>> / <<EndFields>>'
                    stop 'ReadFieldsToConvert - ModuleUpscaleHDF5 - ERR20'
                end if

                Me%FieldsToUpscale(Count) = PropertyName

                Count = Count + 1
            end do

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModuleUpscaleHDF5 - ERR30'
        endif

    end subroutine ReadFieldsToConvert

    !--------------------------------------------------------------------------

    subroutine ReadChildKeywords(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer                                     :: ClientNumber
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: BlockFound
        logical                                     :: AtLeastOneBlock = .false.
        !Begin-----------------------------------------------------------------
        do
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,           &
                                        '<<begin_child>>','<<end_child>>',      &
                                        BlockFound, STAT = STAT_CALL)
           if(STAT_CALL .EQ. SUCCESS_)then
                if (BlockFound) then                                                  
                    AtLeastOneBlock = .true.
                    call AddChildFile
                else
                    if (STAT_CALL .NE. SUCCESS_) then
                        write(*,*) 'Could not find any <<begin_child>>/<<end_child>> block. '
                        stop 'ReadChildKeywords - ModuleUpscaleHDF5 - ERR10'
                    end if
                    exit
                end if
            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then
                write(*,*) 'Error calling ExtractBlockFromBuffer on <<begin_child>>/<<end_child>> blocks. '
                stop 'ReadChildKeywords - ModuleUpscaleHDF5 - ERR20'
            else
                write(*,*) 'Problem reading a <<begin_child>>/<<end_child>> block. make sure you have at least one '
                stop 'ReadChildKeywords - ModuleUpscaleHDF5 - ERR30'
            end if
        end do

        if (.not. AtLeastOneBlock) then                                            
            write(*,*) 'No child file block is indicated in input file. '
            stop 'ReadChildKeywords - ModuleUpscaleHDF5 - ERR40'
        end if

    end subroutine ReadChildKeywords
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ConstructChildFile(NewChildFile)

        !Arguments-------------------------------------------------------------
        type (T_Child),      pointer           :: NewChildFile
        !External--------------------------------------------------------------
        integer                                   :: iflag, STAT_CALL
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        call GetData(NewChildFile%Info%InputFileName,                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'CHILD_FILENAME',                   &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)'Could not read keyword CHILD_FILENAME, check if you are writing a path for the file'
            stop 'ConstructChildFile - ModuleUpscaleHDF5 - ERR10'
        end if
        if (iflag == 0)then
            write(*,*)'Must specify name of child file to convert (CHILD_FILENAME)'
            stop 'ConstructChildFile - ModuleUpscaleHDF5 - ERR20'
        end if

        call GetData(NewChildFile%Info%GridFileName,                    &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'CHILD_GRID_FILENAME',              &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)'Could not read keyword CHILD_GRID_FILENAME'
            stop 'ConstructChildFile - ModuleUpscaleHDF5 - ERR30'
        end if
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ConstructChildFile - ModuleUpscaleHDF5 - ERR40'
        end if
        
        call GetData(NewChildFile%Info%NRemoveFrame,                    &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'N_REMOVE_FRAME',                   &
                     default      = 0,                                  &                        
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructChildFile - ModuleUpscaleHDF5 - ERR50'
        
        call GetData(NewChildFile%Info%GeometryFileName,                           &
                        Me%ObjEnterData, iflag,                            &
                        SearchType   = FromBlock,                          &
                        keyword      = 'CHILD_GEOMETRY',                  &
                        ClientModule = 'ConvertToHDF5', STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)'Error reading CHILD_GEOMETRY keyword'
            stop 'ConstructChildFile - ModuleUpscaleHDF5 - ERR60'
        end if
        
        if (Me%Upscale3D) then
            if (iflag == 0)then
                write(*,*)'Must specify name of CHILD_GEOMETRY file'
                stop 'ConstructChildFile - ModuleUpscaleHDF5 - ERR61'
            end if
        end if
        
    end subroutine ConstructChildFile

    !--------------------------------------------------------------------------

    subroutine AddChildFile

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        type (T_Child),     pointer                  :: PreviousChildFile, NewChildFile
        !Begin-----------------------------------------------------------------

        !Allocates new child file to be upscaled
        allocate (NewChildFile)
        nullify  (NewChildFile%Next)
        call ConstructChildFile(NewChildFile)

        !Put Child files in list by order of priority
        if (.not. associated(Me%FirstChildFile)) then
            call CopyChildFile(Me%FirstChildFile, NewChildFile)
            !deallocate(Me%FirstChildFile%Next) 
            !nullify(Me%FirstChildFile%Next)
        else
            !check next files in list
            !locate previous file in the first file
            allocate(PreviousChildFile)
            PreviousChildFile => Me%FirstChildFile                   
            !Pointed to the first file and will now iterate 
            !through all allocated files and add the new file in the end of the list
            do while(associated(PreviousChildFile))
                if (.not. associated(PreviousChildFile%Next)) then
                    !current file is the last file in the list of relevant files
                    !pointing the last item in the list to the newest user given file.
                    call CopyChildFile(PreviousChildFile%Next, NewChildFile)
                    !current file was added to list
                    exit
                end if

                !check next file in list
                PreviousChildFile => PreviousChildFile%Next
            end do
        end if

        nullify(NewChildFile)

    end subroutine AddChildFile

    !--------------------------------------------------------------------------

    subroutine CopyChildFile(ChildFileNew, ChildFileX)

        !Arguments-------------------------------------------------------------
        type(T_Child), pointer                     :: ChildFileNew
        type(T_Child), pointer                     :: ChildFileX

        !Local-----------------------------------------------------------------
      
        !Begin-----------------------------------------------------------------

        !Copies the values of the fields of a FatherFile to another FatherFile
        allocate(ChildFileNew)

        ChildFileNew%Info%InputFileName     =  ChildFileX%Info%InputFileName
        ChildFileNew%Info%GridFileName      =  ChildFileX%Info%GridFileName
        ChildFileNew%Info%Level             =  ChildFileX%Info%Level
        ChildFileNew%Info%NRemoveFrame      =  ChildFileX%Info%NRemoveFrame

    end subroutine CopyChildFile

    !--------------------------------------------------------------------------

    subroutine ConstructChildGrid(ObjChildFile)

        !Arguments-------------------------------------------------------------
        type (T_Child),      pointer               :: ObjChildFile
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer, dimension(:,:,:), pointer          :: WaterPoints3D
        real, dimension(:, :), pointer              :: SurfaceElevation
        !Local-----------------------------------------------------------------
        logical                                     :: exist
        !Begin-----------------------------------------------------------------
       
        write(*,*)'Constructing child grid...', trim(ObjChildFile%Info%GridFileName)

        !Verifies if file exists
        inquire(FILE = ObjChildFile%Info%GridFileName, EXIST = exist)
        if (.not. exist) then
            write(*,*)'Child Grid file does not exist'
            stop 'ConstructChildGrid - ModuleUpscaleHDF5 - ERR10'
        endif

        call ConstructHorizontalGrid(ObjChildFile%Info%ObjHorizontalGrid,              & 
                                     ObjChildFile%Info%GridFileName, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructChildGrid - ModuleUpscaleHDF5 - ERR20'

        call GetHorizontalGridSize(ObjChildFile%Info%ObjHorizontalGrid,                &
                                   WorkSize = ObjChildFile%Info%WorkSize2D,             &
                                   Size     = ObjChildFile%Info%Size2D,                 &
                                   STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructChildGrid -  ModuleUpscaleHDF5 - ERR30'

        call ConstructGridData(GridDataID       = ObjChildFile%Info%ObjBathymetry,      &
                               HorizontalGridID = ObjChildFile%Info%ObjHorizontalGrid,  &
                               FileName         = ObjChildFile%Info%GridFileName,       &
                               STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructChildGrid -  ModuleUpscaleHDF5 - ERR40'

        call ConstructHorizontalMap(HorizontalMapID  = ObjChildFile%Info%ObjHorizontalMap,     &
                                    GridDataID       = ObjChildFile%Info%ObjBathymetry,        &
                                    HorizontalGridID = ObjChildFile%Info%ObjHorizontalGrid,    &
                                    ActualTime       = Me%BeginTime,                           & 
                                    STAT             = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructChildGrid -  ModuleUpscaleHDF5 - ERR50'

        if (Me%Upscale3D) then
            call ConstructGeometry      (GeometryID       = ObjChildFile%Info%ObjGeometry,         &
                                            GridDataID       = ObjChildFile%Info%ObjBathymetry,       &
                                            HorizontalGridID = ObjChildFile%Info%ObjHorizontalGrid,   &
                                            HorizontalMapID  = ObjChildFile%Info%ObjHorizontalMap,    &
                                            ActualTime       = Me%BeginTime,               &
                                            NewDomain        = ObjChildFile%Info%GeometryFileName,    &
                                            STAT             = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructChildGrid -  ModuleUpscaleHDF5 - ERR60'

            call GetGeometrySize(GeometryID     = ObjChildFile%Info%ObjGeometry,                   &
                                    Size        = ObjChildFile%Info%Size3D,                        &
                                    WorkSize    = ObjChildFile%Info%WorkSize3D,                    &
                                    STAT        = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructChildGrid -  ModuleUpscaleHDF5 - ERR70'

            call ConstructMap ( Map_ID          = ObjChildFile%Info%ObjMap,                        &
                                GeometryID      = ObjChildFile%Info%ObjGeometry,                   &
                                HorizontalMapID = ObjChildFile%Info%ObjHorizontalMap,              &
                                TimeID          = Me%ObjTime,                           &
                                STAT            = STAT_CALL)  
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructChildGrid -  ModuleUpscaleHDF5 - ERR80'
        
            allocate(SurfaceElevation (Me%Father%Size3D%ILB:Me%Father%Size3D%IUB,             &
                                       Me%Father%Size3D%JLB:Me%Father%Size3D%JUB))
            SurfaceElevation(:,:) = 0
            
            call GetWaterPoints3D(ObjChildFile%Info%ObjMap, WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'ConstructChildGrid -  ModuleUpscaleHDF5 - ERR90'

            call ComputeInitialGeometry(GeometryID      = ObjChildFile%Info%ObjGeometry,           &
                                        WaterPoints3D   = WaterPoints3D,                &
                                        SurfaceElevation= SurfaceElevation,             &
                                        ActualTime      = Me%BeginTime,                 &
                                        STAT            = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructChildGrid -  ModuleUpscaleHDF5 - ERR100'
            
            call UngetMap(ObjChildFile%Info%ObjMap, WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructChildGrid - ModuleUpscaleHDF5 - ERR110'
            
            deallocate(SurfaceElevation)
        end if

    end subroutine ConstructChildGrid

    !------------------------------------------------------------------------

    subroutine ConstructGridPolygon(ObjGrid)

        !Arguments-------------------------------------------------------------
        type (T_Grid)                               :: ObjGrid
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer               :: GridLongitude, GridLatitude
        integer                                     :: CurrentVertix
        integer                                     :: i, j, ILB,IUB,JLB,JUB
        !Begin-----------------------------------------------------------------

        !Obtain coordinates of cells center:
        call GetGridLatitudeLongitude(ObjGrid%ObjHorizontalGrid,                    &
                                      GridLatitude = GridLatitude,                  & 
                                      GridLongitude = GridLongitude, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGridPolygon -  ModuleUpscaleHDF5 - ERR10'

        allocate(ObjGrid%Polygon)
        nullify(ObjGrid%Polygon%Next)

        if (ObjGrid%NRemoveFrame > 0) write(*,*) 'Removing', ObjGrid%NRemoveFrame, ' boundary points'

        ILB = ObjGrid%WorkSize2D%ILB + ObjGrid%NRemoveFrame
        JLB = ObjGrid%WorkSize2D%JLB + ObjGrid%NRemoveFrame
        IUB = ObjGrid%WorkSize2D%IUB - ObjGrid%NRemoveFrame
        JUB = ObjGrid%WorkSize2D%JUB - ObjGrid%NRemoveFrame

        !Construct the vertixes of polygon (should be in center of outmost cells):
        !get polygon number of vertixes
        ObjGrid%Polygon%Count = 2 * (IUB - ILB + 1) + 2 * (JUB - JLB + 1 - 2) + 1

        allocate(ObjGrid%Polygon%VerticesF(1:ObjGrid%Polygon%Count))

        CurrentVertix = 1
        
        do i = ILB, IUB
            ObjGrid%Polygon%VerticesF(CurrentVertix)%X  = GridLongitude(i,JLB)
            ObjGrid%Polygon%VerticesF(CurrentVertix)%Y  = GridLatitude(i,JLB)
            CurrentVertix = CurrentVertix + 1
        end do

        do j = JLB + 1, JUB-1
            ObjGrid%Polygon%VerticesF(CurrentVertix)%X  = GridLongitude(IUB,j)
            ObjGrid%Polygon%VerticesF(CurrentVertix)%Y  = GridLatitude(IUB,j)
            CurrentVertix = CurrentVertix + 1
        end do
        
        do i = ILB, IUB
            ObjGrid%Polygon%VerticesF(CurrentVertix)%X  = GridLongitude(IUB+ILB-i,JUB)
            ObjGrid%Polygon%VerticesF(CurrentVertix)%Y  = GridLatitude(IUB+ILB-i,JUB)
            CurrentVertix = CurrentVertix + 1
        end do

        do j = JLB + 1, JUB-1
            ObjGrid%Polygon%VerticesF(CurrentVertix)%X  = GridLongitude(ILB,JUB+JLB-j)
            ObjGrid%Polygon%VerticesF(CurrentVertix)%Y  = GridLatitude(ILB,JUB+JLB-j)
            CurrentVertix = CurrentVertix + 1
        end do

        !close polygon
        ObjGrid%Polygon%VerticesF(CurrentVertix)%X  = ObjGrid%Polygon%VerticesF(1)%X
        ObjGrid%Polygon%VerticesF(CurrentVertix)%Y  = ObjGrid%Polygon%VerticesF(1)%Y

        call SetLimitsPolygon(ObjGrid%Polygon)
       
    end subroutine ConstructGridPolygon

    !--------------------------------------------------------------------------
    
    subroutine SetLimitsPolygon(Polygon)

        !Arguments-------------------------------------------------------------
        type (T_Polygon)                  :: Polygon
        !Begin-----------------------------------------------------------------

        Polygon%Limits%Left   = minval(Polygon%VerticesF%X)
        Polygon%Limits%Right  = maxval(Polygon%VerticesF%X)
        Polygon%Limits%Bottom = minval(Polygon%VerticesF%Y)
        Polygon%Limits%Top    = maxval(Polygon%VerticesF%Y)

    end subroutine SetLimitsPolygon
    
    !--------------------------------------------------------------------------

    subroutine ConstructNewGrid
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        !Local-----------------------------------------------------------------
        logical                                     :: exist
        integer                                     :: HDF5_READ
        integer, dimension(:,:,:), pointer          :: WaterPoints3D
        real, dimension(:, :), pointer              :: SurfaceElevation
        !Begin-----------------------------------------------------------------

        write(*,*)'Constructing father grid...'
       
        !Verifies if file exists
        inquire(FILE = Me%Father%InputFileName, EXIST = exist)
        if (.not. exist) then
            write(*,*)'Grid file does not exist'
            stop 'ConstructNewGrid - ModuleUpscaleHDF5 - ERR10'
        endif
        
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
        
        !Open HDF5 file
        call ConstructHDF5 (Me%Father%ObjHDF5, trim(Me%Father%InputFileName), HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructNewGrid - ModuleUpscaleHDF5 - ERR20' 

        if (Me%Father%GridFromHDF5) then
            call ConstructHorizontalGrid(HorizontalGridID = Me%Father%ObjHorizontalGrid,&
                                         HDF5ID           = Me%Father%ObjHDF5,          &
                                         STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFatherGrid - ModuleInterpolateGrids - ERR30'
        else
            call ConstructHorizontalGrid(Me%Father%ObjHorizontalGrid, Me%Father%InputFileGrid, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid - ModuleUpscaleHDF5 - ERR30'
        end if

        call GetHorizontalGridSize  (Me%Father%ObjHorizontalGrid,                      &
                                     WorkSize = Me%Father%WorkSize2D,                  &
                                     Size     = Me%Father%Size2D, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleUpscaleHDF5 - ERR40'
        
        if (Me%Father%Batim_from_file) then
            
            allocate(Me%Father%Bathym(Me%Father%Size2D%ILB:Me%Father%Size2D%IUB,Me%Father%Size2D%JLB:Me%Father%Size2D%JUB))
            
            call HDF5SetLimits  (Me%Father%ObjHDF5, Me%Father%WorkSize2D%ILB,                     & 
                                 Me%Father%WorkSize2D%IUB,Me%Father%WorkSize2D%JLB,               & 
                                 Me%Father%WorkSize2D%JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructNewGrid - ModuleUpscaleHDF5 - ERR41'
        
            call HDF5ReadData(HDF5ID       = Me%Father%ObjHDF5,             &
                              GroupName    = "/Grid",                       &
                              Name         = "Bathymetry", Array2D = Me%Father%Bathym, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructNewGrid - ModuleUpscaleHDF5 - ERR42'
            
            call ConstructGridData(GridDataID       = Me%Father%ObjBathymetry,          &
                                   HorizontalGridID = Me%Father%ObjHorizontalGrid,      &
                                   InMatrix2D       = Me%Father%Bathym,                 &
                                   STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructNewGrid - ModuleUpscaleHDF5 - ERR43'
        else
            call ConstructGridData      (GridDataID       = Me%Father%ObjBathymetry,        &
                                         HorizontalGridID = Me%Father%ObjHorizontalGrid,    &
                                         FileName         = Me%Father%InputFileGrid, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleUpscaleHDF5 - ERR50'
        end if
        
        call set_module_time
        
        call ConstructHorizontalMap (HorizontalMapID  = Me%Father%ObjHorizontalMap,     &
                                     GridDataID       = Me%Father%ObjBathymetry,        &
                                     HorizontalGridID = Me%Father%ObjHorizontalGrid,    &
                                     ActualTime       = Me%BeginTime, STAT = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleUpscaleHDF5 - ERR60'

        if (Me%Upscale3D) then
            call ConstructGeometry (GeometryID       = Me%Father%ObjGeometry,         &
                                    GridDataID       = Me%Father%ObjBathymetry,       &
                                    HorizontalGridID = Me%Father%ObjHorizontalGrid,   &
                                    HorizontalMapID  = Me%Father%ObjHorizontalMap,    &
                                    ActualTime       = Me%BeginTime,               &
                                    NewDomain        = Me%Father%GeometryFileName, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleUpscaleHDF5 - ERR70'

            call GetGeometrySize (GeometryID = Me%Father%ObjGeometry,                   &
                                    Size     = Me%Father%Size3D,                        &
                                    WorkSize = Me%Father%WorkSize3D, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleUpscaleHDF5 - ERR80'

            call ConstructMap ( Map_ID          = Me%Father%ObjMap,                        &
                                GeometryID      = Me%Father%ObjGeometry,                   &
                                HorizontalMapID = Me%Father%ObjHorizontalMap,              &
                                TimeID          = Me%ObjTime, STAT = STAT_CALL)  
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleUpscaleHDF5 - ERR90'
        
            call GetWaterPoints3D(Me%Father%ObjMap, WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'ConstructNewGrid -  ModuleUpscaleHDF5 - ERR100'

            allocate(SurfaceElevation (Me%Father%Size3D%ILB:Me%Father%Size3D%IUB,             &
                                       Me%Father%Size3D%JLB:Me%Father%Size3D%JUB))
            SurfaceElevation(:,:) = 0
            
            call ComputeInitialGeometry (GeometryID          = Me%Father%ObjGeometry,           &
                                            WaterPoints3D    = WaterPoints3D,                &
                                            SurfaceElevation = SurfaceElevation,             &
                                            ActualTime       = Me%BeginTime,                 &
                                            STAT             = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleUpscaleHDF5 - ERR110'
            
            call UngetMap(Me%Father%ObjMap, WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleUpscaleHDF5 - ERR110'
        end if
        
    end subroutine ConstructNewGrid

    !------------------------------------------------------------------------

    !subroutine ConstructRelevantChild
    !
    !    !Arguments-------------------------------------------------------------
    !    
    !    !External--------------------------------------------------------------
    !    integer                                     :: STAT_CALL
    !    !Local-----------------------------------------------------------------
    !    real, dimension(:,:), pointer               :: GridLongitude, GridLatitude
    !    integer                                     :: i, j
    !    type (T_Child),     pointer                 :: ObjChildFile
    !    !Begin-----------------------------------------------------------------
    !
    !    !Obtain coordinates of cells center of new grid:
    !    call GetGridLatitudeLongitude(Me%New%ObjHorizontalGrid,                     &
    !                                  GridLatitude = GridLatitude,                  & 
    !                                  GridLongitude = GridLongitude, STAT = STAT_CALL)
    !    if(STAT_CALL .ne. SUCCESS_) stop 'ConstructRelevantChild -  ModuleUpscaleHDF5 - ERR10'
    !    
    !    if (Me%Upscale3D) then
    !        !Allocate father matrix based on the inputfilename hdf5:
    !        allocate(Me%RelevantFather(Me%New%WorkSize3D%ILB:Me%New%WorkSize3D%IUB,     &
    !                                   Me%New%WorkSize3D%JLB:Me%New%WorkSize3D%JUB,     &
    !                                   Me%New%WorkSize3D%KLB:Me%New%WorkSize3D%KUB))
    !    else
    !        allocate(Me%RelevantFather(Me%New%WorkSize3D%ILB:Me%New%WorkSize3D%IUB,     &
    !                                   Me%New%WorkSize3D%JLB:Me%New%WorkSize3D%JUB,     &
    !                                   1:1))
    !    end if
    !    
    !
    !    !Constructs the connection between child and father domains
    !    
    !    
    !    
    !    
    !
    !    call UnGetHorizontalGrid(Me%New%ObjHorizontalGrid, GridLatitude, STAT = STAT_CALL)
    !    if(STAT_CALL .ne. SUCCESS_) stop 'ConstructRelevantChild -  ModuleUpscaleHDF5 - ERR20'
    !
    !    call UnGetHorizontalGrid(Me%New%ObjHorizontalGrid, GridLongitude, STAT = STAT_CALL)
    !    if(STAT_CALL .ne. SUCCESS_) stop 'ConstructRelevantChild -  ModuleUpscaleHDF5 - ERR30'
    !    
    !end subroutine ConstructRelevantChild

    !--------------------------------------------------------------------------
    
    subroutine FatherSonCommunication
    
        !Arguments-------------------------------------------------------------    
        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        !Local-----------------------------------------------------------------
        type (T_Child),     pointer         :: ObjChildFile
        logical, dimension(:), allocatable  :: pointNotUsed
        logical                             :: Distorted
        integer                             :: CoordType, CoordTypeSon
        integer                             :: ProjType, I_min, I_max, J_min, J_max, count, i, j, p
        logical                             :: ReadCartCorners
        integer                             :: UTM, MIL_PORT, SIMPLE_GEOG, GEOG
        integer                             :: GRID_COORD, NLRD, ILB, IUB, JLB, JUB, NumberOfNodes
        real                                :: XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW, XC, YC, xcenter, ycenter
        !Begin-----------------------------------------------------------------
    
        write(*,*)'Constructing communication between grids...'
        
        call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,             &
                               SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD, NLRD = NLRD)
        
        !Gets Coordinates in use
        call GetGridCoordType(Me%Father%ObjHorizontalGrid, CoordType, ReadCartCorners, ProjType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleUpscaleHDF5 - ERR10'

        call GetCheckDistortion(Me%Father%ObjHorizontalGrid, Distorted, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR40'
        
        if (Distorted) then
            Write (*,*) 'Tool is not yet ready for father distorted grids.'
            Write (*,*) 'Check if the distortion is real, or if it is but a rounding issue'
            !stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR30' 
        end if
        
        !Cycle all the childs coordinates - they must be the same as the parent domain
        ! check for distortion as well
        ObjChildFile => Me%FirstChildFile
    
        do while (associated(ObjChildFile))
            call GetGridCoordType(ObjChildFile%Info%ObjHorizontalGrid, CoordTypeSon, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR20'

            if (CoordType /= CoordTypeSon) then
                Write (*,*) 'Fathergrid coordinate type is different than son grid coordinate type'
                stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR30'
            endif
            ObjChildFile => ObjChildFile%next
        end do

        if(CoordType == UTM .or. CoordType == MIL_PORT .or.                             &
           CoordType == GRID_COORD .or. CoordType == NLRD)then
            call GetHorizontalGrid(Me%Father%ObjHorizontalGrid,                         & 
                                   XX_IE = Me%Father%ConnectionX,                       &
                                   YY_IE = Me%Father%ConnectionY, STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR140' 
        else
            call GetGridLatitudeLongitude(Me%Father%ObjHorizontalGrid,                  &
                                          GridLatitudeConn  = Me%Father%ConnectionY,    &
                                          GridLongitudeConn = Me%Father%ConnectionX, STAT  = STAT_CALL)
            if(STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR160'
        end if
        
        ObjChildFile => Me%FirstChildFile
        
        do while (associated(ObjChildFile))
            call GetGridCoordType(ObjChildFile%Info%ObjHorizontalGrid, CoordTypeSon, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR20'

            if(CoordTypeSon == UTM .or. CoordTypeSon == MIL_PORT .or.                             &
               CoordTypeSon == GRID_COORD .or. CoordTypeSon == NLRD)then
                call GetHorizontalGrid(ObjChildFile%Info%ObjHorizontalGrid,                         & 
                                       XX_IE = ObjChildFile%Info%ConnectionX,                       &
                                       YY_IE = ObjChildFile%Info%ConnectionY, STAT = STAT_CALL)
                if(STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR140'
            else
                call GetGridLatitudeLongitude(ObjChildFile%Info%ObjHorizontalGrid,                  &
                                              GridLatitudeConn  = ObjChildFile%Info%ConnectionY,    &
                                              GridLongitudeConn = ObjChildFile%Info%ConnectionX, STAT  = STAT_CALL)
                if(STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR160'
            end if
            ObjChildFile => ObjChildFile%next
        end do
        
        !-------------------------------------- construct connections vectors ----------------------------------------
        !build child points to be used in the average procedure
        
        ObjChildFile => Me%FirstChildFile
        do while (associated(ObjChildFile))   
            count = 0

            I_min = ObjChildFile%Info%WorkSize2D%ILB + ObjChildFile%Info%NRemoveFrame
            J_min = ObjChildFile%Info%WorkSize2D%JLB + ObjChildFile%Info%NRemoveFrame
            I_max = ObjChildFile%Info%WorkSize2D%IUB - ObjChildFile%Info%NRemoveFrame
            J_max = ObjChildFile%Info%WorkSize2D%JUB - ObjChildFile%Info%NRemoveFrame
            
            NumberOfNodes =  (I_max - I_min + 1) * (J_max - J_min + 1)
            
            allocate(ObjChildFile%Info%NodeX(NumberOfNodes))
            allocate(ObjChildFile%Info%NodeY(NumberOfNodes))
            allocate(ObjChildFile%Info%NodeZ(NumberOfNodes))

            !allocate(Me%State_open(NumberOfNodes))
            ObjChildFile%Info%NodeX(:) = 0
            ObjChildFile%Info%NodeY(:) = 0
            ObjChildFile%Info%NodeZ(:) = 0
            !ObjChildFile%Info%State_open(:) = 0
            
            do j = J_min, J_max
            do i = I_min, I_max
                count = count + 1
                xcenter = ((ObjChildFile%Info%ConnectionX(i, j  ) + ObjChildFile%Info%ConnectionX(i+1, j  ))/2. + &
                           (ObjChildFile%Info%ConnectionX(i, j+1) + ObjChildFile%Info%ConnectionX(i+1, j+1))/2.)/2.
                ObjChildFile%Info%NodeX(Count) = xcenter
    
                ycenter = ((ObjChildFile%Info%ConnectionY(i, j  ) + ObjChildFile%Info%ConnectionY(i+1, j  ))/2. + &
                           (ObjChildFile%Info%ConnectionY(i, j+1) + ObjChildFile%Info%ConnectionY(i+1, j+1))/2.)/2.             

                ObjChildFile%Info%NodeY(Count) = ycenter
                
                !ObjChildFile%Info%State_open(Count) = ObjChildFile%Info%WaterPoints2D(i, j)
            enddo
            enddo
            
            !build father link arrays--------------------------------------------------------------------------
            ILB = Me%Father%Size2D%ILB
            IUB = Me%Father%Size2D%IUB
            JLB = Me%Father%Size2D%JLB
            JUB = Me%Father%Size2D%JUB
            
            allocate(pointNotUsed(1:NumberOfNodes))
            allocate(ObjChildFile%Info%nPointsInside(ILB:IUB, JLB:JUB))
            allocate(ObjChildFile%Info%FatherPoint_ID_I(1:NumberOfNodes))
            allocate(ObjChildFile%Info%FatherPoint_ID_J(1:NumberOfNodes))
            ObjChildFile%Info%nPointsInside(:,:) = 0
            ObjChildFile%Info%FatherPoint_ID_I(:) = 0
            ObjChildFile%Info%FatherPoint_ID_J(:) = 0
            pointNotUsed(:) = .true.
            
            do j = Me%Father%WorkSize2D%JLB, Me%Father%WorkSize2D%JUB
            do i = Me%Father%WorkSize2D%ILB, Me%Father%WorkSize2D%IUB
                
                XSW = Me%Father%ConnectionX(i, j)
                YSW = Me%Father%ConnectionY(i, j)
                XSE = Me%Father%ConnectionX(i, j + 1)
                YSE = Me%Father%ConnectionY(i, j + 1)
                XNE = Me%Father%ConnectionX(i + 1, j + 1)
                YNE = Me%Father%ConnectionY(i + 1, j + 1)
                XNW = Me%Father%ConnectionX(i + 1, j)
                YNW = Me%Father%ConnectionY(i + 1, j)
                
                XC  = (XSW + XSE + XNE + XNW) / 4.
                YC  = (YSW + YSE + YNE + YNW) / 4.

                do p = 1, NumberOfNodes
                    if (pointNotUsed(p)) then
                        if ((ObjChildFile%Info%NodeX(p) < XSE) .and. (ObjChildFile%Info%NodeX(p) > XSW)) then
                            if ((ObjChildFile%Info%NodeY(p) < YNW) .and. (ObjChildFile%Info%NodeY(p) > YSW)) then
                                !point is inside the grid cell
                                ObjChildFile%Info%FatherPoint_ID_I(p) = i
                                ObjChildFile%Info%FatherPoint_ID_J(p) = j
                                pointNotUsed(p) = .false.
                            end if
                        end if
                    end if
                end do
            end do
            end do 
            
            deallocate(pointNotUsed)
            ObjChildFile => ObjChildFile%next
        end do     

    end subroutine FatherSonCommunication

    !------------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: HDF5_CREATE, ILB, IUB, JLB, JUB
        real,       dimension(:,:  ), pointer       :: Bathymetry
        integer,    dimension(:,:,:), pointer       :: WaterPoints3D
        integer,    dimension(:,:  ), pointer       :: WaterPoints2D
        !----------------------------------------------------------------------

        ILB = Me%Father%Size2D%ILB
        IUB = Me%Father%Size2D%IUB
        JLB = Me%Father%Size2D%JLB
        JUB = Me%Father%Size2D%JUB
        
        allocate(WaterPoints3D(ILB:IUB,JLB:JUB, 1))
        
        call GetWaterPoints2D(Me%Father%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleUpscaleHDF5 - ERR10'
        
        WaterPoints3D(ILB:IUB,JLB:JUB,1) = WaterPoints2D(ILB:IUB,JLB:JUB)

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        !Opens HDF5 File
        call ConstructHDF5(Me%New%ObjHDF5, Me%New%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleUpscaleHDF5 - ERR30'
        
        
        call HDF5SetLimits  (Me%New%ObjHDF5, Me%Father%WorkSize2D%ILB,                     & 
                             Me%Father%WorkSize2D%IUB,Me%Father%WorkSize2D%JLB,               & 
                             Me%Father%WorkSize2D%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleUpscaleHDF5 - ERR40'

        if (Me%Father%Batim_from_file) then
            Bathymetry => Me%Father%Bathym
        else
            call GetGridData(Me%Father%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleUpscaleHDF5 - ERR50'
        end if
        
        call HDF5WriteData   (Me%New%ObjHDF5, "/Grid", "Bathymetry", "-",               &
                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleUpscaleHDF5 - ERR60'            

        if (.not. Me%Father%Batim_from_file) then
            call UngetGridData(Me%Father%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleUpscaleHDF5 - ERR70'
        end if
        
        call WriteHorizontalGrid (Me%Father%ObjHorizontalGrid, Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleUpscaleHDF5 - ERR80'

        call HDF5SetLimits  (Me%New%ObjHDF5, Me%Father%WorkSize2D%ILB,                     & 
                             Me%Father%WorkSize2D%IUB, Me%Father%WorkSize2D%JLB,              &
                             Me%Father%WorkSize2D%JUB, 1, 1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleUpscaleHDF5 - ERR90'

        call HDF5WriteData   (Me%New%ObjHDF5, "/Grid", "WaterPoints3D", "-",              &
                              Array3D = WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleUpscaleHDF5 - ERR100'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%New%ObjHDF5, STAT = STAT_CALL)                         
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleUpscaleHDF5 - ERR110'

        call UngetHorizontalMap(Me%Father%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleUpscaleHDF5 - ERR120'

        deallocate(WaterPoints3D)
        nullify   (WaterPoints3D)

    end subroutine Open_HDF5_OutPut_File
    
    !--------------------------------------------------------------------------
    
    subroutine Open_HDF5_OutPut_File3D

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE
        real,       dimension(:,:  ), pointer       :: Bathymetry
        !----------------------------------------------------------------------

        call GetWaterPoints3D(Me%Father%ObjMap, Me%Father%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleUpscaleHDF5 - ERR10'

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        !Opens HDF5 File
        call ConstructHDF5(Me%New%ObjHDF5, Me%New%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleUpscaleHDF5 - ERR20'
        
        call HDF5SetLimits  (Me%New%ObjHDF5, Me%Father%WorkSize2D%ILB, Me%Father%WorkSize2D%IUB,&
                             Me%Father%WorkSize2D%JLB, Me%Father%WorkSize2D%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleUpscaleHDF5 - ERR30'
        
        if (Me%Father%Batim_from_file) then
            Bathymetry => Me%Father%Bathym
        else
            call GetGridData(Me%Father%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File3D - ModuleUpscaleHDF5 - ERR40'
        end if
        
        call HDF5WriteData   (Me%New%ObjHDF5, "/Grid", "Bathymetry", "-",               &
                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File3D - ModuleUpscaleHDF5 - ERR50'            

        if (.not. Me%Father%Batim_from_file) then
            call UngetGridData(Me%Father%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File3D - ModuleUpscaleHDF5 - ERR60'
        end if

        call WriteHorizontalGrid (Me%Father%ObjHorizontalGrid, Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleUpscaleHDF5 - ERR70'

        call HDF5SetLimits  (Me%New%ObjHDF5, Me%Father%WorkSize3D%ILB, Me%Father%WorkSize3D%IUB,&
                             Me%Father%WorkSize3D%JLB, Me%Father%WorkSize3D%JUB, Me%Father%WorkSize3D%KLB, &
                             Me%Father%WorkSize3D%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleUpscaleHDF5 - ERR80'

        call HDF5WriteData   (Me%New%ObjHDF5, "/Grid", "WaterPoints3D", "-",    &
                              Array3D = Me%Father%WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleUpscaleHDF5 - ERR90'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleUpscaleHDF5 - ERR100'

        call UngetMap(Me%Father%ObjMap, Me%Father%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleUpscaleHDF5 - ERR110'

    end subroutine Open_HDF5_OutPut_File3D

    !--------------------------------------------------------------------------

    subroutine OpenAndReadChildHDF5Files

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: CurrentInstant, CurrentInstant_child
        integer                                     :: StartInstant, EndInstant, StartInstant_child, EndInstant_child
        type(T_Time)                                :: CurrentDate, CurrentDate_child
        integer                                     :: CurrentProperty 
        type(T_Field), pointer                      :: NewField, AuxField
        character(len=StringLength)                 :: PropertyName
        integer                                     :: Rank
        integer, dimension(7)                       :: Dimensions
        integer                                     :: Count, CurrentProperty_child, HDF5_READ
        integer                                     :: i, j, ILB, IUB, JLB, JUB, KLB, KUB, p, k, n, KUB_child
        integer                                     :: ILB_child, IUB_child, JLB_child, JUB_child, KLB_child
        type (T_Child),     pointer                 :: ObjChildFile
        logical                                     :: ConvertThisField
        !Begin-----------------------------------------------------------------

        Count = 1

        nullify(Me%FirstChildFile%NewField, NewField)
        
        !1 - ABRIR e construir INPUT FILE - DONE
        !Get number of instants
        !call get_timeinstants(Me%Father, StartInstant, EndInstant)
        
        StartInstant = 1
        EndInstant = Me%Father%NumberOfInstants
        
        !2 - Fazer o ciclo de propriedades
        !Gets father HDF number of properties
        call get_nproperties(Me%Father)
        
        !write time in output hdf5
        do CurrentInstant = StartInstant, EndInstant
            call OutputInstants(CurrentInstant)
        end do
        
        !Allocates auxiliar fields---------------------------------------
        allocate (NewField)
        allocate (AuxField)
        
        !Allocates new instance for every child------------------------
        ObjChildFile => Me%FirstChildFile 

        do while (associated(ObjChildFile))
            
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
        
            !Open HDF5 file
            call ConstructHDF5 (ObjChildFile%Info%ObjHDF5, trim(ObjChildFile%Info%InputFileName), HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR01'
            
            allocate (ObjChildFile%NewField)
            ObjChildFile => ObjChildFile%next
        end do
        
        !Cycle father properties
        do CurrentProperty = 1, Me%Father%NumberOfProperties
        
            call GetHDF5GroupID(Me%Father%ObjHDF5, "/Results",        &
                                CurrentProperty, PropertyName, STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR10'

            if(.not. Me%ConvertAllFields)then
                ConvertThisField = .false.
                do n = 1, Me%nFieldsToUpscale
                    if(Me%FieldsToUpscale(n) == PropertyName) ConvertThisField = .true.
                end do
                if(.not. ConvertThisField) cycle
            end if

            write(*,*)'Reading '//trim(PropertyName)//' fields'
            
            !Father time cycle
            do CurrentInstant = StartInstant, EndInstant
                
                CurrentDate = Me%Father%InstantsArray(CurrentInstant)
                
                !Get father ID field
                call GetHDF5GroupID(Me%Father%ObjHDF5, "/Results/"//trim(PropertyName), &
                                CurrentInstant, NewField%Name, NewField%Units, Rank, Dimensions, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR20'
                
                NewField%Name = trim(PropertyName)
                
                select case (Rank)
                    
                !Father field is 2D
                case(2)
                    !check dimensions
                    if(Dimensions(1) .ne. Me%Father%WorkSize2D%IUB) then
                        write(*,*)'Fields size is not consistent with grid size : '//trim(Me%Father%InputFileName)
                        stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR30'
                    end if

                    if(Dimensions(2) .ne. Me%Father%WorkSize2D%JUB) then
                        write(*,*)'Fields size is not consistent with grid size : '//trim(Me%Father%InputFileName)
                        stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR40'
                    end if 
                    
                    ILB = Me%Father%Size2D%ILB
                    IUB = Me%Father%Size2D%IUB
                    JLB = Me%Father%Size2D%JLB
                    JUB = Me%Father%Size2D%JUB
                
                    !allocate new field
                    nullify (NewField%Values2D)
                    allocate(NewField%Values2D(ILB:IUB, JLB:JUB))
                    allocate(AuxField%Values2D(ILB:IUB, JLB:JUB))
                    
                    call HDF5SetLimits (Me%Father%ObjHDF5,                 &
                                    Me%Father%WorkSize2D%ILB, Me%Father%WorkSize2D%IUB, &
                                    Me%Father%WorkSize2D%JLB, Me%Father%WorkSize2D%JUB, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR50'

                    !read father field and save it in newfield 2D matrix
                    call HDF5ReadData(Me%Father%ObjHDF5, "/Results/"//trim(PropertyName), trim(PropertyName), &
                                      Array2D      = NewField%Values2D,   &
                                      OutputNumber = CurrentInstant, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR60'
                    
                    !Start computations
                    
                    !Cycle all child domains----------------------------------------------------------
                    ObjChildFile => Me%FirstChildFile 

                    do while (associated(ObjChildFile))
                        
                        call get_nproperties(ObjChildFile%Info)
                                
                        do CurrentProperty_child = 1, ObjChildFile%Info%NumberOfProperties
                                
                            call GetHDF5GroupID(ObjChildFile%Info%ObjHDF5, "/Results", CurrentProperty_child, PropertyName, &
                                                STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR70'

                            if(.not. Me%ConvertAllFields)then
                                ConvertThisField = .false.
                                do n = 1, Me%nFieldsToUpscale
                                    if(Me%FieldsToUpscale(n) == PropertyName) ConvertThisField = .true.
                                end do
                                if(.not. ConvertThisField) cycle
                            end if

                            write(*,*)'Reading '//trim(PropertyName)//' child fields'
                            
                            !Get time instants
                            call get_timeinstants(ObjChildFile%Info, StartInstant_child, EndInstant_child)
                        
                            do CurrentInstant_child = StartInstant_child, EndInstant_child
                                !find current full date+time of child hdf5
                                CurrentDate_child = ObjChildFile%Info%InstantsArray(CurrentInstant_child)
                            
                                !Check if time in child HDF is equal to that of the father hdf5
                                if (CurrentDate == CurrentDate_child) then
                                    !start computation for current property
                                    !Get child ID field
                                    call GetHDF5GroupID(ObjChildFile%Info%ObjHDF5, "/Results/"//trim(PropertyName), &
                                                    CurrentInstant, ObjChildFile%Info%Name, ObjChildFile%Info%Units, Rank, &
                                                    Dimensions, STAT = STAT_CALL)                                
                                    if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR80'
                                    
                                    if(Dimensions(1) .ne. ObjChildFile%Info%WorkSize2D%IUB) then
                                        write(*,*)'Fields size is not consistent with grid size :'
                                        write(*,*)' '//trim(ObjChildFile%Info%InputFileName)
                                        stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR90'
                                    end if

                                    if(Dimensions(2) .ne. ObjChildFile%Info%WorkSize2D%JUB) then
                                        write(*,*)'Fields size is not consistent with grid size :'
                                        write(*,*)' '//trim(ObjChildFile%Info%InputFileName)
                                        stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR100'
                                    end if 
                    
                                    ILB_child = ObjChildFile%Info%Size2D%ILB
                                    IUB_child = ObjChildFile%Info%Size2D%IUB
                                    JLB_child = ObjChildFile%Info%Size2D%JLB
                                    JUB_child = ObjChildFile%Info%Size2D%JUB
                
                                    !allocate new field
                                    nullify (ObjChildFile%NewField%Values2D)
                                    allocate(ObjChildFile%NewField%Values2D(ILB:IUB, JLB:JUB))
                    
                                    call HDF5SetLimits (ObjChildFile%Info%ObjHDF5,                 &
                                                    ObjChildFile%Info%WorkSize2D%ILB, ObjChildFile%Info%WorkSize2D%IUB, &
                                                    ObjChildFile%Info%WorkSize2D%JLB, ObjChildFile%Info%WorkSize2D%JUB, &
                                                    STAT = STAT_CALL)
                                    if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR120'

                                    !read father field and save it in newfield 2D matrix
                                    call HDF5ReadData(ObjChildFile%Info%ObjHDF5, "/Results/"//trim(PropertyName), &
                                                    trim(PropertyName), Array2D = ObjChildFile%NewField%Values2D, &
                                                    OutputNumber = CurrentInstant_child, STAT = STAT_CALL)
                                    if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR20130'
                                    
                                    !AuxField is used only to account for the overlapped area between father and child domains
                                    AuxField%Values2D(:,:) = 0
                                    ObjChildFile%Info%nPointsInside(:,:) = 0
                                    
                                    !Fill NodeZ array of child HDF5
                                    ILB_child = ObjChildFile%Info%WorkSize2D%ILB + ObjChildFile%Info%NRemoveFrame
                                    JLB_child = ObjChildFile%Info%WorkSize2D%JLB + ObjChildFile%Info%NRemoveFrame
                                    IUB_child = ObjChildFile%Info%WorkSize2D%IUB - ObjChildFile%Info%NRemoveFrame
                                    JUB_child = ObjChildFile%Info%WorkSize2D%JUB - ObjChildFile%Info%NRemoveFrame
                                    Count = 0
                                    do j = JLB_child, JUB_child
                                    do i = ILB_child, IUB_child
                                        Count = Count + 1
                                        ObjChildFile%Info%NodeZ(Count) = ObjChildFile%NewField%Values2D(i, j)
                                    enddo
                                    enddo
                                    
                                    !Fill every Father grid cell with available child points
                                    do p = 1, size(ObjChildFile%Info%NodeZ)
                                        i = ObjChildFile%Info%FatherPoint_ID_I(p)
                                        j = ObjChildFile%Info%FatherPoint_ID_J(p)
                                        if (i + j > 1) then
                                            if (ObjChildFile%Info%NodeZ(p) > HalfFillValueReal) then
                                                if (NewField%Values2D(i,j) > HalfFillValueReal) then
                                                    ObjChildFile%Info%nPointsInside(i, j) = ObjChildFile%Info%nPointsInside(i, j) + 1
                                                    AuxField%Values2D(i, j) = AuxField%Values2D(i, j) + ObjChildFile%Info%NodeZ(p)
                                                end if
                                            end if
                                        end if
                                    end do
                
                                    do j = Me%Father%WorkSize2D%JLB, Me%Father%WorkSize2D%JUB
                                    do i = Me%Father%WorkSize2D%ILB, Me%Father%WorkSize2D%IUB
                                        if (ObjChildFile%Info%nPointsInside(i, j) > 0) then
                                            NewField%Values2D(i, j) =  AuxField%Values2D(i, j) &
                                                                    / ObjChildFile%Info%nPointsInside(i, j)
                                        end if
                                    end do
                                    end do
                                end if
                            end do
                        end do
                        
                        !Write into output HDF5 output file
                        call OutputFields(NewField, CurrentInstant)
                        
                        ObjChildFile => ObjChildFile%next
                    end do
                    
                !Father field is 3D
                case(3)
                    !Correct worksize dimensions when fields are 2D but matrixes are 3D (only one vertical layer)
                    if (Me%Father%WorkSize3D%ILB < 1) Me%Father%WorkSize3D%ILB = Me%Father%WorkSize2D%ILB
                    if (Me%Father%WorkSize3D%JLB < 1) Me%Father%WorkSize3D%JLB = Me%Father%WorkSize2D%JLB
                    if (Me%Father%WorkSize3D%KLB < 1) Me%Father%WorkSize3D%KLB = 1
                    if (Me%Father%WorkSize3D%IUB < 1) Me%Father%WorkSize3D%IUB = Me%Father%WorkSize2D%IUB
                    if (Me%Father%WorkSize3D%JUB < 1) Me%Father%WorkSize3D%JUB = Me%Father%WorkSize2D%JUB
                    if (Me%Father%WorkSize3D%KUB < 1) Me%Father%WorkSize3D%KUB = 1
                    
                    !check dimensions
                    if(Dimensions(1) .ne. Me%Father%WorkSize3D%IUB) then
                        write(*,*)'Fields size is not consistent with grid size : '//trim(Me%Father%InputFileName)
                        stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR140'
                    end if

                    if(Dimensions(2) .ne. Me%Father%WorkSize3D%JUB) then
                        write(*,*)'Fields size is not consistent with grid size : '//trim(Me%Father%InputFileName)
                        stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR150'
                    end if
                    
                    if(Dimensions(3) .ne. Me%Father%WorkSize3D%KUB) then
                        write(*,*)'Fields size is not consistent with grid size : '//trim(Me%Father%InputFileName)
                        stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR160'
                    end if 
                    
                    ILB = Me%Father%Size3D%ILB
                    IUB = Me%Father%Size3D%IUB
                    JLB = Me%Father%Size3D%JLB
                    JUB = Me%Father%Size3D%JUB
                    KLB = Me%Father%Size3D%KLB
                    KUB = Me%Father%Size3D%KUB
                    
                    if (ILB < 0) ILB = Me%Father%Size2D%ILB
                    if (JLB < 0) JLB = Me%Father%Size2D%JLB
                    if (KLB < 0) KLB = 0
                    if (IUB < 0) IUB = Me%Father%Size2D%IUB
                    if (JUB < 0) JUB = Me%Father%Size2D%JUB
                    if (KUB < 0) KUB = 2
                    
                    !allocate new field
                    nullify (NewField%Values3D)
                    allocate(NewField%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate(AuxField%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))
            
                    call HDF5SetLimits (Me%Father%ObjHDF5,                 &
                                    Me%Father%WorkSize3D%ILB, Me%Father%WorkSize3D%IUB, &
                                    Me%Father%WorkSize3D%JLB, Me%Father%WorkSize3D%JUB, &
                                    Me%Father%WorkSize3D%KLB, Me%Father%WorkSize3D%KUB, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR170'

                    !read father field and save it in newfield 3D matrix with just 1 layer
                    call HDF5ReadData(Me%Father%ObjHDF5, "/Results/"//trim(PropertyName), trim(PropertyName), &
                                      Array3D      = NewField%Values3D,   &
                                      OutputNumber = CurrentInstant, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR180'
                    
                    !Start computations
                    !Cycle all child domains----------------------------------------------------------
                    ObjChildFile => Me%FirstChildFile 

                    do while (associated(ObjChildFile))
                        
                        call get_nproperties(ObjChildFile%Info)
                                
                        do CurrentProperty_child = 1, ObjChildFile%Info%NumberOfProperties
                                
                            call GetHDF5GroupID(ObjChildFile%Info%ObjHDF5, "/Results", CurrentProperty_child, PropertyName, &
                                                STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR190'

                            if(.not. Me%ConvertAllFields)then
                                ConvertThisField = .false.
                                do n = 1, Me%nFieldsToUpscale
                                    if(Me%FieldsToUpscale(n) == PropertyName) ConvertThisField = .true.
                                end do
                                if(.not. ConvertThisField) cycle
                            end if

                            write(*,*)'Reading '//trim(PropertyName)//' child fields'
                            
                            !Get time instants
                            call get_timeinstants(ObjChildFile%Info, StartInstant_child, EndInstant_child)
                        
                            do CurrentInstant_child = StartInstant_child, EndInstant_child
                                !find current full date+time of child hdf5
                                CurrentDate_child = ObjChildFile%Info%InstantsArray(CurrentInstant_child)
                            
                                !Check if time in child HDF is equal to that of the father hdf5
                                if (CurrentDate == CurrentDate_child) then
                                    !start computation for current property
                                    !Get child ID field
                                    call GetHDF5GroupID(ObjChildFile%Info%ObjHDF5, "/Results/"//trim(PropertyName), &
                                                    CurrentInstant, ObjChildFile%Info%Name, ObjChildFile%Info%Units, &
                                                    Rank, Dimensions, STAT = STAT_CALL)                                
                                    if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR200'
                                    
                                    if (ObjChildFile%Info%WorkSize3D%ILB < 1) &
                                        ObjChildFile%Info%WorkSize3D%ILB = ObjChildFile%Info%WorkSize2D%ILB
                                    if (ObjChildFile%Info%WorkSize3D%JLB < 1) &
                                        ObjChildFile%Info%WorkSize3D%JLB = ObjChildFile%Info%WorkSize2D%JLB
                                    if (ObjChildFile%Info%WorkSize3D%KLB < 1) &
                                        ObjChildFile%Info%WorkSize3D%KLB = 1
                                    if (ObjChildFile%Info%WorkSize3D%IUB < 1) &
                                        ObjChildFile%Info%WorkSize3D%IUB = ObjChildFile%Info%WorkSize2D%IUB
                                    if (ObjChildFile%Info%WorkSize3D%JUB < 1) &
                                        ObjChildFile%Info%WorkSize3D%JUB = ObjChildFile%Info%WorkSize2D%JUB
                                    if (ObjChildFile%Info%WorkSize3D%KUB < 1) &
                                        ObjChildFile%Info%WorkSize3D%KUB = 1
                                    
                                    if(Dimensions(1) .ne. ObjChildFile%Info%WorkSize3D%IUB) then
                                        write(*,*)'Fields size not consistent with grid size :'
                                        write(*,*)' '//trim(ObjChildFile%Info%InputFileName)
                                        stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR210'
                                    end if

                                    if(Dimensions(2) .ne. ObjChildFile%Info%WorkSize3D%JUB) then
                                        write(*,*)'Fields size not consistent with grid size :'
                                        write(*,*)' '//trim(ObjChildFile%Info%InputFileName)
                                        stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR220'
                                    end if 
                                    
                                    if(Dimensions(3) .ne. ObjChildFile%Info%WorkSize3D%KUB) then
                                        write(*,*)'Fields size not consistent with grid size :'
                                        write(*,*)' '//trim(ObjChildFile%Info%InputFileName)
                                        stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR230'
                                    end if
                                    
                                    ILB_child = ObjChildFile%Info%Size3D%ILB
                                    IUB_child = ObjChildFile%Info%Size3D%IUB
                                    JLB_child = ObjChildFile%Info%Size3D%JLB
                                    JUB_child = ObjChildFile%Info%Size3D%JUB
                                    KLB_child = ObjChildFile%Info%Size3D%KLB
                                    KUB_child = ObjChildFile%Info%Size3D%KUB
                                    
                                    if (ILB_child < 0) ILB_child = ObjChildFile%Info%Size2D%ILB
                                    if (JLB_child < 0) JLB_child = ObjChildFile%Info%Size2D%JLB
                                    if (KLB_child < 0) KLB_child = 0
                                    if (IUB_child < 0) IUB_child = ObjChildFile%Info%Size2D%IUB
                                    if (JUB_child < 0) JUB_child = ObjChildFile%Info%Size2D%JUB
                                    if (KUB_child < 0) KUB_child = 2
                                    
                                    !allocate new field
                                    nullify (ObjChildFile%NewField%Values3D)
                                    allocate(ObjChildFile%NewField%Values3D(ILB_child:IUB_child, &
                                                                            JLB_child:JUB_child, &
                                                                            KLB_child:KUB_child))
                    
                                    call HDF5SetLimits (ObjChildFile%Info%ObjHDF5,                 &
                                                    ObjChildFile%Info%WorkSize3D%ILB, ObjChildFile%Info%WorkSize3D%IUB, &
                                                    ObjChildFile%Info%WorkSize3D%JLB, ObjChildFile%Info%WorkSize3D%JUB, &
                                                    ObjChildFile%Info%WorkSize3D%KLB, ObjChildFile%Info%WorkSize3D%KUB, &
                                                    STAT = STAT_CALL)
                                    if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR240'

                                    !read father field and save it in newfield 3D matrix with just 1 layer
                                    call HDF5ReadData(ObjChildFile%Info%ObjHDF5, "/Results/"//trim(PropertyName), &
                                                    trim(PropertyName), Array3D = ObjChildFile%NewField%Values3D, &
                                                    OutputNumber = CurrentInstant_child, STAT = STAT_CALL)
                                    if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR250'
                                    
                                    !AuxField is used only to account for the overlapped area between father and child domains
                                    AuxField%Values3D(:,:,:) = 0
                                    
                                    do k = ObjChildFile%Info%WorkSize3D%KLB, ObjChildFile%Info%WorkSize3D%KUB
                                        
                                        ObjChildFile%Info%nPointsInside(:,:) = 0
                                        
                                        ILB_child = ObjChildFile%Info%WorkSize3D%ILB + ObjChildFile%Info%NRemoveFrame
                                        JLB_child = ObjChildFile%Info%WorkSize3D%JLB + ObjChildFile%Info%NRemoveFrame
                                        IUB_child = ObjChildFile%Info%WorkSize3D%IUB - ObjChildFile%Info%NRemoveFrame
                                        JUB_child = ObjChildFile%Info%WorkSize3D%JUB - ObjChildFile%Info%NRemoveFrame
                                        !Fill NodeZ array of child HDF5
                                        Count = 0
                                        do j = JLB_child, JUB_child
                                        do i = ILB_child, IUB_child
                                            Count = Count + 1
                                            ObjChildFile%Info%NodeZ(Count) = ObjChildFile%NewField%Values3D(i, j, k)
                                        enddo
                                        enddo
                                    
                                        !Fill every Father grid cell with available child points
                                        do p = 1, size(ObjChildFile%Info%NodeZ)
                                            i = ObjChildFile%Info%FatherPoint_ID_I(p)
                                            j = ObjChildFile%Info%FatherPoint_ID_J(p)
                                            if (i + j > 1) then
                                                if (ObjChildFile%Info%NodeZ(p) > HalfFillValueReal) then
                                                    if (NewField%Values3D(i,j,k) > HalfFillValueReal) then
                                                        ObjChildFile%Info%nPointsInside(i,j) = &
                                                            ObjChildFile%Info%nPointsInside(i,j) + 1
                                                        AuxField%Values3D(i,j,k) = AuxField%Values3D(i,j,k) &
                                                                                 + ObjChildFile%Info%NodeZ(p)
                                                    end if
                                                end if
                                            end if
                                        end do
                
                                        do j = Me%Father%WorkSize3D%JLB, Me%Father%WorkSize3D%JUB
                                        do i = Me%Father%WorkSize3D%ILB, Me%Father%WorkSize3D%IUB
                                            if (ObjChildFile%Info%nPointsInside(i, j) > 0) then
                                                NewField%Values3D(i, j, k) = AuxField%Values3D(i, j, k) &
                                                                           / ObjChildFile%Info%nPointsInside(i, j)
                                            end if
                                        end do
                                        end do
                                    end do
                                end if
                            end do
                        end do
                        
                        !Write into output HDF5 output file
                        call OutputFields3D(NewField, CurrentInstant)
                        ObjChildFile => ObjChildFile%next
                    end do
                    
                case default
                    stop 'OpenAndReadChildHDF5Files - ModuleUpscaleHDF5 - ERR260'
                end select
            end do
        end do
        
        !deAllocates new instance
        if(associated(NewField%Values2D)) deallocate (NewField%Values2D)
        if(associated(NewField%Values3D)) deallocate (NewField%Values3D)
        deallocate (NewField)
        nullify    (NewField)
        
        ObjChildFile => Me%FirstChildFile
        do while (associated(ObjChildFile))
            if(associated(ObjChildFile%NewField%Values2D)) deallocate (ObjChildFile%NewField%Values2D)
            if(associated(ObjChildFile%NewField%Values3D)) deallocate (ObjChildFile%NewField%Values3D)
            deallocate(ObjChildFile%NewField)
            nullify(ObjChildFile%NewField)
            ObjChildFile => ObjChildFile%next
        end do

    end subroutine OpenAndReadChildHDF5Files

    !--------------------------------------------------------------------------
    
    subroutine set_module_time
        !External--------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, CurrentInstant
        !Begin---------------------------------------------------------------
        call GetHDF5GroupNumberOfItems(Me%Father%ObjHDF5, "/Time", &
                                       Me%Father%NumberOfInstants, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'set_module_time - ModuleUpscaleHDF5 - ERR10'

        allocate(Me%Father%InstantsArray(1:Me%Father%NumberOfInstants))
        
        !fill array with instants
        do CurrentInstant = 1, Me%Father%NumberOfInstants
            Me%Father%InstantsArray(CurrentInstant) = HDF5TimeInstant(CurrentInstant, Me%Father)
        end do

        Me%BeginTime = Me%Father%InstantsArray(1)
        Me%EndTime   = Me%Father%InstantsArray(Me%Father%NumberOfInstants)
        
        call StartComputeTime(Me%ObjTime, Me%BeginTime, Me%BeginTime, Me%EndTime,     & 
                              60., .false., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)'Error on StartComputeTime routine'
            stop 'set_module_time - ModuleUpscaleHDF5 - ERR20'
        end if
        
    end subroutine set_module_time
    
    !--------------------------------------------------------------------------
    
    subroutine get_timeinstants(Grid, StartInstant, EndInstant)
        !External--------------------------------------------------------------
        type(T_Grid), intent(INOUT)  :: Grid
        integer, intent(OUT)         :: StartInstant, EndInstant
        !Local-----------------------------------------------------------------
        integer                      :: STAT_CALL, CurrentInstant
        !Begin---------------------------------------------------------------
        
        call GetHDF5GroupNumberOfItems(Grid%ObjHDF5, "/Time", &
                                       Grid%NumberOfInstants, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'father_timeinstants - ModuleUpscaleHDF5 - ERR10'
        
        write(*,*)'Number of instants in available data: ', Grid%NumberOfInstants

        allocate(Grid%InstantsArray(1:Grid%NumberOfInstants))

        !fill array with instants
        do CurrentInstant = 1, Grid%NumberOfInstants
            Grid%InstantsArray(CurrentInstant) = HDF5TimeInstant(CurrentInstant, Grid)
        end do
        
        StartInstant = 1
        EndInstant = Grid%NumberOfInstants
        
    end subroutine get_timeinstants
    
    !--------------------------------------------------------------------------
    
    subroutine get_nproperties(Grid)
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_Grid)                                :: Grid
        logical                                     :: ExistGroup
        !Local-----------------------------------------------------------------
        !Begin------------------------------------------------------------------
        
        call GetHDF5GroupExist(Grid%ObjHDF5, "/Results", ExistGroup, &
                              nGroup = Grid%NumberOfProperties, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'get_father_nproperties - ModuleUpscaleHDF5 - ERR10'

        if(Grid%NumberOfProperties == 0)then
            write(*,*)'No data available in file: '//trim(Grid%InputFileName)
            stop 'get_father_nproperties - ModuleInterpolateGrids - ERR20'
        end if
        
    end subroutine get_nproperties
    
    !------------------------------------------------------------------------------------ 

    !--------------------------------------------------------------------------

    type(T_Time) function HDF5TimeInstant(Instant,ObjGridFile)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        type (T_Grid)                           :: ObjGridFile
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        !Local-----------------------------------------------------------------
        real, dimension(:), pointer             :: TimeVector
        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (ObjGridFile%ObjHDF5, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID     = ObjGridFile%ObjHDF5,               &
                             GroupName  = "/Time", Name = "Time",                   &
                             Array1D    = TimeVector, OutputNumber   = Instant, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleUpscaleHDF5 - ERR10'

        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))
        deallocate(TimeVector)
        nullify   (TimeVector)

    end function HDF5TimeInstant

    !--------------------------------------------------------------------------

    subroutine OutputInstants(CurrentInstant)
        !Arguments-------------------------------------------------------------
        integer                                         :: CurrentInstant
        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL
        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        type(T_Time)                                    :: CurrentDate
        !Begin-----------------------------------------------------------------
        
        CurrentDate = Me%Father%InstantsArray(CurrentInstant)

        call ExtractDate   (CurrentDate,                                            &
                            AuxTime(1), AuxTime(2), AuxTime(3),                     &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime

        call HDF5SetLimits  (Me%New%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputInstants - ModuleUpscaleHDF5 - ERR10'

        call HDF5WriteData  (Me%New%ObjHDF5, "/Time",                               &
                             "Time", "YYYY/MM/DD HH:MM:SS",                         &
                             Array1D = TimePtr,                                     &
                             OutputNumber = CurrentInstant, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputInstants - ModuleUpscaleHDF5 - ERR20'

    end subroutine OutputInstants

    !------------------------------------------------------------------------

    subroutine OutputFields(NewField, OutputNumber)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer, intent(IN)              :: NewField
        integer                                         :: OutputNumber
        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL
        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (Me%New%ObjHDF5, Me%Father%WorkSize2D%ILB, Me%Father%WorkSize2D%IUB, &
                             Me%Father%WorkSize2D%JLB, Me%Father%WorkSize2D%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleUpscaleHDF5 - ERR10'

        call HDF5WriteData(Me%New%ObjHDF5,                                  &
                           "/Results/"//NewField%Name,                      &
                           NewField%Name, NewField%Units,                   &
                           Array2D      = NewField%Values2D,                &
                           OutputNumber = OutputNumber, STAT= STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleUpscaleHDF5 - ERR20'

        !Writes everything to disk
        call HDF5FlushMemory (Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleUpscaleHDF5 - ERR30'

    end subroutine OutputFields
    
    !--------------------------------------------------------------------------
    
    subroutine OutputFields3D(NewField, OutputNumber)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer, intent(IN)              :: NewField
        integer                                         :: OutputNumber
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (Me%New%ObjHDF5,Me%Father%WorkSize3D%ILB, Me%Father%WorkSize3D%IUB,&
                             Me%Father%WorkSize3D%JLB, Me%Father%WorkSize3D%JUB, &
                             Me%Father%WorkSize3D%KLB, Me%Father%WorkSize3D%KUB, &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields3D - ModuleUpscaleHDF5 - ERR10'

        call HDF5WriteData(Me%New%ObjHDF5, "/Results/"//NewField%Name, &
                            NewField%Name, NewField%Units,             &
                            Array3D      = NewField%Values3D,          &
                            OutputNumber = OutputNumber,               &
                            STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields3D - ModuleUpscaleHDF5 - ERR20'

        !Writes everything to disk
        call HDF5FlushMemory (Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields3D - ModuleUpscaleHDF5 - ERR30'

    end subroutine OutputFields3D
    
    !--------------------------------------------------------------------------

    subroutine KillUpscaleHDF5

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL
        !Local-----------------------------------------------------------------
        integer                                         :: nUsers
        type (T_Child),     pointer                     :: ObjChildFile
        !Begin-----------------------------------------------------------------

        ObjChildFile => Me%FirstChildFile 

        do while (associated(ObjChildFile))
            deallocate(ObjChildFile%Info%NodeX)
            deallocate(ObjChildFile%Info%NodeY)
            deallocate(ObjChildFile%Info%NodeZ)

            call KillChildGrid(ObjChildFile)

            ObjChildFile => ObjChildFile%Next
        enddo

        call UnGetHorizontalGrid(Me%Father%ObjHorizontalGrid, Me%Father%ConnectionX,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillUpscaleHDF5 - ModuleUpscaleHDF5 - ERR20'

        call UnGetHorizontalGrid(Me%Father%ObjHorizontalGrid, Me%Father%ConnectionY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillUpscaleHDF5 - ModuleUpscaleHDF5 - ERR30'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillUpscaleHDF5 - ModuleUpscaleHDF5 - ERR40' 

        call KillHorizontalMap(Me%Father%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillUpscaleHDF5 - ModuleUpscaleHDF5 - ERR50'

        call KillGridData(Me%Father%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillUpscaleHDF5 - ModuleUpscaleHDF5 - ERR60'

        call KillHorizontalGrid(Me%Father%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillUpscaleHDF5 - ModuleUpscaleHDF5 - ERR70'
        
        call KillHDF5(Me%Father%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillUpscaleHDF5 - ModuleUpscaleHDF5 - ERR80'
        
        if(Me%Upscale3D) then
            call KillMap(Me%Father%ObjMap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillUpscaleHDF5 - ModuleUpscaleHDF5 - ERR90'

            call KillGeometry(Me%Father%ObjGeometry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillUpscaleHDF5 - ModuleUpscaleHDF5 - ERR100'
        endif
   
        deallocate(Me)
        nullify(Me)

    end subroutine KillUpscaleHDF5     

    !--------------------------------------------------------------------------
 
    subroutine KillChildGrid(ObjChildFile)
        !Arguments-------------------------------------------------------------
        type (T_Child),     pointer                :: ObjChildFile
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------
       
        write(*,*)'Killing child grid...'

        call UnGetHorizontalGrid(ObjChildFile%Info%ObjHorizontalGrid,                  & 
                                 ObjChildFile%Info%ConnectionX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillChildGrid - ModuleUpscaleHDF5 - ERR10'

        call UnGetHorizontalGrid(ObjChildFile%Info%ObjHorizontalGrid,                  &
                                 ObjChildFile%Info%ConnectionY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillChildGrid - ModuleUpscaleHDF5 - ERR20'

        call KillHorizontalMap(ObjChildFile%Info%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillChildGrid - ModuleUpscaleHDF5 - ERR30'

        call KillGridData(ObjChildFile%Info%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillChildGrid - ModuleUpscaleHDF5 - ERR40'

        call KillHDF5(ObjChildFile%Info%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillChildGrid - ModuleUpscaleHDF5 - ERR50'

        call KillHorizontalGrid(ObjChildFile%Info%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillChildGrid - ModuleUpscaleHDF5 - ERR60'
        
        if(Me%Upscale3D) then
            call KillMap(ObjChildFile%Info%ObjMap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillChildGrid - ModuleUpscaleHDF5 - ERR90'

            call KillGeometry(ObjChildFile%Info%ObjGeometry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillChildGrid - ModuleUpscaleHDF5 - ERR100'
        endif

    end subroutine KillChildGrid


    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !------------------------------------------------------------------------
      

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

end module ModuleUpscaleHDF5









