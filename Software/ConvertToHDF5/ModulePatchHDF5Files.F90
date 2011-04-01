!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid ConvertToHDF5
! MODULE        : PatchHDF5Files
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : August 2005
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Module to fill grid with several HDF5 files with specific grids/areas
!
!------------------------------------------------------------------------------
!DataFile
!
!   <begin_file>
!
!   ACTION                      : char                  [-]     !'PATCH HDF5 FILES' to use this module
!   TYPE_OF_INTERPOLATION       : int                   [-]     !Type of interpolation to translate grids
!                                                               !(only implemented 3 - Triangulation)
!   START                       : YYYY MM DD HH MM SS   [-]     !Start date of new file
!   END                         : YYYY MM DD HH MM SS   [-]     !End date of new file
!
!   OUTPUTFILENAME              : char                  [-]     !Path to HDF5 file to be created
!   NEW_GRID_FILENAME           : char                  [-]     !Path to grid file to be considered for HDF5
!                                                               !file to be created
!
!   <<begin_father>>                                            
!   LEVEL                       : int                   [-]     !Priority level of this father file (1= max.)
!   FATHER_FILENAME             : char                  [-]     !Path to father HDF5 file
!   FATHER_GRID_FILENAME        : char                  [-]     !Path to father grid file     
!   <<end_father>> 
!
!   <end_file>                                             

Module ModulePatchHDF5Files

    use ModuleGlobalData
    use ModuleTime
    use ModuleDrawing
    use ModuleHDF5,             only : GetHDF5FileAccess, ConstructHDF5, HDF5SetLimits, &
                                       HDF5WriteData, HDF5ReadData, HDF5FlushMemory,    &
                                       GetHDF5GroupNumberOfItems, KillHDF5, GetHDF5GroupID
    use ModuleEnterData,        only : GetData, ExtractBlockFromBlock, Block_Unlock,    &
                                       GetBlockSize, RewindBlock
    use ModuleHorizontalGrid,   only : ConstructHorizontalGrid, GetHorizontalGridSize,  & 
                                       GetGridLatitudeLongitude, GetCheckDistortion,    &
                                       ConstructFatherGridLocation,                     & 
                                       WriteHorizontalGrid, UnGetHorizontalGrid,        &
                                       KillHorizontalGrid    
    use ModuleGridData,         only : ConstructGridData, GetGridData, UngetGridData,   &
                                       KillGridData 
    use ModuleHorizontalMap,    only : ConstructHorizontalMap, GetWaterPoints2D,        &
                                       UngetHorizontalMap, KillHorizontalMap
    use ModuleTriangulation,    only : ConstructTriangulation, SetHeightValues,         & 
                                       InterPolation, KillTriangulation

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartPatchHDF5Files
    private ::      ReadOptions
    private ::          ReadFatherKeywords
    private ::              AddFatherFile
    private ::                  ConstructFatherFile
    private ::                  CopyFatherFile
    private ::      ConstructFatherGrid
    private ::      ConstructGridPolygon
    private ::          SetLimitsPolygon
    private ::      ConstructNewGrid
    private ::      ConstructRelevantFather
    private ::      FatherSonCommunication
    private ::      Open_HDF5_OutPut_File
    private ::      OpenAndReadFatherHDF5Files
    private ::          HDF5TimeInstant
    private ::          GetFatherTriangParam
    private ::          Triangulator
    private ::          OutputFields
    private ::      KillPatchHDF5Files
    private ::          KillFatherGrid

    !Selector                    
    
    !Modifier

    !Destructor

    !Management
    
    !Interfaces----------------------------------------------------------------

    !Parameters----------------------------------------------------------------
    integer, parameter                                      :: Triangulation = 3

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
        character(len=PathLength)                           :: FileName
        character(len=PathLength)                           :: GridFileName
        type(T_Time), dimension(:), pointer                 :: InstantsArray
        integer                                             :: NumberOfInstants
        integer                                             :: NumberOfProperties
        type(T_Size2D)                                      :: Size2D
        type(T_Size2D)                                      :: WorkSize2D
        integer, dimension(:,:),    pointer                 :: WaterPoints2D
        real,    dimension(:,:),    pointer                 :: ConnectionX, ConnectionY
        integer                                             :: Level
        type(T_Polygon), pointer                            :: Polygon
        integer                                             :: NRemoveFrame
    end type  T_Grid

    type       T_Father
        type(T_Grid)                                        :: Info
        real,    dimension(:  ), pointer                    :: NodeX, NodeY, NodeZ
        integer                                             :: NumberOfNodes
        integer                                             :: ObjTriangulation     = 0
        type(T_Field), pointer                              :: NewField
        type(T_Father), pointer                             :: Next
    end type  T_Father
    
    private :: T_PatchHDF5Files
    type       T_PatchHDF5Files
        integer                                             :: ObjEnterData         = 0
        integer                                             :: ObjTime              = 0
        integer                                             :: TypeOfInterpolation
        type(T_Grid )                                       :: New
        type(T_Time)                                        :: BeginTime, EndTime
        logical                                             :: TimeWindow = .true.
        real,    dimension(:,:), pointer                    :: SonCenterX, SonCenterY
        type(T_Father), pointer                             :: FirstFatherFile
        integer, dimension(:,:), pointer                    :: RelevantFather
        character(len=StringLength), dimension(:), pointer  :: FieldsToInterpolate
        integer                                             :: nFieldsToInterpolate
        logical                                             :: ConvertAllFields
    end type  T_PatchHDF5Files

    !Global Module Variables
    type (T_PatchHDF5Files), pointer                        :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartPatchHDF5Files(ObjPatchHDF5FilesID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: ObjPatchHDF5FilesID
        integer,           intent(IN )                  :: ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !External----------------------------------------------------------------
        
        
        !Local-------------------------------------------------------------------
        type (T_Father),     pointer                    :: ObjFatherFile

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, ObjPatchHDF5FilesID)

        !Read instructions from user
        call ReadOptions(ClientNumber)

        !Cycle all the fathers to construct grid
        ObjFatherFile => Me%FirstFatherFile

        do while (associated(ObjFatherFile))
            call ConstructFatherGrid(ObjFatherFile)

            call ConstructGridPolygon(ObjFatherFile%Info)

            ObjFatherFile => ObjFatherFile%Next
        enddo

        call ConstructNewGrid
 
        !Construct array of relevant father files (one for each cell)
        call ConstructRelevantFather 

        !Construct relation between the fathers and the new grid
        call FatherSonCommunication(Me%New)

        !Open Output File (construct header data)
        call Open_HDF5_OutPut_File

        !Read data from father files
        call OpenAndReadFatherHDF5Files

        call KillPatchHDF5Files

        STAT = SUCCESS_

        !----------------------------------------------------------------------

    end subroutine StartPatchHDF5Files

    !------------------------------------------------------------------------

    subroutine ReadOptions(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer, intent(in)                         :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag, iflag1

        !Begin-----------------------------------------------------------------
        
        write(*,*)'Reading instructions...'

        call GetData(Me%New%FileName,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUTFILENAME',                   &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModulePatchHDF5Files - ERR10'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModulePatchHDF5Files - ERR20'
        end if

        call GetData(Me%New%GridFileName,                               &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'NEW_GRID_FILENAME',                &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModulePatchHDF5Files - ERR30'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModulePatchHDF5Files - ERR40'
        end if

        call GetData(Me%TypeOfInterpolation,                            &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'TYPE_OF_INTERPOLATION',            &
                     Default      = Triangulation,                      &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModulePatchHDF5Files - ERR50'

        if (iflag == 0)then
            write(*,*)'Must specify type of interpolation'
            stop 'ReadOptions - ModulePatchHDF5Files - ERR60'
        end if

        if(Me%TypeOfInterpolation == Triangulation)then
            write(*,*)
            write(*,*)'Type of interpolation : Triangulation'
            write(*,*)
            !(to start only triangulation is considered)
        else
            write(*,*) 'Unknown type of interpolation'
            stop       'ReadOptions - ModulePatchHDF5Files - ERR70' 
        end if

        call GetData(Me%BeginTime,                                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'START',                            &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModulePatchHDF5Files - ERR80'

        call GetData(Me%EndTime,                                        &
                     Me%ObjEnterData, iflag1,                           &
                     SearchType   = FromBlock,                          &
                     keyword      = 'END',                              &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModulePatchHDF5Files - ERR90'

        if (iflag==0 .and. iflag1==0) Me%TimeWindow = .FALSE.

        call StartComputeTime(Me%ObjTime, Me%BeginTime, Me%BeginTime, Me%EndTime,     & 
                              60., .false., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModulePatchHDF5Files - ERR100'

        call ReadFatherKeywords (ClientNumber) 

        call ReadFieldsToConvert(ClientNumber)
        
        !(to start only 2D interpolation is considered)

    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    subroutine ReadFieldsToConvert(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer, intent(in)                 :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL
        integer                             :: StartLine, EndLine, Count
        integer                             :: CurrentLineNumber
        logical                             :: BlockFound
        character(len=StringLength)         :: PropertyName

        !Begin-----------------------------------------------------------------
        
        call RewindBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModulePatchHDF5Files - ERR00'

        call ExtractBlockFromBlock(Me%ObjEnterData,             &
                                    ClientNumber,               &
                                    '<<BeginFields>>',          &
                                    '<<EndFields>>',            &
                                    BlockFound,                 &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModulePatchHDF5Files - ERR01'

        if(BlockFound)then
    
            Me%ConvertAllFields = .false.
    
            call GetBlockSize(Me%ObjEnterData,                      &
                              ClientNumber,                         &
                              StartLine,                            &
                              EndLine,                              &
                              FromBlockInBlock,                     &
                              STAT = STAT_CALL)

            Me%nFieldsToInterpolate = EndLine - StartLine - 1

            allocate(Me%FieldsToInterpolate(1:Me%nFieldsToInterpolate))
        
            Count = 1

            do CurrentLineNumber = StartLine + 1 , EndLine - 1

                call GetData(PropertyName,                          &
                             Me%ObjEnterData,                       &
                             flag,                                  &
                             SearchType  = FromBlock_,              &
                             Buffer_Line = CurrentLineNumber,       &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModulePatchHDF5Files - ERR02'

                Me%FieldsToInterpolate(Count) = PropertyName

                Count = Count + 1

            end do

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModulePatchHDF5Files - ERR03'
        
        else

            Me%ConvertAllFields = .true. 

        endif

    end subroutine ReadFieldsToConvert

    !--------------------------------------------------------------------------

    subroutine ReadFatherKeywords(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer                                     :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: BlockFound
        logical                                     :: AtLeastOneBlock = .false.
        integer                                     :: NumberFathers = 0

        !Begin-----------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBlock(Me%ObjEnterData,                         &
                                        ClientNumber,                           &
                                        '<<begin_father>>',                     &
                                        '<<end_father>>',                       &
                                        BlockFound,                             &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                    
                    AtLeastOneBlock = .true.

                    call AddFatherFile 

                    NumberFathers = NumberFathers + 1

                else cd2

                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ReadFatherKeywords - ModulePatchHDF5Files - ERR10'

                    exit do1

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ReadFatherKeywords - ModulePatchHDF5Files - ERR20'
            else cd1
                stop 'ReadFatherKeywords - ModulePatchHDF5Files - ERR30'
            end if cd1

        end do do1

        if (.not. AtLeastOneBlock) then                                            
            write(*,*) 'No father file block is indicated in input file. '
            stop 'ReadFatherKeywords - ModulePatchHDF5Files - ERR40'
        end if

        !If only one block then program breaks: better to use InterpolateGrids
        if (NumberFathers < 2) then                                            
            write(*,*) 'At least two father file blocks needed in input file. '
            write(*,*) 'Only one father file: use INTERPOLATE GRIDS action instead. '
            stop 'ReadFatherKeywords - ModulePatchHDF5Files - ERR50'
        end if

    end subroutine ReadFatherKeywords
    
    !--------------------------------------------------------------------------

    subroutine AddFatherFileOld(ObjFatherFile)

        !Arguments-------------------------------------------------------------
        type (T_Father),     pointer                  :: ObjFatherFile

        !Local-----------------------------------------------------------------
        type (T_Father),     pointer                  :: PreviousFatherFile
        type (T_Father),     pointer                  :: NewFatherFile

        !Begin-----------------------------------------------------------------

        !Allocates new HDF5File
        allocate (NewFatherFile)
        nullify  (NewFatherFile%Next)

        !Insert new file into list and makes current ?? point to it
        if (.not. associated(Me%FirstFatherFile)) then
            Me%FirstFatherFile         => NewFatherFile
            ObjFatherFile              => NewFatherFile
        else
            PreviousFatherFile         => Me%FirstFatherFile
            ObjFatherFile              => Me%FirstFatherFile%Next
            do while (associated(ObjFatherFile))
                PreviousFatherFile     => ObjFatherFile
                ObjFatherFile          => ObjFatherFile%Next
            enddo
            ObjFatherFile              => NewFatherFile
            PreviousFatherFile%Next    => NewFatherFile
        end if

    end subroutine AddFatherFileOld

    !--------------------------------------------------------------------------

    subroutine ConstructFatherFile(NewFatherFile)

        !Arguments-------------------------------------------------------------
        type (T_Father),      pointer           :: NewFatherFile

        !External--------------------------------------------------------------
        integer                                   :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        call GetData(NewFatherFile%Info%FileName,                       &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'FATHER_FILENAME',                  &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructFatherFile - ModulePatchHDF5Files - ERR10'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ConstructFatherFile - ModulePatchHDF5Files - ERR20'
        end if

        call GetData(NewFatherFile%Info%GridFileName,                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'FATHER_GRID_FILENAME',             &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructFatherFile - ModulePatchHDF5Files - ERR30'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ConstructFatherFile - ModulePatchHDF5Files - ERR40'
        end if

        call GetData(NewFatherFile%Info%Level,                          &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'LEVEL',                            &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructFatherFile - ModulePatchHDF5Files - ERR50'
        
        if (iflag == 0)then
            write(*,*)'Must specify priority level of file to convert'
            stop 'ConstructFatherFile - ModulePatchHDF5Files - ERR60'
        end if          

        call GetData(NewFatherFile%Info%NRemoveFrame,                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'N_REMOVE_FRAME',                   &
                     default      = 0,                                  &                        
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructFatherFile - ModulePatchHDF5Files - ERR70'

        
    end subroutine ConstructFatherFile

    !--------------------------------------------------------------------------

    subroutine AddFatherFile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Father),     pointer                  :: PreviousFatherFile, LastFatherFile
        type (T_Father),     pointer                  :: NewFatherFile, FatherFileAux

        !Begin-----------------------------------------------------------------

        !Allocates new father file
        allocate (NewFatherFile)
        nullify  (NewFatherFile%Next)
        call ConstructFatherFile(NewFatherFile)

        !Put Father files in list by order of priority
        if (.not. associated(Me%FirstFatherFile)) then

            call CopyFatherFile(Me%FirstFatherFile, NewFatherFile)
            deallocate(Me%FirstFatherFile%Next) 
            nullify(Me%FirstFatherFile%Next)

        else

            if (NewFatherFile%Info%Level < Me%FirstFatherFile%Info%Level) then
                !current file should be the first file in list

                !save the previous list 
                allocate(FatherFileAux)
                call CopyFatherFile(FatherFileAux, Me%FirstFatherFile)
                FatherFileAux%Next => Me%FirstFatherFile%Next              

                !make the first element in the list of relevant files equal to current file
                call CopyFatherFile(Me%FirstFatherFile, NewFatherFile)
                Me%FirstFatherFile%Next => FatherFileAux

            elseif (NewFatherFile%Info%Level == Me%FirstFatherFile%Info%Level) then
            
               write(*,*)'Two father files with same priority level:'
               write(*,*) trim(NewFatherFile%Info%FileName)
               write(*,*) trim(Me%FirstFatherFile%Info%FileName)
               stop 'AddFatherFile - ModulePatchHDF5Files - ERR10'
            
            else 
                !check next files in list

                !locate previous file in the first file
                allocate(PreviousFatherFile)
                PreviousFatherFile => Me%FirstFatherFile                   

                do while(associated(PreviousFatherFile))

                    if (.not. associated(PreviousFatherFile%Next)) then
        
                        !current file is the last file in the list of relevant files
                        call CopyFatherFile(PreviousFatherFile%Next, NewFatherFile)

                        !current file was added to list
                        exit

                    else

                        !check if current file should be located before the next file
                        if (NewFatherFile%Info%Level < PreviousFatherFile%Info%Level) then
                            
                            !current file should be located before next file

                            !save the previous list begining in PreviousHDF5File%Next
                            allocate(LastFatherFile)
                            LastFatherFile => PreviousFatherFile%Next
                            allocate(FatherFileAux)

                            call CopyFatherFile(FatherFileAux, LastFatherFile)
                            FatherFileAux%Next => LastFatherFile%Next

                            call CopyFatherFile(LastFatherFile, NewFatherFile)

                            PreviousFatherFile%Next => LastFatherFile
                            LastFatherFile%Next => FatherFileAux

                            !current file was added to list
                            exit

                        elseif (NewFatherFile%Info%Level == PreviousFatherFile%Info%Level) then

                           write(*,*)'Two father files with same priority level:'
                           write(*,*) trim(NewFatherFile%Info%FileName)
                           write(*,*) trim(PreviousFatherFile%Info%FileName)
                           stop 'AddFatherFile - ModulePatchHDF5Files - ERR20'

                        end if

                    end if

                    !check next file in list
                    PreviousFatherFile => PreviousFatherFile%Next     

                end do

            end if

            nullify(NewFatherFile)

        end if

    end subroutine AddFatherFile

    !--------------------------------------------------------------------------

    subroutine CopyFatherFile(FatherFileNew, FatherFileX)

        !Arguments-------------------------------------------------------------
        type(T_Father), pointer                     :: FatherFileNew
        type(T_Father), pointer                     :: FatherFileX

        !Local-----------------------------------------------------------------
      
        !Begin-----------------------------------------------------------------

        !Copies the values of the fields of a FatherFile to another FatherFile
        allocate(FatherFileNew)

        FatherFileNew%Info%FileName          =  FatherFileX%Info%FileName
        FatherFileNew%Info%GridFileName      =  FatherFileX%Info%GridFileName
        FatherFileNew%Info%Level             =  FatherFileX%Info%Level
        FatherFileNew%Info%NRemoveFrame      =  FatherFileX%Info%NRemoveFrame

    end subroutine CopyFatherFile

    !--------------------------------------------------------------------------

    subroutine ConstructFatherGrid(ObjFatherFile)

        !Arguments-------------------------------------------------------------
        type (T_Father),      pointer               :: ObjFatherFile
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        logical                                     :: exist

        !Begin-----------------------------------------------------------------
       
        write(*,*)'Constructing father grid...'

        !Verifies if file exists
        inquire(FILE = ObjFatherFile%Info%GridFileName, EXIST = exist)
        if (.not. exist) then
            write(*,*)'Grid file does not exist'
            stop 'ConstructFatherGrid - ModulePatchHDF5Files - ERR10'
        endif

        call ConstructHorizontalGrid(ObjFatherFile%Info%ObjHorizontalGrid,              & 
                                     ObjFatherFile%Info%GridFileName, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                                     &
            stop 'ConstructFatherGrid - ModulePatchHDF5Files - ERR20'

        call GetHorizontalGridSize(ObjFatherFile%Info%ObjHorizontalGrid,                &
                                   WorkSize = ObjFatherFile%Info%WorkSize2D,            &
                                   Size     = ObjFatherFile%Info%Size2D,                &
                                   STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                                     &
            stop 'ConstructFatherGrid -  ModulePatchHDF5Files - ERR30'


        call ConstructGridData(GridDataID       = ObjFatherFile%Info%ObjBathymetry,     &
                               HorizontalGridID = ObjFatherFile%Info%ObjHorizontalGrid, &
                               FileName         = ObjFatherFile%Info%GridFileName,      &
                               STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                                     &
            stop 'ConstructFatherGrid -  ModulePatchHDF5Files - ERR40'

        call ConstructHorizontalMap(HorizontalMapID  = ObjFatherFile%Info%ObjHorizontalMap,     &
                                    GridDataID       = ObjFatherFile%Info%ObjBathymetry,        &
                                    HorizontalGridID = ObjFatherFile%Info%ObjHorizontalGrid,    &
                                    ActualTime       = Me%BeginTime,                            & 
                                    STAT             = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_)                                                     & 
            stop 'ConstructFatherGrid -  ModulePatchHDF5Files - ERR50'


    end subroutine ConstructFatherGrid

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
        if(STAT_CALL .ne. SUCCESS_)                                                 & 
            stop 'ConstructGridPolygon -  ModulePatchHDF5Files - ERR10'

        allocate(ObjGrid%Polygon)
        nullify(ObjGrid%Polygon%Next)

        if (ObjGrid%NRemoveFrame > 0) &
            write(*,*) 'Removing', ObjGrid%NRemoveFrame, ' boundary points'

        ILB = ObjGrid%WorkSize2D%ILB + ObjGrid%NRemoveFrame
        JLB = ObjGrid%WorkSize2D%JLB + ObjGrid%NRemoveFrame
        IUB = ObjGrid%WorkSize2D%IUB - ObjGrid%NRemoveFrame
        JUB = ObjGrid%WorkSize2D%JUB - ObjGrid%NRemoveFrame

        !Construct the vertixes of polygon (should be in center of outmost cells):
        !get polygon number of vertixes
        ObjGrid%Polygon%Count = 2 * (IUB - ILB + 1) + 2 * (JUB - JLB + 1 - 2) + 1

        allocate(ObjGrid%Polygon%VerticesF(1:ObjGrid%Polygon%Count))

        CurrentVertix = 1
        
        j = JLB

        do i = ILB, IUB
            
            ObjGrid%Polygon%VerticesF(CurrentVertix)%X  = GridLongitude(i,j)
            ObjGrid%Polygon%VerticesF(CurrentVertix)%Y  = GridLatitude(i,j)

            CurrentVertix = CurrentVertix + 1

        end do

        i = IUB

        do j = JLB + 1, JUB-1
            
            ObjGrid%Polygon%VerticesF(CurrentVertix)%X  = GridLongitude(i,j)
            ObjGrid%Polygon%VerticesF(CurrentVertix)%Y  = GridLatitude(i,j)

            CurrentVertix = CurrentVertix + 1

        end do

        j = JUB

        do i = ILB, IUB
            
            ObjGrid%Polygon%VerticesF(CurrentVertix)%X  =                           &
                GridLongitude(IUB+ILB-i,j)
            ObjGrid%Polygon%VerticesF(CurrentVertix)%Y  =                           &
                GridLatitude(IUB+ILB-i,j)

            CurrentVertix = CurrentVertix + 1

        end do

        i = ILB

        do j = JLB + 1, JUB-1
            
            ObjGrid%Polygon%VerticesF(CurrentVertix)%X  =                           &
                GridLongitude(i,JUB+JLB-j)
            ObjGrid%Polygon%VerticesF(CurrentVertix)%Y  =                           &
                GridLatitude(i,JUB+JLB-j)

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

        !Begin-----------------------------------------------------------------

        write(*,*)'Constructing new grid...'
       
        !Verifies if file exists
        inquire(FILE = Me%New%GridFileName, EXIST = exist)
        if (.not. exist) then
            write(*,*)'Grid file does not exist'
            stop 'ConstructNewGrid - ModulePatchHDF5Files - ERR10'
        endif

        call ConstructHorizontalGrid(Me%New%ObjHorizontalGrid, Me%New%GridFileName, & 
                                     STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid - ModulePatchHDF5Files - ERR20'


        call GetHorizontalGridSize  (Me%New%ObjHorizontalGrid,                      &
                                     WorkSize = Me%New%WorkSize2D,                  &
                                     Size     = Me%New%Size2D,                      &
                                     STAT     = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModulePatchHDF5Files - ERR30'

        call ConstructGridData      (GridDataID       = Me%New%ObjBathymetry,        &
                                     HorizontalGridID = Me%New%ObjHorizontalGrid,    &
                                     FileName         = Me%New%GridFileName,         &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModulePatchHDF5Files - ERR40'

        call ConstructHorizontalMap (HorizontalMapID  = Me%New%ObjHorizontalMap,     &
                                     GridDataID       = Me%New%ObjBathymetry,        &
                                     HorizontalGridID = Me%New%ObjHorizontalGrid,    &
                                     ActualTime       = Me%BeginTime,                & 
                                     STAT             = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModulePatchHDF5Files - ERR50'


    end subroutine ConstructNewGrid

    !------------------------------------------------------------------------

    subroutine ConstructRelevantFather

        !Arguments-------------------------------------------------------------
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer               :: GridLongitude, GridLatitude
        integer                                     :: i, j
        type (T_PointF),   pointer                  :: Point
        type (T_Father),     pointer                :: ObjFatherFile

        !Begin-----------------------------------------------------------------

        !Obtain coordinates of cells center of new grid:
        call GetGridLatitudeLongitude(Me%New%ObjHorizontalGrid,                     &
                                      GridLatitude = GridLatitude,                  & 
                                      GridLongitude = GridLongitude, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                                 & 
            stop 'ConstructRelevantFather -  ModulePatchHDF5Files - ERR10'

        !Allocate array of relevant fathers:
        allocate(Me%RelevantFather(Me%New%WorkSize2D%ILB:Me%New%WorkSize2D%IUB,     &
                                   Me%New%WorkSize2D%JLB:Me%New%WorkSize2D%JUB))

        allocate(Point)

        !Cycle every cell of the new grid and check which is the highest priority 
        !father polygon where it is located.
        do j = Me%New%WorkSize2D%JLB, Me%New%WorkSize2D%JUB
        do i = Me%New%WorkSize2D%ILB, Me%New%WorkSize2D%IUB

            Point%X = GridLongitude(i,j)
            Point%Y = GridLatitude(i,j)

            ObjFatherFile => Me%FirstFatherFile

do2 :       do while (associated(ObjFatherFile))

                if (IsPointInsidePolygon(Point, ObjFatherFile%Info%Polygon)) then
                    Me%RelevantFather(i,j) = ObjFatherFile%Info%Level
                    !(this because father list is ordered by decreasing priority)
                    exit do2
                else
                    ObjFatherFile => ObjFatherFile%Next
                endif
            enddo do2

        enddo
        enddo

        call UnGetHorizontalGrid(Me%New%ObjHorizontalGrid,                              &
                                 GridLatitude, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                                     & 
            stop 'ConstructRelevantFather -  ModulePatchHDF5Files - ERR20'

        call UnGetHorizontalGrid(Me%New%ObjHorizontalGrid,                              &
                                 GridLongitude, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                                     & 
            stop 'ConstructRelevantFather -  ModulePatchHDF5Files - ERR30'
        
    end subroutine ConstructRelevantFather

    !--------------------------------------------------------------------------
    
    subroutine FatherSonCommunication(NewGrid)

        !Arguments-------------------------------------------------------------
        type (T_Grid)                               :: NewGrid

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        type (T_Father),     pointer                :: ObjFatherFile

        !Begin-----------------------------------------------------------------

        write(*,*)'Constructing communication between grids...'

        !Cycle all the fathers to construct grid
        ObjFatherFile => Me%FirstFatherFile

        do while (associated(ObjFatherFile))

            call GetGridLatitudeLongitude(ObjFatherFile%Info%ObjHorizontalGrid,                  &
                                          GridLatitudeConn  = ObjFatherFile%Info%ConnectionY,    &
                                          GridLongitudeConn = ObjFatherFile%Info%ConnectionX,    &
                                          STAT  = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)                                                 &
                stop 'FatherSonCommunication - ModulePatchHDF5Files - ERR30'            

            ObjFatherFile => ObjFatherFile%Next
        enddo

        call GetGridLatitudeLongitude(NewGrid%ObjHorizontalGrid,                        &
                                      GridLatitudeConn  = NewGrid%ConnectionY,          &
                                      GridLongitudeConn = NewGrid%ConnectionX,          &
                                      STAT  = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                                     &
            stop 'FatherSonCommunication - ModulePatchHDF5Files - ERR40'


    end subroutine FatherSonCommunication

    !------------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: HDF5_CREATE
        real,       dimension(:,:  ), pointer       :: Bathymetry
        integer,    dimension(:,:,:), pointer       :: WaterPoints3D
        integer,    dimension(:,:  ), pointer       :: WaterPoints2D

        !----------------------------------------------------------------------

        allocate(WaterPoints3D(Me%New%Size2D%ILB:Me%New%Size2D%IUB,                     &
                               Me%New%Size2D%JLB:Me%New%Size2D%JUB, 1))


        call GetWaterPoints2D(Me%New%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR20'
        
        WaterPoints3D(Me%New%WorkSize2D%ILB:Me%New%WorkSize2D%IUB,                      &
                      Me%New%WorkSize2D%JLB:Me%New%WorkSize2D%JUB,1) =                  &
            WaterPoints2D(Me%New%WorkSize2D%ILB:Me%New%WorkSize2D%IUB,                  & 
                          Me%New%WorkSize2D%JLB:Me%New%WorkSize2D%JUB)

        
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        
        !Opens HDF5 File
        call ConstructHDF5(Me%New%ObjHDF5, Me%New%FileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR30'
        
        
        call HDF5SetLimits  (Me%New%ObjHDF5, Me%New%WorkSize2D%ILB,                     & 
                             Me%New%WorkSize2D%IUB,Me%New%WorkSize2D%JLB,               & 
                             Me%New%WorkSize2D%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR40'

        call GetGridData(Me%New%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR50'
        
        call HDF5WriteData   (Me%New%ObjHDF5, "/Grid", "Bathymetry", "-",               &
                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR60'            

        call UngetGridData(Me%New%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR70'

        call WriteHorizontalGrid (Me%New%ObjHorizontalGrid, Me%New%ObjHDF5,             & 
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR80'

        call HDF5SetLimits  (Me%New%ObjHDF5, Me%New%WorkSize2D%ILB,                     & 
                             Me%New%WorkSize2D%IUB, Me%New%WorkSize2D%JLB,              &
                             Me%New%WorkSize2D%JUB, 1,1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR90'

        call HDF5WriteData   (Me%New%ObjHDF5, "/Grid", "WaterPoints", "-",              &
                              Array3D = WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR100'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%New%ObjHDF5, STAT = STAT_CALL)                         
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR110'

        call UngetHorizontalMap(Me%New%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR120'

        deallocate(WaterPoints3D)
        nullify   (WaterPoints3D)

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

    subroutine OpenAndReadFatherHDF5Files

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        type(T_Grid ), pointer                      :: AuxGrid
        integer                                     :: CurrentInstant
        integer                                     :: StartInstant, EndInstant
        integer                                     :: AuxStartInstant, AuxEndInstant
        integer                                     :: CurrentProperty 
        type(T_Field), pointer                      :: NewField, AuxField
        character(len=StringLength)                 :: PropertyName
        integer                                     :: Rank
        integer, dimension(7)                       :: Dimensions
        integer                                     :: Count, Count2, NewCurrentInstant
        integer                                     :: i, j, n
        integer, dimension(:,:  ), pointer          :: FatherWP2D, AuxWP2D
        type (T_Father),     pointer                :: ObjFatherFile
        logical                                     :: FirstFather
        integer                                     :: AuxNumberOfProperties
        logical                                     :: exist, ConvertThisField, FirstTime
        integer                                     :: HDF5_READ

        !Begin-----------------------------------------------------------------
        FirstTime   = .true.
        FirstFather = .true.

        Count = 1

        nullify(Me%FirstFatherFile%NewField, NewField, AuxField)

        !Initial check of father files
        !Cycle every father file: --> afterwards implement to check only relevant fathers
        ObjFatherFile => Me%FirstFatherFile

        do while (associated(ObjFatherFile))

            write(*,*)'Checking father file:', trim(ObjFatherFile%Info%FileName)

            !Verifies if file exists
            inquire(FILE = ObjFatherFile%Info%FileName, EXIST = exist)
            if (.not. exist) then
                write(*,*)'HDF5 file does not exist.'
                stop 'OpenAndReadHDF5Files - ModulePatchHDF5Files - ERR10'
            endif

            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

            !Open HDF5 file
            call ConstructHDF5 (ObjFatherFile%Info%ObjHDF5,                         & 
                trim(ObjFatherFile%Info%FileName), HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            & 
                stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR20'

            !Get number of instants
            call GetHDF5GroupNumberOfItems(ObjFatherFile%Info%ObjHDF5, "/Time",     &
                                           ObjFatherFile%Info%NumberOfInstants,     &
                                           STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            & 
                stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR30'

            write(*,*)'Number of instants in available data: ',                     &
                ObjFatherFile%Info%NumberOfInstants

            allocate(ObjFatherFile%Info%InstantsArray(1:                            &
                ObjFatherFile%Info%NumberOfInstants))

            !fill array with instants
            do CurrentInstant = 1, ObjFatherFile%Info%NumberOfInstants
                ObjFatherFile%Info%InstantsArray(CurrentInstant) =                  &
                    HDF5TimeInstant(CurrentInstant,ObjFatherFile)
            end do
       

            !check time window
            if (Me%TimeWindow) then
                if(ObjFatherFile%Info%InstantsArray(1) .gt. Me%BeginTime)then
                    write(*,*)'Data available starts after speficied date.'
                    stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR40'
                end if

                if(ObjFatherFile%Info%InstantsArray(ObjFatherFile%Info%NumberOfInstants)    &
                    .lt. Me%EndTime)then
                    write(*,*)'Data available ends before speficied date.'
                    stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR50'
                end if
            else
                Me%BeginTime = ObjFatherFile%Info%InstantsArray(1)
                Me%EndTime   = ObjFatherFile%Info%InstantsArray(ObjFatherFile%Info%NumberOfInstants)
            endif


            !select time window begin
            do CurrentInstant = 1, ObjFatherFile%Info%NumberOfInstants
                if(ObjFatherFile%Info%InstantsArray(CurrentInstant) .ge.            &
                    Me%BeginTime)then

                    AuxStartInstant = CurrentInstant
                    exit

                end if
            end do

            if (FirstFather) then
                StartInstant = AuxStartInstant
            else
                !Check the time window begin against the previous father
                if (AuxStartInstant .ne. StartInstant) then
                    write(*,*)'Father file has not start time consistent with previous father'
                    stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR60'
                endif
            endif
    
            !select time window end
            do CurrentInstant = StartInstant, ObjFatherFile%Info%NumberOfInstants
                if(ObjFatherFile%Info%InstantsArray(CurrentInstant) .ge.            & 
                    Me%EndTime)then

                    AuxEndInstant = CurrentInstant
                    exit 

                end if
            end do

            if (FirstFather) then
                EndInstant = AuxEndInstant
            else
                !Check the time window end against the previous father
                if (AuxEndInstant .ne. EndInstant) then
                    write(*,*)'Father file has not end time consistent with previous father'
                    stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR70'
                endif
            endif          

            !(Assumed an equal DT in all fathers)

            !check number of properties in file
            call GetHDF5GroupNumberOfItems(ObjFatherFile%Info%ObjHDF5, "/Results",  &
                                          ObjFatherFile%Info%NumberOfProperties,    &
                                          STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            & 
                stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR80'

            if(ObjFatherFile%Info%NumberOfProperties == 0)then
                write(*,*)'No data available in file: '                             &
                    //trim(ObjFatherFile%Info%FileName)
                stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR90'
            end if

            if ( .not. FirstFather) then
                !Check number of properties against the previous father
                if (AuxNumberOfProperties .ne.                                      & 
                    ObjFatherFile%Info%NumberOfProperties) then
                    write(*,*)'Father file has not number of properties'
                    write(*,*)'consistent with previous father.'
                    stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR100'
                endif
            endif          
            AuxNumberOfProperties = ObjFatherFile%Info%NumberOfProperties

            allocate(ObjFatherFile%Info%WaterPoints2D(ObjFatherFile%Info%Size2D%ILB:    &
                                         ObjFatherFile%Info%Size2D%IUB,                 &
                                         ObjFatherFile%Info%Size2D%JLB:                 &
                                         ObjFatherFile%Info%Size2D%JUB))

            ObjFatherFile%Info%WaterPoints2D(:,:) = 0

            !Construct triangulation parameters for father
            call GetWaterPoints2D(ObjFatherFile%Info%ObjHorizontalMap,              & 
                 FatherWP2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR110'

            ObjFatherFile%Info%WaterPoints2D(:,:) = FatherWP2D(:,:)

            call UnGetHorizontalMap(ObjFatherFile%Info%ObjHorizontalMap,            & 
                 FatherWP2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR120'
            
            if(Me%TypeOfInterpolation == Triangulation) then
                call GetFatherTriangParam(ObjFatherFile)

                !Constructs Triangulation
                call ConstructTriangulation (ObjFatherFile%ObjTriangulation,        &
                                             ObjFatherFile%NumberOfNodes,           &
                                             ObjFatherFile%NodeX,                   &
                                             ObjFatherFile%NodeY,                   &
                                             STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                          & 
                    stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR130'
            endif

            if (FirstFather) FirstFather = .false.

            ObjFatherFile => ObjFatherFile%Next
        enddo

        Me%New%NumberOfInstants = EndInstant - StartInstant + 1

        allocate(Me%New%InstantsArray(1:Me%New%NumberOfInstants))

        Me%New%InstantsArray =                                                      & 
            Me%FirstFatherFile%Info%InstantsArray(StartInstant:EndInstant)

        Me%New%NumberOfProperties = Me%FirstFatherFile%Info%NumberOfProperties

        AuxGrid => Me%New

        if(Me%TypeOfInterpolation == Triangulation)then

            call GetWaterPoints2D(AuxGrid%ObjHorizontalMap,                         & 
                                  AuxWP2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR140'

            allocate(Me%SonCenterX(AuxGrid%Size2D%ILB:AuxGrid%Size2D%IUB,           &
                                   AuxGrid%Size2D%JLB:AuxGrid%Size2D%JUB))

            allocate(Me%SonCenterY(AuxGrid%Size2D%ILB:AuxGrid%Size2D%IUB,           &
                                   AuxGrid%Size2D%JLB:AuxGrid%Size2D%JUB))

            do j = AuxGrid%WorkSize2D%JLB, AuxGrid%WorkSize2D%JUB
            do i = AuxGrid%WorkSize2D%ILB, AuxGrid%WorkSize2D%IUB
                
                !Find Son cell center
                Me%SonCenterX(i,j) = ((AuxGrid%ConnectionX(i, j  ) +                &
                                       AuxGrid%ConnectionX(i+1, j  ))/2. +          &
                                      (AuxGrid%ConnectionX(i, j+1) +                &
                                       AuxGrid%ConnectionX(i+1, j+1))/2.)/2.
    
                Me%SonCenterY(i,j) = ((AuxGrid%ConnectionY(i, j  ) +                &
                                       AuxGrid%ConnectionY(i+1, j  ))/2. +          &
                                      (AuxGrid%ConnectionY(i, j+1) +                &
                                       AuxGrid%ConnectionY(i+1, j+1))/2.)/2.

            enddo
            enddo

        endif

prop:   do CurrentProperty = 1, Me%New%NumberOfProperties


            !(Assuming that every father has same properties then read name from first)

            !get property name
            call GetHDF5GroupID(Me%FirstFatherFile%Info%ObjHDF5, "/Results",        &
                                CurrentProperty, PropertyName,                      &
                                STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR150'


            if(.not. Me%ConvertAllFields)then

                ConvertThisField = .false.

                do n = 1, Me%nFieldsToInterpolate

                    if(Me%FieldsToInterpolate(n) == PropertyName)then

                        ConvertThisField = .true.

                    endif

                end do

                if(.not. ConvertThisField) cycle prop

            end if



            write(*,*)'Reading '//trim(PropertyName)//' fields'

            NewCurrentInstant = 0

            do CurrentInstant = StartInstant, EndInstant

                NewCurrentInstant = NewCurrentInstant + 1

                !Allocates new instance
                allocate (NewField)

                !Allocates aux instance
                allocate (AuxField)

                if (FirstTime) then 
                    call OutputInstants(NewCurrentInstant)
                endif

                !Get Father Fields 
                ObjFatherFile => Me%FirstFatherFile 

                do while (associated(ObjFatherFile))

                    !Allocates new instance
                    allocate (ObjFatherFile%NewField)

                    !Get field ID
                    call GetHDF5GroupID(ObjFatherFile%Info%ObjHDF5,                 & 
                                    "/Results/"//trim(PropertyName),                &
                                    CurrentInstant, ObjFatherFile%NewField%Name,    &
                                    ObjFatherFile%NewField%Units, Rank, Dimensions, &
                                    STAT = STAT_CALL)                                
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR160'

                    !check dimensions
                    if(Dimensions(1) .ne. ObjFatherFile%Info%WorkSize2D%IUB) then
                        write(*,*)'Fields size is not consistent with grid size : ' &
                        //trim(ObjFatherFile%Info%FileName)
                        stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR170'
                    end if

                    if(Dimensions(2) .ne. ObjFatherFile%Info%WorkSize2D%JUB) then
                        write(*,*)'Fields size is not consistent with grid size : ' & 
                        //trim(ObjFatherFile%Info%FileName)
                        stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR180'
                    end if 
            
                    !allocate field
                    nullify (ObjFatherFile%NewField%Values2D)
                    allocate(ObjFatherFile%NewField%Values2D                        &
                            (ObjFatherFile%Info%Size2D%ILB:                         &
                             ObjFatherFile%Info%Size2D%IUB,                         &
                             ObjFatherFile%Info%Size2D%JLB:                         &
                             ObjFatherFile%Info%Size2D%JUB))
            
                    call HDF5SetLimits (ObjFatherFile%Info%ObjHDF5,                 &
                                    ObjFatherFile%Info%WorkSize2D%ILB,              &
                                    ObjFatherFile%Info%WorkSize2D%IUB,              &
                                    ObjFatherFile%Info%WorkSize2D%JLB,              &
                                    ObjFatherFile%Info%WorkSize2D%JUB,              &
                                    STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    & 
                    stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR190'

                    !read field
                    call HDF5ReadData(ObjFatherFile%Info%ObjHDF5,                   &
                                  "/Results/"//trim(PropertyName),                  &
                                  trim(PropertyName),                               &
                                  Array2D      = ObjFatherFile%NewField%Values2D,   &
                                  OutputNumber = CurrentInstant, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR200'

                    ObjFatherFile%NewField%IDNumber = Count
                  
                    if(Me%TypeOfInterpolation == Triangulation) then

                        Count2 = 0

                        do j = ObjFatherFile%Info%WorkSize2D%JLB,                   &
                                    ObjFatherFile%Info%WorkSize2D%JUB
                        do i = ObjFatherFile%Info%WorkSize2D%ILB,                   &
                                    ObjFatherFile%Info%WorkSize2D%IUB

                            if (ObjFatherFile%Info%WaterPoints2D(i, j) ==           &
                                WaterPoint) then

                                Count2           = Count2 + 1

                                ObjFatherFile%NodeZ(Count2) =                       &
                                    ObjFatherFile%NewField%Values2D(i, j)
                

                            endif
    
                        enddo
                        enddo

                    endif

                    ObjFatherFile => ObjFatherFile%Next
                enddo
                
                !Construct new fields for father and son
                Count                   = Count + 1

                NewField%IDNumber       = Me%FirstFatherFile%NewField%IDNumber
                NewField%Name           = trim(PropertyName)

                AuxField%IDNumber       = Me%FirstFatherFile%NewField%IDNumber
                AuxField%Name           = trim(PropertyName)

                NewField%Units = trim(Me%FirstFatherFile%NewField%Units)
                NewField%Date  = Me%New%InstantsArray(NewCurrentInstant)

                AuxField%Units = trim(Me%FirstFatherFile%NewField%Units)
                AuxField%Date  = Me%New%InstantsArray(NewCurrentInstant)          

                select case (Rank)

                    case(2)

                        !allocate field
                        nullify (NewField%Values2D)
                        allocate(NewField%Values2D(AuxGrid%Size2D%ILB:              &
                                                   AuxGrid%Size2D%IUB,              &
                                                   AuxGrid%Size2D%JLB:              &
                                                   AuxGrid%Size2D%JUB))

                        NewField%Values2D(:,:) = FillValueReal

                        write(*,*)'Interpolating : '//trim(Me%FirstFatherFile%NewField%Name), &
                                                    Me%FirstFatherFile%NewField%IDNumber

                        !Cell cycle here
                        do j = Me%New%WorkSize2D%JLB, Me%New%WorkSize2D%JUB                                               
                        do i = Me%New%WorkSize2D%ILB, Me%New%WorkSize2D%IUB

                            !Get father ID
                            ObjFatherFile => Me%FirstFatherFile 

        do3 :               do while (associated(ObjFatherFile))
                                if (Me%RelevantFather(i,j) == ObjFatherFile%Info%Level) then
                                    exit do3
                                 else
                                    ObjFatherFile => ObjFatherFile%Next
                                endif
                            enddo do3

                            if(Me%TypeOfInterpolation == Triangulation)then

                                Me%New%WaterPoints2D => AuxWP2D
                               
                                call Triangulator (NewField, AuxGrid, ObjFatherFile, i, j)

                            else

                                write(*,*) 'Unknown type of interpolation'
                                stop       'StartInterpolateGrids - ModulePatchHDF5Files - ERR210' 

                            end if

                        enddo
                        enddo

                        call OutputFields(NewField, NewCurrentInstant)                  

                    case default 
                    
                        write(*,*)'Interpolation only available for 2D fields.'
                        stop 'OpenAndReadFatherHDF5Files - ModulePatchHDF5Files - ERR220'
                    
                end select

                !deAllocates new instance
                if(associated(NewField%Values2D)) deallocate (NewField%Values2D)
                if(associated(NewField%Values3D)) deallocate (NewField%Values3D)
                deallocate (NewField)
                nullify    (NewField)
           
                !deAllocates new instance
                if(associated(AuxField%Values2D)) deallocate (AuxField%Values2D)
                if(associated(AuxField%Values3D)) deallocate (AuxField%Values3D)
                deallocate (AuxField)
                nullify    (AuxField)

                ObjFatherFile => Me%FirstFatherFile 

                do while (associated(ObjFatherFile))

                    !deAllocates new instance
                    if(associated(ObjFatherFile%NewField%Values2D))                 & 
                        deallocate (ObjFatherFile%NewField%Values2D)
                    if(associated(ObjFatherFile%NewField%Values3D))                 & 
                        deallocate (ObjFatherFile%NewField%Values3D)
                    deallocate (ObjFatherFile%NewField)
                    nullify    (ObjFatherFile%NewField)

                    ObjFatherFile => ObjFatherFile%Next
                enddo

            end do

            FirstTime = .false.

        end do prop

        if (Me%TypeOfInterpolation == Triangulation) then

            deallocate(Me%SonCenterX)
            deallocate(Me%SonCenterY)

            call UnGetHorizontalMap(AuxGrid%ObjHorizontalMap,                       & 
                                    AuxWP2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'Open_HDF5_OutPut_File - ModulePatchHDF5Files - ERR230'
        endif

        nullify(AuxGrid)


    end subroutine OpenAndReadFatherHDF5Files

    !--------------------------------------------------------------------------

    subroutine Triangulator (NewField, NewGrid, ObjFatherFile, i, j)

        !Arguments-------------------------------------------------------------
        type(T_Field),           pointer            :: NewField
        type(T_Grid )                               :: NewGrid
        type (T_Father),     pointer                :: ObjFatherFile
        integer                                     :: i, j

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        logical                                     :: FillOutsidePoints   = .false.
        
        !Begin-----------------------------------------------------------------


iN1:     if (ObjFatherFile%NumberOfNodes >= 3) then

        call SetHeightValues(ObjFatherFile%ObjTriangulation, ObjFatherFile%NodeZ,       & 
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Triangulator - ModulePatchHDF5Files - ERR10'

        if(NewGrid%WaterPoints2D(i, j) == WaterPoint) then

                NewField%Values2D(i, j) = InterPolation(ObjFatherFile%ObjTriangulation, &
                                                        Me%SonCenterX(i,j),             &
                                                        Me%SonCenterY(i,j),             &
                                                        FillOutsidePoints,              &
                                                        Default = null_real,            &
                                                        STAT    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              & 
                    stop 'Triangulator - ModuleInterpolateGrids - ERR20'

        end if

        endif iN1

    end subroutine Triangulator   

    !--------------------------------------------------------------------------

    subroutine GetFatherTriangParam (ObjFatherFile)

        !Arguments-------------------------------------------------------------
        type (T_Father),     pointer                :: ObjFatherFile

        !Local-----------------------------------------------------------------
        integer                                     :: Count, iaux, jaux
        
        !Begin-----------------------------------------------------------------

        ObjFatherFile%NumberOfNodes =  Sum(ObjFatherFile%Info%WaterPoints2D             &
                                           (ObjFatherFile%Info%WorkSize2D%ILB:          &
                                            ObjFatherFile%Info%WorkSize2D%IUB,          &
                                            ObjFatherFile%Info%WorkSize2D%JLB:          &
                                            ObjFatherFile%Info%WorkSize2D%JUB))

iN2:     if (ObjFatherFile%NumberOfNodes >= 3) then

        allocate(ObjFatherFile%NodeX(ObjFatherFile%NumberOfNodes))
        allocate(ObjFatherFile%NodeY(ObjFatherFile%NumberOfNodes))
        allocate(ObjFatherFile%NodeZ(ObjFatherFile%NumberOfNodes))

        Count = 0

        do jaux = ObjFatherFile%Info%WorkSize2D%JLB, ObjFatherFile%Info%WorkSize2D%JUB
        do iaux = ObjFatherFile%Info%WorkSize2D%ILB, ObjFatherFile%Info%WorkSize2D%IUB

            if (ObjFatherFile%Info%WaterPoints2D(iaux, jaux) == WaterPoint) then

                Count           = Count + 1

                ObjFatherFile%NodeX(Count) =                                               &
                                   ((ObjFatherFile%Info%ConnectionX(iaux, jaux  ) +        &
                                    ObjFatherFile%Info%ConnectionX(iaux+1, jaux  ))/2. +   &
                                   (ObjFatherFile%Info%ConnectionX(iaux, jaux+1) +         &
                                    ObjFatherFile%Info%ConnectionX(iaux+1, jaux+1))/2.)/2.
        
                ObjFatherFile%NodeY(Count) =                                               &
                                   ((ObjFatherFile%Info%ConnectionY(iaux, jaux  ) +        &
                                    ObjFatherFile%Info%ConnectionY(iaux+1, jaux  ))/2. +   &
                                   (ObjFatherFile%Info%ConnectionY(iaux, jaux+1) +         &
                                    ObjFatherFile%Info%ConnectionY(iaux+1, jaux+1))/2.)/2.              

            endif

        enddo
        enddo

        !else?

        endif iN2

    end subroutine GetFatherTriangParam   

    !--------------------------------------------------------------------------

    type(T_Time) function HDF5TimeInstant(Instant,ObjFatherFile)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        type (T_Father),     pointer            :: ObjFatherFile
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer             :: TimeVector

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (ObjFatherFile%Info%ObjHDF5, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = ObjFatherFile%Info%ObjHDF5,               &
                             GroupName      = "/Time",                                  &
                             Name           = "Time",                                   &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModulePatchHDF5Files - ERR10'

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
        
        CurrentDate = Me%New%InstantsArray(CurrentInstant)

        call ExtractDate   (CurrentDate,                                            &
                            AuxTime(1), AuxTime(2), AuxTime(3),                     &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime

        call HDF5SetLimits  (Me%New%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputInstants - ModulePatchHDF5Files - ERR10'


        call HDF5WriteData  (Me%New%ObjHDF5, "/Time",                               &
                             "Time", "YYYY/MM/DD HH:MM:SS",                         &
                             Array1D = TimePtr,                                     &
                             OutputNumber = CurrentInstant, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputInstants - ModulePatchHDF5Files - ERR20'


    end subroutine OutputInstants

    !------------------------------------------------------------------------

    subroutine OutputFields(NewField, OutputNumber)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                          :: NewField
        integer                                         :: OutputNumber

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (Me%New%ObjHDF5, Me%New%WorkSize2D%ILB,         & 
                             Me%New%WorkSize2D%IUB,                         &
                             Me%New%WorkSize2D%JLB, Me%New%WorkSize2D%JUB,  &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModulePatchHDF5Files - ERR10'

        call HDF5WriteData(Me%New%ObjHDF5,                                  &
                           "/Results/"//NewField%Name,                      &
                           NewField%Name,                                   &
                           NewField%Units,                                  &
                           Array2D      = NewField%Values2D,                &
                           OutputNumber = OutputNumber,                     &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModulePatchHDF5Files - ERR20'

        !Writes everything to disk
        call HDF5FlushMemory (Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModulePatchHDF5Files - ERR30'

    end subroutine OutputFields
    
    !--------------------------------------------------------------------------

    subroutine KillPatchHDF5Files

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                         :: nUsers
        type (T_Father),     pointer                    :: ObjFatherFile
        
        !Begin-----------------------------------------------------------------

        ObjFatherFile => Me%FirstFatherFile 

        do while (associated(ObjFatherFile))
            if (Me%TypeOfInterpolation == Triangulation) then

                call KillTriangulation (ObjFatherFile%ObjTriangulation, STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'KillPatchHDF5Files - ModulePatchHDF5Files - ERR10'
                
                !Deallocate variables needed triangulation
                deallocate(ObjFatherFile%NodeX)
                deallocate(ObjFatherFile%NodeY)
                deallocate(ObjFatherFile%NodeZ)

            endif

            call KillFatherGrid(ObjFatherFile)

            ObjFatherFile => ObjFatherFile%Next
        enddo

        call UnGetHorizontalGrid(Me%New%ObjHorizontalGrid, Me%New%ConnectionX,  & 
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillPatchHDF5Files - ModulePatchHDF5Files - ERR20'


        call UnGetHorizontalGrid(Me%New%ObjHorizontalGrid, Me%New%ConnectionY,  & 
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModulePatchHDF5Files - ERR30'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillInterpolateGrids - ModulePatchHDF5Files - ERR40' 

        call KillHorizontalMap(Me%New%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModulePatchHDF5Files - ERR50'

        call KillGridData(Me%New%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModulePatchHDF5Files - ERR60'

        call KillHorizontalGrid(Me%New%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModulePatchHDF5Files - ERR70'
        
        call KillHDF5(Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModulePatchHDF5Files - ERR80'     
   
        deallocate(Me)
        nullify(Me)

    end subroutine KillPatchHDF5Files     

    !--------------------------------------------------------------------------
 
    subroutine KillFatherGrid(ObjFatherFile)
        !Arguments-------------------------------------------------------------
        type (T_Father),     pointer                :: ObjFatherFile
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
       
        write(*,*)'Killing father grid...'


        call UnGetHorizontalGrid(ObjFatherFile%Info%ObjHorizontalGrid,                  & 
                                 ObjFatherFile%Info%ConnectionX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillFatherGrid - ModulePatchHDF5Files - ERR10'


        call UnGetHorizontalGrid(ObjFatherFile%Info%ObjHorizontalGrid,                  &
                                 ObjFatherFile%Info%ConnectionY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillFatherGrid - ModulePatchHDF5Files - ERR20'

        call KillHorizontalMap( ObjFatherFile%Info%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModulePatchHDF5Files - ERR30'

        call KillGridData(ObjFatherFile%Info%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModulePatchHDF5Files - ERR40'

        call KillHDF5(ObjFatherFile%Info%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModulePatchHDF5Files - ERR50'

        call KillHorizontalGrid(ObjFatherFile%Info%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModulePatchHDF5Files - ERR60'


    end subroutine KillFatherGrid


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

end module ModulePatchHDF5Files









