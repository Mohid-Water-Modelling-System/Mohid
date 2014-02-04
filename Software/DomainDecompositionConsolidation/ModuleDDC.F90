!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model - Results consolidation
! PROJECT       : DDC - Domain Decomposition Consolidation
! PROGRAM       : MainDDC
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          :  2013
! REVISION      : Ricardo Miranda/Joao Ribeiro
! DESCRIPTION   : Program to consolidate results from MPI run with domain decomposition
!
!
!------------------------------------------------------------------------------


Module ModuleDDC

    use ModuleGlobalData
    use ModuleHashTable
    use ModuleHDF5
    use HDF5

    implicit none

    private

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructDDC
    private ::      ReadTreeFile
    private ::          ReadLines
    private ::              ModelLevel
    private ::              ModelPath
    private ::              AllocateDirectoryList
    
    !Selector

    !Modifier
    public  :: ModifyDDC
    private ::      ScanDirectoryList
    private ::          ScanFileList
    private ::              OpenDecomposedFiles
    private ::                  OpenDecomposedFiles2
    private ::                      ScanDecomposedFiles
    private ::                          GetHDF5ReadWriteFileAccessCode
    private ::                          GetDomainDecomposition
    private ::                              GetHDF5ReadFileAccessCode
    private ::                                 GetMappingValues
    private ::          WriteConsolidatedHDF
    private ::              SetArrayLimits
    private ::              ReadHDFDataSetHyperslab
    private ::              ReadHDFDataSet
    private ::              DefineHDF5Limits
    private ::              CreateGlobalSize
    private ::              WriteHDFDataSet
    private ::              WriteTimeValues
    
    !Destructor
    public  :: KillDDC
    private ::      DeAllocateInstance
    private ::          deallocateDirectoryList
    private ::               KillConsolidatedFiles
    private ::                    KillConsolidatedFiles2

    !Management
    private ::      Ready
    private ::      Read_Lock_DDC
    private ::      Read_Unlock_DDC
    private ::      VerifyReadLock_DDC

    !Interfaces----------------------------------------------------------------

    !Types---------------------------------------------------------------------

    public  :: T_DDC
    type       T_DDC
        private

        !nbrModels stores total number of models
        integer                                     :: nbrModels        = NULL_INT
        type(T_DirectoryList), pointer              :: DirectoryList
    end type  T_DDC

    private :: T_DirectoryList
    type       T_DirectoryList
        character(PathLength)                       :: Directory        = NULL_STR
        type(T_HashTable), pointer                  :: hash_map_out !consolidated HDF5 output files
        type(T_HashTable), pointer                  :: hash_map_in  !decomposed HDF5 input files
        type(T_DirectoryList), pointer              :: Next
    end type  T_DirectoryList

    !--------------------------------------------------------------------------

    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    function ConstructDDC()

        !Function----------------------------------------------------------------
        type (T_DDC), pointer                       :: ConstructDDC

        !External----------------------------------------------------------------
        type (T_DDC), pointer                       :: Me

        !Local-------------------------------------------------------------------
        type (T_DDC), pointer                       :: NewObjDDC
        integer                                     :: STAT_CALL    = NULL_INT

        !------------------------------------------------------------------------

        allocate(NewObjDDC)
        NewObjDDC%nbrModels = 0         !Initializes model count

        Me => NewObjDDC

        STAT_CALL = ReadTreeFile(Me)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDDC - ModuleDDC - ERR01'

        ConstructDDC => Me

        !----------------------------------------------------------------------

    end function ConstructDDC

    !--------------------------------------------------------------------------

    integer function ReadTreeFile(Me)

        !Arguments-------------------------------------------------------------
        type (T_DDC), pointer                       :: Me

        !Local-----------------------------------------------------------------
        logical                                     :: TreeExists
        integer                                     :: STAT_        
        integer                                     :: STAT_CALL    
        integer                                     :: iTree        
        character(StringLength)                     :: Coment1      
        character(StringLength)                     :: Coment2      

        !------------------------------------------------------------------------

        STAT_        = NULL_INT
        STAT_CALL    = NULL_INT
        iTree        = NULL_INT
        Coment1      = NULL_STR
        Coment2      = NULL_STR

        !Verifies if Tree file exists and allocates the list of models
        inquire(file='Tree.dat', EXIST = TreeExists)
        if (.NOT. TreeExists) STAT_ = FILE_NOT_FOUND_ERR_

if1:    if (TreeExists) then
            call UnitsManager(iTree, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadTreeFile - ModuleDDC - ERR01'

            open(UNIT = iTree, FILE = 'Tree.dat', status = 'OLD', IOSTAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadTreeFile - ModuleDDC - ERR02'

            read(unit=iTree, fmt=*) Coment1
            read(unit=iTree, fmt=*) Coment2

            call ReadLines(Me,                                                  &
                           iTree         = iTree,                               &
                           DirectoryList = Me%DirectoryList)

            call UnitsManager(iTree, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadTreeFile - ModuleDDC - ERR03'

            STAT_ = SUCCESS_
        endif if1

        ReadTreeFile = STAT_

        !------------------------------------------------------------------------

    end function ReadTreeFile

    !--------------------------------------------------------------------------

    recursive subroutine ReadLines(Me, iTree, DirectoryList)

        !Arguments-------------------------------------------------------------
        type (T_DDC), pointer                                   :: Me
        type(T_DirectoryList), pointer                          :: DirectoryList
        integer, intent(IN)                                     :: iTree

        !Local-----------------------------------------------------------------
        integer                                                 :: iMPI         
        integer                                                 :: STAT_CALL    
        integer                                                 :: ModelLevel_  
        character(PathLength)                                   :: AuxString    
        character(PathLength)                                   :: AuxString2   
        character(PathLength)                                   :: ModelPath_   
        integer                                                 :: nbrModels_   
        integer                                                 :: int

        !------------------------------------------------------------------------
        
        iMPI         = NULL_INT
        STAT_CALL    = NULL_INT
        ModelLevel_  = NULL_INT
        AuxString    = null_str
        AuxString2   = null_str
        ModelPath_   = null_str
        nbrModels_   = NULL_INT
        
        read(unit = iTree, fmt='(a256)', IOSTAT = STAT_CALL) AuxString
if1 :   if (STAT_CALL .EQ. SUCCESS_) then
            AuxString2 = trim(adjustl(AuxString))
if2 :       if (AuxString2(1:1) == '+') then
                iMPI = scan(AuxString2,":")
if3 :           if (iMPI > 0) then
                     read(AuxString2(iMPI+1:),'(I)') nbrModels_
!                    nbrModels_ = int(AuxString2(iMPI+1:))
!                    read(AuxString2(iMPI+1:), '(i10)') nbrModels_
                    Me%nbrModels = Me%nbrModels + nbrModels_

                    ModelLevel_  = ModelLevel(AuxString  = AuxString2,          &
                                              Level      = 0)   !Level=0 because it is the '+' counting start
                    ModelPath_   = ModelPath (AuxString  = AuxString2,          &
                                              Level      = ModelLevel_)

                    DirectoryList => AllocateDirectoryList(ModelPath = ModelPath_)
                    call ReadLines(Me,                                          &
                                   iTree         = iTree,                       &
                                   DirectoryList = DirectoryList%Next)
                else if3
                    !This line does not have domain decomposition, reads next Tree.dat's line
                    Me%nbrModels = Me%nbrModels + 1

                    call ReadLines(Me,                                          &
                                   iTree         = iTree,                       &
                                   DirectoryList = DirectoryList)
                endif if3
            endif if2
        endif if1

        !------------------------------------------------------------------------

    end subroutine ReadLines


    !--------------------------------------------------------------------------

    integer pure recursive function ModelLevel(AuxString, Level)

        !Arguments-------------------------------------------------------------
        Character(len=*), intent(in)                                :: AuxString

        !Local-----------------------------------------------------------------
        integer, intent(in)                                         :: Level

        !------------------------------------------------------------------------

 if1 :  if (AuxString((Level+1):(Level+1)) == '+') then
            ModelLevel = ModelLevel(AuxString = AuxString,                      &
                                    Level     = Level + 1)
        else if1
            ModelLevel = Level
        endif if1

        !------------------------------------------------------------------------

    end function ModelLevel

    !--------------------------------------------------------------------------

    character(len=PathLength) pure function ModelPath (AuxString, Level)

        !Arguments-------------------------------------------------------------
        character(len=*), intent(in)                                :: AuxString
        integer,          intent(in)                                :: Level

        !Local-----------------------------------------------------------------
        integer                                                     :: position

        !------------------------------------------------------------------------

        position   = scan(AuxString, ":")
        ModelPath  = AuxString(Level+1:position-5)//"res"

        !------------------------------------------------------------------------

    end function ModelPath


    !--------------------------------------------------------------------------

    function AllocateDirectoryList(ModelPath)

        !Function--------------------------------------------------------------
        type(T_DirectoryList), pointer                  :: AllocateDirectoryList

        !Arguments-------------------------------------------------------------
        character(PathLength), intent(IN)               :: ModelPath

        !Local-------------------------------------------------------------------
        type (T_DirectoryList), pointer                 :: NewDirectoryList

        !------------------------------------------------------------------------

        allocate(NewDirectoryList)
        NewDirectoryList%Directory = adjustl(trim(ModelPath))
        NewDirectoryList%hash_map_out => hash_init()
        NewDirectoryList%hash_map_in => hash_init()
        nullify(NewDirectoryList%Next)

        AllocateDirectoryList => NewDirectoryList

        !------------------------------------------------------------------------

    end function AllocateDirectoryList

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

    integer function ModifyDDC(Me)

        !Arguments-------------------------------------------------------------
        type (T_DDC), pointer                                       :: Me

        !Local-----------------------------------------------------------------
        integer                                                     :: STAT_        
        integer                                                     :: ready_       

        !----------------------------------------------------------------------

        STAT_        = NULL_INT
        ready_       = NULL_INT

        STAT_ = UNKNOWN_

        ready_ = Ready(Me)

!        if (ready_ .EQ. IDLE_ERR_) then
if1 :       if (associated(Me%DirectoryList)) then
                call ScanDirectoryList(DirectoryList = Me%DirectoryList,    &
                                       nbrModels     = Me%nbrModels)
            endif if1

            STAT_ = SUCCESS_
!        else
!            STAT_ = ready_
!        end if

        ModifyDDC = STAT_

        !------------------------------------------------------------------------

    end function ModifyDDC

    !--------------------------------------------------------------------------

    recursive subroutine ScanDirectoryList(DirectoryList, nbrModels)

        !Arguments-------------------------------------------------------------
        type(T_DirectoryList), pointer                          :: DirectoryList
        integer, intent(IN)                                     :: nbrModels

        !Local--------------------------------------------------------------

        !External--------------------------------------------------------------
        CHARACTER(PathLength)                   :: FirstHDFFile
        INTEGER                                 :: TID
        INTEGER                                 :: OMP_GET_THREAD_NUM  , i, j, k

        integer omp_get_nested
        !------------------------------------------------------------------------

!$OMP PARALLEL
!$OMP SECTIONS
!$OMP SECTION
!        TID = OMP_GET_THREAD_NUM()
!        k= omp_get_nested()
!        PRINT *, 'ScanDirectoryList 1, OpenMP thread = ', TID
!        PRINT *, 'omp_get_nested = ', k
        call OpenDecomposedFiles(hash_map_out = DirectoryList%hash_map_out,        &
                                 hash_map_in  = DirectoryList%hash_map_in,         &
                                 Directory = DirectoryList%Directory,              &
                                 nbrModels = nbrModels)

if1 :   if (hash_get_first_exists(DirectoryList%hash_map_in)) then

            FirstHDFFile = hash_get_first_key(DirectoryList%hash_map_in)
            call ScanFileList (HDFFile     = FirstHDFFile,                          &
                               hash_map_in = DirectoryList%hash_map_in)
        endif if1

        if (associated(DirectoryList%Next)) call ScanDirectoryList(DirectoryList = DirectoryList%Next,  &
                                                                   nbrModels     = nbrModels)
!$OMP END SECTIONS NOWAIT
!$OMP END PARALLEL


        !------------------------------------------------------------------------

    end subroutine ScanDirectoryList

    !--------------------------------------------------------------------------

    recursive subroutine ScanFileList(HDFFile,                                  &
                                      hash_map_in)

        !Arguments-------------------------------------------------------------
        CHARACTER(PathLength), intent(IN)                       :: HDFFile
        type(T_HashTable), pointer                              :: hash_map_in 
         
        !Local--------------------------------------------------------------
        character(StringLength)                                 :: GroupName    
        integer                                                 :: IDIn              
        integer                                                 :: ObjHDF5_In
        integer                                                 :: IDOut           
        integer                                                 :: ObjHDF5_Out  
        integer                                                 :: STAT_CALL 
        CHARACTER(PathLength)                                   :: HDFFileNext 
        integer, dimension(4)                                   :: DomainSize
        integer, dimension(4)                                   :: WindowPosition
        integer, dimension(4)                                   :: WindowFrame

        !------------------------------------------------------------------------

        GroupName    = ""
        IDIn         = NULL_INT
        ObjHDF5_In   = NULL_INT 
        IDOut        = NULL_INT      
        ObjHDF5_Out  = NULL_INT
        STAT_CALL    = NULL_INT
        HDFFileNext  = NULL_STR       
    
        IDOut       = hash_get     (hash_map_in, HDFFile)
        ObjHDF5_Out = hash_getObjID(hash_map_in, HDFFile)
                    
        call GetHDF5ReadFileAccessCode(HDFFile    = HDFFile,                        &
                                       IDIn       = IDIn,                           &
                                       ObjHDF5_In = ObjHDF5_In)
        
        DomainSize     = hash_getDomainSize(hash_map_in, HDFFile)
        WindowPosition = hash_getWindowPosition(hash_map_in, HDFFile)
        WindowFrame = hash_getWindowFrame(hash_map_in, HDFFile)
        
        call WriteConsolidatedHDF ( IDOut           = IDOut,                         &
                                    ObjHDF5_Out     = ObjHDF5_Out,                   &
                                    IDIn            = IDIn,                          &
                                    ObjHDF5_In      = ObjHDF5_In,                    &
                                    GroupName       =  adjustl(trim(GroupName)),     &
                                    DomainSize      = DomainSize,                    &
                                    WindowPosition  = WindowPosition,                &
                                    WindowFrame     = WindowFrame,                   &
                                    HDFFile         = HDFFile)
        
        call KillHDF5(  HDF5ID   = ObjHDF5_In,                                      &
                        STAT     = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ScanFileList - ModuleDDC - ERR01'
        
if1 :   if (hash_get_next_exists(hash_map_in, HDFFile)) then

            HDFFileNext = hash_get_next_key(hash_map_in, HDFFile)
            call ScanFileList(HDFFileNext, hash_map_in)
            
        endif if1
        
        !------------------------------------------------------------------------

    end subroutine ScanFileList

    !--------------------------------------------------------------------------

    recursive subroutine WriteConsolidatedHDF(IDOut, ObjHDF5_Out, IDIn,     &
                                              ObjHDF5_In, GroupName,        &
                                              DomainSize, WindowPosition,   &
                                              WindowFrame, HDFFile)

        !Arguments-------------------------------------------------------------
        character(*), intent(IN)                                :: GroupName  
        integer, intent(IN)                                     :: IDIn               
        integer, intent(IN)                                     :: ObjHDF5_In    
        integer, intent(IN)                                     :: IDOut              
        integer, intent(IN)                                     :: ObjHDF5_Out
        integer, dimension(4), intent(IN)                       :: DomainSize
        integer, dimension(4), intent(IN)                       :: WindowPosition
        integer, dimension(4), intent(IN)                       :: WindowFrame
        character(*), intent(IN)                                :: HDFFile
        
        !Local-------------------------------------------------------------------
        logical                                                 :: Exist
        integer, dimension(:),       pointer                    :: IArray1D
        integer, dimension(:, :),    pointer                    :: IArray2D
        integer, dimension(:, :, :), pointer                    :: IArray3D
        real(4), dimension(:),       pointer                    :: R4Array1D
        real(4), dimension(:, :),    pointer                    :: R4Array2D
        real(4), dimension(:, :, :), pointer                    :: R4Array3D
        real(8), dimension(:),       pointer                    :: R8Array1D
        real(8), dimension(:, :),    pointer                    :: R8Array2D
        real(8), dimension(:, :, :), pointer                    :: R8Array3D
        character(StringLength)                                 :: Name       
        character(StringLength)                                 :: Units
        character(StringLength)                                 :: NewGroupName 
        integer(HID_T)                                          :: GroupType
        integer(HID_T)                                          :: DataType
        integer(HSIZE_T), dimension(7)                          :: dims  
        integer                                                 :: nItems       
        integer                                                 :: Rank
        integer                                                 :: idx 
        integer                                                 :: STAT_CALL         
        integer, dimension(7)                                   :: Dimensions
        integer, dimension(6)                                   :: LimitsArray
        integer                                                 :: LimitsArrayFactor
        real                                                    :: Minimum
        real                                                    :: Maximum
        
        !------------------------------------------------------------------------
            
        DataType        = NULL_INT
        NewGroupName    = NULL_STR
        idx             = NULL_INT
        LimitsArray     = 1

        !Get the number of members in the Group  
        call GetHDF5GroupNumberOfItems (HDF5ID = ObjHDF5_In,                         &
                                        GroupName = adjustl(trim(GroupName))//"/",   &
                                        nItems = nItems, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteConsolidatedHDF - ModuleDDC - ERR01'

do1 :   do idx = 1, nItems

            LimitsArrayFactor = 0
        
            !Gets information about the group
            call GetHDF5ObjectInfo (HDF5ID = ObjHDF5_In, FatherGroupName = adjustl(trim(GroupName))//"/",           &
                                    GroupPosition = idx, GroupName = Name,GroupType = GroupType,                    &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteConsolidatedHDF - ModuleDDC - ERR02'
                          
if1 :       if (GroupType == H5G_DATASET_F) then

                !Check if Dataset Exists
                call GetHDF5DataSetExist(   ObjHDF5_Out, adjustl(trim(GroupName))//"/"//adjustl(trim(Name)),        &
                                            Exist, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteConsolidatedHDF - ModuleDDC - ERR03'
                
                !Get Dataset Exists
                call GetHDF5GroupID(    HDF5ID = ObjHDF5_In, FatherGroupName = adjustl(trim(GroupName))//"/",       &
                                        GroupPosition = idx, GroupName = Name,                                      &
                                        Units = Units, Rank = Rank, Dimensions = Dimensions,                        &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteConsolidatedHDF - ModuleDDC - ERR04'                       
                
                !Get Dataset DataType
                call GetHDF5DataTypeID( HDF5ID = ObjHDF5_In, FatherGroupName = adjustl(trim(GroupName))//"/",       &
                                        GroupPosition = idx, GroupName = Name,                                      &
                                        DataType = DataType, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteConsolidatedHDF - ModuleDDC - ERR05'                                               
                
                !Get Dataset Value Real Minimum Attribute
                call HDF5ReadGenericRealAttribute(  ObjHDF5_In, adjustl(trim(GroupName)), adjustl(trim(Name)), 3,   &
                                                    "Minimum", Minimum, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteConsolidatedHDF - ModuleDDC - ERR06'
                
                !Get Dataset Value Real Maximum Attribute
                call HDF5ReadGenericRealAttribute(  ObjHDF5_In, adjustl(trim(GroupName)), adjustl(trim(Name)), 3,   &
                                                    "Maximum", Maximum, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteConsolidatedHDF - ModuleDDC - ERR07'
         
                !Check if it is Coordenates Feature 
if12 :          if (adjustl(trim(Name)) .EQ. "ConnectionX" .OR. adjustl(trim(Name)) .EQ. "ConnectionY" .OR.         &
                    adjustl(trim(Name)) .EQ. "Latitude" .OR. adjustl(trim(Name)) .EQ. "Longitude") then
                    
                    LimitsArrayFactor = 1
                    
                endif if12
                
                !Set Array Limits
                call SetArrayLimits(    Rank = Rank, GroupName = adjustl(trim(GroupName)),                          &
                                        WindowPosition = WindowPosition, Dimensions = Dimensions,                   &
                                        LimitsArrayFactor = LimitsArrayFactor, LimitsArray = LimitsArray)
                
                !Read Array 
if17 :          if(adjustl(trim(GroupName)) .NE. "/Time") then
                
                    call ReadHDFDataSetHyperslab(   ObjHDF5_In = ObjHDF5_In,                       &
                                                    GroupName   = adjustl(trim(GroupName)),                                 &
                                                    obj_name    = adjustl(trim(Name)),                                      &
                                                    Rank        = Rank,                                                     &
                                                    InnerWindow = LimitsArray,                                              &
                                                    WindowFrame = WindowFrame,                                              &
                                                    DataType    = DataType,                                                 &
                                                    IArray1D    = IArray1D,  IArray2D  = IArray2D,  IArray3D  = IArray3D,   &
                                                    R4Array1D   = R4Array1D, R4Array2D = R4Array2D, R4Array3D = R4Array3D,  &
                                                    R8Array1D   = R8Array1D, R8Array2D = R8Array2D, R8Array3D = R8Array3D)
                
                else if17
                
                    call ReadHDFDataSet(ObjHDF5_In = ObjHDF5_In,                  &
                                        GroupName  = adjustl(trim(GroupName)),    &
                                        obj_name   = adjustl(trim(Name)),         &
                                        Rank       = Rank,                        &
                                        dims       = Dimensions,                  &
                                        DataType   = DataType,                    &
                                        IArray1D   = IArray1D,  IArray2D  = IArray2D,  IArray3D  = IArray3D,            &
                                        R4Array1D  = R4Array1D, R4Array2D = R4Array2D, R4Array3D = R4Array3D,           &
                                        R8Array1D  = R8Array1D, R8Array2D = R8Array2D, R8Array3D = R8Array3D)      
                endif if17
                
if5 :           if(Exist .AND. adjustl(trim(GroupName)) .NE. "/Time") then

                    !HDF5 SetLimits 
                    call DefineHDF5Limits  (ObjHDF5_Out = ObjHDF5_Out,                                                  &
                                            GroupName   = adjustl(trim(GroupName)),                                     &
                                            obj_name    = adjustl(trim(Name)),                                          &
                                            Rank        = Rank,                                                         &
                                            LimitsArray = LimitsArray,                                                  &
                                            WindowFrame = WindowFrame,                                                  &
                                            DomainSize  = DomainSize,                                                   &
                                            InnerWindow = WindowPosition,                                               &
                                            LimitsArrayFactor = LimitsArrayFactor)
                    !HDF5 Write
                    call WriteHDFDataSet   (ObjHDF5_Out = ObjHDF5_Out,                                                  &
                                            GroupName  = adjustl(trim(GroupName)),                                      &
                                            obj_name   = adjustl(trim(Name)),                                           &
                                            Rank       = Rank,                                                          &
                                            DataType   = DataType,                                                      &
                                            Minimum    = Minimum,                                                       &
                                            Maximum    = Maximum,                                                       &
                                            IArray1D   = IArray1D,  IArray2D  = IArray2D,  IArray3D  = IArray3D,        &
                                            R4Array1D  = R4Array1D, R4Array2D = R4Array2D, R4Array3D = R4Array3D,       &
                                            R8Array1D  = R8Array1D, R8Array2D = R8Array2D, R8Array3D = R8Array3D)

                elseif (adjustl(trim(GroupName)) .NE. "/Time") then if5
                    
                    !Create Global HDF5
                    call CreateGlobalSize  (ObjHDF5_Out  = ObjHDF5_Out,                                                 &
                                            GroupName  = adjustl(trim(GroupName)),                                      &
                                            obj_name   = adjustl(trim(Name)),                                           &
                                            Units      = Units,                                                         &
                                            Rank       = Rank,                                                          &
                                            DomainSize = DomainSize,                                                    &
                                            dims       = Dimensions,                                                    &
                                            DataType   = DataType)
                    !HDF5 SetLimits 
                    call DefineHDF5Limits  (ObjHDF5_Out = ObjHDF5_Out,                                                  &
                                            GroupName   = adjustl(trim(GroupName)),                                     &
                                            obj_name    = adjustl(trim(Name)),                                          &
                                            Rank        = Rank,                                                         &
                                            LimitsArray = LimitsArray,                                                  &
                                            WindowFrame = WindowFrame,                                                  &
                                            DomainSize  = DomainSize,                                                   &
                                            InnerWindow = WindowPosition,                                               &
                                            LimitsArrayFactor = LimitsArrayFactor)
                    !HDF5 Write                        
                    call WriteHDFDataSet   (ObjHDF5_Out = ObjHDF5_Out,                                                  &
                                            GroupName  = adjustl(trim(GroupName)),                                      &
                                            obj_name   = adjustl(trim(Name)),                                           &
                                            Rank       = Rank,                                                          &
                                            DataType   = DataType,                                                      &
                                            Minimum    = Minimum,                                                          &
                                            Maximum    = Maximum,                                                          &
                                            IArray1D   = IArray1D,  IArray2D  = IArray2D,  IArray3D  = IArray3D,        &
                                            R4Array1D  = R4Array1D, R4Array2D = R4Array2D, R4Array3D = R4Array3D,       &
                                            R8Array1D  = R8Array1D, R8Array2D = R8Array2D, R8Array3D = R8Array3D)
                    
                elseif (adjustl(trim(GroupName)) .EQ. "/Time" .AND. Exist .EQ. .FALSE.) then if5
                
                    !HDF5 Write Time
                    call WriteTimeValues(   HDF5ID = ObjHDF5_Out, Dimensions = Dimensions, DataType = DataType,     &
                                            GroupName = GroupName, Name = adjustl(trim(Name)), Minimum = Minimum,   &
                                            Maximum = Maximum, Units  = Units, IArray1D = IArray1D,                 &
                                            R4Array1D = R4Array1D, R8Array1D = R8Array1D)

                end if if5
               
            elseif (GroupType ==H5G_GROUP_F) then if1
             
if3 :           if(adjustl(trim(GroupName)) .EQ. "/") then
                    
                    NewGroupName = adjustl(trim(GroupName))//adjustl(trim(Name))
                    
                elseif (adjustl(trim(GroupName))//"/"//adjustl(trim(Name)) .NE.                                         &
                        "/Grid/Decomposition") then if3
                    
                    NewGroupName = adjustl(trim(GroupName))//"/"//adjustl(trim(Name))
                    
                endif if3
                
if2 :           if (adjustl(trim(GroupName))//"/"//adjustl(trim(Name)) .NE. &
                    "/Grid/Decomposition") then
                    
                    call WriteConsolidatedHDF  (IDOut           = IDOut,                                                    &
                                                ObjHDF5_Out     = ObjHDF5_Out,                                              &
                                                IDIn            = IDIn,                                                     &
                                                ObjHDF5_In      = ObjHDF5_In,                                               &
                                                GroupName       = NewGroupName,                                             &
                                                DomainSize      = DomainSize,                                               &
                                                WindowPosition  = WindowPosition,                                           &
                                                WindowFrame     = WindowFrame,                                              &
                                                HDFFile         = HDFFile)
                                  
                endif if2

            endif if1
             
        end do do1

        !------------------------------------------------------------------------

    end subroutine WriteConsolidatedHDF

        !--------------------------------------------------------------------------

    subroutine WriteTimeValues( HDF5ID, Dimensions, DataType, GroupName, Name,  & 
                                Minimum, Maximum, Units, IArray1D, R4Array1D,   &
                                R8Array1D)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                     :: HDF5ID
        integer, dimension(7), intent(IN)                       :: Dimensions
        integer(HID_T), intent(IN)                              :: DataType
        character(*), intent(IN)                                :: GroupName
        character(*), intent(IN)                                :: Name
        character(*), intent(IN)                                :: Units
        real, intent(IN)                                        :: Minimum
        real, intent(IN)                                        :: Maximum
        integer, dimension(:), pointer, intent(IN), optional    :: IArray1D
        real(4), dimension(:), pointer, intent(IN), optional    :: R4Array1D
        real(8), dimension(:), pointer, intent(IN), optional    :: R8Array1D
        
        !Local-------------------------------------------------------------------
        integer                                                 :: STAT_CALL
        
        !------------------------------------------------------------------------
                    
        call HDF5SetLimits  (HDF5ID = HDF5ID, ILB = 1, IUB = Dimensions(1), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteTimeValues - ModuleDDC - ERR01'
                    
if1 :   if(DataType .EQ. H5T_NATIVE_INTEGER)   then

            call HDF5WriteData(HDF5ID = HDF5ID, GroupName = GroupName, Name = adjustl(trim(Name)),     & 
                               Units  = Units, Array1D   = IArray1D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteTimeValues - ModuleDDC - ERR02'

       else if (DataType .EQ. H5T_NATIVE_REAL) then if1
        
            call HDF5WriteData( HDF5ID = HDF5ID, GroupName = GroupName, Name = adjustl(trim(Name)),    & 
                                Units  = Units, Array1D   = R4Array1D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteTimeValues - ModuleDDC - ERR03'

       else if (DataType .EQ. H5T_NATIVE_DOUBLE) then if1

            call HDF5WriteData( HDF5ID = HDF5ID, GroupName = GroupName, Name = adjustl(trim(Name)),    & 
                                Units  = Units, Array1D   = R8Array1D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteTimeValues - ModuleDDC - ERR04'

       end if if1
       
       call HDF5UpdateGenericRealAttribute(HDF5ID = HDF5ID, GroupName = GroupName,    &
                                            ItemName = adjustl(trim(Name)),                             &
                                            ItemType = H5G_DATASET_F,                                   &
                                            AttributeName = "Minimum" ,                                 &
                                            ValueReal = Minimum, STAT = STAT_CALL)
                                            
       call HDF5UpdateGenericRealAttribute(HDF5ID = HDF5ID, GroupName = GroupName,    &
                                            ItemName = adjustl(trim(Name)),                             &
                                            ItemType = H5G_DATASET_F,                                   &
                                            AttributeName = "Maximum" ,                                 &
                                            ValueReal = Maximum, STAT = STAT_CALL)

       !Writes everything to disk
       call HDF5FlushMemory (HDF5ID = HDF5ID, STAT = STAT_CALL)
       if (STAT_CALL /= SUCCESS_) stop 'WriteTimeValues - ModuleDDC - ERR05'

          
        !------------------------------------------------------------------------
    
    end subroutine WriteTimeValues  

        !--------------------------------------------------------------------------

    subroutine SetArrayLimits(  Rank, GroupName, WindowPosition, Dimensions,    &
                                LimitsArrayFactor, LimitsArray)

        !Arguments-------------------------------------------------------------
        character(*), intent(IN)                                :: GroupName
        integer, intent(IN)                                     :: Rank
        integer, dimension(4), intent(IN)                       :: WindowPosition
        integer, dimension(4), intent(IN)                       :: Dimensions
        integer, intent(IN)                                     :: LimitsArrayFactor
        integer, dimension(6), intent(OUT)                      :: LimitsArray
        
        !Local-------------------------------------------------------------------
        
        
        !------------------------------------------------------------------------
                
if1 :   if (Rank .EQ. 1 .AND. adjustl(trim(GroupName)) .NE. "/Time") then
                    
            LimitsArray(1)= WindowPosition(1)
            LimitsArray(2)= WindowPosition(2)
                
        elseif (Rank .EQ. 2) then if1

            LimitsArray(1)= WindowPosition(1)
            LimitsArray(2)= WindowPosition(2)+LimitsArrayFactor
            LimitsArray(3)= WindowPosition(3)
            LimitsArray(4)= WindowPosition(4)+LimitsArrayFactor
                        
        elseif (Rank .EQ. 3) then if1
    
            LimitsArray(1)= WindowPosition(1)
            LimitsArray(2)= WindowPosition(2)
            LimitsArray(3)= WindowPosition(3)
            LimitsArray(4)= WindowPosition(4)
            LimitsArray(5)= 1
            LimitsArray(6)= Dimensions(3)
    
        endif if1
          
        !------------------------------------------------------------------------
    
    end subroutine SetArrayLimits  

    !--------------------------------------------------------------------------

    subroutine ReadHDFDataSetHyperslab(ObjHDF5_In, GroupName, obj_name,      &
                              Rank, InnerWindow, WindowFrame, DataType,      &
                              IArray1D, IArray2D, IArray3D,                  &
                              R4Array1D, R4Array2D, R4Array3D,               &
                              R8Array1D, R8Array2D, R8Array3D)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                     :: ObjHDF5_In
        character(*), intent(IN)                                :: GroupName
        character(*), intent(IN)                                :: obj_name 
        integer, intent(IN)                                     :: Rank
        integer, dimension(6), intent(IN)                       :: InnerWindow 
        integer, dimension(4), intent(IN)                       :: WindowFrame 
        integer(HID_T), intent(IN)                              :: DataType
        integer, dimension(:),       pointer, intent(OUT)       :: IArray1D
        integer, dimension(:, :),    pointer, intent(OUT)       :: IArray2D
        integer, dimension(:, :, :), pointer, intent(OUT)       :: IArray3D
        real(4), dimension(:),       pointer, intent(OUT)       :: R4Array1D
        real(4), dimension(:, :),    pointer, intent(OUT)       :: R4Array2D
        real(4), dimension(:, :, :), pointer, intent(OUT)       :: R4Array3D
        real(8), dimension(:),       pointer, intent(OUT)       :: R8Array1D
        real(8), dimension(:, :),    pointer, intent(OUT)       :: R8Array2D
        real(8), dimension(:, :, :), pointer, intent(OUT)       :: R8Array3D
        
        !Local-------------------------------------------------------------------
        integer                                                 :: STAT_CALL  
        integer                                                 :: ILB
        integer                                                 :: IUB  
        integer                                                 :: JLB
        integer                                                 :: JUB  
        integer                                                 :: KLB
        integer                                                 :: KUB     
        
        !------------------------------------------------------------------------

        STAT_CALL    = NULL_INT
        
        ILB = InnerWindow(1) - WindowFrame(1) + 1 
        IUB = InnerWindow(2) - WindowFrame(1) + 1
        JLB = InnerWindow(3) - WindowFrame(3) + 1
        JUB = InnerWindow(4) - WindowFrame(3) + 1
        KLB = InnerWindow(5)
        KUB = InnerWindow(6)

if1 :   if (Rank .EQ. 1) then

            call HDF5SetLimits  (ObjHDF5_In, ILB, IUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDFDataSetHyperslab - ModuleDDC - ERR01'

if12 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                allocate (IArray1D(ILB:IUB))
                
                call HDF5ReadWindow(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                    Array1D = IArray1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSetHyperslab - ModuleDDC - ERR02'
            
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if12
            
                allocate (R4Array1D(ILB:IUB))
                
                call HDF5ReadWindow(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                    Array1D = R4Array1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSetHyperslab - ModuleDDC - ERR03'
            
            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if12
            
                allocate (R8Array1D(ILB:IUB))
              
                call HDF5ReadWindow(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                    Array1D = R8Array1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSetHyperslab - ModuleDDC - ERR04'
            
            endif if12

        elseif (Rank .EQ. 2) then if1
        
            call HDF5SetLimits  (ObjHDF5_In, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDFDataSetHyperslab - ModuleDDC - ERR05'
        
if13 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                allocate (IArray2D(ILB:IUB,JLB:JUB))
                          
                call HDF5ReadWindow(ObjHDF5_In, GroupName, adjustl(trim(obj_name)),  & 
                                    Array2D = IArray2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSetHyperslab - ModuleDDC - ERR06'
                
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if13
            
                allocate (R4Array2D(ILB:IUB,JLB:JUB))
                
                call HDF5ReadWindow(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                    Array2D = R4Array2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSetHyperslab - ModuleDDC - ERR07'

            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if13
            
                allocate (R8Array2D(ILB:IUB,JLB:JUB))
                          
                call HDF5ReadWindow(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                    Array2D = R8Array2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSetHyperslab - ModuleDDC - ERR08'

            endif if13
        
        elseif (Rank .EQ. 3) then if1
        
            call HDF5SetLimits  (ObjHDF5_In, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDFDataSetHyperslab - ModuleDDC - ERR09'
        
if14 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                allocate (IArray3D(ILB:IUB,JLB:JUB,KLB:KUB))
                
                call HDF5ReadWindow(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                    Array3D = IArray3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSetHyperslab - ModuleDDC - ERR10'
                
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if14
            
                allocate (R4Array3D(ILB:IUB,JLB:JUB,KLB:KUB))
                          
                call HDF5ReadWindow(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                    Array3D = R4Array3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSetHyperslab - ModuleDDC - ERR11'
                
            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if14
            
                allocate (R8Array3D(ILB:IUB,JLB:JUB,KLB:KUB))
                          
                call HDF5ReadWindow(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                    Array3D = R8Array3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSetHyperslab - ModuleDDC - ERR12'
                
            endif if14
        
        endif if1
          
        !------------------------------------------------------------------------
    
    end subroutine ReadHDFDataSetHyperslab  
      
        !--------------------------------------------------------------------------

    subroutine WriteHDFDataSet (ObjHDF5_Out, GroupName, obj_name,                   &
                                Rank, DataType, Minimum, Maximum,                   &
                                IArray1D, IArray2D, IArray3D,                       &
                                R4Array1D, R4Array2D, R4Array3D,                    &
                                R8Array1D, R8Array2D, R8Array3D)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                     :: ObjHDF5_Out
        character(*), intent(IN)                                :: GroupName
        character(*), intent(IN)                                :: obj_name 
        integer, intent(IN)                                     :: Rank
        integer(HID_T), intent(IN)                              :: DataType
        real, intent(IN)                                        :: Minimum
        real, intent(IN)                                        :: Maximum
        
        !Local-------------------------------------------------------------------
        integer, dimension(:),       pointer                    :: IArray1D
        integer, dimension(:, :),    pointer                    :: IArray2D
        integer, dimension(:, :, :), pointer                    :: IArray3D
        real(4), dimension(:),       pointer                    :: R4Array1D
        real(4), dimension(:, :),    pointer                    :: R4Array2D
        real(4), dimension(:, :, :), pointer                    :: R4Array3D
        real(8), dimension(:),       pointer                    :: R8Array1D
        real(8), dimension(:, :),    pointer                    :: R8Array2D
        real(8), dimension(:, :, :), pointer                    :: R8Array3D
        integer                                                 :: STAT_CALL    
        
        !------------------------------------------------------------------------
        
        STAT_CALL = NULL_INT
        
if1 :   if (Rank .EQ. 1) then

if12 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                call HDF5WriteWindow(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                     Array1D = IArray1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR01'
                
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if12
            
                call HDF5WriteWindow(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                     Array1D = R4Array1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR02'

            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if12
            
                call HDF5WriteWindow(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                     Array1D = R8Array1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR03'

            endif if12

        elseif (Rank .EQ. 2) then if1
        
if13 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 
                          
                call HDF5WriteWindow(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                     Array2D = IArray2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR04'

            elseif (DataType .EQ. H5T_NATIVE_REAL) then if13
            
                call HDF5WriteWindow(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                     Array2D = R4Array2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR05'

            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if13
            
                call HDF5WriteWindow(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                     Array2D = R8Array2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR06'
            
            endif if13
        
        elseif (Rank .EQ. 3) then if1
        
if14 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 
                 
                call HDF5WriteWindow(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                     Array3D = IArray3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR07'
            
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if14
            
                call HDF5WriteWindow(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                     Array3D = R4Array3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR08'

            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if14
            
                call HDF5WriteWindow(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                     Array3D = R8Array3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR09'
            
            endif if14
        
        endif if1
        
        call HDF5UpdateGenericRealAttribute(HDF5ID = ObjHDF5_Out, GroupName = GroupName,    &
                                            ItemName = adjustl(trim(obj_name)),                             &
                                            ItemType = H5G_DATASET_F,                                   &
                                            AttributeName = "Minimum" ,                                 &
                                            ValueReal = Minimum, STAT = STAT_CALL)
                                            
        call HDF5UpdateGenericRealAttribute(HDF5ID = ObjHDF5_Out, GroupName = GroupName,    &
                                            ItemName = adjustl(trim(obj_name)),                             &
                                            ItemType = H5G_DATASET_F,                                   &
                                            AttributeName = "Maximum" ,                                 &
                                            ValueReal = Maximum, STAT = STAT_CALL)
        
        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5_Out, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR10'
          
        !------------------------------------------------------------------------
    
    end subroutine WriteHDFDataSet
    
        !--------------------------------------------------------------------------

    recursive subroutine DefineHDF5Limits  (ObjHDF5_Out, GroupName, obj_name,           &
                                            Rank, LimitsArray, WindowFrame, DomainSize, &
                                            InnerWindow, LimitsArrayFactor)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                     :: ObjHDF5_Out
        character(*), intent(IN)                                :: GroupName
        character(*), intent(IN)                                :: obj_name 
        integer, intent(IN)                                     :: Rank
        integer, intent(IN)                                     :: LimitsArrayFactor
        integer, dimension(6), intent(IN)                       :: LimitsArray
        integer, dimension(4), intent(IN)                       :: WindowFrame
        integer, dimension(4), intent(IN)                       :: DomainSize
        integer, dimension(4), intent(IN)                       :: InnerWindow
        
        !Local-------------------------------------------------------------------
        integer                                                 :: STAT_CALL
        integer                                                 :: ILB
        integer                                                 :: IUB  
        integer                                                 :: JLB
        integer                                                 :: JUB  
        integer                                                 :: KLB
        integer                                                 :: KUB    
        
        !------------------------------------------------------------------------
        
if1 :   if (DomainSize(1) .GT. 1) then

            ILB = InnerWindow(1)-DomainSize(1)+1
            IUB = InnerWindow(2)-DomainSize(1)+1+LimitsArrayFactor
            JLB = InnerWindow(3)-DomainSize(3)+1
            JUB = InnerWindow(4)-DomainSize(3)+1+LimitsArrayFactor
            KLB = LimitsArray(5)
            KUB = LimitsArray(6)
            
        else if1
        
            ILB = LimitsArray(1)
            IUB = LimitsArray(2)  
            JLB = LimitsArray(3)
            JUB = LimitsArray(4)  
            KLB = LimitsArray(5) 
            KUB = LimitsArray(6)
        
        endif if1
      
        call HDF5SetLimits( HDF5ID = ObjHDF5_Out,       & 
                            ILB = ILB,   IUB = IUB,     &
                            JLB = JLB,   JUB = JUB,     &
                            KLB = KLB,   KUB = KUB,     &
                            STAT   = STAT_CALL) 
       if (STAT_CALL /= SUCCESS_) stop 'DefineHDF5Limits - ModuleDDC - ERR01'
        
        !------------------------------------------------------------------------

    end subroutine DefineHDF5Limits
    
        !--------------------------------------------------------------------------

    recursive subroutine CreateGlobalSize  (ObjHDF5_Out, GroupName, obj_name,     &
                                            Units, Rank, DomainSize,              &
                                            dims, DataType)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                     :: ObjHDF5_Out
        character(*), intent(IN)                                :: GroupName
        character(*), intent(IN)                                :: obj_name 
        integer, intent(IN)                                     :: Rank
        integer, dimension(4), intent(IN)                       :: DomainSize
        integer, dimension(7), intent(IN)                       :: dims 
        integer(HID_T), intent(IN)                              :: DataType
        character(*), intent(IN)                                :: Units   
        
        !Local-------------------------------------------------------------------
        integer, dimension(:),       pointer                    :: IArray1D
        integer, dimension(:, :),    pointer                    :: IArray2D
        integer, dimension(:, :, :), pointer                    :: IArray3D
        real(4), dimension(:),       pointer                    :: R4Array1D
        real(4), dimension(:, :),    pointer                    :: R4Array2D
        real(4), dimension(:, :, :), pointer                    :: R4Array3D
        real(8), dimension(:),       pointer                    :: R8Array1D
        real(8), dimension(:, :),    pointer                    :: R8Array2D
        real(8), dimension(:, :, :), pointer                    :: R8Array3D
        integer                                                 :: STAT_CALL    
        integer                                                 :: LimitsArrayFactor
        integer                                                 :: ILB
        integer                                                 :: IUB  
        integer                                                 :: JLB
        integer                                                 :: JUB  
        integer                                                 :: KLB
        integer                                                 :: KUB
        
        !------------------------------------------------------------------------
        
        LimitsArrayFactor   = 0
        STAT_CALL           = NULL_INT
        
if11 :  if (adjustl(trim(obj_name)) .EQ. "ConnectionX" .OR. adjustl(trim(obj_name)) .EQ. "ConnectionY" .OR.  &
            adjustl(trim(obj_name)) .EQ. "Latitude" .OR. adjustl(trim(obj_name)) .EQ. "Longitude") then
                    
            LimitsArrayFactor = 1
                    
        endif if11
        
if21 :  if (Rank .EQ. 1) then

if22 :      if (DomainSize(1) .GT. 1) then
            
                ILB = DomainSize(1)-DomainSize(1)+1
                IUB = DomainSize(2)-DomainSize(1)+1
                JLB = 1
                JUB = 1  
                KLB = 1 
                KUB = 1
            
            else if22
                      
                ILB= DomainSize(1)
                IUB= DomainSize(2)
                JLB = 1
                JUB = 1  
                KLB = 1 
                KUB = 1 
            
            endif if22
                
        elseif (Rank .EQ. 2) then if21
            
if23 :      if (DomainSize(1) .GT. 1) then
            
                ILB = DomainSize(1)-DomainSize(1)+1
                IUB = DomainSize(2)-DomainSize(1)+1+LimitsArrayFactor
                JLB = DomainSize(3)-DomainSize(3)+1
                JUB = DomainSize(4)-DomainSize(3)+1+LimitsArrayFactor 
                KLB = 1 
                KUB = 1
            
            else if23    
            
                ILB = DomainSize(1)
                IUB = DomainSize(2)+LimitsArrayFactor
                JLB = DomainSize(3)
                JUB = DomainSize(4)+LimitsArrayFactor
                KLB = 1 
                KUB = 1
                
            endif if23     
                
        elseif (Rank .EQ. 3) then if21
        
if24 :      if (DomainSize(1) .GT. 1) then
            
                ILB = DomainSize(1)-DomainSize(1)+1
                IUB = DomainSize(2)-DomainSize(1)+1+LimitsArrayFactor
                JLB = DomainSize(3)-DomainSize(3)+1
                JUB = DomainSize(4)-DomainSize(3)+1+LimitsArrayFactor 
                KLB = 1 
                KUB = dims(3)
            
            else if24  
                
                ILB = DomainSize(1)
                IUB = DomainSize(2)
                JLB = DomainSize(3)
                JUB = DomainSize(4)
                KLB = 1
                KUB = dims(3)
            
            endif if24
                
        endif if21
        
        call HDF5SetLimits (HDF5ID = ObjHDF5_Out,                                   & 
                            ILB = ILB,   IUB = IUB,                                 &
                            JLB = JLB,   JUB = JUB,                                 &
                            KLB = KLB,   KUB = KUB,                                 &
                            STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateGlobalSize - ModuleDDC - ERR01'
        
if1 :   if (Rank .EQ. 1) then

if12 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                allocate (IArray1D(ILB:IUB))
                IArray1D(:) = null_int
            
                call HDF5WriteData(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                   Units, Array1D = IArray1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CreateGlobalSize - ModuleDDC - ERR02'

                deallocate(IArray1D)
            
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if12
            
                allocate (R4Array1D(ILB:IUB))
                R4Array1D(:) = null_real
                
                call HDF5WriteData(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                   Units, Array1D = R4Array1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CreateGlobalSize - ModuleDDC - ERR03'

                deallocate(R4Array1D)
            
            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if12
            
                allocate (R8Array1D(ILB:IUB))
                R8Array1D(:) = null_real
             
                call HDF5WriteData(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                   Units, Array1D = R8Array1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CreateGlobalSize - ModuleDDC - ERR04'

                deallocate(R8Array1D)
            
            endif if12

        elseif (Rank .EQ. 2) then if1
        
if13 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                allocate (IArray2D (ILB:IUB,                                        &
                                    JLB:JUB))
                IArray2D(:,:) = null_int
                         
                call HDF5WriteData(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                  Units, Array2D = IArray2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CreateGlobalSize - ModuleDDC - ERR05'

                deallocate(IArray2D)
            
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if13
            
                allocate (R4Array2D(ILB:IUB,                                        &
                                    JLB:JUB))
                R4Array2D(:,:) = null_real
                          
                call HDF5WriteData(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                  Units, Array2D = R4Array2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CreateGlobalSize - ModuleDDC - ERR06'

                deallocate(R4Array2D)                          

            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if13
            
                allocate (R8Array2D(ILB:IUB,                                       &
                                    JLB:JUB))
                R8Array2D(:,:) = null_real
                          
                call HDF5WriteData(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                  Units, Array2D = R8Array2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CreateGlobalSize - ModuleDDC - ERR07'

                deallocate(R8Array2D) 
            
            endif if13
        
        elseif (Rank .EQ. 3) then if1
        
if14 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                allocate (IArray3D (ILB:IUB,                                         &
                                    JLB:JUB,                                         &
                                    kLB:kUB))
                IArray3D(:,:,:) = null_int
                          
                call HDF5WriteData(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                   Units, Array3D = IArray3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CreateGlobalSize - ModuleDDC - ERR08'
                            
                deallocate(IArray3D)
            
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if14
            
                allocate (R4Array3D(ILB:IUB,                                         &
                                    JLB:JUB,                                         &
                                    kLB:kUB))
                R4Array3D(:,:,:) = null_real
                          
                call HDF5WriteData(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                   Units, Array3D = R4Array3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CreateGlobalSize - ModuleDDC - ERR09'
                            
                deallocate(R4Array3D)                          
            
            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if14
            
                allocate (R8Array3D(ILB:IUB,                                         &
                                    JLB:JUB,                                         &
                                    kLB:kUB))
                R8Array3D(:,:,:) = null_real
                          
                call HDF5WriteData(ObjHDF5_Out, GroupName, adjustl(trim(obj_name)), & 
                                   Units, Array3D = R8Array3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CreateGlobalSize - ModuleDDC - ERR10'
                            
                deallocate(R8Array3D)                             
            
            endif if14
        
        endif if1
        
        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5_Out, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CreateGlobalSize - ModuleDDC - ERR11'
          
        !------------------------------------------------------------------------

    end subroutine CreateGlobalSize
    
    !--------------------------------------------------------------------------

    subroutine ReadHDFDataSet(ObjHDF5_In, GroupName, obj_name,               &
                              Rank, dims, DataType,                          &
                              IArray1D, IArray2D, IArray3D,                  &
                              R4Array1D, R4Array2D, R4Array3D,               &
                              R8Array1D, R8Array2D, R8Array3D)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                     :: ObjHDF5_In
        character(*), intent(IN)                                :: GroupName
        character(*), intent(IN)                                :: obj_name 
        integer, intent(IN)                                     :: Rank
        integer, dimension(7), intent(IN)                       :: dims 
        integer(HID_T), intent(IN)                              :: DataType
        integer, dimension(:),       pointer, intent(OUT)       :: IArray1D
        integer, dimension(:, :),    pointer, intent(OUT)       :: IArray2D
        integer, dimension(:, :, :), pointer, intent(OUT)       :: IArray3D
        real(4), dimension(:),       pointer, intent(OUT)       :: R4Array1D
        real(4), dimension(:, :),    pointer, intent(OUT)       :: R4Array2D
        real(4), dimension(:, :, :), pointer, intent(OUT)       :: R4Array3D
        real(8), dimension(:),       pointer, intent(OUT)       :: R8Array1D
        real(8), dimension(:, :),    pointer, intent(OUT)       :: R8Array2D
        real(8), dimension(:, :, :), pointer, intent(OUT)       :: R8Array3D
        
        !Local-------------------------------------------------------------------
        integer                                                 :: STAT_CALL   
        
        !------------------------------------------------------------------------

        STAT_CALL    = NULL_INT

if1 :   if (Rank .EQ. 1) then

            call HDF5SetLimits  (ObjHDF5_In, 1, dims(1), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDFDataSet - ModuleDDC - ERR01'

if12 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                allocate (IArray1D(1:dims(1)))
                
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array1D = IArray1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR02'
            
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if12
            
                allocate (R4Array1D(1:dims(1)))
                
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array1D = R4Array1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR03'
            
            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if12
            
                allocate (R8Array1D(1:dims(1)))
              
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array1D = R8Array1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR04'
            
            endif if12

        elseif (Rank .EQ. 2) then if1
        
            call HDF5SetLimits  (ObjHDF5_In, 1, dims(1),1, dims(2), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDFDataSet - ModuleDDC - ERR05'
        
if13 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                allocate (IArray2D(1:dims(1),1:dims(2)))
                          
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)),  & 
                                  Array2D = IArray2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR06'
                
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if13
            
                allocate (R4Array2D(1:dims(1),1:dims(2)))
                
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array2D = R4Array2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR07'

            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if13
            
                allocate (R8Array2D(1:dims(1),1:dims(2)))
                          
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array2D = R8Array2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR08'

            endif if13
        
        elseif (Rank .EQ. 3) then if1
        
            call HDF5SetLimits  (ObjHDF5_In, 1, dims(1),1, dims(2),1,dims(3), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDFDataSet - ModuleDDC - ERR09'
        
if14 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                allocate (IArray3D(1:dims(1),1:dims(2),1:dims(3)))
                
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array3D = IArray3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR10'
                
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if14
            
                allocate (R4Array3D(1:dims(1),1:dims(2),1:dims(3)))
                          
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array3D = R4Array3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR11'
                
            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if14
            
                allocate (R8Array3D(1:dims(1),1:dims(2),1:dims(3)))
                          
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array3D = R8Array3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR12'
                
            endif if14
        
        endif if1
          
        !------------------------------------------------------------------------
    
    end subroutine ReadHDFDataSet
    
    !--------------------------------------------------------------------------
    
    recursive subroutine ScanDecomposedFiles(iHDFFile, hash_map_out, hash_map_in, Directory)

        !Arguments-------------------------------------------------------------
        character(PathLength)                                   :: Directory
        integer, intent(IN)                                     :: iHDFFile
        type(T_HashTable), pointer                              :: hash_map_out
        type(T_HashTable), pointer                              :: hash_map_in

        !Local-----------------------------------------------------------------
        integer                                                 :: iUnderScore      
        integer                                                 :: IDOut            
        integer                                                 :: iUnderScore2     
        integer                                                 :: STAT_            
        integer                                                 :: STAT_CALL        
        character(PathLength)                                   :: HDFinfile        
        character(PathLength)                                   :: MPIResultsFile   
        character(PathLength)                                   :: ConsolidatedFile 
        integer                                                 :: ObjHDF5_Out 
        integer                                                 :: FirstTime     
               
        !------------------------------------------------------------------------

        FirstTime        =  1
        iUnderScore      =  NULL_INT
        IDOut            =  NULL_INT
        iUnderScore2     =  NULL_INT
        STAT_            =  NULL_INT
        STAT_CALL        =  NULL_INT
        HDFinfile        =  null_str
        MPIResultsFile   =  null_str
        ConsolidatedFile =  null_str
        ObjHDF5_Out      =  NULL_INT


        read(unit = iHDFFile, fmt='(a256)', IOSTAT = STAT_CALL) HDFinfile
if1 :   if (STAT_CALL .EQ. SUCCESS_) then

            MPIResultsFile = trim(adjustl(HDFinfile))

            !Checks if it is a Domain Decomposition HDF5 results file type
            iUnderScore = scan(MPIResultsFile,"MPI_")
if2 :       if (iUnderScore > 0) then
                iUnderScore2 = scan(MPIResultsFile(5:),"_")

                
if3 :           if (iUnderScore2 > 0) then
                    ConsolidatedFile = adjustl(trim(Directory)) // '\' // MPIResultsFile(5+iUnderScore2:)
                    
                    IDOut = hash_get(hash_map_out, ConsolidatedFile)
                    
if4 :               if (IDOut <0) then

                        IDOut = NULL_INT
                        
                        call GetHDF5ReadWriteFileAccessCode(ConsolidatedFile = ConsolidatedFile, &
                                                            IDOut           = IDOut,             &
                                                            ObjHDF5_Out     = ObjHDF5_Out)

                        call hash_set(hash_map_out, key = ConsolidatedFile, value_ = IDOut)
                        
                        call hash_setObjID(hash_map_out, key = ConsolidatedFile, ObjID = ObjHDF5_Out)
                        
                    else if4
                    
                        ObjHDF5_Out = hash_getObjID(Me = hash_map_out, key = ConsolidatedFile)
                        
                    endif if4

                    call hash_set(hash_map_in, key = adjustl(trim(Directory)) // '\' // HDFinfile, value_ = IDOut)
                
                    call hash_setObjID(hash_map_in, key = adjustl(trim(Directory)) // '\' // HDFinfile, ObjID = ObjHDF5_Out)
                
                    call GetDomainDecomposition(HDFFile = adjustl(trim(Directory)) // '\' // HDFinfile,   &
                                                hash_map_in = hash_map_in,                                &
                                                FirstTime = FirstTime)
                                    
                endif if3
            endif if2
            
            call ScanDecomposedFiles(iHDFFile = iHDFFile,                       &
                                     hash_map_out = hash_map_out,               &
                                     hash_map_in = hash_map_in,                 &
                                     Directory = Directory)
        endif if1

        !------------------------------------------------------------------------

    end subroutine ScanDecomposedFiles
    
    !--------------------------------------------------------------------------

    subroutine GetHDF5ReadWriteFileAccessCode(ConsolidatedFile,IDOut,ObjHDF5_Out)

        !Arguments-------------------------------------------------------------
        character(PathLength), intent(IN)               :: ConsolidatedFile 
        integer, intent(OUT)                            :: IDOut              
        integer, intent(OUT)                            :: ObjHDF5_Out              
                
        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL        
        integer                                         :: HDF5_READWRITE   
        integer                                         :: HDF5_CREATE      
        integer                                         :: Status           
        logical                                         :: HDF5FileExists
        
        !------------------------------------------------------------------------

        STAT_CALL        = NULL_INT
        HDF5_READWRITE   = NULL_INT
        HDF5_CREATE      = NULL_INT
        Status           = NULL_INT
        ObjHDF5_Out      = NULL_INT

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
                
        call ConstructHDF5 (HDF5ID   = ObjHDF5_Out,                             &
                            FileName = adjustl(trim(ConsolidatedFile)),         &
                            Access   = HDF5_CREATE,                             &
                            STAT     = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ReadWriteFileAccessCode - ModuleDDC - ERR01'
        
        call GetHDF5FileID (HDF5ID = ObjHDF5_Out,                               &
                            FileID = IDOut,                                     &
                            STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ReadWriteFileAccessCode - ModuleDDC - ERR02'
        
        !------------------------------------------------------------------------

    end subroutine GetHDF5ReadWriteFileAccessCode
    
    !--------------------------------------------------------------------------
    
    recursive subroutine GetDomainDecomposition(HDFFile,hash_map_in,FirstTime)

        !Arguments-------------------------------------------------------------
        CHARACTER(PathLength)                       :: HDFFile
        type(T_HashTable), pointer                  :: hash_map_in
        integer                                     :: FirstTime            
            
        !Local-------------------------------------------------------------------
        CHARACTER(PathLength)                       :: DecompositionGroup 
        integer                                     :: Group               
        character(PathLength)                       :: GroupName           
        integer                                     :: IDIn               
        integer                                     :: ObjHDF5_In         
        integer, pointer,   dimension(:)            :: DataVal1D   
        integer                                     :: STAT_CALL  
        integer, dimension(4)                       :: DomainSize         
        integer, dimension(4)                       :: WindowPosition
        integer                                     :: i 

        !------------------------------------------------------------------------

        IDIn       = NULL_INT
        ObjHDF5_In = NULL_INT
        Group      = NULL_INT
        GroupName  = NULL_STR

if1 :   if(FirstTime .EQ. 1) then
            DecompositionGroup = "/Grid/Decomposition/Global/"            
        elseif (FirstTime .EQ. 2) then if1
            DecompositionGroup = "/Grid/Decomposition/InnerMapping/"
        elseif (FirstTime .EQ. 3) then if1
            DecompositionGroup = "/Grid/Decomposition/Mapping/"
        endif if1
        
        call GetHDF5ReadFileAccessCode(HDFFile = HDFFile,               &
                                       IDIn = IDIn,                     &
                                       ObjHDF5_In = ObjHDF5_In)

        call GetMappingValues(ObjHDF5_In = ObjHDF5_In,                  &
                              IDIn = IDIn,                              &
                              GroupName = DecompositionGroup,           &
                              DataVal1D = DataVal1D) 
         
if2 :   if(FirstTime .EQ. 1) then
            call hash_setDomainSize(Me = hash_map_in,                   &
                                    key = HDFFile,                      &  
                                    DomainSize = DataVal1D)            
        
            call KillHDF5(HDF5ID   = ObjHDF5_In,                                                  &
                          STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetDomainDecomposition - ModuleDDC - ERR01'
       
        elseif (FirstTime .EQ. 2) then if2

            call hash_setWindowPosition(Me = hash_map_in,               &
                                        key = HDFFile,                  &  
                                        WindowPosition = DataVal1D)
                                        
        elseif (FirstTime .EQ. 3) then if2
        
            call hash_setWindowFrame(Me = hash_map_in,                  &
                                        key = HDFFile,                  &  
                                        WindowFrame = DataVal1D)
                                        
            call KillHDF5(HDF5ID   = ObjHDF5_In,                                                  &
                          STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetDomainDecomposition - ModuleDDC - ERR02'
        
        endif if2
        
if3 :   if (FirstTime .LE. 2) then

            call GetDomainDecomposition(HDFFile = HDFFile,   &
                                            hash_map_in = hash_map_in,                                &
                                            FirstTime = FirstTime+1)
        endif if3
        
        !------------------------------------------------------------------------

        end subroutine GetDomainDecomposition
        !--------------------------------------------------------------------------

        subroutine GetMappingValues(ObjHDF5_In, IDIn, GroupName, DataVal1D)

       !Arguments-------------------------------------------------------------
        integer(HID_T)                                      :: IDIn
        character(len=*)                                    :: GroupName
        integer                                             :: ObjHDF5_In
        integer(4), pointer, dimension(:)                   :: DataVal1D 
        
        !Local-----------------------------------------------------------------
        character(StringLength)                             :: obj_name
        integer                                             :: idx
        integer(HID_T)                                      :: GroupType
        integer                                             :: STAT_CALL
        
        !------------------------------------------------------------------------

        idx      = 1
        obj_name = adjustl(trim("ILB_IUB_JLB_JUB"))

        !Gets information about the group
        call GetHDF5ObjectInfo (HDF5ID = ObjHDF5_In, FatherGroupName = adjustl(trim(GroupName)),  &
                                GroupPosition = idx, GroupName = obj_name,GroupType = GroupType,      &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetMappingValues - ModuleDDC - ERR01'

if1 :   if (GroupType == H5G_DATASET_F) then

            allocate(DataVal1D(1:4))
                
            call HDF5SetLimits  (HDF5ID = ObjHDF5_In, ILB = 1, IUB = 4, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'GetMappingValues - ModuleDDC - ERR02'
                
            call HDF5ReadData(HDF5ID = ObjHDF5_In, GroupName = adjustl(trim(GroupName)), Name = obj_name,                          &
                                  Array1D = DataVal1D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetMappingValues - ModuleDDC - ERR03'

        elseif (GroupType ==H5G_GROUP_F) then

            STAT_CALL = -1 * H5G_GROUP_F
            if (STAT_CALL /= SUCCESS_) stop 'GetMappingValues - ModuleDDC - ERR04'

        endif if1
         
        !------------------------------------------------------------------------

    end subroutine GetMappingValues

    !--------------------------------------------------------------------------

     subroutine GetHDF5ReadFileAccessCode(HDFFile,IDIn, ObjHDF5_In)

        !Arguments-------------------------------------------------------------
        character(PathLength), intent(IN)                        :: HDFFile 
        integer, intent(OUT)                                     :: IDIn  
        integer, intent(OUT)                                     :: ObjHDF5_In                  
              
        !Local-------------------------------------------------------------------
        integer                                                  :: STAT_CALL       
        integer                                                  :: HDF5_READ       
        
        !------------------------------------------------------------------------

        STAT_CALL  = NULL_INT
        HDF5_READ  = NULL_INT
        ObjHDF5_In = NULL_INT

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
                
        call ConstructHDF5 (HDF5ID   = ObjHDF5_In,                             &
                            FileName = adjustl(trim(HDFFile)),                 &
                            Access   = HDF5_READ,                              &
                            STAT     = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ReadFileAccessCode - ModuleDDC - ERR01'
                            
        call GetHDF5FileID (HDF5ID = ObjHDF5_In,                               &
                            FileID = IDIn,                                     &
                            STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetHDF5ReadFileAccessCode - ModuleDDC - ERR02'  
          
        !------------------------------------------------------------------------

    end subroutine GetHDF5ReadFileAccessCode
        
    !--------------------------------------------------------------------------
    
    subroutine OpenDecomposedFiles(hash_map_out, hash_map_in, Directory, nbrModels)

        !Arguments-------------------------------------------------------------
        type(T_HashTable), pointer                  :: hash_map_out
        type(T_HashTable), pointer                  :: hash_map_in
        character(PathLength), intent(IN)           :: Directory
        integer, intent(IN)                         :: nbrModels

        !Local-----------------------------------------------------------------
        character(StringLength)                     :: DecomposedFiles  = NULL_STR

        !----------------------------------------------------------------------

        DecomposedFiles = adjustl(trim('DecomposedFiles'))

        !Starts looking for first DecomposedFiles: DecomposedFiles_0.dat
        call OpenDecomposedFiles2(hash_map_out    = hash_map_out,               &
                                  hash_map_in     = hash_map_in,                &
                                  Directory       = Directory,                  &
                                  DecomposedFiles = DecomposedFiles,            &
                                  iModel          = 0,                          &
                                  nbrModels       = nbrModels)

        !------------------------------------------------------------------------

    end subroutine OpenDecomposedFiles

    !--------------------------------------------------------------------------

    recursive subroutine OpenDecomposedFiles2(hash_map_out, hash_map_in, Directory, DecomposedFiles, iModel, nbrModels)

        !Arguments-------------------------------------------------------------
        type(T_HashTable), pointer                  :: hash_map_out
        type(T_HashTable), pointer                  :: hash_map_in
        character(StringLength), intent(IN)         :: DecomposedFiles
        character(PathLength), intent(IN)           :: Directory
        integer, intent(IN)                         :: iModel
        integer, intent(IN)                         :: nbrModels

        !Local-----------------------------------------------------------------
        logical                                     :: DDFileExists
        integer                                     :: STAT_CALL        = NULL_INT
        integer                                     :: iDDFile          = NULL_INT
        integer                                     :: ready_           = NULL_INT
        character(StringLength)                     :: Coment1          = NULL_STR
        character(StringLength)                     :: Coment2          = NULL_STR
        character(StringLength)                     :: DDFile           = NULL_STR
        character(StringLength)                     :: AuxString        = null_str

        !----------------------------------------------------------------------

if2 :   if (iModel .LE. nbrModels) then
            write (AuxString, '(i10)') iModel
            DDFile = adjustl(trim(                                              &
                         adjustl(trim(Directory))//'\MPI_'//                    &
                         adjustl(trim(AuxString))))//'_'//                      &
                         adjustl(trim(DecomposedFiles))//'.dat'
            inquire(file = DDFile, EXIST = DDFileExists)

if1 :       if (DDFileExists) then
                call UnitsManager(iDDFile, OPEN_FILE, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OpenDecomposedFiles2 - ModuleDDC - ERR01'

                open(UNIT = iDDFile, FILE = DDFile, status = 'OLD', IOSTAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OpenDecomposedFiles2 - ModuleDDC - ERR02'

!                read(unit=iDDFile, fmt=*) Coment1
!                read(unit=iDDFile, fmt=*) Coment2

                call ScanDecomposedFiles(iHDFFile      = iDDFile,               &
                                         hash_map_out  = hash_map_out,          &
                                         hash_map_in   = hash_map_in,           &
                                         Directory     = Directory)
                                         
                call UnitsManager(iDDFile, CLOSE_FILE, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OpenDecomposedFiles2 - ModuleDDC - ERR03'
            endif if1

            call OpenDecomposedFiles2(hash_map_out  = hash_map_out,             &
                                      hash_map_in   = hash_map_in,              &
                                      Directory       = Directory,              &
                                      DecomposedFiles = DecomposedFiles,        &
                                      iModel          = iModel + 1,             &
                                      nbrModels       = nbrModels)
        endif if2

        !------------------------------------------------------------------------

    end subroutine OpenDecomposedFiles2

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    integer function KillDDC(Me)

        !Arguments---------------------------------------------------------------
        type (T_DDC), pointer                       :: Me

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                             :: ready_
        integer                             :: STAT_        = NULL_INT

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        ready_ = Ready(Me)

cd1 :   if (ready_ .NE. OFF_ERR_) then

!            nUsers = DeassociateInstance(mDDC_,  Me%InstanceID)

!            if (nUsers == 0) then

                !Deallocates Instance
                call DeallocateInstance(Me)

                STAT_      = SUCCESS_

!            end if
        else
            STAT_ = ready_
        end if cd1

        KillDDC = STAT_

        !------------------------------------------------------------------------

    end function KillDDC

    !------------------------------------------------------------------------

    subroutine DeallocateInstance(Me)

        !Arguments-------------------------------------------------------------
        type (T_DDC), pointer                       :: Me

        !Local-----------------------------------------------------------------

        !------------------------------------------------------------------------

        if (associated(Me%DirectoryList)) call deallocateDirectoryList(Me%DirectoryList)

        deallocate (Me)
        nullify    (Me)

        !------------------------------------------------------------------------

    end subroutine DeallocateInstance

    !------------------------------------------------------------------------

    recursive subroutine deallocateDirectoryList (DirectoryList)

        !Arguments-------------------------------------------------------------
        type(T_DirectoryList), pointer              :: DirectoryList

        !Local-----------------------------------------------------------------

        !------------------------------------------------------------------------

if1 :   if (associated(DirectoryList%Next)) then
            call KillConsolidatedFiles(DirectoryList%hash_map_out)

            call KillHash_map           (DirectoryList%hash_map_out)
            nullify   (DirectoryList%hash_map_out)
            call KillHash_map           (DirectoryList%hash_map_in)
            nullify   (DirectoryList%hash_map_in)
            call deallocateDirectoryList(DirectoryList%Next)
        endif if1

        deallocate (DirectoryList)
        nullify    (DirectoryList)

        !------------------------------------------------------------------------

    end subroutine deallocateDirectoryList
    
    !------------------------------------------------------------------------

    recursive subroutine KillConsolidatedFiles (hash_map_out)

        !Arguments-------------------------------------------------------------
        type(T_HashTable), pointer                      :: hash_map_out
        
        !Local-----------------------------------------------------------------
        CHARACTER(StringLength)                         :: key

        !------------------------------------------------------------------------

if1 :   if (hash_get_first_exists(hash_map_out)) then

            key = hash_get_first_key(hash_map_out)
            
            call KillConsolidatedFiles2(hash_map_out, key)

        endif if1

        !------------------------------------------------------------------------

    end subroutine KillConsolidatedFiles
    
    !------------------------------------------------------------------------

    recursive subroutine KillConsolidatedFiles2 (hash_map_out, key)

        !Arguments-------------------------------------------------------------
        type(T_HashTable), pointer                      :: hash_map_out
        CHARACTER(*), INTENT(IN)                        :: key

        !Local-----------------------------------------------------------------
        integer                                         :: ObjHDF5_Out      = NULL_INT
        integer                                         :: STAT_CALL        = NULL_INT
        CHARACTER(StringLength)                         :: key2             = NULL_STR
        !------------------------------------------------------------------------

            ObjHDF5_Out = hash_getObjID(hash_map_out, key)
            
            call KillHDF5(HDF5ID   = ObjHDF5_Out,                                    &
                          STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillConsolidatedFiles2 - ModuleDDC - ERR001'

if1 :       if (hash_get_next_exists(hash_map_out, key)) then
                key2 = hash_get_next_key(hash_map_out, Key)
                call KillConsolidatedFiles2(hash_map_out, key2)
            endif if1
        !------------------------------------------------------------------------

    end subroutine KillConsolidatedFiles2

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    integer function Ready(Me)

        !Arguments-------------------------------------------------------------
        type (T_DDC), pointer                       :: Me
        integer                                     :: ready_ = NULL_INT

        !----------------------------------------------------------------------

!        nullify (Me)

!cd1:    if (ObjDDC_ID > 0) then
!            call LocateObjDDC (ObjDDC_ID)
!            ready_ = VerifyReadLock (mDDC_, Me%InstanceID)
!        else cd1
!            ready_ = OFF_ERR_
!        end if cd1

        Ready = ready_

        !----------------------------------------------------------------------

    end function Ready

    !--------------------------------------------------------------------------

    subroutine Read_Lock_DDC(Me)

        !Arguments---------------------------------------------------------------
        type (T_DDC), pointer                       :: Me

        !----------------------------------------------------------------------

!        ObjCollector(iModule, iInstance)%Read_Lock = ACTIVE
!        ObjCollector(iModule, iInstance)%readers   = ObjCollector(iModule, iInstance)%readers + 1

        !----------------------------------------------------------------------

    end subroutine Read_Lock_DDC

    !--------------------------------------------------------------------------

    subroutine Read_Unlock_DDC()

        !Arguments-------------------------------------------------------------

        !----------------------------------------------------------------------

        !Cannot Read_unlock if Instance is not Read_lock
!        if (.not. ObjCollector(iModule, iInstance)%Read_Lock) then
!            write(*, *) 'Number of Readers error.'
!            write(*, *) 'Routine Name    : ', trim(RoutineName)
!            write(*, *) 'Module  Name    : ', trim(MohidModules(iModule)%Name)
!            write(*, *) 'Instance  ID    : ', iInstance
!            stop        'Read_Unlock - ModuleGlobalData - ERR01'
!        end if

        !Decreases number of readers
!        ObjCollector(iModule, iInstance)%readers = ObjCollector(iModule, iInstance)%readers - 1

        !If number of reades equal to zero set Read_lock to IDLE
!        if (ObjCollector(iModule, iInstance)%readers == 0) ObjCollector(iModule, iInstance)%Read_Lock = IDLE

        !if Number of readers is negative, somethink is wrong
!        if (ObjCollector(iModule, iInstance)%readers < 0) then
!            write(*, *) 'Negative number of readers'
!            write(*, *) 'Routine Name    : ', trim(RoutineName)
!            write(*, *) 'Module  Name    : ', trim(MohidModules(iModule)%Name)
!            write(*, *) 'Instance  ID    : ', iInstance
!            stop        'Read_Unlock - ModuleGlobalData - ERR02'
!        end if

        !----------------------------------------------------------------------

    end subroutine Read_Unlock_DDC

    !--------------------------------------------------------------------------

    integer function VerifyReadLock_DDC()

        !Arguments-------------------------------------------------------------

        !------------------------------------------------------------------------

!        if  (ObjCollector(iModule, iInstance)%Read_Lock) then
!            VerifyReadLock = READ_LOCK_ERR_
!        else
            VerifyReadLock_DDC = IDLE_ERR_
!        end if

        !------------------------------------------------------------------------

    end function VerifyReadLock_DDC

    !-------------------------------------------------------------------------

end module ModuleDDC

