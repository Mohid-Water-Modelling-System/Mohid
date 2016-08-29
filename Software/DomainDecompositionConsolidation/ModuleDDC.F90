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
    private ::              ReadHDFDataSet
    private ::              AuxArrayAllocation
    private ::              GetLimitsHDF
    private ::              GetIndexMPIHDF
    private ::              GetIndexAUXHDF
    private ::              MergeArraysMPI
    private ::              MPIArrayDeallocate
    private ::              WriteHDFDataSetLimits
    private ::              WriteHDFDataSet
    private ::              AuxArrayDeallocate
    
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

    type T_MPIMatrixes
        real(4),    dimension(:    ), pointer :: MPIDataR4_1D
        real(4),    dimension(:,:  ), pointer :: MPIDataR4_2D
        real(4),    dimension(:,:,:), pointer :: MPIDataR4_3D
        real(8),    dimension(:    ), pointer :: MPIDataR8_1D
        real(8),    dimension(:,:  ), pointer :: MPIDataR8_2D
        real(8),    dimension(:,:,:), pointer :: MPIDataR8_3D
        integer,    dimension(:    ), pointer :: MPIDataI4_1D
        integer,    dimension(:,:  ), pointer :: MPIDataI4_2D
        integer,    dimension(:,:,:), pointer :: MPIDataI4_3D        
    end type T_MPIMatrixes
    
    type T_AuxMatrixes
        real(4),    dimension(:    ), pointer :: AuxDataR4_1D
        real(4),    dimension(:,:  ), pointer :: AuxDataR4_2D
        real(4),    dimension(:,:,:), pointer :: AuxDataR4_3D
        real(8),    dimension(:    ), pointer :: AuxDataR8_1D
        real(8),    dimension(:,:  ), pointer :: AuxDataR8_2D
        real(8),    dimension(:,:,:), pointer :: AuxDataR8_3D
        integer,    dimension(:    ), pointer :: AuxDataI4_1D
        integer,    dimension(:,:  ), pointer :: AuxDataI4_2D
        integer,    dimension(:,:,:), pointer :: AuxDataI4_3D        
    end type T_AuxMatrixes 
    
    public  :: T_DDC
    type       T_DDC
        private

        !nbrModels stores total number of models
        integer                                         :: nbrModels        = NULL_INT
        type(T_DirectoryList), pointer                  :: DirectoryList
    end type  T_DDC

    private :: T_DirectoryList
    type       T_DirectoryList
        character(PathLength)                           :: Directory        = NULL_STR
        type (T_HashTable), pointer                     :: hash_map_out !consolidated HDF5 output files
        type (T_HashTable), pointer                     :: hash_map_in  !decomposed HDF5 input files
        type (T_DirectoryList), pointer                 :: Next
        type (T_AuxMatrixes), pointer                   :: AuxMatrixes
        type (T_MPIMatrixes), pointer                   :: MPIMatrixes
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
        
        nullify (Me%DirectoryList)

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

        position  = scan(AuxString, "/", back = .true.)
        if (position == 0) then
            position = scan(AuxString, "\", back = .true.)
        endif
        if (position == 0) then
            ModelPath = "../res"
        else
            ModelPath  = AuxString(Level+1:position)//"res"
        endif
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
        
        allocate(NewDirectoryList%AuxMatrixes)
        allocate(NewDirectoryList%MPIMatrixes)
        
        nullify(NewDirectoryList%AuxMatrixes%AuxDataR4_1D)
        nullify(NewDirectoryList%AuxMatrixes%AuxDataR4_2D)
        nullify(NewDirectoryList%AuxMatrixes%AuxDataR4_3D)
        nullify(NewDirectoryList%AuxMatrixes%AuxDataR8_1D)
        nullify(NewDirectoryList%AuxMatrixes%AuxDataR8_2D)
        nullify(NewDirectoryList%AuxMatrixes%AuxDataR8_3D)
        nullify(NewDirectoryList%AuxMatrixes%AuxDataI4_1D)
        nullify(NewDirectoryList%AuxMatrixes%AuxDataI4_2D)
        nullify(NewDirectoryList%AuxMatrixes%AuxDataI4_3D)
        
        nullify(NewDirectoryList%MPIMatrixes%MPIDataR4_1D)
        nullify(NewDirectoryList%MPIMatrixes%MPIDataR4_2D)
        nullify(NewDirectoryList%MPIMatrixes%MPIDataR4_3D)
        nullify(NewDirectoryList%MPIMatrixes%MPIDataR8_1D)
        nullify(NewDirectoryList%MPIMatrixes%MPIDataR8_2D)
        nullify(NewDirectoryList%MPIMatrixes%MPIDataR8_3D)
        nullify(NewDirectoryList%MPIMatrixes%MPIDataI4_1D)
        nullify(NewDirectoryList%MPIMatrixes%MPIDataI4_2D)
        nullify(NewDirectoryList%MPIMatrixes%MPIDataI4_3D)
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
        CHARACTER(PathLength)                   :: FirstHDFFileIn
        CHARACTER(PathLength)                   :: FirstHDFFileOut
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

            FirstHDFFileIn = hash_get_first_key(DirectoryList%hash_map_in)
            FirstHDFFileOut = hash_get_first_key(DirectoryList%hash_map_out)
            
            call ScanFileList (HDFFileIn    = FirstHDFFileIn,                   &
                               hash_map_in  = DirectoryList%hash_map_in,        &
                               HDFFileOut   = FirstHDFFileOut,                  &
                               hash_map_out = DirectoryList%hash_map_out,       &
                               AuxMatrixes  = DirectoryList%AuxMatrixes,        &
                               MPIMatrixes  = DirectoryList%MPIMatrixes)
        endif if1

        if (associated(DirectoryList%Next)) call ScanDirectoryList(DirectoryList = DirectoryList%Next,  &
                                                                   nbrModels     = nbrModels)
!$OMP END SECTIONS NOWAIT
!$OMP END PARALLEL


        !------------------------------------------------------------------------

    end subroutine ScanDirectoryList

    !--------------------------------------------------------------------------

    recursive subroutine ScanFileList(HDFFileIn, hash_map_in, HDFFileOut,       &
                                      hash_map_out, AuxMatrixes, MPIMatrixes)

        !Arguments-------------------------------------------------------------
        CHARACTER(PathLength), intent(IN)                       :: HDFFileIn
        CHARACTER(PathLength), intent(IN)                       :: HDFFileOut
        type(T_HashTable), pointer                              :: hash_map_in 
        type(T_HashTable), pointer                              :: hash_map_out
        type(T_AuxMatrixes), pointer                            :: AuxMatrixes
        type(T_MPIMatrixes), pointer                            :: MPIMatrixes
         
        !Local--------------------------------------------------------------
        integer                                                 :: STAT_CALL 
        CHARACTER(PathLength)                                   :: HDFFileNext
        CHARACTER(PathLength)                                   :: HDFFileNextOut
        integer                                                 :: IDOut                   
        integer                                                 :: ObjHDF5_Out
        integer                                                 :: CountInFiles                   
        integer                                                 :: TotalFiles
        integer                                                 :: IDin                   
        integer                                                 :: ObjHDF5_In
        integer                                                 :: AuxIDOut                            
        integer                                                 :: AuxObjHDF5_Out
        CHARACTER(PathLength), dimension(:), pointer            :: FileArrayIn
        integer, dimension(:), pointer                          :: IdInArray
        integer, dimension(:), pointer                          :: ObjHDF5_InArray
        integer                                                 :: TotalFilesIn                   
        integer                                                 :: AuxLoop                   
        integer                                                 :: i
        character(StringLength)                                 :: GroupName  
        
        !------------------------------------------------------------------------
        
        allocate(FileArrayIn(1000))
        allocate(IdInArray(1000))
        allocate(ObjHDF5_InArray(1000))
        
        FileArrayIn     = null_str
        IdInArray       = null_int
        ObjHDF5_InArray = null_int
        TotalFilesIn    = null_int
        
        GroupName    = ""
        CountInFiles=1.
        TotalFiles = 0.
        AuxLoop = 0.
        HDFFileNext=HDFFileIn
        
        IDOut       = hash_get     (hash_map_out, HDFFileOut)
        ObjHDF5_Out = hash_getObjID(hash_map_out, HDFFileOut)
        
        AuxIDOut       = hash_get     (hash_map_in, HDFFileNext)
        AuxObjHDF5_Out = hash_getObjID(hash_map_in, HDFFileNext)
                                        
if1 :   if (AuxObjHDF5_Out == ObjHDF5_Out) then

            TotalFiles = TotalFiles + 1. 
        
            !Get access to first MPI file 
            call GetHDF5ReadFileAccessCode( HDFFile    = HDFFileNext,                       &
                                            IDIn       = IDIn,                              &
                                            ObjHDF5_In = ObjHDF5_In)
            
            FileArrayIn(TotalFiles)     = HDFFileNext
            IdInArray(TotalFiles)       = IDIn
            ObjHDF5_InArray(TotalFiles) = ObjHDF5_In
            TotalFilesIn                = TotalFiles
            
        endif if1
        
do1 :   do while(AuxLoop == 0)
            
if2 :       if (hash_get_next_exists(hash_map_in, HDFFileNext)) then

                HDFFileNext = hash_get_next_key(hash_map_in, HDFFileNext)
                
                AuxIDOut       = hash_get     (hash_map_in, HDFFileNext)
                AuxObjHDF5_Out = hash_getObjID(hash_map_in, HDFFileNext)
                                                
if3 :           if (AuxObjHDF5_Out == ObjHDF5_Out) then 
                
                    !Get access to the MPI file 
                    call GetHDF5ReadFileAccessCode( HDFFile    = HDFFileNext,               &
                                                    IDIn       = IDIn,                      &
                                                    ObjHDF5_In = ObjHDF5_In)
                
                    TotalFiles = TotalFiles + 1.
                    
                    FileArrayIn(TotalFiles)     = HDFFileNext
                    IdInArray(TotalFiles)       = IDIn
                    ObjHDF5_InArray(TotalFiles) = ObjHDF5_In
                    TotalFilesIn                = TotalFiles
                    
                endif if3
                
            else if2
            
                AuxLoop = 1
            
            endif if2
            
        enddo do1
        
        print *, adjustl(trim(HDFFileOut))
        print *, ""
        
        call WriteConsolidatedHDF(  IDOut           = IDOut,                            &
                                    ObjHDF5_Out     = ObjHDF5_Out,                      & 
                                    GroupName       = adjustl(trim(GroupName)),         &
                                    IdInArray       = IdInArray,                        &
                                    ObjHDF5_InArray = ObjHDF5_InArray,                  &
                                    FileArrayIn     = FileArrayIn,                      &
                                    NumberOfFiles   = TotalFilesIn,                     &
                                    hash_map_in     = hash_map_in,                      &
                                    AuxMatrixes     = AuxMatrixes,                      &
                                    MPIMatrixes     = MPIMatrixes)
                        
do2 :   do i=1,TotalFilesIn

            call KillHDF5(  HDF5ID   = ObjHDF5_InArray(i),              &
                            STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ScanFileList - ModuleDDC - ERR01'


        enddo do2 
        
        deallocate(FileArrayIn)
        deallocate(IdInArray)
        deallocate(ObjHDF5_InArray)      

if4 :   if (hash_get_next_exists(hash_map_out, HDFFileOut)) then

            HDFFileNextOut = hash_get_next_key(hash_map_out, HDFFileOut)
            call ScanFileList(  HDFFileIn    = HDFFileIn,                   &
                                hash_map_in  = hash_map_in,                 &
                                HDFFileOut   = HDFFileNextOut,              &
                                hash_map_out = hash_map_out,                &
                                AuxMatrixes  = AuxMatrixes,                 &
                                MPIMatrixes  = MPIMatrixes)
            
        endif if4
        
        !------------------------------------------------------------------------

    end subroutine ScanFileList

    !--------------------------------------------------------------------------

    recursive subroutine WriteConsolidatedHDF(  IDOut, ObjHDF5_Out, GroupName,              &
                                                IdInArray, ObjHDF5_InArray, FileArrayIn,    &
                                                NumberOfFiles, hash_map_in, AuxMatrixes,    &
                                                MPIMatrixes)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                     :: IDOut              
        integer, intent(IN)                                     :: ObjHDF5_Out
        character(*), intent(IN)                                :: GroupName  
        CHARACTER(*), dimension(:), pointer                     :: FileArrayIn
        integer, dimension(:), pointer                          :: IdInArray
        integer, dimension(:), pointer                          :: ObjHDF5_InArray
        integer                                                 :: NumberOfFiles
        type(T_HashTable), pointer                              :: hash_map_in
        type(T_AuxMatrixes), pointer                            :: AuxMatrixes        
        type(T_MPIMatrixes), pointer                            :: MPIMatrixes
        
        !Local-------------------------------------------------------------------
        integer                                                 :: nItems
        integer                                                 :: idx  
        integer                                                 :: STAT_CALL
        character(StringLength)                                 :: Name 
        character(StringLength)                                 :: NewGroupName 
        integer(HID_T)                                          :: GroupType
        integer                                                 :: i       
        integer                                                 :: Rank       
        character(StringLength)                                 :: Units
        integer, dimension(7)                                   :: Dimensions 
        integer, dimension(4)                                   :: DomainSize
        integer, dimension(4)                                   :: WindowPosition
        integer, dimension(4)                                   :: WindowFrame
        integer                                                 :: LimitsArrayFactor
        integer, dimension(6)                                   :: LimitsArray
        integer(HID_T)                                          :: DataType
        integer                                                 :: ILB_MPI
        integer                                                 :: IUB_MPI  
        integer                                                 :: JLB_MPI
        integer                                                 :: JUB_MPI  
        integer                                                 :: KLB_MPI
        integer                                                 :: KUB_MPI 
        integer                                                 :: ILB
        integer                                                 :: IUB  
        integer                                                 :: JLB
        integer                                                 :: JUB  
        integer                                                 :: KLB
        integer                                                 :: KUB 
               
        !------------------------------------------------------------------------
            
        DataType            = NULL_INT
        NewGroupName        = NULL_STR
        Name                = NULL_STR
        idx                 = NULL_INT
        nItems              = NULL_INT
        GroupType           = NULL_INT
        Dimensions          = NULL_INT
        ILB_MPI             = NULL_INT
        IUB_MPI             = NULL_INT  
        JLB_MPI             = NULL_INT
        JUB_MPI             = NULL_INT  
        KLB_MPI             = NULL_INT
        KUB_MPI             = NULL_INT 
        ILB                 = NULL_INT
        IUB                 = NULL_INT  
        JLB                 = NULL_INT
        JUB                 = NULL_INT  
        KLB                 = NULL_INT
        KUB                 = NULL_INT 
        LimitsArrayFactor   = 0
        LimitsArray         = 1

        !Get the number of members in the Group  
        call GetHDF5GroupNumberOfItems (HDF5ID = ObjHDF5_InArray(1),                         &
                                        GroupName = adjustl(trim(GroupName))//"/",   &
                                        nItems = nItems, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteConsolidatedHDF - ModuleDDC - ERR01'

do1 :   do idx = 1, nItems

            !Gets information about the group
            call GetHDF5ObjectInfo (HDF5ID = ObjHDF5_InArray(1), FatherGroupName = adjustl(trim(GroupName))//"/",           &
                                    GroupPosition = idx, GroupName = Name,GroupType = GroupType,                    &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteConsolidatedHDF - ModuleDDC - ERR02'

                          
if1 :       if (GroupType == H5G_DATASET_F.AND.adjustl(trim(GroupName))//"/"//adjustl(trim(Name)) .NE. &
                    "/Grid/Decomposition") then

                !Check if it is Coordenates Feature 
if12 :          if (adjustl(trim(Name)) .EQ. "ConnectionX" .OR. adjustl(trim(Name)) .EQ. "ConnectionY" .OR.         &
                    adjustl(trim(Name)) .EQ. "Latitude" .OR. adjustl(trim(Name)) .EQ. "Longitude") then
                    
                    LimitsArrayFactor = 1
                    
                else
                
                    LimitsArrayFactor = 0.    
                    
                endif if12
                
do2 :           do i=1,NumberOfFiles

                    Dimensions          = NULL_INT
                    ILB_MPI             = NULL_INT
                    IUB_MPI             = NULL_INT  
                    JLB_MPI             = NULL_INT
                    JUB_MPI             = NULL_INT  
                    KLB_MPI             = NULL_INT
                    KUB_MPI             = NULL_INT 
                    ILB                 = NULL_INT
                    IUB                 = NULL_INT  
                    JLB                 = NULL_INT
                    JUB                 = NULL_INT  
                    KLB                 = NULL_INT
                    KUB                 = NULL_INT
                    LimitsArray         = NULL_INT
                    
                    DomainSize      = hash_getDomainSize(hash_map_in, FileArrayIn(i))
                    WindowPosition  = hash_getWindowPosition(hash_map_in, FileArrayIn(i))
                    WindowFrame     = hash_getWindowFrame(hash_map_in, FileArrayIn(i))
                    
                    !Get Dataset Information
                    call GetHDF5GroupID(HDF5ID = ObjHDF5_InArray(i), FatherGroupName = adjustl(trim(GroupName))//"/",       &
                                        GroupPosition = idx, GroupName = Name,                                              &
                                        Units = Units, Rank = Rank, Dimensions = Dimensions,                                &
                                        STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteConsolidatedHDF - ModuleDDC - ERR03'
                    
                    !Get Dataset DataType
                    call GetHDF5DataTypeID( HDF5ID = ObjHDF5_InArray(i), FatherGroupName = adjustl(trim(GroupName))//"/",   &
                                            GroupPosition = idx, GroupName = Name,                                          &
                                            DataType = DataType, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteConsolidatedHDF - ModuleDDC - ERR04'  
                    
                    call ReadHDFDataSet(ObjHDF5_In  = ObjHDF5_InArray(i),   GroupName   = adjustl(trim(GroupName)),           &
                                        obj_name    = adjustl(trim(Name)),  Rank        = Rank,                               &
                                        dims        = Dimensions,           DataType    = DataType,                           & 
                                        MPIMatrixes = MPIMatrixes)
                                        
if13 :              if(i == 1) then

                        call AuxArrayAllocation(DomainSize  = DomainSize,   Dimensions          = Dimensions,           &
                                                AuxMatrixes = AuxMatrixes,  LimitsArrayFactor   = LimitsArrayFactor,    & 
                                                MPIMatrixes = MPIMatrixes)

                    endif if13
                    
                    call GetLimitsHDF(  WindowPosition      = WindowPosition,           &
                                        Dimensions          = Dimensions,               &
                                        Rank                = Rank,                     &
                                        LimitsArrayFactor   = LimitsArrayFactor,        &
                                        GroupName           = GroupName,                & 
                                        LimitsArray         = LimitsArray)
                                        
                    call GetIndexMPIHDF(WindowFrame = WindowFrame,                      &
                                        LimitsArray = LimitsArray,                      &
                                        GroupName   = GroupName,                        &
                                        ILB_MPI     = ILB_MPI, IUB_MPI = IUB_MPI,       &
                                        JLB_MPI     = JLB_MPI, JUB_MPI = JUB_MPI,       &
                                        KLB_MPI     = KLB_MPI, KUB_MPI = KUB_MPI)
!                                        
                    call GetIndexAUXHDF(WindowPosition      = WindowPosition,          &
                                        DomainSize          = DomainSize,               &
                                        LimitsArray         = LimitsArray,              &
                                        LimitsArrayFactor   = LimitsArrayFactor,        &
                                        GroupName           = GroupName,                &
                                        ILB                 = ILB , IUB  = IUB,         &
                                        JLB                 = JLB , JUB  = JUB,         &
                                        KLB                 = KLB , KUB  = KUB)
                    
                    call MergeArraysMPI(ILB_MPI = ILB_MPI, IUB_MPI = IUB_MPI,                       &
                                        JLB_MPI = JLB_MPI, JUB_MPI = JUB_MPI,                       &
                                        KLB_MPI = KLB_MPI, KUB_MPI = KUB_MPI,                       &
                                        ILB = ILB, IUB = IUB,                                       &
                                        JLB = JLB, JUB = JUB,                                       &
                                        KLB = KLB, KUB = KUB,                                       &
                                        AuxMatrixes = AuxMatrixes,                                  &
                                        MPIMatrixes = MPIMatrixes)
                
                    call MPIArrayDeallocate(MPIMatrixes= MPIMatrixes)
                
                enddo do2
               
                call WriteHDFDataSetLimits( HDF5ID      = ObjHDF5_Out,  Rank        = Rank,         &
                                            DomainSize  = DomainSize,   Dimensions  = Dimensions,   &
                                            LimitsArrayFactor = LimitsArrayFactor,  GroupName   = GroupName)
                
                call WriteHDFDataSet(   HDF5ID      = ObjHDF5_Out,  GroupName   = GroupName,        &
                                        Name        = Name,         Units       = Units,            &
                                        Rank        = Rank,         DataType    = DataType,         &
                                        AuxMatrixes = AuxMatrixes)
                
                call AuxArrayDeallocate( AuxMatrixes= AuxMatrixes)
                                        
            elseif (GroupType ==H5G_GROUP_F) then if1
             
if3 :           if(adjustl(trim(GroupName)) .EQ. "/") then
                    
                    NewGroupName = adjustl(trim(GroupName))//adjustl(trim(Name))
                    
                elseif (adjustl(trim(GroupName))//"/"//adjustl(trim(Name)) .NE.                                         &
                        "/Grid/Decomposition") then if3
                    
                    NewGroupName = adjustl(trim(GroupName))//"/"//adjustl(trim(Name))
                    
                endif if3
                
if2 :           if (adjustl(trim(GroupName))//"/"//adjustl(trim(Name)) .NE. &
                    "/Grid/Decomposition") then
                    
                    call WriteConsolidatedHDF(  IDOut           = IDOut,                            &
                                                ObjHDF5_Out     = ObjHDF5_Out,                      & 
                                                GroupName       = adjustl(trim(NewGroupName)),      &
                                                IdInArray       = IdInArray,                        &
                                                ObjHDF5_InArray = ObjHDF5_InArray,                  &
                                                FileArrayIn     = FileArrayIn,                      &
                                                NumberOfFiles   = NumberOfFiles,                    &
                                                hash_map_in     = hash_map_in,                      &
                                                AuxMatrixes     = AuxMatrixes,                      &
                                                MPIMatrixes     = MPIMatrixes)
                                  
                endif if2

            endif if1
             
        end do do1

        !------------------------------------------------------------------------

    end subroutine WriteConsolidatedHDF
        
        !--------------------------------------------------------------------------
    
    subroutine GetIndexAuxHDF(  WindowPosition, DomainSize, LimitsArray,        &
                                LimitsArrayFactor, GroupName,                   &
                                ILB, IUB,                                       &
                                JLB, JUB,                                       &
                                KLB, KUB)

        !Arguments-------------------------------------------------------------
        integer, dimension(4),intent(IN)                        :: WindowPosition
        integer, dimension(4),intent(IN)                        :: DomainSize 
        integer, dimension(6),intent(IN)                        :: LimitsArray
        integer,intent(IN)                                      :: LimitsArrayFactor
        character(*), intent(IN)                                :: GroupName
        integer,intent(out)                                     :: ILB
        integer,intent(out)                                     :: IUB  
        integer,intent(out)                                     :: JLB
        integer,intent(out)                                     :: JUB  
        integer,intent(out)                                     :: KLB
        integer,intent(out)                                     :: KUB 
          
        !Local-------------------------------------------------------------------

        
        !------------------------------------------------------------------------
        
if16 :  if (DomainSize(1) .GT. 1 .AND. adjustl(trim(GroupName)) .NE. "/Time") then

            ILB = WindowPosition(1)-DomainSize(1)+1
            IUB = WindowPosition(2)-DomainSize(1)+1+LimitsArrayFactor
            JLB = WindowPosition(3)-DomainSize(3)+1
            JUB = WindowPosition(4)-DomainSize(3)+1+LimitsArrayFactor
            KLB = LimitsArray(5)
            KUB = LimitsArray(6)
            
        else if16
        
            ILB = LimitsArray(1)
            IUB = LimitsArray(2)  
            JLB = LimitsArray(3)
            JUB = LimitsArray(4)  
            KLB = LimitsArray(5) 
            KUB = LimitsArray(6)
        
        endif if16
               
        !------------------------------------------------------------------------
    
    end subroutine GetIndexAuxHDF
        
        !--------------------------------------------------------------------------
    
    subroutine GetIndexMPIHDF(  WindowFrame, LimitsArray, GroupName,        &
                                ILB_MPI, IUB_MPI,                           &
                                JLB_MPI, JUB_MPI,                           &
                                KLB_MPI, KUB_MPI)

        !Arguments-------------------------------------------------------------
        integer, dimension(4),intent(IN)                        :: WindowFrame   
        integer, dimension(6),intent(IN)                        :: LimitsArray
        character(*), intent(IN)                                :: GroupName
        integer,intent(out)                                     :: ILB_MPI
        integer,intent(out)                                     :: IUB_MPI  
        integer,intent(out)                                     :: JLB_MPI
        integer,intent(out)                                     :: JUB_MPI  
        integer,intent(out)                                     :: KLB_MPI
        integer,intent(out)                                     :: KUB_MPI 
          
        !Local-------------------------------------------------------------------

        
        !------------------------------------------------------------------------
        
if15 :  if (adjustl(trim(GroupName)) .NE. "/Time") then 
       
            ILB_MPI = LimitsArray(1) - WindowFrame(1) + 1 
            IUB_MPI = LimitsArray(2) - WindowFrame(1) + 1
            JLB_MPI = LimitsArray(3) - WindowFrame(3) + 1
            JUB_MPI = LimitsArray(4) - WindowFrame(3) + 1
            KLB_MPI = LimitsArray(5)
            KUB_MPI = LimitsArray(6)
        
        else if15
        
            ILB_MPI = LimitsArray(1)  
            IUB_MPI = LimitsArray(2) 
            JLB_MPI = LimitsArray(3) 
            JUB_MPI = LimitsArray(4) 
            KLB_MPI = LimitsArray(5)
            KUB_MPI = LimitsArray(6)
        
        endif if15
          
        !------------------------------------------------------------------------
    
    end subroutine GetIndexMPIHDF
        
        !--------------------------------------------------------------------------
    
    subroutine GetLimitsHDF(WindowPosition, Dimensions, Rank,                     &
                            LimitsArrayFactor, GroupName, LimitsArray)

        !Arguments-------------------------------------------------------------
        integer, dimension(4),intent(IN)                        :: WindowPosition
        integer, dimension(7),intent(IN)                        :: Dimensions 
        integer,intent(IN)                                      :: Rank
        integer,intent(IN)                                      :: LimitsArrayFactor
        character(*), intent(IN)                                :: GroupName   
        integer, dimension(6),intent(out)                       :: LimitsArray
          
        !Local-------------------------------------------------------------------

        
        !------------------------------------------------------------------------
        
        LimitsArray = 1
        
if14 :  if (Rank .EQ. 1 .AND. adjustl(trim(GroupName)) .NE. "/Time") then
                    
            LimitsArray(1)= WindowPosition(1)
            LimitsArray(2)= WindowPosition(2)
            
        elseif (Rank .EQ. 1 .AND. adjustl(trim(GroupName)) .EQ. "/Time") then if14
        
            LimitsArray(1)= 1
            LimitsArray(2)= Dimensions(1)    
                
        elseif (Rank .EQ. 2) then if14

            LimitsArray(1)= WindowPosition(1)
            LimitsArray(2)= WindowPosition(2)+LimitsArrayFactor
            LimitsArray(3)= WindowPosition(3)
            LimitsArray(4)= WindowPosition(4)+LimitsArrayFactor
                        
        elseif (Rank .EQ. 3) then if14
    
            LimitsArray(1)= WindowPosition(1)
            LimitsArray(2)= WindowPosition(2)
            LimitsArray(3)= WindowPosition(3)
            LimitsArray(4)= WindowPosition(4)
            LimitsArray(5)= 1
            LimitsArray(6)= Dimensions(3)
    
        endif if14
          
        !------------------------------------------------------------------------
    
    end subroutine GetLimitsHDF
    
        !------------------------------------------------------------------------
    
    subroutine MergeArraysMPI(ILB_MPI, IUB_MPI, JLB_MPI, JUB_MPI, KLB_MPI, KUB_MPI,     &
                              ILB, IUB, JLB, JUB, KLB, KUB,                             &
                              AuxMatrixes, MPIMatrixes)

        !Arguments-------------------------------------------------------------
        type(T_AuxMatrixes), pointer                            :: AuxMatrixes
        type(T_MPIMatrixes), pointer                            :: MPIMatrixes
        integer,intent(IN)                                      :: ILB_MPI
        integer,intent(IN)                                      :: IUB_MPI  
        integer,intent(IN)                                      :: JLB_MPI
        integer,intent(IN)                                      :: JUB_MPI  
        integer,intent(IN)                                      :: KLB_MPI
        integer,intent(IN)                                      :: KUB_MPI 
        integer,intent(IN)                                      :: ILB
        integer,intent(IN)                                      :: IUB  
        integer,intent(IN)                                      :: JLB
        integer,intent(IN)                                      :: JUB  
        integer,intent(IN)                                      :: KLB
        integer,intent(IN)                                      :: KUB 
           
        !Local-------------------------------------------------------------------

        
        !------------------------------------------------------------------------
        
if1 :   if(associated(MPIMatrixes%MPIDataI4_1D) .EQ. .TRUE.)   then

            AuxMatrixes%AuxDataI4_1D(ILB:IUB)=MPIMatrixes%MPIDataI4_1D(ILB_MPI:IUB_MPI)
            
        elseif (associated(MPIMatrixes%MPIDataI4_2D) .EQ. .TRUE.) then if1
        
            AuxMatrixes%AuxDataI4_2D(ILB:IUB,JLB:JUB)=MPIMatrixes%MPIDataI4_2D(ILB_MPI:IUB_MPI,JLB_MPI:JUB_MPI)

        elseif (associated(MPIMatrixes%MPIDataI4_3D) .EQ. .TRUE.) then if1

            AuxMatrixes%AuxDataI4_3D(ILB:IUB,JLB:JUB,KLB:KUB)=MPIMatrixes%MPIDataI4_3D(ILB_MPI:IUB_MPI,JLB_MPI:JUB_MPI,KLB_MPI:KUB_MPI)
            
        end if if1

if2 :   if(associated(MPIMatrixes%MPIDataR4_1D) .EQ. .TRUE.) then 

            AuxMatrixes%AuxDataR4_1D(ILB:IUB)=MPIMatrixes%MPIDataR4_1D(ILB_MPI:IUB_MPI)
            
        elseif (associated(MPIMatrixes%MPIDataR4_2D) .EQ. .TRUE.) then if2

            AuxMatrixes%AuxDataR4_2D(ILB:IUB,JLB:JUB)=MPIMatrixes%MPIDataR4_2D(ILB_MPI:IUB_MPI,JLB_MPI:JUB_MPI)

        elseif (associated(MPIMatrixes%MPIDataR4_3D) .EQ. .TRUE.) then if2

            AuxMatrixes%AuxDataR4_3D(ILB:IUB,JLB:JUB,KLB:KUB)=MPIMatrixes%MPIDataR4_3D(ILB_MPI:IUB_MPI,JLB_MPI:JUB_MPI,KLB_MPI:KUB_MPI)
            
        end if if2  

if3 :   if(associated(MPIMatrixes%MPIDataR8_1D) .EQ. .TRUE.) then 

            AuxMatrixes%AuxDataR8_1D(ILB:IUB)=MPIMatrixes%MPIDataR8_1D(ILB_MPI:IUB_MPI)

        elseif (associated(MPIMatrixes%MPIDataR8_2D) .EQ. .TRUE.) then if3

            AuxMatrixes%AuxDataR8_2D(ILB:IUB,JLB:JUB)=MPIMatrixes%MPIDataR8_2D(ILB_MPI:IUB_MPI,JLB_MPI:JUB_MPI)

        elseif (associated(MPIMatrixes%MPIDataR8_3D) .EQ. .TRUE.) then if3

            AuxMatrixes%AuxDataR8_3D(ILB:IUB,JLB:JUB,KLB:KUB)=MPIMatrixes%MPIDataR8_3D(ILB_MPI:IUB_MPI,JLB_MPI:JUB_MPI,KLB_MPI:KUB_MPI)
            
        end if if3
          
        !------------------------------------------------------------------------
    
    end subroutine MergeArraysMPI  

        !--------------------------------------------------------------------------
    
    subroutine AuxArrayAllocation(  DomainSize, Dimensions, AuxMatrixes,        &
                                    LimitsArrayFactor, MPIMatrixes)

        !Arguments-------------------------------------------------------------
        integer, dimension(7)                                   :: Dimensions 
        integer, dimension(4)                                   :: DomainSize
        type(T_AuxMatrixes), pointer                            :: AuxMatrixes
        integer                                                 :: LimitsArrayFactor
        type(T_MPIMatrixes), pointer                            :: MPIMatrixes
           
        !Local-------------------------------------------------------------------

        
        !------------------------------------------------------------------------
        
if1 :   if(associated(MPIMatrixes%MPIDataI4_1D) .EQ. .TRUE.)   then

            allocate(AuxMatrixes%AuxDataI4_1D(  Dimensions(1) + LimitsArrayFactor))
            
        elseif (associated(MPIMatrixes%MPIDataI4_2D) .EQ. .TRUE.) then if1
        
            allocate(AuxMatrixes%AuxDataI4_2D(  DomainSize(2) + LimitsArrayFactor,          &
                                                DomainSize(4) + LimitsArrayFactor))

        elseif (associated(MPIMatrixes%MPIDataI4_3D) .EQ. .TRUE.) then if1

            allocate(AuxMatrixes%AuxDataI4_3D(  DomainSize(2) + LimitsArrayFactor,          &
                                                DomainSize(4) + LimitsArrayFactor,          &
                                                Dimensions(3)))
            
        end if if1

if2 :   if(associated(MPIMatrixes%MPIDataR4_1D) .EQ. .TRUE.) then 

            allocate(AuxMatrixes%AuxDataR4_1D(  Dimensions(1) + LimitsArrayFactor))
            
        elseif (associated(MPIMatrixes%MPIDataR4_2D) .EQ. .TRUE.) then if2

            allocate(AuxMatrixes%AuxDataR4_2D(  DomainSize(2) + LimitsArrayFactor,          &
                                                DomainSize(4) + LimitsArrayFactor))

        elseif (associated(MPIMatrixes%MPIDataR4_3D) .EQ. .TRUE.) then if2

            allocate(AuxMatrixes%AuxDataR4_3D(  DomainSize(2) + LimitsArrayFactor,          &
                                                DomainSize(4) + LimitsArrayFactor,          &
                                                Dimensions(3)))
            
        end if if2  

if3 :   if(associated(MPIMatrixes%MPIDataR8_1D) .EQ. .TRUE.) then 

            allocate(AuxMatrixes%AuxDataR8_1D(  Dimensions(1) + LimitsArrayFactor))

        elseif (associated(MPIMatrixes%MPIDataR8_2D) .EQ. .TRUE.) then if3

            allocate(AuxMatrixes%AuxDataR8_2D(  DomainSize(2) + LimitsArrayFactor,          &
                                                DomainSize(4) + LimitsArrayFactor))

        elseif (associated(MPIMatrixes%MPIDataR8_3D) .EQ. .TRUE.) then if3

            allocate(AuxMatrixes%AuxDataR8_3D(  DomainSize(2) + LimitsArrayFactor,          &
                                                DomainSize(4) + LimitsArrayFactor,          &
                                                Dimensions(3)))

        end if if3
          
        !------------------------------------------------------------------------
    
    end subroutine AuxArrayAllocation  

        !--------------------------------------------------------------------------

    recursive subroutine WriteHDFDataSetLimits( HDF5ID, Rank, DomainSize,               &
                                                Dimensions, LimitsArrayFactor,          &  
                                                GroupName)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                     :: HDF5ID
        integer, intent(IN)                                     :: Rank
        integer, dimension(4), intent(IN)                       :: DomainSize
        integer, dimension(7), intent(IN)                       :: Dimensions    
        integer, intent(IN)                                     :: LimitsArrayFactor
        CHARACTER(*), intent(IN)                                :: GroupName
        
        !Local-------------------------------------------------------------------
        integer                                                 :: STAT_CALL
        integer                                                 :: ILB
        integer                                                 :: IUB  
        integer                                                 :: JLB
        integer                                                 :: JUB  
        integer                                                 :: KLB
        integer                                                 :: KUB
        
        !------------------------------------------------------------------------
        
        STAT_CALL           = NULL_INT

if21 :  if (Rank .EQ. 1) then

if22 :      if (DomainSize(1) .GT. 1 .AND. adjustl(trim(GroupName)) .NE. "/Time") then
            
                ILB = DomainSize(1)-DomainSize(1)+1
                IUB = DomainSize(2)-DomainSize(1)+1
                JLB = 1
                JUB = 1  
                KLB = 1 
                KUB = 1
            
            elseif (adjustl(trim(GroupName)) .EQ. "/Time") then if22
                      
                ILB= 1
                IUB= Dimensions(1)
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
                KUB = Dimensions(3)
            
            else if24  
                
                ILB = DomainSize(1)
                IUB = DomainSize(2)
                JLB = DomainSize(3)
                JUB = DomainSize(4)
                KLB = 1
                KUB = Dimensions(3)
            
            endif if24
                
        endif if21
        
        call HDF5SetLimits (HDF5ID = HDF5ID,                                        & 
                            ILB = ILB,   IUB = IUB,                                 &
                            JLB = JLB,   JUB = JUB,                                 &
                            KLB = KLB,   KUB = KUB,                                 &
                            STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSetLimits - ModuleDDC - ERR01'
        
        !------------------------------------------------------------------------

    end subroutine WriteHDFDataSetLimits

        !--------------------------------------------------------------------------
    
    subroutine AuxArrayDeallocate(AuxMatrixes)

        !Arguments-------------------------------------------------------------
        type(T_AuxMatrixes), pointer                            :: AuxMatrixes
        
        !Local-------------------------------------------------------------------
        integer                                                 :: STAT_CALL
        
        !------------------------------------------------------------------------
        
if1 :   if(associated(AuxMatrixes%AuxDataI4_1D) .EQ. .TRUE.)   then

            deallocate(AuxMatrixes%AuxDataI4_1D)
            
        elseif (associated(AuxMatrixes%AuxDataI4_2D) .EQ. .TRUE.) then if1
        
            deallocate(AuxMatrixes%AuxDataI4_2D)

        elseif (associated(AuxMatrixes%AuxDataI4_3D) .EQ. .TRUE.) then if1

            deallocate(AuxMatrixes%AuxDataI4_3D)
            
        end if if1

if2 :   if(associated(AuxMatrixes%AuxDataR4_1D) .EQ. .TRUE.) then 

            deallocate(AuxMatrixes%AuxDataR4_1D)
            
        elseif (associated(AuxMatrixes%AuxDataR4_2D) .EQ. .TRUE.) then if2

            deallocate(AuxMatrixes%AuxDataR4_2D)

        elseif (associated(AuxMatrixes%AuxDataR4_3D) .EQ. .TRUE.) then if2

            deallocate(AuxMatrixes%AuxDataR4_3D)
            
        end if if2  

if3 :   if(associated(AuxMatrixes%AuxDataR8_1D) .EQ. .TRUE.) then 

            deallocate(AuxMatrixes%AuxDataR8_1D)

        elseif (associated(AuxMatrixes%AuxDataR8_2D) .EQ. .TRUE.) then if3

            deallocate(AuxMatrixes%AuxDataR8_2D)

        elseif (associated(AuxMatrixes%AuxDataR8_3D) .EQ. .TRUE.) then if3

            deallocate(AuxMatrixes%AuxDataR8_3D)

        end if if3
          
        !------------------------------------------------------------------------
    
    end subroutine AuxArrayDeallocate  

        !--------------------------------------------------------------------------
        
    subroutine MPIArrayDeallocate(MPIMatrixes)

        !Arguments-------------------------------------------------------------
        type(T_MPIMatrixes), pointer                            :: MPIMatrixes
        
        !Local-------------------------------------------------------------------
        integer                                                 :: STAT_CALL
        
        !------------------------------------------------------------------------
        
if1 :   if(associated(MPIMatrixes%MPIDataI4_1D) .EQ. .TRUE.)   then

            deallocate(MPIMatrixes%MPIDataI4_1D)
            
        elseif (associated(MPIMatrixes%MPIDataI4_2D) .EQ. .TRUE.) then if1
        
            deallocate(MPIMatrixes%MPIDataI4_2D)

        elseif (associated(MPIMatrixes%MPIDataI4_3D) .EQ. .TRUE.) then if1

            deallocate(MPIMatrixes%MPIDataI4_3D)
            
        end if if1

if2 :   if(associated(MPIMatrixes%MPIDataR4_1D) .EQ. .TRUE.) then 

            deallocate(MPIMatrixes%MPIDataR4_1D)
            
        elseif (associated(MPIMatrixes%MPIDataR4_2D) .EQ. .TRUE.) then if2

            deallocate(MPIMatrixes%MPIDataR4_2D)

        elseif (associated(MPIMatrixes%MPIDataR4_3D) .EQ. .TRUE.) then if2

            deallocate(MPIMatrixes%MPIDataR4_3D)
            
        end if if2  

if3 :   if(associated(MPIMatrixes%MPIDataR8_1D) .EQ. .TRUE.) then 

            deallocate(MPIMatrixes%MPIDataR8_1D)

        elseif (associated(MPIMatrixes%MPIDataR8_2D) .EQ. .TRUE.) then if3

            deallocate(MPIMatrixes%MPIDataR8_2D)

        elseif (associated(MPIMatrixes%MPIDataR8_3D) .EQ. .TRUE.) then if3

            deallocate(MPIMatrixes%MPIDataR8_3D)

        end if if3
          
        !------------------------------------------------------------------------
    
    end subroutine MPIArrayDeallocate  

        !--------------------------------------------------------------------------

    subroutine WriteHDFDataSet( HDF5ID, GroupName, Name, Units,            &
                                Rank, DataType, AuxMatrixes)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                     :: HDF5ID
        character(*), intent(IN)                                :: GroupName
        character(*), intent(IN)                                :: Name
        character(*), intent(IN)                                :: Units  
        integer, intent(IN)                                     :: Rank
        integer(HID_T), intent(IN)                              :: DataType
        type(T_AuxMatrixes), pointer                            :: AuxMatrixes
        
        !Local-------------------------------------------------------------------
        integer                                                 :: STAT_CALL    
        
        !------------------------------------------------------------------------
        
        STAT_CALL = NULL_INT
        
if1 :   if (Rank .EQ. 1) then

if12 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                call HDF5WriteData( HDF5ID  = HDF5ID,                   GroupName   = adjustl(trim(GroupName)),     & 
                                    Name    = adjustl(trim(Name)),      Units       = adjustl(trim(Units)),         &
                                    Array1D = AuxMatrixes%AuxDataI4_1D, STAT        = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR01'
                
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if12
            
                call HDF5WriteData( HDF5ID  = HDF5ID,                   GroupName   = adjustl(trim(GroupName)),     & 
                                    Name    = adjustl(trim(Name)),      Units       = adjustl(trim(Units)),         &
                                    Array1D = AuxMatrixes%AuxDataR4_1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR02'

            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if12
            
                call HDF5WriteData( HDF5ID  = HDF5ID,                   GroupName   = adjustl(trim(GroupName)),     & 
                                    Name    = adjustl(trim(Name)),      Units       = adjustl(trim(Units)),         &
                                    Array1D = AuxMatrixes%AuxDataR8_1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR03'

            endif if12

        elseif (Rank .EQ. 2) then if1
        
if13 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 
                          
                call HDF5WriteData( HDF5ID  = HDF5ID,                   GroupName   = adjustl(trim(GroupName)),     & 
                                    Name    = adjustl(trim(Name)),      Units       = adjustl(trim(Units)),         &
                                    Array2D = AuxMatrixes%AuxDataI4_2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR04'

            elseif (DataType .EQ. H5T_NATIVE_REAL) then if13
            
                call HDF5WriteData(HDF5ID  = HDF5ID,                   GroupName   = adjustl(trim(GroupName)),     & 
                                    Name    = adjustl(trim(Name)),      Units       = adjustl(trim(Units)),         &
                                    Array2D = AuxMatrixes%AuxDataR4_2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR05'

            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if13
            
                call HDF5WriteData( HDF5ID  = HDF5ID,                   GroupName   = adjustl(trim(GroupName)),     & 
                                    Name    = adjustl(trim(Name)),      Units       = adjustl(trim(Units)),         &
                                    Array2D = AuxMatrixes%AuxDataR8_2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR06'
            
            endif if13
        
        elseif (Rank .EQ. 3) then if1
        
if14 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 
                 
                call HDF5WriteData( HDF5ID  = HDF5ID,                   GroupName   = adjustl(trim(GroupName)),     & 
                                    Name    = adjustl(trim(Name)),      Units       = adjustl(trim(Units)),         &
                                    Array3D = AuxMatrixes%AuxDataI4_3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR07'
            
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if14
            
                call HDF5WriteData( HDF5ID  = HDF5ID,                   GroupName   = adjustl(trim(GroupName)),     & 
                                    Name    = adjustl(trim(Name)),      Units       = adjustl(trim(Units)),         &
                                    Array3D = AuxMatrixes%AuxDataR4_3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR08'

            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if14
            
                call HDF5WriteData( HDF5ID  = HDF5ID,                   GroupName   = adjustl(trim(GroupName)),     & 
                                    Name    = adjustl(trim(Name)),      Units       = adjustl(trim(Units)),         &
                                    Array3D = AuxMatrixes%AuxDataR8_3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR09'
            
            endif if14
        
        endif if1
        
        !Writes everything to disk
        call HDF5FlushMemory (HDF5ID, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDFDataSet - ModuleDDC - ERR10'
          
        !------------------------------------------------------------------------
    
    end subroutine WriteHDFDataSet

    !--------------------------------------------------------------------------

    subroutine ReadHDFDataSet(ObjHDF5_In, GroupName, obj_name,               &
                              Rank, dims, DataType, MPIMatrixes)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                     :: ObjHDF5_In
        character(*), intent(IN)                                :: GroupName
        character(*), intent(IN)                                :: obj_name 
        integer, intent(IN)                                     :: Rank
        integer, dimension(7), intent(IN)                       :: dims 
        integer(HID_T), intent(IN)                              :: DataType        
        type(T_MPIMatrixes), pointer                            :: MPIMatrixes
        
        !Local-------------------------------------------------------------------
        integer                                                 :: STAT_CALL   
        
        !------------------------------------------------------------------------

        STAT_CALL    = NULL_INT

if1 :   if (Rank .EQ. 1) then

            call HDF5SetLimits  (ObjHDF5_In, 1, dims(1), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDFDataSet - ModuleDDC - ERR01'

if12 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                allocate (MPIMatrixes%MPIDataI4_1D(1:dims(1)))
                
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array1D = MPIMatrixes%MPIDataI4_1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR02'
            
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if12
            
                allocate (MPIMatrixes%MPIDataR4_1D(1:dims(1)))
                
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array1D = MPIMatrixes%MPIDataR4_1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR03'
            
            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if12
            
                allocate (MPIMatrixes%MPIDataR8_1D(1:dims(1)))
              
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array1D = MPIMatrixes%MPIDataR8_1D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR04'
            
            endif if12

        elseif (Rank .EQ. 2) then if1
        
            call HDF5SetLimits  (ObjHDF5_In, 1, dims(1),1, dims(2), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDFDataSet - ModuleDDC - ERR05'
        
if13 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                allocate (MPIMatrixes%MPIDataI4_2D(1:dims(1),1:dims(2)))
                          
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)),  & 
                                  Array2D = MPIMatrixes%MPIDataI4_2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR06'
                
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if13
            
                allocate (MPIMatrixes%MPIDataR4_2D(1:dims(1),1:dims(2)))
                
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array2D = MPIMatrixes%MPIDataR4_2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR07'

            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if13
            
                allocate (MPIMatrixes%MPIDataR8_2D(1:dims(1),1:dims(2)))
                          
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array2D = MPIMatrixes%MPIDataR8_2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR08'

            endif if13
        
        elseif (Rank .EQ. 3) then if1
        
            call HDF5SetLimits  (ObjHDF5_In, 1, dims(1),1, dims(2),1,dims(3), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDFDataSet - ModuleDDC - ERR09'
        
if14 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then 

                allocate (MPIMatrixes%MPIDataI4_3D(1:dims(1),1:dims(2),1:dims(3)))
                
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array3D = MPIMatrixes%MPIDataI4_3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR10'
                
            elseif (DataType .EQ. H5T_NATIVE_REAL) then if14
            
                allocate (MPIMatrixes%MPIDataR4_3D(1:dims(1),1:dims(2),1:dims(3)))
                          
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array3D = MPIMatrixes%MPIDataR4_3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDFDataSet - ModuleDDC - ERR11'
                
            elseif (DataType .EQ. H5T_NATIVE_DOUBLE) then if14
            
                allocate (MPIMatrixes%MPIDataR8_3D(1:dims(1),1:dims(2),1:dims(3)))
                          
                call HDF5ReadData(ObjHDF5_In, GroupName, adjustl(trim(obj_name)), & 
                                  Array3D = MPIMatrixes%MPIDataR8_3D, STAT = STAT_CALL)
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
                    ConsolidatedFile = adjustl(trim(Directory)) // '/' // MPIResultsFile(5+iUnderScore2:)
                    
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

                    call hash_set(hash_map_in, key = adjustl(trim(Directory)) // '/' // HDFinfile, value_ = IDOut)
                
                    call hash_setObjID(hash_map_in, key = adjustl(trim(Directory)) // '/' // HDFinfile, ObjID = ObjHDF5_Out)
                
                    call GetDomainDecomposition(HDFFile = adjustl(trim(Directory)) // '/' // HDFinfile,   &
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
        allocate(DataVal1D(1:4))

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

            deallocate(DataVal1D)
    
            call GetDomainDecomposition(HDFFile = HDFFile,   &
                                            hash_map_in = hash_map_in,                                &
                                            FirstTime = FirstTime+1)
        endif if3
        
if4 :   if (associated (DataVal1D)) then 
            deallocate(DataVal1D) 
        endif if4
        
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
                         adjustl(trim(Directory))//'/MPI_'//                    &
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

