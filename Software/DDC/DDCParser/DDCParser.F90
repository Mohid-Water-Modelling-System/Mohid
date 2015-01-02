!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model - Results consolidation
! PROJECT       : DDC - Domain Decomposition Consolidation
! PROGRAM       : MainDDC
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : 2014
! REVISION      : Ricardo Miranda
! DESCRIPTION   : Program to consolidate results from MPI run with domain decomposition
!
!------------------------------------------------------------------------------
! Run:
! /opt/mpich/bin/mpiexec -n 1 /myDirectory/projects/lang/fortran/DDC/src/Solutions/Linux/DDCParser/DDCParser : -n 2  /myDirectory/projects/lang/fortran/DDC/src/Solutions/Linux/DCCWorker/DCCWorker

program DCCParser

    use mpi
    use moduleMPImanagement

    use ModuleGlobalData
    use ModuleHashTable

    implicit none

    !Constructor
!    public  ::      ConstructDDC
!    private ::          startMPI
!    private ::          readTreeFile
!    private ::              readLines
!    private ::                  modelLevel
!    private ::                  modelPath
!    private ::                  allocateDirectoryList
!    public  ::      createTasks
!    private ::          scanDirectoryList
!    private ::              scanFileList
!    private ::              openDecomposedFiles
!    private ::                  openDecomposedFiles2
!    private ::      barrier

    !Sets & Gets
!    private ::      setMyMPI_id
!    private ::      getMyMPI_id
!    private ::      setNumprocs
!    private ::      getNumprocs
!    private ::      setSlash
!    private ::      getSlash

    !Modifier
!    private :: main
!    public  ::      Loop

    !Recv
!    private ::      recvIdleWorker
!    private ::      recvTaskCompleted
!    private ::          isLastTask
!    private ::              endPrg ««-- Point of EXIT

    !Send
!    private ::          sendTask
!    private ::              switchTask
!    private ::          sendWokersPoisonPill
!    private ::          sendMyMPI_id

    !Destructor
!    public  ::      killDDC
!    private ::          deallocateInstance
!    private ::              deallocateDirectoryList
!    private ::          stopMPI

    !Types---------------------------------------------------------------------
    type       T_DirectoryList        !    private
        character(PathLength)                 :: Directory      = NULL_STR
        type(T_HashTable),            pointer :: hash_map_out   !consolidated HDF5 output files
        type(T_HashTable),            pointer :: hash_map_in    !decomposed HDF5 input files
        type(T_DirectoryList),        pointer :: Next
    end type  T_DirectoryList

    type       T_DDC                  !    public
        integer                               :: myMPI_id       = NULL_INT
        integer                               :: numprocs       = NULL_INT
        integer                               :: nbrModels      = NULL_INT
        type(T_DirectoryList),        pointer :: DirectoryList
        type(T_HashTable),            pointer :: hash_tasks     !list of tasks (HDF files to consolidate)

        character(1)                          :: slash          = '/'
    end type  T_DDC

    !--------------------------------------------------------------------------

    call main()

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN M

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine main()
        type (T_DDC), pointer       :: Me
        integer                     :: STAT_CALL  = UNKNOWN_
        character(len=1)            :: slash      = '*'

        nullify(Me)
        Me => constructDDC()
        if (.NOT. associated(Me))                                             &
            stop "subroutine main, program DDCParser, error calling ConstructDDC, ERR01"

if7 :   if (command_argument_count() .GT. 0) then
            call get_command_argument(1,                                      &
                                      slash,                                  &
                                      STATUS = STAT_CALL)
if2 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "subroutine main, program DDCParser, error calling setSlash, ERR07"
            end if if2

            STAT_CALL = setSlash(Me, slash)
if3 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "subroutine main, program DDCParser, error calling setSlash, ERR03"
            end if if3
        end if if7

        STAT_CALL = createTasks(Me)
if1 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine main, program DDCParser, error calling createTasks, ERR02"
        end if if1

        STAT_CALL = barrier()
if5 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine main, program DDCParser, error calling barrier, ERR05"
        end if if5

        STAT_CALL = sendMyMPI_id(Me)
if6 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine main, program DDCParser, error calling sendMyMPI_id, ERR06"
        end if if6

        call Loop(Me)

    end subroutine main

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    function constructDDC()
        type (T_DDC), pointer       :: constructDDC
        type (T_DDC), pointer       :: NewObjDDC
        integer                     :: STAT_CALL  = UNKNOWN_

        allocate(NewObjDDC)
        NewObjDDC%nbrModels  =  0         !Initializes model count
        NewObjDDC%hash_tasks => hash_init()

        nullify (NewObjDDC%DirectoryList)

        STAT_CALL = startMPI    (NewObjDDC)
if1 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function constructDDC, program DDCParser, error calling startMPI, ERR01"
        end if if1

        STAT_CALL = readTreeFile(NewObjDDC)
if2 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function constructDDC, program DDCParser, error calling readTreeFile, ERR02"
        end if if2

        constructDDC => NewObjDDC

    end function constructDDC

    !---------------------------------------------------------------------------

    integer function startMPI(Me)
        type (T_DDC), pointer       :: Me
        integer                     :: STAT_CALL  = UNKNOWN_
        integer                     :: myMPI_id   = NULL_INT
        integer                     :: numprocs   = NULL_INT

        call MPI_INIT(IERROR = STAT_CALL)
if1 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function startMPI, program DDCParser, error calling MPI_INIT, ERR01"
        end if if1

        call MPI_COMM_RANK(MPI_COMM_WORLD,                             &
                           myMPI_id,                                   &
                           IERROR = STAT_CALL)
if2 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function startMPI, program DDCParser, error calling MPI_COMM_RANK, ERR02.1"
        end if if2

        STAT_CALL = setMyMPI_id(Me, myMPI_id = myMPI_id)
if21 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function startMPI, program DDCParser, error calling setMyMPI_id, ERR02.2"
        end if if21

        call MPI_COMM_SIZE(MPI_COMM_WORLD,                               &
                           numprocs,                                     &
                           IERROR = STAT_CALL)
if3 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function startMPI, program DDCParser, error calling MPI_COMM_SIZE, ERR03"
        end if if3

        STAT_CALL = setNumprocs(Me, numprocs)
if5 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function startMPI, program DDCParser, error calling MPI_COMM_SIZE, ERR03"
        end if if5

        print *, 'Process ', myMPI_id, ' of ', numprocs, ' is alive -> Program DDCParser'

        startMPI = SUCCESS_

    end function startMPI

    !---------------------------------------------------------------------------

    integer function barrier()
        integer                     :: STAT_CALL  = UNKNOWN_

        call MPI_BARRIER(MPI_COMM_WORLD, STAT_CALL)
if4 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function startMPI, program DDCParser, error calling MPI_BARRIER, ERR04"
        end if if4

        barrier = SUCCESS_

    end function barrier

    !--------------------------------------------------------------------------

    integer function readTreeFile(Me)
        type (T_DDC), pointer       :: Me

        logical                     :: TreeExists           = .false.
        logical                     :: TreeExistsCapitall   = .false. 
        integer                     :: STAT_                = UNKNOWN_
        integer                     :: STAT_CALL            = UNKNOWN_
        integer                     :: iTree                = NULL_INT        
        character(StringLength)     :: Coment1              = NULL_STR
        character(StringLength)     :: Coment2              = NULL_STR

        STAT_ = SUCCESS_        
        
        call UnitsManager(iTree, OPEN_FILE, STAT = STAT_CALL)          
        
if1 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function readTreeFile, program DDCParser, error calling UnitsManager, ERR01"
        end if if1

        !Verifies if Tree file exists and allocates the list of models
        inquire(file='tree.dat', EXIST = TreeExists)

if2:    if (TreeExists) then

            open(UNIT = iTree, FILE = 'tree.dat', status = 'OLD', IOSTAT = STAT_CALL)
if3 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function readTreeFile, program DDCParser, error calling UnitsManager, ERR02"
            end if if3
            
        else if2          

            inquire(file='Tree.dat', EXIST = TreeExistsCapitall)
            
if4:        if (TreeExistsCapitall) then

                open(UNIT = iTree, FILE = 'Tree.dat', status = 'OLD', IOSTAT = STAT_CALL)
                
if5 :           if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "function readTreeFile, program DDCParser, error calling open, ERR03"
                end if if5
                
            else  if4
                            
                STAT_ = FILE_NOT_FOUND_ERR_
                
            endif if4
            
        endif if2            
        
        rename

        read(unit=iTree, fmt=*) Coment1
        read(unit=iTree, fmt=*) Coment2

        call readLines(Me,                                                  &
                       iTree         = iTree,                               &
                       DirectoryList = Me%DirectoryList)

        call UnitsManager(iTree, CLOSE_FILE, STAT = STAT_CALL)
if6 :       if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function readTreeFile, program DDCParser, error calling UnitsManager, ERR03"
        end if if6

        STAT_ = STAT_CALL

        readTreeFile = STAT_

    end function readTreeFile

    !--------------------------------------------------------------------------

    recursive subroutine readLines(Me, iTree, DirectoryList)
        type (T_DDC), pointer           :: Me
        type(T_DirectoryList), pointer  :: DirectoryList
        integer, intent(IN)             :: iTree

        !Local-----------------------------------------------------------------
        integer                         :: iMPI         = NULL_INT
        integer                         :: STAT_CALL    = UNKNOWN_
        integer                         :: ModelLevel_  = NULL_INT
        character(PathLength)           :: AuxString    = NULL_STR
        character(PathLength)           :: AuxString2   = NULL_STR
        character(PathLength)           :: ModelPath_   = NULL_STR
        integer                         :: nbrModels_   = NULL_INT

        read(unit = iTree, fmt='(a256)', IOSTAT = STAT_CALL) AuxString

if1 :   if (STAT_CALL .EQ. SUCCESS_) then
            AuxString2 = trim(adjustl(AuxString))
if2 :       if (AuxString2(1:1) == '+') then
                iMPI = scan(AuxString2,":")
if3 :           if (iMPI > 0) then
                    read(AuxString2(iMPI+1:),'(I)') nbrModels_
                    Me%nbrModels = Me%nbrModels + nbrModels_

                    ModelLevel_  = modelLevel(AuxString  = AuxString2,          &
                                              Level      = 0)   !Level=0 because it is the '+' counting start
                    ModelPath_   = modelPath (AuxString  = AuxString2,          &
                                              Level      = ModelLevel_)

                    DirectoryList => allocateDirectoryList(modelPath = ModelPath_)

                    call readLines(Me,                                          &
                                   iTree         = iTree,                       &
                                   DirectoryList = DirectoryList%Next)

                else if3
                    !This line does not have domain decomposition, reads next Tree.dat's line
                    Me%nbrModels = Me%nbrModels + 1

                    call readLines(Me,                                          &
                                   iTree         = iTree,                       &
                                   DirectoryList = DirectoryList)
                endif if3
            endif if2
        endif if1
    end subroutine readLines

    !--------------------------------------------------------------------------

    integer pure recursive function modelLevel(AuxString, Level)
        Character(len=*), intent(in)  :: AuxString
        integer, intent(in)           :: Level

if1 :   if (AuxString((Level+1):(Level+1)) == '+') then
            modelLevel = ModelLevel(AuxString = AuxString,                      &
                                    Level     = Level + 1)
        else if1
            modelLevel = Level
        endif if1

    end function modelLevel

    !--------------------------------------------------------------------------

    character(len=PathLength) pure function modelPath(AuxString, Level)
        character(len=*), intent(in)  :: AuxString
        integer,          intent(in)  :: Level
        integer                       :: position_

        position_  = scan(AuxString, ":")
        modelPath  = AuxString(Level+1:position_-5)//"res"

    end function modelPath

    !--------------------------------------------------------------------------

    function allocateDirectoryList(modelPath)
        type(T_DirectoryList), pointer    :: AllocateDirectoryList
        character(PathLength), intent(IN) :: modelPath
        type (T_DirectoryList), pointer   :: NewDirectoryList

        allocate(NewDirectoryList)
        NewDirectoryList%Directory = adjustl(trim(modelPath))
        NewDirectoryList%hash_map_out => hash_init()
        NewDirectoryList%hash_map_in  => hash_init()
        nullify(NewDirectoryList%Next)

        allocateDirectoryList => NewDirectoryList

    end function allocateDirectoryList

    !--------------------------------------------------------------------------

    integer function createTasks(Me)
        type(T_DDC),       pointer  :: Me

if1 :   if (associated(Me%DirectoryList)) then
            print*, 'Using slash -> ', getSlash(Me)

            call scanDirectoryList(DirectoryList = Me%DirectoryList,          &
                                   nbrModels     = Me%nbrModels,              &
                                   hash_tasks    = Me%hash_tasks,             &
                                   slash         = getSlash(Me))
        endif if1

        createTasks = SUCCESS_

        !------------------------------------------------------------------------

    end function createTasks

    !--------------------------------------------------------------------------

    recursive subroutine scanDirectoryList(DirectoryList,                     &
                                           nbrModels,                         &
                                           hash_tasks,                        &
                                           slash)

        type(T_DirectoryList), pointer  :: DirectoryList
        integer, intent(IN)             :: nbrModels
        character(PathLength)           :: FirstHDFFileOut  = NULL_STR
        type (T_HashTable),    pointer  :: hash_tasks
        character(1), intent(IN)        :: slash

        call openDecomposedFiles(hash_map_out = DirectoryList%hash_map_out,   &
                                 hash_map_in  = DirectoryList%hash_map_in,    &
                                 Directory    = DirectoryList%Directory,      &
                                 nbrModels    = nbrModels,                    &
                                 slash        = slash)

if1 :   if (hash_get_first_exists(DirectoryList%hash_map_in)) then
            FirstHDFFileOut = hash_get_first_key(DirectoryList%hash_map_out)

            call scanFileList(hash_map_in  = DirectoryList%hash_map_in,       &
                              HDFFileOut   = FirstHDFFileOut,                 &
                              hash_map_out = DirectoryList%hash_map_out,      &
                              hash_tasks   = hash_tasks)
        endif if1

        if (associated(DirectoryList%Next))                                   &
            call scanDirectoryList(DirectoryList = DirectoryList%Next,        &
                                   nbrModels     = nbrModels,                 &
                                   hash_tasks    = hash_tasks,                &
                                   slash         = slash)

    end subroutine scanDirectoryList

    !--------------------------------------------------------------------------

    subroutine openDecomposedFiles(hash_map_out,                              &
                                   hash_map_in,                               &
                                   Directory,                                 &
                                   nbrModels,                                 &
                                   slash)

        type(T_HashTable), pointer        :: hash_map_out
        type(T_HashTable), pointer        :: hash_map_in
        character(PathLength), intent(IN) :: Directory
        integer, intent(IN)               :: nbrModels
        character(1), intent(IN)          :: slash
        character(StringLength)           :: DecomposedFiles  = NULL_STR
        DecomposedFiles = adjustl(trim('DecomposedFiles'))

        !Starts looking for first DecomposedFiles: DecomposedFiles_0.dat
        call openDecomposedFiles2(hash_map_out    = hash_map_out,               &
                                  hash_map_in     = hash_map_in,                &
                                  Directory       = Directory,                  &
                                  DecomposedFiles = DecomposedFiles,            &
                                  iModel          = 0,                          &
                                  nbrModels       = nbrModels,                  &
                                  slash           = slash)

    end subroutine openDecomposedFiles

    !--------------------------------------------------------------------------

    recursive subroutine openDecomposedFiles2(hash_map_out,                     &
                                              hash_map_in,                      &
                                              Directory,                        &
                                              DecomposedFiles,                  &
                                              iModel,                           &
                                              nbrModels,                        &
                                              slash)

        type(T_HashTable), pointer          :: hash_map_out
        type(T_HashTable), pointer          :: hash_map_in
        character(StringLength), intent(IN) :: DecomposedFiles
        character(PathLength), intent(IN)   :: Directory
        integer, intent(IN)                 :: iModel
        integer, intent(IN)                 :: nbrModels
        character(1), intent(IN)            :: slash

        logical                             :: DDFileExists
        integer                             :: STAT_CALL        = NULL_INT
        integer                             :: iDDFile          = NULL_INT
        character(StringLength)             :: DDFile           = NULL_STR
        character(StringLength)             :: AuxString        = NULL_STR

if2 :   if (iModel .LT. nbrModels) then
            write (AuxString, '(i10)') iModel
            DDFile = adjustl(trim(                                            &
                         adjustl(trim(Directory))//'/MPI_'//                  &
                         adjustl(trim(AuxString))))//'_'//                    &
                         adjustl(trim(DecomposedFiles))//'.dat'
            inquire(file = DDFile, EXIST = DDFileExists)

if1 :       if (DDFileExists) then
                call UnitsManager(iDDFile, OPEN_FILE, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                    &
                    stop "function openDecomposedFiles2, program DDCParser, error calling UnitsManager, ERR01"

                open(UNIT = iDDFile, FILE = DDFile, status = 'OLD', IOSTAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                    &
                    stop "function openDecomposedFiles2, program DDCParser, error calling UnitsManager, ERR02"

                call scanDecomposedFiles(iHDFFile      = iDDFile,             &
                                         hash_map_out  = hash_map_out,        &
                                         hash_map_in   = hash_map_in,         &
                                         Directory     = Directory,           &
                                         slash         = slash)

                call UnitsManager(iDDFile, CLOSE_FILE, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                    &
                    stop "function openDecomposedFiles2, program DDCParser, error calling UnitsManager, ERR03"
            endif if1

            call openDecomposedFiles2(hash_map_out    = hash_map_out,         &
                                      hash_map_in     = hash_map_in,          &
                                      Directory       = Directory,            &
                                      DecomposedFiles = DecomposedFiles,      &
                                      iModel          = iModel + 1,           &
                                      nbrModels       = nbrModels,            &
                                      slash           = slash)
        endif if2

    end subroutine openDecomposedFiles2

    !--------------------------------------------------------------------------

    recursive subroutine scanDecomposedFiles(iHDFFile,                        &
                                             hash_map_out,                    &
                                             hash_map_in,                     &
                                             Directory,                       &
                                             slash)

        character(PathLength)       :: Directory
        integer, intent(IN)         :: iHDFFile
        type(T_HashTable), pointer  :: hash_map_out
        type(T_HashTable), pointer  :: hash_map_in
        character(1), intent(IN)    :: slash

        integer                     :: iUnderScore
        integer                     :: iUnderScore2
        integer                     :: STAT_CALL
        character(PathLength)       :: HDFinfile
        character(PathLength)       :: MPIResultsFile
        character(PathLength)       :: ConsolidatedFile

        iUnderScore        =  NULL_INT
        iUnderScore2       =  NULL_INT
        STAT_CALL          =  UNKNOWN_
        HDFinfile          =  NULL_STR
        MPIResultsFile     =  NULL_STR
        ConsolidatedFile   =  NULL_STR

        read(unit = iHDFFile, fmt='(a256)', IOSTAT = STAT_CALL) HDFinfile
if1 :   if (STAT_CALL .EQ. SUCCESS_) then
            !Checks if it is a Domain Decomposition HDF5 results file type
            iUnderScore = scan(HDFinfile,"MPI_")
if2 :       if (iUnderScore > 0) then
                iUnderScore2 = scan(HDFinfile(5:),"_")

if3 :           if (iUnderScore2 > 0) then
                    ConsolidatedFile = trim(adjustl(Directory)) // slash //   &
                                       trim(adjustl(HDFinfile(5+iUnderScore2:)))

                    if (hash_get(hash_map_out, ConsolidatedFile) < 0)         &
                        call hash_set(hash_map_out, key = ConsolidatedFile)

                    call hash_set(hash_map_in,                                  &
                                  key    = trim(adjustl(Directory)) // slash // &
                                           trim(adjustl(HDFinfile)),            &
                                  value_ = hash_get(hash_map_out, key = ConsolidatedFile))
                endif if3
            endif if2

            call scanDecomposedFiles(iHDFFile     = iHDFFile,                 &
                                     hash_map_out = hash_map_out,             &
                                     hash_map_in  = hash_map_in,              &
                                     Directory    = Directory,                &
                                     slash        = slash)
        endif if1

    end subroutine scanDecomposedFiles

    !--------------------------------------------------------------------------

    recursive subroutine scanFileList(hash_map_in,                            &
                                      HDFFileOut,                             &
                                      hash_map_out,                           &
                                      hash_tasks)

        CHARACTER(PathLength), intent(IN) :: HDFFileOut
        type(T_HashTable), pointer        :: hash_map_in
        type(T_HashTable), pointer        :: hash_map_out
        type(T_HashTable), pointer        :: hash_tasks

        CHARACTER(PathLength)             :: HDFFileIn
        CHARACTER(PathLength)             :: TaskFile
        integer                           :: iTaskFile
        integer                           :: nTaskFile
        integer                           :: IDOut
        integer                           :: IDIn
        integer                           :: STAT_CALL
        logical                           :: continue_
        CHARACTER(PathLength)             :: auxStr

        integer                           :: I
        integer,dimension(8)              :: values       !DATE_AND_TIME

        HDFFileIn    = NULL_STR
        iTaskFile    = NULL_INT
        nTaskFile    = 0
        IDOut        = NULL_INT
        IDIn         = NULL_INT
        continue_    =.TRUE.

        HDFFileIn = hash_get_first_key(hash_map_in)
        IDOut     = hash_get          (hash_map_out, HDFFileOut)

        call sleep(2)     !pauses to enforce different file name

        call DATE_AND_TIME(VALUES=values)
do2 :   do I = 1, 8
            nTaskFile = nTaskFile * 100 +  values(I)
        end do do2

        nTaskFile = abs(nTaskFile)
        write (auxStr, '(i10)') nTaskFile
        TaskFile = trim(adjustl(auxStr)) // '.tmp'

        !puts a new task in the task list
        call hash_set(hash_tasks,                                             &
                      key    = TaskFile,                                      &
                      value_ = nTaskFile)

        call UnitsManager(iTaskFile, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'scanFileList - ModuleDDC - ERR01'

        open(UNIT = iTaskFile, FILE = TaskFile, status = 'NEW', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'scanFileList - ModuleDDC - ERR02'

        !The first line is the consolidated HDF file name
        write(unit=iTaskFile, fmt=*) trim(adjustl(HDFFileOut))

do1 :   do while(continue_)
            IDIn  = hash_get(hash_map_in,  HDFFileIn )

            if (IDOut .EQ. IDIn)                                              &
                write(unit=iTaskFile, fmt=*) trim(adjustl(HDFFileIn))

            continue_ = hash_get_next_exists(hash_map_in, HDFFileIn)
            if (continue_)                                                    &
                HDFFileIn = hash_get_next_key(hash_map_in, HDFFileIn)
        enddo do1

        if (hash_get_next_exists(hash_map_out, HDFFileOut))                               &
            call scanFileList(hash_map_in  = hash_map_in,                                 &
                              HDFFileOut   = hash_get_next_key(hash_map_out, HDFFileOut), &
                              hash_map_out = hash_map_out,                                &
                              hash_tasks   = hash_tasks)

        call UnitsManager(iTaskFile, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'scanFileList - ModuleDDC - ERR03'

    end subroutine scanFileList

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SETS & GETS SETS & GETS SETS & GETS SETS & GETS SETS & GETS SETS & GETS SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    integer function setMyMPI_id(Me, myMPI_id)
        type (T_DDC), pointer       :: Me
        integer, intent(IN)         :: myMPI_id

        Me%myMPI_id = myMPI_id

        setMyMPI_id = SUCCESS_

    end function setMyMPI_id

    !---------------------------------------------------------------------------

    integer function getMyMPI_id(Me)
        type (T_DDC), pointer       :: Me

        getMyMPI_id = Me%myMPI_id

    end function getMyMPI_id

    !---------------------------------------------------------------------------

    integer function setNumprocs(Me, numprocs)
        type (T_DDC), pointer       :: Me
        integer, intent(IN)         :: numprocs

        Me%numprocs = numprocs

        setNumprocs = SUCCESS_

    end function setNumprocs

    !---------------------------------------------------------------------------

    integer function getNumprocs(Me)
        type (T_DDC), pointer       :: Me

        getNumprocs = Me%numprocs

    end function getNumprocs

    !---------------------------------------------------------------------------

    integer function setSlash(Me, slash)
        type (T_DDC), pointer         :: Me
        character(len=1), intent(IN)  :: slash

        Me%slash = slash

        setSlash = SUCCESS_

    end function setSlash

    !---------------------------------------------------------------------------

    character(1) function getSlash(Me)
        type (T_DDC), pointer       :: Me

        getSlash = Me%slash

    end function getSlash

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    recursive subroutine Loop(Me)
        type (T_DDC), pointer       :: Me
        integer                     :: STAT_CALL   = UNKNOWN_
        integer                     :: STATUS(MPI_STATUS_SIZE)

        call MPI_PROBE(MPI_ANY_SOURCE,                                          &
                       MPI_ANY_TAG,                                             &
                       MPI_COMM_WORLD,                                          &
                       STATUS,                                                  &
                       STAT_CALL)
if25 :  if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine Loop, program DCCParser, error calling MPI_PROBE, ERR01"
        end if if25

if42 :  if (STATUS(MPI_TAG) .EQ. getMsgIdleWorkerTag()) then

            STAT_CALL = recvIdleWorker(Me%hash_tasks, workerRank = STATUS(MPI_SOURCE))
if4 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "subroutine Loop, program DCCParser, error calling recvIdleWorker, ERR02"
            end if if4

        else if (STATUS(MPI_TAG) .EQ. getMsgTaskCompletedTag()) then if42

            STAT_CALL = recvTaskCompleted(Me, workerRank = STATUS(MPI_SOURCE))
if2 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "subroutine Loop, program DCCParser, error calling recvTaskCompleted, ERR02"
            end if if2

        else if42

            print*, "received message tag = ", STATUS(MPI_TAG)
            stop "subroutine Loop, program DCCParser, message not recognized, ERR99"
        end if if42

        call Loop(Me)
    end subroutine Loop

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SEND SEND SEND SEND SEND SEND SEND SEND SEND SEND SEND SEND SEND SEND SEND S

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    integer function sendTask(hash_tasks, workerRank)
        type(T_HashTable), pointer  :: hash_tasks
        integer, intent(IN)         :: workerRank
        integer                     :: nTaskFile    = NULL_INT
        integer                     :: STAT_        = UNKNOWN_
        integer                     :: STAT_CALL    = NULL_INT

        nTaskFile = hash_get_first(hash_tasks)  !Picks a task to send to a worker
if1 :   if (nTaskFile .GE. 0) then
            call MPI_SEND(nTaskFile,                                          &
                          1,                                                  &
                          MPI_INTEGER,                                        &
                          workerRank,                                         &
                          getMsgSendTaskTag(),                                &
                          MPI_COMM_WORLD,                                     &
                          STAT_CALL)
if16 :      if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function sendTask, error calling MPI_SEND, ERR01"
            end if if16

            STAT_CALL = switchTask(hash_tasks)
if17 :      if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function sendTask, error calling MPI_SEND, ERR02"
            end if if17

            STAT_ = SUCCESS_

        else if (nTaskFile .EQ. 0) then if1

            !No more tasks to distribute

        else if1
            STAT_ = OUT_OF_BOUNDS_ERR_
        end if if1

        sendTask = STAT_

    end function sendTask

    !---------------------------------------------------------------------------

    integer function switchTask(hash_tasks)
        !The 1st element on the list is moved to the last position
        type(T_HashTable), pointer  :: hash_tasks
        CHARACTER(PathLength)       ::  TaskFile  = NULL_STR
        integer                     :: nTaskFile  = NULL_INT
        integer                     :: STAT_CALL  = UNKNOWN_

         TaskFile = hash_get_first_key(hash_tasks)
        nTaskFile = hash_get_first    (hash_tasks)

        STAT_CALL = hash_pop(hash_tasks, key = TaskFile)
if16 :  if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function switchTask, error calling hash_pop, ERR01"
        end if if16

        call hash_set(hash_tasks,                                             &
                      key    =  TaskFile,                                     &
                      value_ = nTaskFile)

        switchTask = SUCCESS_

    end function switchTask

    !---------------------------------------------------------------------------

    integer function sendWorkersPoisonPill(Me)
        type (T_DDC), pointer       :: Me
        integer                     :: STAT_CALL    = NULL_INT
        integer                     :: N            = NULL_INT

do1 :   do N = 0, getNumprocs(Me)-1
if2 :   if (N .NE. getMyMPI_id(Me)) then

            call MPI_SEND(getMsgPoisonPillTag(),                              &
                          1,                                                  &
                          MPI_INTEGER,                                        &
                          N,                                                  &
                          getMsgPoisonPillTag(),                              &
                          MPI_COMM_WORLD,                                     &
                          STAT_CALL)
if1 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function sendWorkersPoisonPill, error calling MPI_SEND, ERR01"
            end if if1
        end if if2
        end do do1

        sendWorkersPoisonPill = SUCCESS_

    end function sendWorkersPoisonPill


    !---------------------------------------------------------------------------

    integer function sendMyMPI_id(Me)
        type (T_DDC), pointer       :: Me
        integer                     :: STAT_CALL    = NULL_INT
        integer                     :: N            = NULL_INT

do1 :   do N = 0, getNumprocs(Me)-1
if2 :   if (N .NE. getMyMPI_id(Me)) then

            call MPI_SEND(getMyMPI_id(Me),                                    &
                          1,                                                  &
                          MPI_INTEGER,                                        &
                          N,                                                  &
                          getMsgRankTag(),                                    &
                          MPI_COMM_WORLD,                                     &
                          STAT_CALL)
if1 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function sendMyMPI_id, error calling MPI_SEND, ERR01"
            end if if1
        end if if2
        end do do1

        sendMyMPI_id = SUCCESS_

    end function sendMyMPI_id

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !RECV RECV RECV RECV RECV RECV RECV RECV RECV RECV RECV RECV RECV RECV RECV R

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    integer function recvIdleWorker(hash_tasks, workerRank)
        type(T_HashTable), pointer  :: hash_tasks
        integer,intent(IN)          :: workerRank
        integer                     :: STAT_CALL    = NULL_INT
        integer                     :: msg          = NULL_INT
        integer                     :: STATUS(MPI_STATUS_SIZE)

        call MPI_RECV(msg,                                                      &
                      1,                                                        &
                      MPI_INTEGER,                                              &
                      workerRank,                                               &
                      getMsgIdleWorkerTag(),                                    &
                      MPI_COMM_WORLD,                                           &
                      STATUS,                                                   &
                      STAT_CALL)
if1 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function recvIdleWorker, error calling MPI_RECV, ERR01"
        end if if1

        STAT_CALL = sendTask(hash_tasks, workerRank)
if2 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function recvIdleWorker, error calling sendTask, ERR02"
        end if if2

        recvIdleWorker = SUCCESS_

    end function recvIdleWorker

    !------------------------------------------------------------------------

    integer function recvTaskCompleted(Me, workerRank)
        type (T_DDC), pointer       :: Me
        integer,intent(IN)          :: workerRank
        integer                     :: STAT_CALL    = NULL_INT
        integer                     :: nTaskFile    = NULL_INT
        integer                     :: STATUS(MPI_STATUS_SIZE)
        CHARACTER(PathLength)       :: TaskFile     = NULL_STR
        integer                     :: iFile        = NULL_INT

        call MPI_RECV(nTaskFile,                                                &
                      1,                                                        &
                      MPI_INTEGER,                                              &
                      workerRank,                                               &
                      getMsgTaskCompletedTag(),                                 &
                      MPI_COMM_WORLD,                                           &
                      STATUS,                                                   &
                      STAT_CALL)
if1 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function recvTaskCompleted, error calling MPI_SEND, ERR01"
        end if if1

        write (TaskFile, '(i10)') nTaskFile
        STAT_CALL = hash_pop(Me%hash_tasks, trim(adjustl(TaskFile)) // '.tmp')
if2 :   if ((STAT_CALL .NE.       SUCCESS_) .AND.                                   &
            (STAT_CALL .NE. (-1 * NOT_FOUND_ERR_))) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function recvTaskCompleted, error calling isLastTask, ERR02"
        end if if2

        !Delete task file
        call UnitsManager(iFile, OPEN_FILE, STAT = STAT_CALL)
if5 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function recvTaskCompleted, error calling UnitsManager, ERR03"
        end if if5

        open(unit   = iFile,                                                  &
             file   = trim(adjustl(TaskFile)) // '.tmp',                      &
             status = 'old',                                                  &
             iostat = STAT_CALL)
        if (STAT_CALL .EQ. SUCCESS_)                                          &
            close(unit   = iFile,                                             &
                  status = 'delete')

        STAT_CALL = isLastTask(Me)
if3 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
        end if if3

        recvTaskCompleted = SUCCESS_

    end function recvTaskCompleted

    !---------------------------------------------------------------------------

    integer function isLastTask(Me)
        !When there ar no more tasks the program will stop
        type (T_DDC), pointer       :: Me
        integer                     :: STAT_CALL  = UNKNOWN_
        integer                     :: countTasks = NULL_INT

        countTasks = hash_count(Me%hash_tasks)
if1 :   if (countTasks .EQ. 0) then
            STAT_CALL = sendWorkersPoisonPill(Me)
if16 :      if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function isLastTask, error calling sendWorkersPoisonPill, ERR02"
            end if if16

            call endPrg(Me)
        end if if1

        isLastTask = SUCCESS_

    end function isLastTask

    !---------------------------------------------------------------------------

    subroutine endPrg(Me)
        type (T_DDC), pointer       :: Me
        integer                     :: STAT_CALL  = UNKNOWN_

        STAT_CALL = killDDC(Me)
if2 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine endPrg, program DDCParser, error calling killDDC, ERR03"
        end if if2

        print*, "Program DDCParser terminated"

        call EXIT(SUCCESS_)

    end subroutine endPrg

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    integer function killDDC(Me)
        type (T_DDC), pointer       :: Me
        integer                     :: STAT_CALL = NULL_INT

        call deallocateInstance(Me)

!        STAT_CALL = barrier()
!if5 :   if (STAT_CALL .NE. SUCCESS_) then
!            print*, "STAT_CALL = ", STAT_CALL
!            stop "function killDDC, program DDCParser, error calling barrier, ERR05"
!        end if if5

        STAT_CALL = stopMPI()
if2 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function killDDC, error calling stopMPI, ERR02"
        end if if2

        killDDC = SUCCESS_

    end function killDDC

    !------------------------------------------------------------------------

    subroutine deallocateInstance(Me)
        type (T_DDC), pointer                       :: Me

        if (associated(Me%DirectoryList)) call deallocateDirectoryList(Me%DirectoryList)

        deallocate (Me)
        nullify    (Me)

    end subroutine deallocateInstance

    !------------------------------------------------------------------------

    recursive subroutine deallocateDirectoryList (DirectoryList)
        type(T_DirectoryList), pointer              :: DirectoryList

if1 :   if (associated(DirectoryList%Next)) then
            call KillHash_map           (DirectoryList%hash_map_out)
            nullify   (DirectoryList%hash_map_out)
            call KillHash_map           (DirectoryList%hash_map_in)
            nullify   (DirectoryList%hash_map_in)
            call deallocateDirectoryList(DirectoryList%Next)
        endif if1

        deallocate (DirectoryList)
        nullify    (DirectoryList)

    end subroutine deallocateDirectoryList

    !---------------------------------------------------------------------------

    integer function stopMPI()
        integer                     :: STAT_CALL = NULL_INT

        call MPI_FINALIZE(STAT_CALL)
if2 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function stopMPI, error calling MPI_FINALIZE, ERR02"
        end if if2

        stopMPI =  SUCCESS_

    end function stopMPI

end program DCCParser


