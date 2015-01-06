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
! /opt/mpich/bin/mpiexec -n 1 /myDirectory/projects/lang/fortran/DDC/src/Solutions/Linux/DCCWorker/DCCWorker : -n 2  /myDirectory/projects/lang/fortran/DDC/src/Solutions/Linux/DCCWorker/DCCWorker

program DDCWorker

    use IFPORT
    use mpi
    use moduleMPImanagement

    use ModuleGlobalData
    use ModuleHDF5_OO
    use HDF5

    implicit none

    !Constructor
!    private ::      constructDDC
!    private ::          startMPI
!    private ::          barrier

    !Sets & Gets
!    private ::      setMyMPI_id
!    private ::      getMyMPI_id
!    private ::      setParser_id
!    private ::      getParser_id
!    private ::      setConsolidatedFile
!    private ::      getConsolidatedFile

    !Modifier
!    private :: main
!    public  ::      Loop

    !Recv
!    private ::      recvIdle
!    private ::      recvTask
!    private ::          processTask
!    private ::              allocateTaskAuxMatrices
!    private ::              readInputFile
!    private ::                  readDecomposedFilesNames
!    private ::                  deallocateTaskAuxMatrices
!    private ::                      deallocateAuxMatrices
!    private ::              writeConsolidatedHDF
!    private ::                  allocateTaskMPIMatrices
!    private ::                  writeHDFDataSetLimits
!    private ::                  writeConsolidatedHDF2
!    private ::                      getDomainSize
!    private ::                          GetMappingValues
!    private ::                      readHDFDataSet
!    private ::                      auxArrayAllocation
!    private ::                      deallocateTaskMPIMatrices
!    private ::                          deallocateMPIMatrices
!    private ::                      getLimitsHDF
!    private ::                      getIndexMPIHDF
!    private ::                      getIndexAuxHDF
!    private ::                      mergeArraysMPI
!    private ::      recvParser_id
!    private ::      recvPoisonPill
!    private ::          endPrg ««-- Point of EXIT

    !Send
!    private ::          sendIdle
!    private ::          sendTaskCompleted

    !Destructor
!    private ::      killDDC
!    private ::          deallocateInstance
!    private ::          deallocateInputFile
!    private ::          stopMPI

    !Types---------------------------------------------------------------------
    type T_MPIMatrices
        real(4),    dimension(:    ), pointer :: MPIDataR4_1D
        real(4),    dimension(:,:  ), pointer :: MPIDataR4_2D
        real(4),    dimension(:,:,:), pointer :: MPIDataR4_3D
        real(8),    dimension(:    ), pointer :: MPIDataR8_1D
        real(8),    dimension(:,:  ), pointer :: MPIDataR8_2D
        real(8),    dimension(:,:,:), pointer :: MPIDataR8_3D
        integer,    dimension(:    ), pointer :: MPIDataI4_1D
        integer,    dimension(:,:  ), pointer :: MPIDataI4_2D
        integer,    dimension(:,:,:), pointer :: MPIDataI4_3D
    end type T_MPIMatrices

    type T_AuxMatrices
        real(4),    dimension(:    ), pointer :: AuxDataR4_1D
        real(4),    dimension(:,:  ), pointer :: AuxDataR4_2D
        real(4),    dimension(:,:,:), pointer :: AuxDataR4_3D
        real(8),    dimension(:    ), pointer :: AuxDataR8_1D
        real(8),    dimension(:,:  ), pointer :: AuxDataR8_2D
        real(8),    dimension(:,:,:), pointer :: AuxDataR8_3D
        integer,    dimension(:    ), pointer :: AuxDataI4_1D
        integer,    dimension(:,:  ), pointer :: AuxDataI4_2D
        integer,    dimension(:,:,:), pointer :: AuxDataI4_3D
    end type T_AuxMatrices

    type       T_Task
        type (T_AuxMatrices), pointer         :: AuxMatrices
        type (T_MPIMatrices), pointer         :: MPIMatrices
    end type  T_Task

    type       T_InputFile
        type (T_HDF5), pointer                :: inputFileHDF5

        type (T_InputFile), pointer           :: Next
    end type  T_InputFile

    type       T_DDC
        integer                               :: myMPI_id         = NULL_INT
        integer                               :: parser_id        = NULL_INT

        character(PathLength)                 :: ConsolidatedFileTmp = NULL_STR
        character(PathLength)                 :: ConsolidatedFile    = NULL_STR
        type (T_HDF5), pointer                :: outputFileHDF5
        type (T_InputFile), pointer           :: inputFile

        type (T_Task), pointer                :: task
    end type  T_DDC

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

        nullify(Me)
        Me => ConstructDDC()
        if (.NOT. associated(Me))                                             &
            stop "subroutine main, program DCCWorker, error calling ConstructDDC, ERR01"

        STAT_CALL = barrier()
if5 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine main, program DCCWorker, error calling barrier, ERR05"
        end if if5

        call Loop(Me)

    end subroutine main

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    function constructDDC()
        type (T_DDC), pointer       :: ConstructDDC
        type (T_DDC), pointer       :: NewObjDDC
        integer                     :: STAT_CALL  = UNKNOWN_

        allocate(NewObjDDC)
        nullify (NewObjDDC%task)
        nullify (NewObjDDC%inputFile)

        STAT_CALL = startMPI(NewObjDDC)
if1 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function ConstructDDC, program DCCWorker, error calling startMPI, ERR01"
        end if if1

        ConstructDDC => NewObjDDC

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
            stop "function startMPI, program DCCWorker, error calling MPI_INIT, ERR01"
        end if if1

        call MPI_COMM_RANK(MPI_COMM_WORLD,                             &
                           myMPI_id,                                   &
                           IERROR = STAT_CALL)
if2 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function startMPI, program DCCWorker, error calling MPI_COMM_RANK, ERR02.1"
        end if if2

        STAT_CALL = setMyMPI_id(Me, myMPI_id = myMPI_id)
if21 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function startMPI, program DCCWorker, error calling setMyMPI_id, ERR02.2"
        end if if21

        call MPI_COMM_SIZE(MPI_COMM_WORLD,                               &
                           numprocs,                                     &
                           IERROR = STAT_CALL)
if3 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function startMPI, program DDCWorker, error calling MPI_COMM_SIZE, ERR03"
        end if if3

        print *, 'Process ', myMPI_id, ' of ', numprocs, ' is alive -> Program DDCWorker'

        startMPI = SUCCESS_

    end function startMPI

    !---------------------------------------------------------------------------

    integer function barrier()
        integer                     :: STAT_CALL  = UNKNOWN_

        call MPI_BARRIER(MPI_COMM_WORLD, STAT_CALL)
if4 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function startMPI, program DCCWorker, error calling MPI_BARRIER, ERR04"
        end if if4

        barrier = SUCCESS_

    end function barrier


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

    integer function setParser_id(Me, parser_id)
        type (T_DDC), pointer       :: Me
        integer, intent(IN)         :: parser_id

        Me%parser_id = parser_id

        setParser_id = SUCCESS_

    end function setParser_id

    !---------------------------------------------------------------------------

    integer function getParser_id(Me)
        type (T_DDC), pointer       :: Me

        getParser_id = Me%parser_id

    end function getParser_id

    !---------------------------------------------------------------------------

    integer function setConsolidatedFile(Me, ConsolidatedFile)
        type (T_DDC), pointer             :: Me
        character(PathLength), intent(IN) :: ConsolidatedFile

        Me%ConsolidatedFile = ConsolidatedFile

        setConsolidatedFile = SUCCESS_

    end function setConsolidatedFile

    !---------------------------------------------------------------------------

    integer function setConsolidatedFileTmp(Me, ConsolidatedFileTmp)
        type (T_DDC), pointer             :: Me
        character(PathLength), intent(IN) :: ConsolidatedFileTmp

        Me%ConsolidatedFileTmp = ConsolidatedFileTmp

        setConsolidatedFileTmp = SUCCESS_

    end function setConsolidatedFileTmp

    !---------------------------------------------------------------------------

    function getConsolidatedFile(Me)
        character(PathLength)       :: getConsolidatedFile
        type (T_DDC), pointer       :: Me

        getConsolidatedFile = Me%ConsolidatedFile

    end function getConsolidatedFile

    !---------------------------------------------------------------------------

    function getConsolidatedFileTmp(Me)
        character(PathLength)       :: getConsolidatedFileTmp
        type (T_DDC), pointer       :: Me

        getConsolidatedFileTmp = Me%ConsolidatedFileTmp

    end function getConsolidatedFileTmp

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
            stop "subroutine Loop, program DCCWorker, error calling MPI_PROBE, ERR01"
        end if if25

if42 :  if      (STATUS(MPI_TAG) .EQ. getMsgSendTaskTag()) then

            STAT_CALL = RecvTask(Me)
if4 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "subroutine Loop, program DCCWorker, error calling RecvIdleWorker, ERR02"
            end if if4

        else if (STATUS(MPI_TAG) .EQ. getMsgPoisonPillTag()) then if42

            STAT_CALL = recvPoisonPill(Me)
if2 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "subroutine Loop, program DCCWorker, error calling recvPoisonPill, ERR03"
            end if if2

        else if (STATUS(MPI_TAG) .EQ. getMsgRankTag()) then if42

            STAT_CALL = recvParser_id(Me, STATUS(MPI_SOURCE))
if8 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "subroutine Loop, program DCCWorker, error calling recvParser_id, ERR04"
            end if if8

        else if42

            print*, "received message tag = ", STATUS(MPI_TAG)
            stop "subroutine Loop, program DCCWorker, message not recognized, ERR99"
        end if if42

        call Loop(Me)
    end subroutine Loop

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SEND SEND SEND SEND SEND SEND SEND SEND SEND SEND SEND SEND SEND SEND SEND S

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    integer function sendIdle(Me)
        type (T_DDC), pointer       :: Me
        integer                     :: STAT_CALL    = NULL_INT

        call MPI_SEND(getMyMPI_id(Me),                                        &
                      1,                                                      &
                      MPI_INTEGER,                                            &
                      getParser_id(Me),                                       &
                      getMsgIdleWorkerTag(),                                  &
                      MPI_COMM_WORLD,                                         &
                      STAT_CALL)
if16 :  if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function sendIdle, error calling MPI_SEND, ERR01"
        end if if16

        sendIdle = SUCCESS_

    end function sendIdle

    !---------------------------------------------------------------------------

    integer function sendTaskCompleted(Me, nTaskFile)
        type (T_DDC), pointer       :: Me
        integer                     :: nTaskFile
        integer                     :: STAT_CALL    = NULL_INT

        call MPI_SEND(nTaskFile,                                              &
                      1,                                                      &
                      MPI_INTEGER,                                            &
                      getParser_id(Me),                                       &
                      getMsgTaskCompletedTag(),                               &
                      MPI_COMM_WORLD,                                         &
                      STAT_CALL)
if16 :  if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function sendTaskCompleted, error calling MPI_SEND, ERR01"
        end if if16

        sendTaskCompleted = SUCCESS_

    end function sendTaskCompleted

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !RECV RECV RECV RECV RECV RECV RECV RECV RECV RECV RECV RECV RECV RECV RECV R

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    integer function recvTask(Me)
        type (T_DDC), pointer       :: Me
        integer                     :: STAT_CALL    = NULL_INT
        integer                     :: nTaskFile    = NULL_INT
        integer                     :: STATUS(MPI_STATUS_SIZE)

        call MPI_RECV(nTaskFile,                                                &
                      1,                                                        &
                      MPI_INTEGER,                                              &
                      getParser_id(Me),                                         &
                      getMsgSendTaskTag(),                                      &
                      MPI_COMM_WORLD,                                           &
                      STATUS,                                                   &
                      STAT_CALL)
if1 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function RecvTask, program DCCWorker, error calling MPI_SEND, ERR01"
        end if if1

        STAT_CALL = processTask(Me, nTaskFile = nTaskFile)
if9 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function RecvTask, program DCCWorker, error calling processTask, ERR09"
        end if if9

        STAT_CALL = renameOutputFile(consolidatedFile    = getConsolidatedFile(Me), &
                                     consolidatedFileTmp = getConsolidatedFileTmp(Me))
if6 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function RecvTask, program DCCWorker, error calling renameOutputFile, ERR06"
        end if if6

        STAT_CALL = SendTaskCompleted(Me, nTaskFile)
if4 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function RecvTask, program DCCWorker, error calling SendTaskCompleted, ERR02"
        end if if4

        print*, "File ", adjustl(trim(getConsolidatedFile(Me))), " is ready"

        if (associated(Me%inputFile))                                         &
            call deallocateInputFile(Me%inputFile)

        STAT_CALL = sendIdle(Me)
if5 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function RecvTask, program DCCWorker, error calling sendIdle, ERR04"
        end if if5

        recvTask = SUCCESS_

    end function recvTask

    !---------------------------------------------------------------------------

    integer function processTask(Me, nTaskFile)
        type (T_DDC), pointer       :: Me
        integer, intent(IN)         :: nTaskFile
        integer                     :: STAT_CALL    = NULL_INT
        character(PathLength)       :: HDF5File     = NULL_STR

        !Opens Me%outputFileHDF5
        HDF5File = readInputFile(Me, nTaskFile = nTaskFile)
        if (HDF5File .EQ. NULL_STR)                                           &
            stop "function processTask, program DCCWorker, error calling readInputFile, ERR0"


        STAT_CALL = writeConsolidatedHDF(outputFileHDF5 = Me%outputFileHDF5,  &
                                         inputFile      = Me%inputFile,       &
                                         task           = Me%task,            &
                                         GroupName      = "")
if3 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function processTask, program DCCWorker, error calling writeConsolidatedHDF, ERR02"
        end if if3

        STAT_CALL= KillHDF5(Me%outputFileHDF5)
if6 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function processTask, program DCCWorker, error calling KillHDF5, ERR06"
        end if if6

        processTask = SUCCESS_

    end function processTask

    !---------------------------------------------------------------------------

    function readInputFile(Me, nTaskFile)
    !returns the name of the consolidated HDF5 file
        character(PathLength)       :: readInputFile
        type (T_DDC), pointer       :: Me
        integer, intent(IN)         :: nTaskFile
        character(PathLength)       ::  TaskFile
        character(PathLength)       :: auxStr, auxStr2, auxStr3, auxStr4, auxStr5
        character(PathLength)       :: HDF5File     = NULL_STR
        integer                     :: STAT_CALL    = NULL_INT
        integer                     :: iFile        = NULL_INT
        integer                     :: HDF5_CREATE  = NULL_INT
        integer                     :: iPoint       = NULL_INT

        call UnitsManager(iFile, OPEN_FILE, STAT = STAT_CALL)
if1 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function readInputFile, program DCCWorker, error calling UnitsManager, ERR01"
        end if if1

        write (auxStr, '(i10)') nTaskFile
        TaskFile = trim(adjustl(auxStr)) // '.tmp'
        open(UNIT = iFile, FILE = TaskFile, status = 'OLD', IOSTAT = STAT_CALL)
if2 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function readInputFile, program DCCWorker, error calling open, ERR02"
        end if if2

        read(unit = iFile, fmt='(a256)', IOSTAT = STAT_CALL) auxStr2
if3 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function readInputFile, program DCCWorker, error calling read, ERR04"
        end if if3

        iPoint  = scan(auxStr2, ".", .TRUE.)
        auxStr4 = trim(adjustl(auxStr2(iPoint:)))
        auxStr5 = trim(adjustl(auxStr2(:iPoint)))

        write (auxStr3, '(i10)') getMyMPI_id(Me)
        HDF5File = adjustl(trim(auxStr5)) // 'WRK_' // trim(adjustl(auxStr3)) // trim(adjustl(auxStr4))
        call GetHDF5FileAccess(HDF5_CREATE = HDF5_CREATE)
        Me%outputFileHDF5 => ConstructHDF5(FileName = adjustl(trim(HDF5File)),  &
                                           Access   = HDF5_CREATE)

        STAT_CALL = setConsolidatedFile(Me, auxStr2)
if5 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function readInputFile, program DCCWorker, error calling setConsolidated, ERR05"
        end if if5

        STAT_CALL = setConsolidatedFileTmp(Me, HDF5File)
if8 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function readInputFile, program DCCWorker, error calling setConsolidated, ERR08"
        end if if8

        STAT_CALL = readDecomposedFilesNames(Me%inputFile, iFile = iFile)
if4 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function readInputFile, program DCCWorker, error calling readDecomposedFilesNames, ERR05"
        end if if4

        call UnitsManager(iFile, CLOSE_FILE, STAT = STAT_CALL)
if7 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function readInputFile, program DCCWorker, error calling UnitsManager, ERR03"
        end if if7

        readInputFile = HDF5File

    end function readInputFile

    !---------------------------------------------------------------------------

    recursive integer function readDecomposedFilesNames(inputFile, iFile)
        type (T_InputFile), pointer :: inputFile
        integer, intent(IN)         :: iFile
        integer                     :: STAT_CALL
        character(PathLength)       :: HDF5File     = NULL_STR
        integer                     :: HDF5_READ    = NULL_INT

        STAT_CALL = NULL_INT

        read(unit = iFile, fmt='(a256)', IOSTAT = STAT_CALL) HDF5File
if1 :   if (STAT_CALL .EQ. SUCCESS_) then
            allocate(inputFile)
            nullify (inputFile%Next)
            nullify (inputFile%inputFileHDF5)

            call GetHDF5FileAccess(HDF5_READ = HDF5_READ)
            inputFile%inputFileHDF5 => ConstructHDF5(FileName = adjustl(trim(HDF5File)),  &
                                                     Access   = HDF5_READ)

            STAT_CALL = readDecomposedFilesNames(inputFile%Next, iFile)
if2 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function readDecomposedFilesNames, error calling readDecomposedFilesNames, ERR02"
            end if if2
        end if if1

        readDecomposedFilesNames = SUCCESS_

    end function readDecomposedFilesNames

    !--------------------------------------------------------------------------

    recursive integer function writeConsolidatedHDF(outputFileHDF5,           &
                                                    inputFile,                &
                                                    task,                     &
                                                    GroupName)

        type (T_HDF5), pointer      :: outputFileHDF5
        type (T_InputFile), pointer :: inputFile
        type (T_Task), pointer      :: task
        character(*), intent(IN)    :: GroupName
        integer                     :: nItems
        integer                     :: idx
        character(StringLength)     :: NewGroupName
        character(StringLength)     :: GroupName2
        integer(HID_T)              :: GroupType
        integer                     :: LimitsArrayFactor
        integer                     :: STAT_CALL
        integer                     :: Rank
        integer, dimension(4)       :: DomainSize
        integer, dimension(7)       :: Dimensions
        character(StringLength)     :: Units
        integer(HID_T)              :: DataType

        nItems              = NULL_INT
        idx                 = NULL_INT
        Rank                = NULL_INT
        DomainSize          = NULL_INT
        Dimensions          = NULL_INT
        GroupName2          = NULL_STR
        GroupType           = NULL_INT
        LimitsArrayFactor   = NULL_INT
        NewGroupName        = NULL_STR
        Units               = NULL_STR
        STAT_CALL           = NULL_INT

        !Get the number of members in the Group
        nItems = GetHDF5GroupNumberOfItems(inputFile%inputFileHDF5,           &
                                           GroupName = adjustl(trim(GroupName))//"/")

do1 :   do idx = 1, nItems
            !Gets information about the group
            call GetHDF5ObjectInfo(inputFile%inputFileHDF5,                         &
                                   FatherGroupName = adjustl(trim(GroupName))//"/", &
                                   GroupPosition   = idx,                           &
                                   GroupName       = GroupName2,                    &
                                   GroupType       = GroupType,                     &
                                   STAT            = STAT_CALL)
if7 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "subroutine writeConsolidatedHDF, program DCCWorker, error calling GetHDF5GroupNumberOfItems, ERR01"
            end if if7

if1 :       if ((GroupType .EQ. H5G_DATASET_F) .AND.                          &
                (adjustl(trim(GroupName))//"/"//adjustl(trim(GroupName2)) .NE. "/Grid/Decomposition")) then

                STAT_CALL = allocateTaskAuxMatrices(task)
if2 :           if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "function writeConsolidatedHDF, program DCCWorker, error calling allocateTaskAuxMatrices, ERR06"
                end if if2

                !Check if it is Coordenates Feature
if12 :          if ((adjustl(trim(GroupName2)) .EQ. "ConnectionX") .OR.       &
                    (adjustl(trim(GroupName2)) .EQ. "ConnectionY") .OR.       &
                    (adjustl(trim(GroupName2)) .EQ. "Latitude"   ) .OR.       &
                    (adjustl(trim(GroupName2)) .EQ. "Longitude"  )) then

                    LimitsArrayFactor = 1
                else if12
                    LimitsArrayFactor = 0
                endif if12

                call writeConsolidatedHDF2(outputFileHDF5,                    &
                                           inputFile         = inputFile,     &
                                           task              = task,          &
                                           GroupName         = GroupName,     &
                                           GroupName2        = GroupName2,    &
                                           idx               = idx,           &
                                           Rank              = Rank,          &
                                           DomainSize        = DomainSize,    &
                                           Dimensions        = Dimensions,    &
                                           Units             = Units,         &
                                           DataType          = DataType,      &
                                           LimitsArrayFactor = LimitsArrayFactor)

                call writeHDFDataSetLimits(outputFileHDF5,                        &
                                           Rank              = Rank,              &
                                           DomainSize        = DomainSize,        &
                                           Dimensions        = Dimensions,        &
                                           LimitsArrayFactor = LimitsArrayFactor, &
                                           GroupName         = GroupName)

                call writeHDFDataSet(outputFileHDF5,                          &
                                     GroupName   = GroupName,                 &
                                     Name        = GroupName2,                &
                                     Units       = Units,                     &
                                     Rank        = Rank,                      &
                                     DataType    = DataType,                  &
                                     AuxMatrices = task%AuxMatrices)

                STAT_CALL = deallocateTaskAuxMatrices(task)
if6 :           if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "function writeConsolidatedHDF, program DCCWorker, error calling deallocateTaskAuxMatrices, ERR02"
                end if if6

            elseif (GroupType .EQ. H5G_GROUP_F) then if1
if3 :           if      (adjustl(trim(GroupName))                                 .EQ. "/") then
                    NewGroupName = adjustl(trim(GroupName))//adjustl(trim(GroupName2))
                else if (adjustl(trim(GroupName))//"/"//adjustl(trim(GroupName2)) .NE. "/Grid/Decomposition") then if3
                    NewGroupName = adjustl(trim(GroupName))//"/"//adjustl(trim(GroupName2))

                    STAT_CALL = writeConsolidatedHDF(outputFileHDF5 = outputFileHDF5, &
                                                     inputFile      = inputFile,      &
                                                     task           = task,           &
                                                     GroupName      = adjustl(trim(NewGroupName)))
                endif if3
            endif if1
        end do do1

        WriteConsolidatedHDF = SUCCESS_

    end function writeConsolidatedHDF

    !--------------------------------------------------------------------------

    subroutine writeHDFDataSet(outputFileHDF5,                                &
                               GroupName,                                     &
                               Name ,                                         &
                               Units,                                         &
                               Rank,                                          &
                               DataType,                                      &
                               AuxMatrices)

        type (T_HDF5), pointer        :: outputFileHDF5
        character(*), intent(IN)      :: GroupName
        character(*), intent(IN)      :: Name
        character(*), intent(IN)      :: Units
        integer, intent(IN)           :: Rank
        integer(HID_T), intent(IN)    :: DataType
        type(T_AuxMatrices), pointer  :: AuxMatrices
        integer                       :: STAT_CALL      = NULL_INT

if1 :   if (Rank .EQ. 1) then
if12 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then
                call HDF5WriteData(outputFileHDF5,                            &
                                   GroupName   = adjustl(trim(GroupName)),    &
                                   Name        = adjustl(trim(Name)),         &
                                   Units       = adjustl(trim(Units)),        &
                                   Array1D     = AuxMatrices%AuxDataI4_1D,    &
                                   STAT        = STAT_CALL)
if8 :           if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "subroutine writeHDFDataSet, program DCCWorker, error calling HDF5WriteData, ERR03"
                end if if8

            else if (DataType .EQ. H5T_NATIVE_REAL) then if12
                call HDF5WriteData(outputFileHDF5,                            &
                                   GroupName   = adjustl(trim(GroupName)),    &
                                   Name        = adjustl(trim(Name)),         &
                                   Units       = adjustl(trim(Units)),        &
                                   Array1D     = AuxMatrices%AuxDataR4_1D,    &
                                   STAT        = STAT_CALL)
if9 :           if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "subroutine writeHDFDataSet, program DCCWorker, error calling HDF5WriteData, ERR04"
                end if if9

            else if (DataType .EQ. H5T_NATIVE_DOUBLE) then if12
                call HDF5WriteData(outputFileHDF5,                            &
                                   GroupName   = adjustl(trim(GroupName)),    &
                                   Name        = adjustl(trim(Name)),         &
                                   Units       = adjustl(trim(Units)),        &
                                   Array1D     = AuxMatrices%AuxDataR8_1D,    &
                                   STAT        = STAT_CALL)
if10 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "subroutine writeHDFDataSet, program DCCWorker, error calling HDF5WriteData, ERR05"
                end if if10
            end if if12

        else if (Rank .EQ. 2) then if1
if13 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then
                call HDF5WriteData(outputFileHDF5,                            &
                                   GroupName   = adjustl(trim(GroupName)),    &
                                   Name        = adjustl(trim(Name)),         &
                                   Units       = adjustl(trim(Units)),        &
                                   Array2D     = AuxMatrices%AuxDataI4_2D,    &
                                   STAT        = STAT_CALL)
if11 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "subroutine writeHDFDataSet, program DCCWorker, error calling HDF5WriteData, ERR06"
                end if if11

            else if (DataType .EQ. H5T_NATIVE_REAL) then if13
                call HDF5WriteData(outputFileHDF5,                            &
                                   GroupName   = adjustl(trim(GroupName)),    &
                                   Name        = adjustl(trim(Name)),         &
                                   Units       = adjustl(trim(Units)),        &
                                   Array2D     = AuxMatrices%AuxDataR4_2D,    &
                                   STAT        = STAT_CALL)
if31 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "subroutine writeHDFDataSet, program DCCWorker, error calling HDF5WriteData, ERR07"
                end if if31

            else if (DataType .EQ. H5T_NATIVE_DOUBLE) then if13
                call HDF5WriteData(outputFileHDF5,                            &
                                   GroupName   = adjustl(trim(GroupName)),    &
                                   Name        = adjustl(trim(Name)),         &
                                   Units       = adjustl(trim(Units)),        &
                                   Array2D     = AuxMatrices%AuxDataR8_2D,    &
                                   STAT        = STAT_CALL)
if38 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "subroutine writeHDFDataSet, program DCCWorker, error calling HDF5WriteData, ERR08"
                end if if38
            end if if13
        else if (Rank .EQ. 3) then if1
if14 :      if (DataType .EQ. H5T_NATIVE_INTEGER) then
                call HDF5WriteData(outputFileHDF5,                            &
                                   GroupName   = adjustl(trim(GroupName)),    &
                                   Name        = adjustl(trim(Name)),         &
                                   Units       = adjustl(trim(Units)),        &
                                   Array3D     = AuxMatrices%AuxDataI4_3D,    &
                                   STAT        = STAT_CALL)
if30 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "subroutine writeHDFDataSet, program DCCWorker, error calling HDF5WriteData, ERR09"
                end if if30

            else if (DataType .EQ. H5T_NATIVE_REAL) then if14
                call HDF5WriteData(outputFileHDF5,                            &
                                   GroupName   = adjustl(trim(GroupName)),    &
                                   Name        = adjustl(trim(Name)),         &
                                   Units       = adjustl(trim(Units)),        &
                                   Array3D     = AuxMatrices%AuxDataR4_3D,    &
                                   STAT        = STAT_CALL)
if32 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "subroutine writeHDFDataSet, program DCCWorker, error calling HDF5WriteData, ERR10"
                end if if32

            else if (DataType .EQ. H5T_NATIVE_DOUBLE) then if14
                call HDF5WriteData(outputFileHDF5,                            &
                                   GroupName   = adjustl(trim(GroupName)),    &
                                   Name        = adjustl(trim(Name)),         &
                                   Units       = adjustl(trim(Units)),        &
                                   Array3D     = AuxMatrices%AuxDataR8_3D,    &
                                   STAT        = STAT_CALL)
if7 :           if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "subroutine writeHDFDataSet, program DCCWorker, error calling HDF5WriteData, ERR02"
                end if if7
            end if if14
        end if if1

        !Writes everything to disk
        call HDF5FlushMemory(outputFileHDF5, STAT = STAT_CALL)
if6 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine writeHDFDataSet, program DCCWorker, error calling HDF5FlushMemory, ERR01"
        end if if6

    end subroutine writeHDFDataSet

    !--------------------------------------------------------------------------

    recursive subroutine writeHDFDataSetLimits(outputFileHDF5,                &
                                               Rank,                          &
                                               DomainSize,                    &
                                               Dimensions,                    &
                                               LimitsArrayFactor,             &
                                               GroupName)

        type (T_HDF5), pointer            :: outputFileHDF5
        integer, intent(IN)               :: Rank
        integer, dimension(4), intent(IN) :: DomainSize
        integer, dimension(7), intent(IN) :: Dimensions
        integer, intent(IN)               :: LimitsArrayFactor
        CHARACTER(*), intent(IN)          :: GroupName
        integer                                                 :: STAT_CALL
        integer                                                 :: ILB, IUB
        integer                                                 :: JLB, JUB
        integer                                                 :: KLB, KUB

if21 :  if (Rank .EQ. 1) then
if22 :      if (DomainSize(1) .GT. 1 .AND. adjustl(trim(GroupName)) .NE. "/Time") then
                ILB = DomainSize(1)-DomainSize(1)+1
                IUB = DomainSize(2)-DomainSize(1)+1
                JLB = 1
                JUB = 1
                KLB = 1
                KUB = 1

            else if (adjustl(trim(GroupName)) .EQ. "/Time") then if22
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
        else if (Rank .EQ. 2) then if21
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
        else if (Rank .EQ. 3) then if21
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

        call HDF5SetLimits(outputFileHDF5,                                    &
                           ILB = ILB, IUB = IUB,                              &
                           JLB = JLB, JUB = JUB,                              &
                           KLB = KLB, KUB = KUB,                              &
                           STAT = STAT_CALL)
if6 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine writeHDFDataSetLimits, program DCCWorker, error calling HDF5SetLimits, ERR01"
        end if if6

    end subroutine writeHDFDataSetLimits

    !--------------------------------------------------------------------------

    recursive subroutine writeConsolidatedHDF2(outputFileHDF5,                &
                                               inputFile,                     &
                                               task,                          &
                                               GroupName,                     &
                                               GroupName2,                    &
                                               idx,                           &
                                               Rank,                          &
                                               DomainSize,                    &
                                               Dimensions,                    &
                                               Units,                         &
                                               DataType,                      &
                                               LimitsArrayFactor)
    !Iterates input files

        type (T_HDF5), pointer                :: outputFileHDF5
        type (T_InputFile), pointer           :: inputFile
        type (T_Task), pointer                :: task
        character(*), intent(IN)              :: GroupName
        character(*), intent(IN)              :: GroupName2
        integer, intent(IN)                   :: idx
        integer, intent(IN)                   :: LimitsArrayFactor
        integer, intent(OUT)                  :: Rank
        integer, dimension(4), intent(OUT)    :: DomainSize
        integer, dimension(7), intent(OUT)    :: Dimensions
        character(StringLength), intent(OUT)  :: Units
        integer(HID_T), intent(OUT)           :: DataType
        integer, dimension(4)                 :: WindowPosition
        integer, dimension(4)                 :: WindowFrame
        integer                               :: STAT_CALL

        integer                               :: ILB_MPI, IUB_MPI
        integer                               :: JLB_MPI, JUB_MPI
        integer                               :: KLB_MPI, KUB_MPI
        integer                               :: ILB, IUB
        integer                               :: JLB, JUB
        integer                               :: KLB, KUB
        integer, dimension(6)                 :: LimitsArray


        DomainSize      = NULL_INT
        WindowPosition  = NULL_INT
        WindowFrame     = NULL_INT
        Dimensions      = NULL_INT
        STAT_CALL       = NULL_INT
        Rank            = NULL_INT
        Units           = NULL_STR
        DataType        = NULL_INT

        ILB_MPI         = NULL_INT
        IUB_MPI         = NULL_INT
        JLB_MPI         = NULL_INT
        JUB_MPI         = NULL_INT
        KLB_MPI         = NULL_INT
        KUB_MPI         = NULL_INT
        ILB             = NULL_INT
        IUB             = NULL_INT
        JLB             = NULL_INT
        JUB             = NULL_INT
        KLB             = NULL_INT
        KUB             = NULL_INT
        LimitsArray     = NULL_INT

        DomainSize      = getDomainSize(inputFile%inputFileHDF5, key = 1)
        WindowPosition  = getDomainSize(inputFile%inputFileHDF5, key = 2)
        WindowFrame     = getDomainSize(inputFile%inputFileHDF5, key = 3)

        !Get Dataset Information
        call GetHDF5GroupID(inputFile%inputFileHDF5,                          &
                            FatherGroupName = adjustl(trim(GroupName))//"/",  &
                            GroupPosition   = idx,                            &
                            GroupName       = GroupName2,                     &
                            Units           = Units,                          &
                            Rank            = Rank,                           &
                            Dimensions      = Dimensions,                     &
                            STAT            = STAT_CALL)
if7 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine writeConsolidatedHDF2, program DCCWorker, error calling GetHDF5GroupID, ERR01"
        end if if7

        !Get Dataset DataType
        call GetHDF5DataTypeID(inputFile%inputFileHDF5,                         &
                               FatherGroupName = adjustl(trim(GroupName))//"/", &
                               GroupPosition   = idx,                           &
                               GroupName       = GroupName2,                    &
                               DataType        = DataType,                      &
                               STAT            = STAT_CALL)
if8 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine writeConsolidatedHDF2, program DCCWorker, error calling GetHDF5DataTypeID, ERR02"
        end if if8

        STAT_CALL = allocateTaskMPIMatrices(task)
if81 :  if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine writeConsolidatedHDF2, program DCCWorker, error calling allocateTaskMPIMatrices, ERR81"
        end if if81

        STAT_CALL = readHDFDataSet(inputFile%inputFileHDF5,                   &
                                   GroupName   = adjustl(trim(GroupName)),    &
                                   obj_name    = adjustl(trim(GroupName2)),   &
                                   Rank        = Rank,                        &
                                   dims        = Dimensions,                  &
                                   DataType    = DataType,                    &
                                   MPIMatrices = task%MPIMatrices)
if9 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine writeConsolidatedHDF2, program DCCWorker, error calling ReadHDFDataSet, ERR03"
        end if if9

        STAT_CALL = auxArrayAllocation(DomainSize        = DomainSize,        &
                                       Dimensions        = Dimensions,        &
                                       AuxMatrices       = task%AuxMatrices,  &
                                       LimitsArrayFactor = LimitsArrayFactor, &
                                       MPIMatrices       = task%MPIMatrices)
if10 :  if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine writeConsolidatedHDF2, program DCCWorker, error calling AuxArrayAllocation, ERR04"
        end if if10

        LimitsArray =  getLimitsHDF(WindowPosition      = WindowPosition,     &
                                    Dimensions          = Dimensions,         &
                                    Rank                = Rank,               &
                                    LimitsArrayFactor   = LimitsArrayFactor,  &
                                    GroupName           = GroupName)


        call getIndexMPIHDF(WindowFrame = WindowFrame,                        &
                            LimitsArray = LimitsArray,                        &
                            GroupName   = GroupName,                          &
                            ILB_MPI     = ILB_MPI, IUB_MPI = IUB_MPI,         &
                            JLB_MPI     = JLB_MPI, JUB_MPI = JUB_MPI,         &
                            KLB_MPI     = KLB_MPI, KUB_MPI = KUB_MPI)

        call getIndexAuxHDF(WindowPosition      = WindowPosition,             &
                            DomainSize          = DomainSize,                 &
                            LimitsArray         = LimitsArray,                &
                            LimitsArrayFactor   = LimitsArrayFactor,          &
                            GroupName           = GroupName,                  &
                            ILB                 = ILB , IUB  = IUB,           &
                            JLB                 = JLB , JUB  = JUB,           &
                            KLB                 = KLB , KUB  = KUB)

        call mergeArraysMPI(ILB_MPI = ILB_MPI, IUB_MPI = IUB_MPI,             &
                            JLB_MPI = JLB_MPI, JUB_MPI = JUB_MPI,             &
                            KLB_MPI = KLB_MPI, KUB_MPI = KUB_MPI,             &
                            ILB = ILB, IUB = IUB,                             &
                            JLB = JLB, JUB = JUB,                             &
                            KLB = KLB, KUB = KUB,                             &
                            AuxMatrices = task%AuxMatrices,                   &
                            MPIMatrices = task%MPIMatrices)

        STAT_CALL = deallocateTaskMPIMatrices(task)
if82 :  if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine writeConsolidatedHDF2, program DCCWorker, error calling deallocateTaskMPIMatrices, ERR82"
        end if if82

        if (associated(inputFile%Next))                                       &
            call WriteConsolidatedHDF2(outputFileHDF5,                        &
                                       inputFile         = inputFile%NExt,    &
                                       task              = task,              &
                                       GroupName         = GroupName,         &
                                       GroupName2        = GroupName2,        &
                                       idx               = idx,               &
                                       Rank              = Rank,              &
                                       DomainSize        = DomainSize,        &
                                       Dimensions        = Dimensions,        &
                                       Units             = Units,             &
                                       DataType          = DataType,          &
                                       LimitsArrayFactor = LimitsArrayFactor)

    end subroutine WriteConsolidatedHDF2

    !------------------------------------------------------------------------

    integer function allocateTaskAuxMatrices(task)
        type(T_Task), pointer       :: task

        if (.NOT. associated(task))                                           &
            allocate(task)

        allocate(task%AuxMatrices)

        nullify(task%AuxMatrices%AuxDataR4_1D)
        nullify(task%AuxMatrices%AuxDataR4_2D)
        nullify(task%AuxMatrices%AuxDataR4_3D)
        nullify(task%AuxMatrices%AuxDataR8_1D)
        nullify(task%AuxMatrices%AuxDataR8_2D)
        nullify(task%AuxMatrices%AuxDataR8_3D)
        nullify(task%AuxMatrices%AuxDataI4_1D)
        nullify(task%AuxMatrices%AuxDataI4_2D)
        nullify(task%AuxMatrices%AuxDataI4_3D)

        allocateTaskAuxMatrices = SUCCESS_

    end function allocateTaskAuxMatrices

    !--------------------------------------------------------------------------

    function getLimitsHDF(WindowPosition,                                     &
                          Dimensions,                                         &
                          Rank,                                               &
                          LimitsArrayFactor,                                  &
                          GroupName)

        integer, dimension(6)             :: getLimitsHDF
        integer, dimension(4),intent(IN)  :: WindowPosition
        integer, dimension(7),intent(IN)  :: Dimensions
        integer,intent(IN)                :: Rank
        integer,intent(IN)                :: LimitsArrayFactor
        character(*), intent(IN)          :: GroupName
        integer, dimension(6)             :: auxLimitsArray

        auxLimitsArray = 1

if14 :  if      (Rank .EQ. 1 .AND. adjustl(trim(GroupName)) .NE. "/Time") then
            auxLimitsArray(1)= WindowPosition(1)
            auxLimitsArray(2)= WindowPosition(2)

        else if (Rank .EQ. 1 .AND. adjustl(trim(GroupName)) .EQ. "/Time") then if14
            auxLimitsArray(1)= 1
            auxLimitsArray(2)= Dimensions(1)

        else if (Rank .EQ. 2) then if14
            auxLimitsArray(1)= WindowPosition(1)
            auxLimitsArray(2)= WindowPosition(2)+LimitsArrayFactor
            auxLimitsArray(3)= WindowPosition(3)
            auxLimitsArray(4)= WindowPosition(4)+LimitsArrayFactor

        else if (Rank .EQ. 3) then if14
            auxLimitsArray(1)= WindowPosition(1)
            auxLimitsArray(2)= WindowPosition(2)
            auxLimitsArray(3)= WindowPosition(3)
            auxLimitsArray(4)= WindowPosition(4)
            auxLimitsArray(5)= 1
            auxLimitsArray(6)= Dimensions(3)
        end if if14

        getLimitsHDF = auxLimitsArray

    end function getLimitsHDF

    !------------------------------------------------------------------------

    integer function allocateTaskMPIMatrices(task)
        type(T_Task), pointer       :: task

        if (.NOT. associated(task))                                           &
            allocate(task)

        allocate(task%MPIMatrices)

        nullify(task%MPIMatrices%MPIDataR4_1D)
        nullify(task%MPIMatrices%MPIDataR4_2D)
        nullify(task%MPIMatrices%MPIDataR4_3D)
        nullify(task%MPIMatrices%MPIDataR8_1D)
        nullify(task%MPIMatrices%MPIDataR8_2D)
        nullify(task%MPIMatrices%MPIDataR8_3D)
        nullify(task%MPIMatrices%MPIDataI4_1D)
        nullify(task%MPIMatrices%MPIDataI4_2D)
        nullify(task%MPIMatrices%MPIDataI4_3D)

        allocateTaskMPIMatrices = SUCCESS_

    end function allocateTaskMPIMatrices

    !--------------------------------------------------------------------------

    subroutine getIndexAuxHDF(WindowPosition,                                 &
                              DomainSize,                                     &
                              LimitsArray,                                    &
                              LimitsArrayFactor,                              &
                              GroupName,                                      &
                              ILB, IUB,                                       &
                              JLB, JUB,                                       &
                              KLB, KUB)

        !Arguments-------------------------------------------------------------
        integer, dimension(4),intent(IN)  :: WindowPosition
        integer, dimension(4),intent(IN)  :: DomainSize
        integer, dimension(6),intent(IN)  :: LimitsArray
        integer,intent(IN)                :: LimitsArrayFactor
        character(*), intent(IN)          :: GroupName
        integer,intent(out)               :: ILB, IUB
        integer,intent(out)               :: JLB, JUB
        integer,intent(out)               :: KLB, KUB

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
        end if if16

    end subroutine getIndexAuxHDF

    !------------------------------------------------------------------------

    subroutine mergeArraysMPI(ILB_MPI, IUB_MPI,                               &
                              JLB_MPI, JUB_MPI,                               &
                              KLB_MPI, KUB_MPI,                               &
                              ILB,     IUB,                                   &
                              JLB,     JUB,                                   &
                              KLB,     KUB,                                   &
                              AuxMatrices,                                    &
                              MPIMatrices)

        !Arguments-------------------------------------------------------------
        type(T_AuxMatrices), pointer  :: AuxMatrices
        type(T_MPIMatrices), pointer  :: MPIMatrices
        integer,intent(IN)            :: ILB_MPI, IUB_MPI
        integer,intent(IN)            :: JLB_MPI, JUB_MPI
        integer,intent(IN)            :: KLB_MPI, KUB_MPI
        integer,intent(IN)            :: ILB,     IUB
        integer,intent(IN)            :: JLB,     JUB
        integer,intent(IN)            :: KLB,     KUB

        if(associated(MPIMatrices%MPIDataI4_1D))                                          &
            AuxMatrices%AuxDataI4_1D(ILB:IUB) = MPIMatrices%MPIDataI4_1D(ILB_MPI:IUB_MPI)

        if (associated(MPIMatrices%MPIDataI4_2D))                                         &
            AuxMatrices%AuxDataI4_2D(ILB:IUB,                                             &
                                     JLB:JUB) = MPIMatrices%MPIDataI4_2D(ILB_MPI:IUB_MPI, &
                                                                         JLB_MPI:JUB_MPI)

        if (associated(MPIMatrices%MPIDataI4_3D))                                         &
            AuxMatrices%AuxDataI4_3D(ILB:IUB,                                             &
                                    JLB:JUB,                                              &
                                    KLB:KUB) = MPIMatrices%MPIDataI4_3D(ILB_MPI:IUB_MPI,  &
                                                                        JLB_MPI:JUB_MPI,  &
                                                                        KLB_MPI:KUB_MPI)

        if(associated(MPIMatrices%MPIDataR4_1D))                                          &
            AuxMatrices%AuxDataR4_1D(ILB:IUB) = MPIMatrices%MPIDataR4_1D(ILB_MPI:IUB_MPI)

        if (associated(MPIMatrices%MPIDataR4_2D))                                         &
           AuxMatrices%AuxDataR4_2D(ILB:IUB,                                              &
                                    JLB:JUB) = MPIMatrices%MPIDataR4_2D(ILB_MPI:IUB_MPI,  &
                                                                        JLB_MPI:JUB_MPI)

        if (associated(MPIMatrices%MPIDataR4_3D))                                         &
            AuxMatrices%AuxDataR4_3D(ILB:IUB,                                             &
                                     JLB:JUB,                                             &
                                     KLB:KUB) = MPIMatrices%MPIDataR4_3D(ILB_MPI:IUB_MPI, &
                                                                         JLB_MPI:JUB_MPI, &
                                                                         KLB_MPI:KUB_MPI)

        if (associated(MPIMatrices%MPIDataR8_1D))                                         &
            AuxMatrices%AuxDataR8_1D(ILB:IUB) = MPIMatrices%MPIDataR8_1D(ILB_MPI:IUB_MPI)

        if (associated(MPIMatrices%MPIDataR8_2D))                                         &
            AuxMatrices%AuxDataR8_2D(ILB:IUB,                                             &
                                     JLB:JUB) = MPIMatrices%MPIDataR8_2D(ILB_MPI:IUB_MPI, &
                                                                         JLB_MPI:JUB_MPI)

        if (associated(MPIMatrices%MPIDataR8_3D))                                         &
            AuxMatrices%AuxDataR8_3D(ILB:IUB,                                             &
                                    JLB:JUB,                                              &
                                    KLB:KUB) = MPIMatrices%MPIDataR8_3D(ILB_MPI:IUB_MPI,  &
                                                                        JLB_MPI:JUB_MPI,  &
                                                                        KLB_MPI:KUB_MPI)


        !------------------------------------------------------------------------

    end subroutine mergeArraysMPI

    !--------------------------------------------------------------------------

    integer function readHDFDataSet(inputFileHDF5,                            &
                                    GroupName,                                &
                                    obj_name,                                 &
                                    Rank,                                     &
                                    dims,                                     &
                                    DataType,                                 &
                                    MPIMatrices)

        type (T_HDF5), pointer            :: inputFileHDF5
        character(*), intent(IN)          :: GroupName
        character(*), intent(IN)          :: obj_name
        integer, intent(IN)               :: Rank
        integer, dimension(7), intent(IN) :: dims
        integer(HID_T), intent(IN)        :: DataType
        type (T_MPIMatrices), pointer     :: MPIMatrices
        integer                           :: STAT_CALL          = NULL_INT

if1 :   if      (Rank .EQ. 1) then
            call HDF5SetLimits(inputFileHDF5,                                 &
                               ILB = 1, IUB = dims(1),                        &
                               STAT = STAT_CALL)
if20 :      if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function ReadHDFDataSet, program DCCWorker, error calling HDF5SetLimits, ERR01"
            end if if20

if12 :      if      (DataType .EQ. H5T_NATIVE_INTEGER) then
                if (.NOT. associated(MPIMatrices%MPIDataI4_1D))               &
                    allocate(MPIMatrices%MPIDataI4_1D(1:dims(1)))

                call HDF5ReadData(inputFileHDF5,                              &
                                  GroupName = GroupName,                      &
                                  Name      = adjustl(trim(obj_name)),        &
                                  Array1D   = MPIMatrices%MPIDataI4_1D,       &
                                  STAT      = STAT_CALL)
if21 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "function ReadHDFDataSet, program DCCWorker, error calling HDF5ReadData, ERR02"
                end if if21

            else if (DataType .EQ. H5T_NATIVE_REAL) then if12
                if (.NOT. associated(MPIMatrices%MPIDataR4_1D))               &
                    allocate(MPIMatrices%MPIDataR4_1D(1:dims(1)))

                call HDF5ReadData(inputFileHDF5,                              &
                                  GroupName = GroupName,                      &
                                  Name      = adjustl(trim(obj_name)),        &
                                  Array1D   = MPIMatrices%MPIDataR4_1D,       &
                                  STAT      = STAT_CALL)
if22 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "function ReadHDFDataSet, program DCCWorker, error calling HDF5ReadData, ERR03"
                end if if22

            else if (DataType .EQ. H5T_NATIVE_DOUBLE) then if12
                if (.NOT. associated(MPIMatrices%MPIDataR8_1D))               &
                    allocate(MPIMatrices%MPIDataR8_1D(1:dims(1)))

                call HDF5ReadData(inputFileHDF5,                              &
                                  GroupName = GroupName,                      &
                                  Name      = adjustl(trim(obj_name)),        &
                                  Array1D   = MPIMatrices%MPIDataR8_1D,       &
                                  STAT      = STAT_CALL)
if32 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "function ReadHDFDataSet, program DCCWorker, error calling HDF5ReadData, ERR04"
                end if if32
            end if if12

        else if (Rank .EQ. 2) then if1
            call HDF5SetLimits(inputFileHDF5,                                 &
                               ILB = 1, IUB = dims(1),                        &
                               JLB = 1, JUB = dims(2),                        &
                               STAT = STAT_CALL)
if23 :      if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function ReadHDFDataSet, program DCCWorker, error calling HDF5SetLimits, ERR05"
            end if if23

if13 :      if      (DataType .EQ. H5T_NATIVE_INTEGER) then
                if (.NOT. associated(MPIMatrices%MPIDataI4_2D))               &
                    allocate(MPIMatrices%MPIDataI4_2D(1:dims(1),              &
                                                      1:dims(2)))

                call HDF5ReadData(inputFileHDF5,                              &
                                  GroupName = GroupName,                      &
                                  Name      = adjustl(trim(obj_name)),        &
                                  Array2D   = MPIMatrices%MPIDataI4_2D,       &
                                  STAT      = STAT_CALL)
if24 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "function ReadHDFDataSet, program DCCWorker, error calling HDF5ReadData, ERR06"
                end if if24

            else if (DataType .EQ. H5T_NATIVE_REAL) then if13
                if (.NOT. associated(MPIMatrices%MPIDataR4_2D))               &
                    allocate(MPIMatrices%MPIDataR4_2D(1:dims(1),              &
                                                      1:dims(2)))

                call HDF5ReadData(inputFileHDF5,                              &
                                  GroupName = GroupName,                      &
                                  Name      = adjustl(trim(obj_name)),        &
                                  Array2D   = MPIMatrices%MPIDataR4_2D,       &
                                  STAT      = STAT_CALL)
if25 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "function ReadHDFDataSet, program DCCWorker, error calling HDF5ReadData, ERR07"
                end if if25

            else if (DataType .EQ. H5T_NATIVE_DOUBLE) then if13
                if (.NOT. associated(MPIMatrices%MPIDataR8_2D))               &
                    allocate(MPIMatrices%MPIDataR8_2D(1:dims(1),              &
                                                      1:dims(2)))

                call HDF5ReadData(inputFileHDF5,                              &
                                  GroupName = GroupName,                      &
                                  Name      = adjustl(trim(obj_name)),        &
                                  Array2D   = MPIMatrices%MPIDataR8_2D,       &
                                  STAT      = STAT_CALL)
if26 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "function ReadHDFDataSet, program DCCWorker, error calling HDF5ReadData, ERR08"
                end if if26
            end if if13
        else if (Rank .EQ. 3) then if1
            call HDF5SetLimits(inputFileHDF5,                                 &
                               ILB = 1, IUB = dims(1),                        &
                               JLB = 1, JUB = dims(2),                        &
                               KLB = 1, KUB = dims(3),                        &
                               STAT = STAT_CALL)
if27 :      if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function ReadHDFDataSet, program DCCWorker, error calling HDF5SetLimits, ERR09"
            end if if27

if14 :      if      (DataType .EQ. H5T_NATIVE_INTEGER) then
                if (.NOT. associated(MPIMatrices%MPIDataI4_3D))               &
                    allocate(MPIMatrices%MPIDataI4_3D(1:dims(1),              &
                                                      1:dims(2),              &
                                                      1:dims(3)))

                call HDF5ReadData(inputFileHDF5,                              &
                                  GroupName = GroupName,                      &
                                  Name      = adjustl(trim(obj_name)),        &
                                  Array3D   = MPIMatrices%MPIDataI4_3D,       &
                                  STAT      = STAT_CALL)
if28 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "function ReadHDFDataSet, program DCCWorker, error calling HDF5ReadData, ERR10"
                end if if28

            else if (DataType .EQ. H5T_NATIVE_REAL) then if14
                if (.NOT. associated(MPIMatrices%MPIDataR4_3D))               &
                    allocate(MPIMatrices%MPIDataR4_3D(1:dims(1),              &
                                                      1:dims(2),              &
                                                      1:dims(3)))

                call HDF5ReadData(inputFileHDF5,                              &
                                  GroupName = GroupName,                      &
                                  Name      = adjustl(trim(obj_name)),        &
                                  Array3D   = MPIMatrices%MPIDataR4_3D,       &
                                  STAT      = STAT_CALL)
if29 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "function ReadHDFDataSet, program DCCWorker, error calling HDF5ReadData, ERR11"
                end if if29

            else if (DataType .EQ. H5T_NATIVE_DOUBLE) then if14
                if (.NOT. associated(MPIMatrices%MPIDataR8_3D))               &
                    allocate(MPIMatrices%MPIDataR8_3D(1:dims(1),              &
                                                      1:dims(2),              &
                                                      1:dims(3)))


                call HDF5ReadData(inputFileHDF5,                              &
                                  GroupName = GroupName,                      &
                                  Name      = adjustl(trim(obj_name)),        &
                                  Array3D   = MPIMatrices%MPIDataR8_3D,       &
                                  STAT      = STAT_CALL)
if30 :          if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    stop "function ReadHDFDataSet, program DCCWorker, error calling HDF5ReadData, ERR12"
                end if if30
            end if if14
        end if if1

        readHDFDataSet = SUCCESS_

    end function readHDFDataSet

    !--------------------------------------------------------------------------

    integer function auxArrayAllocation(DomainSize,                           &
                                        Dimensions,                           &
                                        AuxMatrices,                          &
                                        LimitsArrayFactor,                    &
                                        MPIMatrices)

        integer, dimension(7), intent(IN) :: Dimensions
        integer, dimension(4), intent(IN) :: DomainSize
        type(T_AuxMatrices), pointer      :: AuxMatrices
        integer, intent(IN)               :: LimitsArrayFactor
        type(T_MPIMatrices), pointer      :: MPIMatrices

        if ((      associated(MPIMatrices%MPIDataI4_1D)) .AND.                    &
            (.NOT. associated(AuxMatrices%AuxDataI4_1D)))                         &
            allocate(AuxMatrices%AuxDataI4_1D(  Dimensions(1) + LimitsArrayFactor))

        if ((      associated(MPIMatrices%MPIDataI4_2D)) .AND.                    &
            (.NOT. associated(AuxMatrices%AuxDataI4_2D)))                         &
            allocate(AuxMatrices%AuxDataI4_2D(DomainSize(2) + LimitsArrayFactor,  &
                                              DomainSize(4) + LimitsArrayFactor))

        if ((      associated(MPIMatrices%MPIDataI4_3D)) .AND.                    &
            (.NOT. associated(AuxMatrices%AuxDataI4_3D)))                         &
            allocate(AuxMatrices%AuxDataI4_3D(DomainSize(2) + LimitsArrayFactor,  &
                                              DomainSize(4) + LimitsArrayFactor,  &
                                              Dimensions(3)))

        if ((      associated(MPIMatrices%MPIDataR4_1D)) .AND.                    &
            (.NOT. associated(AuxMatrices%AuxDataR4_1D)))                         &
            allocate(AuxMatrices%AuxDataR4_1D(Dimensions(1) + LimitsArrayFactor))

        if ((      associated(MPIMatrices%MPIDataR4_2D)) .AND.                    &
            (.NOT. associated(AuxMatrices%AuxDataR4_2D)))                         &
            allocate(AuxMatrices%AuxDataR4_2D(DomainSize(2) + LimitsArrayFactor,  &
                                              DomainSize(4) + LimitsArrayFactor))

        if ((      associated(MPIMatrices%MPIDataR4_3D)) .AND.                    &
            (.NOT. associated(AuxMatrices%AuxDataR4_3D)))                         &
            allocate(AuxMatrices%AuxDataR4_3D(DomainSize(2) + LimitsArrayFactor,  &
                                              DomainSize(4) + LimitsArrayFactor,  &
                                              Dimensions(3)))

        if ((      associated(MPIMatrices%MPIDataR8_1D)) .AND.                    &
            (.NOT. associated(AuxMatrices%AuxDataR8_1D)))                         &
            allocate(AuxMatrices%AuxDataR8_1D(Dimensions(1) + LimitsArrayFactor))

        if ((      associated(MPIMatrices%MPIDataR8_2D)) .AND.                    &
            (.NOT. associated(AuxMatrices%AuxDataR8_2D)))                         &
            allocate(AuxMatrices%AuxDataR8_2D(DomainSize(2) + LimitsArrayFactor,  &
                                              DomainSize(4) + LimitsArrayFactor))

        if ((      associated(MPIMatrices%MPIDataR8_3D)) .AND.                    &
            (.NOT. associated(AuxMatrices%AuxDataR8_3D)))                         &
            allocate(AuxMatrices%AuxDataR8_3D(DomainSize(2) + LimitsArrayFactor,  &
                                              DomainSize(4) + LimitsArrayFactor,  &
                                              Dimensions(3)))

        auxArrayAllocation = SUCCESS_

    end function auxArrayAllocation

    !--------------------------------------------------------------------------

    subroutine getIndexMPIHDF(WindowFrame,                                    &
                              LimitsArray,                                    &
                              GroupName,                                      &
                              ILB_MPI, IUB_MPI,                               &
                              JLB_MPI, JUB_MPI,                               &
                              KLB_MPI, KUB_MPI)

        integer, dimension(4),intent(IN)  :: WindowFrame
        integer, dimension(6),intent(IN)  :: LimitsArray
        character(*), intent(IN)          :: GroupName
        integer,intent(out)               :: ILB_MPI, IUB_MPI
        integer,intent(out)               :: JLB_MPI, JUB_MPI
        integer,intent(out)               :: KLB_MPI, KUB_MPI

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

    end subroutine getIndexMPIHDF

    !------------------------------------------------------------------------

    function getDomainSize(inputFileHDF5, key)
        integer, dimension(4)       :: getDomainSize
        type (T_HDF5), pointer      :: inputFileHDF5
        integer, intent(IN)         :: key
        character(PathLength)       :: DecompositionGroup = NULL_STR
        integer, dimension(4)       :: DataVal1D          = NULL_INT

if1 :   if      (key .EQ. 1) then
            DecompositionGroup = "/Grid/Decomposition/Global/"
        else if (key .EQ. 2) then if1
            DecompositionGroup = "/Grid/Decomposition/InnerMapping/"
        else if (key .EQ. 3) then if1
            DecompositionGroup = "/Grid/Decomposition/Mapping/"
        endif if1

        DataVal1D = GetMappingValues(inputFileHDF5,                           &
                                     GroupName = DecompositionGroup)

        getDomainSize = DataVal1D

    end function getDomainSize

    !------------------------------------------------------------------------

    function getMappingValues(inputFileHDF5, GroupName)

        integer, dimension(4)              :: getMappingValues
        type (T_HDF5), pointer             :: inputFileHDF5
        character(PathLength), intent(IN)  :: GroupName
        character(PathLength)              :: obj_name
        integer, dimension(:), pointer     :: DataVal1D
        integer, dimension(4)              :: DataVal1D_aux
        integer(HID_T)                     :: GroupType
        integer                            :: STAT_CALL   = NULL_INT

        obj_name = trim(adjustl("ILB_IUB_JLB_JUB"))

        allocate(DataVal1D(4))

        !Gets information about the group
        call GetHDF5ObjectInfo(inputFileHDF5,                                 &
                               FatherGroupName = adjustl(trim(GroupName)),    &
                               GroupPosition   = 1,                           &
                               GroupName       = adjustl(trim(obj_name)),     &
                               GroupType       = GroupType,                   &
                               STAT = STAT_CALL)
if7 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function getMappingValues, program DCCWorker, error calling GetHDF5ObjectInfo, ERR01"
        end if if7

if1 :   if (GroupType == H5G_DATASET_F) then
            call HDF5SetLimits(inputFileHDF5,                                 &
                               ILB = 1, IUB = 4,                              &
                               STAT = STAT_CALL)
if3 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function getMappingValues, program DCCWorker, error calling HDF5SetLimits, ERR02"
            end if if3

            call HDF5ReadData(inputFileHDF5,                                  &
                              GroupName = adjustl(trim(GroupName)),           &
                              Name      = adjustl(trim(obj_name)),            &
                              Array1D   = DataVal1D,                          &
                              STAT      = STAT_CALL)
if4 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function getMappingValues, program DCCWorker, error calling HDF5ReadData, ERR03"
            end if if4

        elseif (GroupType .EQ. H5G_GROUP_F) then if1
            STAT_CALL = -1 * H5G_GROUP_F
if2 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function getMappingValues, program DCCWorker, error calling HDF5ReadData, ERR04"
            end if if2
        end if if1

        DataVal1D_aux = DataVal1D

        deallocate(DataVal1D)

        getMappingValues = DataVal1D_aux

    end function getMappingValues

    !------------------------------------------------------------------------

    integer function renameOutputFile(consolidatedFile, consolidatedFileTmp)
        character(PathLength),      intent(IN)  :: consolidatedFile
        character(PathLength),      intent(IN)  :: consolidatedFileTmp
        character(StringLength)                 :: OnlyFileName
        integer                                 :: STAT_CALL    = NULL_INT
        integer                                 :: iFile        = NULL_INT
        integer                                 :: iFile2       = NULL_INT

        call UnitsManager(iFile, OPEN_FILE, STAT = STAT_CALL)
if1 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function renameOutputFile, error calling UnitsManager, ERR01"
        end if if1

        open(unit   = iFile,                                                            &
             file   = trim(adjustl(consolidatedFile)),                                  &
             status = 'old',                                                            &
             iostat = STAT_CALL)

if2 :   if (STAT_CALL .EQ. SUCCESS_) then
            call UnitsManager(iFile, CLOSE_FILE, STAT = STAT_CALL)
if4 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function renameOutputFile, error calling UnitsManager, ERR02"
            end if if4

            call UnitsManager(iFile2, OPEN_FILE, STAT = STAT_CALL)
if3 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "function renameOutputFile, error calling UnitsManager, ERR03"
            end if if3

            open(unit   = iFile2,                                                       &
                 file   = trim(adjustl(consolidatedFileTmp)),                           &
                 status = 'old',                                                        &
                 iostat = STAT_CALL)
            if (STAT_CALL .EQ. SUCCESS_)                                                &
                close(unit   = iFile2,                                                  &
                      status = 'delete')
        else if2
            !Linux
            STAT_CALL = SYSTEM("mv "                                                    &
                               // trim(adjustl(consolidatedFileTmp)) // " "             &
                               // trim(adjustl(consolidatedFile)))

if5 :       if (STAT_CALL .NE. SUCCESS_) then
                !Windows
                STAT_CALL = ReturnOnlyFileName(consolidatedFile, OnlyFileName)
                
if7:            if (STAT_CALL /= SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    print*, "Error in extracting the filename  from",                   &
                            trim(adjustl(consolidatedFile))
                    stop "function renameOutputFile, error calling ReturnOnlyFileName, ERR04"
                    
                end if if7                
                
                STAT_CALL = SYSTEM("ren "                                               &
                                   // trim(adjustl(consolidatedFileTmp)) // " "         &
                                   // trim(adjustl(OnlyFileName)))
if6 :           if (STAT_CALL .NE. SUCCESS_) then
                    print*, "STAT_CALL = ", STAT_CALL
                    print*, "Error renaming file",                                      &
                            trim(adjustl(consolidatedFileTmp)), " to ",                 &
                            trim(adjustl(consolidatedFile))
                    stop "function renameOutputFile, error calling SYSTEM, ERR05"
                end if if6
            end if if5
        end if if2

        renameOutputFile = SUCCESS_

    end function renameOutputFile

    !------------------------------------------------------------------------
    
    integer function ReturnOnlyFileName(FileNamePlusPath, OnlyFileName)
        character(*),   intent(IN)  :: FileNamePlusPath
        character(*),   intent(OUT) :: OnlyFileName
        integer                     :: STAT_CALL    = NULL_INT
        integer                     :: n            = NULL_INT
        integer                     :: j            = NULL_INT
        integer                     :: i            = NULL_INT        


        n = len_trim(FileNamePlusPath)
        
        j = 1
        
        do i=n,1,-1
            if (FileNamePlusPath(i:i)=='/' .or. FileNamePlusPath(i:i)=='\') then
                j = i+1
                exit
            endif
        enddo
        
if1:    if (j < n) then
        
            OnlyFileName = FileNamePlusPath(j:n)
            
            STAT_CALL = SUCCESS_
            
        else if1
        
            write(*,*) "STAT_CALL = ", STAT_CALL
            stop "function ReturnOnlyFileName, DDCWorker, ERR06"

        end if if1
        
        ReturnOnlyFileName = STAT_CALL


    end function ReturnOnlyFileName

    !------------------------------------------------------------------------
    

    integer function recvPoisonPill(Me)
        type (T_DDC), pointer       :: Me
        integer                     :: STAT_CALL    = NULL_INT
        integer                     :: msg          = NULL_INT
        integer                     :: STATUS(MPI_STATUS_SIZE)

        call MPI_RECV(msg,                                                      &
                      1,                                                        &
                      MPI_INTEGER,                                              &
                      getParser_id(Me),                                         &
                      getMsgPoisonPillTag(),                                    &
                      MPI_COMM_WORLD,                                           &
                      STATUS,                                                   &
                      STAT_CALL)
if1 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function recvPoisonPill, program DCCWorker, error calling MPI_RECV, ERR01"
        end if if1

        call endPrg(Me)

        recvPoisonPill = SUCCESS_

    end function recvPoisonPill

    !---------------------------------------------------------------------------

    subroutine endPrg(Me)
        type (T_DDC), pointer       :: Me
        integer                     :: STAT_CALL  = UNKNOWN_

        STAT_CALL = killDDC(Me)
if2 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "subroutine endPrg, program DCCWorker, error calling killDDC, ERR03"
        end if if2

        call EXIT(SUCCESS_)

    end subroutine endPrg

    !---------------------------------------------------------------------------

    integer function recvParser_id(Me, parser_id)
        type (T_DDC), pointer       :: Me
        integer, intent(IN)         :: parser_id
        integer                     :: STAT_CALL    = NULL_INT
        integer                     :: msg          = NULL_INT
        integer                     :: STATUS(MPI_STATUS_SIZE)

        call MPI_RECV(msg,                                                      &
                      1,                                                        &
                      MPI_INTEGER,                                              &
                      parser_id,                                                &
                      getMsgRankTag(),                                          &
                      MPI_COMM_WORLD,                                           &
                      STATUS,                                                   &
                      STAT_CALL)
if1 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function recvParser_id, program DCCWorker, error calling MPI_RECV, ERR01"
        end if if1

        STAT_CALL = setParser_id(Me, parser_id)
if6 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function recvParser_id, program DCCWorker, error calling setParser_id, ERR05"
        end if if6

        STAT_CALL = sendIdle(Me)
if5 :   if (STAT_CALL .NE. SUCCESS_) then
            print*, "STAT_CALL = ", STAT_CALL
            stop "function recvParser_id, program DCCWorker, error calling sendIdle, ERR02"
        end if if5

        recvParser_id = SUCCESS_

    end function recvParser_id

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    integer function killDDC(Me)
        type (T_DDC), pointer                       :: Me
        integer                     :: STAT_CALL = NULL_INT

        call DeallocateInstance(Me)

!        STAT_CALL = barrier()
!if5 :   if (STAT_CALL .NE. SUCCESS_) then
!            print*, "STAT_CALL = ", STAT_CALL
!            stop "function killDDC, program DCCWorker, error calling barrier, ERR05"
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
        type (T_DDC), pointer       :: Me

        deallocate(Me%task)
        deallocate(Me%inputFile)
        deallocate(Me)
        nullify   (Me)

    end subroutine deallocateInstance

    !------------------------------------------------------------------------

    integer function deallocateTaskMPIMatrices(task)
        type(T_Task), pointer       :: task

        if (associated(task%MPIMatrices))                                     &
            call deallocateMPIMatrices(task%MPIMatrices)

        deallocate(task%MPIMatrices)
        nullify   (task%MPIMatrices)

        deallocateTaskMPIMatrices = SUCCESS_

    end function deallocateTaskMPIMatrices

    !------------------------------------------------------------------------

    integer function deallocateTaskAuxMatrices(task)
        type(T_Task), pointer       :: task

        if (associated(task%AuxMatrices))                                     &
            call deallocateAuxMatrices(task%AuxMatrices)

        deallocate(task%AuxMatrices)
        nullify   (task%AuxMatrices)

        deallocateTaskAuxMatrices = SUCCESS_

    end function deallocateTaskAuxMatrices

    !------------------------------------------------------------------------

    recursive subroutine deallocateInputFile(inputFile)
        type(T_InputFile), pointer  :: inputFile
        integer                     :: STAT_CALL    = NULL_INT

        if (associated(inputFile%Next))                                       &
            call deallocateInputFile(inputFile%Next)

if1 :   if (associated(inputFile%inputFileHDF5)) then
            STAT_CALL= KillHDF5(inputFile%inputFileHDF5)
if2 :       if (STAT_CALL .NE. SUCCESS_) then
                print*, "STAT_CALL = ", STAT_CALL
                stop "subroutine deallocateInputFile, program DCCWorker, error calling KillHDF5, ERR02"
            end if if2
        end if if1

        deallocate(inputFile)
        nullify   (inputFile)

    end subroutine deallocateInputFile

    !------------------------------------------------------------------------

    subroutine deallocateAuxMatrices(AuxMatrices)
        type (T_AuxMatrices), pointer         :: AuxMatrices

if1 :   if (associated(AuxMatrices%AuxDataR4_1D)) then
            deallocate(AuxMatrices%AuxDataR4_1D)
            nullify   (AuxMatrices%AuxDataR4_1D)
        end if if1

if2 :   if (associated(AuxMatrices%AuxDataR4_2D)) then
            deallocate(AuxMatrices%AuxDataR4_2D)
            nullify   (AuxMatrices%AuxDataR4_2D)
        end if if2

if3 :   if (associated(AuxMatrices%AuxDataR4_3D)) then
            deallocate(AuxMatrices%AuxDataR4_3D)
            nullify   (AuxMatrices%AuxDataR4_3D)
        end if if3

if4 :   if (associated(AuxMatrices%AuxDataR8_1D)) then
            deallocate(AuxMatrices%AuxDataR8_1D)
            nullify   (AuxMatrices%AuxDataR8_1D)
        end if if4

if5 :   if (associated(AuxMatrices%AuxDataR8_2D)) then
            deallocate(AuxMatrices%AuxDataR8_2D)
            nullify   (AuxMatrices%AuxDataR8_2D)
        end if if5

if6 :   if (associated(AuxMatrices%AuxDataR8_3D)) then
            deallocate(AuxMatrices%AuxDataR8_3D)
            nullify   (AuxMatrices%AuxDataR8_3D)
        end if if6

if7 :   if (associated(AuxMatrices%AuxDataI4_1D)) then
            deallocate(AuxMatrices%AuxDataI4_1D)
            nullify   (AuxMatrices%AuxDataI4_1D)
        end if if7

if8 :   if (associated(AuxMatrices%AuxDataI4_2D)) then
            deallocate(AuxMatrices%AuxDataI4_2D)
            nullify   (AuxMatrices%AuxDataI4_2D)
        end if if8

if9 :   if (associated(AuxMatrices%AuxDataI4_3D)) then
            deallocate(AuxMatrices%AuxDataI4_3D)
            nullify   (AuxMatrices%AuxDataI4_3D)
        end if if9

    end subroutine deallocateAuxMatrices

    !------------------------------------------------------------------------

    recursive subroutine deallocateMPIMatrices(MPIMatrices)
        type (T_MPIMatrices), pointer         :: MPIMatrices

if1 :   if (associated(MPIMatrices%MPIDataR4_1D)) then
            deallocate(MPIMatrices%MPIDataR4_1D)
            nullify   (MPIMatrices%MPIDataR4_1D)
        end if if1

if2 :   if (associated(MPIMatrices%MPIDataR4_2D)) then
            deallocate(MPIMatrices%MPIDataR4_2D)
            nullify   (MPIMatrices%MPIDataR4_2D)
        end if if2

if3 :   if (associated(MPIMatrices%MPIDataR4_3D)) then
            deallocate(MPIMatrices%MPIDataR4_3D)
            nullify   (MPIMatrices%MPIDataR4_3D)
        end if if3

if4 :   if (associated(MPIMatrices%MPIDataR8_1D)) then
            deallocate(MPIMatrices%MPIDataR8_1D)
            nullify   (MPIMatrices%MPIDataR8_1D)
        end if if4

if5 :   if (associated(MPIMatrices%MPIDataR8_2D)) then
            deallocate(MPIMatrices%MPIDataR8_2D)
            nullify   (MPIMatrices%MPIDataR8_2D)
        end if if5

if6 :   if (associated(MPIMatrices%MPIDataR8_3D)) then
            deallocate(MPIMatrices%MPIDataR8_3D)
            nullify   (MPIMatrices%MPIDataR8_3D)
        end if if6

if7 :   if (associated(MPIMatrices%MPIDataI4_1D)) then
            deallocate(MPIMatrices%MPIDataI4_1D)
            nullify   (MPIMatrices%MPIDataI4_1D)
        end if if7

if8 :   if (associated(MPIMatrices%MPIDataI4_2D)) then
            deallocate(MPIMatrices%MPIDataI4_2D)
            nullify   (MPIMatrices%MPIDataI4_2D)
        end if if8

if9 :   if (associated(MPIMatrices%MPIDataI4_3D)) then
            deallocate(MPIMatrices%MPIDataI4_3D)
            nullify   (MPIMatrices%MPIDataI4_3D)
        end if if9

    end subroutine deallocateMPIMatrices

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

end program DDCWorker

