!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : GlueWW3_OBC
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to serve as GlueWW3_OBC to create new modules
!
!------------------------------------------------------------------------------


Module ModuleGlueWW3_OBC

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructGlueWW3_OBC
    private ::      AllocateInstance

    !Selector
    public  :: GetGlueWW3_OBCPointer
    public  :: GetGlueWW3_OBCInteger
    public  :: UnGetGlueWW3_OBC
                     
    
    !Modifier
    public  :: ModifyGlueWW3_OBC

    !Destructor
    public  :: KillGlueWW3_OBC                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjGlueWW3_OBC 
    
    !Parameters----------------------------------------------------------------
    integer,    parameter                   :: NmaxInst = 10000
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetGlueWW3_OBC3D_I
    private :: UnGetGlueWW3_OBC3D_R8
    interface  UnGetGlueWW3_OBC
        module procedure UnGetGlueWW3_OBC3D_I
        module procedure UnGetGlueWW3_OBC3D_R8
    end interface  UnGetGlueWW3_OBC

    !Types---------------------------------------------------------------------

    type T_File     
        character(LEN=10)                   :: VERBPTBC = 'III  1.03 '
        character(LEN=32)                   :: IDSTRBC  = 'WAVEWATCH III BOUNDARY DATA FILE'  
        integer, dimension(2)               :: TIME1
        integer                             :: NPTS  
        integer,    pointer                 :: NK, NTH, NBO
        integer                             :: NSPEC
        real,       pointer                 :: XFR, FR1
        real,       pointer, dimension(:)   :: TH, XBPO, YBPO
        integer,    pointer, dimension(:,:) :: IPBPO
        real,       pointer, dimension(:,:) :: RDBPO, ABPOS
        integer                             :: Ninstants
        type(T_Time), pointer, dimension(:) :: DateTime
        integer                             :: Unit, UnitStored
        character(len=PathLength)           :: Name, NameStored
        type (T_File)  , pointer            :: Next, Previous
    end type T_File
    
    
    private :: T_GlueWW3_OBC
    type       T_GlueWW3_OBC
        integer                             :: InstanceID
        type (T_Size3D)                     :: Size, WorkSize
        type (T_File)       , pointer       :: FirstFile, LastFile
        integer                             :: FilesNumber
        character(len=PathLength)           :: FileOutput
        integer                             :: UnitOutput
        type (T_GlueWW3_OBC), pointer       :: Next
    end type  T_GlueWW3_OBC

    !Global Module Variables
    type (T_GlueWW3_OBC), pointer           :: FirstObjGlueWW3_OBC
    type (T_GlueWW3_OBC), pointer           :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructGlueWW3_OBC(ObjGlueWW3_OBCID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjGlueWW3_OBCID 
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mGlueWW3_OBC_)) then
            nullify (FirstObjGlueWW3_OBC)
            call RegisterModule (mGlueWW3_OBC_) 
        endif

        call Ready(ObjGlueWW3_OBCID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            call ReadKeywords

            !Returns ID
            ObjGlueWW3_OBCID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleGlueWW3_OBC - ConstructGlueWW3_OBC - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructGlueWW3_OBC
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_GlueWW3_OBC), pointer                         :: NewObjGlueWW3_OBC
        type (T_GlueWW3_OBC), pointer                         :: PreviousObjGlueWW3_OBC


        !Allocates new instance
        allocate (NewObjGlueWW3_OBC)
        nullify  (NewObjGlueWW3_OBC%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjGlueWW3_OBC)) then
            FirstObjGlueWW3_OBC         => NewObjGlueWW3_OBC
            Me                    => NewObjGlueWW3_OBC
        else
            PreviousObjGlueWW3_OBC      => FirstObjGlueWW3_OBC
            Me                    => FirstObjGlueWW3_OBC%Next
            do while (associated(Me))
                PreviousObjGlueWW3_OBC  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjGlueWW3_OBC
            PreviousObjGlueWW3_OBC%Next => NewObjGlueWW3_OBC
        endif

        Me%InstanceID = RegisterNewInstance (mGlueWW3_OBC_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        type (T_File), pointer                      :: NewFile
        integer                                     :: STAT_CALL, ClientNumber, iflag
        integer                                     :: Line, FirstLine, LastLine
        integer                                     :: ObjEnterData = 0
        integer                                     :: FromFile, FromBlock
        character(LEN = StringLength), parameter    :: input_files_begin   = '<begin_input_files>'
        character(LEN = StringLength), parameter    :: input_files_end     = '<end_input_files>'
        logical                                     :: BlockFound, ReadError
        character(LEN = PathLength  )               :: InputFile
        
        !Begin-----------------------------------------------------------------

        call ConstructEnterData (ObjEnterData, FileName = 'GlueWW3_OBC.dat', STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ModuleGlueWW3_OBC - ERR10'
        
        call GetExtractType     (FromFile = FromFile)        
        
        call GetData(Me%FileOutput,                                                     &
                     ObjEnterData,iflag,                                                &
                     SearchType     = FromFile,                                         &
                     keyword        = 'FILE_OUTPUT',                                    &
                     ClientModule   = 'MainGlueWW3_OBC',                                &
                     STAT           = STAT_CALL)              
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ModuleGlueWW3_OBC - ERR20'
        
        call UnitsManager(UNIT = Me%UnitOutput, OPENCLOSE = OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ModuleGlueWW3_OBC - ERR30'
    
        open (UNIT=Me%UnitOutput, FILE=Me%FileOutput, FORM='UNFORMATTED', STATUS='UNKNOWN', IOSTAT=STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ModuleGlueWW3_OBC - ERR40'
        
        call GetExtractType     (FromBlock = FromBlock)
        
        call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                         &
                                   input_files_begin, input_files_end,                  &
                                   BlockFound = BlockFound,                             &
                                   FirstLine = FirstLine, LastLine = LastLine,          &
                                   STAT = STAT_CALL)

IS:     if(STAT_CALL .EQ. SUCCESS_) then

            !The block is found to exist before when reading depth
BF:         if (BlockFound) then            

                do line = FirstLine + 1, LastLine - 1

                    call GetData(InputFile, EnterDataID = ObjEnterData, flag = iflag,&
                                 Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ModuleGlueWW3_OBC - ERR50'
                    
                    nullify (NewFile)
                    allocate(NewFile)
                    
                    NewFile%Name = InputFile
                    
                    call ReadWW3_OBC_Header(NewFile, ReadError)
                    
                    if (ReadError) then
                        stop 'ReadKeywords - ModuleGlueWW3_OBC - ERR60'
                    else
                        call Add_File   (NewFile)
                        
                        call ReadWW3_OBC(NewFile)
                    endif        
                    
                enddo                                    
                
            else  BF
                stop 'ReadKeywords - ModuleGlueWW3_OBC - ERR70'
            endif BF

            call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ModuleGlueWW3_OBC - ERR80'

        else   IS

            stop 'ReadKeywords - ModuleGlueWW3_OBC - ERR90'

        end if IS

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ModuleGlueWW3_OBC - ERR100'

    end subroutine ReadKeywords


    !--------------------------------------------------------------------------
    ! This subroutine adds a new file to the File List  

    subroutine Add_File(NewFile)

        !Arguments-------------------------------------------------------------
        type(T_File),           pointer     :: NewFile

        !----------------------------------------------------------------------

        ! Add to the List a new File
        if (.not.associated(Me%FirstFile)) then
            Me%FilesNumber       = 1
            Me%FirstFile         => NewFile
            Me%LastFile          => NewFile
        else
            NewFile%Previous     => Me%LastFile
            Me%LastFile%Next     => NewFile
            Me%LastFile          => NewFile
            Me%FilesNumber       =  Me%FilesNumber + 1
        end if 

        !----------------------------------------------------------------------

    end subroutine Add_File 

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------    
    subroutine ReadWW3_OBC_Header(File, ReadError)

        !Arguments-------------------------------------------------------------            
        type (T_File), pointer :: File
        logical                :: ReadError
        
        !Local-----------------------------------------------------------------            
        integer                :: i, j, STAT_CALL, IERR

        !Begin-----------------------------------------------------------------             
        
        call UnitsManager(UNIT = File%UnitStored, OPENCLOSE = OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadWW3_OBC_Header - ERR10'                        
        
        write(File%NameStored,*) File%UnitStored
        File%NameStored = trim(adjustl(File%NameStored))
        
        open (UNIT=File%UnitStored, FILE=File%NameStored, FORM='UNFORMATTED',           &
              STATUS='UNKNOWN', ACCESS = 'STREAM', IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadWW3_OBC_Header - ERR20' 

        call UnitsManager(UNIT = File%Unit, OPENCLOSE = OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadWW3_OBC_Header - ERR30'                        
    
        open (UNIT=File%Unit, FILE=File%Name, FORM='UNFORMATTED', STATUS='OLD', IOSTAT=IERR)
        
        if (IERR == SUCCESS_) then
        
            allocate (File%TH(1))
            allocate (File%NK, File%NTH, File%XFR, File%FR1, File%NBO)
            
            read (File%Unit) File%IDSTRBC, File%VERBPTBC, File%NK, File%NTH, File%XFR,  &
                             File%FR1, File%TH, File%NBO
            
            write(*,*) File%IDSTRBC, trim(File%Name)
            write(*,*) 'Version ', File%VERBPTBC
            
            File%NSPEC  = File%NK * File%NTH              
            
            allocate (File%XBPO (1:File%NBO    ), File%YBPO (1:File%NBO    ))        
            allocate (File%IPBPO(1:File%NBO,1:4), File%RDBPO(1:File%NBO,1:4))
    !
            read (File%Unit)                                                                &
                 (File%XBPO(I),I=1,File%NBO),                                               &
                 (File%YBPO(I),I=1,File%NBO),                                               &
                 ((File%IPBPO(I,J),I=1,File%NBO),J=1,4),                                    &
                 ((File%RDBPO(I,J),I=1,File%NBO),J=1,4) 
                 
            ReadError = .false.

        else
        
            ReadError = .true.
        
        endif        

    end subroutine ReadWW3_OBC_Header
    
    !--------------------------------------------------------------------------    
    subroutine ReadWW3_OBC(File)

        !Arguments-------------------------------------------------------------            
        type (T_File), pointer :: File
        
        !Local----------------------------------------------------------------- 
        type (T_Time), dimension(:), pointer :: AuxTime
        integer                              :: IS,ISOUT, N, STAT_CALL          
        real                                 :: Year, Month, Day, Hour, Minutes, Seconds
        real(8)                              :: AuxT1, AuxT2
        logical                              :: EndRead
        
        !Begin-----------------------------------------------------------------  
        
        allocate(AuxTime(1:NmaxInst))
        
        do N = 1, NmaxInst
            call null_time(AuxTime(N))
        enddo            
 
        EndRead = .false.
        
        N = 0
        
        Do While(.not. EndRead)
        
            read (File%Unit, IOSTAT = STAT_CALL) File%TIME1, File%NPTS  

            if (STAT_CALL == SUCCESS_) then

                N = N + 1
                
                if (N > NmaxInst) then
                
                    stop 'ReadWW3_OBC - ERR10'                        
                    
                endif

                AuxT1 = File%TIME1(1)
                AuxT2 = File%TIME1(2)

                Year    = int(AuxT1 / 1e4)
                Month   = int(AuxT1 -Year * 1e4) / 100
                Day     = AuxT1 -Year * 1e4 - Month * 100
                
                Hour    = int(AuxT2 / 1e4)
                Minutes = int(AuxT2 -Hour * 1e4) / 100
                Seconds = AuxT2 -Hour * 1e4 - Minutes * 100
                
                write (File%UnitStored) Year, Month, Day, Hour, Minutes, Seconds
                
                call SetDate(AuxTime(N), Year, Month, Day, Hour, Minutes, Seconds)
                
                if (.not.associated(File%ABPOS)) then
                    allocate(File%ABPOS(1:File%NSPEC,1:File%NPTS))
                endif
                
                do ISOUT=1, File%NPTS
                    read  (File%Unit      ) (File%ABPOS(IS,ISOUT),IS=1,File%NSPEC)
                    write (File%UnitStored) (File%ABPOS(IS,ISOUT),IS=1,File%NSPEC)
                enddo
                
                EndRead = .false. 
            else
                EndRead = .true.     
            endif
        enddo            
        
        File%Ninstants = N 
        
        allocate(File%DateTime(1:File%Ninstants))        
        File%DateTime(1:File%Ninstants) = AuxTime(1:File%Ninstants)
        
        deallocate(AuxTime)
        

    end subroutine ReadWW3_OBC


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    
    !--------------------------------------------------------------------------
    subroutine GetGlueWW3_OBCPointer (ObjGlueWW3_OBCID, Matrix, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjGlueWW3_OBCID
        real(8), dimension(:, :, :),  pointer           :: Matrix
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjGlueWW3_OBCID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mGlueWW3_OBC_, Me%InstanceID)

            !Matrix => Me%Matrix

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetGlueWW3_OBCPointer
    
    !--------------------------------------------------------------------------
    
    subroutine GetGlueWW3_OBCInteger (ObjGlueWW3_OBCID, Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjGlueWW3_OBCID
        real                                            :: Int
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjGlueWW3_OBCID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Int = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetGlueWW3_OBCInteger

    !--------------------------------------------------------------------------

    subroutine UnGetGlueWW3_OBC3D_I(ObjGlueWW3_OBCID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjGlueWW3_OBCID
        integer, dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjGlueWW3_OBCID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mGlueWW3_OBC_, Me%InstanceID, "UnGetGlueWW3_OBC3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetGlueWW3_OBC3D_I

    !--------------------------------------------------------------------------

    subroutine UnGetGlueWW3_OBC3D_R8(ObjGlueWW3_OBCID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjGlueWW3_OBCID
        real(8), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjGlueWW3_OBCID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mGlueWW3_OBC_, Me%InstanceID,  "UnGetGlueWW3_OBC3D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetGlueWW3_OBC3D_R8

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyGlueWW3_OBC(ObjGlueWW3_OBCID, StartTime, EndTime, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjGlueWW3_OBCID
        type (T_Time)                               :: StartTime, EndTime
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjGlueWW3_OBCID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
            
            call WriteFinalFile(StartTime, EndTime)
                
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyGlueWW3_OBC
    
    !----------------------------------------------------------------------------
    
    subroutine WriteFinalFile(StartTime, EndTime)

        !Arguments-------------------------------------------------------------
        type (T_Time)                               :: StartTime, EndTime
        !Local-----------------------------------------------------------------
        type (T_File), pointer                      :: CurrentFile, FirstFile
        real                                        :: Year,  Month,  Day,  Hour,  Minutes,  Seconds 
        real                                        :: YearX, MonthX, DayX, HourX, MinutesX, SecondsX
        integer, dimension(2)                       :: DateX
        integer                                     :: n, nStart, p, ISOUT, IS, i, j, Precision
        logical                                     :: NotLastDate
        real(8)                                     :: TimeD, DateD
        !----------------------------------------------------------------------
        
        if (Me%FirstFile%DateTime(1) >  StartTime) then
            write(*,*) 'Start date older than the oldest date available in files'
            stop 'WriteFinalFile - ModuleGlueWW3_OBC - ERR100' 
        endif

        
        if (Me%LastFile%DateTime(Me%LastFile%Ninstants) <  EndTime) then
            write(*,*) 'End date newer than the newest date available in files'
            stop 'WriteFinalFile - ModuleGlueWW3_OBC - ERR200' 
        endif

        
        CurrentFile => Me%LastFile

       !Find first file and first instant
        do while (associated(CurrentFile))
            if (CurrentFile%DateTime(1) <= StartTime) then
                FirstFile => CurrentFile
                do n=1,FirstFile%Ninstants-1
                    
                    if (CurrentFile%DateTime(n  )<= StartTime .and.                     &
                        CurrentFile%DateTime(n+1)>  StartTime) then
                        nStart = n
                        exit
                    endif                        
                
                enddo
                exit
            endif
            CurrentFile => CurrentFile%Previous
        enddo 
        
        if (.not. associated(FirstFile))  then
            stop 'WriteFinalFile - ModuleGlueWW3_OBC - ERR10'          
        endif
        
        !Write header
        write (Me%UnitOutput) FirstFile%IDSTRBC, FirstFile%VERBPTBC, FirstFile%NK,      &
                              FirstFile%NTH, FirstFile%XFR,                             &
                              FirstFile%FR1, FirstFile%TH, FirstFile%NBO
        
        write(*,*) FirstFile%IDSTRBC, trim(Me%FileOutput)
        write(*,*) 'Version ', FirstFile%VERBPTBC
        
        FirstFile%NSPEC  = FirstFile%NK * FirstFile%NTH              
        
!
        write (Me%UnitOutput)                                                           &
             (FirstFile%XBPO(I),I=1,FirstFile%NBO),                                     &
             (FirstFile%YBPO(I),I=1,FirstFile%NBO),                                     &
             ((FirstFile%IPBPO(I,J),I=1,FirstFile%NBO),J=1,4),                          &
             ((FirstFile%RDBPO(I,J),I=1,FirstFile%NBO),J=1,4)         

        n           =  nStart 
        NotLastDate = .true. 
        
        CurrentFile => FirstFile
        
        do while (associated(CurrentFile) .and. NotLastDate)
        
            Precision = sizeof(Year)
        
            p = 1 + (n-1) * (6 + CurrentFile%NSPEC * CurrentFile%NPTS) * Precision
            
            read (CurrentFile%UnitStored, POS = p) Year, Month, Day, Hour, Minutes, Seconds
            
            p = p + 6 * Precision
            
            call ExtractDate(CurrentFile%DateTime(n), YearX, MonthX, DayX, HourX, MinutesX, SecondsX)
            
            if (YearX /= Year .or. MonthX   /= Month   .or. DayX     /= Day .or.        &
                HourX /= Hour .or. MinutesX /= Minutes .or. SecondsX /= Seconds) then
                stop 'WriteFinalFile - ModuleGlueWW3_OBC - ERR20'
            endif

            !Date
            DateD    = dble(Year)*1e4 + dble(Month)  *100 + dble(Day)
            DateX(1) = int(DateD)
            !Time
            TimeD    = dble(Hour)*1e4 + dble(Minutes)*100 + dble(Seconds)
            DateX(2) = int(TimeD)
            
            write(Me%UnitOutPut) DateX, CurrentFile%NPTS

            do ISOUT=1, CurrentFile%NPTS
                read(CurrentFile%UnitStored, POS = p) (CurrentFile%ABPOS(IS,ISOUT),IS=1,CurrentFile%NSPEC)
                p = p + CurrentFile%NSPEC * Precision
                
                write(Me%UnitOutput) (CurrentFile%ABPOS(IS,ISOUT),IS=1,CurrentFile%NSPEC)
                
            enddo
            
            if (CurrentFile%DateTime(n) >= EndTime) then
                NotLastDate = .false. 
            else
                if (associated(CurrentFile%Next)) then
                    if (CurrentFile%DateTime(n+1) >= CurrentFile%Next%DateTime(1)) then
                        CurrentFile => CurrentFile%Next
                        n = 1
                    else
                        n = n + 1                
                        if (n > CurrentFile%Ninstants) stop 'WriteFinalFile - ModuleGlueWW3_OBC - ERR30'
                    endif
                else
                    n = n + 1                
                    if (n > CurrentFile%Ninstants) stop 'WriteFinalFile - ModuleGlueWW3_OBC - ERR40'
                endif
            endif            
        enddo 

    end subroutine WriteFinalFile
    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillGlueWW3_OBC(ObjGlueWW3_OBCID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjGlueWW3_OBCID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjGlueWW3_OBCID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mGlueWW3_OBC_,  Me%InstanceID)

            if (nUsers == 0) then

                call KillVariables

                !Deallocates Instance
                call DeallocateInstance ()
                
                ObjGlueWW3_OBCID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillGlueWW3_OBC
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_GlueWW3_OBC), pointer          :: AuxObjGlueWW3_OBC
        type (T_GlueWW3_OBC), pointer          :: PreviousObjGlueWW3_OBC

        !Updates pointers
        if (Me%InstanceID == FirstObjGlueWW3_OBC%InstanceID) then
            FirstObjGlueWW3_OBC => FirstObjGlueWW3_OBC%Next
        else
            PreviousObjGlueWW3_OBC => FirstObjGlueWW3_OBC
            AuxObjGlueWW3_OBC      => FirstObjGlueWW3_OBC%Next
            do while (AuxObjGlueWW3_OBC%InstanceID /= Me%InstanceID)
                PreviousObjGlueWW3_OBC => AuxObjGlueWW3_OBC
                AuxObjGlueWW3_OBC      => AuxObjGlueWW3_OBC%Next
            enddo

            !Now update linked list
            PreviousObjGlueWW3_OBC%Next => AuxObjGlueWW3_OBC%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance
    

    !------------------------------------------------------------------------

    subroutine KillVariables()

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        type (T_File), pointer              :: CurrentFile
        integer                             :: STAT_CALL      

        !------------------------------------------------------------------------
    
        call UnitsManager(UNIT = Me%UnitOutput, OPENCLOSE = CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'KillVariables - ModuleGlueWW3_OBC - ERR10'    
        
        !Kill Files
        CurrentFile => Me%FirstFile
        do while (associated(CurrentFile))       
            
            call UnitsManager(UNIT = CurrentFile%UnitStored, STATUS = 'DELETE',         &
                              OPENCLOSE = CLOSE_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'KillVariables - ModuleGlueWW3_OBC - ERR20'    
            
            CurrentFile => CurrentFile%Next
        enddo
                
        !------------------------------------------------------------------------

    end subroutine KillVariables
        

  

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjGlueWW3_OBC_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjGlueWW3_OBC_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjGlueWW3_OBC_ID > 0) then
            call LocateObjGlueWW3_OBC (ObjGlueWW3_OBC_ID)
            ready_ = VerifyReadLock (mGlueWW3_OBC_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjGlueWW3_OBC (ObjGlueWW3_OBCID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjGlueWW3_OBCID

        !Local-----------------------------------------------------------------

        Me => FirstObjGlueWW3_OBC
        do while (associated (Me))
            if (Me%InstanceID == ObjGlueWW3_OBCID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleGlueWW3_OBC - LocateObjGlueWW3_OBC - ERR01'

    end subroutine LocateObjGlueWW3_OBC

    !--------------------------------------------------------------------------

end module ModuleGlueWW3_OBC









