!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : EnterData
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to Write/Read ASCII Files
!
!------------------------------------------------------------------------------
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License 
!version 2, as published by the Free Software Foundation.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!
!------------------------------------------------------------------------------

Module ModuleEnterData

    use ModuleGlobalData
    use ModuleTime               

    implicit none

    private

    !Subroutine Tree----------------------------------------------------------------

    !Constructor
    public  ::  ConstructEnterData
    private ::      AllocateInstance

    !Modifier
    public  :: RewindBuffer     !Sets buffer for new search
    public  :: RewindBlock      !Sets block for new search
    public  :: ReplaceFullBufferLine

    !Selector
    private ::      LookForLineInBlock
    private ::          KeywordsEqual !Logical Function
    private ::      ScanLine
    public  :: GetFileName
    public  :: GetData
    public  :: GetExtractType
    public  :: GetBufferSize
    public  :: GetFullBufferLine
    public  :: GetBlockSize
    public  :: GetKeywordFromLine
    public  :: ReadFileName
    private ::      CreateName
    public  :: ExtractBlockFromBuffer
    private ::      ExtractBlockFromBuffer1
    public  :: ExtractBlockFromBlock
    public  :: ExtractBlockFromBlockFromBlock


    public  :: GetOutPutTime
    private ::      OutPutTimeInternal
    public  :: GetOutPutTimeWindows
    private ::      CountNumberOfWindows
    private ::      ReadOutPutWindows   

    !Output
    public  :: WriteDataLine                !ASCII Output
    

    !Destructor
    public  :: KillEnterData
    private ::      DeallocateInstance

    !Management
    public  :: Block_Unlock              
    private ::      Block_Lock              

    private :: Ready
    private ::      LocateObjEnterData
    
    !Interface-----------------------------------------------------------------

    private :: ReadIntegerVector_New
    private :: ReadInteger_New

    private :: ReadRealVector_New
    private :: ReadReal_New

    private :: ReadDbleVector_New
    private :: ReadDble_New

    private :: ReadLogical_New

    private :: ReadString_New
    private :: ReadStringVector_New

    private :: ReadTime

    interface  GetData
        module procedure ReadIntegerVector_New
        module procedure ReadInteger_New
        module procedure ReadRealVector_New
        module procedure ReadReal_New
        module procedure ReadDbleVector_New
        module procedure ReadDble_New
        module procedure ReadLogical_New
        module procedure ReadString_New
        module procedure ReadStringVector_New
        module procedure ReadTime
    end interface GetData

    interface WriteDataLine
        module procedure WriteDataLineLogical
        module procedure WriteDataLineInteger
        module procedure WriteDataLineString
        module procedure WriteDataLineBlock
        module procedure WriteDataLineReal
        module procedure WriteDataLineIntVector
        module procedure WriteDataLineRealVector
        module procedure WriteDataLineTime
        module procedure WriteDataLineOutputTime
    end interface

   !Parameter-----------------------------------------------------------------

    integer, parameter :: FORMATTED_     = 1
    integer, parameter :: UNFORMATTED_   = 2

    !Type----------------------------------------------------------------------
    Type       T_Dataline
        character(LEN = line_length) :: full_line
        integer                      :: delimiter_pos = null_int
    end type T_Dataline

    type       T_Block
        integer :: BeginBlock = null_int
        integer :: EndBlock   = null_int
    end type T_Block


    Type      T_EnterData
        integer                         :: InstanceID            
        type (T_Block   )               :: Block
        type (T_Block   )               :: BlockFromBlock
        type (T_Block   )               :: BlockFromBlockFromBlock
      
        type (T_Dataline), pointer      :: BufferLines(:)
        logical                         :: BLOCK_LOCK                   = IDLE
        integer                         :: BlockClientIDnumber          = null_int
        integer                         :: unit                         = null_int
        integer                         :: BufferSize
        character(LEN = PathLength)     :: FileName
        type(T_EnterData), pointer      :: Next
    end type T_EnterData

    !Global Variables
    type (T_EnterData), pointer         :: FirstEnterData
    type (T_EnterData), pointer         :: Me

    !--------------------------------------------------------------------------

    contains



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    ! This soubroutine extracts, the whole data from unit "Me%unit"          
    ! and stores them in sring "ObjEnterData" !

    subroutine ConstructEnterData(EnterDataID, FileName, FORM, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataID
        character(LEN = *), intent(IN )             :: FileName
        integer, optional,  intent(OUT)             :: STAT    
        integer, optional,  intent(IN )             :: FORM

        !Local-----------------------------------------------------------------
        integer                                     :: ready_     
        integer                                     :: STAT_CALL    
        logical                                     :: exists
        character(LEN = line_length)                :: string
        character(LEN = 10000)                      :: auxstring
        character(LEN = 1)                          :: one_char
        integer                                     :: STAT_
        integer                                     :: I, line
        integer                                     :: FORM_

        !----------------------------------------------------------------------

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mEnterData_)) then
            nullify (FirstEnterData)
            call RegisterModule (mEnterData_) 
        endif

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)

cd0 :   if (ready_ .EQ. OFF_ERR_) then


            !Verifies if file exits
            inquire(FILE = trim(adjustl(FileName)), EXIST = exists)

            if (.NOT. exists) then
                write (*,*)'Data File does not exists : ', trim(adjustl(FileName))
                STAT_ = FILE_NOT_FOUND_ERR_
            endif

            call AllocateInstance         

            Me%FileName = FileName

            call UnitsManager(Me%unit, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEnterData - ModuleEnterData - ERR01'


cd16 :      if (present(FORM)) then
                if      (FORM .EQ. FORMATTED_  ) then
                    FORM_ = FORMATTED_
                else if (FORM .EQ. UNFORMATTED_) then    
                    FORM_ = UNFORMATTED_
                else    
                    stop 'ModuleEnterData - ConstructEnterData - ERR06' 
                end if 
            else cd16
               FORM_ = FORMATTED_
            end if cd16

if8 :       if      (FORM_ .EQ. FORMATTED_  ) then
                open(UNIT = Me%unit, FILE = trim(adjustl(Me%FileName)), &
                     FORM = "FORMATTED",   STATUS = "OLD", ACTION = "READ", IOSTAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleEnterData - ConstructEnterData - ERR02' 

                rewind(Me%unit) 
                             
            else if (FORM_ .EQ. UNFORMATTED_) then if8

                open(UNIT = Me%unit, FILE = trim(adjustl(Me%FileName)), &
                     FORM = "UNFORMATTED", STATUS = "OLD", ACTION = "READ", IOSTAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleEnterData - ConstructEnterData - ERR03' 

                rewind(Me%unit)              
                             
            end if if8


            !Counts the number of lines
            Me%BufferSize = 0
            rewind(Me%unit)
do2 :       do 
                read (Me%unit, "(A)", end=100) auxstring
                if(len_trim(trim(auxstring)) > line_length)then
                    write(*,*) 'Maximum of ', len(string),' characters is supported.'
                    write(*,*) 'String: '//trim(auxstring)
                    write(*,*) 'File  : '//trim(adjustl(Me%FileName))
                    write(*,*) 'Line  : ', Me%BufferSize + 1
                    stop 'ModuleEnterData - ConstructEnterData - ERR05' 
                end if
                Me%BufferSize = Me%BufferSize + 1
            end do do2

100         continue


            !Allocates buffer
cd3 :       if (Me%BufferSize .GT. 0) then
                allocate(Me%BufferLines(Me%BufferSize))
            else cd3
                nullify(Me%BufferLines)
            end if cd3

            !Copy file to buffer
            rewind(Me%unit)
do3 :       do I = 1, Me%BufferSize
                read (Me%unit, "(A)") string
                string = adjustl(string)

                do line = 1, line_length

                    one_char = string(line:line)

                    if(one_char == tab)then
                        string(line:line) = space
                    end if

                end do

                Me%BufferLines(I)%full_line = string
            end do do3

            call UnitsManager          (Me%unit, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleEnterData - ConstructEnterData - ERR06' 


            !Returns ID
            EnterDataID     = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0

            stop 'ModuleEnterData - ConstructEnterData - ERR99' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructEnterData

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_EnterData), pointer                 :: NewEnterData
        type (T_EnterData), pointer                 :: PreviousEnterData


        !Allocates new instance
        allocate (NewEnterData)
        nullify  (NewEnterData%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstEnterData)) then
            FirstEnterData          => NewEnterData
            Me                      => NewEnterData
        else
            PreviousEnterData       => FirstEnterData
            Me                      => FirstEnterData%Next
            do while (associated(Me))
                PreviousEnterData   => Me
                Me                  => Me%Next
            enddo
            Me                      => NewEnterData
            PreviousEnterData%Next  => NewEnterData
        endif

        Me%InstanceID = RegisterNewInstance (mENTERDATA_)

    end subroutine AllocateInstance

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine RewindBuffer(EnterDataID, STAT)

        !Parameter-------------------------------------------------------------
        integer                                     :: EnterDataID
        integer, optional, intent(OUT)              :: STAT    

        !Local-----------------------------------------------------------------
        integer                                     :: ready_     
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            Me%Block%BeginBlock                   = null_int
            Me%Block%EndBlock                     = null_int
            Me%BlockFromBlock%BeginBlock          = null_int
            Me%BlockFromBlock%EndBlock            = null_int
            Me%BlockFromBlockFromBlock%BeginBlock = null_int
            Me%BlockFromBlockFromBlock%EndBlock   = null_int
   
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------
    
    end subroutine RewindBuffer

    !--------------------------------------------------------------------------

    subroutine RewindBlock(EnterDataID, ClientNumber, STAT)

        !Parameter-------------------------------------------------------------
        integer                                     :: EnterDataID
        integer, optional                           :: ClientNumber
        integer, optional, intent(OUT)              :: STAT    

        !Local-----------------------------------------------------------------
        integer                                     :: ready_     
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)

cd4 :   if (.NOT. present(ClientNumber)) then
cd1 :       if ( ready_ .EQ. IDLE_ERR_       ) then
            
                Me%BlockFromBlock%BeginBlock          = null_int
                Me%BlockFromBlock%EndBlock            = null_int
                Me%BlockFromBlockFromBlock%BeginBlock = null_int
                Me%BlockFromBlockFromBlock%EndBlock   = null_int
           
                STAT_ = SUCCESS_
            else               
                STAT_ = ready_
            end if cd1
        else 
cd5 :       if (ready_ .NE. OFF_ERR_       ) then  

cd2 :           if (Me%BLOCK_LOCK) then
cd6 :               if (Me%BlockClientIDnumber .EQ. ClientNumber) then
       
                        Me%BlockFromBlock%BeginBlock          = null_int
                        Me%BlockFromBlock%EndBlock            = null_int
                        Me%BlockFromBlockFromBlock%BeginBlock = null_int
                        Me%BlockFromBlockFromBlock%EndBlock   = null_int

                        STAT_ = SUCCESS_
                    else
                        STAT_ = CLIENT_NB_ERR_
                    end if cd6
                else
                    STAT_ = BLOCK_UNLOCK_ERR_
                end if cd2
            else               
                STAT_ = ready_
            end if cd5
        end if cd4


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------
    
    end subroutine RewindBlock


    !--------------------------------------------------------------------------

    subroutine ReplaceFullBufferLine(EnterDataID, BufferLine, FullBufferLine, STAT)
        

        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataID
        integer,            intent (IN)             :: BufferLine
        character(*),       intent (IN)             :: FullBufferLine
        integer,            intent(OUT), optional   :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            Me%BufferLines (BufferLine)%full_line = trim(FullBufferLine)    

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


        !----------------------------------------------------------------------

    end subroutine ReplaceFullBufferLine


    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine CopyBlockSize(EnterDataIDTo, EnterDataIDFrom, ClientNumberFrom, STAT)
        

        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataIDTo
        integer                                     :: EnterDataIDFrom, ClientNumberFrom
        integer,            intent(OUT), optional   :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: StartLine, EndLine, STAT_CALL
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataIDTo, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            call GetBlockSize(EnterDataIDFrom, ClientNumberFrom, StartLine, EndLine,    &
                              SearchType = FromBlock_, STAT = STAT_CALL)

            if(STAT_CALL /= SUCCESS_) stop "GetBlockSize - ModuleEnterData - ERR10"

            Me%Block%BeginBlock = StartLine
            Me%Block%EndBlock   = EndLine

            call GetBlockSize(EnterDataIDFrom, ClientNumberFrom, StartLine, EndLine,    &
                              SearchType = FromBlockInBlock_, STAT = STAT_CALL)

            if(STAT_CALL /= SUCCESS_) stop "GetBlockSize - ModuleEnterData - ERR20"

            Me%BlockFromBlock%BeginBlock = StartLine
            Me%BlockFromBlock%EndBlock   = EndLine


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


        !----------------------------------------------------------------------

    end subroutine CopyBlockSize


    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    !     Esta Rotina destina-se a leitura do NOME do ficheiro FILE_NAME        
    !     definido pela palavra chave KEYWORD no ficheiro nomfich.dat.     
    !                                                                      
    !     Sera emitida a mensagem Message no caso de ser criado de forma   
    !     automatica um ficheiro de output ou de nao existir um ficheiro   
    !     de input.                                                        
    !                                                                      
    !     No caso de criacao automatica do nome do ficheiro, este tem a    
    !     extensao Extensao composta por 3 caracteres.                     
    !                                                                      
    !     TIME_END  -> instante da simulacao para gerar nome.                     

    subroutine ReadFileName(KEYWORD, FILE_NAME, Message, TIME_END, Extension, FilesInput,STAT)

        !Arguments-------------------------------------------------------------

        integer,            optional, intent(OUT) :: STAT     

        character(LEN = *),           intent(OUT) :: FILE_NAME
        character(LEN = *),           intent(IN ) :: KEYWORD
        character(LEN = *), optional              :: Extension
        character(LEN = *), optional              :: Message
        character(LEN = *), optional              :: FilesInput

        type(T_Time),       optional, intent(IN ) :: TIME_END

        !External--------------------------------------------------------------

        integer :: STAT_CALL
        integer :: flag = 0

        character(LEN = line_length) :: Root
        character(LEN = line_length) :: RootMessage
        character(LEN  = PathLength) :: FileNomfich       
        

        integer                      :: ObjEnterData = 0

        !Local-----------------------------------------------------------------

        integer :: STAT_             
        integer :: I

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_
        
        !FilesName inital value is defined in GlobalData 
        !The name is change in main.f90 calling SetFilesName subroutine
        FileNomfich = trim(FilesName)
        
        if (present(FilesInput)) FileNomfich = trim(FilesInput)

        !Files name is defined in the ModuleGlobalData
        call ConstructEnterData(ObjEnterData, trim(FileNomfich), STAT = STAT_CALL)
cd4 :   if      (STAT_CALL .EQ. FILE_NOT_FOUND_ERR_) then
            STAT_ = FILE_NOT_FOUND_ERR_
        
        else if (STAT_CALL .EQ. SUCCESS_           ) then cd4

cd1 :       if (present(Message)) then

                call GetData(FILE_NAME, ObjEnterData, flag, SearchType = FromFile_, text = Message, keyword = KEYWORD)

            else cd1

                call GetData(FILE_NAME, ObjEnterData, flag, SearchType = FromFile_, keyword = KEYWORD)

            end if cd1


cd2 :       if      (flag .EQ. 0) then
cd3 :           if (present(Extension) .AND. present(TIME_END)) then
do1 :               do I = 1, line_length
                        Root(I:I) = space
                    end do do1

                    RootMessage = "Root path not found in nomfich."
                    RootMessage = trim(RootMessage)


                    call GetData(Root, ObjEnterData, flag, SearchType = FromFile_,       &
                                                           text       = RootMessage,     &
                                                           keyword    = "ROOT")
                    if (flag .EQ. 0)                                                     &
                        call GetData(Root, ObjEnterData, flag, SearchType = FromFile_,   &
                                                               text       = RootMessage, &
                                                               keyword    = "RAIZ")
cd5 :               if (flag .EQ. 0) then
                        write(*,*) 
                        write(*,*) "Keyword ROOT not found in nomfich.dat"
                        stop       "SUBROUTINE ReadFileName; Module ModuleEnterData. ERR01"
                    end if cd5


                    call CreateName(FILE_NAME, Root, Extension, TIME_END)

                    FILE_NAME = trim(FILE_NAME)

                    STAT_ = SUCCESS_

                else cd3

                    STAT_ = KEYWORD_NOT_FOUND_ERR_
                end if cd3

            else if (flag .EQ. 1) then cd2

                FILE_NAME = trim(FILE_NAME)

                STAT_ = SUCCESS_
            end if cd2


        end if cd4

        call KillEnterData(ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop "ReadFileName - ModuleEnterData - ERR02"


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ReadFileName

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine CreateName(FILE_NAME, Root, Extension, TIME_END)
        implicit none

        !Arguments--------------------------------------------------------------

        character(LEN = *), intent(OUT) :: FILE_NAME
        character(LEN = *), intent(IN ) :: Root
        character(LEN = *), intent(IN ) :: Extension

        type(T_Time),       intent(IN ) :: TIME_END

        !External--------------------------------------------------------------

        real :: Year, Month, Day, Hour

        !Local-----------------------------------------------------------------

        integer :: length1, length2, length3
        integer :: I

        character(LEN = 12) :: VALOR2
        character(LEN = Line_Length) :: File_Aux

        !----------------------------------------------------------------------

        call ExtractDate(TIME_END, Year = Year, Month = Month, Day = Day, Hour = Hour)

        !Se o ano tem 4 digitos escolhe somente dois (e.g. 1998 -> 98)
         if (Year .GT. 100) &
            Year = Year - int(Year / 100) * 100


        write(VALOR2(1:8),"(4I2)") int(Year), int(Month), int(Day), int(Hour)


do1 :   do I = 1, 8
            IF (VALOR2(I:I) .EQ. space) VALOR2(I:I)="0"
            FILE_NAME(I:I) = VALOR2(I:I)
        end do do1

do2 :   do I = 9, len(FILE_NAME)
            FILE_NAME(I:I) = space
        end do do2

        FILE_NAME(9:9) = dot



        length1 = len_trim(Root     )
        length2 = len_trim(FILE_NAME)
        length3 = len_trim(Extension)
cd3 :   if (len(FILE_NAME) .LT. (length1 + length2 + length3)) then
            write(*,*)  
            write(*,*) "Root too long. "
            stop       "SUBROUTINE CreateName; Module ModuleEnterData. ERR01."
        end if cd3

     
        FILE_AUX  = trim(Root)//trim(FILE_NAME)//trim(Extension)
        FILE_NAME = trim(FILE_AUX)

        !----------------------------------------------------------------------

    end subroutine CreateName

    !--------------------------------------------------------------------------

    subroutine GetFileName(EnterDataID, FileName, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataID
        character(LEN = *), intent(OUT)             :: FileName
        integer,            intent(OUT), optional   :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------
        
        STAT_       = UNKNOWN_

        call Ready(EnterDataID, ready_)
        
cd0 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            FileName = Me%FileName

            STAT_ = SUCCESS_

        else cd0

            STAT_ = ready_
        end if cd0


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetFileName

    !--------------------------------------------------------------------------

    subroutine GetExtractType(FromFile, FromBlock, FromBlockInBlock)

        !Arguments-------------------------------------------------------------

        integer, optional, intent(OUT) :: FromFile      
        integer, optional, intent(OUT) :: FromBlock       
        integer, optional, intent(OUT) :: FromBlockInBlock       

        !----------------------------------------------------------------------

        if (present(FromFile        )) FromFile         = FromFile_
        if (present(FromBlock       )) FromBlock        = FromBlock_
        if (present(FromBlockInBlock)) FromBlockInBlock = FromBlockInBlock_

        !----------------------------------------------------------------------

    end subroutine GetExtractType

    !--------------------------------------------------------------------------

    subroutine GetBufferSize(EnterDataID, BufferSize, STAT)
        
        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataID
        integer,            intent (OUT)            :: BufferSize
        integer,            intent(OUT), optional   :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            BufferSize = Me%BufferSize

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetBufferSize

    !--------------------------------------------------------------------------

    subroutine GetFullBufferLine(EnterDataID, BufferLine, FullBufferLine, STAT)
        

        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataID
        integer,            intent (IN)             :: BufferLine
        character(*),       intent (OUT)            :: FullBufferLine
        integer,            intent(OUT), optional   :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            FullBufferLine = trim(Me%BufferLines (BufferLine)%full_line)

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetFullBufferLine


    !--------------------------------------------------------------------------

    subroutine GetBlockSize(EnterDataID, ClientNumber, StartLine, EndLine, SearchType, STAT)
        
        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataID
        integer,           intent (IN)              :: ClientNumber
        integer,           intent (OUT)             :: StartLine
        integer,           intent (OUT)             :: Endline
        integer, optional, intent (IN)              :: SearchType
        integer, optional, intent (OUT)             :: STAT


        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)
        
cd1 :   if (ready_ .NE. OFF_ERR_) then

cd2 :       if (Me%BLOCK_LOCK) then
cd3 :           if (ClientNumber .NE. Me%BlockClientIDnumber) then

                    write(*,*) 
                    write(*,*) "Client ID mismatch."
                    stop       "Subroutine GetBlockSize; module ModuleEnterData. ERR01."

                else

                    if (present(SearchType)) then
                        
                        if (SearchType == FromBlock_       ) then

                            StartLine = Me%Block%BeginBlock
                            EndLine   = Me%Block%EndBlock

                        else if (SearchType .EQ. FromFile_        ) then

                            StartLine = 1
                            EndLine   = Me%BufferSize

                        else if (SearchType .EQ. FromBlockInBlock_) then

                            StartLine = Me%BlockFromBlock%BeginBlock
                            EndLine   = Me%BlockFromBlock%EndBlock

                        endif

                    else
                    
                        StartLine = Me%Block%BeginBlock
                        EndLine   = Me%Block%EndBlock

                    endif

                    STAT_ = SUCCESS_
         
                end if cd3

            else

                STAT_ = BLOCK_UNLOCK_ERR_

            end if cd2
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetBlockSize

    !--------------------------------------------------------------------------     
    !This soubroutine receives the buffer line and check if there is a keywrod 
    !and if there is a keyword return a string with the keyword.

    subroutine GetKeywordFromLine(EnterDataID, BufferLine, ExistKeyword, Keyword, STAT)
        

        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataID
        integer,            intent (IN)             :: BufferLine
        logical,            intent (OUT)            :: ExistKeyword
        character(*),       intent (OUT)            :: Keyword
        integer,            intent(OUT), optional   :: STAT

        !Local-----------------------------------------------------------------
        type(T_Dataline)                            :: string
        integer                                     :: ready_          
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            read(Me%BufferLines(BufferLine)%full_line,"(A)") string%full_line

            !locate the symbol delimiting the keyword
            call ScanLine(string, ExistKeyword)

            if (ExistKeyword) then

                Keyword = string%full_line(1:string%delimiter_pos-1)

            endif

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetKeywordFromLine                                              


    !--------------------------------------------------------------------------     
    !This soubroutine receives a string with Me%BufferSize (Me%BufferLines) and extracts 
    ! the line started by "Keyword". If a line is found "flag is "1".       
    ! Otherwise "flag is "0".                                              

    subroutine LookForLineInBlock(keyword, string, SearchType_, flag, CaseSensitive)

        !Parameter-------------------------------------------------------------

        character(LEN = *), intent(IN ) :: keyword

        integer,            intent(IN ) :: SearchType_
        integer,            intent(OUT) :: flag
        logical, optional,  intent(IN ) :: CaseSensitive           

        !External--------------------------------------------------------------

        type(T_Dataline), intent(OUT)   :: string

        !Local-----------------------------------------------------------------

        integer                         :: I
        integer                         :: StartSearch, EndSearch
        character(LEN=line_length)      :: ReadKeyWord
        logical                         :: CaseSensitiveAux, ValidDelimiter

        !----------------------------------------------------------------------

        flag   = 0

cd1 :   if      (SearchType_ .EQ. FromBlock_       ) then

            StartSearch = Me%Block%BeginBlock
            EndSearch   = Me%Block%EndBlock

        else if (SearchType_ .EQ. FromFile_        ) then

            StartSearch = 1
            EndSearch   = Me%BufferSize

        else if (SearchType_ .EQ. FromBlockInBlock_) then

            StartSearch = Me%BlockFromBlock%BeginBlock
            EndSearch   = Me%BlockFromBlock%EndBlock

        else
            write(*,*)  
            write(*,*) "Error SearchType_ = 0. "
            stop       "Subroutine LookForLineInBlock; Module ModuleEnterData. ERR01."
        end if cd1


do1 :   do I = StartSearch, EndSearch
            read(Me%BufferLines(I)%full_line,"(A)") string%full_line

            !locate the symbol delimiting the keyword
            call ScanLine(string, ValidDelimiter)

            ReadKeyWord = string%full_line(1:string%delimiter_pos-1)

            if (present(CaseSensitive)) then

                CaseSensitiveAux = CaseSensitive

            else
                !By default the enter data module is case sensitive
                CaseSensitiveAux = .true.


            endif
            

            if (KeywordsEqual(trim(ReadKeyWord), trim(keyword), CaseSensitiveAux)) then
                if(ValidDelimiter)then
                    flag = 1
                    exit do1
                else
                    write(*,*)'Invalid delimiter between keyword and value.'
                    write(*,*)'Line         = ', trim(string%full_line)
                    write(*,*)'Line Number  = ', i
                    write(*,*)'File         = ', trim(Me%FileName)
                    stop      'Subroutine LookForLineInBlock; Module ModuleEnterData. ERR02.'
                end if
            end if 
        end do do1

        !----------------------------------------------------------------------

    end subroutine LookForLineInBlock

    !--------------------------------------------------------------------------     

    function KeywordsEqual(ReadKey, ReferenceKey, CaseSentive)
    logical :: KeywordsEqual

    !Arguments --------------------------------------------------------------------------

    character (Len = *) :: ReadKey, ReferenceKey
    logical             :: CaseSentive

    !Local ------------------------------------------------------------------------------
    character (Len = 1) :: a, b
    integer             :: i, ia, ib, DimRef

   
    if (.not. CaseSentive) then

        DimRef = LEN(ReferenceKey)

        if (DimRef == LEN(ReadKey)) then

            do i = 1, DimRef

                a    = ReadKey     (i:i)
                b    = ReferenceKey(i:i)

                ia   = iachar(a)
                ib   = iachar(b)

                if ( ia >= 65 .and. ia <= 90) ReadKey     (i:i) = Achar(ia + 32)
                if ( ib >= 65 .and. ib <= 90) ReferenceKey(i:i) = Achar(ib + 32)

            enddo

        endif

    endif

    if (ReadKey == ReferenceKey) then 

        KeywordsEqual =.true.

    else 

        KeywordsEqual =.false.


    endif

    end function KeywordsEqual



    !--------------------------------------------------------------------------     
    !subroutine to locate the symbol delimiting the keyword from data in 
    ! the string extracted from an input file                            

    subroutine ScanLine(string, ValidDelimiter)

        !External--------------------------------------------------------------
        type(T_Dataline)    :: string
        logical, optional   :: ValidDelimiter

        !Local-----------------------------------------------------------------

        integer :: delimiter_pos
        integer :: IUB

        character(LEN = 1) :: space     = char(32)
        character(LEN = 1) :: slash     = char(47)
        character(LEN = 1) :: backslash = char(92)

        !----------------------------------------------------------------------
        if(present(ValidDelimiter))ValidDelimiter = .true.

        delimiter_pos = scan(string%full_line, delimiter)

        IUB = len_trim(string%full_line)

cd1 :   if       (delimiter_pos .EQ. 0) then
            delimiter_pos = scan(string%full_line, space)
            if(present(ValidDelimiter))ValidDelimiter = .false.

        else if ((delimiter_pos .GT. 0) .AND. (delimiter_pos .LT. IUB)) then
            if ((string%full_line(delimiter_pos+1:delimiter_pos+1) .EQ. slash    ) .OR. &   
                (string%full_line(delimiter_pos+1:delimiter_pos+1) .EQ. backslash))     &   
                delimiter_pos = scan(string%full_line, space)

        else if                               (delimiter_pos .EQ. IUB)  then
            write(*,*)  
            write(*,*) "Found : at the end of line." 
            write(*,*) "Line: "//trim(string%full_line) 
            stop       "Subroutine ScanLine; Module ModuleEnterData. ERR01." 
        end if cd1

        string%delimiter_pos = delimiter_pos

        !----------------------------------------------------------------------

    end subroutine ScanLine

    !--------------------------------------------------------------------------     

    subroutine ReadInteger_New(value,                                         &
                               EnterDataID,                                   &
                               flag,                                          &
                               SearchType,                                    &
                               text,                                          &
                               Buffer_Line,                                   &
                               keyword,                                       &
                               Default,                                       &
                               ClientModule,                                  & !Future use for interface
                               CaseSensitive,                                 &
                               STAT)
      
        !Parameter-------------------------------------------------------------
        integer                                     :: EnterDataID 
        character(LEN = *), optional, intent(IN )   :: keyword
        character(LEN = *), optional, intent(IN )   :: text
        character(LEN = *), optional, intent(IN )   :: ClientModule
        logical           , optional, intent(IN )   :: CaseSensitive
        integer, optional,       intent(IN )        :: Buffer_Line
        integer, optional,       intent(IN )        :: Default
        integer,                 intent(OUT)        :: flag 
        integer, optional,       intent(OUT)        :: STAT
        integer, optional,       intent(IN )        :: SearchType
        integer,                 intent(OUT)        :: value

        !External--------------------------------------------------------------
        integer                                     :: SearchTypeOUT
        character(StringLength)                     :: ClientModuleOUT, DefaultOUT, ValueOUT
        logical                                     :: CaseSensitiveOUT
        integer                                     :: STAT_CALL
        integer                                     :: ready_

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable
        integer :: length          
                                        
        character (len=line_length) :: aux_str = ""

        Type (T_Dataline) :: string
          
        !----------------------------------------------------------------------

        STAT_       = UNKNOWN_

        length      = null_int

        call Ready(EnterDataID, ready_)
        
cd0 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            value = null_int
            flag   = 0


if9 :       if (present(SearchType)) then
                SearchTypeOut = SearchType
            else
                SearchTypeOut = FromFile_
            end if if9


            if (.not.present(keyword).and..not.present(Buffer_Line))          &
              stop "GetData (ReadInteger_New); ModuleEnter_Data. ERR01."
            
            if (present(keyword).and.present(Buffer_Line))                    &
              stop "GetData (ReadInteger_New); ModuleEnter_Data. ERR02."

            if (present(keyword)) then

                if (present (CaseSensitive)) then

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag, CaseSensitive)

                else

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag)

                endif

            elseif (present(Buffer_Line)) then

                 string%full_line     = Me%BufferLines(Buffer_Line)%full_line
                 string%delimiter_pos = 0
                 flag                 = 1

            endif

cd2 :       if (flag == 1) then
                read (string%full_line(string%delimiter_pos + 1 : line_length),"(A)") aux_str
                aux_str = adjustl (aux_str)
                length  = len_trim(aux_str)
                

                read (aux_str(1 : length),*, IOSTAT = STAT_CALL, err=100) value
100             continue
cd3 :           if (STAT_CALL .EQ. SUCCESS_) then
                    STAT_ = SUCCESS_
                else
                    STAT_ = SIZE_ERR_
                end if cd3

            else
cd6 :           if (present(Default)) then
                    value = Default

                    STAT_ = SUCCESS_
                else
cd4 :               if (present(text)) then
                        call WriteErrorMessage(keyword, text, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                              &
                          stop "GetData (ReadInteger_New); ModuleEnter_Data. ERR03."
                    end if cd4

                    STAT_ = SUCCESS_
                end if cd6
            end if cd2

            if (present(keyword)) then
                
                if (present(ClientModule)) then
                    ClientModuleOUT = ClientModule
                else
                    ClientModuleOUT = "Unknown Client"
                endif

                if (present(Default)) then
                    write(DefaultOUT, fmt=*) Default
                else
                    DefaultOUT = "Do not have"
                endif

                if (present(CaseSensitive)) then
                    CaseSensitiveOUT = CaseSensitive
                else
                    CaseSensitiveOUT = .true.
                endif

                write (ValueOUT, fmt=*) value
                
                call LogKeyWord (keyword, SearchTypeOUT, ClientModuleOUT, DefaultOUT, CaseSensitiveOUT, ValueOUT)

            endif

        else 
            STAT_ = ready_
        end if cd0


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ReadInteger_New

    !--------------------------------------------------------------------------     

    subroutine ReadString_New(value,                                          &
                              EnterDataID,                                    &
                              flag,                                           &
                              SearchType,                                     &
                              text,                                           &
                              Buffer_Line,                                    &
                              keyword,                                        &
                              Default,                                        &
                              ClientModule,                                   & !Future use for interface
                              CaseSensitive,                                  &
                              STAT)
      
        !Parameter-------------------------------------------------------------
        integer                                     :: EnterDataID
        character(LEN = *), optional, intent(IN )   :: keyword
        character(LEN = *), optional, intent(IN )   :: text
        character(LEN = *), optional, intent(IN )   :: ClientModule
        character(LEN = *), optional, intent(IN )   :: Default
        logical           , optional, intent(IN )   :: CaseSensitive
        integer, optional,       intent(IN )        :: Buffer_Line
        integer,                 intent(OUT)        :: flag
        integer, optional,       intent(OUT)        :: STAT
        integer, optional,       intent(IN )        :: SearchType
        character(LEN = *),      intent(OUT)        :: value

        !External--------------------------------------------------------------
        integer                                     :: SearchTypeOUT
        character(StringLength)                     :: ClientModuleOUT, DefaultOUT, ValueOUT
        logical                                     :: CaseSensitiveOUT
        integer                                     :: STAT_CALL
        integer                                     :: ready_

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable
        integer :: length          
        integer :: length_Aux, length_Out, I
                                        
        character (len=line_length)                 :: aux_str = ""     
        Type (T_Dataline)                           :: string
          
        !----------------------------------------------------------------------

        STAT_       = UNKNOWN_
        length      = null_int

        call Ready(EnterDataID, ready_)
        
cd0 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            value = trim(adjustl("************"))
            flag  = 0


if9 :       if (present(SearchType)) then
                SearchTypeOut = SearchType
            else
                SearchTypeOut = FromFile_
            end if if9


            if (.not.present(keyword).and..not.present(Buffer_Line))          &
              stop "GetData (ReadString_New); ModuleEnter_Data. ERR01."
            
            if (present(keyword).and.present(Buffer_Line))                    &
              stop "GetData (ReadString_New); ModuleEnter_Data. ERR02."

            if (present(keyword)) then

                if (present (CaseSensitive)) then

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag, CaseSensitive)

                else

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag)

                endif

            elseif (present(Buffer_Line)) then

                 string%full_line     = Me%BufferLines(Buffer_Line)%full_line
                 string%delimiter_pos = 0
                 flag                 = 1

            endif

cd1 :       if (flag == 1) then

                read (string%full_line(string%delimiter_pos + 1 : line_length),"(A)",err=100, iostat=STAT_CALL) aux_str
100             continue
cd3 :           if (STAT_CALL .EQ. SUCCESS_) then
                    aux_str = adjustl(aux_str)
                
                    length_Aux = len_trim(aux_str)

                    length_Out=len (value)

                    if (length_Aux>length_Out) then
                        write(*,*) trim(aux_str)
                        stop "GetData (ReadString_New); ModuleEnter_Data. ERR04."
                    endif

                    do i=1,length_Out
                       value(i:i) = Space
                    enddo
 
                    do i = 1, length_Aux
                       read (aux_str(i:i),"(A1)") value(i:i)
                    end do 

                    STAT_ = SUCCESS_
                else
                    STAT_ = SIZE_ERR_
                end if cd3

            else
cd6 :           if (present(Default)) then
                    value = Default

                    STAT_ = SUCCESS_
                else
if8 :               if (present(text)) then
                        call WriteErrorMessage(keyword, text, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                              &
                            stop "GetData (ReadString_New); ModuleEnter_Data. ERR03."
                    end if if8

                    STAT_ = SUCCESS_
                end if cd6
            end if cd1

            if (present(keyword)) then
                
                if (present(ClientModule)) then
                    ClientModuleOUT = ClientModule
                else
                    ClientModuleOUT = "Unknown Client"
                endif

                if (present(Default)) then
                    DefaultOUT = Default
                else
                    DefaultOUT = "Do not have"
                endif

                if (present(CaseSensitive)) then
                    CaseSensitiveOUT = CaseSensitive
                else
                    CaseSensitiveOUT = .true.
                endif

                ValueOUT = value
                
                call LogKeyWord (keyword, SearchTypeOUT, ClientModuleOUT, DefaultOUT, CaseSensitiveOUT, ValueOUT)

            endif

        else 
            STAT_ = ready_
        end if cd0


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ReadString_New


    !--------------------------------------------------------------------------     

    subroutine ReadStringVector_New(vector,                                   &
                                  EnterDataID,                                &
                                  flag,                                       &
                                  SearchType,                                 &
                                  text,                                       &
                                  Buffer_Line,                                &
                                  keyword,                                    &
                                  Default,                                    &
                                  ClientModule,                               & !Future use for interface
                                  CaseSensitive,                              &
                                  STAT)
      
        !Parameter-------------------------------------------------------------
        integer                                     :: EnterDataID    
        character(LEN = *), optional, intent(IN )   :: keyword
        character(LEN = *), optional, intent(IN )   :: text
        character(LEN = *), optional, intent(IN )   :: ClientModule
        logical           , optional, intent(IN )   :: CaseSensitive

        integer, optional,       intent(IN ) :: Buffer_Line

        integer,                 intent(OUT) :: flag 
        integer, optional,       intent(OUT) :: STAT
        integer, optional,       intent(IN ) :: SearchType

        character(LEN = *), dimension(:), intent(OUT) :: vector
        character(LEN = *), optional,     intent(IN ) :: Default


        !External--------------------------------------------------------------
        character(StringLength)                     :: ClientModuleOUT, DefaultOUT        
        logical                                     :: CaseSensitiveOUT
        integer                                     :: SearchTypeOUT
        integer                                     :: STAT_CALL
        integer                                     :: ready_          

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable
        integer :: vector_size      
        integer :: length, i, n_data            
                                        
        character (len=line_length) :: aux_str      = ""

        Type (T_Dataline) :: string
          
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        length      = null_int
        vector_size = null_int
        n_data      = null_int


        call Ready(EnterDataID, ready_)
        
cd0 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            vector = null_str
            flag   = 0


if9 :       if (present(SearchType)) then
                SearchTypeOut = SearchType
            else
                SearchTypeOut = FromFile_
            end if if9


            if (.not.present(keyword).and..not.present(Buffer_Line))          &
              stop "Error 01 ReadDbleVector_New - ModuleEnter_Data"
            
            if (present(keyword).and.present(Buffer_Line))                    &
              stop "Error 02 ReadDbleVector_New - ModuleEnter_Data"

            if (present(keyword)) then

                if (present (CaseSensitive)) then

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag, CaseSensitive)

                else

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag)

                endif

            elseif (present(Buffer_Line)) then

                 string%full_line     = Me%BufferLines(Buffer_Line)%full_line
                 string%delimiter_pos = 0
                 flag                 = 1

            endif

cd2 :       if (flag == 1) then
                vector_size = size(vector)

                read (string%full_line(string%delimiter_pos + 1 : line_length),"(A)") aux_str
                aux_str = adjustl (aux_str)
                length  = len_trim(aux_str)
                
                n_data  = 0

cd1 :           if (length > 0) then ! counts the number of data
                    n_data = 1
do1 :                   do i = 2, length
                            if((aux_str (i:i) == " ") .AND.                   &
                               (aux_str (i-1:i-1) /= " "))                    &
                                n_data = n_data + 1 
                        end do do1
                end if cd1

                read (aux_str(1 : length),*, err=100) (vector (i), i = 1, min (n_data, vector_size))
                flag = min (n_data, vector_size) 
            else
cd6 :           if (present(Default)) then
                    vector = Default

                    STAT_ = SUCCESS_
                else
cd4 :               if (present(text)) then
                        call WriteErrorMessage(keyword, text, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                              &
                            stop "Subroutine ReadDbleVector_New; module ModuleEnterData. ERR01."
                    end if cd4

                    STAT_ = KEYWORD_NOT_FOUND_ERR_
                end if cd6
            end if cd2

cd3 :       if (vector_size .NE. n_data) then
                STAT_ = SIZE_ERR_
            else 
                STAT_ = SUCCESS_
            end if cd3


            if (present(keyword)) then
                
                if (present(ClientModule)) then
                    ClientModuleOUT = ClientModule
                else
                    ClientModuleOUT = "Unknown Client"
                endif

                if (present(Default)) then
                    write(DefaultOUT, fmt=*) Default
                else
                    DefaultOUT = "Do not have"
                endif

                if (present(CaseSensitive)) then
                    CaseSensitiveOUT = CaseSensitive
                else
                    CaseSensitiveOUT = .true.
                endif

!                if (vector_size > 0) then
!                    write (ValueOUT, fmt=*, err=10) vector (1:min (n_data, vector_size))
!                     write (ValueOUT, *) vector (1:min (n_data, vector_size))                                           
!                else
!                    ValueOUT = ' '
!                endif

!10              continue
                
!                call LogKeyWord (keyword, SearchTypeOUT, ClientModuleOUT, DefaultOUT, CaseSensitiveOUT, ValueOUT)

            endif

        else 
            STAT_ = ready_
        end if cd0


100     continue


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ReadStringVector_New

    !--------------------------------------------------------------------------     

    subroutine ReadLogical_New(value,                                         &
                               EnterDataID,                                   &
                               flag,                                          &
                               SearchType,                                    &
                               text,                                          &
                               Buffer_Line,                                   &
                               keyword,                                       &
                               Default,                                       &
                               ClientModule,                                  & !Future use for interface
                               CaseSensitive,                                 &
                               STAT)
      
        !Parameter-------------------------------------------------------------
        integer                                     :: EnterDataID
        character(LEN = *), optional, intent(IN )   :: keyword
        character(LEN = *), optional, intent(IN )   :: text
        character(LEN = *), optional, intent(IN )   :: ClientModule
        logical           , optional, intent(IN )   :: CaseSensitive
        integer, optional,       intent(IN )        :: Buffer_Line
        integer,                 intent(OUT)        :: flag 
        integer, optional,       intent(OUT)        :: STAT
        integer, optional,       intent(IN )        :: SearchType
        logical,                 intent(OUT)        :: value
        logical, optional,       intent(IN )        :: Default

        !External--------------------------------------------------------------
        character(StringLength)                     :: ClientModuleOUT, DefaultOUT, ValueOUT
        logical                                     :: CaseSensitiveOUT
        integer                                     :: SearchTypeOUT
        integer                                     :: STAT_CALL
        integer                                     :: ready_          

        !Local-----------------------------------------------------------------

        integer :: STAT_              
        integer :: length          
        integer :: aux_int
                                        
        Type (T_Dataline) :: string
          
        !----------------------------------------------------------------------

        STAT_       = UNKNOWN_

        length      = null_int

        flag   = 0


        call Ready(EnterDataID, ready_)
        
cd0 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

if9 :       if (present(SearchType)) then
                SearchTypeOut = SearchType
            else
                SearchTypeOut = FromFile_
            end if if9


            if (.not.present(keyword).and..not.present(Buffer_Line))          &
              stop "GetData (ReadLogical_New); ModuleEnter_Data. ERR01."
            
            if (present(keyword).and.present(Buffer_Line))                    &
              stop "GetData (ReadLogical_New); ModuleEnter_Data. ERR02."

            if (present(keyword)) then

                if (present (CaseSensitive)) then

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag, CaseSensitive)

                else

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag)

                endif

            elseif (present(Buffer_Line)) then

                 string%full_line     = Me%BufferLines(Buffer_Line)%full_line
                 string%delimiter_pos = 0
                 flag                 = 1

            endif

cd1 :       if (flag == 1) then
                read (string%full_line(string%delimiter_pos+1 : line_length), *, IOSTAT = STAT_CALL, &
                                                                                 ERR    = 100) aux_int
100             continue
cd4 :           if (STAT_CALL .NE. SUCCESS_) then
cd5 :           if (present(Buffer_Line)   ) then   !If reading from a specified line it is necessary to have 
                                                    !  a logical value explictly written in data file
                    aux_int = 0
                else
                    aux_int = 1
                end if cd5
                end if cd4


cd2 :           if      (aux_int .EQ. 1) then

                    value = .TRUE.

                else if (aux_int .EQ. 0) then cd2 

                    value = .FALSE.
                else

                    write(*,*    ) 
                    write(*,*    ) "Value in front of a logical keyword must be 0 or 1."
                    write(*,"(A)") "Keyword: ", keyword
                    stop           "Subroutine GetData (ReadLogical_New); module ModuleEnterData. ERR03."
                end if cd2

                STAT_ = SUCCESS_
            else
cd6 :           if (present(Default)) then

                    value = Default

                    STAT_ = SUCCESS_

                else
                    value = .FALSE.

cd3 :               if (present(text)) then
                        call WriteErrorMessage(keyword, text, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                              &
                            stop "Subroutine GetData (ReadLogical_New); module ModuleEnterData. ERR04."
                    end if cd3

                    STAT_ = SUCCESS_
                end if cd6
            end if cd1

            if (present(keyword)) then
                
                if (present(ClientModule)) then
                    ClientModuleOUT = ClientModule
                else
                    ClientModuleOUT = "Unknown Client"
                endif

                if (present(Default)) then
                    if (Default) then
                        DefaultOUT = "TRUE"
                    else
                        DefaultOUT = "FALSE"
                    endif
                else
                    DefaultOUT = "Do not have"
                endif

                if (present(CaseSensitive)) then
                    CaseSensitiveOUT = CaseSensitive
                else
                    CaseSensitiveOUT = .true.
                endif

                if (value) then
                    ValueOUT = "TRUE"
                else
                    ValueOUT = "FALSE"
                endif
                
                call LogKeyWord (keyword, SearchTypeOUT, ClientModuleOUT, DefaultOUT, CaseSensitiveOUT, ValueOUT)

            endif

        else 
            STAT_ = ready_
        end if cd0


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ReadLogical_New

    !--------------------------------------------------------------------------     

    !--------------------------------------------------------------------------     

    subroutine ReadDble_New(value,                                            &
                            EnterDataID,                                      &
                            flag,                                             &
                            SearchType,                                       &
                            text,                                             &
                            Buffer_Line,                                      &
                            keyword,                                          &
                            Default,                                          &
                            ClientModule,                                     & !Future use for interface
                            CaseSensitive,                                    &
                            STAT)
      
        !Parameter-------------------------------------------------------------
        integer                                     :: EnterDataID
        character(LEN = *), optional, intent(IN )   :: keyword
        character(LEN = *), optional, intent(IN )   :: text
        character(LEN = *), optional, intent(IN )   :: ClientModule
        logical           , optional, intent(IN )   :: CaseSensitive
        integer, optional,       intent(IN )        :: Buffer_Line
        integer,                 intent(OUT)        :: flag 
        integer, optional,       intent(OUT)        :: STAT
        integer, optional,       intent(IN )        :: SearchType
        real(8),                 intent(OUT)        :: value
        real(8), optional,       intent(IN )        :: Default


        !External--------------------------------------------------------------
        character(StringLength)                     :: ClientModuleOUT, DefaultOUT, ValueOUT
        logical                                     :: CaseSensitiveOUT
        integer                                     :: SearchTypeOUT
        integer                                     :: STAT_CALL
        integer                                     :: ready_          

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable
        integer :: length          
                                        
        character (len=line_length) :: aux_str = ""

        Type (T_Dataline) :: string
          
        !----------------------------------------------------------------------

        STAT_       = UNKNOWN_

        length      = null_int



        call Ready(EnterDataID, ready_)
        
cd0 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            
            value = null_real
            flag   = 0


if9 :       if (present(SearchType)) then
                SearchTypeOut = SearchType
            else
                SearchTypeOut = FromFile_
            end if if9


            if (.not.present(keyword).and..not.present(Buffer_Line))          &
              stop "GetData (ReadDble_New); ModuleEnter_Data. ERR01."
            
            if (present(keyword).and.present(Buffer_Line))                    &
              stop "GetData (ReadDble_New); ModuleEnter_Data. ERR02."

            if (present(keyword)) then

                if (present (CaseSensitive)) then

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag, CaseSensitive)

                else

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag)

                endif

            elseif (present(Buffer_Line)) then

                 string%full_line     = Me%BufferLines(Buffer_Line)%full_line
                 string%delimiter_pos = 0
                 flag                 = 1

            endif

cd2 :       if (flag == 1) then
                read (string%full_line(string%delimiter_pos + 1 : line_length),"(A)") aux_str
                aux_str = adjustl (aux_str)
                length  = len_trim(aux_str)
                

                read (aux_str(1 : length),*, IOSTAT = STAT_CALL, err=100) value
100             continue
cd3 :           if (STAT_CALL .EQ. SUCCESS_) then
                    STAT_ = SUCCESS_
                else
                    STAT_ = SIZE_ERR_
                end if cd3

            else
cd6 :           if (present(Default)) then
                    value = Default

                    STAT_ = SUCCESS_
                else
cd4 :               if (present(text)) then
                        call WriteErrorMessage(keyword, text, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                              &
                          stop "GetData (ReadDble_New); ModuleEnter_Data. ERR03."
                    end if cd4

                    STAT_ = SUCCESS_
                end if cd6
            end if cd2

            if (present(keyword)) then
                
                if (present(ClientModule)) then
                    ClientModuleOUT = ClientModule
                else
                    ClientModuleOUT = "Unknown Client"
                endif

                if (present(Default)) then
                    write(DefaultOUT, fmt=*) Default
                else
                    DefaultOUT = "Do not have"
                endif

                if (present(CaseSensitive)) then
                    CaseSensitiveOUT = CaseSensitive
                else
                    CaseSensitiveOUT = .true.
                endif

                write (ValueOUT, fmt=*) value
                
                call LogKeyWord (keyword, SearchTypeOUT, ClientModuleOUT, DefaultOUT, CaseSensitiveOUT, ValueOUT)

            endif

        else 
            STAT_ = ready_
        end if cd0


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ReadDble_New

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine ReadIntegerVector_New(Vector,                                  &
                                     EnterDataID,                             &
                                     flag,                                    &
                                     SearchType,                              &
                                     text,                                    &
                                     Buffer_Line,                             &
                                     keyword,                                 &
                                     Default,                                 &
                                     ClientModule,                            & !Future use for interface
                                     CaseSensitive,                           &
                                     STAT)

        !Arguments-------------------------------------------------------------
        integer, dimension(:), intent(OUT)            :: vector
        integer                                       :: EnterDataID
        integer           , intent(OUT)               :: flag         
        integer           , intent(IN ), optional     :: SearchType
        character(LEN = *), intent(IN ), optional     :: text
        integer           , intent(IN ), optional     :: Buffer_Line        
        character(LEN = *), intent(IN ), optional     :: keyword
        integer           , intent(IN ), optional     :: Default        
        character(LEN = *), intent(IN ), optional     :: ClientModule
        logical           , intent(IN ), optional     :: CaseSensitive
        integer           , intent(OUT), optional     :: STAT
        

        !Local-----------------------------------------------------------------
        character(StringLength)                       :: ClientModuleOUT, DefaultOUT, ValueOUT
        logical                                       :: CaseSensitiveOUT
        integer                                       :: SearchTypeOUT
        integer                                       :: STAT_CALL
        integer                                       :: ready_            

        Type(T_Dataline) :: string


        integer :: vector_size        
        integer :: length, i, n_data       
        integer :: STAT_              !Auxiliar local variable
                                        
        character (len=line_length) :: aux_str = ""
                    
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_


        call Ready(EnterDataID, ready_)
        
cd0 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            flag   = 0
            Vector = null_int
            n_data = null_int


if9 :       if (present(SearchType)) then
                SearchTypeOut = SearchType
            else
                SearchTypeOut = FromFile_
            end if if9

            if (.not.present(keyword).and..not.present(Buffer_Line))          &
              stop "Error 01 ReadIntegerVector_New - ModuleEnter_Data"
            
            if (present(keyword).and.present(Buffer_Line))                    &
              stop "Error 02 ReadIntegerVector_New - ModuleEnter_Data"

            if (present(keyword)) then

                if (present (CaseSensitive)) then

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag, CaseSensitive)

                else

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag)

                endif

            elseif (present(Buffer_Line)) then

                 string%full_line     = Me%BufferLines(Buffer_Line)%full_line
                 string%delimiter_pos = 0
                 flag                 = 1

            endif
    
            vector_size = size   (vector)
cd1 :       if (flag == 1) then
                ! vector_size = size   (vector) ! GR : This line is wrong here as the vector_size 
                                                ! is required whether flag is 1 or 0
 
                read (string%full_line(string%delimiter_pos + 1 : line_length),"(A)") aux_str
                aux_str = adjustl (aux_str)
                length  = len_trim(aux_str)

cd2 :           if (length > 0) then ! counts the number of data
                    n_data = 1

do1 :               do i = 2, length
                        if((aux_str (i:i) == " ") .and. (aux_str (i-1:i-1) /= " ")) &
                            n_data = n_data + 1 
                    end do do1
                end if cd2

                read (aux_str(1 : length),*,err=100) (vector (i), i = 1, min (n_data, vector_size))
                flag = min (n_data, vector_size) 



cd3 :           if (vector_size .NE. n_data) then
                    STAT_ = SIZE_ERR_
                else 
                    STAT_ = SUCCESS_
                end if cd3

            else
cd6 :           if (present(Default)) then
                    vector = Default

                    STAT_ = SUCCESS_
                else
cd4 :               if (present(text)) then
                        call WriteErrorMessage(keyword, text, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                              &
                            stop "Subroutine ReadIntegerVector_New; module ModuleEnterData. ERR01."
                    end if cd4

                    STAT_ = KEYWORD_NOT_FOUND_ERR_
                end if cd6
            end if cd1 

            if (present(keyword)) then
                
                if (present(ClientModule)) then
                    ClientModuleOUT = ClientModule
                else
                    ClientModuleOUT = "Unknown Client"
                endif

                if (present(Default)) then
                    write(DefaultOUT, fmt=*) Default
                else
                    DefaultOUT = "Do not have"
                endif

                if (present(CaseSensitive)) then
                    CaseSensitiveOUT = CaseSensitive
                else
                    CaseSensitiveOUT = .true.
                endif
                
                ! GR : vector_size must be properly initialized before this line
                if (vector_size > 0) then
                    write (ValueOUT, fmt=*) vector (1:min (n_data, vector_size))
                else
                    ValueOUT = ' '
                endif
                
                call LogKeyWord (keyword, SearchTypeOUT, ClientModuleOUT, DefaultOUT, CaseSensitiveOUT, ValueOUT)

            endif

        else 
            STAT_ = ready_
        end if cd0


100     continue


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ReadIntegerVector_New

    !--------------------------------------------------------------------------     



    !--------------------------------------------------------------------------     

    subroutine ReadRealVector_New(vector,                                     &
                                  EnterDataID,                                &
                                  flag,                                       &
                                  SearchType,                                 &
                                  text,                                       &
                                  Buffer_Line,                                &
                                  keyword,                                    &
                                  Default,                                    &
                                  ClientModule,                               & !Future use for interface
                                  CaseSensitive,                              &
                                  STAT)

        !Parameter-------------------------------------------------------------
        integer                                     :: EnterDataID
        character(LEN = *), optional, intent(IN )   :: keyword
        character(LEN = *), optional, intent(IN )   :: text
        character(LEN = *), optional, intent(IN )   :: ClientModule
        logical           , optional, intent(IN )   :: CaseSensitive

        integer, optional,       intent(IN )        :: Buffer_Line

        integer,                 intent(OUT)        :: flag
        integer, optional,       intent(OUT)        :: STAT
        integer, optional,       intent(IN )        :: SearchType

        real(4), dimension(:),   intent(OUT)        :: vector
        real(4), optional,       intent(IN )        :: Default


        !External--------------------------------------------------------------
        character(StringLength)                     :: ClientModuleOUT, DefaultOUT, ValueOUT
        logical                                     :: CaseSensitiveOUT
        integer                                     :: SearchTypeOUT
        integer                                     :: STAT_CALL
        integer                                     :: ready_          

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable
        integer :: vector_size    
        integer :: length, i, n_data          
                                        
        character (len=line_length) :: aux_str = "" 

        Type (T_Dataline) :: string
          
        !----------------------------------------------------------------------

        STAT_       = UNKNOWN_

        length      = null_int
        vector_size = null_int
        n_data      = null_int

        call Ready(EnterDataID, ready_)
        
cd0 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            
            vector = null_real
            flag   = 0


if9 :       if (present(SearchType)) then
                SearchTypeOut = SearchType
            else
                SearchTypeOut = FromFile_
            end if if9


            if (.not.present(keyword).and..not.present(Buffer_Line))          &
              stop "Error 01 ReadRealVector_New - ModuleEnter_Data"
            
            if (present(keyword).and.present(Buffer_Line))                    &
              stop "Error 02 ReadRealVector_New - ModuleEnter_Data"

            if (present(keyword)) then

                if (present (CaseSensitive)) then

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag, CaseSensitive)

                else

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag)

                endif

            elseif (present(Buffer_Line)) then

                 string%full_line     = Me%BufferLines(Buffer_Line)%full_line
                 string%delimiter_pos = 0
                 flag                 = 1

            endif

cd2 :       if (flag == 1) then
                vector_size = size(vector)

                read (string%full_line(string%delimiter_pos + 1 : line_length),"(A)") aux_str
                aux_str = adjustl (aux_str)
                length  = len_trim(aux_str)
                
                n_data  = 0

cd1 :           if (length > 0) then ! counts the number of data
                    n_data = 1
do1 :                   do i = 2, length
                            if((aux_str (i:i) == " ") .AND.                   &
                               (aux_str (i-1:i-1) /= " "))                    &
                                n_data = n_data + 1 
                        end do do1
                end if cd1

                read (aux_str(1 : length),*, err=100) (vector (i), i = 1, min (n_data, vector_size))
                flag = min (n_data, vector_size) 
            else
cd6 :           if (present(Default)) then
                    vector = Default

                    STAT_ = SUCCESS_
                else
cd4 :               if (present(text)) then
                        call WriteErrorMessage(keyword, text, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                              &
                            stop "Subroutine ReadRealVector_New; module ModuleEnterData. ERR01."
                    end if cd4

                    STAT_ = KEYWORD_NOT_FOUND_ERR_
                end if cd6
            end if cd2

cd3 :       if (vector_size .NE. n_data) then
                STAT_ = SIZE_ERR_
            else 
                STAT_ = SUCCESS_
            end if cd3

            if (present(keyword)) then
                
                if (present(ClientModule)) then
                    ClientModuleOUT = ClientModule
                else
                    ClientModuleOUT = "Unknown Client"
                endif

                if (present(Default)) then
                    write(DefaultOUT, fmt=*) Default
                else
                    DefaultOUT = "Do not have"
                endif

                if (present(CaseSensitive)) then
                    CaseSensitiveOUT = CaseSensitive
                else
                    CaseSensitiveOUT = .true.
                endif

                if (vector_size > 0) then
                    write (ValueOUT, fmt=*, err=10) vector (1:min (n_data, vector_size))
                else
                    ValueOUT = ' '
                endif
                
10              continue

                call LogKeyWord (keyword, SearchTypeOUT, ClientModuleOUT, DefaultOUT, CaseSensitiveOUT, ValueOUT)

            endif

        else 
            STAT_ = ready_
        end if cd0


100     continue


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ReadRealVector_New

    !--------------------------------------------------------------------------     



    !--------------------------------------------------------------------------     

    subroutine ReadReal_New(value,                                            &
                            EnterDataID,                                      &
                            flag,                                             &
                            SearchType,                                       &
                            text,                                             &
                            Buffer_Line,                                      &
                            keyword,                                          &
                            Default,                                          &
                            ClientModule,                                     & !Future use for interface
                            CaseSensitive,                                    &
                            STAT)
      
        !Parameter-------------------------------------------------------------
        integer                                     :: EnterDataID    
        character(LEN = *), optional, intent(IN )   :: keyword
        character(LEN = *), optional, intent(IN )   :: text
        character(LEN = *), optional, intent(IN )   :: ClientModule
        logical           , optional, intent(IN )   :: CaseSensitive
        integer, optional,       intent(IN )        :: Buffer_Line

        integer,                 intent(OUT)        :: flag 
        integer, optional,       intent(OUT)        :: STAT
        integer, optional,       intent(IN )        :: SearchType

        real(4),                 intent(OUT)        :: value
        real(4), optional,       intent(IN )        :: Default


        !External--------------------------------------------------------------
        character(StringLength)                     :: ClientModuleOUT, DefaultOUT, ValueOUT
        logical                                     :: CaseSensitiveOUT
        integer                                     :: SearchTypeOUT
        integer                                     :: STAT_CALL
        integer                                     :: ready_          

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable
        integer :: length          
                                        
        character (len=line_length) :: aux_str = ""      

        Type (T_Dataline) :: string
          
        !----------------------------------------------------------------------

        STAT_       = UNKNOWN_

        length      = null_int



        call Ready(EnterDataID, ready_)
        
cd0 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then



            value = null_real
            flag   = 0



if9 :       if (present(SearchType)) then
                SearchTypeOut = SearchType
            else
                SearchTypeOut = FromFile_
            end if if9


            if (.not.present(keyword).and..not.present(Buffer_Line))          &
              stop "GetData (ReadReal_New); ModuleEnter_Data. ERR01."
            
            if (present(keyword).and.present(Buffer_Line))                    &
              stop "GetData (ReadReal_New); ModuleEnter_Data. ERR02."

            if (present(keyword)) then

                if (present (CaseSensitive)) then

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag, CaseSensitive)

                else

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag)

                endif

            elseif (present(Buffer_Line)) then

                 string%full_line     = Me%BufferLines(Buffer_Line)%full_line
                 string%delimiter_pos = 0
                 flag                 = 1

            endif

cd2 :       if (flag == 1) then
                read (string%full_line(string%delimiter_pos + 1 : line_length),"(A)") aux_str
                aux_str = adjustl (aux_str)
                length  = len_trim(aux_str)
                

                read (aux_str(1 : length),*, IOSTAT = STAT_CALL, err=100) value
100             continue
cd3 :           if (STAT_CALL .EQ. SUCCESS_) then
                    STAT_ = SUCCESS_
                else
                    STAT_ = SIZE_ERR_
                end if cd3

            else
cd6 :           if (present(Default)) then
                    value = Default

                    STAT_ = SUCCESS_
                else
cd4 :               if (present(text)) then
                        call WriteErrorMessage(keyword, text, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                              &
                          stop "GetData (ReadReal_New); ModuleEnter_Data. ERR03."
                    end if cd4

                    STAT_ = SUCCESS_
                end if cd6
            end if cd2

            if (present(keyword)) then
                
                if (present(ClientModule)) then
                    ClientModuleOUT = ClientModule
                else
                    ClientModuleOUT = "Unknown Client"
                endif

                if (present(Default)) then
                    write(DefaultOUT, fmt=*) Default
                else
                    DefaultOUT = "Do not have"
                endif

                if (present(CaseSensitive)) then
                    CaseSensitiveOUT = CaseSensitive
                else
                    CaseSensitiveOUT = .true.
                endif

                write (ValueOUT, fmt=*) value
                
                call LogKeyWord (keyword, SearchTypeOUT, ClientModuleOUT, DefaultOUT, CaseSensitiveOUT, ValueOUT)

            endif


        else 
            STAT_ = ready_
        end if cd0


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ReadReal_New

    !--------------------------------------------------------------------------     



    !--------------------------------------------------------------------------     

    subroutine ReadDbleVector_New(vector,                                     &
                                  EnterDataID,                                &
                                  flag,                                       &
                                  SearchType,                                 &
                                  text,                                       &
                                  Buffer_Line,                                &
                                  keyword,                                    &
                                  Default,                                    &
                                  ClientModule,                               & !Future use for interface
                                  CaseSensitive,                              &
                                  STAT)
      
        !Parameter-------------------------------------------------------------
        integer                                     :: EnterDataID    
        character(LEN = *), optional, intent(IN )   :: keyword
        character(LEN = *), optional, intent(IN )   :: text
        character(LEN = *), optional, intent(IN )   :: ClientModule
        logical           , optional, intent(IN )   :: CaseSensitive

        integer, optional,       intent(IN ) :: Buffer_Line

        integer,                 intent(OUT) :: flag 
        integer, optional,       intent(OUT) :: STAT
        integer, optional,       intent(IN ) :: SearchType

        real(8), dimension(:),   intent(OUT) :: vector
        real(8), optional,       intent(IN ) :: Default


        !External--------------------------------------------------------------
        character(StringLength)                     :: ClientModuleOUT, DefaultOUT, ValueOUT
        logical                                     :: CaseSensitiveOUT
        integer                                     :: SearchTypeOUT
        integer                                     :: STAT_CALL
        integer                                     :: ready_          

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable
        integer :: vector_size      
        integer :: length, i, n_data            
                                        
        character (len=line_length) :: aux_str      = ""

        Type (T_Dataline) :: string
          
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        length      = null_int
        vector_size = null_int
        n_data      = null_int


        call Ready(EnterDataID, ready_)
        
cd0 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            vector = null_real
            flag   = 0


if9 :       if (present(SearchType)) then
                SearchTypeOut = SearchType
            else
                SearchTypeOut = FromFile_
            end if if9


            if (.not.present(keyword).and..not.present(Buffer_Line))          &
              stop "Error 01 ReadDbleVector_New - ModuleEnter_Data"
            
            if (present(keyword).and.present(Buffer_Line))                    &
              stop "Error 02 ReadDbleVector_New - ModuleEnter_Data"

            if (present(keyword)) then

                if (present (CaseSensitive)) then

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag, CaseSensitive)

                else

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag)

                endif

            elseif (present(Buffer_Line)) then

                 string%full_line     = Me%BufferLines(Buffer_Line)%full_line
                 string%delimiter_pos = 0
                 flag                 = 1

            endif

cd2 :       if (flag == 1) then
                vector_size = size(vector)

                read (string%full_line(string%delimiter_pos + 1 : line_length),"(A)") aux_str
                aux_str = adjustl (aux_str)
                length  = len_trim(aux_str)
                
                n_data  = 0

cd1 :           if (length > 0) then ! counts the number of data
                    n_data = 1
do1 :                   do i = 2, length
                            if((aux_str (i:i) == " ") .AND.                   &
                               (aux_str (i-1:i-1) /= " "))                    &
                                n_data = n_data + 1 
                        end do do1
                end if cd1

                read (aux_str(1 : length),*, err=100) (vector (i), i = 1, min (n_data, vector_size))
                flag = min (n_data, vector_size) 
            else
cd6 :           if (present(Default)) then
                    vector = Default

                    STAT_ = SUCCESS_
                else
cd4 :               if (present(text)) then
                        call WriteErrorMessage(keyword, text, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                              &
                            stop "Subroutine ReadDbleVector_New; module ModuleEnterData. ERR01."
                    end if cd4

                    STAT_ = KEYWORD_NOT_FOUND_ERR_
                end if cd6
            end if cd2

cd3 :       if (vector_size .NE. n_data) then
                STAT_ = SIZE_ERR_
            else 
                STAT_ = SUCCESS_
            end if cd3


            if (present(keyword)) then
                
                if (present(ClientModule)) then
                    ClientModuleOUT = ClientModule
                else
                    ClientModuleOUT = "Unknown Client"
                endif

                if (present(Default)) then
                    write(DefaultOUT, fmt=*) Default
                else
                    DefaultOUT = "Do not have"
                endif

                if (present(CaseSensitive)) then
                    CaseSensitiveOUT = CaseSensitive
                else
                    CaseSensitiveOUT = .true.
                endif

                if (vector_size > 0) then
                    write (ValueOUT, fmt=*, err=10) vector (1:min (n_data, vector_size))
                else
                    ValueOUT = ' '
                endif

10              continue
                
                call LogKeyWord (keyword, SearchTypeOUT, ClientModuleOUT, DefaultOUT, CaseSensitiveOUT, ValueOUT)

            endif

        else 
            STAT_ = ready_
        end if cd0


100     continue


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ReadDbleVector_New

    !--------------------------------------------------------------------------     



    !--------------------------------------------------------------------------     

    subroutine ReadTime(value,                                                &
                        EnterDataID,                                          &
                        flag,                                                 &
                        SearchType,                                           &
                        text,                                                 &
                        Buffer_Line,                                          &
                        keyword,                                              &
                        Default,                                              &
                        ClientModule,                                         & !Future use for interface
                        CaseSensitive,                                        &
                        STAT)
      
        !Parameter-------------------------------------------------------------
        integer                                     :: EnterDataID    
        character(LEN = *), optional, intent(IN )   :: keyword
        character(LEN = *), optional, intent(IN )   :: text
        character(LEN = *), optional, intent(IN )   :: ClientModule
        logical           , optional, intent(IN )   :: CaseSensitive

        integer, optional,       intent(IN ) :: Buffer_Line

        integer,                 intent(OUT) :: flag 
        integer, optional,       intent(OUT) :: STAT
        integer, optional,       intent(IN ) :: SearchType

        type(T_Time),            intent(OUT) :: value
        type(T_Time), optional,  intent(IN ) :: Default
     

        !External--------------------------------------------------------------
        character(len=256)                          :: ClientModuleOUT, DefaultOUT, ValueOUT
        logical                                     :: CaseSensitiveOUT
        integer                                     :: SearchTypeOUT
        integer                                     :: STAT_CALL
        integer                                     :: ready_          

        real    :: Year, Month, Day, Hour, Minute, Second

        !Local-----------------------------------------------------------------

        integer :: STAT_             
        integer :: length            
                                        
        character (len=line_length) :: aux_str = ""      

        Type (T_Dataline) :: string
          
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        flag   = 0

        call Ready(EnterDataID, ready_)
        
cd0 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then


if9 :       if (present(SearchType)) then
                SearchTypeOut = SearchType
            else
                SearchTypeOut = FromFile_
            end if if9


            if (.not.present(keyword).and..not.present(Buffer_Line))          &
                stop "Subroutine GetData (ReadTime); module ModuleEnterData. ERR01."
            
            if (present(keyword).and.present(Buffer_Line))                    &
                stop "Subroutine GetData (ReadTime); module ModuleEnterData. ERR02."

            if (present(keyword)) then

                if (present (CaseSensitive)) then

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag, CaseSensitive)

                else

                    call LookForLineInBlock(keyword, string, SearchTypeOUT, flag)

                endif

            elseif (present(Buffer_Line)) then

                 string%full_line     = Me%BufferLines(Buffer_Line)%full_line
                 string%delimiter_pos = 0
                 flag                 = 1

            endif

cd2 :       if (flag == 1) then
                read (string%full_line(string%delimiter_pos + 1 : line_length),"(A)") aux_str
                aux_str = adjustl (aux_str)
                length  = len_trim(aux_str)
                
                read (aux_str(1 : length),*, IOSTAT = STAT_CALL, err=100) Year, Month, Day, Hour, Minute, Second
100             continue
cd3 :           if (STAT_CALL .EQ. SUCCESS_) then
                    call SetDate(value, Year, Month, Day, Hour, Minute, Second)

                    STAT_ = SUCCESS_
                else
                    STAT_ = SIZE_ERR_
                end if cd3

            else
cd6 :           if (present(Default)) then
                    value = Default

                    STAT_ = SUCCESS_
                else
cd4 :               if (present(text)) then
                        call WriteErrorMessage(keyword, text, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                              &
                            stop "Subroutine GetData (ReadTime); module ModuleEnterData. ERR03."
                    end if cd4

                    STAT_ = SUCCESS_
                end if cd6
            end if cd2

            if (present(keyword)) then
                
                if (present(ClientModule)) then
                    ClientModuleOUT = ClientModule
                else
                    ClientModuleOUT = "Unknown Client"
                endif

                if (present(Default)) then
                    call ExtractDate (Default, Year, Month, Day, Hour, Minute, Second)
                    write(DefaultOUT, fmt=*) Year, Month, Day, Hour, Minute, Second
                else
                    DefaultOUT = "Do not have"
                endif

                if (present(CaseSensitive)) then
                    CaseSensitiveOUT = CaseSensitive
                else
                    CaseSensitiveOUT = .true.
                endif

                call ExtractDate (value, Year, Month, Day, Hour, Minute, Second)
                write (ValueOUT, fmt=*) Year, Month, Day, Hour, Minute, Second
                
                call LogKeyWord (keyword, SearchTypeOUT, ClientModuleOUT, DefaultOUT, CaseSensitiveOUT, ValueOUT)

            endif


        else 
            STAT_ = ready_
        end if cd0


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ReadTime

    !--------------------------------------------------------------------------     

    !--------------------------------------------------------------------------
    !This soubroutine extracts, from unit "Me%unit", the lines 
    ! between keywords "block_begin" and "block_end" 
    ! and stores them in Me%BufferLines 

    subroutine ExtractBlockFromBuffer(EnterDataID, ClientNumber,               &
                                                    block_begin,  block_end,   &
                                                    BlockFound,                &
                                                    FirstLine, LastLine, STAT) 
                                 
        !External--------------------------------------------------------------
        integer                                     :: EnterDataID    
        integer,            intent(INOUT)           :: ClientNumber
        integer, optional,  intent(OUT  )           :: STAT    
        integer, optional,  intent(OUT  )           :: FirstLine    
        integer, optional,  intent(OUT  )           :: LastLine    
        logical,            intent(OUT  )           :: BlockFound
        character(LEN = *), intent(IN   )           :: block_begin
        character(LEN = *), intent(IN   )           :: block_end

                     
        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT1    
        integer                                     :: FirstLine1, LastLine1    
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then
if8 :       if (.NOT. Me%BLOCK_LOCK) then                       !First block
                call Block_Lock

                call ExtractBlockFromBuffer1(block_begin,  block_end,  &
                                             BlockFound,               &
                                             FirstLine1, LastLine1, STAT1)

                ClientNumber = Me%BlockClientIDnumber 

                if (present(FirstLine)) FirstLine = FirstLine1
                if (present(LastLine )) LastLine  = LastLine1

                STAT_ = STAT1
            else                                                
if9 :           if (Me%BlockClientIDnumber .EQ. ClientNumber) then     !Second or higher block
                    call ExtractBlockFromBuffer1(block_begin,  block_end,  &
                                                 BlockFound,               &
                                                 FirstLine1, LastLine1, STAT1)

                    if (present(FirstLine)) FirstLine = FirstLine1
                    if (present(LastLine )) LastLine  = LastLine1

                    STAT_ = STAT1
                else
                    STAT_ = BLOCK_LOCK_ERR_     !ExtractBlockFromBuffer in use by other client
                end if if9
            end if if8
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ExtractBlockFromBuffer

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine ExtractBlockFromBuffer1(block_begin,  block_end,  &
                                       BlockFound,               &
                                       FirstLine1, LastLine1, STAT1) 
                                 
        !External--------------------------------------------------------------
           
        integer, intent(OUT) :: STAT1
        integer, intent(OUT) :: FirstLine1    
        integer, intent(OUT) :: LastLine1    
            
        logical,intent(OUT) :: BlockFound

        character(LEN = *), intent(IN ) :: block_begin
        character(LEN = *), intent(IN ) :: block_end

        !Local-----------------------------------------------------------------

        integer :: I, start

        logical :: FoundBegin, FoundEnd

        character(LEN = line_length) :: ReadKeyWord

        !----------------------------------------------------------------------

        start = max(Me%Block%EndBlock, 1)   !Verifyes if a previous block was found

        !Looks for block_begin 
        FoundBegin = .FALSE.
        FoundEnd   = .FALSE.
do1 :   do I = start, Me%BufferSize
            call ScanLine(Me%BufferLines(I))

            ReadKeyWord = Me%BufferLines(I)%full_line(1:Me%BufferLines(I)%delimiter_pos-1)
            if (trim(adjustl(ReadKeyWord)) .EQ. trim(adjustl(block_begin))) then
                Me%Block%BeginBlock = I
                FoundBegin = .TRUE.
                exit do1
            end if 
        end do do1


cd3 :   if (FoundBegin) then
do2 :       do I = Me%Block%BeginBlock+1, Me%BufferSize

                call ScanLine(Me%BufferLines(I))
                ReadKeyWord = Me%BufferLines(I)%full_line(1:Me%BufferLines(I)%delimiter_pos-1)

if8 :           if (trim(adjustl(ReadKeyWord)) .EQ. trim(adjustl(block_begin))) then
                    FoundEnd = .FALSE.
                    exit do2
                end if if8



cd2 :           if (trim(adjustl(ReadKeyWord)) .EQ. trim(adjustl(block_end))) then
                    Me%Block%EndBlock = I
                    FoundEnd = .TRUE.
                    exit do2
                end if cd2

            end do do2
        end if cd3


        FirstLine1 = Me%Block%BeginBlock
        LastLine1  = Me%Block%EndBlock


cd5 :   if      (       FoundBegin  .AND.        FoundEnd ) then
            BlockFound = .TRUE.
            STAT1      = SUCCESS_

        else if (       FoundBegin  .AND. (.NOT. FoundEnd)) then
            BlockFound = .FALSE.
            STAT1      = BLOCK_END_ERR_

        else if ((.NOT. FoundBegin) .AND. (.NOT. FoundEnd)) then
            BlockFound = .FALSE.
            STAT1      = SUCCESS_

            Me%Block%BeginBlock = null_int    !Ready to scan the file again
            Me%Block%EndBlock   = null_int    !Ready to scan the file again
        end if cd5


        STAT1 = SUCCESS_

    !--------------------------------------------------------------------------

    end subroutine ExtractBlockFromBuffer1

    !--------------------------------------------------------------------------

    subroutine ExtractBlockFromBlock(EnterDataID,  ClientNumber,              &
                                                   block_begin,  block_end,   &
                                                   BlockInBlockFound,         &
                                                   FirstLine, LastLine, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataID    
        integer,            intent(IN )             :: ClientNumber
        integer, optional,  intent(OUT)             :: FirstLine    
        integer, optional,  intent(OUT)             :: LastLine    
        integer, optional,  intent(OUT)             :: STAT    
        character(LEN = *), intent(IN )             :: block_begin
        character(LEN = *), intent(IN )             :: block_end
        logical,            intent(OUT)             :: BlockInBlockFound


        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_
        integer                                     :: I, start
        logical                                     :: FoundBegin, FoundEnd
        character(LEN = line_length)                :: ReadKeyWord

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

cd22 :      if ((Me%BLOCK_LOCK) .AND.                         &
                (Me%BlockClientIDnumber .EQ. ClientNumber)) then

                start = max(Me%BlockFromBlock%EndBlock,             & !Verifyes if a previous block was found
                            Me%Block%BeginBlock)   

                !Looks for block_begin 
                FoundBegin = .FALSE.
                FoundEnd   = .FALSE.
do1 :           do I = start+1, Me%Block%EndBlock-1

                    call ScanLine(Me%BufferLines(I))
                    ReadKeyWord = Me%BufferLines(I)%full_line(1:Me%BufferLines(I)%delimiter_pos-1)

                    if (trim(adjustl(ReadKeyWord)) .EQ. trim(adjustl(block_begin))) then
                        Me%BlockFromBlock%BeginBlock = I
                        FoundBegin = .TRUE.
                        exit do1
                    end if    
                end do do1


cd3 :           if (FoundBegin) then
do2 :               do I = Me%BlockFromBlock%BeginBlock+1, Me%Block%EndBlock-1

                        call ScanLine(Me%BufferLines(I))
                        ReadKeyWord = Me%BufferLines(I)%full_line(1:Me%BufferLines(I)%delimiter_pos-1)

if8 :                   if (trim(adjustl(ReadKeyWord)) .EQ. trim(adjustl(block_begin))) then
                            FoundEnd = .FALSE.
                            exit do2
                        end if if8



cd2 :                   if (trim(adjustl(ReadKeyWord)) .EQ. trim(adjustl(block_end))) then
                            Me%BlockFromBlock%EndBlock = I
                            FoundEnd = .TRUE.
                            exit do2
                        end if cd2

                    end do do2
                end if cd3


                if (present(FirstLine)) FirstLine = Me%BlockFromBlock%BeginBlock
                if (present(LastLine))  LastLine  = Me%BlockFromBlock%EndBlock


cd5 :           if      (       FoundBegin  .AND.        FoundEnd ) then
                    BlockInBlockFound = .TRUE.
                    STAT_             = SUCCESS_

                else if (       FoundBegin  .AND. (.NOT. FoundEnd)) then
                    BlockInBlockFound = .FALSE.
                    STAT_             = BLOCK_END_ERR_

                else if ((.NOT. FoundBegin) .AND. (.NOT. FoundEnd)) then
                    BlockInBlockFound = .FALSE.
                    STAT_             = SUCCESS_

                    Me%BlockFromBlock%BeginBlock = null_int   !Ready to scan the file again
                    Me%BlockFromBlock%EndBlock   = null_int   !Ready to scan the file again
                end if cd5


            STAT_ = SUCCESS_
            end if cd22
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ExtractBlockFromBlock

    !--------------------------------------------------------------------------



    !----------------------------------------------------------------------------

    subroutine ExtractBlockFromBlockFromBlock(EnterDataID,                      &
                                              ClientNumber,                     &
                                              block_begin, block_end,           &
                                              BlockInBlockInBlockFound,         &
                                              FirstLine, LastLine,              &
                                              STAT) 

        !Arguments---------------------------------------------------------------
        integer                                     :: EnterDataID
        integer,            intent(IN )             :: ClientNumber
        integer, optional,  intent(OUT)             :: FirstLine    
        integer, optional,  intent(OUT)             :: LastLine    
        integer, optional,  intent(OUT)             :: STAT    
                                                    
        character(LEN = *), intent(IN )             :: block_begin
        character(LEN = *), intent(IN )             :: block_end
                                                    
        logical,            intent(OUT)             :: BlockInBlockInBlockFound

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_
        integer                                     :: I, start
        logical                                     :: FoundBegin, FoundEnd
        character(LEN = line_length)                :: ReadKeyWord

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

cd22 :      if ((Me%BLOCK_LOCK) .AND.                           &
                (Me%BlockClientIDnumber .EQ. ClientNumber)) then

                start = max(Me%BlockFromBlockFromBlock%EndBlock,      & !Verifyes if a previous BlockFromBlock was found
                            Me%BlockFromBlock%BeginBlock)   

                !Looks for block_begin 
                FoundBegin = .FALSE.
                FoundEnd   = .FALSE.
do1 :           do I = start+1, Me%BlockFromBlock%EndBlock-1

                    call ScanLine(Me%BufferLines(I))
                    ReadKeyWord = Me%BufferLines(I)%full_line(1:Me%BufferLines(I)%delimiter_pos-1)

                    if (trim(adjustl(ReadKeyWord)) .EQ. trim(adjustl(block_begin))) then
                        Me%BlockFromBlockFromBlock%BeginBlock = I
                        FoundBegin = .TRUE.
                        exit do1
                    end if    
                end do do1


cd3 :           if (FoundBegin) then
do2 :               do I = Me%BlockFromBlockFromBlock%BeginBlock+1, Me%BlockFromBlock%EndBlock-1

                        call ScanLine(Me%BufferLines(I))
                        ReadKeyWord = Me%BufferLines(I)%full_line(1:Me%BufferLines(I)%delimiter_pos-1)

if8 :                   if (trim(adjustl(ReadKeyWord)) .EQ. trim(adjustl(block_begin))) then
                            FoundEnd = .FALSE.
                            exit do2
                        end if if8



cd2 :                   if (trim(adjustl(ReadKeyWord)) .EQ. trim(adjustl(block_end))) then
                            Me%BlockFromBlockFromBlock%EndBlock = I
                            FoundEnd = .TRUE.
                            exit do2
                        end if cd2

                    end do do2
                end if cd3


                if (present(FirstLine)) FirstLine = Me%BlockFromBlockFromBlock%BeginBlock
                if (present(LastLine )) LastLine  = Me%BlockFromBlockFromBlock%EndBlock


cd5 :           if      (       FoundBegin  .AND.        FoundEnd ) then
                    BlockInBlockInBlockFound = .TRUE.
                    STAT_                    = SUCCESS_

                else if (       FoundBegin  .AND. (.NOT. FoundEnd)) then
                    BlockInBlockInBlockFound = .FALSE.
                    STAT_                    = BLOCK_END_ERR_

                else if ((.NOT. FoundBegin) .AND. (.NOT. FoundEnd)) then
                    BlockInBlockInBlockFound = .FALSE.
                    STAT_                    = SUCCESS_

                    Me%BlockFromBlockFromBlock%BeginBlock = null_int   !Ready to scan the file again
                    Me%BlockFromBlockFromBlock%EndBlock   = null_int   !Ready to scan the file again
                end if cd5


            STAT_ = SUCCESS_
            end if cd22
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ExtractBlockFromBlockFromBlock


    !--------------------------------------------------------------------------

    subroutine GetOutPutTime(EnterDataID, CurrentTime, EndTime, Keyword, SearchType,    &
                             OutPutsTime, OutPutsON, OutPutsNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: EnterDataID
        type (T_Time),     intent (IN)              :: CurrentTime, EndTime
        character(Len= *), intent (IN)              :: Keyword
        integer,           intent (IN)              :: SearchType
        type(T_Time), dimension(:), pointer         :: OutPutsTime
        logical,           intent(OUT)              :: OutPutsON            
        integer, optional, intent(OUT)              :: OutPutsNumber
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        type (T_Time)                               :: AuxTime
        integer                                     :: ready_          
        integer                                     :: STAT_

        STAT_ = UNKNOWN_


        call Ready(EnterDataID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call SetDate(AuxTime, 0, 0, 0, 0, 0, 0)

            if (CurrentTime < AuxTime .or. EndTime < AuxTime .or. CurrentTime > EndTime) then

                STAT_ = UNKNOWN_

            else

                if (present(OutPutsNumber)) then 
        
                    call OutPutTimeInternal   (EnterDataID, CurrentTime, EndTime,            &
                                               Keyword, SearchType,                          &
                                               OutPutsTime, OutPutsON, OutPutsNumber)
                else

                    call OutPutTimeInternal   (EnterDataID, CurrentTime, EndTime,            &
                                               Keyword, SearchType,                          &
                                               OutPutsTime, OutPutsON)

                endif


                STAT_ = SUCCESS_

            endif

        else
         
            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetOutPutTime

    !--------------------------------------------------------------------------

    subroutine OutPutTimeInternal(EnterDataID, CurrentTime, EndTime,                     &
                                  Keyword, SearchType, OutPutsTime,                      &
                                  OutPutsON, OutPutsNumber)

        !Arguments---------------------------------------------------------------
        integer                                     :: EnterDataID
        type (T_Time),     intent (IN)              :: CurrentTime, EndTime
        character(Len= *), intent (IN)              :: Keyword
        integer,           intent (IN)              :: SearchType
        type(T_Time), dimension(:), pointer         :: OutPutsTime
        logical,           intent (OUT)             :: OutPutsON            
        integer, optional, intent (OUT)             :: OutPutsNumber

        !Local-----------------------------------------------------------------
        integer                                     :: Status
        integer                                     :: iflag, ReadValuesNumber
        integer                                     :: i, ExtraOutPuts, OutPutsNumber_ =0
        real                                        :: ExtraTime_Seconds
        real, allocatable, dimension (:)            :: AuxDT
        type(T_Time), dimension(:), pointer         :: AuxTime 


        !----------------------------------------------------------------------

        ! Output Time 
        ! The first (iflag-1) values are consider specific outputs in time
        ! this first values must be in ascending order,  
        ! the values are given in seconds and the TimeI is consider to be the zero.
        ! The (iflag) value is consider to be the interval between values 
        ! from the (iglag-1) output forward, example: 

        ! TimeI = 1998 1 1  0 0 0  
        ! TimeE = 1998 1 1 12 0 0
        ! OUTPUT_TIME : 0. 7200.  14400. 14400.

        ! The result of this are the follow outputs:
        ! 1998 1 1 0 0 0  
        ! 1998 1 1 2 0 0  
        ! 1998 1 1 4 0 0  
        ! 1998 1 1 8 0 0  
        ! 1998 1 1 12 0 0  


        !Allocate auxliar variable
        allocate(AuxDT(240), STAT = Status)

        if (Status /= SUCCESS_)                                                          &
            call SetError(FATAL_, INTERNAL_, "OutPutTimeInternal - EnterData - ERR01") 

        

        AuxDT = FillValueReal

        call GetData(AuxDT, EnterDataID, iflag,                                          &
                     Keyword    = Keyword,                                               &
                     SearchType = SearchType,                                            &
                     STAT       = Status)

        if (Status /= SUCCESS_)   then

            !By default 240 values are read so this error always exist.
            if (Status /= SIZE_ERR_) then 

                call SetError(FATAL_, INTERNAL_, "OutPutTimeInternal - EnterData - ERR02") 

            endif

        endif

        OutPutsOn = .true.

        ! If the OUTPUT_TIME is not defined then by default is consider
        ! that this property don"t have outputs
        if (iflag==0) OutPutsOn = .false.


        ! If iflag=1 it can means :
        ! 1 - only one value is read  
        ! 2 - the keyword was find but don"t have any values AuxDT(1)=FillValueReal
cd1 :   if (iflag==1      ) then
cd2 :   if (AuxDT(1) > 0) then 
            ! When AuxDT(1)> 0 then is consider that the interval between outputs
            ! is equal to this value and the first output equal to timeI
            AuxDT(2)       = AuxDT(1)
            AuxDT(1)       = 0
            ReadValuesNumber = 2
        else 
            !When AuxDT(1) <= 0 is consider that this property don"t have any OutPuts
            OutPutsOn = .false.
        end if cd2
        end if cd1

cd3:    if (OutPutsOn) then
            ! iflag returns also the number of values read                   
            ReadValuesNumber = iflag

            ! We need to now the interval between the last specific output and the end of the simulation (TimeE).
            ! This information is important to compute the number of outputs after the last specific 
            ! output 


            ExtraTime_Seconds = EndTime - (CurrentTime + AuxDT(ReadValuesNumber - 1))

            if (ExtraTime_Seconds < 0)                                                   &
                call SetError(FATAL_, INTERNAL_, "OutPutTimeInternal - EnterData - ERR03")        

            ! ExtraOutPuts is equal to the total number of outputs after the last specific output         
            ExtraOutPuts = int(ExtraTime_Seconds/AuxDT(ReadValuesNumber)) 
            ! The total number of output is the number of vaulea read less one plus the 
            ! "ExtraOutPuts"
            OutPutsNumber_ = ReadValuesNumber + ExtraOutPuts + 1

            allocate (AuxTime(1:OutPutsNumber_), STAT = Status)

            if (Status /= SUCCESS_)                                                      &
                call SetError(FATAL_, INTERNAL_, "OutPutTimeInternal - EnterData - ERR04") 

do2 :       do i=1,ReadValuesNumber-1
                !AppHydrodynamic%OutPut%OutTime(i) = AppHydrodynamic%CurrentTime + AuxDT(i)
                AuxTime(i) = CurrentTime + AuxDT(i)
            end do do2

do1 :       do i=ReadValuesNumber,OutPutsNumber_
                !AppHydrodynamic%OutPut%OutTime(i) = AppHydrodynamic%OutPut%OutTime(i-1) + AuxDT(ReadValuesNumber)
                AuxTime(i) = AuxTime(i-1) +  AuxDT(ReadValuesNumber)

                if (AuxTime(i) >= EndTime) then

                    AuxTime(i) = EndTime

                    OutPutsNumber_ = i
                    exit

                endif

            end do do1

            allocate(OutPutsTime(1:OutPutsNumber_), STAT = Status)

            if (Status /= SUCCESS_)                                                      &
                call SetError(FATAL_, INTERNAL_, "OutPutTimeInternal - EnterData - ERR05") 


            OutPutsTime(1:OutPutsNumber_) = AuxTime(1:OutPutsNumber_)

            deallocate (AuxTime, STAT = Status)

            if (Status /= SUCCESS_)                                                      &
                call SetError(FATAL_, INTERNAL_, "OutPutTimeInternal - EnterData - ERR06") 

        end if cd3

        !Deallocate auxliar variable
        Deallocate(AuxDT, STAT = Status)

        if (Status /= SUCCESS_)                                                          &
                call SetError(FATAL_, INTERNAL_, "OutPutTimeInternal - EnterData - ERR07") 


        if (present(OutPutsNumber)) OutPutsNumber = OutPutsNumber_

    end subroutine OutPutTimeInternal


    !--------------------------------------------------------------------------

    subroutine GetOutPutTimeWindows(EnterDataID, CurrentTime, EndTime,             &
                                    OutPutWindows, OutPutWindowsON, WindowsNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: EnterDataID
        type (T_Time),     intent (IN)              :: CurrentTime, EndTime
        type(T_OutPutTime), dimension(:), pointer   :: OutPutWindows
        logical,           intent(OUT)              :: OutPutWindowsON            
        integer,           intent(OUT)              :: WindowsNumber
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_, ClientNumber          
        integer                                     :: STAT_, STAT_CALL
        logical                                     :: FoundBlock

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call ExtractBlockFromBuffer(Me%InstanceID, ClientNumber,                    &
                                        '<beginoutput>', '<endoutput>',                 &
                                        FoundBlock, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetOutPutTimeWindows - EnterData - ERR10'

if1:        if (FoundBlock) then

                call CountNumberOfWindows(WindowsNumber, ClientNumber)
                
                allocate(OutPutWindows   (WindowsNumber))
                
                if (WindowsNumber > 0) OutPutWindowsON = .true.
                
                call ReadOutPutWindows   (CurrentTime, EndTime, ClientNumber, OutPutWindows)
                
            endif if1
            
            call RewindBuffer(Me%InstanceID, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetOutPutTimeWindows - EnterData - ERR20'   
            
            call Block_Unlock(Me%InstanceID, ClientNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetOutPutTimeWindows - EnterData - ERR30'
                                                 
            STAT_ = SUCCESS_

        else cd1
         
            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetOutPutTimeWindows

    !--------------------------------------------------------------------------

    subroutine CountNumberOfWindows(WindowsNumber, ClientNumber)
    
    
        !Arguments-------------------------------------------------------------
        integer                     :: WindowsNumber
        
        !Local-----------------------------------------------------------------
        integer                     :: ClientNumber, STAT_CALL
        logical                     :: FoundWindow

        !Begin-----------------------------------------------------------------
        
        WindowsNumber = 0
                                                                                            
d1:     do 
    
            call ExtractBlockFromBlock(Me%InstanceID, ClientNumber,                     &
                                       '<<beginoutwindow>>', '<<endoutwindow>>',        &
                                       FoundWindow, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CountNumberOfWindows - EnterData - ERR10'

if1:        if (FoundWindow) then
                             
                WindowsNumber = WindowsNumber + 1                 
                
            else  if1
            
                call RewindBlock(Me%InstanceID, ClientNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CountNumberOfWindows - EnterData - ERR20'
                
                exit

            endif if1
                        
        enddo d1
        
   
    end subroutine CountNumberOfWindows
    
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    subroutine ReadOutPutWindows(CurrentTime, EndTime, ClientNumber, OutPutWindows)
    
    
        !Arguments-------------------------------------------------------------
        type (T_Time),     intent (IN)              :: CurrentTime, EndTime
        integer,           intent (IN)              :: ClientNumber
        type(T_OutPutTime), dimension(:), pointer   :: OutPutWindows
                
        !Local-----------------------------------------------------------------
        integer, dimension(2)                       :: Aux
        integer                                     :: STAT_CALL, WN, iflag
        logical                                     :: FoundWindow

        !Begin-----------------------------------------------------------------
        
        WN = 0
                                                                                            
d1:     do
    
            call ExtractBlockFromBlock(Me%InstanceID, ClientNumber,                     &
                                       '<<beginoutwindow>>', '<<endoutwindow>>',        &
                                       FoundWindow, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutPutWindows - EnterData - ERR10'

if1:        if (FoundWindow) then

                WN = WN + 1
                                                                
                call GetOutPutTime(Me%InstanceID,                                       &
                    CurrentTime   = CurrentTime,                                        &
                    EndTime       = EndTime,                                            &
                    keyword       = 'OUTPUT_TIME_W',                                    &
                    SearchType    = FromBlockinBlock,                                   &
                    OutPutsTime   = OutPutWindows(WN)%OutTime,                          &
                    OutPutsOn     = OutPutWindows(WN)%ON,                               &
                    OutPutsNumber = OutPutWindows(WN)%Number,                           &
                    STAT          = STAT_CALL)

                if (STAT_CALL /= SUCCESS_)                                              &
                    call SetError(FATAL_, INTERNAL_, 'ReadOutPutWindows - EnterData - ERR20') 


                call GetDAta    (Aux,                                                   &
                                 Me%InstanceID, iflag,                                  &
                                 Keyword        = 'KLB_KUB_W',                          &
                                 SearchType     = FromBlockinBlock,                     &
                                 ClientModule   = 'ModuleEnterData',                    &
                                 STAT           = STAT_CALL)
                                 
                if (STAT_CALL /= SUCCESS_)                                              &
                    call SetError(FATAL_, INTERNAL_, 'ReadOutPutWindows - EnterData - ERR30') 

                if (iflag == 0)                                                         &
                    call SetError(FATAL_, INTERNAL_, 'ReadOutPutWindows - EnterData - ERR40') 
                    
                OutPutWindows(WN)%KLB = Aux(1)
                OutPutWindows(WN)%KUB = Aux(2)

                call GetDAta    (Aux,                                                   &
                                 Me%InstanceID, iflag,                                  &
                                 Keyword        = 'ILB_IUB_W',                          &
                                 SearchType     = FromBlockinBlock,                     &
                                 ClientModule   = 'ModuleEnterData',                    &
                                 STAT           = STAT_CALL)
                                 
                if (STAT_CALL /= SUCCESS_)                                              &
                    call SetError(FATAL_, INTERNAL_, 'ReadOutPutWindows - EnterData - ERR50') 

                if (iflag == 0)                                                         &
                    call SetError(FATAL_, INTERNAL_, 'ReadOutPutWindows - EnterData - ERR60') 

                OutPutWindows(WN)%ILB = Aux(1)
                OutPutWindows(WN)%IUB = Aux(2)

                call GetDAta    (Aux,                                                   &
                                 Me%InstanceID, iflag,                                  &
                                 Keyword        = 'JLB_JUB_W',                          &
                                 SearchType     = FromBlockinBlock,                     &
                                 ClientModule   = 'ModuleEnterData',                    &
                                 STAT           = STAT_CALL)
                                 
                if (STAT_CALL /= SUCCESS_)                                              &
                    call SetError(FATAL_, INTERNAL_, 'ReadOutPutWindows - EnterData - ERR70') 

                if (iflag == 0)                                                         &
                    call SetError(FATAL_, INTERNAL_, 'ReadOutPutWindows - EnterData - ERR80') 

                OutPutWindows(WN)%JLB = Aux(1)
                OutPutWindows(WN)%JUB = Aux(2)

                
            else  if1
            
                exit

            endif if1
                        
        enddo d1

    
    end subroutine ReadOutPutWindows
    
    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillEnterData(EnterDataID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_
        integer                                     :: STAT_CALL
        integer                                     :: nUsers

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)

cd1 :   if (ready_ .NE. OFF_ERR_) then
            

            nUsers = DeassociateInstance(mENTERDATA_,  Me%InstanceID)

            if (nUsers == 0) then

                !Deallocates buffer
                if (associated(Me%BufferLines)) then
                    deallocate(Me%BufferLines, STAT = STAT_CALL)
                end if

                !Deallocates Instance
                call DeallocateInstance

                EnterDataID = 0


                STAT_ = SUCCESS_

            end if

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine KillEnterData

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_EnterData), pointer             :: AuxEnterData
        type (T_EnterData), pointer             :: PreviousEnterData

        !Updates pointers
        if (Me%InstanceID == FirstEnterData%InstanceID) then
            FirstEnterData => FirstEnterData%Next
        else
            PreviousEnterData => FirstEnterData
            AuxEnterData      => FirstEnterData%Next
            do while (AuxEnterData%InstanceID /= Me%InstanceID)
                PreviousEnterData => AuxEnterData
                AuxEnterData      => AuxEnterData%Next
            enddo

            !Now update linked list
            PreviousEnterData%Next => AuxEnterData%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 
            
    end subroutine DeallocateInstance


    !--------------------------------------------------------------------------     


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine Ready (EnterDataID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (EnterDataID > 0) then
            call LocateObjEnterData (EnterDataID)
            ready_ = VerifyReadLock (mENTERDATA_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1


        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjEnterData (EnterDataID)

        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataID

        !Local-----------------------------------------------------------------

        Me => FirstEnterData
        do while (associated (Me))
            if (Me%InstanceID == EnterDataID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                                          &
            stop 'ModuleEnterData - LocateObjEnterData - ERR01'

    end subroutine LocateObjEnterData

    !--------------------------------------------------------------------------

    subroutine Block_Lock

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer :: number
        real    :: x

       !-----------------------------------------------------------------------

cd1 :   if (.NOT. Me%BLOCK_LOCK) then     
            Me%BLOCK_LOCK = ACTIVE

            call RANDOM_NUMBER(x)
            number = int(100000.0 * x)

            Me%BlockClientIDnumber = number
        else
            stop "Subroutine Block_Lock; module ModuleEnterData. ERR01."
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Block_Lock

    !--------------------------------------------------------------------------

    subroutine Block_Unlock(EnterDataID, ClientNumber, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataID
        integer,           intent(IN )              :: ClientNumber
        integer, optional, intent(OUT)              :: STAT    

        !Local-----------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(EnterDataID, ready_)
        
cd1 :   if ((ready_ .NE. OFF_ERR_       )) then

cd2 :       if (Me%BLOCK_LOCK) then
cd3 :           if (ClientNumber .NE. Me%BlockClientIDnumber) then

                    write(*,*) 
                    write(*,*) "Client ID mismatch."
                    stop       "Subroutine Block_Unlock; module ModuleEnterData. ERR01."

                else

                    Me%BlockClientIDnumber = null_int
                    Me%BLOCK_LOCK          = IDLE
                    STAT_                                  = SUCCESS_
         
                end if cd3
            else

                STAT_ = BLOCK_UNLOCK_ERR_

            end if cd2
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine Block_Unlock

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !ASCII OUTPUT - ASCII OUTPUT - ASCII OUTPUT - ASCII OUTPUT - ASCII OUTPUT - 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine WriteDataLineLogical(UnitID, KeyWord, LogicalValue)

        !Arguments-------------------------------------------------------------
        integer                     :: UnitID
        character(len=*)            :: KeyWord
        logical                     :: LogicalValue

        !Local-----------------------------------------------------------------
        character(len=25)           :: KeyWord_Delimiter
        character(len=27)           :: Line

        if (len_trim(adjustl(KeyWord)) < 25) then

            KeyWord_Delimiter        = trim(adjustl(KeyWord))
            KeyWord_Delimiter(25:25) = ":"
            if (LogicalValue) then
                Line = KeyWord_Delimiter//" 1"
            else
                Line = KeyWord_Delimiter//" 0"
            endif
    
            write(UnitID, "(A)") trim(Line)
         

        else
            if (LogicalValue) then
                write(UnitID, *)trim(adjustl(KeyWord))//" : 1"
            else
                write(UnitID, *)trim(adjustl(KeyWord))//" : 0"
            endif
        endif

    end subroutine WriteDataLineLogical

    !--------------------------------------------------------------------------

    !WriteDataLineInteger
    !   Writes an integer value
    !
    subroutine WriteDataLineInteger(UnitID, KeyWord, IntegerValue)

        !Arguments-------------------------------------------------------------
        integer                     :: UnitID
        character(len=*)            :: KeyWord
        integer                     :: IntegerValue

        !Local-----------------------------------------------------------------
        character(len=25)           :: KeyWord_Delimiter = " "
        character(len=10)           :: KeyWordValue      = " "
        character(len=36)           :: Line

        write(KeyWordValue, "(i10)")IntegerValue

        if (len_trim(adjustl(KeyWord)) < 25) then

            KeyWord_Delimiter        = trim(adjustl(KeyWord))
            KeyWord_Delimiter(25:25) = ":"
        
            Line = KeyWord_Delimiter//" "//trim(adjustl(KeyWordValue))
            write(UnitID, "(A)") trim(Line)
         
        else

            write(UnitID, *)trim(adjustl(KeyWord))//" : "//trim(adjustl(KeyWordValue))

        endif


    end subroutine WriteDataLineInteger

    !--------------------------------------------------------------------------

    !WriteDataLineReal
    !   Writes a real value
    !
    subroutine WriteDataLineReal(UnitID, KeyWord, RealValue, Form)

        !Arguments-------------------------------------------------------------
        integer                     :: UnitID
        character(len=*)            :: KeyWord
        real                        :: RealValue
        character(len=*), optional  :: Form

        !Local-----------------------------------------------------------------
        character(len=25)           :: KeyWord_Delimiter = " "
        character(len=15)           :: KeyWordValue      = " "
        character(len=41)           :: Line

        if (present (Form)) then
            write(KeyWordValue, Form)RealValue
        else
            KeyWordValue =  FormReal (RealValue)
        endif

        if (len_trim(adjustl(KeyWord)) < 25) then

            KeyWord_Delimiter        = trim(adjustl(KeyWord))
            KeyWord_Delimiter(25:25) = ":"

            Line = KeyWord_Delimiter//" "//trim(adjustl(KeyWordValue))
            write(UnitID, "(A)") trim(Line)

        else

            write(UnitID, *)trim(adjustl(KeyWord))//" : "//trim(adjustl(KeyWordValue))

        endif

    end subroutine WriteDataLineReal

    !--------------------------------------------------------------------------

     !WriteDataLineString
    !   Writes a string value
    !
    subroutine WriteDataLineString(UnitID, KeyWord, StringValue)

        !Arguments-------------------------------------------------------------
        integer                     :: UnitID
        character(len=*)            :: KeyWord
        character(len=*)            :: StringValue

        !Local-----------------------------------------------------------------
        character(len=25)           :: KeyWord_Delimiter = " "
        character(len=255)          :: Line

        if (len_trim(StringValue) == 0) StringValue = "No data"

        if (len_trim(adjustl(KeyWord)) < 25) then

            KeyWord_Delimiter        = trim(adjustl(KeyWord))
            KeyWord_Delimiter(25:25) = ":"

            Line = KeyWord_Delimiter//" "//trim(adjustl(StringValue))

            write(UnitID, "(A)") trim(Line)

        else

            write(UnitID, *)trim(adjustl(KeyWord))//" : "//trim(adjustl(StringValue))

        endif

    end subroutine WriteDataLineString

    !WriteDataLineBlock
    !   Writes a simple line (Blockend / Blockbegin)
    !
    subroutine WriteDataLineBlock(UnitID, Block)

        !Arguments-------------------------------------------------------------
        integer                     :: UnitID
        character(len=*)            :: Block

        !Local-----------------------------------------------------------------

        write(UnitID, "(A)") trim(Block)


    end subroutine WriteDataLineBlock


    !WriteDataLineRealVector
    !   Writes a real vector
    !
    subroutine WriteDataLineRealVector(UnitID, KeyWord, nValues, VectorValues)

        !Arguments-------------------------------------------------------------
        integer                     :: UnitID
        character(len=*)            :: KeyWord
        integer                     :: nValues
        real, dimension(:), pointer :: VectorValues

        !Local-----------------------------------------------------------------
        character(len=25)           :: KeyWord_Delimiter = " "
        character(len=255)          :: Line              = " "
        integer                     :: i
        character(10), allocatable, dimension(:) :: CharVectorValues


        allocate(CharVectorValues(nValues))
        do i = 1, nValues
            CharVectorValues(i) = FormReal (VectorValues(i))
        enddo

        Line = " "//trim(adjustl(CharVectorValues(1)))
        do i = 2, nValues
            Line = trim(Line)//" "//trim(adjustl(CharVectorValues(i)))
        enddo

        if (len_trim(adjustl(KeyWord)) < 25) then

            KeyWord_Delimiter        = trim(adjustl(KeyWord))
            KeyWord_Delimiter(25:25) = ":"

            write(UnitID, "(A)")KeyWord_Delimiter//" "//trim(Line)

        else

            write(UnitID, *)trim(adjustl(KeyWord))//" : "//trim(Line)

        endif
    
    end subroutine WriteDataLineRealVector

    !--------------------------------------------------------------------------

    !WriteDataLineVectorInt
    !   Writes a real vector
    !
    subroutine WriteDataLineIntVector(UnitID, KeyWord, nValues, VectorValues)

        !Arguments-------------------------------------------------------------
        integer                                     :: UnitID
        character(len=*)                            :: KeyWord
        integer                                     :: nValues
        integer, dimension(:), pointer              :: VectorValues

        !Local-----------------------------------------------------------------
        character(len=25)           :: KeyWord_Delimiter = " "
        character(len=255)          :: Line
        integer                     :: i
        character(10), allocatable, dimension(:) :: CharVectorValues


        allocate(CharVectorValues(nValues))
        do i = 1, nValues
            write(CharVectorValues(i), "(i10)")VectorValues(i)
        enddo

        do i = 1, nValues
            Line((i-1)*10+1:i*10) = trim(adjustl(CharVectorValues(i)))
        enddo

          if (len_trim(adjustl(KeyWord)) < 25) then

            KeyWord_Delimiter        = trim(adjustl(KeyWord))
            KeyWord_Delimiter(25:25) = ":"

            write(UnitID, "(A)")KeyWord_Delimiter//" "//trim(Line)

        else

            write(UnitID, *)trim(adjustl(KeyWord))//" : "//trim(Line)

        endif
    
    end subroutine WriteDataLineIntVector

    !--------------------------------------------------------------------------

    !WriteDataLineTime
    !   Writes a time value (YYYY MM DD HH MM SS)
    !
    subroutine WriteDataLineTime(UnitID, KeyWord, Time)

        !Arguments-------------------------------------------------------------
        integer                     :: UnitID
        character(len=*)            :: KeyWord
        type (T_Time)               :: Time

        !Local-----------------------------------------------------------------
        character(len=25)           :: KeyWord_Delimiter = " "
        character(len=200)          :: KeyWordValue      = " "
        character(len=255)          :: Line
        
        real                        :: Year, Month, Day, Hour, Minute, Second

        call ExtractDate(Time, Year, Month, Day, Hour, Minute, Second)            

        write(KeyWordValue, fmt=10)Year, Month, Day, Hour, Minute, Second
        10 format(f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)
 
        if (len_trim(adjustl(KeyWord)) < 25) then

            KeyWord_Delimiter        = trim(adjustl(KeyWord))
            KeyWord_Delimiter(25:25) = ":"

            Line = KeyWord_Delimiter//" "//trim(adjustl(KeyWordValue))
            write(UnitID, fmt=*)trim(Line)
!         20 format(a<len_trim(Line)>)

        else

            write(UnitID, *)trim(adjustl(KeyWord))//" : "//trim(adjustl(KeyWordValue))

        endif
    
    end subroutine WriteDataLineTime

    !--------------------------------------------------------------------------

    subroutine WriteDataLineOutputTime(UnitID, KeyWord, Time1, Time2)

        !Arguments-------------------------------------------------------------
        integer                     :: UnitID
        character(len=*)            :: KeyWord
        integer                     :: Time1, Time2

        !Local-----------------------------------------------------------------
        character(len=25)           :: KeyWord_Delimiter = " "
        character(len=10)           :: KeyWordValue1     = " "
        character(len=10)           :: KeyWordValue2     = " "
        character(len=36)           :: Line

        write(KeyWordValue1, "(i10)")Time1
        write(KeyWordValue2, "(i10)")Time2

        if (len_trim(adjustl(KeyWord)) < 25) then

            KeyWord_Delimiter        = trim(adjustl(KeyWord))
            KeyWord_Delimiter(25:25) = ":"
        
            Line = KeyWord_Delimiter//" "//trim(adjustl(KeyWordValue1))//" "//trim(adjustl(KeyWordValue2))
            write(UnitID, *)trim(Line)
!         10 format(a<len_trim(Line)>)

        else

            write(UnitID, *)trim(adjustl(KeyWord))//" : "//trim(adjustl(KeyWordValue1))//" "//trim(adjustl(KeyWordValue2))

        endif


    end subroutine WriteDataLineOutputTime

    !--------------------------------------------------------------------------

    character(len=StringLength) function FormReal(RealValue)

        !Arguments-------------------------------------------------------------
        real, intent(in)                            :: RealValue

        !Local-----------------------------------------------------------------
        character(len=8)                            :: AuxString
        character(len=StringLength)                 :: AuxString2
        character(len=6)                            :: FormatStatment
        logical                                     :: HaveFormat

           !-123456.  f8.0
           !-12345.6  f8.1
           !-1234.56  f8.2
           !-123.456  f8.3
           !-12.4567  f8.4
           !-1.23456  f8.5
           !-.123456  f8.6
           !-.000012  f8.6
        HaveFormat = .false.
        if     (abs(RealValue) < 1000000. .and. abs(RealValue) >= 100000.) then
            
            FormatStatment = "(f8.0)"
            HaveFormat     = .true.

        elseif (abs(RealValue) < 100000.  .and. abs(RealValue) >= 10000. ) then

            FormatStatment = "(f8.1)"
            HaveFormat     = .true.

        elseif (abs(RealValue) < 10000.   .and. abs(RealValue) >= 1000.  ) then

            FormatStatment = "(f8.1)"
            HaveFormat     = .true.

        elseif (abs(RealValue) < 1000.    .and. abs(RealValue) >= 100.   ) then

            FormatStatment = "(f7.2)"
            HaveFormat     = .true.

        elseif (abs(RealValue) < 100.     .and. abs(RealValue) >= 10.    ) then

            FormatStatment = "(f6.2)"
            HaveFormat     = .true.

        elseif (abs(RealValue) < 10.      .and. abs(RealValue) > 1.      ) then

            FormatStatment = "(f5.2)"
            HaveFormat     = .true.

        elseif (abs(RealValue) <= 1.       .and. abs(RealValue) >= 0.1    ) then

            FormatStatment = "(f5.3)"
            HaveFormat     = .true.

        elseif (abs(RealValue) < 0.1       .and. abs(RealValue) >= 0.01   ) then

            FormatStatment = "(f7.3)"
            HaveFormat     = .true.

        elseif (abs(RealValue) < 0.01) then

            FormatStatment = "(f6.3)"
            HaveFormat     = .true.

        endif

        if (RealValue== -1.00) FormatStatment = "(f5.2)"

        if (HaveFormat) then
            write(AuxString, fmt=FormatStatment)RealValue
            FormReal = trim(adjustl(AuxString))
            return
        else

            write(AuxString2, *)RealValue
            FormReal = trim(adjustl(AuxString2))
            return
        endif


    end function FormReal


end module ModuleEnterData

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------


