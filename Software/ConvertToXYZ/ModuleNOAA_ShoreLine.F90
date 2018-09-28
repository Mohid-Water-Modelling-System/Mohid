!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : ConvertToXYZ
! MODULE        : NOAA_ShoreLine
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2004
! REVISION      : Luis Fernandes
! DESCRIPTION   : Module to convert NOAA ShoreLines file into MOHID GIS Line format 
!                   - To download NOAA ShoreLines files go to 
!                     http://rimmer.ngdc.noaa.gov/mgg/coast/getcoast.html
!
!------------------------------------------------------------------------------

!   INPUT_FILENAME              : char              -           !Path to input file to convert
!   OUTPUT_FOLDER               : char              -           !Path to folder where files are generated
!   OUTPUT_NAME                 : char              -           !Name to give to output file(s)
!   WRITE_AS_POLYGONS           : char              -           !Write as polygons. If false lines are written
!   MAX_BUFFER_SIZE             : int           [1000000]       !Aprox. output file size in bytes

module ModuleNOAA_ShoreLine

    use ModuleGlobalData
    use ModuleEnterData

    implicit none

    private

    public  :: ConvertNOAA_ShoreLine
    private ::      ReadOptions
    private ::      Open_New_Output_File 

    !Types--------------------------------------------------------------
    type     T_NOAA_ShoreLine
        integer                                 :: InputUnit
        integer                                 :: OutputUnit           = 0
        character(len=StringLength)             :: FileName
        character(len=StringLength)             :: OutputFolder
        character(len=StringLength)             :: OutputName
        character(len=StringLength)             :: FirstPoint
        logical                                 :: WriteAsPolygons      = .true.
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjEnterDataNOAA     = 0
        integer                                 :: BufferSize           = null_int
        integer                                 :: MaxBufferSize
        logical                                 :: GettingFirstPoint    = .false.
    end type T_NOAA_ShoreLine

    type(T_NOAA_ShoreLine), pointer             :: Me
 
    !Parameters---------------------------------------------------------
    
    contains
    
    !------------------------------------------------------------------
    
    subroutine ConvertNOAA_ShoreLine(EnterDataID, STAT)
        
        !Arguments---------------------------------------------------------------
        integer,           intent(IN )              :: EnterDataID
        integer, optional, intent(OUT)              :: STAT

        !Local---------------------------------------------------------
        integer                                     :: STAT_CALL, nUsers
        integer                                     :: BlockSize, BufferLine
        character(len=25)                           :: FullBufferLine

        !Begin---------------------------------------------------------


        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData     = AssociateInstance(mEnterData_, EnterDataID)

        call ReadOptions

        write(*,*)"Reading file. Please wait..."

        call RemoveTabs

        call ConstructEnterData(Me%ObjEnterDataNOAA, Me%FileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ModuleNOAA_ShoreLine - ERR01'

        call GetBufferSize(Me%ObjEnterDataNOAA, BlockSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ModuleNOAA_ShoreLine - ERR02'

        write(*,*)"Converting file. Please wait..."

        Me%BufferSize = Me%MaxBufferSize + 1 

        do BufferLine = 1, BlockSize

            call GetFullBufferLine(Me%ObjEnterDataNOAA, BufferLine, FullBufferLine)

            if(scan(FullBufferLine, "#") .ne. 0)then

                if(Me%OutputUnit .ne. 0)then

                    !write(Me%OutputUnit,*) trim(Me%FirstPoint)

                    if(Me%WriteAsPolygons)then
                        write(Me%OutputUnit,*)"<endpolygon>"
                    else
                        write(Me%OutputUnit,*)"<end_line>"
                    endif 

                end if

                call Open_New_Output_File

                if(Me%WriteAsPolygons)then
                    write(Me%OutputUnit,*)"<beginpolygon>"
                else
                    write(Me%OutputUnit,*)"<begin_line>"
                endif 


                
                Me%GettingFirstPoint = .true.

            else
                if(Me%GettingFirstPoint)then
                    Me%FirstPoint        = FullBufferLine
                    Me%GettingFirstPoint = .false.
                end if               

                write(Me%OutputUnit,*)FullBufferLine

                Me%BufferSize = Me%BufferSize + 25

            end if

        end do
        
        !write(Me%OutputUnit,*) trim(Me%FirstPoint)

        if(Me%WriteAsPolygons)then
            write(Me%OutputUnit,*)"<endpolygon>"
        else
            write(Me%OutputUnit,*)"<end_line>"
        endif 


        call UnitsManager(Me%OutputUnit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ModuleNOAA_ShoreLine - ERR03'

        write(*,*)"Finished converting NOAA shoreline file. Please wait..."

        call KillEnterData(Me%ObjEnterDataNOAA, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ModuleNOAA_ShoreLine - ERR04'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'ConvertNOAA_ShoreLine - ModuleNOAA_ShoreLine - ERR05'
        
        STAT = SUCCESS_


    end subroutine ConvertNOAA_ShoreLine
    

    !--------------------------------------------------------------------------


    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                 :: iflag, STAT_CALL

        !Begin---------------------------------------------------------

        call GetData(Me%FileName,                                   &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'INPUT_FILENAME',               &
                     ClientModule = 'ConvertToXYZ',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNOAA_ShoreLine - ERR50'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleEtopo2 - ERR60'
        end if

        call GetData(Me%OutputFolder,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUT_FOLDER',                    &
                     ClientModule = 'ConvertToXYZ',                     &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNOAA_ShoreLine - ERR80'


        call GetData(Me%OutputName,                                     &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUT_NAME',                      &
                     Default      = 'Coast_Line',                       &
                     ClientModule = 'ConvertToXYZ',                     &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNOAA_ShoreLine - ERR80'



        call GetData(Me%WriteAsPolygons,                                &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'WRITE_AS_POLYGONS',                &
                     ClientModule = 'ConvertToXYZ',                     &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNOAA_ShoreLine - ERR80'

        call GetData(Me%MaxBufferSize,                                  &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'MAX_BUFFER_SIZE',                  &
                     Default      = 1000000,                            &
                     ClientModule = 'ConvertToXYZ',                     &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNOAA_ShoreLine - ERR80'





    end subroutine ReadOptions


    !------------------------------------------------------------------
    
    
    subroutine Open_New_Output_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        character(len=10)                           :: char_unit
        integer, save                               :: count
        character(len=StringLength)                 :: OutputFileName

        !Begin---------------------------------------------------------

        if(Me%BufferSize > Me%MaxBufferSize)then
            
            Me%BufferSize = 0

            if(Me%OutputUnit .ne. 0)then
                call UnitsManager(Me%OutputUnit, CLOSE_FILE, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'Open_New_Output_File - ModuleNOAA_ShoreLine - ERR01'
            end if

            count = count + 1

            write(char_unit, 10)count

            if(count>9)then
                char_unit =  '00'//trim(adjustl(char_unit))
            else
                char_unit = '000'//trim(adjustl(char_unit))
            end if

            write(*,*)"Converting file "//trim(char_unit)//'...'

            call UnitsManager(Me%OutputUnit, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_New_Output_File - ModuleNOAA_ShoreLine - ERR02'
            
            OutputFileName = trim(Me%OutputFolder)//trim(Me%OutputName)
            OutputFileName = trim(adjustl(OutputFileName))//trim(adjustl(char_unit))

            if(Me%WriteAsPolygons)then
                OutputFileName = trim(adjustl(OutputFileName))//'.xy'
            else
                OutputFileName = trim(adjustl(OutputFileName))//'.lin'
            endif 

            open(Unit   = Me%OutputUnit,        &
                 File   = OutputFileName,       &
                 STATUS = 'UNKNOWN',            &
                 IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Open_New_Output_File - ModuleNOAA_ShoreLine - ERR03'

            10 format(i5)

        end if

    end subroutine Open_New_Output_File

    !------------------------------------------------------------------


    subroutine RemoveTabs

        !Local-----------------------------------------------------------------
        character(len = StringLength)                           :: auxstring
        character(len = 1)                                      :: lchar
        integer                                                 :: i, BufferSize = 0, line
        character(len = StringLength), dimension(:), pointer    :: buffer
        logical                                                 :: FoundTabs = .false.
        integer                                                 :: Unit, STAT_CALL
        !Begin---------------------------------------------------------


        call UnitsManager(Unit, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RemoveTabs - ModuleNOAA_ShoreLine - ERR10'

        open(Unit, file = Me%FileName, action = 'read', status='old')

        write(*,*)
        write(*,*)"Getting buffer size..."

        do 
            read (Unit, "(A)", end=100) auxstring
            BufferSize = BufferSize + 1
        end do 

        100 continue

        allocate(buffer (BufferSize))

        rewind(Unit)

        write(*,*)
        write(*,*)"Searching tabs. Please wait..."

        do line = 1, BufferSize

            read(Unit, '(A)')auxstring

            do i = 2, StringLength

                lchar = auxstring(i-1:i-1)

                if(lchar == char(9))then
                    FoundTabs = .true.
                    auxstring(i-1:i-1) = char(32)
                end if

            end do

            buffer(line) = auxstring

        end do


        call UnitsManager(Unit, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RemoveTabs - ModuleNOAA_ShoreLine - ERR20'

        if(FoundTabs)then

            write(*,*)
            write(*,*)"Found tabs. Re-writing file..."

            call UnitsManager(Unit, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RemoveTabs - ModuleNOAA_ShoreLine - ERR30'



            open(Unit, file = Me%FileName, status = 'unknown')

            do line = 1, BufferSize
                write(Unit,*)trim(adjustl(buffer(line)))
            end do

            write(*,*)
            write(*,*)"Tab formatting successfully removed..."
            write(*,*)

            call UnitsManager(Unit, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RemoveTabs - ModuleNOAA_ShoreLine - ERR40'

        end if


    end subroutine RemoveTabs
    

end module ModuleNOAA_ShoreLine





