!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : NASA
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2004
! REVISION      : Pedro Galvao, Pedro Chambel Leitao, Luis Fernandes
! DESCRIPTION   : Module to convert NASA file into XYZ format within a defined window.
!                 Compile unformatted file conversion in BIG ENDIAN
!                 Source on the internet: ftp://edcsgs9.cr.usgs.gov/pub/data/srtm/
!
!------------------------------------------------------------------------------


!   INPUT_FOLDER                : char              -           !Path to input file to convert
!   XYZ_OUTPUT_FILENAME         : char              -           !Path to XYZ file generated


Module ModuleNASA

    use ModuleGlobalData
    use ModuleDrawing
    use ModuleEnterData
    use ModuleFunctions, only : LatLonToUTM

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertNASA
    private ::      ReadOptions
    private ::      OpenFile
    private ::      ConvertToXYZ
    private ::      KillNASA

    !Parameters----------------------------------------------------------------
    integer, parameter                  :: PointNoValue     = -32768

    !Types---------------------------------------------------------------------

    type    T_Map
        integer                         :: OriginX
        integer                         :: OriginY
        integer                         :: Unit
        character(len=StringLength)     :: FileName
        type(T_Map), pointer            :: Next
    end type T_Map

    type      T_NASA
        integer                         :: OutputUnit
        integer                         :: intminX
        integer                         :: intminY
        integer                         :: intmaxX
        integer                         :: intmaxY
        integer                         :: Sampling
        logical                         :: ConvertToUTM
        integer                         :: Datum
        type(T_Limits)                  :: Window
        character(len=StringLength)     :: Folder, OutputFileName
        type(T_Map),           pointer  :: FirstMap
        integer                         :: ObjEnterData = 0
    end type  T_NASA

    type(T_NASA), pointer               :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertNASA(Window, EnterDataID, STAT)

        !Arguments---------------------------------------------------------------
        type(T_Limits),    intent(IN )                  :: Window
        integer                                         :: EnterDataID
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_

        nullify (Me)
        allocate(Me)

        Me%Window       = Window
        Me%ObjEnterData = AssociateInstance(mEnterData_, EnterDataID)

        call ReadOptions

        call SelectFiles

        call ConvertToXYZ

        call KillNASA

        STAT = SUCCESS_

    end subroutine ConvertNASA

    !--------------------------------------------------------------------------


    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                 :: iflag, STAT_CALL
        
        !Begin---------------------------------------------------------


        call GetData(Me%Folder,                                     &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'INPUT_FOLDER',                 &
                     ClientModule = 'ModuleNASA',                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNASA - ERR01'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleNASA - ERR02'
        end if


        call GetData(Me%OutputFileName,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'XYZ_OUTPUT_FILENAME',          &
                     ClientModule = 'ModuleNASA',                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNASA - ERR50'

        !SAMPLE_POINTS=1 => Topography is read every point
        !SAMPLE_POINTS=2 =>Topography is read every one point from 2 points  
                !(1/(2^2) of the points=25%)
        !SAMPLE_POINTS=n =>Topography is read every one point from n points
                ! (1/(n^2) of the points)
        call GetData(Me%Sampling,                                   &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'SAMPLE_POINTS',                &
                     Default      =  1 ,                            &
                     ClientModule = 'ModuleNASA',                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNASA - ERR50'

        !Convert To UTM
        call GetData(Me%ConvertToUTM,                               &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'CONVERT_TO_UTM',               &
                     Default      =  .false.,                       &
                     ClientModule = 'ModuleNASA',                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNASA - ERR50'

        if (Me%ConvertToUTM) then

            call GetData(Me%Datum,                                      &
                         Me%ObjEnterData, iflag,                        &
                         SearchType   = FromBlock,                      &
                         keyword      = 'DATUM',                        &
                         Default      =  WGS_84_DATUM,                  &
                         ClientModule = 'ModuleNASA',                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNASA - ERR50'

        end if 


    end subroutine ReadOptions


    !--------------------------------------------------------------------------


    subroutine ConvertToXYZ

        !Local-----------------------------------------------------------------
        type(T_Map), pointer                        :: MapFile
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call UnitsManager(Me%OutputUnit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR02'

        open(Unit   = Me%OutputUnit,                &
             File   = Me%OutputFileName,            &
             STATUS = 'UNKNOWN',                    &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR03'



        write(Me%OutputUnit,*)"<begin_xyz>"

        MapFile => Me%FirstMap

        do while(associated(MapFile))

            call OpenFile(MapFile%FileName, MapFile%Unit)

            call ReadFileAndWriteOutput(MapFile)

            MapFile => MapFile%Next

        end do

        write(Me%OutputUnit,*)"<end_xyz>"

    end subroutine ConvertToXYZ

    !--------------------------------------------------------------------------

    subroutine SelectFiles

        !Local-----------------------------------------------------------------
        type(T_Map), pointer                        :: NewMap
        integer                                     :: i, j
        integer                                     :: nx, ny
        integer                                     :: latmin, latmax, longmin, longmax
        character(len=10)                           :: char_i, char_j, char_EW, char_NS
        character(len=StringLength)                 :: FileName
        logical                                     :: exist

        !Begin---------------------------------------------------------

        write(*,*)
        write(*,*)"Selecting NASA file to read..."

        nullify(Me%FirstMap)

        if (Me%Window%Bottom < 0.) then
            latmin = int(Me%Window%Bottom) - 1
        else
            latmin = int(Me%Window%Bottom)
        endif

        if (Me%Window%Top < 0.) then
            latmax = int(Me%Window%Top)
        else
            latmax = int(Me%Window%Top) + 1
        endif

        if (Me%Window%Left < 0.) then
            longmin = int(Me%Window%Left) - 1
        else
            longmin = int(Me%Window%Left)
        endif

        if (Me%Window%Right < 0.) then
            longmax = int(Me%Window%Right)
        else
            longmax = int(Me%Window%Right) + 1
        endif

        !Each cell has one degree so subtract one
        longmax = longmax - 1
        latmax  = latmax  - 1

        !maps allongx
        nx = longmax - longmin + 1
    
        !maps allongy
        ny = latmax - latmin + 1

        do i = longmin, longmax
        do j = latmin,  latmax 
            
            !East West component
            if (i < 0) then
                write(char_i, fmt = 10) -i 
                if (i > -10) then
                    char_EW = 'W00'
                elseif (i > -100) then
                    char_EW = 'W0'
                else
                    char_EW = 'W'
                end if    
            else
                write(char_i, fmt = 10)  i
                if (i < 10) then
                    char_EW = 'E00'
                elseif (i < 100) then
                    char_EW = 'E0'
                else 
                    char_EW = 'E'
                end if    
            endif

            !North South Component
            if (j < 0) then
                write(char_j, fmt = 10) -j
                if (j > -10) then
                    char_NS = 'S0'
                else
                    char_NS = 'S'
                end if    
            else
                write(char_j, fmt = 10)  j
                if (j < 10) then
                    char_NS = 'N0'
                else
                    char_NS = 'N'
                end if    
            endif

            Filename = trim(Me%Folder)//trim(adjustl(char_NS))// &
                       trim(adjustl(char_j))//trim(adjustl(char_EW))&
                       //trim(adjustl(char_i))//'.hgt'

            write(*,*)'Looking for :', trim(FileName)

            !Verifies if file exists
            inquire(file = FileName, exist = exist)
            if (.not. exist) then
                write(*,*) trim(FileName)
                write(*,*)'NASA file does not exist'
            else
                call AddMap(Me%FirstMap, NewMap)
                
                if (i < 0) then
                    NewMap%OriginX = i
                else
                    NewMap%OriginX = i
                endif
                if (j < 0) then
                    NewMap%OriginY = j
                else
                    NewMap%OriginY = j
                endif

                NewMap%Filename = Filename
            endif

        enddo
        enddo

        10 format(i10)

        write(*,*)
        write(*,*)"Selected ",nx*ny, "files..."


    end subroutine SelectFiles
    
    !------------------------------------------------------------------


    subroutine OpenFile(FileName, Unit)

        !Arguments-----------------------------------------------------
        character(len=*)                            :: FileName
        integer, intent(out)                        :: Unit

        !Local---------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin---------------------------------------------------------

        call UnitsManager(Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenFile - ModuleNASA - ERR02'

        open(Unit   = Unit,                 &
             File   = FileName,             &
             Form   = 'BINARY',             &
             STATUS = 'OLD',                &
             Action = 'READ',               &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenFile - ModuleNASA - ERR03'

        write(*,*)
        write(*,*)"Opened NASA file..."


    end subroutine OpenFile
    
    !------------------------------------------------------------------

    subroutine ReadFileAndWriteOutput(MapFile)

        !Arguments-----------------------------------------------------
        type(T_Map),      pointer           :: MapFile

        !Local---------------------------------------------------------
        integer                             :: i, j, STAT_CALL
        integer(2)                          :: val(1201, 1201), maxval
        real                                :: pointx, pointy 
        integer                             :: count_i, count_j
        real(8)                             :: utm_x, utm_y
        integer                             :: grid_zone(2)        

        !Begin---------------------------------------------------------

        write(*,*)
        write(*,*)"Reading NASA files. Please wait..."

        do i=1201,1,-1
        do j=1,1201

            read(MapFile%Unit) val(i,j)
            maxval = max(maxval, val(i,j))

        enddo
        enddo

        call UnitsManager(MapFile%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'KillNASA - ModuleNASA - ERR01'

        count_i = 1

        do i=1,1201
            count_j = 1
        do j=1,1201

            if (i == count_i .and. j == count_j) then

                PointX = MapFile%OriginX + (j-1) * 3. / 3600.
                PointY = MapFile%OriginY + (i-1) * 3. / 3600.

                if ((PointX .GE. Me%Window%Left  ) .and. (PointX .LE. Me%Window%Right)) then

                    if ((PointY .GE. Me%Window%Bottom) .and. (PointY .LE. Me%Window%Top)) then

                        if (val(i,j) .ne. PointNoValue) then

                            if (Me%ConvertToUTM) then
                                call LatLonToUTM (dble(pointx), dble(pointy), utm_x, utm_y, grid_zone, Me%Datum)                                
                                write(Me%OutputUnit, '(F20.6, 1X, F20.6, 1X, I10)') utm_x, utm_y, val(i,j)
                            else
                                write(Me%OutputUnit, '(F11.6, 1X, F10.6, 1X, I10)') pointx, pointy,  val(i,j)
                            endif


                        endif

                    end if

                end if

            endif

            if (j == count_j) then
                count_j = count_j + Me%Sampling
            endif


        enddo
            if (i == count_i ) then
                count_i = count_i + Me%Sampling
            endif
        enddo

        write(*,*)
        write(*,*)"Finished reading..."

    end subroutine ReadFileAndWriteOutput


    !------------------------------------------------------------------------
    
    subroutine AddMap (FirstMap, ObjMap)

        !Arguments-------------------------------------------------------------
        type (T_Map), pointer                   :: FirstMap
        type (T_Map), pointer                   :: ObjMap

        !Local-----------------------------------------------------------------
        type (T_Map), pointer                   :: NewMap
        type (T_Map), pointer                   :: PreviousMap
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewMap)
        nullify  (NewMap%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstMap)) then
            FirstMap         => NewMap
            ObjMap           => NewMap
        else
            PreviousMap      => FirstMap
            ObjMap           => FirstMap%Next
            do while (associated(ObjMap))
                PreviousMap  => ObjMap
                ObjMap       => ObjMap%Next
            enddo
            ObjMap           => NewMap
            PreviousMap%Next => NewMap
        endif


    end subroutine AddMap
    
    
    !------------------------------------------------------------------------
  
   
    subroutine KillNASA
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call UnitsManager(Me%OutputUnit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'KillNASA - ModuleNASA - ERR01'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0)           stop 'KillNASA - ModuleNASA - ERR02'

        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillNASA

    !--------------------------------------------------------------------------

 
end module ModuleNASA









