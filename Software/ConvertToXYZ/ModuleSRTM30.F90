!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : SRTM30
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2004
! REVISION      : Luis Fernandes
! DESCRIPTION   : Module to convert SRTM30 file into XYZ format within a defined window.
!                 Compile unformatted file conversion in BIG ENDIAN
!                 Source on the internet: ftp://topex.ucsd.edu/pub/srtm30_plus/data/
!
!------------------------------------------------------------------------------


!   INPUT_FILE                  : char              -           !Path to input file to convert
!   XYZ_OUTPUT_FILENAME         : char              -           !Path to XYZ file generated
!   WRITE_AS_BATHYMETRY         : 0/1               -           !Write XYZ values multiplied by -1.
!   ORIGIN_X                    : real              -           !Lower left corner XX coordinate
!   ORIGIN_Y                    : real              -           !Lower left corner YY coordinate

Module ModuleSRTM30

    use ModuleGlobalData
    use ModuleDrawing
    use ModuleEnterData
    use ModuleFunctions

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertSRTM30
    private ::      ReadOptions
    private ::      OpenFile
    private ::      ConvertToXYZ
    private ::      KillSRTM30

    !Parameters----------------------------------------------------------------
    integer, parameter                  :: PointNoValue     = -32768

    !Types---------------------------------------------------------------------

    type    T_Map
        type(T_Map), pointer            :: Next
    end type T_Map

    type      T_SRTM30
        integer                         :: OriginX
        integer                         :: OriginY
        integer                         :: InputUnit
        integer                         :: OutputUnit
        integer                         :: intminX
        integer                         :: intminY
        integer                         :: intmaxX
        integer                         :: intmaxY
        integer                         :: Datum
        type(T_Limits)                  :: Window
        character(len=PathLength)       :: InputFile
        character(len=PathLength)       :: OutputFile
        integer                         :: ObjEnterData = 0
    end type  T_SRTM30

    type(T_SRTM30), pointer             :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertSRTM30(Window, EnterDataID, STAT)

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

        call OpenFile

        call ConvertToXYZ

        call KillSRTM30

        STAT = SUCCESS_

    end subroutine ConvertSRTM30

    !--------------------------------------------------------------------------


    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                 :: iflag, STAT_CALL
        
        !Begin---------------------------------------------------------


        call GetData(Me%InputFile,                                  &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'INPUT_FILE',                   &
                     ClientModule = 'ModuleSRTM30',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleSRTM30 - ERR01'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleSRTM30 - ERR02'
        end if


        call GetData(Me%OutputFile,                                 &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'XYZ_OUTPUT_FILENAME',          &
                     ClientModule = 'ModuleSRTM30',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleSRTM30 - ERR10'


        call GetData(Me%OriginX,                                    &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'ORIGIN_X',                     &
                     ClientModule = 'ModuleSRTM30',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleSRTM30 - ERR20'

        call GetData(Me%OriginY,                                    &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'ORIGIN_Y',                     &
                     ClientModule = 'ModuleSRTM30',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleSRTM30 - ERR30'

    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    subroutine ConvertToXYZ

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call UnitsManager(Me%OutputUnit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR02'

        open(Unit   = Me%OutputUnit,                &
             File   = Me%OutputFile,                &
             STATUS = 'UNKNOWN',                    &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR03'

        write(Me%OutputUnit,*)"<begin_xyz>"

        call ReadFileAndWriteOutput

        write(Me%OutputUnit,*)"<end_xyz>"

    end subroutine ConvertToXYZ

    !--------------------------------------------------------------------------

    subroutine OpenFile

        !Local---------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin---------------------------------------------------------

        call UnitsManager(Me%InputUnit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenFile - ModuleSRTM30 - ERR02'

        open(Unit   = Me%InputUnit,         &
             File   = Me%InputFile,         &
             Form   = 'BINARY',             &
             STATUS = 'OLD',                &
             Action = 'READ',               &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenFile - ModuleSRTM30 - ERR03'

        write(*,*)
        write(*,*)"Opened SRTM30 file..."


    end subroutine OpenFile
    
    !------------------------------------------------------------------

    subroutine ReadFileAndWriteOutput

        !Local---------------------------------------------------------
        integer                             :: i, j, STAT_CALL
        integer(2)                          :: val(6000, 4800)
        real                                :: pointx, pointy 

        !Begin---------------------------------------------------------

        write(*,*)
        write(*,*)"Reading SRTM30 files. Please wait..."
        
        do i=6000,1,-1
        do j=1,4800
            read(Me%InputUnit)val(i,j)
        enddo
        enddo

        call UnitsManager(Me%InputUnit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'KillSRTM30 - ModuleSRTM30 - ERR01'

        do i=1,6000
        do j=1,4800

            PointX = Me%OriginX + (j-1) * 30. / 3600.
            PointY = Me%OriginY + (i-1) * 30. / 3600.

            if ((PointX .GE. Me%Window%Left) .and. (PointX .LE. Me%Window%Right)) then
                
                if ((PointY .GE. Me%Window%Bottom) .and. (PointY .LE. Me%Window%Top)) then
                    
                    if (val(i,j) .ne. PointNoValue) then
                        
                        if(val(i,j) .ge. Me%Window%MinimumValue .and. &
                           val(i,j) .le. Me%Window%MaximumValue)then

                            write(Me%OutputUnit, 999) PointX, PointY, val(i,j) * (-1)

                        end if

                    endif

                end if

            end if

        enddo
        enddo

        999 format(F11.6, 1X, F10.6, 1X, I10)

        write(*,*)
        write(*,*)"Finished reading..."

    end subroutine ReadFileAndWriteOutput


    !------------------------------------------------------------------------
  
   
    subroutine KillSRTM30
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call UnitsManager(Me%OutputUnit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'KillSRTM30 - ModuleSRTM30 - ERR01'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0)           stop 'KillSRTM30 - ModuleSRTM30 - ERR02'

        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillSRTM30

    !--------------------------------------------------------------------------

 
end module ModuleSRTM30









