!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Etopo5
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2004
! REVISION      : Henrique Coelho, Paulo Chambel Leitao, Madalena Santos, Luis Fernandes
! DESCRIPTION   : Module to convert Etopo5 file into XYZ format within
!                 a defined window.
!
!------------------------------------------------------------------------------


!   INPUT_FILENAME              : char              -           !Path to input file to convert
!   OUTPUT_FILENAME             : char              -           !Path to XYZ file generated
!   WRITE_AS_BATHYMETRY         : 0/1               -           !Write XYZ values multiplied by -1.


Module ModuleEtopo5

    use ModuleGlobalData
    use ModuleDrawing
    use ModuleEnterData

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertEtopo5
    private ::      ReadOptions
    private ::      OpenFile
    private ::      PrepareReading
    private ::      KillEtopo5

    !Parameters----------------------------------------------------------------
    integer, parameter                  :: iLongitude   = 4320
    integer, parameter                  :: iLatitude    = 2160   

    !Types---------------------------------------------------------------------
    type      T_Etopo5
        integer                         :: Unit
        integer                         :: OutputUnit
        type(T_Limits)                  :: Window
        logical                         :: WriteAsBathymetry
        character(len=StringLength)     :: FileName, OutputFileName
        integer, dimension(:,:), pointer:: Topography
        real,    dimension(:  ), pointer:: Latitude
        real,    dimension(:  ), pointer:: Longitude
        integer                         :: ObjEnterData = 0
    end type  T_Etopo5

    type(T_Etopo5), pointer             :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertEtopo5(Window, EnterDataID, STAT)

        !Arguments---------------------------------------------------------------
        type(T_Limits),    intent(IN )                  :: Window
        integer,           intent(IN )                  :: EnterDataID
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_

        nullify (Me)
        allocate(Me)

        Me%Window               = Window
        Me%ObjEnterData         = AssociateInstance(mEnterData_, EnterDataID)

        call ReadOptions

        call OpenFile

        call PrepareReading

        call ReadFileAndWriteOutput

        call KillEtopo5

        STAT = SUCCESS_

    end subroutine ConvertEtopo5

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
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleEtopo5 - ERR50'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleEtopo5 - ERR60'
        end if


        call GetData(Me%OutputFileName,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'OUTPUT_FILENAME',              &
                     ClientModule = 'ModuleEtopo5',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleEtopo5 - ERR50'

        call GetData(Me%WriteAsBathymetry,                          &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'WRITE_AS_BATHYMETRY',          &
                     Default      = .false.,                        &
                     ClientModule = 'ModuleEtopo5',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleEtopo5 - ERR122'


    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    subroutine PrepareReading

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, j2

        !Begin---------------------------------------------------------


        write(*,*)
        write(*,*)"Preparing to read..."

        if(Me%Window%Left .lt. -180. .or. Me%Window%Right .gt. 180.)then
            write(*,*)'Window wrongly defined'
            write(*,*)'Window Left  : ', Me%Window%Left 
            write(*,*)'Window Right : ', Me%Window%Right 
            stop 'PrepareReading - ModuleEtopo5 - ERR01'
        end if

        if(Me%Window%Top .gt. 90.   .or. Me%Window%Bottom .lt. -90.)then
            write(*,*)'Window wrongly defined'
            write(*,*)'Window Top    : ', Me%Window%Top 
            write(*,*)'Window Bottom : ', Me%Window%Bottom 
            stop 'PrepareReading - ModuleEtopo5 - ERR02'
        end if

        allocate(Me%Topography(iLatitude, iLongitude))
        allocate(Me%Longitude (iLongitude))
        allocate(Me%Latitude  (iLatitude ))


        j2=iLongitude /2 

        do j=1,j2
            Me%Longitude(j)     = 0.    + float(j-1)/12
            Me%Longitude(j2+j)  = -180. + float(j-1)/12
        enddo

        do i=1,iLatitude
            Me%Latitude(i)      = 90    - float(i-1)/12
        enddo  


    end subroutine PrepareReading
    
    
    !--------------------------------------------------------------------------

    
    subroutine OpenFile

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: exist

        !Begin---------------------------------------------------------

        !Verifies if file exists
        inquire(file = Me%FileName, exist = exist)
        if (.not. exist) then
            write(*,*)'ETOPO file does not exist'
            stop 'OpenFile - ConvertToXYZ - ERR01'
        endif

        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenFile - ModuleEtopo5 - ERR02'

        open(Unit   = Me%Unit,              &
             File   = Me%FileName,          &
             STATUS = 'OLD',                &
             Action = 'READ',               &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenFile - ModuleEtopo5 - ERR03'

        rewind(Me%Unit)

        write(*,*)
        write(*,*)"Opened Etopo5 file..."


    end subroutine OpenFile
    
    !------------------------------------------------------------------


    subroutine ReadFileAndWriteOutput


        !Local-----------------------------------------------------------------
        integer                                     :: ilat, ilon, i, j, STAT_CALL

        !Begin---------------------------------------------------------


        call UnitsManager(Me%OutputUnit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileAndWriteOutput - ModuleEtopo5 - ERR02'

        open(Unit   = Me%OutputUnit,                &
             File   = Me%OutputFileName,            &
             STATUS = 'UNKNOWN',                    &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileAndWriteOutput - ModuleEtopo5 - ERR03'



        write(*,*)
        write(*,*)"Reading Etopo5 file. Please wait..."

        do ilat=1,ilatitude
        do ilon=1,216
            read(Me%Unit,20) (Me%Topography(ilat,j),j=(ilon-1)*20+1,ilon*20)
        enddo
        enddo

        write(*,*)
        write(*,*)"Writing XYZ file. Please wait..."


        write(Me%OutputUnit,*)"<begin_xyz>"

        do i = 1, iLatitude             
        do j = 1, iLongitude

            if(Me%Latitude  (i  ) .ge. Me%Window%Bottom         .and. &
               Me%Latitude  (i  ) .le. Me%Window%Top            .and. &
               Me%Longitude (j  ) .ge. Me%Window%Left           .and. &
               Me%Longitude (j  ) .le. Me%Window%Right          .and. &
               Me%Topography(i,j) .ge. Me%Window%MinimumValue   .and. &
               Me%Topography(i,j) .le. Me%Window%MaximumValue)then

                if(Me%WriteAsBathymetry)then
                    write(Me%OutputUnit,10) Me%Longitude(j),Me%Latitude(i),-Me%Topography(i,j)
                else
                    write(Me%OutputUnit,10) Me%Longitude(j),Me%Latitude(i),Me%Topography(i,j)
                end if

            end if

        enddo
        enddo


        write(Me%OutputUnit,*)"<end_xyz>"
        
        write(*,*)
        write(*,*)"Finished reading..."

        10 format(F10.5, 1X, F9.5, 1X, I5)
        20 format(20I6)

    end subroutine ReadFileAndWriteOutput
    
    
    subroutine KillEtopo5
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        
        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'KillEtopo5 - ModuleEtopo5 - ERR01'

        call UnitsManager(Me%OutputUnit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'KillEtopo5 - ModuleEtopo5 - ERR02'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0)           stop 'KillEtopo2 - ModuleEtopo5 - ERR03'

        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillEtopo5

    !--------------------------------------------------------------------------

 
end module ModuleEtopo5









