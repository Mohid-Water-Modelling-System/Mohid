!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Etopo2
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2004
! REVISION      : Henrique Coelho, Paulo Chambel Leitao, Madalena Santos, Luis Fernandes
! DESCRIPTION   : Module to convert Etopo2 file into XYZ format within a defined window.
!                 Source on the Internet: http://dss.ucar.edu/datasets/ds759.3/
!
!------------------------------------------------------------------------------



!   INPUT_FILENAME              : char              -           !Path to input file to convert
!   OUTPUT_FILENAME             : char              -           !Path to XYZ file generated
!   WRITE_AS_BATHYMETRY         : 0/1               -           !Write XYZ values multiplied by -1.



Module ModuleEtopo2

    use ModuleGlobalData
    use ModuleDrawing
    use ModuleEnterData

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertEtopo2
    private ::      ReadOptions
    private ::      OpenFile
    private ::      PrepareReading
    private ::      KillEtopo2

    !Parameters----------------------------------------------------------------
    integer, parameter                  :: iLongitude   = 10800
    integer, parameter                  :: iLatitude    = 5400   

    !Types---------------------------------------------------------------------
    type      T_Etopo2
        integer                         :: Unit
        integer                         :: OutputUnit
        type(T_Limits)                  :: Window
        logical                         :: WriteAsBathymetry
        character(len=StringLength)     :: FileName
        character(len=StringLength)     :: OutputFileName
        integer, dimension(:), pointer  :: Topography
        real,    dimension(:), pointer  :: Latitude
        real,    dimension(:), pointer  :: Longitude
        integer                         :: ObjEnterData = 0
    end type  T_Etopo2

    type(T_Etopo2), pointer             :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertEtopo2(Window, EnterDataID, STAT)

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

        call KillEtopo2

        STAT = SUCCESS_

    end subroutine ConvertEtopo2


    !--------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                 :: iflag, STAT_CALL

        !Begin---------------------------------------------------------

        call GetData(Me%FileName,                                   &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'INPUT_FILENAME',               &
                     ClientModule = 'ModuleEtopo2',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleEtopo2 - ERR50'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleEtopo2 - ERR60'
        end if

        call GetData(Me%OutputFileName,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'OUTPUT_FILENAME',              &
                     ClientModule = 'ModuleEtopo2',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleEtopo2 - ERR50'


        call GetData(Me%WriteAsBathymetry,                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'WRITE_AS_BATHYMETRY',              &
                     Default      = .false.,                            &
                     ClientModule = 'ModuleEtopo2',                     &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleEtopo2 - ERR122'


    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    
    subroutine PrepareReading

        !Local-----------------------------------------------------------------
        integer                                     :: i, j

        !Begin---------------------------------------------------------


        write(*,*)
        write(*,*)"Preparing to read..."

        if(Me%Window%Left .lt. -180. .or. Me%Window%Right .gt. 180.)then
            write(*,*)'Window wrongly defined'
            write(*,*)'Window Left  : ', Me%Window%Left 
            write(*,*)'Window Right : ', Me%Window%Right 
            stop 'PrepareReading - ModuleEtopo2 - ERR01'
        end if

        if(Me%Window%Top .gt. 90.   .or. Me%Window%Bottom .lt. -90.)then
            write(*,*)'Window wrongly defined'
            write(*,*)'Window Top    : ', Me%Window%Top 
            write(*,*)'Window Bottom : ', Me%Window%Bottom 
            stop 'PrepareReading - ModuleEtopo2 - ERR02'
        end if

        allocate(Me%Topography(iLongitude))
        allocate(Me%Longitude (iLongitude))
        allocate(Me%Latitude  (iLatitude ))

        do i = 1, iLongitude
            Me%Longitude(i) = -180.0 + (float(i)-1.0)*2.0/60.0 
        enddo

        do j = 1, iLatitude
            Me%Latitude(j)  =  90.0  - (float(j)-1.0)*2.0/60.0
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
            stop 'OpenFile - ModuleEtopo2 - ERR01'
        endif

        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenFile - ModuleEtopo2 - ERR02'

        open(Unit   = Me%Unit,              &
             File   = Me%FileName,          &
             STATUS = 'OLD',                &
             Action = 'READ',               &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenFile - ModuleEtopo2 - ERR03'

        rewind(Me%Unit)

        write(*,*)
        write(*,*)"Opened Etopo2 file..."


    end subroutine OpenFile
    
    !------------------------------------------------------------------


    subroutine ReadFileAndWriteOutput


        !Local-----------------------------------------------------------------
        integer                                     :: ilonl, ilonr
        integer                                     :: nlon,  nlat
        integer                                     :: ilatb, ilatt
        integer                                     :: ndlon, ndlat
        real                                        :: rlat,  rlon
        integer                                     :: i, j, jnext
        integer                                     :: STAT_CALL

        !Begin---------------------------------------------------------


        call UnitsManager(Me%OutputUnit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileAndWriteOutput - ModuleEtopo2 - ERR02'

        open(Unit   = Me%OutputUnit,                &
             File   = Me%OutputFileName,            &
             STATUS = 'UNKNOWN',                    &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileAndWriteOutput - ModuleEtopo2 - ERR03'



        ilonl = nint(Me%Window%Left  + 180.)*30 + 1
        ilonr = nint(Me%Window%Right + 180.)*30 + 1

        nlon = ilonr - ilonl + 1

        ilatb = nint(90.- Me%Window%Bottom)*30 + 1
        ilatt = nint(90.- Me%Window%Top   )*30 + 1

        if(ilonr .eq. iLongitude + 1)ilonr = iLongitude  !boundary check
        if(ilatb .eq. iLatitude  + 1)ilatb = iLatitude   !boundary check

        nlon = ilonr - ilonl + 1
        nlat = ilatb - ilatt + 1

        ndlon = nlon + 1
        ndlat = nlat + 1

        if(ndlon.gt.iLongitude)ndlon = ndlon - 1
        if(ndlat.gt.iLatitude )ndlat = ndlat - 1

        write(*,*)
        write(*,*)"Reading Etopo2 file. Please wait..."

        write(Me%OutputUnit,*)"<begin_xyz>"

        jnext = 0
        rlat  = Me%Window%Top
        rlon  = Me%Window%Left

        do j = 1, iLatitude
        
            if(j .gt. ilatb) exit
        
            read(Me%Unit,'(15i6)')(Me%Topography(i),i=1, iLongitude)
        
            if(j .lt. ilatt .or. j .lt. jnext) cycle
        
            do i = ilonl, ilonr

                if(Me%Topography(i) .ge. Me%Window%MinimumValue .and. &
                   Me%Topography(i) .le. Me%Window%MaximumValue)then
                    
                    if(Me%WriteAsBathymetry)then
                        write(Me%OutputUnit,10)rlon,rlat,-Me%Topography(i)
                    else
                        write(Me%OutputUnit,10)rlon,rlat,Me%Topography(i)
                    endif

                end if

                rlon=rlon+2./60.

            enddo

            rlat    = rlat-2./60.
            rlon    = Me%Window%Left
            jnext   = j + 1
	  	     
        enddo

        write(Me%OutputUnit,*)"<end_xyz>"
        
        write(*,*)
        write(*,*)"Finished reading..."

        10 format(F10.5, 1X, F9.5, 1X, I5)

    end subroutine ReadFileAndWriteOutput
    
    
    subroutine KillEtopo2
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        
        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'KillEtopo2 - ModuleEtopo2 - ERR01'

        call UnitsManager(Me%OutputUnit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'KillEtopo2 - ModuleEtopo2 - ERR02'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0)           stop 'KillEtopo2 - ModuleEtopo2 - ERR03'

        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillEtopo2

    !--------------------------------------------------------------------------

 
end module ModuleEtopo2









