!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid ConvertToHDF5
! MODULE        : HYCOMFormat
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : July 2003
! REVISION      : Luis Fernandes
! DESCRIPTION   : Module to convert HYCOM format files into HDF5 format.
!
!------------------------------------------------------------------------------


Module ModuleHYCOMFormat

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalGrid
    use ModuleDrawing
    use netcdf

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertHYCOMFormat
    private ::      ReadOptions
    private ::      OpenAndReadHYCOMFiles
    private ::          OpenHDF5OutputFile
    private ::      KillHYCOMFormat
    

    !Parameters----------------------------------------------------------------



    !Types---------------------------------------------------------------------
    
    private :: T_HYCOMFormat
    type       T_HYCOMFormat
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
        character(len=PathLength)               :: HYCOMGridFileName
        character(len=PathLength)               :: HYCOMTempFileName
        character(len=PathLength)               :: HYCOMSaltFileName
        character(len=PathLength)               :: HYCOMUvelFileName
        character(len=PathLength)               :: HYCOMVvelFileName
        character(len=PathLength)               :: HYCOMUbaroFileName
        character(len=PathLength)               :: HYCOMVbaroFileName
        character(len=PathLength)               :: HYCOMSSHFileName



        logical                                 :: TempON   = .false.
        logical                                 :: SaltON   = .false.
        logical                                 :: UvelON   = .false.
        logical                                 :: VvelON   = .false.
        logical                                 :: UbaroON  = .false.
        logical                                 :: VbaroON  = .false.
        logical                                 :: SSHON    = .false.

        character(len=PathLength)               :: OutputGridFileName
        character(len=PathLength)               :: OutputFileName

        integer                                 :: HYCOMGridID
        integer                                 :: HYCOMTempID
        integer                                 :: HYCOMSaltID
        integer                                 :: HYCOMUvelID
        integer                                 :: HYCOMVvelID
        integer                                 :: HYCOMUbaroID
        integer                                 :: HYCOMVbaroID
        integer                                 :: HYCOMSSHID


        integer                                 :: i_top, i_bottom
        integer                                 :: j_left, j_right
        integer,dimension(:,:,:),   pointer     :: WaterPoints3D
        integer,dimension(:,:  ),   pointer     :: WaterPoints2D
        real,   dimension(:,:  ),   pointer     :: Bathymetry
        real,   dimension(:,:  ),   pointer     :: ConnectionX
        real,   dimension(:,:  ),   pointer     :: ConnectionY
        real                                    :: Xorig
        real                                    :: Yorig
        logical                                 :: InputGridInMohidFormat = .false.
        type(T_Limits)                          :: Window
        type(T_Size3D)                          :: Size, WorkSize, RealSize
        type(T_Time)                            :: Time
    end type  T_HYCOMFormat

    type(T_HYCOMFormat), pointer                :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertHYCOMFormat(EnterDataID, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        call ReadOptions

        call OpenHDF5OutputFile

        call OpenAndReadHYCOMFiles

        call KillHYCOMFormat

        STAT = SUCCESS_


    end subroutine ConvertHYCOMFormat

    !------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------
       
        call GetData(Me%HYCOMGridFileName,                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'HYCOM_GRID_FILENAME',              &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR10'


        call GetData(Me%HYCOMTempFileName,                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'HYCOM_TEMP_FILENAME',              &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR20'

        if (iflag == 1) Me%TempON = .true.

        call GetData(Me%HYCOMSaltFileName,                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'HYCOM_SALT_FILENAME',              &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR30'

        if (iflag == 1) Me%SaltON = .true.


        call GetData(Me%HYCOMUvelFileName,                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'HYCOM_UVEL_FILENAME',              &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR40'

        if (iflag == 1) Me%UvelON = .true.

        call GetData(Me%HYCOMVvelFileName,                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'HYCOM_VVEL_FILENAME',              &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR50'

        if (iflag == 1) Me%VvelON = .true.

        call GetData(Me%HYCOMUbaroFileName,                             &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'HYCOM_UBARO_FILENAME',             &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR60'

        if (iflag == 1) Me%UbaroON = .true.


        call GetData(Me%HYCOMVbaroFileName,                             &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'HYCOM_VBARO_FILENAME',             &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR70'

        if (iflag == 1) Me%VbaroON = .true.

        call GetData(Me%HYCOMSSHFileName,                               &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'HYCOM_SSH_FILENAME',               &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR80'

        if (iflag == 1) Me%SSHON = .true.


        call GetData(Me%OutputGridFileName,                             &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',             &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR90'



        call GetData(Me%OutputFileName,                                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUTFILENAME',                   &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR100'


        call GetData(Me%Window%Left,                                    &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'LEFT',                             &
                     Default      = -99999.,                            &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR110'

        call GetData(Me%Window%Right,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'RIGHT',                            &
                     Default      =  99999.,                            &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR120'

        call GetData(Me%Window%Top,                                     &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'TOP',                              &
                     Default      =  99999.,                            &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR130'

        call GetData(Me%Window%Bottom,                                  &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'BOTTOM',                           &
                     Default      =  -99999.,                           &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHYCOMFormat - ERR140'


    end subroutine ReadOptions


    !--------------------------------------------------------------------------

    subroutine OpenAndReadHYCOMFiles


        call ConstructTime

        call ConstructGrid


        if (Me%SaltON)  call ConstructSalt
        if (Me%TempON)  call ConstructTemp

        if (Me%UvelON)  call ConstructUvel
        if (Me%VvelON)  call ConstructVvel

        if (Me%UbaroON) call ConstructUbaro
        if (Me%VbaroON) call ConstructVbaro

        if (Me%SSHON)   call ConstructSSH


    end subroutine OpenAndReadHYCOMFiles
    
    
    !------------------------------------------------------------------------

    
    subroutine OpenHDF5OutputFile

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, trim(Me%OutputFileName), HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OpenHDF5OutputFile - ModuleHYCOMFormat - ERR02'

    end subroutine OpenHDF5OutputFile



    subroutine ConstructTime
        
        !Local-----------------------------------------------------------------
        integer                                 :: stat, STAT_CALL
        integer                                 :: DateID, DateLength, DateType
        character(len=80)                       :: DateStr
        logical                                 :: exist
        real                                    :: Year, Month, Day
        real                                    :: Hour, Minute, Second
        real,    dimension(6), target           :: AuxTime
        real,    dimension(:), pointer          :: TimePtr

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing time...'

        !Verifies if file exists
        inquire(file = Me%HYCOMGridFileName, exist = exist)
        if (.not. exist) then
            write(*,*)'HYCOM grid NETCDF file does not exist'
            stop 'ConstructTime - ModuleHYCOMFormat - ERR01'
        endif

        stat = nf90_open(Me%HYCOMGridFileName, NF90_NOWRITE, Me%HYCOMGridID)
        if (stat /= nf90_noerr) stop 'ConstructTime - ModuleHYCOMFormat - ERR02'

        stat = nf90_inquire_attribute(Me%HYCOMGridID, NF90_GLOBAL, 'field_date', DateType,      &
                                      DateLength, DateID)
        if (stat /= nf90_noerr) stop 'ConstructTime - ModuleHYCOMFormat - ERR03'

        stat =  nf90_get_att(Me%HYCOMGridID, NF90_GLOBAL, 'field_date', DateStr(1:DateLength))
        if (stat /= nf90_noerr) stop 'ConstructTime - ModuleHYCOMFormat - ERR04'

        read(DateStr(1 :4 ),*)Year
        read(DateStr(5 :6 ),*)Month
        read(DateStr(7 :8 ),*)Day
        read(DateStr(10:11),*)Hour
        read(DateStr(13:14),*)Minute
        read(DateStr(16:17),*)Second

        call SetDate(Me%Time, Year, Month, Day, Hour, Minute, Second)


        call ExtractDate(Me%Time, Year, Month, Day, Hour, Minute, Second) 

        AuxTime(1) = Year
        AuxTime(2) = Month
        AuxTime(3) = Day
        AuxTime(4) = Hour
        AuxTime(5) = Minute
        AuxTime(6) = Second

        TimePtr => AuxTime

        
        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructTime - ModuleHYCOMFormat - ERR05'

        
        call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                             "Time", "YYYY/MM/DD HH:MM:SS",             &
                             Array1D = TimePtr,                         &
                             OutputNumber = 1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructTime - ModuleHYCOMFormat - ERR06'


    end subroutine ConstructTime

    !------------------------------------------------------------------------
    
    subroutine ConstructGrid
        
        !Local-----------------------------------------------------------------
        integer                                 :: stat, nDimensions, nAttributes
        integer                                 :: lon_id, lat_id, lthk_id, var_type
        integer, dimension(4)                   :: DimsID
        type(T_Size2D)                          :: WorkSize2D
        real, dimension(:), pointer             :: Longitude, XX
        real, dimension(:), pointer             :: Latitude,  YY
        character(len=80)                       :: LonStr, LatStr, lthkStr
        integer                                 :: i, j, k, STAT_CALL


        real, dimension(:,:,:), pointer         :: Layers, WorkLayers, SZZ
        real, dimension(:,:  ), pointer         :: Layers2D
        real, dimension(:,:  ), pointer         :: LayersT

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing grid...'

        stat = nf90_inq_varid(Me%HYCOMGridID, 'Longitude', lon_id)
        if (stat /= nf90_noerr) stop 'ConstructGrid - ModuleHYCOMFormat - ERR04'


        stat = nf90_inquire_variable(Me%HYCOMGridID, lon_id, LonStr,        &
                                     var_type, nDimensions,                 &
                                     DimsID, nAttributes)
        if (stat /= nf90_noerr) stop 'OpenReadAndConvertFile - ModuleHYCOMFormat - ERR06'


        Me%RealSize%ILB = 1
        Me%RealSize%JLB = 1
        Me%RealSize%IUB = 1
        Me%RealSize%JUB = 1
        Me%RealSize%KLB = 1
        Me%RealSize%KUB = 1

        do i = 1, nDimensions

            stat = nf90_inquire_dimension(Me%HYCOMGridID, DimsID(i), LonStr, Me%RealSize%JUB)

        enddo

        allocate(Longitude(1:Me%RealSize%JUB))
        
        stat = nf90_get_var(Me%HYCOMGridID, lon_id, Longitude)
        if (stat /= nf90_noerr) stop 'GetVariable - ModuleHYCOMFormat - ERR01'


        stat = nf90_inq_varid(Me%HYCOMGridID, 'Latitude', lat_id)
        if (stat /= nf90_noerr) stop 'ConstructGrid - ModuleHYCOMFormat - ERR04'


        stat = nf90_inquire_variable(Me%HYCOMGridID, lat_id, LatStr,        &
                                     var_type, nDimensions,                 &
                                     DimsID, nAttributes)
        if (stat /= nf90_noerr) stop 'OpenReadAndConvertFile - ModuleHYCOMFormat - ERR06'

        do i = 1, nDimensions

            stat = nf90_inquire_dimension(Me%HYCOMGridID, DimsID(i), LatStr, Me%RealSize%IUB)

        enddo

        allocate(Latitude(1:Me%RealSize%IUB))
        
        stat = nf90_get_var(Me%HYCOMGridID, lat_id, Latitude)
        if (stat /= nf90_noerr) stop 'GetVariable - ModuleHYCOMFormat - ERR01'


        do j = 1, Me%RealSize%JUB
            if(Longitude(j) .ge. Me%Window%Left)then
                Me%j_left = j
                exit
            endif
        enddo   


        do j = Me%RealSize%JUB, Me%j_left, -1
            if(Longitude(j) .le. Me%Window%Right)then
                Me%j_right = j
                exit
            endif
        enddo   

        do i = 1, Me%RealSize%IUB
            if(Latitude(i) .ge. Me%Window%Bottom)then
                Me%i_bottom = i
                exit
            endif
        enddo   


        do i = Me%RealSize%IUB, Me%i_bottom, -1
            if(Latitude(i) .le. Me%Window%Top)then
                Me%i_top = i
                exit
            endif
        enddo


        Me%Xorig = Longitude(Me%j_left  ) + (Longitude(Me%j_left+1  ) - Longitude(Me%j_left  ))/2.
        Me%Yorig = Latitude (Me%i_bottom) + (Latitude (Me%i_bottom+1) - Latitude (Me%i_bottom))/2.


        Me%WorkSize%ILB = 1
        Me%WorkSize%JUB = Me%j_right - Me%j_left
        Me%WorkSize%JLB = 1
        Me%WorkSize%IUB = Me%i_top - Me%i_bottom

        allocate(XX(1:Me%WorkSize%JUB+1))
        allocate(YY(1:Me%WorkSize%IUB+1))

        XX(:) = 0.
        YY(:) = 0.

        do j = Me%WorkSize%JLB+1, Me%WorkSize%JUB

            XX(j) = XX(j-1) + (Longitude(Me%j_left + j) - Longitude(Me%j_left + j - 1))

        enddo

        XX(Me%WorkSize%JUB+1) = XX(Me%WorkSize%JUB) +                      &
                                (Longitude(Me%j_left + Me%WorkSize%JUB) - &
                                 Longitude(Me%j_left + Me%WorkSize%JUB-1))
        
        do i = Me%WorkSize%ILB+1, Me%WorkSize%IUB

            YY(i) = YY(i-1) + (Latitude(Me%i_bottom + i) - Latitude(Me%i_bottom + i-1))

        enddo
        
        YY(Me%WorkSize%IUB+1) = YY(Me%WorkSize%IUB) + (Latitude(Me%i_bottom + Me%WorkSize%IUB) &
                                                     - Latitude(Me%i_bottom + Me%WorkSize%IUB-1))






        stat = nf90_inq_varid(Me%HYCOMGridID, 'lthk', lthk_id)
        if (stat /= nf90_noerr) stop 'ConstructGrid - ModuleHYCOMFormat - ERR04'



        stat = nf90_inquire_variable(Me%HYCOMGridID, lthk_id, lthkStr,          &
                                     var_type, nDimensions,                     &
                                     DimsID, nAttributes)
        if (stat /= nf90_noerr) stop 'OpenReadAndConvertFile - ModuleHYCOMFormat - ERR06'

        stat = nf90_inquire_dimension(Me%HYCOMGridID, DimsID(3), lthkStr, Me%RealSize%KUB)
        if (stat /= nf90_noerr) stop 'OpenReadAndConvertFile - ModuleHYCOMFormat - ERR06'


        Me%WorkSize%KLB = 1
        Me%WorkSize%KUB = Me%RealSize%KUB
        

        allocate(Layers(1:Me%RealSize%JUB, 1:Me%RealSize%IUB, 1:Me%RealSize%KUB))

        stat = nf90_get_var(Me%HYCOMGridID, lthk_id, Layers)
        if (stat /= nf90_noerr) stop 'GetVariable - ModuleHYCOMFormat - ERR01'

        allocate(LayersT(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))

        allocate(WorkLayers(Me%WorkSize%ILB:Me%WorkSize%IUB, &
                            Me%WorkSize%JLB:Me%WorkSize%JUB, &
                            Me%WorkSize%KLB:Me%WorkSize%KUB))
        
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

            Layers2D => Layers(Me%j_left:Me%j_right, Me%i_bottom:Me%i_top, k)

            LayersT = transpose(Layers2D)

            WorkLayers(:,:,Me%WorkSize%KUB - k + 1) = LayersT

        enddo

        
        allocate(Me%Bathymetry(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%ILB:Me%WorkSize%JUB))
        Me%Bathymetry = 0.


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            
            Me%Bathymetry(i,j) = Me%Bathymetry(i,j) + WorkLayers(i, j, k)

        enddo
        enddo
        enddo

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if(Me%Bathymetry(i,j) > 10000.)then
                Me%Bathymetry(i,j) = -99
            end if

        enddo
        enddo


        WorkSize2D%ILB = Me%WorkSize%ILB
        WorkSize2D%IUB = Me%WorkSize%IUB
        WorkSize2D%JLB = Me%WorkSize%JLB
        WorkSize2D%JUB = Me%WorkSize%JUB

        call WriteGridData (FileName        = Me%OutputGridFileName,        &
                            XX              = XX,                           &
                            YY              = YY,                           &
                            COMENT1         = " HYCOM Grid based on file :",&
                            COMENT2         = Me%HYCOMGridFileName,         &
                            WorkSize        = WorkSize2D,                   & 
                            CoordType       = 1,                            &
                            Xorig           = Me%Xorig,                     &
                            Yorig           = Me%Yorig,                     &
                            Zone            = 0,                            &
                            Grid_Angle      = 0.,                           &
                            Latitude        = 0.,                           &
                            Longitude       = 0.,                           &
                            GridData2D_Real = Me%Bathymetry,                &
                            Overwrite       = ON,                           &
                            FillValue       = -99.,                         &
                            STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleHYCOMFormat - ERR01'

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%OutputGridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleHYCOMFormat - ERR02'


        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleHYCOMFormat - ERR02'


        allocate(SZZ(Me%WorkSize%ILB  :Me%WorkSize%IUB,   &
                     Me%WorkSize%JLB  :Me%WorkSize%JUB,   &
                     Me%WorkSize%KLB-1:Me%WorkSize%KUB))


        SZZ(:,:, Me%WorkSize%KLB-1) = Me%Bathymetry(:,:)

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            
            SZZ(i,j, k) = SZZ(i,j, k-1) - WorkLayers(i, j, k)

        enddo
        enddo
        enddo


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do k = Me%WorkSize%KLB-1, Me%WorkSize%KUB
            
            if(SZZ(i,j,k) < -10000. .or. SZZ(i,j,k) > 10000.)then
                
                SZZ(i,j,k) = null_real

            endif

        enddo
        enddo
        enddo


        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, &
                           Me%WorkSize%JLB, Me%WorkSize%JUB, &
                           Me%WorkSize%KLB-1, Me%WorkSize%KUB,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructSalt - ModuleHYCOMFormat - ERR06'


        call HDF5WriteData   (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical", "m",      &
                              Array3D = SZZ, OutputNumber = 1,                      &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructSalt - ModuleHYCOMFormat - ERR07'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",           &
                              Array2D = Me%Bathymetry,                          &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructSalt - ModuleHYCOMFormat - ERR09'

   

        deallocate(SZZ, Layers, WorkLayers, LayersT, XX, YY, Latitude, Longitude)


    end subroutine ConstructGrid

    

    subroutine ConstructSalt
        
        !Local-----------------------------------------------------------------
        integer                                 :: stat, STAT_CALL, i, j, k, Salt_id
        logical                                 :: exist
        real                                    :: FillValue
        character (Len = StringLength)          :: PropName
        real, dimension(:,:,:), pointer         :: Salt, WorkSalt
        real, dimension(:,:  ), pointer         :: Salt2D
        real, dimension(:,:  ), pointer         :: SaltT

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing salt...'

        !Verifies if file exists
        inquire(file = Me%HYCOMSaltFileName, exist = exist)
        if (.not. exist) then
            write(*,*)'HYCOM NETCDF salt file does not exist'
            stop 'ConstructSalt - ModuleHYCOMFormat - ERR01'
        endif

        stat = nf90_open(Me%HYCOMSaltFileName, NF90_NOWRITE, Me%HYCOMSaltID)
        if (stat /= nf90_noerr) stop 'ConstructSalt - ModuleHYCOMFormat - ERR02'


        stat = nf90_inq_varid(Me%HYCOMSaltID, 'salt', Salt_id)
        if (stat /= nf90_noerr) stop 'ConstructSalt - ModuleHYCOMFormat - ERR03'

        stat =  nf90_get_att(Me%HYCOMSaltID, Salt_id, '_FillValue', FillValue)
        if (stat /= nf90_noerr) stop 'ConstructSalt - ModuleHYCOMFormat - ERR04'

        allocate(Salt(1:Me%RealSize%JUB, 1:Me%RealSize%IUB, 1:Me%RealSize%KUB))

        stat = nf90_get_var(Me%HYCOMSaltID, Salt_id, Salt)
        if (stat /= nf90_noerr) stop 'ConstructSalt - ModuleHYCOMFormat - ERR05'


        allocate(SaltT(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))

        allocate(WorkSalt(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB, 1:Me%RealSize%KUB))
        
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

            Salt2D => Salt(Me%j_left:Me%j_right, Me%i_bottom:Me%i_top, k)

            SaltT = transpose(Salt2D)

            WorkSalt(:,:, Me%WorkSize%KUB - k + 1) = SaltT
            
        enddo


        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, &
                           Me%WorkSize%JLB, Me%WorkSize%JUB, Me%WorkSize%KLB, Me%WorkSize%KUB,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructSalt - ModuleHYCOMFormat - ERR06'

        PropName = GetPropertyName(Salinity_) 


        call HDF5WriteData   (Me%ObjHDF5, trim("/Results/"//PropName), trim(PropName), "psu",   &
                              Array3D = WorkSalt, OutputNumber = 1,                 &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructSalt - ModuleHYCOMFormat - ERR07'

        allocate(Me%WaterPoints3D(Me%WorkSize%ILB:Me%WorkSize%IUB, &
                                  Me%WorkSize%JLB:Me%WorkSize%JUB, &
                                  Me%WorkSize%KLB:Me%RealSize%KUB))


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

            if(abs(WorkSalt(i,j,k)) >= abs(FillValue)/2.) then

                Me%WaterPoints3D(i,j,k) = 0

            else

                Me%WaterPoints3D(i,j,k) = 1

            end if

        enddo
        enddo
        enddo

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints", "-",          &
                              Array3D = Me%WaterPoints3D,                               &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructSalt - ModuleHYCOMFormat - ERR08'

        deallocate(Me%WaterPoints3D)
        
        deallocate(WorkSalt, SaltT, Salt)

    end subroutine ConstructSalt


    subroutine ConstructUvel
        
        !Local-----------------------------------------------------------------
        integer                                 :: stat, STAT_CALL, i, j, k, Uvel_id, f1, f2, ff
        logical                                 :: exist
        real                                    :: FillValue
        character (Len = StringLength)          :: PropName
        real, dimension(:,:,:), pointer         :: Uvel, WorkUvel, Uvelcenter
        real, dimension(:,:  ), pointer         :: Uvel2D
        real, dimension(:,:  ), pointer         :: UvelT

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing uvel...'

        !Verifies if file exists
        inquire(file = Me%HYCOMUvelFileName, exist = exist)
        if (.not. exist) then
            write(*,*)'HYCOM NETCDF Uvel file does not exist'
            stop 'ConstructUvel - ModuleHYCOMFormat - ERR01'
        endif

        stat = nf90_open(Me%HYCOMUvelFileName, NF90_NOWRITE, Me%HYCOMUvelID)
        if (stat /= nf90_noerr) stop 'ConstructUvel - ModuleHYCOMFormat - ERR02'


        stat = nf90_inq_varid(Me%HYCOMUvelID, 'uvel', Uvel_id)
        if (stat /= nf90_noerr) stop 'ConstructUvel - ModuleHYCOMFormat - ERR03'

        stat =  nf90_get_att(Me%HYCOMUvelID, Uvel_id, '_FillValue', FillValue)
        if (stat /= nf90_noerr) stop 'ConstructUvel - ModuleHYCOMFormat - ERR04'

        allocate(Uvel(1:Me%RealSize%JUB, 1:Me%RealSize%IUB, 1:Me%RealSize%KUB))

        stat = nf90_get_var(Me%HYCOMUvelID, Uvel_id, Uvel)
        if (stat /= nf90_noerr) stop 'ConstructUvel - ModuleHYCOMFormat - ERR05'


        allocate(UvelT(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB+1))

        allocate(WorkUvel  (Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB+1, 1:Me%RealSize%KUB))
        allocate(Uvelcenter(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB  , 1:Me%RealSize%KUB))
        
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

            Uvel2D => Uvel(Me%j_left:Me%j_right+1, Me%i_bottom:Me%i_top, k)

            UvelT = transpose(Uvel2D)

            WorkUvel(:,:, Me%WorkSize%KUB - k + 1) = UvelT
            
        enddo


        allocate(Me%WaterPoints3D(Me%WorkSize%ILB:Me%WorkSize%IUB  , &
                                  Me%WorkSize%JLB:Me%WorkSize%JUB+1, &
                                  Me%WorkSize%KLB:Me%RealSize%KUB))


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB+1
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

            if(abs(WorkUvel(i,j,k)) >= abs(FillValue)/2.) then

                Me%WaterPoints3D(i,j,k) = 0

            else

                Me%WaterPoints3D(i,j,k) = 1

            end if

        enddo
        enddo
        enddo

!        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3DU", "psu",        &
!                              Array3D = Me%WaterPoints3D,                         &
!                              STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'ConstructUvel - ModuleHYCOMFormat - ERR08'

        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, &
                           Me%WorkSize%JLB, Me%WorkSize%JUB, Me%WorkSize%KLB, Me%WorkSize%KUB,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructUvel - ModuleHYCOMFormat - ERR06'

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

            f1 = Me%WaterPoints3D(i,j,k)
            f2 = Me%WaterPoints3D(i,j+1,k)
            ff = f1 + f2

            if (ff > 0) then

                UvelCenter(i, j, k) = (WorkUvel(i,j,k)*f1 + WorkUvel(i,j+1,k)*f2) / ff

            else

                UvelCenter(i, j, k) = FillValueReal

            end if

        enddo
        enddo
        enddo

        PropName = GetPropertyName(VelocityU_) 


        call HDF5WriteData   (Me%ObjHDF5, trim("/Results/"//PropName), trim(PropName), "m/s",   &
                              Array3D = UvelCenter, OutputNumber = 1,                 &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructUvel - ModuleHYCOMFormat - ERR07'

        

        deallocate(WorkUvel, UvelT, Uvel, UvelCenter, Me%WaterPoints3D)

    end subroutine ConstructUvel


    subroutine ConstructVvel
        
        !Local-----------------------------------------------------------------
        integer                                 :: stat, STAT_CALL, i, j, k, Vvel_id, f1, f2, ff
        logical                                 :: exist
        real                                    :: FillValue
        character (Len = StringLength)          :: PropName
        real, dimension(:,:,:), pointer         :: Vvel, WorkVvel, VvelCenter
        real, dimension(:,:  ), pointer         :: Vvel2D
        real, dimension(:,:  ), pointer         :: VvelT

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing vvel...'

        !Verifies if file exists
        inquire(file = Me%HYCOMVvelFileName, exist = exist)
        if (.not. exist) then
            write(*,*)'HYCOM NETCDF Vvel file does not exist'
            stop 'ConstructVvel - ModuleHYCOMFormat - ERR01'
        endif

        stat = nf90_open(Me%HYCOMVvelFileName, NF90_NOWRITE, Me%HYCOMVvelID)
        if (stat /= nf90_noerr) stop 'ConstructVvel - ModuleHYCOMFormat - ERR02'


        stat = nf90_inq_varid(Me%HYCOMVvelID, 'vvel', Vvel_id)
        if (stat /= nf90_noerr) stop 'ConstructVvel - ModuleHYCOMFormat - ERR03'

        stat =  nf90_get_att(Me%HYCOMVvelID, Vvel_id, '_FillValue', FillValue)
        if (stat /= nf90_noerr) stop 'ConstructVvel - ModuleHYCOMFormat - ERR04'

        allocate(Vvel(1:Me%RealSize%JUB, 1:Me%RealSize%IUB, 1:Me%RealSize%KUB))

        stat = nf90_get_var(Me%HYCOMVvelID, Vvel_id, Vvel)
        if (stat /= nf90_noerr) stop 'ConstructVvel - ModuleHYCOMFormat - ERR05'


        allocate(VvelT(Me%WorkSize%ILB:Me%WorkSize%IUB+1, Me%WorkSize%JLB:Me%WorkSize%JUB))

        allocate(WorkVvel  (Me%WorkSize%ILB:Me%WorkSize%IUB+1, Me%WorkSize%JLB:Me%WorkSize%JUB, 1:Me%RealSize%KUB))
        allocate(VvelCenter(Me%WorkSize%ILB:Me%WorkSize%IUB  , Me%WorkSize%JLB:Me%WorkSize%JUB, 1:Me%RealSize%KUB))
        
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

            Vvel2D => Vvel(Me%j_left:Me%j_right, Me%i_bottom:Me%i_top+1, k)

            VvelT = transpose(Vvel2D)

            WorkVvel(:,:, Me%WorkSize%KUB - k + 1) = VvelT
            
        enddo

        allocate(Me%WaterPoints3D(Me%WorkSize%ILB:Me%WorkSize%IUB+1, &
                                  Me%WorkSize%JLB:Me%WorkSize%JUB,   &
                                  Me%WorkSize%KLB:Me%RealSize%KUB))


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB+1
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

            if(abs(WorkVvel(i,j,k)) >= abs(FillValue)/2.) then

                Me%WaterPoints3D(i,j,k) = 0

            else

                Me%WaterPoints3D(i,j,k) = 1

            end if

        enddo
        enddo
        enddo


        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, &
                           Me%WorkSize%JLB, Me%WorkSize%JUB, Me%WorkSize%KLB, Me%WorkSize%KUB,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructVvel - ModuleHYCOMFormat - ERR06'

!        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3DV", "-",               &
!                              Array3D = Me%WaterPoints3D,                               &
!                              STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'ConstructVvel - ModuleHYCOMFormat - ERR08'

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

            f1 = Me%WaterPoints3D(i,j,k)
            f2 = Me%WaterPoints3D(i+1,j,k)
            ff = f1 + f2

            if (ff > 0) then

                VvelCenter(i, j, k) = (WorkVvel(i,j,k)*f1 + WorkVvel(i+1,j,k)*f2) / ff

            else
                
                VvelCenter(i, j, k) = FillValueReal

            end if

        enddo
        enddo
        enddo




        PropName = GetPropertyName(VelocityV_) 


        call HDF5WriteData   (Me%ObjHDF5, trim("/Results/"//PropName), trim(PropName), "m/s",   &
                              Array3D = VvelCenter, OutputNumber = 1,                 &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructVvel - ModuleHYCOMFormat - ERR07'


        
        deallocate(WorkVvel, VvelT, Vvel, Me%WaterPoints3D, VvelCenter)

    end subroutine ConstructVvel


    subroutine ConstructUbaro
        
        !Local-----------------------------------------------------------------
        integer                                 :: stat, STAT_CALL, i, j, Ubaro_id, f1, f2, ff
        logical                                 :: exist
        real                                    :: FillValue
        character (Len = StringLength)          :: PropName
        real, dimension(:,:  ), pointer         :: Ubaro, WorkUbaro, UbaroCenter
        real, dimension(:,:  ), pointer         :: Ubaro2D
        real, dimension(:,:  ), pointer         :: UbaroT

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing ubaro...'

        !Verifies if file exists
        inquire(file = Me%HYCOMUbaroFileName, exist = exist)
        if (.not. exist) then
            write(*,*)'HYCOM NETCDF Ubaro file does not exist'
            stop 'ConstructUbaro - ModuleHYCOMFormat - ERR01'
        endif

        stat = nf90_open(Me%HYCOMUbaroFileName, NF90_NOWRITE, Me%HYCOMUbaroID)
        if (stat /= nf90_noerr) stop 'ConstructUbaro - ModuleHYCOMFormat - ERR02'


        stat = nf90_inq_varid(Me%HYCOMUbaroID, 'ubaro', Ubaro_id)
        if (stat /= nf90_noerr) stop 'ConstructUbaro - ModuleHYCOMFormat - ERR03'

        stat =  nf90_get_att(Me%HYCOMUbaroID, Ubaro_id, '_FillValue', FillValue)
        if (stat /= nf90_noerr) stop 'ConstructUbaro - ModuleHYCOMFormat - ERR04'

        allocate(Ubaro(1:Me%RealSize%JUB, 1:Me%RealSize%IUB))

        stat = nf90_get_var(Me%HYCOMUbaroID, Ubaro_id, Ubaro)
        if (stat /= nf90_noerr) stop 'ConstructUbaro - ModuleHYCOMFormat - ERR05'


        allocate(UbaroT(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB+1))

        allocate(WorkUbaro  (Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB+1))
        allocate(UbaroCenter(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
        
        Ubaro2D => Ubaro(Me%j_left:Me%j_right+1, Me%i_bottom:Me%i_top)

        UbaroT = transpose(Ubaro2D)

        WorkUbaro(:,:) = UbaroT
            

        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, &
                           Me%WorkSize%JLB, Me%WorkSize%JUB,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructUbaro - ModuleHYCOMFormat - ERR06'


        allocate(Me%WaterPoints2D(Me%WorkSize%ILB:Me%WorkSize%IUB, &
                                  Me%WorkSize%JLB:Me%WorkSize%JUB+1))


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB+1
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if(abs(WorkUbaro(i,j)) >= abs(FillValue)/2.) then

                Me%WaterPoints2D(i,j) = 0

            else

                Me%WaterPoints2D(i,j) = 1

            end if

        enddo
        enddo

!        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2DU", "-",               &
!                              Array2D = Me%WaterPoints2D,                               &
!                              STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'ConstructUbaro - ModuleHYCOMFormat - ERR08'

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            f1 = Me%WaterPoints2D(i,j)
            f2 = Me%WaterPoints2D(i,j+1)
            ff = f1 + f2

            if (ff > 0) then

                UbaroCenter(i, j) = (WorkUbaro(i,j)*f1 + WorkUbaro(i,j+1)*f2) / ff

            else
                
                UbaroCenter(i, j) = FillValueReal

            end if

        enddo
        enddo


        PropName = GetPropertyName(BarotropicVelocityU_) 


        call HDF5WriteData   (Me%ObjHDF5, trim("/Results/"//PropName), trim(PropName), "m/s",   &
                              Array2D = UbaroCenter, OutputNumber = 1,                 &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructUbaro - ModuleHYCOMFormat - ERR07'        

        deallocate(WorkUbaro, UbaroT, Ubaro, Me%WaterPoints2D,UbaroCenter)

    end subroutine ConstructUbaro


    subroutine ConstructVbaro
        
        !Local-----------------------------------------------------------------
        integer                                 :: stat, STAT_CALL, i, j, Vbaro_id, f1, f2, ff
        logical                                 :: exist
        real                                    :: FillValue
        character (Len = StringLength)          :: PropName
        real, dimension(:,:  ), pointer         :: Vbaro, WorkVbaro, VbaroCenter
        real, dimension(:,:  ), pointer         :: Vbaro2D
        real, dimension(:,:  ), pointer         :: VbaroT

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing vbaro...'

        !Verifies if file exists
        inquire(file = Me%HYCOMVbaroFileName, exist = exist)
        if (.not. exist) then
            write(*,*)'HYCOM NETCDF Vbaro file does not exist'
            stop 'ConstructVbaro - ModuleHYCOMFormat - ERR01'
        endif

        stat = nf90_open(Me%HYCOMVbaroFileName, NF90_NOWRITE, Me%HYCOMVbaroID)
        if (stat /= nf90_noerr) stop 'ConstructVbaro - ModuleHYCOMFormat - ERR02'


        stat = nf90_inq_varid(Me%HYCOMVbaroID, 'vbaro', Vbaro_id)
        if (stat /= nf90_noerr) stop 'ConstructVbaro - ModuleHYCOMFormat - ERR03'

        stat =  nf90_get_att(Me%HYCOMVbaroID, Vbaro_id, '_FillValue', FillValue)
        if (stat /= nf90_noerr) stop 'ConstructVbaro - ModuleHYCOMFormat - ERR04'

        allocate(Vbaro(1:Me%RealSize%JUB, 1:Me%RealSize%IUB))

        stat = nf90_get_var(Me%HYCOMVbaroID, Vbaro_id, Vbaro)
        if (stat /= nf90_noerr) stop 'ConstructVbaro - ModuleHYCOMFormat - ERR05'


        allocate(VbaroT(Me%WorkSize%ILB:Me%WorkSize%IUB+1, Me%WorkSize%JLB:Me%WorkSize%JUB))

        allocate(WorkVbaro  (Me%WorkSize%ILB:Me%WorkSize%IUB+1, Me%WorkSize%JLB:Me%WorkSize%JUB))

        allocate(VbaroCenter(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
        
        Vbaro2D => Vbaro(Me%j_left:Me%j_right, Me%i_bottom:Me%i_top+1)

        VbaroT = transpose(Vbaro2D)

        WorkVbaro(:,:) = VbaroT
        



        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, &
                           Me%WorkSize%JLB, Me%WorkSize%JUB,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructVbaro - ModuleHYCOMFormat - ERR06'


        allocate(Me%WaterPoints2D(Me%WorkSize%ILB:Me%WorkSize%IUB+1,                    &
                                  Me%WorkSize%JLB:Me%WorkSize%JUB))


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB + 1

            if(abs(WorkVbaro(i,j)) >= abs(FillValue)/2.)then

                Me%WaterPoints2D(i,j) = 0

            else

                Me%WaterPoints2D(i,j) = 1

            end if

        enddo
        enddo

!        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2DV", "-",               &
!                              Array2D = Me%WaterPoints2D,                               &
!                              STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'ConstructVbaro - ModuleHYCOMFormat - ERR08'


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            f1 = Me%WaterPoints2D(i  ,j)
            f2 = Me%WaterPoints2D(i+1,j)
            ff = f1 + f2

            if (ff > 0) then

                VbaroCenter(i, j) = (WorkVbaro(i,j)*f1 + WorkVbaro(i+1,j)*f2) / ff

            else
                
                VbaroCenter(i, j) = FillValueReal

            end if

        enddo
        enddo


       PropName = GetPropertyName(BarotropicVelocityV_) 


        call HDF5WriteData   (Me%ObjHDF5, trim("/Results/"//PropName), trim(PropName), "m/s",   &
                              Array2D = VbaroCenter, OutputNumber = 1,                 &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructVbaro - ModuleHYCOMFormat - ERR07'

 
        

        deallocate(WorkVbaro, VbaroT, Vbaro, Me%WaterPoints2D, VbaroCenter)

    end subroutine ConstructVbaro


    subroutine ConstructSSH
        
        !Local-----------------------------------------------------------------
        integer                                 :: stat, STAT_CALL, i, j, SSH_id
        logical                                 :: exist
        real                                    :: FillValue
        character (Len = StringLength)          :: PropName
        real, dimension(:,:  ), pointer         :: SSH, WorkSSH
        real, dimension(:,:  ), pointer         :: SSH2D
        real, dimension(:,:  ), pointer         :: SSHT

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing ssh...'

        !Verifies if file exists
        inquire(file = Me%HYCOMSSHFileName, exist = exist)
        if (.not. exist) then
            write(*,*)'HYCOM NETCDF SSH file does not exist'
            stop 'ConstructSSH - ModuleHYCOMFormat - ERR01'
        endif

        stat = nf90_open(Me%HYCOMSSHFileName, NF90_NOWRITE, Me%HYCOMSSHID)
        if (stat /= nf90_noerr) stop 'ConstructSSH - ModuleHYCOMFormat - ERR02'


        stat = nf90_inq_varid(Me%HYCOMSSHID, 'ssh', SSH_id)
        if (stat /= nf90_noerr) stop 'ConstructSSH - ModuleHYCOMFormat - ERR03'

        stat =  nf90_get_att(Me%HYCOMSSHID, SSH_id, '_FillValue', FillValue)
        if (stat /= nf90_noerr) stop 'ConstructSSH - ModuleHYCOMFormat - ERR04'

        allocate(SSH(1:Me%RealSize%JUB, 1:Me%RealSize%IUB))

        stat = nf90_get_var(Me%HYCOMSSHID, SSH_id, SSH)
        if (stat /= nf90_noerr) stop 'ConstructSSH - ModuleHYCOMFormat - ERR05'


        allocate(SSHT(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))

        allocate(WorkSSH(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
        
        SSH2D => SSH(Me%j_left:Me%j_right, Me%i_bottom:Me%i_top)

        SSHT = transpose(SSH2D)

        WorkSSH(:,:) = SSHT
        



        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, &
                           Me%WorkSize%JLB, Me%WorkSize%JUB,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructSSH - ModuleHYCOMFormat - ERR06'


        PropName = GetPropertyName(WaterLevel_) 


        call HDF5WriteData   (Me%ObjHDF5, trim("/Results/"//PropName), trim(PropName), "m",   &
                              Array2D = WorkSSH, OutputNumber = 1,                 &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructSSH - ModuleHYCOMFormat - ERR07'

        allocate(Me%WaterPoints2D(Me%WorkSize%ILB:Me%WorkSize%IUB, &
                                  Me%WorkSize%JLB:Me%WorkSize%JUB))


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if(abs(WorkSSH(i,j)) >= abs(FillValue)/2.) then

                Me%WaterPoints2D(i,j) = 0

            else

                Me%WaterPoints2D(i,j) = 1

            end if

        enddo
        enddo

!        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2DZ", "-",         &
!                              Array2D = Me%WaterPoints2D,                         &
!                              STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'ConstructSSH - ModuleHYCOMFormat - ERR08'
        
        deallocate(WorkSSH, SSHT, SSH, Me%WaterPoints2D)

    end subroutine ConstructSSH


    subroutine ConstructTemp
        
        !Local-----------------------------------------------------------------
        integer                                 :: stat, STAT_CALL, i, j, k, Temp_id
        logical                                 :: exist
        real                                    :: FillValue
        character (Len = StringLength)          :: PropName
        real, dimension(:,:,:), pointer         :: Temp, WorkTemp
        real, dimension(:,:  ), pointer         :: Temp2D
        real, dimension(:,:  ), pointer         :: TempT

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing temperature...'

        !Verifies if file exists
        inquire(file = Me%HYCOMTempFileName, exist = exist)
        if (.not. exist) then
            write(*,*)'HYCOM NETCDF Temp file does not exist'
            stop 'ConstructTemp - ModuleHYCOMFormat - ERR01'
        endif

        stat = nf90_open(Me%HYCOMTempFileName, NF90_NOWRITE, Me%HYCOMTempID)
        if (stat /= nf90_noerr) stop 'ConstructTemp - ModuleHYCOMFormat - ERR02'

        stat = nf90_inq_varid(Me%HYCOMTempID, 'temp', Temp_id)
        if (stat /= nf90_noerr) stop 'ConstructTemp - ModuleHYCOMFormat - ERR04'

        stat =  nf90_get_att(Me%HYCOMTempID, Temp_id, '_FillValue', FillValue)
        if (stat /= nf90_noerr) stop 'ConstructTemp - ModuleHYCOMFormat - ERR04'

        allocate(Temp(1:Me%RealSize%JUB, 1:Me%RealSize%IUB, 1:Me%RealSize%KUB))

        stat = nf90_get_var(Me%HYCOMTempID, Temp_id, Temp)
        if (stat /= nf90_noerr) stop 'GetVariable - ModuleHYCOMFormat - ERR01'


        allocate(TempT(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))

        allocate(WorkTemp(Me%WorkSize%ILB:Me%WorkSize%IUB, &
                          Me%WorkSize%JLB:Me%WorkSize%JUB, &
                          Me%WorkSize%KLB:Me%WorkSize%KUB))
        
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

            Temp2D => Temp(Me%j_left:Me%j_right, Me%i_bottom:Me%i_top, k)

            TempT = transpose(Temp2D)

            WorkTemp(:,:,Me%WorkSize%KUB - k + 1) = TempT

        enddo



        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, &
                           Me%WorkSize%JLB, Me%WorkSize%JUB, 1, Me%RealSize%KUB,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructTemp - ModuleHYCOMFormat - ERR04'

        PropName = GetPropertyName(Temperature_) 

        call HDF5WriteData   (Me%ObjHDF5, trim("/Results/"//PropName), trim(PropName), "ºC",   &
                              Array3D = WorkTemp, OutputNumber = 1,                         &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructTemp - ModuleHYCOMFormat - ERR04'


        allocate(Me%WaterPoints3D(Me%WorkSize%ILB:Me%WorkSize%IUB, &
                                  Me%WorkSize%JLB:Me%WorkSize%JUB, &
                                  Me%WorkSize%KLB:Me%RealSize%KUB))


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

            if(abs(WorkTemp(i,j,k)) >= abs(FillValue)/2.) then

                Me%WaterPoints3D(i,j,k) = 0

            else

                Me%WaterPoints3D(i,j,k) = 1

            end if

        enddo
        enddo
        enddo

!        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3Dtemp", "-",            &
!                              Array3D = Me%WaterPoints3D,                               &
!!                              STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'ConstructSalt - ModuleHYCOMFormat - ERR08'

        deallocate(Me%WaterPoints3D)

        deallocate(WorkTemp, TempT, Temp)

    end subroutine ConstructTemp

    !------------------------------------------------------------------------
    
    subroutine KillHYCOMFormat
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillHYCOMFormat - ModuleHYCOMFormat - ERR02'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillHYCOMFormat - ModuleHYCOMFormat - ERR03'

        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillHYCOMFormat

    !--------------------------------------------------------------------------

end module ModuleHYCOMFormat
