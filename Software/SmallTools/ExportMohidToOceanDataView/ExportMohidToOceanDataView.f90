program ExportMohidToOceanDataView

    use ModuleGlobalData
    use ModuleHDF5


    implicit none


    integer                             :: ObjHDF5          = 0
    integer                             :: HDF5_READWRITE, STAT_CALL, i,j,k, n, imax, jmax, kmax, Ninstants
    real,    dimension(:,:,:), pointer  :: Temperature, Vertical
    integer,    dimension(:,:,:), pointer  :: OpenPoints3D
    real,    dimension(:,:  ), pointer  :: ConnectionX, ConnectionY
        real, dimension(6), target          :: AuxTime
        real, dimension(:), pointer         :: TimePtr
    character(Len=4)                        :: ic, jc,Year,Month,day,minute
    character(Len=256)  :: CruiseDate, Station, CruiseName
    character(Len=8)    :: Hour
    real                :: Lat, Long,Depth

    call StartUpMohid("ReadHDF5")

    kmax = 8
    imax = 116
    jmax = 143
    Ninstants = 18

    allocate(Temperature(1:imax,1:jmax,1:kmax))
    allocate(OpenPoints3D(1:imax,1:jmax,1:kmax))
    allocate(Vertical(1:imax,1:jmax,0:kmax))
    allocate(ConnectionX(1:imax+1,1:jmax+1))
    allocate(ConnectionY(1:imax+1,1:jmax+1))

    call GetHDF5FileAccess  (HDF5_READWRITE = HDF5_READWRITE)


    call ConstructHDF5 (ObjHDF5, trim("WaterProperties_2.hdf5"), &
                        HDF5_READWRITE, STAT = STAT_CALL)
    if (STAT_CALL .NE. SUCCESS_) stop 'ReadHDF5 - ERR01'

    !Sets limits for next write operations
    !call HDF5SetLimits   (ObjHDF5, 1, imax+1, 1, jmax+1,           &
    !                      STAT = STAT_CALL)
    !if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR02'


    !call HDF5ReadData   (ObjHDF5, "/Grid", "ConnectionX",                       &
    !                      Array2D = ConnectionX,                                        &
    !                      STAT = STAT_CALL)
    !if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR03'

    !call HDF5ReadData   (ObjHDF5, "/Grid", "ConnectionY",                       &
    !                      Array2D = ConnectionY,                                        &
    !                      STAT = STAT_CALL)
    !if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR04'

    open(34,file='InputLongLatMohid.txt')

    do i=1,imax
    do j=1,jmax
        read(34,*) ConnectionX(i,j), ConnectionY(i,j)
    enddo
    enddo

    close(34)

    open(33,file='outputMohid.txt')

    write(33,'(A178)') 'Cruise;Station;Type;mon/day/yr;hh:mm ; Longitude [degrees_east] ;  Latitude [degrees_north];Depth [m] ;  Temperature [°C] ;'

    do n=1,Ninstants

        call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleHydrodynamic - ERR01'

        TimePtr => AuxTime

        call HDF5ReadData  (ObjHDF5, "/Time", "Time",       &
                             Array1D = TimePtr, OutputNumber = n, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleHydrodynamic - ERR02'


        call HDF5SetLimits  (ObjHDF5,1,imax,1,jmax,1,kmax,STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadHDF5 - ERR02'
    
        call HDF5ReadData(ObjHDF5, "/Results/"//"temperature",  &
                          "temperature", Array3D = Temperature, outputNumber = n, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadHDF5 - ERR03'

        call HDF5ReadData  (ObjHDF5, "/Grid/OpenPoints", "OpenPoints",              &
                             Array3D = OpenPoints3D, OutputNumber = n,         &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModuleHydrodynamic - ERR06'


        call HDF5SetLimits  (ObjHDF5,1,imax,1,jmax,0,kmax,STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadHDF5 - ERR02'

        call HDF5ReadData(ObjHDF5, "/Grid/"//"VerticalZ",  &
                          "Vertical", Array3D = Vertical, outputNumber = n, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadHDF5 - ERR04'

        do i = 1, imax
        do j = 1, jmax
        do k = 1, kmax

            if(OpenPoints3D(i,j,k) == 1)then
                CruiseName = 'Mohid'
                write(ic,'(I4)') i
                write(jc,'(I4)') j
                Station = trim(adjustl(ic))//'_'//trim(adjustl(jc))
                write(Year,'(I4)') int(TimePtr(1))
                write(Month,'(I4)') int(TimePtr(2))
                write(Day,'(I4)') int(TimePtr(3))
                CruiseDate = trim(adjustl(Month))//'/'//trim(adjustl(day))//'/'//trim(adjustl(Year))
                write(Hour,'(I4)') int(TimePtr(4))
                write(Minute,'(I4)') int(TimePtr(5))
                Hour = trim(adjustl(Hour))//':'//trim(adjustl(minute))                             
                Long = (ConnectionX(i,j)+ConnectionX(i+1,j)+ConnectionX(i,j+1)+ConnectionX(i+1,j+1))/4.
                Lat  = (ConnectionY(i,j)+ConnectionY(i+1,j)+ConnectionY(i,j+1)+ConnectionY(i+1,j+1))/4.
                Depth = (Vertical(i,j,k) + Vertical(i,j,k-1))/2. - Vertical(i,j,kmax)
                write(33,'(A30,1A,A30,1A,A1,1A,A10,1A,A5,4(1A,f10.6))') trim(CruiseName), ';',trim(Station),';', 'B',';', &
                trim(CruiseDate), ';',trim(Hour(1:5)), ';',Long, ';',Lat,';',Depth,';',Temperature(i,j,k)


            end if
        enddo
        enddo
        enddo

        exit

    enddo

    call KillHDF5(ObjHDF5, STAT = STAT_CALL)
    if (STAT_CALL .NE. SUCCESS_) stop 'ReadHDF5 - ERR04'

    close(33)

    deallocate(Temperature)


end program ExportMohidToOceanDataView

