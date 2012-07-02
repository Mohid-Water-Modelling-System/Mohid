!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : ConvertToHDF5
! MODULE        : Aladin Format
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Januray 2008
! REVISION      : Guillaume Riflet
! DESCRIPTION   : Module to convert Aldin files into HDF5 format.
!                 Optimized for CVF6.6. Written in ANSI FORTRAN 95
!
!------------------------------------------------------------------------------
!DataFile
!
!   OUTPUT_GRID_FILENAME        : char              -           !Path to grid data file generated from Aladin file
!   OUTPUTFILENAME              : char              -           !Path to HDF5 file generated from Aladin file
!   INPUT_GRID_FILENAME         : char              -           !Path to netcdf griddata file
!
!   <<begin_input_files>>
!   ALADIN_BULKCLOUD_OPASYMP_19723_20088.nc
!   ALADIN_BULKHUMI_OPASYMP_19723_20088.nc
!   ...
!   ... (see below for available fields)
!   <<end_input_files>>

!ALADIN_BULKCLOUD_OPASYMP_19723_20088 {
!dimensions:
!        time_counter = UNLIMITED ; // (2928 currently)
!        lat = 245 ;
!        lon = 251 ;
!variables:
!        float soclotot(time_counter, lat, lon) ;
!                soclotot:long_name = "Couverture nuageuse" ;
!                soclotot:unit = "%" ;
!                soclotot:missing_value = 0.f ;
!        double time_counter(time_counter) ;
!                time_counter:units = "seconds since 1950-01-01 00:00:00" ;
!                time_counter:calendar = "gregorian" ;
!                time_counter:title = "Time" ;
!                time_counter:long_name = "Time axis" ;
!                time_counter:time_origin = "1950-01-01 00:00:00" ;
!        float lat(lat) ;
!                lat:units = "degrees_north" ;
!                lat:valid_min = 0.f ;
!                lat:valid_max = 0.f ;
!                lat:long_name = "Latitude" ;
!                lat:nav_model = "Default grid" ;
!        float lon(lon) ;
!                lon:units = "degrees_east" ;
!                lon:valid_min = 0.f ;
!                lon:valid_max = 0.f ;
!                lon:long_name = "Longitude" ;
!                lon:nav_model = "Default grid" ;
!}

!--- Available properties in Aladin output file ---

!Aladin:
!variables:
!        float soclotot(time_counter, lat, lon) ;
!                soclotot:long_name = "Couverture nuageuse" ;
!                soclotot:unit = "%" ;
!                soclotot:missing_value = 0.f ;
!Mohid: cloud cover

!Aladin:
!variables:
!        float sohumrel(time_counter, lat, lon) ;
!                sohumrel:long_name = "Humidite specifique" ;
!                sohumrel:unit = "%" ;
!                sohumrel:missing_value = 0.f ;
!Mohid: relative humidity

!Aladin:
!       float sofluxir(time_counter, lat, lon) ;
!                sofluxir:long_name = "Flux IR" ;
!                sofluxir:unit = "Wm-2" ;
!Mohid: Downward long wave radiation ---> probably it's this.
!Mohid: #####ERROR:non-solar flux changes#####

!Aladin:
!        float sosspres(time_counter, lat, lon) ;
!                sosspres:long_name = "Pression de surface" ;
!                sosspres:unit = "Pa" ;
!                sosspres:missing_value = 0.f ;
!Mohid: mean sea level pressure OR atmospheric pressure

!Aladin:
!variables:
!        float sosolarf(time_counter, lat, lon) ;
!                sosolarf:long_name = "Flux solaire" ;
!                sosolarf:unit = "Wm-2" ;
!                sosolarf:missing_value = 0.f ;
!Mohid: solar radiation

!Aladin:
!        float sotemair(time_counter, lat, lon) ;
!                sotemair:long_name = "Temperature air 2m" ;
!                sotemair:unit = "Kelvin" ;
!Mohid: air temperature

!Aladin:
!variables:
!        float sowinmod(time_counter, lat, lon) ;
!                sowinmod:long_name = "Module du vent +á 10m" ;
!                sowinmod:unit = "ms-1" ;
!                sowinmod:missing_value = 0.f ;
!Mohid: wind modulus

!Aladin:
!variables:
!        float sowaprec(time_counter, lat, lon) ;
!                sowaprec:unit = "mmday-1" ;
!                sowaprec:missing_value = 0.f ;
!Mohid: precipitation

!Aladin:
!variables:
!        float sowindu10(time_counter, lat, lon) ;
!                sowindu10:long_name = "Vitesse du vent longitudinale 10m" ;
!                sowindu10:unit = "ms-1" ;
!                sowindu10:missing_value = 0.f ;
!Mohid: wind velocity X

!Aladin
!variables:
!        float sowindv10(time_counter, lat, lon) ;
!                sowindv10:long_name = "Vitesse du vent meridienne 10m" ;
!                sowindv10:unit = "ms-1" ;
!                sowindv10:missing_value = 0.f ;
!Mohid: wind velocity Y

!Aladin:
!variables:
!        float sozotaux(time_counter, lat, lon) ;
!                sozotaux:unit = "Pa" ;
!                sozotaux:missing_value = 0.f ;
!Mohid: wind stress X

!Aladin:
!variables:
!        float sometauy(time_counter, lat, lon) ;
!                sometauy:unit = "Pa" ;
!                sometauy:missing_value = 0.f ;
!Mohid: wind stress Y
Module ModuleAladinFormat

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalGrid
    use ModuleHorizontalMap

#ifdef _USE_NIX
    ! Manages NetCDF files
    use netcdf
#else
    use netcdf90
#endif
    
    implicit none

    private 
    
    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertAladinFormat
    private ::      ReadOptions
    private ::      OpenAndReadAladinFile
    private ::          OpenAndReadGridAladin
    private ::          ConstructGrid
    private ::      WriteHDFALadinFile
    private ::          CheckName
    private ::          Open_HDF5_OutPut_File
    private ::          Write_HDF5_Grid_Data
    private ::          OpenAndReadAladinFields
    private ::              ReadAladinFile
    private ::                WriteHDF5Field
    private ::          Close_HDF5_OutPut_File
    private ::      ComputeRelativeHumidity
    private ::          AllocateVariables
    private ::          GetVariables
    private ::          CalculateRelativeHumidity
    private ::          WriteRelativeHumidity
    private ::          ClearVariables
    private ::      KillAladinFormat
    
    !Parameters---------------------------------------------------------------
    character(LEN = StringLength), parameter    :: input_files_begin   = '<<begin_input_files>>'
    character(LEN = StringLength), parameter    :: input_files_end     = '<<end_input_files>>'

    !Types---------------------------------------------------------------------
    
    private :: T_Date
    type       T_Date
        type(T_Time)                            :: Date
        type(T_Date), pointer                   :: Next
    end type  T_Date


    private :: T_Field
    type       T_Field
        character(len=StringLength)             :: Name
        character(len=StringLength)             :: Units
        integer                                 :: GridLocation
        type(T_Time)                            :: Date
        integer                                 :: nDimensions
        real, dimension(:,:),       pointer     :: Values2D
        integer                                 :: OutputNumber         = 1
        type(T_Size2D)                          :: Size, WorkSize
        type(T_Field),              pointer     :: Next
    end type  T_Field

    private :: T_AladinFormat
    type T_AladinFormat
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjHorizontalMap     = 0
        integer                                 :: ObjTime              = 0
        integer                                 :: Unit, ClientNumber
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: InputGridFile
        character(len=PathLength)               :: OutputFileName
        integer                                 :: ReadOptionType
        integer                                 :: imax, jmax
        real, dimension(:,:),       pointer     :: Bathymetry
        real, dimension(:  ),       pointer     :: XX, YY
        integer, dimension(:,:),    pointer     :: WaterPoints
        real, dimension(:,:),       pointer     :: Temperature 
        real, dimension(:,:),       pointer     :: Pressure 
        real, dimension(:,:),       pointer     :: SpecificHumidity
        real, dimension(:,:),       pointer     :: RelativeHumidity
        type(T_Size2D)                          :: Size, WorkSize
        type(T_Field),              pointer     :: FirstField         
        type(T_Date),               pointer     :: FirstDate 
        type(T_Time)                            :: RefDateTime    
        type(T_Time), dimension(:), pointer     :: Times  
        integer                                 :: nInstants
        integer, dimension(12)                  :: Instants(1:13) = 0
        logical                                 :: KelvinToCelsius = .false.
        logical                                 :: PercentageToDecimal = .false.
        logical                                 :: ComputeRelativeHumidity = .false.
        logical                                 :: Convert = .false.
    end type

    type(T_AladinFormat), pointer                              :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine ConvertAladinFormat(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !Local-------------------------------------------------------------------
        integer                                 STAT_CALL

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        Me%ClientNumber = ClientNumber

        !The time in Aladin is compute in seconds from 1950/1/1 : 0h:0m:0s (???)
        call SetDate (Me%RefDateTime, Year=1950, Month=1, Day=1, Hour=0, Minute=0, Second=0) 

        call StartComputeTime(Me%ObjTime, Me%RefDateTime, Me%RefDateTime, Me%RefDateTime, DT = 0.0,    &
                                 VariableDT = .false., STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleAladinFormat - ERR02a'

        call ReadOptions

        call OpenAndReadAladinFile

        if (Me%Convert) then

            write(*,*) 'Converting Aladin to Mohid HDF5 format...'

            call WriteHDFALadinFile

            write(*,*) 'Done converting Aladin to Mohid HDF5 format.'

        endif

        if (Me%ComputeRelativeHumidity) then

            write(*,*) 'Computing relative humidity...'

            call ComputeRelativeHumidity

            write(*,*) 'Done computing relative humidity.'

        endif

        call KillAladinFormat

        STAT = SUCCESS_


    end subroutine ConvertAladinFormat

    !------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------

        !Read output griddata filename
        call GetData(Me%GridFileName,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',                             &
                     ClientModule = 'ModuleAladinFormat',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleAladinFormat - ERR10'

        !Read output filename
        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'ModuleAladinFormat',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleAladinFormat - ERR20'

        !Read input Aladin netcdf gridded data file to generate the griddata
        call GetData(Me%InputGridFile,                                              &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'INPUT_GRID_FILENAME',                          &
                     ClientModule = 'ModuleAladinFormat',                         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleAladinFormat - ERR80'

        !Read convert param type bool
        call GetData(Me%Convert,                                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'CONVERT',                                          &
                     default      = .true.,                                             &         
                     ClientModule = 'ModuleAladinFormat',                               &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleAladinFormat - ERR90'

        !Read ComputeRelativeHumidity param bool,
        !Requires that the output file already exists with AirTemperature, Pressure and
        !Specific Humidity.
        call GetData(Me%ComputeRelativeHumidity,                                      &
                     Me%ObjEnterData, iflag,                                          &
                     SearchType   = FromBlock,                                        &
                     keyword      = 'COMPUTE_RELATIVE_HUMIDITY',                      &
                     default      = .true.,                                           &
                     ClientModule = 'ModuleAladinFormat',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleAladinFormat - ERR100'

    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    subroutine OpenAndReadAladinFile

        call OpenAndReadGridAladin

        call ConstructGrid

    end subroutine OpenAndReadAladinFile    
    
    !------------------------------------------------------------------------

    subroutine WriteHDFALadinFile

        !Local-----------------------------------------------------------------
        integer                 :: HDF5_CREATE    

        !Begin----------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        call Open_HDF5_OutPut_File (HDF5_CREATE)

        call Write_HDF5_Grid_Data

        call OpenAndReadAladinFields

        call Close_HDF5_OutPut_File

    end subroutine WriteHDFALadinFile    

    !------------------------------------------------------------------------

    subroutine ComputeRelativeHumidity

        !Local-----------------------------------------------------------------
        integer                             :: HDF5_IO
        integer                             :: ni, iOut
        character(Len=StringLength)         :: MohidName

        call AllocateVariables

d1:     do ni = 1, Me%nInstants

            MohidName = GetPropertyName(RelativeHumidity_)

            iOut = OutputInstants(MohidName)

            call GetHDF5FileAccess  (HDF5_READ = HDF5_IO)
            call Open_HDF5_OutPut_File (HDF5_IO)
            call GetVariables(iOut)   
            call Close_HDF5_OutPut_File
    
            call CalculateRelativeHumidity

            call GetHDF5FileAccess  (HDF5_READWRITE = HDF5_IO)
            call Open_HDF5_OutPut_File (HDF5_IO)
            call WriteRelativeHumidity (iOut)
            call Close_HDF5_OutPut_File

        enddo d1

        call ClearVariables

    end subroutine ComputeRelativeHumidity
    
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------

    subroutine AllocateVariables

        !Local-----------------------------------------------------------------
        allocate(Me%Temperature(Me%WorkSize%ILB:Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB:Me%WorkSize%JUB))

        allocate(Me%Pressure(Me%WorkSize%ILB:Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB:Me%WorkSize%JUB))

        allocate(Me%SpecificHumidity(Me%WorkSize%ILB:Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB:Me%WorkSize%JUB))

        allocate(Me%RelativeHumidity(Me%WorkSize%ILB:Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB:Me%WorkSize%JUB))

        Me%RelativeHumidity = 0.

    end subroutine AllocateVariables
    
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------

    subroutine ReadTime

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: CurrentInstant

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)"Reading time..."
        write(*,*)

        !Gets number of time instants
        call GetHDF5GroupNumberOfItems(Me%ObjHDF5, '/Time', &
                                       Me%nInstants, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadTime - ModuleAladinFormat - ERR01'

        write(*,*)'Number of instants in hdf5 file : ', Me%nInstants

        allocate(Me%Times(1:Me%nInstants))

        do CurrentInstant = 1, Me%nInstants

            Me%Times(CurrentInstant) = HDF5TimeInstant(CurrentInstant)

        end do

    end subroutine ReadTime

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    type(T_Time) function HDF5TimeInstant(Instant)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer             :: TimeVector

        !Begin-----------------------------------------------------------------

        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = Me%ObjHDF5,                       &
                             GroupName      = '/Time',            &
                             Name           = 'Time',                 &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'HDF5TimeInstant - ModuleAladinFormat - ERR01'

        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))

        deallocate(TimeVector)
        nullify   (TimeVector)

    end function HDF5TimeInstant
    
    !--------------------------------------------------------------------------

    !------------------------------------------------------------------------

    subroutine GetVariables (iOut)

        !Arguments------------------------------------------------
        integer                             :: iOut

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL

        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'GetVariables - ModuleAladinFormat - ERR40'
            
        call HDF5ReadData(HDF5ID       = Me%ObjHDF5,                    &
                          GroupName    = "/Results/specific humidity",  &
                          Name         = "specific humidity",           &
                          Array2D      = Me%SpecificHumidity,           &
                          OutputNumber = iOut,                          &
                          STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'GetVariables - ModuleAladinFormat - ERR50'            

        call HDF5ReadData(HDF5ID       = Me%ObjHDF5,                    &
                          GroupName    = "/Results/air temperature",    &
                          Name         = "air temperature",             &
                          Array2D      = Me%Temperature,                &
                          OutputNumber = iOut,                          &
                          STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'GetVariables - ModuleAladinFormat - ERR51'            

        call HDF5ReadData(HDF5ID       = Me%ObjHDF5,                        &
                          GroupName    = "/Results/atmospheric pressure",   &
                          Name         = "atmospheric pressure",            &
                          Array2D      = Me%Pressure,                       &
                          OutputNumber = iOut,                              &
                          STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'GetVariables - ModuleAladinFormat - ERR52'            

    end subroutine GetVariables
    
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------

    subroutine CalculateRelativeHumidity

        !Locals-------------------------------------------------------
        ! mixture(%), dewpoint mixture (%), dewpoint water pressure (Pa), specific humidity (%)
        real                                        :: w, ws, es, q
        integer                                     :: i, j

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        
            q = Me%SpecificHumidity(i,j)

            w = q / (1 - q)

            ! from mm5tograds
            if (Me%Temperature(i,j) .le. 0.) then
                es = 6.11 * EXP (22.514 - 6150./(Me%Temperature(i,j) + AbsoluteZero))
            else
                !es = 6.11 * EXP (19.83 - 5417./(Me%Temperature(i,j) + AbsoluteZero))
                es = 6.112* EXP (17.67*((Me%Temperature(i,j) - 273.15 + AbsoluteZero)    &
                    /(Me%Temperature(i,j) - 29.65 + AbsoluteZero)))               
            endif

            ws = 0.622 * es / ((Me%Pressure(i,j)/100.) - es)

            ! 5% < Rel. Hum. < 100%
            Me%RelativeHumidity(i,j) =  min(max(100. * w / ws,5.), 100.) * 0.01
            !Me%RelativeHumidity(i,j) =  w / ws

        enddo
        enddo

    end subroutine CalculateRelativeHumidity

    !------------------------------------------------------------------------

    !------------------------------------------------------------------------

    subroutine WriteRelativeHumidity (iOut)

        !Arguments------------------------------------------------
        integer                             :: iOut

        !Local----------------------------------------------------
        integer                             :: STAT_CALL
        character(len=20)                   :: PropUnits
        character(len=StringLength)         :: MohidName

        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteRelativeHumidity - ModuleMecatorFormat - ERR02'
        
        MohidName = GetPropertyName(RelativeHumidity_)

        PropUnits = MohidUnits(MohidName)

        call HDF5WriteData (Me%ObjHDF5, "/Results/relative humidity",              &
                            "relative humidity",trim(PropUnits), Array2D = Me%RelativeHumidity,      &
                            OutputNumber = iOut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleAladinFormat - ERR20'

    end subroutine WriteRelativeHumidity

    !------------------------------------------------------------------------

    subroutine ClearVariables

        deallocate(Me%Temperature)

        deallocate(Me%Pressure)

        deallocate(Me%SpecificHumidity)

        deallocate(Me%RelativeHumidity)

    end subroutine ClearVariables
    
    !------------------------------------------------------------------------

    !Here we're fooling converttohdf5 in thinking there's an all-water bathym
    subroutine OpenAndReadGridAladin

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:  ), pointer     :: Aux2D
        real, dimension(:    ), allocatable     :: Aux1DLon, Aux1DLat
        logical                                 :: exist
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: WILB, WIUB, WJLB,WJUB
        integer                                 :: i, j, n
        integer                                 :: ncid, status, dimid


        !Begin----------------------------------------------------------------

        !Verifies if file exists
        inquire(file = Me%InputGridFile, exist = exist)
i1:     if (exist) then

            status=NF90_OPEN(trim(Me%InputGridFile),NF90_NOWRITE,ncid)
            if (status /= nf90_noerr) stop 'OpenAndReadGridAladin - ModuleAladinFormat - ERR10'

            status=NF90_INQ_DIMID(ncid,"lat",dimid)
            if (status /= nf90_noerr) stop 'OpenAndReadGridAladin - ModuleAladinFormat - ERR20'

            status=NF90_INQUIRE_DIMENSION(ncid,dimid,len = Me%imax)
            if (status /= nf90_noerr) stop 'OpenAndReadGridAladin - ModuleAladinFormat - ERR30'

            status=NF90_INQ_DIMID(ncid,"lon",dimid)
            if (status /= nf90_noerr) stop 'OpenAndReadGridAladin - ModuleAladinFormat - ERR40'

            status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%jmax)
            if (status /= nf90_noerr) stop 'OpenAndReadGridAladin - ModuleAladinFormat - ERR50'
          
            !The border cells are not considered because with the available data
            !is not possible to compute the horizontal position of the cells corners
            Me%WorkSize%ILB = 1
            Me%WorkSize%IUB = Me%imax - 2
                
            Me%WorkSize%JLB = 1
            Me%WorkSize%JUB = Me%jmax - 2

            Me%Size%ILB     = Me%WorkSize%ILB - 1
            Me%Size%IUB     = Me%WorkSize%IUB + 1
            Me%Size%JLB     = Me%WorkSize%JLB - 1
            Me%Size%JUB     = Me%WorkSize%JUB + 1

            WILB            = Me%WorkSize%ILB 
            WIUB            = Me%WorkSize%IUB 
            WJLB            = Me%WorkSize%JLB 
            WJUB            = Me%WorkSize%JUB 

            ILB             = Me%Size%ILB
            IUB             = Me%Size%IUB
            JLB             = Me%Size%JLB
            JUB             = Me%Size%JUB

            allocate(Me%WaterPoints (ILB:IUB, JLB:JUB))
            allocate(Me%Bathymetry  (ILB:IUB, JLB:JUB))
            
            !Longitude == XX == j
            !Latitude == YY == i
            allocate(Me%XX          (JLB:JUB))
            allocate(Me%YY          (ILB:IUB))

            allocate(Aux2D(Me%jmax, Me%imax      ))
            allocate(Aux1DLon(Me%jmax ))
            allocate(Aux1DLat(Me%imax ))

            !There's no mask in the Aladin netcdf variables so it's all water!
            do i=WILB,WIUB
            do j=WJLB,WJUB
               Me%WaterPoints(i,j) = WaterPoint
            enddo
            enddo

            !There's no orography in the Aladin outputs so there's no real bathymetry
            do i=WILB,WIUB
            do j=WJLB,WJUB
                Me%Bathymetry(i,j) = 100.
            enddo
            enddo

            status = nf90_inq_varid(ncid, 'lat', n)
            if (status /= nf90_noerr) stop 'OpenAndReadGridAladin - ModuleAladinFormat - ERR100'

            status = NF90_GET_VAR(ncid,n,Aux1DLat)
            if (status /= nf90_noerr) stop 'OpenAndReadGridAladin - ModuleAladinFormat - ERR110'

            !Calculating Latitude on the cell faces
            do i=WILB,WIUB + 1
               Me%YY(i) = (Aux1DLat(i) + Aux1DLat(i + 1))/2.
            enddo
  
            status = nf90_inq_varid(ncid, 'lon', n)
            if (status /= nf90_noerr) stop 'OpenAndReadGridAladin - ModuleAladinFormat - ERR120'

            status = NF90_GET_VAR(ncid,n,Aux1DLon)
            if (status /= nf90_noerr) stop 'OpenAndReadGridAladin - ModuleAladinFormat - ERR130'

            !Calculate Longitude on the cell faces
            do j=WJLB,WJUB + 1
               Me%XX(j) = (Aux1DLon(j) + Aux1DLon(j + 1))/2.
            enddo

            status=NF90_INQ_DIMID(ncid,"time_counter",dimid)
            if (status /= nf90_noerr) stop 'ReadAladinFile - ModuleAladinFormat - ERR20'

            status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%nInstants)
            if (status /= nf90_noerr) stop 'ReadAladinFile - ModuleAladinFormat - ERR30'

           status=NF90_CLOSE(ncid)
            if (status /= nf90_noerr)  stop 'OpenAndReadGridAladin - ModuleAladinFormat - ERR140'

            deallocate(Aux1DLat)
            deallocate(Aux1DLon)
            deallocate(Aux2D)

        else i1

            write (*,*) "The input grid file do not exist : ",trim(Me%InputGridFile)
            stop 'OpenAndReadGridAladin - ModuleAladinFormat - ERR150'

        endif i1

    end subroutine OpenAndReadGridAladin
    
    !------------------------------------------------------------------------

    
    !------------------------------------------------------------------------

    subroutine ConstructGrid
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL !, UnitID, i, j
        type(T_Size2D)              :: WorkSize2D
        
        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing grid...'


        WorkSize2D%ILB = Me%WorkSize%ILB
        WorkSize2D%JLB = Me%WorkSize%JLB

        WorkSize2D%IUB = Me%WorkSize%IUB
        WorkSize2D%JUB = Me%WorkSize%JUB

        call WriteGridData  (FileName       = Me%GridFileName,                              &
                             XX             = Me%XX,                                        &
                             YY             = Me%YY,                                        &
                             COMENT1        = 'Grid Data created from Aladin NetCDF file',  &
                             COMENT2        = trim(Me%InputGridFile),                       &
                             WorkSize       = WorkSize2D,                                   &
                             CoordType      = 4,                                            &
                             Xorig          = 0.,                                           &
                             Yorig          = 0.,                                           &
                             Zone           = 29,                                           &
                             GRID_ANGLE     = 0.,                                           &
                             Latitude       = Me%XX(1),                                     &
                             Longitude      = Me%YY(1),                                     &
                             FillValue      = -99.,                                         &
                             Overwrite      = .true.,                                       &
                             GridData2D_Real= Me%Bathymetry,                                &
                             STAT           = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructGrid - ModuleAladinFormat - ERR10'


        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleAladinFormat - ERR20'


    end subroutine ConstructGrid

    
    !------------------------------------------------------------------------


    !----------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File (HDF5_IO_CODE)

        !Arguments-------------------------------------------------------------
        integer                                     :: HDF5_IO_CODE

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------
       
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_IO_CODE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR01'
        
        
    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine Write_HDF5_Grid_Data

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR02'
        
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR03'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR04'            
   
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR07'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",    &
                              Array2D = Me%WaterPoints,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR08'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR09'

    end subroutine Write_HDF5_Grid_Data

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine Close_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Close_HDF5_OutPut_File - ModuleAladinFormat - ERR60'

    end subroutine Close_HDF5_OutPut_File

    !--------------------------------------------------------------------------

    !------------------------------------------------------------------------

    subroutine OpenAndReadAladinFields

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        character(len=PathLength)               :: InPutFile
        logical                                 :: exist, BlockFound
        integer                                 :: iflag, line, FirstLine, LastLine,    &
                                                   STAT_CALL


        !Begin----------------------------------------------------------------


        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   input_files_begin, input_files_end,                  &
                                   BlockInBlockFound = BlockFound,                      &
                                   FirstLine = FirstLine, LastLine = LastLine,          &
                                   STAT = STAT_CALL)

IS:     if(STAT_CALL .EQ. SUCCESS_) then

BF:         if (BlockFound) then

                do line = FirstLine + 1, LastLine - 1

                    call GetData(InputFile, EnterDataID = Me%ObjEnterData, flag = iflag, &
                                 Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadAladinFields - ModuleAladinFormat - ERR10'

                    inquire(file = InputFile, exist = exist)
       
i1:                 if (exist) then
                                                    
                        call ReadAladinFile (InputFile) 

                    endif i1
                enddo

            else BF

                stop 'OpenAndReadAladinFields - ModuleAladinFormat - ERR20'

            end if BF

            call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OpenAndReadAladinFields - ModuleAladinFormat - ERR30'

        else   IS

            stop 'OpenAndReadAladinFields - ModuleAladinFormat - ERR40'

        end if IS



    end subroutine OpenAndReadAladinFields


    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    subroutine ReadAladinFile(InputFile)
    !
    !Local subsubroutines : 
    !   CheckName, check
    !   OutputInstants, check
    !   GetPropertyName, 
    !   GetPropertyIDNumber
    !   WriteHDF5Field

        !Arguments-------------------------------------------------------------
        character (Len=*)                       :: InputFile
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:  ), pointer       :: Aux3D
        real, dimension(:,:    ), pointer       :: Aux2D
        real(8), dimension(:    ), allocatable     :: AuxSecs
        integer                                 :: ni, n, iOut, i, j
        integer                                 :: ncid, status, dimid
        integer                                 :: nDimensions
        integer                                 :: nDims, nVars, nAtrr, xtype
        character (len=80)                      :: nameAux
        character(Len=StringLength)             :: MohidName
        type (T_Time)                           :: FieldTime


        !Begin----------------------------------------------------------------

        !Verifies if file exists

        status=NF90_OPEN(trim(InputFile),NF90_NOWRITE,ncid)
        if (status /= nf90_noerr) stop 'ReadAladinFile - ModuleAladinFormat - ERR10'

        status=NF90_INQ_DIMID(ncid,"time_counter",dimid)
        if (status /= nf90_noerr) stop 'ReadAladinFile - ModuleAladinFormat - ERR20'

        status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%nInstants)
        if (status /= nf90_noerr) stop 'ReadAladinFile - ModuleAladinFormat - ERR30'

        allocate(AuxSecs(1:Me%nInstants))

        status = nf90_inq_varid(ncid, 'time_counter', n)
        if (status /= nf90_noerr) stop 'ReadAladinFile - ModuleAladinFormat - ERR40'

        status = NF90_GET_VAR(ncid,n,AuxSecs)
        if (status /= nf90_noerr) stop 'ReadAladinFile - ModuleAladinFormat - ERR50'

        status=NF90_INQUIRE(ncid, nDims, nVars, nAtrr)
        if (status /= nf90_noerr) stop 'OpenAndReadAladinFileV1 - ModuleAladinFormat - ERR60'

d0:     do n=1,nVars

            status=NF90_INQUIRE_VARIABLE(ncid, n, nameAux, xtype, nDimensions)
            if (status /= nf90_noerr) stop 'OpenAndReadAladinFileV1 - ModuleAladinFormat - ERR70'

            if (nDimensions == 3) then

                !Allocate memory
                allocate(Aux3D(Me%jmax, Me%imax, Me%nInstants))
                allocate(Aux2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                Aux2D(:,:) = FillValueReal
                
                !Fill memory
                status = NF90_GET_VAR(ncid,n,Aux3D)
                if (status /= nf90_noerr) stop 'ReadAladinFile - ModuleAladinFormat - ERR80'

            elseif (nDimensions == 2 .or. nDimensions == 1) then
                !Grid properties already written
                cycle
            endif


i2:         if (CheckName(nameAux, MohidName)) then
                
d1:             do ni = 1, Me%nInstants

                    iOut = OutputInstants(MohidName)

                    !Aladin returns its results every 3 hours in seconds since 1950
                    FieldTime = Me%RefDateTime + AuxSecs(ni)

                    if (nDimensions == 3) then

                        !The boundary cells are not read
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                            if (Me%WaterPoints(i, j) == WaterPoint) then
                                
                                Aux2D(i, j) = Aux3D(j+1,i+1,ni)

                                if (Me%PercentageToDecimal) then
                                    Aux2D(i, j) = Aux2D(i, j) / 100
                                endif

                                if (Me%KelvinToCelsius) then

                                    if (Aux2D(i,j) .gt. 100) then
                                        Aux2D(i, j) = Aux2D(i, j) - AbsoluteZero
                                    endif

                                endif

                            endif

                        enddo
                        enddo

                        call WriteHDF5Field(FieldTime, MohidName, iOut, Aux2D = Aux2D)

                    endif


                enddo d1

                !Turn it off for next property before exiting "if"
                if (Me%PercentageToDecimal) then
                    Me%PercentageToDecimal = .false.
                endif

                if (Me%KelvinToCelsius) then
                    Me%KelvinToCelsius = .false.
                endif

            else i2

                write (*,*) 'This property name was not defined yet'
                stop        'ReadAladinFile - ModuleAladinFormat - ERR90'        

            endif i2

            if (nDimensions == 3) then
                deallocate(Aux3D)
                deallocate(Aux2D)
            endif

        enddo d0

        deallocate(AuxSecs)

        !Closes Netcdf file handle
        status=NF90_CLOSE(ncid)
        if (status /= nf90_noerr) stop 'ReadAladinFile - ModuleAladinFormat - ERR91'

    end subroutine ReadAladinFile
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    logical function CheckName(AladinName, MohidName)
        
        !Arguments-----------------------------------------------------------
        character(Len=*)                :: AladinName
        character(Len=StringLength)     :: MohidName
        
        !Begin-----------------------------------------------------------------


        select case(trim(AladinName))

            case('soclotot')

                MohidName = GetPropertyName(CloudCover_)
                Me%PercentageToDecimal = .true.
                CheckName = .true.

            case('sohumrel')

                MohidName = GetPropertyName(SpecificHumidity_)
                Me%PercentageToDecimal = .false.
                CheckName = .true.

            case('sofluxir')

                !MohidName = GetPropertyName(NonSolarFlux_)
                MohidName = GetPropertyName(DownwardLongWaveRadiation_)
                CheckName = .true.

            case('sosspres')

                MohidName = GetPropertyName(AtmosphericPressure_)
                CheckName = .true.

            case('sosolarf')

                MohidName = GetPropertyName(SolarRadiation_)
                CheckName = .true.


            case('sotemair')

                MohidName = GetPropertyName(AirTemperature_)
                Me%KelvinToCelsius = .true.
                CheckName = .true.


            case('sowinmod')

                MohidName = GetPropertyName(WindModulus_)
                CheckName = .true.

            case('sowaprec')

                MohidName = GetPropertyName(Precipitation_)
                CheckName = .true.


            case('sozotaux')

                MohidName = GetPropertyName(WindStressX_)
                CheckName = .true.


            case('sometauy')

                MohidName = GetPropertyName(WindStressY_)
                CheckName = .true.


            case('sowindu10')

                MohidName = GetPropertyName(WindVelocityX_)
                CheckName = .true.


            case('sowindv10')

                MohidName = GetPropertyName(WindVelocityY_)
                CheckName = .true.


            case default
                
                CheckName = .false.

        end select


    end function CheckName
    
    
    !--------------------------------------------------------------------------

    integer function OutputInstants(MohidName)
        
        !Arguments-----------------------------------------------------------
        character(Len=StringLength)     :: MohidName
        !Local-----------------------------------------------------------------
        integer                         :: MohidID
        
        !Begin-----------------------------------------------------------------


        MohidID = GetPropertyIDNumber(MohidName)

        select case(MohidID)

            case(CloudCover_)
                
                Me%Instants(1) = Me%Instants(1) + 1
                OutputInstants = Me%Instants(1)

            case(SpecificHumidity_)
                
                Me%Instants(2) = Me%Instants(2) + 1
                OutputInstants = Me%Instants(2)

            case(DownwardLongWaveRadiation_)
                
                Me%Instants(3) = Me%Instants(3) + 1
                OutputInstants = Me%Instants(3)

            case(AtmosphericPressure_)
                
                Me%Instants(4) = Me%Instants(4) + 1
                OutputInstants = Me%Instants(4)

            case(SolarRadiation_)
                
                Me%Instants(5) = Me%Instants(5) + 1
                OutputInstants = Me%Instants(5)

            case(AirTemperature_)

                Me%Instants(6) = Me%Instants(6) + 1
                OutputInstants = Me%Instants(6)

            case(WindModulus_)

                Me%Instants(7) = Me%Instants(7) + 1
                OutputInstants = Me%Instants(7)

            case(Precipitation_)

                Me%Instants(8) = Me%Instants(8) + 1
                OutputInstants = Me%Instants(8)

            case(WindStressX_)

                Me%Instants(9) = Me%Instants(9) + 1
                OutputInstants = Me%Instants(9)

            case(WindStressY_)

                Me%Instants(10) = Me%Instants(10) + 1
                OutputInstants = Me%Instants(10)

            case(WindVelocityX_)

                Me%Instants(11) = Me%Instants(11) + 1
                OutputInstants = Me%Instants(11)

            case(WindVelocityY_)

                Me%Instants(12) = Me%Instants(12) + 1
                OutputInstants = Me%Instants(12)

            case(RelativeHumidity_)
                
                Me%Instants(13) = Me%Instants(13) + 1
                OutputInstants = Me%Instants(13)

            case default

                write(*,*) 'Warning - ModuleAladinFormat - WNG01'

        end select


    end function OutputInstants

    !------------------------------------------------------------------------
    
    
    !------------------------------------------------------------------------
    subroutine WriteHDF5Field(FieldTime, MohidName, iOut, Aux2D, Aux3D)

        !Arguments-------------------------------------------------------------
        type (T_Time)                                   :: FieldTime
        character(Len=StringLength)                     :: MohidName
        integer                                         :: iOut

        real, dimension(:,:  ), pointer, optional       :: Aux2D
        real, dimension(:,:,:), pointer, optional       :: Aux3D

        !Local-----------------------------------------------------------------
        type (T_Time)                                   :: AuxFieldTime
        character(Len=StringLength)                     :: PropUnits, AuxName
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        integer                                         :: STAT_CALL, nItems
        integer                                         :: WorkILB, WorkJLB, WorkKLB
        integer                                         :: WorkIUB, WorkJUB, WorkKUB
        logical                                         :: Exist

        !Begin-----------------------------------------------------------------
        
        !Bounds

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        PropUnits = MohidUnits(MohidName)

        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                     &
                             WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleAladinFormat - ERR10'

        if      (present(Aux2D)) then

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(MohidName),              &
                                 trim(MohidName),trim(PropUnits), Array2D = Aux2D,      &
                                 OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleAladinFormat - ERR20'

        else if (present(Aux3D)) then

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(MohidName),              &
                                 trim(MohidName),trim(PropUnits), Array3D = Aux3D,      &
                                 OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleAladinFormat - ERR30'

        else 

            stop 'WriteHDF5Field - ModuleAladinFormat - ERR40'

        endif

        call GetHDF5GroupExist (Me%ObjHDF5, "/Time", Exist, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleAladinFormat - ERR35'

        if (Exist) then
        
            call GetHDF5GroupNumberOfItems (Me%ObjHDF5,  "/Time", nItems, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleAladinFormat - ERR37'
        
        endif    


        if (.not. Exist .or. iOut > nItems) then

            !Writes current time
            call ExtractDate   (FieldTime, AuxTime(1), AuxTime(2), AuxTime(3),              &
                                           AuxTime(4), AuxTime(5), AuxTime(6))
            TimePtr => AuxTime
            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleAladinFormat - ERR50'

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",        &
                                 Array1D = TimePtr, OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleAladinFormat - ERR60'

            AuxName = MohidName

        else if (iOut == nItems) then

            TimePtr => AuxTime
            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleAladinFormat - ERR100'

            call HDF5ReadData  (Me%ObjHDF5, "/Time", "Time",                               &
                                 Array1D = TimePtr, OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleAladinFormat - ERR110'


            call SetDate   (AuxFieldTime, AuxTime(1), AuxTime(2), AuxTime(3),              &
                                          AuxTime(4), AuxTime(5), AuxTime(6))


            if (FieldTime /= AuxFieldTime) then
!                write(*,*) 'The time instants of property ',trim(MohidName)
!                write(*,*) 'are not consistent with property ',trim(AuxName)
                stop 'WriteHDF5Field - ModuleAladinFormat - ERR120'   
            endif             
            
        endif

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleAladinFormat - ERR90'


    end subroutine WriteHDF5Field
 
    !------------------------------------------------------------------------
 

    !--------------------------------------------------------------------------

    character(Len=StringLength) function MohidUnits(MohidName)
        
        !Arguments-----------------------------------------------------------
        character(Len=StringLength)     :: MohidName
        !Local-----------------------------------------------------------------
        integer                         :: MohidID
        
        !Begin-----------------------------------------------------------------


        MohidID = GetPropertyIDNumber(MohidName)

        select case(MohidID)

            case(CloudCover_, SpecificHumidity_, RelativeHumidity_)
               
                MohidUnits = '%'

            case(DownwardLongWaveRadiation_, SolarRadiation_)
               
                MohidUnits = 'W/m2'

            case(AtmosphericPressure_, WindStressX_, WindStressY_)
               
                MohidUnits = 'Pa'

            case(AirTemperature_)
               
                MohidUnits = 'ºC'

            case(WindModulus_, WindVelocityX_, WindVelocityY_)
               
                MohidUnits = 'm/s'

            case(Precipitation_)
               
                MohidUnits = 'mm/day'

            !Legacy code units
            case(VelocityU_, VelocityV_, BarotropicVelocityU_, BarotropicVelocityV_)
               
                MohidUnits = 'm/s'

            case(Temperature_)
                
                MohidUnits = 'ºC'

            case(Salinity_)
                
                MohidUnits = 'psu'

            case(WaterLevel_)
                
                MohidUnits = 'm'

        end select


    end function MohidUnits
    
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine KillAladinFormat
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillAladinFormat - ModuleAladinFormat - ERR50'
        
!        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'KillAladinFormat - ModuleAladinFormat - ERR60'

        deallocate(Me%XX)
        deallocate(Me%YY)
        deallocate(Me%Bathymetry)
        deallocate(Me%WaterPoints)

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillAladinFormat - ModuleAladinFormat - ERR70'

        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillAladinFormat

    !--------------------------------------------------------------------------


end module ModuleAladinFormat









