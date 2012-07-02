!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : ConvertToHDF5
! MODULE        : IHRadar Format
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2012
! REVISION      : v1
! DESCRIPTION   : Module to convert SeaSonde HF radar ascii files into HDF5 format.
!
!------------------------------------------------------------------------------
!DataFile
!
!   IH_GRID_VERSION             : integer           -           !Grid version from IH radar file.
!   INPUT_GRID_FILENAME         : char              -           !Path to IH radar grid file
!   OUTPUT_GRID_FILENAME        : char              -           !!Path to grid data file generated from seasonde radar file
!   OUTPUT_HDF5_FILENAME        : char              -           !Path to HDF5 file generated from seasonde radar file
!   OUTPUT_NETCDF_FILENAME      : char              -           !Path to Netcdf file generated from seasonde radar file
!   NETCDF_TITLE                : char              -           !Netcdf file global attribute
!   NETCDF_CONVENTION           : char              CF-1.0      !Netcdf file global attribute
!   NETCDF_VERSION              : char              3.6.1       !Netcdf file global attribute
!   NETCDF_HISTORY              : char              -           !Netcdf file global attribute
!   NETCDF_INSTITUTION          : char              Instituto Superior Tecnico           !Netcdf file global attribute
!   NETCDF_REFERENCES           : char              http://www.mohid.com           !Netcdf file global attribute
!   NETCDF_DATE                 : char              -           !Netcdf file global attribute
!
!   <<begin_input_files>>
!   TOTL_IHOC_2012_05_01_1400.tuv
!   TOTL_IHOC_2012_05_01_1500.tuv
!   ...
!   ... (see below for available fields)
!   <<end_input_files>>

!IH Radar Ascii File
!http://websig.hidrografico.pt/www/content/produtos/simoc/TOTL_IHOC_2012_05_01_1400.tuv
!
!%CTF: 1.00
!%FileType: LLUV tots "CurrentMap"
!%LLUVSpec: 1.03  2006 10 05
!%Manufacturer: CODAR Ocean Sensors. SeaSonde
!%Site: IHOC ""
!%TimeStamp: 2012 05 27  14 00 00
!%TimeZone: "UTC" +0.000 0
!%TimeCoverage: 75 Minutes
!%Origin:  38.4333333   -9.8666667
!%GreatCircle: "WGS84" 6378137.000  298.257223562997
!%GeodVersion: "CGEO" 1.50  2006 05 07
!%CombineMethod: 1
!%GridCreatedBy: SeaDisplay 5.1.5
!%GridVersion: 4
!%GridTimeStamp: 4  2011 12 01  00 00 00
!%GridLastModified: 2012 04 09  14 48 04
!%GridAxisOrientation: 0.0 DegNCW
!%GridAxisType: 6
!%GridSpacing: 1.400 km
!%AveragingRadius: 8.000 km
!%DistanceAngularLimit: 20.0
!%CurrentVelocityLimit: 100.0 cm/s
!%% SiteSource # Name  Lat           Lon    Coverage(s) RngStep(km)  Pattern AntBearing(NCW)
!%SiteSource:  1 JLSM   38.6745000   -9.3264167   75.00     0.887 Ideal    269.0
!%SiteSource:  2 EPSM   38.4154667   -9.2166833   75.00     0.887 Meas     286.0
!%TableType: LLUV TOT4
!%TableColumns: 20
!%TableColumnTypes: LOND LATD VELU VELV VFLG UQAL VQAL CQAL XDST YDST RNGE BEAR VELO HEAD S1CN S2CN S3CN S4CN S5CN S6CN 
!%TableRows: 1165
!%TableStart:
!%%   Longitude   Latitude    U comp   V comp  VectorFlag   U StdDev    V StdDev   Covariance  X Distance  Y Distance   Range   Bearing   Velocity  Direction  Site Contributers       
!%%     (deg)       (deg)     (cm/s)   (cm/s)  (GridCode)    Quality     Quality     Quality      (km)        (km)       (km)  (deg NCW)   (cm/s)   (deg NCW)  (#1)(#2)(#3)(#4)(#5)(#6)
!    -9.8506463  38.3702715    1.834  -15.769          0       1.040       3.840       2.560      1.4000     -7.0000    7.1386   168.7     15.880     173.4      9   2   0   0   0   0
!    -9.8506436  38.3828837    2.283  -23.838          0       0.990       4.200       2.650      1.4000     -5.6000    5.7723   166.0     23.950     174.5      7   2   0   0   0   0
!    -9.8506408  38.3954959    2.612  -27.267          0       0.820       3.910       2.220      1.4000     -4.2000    4.4272   161.6     27.390     174.5      7   2   0   0   0   0
!    -9.8506380  38.4081080    2.138  -25.013          0       0.900       4.210       2.520      1.4000     -2.8000    3.1305   153.4     25.100     175.1      9   2   0   0   0   0
!    ...
!    -9.1297477  38.3679545   20.221  -23.394          0       1.050       1.050      -0.820     64.4000     -7.0000   64.7793    96.2     30.920     139.2      9 103   0   0   0   0
!%TableEnd:
!%%
!%ProcessedTimeStamp: 2012 05 27  15 13 17
!%ProcessingTool: "RadialsToCurrents" 10.2.5
!%ProcessingTool: "CheckForRadials" 11.1.5
!%ProcessingTool: "TotalArchiver" 11.0.4
!%End:

!--- Available properties in IH Radar output file ---

!IH-Radar:
!Longitude   Latitude    U comp   V comp  VectorFlag   U StdDev    V StdDev   Covariance  X Distance  Y Distance   Range   Bearing   Velocity  Direction  Site Contributers
!  (deg)       (deg)     (cm/s)   (cm/s)  (GridCode)    Quality     Quality     Quality      (km)        (km)       (km)  (deg NCW)   (cm/s)   (deg NCW)  (#1)(#2)(#3)(#4)(#5)(#6)
!  CurrentVelocityLimit: 100.0 cm/s
!Mohid:
!Longitude   Latitude    Velocity U  Velocity V  U StdDev    V StdDev   Covariance  X Distance  Y Distance   Range   Bearing   Velocity  Direction
!  (º)         (º)          (m/s)       (m/s)                                          (m)          (m)        (m)     (º)       (m/s)     (º)  

Module ModuleIHRadarFormat

    use ModuleGlobalData
    use ModuleFunctions,   only: SetMatrixValue
    use ModuleHDF5
    use ModuleNetcdf
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData,    only: WriteGridData
    use ModuleHorizontalGrid, only: ConstructHorizontalGrid, GetHorizontalGridSize, &
                                    LocateCellPolygons, GetCornersCoordinates,  &
                                    UngetHorizontalGrid, KillHorizontalGrid, WriteHorizontalGrid
    use ModuleHorizontalMap
    
    implicit none

    private 
    
    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertIHRadarFormat
    private ::      ReadOptions
    private ::      CreateGridIHRadarFile
    private ::          LoadIHRadarGrid
    private ::          ConstructGridDataFile
    private ::      WriteIHRadarFile
    private ::          Open_HDF5_OutPut_File
    private ::          OpenAndReadIHRadarFields
    private ::              ReadIHRadarFile
    private ::                WriteHDF5Field
    private ::          Close_HDF5_OutPut_File
    private ::          AllocateVariables
    private ::          ClearVariables
    private ::      KillIHRadarFormat
    
    !Parameters---------------------------------------------------------------
    character(LEN = StringLength), parameter    :: input_files_begin   = '<<begin_input_files>>'
    character(LEN = StringLength), parameter    :: input_files_end     = '<<end_input_files>>'
    character(LEN = StringLength), parameter    :: data_begin   = '%TableStart'
    character(LEN = StringLength), parameter    :: data_end     = '%TableEnd'
    integer, parameter                          :: no_output_       = -99
    integer, parameter                          :: hdf5_       = 1
    integer, parameter                          :: netcdf_     = 2

    !Types---------------------------------------------------------------------
    private :: T_Netcdf_Global_Attributes
    type T_Netcdf_Global_Attributes
        character(len=StringLength)                         :: Title
        character(len=StringLength)                         :: Convention
        character(len=StringLength)                         :: Version
        character(len=StringLength)                         :: History
        character(len=StringLength)                         :: Source
        character(len=StringLength)                         :: Institution
        character(len=StringLength)                         :: References
        integer                                             :: iDate
    end type T_Netcdf_Global_Attributes
    
    private :: T_IHRadarFormat
    type T_IHRadarFormat
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjNetcdf            = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjHorizontalMap     = 0
        integer                                 :: ObjTime              = 0
        integer                                 :: Unit, ClientNumber
        logical                                 :: OutputHDF5
        character(len=PathLength)               :: OutputHDF5FileName
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: InputGridFile
        integer                                 :: GridVersion
        integer                                 :: ReadOptionType
        integer                                 :: imax, jmax
        real, dimension(:,:),       pointer     :: Bathymetry
        real, dimension(:,:),       pointer     :: XX_IE, YY_IE
               
        integer, dimension(:,:),    pointer     :: OpenPoints
        real, dimension(:,:),       pointer     :: velu
        real, dimension(:,:),       pointer     :: velv
        real, dimension(:,:),       pointer     :: velu_stdv
        real, dimension(:,:),       pointer     :: velv_stdv
        real, dimension(:,:),       pointer     :: covariance
        real, dimension(:,:),       pointer     :: xdistance
        real, dimension(:,:),       pointer     :: ydistance
        real, dimension(:,:),       pointer     :: range
        real, dimension(:,:),       pointer     :: bearing
        real, dimension(:,:),       pointer     :: velocity
        real, dimension(:,:),       pointer     :: direction

        logical                                 :: OutputNetcdf
        character(len=PathLength)               :: OutputNetcdfFileName
        type(T_Netcdf_Global_Attributes)        :: OutputNetcdfAttr
                
        type(T_Size2D)                          :: Size, WorkSize
        type(T_Time)                            :: Time
    end type

    type(T_IHRadarFormat), pointer                              :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine ConvertIHRadarFormat(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        Me%ClientNumber = ClientNumber

        !GRiflet: time needs to be rethought completely from scratch ...
        !The time in IHRadar is compute in seconds from 1950/1/1 : 0h:0m:0s (???)
        !call SetDate (Me%RefDateTime, Year=1950, Month=1, Day=1, Hour=0, Minute=0, Second=0) 
        !call StartComputeTime(Me%ObjTime, Me%RefDateTime, Me%RefDateTime, Me%RefDateTime, DT = 0.0,    &
        !                         VariableDT = .false., STAT = STAT_CALL)   
        !if (STAT_CALL /= SUCCESS_) stop 'ConstructGridDataFile - ModuleIHRadarFormat - ERR02a'

        call ReadOptions
        
        if (Me%OutputNetcdf) then
            call ReadNetcdfOptions
        endif

        !Open and load the grid file; then create a new griddata file
        call CreateGridIHRadarFile

        write(*,*) 'Converting IHRadar to Mohid HDF5 format...'

        call AllocateVariables

        !open the list os ascii IH radar files, read and write 
        !in the HDF5 file in MOHID format
        call WriteIHRadarFile
        
        write(*,*) 'Done converting IHRadar to Mohid HDF5 format.'

        call KillIHRadarFormat

        STAT = SUCCESS_

    end subroutine ConvertIHRadarFormat

    !------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------

        !Read output filename
        call GetData(Me%OutputHDF5FileName,                                                &
                     Me%ObjEnterData, iflag,                                           &
                     SearchType   = FromBlock,                                         &
                     keyword      = 'OUTPUTFILENAME',                                 &
                     ClientModule = 'ModuleIHRadarFormat',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleIHRadarFormat - ERR10'
        if (iflag) then
            Me%OutputHDF5 = .true.
        else
            Me%OutputHDF5 = .false.
        endif

        !Read output filename
        call GetData(Me%OutputHDF5FileName,                                                &
                     Me%ObjEnterData, iflag,                                           &
                     SearchType   = FromBlock,                                         &
                     keyword      = 'OUTPUT_HDF5_FILENAME',                            &
                     ClientModule = 'ModuleIHRadarFormat',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleIHRadarFormat - ERR11'
        if (iflag) then
            Me%OutputHDF5 = .true.
        else
            Me%OutputHDF5 = .false.
        endif

        !Read output netcdf filename
        call GetData(Me%OutputNetcdfFileName,                                          &
                     Me%ObjEnterData, iflag,                                           &
                     SearchType   = FromBlock,                                         &
                     keyword      = 'OUTPUT_NETCDF_FILENAME',                          &
                     ClientModule = 'ModuleIHRadarFormat',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleIHRadarFormat - ERR12'        
        if (iflag) then
            Me%OutputNetcdf = .true.
        else
            Me%OutputNetcdf = .false.
        endif
        
        if (.not.Me%OutputNetcdf .or. .not.Me%OutputHDF5) then
            write(*,*) 'Please define one of OUTPUT_NETCDF_FILENAME or OUTPUT_HDF5_FILENAME (aka OUTPUTFILENAME) keywords.'
            stop 'ReadOptions - ModuleIHRadarFormat - ERR13'        
        endif  
        
        !Read output grid filename
        call GetData(Me%GridFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',                                   &
                     ClientModule = 'ModuleIHRadarFormat',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleIHRadarFormat - ERR20'
        if (.not.iflag) then
            write(*,*) 'Please define OUTPUT_GRID_FILENAME keyword'
            stop 'ReadOptions - ModuleIHRadarFormat - ERR21'
        endif

        !Read input IHRadar netcdf gridded data file to generate the griddata
        call GetData(Me%InputGridFile,                                              &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'INPUT_GRID_FILENAME',                          &
                     ClientModule = 'ModuleIHRadarFormat',                         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleIHRadarFormat - ERR80'
        if (.not.iflag) then
            write(*,*) 'Please define INPUT_GRID_FILENAME keyword'
            stop 'ReadOptions - ModuleIHRadarFormat - ERR22'
        endif

        call GetData(Me%GridVersion,                                                &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'IH_GRID_VERSION',                       &
                     ClientModule = 'ModuleIHRadarFormat',                          &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleIHRadarFormat - ERR81'
        if (.not.iflag) then
            write(*,*) 'Please define IH_GRID_VERSION keyword'
            stop 'ReadOptions - ModuleIHRadarFormat - ERR23'
        endif

    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    subroutine ReadNetcdfOptions
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------

        call GetData(Me%OutputNetcdfAttr%Title,                                                &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_TITLE',                                     &
                     ClientModule = 'ModuleIHRadarFormat',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadNetcdfOptions - ModuleIHRadarFormat - ERR04'

        call GetData(Me%OutputNetcdfAttr%Convention,                                           &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_CONVENTION',                                &
                     Default      = 'CF-1.0',                                           &
                     ClientModule = 'ModuleIHRadarFormat',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadNetcdfOptions - ModuleIHRadarFormat - ERR05'

        call GetData(Me%OutputNetcdfAttr%Version,                                              &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_VERSION',                                   &
                     Default      = '3.6.1',                                            &
                     ClientModule = 'ModuleIHRadarFormat',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadNetcdfOptions - ModuleIHRadarFormat - ERR06'

        call GetData(Me%OutputNetcdfAttr%History,                                              &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_HISTORY',                                   &
                     ClientModule = 'ModuleIHRadarFormat',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadNetcdfOptions - ModuleIHRadarFormat - ERR07'
        
        call GetData(Me%OutputNetcdfAttr%Source,                                               &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_SOURCE',                                    &
                     ClientModule = 'ModuleIHRadarFormat',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadNetcdfOptions - ModuleIHRadarFormat - ERR08'

        
        call GetData(Me%OutputNetcdfAttr%Institution,                                          &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_INSTITUTION',                               &
                     Default      = 'Instituto Superior Tecnico',                       &
                     ClientModule = 'ModuleIHRadarFormat',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadNetcdfOptions - ModuleIHRadarFormat - ERR09'
        
        call GetData(Me%OutputNetcdfAttr%References,                                           &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_REFERENCES',                                &
                     Default      = 'http://www.mohid.com',                             &
                     ClientModule = 'ModuleIHRadarFormat',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadNetcdfOptions - ModuleIHRadarFormat - ERR10'

        call GetData(Me%OutputNetcdfAttr%iDate,                                                &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_DATE',                                      &
                     ClientModule = 'ModuleIHRadarFormat',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadNetcdfOptions - ModuleIHRadarFormat - ERR11'
    
    end subroutine ReadNetcdfOptions

    !------------------------------------------------------------------------

    subroutine CreateGridIHRadarFile

        !Load grid from file
        call LoadIHRadarGrid
        
        !Construct fake bathymetry and openpoints and write to output griddata file
        call ConstructGridDataFile

    end subroutine CreateGridIHRadarFile    
    
    !------------------------------------------------------------------------

    subroutine WriteIHRadarFile

        !Begin----------------------------------------------------------------

        !Gets File Access Code
        call Open_Output_File

        call OpenAndReadIHRadarFields

        call Close_Output_File

    end subroutine WriteIHRadarFile    
   
    !------------------------------------------------------------------------

    subroutine AllocateVariables

        !Local-----------------------------------------------------------------
        integer                 :: ILB, IUB, JLB, JUB

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        allocate(Me%velu(ILB:IUB,JLB:JUB))
        allocate(Me%velv(ILB:IUB,JLB:JUB))
        allocate(Me%velu_stdv(ILB:IUB,JLB:JUB))
        allocate(Me%velv_stdv(ILB:IUB,JLB:JUB))
        allocate(Me%covariance(ILB:IUB,JLB:JUB))
        allocate(Me%xdistance(ILB:IUB,JLB:JUB))
        allocate(Me%ydistance(ILB:IUB,JLB:JUB))
        allocate(Me%range(ILB:IUB,JLB:JUB))
        allocate(Me%bearing(ILB:IUB,JLB:JUB))
        allocate(Me%velocity(ILB:IUB,JLB:JUB))
        allocate(Me%direction(ILB:IUB,JLB:JUB))

    end subroutine AllocateVariables

    !------------------------------------------------------------------------

    !------------------------------------------------------------------------

    subroutine ClearVariables

        deallocate(Me%velu)
        deallocate(Me%velv)
        deallocate(Me%velu_stdv)
        deallocate(Me%velv_stdv)
        deallocate(Me%covariance)
        deallocate(Me%xdistance)
        deallocate(Me%ydistance)
        deallocate(Me%range)
        deallocate(Me%bearing)
        deallocate(Me%velocity)
        deallocate(Me%direction)       

    end subroutine ClearVariables
    
    !------------------------------------------------------------------------

    !Here we're fooling converttohdf5 in thinking there's an all-water bathym
    !GRiflet: This subroutine must be completely rewritten ...
    !It's where we create a fake bathymetry from an input "grid" 
    !file in the grd MOHID format (the template reads it from a netcdf).
    subroutine LoadIHRadarGrid

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        logical                                 :: exist
        integer                                 :: STAT_CALL


        !Begin----------------------------------------------------------------

        !Verifies if file exists
        inquire(file = Me%InputGridFile, exist = exist)
i1:     if (exist) then

            !Loads grid from file          
            call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%InputGridFile, STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'LoadIHRadarGrid - ModuleIHRadarFormat - ERR10'
            
            !Reads size of grid
            call GetHorizontalGridSize(Me%ObjHorizontalGrid, Size = Me%Size, &
                                        WorkSize = Me%Worksize, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'LoadIHRadarGrid - ModuleIHRadarFormat - ERR20'
            
            !Load the grid corners coordinates XX and YY
            call GetCornersCoordinates(Me%ObjHorizontalGrid, Me%XX_IE, Me%YY_IE, STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'LoadIHRadarGrid - ModuleIHRadarFormat - ERR30'

        else i1

            write (*,*) "The input grid file do not exist : ",trim(Me%InputGridFile)
            stop 'LoadIHRadarGrid - ModuleIHRadarFormat - ERR150'

        endif i1

    end subroutine LoadIHRadarGrid
    
    !------------------------------------------------------------------------

    
    !------------------------------------------------------------------------
    
    subroutine ConstructGridDataFile()
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL !, UnitID, i, j
        
        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing grid data...'

        !Create a new bathymetry
        allocate(Me%OpenPoints(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
        allocate(Me%Bathymetry(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
        
        call SetMatrixValue(Me%OpenPoints, Me%WorkSize, WaterPoint)
        call SetMatrixValue(Me%Bathymetry, Me%WorkSize, 100., MapMatrix = Me%OpenPoints)

        !Write new grid data file with bathymetry
        call WriteGridData  (FileName       = Me%GridFileName,                              &
                             COMENT1        = 'Grid Data created from IHRadar grid file',   &
                             COMENT2        = trim(Me%InputGridFile),                       &
                             HorizontalGridID = Me%ObjHorizontalGrid,                          &
                             FillValue      = 100.,                                         &
                             Overwrite      = .true.,                                       &
                             GridData2D_Real= Me%Bathymetry,                                &
                             STAT           = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructGridDataFile - ModuleIHRadarFormat - ERR100'

    end subroutine ConstructGridDataFile

    
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------

    subroutine OpenAndReadIHRadarFields

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

                    if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadIHRadarFields - ModuleIHRadarFormat - ERR10'

                    inquire(file = InputFile, exist = exist)
       
i1:                 if (exist) then

                        write(*,*) '... processing file ', line - FirstLine
                                                    
                        call ReadIHRadarFile (InputFile)
                        
                        call WriteIHRadar (line - FirstLine)

                    endif i1
                    
                enddo

            else BF

                stop 'OpenAndReadIHRadarFields - ModuleIHRadarFormat - ERR20'

            end if BF

            call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OpenAndReadIHRadarFields - ModuleIHRadarFormat - ERR30'

        else   IS

            stop 'OpenAndReadIHRadarFields - ModuleIHRadarFormat - ERR40'

        end if IS

    end subroutine OpenAndReadIHRadarFields


    !------------------------------------------------------------------------

    subroutine ReadIHRadarFile(InputFile)
    !GRiflet: this one must be completely rewritten.
    !
    !Local subsubroutines : 
    !   CheckName, check
    !   OutputInstants, check
    !   GetPropertyName, 
    !   GetPropertyIDNumber
    !   WriteHDF5Field

        !Arguments-------------------------------------------------------------
        character(len=PathLength)               :: InputFile        

        !Local-----------------------------------------------------------------
        integer                                 :: ObjEnterData, ClientNumber
        integer                                 :: GridVersion
        character(len=PathLength)               :: dataline
        logical                                 :: BlockFound
        integer                                 :: iflag, line, FirstLine, LastLine, lines, length, i, j
        real, dimension(:), pointer             :: lon, lat, u, v, stdu, stdv, cov, xdist, ydist, rang, bear, vel, dir
        real                                    :: aux
        integer                                 :: STAT_CALL
        real, dimension(:), pointer             :: axis_ptr

        !Begin-----------------------------------------------------------------

        !Griflet: Let's read the IH Radar ascii file
        call ConstructEnterData (ObjEnterData, InputFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadIHRadarFile - ModuleIHRadarFormat - ERR01'

        !Griflet: Check if grid version is consistent (make this GetData call its own logical function)
        call GetData(GridVersion,                                              &
                     ObjEnterData, iflag,                                           &
                     SearchType   = fromFile,                                       &
                     keyword      = '%GridVersion',                            &
                     ClientModule = 'ModuleIHRadarFormat',                          &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadIHRadarFile - ModuleIHRadarFormat - ERR02'
        
        if (GridVersion == Me%GridVersion) then

            !Griflet: Get time
            call GetData(Me%Time,                                             &
                         ObjEnterData, iflag,                                   &
                         SearchType   = FromFile,                               &
                         keyword      = '%TimeStamp',                             &
                         ClientModule = 'ModuleIHRadarFormat',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadIHRadarFile - ModuleIHRadarFormat - ERR03'
            
            !Griflet: Get TableRows or lines
            call GetData(lines,                                                 &
                         ObjEnterData, iflag,                                   &
                         SearchType   = FromFile,                               &
                         keyword      = '%TableRows',                           &
                         ClientModule = 'ModuleIHRadarFormat',                  &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadIHRadarFile - ModuleIHRadarFormat - ERR04'

            !Griflet: Look for block, loop for each line in block
            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                    &
                                       data_begin, data_end,                  &
                                       BlockFound = BlockFound,                      &
                                       FirstLine = FirstLine, LastLine = LastLine,          &
                                       STAT = STAT_CALL)
            if(STAT_CALL .eq. SUCCESS_) then
            
                if (BlockFound) then

                    !Griflet: allocate arrays and clean matrices
                    call SetMatrixValue(Me%Openpoints, Me%WorkSize, WaterPoint)
                    call SetMatrixValue(Me%velu, Me%WorkSize, FillValueReal)
                    call SetMatrixValue(Me%velv, Me%WorkSize, FillValueReal)
                    call SetMatrixValue(Me%velu_stdv, Me%WorkSize, FillValueReal)
                    call SetMatrixValue(Me%velv_stdv, Me%WorkSize, FillValueReal)
                    call SetMatrixValue(Me%covariance, Me%WorkSize, FillValueReal)
                    call SetMatrixValue(Me%xdistance, Me%WorkSize, FillValueReal)
                    call SetMatrixValue(Me%ydistance, Me%WorkSize, FillValueReal)
                    call SetMatrixValue(Me%range, Me%WorkSize, FillValueReal)
                    call SetMatrixValue(Me%bearing, Me%WorkSize, FillValueReal)
                    call SetMatrixValue(Me%velocity, Me%WorkSize, FillValueReal)
                    call SetMatrixValue(Me%direction, Me%WorkSize, FillValueReal)
                    
                    allocate(lon(1:lines))
                    allocate(lat(1:lines))
                    allocate(u(1:lines))
                    allocate(v(1:lines))
                    allocate(stdu(1:lines))
                    allocate(stdv(1:lines))
                    allocate(cov(1:lines))
                    allocate(xdist(1:lines))
                    allocate(ydist(1:lines))
                    allocate(rang(1:lines))
                    allocate(bear(1:lines))
                    allocate(vel(1:lines))
                    allocate(dir(1:lines))

                    !Griflet: skip the headlines (first two lines)
                    !and read the ascii file data
                    do line = FirstLine + 3, LastLine - 1
                    
                        !Griflet: compute corretct line index
                        i = line - Firstline -2                    

                        !Griflet: get one data line ...
                        call GetData(dataline, EnterDataID = ObjEnterData, flag = iflag, &
                                     Buffer_Line = line, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ReadIHRadarFile - ModuleIHRadarFormat - ERR10'                    
                        
                        dataline = adjustl (dataline)
                        length  = len_trim(dataline)                
                        read (dataline(1 : length),*) &
                            lon(i), lat(i), u(i), v(i), stdu(i), &
                            stdv(i), aux, cov(i), xdist(i), ydist(i), &
                            rang(i), bear(i), vel(i), dir(i)  
                        continue
                        if (STAT_CALL .EQ. SUCCESS_) then
                        end if
                                            
                    end do
                    
                    !Griflet: convert the ascii file data from arrays to an hdf5 matrix
                    !Must invocate an algorithm of point in polygon,
                    do line = 1, lines
                        i = 0
                        j = 0
                        !call LocateCellPolygons(Me%XX_IE, Me%YY_IE,                         &
                        !                    lon(line), lat(line), Me%OpenPoints,            &
                        !                    Me%WorkSize%ILB, Me%WorkSize%IUB,           &
                        !                    Me%WorkSize%JLB, Me%WorkSize%JUB,           &
                        !                    i, j)
                        axis_ptr => Me%XX_IE(1,Me%WorkSize%JLB:Me%WorkSize%JUB+1)
                        j = LocateCellIn1DAxis(axis_ptr, lon(line), Me%WorkSize%JLB, Me%WorkSize%JUB+1)
                        axis_ptr => Me%YY_IE( Me%WorkSize%ILB:Me%WorkSize%IUB+1,1)
                        i = LocateCellIn1DAxis(axis_ptr, lat(line), Me%WorkSize%ILB, Me%WorkSize%IUB+1)
                        !Griflet: set the confirmed water point as zero
                        Me%OpenPoints(i,j) = 0                   
                        Me%velu(i,j) = u(line) * 1E-2
                        Me%velv(i,j) = v(line) * 1E-2
                        Me%velu_stdv(i,j) = stdu(line) * 1E-2
                        Me%velv_stdv(i,j) = stdv(line) * 1E-2
                        Me%covariance(i,j) = cov(line)
                        Me%xdistance(i,j) = xdist(line) * 1E3
                        Me%ydistance(i,j) = ydist(line) * 1E3
                        Me%range(i,j) = rang(line) * 1E3
                        Me%bearing(i,j) = bear(line)
                        Me%velocity(i,j) = vel(line) * 1E-2
                        Me%direction(i,j) = dir(line)
                    end do
                    
                    !Griflet: switch the mask values (zero become one and one become zero)
                    do j = Me%WorkSize%JLB,Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB,Me%WorkSize%IUB
                        Me%OpenPoints(i,j) = WaterPoint - Me%OpenPoints(i,j)
                    enddo
                    enddo
                    
                    deallocate(lon)
                    deallocate(lat)
                    deallocate(u)
                    deallocate(v)
                    deallocate(stdu)
                    deallocate(stdv)
                    deallocate(cov)
                    deallocate(xdist)
                    deallocate(ydist)
                    deallocate(rang)
                    deallocate(bear)
                    deallocate(vel)
                    deallocate(dir)
                    
                end if
                
            end if        
        
        else
        
            write(*,*) 'Inconsistent IH Seasonde Radar grid version with declared IH_GRID_VERSION in ConvertToHdf5Action.dat'
            stop 'ReadIHRadarFile - ModuleIHRadarFormat - ERR230'
        
        end if            

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadIHRadarFile - ModuleIHRadarFormat - ERR240'

    end subroutine ReadIHRadarFile

    !--------------------------------------------------------------------------

    integer function LocateCellIn1DAxis(axis, pos, ILB, IUB)

        !Arguments -------------------------------------------------
        real, dimension(:), pointer                     :: axis
        real                                            :: pos
        integer                                         :: ILB, IUB        

        !Local -----------------------------------------------------
        integer                                         :: ICenter, IL, IU
        logical                                         :: cellfound        

        !Begin -----------------------------------------------------        

        !Initialize algorithm
        cellfound = .false.        
        IL = ILB
        IU = IUB        

        !Iterative algorithm
        do while (.not.cellfound)
                
            ICenter = (IL + IU)/2     
                                           
            if ( pos > axis(ICenter) ) then
                !Cell is in [ICenter IU]
                IL = ICenter            
            else if ( pos <= axis(ICenter) ) then
                !Cell is in [IL ICenter]
                IU = ICenter
            else
                 stop 'LocateCellIn1DAxis - ModuleIHRadarFormat - ERR10'
            end if           
             
            if (IU-IL == 1) cellfound = .true. 
                   
        end do

        LocateCellIn1DAxis = IL

    end function LocateCellIn1DAxis

    !--------------------------------------------------------------------------
      
    subroutine Open_Output_File

        !Construct file
        !Write Bathymetry
        !Write Horizontal Grid
        !Write WaterPoints2D

        if (Me%OutputHDF5) then
            call Open_HDF5_Output_File
        endif
        
        if (Me%OutputNetcdf) then
            !Griflet: todo
        endif
                
    end subroutine

    !------------------------------------------------------------------------

    subroutine WriteIHRadar(iOut)
    
        !Arguments -------------------------------------------------
        integer                                         :: iOut
        
        !Begin -----------------------------------------------------

        !Write OpenPoints
        !Write Fields
        !Write Time

        if (Me%OutputHDF5) then        
            call Write_HDF5_Fields(iOut)
        endif
                
        if (Me%OutputNetcdf) then        
            !Griflet: todo
        endif
        
    end subroutine
  
    !--------------------------------------------------------------------------
    
    subroutine Close_Output_File

        !Kill file

        if (Me%OutputHDF5) then        
            call Close_HDF5_Output_File
        endif
                
        if (Me%OutputNetcdf) then        
            !Griflet: todo
        endif

    end subroutine

    !--------------------------------------------------------------------------

    subroutine Open_HDF5_Output_File

        !Local-----------------------------------------------------------------
        integer                                     :: HDF5_CREATE
        integer                                     :: HDF5_IO_CODE
        integer                                     :: STAT_CALL

        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputHDF5FileName, HDF5_IO_CODE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_Output_File - ModuleMecatorFormat - ERR01'

        call WriteHDF5GridData

    end subroutine

    !------------------------------------------------------------------------
    
    subroutine WriteHDF5GridData

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDF5GridData - ModuleMecatorFormat - ERR02'
        
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDF5GridData - ModuleMecatorFormat - ERR03'            

        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDF5GridData - ModuleMecatorFormat - ERR04'            
           
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDF5GridData - ModuleMecatorFormat - ERR07'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",    &
                              Array2D = Me%OpenPoints,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDF5GridData - ModuleMecatorFormat - ERR08'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDF5GridData - ModuleMecatorFormat - ERR09'

    end subroutine WriteHDF5GridData

    !--------------------------------------------------------------------------

    subroutine Write_HDF5_Fields(iOut)

        !Arguments -------------------------------------------------
        integer                                         :: iOut

        !Locals ----------------------------------------------------
        character(StringLength)                         :: MohidName
        character(StringLength)                         :: units
        integer                                         :: STAT_CALL

        !Begin------------------------------------------------------

        !Write OpenPoints
        call WriteHDF5Openpoints(iOut)

        !Write fields
        MohidName = "velocity U"
        units = "m/s"
        call WriteHDF5Field(MohidName, Me%velu, units, iOut)
        MohidName = "velocity V"
        units = "m/s"
        call WriteHDF5Field(MohidName, Me%velv, units, iOut)
        MohidName = "velU std dev"
        units = "m/s"
        call WriteHDF5Field(MohidName, Me%velu_stdv, units, iOut)
        MohidName = "velV std dev"
        units = "m/s"
        call WriteHDF5Field(MohidName, Me%velv_stdv, units, iOut)
        MohidName = "covariance"
        units = "%"
        call WriteHDF5Field(MohidName, Me%covariance, units, iOut)
        MohidName = "x distance"
        units = "m"
        call WriteHDF5Field(MohidName, Me%xdistance, units, iOut)
        MohidName = "y distance"
        units = "m"
        call WriteHDF5Field(MohidName, Me%ydistance, units, iOut)
        MohidName = "range"
        units = "m"
        call WriteHDF5Field(MohidName, Me%range, units, iOut)
        MohidName = "bearing"
        units = "º"
        call WriteHDF5Field(MohidName, Me%bearing, units, iOut)
        MohidName = "velocity modulus"
        units = "m/s"
        call WriteHDF5Field(MohidName, Me%velocity, units, iOut)
        MohidName = "velocity direction"
        units = "º"
        call WriteHDF5Field(MohidName, Me%direction, units, iOut)

        !Write time
        call WriteHDF5Time(iOut)

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Fields - ModuleIHRadarFormat - ERR10'

    end subroutine Write_HDF5_Fields

    !------------------------------------------------------------------------

    subroutine WriteHDF5Openpoints(iOut)

        !Arguments-------------------------------------------------------------
        integer                                         :: iOut

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: WorkILB, WorkJLB
        integer                                         :: WorkIUB, WorkJUB

        !Begin-----------------------------------------------------------------
        
        !Bounds
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        !Begin-----------------------------------------------------------------
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                     &
                             WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Openpoints - ModuleIHRadarFormat - ERR10'

        call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints",                  &
                                 "OpenPoints","-", Array2D = Me%OpenPoints,      &
                                 OutputNumber = iOut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Openpoints - ModuleIHRadarFormat - ERR20'
        
    end subroutine WriteHDF5Openpoints
 
    !------------------------------------------------------------------------

    subroutine WriteHDF5Field(MohidName, Field, PropUnits, iOut)

        !Arguments-------------------------------------------------------------
        character(Len=StringLength)                     :: MohidName
        real, dimension(:,:  ), pointer                 :: Field
        character(Len=StringLength)                     :: PropUnits
        integer                                         :: iOut

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: WorkILB, WorkJLB
        integer                                         :: WorkIUB, WorkJUB

        !Begin-----------------------------------------------------------------
        
        !Bounds
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                     &
                             WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleIHRadarFormat - ERR10'

        call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(MohidName),                  &
                                 trim(MohidName),trim(PropUnits), Array2D = Field,      &
                                 OutputNumber = iOut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleIHRadarFormat - ERR20'

    end subroutine WriteHDF5Field 

    !------------------------------------------------------------------------

    subroutine WriteHDF5Time(iOut)

        !Arguments-------------------------------------------------------------
        integer                                         :: iOut

        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------
        !Writes current time
        call ExtractDate   (Me%Time, AuxTime(1), AuxTime(2), AuxTime(3),              &
                                           AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime
        
        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Time - ModuleIHRadarFormat - ERR10'

        call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",        &
                               Array1D = TimePtr, OutputNumber = iOut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Time - ModuleIHRadarFormat - ERR20'
        
    end subroutine WriteHDF5Time
 
    !------------------------------------------------------------------------

    subroutine Close_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Close_HDF5_OutPut_File - ModuleIHRadarFormat - ERR60'

    end subroutine Close_HDF5_OutPut_File

    !--------------------------------------------------------------------------
    
    subroutine Open_Netcdf_Output_File

        !Local-----------------------------------------------------------------
        integer                                     :: NCDF_CREATE, STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        call GetNCDFFileAccess(NCDF_CREATE = NCDF_CREATE)
        
        call ConstructNETCDF(Me%ObjNETCDF, Me%OutputNetcdfFileName, NCDF_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OpenNCDFFile - ModuleIHRadarFormat - ERR01'

        call NETCDFWriteHeader(NCDFID         = Me%ObjNETCDF,                 &
                               Title          = Me%OutputNetcdfAttr%Title,           &
                               Convention     = Me%OutputNetcdfAttr%Convention,      &
                               Version        = Me%OutputNetcdfAttr%Version,         &
                               History        = Me%OutputNetcdfAttr%History,         &
                               iDate          = Me%OutputNetcdfAttr%iDate,           &
                               Source         = Me%OutputNetcdfAttr%Source,          &
                               Institution    = Me%OutputNetcdfAttr%Institution,     &
                               References     = Me%OutputNetcdfAttr%References,      &
                               STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OpenNCDFFile - ModuleIHRadarFormat - ERR02'
        
        write(*,*)
        write(*,*)'Opened ncdf file                : ', trim(Me%OutputNetcdfFileName)

        call NETCDFSetDimensions(Me%ObjNETCDF, Me%Size%IUB, Me%Size%JUB, 0, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadSetDimensions - Convert2netcdf - ERR02'

    end subroutine Open_Netcdf_Output_File
    
    !--------------------------------------------------------------------------

    subroutine KillIHRadarFormat
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%XX_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillIHRadarFormat - ModuleIHRadarFormat - ERR10'
        
        call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillIHRadarFormat - ModuleIHRadarFormat - ERR20'
        
        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillIHRadarFormat - ModuleIHRadarFormat - ERR50'

!        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'KillIHRadarFormat - ModuleIHRadarFormat - ERR60'

        deallocate(Me%Bathymetry)
        deallocate(Me%OpenPoints)

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillIHRadarFormat - ModuleIHRadarFormat - ERR70'

        deallocate(Me)
        nullify   (Me)

    end subroutine KillIHRadarFormat

    !--------------------------------------------------------------------------

end module ModuleIHRadarFormat