!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Convert2netcdf
! PROGRAM       : Convert2netcdf
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : October 2006
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Program to convert MOHID HDF5 files to NETCDF format
!
!------------------------------------------------------------------------------

!Input data file must be named: "Convert2netcdf.dat"

! HDF_FILE             : [char]        [-]                      !Path to the HDF5 file to be converted
! HDF_SIZE_GROUP       : [char]     ["/Grid"]                   !HDF group containing the HDF dataset from where 
                                                                !the dimensions will be read
! HDF_SIZE_DATASET     : [char]   ["WaterPoints3D"]             !Name of the HDF dataset from where the dimensions
                                                                !will be read
! HDF_TIME_VAR         : [char]      ["Time"]                   !Name of the HDF group containing the 
                                                                !time values
! HDF_READ_LATLON      : int           [1]                      !Type of projection (0 - metric, 1 - lat/lon)
! HDF_READ_SIGMA       : int           [0]                      !Type of vertical coordinate (0 - cartesian, 1 - sigma)
!
! IMPOSE_MASK          : bool          [0]                      !In case we want to specify the mask variable.
! HDF_MASK             : [char]   ["WaterPoints3D"]             !If we chose to specify a mask variable, then we define here
                                                                !the variable name
! HDF_MASK_IS_3D       : bool          [1]                      !Is the mask variable 3D? yes or no.
! RESULTS_ARE_2D       : bool          [0]                      !Sometimes the mask variable is 3D, but the results
                                                                !are 2D (ex: surface results). Use this keyword to enforce
                                                                !this condition
!
! NETCDF_FILE          : [char]        [-]                      !netcdf file to be created
! NETCDF_TITLE         : [char]        [-]                      !netcdf file title 
! NETCDF_CONVENTION    : [char]     ["CF-1.0"]                  !netcdf naming convention   
! NETCDF_VERSION       : [char]      [3.6.1]                    !netcdf library version
! NETCDF_HISTORY       : [char]        [-]                      !netcdf file history
! NETCDF_SOURCE        : [char]        [-]                      !netcdf file source 
! NETCDF_INSTITUTION   : [char] ["Instituto Superior T�cnico"]  !netcdf institution
! NETCDF_REFERENCES    : [char]   [http://www.mohid.com/]       !netcdf references
! NETCDF_DATE          : int           [-]                      !current year

! DEPTH_OFFSET         : real          [0.]                     !add_offset attribute for depth
! REFERENCE_TIME       : date    [2004 1 1 0 0 0]               !initial date for time dimension
! CONVERT_EVERYTHING   : bool          [1]                      !convert all info in HDF file or not

!<begin_groups>
!/Results/salinity     : e.g. converts only the salinity results
!<end_groups>

! MULTIPLY_FACTOR      : real          [1.]                     !multiply by fACTOR 
! ADD_FACTOR


program Convert2netcdf

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleHDF5
    use ModuleNETCDF
    use ModuleFunctions
    use HDF5

    implicit none

    type T_HDFFile
        integer                                             :: FileID       = 0
        character(len=PathLength)                           :: Name         = null_str
        integer                                             :: ObjHDF5      = 0
        integer                                             :: nInstants    = FillValueInt
        real(8), dimension(:), pointer                      :: Times        => null()
        type(T_Time)                                        :: InitialDate
        type(T_Size3D)                                      :: Size
        character(len=StringLength)                         :: SizeGroup, SizeDataSet
        character(len=StringLength)                         :: HdfMask          = null_str
        character(len=StringLength)                         :: TimeVar          = null_str
        character(len=StringLength)                         :: TimeGroupName    = null_str
        character(len=StringLength)                         :: VertVar          = null_str
        logical                                             :: ReadLatLon       = .true.
        logical                                             :: Sigma            = .false.
        logical                                             :: HdfMaskIs3D      = .true.
        logical                                             :: ImposeMask       = .false.
        logical                                             :: ResultsAre2D     = .false.
        logical                                             :: OutputIs2D       = .false. 
    end type T_HDFFile                                      
                                                            
    type T_NCDFFile                                         
        character(len=PathLength)                           :: Name                 = null_str
        integer                                             :: ObjNETCDF            = 0
        character(len=StringLength)                         :: Title                = null_str
        character(len=StringLength)                         :: Convention           = null_str
        character(len=StringLength)                         :: Version              = null_str
        character(len=StringLength)                         :: History              = null_str  
        character(len=StringLength)                         :: Source               = null_str
        character(len=StringLength)                         :: Institution          = null_str
        character(len=StringLength)                         :: References           = null_str
        integer                                             :: iDate                = FillValueInt
        real                                                :: geospatial_lat_min   = FillValueReal
        real                                                :: geospatial_lat_max   = FillValueReal
        real                                                :: geospatial_lon_min   = FillValueReal
        real                                                :: geospatial_lon_max   = FillValueReal
        character(len=StringLength)                         :: CoordSysBuilder      = null_str
        character(len=StringLength)                         :: contact              = null_str
        character(len=StringLength)                         :: field_type           = null_str
        character(len=StringLength)                         :: bulletin_date        = null_str
        character(len=StringLength)                         :: bulletin_type        = null_str
        character(len=StringLength)                         :: comment              = null_str
        character(len=StringLength)                         :: MetadataAtt          = null_str        
        character(len=LinkLength  )                         :: MetadataLink         = null_str
    end type T_NCDFFile

    type T_Conv2netcdf

        character(len=PathLength)                           :: DataFile         = null_str
        character(len=PathLength)                           :: netcdfFile
                                                            
        type(T_Time)                                        :: InitialSystemTime
        type(T_Time)                                        :: FinalSystemTime
        real                                                :: TotalCPUTime
        real                                                :: ElapsedSeconds
        integer, dimension(8)                               :: F95Time
        logical                                             :: MohidStandardInOutUnits = .false. 
        logical                                             :: OdysseaProject          = .false. 
                                                            
        integer                                             :: ObjEnterData     = 0
                                                            
        real,    dimension(:,:  ), pointer                  :: Float2DIn           => null()
        real,    dimension(:,:,:), pointer                  :: Float3DIn           => null()
        integer, dimension(:,:  ), pointer                  :: Int2DIn             => null()
        integer, dimension(:,:,:), pointer                  :: Int3DIn             => null()
                                                            
        real,    dimension(:,:  ), pointer                  :: Float2DOut          => null()
        real,    dimension(:,:,:), pointer                  :: Float3DOut          => null()
        integer, dimension(:,:  ), pointer                  :: Int2DOut            => null()
        integer, dimension(:,:,:), pointer                  :: Int3DOut            => null()
        
        logical                                             :: IsMapping           = .false. 


        type(T_HDFFile)                                     :: HDFFile
        type(T_NCDFFile)                                    :: NCDF_File

        real                                                :: DepthAddOffSet   = 0.
        type(T_Time)                                        :: ReferenceTime
        integer                                             :: DecimalPlaces    = FillValueInt
        logical                                             :: StagGoogleOut    = .true. 

        logical                                             :: ConvertEverything=ON
        integer                                             :: nGroupsToConvert
        character(len=StringLength), dimension(:), pointer  :: GroupsToConvert => null()

        real,   dimension(:,:,:), pointer                   :: Depth3DIN        => null()   
        real,   dimension(:), pointer                       :: DepthVector      => null()
        integer                                             :: DepthLayers      =  FillValueInt
        logical                                             :: DepthLayersON    = .false. 
        
        real                                                :: Add_Factor       = FillValueReal
        real                                                :: Multiply_Factor  = FillValueReal        
        
        logical                                             :: SimpleGrid       = .false. 
        
        real                                                :: MissingValue     = FillValueReal

    end type T_Conv2netcdf

    type(T_Conv2netcdf)                     :: Me

    call ConstructConvert2netcdf
    call ModifyConvert2netcdf
    call KillConvert2netcdf

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructConvert2netcdf

        call StartUpMohid("Convert2netcdf")

        call StartCPUTime

        call ReadKeywords

    end subroutine ConstructConvert2netcdf
    
    !--------------------------------------------------------------------------

    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag
        logical                                     :: exist, Exist2
        real, dimension(:), allocatable             :: Aux1D
        
        !Begin-----------------------------------------------------------------
        
        !Read input file name from nomfich file
        call ReadFileName('IN_MODEL', Me%DataFile, "Convert2netcdf", STAT = STAT_CALL)
        
        if     (STAT_CALL == FILE_NOT_FOUND_ERR_) then
            Me%DataFile = "Convert2netcdf.dat"
        elseif (STAT_CALL /= SUCCESS_              ) then
            stop 'ReadKeywords - ModuleValida4D - ERR10'
        endif                        
                       

        call ConstructEnterData (Me%ObjEnterData, Me%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR20'

        call GetData(Me%NCDF_File%Name,                                                 &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_FILE',                                      &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR30'

        call GetData(Me%NCDF_File%Title,                                                &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_TITLE',                                     &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR40'

        call GetData(Me%NCDF_File%Convention,                                           &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_CONVENTION',                                &
                     Default      = 'CF-1.0',                                           &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR50'

        call GetData(Me%NCDF_File%Version,                                              &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_VERSION',                                   &
                     Default      = '3.6.1',                                            &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR60'

        call GetData(Me%NCDF_File%History,                                              &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_HISTORY',                                   &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR70'
        
        call GetData(Me%NCDF_File%Source,                                               &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_SOURCE',                                    &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR80'

        
        call GetData(Me%NCDF_File%Institution,                                          &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_INSTITUTION',                               &
                     Default      = 'Instituto Superior Tecnico',                       &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR90'
        
        call GetData(Me%NCDF_File%References,                                           &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_REFERENCES',                                &
                     Default      = 'http://www.mohid.com',                             &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR100'

        call GetData(Me%NCDF_File%iDate,                                                &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_DATE',                                      &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR110'
        
        
        call GetData(Me%NCDF_File%geospatial_lat_min,                                   &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_LAT_MIN',                                   &
                     Default      = -90.,                                               &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR120'        
        
        call GetData(Me%NCDF_File%geospatial_lat_max,                                   &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_LAT_MAX',                                   &
                     Default      = +90.,                                               &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR130'          
  
        call GetData(Me%NCDF_File%geospatial_lon_min,                                   &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_LON_MIN',                                   &
                     Default      = -180.,                                              &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR140'           

        call GetData(Me%NCDF_File%geospatial_lon_max,                                   &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_LON_MAX',                                   &
                     Default      = +180.,                                              &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR150'         
        
        call GetData(Me%NCDF_File%CoordSysBuilder,                                      &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_COORD_SYSTEM',                              &
                     Default      = "ucar.nc2.dataset.conv.CF1Convention",              &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR160'
        
        call GetData(Me%NCDF_File%contact,                                              &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_CONTACT',                                   &
                     Default      = "general@mohid.com",                                &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR170'        

        call GetData(Me%NCDF_File%field_type,                                           &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_FIELD_TYPE',                                &
                     Default      = "mean",                                             &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR180'        
                
        call GetData(Me%NCDF_File%bulletin_date,                                        &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_BULLETIN_DATE',                             &
                     Default      = "2018-01-01 00:00:00",                              &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR190'            

        call GetData(Me%NCDF_File%bulletin_type,                                        &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_BULLETIN_TYPE',                             &
                     Default      = "operational",                                      &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR200'
       
        call GetData(Me%NCDF_File%comment,                                              &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_COMMENT',                                   &
                     Default      = "MOHID product",                                    &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR210'       
        
        call GetData(Me%NCDF_File%MetadataAtt,                                          &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NETCDF_METADATA_ATTRIBUTE',                        &
                     Default      = trim(null_str),                                     &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR220'           
        
        if (iflag /= 0) then
            
            call GetData(Me%NCDF_File%MetadataLink,                                     &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromFile,                                       &
                         keyword      = 'NETCDF_METADATA_LINK',                         &
                         Default      = trim(null_str),                                 &
                         ClientModule = 'Convert2netcdf',                               &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR230'           
        
        endif

        call GetData(Me%HDFFile%Name,                                                   &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF_FILE',                                         &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR240'

        !Verifies if file exists
        inquire(FILE = trim(Me%HDFFile%Name), EXIST = exist)
        if (exist) then
            call OpenHDF5File
        else
            write(*,*)'HDF5 file does not exist'
            stop 'ReadKeywords - Convert2netcdf - ERR250'
        endif


        call GetData(Me%HDFFile%TimeVar,                                                &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF_TIME_VAR',                                     &
                     Default      = 'Time',                                             &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR260'
        
        call GetData(Me%HDFFile%TimeGroupName,                                          &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF_TIME_GROUP_NAME',                              &
                     Default      = trim(Me%HDFFile%TimeVar),                           &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR265'        

        call GetData(Me%HDFFile%VertVar,                                                &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF_VERT_VAR',                                     &
                     Default      = 'VerticalZ/Vertical',                               &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR270'


        call GetData(Me%HDFFile%SizeGroup,                                              &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF_SIZE_GROUP',                                   &
                     Default      = '/Grid',                                            &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR280'

        call GetData(Me%HDFFile%ImposeMask,                                             &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'IMPOSE_MASK',                                      &
                     Default      = .false.,                                            &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR290'

        call GetData(Me%HDFFile%HdfMask,                                                &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF_MASK',                                         &
                     Default      = 'WaterPoints3D',                                    &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR300'
        
        call GetHDF5DataSetExist (HDF5ID        = Me%HDFFile%ObjHDF5,           &
                                  DataSetName   = trim(Me%HDFFile%SizeGroup)//  &
                                                  "/"//trim(Me%HDFFile%HdfMask),&
                                  Exist         = Exist,                        & 
                                  STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR310'
        
        if (.not.Exist) then
        
            call GetHDF5DataSetExist (HDF5ID        = Me%HDFFile%ObjHDF5,           &
                                      DataSetName   = trim(Me%HDFFile%SizeGroup)//  &
                                                      "/"//"WaterPoints",           &
                                      Exist         = Exist2,                       & 
                                      STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR320'
            
            if (Exist2) then
                Me%HDFFile%HdfMask = "WaterPoints"
            else
                write(*,*) 'Name define in keyword HDF_MASK not valid' 
                stop 'ReadMask - Convert2netcdf - ERR330'
            endif

        endif
        

        call GetData(Me%HDFFile%HdfMaskIs3D,                                            &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF_MASK_IS_3D',                                   &
                     Default      = .true.,                                             &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR340'

        call GetData(Me%HDFFile%ResultsAre2D,                                           &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESULTS_ARE_2D',                                   &
                     Default      = .false.,                                            &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR350'

        call GetData(Me%HDFFile%OutputIs2D,                                             &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'OUTPUT_IS_2D',                                     &
                     Default      = .false.,                                            &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR360'


        call GetData(Me%HDFFile%SizeDataSet,                                            &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF_SIZE_DATASET',                                 &
                     Default      = trim(Me%HDFFile%HdfMask),                           &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR370'
        
        call GetData(Me%HDFFile%ReadLatLon,                                             &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF_READ_LATLON',                                  &
                     Default      = .true.,                                             &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR380'

        call GetData(Me%HDFFile%Sigma,                                                  &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF_READ_SIGMA',                                   &
                     Default      = .false.,                                            &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR390'

        call GetData(Me%ConvertEverything,                                              &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'CONVERT_EVERYTHING',                               &
                     Default      = .true.,                                             &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR400'


        call GetData(Me%DepthAddOffSet,                                                 &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DEPTH_OFFSET',                                     &
                     Default      = 0.,                                                 &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR410'


        call GetData(Me%ReferenceTime,                                                  &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'REFERENCE_TIME',                                   &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR420'

        if(iflag == 0)then
            call SetDate(Me%ReferenceTime, 2004, 1, 1, 0, 0, 0)
        end if

        call GetData(Me%DecimalPlaces,                                                  &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DECIMAL_PLACES',                                   &
                     ClientModule = 'Convert2netcdf',                                   &
                     Default      = FillValueInt,                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR430'
        
        if (iflag /= 0) then
            if (Me%DecimalPlaces < 0) then
                stop 'ReadKeywords - Convert2netcdf - ERR440'
            endif
            
            if (Me%DecimalPlaces > 16) then
                stop 'ReadKeywords - Convert2netcdf - ERR450'
            endif            
            
        endif
        
       

        call GetData(Me%StagGoogleOut ,                                                 &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'STAG_GOOGLE_OUT',                                  &
                     ClientModule = 'Convert2netcdf',                                   &
                     Default      = .true.,                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR460'
        
        iflag = 0
        
        allocate(Aux1D(1:100))
        
        call GetData(Aux1D,                                                             &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DEPTH_LAYERS',                                     &
                     ClientModule = 'Convert2netcdf',                                   &
                     STAT         = STAT_CALL)
        
        if (iflag > 0) then
        
            allocate(Me%DepthVector(1:iflag))
            
            Me%DepthVector(1:iflag) = Aux1D(1:iflag)
            
            deallocate(Aux1D)
            
!            call GetData(Me%DepthVector,                                                &
!                         Me%ObjEnterData,iflag,                                         &
!                         SearchType   = FromFile,                                       &
!                         keyword      = 'DEPTH_LAYERS',                                 &
!                         ClientModule = 'Convert2netcdf',                               &
!                         STAT         = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR350'              
            
            Me%DepthLayers          = iflag    
            Me%DepthLayersON        = .true.            
            
        else
            Me%DepthLayersON        = .false.
        endif            
        
        call GetData(Me%Add_Factor,                                                     &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'ADD_FACTOR',                                       &
                     ClientModule = 'Convert2netcdf',                                   &
                     Default      = 0.,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR470'
        
        call GetData(Me%Multiply_Factor,                                                &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MULTIPLY_FACTOR',                                  &
                     ClientModule = 'Convert2netcdf',                                   &
                     Default      = 1.,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR480'
        
        call GetData(Me%SimpleGrid,                                                     &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'SIMPLE_GRID',                                      &
                     ClientModule = 'Convert2netcdf',                                   &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR490'
        
     
        if(.not. Me%ConvertEverything)then

             call ReadVGroupsToConvert

        end if
        
        call GetData(Me%MissingValue,                                                   &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MISSING_VALUE',                                    &
                     ClientModule = 'Convert2netcdf',                                   &
                     Default      = FillValueReal,                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR500'
        
        call GetData(Me%MohidStandardInOutUnits,                                        &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MOHID_STANDARD_IN_OUT_UNITS',                      &
                     ClientModule = 'Convert2netcdf',                                   &
                     Default      = .true.,                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR510'        
                    

        call GetData(Me%OdysseaProject,                                                 &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'ODYSSEA_PROJECT',                                  &
                     ClientModule = 'Convert2netcdf',                                   &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR520'        
        

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - Convert2netcdf - ERR990'

    end subroutine ReadKeywords

    !--------------------------------------------------------------------------

    subroutine ReadVGroupsToConvert

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ClientNumber, iflag, CurrentLineNumber
        logical                                     :: BlockFound
        integer                                     :: StartLine, EndLine, Count
        character(len=StringLength)                 :: PropertyName

        !Begin-----------------------------------------------------------------

        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                      &
                                    '<begin_groups>', '<end_groups>', BlockFound,       &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadVGroupsToConvert - ConvertToNETCDF - ERR01'

        if(.not. BlockFound)then
            
            write(*,*)"Could not find <begin_groups>...<end_groups> block"
            stop 'ReadVGroupsToConvert - ConvertToNETCDF - ERR10'
        
        else

            call GetBlockSize(Me%ObjEnterData,                      &
                              ClientNumber,                         &
                              StartLine,                            &
                              EndLine,                              &
                              FromBlock,                            &
                              STAT = STAT_CALL)

            Me%nGroupsToConvert = EndLine - StartLine - 1

            allocate(Me%GroupsToConvert(1:Me%nGroupsToConvert))

            Count = 1

            do CurrentLineNumber = StartLine + 1 , EndLine - 1

                call GetData(PropertyName,                          &
                             Me%ObjEnterData,                       &
                             iflag,                                 &
                             SearchType  = FromBlock_,              &
                             Buffer_Line = CurrentLineNumber,       &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadVGroupsToConvert - ConvertToNETCDF - ERR20'

                Me%GroupsToConvert(Count) = trim(PropertyName)

                Count = Count + 1

            end do

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModuleInterpolateGrids - ERR03'


        end if

    end subroutine ReadVGroupsToConvert

    !--------------------------------------------------------------------------

    subroutine ModifyConvert2netcdf
        
        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL
        !Begin-----------------------------------------------------------------
    
        !call OpenHDF5File
    

        call GetHDF5FileID (Me%HDFFile%ObjHDF5, Me%HDFFile%FileID, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyConvert2netcdf - Convert2netcdf - ERR10'

        call ReadSetDimensions    
    
        call ReadLatLonWindow    
        
        call OpenNCDFFile
        
        call ReadWriteTime

        call ReadWriteGrid

        call ReadWriteData

        call CloseFiles

    end subroutine ModifyConvert2netcdf
    
    !--------------------------------------------------------------------------
   
    subroutine OpenHDF5File

        !Local-----------------------------------------------------------------
        integer                             :: HDF5_READ, STAT_CALL

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        call ConstructHDF5 (Me%HDFFile%ObjHDF5, Me%HDFFile%Name, Access = HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5File - Convert2netcdf - ERR01'

        write(*,*)'Opened hdf5 file                : ', trim(Me%HDFFile%Name)


    end subroutine OpenHDF5File

    !--------------------------------------------------------------------------
    
    subroutine CloseFiles
        
        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call KillHDF5 (Me%HDFFile%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CloseFiles - Convert2netcdf - ERR01'
       
        call KillNETCDF(Me%NCDF_File%ObjNETCDF, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CloseFiles - Convert2netcdf - ERR02'
        
        write(*,*)

    end subroutine CloseFiles
    
    !--------------------------------------------------------------------------
       
    subroutine OpenNCDFFile

        !Local-----------------------------------------------------------------
        integer                                     :: NCDF_CREATE, STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call GetNCDFFileAccess(NCDF_CREATE = NCDF_CREATE)
        
        call ConstructNETCDF(Me%NCDF_File%ObjNETCDF, Me%NCDF_File%Name, NCDF_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OpenNCDFFile - Convert2netcdf - ERR01'

        call NETCDFWriteHeader (NCDFID              = Me%NCDF_File%ObjNETCDF,           &
                                Title               = Me%NCDF_File%Title,               &
                                Convention          = Me%NCDF_File%Convention,          &
                                Version             = Me%NCDF_File%Version,             &
                                History             = Me%NCDF_File%History,             &
                                iDate               = Me%NCDF_File%iDate,               &
                                Source              = Me%NCDF_File%Source,              &
                                Institution         = Me%NCDF_File%Institution,         &
                                References          = Me%NCDF_File%References,          &
                                geospatial_lat_min  = Me%NCDF_File%geospatial_lat_min,  & 
                                geospatial_lat_max  = Me%NCDF_File%geospatial_lat_max,  &
                                geospatial_lon_min  = Me%NCDF_File%geospatial_lon_min,  & 
                                geospatial_lon_max  = Me%NCDF_File%geospatial_lon_max,  &
                                CoordSysBuilder     = Me%NCDF_File%CoordSysBuilder,     &
                                contact             = Me%NCDF_File%contact,             & 
                                field_type          = Me%NCDF_File%field_type,          &     
                                bulletin_date       = Me%NCDF_File%bulletin_date,       &
                                bulletin_type       = Me%NCDF_File%bulletin_type,       & 
                                comment             = Me%NCDF_File%comment,             &   
                                MetadataAtt         = Me%NCDF_File%MetadataAtt,         &   
                                MetadataLink        = Me%NCDF_File%MetadataLink,        &   
                                STAT                = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OpenNCDFFile - Convert2netcdf - ERR02'
        
        if (Me%DepthLayersON) then
            call NETCDFSetDimensions(Me%NCDF_File%ObjNETCDF, int(Me%HDFFile%Size%IUB,4),int(Me%HDFFile%Size%JUB,4), &
                                     int(Me%DepthLayers,4), SimpleGrid = Me%SimpleGrid, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'OpenNCDFFile - Convert2netcdf - ERR20'
        else
            call NETCDFSetDimensions(Me%NCDF_File%ObjNETCDF, int(Me%HDFFile%Size%IUB,4), int(Me%HDFFile%Size%JUB,4),&
                                     int(Me%HDFFile%Size%KUB,4), SimpleGrid = Me%SimpleGrid, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'OpenNCDFFile - Convert2netcdf - ERR30'
        endif        
        
        
        write(*,*)
        write(*,*)'Opened ncdf file                : ', trim(Me%NCDF_File%Name)

    end subroutine OpenNCDFFile
    
    !--------------------------------------------------------------------------

    subroutine ReadWriteTime

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: CurrentInstant
        character(len=19)                                   :: CharInitialDate  

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)"Reading time..."
        write(*,*)

        !Gets number of time instants
        call GetHDF5GroupNumberOfItems(Me%HDFFile%ObjHDF5, '/'//trim(Me%HDFFile%TimeGroupName), &
                                       Me%HDFFile%nInstants, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadWriteTime - Convert2netcdf - ERR01'

        write(*,*)'Number of instants in hdf5 file : ', Me%HDFFile%nInstants

        Me%HDFFile%InitialDate = HDF5TimeInstant(1)

        allocate(Me%HDFFile%Times(1:Me%HDFFile%nInstants))

        do CurrentInstant = 1, Me%HDFFile%nInstants

            Me%HDFFile%Times(CurrentInstant) = HDF5TimeInstant(CurrentInstant) - Me%ReferenceTime

        end do

        CharInitialDate = TimeToStringV2(Me%ReferenceTime)

        call NETCDFWriteTime(NCDFID         = Me%NCDF_File%ObjNETCDF,           &
                             InitialDate    = CharInitialDate,                  &
                             nInstants      = Me%HDFFile%nInstants,             &
                             Times          = Me%HDFFile%Times,                 &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadWriteTime - Convert2netcdf - ERR02'

    end subroutine ReadWriteTime

    !--------------------------------------------------------------------------

    subroutine ReadWriteGrid

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
       
        write(*,*)
        write(*,*)"Reading grid..."
        write(*,*)

        call ReadWriteLatLon

        if (Me%HDFFile%Size%KUB .gt. 0 .and. .not. Me%HDFFile%OutputIs2D) then
            if ( .not. Me%HDFFile%ResultsAre2D ) then
                if (Me%DepthLayersON) then
                    call WriteDepthLayers
                else
                    call ReadWriteVertical
                endif                                        
            else
                call WriteVerticalNullDepth
            endif
        endif

        if (.not.Me%OdysseaProject) then

            call ReadWriteBathymetry

            call ReadMask
            
        endif                        

    end subroutine ReadWriteGrid

    !--------------------------------------------------------------------------

    subroutine ReadSetDimensions

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: gr_id, dset_id, class_id
        !integer(HID_T)                              :: ssize
        integer(SIZE_T)                             :: ssize        
        integer(HID_T)                              :: space_id, datatype_id
        integer(HID_T)                              :: rank
        integer(HSIZE_T), dimension(7)              :: dims, maxdims
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        dims = 0        
        write(*,*)"Reading sizes..."

        call h5gopen_f(Me%HDFFile%FileID, Me%HDFFile%SizeGroup, gr_id, STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadSetDimensions - Convert2netcdf - ERR10'

        !Opens data set
        call h5dopen_f(gr_id, trim(adjustl(Me%HDFFile%SizeDataSet)), dset_id, STAT_CALL)
        call h5dget_space_f                 (dset_id, space_id, STAT_CALL)
        call h5sget_simple_extent_ndims_f   (space_id, rank, STAT_CALL)
        call h5dget_type_f                  (dset_id, datatype_id,   STAT_CALL)
        call h5tget_size_f                  (datatype_id, ssize,      STAT_CALL)
        call h5tget_class_f                 (datatype_id, class_id,  STAT_CALL) 
        call h5tclose_f                     (datatype_id, STAT_CALL) 
        call h5sget_simple_extent_dims_f    (space_id, dims, maxdims, STAT_CALL) 
        call h5sclose_f                     (space_id, STAT_CALL)
        call h5dclose_f                     (dset_id, STAT_CALL)
        call h5gclose_f                     (gr_id, STAT_CALL)

        allocate(Me%Int2DIn   (1:dims(1), 1:dims(2)))
        allocate(Me%Int3DIn   (1:dims(1), 1:dims(2), 1:dims(3)))
        allocate(Me%Float2DIn (1:dims(1), 1:dims(2)))
        allocate(Me%Float3DIn (1:dims(1), 1:dims(2), 1:dims(3)))
        
        allocate(Me%Int2DOut  (1:dims(2), 1:dims(1)))
        allocate(Me%Float2DOut(1:dims(2), 1:dims(1)))
        
        if (Me%DepthLayersON) then        
            allocate(Me%Float3DOut(1:dims(2), 1:dims(1), 1:Me%DepthLayers))        
            allocate(Me%Int3DOut  (1:dims(2), 1:dims(1), 1:Me%DepthLayers))
            allocate(Me%Depth3DIN (1:dims(1), 1:dims(2), 1:dims(3)))                    
        else
            allocate(Me%Float3DOut(1:dims(2), 1:dims(1), 1:dims(3)))        
            allocate(Me%Int3DOut  (1:dims(2), 1:dims(1), 1:dims(3)))            
        endif            

        Me%Int2DIn    = null_int
        Me%Int3DIn    = null_int  
        Me%Float2DIn  = null_real
        Me%Float3DIn  = null_real  

        Me%Int2DOut    = null_int
        Me%Int3DOut    = null_int  
        Me%Float2DOut  = null_real
        Me%Float3DOut  = null_real  

       
        Me%HDFFile%Size%ILB = 1
        Me%HDFFile%Size%JLB = 1
        Me%HDFFile%Size%KLB = 1

        Me%HDFFile%Size%IUB = dims(1)
        Me%HDFFile%Size%JUB = dims(2)
        Me%HDFFile%Size%KUB = dims(3)
        
        write(*,*)
        write(*,*)"IUB", Me%HDFFile%Size%IUB
        write(*,*)"JUB", Me%HDFFile%Size%JUB
        write(*,*)"KUB", Me%HDFFile%Size%KUB
        write(*,*)

    end subroutine ReadSetDimensions

    !--------------------------------------------------------------------------

    subroutine ReadWriteLatLon

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, i, j
        real, dimension(:,:), pointer       :: Lat, Lon, Lat_Stag, Lon_Stag, Aux4
        real(8), dimension(:,:), pointer    :: SphericX, SphericY               
        character(len=StringLength)         :: LatVar, LonVar

        !Begin-----------------------------------------------------------------

        write(*,*)"Reading and writing latitude and longitude..."

        allocate(Aux4        (1:Me%HDFFile%Size%IUB+1, 1:Me%HDFFile%Size%JUB+1))


        allocate(Lat        (1:Me%HDFFile%Size%JUB, 1:Me%HDFFile%Size%IUB))
        allocate(Lon        (1:Me%HDFFile%Size%JUB, 1:Me%HDFFile%Size%IUB))

        allocate(Lat_Stag   (1:Me%HDFFile%Size%JUB+1, 1:Me%HDFFile%Size%IUB+1))
        allocate(Lon_Stag   (1:Me%HDFFile%Size%JUB+1, 1:Me%HDFFile%Size%IUB+1))

        call HDF5SetLimits(Me%HDFFile%ObjHDF5, ILB = 1, IUB = Me%HDFFile%Size%IUB+1, &
                                               JLB = 1, JUB = Me%HDFFile%Size%JUB+1, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteLatLon - Convert2netcdf - ERR10'

        if(Me%HDFFile%ReadLatLon)then
            LatVar = "Latitude"
            LonVar = "Longitude"
        else
            LatVar = "ConnectionY"
            LonVar = "ConnectionX"
        end if

        call HDF5ReadData(HDF5ID       = Me%HDFFile%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(LatVar),                  &
                          Array2D      = Aux4,                          &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteLatLon - Convert2netcdf - ERR20'

        do j = 1, Me%HDFFile%Size%JUB+1
        do i = 1, Me%HDFFile%Size%IUB+1
            Lat_Stag(j,i) = Aux4(i,j)
        enddo
        enddo        

        call HDF5ReadData(HDF5ID       = Me%HDFFile%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(LonVar),                  &
                          Array2D      = Aux4,                          &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteLatLon - Convert2netcdf - ERR30'


        do j = 1, Me%HDFFile%Size%JUB+1
        do i = 1, Me%HDFFile%Size%IUB+1
            Lon_Stag(j,i) = Aux4(i,j)
        enddo
        enddo        

        do j = 1, Me%HDFFile%Size%JUB 
        do i = 1, Me%HDFFile%Size%IUB 
            Lat(j,i) = (Lat_Stag(j,i) + Lat_Stag(j,i+1)+Lat_Stag(j+1,i) + Lat_Stag(j+1,i+1))/4.
            Lon(j,i) = (Lon_Stag(j,i) + Lon_Stag(j,i+1)+Lon_Stag(j+1,i) + Lon_Stag(j+1,i+1))/4.            
        enddo
        enddo
        
        if (Me%SimpleGrid) then
        
            call NETCDFWriteLatLon1D(NCDFID           = Me%NCDF_File%ObjNETCDF,         &
                                     Lat              = Lat,                            &
                                     Lon              = Lon,                            &
                                     STAT             = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteLatLon - Convert2netcdf - ERR40'
        
        else
        
            allocate(SphericX(1:Me%HDFFile%Size%JUB+1, 1:Me%HDFFile%Size%IUB+1))
            allocate(SphericY(1:Me%HDFFile%Size%JUB+1, 1:Me%HDFFile%Size%IUB+1))
            
            if (Me%StagGoogleOut) then

                call WGS84toGoogleMaps(lon_stag, lat_stag,                              &
                                       1, Me%HDFFile%Size%JUB, 1, Me%HDFFile%Size%IUB,  &
                                       SphericX, SphericY)

                call NETCDFWriteLatLon(NCDFID           = Me%NCDF_File%ObjNETCDF,       &
                                       Lat              = Lat,                          &
                                       Lon              = Lon,                          &
                                       Lat_Stag         = Lat_Stag,                     &
                                       Lon_Stag         = Lon_Stag,                     &
                                       SphericX         = SphericX,                     &
                                       SphericY         = SphericY,                     &
                                       STAT             = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteLatLon - Convert2netcdf - ERR50'

                deallocate(SphericX, SphericY)
                
            else

                call NETCDFWriteLatLon(NCDFID           = Me%NCDF_File%ObjNETCDF,       &
                                       Lat              = Lat,                          &
                                       Lon              = Lon,                          &
                                       Lat_Stag         = Lat_Stag,                     &
                                       Lon_Stag         = Lon_Stag,                     &
                                       STAT             = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteLatLon - Convert2netcdf - ERR60'        
            
            endif
                    
        endif

        deallocate(Lat, Lon, Lat_Stag, Lon_Stag)
        nullify   (Lat, Lon, Lat_Stag, Lon_Stag)

        write(*,*)"Done!"
        write(*,*)

    end subroutine ReadWriteLatLon
    
    
    !--------------------------------------------------------------------------

    subroutine ReadLatLonWindow

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, i, j
        real, dimension(:,:), pointer       :: Lat, Lon, Lat_Stag, Lon_Stag, Aux4
        real(8), dimension(:,:), pointer    :: SphericX, SphericY               
        character(len=StringLength)         :: LatVar, LonVar

        !Begin-----------------------------------------------------------------

        write(*,*)"Reading latitude and longitude window..."

        allocate(Aux4        (1:Me%HDFFile%Size%IUB+1, 1:Me%HDFFile%Size%JUB+1))


        allocate(Lat        (1:Me%HDFFile%Size%JUB, 1:Me%HDFFile%Size%IUB))
        allocate(Lon        (1:Me%HDFFile%Size%JUB, 1:Me%HDFFile%Size%IUB))

        call HDF5SetLimits(Me%HDFFile%ObjHDF5, ILB = 1, IUB = Me%HDFFile%Size%IUB+1, &
                                               JLB = 1, JUB = Me%HDFFile%Size%JUB+1, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadLatLonWindow - Convert2netcdf - ERR10'

        if(Me%HDFFile%ReadLatLon)then
            LatVar = "Latitude"
            LonVar = "Longitude"
        else
            LatVar = "ConnectionY"
            LonVar = "ConnectionX"
        end if

        call HDF5ReadData(HDF5ID       = Me%HDFFile%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(LatVar),                  &
                          Array2D      = Aux4,                          &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadLatLonWindow - Convert2netcdf - ERR20'
        
        Me%NCDF_File%geospatial_lat_min = Aux4(1,1)
        Me%NCDF_File%geospatial_lat_max = Aux4(Me%HDFFile%Size%IUB+1,Me%HDFFile%Size%JUB+1)

        call HDF5ReadData(HDF5ID       = Me%HDFFile%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(LonVar),                  &
                          Array2D      = Aux4,                          &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadLatLonWindow - Convert2netcdf - ERR30'

        Me%NCDF_File%geospatial_lon_min = Aux4(1,1)
        Me%NCDF_File%geospatial_lon_max = Aux4(Me%HDFFile%Size%IUB+1,Me%HDFFile%Size%JUB+1)       

        deallocate(Lat, Lon)
        nullify   (Lat, Lon)

        write(*,*)"Done!"
        write(*,*)

    end subroutine ReadLatLonWindow
    
    
    !--------------------------------------------------------------------------

    subroutine ReadDepthIn3D(OutputNumber)
    
        !Arguments-------------------------------------------------------------    
        integer                             :: OutputNumber

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, i, j, k, KUBin
        real,    dimension(:,:,:), pointer  :: Vert3D
        integer, dimension(:,:,:), pointer  :: WaterPoints3D
        character(len=StringLength)         :: PointsVar

        !Begin-----------------------------------------------------------------

        write(*,*)"Reading and writing Vertical Coordinate..."
        
        KUBin   =  Me%HDFFile%Size%KUB

        !Same as SZZ
        allocate(Vert3D      (1:Me%HDFFile%Size%IUB, &
                              1:Me%HDFFile%Size%JUB, &
                              1:KUBin+1))

        allocate(WaterPoints3D (1:Me%HDFFile%Size%IUB, &
                                1:Me%HDFFile%Size%JUB, &
                                1:KUBin))
                                
            !Read VerticalZ code
        call HDF5SetLimits(Me%HDFFile%ObjHDF5, ILB = 1, IUB = Me%HDFFile%Size%IUB, &
                                               JLB = 1, JUB = Me%HDFFile%Size%JUB, &
                                               KLB = 1, KUB = KUBin+1,             &
                                               STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDepthIn3D - Convert2netcdf - ERR10'

        !The first time instant is considered only.
        call HDF5ReadData(HDF5ID       = Me%HDFFile%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(Me%HDFFile%VertVar),      &
                          Array3D      = Vert3D,                        &
                          OutputNumber = OutputNumber,                  &
                          STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDepthIn3D - Convert2netcdf - ERR20'
                
        !Read Waterpoints code
        call HDF5SetLimits(Me%HDFFile%ObjHDF5, ILB = 1, IUB = Me%HDFFile%Size%IUB, &
                                               JLB = 1, JUB = Me%HDFFile%Size%JUB, &
                                               KLB = 1, KUB = KUBin,               &
                                               STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDepthIn3D - Convert2netcdf - ERR30'

        PointsVar = "WaterPoints3D"

        !Get the land mask
        call HDF5ReadData(HDF5ID       = Me%HDFFile%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(PointsVar),               &
                          Array3D      = WaterPoints3D,                 &
                          STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDepthIn3D - Convert2netcdf - ERR40'


            !Returns the mean depth of the cell's top and down faces in a cartesian referential

        do k = 1, KUBin

            do j = 1, Me%HDFFile%Size%JUB
            do i = 1, Me%HDFFile%Size%IUB
            
                Me%Depth3DIN(i,j,k) = FillValueReal
                
                if(WaterPoints3D(i, j, k) == Waterpoint)then
                    if(Vert3D(i, j, k) > FillValueReal/2. .and. Vert3D(i,j,k+1) > FillValueReal/2.)then
                        Me%Depth3DIN(i,j,k) = (Vert3D(i, j, k) + Vert3D(i, j, k+1)) / 2.
                    endif
                endif
                
            enddo
            enddo
        
        enddo            

        deallocate(WaterPoints3D)
        nullify   (WaterPoints3D)

        deallocate(Vert3D)
        nullify   (Vert3D)


        write(*,*)"Done!"
        write(*,*)

    end subroutine ReadDepthIn3D

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine WriteDepthLayers

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, i, j, k, KUBout
        real, dimension(:  ), pointer       :: Vert1DStag


        !Begin-----------------------------------------------------------------

        write(*,*)"writing Vertical Coordinate..."
        
        KUBout  =  Me%DepthLayers
        
        allocate(Vert1DStag  (1:KUBout+1))        

            
        Vert1DStag(1:KUBout) = Me%DepthVector(1:KUBout)
        Vert1DStag(KUBout+1) = 0.
        

        call NETCDFWriteVert(NCDFID           = Me%NCDF_File%ObjNETCDF,                 &
                             Vert             = Me%DepthVector,                         &
                             VertCoordinate   = .false.,                                &
                             SimpleGrid       = Me%SimpleGrid,                          &
                             OffSet           = Me%DepthAddOffSet,                      &                             
                             STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteDepthLayers - Convert2netcdf - ERR10'
        
        if (.not.Me%SimpleGrid) then
        
        call NETCDFWriteVertStag(NCDFID         = Me%NCDF_File%ObjNETCDF,               &
                                 VertStag       = Vert1DStag,                           &
                                 STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteDepthLayers - Convert2netcdf - ERR20'                

        endif 
        
        deallocate(Vert1DStag)
        nullify   (Vert1DStag)
        
        write(*,*)"Done!"
        write(*,*)

    end subroutine WriteDepthLayers

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine ReadWriteVertical

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, i, j, k, kk
        real, dimension(:  ), pointer       :: Vert1D, Vert1DStag
        real, dimension(:,:,:), pointer     :: Vert3D
        integer, dimension(:,:,:),pointer   :: WaterPoints3D
        real                                :: Vert, area
        character(len=StringLength)         :: PointsVar


        !Begin-----------------------------------------------------------------

        write(*,*)"Reading and writing Vertical Coordinate..."
        
        allocate(Vert1D      (1:Me%HDFFile%Size%KUB))
        allocate(Vert1DStag  (1:Me%HDFFile%Size%KUB+1))        

        !Same as SZZ
        allocate(Vert3D      (1:Me%HDFFile%Size%IUB, &
                              1:Me%HDFFile%Size%JUB, &
                              1:Me%HDFFile%Size%KUB+1))

        allocate(WaterPoints3D (1:Me%HDFFile%Size%IUB, &
                                1:Me%HDFFile%Size%JUB, &
                                1:Me%HDFFile%Size%KUB))
                                
        if (.not. Me%HDFFile%Sigma) then                                

            !Read VerticalZ code
            call HDF5SetLimits(Me%HDFFile%ObjHDF5, ILB = 1, IUB = Me%HDFFile%Size%IUB, &
                                                   JLB = 1, JUB = Me%HDFFile%Size%JUB, &
                                                   KLB = 1, KUB = Me%HDFFile%Size%KUB+1, &
                                                   STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - Convert2netcdf - ERR01'

            !The first time instant is considered only.
            call HDF5ReadData(HDF5ID       = Me%HDFFile%ObjHDF5,            &
                              GroupName    = "/Grid",                       &
                              Name         = trim(Me%HDFFile%VertVar),      &
                              Array3D      = Vert3D,                        &
                              OutputNumber = 1,                             &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - Convert2netcdf - ERR02'

        endif
                
        !Read Waterpoints code
        call HDF5SetLimits(Me%HDFFile%ObjHDF5, ILB = 1, IUB = Me%HDFFile%Size%IUB, &
                                               JLB = 1, JUB = Me%HDFFile%Size%JUB, &
                                               KLB = 1, KUB = Me%HDFFile%Size%KUB, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - Convert2netcdf - ERR03'

        PointsVar = "WaterPoints3D"

        !Get the land mask
        call HDF5ReadData(HDF5ID       = Me%HDFFile%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(PointsVar),               &
                          Array3D      = WaterPoints3D,                 &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - Convert2netcdf - ERR04'

        !Is it a sigma vertical coordinate ?
        if (Me%HDFFile%Sigma) then
            
            !Returns the depth of the cell's bottom face in a sigma referential i.e. [0 1]
            do i = 1, Me%HDFFile%Size%KUB
                Vert1D(i) = i/Me%HDFFile%Size%KUB - 0.5  !assuming constant spacing
            enddo

            !Returns the depth of the cell's bottom face in a sigma referential i.e. [0 1]
            do i = 0, Me%HDFFile%Size%KUB
                Vert1DStag(i+1) = i/Me%HDFFile%Size%KUB  !assuming constant spacing
            enddo

        !No? Then it must be a z-vertical coordinate.
        else
        
            !Returns the mean depth of the cell's top and down faces in a cartesian referential

            do k = 1, Me%HDFFile%Size%KUB

                Vert = 0.
                area = 0.

                do j = 1, Me%HDFFile%Size%JUB
                do i = 1, Me%HDFFile%Size%IUB
                    if(WaterPoints3D(i, j, k) == Waterpoint)then
                        if(Vert3D(i, j, k) > FillValueReal/2. .and. Vert3D(i,j,k+1) > FillValueReal/2.)then
                            Vert = Vert + Vert3D(i, j, k) + Vert3D(i, j, k+1)
                            area = area + 1 
                        endif
                    endif
                enddo
                enddo
            
                if(area > 0.)then
                    Vert1D(k) = Vert / area * 0.5
                else
                    Vert1D(k) = maxval(Vert3D)
                endif

            enddo
            
            do k = 1, Me%HDFFile%Size%KUB+1

                Vert = 0.
                area = 0.

                do j = 1, Me%HDFFile%Size%JUB
                do i = 1, Me%HDFFile%Size%IUB
                    if (k==Me%HDFFile%Size%KUB+1) then
                        kk = k-1
                    else
                        kk = k
                    endif
                    if(WaterPoints3D(i, j, kk) == Waterpoint)then
                        if(Vert3D(i, j, k) > FillValueReal/2.)then
                            Vert = Vert + Vert3D(i, j, k) 
                            area = area + 1 
                        endif
                    endif
                enddo
                enddo
            
                if(area > 0.)then
                    Vert1DStag(k) = Vert / area 
                else
                    Vert1DStag(k) = maxval(Vert3D)
                endif

            enddo            
       
        endif

        call NETCDFWriteVert(NCDFID           = Me%NCDF_File%ObjNETCDF,                 &
                             Vert             = Vert1D,                                 &
                             VertCoordinate   = Me%HDFFile%Sigma,                       &
                             SimpleGrid       = Me%SimpleGrid,                          &
                             OffSet           = Me%DepthAddOffSet,                      &
                             STAT             = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - Convert2netcdf - ERR05'
        
        if (.not.Me%SimpleGrid) then        
        
        call NETCDFWriteVertStag(NCDFID         = Me%NCDF_File%ObjNETCDF,               &
                                 VertStag       = Vert1DStag,                           &
                                 STAT           = STAT_CALL)
        endif

        deallocate(WaterPoints3D)
        nullify   (WaterPoints3D)

        deallocate(Vert3D)
        nullify   (Vert3D)

        deallocate(Vert1D)
        nullify   (Vert1D)
        
        deallocate(Vert1DStag)
        nullify   (Vert1DStag)        

        write(*,*)"Done!"
        write(*,*)

    end subroutine ReadWriteVertical

    !--------------------------------------------------------------------------
    subroutine WriteVerticalNullDepth

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, i, j, k, kk
        real, dimension(:  ), pointer       :: Vert1D, Vert1DStag


        !Begin-----------------------------------------------------------------
    
    
        write(*,*)"Writing Vertical Coordinate..."
        
        if (Me%HDFFile%Size%KUB /= 1) then
            stop "WriteVerticalNullDepth - ConvertNetCDF - ERR10"
        endif
        
        allocate(Vert1D      (1:Me%HDFFile%Size%KUB))
        allocate(Vert1DStag  (1:Me%HDFFile%Size%KUB+1))          
        
        Vert1D(1) = 0.5
        Vert1DStag(1) = 1
        Vert1DStag(2) = 0
        

        call NETCDFWriteVert(NCDFID           = Me%NCDF_File%ObjNETCDF,                 &
                             Vert             = Vert1D,                                 &
                             VertCoordinate   = Me%HDFFile%Sigma,                       &
                             SimpleGrid       = Me%SimpleGrid,                          &
                             OffSet           = Me%DepthAddOffSet,                      &
                             STAT             = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteVerticalNullDepth - Convert2netcdf - ERR10'
        
        if (.not.Me%SimpleGrid) then        
        
        call NETCDFWriteVertStag(NCDFID         = Me%NCDF_File%ObjNETCDF,               &
                                 VertStag       = Vert1DStag,                           &
                                 STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteVerticalNullDepth - Convert2netcdf - ERR20'                

        endif
        
        deallocate(Vert1D      )
        deallocate(Vert1DStag  )
        
        
    end subroutine WriteVerticalNullDepth
    
    !--------------------------------------------------------------------------    

    subroutine ReadWriteBathymetry

        !Local-----------------------------------------------------------------
        character(len=StringLength)         :: NCDFName, LongName, StandardName, Units, Positive
        real                                :: MinValue, MaxValue, ValidMin, ValidMax, MissingValue
        integer                             :: STAT_CALL, i, j

        !Begin-----------------------------------------------------------------

        write(*,*)"Reading and writing bathymetry..."
        
        call HDF5SetLimits(Me%HDFFile%ObjHDF5, ILB = 1, IUB = Me%HDFFile%Size%IUB, &
                                               JLB = 1, JUB = Me%HDFFile%Size%JUB, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteBathymetry - Convert2netcdf - ERR01'

        call HDF5ReadData(HDF5ID       = Me%HDFFile%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = "Bathymetry",                  &
                          Array2D      = Me%Float2DIn,                  &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteBathymetry - Convert2netcdf - ERR01'
        
        do i=1,Me%HDFFile%Size%IUB
        do j=1,Me%HDFFile%Size%JUB
            Me%Float2DOut(j,i) = Me%Float2DIn(i, j)
        enddo
        enddo
        
       Units = "m"

        call BuildAttributes("Bathymetry", NCDFName, LongName, StandardName, &
                                           Units, ValidMin, ValidMax,        &
                                           MinValue, MaxValue, MissingValue, &
                                           Positive, Me%Float2DOut)
        
        call NETCDFWriteData (NCDFID        = Me%NCDF_File%ObjNETCDF,   &
                              Name          = trim(NCDFName),           &
                              LongName      = trim(LongName),           &
                              StandardName  = trim(StandardName),       & 
                              Units         = trim(Units),              &
                              ValidMin      = ValidMin,                 &
                              ValidMax      = ValidMax,                 &
                              MinValue      = MinValue,                 &
                              MaxValue      = MaxValue,                 &
                              MissingValue  = MissingValue,             &
                              Array2D       = Me%Float2DOut,            &
                              STAT          = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteBathymetry - Convert2netcdf - ERR01'

        write(*,*)"Done!"
        write(*,*)

    end subroutine ReadWriteBathymetry

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadMask

        !Local-----------------------------------------------------------------
        character(len=StringLength)         :: NCDFName, LongName, StandardName, Units
        real                                :: MinValue, MaxValue, ValidMin, ValidMax, MissingValue
        integer                             :: STAT_CALL, i, j, k, KUB, kfloor

        !Begin-----------------------------------------------------------------

        write(*,*)"Reading and writing mask..."
        
        Units ="-"
        
        KUB = Me%HDFFile%Size%KUB

        if ( Me%HDFFile%ImposeMask ) then
        
            if ( Me%HDFFile%HdfMaskIs3D ) then

                call HDF5SetLimits(Me%HDFFile%ObjHDF5, ILB = 1, IUB = Me%HDFFile%Size%IUB,  &
                                               JLB = 1, JUB = Me%HDFFile%Size%JUB,          &
                                               KLB = 1, KUB = Me%HDFFile%Size%KUB,          &
                                               STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR10'
                
                call HDF5ReadData(HDF5ID       = Me%HDFFile%ObjHDF5,                        &
                                  GroupName    = Me%HDFFile%SizeGroup,                      &
                                  Name         = trim(Me%HDFFile%HdfMask),                  &
                                  Array3D      = Me%Int3DIn,                                &
                                  STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR50'
                
                if (Me%HDFFile%OutputIs2D) then
                
                    do i=1,Me%HDFFile%Size%IUB
                    do j=1,Me%HDFFile%Size%JUB
                        k                = 1
                        Me%Int2DIn (i,j) = Me%Int3DIn(i,j,k)                        
                        Me%Int2DOut(j,i) = Me%Int3DIn(i,j,k)
                    enddo                
                    enddo
            
                    call BuildAttributes(trim(Me%HDFFile%HdfMask), NCDFName,                    &
                                         LongName, StandardName,                                &
                                         Units, ValidMin, ValidMax,                             &
                                         MinValue, MaxValue, MissingValue, Int2D = Me%Int2DOut)
            
                    call NETCDFWriteData (NCDFID        = Me%NCDF_File%ObjNETCDF,               &
                                          Name          = trim(NCDFName),                       &
                                          LongName      = trim(LongName),                       &
                                          StandardName  = trim(StandardName),                   & 
                                          Units         = trim(Units),                          &
                                          ValidMin      = ValidMin,                             &
                                          ValidMax      = ValidMax,                             &
                                          MinValue      = MinValue,                             &
                                          MaxValue      = MaxValue,                             &
                                          MissingValue  = MissingValue,                         &
                                          Array2D       = Me%Int2DOut,                          &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR60'
                    
                else
                
                    if (Me%DepthLayersON) then
                    
                        call ReadDepthIn3D(1)
                    
                        do i=1,Me%HDFFile%Size%IUB
                        do j=1,Me%HDFFile%Size%JUB

                            kfloor = FillValueInt
                            do k=1,Me%HDFFile%Size%KUB                                        
                                if (Me%Depth3DIN(i,j,k) > HalfFillValueReal) then
                                    kfloor = k
                                    exit
                                endif                                    
                            enddo                                    
                            
                            Me%Int3DOut(j,i,:)=0
                            
                            if (kfloor > FillValueInt) then
                            
                                do k=1,Me%DepthLayers
                                    if (Me%DepthVector(k) < Me%Depth3DIN(i,j,kfloor)) then
                                        Me%Int3DOut(j,i,k) = 1
                                    endif                                        
                                enddo
                            endif                                
                        enddo                
                        enddo
                        
                    else                        
                        do i=1,Me%HDFFile%Size%IUB
                        do j=1,Me%HDFFile%Size%JUB
                        do k=1,Me%HDFFile%Size%KUB                
                            Me%Int3DOut(j,i,k) = Me%Int3DIn(i,j,k)
                        enddo
                        enddo                
                        enddo                    
                    endif
                    
                    call BuildAttributes(trim(Me%HDFFile%HdfMask), NCDFName,                    &
                                         LongName, StandardName,                                &
                                         Units, ValidMin, ValidMax,                             &
                                         MinValue, MaxValue, MissingValue, Int3D = Me%Int3DOut)
            
                    call NETCDFWriteData (NCDFID        = Me%NCDF_File%ObjNETCDF,               &
                                          Name          = trim(NCDFName),                       &
                                          LongName      = trim(LongName),                       &
                                          StandardName  = trim(StandardName),                   & 
                                          Units         = trim(Units),                          &
                                          ValidMin      = ValidMin,                             &
                                          ValidMax      = ValidMax,                             &
                                          MinValue      = MinValue,                             &
                                          MaxValue      = MaxValue,                             &
                                          MissingValue  = MissingValue,                         &
                                          Array3D       = Me%Int3DOut,                          &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR60'
                
                
                endif                    

            else

                call HDF5SetLimits(Me%HDFFile%ObjHDF5, ILB = 1, IUB = Me%HDFFile%Size%IUB,  &
                                               JLB = 1, JUB = Me%HDFFile%Size%JUB,          &
                                               STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR70'

                call HDF5ReadData(HDF5ID   = Me%HDFFile%ObjHDF5,                            &
                                  GroupName    = Me%HDFFile%SizeGroup,                      &
                                  Name         = Me%HDFFile%HdfMask,                        &
                                  Array2D      = Me%Int2DIn,                                &
                                  STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR80'
                
                do i=1,Me%HDFFile%Size%IUB
                do j=1,Me%HDFFile%Size%JUB
                    Me%Int2DOut(j,i) = Me%Int2DIn(i,j)
                enddo                
                enddo                

                call BuildAttributes(trim(Me%HDFFile%HdfMask), NCDFName, LongName,          &
                                     StandardName, Units, ValidMin, ValidMax,               &
                                     MinValue, MaxValue, MissingValue, Int2D = Me%Int2DOut)
        
                call NETCDFWriteData (NCDFID        = Me%NCDF_File%ObjNETCDF,               &
                                      Name          = trim(NCDFName),                       &
                                      LongName      = trim(LongName),                       &
                                      StandardName  = trim(StandardName),                   & 
                                      Units         = trim(Units),                          &
                                      ValidMin      = ValidMin,                             &
                                      ValidMax      = ValidMax,                             &
                                      MinValue      = MinValue,                             &
                                      MaxValue      = MaxValue,                             &
                                      MissingValue  = MissingValue,                         &
                                      Array2D       = Me%Int2DOut,                          &
                                      STAT          = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR90'

            endif

        else

            call ReadMaskOld

        endif

        write(*,*)"Done!"
        write(*,*)

    end subroutine ReadMask

    !--------------------------------------------------------------------------

    subroutine ReadMaskOld

        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: NCDFName, LongName, StandardName, Units
        real                                        :: MinValue, MaxValue, ValidMin, ValidMax, MissingValue
        integer                                     :: STAT_CALL, nItems
        integer                                     :: item
        character(len=StringLength)                 :: obj_name
        integer                                     :: rank, obj_type
        integer(HSIZE_T), dimension(7)              :: dims, maxdims
        integer(HID_T)                              :: gr_id, dset_id, class_id
        !integer(HID_T)                              :: ssize
        integer(SIZE_T)                             :: ssize        
        integer(HID_T)                              :: space_id, datatype_id
        !logical                                     :: IsMapping
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB, i, j, k

        !Begin-----------------------------------------------------------------

        ILB = Me%HDFFile%Size%ILB 
        IUB = Me%HDFFile%Size%IUB 
        JLB = Me%HDFFile%Size%JLB 
        JUB = Me%HDFFile%Size%JUB 
        KLB = Me%HDFFile%Size%KLB 
        KUB = Me%HDFFile%Size%KUB

        Me%IsMapping = .false.

        write(*,*)"Reading and writing mask old..."

        call h5gopen_f(Me%HDFFile%FileID, "/Grid", gr_id, STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR01'

        !Get number of items in group
        call h5gn_members_f(gr_id, "/Grid", nItems, STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR02'

        Units ="-"

        do item = 1, nItems

            !Get info on object
            call h5gget_obj_info_idx_f(gr_id, "/Grid", item-1, obj_name, obj_type, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR03'

            if (obj_type == H5G_DATASET_F) then

                call h5dopen_f      (gr_id, trim(adjustl(obj_name)), dset_id, STAT_CALL)
                call h5dget_space_f (dset_id, space_id, STAT_CALL)
                call h5sget_simple_extent_ndims_f (space_id, rank, STAT_CALL)
                call h5dget_type_f (dset_id, datatype_id,   STAT_CALL)
                call h5tget_size_f (datatype_id, ssize,      STAT_CALL)
                call h5tget_class_f(datatype_id, class_id,  STAT_CALL) 
                call h5tclose_f    (datatype_id, STAT_CALL) 
                call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT_CALL) 
                call h5sclose_f     (space_id, STAT_CALL)
                call h5dclose_f (dset_id, STAT_CALL)

                !This isn't the right way to look for a Mask.
                !The mask variable should be explicited in the configuration file.
                !And the defaults should be WaterPoints2D and WaterPoints3D.
                if(class_id == H5T_INTEGER_F .and. trim(obj_name)/= "Define Cells") then
                    Me%IsMapping = .true.
                else
                    Me%IsMapping = .false.
                endif

                if(Me%IsMapping)then

                    select case(Rank)

                        case(2)

                            call HDF5SetLimits(Me%HDFFile%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR04'

                            call HDF5ReadData(HDF5ID       = Me%HDFFile%ObjHDF5,            &
                                              GroupName    = "/Grid/",                      &
                                              Name         = trim(obj_name),                &
                                              Array2D      = Me%Int2DIn,                    &
                                              STAT         = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR05'
                            
                            do i=ILB, IUB
                            do j=JLB, JUB
                                Me%Int2DOut(j,i) = Me%Int2DIn(i,j)
                            enddo                
                            enddo                            

                            call BuildAttributes(obj_name, NCDFName, LongName, StandardName, &
                                                 Units, ValidMin, ValidMax,             &
                                                 MinValue, MaxValue, MissingValue,      &
                                                 Int2D = Me%Int2DOut)

                            if (trim(NCDFName) /= "mask") then
                            call CheckAndCorrectVarName(obj_name, NCDFName)
                            endif

                            call NETCDFWriteData (NCDFID        = Me%NCDF_File%ObjNETCDF,   &
                                                  Name          = trim(NCDFName),           &
                                                  LongName      = trim(LongName),           &
                                                  StandardName  = trim(StandardName),       & 
                                                  Units         = trim(Units),              &
                                                  ValidMin      = ValidMin,                 &
                                                  ValidMax      = ValidMax,                 &
                                                  MinValue      = MinValue,                 &
                                                  MaxValue      = MaxValue,                 &
                                                  MissingValue  = MissingValue,             &
                                                  Array2D       = Me%Int2DOut,              &
                                                  STAT          = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR06'


                        case(3)

                            call HDF5SetLimits(Me%HDFFile%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR07'

                            call HDF5ReadData(HDF5ID       = Me%HDFFile%ObjHDF5,            &
                                              GroupName    = "/Grid",                       &
                                              Name         = trim(obj_name),                &
                                              Array3D      = Me%Int3DIn,                    &
                                              STAT         = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR08'
                            
                            do i=ILB, IUB
                            do j=JLB, JUB
                            do k=KLB, KUB
                                Me%Int3DOut(j,i,k) = Me%Int3DIn(i,j,k)
                            enddo                
                            enddo                            
                            enddo

                            call BuildAttributes(obj_name, NCDFName, LongName, StandardName, &
                                                           Units, ValidMin, ValidMax,        &
                                                           MinValue, MaxValue, MissingValue, &
                                                           Int3D = Me%Int3DOut)

                            call NETCDFWriteData (NCDFID        = Me%NCDF_File%ObjNETCDF,   &
                                                  Name          = trim(NCDFName),           &
                                                  LongName      = trim(LongName),           &
                                                  StandardName  = trim(StandardName),       & 
                                                  Units         = trim(Units),              &
                                                  ValidMin      = ValidMin,                 &
                                                  ValidMax      = ValidMax,                 &
                                                  MinValue      = MinValue,                 &
                                                  MaxValue      = MaxValue,                 &
                                                  MissingValue  = MissingValue,             &
                                                  Array3D       = Me%Int3DOut,              &
                                                  STAT          = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadMask - Convert2netcdf - ERR09'

                    end select

                end if

            elseif (obj_type == H5G_GROUP_F) then


            endif

        enddo

        call h5gclose_f(gr_id, STAT_CALL)

        write(*,*)"Done!"
        write(*,*)

    end subroutine ReadMaskOld

    !-------------------------------------------------------------------------
    
    subroutine BuildAttributes(Name, NCDFName, LongName, StandardName, Units,           &
                               ValidMin, ValidMax, Min, Max, MissingValue, Positive,    &
                               Float2D, Float3D, Int2D, Int3D, Add_Factor, Multiply_Factor)

        !Arguments-------------------------------------------------------------
        character(len=*), intent(in )                       :: Name
        character(len=*), intent(out)                       :: NCDFName
        character(len=*), intent(out)                       :: LongName
        character(len=*), intent(out)                       :: StandardName
        character(len=*), intent(inout)                     :: Units
        character(len=*), intent(out), optional             :: Positive
        real,             intent(out)                       :: ValidMin, Min
        real,             intent(out)                       :: ValidMax, Max
        real,             intent(out)                       :: MissingValue
        real,    dimension(:,:  ), pointer, optional        :: Float2D 
        real,    dimension(:,:,:), pointer, optional        :: Float3D 
        integer, dimension(:,:  ), pointer, optional        :: Int2D   
        integer, dimension(:,:,:), pointer, optional        :: Int3D   
        real,                               optional        :: Add_Factor
        real,                               optional        :: Multiply_Factor

        !Local-----------------------------------------------------------------
        integer                                             :: i, j, k, KUBout
        real                                                :: Add_Factor_
        real                                                :: Multiply_Factor_
        character(len=StringLength)                         :: Units_
        !Begin-----------------------------------------------------------------
        
        if (present(Add_Factor)) then
            Add_Factor_      = Add_Factor
        else
            Add_Factor_      = 0.
        endif

        if (present(Multiply_Factor)) then
            Multiply_Factor_      = Multiply_Factor
        else
            Multiply_Factor_      = 1.
        endif

        if (Me%DepthLayersON) then
            KUBout = Me%DepthLayers   
        else
            KUBout = Me%HDFFile%Size%KUB
        endif          
        

        select case(trim(adjustl(Name)))

            case("Bathymetry")
                NCDFName        = "bathymetry"
                LongName        = "bathymetry"
                StandardName    = "sea_floor_depth_below_geoid"
                Units_          = "m"
                ValidMin        = -50.
                ValidMax        = 11000.
                Positive        = "down"
                MissingValue    = -99.0

            case("WaterPoints2D", "WaterPoints3D", "MappingPoints2D", "WaterPoints")
                NCDFName        = "mask"
                LongName        = "mask of potential water points"
                StandardName    = "land_binary_mask"
                Units_          = null_str
                ValidMin        = 0
                ValidMax        = 1
                MissingValue    = -99
            
            case("OpenPoints2D", "OpenPoints3D")
                NCDFName        = "mask"
                LongName        = "mask of effective water points at one given instant"
                StandardName    = "mask"
                Units_          = null_str
                ValidMin        = 0
                ValidMax        = 1
                MissingValue    = -99

            case("temperature")
                NCDFName        = "temperature"
                LongName        = "sea water temperature"
                StandardName    = "sea_water_temperature"
                Units_          = "degC"
                ValidMin        = 0.
                ValidMax        = 50.
                MissingValue    = Me%MissingValue

            case("salinity")
                NCDFName        = "salinity"
                LongName        = "sea water salinity"
                StandardName    = "sea_water_salinity"
                Units_          = "1e-3"
                ValidMin        = 0.
                ValidMax        = 40.
                MissingValue    = Me%MissingValue

            case("density")
                NCDFName        = "density"
                LongName        = "sea water density"
                StandardName    = "sea_water_density"
                Units_          = "kg m-3"
                ValidMin        = 900.
                ValidMax        = 1200.
                MissingValue    = Me%MissingValue

            case("oxygen")
                NCDFName        = "dissolved_oxygen"
                LongName        = "mass concentration of oxygen in sea water"
                StandardName    = "mass_concentration_of_oxygen_in_sea_water"
                Units_          = "mg l-1"
                ValidMin        = 0.
                ValidMax        = 30.
                MissingValue    = Me%MissingValue

            case("dissolved_oxygen_percent_saturation")
                NCDFName        = "dissolved_oxygen_percent_saturation"
                LongName        = "dissolved oxygen percent saturation"
                StandardName    = "dissolved_oxygen_percent_saturation"
                Units_          = "%"
                ValidMin        = 0.
                ValidMax        = 200.
                MissingValue    = Me%MissingValue

            case("velocity_U")
                NCDFName        = "u"
                LongName        = "eastward sea water velocity"
                StandardName    = "eastward_sea_water_velocity"
                Units_          = "m s-1"
                ValidMin        = -5.
                ValidMax        = 5.
                MissingValue    = Me%MissingValue

            case("velocity_V")
                NCDFName        = "v"
                LongName        = "northward sea water velocity"
                StandardName    = "northward_sea_water_velocity"
                Units_          = "m s-1"
                ValidMin        = -5.
                ValidMax        = 5.
                MissingValue    = Me%MissingValue

            case("velocity_W")
                NCDFName        = "velocity_W"
                LongName        = "upward sea water velocity"
                StandardName    = "upward_sea_water_velocity"
                Units_          = "m s-1"
                ValidMin        = -2.
                ValidMax        = 2.
                MissingValue    = Me%MissingValue

            case("velocity_modulus")
                NCDFName        = "vm"
                LongName        = "sea water speed"
                StandardName    = "sea_water_speed"
                Units_          = "m s-1"
                ValidMin        = -5.
                ValidMax        = 5.
                MissingValue    = Me%MissingValue

            case("water_level")
                NCDFName        = "ssh"
                LongName        = "sea surface height"
                StandardName    = "sea_surface_height"
                Units_          = "m"
                ValidMin        = -20.
                ValidMax        = 20.
                MissingValue    = Me%MissingValue

            case("wind_modulus")
                NCDFName        = "wind_modulus"
                LongName        = "wind speed"
                StandardName    = "wind_speed"
                Units_          = "m s-1"
                ValidMin        = 0.0
                ValidMax        = 100.
                MissingValue    = Me%MissingValue
            case("wind_gust")
                NCDFName        = "wind_speed_of_gust"
                LongName        = "wind speed of gust"
                StandardName    = "wind_speed_of_gust"
                Units_          = "m s-1"
                ValidMin        = 0.0
                ValidMax        = 200.
                MissingValue    = Me%MissingValue                
            case("wind_velocity_X")
                
                if (Me%OdysseaProject) then
                    NCDFName        = "grid_eastward_wind"
                    LongName        = "grid_eastward_wind"
                    StandardName    = "grid_eastward_wind"
                else
                NCDFName        = "x_wind"
                LongName        = "x wind"
                StandardName    = "x_wind"
                endif
                
                Units_          = "m s-1"
                ValidMin        = -100.
                ValidMax        = 100.
                MissingValue    = Me%MissingValue

            case("wind_velocity_Y")
                if (Me%OdysseaProject) then
                    NCDFName        = "grid_northward_wind"
                    LongName        = "grid_northward_wind"
                    StandardName    = "grid_northward_wind"
                else                
                NCDFName        = "y_wind"
                LongName        = "y wind"
                StandardName    = "y_wind"
                endif
                Units_          = "m s-1"
                ValidMin        = -100.
                ValidMax        = 100.
                MissingValue    = Me%MissingValue

            case("air_temperature")
                NCDFName        = "air_temperature"
                LongName        = "air temperature"
                StandardName    = "air_temperature"
                Units_          = "degC"
                ValidMin        = -90.
                ValidMax        = 60.
                MissingValue    = Me%MissingValue

            case("atmospheric_pressure")
                NCDFName        = "air_pressure_at_mean_sea_level"
                LongName        = "air_pressure_at_mean_sea_level"
                StandardName    = "air_pressure_at_mean_sea_level"
                Units_          = "Pa"
                ValidMin        = 85000.
                ValidMax        = 110000.
                MissingValue    = Me%MissingValue

            case("mean_sea_level_pressure")
                NCDFName        = "mean_sea_level_pressure"
                LongName        = "air pressure at sea level"
                StandardName    = "air_pressure_at_sea_level"
                Units_          = "Pa"
                ValidMin        = 85000.
                ValidMax        = 110000.
                MissingValue    = Me%MissingValue

            case("relative_humidity")
                NCDFName        = "relative_humidity"
                LongName        = "relative humidity"
                StandardName    = "relative_humidity"
                Units_          = "1"
                ValidMin        = 0.
                ValidMax        = 1.
                MissingValue    = Me%MissingValue

            case("short_wave_solar_radiation_extinction")
                NCDFName        = "volume_absorption_coefficient_of_radiative_flux_in_sea_water"
                LongName        = "short wave solar radiation light extinction coefficient"
                StandardName    = "volume_absorption_coefficient_of_radiative_flux_in_sea_water"
                Units_          = "m-1"
                ValidMin        = 0.
                ValidMax        = 100.
                MissingValue    = Me%MissingValue

            case("solar_radiation")
                NCDFName        = "solar_radiation"
                LongName        = "downwelling shortwave flux in air"
                StandardName    = "downwelling_shortwave_flux_in_air"
                Units_          = "W m-2"
                ValidMin        = 0.
                ValidMax        = 1400.
                MissingValue    = Me%MissingValue

            case("downward_long_wave_radiation")
                NCDFName        = "downward_long_wave_radiation"
                LongName        = "downwelling longwave flux in air"
                StandardName    = "downwelling_longwave_flux_in_air"
                Units_          = "W m-2"
                ValidMin        = 0.
                ValidMax        = 1400.
                MissingValue    = Me%MissingValue

            case("phytoplankton")
                NCDFName        = "phytoplankton"
                LongName        = "mole concentration of phytoplankton expressed as carbon in sea water"
                StandardName    = "mole_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water"
                Units_          = "millimol m-3"
                Multiply_Factor_ = 1000./12.0107
                ValidMin        = 0. * Multiply_Factor_
                ValidMax        = 10. * Multiply_Factor_
                MissingValue    = Me%MissingValue

            case("zooplankton")
                NCDFName        = "zooplankton"
                LongName        = "mole concentration of zooplankton expressed as carbon in sea water"
                StandardName    = "mole_concentration_of_zooplankton_expressed_as_carbon_in_sea_water"
                Units_          = "millimol m-3"
                Multiply_Factor_ = 1000./12.0107
                ValidMin        = 0. * Multiply_Factor_
                ValidMax        = 10. * Multiply_Factor_
                MissingValue    = Me%MissingValue

            case("nitrate")
                NCDFName        = "nitrate"
                LongName        = "mole concentration of nitrate in sea water"
                StandardName    = "mole_concentration_of_nitrate_in_sea_water"
                Units_          = "millimol m-3"
                Multiply_Factor_ = 1000./14.0067
                ValidMin        = 0. * Multiply_Factor_
                ValidMax        = 10. * Multiply_Factor_
                MissingValue    = Me%MissingValue

            case("ammonia")
                NCDFName        = "ammonia"
                LongName        = "mole concentration of ammonium in sea water"
                StandardName    = "mole_concentration_of_ammonium_in_sea_water"
                Units_          = "millimol m-3"
                Multiply_Factor_ = 1000./14.0067
                ValidMin        = 0. * Multiply_Factor_
                ValidMax        = 10. * Multiply_Factor_
                MissingValue    = Me%MissingValue

            case("cohesive_sediment")
                NCDFName        = "cohesive_sediment"
                LongName        = "mass concentration of suspended matter in sea water"
                StandardName    = "mass_concentration_of_suspended_matter_in_sea_water"
                Units_          = "mg l-1"
                ValidMin        = 0.
                ValidMax        = 1000.
                MissingValue    = Me%MissingValue

            case("inorganic_phosphorus")
                NCDFName        = "inorganic_phosphorus"
                LongName        = "mole concentration of phosphate in sea water"
                StandardName    = "mole_concentration_of_phosphate_in_sea_water"
                Units_          = "millimol m-3"
                Multiply_Factor_ = 1000./30.974
                ValidMin        = 0. * Multiply_Factor_
                ValidMax        = 10. * Multiply_Factor_
                MissingValue    = Me%MissingValue

            case("particulate_organic_nitrogen")
                NCDFName        = "particulate_organic_nitrogen"
                LongName        = "mole concentration of particulate organic matter expressed as nitrogen in sea water"
                StandardName    = "mole_concentration_of_particulate_organic_matter_expressed_as_nitrogen_in_sea_water"
                Units_          = "millimol m-3"
                Multiply_Factor_ = 1000./14.0067
                ValidMin        = 0. * Multiply_Factor_
                ValidMax        = 10. * Multiply_Factor_
                MissingValue    = Me%MissingValue

            case("particulate_organic_phosphorus")
                NCDFName        = "particulate_organic_phosphorus"
                LongName        = "mole concentration of particulate organic matter expressed as phosphorus in sea water"
                StandardName    = "mole_concentration_of_particulate_organic_matter_expressed_as_phosphorus_in_sea_water"
                Units_          = "millimol m-3"
                Multiply_Factor_ = 1000./30.974
                ValidMin        = 0. * Multiply_Factor_
                ValidMax        = 10. * Multiply_Factor_
                MissingValue    = Me%MissingValue

            case("mean_wave_direction")
                NCDFName        = "mean_wave_direction"
                LongName        = "sea surface wave to direction"
                StandardName    = "sea_surface_wave_to_direction"
                Units_          = "degree"
                ValidMin        = 0. 
                ValidMax        = 360.
                MissingValue    = Me%MissingValue
                
            case("swell_wave_direction")
                NCDFName        = "swell_wave_direction"
                LongName        = "sea surface swell wave to direction"
                StandardName    = "sea_surface_swell_wave_to_direction"
                Units_          = "degree"
                ValidMin        = 0. 
                ValidMax        = 360.
                MissingValue    = Me%MissingValue                
                
! Next 2 keywords are not CF compliant
            case("mean_wave_direction_X")
                NCDFName        = "wave_X"
                LongName        = "mean wave direction X"
                StandardName    = "mean_wave_direction_X"
                Units_          = "-"
                ValidMin        = -1. 
                ValidMax        = 1.
                MissingValue    = Me%MissingValue

            case("mean_wave_direction_Y")
                NCDFName        = "wave_Y"
                LongName        = "mean wave direction Y"
                StandardName    = "mean_wave_direction_Y"
                Units_          = "-"
                ValidMin        = -1. 
                ValidMax        = 1.
                MissingValue    = Me%MissingValue

            case("significant_wave_height")
                NCDFName        = "significant_wave_height"
                LongName        = "sea surface wave significant height"
                StandardName    = "sea_surface_wave_significant_height"
                Units_          = "m"
                ValidMin        = 0. 
                ValidMax        = 40.
                MissingValue    = Me%MissingValue

            case("swell_wave_height")
                NCDFName        = "swell_wave_height"
                LongName        = "sea surface swell wave height"
                StandardName    = "sea_surface_swell_wave_height"
                Units_          = "m"
                ValidMin        = 0. 
                ValidMax        = 40.
                MissingValue    = Me%MissingValue

            case("wind_wave_height")
                NCDFName        = "wind_wave_height"
                LongName        = "sea surface wind wave height"
                StandardName    = "sea_surface_wind_wave_height"
                Units_          = "m"
                ValidMin        = 0. 
                ValidMax        = 40.
                MissingValue    = Me%MissingValue

            case("mean_wave_period")
                NCDFName        = "mean_wave_period"
                LongName        = "sea surface wave zero upcrossing period"
                StandardName    = "sea_surface_wave_zero_upcrossing_period"
                Units_          = "s"
                ValidMin        = 0. 
                ValidMax        = 30.
                MissingValue    = Me%MissingValue


            case("swell_wave_period")
                NCDFName        = "swell_wave_period"
                LongName        = "sea surface swell wave period"
                StandardName    = "sea_surface_swell_wave_period"
                Units_          = "s"
                ValidMin        = 0. 
                ValidMax        = 30.
                MissingValue    = Me%MissingValue

            case("wind_wave_period")
                NCDFName        = "wind_wave_period"
                LongName        = "sea surface wind wave period"
                StandardName    = "sea_surface_wind_wave_period"
                Units_          = "s"
                ValidMin        = 0. 
                ValidMax        = 30.
                MissingValue    = Me%MissingValue
                
            case default

                NCDFName        = trim(adjustl(Name))
                LongName        = trim(adjustl(Name))
                StandardName    = trim(adjustl(Name))
                Units_          = "unknown"
                ValidMin        =   Me%MissingValue / 10.
                ValidMax        = - Me%MissingValue / 10. 
                MissingValue    = Me%MissingValue

        end select

        Min = ValidMax
        Max = ValidMin

        
        if (Me%MohidStandardInOutUnits) then                
            
            Units = Units_
            
        else
            
            NCDFName        = trim(adjustl(Name))
            LongName        = trim(adjustl(Name))
            StandardName    = trim(adjustl(Name))
            ValidMin        =   Me%MissingValue / 10.
            ValidMax        = - Me%MissingValue / 10. 
            MissingValue    = Me%MissingValue
            !Units attributte is read before 

            if (present(Add_Factor)) then
                Add_Factor_      = Add_Factor
            else
                Add_Factor_      = 0.
            endif

            if (present(Multiply_Factor)) then
                Multiply_Factor_      = Multiply_Factor
            else
                Multiply_Factor_      = 1.
            endif  
        
            
        endif        

if1:   if(present(Int2D) .or. present(Int3D))then
           
            Min = 0
            Max = 1

        elseif(present(Float2D))then if1

            do j = 1, Me%HDFFile%Size%JUB
            do i = 1, Me%HDFFile%Size%IUB
                Float2D(j,i) = Float2D(j,i) * Multiply_Factor_ + Add_Factor_

                if(Float2D(j,i) .gt. FillValueReal/2. .and. Float2D(j,i) .lt.  Min .and. &
                                                            Float2D(j,i) .ge. ValidMin)then

                    Min = Float2D(j,i)

                end if

                if(Float2D(j,i) > Max .and. Float2D(j,i) .le. ValidMax) then

                    Max = Float2D(j,i) 

                endif
                
                if(Float2D(j,i) < FillValueReal/2.) then
                    Float2D(j,i) = Me%MissingValue
                endif
                
                if (Me%HDFFile%ImposeMask .or. Me%IsMapping) then
                    if (Me%Int2DOut(j,i) == 0) Float2D(j,i) = Me%MissingValue             
                endif                    

            enddo
            enddo

        elseif(present(Float3D))then if1

            do j = 1, Me%HDFFile%Size%JUB
            do i = 1, Me%HDFFile%Size%IUB
            do k = 1, KUBout
                Float3D(j,i,k) = Float3D(j, i, k) * Multiply_Factor_ + Add_Factor_

                if(Float3D(j,i,k) .gt. FillValueReal/2. .and. Float3D(j,i,k) .lt.  Min .and. &
                                                              Float3D(j,i,k) .ge. ValidMin)then

                    Min = Float3D(j,i,k)

                end if

                if(Float3D(j,i,k) > Max .and. Float3D(j,i,k) .le. ValidMax) then

                    Max = Float3D(j,i,k) 

                endif
                
                
                if(Float3D(j,i,k) < FillValueReal/2.) then
                    Float3D(j,i,k) = Me%MissingValue
                endif
                
                if (Me%HDFFile%ImposeMask .or. Me%IsMapping) then
                    if (Me%Int3DOut(j,i,k) == 0) Float3D(j,i,k) = Me%MissingValue             
                endif                                    

            enddo
            enddo
            enddo

        endif if1

    end subroutine BuildAttributes
    
    
    !--------------------------------------------------------------------------

    type(T_Time) function HDF5TimeInstant(Instant)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer             :: TimeVector

        !Begin-----------------------------------------------------------------

        call HDF5SetLimits  (Me%HDFFile%ObjHDF5, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = Me%HDFFile%ObjHDF5,                       &
                             GroupName      = '/'//trim(Me%HDFFile%TimeGroupName),      &
                             Name           = trim(Me%HDFFile%TimeVar),                 &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'HDF5TimeInstant - Convert2netcdf - ERR01'

        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))

        deallocate(TimeVector)
        nullify   (TimeVector)

    end function HDF5TimeInstant
    
    !--------------------------------------------------------------------------

    subroutine ReadWriteData

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: item
        integer                                             :: nGItems
        character(len=StringLength)                         :: obj_name
        integer                                             :: obj_type
        integer(HID_T)                                      :: gr_id 

        !Begin-----------------------------------------------------------------

        call h5gopen_f(Me%HDFFile%FileID, "/", gr_id, STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteData - Convert2netcdf - ERR01'

        !Get number of items in group
        call h5gn_members_f(gr_id, "/", nGItems, STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteData - Convert2netcdf - ERR02'

        if(Me%ConvertEverything)then

            do item = 1, nGItems

                !Get info on object
                call h5gget_obj_info_idx_f(gr_id, "/", item-1, obj_name, obj_type, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteData - Convert2netcdf - ERR03'
            
                if (obj_type == H5G_GROUP_F) then

                    call BrowseGroup(obj_name)

                end if

            enddo

        else

            do item = 1, Me%nGroupsToConvert
                call BrowseGroup(trim(Me%GroupsToConvert(item)))
            enddo

        end if

        call h5gclose_f(gr_id, STAT_CALL)

        write(*,*)"Done!"
        write(*,*)

    end subroutine ReadWriteData
   
    !--------------------------------------------------------------------------

    recursive subroutine BrowseGroup(FatherGroupName)

        !Arguments-------------------------------------------------------------
        character(len=*)                                    :: FatherGroupName

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: item
        integer                                             :: nGItems
        character(len=StringLength)                         :: obj_name = null_str
        integer                                             :: obj_type = null_int
        integer(HID_T)                                      :: gr_id

        !Begin-----------------------------------------------------------------

        call h5gopen_f(Me%HDFFile%FileID, "/"//trim(FatherGroupName), gr_id, STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'BrowseGroup - Convert2netcdf - ERR01'

        !Get number of items in group
        call h5gn_members_f(gr_id, "/"//trim(FatherGroupName), nGItems, STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'BrowseGroup - Convert2netcdf - ERR02'

        write(*,*) nGItems, " items -> "//"/"//trim(FatherGroupName)

        do item = 1, nGItems

            !Get info on object
            call h5gget_obj_info_idx_f(gr_id, "/"//trim(FatherGroupName), item-1, obj_name, obj_type, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'BrowseGroup - Convert2netcdf - ERR03'

            if    (obj_type == H5G_GROUP_F)then
                
                !Only reads groups and datasets not included in Grid and Time
                if (trim(FatherGroupName) .ne. "Grid" .and. trim(FatherGroupName) .ne. "Time") then
                    call BrowseGroup(trim(FatherGroupName)//"/"//trim(obj_name))
                endif

            elseif(obj_type == H5G_DATASET_F)then
                
                !Only reads groups and datasets not included in Grid and Time
                if (trim(FatherGroupName) .ne. "Grid" .and. trim(FatherGroupName) .ne. "Time") then
                    call ReadDataSet(obj_name, gr_id, item, nGItems) 
                endif

            end if

        enddo

    end subroutine BrowseGroup

    !--------------------------------------------------------------------------
    
    subroutine ReadDataSet(obj_name, gr_id, item, nGItems)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: obj_name
        integer(HID_T)                              :: gr_id
        integer                                     :: item, nGItems

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer(HID_T)                              :: class_id, space_id, dset_id
        integer(HID_T)                              :: datatype_id, rank, NumType
        integer(HID_T)                              :: attr_id, type_id
        integer(SIZE_T)                             :: ssize        
        integer(HSIZE_T), dimension(7)              :: dims
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        character(len=StringLength)                 :: Name, NCDFName, LongName, StandardName, Units
        real                                        :: MinValue, MaxValue, ValidMin, ValidMax, MissingValue
        integer                                     :: i, j, k, KUBout, kfloor
        real                                        :: DecimalPlaces, Z
        real(8),    dimension(:), pointer           :: Depth1D, Matrix1D
        
        
        
        
        

        !Begin-----------------------------------------------------------------
        
        ILB = Me%HDFFile%Size%ILB 
        IUB = Me%HDFFile%Size%IUB 
        JLB = Me%HDFFile%Size%JLB 
        JUB = Me%HDFFile%Size%JUB 
        KLB = Me%HDFFile%Size%KLB 
        KUB = Me%HDFFile%Size%KUB

        dims(1) = Me%HDFFile%Size%IUB - Me%HDFFile%Size%ILB + 1
        dims(2) = Me%HDFFile%Size%JUB - Me%HDFFile%Size%JLB + 1
        dims(3) = Me%HDFFile%Size%KUB - Me%HDFFile%Size%KLB + 1


        call h5dopen_f (gr_id, trim(adjustl(obj_name)), dset_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadDataSet - Convert2netcdf - ERR10'
        endif

        call h5dget_space_f                 (dset_id,       space_id,       STAT_CALL)
        call h5sget_simple_extent_ndims_f   (space_id,      rank,           STAT_CALL)
        call h5dget_type_f                  (dset_id,       datatype_id,    STAT_CALL)
        call h5tget_size_f                  (datatype_id, ssize,           STAT_CALL)
        call h5tget_class_f                 (datatype_id,   class_id,       STAT_CALL) 
        call h5tclose_f                     (datatype_id,                   STAT_CALL) 
        call h5sclose_f                     (space_id,                      STAT_CALL)
            
        if      (class_id == H5T_FLOAT_F  ) then
            NumType = H5T_NATIVE_REAL
        elseif  (class_id == H5T_INTEGER_F) then
            NumType = H5T_NATIVE_INTEGER
        end if

        call CheckAndCorrectVarName(obj_name, Name)

        !Reads Units
        call h5aopen_name_f     (dset_id, "Units", attr_id,     STAT_CALL)
        call h5Tcopy_f          (H5T_NATIVE_CHARACTER, type_id, STAT_CALL)
        call h5Tset_size_f      (type_id, Int8(StringLength),   STAT_CALL)
        call h5aread_f          (attr_id, type_id, Units, dims, STAT_CALL)
        call h5aclose_f         (attr_id,                       STAT_CALL) 
        call h5Tclose_f         (type_id,                       STAT_CALL)
        

        select case(rank) 

            case(2)

                call h5dread_f   (dset_id, NumType,                                         &
                                  Me%Float2DIn(Me%HDFFile%Size%ILB:Me%HDFFile%Size%IUB,     &
                                               Me%HDFFile%Size%JLB:Me%HDFFile%Size%JUB),    &
                                  dims, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'ReadDataSet - Convert2netcdf - ERR50'
                endif


                if (Me%DecimalPlaces > FillValueInt) then
                
                    DecimalPlaces = 10.**Me%DecimalPlaces

                    do i=1,Me%HDFFile%Size%IUB
                    do j=1,Me%HDFFile%Size%JUB
                        if (Me%Float2DIn(i, j) /= FillValueReal) then
                            Me%Float2DOut(j,i) = int(Me%Float2DIn(i, j) * DecimalPlaces)/ DecimalPlaces
                        else
                            Me%Float2DOut(j,i) =  FillValueReal                                                       
                        endif                            
                    enddo
                    enddo                
                    
                else                    
                    
                    do i=1,Me%HDFFile%Size%IUB
                    do j=1,Me%HDFFile%Size%JUB
                        if (Me%Float2DIn(i, j) /= FillValueReal) then
                            Me%Float2DOut(j,i) = Me%Float2DIn(i, j)
                        else
                            Me%Float2DOut(j,i) =  FillValueReal                                                       
                        endif                            
                    enddo
                    enddo                
                    
                endif                                        

                call BuildAttributes(Name, NCDFName, LongName, StandardName, &
                                           Units, ValidMin, ValidMax,        &
                                           MinValue, MaxValue, MissingValue, &
                                           Float2D = Me%Float2DOut,          &
                                           Add_Factor = Me%Add_Factor,       &
                                           Multiply_Factor = Me%Multiply_Factor)

                if(nGItems > 1)then

                    call NETCDFWriteData (NCDFID        = Me%NCDF_File%ObjNETCDF,   &
                                          Name          = trim(NCDFName),           &
                                          LongName      = trim(LongName),           &
                                          StandardName  = trim(StandardName),       & 
                                          Units         = trim(Units),              &
                                          ValidMin      = ValidMin,                 &
                                          ValidMax      = ValidMax,                 &
                                          MinValue      = MinValue,                 &
                                          MaxValue      = MaxValue,                 &
                                          MissingValue  = MissingValue,             &
                                          OutputNumber  = item,                     &
                                          Array2D       = Me%Float2DOut,            &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) then
                        stop 'ReadDataSet - Convert2netcdf - ERR60'
                    endif
                
                else
                    
                    call NETCDFWriteData (NCDFID        = Me%NCDF_File%ObjNETCDF,   &
                                          Name          = trim(NCDFName),           &
                                          LongName      = trim(LongName),           &
                                          StandardName  = trim(StandardName),       & 
                                          Units         = trim(Units),              &
                                          ValidMin      = ValidMin,                 &
                                          ValidMax      = ValidMax,                 &
                                          MinValue      = MinValue,                 &
                                          MaxValue      = MaxValue,                 &
                                          MissingValue  = MissingValue,             &
                                          Array2D       = Me%Float2DOut,             &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) then
                        stop 'ReadDataSet - Convert2netcdf - ERR70'
                    endif
                end if


            case(3)

                call h5dread_f   (dset_id, NumType,                                         &
                                  Me%Float3DIn(Me%HDFFile%Size%ILB:Me%HDFFile%Size%IUB,     &
                                               Me%HDFFile%Size%JLB:Me%HDFFile%Size%JUB,     &
                                               Me%HDFFile%Size%KLB:Me%HDFFile%Size%KUB),    &
                                  dims, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'ReadDataSet - Convert2netcdf - ERR80'
                endif

                if (Me%HDFFile%OutputIs2D) then                    

                    if (Me%DecimalPlaces > FillValueInt) then
                    
                        DecimalPlaces = 10.**Me%DecimalPlaces
                        
                        do i=1,Me%HDFFile%Size%IUB
                        do j=1,Me%HDFFile%Size%JUB
                            k                  = 1
                            Me%Float2DIn (i,j) = Me%Float3DIn(i, j, k)                            
                            if (Me%Float2DIn(i, j) /= FillValueReal) then                            
                                Me%Float2DOut(j,i) = int(Me%Float2DIn(i, j) * DecimalPlaces) / DecimalPlaces
                            else
                                Me%Float2DOut(j,i) =  FillValueReal                                                       
                            endif                                                                
                        enddo
                        enddo

                    else

                        do i=1,Me%HDFFile%Size%IUB
                        do j=1,Me%HDFFile%Size%JUB
                            k                  = 1
                            Me%Float2DIn (i,j) = Me%Float3DIn(i, j, k)
                            if (Me%Float2DIn(i, j) /= FillValueReal) then                                                        
                                Me%Float2DOut(j,i) = Me%Float2DIn(i, j)
                            else
                                Me%Float2DOut(j,i) = FillValueReal
                            endif                                
                        enddo
                        enddo
                        
                    endif                    

                    call BuildAttributes(Name, NCDFName, LongName, StandardName, &
                                               Units, ValidMin, ValidMax,        &
                                               MinValue, MaxValue, MissingValue, &
                                               Float2D = Me%Float2DOut,           &
                                               Add_Factor = Me%Add_Factor,       &
                                               Multiply_Factor = Me%Multiply_Factor)

                    if(nGItems > 1)then

                        call NETCDFWriteData (NCDFID        = Me%NCDF_File%ObjNETCDF,   &
                                              Name          = trim(NCDFName),           &
                                              LongName      = trim(LongName),           &
                                              StandardName  = trim(StandardName),       & 
                                              Units         = trim(Units),              &
                                              ValidMin      = ValidMin,                 &
                                              ValidMax      = ValidMax,                 &
                                              MinValue      = MinValue,                 &
                                              MaxValue      = MaxValue,                 &
                                              OutputNumber  = item,                     &
                                              Array2D       = Me%Float2DOut,            &
                                              STAT          = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) then
                            stop 'ReadDataSet - Convert2netcdf - ERR90'
                        endif
                    
                    else
                        
                        call NETCDFWriteData (NCDFID        = Me%NCDF_File%ObjNETCDF,   &
                                              Name          = trim(NCDFName),           &
                                              LongName      = trim(LongName),           &
                                              StandardName  = trim(StandardName),       & 
                                              Units         = trim(Units),              &
                                              ValidMin      = ValidMin,                 &
                                              ValidMax      = ValidMax,                 &
                                              MinValue      = MinValue,                 &
                                              MaxValue      = MaxValue,                 &
                                              Array2D       = Me%Float2DOut,            &
                                              STAT          = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) then
                            stop 'ReadDataSet - Convert2netcdf - ERR100'
                        endif
                    end if
                    
                else                 
                
                    if (Me%DepthLayersON) then
                        KUBout = Me%DepthLayers   
                    else
                        KUBout = Me%HDFFile%Size%KUB
                    endif                     
                    

                    if (Me%DepthLayersON) then                       
                    
                        call ReadDepthIn3D(item)
                        
                        allocate(Depth1D (Me%HDFFile%Size%KLB:Me%HDFFile%Size%KUB))
                        allocate(Matrix1D(Me%HDFFile%Size%KLB:Me%HDFFile%Size%KUB))
                                            
                        do i=1,Me%HDFFile%Size%IUB
                        do j=1,Me%HDFFile%Size%JUB

                            kfloor = FillValueInt
                            
                            do k=1,Me%HDFFile%Size%KUB                                        
                                if (Me%Depth3DIN(i,j,k) > HalfFillValueReal) then
                                    kfloor = k
                                    exit
                                endif                                    
                            enddo                                    
                            
                            Me%Float3DOut(j,i,:) = FillValueReal
                            
                            if (kfloor > FillValueInt) then
                                
                                !Depth = distance to surface
                                Depth1D (kfloor:Me%HDFFile%Size%KUB) = Me%Depth3DIN(i,j,kfloor:Me%HDFFile%Size%KUB)- Me%Depth3DIN(i,j,Me%HDFFile%Size%KUB)
                                Matrix1D(kfloor:Me%HDFFile%Size%KUB) = Me%Float3DIn(i,j,kfloor:Me%HDFFile%Size%KUB)
                                
                                do k=1,Me%DepthLayers
                                    Z = Me%DepthVector(k)
                                    if (Z < Me%Depth3DIN(i,j,kfloor)) then
                                        Me%Float3DOut(j,i,k) = ValueAtDepthZ(Z       = Z,            &
                                                                             KLB     = kfloor,       &
                                                                             KUB     = Me%HDFFile%Size%KUB, &
                                                                             Depth1D = Depth1D,      &
                                                                             Matrix1D= Matrix1D)
                                    endif                                        
                                enddo
                            endif                                
                        enddo                
                        enddo
                        
                        deallocate(Depth1D )
                        deallocate(Matrix1D)
                        

                    else                    

                        do i=1,Me%HDFFile%Size%IUB
                        do j=1,Me%HDFFile%Size%JUB
                        do k=1,Me%HDFFile%Size%KUB                
                            if (Me%Float3DIn(i, j, k) /= FillValueReal) then                                                              
                                Me%Float3DOut(j,i,k) = Me%Float3DIn(i, j, k)
                            else
                                Me%Float3DOut(j,i,k) = FillValueReal
                            endif                                    
                        enddo
                        enddo
                        enddo

                    endif

                    if (Me%DecimalPlaces > FillValueInt) then
                    
                        DecimalPlaces = 10.**Me%DecimalPlaces
                        
                        do i=1,Me%HDFFile%Size%IUB
                        do j=1,Me%HDFFile%Size%JUB
                        do k=1,KUBout  
                            if (Me%Float3DOut(j,i,k) /= FillValueReal) then                                      
                                Me%Float3DOut(j,i,k) = int(Me%Float3DOut(j,i,k) * DecimalPlaces) / DecimalPlaces
                            else                       
                                Me%Float3DOut(j,i,k) = FillValueReal         
                            endif                                
                        enddo
                        enddo
                        enddo

                    endif

                    call BuildAttributes(Name, NCDFName, LongName, StandardName, &
                                               Units, ValidMin, ValidMax,        &
                                               MinValue, MaxValue, MissingValue, &
                                               Float3D = Me%Float3DOut,           &
                                               Add_Factor = Me%Add_Factor,       &
                                               Multiply_Factor = Me%Multiply_Factor)                                               

                    if(nGItems > 1)then

                        call NETCDFWriteData (NCDFID        = Me%NCDF_File%ObjNETCDF,   &
                                              Name          = trim(NCDFName),           &
                                              LongName      = trim(LongName),           &
                                              StandardName  = trim(StandardName),       & 
                                              Units         = trim(Units),              &
                                              ValidMin      = ValidMin,                 &
                                              ValidMax      = ValidMax,                 &
                                              MinValue      = MinValue,                 &
                                              MaxValue      = MaxValue,                 &
                                              OutputNumber  = item,                     &
                                              Array3D       = Me%Float3DOut,            &
                                              STAT          = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) then
                            stop 'ReadDataSet - Convert2netcdf - ERR110'
                        endif
                    
                    else
                        
                        call NETCDFWriteData (NCDFID        = Me%NCDF_File%ObjNETCDF,   &
                                              Name          = trim(NCDFName),           &
                                              LongName      = trim(LongName),           &
                                              StandardName  = trim(StandardName),       & 
                                              Units         = trim(Units),              &
                                              ValidMin      = ValidMin,                 &
                                              ValidMax      = ValidMax,                 &
                                              MinValue      = MinValue,                 &
                                              MaxValue      = MaxValue,                 &
                                              Array3D       = Me%Float3DOut,            &
                                              STAT          = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) then
                            stop 'ReadDataSet - Convert2netcdf - ERR120'
                        endif
                    end if

                endif
            
         end select

        call h5dclose_f  (dset_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadDataSet - Convert2netcdf - ERR990'
        endif

    end subroutine ReadDataSet

    !--------------------------------------------------------------------------
    
    real(8) function ValueAtDepthZ(Z, KLB, KUB, Depth1D, Matrix1D)
    
        !Arguments-------------------------------------------------------------
        real(8), dimension(:), pointer :: Depth1D, Matrix1D
        real                           :: Z
        integer                        :: KLB, KUB

        !Local-----------------------------------------------------------------        
        integer     ::kb, Ndepths, k

        !Begin-----------------------------------------------------------------            

        do k = KUB,KLB+1,-1
            if (Depth1D (k-1)<Depth1D (k)) exit
        enddo
        
        kb = k
        
        Ndepths       = KUB - kb + 1

        ValueAtDepthZ = InterpolateProfileR8(dble(Z), Ndepths, Depth1D (kb:KUB), Matrix1D(kb:KUB))
        
    end function ValueAtDepthZ

    !--------------------------------------------------------------------------        
    
    subroutine CheckAndCorrectVarName(obj_name, Name)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: obj_name

        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: Name
        integer                                     :: i

        !Begin-----------------------------------------------------------------

        if(scan(obj_name, "_") .ne. 0 .and. scan(obj_name, "0") .ne. 0)then
            Name = obj_name(1:len_trim(obj_name)-6)
        else
            Name = trim(obj_name)
        endif

        do i = 1, len_trim(Name)
            if (Name(i:i) == ' ') then
                Name(i:i) =  '_'
            endif
        enddo

    end subroutine CheckAndCorrectVarName
    
    !--------------------------------------------------------------------------

    subroutine KillConvert2netcdf

        call StopCPUTime

        call ShutdownMohid ("Convert2netcdf", Me%ElapsedSeconds, Me%TotalCPUTime)

    end subroutine KillConvert2netcdf
    
    !--------------------------------------------------------------------------

    subroutine StartCPUTime

        call date_and_time(Values = Me%F95Time)

        call SetDate      (Me%InitialSystemTime, float(Me%F95Time(1)), float(Me%F95Time(2)),      &
                                                 float(Me%F95Time(3)), float(Me%F95Time(5)),      &
                                                 float(Me%F95Time(6)), float(Me%F95Time(7))+      &
                                                 Me%F95Time(8)/1000.)

    end subroutine StartCPUTime

    !--------------------------------------------------------------------------

    subroutine StopCPUTime

        call date_and_time(Values = Me%F95Time)

        call SetDate      (Me%FinalSystemTime,   float(Me%F95Time(1)), float(Me%F95Time(2)),      &
                                                 float(Me%F95Time(3)), float(Me%F95Time(5)),      &
                                                 float(Me%F95Time(6)), float(Me%F95Time(7))+      &
                                                 Me%F95Time(8)/1000.)

        call cpu_time(Me%TotalCPUTime)

        Me%ElapsedSeconds = Me%FinalSystemTime - Me%InitialSystemTime

    end subroutine StopCPUTime
    
    !--------------------------------------------------------------------------
    
end program Convert2netcdf

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior T�cnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------