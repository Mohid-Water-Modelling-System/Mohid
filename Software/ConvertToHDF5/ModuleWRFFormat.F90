!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : WRFFormat
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2008
! REVISION      : Rosa Trancoso
! DESCRIPTION   : Module to convert WRF netcdf files into HDF5 format.
!
!------------------------------------------------------------------------------
!DataFile
!
!   FILENAME                    : char              -           !Path to WRF original file
!   OUTPUT_GRID_FILENAME        : char              -           !Path to grid data file generated from WRF file
!   OUTPUTFILENAME              : char              -           !Path to HDF5 file generated from WRF file
!   WRITE_XYZ                   : 0/1              [0]          !Flag to write xyz center grid cells
!   WRITE_TERRAIN               : 0/1              [0]          !Flag to write WRF TERRAIN fields
!   WRITE_HEADERS               : 0/1              [0]          !Flag to output WRF ans TERRAIN headers in final HDF5 file
!   COMPUTE_WINDSTRESS          : 0/1              [0]          !Flag to compute and write wind shear stress fields
!   COMPUTE_RELATIVE_HUMIDITY   : 0/1              [0]          !Flag to compute and write relative humidity fields at 2-m
!   COMPUTE_RELATIVE_HUMIDITY_3D: 0/1              [0]          !Flag to compute and write relative humidity 3D fields
!   COMPUTE_PRECIPITATION       : 0/1              [0]          !Flag to compute and write precipitation field

!
!   START                       : YYYY MM DD HH MM SS   [-]     !Start date of new file
!   END                         : YYYY MM DD HH MM SS   [-]     !End date of new file


!   <<BeginFields>>
!   air temperature
!   atmospheric pressure
!   wind velocity X
!   ...
!   ... (see below for available fields)
!   <<EndFields>>


!Some notes on how WRF file is organized and read
!dimensions:
!	Time = UNLIMITED ; // (25 currently)
!	DateStrLen       = 19 ;
!	west_east        = 39 ;
!	south_north      = 54 ;
!	bottom_top       = 27 ;
!	bottom_top_stag  = 28 ;
!	soil_layers_stag = 5 ;
!	west_east_stag   = 40 ;
!	south_north_stag = 55 ;

!	float U(Time, bottom_top, south_north, west_east_stag) ;
!		U:FieldType = 104 ;
!		U:MemoryOrder = "XYZ" ;
!		U:description = "x-wind component" ;
!		U:units = "m s-1" ;
!		U:stagger = "X" ;
!		U:coordinates = "XLONG_U XLAT_U" ;

Module ModuleWRFFormat

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalGrid    
    use netcdf
    use proj4


    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertWRFFormat
    private ::      ReadOptions
    private ::          ReadFieldsToConvert
    private ::          AtLeastOne3DFieldIsToConvert
    private ::          FieldIsToConvert
    private ::      Open_HDF5_OutPut_File
    private ::      OpenAndReadWRFFile
    private ::          FillDimInformation
    private ::          ReadWriteGridData
    private ::              InitializeProjection
    private ::              CorrectCorners
    private ::              ConstructGrid 
    private ::              WriteCenterGridXYZ
    private ::              WriteGridToHDF5File
    private ::          SetNewDate    
    private ::          VariableIsToRead
    private ::          FieldIsToRead
    private ::          AddField
    private ::          SetNewFieldAttributes
    private ::      BeginEndTime
    private ::      ComputeVerticalCoordinate        
    private ::      ComputeAtmosphericPressure3D
    private ::      ComputeAirTemperature3D
    private ::      ComputeMeanSeaLevelPressureMM5
    private ::      ComputeMeanSeaLevelPressureWRF
    

    !copy from ModuleMM5Format - if change here, please change there. This points to a ModuleMeteoFormat!

    private ::      ComputeRelativeHumidity
    private ::      ComputeRelativeHumidity3D
    private ::      ComputePrecipitation
    private ::      ComputeWindShearStress
    private ::      ComputeWindModulus

    private ::      OutputFields
    private ::      KillWRFFormat


    !Parameters----------------------------------------------------------------
    integer, parameter                                      :: DotGrid              = 1          !Cell corners
    integer, parameter                                      :: CrossGrid            = 2          !Cell center
    integer, parameter                                      :: FullLevels           = 1          !layer face
    integer, parameter                                      :: HalfLevels           = 2          !layer center

    !Perfect gas constant in Pa m3 kg-1 K-1
    real, parameter                                      :: R_dry    = 287.04       !J kg-1 K-1
    real, parameter                                      :: Cp_dry   = 7.*R_dry/2.  !1004.64 J kg-1 K-1
    real, parameter                                      :: p1000mb  = 100000.      !Pa
    real, parameter                                      :: GAMMA    = 6.5E-3       !K/m (neutral stability)

    character(StringLength), parameter  :: AccConvPrecipitation     = "accumulated convective precipitation"
    character(StringLength), parameter  :: AccNonConvPrecipitation  = "accumulated non-convective precipitation"


    !Types---------------------------------------------------------------------

      
    type       T_Dim
        integer                                             :: ID
        character(len=256)                                  :: Name
        integer                                             :: Size
    end type  T_Dim


    type       T_Date
        type(T_Time)                                        :: Date
        type(T_Date), pointer                               :: Next
    end type  T_Date

    type       T_Field
        character(len=StringLength)                         :: Name
        character(len=StringLength)                         :: Units
        logical                                             :: Convert                  = .false.
        integer                                             :: nDimensions
        logical                                             :: StaggeredX               = .false.
        logical                                             :: StaggeredY               = .false.
        logical                                             :: StaggeredZ               = .false.        
        type(T_Time)                                        :: Date
        real, dimension(:,:  ),     pointer                 :: Values2D
        real, dimension(:,:,:),     pointer                 :: Values3D
        integer                                             :: OutputNumber             = 1        
        type(T_Size3D)                                      :: Size, WorkSize
        type(T_Field),              pointer                 :: Next
    end type  T_Field

 
    type       T_WRFFormat
        integer                                             :: ObjEnterData             = 0
        integer                                             :: ObjHDF5                  = 0
        integer                                             :: ObjHorizontalGrid        = 0
        integer                                             :: Unit
        character(len=PathLength)                           :: FileName
        character(len=PathLength)                           :: TerrainFileName
        character(len=PathLength)                           :: GridFileName
        character(len=PathLength)                           :: OutputFileName
        real, dimension(:,:  ),       pointer               :: Bathymetry
        real, dimension(:,:  ),       pointer               :: CenterX
        real, dimension(:,:  ),       pointer               :: CenterY
        real, dimension(:,:  ),       pointer               :: ConnectionX
        real, dimension(:,:  ),       pointer               :: ConnectionY
        real, dimension(:,:  ),       pointer               :: LandUse
        real, dimension(:    ),       pointer               :: Sigma
        logical                                             :: WriteXYZ                 = .false.
        logical                                             :: ComputeWindStress        = .false.
        logical                                             :: ComputeWindModulus       = .false.
        logical                                             :: ComputeRelativeHumidity  = .false.
        logical                                             :: ComputeRelativeHumidity3D= .false.
        logical                                             :: ComputePrecipitation     = .false.
        logical                                             :: ComputeMeanSeaLevelPressureMM5 = .false.
        logical                                             :: ComputeMeanSeaLevelPressureWRF = .false.
        logical                                             :: ComputePressure3D        = .false.
        logical                                             :: ComputeTemperature3D     = .false.

        logical                                             :: WriteTerrain             = .false.
        logical                                             :: TerrainFileNotSpecified  = .false.
        logical                                             :: OutputGridNotSpecified   = .false.
        real                                                :: PTop                     = null_real
        integer                                             :: IfSnow                   = null_int
        type(T_Size3D)                                      :: Size, WorkSize
        type(T_Field),                pointer               :: FirstField
        type(T_Date),                 pointer               :: FirstDate
        character(len=StringLength), dimension(:), pointer  :: FieldsToConvert
        character(len=StringLength), dimension(:), pointer  :: FieldsToRead
        logical                                             :: TimeWindow
        type(T_Time)                                        :: StartTime, EndTime
        real                                                :: OutputDTInterval = null_real
        
        integer                                             :: ProjType                 = null_int
        real                                                :: CoarseDomainCenterLat    = null_real
        real                                                :: CoarseDomainCenterLon    = null_real
        real                                                :: TrueLatUpper             = null_real
        real                                                :: TrueLatLower             = null_real
        real                                                :: DY                       = null_real
      
        type(prj90_projection)                              :: Proj
        character(len=20), dimension(8)                     :: Params


    end type  T_WRFFormat

    type(T_WRFFormat), pointer                              :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertWRFFormat(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        call ReadOptions(ClientNumber)

        call Open_HDF5_OutPut_File

        call OpenAndReadWRFFile

        call BeginEndTime

        if(AtLeastOne3DFieldIsToConvert())      call ComputeVerticalCoordinate        

        if(Me%ComputePressure3D)                call ComputeAtmosphericPressure3D        

        if(Me%ComputeTemperature3D)             call ComputeAirTemperature3D        

        if(Me%ComputeMeanSeaLevelPressureMM5)   call ComputeMeanSeaLevelPressureMM5
        if(Me%ComputeMeanSeaLevelPressureWRF)   call ComputeMeanSeaLevelPressureWRF

        !like in ModuleMM5Format

        if(Me%ComputeWindStress)                call ComputeWindShearStress

        if (Me%ComputeWindModulus)              call ComputeWindModulus

        if(Me%ComputeRelativeHumidity)          call ComputeRelativeHumidity

        if(Me%ComputeRelativeHumidity3D)        call ComputeRelativeHumidity3D

        if(Me%ComputePrecipitation)             call ComputePrecipitation

        call OutputFields

        call KillWRFFormat


        STAT = SUCCESS_


    end subroutine ConvertWRFFormat

    !------------------------------------------------------------------------

    subroutine ReadOptions(ClientNumber)
        
        !Arguments-------------------------------------------------------------
        integer                             :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag, iflag1

        !Begin-----------------------------------------------------------------
       
        call GetData(Me%FileName,                                       &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'FILENAME',                         &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR01'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleWRFFormat - ERR02'
        end if

        call GetData(Me%GridFileName,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',             &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR03'

        if (iflag == 0)then
            Me%OutputGridNotSpecified = .true.
            write(*,*)
            write(*,*)'Output grid file (Mohid format) was not specified.'
            write(*,*)'    - Mohid grid will not be constructed.'
            write(*,*)'    - HDF5 file will not be usable to interpolate'
            write(*,*)'      to other grids.'
            write(*,*)
        end if


        call GetData(Me%OutputFileName,                                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUTFILENAME',                   &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR04'


        call GetData(Me%WriteXYZ,                                       &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'WRITE_XYZ',                        &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR06'

        call GetData(Me%ComputeWindStress,                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'COMPUTE_WINDSTRESS',               &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR07'
        

        call GetData(Me%ComputeWindModulus,                             &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'COMPUTE_WINDMODULUS',              &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR07a'


        call GetData(Me%ComputeRelativeHumidity,                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'COMPUTE_RELATIVE_HUMIDITY',        &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR08'

        call GetData(Me%ComputeRelativeHumidity3D,                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'COMPUTE_RELATIVE_HUMIDITY_3D',     &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR08'


        call GetData(Me%ComputePrecipitation,                           &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'COMPUTE_PRECIPITATION',            &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR08a'

        call GetData(Me%ComputeMeanSeaLevelPressureMM5,                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'COMPUTE_MSLP_MM5',                 &
                     Default      = ON,                                 &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR08b'

        call GetData(Me%ComputeMeanSeaLevelPressureWRF,                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'COMPUTE_MSLP_WRF',                 &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR08b1'

        call GetData(Me%WriteTerrain,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'WRITE_TERRAIN',                    &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR09'
        

        call GetData(Me%StartTime,                                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'START',                            &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR10'


        call GetData(Me%EndTime,                                        &
                     Me%ObjEnterData, iflag1,                           &
                     SearchType   = FromBlock,                          &
                     keyword      = 'END',                              &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR11'

        if (iflag==1 .AND. iflag1==1) Me%TimeWindow = .TRUE.
        
        if (Me%TimeWindow) then

            if (Me%StartTime .GE. Me%EndTime) then
                write (*,*)
                write (*,*) 'START greater or equal than END'
                stop 'ReadOptions - ModuleWRFFormat - ERR12'
            endif

        endif

        call GetData(Me%OutputDTInterval,                               &
                     Me%ObjEnterData, iflag1,                           &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUT_DT',                        & 
                     Default      = 0.0,                                &
                     ClientModule = 'ModuleWRFFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWRFFormat - ERR11'

        call ReadFieldsToConvert(ClientNumber)

    end subroutine ReadOptions


    !--------------------------------------------------------------------------


    subroutine ReadFieldsToConvert(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer                             :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL
        integer                             :: StartLine, EndLine, Count
        integer                             :: CurrentLineNumber
        integer                             :: nPropertiesToConvert
        logical                             :: BlockFound
        character(len=StringLength)         :: PropertyName

        logical                             :: NeedsPressureBaseState           = .false.
        logical                             :: NeedsPressurePerturbation        = .false.
        logical                             :: NeedsGeopotenticalBaseState      = .false.
        logical                             :: NeedsGeopotenticalPerturbation   = .false.

        logical                             :: NeedsTemperature3D               = .false.
        logical                             :: NeedsPressure3D                  = .false.
        logical                             :: NeedsPotTemp3D                   = .false.
        logical                             :: NeedsMixingRatio3D               = .false.

        logical                             :: NeedsSurfaceTemperature          = .false.
        logical                             :: NeedsSurfacePressure             = .false.
        logical                             :: NeedsSurfaceMixingRatio          = .false.
        logical                             :: NeedsWindVelocityX               = .false.
        logical                             :: NeedsWindVelocityY               = .false.
        logical                             :: NeedsWindShearVelocity           = .false.
        logical                             :: NeedsConvectivePCP               = .false.
        logical                             :: NeedsNonConvectivePCP            = .false.
        logical                             :: NeedsTerrain                     = .false.

        !Begin-----------------------------------------------------------------

        call ExtractBlockFromBlock(Me%ObjEnterData,             &
                                    ClientNumber,               &
                                    '<<BeginFields>>',          &
                                    '<<EndFields>>',            &
                                    BlockFound,                 &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModuleWRFFormat - ERR01'

        call GetBlockSize(Me%ObjEnterData,                      &
                          ClientNumber,                         &
                          StartLine,                            &
                          EndLine,                              &
                          FromBlockInBlock,                     &
                          STAT = STAT_CALL)

        nPropertiesToConvert = EndLine - StartLine - 1

        allocate(Me%FieldsToConvert(1:nPropertiesToConvert))

        Count = 1

        do CurrentLineNumber = StartLine + 1 , EndLine - 1

            call GetData(PropertyName,                          &
                         Me%ObjEnterData,                       &
                         flag,                                  &
                         SearchType  = FromBlock_,              &
                         Buffer_Line = CurrentLineNumber,       &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModuleWRFFormat - ERR02'

            Me%FieldsToConvert(Count) = PropertyName

            Count = Count + 1

        end do

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModuleWRFFormat - ERR03'

        allocate(Me%FieldsToRead(1:100))

        Me%FieldsToRead(:) = null_str

        Me%FieldsToRead(1:nPropertiesToConvert) = Me%FieldsToConvert(1:nPropertiesToConvert)

        Count = nPropertiesToConvert + 1

!         if(FieldIsToConvert(trim(GetPropertyName(AtmosphericPressure_)))) then
!             NeedsPressure3D             = .true.
!         end if

        if(FieldIsToConvert(trim(GetPropertyName(AtmosphericPressure_))//"_3D"))then
            NeedsPressure3D             = .true.            
        end if

        if(FieldIsToConvert(trim(GetPropertyName(AirTemperature_))//"_3D"))then
            NeedsTemperature3D          = .true.
        end if

        if(Me%ComputeRelativeHumidity)then
            NeedsPressure3D             = .true.
            NeedsSurfaceTemperature     = .true.
            NeedsSurfaceMixingRatio     = .true.
        end if

        if(Me%ComputeRelativeHumidity3D)then
            NeedsPressure3D             = .true.
            NeedsTemperature3D          = .true.
            NeedsMixingRatio3D          = .true.
        end if

        if(Me%ComputePrecipitation)then
            NeedsConvectivePCP          = .true.
            NeedsNonConvectivePCP       = .true.
        end if

        if(Me%ComputeMeanSeaLevelPressureMM5)then
            NeedsSurfacePressure        = .true.
            NeedsPressure3D             = .true.
            NeedsTemperature3D          = .true.
            NeedsTerrain                = .true.
        end if

        if(Me%ComputeMeanSeaLevelPressureWRF)then
            NeedsPressure3D             = .true.
            NeedsPotTemp3D              = .true.
            NeedsMixingRatio3D          = .true.
        end if


        if(Me%ComputeWindStress)then
            NeedsWindVelocityX          = .true.
            NeedsWindVelocityY          = .true.
            NeedsWindShearVelocity      = .true. 
        end if

        if (Me%ComputeWindModulus) then 
            NeedsWindVelocityX          = .true.
            NeedsWindVelocityY          = .true.
        endif

        if(AtLeastOne3DFieldIsToConvert())then
            !Needs VerticalZ
            NeedsGeopotenticalPerturbation  = .true.
            NeedsGeopotenticalBaseState     = .true.
        end if


        if (NeedsTemperature3D) then
            Me%ComputeTemperature3D     = .true.
            NeedsPotTemp3D              = .true.
            NeedsPressure3D             = .true.  
        endif

        if (NeedsPressure3D) then
            Me%ComputePressure3D        = .true.
            NeedsPressureBaseState      = .true.
            NeedsPressurePerturbation   = .true.
        endif


        !allocate Me%FieldsToRead ------------------------------------------------

        if (NeedsPressureBaseState .and. .not. FieldIsToConvert(trim('base state pressure_3D'))) then
            Me%FieldsToRead(Count) = trim('base state pressure_3D')
            Count = Count + 1
        end if

        if(NeedsPressurePerturbation .and. .not. FieldIsToConvert(trim(GetPropertyName(AtmosphericPressure_))//"_3D"))then
            Me%FieldsToRead(Count) = trim(GetPropertyName(AtmosphericPressure_))//"_3D"
            Count = Count + 1
        end if

        if(NeedsGeopotenticalBaseState .and. .not. FieldIsToConvert(trim('geopotential base state_3D'))) then
            Me%FieldsToRead(Count) = trim('geopotential base state_3D')
            Count = Count + 1
        end if

       if(NeedsGeopotenticalPerturbation .and. .not. FieldIsToConvert(trim('geopotential perturbation_3D'))) then
            Me%FieldsToRead(Count) = trim('geopotential perturbation_3D')
            Count = Count + 1
        end if

        if(NeedsTemperature3D .and. .not. FieldIsToConvert(trim(GetPropertyName(AirTemperature_))//"_3D"))then
            Me%FieldsToRead(Count) = trim(GetPropertyName(AirTemperature_))//"_3D"
            Count = Count + 1
        end if

        if(NeedsPotTemp3D .and. .not. FieldIsToConvert(trim("potential temperature_3D")))then
            Me%FieldsToRead(Count) = trim("potential temperature_3D")
            Count = Count + 1
        end if

        if(NeedsMixingRatio3D .and. .not. FieldIsToConvert(trim("mixing ratio_3D")))then
            Me%FieldsToRead(Count) = trim("mixing ratio_3D")
            Count = Count + 1
        end if

        !2D -------------------------------------------------------------------

        if(NeedsSurfaceTemperature .and. .not. FieldIsToConvert(trim(GetPropertyName(AirTemperature_))))then
            Me%FieldsToRead(Count) = trim(GetPropertyName(AirTemperature_))
            Count = Count + 1
        end if

        if(NeedsSurfacePressure .and. .not. FieldIsToConvert(trim(GetPropertyName(AtmosphericPressure_))))then
            Me%FieldsToRead(Count) = trim(GetPropertyName(AtmosphericPressure_))
            Count = Count + 1
        end if

        if(NeedsSurfaceMixingRatio .and. .not. FieldIsToConvert(trim("2-meter mixing ratio")))then
            Me%FieldsToRead(Count) = trim("2-meter mixing ratio")
            Count = Count + 1
        end if


        if(NeedsWindVelocityX .and. .not. FieldIsToConvert(trim(GetPropertyName(WindVelocityX_))))then
            Me%FieldsToRead(Count) = trim(GetPropertyName(WindVelocityX_))
            Count = Count + 1
        end if

        if(NeedsWindVelocityY .and. .not. FieldIsToConvert(trim(GetPropertyName(WindVelocityY_))))then
            Me%FieldsToRead(Count) = trim(GetPropertyName(WindVelocityY_))
            Count = Count + 1
        end if


        if(NeedsWindShearVelocity .and. .not. FieldIsToConvert(trim(GetPropertyName(WindShearVelocity_))))then
            Me%FieldsToRead(Count) = trim(GetPropertyName(WindShearVelocity_))
            Count = Count + 1
        end if

        if(NeedsConvectivePCP .and. .not. FieldIsToConvert(trim(AccConvPrecipitation)))then
            Me%FieldsToRead(Count) = trim(AccConvPrecipitation)
            Count = Count + 1 
        endif

        if(NeedsNonConvectivePCP .and. .not. FieldIsToConvert(trim(AccNonConvPrecipitation)))then
            Me%FieldsToRead(Count) = trim(AccNonConvPrecipitation)
            Count = Count + 1 
        endif

        if(NeedsTerrain .and. .not. FieldIsToConvert(trim('terrain')))then
            Me%FieldsToRead(Count) = trim('terrain')
            Count = Count + 1 
        endif

    end subroutine ReadFieldsToConvert


    !--------------------------------------------------------------------------

    logical function FieldIsToConvert(FieldName)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: FieldName

        !Local----------------------------------------------------------------
        integer                                     :: i

        !Begin----------------------------------------------------------------
        
        FieldIsToConvert = .false.

        do i = 1, size(Me%FieldsToConvert)

            if(Me%FieldsToConvert(i) == trim(FieldName))then
                FieldIsToConvert = .true.
            end if

        enddo

    end function FieldIsToConvert 

    !--------------------------------------------------------------------------
    
    logical function AtLeastOne3DFieldIsToConvert()

        !Local----------------------------------------------------------------
        integer                                     :: i

        !Begin----------------------------------------------------------------
        
        AtLeastOne3DFieldIsToConvert = .false.

        do i = 1, size(Me%FieldsToConvert)

            if(scan(Me%FieldsToConvert(i), "_3D") .ne. 0)then
                AtLeastOne3DFieldIsToConvert = .true.
                return
            end if

        enddo

    end function AtLeastOne3DFieldIsToConvert 

    !----------------------------------------------------------------------

    subroutine FillDimInformation(DimName, Size)
        
        !Arguments-------------------------------------------------------------
        character (len= *)                          :: DimName
        integer                                     :: Size

        select case (DimName)

            case ('west_east')                 

                Me%Size%JLB = 0
                Me%Size%JUB = Size + 1
                
                Me%WorkSize%JLB = Me%Size%JLB + 1
                Me%WorkSize%JUB = Size

            case ('south_north')

                Me%Size%ILB = 0
                Me%Size%IUB = Size + 1
                
                Me%WorkSize%ILB = Me%Size%ILB + 1
                Me%WorkSize%IUB = Size

            case ('bottom_top')

                Me%Size%KLB = 0           ; Me%WorkSize%KLB = Me%Size%KLB + 1
                Me%Size%KUB = size + 1
                
                Me%WorkSize%KLB = Me%Size%KLB + 1
                Me%WorkSize%KUB = Size

            case default

        end select

    end subroutine FillDimInformation

   
    !--------------------------------------------------------------------------
    
    subroutine OpenAndReadWRFFile

        !Local-----------------------------------------------------------------
        logical                                         :: exist
        integer                                         :: status
        integer                                         :: ncid, ndims, nvars, NGlobalAtts, unlimitedDimId
        type(T_Dim), pointer, dimension(:)              :: Dims
        integer                                         :: i, j, k, it, len
        character(Len=StringLength)                     :: name
        integer                                         :: varid, xtype, nDimensions, natts
        integer, pointer, dimension(:)                  :: VarDimIds
        type(T_Date), pointer                           :: CurrentDate
        character(1)    , pointer, dimension(:,:,:,:)   :: DataAuxChar
        real            , pointer, dimension(:,:,:,:)   :: DataAux
        character(len=19)                               :: current_date
        character(Len=StringLength)                     :: MohidName

        type(T_Field), pointer                          :: NewField
        character(Len=StringLength)                     :: Units
        character(1)                                    :: stagger
        
        integer                                         :: ILB , IUB , JLB , JUB , KLB , KUB
        integer                                         :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        integer                                         :: NX,NY,NZ, DateStrLen, TimeSize
        real                                            :: aux1, aux2
       
        

        !Begin-----------------------------------------------------------------
        
        write(*,*)'---------------------------'
        write(*,*)
        write(*,*)'Reading WRF output file...'

        nullify(Me%FirstDate  )

        !Verifies if file exists
        inquire(file = Me%FileName, exist = exist)
        if (.not. exist) then
            write(*,*)'WRF output file does not exist'
            stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR01'
        endif

        status = nf90_open(trim(Me%FileName),NF90_NOWRITE,ncid)
        call handle_error(status); if (status /= nf90_noerr) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR01'

        status = nf90_inquire(ncid, ndims, nvars, NGlobalAtts, unlimitedDimId)
        call handle_error(status); if (status /= NF90_NOERR) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR02'

        write(*,*)'ndims =', ndims
        write(*,*)'nvars =', nvars 
        write(*,*)'natts =', NGlobalAtts
        write(*,*)'unlimitedDimId = ', unlimitedDimId

        ! Size ----------------------------------------------------------------

        allocate(Dims(ndims))
        do i=1,ndims
    
            status = nf90_inquire_dimension(ncid, i, name, len)
            call handle_error(status); if (status /= NF90_NOERR) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR03'

            write(*,'(1x,i5,a20,i5)') i, trim(name), len
            Dims(i)%ID   = i
            Dims(i)%Name = name
            Dims(i)%Size = len
        
            call FillDimInformation(name, len)

        enddo

        ! GEOMETRY - time invariant  ------------------------------------------

        write(*,*) 'Reading Geometry...'
        call ReadWriteGridData (ncid, nvars)

        ! Time ----------------------------------------------------------------

        write(*,*) 'Reading Times...'

        status = nf90_inq_varid(ncid, 'Times', varid)
        call handle_error(status); if (status /= NF90_NOERR) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR04'

        status = nf90_inquire_variable(ncid, varid, name, xtype, ndims=nDimensions, natts=natts)
        call handle_error(status); if (status /= NF90_NOERR) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR05'

        allocate(VarDimIds(nDimensions))

        status = nf90_inquire_variable(ncid, varid, name, dimids=VarDimIds)
        call handle_error(status); if (status /= NF90_NOERR) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR06'
        
        DateStrLen = Dims(VarDimIds(1))%Size !DateStrLen = 19
        TimeSize   = Dims(VarDimIds(2))%Size !Time

        allocate(DataAuxChar(DateStrLen,TimeSize, 1,1))

        status = nf90_get_var(ncid, varid, DataAuxChar, start=(/1,1,1/), count=(/DateStrLen, TimeSize, 1/))
        call handle_error(status); if (status /= NF90_NOERR) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR07'
    
        do j=1,TimeSize
            current_date = ''
            do i=1,DateStrLen
                current_date = trim(adjustl(current_date))//trim(adjustl(DataAuxChar(i,j,1,1)))
            enddo
            call SetNewDate(current_date)
        enddo

        deallocate(DataAuxChar)
        deallocate(VarDimIds)
            
        !VARS - time dependent ------------------------------------------------

        write(*,*) 'Reading Variables...'

        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB
        KLB = Me%Size%KLB; WKLB = Me%WorkSize%KLB
        KUB = Me%Size%KUB; WKUB = Me%WorkSize%KUB

        nullify(Me%FirstField)

        CurrentDate => Me%FirstDate
        it = 1
doT:    do while(associated(CurrentDate))                

doV:        do varid = 1, nvars
            
            status = nf90_inquire_variable(ncid, varid, name, xtype, ndims=nDimensions, natts=natts)
            call handle_error(status); if (status /= NF90_NOERR) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR08'

if0:        if(VariableIsToRead(name, MohidName)) then

                allocate(VarDimIds(nDimensions))

                status = nf90_inquire_variable(ncid, varid, name, dimids=VarDimIds)
                call handle_error(status); if (status /= NF90_NOERR) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR09'

                call AddField(Me%FirstField, NewField)

                NewField%nDimensions    = nDimensions - 1

                status = nf90_get_att(ncid, varid, 'units', Units)
                call handle_error(status); if (status /= NF90_NOERR) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR10'                

                call SetNewFieldAttributes(Field    = NewField,                &
                                           Name     = trim(MohidName),         &
                                           Units    = trim(Units),             &
                                           Date     = CurrentDate%Date,        &
                                           WorkSize = Me%WorkSize,             &
                                           Convert  = FieldIsToConvert(MohidName))

                status = nf90_get_att(ncid, varid, 'stagger', stagger)
                call handle_error(status); if (status /= NF90_NOERR) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR11'

                if (trim(stagger) == 'X') NewField%StaggeredX = .TRUE.
                if (trim(stagger) == 'Y') NewField%StaggeredY = .TRUE.
                if (trim(stagger) == 'Z') NewField%StaggeredZ = .TRUE.


                select case(NewField%nDimensions)

                    case(2)

                        !Memory Order is always XY => JI

                        NX = Dims(VarDimIds(1))%Size
                        NY = Dims(VarDimIds(2))%Size
                        allocate(DataAux(NX,NY, 1,1))

                        status = nf90_get_var(ncid, varid, DataAux, start=(/1,1,it/), count=(/NX, NY, 1/))
                        call handle_error(status); if (status /= NF90_NOERR) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR12'

                        do j=1,NY
                        do i=1,NX

                            if (isnan(DataAux(i,j, 1, 1))) then
                                write(*,*) 'NaN values in WRF property ', trim(name)
                                write(*,*) 'at (i,j)    = ', i,j
                                stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR13'
                            endif 

                       enddo
                       enddo


                        allocate(NewField%Values2D(ILB:IUB, JLB:JUB))

                        if (NewField%StaggeredX .and. NewField%StaggeredY) then !!!!!!!!!!isto nao existe!!!!!!!!!!!

                            do j = WJLB,WJUB
                            do i = WILB,WIUB
                                NewField%Values2D(i, j) =   &
                                (DataAux(j  , i  , 1, 1) +  &
                                 DataAux(j+1, i+1, 1, 1) +  &
                                 DataAux(j+1, i  , 1, 1) +  &
                                 DataAux(j  , i+1, 1, 1))/4.
                            enddo
                            enddo

!                             NewField%Values2D(WILB:WIUB, WJLB:WJUB) =               &
!                             (DataAux(WILB    :WIUB    , WJLB    :WJUB    , 1, 1) +  &
!                              DataAux(WILB + 1:WIUB + 1, WJLB + 1:WJUB + 1, 1, 1) +  &
!                              DataAux(WILB + 1:WIUB + 1, WJLB    :WJUB    , 1, 1) +  &
!                              DataAux(WILB    :WIUB    , WJLB + 1:WJUB + 1, 1, 1))/ 4.

                        else if (NewField%StaggeredX) then

                            do j = WJLB,WJUB
                            do i = WILB,WIUB
                                NewField%Values2D(i, j) =   &
                                (DataAux(j  , i, 1, 1) +    &
                                 DataAux(j+1, i, 1, 1))/2.
                            enddo
                            enddo

!                             NewField%Values2D(WILB:WIUB, WJLB:WJUB) =               &
!                             (DataAux(WILB    :WIUB    , WJLB    :WJUB    , 1, 1) +  &
!                              DataAux(WILB + 1:WIUB + 1, WJLB    :WJUB,     1, 1))/2

                        elseif (NewField%StaggeredY) then

                            do j = WJLB,WJUB
                            do i = WILB,WIUB
                                NewField%Values2D(i, j) =   &
                                (DataAux(j, i  , 1, 1) +    &
                                 DataAux(j, i+1, 1, 1))/2.
                            enddo
                            enddo

!                             NewField%Values2D(WILB:WIUB, WJLB:WJUB) =               &
!                             (DataAux(WILB    :WIUB    , WJLB    :WJUB    , 1, 1) +  &
!                              DataAux(WILB    :WIUB    , WJLB + 1:WJUB + 1, 1, 1))/2

                        else

                            do j = WJLB,WJUB
                            do i = WILB,WIUB
                                NewField%Values2D(i, j) = DataAux(j, i, 1,1)
                            enddo
                            enddo

!                            NewField%Values2D(WILB:WIUB,WJLB:WJUB) = DataAux(WILB:WIUB,WJLB:WJUB, 1, 1)

                        end if

                    case(3) !3d + time

                        !Memory Order is always XYZ => JIK
                        NX = Dims(VarDimIds(1))%Size
                        NY = Dims(VarDimIds(2))%Size
                        NZ = Dims(VarDimIds(3))%Size
                        allocate(DataAux(NX,NY,NZ,1))

                        status = nf90_get_var(ncid, varid, DataAux, start=(/1,1,1,it/), count=(/NX,NY,NZ, 1/))
                        call handle_error(status); if (status /= NF90_NOERR) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR14'

                        do k=1,NZ
                        do j=1,NY
                        do i=1,NX

                            if (isnan(DataAux(i,j,k,1))) then
                                write(*,*) 'NaN values in WRF property ', trim(name)
                                write(*,*) 'at (i,j,k)    = ', i,j,k
                                stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR15'
                            endif 

                        enddo
                        enddo
                        enddo
            
                        allocate(NewField%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))

                        if (NewField%StaggeredX .and. NewField%StaggeredY) then

                            if (NewField%StaggeredZ) then

                                do j = WJLB, WJUB
                                do i = WILB, WIUB
                                do k = WKLB, WKUB

                                    aux1 =  (DataAux(j  , i  , k, 1) + DataAux(j+1, i+1, k, 1) +        &
                                             DataAux(j+1, i  , k, 1) + DataAux(j  , i+1, k, 1))/4.

                                    aux1 =  (DataAux(j  , i  , k+1, 1) + DataAux(j+1, i+1, k+1, 1) +    &
                                             DataAux(j+1, i  , k+1, 1) + DataAux(j  , i+1, k+1, 1))/4.

                                    NewField%Values3D(i, j, k) = (aux1 + aux2) / 2.

                                enddo
                                enddo
                                enddo

                            else 

                                do j = WJLB, WJUB
                                do i = WILB, WIUB
                                do k = WKLB, WKUB

                                    NewField%Values3D(i, j, k) =                            &
                                    (DataAux(j  , i  , k, 1) + DataAux(j+1, i+1, k, 1) +    &
                                     DataAux(j+1, i  , k, 1) + DataAux(j  , i+1, k, 1))/4.

                                enddo
                                enddo
                                enddo

                            endif

                        else if (NewField%StaggeredX) then


                            if (NewField%StaggeredZ) then

                                do j = WJLB, WJUB
                                do i = WILB, WIUB
                                do k = WKLB, WKUB

                                    aux1 = (DataAux(j  , i, k  , 1) + DataAux(j+1, i, k  , 1))/2.
                                    aux2 = (DataAux(j  , i, k+1, 1) + DataAux(j+1, i, k+1, 1))/2.

                                    NewField%Values3D(i, j, k) = (aux1 + aux2) / 2.

                                enddo
                                enddo
                                enddo

                            else 

                                do j = WJLB, WJUB
                                do i = WILB, WIUB
                                do k = WKLB, WKUB

                                    NewField%Values3D(i, j, k) =                                              &
                                    (DataAux(j,   i, k, 1) + DataAux(j+1, i, k, 1))/2

                                enddo
                                enddo
                                enddo

                            endif


                        else if (NewField%StaggeredY) then


                            if (NewField%StaggeredZ) then

                                do j = WJLB, WJUB
                                do i = WILB, WIUB
                                do k = WKLB, WKUB

                                    aux1 = (DataAux(j, i, k  , 1) + DataAux(j, i+1, k  , 1))/2.
                                    aux2 = (DataAux(j, i, k+1, 1) + DataAux(j, i+1, k+1, 1))/2.

                                    NewField%Values3D(i, j, k) = (aux1 + aux2) / 2.

                                enddo
                                enddo
                                enddo

                            else 

                                do j = WJLB, WJUB
                                do i = WILB, WIUB
                                do k = WKLB, WKUB

                                    NewField%Values3D(i, j, k) =                                              &
                                    (DataAux(j, i, k, 1) + DataAux(j, i+1, k, 1))/2.

                                enddo
                                enddo
                                enddo

                            endif

                        else

                            if (NewField%StaggeredZ) then

                                do j = WJLB, WJUB
                                do i = WILB, WIUB
                                do k = WKLB, WKUB

                                    NewField%Values3D(i, j, k) = (DataAux(j,i, k  , 1) + &
                                                                  DataAux(j,i, k+1, 1)) / 2.

                                enddo
                                enddo
                                enddo

                            else 

                                do j = WJLB, WJUB
                                do i = WILB, WIUB
                                do k = WKLB, WKUB

                                    NewField%Values3D(i, j, k) = DataAux(j,i, k, 1)

                                enddo
                                enddo
                                enddo

                            endif

                        end if

                    case default

                        write(*,*) 'Var ', trim(adjustl(name)), ' has ', nDimensions
                        stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR14'

                    end select                

                    deallocate(DataAux)
                    deallocate(VarDimIds)

                endif if0

            end do doV

            CurrentDate => CurrentDate%Next
            it = it + 1

        enddo doT
                    
100     write(*,*)'Finished reading WRF output file.'

        status = nf90_close(ncid)
        call handle_error(status); if (status /= NF90_NOERR) stop 'OpenAndReadWRFFile - ModuleWRFFormat - ERR99'

    end subroutine OpenAndReadWRFFile
    
    !----------------------------------------------------------------------

    subroutine BeginEndTime

        !Local-------------------------------------------------------------
        type(T_Date), pointer                   :: CurrentDate
        type(T_Time)                            :: StartTime, EndTime
        real                                    :: Year, Month, Day
        real                                    :: Hour, Minute, Second

100 format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)
      
        CurrentDate => Me%FirstDate
            
        StartTime = CurrentDate%Date

        do while(associated(CurrentDate))                

            EndTime = CurrentDate%Date
            CurrentDate => CurrentDate%Next

        end do 
        
if1:    if (Me%TimeWindow) then


            if (Me%StartTime .GT. EndTime .OR. Me%EndTime .LT. StartTime) then
                
                write (*,*) 
                write (*,*) 'Time Window not in file'
                write (*,*) 'File START at:'
                call ExtractDate(StartTime, Year, Month, Day, Hour, Minute, Second)        
                write (*,*) 'File END at:'
                call ExtractDate(EndTime, Year, Month, Day, Hour, Minute, Second)        
                write (*,fmt=100) Year, Month, Day, Hour, Minute, Second
                stop 'BeginEndTime - ModuleMM5Format - ERR01'                

            endif


            if (Me%StartTime .LT. StartTime) then
                
                write (*,*) 
                write (*,*) 'START is before the start time of file'
                write (*,*) 'File begins at:'
                call ExtractDate(StartTime, Year, Month, Day, Hour, Minute, Second)        
                write (*,fmt=100) Year, Month, Day, Hour, Minute, Second
                write (*,*) 'BeginEndTime - ModuleMM5Format - WARN01'                

            endif

            if (Me%EndTime .GT. EndTime) then                
                
                write (*,*) 
                write (*,*) 'END is after the end time of file'
                write (*,*) 'File ends at:'
                call ExtractDate(EndTime, Year, Month, Day, Hour, Minute, Second)        
                write (*,fmt=100) Year, Month, Day, Hour, Minute, Second
                write (*,*) 'BeginEndTime - ModuleMM5Format - WARN02'                
            
            endif
        
        else !if1

            Me%StartTime = StartTime
            Me%EndTime   = EndTime

        endif if1


    end subroutine BeginEndTime
    
    !--------------------------------------------------------------------------

    subroutine ReadWriteGridData (ncid, nvars)

        !Arguments-------------------------------------------------------------
        integer, intent(in)                             :: ncid, nvars      

        !Local-----------------------------------------------------------------
        integer                                         :: status, varid
        character(Len=StringLength)                     :: name
        integer                                         :: ILB , IUB , JLB , JUB , KLB , KUB
        integer                                         :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        real            , pointer, dimension(:,:,:,:)   :: DataAux
        logical                                         :: Bathymetry_OK  = .false.
        !logical                                         :: LandUse_OK     = .false.        
        logical                                         :: CenterX_OK     = .false.
        logical                                         :: CenterY_OK     = .false.
        logical                                         :: ConnectionX_OK = .false.
        logical                                         :: ConnectionY_OK = .false.
        


        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB
        KLB = Me%Size%KLB; WKLB = Me%WorkSize%KLB
        KUB = Me%Size%KUB; WKUB = Me%WorkSize%KUB


        call InitializeProjection(ncid)

        ! get HGT   -> Batymetry
        ! get XLAT  -> CenterY
        ! get XLONG -> CenterX
        ! get XLAT_U, XLONG_U -> passar para XY -> 
        !        diminuir o Y de DY/2 e aumentar o ultimo Y de DY/2 
        !        passar a Lat,Lon -> ConnectionX, ConnectionY.

        status = nf90_get_att(ncid, NF90_GLOBAL, 'DY', Me%DY)
        call handle_error(status); if (status /= NF90_NOERR) stop 'ReadWriteGridData - ModuleWRFFormat - ERR00' 

        do varid = 1, nvars

            status = nf90_inquire_variable(ncid, varid, name)
            call handle_error(status); if (status /= NF90_NOERR) stop 'ReadWriteGridData - ModuleWRFFormat - ERR01'

            select case (trim(name))

                case ('HGT')

                    !Memory Order = XY => JI
                    allocate(DataAux(WJUB,WIUB, 1,1))

                    status = nf90_get_var(ncid, varid, DataAux, start=(/1,1,1/), count=(/WJUB, WIUB, 1/))
                    call handle_error(status); if (status /= NF90_NOERR) stop 'ReadWriteGridData - ModuleWRFFormat - ERR02'

                    allocate(Me%Bathymetry(ILB:IUB, JLB:JUB))

                    Me%Bathymetry(WILB:WIUB, WJLB:WJUB) = transpose(DataAux(WJLB:WJUB, WILB:WIUB, 1, 1))

                    deallocate(DataAux)

                    Bathymetry_OK = .true.

                case ('XLONG')

                    !Memory Order = XY => JI
                    allocate(DataAux(WJUB,WIUB, 1,1))

                    status = nf90_get_var(ncid, varid, DataAux, start=(/1,1,1/), count=(/WJUB, WIUB, 1/))
                    call handle_error(status); if (status /= NF90_NOERR) stop 'ReadWriteGridData - ModuleWRFFormat - ERR02'

                    allocate(Me%CenterX(ILB:IUB, JLB:JUB))

                    Me%CenterX(WILB:WIUB, WJLB:WJUB) = transpose(DataAux(WJLB:WJUB, WILB:WIUB, 1, 1))

                    deallocate(DataAux)

                    CenterX_OK = .true.

                case ('XLAT')

                    !Memory Order = XY => JI
                    allocate(DataAux(WJUB,WIUB, 1,1))

                    status = nf90_get_var(ncid, varid, DataAux, start=(/1,1,1/), count=(/WJUB, WIUB, 1/))
                    call handle_error(status); if (status /= NF90_NOERR) stop 'ReadWriteGridData - ModuleWRFFormat - ERR02'

                    allocate(Me%CenterY(ILB:IUB, JLB:JUB))

                    Me%CenterY(WILB:WIUB, WJLB:WJUB) = transpose(DataAux(WJLB:WJUB, WILB:WIUB, 1, 1))

                    deallocate(DataAux)

                    CenterY_OK = .true.

                case ('XLONG_U')

                    !Memory Order = XY => JI
                    !Staggered = X

                    allocate(DataAux(JUB,WIUB, 1,1))

                    status = nf90_get_var(ncid, varid, DataAux, start=(/1,1,1/), count=(/JUB, WIUB, 1/))
                    call handle_error(status); if (status /= NF90_NOERR) stop 'ReadWriteGridData - ModuleWRFFormat - ERR02'

                    allocate(Me%ConnectionX(ILB:IUB, JLB:JUB))

                    Me%ConnectionX(WILB:WIUB, WJLB:WJUB+1) = transpose(DataAux(WJLB:WJUB+1, WILB:WIUB, 1, 1))

                    deallocate(DataAux)

                    ConnectionX_OK = .true.

                case ('XLAT_U')

                    !Memory Order = XY => JI
                    !Staggered = X

                    allocate(DataAux(JUB,WIUB, 1,1))

                    status = nf90_get_var(ncid, varid, DataAux, start=(/1,1,1/), count=(/JUB, WIUB, 1/))
                    call handle_error(status); if (status /= NF90_NOERR) stop 'ReadWriteGridData - ModuleWRFFormat - ERR02'

                    allocate(Me%ConnectionY(ILB:IUB, JLB:JUB))

                    Me%ConnectionY(WILB:WIUB, WJLB:WJUB+1) = transpose(DataAux(WJLB:WJUB+1, WILB:WIUB, 1, 1))

                    deallocate(DataAux)

                    ConnectionY_OK = .true.

                case default

            end select

            if (ConnectionX_OK .and. ConnectionY_OK) call CorrectCorners

            if (Bathymetry_OK  .and. & !LandUse_OK .and. &
                CenterX_OK     .and. CenterY_OK .and. &
                ConnectionX_OK .and. ConnectionY_OK) exit

        enddo
        
        if (.not. Bathymetry_OK)    stop 'Missing HGT - ReadWriteGridData - ModuleWRFFormat - ERR08'
        !if (.not. LandUse_OK)       stop 'Missing LU_INDEX - ReadWriteGridData - ModuleWRFFormat - ERR09'
        if (.not. CenterX_OK)       stop 'Missing XLAT - ReadWriteGridData - ModuleWRFFormat - ERR10'
        if (.not. CenterY_OK)       stop 'Missing XLONG - ReadWriteGridData - ModuleWRFFormat - ERR11'
        if (.not. ConnectionX_OK)   stop 'Missing XLAT_U - ReadWriteGridData - ModuleWRFFormat - ERR12'
        if (.not. ConnectionY_OK)   stop 'Missing ZLONG_U - ReadWriteGridData - ModuleWRFFormat - ERR13'
            
        if(                  .not. Me%OutputGridNotSpecified) call ConstructGrid 

        if(Me%WriteXYZ .and. .not. Me%OutputGridNotSpecified) call WriteCenterGridXYZ

        call WriteGridToHDF5File

        !Kill Projection
        status = prj90_free(Me%Proj)

    end subroutine ReadWriteGridData

    !--------------------------------------------------------------------------

    subroutine InitializeProjection (ncid)

        !Local-----------------------------------------------------------------
        integer, intent(in)                 :: ncid
        integer                             :: status

        !MAP_PROJ. 1: LAMBERT CONFORMAL, 2: POLAR STEREOGRAPHIC, 3: MERCATOR
                
        status = nf90_get_att(ncid, NF90_GLOBAL, 'MAP_PROJ', Me%ProjType)
        call handle_error(status); if (status /= NF90_NOERR) stop 'InitializeProjection - ModuleWRFFormat - ERR01'


        if (Me%ProjType == 1) then
            
            Me%ProjType = LAMB_CONF_CONIC_
        
            status = nf90_get_att(ncid, NF90_GLOBAL, 'CEN_LAT', Me%CoarseDomainCenterLat)
            call handle_error(status); if (status /= NF90_NOERR) stop 'InitializeProjection - ModuleWRFFormat - ERR02'

            status = nf90_get_att(ncid, NF90_GLOBAL, 'CEN_LON', Me%CoarseDomainCenterLon)
            call handle_error(status); if (status /= NF90_NOERR) stop 'InitializeProjection - ModuleWRFFormat - ERR03'

            status = nf90_get_att(ncid, NF90_GLOBAL, 'TRUELAT1', Me%TrueLatLower) !sp1
            call handle_error(status); if (status /= NF90_NOERR) stop 'InitializeProjection - ModuleWRFFormat - ERR04'

            status = nf90_get_att(ncid, NF90_GLOBAL, 'TRUELAT2', Me%TrueLatUpper) !sp2
            call handle_error(status); if (status /= NF90_NOERR) stop 'InitializeProjection - ModuleWRFFormat - ERR05'

            Me%Params(1) = 'proj=lcc'
            Me%Params(2) = 'ellps=sphere'
            write(Me%Params(3),'(a6,f6.3)') 'lat_1=',Me%TrueLatLower
            write(Me%Params(4),'(a6,f6.3)') 'lat_2=',Me%TrueLatUpper
            write(Me%Params(5),'(a6,f6.3)') 'lon_0=',Me%CoarseDomainCenterLon
            write(Me%Params(6),'(a6,f6.3)') 'lat_0=',Me%CoarseDomainCenterLat
!            Me%Params(7) = 'x_0=1903970.98145531'
!            Me%Params(8) = 'y_0=898179.31322811'


        elseif (Me%ProjType == 2) then

            write(*,*) 'Not ready for Polar Stereographic Map Projection'
            stop 'InitializeProjection - ModuleWRFFormat - ERR99'

        elseif (Me%ProjType == 3) then
            Me%ProjType = SPHERE_MERCATOR_            
            Me%params(1) = '+init=esri:53004'     
        
            !+proj=merc +lat_ts=0 +lon_0=0 +k=1.000000 +x_0=0 +y_0=0
            !+a=6371000 +b=6371000 +units=m +no_defs
                        
        endif

        status=prj90_init(Me%Proj,Me%Params)
        call handle_proj_error(status); if (status /= PRJ90_NOERR) stop 'InitializeProjection - ModuleWRFFormat - ERR07'

    end subroutine InitializeProjection

    !--------------------------------------------------------------------------

    subroutine CorrectCorners
    
        !Local ----------------------------------------------------------------
        !real(kind=kind(1.0d0))                      :: x,y,lat,lon
        real(8)                                     :: x,y,lat,lon
        integer                                     :: ILB,IUB, WILB, WIUB
        integer                                     :: JLB,JUB, WJLB, WJUB
        integer                                     :: i, j, status        

        !Begin ----------------------------------------------------------------

        ! get XLAT_U, XLONG_U -> passar para XY -> 
        !        diminuir o Y de DY/2 e aumentar o ultimo Y de DY/2 
        !        passar a Lat,Lon -> ConnectionX, ConnectionY.


        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB   !south_north
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB   
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB   !west_east
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB

        do j = WJLB, JUB
        do i = WILB, WIUB

            lon = dble(Me%ConnectionX(i,j))
            lat = dble(Me%ConnectionY(i,j))

            status = prj90_fwd(Me%Proj, lon, lat, x, y)
            call handle_proj_error(status); if (status /= PRJ90_NOERR) stop 'CorrectCorners - ModuleWRFFormat - ERR01'

            y = y - dble(Me%DY)/2.

            status = prj90_inv(Me%Proj, x, y, lon, lat)
            call handle_proj_error(status); if (status /= PRJ90_NOERR) stop 'CorrectCorners - ModuleWRFFormat - ERR02'

            Me%ConnectionX(i,j) = lon
            Me%ConnectionY(i,j) = lat

        enddo
        enddo

        do j = WJLB, JUB

            lon = dble(Me%ConnectionX(WIUB, j))
            lat = dble(Me%ConnectionY(WIUB, j))

            status = prj90_fwd(Me%Proj, lon, lat, x, y)
            call handle_proj_error(status); if (status /= PRJ90_NOERR) stop 'CorrectCorners - ModuleWRFFormat - ERR03'

            y = y + dble(Me%DY)/2.

            status = prj90_inv(Me%Proj, x, y, lon, lat)
            call handle_proj_error(status); if (status /= PRJ90_NOERR) stop 'CorrectCorners - ModuleWRFFormat - ERR04'

            Me%ConnectionX(IUB, j) = lon
            Me%ConnectionY(IUB, j) = lat

        enddo

    end subroutine CorrectCorners

    !------------------------------------------------------------------------

    subroutine SetNewDate(current_date)

        !Arguments-------------------------------------------------------------
        character (len=*)                           :: current_date
        type(T_Date), pointer                       :: NewDate
        
        !Local-----------------------------------------------------------------
        real,              dimension(6)             :: TimeReal
        integer                                     :: i
        
        !Begin-----------------------------------------------------------------

        do i=1,len_trim(current_date)

            if (current_date(i:i) =='_'.or.current_date(i:i) ==':'.or. current_date(i:i) =='-') then
                current_date(i:i) = ' '
            endif

        enddo

        read(current_date,*) TimeReal

        call AddDate(Me%FirstDate, NewDate)

        call SetDate(NewDate%Date, Year = TimeReal(1), Month  = TimeReal(2), Day    = TimeReal(3), &
                                   Hour = TimeReal(4), Minute = TimeReal(5), Second = TimeReal(6))



    end subroutine SetNewDate

    !------------------------------------------------------------------------

    subroutine ConstructGrid
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL !, UnitID, i, j
        type(T_Size2D)              :: Size2D
        
        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing grid...'

        Size2D%ILB = Me%WorkSize%ILB
        Size2D%IUB = Me%WorkSize%IUB
        Size2D%JLB = Me%WorkSize%JLB
        Size2D%JUB = Me%WorkSize%JUB


        call WriteGridData  (FileName       = Me%GridFileName,                              &
                             ConnectionX    = Me%ConnectionX,                               &
                             ConnectionY    = Me%ConnectionY,                               &
                             COMENT1        = 'Grid Data created from WRF file',            &
                             COMENT2        = trim(Me%FileName),                            &
                             WorkSize       = Size2D,                                       &
                             CoordType      = SIMPLE_GEOG_,                                 &
                             Xorig          = Me%ConnectionX(1,1),                          &
                             Yorig          = Me%ConnectionY(1,1),                          &
                             Zone           = -99,                                          &
                             GRID_ANGLE     = 0.,                                           &
                             Latitude       = Me%CoarseDomainCenterLat,                     &
                             Longitude      = Me%CoarseDomainCenterLon,                     &
                             FillValue      = -99.,                                         &
                             Overwrite      = .true.,                                       &
                             GridData2D_Real= Me%Bathymetry,                                &
                             ProjType       = Me%ProjType,                                  &
                             Datum          = SPHERE_DATUM,                                 &
                             SP1            = Me%TrueLatLower,                              &
                             SP2            = Me%TrueLatUpper,                              &
                             STAT           = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructGrid - ModuleWRFFormat - ERR01'

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleWRFFormat - ERR02'


    end subroutine ConstructGrid

    !------------------------------------------------------------------------

    subroutine WriteCenterGridXYZ
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, i, j, UnitID
        integer                     :: WILB, WIUB, WJLB, WJUB
        
        !Begin-----------------------------------------------------------------

        WILB = Me%WorkSize%ILB 
        WIUB = Me%WorkSize%IUB 
        WJLB = Me%WorkSize%JLB 
        WJUB = Me%WorkSize%JUB 

        write(*,*)
        write(*,*)'Constructing xyz terrain file...'

        call UnitsManager(UnitID, OPEN_FILE, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'WriteCenterGridXYZ - ModuleWRFFormat - ERR01'

        open(unit=UnitID, status = "unknown", file = trim(Me%GridFileName)//".xyz")
        rewind(UnitID)

        write(UnitID, *)"<begin_xyz>"
        do j = WJLB, WJUB
        do i = WILB, WIUB
            write(UnitID, *)Me%CenterX(i,j), Me%CenterY(i,j), Me%Bathymetry(i,j)
        enddo
        enddo
        write(UnitID, *)"<end_xyz>"


        call UnitsManager(UnitID, CLOSE_FILE, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'WriteCenterGridXYZ - ModuleWRFFormat - ERR02'


    end subroutine WriteCenterGridXYZ

    !----------------------------------------------------------------------

    subroutine WriteGridToHDF5File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer,    dimension(:,:), pointer         :: WaterPoints2D

        !----------------------------------------------------------------------

        allocate(WaterPoints2D(Me%Size%ILB:Me%Size%IUB,&
                               Me%Size%JLB:Me%Size%JUB))
     

        WaterPoints2D = 1
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleWRFFormat - ERR02'

        Me%Bathymetry = -1 * Me%Bathymetry

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleWRFFormat - ERR03'

        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleWRFFormat - ERR02a'

        
        !call HDF5WriteData   (Me%ObjHDF5, "/Grid", "LandUse", "-",       &
        !                      Array2D =  Me%LandUse, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleWRFFormat - ERR03a'

        if(Me%OutputGridNotSpecified)then

            call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB+1,&
                                 Me%WorkSize%JLB, Me%WorkSize%JUB+1, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleWRFFormat - ERR04'


            call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionX", "-",       &
                                  Array2D =  Me%ConnectionX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleWRFFormat - ERR05'

            call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionY", "-",       &
                                  Array2D =  Me%ConnectionY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleWRFFormat - ERR06'

        else

            call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleWRFFormat - ERR07'

        end if
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, &
                             Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleWRFFormat - ERR08'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",    &
                              Array2D = WaterPoints2D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleWRFFormat - ERR09'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleWRFFormat - ERR10'

        deallocate(WaterPoints2D)
        nullify   (WaterPoints2D)

    end subroutine WriteGridToHDF5File

    !------------------------------------------------------------------------

    subroutine ComputeVerticalCoordinate
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer                  :: VerticalZ        
        real, dimension(:,:,:), pointer         :: GeoPP, GeoBS       
        logical                                 :: GeopotentialPerturbation_OK
        logical                                 :: GeopotentialBaseState_OK
        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB, k
        type(T_Field), pointer                  :: Field
        type(T_Date), pointer                   :: CurrentDate


        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Computing Vertical Coordinate...'

        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 
        KLB = Me%Size%KLB; WKLB = Me%WorkSize%KLB 
        KUB = Me%Size%KUB; WKUB = Me%WorkSize%KUB 


        GeopotentialPerturbation_OK  = .false.
        GeopotentialBaseState_OK     = .false.


        CurrentDate => Me%FirstDate

        do while(associated(CurrentDate))

            Field => Me%FirstField

            do while(associated(Field))
            
                if(Field%Name == trim('geopotential base state_3D')         .and. &
                   Field%Date == CurrentDate%Date)then
                    
                    
                    GeoPP => Field%Values3D
                    GeopotentialPerturbation_OK = .true.

                end if

                if(Field%Name == trim('geopotential base state_3D')          .and. &
                   Field%Date == CurrentDate%Date) then
                                        
                    GeoBS => Field%Values3D
                    GeopotentialBaseState_OK    = .true.

                end if

                if(GeopotentialPerturbation_OK .and. GeopotentialBaseState_OK) then

                    call AddField(Me%FirstField, VerticalZ)

                    call SetNewFieldAttributes(Field         =  VerticalZ,          &
                                               Name          = 'VerticalZ',         &
                                               Units         = 'm',                 &
                                               Date          = CurrentDate%Date,    &
                                               nDimensions   = 3,                   &
                                               WorkSize      = Me%WorkSize,         &
                                               Convert       = .true.)

                    allocate(VerticalZ%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))
                    VerticalZ%Values3D = 0.

                    VerticalZ%Values3D(WILB:WIUB, WJLB:WJUB, 0) = - Me%Bathymetry(WILB:WIUB, WJLB:WJUB)

                    VerticalZ%Values3D(WILB:WIUB, WJLB:WJUB, WKLB:WKUB+1) =         &
                                        (GeoPP(WILB:WIUB, WJLB:WJUB, WKLB:WKUB+1) + &
                                         GeoBS(WILB:WIUB, WJLB:WJUB, WKLB:WKUB+1)) / Gravity

!                    write(*,*) 'VerticalZ(20,20,0:28) =', (VerticalZ%Values3D(20,20,k), k=KLB,KUB)

                    GeopotentialPerturbation_OK  = .false.
                    GeopotentialBaseState_OK     = .false.

                    nullify(GeoPP, GeoBS) 
                    exit


               end if

                Field => Field%Next

            end do


            CurrentDate => CurrentDate%Next

        end do 


    end subroutine ComputeVerticalCoordinate


    !------------------------------------------------------------------------

    subroutine ComputeAtmosphericPressure3D
        
        !Local-----------------------------------------------------------------        
        type(T_Field), pointer                  :: Pressure3D
        real, dimension(:,:,:), pointer         :: PB
        logical                                 :: PressureBaseState_OK, Pressure3D_OK
        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        type(T_Field), pointer                  :: Field !, SurfacePressure
        type(T_Date), pointer                   :: CurrentDate


        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Computing Atmospheric Pressure 3D...'


        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 
        KLB = Me%Size%KLB; WKLB = Me%WorkSize%KLB 
        KUB = Me%Size%KUB; WKUB = Me%WorkSize%KUB

        CurrentDate => Me%FirstDate

        do while(associated(CurrentDate))

            ! Compute Pressure 3D

            PressureBaseState_OK = .false.

            Field => Me%FirstField

            do while(associated(Field))

                if(Field%Name == trim('base state pressure_3D').and. Field%Date == CurrentDate%Date) then

                    PB => Field%Values3D
                    PressureBaseState_OK = .true.

                end if

                if(Field%Name == trim(GetPropertyName(AtmosphericPressure_))//"_3D" .and.   &
                   Field%Date == CurrentDate%Date) then

                    Pressure3D => Field
                    Pressure3D_OK = .true.

                endif

                if (PressureBaseState_OK .and. Pressure3D_OK) then

                    Pressure3D%Values3D(WILB:WIUB, WJLB:WJUB, WKLB:WKUB) =                              &
                                                Pressure3D%Values3D(WILB:WIUB, WJLB:WJUB, WKLB:WKUB) +  &
                                                PB(WILB:WIUB, WJLB:WJUB, WKLB:WKUB)

                    Pressure3D%Units    = 'Pa'

                    PressureBaseState_OK = .false.
                    Pressure3D_OK        = .false.

                    nullify(PB)
                    exit

                endif

                Field => Field%Next

            end do

            CurrentDate => CurrentDate%Next

        end do

    end subroutine ComputeAtmosphericPressure3D

    !------------------------------------------------------------------------

    subroutine ComputeAirTemperature3D

        !Local-----------------------------------------------------------------        
        type(T_Field), pointer                  :: Temperature3D
        real, dimension(:,:,:), pointer         :: P, Theta
        logical                                 :: PotTemp3D_OK, Pressure3D_OK
        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        type(T_Field), pointer                  :: Field
        type(T_Date), pointer                   :: CurrentDate
        real                                    :: rcp


        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Computing Air Temperature 3D...'

        rcp = R_dry / Cp_dry

        WILB = Me%WorkSize%ILB 
        WIUB = Me%WorkSize%IUB 
        WJLB = Me%WorkSize%JLB 
        WJUB = Me%WorkSize%JUB 
        WKLB = Me%WorkSize%KLB 
        WKUB = Me%WorkSize%KUB 

        CurrentDate => Me%FirstDate

        PotTemp3D_OK  = .false.
        Pressure3D_OK = .false.


        do while(associated(CurrentDate))

            Field => Me%FirstField

            do while(associated(Field))

                if(Field%Name == trim('potential temperature_3D')                   .and.   &
                   Field%Date == CurrentDate%Date)then

                    Theta => Field%Values3D
                    PotTemp3D_OK = .true.

                end if


                if(Field%Name == trim(GetPropertyName(AtmosphericPressure_))//"_3D" .and.   &
                   Field%Date == CurrentDate%Date)then

                    P => Field%Values3D
                    Pressure3D_OK = .true.  

                end if

                if(PotTemp3D_OK .and. Pressure3D_OK)then

                    call AddField(Me%FirstField, Temperature3D)

                    call SetNewFieldAttributes(Field         =  Temperature3D,                          &
                                               Name          =  trim(GetPropertyName(AirTemperature_))//"_3D", &
                                               Units         = 'K',                                     &
                                               Date          = CurrentDate%Date,                        &
                                               WorkSize      = Me%WorkSize,                             &
                                               nDimensions   = 3,                                       &
                                               Convert       = .true.)

                    allocate(Temperature3D%Values3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
                    Temperature3D%Values3D = null_real

                    Temperature3D%Values3D(WILB:WIUB, WJLB:WJUB, WKLB:WKUB) =               &
                                ( Theta(WILB:WIUB, WJLB:WJUB, WKLB:WKUB) + 300. ) *         &
                                ( P(WILB:WIUB, WJLB:WJUB, WKLB:WKUB) / p1000mb  )** rcp

                    PotTemp3D_OK    = .false.
                    Pressure3D_OK   = .false.

                    nullify(Theta, P)
                    exit

                endif

                Field => Field%Next

            end do

            CurrentDate => CurrentDate%Next

        end do 

        nullify(P,Field)

    end subroutine ComputeAirTemperature3D

    !------------------------------------------------------------------------

    subroutine ComputeMeanSeaLevelPressureMM5

        !Local-----------------------------------------------------------------
        type(T_Field), pointer                  :: MeanSeaLevelPressure
        real, dimension(:,:  ), pointer         :: ReferenceSurfacePressure, Terrain
        real, dimension(:,:,:), pointer         :: Pressure3D, Temperature3D
        logical                                 :: Pref_OK , Pressure_OK, Temperature_OK, Terrain_OK
        type(T_Field), pointer                  :: Field
        type(T_Date), pointer                   :: CurrentDate

        !Begin-----------------------------------------------------------------

        write(*,*) 'Computinng MeanSeaLevelPressure with MM5 algorithm...'

        Temperature_OK  = .false.
        Pref_OK         = .false.
        Pressure_OK     = .false.
        Terrain_OK      = .false.

        CurrentDate => Me%FirstDate

        do while(associated(CurrentDate))

            Field => Me%FirstField

            do while(associated(Field))

                if(Field%Name == GetPropertyName(AtmosphericPressure_)              .and. &
                   Field%Date == CurrentDate%Date)then

                    ReferenceSurfacePressure    => Field%Values2D
                    Pref_OK                     = .true.


                end if

                if(Field%Name == trim(GetPropertyName(AtmosphericPressure_))//"_3D" .and. &
                   Field%Date == CurrentDate%Date)then

                    Pressure3D                  => Field%Values3D
                    Pressure_OK                 = .true.

                end if

                if(Field%Name == trim(GetPropertyName(AirTemperature_))//"_3D"      .and. &
                   Field%Date == CurrentDate%Date)then

                    Temperature3D               => Field%Values3D
                    Temperature_OK              = .true.

                end if

                if(Field%Name == trim('terrain')                                    .and. &
                   Field%Date == CurrentDate%Date)then

                    !nao pode ser Me%Bathymetry pq o terreno muda com 2-way nesting                    
                    Terrain                 => Field%Values2D
                    Terrain_OK              = .true.

                end if

                if(Pref_OK .and. Pressure_OK .and. Temperature_OK .and. Terrain_OK)then

                    call AddField(Me%FirstField, MeanSeaLevelPressure)

                    call SetNewFieldAttributes(Field         =  MeanSeaLevelPressure,                   &
                                               Name          =  GetPropertyName(MeanSeaLevelPressure_), &
                                               Units         = 'Pa',                                    &
                                               Date          = CurrentDate%Date,                        &
                                               WorkSize      = Me%WorkSize,                             &
                                               nDimensions   = 2,                                       &
                                               Convert       = .true.)

                    allocate(MeanSeaLevelPressure%Values2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                    MeanSeaLevelPressure%Values2D = 0.

                    call ComputeSeaPressure_MM5 (MeanSeaLevelPressure, ReferenceSurfacePressure,        &
                                                 Pressure3D, Temperature3D,  Terrain)

                    Temperature_OK  = .false.
                    Pref_OK         = .false.
                    Pressure_OK     = .false.
                    Terrain_OK      = .false.             

                    nullify(ReferenceSurfacePressure, Pressure3D, Temperature3D, Terrain)
                    exit

               end if


                Field => Field%Next

            end do


            CurrentDate => CurrentDate%Next

        end do 

        nullify(Field)

    end subroutine ComputeMeanSeaLevelPressureMM5

    !------------------------------------------------------------------------

    subroutine ComputeMeanSeaLevelPressureWRF

        !Local-----------------------------------------------------------------
        type(T_Field), pointer                  :: MeanSeaLevelPressure
        real, dimension(:,:,:), pointer         :: VerticalZ, Pressure3D, Temperature3D, MixRatio3D
        logical                                 :: VerticalZ_OK, Pressure_OK, Temperature_OK, MixRatio_OK
        type(T_Field), pointer                  :: Field
        type(T_Date), pointer                   :: CurrentDate

        !Begin-----------------------------------------------------------------

        write(*,*) 'Computinng MeanSeaLevelPressure with WRF algorithm...'

        Temperature_OK  = .false.
        Pressure_OK     = .false.
        MixRatio_OK     = .false.

        CurrentDate => Me%FirstDate

        do while(associated(CurrentDate))

            Field => Me%FirstField

            do while(associated(Field))

                if(Field%Name == trim('VerticalZ')                                  .and. &
                   Field%Date == CurrentDate%Date)then

                    VerticalZ                   => Field%Values3D
                    VerticalZ_OK                = .true.

                end if

                if(Field%Name == trim(GetPropertyName(AtmosphericPressure_))//"_3D" .and. &
                   Field%Date == CurrentDate%Date)then

                    Pressure3D                  => Field%Values3D
                    Pressure_OK                 = .true.

                end if

                if(Field%Name == trim(GetPropertyName(AirTemperature_))//"_3D"      .and. &
                   Field%Date == CurrentDate%Date)then

                    Temperature3D               => Field%Values3D
                    Temperature_OK              = .true.

                end if

                if(Field%Name == trim('mixing ratio_3D')                            .and. &
                   Field%Date == CurrentDate%Date)then

                    MixRatio3D               => Field%Values3D
                    MixRatio_OK              = .true.

                end if

                if(VerticalZ_OK .and. Pressure_OK .and. Temperature_OK .and. MixRatio_OK) then

                    call AddField(Me%FirstField, MeanSeaLevelPressure)

                    call SetNewFieldAttributes(Field         =  MeanSeaLevelPressure,                   &
                                               Name          =  'mslp_WRF',                             &
                                               Units         = 'Pa',                                    &
                                               Date          = CurrentDate%Date,                        &
                                               WorkSize      = Me%WorkSize,                             &
                                               nDimensions   = 2,                                       &
                                               Convert       = .true.)

                    allocate(MeanSeaLevelPressure%Values2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                    MeanSeaLevelPressure%Values2D = 0.

                    call compute_seaprs (VerticalZ, Temperature3D, Pressure3D, MixRatio3D, MeanSeaLevelPressure)

                    VerticalZ_OK    = .false.
                    Temperature_OK  = .false.
                    Pressure_OK     = .false.
                    MixRatio_OK     = .false.

                    nullify(VerticalZ, Pressure3D, Temperature3D, MixRatio3D)
                    exit

               end if


                Field => Field%Next

            end do

            CurrentDate => CurrentDate%Next

        end do 

        nullify(Field)

    end subroutine ComputeMeanSeaLevelPressureWRF

    !--------------------------------------------------------------------------

    subroutine ComputeSeaPressure_MM5 (MeanSeaLevelPressure, ReferenceSurfacePressure,        &
                                       Pressure3D, Temperature3D,  Terrain)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                  :: MeanSeaLevelPressure
        real, dimension(:,:  ), pointer         :: ReferenceSurfacePressure, Terrain
        real, dimension(:,:,:), pointer         :: Pressure3D, Temperature3D

        !Local-----------------------------------------------------------------        
        real, dimension(:,:  ), pointer         :: PL, XKLEV, TS, T0               
        logical                                 :: Restart, L1,L2,L3        
        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB, KUPTO
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: i,j,k
        real                                    :: PCONST = 10000.    !Pa (find layer 100 hPa above surface)
        real                                    :: TC = 273.15+17.5 !critical temperature
        integer                                 :: KLO, KHI
        real                                    :: PLO, PHI, TLO, THI, PSFC
        real                                    :: TL, TBAR, HL, XTERM, Tmed



        !Begin-----------------------------------------------------------------
 
        !write(*,*)
        !write(*,*)'Computing mean sea level pressure with MM5 subroutine...'

        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 
        KLB = Me%Size%KLB; WKLB = Me%WorkSize%KLB 
        KUB = Me%Size%KUB; WKUB = Me%WorkSize%KUB 

        if (mod(real(WKUB), 2.) .eq. 0.) then
            KUPTO = WKUB / 2
        else
            KUPTO = WKUB - WKUB/2 + 1
        endif

        XTERM = GAMMA * R_dry / Gravity        
!        Terrain => Me%Bathymetry(WILB:WIUB, WJLB:WJUB)  !Nao pode ser pq o terreno muda com 2-way nesting

        allocate(PL    (WILB:WIUB, WJLB:WJUB))
        allocate(XKLEV (WILB:WIUB, WJLB:WJUB))
        allocate(TS    (WILB:WIUB, WJLB:WJUB))
        allocate(T0    (WILB:WIUB, WJLB:WJUB))

        PL      = null_real
        XKLEV   = null_real
        TS      = null_real
        T0      = null_real

        Restart= .true.

        do while (Restart) 

            Restart = .false.

            !COMPUTE PRESSURE AT PCONST MB ABOVE SURFACE (PL) 

            do j = WJLB, WJUB
            do i = WILB, WIUB
                PL (i,j) = ReferenceSurfacePressure(i,j) - PCONST
                XKLEV(i,j) = 0
            enddo
            enddo

            !FIND 2 LEVELS ON SIGMA SURFACES SURROUNDING PL AT EACH I,J


            do j = WJLB, WJUB
            do i = WILB, WIUB

                do k = WKLB+1, KUPTO
                    if ((Pressure3D(i,j,k)   .lt. PL(i,j)) .and. &
                        (Pressure3D(i,j,k-1) .ge. PL(i,j)))  then
                        XKLEV(i,j) = float(k)
                        exit
                    endif
                enddo

                if (XKLEV(i,j) .lt. 1) then

                    write(*,*) 'Error finding pressure level', PCONST, 'hPa above the surface'
                    write(*,*) 'Last k level =', KUPTO

                    if (KUPTO .ne. 1) then

                        write (*,*) 'Try again with KUPTO=1'
                        KUPTO=1
                        Restart=.true.
                        exit 
                    else

                        write (*,*) 'i,j = ',i,j
                        write (*,*) 'PL = ', PL(i,j)
                        write (*,*) 'ReferenceSurfacePressure = ', ReferenceSurfacePressure(i,j)
                        stop 'ComputeMeanSeaLevelPressure - ModuleMM5Format - ERR01'

                    endif
                endif

            enddo
            enddo
            !if (Restart) exit
        enddo


        !GET TEMPERATURE AT PL (TL), EXTRAPOLATE T AT SURFACE (TS)
        !AND T AT SEA LEVEL (T0) WITH 6.5 K/KM LAPSE RATE
        PL = PL * 0.01

        do j= WJLB, WJUB
        do i= WILB, WIUB

            KLO = NINT(XKLEV(i,j))
            KHI = NINT(XKLEV(i,j))+1

            PLO  = Pressure3D(i,j,KLO) * 0.01
            PHI  = Pressure3D(i,j,KHI) * 0.01
            PSFC = ReferenceSurfacePressure(i,j) * 0.01

            TLO  = Temperature3D(i,j,KLO)
            THI  = Temperature3D(i,j,KHI)

            if ((PLO .le. 1.e-3) .or. (PHI .le. 1.e-3) .or. (PSFC .le. 1.e-3) .or. &
                (PL(i,j) .le. 1.e-3)) then
                write (*,*) 'PLO    =', PLO
                write (*,*) 'PHI    =', PHI
                write (*,*) 'PSFC   =', PSFC
                write (*,*) 'PL(i,j)=', PL(i,j)
                stop
            end if

            TL  = THI-(THI-TLO)*ALOG(PL(i,j)/PHI)/ALOG(PLO/PHI)
            TS(i,j) = TL *(PSFC/PL(i,j))**XTERM
            TBAR = (TS(i,j)+TL)*0.5
            HL   = Terrain(i,j) - R_dry / Gravity * ALOG(PL(i,j)/PSFC)*TBAR
            T0(i,j)=TL+GAMMA*HL    

        enddo
        enddo

        !CORRECT SEA LEVEL TEMPERATURE IF TOO HOT

        do j= WJLB, WJUB
        do i= WILB, WIUB 

            L1=T0(i,j).LT.TC
            L2=TS(i,j).LE.TC
            L3=.NOT.L1

            if (L2.AND.L3) then
                T0(i,j)=TC
            else if ((.NOT. L2).AND. L3) then
                T0(i,j)=TC-0.005*(TS(i,j)-TC)**2
            endif

        enddo
        enddo

        !COMPUTE SEA LEVEL PRESSURE

        do j= WJLB, WJUB
        do i= WILB, WIUB

            Tmed = (TS(i,j)+T0(i,j))/ 2.
            if (Tmed .le. 1.e-3) then
                write (*,*) '(i,j) =', i,j
                write (*,*) 'TS    =', TS(i,j)
                write (*,*) 'T0    =', T0(i,j)
                stop
            endif
            MeanSeaLevelPressure%Values2D (i,j) = ReferenceSurfacePressure(i,j) *       &
                                                  exp(2.* Gravity * Terrain(i,j)/       &
                                                  (R_dry * (TS(i,j)+T0(i,j))))
        enddo
        enddo
        deallocate(PL, XKLEV, TS, T0)

!      print *,'sea pres input at weird location i=20,j=1,k=1'
!      print *,'t=',Temperature3D(20,1,1),Temperature3D(20,2,1),Temperature3D(20,3,1)
!      print *,'terrain=',Terrain(20,1),Terrain(20,2),Terrain(20,3)

!      print *,'psfc=',ReferenceSurfacePressure(20,1), &
!                      ReferenceSurfacePressure(20,2), &
!                      ReferenceSurfacePressure(20,3)

!      print *,'p3d=',Pressure3D(20,1,1),Pressure3D(20,2,1),Pressure3D(20,3,1)

!      print *,'slp=',MeanSeaLevelPressure%Values2D(20,1),     &
!                     MeanSeaLevelPressure%Values2D(20,2),     &
!                     MeanSeaLevelPressure%Values2D(20,3)

    end subroutine ComputeSeaPressure_MM5

!-------------------------------------------------------------------------
!
! This routines has been taken "as is" from wrf_user_fortran_util_0.f
!
! This routine assumes
!    index order is (i,j,k)
!    wrf staggering
!    units: pressure (Pa), temperature(K), height (m), mixing ratio (kg kg{-1})
!    availability of 3d p, t, and qv; 2d terrain; 1d half-level zeta string
!    output units of SLP are Pa, but you should divide that by 100 for the
!          weather weenies.
!    virtual effects are included
!
! Dave

    subroutine compute_seaprs (z, t , p , q , MeanSeaLevelPressure)

! where
!           p = p+pb                                !half levels
!           t = temperature3D
!           q = QVAPOR

! also needs
!           ph = (ph+phb)/9.81                      !full levels
!           z = 0.5*(ph(:,:,1:nz)+ph(:,:,2:nz+1))   !half levels

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                  :: MeanSeaLevelPressure
        real, dimension(:,:,:), pointer         :: z, t, p, q

        !Local-----------------------------------------------------------------        
        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer, dimension(:,:), pointer        :: level
        real   , dimension(:,:), pointer        :: t_surf, t_sea_level

        real                                    :: PCONST = 10000.    !Pa (find layer 100 hPa above surface)
        real                                    :: TC = 273.15+17.5 !critical temperature
        logical                                 :: ridiculous_mm5_test = .TRUE.

        integer                                 :: i,j,k, klo, khi

        REAL plo , phi , tlo, thi , zlo , zhi
        REAL p_at_pconst , t_at_pconst , z_at_pconst
        REAL z_half_lowest

        REAL, PARAMETER :: rcp          = R_dry / Cp_dry

        LOGICAL  l1 , l2 , l3, found

        !Begin-----------------------------------------------------------------
 
        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 
        KLB = Me%Size%KLB; WKLB = Me%WorkSize%KLB 
        KUB = Me%Size%KUB; WKUB = Me%WorkSize%KUB 

        allocate(level      (ILB:IUB,JLB:JUB))
        allocate(t_surf     (ILB:IUB,JLB:JUB))
        allocate(t_sea_level(ILB:IUB,JLB:JUB))

        z = 0.5 * (z(:,:,WKLB:WKUB) + z(:,:,WKLB+1:WKUB+1))

!     Find least zeta level that is PCONST Pa above the surface.  We later use this
!     level to extrapolate a surface pressure and temperature, which is supposed
!     to reduce the effect of the diurnal heating cycle in the pressure field.

        DO j = WJLB, WJUB
        DO i = WILB, WIUB

            level(i,j) = -1

            k = 1
            found = .false.
            do while( (.not. found) .and. (k.le.WKUB))
               IF ( p(i,j,k) .LT. p(i,j,1)-PCONST ) THEN
                  level(i,j) = k
                  found = .true.
               END IF
               k = k+1
            END DO

            IF ( level(i,j) .EQ. -1 ) THEN
            PRINT '(A,I4,A)','Troubles finding level ',   &
                        NINT(PCONST)/100,' above ground.'
            PRINT '(A,I4,A,I4,A)',                        &
                  'Problems first occur at (',i,',',j,')'
            PRINT '(A,F6.1,A)',                           &
                  'Surface pressure = ',p(i,j,1)/100,' hPa.'
            STOP 'Error_in_finding_100_hPa_up'
         END IF


        END DO
        END DO

!     Get temperature PCONST Pa above surface.  Use this to extrapolate
!     the temperature at the surface and down to sea level.

        DO j = WJLB, WJUB
        DO i = WILB, WIUB

            klo = MAX ( level(i,j) - 1 , 1      )
            khi = MIN ( klo + 1        , WKUB - 1 )

            IF ( klo .EQ. khi ) THEN
               PRINT '(A)','Trapping levels are weird.'
               PRINT '(A,I3,A,I3,A)','klo = ',klo,', khi = ',khi, &
                            ': and they should not be equal.'
               STOP 'Error_trapping_levels'
            END IF

            plo = p(i,j,klo)
            phi = p(i,j,khi)
            tlo = t(i,j,klo)*(1. + 0.608 * q(i,j,klo) )
            thi = t(i,j,khi)*(1. + 0.608 * q(i,j,khi) )

            zlo = z(i,j,klo)
            zhi = z(i,j,khi)

            p_at_pconst = p(i,j,1) - pconst
            t_at_pconst = thi-(thi-tlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)
            z_at_pconst = zhi-(zhi-zlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)

            t_surf(i,j) = t_at_pconst*(p(i,j,1)/p_at_pconst)**(GAMMA*R_dry/Gravity)
            t_sea_level(i,j) = t_at_pconst+GAMMA*z_at_pconst

        END DO
        END DO

!     If we follow a traditional computation, there is a correction to the sea level
!     temperature if both the surface and sea level temnperatures are *too* hot.

        IF ( ridiculous_mm5_test ) THEN
            DO j = WJLB, WJUB
            DO i = WILB, WIUB
               l1 = t_sea_level(i,j) .LT. TC
               l2 = t_surf     (i,j) .LE. TC
               l3 = .NOT. l1
               IF ( l2 .AND. l3 ) THEN
                  t_sea_level(i,j) = TC
               ELSE
                  t_sea_level(i,j) = TC - 0.005*(t_surf(i,j)-TC)**2
               END IF
            END DO
            END DO
        END IF

!     The grand finale: ta da!

        DO j = WJLB, WJUB
        DO i = WILB, WIUB
            z_half_lowest=z(i,j,1)


            MeanSeaLevelPressure%Values2D(i,j) = p(i,j,1) *                             &
                                            EXP((2.*Gravity*z_half_lowest)/             &
                                               (R_dry*(t_sea_level(i,j)+t_surf(i,j))))

         !sea_level_pressure(i,j) = sea_level_pressure(i,j)*0.01

        END DO
        END DO

!        print *,'sea pres input at weird location i=20,j=1,k=1'
!        print *,'t=',t(20,1,1),t(20,2,1),t(20,3,1)
!        print *,'z=',z(20,1,1),z(20,2,1),z(20,3,1)
!        print *,'p=',p(20,1,1),p(20,2,1),p(20,3,1)
!        print *,'slp=',MeanSeaLevelPressure%Values2D(20,1),     &
!                        MeanSeaLevelPressure%Values2D(20,2),     &
!                        MeanSeaLevelPressure%Values2D(20,3)

        deallocate(level,t_surf,t_sea_level)

    end subroutine compute_seaprs

    !--------------------------------------------------------------------------
    
    subroutine ComputeRelativeHumidity


        !Local-----------------------------------------------------------------
        type(T_Field), pointer                          :: RelativeHumidity
        type(T_Field), pointer                          :: Field
        type(T_Date), pointer                           :: CurrentDate
        real, dimension(:,:), pointer                   :: SurfacePressure
        real, dimension(:,:), pointer                   :: SurfaceTemperature
        real, dimension(:,:), pointer                   :: SurfaceMixingRatio
        logical                                         :: Pressure_OK      = .false.
        logical                                         :: Temperature_OK   = .false.
        logical                                         :: MixingRatio_OK   = .false.
        integer                                         :: ILB, IUB, JLB, JUB
        integer                                         :: WILB, WIUB, WJLB, WJUB
        integer                                         :: i, j
        real                                            :: es1, qs1

        !Begin-----------------------------------------------------------------
    
        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 

        Pressure_OK      = .false.
        Temperature_OK   = .false.
        MixingRatio_OK   = .false.

        write(*,*)
        write(*,*)'Computing relative humidity...'

        CurrentDate => Me%FirstDate

        do while(associated(CurrentDate))

            Field => Me%FirstField

            do while(associated(Field))

                if(Field%Name == GetPropertyName(AirTemperature_)    .and. &
                   Field%Date == CurrentDate%Date)then

                    SurfaceTemperature  => Field%Values2D
                    Temperature_OK      = .true.

                end if

                if(Field%Name == GetPropertyName(AtmosphericPressure_)    .and. &
                   Field%Date == CurrentDate%Date)then

                    SurfacePressure     => Field%Values2D
                    Pressure_OK         = .true.

                end if

                if(Field%Name == trim("2-meter mixing ratio").and. &
                   Field%Date == CurrentDate%Date)then

                    SurfaceMixingRatio  => Field%Values2D
                    MixingRatio_OK      = .true.

                end if

                if(Temperature_OK .and. Pressure_OK .and. MixingRatio_OK)then

                    call AddField(Me%FirstField, RelativeHumidity)


                    call SetNewFieldAttributes(Field         = RelativeHumidity,                    &
                                               Name          = GetPropertyName(RelativeHumidity_),  &
                                               Units         = '-',                                 &
                                               Date          = CurrentDate%Date,                    &
                                               WorkSize      = Me%WorkSize,                         &
                                               nDimensions   = 2,                                   &
                                               Convert       = .true.)

                    allocate(RelativeHumidity%Values2D(ILB:IUB, JLB:JUB))

                    do j = WJLB, WJUB
                    do i = WILB, WIUB


                       ! from mm5tograds
                       if (SurfaceTemperature(i,j) .le. 273.15 .and. Me%IfSnow .ne.0) then
                            es1 = 6.11 * EXP (22.514-(6150./SurfaceTemperature(i,j)))
                        else
                            es1 = 6.112* EXP (17.67*((SurfaceTemperature(i,j)-273.15)/(SurfaceTemperature(i,j)-29.65)))
                        endif

                        qs1 = 0.622 * es1 / ((SurfacePressure(i,j)/100.) - es1)

                        RelativeHumidity%Values2D(i,j) =  min(max(100. * SurfaceMixingRatio(i,j) / qs1, 5.), 100.)

                        RelativeHumidity%Values2D(i,j) = RelativeHumidity%Values2D(i,j)  * 0.01

                    enddo
                    enddo

                    nullify(SurfaceTemperature)
                    nullify(SurfacePressure   )
                    nullify(SurfaceMixingRatio)

                    Temperature_OK = .false.
                    Pressure_OK    = .false.
                    MixingRatio_OK = .false.


                end if

                Field => Field%Next

            end do


            CurrentDate => CurrentDate%Next

        end do


    end subroutine ComputeRelativeHumidity


    !--------------------------------------------------------------------------
    
    subroutine ComputeRelativeHumidity3D


        !Local-----------------------------------------------------------------
        type(T_Field), pointer                          :: RelativeHumidity
        type(T_Field), pointer                          :: Field
        type(T_Date), pointer                           :: CurrentDate
        real, dimension(:,:,:), pointer                 :: Pressure, Temperature, MixingRatio        
        logical                                         :: Pressure_OK      = .false.
        logical                                         :: Temperature_OK   = .false.
        logical                                         :: MixingRatio_OK   = .false.
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                         :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        integer                                         :: i, j, k
        real                                            :: es1, qs1

        !Begin-----------------------------------------------------------------
    
        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 
        KLB = Me%Size%KLB; WKLB = Me%WorkSize%KLB 
        KUB = Me%Size%KUB; WKUB = Me%WorkSize%KUB 

        Pressure_OK      = .false.
        Temperature_OK   = .false.
        MixingRatio_OK   = .false.

        write(*,*)
        write(*,*)'Computing relative humidity 3D...'

        CurrentDate => Me%FirstDate

        do while(associated(CurrentDate))

            Field => Me%FirstField

            do while(associated(Field))
            
                if(Field%Name == trim(GetPropertyName(AirTemperature_))//"_3D"    .and. &                
                   Field%Date == CurrentDate%Date)then
                    
                    Temperature  => Field%Values3D
                    Temperature_OK      = .true.

                end if

                if(Field%Name == trim(GetPropertyName(AtmosphericPressure_))//"_3D"    .and. &
                   Field%Date == CurrentDate%Date)then
                    
                    Pressure     => Field%Values3D
                    Pressure_OK         = .true.

                end if

                if(Field%Name == trim("mixing ratio")//"_3D" .and. &
                   Field%Date == CurrentDate%Date)then

                    MixingRatio  => Field%Values3D
                    MixingRatio_OK      = .true.

                end if

                if(Temperature_OK .and. Pressure_OK .and. MixingRatio_OK)then

                    
                    call AddField(Me%FirstField, RelativeHumidity)


                    call SetNewFieldAttributes(Field         = RelativeHumidity,                                &
                                               Name          = trim(GetPropertyName(RelativeHumidity_))//"_3D", &
                                               Units         = '-',                                             &
                                               Date          = CurrentDate%Date,                                &
                                               WorkSize      = Me%WorkSize,                                     &
                                               nDimensions   = 3,                                               &
                                               Convert       = .true.)

                    allocate(RelativeHumidity%Values3D(ILB:IUB, JLB:JUB,KLB:KUB))

                    do j = WJLB, WJUB
                    do i = WILB, WIUB
                    do k = WKLB, WKUB

                        ! from mm5tograds                        
                        if (Temperature(i,j,k) .le. 273.15 .and. Me%IfSnow .ne.0) then
                            es1 = 6.11 * EXP (22.514-(6150./Temperature(i,j,k)))
                        else
                            es1 = 6.112* EXP (17.67*((Temperature(i,j,k)-273.15)/(Temperature(i,j,k)-29.65)))
                        endif

                        qs1 = 0.622 * es1 / ((Pressure(i,j,k)/100.) - es1)

                        RelativeHumidity%Values3D(i,j, k) =  min(max(100. * MixingRatio(i,j,k) / qs1, 5.),100.)

                        RelativeHumidity%Values3D(i,j, k) =  RelativeHumidity%Values3D(i,j, k) * 0.01

       
                    enddo
                    enddo
                    enddo

                    nullify(Temperature)
                    nullify(Pressure   )
                    nullify(MixingRatio)

                    Temperature_OK = .false.
                    Pressure_OK    = .false.
                    MixingRatio_OK = .false.


                end if

                Field => Field%Next

            end do


            CurrentDate => CurrentDate%Next

        end do           


    end subroutine ComputeRelativeHumidity3D

     !------------------------------------------------------------------------
    
    subroutine ComputePrecipitation


        !Local-----------------------------------------------------------------
        type(T_Field), pointer                          :: Precipitation
        type(T_Field), pointer                          :: Field
        type(T_Date), pointer                           :: CurrentDate, PreviousDate
        real, dimension(:,:), pointer                   :: AccConvectivePCP
        real, dimension(:,:), pointer                   :: AccNonConvectivePCP
        real, dimension(:,:), pointer                   :: OldAccConvectivePCP
        real, dimension(:,:), pointer                   :: OldAccNonConvectivePCP
        logical                                         :: AccConvectivePCP_OK          = .false.
        logical                                         :: AccNonConvectivePCP_OK       = .false.
        logical                                         :: OldAccConvectivePCP_OK       = .false.
        logical                                         :: OldAccNonConvectivePCP_OK    = .false.
        integer                                         :: ILB, IUB, JLB, JUB
        integer                                         :: WILB, WIUB, WJLB, WJUB
        integer                                         :: i, j
        real                                            :: CurrTotalAccPrecipitation
        real                                            :: OldTotalAccPrecipitation
        real                                            :: TotalPrecipitation        

        !Begin-----------------------------------------------------------------

        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 

        AccConvectivePCP_OK         = .false.
        OldAccConvectivePCP_OK      = .false.     
        AccNonConvectivePCP_OK      = .false.
        OldAccNonConvectivePCP_OK   = .false.
        
        write(*,*)
        write(*,*)'Computing precipitation...'

        CurrentDate     => Me%FirstDate
        PreviousDate    => CurrentDate
                

do1:    do while(associated(CurrentDate))

            Field => Me%FirstField

do2:        do while(associated(Field))
            
                if(Field%Name == trim(AccConvPrecipitation) .and. &
                   Field%Date == CurrentDate%Date)then
                    
                    AccConvectivePCP        => Field%Values2D
                    AccConvectivePCP_OK     = .true.

                end if

                if(Field%Name == trim(AccConvPrecipitation) .and. &
                   Field%Date == PreviousDate%Date)then
                    
                    OldAccConvectivePCP        => Field%Values2D
                    OldAccConvectivePCP_OK     = .true.

                end if

                if(Field%Name == trim(AccNonConvPrecipitation) .and. &
                   Field%Date == CurrentDate%Date)then
                    
                    AccNonConvectivePCP     => Field%Values2D
                    AccNonConvectivePCP_OK  = .true.

                end if

                if(Field%Name == trim(AccNonConvPrecipitation) .and. &
                   Field%Date == PreviousDate%Date)then
                    
                    OldAccNonConvectivePCP     => Field%Values2D
                    OldAccNonConvectivePCP_OK  = .true.

                end if

                if(AccConvectivePCP_OK    .and. AccNonConvectivePCP_OK .and. &
                   OldAccConvectivePCP_OK .and. OldAccNonConvectivePCP_OK )then
                    
                    call AddField(Me%FirstField, Precipitation)

                    call SetNewFieldAttributes(Field         = Precipitation,                       &
                                               Name          = GetPropertyName(Precipitation_),     &
                                               Units         = 'mm',                                &
                                               Date          = CurrentDate%Date,                    &
                                               WorkSize      = Me%WorkSize,                         &
                                               nDimensions   = 2,                                   &
                                               Convert       = .true.)

                    allocate(Precipitation%Values2D(ILB:IUB, JLB:JUB))

                    if (CurrentDate%Date == Me%FirstDate%Date) then

                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                                                    !mm
                            TotalPrecipitation   =  10.0 * ( AccConvectivePCP(i,j) + AccNonConvectivePCP (i,j) ) 

                            if (TotalPrecipitation .LE. AllmostZero) TotalPrecipitation = 0.0

                            Precipitation%Values2D(i,j) = TotalPrecipitation


                        enddo
                        enddo
                            
                    else

                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                            CurrTotalAccPrecipitation   =  AccConvectivePCP(i,j)    + AccNonConvectivePCP (i,j)
                            OldTotalAccPrecipitation    =  OldAccConvectivePCP(i,j) + OldAccNonConvectivePCP (i,j)

                                                !mm
                            TotalPrecipitation = 10.0 * ( CurrTotalAccPrecipitation - OldTotalAccPrecipitation )

                            if (TotalPrecipitation .LE. AllmostZero) TotalPrecipitation = 0.0

                            Precipitation%Values2D(i,j) = TotalPrecipitation


                        enddo
                        enddo

                    endif

                    nullify(AccConvectivePCP        )
                    nullify(OldAccConvectivePCP     )
                    nullify(AccNonConvectivePCP     )
                    nullify(OldAccNonConvectivePCP  )

                    AccConvectivePCP_OK         = .false.
                    OldAccConvectivePCP_OK      = .false.     
                    AccNonConvectivePCP_OK      = .false.
                    OldAccNonConvectivePCP_OK   = .false.

                end if

                Field => Field%Next

            end do do2

            PreviousDate    => CurrentDate
            CurrentDate     => CurrentDate%Next

        end do  do1   


    end subroutine ComputePrecipitation

    !------------------------------------------------------------------------

    subroutine ComputeWindShearStress
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer                          :: Field
        type(T_Date), pointer                           :: CurrentDate
        real, dimension(:,:), pointer                   :: WindVelocityX
        real, dimension(:,:), pointer                   :: WindVelocityY
        real, dimension(:,:), pointer                   :: WindShearVelocity
        logical                                         :: VelocityX_OK     = .false.
        logical                                         :: VelocityY_OK     = .false.
        logical                                         :: ShearVelocity_OK = .false.
        type(T_Field), pointer                          :: WindStressX
        type(T_Field), pointer                          :: WindStressY
        integer                                         :: ILB, IUB, JLB, JUB
        integer                                         :: WILB, WIUB, WJLB, WJUB
        integer                                         :: i, j
        real                                            :: WindModulus

        !Begin-----------------------------------------------------------------
        
        WILB = Me%WorkSize%ILB 
        WIUB = Me%WorkSize%IUB 
        WJLB = Me%WorkSize%JLB 
        WJUB = Me%WorkSize%JUB 

        ILB  = Me%Size%ILB
        IUB  = Me%Size%IUB
        JLB  = Me%Size%JLB
        JUB  = Me%Size%JUB

        write(*,*)
        write(*,*)'Computing wind shear stress...'

        CurrentDate => Me%FirstDate

        do while(associated(CurrentDate))

            Field => Me%FirstField

            do while(associated(Field))
            
                if(Field%Name == GetPropertyName(WindVelocityX_)    .and. &
                   Field%Date == CurrentDate%Date)then
                    
                    WindVelocityX       => Field%Values2D
                    VelocityX_OK        = .true.

                end if

                if(Field%Name == GetPropertyName(WindVelocityY_)    .and. &
                   Field%Date == CurrentDate%Date)then
                    
                    WindVelocityY       => Field%Values2D
                    VelocityY_OK        = .true.

                end if

                if(Field%Name == GetPropertyName(WindShearVelocity_).and. &
                   Field%Date == CurrentDate%Date)then

                    WindShearVelocity   => Field%Values2D
                    ShearVelocity_OK    = .true.

                end if

                if(VelocityX_OK .and. VelocityY_OK .and. ShearVelocity_OK)then

                    
                    call AddField(Me%FirstField, WindStressX)
                    call AddField(Me%FirstField, WindStressY)

                    call SetNewFieldAttributes(Field         = WindStressX,                      &
                                               Name          = GetPropertyName(WindStressX_),    &
                                               Units         = 'N/m2',                           &
                                               Date          = CurrentDate%Date,                 &
                                               WorkSize      = Me%WorkSize,                      &
                                               nDimensions   = 2,                                &
                                               Convert       = .true.)

                    call SetNewFieldAttributes(Field         = WindStressY,                      &
                                               Name          = GetPropertyName(WindStressY_),    &
                                               Units         = 'N/m2',                           &
                                               Date          = CurrentDate%Date,                 &
                                               WorkSize      = Me%WorkSize,                      &
                                               nDimensions   = 2,                                &
                                               Convert       = .true.)


                    allocate(WindStressX%Values2D(ILB:IUB, JLB:JUB))
                    allocate(WindStressY%Values2D(ILB:IUB, JLB:JUB))


                    do i = WILB, WIUB
                    do j = WJLB, WJUB

                        WindModulus = sqrt(WindVelocityX(i,j)**2+WindVelocityY(i,j)**2.)

                        if (WindModulus > 0.) then

                            !N/m2      = [kg/m^3] * [m2/s2]       * []
                            WindStressX%Values2D(i,j) = Air_Density * WindShearVelocity(i,j)**2. * WindVelocityX(i,j)/WindModulus


                            !N/m2      = [kg/m^3] * [m2/s2]       * []
                            WindStressY%Values2D(i,j) = Air_Density * WindShearVelocity(i,j)**2. * WindVelocityY(i,j)/WindModulus
                        else

                            WindStressX%Values2D(i,j) = 0.
                            WindStressY%Values2D(i,j) = 0.

                        endif

                    enddo
                    enddo

                    nullify(WindVelocityX    )
                    nullify(WindVelocityY    )
                    nullify(WindShearVelocity)

                    VelocityX_OK     = .false.
                    VelocityY_OK     = .false.
                    ShearVelocity_OK = .false.


                end if

                Field => Field%Next

            end do


            CurrentDate => CurrentDate%Next

        end do           
            
    end subroutine ComputeWindShearStress
    

    !------------------------------------------------------------------------

    subroutine ComputeWindModulus
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer                          :: Field
        type(T_Date), pointer                           :: CurrentDate
        real, dimension(:,:), pointer                   :: WindVelocityX
        real, dimension(:,:), pointer                   :: WindVelocityY
        logical                                         :: VelocityX_OK     = .false.
        logical                                         :: VelocityY_OK     = .false.
        type(T_Field), pointer                          :: WindModulus, WindDirection
        integer                                         :: ILB, IUB, JLB, JUB
        integer                                         :: WILB, WIUB, WJLB, WJUB
        integer                                         :: i, j
        real                                            :: wd

        !Begin-----------------------------------------------------------------
        
        WILB = Me%WorkSize%ILB 
        WIUB = Me%WorkSize%IUB 
        WJLB = Me%WorkSize%JLB 
        WJUB = Me%WorkSize%JUB 

        ILB  = Me%Size%ILB
        IUB  = Me%Size%IUB
        JLB  = Me%Size%JLB
        JUB  = Me%Size%JUB

        write(*,*)
        write(*,*)'Computing wind modulus ...'

        CurrentDate => Me%FirstDate

        do while(associated(CurrentDate))

            Field => Me%FirstField

            do while(associated(Field))
            
                if(Field%Name == GetPropertyName(WindVelocityX_)    .and. &
                   Field%Date == CurrentDate%Date)then
                    
                    WindVelocityX       => Field%Values2D
                    VelocityX_OK        = .true.

                end if

                if(Field%Name == GetPropertyName(WindVelocityY_)    .and. &
                   Field%Date == CurrentDate%Date)then
                    
                    WindVelocityY       => Field%Values2D
                    VelocityY_OK        = .true.

                end if


                if(VelocityX_OK .and. VelocityY_OK) then

                    
                    call AddField(Me%FirstField, WindModulus)
                    call AddField(Me%FirstField, WindDirection)                    

                    call SetNewFieldAttributes(Field         = WindModulus,                         &
                                               Name          = GetPropertyName(WindModulus_),       &
                                               Units         = 'm/s',                               &
                                               Date          = CurrentDate%Date,                    &
                                               WorkSize      = Me%WorkSize,                         &
                                               nDimensions   = 2,                                   &
                                               Convert       = .true.)

                    call SetNewFieldAttributes(Field         = WindDirection,                       &
                                               Name          = GetPropertyName(WindDirection_),     &
                                               Units         = 'degrees',                           &
                                               Date          = CurrentDate%Date,                    &
                                               WorkSize      = Me%WorkSize,                         &
                                               nDimensions   = 2,                                   &
                                               Convert       = .true.)


                    allocate(WindModulus%Values2D   (ILB:IUB, JLB:JUB))
                    allocate(WindDirection%Values2D (ILB:IUB, JLB:JUB))


                    do i = WILB, WIUB
                    do j = WJLB, WJUB

                        WindModulus%Values2D(i,j)= sqrt(WindVelocityX(i,j)**2+WindVelocityY(i,j)**2.)

                        wd = 3.* Pi / 2. - atan2(WindVelocityY(i,j), WindVelocityX(i,j))
                        wd = wd * 180 / Pi

                        if (wd .gt. 360.) wd = wd - 360.

                        WindDirection%Values2D(i,j) = wd
                        
                    enddo
                    enddo

                    nullify(WindVelocityX    )
                    nullify(WindVelocityY    )
                    
                    VelocityX_OK     = .false.
                    VelocityY_OK     = .false.
                    

                end if

                Field => Field%Next

            end do


            CurrentDate => CurrentDate%Next

        end do           
            
    end subroutine ComputeWindModulus

    !------------------------------------------------------------------------
          
    subroutine SetNewFieldAttributes(Field, Name, Units, Date, WorkSize, nDimensions, Convert)

        !Arguments-------------------------------------------------------------
        type(T_Field),    pointer                 :: Field
        character(len=*),           intent(IN)    :: Name
        character(len=*),           intent(IN)    :: Units
        type(T_Time),               intent(IN)    :: Date
        logical,                    intent(IN)    :: Convert
        integer,       optional,    intent(IN)    :: nDimensions
        type(T_Size3D),optional,    intent(IN)    :: WorkSize

        !Begin-----------------------------------------------------------------

        Field%Name        = Name ;  Field%Date        = Date

        Field%Units       = Units;  Field%Convert     = Convert

        if(present(WorkSize   ))    Field%WorkSize    = WorkSize

        if(present(nDimensions))    Field%nDimensions = nDimensions


    end subroutine SetNewFieldAttributes

    !------------------------------------------------------------------------
    
    subroutine AddField (FirstField, ObjField)

        !Arguments-------------------------------------------------------------
        type (T_Field), pointer                   :: FirstField
        type (T_Field), pointer                   :: ObjField

        !Local-----------------------------------------------------------------
        type (T_Field), pointer                   :: NewField
        type (T_Field), pointer                   :: PreviousField
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewField)
        nullify  (NewField%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstField)) then
            FirstField         => NewField
            ObjField           => NewField
        else
            PreviousField      => FirstField
            ObjField           => FirstField%Next
            do while (associated(ObjField))
                PreviousField  => ObjField
                ObjField       => ObjField%Next
            enddo
            ObjField           => NewField
            PreviousField%Next => NewField
        endif


    end subroutine AddField
    
    
    !------------------------------------------------------------------------


    subroutine AddDate (FirstDate, ObjDate)

        !Arguments-------------------------------------------------------------
        type (T_Date), pointer                   :: FirstDate
        type (T_Date), pointer                   :: ObjDate

        !Local-----------------------------------------------------------------
        type (T_Date), pointer                   :: NewDate
        type (T_Date), pointer                   :: PreviousDate
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewDate)
        nullify  (NewDate%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstDate)) then
            FirstDate         => NewDate
            ObjDate           => NewDate
        else
            PreviousDate      => FirstDate
            ObjDate           => FirstDate%Next
            do while (associated(ObjDate))
                PreviousDate  => ObjDate
                ObjDate       => ObjDate%Next
            enddo
            ObjDate           => NewDate
            PreviousDate%Next => NewDate
        endif


    end subroutine AddDate
    
   !----------------------------------------------------------------------    

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleWRFFormat - ERR01'
        

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

    logical function VariableIsToRead(WRFName, MohidName)
        
        !Arguments-----------------------------------------------------------
        character(Len=*)                :: WRFName
        character(Len=StringLength)     :: MohidName
        integer                         :: i
        
        !Begin-----------------------------------------------------------------

        VariableIsToRead = .false.

        select case(trim(WRFName))


            !2D ---------------------------------------------------------------

            case('U10')

                MohidName = GetPropertyName(WindVelocityX_)

            case('V10')

                MohidName = GetPropertyName(WindVelocityY_)

            case('UST')

                MohidName = GetPropertyName(WindShearVelocity_)

            case('PSFC')

                MohidName = GetPropertyName(AtmosphericPressure_)


            case('T2') 

                MohidName = GetPropertyName(AirTemperature_)

            case('Q2')

                MohidName = "2-meter mixing ratio"


            case('GLW')

                MohidName = GetPropertyName(DownwardLongWaveRadiation_)

            case('SWDOWN')

                MohidName = GetPropertyName(SolarRadiation_)

            case('LH')

                MohidName = GetPropertyName(LatentHeat_)

            case('HFX')

                MohidName = GetPropertyName(SensibleHeat_)

            case('OLR')

                MohidName = GetPropertyName(UpwardLongWaveRadiation_)


            case('SST')

                MohidName = "sea water temperature"

            case('PBLH')

                MohidName = "PBL height"

            case('ALBEDO')

                MohidName = GetPropertyName(Albedo_)

            case('RAINC')      

                MohidName = AccConvPrecipitation


            case('RAINNC')      

                MohidName = AccNonConvPrecipitation


            case('TSK') !SURFACE SKIN TEMPERATURE ou TSLB (SOIL TEMPERATURE)?

                MohidName = "ground temperature"

            !3D ---------------------------------------------------------------

            case('U')

                MohidName = GetPropertyName(WindVelocityX_)
                MohidName = trim(MohidName)//"_3D"

            case('V')

                MohidName = GetPropertyName(WindVelocityY_)
                MohidName = trim(MohidName)//"_3D"

            case('W')

                MohidName = 'wind velocity Z'
                MohidName = trim(MohidName)//"_3D"

            case('T')

                !Potential Temperature
                MohidName = "potential temperature_3D"
!                 MohidName = GetPropertyName(AirTemperature_)
!                 MohidName = trim(MohidName)//"_3D"

            case('P')

                !Pressure Perturbation
                MohidName = GetPropertyName(AtmosphericPressure_)
                MohidName = trim(MohidName)//"_3D"

            case('PB')

                MohidName = 'base state pressure_3D'

            case('PH')

                MohidName = 'geopotential perturbation_3D'

            case('PHB')

                MohidName = 'geopotential base state_3D'

            case('QVAPOR')

                MohidName = 'mixing ratio_3D'

            case('QCLOUD')

                MohidName = 'cloud water mixing ratio_3D'

            case('QRAIN')

                MohidName = 'rain water mixing ratio_3D'

            case('HGT')

                MohidName = 'terrain'

            case default

                MohidName = ''

        end select


        do i = 1, size(Me%FieldsToRead)

            if(trim(MohidName) == trim(Me%FieldsToRead(i)))then
                VariableIsToRead = .true.
            end if

        enddo


    end function VariableIsToRead

    !------------------------------------------------------------------------    

    logical function FieldIsToRead(FieldName)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: FieldName

        !Local----------------------------------------------------------------
        integer                                     :: i

        !Begin----------------------------------------------------------------
        
        FieldIsToRead = .false.

        do i = 1, size(Me%FieldsToRead)

            if(Me%FieldsToRead(i) == trim(FieldName))then
                FieldIsToRead = .true.
            end if

        enddo

    end function FieldIsToRead

    !--------------------------------------------------------------------------

    subroutine ConvertToMohidUnits(Field)

        !Local-----------------------------------------------------------------
        type (T_Field), pointer                     :: Field
        integer                                     :: i,j,k
        
        !Begin-----------------------------------------------------------------

        ! tem que ser assim pq as unidades do netcdf sao todas diferentes - formatacao 

        select case(trim(Field%Name))


            case('air temperature')

                Field%Units     = '�C'
                Field%Values2D  = Field%Values2D - AbsoluteZero

            case('air temperature_3D')

                Field%Units     = '�C'
                Field%Values3D  = Field%Values3D - AbsoluteZero
            
            case('wind velocity X')

                Field%Units     = 'm/s'

            case('wind velocity Y')

                Field%Units     = 'm/s'

            case('wind velocity X_3D')

                Field%Units     = 'm/s'

            case('wind velocity Y_3D')

                Field%Units     = 'm/s'

            case('wind velocity Z_3D')

                Field%Units     = 'm/s'

            case('cloud water mixing ratio_3D')

                Field%Units     = 'kg/kg'       
                
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                    if (Field%Values3D(i,j,k) .le. AlmostZero) Field%Values3D(i,j,k) = 0.0
                enddo
                enddo
                enddo

            case('mixing ratio_3D')

                Field%Units     = 'kg/kg'

                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                    if (Field%Values3D(i,j,k) .le. AlmostZero) Field%Values3D(i,j,k) = 0.0
                enddo
                enddo
                enddo

            case default


        end select



    end subroutine ConvertToMohidUnits
    
    !------------------------------------------------------------------------
    
    subroutine WriteMapping
        
        !Local-----------------------------------------------------------------
        integer,    dimension(:,:,:), pointer           :: WaterPoints3D
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------
        

        allocate(WaterPoints3D(Me%Size%ILB:Me%Size%IUB,&
                               Me%Size%JLB:Me%Size%JUB,&
                               Me%Size%KLB:Me%Size%KUB))

        WaterPoints3D = 1


        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, &
                             Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteMapping - ModuleWRFFormat - ERR10'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3D", "-",    &
                              Array3D = WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteMapping - ModuleWRFFormat - ERR20'

        deallocate(WaterPoints3D)
        nullify   (WaterPoints3D)

    end subroutine WriteMapping
    
    !--------------------------------------------------------------------------
   
    subroutine OutputFields

        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        integer                                         :: STAT_CALL, OutputNumber
        type(T_Field), pointer                          :: Field
        type(T_Date), pointer                           :: CurrentDate, PreviousOutputDate
        real                                            :: dt

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Writing HDF5 file...'

        call WriteMapping

        OutputNumber = 1
        CurrentDate        => Me%FirstDate
        PreviousOutputDate => CurrentDate
        dt = 0.0

        do while(associated(CurrentDate))

            call ExtractDate   (CurrentDate%Date,                           &
                                AuxTime(1), AuxTime(2), AuxTime(3),         &
                                AuxTime(4), AuxTime(5), AuxTime(6))


ifT:        if (CurrentDate%Date .GE. Me%StartTime .AND. CurrentDate%Date .LE. Me%EndTime) then            
            dt = CurrentDate%Date - PreviousOutputDate%Date
ifDT:       if (CurrentDate%Date .EQ. Me%FirstDate%Date .OR. dt >= Me%OutputDTInterval) then

            TimePtr => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleWRFFormat - ERR30'


            call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                 "Time", "YYYY/MM/DD HH:MM:SS",             &
                                 Array1D = TimePtr,                         &
                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleWRFFormat - ERR40'


            Field => Me%FirstField

            do while(associated(Field))

                if(Field%Convert)then
            
                    if(Field%Date == CurrentDate%Date)then
                
                        Field%OutputNumber = OutputNumber

                        call ConvertToMohidUnits(Field)

                        if(Field%nDimensions == 2)then
                           
                            call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                               Field%WorkSize%JLB, Field%WorkSize%JUB, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleWRFFormat - ERR50'

                            call HDF5WriteData(Me%ObjHDF5,                                      &
                                               "/Results/"//Field%Name,                         &
                                               Field%Name,                                      &
                                               Field%Units,                                     &
                                               Array2D      = Field%Values2D,                   &
                                               OutputNumber = Field%OutputNumber,               &
                                               STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleWRFFormat - ERR60'


                        elseif(Field%nDimensions == 3)then

                            if(Field%Name == 'VerticalZ')then

                                call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                                   Field%WorkSize%JLB, Field%WorkSize%JUB,            &
                                                   Field%WorkSize%KLB - 1, Field%WorkSize%KUB, STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleWRFFormat - ERR70'

                                !Invert vertical axis
                                Field%Values3D = -1. * Field%Values3D

                                call HDF5WriteData(Me%ObjHDF5,                                  &
                                               "/Grid/"//Field%Name,                            &
                                               'Vertical',                                      &
                                               Field%Units,                                     &
                                               Array3D      = Field%Values3D,                   &
                                               OutputNumber = Field%OutputNumber,               &
                                               STAT         = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleWRFFormat - ERR80'

                            else
                    
                                call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                                   Field%WorkSize%JLB, Field%WorkSize%JUB,            &
                                                   Field%WorkSize%KLB, Field%WorkSize%KUB, STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleWRFFormat - ERR90'

                                call HDF5WriteData(Me%ObjHDF5,                                      &
                                                   "/Results3D/"//Field%Name,                       &
                                                   Field%Name,                                      &
                                                   Field%Units,                                     &
                                                   Array3D      = Field%Values3D,                   &
                                                   OutputNumber = Field%OutputNumber,               &
                                                   STAT         = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleWRFFormat - ERR100'

                            end if
                        end if
                    end if
                end if

                Field => Field%Next

            end do

            OutputNumber = OutputNumber + 1
            PreviousOutputDate => CurrentDate

            endif ifDT            
            endif ifT

            CurrentDate  => CurrentDate%Next

        end do

        write(*,*)
        write(*,*)'Closing HDF5 file...'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleWRFFormat - ERR110'



    end subroutine OutputFields

    !--------------------------------------------------------------------------

    subroutine KillWRFFormat
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        if(.not. Me%OutputGridNotSpecified)then
            call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillWRFFormat - ModuleWRFFormat - ERR01'

            deallocate(Me%Bathymetry )
            deallocate(Me%CenterX    )
            deallocate(Me%ConnectionX)
            deallocate(Me%CenterY    )
            deallocate(Me%ConnectionY)

        end if

        
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillWRFFormat - ModuleWRFFormat - ERR02'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillWRFFormat - ModuleWRFFormat - ERR03'



        deallocate(Me%FirstField)
        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillWRFFormat

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine handle_error (status)

        integer, intent(in)         :: status

        if (status /= NF90_NOERR) write(*,*) trim(nf90_strerror(status))        

    end subroutine handle_error

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine handle_proj_error (status)

        integer, intent(in)         :: status

        if (status /= PRJ90_NOERR) write(*,*) trim(prj90_strerrno(status))

    end subroutine handle_proj_error    

    !--------------------------------------------------------------------------
 
end module ModuleWRFFormat

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior T�cnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------










