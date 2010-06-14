!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : MM5Format
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2004
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to convert MM5Format files into HDF5 format.
!                 Compiler settings: Unformatted File Conversion = BIG ENDIAN
!                 Optimized for CVF6.6. Written in ANSI FORTRAN 95
!
!------------------------------------------------------------------------------
!DataFile
!
!   FILENAME                    : char              -           !Path to MM5 original file
!   TERRAIN_FILENAME            : char              -           !Path to MM5 TERRAIN file
!   OUTPUT_GRID_FILENAME        : char              -           !Path to grid data file generated from MM5 file
!   OUTPUTFILENAME              : char              -           !Path to HDF5 file generated from MM5 file
!   WRITE_XYZ                   : 0/1              [0]          !Flag to write xyz center grid cells
!   WRITE_TERRAIN               : 0/1              [0]          !Flag to write MM5 TERRAIN fields
!   WRITE_HEADERS               : 0/1              [0]          !Flag to output MM5 ans TERRAIN headers in final HDF5 file
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

!---Available properties in MM5 output file (copy from MOHID NAME column to input data file)

!         ---MM5 NAME---    ---MOHID NAME---
!
!            T2             air temperature
!            PSTARCRS       atmospheric pressure
!            U10            wind velocity X
!            V10            wind velocity Y
!            UST            wind shear velocity
!            LHFLUX         latent heat
!            SWDOWN         sensible heat
!            SWDOWN         solar radiation
!            LWDOWN         infrared radiation
!            SWOUT          top outgoing shortwave radiation
!            LWOUT          top outgoing longwave radiation
!            SOIL T 1       soil temperature layer 1
!            SOIL T 1       soil temperature layer 2
!            SOIL T 1       soil temperature layer 3
!            SOIL T 1       soil temperature layer 4
!            SOIL T 1       soil temperature layer 5
!            SOIL T 1       soil temperature layer 6
!            Q2             2-meter mixing ratio
!            TSEASFC        sea water temperature
!            PBL HGT        PBL height
!            PBL REGIME     PBL regime
!            RAIN CON       accumulated convective precipitation        (cm)
!            RAIN NON       accumulated non-convective precipitation    (cm)
!            GROUND T       ground temperature
!            RES TEMP       infinite reservoir slab temperature
!            U              wind velocity X_3D
!            V              wind velocity Y_3D
!            W              wind velocity Z_3D
!            T              air temperature_3D
!            PP             atmospheric pressure_3D
!            Q              mixing ratio_3D
!            CLW            cloud water mixing ratio_3D
!            RNW            rain water mixing ratio_3D
!            ICE            cloud ice mixing ratio_3D
!            SNOW           snow mixing ratio_3D
!            RAD TEND       atmospheric radiation tendency_3D


!Some notes on how MM5 file is organized and read
!   nDimensions    = dimension of the field
!   start_index    = starting indices of the field array
!   end_index      = ending indices of the field array
!   time           = the integration or forecast time for this field
!   staggering     = whether the field is at dot('D') or cross point('C')
!   ordering       = the order of the field array dimension
!                      - YXP: 3-D field, pressure data dimensioned by (IX,JX,KXP)
!                      - YXS: 3-D field, sigma data dimensioned by (IX,JX,KXS)
!                      - YXW: 3-D field, sigma data dimensioned by (IX,JX,KXS+1)
!                      - YX : 2-D field, with array dimensioned by (IX,JX) with IX in Y direction)
!                      - CA : 2-D field, with array dimensioned by (landuse-categories,2).
!                             Arrays to store land property values, such as albedo, roughness length, etc.
!                      - XSB: 3-D field, containing north and south boundary arrays, dimensioned by (JX,KXS,5) 
!                      - YSB: 3-D field, containing west and east boundary arrays, dimensioned by (IX,KXS,5) 
!                      - XWB: 3-D field, containing north and south boundary arrays for vertical motion, 
!                             dimensioned by (JX,KXS+1,5)
!                      - YWB: 3-D field, containing west and east boundary arrays for vertical motion, 
!                             dimensioned by (IX,KXS+1,5) 
!                      - P  : 1-D field, pressure level array
!                      - S  : 1-D field, sigma level array
!   current_date   = representation of date valid for this field
!   MM5Name        = MM5 field name
!   Units          = field units
!   description    = field description


Module ModuleMM5Format

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalGrid

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertMM5Format
    private ::      ReadOptions
    private ::          ReadFieldsToConvert
    private ::          FieldIsToConvert
    private ::      Open_HDF5_OutPut_File
    private ::      OpenAndReadMM5TerrainFile
    private ::          ConstructGrid
    private ::          WriteCenterGridXYZ
    private ::          WriteGridToHDF5File
    private ::      OpenAndReadMM5File
    private ::          AddField
    private ::          AddDate
    private ::          SetNewDate
    private ::          AttributeDimensions
    private ::          FieldIsToRead
    private ::          CheckGridInformation
    private ::          SetNewFieldAttributes
    private ::          ConvertToMohidUnits2D
    private ::          ConvertToMohidUnits3D
    private ::          FillGridInformation
    private ::      BeginEndTime
    private ::      ComputeWindShearStress
    private ::      ComputeVerticalCoordinate
    private ::      ComputeMeanSeaLevelPressure
    private ::      ComputeRelativeHumidity
    private ::      ComputeRelativeHumidity3D
    private ::      ComputePrecipitation
    private ::      OutputFields
    private ::      KillMM5Format
    

    !MM5-----------------------------------------------------------------------
    private ::          printout_big_header

    !Parameters----------------------------------------------------------------
    integer, parameter                                      :: DotGrid              = 1          !Cell corners
    integer, parameter                                      :: CrossGrid            = 2          !Cell center
    integer, parameter                                      :: FullLevels           = 1          !layer face
    integer, parameter                                      :: HalfLevels           = 2          !layer center

    !Perfect gas constant in Pa m3 kg-1 K-1
    integer, parameter                                      :: Perfect_gas_R        = 287.04 !Pa m3 kg-1 K-1

    character(StringLength), parameter  :: AccConvPrecipitation     = "accumulated convective precipitation"
    character(StringLength), parameter  :: AccNonConvPrecipitation  = "accumulated non-convective precipitation"


    !Types---------------------------------------------------------------------
    type       T_Date
        type(T_Time)                                        :: Date
        type(T_Date), pointer                               :: Next
    end type  T_Date

    type       T_Field
        character(len=StringLength)                         :: Name
        character(len=StringLength)                         :: Units
        logical                                             :: Convert                  = .false.
        integer                                             :: nDimensions
        integer                                             :: GridLocation
        integer                                             :: GridLocationZ
        type(T_Time)                                        :: Date
        real, dimension(:,:  ),     pointer                 :: Values2D
        real, dimension(:,:,:),     pointer                 :: Values3D
        integer                                             :: OutputNumber             = 1
        type(T_Size3D)                                      :: Size, WorkSize
        type(T_Field),              pointer                 :: Next
    end type  T_Field

 
    type       T_MM5Format
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
        logical                                             :: ComputeMeanSeaLevelPressure = .false.
        logical                                             :: WriteTerrain             = .false.
        logical                                             :: WriteHeaders             = .false.
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

    end type  T_MM5Format

    type(T_MM5Format), pointer                              :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertMM5Format(EnterDataID, ClientNumber, STAT)

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

        if(.not. Me%TerrainFileNotSpecified)    call OpenAndReadMM5TerrainFile

        call OpenAndReadMM5File

        call BeginEndTime

        call CorrectInitialCondition

        call ComputeVerticalCoordinate

        if(Me%ComputeMeanSeaLevelPressure)      call ComputeMeanSeaLevelPressure

        if(Me%ComputeWindStress)                call ComputeWindShearStress

        if (Me%ComputeWindModulus)              call ComputeWindModulus

        if(Me%ComputeRelativeHumidity)          call ComputeRelativeHumidity
        if(Me%ComputeRelativeHumidity3D)        call ComputeRelativeHumidity3D

        if(Me%ComputePrecipitation)             call ComputePrecipitation

        call OutputFields

        call KillMM5Format


        STAT = SUCCESS_


    end subroutine ConvertMM5Format

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
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR01'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleMM5Format - ERR02'
        end if

        call GetData(Me%TerrainFileName,                                &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'TERRAIN_FILENAME',                 &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR02'
        
        if (iflag == 0)then
            Me%TerrainFileNotSpecified = .true.
            write(*,*)
            write(*,*)'MM5 TERRAIN file was not specified.'
            write(*,*)'    - Mohid grid will not be constructed.'
            write(*,*)'    - HDF5 file will not be viewable using'
            write(*,*)'      Mohid HDF5 Post Processor or Mohid GIS.'
            write(*,*)
        end if

        call GetData(Me%GridFileName,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',             &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR03'

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
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR04'


        call GetData(Me%WriteXYZ,                                       &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'WRITE_XYZ',                        &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR06'

        call GetData(Me%ComputeWindStress,                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'COMPUTE_WINDSTRESS',               &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR07'
        

        call GetData(Me%ComputeWindModulus,                             &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'COMPUTE_WINDMODULUS',              &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR07a'


        call GetData(Me%ComputeRelativeHumidity,                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'COMPUTE_RELATIVE_HUMIDITY',        &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR08'

        call GetData(Me%ComputeRelativeHumidity3D,                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'COMPUTE_RELATIVE_HUMIDITY_3D',     &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR08'


        call GetData(Me%ComputePrecipitation,                           &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'COMPUTE_PRECIPITATION',            &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR08a'

        call GetData(Me%ComputeMeanSeaLevelPressure,                    &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'COMPUTE_MSLP',                     &
                     Default      = ON,                                 &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR08b'

        call GetData(Me%WriteTerrain,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'WRITE_TERRAIN',                    &
                     Default      = OFF,                                &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR09'
        
        call GetData(Me%WriteHeaders,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'WRITE_HEADERS',                    &
                     Default      = ON,                                 &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR09'

        call GetData(Me%StartTime,                                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'START',                            &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR10'


        call GetData(Me%EndTime,                                        &
                     Me%ObjEnterData, iflag1,                           &
                     SearchType   = FromBlock,                          &
                     keyword      = 'END',                              &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR11'

        if (iflag==1 .AND. iflag1==1) Me%TimeWindow = .TRUE.
        
        if (Me%TimeWindow) then

            if (Me%StartTime .GE. Me%EndTime) then
                write (*,*)
                write (*,*) 'START greater or equal than END'
                stop 'ReadOptions - ModuleMM5Format - ERR12'
            endif

        endif

        call GetData(Me%OutputDTInterval,                               &
                     Me%ObjEnterData, iflag1,                           &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUT_DT',                        & 
                     Default      = 0.0,                                &
                     ClientModule = 'ModuleMM5Format',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR11'

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
        logical                             :: NeedsSurfacePressure     = .false.
        logical                             :: NeedsSurfaceTemperature  = .false.
        logical                             :: NeedsPressure3D          = .false.
        logical                             :: NeedsTemperature3D       = .false.
        logical                             :: NeedsWindVelocityX       = .false.
        logical                             :: NeedsWindVelocityY       = .false.
        logical                             :: NeedsSurfaceMixingRatio  = .false.
        logical                             :: NeedsMixingRatio         = .false.
        logical                             :: NeedsWindShearVelocity   = .false.
        logical                             :: NeedsConvectivePCP       = .false.
        logical                             :: NeedsNonConvectivePCP    = .false.
        logical                             :: NeedsTerrain             = .false.


        !Begin-----------------------------------------------------------------
        
        call ExtractBlockFromBlock(Me%ObjEnterData,             &
                                    ClientNumber,               &
                                    '<<BeginFields>>',          &
                                    '<<EndFields>>',            &
                                    BlockFound,                 &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModuleMM5Format - ERR01'
    
    
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
            if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModuleMM5Format - ERR02'

            Me%FieldsToConvert(Count) = PropertyName

            Count = Count + 1

        end do

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModuleMM5Format - ERR03'

        allocate(Me%FieldsToRead(1:100))

        Me%FieldsToRead(:) = null_str

        Me%FieldsToRead(1:nPropertiesToConvert) = Me%FieldsToConvert(1:nPropertiesToConvert)

        Count = nPropertiesToConvert + 1


        if(FieldIsToConvert(trim(GetPropertyName(AtmosphericPressure_))))then
            NeedsPressure3D         = .true. ; NeedsTemperature3D       = .true.
        end if

        if(FieldIsToConvert(trim(GetPropertyName(AtmosphericPressure_))//"_3D"))then
            NeedsTemperature3D      = .true.
        end if

        if(Me%ComputeRelativeHumidity)then
!            NeedsPressure3D         = .true. ; NeedsTemperature3D       = .true.
            NeedsSurfacePressure    = .true. ; NeedsSurfaceTemperature  = .true.
            NeedsSurfaceMixingRatio = .true.
        end if

        if(Me%ComputeRelativeHumidity3D)then
            NeedsPressure3D         = .true. ; NeedsTemperature3D       = .true.            
            NeedsMixingRatio        = .true.
        end if

        if(Me%ComputePrecipitation)then
            NeedsConvectivePCP          = .true.
            NeedsNonConvectivePCP       = .true.
        end if

        if(Me%ComputeMeanSeaLevelPressure)then
            NeedsSurfacePressure        = .true.
            NeedsPressure3D             = .true.  
            NeedsTemperature3D          = .true.
            NeedsTerrain                = .true.
        end if

        if(Me%ComputeWindStress)then
            NeedsWindVelocityX      = .true. ; NeedsWindVelocityY       = .true.
            NeedsWindShearVelocity  = .true. 
        end if

        if (Me%ComputeWindModulus) then 
            NeedsWindVelocityX          = .true.
            NeedsWindVelocityY          = .true.
        endif

        if(AtLeastOne3DFieldIsToConvert())then
            NeedsPressure3D         = .true. ; NeedsTemperature3D       = .true.
            NeedsSurfacePressure    = .true.
        end if

        if(FieldIsToConvert(trim(GetPropertyName(AtmosphericPressure_))))then
            NeedsPressure3D         = .true. ; NeedsTemperature3D       = .true.
        end if

        !MM5 doesn't have Surface temperature initial condition

        if(FieldIsToConvert(trim(GetPropertyName(AirTemperature_))))then
            NeedsTemperature3D       = .true.
        end if

        !allocate Me%FieldsToRead ------------------------------------------------

        if(NeedsPressure3D .and. .not. FieldIsToConvert(trim(GetPropertyName(AtmosphericPressure_))//"_3D"))then
            Me%FieldsToRead(Count) = trim(GetPropertyName(AtmosphericPressure_))//"_3D"
            Count = Count + 1
        end if

        if(NeedsTemperature3D .and. .not. FieldIsToConvert(trim(GetPropertyName(AirTemperature_))//"_3D"))then
            Me%FieldsToRead(Count) = trim(GetPropertyName(AirTemperature_))//"_3D"
            Count = Count + 1
        end if

        if(NeedsSurfacePressure .and. .not. FieldIsToConvert(trim(GetPropertyName(AtmosphericPressure_))))then
            Me%FieldsToRead(Count) = trim(GetPropertyName(AtmosphericPressure_))
            Count = Count + 1
        end if

        if(NeedsSurfaceTemperature .and. .not. FieldIsToConvert(trim(GetPropertyName(AirTemperature_))))then
            Me%FieldsToRead(Count) = trim(GetPropertyName(AirTemperature_))
            Count = Count + 1
        end if

        if(NeedsSurfaceMixingRatio .and. .not. FieldIsToConvert(trim("2-meter mixing ratio")))then
            Me%FieldsToRead(Count) = trim("2-meter mixing ratio")
            Count = Count + 1
        end if

        if(NeedsMixingRatio .and. .not. FieldIsToConvert(trim("mixing ratio")//"3D"))then
            Me%FieldsToRead(Count) = trim("mixing ratio")//"_3D"
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

    
    !--------------------------------------------------------------------------

    
    subroutine OpenAndReadMM5TerrainFile

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: exist

        !Others----------------------------------------------------------------
        integer                                     :: TerrainUnit
        integer                                     :: flag
        integer,            dimension(50,20)        :: bhi
        real,               dimension(20,20)        :: bhr
        character(len=80),  dimension(50,20)        :: bhic
        character(len=80),  dimension(20,20)        :: bhrc
        integer                                     :: nDimensions
        real                                        :: time
        integer,            dimension(4)            :: start_index, end_index, expanded_end_index
        character (len= 4)                          :: staggering
        character (len= 4)                          :: ordering
        character (len=24)                          :: current_date
        character (len= 9)                          :: MM5Name
        character (len=25)                          :: units
        character (len=46)                          :: description
        real, pointer,     dimension(:,:,:,:)       :: DataAux, DataAux2
        integer                                     :: ierr, ier
        character (len=StringLength)                :: GridInformationName
        integer                                     :: MotherDomainID, DomainID
        logical                                     :: IsCoarseDomainExpanded
        integer                                     :: GridOffSetI, GridOffSetJ
        integer                                     :: i,j,k, i2,j2
        

        !Begin-----------------------------------------------------------------
        
        ierr = 0
        IsCoarseDomainExpanded = .false.
        MotherDomainID  = null_int
        DomainID        = null_int
        GridOffSetI     = null_int
        GridOffSetJ     = null_int

        nullify(Me%ConnectionX)
        nullify(Me%ConnectionY)
        nullify(Me%Bathymetry )
        nullify(Me%CenterX    )
        nullify(Me%CenterY    )
        nullify(Me%LandUse    )

        !Verifies if file exists
        inquire(file = Me%TerrainFileName, exist = exist)
        if (.not. exist) then
            write(*,*)'MM5 TERRAIN file does not exist'
            stop 'OpenAndReadMM5TerrainFile - ModuleMM5Format - ERR01'
        endif

        write(*,*)'---------------------------'
        write(*,*)
        write(*,*)'Reading MM5 TERRAIN file...'
        write(*,*)

        call UnitsManager(TerrainUnit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMM5TerrainFile - ModuleMM5Format - ERR02'

        open(Unit   = TerrainUnit,          &
             File   = Me%TerrainFileName,   &
             Form   = 'UNFORMATTED',        &
             STATUS = 'OLD',                &
             Action = 'READ',               &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMM5TerrainFile - ModuleMM5Format - ERR03'

        rewind(TerrainUnit)

        read  (TerrainUnit, IOSTAT = STAT_CALL) flag
        if(STAT_CALL/=0) stop 'OpenAndReadMM5TerrainFile - ModuleMM5Format - ERR03b - Check big-endian compatibility flag'

do0:    do while (ierr == 0)

cd0:        if (flag == 0) then

                read(TerrainUnit,iostat=ier) bhi, bhr, bhic, bhrc 
                if(ier/=0) stop 'OpenAndReadMM5TerrainFile - ModuleMM5Format - ERR04 - Check big-endian compatibility flag'

                !BHI(  7, 1): MAP PROJECTION. 1: LAMBERT CONFORMAL, 2: POLAR STEREOGRAPHIC, 3: MERCATOR
                
                Me%ProjType = bhi(7,1)

                if (Me%ProjType == 1) then
                    Me%ProjType = LAMB_CONF_CONIC_


                elseif (Me%ProjType == 2) then
                    write(*,*) 'Not ready for Polar Stereographic Map Projection'
                    stop 'OpenAndReadMM5TerrainFile - ModuleMM5Format - ERR04a'
#ifdef _USE_PROJ4
                elseif (Me%ProjType == 3) then
                    Me%ProjType = SPHERE_MERCATOR_
#else
                elseif (Me%ProjType == 3) then
                    write(*,*) 'Not ready for Mercator Map Projection'
                    stop 'OpenAndReadMM5TerrainFile - ModuleMM5Format - ERR04b'  
#endif
                endif


                Me%CoarseDomainCenterLat = bhr(2,1)
                Me%CoarseDomainCenterLon = bhr(3,1)
                Me%TrueLatUpper          = bhr(5,1)
                Me%TrueLatLower          = bhr(6,1)

                if (bhi(8,1) == 1) IsCoarseDomainExpanded = .true.
                DomainID        = bhi(13,1)
                MotherDomainID  = bhi(14,1)                
                GridOffSetI     = bhi(11,1)
                GridOffSetJ     = bhi(12,1)

                
           
            elseif (flag == 1) then cd0

                read (TerrainUnit,iostat=ier) nDimensions, start_index, end_index, time,   &
                                              staggering, ordering,current_date,           &
                                              MM5Name, units, description
                
                if(ier/=0) stop 'OpenAndReadMM5TerrainFile - ModuleMM5Format - ERR05'

                

                if (nDimensions == 1) then
                    allocate(DataAux2(end_index(1), 1, 1, 1))
                elseif (nDimensions == 2) then
                    allocate(DataAux2(end_index(1), end_index(2), 1, 1))
                elseif (nDimensions == 3) then
                    allocate(DataAux2(end_index(1), end_index(2), end_index(3), 1))
                endif

                read(TerrainUnit, end = 100, err=100, iostat = ier) DataAux2

                if ((MotherDomainID .eq. DomainID) .and. IsCoarseDomainExpanded) then

                    expanded_end_index = end_index

                    end_index(1) = end_index(1) - 2*GridOffSetI
                    end_index(2) = end_index(2) - 2*GridOffSetJ

                    if (nDimensions == 1) then
                        allocate(DataAux(end_index(1), 1, 1, 1))
                        DataAux = null_real
                        i=1
                        do i2 = GridOffSetI, expanded_end_index(1) - GridOffSetI
                            DataAux(i,1,1,1) = DataAux(i2,1,1,1)
                            i=i+1
                        enddo

                    elseif (nDimensions == 2) then
                        allocate(DataAux(end_index(1),end_index(2), 1, 1))
                        DataAux = null_real

                        i=1
                        j=1

                        do j2 = GridOffSetJ+1, expanded_end_index(2) - GridOffSetJ                    
                        i=1
                        do i2 = GridOffSetI+1, expanded_end_index(1) - GridOffSetI
                            DataAux(i,j,1,1) = DataAux2(i2,j2,1,1)
                            i=i+1
                        enddo
                        j=j+1
                        enddo

                    elseif (nDimensions == 3) then
                        allocate(DataAux(end_index(1), end_index(2), end_index(3), 1))
                        DataAux = null_real
                        
                        do k = 1, end_index(3)
                        i=1
                        j=1
                        do j2 = GridOffSetJ+1, expanded_end_index(2) - GridOffSetJ                    
                        i=1
                        do i2 = GridOffSetI+1, expanded_end_index(1) - GridOffSetI
                            DataAux(i,j,k,1) = DataAux2(i2,j2,k,1)
                            i=i+1
                        enddo
                        j=j+1
                        enddo                        
                        enddo

                    endif
                
                                        
                else
                    DataAux => DataAux2
                endif

                if(Me%WriteTerrain)then

                    call WriteTerrainField(MM5Name, DataAux, units, start_index, end_index)

                end if

                if(CheckGridInformation(MM5Name, GridInformationName))then

                    call FillGridInformation(GridInformationName, DataAux, end_index)

                end if

                deallocate(DataAux)

            end if cd0

            read(TerrainUnit, iostat=ierr) flag

        enddo do0

100     write(*,*)
        write(*,*)'Done.'
        write(*,*)

        if(                  .not. Me%OutputGridNotSpecified) call ConstructGrid

        if(Me%WriteXYZ .and. .not. Me%OutputGridNotSpecified) call WriteCenterGridXYZ

        call WriteGridToHDF5File


    end subroutine OpenAndReadMM5TerrainFile
    
    !----------------------------------------------------------------------

    subroutine FillGridInformation(Name, DataAux, end_index)
        
        !Arguments-------------------------------------------------------------
        character (len= *)                          :: Name
        real, pointer,     dimension(:,:,:,:)       :: DataAux
        integer,           dimension(4)             :: end_index

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: WILB, WIUB, WJLB, WJUB

        !Begin-----------------------------------------------------------------

        Me%Size%ILB = 0           ; Me%WorkSize%ILB = Me%Size%ILB + 1
        Me%Size%IUB = end_index(1); Me%WorkSize%IUB = Me%Size%IUB - 1
        Me%Size%JLB = 0           ; Me%WorkSize%JLB = Me%Size%JLB + 1
        Me%Size%JUB = end_index(2); Me%WorkSize%JUB = Me%Size%JUB - 1

        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 

        select case(Name)

            case('Bathymetry')

                allocate(Me%Bathymetry(ILB:IUB, JLB:JUB))

                Me%Bathymetry(WILB:WIUB, WJLB:WJUB) = DataAux(WILB:WIUB, WJLB:WJUB, 1, 1)

            case('CenterX')

                allocate(Me%CenterX    (ILB:IUB, JLB:JUB))

                Me%CenterX(WILB:WIUB, WJLB:WJUB) = DataAux(WILB:WIUB, WJLB:WJUB, 1, 1)

            case('CenterY')

                allocate(Me%CenterY    (ILB:IUB, JLB:JUB))

                Me%CenterY(WILB:WIUB, WJLB:WJUB) = DataAux(WILB:WIUB, WJLB:WJUB, 1, 1)


            case('ConnectionX')

                allocate(Me%ConnectionX(ILB:IUB, JLB:JUB))

                Me%ConnectionX(WILB:WIUB+1, WJLB:WJUB+1) = DataAux(WILB:WIUB+1, WJLB:WJUB+1,1,1)

            case('ConnectionY')

                allocate(Me%ConnectionY(ILB:IUB, JLB:JUB))

                Me%ConnectionY(WILB:WIUB+1, WJLB:WJUB+1) = DataAux(WILB:WIUB+1, WJLB:WJUB+1,1,1)

            case('LandUse')

                allocate(Me%LandUse(ILB:IUB, JLB:JUB))

                Me%LandUse(WILB:WIUB, WJLB:WJUB) = DataAux(WILB:WIUB, WJLB:WJUB,1,1)


            case default

        end select

    end subroutine FillGridInformation
    
    
    !----------------------------------------------------------------------

    
    subroutine WriteTerrainField(MM5Name, DataAux, units, start_index, end_index)

        !Arguments-------------------------------------------------------------
        character (len= 9)                          :: MM5Name
        character (len=25)                          :: units
        real, pointer,     dimension(:,:,:,:)       :: DataAux
        integer,           dimension(4)             :: start_index, end_index

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        real,    dimension(:,:), pointer            :: Field

        !----------------------------------------------------------------------

        call HDF5SetLimits  (Me%ObjHDF5, start_index(1), end_index(1)-1,&
                             start_index(2), end_index(2)-1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteTerrainField - ModuleMM5Format - ERR01'

        Field => DataAux(:, :, 1, 1)
       
        call HDF5WriteData   (Me%ObjHDF5, "/Terrain", MM5Name, units,       &
                              Array2D =  Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteTerrainField - ModuleMM5Format - ERR02' 

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteTerrainField - ModuleMM5Format - ERR03'


    end subroutine WriteTerrainField


    !--------------------------------------------------------------------------

    
    subroutine OpenAndReadMM5File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: exist

        !Others----------------------------------------------------------------
        integer,            dimension(50,20)        :: bhi
        real,               dimension(20,20)        :: bhr
        character(len=80),  dimension(50,20)        :: bhic
        character(len=80),  dimension(20,20)        :: bhrc
        integer                                     :: flag
        integer                                     :: nDimensions
        real                                        :: time
        integer,            dimension(4)            :: start_index, end_index
        character (len= 4)                          :: staggering
        character (len= 4)                          :: ordering
        character (len=24)                          :: current_date
        character (len= 9)                          :: MM5Name
        character (len=25)                          :: units
        character (len=46)                          :: description
        real, pointer,     dimension(:,:,:,:)       :: DataAux, LatCrs, LonCrs
        logical                                     :: PrintLatCrs = .true.
        logical                                     :: PrintLonCrs = .true.
        integer                                     :: ierr, ier, i, j, k
        logical                                     :: newtime      = .true.
        logical                                     :: TerrainOK    = .false.
        character(len=StringLength)                 :: MohidName, GridInformationName
        type(T_Field), pointer                      :: NewField
        type(T_Date ), pointer                      :: NewDate
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        real                                        :: aux1,aux2

        !Begin-----------------------------------------------------------------
        
        write(*,*)'---------------------------'
        write(*,*)
        write(*,*)'Reading MM5 output file...'

        nullify(NewField      )
        nullify(Me%FirstField )
        nullify(Me%FirstDate  )

        ierr = 0

        !Verifies if file exists
        inquire(file = Me%FileName, exist = exist)
        if (.not. exist) then
            write(*,*)'MM5 output file does not exist'
            stop 'OpenAndReadMM5File - ModuleMM5Format - ERR01'
        endif


        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMM5File - ModuleMM5Format - ERR02'

        open(Unit   = Me%Unit,          &
             File   = Me%FileName,      &
             Form   = 'UNFORMATTED',    &
             STATUS = 'OLD',            &
             Action = 'READ',           &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMM5File - ModuleMM5Format - ERR03'

        
        rewind(Me%Unit)
        read  (Me%Unit, IOSTAT = STAT_CALL) flag

do0:    do while (ierr == 0)

cd0:        if (flag == 0) then

                read(Me%Unit,iostat=ier) bhi, bhr, bhic, bhrc 

                Me%Ptop   = bhr(2,2)
                Me%IfSnow = bhi(16,13)

                if(ier/=0) stop 'OpenAndReadMM5File - ModuleMM5Format - ERR04'

    
                if (Me%WriteHeaders) call WriteGlobalHeaders(bhi, bhr, bhic, bhrc)
!                call printout_big_header(bhi, bhr, bhic, bhrc)
            
            elseif (flag == 1) then cd0
                
                read (Me%Unit,iostat=ier) nDimensions, start_index, end_index, time,   &
                                          staggering, ordering, current_date,          &
                                          MM5Name, units, description

                if(ier/=0) stop 'OpenAndReadMM5File - ModuleMM5Format - ERR05'

                if (newtime) then

                    call SetNewDate(current_date, NewDate)

                    newtime = .false.

                endif

                if     (nDimensions == 1) then
                    allocate(DataAux(end_index(1), 1, 1, 1))
                elseif (nDimensions == 2) then
                    allocate(DataAux(end_index(1), end_index(2), 1, 1))
                elseif (nDimensions == 3) then
                    allocate(DataAux(end_index(1), end_index(2), end_index(3), 1))
                endif

                read(Me%Unit, END = 100, ERR=100) DataAux

                !Temporario - So para saber ao certo as dimensoes do dominio

                if(MM5Name == 'LATITCRS' .and. PrintLatCrs )then
                    allocate(LatCrs(end_index(1), end_index(2), 1, 1))
                    LatCrs = DataAux

                    write (*,'(/,"LatCrs( 1, 1, 1, 1)", F10.4)') LatCrs(1,1,1,1)
                    write (*,'(/,"LatCrs(", I2,",1,1,1)", F10.4)') end_index(1), LatCrs(end_index(1),1,1,1)
                    write (*,'(/,"LatCrs( 1,", I2, ", 1, 1)", F10.4)') end_index(2), LatCrs(1,end_index(2),1,1)
                    write (*,'(/,"LatCrs(",I2,",", I2, ", 1, 1)", F10.4)') end_index(1), end_index(2), &
                                                                           LatCrs(end_index(1),end_index(2),1,1)
                    write (*,*)

                    PrintLatCrs = .false.

                endif

                
                if(MM5Name == 'LONGICRS' .and. PrintLonCrs)then
                    allocate(LonCrs(end_index(1), end_index(2), 1, 1))
                    LonCrs = DataAux

                    write (*,'(/,"LonCrs( 1, 1, 1, 1)", F10.4)') LonCrs(1,1,1,1)
                    write (*,'(/,"LonCrs(", I2,",1,1,1)", F10.4)') end_index(1), LonCrs(end_index(1),1,1,1)
                    write (*,'(/,"LonCrs( 1,", I2, ", 1, 1)", F10.4)') end_index(2), LonCrs(1,end_index(2),1,1)
                    write (*,'(/,"LonCrs(",I2,",", I2, ", 1, 1)", F10.4)') end_index(1), end_index(2), &
                                                                           LonCrs(end_index(1),end_index(2),1,1)
                    write (*,*)

                    PrintLonCrs = .false.

                endif

                if(nDimensions == 1)then

                    if(MM5Name == 'SIGMAH')then

                        allocate(Me%Sigma (start_index(1):end_index(1)))
                        Me%Sigma = DataAux(start_index(1):end_index(1), 1, 1, 1)

                    end if

                end if

                if(.not. TerrainOK .and. (Me%WriteTerrain .and. Me%TerrainFileNotSpecified))then

                    if(IsTerrainField(MM5Name))then

                        call WriteTerrainField(MM5Name, DataAux, units, start_index, end_index)

                    end if

                end if

                if(.not. TerrainOK .and. Me%TerrainFileNotSpecified)then

                    if(CheckGridInformation(MM5Name, GridInformationName))then

                        call FillGridInformation(GridInformationName, DataAux, end_index)

                    end if

                end if                    

                if(FieldIsToRead(MM5Name, MohidName))then

                    call AddField(Me%FirstField, NewField)

                    call AttributeDimensions(NewField, nDimensions, start_index, end_index, staggering, ordering)

                    call SetNewFieldAttributes(Field    = NewField,                &
                                               Name     = trim(MohidName),         &
                                               Units    = trim(Units),             &
                                               Date     = NewDate%Date,            &
                                               Convert  = FieldIsToConvert(MohidName))

                    WILB = NewField%WorkSize%ILB; ILB  = NewField%Size%ILB
                    WIUB = NewField%WorkSize%IUB; IUB  = NewField%Size%IUB
                    WJLB = NewField%WorkSize%JLB; JLB  = NewField%Size%JLB
                    WJUB = NewField%WorkSize%JUB; JUB  = NewField%Size%JUB
                    WKLB = NewField%WorkSize%KLB; KLB  = NewField%Size%KLB
                    WKUB = NewField%WorkSize%KUB; KUB  = NewField%Size%KUB


                    if (IUB .ne. Me%Size%IUB) then
                        write (*,*) 
                        write (*,*) 'Error in field ', trim(MohidName)
                        write (*,*) 'Field IUB = ', IUB
                        write (*,*) 'Grid IUB  = ', Me%Size%IUB
                        stop 'OpenAndReadMM5File - ModuleMM5Format - ERR06'
                    endif

                    if ((nDimensions .gt. 1) .and. (JUB .ne. Me%Size%JUB)) then
                        write (*,*) 
                        write (*,*) 'Error in field ', trim(MohidName)
                        write (*,*) 'Field JUB = ', JUB
                        write (*,*) 'Grid JUB  = ', Me%Size%JUB
                        stop 'OpenAndReadMM5File - ModuleMM5Format - ERR07'
                    endif

                    if(Me%Size%KUB     .ne. null_int) then
                        if ((nDimensions .gt. 2) .and. (KUB .ne. Me%Size%KUB)) then
                            write (*,*) 
                            write (*,*) 'Error in field ', trim(MohidName)
                            write (*,*) 'Field KUB = ', KUB
                            write (*,*) 'Grid KUB  = ', Me%Size%KUB
                            stop 'OpenAndReadMM5File - ModuleMM5Format - ERR08'
                        endif
                    endif

                    select case(NewField%nDimensions)

                        case(2)
                            
                            allocate(NewField%Values2D(ILB:IUB, JLB:JUB))

                            if(NewField%GridLocation == DotGrid)then

                                NewField%Values2D(WILB:WIUB, WJLB:WJUB) =               &
                                (DataAux(WILB    :WIUB    , WJLB    :WJUB    , 1, 1) +  &
                                 DataAux(WILB + 1:WIUB + 1, WJLB + 1:WJUB + 1, 1, 1) +  &
                                 DataAux(WILB + 1:WIUB + 1, WJLB    :WJUB    , 1, 1) +  &
                                 DataAux(WILB    :WIUB    , WJLB + 1:WJUB + 1, 1, 1))/ 4.

                            elseif(NewField%GridLocation == CrossGrid)then

                                NewField%Values2D(WILB:WIUB,WJLB:WJUB) = DataAux(WILB:WIUB,WJLB:WJUB, 1, 1)

                            end if

                        case(3)

                            if(Me%WorkSize%KLB .eq. null_int) Me%WorkSize%KLB = NewField%WorkSize%KLB
                            if(Me%WorkSize%KUB .eq. null_int) Me%WorkSize%KUB = NewField%WorkSize%KUB
                            if(Me%Size%KLB     .eq. null_int) Me%Size%KLB     = NewField%Size%KLB
                            if(Me%Size%KUB     .eq. null_int) Me%Size%KUB     = NewField%Size%KUB

                            allocate(NewField%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))


                            if(NewField%GridLocation == DotGrid)then
                               
                                if (NewField%GridLocationZ == FullLevels) then

                                    do j = WJLB, WJUB
                                    do i = WILB, WIUB
                                    do k = WKLB, WKUB

                                        aux1 = (DataAux(i,   j, WKUB + 1 - k, 1) + DataAux(i+1, j+1, WKUB + 1 - k, 1) +  &
                                                DataAux(i+1, j, WKUB + 1 - k, 1) + DataAux(i  , j  , WKUB + 1 - k, 1))/ 4.


                                        aux2 = (DataAux(i,   j, WKUB + 2 - k, 1) + DataAux(i+1, j+1, WKUB + 2 - k, 1) +  &
                                                DataAux(i+1, j, WKUB + 2 - k, 1) + DataAux(i  , j  , WKUB + 2 - k, 1))/ 4.


                                        NewField%Values3D(i, j, k) = (aux1 + aux2) / 2.

                                    enddo
                                    enddo
                                    enddo

                                else 

                                    do j = WJLB, WJUB
                                    do i = WILB, WIUB
                                    do k = WKLB, WKUB

                                        NewField%Values3D(i, j, k) =                                              &
                                        (DataAux(i,   j, WKUB + 1 - k, 1) + DataAux(i+1, j+1, WKUB + 1 - k, 1) +  &
                                         DataAux(i+1, j, WKUB + 1 - k, 1) + DataAux(i  , j  , WKUB + 1 - k, 1))/ 4.

                                    enddo
                                    enddo
                                    enddo

                                endif

                            elseif(NewField%GridLocation == CrossGrid)then
    
                                if (NewField%GridLocationZ == FullLevels) then

                                    do j = WJLB, WJUB
                                    do i = WILB, WIUB
                                    do k = WKLB, WKUB

                                        NewField%Values3D(i, j, k) = (DataAux(i,j, WKUB + 1 - k, 1) + &
                                                                      DataAux(i,j, WKUB + 2 - k, 1)) / 2.

                                    enddo
                                    enddo
                                    enddo

                                else


                                    do j = WJLB, WJUB
                                    do i = WILB, WIUB
                                    do k = WKLB, WKUB

                                        NewField%Values3D(i, j, k) = DataAux(i,j, WKUB + 1 - k, 1)

                                    enddo
                                    enddo
                                    enddo

                                end if

                            end if
                            

                        case default

                    end select

                end if

                deallocate(DataAux)

            elseif (flag == 2) then cd0

                newtime     = .true.
                TerrainOK   = .true.

            else cd0

                stop

            endif cd0

            read(Me%Unit, iostat=ierr) flag

        enddo do0

        if(Me%WriteXYZ .and. .not. Me%OutputGridNotSpecified) call WriteCenterGridXYZ

100     write(*,*)'Finished reading MM5 output file.'

        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMM5File - ModuleMM5Format - ERR05'

    end subroutine OpenAndReadMM5File
    
    
    !------------------------------------------------------------------------

    
    subroutine SetNewDate(current_date, NewDate)

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

    !------------------------------------------------------------------------

    subroutine CorrectInitialCondition


        type(T_Field), pointer                  :: Field
        real, dimension(:,:  ), pointer         :: SurfaceTemperature

        if (Me%StartTime == Me%FirstDate%Date) then

            !Find T3D

            Field => Me%FirstField        

            do while(associated(Field))
        
                if(Field%Name == trim(GetPropertyName(AirTemperature_))//"_3D"      .and. &
                   Field%Date == Me%FirstDate%Date)then
                    
                    SurfaceTemperature  => Field%Values3D(:,:,1)
            
                end if

                Field => Field%Next
            
            end do

            !Correct surface temperature

            Field => Me%FirstField

            do while(associated(Field))
        
                if(Field%Name == trim(GetPropertyName(AirTemperature_))      .and. &
                   Field%Date == Me%FirstDate%Date)then
                    
                    Field%Values2D = SurfaceTemperature
            
                end if
    
                Field => Field%Next

            end do

        endif


    end subroutine CorrectInitialCondition

    !------------------------------------------------------------------------

    subroutine AttributeDimensions(Field, nDimensions, start_index, end_index, staggering, ordering)
        
        !Arguments-------------------------------------------------------------
        type(T_Field), pointer              :: Field
        integer                             :: nDimensions
        integer, dimension(4)               :: start_index
        integer, dimension(4)               :: end_index
        character (len= 4)                  :: staggering, ordering

        !Begin-----------------------------------------------------------------

        Field%WorkSize%ILB   = start_index(1)
        Field%WorkSize%IUB   = end_index  (1) - 1
        Field%WorkSize%JLB   = start_index(2)
        Field%WorkSize%JUB   = end_index  (2) - 1
        Field%WorkSize%KLB   = start_index(3)

        if (trim(ordering) == 'YXW') then
            Field%WorkSize%KUB   = end_index  (3) - 1
            Field%GridLocationZ = FullLevels
        else 
            Field%WorkSize%KUB   = end_index  (3)
            Field%GridLocationZ = HalfLevels  
        endif 

        Field%Size%ILB       = Field%WorkSize%ILB - 1
        Field%Size%IUB       = Field%WorkSize%IUB + 1
        Field%Size%JLB       = Field%WorkSize%JLB - 1
        Field%Size%JUB       = Field%WorkSize%JUB + 1
        Field%Size%KLB       = Field%WorkSize%KLB - 1
        Field%Size%KUB       = Field%WorkSize%KUB + 1

        Field%nDimensions    = nDimensions
        
        if    (trim(staggering) == 'D')then
            Field%GridLocation  = DotGrid
        elseif(trim(staggering) == 'C')then
            Field%GridLocation  = CrossGrid
        end if

    end subroutine AttributeDimensions

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
        if(STAT_CALL .ne. SUCCESS_)stop 'WriteCenterGridXYZ - ModuleMM5Format - ERR01'

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
        if(STAT_CALL .ne. SUCCESS_)stop 'WriteCenterGridXYZ - ModuleMM5Format - ERR02'


    end subroutine WriteCenterGridXYZ

    
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
                             COMENT1        = 'Grid Data created from MM5 TERRAIN file',    &
                             COMENT2        = trim(Me%TerrainFileName),                     &
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
                             Datum          = SPHERE_DATUM,                                 &
                             ProjType       = Me%ProjType,                                  & 
                             SP1            = Me%TrueLatLower,                              &
                             SP2            = Me%TrueLatUpper,                              &
                             STAT           = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructGrid - ModuleMM5Format - ERR01'

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleMM5Format - ERR02'


    end subroutine ConstructGrid

    
    !------------------------------------------------------------------------

    subroutine ConvertToMohidUnits2D(Field)

        !Local-----------------------------------------------------------------
        type (T_Field), pointer                   :: Field
        
        !Begin-----------------------------------------------------------------

        select case(trim(Field%Units))


            case('K')

                Field%Units     = 'ºC'
                Field%Values2D  = Field%Values2D - AbsoluteZero
            
            case('m s{-1}')

                Field%Units     = 'm/s'

            case('cm')

                Field%Units     = 'm/s'
                Field%Values2D  = Field%Values2D / 3600. / 100.

            case('W/m^2')

                Field%Units     = 'W/m2'

            case default



        end select



    end subroutine ConvertToMohidUnits2D

    !------------------------------------------------------------------------

    subroutine ConvertToMohidUnits3D(Field)

        !Local-----------------------------------------------------------------
        type (T_Field), pointer                   :: Field
        
        !Begin-----------------------------------------------------------------

        select case(trim(Field%Units))

            case('K')

                Field%Units     = 'ºC'
                Field%Values3D  = Field%Values3D - AbsoluteZero
            
            case('m s{-1}')

                Field%Units     = 'm/s'

            case('cm')

                Field%Units     = 'm/s'
                Field%Values3D  = Field%Values3D / 3600. / 100.

            case('W/m^2')

                Field%Units     = 'W/m2'

            case default

        end select

    end subroutine ConvertToMohidUnits3D

    !------------------------------------------------------------------------

    subroutine ComputeVerticalCoordinate
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer                  :: VerticalZ
        real, dimension(:,:  ), pointer         :: ReferenceSurfacePressure
        real, dimension(:,:,:), pointer         :: Pressure3D
        real, dimension(:,:,:), pointer         :: Temperature
        logical                                 :: Pref_OK , Pressure_OK, Temperature_OK            
        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: i,j,k
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

        Temperature_OK  = .false.
        Pref_OK         = .false.
        Pressure_OK     = .false.


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
                    
                    Temperature                 => Field%Values3D
                    Temperature_OK              = .true.

                end if


                if(Pref_OK .and. Pressure_OK .and. Temperature_OK)then

                    call AddField(Me%FirstField, VerticalZ)

                    call SetNewFieldAttributes(Field         =  VerticalZ,          &
                                               Name          = 'VerticalZ',         &
                                               Units         = 'm',                 &
                                               Date          = CurrentDate%Date,    &
                                               WorkSize      = Me%WorkSize,         &
                                               nDimensions   = 3,                   &
                                               Convert       = .true.)

                    allocate(VerticalZ%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))
                    VerticalZ%Values3D = 0.
                    
                    do j = WJLB, WJUB
                    do i = WILB, WIUB
                    do k = WKLB, WKUB
                        
                        !Pressure(i,j,k) = Sigma(k) * PSTARCRS(i,j) + PTOP + PP(i,j,k)
                        Pressure3D(i,j,k) = Me%Sigma(WKUB + 1 - k) * ReferenceSurfacePressure(i,j) + &
                                            Me%PTop + Pressure3D(i,j,k)

                    enddo
                    enddo
                    enddo


                    Pressure3D(ILB:IUB, JLB:JUB,WKUB + 1) = Me%Ptop

                    VerticalZ%Values3D(WILB:WIUB, WJLB:WJUB, 0) = - Me%Bathymetry(WILB:WIUB, WJLB:WJUB)

                    do j = WJLB, WJUB
                    do i = WILB, WIUB
                    do k = WKLB+1, WKUB+1

                        VerticalZ%Values3D(i,j,k-1) = VerticalZ%Values3D(i,j,k-2) + &
                                                      (Pressure3D(i,j,k) - Pressure3D(i,j,k-1))/ &
                                                      (-Gravity * Pressure3D(i,j,k-1)     / &
                                                      (Perfect_gas_R * Temperature(i,j,k-1)))

                    enddo
                    enddo
                    enddo


                    ReferenceSurfacePressure(WILB:WIUB, WJLB:WJUB) =                            &
                                                Pressure3D(WILB:WIUB, WJLB:WJUB, WKLB) +        &
                                                (1. - Me%Sigma(WKUB + 1 - WKLB)) *              &
                                                ReferenceSurfacePressure(WILB:WIUB, WJLB:WJUB)

                    !ReferenceSurfacePressure(WILB:WIUB, WJLB:WJUB) = Pressure3D(WILB:WIUB, WJLB:WJUB, WKLB) 


                    Temperature_OK  = .false.
                    Pref_OK         = .false.
                    Pressure_OK     = .false.

                    nullify(Pressure3D, ReferenceSurfacePressure) 


               end if

                Field => Field%Next

            end do


            CurrentDate => CurrentDate%Next

        end do 


    end subroutine ComputeVerticalCoordinate
    
    
    !------------------------------------------------------------------------

    subroutine ComputeMeanSeaLevelPressure 

    
        !Local-----------------------------------------------------------------        
        type(T_Field), pointer                  :: MeanSeaLevelPressure
        real, dimension(:,:  ), pointer         :: ReferenceSurfacePressure, Terrain
        real, dimension(:,:  ), pointer         :: PL, XKLEV, TS, T0       
        real, dimension(:,:,:), pointer         :: Pressure3D
        real, dimension(:,:,:), pointer         :: Temperature
        logical                                 :: Pref_OK , Pressure_OK, Temperature_OK, Terrain_OK
        logical                                 :: Restart, L1,L2,L3        
        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB, KUPTO
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: i,j,k
        type(T_Field), pointer                  :: Field
        type(T_Date), pointer                   :: CurrentDate
        real                                    :: PCONST = 10000.    !Pa (find layer 100 hPa above surface)
        real                                    :: GAMMA = 6.5E-3   !k/m (neutral stability)
        real                                    :: TC = 273.15+17.5 !critical temperature
        integer                                 :: KLO, KHI
        real                                    :: PLO, PHI, TLO, THI, PSFC
        real                                    :: TL, TBAR, HL, XTERM, Tmed



        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Computing mean sea level pressure ...'

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

        XTERM = GAMMA * Perfect_gas_R / Gravity        
!        Terrain => Me%Bathymetry(WILB:WIUB, WJLB:WJUB)  !Nao pode ser pq o terreno muda com 2-way nesting

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
                    
                    Temperature                 => Field%Values3D
                    Temperature_OK              = .true.

                end if

                if(Field%Name == trim('terrain')                                    .and. &
                   Field%Date == CurrentDate%Date)then
                    
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

                    allocate(MeanSeaLevelPressure%Values2D(ILB:IUB, JLB:JUB))
                    MeanSeaLevelPressure%Values2D = 0.

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

                        if (Restart) exit

                        enddo
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
                        
                        TLO  = Temperature(i,j,KLO)                                               
                        THI  = Temperature(i,j,KHI)                                               

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
                        HL   = Terrain(i,j) - Perfect_Gas_R / Gravity * ALOG(PL(i,j)/PSFC)*TBAR                 
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
                                                              (Perfect_Gas_R * (TS(i,j)+T0(i,j))))
                    enddo
                    enddo

                    Temperature_OK  = .false.
                    Pref_OK         = .false.
                    Pressure_OK     = .false.
                    Terrain_OK      = .false.             

                    nullify(ReferenceSurfacePressure, Pressure3D, Temperature, Terrain)

                    deallocate(PL, XKLEV, TS, T0)
                    
                    exit

               end if


                Field => Field%Next

            end do


            CurrentDate => CurrentDate%Next

        end do 

        nullify(Terrain)


    end subroutine ComputeMeanSeaLevelPressure

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

                    ILB = Field%Size%ILB; WILB = Field%WorkSize%ILB 
                    IUB = Field%Size%IUB; WIUB = Field%WorkSize%IUB 
                    JLB = Field%Size%JLB; WJLB = Field%WorkSize%JLB 
                    JUB = Field%Size%JUB; WJUB = Field%WorkSize%JUB 
                    
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
        if (STAT_CALL /= SUCCESS_)stop 'WriteMapping - ModuleMM5Format - ERR10'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3D", "-",    &
                              Array3D = WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteMapping - ModuleMM5Format - ERR20'

        deallocate(WaterPoints3D)
        nullify   (WaterPoints3D)

    end subroutine WriteMapping
    
    
    !------------------------------------------------------------------------



    
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
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR30'


            call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                 "Time", "YYYY/MM/DD HH:MM:SS",             &
                                 Array1D = TimePtr,                         &
                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR40'


            Field => Me%FirstField

            do while(associated(Field))

                if(Field%Convert)then
            
                    if(Field%Date == CurrentDate%Date)then
                
                        Field%OutputNumber = OutputNumber

                        if(Field%nDimensions == 2)then

                            call ConvertToMohidUnits2D(Field)

                            call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                               Field%WorkSize%JLB, Field%WorkSize%JUB, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR50'

                            call HDF5WriteData(Me%ObjHDF5,                                      &
                                               "/Results/"//Field%Name,                         &
                                               Field%Name,                                      &
                                               Field%Units,                                     &
                                               Array2D      = Field%Values2D,                   &
                                               OutputNumber = Field%OutputNumber,               &
                                               STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR60'


                        elseif(Field%nDimensions == 3)then

                            if(Field%Name == 'VerticalZ')then

                                call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                                   Field%WorkSize%JLB, Field%WorkSize%JUB,            &
                                                   Field%WorkSize%KLB - 1, Field%WorkSize%KUB, STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR70'

                                !Invert vertical axis
                                Field%Values3D = -1. * Field%Values3D

                                call HDF5WriteData(Me%ObjHDF5,                                  &
                                               "/Grid/"//Field%Name,                            &
                                               'Vertical',                                      &
                                               Field%Units,                                     &
                                               Array3D      = Field%Values3D,                   &
                                               OutputNumber = Field%OutputNumber,               &
                                               STAT         = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR80'

                            else
                    
                                call ConvertToMohidUnits3D(Field)

                                call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                                   Field%WorkSize%JLB, Field%WorkSize%JUB,            &
                                                   Field%WorkSize%KLB, Field%WorkSize%KUB, STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR90'

                                call HDF5WriteData(Me%ObjHDF5,                                      &
                                                   "/Results3D/"//Field%Name,                       &
                                                   Field%Name,                                      &
                                                   Field%Units,                                     &
                                                   Array3D      = Field%Values3D,                   &
                                                   OutputNumber = Field%OutputNumber,               &
                                                   STAT         = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR100'

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
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR110'



    end subroutine OutputFields


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
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR01'
        

    end subroutine Open_HDF5_OutPut_File


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
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleMM5Format - ERR02'

        Me%Bathymetry = -1 * Me%Bathymetry

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleMM5Format - ERR03'

        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleMM5Format - ERR02a'

        
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "LandUse", "-",       &
                              Array2D =  Me%LandUse, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleMM5Format - ERR03a'

        if(Me%OutputGridNotSpecified)then

            call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB+1,&
                                 Me%WorkSize%JLB, Me%WorkSize%JUB+1, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleMM5Format - ERR04'


            call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionX", "-",       &
                                  Array2D =  Me%ConnectionX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleMM5Format - ERR05'

            call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionY", "-",       &
                                  Array2D =  Me%ConnectionY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleMM5Format - ERR06'

        else

            call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleMM5Format - ERR07'

        end if
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, &
                             Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleMM5Format - ERR08'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",    &
                              Array2D = WaterPoints2D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleMM5Format - ERR09'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleMM5Format - ERR10'

        deallocate(WaterPoints2D)
        nullify   (WaterPoints2D)

    end subroutine WriteGridToHDF5File


    !--------------------------------------------------------------------------


    logical function IsTerrainField(FieldName)

        !Arguments-----------------------------------------------------------
        character(Len=*)                :: FieldName
        
        !Begin-----------------------------------------------------------------

        IsTerrainField = .false.

        select case(trim(FieldName))

            case('TERRAIN')
            
                IsTerrainField = .true.

            case('LONGICRS')
            
                IsTerrainField = .true.

            case('LATITCRS')
            
                IsTerrainField = .true.

            case('LATITDOT')
                
                IsTerrainField = .true.

            case('LONGIDOT')

                IsTerrainField = .true.
            
            case('LAND USE')

                IsTerrainField = .true.

            case('MAPFACCR')

                IsTerrainField = .true.

            case('MAPFACDT')

                IsTerrainField = .true.

            case('CORIOLIS')

                IsTerrainField = .true.
                
            case default

        end select 


    end function IsTerrainField
    
    !--------------------------------------------------------------------------


    logical function CheckGridInformation(MM5Name, GridInformationName)
        
        !Arguments-----------------------------------------------------------
        character(Len=*)                :: MM5Name
        character(Len=StringLength)     :: GridInformationName
        
        !Begin-----------------------------------------------------------------
        
        CheckGridInformation = .false.

        select case(trim(MM5Name))

            case('TERRAIN')
            
                GridInformationName  = 'Bathymetry'
                CheckGridInformation = .true.

            case('LONGICRS')
            
                GridInformationName  = 'CenterX'
                CheckGridInformation = .true.

            case('LATITCRS')
            
                GridInformationName  = 'CenterY'
                CheckGridInformation = .true.

            case('LATITDOT')
                
                GridInformationName  = 'ConnectionY'
                CheckGridInformation = .true.

            case('LONGIDOT')

                GridInformationName  = 'ConnectionX'
                CheckGridInformation = .true.


            case('LAND USE')

                GridInformationName  = 'LandUse'
                CheckGridInformation = .true.

            case default

        end select 
    
    end function CheckGridInformation
    
    !--------------------------------------------------------------------------

    logical function FieldIsToRead(MM5Name, MohidName)
        
        !Arguments-----------------------------------------------------------
        character(Len=*)                :: MM5Name
        character(Len=StringLength)     :: MohidName
        integer                         :: i
        
        !Begin-----------------------------------------------------------------

        FieldIsToRead = .false.

        select case(trim(MM5Name))

            case('U10')

                MohidName = GetPropertyName(WindVelocityX_)

            case('V10')

                MohidName = GetPropertyName(WindVelocityY_)

            case('UST')

                MohidName = GetPropertyName(WindShearVelocity_)

            case('LWDOWN')

                MohidName = GetPropertyName(DownwardLongWaveRadiation_)

            case('SWDOWN')

                MohidName = GetPropertyName(SolarRadiation_)

            case('LHFLUX')

                MohidName = GetPropertyName(LatentHeat_)

            case('SHFLUX')

                MohidName = GetPropertyName(SensibleHeat_)

            case('SWOUT')

                MohidName = "top outgoing shortwave radiation"

            case('LWOUT')

                MohidName = GetPropertyName(UpwardLongWaveRadiation_)

            case('SOIL T 1')

                MohidName = "soil temperature layer 1"

            case('SOIL T 2')

                MohidName = "soil temperature layer 2"

            case('SOIL T 3')

                MohidName = "soil temperature layer 3"
            
            case('SOIL T 4')

                MohidName = "soil temperature layer 4"
            
            case('SOIL T 5')

                MohidName = "soil temperature layer 5"

            case('SOIL T 6')

                MohidName = "soil temperature layer 6"

            case('Q2')

                MohidName = "2-meter mixing ratio"
 

            case('TSEASFC')

                MohidName = "sea water temperature"

            case('PBL HGT')

                MohidName = "PBL height"

            case('PBL REGIME')

                MohidName = "PBL regime"


            !case('ALBD')
                
            !    MohidName = GetPropertyName(Albedo_)
            !    FieldIsToRead = .true.


            case('RAIN CON')      

                MohidName = AccConvPrecipitation


            case('RAIN NON')      

                MohidName = AccNonConvPrecipitation


            case('GROUND T')

                MohidName = "ground temperature"

            case('RES TEMP')

                MohidName = "infinite reservoir slab temperature"


            case('PSTARCRS')  
                !This is the reference pressure at the surface
                !which is converted to surface pressure
                MohidName = GetPropertyName(AtmosphericPressure_)

            case('T2') !it used to be'GROUND T'. T2 is temperature at 2 meters altitude

                MohidName = GetPropertyName(AirTemperature_)
            
            case('U')

                MohidName = GetPropertyName(WindVelocityX_)
                MohidName = trim(MohidName)//"_3D"

            case('V')

                MohidName = GetPropertyName(WindVelocityY_)
                MohidName = trim(MohidName)//"_3D"

            case('T')

                MohidName = GetPropertyName(AirTemperature_)
                MohidName = trim(MohidName)//"_3D"
            
            case('PP')
                !This is the pressure perturbation field which is
                !converted to pressure values
                MohidName = GetPropertyName(AtmosphericPressure_)
                MohidName = trim(MohidName)//"_3D"

            case('W')

                MohidName = 'wind velocity Z'
                MohidName = trim(MohidName)//"_3D"

            case('Q')
                
                MohidName = 'mixing ratio'
                MohidName = trim(MohidName)//"_3D"
            
            case('CLW')
                
                MohidName = 'cloud water mixing ratio'
                MohidName = trim(MohidName)//"_3D"

            case('RNW')
                
                MohidName = 'rain water mixing ratio'
                MohidName = trim(MohidName)//"_3D"

            case('ICE')
                
                MohidName = 'cloud ice mixing ratio'
                MohidName = trim(MohidName)//"_3D"

            case('SNOW')
                
                MohidName = 'snow mixing ratio'
                MohidName = trim(MohidName)//"_3D"

            case('RAD TEND')
                
                MohidName = 'atmospheric radiation tendency'
                MohidName = trim(MohidName)//"_3D"

            case('TERRAIN')

                MohidName = 'terrain'

!            case('LAND USE')

!                MohidName = 'land use'
            case default
                
                MohidName = ''

        end select

        do i = 1, size(Me%FieldsToRead)

            if(trim(MohidName) == trim(Me%FieldsToRead(i)))then
                FieldIsToRead = .true.
            end if

        enddo


    end function FieldIsToRead
    
    
    !--------------------------------------------------------------------------

    
    subroutine KillMM5Format
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        if(.not. Me%OutputGridNotSpecified)then
            call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillMM5Format - ModuleMM5Format - ERR01'

            deallocate(Me%Bathymetry )
            deallocate(Me%CenterX    )
            deallocate(Me%ConnectionX)
            deallocate(Me%CenterY    )
            deallocate(Me%ConnectionY)

        end if

        
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillMM5Format - ModuleMM5Format - ERR02'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillMM5Format - ModuleMM5Format - ERR03'



        deallocate(Me%FirstField)
        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillMM5Format

    !--------------------------------------------------------------------------

    subroutine printout_big_header(bhi, bhr, bhic, bhrc)


        !Arguments-----------------------------------------------------------
        integer,            dimension(50,20)    :: bhi
        real,               dimension(20,20)    :: bhr
        character(len=80),  dimension(50,20)    :: bhic
        character(len=80),  dimension(20,20)    :: bhrc
        
        !Local-----------------------------------------------------------------
        integer                                 :: i, j, v3j
        
        !Begin-----------------------------------------------------------------

        write(*,'(/)')
        v3j = bhi(1,1)
        
        if (bhi(1,1) == 11) v3j = v3j+5
        
        do j = 1, v3j

            if (j < 8 .or. j>10) then

                if (j == 1)     write(*, '(  "TERRAIN Portion of big header:"         )')
                if (j == 2)     write(*, '(/,"REGRID Portion of big header:"          )')
                if (j == 3)     write(*, '(/,"RAWINS Portion of big header:"          )')
                if (j == 4)     write(*, '(/,"SFC RAWINS Portion of big header:"      )')
                if (j == 5)     write(*, '(/,"INTERP Portion of big header:"          )')
                if (j == 11)    write(*, '(/,"MM5 Portion of big header:"             )')
                if (j == 6)     write(*, '(/,"MM5 Substrate Temp File big header:"    )')
                if (j == 7)     write(*, '(/,"MM5 Boundary File big header:"          )')
                if (j == 8)     write(*, '(/,"Interpolated MM5 Portion of big header:")')
                                write(*, '(/,"***Integers:"/)                          ')

                
                do i = 1, size(bhi,1)
                
                    if (bhi(i,j) /= -999) then
                        write(*,'("BHI(",I3,",",I2,"):",I8," : ",A)')i, j, bhi(i,j),trim(bhic(i,j))
                    endif

                enddo

                write(*,'(/,"***Floats:"/)')
        
                do i = 1, size(bhr,1)
                    if (bhr(i,j) /= -999.) then
                    write(*,'("BHR(",I3,",",I2,"):",F9.2," : ",A)')&
                    i, j, bhr(i,j),trim(bhrc(i,j))
                    endif
                enddo

                write(*,'(/)')
            endif

        enddo

    end subroutine printout_big_header

  
    !--------------------------------------------------------------------------

    subroutine WriteGlobalHeaders(bhi, bhr, bhic, bhrc)


        !Arguments-----------------------------------------------------------
        integer,            dimension(50,20)    :: bhi
        real,               dimension(20,20)    :: bhr
        character(len=80),  dimension(50,20)    :: bhic
        character(len=80),  dimension(20,20)    :: bhrc
        
        !Local-----------------------------------------------------------------
        integer                                 :: i, j, v3j
        integer                                 :: STAT_CALL
        character(len=StringLength)             :: GroupName
        integer                                 :: IMOVE
        
        !Begin-----------------------------------------------------------------

        v3j = bhi(1,1)
        IMOVE = bhi(4,14)
        
        if (bhi(1,1) == 11) v3j = v3j+5
        
do1:    do j = 1, v3j

if1:        if (j < 8 .or. j>10) then

                if (j == 1)     GroupName = "HeaderTerrain"
                if (j == 2)     cycle !write(*, '(/,"REGRID Portion of big header:"          )')
                if (j == 3)     cycle !write(*, '(/,"RAWINS Portion of big header:"          )')
                if (j == 4)     cycle !write(*, '(/,"SFC RAWINS Portion of big header:"      )')
                if (j == 5)     cycle !write(*, '(/,"INTERP Portion of big header:"          )')
                if (j == 11)    GroupName = "HeaderMM5" 
                if (j == 6)     cycle !write(*, '(/,"MM5 Substrate Temp File big header:"    )')
                if (j == 7)     cycle !write(*, '(/,"MM5 Boundary File big header:"          )')
                if (j == 8)     cycle !write(*, '(/,"Interpolated MM5 Portion of big header:")')

                
                do i = 1, size(bhi,1)
                    if (bhi(i,j) /= -999) then

                        if (IMOVE == 0) then
                            if (j == 14 .and. i >= 5  .and. i <= 35) cycle
                        endif

                        call HDF5WriteGlobalAttribute(HDF5ID        = Me%ObjHDF5,       &
                                                      GroupName     = GroupName,        &
                                                      AttributeName = trim(bhic(i,j)),  &
                                                      Att_Int       = bhi(i,j),         &
                                                      STAT          = STAT_CALL)
        
                        if (STAT_CALL /= 0) stop 'WriteTerrainHeader - ModuleMM5Format - ERR01'                    

        
                    endif

                enddo

           
                do i = 1, size(bhr,1)
                    if (bhr(i,j) /= -999.) then
                        call HDF5WriteGlobalAttribute(HDF5ID        = Me%ObjHDF5,       &
                                                      GroupName     = GroupName,        &
                                                      AttributeName = trim(bhrc(i,j)),  &
                                                      Att_Real      = bhr(i,j),         &
                                                      STAT          = STAT_CALL)

                        if (STAT_CALL /= 0) stop 'WriteTerrainHeader - ModuleMM5Format - ERR02'                    
                 
                    endif
                enddo

            endif if1
        enddo do1


    end subroutine WriteGlobalHeaders

    !--------------------------------------------------------------------------

 
end module ModuleMM5Format









