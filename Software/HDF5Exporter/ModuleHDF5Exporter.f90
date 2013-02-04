!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : HDF5Exporter files
! MODULE        : ExportHDF5ToTimeSerie
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : April 2004
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Module to create Time Series from HDF5 files
!                 HDF5 files for same time serie must have the same DT
!
!------------------------------------------------------------------------------

!DataFile
!
!   <BeginHDF5File>                                                                 
!   NAME                    : char                  [-]         !Name of HDF5 file 
!   <EndHDF5File>                                               !(specify one block for each HDF5 file to use)
!
!   EXPORT_TYPE             : int                   [1]         !Define the export type (if from cell or from area)
!                                                               !1 - Export CELL to Timeseries ; 2 - Export AREA to Timeseries                      
!   MASK_GRID               : char                  [-]         !Name of the MASK grid data file (if used)
!   AREA_FILL_VALUE         : int                   [-]         !FILL_VALUE that do not enter the calculations
!
!   START_TIME              : YYYY MM DD HH MM SS   [-]         !Start date of time series
!   END_TIME                : YYYY MM DD HH MM SS   [-]         !End date of time series
!
!   COMPUTE_RESIDUAL        : 0/1                   [1]         !Residual appearance in time series: 
!                                                               !0 - Appears not ; 1 - Appears   
!                                                               !(residual is meanless in time series produced: 
!                                                               !should always be 0)
!
!   VARIABLE_GRID           : 0/1                   [0]         !Use grid variable in time , usually
!                                                               !to satellite images 
!                                                               !0 - normal grid 1-variable grid
!   <BeginParameter>                                            
!   HDF_GROUP               : char                  [-]         !Path of the HDF5 group in HDF5 file for the property
!   PROPERTY                : char                  [-]         !Property name (should be equal to the one specified 
!   <EndParameter>                                              !in ModuleGlobalData)
!                                                               !(specify one block for each property)                                                                                 
!   <BeginTimeSerie>
!   NAME                    : char                  [-]         !Name for output time serie file
!   MASK_ID                 : int                   [-]         !ID on the MASK grid data file (if used MASK_GRID)
!   LOCALIZATION_I          : int                   [-]         !Horizontal - latitude location for time serie
!   LOCALIZATION_J          : int                   [-]         !Horizontal - longitude location for time serie
!   LOCALIZATION_K          : int                   [-]         !Vertical location for time serie
!   LATITUDE                : real                  [-]                           
!   LONGITUDE               : real                  [-]
!   <EndTimeSerie>                                              !(specify one block for each output time serie file)
!
!   WATERPOINTS_NAME        : char                  [-]         !Name of HDF5 item containing water points
!
!   GRID_FILENAME           : char                  [-]         !Name of the grid file name 
!   

Module ModuleExportHDF5ToTimeSerie

    use ModuleGlobalData         
    use ModuleTime               
    use ModuleEnterData,         only : ConstructEnterData, KillEnterData,             &
                                        GetData, ExtractBlockFromBuffer, Block_Unlock
    use ModuleTimeSerie,         only : StartTimeSerie, WriteTimeSerie, KillTimeSerie, &
                                        GetTimeSerieLocation, GetNumberOfTimeSeries,   &
                                        CorrectsCellsTimeSerie,                        &
                                        WriteSpecificTimeSerieLine
    use ModuleFunctions         
    use ModuleHDF5,              only : GetHDF5FileAccess, ConstructHDF5,              &
                                        GetHDF5GroupNumberOfItems, HDF5SetLimits,      &
                                        HDF5ReadData, GetHDF5GroupID, KillHDF5,        &
                                        GetHDF5DataSetExist     
    use ModuleDrawing,           only : IsPointInsidePolygon, SetLimits, T_Point,      &
                                        T_PointF, T_Polygon, GetSpecificPolygon,       &
                                        GetPolygonsNumber, New
    use ModuleHorizontalGrid,    only : ConstructHorizontalGrid, GetXYCellZ,           &
                                        GetHorizontalGridSize, KillHorizontalGrid,     &
                                        GetZCoordinates, UnGetHorizontalGrid
    use ModuleGridData,          only : ConstructGridData, GetGridData, UngetGridData, &
                                        GetIsGridData3D, KillGridData           

    ! To lower compilation costs are specified which subrotines are used in each module

    implicit none 

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartExportHDF5ToTimeSerie
    private ::          ReadKeyWords    
    private ::              ReadGlobalData
    private ::              ReadParameters
    private ::                  AddParameter
    private ::                  ConstructTimeSerieParameters
    private ::              ReadHDF5FileName            
    private ::                  AddHDF5File
    private ::                  ConstructHDF5File
    private ::          OpenAndDateHDF5Files
    private ::              HDF5TimeInstant
    private ::              HDF5Evaluator
    private ::                  SamePeriodHDF5 
    private ::              ObtainInstantsTimes
    private ::              AddTSHDF5File
    private ::                  CreateTSHDF5File
    private ::              ObtainWaterPoints
    private ::              AddTSInstantsTimes
    private ::          ConstructParameterList
    private ::          OpenOutputFiles    

    !Selector
    
    !Modifier
    public  :: ModifyExportHDF5ToTimeSerie
    private ::          OpenAndReadHDF5File
    private ::              ReadParameterFields
    private ::                  AddField
    private ::          ComputeIJonVariableGrids
    private ::              ComputeIJ

    !Destructor
    public  :: KillExportHDF5ToTimeSerie
    private ::      KillIndividualHDF5File
    private ::      KillIndividualParameter


    !Management
    
    !Interfaces----------------------------------------------------------------

    !Parameters----------------------------------------------------------------
    integer, parameter :: ExportCellToTimeseries = 1
    integer, parameter :: ExportAreaToTimeseries = 2
    
    !Types---------------------------------------------------------------------

    ! Definition of type T_Field
    type       T_Field
        character(len=StringLength)                 :: Name
        character(len=StringLength)                 :: Units
        integer                                     :: IDNumber
        type(T_Time)                                :: Date
        real, dimension(:,:  ),     pointer         :: Values2D
        real, dimension(:,:,: ),    pointer         :: Values3D 
        type(T_Field),              pointer         :: Next   => null()
        integer, dimension(:),      pointer         :: PositionI
        integer, dimension(:),      pointer         :: PositionJ
        integer, dimension(:),      pointer         :: PositionK
        real                                        :: Latitude
        real                                        :: Longitude
        logical                                     :: OutPoint = .false.
    end type  T_Field

    ! Definition of type T_Parameter
    type       T_Parameter
        character(len=StringLength)                 :: Name
        character(len=PathLength)                   :: Group
        type(T_Field),              pointer         :: FirstField
        type(T_Field),              pointer         :: CurrentField
        integer                                     :: Rank
        character(len=StringLength)                 :: Units
        integer                                     :: LastRank
        character(len=StringLength)                 :: LastName
        type(T_Parameter),          pointer         :: Next   => null()
    end type  T_Parameter

    type      T_TimeSeriesData
        integer                                     :: ID
        integer                                     :: Layer = -1
        real                                        :: Sum  !Sum
        real                                        :: QSum
        real                                        :: Mean !Mean
        real                                        :: Max  !Maximun value found
        real                                        :: Min  !Minimun value found
        real                                        :: SD   !Standard Deviation
        real                                        :: CV   !Coefficient of Variation
        logical                                     :: First
        integer                                     :: Count = 0
        real, dimension(:), pointer                 :: ParamsData
        type(T_TimeSeriesData), pointer             :: Next => null()
    end type  T_TimeSeriesData

    ! Definition of type T_TimeSerieTime
    type       T_TimeSerieTime
        type(T_Time)                                :: Time
        type(T_TimeSerieTime),      pointer         :: Next   => null()
    end type  T_TimeSerieTime

    ! Definition of type T_HDF5File
    type       T_HDF5File
        integer                                     :: HDFID = 0
        character(len=StringLength)                 :: Name
        type(T_Time)                                :: StartTime
        type(T_Time)                                :: EndTime
        type(T_Time)                                :: StartFieldTime
        type(T_Time)                                :: EndFieldTime
        integer                                     :: Rank
        integer                                     :: NumberOfInstants
        type(T_Time), dimension(:), pointer         :: InstantsArray 
        integer                                     :: StartInstant, EndInstant
        type(T_TimeSerieTime),      pointer         :: FirstInstantTime
        type(T_HDF5File),           pointer         :: Next   => null()
    end type  T_HDF5File

    ! Definition of type T_ExportHDF5ToTimeSerie
    private :: T_ExportHDF5ToTimeSerie
    type       T_ExportHDF5ToTimeSerie
        integer                                     :: InstanceID
        type(T_Time)                                :: StartTSTime, EndTSTime
        logical                                     :: UseStartTSTime, UseEndTSTime
        type(T_Time)                                :: ActualEndTSTime 
        logical                                     :: VariableGrid = .false.
        real                                        :: DT
        logical                                     :: VariableDT
        type(T_Size2D)                              :: Size2D
        type(T_Size3D)                              :: Size3D
        type (T_Parameter),            pointer      :: FirstParameter
        integer                                     :: ParameterNumber = 0
        character(len=StringLength), dimension(:), pointer :: ParameterList
        integer                                     :: ObjEnterData          = 0
        integer                                     :: ObjTime               = 0
        integer                                     :: ObjTimeSerie          = 0
        integer                                     :: ObjHorizontalGrid     = 0
        integer                                     :: ObjMaskHorizontalGrid = 0
        integer                                     :: ObjMaskGrid           = 0
        type(T_TimeSeriesData), pointer             :: FirstTimeSeriesData   => null()
        integer                                     :: TimeSeriesDataNumber  = 0
        type(T_Time)                                :: CurrentTime
        character(PathLength)                       :: DataFile
        character(PathLength)                       :: GridFileName
        logical                                     :: GridFileNameON
        logical                                     :: MaskIs3D
        real                                        :: AreaFillValue
        character(PathLength)                       :: PolygonsFile
        logical                                     :: PolygonON
        real, dimension(:,:), pointer               :: mask_2D
        type(T_HDF5File),              pointer      :: FirstHDF5File
        type(T_HDF5File),              pointer      :: FirstTSHDF5File
        type(T_HDF5File),              pointer      :: LastTSHDF5File
        type(T_TimeSerieTime),         pointer      :: FirstTimeSerieTime
        character(len=StringLength)                 :: WaterPointsName, WaterPointsGroup, TimeGroup
        integer                                     :: WaterPointsRank
        integer, dimension(:,:  ),     pointer      :: WaterPoints2D
        integer, dimension(:,:,: ),    pointer      :: WaterPoints3D
        integer                                     :: DecimationFactor
        integer                                     :: ExportType            = ExportCellToTimeseries
        logical                                     :: CheckPropertyName = .true.
        logical                                     :: UsePointsMatrix   = .true.
        type(T_ExportHDF5ToTimeSerie), pointer      :: Next   => null()
    end type  T_ExportHDF5ToTimeSerie

    !Global Module Variables
    type (T_ExportHDF5ToTimeSerie),    pointer      :: Me

    !--------------------------------------------------------------------------
    
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartExportHDF5ToTimeSerie(ObjExportHDF5ToTimeSerieID, DataFile)

        !Arguments---------------------------------------------------------------
        integer                                     :: ObjExportHDF5ToTimeSerieID 
        character(PathLength), intent(IN)           :: DataFile

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        !Assures nullification of the global variable
        allocate(Me)

        !Returns ID
        ObjExportHDF5ToTimeSerieID          = 1

        !Atribute the name of data file            
        Me%DataFile = DataFile

        ! Read keyword file and HDF5 file
        call ReadKeywords

        if (Me%ParameterNumber > 0) then

            ! Open and obtain key features of the HDF5 files
            call OpenAndDateHDF5Files 

            ! Construct parameter list
            call ConstructParameterList

            ! Open and write header of Output Files
            call OpenOutputFiles

        end if

        !----------------------------------------------------------------------

    end subroutine StartExportHDF5ToTimeSerie
 
    !--------------------------------------------------------------------------

    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call ConstructEnterData (Me%ObjEnterData, Me%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ReadKeywords - ModuleExportHDF5ToTimeSerie - ERR10'

        call ReadGlobalData

        call ReadParameters

        call ReadHDF5FileName

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadKeywords - ModuleExportHDF5ToTimeSerie - ERR20'
        end if

    end subroutine ReadKeywords

    !--------------------------------------------------------------------------

    subroutine AddParameter (ObjParameter)

        !Arguments-------------------------------------------------------------
        type (T_Parameter),     pointer             :: ObjParameter

        !Local-----------------------------------------------------------------
        type (T_Parameter),     pointer             :: PreviousParameter
        type (T_Parameter),     pointer             :: NewParameter

        !Begin-----------------------------------------------------------------

        !Allocates new Parameter
        allocate (NewParameter)
        nullify  (NewParameter%Next)

        !Insert new Parameter into list and makes current ?? point to it
        if (.not. associated(Me%FirstParameter)) then
            Me%FirstParameter         => NewParameter
            ObjParameter              => NewParameter
            Me%ParameterNumber = 1
        else
            PreviousParameter         => Me%FirstParameter
            ObjParameter              => Me%FirstParameter%Next
            do while (associated(ObjParameter))
                PreviousParameter     => ObjParameter
                ObjParameter          => ObjParameter%Next
            enddo
            ObjParameter              => NewParameter
            PreviousParameter%Next    => NewParameter
            ! Count number of parameters in list
            Me%ParameterNumber = Me%ParameterNumber + 1
        end if

    end subroutine AddParameter

    !--------------------------------------------------------------------------

    subroutine ConstructTimeSerieParameters (NewParameter)

        !Arguments-------------------------------------------------------------
        type (T_Parameter),      pointer            :: NewParameter

        !External--------------------------------------------------------------
        integer                                     :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        
        ! Obtain parameter name
        call GetData(NewParameter%Name,                         &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'PROPERTY',                 &
                     ClientModule = 'ExportToTimeSerie',        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_ .OR. iflag .EQ. 0)          &
        stop 'ConstructTimeSerieParameters - ModuleExportHDF5ToTimeSerie - ERR10'
        
        if (Me%CheckPropertyName) then
            if (.not.CheckPropertyName(NewParameter%Name)) then
                write(*,*)
                write(*,*) 'The property name is not recognised by the model.'
                write(*,*) 'ConstructTimeSerieParameters - ModuleExportHDF5ToTimeSerie - WARN10' 
            end if
        endif           
        
        ! Obtain parameter group
        call GetData(NewParameter%Group,                        &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'HDF_GROUP',                &
                     ClientModule = 'ExportToTimeSerie',        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_ .OR. iflag .EQ. 0)           &
        stop 'ConstructTimeSerieParameters - ModuleExportHDF5ToTimeSerie - ERR20'

    end subroutine ConstructTimeSerieParameters 

    !--------------------------------------------------------------------------

    type(T_Time) function HDF5TimeInstant(Instant, ObjHDF5File)

        !Arguments-------------------------------------------------------------
        integer                                     :: Instant
        type(T_HDF5File), pointer                   :: ObjHDF5File
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer                 :: TimeVector

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (ObjHDF5File%HDFID, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = ObjHDF5File%HDFID,        &
                             GroupName      = "/"//trim(Me%TimeGroup),  &
                             Name           = trim(Me%TimeGroup),       &
                             Array1D        = TimeVector,               &
                             OutputNumber   = Instant,                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                      &
        stop 'HDF5TimeInstant - ModuleExportHDF5ToTimeSerie - ERR10'

        call SetDate(HDF5TimeInstant, Year  = TimeVector(1),            &
                     Month  = TimeVector(2), Day      = TimeVector(3),  &
                     Hour   = TimeVector(4), Minute   = TimeVector(5),  &
                     Second = TimeVector(6))


        deallocate(TimeVector)
        nullify   (TimeVector)

    end function HDF5TimeInstant

    !--------------------------------------------------------------------------

    subroutine AddField (FirstField, ObjField)

        !Arguments-------------------------------------------------------------
        type (T_Field), pointer                     :: FirstField
        type (T_Field), pointer                     :: ObjField
        
        !Local-----------------------------------------------------------------
        type (T_Field), pointer                     :: NewField
        type (T_Field), pointer                     :: PreviousField
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewField)
        nullify  (NewField%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstField)) then
        !FirstField should be the same for all HDF5 files 
            FirstField              => NewField
            ObjField                => NewField
        else
            PreviousField           => FirstField
            ObjField                => FirstField%Next
            do while (associated(ObjField))
                PreviousField       => ObjField
                ObjField            => ObjField%Next
            enddo
            ObjField                => NewField
            PreviousField%Next      => NewField
        endif

    end subroutine AddField

    !--------------------------------------------------------------------------

    subroutine ConstructParameterList

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: nParameters
        type (T_Parameter),       pointer           :: ParameterX
        
        !Begin-----------------------------------------------------------------

        !Allocates Parameter list and Parameter group list
        if (Me%ExportType == ExportCellToTimeSeries) then
            allocate(Me%ParameterList(Me%ParameterNumber))
        else
            allocate(Me%ParameterList(Me%ParameterNumber * 7))
        endif

        !Fills up ParameterList
        ParameterX   => Me%FirstParameter
        nParameters  = 0

        do while (associated(ParameterX))
        
            if (Me%ExportType == ExportCellToTimeSeries) then
                nParameters = nParameters + 1
                Me%ParameterList(nParameters) = trim(adjustl(ParameterX%Name))
            else
                nParameters = nParameters + 1
                Me%ParameterList(nParameters) = trim(adjustl(ParameterX%Name))//'_sum'
                nParameters = nParameters + 1
                Me%ParameterList(nParameters) = trim(adjustl(ParameterX%Name))//'_mean'
                nParameters = nParameters + 1
                Me%ParameterList(nParameters) = trim(adjustl(ParameterX%Name))//'_min'
                nParameters = nParameters + 1
                Me%ParameterList(nParameters) = trim(adjustl(ParameterX%Name))//'_max'
                nParameters = nParameters + 1
                Me%ParameterList(nParameters) = trim(adjustl(ParameterX%Name))//'_sd'
                nParameters = nParameters + 1
                Me%ParameterList(nParameters) = trim(adjustl(ParameterX%Name))//'_cv'
                nParameters = nParameters + 1
                Me%ParameterList(nParameters) = trim(adjustl(ParameterX%Name))//'_count'
            endif

            ParameterX=>ParameterX%Next

        enddo

    end subroutine ConstructParameterList

    !--------------------------------------------------------------------------

    subroutine ReadParameterFields(ObjHDF5File, ObjParameter)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: ObjHDF5File
        type(T_Parameter), pointer                  :: ObjParameter
          
        !Local-----------------------------------------------------------------
        integer                                     :: NewCurrentInstant
        integer                                     :: CurrentInstant
        type (T_Field), pointer                     :: NewField
        integer                                     :: Count = 1
        integer, dimension(7)                       :: Dimensions
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !read/copy parameter fields
        write(*,*)'Reading '//trim(ObjParameter%Name)//' fields'

        NewCurrentInstant = 0

        do CurrentInstant = ObjHDF5File%StartInstant, ObjHDF5File%EndInstant,           &
                            (Me%DecimationFactor + 1) 

            NewCurrentInstant = NewCurrentInstant + 1
                
            !construct new fields
            call AddField(ObjParameter%FirstField, NewField)
            NewField%IDNumber = Count
            Count = Count + 1

            !get field ID, Rank and Dimensions
            !(this needs to be done for each instant)
            call GetHDF5GroupID(ObjHDF5File%HDFID, ObjParameter%Group,                  &
                                CurrentInstant, NewField%Name,                          &
                                NewField%Units, ObjParameter%Rank,                      &
                                Dimensions,                                             &
                                STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_)                                                &  
            stop 'ReadParameterFields - ModuleExportHDF5ToTimeSerie - ERR10'
               
            NewField%Units = trim(NewField%Units)
            NewField%Date  = ObjHDF5File%InstantsArray(NewCurrentInstant)
            
            !get field values
            select case (ObjParameter%Rank)

                case(2)
                ! The HDF5 file contains 2D data

                    !calculate Size2D
                    Me%Size2D%ILB = 1
                    Me%Size2D%IUB = Dimensions(1)
                    Me%Size2D%JLB = 1
                    Me%Size2D%JUB = Dimensions(2)

                    !allocate field
                    nullify (NewField%Values2D)
                    allocate(NewField%Values2D(Me%Size2D%ILB:Me%Size2D%IUB,             &
                                               Me%Size2D%JLB:Me%Size2D%JUB))
                       
                    call HDF5SetLimits (ObjHDF5File%HDFID, Me%Size2D%ILB,               &
                                        Me%Size2D%IUB, Me%Size2D%JLB,Me%Size2D%JUB,     &
                                        STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        & 
                    stop 'ReadParameterFields - ModuleExportHDF5ToTimeSerie - ERR20'
        
                    !read field
                    call HDF5ReadData(ObjHDF5File%HDFID, ObjParameter%Group,            &
                                      trim(NewField%Name),                              &
                                      Array2D      = NewField%Values2D,                 &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ReadParameterFields - ModuleExportHDF5ToTimeSerie - ERR30'

                case(3)
                ! The HDF5 file contains 3D data

                    !calculate Size3D
                    Me%Size3D%ILB = 1
                    Me%Size3D%IUB = Dimensions(1)
                    Me%Size3D%JLB = 1
                    Me%Size3D%JUB = Dimensions(2)
                    Me%Size3D%KLB = 1
                    Me%Size3D%KUB = Dimensions(3)
                        
                    !allocate field
                    nullify (NewField%Values3D)
                    allocate(NewField%Values3D(Me%Size3D%ILB:Me%Size3D%IUB,             &
                                               Me%Size3D%JLB:Me%Size3D%JUB,             &
                                               Me%Size3D%KLB:Me%Size3D%KUB))
                        
                    call HDF5SetLimits  (ObjHDF5File%HDFID, Me%Size3D%ILB,              &
                                         Me%Size3D%IUB, Me%Size3D%JLB,Me%Size3D%JUB,    &
                                         Me%Size3D%KLB,Me%Size3D%KUB, STAT = STAT_CALL)                
                    if (STAT_CALL .NE. SUCCESS_)                                        & 
                    stop 'ReadParameterFields - ModuleExportHDF5ToTimeSerie - ERR40'
        
                    !read field
                    call HDF5ReadData(ObjHDF5File%HDFID, ObjParameter%Group,            &
                                      trim(NewField%Name),                              &
                                      Array3D      = NewField%Values3D,                 &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ReadParameterFields - ModuleExportHDF5ToTimeSerie - ERR50'
 
            case default 
                    
                write(*,*)'Time Serie created only for 2D or 3D HDF5 files.'
                stop 'ReadParameterFields - ModuleExportHDF5ToTimeSerie - ERR60'
                    
            end select

        end do

    end subroutine ReadParameterFields

    !--------------------------------------------------------------------------

    subroutine OpenOutputFiles

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, TimeSerieNumber, dn, Id, Jd
        real                                        :: CoordX, CoordY
        logical                                     :: CoordON

        !Begin-----------------------------------------------------------------

        !Call StartComputeTime for the whole Time Serie

        call StartComputeTime(Me%ObjTime, Me%FirstTSHDF5File%InstantsArray(1), Me%FirstTSHDF5File%InstantsArray(1),      &
                              Me%LastTSHDF5File%InstantsArray                       &
                              (Me%LastTSHDF5File%NumberOfInstants),                 &
                              DT = Me%DT, VariableDT = Me%VariableDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
        stop 'OpenOutputFiles - ModuleExportHDF5ToTimeSerie - ERR010'

        if (Me%ExportType == ExportAreaToTimeSeries) then
        
            call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime, &
                                trim(Me%DataFile),           &
                                Me%ParameterList, "ets",     &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'OpenOutputFiles - ModuleExportHDF5ToTimeSerie - ERR020'
        
        else

            select case(Me%WaterPointsRank)

                case(2)
                       
                    !Constructs Time Serie Header
                    call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                    &
                                        trim(Me%DataFile),                              &
                                        Me%ParameterList, "ets",                        &
                                        WaterPoints2D = Me%WaterPoints2D,               &
                                        STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                    stop 'OpenOutputFiles - ModuleExportHDF5ToTimeSerie - ERR030'

                    deallocate(Me%WaterPoints2D)
                    nullify(Me%WaterPoints2D)

                case(3)

                    !Constructs Time Serie Header
                    call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                    &
                                        trim(Me%DataFile),                              &
                                        Me%ParameterList, "ets",                        &
                                        WaterPoints3D = Me%WaterPoints3D,               &
                                        STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                    stop 'OpenOutputFiles - ModuleExportHDF5ToTimeSerie - ERR040'

                    deallocate(Me%WaterPoints3D)
                    nullify(Me%WaterPoints3D)

            end select

        endif

        if(.not. Me%VariableGrid)then

            !Corrects if necessary the cell of the time serie based in the time serie coordinates
            call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OpenOutputFiles - ModuleExportHDF5ToTimeSerie - ERR050'

            do dn = 1, TimeSerieNumber

                call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                          CoordX   = CoordX,                                &
                                          CoordY   = CoordY,                                & 
                                          CoordON  = CoordON,                               &
                                          STAT     = STAT_CALL)
                if (CoordON) then

                    if (.not. Me%GridFileNameON) stop 'OpenOutputFiles - ModuleExportHDF5ToTimeSerie - ERR060'

                    call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OpenOutputFiles - ModuleExportHDF5ToTimeSerie - ERR070'
                    
                    if (.not.Me%PolygonON) then

                        if (Id < 0 .or. Jd < 0) stop 'OpenOutputFiles - ModuleExportHDF5ToTimeSerie - ERR080'

                        call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OpenOutputFiles - ModuleExportHDF5ToTimeSerie - ERR090'
                        
                        write(*,*) 'I = ', Id
                        write(*,*) 'J = ', Jd
                        
                    endif
                endif
            enddo


            if (Me%GridFileNameON) then
                call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OpenOutputFiles - ModuleExportHDF5ToTimeSerie - ERR100'
            endif

        end if


        !Deallocates Parameter list
        deallocate(Me%ParameterList)
        nullify(Me%ParameterList)

    end subroutine OpenOutputFiles
  
    !--------------------------------------------------------------------------

    subroutine ReadGlobalData

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        type(T_Size2D)                              :: WorkSize2D, Size2D
        integer                                     :: STAT_CALL, iflag, i, j, n, PolygonsNumber
        logical                                     :: exist
        character(len=PathLength)                   :: mask_grid_filename
        type (T_PointF),   pointer                  :: Point
        type (T_Polygon),  pointer                  :: MaskPolygons, CurrentPolygon
        real,   dimension(:,:), pointer             :: CoordX, CoordY

        !Begin-----------------------------------------------------------------

        !Check if will do the export of an area to a timeseries
        call GetData(Me%ExportType, Me%ObjEnterData, iflag, &
                     keyword      = 'EXPORT_TYPE',          &
                     SearchType   = FromFile,               &
                     Default      = 1,                      &
                     ClientModule = 'ExportToTimeSerie',    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR010'

        Me%UsePointsMatrix = .true.
        
        if (Me%ExportType == ExportAreaToTimeseries) then
            
            call GetData(Me%UsePointsMatrix, Me%ObjEnterData, iflag,  &
                         keyword      = 'USE_POINTS',                 &
                         SearchType   = FromFile,                     &
                         ClientModule = 'ExportToTimeSerie',          &
                         default      = .true.,                       &
                         STAT         = STAT_CALL)
            if ((STAT_CALL /= SUCCESS_) .or. (iflag < 1)) &
                stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR019'                          
            
            
            call GetData(Me%PolygonsFile, Me%ObjEnterData, iflag,                       &
                         keyword      = 'POLYGONS_FILE',                                &
                         SearchType   = FromFile,                                       &
                         Default      = "******.***",                                   &
                         ClientModule = 'ExportToTimeSerie',                            &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR040'               
                                                            
            if (iflag == 1) then
                Me%PolygonON = .true.
            else
                Me%PolygonON = .false.
            endif
            
            if (Me%PolygonON) then
            
                if (.not. Me%GridFileNameON) then
                
                    write(*,*)'Grid file not define'
                    stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR045'
                                
                endif
            
                call GetHorizontalGridSize(Me%ObjHorizontalGrid, Size = Size2D,         &
                                           WorkSize = WorkSize2D, STAT = STAT_CALL)           
                if (STAT_CALL /= SUCCESS_) &
                   stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR050'            
            
                Me%MaskIs3D = .false. 
                
                allocate(Me%mask_2D(Size2D%ILB:Size2D%IUB, Size2D%JLB:Size2D%JUB))
                
                Me%mask_2D(:,:) = -99
                
                call New(MaskPolygons, Me%PolygonsFile)
                
                PolygonsNumber = GetPolygonsNumber(MaskPolygons)
                
                call GetZCoordinates(Me%ObjHorizontalGrid, CoordX, CoordY, STAT= STAT_CALL)           
                if (STAT_CALL /= SUCCESS_) &
                   stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR060'   
                   
                allocate(Point)                   
                
                do n=1, PolygonsNumber
                    do j=WorkSize2D%JLB, WorkSize2D%JUB
                    do i=WorkSize2D%ILB, WorkSize2D%IUB
                        Point%X = CoordX(i, j)
                        Point%Y = CoordY(i, j)
                        call GetSpecificPolygon(MaskPolygons, n, CurrentPolygon)
                        
                        if (IsPointInsidePolygon(Point, CurrentPolygon)) then
                            Me%mask_2D(i,j) = n
                        endif
                    enddo
                    enddo
                
                enddo
                
                deallocate(Point)

                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordX, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR070'                 

                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordY, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR080'                 

            
            else                                   
            
                call GetData(mask_grid_filename, Me%ObjEnterData, iflag,  &
                             keyword      = 'MASK_GRID',                  &
                             SearchType   = FromFile,                     &
                             ClientModule = 'ExportToTimeSerie',          &
                             STAT         = STAT_CALL)
                if ((STAT_CALL /= SUCCESS_) .or. (iflag < 1)) &
                    stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR020'
                    
                call ConstructHorizontalGrid(Me%ObjMaskHorizontalGrid, &
                                             mask_grid_filename,       &
                                             STAT = STAT_CALL)           
                if (STAT_CALL /= SUCCESS_) &
                   stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR030'
            
                call ConstructGridData (Me%ObjMaskGrid,                 &
                                        Me%ObjMaskHorizontalGrid,       &
                                        FileName = mask_grid_filename,  &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR040'            
                
                call GetIsGridData3D(Me%ObjMaskGrid, Me%MaskIs3D, STAT = STAT_CALL)            
                if (STAT_CALL /= SUCCESS_) &
                    stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR050'                                
                
            endif
                        
            call LoadTimeSeriesMaskID
                            
            call GetData(Me%AreaFillValue, Me%ObjEnterData, iflag,  &
                         keyword      = 'AREA_FILL_VALUE',            &
                         SearchType   = FromFile,                     &
                         ClientModule = 'ExportToTimeSerie',          &
                         Default      = -99.0,                        & 
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR051'                                     

        endif
            
        call GetData(Me%CheckPropertyName,                      &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromFile,                   &
                     keyword      = 'CHECK_PROPERTY',           &
                     ClientModule = 'ExportToTimeSerie',        &
                     Default      = .true.,                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)          &
        stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR52'               
        
        
        ! Obtain the start and end times for the Time Serie
        ! Start Time
        call GetData(Me%StartTSTime, Me%ObjEnterData, iflag,                &
                     keyword      = 'START_TIME',                           &
                     SearchType   = FromFile,                               &
                     ClientModule = 'ExportToTimeSerie',                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR060'   
        if (iflag == 0) then
            Me%UseStartTSTime = .false.
        else
            Me%UseStartTSTime = .true.
        endif

        ! End Time 
        call GetData(Me%EndTSTime,   Me%ObjEnterData, iflag,                &
                     keyword      = 'END_TIME',                             &
                     SearchType   = FromFile,                               &
                     ClientModule = 'ExportToTimeSerie',                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR070'   
        if (iflag == 0) then
            Me%UseEndTSTime = .false.
        else
            Me%UseEndTSTime = .true.
        endif

        ! Verifies Time Variables
        if (Me%UseStartTSTime .and. Me%UseEndTSTime) then
            if (Me%EndTSTime .lt. Me%StartTSTime) then
                write (*,*) 'Time Serie End Time is BEFORE Time Serie Start Time'
                write (*,*) 'Module :','ExportHDF5ToTimeSerie'
                stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR080'
            endif

            if (Me%EndTSTime .eq. Me%StartTSTime) then
                write (*,*) 'Time Serie End Time is EQUAL Time Serie Start Time'
                write (*,*) 'Module :','ExportHDF5ToTimeSerie'
                stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR090'
            endif
        endif
            
        call GetData(Me%VariableGrid,Me%ObjEnterData, iflag,                &
                     keyword      = 'VARIABLE_GRID',                        &
                     SearchType   = FromFile,                               &
                     ClientModule = 'ExportToTimeSerie',                    &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR100'

        if (Me%UsePointsMatrix) then
        
            call GetData(Me%WaterPointsName,   Me%ObjEnterData, iflag,          &
                         keyword      = 'WATERPOINTS_NAME',                     &
                         SearchType   = FromFile,                               &
                         ClientModule = 'ExportToTimeSerie',                    &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                          &
                stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR110'
                
            call GetData(Me%WaterPointsGroup,   Me%ObjEnterData, iflag,         &
                         keyword      = 'WATERPOINTS_GROUP',                    &
                         SearchType   = FromFile,                               &
                         ClientModule = 'ExportToTimeSerie',                    &
                         default      = "/Grid",                                &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR120'
        
        endif
        
        call GetData(Me%GridFileName,                                       &
                     Me%ObjEnterData, iflag,                                &
                     SearchType   = FromFile,                               &
                     keyword      = 'GRID_FILENAME',                        &
                     ClientModule = 'ConvertToHDF5',                        &
                     STAT         = STAT_CALL)                                          
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR130'
                       
        if (iflag==0) then
            Me%GridFileNameON = .false.
        else
            Me%GridFileNameON = .true.

            inquire(FILE = Me%GridFileName, EXIST = exist)
            if (.not. exist) then
                write(*,*)'Grid file does not exist'
                stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR140'
            endif

            call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR150'
        
        endif

        call GetData(Me%TimeGroup,   Me%ObjEnterData, iflag,                &
                     keyword      = 'TIME_GROUP',                           &
                     SearchType   = FromFile,                               &
                     ClientModule = 'ExportToTimeSerie',                    &
                     Default      = "Time",                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ )                                         &
            stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR160'   


        ! Get decimation factor
        call GetData(Me%DecimationFactor,                                   &
                     Me%ObjEnterData, iflag,                                &
                     SearchType = FromFile,                                 &
                     keyword    = 'DECIMATION_FACTOR',                      &
                     Default    = 0,                                        &
                     ClientModule = 'ExportToTimeSerie',                    &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ReadGlobalData - ModuleExportHDF5ToTimeSerie - ERR170'

    end subroutine ReadGlobalData

    !--------------------------------------------------------------------------
        
    subroutine LoadTimeSeriesMaskID
    
        type (T_TimeSeriesData), pointer            :: NewProperty  => null()
        type (T_TimeSeriesData), pointer            :: LastProperty => null()
        integer                                     :: ClientNumber
        integer                                     :: STAT_CALL
        logical                                     :: BlockFound
        integer                                     :: iflag
            
do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,    &
                                        ClientNumber,       &
                                        '<BeginTimeSerie>', &
                                        '<EndTimeSerie>',   &
                                        BlockFound,         &
                                        STAT = STAT_CALL)

cd1 :       if (STAT_CALL .EQ. SUCCESS_) then   
 
cd2 :           if (BlockFound) then                                                  

                    ! Construct a New Property 
                    allocate (NewProperty)
                    NewProperty%Next => null()
                    
                    call GetData(NewProperty%ID, Me%ObjEnterData, iflag, &
                                 keyword      = 'MASK_ID',               &
                                 SearchType   = FromBlock,               &
                                 ClientModule = 'ExportToTimeSerie',     &
                                 STAT         = STAT_CALL)
                    if ((STAT_CALL /= SUCCESS_) .or. (iflag /= 1)) &
                        stop 'LoadTimeSeriesMaskID - ModuleExportHDF5ToTimeSerie - ERR010'

                    call GetData(NewProperty%Layer, Me%ObjEnterData, iflag, &
                                 keyword      = 'LAYER',                    &
                                 SearchType   = FromBlock,                  &
                                 Default      = -1,                         &
                                 ClientModule = 'ExportToTimeSerie',        &
                                 STAT         = STAT_CALL)
                    if ((STAT_CALL /= SUCCESS_)) &
                        stop 'LoadTimeSeriesMaskID - ModuleExportHDF5ToTimeSerie - ERR020'
                        
                    ! Add new Timeseries to the List 
                    if (.not. associated(Me%FirstTimeSeriesData)) then
                        Me%FirstTimeSeriesData => NewProperty
                        LastProperty           => NewProperty
                    else                        
                        LastProperty%Next => NewProperty
                        LastProperty      => NewProperty
                    end if
                    Me%TimeSeriesDataNumber = Me%TimeSeriesDataNumber + 1
                    
                else
                    
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'LoadTimeSeriesMaskID - ModuleExportHDF5ToTimeSerie - ERR030'

                    exit do1    !No more blocks

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'LoadTimeSeriesMaskID - ModuleExportHDF5ToTimeSerie - ERR040'
            end if cd1
        end do do1    
    
    end subroutine LoadTimeSeriesMaskID
    
    !--------------------------------------------------------------------------

    subroutine ReadParameters

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, ClientNumber
        type (T_Parameter),       pointer           :: NewParameter
        logical                                     :: BlockFound
        logical                                     :: AtLeastOneBlock = .false.

        !Begin-----------------------------------------------------------------

        ! Obtain Parameters for the Time Serie
do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<BeginParameter>',   &
                                        block_end       = '<EndParameter>',     &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  

                    AtLeastOneBlock = .true.
                    
                    call AddParameter (NewParameter)

                    call ConstructTimeSerieParameters (NewParameter)

                    nullify(NewParameter)

                else cd2
                    call Block_Unlock(Me%ObjEnterData,                          & 
                                      ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ReadParameters - ModuleExportHDF5ToTimeSerie - ERR10'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ReadParameters - ModuleExportHDF5ToTimeSerie - ERR20'
            else cd1
                stop 'ReadParameters - ModuleExportHDF5ToTimeSerie - ERR30'
            end if cd1

        end do do1

        if (.not. AtLeastOneBlock) then                                            
            write(*,*) 'No property block is indicated in input file. '
            stop 'ReadParameters - ModuleExportHDF5ToTimeSerie - ERR40'
        end if

    end subroutine ReadParameters

    !--------------------------------------------------------------------------

    subroutine ReadHDF5FileName

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, ClientNumber
        type (T_HDF5File),       pointer            :: NewHDF5File
        logical                                     :: BlockFound
        logical                                     :: AtLeastOneBlock = .false.

        !Begin-----------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<BeginHDF5File>',    &
                                        block_end       = '<EndHDF5File>',      &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                    
                    AtLeastOneBlock = .true.

                    call AddHDF5File                     (NewHDF5File)

                    call ConstructHDF5File               (NewHDF5File)

                    nullify(NewHDF5File)

                else cd2
                    call Block_Unlock(Me%ObjEnterData,                          & 
                                      ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ReadHDF5FileName - ModuleExportHDF5ToTimeSerie - ERR10'

                    exit do1

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ReadHDF5FileName - ModuleExportHDF5ToTimeSerie - ERR20'
            else cd1
                stop 'ReadHDF5FileName - ModuleExportHDF5ToTimeSerie - ERR30'
            end if cd1

        end do do1

        if (.not. AtLeastOneBlock) then                                            
            write(*,*) 'No HDF5 file block is indicated in input file. '
            stop 'ReadHDF5FileName - ModuleExportHDF5ToTimeSerie - ERR40'
        end if

    end subroutine ReadHDF5FileName

    !--------------------------------------------------------------------------

    subroutine AddHDF5File(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type (T_HDF5File),     pointer              :: ObjHDF5File

        !Local-----------------------------------------------------------------
        type (T_HDF5File),     pointer              :: PreviousHDF5File
        type (T_HDF5File),     pointer              :: NewHDF5File

        !Begin-----------------------------------------------------------------

        !Allocates new HDF5File
        allocate (NewHDF5File)
        nullify  (NewHDF5File%Next)

        !Insert new HDF5 file into list and makes current point to it
        if (.not. associated(Me%FirstHDF5File)) then
            Me%FirstHDF5File         => NewHDF5File
            ObjHDF5File              => NewHDF5File
        else
            PreviousHDF5File         => Me%FirstHDF5File
            ObjHDF5File              => Me%FirstHDF5File%Next
            do while (associated(ObjHDF5File))
                PreviousHDF5File     => ObjHDF5File
                ObjHDF5File          => ObjHDF5File%Next
            enddo
            ObjHDF5File              => NewHDF5File
            PreviousHDF5File%Next    => NewHDF5File
        end if

    end subroutine AddHDF5File

    !--------------------------------------------------------------------------

    subroutine ConstructHDF5File(NewHDF5File)

        !Arguments-------------------------------------------------------------
        type (T_HDF5File),      pointer             :: NewHDF5File

        !External--------------------------------------------------------------
        integer                                     :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        
        ! Obtain HDF5 file name
        call GetData(NewHDF5File%Name,                                  &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'NAME',                             &
                     ClientModule = 'ExportToTimeSerie',                &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                    &
        stop 'ConstructHDF5File - ModuleExportHDF5ToTimeSerie - ERR10'

    end subroutine ConstructHDF5File

    !--------------------------------------------------------------------------

    subroutine OpenAndDateHDF5Files

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type (T_HDF5File), pointer                  :: HDF5FileX
        character(len=StringLength)                 :: ParameterName
        logical                                     :: exist, FirstTime
        integer                                     :: HDF5_READ
        type (T_Parameter), pointer                 :: ObjParameter
        type(T_Time), dimension(:), pointer         :: AuxInstantsArray 
        integer                                     :: CurrentInstant, AuxNumberInstants
        logical                                     :: Relevant
        real                                        :: LastDT, AuxDT
        real                                        :: Year, Month, Day, Hour 
        real                                        :: Minute, Second
        integer                                     :: Count
      
        !Begin-----------------------------------------------------------------

        FirstTime = .true.

        HDF5FileX => Me%FirstHDF5File
        
        !In a DO cycle open all HDF5 files provided by the user
        do while (associated(HDF5FileX))

            !Verifies if file exists
            inquire(FILE = HDF5FileX%Name, EXIST = exist)
            if (.not. exist) then
                write(*,*)'HDF5 file does not exist:'//trim(HDF5FileX%Name)
                stop 'OpenAndDateHDF5Files - ModuleExportHDF5ToTimeSerie - ERR10'
            endif

            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

            !Open HDF5 file
            call ConstructHDF5 (HDF5FileX%HDFID, trim(HDF5FileX%Name),              &
                                HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'OpenAndDateHDF5Files - ModuleExportHDF5ToTimeSerie - ERR20'

            !Obtain start and end times of HDF5 file
            !(obtain number of instants) 
            call GetHDF5GroupNumberOfItems(HDF5FileX%HDFID, "/"//(Me%TimeGroup),    &
                                           HDF5FileX%NumberOfInstants, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            & 
            stop 'OpenAndDateHDF5Files - ModuleExportHDF5ToTimeSerie - ERR30'

            !(obtain HDF5 start time)
            HDF5FileX%StartTime = HDF5TimeInstant(1, HDF5FileX)

            !(obtain HDF5 end time)
            HDF5FileX%EndTime = HDF5TimeInstant(HDF5FileX%NumberOfInstants, HDF5FileX)
            
            if (FirstTime .and. (.not. Me%UseStartTSTime)) &
                Me%StartTSTime = HDF5FileX%StartTime

            !Get info about the rank and variables present
            !(only data needed for checking are obtained from each file)
            ObjParameter => Me%FirstParameter

            do while(associated(ObjParameter))

                !get field ID, Rank
                call GetHDF5GroupID(HDF5FileX%HDFID, ObjParameter%Group,            &
                                1, ParameterName, ObjParameter%Units,               &
                                ObjParameter%Rank,                                  &
                                STAT = STAT_CALL)                                

                !check if file contains parameter required
                if (STAT_CALL .NE. SUCCESS_) then  
                    write(*,*)'HDF5 file do not contain parameter required:'//trim(HDF5FileX%Name)
                    write(*,*)'Parameter required:'//trim(ObjParameter%Name)
                    stop 'OpenAndDateHDF5Files - ModuleExportHDF5ToTimeSerie - ERR40'
                end if

                if (.not. FirstTime) then
                    !check if name of parameter is consistent with the previous
                    if (ParameterName .NE. ObjParameter%LastName) then
                        write(*,*)'HDF5 file do not contain parameter required:'//trim(HDF5FileX%Name)
                        write(*,*)'Parameter required:'//trim(ObjParameter%Name)
                        stop 'OpenAndDateHDF5Files - ModuleExportHDF5ToTimeSerie - ERR50'
                    end if

                    !check if rank of parameter is consistent with the previous
                    if (ObjParameter%Rank .NE. ObjParameter%LastRank) then
                        write(*,*)'File rank not consistent with previous rank:'//trim(HDF5FileX%Name)
                        write(*,*)'Parameter:'//trim(ObjParameter%Name)
                        stop 'OpenAndDateHDF5Files - ModuleExportHDF5ToTimeSerie - ERR60'
                    end if

                end if 

                ObjParameter%LastRank = ObjParameter%Rank
                ObjParameter%LastName = ParameterName               
                
                ObjParameter => ObjParameter%Next

            end do

            !Check if the HDF5 file is relevant for the time series
            call HDF5Evaluator(HDF5FileX, Relevant)

            !If HDF5 file is relevant then obtain key parameters and instants
            if (Relevant) then

                !Get useful time information from file:
                !Set instant array for this file 
                allocate(AuxInstantsArray(1:HDF5FileX%NumberOfInstants))

                !Fill array with instants
                do CurrentInstant = 1, HDF5FileX%NumberOfInstants

                   AuxInstantsArray(CurrentInstant) =                               &
                        HDF5TimeInstant(CurrentInstant, HDF5FileX)

                end do

                !Get start and end instants for this file
                !select time window begin
                do CurrentInstant = 1, HDF5FileX%NumberOfInstants

                    if (AuxInstantsArray(CurrentInstant)                            &
                        .ge. HDF5FileX%StartFieldTime) then

                        HDF5FileX%StartInstant = CurrentInstant
                        HDF5FileX%StartFieldTime = HDF5TimeInstant(CurrentInstant,  &
                                                                   HDF5FileX)
                
                        exit

                    end if

                end do
    
                !select time window end
                do CurrentInstant = HDF5FileX%StartInstant,                         &
                        HDF5FileX%NumberOfInstants, (Me%DecimationFactor + 1)

                    if (AuxInstantsArray(CurrentInstant)                            &
                        .eq. HDF5FileX%EndFieldTime) then

                        HDF5FileX%EndInstant = (CurrentInstant)
                        HDF5FileX%EndFieldTime = HDF5TimeInstant(CurrentInstant, HDF5FileX)

                        exit

                    end if

                    if (AuxInstantsArray(CurrentInstant)                            &
                        .ge. HDF5FileX%EndFieldTime) then

                        HDF5FileX%EndInstant = (CurrentInstant-1)
                        HDF5FileX%EndFieldTime = HDF5TimeInstant(CurrentInstant-1,  &
                                                                 HDF5FileX)

                        exit 

                    end if

                end do

                !Check to see the presence of only one time in the time serie
                if (HDF5FileX%StartFieldTime .eq. HDF5FileX%EndFieldTime) then
                    if (Me%StartTSTime >= HDF5FileX%StartTime .AND.                 &
                        Me%EndTSTime <= HDF5FileX%EndTime) then          
                        write(*,*) 'Time series has only one time:' 
100                 format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)
                        call ExtractDate(HDF5FileX%StartFieldTime, Year, Month,     & 
                             Day, Hour, Minute, Second)
                        write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
                        write(*,*) 'This is not allowed.'
                        stop 'OpenAndDateHDF5Files - ModuleExportHDF5ToTimeSerie - ERR70'
                    end if
                end if

                AuxNumberInstants = HDF5FileX%NumberOfInstants

                HDF5FileX%NumberOfInstants = 1 + (HDF5FileX%EndInstant -            &
                                            HDF5FileX%StartInstant)/                &
                                            (Me%DecimationFactor + 1)

                allocate(HDF5FileX%InstantsArray(1:HDF5FileX%NumberOfInstants))

                !HDF5FileX%InstantsArray =                                           &
                !    AuxInstantsArray(HDF5FileX%StartInstant:HDF5FileX%EndInstant)

                Count = 0
                do CurrentInstant = HDF5FileX%StartInstant, HDF5FileX%EndInstant,   &
                                       (Me%DecimationFactor + 1)
                    Count = Count + 1
                    HDF5FileX%InstantsArray(Count) = AuxInstantsArray(CurrentInstant)
                enddo

                !!get instant times 
                !call ObtainInstantsTimes(HDF5FileX)               

                !Calculate DT for time serie creation
                if (HDF5FileX%NumberOfInstants .ge. 2) then
                    Me%DT = HDF5FileX%InstantsArray(2) - HDF5FileX%InstantsArray(1)
                else 
                    write(*,*) 'HDF5 file with less than 2 time instants:'                
                    write(*,*) HDF5FileX%Name
                end if

                !Calculate DT for consistency checking
                if (AuxNumberInstants .ge. 3) then
                    AuxDT = AuxInstantsArray(3) - AuxInstantsArray(2)
                else
                    !Check if this DT is equal to the DT of last file: they must equal!
                    if ((AuxDT .NE. LastDT) .and. (.not. FirstTime) .and.           &
                        (AuxNumberInstants .ge. 3)) then
                        write(*,*) 'HDF5 files do not have the same DT'                
                        write(*,*) Me%FirstTSHDF5File%Name
                        write(*,*) HDF5FileX%Name
                        stop 'OpenAndDateHDF5Files - ModuleExportHDF5ToTimeSerie - ERR80'                
                    end if 
                    LastDT = AuxDT
                end if

                !Add file to list of relevant files
                call AddTSHDF5File(HDF5FileX)

                if (FirstTime) then

                    !Get water points to check time serie location after
                    !(these are assumed equal in all relevant HDF5 files: for economy no
                    !check is made)
                    if (Me%UsePointsMatrix) &
                    call ObtainWaterPoints

                    !next run of cycle (next parameter) is not the first one                
                    FirstTime = .false.
                
                endif
          
                deallocate(AuxInstantsArray)
                nullify(AuxInstantsArray)

            end if

            call KillHDF5(HDF5FileX%HDFID)           

            HDF5FileX => HDF5FileX%Next

        end do

        if (FirstTime) then

            !No HDF5 file is provided with suitable data
            write(*,*)
            call ExtractDate(Me%StartTSTime, Year, Month, Day, Hour,                &  
                             Minute, Second)        
            write(*,*) 'Data lacking from'      
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
            write(*,*) 'to'      
            call ExtractDate(Me%EndTSTime, Year, Month, Day, Hour, Minute, Second)
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
            write(*,*) 'There are not requested data in the HDF5 files provided.'      
            stop 'OpenAndDateHDF5Files - ModuleExportHDF5ToTimeSerie - ERR90'                

        else

            !assume VariableDT
            Me%VariableDT = .True.

            !For each HDF5 file needed for the time serie
            HDF5FileX => Me%FirstTSHDF5File

            do while(associated(HDF5FileX))

                !get instant times 
                call ObtainInstantsTimes(HDF5FileX)               

                !Get instants' times and put them in list 
                call AddTSInstantsTimes(HDF5FileX)

                HDF5FileX => HDF5FileX%Next            

            end do

        end if

    end subroutine OpenAndDateHDF5Files

    !--------------------------------------------------------------------------

    subroutine HDF5Selector

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX
        type(T_HDF5File), pointer                   :: OtherFile       
        integer                                     :: STAT_CALL
      
        !Begin-----------------------------------------------------------------

        !Accordingly to StartTSTime and EndTSTime select the HDF5 files to use
        HDF5FileX => Me%FirstHDF5File

        do while (associated(HDF5FileX))

            !Check if there is a file with same period
            call SamePeriodHDF5(HDF5FileX, OtherFile, STAT_CALL)

            if (STAT_CALL .EQ. SUCCESS_) then 
                write(*,*)'Same period is contained in two HDF5 files:'
                write(*,*) HDF5FileX%Name
                write(*,*) OtherFile%Name
                stop 'HDF5Selector - ModuleExportHDF5ToTimeSerie - ERR10'
            end if

            !See if the file is to be considered for time serie 
            !Check start and end time
            if ((HDF5FileX%EndTime >= Me%StartTSTime) .and.                 &
                (HDF5FileX%EndTime <= Me%EndTSTime)) then

                if (HDF5FileX%StartTime < Me%StartTSTime) then

                    !Start time serie time is after start time of file
                    HDF5FileX%StartFieldTime = Me%StartTSTime
                    HDF5FileX%EndFieldTime = HDF5FileX%EndTime

                    call AddTSHDF5File(HDF5FileX)

                else 

                    !Start time serie time is before start time of file
                    HDF5FileX%StartFieldTime = HDF5FileX%StartTime
                    HDF5FileX%EndFieldTime = HDF5FileX%EndTime

                    call AddTSHDF5File(HDF5FileX)

                end if

            else if ((HDF5FileX%StartTime >= Me%StartTSTime) .and.          &
                     (HDF5FileX%StartTime <= Me%EndTSTime)) then

                !End time serie time is before end time of file
                HDF5FileX%StartFieldTime = HDF5FileX%StartTime
                HDF5FileX%EndFieldTime = Me%EndTSTime

                call AddTSHDF5File(HDF5FileX)

            end if

            HDF5FileX => HDF5FileX%Next

        end do 

    end subroutine HDF5Selector

    !--------------------------------------------------------------------------

    subroutine HDF5Evaluator(HDF5FileX, Relevant)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX
        logical                                     :: Relevant

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_HDF5File), pointer                   :: OtherFile      

        !Begin-----------------------------------------------------------------

        Relevant = .FALSE.

        !Check if there is a file with same period
        call SamePeriodHDF5(HDF5FileX, OtherFile, STAT_CALL)

        if (STAT_CALL .EQ. SUCCESS_) then 
            write(*,*)'Same period is contained in two HDF5 files:'
            write(*,*) HDF5FileX%Name
            write(*,*) OtherFile%Name
            stop 'HDF5Evaluator - ModuleExportHDF5ToTimeSerie - ERR10'
        end if

        if (.not. Me%UseEndTSTime) then
            Me%EndTSTime = HDF5FileX%EndTime
        endif
        
        !See if the file is to be considered for time serie 
        !Check start and end time
        if ((HDF5FileX%EndTime >= Me%StartTSTime) .and.                     &
            (HDF5FileX%EndTime <= Me%EndTSTime)) then

            !End time serie time is between start and end times of file
            if (HDF5FileX%StartTime < Me%StartTSTime) then

                !Start time serie time is after start time of file
                HDF5FileX%StartFieldTime = Me%StartTSTime
                HDF5FileX%EndFieldTime = HDF5FileX%EndTime

                Relevant = .TRUE.

            else 

                !Start time serie time is before start time of file
                HDF5FileX%StartFieldTime = HDF5FileX%StartTime
                HDF5FileX%EndFieldTime = HDF5FileX%EndTime

                Relevant = .TRUE.

            end if

        else if ((HDF5FileX%StartTime >= Me%StartTSTime) .and.              &
                 (HDF5FileX%StartTime <= Me%EndTSTime)) then

            !End time serie time is before end time of file
            HDF5FileX%StartFieldTime = HDF5FileX%StartTime
            HDF5FileX%EndFieldTime = Me%EndTSTime

            Relevant = .TRUE.

        else if ((HDF5FileX%StartTime < Me%StartTSTime) .and.               &
                 (HDF5FileX%EndTime > Me%EndTSTime)) then

            !Time serie is contained in file
            HDF5FileX%StartFieldTime = Me%StartTSTime
            HDF5FileX%EndFieldTime = Me%EndTSTime

            Relevant = .TRUE.

        end if

    end subroutine HDF5Evaluator

    !--------------------------------------------------------------------------

    subroutine AddTSHDF5File(HDF5FileX)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX

        !Local-----------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileAux
        type(T_HDF5File), pointer                   :: PreviousHDF5File, LastHDF5File
      
        !Begin-----------------------------------------------------------------

        if (.not. associated(Me%FirstTSHDF5File)) then

            call CreateTSHDF5File(Me%FirstTSHDF5File, HDF5FileX)
            call CreateTSHDF5File(Me%LastTSHDF5File, HDF5FileX)
            deallocate(Me%FirstTSHDF5File%Next) 
            nullify(Me%FirstTSHDF5File%Next)

        else

            if (HDF5FileX%StartFieldTime < Me%FirstTSHDF5File%StartFieldTime) then
                !current file should be the first file in list

                !save the previous list 
                allocate(HDF5FileAux)
                call CreateTSHDF5File(HDF5FileAux, Me%FirstTSHDF5File)
                HDF5FileAux%Next => Me%FirstTSHDF5File%Next              

                !make the first element in the list of relevant files equal to current file
                call CreateTSHDF5File(Me%FirstTSHDF5File, HDF5FileX)
                Me%FirstTSHDF5File%Next => HDF5FileAux

            else
                !check next files in list

                !locate previous file in the first file
                allocate(PreviousHDF5File)
                PreviousHDF5File => Me%FirstTSHDF5File                   

                do while(associated(PreviousHDF5File))

                    if (.not. associated(PreviousHDF5File%Next)) then
        
                        !current file is the last file in the list of relevant files
                        call CreateTSHDF5File(PreviousHDF5File%Next, HDF5FileX)
                        allocate(LastHDF5File)
                        LastHDF5File => PreviousHDF5File%Next
                        deallocate(LastHDF5File%Next)
                        nullify(LastHDF5File%Next)
                        call CreateTSHDF5File(Me%LastTSHDF5File, HDF5FileX)

                        !current file was added to list
                        !Adjust end field time if equal to next file start field time
                        if ((Me%LastTSHDF5File%StartFieldTime -                     &
                            PreviousHDF5File%EndFieldTime) .lt. Me%DT) then

                            call AdjustHDF5EndInstant(PreviousHDF5File)

                        endif
                        
                        exit

                    else

                        !check if current file should be located before the next file
                        if (HDF5FileX%StartFieldTime <                              &
                            PreviousHDF5File%Next%StartFieldTime) then
                            !current file should be located before next file

                            !save the previous list begining in PreviousHDF5File%Next
                            allocate(LastHDF5File)
                            LastHDF5File => PreviousHDF5File%Next
                            allocate(HDF5FileAux)

                            call CreateTSHDF5File(HDF5FileAux, LastHDF5File)
                            HDF5FileAux%Next => LastHDF5File%Next

                            call CreateTSHDF5File(LastHDF5File, HDF5FileX)

                            PreviousHDF5File%Next => LastHDF5File
                            LastHDF5File%Next => HDF5FileAux

                            !current file was added to list

                            !Adjust end field time if equal to next file start field time
                            if ((LastHDF5File%StartFieldTime -                      &
                                HDF5FileAux%EndFieldTime) .lt. Me%DT) then

                                call AdjustHDF5EndInstant(HDF5FileAux)

                            endif

                            if ((HDF5FileAux%StartFieldTime -                       &
                                PreviousHDF5File%EndFieldTime) .lt. Me%DT) then

                                call AdjustHDF5EndInstant(PreviousHDF5File)

                            endif

                            exit

                        end if

                    end if

                    !check next file in list
                    PreviousHDF5File => PreviousHDF5File%Next     

                end do

            end if

        end if

    end subroutine AddTSHDF5File

    !--------------------------------------------------------------------------

    subroutine SamePeriodHDF5(HDF5FileX, HDF5FileAux, STAT)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX
        type(T_HDF5File), pointer                   :: HDF5FileAux
        integer, intent(OUT)                        :: STAT

        !Local-----------------------------------------------------------------
      
        !Begin-----------------------------------------------------------------

        STAT = UNKNOWN_

        HDF5FileAux => Me%FirstHDF5File

        do while (associated(HDF5FileAux))

            if (HDF5FileX%Name .NE. HDF5FileAux%Name) then 
            !(not the same HDF5 file)

                !Check if the same period is in more than one file
                if (((HDF5FileX%StartTime >= HDF5FileAux%StartTime)             &
                    .and. (HDF5FileX%StartTime < HDF5FileAux%EndTime))          &
                    .or. ((HDF5FileX%EndTime > HDF5FileAux%StartTime)           &
                    .and. (HDF5FileX%EndTime <= HDF5FileAux%EndTime))           &
                    .or. ((HDF5FileX%StartTime <= HDF5FileAux%StartTime)        &
                    .and. (HDF5FileX%EndTime >= HDF5FileAux%EndTime))) then
                    !It is allowed that the end instant of a file is equal to 
                    !the start instant of another file (for continuous runs)

                    STAT = SUCCESS_

                end if

            end if

            HDF5FileAux => HDF5FileAux%Next

        end do

    end subroutine SamePeriodHDF5

    !--------------------------------------------------------------------------

    subroutine CreateTSHDF5File(HDF5FileNew, HDF5FileX)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX
        type(T_HDF5File), pointer                   :: HDF5FileNew

        !Local-----------------------------------------------------------------
      
        !Begin-----------------------------------------------------------------

        !This subroutine atributes the values of the fields of a HDF5File to another HDF5File
        allocate(HDF5FileNew)

        HDF5FileNew%Name             =  HDF5FileX%Name
        HDF5FileNew%StartTime        =  HDF5FileX%StartTime
        HDF5FileNew%EndTime          =  HDF5FileX%EndTime
        HDF5FileNew%StartFieldTime   =  HDF5FileX%StartFieldTime
        HDF5FileNew%EndFieldTime     =  HDF5FileX%EndFieldTime
        HDF5FileNew%Rank             =  HDF5FileX%Rank
        HDF5FileNew%NumberOfInstants =  HDF5FileX%NumberOfInstants
        HDF5FileNew%StartInstant     =  HDF5FileX%StartInstant
        HDF5FileNew%EndInstant       =  HDF5FileX%EndInstant
        HDF5FileNew%InstantsArray    => HDF5FileX%InstantsArray
        HDF5FileNew%FirstInstantTime => HDF5FileX%FirstInstantTime

    end subroutine CreateTSHDF5File

    !--------------------------------------------------------------------------

    subroutine AdjustHDF5EndInstant(HDF5FileX)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX

        !Local-----------------------------------------------------------------
        integer                                     :: HDF5_READ
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        HDF5FileX%EndInstant = HDF5FileX%EndInstant - 1
        
        if (HDF5FileX%EndInstant .lt. HDF5FileX%StartInstant) then

            write(*,*)'HDF5 file with no relevant data:'
            write(*,*) trim(HDF5FileX%Name)
            write(*,*)'Remove correspondent HDF5 block from input file.'
            stop 'AdjustHDF5EndInstant - ModuleHDF5Exporter - ERR01'
        else

            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

            !Open HDF5 file
            call ConstructHDF5 (HDF5FileX%HDFID, trim(HDF5FileX%Name), HDF5_READ,   &
                                STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'AdjustHDF5EndInstant - ModuleHDF5Exporter - ERR02'

            HDF5FileX%EndFieldTime = HDF5TimeInstant(HDF5FileX%EndInstant, HDF5FileX)
            
            call killhdf5(HDF5FileX%HDFID) 

            HDF5FileX%NumberOfInstants = HDF5FileX%NumberOfInstants - 1
        endif

    end subroutine AdjustHDF5EndInstant

    !--------------------------------------------------------------------------

    subroutine AddTSInstantsTimes(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                 :: ObjHDF5File
        
        !Local-----------------------------------------------------------------
        type (T_TimeSerieTime), pointer           :: NewTime
        type (T_TimeSerieTime), pointer           :: ObjTime
        type (T_TimeSerieTime), pointer           :: PreviousTime
        type (T_TimeSerieTime), pointer           :: ObjInstantTime
        
        !Begin-----------------------------------------------------------------

        ObjInstantTime => ObjHDF5File%FirstInstantTime
        
        do while (associated(ObjInstantTime))

            !Allocates new instance
            allocate (NewTime)
            nullify  (NewTime%Next)

            NewTime%Time = ObjInstantTime%Time

            !Insert New Instance into list and makes Current point to it
            if (.not. associated(Me%FirstTimeSerieTime)) then
            !FirstField should be the same for all HDF5 files 
                Me%FirstTimeSerieTime   => NewTime
                ObjTime                 => NewTime
            else
                PreviousTime            => Me%FirstTimeSerieTime
                ObjTime                 => Me%FirstTimeSerieTime%Next
                do while (associated(ObjTime))
                    PreviousTime        => ObjTime
                    ObjTime             => ObjTime%Next
                enddo

                ObjTime                 => NewTime
                PreviousTime%Next       => NewTime

            endif

            ObjInstantTime => ObjInstantTime%Next          

        end do 

    end subroutine AddTSInstantsTimes

    !--------------------------------------------------------------------------

    subroutine ObtainInstantsTimes(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                 :: ObjHDF5File
        
        !Local-----------------------------------------------------------------
        type (T_TimeSerieTime), pointer           :: NewTime
        type (T_TimeSerieTime), pointer           :: ObjTime
        type (T_TimeSerieTime), pointer           :: PreviousTime
        integer                                   :: CurrentInstant
        integer                                   :: HDF5_READ
        integer                                   :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Open HDF5 file
        call ConstructHDF5 (ObjHDF5File%HDFID,                                  &
                            trim(ObjHDF5File%Name),                             &
                            HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ObtainInstantsTimes - ModuleAssimilationPreProcessor - ERR01'

        do CurrentInstant = ObjHDF5File%StartInstant, ObjHDF5File%EndInstant,   &
                            (Me%DecimationFactor + 1)

            !Allocates new instance
            allocate (NewTime)
            nullify  (NewTime%Next)

            NewTime%Time = HDF5TimeInstant(CurrentInstant, ObjHDF5File)

            !Insert New Instance into list and makes Current point to it
            if (.not. associated(ObjHDF5File%FirstInstantTime)) then
            !FirstField should be the same for all HDF5 files 
                ObjHDF5File%FirstInstantTime    => NewTime
                ObjTime                         => NewTime
            else
                PreviousTime                    => ObjHDF5File%FirstInstantTime
                ObjTime                         => ObjHDF5File%FirstInstantTime%Next
                do while (associated(ObjTime))
                    PreviousTime                => ObjTime
                    ObjTime                     => ObjTime%Next
                enddo
                ObjTime                         => NewTime
                PreviousTime%Next               => NewTime
            endif

        end do 

        call killhdf5(ObjHDF5File%HDFID) 

    end subroutine ObtainInstantsTimes

    !--------------------------------------------------------------------------

    subroutine ObtainWaterPoints

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer  , dimension(7)                     :: Dimensions
        integer                                     :: nItems, AuxPosition, item
        character(len=StringLength)                 :: AuxName
        logical                                     :: exist = ON
        integer                                     :: HDF5_READ, STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !open HDF5 file           
        call ConstructHDF5 (Me%FirstTSHDF5File%HDFID,                               &
                            trim(Me%FirstTSHDF5File%Name),                          &
                            HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ObtainWaterPoints - ModuleExportHDF5ToTimeSerie - ERR10'

        !check if provided name for water points data set exist
        call GetHDF5DataSetExist(Me%FirstTSHDF5File%HDFID,                          & 
                        trim(Me%WaterPointsGroup)//"/"//trim(Me%WaterPointsName),  &
                        exist, STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ObtainWaterPoints - ModuleExportHDF5ToTimeSerie - ERR20'

        if (.not. exist) then
        
            write(*,*)
            write(*,*)'Expected Water Points variable does not exist in group:'
            write(*,*) trim(Me%WaterPointsGroup)//"/"//trim(Me%WaterPointsName)
            stop 'ObtainWaterPoints - ModuleExportHDF5ToTimeSerie - ERR30'

        else
    
            !find position for data set
            call GetHDF5GroupNumberOfItems(Me%FirstTSHDF5File%HDFID,                &
                                           trim(Me%WaterPointsGroup), nItems, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            & 
                stop 'ObtainWaterPoints - ModuleExportHDF5ToTimeSerie - ERR40'

do2:        do item = 1, nItems

                call GetHDF5GroupID(Me%FirstTSHDF5File%HDFID,                       &
                                    trim(Me%WaterPointsGroup), item,                &
                                    AuxName, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        & 
                    stop 'ObtainWaterPoints - ModuleExportHDF5ToTimeSerie - ERR50'

                if (trim(AuxName) == trim(Me%WaterPointsName)) then

                    AuxPosition = item
                    exit do2

                endif
            
            end do do2                        

            !get field rank and dimensions
            call GetHDF5GroupID(Me%FirstTSHDF5File%HDFID, trim(Me%WaterPointsGroup),&
                                AuxPosition, Me%WaterPointsName,                    &
                                Rank = Me%WaterPointsRank,                          & 
                                Dimensions = Dimensions,                            &
                                STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_)                                            & 
                stop 'ObtainWaterPoints - ModuleExportHDF5ToTimeSerie - ERR60'
    
            !get water points values
            select case (Me%WaterPointsRank)

                case(2)
                ! 2D WaterPoints

                    !calculate Size2D
                    Me%Size2D%ILB = 1
                    Me%Size2D%IUB = Dimensions(1)
                    Me%Size2D%JLB = 1
                    Me%Size2D%JUB = Dimensions(2)

                    !allocate field
                    nullify (Me%WaterPoints2D)
                    allocate(Me%WaterPoints2D(Me%Size2D%ILB:Me%Size2D%IUB,          &
                                              Me%Size2D%JLB:Me%Size2D%JUB))
               
                    call HDF5SetLimits (Me%FirstTSHDF5File%HDFID, Me%Size2D%ILB,    &
                                        Me%Size2D%IUB, Me%Size2D%JLB,Me%Size2D%JUB, &
                                        STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    & 
                        stop 'ObtainWaterPoints - ModuleExportHDF5ToTimeSerie - ERR70'

                    !read field
                    call HDF5ReadData(Me%FirstTSHDF5File%HDFID,                     &
                                      trim(Me%WaterPointsGroup),                    &
                                      trim(Me%WaterPointsName),                     &
                                      Array2D      = Me%WaterPoints2D,              &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'ObtainWaterPoints - ModuleExportHDF5ToTimeSerie - ERR80'

                case(3)
                ! 3D WaterPoints

                    !calculate Size3D
                    Me%Size3D%ILB = 1
                    Me%Size3D%IUB = Dimensions(1)
                    Me%Size3D%JLB = 1
                    Me%Size3D%JUB = Dimensions(2)
                    Me%Size3D%KLB = 1
                    Me%Size3D%KUB = Dimensions(3)
                
                    !allocate field
                    nullify (Me%WaterPoints3D)
                    allocate(Me%WaterPoints3D(Me%Size3D%ILB:Me%Size3D%IUB,          &
                                               Me%Size3D%JLB:Me%Size3D%JUB,         &
                                               Me%Size3D%KLB:Me%Size3D%KUB))
                
                    call HDF5SetLimits  (Me%FirstTSHDF5File%HDFID, Me%Size3D%ILB,   &
                                         Me%Size3D%IUB, Me%Size3D%JLB,Me%Size3D%JUB, &
                                         Me%Size3D%KLB,Me%Size3D%KUB, STAT = STAT_CALL)                
                    if (STAT_CALL .NE. SUCCESS_)                                    & 
                        stop 'ObtainWaterPoints - ModuleExportHDF5ToTimeSerie - ERR90'

                    !read field
                    call HDF5ReadData(Me%FirstTSHDF5File%HDFID,                     &
                                      trim(Me%WaterPointsGroup),                    &
                                      trim(Me%WaterPointsName),                     &
                                      Array3D      = Me%WaterPoints3D,              &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'ObtainWaterPoints - ModuleExportHDF5ToTimeSerie - ERR100'

            case default 
            
                write(*,*)'Water Points variable must be 2D or 3D.'
                stop 'ObtainWaterPoints - ModuleExportHDF5ToTimeSerie - ERR110'
            
            end select

            call KillHDF5(Me%FirstTSHDF5File%HDFID)           

        end if

    end subroutine ObtainWaterPoints

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine ModifyExportHDF5ToTimeSerie

        if (Me%ExportType == ExportCellToTimeSeries) then
        
            call ModifyExportHDF5CellToTimeSerie
        
        else
        
            call ModifyExportHDF5AreaToTimeSerie
        
        endif

    end subroutine ModifyExportHDF5ToTimeSerie

    !--------------------------------------------------------------------------
    
    subroutine ModifyExportHDF5CellToTimeSerie
    
        !Local-----------------------------------------------------------------
        logical                                     :: Running
        type (T_Parameter), pointer                 :: ObjParameter, ParameterToKill
        type (T_TimeSerieTime), pointer             :: ObjTimeSerieTime
        integer                                     :: STAT_CALL
        real                                        :: DT
        integer                                     :: Count = 0
        logical                                     :: FirstFile
        type(T_HDF5File), pointer                   :: ObjHDF5File
        type(T_Time)                                :: LastStartTime
        real                                        :: LastDT

        !Begin-----------------------------------------------------------------

        ObjTimeSerieTime => Me%FirstTimeSerieTime

        !Cycle each relevant HDF5 file
        FirstFile =.true.

        Me%ActualEndTSTime = Me%EndTSTime

        !For each HDF5 file needed for the time serie
        ObjHDF5File => Me%FirstTSHDF5File
        
        !Initialize variables needed to check lacks in data
        LastStartTime = Me%StartTSTime
        LastDT = 0

        do while(associated(ObjHDF5File))

            !Open and read relevant data
            call OpenAndReadHDF5File(FirstFile, ObjHDF5File, LastDT, LastStartTime)
            !(HDF5 parameter data must be read in this cycle for memory to be released
            !after writing in the time serie)

            if (FirstFile) FirstFile = .false.

            !Cycle each parameter
            !For each parameter write Time Serie
            Running      = .true.

            Me%CurrentTime  = ObjTimeSerieTime%Time

            write(*,*)     
            write(*,*)'Writing time serie files...'

            do while (Running)

                !For each parameter write Time Serie 
                ObjParameter => Me%FirstParameter
            
                Count = Count + 1 
            
                if (Count .EQ. 1) then 

                    ObjParameter%CurrentField => ObjParameter%FirstField 

                end if

                do while(associated(ObjParameter))

                    if (.not. Me%VariableGrid) then

                        select case (ObjParameter%Rank)

                            !Values in the begining of files in middle of time serie are not 
                            !writen. These values are very close in time to last value of 
                            !previous file.

                            case(2)

                                call WriteTimeSerie(Me%ObjTimeSerie,                            &
                                                Data2D = ObjParameter%CurrentField%Values2D,    &
                                                STAT = STAT_CALL)
 
                            case(3)

                                call WriteTimeSerie(Me%ObjTimeSerie,                            &
                                                Data3D = ObjParameter%CurrentField%Values3D,    &
                                                STAT = STAT_CALL)

                            case default 
            
                                write(*,*)'Time Serie created only for 2D or 3D HDF5 parameters.'
                                stop 'ModifyExportHDF5ToTimeSerie - ModifyExportHDF5CellToTimeSerie - ERR010'
            
                        end select
                    
                    
                    elseif(Me%VariableGrid) then
                        
                        select case(ObjParameter%Rank)
                            
                            case(2)

                                if (.not.ObjParameter%CurrentField%OutPoint) then                            

                                call WriteTimeSerie(Me%ObjTimeSerie,                                &
                                                    Data2D = ObjParameter%CurrentField%Values2D,    &
                                                    Icell  = ObjParameter%CurrentField%PositionI,   &
                                                    Jcell  = ObjParameter%CurrentField%PositionJ,   &
                                                    Kcell  = ObjParameter%CurrentField%PositionK,   &
                                                    STAT = STAT_CALL)
                                 endif
                            case(3)

                                call WriteTimeSerie(Me%ObjTimeSerie,                                &
                                                    Data3D = ObjParameter%CurrentField%Values3D,    &
                                                    STAT = STAT_CALL)

                            case default 
            
                                write(*,*)'Time Serie created only for 2D or 3D HDF5 parameters.'
                                stop 'ModifyExportHDF5ToTimeSerie - ModifyExportHDF5CellToTimeSerie - ERR020'
            
                        end select

                    endif

                    ObjParameter%CurrentField => ObjParameter%CurrentField%Next               

                    ObjParameter              => ObjParameter%Next               

                end do

                if (associated(ObjTimeSerieTime%Next)) then

                    DT = ObjTimeSerieTime%Next%Time - Me%CurrentTime

                    ObjTimeSerieTime => ObjTimeSerieTime%Next

                end if

                !Actualization of time            
                Me%CurrentTime = Me%CurrentTime + DT

                call ActualizeCurrentTime(Me%ObjTime, DT, STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)                                         &
                stop 'ModifyExportHDF5ToTimeSerie - ModifyExportHDF5CellToTimeSerie - ERR030'

                !Running dependent of the last time of file
                if ((Me%CurrentTime <= ObjHDF5File%EndFieldTime) .and. (DT .ne. 0)) then
                !(if DT = 0 then Running = .false. because is the linking instant between
                !continuous run adjacent files)
                    Running = .true.
                else
                    Running = .false.
                end if

            end do

            Count = 0

            !Kill field space for each parameter 
            ObjParameter => Me%FirstParameter

            do while(associated(ObjParameter))  

                ParameterToKill => ObjParameter
                call KillIndividualParameterFields(ObjParameter)
                ObjParameter    => ObjParameter%Next

            end do

            ObjHDF5File => ObjHDF5File%Next           

        end do 

        call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                                 &
        stop 'ModifyExportHDF5ToTimeSerie - ModifyExportHDF5CellToTimeSerie - ERR040'

    end subroutine ModifyExportHDF5CellToTimeSerie
    
    !--------------------------------------------------------------------------
    
    subroutine ModifyExportHDF5AreaToTimeSerie
    
        !Local-----------------------------------------------------------------
        logical                                     :: Running
        type (T_Parameter), pointer                 :: ObjParameter, ParameterToKill
        type (T_TimeSerieTime), pointer             :: ObjTimeSerieTime
        integer                                     :: STAT_CALL
        real                                        :: DT
        integer                                     :: Count = 0
        logical                                     :: FirstFile
        type(T_HDF5File), pointer                   :: ObjHDF5File
        type(T_Time)                                :: LastStartTime
        real                                        :: LastDT
        
        real, dimension(:,:), pointer               :: mask_2d
        real, dimension(:,:,:), pointer             :: mask_3d
        type (T_TimeSeriesData), pointer            :: TS
        integer                                     :: i, j, k
        integer                                     :: parameter_index, ts_index, p_index
        real                                        :: val   
        logical                                     :: go

        !Begin-----------------------------------------------------------------

        ObjTimeSerieTime => Me%FirstTimeSerieTime

        !Cycle each relevant HDF5 file
        FirstFile =.true.

        Me%ActualEndTSTime = Me%EndTSTime

        !For each HDF5 file needed for the time serie
        ObjHDF5File => Me%FirstTSHDF5File
        
        !Initialize variables needed to check lacks in data
        LastStartTime = Me%StartTSTime
        LastDT = 0
        
        TS => Me%FirstTimeSeriesData
        do while (associated(TS))
            allocate (TS%ParamsData (Me%ParameterNumber * 7))
            
            TS%ParamsData = 0.0
            TS%Sum        = 0.0
            TS%QSum       = 0.0
            TS%Mean       = 0.0
            TS%Max        = 0.0
            TS%Min        = 0.0
            TS%SD         = 0.0
            TS%CV         = 0.0
            TS%Count      = 0
            
            TS => TS%Next
        enddo
                
        if (Me%PolygonON) then
        
            mask_2d => Me%mask_2D

        else
        
            if (Me%MaskIs3D) then
                call GetGridData (Me%ObjMaskGrid, mask_3d, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyExportHDF5ToTimeSerie - ModuleExportHDF5ToTimeSerie - ERR010'
            else
                call GetGridData (Me%ObjMaskGrid, mask_2d, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyExportHDF5ToTimeSerie - ModuleExportHDF5ToTimeSerie - ERR020'
            endif

        endif

        do while(associated(ObjHDF5File))

            !Open and read relevant data
            call OpenAndReadHDF5File(FirstFile, ObjHDF5File, LastDT, LastStartTime)
            !(HDF5 parameter data must be read in this cycle for memory to be released
            !after writing in the time serie)

            if (FirstFile) FirstFile = .false.

            !Cycle each parameter
            !For each parameter write Time Serie
            Running      = .true.

            Me%CurrentTime  = ObjTimeSerieTime%Time

            write(*,*)     
            write(*,*)'Writing time serie files...'

            do while (Running)

                !For each parameter write Time Serie 
                ObjParameter => Me%FirstParameter
                
                TS => Me%FirstTimeSeriesData
                do while (associated(TS))                
                    TS%First = .true.
                    TS       => TS%Next
                enddo
                
                parameter_index = 1            
                Count           = Count + 1  
                           
                if (Count .EQ. 1) &
                    ObjParameter%CurrentField => ObjParameter%FirstField
                
                TS => Me%FirstTimeSeriesData
                do while (associated(TS))
                    TS%ParamsData = 0.0
                    TS            => TS%Next
                enddo                
                
                do while(associated(ObjParameter))               
                              
                    !Create a timeseries for an area, instead of a point, doing the sum, mean, max, min, number of points and standard deviation                        
                    select case (ObjParameter%Rank)
                    case (2)
                        if (Me%MaskIs3D) &
                            stop 'ModifyExportHDF5ToTimeSerie - ModuleExportHDF5ToTimeSerie - ERR00'
                    
                        !First do the sum of the values, by ID of TimeSeries
                        do i = 1, Me%Size2D%IUB
                        do j = 1, Me%Size2D%JUB
                        
                            go = .false.
                            
                            if (.not. Me%UsePointsMatrix) then
                                go = .true.
                            else
                                if (Me%WaterPoints2D(i,j) == 1) then
                                    go = .true.
                                endif
                            endif
                            
                            if (go) then
                                if ((mask_2d(i,j) > 0) .and. (ObjParameter%CurrentField%Values2D(i, j) .ne. Me%AreaFillValue)) then
                                    TS => Me%FirstTimeSeriesData
do_ts1:                             do while (associated(TS))
                                        if (mask_2d(i, j) == TS%ID) then
                                            val      = ObjParameter%CurrentField%Values2D(i, j)
                                            TS%Sum   = TS%Sum + val
                                            TS%Count = TS%Count + 1
                                            if (TS%First) then
                                                TS%Max = val
                                                TS%Min = val
                                                TS%First = .false.
                                            else
                                                if (TS%Max < val) TS%Max = val
                                                if (TS%Min > val) TS%Min = val
                                            endif
                                            exit do_ts1
                                        
                                        endif
                                
                                        TS => TS%Next
                                    enddo do_ts1
                                endif
                            endif
                        enddo
                        enddo
                        
                        !Compute mean for each TS
                        TS => Me%FirstTimeSeriesData
                        do while (associated(TS))
                            if (TS%Count > 0) then
                                TS%Mean = TS%Sum / TS%Count
                            else
                                TS%Mean = 0.0
                            endif
                            
                            TS => TS%Next
                        enddo
                        
                        !Compute (x - xm)^2
                        do i = 1, Me%Size2D%IUB
                        do j = 1, Me%Size2D%JUB
                            
                            go = .false.
                            
                            if (.not. Me%UsePointsMatrix) then
                                go = .true.
                            else
                                if (Me%WaterPoints2D(i,j) == 1) then
                                    go = .true.
                                endif
                            endif
                            
                            if (go) then                        
                                if ((mask_2d(i,j) > 0)  .and. (ObjParameter%CurrentField%Values2D(i, j) .ne. Me%AreaFillValue)) then
                                    TS => Me%FirstTimeSeriesData
do_ts2:                             do while (associated(TS))
                                        if (mask_2d(i, j) == TS%ID) then
                                            val      = ObjParameter%CurrentField%Values2D(i, j)
                                            TS%QSum  = TS%QSum + (val - TS%Mean)**2
                                            exit do_ts2
                                        
                                        endif
                                
                                        TS => TS%Next
                                    enddo do_ts2
                                endif
                            endif                            
                        enddo
                        enddo            

                        !Compute SD for each TS
                        TS => Me%FirstTimeSeriesData
                        do while (associated(TS))
                            if (TS%Count > 0) then
                                TS%SD = sqrt(TS%QSum / (TS%Count - 1))
                                TS%CV = TS%SD / abs(TS%Mean)
                            else
                                TS%SD = 0.0
                                TS%CV = 0.0
                            endif
                            
                            TS => TS%Next
                        enddo
                        
                    case (3)
                    
                        if (Me%MaskIs3D) then

                            !First do the sum of the values, by ID of TimeSeries
                            do i = 1, Me%Size3D%IUB
                            do j = 1, Me%Size3D%JUB
                            do k = 1, Me%Size3D%KUB
                                
                                go = .false.
                            
                                if (.not. Me%UsePointsMatrix) then
                                    go = .true.
                                else
                                    if (Me%WaterPoints3D(i,j,k) == 1) then
                                        go = .true.
                                    endif
                                endif
                            
                                if (go) then                                  
                                    if ((mask_3d(i,j,k) > 0)  .and. (ObjParameter%CurrentField%Values3D(i, j, k) .ne. Me%AreaFillValue)) then
                                        TS => Me%FirstTimeSeriesData
do_ts3:                                 do while (associated(TS))
                                            if (mask_3d(i,j,k) == TS%ID) then
                                                val      = ObjParameter%CurrentField%Values3D(i,j,k)
                                                TS%Sum   = TS%Sum + val
                                                TS%Count = TS%Count + 1
                                                if (TS%First) then
                                                    TS%Max = val
                                                    TS%Min = val
                                                    TS%First = .false.
                                                else
                                                    if (TS%Max < val) TS%Max = val
                                                    if (TS%Min > val) TS%Min = val
                                                endif
                                                exit do_ts3                                        
                                            endif
                                    
                                            TS => TS%Next
                                        enddo do_ts3
                                    endif
                                endif                                
                            enddo
                            enddo
                            enddo
                            
                            !Compute mean for each TS
                            TS => Me%FirstTimeSeriesData
                            do while (associated(TS))
                                if (TS%Count > 0) then
                                    TS%Mean = TS%Sum / TS%Count
                                else
                                    TS%Mean = 0.0
                                endif
                                
                                TS => TS%Next
                            enddo
                            
                            !Compute (x - xm)^2
                            do i = 1, Me%Size3D%IUB
                            do j = 1, Me%Size3D%JUB
                            do k = 1, Me%Size3D%KUB
                                go = .false.
                            
                                if (.not. Me%UsePointsMatrix) then
                                    go = .true.
                                else
                                    if (Me%WaterPoints3D(i,j,k) == 1) then
                                        go = .true.
                                    endif
                                endif
                            
                                if (go) then                                  
                                    if ((mask_3d(i,j,k) > 0) .and. (ObjParameter%CurrentField%Values3D(i, j, k) .ne. Me%AreaFillValue)) then
                                        TS => Me%FirstTimeSeriesData
do_ts4:                                 do while (associated(TS))
                                            if (mask_3d(i,j,k) == TS%ID) then
                                                val      = ObjParameter%CurrentField%Values3D(i,j,k)
                                                TS%QSum  = TS%QSum + (val - TS%Mean)**2
                                                exit do_ts4
                                            
                                            endif
                                    
                                            TS => TS%Next
                                        enddo do_ts4
                                    endif
                                endif                                
                            enddo
                            enddo
                            enddo            

                            !Compute SD for each TS
                            TS => Me%FirstTimeSeriesData
                            do while (associated(TS))
                                if (TS%Count > 0) then
                                    TS%SD = sqrt(TS%QSum / (TS%Count - 1))
                                else
                                    TS%QSum = 0.0
                                endif
                                
                                TS => TS%Next
                            enddo
                                                    
                        else
                        
                            !First do the sum of the values, by ID of TimeSeries
                                                               
                            TS => Me%FirstTimeSeriesData
                            do while (associated(TS))

                                if (TS%Layer < 0) then
                                    k = Me%Size3D%KUB
                                else
                                    k = TS%Layer
                                endif
                                           
                                do i = 1, Me%Size3D%IUB
                                do j = 1, Me%Size3D%JUB  
                                
                                    go = .false.
                                    if (.not. Me%UsePointsMatrix) then
                                        go = .true.
                                    elseif (Me%WaterPoints3D(i,j,k) == 1) then
                                        go = .true.
                                    endif                                            
                                
                                    if (go) then                                                              
                                        if ((mask_2d(i,j) > 0) .and. (ObjParameter%CurrentField%Values3D(i, j, k) .ne. Me%AreaFillValue)) then
                                            if (mask_2d(i,j) == TS%ID) then
                                                val      = ObjParameter%CurrentField%Values3D(i,j,k)
                                                TS%Sum   = TS%Sum + val
                                                TS%Count = TS%Count + 1
                                                if (TS%First) then
                                                    TS%Max = val
                                                    TS%Min = val
                                                    TS%First = .false.
                                                else
                                                    if (TS%max < val) TS%Max = val
                                                    if (TS%min > val) TS%Min = val
                                                endif                                       
                                            endif
                                        endif
                                    endif
                                                                            
                                enddo
                                enddo
                                

                                !Compute mean for each TS
                                if (TS%Count > 0) then
                                    TS%Mean = TS%Sum / TS%Count
                                else
                                    TS%Mean = 0.0
                                endif
                                  
                                !Compute (x - xm)^2
                                do i = 1, Me%Size3D%IUB
                                do j = 1, Me%Size3D%JUB  
                                
                                    go = .false.
                                    if (.not. Me%UsePointsMatrix) then
                                        go = .true.
                                    elseif (Me%WaterPoints3D(i,j,k) == 1) then
                                        go = .true.
                                    endif
                            
                                    if (go) then                                                              
                                
                                        if ((mask_2d(i,j) > 0) .and. (ObjParameter%CurrentField%Values3D(i, j, k) .ne. Me%AreaFillValue)) then

                                            if (mask_2d(i,j) == TS%ID) then
                                                val      = ObjParameter%CurrentField%Values3D(i,j,k)
                                                TS%QSum  = TS%QSum + (val - TS%Mean)**2
                                            endif                                            
                                        
                                        endif
                                    endif    
                                    
                                enddo
                                enddo  
                                
                                if (TS%Count > 0) then
                                    TS%SD = sqrt(TS%QSum / (TS%Count - 1))
                                    TS%CV = TS%SD / abs(TS%Mean)
                                else
                                    TS%SD = 0.0
                                    TS%CV = 0.0                                    
                                endif
                                
                                TS => TS%Next
                            enddo                        
                        endif
                                        
                    end select

                    ts_index = 1
                    p_index  = (parameter_index - 1) * 7
                    TS => Me%FirstTimeSeriesData
                    do while (associated(TS))
                        TS%ParamsData (p_index + 1) = TS%Sum
                        TS%ParamsData (p_index + 2) = TS%Mean
                        TS%ParamsData (p_index + 3) = TS%Min
                        TS%ParamsData (p_index + 4) = TS%Max
                        TS%ParamsData (p_index + 5) = TS%SD
                        TS%ParamsData (p_index + 6) = TS%CV
                        TS%ParamsData (p_index + 7) = TS%Count
                            
                        TS%First = .true.
                        TS%Sum   = 0.0
                        TS%QSum  = 0.0
                        TS%Mean  = 0.0
                        TS%Max   = 0.0
                        TS%Min   = 0.0
                        TS%SD    = 0.0
                        TS%CV    = 0.0
                        TS%Count = 0
                        TS       => TS%Next
                        ts_index = ts_index + 1
                    enddo
                   

                    ObjParameter%CurrentField => ObjParameter%CurrentField%Next               
                    ObjParameter              => ObjParameter%Next               
                    
                    parameter_index = parameter_index + 1
                                                                   
                enddo

                TS       => Me%FirstTimeSeriesData                
                ts_index = 1
                do while (associated(TS))
                    call WriteSpecificTimeSerieLine(Me%ObjTimeSerie, ts_index, TS%ParamsData, CheckTime = .false., STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ModifyExportHDF5ToTimeSerie - ModuleExportHDF5ToTimeSerie - ERR030'
                        
                    TS%ParamsData = 0.0
                    TS            => TS%Next
                    ts_index      = ts_index + 1
                enddo                

                if (associated(ObjTimeSerieTime%Next)) then
                    DT               = ObjTimeSerieTime%Next%Time - Me%CurrentTime
                    ObjTimeSerieTime => ObjTimeSerieTime%Next
                end if

                !Actualization of time            
                Me%CurrentTime = Me%CurrentTime + DT

                call ActualizeCurrentTime(Me%ObjTime, DT, STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)                                         &
                stop 'ModifyExportHDF5ToTimeSerie - ModuleExportHDF5ToTimeSerie - ERR40'

                !Running dependent of the last time of file
                if ((Me%CurrentTime <= ObjHDF5File%EndFieldTime) .and. (DT .ne. 0)) then
                !(if DT = 0 then Running = .false. because is the linking instant between
                !continuous run adjacent files)
                    Running = .true.
                else
                    Running = .false.
                end if

            end do

            Count = 0

            !Kill field space for each parameter 
            ObjParameter => Me%FirstParameter
            do while(associated(ObjParameter))  
                ParameterToKill => ObjParameter
                call KillIndividualParameterFields(ObjParameter)
                ObjParameter    => ObjParameter%Next
            end do

            ObjHDF5File => ObjHDF5File%Next           

        end do 
        
        if (Me%PolygonON) then
        
            nullify(mask_2d)
            deallocate(Me%mask_2d)
            nullify(mask_2d)

        else

            if (Me%MaskIs3D) then
                call UngetGridData  (Me%ObjMaskGrid, mask_3d, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyExportHDF5ToTimeSerie - ModuleExportHDF5ToTimeSerie - ERR050'
            else
                call UngetGridData  (Me%ObjMaskGrid, mask_2d, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyExportHDF5ToTimeSerie - ModuleExportHDF5ToTimeSerie - ERR060'
            endif

        endif

        TS => Me%FirstTimeSeriesData
        do while (associated(TS))
            deallocate (TS%ParamsData)
            TS => TS%Next
        enddo

        call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                                 &
        stop 'ModifyExportHDF5ToTimeSerie - ModuleExportHDF5ToTimeSerie - ERR070'
    
    end subroutine ModifyExportHDF5AreaToTimeSerie    
    
    !--------------------------------------------------------------------------

    subroutine  OpenAndReadHDF5File(FirstFile, ObjHDF5File, LastDT, LastStartTime)

        !Arguments-------------------------------------------------------------
        logical                                     :: FirstFile
        type(T_HDF5File), pointer                   :: ObjHDF5File
        real                                        :: LastDT
        type(T_Time)                                :: LastStartTime

        !Local-----------------------------------------------------------------
        type(T_Parameter), pointer                  :: ObjParameter
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_READ
        real                                        :: Year, Month, Day, Hour 
        real                                        :: Minute, Second 

        !Begin-----------------------------------------------------------------

        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Open HDF5 file
        call ConstructHDF5 (ObjHDF5File%HDFID, trim(ObjHDF5File%Name),              & 
                            HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) then
            write(*,*)
            write(*,*) 'HDF5 file cannot be opened'//ObjHDF5File%Name                
            stop 'OpenAndReadHDF5File - ModuleExportHDF5ToTimeSerie - ERR10'
        end if

        !Check if there are time periods lacking
        if (ObjHDF5File%StartFieldTime > (LastStartTime + LastDT)) then
            !which DT should be? The new or the last? Suppose the last.
100         format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)
            call ExtractDate(LastStartTime, Year, Month, Day, Hour, Minute, Second)
            write(*,*)
            write(*,*) 'Data lacking from'      
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
            write(*,*) 'to'      
            call ExtractDate(ObjHDF5File%StartFieldTime, Year, Month, Day, Hour,    & 
                             Minute, Second)
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
        end if

        !Read fields data
        write(*,*)'Reading data from HDF5 file:'//trim(ObjHDF5File%Name) 

        !getting parameter fields data
        ObjParameter => Me%FirstParameter

        do while(associated(ObjParameter))

            if (FirstFile) then
            !it is the first file accessed

                nullify(ObjParameter%FirstField)
                !nullify have to be here because for each new parameter has to be 
                !nullified for the first file
            
            end if

            !Allocates first field for this parameter
            allocate (ObjParameter%FirstField)
            nullify (ObjParameter%FirstField)

            call ReadParameterFields(ObjHDF5File, ObjParameter)

            ObjParameter%CurrentField => ObjParameter%FirstField

            ObjParameter => ObjParameter%Next

        end do

        !in variable grids it needs to calculate i,j position for each time instant
        !which is the same for each parameter thus it can do it just once using the 
        !first parameter
        if (Me%VariableGrid) then
            call ComputeIJonVariableGrids(ObjHDF5File, Me%FirstParameter)
        endif


        LastStartTime = ObjHDF5File%EndFieldTime 
        LastDT = Me%DT

       !check if there are data lacking after the last data from last file 
       if (.not. associated(ObjHDF5File%Next)) then
           !(last file)
           if (ObjHDF5File%EndFieldTime < (Me%EndTSTime)) then
                !no DT is considered here
                call ExtractDate(ObjHDF5File%EndFieldTime, Year, Month, Day, Hour,  &  
                                 Minute, Second)
                write(*,*)
                write(*,*) 'Data lacking from'      
                write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
                write(*,*) 'to'      
                call ExtractDate(Me%EndTSTime, Year, Month, Day, Hour, Minute, Second)
                write(*,fmt=100) Year, Month, Day, Hour, Minute, Second                    
            end if
        end if

        call killhdf5(ObjHDF5File%HDFID)

    end subroutine OpenAndReadHDF5File

    !--------------------------------------------------------------------------
    
    subroutine  ComputeIJonVariableGrids (ObjHDF5File, ObjParameter)

         !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: ObjHDF5File
        type(T_Parameter), pointer                  :: ObjParameter
          
        !Local-----------------------------------------------------------------
        integer                                     :: CurrentInstant
        type (T_Field), pointer                     :: CurrentField
        integer, dimension(7)                       :: Dimensions
        integer                                     :: STAT_CALL
        character(len=StringLength)                 :: LatName, LonName
        character(len=StringLength)                 :: Units
        integer                                     :: Rank
        real(4), dimension(:,:), pointer            :: LatArray
        real(4), dimension(:,:), pointer            :: LonArray
        integer                                     :: ILB,IUB,JLB,JUB
        integer                                     :: NumberOfTimeSeries
        integer                                     :: i
        real                                        :: Lat, Long
        integer                                     :: positionI, positionJ, positionK
        
        !Begin-----------------------------------------------------------------

        CurrentField => ObjParameter%FirstField
        CurrentInstant =0

        do while(associated(CurrentField))

            CurrentInstant = CurrentInstant + 1

            !get field ID, Rank and Dimensions
            !(this needs to be done for each instant)
            call GetHDF5GroupID(ObjHDF5File%HDFID, ObjParameter%Group,                  &
                            CurrentInstant, CurrentField%Name,                          &
                            CurrentField%Units, ObjParameter%Rank,                      &
                            Dimensions,                                                 &
                            STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_)                                                &  
                stop 'ComputeIJonVariableGrids - ModuleExportHDF5ToTimeSerie - ERR10'

            call GetHDF5GroupID(ObjHDF5File%HDFID, '/Grid/Latitude',                    &
                                CurrentInstant, LatName,Units, Rank,                    &
                                Dimensions,                                             &
                                STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_)                                                &  
                stop 'ComputeIJonVariableGrids - ModuleExportHDF5ToTimeSerie - ERR20'

            call GetHDF5GroupID(ObjHDF5File%HDFID, '/Grid/Longitude',                   &
                                CurrentInstant, LonName,Units, Rank,                    &
                                Dimensions,                                             &
                                STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_)                                                &  
                stop 'ComputeIJonVariableGrids - ModuleExportHDF5ToTimeSerie - ERR30'           

            ILB = 1
            IUB = Dimensions(1)
            JLB = 1
            JUB = Dimensions(2)

            nullify (LatArray)
            nullify (LonArray)
            allocate(LatArray(ILB:IUB,JLB:JUB))
            allocate(LonArray(ILB:IUB,JLB:JUB))
    
            call HDF5SetLimits (ObjHDF5File%HDFID,ILB,IUB,JLB,JUB, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                & 
                stop 'ComputeIJonVariableGrids - ModuleExportHDF5ToTimeSerie - ERR40'

            !read  Latitude field
            call HDF5ReadData(ObjHDF5File%HDFID, '/Grid/Latitude',                      &
                            trim(LatName),                                              &
                            Array2D      = LatArray,                                    &
                            STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'ComputeIJonVariableGrids - ModuleExportHDF5ToTimeSerie - ERR505'

            call HDF5SetLimits (ObjHDF5File%HDFID,ILB,IUB,JLB,JUB, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                & 
                stop 'ComputeIJonVariableGrids - ModuleExportHDF5ToTimeSerie - ERR60'
                
            !read  Longitude field
            call HDF5ReadData(ObjHDF5File%HDFID, '/Grid/Longitude',                     &
                            trim(LonName),                                              &
                            Array2D      = LonArray,                                    &
                            STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'ComputeIJonVariableGrids - ModuleExportHDF5ToTimeSerie - ERR70'

                        
            call GetNumberOfTimeSeries(Me%ObjTimeSerie, NumberOfTimeSeries,             &
                STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'ComputeIJonVariableGrids - ModuleExportHDF5ToTimeSerie - ERR80'
                            
            nullify (CurrentField%PositionJ)
            nullify (CurrentField%PositionI)
            nullify (CurrentField%PositionK)

            allocate(CurrentField%PositionJ(1:NumberOfTimeSeries))
            allocate(CurrentField%PositionI(1:NumberOfTimeSeries))
            allocate(CurrentField%PositionK(1:NumberOfTimeSeries))

            CurrentField%PositionK = 1

            do i=1,NumberOfTimeSeries

                call GetTimeSerieLocation(Me%ObjTimeSerie, i, positionI,                &
                                               positionJ, positionK,                    &
                                               Latitude = Lat, Longitude = Long,        &
                                               STAT = STAT_CALL )
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'ComputeIJonVariableGrids - ModuleExportHDF5ToTimeSerie - ERR90'
                                    
                call ComputeIJ (LatArray, LonArray, Lat, Long ,                         &
                                ILB, IUB, JLB, JUB,                                     &
                                CurrentField%PositionI(i), CurrentField%PositionJ(i))
            
                ! se uma das series temporais tiver fora do dominio então não
                ! escreve nada deste fiels incluindo as outras series que podem estar dentro)
            
                 if ((CurrentField%PositionI(i).eq.-1).or.(CurrentField%PositionJ(i).eq.-1)) then
                    CurrentField%OutPoint=.true.
                 endif
            end do

            deallocate(LatArray)
            deallocate(LonArray)
            nullify   (LatArray)
            nullify   (LonArray)

            CurrentField => CurrentField%Next

        end do

    end subroutine ComputeIJonVariableGrids

    !--------------------------------------------------------------------------

    subroutine ComputeIJ (LatArray,   LonArray,  Latitude,   Longitude, &
                          ILB, IUB,  JLB,  JUB, PositionI,  PositionJ)

        !Arguments-------------------------------------------------------------
        integer, intent (OUT)                           :: PositionI, PositionJ
        real, intent (IN)                               :: Latitude, Longitude
        integer, intent (IN)                            :: ILB, IUB,  JLB,  JUB
        real(4), dimension(:,:), pointer                :: LatArray
        real(4), dimension(:,:), pointer                :: LonArray
 
        !Local-----------------------------------------------------------------
        type(T_Polygon),  pointer                       :: OriginalArea
        type(T_PointF),   pointer                       :: AuxPoint
        
        type     T_SelectedArea
            type(T_PointF)                              :: PointF
            type(T_Point)                               :: Point
        end type T_SelectedArea 

        type(T_SelectedArea), dimension(:),  pointer    :: UserArea

        logical                                         :: pointin
        real                                            :: lllat,ullat,urlat,lrlat
        real                                            :: urlon,ullon,lrlon,lllon 
        !real                                            :: aux
        integer                                         :: column, line

        !Begin-----------------------------------------------------------------

        allocate(OriginalArea)
        allocate(OriginalArea%VerticesF(1:5))
        nullify (UserArea)
        allocate(UserArea(1:5))
        nullify (Auxpoint)
        allocate(Auxpoint)
        
        OriginalArea%Count = 5

        OriginalArea%VerticesF(1)%X = LonArray(ilb,jlb) 
        OriginalArea%VerticesF(2)%X = LonArray(iub,jlb) 
        OriginalArea%VerticesF(3)%X = LonArray(iub,jub) 
        OriginalArea%VerticesF(4)%X = LonArray(ilb,jub)
        OriginalArea%VerticesF(5)%X = LonArray(ilb,jlb) 

        OriginalArea%VerticesF(1)%Y = LatArray(ilb,jlb)
        OriginalArea%VerticesF(2)%Y = LatArray(iub,jlb)
        OriginalArea%VerticesF(3)%Y = LatArray(iub,jub)  
        OriginalArea%VerticesF(4)%Y = LatArray(ilb,jub)
        OriginalArea%VerticesF(5)%Y = LatArray(ilb,jlb)

        lllat=LatArray(iub,1)
        ullat=LatArray(1,1)
        urlat=LatArray(1,jub)
        lrlat=LatArray(iub,jub)

        urlon=LonArray(1,jub)
        ullon=LonArray(1,1)
        lrlon=LonArray(iub,jub)
        lllon=LonArray(iub,1)
        
    
    
        ! --- Confirm that the selected point is inside the original one

        call SetLimits(OriginalArea)
        
        Auxpoint%X = Longitude
        Auxpoint%Y = Latitude
   
        pointin = IsPointInsidePolygon(Auxpoint, OriginalArea)
        
        if (.not.pointin) then
           write(*,*) "Point outside the Grid!"
        endif 
                     
        !isto é para ver qual é a coluna tento em conta que o j ta dentro dos limites
        !a linha 1 é a mais a norte!!
        
        if (pointin) then

            Auxpoint%Y=(min(ullat,urlat)+max(lllat,lrlat))/2   
        
            do column = jlb,jub-1
 
                OriginalArea%VerticesF(1)%X = LonArray(ilb, column) 
                OriginalArea%VerticesF(2)%X = LonArray(ilb, column +1) 
                OriginalArea%VerticesF(3)%X = LonArray(iub, column +1) 
                OriginalArea%VerticesF(4)%X = LonArray(iub, column)
                OriginalArea%VerticesF(5)%X = LonArray(ilb, column)            

                OriginalArea%VerticesF(1)%Y = LatArray(ilb, column) 
                OriginalArea%VerticesF(2)%Y = LatArray(ilb, column +1)
                OriginalArea%VerticesF(3)%Y = LatArray(iub, column +1)  
                OriginalArea%VerticesF(4)%Y = LatArray(iub, column)
                OriginalArea%VerticesF(5)%Y = LatArray(ilb, column)   

               call SetLimits(OriginalArea)

                           
               if (IsPointInsidePolygon(Auxpoint, OriginalArea)) then
                       PositionJ = column
                       exit
               endif

            enddo

            !isto é para ver qual é a linha tendo em conta que o i ta dentro dos limites
            Auxpoint%X=(min(urlon,lrlon)+max(lllon,ullon))/2

            do line = ilb, iub-1
 
                OriginalArea%VerticesF(1)%X = LonArray(line  ,jlb) 
                OriginalArea%VerticesF(2)%X = LonArray(line+1,jlb) 
                OriginalArea%VerticesF(3)%X = LonArray(line+1,jub) 
                OriginalArea%VerticesF(4)%X = LonArray(line  ,jub)
                OriginalArea%VerticesF(5)%X = LonArray(line  ,jlb) 

                OriginalArea%VerticesF(1)%Y = LatArray(line  ,jlb) 
                OriginalArea%VerticesF(2)%Y = LatArray(line+1,jlb)
                OriginalArea%VerticesF(3)%Y = LatArray(line+1,jub)  
                OriginalArea%VerticesF(4)%Y = LatArray(line  ,jub)
                OriginalArea%VerticesF(5)%Y = LatArray(line  ,jlb) 


                call SetLimits(OriginalArea)

                if (IsPointInsidePolygon(Auxpoint, OriginalArea)) then
                    PositionI= line
                    exit
                endif

            enddo

        else

            PositionI = -1
            PositionJ = -1

        endif

    end subroutine ComputeIJ

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine KillExportHDF5ToTimeSerie

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX, HDF5FileToKill
        type(T_Parameter), pointer                  :: ParameterX, ParameterToKill
        integer                                     :: STAT_CALL
        type (T_TimeSeriesData), pointer            :: TSX, TSToKill
        
        !Begin-----------------------------------------------------------------

        !Kill all HDF5Files from the list of the ones relevant for Time Series
        HDF5FileX=> Me%FirstTSHDF5File

        do while(associated(HDF5FileX))  

            HDF5FileToKill  => HDF5FileX
            HDF5FileX       => HDF5FileX%Next
            call KillIndividualHDF5File(HDF5FileToKill)

        end do
        nullify(Me%FirstTSHDF5File)
 
        if (Me%ExportType == ExportAreaToTimeseries) then

            if (.not.Me%PolygonON) then
                
                call KillGridData(Me%ObjMaskGrid, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    print *, 'WARN: KillExportHDF5ToTimeSerie - ModuleExportHDF5ToTimeSerie - WRN010' 
                    
                call KillHorizontalGrid(Me%ObjMaskHorizontalGrid, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    print *, 'WARN: KillExportHDF5ToTimeSerie - ModuleExportHDF5ToTimeSerie - WRN020'
            
            endif
            
            if (associated(Me%WaterPoints2D)) then    
                deallocate(Me%WaterPoints2D)
                nullify(Me%WaterPoints2D)
            endif
            
            if (associated(Me%WaterPoints3D)) then                
                deallocate(Me%WaterPoints3D)
                nullify(Me%WaterPoints3D)
            endif
            
            TSX => Me%FirstTimeSeriesData
            do while (associated(TSX))
                TSToKill => TSX
                TSX      => TSX%Next
                
                if (associated(TSToKill%ParamsData)) then
                    deallocate (TSToKill%ParamsData)
                    nullify (TSToKill%ParamsData)
                endif
                
                deallocate (TSToKill)
                nullify (TSToKill)
            enddo 
            
        endif

        !Kill all parameters from the list for Time Series
        ParameterX=> Me%FirstParameter
        do while(associated(ParameterX))  
            ParameterToKill => ParameterX
            ParameterX      => ParameterX%Next
            call KillIndividualParameter(ParameterToKill)
        end do
        nullify(Me%FirstParameter)

        deallocate(Me)
        nullify(Me)

        !----------------------------------------------------------------------

    end subroutine KillExportHDF5ToTimeSerie    

    !--------------------------------------------------------------------------

    subroutine KillIndividualHDF5File(HDF5ToDispose)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                       :: HDF5ToDispose

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Deallocate the InstantsArray
        if (associated(HDF5ToDispose%InstantsArray)) then 
            deallocate(HDF5ToDispose%InstantsArray, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'KillIndividualHDF5File - ModuleExportHDF5ToTimeSerie - ERR10'               
            nullify(HDF5ToDispose%InstantsArray)
        end if

        deallocate(HDF5ToDispose, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
        stop 'KillIndividualHDF5File - ModuleExportHDF5ToTimeSerie - ERR20'
        nullify(HDF5ToDispose)

    end subroutine KillIndividualHDF5File

    !--------------------------------------------------------------------------

    subroutine KillIndividualParameter(ParameterToDispose)

        !Arguments-------------------------------------------------------------
        type(T_Parameter), pointer                  :: ParameterToDispose

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !The fields of the parameters have already been killed
        deallocate(ParameterToDispose, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
        stop 'KillIndividualParameter - ModuleExportHDF5ToTimeSerie - ERR10'
        nullify(ParameterToDispose)

    end subroutine KillIndividualParameter

    !--------------------------------------------------------------------------

    subroutine KillIndividualParameterFields(ParameterToDispose)

        !Arguments-------------------------------------------------------------
        type(T_Parameter), pointer                  :: ParameterToDispose

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_Field), pointer                      :: FieldToKill, CurrentField

        !Begin-----------------------------------------------------------------

        CurrentField => ParameterToDispose%FirstField

        do while(associated(CurrentField))

            FieldToKill => CurrentField
            CurrentField => CurrentField%Next

            if (associated(FieldToKill%Values2D)) then
                deallocate(FieldToKill%Values2D, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'KillIndividualParameterFields - ModuleExportHDF5ToTimeSerie - ERR10'  
                nullify(FieldToKill%Values2D)
            end if

            if (associated(FieldToKill%Values3D)) then
                deallocate(FieldToKill%Values3D, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'KillIndividualParameterFields - ModuleExportHDF5ToTimeSerie - ERR20'  
                nullify(FieldToKill%Values3D)
            end if

            if (associated(FieldToKill%PositionJ)) then
                deallocate(FieldToKill%PositionJ, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'KillIndividualParameterFields - ModuleExportHDF5ToTimeSerie - ERR30'  
                nullify(FieldToKill%PositionJ)
            end if
            
            if (associated(FieldToKill%PositionI)) then
                deallocate(FieldToKill%PositionI, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'KillIndividualParameterFields - ModuleExportHDF5ToTimeSerie - ERR40'  
                nullify(FieldToKill%PositionJ)
            end if    

        end do 

    end subroutine KillIndividualParameterFields

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

end module ModuleExportHDF5ToTimeSerie









