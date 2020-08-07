!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : FillMatrix
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to Fill Matrixes
!
!------------------------------------------------------------------------------
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License 
!version 2, as published by the Free Software Foundation.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!
!------------------------------------------------------------------------------

Module ModuleFillMatrix
    
#ifndef _NOT_IEEE_ARITHMETIC
    use ieee_arithmetic
#endif
  
    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions,        only : InterpolateValueInTime, InterpolateProfile,      &
                                       SetMatrixValue, InterpolateMatrix2DInTime,       &
                                       InterpolateMatrix3DInTime,                       &
                                       InterpolateAngle2DInTime,                        &
                                       InterpolateAngle3DInTime,                        &
                                       InterpolateLinearyMatrix2D,                      &
                                       InterpolateLinearyMatrix3D,                      &
                                       AngleFromFieldToGrid,                            &
                                       WaveLengthHuntsApproximation
    use ModuleDrawing
    use ModuleBoxDif,           only : StartBoxDif, GetBoxes, GetNumberOfBoxes,         &
                                       UngetBoxDif, KillBoxDif
    use ModuleGridData,         only : ConstructGridData, GetGridData, UnGetGridData,   &
                                       KillGridData, GetGridDataType 
    use ModuleHorizontalGrid,   only : GetGridAngle, GetHorizontalGridSize,             &
                                       GetGridBorderLimits, GetLatitudeLongitude,       &
                                       GetDDecompOpenBorders,                           &
                                       GetDDecompParameters,                            &
                                       GetDDecompWorkSize2D, GetZCoordinates,           &
                                       UnGetHorizontalGrid, RotateAngleFieldToGrid,     &
                                       RotateVectorFieldToGrid, GetCheckDistortion,     &
                                       GetGridOutBorderCartLimits,                      &
                                       GetGridBorderCartPolygon,                        &
                                       GetHorizontalGrid, ConstructHorizontalGrid,      &
                                       KillHorizontalGrid, GetDDecompON, GetCellRotation
    use ModuleTimeSerie,        only : StartTimeSerieInput, GetTimeSerieValue,          &
                                       GetTimeSerieDTForNextEvent,                      &
                                       GetTimeSerieTimeOfNextDataset,                   & 
                                       GetTimeSerieTimeOfDataset,                       &
                                       GetTimeSerieDataValues,GetTimeSerieValueForIndex,&
                                        KillTimeSerie
    use ModuleGeometry,         only : GetGeometryDistances, UnGetGeometry, GetGeometrySize
    use ModuleHDF5,             only : ConstructHDF5, HDF5ReadData, GetHDF5GroupID,     &
                                       GetHDF5FileAccess, GetHDF5GroupNumberOfItems,    &
                                       HDF5SetLimits, GetHDF5ArrayDimensions, KillHDF5, &
                                       GetHDF5GroupExist, GetHDF5DataSetExist,          &
                                       GetHDF5AllDataSetsOK
                                       
    use ModuleField4D,          only : ConstructField4D, GetField4DNumberOfInstants,    &
                                       GetField4DInstant, ModifyField4D,                &
                                       ModifyField4DXYZ, GetField4DGeneric4DValue,      &
                                       KillField4D
    use ModuleStopWatch,        only : StartWatch, StopWatch
                                       

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructFillMatrix
    private ::      AllocateInstance
    private ::      ReadOptions
    private ::          ConstructProfileTSDefault
    private ::          ConstructSpaceLayers
    private ::          ConstructSpaceBox
    private ::          ConstructSpaceASCIIFile
    private ::          ConstructSpaceProfile
    private ::          ConstructAnalyticProfile
    private ::          ConstructSpaceTimeSerie
    private ::          ConstructProfileTimeSerie
    private ::          ConstructHDFInput
    private ::              HDF5TimeInstant 
    private ::              ReadHDF5Values2D
    private ::              ReadHDF5Values3D
    private ::              ProfileTimeSerieField

    !Selector
    public  :: GetIfMatrixRemainsConstant
    public  :: GetDefaultValue
    public  :: GetFillMatrixDTPrediction
    public  :: GetHDFTimeLimits
    public  :: GetNumberOfInstants
    public  :: GetTimeInstant
    public  :: GetNextValueForDTPred
    public  :: GetValuesProcessingOptions
    public  :: GetMultiTimeSeries
    public  :: UngetFillMatrix
    public  :: GetAngleField
    public  :: GetVectorialField
    public  :: GetAnalyticCelerityON
    public  :: GetAnalyticCelerity
    public  :: GetFilenameHDF
                     
    
    !Modifier
    public  :: ModifyFillMatrix
    public  :: ModifyFillMatrixVectorial
    private ::      ModifySpaceTimeSerie
    private ::      ModifyHDFInput2D
    private ::          ModifyHDFInput2DTime
    private ::          ModifyHDFInput2DGeneric4D
    private ::      ModifyHDFInput3D
    private ::          ModifyHDFInput3DTime
    private ::          ModifyHDFInput3DGeneric4D
    private ::      ModifyProfileTimeSerie
   

    !Destructor
    public  :: KillFillMatrix                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjFillMatrix 
    
    !Data
    public  ::  T_Station
    
    !Interfaces----------------------------------------------------------------
    interface ConstructFillMatrix
        module procedure ConstructFillMatrix2D
        module procedure ConstructFillMatrix3D
        module procedure ConstructFillMatrix2DVectorial
        module procedure ConstructFillMatrix3DVectorial        
    end interface ConstructFillMatrix
    
    !interface ModifyFillMatrix
    !    module procedure ModifyFillMatrix
    !    module procedure ModifyFillMatrixVectorial       
    !end interface ModifyFillMatrix    

    interface  UngetFillMatrix
        module procedure UngetFillMatrix2D
        module procedure UngetFillMatrix3D
    end interface UngetFillMatrix

    interface  GetVectorialField
        module procedure GetVectorialField2D
        module procedure GetVectorialField3D
    end interface GetVectorialField
    
    interface  GetAngleField
        module procedure GetAngleField2D
        module procedure GetAngleField3D
    end interface GetAngleField    
    
    interface  GetDefaultValue
        module procedure GetDefaultValueScalar
        module procedure GetDefaultValueVectorial
    end interface  GetDefaultValue
    
    !Parameter-----------------------------------------------------------------
    integer, parameter                              :: Time_4D_         = 1
    integer, parameter                              :: Generic_4D_      = 2

    integer, parameter                              :: Dim2D            = 2
    integer, parameter                              :: Dim3D            = 3
    
    !Initialization Methods
    integer, parameter                              :: Constant         = 1
    integer, parameter                              :: Layers           = 2
    integer, parameter                              :: Boxes            = 3
    integer, parameter                              :: ASCIIFile        = 4
    integer, parameter                              :: Profile          = 5
    integer, parameter                              :: ReadTimeSerie    = 6
    integer, parameter                              :: ReadHDF          = 7
    integer, parameter                              :: AnalyticProfile  = 8
    integer, parameter                              :: ProfileTimeSerie = 9
    integer, parameter                              :: Sponge           = 10
    integer, parameter                              :: MultiTimeserie   = 11
    integer, parameter                              :: AnalyticWave     = 12
    
    !Data Processing Type (used with MultiTimeserie)    
    integer, parameter                              :: Interpolate      = 1 !Will interpolate the value between two times
    integer, parameter                              :: Accumulate       = 2 !Will distribute the value in time
    integer, parameter                              :: NoProcessing     = 3 !Will use original value (exact or most right one)
    
    !Variable from file
    integer, parameter                              :: None             = 1

    !Type of analytic profiles 
    integer, parameter                              :: Linear           = 1
    integer, parameter                              :: Exponential      = 2

    !type of values
    integer, parameter                              :: InterpolatedValues = 1
    integer, parameter                              :: AccumulatedValues  = 2
    integer, parameter                              :: OriginalValues     = 3
    
    !type of values 
    integer, parameter                              :: sponge_exp_                = 1
    integer, parameter                              :: sponge_linear_             = 2
    integer, parameter                              :: sponge_wave_stress_dump_   = 3
    
    !wave types
    integer, parameter                              :: SineWaveSeaLevel_          = 1
    integer, parameter                              :: CnoidalWaveSeaLevel_       = 2
    integer, parameter                              :: SolitartyWaveSeaLevel_     = 3
    integer, parameter                              :: SineWaveVelX_              = 4
    integer, parameter                              :: CnoidalWaveVelX_           = 5
    integer, parameter                              :: SolitartyWaveVelX_         = 6
    integer, parameter                              :: SineWaveVelY_              = 7
    integer, parameter                              :: CnoidalWaveVelY_           = 8
    integer, parameter                              :: SolitartyWaveVelY_         = 9

    
    !wave cell types
    integer, parameter                              :: EnteringWaveCell_  = 1
    integer, parameter                              :: LeavingWaveCell_   = 2
    integer, parameter                              :: InteriorWaveCell_  = 3
    integer, parameter                              :: ExteriorWaveCell_  = 4    
    
    
    !Parameter-----------------------------------------------------------------
    character(LEN = StringLength), parameter :: BeginProfile      = '<BeginProfile>'
    character(LEN = StringLength), parameter :: EndProfile        = '<EndProfile>'
    character(LEN = StringLength), parameter :: BeginLayers       = '<BeginLayers>'
    character(LEN = StringLength), parameter :: EndLayers         = '<EndLayers>'
    character(LEN = StringLength), parameter :: BeginDepth        = '<BeginDepth>'
    character(LEN = StringLength), parameter :: EndDepth          = '<EndDepth>'
    character(LEN = StringLength), parameter :: BeginTimes        = '<BeginTime>'
    character(LEN = StringLength), parameter :: EndTimes          = '<EndTime>'
    character(LEN = StringLength), parameter :: BeginProfileValues= '<BeginProfileValues>'
    character(LEN = StringLength), parameter :: EndProfileValues  = '<EndProfileValues>'
    character(LEN = StringLength), parameter :: BeginMTSBlock     = '<BeginStation>'
    character(LEN = StringLength), parameter :: EndMTSBlock       = '<EndStation>'
    
    character(LEN = StringLength), parameter :: Begin_files       = '<<BeginFiles>>'
    character(LEN = StringLength), parameter :: End_files         = '<<EndFiles>>'
    character(LEN = StringLength), parameter :: Begin_files_2     = '<<<BeginFiles>>>'
    character(LEN = StringLength), parameter :: End_files_2       = '<<<EndFiles>>>'
    

    !Types---------------------------------------------------------------------

    type T_Layers
        real, dimension(:), pointer                 :: Values   => null()
    end type T_Layers

    type T_Boxes
        character(PathLength)                       :: FileName     = null_str 
        integer                                     :: ObjBoxDif    = null_int 
        real, dimension(:), pointer                 :: Values       => null()
    end type T_Boxes

    type T_ASCIIFile
        character(PathLength)                       :: FileName     = null_str 
        integer                                     :: GridDataID   = null_int 
        type (T_ASCIIFile), pointer                 :: Next             => null()
        type (T_ASCIIFile), pointer                 :: Prev             => null()         
    end type T_ASCIIFile

    type T_Sponge
        real                                        :: OutValue     = null_real
        integer                                     :: Cells        = null_int 
        logical                                     :: Growing      = .false.  
        integer                                     :: Evolution    = null_int 
        !1 - South; 2 - North; 3 - West; 4 - East        
        logical, dimension(1:4)                     :: OpenBordersON = .true.
    end type T_Sponge
    
    type T_EnteringCell
        integer, dimension(:),   pointer            :: i             => null()
        integer, dimension(:),   pointer            :: j             => null()
        real,    dimension(:,:), pointer            :: dx            => null()
        real,    dimension(:,:), pointer            :: dy            => null()    
        real,    dimension(:),   pointer            :: TimeLag       => null()
        integer                                     :: nCells        =  null_int
        integer                                     :: iStart        = null_int
        integer                                     :: jStart        = null_int
        integer                                     :: n             = null_int
    end type T_EnteringCell

    type T_AnalyticWave
        logical                                     :: ON               = .false. 
        real                                        :: Amplitude        = null_real
        real                                        :: Direction        = null_real
        real                                        :: Period           = null_real
        real                                        :: AverageValue     = null_real
        real                                        :: DepthValue       = null_real
        real                                        :: CoefValue        = null_real
        real                                        :: Dif              = null_real
        logical                                     :: SlowStartON      = .false. 
        real                                        :: SlowStartPeriod  = null_real
       ! WaveType = 1 (Sine), WaveType = 2 (Cnoidal), WaveType = 3 (solitary)
        integer                                     :: WaveType      = null_int
        real,    dimension(:,:), pointer            :: X2D           => null()
        real,    dimension(:,:), pointer            :: Celerity      => null()
        real,    dimension(:,:), pointer            :: AmpAux        => null()
!        EnteringWaveCell_  = 1
!        LeavingWaveCell_   = 2
!        InteriorWaveCell_  = 3
!        ExteriorWaveCell_  = 4
        integer, dimension(:,:), pointer            :: CellType     => null()
        integer, dimension(:,:), pointer            :: TlagMissing  => null()        
        type (T_EnteringCell)                       :: EnteringCell 
    end type T_AnalyticWave


    type T_TimeSerie
        character(PathLength)                       :: FileName     = null_str 
        integer                                     :: ObjTimeSerie = 0
        integer                                     :: Column       = null_int 
        type (T_Time)                               :: NextTime, &
                                                       PreviousTime
        integer                                     :: NextInstant      = 0,     &
                                                       PreviousInstant  = 0
        real                                        :: PreviousValue    = 0.,    & 
                                                       NextValue        = 0.,    &
                                                       CurrentValue     = 0.
        integer                                     :: NumberOfInstants = 0
        
        logical                                     :: RemainsConstant  = .false.
        
        real                                        :: PredictedDT          = -null_real
        real                                        :: DTForNextEvent       = -null_real
        real                                        :: DTForNextDataset     = -null_real
        type(T_Time)                                :: NextEventStart
        type(T_Time)                                :: NextEventEnd
        type(T_Time)                                :: TimeOfNextDataset
        integer                                     :: InstantOfNextDataset  = null_int
        real                                        :: NextValueForDTPred    = -null_real
        
        type (T_TimeSerie), pointer                 :: Next             => null()
        type (T_TimeSerie), pointer                 :: Prev             => null()        
        
    end type T_TimeSerie

    type T_ProfileTimeSerie
        character(PathLength)                       :: FileName = null_str 
        type (T_Time)                               :: NextTime, PreviousTime
        integer                                     :: NextInstant      = null_int, & 
                                                       PreviousInstant  = null_int    
        real,           dimension(:,:,:), pointer   :: PreviousField3D  => null(), &
                                                       NextField3D      => null()
        real,           dimension(:,:  ), pointer   :: Values           => null(), &
                                                       Depths           => null()
        type(T_Time),   dimension(:    ), pointer   :: TimeInstants     => null()
        integer                                     :: NumberOfInstants = null_int, & 
                                                       nValues          = null_int, & 
                                                       nDepths          = null_int, &      
                                                       FirstInstant     = null_int, & 
                                                       LastInstant      = null_int    
        logical                                     :: CyclicTimeON     = .false.
    end type T_ProfileTimeSerie

    type T_Station   
        character(PathLength)                       :: FileName         = null_str 
        integer                                     :: ObjTimeSerie     = 0        
        integer                                     :: Column           = null_int 
        integer                                     :: FillID           = null_int 
        logical                                     :: RemainConstant   = .false.        
        logical                                     :: ValueIsDefined   = .false.
        real                                        :: NewValue         = -null_real                      
        integer                                     :: NumberOfInstants = 0 
        type(T_Time)                                :: NextEventStart
        type(T_Time)                                :: NextEventEnd
        type(T_Time)                                :: NextTime
        type(T_Time)                                :: PreviousTime
        integer                                     :: NextInstant      = 0
        integer                                     :: PreviousInstant  = 0
        real                                        :: PredictedDT      = -null_real
        real                                        :: DTForNextEvent   = 0.
        real                                        :: DTForNextDataset = -null_real
        real                                        :: NextValue        = 0.
        real                                        :: PreviousValue    = 0.
        real                                        :: NextValueForDTPred = 0.        
    end type T_Station

    type T_MultiTimeSerie
        integer, dimension(:,:), allocatable        :: FillGrid2D
        integer, dimension(:,:,:), allocatable      :: FillGrid3D
        integer                                     :: DataProcessing   = null_int 
        type(T_Station), dimension(:), allocatable  :: StationsList
        integer                                     :: NumberOfSources  = 0
    end type T_MultiTimeserie

    !Generic 4D
    type T_Generic4D
        logical                                     :: ON                   = .false.   
        logical                                     :: ReadFromTimeSerie    = .false.   
        integer                                     :: ObjTimeSerie         = null_int  
        integer                                     :: TimeSerieColumn      = null_int  
        real                                        :: CurrentValue         = null_real 

    end type 

    type T_Field4D
        character(PathLength)                       :: FileName             = null_str, & 
                                                       VGroupPath           = null_str, & 
                                                       FieldName            = null_str 
        real                                        :: MultiplyingFactor    = null_real 
        logical                                     :: HasMultiplyingFactor = .false.
        real                                        :: AddingFactor         = null_real 
        logical                                     :: HasAddingFactor      = .false.
        type (T_Time)                               :: NextTime,  PreviousTime
        type (T_Time)                               :: StartTime,  EndTime        
        real                                        :: Next4DValue          = FillValueReal
        real                                        :: Previous4DValue      = FillValueReal
        integer                                     :: NextInstant          = null_int, & 
                                                       PreviousInstant      = null_int    
        real, dimension(:,:  ), pointer             :: PreviousField2D      => null(), &
                                                       NextField2D          => null(), &
                                                       Array2D              => null()
        real, dimension(:,:,:), pointer             :: PreviousField3D      => null(), &
                                                       NextField3D          => null(), &
                                                       ReadField3D          => null(), &
                                                       Array3D              => null()
        
        
        real                                        :: PredictedDT          = -null_real
        real                                        :: DTForNextEvent       = -null_real
        real                                        :: DTForNextDataset     = -null_real
        type(T_Time)                                :: NextEventStart
        type(T_Time)                                :: NextEventEnd
        type(T_Time)                                :: TimeOfNextDataset
        integer                                     :: InstantOfNextDataset  = null_int
        real                                        :: NextValueForDTPred    = -null_real        
        
        integer                                     :: ObjHDF5              =  0
        logical                                     :: RemainsConstant      = .false.
        integer                                     :: NumberOfInstants     = null_int 
        logical                                     :: CyclicTimeON         = .false.
        logical                                     :: From2Dto3D           = .false.
        logical                                     :: From3Dto2D           = .false.        
        type(T_Generic4D)                           :: Generic4D
        !logical                                     :: ArgumentFileName     = .false. 
        integer                                     :: ObjField4D           = 0
        logical                                     :: Field4D              = .false.
        logical                                     :: HarmonicsON          = .false.
        real                                        :: HarmonicsDT          = null_real
        logical                                     :: SpatialInterpolON    = .false.
        logical                                     :: InterpolOnlyVertically = .false.        
        logical                                     :: GenericYear          = .false.
        integer                                     :: Ncells
        real,    dimension(:), pointer              :: X                    => null()
        real,    dimension(:), pointer              :: Y                    => null()        
        real,    dimension(:), pointer              :: Z                    => null()
        real,    dimension(:), pointer              :: Prop                 => null()
        logical, dimension(:), pointer              :: NoData               => null()  
        logical                                     :: Extrapolate          = .false.       
        integer                                     :: ExtrapolateMethod    = null_int
        character(len=PathLength), dimension(:), pointer :: FileNameList    => null()
        type (T_Field4D), pointer                   :: Next                 => null()
        type (T_Field4D), pointer                   :: Prev                 => null()         
    end type T_Field4D


    private :: T_FillMatrix
    type       T_FillMatrix
        integer                                     :: InstanceID
        type (T_Size2D)                             :: Size2D, WorkSize2D
        type (T_Size3D)                             :: Size3D, WorkSize3D
        type (T_PropertyID)                         :: PropertyID
        integer                                     :: Dim                      = null_int 
        integer                                     :: TypeZUV                  = null_int 
        integer                                     :: TimeEvolution            = null_int 
        integer                                     :: SpaceEvolution           = null_int 
        integer                                     :: InitializationMethod     = null_int 
        integer                                     :: InitializationDefault    = null_int 
        logical                                     :: RemainsConstant          = .false.
        
        logical                                     :: VectorialProp            = .false.  ! w/ x and y fields (orig & rotated)
        logical                                     :: RotateAngleToGrid        = .false.  !scalar w/ original and rotated field

!        logical                                     :: AccumulatedValue     = .false.
!        logical                                     :: NoInterpol           = .false.
                
!        integer                                     :: ValuesType           
        logical                                     :: InterpolateValues        = .false. 
        logical                                     :: AccumulateValues         = .false. 
        logical                                     :: UseOriginalValues        = .false. 
        logical                                     :: PreviousInstantValues    = .false. 
        logical                                     :: IgnoreNoDataPoint        = .false. 
        integer                                     :: PredictDTMethod !1 for old method, 2 for new method (for rain, mainly)
        real                                        :: NoDataValue              = null_real 
        
        logical                                     :: Backtracking             = .false.
        
        character(len=StringLength)                 :: OverrideValueKeyword     = null_str
        logical                                     :: OverrideValueKeywordON   = .false.
         
        real                                        :: MinForDTDecrease         = AllmostZero
        !real                                        :: DefaultValue             = null_real 
        real, dimension(3)                          :: DefaultValue             = null_real 
        real                                        :: PredictedDT              = -null_real
        real                                        :: DTForNextEvent           = -null_real
        real                                        :: DTForNextDataset     = -null_real
        !type(T_Time)                                :: NextEventStart
        !type(T_Time)                                :: NextEventEnd
        !type(T_Time)                                :: TimeOfNextDataset
        !integer                                     :: InstantOfNextDataset
        logical                                     :: ValueIsUsedForDTPrediction = .false.
        real                                        :: NextValueForDTPred
        real,    dimension(:, :   ), pointer        :: Matrix2D         => null()
        real,    dimension(:, :   ), pointer        :: Matrix2DX        => null()   !input vectorial field x (meridian value)
        real,    dimension(:, :   ), pointer        :: Matrix2DY        => null()   !input vectorial field y (parallel value)
        real,    dimension(:, :   ), pointer        :: Matrix2DU        => null()   !model grid vectorial field U (cell ref)
        real,    dimension(:, :   ), pointer        :: Matrix2DV        => null()   !model grid vectorial field V (cell ref)
        real,    dimension(:, :   ), pointer        :: Matrix2DFieldAngle => null() !input angle field in nautic or current ref
        real,    dimension(:, :   ), pointer        :: Matrix2DCellAngle  => null() !cell angle field in cell ref
        real,    dimension(:, :, :), pointer        :: Matrix3D         => null()
        real,    dimension(:, :, :), pointer        :: Matrix3DX        => null()
        real,    dimension(:, :, :), pointer        :: Matrix3DY        => null()
        real,    dimension(:, :, :), pointer        :: Matrix3DZ        => null()
        real,    dimension(:, :, :), pointer        :: Matrix3DU        => null()
        real,    dimension(:, :, :), pointer        :: Matrix3DV        => null()
        real,    dimension(:, :, :), pointer        :: Matrix3DW        => null()
        real,    dimension(:, :, :), pointer        :: Matrix3DFieldAngle => null()
        real,    dimension(:, :, :), pointer        :: Matrix3DCellAngle  => null()
        integer, dimension(:, :   ), pointer        :: PointsToFill2D   => null()        
        integer, dimension(:, :, :), pointer        :: PointsToFill3D   => null()
        
        logical                                     :: UseZ             = .false.   !In case of 3D,use Z component?
        
        type(T_Time)                                :: BeginTime, EndTime
        
        logical                                     :: ArgumentFileName     = .false. 
        character(len=PathLength), dimension(2)     :: FileNameHDF          = null_str
        type(T_ASCIIFile), pointer                  :: FirstASCIIFile      
        type(T_TimeSerie), pointer                  :: FirstTimeSerie      
        type(T_Field4D), pointer                    :: FirstHDF
        integer                                     :: nASCIIFiles      = null_int
        integer                                     :: nTimeSeries      = null_int
        integer                                     :: nHDFs            = null_int
        
        logical                                     :: CheckDates = .true.
        character(len=PathLength)                   :: SpongeFILE_DT    = null_str

        logical                                     :: CheckHDF5_File   = .false. 

        !Initialization Methods
        type (T_Layers   )                          :: Layers
        type (T_Boxes    )                          :: Boxes
        type (T_TimeSerie)                          :: TimeSerie
        type (T_ASCIIFile)                          :: ASCIIFile
        type (T_Sponge   )                          :: Sponge
        type (T_AnalyticWave)                       :: AnalyticWave
        type (T_Field4D  )                          :: HDF
        type (T_ProfileTimeSerie)                   :: ProfileTimeSerie
        type (T_MultiTimeSerie)                     :: MultiTimeSerie
        integer                                     :: ObjEnterData         = 0 
        integer                                     :: ObjTime              = 0
        integer                                     :: ObjHorizontalGrid    = 0    
        integer                                     :: ObjGeometry          = 0    

        type(T_FillMatrix), pointer                 :: Next => null()
    end type  T_FillMatrix

    !Global Module Variables
    type (T_FillMatrix), pointer                    :: FirstObjFillMatrix   => null()
    type (T_FillMatrix), pointer                    :: Me                   => null()

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructFillMatrix2D(PropertyID, EnterDataID, TimeID,                   &
                                     HorizontalGridID, ExtractType, PointsToFill2D,     &
                                     Matrix2D, TypeZUV, Matrix2DInputRef, FileNameHDF,  &    
                                     ObjFillMatrix, OverrideValueKeyword, ClientID,     &
                                     PredictDTMethod, MinForDTDecrease,                 &
                                     ValueIsUsedForDTPrediction, CheckDates,            &
                                     RotateAngleToGrid, SpongeFILE_DT, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: EnterDataID
        integer                                         :: TimeID
        integer                                         :: HorizontalGridID
        integer                                         :: ExtractType
        integer, dimension(:, :), pointer               :: PointsToFill2D
        real, dimension(:, :), pointer                  :: Matrix2D
        real, dimension(:, :), pointer, optional        :: Matrix2DInputRef      !original field (e.g. angle)
        integer                                         :: TypeZUV
        type (T_PropertyID)                             :: PropertyID
        integer, optional, intent(IN)                   :: ClientID
        character(*), optional, intent(IN )             :: FileNameHDF, OverrideValueKeyword
        integer,      optional, intent(INOUT)           :: ObjFillMatrix
        logical,      optional, intent(IN)              :: CheckDates
        integer,      optional, intent(OUT)             :: STAT     
        integer,      optional, intent(IN )             :: PredictDTMethod 
        real,         optional, intent(IN )             :: MinForDTDecrease  
        logical,      optional, intent(IN )             :: ValueIsUsedForDTPrediction
        logical,      optional, intent(IN )             :: RotateAngleToGrid
        character(*), optional, intent(IN )             :: SpongeFILE_DT

        !Local-------------------------------------------------------------------
        integer                                         :: ready_, STAT_, STAT_CALL, nUsers, ObjFillMatrix_
        integer                                         :: PredictDTMethod_, Referential
        type (T_PropertyID), pointer                    :: Prop

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mFillMatrix_)) then
            nullify (FirstObjFillMatrix)
            call RegisterModule (mFillMatrix_) 
        endif
        
        if (present(ObjFillMatrix)) then
            ObjFillMatrix_ = ObjFillMatrix
        else
            ObjFillMatrix_ = PropertyID%ObjFillMatrix
        endif

        call Ready(ObjFillMatrix_, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            if (present(CheckDates)) then
                Me%CheckDates = CheckDates
            endif
            
            if (present(SpongeFILE_DT)) then
                Me%SpongeFILE_DT = trim(SpongeFILE_DT)
            else
                Me%SpongeFILE_DT = null_str
            endif
            
!~             if (Check_Vectorial_Property(PropertyID%IDNumber)) then
            if (PropertyID%IsVectorial) then
                write(*,*) 'Constructing vectorial property but expected scalar'
                stop 'ConstructFillMatrix2D - ModuleFillMatrix - ERR00'
            endif
            
            if (present(PredictDTMethod)) then
                PredictDTMethod_ = PredictDTMethod
            else
                PredictDTMethod_ = 1
            endif
            
            if (present(MinForDTDecrease)) then
                Me%MinForDTDecrease = MinForDTDecrease
            endif
            
            if (present(ValueIsUsedForDTPrediction)) then
                Me%ValueIsUsedForDTPrediction = ValueIsUsedForDTPrediction
            else
                Me%ValueIsUsedForDTPrediction = .false.
            endif
            
            Me%ObjEnterData      = AssociateInstance (mENTERDATA_,      EnterDataID     )
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)

            ! JPW 2012-01-28: Use HorizontalGridSize to set size of the matrix.
            ! Using ubound may cause inconsistency with padded matrices (making IUB larger than the actual grid bound).
            ! PointsToFill2D should have the same dimensions as HorizontalGrid anyway.
            ! Matrices that are used for the Thomas algorithm are padded if _USE_PAGELOCKED (CUDA only) or _PAD_MATRICES (Fortran) is defined
            call GetHorizontalGridSize(HorizontalGridID, Me%Size2D, Me%WorkSize2D, STAT = STAT_)
!            Me%Size2D%ILB       = lbound(PointsToFill2D, dim = 1)
!            Me%Size2D%IUB       = ubound(PointsToFill2D, dim = 1)
!            Me%Size2D%JLB       = lbound(PointsToFill2D, dim = 2)
!            Me%Size2D%JUB       = ubound(PointsToFill2D, dim = 2)
!            Me%WorkSize2D%ILB   = Me%Size2D%ILB + 1
!            Me%WorkSize2D%IUB   = Me%Size2D%IUB - 1
!            Me%WorkSize2D%JLB   = Me%Size2D%JLB + 1
!            Me%WorkSize2D%JUB   = Me%Size2D%JUB - 1


            Me%Size3D       = T_Size3D(null_int, null_int, null_int, null_int, null_int, null_int)
            Me%Dim          = Dim2D
            Me%TypeZUV      = TypeZUV
            Me%PropertyID   = PropertyID
            
            if (present(RotateAngleToGrid)) then
                Me%RotateAngleToGrid = RotateAngleToGrid
            else
                if (Me%PropertyID%IsAngle) then
                    Me%RotateAngleToGrid = .true.                    
                endif 
            endif            
            
            Me%Matrix2D       => Matrix2D
            Me%PointsToFill2D => PointsToFill2D
            
            !!get the orginal field. will be given to user to output
!~             if (Check_Angle_Property(Me%PropertyID%IDNumber)) then
            if (Me%RotateAngleToGrid) then
                if (.not. present(Matrix2DInputRef)) then
                    write(*,*) 'Constructing angle property but not given original field'
                    stop 'ConstructFillMatrix2D - ModuleFillMatrix - ERR10'                    
                endif
                
                !point to input angle and cell angle (rotated to cell ref)
                Me%Matrix2DFieldAngle => Matrix2DInputRef
                Me%Matrix2DCellAngle  => Matrix2D
                where (PointsToFill2D == WaterPoint) Me%Matrix2DFieldAngle = null_real
                where (PointsToFill2D == WaterPoint) Me%Matrix2DCellAngle = null_real
            endif            
            
            where (PointsToFill2D == WaterPoint) Me%Matrix2D      = null_real
            
            if (Me%TypeZUV == TypeU_) then
                Me%Size2D%JUB       = Me%Size2D%JUB + 1
                Me%WorkSize2D%JUB   = Me%WorkSize2D%JUB + 1
            endif

            if (Me%TypeZUV == TypeV_) then
                Me%Size2D%IUB       = Me%Size2D%IUB + 1
                Me%WorkSize2D%IUB   = Me%WorkSize2D%IUB + 1
            endif
            
            if (present(FileNameHDF)) then
            
                Me%ArgumentFileName = .true.
                Me%FileNameHDF(1)   = trim(FileNameHDF)
                
            else
            
                Me%ArgumentFileName = .false.
            
            endif
            
            if(present(OverrideValueKeyword))then
                Me%OverrideValueKeyword   = trim(adjustl(OverrideValueKeyword))
                Me%OverrideValueKeywordON = .true.
            else
                Me%OverrideValueKeywordON = .false.
            end if
            
            if (present(ClientID)) then
                call ReadOptions (ExtractType,                          &
                                  PointsToFill2D = PointsToFill2D,      &
                                  PredictDTMethod = PredictDTMethod_,   &
                                  ClientID = ClientID)
            else
                call ReadOptions (ExtractType,                          &
                                  PointsToFill2D = PointsToFill2D,      &
                                  PredictDTMethod = PredictDTMethod_)
            endif
            
            
            !Is this a angle property? convert to cell referential angle
            if (Me%RotateAngleToGrid) then

                !angle referential
                Prop => Me%PropertyID
                Referential = Get_Angle_Referential(Prop)                                    
                
                !!Need to rotate input field            
                call RotateAngleFieldToGrid(HorizontalGridID      = Me%ObjHorizontalGrid,                 &
                                                AngleIn           = Me%Matrix2DFieldAngle,                &
                                                InReferential     = Referential,                          &
                                                AngleOut          = Me%Matrix2DCellAngle,                 &
                                                WaterPoints2D     = PointsToFill2D,                       &
                                                Rotate            = .true.,                               &
                                                STAT              = STAT_CALL)                                       
                
            endif
            
            nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
            if (nUsers == 0) stop 'ConstructFillMatrix2D - ModuleFillMatrix - ERR20' 
            
            if(Me%TimeEvolution == None)then
                PropertyID%SolutionFromFile  = .false.
            else 
                PropertyID%SolutionFromFile  = .true.
            end if

            !Returns ID
            ObjFillMatrix_ = Me%InstanceID
            
            
            if (present(ObjFillMatrix)) then
                ObjFillMatrix            = ObjFillMatrix_ 
            else
                PropertyID%ObjFillMatrix = ObjFillMatrix_
            endif       
            
            nullify(Me%Matrix2D      )
            nullify(Me%Matrix2DFieldAngle )
            nullify(Me%Matrix2DCellAngle )
            nullify(Me%PointsToFill2D)
                 

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ConstructFillMatrix2D - ModuleFillMatrix - ERR02'  

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructFillMatrix2D

    !----------------------------------------------------------------------
 
    subroutine ConstructFillMatrix2DVectorial(PropertyID, EnterDataID, TimeID,          &
                                     HorizontalGridID, ExtractType, PointsToFill2D,     &
                                     Matrix2DU, Matrix2DV, Matrix2DX, Matrix2DY,        &
                                     TypeZUV, FileNameHDF, ObjFillMatrix,     &
                                     OverrideValueKeyword, ClientID, PredictDTMethod,   &
                                     MinForDTDecrease, ValueIsUsedForDTPrediction, CHeckDates, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: EnterDataID
        integer                                         :: TimeID
        integer                                         :: HorizontalGridID
        integer                                         :: ExtractType
        integer, dimension(:, :), pointer               :: PointsToFill2D
        real, dimension(:, :), pointer                  :: Matrix2DU
        real, dimension(:, :), pointer                  :: Matrix2DV
        real, dimension(:, :), pointer                  :: Matrix2DX
        real, dimension(:, :), pointer                  :: Matrix2DY        
        integer                                         :: TypeZUV
        type (T_PropertyID)                             :: PropertyID
        integer, optional, intent(IN)                   :: ClientID
        character(*), dimension(2), optional, intent(IN ) :: FileNameHDF
        character(*), optional, intent(IN )             :: OverrideValueKeyword
        integer,      optional, intent(INOUT)           :: ObjFillMatrix
        logical,      optional, intent(IN)              :: CheckDates
        integer,      optional, intent(OUT)             :: STAT     
        integer,      optional, intent(IN )             :: PredictDTMethod 
        real,         optional, intent(IN )             :: MinForDTDecrease  
        logical,      optional, intent(IN )             :: ValueIsUsedForDTPrediction

        !Local-------------------------------------------------------------------
        integer                                         :: ready_, STAT_, STAT_CALL, nUsers, ObjFillMatrix_
        integer                                         :: PredictDTMethod_
 
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mFillMatrix_)) then
            nullify (FirstObjFillMatrix)
            call RegisterModule (mFillMatrix_) 
        endif
        
        if (present(ObjFillMatrix)) then
            ObjFillMatrix_ = ObjFillMatrix
        else
            ObjFillMatrix_ = PropertyID%ObjFillMatrix
        endif

        call Ready(ObjFillMatrix_, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            if (present(CheckDates)) then
                Me%CheckDates = CheckDates
            endif
            
            Me%VectorialProp = .true.
            
!~             if (.not. Check_Vectorial_Property(PropertyID%IDNumber)) then
            if (.not. PropertyID%IsVectorial) then
                write(*,*) 'Constructing scalar property but expected vectorial'
                stop 'ConstructFillMatrix2DVectorial - ModuleFillMatrix - ERR00'
            endif            
            
            if (present(PredictDTMethod)) then
                PredictDTMethod_ = PredictDTMethod
            else
                PredictDTMethod_ = 1
            endif
            
            if (present(MinForDTDecrease)) then
                Me%MinForDTDecrease = MinForDTDecrease
            endif
            
            if (present(ValueIsUsedForDTPrediction)) then
                Me%ValueIsUsedForDTPrediction = ValueIsUsedForDTPrediction
            else
                Me%ValueIsUsedForDTPrediction = .false.
            endif
            
            Me%ObjEnterData      = AssociateInstance (mENTERDATA_,      EnterDataID     )
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)

            ! JPW 2012-01-28: Use HorizontalGridSize to set size of the matrix.
            ! Using ubound may cause inconsistency with padded matrices (making IUB larger than the actual grid bound).
            ! PointsToFill2D should have the same dimensions as HorizontalGrid anyway.
            ! Matrices that are used for the Thomas algorithm are padded if _USE_PAGELOCKED (CUDA only) or _PAD_MATRICES (Fortran) is defined
            call GetHorizontalGridSize(HorizontalGridID, Me%Size2D, Me%WorkSize2D, STAT = STAT_)
!            Me%Size2D%ILB       = lbound(PointsToFill2D, dim = 1)
!            Me%Size2D%IUB       = ubound(PointsToFill2D, dim = 1)
!            Me%Size2D%JLB       = lbound(PointsToFill2D, dim = 2)
!            Me%Size2D%JUB       = ubound(PointsToFill2D, dim = 2)
!            Me%WorkSize2D%ILB   = Me%Size2D%ILB + 1
!            Me%WorkSize2D%IUB   = Me%Size2D%IUB - 1
!            Me%WorkSize2D%JLB   = Me%Size2D%JLB + 1
!            Me%WorkSize2D%JUB   = Me%Size2D%JUB - 1


            Me%Size3D       = T_Size3D(null_int, null_int, null_int, null_int, null_int, null_int)
            Me%Dim          = Dim2D
            Me%TypeZUV      = TypeZUV
            Me%PropertyID   = PropertyID
            
            !Point result matrixes to argument (user wants rotated to grid result)
            Me%Matrix2DU      => Matrix2DU
            Me%Matrix2DV      => Matrix2DV
            Me%PointsToFill2D => PointsToFill2D
            
            !Point input matrixes to argument (user inpt field)
            Me%Matrix2DX      => Matrix2DX
            Me%Matrix2DY      => Matrix2DY      
     
            
            !Start matrice
            !where (PointsToFill2D == WaterPoint) Me%Matrix2D = null_real
            where (PointsToFill2D == WaterPoint) Me%Matrix2DU = null_real
            where (PointsToFill2D == WaterPoint) Me%Matrix2DV = null_real
            where (PointsToFill2D == WaterPoint) Me%Matrix2DX = null_real
            where (PointsToFill2D == WaterPoint) Me%Matrix2DY = null_real            

            if (Me%TypeZUV == TypeU_) then
                Me%Size2D%JUB       = Me%Size2D%JUB + 1
                Me%WorkSize2D%JUB   = Me%WorkSize2D%JUB + 1
            endif

            if (Me%TypeZUV == TypeV_) then
                Me%Size2D%IUB       = Me%Size2D%IUB + 1
                Me%WorkSize2D%IUB   = Me%WorkSize2D%IUB + 1
            endif
            
            if (present(FileNameHDF)) then
            
                Me%ArgumentFileName = .true.
                Me%FileNameHDF(1)      = trim(FileNameHDF(1))
                Me%FileNameHDF(2)      = trim(FileNameHDF(2))
                
            else
            
                Me%ArgumentFileName = .false.
            
            endif
            
            if(present(OverrideValueKeyword))then
                Me%OverrideValueKeyword   = trim(adjustl(OverrideValueKeyword))
                Me%OverrideValueKeywordON = .true.
            else
                Me%OverrideValueKeywordON = .false.
            end if
            
            !Here the Me%Matrix2DX and Me%Matrix2DY from user input field are filled
            if (present(ClientID)) then
                call ReadOptions (ExtractType,                          &
                                  PointsToFill2D = PointsToFill2D,      &
                                  PredictDTMethod = PredictDTMethod_,   &
                                  ClientID = ClientID)
            else
                call ReadOptions (ExtractType,                          &
                                  PointsToFill2D = PointsToFill2D,      &
                                  PredictDTMethod = PredictDTMethod_)
            endif
            
            nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
            if (nUsers == 0) stop 'ConstructFillMatrix2D - ModuleFillMatrix - ERR20' 
            
            !!Need to rotate input field (Me%Matrix2DX and Me%Matrix2DY) to grid (Me%Matrix2DU and Me%Matrix2DV)     
            call RotateVectorFieldToGrid(HorizontalGridID  = Me%ObjHorizontalGrid, &
                                            VectorInX         = Me%Matrix2DX,                         &
                                            VectorInY         = Me%Matrix2DY,                         &
                                            VectorOutX        = Me%Matrix2DU,                         &
                                            VectorOutY        = Me%Matrix2DV,                         &   
                                            WaterPoints2D     = PointsToFill2D,                       &
                                            RotateX           = .true.,                               &
                                            RotateY           = .true.,                               &
                                            STAT              = STAT_CALL)                
            

            
            if(Me%TimeEvolution == None)then
                PropertyID%SolutionFromFile  = .false.
            else 
                PropertyID%SolutionFromFile  = .true.
            end if

            !Returns ID
            ObjFillMatrix_ = Me%InstanceID
            
            
            if (present(ObjFillMatrix)) then
                ObjFillMatrix            = ObjFillMatrix_ 
            else
                PropertyID%ObjFillMatrix = ObjFillMatrix_
            endif       
            
            !nullify will remove pointer to Matrix2DU, Matrix2DV and X and Y (Module fields) since they were already modifed
            nullify(Me%Matrix2D       )
            nullify(Me%Matrix2DU      )
            nullify(Me%Matrix2DV      )
            nullify(Me%Matrix2DX      )
            nullify(Me%Matrix2DY      )
            
            nullify(Me%PointsToFill2D)

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ConstructFillMatrix2D - ModuleFillMatrix - ERR02'  

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructFillMatrix2DVectorial

    !----------------------------------------------------------------------

    subroutine ConstructFillMatrix3D(PropertyID, EnterDataID, TimeID, HorizontalGridID, &
                                     GeometryID, ExtractType, PointsToFill3D, Matrix3D, &
                                     TypeZUV, Matrix3DInputRef, FillMatrix, FileNameHDF,&
                                     ObjFillMatrix, OverrideValueKeyword, ClientID,     &
                                     predictDTMethod, MinForDTDecrease,                 &
                                     ValueIsUsedForDTPrediction, CheckDates,            &
                                     RotateAngleToGrid, SpongeFILE_DT, STAT)

        !Arguments---------------------------------------------------------------
        type (T_PropertyID)                             :: PropertyID
        integer                                         :: EnterDataID
        integer                                         :: TimeID
        integer                                         :: HorizontalGridID
        integer                                         :: GeometryID
        integer                                         :: ExtractType
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        real, dimension(:, :, :), pointer               :: Matrix3D
        real, dimension(:, :, :), pointer, optional     :: Matrix3DInputRef    !original field (e.g. angle)
        integer                                         :: TypeZUV
        integer, optional, intent(IN)                   :: ClientID
        real        , optional, intent(IN )             :: FillMatrix
        character(*), optional, intent(IN )             :: FileNameHDF, OverrideValueKeyword
        integer,      optional, intent(INOUT)           :: ObjFillMatrix
        logical,      optional, intent(IN)              :: CheckDates
        integer,      optional, intent(OUT)             :: STAT     
        integer,      optional, intent(IN )             :: PredictDTMethod  
        real,         optional, intent(IN )             :: MinForDTDecrease 
        logical,      optional, intent(IN )             :: ValueIsUsedForDTPrediction 
        logical,      optional, intent(IN )             :: RotateAngleToGrid
        character(*), optional, intent(IN )             :: SpongeFILE_DT
        !Local-------------------------------------------------------------------
        real                                            :: FillMatrix_
        integer                                         :: ready_, STAT_, STAT_CALL, nUsers, ObjFillMatrix_
        integer                                         :: PredictDTMethod_, Referential
        type (T_PropertyID), pointer                    :: Prop
 
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mFillMatrix_)) then
            nullify (FirstObjFillMatrix)
            call RegisterModule (mFillMatrix_) 
        endif

        if (present(ObjFillMatrix)) then
            ObjFillMatrix_ = ObjFillMatrix
        else
            ObjFillMatrix_ = PropertyID%ObjFillMatrix
        endif

        call Ready(ObjFillMatrix_, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            if (present(RotateAngleToGrid)) then
                Me%RotateAngleToGrid = RotateAngleToGrid
            else
                if (Me%PropertyID%IsAngle) then
                    Me%RotateAngleToGrid = .true.                    
                endif 
            endif
            
            if (present(CheckDates)) then
                Me%CheckDates = CheckDates
            endif
            
            if (present(SpongeFILE_DT)) then
                Me%SpongeFILE_DT = trim(SpongeFILE_DT)
            else
                Me%SpongeFILE_DT = null_str
            endif            

!~             if (Check_Vectorial_Property(PropertyID%IDNumber)) then
            if (PropertyID%IsVectorial) then
                write(*,*) 'Constructing vectorial property but expected scalar'
                stop 'ConstructFillMatrix3D - ModuleFillMatrix - ERR00'
            endif               
            
            if (present(PredictDTMethod)) then
                PredictDTMethod_ = PredictDTMethod
            else
                PredictDTMethod_ = 1
            endif

            if (present(MinForDTDecrease)) then
                Me%MinForDTDecrease = MinForDTDecrease
            endif

            if (present(ValueIsUsedForDTPrediction)) then
                Me%ValueIsUsedForDTPrediction = ValueIsUsedForDTPrediction
            else
                Me%ValueIsUsedForDTPrediction = .false.
            endif

            Me%ObjEnterData      = AssociateInstance (mENTERDATA_,      EnterDataID     )
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )

            ! JPW 2012-01-28: Use GeometrySize to set size of the matrix.
            ! Using ubound may cause inconsistency with padded matrices (making IUB larger than the actual grid bound).
            ! PointsToFill3D should have the same dimensions as Geometry anyway.
            ! Matrices that are used for the Thomas algorithm are padded if _USE_PAGELOCKED (CUDA only) or _PAD_MATRICES (Fortran) is defined
            call GetGeometrySize(GeometryID , Me%Size3D, Me%WorkSize3D, STAT_)
!            Me%Size3D%ILB       = lbound(PointsToFill3D, dim = 1)
!            Me%Size3D%IUB       = ubound(PointsToFill3D, dim = 1)
!            Me%Size3D%JLB       = lbound(PointsToFill3D, dim = 2)
!            Me%Size3D%JUB       = ubound(PointsToFill3D, dim = 2)
!            Me%Size3D%KLB       = lbound(PointsToFill3D, dim = 3)
!            Me%Size3D%KUB       = ubound(PointsToFill3D, dim = 3)
!            Me%WorkSize3D%ILB   = Me%Size3D%ILB + 1
!            Me%WorkSize3D%IUB   = Me%Size3D%IUB - 1
!            Me%WorkSize3D%JLB   = Me%Size3D%JLB + 1
!            Me%WorkSize3D%JUB   = Me%Size3D%JUB - 1
!            Me%WorkSize3D%KLB   = Me%Size3D%KLB + 1
!            Me%WorkSize3D%KUB   = Me%Size3D%KUB - 1


            Me%Size2D       = T_Size2D(null_int, null_int, null_int, null_int)
            Me%Dim          = Dim3D
            Me%TypeZUV      = TypeZUV
            Me%PropertyID   = PropertyID

            Me%Matrix3D       => Matrix3D
            Me%PointsToFill3D => PointsToFill3D
            
            !!get the orginal field. will be given to user to output
!~             if (Check_Angle_Property(Me%PropertyID%IDNumber)) then
            if (Me%RotateAngleToGrid) then
                if (.not. present(Matrix3DInputRef)) then
                    write(*,*) 'Constructing angle property but not given original field'
                    stop 'ConstructFillMatrix3D - ModuleFillMatrix - ERR10'                    
                endif
                
                !point to input angle and cell angle (rotated to cell ref)
                Me%Matrix3DFieldAngle => Matrix3DInputRef
                Me%Matrix3DCellAngle  => Matrix3D
                where (PointsToFill3D == WaterPoint) Me%Matrix3DFieldAngle = null_real
                where (PointsToFill3D == WaterPoint) Me%Matrix3DCellAngle  = null_real
            endif 

            if (present(FillMatrix)) then
                FillMatrix_ = FillMatrix
            else
                FillMatrix_ = null_real
            endif

            where (PointsToFill3D == WaterPoint) Me%Matrix3D = FillMatrix_

            if (Me%TypeZUV == TypeU_) then
                Me%Size3D%JUB       = Me%Size3D%JUB + 1
                Me%WorkSize3D%JUB   = Me%WorkSize3D%JUB + 1
            endif

            if (Me%TypeZUV == TypeV_) then
                Me%Size3D%IUB       = Me%Size3D%IUB + 1
                Me%WorkSize3D%IUB   = Me%WorkSize3D%IUB + 1
            endif

            !Specific of the Vertical_Z matrix see ModuleGeometry 
            if (Me%TypeZUV == TypeW_) then
                Me%WorkSize3D%KLB   = Me%WorkSize3D%KLB - 1
            endif


            if (present(FileNameHDF)) then
            
                Me%ArgumentFileName = .true.
                Me%FileNameHDF(1)  = trim(FileNameHDF)
                
            else
            
                Me%ArgumentFileName = .false.
            
            endif
            
            if(present(OverrideValueKeyword))then
                Me%OverrideValueKeyword   = trim(adjustl(OverrideValueKeyword))
                Me%OverrideValueKeywordON = .true.
            else
                Me%OverrideValueKeywordON = .false.
            end if


            if (present(ClientID)) then
                call ReadOptions (ExtractType,                          &
                                  PointsToFill3D = PointsToFill3D,      &
                                  PredictDTMethod = PredictDTMethod_,   &
                                  ClientID = ClientID)
            else
                call ReadOptions (ExtractType,                          &
                                  PointsToFill3D = PointsToFill3D,      &
                                  PredictDTMethod = PredictDTMethod_)
            endif


            !Is this a angle property? convert to cell referential angle
            if (Me%RotateAngleToGrid) then
                    
                !angle referential
                Prop => Me%PropertyID
                Referential = Get_Angle_Referential(Prop)
                
                !!Need to rotate input field         
                call RotateAngleFieldToGrid(HorizontalGridID      = Me%ObjHorizontalGrid,               &
                                                AngleIn           = Me%Matrix3DFieldAngle,                &
                                                InReferential     = Referential,                          &
                                                AngleOut          = Me%Matrix3DCellAngle,                 &
                                                WaterPoints3D     = PointsToFill3D,                       &
                                                Rotate            = .true.,                               &
                                                KLB               = Me%WorkSize3D%KLB,                    &
                                                KUB               = Me%WorkSize3D%KUB,                    &
                                                STAT              = STAT_CALL)                  
              
                
            endif            
                        
            nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
            if (nUsers == 0) stop 'ConstructFillMatrix3D - ModuleFillMatrix - ERR01' 

            if(Me%TimeEvolution == None)then
                PropertyID%SolutionFromFile  = .false.
            else 
                PropertyID%SolutionFromFile  = .true.
            end if

            !Returns ID
            ObjFillMatrix_ = Me%InstanceID
                        
            if (present(ObjFillMatrix)) then
                ObjFillMatrix            = ObjFillMatrix_ 
            else
                PropertyID%ObjFillMatrix = ObjFillMatrix_
            endif            
            
            nullify(Me%Matrix3D           )
            nullify(Me%Matrix3DFieldAngle )
            nullify(Me%Matrix3DCellAngle  )
            nullify(Me%PointsToFill3D)

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ConstructFillMatrix3D - ModuleFillMatrix - ERR02' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructFillMatrix3D

    !--------------------------------------------------------------------------

    subroutine ConstructFillMatrix3DVectorial(PropertyID, EnterDataID, TimeID,          &
                                     HorizontalGridID, GeometryID, ExtractType, PointsToFill3D,     &
                                     Matrix3DU, Matrix3DV, Matrix3DW, Matrix3DX, Matrix3DY,        &
                                     TypeZUV, FillMatrix, FileNameHDF, ObjFillMatrix,     &
                                     OverrideValueKeyword, ClientID, PredictDTMethod,   &
                                     MinForDTDecrease, ValueIsUsedForDTPrediction, CheckDates, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: EnterDataID
        integer                                         :: TimeID
        integer                                         :: GeometryID
        integer                                         :: HorizontalGridID
        integer                                         :: ExtractType
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        real, dimension(:, :, :), pointer               :: Matrix3DU
        real, dimension(:, :, :), pointer               :: Matrix3DV
        real, dimension(:, :, :), pointer               :: Matrix3DX
        real, dimension(:, :, :), pointer               :: Matrix3DY        
        integer                                         :: TypeZUV
        type (T_PropertyID)                             :: PropertyID
        real, dimension(:, :, :), pointer, optional     :: Matrix3DW        
        integer, optional, intent(IN)                   :: ClientID
        real        , optional, intent(IN )             :: FillMatrix
        character(*), dimension(2), optional, intent(IN ) :: FileNameHDF
        character(*), optional, intent(IN )             :: OverrideValueKeyword
        integer,      optional, intent(INOUT)           :: ObjFillMatrix
        logical,      optional, intent(IN)              :: CheckDates
        integer,      optional, intent(OUT)             :: STAT     
        integer,      optional, intent(IN )             :: PredictDTMethod 
        real,         optional, intent(IN )             :: MinForDTDecrease  
        logical,      optional, intent(IN )             :: ValueIsUsedForDTPrediction

        !Local-------------------------------------------------------------------
        real                                            :: FillMatrix_
        integer                                         :: ready_, STAT_, STAT_CALL, nUsers, ObjFillMatrix_
        integer                                         :: PredictDTMethod_
 
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mFillMatrix_)) then
            nullify (FirstObjFillMatrix)
            call RegisterModule (mFillMatrix_) 
        endif
        
        if (present(ObjFillMatrix)) then
            ObjFillMatrix_ = ObjFillMatrix
        else
            ObjFillMatrix_ = PropertyID%ObjFillMatrix
        endif

        call Ready(ObjFillMatrix_, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            if (present(CheckDates)) then
                Me%CheckDates = CheckDates
            endif
            
            Me%VectorialProp = .true.

!~             if (.not. Check_Vectorial_Property(PropertyID%IDNumber)) then
            if (.not. PropertyID%IsVectorial) then
                write(*,*) 'Constructing scalar property but expected vectorial'
                stop 'ConstructFillMatrix3DVectorial - ModuleFillMatrix - ERR00'
            endif                                     
            
            if (present(PredictDTMethod)) then
                PredictDTMethod_ = PredictDTMethod
            else
                PredictDTMethod_ = 1
            endif

            if (present(MinForDTDecrease)) then
                Me%MinForDTDecrease = MinForDTDecrease
            endif

            if (present(ValueIsUsedForDTPrediction)) then
                Me%ValueIsUsedForDTPrediction = ValueIsUsedForDTPrediction
            else
                Me%ValueIsUsedForDTPrediction = .false.
            endif

            Me%ObjEnterData      = AssociateInstance (mENTERDATA_,      EnterDataID     )
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )

            ! JPW 2012-01-28: Use GeometrySize to set size of the matrix.
            ! Using ubound may cause inconsistency with padded matrices (making IUB larger than the actual grid bound).
            ! PointsToFill3D should have the same dimensions as Geometry anyway.
            ! Matrices that are used for the Thomas algorithm are padded if _USE_PAGELOCKED (CUDA only) or _PAD_MATRICES (Fortran) is defined
            call GetGeometrySize(GeometryID , Me%Size3D, Me%WorkSize3D, STAT_)
!            Me%Size3D%ILB       = lbound(PointsToFill3D, dim = 1)
!            Me%Size3D%IUB       = ubound(PointsToFill3D, dim = 1)
!            Me%Size3D%JLB       = lbound(PointsToFill3D, dim = 2)
!            Me%Size3D%JUB       = ubound(PointsToFill3D, dim = 2)
!            Me%Size3D%KLB       = lbound(PointsToFill3D, dim = 3)
!            Me%Size3D%KUB       = ubound(PointsToFill3D, dim = 3)
!            Me%WorkSize3D%ILB   = Me%Size3D%ILB + 1
!            Me%WorkSize3D%IUB   = Me%Size3D%IUB - 1
!            Me%WorkSize3D%JLB   = Me%Size3D%JLB + 1
!            Me%WorkSize3D%JUB   = Me%Size3D%JUB - 1
!            Me%WorkSize3D%KLB   = Me%Size3D%KLB + 1
!            Me%WorkSize3D%KUB   = Me%Size3D%KUB - 1


            Me%Size2D       = T_Size2D(null_int, null_int, null_int, null_int)
            Me%Dim          = Dim3D
            Me%TypeZUV      = TypeZUV
            Me%PropertyID   = PropertyID            
            
            !Point to result matrixes
            Me%Matrix3DU      => Matrix3DU
            Me%Matrix3DV      => Matrix3DV
            if(present(Matrix3DW)) Me%Matrix3DW      => Matrix3DW
            Me%PointsToFill3D => PointsToFill3D
            
            !Point to input matrixes
            Me%Matrix3DX      => Matrix3DX
            Me%Matrix3DY      => Matrix3DY
            !No need for Z that is the same as W (not transformation)

            if (present(FillMatrix)) then
                FillMatrix_ = FillMatrix
            else
                FillMatrix_ = null_real
            endif

            where (PointsToFill3D == WaterPoint) Me%Matrix3DU = FillMatrix_
            where (PointsToFill3D == WaterPoint) Me%Matrix3DV = FillMatrix_
            if(present(Matrix3DW)) where (PointsToFill3D == WaterPoint) Me%Matrix3DW = FillMatrix_
            where (PointsToFill3D == WaterPoint) Me%Matrix3DX = FillMatrix_
            where (PointsToFill3D == WaterPoint) Me%Matrix3DY = FillMatrix_            

            if (Me%TypeZUV == TypeU_) then
                Me%Size3D%JUB       = Me%Size3D%JUB + 1
                Me%WorkSize3D%JUB   = Me%WorkSize3D%JUB + 1
            endif

            if (Me%TypeZUV == TypeV_) then
                Me%Size3D%IUB       = Me%Size3D%IUB + 1
                Me%WorkSize3D%IUB   = Me%WorkSize3D%IUB + 1
            endif

            !Specific of the Vertical_Z matrix see ModuleGeometry 
            if (Me%TypeZUV == TypeW_) then
                Me%WorkSize3D%KLB   = Me%WorkSize3D%KLB - 1
            endif


            if (present(FileNameHDF)) then
            
                Me%ArgumentFileName = .true.
                Me%FileNameHDF(1)      = trim(FileNameHDF(1))
                Me%FileNameHDF(2)      = trim(FileNameHDF(2))
                
            else
            
                Me%ArgumentFileName = .false.
            
            endif
            
            if(present(OverrideValueKeyword))then
                Me%OverrideValueKeyword   = trim(adjustl(OverrideValueKeyword))
                Me%OverrideValueKeywordON = .true.
            else
                Me%OverrideValueKeywordON = .false.
            end if

            !Here matrixes Me%Matrix3DX and Me%Matrix3DY are filled from user info
            if (present(ClientID)) then
                call ReadOptions (ExtractType,                          &
                                  PointsToFill3D = PointsToFill3D,      &
                                  PredictDTMethod = PredictDTMethod_,   &
                                  ClientID = ClientID)
            else
                call ReadOptions (ExtractType,                          &
                                  PointsToFill3D = PointsToFill3D,      &
                                  PredictDTMethod = PredictDTMethod_)
            endif                               
            
            !!Need to rotate input field (Me%Matrix3DX and Me%Matrix3DY) to grid (Me%Matrix3DU and Me%Matrix3DV))            
            call RotateVectorFieldToGrid(HorizontalGridID  = Me%ObjHorizontalGrid, &
                                            VectorInX         = Me%Matrix3DX,                         &
                                            VectorInY         = Me%Matrix3DY,                         &
                                            VectorOutX        = Me%Matrix3DU,                         &
                                            VectorOutY        = Me%Matrix3DV,                         &   
                                            WaterPoints3D     = PointsToFill3D,                       &
                                            RotateX           = .true.,                               &
                                            RotateY           = .true.,                               &
                                            KLB               = Me%WorkSize3D%KLB,                    &
                                            KUB               = Me%WorkSize3D%KUB,                    &            
                                            STAT              = STAT_CALL)                
            

            
            nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
            if (nUsers == 0) stop 'ConstructFillMatrix3D - ModuleFillMatrix - ERR01' 

            if(Me%TimeEvolution == None)then
                PropertyID%SolutionFromFile  = .false.
            else 
                PropertyID%SolutionFromFile  = .true.
            end if

            !Returns ID
            ObjFillMatrix_ = Me%InstanceID
                        
            if (present(ObjFillMatrix)) then
                ObjFillMatrix            = ObjFillMatrix_ 
            else
                PropertyID%ObjFillMatrix = ObjFillMatrix_
            endif            
            
            !nullify will remove pointer to Matrix3DU and Matrix3DV (Module field) since they were alrady modifed
            !when modifying Me%Matrix3DU and Me%Matrix3DV    
            nullify(Me%Matrix3D      )
            nullify(Me%Matrix3DU     )
            nullify(Me%Matrix3DV     )
            nullify(Me%Matrix3DW     )
            nullify(Me%Matrix3DX     )
            nullify(Me%Matrix3DY     )
            
            nullify(Me%PointsToFill3D)          

            STAT_ = SUCCESS_                 


        else cd0
            
            stop 'ConstructFillMatrix2D - ModuleFillMatrix - ERR02'  

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructFillMatrix3DVectorial

    !----------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_FillMatrix), pointer                         :: NewObjFillMatrix
        type (T_FillMatrix), pointer                         :: PreviousObjFillMatrix


        !Allocates new instance
        allocate (NewObjFillMatrix)
        nullify  (NewObjFillMatrix%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjFillMatrix)) then
            FirstObjFillMatrix    => NewObjFillMatrix
            Me                    => NewObjFillMatrix
        else
            PreviousObjFillMatrix => FirstObjFillMatrix
            Me                    => FirstObjFillMatrix%Next
            do while (associated(Me))
                PreviousObjFillMatrix  => Me
                Me                     => Me%Next
            enddo
            Me                          => NewObjFillMatrix
            PreviousObjFillMatrix%Next  => NewObjFillMatrix
        endif

        Me%InstanceID = RegisterNewInstance (mFILLMATRIX_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ReadOptions (ExtractType, PointsToFill2D, PointsToFill3D, ClientID, PredictDTMethod)

        !Arguments-------------------------------------------------------------
        integer                                         :: ExtractType
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D
        integer, optional, intent(IN)                   :: ClientID
        integer, optional, intent(IN)                   :: PredictDTMethod

        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL, ClientID_
        integer                                         :: iflag
        character(len=StringLength)                     :: AuxString
        logical                                         :: AuxBoolean
        integer                                         :: AuxInteger
        logical                                         :: CanContinue
        !---------------------------------------------------------------------
        
        if (present(ClientID)) then
            ClientID_ = ClientID
        else
            ClientID_ = FillValueInt           
        endif

        !Reads Time Evolution
        call GetData(AuxString, Me%ObjEnterData,  iflag,                            &
                     SearchType     = ExtractType,                                  &
                     keyword        = 'FILE_IN_TIME',                               &
                     default        = "None",                                       &
                     ClientModule   = 'ModuleFillMatrix',                           &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR010'

        !write(*,*) trim(adjustl(AuxString))

        select case (trim(adjustl(AuxString)))
            case ("None",       "NONE", "none")
                Me%TimeEvolution    = None
            case ("Hdf",        "HDF",          "hdf")
                Me%TimeEvolution    = ReadHDF
            case ("Timeserie",  "TIMESERIE",    "timeserie",    "TimeSerie")
                Me%TimeEvolution    = ReadTimeSerie
            case ("Profile_Timeserie",  "PROFILE_TIMESERIE",    "profile_timeserie",    "Profile_TimeSerie")
                Me%TimeEvolution    = ProfileTimeSerie
            case ("Multitimeserie", "MULTITIMESERIE", "multitimeserie", "MultiTimeserie")
                Me%TimeEvolution    = MultiTimeserie
            case ("Analytic Wave",  "ANALYTIC WAVE",    "analytic wave")
                Me%TimeEvolution    = AnalyticWave                     
            case default
                write(*,*)'Invalid option for keyword FILE_IN_TIME'
                stop 'ReadOptions - ModuleFillMatrix - ERR020'
        end select

        !write (*,*) Me%TimeEvolution
        
        if (Me%ArgumentFileName) Me%TimeEvolution    = ReadHDF

        !write (*,*) Me%TimeEvolution


        if(Me%TimeEvolution == None)then

            call GetData(AuxString, Me%ObjEnterData,  iflag,                            &
                         SearchType     = ExtractType,                                  &
                         keyword        = 'INITIALIZATION_METHOD',                      &
                         default        = "Constant",                                   &
                         ClientModule   = 'ModuleFillMatrix',                           &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR030'

            !write(*,*) trim(adjustl(AuxString))

            select case (trim(adjustl(AuxString)))
                case ("Constant",   "CONSTANT",   "constant")
                    Me%InitializationMethod = Constant
                case ("Layers",     "LAYERS",     "layers")
                    Me%InitializationMethod = Layers
                case ("Boxes",      "BOXES",      "boxes")
                    Me%InitializationMethod = Boxes
                case ("ASCII_File", "ASCII_FILE", "ascii_file",   "Ascii_file")
                    Me%InitializationMethod = AsciiFile
                case ("Profile",    "PROFILE",    "profile")
                    Me%InitializationMethod = Profile
                case ("Analytic Profile",    "ANALYTIC PROFILE",    "analytic profile")
                    Me%InitializationMethod = AnalyticProfile
                case ("Hdf",        "HDF",          "hdf")
                    Me%InitializationMethod = ReadHDF
                case ("Timeserie",  "TIMESERIE",    "timeserie",    "TimeSerie")
                    Me%InitializationMethod = ReadTimeSerie
                case ("Profile_Timeserie",  "PROFILE_TIMESERIE",    "profile_timeserie",    "Profile_TimeSerie")
                    Me%InitializationMethod = ProfileTimeSerie
                case ("Sponge",  "SPONGE",    "sponge")
                    Me%InitializationMethod = Sponge
                case default
                    write(*,*)'Invalid option for keyword INITIALIZATION_METHOD'
                    stop 'ReadOptions - ModuleFillMatrix - ERR040'
            end select


            call GetData(Me%RemainsConstant, Me%ObjEnterData,  iflag,                   &
                         SearchType     = ExtractType,                                  &
                         keyword        = 'REMAIN_CONSTANT',                            &
                         default        = .false.,                                      &
                         ClientModule   = 'ModuleFillMatrix',                           &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR050'

            call GetData(AuxString, Me%ObjEnterData,  iflag,                            &
                         SearchType     = ExtractType,                                  &
                         keyword        = 'INITIALIZATION_DEFAULT',                     &
                         default        = "Constant",                                   &
                         ClientModule   = 'ModuleFillMatrix',                           &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR060'

            select case (trim(adjustl(AuxString)))
                case ("Constant",   "CONSTANT",   "constant")
                    Me%InitializationDefault = Constant
                case ("Profile_Timeserie",  "PROFILE_TIMESERIE",    "profile_timeserie",    "Profile_TimeSerie")
                    Me%InitializationDefault = ProfileTimeSerie
            end select
        
        end if
        

        ! Check if the simulation goes backward in time or forward in time (default mode)
        call GetBackTracking(Me%ObjTime, Me%BackTracking, STAT = STAT_CALL)                    
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR070' 

        call GetComputeTimeLimits(Me%ObjTime, BeginTime = Me%BeginTime,                 &
                                  EndTime = Me%EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR080'        

        !Gets the default value        
        call GetData(Me%DefaultValue, Me%ObjEnterData,  iflag,                          &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'DEFAULTVALUE',                                   &
                     ClientModule   = 'ModuleFillMatrix',                               &
                     default        = null_real,                                        &
                     STAT           = STAT_CALL)        
        if (STAT_CALL .NE. SUCCESS_) then            
            !normal case, found only one value (STAT_CALL size err and iflag = 1). ignore
            !vectorial 2D case, found only two values (STAT_CALL size err and iflag = 2). ignore
            if (STAT_CALL /= SIZE_ERR_) then
                write(*,*) "Property ", trim(Me%PropertyID%Name)
                write(*,*) "Wrong number of arguments for DEFAULTVALUE keyword"
                stop 'ReadOptions - ModuleFillMatrix - ERR090'                
            endif          
            
            !need to define two or three values, not only one
            if (Me%VectorialProp .and. iflag == 1) then
                write(*,*) 'DEFAULTVALUE for vectorial prop needs to define the x and y'
                write(*,*) '(and optionally z) components in 2/3 values'
                write(*,*) 'separated by space. Property : ', trim(Me%PropertyID%Name)
                stop 'ReadOptions - ModuleFillMatrix - ERR095'
            endif
        endif
        
        if (iflag == 0 .and. .not. Me%ArgumentFileName) then
        
            if(Me%OverrideValueKeywordON .and. Me%InitializationMethod == Constant)then
            
                call GetData(Me%DefaultValue, Me%ObjEnterData,  iflag,                 &
                             SearchType     = ExtractType,                             &
                             keyword        = trim(Me%OverrideValueKeyword),           &
                             ClientModule   = 'ModuleFillMatrix',                      &
                             default        = null_real,                               &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) then
                    
                    if(STAT_CALL == SIZE_ERR_)then 
                        if(.not. Me%VectorialProp .and. iflag .ne. 1)then
                            write(*,*) "Property ", trim(Me%PropertyID%Name)
                            write(*,*) "DEFAULTVALUE must be a scalar"
                            stop 'ReadOptions - ModuleFillMatrix - ERR99'
                        elseif(iflag > 3)then
                            write(*,*) "Property ", trim(Me%PropertyID%Name)
                            write(*,*) "Too many arguments in keyword "//trim(Me%OverrideValueKeyword)
                            stop 'ReadOptions - ModuleFillMatrix - ERR98'
                        endif
                    else
                        stop 'ReadOptions - ModuleFillMatrix - ERR100'
                    endif 
                    
                
                end if
                
                if(iflag == 0)then
                    
                    write(*,*)'Please define override keyword '//trim(Me%OverrideValueKeyword)
                    write(*,*)'to give a default value for property '//trim(Me%PropertyID%Name)
                    stop 'ReadOptions - ModuleFillMatrix - ERR101'
                    
                end if
                
            elseif(Me%OverrideValueKeywordON .and. Me%InitializationMethod == Boxes)then
            
                Me%DefaultValue = null_real
                
            elseif(Me%OverrideValueKeywordON .and. Me%InitializationMethod == AsciiFile)then
            
                Me%DefaultValue = null_real
                
            elseif((Me%OverrideValueKeywordON .and. .not. Me%InitializationMethod == Boxes    ) .or. &
                   (Me%OverrideValueKeywordON .and. .not. Me%InitializationMethod == AsciiFile) .or. &
                   (Me%OverrideValueKeywordON .and. .not. Me%InitializationMethod == Constant ))then
                   
                write(*,*)'Initialization method for property '//trim(Me%PropertyID%Name)
                write(*,*)'can only be CONSTANT, BOXES or ASCII'
                stop 'ReadOptions - ModuleFillMatrix - ERR102'

            
            else

                write(*,*)'Please define default value for property '//trim(Me%PropertyID%Name)
                stop 'ReadOptions - ModuleFillMatrix - ERR103'

            end if
        
        
        end if
       
        !Keywords that define how to "work" values (interpolate, accumulate, nothing)
        call GetData(Me%InterpolateValues,                                              &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'INTERPOLATE_VALUES',                             &
                     default        = .false.,                                          &
                     ClientModule   = 'ModuleFillMatrix',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR104'       
        
        call GetData(Me%AccumulateValues,                                               &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'ACCUMULATE_VALUES',                              &
                     default        = .false.,                                          &
                     ClientModule   = 'ModuleFillMatrix',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR105'     
          
        call GetData(Me%UseOriginalValues,                                              &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'USE_ORIGINAL_VALUES',                            &
                     default        = .false.,                                          &
                     ClientModule   = 'ModuleFillMatrix',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR110'       

        if ((Me%InterpolateValues .AND. Me%AccumulateValues)   .OR. &
            (Me%InterpolateValues .AND. Me%UseOriginalValues)  .OR. &
            (Me%AccumulateValues  .AND. Me%UseOriginalValues)) then            
            write (*,*) 'The keywords INTERPOLATE_VALUES, ACCUMULATE_VALUES and'
            write (*,*) 'USE_ORIGINAL_VALUES are mutually exclusives.'
            write (*,*) 'Only one can be set to true.'
            stop 'ReadOptions - ModuleFillMatrix - ERR120'            
        endif

        !Property ValuesType (Interpolated, accumulate, original value) 
        call GetData(AuxInteger,                                                        &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'VALUES_TYPE',                                    &
                     default        = InterpolatedValues,                               &
                     ClientModule   = 'ModuleFillMatrix',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR120-A'       
        if (iflag .NE. 0) then
            write(*,*) 
            write(*,*) 'ModuleFillMatrix WARNING:'
            write(*,*) 'VALUES_TYPE keyword was removed.'
            write(*,*) 'Use these instead: '
            write(*,*) '   INTERPOLATE_VALUES  => to interpolate values'
            write(*,*) '   ACCUMULATE_VALUES   => to accumulate values'
            write(*,*) '   USE_ORIGINAL_VALUES => to use original values'
            write(*,*)         
            stop 'ReadOptions - ModuleFillMatrix - ERR120-B' 
        endif
        
        !Property not interpolated (e.g. Concentration on rain)
        call GetData(AuxBoolean, Me%ObjEnterData,  iflag,                               &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'NO_INTERPOLATION_OR_ACCUMULATION',               &
                     default        = .false.,                                          &
                     ClientModule   = 'ModuleFillMatrix',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR130'
        if (iflag .NE. 0) then
            write(*,*) 
            write(*,*) 'ModuleFillMatrix WARNING:'
            write(*,*) 'NO_INTERPOLATION_OR_ACCUMULATION keyword is deprecated.'
            write(*,*) 'To use values without interpolation or accumulation use instead: '
            write(*,*) '   "USE_ORIGINAL_VALUES : 1"'
            write(*,*) 
        
            if (AuxBoolean) then            
                Me%InterpolateValues = .false.
                Me%AccumulateValues  = .false.
                Me%UseOriginalValues = .true.
            endif
        endif
        
        if (.NOT. AuxBoolean) then        
            !Accumulitve Property (e.g. Rain from gauges) !this keyword name should be changed because is
            !not describing well what is the computation. However a lot of people use it already...
            call GetData(AuxBoolean, Me%ObjEnterData,  iflag,                               &
                         SearchType     = ExtractType,                                      &
                         keyword        = 'NO_INTERPOLATION',                               &
                         default        = .false.,                                          &
                         ClientModule   = 'ModuleFillMatrix',                               &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR140'
            if (iflag .NE. 0) then
                write(*,*) 
                write(*,*) 'ModuleFillMatrix WARNING:'
                write(*,*) 'NO_INTERPOLATION keyword is deprecated.'
                write(*,*) 'To use interpolated values use instead : '
                write(*,*) '   "INTERPOLATE_VALUES : 1"'
                write(*,*) 'To use accumulated values use instead :'
                write(*,*) '   "ACCUMULATE_VALUES : 1" '
                write(*,*) 
            
                if (AuxBoolean) then                    
                    Me%InterpolateValues = .false.
                    Me%AccumulateValues  = .true.
                    Me%UseOriginalValues = .false.                   
                else       
                    Me%InterpolateValues = .true.
                    Me%AccumulateValues  = .false.
                    Me%UseOriginalValues = .false.                           
                endif            
            endif
        endif  
        
        if ((.NOT. Me%InterpolateValues) .AND. &
            (.NOT. Me%AccumulateValues)  .AND. &
            (.NOT. Me%UseOriginalValues)) then  
            Me%InterpolateValues = .true.     
        endif
        

        !When to shut down DT?
        call GetData(Me%MinForDTDecrease, Me%ObjEnterData,  iflag,                      &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'MIN_FOR_DT_DECREASE',                            &
                     default        = AllMostZero,                                      &
                     ClientModule   = 'ModuleFillMatrix',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR150'
        
        !Assumes the last instant value to avoid linear interpolation
        call GetData(Me%PreviousInstantValues, Me%ObjEnterData,  iflag,                 &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'PREVIOUS_INSTANT_VALUES',                        &
                     default        = .false.,                                          &
                     ClientModule   = 'ModuleFillMatrix',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR155'

        call GetData(Me%IgnoreNoDataPoint,                                              &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'IGNORE_NODATA_POINT',                            &
                     default        = .true.,                                           &
                     ClientModule   = 'ModuleFillMatrix',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR156' 
        
        if (.NOT. Me%IgnoreNoDataPoint) then
            call GetData(Me%NoDataValue,                                                &
                         Me%ObjEnterData,  iflag,                                       &
                         SearchType     = ExtractType,                                  &
                         keyword        = 'NODATA_VALUE',                               &
                         default        = -99.0,                                        &
                         ClientModule   = 'ModuleFillMatrix',                           &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR157'                  
        endif
        

        call GetData(Me%PredictDTMethod,                                            &
                     Me%ObjEnterData,  iflag,                                       &
                     SearchType     = ExtractType,                                  &
                     keyword        = 'PREDICT_DT_METHOD',                          &
                     default        = PredictDTMethod,                              &
                     ClientModule   = 'ModuleFillMatrix',                           &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR158'

        !Fill Matrix with default Value. This is the user field (not rotated to grid)
        if (Me%VectorialProp) then
            if (Me%Dim == Dim2D) then
                where (PointsToFill2D == WaterPoint) Me%Matrix2DX = Me%DefaultValue(1)
                where (PointsToFill2D == WaterPoint) Me%Matrix2DY = Me%DefaultValue(2)
            else
                where (PointsToFill3D == WaterPoint) Me%Matrix3DX = Me%DefaultValue(1)
                where (PointsToFill3D == WaterPoint) Me%Matrix3DY = Me%DefaultValue(2)
            endif            
        
        !User field angle (not rotated to grid)
        else if (Me%RotateAngleToGrid) then
            if (Me%Dim == Dim2D) then
                where (PointsToFill2D == WaterPoint) Me%Matrix2DFieldAngle = Me%DefaultValue(1)
            else
                where (PointsToFill3D == WaterPoint) Me%Matrix3DFieldAngle = Me%DefaultValue(1)
            endif            
        
        else
            if (Me%Dim == Dim2D) then
                where (PointsToFill2D == WaterPoint) Me%Matrix2D = Me%DefaultValue(1)
            else
                where (PointsToFill3D == WaterPoint) Me%Matrix3D = Me%DefaultValue(1)
            endif
        endif
        
        CanContinue  = .true.
        if (Me%VectorialProp) then
            if (Me%TimeEvolution == None) then            
                if ((Me%InitializationDefault /= Constant) .or. (Me%InitializationMethod /= Constant            &
                    .and. Me%InitializationMethod /= ASCIIFile .and. Me%InitializationMethod /= ReadTimeSerie  &
                    .and. Me%InitializationMethod /= ReadHDF)) then
                    CanContinue = .false.
                endif
            elseif (Me%TimeEvolution /= ReadTimeSerie .and. Me%TimeEvolution /= ReadHDF) then
                CanContinue = .false.
            endif
        endif
        if (.not. CanContinue) then
            write(*,*) 'For now vectorial property can only be defined by constant, asciifile, timeserie or HDF'
            write(*,*) 'Property : ', trim(Me%PropertyID%Name)
            stop ' ReadOptios - ModuleFillMatrix - ERR159'
        endif
    
    
        select case (Me%TimeEvolution)

            case(None)

                !Fill Matrix with an initialization default
                select case(Me%InitializationDefault)

                    case(Constant)

                        !Do nothing. Matrix values equal default value

                    case(ProfileTimeSerie)

                        if (Me%Dim == Dim2D) then
                            write(*,*)'Cannot Initialize 2D matrix by profile time serie'
                            stop 'ReadOptions - ModuleFillMatrix - ERR160'
                        else
                            call ConstructProfileTSDefault (PointsToFill3D, ExtractType)
                        endif
    
                end select

                select case(Me%InitializationMethod)

                    case(Constant)

                        !Do nothing. Matrix values equal default value

                    case(Layers)

                        if (Me%Dim == Dim2D) then
                            write(*,*)'Cannot Initialise 2D Matrix by Layers'
                            stop 'ReadOptions - ModuleFillMatrix - ERR170'
                        else
                            call ConstructSpaceLayers (ExtractType, PointsToFill3D)
                        endif

                    case(Boxes)

                        if (Me%Dim == Dim2D) then
                            call ConstructSpaceBox (ExtractType, PointsToFill2D = PointsToFill2D)
                        else
                            call ConstructSpaceBox (ExtractType, PointsToFill3D = PointsToFill3D)
                        endif

                    case(AsciiFile)
                        
                        if (Me%Dim == Dim2D) then
                            call ConstructSpaceASCIIFile (ExtractType, PointsToFill2D = PointsToFill2D)
                        else
                            call ConstructSpaceASCIIFile (ExtractType, PointsToFill3D = PointsToFill3D)
                        endif


                    case(Sponge)
                        
                        if (Me%Dim == Dim2D) then
                            call ConstructSponge (ExtractType, PointsToFill2D = PointsToFill2D)
                        else
                            call ConstructSponge (ExtractType, PointsToFill3D = PointsToFill3D)
                        endif

                    case(Profile)

                        if (Me%Dim == Dim2D) then
                            write(*,*)'Cannot Initialise 2D Matrix by Profile'
                            stop 'ReadOptions - ModuleFillMatrix - ERR180'
                        else
                            call ConstructSpaceProfile(ExtractType, PointsToFill3D)
                        endif


                    case(AnalyticProfile)

                        if (Me%Dim == Dim2D) then
                            write(*,*)'Cannot Initialise 2D Matrix by Profile'
                            stop 'ReadOptions - ModuleFillMatrix - ERR190'
                        else
                            call ConstructAnalyticProfile(ExtractType, PointsToFill3D)
                        endif

                    case(ReadTimeSerie)

                        if (Me%Dim == Dim2D) then
                            call ConstructSpaceTimeSerie    (ExtractType, PointsToFill2D = PointsToFill2D) 
                        else
                            call ConstructSpaceTimeSerie    (ExtractType, PointsToFill3D = PointsToFill3D) 
                        endif                        

                    case(ReadHDF)

                        if (Me%Dim == Dim2D) then
                            call ConstructHDFInput (ExtractType, ClientID_, PointsToFill2D = PointsToFill2D)
                        else    
                            call ConstructHDFInput (ExtractType, ClientID_, PointsToFill3D = PointsToFill3D)
                        endif

                    case(ProfileTimeSerie)

                        call ConstructProfileTimeSerie (PointsToFill3D, ExtractType)

                    end select

            case(ReadTimeSerie)

                if (Me%Dim == Dim2D) then
                    call ConstructSpaceTimeSerie    (ExtractType, PointsToFill2D = PointsToFill2D) 
                else
                    call ConstructSpaceTimeSerie    (ExtractType, PointsToFill3D = PointsToFill3D) 
                endif

            case(ReadHDF)

                if (Me%Dim == Dim2D) then
                    call ConstructHDFInput (ExtractType, ClientID_, PointsToFill2D = PointsToFill2D)
                else    
                    call ConstructHDFInput (ExtractType, ClientID_, PointsToFill3D = PointsToFill3D)
                endif

            case(ProfileTimeSerie)

                if (Me%Dim == Dim2D) then
                    write(*,*)'Cannot read 2D matrix by profile time serie'
                    stop 'ReadOptions - ModuleFillMatrix - ERR200'
                else
                    call ConstructProfileTimeSerie (PointsToFill3D, ExtractType)
                endif
                
            case(MultiTimeserie)
            
                if (.not. present(ClientID)) &
                    stop 'ReadOptions - ModuleFillMatrix - ERR210'                            

                if (Me%Dim == Dim2D) then
                    call ConstructMultiTimeSerie (ClientID, PointsToFill2D = PointsToFill2D)
                    call ModifyMultiTimeSerie    (PointsToFill2D = PointsToFill2D) 
                else
                    call ConstructMultiTimeSerie (ClientID, PointsToFill3D = PointsToFill3D)
                    call ModifyMultiTimeSerie    (PointsToFill3D = PointsToFill3D) 
                endif     
                
            case(AnalyticWave)
                call ConstructAnalyticWave (ExtractType)
                if (Me%Dim == Dim2D) then
                    call ModifyAnalyticWave    (PointsToFill2D = PointsToFill2D) 
                else
                    call ModifyAnalyticWave    (PointsToFill3D = PointsToFill3D) 
                endif     
                    
                           
                
        end select

    end subroutine ReadOptions

    !--------------------------------------------------------------------------
    
    subroutine InsertTimeSerieToList(TimeSerie)
        !Arguments-------------------------------------------------------------
        type (T_TimeSerie), pointer                    :: TimeSerie
        !Local-----------------------------------------------------------------
        type (T_TimeSerie), pointer                    :: CurrentTimeSerie  => null()
        type (T_TimeSerie), pointer                    :: PreviousTimeSerie => null()
        

        !Inserts a new property to the list of properties
        if (.not. associated(Me%FirstTimeSerie)) then
            Me%FirstTimeSerie => TimeSerie
            CurrentTimeSerie  => TimeSerie
        else
            PreviousTimeSerie => Me%FirstTimeSerie
            CurrentTimeSerie  => PreviousTimeSerie%Next
            do while (associated(CurrentTimeSerie))
                PreviousTimeSerie => CurrentTimeSerie
                CurrentTimeSerie  => PreviousTimeSerie%Next
            enddo
            PreviousTimeSerie%Next => TimeSerie
            TimeSerie%Prev      => PreviousTimeSerie
        endif        
        
        Me%nTimeSeries = Me%nTimeSeries + 1        
        
    end subroutine InsertTimeSerieToList
    
    !--------------------------------------------------------------------------
    
    
    subroutine InsertASCIIFileToList(ASCIIFile)
        !Arguments-------------------------------------------------------------
        type (T_ASCIIFile), pointer                    :: ASCIIFile
        !Local-----------------------------------------------------------------
        type (T_ASCIIFile), pointer                    :: CurrentASCIIFile  => null()
        type (T_ASCIIFile), pointer                    :: PreviousASCIIFile => null()
        
        
        !Inserts a new property to the list of properties
        if (.not. associated(Me%FirstASCIIFile)) then
            Me%FirstASCIIFile => ASCIIFile
            CurrentASCIIFile  => ASCIIFile
        else
            PreviousASCIIFile => Me%FirstASCIIFile
            CurrentASCIIFile  => PreviousASCIIFile%Next
            do while (associated(CurrentASCIIFile))
                PreviousASCIIFile => CurrentASCIIFile
                CurrentASCIIFile  => PreviousASCIIFile%Next
            enddo
            PreviousASCIIFile%Next => ASCIIFile
            ASCIIFile%Prev      => PreviousASCIIFile
        endif        

        
        Me%nASCIIFiles = Me%nASCIIFiles + 1        
        
    end subroutine InsertASCIIFileToList
    
    !--------------------------------------------------------------------------    

    subroutine InsertHDFToList(HDF)
        !Arguments-------------------------------------------------------------
        type (T_Field4D), pointer                      :: HDF
        !Local-----------------------------------------------------------------
        type (T_Field4D), pointer                      :: CurrentHDF  => null()
        type (T_Field4D), pointer                      :: PreviousHDF => null()
        
        !Inserts a new property to the list of properties
        if (.not. associated(Me%FirstHDF)) then
            Me%FirstHDF => HDF
            CurrentHDF  => HDF
        else
            PreviousHDF => Me%FirstHDF
            CurrentHDF  => PreviousHDF%Next
            do while (associated(CurrentHDF))
                PreviousHDF => CurrentHDF
                CurrentHDF  => PreviousHDF%Next
            enddo
            PreviousHDF%Next => HDF
            HDF%Prev         => PreviousHDF
        endif
        
        
        Me%nHDFs = Me%nHDFs + 1        
        
    end subroutine InsertHDFToList
    
    !--------------------------------------------------------------------------
        
    
    subroutine ConstructMultiTimeSerie (ClientID, PointsToFill2D, PointsToFill3D)
    
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: ClientID   
        integer, dimension(:, :), pointer, optional     :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D     

        !Local----------------------------------------------------------------
        character(PathLength)                           :: FileName
        integer                                         :: MaskGridDataID = 0
        integer                                         :: TypeZUV
        real, dimension(:, :),    pointer               :: GridData2D
        real, dimension(:, :, :), pointer               :: GridData3D
        integer                                         :: STAT_CALL
        integer                                         :: iflag
        integer                                         :: ilb, iub, jlb, jub, klb, kub
        integer                                         :: i, j, k
        integer                                         :: index
        logical                                         :: found
        type(T_Station), dimension(:), pointer          :: sl
        integer                                         :: ClientNumber
        Type(T_Time)                                    :: CurrentTime

        !----------------------------------------------------------------------    
        
        ClientNumber = ClientID
        
        !Gets the number of Source blocks
        call GetNumberOfBlocks(Me%ObjEnterData, BeginMTSBlock, EndMTSBlock,   &
                               FromBlock_, Me%MultiTimeSerie%NumberOfSources,  &
                               ClientNumber, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR010'        
        
        allocate (Me%MultiTimeSerie%StationsList(Me%MultiTimeSerie%NumberOfSources))
        
        sl => Me%MultiTimeSerie%StationsList        
        do index = 1, Me%MultiTimeSerie%NumberOfSources

            call ExtractBlockFromBlock(Me%ObjEnterData,                      &
                                       ClientNumber      = ClientNumber,     &
                                       block_begin       = BeginMTSBlock,    &
                                       block_end         = EndMTSBlock,      &
                                       BlockInBlockFound = found,            &
                                       STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR020'

            if (found) then
            
                call GetData(sl(index)%FillID,                             &
                             Me%ObjEnterData , iflag,                      &
                             SearchType   = FromBlockInBlock,              &
                             keyword      = 'FILL_ID',                     &
                             ClientModule = 'ModuleFillMatrix',            &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL /= SUCCESS_) &
                    stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR021'
                if (iflag /= 1) &
                    stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR022'
                    
                call GetData(sl(index)%RemainConstant,                     &
                             Me%ObjEnterData, iflag,                       &
                             SearchType   = FromBlockInBlock,              &
                             keyword      = 'REMAIN_CONSTANT',             &
                             default      = .false.,                       &
                             ClientModule = 'ModuleFillMatrix',            &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR023'
                    
                if (sl(index)%RemainConstant) then
                
                    call GetData(sl(index)%NewValue,                       &
                                 Me%ObjEnterData, iflag,                   &
                                 SearchType   = FromBlockInBlock,          &
                                 keyword      = 'DEFAULTVALUE',            &                                 
                                 ClientModule = 'ModuleFillMatrix',        &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR024'                
                
                else

                    call GetData(sl(index)%FileName,                           &
                                 Me%ObjEnterData, iflag,                       &
                                 SearchType   = FromBlockInBlock,              &
                                 keyword      = 'FILE_NAME',                   &                             
                                 ClientModule = 'ModuleFillMatrix',            &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR025'
                        
                    call GetData(sl(index)%Column,                             &
                                 Me%ObjEnterData , iflag,                      &
                                 SearchType   = FromBlockInBlock,              &
                                 keyword      = 'DATA_COLUMN',                 &
                                 ClientModule = 'ModuleFillMatrix',            &
                                 STAT         = STAT_CALL)                                      
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR026'
                    if (iflag /= 1) &
                        stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR027'                

                    !Starts Time Serie
                    call StartTimeSerieInput(sl(index)%ObjTimeSerie,            &
                                             sl(index)%FileName,                &
                                             Me%ObjTime,                        &
                                             CheckDates = Me%CheckDates,        &
                                             STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR028'
                endif
                    
            else

                stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR030'

            endif

        enddo   

        !call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
        !if (STAT_CALL /= SUCCESS_) &
        !    stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR040'
                 
        !Gets the name of the mask file
        call GetData(FileName,                                                   &
                     Me%ObjEnterData, iflag,                                     &
                     SearchType   = FromBlock,                                   &
                     keyword      = 'ASCII_MASK_FILENAME',                       &
                     ClientModule = 'ModuleFillMatrix',                          &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR050'
        if (iflag /= 1)then
            write(*,*)'ASCII_MASK_FILENAME is missing'
            stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR051'
        end if                
        
        call GetData(Me%MultiTimeSerie%DataProcessing,             &
                     Me%ObjEnterData , iflag,                      &
                     SearchType   = FromBlock,                     &                             
                     keyword      = 'DATA_PROCESSING',             &
                     Default      = Interpolate,                   & 
                     ClientModule = 'ModuleFillMatrix',            &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR060'
        if (iflag /= 1) &
            stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR070'         
                
        if (Me%Dim == Dim2D) then

            ilb = Me%WorkSize2D%ILB
            iub = Me%WorkSize2D%IUB
            jlb = Me%WorkSize2D%JLB
            jub = Me%WorkSize2D%JUB

            allocate (Me%MultiTimeSerie%FillGrid2D(ilb:iub, jlb:jub))

            call ConstructGridData(MaskGridDataID, Me%ObjHorizontalGrid,         &
                                   FileName     = FileName,                      &
                                   DefaultValue = Me%DefaultValue(1),            &
                                   STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR080'

            call GetGridDataType(MaskGridDataID, TypeZUV, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR090'

            if(TypeZUV .ne. Me%TypeZUV)then
                write(*,*)'Inconsistency found in type ZUV'
                write(*,*)'Grid data: '//trim(FileName)
                stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR100'
            end if

            call GetGridData (MaskGridDataID, GridData2D, STAT = STAT_CALL)  
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR110'
            
            do j = jlb, jub
            do i = ilb, iub
            
                if (PointsToFill2D(i,j) == WaterPoint) then
                
                    found = .false.
                    
                    do index = 1, Me%MultiTimeSerie%NumberOfSources
                        if (sl(index)%FillID == GridData2D(i,j)) then
                            found = .true.
                            exit
                        endif
                    enddo
                    
                    if (found) then
                        Me%MultiTimeSerie%FillGrid2D(i,j) = index
                    else
                        stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR120'
                    endif 
                    
                else
                
                    Me%MultiTimeSerie%FillGrid2D(i,j) = -99
                    
                endif
            
            enddo
            enddo

            call UnGetGridData    (MaskGridDataID, GridData2D, STAT = STAT_CALL)  
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR130'

        else

            ilb = Me%WorkSize3D%ILB
            iub = Me%WorkSize3D%IUB
            jlb = Me%WorkSize3D%JLB
            jub = Me%WorkSize3D%JUB
            klb = Me%WorkSize3D%KLB
            kub = Me%WorkSize3D%KUB
            
            allocate (Me%MultiTimeSerie%FillGrid3D(ilb:iub, jlb:jub, klb:kub))

            call ConstructGridData(MaskGridDataID, Me%ObjHorizontalGrid,         &
                                   FileName     = FileName,                      &
                                   KLB          = Me%Worksize3D%KLB,             &
                                   KUB          = Me%Worksize3D%KUB,             &
                                   DefaultValue = Me%DefaultValue(1),            &
                                   STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR140'

            call GetGridDataType(Me%ASCIIFile%GridDataID, TypeZUV, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR150'
            
            if(TypeZUV .ne. Me%TypeZUV)then
                write(*,*)'Inconsistency found in type ZUV'
                write(*,*)'Grid data: '//trim(FileName)
                stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR160'
            end if

            call GetGridData (MaskGridDataID, GridData3D, STAT = STAT_CALL)  
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR170'

            do j = jlb, jub
            do i = ilb, iub
            do k = klb, kub
            
                if (PointsToFill3D(i,j,k) == WaterPoint) then
                
                    found = .false.
                    
                    do index = 1, Me%MultiTimeSerie%NumberOfSources
                        if (sl(index)%FillID == GridData3D(i,j,k)) then
                            found = .true.
                            exit
                        endif
                    enddo
                    
                    if (found) then
                        Me%MultiTimeSerie%FillGrid3D(i,j,k) = index
                    else
                        stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR180'
                    endif 
                    
                else
                
                    Me%MultiTimeSerie%FillGrid3D(i,j,k) = -99
                    
                endif
            
            enddo
            enddo
            enddo

            call UnGetGridData    (MaskGridDataID, GridData3D, STAT = STAT_CALL)  
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR190'

        endif

        call KillGridData (MaskGridDataID, STAT = STAT_CALL)  
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR200'  
                
        if (Me%PredictDTMethod == 2) then                    

            !Gets Current Time
            call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR210'                      
        
            Me%DTForNextEvent       = -null_real
            Me%PredictedDT          = -null_real
            Me%DTForNextDataset     = -null_real
            Me%NextValueForDTPred   = 0.0
        
            sl => Me%MultiTimeSerie%StationsList
            do index = 1, Me%MultiTimeSerie%NumberOfSources
            
                if (.not. sl(index)%RemainConstant) then
            
                    call GetTimeSerieDataValues(sl(index)%ObjTimeSerie,     &
                                                sl(index)%NumberOfInstants, &
                                                STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR219'            
                
                    if(sl(index)%NumberOfInstants > 1)then                
                    
                        sl(index)%PreviousInstant = 1                
                        call GetTimeSerieTimeOfDataset(sl(index)%ObjTimeSerie,      &
                                                       sl(index)%PreviousInstant,   &
                                                       sl(index)%PreviousTime,      &
                                                       STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)    &
                            stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR220'
                                                         
                        sl(index)%NextInstant = 1                
                        sl(index)%NextTime    = sl(index)%PreviousTime
                                                                                        
                        call ActualizeMultiTimeSerieValues (sl(index))
                        
                        sl(index)%NextEventStart     = sl(index)%PreviousTime
                        sl(index)%NextEventEnd       = sl(index)%PreviousTime
                        sl(index)%NextValueForDTPred = sl(index)%NextValue
                        sl(index)%DTForNextEvent     = 0.0
                         
!                        
!                        if (Me%AccumulateValues) then
!                            if (sl(index)%NextValue > Me%MinForDTDecrease) then        
!                                sl(index)%NextEventStart = sl(index)%PreviousTime
!                                sl(index)%NextEventEnd   = sl(index)%NextTime
!                            else
!                                call FindNextEventInMultiTimeSerie(CurrentTime, sl(index))
!                            endif
!                        else
!                            sl(index)%DTForNextEvent = 0.0
!                        endif
!                        
!                        if (Me%ValueIsUsedForDTPrediction .and. (sl(index)%DTForNextEvent == 0.0)) then
!                            if (Me%UseOriginalValues) then
!                                if (sl(index)%NextValue > sl(index)%NextValueForDTPred) &
!                                    sl(index)%NextValueForDTPred = sl(index)%NextValue
!                            elseif (Me%AccumulateValues) then
!                                value = sl(index)%NextValue / (sl(index)%PreviousTime - sl(index)%NextTime)
!                                if (value > Me%NextValueForDTPred) &
!                                    sl(index)%NextValueForDTPred = value
!                            else
!                                stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR240'
!                            endif
!                        endif

                    else
                    
                        stop 'ConstructMultiTimeSerie - ModuleFillMatrix - ERR241'

                    endif
                    
                else
                    
                    if (sl(index)%NewValue > 0.0) then
                        sl(index)%DTForNextEvent        = 0.0
                        sl(index)%PredictedDT           = -null_real
                        sl(index)%DTForNextDataset      = -null_real
                        sl(index)%NextValueForDTPred    = sl(index)%NewValue
                    else
                        sl(index)%DTForNextEvent        = -null_real
                        sl(index)%PredictedDT           = -null_real
                        sl(index)%DTForNextDataset      = -null_real
                        sl(index)%NextValueForDTPred    = 0.0
                    endif
                    
                endif
!  
!                if (Me%DTForNextEvent > sl(index)%DTForNextEvent) &
!                    Me%DTForNextEvent = sl(index)%DTForNextEvent
!                
!                if (Me%PredictedDT > sl(index)%PredictedDT) &
!                    Me%PredictedDT = sl(index)%PredictedDT 
!
!                if (Me%DTForNextDataset > sl(index)%DTForNextDataset) &
!                    Me%DTForNextDataset = sl(index)%DTForNextDataset
!
!                if ((sl(index)%DTForNextEvent == 0.0) .and. &
!                    (Me%NextValueForDTPred < sl(index)%NextValueForDTPred)) then
!                    Me%NextValueForDTPred = sl(index)%NextValueForDTPred
!                endif
  
            enddo
        
        endif
                
        !----------------------------------------------------------------------
    
    end subroutine ConstructMultiTimeSerie
    
    !--------------------------------------------------------------------------

    subroutine ConstructSpaceLayers (ExtractType, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        integer                                     :: ExtractType
        integer, dimension(:, :, :), pointer        :: PointsToFill3D

        !Local----------------------------------------------------------------
        character(PathLength)                       :: FileName
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer                                     :: k, NLayers
        integer                                     :: ObjEnterData = 0
        integer                                     :: ClientNumber
        integer                                     :: FirstLine, LastLine
        logical                                     :: BlockFound
        integer                                     :: line, l
        real, dimension(:), allocatable             :: Aux 

        !----------------------------------------------------------------------
        !Begin----------------------------------------------------------------

        allocate (Me%Layers%Values(Me%WorkSize3D%KUB))


        call GetData(FileName,                                                   &
                     Me%ObjEnterData , iflag,                                    &
                     SearchType   = ExtractType,                                 &
                     keyword      = 'FILENAME',                                  &
                     ClientModule = 'ModuleFillMatrix',                          &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceLayers - ModuleFillMatrix - ERR01'


flag0:  if (iflag /= 0) then
            
    
            !Opens File
            call ConstructEnterData(ObjEnterData, FileName, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ConstructSpaceLayers - ModuleFillMatrix - ERR02'


            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                          &
                                        BeginLayers, EndLayers, BlockFound,                  &
                                        FirstLine = FirstLine, LastLine = LastLine,          &
                                        STAT = STAT_CALL)
        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceLayers - ModuleFillMatrix - ERR03'

BF:         if (BlockFound) then

                NLayers =  LastLine - FirstLine - 1

                if (NLayers > Me%WorkSize3D%KUB) stop 'ConstructSpaceLayers - ModuleFillMatrix - ERR04'


                !Allocates auxiliar variables
                allocate (Aux(2))
            
                l = 1
                do line = FirstLine + 1, LastLine - 1

                    call GetData(Aux, EnterDataID = ObjEnterData, flag = iflag, &
                                    SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceLayers - ModuleFillMatrix - ERR05'

                    if (Aux(1) < Me%WorkSize3D%KLB .or. Aux(1) > Me%WorkSize3D%KUB)      &
                        stop 'ConstructSpaceLayers - ModuleFillMatrix - ERR06'

                    Me%Layers%Values(int(Aux(1))) = Aux(2)
                    l = l + 1

                enddo

                deallocate(Aux)

                call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceLayers - ModuleFillMatrix - ERR07'

                call KillEnterData(ObjEnterData, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceLayers - ModuleFillMatrix - ERR08'


            else 

                write(*,*) 'Block <BeginLayers>, <EndLayers> not found'
                write(*,*) 'FileName = ', trim(FileName)
                stop 'ConstructSpaceLayers - ModuleFillMatrix - ERR09'

            endif BF                            

        else if (iflag == 0) then flag0

            call GetData(Me%Layers%Values,                                           &
                         Me%ObjEnterData, iflag,                                     &
                         SearchType    = ExtractType,                                &
                         keyword       = 'LAYERS_VALUES',                            &
                         ClientModule  = 'ModuleFillMatrix',                         &
                         STAT          = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceLayers - ModuleFillMatrix - ERR10'

            if (iflag /= Me%WorkSize3D%KUB) then
                write(*,*)'Invalid Number of Layers'
                stop 'ConstructSpaceLayers - ModuleFillMatrix - ERR11'
            endif
            
            
        endif flag0

        !Fill Matrix with default Value
        do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            where (PointsToFill3D(:,:,k) == WaterPoint) Me%Matrix3D(:,:,k) = Me%Layers%Values(k)
        enddo   

    end subroutine ConstructSpaceLayers

    !--------------------------------------------------------------------------

    subroutine ConstructSpaceProfile (ExtractType, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        integer                                     :: ExtractType
        integer, dimension(:, :, :), pointer        :: PointsToFill3D

        !Local----------------------------------------------------------------
        real, dimension(:,:,:), pointer             :: SZZ
        real, dimension(:),     pointer             :: Values, Depth
        real                                        :: CellDepth
        character(PathLength)                       :: FileName
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer                                     :: i, j, k, NDEPTHS
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: ObjEnterData = 0
        integer                                     :: ClientNumber
        integer                                     :: FirstLine, LastLine
        logical                                     :: BlockFound
        integer                                     :: line, l
        real, dimension(:), allocatable             :: Aux 

        !----------------------------------------------------------------------
        !Begin----------------------------------------------------------------

        ILB = Me%WorkSize3D%ILB
        IUB = Me%WorkSize3D%IUB

        JLB = Me%WorkSize3D%JLB
        JUB = Me%WorkSize3D%JUB

        KLB = Me%WorkSize3D%KLB
        KUB = Me%WorkSize3D%KUB

        !Gets Center of the cells
        call GetGeometryDistances(Me%ObjGeometry, SZZ = SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR01'


        call GetData(FileName,                                                   &
                     Me%ObjEnterData , iflag,                                    &
                     SearchType   = ExtractType,                                 &
                     keyword      = 'FILENAME',                                  &
                     ClientModule = 'ModuleFillMatrix',                          &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR02'


flag0:  if (iflag /= 0) then

            !Opens File
            call ConstructEnterData(ObjEnterData, FileName, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR03'


            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                          &
                                        BeginProfile, EndProfile, BlockFound,                &
                                        FirstLine = FirstLine, LastLine = LastLine,          &
                                        STAT = STAT_CALL)
        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR04'

BF:         if (BlockFound) then

                NDEPTHS =  LastLine - FirstLine - 1

                !Allocates auxiliar variables
                allocate (Values        (NDEPTHS))
                allocate (Depth         (NDEPTHS))
                allocate (Aux(2))
            
                l = 1
                do line = FirstLine + 1, LastLine - 1

                    call GetData(Aux, EnterDataID = ObjEnterData, flag = iflag, &
                                 SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR05'

                    Depth (l) = Aux(1)
                    Values(l) = Aux(2)
                    l = l + 1

                enddo

                deallocate(Aux)

                call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR06'

                call KillEnterData(ObjEnterData, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR07'


            else 

                write(*,*) 'Block <BeginProfile>, <EndProfile> not found'
                write(*,*) 'FileName = ', trim(FileName)
                stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR08'

            endif BF

        else if (iflag == 0) then flag0

            !Get the number of values
            call GetData(NDEPTHS, Me%ObjEnterData, iflag,                                    &
                         SearchType   = ExtractType,                                         &
                         keyword      = 'NDEPTHS',                                           &
                         ClientModule = 'FillMatrix',                                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR09'

            !Allocates auxiliar variables
            allocate (Values        (NDEPTHS))
            allocate (Depth         (NDEPTHS))

            !Gets Depth
            call GetData (Depth, Me%ObjEnterData, iflag,                                     &
                          SearchType   = ExtractType,                                        &
                          keyword      = 'DEPTH_PROFILE',                                    &
                          ClientModule = 'FillMatrix',                                       &
                          STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR10'


            if (iflag /= NDEPTHS) then
                write (*,*) 'Invalid Number of Depth Values. Keyword DEPTH_PROFILE'
                if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR11'
            endif

        
            !Gets Profile
            call GetData (Values, Me%ObjEnterData, iflag,                                    &
                          SearchType   = ExtractType,                                        &
                          keyword      = 'PROFILE_VALUES',                                   &
                          ClientModule = 'FillMatrix',                                       &
                          STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR12'

            if (iflag /= NDEPTHS) then
                write (*,*) 'Invalid Number of Values. Keyword PROFILE_VALUES'
                if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR13'
            endif

        endif flag0

        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB

            if (PointsToFill3D(i, j, k) == WaterPoint) then

                CellDepth          = (SZZ(i, j, k) + SZZ(i, j, k - 1)) / 2 - SZZ(i, j, KUB) 
                Me%Matrix3D(i,j,k) = InterpolateProfile (CellDepth, NDEPTHS, Depth, Values)

            endif

        enddo
        enddo
        enddo


        !Deallocates auxiliar variables
        deallocate (Values)
        deallocate (Depth )

        call UnGetGeometry (Me%ObjGeometry, SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - ModuleFillMatrix - ERR14'


    end subroutine ConstructSpaceProfile

    !--------------------------------------------------------------------------

    subroutine ConstructProfileTimeSerie (PointsToFill3D, ExtractType)

        !Arguments-------------------------------------------------------------
        integer                                     :: ExtractType
        integer, dimension(:, :, :), pointer        :: PointsToFill3D

        !Local----------------------------------------------------------------
        real,           dimension(:,:,:), pointer   :: SZZ
        type(T_Time)                                :: AuxTime
        character(PathLength)                       :: FileName
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: ObjEnterData = 0
        integer                                     :: ClientNumber
        integer                                     :: FirstLine, LastLine
        logical                                     :: BlockFound
        logical                                     :: FoundSecondInstant
        integer                                     :: line, l
        real, dimension(:), allocatable             :: Aux 
        type(T_Time)                                :: NextTime, PreviousTime, Now
        type(T_Time)                                :: LastInstantTime, EndTime

        !Begin----------------------------------------------------------------

        ILB = Me%Size3D%ILB
        IUB = Me%Size3D%IUB
        JLB = Me%Size3D%JLB
        JUB = Me%Size3D%JUB
        KLB = Me%Size3D%KLB
        KUB = Me%Size3D%KUB

        !Gets Center of the cells
        call GetGeometryDistances(Me%ObjGeometry, SZZ = SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR10'

        call GetData(FileName,                                                   &
                     Me%ObjEnterData , iflag,                                    &
                     SearchType   = ExtractType,                                 &
                     keyword      = 'FILENAME',                                  &
                     ClientModule = 'ModuleFillMatrix',                          &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR20'


flag0:  if (iflag /= 0) then

            !Opens File
            call ConstructEnterData(ObjEnterData, FileName, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR30'

            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                             &
                                        BeginTimes, EndTimes, BlockFound,                       &
                                        FirstLine = FirstLine, LastLine = LastLine,             &
                                        STAT = STAT_CALL)
        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR40'

            if (BlockFound) then

                Me%ProfileTimeSerie%NumberOfInstants =  LastLine - FirstLine - 1

               !Allocates auxiliar variables
                allocate (Me%ProfileTimeSerie%TimeInstants (Me%ProfileTimeSerie%NumberOfInstants))

                l = 1
                do line = FirstLine + 1, LastLine - 1

                    call GetData(AuxTime, EnterDataID = ObjEnterData, flag = iflag, &
                                 SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR50'

                    Me%ProfileTimeSerie%TimeInstants (l) = AuxTime

                    l = l + 1

                enddo

                call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR60'

            else 

                write(*,*) 'Block <BeginTime>, <EndTime> not found'
                write(*,*) 'FileName = ', trim(FileName)
                stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR70'

            endif

            call RewindBuffer(ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR80'



            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                          &
                                        BeginProfileValues, EndProfileValues, BlockFound,    &
                                        FirstLine = FirstLine, LastLine = LastLine,          &
                                        STAT = STAT_CALL)
        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR90'

BF:         if (BlockFound) then

                allocate(Aux(1:Me%ProfileTimeSerie%NumberOfInstants))

                call GetData(Aux, EnterDataID = ObjEnterData, flag = iflag, &
                                  SearchType = FromBlock, Buffer_Line = FirstLine+1, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR100'

                Me%ProfileTimeSerie%nValues =  LastLine - FirstLine - 1

                if(iflag .ne. Me%ProfileTimeSerie%NumberOfInstants)then
                    write(*,*) 'Number of time instants is not equal to number of values'
                    write(*,*) 'FileName = ', trim(FileName)
                    stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR110'

                endif

 
                !Allocates auxiliar variables
                allocate (Me%ProfileTimeSerie%Values (1:Me%ProfileTimeSerie%nValues,         &
                                                      1:Me%ProfileTimeSerie%NumberOfInstants))

                l = 1
                do line = FirstLine + 1, LastLine - 1

                    call GetData(Aux, EnterDataID = ObjEnterData, flag = iflag, &
                                 SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR120'

                    Me%ProfileTimeSerie%Values(l,:) = Aux(:)

                    l = l + 1

                enddo

                deallocate(Aux)

                call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR130'

                call RewindBuffer(ObjEnterData, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR140'

            else 

                write(*,*) 'Block <BeginProfileValues>, <EndProfileValues> not found'
                write(*,*) 'FileName = ', trim(FileName)
                stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR150'

            endif BF

            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                          &
                                        BeginDepth, EndDepth, BlockFound,                    &
                                        FirstLine = FirstLine, LastLine = LastLine,          &
                                        STAT = STAT_CALL)
        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR160'

            if (BlockFound) then

                allocate(Aux(1:Me%ProfileTimeSerie%NumberOfInstants))

                Me%ProfileTimeSerie%nDepths =  LastLine - FirstLine - 1

                call GetData(Aux, EnterDataID = ObjEnterData, flag = iflag, &
                                 SearchType = FromBlock, Buffer_Line = FirstLine+1, STAT = STAT_CALL)

                if(iflag .ne. Me%ProfileTimeSerie%NumberOfInstants)then
                    write(*,*) 'Number of depths is not equal to number of values'
                    write(*,*) 'FileName = ', trim(FileName)
                    stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR170'

                endif

                !Allocates auxiliar variables
                allocate (Me%ProfileTimeSerie%Depths (1:Me%ProfileTimeSerie%nDepths,         &
                                                      1:Me%ProfileTimeSerie%NumberOfInstants))

                l = 1
                do line = FirstLine + 1, LastLine - 1

                    call GetData(Aux, EnterDataID = ObjEnterData, flag = iflag, &
                                 SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR180'

                    Me%ProfileTimeSerie%Depths (l,:) = Aux(:)

                    l = l + 1

                enddo

                deallocate(Aux)

                call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR190'

            else 

                write(*,*) 'Block <BeginProfileValues>, <EndProfileValues> not found'
                write(*,*) 'FileName = ', trim(FileName)
                stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR200'

            endif



            call KillEnterData(ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR210'


        else if (iflag == 0) then flag0


                write(*,*) 'Could not find profile time serie file.'
                write(*,*) 'FileName = ', trim(FileName)
                stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR220'

        endif flag0

        call GetComputeTimeLimits(Me%ObjTime, EndTime = EndTime, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR230'

        call GetComputeCurrentTime(Me%ObjTime, Now, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR240'

        if(Me%ProfileTimeSerie%NumberOfInstants > 1)then

            Me%ProfileTimeSerie%PreviousInstant  = 1
            Me%ProfileTimeSerie%NextInstant      = Me%ProfileTimeSerie%PreviousInstant
            PreviousTime                         = Me%ProfileTimeSerie%TimeInstants(Me%ProfileTimeSerie%PreviousInstant)
            call CheckCyclicMonths(PreviousTime, RefTime = Now, CyclicTimeON = Me%ProfileTimeSerie%CyclicTimeON)


i23:        if (Me%ProfileTimeSerie%CyclicTimeON) then

                if (Me%ProfileTimeSerie%NumberOfInstants /= 12) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR245'

            else i23

                if(PreviousTime .gt. Now)then
                    write(*,*)
                    write(*,*)'Could not create solution from profile time serie file'
                    write(*,*)'First file instant greater than current time'
                    write(*,*)'Matrix name: '//trim(Me%PropertyID%Name)
                    stop      'ConstructProfileTimeSerie - ModuleFillMatrix - ERR250'
                end if
        
                LastInstantTime         = Me%ProfileTimeSerie%TimeInstants(Me%ProfileTimeSerie%NumberOfInstants)

                if (Me%ProfileTimeSerie%CyclicTimeON)  call CheckCyclicMonths(LastInstantTime, RefTime = EndTime)

                if(LastInstantTime .lt. EndTime)then
                    write(*,*)
                    write(*,*)'Could not create solution from profile time serie file'
                    write(*,*)'Last instant in file lower than simulation ending time'
                    write(*,*)'Matrix name: '//trim(Me%PropertyID%Name)
                    stop      'ConstructProfileTimeSerie - ModuleFillMatrix - ERR260'
                end if

            endif i23

            FoundSecondInstant = .false.
            
            !if number of instants greater than 1 then 
            !find first and second instants
            do while(.not. FoundSecondInstant)
                
                Me%ProfileTimeSerie%PreviousInstant  = Me%ProfileTimeSerie%NextInstant
                Me%ProfileTimeSerie%NextInstant      = Me%ProfileTimeSerie%NextInstant + 1

                if (Me%ProfileTimeSerie%CyclicTimeON .and. Me%ProfileTimeSerie%NextInstant &
                    .gt. Me%ProfileTimeSerie%NumberOfInstants) then
                    Me%ProfileTimeSerie%NextInstant  = 1
                end if

                NextTime                = Me%ProfileTimeSerie%TimeInstants(Me%ProfileTimeSerie%NextInstant)
                
                if (Me%ProfileTimeSerie%CyclicTimeON) call CheckCyclicMonths(NextTime, PreviousTime = PreviousTime)

                if(PreviousTime .le. Now .and. NextTime .ge. Now) then
                    FoundSecondInstant  = .true.
                    exit
                end if

                PreviousTime            = NextTime

                if(Me%ProfileTimeSerie%NextInstant .gt. Me%ProfileTimeSerie%NumberOfInstants &
                   .and. .not. Me%ProfileTimeSerie%CyclicTimeON)then
                    write(*,*)
                    write(*,*)'Could not read solution from Profile Time Serie file'
                    write(*,*)'Could not find second instant in file'
                    write(*,*)'Matrix name: '//trim(Me%PropertyID%Name)
                    stop      'ConstructProfileTimeSerie - ModuleFillMatrix - ERR270'
                end if

            end do

            Me%ProfileTimeSerie%PreviousTime     = PreviousTime
            Me%ProfileTimeSerie%NextTime         = NextTime

        endif

        call UnGetGeometry (Me%ObjGeometry, SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR280'

        allocate(Me%ProfileTimeSerie%PreviousField3D (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%ProfileTimeSerie%NextField3D     (ILB:IUB, JLB:JUB, KLB:KUB))

        Me%ProfileTimeSerie%PreviousField3D(:,:,:) = FillValueReal
        Me%ProfileTimeSerie%NextField3D    (:,:,:) = FillValueReal

        call ProfileTimeSerieField(PointsToFill3D, Me%ProfileTimeSerie%PreviousInstant, &
                                                   Me%ProfileTimeSerie%PreviousField3D)
        call ProfileTimeSerieField(PointsToFill3D, Me%ProfileTimeSerie%NextInstant,     &
                                                   Me%ProfileTimeSerie%NextField3D    )

        if (Me%PropertyID%IsAngle) then            

            call InterpolateAngle3DInTime (ActualTime       = Now,                                   &
                                           Size             = Me%WorkSize3D,                         &
                                           Time1            = Me%ProfileTimeSerie%PreviousTime,      &
                                           Matrix1          = Me%ProfileTimeSerie%PreviousField3D,   &
                                           Time2            = Me%ProfileTimeSerie%NextTime,          &
                                           Matrix2          = Me%ProfileTimeSerie%NextField3D,       &
                                           MatrixOut        = Me%Matrix3D,                           &
                                           PointsToFill3D   = PointsToFill3D)

        else

            call InterpolateMatrix3DInTime(ActualTime       = Now,                                   &
                                           Size             = Me%WorkSize3D,                         &
                                           Time1            = Me%ProfileTimeSerie%PreviousTime,      &
                                           Matrix1          = Me%ProfileTimeSerie%PreviousField3D,   &
                                           Time2            = Me%ProfileTimeSerie%NextTime,          &
                                           Matrix2          = Me%ProfileTimeSerie%NextField3D,       &
                                           MatrixOut        = Me%Matrix3D,                           &
                                           PointsToFill3D   = PointsToFill3D)
                                           
        endif                                           

    end subroutine ConstructProfileTimeSerie

    !--------------------------------------------------------------------------

    subroutine ConstructProfileTSDefault (PointsToFill3D, ExtractType)

        !Arguments-------------------------------------------------------------
        integer                                     :: ExtractType
        integer, dimension(:, :, :), pointer        :: PointsToFill3D

        !Local----------------------------------------------------------------
        real,           dimension(:,:,:), pointer   :: SZZ
        type(T_Time)                                :: AuxTime
        character(PathLength)                       :: FileName
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: ObjEnterData = 0
        integer                                     :: ClientNumber
        integer                                     :: FirstLine, LastLine
        logical                                     :: BlockFound
        logical                                     :: FoundSecondInstant
        integer                                     :: line, l
        real, dimension(:), allocatable             :: Aux 
        type(T_Time)                                :: NextTime, PreviousTime, Now
        type(T_Time)                                :: LastInstantTime, EndTime

        !Begin----------------------------------------------------------------

        ILB = Me%Size3D%ILB
        IUB = Me%Size3D%IUB
        JLB = Me%Size3D%JLB
        JUB = Me%Size3D%JUB
        KLB = Me%Size3D%KLB
        KUB = Me%Size3D%KUB

        !Gets Center of the cells
        call GetGeometryDistances(Me%ObjGeometry, SZZ = SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                               & 
            stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR10'

        call GetData(FileName,                                                   &
                     Me%ObjEnterData , iflag,                                    &
                     SearchType   = ExtractType,                                 &
                     keyword      = 'FILENAME_DEFAULT',                          &
                     ClientModule = 'ModuleFillMatrix',                          &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_)                                             &
            stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR20'


flag0:  if (iflag /= 0) then

            !Opens File
            call ConstructEnterData(ObjEnterData, FileName, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR30'

            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                             &
                                        BeginTimes, EndTimes, BlockFound,                       &
                                        FirstLine = FirstLine, LastLine = LastLine,             &
                                        STAT = STAT_CALL)
        
            if (STAT_CALL /= SUCCESS_)                                           &
                stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR40'

            if (BlockFound) then

                Me%ProfileTimeSerie%NumberOfInstants =  LastLine - FirstLine - 1

               !Allocates auxiliar variables
                allocate (Me%ProfileTimeSerie%TimeInstants (Me%ProfileTimeSerie%NumberOfInstants))

                l = 1
                do line = FirstLine + 1, LastLine - 1

                    call GetData(AuxTime, EnterDataID = ObjEnterData, flag = iflag, &
                                 SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_)                                   & 
                        stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR50'

                    Me%ProfileTimeSerie%TimeInstants (l) = AuxTime

                    l = l + 1

                enddo

                call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                       & 
                    stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR60'

            else 

                write(*,*) 'Block <BeginTime>, <EndTime> not found'
                write(*,*) 'FileName = ', trim(FileName)
                stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR70'

            endif

            call RewindBuffer(ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                           & 
                stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR80'



            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                          &
                                        BeginProfileValues, EndProfileValues, BlockFound,    &
                                        FirstLine = FirstLine, LastLine = LastLine,          &
                                        STAT = STAT_CALL)
        
            if (STAT_CALL /= SUCCESS_)                                           & 
                stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR90'

BF:         if (BlockFound) then

                allocate(Aux(1:Me%ProfileTimeSerie%NumberOfInstants))

                call GetData(Aux, EnterDataID = ObjEnterData, flag = iflag,      &
                                  SearchType = FromBlock, Buffer_Line = FirstLine+1, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                       & 
                    stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR100'

                Me%ProfileTimeSerie%nValues =  LastLine - FirstLine - 1

                if(iflag .ne. Me%ProfileTimeSerie%NumberOfInstants)then
                    write(*,*) 'Number of time instants is not equal to number of values'
                    write(*,*) 'FileName = ', trim(FileName)
                    stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR110'

                endif

 
                !Allocates auxiliar variables
                allocate (Me%ProfileTimeSerie%Values (1:Me%ProfileTimeSerie%nValues,         &
                                                      1:Me%ProfileTimeSerie%NumberOfInstants))

                l = 1
                do line = FirstLine + 1, LastLine - 1

                    call GetData(Aux, EnterDataID = ObjEnterData, flag = iflag,  &
                                 SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_)                                   &                              
                        stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR120'

                    Me%ProfileTimeSerie%Values(l,:) = Aux(:)

                    l = l + 1

                enddo

                deallocate(Aux)

                call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                       & 
                    stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR130'

                call RewindBuffer(ObjEnterData, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                       &
                    stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR140'

            else 

                write(*,*) 'Block <BeginProfileValues>, <EndProfileValues> not found'
                write(*,*) 'FileName = ', trim(FileName)
                stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR150'

            endif BF

            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                          &
                                        BeginDepth, EndDepth, BlockFound,                    &
                                        FirstLine = FirstLine, LastLine = LastLine,          &
                                        STAT = STAT_CALL)
        
            if (STAT_CALL /= SUCCESS_)                                           &
                stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR160'

            if (BlockFound) then

                allocate(Aux(1:Me%ProfileTimeSerie%NumberOfInstants))

                Me%ProfileTimeSerie%nDepths =  LastLine - FirstLine - 1

                call GetData(Aux, EnterDataID = ObjEnterData, flag = iflag,      &
                                 SearchType = FromBlock, Buffer_Line = FirstLine+1, STAT = STAT_CALL)

                if(iflag .ne. Me%ProfileTimeSerie%NumberOfInstants)then
                    write(*,*) 'Number of depths is not equal to number of values'
                    write(*,*) 'FileName = ', trim(FileName)
                    stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR170'

                endif

                !Allocates auxiliar variables
                allocate (Me%ProfileTimeSerie%Depths (1:Me%ProfileTimeSerie%nDepths,         &
                                                      1:Me%ProfileTimeSerie%NumberOfInstants))

                l = 1
                do line = FirstLine + 1, LastLine - 1

                    call GetData(Aux, EnterDataID = ObjEnterData, flag = iflag,  &
                                 SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_)                                   &
                        stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR180'

                    Me%ProfileTimeSerie%Depths (l,:) = Aux(:)

                    l = l + 1

                enddo

                deallocate(Aux)

                call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                       & 
                    stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR190'

            else 

                write(*,*) 'Block <BeginProfileValues>, <EndProfileValues> not found'
                write(*,*) 'FileName = ', trim(FileName)
                stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR200'

            endif



            call KillEnterData(ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                           &
                stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR210'


        else if (iflag == 0) then flag0


                write(*,*) 'Could not find profile time serie file.'
                write(*,*) 'FileName = ', trim(FileName)
                stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR220'

        endif flag0

        call GetComputeTimeLimits(Me%ObjTime, EndTime = EndTime, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                             & 
            stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR230'

        call GetComputeCurrentTime(Me%ObjTime, Now, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                             &
            stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR240'

        if(Me%ProfileTimeSerie%NumberOfInstants > 1)then

            Me%ProfileTimeSerie%PreviousInstant  = 1
            Me%ProfileTimeSerie%NextInstant      = Me%ProfileTimeSerie%PreviousInstant
            PreviousTime                         = Me%ProfileTimeSerie%TimeInstants(Me%ProfileTimeSerie%PreviousInstant)

            call CheckCyclicMonths(PreviousTime, RefTime = Now, CyclicTimeON = Me%ProfileTimeSerie%CyclicTimeON)


i23:        if (Me%ProfileTimeSerie%CyclicTimeON) then
    
                if (Me%ProfileTimeSerie%NumberOfInstants /= 12) stop 'ConstructProfileTimeSerie - ModuleFillMatrix - ERR245'

            else i23

                if(PreviousTime .gt. Now)then
                    write(*,*)
                    write(*,*)'Could not create solution from profile time serie file'
                    write(*,*)'First file instant greater than current time'
                    write(*,*)'Matrix name: '//trim(Me%PropertyID%Name)
                    stop      'ConstructProfileTSDefault - ModuleFillMatrix - ERR250'
                end if
            
                LastInstantTime         = Me%ProfileTimeSerie%TimeInstants(Me%ProfileTimeSerie%NumberOfInstants)

                if (Me%ProfileTimeSerie%CyclicTimeON) call CheckCyclicMonths(LastInstantTime, RefTime = EndTime)

    !            if(LastInstantTime .lt. EndTime)then
                if(LastInstantTime .lt. Now)then
                    write(*,*)
                    write(*,*)'Could not create solution from profile time serie file'
                    write(*,*)'Last instant in file lower than simulation ending time'
                    write(*,*)'Matrix name: '//trim(Me%PropertyID%Name)
                    stop      'ConstructProfileTSDefault - ModuleFillMatrix - ERR260'
                end if

            end if i23

            FoundSecondInstant = .false.
            
            !if number of instants greater than 1 then 
            !find first and second instants
            do while(.not. FoundSecondInstant)
                
                Me%ProfileTimeSerie%PreviousInstant  = Me%ProfileTimeSerie%NextInstant
                Me%ProfileTimeSerie%NextInstant      = Me%ProfileTimeSerie%NextInstant + 1

                if (Me%ProfileTimeSerie%CyclicTimeON .and. Me%ProfileTimeSerie%NextInstant .gt. &
                    Me%ProfileTimeSerie%NumberOfInstants) then
                    Me%ProfileTimeSerie%NextInstant  = 1
                end if

                NextTime                = Me%ProfileTimeSerie%TimeInstants(Me%ProfileTimeSerie%NextInstant)

                if (Me%ProfileTimeSerie%CyclicTimeON) call CheckCyclicMonths(NextTime, PreviousTime = PreviousTime)

                if(PreviousTime .le. Now .and. NextTime .ge. Now) then
                    FoundSecondInstant  = .true.
                    exit
                end if

                PreviousTime            = NextTime

                if(Me%ProfileTimeSerie%NextInstant .gt. Me%ProfileTimeSerie%NumberOfInstants &
                   .and. .not. Me%ProfileTimeSerie%CyclicTimeON)then
                    write(*,*)
                    write(*,*)'Could not read solution from Profile Time Serie file'
                    write(*,*)'Could not find second instant in file'
                    write(*,*)'Matrix name: '//trim(Me%PropertyID%Name)
                    stop      'ConstructProfileTSDefault - ModuleFillMatrix - ERR270'
                end if

            end do

            Me%ProfileTimeSerie%PreviousTime     = PreviousTime
            Me%ProfileTimeSerie%NextTime         = NextTime

        endif

        call UnGetGeometry (Me%ObjGeometry, SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                               &
            stop 'ConstructProfileTSDefault - ModuleFillMatrix - ERR280'

        allocate(Me%ProfileTimeSerie%PreviousField3D (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%ProfileTimeSerie%NextField3D     (ILB:IUB, JLB:JUB, KLB:KUB))

        Me%ProfileTimeSerie%PreviousField3D(:,:,:) = FillValueReal
        Me%ProfileTimeSerie%NextField3D    (:,:,:) = FillValueReal

        call ProfileTimeSerieField(PointsToFill3D, Me%ProfileTimeSerie%PreviousInstant, &
                                                   Me%ProfileTimeSerie%PreviousField3D)
        call ProfileTimeSerieField(PointsToFill3D, Me%ProfileTimeSerie%NextInstant,     &
                                                   Me%ProfileTimeSerie%NextField3D    )

        if (Me%PropertyID%IsAngle) then            
            call InterpolateAngle3DInTime (ActualTime       = Now,                                   &
                                           Size             = Me%WorkSize3D,                         &
                                           Time1            = Me%ProfileTimeSerie%PreviousTime,      &
                                           Matrix1          = Me%ProfileTimeSerie%PreviousField3D,   &
                                           Time2            = Me%ProfileTimeSerie%NextTime,          &
                                           Matrix2          = Me%ProfileTimeSerie%NextField3D,       &
                                           MatrixOut        = Me%Matrix3D,                           &
                                           PointsToFill3D   = PointsToFill3D)
        else
            call InterpolateMatrix3DInTime(ActualTime       = Now,                                   &
                                           Size             = Me%WorkSize3D,                         &
                                           Time1            = Me%ProfileTimeSerie%PreviousTime,      &
                                           Matrix1          = Me%ProfileTimeSerie%PreviousField3D,   &
                                           Time2            = Me%ProfileTimeSerie%NextTime,          &
                                           Matrix2          = Me%ProfileTimeSerie%NextField3D,       &
                                           MatrixOut        = Me%Matrix3D,                           &
                                           PointsToFill3D   = PointsToFill3D)
        endif
    end subroutine ConstructProfileTSDefault

    !--------------------------------------------------------------------------

    subroutine ConstructAnalyticProfile (ExtractType, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        integer                                     :: ExtractType
        integer, dimension(:, :, :), pointer        :: PointsToFill3D

        !Local----------------------------------------------------------------
        real, dimension(:,:,:), pointer             :: SZZ
        real                                        :: CellDepth, CoefA, CoefB
        character(len=StringLength)                 :: AuxString
        integer                                     :: STAT_CALL
        integer                                     :: iflag, ProfileType
        integer                                     :: i, j, k
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB

        !----------------------------------------------------------------------
        !Begin----------------------------------------------------------------

        ILB = Me%WorkSize3D%ILB
        IUB = Me%WorkSize3D%IUB

        JLB = Me%WorkSize3D%JLB
        JUB = Me%WorkSize3D%JUB

        KLB = Me%WorkSize3D%KLB
        KUB = Me%WorkSize3D%KUB

        !Gets Center of the cells
        call GetGeometryDistances(Me%ObjGeometry, SZZ = SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticProfile - ModuleFillMatrix - ERR01'

        call GetData(AuxString, Me%ObjEnterData,  iflag,                            &
                     SearchType     = ExtractType,                                  &
                     keyword        = 'PROFILE_TYPE',                               &
                     default        = "Linear",                                     &
                     ClientModule   = 'ModuleFillMatrix',                           &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructAnalyticProfile - ModuleFillMatrix - ERR03'


        select case (trim(adjustl(AuxString)))
            case ("Linear",   "LINEAR",   "linear")
                ProfileType = Linear
            case ("Exponential",     "EXPONENTIAL",     "exponential")
                ProfileType = Exponential
            case default
                write(*,*)'Invalid option for keyword INITIALIZATION_METHOD'
                stop 'ConstructAnalyticProfile - ModuleFillMatrix - ERR04'
        end select


        call GetData(CoefA, Me%ObjEnterData,  iflag,                                &
                     SearchType     = ExtractType,                                  &
                     keyword        = 'CoefA',                                      &
                     default        = 0.1,                                          &
                     ClientModule   = 'ModuleFillMatrix',                           &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructAnalyticProfile - ModuleFillMatrix - ERR05'

        call GetData(CoefB, Me%ObjEnterData,  iflag,                                &
                     SearchType     = ExtractType,                                  &
                     keyword        = 'CoefB',                                      &
                     default        = 4500.,                                        &
                     ClientModule   = 'ModuleFillMatrix',                           &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructAnalyticProfile - ModuleFillMatrix - ERR06'
        

        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB

            if (PointsToFill3D(i, j, k) == WaterPoint) then

                CellDepth          = (SZZ(i, j, k) + SZZ(i, j, k - 1)) / 2 - SZZ(i, j, KUB) 

                if      (ProfileType == Linear) then
                    !Linear profile
                    Me%Matrix3D(i,j,k) = Me%DefaultValue(1) + CoefA * CellDepth / CoefB
                    !Me%Matrix3D(i,j,k) = Me%DefaultValue + CoefA * CellDepth / CoefB
                else if (ProfileType == Exponential) then
                    !Exponential profile
                    Me%Matrix3D(i,j,k) = Me%DefaultValue(1) - CoefA * exp(- CellDepth / CoefB)
                    !Me%Matrix3D(i,j,k) = Me%DefaultValue - CoefA * exp(- CellDepth / CoefB)
                endif

            endif

        enddo
        enddo
        enddo


        call UnGetGeometry (Me%ObjGeometry, SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticProfile - ModuleFillMatrix - ERR14'


    end subroutine ConstructAnalyticProfile

    !--------------------------------------------------------------------------



    subroutine ConstructSpaceBox (ExtractType, PointsToFill2D, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        integer                                         :: ExtractType
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D

        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: iflag
        integer                                         :: i, j, k
        integer, dimension (:, :   ), pointer           :: Boxes2D
        integer, dimension (:, :, :), pointer           :: Boxes3D
        integer                                         :: BoxesNumber


        if (Me%TypeZUV /= TypeZ_) then
            write(*,*)'Cannot initialize U or V matrixes with boxes'
            stop 'ConstructSpaceBox - ModuleFillMatrix - ERR01'
        endif

        !Gets name of the Box definition file
        call GetData(Me%Boxes%FileName,                                          &
                     Me%ObjEnterData, iflag,                                     &
                     SearchType   = ExtractType,                                 &
                     keyword      = 'FILENAME',                                  &
                     ClientModule = 'ModuleFillMatrix',                          &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceBox - ModuleFillMatrix - ERR02'
        if (iflag==0)then
            write(*,*)'Box File Name not given'
            stop 'ConstructSpaceBox - ModuleFillMatrix - ERR03'
        end if

        !Starts BoxDif / Gets Boxes and number of boxes
        if (Me%Dim == Dim2D) then
            
            call StartBoxDif(BoxDifID           = Me%Boxes%ObjBoxDif,               &
                             TimeID             = Me%ObjTime,                       &
                             HorizontalGridID   = Me%ObjHorizontalGrid,             &
                             BoxesFilePath      = Me%Boxes%FileName,                &
                             WaterPoints2D      = PointsToFill2D,                   &
                             STAT               = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceBox - ModuleFillMatrix - ERR04'

            call GetBoxes(Me%Boxes%ObjBoxDif, Boxes2D = Boxes2D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceBox - ModuleFillMatrix - ERR05'

            call GetNumberOfBoxes(Me%Boxes%ObjBoxDif, NumberOfBoxes2D = BoxesNumber, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceBox - ModuleFillMatrix - ERR06'

        else

            call StartBoxDif(BoxDifID           = Me%Boxes%ObjBoxDif,               &
                             TimeID             = Me%ObjTime,                       &
                             HorizontalGridID   = Me%ObjHorizontalGrid,             &
                             BoxesFilePath      = Me%Boxes%FileName,                &
                             WaterPoints3D      = PointsToFill3D,                   &
                             Size3D             = Me%Size3D,                        &
                             WorkSize3D         = Me%WorkSize3D,                    &
                             STAT               = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceBox - ModuleFillMatrix - ERR07'

            call GetBoxes(Me%Boxes%ObjBoxDif, Boxes3D = Boxes3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceBox - ModuleFillMatrix - ERR08'

            call GetNumberOfBoxes(Me%Boxes%ObjBoxDif, NumberOfBoxes3D = BoxesNumber, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceBox - ModuleFillMatrix - ERR09'

        endif
        
        if(Me%OverrideValueKeywordON)then
            BoxesNumber = BoxesNumber + 1  !to account for box 0 
        endif

        !Gets boxes Values
        allocate (Me%Boxes%Values(BoxesNumber))

        call GetData(Me%Boxes%Values,                                            &
                     Me%ObjEnterData , iflag,                                    &
                     SearchType   = ExtractType,                                 &
                     keyword      = 'BOXES_VALUES',                              &
                     ClientModule = 'ModuleFillMatrix',                          &
                     STAT         = STAT_CALL)                                      

        if       (STAT_CALL .EQ. SIZE_ERR_)  then
            write(*,*) 'Incorrect number of boxes'
            stop 'ConstructSpaceBox - ModuleFillMatrix - ERR11'
        else if ((STAT_CALL .NE. SIZE_ERR_) .AND.  (STAT_CALL .NE. SUCCESS_)) then
            stop 'ConstructSpaceBox - ModuleFillMatrix - ERR12'
        end if 
                  
        if (iflag==0) then
        
            if(Me%OverrideValueKeywordON)then
                
                call GetData(Me%Boxes%Values, Me%ObjEnterData,  iflag,             &
                             SearchType     = ExtractType,                         &
                             keyword        = trim(Me%OverrideValueKeyword),       &
                             ClientModule   = 'ModuleFillMatrix',                  &
                             STAT           = STAT_CALL)
                
                if       (STAT_CALL .EQ. SIZE_ERR_)  then
                    write(*,*) 'Incorrect number of boxes for property '
                    write(*,*)trim(Me%PropertyID%Name)//', '//trim(Me%OverrideValueKeyword) 
                    stop 'ConstructSpaceBox - ModuleFillMatrix - ERR12b'
                else if ((STAT_CALL .NE. SIZE_ERR_) .AND.  (STAT_CALL .NE. SUCCESS_)) then
                    stop 'ConstructSpaceBox - ModuleFillMatrix - ERR12c'
                end if           
                
                if(iflag == 0)then
                    write(*,*)'Please define override keyword '//trim(Me%PropertyID%Name)
                    write(*,*)'to give a boxes values for property '//trim(Me%PropertyID%Name)
                    stop 'ConstructSpaceBox - ModuleFillMatrix - ERR12d'
                end if
                
                if (Me%Dim == Dim2D) then
                    Boxes2D = Boxes2D + 1
                else
                    Boxes3D = Boxes3D + 1
                endif
            
            else

                write(*,*) 'Boxes Values not given for property '//trim(Me%PropertyID%Name)           
                stop       'ConstructSpaceBox - ModuleFillMatrix - ERR13'

            end if
        
        end if

        !Fills Matrix
        if (Me%Dim == Dim2D) then
            do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                if (Boxes2D(i,j) > 0) then
                    Me%Matrix2D (i, j) = Me%Boxes%Values(Boxes2D(i, j))
                end if
            end do
            end do

            call UngetBoxDif(Me%Boxes%ObjBoxDif, Boxes2D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceBox - ModuleFillMatrix - ERR14'

        else

            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
                if (Boxes3D(i,j,k) > 0) then
                    Me%Matrix3D (i, j, k) = Me%Boxes%Values(Boxes3D(i, j, k))
                end if
            end do
            end do
            end do

            call UngetBoxDif(Me%Boxes%ObjBoxDif, Boxes3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceBox - ModuleFillMatrix - ERR15'

        endif

        deallocate (Me%Boxes%Values)

        call KillBoxDif(Me%Boxes%ObjBoxDif, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceBox - ModuleFillMatrix - ERR16'

    end subroutine ConstructSpaceBox

    !--------------------------------------------------------------------------

    subroutine ConstructSpaceASCIIFile (ExtractType, PointsToFill2D, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        integer                                         :: ExtractType
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D

        !Local----------------------------------------------------------------
        integer                                     :: STAT_CALL, i, j, k
        integer                                     :: iflag, TypeZUV
        real, dimension(:, :),    pointer           :: GridData2D
        real, dimension(:, :, :), pointer           :: GridData3D
        integer                                     :: file
        character(len = PathLength), dimension(2)   :: Filename         = null_str
        type(T_ASCIIFile), pointer                  :: NewASCIIFile, CurrentASCIIFile
        integer                                     :: nASCIIFiles     = 1
        logical                                     :: exist
        !Begin----------------------------------------------------------------
        
        Me%nASCIIFiles = 0
        
        !Gets the name of the data file
        !If vectorial use different keywords for x and y since GetData would not allow spaces
        !or would interpret as a new filename        
        if (Me%VectorialProp) then         
            nASCIIFiles = 2
            call GetData(FileName(1),                                                &
                         Me%ObjEnterData , iflag,                                    &
                         SearchType   = ExtractType,                                 &
                         keyword      = 'FILENAME_X',                                &
                         ClientModule = 'ModuleFillMatrix',                          &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR00a'              
            if (iflag==0)then
                write(*,*) 'ASCII FILENAME_X keyword not given not given for vectorial property '//trim(Me%PropertyID%Name)
                stop       'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR00b'                
            endif
            inquire (file=trim(FileName(1)), exist = exist)
            if (.not. exist) then
                write(*,*)'Could not find file '//trim(FileName(1))
                stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR00c'
            endif      
            
            call GetData(FileName(2),                                                &
                         Me%ObjEnterData , iflag,                                    &
                         SearchType   = ExtractType,                                 &
                         keyword      = 'FILENAME_Y',                                &
                         ClientModule = 'ModuleFillMatrix',                          &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR00d'     
            if (iflag==0)then
                write(*,*) 'ASCII FILENAME_Y keyword not given not given for vectorial property '//trim(Me%PropertyID%Name)
                stop       'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR00e'                
            endif
            inquire (file=trim(FileName(2)), exist = exist)
            if (.not. exist) then
                write(*,*)'Could not find file '//trim(FileName(2))
                stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR00f'
            endif      
            
        else
            
            call GetData(FileName(1),                                                &
                         Me%ObjEnterData , iflag,                                    &
                         SearchType   = ExtractType,                                 &
                         keyword      = 'FILENAME',                                  &
                         ClientModule = 'ModuleFillMatrix',                          &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR01a'
        
            if (iflag==0)then

               if(Me%OverrideValueKeywordON)then
                
                    call GetData(FileName(1),                                          &
                                 Me%ObjEnterData , iflag,                              &
                                 SearchType   = ExtractType,                           &
                                 keyword      = trim(Me%OverrideValueKeyword),         &
                                 ClientModule = 'ModuleFillMatrix',                    &
                                 STAT         = STAT_CALL)                                      
                
                    if (STAT_CALL .NE. SUCCESS_)  then
                        write(*,*)trim(Me%PropertyID%Name)//', '//trim(Me%OverrideValueKeyword) 
                        stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR01b'
                    end if           
                
                    if(iflag == 0)then
                        write(*,*)'Please define the ASCII file in the override keyword '//trim(Me%PropertyID%Name)
                        write(*,*)'to give values for property '//trim(Me%PropertyID%Name)
                        stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR01c'
                    end if
                           
                else

                    write(*,*) 'ASCII File Name not given not given for property '//trim(Me%PropertyID%Name)           
                    stop       'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR01d'

                end if     
            else
                inquire (file=trim(FileName(1)), exist = exist)
                if (.not. exist) then
                    write(*,*)'Could not find file '//trim(FileName(1))
                    stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR01e'
                endif                   
            endif
        end if      
        
        do i = 1, nASCIIFiles

            nullify  (NewASCIIFile)
            allocate (NewASCIIFile)
            nullify  (NewASCIIFile%Next)
            nullify  (NewASCIIFile%Prev)

            !in vectorial prop always exist two names
            NewASCIIFile%Filename = Filename(i)  
            
            call InsertASCIIFileToList(NewASCIIFile)
                
        enddo
            
        CurrentASCIIFile => Me%FirstASCIIFile
        file = 1
        do while (associated(CurrentASCIIFile))        
        
            !Associate Matrix2D and Me%Matrix3D to the input field ones
            if (Me%VectorialProp .or. Me%RotateAngleToGrid) call AssociateMatrixes(file)                        
            
            if (Me%Dim == Dim2D) then

                
                call ConstructGridData(CurrentASCIIFile%GridDataID, Me%ObjHorizontalGrid,        &
                                       FileName     = CurrentASCIIFile%FileName,                 &
                                       DefaultValue = Me%DefaultValue(file),                     &
                                       STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR03'

                call GetGridDataType(CurrentASCIIFile%GridDataID, TypeZUV, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR04'

                if(TypeZUV .ne. Me%TypeZUV)then
                    write(*,*)'Inconsistency found in type ZUV'
                    write(*,*)'Grid data: '//trim(CurrentASCIIFile%FileName)
                    stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR05'
                end if

                call GetGridData      (CurrentASCIIFile%GridDataID, GridData2D, STAT = STAT_CALL)  
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR06'
                                
                !Copies data
                where (PointsToFill2D == WaterPoint) Me%Matrix2D = GridData2D
                
                call UnGetGridData    (CurrentASCIIFile%GridDataID, GridData2D, STAT = STAT_CALL)  
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR07'

            else

                call ConstructGridData(CurrentASCIIFile%GridDataID, Me%ObjHorizontalGrid,        &
                                       FileName     = CurrentASCIIFile%FileName,                 &
                                       KLB          = Me%Worksize3D%KLB,                     &
                                       KUB          = Me%Worksize3D%KUB,                     &
                                       DefaultValue = Me%DefaultValue(i),                    &
                                       STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR08'


                call GetGridDataType(CurrentASCIIFile%GridDataID, TypeZUV, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR09'
            
                if(TypeZUV .ne. Me%TypeZUV)then
                    write(*,*)'Inconsistency found in type ZUV'
                    write(*,*)'Grid data: '//trim(CurrentASCIIFile%FileName)
                    stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR10'
                end if

                call GetGridData      (CurrentASCIIFile%GridDataID, GridData3D, STAT = STAT_CALL)  
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR11'
                
                !Copies data
                do k = Me%Worksize3D%KLB, Me%Worksize3D%KUB
                do j = Me%Worksize3D%JLB, Me%Worksize3D%JUB
                do i = Me%Worksize3D%ILB, Me%Worksize3D%IUB
                    if (PointsToFill3D(i, j, k) == WaterPoint) Me%Matrix3D(i, j, k) = GridData3D(i, j, k)
                enddo
                enddo
                enddo
                
                call UnGetGridData    (CurrentASCIIFile%GridDataID, GridData3D, STAT = STAT_CALL)  
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR12'

            endif   
            
            call KillGridData (CurrentASCIIFile%GridDataID, STAT = STAT_CALL)  
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR13'            
                   
            
            file = file + 1
            CurrentASCIIFile => CurrentASCIIFile%Next
        enddo
        


    end subroutine ConstructSpaceASCIIFile

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ConstructSponge (ExtractType, PointsToFill2D, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        integer                                         :: ExtractType
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D

        !Local----------------------------------------------------------------
        real, dimension(:),       pointer           :: AuxT
        integer, dimension(4,4)                     :: dij
        integer, dimension(4)                       :: AuxI
        real                                        :: Aux, default
        integer                                     :: STAT_CALL, i, j, k, l, sp
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: iflag
        real                                        :: DefaultOutVal     

        !Begin----------------------------------------------------------------



        if (Me%Dim == Dim2D) then
        
            JLB = Me%Worksize2D%JLB
            JUB = Me%Worksize2D%JUB
            ILB = Me%Worksize2D%ILB
            IUB = Me%Worksize2D%IUB
            KLB = 1
            KUB = 1
                        
        else
        
            JLB = Me%Worksize3D%JLB
            JUB = Me%Worksize3D%JUB
            ILB = Me%Worksize3D%ILB
            IUB = Me%Worksize3D%IUB
            KLB = Me%Worksize3D%KLB
            KUB = Me%Worksize3D%KUB
            
        endif
        
        DefaultOutVal = 1e5
        if      (present(PointsToFill2D)) then
            call CheckViscSpongeOption(ExtractType = ExtractType, PointsToFill2D = PointsToFill2D, DefaultOutVal = DefaultOutVal)
        else if (present(PointsToFill3D)) then
            call CheckViscSpongeOption(ExtractType = ExtractType, PointsToFill3D = PointsToFill3D, DefaultOutVal = DefaultOutVal)
        endif
        
        call CheckReferenceSolDT(ExtractType = ExtractType, DefaultOutVal = DefaultOutVal)

        !Gets the sponge value in the model open boundary
        call GetData(Me%Sponge%OutValue,                                                &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'SPONGE_OUT',                                       &
                     Default      = DefaultOutVal,                                      &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSponge - ModuleFillMatrix - ERR20'
        
        

        !Gets the number of sponge cells
        call GetData(Me%Sponge%Cells,                                                   &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'SPONGE_CELLS',                                     &
                     Default      = 10,                                                  &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSponge - ModuleFillMatrix - ERR30'
        
        if (Me%Sponge%Cells > IUB - ILB + 1 .and. Me%Sponge%Cells > JUB - JLB + 1) then
            stop 'ConstructSponge - ModuleFillMatrix - ERR36'
        endif        

        !Gets the nsponge evolution
        call GetData(Me%Sponge%Evolution,                                               &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'SPONGE_EVOLUTION',                                 &
                     Default      = sponge_exp_,                                        &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSponge - ModuleFillMatrix - ERR40'

        if      (Me%Sponge%Evolution /= sponge_exp_             .and.                   &
                 Me%Sponge%Evolution /= sponge_linear_          .and.                   &
                 Me%Sponge%Evolution /= sponge_wave_stress_dump_) then        

            write(*,*) 'Sponge evolution can only be linear or exponential'
            stop       'ConstructSponge - ModuleFillMatrix - ERR50'
        
        endif
        
        if (Me%DefaultValue(1) < Me%Sponge%OutValue) then
        !if (Me%DefaultValue < Me%Sponge%OutValue) then
        
            Me%Sponge%Growing = .true.
            
        else
        
            Me%Sponge%Growing = .false.
        
        endif
        
        if (Me%TypeZUV == TypeU_ .or. Me%TypeZUV == TypeV_) then
            allocate(AuxT(Me%Sponge%Cells + 1))
        else
            allocate(AuxT(Me%Sponge%Cells))
        
        endif
        
        default = Me%DefaultValue(1)
        !default = Me%DefaultValue
        
        if      (Me%Sponge%Evolution == sponge_exp_) then
        
            do sp = 1, Me%Sponge%Cells

                AuxT(sp) = log(Me%Sponge%OutValue) * real(Me%Sponge%Cells - sp) /real(Me%Sponge%Cells - 1) + &
                           log(default)  * real(sp - 1)               /real(Me%Sponge%Cells - 1)
                     
                AuxT(sp) = exp(AuxT(sp))
                
            enddo
        
        elseif (Me%Sponge%Evolution == sponge_linear_) then        
        
            do sp = 1, Me%Sponge%Cells

                AuxT(sp) = Me%Sponge%OutValue * real(Me%Sponge%Cells - sp) /real(Me%Sponge%Cells - 1) + &
                           default    * real(sp - 1)               /real(Me%Sponge%Cells - 1)
                     
            enddo

        elseif (Me%Sponge%Evolution == sponge_wave_stress_dump_) then       
        
                AuxT(Me%Sponge%Cells) = 1.
        
                do sp = Me%Sponge%Cells,2,-1

                    AuxT(sp-1) = max(1e-10,AuxT(sp) * 0.5)
                    
                enddo

        endif        
        
        if (Me%TypeZUV == TypeU_ .or. Me%TypeZUV == TypeV_) then
            AuxT(Me%Sponge%Cells+1) = AuxT(Me%Sponge%Cells)
            Me%Sponge%Cells = Me%Sponge%Cells + 1
        endif
        
        
        dij(:,:) = 0
        
        Me%Sponge%OpenBordersON  = GetDDecompOpenBorders(Me%ObjHorizontalGrid, &
                                                                     STAT = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSponge - ModuleFillMatrix - ERR60'

dsp:    do sp = 1, Me%Sponge%Cells
            
            AuxI(1) = min(ILB + sp - 1, IUB)
            AuxI(2) = max(IUB - sp + 1, ILB)
            AuxI(3) = min(JLB + sp - 1, JUB)
            AuxI(4) = max(JUB - sp + 1, JLB)
            
            
            
            dij(1,1) = max (1 - sp, 1 - (IUB - ILB))
            dij(2,2) = min (sp - 1, (IUB - ILB) - 1)
            dij(3,3) = max (1 - sp, 1 - (JUB - JLB))
            dij(4,4) = min (sp - 1, (JUB - JLB) - 1)
            
            
            Aux     = AuxT(sp)

            !1 - South; 2 - North; 3 - West; 4 - East        
            !Me%Sponge%OpenBordersON(:)

dk:         do k =  KLB,  KUB

               !Southern and Northen boundary
    dl:         do l=1,2
    
                    !1 - South; 2 - North; 3 - West; 4 - East        
                    if (.not.Me%Sponge%OpenBordersON(l)) cycle 
                   
                    i = AuxI(l)

    dj:             do j = JLB, JUB
                                
                        if (CheckSponge(PointsToFill2D, PointsToFill3D, sp, dij(l,1), dij(l,2), 0, 0, i, j, k)) then
                            call FillSponge(PointsToFill2D, PointsToFill3D, Aux, i, j, k)
                        endif
                
                    enddo dj
                 
                enddo dl
                
                !Western and Eastern boundary
dl2:            do l=3,4

                    !1 - South; 2 - North; 3 - West; 4 - East        
                    if (.not.Me%Sponge%OpenBordersON(l)) cycle 
                    
                    j = AuxI(l)

di:                 do i = ILB, IUB
                                                
                        if (CheckSponge(PointsToFill2D, PointsToFill3D, sp, 0, 0, dij(l,3), dij(l,4), i, j, k)) then
                            call FillSponge(PointsToFill2D, PointsToFill3D, Aux, i, j, k)
                        endif

                    enddo di
                        
                enddo dl2
                
            enddo dk
            
        enddo dsp

        deallocate(AuxT)

    end subroutine ConstructSponge
    
    !--------------------------------------------------------------------------
    
    subroutine CheckViscSpongeOption (ExtractType, PointsToFill2D, PointsToFill3D, DefaultOutVal)

        !Arguments-------------------------------------------------------------
        integer                                         :: ExtractType
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D
        real                                            :: DefaultOutVal

        !Local----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        real,   dimension(:,:), pointer                 :: DUX, DVY
        real                                            :: DX, DT, ViscOut, ViscIn 
        logical                                         :: ViscTurbSponge

        !Begin----------------------------------------------------------------

    
    
        !Check ifa default viscosity want to be assumed
        call GetData(ViscTurbSponge,                                                    &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'VISC_TURB_SPONGE',                                 &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) then
            stop 'CheckViscSpongeOption - ModuleFillMatrix - ERR10'  
        endif            
        
        if (ViscTurbSponge) then
            
            call GetComputeTimeStep(Me%ObjTime, DT = DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'CheckViscSpongeOption - ModuleFillMatrix - ERR20'  
            endif
            
            call GetHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,             &
                                   DUX              = DUX,                              &
                                   DVY              = DVY,                              &
                                   STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'CheckViscSpongeOption - ModuleFillMatrix - ERR30'  
            endif             
           
            DX = max(DUX(1,1), DVY(1,1))
            
            !dx / 100
            ViscIn = DX/100.
            
            ! 0.1 * dx^2/dt -> stability limit 0.5 * dx^2/dt
            ViscOut       = 0.1*DX**2/DT
            DefaultOutVal = ViscOut
            
            if (present(PointsToFill2D)) then
                where (PointsToFill2D == WaterPoint) Me%Matrix2D = ViscIn
            endif                
            if (present(PointsToFill3D)) then
                where (PointsToFill3D == WaterPoint) Me%Matrix3D = ViscIn
            endif  
            
            call UnGetHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     Array            = DUX,                            &
                                     STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'CheckViscSpongeOption - ModuleFillMatrix - ERR40'  
            endif     
            
            call UnGetHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     Array            = DVY,                            &
                                     STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'CheckViscSpongeOption - ModuleFillMatrix - ERR50'  
            endif              
        endif
            
    end subroutine CheckViscSpongeOption
    
    !--------------------------------------------------------------------------    
    subroutine CheckReferenceSolDT (ExtractType, DefaultOutVal)

        !Arguments-------------------------------------------------------------
        integer                                         :: ExtractType
        real                                            :: DefaultOutVal

        !Local----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        logical                                         :: FILE_DT, exist
        integer, save                                   :: ObjHDF5
        integer                                         :: HDF5_READ        
        real,   dimension(:), pointer                   :: TimeVector
        type (T_Time)                                   :: StartDate, EndDate

        !Begin----------------------------------------------------------------

        !Check ifa default viscosity want to be assumed
        call GetData(FILE_DT,                                                           &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'FILE_DT',                                          &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) then
            stop 'CheckReferenceSolDT - ModuleFillMatrix - ERR10'  
        endif            
        
        if (FILE_DT) then
            
            if (Me%SpongeFILE_DT == null_str) then
                stop 'CheckReferenceSolDT - ModuleFillMatrix - ERR20'  
            endif 

            inquire (file=trim(Me%SpongeFILE_DT), exist = exist)
            if (.not. exist) then
                write(*,*)'Could not find file '//trim(Me%SpongeFILE_DT)
                stop 'CheckReferenceSolDT - ModuleFillMatrix - ERR30'
            endif  
            
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)            
            
            call ConstructHDF5 (ObjHDF5, trim(Me%SpongeFILE_DT), HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'CheckReferenceSolDT - ModuleFillMatrix - ERR40'  
            endif  

            call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'CheckReferenceSolDT - ModuleFillMatrix - ERR50'  
            endif  

            allocate(TimeVector(6))

            call HDF5ReadData   (HDF5ID         = ObjHDF5,                              &
                                 GroupName      = "/Time",                              &
                                 Name           = "Time",                               &
                                 Array1D        = TimeVector,                           &
                                 OutputNumber   = 1,                                    &
                                 STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'CheckReferenceSolDT - ModuleFillMatrix - ERR60'
            endif                

            call SetDate(StartDate, Year     = TimeVector(1), Month  = TimeVector(2),   &
                                    Day      = TimeVector(3), Hour   = TimeVector(4),   &
                                    Minute   = TimeVector(5), Second = TimeVector(6))
            

            call HDF5ReadData   (HDF5ID         = ObjHDF5,                              &
                                 GroupName      = "/Time",                              &
                                 Name           = "Time",                               &
                                 Array1D        = TimeVector,                           &
                                 OutputNumber   = 2,                                    &
                                 STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'CheckReferenceSolDT - ModuleFillMatrix - ERR70'
            endif        

            call SetDate(EndDate, Year     = TimeVector(1), Month  = TimeVector(2),     &
                                  Day      = TimeVector(3), Hour   = TimeVector(4),     &
                                  Minute   = TimeVector(5), Second = TimeVector(6))
                                     
            deallocate(TimeVector)
            
            
            DefaultOutVal = EndDate - StartDate
            
            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'CheckReferenceSolDT - ModuleFillMatrix - ERR80'  
            endif              
            
        endif            
            
    end subroutine CheckReferenceSolDT
    
    !--------------------------------------------------------------------------    

    
    subroutine FillSponge(PointsToFill2D, PointsToFill3D, Aux, i, j, k)
    
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D
        real                                            :: Aux
        integer                                         :: i, j, k

        !Local-----------------------------------------------------------------            

        !Begin-----------------------------------------------------------------    
    
i2:     if (Me%Dim == Dim2D) then

i3:         if (PointsToFill2D(i, j) == WaterPoint) then

ig:             if    (       Me%Sponge%Growing .and. Aux >  Me%Matrix2D(i, j)) then
                
                    Me%Matrix2D(i, j) = Aux
                
                elseif (.not. Me%Sponge%Growing .and. Aux <  Me%Matrix2D(i, j)) then ig
            
                    Me%Matrix2D(i, j) = Aux
                    
                endif ig
                
            endif i3
            
        else i2
                                    
       
i4:         if (PointsToFill3D(i, j, k) == WaterPoint) then

i5:             if (      Me%Sponge%Growing .and. Aux >  Me%Matrix3D(i, j, k)) then
                
                    Me%Matrix3D(i, j, k) = Aux
                
                elseif (.not. Me%Sponge%Growing .and. Aux <  Me%Matrix3D(i, j, k)) then i5
            
                    Me%Matrix3D(i, j, k) = Aux
                
                endif i5
            
            endif i4
                
        endif  i2      
    
    end subroutine FillSponge

    !--------------------------------------------------------------------------    

    !--------------------------------------------------------------------------

    subroutine ConstructAnalyticWave (ExtractType)

        !Arguments-------------------------------------------------------------
        integer                             :: ExtractType

        !Local-----------------------------------------------------------------
        integer,   dimension(:,:), pointer  :: Aux2D
        real,      dimension(:,:), pointer  :: XX2D, YY2D, Bathymetry, DUX, DVY
        integer                             :: STAT_CALL
        integer                             :: iflag, i, j
        integer                             :: ObjHorizontalGridAux
        type (T_Polygon), pointer           :: InteriorCellsArea
        type (T_PointF ), pointer           :: Point, Point0
        type (T_Segment), pointer           :: Segment   
        real(8)                             :: PointX, PointY         
        real(8)                             :: DifX, DifY         
        real(8)                             :: PointDirection, DifDirection    
        real                                :: East, West, North, South
        real(8)                             :: CircleCenterX, CircleCenterY, Radius
        real(8)                             :: WaveDirection, Xorig, Yorig
        real(8)                             :: RunPeriod
        real(8)                             :: Omega, k, Aux, A, H, MinDX, T, L
        character(len = PathLength)         :: BathymetryFile
        real(8)                             :: Dif, AuxL1, AuxL
        integer                             :: ObjBathymetry
        integer                             :: EnterNcells, n
        real                                :: WaveHeight
        !Begin----------------------------------------------------------------

        !Gets the wave amplitude (m)
        call GetData(WaveHeight,                                                        &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'WAVE_HEIGHT',                                      &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (iflag == 0           ) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR10'
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR20'
        
        Me%AnalyticWave%Amplitude = WaveHeight / 2. 
        

        !Gets the wave direction (degrees in meteorological convention)
        call GetData(Me%AnalyticWave%Direction,                                         &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'WAVE_DIRECTION',                                   &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (iflag == 0           ) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR30'
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR40'
        
        !from nautical reference in degrees to cartesian in radians
        Me%AnalyticWave%Direction = (270 - Me%AnalyticWave%Direction) * Pi/180.


        !Gets the wave period (s)
        call GetData(Me%AnalyticWave%Period,                                            &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'WAVE_PERIOD',                                      &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (iflag == 0           ) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR50'
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR60'
        
        RunPeriod   = Me%EndTime - Me%BeginTime
        
        !if(real(int(RunPeriod/Me%AnalyticWave%Period)) /= real(RunPeriod/Me%AnalyticWave%Period)) then
        !    write(*,*) 'Wave Period needs to be a multiple of the run period'
        !    stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR70'
        !endif
        
        if (Me%AnalyticWave%Period <= 0.) then
            write(*,*) 'wave period needs to be greater than 0'
            stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR80'
        endif
        

        !Gets the wave type 
        !wave types
!            SineWaveSeaLevel_          = 1
!            CnoidalWaveSeaLevel_       = 2
!            SolitartyWaveSeaLevel_     = 3
!            SineWaveVelX_              = 4
!            CnoidalWaveVelX_           = 5
!            SolitartyWaveVelX_         = 6
!            SineWaveVelY_              = 7
!            CnoidalWaveVelY_           = 8
!            SolitartyWaveVelY_         = 9
        call GetData(Me%AnalyticWave%WaveType,                                          &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'WAVE_TYPE',                                        &
                     default      = SineWaveSeaLevel_,                                  &                   
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (iflag == 0           ) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR90'
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR100'
        
        if (Me%AnalyticWave%WaveType < 0 .and. Me%AnalyticWave%WaveType > 9) then
            write(*,*) 'Not valid WAVE_TYPE option =', Me%AnalyticWave%WaveType
            stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR110'
        endif
        
        !Gets the wave reference water level (m)
        call GetData(Me%AnalyticWave%AverageValue,                                      &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'AVERAGE_VALUE',                                    &
                     default      = 0.,                                                 &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR120'
        
        !Gets the depth use to compute the imposed wave (m)
        call GetData(Me%AnalyticWave%DepthValue,                                        &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'DEPTH_VALUE',                                      &
                     default      = null_real,                                          &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR123'

        !Gets the coef value to adjust the wave celerity to be imposed in the boundary
        call GetData(Me%AnalyticWave%CoefValue,                                         &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'COEF_VALUE',                                       &
                     default      = 1.,                                                 &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR125' 
        
        
        !Gets if the wave is to imposed gradually or not along a period (by default wave period).
        call GetData(Me%AnalyticWave%SlowStartON,                                       &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'SLOWSTART',                                        &
                     default      = .true.,                                             &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR126'
        
        if (Me%AnalyticWave%SlowStartON) then
        
            call GetData(Me%AnalyticWave%SlowStartPeriod,                               &
                         Me%ObjEnterData , iflag,                                       &
                         SearchType   = ExtractType,                                    &
                         keyword      = 'SLOWSTART_PERIOD',                             &
                         default      = Me%AnalyticWave%Period,                         &
                         ClientModule = 'ModuleFillMatrix',                             &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR126a'         
            
        endif            
        
        !Gets the dif coef.
        call GetData(Me%AnalyticWave%Dif,                                               &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'DIF_VALUE',                                        &
                     default      = 1e-5,                                               &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR127'        
        
        
        !Gets the file name of the Bathymetry
        call ReadFileName('IN_BATIM', BathymetryFile, "Bathymetry File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR130'
        
        ObjBathymetry = 0
        
        !Horizontal Grid Data - Water Column (Bathymetry)
        call ConstructGridData      (GridDataID       = ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,        &
                                     TimeID           = Me%ObjTime,                  &
                                     FileName         = BathymetryFile,              &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR140'
        

        !Gets Bathymetry
        call GetGridData(ObjBathymetry, Bathymetry, STAT = STAT_CALL)
                           
        allocate(Me%AnalyticWave%AmpAux  (Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))
        allocate(Me%AnalyticWave%Celerity(Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))
                          
               
        A = Me%AnalyticWave%Amplitude
        T = Me%AnalyticWave%Period    
                
        if (Me%AnalyticWave%WaveType == SineWaveSeaLevel_  .or.                         &
            Me%AnalyticWave%WaveType == SineWaveVelX_      .or.                         &
            Me%AnalyticWave%WaveType == SineWaveVelY_) then
            
            H = Me%AnalyticWave%DepthValue + Me%AnalyticWave%AverageValue
    
                            
            do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB        

                if (Me%AnalyticWave%DepthValue < HalfFillValueReal) then
                    H = max(Bathymetry(i, j) + Me%AnalyticWave%AverageValue, 0.1)
                endif                    
                
                Me%AnalyticWave%AmpAux(i, j) = A        
            
                L = WaveLengthHuntsApproximation(real(T), real(H))
                
                AuxL = sqrt(Gravity*H) * T
                Dif  = - FillValueReal 
                
                do while (Dif > Me%AnalyticWave%Dif)
                   AuxL1 = Gravity / (2*Pi)* T**2. * tanh(2*Pi*H/AuxL)
                   Dif   = abs(AuxL1 - AuxL)
                   AuxL  = AuxL1
                enddo
                
                L = AuxL
                
                if (L <= 0) then
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR150'
                endif 
                
                Omega = 2. * Pi / T
                k     = 2. * Pi / L
                
                
                Me%AnalyticWave%Celerity(i,j) = Me%AnalyticWave%CoefValue * L / T
                
                !write(*,*) 'Wave celerity - ', Me%AnalyticWave%Celerity(i, j)
                
                
                Aux = L / H
                
                if (Me%AnalyticWave%WaveType == SineWaveVelX_      .or.                     &
                    Me%AnalyticWave%WaveType == SineWaveVelY_) then     
                    
                    if     ( Aux >  20.               ) then
                        
                        Me%AnalyticWave%AmpAux(i, j) = Omega * A / k / H
                        
                    elseif  (Aux >  2. .and. Aux <=20.) then 
                    
                        Me%AnalyticWave%AmpAux(i, j) = Omega * A / k / H      
                    
                    elseif  (Aux <= 2.                ) then
                    
                        Me%AnalyticWave%AmpAux(i, j) = Omega * A / k / H * (1. - exp(-k*H))                
                    
                    endif
                endif            
                    
            enddo
            enddo
                
        endif
        
        !Kills Bathymetry
        call KillGridData(ObjBathymetry, STAT = STAT_CALL)
        
        call GetHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,             & 
                               XX2D_Z           = XX2D,                             &
                               YY2D_Z           = YY2D,                             &
                               DUX              = DUX,                              &
                               DVY              = DVY,                              &
                               STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR160'
        
        
        
        if (GetDDecompON(Me%ObjHorizontalGrid)) then
        
            ObjHorizontalGridAux = 0
    
            !Horizontal Grid
            call ConstructHorizontalGrid(HorizontalGridID = ObjHorizontalGridAux,       &
                                         DataFile         = BathymetryFile,             &
                                         STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR170'
        
        else
        
            ObjHorizontalGridAux = Me%ObjHorizontalGrid
        
        endif
        
        call GetGridOutBorderCartLimits(HorizontalGridID  = ObjHorizontalGridAux,       & 
                                        West              = West,                       &
                                        East              = East,                       &
                                        South             = South,                      &
                                        North             = North,                      &
                                        STAT              = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR180'
        
        call GetGridBorderCartPolygon  (HorizontalGridID  = ObjHorizontalGridAux,       & 
                                        Polygon           = InteriorCellsArea,          &
                                        STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR190'
        
        allocate(Me%AnalyticWave%X2D        (Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))
        allocate(Me%AnalyticWave%CellType   (Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))
        allocate(Me%AnalyticWave%TlagMissing(Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))        
        allocate(Aux2D                      (Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))
           

        Me%AnalyticWave%X2D         (:,:) = null_real
        Me%AnalyticWave%CellType    (:,:) = null_int
        Aux2D                       (:,:) = null_int 
        Me%AnalyticWave%TlagMissing (:,:) = 0
        
        
        CircleCenterX = (West  + East )/2.
        CircleCenterY = (South + North)/2.        
        
        Radius = 100 * sqrt(((East-West)/2.)**2.+((North-South)/2.)**2.) 
        
        WaveDirection = Me%AnalyticWave%Direction       
        
        !point in the circle first reach by the incoming wave 
        Xorig = CircleCenterX + Radius * cos(WaveDirection + Pi)
        Yorig = CircleCenterY + Radius * sin(WaveDirection + Pi)
        
        allocate(Segment)
        allocate(Point)
        allocate(Point0)        
        
        Point0%X    = Xorig
        Point0%Y    = Yorig
        
        Segment%StartAt     => Point0
        
        MinDX       = - FillValueReal
                
        do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
        do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB        
        
            PointX                      = XX2D(i, j)                   
            PointY                      = YY2D(i, j)   

            DifX                        = PointX - Xorig
            DifY                        = PointY - Yorig        

            PointDirection              = atan2(DifY, DifX)
            DifDirection                = PointDirection - WaveDirection

            Me%AnalyticWave%X2D(i, j)   = sqrt(DifX**2.+DifY**2.)*cos(DifDirection)
            
            Point%X                     = PointX
            Point%Y                     = PointY
            
            Segment%EndAt               => Point
            
!            EnteringWaveCell_ 
!            LeavingWaveCell_  
!            InteriorWaveCell_ 
!            ExteriorWaveCell_

!            if (IsPointInsidePolygon(Point, InteriorCellsArea)) then
!                !interior cell
!                Aux2D(i, j) = InteriorWaveCell_
!            !boundary cell
!            else
!                if (Intersect2D_SegPoly(Segment= Segment, Polygon= InteriorCellsArea)) then
!                    !leaving wave cell
!                    Aux2D(i, j) = LeavingWaveCell_
!                else
!                    !entering wave cell
!                    Aux2D(i, j) = EnteringWaveCell_
!                endif
!            endif
!            
        enddo
        enddo
!        
!        !SW corner
!        i = Me%WorkSize2D%ILB
!        j = Me%WorkSize2D%JLB
!        
!        Aux2D(i, j) = ExteriorWaveCell_
!        
!        !NW corner
!        i = Me%WorkSize2D%IUB
!        j = Me%WorkSize2D%JLB
!        
!        Aux2D(i, j) = ExteriorWaveCell_
!        
!        !SE corner
!        i = Me%WorkSize2D%ILB
!        j = Me%WorkSize2D%JUB
!        
!        Aux2D(i, j) = ExteriorWaveCell_
!        
!        !NE corner
!        i = Me%WorkSize2D%IUB
!        j = Me%WorkSize2D%JUB
!        
!        Aux2D(i, j) = ExteriorWaveCell_
!        
!        Me%AnalyticWave%CellType(:,:) = Aux2D(:,:)
        
        
        !leaving wave adjacent cells
!        do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
!        do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB   
!        
!                if (Aux2D(i, j) == ExteriorWaveCell_) cycle
!        
!                if (Aux2D(i-1,j  ) == LeavingWaveCell_ .or.                             &
!                    Aux2D(i+1,j  ) == LeavingWaveCell_ .or.                             &
!                    Aux2D(i  ,j-1) == LeavingWaveCell_ .or.                             &                    
!                    Aux2D(i  ,j+1) == LeavingWaveCell_ ) then
!                    Me%AnalyticWave%CellType(i, j) = LeavingWaveCell_
!                endif                    
!            
!        enddo
!        enddo        
        
        !entering wave adjacent cells - these cells have priority over Leaving cells
!        do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
!        do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB        
!        
!                if (Aux2D(i, j) == ExteriorWaveCell_) cycle        
!        
!                if (Aux2D(i-1,j  ) == EnteringWaveCell_ .or.                            &
!                    Aux2D(i+1,j  ) == EnteringWaveCell_ .or.                            &
!                    Aux2D(i  ,j-1) == EnteringWaveCell_ .or.                            &                    
!                    Aux2D(i  ,j+1) == EnteringWaveCell_ ) then
!                    Me%AnalyticWave%CellType(i, j) = EnteringWaveCell_
!                endif       
!                
!                !Temporary
!                if (Aux2D(i, j) == InteriorWaveCell_) then
!                    Me%AnalyticWave%CellType(i, j) = EnteringWaveCell_
!                endif          
!            
!        enddo
!        enddo     
        
        Me%AnalyticWave%CellType(:,:) = EnteringWaveCell_
        
        do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
        do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB   
        
            if (Me%AnalyticWave%CellType(i, j) == EnteringWaveCell_) then             
        
                if (Me%AnalyticWave%X2D(i, j) < MinDX) then
                    MinDX = Me%AnalyticWave%X2D(i, j)
                    Me%AnalyticWave%EnteringCell%iStart = i
                    Me%AnalyticWave%EnteringCell%jStart = j                
                endif        
                
            endif                

        enddo
        enddo            
        
        do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
        do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB        
            Me%AnalyticWave%X2D(i, j) = Me%AnalyticWave%X2D(i, j) - MinDX
        enddo
        enddo            
        
        Me%AnalyticWave%EnteringCell%nCells = 0
        
        do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
        do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB        
            
            if (Me%AnalyticWave%CellType(i, j) == EnteringWaveCell_) then

                Me%AnalyticWave%EnteringCell%nCells = Me%AnalyticWave%EnteringCell%nCells + 1

            endif                    
        
        enddo
        enddo      
        
        EnterNcells = Me%AnalyticWave%EnteringCell%nCells 
        
        allocate(Me%AnalyticWave%EnteringCell%i      (EnterNcells))
        allocate(Me%AnalyticWave%EnteringCell%j      (EnterNcells))    
        allocate(Me%AnalyticWave%EnteringCell%TimeLag(EnterNcells))            

        
        Me%AnalyticWave%EnteringCell%dx => DUX
        Me%AnalyticWave%EnteringCell%dy => DVY
        
        Me%AnalyticWave%EnteringCell%i      (1) = Me%AnalyticWave%EnteringCell%iStart
        Me%AnalyticWave%EnteringCell%j      (1) = Me%AnalyticWave%EnteringCell%jStart
        Me%AnalyticWave%EnteringCell%n          = 1
        Me%AnalyticWave%EnteringCell%TimeLag(1) = 0.     
        
        Me%AnalyticWave%TlagMissing (Me%AnalyticWave%EnteringCell%iStart,               &
                                     Me%AnalyticWave%EnteringCell%jStart) = 1        
           
        
        !Compute TimeLag
        !call ComputeTimeLag(i    = Me%AnalyticWave%EnteringCell%i      (1),             &
        !                    j    = Me%AnalyticWave%EnteringCell%j      (1),             &
        !                    Tlag = Me%AnalyticWave%EnteringCell%TimeLag(1)) 
        
        n = 0
        
        do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
        do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB        
            
            if (Me%AnalyticWave%CellType(i, j) == EnteringWaveCell_) then
                n = n + 1
                Me%AnalyticWave%EnteringCell%TimeLag(n) =  Me%AnalyticWave%X2D(i, j) / Me%AnalyticWave%Celerity(i, j)
                Me%AnalyticWave%EnteringCell%i      (n) = i
                Me%AnalyticWave%EnteringCell%j      (n) = j
            endif                    
        
        enddo
        enddo              
                            
        nullify(Me%AnalyticWave%EnteringCell%dx)
        nullify(Me%AnalyticWave%EnteringCell%dy)
                            

        call UngetHorizontalGrid(HorizontalGridID = ObjHorizontalGridAux,               &
                                 Polygon          = InteriorCellsArea,                  &
                                 STAT             = STAT_CALL)                          
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR210'
        
        if (GetDDecompON(Me%ObjHorizontalGrid)) then        
            call KillHorizontalGrid(HorizontalGridID = ObjHorizontalGridAux,            &
                                    STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR220'           
        endif            

        call UngetHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,               & 
                                 Array            = XX2D,                               &
                                 STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR230'

        call UngetHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,               & 
                                 Array            = YY2D,                               &
                                 STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR240'

        call UngetHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,               & 
                                 Array            = DUX,                                &
                                 STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR250'

        call UngetHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,               & 
                                 Array            = DVY,                                &
                                 STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAnalyticWave - ModuleFillMatrix - ERR260'


        deallocate(Segment)
        deallocate(Point)
        deallocate(Point0)        
        
        deallocate(Aux2D)
        
        
        Me%AnalyticWave%ON = .true.
        

    end subroutine ConstructAnalyticWave

    !--------------------------------------------------------------------------
    
    recursive subroutine ComputeTimeLag(i, j, Tlag)
            
        !Arguments-------------------------------------------------------------
        integer                                         :: i, j
        real                                            :: Tlag
        !Local-----------------------------------------------------------------            
        integer                                         :: i1, j1
        real                                            :: T
        real                                            :: CellRotationX, dr1, dr2
        real                                            :: dx1, dy1, dx2, dy2, dt1, dt2, NewTime
        real                                            :: wc1, wc2
        integer                                         :: STAT_CALL
        integer                                         :: nn, ns, nw, ne        
        integer                                         :: in, is, iw, ie
        integer                                         :: jn, js, jw, je        

        !Begin-----------------------------------------------------------------    
    

        
        T  = Tlag
        i1 = i
        j1 = j

        dx1 = Me%AnalyticWave%EnteringCell%dx(i1, j1) / 2.
        dy1 = Me%AnalyticWave%EnteringCell%dy(i1, j1) / 2.      
        wc1 = Me%AnalyticWave%Celerity       (i1, j1)
        
        call GetCellRotation(Me%ObjHorizontalGrid, i1, j1, CellRotationX, STAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) then
            stop 'ComputeTimeLag - ModuleFillMatrix - ERR10'
        endif

        dr1 = Me%AnalyticWave%Direction - CellRotationX
        
        !North
        in = i1 + 1
        jn = j1 
        
        if (Me%AnalyticWave%CellType    (in, jn) == EnteringWaveCell_ .and.             &
            Me%AnalyticWave%TlagMissing (in, jn) == 0) then
            
            call GetCellRotation(Me%ObjHorizontalGrid, in, jn, CellRotationX, STAT = STAT_CALL)
            
            if (STAT_CALL /= SUCCESS_) then
                stop 'ComputeTimeLag - ModuleFillMatrix - ERR20'
            endif
            
            dr2 = Me%AnalyticWave%Direction - CellRotationX                 
                 
            dy2 = Me%AnalyticWave%EnteringCell%dy(in, jn) / 2.      
            wc2 = Me%AnalyticWave%Celerity       (in, jn)

            dt1 = dy1 * sin(dr1) / wc1
            dt2 = dy2 * sin(dr2) / wc2
                         
            NewTime = T + dt1 + dt2
            
            Me%AnalyticWave%EnteringCell%n = Me%AnalyticWave%EnteringCell%n + 1
            
            nn  = Me%AnalyticWave%EnteringCell%n

            Me%AnalyticWave%TlagMissing(in, jn) = nn

            Me%AnalyticWave%EnteringCell%i(nn)       = in
            Me%AnalyticWave%EnteringCell%j(nn)       = jn
            Me%AnalyticWave%EnteringCell%TimeLag(nn) = NewTime
            
        endif
        
        !South
        is = i1 - 1
        js = j1 
        
        if (Me%AnalyticWave%CellType    (is, js) == EnteringWaveCell_ .and.             &
            Me%AnalyticWave%TlagMissing (is, js) == 0) then

            
            call GetCellRotation(Me%ObjHorizontalGrid, is, js, CellRotationX, STAT = STAT_CALL)
            
            if (STAT_CALL /= SUCCESS_) then
                stop 'ComputeTimeLag - ModuleFillMatrix - ERR30'
            endif
            
            dr2 = Me%AnalyticWave%Direction - CellRotationX                 
                 
            dy2 = Me%AnalyticWave%EnteringCell%dy(is, js) / 2.      
            wc2 = Me%AnalyticWave%Celerity       (is, js)

            dt1 = - dy1 * sin(dr1) / wc1
            dt2 = - dy2 * sin(dr2) / wc2
                         
            NewTime = T + dt1 + dt2
        
            Me%AnalyticWave%EnteringCell%n = Me%AnalyticWave%EnteringCell%n + 1

            ns  = Me%AnalyticWave%EnteringCell%n

            Me%AnalyticWave%TlagMissing(is, js) = ns

            Me%AnalyticWave%EnteringCell%i(ns)       = is
            Me%AnalyticWave%EnteringCell%j(ns)       = js
            Me%AnalyticWave%EnteringCell%TimeLag(ns) = NewTime

        endif   
        
        !West
        iw = i1
        jw = j1 - 1 
        
        if (Me%AnalyticWave%CellType    (iw, jw) == EnteringWaveCell_ .and.             &
            Me%AnalyticWave%TlagMissing (iw, jw) == 0) then
        
            call GetCellRotation(Me%ObjHorizontalGrid, iw, jw, CellRotationX, STAT = STAT_CALL)
            
            if (STAT_CALL /= SUCCESS_) then
                stop 'ComputeTimeLag - ModuleFillMatrix - ERR40'
            endif
            
            dr2 = Me%AnalyticWave%Direction - CellRotationX                 
                 
            dx2 = Me%AnalyticWave%EnteringCell%dx(iw, jw) / 2.      
            wc2 = Me%AnalyticWave%Celerity       (iw, jw)

            dt1 = - dx1 * cos(dr1) / wc1
            dt2 = - dx2 * cos(dr2) / wc2
                         
            NewTime = T + dt1 + dt2

            Me%AnalyticWave%EnteringCell%n = Me%AnalyticWave%EnteringCell%n + 1
            
            nw  = Me%AnalyticWave%EnteringCell%n
            
            Me%AnalyticWave%TlagMissing(iw, jw) = nw            

            Me%AnalyticWave%EnteringCell%i(nw)       = iw
            Me%AnalyticWave%EnteringCell%j(nw)       = jw
            Me%AnalyticWave%EnteringCell%TimeLag(nw) = NewTime
       
        
        endif             

        !East
        ie = i1
        je = j1 + 1 
        
        if (Me%AnalyticWave%CellType    (ie, je) == EnteringWaveCell_ .and.             &
            Me%AnalyticWave%TlagMissing (ie, je) == 0) then
        
            call GetCellRotation(Me%ObjHorizontalGrid, ie, je, CellRotationX, STAT = STAT_CALL)
            
            if (STAT_CALL /= SUCCESS_) then
                stop 'ComputeTimeLag - ModuleFillMatrix - ERR50'
            endif
            
            dr2 = Me%AnalyticWave%Direction - CellRotationX                 
                 
            dx2 = Me%AnalyticWave%EnteringCell%dx(ie, je) / 2.      
            wc2 = Me%AnalyticWave%Celerity       (ie, je)

            dt1 = dx1 * cos(dr1) / wc1
            dt2 = dx2 * cos(dr2) / wc2
                         
            NewTime = T + dt1 + dt2
            
            Me%AnalyticWave%EnteringCell%n = Me%AnalyticWave%EnteringCell%n + 1
            
            ne  = Me%AnalyticWave%EnteringCell%n
            
            Me%AnalyticWave%TlagMissing(ie, je) = ne
            
            Me%AnalyticWave%EnteringCell%i(ne)       = ie
            Me%AnalyticWave%EnteringCell%j(ne)       = je
            Me%AnalyticWave%EnteringCell%TimeLag(ne) = NewTime
            
        endif       

        !Compute TimeLag
        if (Me%AnalyticWave%CellType    (in, jn) == EnteringWaveCell_ .and.             &
            Me%AnalyticWave%TlagMissing (in, jn) == nn) then
            call ComputeTimeLag(i=in,j=jn,Tlag=Me%AnalyticWave%EnteringCell%TimeLag(nn))
        endif            

        if (Me%AnalyticWave%CellType    (is, js) == EnteringWaveCell_ .and.             &
            Me%AnalyticWave%TlagMissing (is, js) == ns) then
            call ComputeTimeLag(i=is,j=js,Tlag=Me%AnalyticWave%EnteringCell%TimeLag(ns))
        endif            
        
        if (Me%AnalyticWave%CellType    (iw, jw) == EnteringWaveCell_ .and.             &
            Me%AnalyticWave%TlagMissing (iw, jw) == nw) then
            call ComputeTimeLag(i=iw,j=jw,Tlag=Me%AnalyticWave%EnteringCell%TimeLag(nw))
        endif            
        
        if (Me%AnalyticWave%CellType    (ie, je) == EnteringWaveCell_ .and.             &
            Me%AnalyticWave%TlagMissing (ie, je) == ne) then
            call ComputeTimeLag(i=ie,j=je,Tlag=Me%AnalyticWave%EnteringCell%TimeLag(ne))
        endif                    
    
    end subroutine ComputeTimeLag

    !--------------------------------------------------------------------------
    
    function CheckSponge(PointsToFill2D, PointsToFill3D, sp, di1, di2, dj1, dj2, i, j, k)
    
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D
        integer                                         :: sp, di1, di2, dj1, dj2, i, j, k
        logical                                         :: CheckSponge

        !Local-----------------------------------------------------------------            


        !Begin-----------------------------------------------------------------    
    
i2:     if (Me%Dim == Dim2D) then

            if (sum(PointsToFill2D(i+di1:i+di2, j+dj1:j+dj2)) == sp) then
                
                CheckSponge = .true.
                
            else
            
                CheckSponge = .false.
                
            endif
            
        else i2
                                    
       
            if (sum(PointsToFill3D(i+di1:i+di2, j+dj1:j+dj2,k)) == sp) then
                
                CheckSponge = .true.
                
            else
            
                CheckSponge = .false.

            endif        

        endif  i2      
    
    end function CheckSponge

    !--------------------------------------------------------------------------    

    subroutine ConstructSpaceTimeSerie (ExtractType, PointsToFill2D, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        integer                                         :: ExtractType
        integer, dimension(:, :   ), pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D
        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL, i
        integer                                         :: iflag, file, column
        character(len = PathLength), dimension(3)       :: Filename   = " "
        character(len = StringLength), dimension(3)     :: DataColumn = " "
        type(T_TimeSerie), pointer                      :: CurrentTimeSerie, NewTimeSerie
        integer                                         :: nTimeSeries   = 1
        logical                                         :: exist
        !Begin----------------------------------------------------------------
        
        FileName(:) = " "
        DataColumn(:) = " "
        Me%nTimeSeries = 0
        nTimeSeries   = 1
        !Gets the name of the data file
        !If vectorial use different keywords for x and y since GetData would not allow spaces
        !or would interpret as a new filename    
        
        !always search for one filename
        call GetData(FileName(1),                                                &
                        Me%ObjEnterData , iflag,                                    &
                        SearchType   = ExtractType,                                 &
                        keyword      = 'FILENAME',                                  &
                        ClientModule = 'ModuleFillMatrix',                          &
                        STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR01a'       
        if (iflag==0) then
            
            !if scalar needs to have it defined
            if (.not. Me%VectorialProp) then
                write(*,*) 'TIMESERIE File Name not given for property '//trim(Me%PropertyID%Name)           
                stop       'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR01b'     
                
            !if vectorial, and not one file defined, search for two files for each component
            else
                
                call GetData(FileName(1),                                                &
                             Me%ObjEnterData , iflag,                                    &
                             SearchType   = ExtractType,                                 &
                             keyword      = 'FILENAME_X',                                &
                             ClientModule = 'ModuleFillMatrix',                          &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR01c'              
                if (iflag==0)then
                    write(*,*) 'TIMESERIE FILENAME_X keyword not given not given for vectorial property '//trim(Me%PropertyID%Name)
                    stop       'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR01d'                
                endif
                inquire (file=trim(FileName(1)), exist = exist)
                if (.not. exist) then
                    write(*,*)'Could not find file '//trim(FileName(1))
                    stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR01e'
                endif  
                
                call GetData(FileName(2),                                                &
                             Me%ObjEnterData , iflag,                                    &
                             SearchType   = ExtractType,                                 &
                             keyword      = 'FILENAME_Y',                                &
                             ClientModule = 'ModuleFillMatrix',                          &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR01f'     
                if (iflag==0)then
                    write(*,*) 'TIMESERIE FILENAME_Y keyword not given not given for vectorial property '//trim(Me%PropertyID%Name)
                    stop       'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR01g'                
                endif
                inquire (file=trim(FileName(2)), exist = exist)
                if (.not. exist) then
                    write(*,*)'Could not find file '//trim(FileName(2))
                    stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR01h'
                endif 
                
                if (Me%Dim == Dim3D) then
                    
                    call GetData(FileName(3),                                                &
                                 Me%ObjEnterData , iflag,                                    &
                                 SearchType   = ExtractType,                                 &
                                 keyword      = 'FILENAME_Z',                                &
                                 ClientModule = 'ModuleFillMatrix',                          &
                                 STAT         = STAT_CALL)                                      
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR01i'     
                    if (iflag==0)then
                        !write(*,*) 'TIMESERIE FILENAME_Z keyword not given not given for vectorial property '//trim(Me%PropertyID%Name)
                        !stop       'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR01j'                
                    else
                        inquire (file=trim(FileName(3)), exist = exist)
                        if (.not. exist) then
                            write(*,*)'Could not find file '//trim(FileName(3))
                            stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR01k'
                        endif      
                    endif
                endif
                
            endif       
        
        else
            inquire (file=trim(FileName(1)), exist = exist)
            if (.not. exist) then
                write(*,*)'Could not find file '//trim(FileName(1))
                stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR01l'
            endif                 
        endif        
        
        !mandatory two columns
        if (Me%VectorialProp) then            
                        
            call GetData(DataColumn(1),                                                &
                         Me%ObjEnterData , iflag,                                    &
                         SearchType   = ExtractType,                                 &
                         keyword      = 'DATA_COLUMN_X',                             &
                         ClientModule = 'ModuleFillMatrix',                          &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR02a'              
            if (iflag==0)then
                write(*,*) 'TIMESERIE DATA_COLUMN_X keyword not given not given for vectorial property '//trim(Me%PropertyID%Name)
                stop       'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR02b'                
            endif
            
            call GetData(DataColumn(2),                                                &
                         Me%ObjEnterData , iflag,                                    &
                         SearchType   = ExtractType,                                 &
                         keyword      = 'DATA_COLUMN_Y',                             &
                         ClientModule = 'ModuleFillMatrix',                          &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR02c'     
            if (iflag==0)then
                write(*,*) 'TIMESERIE DATA_COLUMN_Y keyword not given not given for vectorial property '//trim(Me%PropertyID%Name)
                stop       'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR02d'                
            endif            
            
            if (Me%Dim == Dim3D) then
                Me%UseZ = .true.
                call GetData(DataColumn(3),                                              &
                             Me%ObjEnterData , iflag,                                    &
                             SearchType   = ExtractType,                                 &
                             keyword      = 'DATA_COLUMN_Z',                             &
                             ClientModule = 'ModuleFillMatrix',                          &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR02e'
                if (iflag==0)then
                    Me%UseZ = .false.
                    !write(*,*) 'TIMESERIE DATA_COLUMN_Z keyword not given not given for vectorial property '//trim(Me%PropertyID%Name)                    
                    !stop       'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR02f'                
                endif                    
                
                !verify that user provided W omponent
                if (Me%UseZ .and. .not. associated(Me%Matrix3DW)) then
                    write(*,*) 'Constructing vectorial property that needs W component to be given'
                    stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR02f'            
                endif                  
                
            endif
            
            !Default
            nTimeSeries = 2            
            !if second filename not defined, is the same as previous
            if (FileName(2) == " ") FileName(2) = FileName(1)            
           
            if (Me%Dim == Dim3D .and. Me%UseZ) then
                nTimeSeries = 3
                if (FileName(3) == " ") FileName(3) = FileName(1)
            endif            
            
            
        else            
                   
            call GetData(DataColumn(1),                                              &
                         Me%ObjEnterData , iflag,                                    &
                         SearchType   = ExtractType,                                 &
                         keyword      = 'DATA_COLUMN',                               &
                         ClientModule = 'ModuleFillMatrix',                          &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR02g'           

            if (iflag==0)then
                write(*,*)'Data Column not given for property'
                write(*,*)trim(Me%PropertyID%Name)
                stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR02h'
            endif            
            
        endif                      
        
        do i = 1, nTimeSeries

            nullify  (NewTimeSerie)
            allocate (NewTimeSerie)
            nullify  (NewTimeSerie%Next)
            nullify  (NewTimeSerie%Prev)
                
            !convert to int
            read(DataColumn(i),*,iostat=STAT_CALL) column        

            NewTimeSerie%Filename = Filename(i)
            NewTimeSerie%Column   = column           
            
            
            call InsertTimeSerieToList(NewTimeSerie)
                
        enddo
        
        if (Me%PredictDTMethod == 2) then
            Me%DTForNextEvent       = -null_real
            Me%PredictedDT          = -null_real
            Me%DTForNextDataset     = -null_real
            Me%NextValueForDTPred   = 0.0        
        endif
        
        CurrentTimeSerie => Me%FirstTimeSerie
        file = 1
        do while (associated(CurrentTimeSerie))              
            
            !Associate Matrix2D and Me%Matrix3D to the input field ones
            if (Me%VectorialProp .or. Me%RotateAngleToGrid) call AssociateMatrixes(file)                        
            
            !Starts Time Serie
            call StartTimeSerieInput(CurrentTimeSerie%ObjTimeSerie, &
                                     CurrentTimeSerie%FileName,     &
                                     Me%ObjTime,                    &
                                     CheckDates = Me%CheckDates,    &
                                     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR050'
        
            call GetTimeSerieDataValues(CurrentTimeSerie%ObjTimeSerie,       &
                                        CurrentTimeSerie%NumberOfInstants,   &
                                        STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR060'
            
            Me%RemainsConstant = .true.
            !if only one instant is found then values remain constant
            if(CurrentTimeSerie%NumberOfInstants == 1) then
            
                CurrentTimeSerie%RemainsConstant = .true.
                call GetTimeSerieValueForIndex (CurrentTimeSerie%ObjTimeSerie,   &
                                                1,                               &
                                                CurrentTimeSerie%Column,         &
                                                CurrentTimeSerie%CurrentValue,   &
                                                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR065'             
            else
                !if any of the series do not remain constant then global is false
                Me%RemainsConstant = .false.
            endif
            
            !this not used in case of vectorial since general properties are written (Me%NextEventStart ...) and the complexity for now
            !is not needed
            if (Me%PredictDTMethod == 2) then
        
                if (Me%PropertyID%IDNumber == WindDirection_) then
                    write(*,*) 'The method 2 to predict DT do not works with WindDirection_ property'
                    stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR070'            
                endif
        
                if(CurrentTimeSerie%NumberOfInstants > 1) then
                             
                    CurrentTimeSerie%PreviousInstant = 1                
                    call GetTimeSerieTimeOfDataset(CurrentTimeSerie%ObjTimeSerie,        &
                                                   CurrentTimeSerie%PreviousInstant,     &
                                                   CurrentTimeSerie%PreviousTime,        &
                                                   STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)    &
                        stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR080'
                                                 
                    CurrentTimeSerie%NextInstant = 1
                    CurrentTimeSerie%NextTime    = CurrentTimeSerie%PreviousTime
                
                    call ActualizeTimeSerieValues(CurrentTimeSerie)
                
                    CurrentTimeSerie%NextEventStart        = CurrentTimeSerie%PreviousTime
                    CurrentTimeSerie%NextEventEnd          = CurrentTimeSerie%PreviousTime
                    CurrentTimeSerie%NextValueForDTPred    = CurrentTimeSerie%NextValue
                    CurrentTimeSerie%DTForNextEvent        = 0.0
            
                endif

            endif
            
            if (Me%Dim == Dim2D) then
                call ModifySpaceTimeSerie    (CurrentTimeSerie, PointsToFill2D = PointsToFill2D) 
            else
                call ModifySpaceTimeSerie    (CurrentTimeSerie, PointsToFill3D = PointsToFill3D) 
            endif            
            
            
            file = file + 1
            CurrentTimeSerie => CurrentTimeSerie%Next
            
        enddo        
        
    end subroutine ConstructSpaceTimeSerie

    !--------------------------------------------------------------------------
    
    subroutine AssociateMatrixes(file) 
    
        !Arguments-------------------------------------------------------------
        integer                                         :: file
    
        !External--------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
    
        !Begin-----------------------------------------------------------------    
    
        !Vectorial prop get user field data by associating Me%Matrix2D and Me%Matrix3D to it
        !Routines change Me%Matrix2D and Me%Matrix3D so they change the correct one that will be 
        !available to origin module (for output only)
        if (Me%VectorialProp) then
            if (file == 1) then
                if (Me%Dim == Dim2D) then
                    Me%Matrix2D => Me%Matrix2DX
                else
                    Me%Matrix3D => Me%Matrix3DX
                endif
            else if (file == 2) then
                if (Me%Dim == Dim2D) then
                    Me%Matrix2D => Me%Matrix2DY
                else
                    Me%Matrix3D => Me%Matrix3DY
                endif
            else if (file == 3) then !3D
                !no need for Z matrix since w is independent of horizotal referential
                !fill directly result matrix
                Me%Matrix3D => Me%Matrix3DW              
            endif                    
        else if (Me%RotateAngleToGrid) then
            if (Me%Dim == Dim2D) then
                Me%Matrix2D => Me%Matrix2DFieldAngle
            else
                Me%Matrix3D => Me%Matrix3DFieldAngle
            endif                    
        endif
    
    end subroutine AssociateMatrixes

    !-------------------------------------------------------------------------- 

    subroutine ConstructHDFInput (ExtractType, ClientID, PointsToFill2D, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        integer                                         :: ExtractType
        integer                                         :: ClientID
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: HDF5_READ
        type(T_Time)                                    :: Now, CurrentTime

        !Local-----------------------------------------------------------------
        type (T_Field4D), pointer                       :: CurrentHDF
        integer                                         :: file
        !Begin-----------------------------------------------------------------
        
        !get general options and list of HDF (one if scalar, two if vectorial prop)
        call ReadOptionsHDFinput(ExtractType, ClientID) 
        
        
        CurrentHDF => Me%FirstHDF
        file = 1
        do while (associated(CurrentHDF))
            
            !Associate Matrix2D and Me%Matrix3D to the input field ones
            if (Me%VectorialProp .or. Me%RotateAngleToGrid) call AssociateMatrixes(file)            
            
            call AllocateHDFInput(CurrentHDF)
             
                
if4D:       if (CurrentHDF%Field4D) then

                call BuildField4D(ExtractType, ClientID, PointsToFill2D, PointsToFill3D, CurrentHDF)    
        
            else if4D
        
                call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

                call ConstructHDF5 (CurrentHDF%ObjHDF5, trim(CurrentHDF%FileName), HDF5_READ, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR180'

                if (Me%CheckHDF5_File) then
        
                    call GetHDF5AllDataSetsOK (HDF5ID   = CurrentHDF%ObjHDF5,           &
                                               STAT     = STAT_CALL)                                      
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR185'
            
                endif                

                call GetHDF5GroupNumberOfItems(CurrentHDF%ObjHDF5, "/Time", &
                                               CurrentHDF%NumberOfInstants, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR190'
        
            endif if4D
        
            CurrentHDF%StartTime = HDF5TimeInstant(1, CurrentHDF)
            CurrentHDF%EndTime   = HDF5TimeInstant(CurrentHDF%NumberOfInstants, CurrentHDF)
            
            !if only one instant is found then values remain constant
            Me%RemainsConstant = .true.
            if(CurrentHDF%NumberOfInstants == 1 .and. .not. CurrentHDF%HarmonicsON) then
                CurrentHDF%RemainsConstant = .true.
            else
                !if any of the hdf do not remain constant then global is false
                Me%RemainsConstant = .false.
            endif            
        

        
            call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR200'
        
            !Backtracking time inversion is also done in the ModuleField4D    
            if (Me%BackTracking .and. .not. CurrentHDF%Field4D) then  
                call BacktrackingTime(Now)
            else   
                Now = CurrentTime
            endif

            !Initial field
    i1:     if (CurrentHDF%Generic4D%ON .or. (CurrentHDF%Field4D .and. CurrentHDF%HarmonicsOn)) then
        
                if (Me%Dim == Dim2D) then
                    call ModifyHDFInput2D (PointsToFill2D, CurrentHDF) 
                else
                    call ModifyHDFInput3D (PointsToFill3D, CurrentHDF)
                endif
            
            else i1
            
    i2:         if(CurrentHDF%NumberOfInstants > 1)then

    i3:             if (Me%PredictDTMethod == 2) then
                
                        call ConstructHDFPredictDTMethod2(Now, PointsToFill3D, PointsToFill2D, CurrentHDF)
                
                    else i3

                        call ConstructHDFPredictDTMethod1(Now, CurrentHDF)

                    endif i3

                elseif(CurrentHDF%NumberOfInstants == 1)then i2


                    call ConstructHDFOneInstant(Now, CurrentHDF) 

                else i2
                    write(*,*)
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'No time information found'
                    stop 'ConstructHDFInput - ModuleFillMatrix - ERR300'
                end if i2

    i4:         if (Me%PredictDTMethod == 1) then
            
    i5:             if(Me%Dim == Dim2D)then

                        call DTMethod_1_PrevNext2D(Now, PointsToFill2D, CurrentHDF)    

                    else i5

                        call DTMethod_1_PrevNext3D(Now, PointsToFill3D, CurrentHDF)    

                    end if i5
                
                endif i4

            endif i1
            
                
            file = file + 1
            CurrentHDF => CurrentHDF%Next
        enddo            
            
    end subroutine ConstructHDFInput

    !-----------------------------------------------------------------------------------

    subroutine AllocateHDFInput(CurrentHDF)
    
         !Arguments------------------------------------------------------------
        type(T_Field4D)                                 :: CurrentHDF
        !Local----------------------------------------------------------------- 
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB

        !Begin-----------------------------------------------------------------           


        
        nullify(CurrentHDF%PreviousField2D, CurrentHDF%NextField2D)
        nullify(CurrentHDF%PreviousField3D, CurrentHDF%NextField3D)
        
i0:     if(Me%Dim == Dim2D)then

            ILB = Me%Size2D%ILB
            IUB = Me%Size2D%IUB
            JLB = Me%Size2D%JLB
            JUB = Me%Size2D%JUB

            allocate(CurrentHDF%PreviousField2D (ILB:IUB, JLB:JUB))
            allocate(CurrentHDF%NextField2D     (ILB:IUB, JLB:JUB))
            allocate(CurrentHDF%Array2D         (ILB:IUB, JLB:JUB))

            CurrentHDF%PreviousField2D(:,:) = FillValueReal
            CurrentHDF%NextField2D    (:,:) = FillValueReal

        else i0

            ILB = Me%Size3D%ILB
            IUB = Me%Size3D%IUB
            JLB = Me%Size3D%JLB
            JUB = Me%Size3D%JUB
            KLB = Me%Size3D%KLB
            KUB = Me%Size3D%KUB

            allocate(CurrentHDF%PreviousField3D (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(CurrentHDF%NextField3D     (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(CurrentHDF%Array3D         (ILB:IUB, JLB:JUB, KLB:KUB))

            CurrentHDF%PreviousField3D(:,:,:) = FillValueReal
            CurrentHDF%NextField3D    (:,:,:) = FillValueReal

        endif i0    
            

    end subroutine AllocateHDFInput

    !-----------------------------------------------------------------------------------

    subroutine ReadOptionsHDFinput(ExtractType, ClientID)
    
         !Arguments------------------------------------------------------------
        integer                                         :: ExtractType   
        integer                                         :: ClientID
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, iflag, i
        integer                                         :: ILB, IUB, JLB, JUB
        logical                                         :: MasterOrSlave, LastGroupEqualField
        character(len = PathLength  ), dimension(3)     :: FileName  = " "        
        character(len = StringLength), dimension(3)     :: FieldName = " "
        type(T_Field4D), pointer                        :: NewHDF, CurrentHDF
        integer                                         :: nHDFs           = 1
        logical                                         :: exist
        real                                            :: DT           
        !Begin-----------------------------------------------------------------           
        
        FileName(:) = " "
        FieldName(:) = " "
        Me%nHDFs = 0
        nHDFs           = 1
        
i0:     if(Me%Dim == Dim2D)then

            ILB = Me%Size2D%ILB
            IUB = Me%Size2D%IUB
            JLB = Me%Size2D%JLB
            JUB = Me%Size2D%JUB

        else i0

            ILB = Me%Size3D%ILB
            IUB = Me%Size3D%IUB
            JLB = Me%Size3D%JLB
            JUB = Me%Size3D%JUB
        
        endif i0
    
!        !interpolation of angles is not done
!        if (Me%PropertyID%IDNumber == WindDirection_) then 
!            write(*,*) 'Trying to construct an HDF for property wind direction. Not available option.'
!            stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR10'              
!        endif
        
        !Always search for one filename
        if (.not. Me%ArgumentFileName) then
            
            call GetData(FileName(1),                                                      &
                            Me%ObjEnterData , iflag,                                       &
                            SearchType   = ExtractType,                                    &
                            keyword      = 'FILENAME',                                     &
                            ClientModule = 'ModuleFillMatrix',                             &
                            STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR20'            

            if (iflag==0)then
                
                !need to have defined file
                if (.not. Me%VectorialProp) then                                
                    write(*,*)'HDF filename not given for property '//trim(Me%PropertyID%Name)
                    stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR30'
                
                !If not one file, search foreach file for each component
                else
                    
                    call GetData(FileName(1),                                                &
                                 Me%ObjEnterData , iflag,                                    &
                                 SearchType   = ExtractType,                                 &
                                 keyword      = 'FILENAME_X',                                &
                                 ClientModule = 'ModuleFillMatrix',                          &
                                 STAT         = STAT_CALL)                                      
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR40'              
                    if (iflag==0)then
                        write(*,*) 'HDF FILENAME_X keyword not given not given for vectorial property '//trim(Me%PropertyID%Name)
                        stop       'ReadOptionsHDFinput - ModuleFillMatrix - ERR50'                
                    endif
                    inquire (file=trim(FileName(1)), exist = exist)
                    if (.not. exist) then
                        write(*,*)'Could not find file '//trim(FileName(1))
                        stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR60'
                    endif    
                    
                    call GetData(FileName(2),                                                &
                                 Me%ObjEnterData , iflag,                                    &
                                 SearchType   = ExtractType,                                 &
                                 keyword      = 'FILENAME_Y',                                &
                                 ClientModule = 'ModuleFillMatrix',                          &
                                 STAT         = STAT_CALL)                                      
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR70'     
                    if (iflag==0)then
                        write(*,*) 'HDF FILENAME_Y keyword not given not given for vectorial property '//trim(Me%PropertyID%Name)
                        stop       'ReadOptionsHDFinput - ModuleFillMatrix - ERR80'                
                    endif                                
                    inquire (file=trim(FileName(2)), exist = exist)
                    if (.not. exist) then
                        write(*,*)'Could not find file '//trim(FileName(2))
                        stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR90'
                    endif  
                    
                    if (Me%Dim == Dim3D) then                        
                        call GetData(FileName(3),                                                &
                                     Me%ObjEnterData , iflag,                                    &
                                     SearchType   = ExtractType,                                 &
                                     keyword      = 'FILENAME_Z',                                &
                                     ClientModule = 'ModuleFillMatrix',                          &
                                     STAT         = STAT_CALL)                                      
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR02g'     
                        if (iflag==0)then
                            !write(*,*) 'HDF FILENAME_Z keyword not given not given for vectorial property '//trim(Me%PropertyID%Name)
                            !stop       'ReadOptionsHDFinput - ModuleFillMatrix - ERR02h'                
                        else                                
                            inquire (file=trim(FileName(3)), exist = exist)
                            if (.not. exist) then
                                write(*,*)'Could not find file '//trim(FileName(3))
                                stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR100'
                            endif    
                        endif
                    endif
                    
                endif
            else
        
                inquire (file=trim(FileName(1)), exist = exist)
                if (.not. exist) then
                    write(*,*)'Could not find file '//trim(FileName(1))
                    stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR110'
                endif                 
            
            endif
        else
            
            if (.not. Me%VectorialProp) then 
                FileName(1) = Me%FileNameHDF(1)                
            else
                FileName(1) = Me%FileNameHDF(1)
                FileName(2) = Me%FileNameHDF(2)                
            endif
        endif           
                
        
        if (Me%VectorialProp) then                                    
            
            call GetData(FieldName(1),                                               &
                         Me%ObjEnterData , iflag,                                    &
                         SearchType   = ExtractType,                                 &
                         keyword      = 'HDF_FIELD_NAME_X',                          &
                         ClientModule = 'ModuleFillMatrix',                          &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR120'              
            if (iflag==0)then
                write(*,*) 'HDF HDF_FIELD_NAME_X keyword not given not given for vectorial property '//trim(Me%PropertyID%Name) 
                stop       'ReadOptionsHDFinput - ModuleFillMatrix - ERR04b'                
            endif
            
            call GetData(FieldName(2),                                               &
                         Me%ObjEnterData , iflag,                                    &
                         SearchType   = ExtractType,                                 &
                         keyword      = 'HDF_FIELD_NAME_Y',                          &
                         ClientModule = 'ModuleFillMatrix',                          &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR130'     
            if (iflag==0)then
                write(*,*) 'HDF HDF_FIELD_NAME_Y keyword not given not given for vectorial property '//trim(Me%PropertyID%Name)
                stop       'ReadOptionsHDFinput - ModuleFillMatrix - ERR04d'                
            endif            
            
            if (Me%Dim == Dim3D) then
                Me%UseZ = .true.
                call GetData(FieldName(3),                                               &
                             Me%ObjEnterData , iflag,                                    &
                             SearchType   = ExtractType,                                 &
                             keyword      = 'HDF_FIELD_NAME_Z',                          &
                             ClientModule = 'ModuleFillMatrix',                          &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR140'     
                if (iflag==0)then
                    Me%UseZ = .false.
                    !write(*,*) 'HDF HDF_FIELD_NAME_Z keyword not given not given for vectorial property '//trim(Me%PropertyID%Name)
                    !stop       'ConstructHDFInput - ModuleFillMatrix - ERR04f'                
                endif    
                
                !verify that user provided W omponent
                if (Me%UseZ .and. .not. associated(Me%Matrix3DW)) then
                    write(*,*) 'Constructing vectorial property that needs W component to be given'
                    stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR150'            
                endif                  
                
            endif
            
            !Default
            nHDFs = 2            
            !if second filename not defined, is the same as previous
            if (FileName(2) == " ") FileName(2) = FileName(1)
            
            if (Me%Dim == Dim3D .and. Me%UseZ) then
                nHDFs = 3
                if (FileName(3) == " ") FileName(3) = FileName(1)
            endif            
            
            
        else        
                
            call GetData(FieldName(1),                                                      &
                         Me%ObjEnterData , iflag,                                           &
                         SearchType   = ExtractType,                                        &
                         keyword      = 'HDF_FIELD_NAME',                                   &
                         default      = trim(Me%PropertyID%Name),                           &
                         ClientModule = 'ModuleFillMatrix',                                 &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR160' 
            
        endif
            
        call GetData(Me%CheckHDF5_File,                                                 &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'CHECK_HDF5_FILE',                                  &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR165'         
        
        do i = 1, nHDFs

            nullify  (NewHDF)
            allocate (NewHDF)
            nullify  (NewHDF%Next)
            nullify  (NewHDF%Prev)

            NewHDF%Filename  = Filename(i)                    
            NewHDF%FieldName = FieldName(i)            
            
            
            call InsertHDFToList(NewHDF)
                
        enddo        
        
        CurrentHDF => Me%FirstHDF
        
d1:     do while (associated(CurrentHDF))         
        
            call GetData(CurrentHDF%Generic4D%ON,                                           &
                         Me%ObjEnterData , iflag,                                           &
                         SearchType   = ExtractType,                                        &
                         keyword      = '4D',                                               &
                         default      = .false.,                                            &
                         ClientModule = 'ModuleFillMatrix',                                 &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR170'


            if (CurrentHDF%Generic4D%ON) call Generic4thDimension(ExtractType, CurrentHDF)
        

            call GetData(CurrentHDF%VGroupPath,                                             &
                         Me%ObjEnterData , iflag,                                           &
                         SearchType   = ExtractType,                                        &
                         keyword      = 'VGROUP_PATH',                                      &
                         default      = "/Results",                                         &
                         ClientModule = 'ModuleFillMatrix',                                 &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR180'

            call GetData(CurrentHDF%MultiplyingFactor,                                      &
                         Me%ObjEnterData , iflag,                                           &
                         SearchType   = ExtractType,                                        &
                         keyword      = 'MULTIPLYING_FACTOR',                               &
                         default      = 1.,                                                 &
                         ClientModule = 'ModuleFillMatrix',                                 &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR190'
        
            if (iflag == 1)then
                CurrentHDF%HasMultiplyingFactor = .true.
            end if

            call GetData(CurrentHDF%AddingFactor,                                               &
                         Me%ObjEnterData , iflag,                                           &
                         SearchType   = ExtractType,                                        &
                         keyword      = 'ADDING_FACTOR',                                    &
                         default      = 0.,                                                 &
                         ClientModule = 'ModuleFillMatrix',                                 &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR200'
        
            if (iflag == 1)then
                CurrentHDF%HasAddingFactor = .true.
            end if


            call GetData(LastGroupEqualField,                                               &
                         Me%ObjEnterData , iflag,                                           &
                         SearchType   = ExtractType,                                        &
                         keyword      = 'LAST_GROUP_EQUAL_FIELD',                           &
                         default      = .true.,                                             &
                         ClientModule = 'ModuleFillMatrix',                                 &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR230'

            if (LastGroupEqualField)                                                        &
                CurrentHDF%VGroupPath=trim(CurrentHDF%VGroupPath)//"/"//trim(CurrentHDF%FieldName)



            
            call GetData(CurrentHDF%From3Dto2D,                                             &
                         Me%ObjEnterData , iflag,                                           &
                         SearchType   = ExtractType,                                        &
                         keyword      = 'FROM_3D_TO_2D',                                    &
                         default      = .false.,                                            &
                         ClientModule = 'ModuleFillMatrix',                                 &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR240'
           
        
            call GetDDecompParameters(HorizontalGridID = Me%ObjHorizontalGrid, &
                                                  MasterOrSlave    = MasterOrSlave,        &
                                                  STAT             = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR260'        


            call GetData(CurrentHDF%Field4D,                                                    &
                         Me%ObjEnterData , iflag,                                           &
                         SearchType   = ExtractType,                                        &
                         keyword      = 'FIELD4D',                                          &
                         default      = .false.,                                            &
                         ClientModule = 'ModuleFillMatrix',                                 &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR270'
            
            call GetData(CurrentHDF%From2Dto3D,                                             &
                         Me%ObjEnterData , iflag,                                           &
                         SearchType   = ExtractType,                                        &
                         keyword      = 'FROM_2D_TO_3D',                                    &
                         default      = .false.,                                            &
                         ClientModule = 'ModuleFillMatrix',                                 &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR240'
            
            if (CurrentHDF%From2Dto3D .and. .not. CurrentHDF%Field4D) then
        
                allocate(CurrentHDF%ReadField3D(ILB:IUB, JLB:JUB, 0:2), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR250'
             
                CurrentHDF%ReadField3D(:,:,:) = FillValueReal  
            
            endif                  
        
        
            if (CurrentHDF%Field4D) then
                call GetData(CurrentHDF%Extrapolate,                                            &
                             Me%ObjEnterData , iflag,                                       &
                             SearchType   = ExtractType,                                    &
                             keyword      = 'EXTRAPOLATE',                                  &
                             default      = .false.,                                        &
                             ClientModule = 'ModuleFillMatrix',                             &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR280'

                !ExtrapolAverage_ = 1, ExtrapolNearstCell_ = 2
                call GetData(CurrentHDF%ExtrapolateMethod,                              &
                             Me%ObjEnterData , iflag,                                   &
                             SearchType   = ExtractType,                                &
                             keyword      = 'EXTRAPOLATE_METHOD',                       &
                             default      = ExtrapolAverage_,                           &
                             ClientModule = 'ModuleFillMatrix',                         &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR290'  
                
                if (CurrentHDF%ExtrapolateMethod /= ExtrapolAverage_ .and.              &
                    CurrentHDF%ExtrapolateMethod /= ExtrapolNearstCell_ ) then
                    stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR300'  
                endif        
                    
                if     (ExtractType == FromBlock)        then
                    call ReadListFilesFromBlock         (CurrentHDF, ClientID)     
                elseif (ExtractType == FromBlockInBlock) then
                    call ReadListFilesFromBlockInBlock  (CurrentHDF, ClientID)     
                endif
                
            endif
        
            !The adding and multiplying functionalities are also available in ModuleField4D
            !This way it is avoid adding and multiplying twice the AddingFactor and the MultiplyingFactor respectively 
            if (CurrentHDF%Field4D) then
                CurrentHDF%HasMultiplyingFactor = .false.
                CurrentHDF%HasAddingFactor      = .false.
            endif            
                
            call GetData(CurrentHDF%HarmonicsON,                                            &
                         Me%ObjEnterData , iflag,                                           &
                         SearchType   = ExtractType,                                        &
                         keyword      = 'HARMONICS',                                        &
                         default      = .false.,                                            &
                         ClientModule = 'ModuleFillMatrix',                                 &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR310'    
            
            if (CurrentHDF%HarmonicsON) then
        
                CurrentHDF%From2Dto3D = .true.    
                
                call GetData(CurrentHDF%HarmonicsDT,                                        &
                                Me%ObjEnterData , iflag,                                    &
                                SearchType   = ExtractType,                                 &
                                keyword      = 'HARMONICS_DT',                              &
                                default      =  900.,                                       &
                                ClientModule = 'ModuleField4D',                             &
                                STAT         = STAT_CALL)                                      
                if (STAT_CALL /= SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR315'   
                
                call GetComputeTimeStep(TimeID = Me%ObjTime, DT = DT, STAT = STAT_CALL)                                      
                if (STAT_CALL /= SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR317' 
                
                if (CurrentHDF%HarmonicsDT < DT) then
                    CurrentHDF%HarmonicsDT = DT
                endif
                
            endif

            call GetData(CurrentHDF%SpatialInterpolON,                                      &
                         Me%ObjEnterData , iflag,                                           &
                         SearchType   = ExtractType,                                        &
                         keyword      = 'SPATIAL_INTERPOL',                                 &
                         default      = .false.,                                            &
                         ClientModule = 'ModuleFillMatrix',                                 &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR320'        
            
            CurrentHDF%InterpolOnlyVertically = .false. 
            
            if (CurrentHDF%SpatialInterpolON) then

                call GetData(CurrentHDF%InterpolOnlyVertically,                             &
                             Me%ObjEnterData , iflag,                                       &
                             SearchType   = ExtractType,                                    &
                             keyword      = 'INTERPOL_ONLY_VERTICALLY',                     &
                             default      = .false.,                                        &
                             ClientModule = 'ModuleFillMatrix',                             &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR330'              
        
            endif
            
            call GetData(CurrentHDF%GenericYear,                                            &
                         Me%ObjEnterData , iflag,                                           &
                         SearchType   = ExtractType,                                        &
                         keyword      = 'GENERIC_YEAR',                                     &
                         default      = .false.,                                            &
                         ClientModule = 'ModuleFillMatrix',                                 &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptionsHDFinput - ModuleFillMatrix - ERR340'  
        
            if (CurrentHDF%GenericYear) then
                CurrentHDF%CyclicTimeON = .true. 
            endif
        
            if (MasterOrSlave) then
                CurrentHDF%Field4D = .true. 
            endif
            
            CurrentHDF => CurrentHDF%Next
            
        enddo d1
    
    end subroutine ReadOptionsHDFinput
    
    !-----------------------------------------------------------------------------    
    
    
    !--------------------------------------------------------------------------

    subroutine ReadListFilesFromBlock(CurrentHDF, ClientNumber)

        !Arguments-------------------------------------------------------------
        type(T_Field4D), pointer                        :: CurrentHDF
        integer                                         :: ClientNumber 

        !Local-----------------------------------------------------------------
        logical                                         :: BlockInBlockFound
        integer                                         :: STAT_CALL, FirstLine, LastLine, iflag
        integer                                         :: n, FileNumber 
        !Begin-----------------------------------------------------------------
        
        call RewindBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ReadListFilesFromBlock - ModuleFillMatrix - ERR10'        

        call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,                       &
                                            Begin_files,                                &
                                            End_files,                                  &
                                            BlockInBlockFound,                          &
                                            FirstLine, LastLine,                        &
                                            STAT = STAT_CALL)                                        
        if (STAT_CALL /= SUCCESS_) stop 'ReadListFilesFromBlock - ModuleFillMatrix - ERR20'
            
        if (BlockInBlockFound) then
                
            FileNumber = LastLine - FirstLine - 1
            
            allocate(CurrentHDF%FileNameList(FileNumber))


            do n = 1, FileNumber

                call GetData(CurrentHDF%FileNameList(n),                                &
                                Me%ObjEnterData,  iflag, Buffer_Line  = FirstLine + n,  &
                                STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadListFilesFromBlock - ModuleFillMatrix - ERR30'

            enddo
        

            call RewindBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'ReadListFilesFromBlock - ModuleFillMatrix - ERR40'
            
        endif            

    end subroutine ReadListFilesFromBlock
                
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    subroutine ReadListFilesFromBlockInBlock(CurrentHDF, ClientNumber)

        !Arguments-------------------------------------------------------------
        type(T_Field4D), pointer                        :: CurrentHDF
        integer                                         :: ClientNumber 

        !Local-----------------------------------------------------------------
        logical                                         :: BlockInBlockInBlockFound
        integer                                         :: STAT_CALL, FirstLine, LastLine, iflag
        integer                                         :: n, FileNumber 
        !Begin-----------------------------------------------------------------
        
        call RewindBlockInBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ReadListFilesFromBlockInBlock - ModuleFillMatrix - ERR10'        

        call ExtractBlockFromBlockFromBlock(Me%ObjEnterData, ClientNumber,              &
                                            Begin_files_2,                              &
                                            End_files_2,                                &
                                            BlockInBlockInBlockFound,                   &
                                            FirstLine, LastLine,                        &
                                            STAT = STAT_CALL)                                        
        if (STAT_CALL /= SUCCESS_) stop 'ReadListFilesFromBlockInBlock - ModuleFillMatrix - ERR20'
            
        if (BlockInBlockInBlockFound) then
                
            FileNumber = LastLine - FirstLine - 1
            
            allocate(CurrentHDF%FileNameList(FileNumber))


            do n = 1, FileNumber

                call GetData(CurrentHDF%FileNameList(n),                                &
                                Me%ObjEnterData,  iflag, Buffer_Line  = FirstLine + n,  &
                                STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadListFilesFromBlockInBlock - ModuleFillMatrix - ERR30'

            enddo
        

            call RewindBlockInBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'ReadListFilesFromBlockInBlock - ModuleFillMatrix - ERR40'
            
        endif            

    end subroutine ReadListFilesFromBlockInBlock
                
    !--------------------------------------------------------------------------
    
    
    subroutine BuildField4D(ExtractType, ClientID, PointsToFill2D, PointsToFill3D, CurrentHDF)    
    
         !Arguments------------------------------------------------------------
        integer                                         :: ExtractType, ClientID
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D 
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D  
        type(T_Field4D)                                 :: CurrentHDF
        !Local-----------------------------------------------------------------
        real                                            :: LatDefault, LongDefault
        type (T_Size2D)                                 :: HaloMap
        logical                                         :: MasterOrSlave
        integer                                         :: STAT_CALL
        !Begin-----------------------------------------------------------------           
        
ifSI:   if (CurrentHDF%SpatialInterpolON) then

            call ConstructField4DInterpol(ExtractType, ClientID, PointsToFill2D, PointsToFill3D, CurrentHDF)

        else ifSI
        
            call GetDDecompParameters(HorizontalGridID = Me%ObjHorizontalGrid,          &
                                      MasterOrSlave    = MasterOrSlave,                 &
                                      STAT             = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'BuildField4D - ModuleFillMatrix - ERR10'        
        
        
ifMS:       if (MasterOrSlave) then
            
                call GetDDecompWorkSize2D(HorizontalGridID = Me%ObjHorizontalGrid,      &
                                          WorkSize         = HaloMap,                   &
                                          STAT             = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'BuildField4D - ModuleFillMatrix - ERR20'
                                                      
                write(*,*) 'With domain decomposition - ILB,IUB, JLB, JUB',             &
                            HaloMap%ILB,HaloMap%IUB, HaloMap%JLB, HaloMap%JUB
                
            else ifMS
            
                if(Me%Dim == Dim2D)then
                    HaloMap     = Me%WorkSize2D
                else                    
                    HaloMap%ILB = Me%WorkSize3D%ILB
                    HaloMap%IUB = Me%WorkSize3D%IUB
                    HaloMap%JLB = Me%WorkSize3D%JLB
                    HaloMap%JUB = Me%WorkSize3D%JUB
                endif
                
                write(*,*) 'No domain decomposition - ILB,IUB, JLB, JUB',               &
                            HaloMap%ILB,HaloMap%IUB, HaloMap%JLB, HaloMap%JUB
            
            endif ifMS                

            call GetLatitudeLongitude(Me%ObjHorizontalGrid, Latitude  = LatDefault,     &
                                                            Longitude = LongDefault,    & 
                                                            STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'BuildField4D - ModuleFillMatrix - ERR30'

            call ConstructField4D(Field4DID         = CurrentHDF%ObjField4D,            &
                                  EnterDataID       = Me%ObjEnterData,                  &
                                  ExtractType       = ExtractType,                      &
                                  FileName          = CurrentHDF%FileName,              &
                                  FieldName         = CurrentHDF%FieldName,             &
                                  TimeID            = Me%ObjTime,                       &   
                                  MaskDim           = Me%Dim,                           &
                                  LatReference      = LatDefault,                       &
                                  LonReference      = LongDefault,                      & 
                                  WindowLimitsJI    = HaloMap,                          &
                                  Extrapolate       = .false.,                          &    
                                  ExtrapolateMethod = ExtrapolAverage_,                 &
                                  PropertyID        = Me%PropertyID,                    &                                  
                                  ClientID          = ClientID,                         &
                                  FileNameList      = CurrentHDF%FileNameList,          &
                                  CheckHDF5_File    = Me%CheckHDF5_File,                &
                                  STAT              = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'BuildField4D - ModuleFillMatrix - ERR40'
        
        endif ifSI
        
        call GetField4DNumberOfInstants(CurrentHDF%ObjField4D, CurrentHDF%NumberOfInstants, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'BuildField4D - ModuleFillMatrix - ERR50'

    end subroutine BuildField4D

    !-----------------------------------------------------------------------------------
    
    subroutine DTMethod_1_PrevNext3D(Now, PointsToFill3D, CurrentHDF)    
    
         !Arguments------------------------------------------------------------
        type(T_Time)                                    :: Now
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D   
        type(T_Field4D)                                 :: CurrentHDF
        !Local-----------------------------------------------------------------
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB, i, j, k
        !Begin-----------------------------------------------------------------           
    
        ILB = Me%Size3D%ILB
        IUB = Me%Size3D%IUB
        JLB = Me%Size3D%JLB
        JUB = Me%Size3D%JUB
        KLB = Me%Size3D%KLB
        KUB = Me%Size3D%KUB

        allocate(CurrentHDF%PreviousField3D (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(CurrentHDF%NextField3D     (ILB:IUB, JLB:JUB, KLB:KUB))

        CurrentHDF%PreviousField3D(:,:,:) = FillValueReal
        CurrentHDF%NextField3D    (:,:,:) = FillValueReal

        call ReadHDF5Values3D(CurrentHDF%PreviousInstant, CurrentHDF%PreviousField3D, CurrentHDF)
        call ReadHDF5Values3D(CurrentHDF%NextInstant,     CurrentHDF%NextField3D, CurrentHDF    )

        !limit maximum values
        do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
        do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
        do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB

#ifndef _NOT_IEEE_ARITHMETIC
            if (ieee_is_nan (CurrentHDF%PreviousField3D(i,j,k)))                    &
                CurrentHDF%PreviousField3D(i,j,k) = FillValueReal 
#endif
        
            if (abs(CurrentHDF%PreviousField3D(i,j,k)) > abs(FillValueReal))        &
                CurrentHDF%PreviousField3D(i,j,k) = FillValueReal

#ifndef _NOT_IEEE_ARITHMETIC
            if (ieee_is_nan (CurrentHDF%NextField3D    (i,j,k)))                    &
                CurrentHDF%NextField3D(i,j,k) = FillValueReal 
#endif
            
            if (abs(CurrentHDF%NextField3D    (i,j,k)) > abs(FillValueReal))        &
                CurrentHDF%NextField3D        (i,j,k) = FillValueReal
        enddo
        enddo
        enddo                


        if (CurrentHDF%PreviousInstant /= CurrentHDF%NextInstant) then
        
           if (Me%PreviousInstantValues) then
            
                Me%Matrix3D = CurrentHDF%PreviousField3D

            else             
            
                if (Me%PropertyID%IsAngle) then                           

                    call InterpolateAngle3DInTime (ActualTime       = Now,                         &
                                                   Size             = Me%WorkSize3D,               &
                                                   Time1            = CurrentHDF%PreviousTime,     &
                                                   Matrix1          = CurrentHDF%PreviousField3D,  &
                                                   Time2            = CurrentHDF%NextTime,         &
                                                   Matrix2          = CurrentHDF%NextField3D,      &
                                                   MatrixOut        = Me%Matrix3D,                 &
                                                   PointsToFill3D   = PointsToFill3D)
                else
                    call InterpolateMatrix3DInTime(ActualTime       = Now,                         &
                                                   Size             = Me%WorkSize3D,               &
                                                   Time1            = CurrentHDF%PreviousTime,     &
                                                   Matrix1          = CurrentHDF%PreviousField3D,  &
                                                   Time2            = CurrentHDF%NextTime,         &
                                                   Matrix2          = CurrentHDF%NextField3D,      &
                                                   MatrixOut        = Me%Matrix3D,                 &
                                                   PointsToFill3D   = PointsToFill3D)
                
                endif                                                   
            endif
            
        else

            !Prev and next are equal (last instant?)
            Me%Matrix3D(:,:,:)  = CurrentHDF%NextField3D(:,:,:)

        endif
                
    end subroutine DTMethod_1_PrevNext3D

    !-----------------------------------------------------------------------------------        

    subroutine DTMethod_1_PrevNext2D(Now, PointsToFill2D, CurrentHDF)    
    
         !Arguments------------------------------------------------------------
        type(T_Time)                                    :: Now
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D  
        type(T_Field4D)                                 :: CurrentHDF
        !Local-----------------------------------------------------------------
        integer                                         :: ILB, IUB, JLB, JUB, i, j
        !Begin-----------------------------------------------------------------          
        
        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        allocate(CurrentHDF%PreviousField2D (ILB:IUB, JLB:JUB))
        allocate(CurrentHDF%NextField2D     (ILB:IUB, JLB:JUB))

        CurrentHDF%PreviousField2D(:,:) = FillValueReal
        CurrentHDF%NextField2D    (:,:) = FillValueReal

        call ReadHDF5Values2D(CurrentHDF%PreviousInstant, CurrentHDF%PreviousField2D, CurrentHDF)
        call ReadHDF5Values2D(CurrentHDF%NextInstant,     CurrentHDF%NextField2D,     CurrentHDF)
        
        !limit maximum values
        do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
        do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
        
#ifndef _NOT_IEEE_ARITHMETIC
            if (ieee_is_nan (CurrentHDF%PreviousField2D(i,j)))                      &
                CurrentHDF%PreviousField2D(i,j) = FillValueReal
#endif
        
            if (abs(CurrentHDF%PreviousField2D   (i,j)) > abs(FillValueReal))       &
                CurrentHDF%PreviousField2D(i,j) = FillValueReal

#ifndef _NOT_IEEE_ARITHMETIC
            if (ieee_is_nan (CurrentHDF%NextField2D    (i,j)))                      &
                CurrentHDF%NextField2D(i,j)     = FillValueReal                
#endif
            
            if (abs(CurrentHDF%NextField2D       (i,j)) > abs(FillValueReal))       &
                CurrentHDF%NextField2D(i,j)     = FillValueReal
        enddo
        enddo
                        
        if (CurrentHDF%PreviousInstant /= CurrentHDF%NextInstant) then
        
            !Interpolates the two matrixes in time
            if (Me%PropertyID%IsAngle) then        
                call InterpolateAngle2DInTime  (ActualTime       = Now,                         &
                                               Size             = Me%WorkSize2D,               &
                                               Time1            = CurrentHDF%PreviousTime,     &
                                               Matrix1          = CurrentHDF%PreviousField2D,  &
                                               Time2            = CurrentHDF%NextTime,         &
                                               Matrix2          = CurrentHDF%NextField2D,      &
                                               MatrixOut        = Me%Matrix2D,                 &
                                               PointsToFill2D   = PointsToFill2D)
            else        
                call InterpolateMatrix2DInTime(ActualTime       = Now,                         &
                                               Size             = Me%WorkSize2D,               &
                                               Time1            = CurrentHDF%PreviousTime,     &
                                               Matrix1          = CurrentHDF%PreviousField2D,  &
                                               Time2            = CurrentHDF%NextTime,         &
                                               Matrix2          = CurrentHDF%NextField2D,      &
                                               MatrixOut        = Me%Matrix2D,                 &
                                               PointsToFill2D   = PointsToFill2D)
            endif                                                                      
        else

            Me%Matrix2D(:,:)  = CurrentHDF%PreviousField2D(:,:)

        endif

    end subroutine DTMethod_1_PrevNext2D

    !-----------------------------------------------------------------------------------    
    
    subroutine ConstructHDFOneInstant(Now, CurrentHDF)    
    
         !Arguments------------------------------------------------------------
        type(T_Time)                                    :: Now
        type(T_Field4D)                                 :: CurrentHDF    
        !Local-----------------------------------------------------------------

        real                                            :: Year, Month, Day, Hour, Minute, Second
        !Begin-----------------------------------------------------------------          
    
        CurrentHDF%PreviousInstant  = 1
        CurrentHDF%NextInstant      = CurrentHDF%PreviousInstant

        CurrentHDF%PreviousTime     = HDF5TimeInstant(CurrentHDF%PreviousInstant, CurrentHDF)
        
        !if(CurrentHDF%GenericYear)then
            !call SetHDFGenericYear(CurrentHDF%PreviousTime, Now)
        !endif
                    
        if (CurrentHDF%CyclicTimeON) then        
            call SetHDFCyclicDates(TimeInstant  = CurrentHDF%PreviousTime,                  &
                                   RefTime      = Now,                                      &
                                   CurrentHDF   = CurrentHDF)
                                   
            CurrentHDF%StartTime = Me%BeginTime
            CurrentHDF%EndTime   = Me%EndTime
        endif

        CurrentHDF%NextTime         = CurrentHDF%PreviousTime

        call ExtractDate(CurrentHDF%PreviousTime, Year, Month, Day, Hour, Minute, Second)

        write(*,*)
        write(*,*)trim(CurrentHDF%FieldName)//' is being read from HDF file.' 
        write(*,*)'Time instant: ', Year, Month, Day, Hour, Minute, Second
        write(*,*)'ConstructHDFOneInstant - ModuleFillMatrix - WRN10'
    
    end subroutine ConstructHDFOneInstant 

    !-----------------------------------------------------------------------------------
   

    subroutine ConstructHDFPredictDTMethod1(Now, CurrentHDF)
    
         !Arguments------------------------------------------------------------
        type(T_Time)                                    :: Now
        type(T_Field4D)                                 :: CurrentHDF
        !Local-----------------------------------------------------------------
        logical                                         :: FoundSecondInstant
        real                                            :: StartTimeYear, EndTimeYear
        !Begin-----------------------------------------------------------------            

        if (Me%Backtracking) then
            CurrentHDF%PreviousInstant  = CurrentHDF%NumberOfInstants
            CurrentHDF%NextInstant      = CurrentHDF%PreviousInstant               
        else
            CurrentHDF%PreviousInstant  = 1
            CurrentHDF%NextInstant      = CurrentHDF%PreviousInstant
        endif
        
        CurrentHDF%PreviousTime         = HDF5TimeInstant(CurrentHDF%PreviousInstant, CurrentHDF)
        
        !if(CurrentHDF%GenericYear)then
        !    call SetHDFGenericYear(CurrentHDF%PreviousTime, RefTime = Now)
        !endif
        
        call SetHDFCyclicDates(TimeInstant  = CurrentHDF%PreviousTime,                  &
                               RefTime      = Now,                                      &
                               CyclicTimeON = CurrentHDF%CyclicTimeON,                  &
                               CurrentHDF   = CurrentHDF)
                               
         if (CurrentHDF%CyclicTimeON .and. .not. CurrentHDF%GenericYear) then
            CurrentHDF%StartTime = Me%BeginTime
            CurrentHDF%EndTime   = Me%EndTime
         endif                               

ib:     if (Me%BackTracking) then  
                
            if(CurrentHDF%PreviousTime .lt. Now)then
                write(*,*)
                write(*,*)'----------Backtracking mode-----------'
                write(*,*)'Could not read solution from HDF5 file'
                write(*,*)'Last file instant greater than current time'
                write(*,*)'Matrix name: '//trim(CurrentHDF%FieldName)
                stop      'ConstructHDFPredictDTMethod1 - ModuleFillMatrix - ERR10'
            end if

            if(Me%TimeEvolution .ne. None)then                        
                if(CurrentHDF%StartTime .gt. Me%BeginTime)then
                    write(*,*)
                    write(*,*)'----------Backtracking mode-----------'                                
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'First instant in file lower than simulation starting time'
                    write(*,*)'Matrix name: '//trim(CurrentHDF%FieldName)
                    stop      'ConstructHDFPredictDTMethod1 - ModuleFillMatrix - ERR20'
                end if
            endif
            
        else   ib
        
            if(CurrentHDF%PreviousTime .gt. Now)then
                write(*,*)
                write(*,*)'Could not read solution from HDF5 file'
                write(*,*)'First file instant greater than current time'
                write(*,*)'Matrix name: '//trim(CurrentHDF%FieldName)
                stop      'ConstructHDFPredictDTMethod1 - ModuleFillMatrix - ERR30'
            end if

            if(CurrentHDF%GenericYear)then
            
                call SetHDFGenericYear(CurrentHDF%EndTime,   RefTime = Me%EndTime)
                call SetHDFGenericYear(CurrentHDF%StartTime, RefTime = Me%BeginTime)
                
                call ExtractDate(CurrentHDF%StartTime, Year = StartTimeYear)
                call ExtractDate(CurrentHDF%EndTime,   Year = EndTimeYear  )
                
                if(StartTimeYear .ne. EndTimeYear)then
                    write(*,*)
                    write(*,*)'When using a generic year HDF5 file'
                    write(*,*)'The year of the start time has to be the same as'
                    write(*,*)'the year of the end time'
                    write(*,*)'Matrix name: '//trim(CurrentHDF%FieldName)
                    stop      'ConstructHDFPredictDTMethod1 - ModuleFillMatrix - ERR40'
                endif
                
            endif


            if(Me%TimeEvolution .ne. None)then
                if(CurrentHDF%EndTime .lt. Me%EndTime)then
                    write(*,*)
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'Last instant in file lower than simulation ending time'
                    write(*,*)'Matrix name: '//trim(CurrentHDF%FieldName)
                    stop      'ConstructHDFPredictDTMethod1 - ModuleFillMatrix - ERR50'
                end if
            end if
        endif ib

        FoundSecondInstant = .false.
    
        !if number of instants greater than 1 then 
        !find first and second instants
d2:      do while(.not. FoundSecondInstant)

            CurrentHDF%PreviousInstant  = CurrentHDF%NextInstant
            if (Me%Backtracking) then
                CurrentHDF%NextInstant      = CurrentHDF%NextInstant - 1
            else                
                CurrentHDF%NextInstant      = CurrentHDF%NextInstant + 1
            endif

            !if (CurrentHDF%CyclicTimeON .and. CurrentHDF%NextInstant .gt. CurrentHDF%NumberOfInstants) then
            !    CurrentHDF%NextInstant  = 1
            !end if


            CurrentHDF%NextTime         = HDF5TimeInstant(CurrentHDF%NextInstant, CurrentHDF)
            
            
            !if(CurrentHDF%GenericYear)then
            !    call SetHDFGenericYear(CurrentHDF%NextTime, Now)
            !endif
            
            if (CurrentHDF%CyclicTimeON) then
                call SetHDFCyclicDates(TimeInstant  = CurrentHDF%NextTime,                  &
                                       RefTime      = Now,                                  &
                                       CurrentHDF   = CurrentHDF)
            endif
            
            if (Me%Backtracking) then
                if(CurrentHDF%PreviousTime .ge. Now .and. CurrentHDF%NextTime .le. Now) then
                    FoundSecondInstant  = .true.
                    exit
                end if
            else
                if(CurrentHDF%PreviousTime .le. Now .and. CurrentHDF%NextTime .ge. Now) then
                    FoundSecondInstant  = .true.
                    exit
                end if
            endif
            CurrentHDF%PreviousTime = CurrentHDF%NextTime

            if (Me%Backtracking) then
                if(CurrentHDF%NextInstant .lt. 1 .and. .not. CurrentHDF%CyclicTimeON) then
                    write(*,*)
                    write(*,*)'----------Backtracking mode-----------------'
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'Could not find second instant in file'
                    write(*,*)'Matrix name: '//trim(CurrentHDF%FieldName)
                    stop      'ConstructHDFPredictDTMethod1 - ModuleFillMatrix - ERR60'
                end if
            
            else
                if(CurrentHDF%NextInstant .gt. CurrentHDF%NumberOfInstants .and. .not. CurrentHDF%CyclicTimeON) then
                    write(*,*)
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'Could not find second instant in file'
                    write(*,*)'Matrix name: '//trim(CurrentHDF%FieldName)
                    stop      'ConstructHDFPredictDTMethod1 - ModuleFillMatrix - ERR70'
                end if
            endif
        end do d2

    end subroutine ConstructHDFPredictDTMethod1

    !--------------------------------------------------------------------------

            

    subroutine ConstructHDFPredictDTMethod2(Now, PointsToFill3D, PointsToFill2D, CurrentHDF)
    
         !Arguments------------------------------------------------------------
        type(T_Time)                                    :: Now
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D
        type(T_Field4D)                                 :: CurrentHDF
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------            

        !This methodology do NOT works with CYCLICTIME or BACKTRACKING (needs revision)
        if (Me%Backtracking) then
            stop 'ConstructHDFPredictDTMethod2 - ModuleFillMatrix - ERR10'
        endif

        CurrentHDF%PreviousInstant = 1
        CurrentHDF%PreviousTime    = HDF5TimeInstant(CurrentHDF%PreviousInstant, CurrentHDF)
        
        if(CurrentHDF%PreviousTime .gt. Now)then
            write(*,*)
            write(*,*)'Could not read solution from HDF5 file'
            write(*,*)'First file instant greater than current time'
            write(*,*)'Matrix name: '//trim(CurrentHDF%FieldName)
            stop      'ConstructHDFPredictDTMethod2 - ModuleFillMatrix - ERR20'
        end if
        
        if(Me%TimeEvolution .ne. None)then
            if(CurrentHDF%EndTime .lt. Me%EndTime)then
                write(*,*)
                write(*,*)'Could not read solution from HDF5 file'
                write(*,*)'Last instant in file lower than simulation ending time'
                write(*,*)'Matrix name: '//trim(CurrentHDF%FieldName)
                stop      'ConstructHDFPredictDTMethod2 - ModuleFillMatrix - ERR30'
            end if
        end if                    
        
        CurrentHDF%NextInstant = 1
        CurrentHDF%NextTime    = CurrentHDF%PreviousTime
        
        call ActualizeHDFValues (CurrentHDF, CurrentHDF%PreviousInstant, CurrentHDF%PreviousField2D)
        call ActualizeHDFValues (CurrentHDF, CurrentHDF%NextInstant, CurrentHDF%NextField2D)
        
        CurrentHDF%NextEventStart        = CurrentHDF%PreviousTime
        CurrentHDF%NextEventEnd          = CurrentHDF%PreviousTime
        
        CurrentHDF%DTForNextEvent        = 0.0               
        
        if (Me%Dim == Dim2D) then
            CurrentHDF%NextValueForDTPred = maxval(CurrentHDF%NextField2D)
            call ModifyHDFInput2D (PointsToFill2D, CurrentHDF) 
        else
            CurrentHDF%NextValueForDTPred = maxval(CurrentHDF%NextField3D)
            call ModifyHDFInput3D (PointsToFill3D, CurrentHDF)
        endif
                    
    end subroutine ConstructHDFPredictDTMethod2
    
    !----------------------------------------------------------------------------------
    
    
    subroutine ConstructField4DInterpol(ExtractType, ClientID, PointsToFill2D, PointsToFill3D, CurrentHDF)
        !Arguments----------------------------------------------------------------------
        integer                                         :: ExtractType, ClientID
        integer,  dimension(:,:,:),   pointer, optional :: PointsToFill3D
        integer,  dimension(:,:),     pointer, optional :: PointsToFill2D
        type(T_Field4D)                                 :: CurrentHDF

        !Local--------------------------------------------------------------------------
        real,       dimension(:,:),     pointer         :: CoordX, CoordY
        real, dimension(1:2,1:2)                        :: WindowLimitsXY
        real                                            :: West, East, South, North  
        real                                            :: LatDefault, LongDefault
        integer                                         :: STAT_CALL, i, j, k, icount, NCells
        real, dimension(4)                              :: Aux4
        integer                                         :: iflag, ObjHorizontalGridAux
        character(len=PathLength)                       :: BathymetryFile
        
        !Begin--------------------------------------------------------------------------      

        Aux4 (:) = FillValueReal

        !West, East, South, North
        call GetData(Aux4,                                                              &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'BORDER_LIMITS',                                    &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR10'
        
        if (iflag < 4) then
            
            call ReadFileName('IN_BATIM', BathymetryFile, "Bathymetry File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR20'
            
            ObjHorizontalGridAux = 0
    
            !Entire grid
            call ConstructHorizontalGrid(HorizontalGridID = ObjHorizontalGridAux,       &
                                         DataFile         = BathymetryFile,             &
                                         STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR30'
            
 
            call GetGridBorderLimits(ObjHorizontalGridAux, West, East, South, North, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR40'

            call KillHorizontalGrid(HorizontalGridID = ObjHorizontalGridAux,            &
                                    STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR50'

        elseif (iflag == 4) then
        
            West = Aux4(1); East = Aux4(2); South = Aux4(3); North = Aux4(4);
        
        else
        
            stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR60'
        
        endif

        write(*,*) 'Border limits', West, East, South, North

        WindowLimitsXY(2,1) = South
        WindowLimitsXY(2,2) = North
        WindowLimitsXY(1,1) = West
        WindowLimitsXY(1,2) = East
        
        if (Me%Dim == Dim2D) then
            
            icount = 0
            do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB        
                if (PointsToFill2D(i,j) == WaterPoint) icount = icount + 1
            enddo
            enddo
            
            Ncells        = icount
            CurrentHDF%Ncells = Ncells
        
        else
        
            icount = 0
            
            if (CurrentHDF%From2Dto3D) then            

                do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB        
                    if (PointsToFill3D(i,j,Me%WorkSize3D%KUB) == WaterPoint) icount = icount + 1
                enddo
                enddo
                
                Ncells        = icount
                CurrentHDF%Ncells = Ncells                
            
            else
            
                do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
                do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB        
                    if (PointsToFill3D(i,j,k) == WaterPoint) icount = icount + 1
                enddo
                enddo
                enddo
                
                Ncells        = icount
                CurrentHDF%Ncells = Ncells                
                
                allocate(CurrentHDF%Z(1:NCells))
         
                CurrentHDF%Z(1:NCells) = FillValueReal                
                
            endif                
            

        

                                                     
        endif                    
            
        allocate(CurrentHDF%X(1:NCells), CurrentHDF%Y(1:NCells), CurrentHDF%Prop(1:NCells), CurrentHDF%NoData(1:NCells))
        
        CurrentHDF%X     (1:NCells) = FillValueReal
        CurrentHDF%Y     (1:NCells) = FillValueReal
        CurrentHDF%Prop  (1:NCells) = FillValueReal
        CurrentHDF%NoData(1:NCells) = .true.
        
        call GetZCoordinates(Me%ObjHorizontalGrid, CoordX, CoordY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR70'
        
        if (Me%Dim == Dim2D) then
            
            icount = 0
            
            do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB        
            
                if (PointsToFill2D(i,j) == WaterPoint) then

                    icount           = icount + 1
                    CurrentHDF%X(icount) = CoordX(i, j)
                    CurrentHDF%Y(icount) = CoordY(i, j)

                endif                    
                
            enddo
            enddo
            
        else
        
            icount = 0
            
            if (CurrentHDF%From2Dto3D) then            

                do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB        
                
                    if (PointsToFill3D(i,j,Me%WorkSize3D%KUB) == WaterPoint) then
                    
                        icount           = icount + 1
                        CurrentHDF%X(icount) = CoordX(i, j)
                        CurrentHDF%Y(icount) = CoordY(i, j)
                
                    endif                    
                
                enddo
                enddo   
            
            else
            
                do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
                do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB        
                
                    if (PointsToFill3D(i,j,k) == WaterPoint) then
                    
                        icount           = icount + 1
                        CurrentHDF%X(icount) = CoordX(i, j)
                        CurrentHDF%Y(icount) = CoordY(i, j)
                
                    endif                    
                
                enddo
                enddo   
                enddo     
                
            endif                
            
        endif       
    
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR80'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR90'
        
        call GetLatitudeLongitude(Me%ObjHorizontalGrid, Latitude  = LatDefault,         &
                                                        Longitude = LongDefault,        & 
                                                        STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR100'
        
        if (CurrentHDF%InterpolOnlyVertically) then
        

            call ConstructField4D(Field4DID         = CurrentHDF%ObjField4D,            &
                                  EnterDataID       = Me%ObjEnterData,                  &
                                  ExtractType       = ExtractType,                      &
                                  FileName          = CurrentHDF%FileName,              &
                                  FieldName         = CurrentHDF%FieldName,             &
                                  TimeID            = Me%ObjTime,                       &   
                                  MaskDim           = Me%Dim,                           &
                                  HorizontalGridID  = Me%ObjHorizontalGrid,             &                                
                                  Extrapolate       = .true.,                           &
                                  ExtrapolateMethod = CurrentHDF%ExtrapolateMethod,     &
                                  PropertyID        = Me%PropertyID,                    &
                                  ClientID          = ClientID,                         &
                                  FileNameList      = CurrentHDF%FileNameList,          &
                                  CheckHDF5_File    = Me%CheckHDF5_File,                &
                                  STAT              = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR110'        
        
        else
        
            call ConstructField4D(Field4DID         = CurrentHDF%ObjField4D,            &
                                  EnterDataID       = Me%ObjEnterData,                  &
                                  ExtractType       = ExtractType,                      &
                                  FileName          = CurrentHDF%FileName,              &
                                  FieldName         = CurrentHDF%FieldName,             &
                                  TimeID            = Me%ObjTime,                       &   
                                  MaskDim           = Me%Dim,                           &
                                  LatReference      = LatDefault,                       &
                                  LonReference      = LongDefault,                      & 
                                  WindowLimitsXY    = WindowLimitsXY,                   &
                                  Extrapolate       = .true.,                           &
                                  ExtrapolateMethod = CurrentHDF%ExtrapolateMethod,     &
                                  PropertyID        = Me%PropertyID,                    &
                                  ClientID          = ClientID,                         &
                                  FileNameList      = CurrentHDF%FileNameList,          &
                                  CheckHDF5_File    = Me%CheckHDF5_File,                &
                                  STAT              = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR120'
        endif
    
    end subroutine ConstructField4DInterpol

    
   !----------------------------------------------------------------------------

    subroutine Generic4thDimension(ExtractType, CurrentHDF)

        !Arguments-------------------------------------------------------------
        integer                            :: ExtractType
        type(T_Field4D)                    :: CurrentHDF
        !Local-----------------------------------------------------------------
        character(len = PathLength)        :: Filename
        integer                            :: STAT_CALL, iflag

        !----------------------------------------------------------------------

        call GetData(CurrentHDF%Generic4D%ReadFromTimeSerie,                                &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = '4D_TIME_SERIE',                                    &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'Generic4thDimension - ModuleFillMatrix - ERR10'


        if (CurrentHDF%Generic4D%ReadFromTimeSerie) then

            call GetData(Filename, Me%ObjEnterData, iflag,                              &
                         keyword        = 'GENERIC_4D_FILENAME',                        &  
                         SearchType     = ExtractType,                                  &
                         ClientModule   = 'ModuleFillMatrix',                           &
                         default        = "******.***",                                 &
                         STAT           = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'Generic4thDimension  - ModuleFillMatrix - ERR20'
            if (iflag == 0) stop 'Generic4thDimension  - ModuleFillMatrix - ERR30'

            call GetData(CurrentHDF%Generic4D%TimeSerieColumn, Me%ObjEnterData, iflag,  &
                         keyword        = 'TIME_SERIE_COLUMN',                          &  
                         SearchType     = ExtractType,                                  &
                         ClientModule   = 'ModuleFillMatrix',                           &
                         default        = 2,                                            &
                         STAT           = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'Generic4thDimension  - ModuleFillMatrix - ERR40'

            !Starts Time Serie
            call StartTimeSerieInput(CurrentHDF%Generic4D%ObjTimeSerie,     &
                                     FileName,                              &
                                     CheckDates =.false.,                   &
                                     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Generic4thDimension - ModuleFillMatrix - ERR50'

        endif


     end subroutine Generic4thDimension
   !----------------------------------------------------------------------------


    type(T_Time) function HDF5TimeInstant(Instant, CurrentHDF)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        type(T_Field4D)                         :: CurrentHDF

        !Local-----------------------------------------------------------------
!        type(T_Time)                            :: TimeInstant
        real,    dimension(:), pointer          :: TimeVector
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------

if4D:   if (CurrentHDF%Field4D) then      
        
            HDF5TimeInstant = GetField4DInstant(Field4DID = CurrentHDF%ObjField4D,          &
                                                Instant   = Instant,                    &
                                                STAT      =  STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleFillMatrix - ERR10'
        
        else if4D
        
            call HDF5SetLimits  (CurrentHDF%ObjHDF5, 1, 6, STAT = STAT_CALL)

            allocate(TimeVector(6))

            call HDF5ReadData   (HDF5ID         = CurrentHDF%ObjHDF5,                           &
                                 GroupName      = "/Time",                                  &
                                 Name           = "Time",                                   &
                                 Array1D        = TimeVector,                               &
                                 OutputNumber   = Instant,                                  &
                                 STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleFillMatrix - ERR20'
            
            call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                          Day      = TimeVector(3), Hour   = TimeVector(4), &
                                          Minute   = TimeVector(5), Second = TimeVector(6))

                                         
            deallocate(TimeVector)
            

        endif if4D


    end function HDF5TimeInstant

    
    !--------------------------------------------------------------------------

    subroutine CheckCyclicMonths(TimeInstant, RefTime, PreviousTime, CyclicTimeON)

        !Arguments-------------------------------------------------------------
        type(T_Time)                            :: TimeInstant
        type(T_Time), optional                  :: RefTime, PreviousTime
        logical     , optional                  :: CyclicTimeON

        !Local-----------------------------------------------------------------
        real                                    :: Year, YearX, Month, MonthX, Day, Hour, Minute, Second

        !Begin-----------------------------------------------------------------
        
        !Cyclic months

        call ExtractDate(TimeInstant, Year = Year, Month  = Month,  Day    = Day,                &
                                      Hour = Hour, Minute = Minute, Second = Second)
        
        !!David
        !!This does not make sense. get info from hdf but this routine is called from profileTimeSerie and ProfileTSDefault
        if(Me%HDF%GenericYear) Year = CyclicTime
        
        if (Year == CyclicTime) then

            if (present(CyclicTimeON)) CyclicTimeON = .true.

            if (present(RefTime) .and. present(PreviousTime)) then 
                stop 'CheckCyclicMonths - ModuleFillMatrix - ERR10'
            endif 

            if      (present(RefTime)) then 

                call ExtractDate(RefTime, Year = Year)

            else if (present(PreviousTime)) then
             
                call ExtractDate(PreviousTime, Year = YearX, Month  = MonthX)
                
                Year = YearX
                if (MonthX > Month) then
                    Year = Year + 1
                endif

            endif

        endif

        call SetDate(TimeInstant, Year = Year, Month  = Month,  Day    = Day,     &
                                  Hour = Hour, Minute = Minute, Second = Second)


    end subroutine CheckCyclicMonths
    
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    subroutine SetHDFCyclicDates(TimeInstant, RefTime, RestartCycle, CyclicTimeON, CurrentHDF)

        !Arguments-------------------------------------------------------------
        type(T_Time)                            :: TimeInstant
        type(T_Time)                            :: RefTime
        logical     , optional, intent(IN)      :: RestartCycle
        logical     , optional, intent(OUT)     :: CyclicTimeON
        type(T_Field4D)                         :: CurrentHDF

        !Local-----------------------------------------------------------------
        logical                                 :: RestartCycle_
        logical                                 :: CyclicTimeON_
        real, dimension(6)                      :: AuxTime, AuxTimeRef 
        integer                                 :: i

        !Begin-----------------------------------------------------------------
        
        call ExtractDate(TimeInstant, Year      = AuxTime(1),                           &
                                      Month     = AuxTime(2),                           &  
                                      Day       = AuxTime(3),                           &
                                      Hour      = AuxTime(4),                           &
                                      Minute    = AuxTime(5),                           &
                                      Second    = AuxTime(6))                                   
                                      
        if(CurrentHDF%GenericYear) AuxTime(1) = CyclicTime
            
        if (present(RestartCycle)) then
            RestartCycle_ = RestartCycle
        else
            RestartCycle_ = .false.                 
        endif            
        
        if (AuxTime(1) == CyclicTime) then
            CyclicTimeON_  = .true. 
        else
            CyclicTimeON_  = .false.         
        endif            

        if (present(CyclicTimeON)) CyclicTimeON = CyclicTimeON_

i1:     if (CyclicTimeON_) then

            call ExtractDate(RefTime, Year      = AuxTimeRef(1),                        &
                                      Month     = AuxTimeRef(2),                        &  
                                      Day       = AuxTimeRef(3),                        &
                                      Hour      = AuxTimeRef(4),                        &
                                      Minute    = AuxTimeRef(5),                        &
                                      Second    = AuxTimeRef(6))  
d1:         do i=1,6
                if (AuxTime(i) == CyclicTime) then
                    AuxTime(i) = AuxTimeRef(i)        
                else
                    if (RestartCycle_) then
                        if (AuxTimeRef(i) > AuxTime(i)) then
                            AuxTime(i-1) = AuxTime(i-1) + 1
                        endif
                    endif                                                
                    exit
                endif
            enddo d1

        endif i1

        call SetDate    (TimeInstant, Year      = AuxTime(1),                           &
                                      Month     = AuxTime(2),                           &  
                                      Day       = AuxTime(3),                           &
                                      Hour      = AuxTime(4),                           &
                                      Minute    = AuxTime(5),                           &
                                      Second    = AuxTime(6))  

    end subroutine SetHDFCyclicDates
    
    !--------------------------------------------------------------------------
    
    
    subroutine SetHDFGenericYear(TimeInstant, RefTime, AddYear)

        !Arguments-------------------------------------------------------------
        type(T_Time)                            :: TimeInstant
        type(T_Time)                            :: RefTime
        logical, optional                       :: AddYear

        !Local-----------------------------------------------------------------
        real                                    :: Year, Month, Day, Hour, Minute, Second
        logical                                 :: AddYear_

        !Begin-----------------------------------------------------------------

        call ExtractDate(TimeInstant, Year = Year, Month  = Month,  Day    = Day,       &
                                      Hour = Hour, Minute = Minute, Second = Second)
 
        call ExtractDate(RefTime, Year = Year)
        
        if(present(AddYear))then
            AddYear_ = AddYear
        else
            AddYear_ = .false.
        end if
        
        if(AddYear_)Year = Year + 1

        call SetDate(TimeInstant, Year = Year, Month  = Month,  Day    = Day,           &
                                  Hour = Hour, Minute = Minute, Second = Second)

    end subroutine SetHDFGenericYear
    
    !--------------------------------------------------------------------------

    real function HDF5Generic4DInstant(Instant, CurrentHDF)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        type(T_Field4D)                         :: CurrentHDF

        !Local-----------------------------------------------------------------
        real,   dimension(:), pointer           :: AuxVector
        type (T_Time)                           :: TimeInstant 
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------

if4D:   if (CurrentHDF%Field4D) then

            if (CurrentHDF%Generic4D%ReadFromTimeSerie) then

                TimeInstant =  GetField4DInstant (CurrentHDF%ObjField4D, Instant, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'HDF5Generic4DInstant - ModuleFillMatrix - ERR10'                

                HDF5Generic4DInstant = TimeSerieValue(CurrentHDF%Generic4D%ObjTimeSerie,        &
                                                      TimeInstant,                              &
                                                      CurrentHDF%Generic4D%TimeSerieColumn) 

            else
            
                HDF5Generic4DInstant = GetField4DGeneric4DValue (Field4DID = CurrentHDF%ObjField4D, &
                                                                 Instant   = Instant,               &
                                                                 STAT      = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'HDF5Generic4DInstant - ModuleFillMatrix - ERR20'

            endif


        else                

            if (CurrentHDF%Generic4D%ReadFromTimeSerie) then

                TimeInstant = HDF5TimeInstant(Instant, CurrentHDF)

                HDF5Generic4DInstant = TimeSerieValue(CurrentHDF%Generic4D%ObjTimeSerie,   &
                                                      TimeInstant,                         &
                                                      CurrentHDF%Generic4D%TimeSerieColumn) 

            else
            
                call HDF5SetLimits  (CurrentHDF%ObjHDF5, 1, 1, STAT = STAT_CALL)

                allocate(AuxVector(1))

                call HDF5ReadData   (HDF5ID         = CurrentHDF%ObjHDF5,               &
                                     GroupName      = "/Generic4D",                     &
                                     Name           = "Generic4D",                      &
                                     Array1D        = AuxVector,                        &
                                     OutputNumber   = Instant,                          &
                                     STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    write(*,*) 'ObjHDF5=', CurrentHDF%ObjHDF5 
                    write(*,*) 'FileName=', trim(CurrentHDF%FileName)
                    write(*,*) 'STAT_CALL= ', STAT_CALL
                    stop 'HDF5Generic4DInstant - ModuleFillMatrix - ERR10'
                endif
                
                HDF5Generic4DInstant = AuxVector(1)
     
                deallocate(AuxVector)

            endif 
            

        endif  if4D          

    end function HDF5Generic4DInstant
    
    !--------------------------------------------------------------------------
    subroutine ReadHDF5Values2D (InstantIn, Field, CurrentHDF)

        !Arguments-------------------------------------------------------------
        integer                                 :: InstantIn
        real, dimension(:,:), pointer           :: Field
        type(T_Field4D)                         :: CurrentHDF
        !Local-----------------------------------------------------------------
        integer                                 :: Instant, STAT_CALL, Imax, Jmax, i, j
        type (T_Time)                           :: CurrentTime

        !Begin-----------------------------------------------------------------
        
        !Backtracking time inversion is also done in the ModuleField4D    
        if (Me%BackTracking .and. .not. CurrentHDF%Field4D) then  
            Instant = CurrentHDF%NumberOfInstants - InstantIn + 1
        else
            Instant = InstantIn
        endif
if4D:   if (CurrentHDF%Field4D) then

            CurrentTime = GetField4DInstant (CurrentHDF%ObjField4D, Instant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR10'
            
            if (CurrentHDF%SpatialInterpolON) then
            
                call ModifyField4DInterpol(CurrentTime      = CurrentTime,              & 
                                           Matrix2D         = Field,                    &
                                           CurrentHDF       = CurrentHDF,               &
                                           Instant          = Instant)

            else
            
                call ModifyField4D(Field4DID        = CurrentHDF%ObjField4D,            &
                                   PropertyIDNumber = Me%PropertyID%IDNumber,           & 
                                   CurrentTime      = CurrentTime,                      &
                                   Matrix2D         = Field,                            &
                                   Instant          = Instant,                          &
                                   STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR10'

            endif                

        else if4D            
        
            call GetHDF5ArrayDimensions(CurrentHDF%ObjHDF5, trim(CurrentHDF%VGroupPath),        &
                              trim(CurrentHDF%FieldName), OutputNumber = Instant,           &
                              Imax = Imax, Jmax = Jmax, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR20'                               
            
            if ((Imax /= Me%WorkSize2D%IUB - Me%WorkSize2D%ILB + 1) .or.                &
                (Jmax /= Me%WorkSize2D%JUB - Me%WorkSize2D%JLB + 1)) then
                
                write (*,*) trim(CurrentHDF%VGroupPath)
                write (*,*) trim(CurrentHDF%FieldName)
                write (*,*) 'miss match between the HDF5 input file and model domain'
                stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR30'                                   

            endif

            call HDF5SetLimits  (CurrentHDF%ObjHDF5, Me%WorkSize2D%ILB, Me%WorkSize2D%IUB,      &
                                 Me%WorkSize2D%JLB, Me%WorkSize2D%JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR40'


            call HDF5ReadData(CurrentHDF%ObjHDF5, trim(CurrentHDF%VGroupPath),                      &
                              trim(CurrentHDF%FieldName),                                       &
                              Array2D = Field, OutputNumber = Instant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR50'
            
            if(CurrentHDF%HasMultiplyingFactor)then
                do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
                do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                    Field(i,j) = Field(i,j) * CurrentHDF%MultiplyingFactor
                enddo
                enddo
            end if

            if(CurrentHDF%HasAddingFactor)then
                do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
                do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                    Field(i,j) = Field(i,j) + CurrentHDF%AddingFactor
                enddo
                enddo
            end if            

        endif if4D



    end subroutine ReadHDF5Values2D
    
    
    !--------------------------------------------------------------------------

    
    subroutine ReadHDF5Values3D (Instant, Field, CurrentHDF)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        real, dimension(:,:,:), pointer         :: Field
        type(T_Field4D)                         :: CurrentHDF
        !Local-----------------------------------------------------------------
        real,    dimension(:,:,:), pointer      :: AuxField, SZZ
        integer, dimension(:,:,:), pointer      :: WaterPoints3D
        real                                    :: HT
        type (T_Time)                           :: CurrentTime
        integer                                 :: Imax, Jmax, Kmax, kbottom
        integer                                 :: STAT_CALL, i, j, k, ILB, IUB, JLB, JUB, KLB, KUB
        character(StringLength)                 :: DataSetVert
        logical                                 :: Exist1, Exist2

        !Begin-----------------------------------------------------------------


         ILB = Me%WorkSize3D%ILB
         IUB = Me%WorkSize3D%IUB
         JLB = Me%WorkSize3D%JLB
         JUB = Me%WorkSize3D%JUB

        if (CurrentHDF%From2Dto3D) then
            KLB = 1
            KUB = 1
            Kmax= 1
        else
            KLB = Me%WorkSize3D%KLB
            KUB = Me%WorkSize3D%KUB
        
            if (.not. CurrentHDF%Field4D) then
                CurrentHDF%ReadField3D => Field
            endif
            
        endif
        
        
if4D:   if (CurrentHDF%Field4D) then

            CurrentTime = GetField4DInstant (CurrentHDF%ObjField4D, Instant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR10'
            
            if (CurrentHDF%SpatialInterpolON) then
            
                call ModifyField4DInterpol(CurrentTime      = CurrentTime,              & 
                                           Matrix3D         = Field,                    &
                                           CurrentHDF       = CurrentHDF,               &
                                           Instant          = Instant)

            else            

                call ModifyField4D(Field4DID        = CurrentHDF%ObjField4D,            &
                                   PropertyIDNumber = Me%PropertyID%IDNumber,           & 
                                   CurrentTime      = CurrentTime,                      & 
                                   Matrix3D         = Field,                            &
                                   Instant          = Instant,                          &                                   
                                   STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR20'
                
            endif                

        else if4D            
        
            if (CurrentHDF%From2Dto3D) then
                call GetHDF5ArrayDimensions(CurrentHDF%ObjHDF5, trim(CurrentHDF%VGroupPath),    &
                                  trim(CurrentHDF%FieldName), OutputNumber = Instant,       &
                                  Imax = Imax, Jmax = Jmax, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR30'
            else
                call GetHDF5ArrayDimensions(CurrentHDF%ObjHDF5, trim(CurrentHDF%VGroupPath),    &
                                  trim(CurrentHDF%FieldName), OutputNumber = Instant,       &
                                  Imax = Imax, Jmax = Jmax, Kmax = Kmax, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR40'                                   
            endif
            
            if (CurrentHDF%From3Dto2D) then

                if ((Imax /= IUB - ILB + 1) .or.                                                &
                    (Jmax /= JUB - JLB + 1)) then
                    
                    write (*,*) trim(CurrentHDF%VGroupPath)
                    write (*,*) trim(CurrentHDF%FieldName)
                    write (*,*) 'miss match between the HDF5 input file and model domain'
                    stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR50'                                   
                    
                endif
                
                if (KUB /= 1) then
                    stop 'When integrating field 3D to 2D KUB number of layers must be 1'
                endif
                
                Field(:,:,:) = 0.
                                
                call HDF5SetLimits  (CurrentHDF%ObjHDF5, ILB, IUB, JLB, JUB, 1, Kmax, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR60'
                
                allocate(AuxField     (ILB-1:IUB+1, JLB-1:JUB+1, 0:Kmax+1))
                allocate(WaterPoints3D(ILB-1:IUB+1, JLB-1:JUB+1, 0:Kmax+1))  
                allocate(SZZ          (ILB-1:IUB+1, JLB-1:JUB+1, 0:Kmax+1))  
                     
                call HDF5ReadData(CurrentHDF%ObjHDF5, trim(CurrentHDF%VGroupPath),                      &
                                  trim(CurrentHDF%FieldName),                                       &
                                  Array3D = AuxField, OutputNumber = Instant, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR70'     
                
                call HDF5ReadData(CurrentHDF%ObjHDF5, "/Grid", "WaterPoints3D",    &
                                  Array3D = WaterPoints3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR80'                      
                
                call HDF5SetLimits  (CurrentHDF%ObjHDF5, ILB, IUB, JLB, JUB, 0, Kmax, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR90'
                

                call GetHDF5DataSetExist (HDF5ID = CurrentHDF%ObjHDF5, DataSetName = "/Grid/VerticalZ/Vertical_00001",&
                                          Exist  = Exist1, STAT = STAT_CALL)                                
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR100'
                

                call GetHDF5DataSetExist (HDF5ID = CurrentHDF%ObjHDF5, DataSetName = "/Grid/VerticalZ/VerticalZ_00001",&
                                          Exist  = Exist2, STAT = STAT_CALL)                                
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR110'
                
               
                if (Exist1) then

                    DataSetVert = "Vertical"            
                else
                
                    if (Exist2) then

                        DataSetVert = "VerticalZ"
                    
                    else
                        write(*,*) 'Missing the follow DataSet /Grid/VerticalZ/Vertical_00001'
                        write(*,*) '                       OR'                            
                        write(*,*) 'Missing the follow DataSet /Grid/VerticalZ/VerticalZ_00001'                        
                        stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR120'
                    endif

                endif
                
                    

                
                call HDF5ReadData(CurrentHDF%ObjHDF5, "/Grid/VerticalZ", trim(DataSetVert),&
                                  Array3D = SZZ, OutputNumber = Instant, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ModuleFillMatrix - ModuleFillMatrix - ERR130'                     

                do j = JLB, JUB
                do i = ILB, IUB

                    kbottom = -99
                    do k = 1, Kmax
                        if (WaterPoints3D (i,j,k)== 1) then
                            if (kbottom< 0) then
                                kbottom = k
                                HT      = (SZZ(i,j,kbottom-1) - SZZ(i,j,Kmax))
                            endif                            
                            if (HT > 0) then
                                Field(i,j,1) = Field(i,j,1) + AuxField(i,j,k) * (SZZ(i,j,k-1)-SZZ(i,j,k)) / HT
                            endif                            
                        endif                        
                    enddo
                
                enddo
                enddo
                     
              
                deallocate(AuxField, WaterPoints3D, SZZ)

            
            else
                                        
                if ((Imax /= IUB - ILB + 1) .or.                                                &
                    (Jmax /= JUB - JLB + 1) .or.                                                &
                    (Kmax /= KUB - KLB + 1)) then
                    
                    if (.not.(Kmax == 0 .and. KUB-KLB == 0)) then
                    
                        write (*,*) trim(CurrentHDF%VGroupPath)
                        write (*,*) trim(CurrentHDF%FieldName)
                        write (*,*) 'miss match between the HDF5 input file and model domain'
                        stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR140'                                   
                    
                    endif

                endif
                
                call HDF5SetLimits  (CurrentHDF%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR150'
                
                     
                call HDF5ReadData(CurrentHDF%ObjHDF5, trim(CurrentHDF%VGroupPath),                      &
                                  trim(CurrentHDF%FieldName),                                       &
                                  Array3D = CurrentHDF%ReadField3D, OutputNumber = Instant, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR160'                
              
            endif
            
            if (CurrentHDF%From2Dto3D) then    
           
                do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
                do j =               JLB,               JUB
                do i =               ILB,               IUB
                    Field(i,j,k) = CurrentHDF%ReadField3D(i,j,1)
                enddo
                enddo
                enddo
            
            else
               nullify(CurrentHDF%ReadField3D)  
            endif        

            if(CurrentHDF%HasMultiplyingFactor)then
            
                do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
                do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
                    Field(i,j,k) = Field(i,j,k) * CurrentHDF%MultiplyingFactor
                enddo
                enddo
                enddo
            
            end if

            if(CurrentHDF%HasAddingFactor)then
            
                do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
                do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
                    Field(i,j,k) = Field(i,j,k) + CurrentHDF%AddingFactor
                enddo
                enddo
                enddo

            end if            
            
        endif if4D            



    end subroutine ReadHDF5Values3D

    !-----------------------------------------------------------------------------------

    subroutine ModifyField4DInterpol(CurrentTime, Matrix3D, Matrix2D, CurrentHDF, Instant)
    
        !Arguments-------------------------------------------------------------
        type(T_Time)                                :: CurrentTime
        real, dimension(:,:,:), pointer, optional   :: Matrix3D
        real, dimension(:,:  ), pointer, optional   :: Matrix2D
        type(T_Field4D)                             :: CurrentHDF
        integer,                         optional   :: Instant        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), pointer             :: ZCellCenter 
        integer                                     :: i, j, k, icount
        integer                                     :: STAT_CALL
        type (T_Size2D)                             :: HaloMap
        integer                                     :: di, dj        
        logical                                     :: MasterOrSlave
        !Begin-----------------------------------------------------------------

        CurrentHDF%NoData   (:) = .true.
        
        call GetDDecompParameters(HorizontalGridID = Me%ObjHorizontalGrid,              &
                                  MasterOrSlave    = MasterOrSlave,                     &
                                  STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR10' 
                
ifMS:   if (MasterOrSlave) then
        
            call GetDDecompWorkSize2D(HorizontalGridID = Me%ObjHorizontalGrid,          &
                                      WorkSize         = HaloMap,                       &
                                      STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR20' 
                                                  
            di = HaloMap%ILB - 1
            dj = HaloMap%JLB - 1
            
        else ifMS
        
            di = 0
            dj = 0
        
        endif ifMS                
        

if2D:   if (Me%Dim == Dim2D) then

            call ModifyField4DXYZ(Field4DID             = CurrentHDF%ObjField4D,            &
                                  PropertyIDNumber      = Me%PropertyID%IDNumber,           &
                                  CurrentTime           = CurrentTime,                      &
                                  X                     = CurrentHDF%X,                     &
                                  Y                     = CurrentHDF%Y,                     &
                                  Field                 = CurrentHDF%Prop,                  &
                                  NoData                = CurrentHDF%NoData,                &
                                  Instant               = Instant,                          &
                                  STAT                  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR30' 
            
            icount = 0
            
            do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB        
                
                if (Me%PointsToFill2D(i,j) == WaterPoint) then
                    
                    icount           = icount + 1
                    if (CurrentHDF%NoData(icount)) then
                        if (CurrentHDF%Extrapolate) then
                            Matrix2D(i, j)   = Me%DefaultValue(1)
                        else
                            write(*,*) GetPropertyName (Me%PropertyID%IDNumber)
                            write(*,*) 'No data in 2D cell I=',i + di, 'J=',j + dj
                            stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR40' 
                        endif                        
                    else                        
                        Matrix2D(i, j)   = CurrentHDF%Prop(icount)
                    endif
                endif                    
                
            enddo   
            enddo      
            
        
        else if2D
                       
F2D3D:      if (CurrentHDF%From2Dto3D) then            
            
                call ModifyField4DXYZ(Field4DID             = CurrentHDF%ObjField4D,            &
                                      PropertyIDNumber      = Me%PropertyID%IDNumber,           &
                                      CurrentTime           = CurrentTime,                      &
                                      X                     = CurrentHDF%X,                     &
                                      Y                     = CurrentHDF%Y,                     &
                                      Field                 = CurrentHDF%Prop,                  &
                                      NoData                = CurrentHDF%NoData,                &
                                      Instant               = Instant,                          &                                  
                                      STAT                  = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR50' 
            
                icount = 0
            
                do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB        
                    if (Me%PointsToFill3D(i,j,Me%WorkSize3D%KUB) == WaterPoint) then                    
                        icount           = icount + 1
                        if (CurrentHDF%NoData(icount)) then
                            if (CurrentHDF%Extrapolate) then 
                                !NeedToExtrapolate = .true.
                                Matrix3D(i, j, k) = Me%DefaultValue(1)
                            else                             
                                write(*,*) 'No data in 2D cell I=',i + di, 'J=',j + dj
                                stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR60' 
                            endif                                
                        else                        
                            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
                                if (Me%PointsToFill3D(i,j,k) == WaterPoint) then                                
                                    Matrix3D(i, j, k)   = CurrentHDF%Prop(icount)
                                endif
                            enddo                                        
                        endif
                    endif                    
                enddo
                enddo   
    
            else F2D3D
        
                call GetGeometryDistances(Me%ObjGeometry, ZCellCenter, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR70' 
        
                icount = 0
            
                do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
                do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB        
                
                    if (Me%PointsToFill3D(i,j,k) == WaterPoint) then
                    
                        icount           = icount + 1
                        CurrentHDF%Z(icount) = ZCellCenter(i, j, k)
                
                    endif                    
                
                enddo
                enddo   
                enddo            
            
                call UnGetGeometry(Me%ObjGeometry, ZCellCenter, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR80'             
        
                call ModifyField4DXYZ(Field4DID             = CurrentHDF%ObjField4D,            &
                                      PropertyIDNumber      = Me%PropertyID%IDNumber,           &
                                      CurrentTime           = CurrentTime,                      &
                                      X                     = CurrentHDF%X,                     &
                                      Y                     = CurrentHDF%Y,                     &
                                      Z                     = CurrentHDF%Z,                     &
                                      Field                 = CurrentHDF%Prop,                  &
                                      NoData                = CurrentHDF%NoData,                &
                                      Instant               = Instant,                          &                                  
                                      STAT                  = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR90' 
            
                icount = 0
            
                do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
                do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB        
                
                    if (Me%PointsToFill3D(i,j,k) == WaterPoint) then
                    
                        icount           = icount + 1
                        if (CurrentHDF%NoData(icount)) then
                            if (CurrentHDF%Extrapolate) then 
                                !NeedToExtrapolate = .true.
                                Matrix3D(i, j, k) = Me%DefaultValue(1)
                            else    
                                write(*,*) 'No data in 3D cell I=',i + di, 'J=',j + dj, 'K=',k
                                stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR100'
                            endif                               
                        else                        
                            Matrix3D(i, j, k)   = CurrentHDF%Prop(icount)
                        endif
                
                    endif                    
                
                enddo
                enddo   
                enddo   
                
            endif F2D3D        
        endif if2D
        
        if (icount > CurrentHDF%Ncells) then
            write(*,*) 'icount =', icount
            write(*,*) 'CurrentHDF%Ncells =', CurrentHDF%Ncells
            stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR110' 
        endif         
        

    end subroutine ModifyField4DInterpol
    
    !-----------------------------------------------------------------------------------    
    

    subroutine ProfileTimeSerieField(PointsToFill3D, Instant, Field)

        !Arguments-------------------------------------------------------------
        integer, dimension(:, :, :), pointer    :: PointsToFill3D
        integer                                 :: Instant
        real, dimension(:,:,:), pointer         :: Field

        !Local-----------------------------------------------------------------
        real                                    :: CellDepth
        integer                                 :: STAT_CALL, i, j, k, NDEPTHS
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        real, dimension(:,:,:), pointer         :: SZZ
        real, dimension(:    ),     pointer     :: Values, Depth

        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize3D%ILB
        IUB = Me%WorkSize3D%IUB
        JLB = Me%WorkSize3D%JLB
        JUB = Me%WorkSize3D%JUB
        KLB = Me%WorkSize3D%KLB
        KUB = Me%WorkSize3D%KUB

        nullify(Values, Depth)

        !Gets Center of the cells
        call GetGeometryDistances(Me%ObjGeometry, SZZ = SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ProfileTimeSerieField - ModuleFillMatrix - ERR01'


        Values => Me%ProfileTimeSerie%Values(:, Instant)
        Depth  => Me%ProfileTimeSerie%Depths(:, Instant)

        NDEPTHS = Me%ProfileTimeSerie%nDepths

        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB

            if (PointsToFill3D(i, j, k) == WaterPoint) then

                CellDepth          = (SZZ(i, j, k) + SZZ(i, j, k - 1)) / 2 - SZZ(i, j, KUB) 
                Field(i,j,k)       = InterpolateProfile (CellDepth, NDEPTHS, Depth, Values)

            endif

        enddo
        enddo
        enddo

        call UnGetGeometry (Me%ObjGeometry, SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ProfileTimeSerieField - ModuleFillMatrix - ERR14'

    end subroutine ProfileTimeSerieField

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetIfMatrixRemainsConstant (FillMatrixID, RemainsConstant, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        logical, intent(OUT)                            :: RemainsConstant
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            RemainsConstant = Me%RemainsConstant

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetIfMatrixRemainsConstant

    !--------------------------------------------------------------------------

    subroutine GetDefaultValueScalar (FillMatrixID, DefaultValue, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        real,    intent(OUT)                            :: DefaultValue
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            DefaultValue = Me%DefaultValue(1)

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetDefaultValueScalar

    !--------------------------------------------------------------------------

    subroutine GetIsVectorialProp (FillMatrixID, IsVectorialProp, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        logical,  intent(OUT)                           :: IsVectorialProp
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            IsVectorialProp = Me%VectorialProp

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetIsVectorialProp

    !--------------------------------------------------------------------------       
    
    subroutine GetDefaultValueVectorial (FillMatrixID, DefaultValue, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        real, dimension(3),  intent(OUT)                :: DefaultValue
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            DefaultValue = Me%DefaultValue

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetDefaultValueVectorial

    !--------------------------------------------------------------------------    
    !Get original field (used for output)
    subroutine GetVectorialField2D (FillMatrixID, FieldX, FieldY, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        real, dimension(:, :),  intent(OUT), pointer    :: FieldX
        real, dimension(:, :),  intent(OUT), pointer    :: FieldY
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_ )) then
            
            call Read_Lock(mFILLMATRIX_, Me%InstanceID)
            call Read_Lock(mFILLMATRIX_, Me%InstanceID)
            
            if (Me%VectorialProp) then
                FieldX = Me%Matrix2DX
                FieldY = Me%Matrix2DY
            else
                FieldX => null()
                FieldY => null()           
            endif

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetVectorialField2D

    !--------------------------------------------------------------------------        
    
    !Get original field (used for output)
    subroutine GetVectorialField3D (FillMatrixID, FieldX, FieldY, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        real, dimension(:, :, :),  intent(OUT), pointer :: FieldX
        real, dimension(:, :, :),  intent(OUT), pointer :: FieldY
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_ )) then
            
            call Read_Lock(mFILLMATRIX_, Me%InstanceID)
            call Read_Lock(mFILLMATRIX_, Me%InstanceID)
            
            if (Me%VectorialProp) then
                FieldX = Me%Matrix3DX
                FieldY = Me%Matrix3DY
            else
                FieldX => null()
                FieldY => null()           
            endif

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetVectorialField3D

    !--------------------------------------------------------------------------        
    
    
    !Get original field (used for output)
    subroutine GetAngleField2D (FillMatrixID, Field, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        real, dimension(:, :),  intent(OUT), pointer    :: Field
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_ )) then
            
            call Read_Lock(mFILLMATRIX_, Me%InstanceID)
            
!~             if (Check_Angle_Property(Me%PropertyID%IDNumber)) then
            if (Me%PropertyID%IsAngle) then
                Field = Me%Matrix2DFieldAngle
            else
                Field => null()          
            endif

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetAngleField2D

    !--------------------------------------------------------------------------        
    
    !Get original field (used for output)
    subroutine GetAngleField3D (FillMatrixID, Field, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        real, dimension(:, :, :),  intent(OUT), pointer :: Field
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_ )) then
            
            call Read_Lock(mFILLMATRIX_, Me%InstanceID)
            
!~             if (Check_Angle_Property(Me%PropertyID%IDNumber)) then
            if (Me%PropertyID%IsAngle) then
                Field = Me%Matrix3DFieldAngle
            else
                Field => null()          
            endif

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetAngleField3D

    !--------------------------------------------------------------------------
    
    subroutine GetFillMatrixDTPrediction (FillMatrixID, PredictedDT, DTForNextEvent, DTForNextDataset, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        real, intent(OUT), optional                     :: PredictedDT,     &
                                                           DTForNextEvent,  &
                                                           DTForNextDataset
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PredictedDT)) &
                PredictedDT = Me%PredictedDT
                
            if (present(DTForNextEvent)) &
                DTForNextEvent = Me%DTForNextEvent

            if (present(DTForNextDataset)) &
                DTForNextDataset = Me%DTForNextDataset

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetFillMatrixDTPrediction


    !--------------------------------------------------------------------------

    subroutine GetHDFTimeLimits (FillMatrixID, StartTime, EndTime, HaveTimeLimits, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        type(T_Time)                                    :: StartTime, EndTime
        logical                                         :: HaveTimeLimits
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            if (Me%TimeEvolution == ReadHDF) then
                HaveTimeLimits = .true.
                
                StartTime = Me%FirstHDF%StartTime
                EndTime   = Me%FirstHDF%EndTime
                
            else
                HaveTimeLimits = .false.
                
                call null_time(StartTime)
                call null_time(EndTime  )
                
            endif


            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetHDFTimeLimits
    
    !--------------------------------------------------------------------------
    
    subroutine GetHDFTimeLimitsVectorial (FillMatrixID, StartTimeX, EndTimeX, StartTimeY, EndTimeY, HaveTimeLimits, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        type(T_Time)                                    :: StartTimeX, EndTimeX, StartTimeY, EndTimeY
        logical                                         :: HaveTimeLimits
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            if (Me%TimeEvolution == ReadHDF) then
                HaveTimeLimits = .true.
                
                StartTimeX = Me%FirstHDF%StartTime
                EndTimeX   = Me%FirstHDF%EndTime
                
                StartTimeY = Me%FirstHDF%Next%StartTime
                EndTimeY   = Me%FirstHDF%Next%EndTime                
                
            else
                HaveTimeLimits = .false.
                
                call null_time(StartTimeX)
                call null_time(StartTimeY)
                call null_time(EndTimeX  )
                call null_time(EndTimeY  )
                
            endif


            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetHDFTimeLimitsVectorial
    
    !--------------------------------------------------------------------------    
    
    subroutine GetNumberOfInstants (FillMatrixID, NumberOfInstants, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        integer, intent(OUT)                            :: NumberOfInstants
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            NumberOfInstants = Me%FirstHDF%NumberOfInstants

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
        
        !----------------------------------------------------------------------
    
    end subroutine GetNumberOfInstants
    
    !--------------------------------------------------------------------------
    
    subroutine GetNumberOfInstantsVectorial (FillMatrixID, NumberOfInstantsX, NumberOfInstantsY, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        integer, intent(OUT)                            :: NumberOfInstantsX, NumberOfInstantsY
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            if (Me%VectorialProp) then
                NumberOfInstantsX = Me%FirstHDF%NumberOfInstants
                NumberOfInstantsY = Me%FirstHDF%Next%NumberOfInstants
            else
                NumberOfInstantsX = null_int
                NumberOfInstantsY = null_int                
            endif
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
        
        !----------------------------------------------------------------------
    
    end subroutine GetNumberOfInstantsVectorial
    
    !--------------------------------------------------------------------------    
    
    subroutine GetTimeInstant(FillMatrixID, Instant, TimeInstant, STAT)

        !Arguments-------------------------------------------------------------
        integer                                 :: FillMatrixID
        integer                                 :: Instant
        type(T_Time), intent(OUT)               :: TimeInstant 
        integer, intent(OUT), optional          :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, ready_
        real,    dimension(:), pointer          :: TimeVector

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call HDF5SetLimits  (Me%FirstHDF%ObjHDF5, 1, 6, STAT = STAT_)

            allocate(TimeVector(6))

            call HDF5ReadData   (HDF5ID         = Me%FirstHDF%ObjHDF5,   &
                                 GroupName      = "/Time",          &
                                 Name           = "Time",           &
                                 Array1D        = TimeVector,       &
                                 OutputNumber   = Instant,          &
                                 STAT           = STAT_)
            if (STAT_ /= SUCCESS_) stop 'GetTimeInstant - ModuleFillMatrix - ERR010'

            call SetDate(TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))
                                     
            deallocate(TimeVector)
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
        
        !----------------------------------------------------------------------

    end subroutine GetTimeInstant

    !----------------------------------------------------------------------
    
    subroutine GetTimeInstantVectorial(FillMatrixID, Instant, TimeInstantX, TimeInstantY, STAT)

        !Arguments-------------------------------------------------------------
        integer                                 :: FillMatrixID
        integer                                 :: Instant
        type(T_Time), intent(OUT)               :: TimeInstantX, TimeInstantY
        integer, intent(OUT), optional          :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, ready_
        real,    dimension(:), pointer          :: TimeVector

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (Me%VectorialProp) then
                
                call HDF5SetLimits  (Me%FirstHDF%ObjHDF5, 1, 6, STAT = STAT_)

                allocate(TimeVector(6))

                call HDF5ReadData   (HDF5ID         = Me%FirstHDF%ObjHDF5,   &
                                     GroupName      = "/Time",          &
                                     Name           = "Time",           &
                                     Array1D        = TimeVector,       &
                                     OutputNumber   = Instant,          &
                                     STAT           = STAT_)
                if (STAT_ /= SUCCESS_) stop 'GetTimeInstant - ModuleFillMatrix - ERR010'

                call SetDate(TimeInstantX, Year     = TimeVector(1), Month  = TimeVector(2), &
                                          Day      = TimeVector(3), Hour   = TimeVector(4), &
                                          Minute   = TimeVector(5), Second = TimeVector(6))
            
            
            
                call HDF5SetLimits  (Me%FirstHDF%Next%ObjHDF5, 1, 6, STAT = STAT_)

                call HDF5ReadData   (HDF5ID         = Me%FirstHDF%Next%ObjHDF5,   &
                                     GroupName      = "/Time",          &
                                     Name           = "Time",           &
                                     Array1D        = TimeVector,       &
                                     OutputNumber   = Instant,          &
                                     STAT           = STAT_)
                if (STAT_ /= SUCCESS_) stop 'GetTimeInstant - ModuleFillMatrix - ERR010'

                call SetDate(TimeInstantY, Year     = TimeVector(1), Month  = TimeVector(2), &
                                          Day      = TimeVector(3), Hour   = TimeVector(4), &
                                          Minute   = TimeVector(5), Second = TimeVector(6))            
                                     
                deallocate(TimeVector)
                
            else
                call  null_time(TimeInstantX)
                call  null_time(TimeInstantY)
            endif
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
        
        !----------------------------------------------------------------------

    end subroutine GetTimeInstantVectorial

    !----------------------------------------------------------------------
    
    subroutine GetNextValueForDTPred(FillMatrixID, max_value, STAT)

        !Arguments-------------------------------------------------------------
        integer                                 :: FillMatrixID
        real, intent(OUT)                       :: max_value
        integer, intent(OUT), optional          :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            max_value = Me%NextValueForDTPred
            
            STAT_ = SUCCESS_

        else 

            STAT_ = ready_

        end if

        if (present(STAT)) STAT = STAT_
        
        !----------------------------------------------------------------------

    end subroutine GetNextValueForDTPred
    
    !--------------------------------------------------------------------------
    
    subroutine GetValuesProcessingOptions(FillMatrixID, Accumulate, Interpolate, UseOriginal, STAT)

        !Arguments-------------------------------------------------------------
        integer                                 :: FillMatrixID
        logical, intent(OUT), optional          :: Accumulate, Interpolate, UseOriginal
        integer, intent(OUT), optional          :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(Accumulate)) &
                Accumulate = Me%AccumulateValues
                
            if (present(Interpolate)) &
                Interpolate = Me%InterpolateValues
            
            if (present(UseOriginal)) &
                UseOriginal = Me%UseOriginalValues
            
            STAT_ = SUCCESS_

        else 

            STAT_ = ready_

        end if

        if (present(STAT)) STAT = STAT_
        
        !----------------------------------------------------------------------

    end subroutine GetValuesProcessingOptions   
    
    !--------------------------------------------------------------------------
    
    subroutine GetMultiTimeSeries (FillMatrixID, List, Counter, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                                             :: FillMatrixID
        type(T_Station), dimension(:), pointer, intent(OUT) :: List
        integer, intent(OUT)                                :: Counter
        integer, intent(OUT), optional                      :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                             :: STAT_, ready_

        !----------------------------------------------------------------------

        !This routine does not "lock", but just because the property do not
        !change after it's construction. if this behaviour changes, it must 
        !implement a "lock" system, like in others modules
        
        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            List => Me%MultiTimeSerie%StationsList
            Counter = Me%MultiTimeSerie%NumberOfSources
            
            STAT_ = SUCCESS_

        else 

            STAT_ = ready_

        end if

        if (present(STAT)) STAT = STAT_
        
        !----------------------------------------------------------------------
        
    end subroutine GetMultiTimeSeries
    
    !--------------------------------------------------------------------------

    !Get if Analytic Celerity is ON
    subroutine GetAnalyticCelerityON (FillMatrixID, AnalyticCelerityON, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        logical                                         :: AnalyticCelerityON
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_ )) then
            
            AnalyticCelerityON = Me%AnalyticWave%ON

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetAnalyticCelerityON

    !--------------------------------------------------------------------------  
        
    !Get Analytic Celerity
    subroutine GetAnalyticCelerity (FillMatrixID, AnalyticCelerity, AnalyticDirection, AnalyticAverageValue, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        real,    dimension(:, :),  intent(OUT), pointer :: AnalyticCelerity
        real                                            :: AnalyticDirection
        real                                            :: AnalyticAverageValue
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_ )) then
            
            call Read_Lock(mFILLMATRIX_, Me%InstanceID)

            AnalyticCelerity => Me%AnalyticWave%Celerity
            
            AnalyticDirection    = Me%AnalyticWave%Direction
            
            AnalyticAverageValue = Me%AnalyticWave%AverageValue

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetAnalyticCelerity

    !--------------------------------------------------------------------------        
    
    


    !--------------------------------------------------------------------------  
        
    !Get hdf filename
    subroutine GetFilenameHDF (FillMatrixID, FilenameHDF, HdfFileExist, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        character(len=*)                                :: FilenameHDF
        logical                                         :: HdfFileExist
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_ )) then
            
            if (associated(Me%FirstHDF)) then
                FilenameHDF     = Me%FirstHDF%FileName
                HdfFileExist    = .true.
            else
                FilenameHDF     = null_str
                HdfFileExist    = .false.                
            endif

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetFilenameHDF

    !--------------------------------------------------------------------------    
   
    
    subroutine UngetFillMatrix2D(FillMatrixID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: FillMatrixID
        real, pointer, dimension(:,:)   :: Array
        integer, optional, intent (OUT) :: STAT
   
        !External--------------------------------------------------------------
        integer                         :: ready_   

        !Local-----------------------------------------------------------------
        integer                         :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mFILLMATRIX_, Me%InstanceID, "UngetFillMatrix2D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetFillMatrix2D

    !--------------------------------------------------------------------------    

    subroutine UngetFillMatrix3D(FillMatrixID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: FillMatrixID
        real, pointer, dimension(:,:,:) :: Array
        integer, optional, intent (OUT) :: STAT
   
        !External--------------------------------------------------------------
        integer                         :: ready_   

        !Local-----------------------------------------------------------------
        integer                         :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mFILLMATRIX_, Me%InstanceID, "UngetFillMatrix3D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetFillMatrix3D

    !--------------------------------------------------------------------------      
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyFillMatrix(FillMatrixID, Matrix2D, Matrix3D, Matrix2DInputRef, Matrix3DInputRef, PointsToFill2D,       &
                                PointsToFill3D, Generic_4D_Value, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        real,    dimension(:, :),    pointer, optional  :: Matrix2D
        real,    dimension(:, :, :), pointer, optional  :: Matrix3D
        real,    dimension(:, :),    pointer, optional  :: Matrix2DInputRef        !original field (e.g. angle)
        real,    dimension(:, :, :), pointer, optional  :: Matrix3DInputRef        
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D
        real,                    intent( IN), optional  :: Generic_4D_Value
        integer,                 intent(OUT), optional  :: STAT

        !Local-----------------------------------------------------------------
        real                                            :: Generic_4D_Value_        
        integer                                         :: STAT_, STAT_CALL, ready_, Referential
        logical                                         :: ModifyError = .false.
        type(T_Field4D), pointer                        :: CurrentHDF
        type(T_TimeSerie), pointer                      :: CurrentTimeSerie
        type (T_PropertyID), pointer                    :: Prop
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleFillMatrix", "ModifyFillMatrix")
        
        
        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            if (present(Matrix2D)) Me%Matrix2D => Matrix2D
            if (present(Matrix3D)) Me%Matrix3D => Matrix3D
            
            !!get the orginal field. will be given to user to output
            if (Me%RotateAngleToGrid) then
                if (.not. present(Matrix2DInputRef) .and. .not. present(Matrix3DInputRef)) then
                    write(*,*) 'Constructing angle property but not given original field'
                    stop 'ModifyFillMatrix - ModuleFillMatrix - ERR10'                    
                endif
                
                if (present(Matrix2D))         Me%Matrix2DCellAngle => Matrix2D
                if (present(Matrix3D))         Me%Matrix3DCellAngle => Matrix3D                
                if (present(Matrix2DInputRef)) Me%Matrix2DFieldAngle => Matrix2DInputRef
                if (present(Matrix3DInputRef)) Me%Matrix3DFieldAngle => Matrix3DInputRef                  
            endif             
          
            
            if (present(PointsToFill2D)) Me%PointsToFill2D => PointsToFill2D
            if (present(PointsToFill3D)) Me%PointsToFill3D => PointsToFill3D



            if (present(Generic_4D_Value)) then
                Generic_4D_Value_ = Generic_4D_Value
            else
                if (associated(Me%FirstHDF)) then
                    if (Me%FirstHDF%Generic4D%ON) then
                        ModifyError = .true. 
                        write(*,*) 'The FillMatrix wants to interpolate along a Generic 4D dimension'
                        write(*,*) 'However, no data is provide for the interpolation'
                    endif
                endif
                Generic_4D_Value_ = FillValueReal
            endif

            !reset
            if (Me%PredictDTMethod == 2) then
                Me%DTForNextEvent       = -null_real
                Me%PredictedDT          = -null_real
                Me%DTForNextDataset     = -null_real
                Me%NextValueForDTPred   = 0.0            
            endif            

            if (.not. ModifyError) then            
                                                
                select case (Me%TimeEvolution)

                    case (ReadTimeSerie)
                        !Not vectorial prop - only one time serie
                        CurrentTimeSerie => Me%FirstTimeSerie

                        !Associate Matrix2D and Me%Matrix3D to the input field ones
                        if (Me%RotateAngleToGrid) call AssociateMatrixes(1)                 
                        
                        if (Me%Dim == Dim2D) then
                            call ModifySpaceTimeSerie    (CurrentTimeSerie, PointsToFill2D = PointsToFill2D) 
                        else
                            call ModifySpaceTimeSerie    (CurrentTimeSerie, PointsToFill3D = PointsToFill3D) 
                        endif
                
                    case (ReadHDF)


                        !Not vectorial prop - only one hdf
                        CurrentHDF => Me%FirstHDF    
                        
                        if(.not. CurrentHDF%RemainsConstant)then
                            
                            !Associate Matrix2D and Me%Matrix3D to the input field ones
                            if (Me%RotateAngleToGrid) call AssociateMatrixes(1)                              
                            
                            if (Me%Dim == Dim2D) then
                                call ModifyHDFInput2D (PointsToFill2D, CurrentHDF, Generic_4D_Value_ ) 
                            else
                                call ModifyHDFInput3D (PointsToFill3D, CurrentHDF, Generic_4D_Value_ )
                            endif

                        end if

                                    
                    case (ProfileTimeSerie)

                        if(.not. Me%RemainsConstant)then

                            call ModifyProfileTimeSerie (PointsToFill3D = PointsToFill3D) 

                        end if
                        
                    case (MultiTimeSerie)

                        if (Me%Dim == Dim2D) then                            
                            call ModifyMultiTimeSerie (PointsToFill2D = PointsToFill2D) 
                        else                            
                            call ModifyMultiTimeSerie (PointsToFill3D = PointsToFill3D) 
                        endif  
                        
                    case(AnalyticWave)

                        if (Me%Dim == Dim2D) then
                            call ModifyAnalyticWave    (PointsToFill2D = PointsToFill2D) 
                        else
                            call ModifyAnalyticWave    (PointsToFill3D = PointsToFill3D) 
                        endif                             

                end select                
                
                !Is this a angle property? convert to cell referential angle
                if (Me%RotateAngleToGrid) then
                    
                    !angle referential
                    Prop => Me%PropertyID
                    Referential = Get_Angle_Referential(Prop)
                    
                    if (Me%Dim == Dim2D) then                          
                        
                        !!Need to rotate input field          
                        call RotateAngleFieldToGrid(HorizontalGridID      = Me%ObjHorizontalGrid,               &
                                                        AngleIn           = Me%Matrix2DFieldAngle,                &
                                                        InReferential     = Referential,                          &
                                                        AngleOut          = Me%Matrix2DCellAngle,                 &
                                                        WaterPoints2D     = PointsToFill2D,                       &
                                                        Rotate            = .true.,                               &
                                                        STAT              = STAT_CALL)                    
                    else
                                             
                        
                        !!Need to rotate input field            
                        call RotateAngleFieldToGrid(HorizontalGridID      = Me%ObjHorizontalGrid,               &
                                                        AngleIn           = Me%Matrix3DFieldAngle,                &
                                                        InReferential     = Referential,                          &
                                                        AngleOut          = Me%Matrix3DCellAngle,                 &
                                                        WaterPoints3D     = PointsToFill3D,                       &
                                                        Rotate            = .true.,                               &
                                                        KLB               = Me%WorkSize3D%KLB,                    &
                                                        KUB               = Me%WorkSize3D%KUB,                    &
                                                        STAT              = STAT_CALL)                  
                    endif                    
                     
                
                endif        
                
                nullify(Me%Matrix2D)
                nullify(Me%Matrix3D)
                nullify(Me%Matrix2DFieldAngle)
                nullify(Me%Matrix3DFieldAngle)        
                nullify(Me%Matrix2DCellAngle)
                nullify(Me%Matrix3DCellAngle)                   

                nullify(Me%PointsToFill2D)
                nullify(Me%PointsToFill3D)

                STAT_ = SUCCESS_
            else

                STAT_ = UNKNOWN_

            endif
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        if (MonitorPerformance) call StopWatch ("ModuleFillMatrix", "ModifyFillMatrix")
        
        
    end subroutine ModifyFillMatrix

    !--------------------------------------------------------------------------

    subroutine ModifyFillMatrixVectorial (FillMatrixID,     &
                                          Matrix2DU,        &
                                          Matrix2DV,        &
                                          Matrix3DU,        &
                                          Matrix3DV,        & 
                                          Matrix3DW,        &
                                          Matrix2DX,        &
                                          Matrix2DY,        &
                                          Matrix3DX,        &
                                          Matrix3DY,        &
                                          PointsToFill2D,   &
                                          PointsToFill3D,   &
                                          Generic_4D_Value, &
                                          STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        real,    dimension(:, :),    pointer, optional  :: Matrix2DU
        real,    dimension(:, :),    pointer, optional  :: Matrix2DV
        real,    dimension(:, :, :), pointer, optional  :: Matrix3DU
        real,    dimension(:, :, :), pointer, optional  :: Matrix3DV
        real,    dimension(:, :, :), pointer, optional  :: Matrix3DW
        real,    dimension(:, :),    pointer, optional  :: Matrix2DX
        real,    dimension(:, :),    pointer, optional  :: Matrix2DY
        real,    dimension(:, :, :), pointer, optional  :: Matrix3DX
        real,    dimension(:, :, :), pointer, optional  :: Matrix3DY        
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D
        real,                    intent( IN), optional  :: Generic_4D_Value
        integer,                 intent(OUT), optional  :: STAT

        !Local-----------------------------------------------------------------
        real                                            :: Generic_4D_Value_        
        integer                                         :: STAT_, STAT_CALL, ready_, file
        logical                                         :: ModifyError = .false.
        type(T_Field4D), pointer                        :: CurrentHDF
        type(T_TimeSerie), pointer                      :: CurrentTimeSerie        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            
            if (Me%Dim == Dim3D .and. Me%UseZ .and. .not. present(Matrix3DW)) then
                write(*,*) 'Modifying vectorial property that needs W component to be given'
                stop 'ModifyFillMatrixVectorial - ModuleFillMatrix - ERR01'            
            endif             
            
            if (present(Matrix2DU)) Me%Matrix2DU => Matrix2DU
            if (present(Matrix2DV)) Me%Matrix2DV => Matrix2DV
            if (present(Matrix3DU)) Me%Matrix3DU => Matrix3DU
            if (present(Matrix3DV)) Me%Matrix3DV => Matrix3DV
            if (present(Matrix3DV)) Me%Matrix3DW => Matrix3DW
            
            if (present(Matrix2DX)) Me%Matrix2DX => Matrix2DX
            if (present(Matrix2DY)) Me%Matrix2DY => Matrix2DY
            if (present(Matrix3DX)) Me%Matrix3DX => Matrix3DX
            if (present(Matrix3DY)) Me%Matrix3DY => Matrix3DY            
            
            if (present(PointsToFill2D)) Me%PointsToFill2D => PointsToFill2D
            if (present(PointsToFill3D)) Me%PointsToFill3D => PointsToFill3D


            if (present(Generic_4D_Value)) then
                Generic_4D_Value_ = Generic_4D_Value
            else
                if (associated(Me%FirstHDF)) then
                    !If vectorial both x and y have same generic 4d property, test only one
                    if (Me%FirstHDF%Generic4D%ON) then
                        ModifyError = .true. 
                        write(*,*) 'The FillMatrix wants to interpolate along a Generic 4D dimension'
                        write(*,*) 'However, no data is provide for the interpolation'
                    endif
                endif
                Generic_4D_Value_ = FillValueReal
            endif
            
            !reset
            if (Me%PredictDTMethod == 2) then
                Me%DTForNextEvent       = -null_real
                Me%PredictedDT          = -null_real
                Me%DTForNextDataset     = -null_real
                Me%NextValueForDTPred   = 0.0            
            endif
        
            if (.not. ModifyError) then
                                
                
                select case (Me%TimeEvolution)

                    case (ReadTimeSerie)

                        !Vectorial prop - two time series
                        CurrentTimeSerie => Me%FirstTimeSerie
                        file = 1
                        do while (associated(CurrentTimeSerie))
                                
                            !Associate Matrix2D and Me%Matrix3D to the input field ones
                            call AssociateMatrixes(file)                                        
                            
                            if (Me%Dim == Dim2D) then
                                call ModifySpaceTimeSerie    (CurrentTimeSerie, PointsToFill2D = PointsToFill2D) 
                            else
                                call ModifySpaceTimeSerie    (CurrentTimeSerie, PointsToFill3D = PointsToFill3D) 
                            endif                                
                                
                            file = file + 1
                            CurrentTimeSerie => CurrentTimeSerie%Next
                        enddo
                
                    case (ReadHDF)

                        !Vectorial prop - two time series
                        CurrentHDF => Me%FirstHDF
                        file = 1
                                                    
                        do while (associated(CurrentHDF))              
                            
                            if(.not. CurrentHDF%RemainsConstant)then                                         
                                
                                !Associate Matrix2D and Me%Matrix3D to the input field ones
                                call AssociateMatrixes(file)                                            
                                
                                if (Me%Dim == Dim2D) then
                                    call ModifyHDFInput2D (PointsToFill2D, CurrentHDF, Generic_4D_Value_ ) 
                                else
                                    call ModifyHDFInput3D (PointsToFill3D, CurrentHDF, Generic_4D_Value_)
                                endif
                                
                                file = file + 1
                                CurrentHDF => CurrentHDF%Next
                                                                
                            end if
                        enddo
                                    
                    case (ProfileTimeSerie)

                        if(.not. Me%RemainsConstant)then

                            call ModifyProfileTimeSerie (PointsToFill3D = PointsToFill3D) 

                        end if
                        
                    case (MultiTimeSerie)

                        if (Me%Dim == Dim2D) then                            
                            call ModifyMultiTimeSerie (PointsToFill2D = PointsToFill2D) 
                        else                            
                            call ModifyMultiTimeSerie (PointsToFill3D = PointsToFill3D) 
                        endif  
    
                end select                                                        
                
                if (Me%Dim == Dim2D) then  
                    !!Need to rotate input field (Me%Matrix2DX and Me%Matrix2DY) to grid (Me%Matrix2DU and Me%Matrix2DV))            
                    call RotateVectorFieldToGrid(HorizontalGridID  = Me%ObjHorizontalGrid, &
                                                    VectorInX         = Me%Matrix2DX,                         &
                                                    VectorInY         = Me%Matrix2DY,                         &
                                                    VectorOutX        = Me%Matrix2DU,                         &
                                                    VectorOutY        = Me%Matrix2DV,                         &   
                                                    WaterPoints2D     = PointsToFill2D,                       &
                                                    RotateX           = .true.,                               &
                                                    RotateY           = .true.,                               &
                                                    STAT              = STAT_CALL)                    
                else
                    !!Need to rotate input field (Me%Matrix3DX and Me%Matrix3DY) to grid (Me%Matrix3DU and Me%Matrix3DV))            
                    call RotateVectorFieldToGrid(HorizontalGridID  = Me%ObjHorizontalGrid, &
                                                    VectorInX         = Me%Matrix3DX,                         &
                                                    VectorInY         = Me%Matrix3DY,                         &
                                                    VectorOutX        = Me%Matrix3DU,                         &
                                                    VectorOutY        = Me%Matrix3DV,                         &   
                                                    WaterPoints3D     = PointsToFill3D,                       &
                                                    RotateX           = .true.,                               &
                                                    RotateY           = .true.,                               &
                                                    KLB               = Me%WorkSize3D%KLB,                    &
                                                    KUB               = Me%WorkSize3D%KUB,                    &                    
                                                    STAT              = STAT_CALL)                    
                endif
                
                nullify(Me%Matrix2D)
                nullify(Me%Matrix3D)
                nullify(Me%Matrix2DU)
                nullify(Me%Matrix2DV)
                nullify(Me%Matrix3DU)
                nullify(Me%Matrix3DV)   
                nullify(Me%Matrix2DX)
                nullify(Me%Matrix2DY)
                nullify(Me%Matrix3DX)
                nullify(Me%Matrix3DY)                     
                nullify(Me%PointsToFill2D)
                nullify(Me%PointsToFill3D)
                
                
                STAT_ = SUCCESS_
            else

                STAT_ = UNKNOWN_

            endif
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyFillMatrixVectorial

    !--------------------------------------------------------------------------                                
    
    subroutine ModifyMultiTimeSerie (PointsToFill2D, PointsToFill3D)
    
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer, optional    :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional :: PointsToFill3D             

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: ilb, iub, jlb, jub, klb, kub
        integer                                         :: i, j, k  
        integer, dimension(:,:), pointer                :: id2d => null()
        integer, dimension(:,:,:), pointer              :: id3d => null()
        type(T_Station), dimension(:), pointer          :: sl   => null()
        type (T_Time)                                   :: Now
        type (T_Time)                                   :: Time1, Time2
        real                                            :: Value1, Value2
        logical                                         :: TimeCycle
        logical                                         :: compare_dt = .false.
        real                                            :: PredictedDT
        real                                            :: DTForNextEvent
        integer                                         :: index
        
        !----------------------------------------------------------------------
        
        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyMultiTimeSerie - ModuleFillMatrix - ERR010'

        if (Me%PredictDTMethod == 2) then
            !This method do not works with TimeCycle        

            Me%DTForNextEvent       = -null_real
            Me%PredictedDT          = -null_real
            Me%DTForNextDataset     = -null_real
            Me%NextValueForDTPred   = 0.0

            sl => Me%MultiTimeSerie%StationsList            
            do index = 1, Me%MultiTimeSerie%NumberOfSources
            
                if (.not. sl(index)%RemainConstant) then
            
                    if (Now > sl(index)%NextTime) then
                        call ActualizeMultiTimeSerieTimes  (Now, sl(index))
                        call ActualizeMultiTimeSerieValues (sl(index))                    
                    endif
                    
                    !avoid evaluate in cosntruct phase where previous and next time are the same
                    if (Now > sl(index)%PreviousTime) then
                        select case (Me%MultiTimeSerie%DataProcessing)                                           
                        case (Interpolate)
                            !Interpolates Value for current instant
                            call InterpolateValueInTime(Now,                     &
                                                        sl(index)%PreviousTime,  &
                                                        sl(index)%PreviousValue, &
                                                        sl(index)%NextTime,      &
                                                        sl(index)%NextValue,     &
                                                        sl(index)%NewValue)
                        case (Accumulate)
                            sl(index)%NewValue = sl(index)%NextValue / (sl(index)%NextTime - sl(index)%PreviousTime)
                        case (NoProcessing) !like USE_ORIGINAL_VALUES                        
                            sl(index)%NewValue = sl(index)%NextValue
                        case default
                            stop 'ModifyMultiTimeSerie - ModuleFillMatrix - ERR020'
                        end select 
                    endif
                                 
                    if (Me%ValueIsUsedForDTPrediction) then
                        if (Now >= sl(index)%NextEventEnd) then
                            call FindNextEventInMultiTimeSerie (Now, sl(index))
                            if (Me%AccumulateValues .and. (sl(index)%NextValueForDTPred > 0.0)) then
                                sl(index)%NextValueForDTPred = (sl(index)%NextValueForDTPred) / &
                                                               (sl(index)%NextEventEnd - sl(index)%NextEventStart)
                            endif
                        endif

                        if (Now >= sl(index)%NextEventStart .and. Now < sl(index)%NextEventEnd) then
                            sl(index)%DTForNextEvent = 0.0
                        else
                            sl(index)%DTForNextEvent = sl(index)%NextEventStart - Now 
                        endif

                        if (sl(index)%DTForNextEvent > 0.0) then
                            sl(index)%PredictedDT = sl(index)%DTForNextEvent
                        else
                            sl(index)%PredictedDT = sl(index)%NextEventEnd - Now
                        endif   
                    endif
                    
!                    
!                    
!                    
!                    
!                        instant = sl(index)%NextInstant + 1
!                        call GetTimeSerieTimeOfDataset(sl(index)%ObjTimeSerie,      &
!                                                       instant,                     &
!                                                       time,                        &
!                                                       STAT = STAT_CALL)
!                        if (STAT_CALL .NE. SUCCESS_)    &
!                            stop 'ModifyMultiTimeSerie - ModuleFillMatrix - ERR020'                 
!                        call GetTimeSerieValueForIndex (sl(index)%ObjTimeSerie,     &
!                                                        instant,                    &
!                                                        sl(index)%Column,           &
!                                                        NextValue,                  &                                            
!                                                        STAT = STAT_CALL)
!                        if (STAT_CALL .NE. SUCCESS_)    &
!                            stop 'ModifyMultiTimeSerie - ModuleFillMatrix - ERR030'
!                            
!                        if (Me%UseOriginalValues) then
!                            if (NextValue > sl(index)%NextValueForDTPred) &
!                                sl(index)%NextValueForDTPred = NextValue
!                        elseif (Me%AccumulateValues) then
!                            value = NextValue / (time - sl(index)%NextTime)
!                            if (value > sl(index)%NextValueForDTPred) &
!                                sl(index)%NextValueForDTPred = value
!                        else
!                            stop 'ModifyMultiTimeSerie - ModuleFillMatrix - ERR040'
!                        endif                
!                    
!                    else
!                    
!                        sl(index)%NextValueForDTPred = 0.0
!                    
!                    endif    
!                    
!                    select case (Me%MultiTimeSerie%DataProcessing)
!                                           
!                        case (Interpolate)
!                        
!                            !Interpolates Value for current instant
!                            call InterpolateValueInTime(Now,                     &
!                                                        sl(index)%PreviousTime,  &
!                                                        sl(index)%PreviousValue, &
!                                                        sl(index)%NextTime,      &
!                                                        sl(index)%NextValue,     &
!                                                        sl(index)%NewValue)
!
!                        case (Accumulate)
!                        
!!                            if(Now > sl(index)%NextEventEnd) then
!!                                call FindNextEventInMultiTimeSerie (Now, sl(index))
!!                            endif
!
!                            if (Now > sl(index)%NextEventStart) then                    
!                                sl(index)%NewValue = sl(index)%NextValue / (sl(index)%NextTime - sl(index)%PreviousTime)                         
!                            else 
!                                sl(index)%NewValue = 0.0                        
!                            endif
!                            
!                            if (Now >= sl(index)%NextEventStart .and. Now < sl(index)%NextEventEnd) then
!                                sl(index)%DTForNextEvent   = 0.0
!                            elseif (Now == sl(index)%NextEventEnd) then
!                                call FindNextEventInMultiTimeSerie (Now, sl(index))
!                                sl(index)%DTForNextEvent   = sl(index)%NextEventStart - Now
!                            else
!                                sl(index)%DTForNextEvent   = sl(index)%NextEventStart - Now
!                            endif
!                            
!                            if (sl(index)%DTForNextEvent == 0.0) then
!                                sl(index)%PredictedDT = sl(index)%NextEventEnd - Now
!                            else
!                                sl(index)%PredictedDT = sl(index)%DTForNextEvent
!                            endif
!                            
!                            sl(index)%DTForNextDataset = sl(index)%NextTime - Now
!                        
!                        case (NoProcessing) !like USE_ORIGINAL_VALUES
!                        
!                            sl(index)%NewValue = sl(index)%NextValue
!                        
!                        case default
!                        
!                            stop 'ModifyMultiTimeSerie - ModuleFillMatrix - ERR030'
!                        
!                    end select                                

                    if (Me%DTForNextEvent > sl(index)%DTForNextEvent) &
                        Me%DTForNextEvent = sl(index)%DTForNextEvent
                    
                    if (Me%PredictedDT > sl(index)%PredictedDT) &
                        Me%PredictedDT = sl(index)%PredictedDT 

                    if (Me%DTForNextDataset > sl(index)%DTForNextDataset) &
                        Me%DTForNextDataset = sl(index)%DTForNextDataset

                    if ((sl(index)%DTForNextEvent == 0.0) .and. &
                        (Me%NextValueForDTPred < sl(index)%NextValueForDTPred)) then
                        Me%NextValueForDTPred = sl(index)%NextValueForDTPred
                    endif    

                endif

            enddo

            if (Me%Dim == Dim2D) then
            
                id2d => Me%MultiTimeSerie%FillGrid2D
                
                ilb = Me%WorkSize2D%ILB
                iub = Me%WorkSize2D%IUB
                jlb = Me%WorkSize2D%JLB
                jub = Me%WorkSize2D%JUB
            
                do j = jlb, jub
                do i = ilb, iub

                    if (PointsToFill2D(i,j) == WaterPoint) then
                    
                        Me%Matrix2D(i,j) = sl(id2d(i,j))%NewValue
                        
                    endif

                enddo
                enddo

            else
            
                id3d => Me%MultiTimeSerie%FillGrid3D
            
                ilb = Me%WorkSize3D%ILB
                iub = Me%WorkSize3D%IUB
                jlb = Me%WorkSize3D%JLB
                jub = Me%WorkSize3D%JUB
                klb = Me%WorkSize3D%KLB
                kub = Me%WorkSize3D%KUB        
            
                do j = jlb, jub
                do i = ilb, iub
                do k = klb, kub
                
                    if (PointsToFill3D(i,j,k) == WaterPoint) then

                        Me%Matrix3D(i,j,k) = sl(id3d(i,j,k))%NewValue

                    endif

                enddo
                enddo
                enddo
                
            endif

        else

            sl => Me%MultiTimeSerie%StationsList
        
            if (Me%Dim == Dim2D) then
            
                id2d => Me%MultiTimeSerie%FillGrid2D
                
                ilb = Me%WorkSize2D%ILB
                iub = Me%WorkSize2D%IUB
                jlb = Me%WorkSize2D%JLB
                jub = Me%WorkSize2D%JUB        
            
                do j = jlb, jub
                do i = ilb, iub
                
                    if (PointsToFill2D(i,j) == WaterPoint) then
                             
                        if (sl(id2d(i,j))%ValueIsDefined .or. sl(id2d(i,j))%RemainConstant) then
                        
                            Me%Matrix2D(i,j) = sl(id2d(i,j))%NewValue
                        
                        else
                               
                            sl(id2d(i,j))%ValueIsDefined = .true.
                                                         
                            !Gets Value for current Time
                            call GetTimeSerieValue (sl(id2d(i,j))%ObjTimeSerie, Now,                &
                                                    sl(id2d(i,j))%Column,                           &
                                                    Time1, Value1, Time2, Value2, TimeCycle,        &
                                                    STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ModifyMultiTimeSerie - ModuleFillMatrix - ERR010'

                            if (TimeCycle) then
                                                    
                                sl(id2d(i,j))%NewValue = Value1
                        
                            else
                            
                                select case (Me%MultiTimeSerie%DataProcessing)
                                
                                    case (Interpolate)
                                    
                                        !Interpolates Value for current instant
                                        call InterpolateValueInTime(Now, Time1, Value1, Time2, Value2, sl(id2d(i,j))%NewValue)
                                    
                                    case (Accumulate)
                                    
                                        sl(id2d(i,j))%NewValue = Value2 / (Time2 - Time1)
                                    
                                        !write (*,*) 'Multi NewValue = ', sl(id2d(i,j))%NewValue
                                    
                                        call GetTimeSerieDTForNextEvent (sl(id2d(i,j))%ObjTimeSerie,                    &
                                                                         sl(id2d(i,j))%NewValue, sl(id2d(i,j))%Column,  &
                                                                         Now, PredictedDT, DTForNextEvent,              &
                                                                         STAT  = STAT_CALL)
                                        if (STAT_CALL /= SUCCESS_) stop 'ModifyMultiTimeSerie - ModuleFillMatrix - ERR020'
                                    
                                        if (compare_dt) then
                                            if (Me%PredictedDT > PredictedDT) then
                                                Me%PredictedDT = PredictedDT
                                            endif
                                            
                                            if (Me%DTForNextEvent > DTForNextEvent) then
                                                Me%DTForNextEvent = DTForNextEvent
                                            endif
                                        else
                                            Me%PredictedDT    = PredictedDT
                                            Me%DTForNextEvent = DTForNextEvent
                                            compare_dt        = .true.
                                        endif
                                    
                                    case (NoProcessing)
                                    
                                        sl(id2d(i,j))%NewValue = Value2
                                    
                                    case default
                                    
                                        stop 'ModifyMultiTimeSerie - ModuleFillMatrix - ERR030'
                                        
                                end select
                            endif
                            
                            Me%Matrix2D(i,j) = sl(id2d(i,j))%NewValue
                            
                        endif
                        
                    endif
                    
                enddo
                enddo        
            
            else
            
                id3d => Me%MultiTimeSerie%FillGrid3D
            
                ilb = Me%WorkSize3D%ILB
                iub = Me%WorkSize3D%IUB
                jlb = Me%WorkSize3D%JLB
                jub = Me%WorkSize3D%JUB
                klb = Me%WorkSize3D%KLB
                kub = Me%WorkSize3D%KUB        
            
                do j = jlb, jub
                do i = ilb, iub
                do k = klb, kub
                
                    if (PointsToFill3D(i,j,k) == WaterPoint) then
                             
                        if (sl(id3d(i,j,k))%ValueIsDefined) then
                        
                            Me%Matrix3D(i,j,k) = sl(id3d(i,j,k))%NewValue
                        
                        else
                               
                            sl(id3d(i,j,k))%ValueIsDefined = .true.
                                                            
                            !Gets Value for current Time
                            call GetTimeSerieValue (sl(id3d(i,j,k))%ObjTimeSerie, Now,              &
                                                    sl(id3d(i,j,k))%Column,                         &
                                                    Time1, Value1, Time2, Value2, TimeCycle,        &
                                                    STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ModifyMultiTimeSerie - ModuleFillMatrix - ERR040'

                            if (TimeCycle) then
                                                    
                                sl(id3d(i,j,k))%NewValue = Value1
                        
                            else
                            
                                select case (Me%MultiTimeSerie%DataProcessing)
                                
                                    case (Interpolate)
                                    
                                        !Interpolates Value for current instant
                                        call InterpolateValueInTime(Now, Time1, Value1, Time2, Value2, sl(id3d(i,j,k))%NewValue)
                                    
                                    case (Accumulate)
                                    
                                        sl(id3d(i,j,k))%NewValue = Value2 / (Time2 - Time1)
                                    
                                        call GetTimeSerieDTForNextEvent (sl(id3d(i,j,k))%ObjTimeSerie,                      &
                                                                         sl(id3d(i,j,k))%NewValue, sl(id3d(i,j,k))%Column,  &
                                                                         Now, PredictedDT, DTForNextEvent,                  &
                                                                         STAT  = STAT_CALL)
                                        if (STAT_CALL /= SUCCESS_) stop 'ModifyMultiTimeSerie - ModuleFillMatrix - ERR050'
                                    
                                        if (compare_dt) then
                                            if (Me%PredictedDT > PredictedDT) then
                                                Me%PredictedDT = PredictedDT
                                            endif
                                            
                                            if (Me%DTForNextEvent > DTForNextEvent) then
                                                Me%DTForNextEvent = DTForNextEvent
                                            endif
                                        else
                                            Me%PredictedDT    = PredictedDT
                                            Me%DTForNextEvent = DTForNextEvent
                                            compare_dt        = .true.
                                        endif
                                    
                                    case (NoProcessing)
                                    
                                        sl(id3d(i,j,k))%NewValue = Value2
                                    
                                    case default
                                    
                                        stop 'ModifyMultiTimeSerie - ModuleFillMatrix - ERR060'
                                        
                                end select
                            endif
                            
                            Me%Matrix3D(i,j,k) = sl(id3d(i,j,k))%NewValue
                            
                        endif
                        
                    endif
                    
                enddo
                enddo
                enddo
            
            endif
            
            do index = 1, Me%MultiTimeSerie%NumberOfSources
                sl(index)%ValueIsDefined = .false.
            enddo

        endif

        !----------------------------------------------------------------------

    end subroutine ModifyMultiTimeSerie
    
    !--------------------------------------------------------------------------
    
    subroutine ActualizeMultiTimeSerieTimes (ActualTime, Station)
    
        !Arguments-------------------------------------------------------------        
        Type (T_Time), intent(IN)                       :: ActualTime
        Type (T_Station)                                :: Station

        !Local-----------------------------------------------------------------      
        integer                                         :: STAT_CALL
        
        !----------------------------------------------------------------------                     
        
        do while (Station%NextTime < ActualTime)  
                          
            !write (*,*) 'Station%NumberOfInstants: ', Station%NumberOfInstants
                          
            if (Station%NextInstant < Station%NumberOfInstants) then
            
                Station%PreviousInstant = Station%NextInstant
                Station%PreviousTime    = Station%NextTime
                Station%NextInstant     = Station%NextInstant + 1
                
            else
            
                stop 'ActualizeMultiTimeSerieTimes - ModuleFillMatrix - ERR010'
                
            endif
                  
            call GetTimeSerieTimeOfDataset(Station%ObjTimeSerie,    &
                                           Station%NextInstant,     &
                                           Station%NextTime,        &
                                           STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)    &
                stop 'ActualizeMultiTimeSerieTimes - ModuleFillMatrix - ERR020'  
           
        enddo
 
        !----------------------------------------------------------------------
        
    end subroutine ActualizeMultiTimeSerieTimes 
    
    !--------------------------------------------------------------------------
    
    subroutine ActualizeMultiTimeSerieValues (Station)
    
        !Arguments-------------------------------------------------------------
        Type (T_Station)                                :: Station

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_
        
        !----------------------------------------------------------------------    
    
        call GetTimeSerieValueForIndex (Station%ObjTimeSerie,      &
                                        Station%PreviousInstant,   &
                                        Station%Column,            &
                                        Station%PreviousValue,     &
                                        STAT = STAT_)
        if (STAT_ /= SUCCESS_) stop 'ActualizeMultiTimeSerieValues - ModuleFillMatrix - ERR010'            

        call GetTimeSerieValueForIndex (Station%ObjTimeSerie,      &
                                        Station%NextInstant,       &
                                        Station%Column,            &
                                        Station%NextValue,         &                                            
                                        STAT = STAT_)
        if (STAT_ /= SUCCESS_) stop 'ActualizeMultiTimeSerieValues - ModuleFillMatrix - ERR020' 
        
        !----------------------------------------------------------------------    
    
    end subroutine ActualizeMultiTimeSerieValues
       
    !-------------------------------------------------------------------------- 
    
    subroutine FindNextEventInMultiTimeSerie(Now, Station)
    
        !Arguments-------------------------------------------------------------        
        Type (T_Time), intent(IN)                       :: Now
        Type (T_Station)                                :: Station

        !Local-----------------------------------------------------------------
        real                                            :: instant_value
        integer                                         :: instant, STAT_CALL

        !----------------------------------------------------------------------

        if (Station%NextInstant < Station%NumberOfInstants) then
        
            if (Now > Station%NextEventEnd) then

                Station%NextEventStart  = Station%PreviousTime
                instant                 = Station%NextInstant
                Station%NextEventEnd    = Station%NextTime                
                instant_value           = Station%NextValue

            else
            
                Station%NextEventStart  = Station%NextTime
                instant                 = Station%NextInstant + 1
                call GetTimeSerieTimeOfDataset(Station%ObjTimeSerie,   &
                                               instant,                &
                                               Station%NextEventEnd,   &
                                               STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'FindNextEventInMultiTimeSerie - ModuleFillMatrix - ERR010'                 
                call GetTimeSerieValueForIndex (Station%ObjTimeSerie,  &
                                                instant,               &
                                                Station%Column,        &
                                                instant_value,         &                                            
                                                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'FindNextEventInMultiTimeSerie - ModuleFillMatrix - ERR020'   
            
            endif
                             
            do while (instant_value <= Me%MinForDTDecrease .and. instant < Station%NumberOfInstants)                
                instant = instant + 1  
                    
                call GetTimeSerieValueForIndex (Station%ObjTimeSerie,  &
                                                instant,               &
                                                Station%Column,        &
                                                instant_value,         &                                            
                                                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'FindNextEventInMultiTimeSerie - ModuleFillMatrix - ERR030'                                
            enddo

            if (instant_value > Me%MinForDTDecrease) then
            
                call GetTimeSerieTimeOfDataset(Station%ObjTimeSerie,   &
                                               Instant - 1,            &
                                               Station%NextEventStart, &
                                               STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'FindNextEventInMultiTimeSerie - ModuleFillMatrix - ERR040' 
                
                call GetTimeSerieTimeOfDataset(Station%ObjTimeSerie,   &
                                               instant,                &
                                               Station%NextEventEnd,   &
                                               STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'FindNextEventInMultiTimeSerie - ModuleFillMatrix - ERR050'              
                    
            else
            
                call GetTimeSerieTimeOfDataset(Station%ObjTimeSerie,   &
                                               instant,                &
                                               Station%NextEventStart, &
                                               STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'FindNextEventInMultiTimeSerie - ModuleFillMatrix - ERR060'                                 
                        
                Station%NextEventEnd = Station%NextEventStart
                instant_value        = 0.0
                                    
            endif 
            
        else
        
            Station%NextEventStart = Station%NextTime
            Station%NextEventEnd   = Station%NextTime                

        endif
        
        Station%NextValueForDTPred = instant_value
!            Station%NextEventStart = Station%NextTime
!        
!            instant = Station%NextInstant + 1
!            call GetTimeSerieTimeOfDataset(Station%ObjTimeSerie,   &
!                                           Instant,                &
!                                           Station%NextEventEnd,   &
!                                           STAT_CALL)
!            if (STAT_CALL .NE. SUCCESS_) &
!                stop 'FindNextEventInMultiTimeSerie - ModuleFillMatrix - ERR010'
!                        
!            call GetTimeSerieValueForIndex (Station%ObjTimeSerie,  &
!                                            instant,               &
!                                            Station%Column,        &
!                                            instant_value,         &                                            
!                                            STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) &
!                stop 'FindNextEventInMultiTimeSerie - ModuleFillMatrix - ERR020'     
!
!            do while (instant_value <= Me%MinForDTDecrease .and. instant < Station%NumberOfInstants)                
!                instant = instant + 1  
!                    
!                call GetTimeSerieValueForIndex (Station%ObjTimeSerie,       &
!                                                instant,                    &
!                                                Station%Column,             &
!                                                instant_value,              &                                            
!                                                STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) &
!                    stop 'FindNextEventInMultiTimeSerie - ModuleFillMatrix - ERR030'                                
!            enddo
!
!            if (instant_value > Me%MinForDTDecrease) then
!            
!                call GetTimeSerieTimeOfDataset(Station%ObjTimeSerie,        &
!                                               instant - 1,                 &
!                                               Station%NextEventStart,      &
!                                               STAT_CALL)
!                if (STAT_CALL .NE. SUCCESS_) &
!                    stop 'FindNextEventInMultiTimeSerie - ModuleFillMatrix - ERR040' 
!                
!                call GetTimeSerieTimeOfDataset(Station%ObjTimeSerie,        &
!                                               instant,                     &
!                                               Station%NextEventEnd,        &
!                                               STAT_CALL)
!                if (STAT_CALL .NE. SUCCESS_) &
!                    stop 'FindNextEventInMultiTimeSerie - ModuleFillMatrix - ERR050' 
!                    
!            else
!            
!                call GetTimeSerieTimeOfDataset(Station%ObjTimeSerie,        &
!                                               instant,                     &
!                                               Station%NextEventStart,      &
!                                               STAT_CALL)
!                if (STAT_CALL .NE. SUCCESS_) &
!                    stop 'FindNextEventInMultiTimeSerie - ModuleFillMatrix - ERR060'                                 
!                        
!                Station%NextEventEnd = Station%NextEventStart
!                                    
!            endif
!
!            Station%DTForNextEvent = Station%NextEventStart - ActualTime
!        
!        endif            

        !----------------------------------------------------------------------
    
    end subroutine FindNextEventInMultiTimeSerie
    
    !--------------------------------------------------------------------------

    real function TimeSerieValue(ObjTimeSerie, Now, TimeSerieColumn)    
        !Arguments--------------------------------------------------------------
        integer                                         :: ObjTimeSerie, TimeSerieColumn
        type (T_Time)                                   :: Now
        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL
        type (T_Time)                                   :: Time1, Time2
        real                                            :: Value1, Value2
        logical                                         :: TimeCycle

        !Begin------------------------------------------------------------------



        !Gets Value for current Time
        call GetTimeSerieValue (ObjTimeSerie, Now, TimeSerieColumn,                     &
                                Time1, Value1, Time2, Value2, TimeCycle,                &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'TimeSerieValue - ModuleFillMatrix - ERR10'
     
        if (TimeCycle) then
            TimeSerieValue = Value1

        else

            !Interpolates Value for current instant
            call InterpolateValueInTime(Now, Time1, Value1, Time2, Value2, TimeSerieValue)

        endif

    end function TimeSerieValue

    !--------------------------------------------------------------------------
    
    subroutine ModifyAnalyticWave(PointsToFill2D, PointsToFill3D)
    
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer, optional    :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional :: PointsToFill3D             

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: ilb, iub, jlb, jub, klb, kub
        integer                                         :: i, j, k  
        type (T_Time)                                   :: Now, RefDate
        real(8)                                         :: Amplitude, Period, AverageValue
        real(8)                                         :: TimeSeconds, T1, T2, T3, Dir, A
        real(8)                                         :: StartPeriod, T4, SlowStartPeriod
        logical                                         :: SlowStartON
        integer                                         :: WaveType, n     
        
        !----------------------------------------------------------------------
        
        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyAnalyticWave - ModuleFillMatrix - ERR010'
        
        Period          = Me%AnalyticWave%Period
        WaveType        = Me%AnalyticWave%WaveType
        AverageValue    = Me%AnalyticWave%AverageValue
        SlowStartON     = Me%AnalyticWave%SlowStartON
        SlowStartPeriod = Me%AnalyticWave%SlowStartPeriod
        
        
        call SetDate(RefDate, Year = 2018, Month = 1, Day = 1, Hour = 0, Minute = 0, Second = 0)        
        
        TimeSeconds   = Now - RefDate
        StartPeriod   = Now - Me%BeginTime 
        
        
        !TimeSeconds   = Now - Me%BeginTime

        if (Me%Dim == Dim2D) then
            
            ILB = Me%WorkSize2D%ILB
            IUB = Me%WorkSize2D%IUB
            JLB = Me%WorkSize2D%JLB
            JUB = Me%WorkSize2D%JUB

            do j = JLB, JUB
            do i = ILB, IUB

                if (PointsToFill2D(i,j) == WaterPoint) then
                
                    Me%Matrix2D(i,j) = AverageValue                    
                
                endif
                
            enddo
            enddo
            
            
            do n = 1, Me%AnalyticWave%EnteringCell%nCells

                                                    
                !X  = Me%AnalyticWave%X2D(i, j)

                i = Me%AnalyticWave%EnteringCell%i(n)
                j = Me%AnalyticWave%EnteringCell%j(n)
                
                if (Me%AnalyticWave%CellType(i, j) /= EnteringWaveCell_) then
                    stop 'ModifyAnalyticWave - ModuleFillMatrix - ERR020'
                endif

                Amplitude     = Me%AnalyticWave%AmpAux  (i, j)
                !WaveCelerity  = Me%AnalyticWave%Celerity(i, j
                
                T1 = TimeSeconds
                !T2 = X / WaveCelerity
                T2 = Me%AnalyticWave%EnteringCell%TimeLag(n)
                T3 = T1 - T2
                T4 = StartPeriod - T2
                A  = Amplitude
                Dir = Me%AnalyticWave%Direction  
            
                
                if (SlowStartON) then

                    if (T4 >= 0) then

                        if (T4 < SlowStartPeriod) then
                            A = A * T4 / SlowStartPeriod

                        endif

                    else
                    
                        A = 0.

                    endif                        

                endif
                
                Me%Matrix2D(i,j) =  AverageValue + ComputeWaveAnalytic1D (A, Period, T3, WaveType, Dir)

            enddo

         else
            
            ILB = Me%WorkSize3D%ILB
            IUB = Me%WorkSize3D%IUB
            JLB = Me%WorkSize3D%JLB
            JUB = Me%WorkSize3D%JUB
            KLB = Me%WorkSize3D%KLB
            KUB = Me%WorkSize3D%KUB        

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
            
                if (PointsToFill3D(i,j,k) == WaterPoint) then
!
!                    X = Me%AnalyticWave%X2D(i, j)
!
!                    if (Me%AnalyticWave%CellType(i, j) == EnteringWaveCell_) then
!                        Me%Matrix3D(i,j, k) = AverageValue + ComputeWaveAnalytic2D (X, Amplitude, Period, TimeSeconds, WaveType)
!                    else
!                        Me%Matrix3D(i,j, k) = AverageValue
!                    endif 
!
                endif

            enddo
            enddo
            enddo
                
        endif
                
    
    end subroutine ModifyAnalyticWave    
    
    !--------------------------------------------------------------------------   
    
    real(8) function ComputeWaveAnalytic1D (Amplitude, Period, T, WaveType, Dir)     

        !Arguments------------------------------------------------------------
        real(8)            :: Amplitude, Period, T, Dir
        integer            :: WaveType

        !Local----------------------------------------------------------------

        !Begin----------------------------------------------------------------
        !wave types
!            SineWaveSeaLevel_          = 1
!            CnoidalWaveSeaLevel_       = 2
!            SolitartyWaveSeaLevel_     = 3
!            SineWaveVelX_              = 4
!            CnoidalWaveVelX_           = 5
!            SolitartyWaveVelX_         = 6
!            SineWaveVelY_              = 7
!            CnoidalWaveVelY_           = 8
!            SolitartyWaveVelY_         = 9

        if      (WaveType == SineWaveSeaLevel_) then
            ComputeWaveAnalytic1D = Amplitude * sin(2.*Pi*T/Period)

        else if (WaveType == SineWaveVelX_    ) then
            ComputeWaveAnalytic1D = Amplitude * sin(2.*Pi*T/Period) * cos(Dir)

        else if (WaveType == SineWaveVelY_    ) then
            ComputeWaveAnalytic1D = Amplitude * sin(2.*Pi*T/Period) * sin(Dir)

        endif

    
    end function ComputeWaveAnalytic1D
    
    !--------------------------------------------------------------------------       

    subroutine ModifyHDFInput3D(PointsToFill3D, CurrentHDF, Generic_4D_Value_)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        type(T_Field4D)                                 :: CurrentHDF        
        real, optional                                  :: Generic_4D_Value_

        !Local----------------------------------------------------------------

        !Begin----------------------------------------------------------------

        if (CurrentHDF%Generic4D%ON) then

            if (.not. present(Generic_4D_Value_)) &
                stop 'ModifyHDFInput3D - ModuleFillMatrix - ERR010'
            call ModifyHDFInput3DGeneric4D(PointsToFill3D, Generic_4D_Value_, CurrentHDF)

        else
            
            call ModifyHDFInput3DTime(PointsToFill3D, CurrentHDF)

        endif

    end subroutine ModifyHDFInput3D

    !----------------------------------------------------------------------------

    subroutine ModifyHDFInput3DHarmonics(CurrentHDF)
        
        !Arguments------------------------------------------------------------
        type(T_Field4D)                                 :: CurrentHDF
        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        type (T_Time)                                   :: Now, CurrentTime

        !Begin----------------------------------------------------------------

        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadNewField - ModuleFillMatrix - ERR10'

        !Backtracking time inversion is also done in the ModuleField4D    
        if (Me%BackTracking .and. .not. CurrentHDF%Field4D) then  
            call BacktrackingTime(Now)
        else   
            Now = CurrentTime
        endif 
    
        if (CurrentHDF%SpatialInterpolON) then
        
            call ModifyField4DInterpol(CurrentTime      = Now,                          & 
                                       Matrix3D         = Me%Matrix3D,                  &
                                       CurrentHDF       = CurrentHDF)
        
        else
            call ModifyField4D(Field4DID        = CurrentHDF%ObjField4D,                &
                               PropertyIDNumber = Me%PropertyID%IDNumber,               & 
                               CurrentTime      = Now,                                  & 
                               Matrix3D         = Me%Matrix3D,                          &
                               STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR20'  
        endif
            
    end subroutine ModifyHDFInput3DHarmonics
        
    !----------------------------------------------------------------------------

    subroutine ModifyHDFInput3DStandard(PointsToFill3D, CurrentHDF)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        type(T_Field4D)                                 :: CurrentHDF
        !Local----------------------------------------------------------------
        integer                                         :: n, i, j, k
        type (T_Time)                                   :: Now

        !Begin----------------------------------------------------------------

i2:     if (ReadNewField(Now,n, CurrentHDF))then 
    
i4:         if (n==1) then
                call SetMatrixValue(CurrentHDF%PreviousField3D, Me%WorkSize3D, CurrentHDF%NextField3D)
            else i4
                call ReadHDF5Values3D(CurrentHDF%PreviousInstant, CurrentHDF%PreviousField3D, CurrentHDF)
                
                !limit maximum values
                do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
                do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB

#ifndef _NOT_IEEE_ARITHMETIC
                    if (ieee_is_nan (CurrentHDF%PreviousField3D(i,j,k)))                    &
                        CurrentHDF%PreviousField3D       (i,j,k) = FillValueReal 
#endif
                
                    if (abs(CurrentHDF%PreviousField3D(i,j,k)) > abs(FillValueReal))        &
                            CurrentHDF%PreviousField3D(i,j,k) = FillValueReal

                enddo
                enddo
                enddo                
                
            endif i4

            call ReadHDF5Values3D(CurrentHDF%NextInstant, CurrentHDF%NextField3D, CurrentHDF)
            
            !limit maximum values
            do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
            
#ifndef _NOT_IEEE_ARITHMETIC
                if (ieee_is_nan (CurrentHDF%NextField3D (i,j,k)))                           &
                    CurrentHDF%NextField3D        (i,j,k) = FillValueReal 
#endif

                if (abs(CurrentHDF%NextField3D    (i,j,k)) > abs(FillValueReal))            &
                        CurrentHDF%NextField3D    (i,j,k) = FillValueReal

            enddo
            enddo
            enddo                
            
        end if i2

i3:     if (CurrentHDF%PreviousInstant /= CurrentHDF%NextInstant) then
                                            
i5:         if (Me%PreviousInstantValues) then
            
                Me%Matrix3D = CurrentHDF%PreviousField3D

            else i5
            
                if (Me%PropertyID%IsAngle) then            
                    call InterpolateAngle3DInTime (ActualTime       = Now,                         &
                                                   Size             = Me%WorkSize3D,               &
                                                   Time1            = CurrentHDF%PreviousTime,     &
                                                   Matrix1          = CurrentHDF%PreviousField3D,  &
                                                   Time2            = CurrentHDF%NextTime,         &
                                                   Matrix2          = CurrentHDF%NextField3D,      &
                                                   MatrixOut        = Me%Matrix3D,                 &
                                                   PointsToFill3D   = PointsToFill3D)
                else
                    call InterpolateMatrix3DInTime(ActualTime       = Now,                         &
                                                   Size             = Me%WorkSize3D,               &
                                                   Time1            = CurrentHDF%PreviousTime,     &
                                                   Matrix1          = CurrentHDF%PreviousField3D,  &
                                                   Time2            = CurrentHDF%NextTime,         &
                                                   Matrix2          = CurrentHDF%NextField3D,      &
                                                   MatrixOut        = Me%Matrix3D,                 &
                                                   PointsToFill3D   = PointsToFill3D)
                endif                                                   
            endif i5
        else i3
            
            !Prev and next are equal (last instant?)
            call SetMatrixValue(Me%Matrix3D, Me%WorkSize3D, CurrentHDF%NextField3D)

        endif i3



    end subroutine ModifyHDFInput3DStandard
            
    !----------------------------------------------------------------------------

    subroutine ModifyHDFInput3DTime(PointsToFill3D, CurrentHDF)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        type(T_Field4D)                                 :: CurrentHDF
        !Local----------------------------------------------------------------

        !Begin----------------------------------------------------------------

i1:     if (CurrentHDF%Field4D .and. CurrentHDF%HarmonicsOn) then

            call ModifyHDFInput3DHarmonics(CurrentHDF)

        else i1 
        
            call ModifyHDFInput3DStandard(PointsToFill3D, CurrentHDF)
        
        endif i1

    end subroutine ModifyHDFInput3DTime


    !--------------------------------------------------------------------------


    subroutine ModifyHDFInput3DGeneric4D(PointsToFill3D, Generic_4D_Value_, CurrentHDF)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        real                                            :: Generic_4D_Value_
        type(T_Field4D)                                 :: CurrentHDF
        !Local----------------------------------------------------------------
        integer                                         :: PrevI, NextI, i, j, k
        !Begin----------------------------------------------------------------
        
i1:     if (.not.(CurrentHDF%Previous4DValue <= Generic_4D_Value_ .and.                     &
                  CurrentHDF%Next4DValue     >= Generic_4D_Value_)) then
            !Found new limits
            PrevI              = 1
            NextI              = 2
            CurrentHDF%Next4DValue = HDF5Generic4DInstant(1, CurrentHDF)
            do 

                CurrentHDF%Previous4DValue  = CurrentHDF%Next4DValue
                CurrentHDF%Next4DValue      = HDF5Generic4DInstant(NextI, CurrentHDF)

                if (CurrentHDF%Previous4DValue <= Generic_4D_Value_ .and.                  &
                    CurrentHDF%Next4DValue     >= Generic_4D_Value_) then
                    exit
                endif

                if (NextI > CurrentHDF%NumberOfInstants) then

                    write(*,*)
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'Time instants inconsistency.'
                    stop      'ModifyHDFInput3DGeneric4D - ModuleFillMatrix - ERR10'

                endif

                PrevI = NextI
                NextI = NextI + 1

            enddo

            CurrentHDF%NextInstant     = NextI
            CurrentHDF%PreviousInstant = PrevI  


            !call SetMatrixValue(CurrentHDF%PreviousField3D, Me%WorkSize3D, CurrentHDF%NextField3D)

            call ReadHDF5Values3D(CurrentHDF%PreviousInstant, CurrentHDF%PreviousField3D, CurrentHDF)

            call ReadHDF5Values3D(CurrentHDF%NextInstant,     CurrentHDF%NextField3D, CurrentHDF)
            
            !limit maximum values
            do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB

#ifndef _NOT_IEEE_ARITHMETIC
                if (ieee_is_nan (CurrentHDF%PreviousField3D(i,j,k)))                        &
                    CurrentHDF%PreviousField3D       (i,j,k) = FillValueReal 
#endif
            
                if (abs(CurrentHDF%PreviousField3D(i,j,k)) > abs(FillValueReal))            &
                        CurrentHDF%PreviousField3D(i,j,k) = FillValueReal
      
#ifndef _NOT_IEEE_ARITHMETIC                  
                if (ieee_is_nan (CurrentHDF%NextField3D    (i,j,k)))                        &
                    CurrentHDF%NextField3D           (i,j,k) = FillValueReal 
#endif                        
                
                if (abs(CurrentHDF%NextField3D    (i,j,k)) > abs(FillValueReal))            &
                        CurrentHDF%NextField3D    (i,j,k) = FillValueReal
            enddo
            enddo
            enddo                
            
                    
        endif i1
        

        if (CurrentHDF%PreviousInstant /= CurrentHDF%NextInstant) then

            call InterpolateLinearyMatrix3D(X                = Generic_4D_Value_,       &
                                            Size             = Me%WorkSize3D,           &
                                            X1               = CurrentHDF%Previous4DValue,  &
                                            Matrix1          = CurrentHDF%PreviousField3D,  &
                                            X2               = CurrentHDF%Next4DValue,      &
                                            Matrix2          = CurrentHDF%NextField3D,      &
                                            MatrixOut        = Me%Matrix3D,             &
                                            PointsToFill3D   = PointsToFill3D)

        else
            
            !Prev and next are equal (last instant?)
            call SetMatrixValue(Me%Matrix3D, Me%WorkSize3D, CurrentHDF%NextField3D)

        endif


    end subroutine ModifyHDFInput3DGeneric4D


    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------


    subroutine ModifyHDFInput2D(PointsToFill2D, CurrentHDF, Generic_4D_Value_)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        type(T_Field4D)                                 :: CurrentHDF        
        real, optional                                  :: Generic_4D_Value_
        !Local----------------------------------------------------------------

        !Begin----------------------------------------------------------------

        if (CurrentHDF%Generic4D%ON) then

            if (present(Generic_4D_Value_)) then
                call ModifyHDFInput2DGeneric4D(PointsToFill2D, Generic_4D_Value_, CurrentHDF)
            endif
        else

            call ModifyHDFInput2DTime(PointsToFill2D, CurrentHDF)

        endif

    end subroutine ModifyHDFInput2D

    !----------------------------------------------------------------------------

    subroutine ModifyHDFInput2DHarmonics(PointsToFill2D, CurrentHDF)
        
        !Arguments------------------------------------------------------------
        integer,    dimension(:,:), pointer             :: PointsToFill2D
        type(T_Field4D)                                 :: CurrentHDF
        !Local----------------------------------------------------------------
        type (T_Time)                                   :: Now, CurrentTime
        integer                                         :: STAT_CALL

        !Begin----------------------------------------------------------------

        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyHDFInput2DHarmonics - ModuleFillMatrix - ERR10'

        !Backtracking time inversion is also done in the ModuleField4D    
        if (Me%BackTracking .and. .not. CurrentHDF%Field4D) then  
            call BacktrackingTime(Now)
        else   
            Now = CurrentTime
        endif         

            
        if (CurrentHDF%SpatialInterpolON) then
        
            if (Me%BackTracking) then  
                stop 'ModifyHDFInput2DHarmonics - ModuleFillMatrix - ERR20'
            endif
            
            if (Now == Me%BeginTime) then            
                
                CurrentHDF%NextTime =  Me%BeginTime
                
                call ModifyField4DInterpol(CurrentTime      = CurrentHDF%NextTime,      & 
                                            Matrix2D        = CurrentHDF%NextField2D,   &
                                            CurrentHDF      = CurrentHDF)            
            
            endif
        
            if (Now >= CurrentHDF%NextTime) then

                CurrentHDF%PreviousTime = CurrentHDF%NextTime
                CurrentHDF%NextTime     = CurrentHDF%NextTime + CurrentHDF%HarmonicsDT            
        
                call ModifyField4DInterpol(CurrentTime      = CurrentHDF%NextTime,      & 
                                            Matrix2D        = Me%Matrix2D,              &
                                            CurrentHDF      = CurrentHDF)
                                               
                CurrentHDF%PreviousField2D(:,:) = CurrentHDF%NextField2D(:,:)
                CurrentHDF%NextField2D    (:,:) = Me%Matrix2D           (:,:)                                               
                    
            endif                    
                                           
            call InterpolateMatrix2DInTime(ActualTime        = Now,                     &
                                            Size             = Me%WorkSize2D,           &
                                            Time1            = CurrentHDF%PreviousTime, &
                                            Matrix1          = CurrentHDF%PreviousField2D,&
                                            Time2            = CurrentHDF%NextTime,     &
                                            Matrix2          = CurrentHDF%NextField2D,  &
                                            MatrixOut        = Me%Matrix2D,             &
                                            PointsToFill2D   = PointsToFill2D)                                           
        
        else
        

            call ModifyField4D(Field4DID        = CurrentHDF%ObjField4D,                &
                                PropertyIDNumber = Me%PropertyID%IDNumber,              & 
                                CurrentTime      = Now,                                 & 
                                Matrix2D         = Me%Matrix2D,                         &
                                STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ModifyHDFInput2DHarmonics - ModuleFillMatrix - ERR30'  

        endif

        
        


    end subroutine ModifyHDFInput2DHarmonics
        
    !----------------------------------------------------------------------------        


    subroutine ModifyHDFInput2DRainType(PointsToFill2D, CurrentHDF)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        type(T_Field4D)                                 :: CurrentHDF
        !Local----------------------------------------------------------------
        type (T_Time)                                   :: Now
        integer                                         :: STAT_CALL
        real                                            :: aux

        !Begin----------------------------------------------------------------

        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyHDFInput2DRainType - ModuleFillMatrix - ERR010'
        
i3:     if (Now > CurrentHDF%NextTime) then
            call ActualizeHDFTimes (Now, CurrentHDF)
            call ActualizeHDFValues (CurrentHDF, CurrentHDF%PreviousInstant, CurrentHDF%PreviousField2D)
            call ActualizeHDFValues (CurrentHDF, CurrentHDF%NextInstant, CurrentHDF%NextField2D)
        endif i3
        
        !avoid evaluate in cosntruct phase where previous and next time are the same
i4:     if (Now > CurrentHDF%PreviousTime) then
i5:         if (Me%UseOriginalValues) then
                Me%Matrix2D = CurrentHDF%NextField2D                    
            else if (Me%AccumulateValues) then i5
                Me%Matrix2D = CurrentHDF%NextField2D / (CurrentHDF%NextTime - CurrentHDF%PreviousTime)
            else i5
                !Interpolates the two matrixes in time
                if (Me%PropertyID%IsAngle) then
                    call InterpolateAngle2DInTime (ActualTime       = Now,                         &
                                                   Size             = Me%WorkSize2D,               &
                                                   Time1            = CurrentHDF%PreviousTime,     &
                                                   Matrix1          = CurrentHDF%PreviousField2D,  &
                                                   Time2            = CurrentHDF%NextTime,         &
                                                   Matrix2          = CurrentHDF%NextField2D,      &
                                                   MatrixOut        = Me%Matrix2D,                 &
                                                   PointsToFill2D   = PointsToFill2D)                              
                else
                    call InterpolateMatrix2DInTime(ActualTime       = Now,                         &
                                                   Size             = Me%WorkSize2D,               &
                                                   Time1            = CurrentHDF%PreviousTime,     &
                                                   Matrix1          = CurrentHDF%PreviousField2D,  &
                                                   Time2            = CurrentHDF%NextTime,         &
                                                   Matrix2          = CurrentHDF%NextField2D,      &
                                                   MatrixOut        = Me%Matrix2D,                 &
                                                   PointsToFill2D   = PointsToFill2D)                
                endif                                                   
            endif i5
        endif i4          
        
i6:     if (Me%ValueIsUsedForDTPrediction) then
i7:         if (Now >= CurrentHDF%NextEventEnd) then                             
                call FindNextEventInHDF (Now, CurrentHDF)
i8:             if (Me%AccumulateValues .and. (CurrentHDF%NextValueForDTPred > 0.0)) then
                    aux = CurrentHDF%NextEventEnd - CurrentHDF%NextEventStart
                    if (aux > 0.) then
                        CurrentHDF%NextValueForDTPred = CurrentHDF%NextValueForDTPred / aux
                    endif                                
                endif i8
            endif i7
            
i9:         if (Now >= CurrentHDF%NextEventStart .and. Now < CurrentHDF%NextEventEnd) then
                CurrentHDF%DTForNextEvent = 0.0
            else i9
                CurrentHDF%DTForNextEvent = CurrentHDF%NextEventStart - Now 
            endif i9
            
i10:        if (CurrentHDF%DTForNextEvent > 0.0) then
                CurrentHDF%PredictedDT = CurrentHDF%DTForNextEvent                    
            else i10
                CurrentHDF%PredictedDT = CurrentHDF%NextEventEnd - Now
            endif i10                
        endif  i6


    end subroutine ModifyHDFInput2DRainType
        
    !----------------------------------------------------------------------------    


    subroutine ModifyHDFInput2DStandard(PointsToFill2D, CurrentHDF)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        type(T_Field4D)                                 :: CurrentHDF
        !Local----------------------------------------------------------------
        integer                                         :: n, i, j
        type (T_Time)                                   :: Now

        !Begin----------------------------------------------------------------

i11:    if (ReadNewField(Now,n, CurrentHDF))then
i12:        if (n==1) then 
                call SetMatrixValue(CurrentHDF%PreviousField2D, Me%WorkSize2D, CurrentHDF%NextField2D, PointsToFill2D)
            else i12
                call ReadHDF5Values2D(CurrentHDF%PreviousInstant, CurrentHDF%PreviousField2D, CurrentHDF)
                
                !limit maximum values
                do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
                do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
            
#ifndef _NOT_IEEE_ARITHMETIC                    
                    if (ieee_is_nan (CurrentHDF%PreviousField2D(i,j)))                      &
                        CurrentHDF%PreviousField2D(i,j) = FillValueReal                
#endif
                    
                    if (abs(CurrentHDF%PreviousField2D(i,j)) > abs(FillValueReal))          &
                        CurrentHDF%PreviousField2D(i,j) = FillValueReal
                enddo
                enddo   
                    
            endif i12

            call ReadHDF5Values2D(CurrentHDF%NextInstant, CurrentHDF%NextField2D, CurrentHDF)
            
            !limit maximum values
            do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
            
#ifndef _NOT_IEEE_ARITHMETIC
                if (ieee_is_nan (CurrentHDF%PreviousField2D(i,j)))                      &
                    CurrentHDF%PreviousField2D(i,j) = FillValueReal                
#endif

                if (abs(CurrentHDF%PreviousField2D(i,j)) > abs(FillValueReal))          &
                    CurrentHDF%PreviousField2D(i,j) = FillValueReal
            enddo
            enddo   
                
        end if i11

i15:    if (CurrentHDF%PreviousInstant /= CurrentHDF%NextInstant) then

i16:        if (Me%UseOriginalValues) then
                
                Me%Matrix2D = CurrentHDF%NextField2D
                
i17:            if (Me%PredictDTMethod == 1) then
                    call PredictDTForHDF (PointsToFill2D, CurrentHDF%PreviousTime, CurrentHDF%NextTime, Now, CurrentHDF)
                elseif (Me%PredictDTMethod == 2) then i17
                    call PredictDTForHDF_New (PointsToFill2D, Now, CurrentHDF)
                else i17
                    stop 'ModifyHDFInput2DStandard - ModuleFillMatrix - ERR010'
                endif i17
            
            else if (Me%AccumulateValues) then  i16     !For Rain
                    
                Me%Matrix2D = CurrentHDF%NextField2D / (CurrentHDF%NextTime - CurrentHDF%PreviousTime)
    
                !This will replace the processed values if the value in NextField2D is a NODATA value
i19:            if (.not. Me%IgnoreNoDataPoint) then
                    do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
                    do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
i18:                            if (CurrentHDF%NextField2D (i,j) == Me%NoDataValue) then
                            Me%Matrix2D (i,j) = Me%NoDataValue
                        endif i18
                    enddo
                    enddo
                endif i19
    
i20:            if (Me%PredictDTMethod == 1) then
                    call PredictDTForHDF (PointsToFill2D, CurrentHDF%PreviousTime, CurrentHDF%NextTime, Now, CurrentHDF)
                elseif (Me%PredictDTMethod == 2) then i20
                    call PredictDTForHDF_New (PointsToFill2D, Now, CurrentHDF)
                else  i20
                    stop 'ModifyHDFInput2DStandard - ModuleFillMatrix - ERR020'
                endif i20
                
            else i16
                !Interpolates the two matrixes in time
                if (Me%PropertyID%IsAngle) then
                    call InterpolateAngle2DInTime (ActualTime       = Now,                         &
                                                   Size             = Me%WorkSize2D,               &
                                                   Time1            = CurrentHDF%PreviousTime,     &
                                                   Matrix1          = CurrentHDF%PreviousField2D,  &
                                                   Time2            = CurrentHDF%NextTime,         &
                                                   Matrix2          = CurrentHDF%NextField2D,      &
                                                   MatrixOut        = Me%Matrix2D,                 &
                                                   PointsToFill2D   = PointsToFill2D)                
                else
                    call InterpolateMatrix2DInTime(ActualTime       = Now,                         &
                                                   Size             = Me%WorkSize2D,               &
                                                   Time1            = CurrentHDF%PreviousTime,     &
                                                   Matrix1          = CurrentHDF%PreviousField2D,  &
                                                   Time2            = CurrentHDF%NextTime,         &
                                                   Matrix2          = CurrentHDF%NextField2D,      &
                                                   MatrixOut        = Me%Matrix2D,                 &
                                                   PointsToFill2D   = PointsToFill2D)
                endif                                                   
                                               
                !This will replace the processed values if the value in NextField2D is a NODATA value
i21:            if (.not. Me%IgnoreNoDataPoint) then
                    do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
                    do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
i22:                            if (CurrentHDF%NextField2D (i,j) == Me%NoDataValue) then
                            Me%Matrix2D (i,j) = Me%NoDataValue
                        endif i22
                    enddo
                    enddo
                endif i21                                               
                
            endif i16
                
        else i15

            !Prev and next are equal (last instant?)
i23:        if (Me%UseOriginalValues .or. Me%InterpolateValues) then

                call SetMatrixValue(Me%Matrix2D, Me%WorkSize2D, CurrentHDF%PreviousField2D, PointsToFill2D)

            else i23
            
                !do nothing
                
            endif i23

        endif i15


    end subroutine ModifyHDFInput2DStandard
        
    !----------------------------------------------------------------------------    
    
    subroutine ModifyHDFInput2DTime(PointsToFill2D, CurrentHDF)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        type(T_Field4D)                                 :: CurrentHDF
        !Local----------------------------------------------------------------


        !Begin----------------------------------------------------------------
        
i1:     if (CurrentHDF%Field4D .and. CurrentHDF%HarmonicsOn) then

            call ModifyHDFInput2DHarmonics(PointsToFill2D, CurrentHDF)
            
        else i1 
        
i2:         if (Me%PredictDTMethod == 2) then

                call ModifyHDFInput2DRainType(PointsToFill2D, CurrentHDF)            
                      
            else i2
            
                call ModifyHDFInput2DStandard(PointsToFill2D, CurrentHDF)
                
            endif i2
                
        endif i1            

    end subroutine ModifyHDFInput2DTime

    !--------------------------------------------------------------------------
    !-------------------------------------------------------------------------
    subroutine BacktrackingTime(Now)

        !Arguments------------------------------------------------------------
        type (T_Time), intent(OUT)                      :: Now
        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        type (T_Time)                                   :: CurrentTime
        real                                            :: TotalTime, AuxPeriod

        !Begin----------------------------------------------------------------

        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'BacktrackingTime - ModuleFillMatrix - ERR10'

        TotalTime = Me%EndTime  - Me%BeginTime                  
        AuxPeriod = CurrentTime - Me%BeginTime
        AuxPeriod = TotalTime   - AuxPeriod
        
        Now = Me%BeginTime + AuxPeriod

    end subroutine BacktrackingTime
    
    !-------------------------------------------------------------------------
    
    logical function ReadNewField(Now,n, CurrentHDF)

        !Arguments------------------------------------------------------------
        type (T_Time), intent(OUT)                      :: Now
        integer      , intent(OUT)                      :: n
        type(T_Field4D)                                 :: CurrentHDF
        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        type (T_Time)                                   :: CurrentTime
        logical                                         :: ReadNewField_

        !Begin----------------------------------------------------------------

        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadNewField - ModuleFillMatrix - ERR10'

        !Backtracking time inversion is also done in the ModuleField4D    
        if (Me%BackTracking .and. .not. CurrentHDF%Field4D) then  
            call BacktrackingTime(Now)
        else   
            Now = CurrentTime
        endif     

        ReadNewField_ = .false.
        
        if (.not. Me%AccumulateValues) then
        !Backtracking time inversion is also done in the ModuleField4D    
        if (Me%BackTracking .and. .not. CurrentHDF%Field4D) then  
                if (Now .le. CurrentHDF%NextTime) ReadNewField_ = .true.
            else            
                if (Now .ge. CurrentHDF%NextTime) ReadNewField_ = .true.
            endif
        else
            !Backtracking time inversion is also done in the ModuleField4D    
            if (Me%BackTracking .and. .not. CurrentHDF%Field4D) then  
                if (Now .lt. CurrentHDF%NextTime) ReadNewField_ = .true.
            else            
                if (Now .gt. CurrentHDF%NextTime) ReadNewField_ = .true.
            endif        
        endif
        
        ReadNewField = ReadNewField_

        if (ReadNewField_)then

            n = 0
            
            do 

                !Backtracking time inversion is also done in the ModuleField4D    
                if (Me%BackTracking .and. .not. CurrentHDF%Field4D) then  
                    if (Now .gt. CurrentHDF%NextTime) exit
                else            
                    if (Now .lt. CurrentHDF%NextTime) exit
                endif            
                
                CurrentHDF%PreviousInstant  = CurrentHDF%NextInstant
                    
                !Backtracking time inversion is also done in the ModuleField4D    
                if (Me%BackTracking .and. .not. CurrentHDF%Field4D) then  
                    if(CurrentHDF%NextInstant .gt. 1)then
                        CurrentHDF%NextInstant  = CurrentHDF%NextInstant - 1
                    else
                        exit
                    endif
                else
                    if(CurrentHDF%NextInstant .lt. CurrentHDF%NumberOfInstants)then
                        CurrentHDF%NextInstant  = CurrentHDF%NextInstant + 1
                    else
                        !if (CurrentHDF%GenericYear) then 
                        if (CurrentHDF%CyclicTimeON) then                         
                            CurrentHDF%NextInstant  = 1
                        else
                            exit
                        endif
                    endif
                endif
                
                CurrentHDF%PreviousTime     = CurrentHDF%NextTime
                CurrentHDF%NextTime         = HDF5TimeInstant(CurrentHDF%NextInstant, CurrentHDF)
                
                !if (CurrentHDF%GenericYear) then
                if (CurrentHDF%CyclicTimeON) then
                    if(CurrentHDF%NextInstant > 1)then
                        !call SetHDFGenericYear(CurrentHDF%NextTime, Now)
                        call SetHDFCyclicDates(TimeInstant  = CurrentHDF%NextTime,          &
                                               RefTime      = Now,                          &
                                               CurrentHDF   = CurrentHDF)
                    else
                        !call SetHDFGenericYear(CurrentHDF%NextTime, Now, AddYear = .true.)
                        call SetHDFCyclicDates(TimeInstant  = CurrentHDF%NextTime,          &
                                               RefTime      = Now,                          & 
                                               RestartCycle = .true.,                       &
                                               CurrentHDF   = CurrentHDF)
                        
                    endif
                endif

                n = n + 1
                
                
            enddo
            
            !Backtracking time inversion is also done in the ModuleField4D    
            if (Me%BackTracking .and. .not. CurrentHDF%Field4D) then  
                if(Now .lt. CurrentHDF%NextTime)then
                    write(*,*)
                    write(*,*)'----------Backtracking mode-----------'
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'Time instants inconsistency.'
                    stop      'ReadNewField - ModuleFillMatrix - ERR20'
                end if
            else
                if(Now .gt. CurrentHDF%NextTime)then
                    write(*,*)
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'Time instants inconsistency.'
                    stop      'ReadNewField - ModuleFillMatrix - ERR30'
                end if
            endif    
    
        endif 
        
    end function ReadNewField

    !-------------------------------------------------------------------------

    subroutine ModifyHDFInput2DGeneric4D(PointsToFill2D, Generic_4D_Value_, CurrentHDF)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        real                                            :: Generic_4D_Value_
        type(T_Field4D)                                 :: CurrentHDF
        !Local----------------------------------------------------------------
        integer                                         :: PrevI, NextI, i, j 
        !Begin----------------------------------------------------------------
        
i1:     if (.not.(CurrentHDF%Previous4DValue <= Generic_4D_Value_ .and.                     &
                  CurrentHDF%Next4DValue     >= Generic_4D_Value_)) then
            !Found new limits
            PrevI              = 1
            NextI              = 2
            CurrentHDF%Next4DValue = HDF5Generic4DInstant(1, CurrentHDF)
            do 

                CurrentHDF%Previous4DValue  = CurrentHDF%Next4DValue
                CurrentHDF%Next4DValue      = HDF5Generic4DInstant(NextI, CurrentHDF)

                if (CurrentHDF%Previous4DValue <= Generic_4D_Value_ .and.                  &
                    CurrentHDF%Next4DValue     >= Generic_4D_Value_) then
                    exit
                endif

                if (NextI > CurrentHDF%NumberOfInstants) then

                    write(*,*)
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'Time instants inconsistency.'
                    stop      'ModifyHDFInput2DGeneric4D - ModuleFillMatrix - ERR10'

                endif

                PrevI = NextI
                NextI = NextI + 1

            enddo

            CurrentHDF%NextInstant     = NextI
            CurrentHDF%PreviousInstant = PrevI  


            !call SetMatrixValue(CurrentHDF%PreviousField2D, Me%WorkSize2D, CurrentHDF%NextField2D)

            call ReadHDF5Values2D(CurrentHDF%PreviousInstant, CurrentHDF%PreviousField2D, CurrentHDF)
            call ReadHDF5Values2D(CurrentHDF%NextInstant,     CurrentHDF%NextField2D, CurrentHDF)
            
            !limit maximum values
            do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB

#ifndef _NOT_IEEE_ARITHMETIC
                if (ieee_is_nan (CurrentHDF%PreviousField2D(i,j)))                          &
                    CurrentHDF%PreviousField2D(i,j) = FillValueReal                    
#endif
            
                if (abs(CurrentHDF%PreviousField2D(i,j)) > abs(FillValueReal))              &
                    CurrentHDF%PreviousField2D(i,j) = FillValueReal

#ifndef _NOT_IEEE_ARITHMETIC
                if (ieee_is_nan (CurrentHDF%NextField2D(i,j)))                              &
                    CurrentHDF%NextField2D    (i,j) = FillValueReal                    
#endif
                
                if (abs(CurrentHDF%NextField2D    (i,j)) > abs(FillValueReal))              &
                    CurrentHDF%NextField2D    (i,j) = FillValueReal
            enddo
            enddo   
            
                    
        endif i1
        

        if (CurrentHDF%PreviousInstant /= CurrentHDF%NextInstant) then

            call InterpolateLinearyMatrix2D(X                = Generic_4D_Value_,       &
                                            Size             = Me%WorkSize2D,           &
                                            X1               = CurrentHDF%Previous4DValue,  &
                                            Matrix1          = CurrentHDF%PreviousField2D,  &
                                            X2               = CurrentHDF%Next4DValue,      &
                                            Matrix2          = CurrentHDF%NextField2D,      &
                                            MatrixOut        = Me%Matrix2D,             &
                                            PointsToFill2D   = PointsToFill2D)

        else
          
            !Prev and next are equal (last instant?)
            call SetMatrixValue(Me%Matrix2D, Me%WorkSize2D, CurrentHDF%NextField2D)

        endif


    end subroutine ModifyHDFInput2DGeneric4D


    !--------------------------------------------------------------------------

    subroutine PredictDTForHDF(PointsToFill2D, Time1, Time2, ActualTime, CurrentHDF)
        
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        type (T_Time)                                   :: Time1, Time2, ActualTime
        type(T_Field4D)                                 :: CurrentHDF
        !Local-----------------------------------------------------------------
        integer                                         :: i, j
        real                                            :: aux1, aux2
        logical                                         :: ValueDifferentZero

        !Searches Maximum Value
        ValueDifferentZero = .false.
doj:    do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
        do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
            if(PointsToFill2D(i,j) == 1)then
                if (Me%Matrix2D(i, j) > Me%MinForDTDecrease) then
                    ValueDifferentZero = .true.
                    exit doj
                endif
            endif
        enddo
        enddo doj

        if (ValueDifferentZero) then

            if (Time2 /= ActualTime) then
                aux1 = Time2 - ActualTime
            else
                aux1 = -FillValueReal
            endif

            !to ensure that DT does not go through two intervals
            aux2 = Time2 - Time1

            CurrentHDF%PredictedDT     = min(aux1, aux2)
            CurrentHDF%DTForNextEvent  = 0.0
        
        else

            !Can run until next Matrix will be read
            !This prediction is different from the one done by the GetTimeSerieDTForNextEvent
            !but I (Frank) assume that the datasets in the HDF always have a quite larger DT
            !then the model will use to run
            CurrentHDF%PredictedDT     = Time2 - ActualTime
            CurrentHDF%DTForNextEvent  = Time2 - ActualTime
        endif
        
        !Compare to global value to update
        if (.not. Me%VectorialProp) then
            if (Me%DTForNextEvent > CurrentHDF%DTForNextEvent) &
                Me%DTForNextEvent = CurrentHDF%DTForNextEvent
                    
            if (Me%PredictedDT > CurrentHDF%PredictedDT) &
                Me%PredictedDT = CurrentHDF%PredictedDT 

            if (Me%DTForNextDataset > CurrentHDF%DTForNextDataset) &
                Me%DTForNextDataset = CurrentHDF%DTForNextDataset

            if ((CurrentHDF%DTForNextEvent == 0.0) .and. &
                (Me%NextValueForDTPred < CurrentHDF%NextValueForDTPred)) then
                Me%NextValueForDTPred = CurrentHDF%NextValueForDTPred
            endif   
        endif                        

    end subroutine PredictDTForHDF

    !--------------------------------------------------------------------------

    subroutine PredictDTForHDF_New(PointsToFill2D, ActualTime, CurrentHDF)
        
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        type (T_Time)                                   :: ActualTime
        type(T_Field4D)                                 :: CurrentHDF
        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------
        
        if (ActualTime >= CurrentHDF%NextTime) then
            call FindDTForNextEventInHDF(PointsToFill2D, ActualTime, CurrentHDF)             
        else
            if (CurrentHDF%DTForNextEvent > 0.0) then
                CurrentHDF%DTForNextEvent = CurrentHDF%PreviousTime - ActualTime
                if (CurrentHDF%DTForNextEvent<0.0) CurrentHDF%DTForNextEvent = 0.0
            endif
            CurrentHDF%DTForNextDataset = CurrentHDF%TimeOfNextDataset - ActualTime
            if (CurrentHDF%DTForNextDataset <= 0.0) then
                CurrentHDF%TimeOfNextDataset = TimeOfNextDatasetInHDF(ActualTime, CurrentHDF) 
                CurrentHDF%DTForNextDataset  = CurrentHDF%TimeOfNextDataset - ActualTime
            endif
        endif
        
        !Compare to global value to update
        if (.not. Me%VectorialProp) then
            if (Me%DTForNextEvent > CurrentHDF%DTForNextEvent) &
                Me%DTForNextEvent = CurrentHDF%DTForNextEvent
                    
            if (Me%PredictedDT > CurrentHDF%PredictedDT) &
                Me%PredictedDT = CurrentHDF%PredictedDT 

            if (Me%DTForNextDataset > CurrentHDF%DTForNextDataset) &
                Me%DTForNextDataset = CurrentHDF%DTForNextDataset

            if ((CurrentHDF%DTForNextEvent == 0.0) .and. &
                (Me%NextValueForDTPred < CurrentHDF%NextValueForDTPred)) then
                Me%NextValueForDTPred = CurrentHDF%NextValueForDTPred
            endif   
        endif         
        
        !----------------------------------------------------------------------

    end subroutine PredictDTForHDF_New

    !--------------------------------------------------------------------------
    
    type (T_Time) function TimeOfNextDatasetInHDF (ActualTime, CurrentHDF)
        !NOT for use with backtracking mode
    
        !Arguments-------------------------------------------------------------        
        Type (T_Time), intent(IN) :: ActualTime
        type(T_Field4D)           :: CurrentHDF
        !Local-----------------------------------------------------------------
         
        !----------------------------------------------------------------------
               
        TimeOfNextDatasetInHDF = CurrentHDF%TimeOfNextDataset
                    
!        if (CurrentHDF%CyclicTimeON) then 
!            TimeOfNextDatasetInHDF = CurrentHDF%NextTime        
!        else  
do1:     do while (TimeOfNextDatasetInHDF <= ActualTime)              
            if (CurrentHDF%InstantOfNextDataset < CurrentHDF%NumberOfInstants) then
                CurrentHDF%InstantOfNextDataset = CurrentHDF%InstantOfNextDataset + 1
                TimeOfNextDatasetInHDF = HDF5TimeInstant(CurrentHDF%InstantOfNextDataset, CurrentHDF)
            else
                exit do1
            endif
        enddo do1
!       endif
            
        !----------------------------------------------------------------------
    
    end function TimeOfNextDatasetInHDF
        
    !--------------------------------------------------------------------------
    
    subroutine FindDTForNextEventInHDF  (PointsToFill2D, ActualTime, CurrentHDF)
    
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer, intent(IN)   :: PointsToFill2D
        Type (T_Time), intent(IN)                       :: ActualTime
        type(T_Field4D)                                 :: CurrentHDF
        !Local-----------------------------------------------------------------
        integer                                         :: i, j        
        logical                                         :: DTForNextDatasetWasSet
        type (T_Time)                                   :: PreviousTime
        
        !----------------------------------------------------------------------
        PreviousTime = CurrentHDF%NextTime
        DTForNextDatasetWasSet = .false.

doF:    do

            CurrentHDF%PreviousInstant  = CurrentHDF%NextInstant
            
            if(CurrentHDF%NextInstant .lt. CurrentHDF%NumberOfInstants)then
                CurrentHDF%NextInstant  = CurrentHDF%NextInstant + 1                
                CurrentHDF%PreviousTime = CurrentHDF%NextTime
                CurrentHDF%NextTime     = HDF5TimeInstant(CurrentHDF%NextInstant, CurrentHDF)
                
                if (.not. DTForNextDatasetWasSet) then
                    DTForNextDatasetWasSet  = .true.
                    CurrentHDF%DTForNextDataset     = CurrentHDF%NextTime - ActualTime
                    CurrentHDF%TimeOfNextDataset    = CurrentHDF%NextTime
                    CurrentHDF%InstantOfNextDataset = CurrentHDF%NextInstant
                endif
            else
                exit doF
            endif

            call ReadHDF5Values2D(CurrentHDF%PreviousInstant, CurrentHDF%PreviousField2D, CurrentHDF)
            do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
#ifndef _NOT_IEEE_ARITHMETIC                    
                if (ieee_is_nan (CurrentHDF%PreviousField2D(i,j))) &
                    CurrentHDF%PreviousField2D(i,j) = FillValueReal
#endif                    
                if (abs(CurrentHDF%PreviousField2D(i,j)) > abs(FillValueReal)) &
                    CurrentHDF%PreviousField2D(i,j) = FillValueReal            
            enddo
            enddo   

            call ReadHDF5Values2D(CurrentHDF%NextInstant, CurrentHDF%NextField2D, CurrentHDF)
            do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
#ifndef _NOT_IEEE_ARITHMETIC            
                if (ieee_is_nan (CurrentHDF%NextField2D(i,j))) &
                    CurrentHDF%NextField2D (i,j) = FillValueReal                      
#endif                    
                if (abs(CurrentHDF%NextField2D(i,j)) > abs(FillValueReal)) &
                    CurrentHDF%NextField2D (i,j) = FillValueReal
            enddo
            enddo  

doM:        do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                if(PointsToFill2D(i,j) == 1)then
                    if (CurrentHDF%NextField2D(i,j) > Me%MinForDTDecrease) then
                        exit doF
                    endif
                endif
            enddo
            enddo doM
            
        enddo doF
        
        if (PreviousTime == CurrentHDF%PreviousTime) then
            CurrentHDF%DTForNextEvent = 0.0
        else
            CurrentHDF%DTForNextEvent = CurrentHDF%PreviousTime - ActualTime
        endif
        !----------------------------------------------------------------------
    
    end subroutine FindDTForNextEventInHDF
    
    !--------------------------------------------------------------------------

    subroutine ModifyProfileTimeSerie(PointsToFill3D)    
    
        !Arguments------------------------------------------------------------
        integer, dimension(:, :, :), pointer            :: PointsToFill3D

        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        type (T_Time)                                   :: Now

        !Begin----------------------------------------------------------------

        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyProfileTimeSerie - ModuleFillMatrix - ERR010'

        if(Now .ge. Me%ProfileTimeSerie%NextTime)then
            
            Me%ProfileTimeSerie%PreviousInstant  = Me%ProfileTimeSerie%NextInstant
            
            if(Me%ProfileTimeSerie%NextInstant .lt. Me%ProfileTimeSerie%NumberOfInstants)then
                Me%ProfileTimeSerie%NextInstant  = Me%ProfileTimeSerie%NextInstant + 1
            else
                if (Me%ProfileTimeSerie%CyclicTimeON) Me%ProfileTimeSerie%NextInstant  = 1
            end if
            
            Me%ProfileTimeSerie%PreviousTime     = Me%ProfileTimeSerie%NextTime
            Me%ProfileTimeSerie%NextTime         = Me%ProfileTimeSerie%TimeInstants(Me%ProfileTimeSerie%NextInstant)

            if (Me%ProfileTimeSerie%CyclicTimeON) then
                call CheckCyclicMonths(Me%ProfileTimeSerie%NextTime, PreviousTime = Me%ProfileTimeSerie%PreviousTime)
            endif

            if(Now .gt. Me%ProfileTimeSerie%NextTime)then
                write(*,*)
                write(*,*)'Could not read solution from profile time file'
                write(*,*)'Time instants inconsistency.'
                stop      'ModifyProfileTimeSerie - ModuleFillMatrix - ERR020'
            end if
            
            call SetMatrixValue(Me%ProfileTimeSerie%PreviousField3D,            &
                                Me%WorkSize3D, Me%ProfileTimeSerie%NextField3D, &
                                PointsToFill3D)

            call ProfileTimeSerieField(PointsToFill3D, Me%ProfileTimeSerie%NextInstant, &
                                                       Me%ProfileTimeSerie%NextField3D)

        end if

        if (Me%PropertyID%IsAngle) then            

            call InterpolateAngle3DInTime (ActualTime       = Now,                                   &
                                           Size             = Me%WorkSize3D,                         &
                                           Time1            = Me%ProfileTimeSerie%PreviousTime,      &
                                           Matrix1          = Me%ProfileTimeSerie%PreviousField3D,   &
                                           Time2            = Me%ProfileTimeSerie%NextTime,          &
                                           Matrix2          = Me%ProfileTimeSerie%NextField3D,       &
                                           MatrixOut        = Me%Matrix3D,                           &
                                           PointsToFill3D   = PointsToFill3D)

        else

            call InterpolateMatrix3DInTime(ActualTime       = Now,                                   &
                                           Size             = Me%WorkSize3D,                         &
                                           Time1            = Me%ProfileTimeSerie%PreviousTime,      &
                                           Matrix1          = Me%ProfileTimeSerie%PreviousField3D,   &
                                           Time2            = Me%ProfileTimeSerie%NextTime,          &
                                           Matrix2          = Me%ProfileTimeSerie%NextField3D,       &
                                           MatrixOut        = Me%Matrix3D,                           &
                                           PointsToFill3D   = PointsToFill3D)
        endif            
    end subroutine ModifyProfileTimeSerie

    !--------------------------------------------------------------------------


    subroutine ModifySpaceTimeSerie (CurrentTimeSerie, PointsToFill2D, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        type(T_TimeSerie)                               :: CurrentTimeSerie        
        integer, dimension(:, :   ), pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D

        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        type (T_Time)                                   :: Now
        type (T_Time)                                   :: Time1, Time2
        real                                            :: Value1, Value2, NewValue
        logical                                         :: TimeCycle
        real                                            :: DT1, DT2, Angle
        real                                            :: u1, u2, v1, v2, uf, vf
        real                                            :: aux 
        
        !Begin----------------------------------------------------------------
        
        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR010'                  
            
        if (.not. CurrentTimeSerie%RemainsConstant) then        
            
            if (Me%PredictDTMethod == 1) then
                !Gets Value for current Time
                call GetTimeSerieValue (CurrentTimeSerie%ObjTimeSerie, Now,               &
                                        CurrentTimeSerie%Column,                          &
                                        Time1, Value1, Time2, Value2, TimeCycle,         &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR020'

        
                if (Me%PropertyID%IDNumber /= WindDirection_) then

                    if (TimeCycle) then
                    
                        NewValue = Value1

                    else

                        if (Me%UseOriginalValues) then
                
                            NewValue = Value2
                    
                        elseif (Me%AccumulateValues) then       !For Rain
                
                            NewValue = Value2 / (Time2 - Time1)
    
                            call GetTimeSerieDTForNextEvent (CurrentTimeSerie%ObjTimeSerie,          &
                                                                NewValue, CurrentTimeSerie%Column, Now, &
                                                                CurrentTimeSerie%PredictedDT, CurrentTimeSerie%DTForNextEvent,  &
                                                                STAT  = STAT_CALL)
                                    if (STAT_CALL /= SUCCESS_) stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR030'
                
                        else
                            !Interpolates Value for current instant
                            call InterpolateValueInTime(Now, Time1, Value1, Time2, Value2, NewValue)
                        endif
                    endif

                else
                    
                    !Gets Grid Angle
                    call GetGridAngle(Me%ObjHorizontalGrid, Angle, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR040'                        
                        
                    if (TimeCycle) then

                        NewValue = Value1  

                    else                        
                
                        DT1      = Now   - Time1
                        DT2      = Time2 - Now

                        !Interpolates angle
                        u1 = DT2 * cos(Value1/180*PI)
                        v1 = DT2 * sin(Value1/180*PI)

                        u2 = DT1 * cos(Value2/180*PI)
                        v2 = DT1 * sin(Value2/180*PI)

                        uf = u1 + u2
                        vf = v1 + v2

                        !1st Quad
                        if (uf > 0. .and. vf > 0.) then
                            NewValue = atan(vf / uf) / PI * 180.
                        !2 e 3 Quad
                        elseif (uf < 0.) then
                            NewValue = atan(vf / uf) / PI * 180. + 180.
                        !4 Quad
                        elseif (uf > 0. .and. vf < 0.) then
                            NewValue = atan(vf / uf) / PI * 180. + 360.
                        !one value zero
                        else
                            NewValue = (DT2 * Value1 + DT1 * Value2) / (DT1 + DT2)
                        endif

                    endif
                        
            
                endif        
            
            else !Me%PredictDTMethod == 2
            
                !This method do not works with TimeCycle
!                if (TimeCycle) then
!                    write(*,*) 'The method 2 to predict DT do not works with TimeCycle.'
!                    stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR050'
!                endif
                
                if (Now > CurrentTimeSerie%NextTime) then
                    call ActualizeTimeSerieTimes (Now, CurrentTimeSerie)
                    call ActualizeTimeSerieValues(CurrentTimeSerie)                   
                endif
                
                !avoid evaluate in cosntruct phase where previous and next time are the same
                if (Now > CurrentTimeSerie%PreviousTime) then
                    if (Me%UseOriginalValues) then
                        CurrentTimeSerie%CurrentValue = CurrentTimeSerie%NextValue
                    elseif (Me%AccumulateValues) then
                        CurrentTimeSerie%CurrentValue = CurrentTimeSerie%NextValue      &
                               / (CurrentTimeSerie%NextTime - CurrentTimeSerie%PreviousTime)
                    elseif (Me%InterpolateValues) then
                        call InterpolateValueInTime(Now,                        &
                                                    CurrentTimeSerie%PreviousTime,  &
                                                    CurrentTimeSerie%PreviousValue, &
                                                    CurrentTimeSerie%NextTime,      &
                                                    CurrentTimeSerie%NextValue,     &
                                                    CurrentTimeSerie%CurrentValue)                    
                    else
                        stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR060'
                    endif
                endif
                
                NewValue = CurrentTimeSerie%CurrentValue
                
                if (Me%ValueIsUsedForDTPrediction) then
                    if (Now >= CurrentTimeSerie%NextEventEnd) then
                        call FindNextEventInTimeSerie (Now, CurrentTimeSerie)
                        if (Me%AccumulateValues .and. (CurrentTimeSerie%NextValueForDTPred > 0.0)) then
                            aux = CurrentTimeSerie%NextEventEnd - CurrentTimeSerie%NextEventStart
                            if (aux > 0.) then
                                CurrentTimeSerie%NextValueForDTPred = CurrentTimeSerie%NextValueForDTPred / aux
                            endif                                
                        endif                    
                    endif
                    
                    if (Now >= CurrentTimeSerie%NextEventStart .and. Now < CurrentTimeSerie%NextEventEnd) then
                        CurrentTimeSerie%DTForNextEvent = 0.0
                    else
                        CurrentTimeSerie%DTForNextEvent = CurrentTimeSerie%NextEventStart - Now 
                    endif
                    
                    if (CurrentTimeSerie%DTForNextEvent > 0.0) then
                        CurrentTimeSerie%PredictedDT = CurrentTimeSerie%DTForNextEvent                    
                    else
                        CurrentTimeSerie%PredictedDT = CurrentTimeSerie%NextEventEnd - Now
                    endif                
                endif            
                
            endif
            
            !Compare to global value to update
            if (.not. Me%VectorialProp) then
                if (Me%DTForNextEvent > CurrentTimeSerie%DTForNextEvent) &
                    Me%DTForNextEvent = CurrentTimeSerie%DTForNextEvent
                    
                if (Me%PredictedDT > CurrentTimeSerie%PredictedDT) &
                    Me%PredictedDT = CurrentTimeSerie%PredictedDT 

                if (Me%DTForNextDataset > CurrentTimeSerie%DTForNextDataset) &
                    Me%DTForNextDataset = CurrentTimeSerie%DTForNextDataset

                if ((CurrentTimeSerie%DTForNextEvent == 0.0) .and. &
                    (Me%NextValueForDTPred < CurrentTimeSerie%NextValueForDTPred)) then
                    Me%NextValueForDTPred = CurrentTimeSerie%NextValueForDTPred
                endif   
            endif            
            
        else
        
            NewValue = CurrentTimeSerie%CurrentValue
            
        endif
          
       
        if (Me%Dim == Dim2D) then
            call SetMatrixValue(Me%Matrix2D, Me%WorkSize2D, NewValue, PointsToFill2D)
        else
            call SetMatrixValue(Me%Matrix3D, Me%WorkSize3D, NewValue, PointsToFill3D)
        endif
     
                        
    end subroutine ModifySpaceTimeSerie
    
    !--------------------------------------------------------------------------
    
!    subroutine PredictDTForTimeSerie(ActualTime)
!        
!        !Arguments-------------------------------------------------------------
!        type (T_Time)                                   :: ActualTime
!
!        !Local-----------------------------------------------------------------
!
!        !----------------------------------------------------------------------
!        
!        if (ActualTime >= Me%TimeSerie%NextTime) then
!        
!            call FindDTForNextEventInTimeSerie(Me%TimeSerie%ObjTimeSerie, ActualTime) 
!                        
!        else
!            if (Me%DTForNextEvent > 0.0) then
!            
!                Me%DTForNextEvent = max(0.0, Me%TimeSerie%PreviousTime - ActualTime)
!                
!            endif
!            
!            Me%DTForNextDataset = Me%TimeOfNextDataset - ActualTime
!            
!            if (Me%DTForNextDataset <= 0.0) then
!            
!                Me%TimeOfNextDataset = TimeOfNextDatasetInTimeSerie(Me%TimeSerie%ObjTimeSerie, ActualTime) 
!                Me%DTForNextDataset  = Me%TimeOfNextDataset - ActualTime
!                
!            endif
!            
!        endif
!        
!        !----------------------------------------------------------------------
!
!    end subroutine PredictDTForTimeSerie
!
!    !--------------------------------------------------------------------------
!    
!    type (T_Time) function TimeOfNextDatasetInTimeSerie (TS_ID, ActualTime)
!    
!        !Arguments-------------------------------------------------------------        
!        integer, intent(IN)         :: TS_ID
!        Type (T_Time), intent(IN)   :: ActualTime
!        
!        !Local-----------------------------------------------------------------
!        integer                     :: index
!        integer                     :: STAT_
!        
!        !----------------------------------------------------------------------
!               
!        TimeOfNextDatasetInTimeSerie = Me%TimeOfNextDataset
!                    
!do1:    do while (TimeOfNextDatasetInTimeSerie <= ActualTime)
!          
!            if (Me%InstantOfNextDataset < Me%TimeSerie%NumberOfInstants) then
!            
!                Me%InstantOfNextDataset = Me%InstantOfNextDataset + 1
!                
!                call GetTimeSerieTimeOfNextDataset(TS_ID, ActualTime, TimeOfNextDatasetInTimeSerie, STAT = STAT_)
!                if (STAT_ /= SUCCESS_) stop 'TimeOfNextDatasetInTimeSerie - ModuleFillMatrix - ERR010'
!                
!            else
!            
!                exit do1
!                
!            endif
!            
!        enddo do1
!            
!        !----------------------------------------------------------------------
!    
!    end function TimeOfNextDatasetInTimeSerie
        
    !--------------------------------------------------------------------------
    
    subroutine ActualizeTimeSerieTimes (ActualTime, CurrentTimeSerie)
    
        !Arguments-------------------------------------------------------------        
        Type (T_Time), intent(IN)                       :: ActualTime
        Type (T_TimeSerie)                              :: CurrentTimeSerie
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        
        !----------------------------------------------------------------------                     
        

        do while (CurrentTimeSerie%NextTime < ActualTime)  
                          
            if (CurrentTimeSerie%NextInstant < CurrentTimeSerie%NumberOfInstants) then
            
                CurrentTimeSerie%PreviousInstant = CurrentTimeSerie%NextInstant
                CurrentTimeSerie%PreviousTime    = CurrentTimeSerie%NextTime
                CurrentTimeSerie%NextInstant     = CurrentTimeSerie%NextInstant + 1
                
            else
            
                stop 'ActualizeTimeSerieTimes - ModuleFillMatrix - ERR010'
                
            endif
                  
            call GetTimeSerieTimeOfDataset(CurrentTimeSerie%ObjTimeSerie,    &
                                            CurrentTimeSerie%NextInstant,     &
                                            CurrentTimeSerie%NextTime,        &
                                            STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)    &
                stop 'ActualizeTimeSerieTimes - ModuleFillMatrix - ERR020'  
           
        enddo
            

        
        !----------------------------------------------------------------------
        
    end subroutine ActualizeTimeSerieTimes 
    
    !--------------------------------------------------------------------------
    
    subroutine ActualizeTimeSerieValues(CurrentTimeSerie)
    
        !Arguments-------------------------------------------------------------
        Type (T_TimeSerie)                              :: CurrentTimeSerie
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_
        
        !----------------------------------------------------------------------    

        
        call GetTimeSerieValueForIndex (CurrentTimeSerie%ObjTimeSerie,      &
                                        CurrentTimeSerie%PreviousInstant,   &
                                        CurrentTimeSerie%Column,            &
                                        CurrentTimeSerie%PreviousValue,     &
                                        STAT = STAT_)
        if (STAT_ /= SUCCESS_) stop 'ActualizeTimeSerieValues - ModuleFillMatrix - ERR010'            

        call GetTimeSerieValueForIndex (CurrentTimeSerie%ObjTimeSerie,      &
                                        CurrentTimeSerie%NextInstant,       &
                                        CurrentTimeSerie%Column,            &
                                        CurrentTimeSerie%NextValue,         &                                            
                                        STAT = STAT_)
        if (STAT_ /= SUCCESS_) stop 'ActualizeTimeSerieValues - ModuleFillMatrix - ERR020' 
            

            
            
        !----------------------------------------------------------------------    
    
    end subroutine ActualizeTimeSerieValues
    
    !--------------------------------------------------------------------------
    
    subroutine FindNextEventInTimeSerie(Now, CurrentTimeSerie)
    
        !Arguments-------------------------------------------------------------        
        Type (T_Time), intent(IN)                       :: Now
        type(T_TimeSerie)                               :: CurrentTimeSerie
        !Local-----------------------------------------------------------------
        real                                            :: instant_value
        integer                                         :: instant, STAT_CALL

        !----------------------------------------------------------------------
     
        
        if (CurrentTimeSerie%NextInstant < CurrentTimeSerie%NumberOfInstants) then
        
            if (Now > CurrentTimeSerie%NextEventEnd) then
        
                CurrentTimeSerie%NextEventStart = CurrentTimeSerie%PreviousTime
                instant                         = CurrentTimeSerie%NextInstant
                CurrentTimeSerie%NextEventEnd   = CurrentTimeSerie%NextTime                
                instant_value                   = CurrentTimeSerie%NextValue

            else
            
                CurrentTimeSerie%NextEventStart = CurrentTimeSerie%NextTime
                instant           = CurrentTimeSerie%NextInstant + 1
                call GetTimeSerieTimeOfDataset(CurrentTimeSerie%ObjTimeSerie,   &
                                                instant,                        &
                                                CurrentTimeSerie%NextEventEnd,  &
                                                STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'FindNextEventInTimeSerie - ModuleFillMatrix - ERR010'                 
                call GetTimeSerieValueForIndex (CurrentTimeSerie%ObjTimeSerie,  &
                                                instant,                    &
                                                CurrentTimeSerie%Column,        &
                                                instant_value,              &                                            
                                                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'FindNextEventInTimeSerie - ModuleFillMatrix - ERR020'   
            
            endif
                             
            do while (instant_value <= Me%MinForDTDecrease .and. instant < CurrentTimeSerie%NumberOfInstants)                
                instant = instant + 1  
                    
                call GetTimeSerieValueForIndex (CurrentTimeSerie%ObjTimeSerie,  &
                                                instant,                    &
                                                CurrentTimeSerie%Column,        &
                                                instant_value,              &                                            
                                                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'FindNextEventInTimeSerie - ModuleFillMatrix - ERR030'                                
            enddo

            if (instant_value > Me%MinForDTDecrease) then
            
                call GetTimeSerieTimeOfDataset(CurrentTimeSerie%ObjTimeSerie,   &
                                                Instant - 1,                    &
                                                CurrentTimeSerie%NextEventStart, &
                                                STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'FindNextEventInTimeSerie - ModuleFillMatrix - ERR040' 
                
                call GetTimeSerieTimeOfDataset(CurrentTimeSerie%ObjTimeSerie,   &
                                                instant,                        &
                                                CurrentTimeSerie%NextEventEnd,  &
                                                STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'FindNextEventInTimeSerie - ModuleFillMatrix - ERR050'              
                    
            else
            
                call GetTimeSerieTimeOfDataset(CurrentTimeSerie%ObjTimeSerie,   &
                                                instant,                         &
                                                CurrentTimeSerie%NextEventStart, &
                                                STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'FindNextEventInTimeSerie - ModuleFillMatrix - ERR060'                                 
                        
                CurrentTimeSerie%NextEventEnd = CurrentTimeSerie%NextEventStart
                instant_value   = 0.0
                                    
            endif 
            
            CurrentTimeSerie%NextValueForDTPred = instant_value
            
        else
        
            CurrentTimeSerie%NextEventStart = CurrentTimeSerie%NextTime
            CurrentTimeSerie%NextEventEnd   = CurrentTimeSerie%NextTime                

        endif
        
        
        !write (*,*) 'instant_value = ', instant_value
     
        
        !----------------------------------------------------------------------
    
    end subroutine FindNextEventInTimeSerie
    
    !--------------------------------------------------------------------------
    
    subroutine ActualizeHDFTimes (ActualTime, CurrentHDF)
    
        !Arguments-------------------------------------------------------------        
        Type (T_Time), intent(IN)                       :: ActualTime
        type(T_Field4D)                                 :: CurrentHDF
        !Local-----------------------------------------------------------------
        
        !----------------------------------------------------------------------                     
        
        do while (CurrentHDF%NextTime < ActualTime)  
                          
            if (CurrentHDF%NextInstant < CurrentHDF%NumberOfInstants) then
            
                CurrentHDF%PreviousInstant = CurrentHDF%NextInstant
                CurrentHDF%PreviousTime    = CurrentHDF%NextTime
                CurrentHDF%NextInstant     = CurrentHDF%NextInstant + 1
                
            else
            
                stop 'ActualizeHDFTimes - ModuleFillMatrix - ERR010'
                
            endif
                  
            CurrentHDF%NextTime = HDF5TimeInstant(CurrentHDF%NextInstant, CurrentHDF)
           
        enddo
 
        !----------------------------------------------------------------------
        
    end subroutine ActualizeHDFTimes 
    
    !--------------------------------------------------------------------------    

    subroutine ActualizeHDFValues (CurrentHDF, instant, field2D, field3D)
    
        !Arguments-------------------------------------------------------------
        type(T_Field4D)                             :: CurrentHDF
        integer                                     :: instant
        real, dimension(:, :), pointer, optional    :: field2D
        real, dimension(:, :, :), pointer, optional :: field3D

        !Local-----------------------------------------------------------------
        integer :: i, j, k
        
        !----------------------------------------------------------------------    
    
        if (present(field2D)) then
    
            call ReadHDF5Values2D(instant, field2D, CurrentHDF)        
        
            !limit maximum values
            do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
        
#ifndef _NOT_IEEE_ARITHMETIC                    
                if (ieee_is_nan (field2D(i,j))) &
                    field2D(i,j) = FillValueReal                
#endif
            
                if (abs(field2D(i,j)) > abs(FillValueReal)) &
                    field2D(i,j) = FillValueReal

            enddo
            enddo   
        
        elseif (present(field3D)) then
        
            call ReadHDF5Values3D(instant, field3D, currentHDF)        
        
            !limit maximum values
            do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
        
#ifndef _NOT_IEEE_ARITHMETIC                    
                if (ieee_is_nan (field3D(i,j,k))) &
                    field3D(i,j,k) = FillValueReal                
#endif
            
                if (abs(field3D(i,j,k)) > abs(FillValueReal)) &
                    field3D(i,j,k) = FillValueReal

            enddo
            enddo 
            enddo
                    
        endif
        
        !----------------------------------------------------------------------    
    
    end subroutine ActualizeHDFValues
    
    !--------------------------------------------------------------------------

    subroutine FindNextEventInHDF(Now, CurrentHDF)
    
        !Arguments-------------------------------------------------------------        
        Type (T_Time), intent(IN)                       :: Now
        type(T_Field4D)                                 :: CurrentHDF
        !Local-----------------------------------------------------------------
        real                                            :: instant_value
        integer                                         :: instant

        !----------------------------------------------------------------------
        
        if (CurrentHDF%NextInstant < CurrentHDF%NumberOfInstants) then

            if (Now > CurrentHDF%NextEventEnd) then
                CurrentHDF%NextEventStart = CurrentHDF%PreviousTime
                instant                   = CurrentHDF%NextInstant
                CurrentHDF%NextEventEnd   = CurrentHDF%NextTime   
                if (Me%Dim == Dim2D) then             
                    CurrentHDF%Array2D    = CurrentHDF%NextField2D
                else
                    CurrentHDF%Array3D    = CurrentHDF%NextField3D
                endif
            else
            
                CurrentHDF%NextEventStart = CurrentHDF%NextTime
                instant                   = CurrentHDF%NextInstant + 1
                CurrentHDF%NextEventEnd   = HDF5TimeInstant(instant, CurrentHDF)
                 
                if (Me%Dim == Dim2D) then
                    call ActualizeHDFValues (CurrentHDF, instant, field2D = CurrentHDF%Array2D)
                else
                    call ActualizeHDFValues (CurrentHDF, instant, field3D = CurrentHDF%Array3D)
                endif
            endif
                 
            if (Me%Dim == Dim2D) then     
                instant_value = maxval(CurrentHDF%Array2D)
            else
                instant_value = maxval(CurrentHDF%Array3D)
            endif
                             
            do while (instant_value <= Me%MinForDTDecrease .and. instant < CurrentHDF%NumberOfInstants)                
                instant = instant + 1  
                    
                if (Me%Dim == Dim2D) then
                    call ActualizeHDFValues (CurrentHDF, instant, field2D = CurrentHDF%Array2D)
                    instant_value = maxval(CurrentHDF%Array2D)
                else
                    call ActualizeHDFValues (CurrentHDF, instant, field3D = CurrentHDF%Array3D)
                    instant_value = maxval(CurrentHDF%Array3D)
                endif                               
            enddo

            if (instant_value > Me%MinForDTDecrease) then
            
                CurrentHDF%NextEventStart   = HDF5TimeInstant(instant - 1, CurrentHDF)
                CurrentHDF%NextEventEnd     = HDF5TimeInstant(instant, CurrentHDF)
                    
            else
            
                CurrentHDF%NextEventStart   = HDF5TimeInstant(instant, CurrentHDF)                        
                CurrentHDF%NextEventEnd     = CurrentHDF%NextEventStart
                instant_value       = 0.0
                                    
            endif 
            
        else
        
            CurrentHDF%NextEventStart = CurrentHDF%NextTime
            CurrentHDF%NextEventEnd   = CurrentHDF%NextTime                

        endif
        
        CurrentHDF%NextValueForDTPred = instant_value        

        !----------------------------------------------------------------------
    
    end subroutine FindNextEventInHDF
    
    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillFillMatrix(FillMatrixID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: FillMatrixID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers, STAT_CALL           
        type(T_Field4D), pointer            :: CurrentHDF
        type(T_Timeserie), pointer          :: CurrentTimeSerie
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mFILLMATRIX_,  Me%InstanceID)

            if (nUsers == 0) then
                                
                CurrentTimeSerie => Me%FirstTimeSerie
                do while (associated(CurrentTimeSerie))
                
                    !Kills Time Serie
                    if (CurrentTimeSerie%ObjTimeSerie /= 0) then
                        call KillTimeSerie (CurrentTimeSerie%ObjTimeSerie, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillFillMatrix - ModuleFillMatrix - ERR10'
                    endif
                    
                    CurrentTimeSerie => CurrentTimeSerie%Next
                
                enddo
                    
                CurrentHDF => Me%FirstHDF
                do while (associated(CurrentHDF))
                    
                    if (CurrentHDF%ObjHDF5 /= 0 .or. CurrentHDF%ObjField4D /= 0 ) then

                        if( associated(CurrentHDF%PreviousField2D))then
                            deallocate(CurrentHDF%PreviousField2D)
                            nullify   (CurrentHDF%PreviousField2D)
                        end if

                        if( associated(CurrentHDF%NextField2D))then
                            deallocate(CurrentHDF%NextField2D)
                            nullify   (CurrentHDF%NextField2D)
                        end if

                        if( associated(CurrentHDF%PreviousField3D))then
                            deallocate(CurrentHDF%PreviousField3D)
                            nullify   (CurrentHDF%PreviousField3D)
                        end if

                        if(associated (CurrentHDF%NextField3D))then
                            deallocate(CurrentHDF%NextField3D)
                            nullify   (CurrentHDF%NextField3D)
                        end if
                    
                        if (associated(CurrentHDF%ReadField3D)) then
                            deallocate(CurrentHDF%ReadField3D, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)stop 'KillFillMatrix - ModuleFillMatrix - ERR20'
                            nullify   (CurrentHDF%ReadField3D)
                        endif
                    
    if4D:               if (CurrentHDF%Field4D) then

                            if (CurrentHDF%SpatialInterpolON) then
    
                                if (Me%Dim == Dim3D .and. .not. CurrentHDF%From2Dto3D) then
                                     deallocate(CurrentHDF%Z)
                                endif                    
                                
                                deallocate(CurrentHDF%X, CurrentHDF%Y, CurrentHDF%Prop, CurrentHDF%NoData)
                        
                            endif
                            
                            if (associated(CurrentHDF%FileNameList)) then
                                deallocate(CurrentHDF%FileNameList)
                            endif
                        
                            call KillField4D(CurrentHDF%ObjField4D, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillFillMatrix - ModuleFillMatrix - ERR30'
                      
                        else if4D

                            call KillHDF5(CurrentHDF%ObjHDF5, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillFillMatrix - ModuleFillMatrix - ERR40'

                        endif if4D
                      
                        if (CurrentHDF%Generic4D%ReadFromTimeSerie) then
                            call KillTimeSerie(CurrentHDF%Generic4D%ObjTimeSerie, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillFillMatrix - ModuleFillMatrix - ERR50'
                        endif                        
                    
                    endif
  
                    CurrentHDF => CurrentHDF%Next
                    
                enddo
                
                if (Me%AnalyticWave%ON) then
                    call KillAnalyticWave
                endif                    

                if (Me%ObjGeometry /= 0) then
                    nUsers = DeassociateInstance (mGEOMETRY_, Me%ObjGeometry)
                    if (nUsers == 0) stop 'KillFillMatrix - ModuleFillMatrix - ERR60'
                endif

                nUsers = DeassociateInstance (mTIME_,           Me%ObjTime           )
                if (nUsers == 0) stop 'KillFillMatrix - ModuleFillMatrix - ERR70'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid )
                if (nUsers == 0) stop 'KillFillMatrix - ModuleFillMatrix - ERR80'
                
                !Deallocates Instance
                call DeallocateInstance ()

                FillMatrixID = 0
                STAT_           = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillFillMatrix
        

    !-----------------------------------------------------------------------------------
    
    subroutine KillAnalyticWave

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        
        deallocate(Me%AnalyticWave%X2D                  )
        deallocate(Me%AnalyticWave%CellType             )
        deallocate(Me%AnalyticWave%TlagMissing          )
        deallocate(Me%AnalyticWave%AmpAux               )
        deallocate(Me%AnalyticWave%Celerity             )
        deallocate(Me%AnalyticWave%EnteringCell%i       )
        deallocate(Me%AnalyticWave%EnteringCell%j       )
        deallocate(Me%AnalyticWave%EnteringCell%TimeLag )        
    
    end subroutine KillAnalyticWave
    
    !-----------------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_FillMatrix), pointer          :: AuxObjFillMatrix
        type (T_FillMatrix), pointer          :: PreviousObjFillMatrix

        !Updates pointers
        if (Me%InstanceID == FirstObjFillMatrix%InstanceID) then
            FirstObjFillMatrix => FirstObjFillMatrix%Next
        else
            PreviousObjFillMatrix => FirstObjFillMatrix
            AuxObjFillMatrix      => FirstObjFillMatrix%Next
            do while (AuxObjFillMatrix%InstanceID /= Me%InstanceID)
                PreviousObjFillMatrix => AuxObjFillMatrix
                AuxObjFillMatrix      => AuxObjFillMatrix%Next
            enddo

            !Now update linked list
            PreviousObjFillMatrix%Next => AuxObjFillMatrix%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjFillMatrix_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjFillMatrix_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjFillMatrix_ID > 0) then
            call LocateObjFillMatrix (ObjFillMatrix_ID)
            ready_ = VerifyReadLock (mFILLMATRIX_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjFillMatrix (ObjFillMatrixID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjFillMatrixID

        !Local-----------------------------------------------------------------

        Me => FirstObjFillMatrix
        do while (associated (Me))
            if (Me%InstanceID == ObjFillMatrixID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleFillMatrix - LocateObjFillMatrix - ERR01'

    end subroutine LocateObjFillMatrix

    !--------------------------------------------------------------------------

end module ModuleFillMatrix

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
