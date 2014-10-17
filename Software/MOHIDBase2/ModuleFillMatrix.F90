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
                                       InterpolateLinearyMatrix2D, InterpolateLinearyMatrix3D
    use ModuleBoxDif,           only : StartBoxDif, GetBoxes, GetNumberOfBoxes,         &
                                       UngetBoxDif, KillBoxDif
    use ModuleGridData,         only : ConstructGridData, GetGridData, UnGetGridData,   &
                                       KillGridData, GetGridDataType 
    use ModuleHorizontalGrid,   only : GetGridAngle, GetHorizontalGridSize,             &
                                       GetGridBorderLimits, GetLatitudeLongitude,       &
                                       GetDDecompOpenBorders,                           &
                                       GetDDecompParameters,                            &
                                       GetDDecompWorkSize2D, GetZCoordinates,           &
                                       UnGetHorizontalGrid
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
                                       GetHDF5GroupExist
                                       
    use ModuleField4D,          only : ConstructField4D, GetField4DNumberOfInstants,    &
                                       GetField4DInstant, ModifyField4D,                &
                                       ModifyField4DXYZ, KillField4D
                                       

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
                     
    
    !Modifier
    public  :: ModifyFillMatrix
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
    end interface ConstructFillMatrix


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
    integer, parameter                              :: sponge_exp_               = 1
    integer, parameter                              :: sponge_linear_            = 2
    
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
    

    !Types---------------------------------------------------------------------

    type T_Layers
        real, dimension(:), pointer                 :: Values   => null()
    end type T_Layers

    type T_Boxes
        character(PathLength)                       :: FileName     = null_str !initialization: Jauch
        integer                                     :: ObjBoxDif    = null_int !initialization: Jauch
        real, dimension(:), pointer                 :: Values       => null()
    end type T_Boxes

    type T_ASCIIFile
        character(PathLength)                       :: FileName     = null_str !initialization: Jauch
        integer                                     :: GridDataID   = null_int !initialization: Jauch
    end type T_ASCIIFile

    type T_Sponge
        real                                        :: OutValue     = null_real !initialization: Jauch
        integer                                     :: Cells        = null_int  !initialization: Jauch
        logical                                     :: Growing      = .false.   !initialization: Jauch
        integer                                     :: Evolution    = null_int  !initialization: Jauch
        !1 - South; 2 - North; 3 - West; 4 - East        
        logical, dimension(1:4)                     :: OpenBordersON = .true.
    end type T_Sponge

    type T_TimeSerie
        character(PathLength)                       :: FileName     = null_str !initialization: Jauch
        integer                                     :: ObjTimeSerie = 0
        integer                                     :: Column       = null_int !initialization: Jauch
        type (T_Time)                               :: NextTime, &
                                                       PreviousTime
        integer                                     :: NextInstant      = 0,     &
                                                       PreviousInstant  = 0
        real                                        :: PreviousValue    = 0.,    & 
                                                       NextValue        = 0.,    &
                                                       CurrentValue     = 0.
        integer                                     :: NumberOfInstants = 0
    end type T_TimeSerie

    type T_ProfileTimeSerie
        character(PathLength)                       :: FileName = null_str !initialization: Jauch
        type (T_Time)                               :: NextTime, PreviousTime
        integer                                     :: NextInstant      = null_int, & !initialization: Jauch
                                                       PreviousInstant  = null_int    !initialization: Jauch
        real,           dimension(:,:,:), pointer   :: PreviousField3D  => null(), &
                                                       NextField3D      => null()
        real,           dimension(:,:  ), pointer   :: Values           => null(), &
                                                       Depths           => null()
        type(T_Time),   dimension(:    ), pointer   :: TimeInstants     => null()
        integer                                     :: NumberOfInstants = null_int, & !initialization: Jauch
                                                       nValues          = null_int, & !initialization: Jauch
                                                       nDepths          = null_int, & !initialization: Jauch     
                                                       FirstInstant     = null_int, & !initialization: Jauch
                                                       LastInstant      = null_int    !initialization: Jauch
        logical                                     :: CyclicTimeON     = .false.
    end type T_ProfileTimeSerie

    type T_Station   
        character(PathLength)                       :: FileName         = null_str !initialization: Jauch
        integer                                     :: ObjTimeSerie     = 0        
        integer                                     :: Column           = null_int !initialization: Jauch
        integer                                     :: FillID           = null_int !initialization: Jauch
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
        integer                                     :: DataProcessing   = null_int !initialization: Jauch
        type(T_Station), dimension(:), allocatable  :: StationsList
        integer                                     :: NumberOfSources  = 0
    end type T_MultiTimeserie

    !Generic 4D
    type T_Generic4D
        logical                                     :: ON                   = .false.   !initialization: Jauch
        logical                                     :: ReadFromTimeSerie    = .false.   !initialization: Jauch
        integer                                     :: ObjTimeSerie         = null_int  !initialization: Jauch
        integer                                     :: TimeSerieColumn      = null_int  !initialization: Jauch
        real                                        :: CurrentValue         = null_real !initialization: Jauch

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
        integer                                     :: ObjHDF5              =  0
        integer                                     :: NumberOfInstants     = null_int 
        logical                                     :: CyclicTimeON         = .false.
        logical                                     :: From2Dto3D           = .false.
        type(T_Generic4D)                           :: Generic4D
        logical                                     :: ArgumentFileName     = .false. 
        integer                                     :: ObjField4D           = 0
        logical                                     :: Field4D              = .false.
        logical                                     :: HarmonicsON          = .false.
        logical                                     :: SpatialInterpolON    = .false.
        logical                                     :: GenericYear          = .false.
        integer                                     :: Ncells
        real,    dimension(:), pointer              :: X                    => null()
        real,    dimension(:), pointer              :: Y                    => null()        
        real,    dimension(:), pointer              :: Z                    => null()
        real,    dimension(:), pointer              :: Prop                 => null()
        logical, dimension(:), pointer              :: NoData               => null()        
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

!        logical                                     :: AccumulatedValue     = .false.
!        logical                                     :: NoInterpol           = .false.
                
!        integer                                     :: ValuesType           
        logical                                     :: InterpolateValues        = .false. !initialization: Jauch
        logical                                     :: AccumulateValues         = .false. !initialization: Jauch
        logical                                     :: UseOriginalValues        = .false. !initialization: Jauch
        logical                                     :: PreviousInstantValues    = .false. !initialization: Jauch
        logical                                     :: IgnoreNoDataPoint        = .false. !initialization: Jauch
        integer                                     :: PredictDTMethod !1 for old method, 2 for new method (for rain, mainly)
        real                                        :: NoDataValue              = null_real !initialization: Jauch
        
        logical                                     :: Backtracking             = .false.
        
        character(len=StringLength)                 :: OverrideValueKeyword     = null_str
        logical                                     :: OverrideValueKeywordON   = .false.
         
        real                                        :: MinForDTDecrease         = AllmostZero
        real                                        :: DefaultValue             = null_real !initialization: Jauch
        real                                        :: PredictedDT              = -null_real
        real                                        :: DTForNextEvent           = -null_real
        real                                        :: DTForNextDataset     = -null_real
        type(T_Time)                                :: NextEventStart
        type(T_Time)                                :: NextEventEnd
        type(T_Time)                                :: TimeOfNextDataset
        integer                                     :: InstantOfNextDataset
        logical                                     :: ValueIsUsedForDTPrediction = .false.
        real                                        :: NextValueForDTPred
        real,    dimension(:, :   ), pointer        :: Matrix2D         => null()
        real,    dimension(:, :, :), pointer        :: Matrix3D         => null()
        integer, dimension(:, :   ), pointer        :: PointsToFill2D   => null()        
        integer, dimension(:, :, :), pointer        :: PointsToFill3D   => null()
        
        type(T_Time)                                :: BeginTime, EndTime

        !Initialization Methods
        type (T_Layers   )                          :: Layers
        type (T_Boxes    )                          :: Boxes
        type (T_TimeSerie)                          :: TimeSerie
        type (T_ASCIIFile)                          :: ASCIIFile
        type (T_Sponge   )                          :: Sponge
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
                                     Matrix2D, TypeZUV, FileNameHDF, ObjFillMatrix,     &
                                     OverrideValueKeyword, ClientID, PredictDTMethod,   &
                                     MinForDTDecrease, ValueIsUsedForDTPrediction, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: EnterDataID
        integer                                         :: TimeID
        integer                                         :: HorizontalGridID
        integer                                         :: ExtractType
        integer, dimension(:, :), pointer               :: PointsToFill2D
        real, dimension(:, :), pointer                  :: Matrix2D
        integer                                         :: TypeZUV
        type (T_PropertyID)                             :: PropertyID
        integer, optional, intent(IN)                   :: ClientID
        character(*), optional, intent(IN )             :: FileNameHDF, OverrideValueKeyword
        integer,      optional, intent(INOUT)           :: ObjFillMatrix
        integer,      optional, intent(OUT)             :: STAT     
        integer,      optional, intent(IN )             :: PredictDTMethod 
        real,         optional, intent(IN )             :: MinForDTDecrease  
        logical,      optional, intent(IN )             :: ValueIsUsedForDTPrediction

        !Local-------------------------------------------------------------------
        integer                                         :: ready_, STAT_, nUsers, ObjFillMatrix_
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
            
            Me%Matrix2D       => Matrix2D
            Me%PointsToFill2D => PointsToFill2D
            
            where (PointsToFill2D == WaterPoint) Me%Matrix2D = null_real

            if (Me%TypeZUV == TypeU_) then
                Me%Size2D%JUB       = Me%Size2D%JUB + 1
                Me%WorkSize2D%JUB   = Me%WorkSize2D%JUB + 1
            endif

            if (Me%TypeZUV == TypeV_) then
                Me%Size2D%IUB       = Me%Size2D%IUB + 1
                Me%WorkSize2D%IUB   = Me%WorkSize2D%IUB + 1
            endif
            
            if (present(FileNameHDF)) then
            
                Me%HDF%ArgumentFileName = .true.
                Me%HDF%FileName         = trim(FileNameHDF)
                
            else
            
                Me%HDF%ArgumentFileName = .false.
            
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
            nullify(Me%PointsToFill2D)
                 

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ConstructFillMatrix2D - ModuleFillMatrix - ERR02'  

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructFillMatrix2D

    !----------------------------------------------------------------------
 
    subroutine ConstructFillMatrix3D(PropertyID, EnterDataID, TimeID,                   &
                                     HorizontalGridID, GeometryID, ExtractType,         &
                                     PointsToFill3D, Matrix3D, TypeZUV, FillMatrix,     &
                                     FileNameHDF, ObjFillMatrix,                        &
                                     OverrideValueKeyword, ClientID, predictDTMethod,   &
                                     MinForDTDecrease, ValueIsUsedForDTPrediction, STAT)

        !Arguments---------------------------------------------------------------
        type (T_PropertyID)                             :: PropertyID
        integer                                         :: EnterDataID
        integer                                         :: TimeID
        integer                                         :: HorizontalGridID
        integer                                         :: GeometryID
        integer                                         :: ExtractType
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        real, dimension(:, :, :), pointer               :: Matrix3D
        integer                                         :: TypeZUV
        integer, optional, intent(IN)                   :: ClientID
        real        , optional, intent(IN )             :: FillMatrix
        character(*), optional, intent(IN )             :: FileNameHDF, OverrideValueKeyword
        integer,      optional, intent(INOUT)           :: ObjFillMatrix
        integer,      optional, intent(OUT)             :: STAT     
        integer,      optional, intent(IN )             :: PredictDTMethod  
        real,         optional, intent(IN )             :: MinForDTDecrease 
        logical,      optional, intent(IN )             :: ValueIsUsedForDTPrediction 

        !Local-------------------------------------------------------------------
        real                                            :: FillMatrix_
        integer                                         :: ready_, STAT_, nUsers, ObjFillMatrix_
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
            
                Me%HDF%ArgumentFileName = .true.
                Me%HDF%FileName         = trim(FileNameHDF)
                
            else
            
                Me%HDF%ArgumentFileName = .false.
            
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
            
            nullify(Me%Matrix3D      )
            nullify(Me%PointsToFill3D)

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ConstructFillMatrix3D - ModuleFillMatrix - ERR02' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructFillMatrix3D

    !--------------------------------------------------------------------------
    
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
            case default
                write(*,*)'Invalid option for keyword FILE_IN_TIME'
                stop 'ReadOptions - ModuleFillMatrix - ERR020'
        end select
        
        if (Me%HDF%ArgumentFileName) Me%TimeEvolution    = ReadHDF


        if(Me%TimeEvolution == None)then

            call GetData(AuxString, Me%ObjEnterData,  iflag,                            &
                         SearchType     = ExtractType,                                  &
                         keyword        = 'INITIALIZATION_METHOD',                      &
                         default        = "Constant",                                   &
                         ClientModule   = 'ModuleFillMatrix',                           &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR030'


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
                     default        = FillValueReal,                                    &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR090'
        
        if (iflag == 0 .and. .not. Me%HDF%ArgumentFileName) then
        
            if(Me%OverrideValueKeywordON .and. Me%InitializationMethod == Constant)then
            
                call GetData(Me%DefaultValue, Me%ObjEnterData,  iflag,                 &
                             SearchType     = ExtractType,                             &
                             keyword        = trim(Me%OverrideValueKeyword),           &
                             ClientModule   = 'ModuleFillMatrix',                      &
                             default        = FillValueReal,                           &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR100'
                
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

        !Fill Matrix with default Value
        if (Me%Dim == Dim2D) then
            where (PointsToFill2D == WaterPoint) Me%Matrix2D = Me%DefaultValue
        else
            where (PointsToFill3D == WaterPoint) Me%Matrix3D = Me%DefaultValue
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

                        call ConstructSpaceTimeSerie (ExtractType)

                        if (Me%Dim == Dim2D) then
                            call ModifySpaceTimeSerie    (PointsToFill2D = PointsToFill2D) 
                        else
                            call ModifySpaceTimeSerie    (PointsToFill3D = PointsToFill3D) 
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

                call ConstructSpaceTimeSerie (ExtractType)

                if (Me%Dim == Dim2D) then
                    call ModifySpaceTimeSerie    (PointsToFill2D = PointsToFill2D) 
                else
                    call ModifySpaceTimeSerie    (PointsToFill3D = PointsToFill3D) 
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
                
        end select

    end subroutine ReadOptions

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
                                   DefaultValue = Me%DefaultValue,               &
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
                                   DefaultValue = Me%DefaultValue,               &
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


        call InterpolateMatrix3DInTime(ActualTime       = Now,                                   &
                                       Size             = Me%WorkSize3D,                         &
                                       Time1            = Me%ProfileTimeSerie%PreviousTime,      &
                                       Matrix1          = Me%ProfileTimeSerie%PreviousField3D,   &
                                       Time2            = Me%ProfileTimeSerie%NextTime,          &
                                       Matrix2          = Me%ProfileTimeSerie%NextField3D,       &
                                       MatrixOut        = Me%Matrix3D,                           &
                                       PointsToFill3D   = PointsToFill3D)

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


        call InterpolateMatrix3DInTime(ActualTime       = Now,                                   &
                                       Size             = Me%WorkSize3D,                         &
                                       Time1            = Me%ProfileTimeSerie%PreviousTime,      &
                                       Matrix1          = Me%ProfileTimeSerie%PreviousField3D,   &
                                       Time2            = Me%ProfileTimeSerie%NextTime,          &
                                       Matrix2          = Me%ProfileTimeSerie%NextField3D,       &
                                       MatrixOut        = Me%Matrix3D,                           &
                                       PointsToFill3D   = PointsToFill3D)

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
                    Me%Matrix3D(i,j,k) = Me%DefaultValue + CoefA * CellDepth / CoefB
                else if (ProfileType == Exponential) then
                    !Exponential profile
                    Me%Matrix3D(i,j,k) = Me%DefaultValue - CoefA * exp(- CellDepth / CoefB)
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

        !Begin----------------------------------------------------------------
        
        !Gets the name of the data file
        call GetData(Me%ASCIIFile%FileName,                                      &
                     Me%ObjEnterData , iflag,                                    &
                     SearchType   = ExtractType,                                 &
                     keyword      = 'FILENAME',                                  &
                     ClientModule = 'ModuleFillMatrix',                          &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR01'
        if (iflag==0)then
        
           if(Me%OverrideValueKeywordON)then
                
                call GetData(Me%ASCIIFile%FileName,                                &
                             Me%ObjEnterData , iflag,                              &
                             SearchType   = ExtractType,                           &
                             keyword      = trim(Me%OverrideValueKeyword),         &
                             ClientModule = 'ModuleFillMatrix',                    &
                             STAT         = STAT_CALL)                                      
                
                if (STAT_CALL .NE. SUCCESS_)  then
                    write(*,*)trim(Me%PropertyID%Name)//', '//trim(Me%OverrideValueKeyword) 
                    stop 'ReadOptions - ModuleFillMatrix - ERR01'
                end if           
                
                if(iflag == 0)then
                    write(*,*)'Please define the ASCII file in the override keyword '//trim(Me%PropertyID%Name)
                    write(*,*)'to give values for property '//trim(Me%PropertyID%Name)
                    stop 'ReadOptions - ModuleFillMatrix - ERR011'
                end if
                           
            else

                write(*,*) 'ASCII File Name not given not given for property '//trim(Me%PropertyID%Name)           
                stop       'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR02'

            end if        
           end if

        if (Me%Dim == Dim2D) then

            call ConstructGridData(Me%ASCIIFile%GridDataID, Me%ObjHorizontalGrid,        &
                                   FileName     = Me%ASCIIFile%FileName,                 &
                                   DefaultValue = Me%DefaultValue,                       &
                                   STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR03'

            call GetGridDataType(Me%ASCIIFile%GridDataID, TypeZUV, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR04'

            if(TypeZUV .ne. Me%TypeZUV)then
                write(*,*)'Inconsistency found in type ZUV'
                write(*,*)'Grid data: '//trim(Me%ASCIIFile%FileName)
                stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR05'
            end if

            call GetGridData      (Me%ASCIIFile%GridDataID, GridData2D, STAT = STAT_CALL)  
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR06'

            !Copies data
            where (PointsToFill2D == WaterPoint) Me%Matrix2D = GridData2D

            call UnGetGridData    (Me%ASCIIFile%GridDataID, GridData2D, STAT = STAT_CALL)  
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR07'

        else

            call ConstructGridData(Me%ASCIIFile%GridDataID, Me%ObjHorizontalGrid,        &
                                   FileName     = Me%ASCIIFile%FileName,                 &
                                   KLB          = Me%Worksize3D%KLB,                     &
                                   KUB          = Me%Worksize3D%KUB,                     &
                                   DefaultValue = Me%DefaultValue,                       &
                                   STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR08'


            call GetGridDataType(Me%ASCIIFile%GridDataID, TypeZUV, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR09'
            
            if(TypeZUV .ne. Me%TypeZUV)then
                write(*,*)'Inconsistency found in type ZUV'
                write(*,*)'Grid data: '//trim(Me%ASCIIFile%FileName)
                stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR10'
            end if

            call GetGridData      (Me%ASCIIFile%GridDataID, GridData3D, STAT = STAT_CALL)  
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR11'

            !Copies data
            do k = Me%Worksize3D%KLB, Me%Worksize3D%KUB
            do j = Me%Worksize3D%JLB, Me%Worksize3D%JUB
            do i = Me%Worksize3D%ILB, Me%Worksize3D%IUB
                if (PointsToFill3D(i, j, k) == WaterPoint) Me%Matrix3D(i, j, k) = GridData3D(i, j, k)
            enddo
            enddo
            enddo

            call UnGetGridData    (Me%ASCIIFile%GridDataID, GridData3D, STAT = STAT_CALL)  
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR12'

        endif

        call KillGridData (Me%ASCIIFile%GridDataID, STAT = STAT_CALL)  
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR13'

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
        real                                        :: Aux
        integer                                     :: STAT_CALL, i, j, k, l, sp
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: iflag

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

        !Gets the sponge value in the model open boundary
        call GetData(Me%Sponge%OutValue,                                                &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'SPONGE_OUT',                                       &
                     Default      = 1.e5,                                               &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSponge - ModuleFillMatrix - ERR20'
        
        

        !Gets the number of sponge cells
        call GetData(Me%Sponge%Cells,                                                   &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'SPONGE_CELLS',                                     &
                     Default      = 10,                                                 &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSponge - ModuleFillMatrix - ERR30'
        

        !Gets the nsponge evolution
        call GetData(Me%Sponge%Evolution,                                               &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'SPONGE_EVOLUTION',                                 &
                     Default      = sponge_exp_,                                        &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSponge - ModuleFillMatrix - ERR40'

        if      (Me%Sponge%Evolution /= sponge_exp_ .and. Me%Sponge%Evolution /= sponge_linear_) then        

            write(*,*) 'Sponge evolution can only be linear or exponential'
            stop       'ConstructSponge - ModuleFillMatrix - ERR50'
        
        endif
        

        if (Me%DefaultValue < Me%Sponge%OutValue) then
        
            Me%Sponge%Growing = .true.
            
        else
        
            Me%Sponge%Growing = .false.
        
        endif
        
        if (Me%TypeZUV == TypeU_ .or. Me%TypeZUV == TypeV_) then
            allocate(AuxT(Me%Sponge%Cells + 1))
        else
            allocate(AuxT(Me%Sponge%Cells))
        
        endif
        
        if      (Me%Sponge%Evolution == sponge_exp_) then
        
            do sp = 1, Me%Sponge%Cells

                AuxT(sp) = log(Me%Sponge%OutValue) * real(Me%Sponge%Cells - sp) /real(Me%Sponge%Cells - 1) + &
                           log(Me%DefaultValue)  * real(sp - 1)               /real(Me%Sponge%Cells - 1)
                     
                AuxT(sp) = exp(AuxT(sp))
                
            enddo
        
        elseif (Me%Sponge%Evolution == sponge_linear_) then        
        
            do sp = 1, Me%Sponge%Cells

                AuxT(sp) = Me%Sponge%OutValue * real(Me%Sponge%Cells - sp) /real(Me%Sponge%Cells - 1) + &
                           Me%DefaultValue    * real(sp - 1)               /real(Me%Sponge%Cells - 1)
                     
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
            
            AuxI(1) = ILB + sp - 1
            AuxI(2) = IUB - sp + 1
            AuxI(3) = JLB + sp - 1
            AuxI(4) = JUB - sp + 1
            
            
            dij(1,1) = 1 - sp
            dij(2,2) = sp - 1
            dij(3,3) = 1 - sp
            dij(4,4) = sp - 1
            
            
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

    subroutine ConstructSpaceTimeSerie (ExtractType)

        !Arguments-------------------------------------------------------------
        integer                                         :: ExtractType

        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: iflag
        
        !Begin----------------------------------------------------------------

        call GetData(Me%TimeSerie%FileName,                                      &
                     Me%ObjEnterData , iflag,                                    &
                     SearchType   = ExtractType,                                 &
                     keyword      = 'FILENAME',                                  &
                     ClientModule = 'ModuleFillMatrix',                          &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR010'

        if (iflag==0)then
            write(*,*)'Time Serie File Name not given for property'
            write(*,*)trim(Me%PropertyID%Name)
            stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR020'
        endif
 
        call GetData(Me%TimeSerie%Column,                                        &
                     Me%ObjEnterData , iflag,                                    &
                     SearchType   = ExtractType,                                 &
                     keyword      = 'DATA_COLUMN',                               &
                     ClientModule = 'ModuleFillMatrix',                          &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR030'

        if (iflag==0)then
            write(*,*)'Data Column not given for property'
            write(*,*)trim(Me%PropertyID%Name)
            stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR040'
        endif

        !Starts Time Serie
        call StartTimeSerieInput(Me%TimeSerie%ObjTimeSerie, &
                                 Me%TimeSerie%FileName,     &
                                 Me%ObjTime,                &
                                 STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR050'
        
        call GetTimeSerieDataValues(Me%TimeSerie%ObjTimeSerie,     &
                                    Me%TimeSerie%NumberOfInstants, &
                                    STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR060'
        
        !if only one instant is found then values remain constant
        if(Me%TimeSerie%NumberOfInstants == 1) then
            
            Me%RemainsConstant = .true.
            call GetTimeSerieValueForIndex (Me%TimeSerie%ObjTimeSerie,      &
                                            1,                              &
                                            Me%TimeSerie%Column,            &
                                            Me%TimeSerie%CurrentValue,      &
                                            STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR065'             
        
        endif
        
        if (Me%PredictDTMethod == 2) then
        
            if (Me%PropertyID%IDNumber == WindDirection_) then
                write(*,*) 'The method 2 to predict DT do not works with WindDirection_ property'
                stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR070'            
            endif
        
            if(Me%TimeSerie%NumberOfInstants > 1)then
                             
                Me%TimeSerie%PreviousInstant = 1                
                call GetTimeSerieTimeOfDataset(Me%TimeSerie%ObjTimeSerie,    &
                                               Me%TimeSerie%PreviousInstant,     &
                                               Me%TimeSerie%PreviousTime,        &
                                               STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)    &
                    stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR080'
                                                 
                Me%TimeSerie%NextInstant = 1
                Me%TimeSerie%NextTime    = Me%TimeSerie%PreviousTime
                
                call ActualizeTimeSerieValues
                
                Me%NextEventStart        = Me%TimeSerie%PreviousTime
                Me%NextEventEnd          = Me%TimeSerie%PreviousTime
                Me%NextValueForDTPred    = Me%TimeSerie%NextValue
                Me%DTForNextEvent        = 0.0
            
            endif

        endif
        
    end subroutine ConstructSpaceTimeSerie

    !--------------------------------------------------------------------------


    subroutine ConstructHDFInput (ExtractType, ClientID, PointsToFill2D, PointsToFill3D)
    !subroutine ConstructHDFInput (ExtractType, HDF5File, PointsToFill2D, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        integer                                         :: ExtractType
        integer                                         :: ClientID
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D
        !character (len = PathLength),         optional  :: HDF5File       

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: iflag, HDF5_READ
        logical                                         :: exist
        type(T_Time)                                    :: Now, CurrentTime

        !Local-----------------------------------------------------------------
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB, i, j, k
        logical                                         :: FoundSecondInstant, LastGroupEqualField
        real                                            :: Year, Month, Day, Hour, Minute, Second
        real                                            :: LatDefault, LongDefault
        type (T_Size2D)                                 :: WindowLimitsJI
        logical                                         :: MasterOrSlave
        real                                            :: StartTimeYear, EndTimeYear
        
        !Begin-----------------------------------------------------------------

        nullify(Me%HDF%PreviousField2D, Me%HDF%NextField2D)
        nullify(Me%HDF%PreviousField3D, Me%HDF%NextField3D)
        
i0:     if(Me%Dim == Dim2D)then

            ILB = Me%Size2D%ILB
            IUB = Me%Size2D%IUB
            JLB = Me%Size2D%JLB
            JUB = Me%Size2D%JUB

            allocate(Me%HDF%PreviousField2D (ILB:IUB, JLB:JUB))
            allocate(Me%HDF%NextField2D     (ILB:IUB, JLB:JUB))
            allocate(Me%HDF%Array2D         (ILB:IUB, JLB:JUB))

            Me%HDF%PreviousField2D(:,:) = FillValueReal
            Me%HDF%NextField2D    (:,:) = FillValueReal

        else i0

            ILB = Me%Size3D%ILB
            IUB = Me%Size3D%IUB
            JLB = Me%Size3D%JLB
            JUB = Me%Size3D%JUB
            KLB = Me%Size3D%KLB
            KUB = Me%Size3D%KUB

            allocate(Me%HDF%PreviousField3D (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(Me%HDF%NextField3D     (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(Me%HDF%Array3D         (ILB:IUB, JLB:JUB, KLB:KUB))

            Me%HDF%PreviousField3D(:,:,:) = FillValueReal
            Me%HDF%NextField3D    (:,:,:) = FillValueReal

        endif i0
        
        call GetData(Me%HDF%Generic4D%ON,                                               &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = '4D',                                               &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR10'


        if (Me%HDF%Generic4D%ON) call Generic4thDimension(ExtractType)
        

        call GetData(Me%HDF%VGroupPath,                                                 &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'VGROUP_PATH',                                      &
                     default      = "/Results",                                         &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR20'

        call GetData(Me%HDF%MultiplyingFactor,                                          &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'MULTIPLYING_FACTOR',                               &
                     default      = 1.,                                                 &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR30'
        
        if (iflag == 1)then
            Me%HDF%HasMultiplyingFactor = .true.
        end if

        call GetData(Me%HDF%AddingFactor,                                               &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'ADDING_FACTOR',                                    &
                     default      = 0.,                                                 &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR40'
        
        if (iflag == 1)then
            Me%HDF%HasAddingFactor = .true.
        end if

        call GetData(Me%HDF%FieldName,                                                  &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'HDF_FIELD_NAME',                                   &
                     default      = trim(Me%PropertyID%Name),                           &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR50'

        call GetData(LastGroupEqualField,                                               &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'LAST_GROUP_EQUAL_FIELD',                           &
                     default      = .true.,                                             &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR60'

        if (LastGroupEqualField)                                                        &
            Me%HDF%VGroupPath=trim(Me%HDF%VGroupPath)//"/"//trim(Me%HDF%FieldName)


        call GetData(Me%HDF%From2Dto3D,                                                 &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'FROM_2D_TO_3D',                                    &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR70'
        
        if (Me%HDF%From2Dto3D) then
        
            allocate(Me%HDF%ReadField3D(ILB:IUB, JLB:JUB, 0:2), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructHDFInput - ModuleFillMatrix - ERR80'
             
            Me%HDF%ReadField3D(:,:,:) = FillValueReal  
            
        endif      

        if (.not. Me%HDF%ArgumentFileName) then

            call GetData(Me%HDF%FileName,                                               &
                         Me%ObjEnterData , iflag,                                       &
                         SearchType   = ExtractType,                                    &
                         keyword      = 'FILENAME',                                     &
                         ClientModule = 'ModuleFillMatrix',                             &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR90'

            if (iflag==0)then
                write(*,*)'HDF filename not given'
                stop 'ConstructHDFInput - ModuleFillMatrix - ERR100'
            endif
        endif

        inquire (file=trim(Me%HDF%FileName), exist = exist)
        if (.not. exist) then
            write(*,*)'Could not find file '//trim(Me%HDF%FileName)
            stop 'ConstructHDFInput - ModuleFillMatrix - ERR110'
        endif
        

        call GetDDecompParameters(HorizontalGridID = Me%ObjHorizontalGrid, &
                                              MasterOrSlave    = MasterOrSlave,        &
                                              STAT             = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR130'        


        call GetData(Me%HDF%Field4D,                                                    &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'FIELD4D',                                          &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR120'
                
        call GetData(Me%HDF%HarmonicsON,                                                &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'HARMONICS',                                        &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR130'        

        call GetData(Me%HDF%SpatialInterpolON,                                          &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'SPATIAL_INTERPOL',                                 &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR135'        
        
        call GetData(Me%HDF%GenericYear,                                                &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'GENERIC_YEAR',                                     &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR135'        

        
        if (MasterOrSlave) then
            Me%HDF%Field4D = .true. 
        endif
        
if4D:   if (Me%HDF%Field4D) then

ifSI:       if (Me%HDF%SpatialInterpolON) then

                call ConstructField4DInterpol(ExtractType, ClientID, PointsToFill2D, PointsToFill3D)

            else ifSI
            
ifMS:           if (MasterOrSlave) then
                
                    call GetDDecompWorkSize2D(HorizontalGridID = Me%ObjHorizontalGrid, &
                                                          WorkSize         = WindowLimitsJI,       &
                                                          STAT             = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR140'
                                                          
                    write(*,*) 'With domain decomposition - ILB,IUB, JLB, JUB',             &
                                WindowLimitsJI%ILB,WindowLimitsJI%IUB, WindowLimitsJI%JLB, WindowLimitsJI%JUB
                    
                else ifMS
                
                    if(Me%Dim == Dim2D)then
                        WindowLimitsJI     = Me%WorkSize2D
                    else                    
                        WindowLimitsJI%ILB = Me%WorkSize3D%ILB
                        WindowLimitsJI%IUB = Me%WorkSize3D%IUB
                        WindowLimitsJI%JLB = Me%WorkSize3D%JLB
                        WindowLimitsJI%JUB = Me%WorkSize3D%JUB
                    endif
                    
                    write(*,*) 'No domain decomposition - ILB,IUB, JLB, JUB',               &
                                WindowLimitsJI%ILB,WindowLimitsJI%IUB, WindowLimitsJI%JLB, WindowLimitsJI%JUB
                
                endif ifMS                

                call GetLatitudeLongitude(Me%ObjHorizontalGrid, Latitude  = LatDefault,     &
                                                                Longitude = LongDefault,    & 
                                                                STAT      = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR150'

                call ConstructField4D(Field4DID         = Me%HDF%ObjField4D,                &
                                      EnterDataID       = Me%ObjEnterData,                  &
                                      ExtractType       = ExtractType,                      &
                                      FileName          = Me%HDF%FileName,                  &
                                      TimeID            = Me%ObjTime,                       &   
                                      MaskDim           = Me%Dim,                           &
                                      LatReference      = LatDefault,                       &
                                      LonReference      = LongDefault,                      & 
                                      WindowLimitsJI    = WindowLimitsJI,                   &
                                      Extrapolate       = .false.,                          &    
                                      PropertyID        = Me%PropertyID,                    &                                  
                                      ClientID          = ClientID,                         &
                                      STAT              = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR160'
            
            endif ifSI
            
            call GetField4DNumberOfInstants(Me%HDF%ObjField4D, Me%HDF%NumberOfInstants, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR170'
            
        
        else if4D
        
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

            call ConstructHDF5 (Me%HDF%ObjHDF5, trim(Me%HDF%FileName), HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR180'

            call GetHDF5GroupNumberOfItems(Me%HDF%ObjHDF5, "/Time", &
                                           Me%HDF%NumberOfInstants, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR190'
        
        endif if4D
        
        Me%HDF%StartTime = HDF5TimeInstant(1)
        Me%HDF%EndTime   = HDF5TimeInstant(Me%HDF%NumberOfInstants)

        !if only one instant is found then values remain constant
        if(Me%HDF%NumberOfInstants == 1 .and. .not. Me%HDF%HarmonicsON) then
            Me%RemainsConstant = .true.
        endif            

        call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR200'
        
        if (Me%BackTracking) then  
            call BacktrackingTime(Now)
        else   
            Now = CurrentTime
        endif

        !Initial field
i1:     if (Me%HDF%Generic4D%ON .or. (Me%HDF%Field4D .and. Me%HDF%HarmonicsOn)) then
            
            if (Me%Dim == Dim2D) then
                call ModifyHDFInput2D (PointsToFill2D) 
            else
                call ModifyHDFInput3D (PointsToFill3D)
            endif

        else i1
            
i2:         if(Me%HDF%NumberOfInstants > 1)then

i2a:            if (Me%PredictDTMethod == 2) then

                    !This methodology do NOT works with CYCLICTIME or BACKTRACKING (needs revision)
                    if (Me%Backtracking) then
                        stop 'ConstructHDFInput - ModuleFillMatrix - ERR210'
                    endif

                    Me%HDF%PreviousInstant = 1
                    Me%HDF%PreviousTime    = HDF5TimeInstant(Me%HDF%PreviousInstant)
                    
                    if(Me%HDF%PreviousTime .gt. Now)then
                        write(*,*)
                        write(*,*)'Could not read solution from HDF5 file'
                        write(*,*)'First file instant greater than current time'
                        write(*,*)'Matrix name: '//trim(Me%HDF%FieldName)
                        stop      'ConstructHDFInput - ModuleFillMatrix - ERR211'
                    end if
                    
                    if(Me%TimeEvolution .ne. None)then
                        if(Me%HDF%EndTime .lt. Me%EndTime)then
                            write(*,*)
                            write(*,*)'Could not read solution from HDF5 file'
                            write(*,*)'Last instant in file lower than simulation ending time'
                            write(*,*)'Matrix name: '//trim(Me%HDF%FieldName)
                            stop      'ConstructHDFInput - ModuleFillMatrix - ERR212'
                        end if
                    end if                    
                    
                    Me%HDF%NextInstant = 1
                    Me%HDF%NextTime    = Me%HDF%PreviousTime
                    
                    call ActualizeHDFValues (Me%HDF%PreviousInstant, Me%HDF%PreviousField2D)
                    call ActualizeHDFValues (Me%HDF%NextInstant, Me%HDF%NextField2D)
                    
                    Me%NextEventStart        = Me%HDF%PreviousTime
                    Me%NextEventEnd          = Me%HDF%PreviousTime
                    
                    Me%DTForNextEvent        = 0.0
                    
                    if (Me%Dim == Dim2D) then
                        Me%NextValueForDTPred = maxval(Me%HDF%NextField2D)
                        call ModifyHDFInput2D (PointsToFill2D) 
                    else
                        Me%NextValueForDTPred = maxval(Me%HDF%NextField3D)
                        call ModifyHDFInput3D (PointsToFill3D)
                    endif
                
                else i2a

                    if (Me%Backtracking) then
                        Me%HDF%PreviousInstant  = Me%HDF%NumberOfInstants
                        Me%HDF%NextInstant      = Me%HDF%PreviousInstant               
                    else
                        Me%HDF%PreviousInstant  = 1
                        Me%HDF%NextInstant      = Me%HDF%PreviousInstant
                    endif
                    
                    Me%HDF%PreviousTime     = HDF5TimeInstant(Me%HDF%PreviousInstant)
                    
                    if(Me%HDF%GenericYear)then
                        call SetHDFGenericYear(Me%HDF%PreviousTime, RefTime = Now)
                    endif

ib:                 if (Me%BackTracking) then  
                            
                        if(Me%HDF%PreviousTime .lt. Now)then
                            write(*,*)
                            write(*,*)'----------Backtracking mode-----------'
                            write(*,*)'Could not read solution from HDF5 file'
                            write(*,*)'Last file instant greater than current time'
                            write(*,*)'Matrix name: '//trim(Me%HDF%FieldName)
                            stop      'ConstructHDFInput - ModuleFillMatrix - ERR220'
                        end if

                        if(Me%TimeEvolution .ne. None)then                        
                            if(Me%HDF%StartTime .gt. Me%BeginTime)then
                                write(*,*)
                                write(*,*)'----------Backtracking mode-----------'                                
                                write(*,*)'Could not read solution from HDF5 file'
                                write(*,*)'First instant in file lower than simulation starting time'
                                write(*,*)'Matrix name: '//trim(Me%HDF%FieldName)
                                stop      'ConstructHDFInput - ModuleFillMatrix - ERR230'
                            end if
                        endif
                        
                    else   ib
                    
                        if(Me%HDF%PreviousTime .gt. Now)then
                            write(*,*)
                            write(*,*)'Could not read solution from HDF5 file'
                            write(*,*)'First file instant greater than current time'
                            write(*,*)'Matrix name: '//trim(Me%HDF%FieldName)
                            stop      'ConstructHDFInput - ModuleFillMatrix - ERR240'
                        end if

                        if(Me%HDF%GenericYear)then
                        
                            call SetHDFGenericYear(Me%HDF%EndTime,   RefTime = Me%EndTime)
                            call SetHDFGenericYear(Me%HDF%StartTime, RefTime = Me%BeginTime)
                            
                            call ExtractDate(Me%HDF%StartTime, Year = StartTimeYear)
                            call ExtractDate(Me%HDF%EndTime,   Year = EndTimeYear  )
                            
                            if(StartTimeYear .ne. EndTimeYear)then
                                write(*,*)
                                write(*,*)'When using a generic year HDF5 file'
                                write(*,*)'The year of the start time has to be the same as'
                                write(*,*)'the year of the end time'
                                write(*,*)'Matrix name: '//trim(Me%HDF%FieldName)
                                stop      'ConstructHDFInput - ModuleFillMatrix - ERR245'
                            endif
                            
                        endif


                        if(Me%TimeEvolution .ne. None)then
                            if(Me%HDF%EndTime .lt. Me%EndTime)then
                                write(*,*)
                                write(*,*)'Could not read solution from HDF5 file'
                                write(*,*)'Last instant in file lower than simulation ending time'
                                write(*,*)'Matrix name: '//trim(Me%HDF%FieldName)
                                stop      'ConstructHDFInput - ModuleFillMatrix - ERR250'
                            end if
                        end if
                    endif ib

                    FoundSecondInstant = .false.
                
                    !if number of instants greater than 1 then 
                    !find first and second instants
d2:                 do while(.not. FoundSecondInstant)

                        Me%HDF%PreviousInstant  = Me%HDF%NextInstant
                        if (Me%Backtracking) then
                            Me%HDF%NextInstant      = Me%HDF%NextInstant - 1
                        else                
                            Me%HDF%NextInstant      = Me%HDF%NextInstant + 1
                        endif

                        if (Me%HDF%CyclicTimeON .and. Me%HDF%NextInstant .gt. Me%HDF%NumberOfInstants) then
                            Me%HDF%NextInstant  = 1
                        end if


                        Me%HDF%NextTime         = HDF5TimeInstant(Me%HDF%NextInstant)
                        
                        if(Me%HDF%GenericYear)then
                            call SetHDFGenericYear(Me%HDF%NextTime, Now)
                        endif

                        if (Me%Backtracking) then
                            if(Me%HDF%PreviousTime .ge. Now .and. Me%HDF%NextTime .le. Now) then
                                FoundSecondInstant  = .true.
                                exit
                            end if
                        else
                            if(Me%HDF%PreviousTime .le. Now .and. Me%HDF%NextTime .ge. Now) then
                                FoundSecondInstant  = .true.
                                exit
                            end if
                        endif
                        Me%HDF%PreviousTime = Me%HDF%NextTime

                        if (Me%Backtracking) then
                            if(Me%HDF%NextInstant .lt. 1 .and. .not. Me%HDF%CyclicTimeON) then
                                write(*,*)
                                write(*,*)'----------Backtracking mode-----------------'
                                write(*,*)'Could not read solution from HDF5 file'
                                write(*,*)'Could not find second instant in file'
                                write(*,*)'Matrix name: '//trim(Me%HDF%FieldName)
                                stop      'ConstructHDFInput - ModuleFillMatrix - ERR260'
                            end if
                        
                        else
                            if(Me%HDF%NextInstant .gt. Me%HDF%NumberOfInstants .and. .not. Me%HDF%CyclicTimeON) then
                                write(*,*)
                                write(*,*)'Could not read solution from HDF5 file'
                                write(*,*)'Could not find second instant in file'
                                write(*,*)'Matrix name: '//trim(Me%HDF%FieldName)
                                stop      'ConstructHDFInput - ModuleFillMatrix - ERR270'
                            end if
                        endif
                    end do d2

                endif i2a

            elseif(Me%HDF%NumberOfInstants == 1)then i2

                Me%HDF%PreviousInstant  = 1
                Me%HDF%NextInstant      = Me%HDF%PreviousInstant

                Me%HDF%PreviousTime     = HDF5TimeInstant(Me%HDF%PreviousInstant)
                
                if(Me%HDF%GenericYear)then
                    call SetHDFGenericYear(Me%HDF%PreviousTime, Now)
                endif

                Me%HDF%NextTime         = Me%HDF%PreviousTime

                call ExtractDate(Me%HDF%PreviousTime, Year, Month, Day, Hour, Minute, Second)

                write(*,*)
                write(*,*)trim(Me%HDF%FieldName)//' is being read from HDF file.' 
                write(*,*)'Time instant: ', Year, Month, Day, Hour, Minute, Second
                write(*,*)'ConstructHDFInput - ModuleFillMatrix - WRN10'
            
            else i2
                write(*,*)
                write(*,*)'Could not read solution from HDF5 file'
                write(*,*)'No time information found'
                stop 'ConstructHDFInput - ModuleFillMatrix - ERR300'

            end if i2

            if (Me%PredictDTMethod == 1) then
i4:         if(Me%Dim == Dim2D)then

                ILB = Me%Size2D%ILB
                IUB = Me%Size2D%IUB
                JLB = Me%Size2D%JLB
                JUB = Me%Size2D%JUB

                allocate(Me%HDF%PreviousField2D (ILB:IUB, JLB:JUB))
                allocate(Me%HDF%NextField2D     (ILB:IUB, JLB:JUB))

                Me%HDF%PreviousField2D(:,:) = FillValueReal
                Me%HDF%NextField2D    (:,:) = FillValueReal

                call ReadHDF5Values2D(Me%HDF%PreviousInstant, Me%HDF%PreviousField2D)
                call ReadHDF5Values2D(Me%HDF%NextInstant,     Me%HDF%NextField2D    )
                
                !limit maximum values
                do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
                do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                
#ifndef _NOT_IEEE_ARITHMETIC
                        if (ieee_is_nan (Me%HDF%PreviousField2D(i,j)))                      &
                            Me%HDF%PreviousField2D(i,j) = FillValueReal
#endif
                
                    if (abs(Me%HDF%PreviousField2D   (i,j)) > abs(FillValueReal))       &
                        Me%HDF%PreviousField2D(i,j) = FillValueReal

#ifndef _NOT_IEEE_ARITHMETIC
                    if (ieee_is_nan (Me%HDF%NextField2D    (i,j)))                      &
                        Me%HDF%NextField2D(i,j)     = FillValueReal                
#endif
                    
                    if (abs(Me%HDF%NextField2D       (i,j)) > abs(FillValueReal))       &
                        Me%HDF%NextField2D(i,j)     = FillValueReal
                enddo
                    enddo
                                
                if (Me%HDF%PreviousInstant /= Me%HDF%NextInstant) then
                
                    !Interpolates the two matrixes in time
                    call InterpolateMatrix2DInTime(ActualTime       = Now,                         &
                                                   Size             = Me%WorkSize2D,               &
                                                   Time1            = Me%HDF%PreviousTime,         &
                                                   Matrix1          = Me%HDF%PreviousField2D,      &
                                                   Time2            = Me%HDF%NextTime,             &
                                                   Matrix2          = Me%HDF%NextField2D,          &
                                                   MatrixOut        = Me%Matrix2D,                 &
                                                   PointsToFill2D   = PointsToFill2D)
                                                                              
                else

                    Me%Matrix2D(:,:)  = Me%HDF%PreviousField2D(:,:)

                endif

            else i4

                ILB = Me%Size3D%ILB
                IUB = Me%Size3D%IUB
                JLB = Me%Size3D%JLB
                JUB = Me%Size3D%JUB
                KLB = Me%Size3D%KLB
                KUB = Me%Size3D%KUB

                allocate(Me%HDF%PreviousField3D (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%HDF%NextField3D     (ILB:IUB, JLB:JUB, KLB:KUB))

                Me%HDF%PreviousField3D(:,:,:) = FillValueReal
                Me%HDF%NextField3D    (:,:,:) = FillValueReal

                call ReadHDF5Values3D(Me%HDF%PreviousInstant, Me%HDF%PreviousField3D)
                call ReadHDF5Values3D(Me%HDF%NextInstant,     Me%HDF%NextField3D    )

                !limit maximum values
                do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
                do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB

#ifndef _NOT_IEEE_ARITHMETIC
                    if (ieee_is_nan (Me%HDF%PreviousField3D(i,j,k)))                    &
                        Me%HDF%PreviousField3D(i,j,k) = FillValueReal 
#endif
                
                    if (abs(Me%HDF%PreviousField3D(i,j,k)) > abs(FillValueReal))        &
                        Me%HDF%PreviousField3D(i,j,k) = FillValueReal

#ifndef _NOT_IEEE_ARITHMETIC
                    if (ieee_is_nan (Me%HDF%NextField3D    (i,j,k)))                    &
                        Me%HDF%NextField3D(i,j,k) = FillValueReal 
#endif
                    
                    if (abs(Me%HDF%NextField3D    (i,j,k)) > abs(FillValueReal))        &
                        Me%HDF%NextField3D        (i,j,k) = FillValueReal
                enddo
                enddo
                enddo                


                if (Me%HDF%PreviousInstant /= Me%HDF%NextInstant) then
                
                   if (Me%PreviousInstantValues) then
                    
                        Me%Matrix3D = Me%HDF%PreviousField3D
  
                    else                

                        call InterpolateMatrix3DInTime(ActualTime       = Now,                         &
                                                       Size             = Me%WorkSize3D,               &
                                                       Time1            = Me%HDF%PreviousTime,         &
                                                       Matrix1          = Me%HDF%PreviousField3D,      &
                                                       Time2            = Me%HDF%NextTime,             &
                                                       Matrix2          = Me%HDF%NextField3D,          &
                                                       MatrixOut        = Me%Matrix3D,                 &
                                                       PointsToFill3D   = PointsToFill3D)
                    endif
                    
                else

                    !Prev and next are equal (last instant?)
                    Me%Matrix3D(:,:,:)  = Me%HDF%NextField3D(:,:,:)

                endif

            end if i4
            endif

        endif i1


    end subroutine ConstructHDFInput

    !-----------------------------------------------------------------------------------
    
    subroutine ConstructField4DInterpol(ExtractType, ClientID, PointsToFill2D, PointsToFill3D)
        !Arguments----------------------------------------------------------------------
        integer                                         :: ExtractType, ClientID
        integer,  dimension(:,:,:),   pointer, optional :: PointsToFill3D
        integer,  dimension(:,:),     pointer, optional :: PointsToFill2D


        !Local--------------------------------------------------------------------------
        real,       dimension(:,:),     pointer         :: CoordX, CoordY
        real, dimension(1:2,1:2)                        :: WindowLimitsXY
        real                                            :: West, East, South, North  
        real                                            :: LatDefault, LongDefault
        integer                                         :: STAT_CALL, i, j, k, icount, NCells
        
        !Begin--------------------------------------------------------------------------      

        call GetGridBorderLimits(Me%ObjHorizontalGrid, West, East, South, North, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR10'

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
            Me%HDF%Ncells = Ncells
        
        else
        
            icount = 0
            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB        
                if (PointsToFill3D(i,j,k) == WaterPoint) icount = icount + 1
            enddo
            enddo
            enddo
            
            Ncells        = icount
            Me%HDF%Ncells = Ncells
        
            allocate(Me%HDF%Z(1:NCells))
         
            Me%HDF%Z(1:NCells) = FillValueReal
                                                     
        endif                    
            
        allocate(Me%HDF%X(1:NCells), Me%HDF%Y(1:NCells), Me%HDF%Prop(1:NCells), Me%HDF%NoData(1:NCells))
        
        Me%HDF%X     (1:NCells) = FillValueReal
        Me%HDF%Y     (1:NCells) = FillValueReal
        Me%HDF%Prop  (1:NCells) = FillValueReal
        Me%HDF%NoData(1:NCells) = .true.
        
        call GetZCoordinates(Me%ObjHorizontalGrid, CoordX, CoordY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR40'
        
        if (Me%Dim == Dim2D) then
            
            icount = 0
            
            do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB        
            
                if (PointsToFill2D(i,j) == WaterPoint) then

                    icount           = icount + 1
                    Me%HDF%X(icount) = CoordX(i, j)
                    Me%HDF%Y(icount) = CoordY(i, j)

                endif                    
                
            enddo
            enddo
            
        else
        
            icount = 0
            
            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB        
                
                if (PointsToFill3D(i,j,k) == WaterPoint) then
                    
                    icount           = icount + 1
                    Me%HDF%X(icount) = CoordX(i, j)
                    Me%HDF%Y(icount) = CoordY(i, j)
                
                endif                    
                
            enddo
            enddo   
            enddo     
            
        endif       
    
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR70'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR80'
        
        call GetLatitudeLongitude(Me%ObjHorizontalGrid, Latitude  = LatDefault,         &
                                                        Longitude = LongDefault,        & 
                                                        STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR90'
        
        call ConstructField4D(Field4DID         = Me%HDF%ObjField4D,                    &
                              EnterDataID       = Me%ObjEnterData,                      &
                              ExtractType       = ExtractType,                          &
                              FileName          = Me%HDF%FileName,                      &
                              TimeID            = Me%ObjTime,                           &   
                              MaskDim           = Me%Dim,                               &
                              LatReference      = LatDefault,                           &
                              LonReference      = LongDefault,                          & 
                              WindowLimitsXY    = WindowLimitsXY,                       &
                              Extrapolate       = .true.,                               &    
                              PropertyID        = Me%PropertyID,                        &                                  
                              ClientID          = ClientID,                             &
                              STAT              = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructField4DInterpol - ModuleFillMatrix - ERR100'
    
    end subroutine ConstructField4DInterpol

    
   !----------------------------------------------------------------------------

    subroutine Generic4thDimension(ExtractType)

        !Arguments-------------------------------------------------------------
        integer                            :: ExtractType

        !Local-----------------------------------------------------------------
        character(len = StringLength)      :: Filename
        integer                            :: STAT_CALL, iflag

        !----------------------------------------------------------------------

        call GetData(Me%HDF%Generic4D%ReadFromTimeSerie,                                &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = '4D_TIME_SERIE',                                    &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'Generic4thDimension - ModuleFillMatrix - ERR10'


        if (Me%HDF%Generic4D%ReadFromTimeSerie) then

            call GetData(Filename, Me%ObjEnterData, iflag,                              &
                         keyword        = 'GENERIC_4D_FILENAME',                        &  
                         SearchType     = ExtractType,                                  &
                         ClientModule   = 'ModuleFillMatrix',                           &
                         default        = "******.***",                                 &
                         STAT           = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'Generic4thDimension  - ModuleFillMatrix - ERR20'
            if (iflag == 0) stop 'Generic4thDimension  - ModuleFillMatrix - ERR30'

            call GetData(Me%HDF%Generic4D%TimeSerieColumn, Me%ObjEnterData, iflag,      &
                         keyword        = 'TIME_SERIE_COLUMN',                          &  
                         SearchType     = ExtractType,                                  &
                         ClientModule   = 'ModuleFillMatrix',                           &
                         default        = 2,                                            &
                         STAT           = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'Generic4thDimension  - ModuleFillMatrix - ERR40'

            !Starts Time Serie
            call StartTimeSerieInput(Me%HDF%Generic4D%ObjTimeSerie,                     &
                                     FileName,                                          &
                                     CheckDates =.false.,                               &
                                     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Generic4thDimension - ModuleFillMatrix - ERR50'

        endif


     end subroutine Generic4thDimension
   !----------------------------------------------------------------------------


    type(T_Time) function HDF5TimeInstant(Instant)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        

        !Local-----------------------------------------------------------------
!        type(T_Time)                            :: TimeInstant
        real,    dimension(:), pointer          :: TimeVector
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------

if4D:   if (Me%HDF%Field4D) then      
        
            HDF5TimeInstant = GetField4DInstant(Field4DID = Me%HDF%ObjField4D,          &
                                                Instant   = Instant,                    &
                                                STAT      =  STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleFillMatrix - ERR10'
        
        else if4D
        
            call HDF5SetLimits  (Me%HDF%ObjHDF5, 1, 6, STAT = STAT_CALL)

            allocate(TimeVector(6))

            call HDF5ReadData   (HDF5ID         = Me%HDF%ObjHDF5,                           &
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
                                      
        if(Me%HDF%GenericYear)Year = CyclicTime
        
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
    
    subroutine SetHDFGenericYear(TimeInstant, RefTime, AddYear)

        !Arguments-------------------------------------------------------------
        type(T_Time)                            :: TimeInstant
        type(T_Time)                            :: RefTime
        logical, optional                       :: AddYear

        !Local-----------------------------------------------------------------
        real                                    :: Year, Month, Day, Hour, Minute, Second
        logical                                 :: AddYear_

        !Begin-----------------------------------------------------------------

        call ExtractDate(TimeInstant, Year = Year, Month  = Month,  Day    = Day,                &
                                      Hour = Hour, Minute = Minute, Second = Second)
 
        call ExtractDate(RefTime, Year = Year)
        
        if(present(AddYear))then
            AddYear_ = AddYear
        else
            AddYear_ = .false.
        end if
        
        if(AddYear_)Year = Year + 1

        call SetDate(TimeInstant, Year = Year, Month  = Month,  Day    = Day,     &
                                  Hour = Hour, Minute = Minute, Second = Second)


    end subroutine SetHDFGenericYear
    
    !--------------------------------------------------------------------------

    real function HDF5Generic4DInstant(Instant)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        

        !Local-----------------------------------------------------------------
        real,   dimension(:), pointer           :: AuxVector
        type (T_Time)                           :: TimeInstant 
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (Me%HDF%Generic4D%ReadFromTimeSerie) then

            TimeInstant = HDF5TimeInstant(Instant)

            HDF5Generic4DInstant = TimeSerieValue(Me%HDF%Generic4D%ObjTimeSerie,        &
                                                  TimeInstant,                          &
                                                  Me%HDF%Generic4D%TimeSerieColumn) 

        else
                
            call HDF5SetLimits  (Me%HDF%ObjHDF5, 1, 1, STAT = STAT_CALL)

            allocate(AuxVector(1))

            call HDF5ReadData   (HDF5ID         = Me%HDF%ObjHDF5,                       &
                                 GroupName      = "/Generic4D",                         &
                                 Name           = "Generic4D",                          &
                                 Array1D        = AuxVector,                            &
                                 OutputNumber   = Instant,                              &
                                 STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Generic4DInstant - ModuleFillMatrix - ERR10'

            HDF5Generic4DInstant = AuxVector(1)
 
            deallocate(AuxVector)

        endif

    end function HDF5Generic4DInstant
    
    !--------------------------------------------------------------------------
    subroutine ReadHDF5Values2D (InstantIn, Field)

        !Arguments-------------------------------------------------------------
        integer                                 :: InstantIn
        real, dimension(:,:), pointer           :: Field
        
        !Local-----------------------------------------------------------------
        integer                                 :: Instant, STAT_CALL, Imax, Jmax, i, j
        type (T_Time)                           :: CurrentTime

        !Begin-----------------------------------------------------------------
        
        if (Me%Backtracking) then
            Instant = Me%HDF%NumberOfInstants - InstantIn + 1
        else
            Instant = InstantIn
        endif
if4D:   if (Me%HDF%Field4D) then

            CurrentTime = GetField4DInstant (Me%HDF%ObjField4D, Instant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR10'
            
            if (Me%HDF%SpatialInterpolON) then
            
                call ModifyField4DInterpol(CurrentTime      = CurrentTime,              & 
                                           Matrix2D         = Field)

            else
            
                call ModifyField4D(Field4DID        = Me%HDF%ObjField4D,                &
                                   PropertyIDNumber = Me%PropertyID%IDNumber,           & 
                                   CurrentTime      = CurrentTime,                      &
                                   Matrix2D         = Field,                            &
                                   STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR10'

            endif                

        else if4D            
        
            call GetHDF5ArrayDimensions(Me%HDF%ObjHDF5, trim(Me%HDF%VGroupPath),        &
                              trim(Me%HDF%FieldName), OutputNumber = Instant,           &
                              Imax = Imax, Jmax = Jmax, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR20'                               
            
            if ((Imax /= Me%WorkSize2D%IUB - Me%WorkSize2D%ILB + 1) .or.                &
                (Jmax /= Me%WorkSize2D%JUB - Me%WorkSize2D%JLB + 1)) then
                
                write (*,*) trim(Me%HDF%VGroupPath)
                write (*,*) trim(Me%HDF%FieldName)
                write (*,*) 'miss match between the HDF5 input file and model domain'
                stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR30'                                   

            endif

            call HDF5SetLimits  (Me%HDF%ObjHDF5, Me%WorkSize2D%ILB, Me%WorkSize2D%IUB,      &
                                 Me%WorkSize2D%JLB, Me%WorkSize2D%JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR40'


            call HDF5ReadData(Me%HDF%ObjHDF5, trim(Me%HDF%VGroupPath),                      &
                              trim(Me%HDF%FieldName),                                       &
                              Array2D = Field, OutputNumber = Instant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR50'

        endif if4D

        if(Me%HDF%HasMultiplyingFactor)then
            do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                Field(i,j) = Field(i,j) * Me%HDF%MultiplyingFactor
            enddo
            enddo
        end if

        if(Me%HDF%HasAddingFactor)then
            do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                Field(i,j) = Field(i,j) + Me%HDF%AddingFactor
            enddo
            enddo
        end if

    end subroutine ReadHDF5Values2D
    
    
    !--------------------------------------------------------------------------

    
    subroutine ReadHDF5Values3D (Instant, Field)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        real, dimension(:,:,:), pointer         :: Field

        !Local-----------------------------------------------------------------
        type (T_Time)                           :: CurrentTime
        integer                                 :: Imax, Jmax, Kmax
        integer                                 :: STAT_CALL, i, j, k, ILB, IUB, JLB, JUB, KLB, KUB

        !Begin-----------------------------------------------------------------


         ILB = Me%WorkSize3D%ILB
         IUB = Me%WorkSize3D%IUB
         JLB = Me%WorkSize3D%JLB
         JUB = Me%WorkSize3D%JUB

        if (Me%HDF%From2Dto3D) then
            KLB = 1
            KUB = 1
            Kmax= 1
        else
            KLB = Me%WorkSize3D%KLB
            KUB = Me%WorkSize3D%KUB
                     
            Me%HDF%ReadField3D => Field
        endif
        
        
if4D:   if (Me%HDF%Field4D) then

            CurrentTime = GetField4DInstant (Me%HDF%ObjField4D, Instant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR10'
            
            if (Me%HDF%SpatialInterpolON) then
            
                call ModifyField4DInterpol(CurrentTime      = CurrentTime,              & 
                                           Matrix3D         = Field)

            else            

                call ModifyField4D(Field4DID        = Me%HDF%ObjField4D,                &
                                   PropertyIDNumber = Me%PropertyID%IDNumber,           & 
                                   CurrentTime      = CurrentTime,                      & 
                                   Matrix3D         = Field,                            &
                                   STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR20'                                   
                
            endif                

        else if4D            
        
            if (Me%HDF%From2Dto3D) then
                call GetHDF5ArrayDimensions(Me%HDF%ObjHDF5, trim(Me%HDF%VGroupPath),    &
                                  trim(Me%HDF%FieldName), OutputNumber = Instant,       &
                                  Imax = Imax, Jmax = Jmax, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR30'
            else
                call GetHDF5ArrayDimensions(Me%HDF%ObjHDF5, trim(Me%HDF%VGroupPath),    &
                                  trim(Me%HDF%FieldName), OutputNumber = Instant,       &
                                  Imax = Imax, Jmax = Jmax, Kmax = Kmax, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR31'                                   
            endif
            
            
            if ((Imax /= IUB - ILB + 1) .or.                                                &
                (Jmax /= JUB - JLB + 1) .or.                                                &
                (Kmax /= KUB - KLB + 1)) then
                
                if (.not.(Kmax == 0 .and. KUB-KLB == 0)) then
                
                    write (*,*) trim(Me%HDF%VGroupPath)
                    write (*,*) trim(Me%HDF%FieldName)
                    write (*,*) 'miss match between the HDF5 input file and model domain'
                    stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR40'                                   
                
                endif

            endif
          


            call HDF5SetLimits  (Me%HDF%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR50'
            
                 
            call HDF5ReadData(Me%HDF%ObjHDF5, trim(Me%HDF%VGroupPath),                      &
                              trim(Me%HDF%FieldName),                                       &
                              Array3D = Me%HDF%ReadField3D, OutputNumber = Instant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR60'
            
        endif if4D            

        if (Me%HDF%From2Dto3D) then    
           
            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j =               JLB,               JUB
            do i =               ILB,               IUB
                Field(i,j,k) = Me%HDF%ReadField3D(i,j,1)
            enddo
            enddo
            enddo
            
        else
           nullify(Me%HDF%ReadField3D)  
        endif        

        if(Me%HDF%HasMultiplyingFactor)then
            
            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
                Field(i,j,k) = Field(i,j,k) * Me%HDF%MultiplyingFactor
            enddo
            enddo
            enddo
            
        end if

        if(Me%HDF%HasAddingFactor)then
            
            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
                Field(i,j,k) = Field(i,j,k) + Me%HDF%AddingFactor
            enddo
            enddo
            enddo

        end if

    end subroutine ReadHDF5Values3D

    !-----------------------------------------------------------------------------------

    subroutine ModifyField4DInterpol(CurrentTime, Matrix3D, Matrix2D)
    
        !Arguments-------------------------------------------------------------
        type(T_Time)                                :: CurrentTime
        real, dimension(:,:,:), pointer, optional   :: Matrix3D
        real, dimension(:,:  ), pointer, optional   :: Matrix2D
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), pointer             :: ZCellCenter 
        integer                                     :: i, j, k, icount
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------

        Me%HDF%NoData   (:) = .true.

if2D:   if (Me%Dim == Dim2D) then

            call ModifyField4DXYZ(Field4DID             = Me%HDF%ObjField4D,            &
                                  PropertyIDNumber      = Me%PropertyID%IDNumber,       &
                                  CurrentTime           = CurrentTime,                  &
                                  X                     = Me%HDF%X,                     &
                                  Y                     = Me%HDF%Y,                     &
                                  Field                 = Me%HDF%Prop,                  &
                                  NoData                = Me%HDF%NoData,                &
                                  STAT                  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR10' 
            
            icount = 0
            
            do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB        
                
                if (Me%PointsToFill2D(i,j) == WaterPoint) then
                    
                    icount           = icount + 1
                    if (Me%HDF%NoData(icount)) then
                        write(*,*) GetPropertyName (Me%PropertyID%IDNumber)
                        write(*,*) 'No data in 2D cell I=',i, 'J=',j
                        stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR20' 
                    else                        
                        Matrix2D(i, j)   = Me%HDF%Prop(icount)
                    endif
                endif                    
                
            enddo   
            enddo      
            
        
        else if2D
        
            call GetGeometryDistances(Me%ObjGeometry, ZCellCenter, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR30' 
        
            icount = 0
            
            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB        
                
                if (Me%PointsToFill3D(i,j,k) == WaterPoint) then
                    
                    icount           = icount + 1
                    Me%HDF%Z(icount) = ZCellCenter(i, j, k)
                
                endif                    
                
            enddo
            enddo   
            enddo            
            
            call UnGetGeometry(Me%ObjGeometry, ZCellCenter, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR40'             
        
            call ModifyField4DXYZ(Field4DID             = Me%HDF%ObjField4D,            &
                                  PropertyIDNumber      = Me%PropertyID%IDNumber,       &
                                  CurrentTime           = CurrentTime,                  &
                                  X                     = Me%HDF%X,                     &
                                  Y                     = Me%HDF%Y,                     &
                                  Z                     = Me%HDF%Z,                     &
                                  Field                 = Me%HDF%Prop,                  &
                                  NoData                = Me%HDF%NoData,                &
                                  STAT                  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR50' 
            
            icount = 0
            
            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB        
                
                if (Me%PointsToFill3D(i,j,k) == WaterPoint) then
                    
                    icount           = icount + 1
                    if (Me%HDF%NoData(icount)) then
                        write(*,*) 'No data in 3D cell I=',i, 'J=',j, 'K=',k
                        stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR60' 
                    else                        
                        Matrix3D(i, j, k)   = Me%HDF%Prop(icount)
                    endif
                
                endif                    
                
            enddo
            enddo   
            enddo            
            
        
        endif if2D
        
        if (icount /= Me%HDF%Ncells) then
            stop 'ModifyField4DInterpol - ModuleFillMatrix - ERR70' 
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

    subroutine GetDefaultValue (FillMatrixID, DefaultValue, STAT)

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

            DefaultValue = Me%DefaultValue

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetDefaultValue

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
                
                StartTime = Me%HDF%StartTime
                EndTime   = Me%HDF%EndTime
                
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

            NumberOfInstants = Me%HDF%NumberOfInstants

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
        
        !----------------------------------------------------------------------
    
    end subroutine GetNumberOfInstants
    
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

            call HDF5SetLimits  (Me%HDF%ObjHDF5, 1, 6, STAT = STAT_)

            allocate(TimeVector(6))

            call HDF5ReadData   (HDF5ID         = Me%HDF%ObjHDF5,   &
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
        
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyFillMatrix(FillMatrixID, Matrix2D, Matrix3D, PointsToFill2D,       &
                                PointsToFill3D, Generic_4D_Value, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        real,    dimension(:, :),    pointer, optional  :: Matrix2D
        real,    dimension(:, :, :), pointer, optional  :: Matrix3D
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D
        real,                    intent( IN), optional  :: Generic_4D_Value
        integer,                 intent(OUT), optional  :: STAT

        !Local-----------------------------------------------------------------
        real                                            :: Generic_4D_Value_        
        integer                                         :: STAT_, ready_
        logical                                         :: ModifyError = .false.
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            if (present(Matrix2D)) Me%Matrix2D => Matrix2D
            if (present(Matrix3D)) Me%Matrix3D => Matrix3D

            if (present(PointsToFill2D)) Me%PointsToFill2D => PointsToFill2D
            if (present(PointsToFill3D)) Me%PointsToFill3D => PointsToFill3D



            if (present(Generic_4D_Value)) then
                Generic_4D_Value_ = Generic_4D_Value
            else
                if (Me%HDF%Generic4D%ON) then
                    ModifyError = .true. 
                    write(*,*) 'The FillMatrix wants to interpolate along a Generic 4D dimension'
                    write(*,*) 'However, no data is provide for the interpolation'
                endif
                Generic_4D_Value_ = FillValueReal
            endif

            if (.not. ModifyError) then

                select case (Me%TimeEvolution)

                    case (ReadTimeSerie)

                        if (Me%Dim == Dim2D) then
                            call ModifySpaceTimeSerie    (PointsToFill2D = PointsToFill2D) 
                        else
                            call ModifySpaceTimeSerie    (PointsToFill3D = PointsToFill3D) 
                        endif
                
                    case (ReadHDF)

                        if(.not. Me%RemainsConstant)then

                            if (Me%Dim == Dim2D) then
                                call ModifyHDFInput2D (PointsToFill2D, Generic_4D_Value_) 
                            else
                                call ModifyHDFInput3D (PointsToFill3D, Generic_4D_Value_)
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

                end select
                
                nullify(Me%Matrix2D)
                nullify(Me%Matrix3D)

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

    end subroutine ModifyFillMatrix

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

    subroutine ModifyHDFInput3D(PointsToFill3D, Generic_4D_Value_)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        real, optional                                  :: Generic_4D_Value_

        !Local----------------------------------------------------------------

        !Begin----------------------------------------------------------------

        if (Me%HDF%Generic4D%ON) then

            if (.not. present(Generic_4D_Value_)) &
                stop 'ModifyHDFInput3D - ModuleFillMatrix - ERR010'
            call ModifyHDFInput3DGeneric4D(PointsToFill3D, Generic_4D_Value_)

        else
            
            call ModifyHDFInput3DTime(PointsToFill3D)

        endif

    end subroutine ModifyHDFInput3D

    !----------------------------------------------------------------------------

    subroutine ModifyHDFInput3DHarmonics
        
        !Arguments------------------------------------------------------------

        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        type (T_Time)                                   :: Now, CurrentTime

        !Begin----------------------------------------------------------------

        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadNewField - ModuleFillMatrix - ERR10'

        if (Me%BackTracking) then  
            call BacktrackingTime(Now)
        else   
            Now = CurrentTime
        endif 
    
        if (Me%HDF%SpatialInterpolON) then
        
            call ModifyField4DInterpol(CurrentTime      = Now,                          & 
                                       Matrix3D         = Me%Matrix3D)
        
        else
            call ModifyField4D(Field4DID        = Me%HDF%ObjField4D,                    &
                               PropertyIDNumber = Me%PropertyID%IDNumber,               & 
                               CurrentTime      = Now,                                  & 
                               Matrix3D         = Me%Matrix3D,                          &
                               STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR20'  
        endif
            
    end subroutine ModifyHDFInput3DHarmonics
        
    !----------------------------------------------------------------------------

    subroutine ModifyHDFInput3DStandard(PointsToFill3D)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        !Local----------------------------------------------------------------
        integer                                         :: n, i, j, k
        type (T_Time)                                   :: Now

        !Begin----------------------------------------------------------------

i2:     if (ReadNewField(Now,n))then 
    
i4:         if (n==1) then
                call SetMatrixValue(Me%HDF%PreviousField3D, Me%WorkSize3D, Me%HDF%NextField3D)
            else i4
                call ReadHDF5Values3D(Me%HDF%PreviousInstant, Me%HDF%PreviousField3D)
                
                !limit maximum values
                do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
                do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB

#ifndef _NOT_IEEE_ARITHMETIC
                    if (ieee_is_nan (Me%HDF%PreviousField3D(i,j,k)))                    &
                        Me%HDF%PreviousField3D       (i,j,k) = FillValueReal 
#endif
                
                    if (abs(Me%HDF%PreviousField3D(i,j,k)) > abs(FillValueReal))        &
                            Me%HDF%PreviousField3D(i,j,k) = FillValueReal

                enddo
                enddo
                enddo                
                
            endif i4

            call ReadHDF5Values3D(Me%HDF%NextInstant, Me%HDF%NextField3D)
            
            !limit maximum values
            do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
            
#ifndef _NOT_IEEE_ARITHMETIC
                if (ieee_is_nan (Me%HDF%NextField3D (i,j,k)))                           &
                    Me%HDF%NextField3D        (i,j,k) = FillValueReal 
#endif

                if (abs(Me%HDF%NextField3D    (i,j,k)) > abs(FillValueReal))            &
                        Me%HDF%NextField3D    (i,j,k) = FillValueReal

            enddo
            enddo
            enddo                
            
        end if i2

i3:     if (Me%HDF%PreviousInstant /= Me%HDF%NextInstant) then
                                            
i5:         if (Me%PreviousInstantValues) then
            
                Me%Matrix3D = Me%HDF%PreviousField3D

            else i5

                call InterpolateMatrix3DInTime(ActualTime       = Now,                         &
                                               Size             = Me%WorkSize3D,               &
                                               Time1            = Me%HDF%PreviousTime,         &
                                               Matrix1          = Me%HDF%PreviousField3D,      &
                                               Time2            = Me%HDF%NextTime,             &
                                               Matrix2          = Me%HDF%NextField3D,          &
                                               MatrixOut        = Me%Matrix3D,                 &
                                               PointsToFill3D   = PointsToFill3D)
            endif i5
        else i3
            
            !Prev and next are equal (last instant?)
            call SetMatrixValue(Me%Matrix3D, Me%WorkSize3D, Me%HDF%NextField3D)

        endif i3



    end subroutine ModifyHDFInput3DStandard
            
    !----------------------------------------------------------------------------

    subroutine ModifyHDFInput3DTime(PointsToFill3D)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        !Local----------------------------------------------------------------

        !Begin----------------------------------------------------------------

i1:     if (Me%HDF%Field4D .and. Me%HDF%HarmonicsOn) then

            call ModifyHDFInput3DHarmonics

        else i1 
        
            call ModifyHDFInput3DStandard(PointsToFill3D)
        
        endif i1

    end subroutine ModifyHDFInput3DTime


    !--------------------------------------------------------------------------


    subroutine ModifyHDFInput3DGeneric4D(PointsToFill3D, Generic_4D_Value_)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        real                                            :: Generic_4D_Value_

        !Local----------------------------------------------------------------
        integer                                         :: PrevI, NextI, i, j, k
        !Begin----------------------------------------------------------------
        
i1:     if (.not.(Me%HDF%Previous4DValue <= Generic_4D_Value_ .and.                     &
                  Me%HDF%Next4DValue     >= Generic_4D_Value_)) then
            !Found new limits
            PrevI              = 1
            NextI              = 2
            Me%HDF%Next4DValue = HDF5Generic4DInstant(1)
            do 

                Me%HDF%Previous4DValue  = Me%HDF%Next4DValue
                Me%HDF%Next4DValue      = HDF5Generic4DInstant(NextI)

                if (Me%HDF%Previous4DValue <= Generic_4D_Value_ .and.                  &
                    Me%HDF%Next4DValue     >= Generic_4D_Value_) then
                    exit
                endif

                if (NextI > Me%HDF%NumberOfInstants) then

                    write(*,*)
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'Time instants inconsistency.'
                    stop      'ModifyHDFInput3DGeneric4D - ModuleFillMatrix - ERR10'

                endif

                PrevI = NextI
                NextI = NextI + 1

            enddo

            Me%HDF%NextInstant     = NextI
            Me%HDF%PreviousInstant = PrevI  


            !call SetMatrixValue(Me%HDF%PreviousField3D, Me%WorkSize3D, Me%HDF%NextField3D)

            call ReadHDF5Values3D(Me%HDF%PreviousInstant, Me%HDF%PreviousField3D)

            call ReadHDF5Values3D(Me%HDF%NextInstant,     Me%HDF%NextField3D)
            
            !limit maximum values
            do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB

#ifndef _NOT_IEEE_ARITHMETIC
                if (ieee_is_nan (Me%HDF%PreviousField3D(i,j,k)))                        &
                    Me%HDF%PreviousField3D       (i,j,k) = FillValueReal 
#endif
            
                if (abs(Me%HDF%PreviousField3D(i,j,k)) > abs(FillValueReal))            &
                        Me%HDF%PreviousField3D(i,j,k) = FillValueReal
      
#ifndef _NOT_IEEE_ARITHMETIC                  
                if (ieee_is_nan (Me%HDF%NextField3D    (i,j,k)))                        &
                    Me%HDF%NextField3D           (i,j,k) = FillValueReal 
#endif                        
                
                if (abs(Me%HDF%NextField3D    (i,j,k)) > abs(FillValueReal))            &
                        Me%HDF%NextField3D    (i,j,k) = FillValueReal
            enddo
            enddo
            enddo                
            
                    
        endif i1
        

        if (Me%HDF%PreviousInstant /= Me%HDF%NextInstant) then

            call InterpolateLinearyMatrix3D(X                = Generic_4D_Value_,       &
                                            Size             = Me%WorkSize3D,           &
                                            X1               = Me%HDF%Previous4DValue,  &
                                            Matrix1          = Me%HDF%PreviousField3D,  &
                                            X2               = Me%HDF%Next4DValue,      &
                                            Matrix2          = Me%HDF%NextField3D,      &
                                            MatrixOut        = Me%Matrix3D,             &
                                            PointsToFill3D   = PointsToFill3D)

        else
            
            !Prev and next are equal (last instant?)
            call SetMatrixValue(Me%Matrix3D, Me%WorkSize3D, Me%HDF%NextField3D)

        endif


    end subroutine ModifyHDFInput3DGeneric4D


    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------


    subroutine ModifyHDFInput2D(PointsToFill2D, Generic_4D_Value_)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        real, optional                                  :: Generic_4D_Value_

        !Local----------------------------------------------------------------

        !Begin----------------------------------------------------------------

        if (Me%HDF%Generic4D%ON) then

            if (.not. present(Generic_4D_Value_)) &
                stop 'ModifyHDFInput2D - ModuleFillMAtrix - ERR010'
            call ModifyHDFInput2DGeneric4D(PointsToFill2D, Generic_4D_Value_)
            
        else

            call ModifyHDFInput2DTime(PointsToFill2D)

        endif

    end subroutine ModifyHDFInput2D

    !----------------------------------------------------------------------------

    subroutine ModifyHDFInput2DHarmonics()
        
        !Arguments------------------------------------------------------------
        !Local----------------------------------------------------------------
        type (T_Time)                                   :: Now, CurrentTime
        integer                                         :: STAT_CALL

        !Begin----------------------------------------------------------------

        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyHDFInput2DHarmonics - ModuleFillMatrix - ERR10'

        if (Me%BackTracking) then  
            call BacktrackingTime(Now)
        else   
            Now = CurrentTime
        endif 
        
        if (Me%HDF%SpatialInterpolON) then
        
            call ModifyField4DInterpol(CurrentTime      = Now,                          & 
                                       Matrix2D         = Me%Matrix2D)
        
        else
        

            call ModifyField4D(Field4DID        = Me%HDF%ObjField4D,                    &
                               PropertyIDNumber = Me%PropertyID%IDNumber,               & 
                               CurrentTime      = Now,                                  & 
                               Matrix2D         = Me%Matrix2D,                          &
                               STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ModifyHDFInput2DHarmonics - ModuleFillMatrix - ERR20'  

        endif

    end subroutine ModifyHDFInput2DHarmonics
        
    !----------------------------------------------------------------------------        


    subroutine ModifyHDFInput2DRainType(PointsToFill2D)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        !Local----------------------------------------------------------------
        type (T_Time)                                   :: Now
        integer                                         :: STAT_CALL

        !Begin----------------------------------------------------------------

        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyHDFInput2DRainType - ModuleFillMatrix - ERR010'
        
i3:             if (Now > Me%HDF%NextTime) then
            call ActualizeHDFTimes (Now)
            call ActualizeHDFValues (Me%HDF%PreviousInstant, Me%HDF%PreviousField2D)
            call ActualizeHDFValues (Me%HDF%NextInstant, Me%HDF%NextField2D)
        endif i3
        
        !avoid evaluate in cosntruct phase where previous and next time are the same
i4:     if (Now > Me%HDF%PreviousTime) then
i5:         if (Me%UseOriginalValues) then
                Me%Matrix2D = Me%HDF%NextField2D                    
            else if (Me%AccumulateValues) then i5
                Me%Matrix2D = Me%HDF%NextField2D / (Me%HDF%NextTime - Me%HDF%PreviousTime)
            else i5
                !Interpolates the two matrixes in time
                call InterpolateMatrix2DInTime(ActualTime       = Now,                         &
                                               Size             = Me%WorkSize2D,               &
                                               Time1            = Me%HDF%PreviousTime,         &
                                               Matrix1          = Me%HDF%PreviousField2D,      &
                                               Time2            = Me%HDF%NextTime,             &
                                               Matrix2          = Me%HDF%NextField2D,          &
                                               MatrixOut        = Me%Matrix2D,                 &
                                               PointsToFill2D   = PointsToFill2D)                
            endif i5
        endif i4          
        
i6:     if (Me%ValueIsUsedForDTPrediction) then
i7:         if (Now >= Me%NextEventEnd) then                             
                call FindNextEventInHDF (Now)
i8:             if (Me%AccumulateValues .and. (Me%NextValueForDTPred > 0.0)) then
                    Me%NextValueForDTPred = Me%NextValueForDTPred / (Me%NextEventEnd - Me%NextEventStart)
                endif i8
            endif i7
            
i9:         if (Now >= Me%NextEventStart .and. Now < Me%NextEventEnd) then
                Me%DTForNextEvent = 0.0
            else i9
                Me%DTForNextEvent = Me%NextEventStart - Now 
            endif i9
            
i10:        if (Me%DTForNextEvent > 0.0) then
                Me%PredictedDT = Me%DTForNextEvent                    
            else i10
                Me%PredictedDT = Me%NextEventEnd - Now
            endif i10                
        endif  i6


    end subroutine ModifyHDFInput2DRainType
        
    !----------------------------------------------------------------------------    


    subroutine ModifyHDFInput2DStandard(PointsToFill2D)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        !Local----------------------------------------------------------------
        integer                                         :: n, i, j
        type (T_Time)                                   :: Now

        !Begin----------------------------------------------------------------

i11:    if (ReadNewField(Now,n))then
i12:        if (n==1) then 
                call SetMatrixValue(Me%HDF%PreviousField2D, Me%WorkSize2D, Me%HDF%NextField2D, PointsToFill2D)
            else i12
                call ReadHDF5Values2D(Me%HDF%PreviousInstant, Me%HDF%PreviousField2D)
                
                !limit maximum values
                do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
                do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
            
#ifndef _NOT_IEEE_ARITHMETIC                    
                    if (ieee_is_nan (Me%HDF%PreviousField2D(i,j)))                      &
                        Me%HDF%PreviousField2D(i,j) = FillValueReal                
#endif
                    
                    if (abs(Me%HDF%PreviousField2D(i,j)) > abs(FillValueReal))          &
                        Me%HDF%PreviousField2D(i,j) = FillValueReal
                enddo
                enddo   
                    
            endif i12

            call ReadHDF5Values2D(Me%HDF%NextInstant, Me%HDF%NextField2D)
            
            !limit maximum values
            do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
            
#ifndef _NOT_IEEE_ARITHMETIC
                if (ieee_is_nan (Me%HDF%PreviousField2D(i,j)))                      &
                    Me%HDF%PreviousField2D(i,j) = FillValueReal                
#endif

                if (abs(Me%HDF%PreviousField2D(i,j)) > abs(FillValueReal))          &
                    Me%HDF%PreviousField2D(i,j) = FillValueReal
            enddo
            enddo   
                
        end if i11

i15:    if (Me%HDF%PreviousInstant /= Me%HDF%NextInstant) then

i16:        if (Me%UseOriginalValues) then
                
                Me%Matrix2D = Me%HDF%NextField2D
                
i17:            if (Me%PredictDTMethod == 1) then
                    call PredictDTForHDF (PointsToFill2D, Me%HDF%PreviousTime, Me%HDF%NextTime, Now)
                elseif (Me%PredictDTMethod == 2) then i17
                    call PredictDTForHDF_New (PointsToFill2D, Now)
                else i17
                    stop 'ModifyHDFInput2DStandard - ModuleFillMatrix - ERR010'
                endif i17
            
            else if (Me%AccumulateValues) then  i16     !For Rain
                    
                Me%Matrix2D = Me%HDF%NextField2D / (Me%HDF%NextTime - Me%HDF%PreviousTime)
    
                !This will replace the processed values if the value in NextField2D is a NODATA value
i19:            if (.not. Me%IgnoreNoDataPoint) then
                    do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
                    do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
i18:                            if (Me%HDF%NextField2D (i,j) == Me%NoDataValue) then
                            Me%Matrix2D (i,j) = Me%NoDataValue
                        endif i18
                    enddo
                    enddo
                endif i19
    
i20:            if (Me%PredictDTMethod == 1) then
                    call PredictDTForHDF (PointsToFill2D, Me%HDF%PreviousTime, Me%HDF%NextTime, Now)
                elseif (Me%PredictDTMethod == 2) then i20
                    call PredictDTForHDF_New (PointsToFill2D, Now)
                else  i20
                    stop 'ModifyHDFInput2DStandard - ModuleFillMatrix - ERR020'
                endif i20
                
            else i16
                !Interpolates the two matrixes in time
                call InterpolateMatrix2DInTime(ActualTime       = Now,                         &
                                               Size             = Me%WorkSize2D,               &
                                               Time1            = Me%HDF%PreviousTime,         &
                                               Matrix1          = Me%HDF%PreviousField2D,      &
                                               Time2            = Me%HDF%NextTime,             &
                                               Matrix2          = Me%HDF%NextField2D,          &
                                               MatrixOut        = Me%Matrix2D,                 &
                                               PointsToFill2D   = PointsToFill2D)
                                               
                !This will replace the processed values if the value in NextField2D is a NODATA value
i21:            if (.not. Me%IgnoreNoDataPoint) then
                    do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
                    do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
i22:                            if (Me%HDF%NextField2D (i,j) == Me%NoDataValue) then
                            Me%Matrix2D (i,j) = Me%NoDataValue
                        endif i22
                    enddo
                    enddo
                endif i21                                               
                
            endif i16
                
        else i15

            !Prev and next are equal (last instant?)
i23:        if (Me%UseOriginalValues .or. Me%InterpolateValues) then

                call SetMatrixValue(Me%Matrix2D, Me%WorkSize2D, Me%HDF%PreviousField2D, PointsToFill2D)

            else i23
            
                !do nothing
                
            endif i23

        endif i15


    end subroutine ModifyHDFInput2DStandard
        
    !----------------------------------------------------------------------------    
    
    subroutine ModifyHDFInput2DTime(PointsToFill2D)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        !Local----------------------------------------------------------------


        !Begin----------------------------------------------------------------
        
i1:     if (Me%HDF%Field4D .and. Me%HDF%HarmonicsOn) then

            call ModifyHDFInput2DHarmonics
            
        else i1 
        
i2:         if (Me%PredictDTMethod == 2) then

                call ModifyHDFInput2DRainType(PointsToFill2D)            
                      
            else i2
            
                call ModifyHDFInput2DStandard(PointsToFill2D)
                
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
    
    logical function ReadNewField(Now,n)

        !Arguments------------------------------------------------------------
        type (T_Time), intent(OUT)                      :: Now
        integer      , intent(OUT)                      :: n
        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        type (T_Time)                                   :: CurrentTime
        logical                                         :: ReadNewField_

        !Begin----------------------------------------------------------------

        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadNewField - ModuleFillMatrix - ERR10'

        if (Me%BackTracking) then  
            call BacktrackingTime(Now)
        else   
            Now = CurrentTime
        endif     

        ReadNewField_ = .false.
        
        if (.not. Me%AccumulateValues) then
            if (Me%BackTracking) then  
                if (Now .le. Me%HDF%NextTime) ReadNewField_ = .true.
            else            
                if (Now .ge. Me%HDF%NextTime) ReadNewField_ = .true.
            endif
        else
            if (Me%BackTracking) then  
                if (Now .lt. Me%HDF%NextTime) ReadNewField_ = .true.
            else            
                if (Now .gt. Me%HDF%NextTime) ReadNewField_ = .true.
            endif        
        endif
        
        ReadNewField = ReadNewField_

        if (ReadNewField_)then

            n = 0
            
            do 

                if (Me%BackTracking) then  
                    if (Now .gt. Me%HDF%NextTime) exit
                else            
                    if (Now .lt. Me%HDF%NextTime) exit
                endif            
                
                Me%HDF%PreviousInstant  = Me%HDF%NextInstant
                    
                if (Me%BackTracking) then
                    if(Me%HDF%NextInstant .gt. 1)then
                        Me%HDF%NextInstant  = Me%HDF%NextInstant - 1
                    else
                        exit
                    endif
                else
                    if(Me%HDF%NextInstant .lt. Me%HDF%NumberOfInstants)then
                        Me%HDF%NextInstant  = Me%HDF%NextInstant + 1
                    else
                        if (Me%HDF%GenericYear) then 
                            Me%HDF%NextInstant  = 1
                        else
                            exit
                        endif
                    endif
                endif
                
                Me%HDF%PreviousTime     = Me%HDF%NextTime
                Me%HDF%NextTime         = HDF5TimeInstant(Me%HDF%NextInstant)
                
                if (Me%HDF%GenericYear) then
                    if(Me%HDF%NextInstant > 1)then
                        call SetHDFGenericYear(Me%HDF%NextTime, Now)
                    else
                        call SetHDFGenericYear(Me%HDF%NextTime, Now, AddYear = .true.)
                    endif
                endif

                n = n + 1
                
                
            enddo
            
            if (Me%BackTracking) then
                if(Now .lt. Me%HDF%NextTime)then
                    write(*,*)
                    write(*,*)'----------Backtracking mode-----------'
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'Time instants inconsistency.'
                    stop      'ReadNewField - ModuleFillMatrix - ERR20'
                end if
            else
                if(Now .gt. Me%HDF%NextTime)then
                    write(*,*)
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'Time instants inconsistency.'
                    stop      'ReadNewField - ModuleFillMatrix - ERR30'
                end if
            endif    
    
        endif 
        
    end function ReadNewField

    !-------------------------------------------------------------------------

    subroutine ModifyHDFInput2DGeneric4D(PointsToFill2D, Generic_4D_Value_)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        real                                            :: Generic_4D_Value_

        !Local----------------------------------------------------------------
        integer                                         :: PrevI, NextI, i, j 
        !Begin----------------------------------------------------------------
        
i1:     if (.not.(Me%HDF%Previous4DValue <= Generic_4D_Value_ .and.                     &
                  Me%HDF%Next4DValue     >= Generic_4D_Value_)) then
            !Found new limits
            PrevI              = 1
            NextI              = 2
            Me%HDF%Next4DValue = HDF5Generic4DInstant(1)
            do 

                Me%HDF%Previous4DValue  = Me%HDF%Next4DValue
                Me%HDF%Next4DValue      = HDF5Generic4DInstant(NextI)

                if (Me%HDF%Previous4DValue <= Generic_4D_Value_ .and.                  &
                    Me%HDF%Next4DValue     >= Generic_4D_Value_) then
                    exit
                endif

                if (NextI > Me%HDF%NumberOfInstants) then

                    write(*,*)
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'Time instants inconsistency.'
                    stop      'ModifyHDFInput2DGeneric4D - ModuleFillMatrix - ERR10'

                endif

                PrevI = NextI
                NextI = NextI + 1

            enddo

            Me%HDF%NextInstant     = NextI
            Me%HDF%PreviousInstant = PrevI  


            !call SetMatrixValue(Me%HDF%PreviousField2D, Me%WorkSize2D, Me%HDF%NextField2D)

            call ReadHDF5Values2D(Me%HDF%PreviousInstant, Me%HDF%PreviousField2D)
            call ReadHDF5Values2D(Me%HDF%NextInstant,     Me%HDF%NextField2D)
            
            !limit maximum values
            do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB

#ifndef _NOT_IEEE_ARITHMETIC
                if (ieee_is_nan (Me%HDF%PreviousField2D(i,j)))                          &
                    Me%HDF%PreviousField2D(i,j) = FillValueReal                    
#endif
            
                if (abs(Me%HDF%PreviousField2D(i,j)) > abs(FillValueReal))              &
                    Me%HDF%PreviousField2D(i,j) = FillValueReal

#ifndef _NOT_IEEE_ARITHMETIC
                if (ieee_is_nan (Me%HDF%NextField2D(i,j)))                              &
                    Me%HDF%NextField2D    (i,j) = FillValueReal                    
#endif
                
                if (abs(Me%HDF%NextField2D    (i,j)) > abs(FillValueReal))              &
                    Me%HDF%NextField2D    (i,j) = FillValueReal
            enddo
            enddo   
            
                    
        endif i1
        

        if (Me%HDF%PreviousInstant /= Me%HDF%NextInstant) then

            call InterpolateLinearyMatrix2D(X                = Generic_4D_Value_,       &
                                            Size             = Me%WorkSize2D,           &
                                            X1               = Me%HDF%Previous4DValue,  &
                                            Matrix1          = Me%HDF%PreviousField2D,  &
                                            X2               = Me%HDF%Next4DValue,      &
                                            Matrix2          = Me%HDF%NextField2D,      &
                                            MatrixOut        = Me%Matrix2D,             &
                                            PointsToFill2D   = PointsToFill2D)

        else
          
            !Prev and next are equal (last instant?)
            call SetMatrixValue(Me%Matrix2D, Me%WorkSize2D, Me%HDF%NextField2D)

        endif


    end subroutine ModifyHDFInput2DGeneric4D


    !--------------------------------------------------------------------------

    subroutine PredictDTForHDF(PointsToFill2D, Time1, Time2, ActualTime)
        
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        type (T_Time)                                   :: Time1, Time2, ActualTime

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

            Me%PredictedDT     = min(aux1, aux2)
            Me%DTForNextEvent  = 0.0
        
        else

            !Can run until next Matrix will be read
            !This prediction is different from the one done by the GetTimeSerieDTForNextEvent
            !but I (Frank) assume that the datasets in the HDF always have a quite larger DT
            !then the model will use to run
            Me%PredictedDT     = Time2 - ActualTime
            Me%DTForNextEvent  = Time2 - ActualTime
        endif
                        

    end subroutine PredictDTForHDF

    !--------------------------------------------------------------------------

    subroutine PredictDTForHDF_New(PointsToFill2D, ActualTime)
        
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        type (T_Time)                                   :: ActualTime

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------
        
        if (ActualTime >= Me%HDF%NextTime) then
            call FindDTForNextEventInHDF(PointsToFill2D, ActualTime)             
        else
            if (Me%DTForNextEvent > 0.0) then
                Me%DTForNextEvent = Me%HDF%PreviousTime - ActualTime
                if (Me%DTForNextEvent<0.0) Me%DTForNextEvent = 0.0
            endif
            Me%DTForNextDataset = Me%TimeOfNextDataset - ActualTime
            if (Me%DTForNextDataset <= 0.0) then
                Me%TimeOfNextDataset = TimeOfNextDatasetInHDF(ActualTime) 
                Me%DTForNextDataset  = Me%TimeOfNextDataset - ActualTime
            endif
        endif
        
        !----------------------------------------------------------------------

    end subroutine PredictDTForHDF_New

    !--------------------------------------------------------------------------
    
    type (T_Time) function TimeOfNextDatasetInHDF (ActualTime)
        !NOT for use with backtracking mode
    
        !Arguments-------------------------------------------------------------        
        Type (T_Time), intent(IN) :: ActualTime
        
        !Local-----------------------------------------------------------------
         
        !----------------------------------------------------------------------
               
        TimeOfNextDatasetInHDF = Me%TimeOfNextDataset
                    
        if (Me%HDF%CyclicTimeON) then 
            TimeOfNextDatasetInHDF = Me%HDF%NextTime        
        else  
do1:        do while (TimeOfNextDatasetInHDF <= ActualTime)              
                if (Me%InstantOfNextDataset < Me%HDF%NumberOfInstants) then
                    Me%InstantOfNextDataset = Me%InstantOfNextDataset + 1
                    TimeOfNextDatasetInHDF = HDF5TimeInstant(Me%InstantOfNextDataset)
                else
                    exit do1
                endif
            enddo do1
        endif
            
        !----------------------------------------------------------------------
    
    end function TimeOfNextDatasetInHDF
        
    !--------------------------------------------------------------------------
    
    subroutine FindDTForNextEventInHDF(PointsToFill2D, ActualTime)
    
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer, intent(IN)   :: PointsToFill2D
        Type (T_Time), intent(IN)                       :: ActualTime

        !Local-----------------------------------------------------------------
        integer                                         :: i, j        
        logical                                         :: DTForNextDatasetWasSet
        type (T_Time)                                   :: PreviousTime
        
        !----------------------------------------------------------------------
        PreviousTime = Me%HDF%NextTime
        DTForNextDatasetWasSet = .false.

doF:    do

            Me%HDF%PreviousInstant  = Me%HDF%NextInstant
            
            if(Me%HDF%NextInstant .lt. Me%HDF%NumberOfInstants)then
                Me%HDF%NextInstant  = Me%HDF%NextInstant + 1                
                Me%HDF%PreviousTime = Me%HDF%NextTime
                Me%HDF%NextTime     = HDF5TimeInstant(Me%HDF%NextInstant)
                
                if (.not. DTForNextDatasetWasSet) then
                    DTForNextDatasetWasSet  = .true.
                    Me%DTForNextDataset     = Me%HDF%NextTime - ActualTime
                    Me%TimeOfNextDataset    = Me%HDF%NextTime
                    Me%InstantOfNextDataset = Me%HDF%NextInstant
                endif
            else
                exit doF
            endif

            call ReadHDF5Values2D(Me%HDF%PreviousInstant, Me%HDF%PreviousField2D)
            do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
#ifndef _NOT_IEEE_ARITHMETIC                    
                if (ieee_is_nan (Me%HDF%PreviousField2D(i,j))) &
                    Me%HDF%PreviousField2D(i,j) = FillValueReal
#endif                    
                if (abs(Me%HDF%PreviousField2D(i,j)) > abs(FillValueReal)) &
                    Me%HDF%PreviousField2D(i,j) = FillValueReal            
            enddo
            enddo   

            call ReadHDF5Values2D(Me%HDF%NextInstant, Me%HDF%NextField2D)
            do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
#ifndef _NOT_IEEE_ARITHMETIC            
                if (ieee_is_nan (Me%HDF%NextField2D(i,j))) &
                    Me%HDF%NextField2D (i,j) = FillValueReal                      
#endif                    
                if (abs(Me%HDF%NextField2D(i,j)) > abs(FillValueReal)) &
                    Me%HDF%NextField2D (i,j) = FillValueReal
            enddo
            enddo  

doM:        do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                if(PointsToFill2D(i,j) == 1)then
                    if (Me%HDF%NextField2D(i,j) > Me%MinForDTDecrease) then
                        exit doF
                    endif
                endif
            enddo
            enddo doM
            
        enddo doF
        
        if (PreviousTime == Me%HDF%PreviousTime) then
            Me%DTForNextEvent = 0.0
        else
            Me%DTForNextEvent = Me%HDF%PreviousTime - ActualTime
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

        call InterpolateMatrix3DInTime(ActualTime       = Now,                                   &
                                       Size             = Me%WorkSize3D,                         &
                                       Time1            = Me%ProfileTimeSerie%PreviousTime,      &
                                       Matrix1          = Me%ProfileTimeSerie%PreviousField3D,   &
                                       Time2            = Me%ProfileTimeSerie%NextTime,          &
                                       Matrix2          = Me%ProfileTimeSerie%NextField3D,       &
                                       MatrixOut        = Me%Matrix3D,                           &
                                       PointsToFill3D   = PointsToFill3D)

    end subroutine ModifyProfileTimeSerie

    !--------------------------------------------------------------------------


    subroutine ModifySpaceTimeSerie (PointsToFill2D, PointsToFill3D)

        !Arguments-------------------------------------------------------------
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
        
        !Begin----------------------------------------------------------------
        
        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR010'

        if (.not. Me%RemainsConstant) then

            if (Me%PredictDTMethod == 1) then
                !Gets Value for current Time
                call GetTimeSerieValue (Me%TimeSerie%ObjTimeSerie, Now,                  &
                                        Me%TimeSerie%Column,                             &
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
    
                    call GetTimeSerieDTForNextEvent (Me%TimeSerie%ObjTimeSerie,          &
                                                     NewValue, Me%TimeSerie%Column, Now, &
                                                     Me%PredictedDT, Me%DTForNextEvent,  &
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

                NewValue = 270. - Value1 + Angle

            else
                
                Value1   = 270. - Value1 + Angle
                Value2   = 270. - Value2 + Angle
                
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
                
                if (Now > Me%Timeserie%NextTime) then
                    call ActualizeTimeSerieTimes (Now)
                    call ActualizeTimeSerieValues                   
        endif
                
                !avoid evaluate in cosntruct phase where previous and next time are the same
                if (Now > Me%Timeserie%PreviousTime) then
                    if (Me%UseOriginalValues) then
                        Me%Timeserie%CurrentValue = Me%Timeserie%NextValue
                    elseif (Me%AccumulateValues) then
                        Me%Timeserie%CurrentValue = Me%Timeserie%NextValue / (Me%Timeserie%NextTime - Me%Timeserie%PreviousTime)
                    elseif (Me%InterpolateValues) then
                        call InterpolateValueInTime(Now,                        &
                                                    Me%Timeserie%PreviousTime,  &
                                                    Me%Timeserie%PreviousValue, &
                                                    Me%Timeserie%NextTime,      &
                                                    Me%Timeserie%NextValue,     &
                                                    Me%Timeserie%CurrentValue)                    
                    else
                        stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR060'
                    endif
                endif
                
                NewValue = Me%Timeserie%CurrentValue
                
                if (Me%ValueIsUsedForDTPrediction) then
                    if (Now >= Me%NextEventEnd) then
                        call FindNextEventInTimeSerie (Now)
                        if (Me%AccumulateValues .and. (Me%NextValueForDTPred > 0.0)) then
                            Me%NextValueForDTPred = Me%NextValueForDTPred / (Me%NextEventEnd - Me%NextEventStart)
                        endif                    
                    endif
                    
                    if (Now >= Me%NextEventStart .and. Now < Me%NextEventEnd) then
                        Me%DTForNextEvent = 0.0
                    else
                        Me%DTForNextEvent = Me%NextEventStart - Now 
                    endif
                    
                    if (Me%DTForNextEvent > 0.0) then
                        Me%PredictedDT = Me%DTForNextEvent                    
                    else
                        Me%PredictedDT = Me%NextEventEnd - Now
                    endif                
                endif
            endif
        else
        
            NewValue = Me%TimeSerie%CurrentValue
            
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
    
    subroutine ActualizeTimeSerieTimes (ActualTime)
    
        !Arguments-------------------------------------------------------------        
        Type (T_Time), intent(IN)                       :: ActualTime

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        
        !----------------------------------------------------------------------                     
        
        do while (Me%TimeSerie%NextTime < ActualTime)  
                          
            if (Me%TimeSerie%NextInstant < Me%TimeSerie%NumberOfInstants) then
            
                Me%TimeSerie%PreviousInstant = Me%TimeSerie%NextInstant
                Me%TimeSerie%PreviousTime    = Me%TimeSerie%NextTime
                Me%TimeSerie%NextInstant     = Me%TimeSerie%NextInstant + 1
                
            else
            
                stop 'ActualizeTimeSerieTimes - ModuleFillMatrix - ERR010'
                
            endif
                  
            call GetTimeSerieTimeOfDataset(Me%TimeSerie%ObjTimeSerie,    &
                                           Me%TimeSerie%NextInstant,     &
                                           Me%TimeSerie%NextTime,        &
                                           STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)    &
                stop 'ActualizeTimeSerieTimes - ModuleFillMatrix - ERR020'  
           
        enddo
 
        !----------------------------------------------------------------------
        
    end subroutine ActualizeTimeSerieTimes 
    
    !--------------------------------------------------------------------------
    
    subroutine ActualizeTimeSerieValues
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_
        
        !----------------------------------------------------------------------    
    
        call GetTimeSerieValueForIndex (Me%TimeSerie%ObjTimeSerie,      &
                                        Me%TimeSerie%PreviousInstant,   &
                                        Me%TimeSerie%Column,            &
                                        Me%TimeSerie%PreviousValue,     &
                                        STAT = STAT_)
        if (STAT_ /= SUCCESS_) stop 'ActualizeTimeSerieValues - ModuleFillMatrix - ERR010'            

        call GetTimeSerieValueForIndex (Me%TimeSerie%ObjTimeSerie,      &
                                        Me%TimeSerie%NextInstant,       &
                                        Me%TimeSerie%Column,            &
                                        Me%TimeSerie%NextValue,         &                                            
                                        STAT = STAT_)
        if (STAT_ /= SUCCESS_) stop 'ActualizeTimeSerieValues - ModuleFillMatrix - ERR020' 
        
        !----------------------------------------------------------------------    
    
    end subroutine ActualizeTimeSerieValues
    
    !--------------------------------------------------------------------------
    
    subroutine FindNextEventInTimeSerie(Now)
    
        !Arguments-------------------------------------------------------------        
        Type (T_Time), intent(IN)                       :: Now

        !Local-----------------------------------------------------------------
        real                                            :: instant_value
        integer                                         :: instant, STAT_CALL

        !----------------------------------------------------------------------
        
        if (Me%TimeSerie%NextInstant < Me%TimeSerie%NumberOfInstants) then
        
            if (Now > Me%NextEventEnd) then
        
                Me%NextEventStart = Me%TimeSerie%PreviousTime
                instant           = Me%TimeSerie%NextInstant
                Me%NextEventEnd   = Me%TimeSerie%NextTime                
                instant_value     = Me%TimeSerie%NextValue

            else
            
                Me%NextEventStart = Me%TimeSerie%NextTime
                instant           = Me%TimeSerie%NextInstant + 1
                call GetTimeSerieTimeOfDataset(Me%TimeSerie%ObjTimeSerie,   &
                                               instant,                     &
                                               Me%NextEventEnd,             &
                                               STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'FindNextEventInTimeSerie - ModuleFillMatrix - ERR010'                 
                call GetTimeSerieValueForIndex (Me%TimeSerie%ObjTimeSerie,  &
                                                instant,                    &
                                                Me%TimeSerie%Column,        &
                                                instant_value,              &                                            
                                                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'FindNextEventInTimeSerie - ModuleFillMatrix - ERR020'   
            
            endif
                             
            do while (instant_value <= Me%MinForDTDecrease .and. instant < Me%TimeSerie%NumberOfInstants)                
                instant = instant + 1  
                    
                call GetTimeSerieValueForIndex (Me%TimeSerie%ObjTimeSerie,  &
                                                instant,                    &
                                                Me%TimeSerie%Column,        &
                                                instant_value,              &                                            
                                                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'FindNextEventInTimeSerie - ModuleFillMatrix - ERR030'                                
            enddo

            if (instant_value > Me%MinForDTDecrease) then
            
                call GetTimeSerieTimeOfDataset(Me%TimeSerie%ObjTimeSerie,   &
                                               Instant - 1,                 &
                                               Me%NextEventStart,           &
                                               STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'FindNextEventInTimeSerie - ModuleFillMatrix - ERR040' 
                
                call GetTimeSerieTimeOfDataset(Me%TimeSerie%ObjTimeSerie,   &
                                               instant,                     &
                                               Me%NextEventEnd,             &
                                               STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'FindNextEventInTimeSerie - ModuleFillMatrix - ERR050'              
                    
            else
            
                call GetTimeSerieTimeOfDataset(Me%TimeSerie%ObjTimeSerie,   &
                                               instant,                     &
                                               Me%NextEventStart,           &
                                               STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'FindNextEventInTimeSerie - ModuleFillMatrix - ERR060'                                 
                        
                Me%NextEventEnd = Me%NextEventStart
                instant_value   = 0.0
                                    
            endif 
            
        else
        
            Me%NextEventStart = Me%TimeSerie%NextTime
            Me%NextEventEnd   = Me%TimeSerie%NextTime                

        endif
        
        Me%NextValueForDTPred = instant_value
        !write (*,*) 'instant_value = ', instant_value

        !----------------------------------------------------------------------
    
    end subroutine FindNextEventInTimeSerie
    
    !--------------------------------------------------------------------------
    
    subroutine ActualizeHDFTimes (ActualTime)
    
        !Arguments-------------------------------------------------------------        
        Type (T_Time), intent(IN)                       :: ActualTime

        !Local-----------------------------------------------------------------
        
        !----------------------------------------------------------------------                     
        
        do while (Me%HDF%NextTime < ActualTime)  
                          
            if (Me%HDF%NextInstant < Me%HDF%NumberOfInstants) then
            
                Me%HDF%PreviousInstant = Me%HDF%NextInstant
                Me%HDF%PreviousTime    = Me%HDF%NextTime
                Me%HDF%NextInstant     = Me%HDF%NextInstant + 1
                
            else
            
                stop 'ActualizeHDFTimes - ModuleFillMatrix - ERR010'
                
            endif
                  
            Me%HDF%NextTime = HDF5TimeInstant(Me%HDF%NextInstant)
           
        enddo
 
        !----------------------------------------------------------------------
        
    end subroutine ActualizeHDFTimes 
    
    !--------------------------------------------------------------------------    

    subroutine ActualizeHDFValues (instant, field2D, field3D)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: instant
        real, dimension(:, :), pointer, optional    :: field2D
        real, dimension(:, :, :), pointer, optional :: field3D

        !Local-----------------------------------------------------------------
        integer :: i, j, k
        
        !----------------------------------------------------------------------    
    
        if (present(field2D)) then
    
            call ReadHDF5Values2D(instant, field2D)        
        
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
        
            call ReadHDF5Values3D(instant, field3D)        
        
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

    subroutine FindNextEventInHDF(Now)
    
        !Arguments-------------------------------------------------------------        
        Type (T_Time), intent(IN)                       :: Now

        !Local-----------------------------------------------------------------
        real                                            :: instant_value
        integer                                         :: instant

        !----------------------------------------------------------------------
        
        if (Me%HDF%NextInstant < Me%HDF%NumberOfInstants) then

            if (Now > Me%NextEventEnd) then
                Me%NextEventStart = Me%HDF%PreviousTime
                instant           = Me%HDF%NextInstant
                Me%NextEventEnd   = Me%HDF%NextTime   
                if (Me%Dim == Dim2D) then             
                    Me%HDF%Array2D    = Me%HDF%NextField2D
                else
                    Me%HDF%Array3D    = Me%HDF%NextField3D
                endif
            else
            
                Me%NextEventStart = Me%HDF%NextTime
                instant           = Me%HDF%NextInstant + 1
                Me%NextEventEnd   = HDF5TimeInstant(instant)
                 
                if (Me%Dim == Dim2D) then
                    call ActualizeHDFValues (instant, field2D = Me%HDF%Array2D)
                else
                    call ActualizeHDFValues (instant, field3D = Me%HDF%Array3D)
                endif
            endif
                 
            if (Me%Dim == Dim2D) then     
                instant_value = maxval(Me%HDF%Array2D)
            else
                instant_value = maxval(Me%HDF%Array3D)
            endif
                             
            do while (instant_value <= Me%MinForDTDecrease .and. instant < Me%TimeSerie%NumberOfInstants)                
                instant = instant + 1  
                    
                if (Me%Dim == Dim2D) then
                    call ActualizeHDFValues (instant, field2D = Me%HDF%Array2D)
                    instant_value = maxval(Me%HDF%Array2D)
                else
                    call ActualizeHDFValues (instant, field3D = Me%HDF%Array3D)
                    instant_value = maxval(Me%HDF%Array3D)
                endif                               
            enddo

            if (instant_value > Me%MinForDTDecrease) then
            
                Me%NextEventStart   = HDF5TimeInstant(instant - 1)
                Me%NextEventEnd     = HDF5TimeInstant(instant)
                    
            else
            
                Me%NextEventStart   = HDF5TimeInstant(instant)                        
                Me%NextEventEnd     = Me%NextEventStart
                instant_value       = 0.0
                                    
            endif 
            
        else
        
            Me%NextEventStart = Me%HDF%NextTime
            Me%NextEventEnd   = Me%HDF%NextTime                

        endif
        
        Me%NextValueForDTPred = instant_value        

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

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mFILLMATRIX_,  Me%InstanceID)

            if (nUsers == 0) then

                !Kills Time Serie
                if (Me%TimeSerie%ObjTimeSerie /= 0) then
                    call KillTimeSerie (Me%TimeSerie%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillFillMatrix - ModuleFillMatrix - ERR10'
                endif

                if (Me%HDF%ObjHDF5 /= 0 .or. Me%HDF%ObjField4D /= 0 ) then

                    if( associated(Me%HDF%PreviousField2D))then
                        deallocate(Me%HDF%PreviousField2D)
                        nullify   (Me%HDF%PreviousField2D)
                    end if

                    if( associated(Me%HDF%NextField2D))then
                        deallocate(Me%HDF%NextField2D)
                        nullify   (Me%HDF%NextField2D)
                    end if

                    if( associated(Me%HDF%PreviousField3D))then
                        deallocate(Me%HDF%PreviousField3D)
                        nullify   (Me%HDF%PreviousField3D)
                    end if

                    if(associated (Me%HDF%NextField3D))then
                        deallocate(Me%HDF%NextField3D)
                        nullify   (Me%HDF%NextField3D)
                    end if
                    
                    if (associated(Me%HDF%ReadField3D)) then
                        deallocate(Me%HDF%ReadField3D, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)stop 'KillFillMatrix - ModuleFillMatrix - ERR20'
                        nullify   (Me%HDF%ReadField3D)
                    endif
                    
if4D:               if (Me%HDF%Field4D) then

                        if (Me%HDF%SpatialInterpolON) then
    
                            if (Me%Dim == Dim3D) then
                                 deallocate(Me%HDF%Z)
                            endif                    
                                
                            deallocate(Me%HDF%X, Me%HDF%Y, Me%HDF%Prop, Me%HDF%NoData)
                        
                        endif
                        
                        call KillField4D(Me%HDF%ObjField4D, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillFillMatrix - ModuleFillMatrix - ERR30'
                      
                    else if4D

                        call KillHDF5(Me%HDF%ObjHDF5, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillFillMatrix - ModuleFillMatrix - ERR40'

                    endif if4D
                      
                    if (Me%HDF%Generic4D%ReadFromTimeSerie) then
                        call KillTimeSerie(Me%HDF%Generic4D%ObjTimeSerie, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillFillMatrix - ModuleFillMatrix - ERR50'
                    endif
  
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
        

    !------------------------------------------------------------------------
    
    
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
