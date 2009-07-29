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

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions,        only : InterpolateValueInTime, InterpolateProfile,       &
                                       SetMatrixValue, InterpolateMatrix2DInTime,        &
                                       InterpolateMatrix3DInTime,                        &
                                       InterpolateLinearyMatrix2D, InterpolateLinearyMatrix3D
    use ModuleBoxDif,           only : StartBoxDif, GetBoxes, GetNumberOfBoxes,          &
                                       UngetBoxDif, KillBoxDif
    use ModuleGridData,         only : ConstructGridData, GetGridData, UnGetGridData,    &
                                       KillGridData, GetGridDataType 
    use ModuleHorizontalGrid,   only : GetGridAngle
    use ModuleTimeSerie,        only : StartTimeSerieInput, GetTimeSerieValue,           &
                                       GetTimeSerieDTForNextEvent, KillTimeSerie
    use ModuleGeometry,         only : GetGeometryDistances, UnGetGeometry
    use ModuleHDF5,             only : ConstructHDF5, HDF5ReadData, GetHDF5GroupID,      &
                                       GetHDF5FileAccess, GetHDF5GroupNumberOfItems,     &
                                       HDF5SetLimits, KillHDF5

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

    !Variable from file
    integer, parameter                              :: None             = 1

    !Type of analytic profiles 
    integer, parameter                              :: Linear           = 1
    integer, parameter                              :: Exponential      = 2


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

    !Types---------------------------------------------------------------------

    type T_Layers
        real, dimension(:), pointer                 :: Values
    end type T_Layers

    type T_Boxes
        character(PathLength)                       :: FileName
        integer                                     :: ObjBoxDif
        real, dimension(:), pointer                 :: Values
    end type T_Boxes

    type T_ASCIIFile
        character(PathLength)                       :: FileName
        integer                                     :: GridDataID
    end type T_ASCIIFile

    type T_TimeSerie
        character(PathLength)                       :: FileName
        integer                                     :: ObjTimeSerie = 0
        integer                                     :: Column
    end type T_TimeSerie

    type T_ProfileTimeSerie
        character(PathLength)                       :: FileName
        type (T_Time)                               :: NextTime,  PreviousTime
        integer                                     :: NextInstant, PreviousInstant 
        real,           dimension(:,:,:), pointer   :: PreviousField3D, NextField3D
        real,           dimension(:,:  ), pointer   :: Values, Depths
        type(T_Time),   dimension(:    ), pointer   :: TimeInstants
        integer                                     :: NumberOfInstants, nValues, nDepths
        logical                                     :: CyclicTimeON = .false.
    end type T_ProfileTimeSerie


    !Generic 4D
    type T_Generic4D
        logical                            :: ON
        logical                            :: ReadFromTimeSerie
        integer                            :: ObjTimeSerie
        integer                            :: TimeSerieColumn
        real                               :: CurrentValue

    end type 

    type T_HDF
        character(PathLength)                       :: FileName, VGroupPath, FieldName
        real                                        :: MultiplyingFactor
        logical                                     :: HasMultiplyingFactor = .false.
        real                                        :: AddingFactor
        logical                                     :: HasAddingFactor = .false.
        type (T_Time)                               :: NextTime,  PreviousTime
        real                                        :: Next4DValue     = FillValueReal
        real                                        :: Previous4DValue = FillValueReal
        integer                                     :: NextInstant, PreviousInstant 
        real, dimension(:,:  ), pointer             :: PreviousField2D, NextField2D
        real, dimension(:,:,:), pointer             :: PreviousField3D, NextField3D, ReadField3D
        integer                                     :: ObjHDF5 = 0
        integer                                     :: NumberOfInstants
        logical                                     :: CyclicTimeON = .false.
        logical                                     :: From2Dto3D   = .false.
        type(T_Generic4D)                           :: Generic4D
    end type T_HDF


    private :: T_FillMatrix
    type       T_FillMatrix
        integer                                     :: InstanceID
        type (T_Size2D)                             :: Size2D, WorkSize2D
        type (T_Size3D)                             :: Size3D, WorkSize3D
        type (T_PropertyID)                         :: PropertyID
        integer                                     :: Dim                  !2D/3D
        integer                                     :: TypeZUV              !Z/U/V
        integer                                     :: TimeEvolution
        integer                                     :: SpaceEvolution
        integer                                     :: InitializationMethod
        integer                                     :: InitializationDefault
        logical                                     :: RemainsConstant      = .false.
        logical                                     :: AccumulatedValue     = .false.
        real                                        :: MinForDTDecrease     = AllmostZero
        real                                        :: DefaultValue
        real                                        :: PredictedDT          = -null_real
        real                                        :: DTForNextEvent       = -null_real
        real, dimension(:, :   ), pointer           :: Matrix2D
        real, dimension(:, :, :), pointer           :: Matrix3D

        !Initialization Methods
        type (T_Layers   )                          :: Layers
        type (T_Boxes    )                          :: Boxes
        type (T_TimeSerie)                          :: TimeSerie
        type (T_ASCIIFile)                          :: ASCIIFile
        type (T_HDF      )                          :: HDF
        type (T_ProfileTimeSerie)                   :: ProfileTimeSerie
        integer                                     :: ObjEnterData         = 0 
        integer                                     :: ObjTime              = 0
        integer                                     :: ObjHorizontalGrid    = 0    
        integer                                     :: ObjGeometry          = 0    

        type(T_FillMatrix), pointer                 :: Next
    end type  T_FillMatrix

    !Global Module Variables
    type (T_FillMatrix), pointer                    :: FirstObjFillMatrix
    type (T_FillMatrix), pointer                    :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructFillMatrix2D(PropertyID, EnterDataID, TimeID,                   &
                                     HorizontalGridID, ExtractType, PointsToFill2D,     &
                                     Matrix2D, TypeZUV, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: EnterDataID
        integer                                         :: TimeID
        integer                                         :: HorizontalGridID
        integer                                         :: ExtractType
        integer, dimension(:, :), pointer               :: PointsToFill2D
        real, dimension(:, :), pointer                  :: Matrix2D
        integer                                         :: TypeZUV
        type (T_PropertyID)                             :: PropertyID
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, nUsers

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mFillMatrix_)) then
            nullify (FirstObjFillMatrix)
            call RegisterModule (mFillMatrix_) 
        endif

        call Ready(PropertyID%ObjFillMatrix, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%ObjEnterData      = AssociateInstance (mENTERDATA_,      EnterDataID     )
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)


            Me%Size2D%ILB       = lbound(PointsToFill2D, dim = 1)
            Me%Size2D%IUB       = ubound(PointsToFill2D, dim = 1)
            Me%Size2D%JLB       = lbound(PointsToFill2D, dim = 2)
            Me%Size2D%JUB       = ubound(PointsToFill2D, dim = 2)
            Me%WorkSize2D%ILB   = Me%Size2D%ILB + 1
            Me%WorkSize2D%IUB   = Me%Size2D%IUB - 1
            Me%WorkSize2D%JLB   = Me%Size2D%JLB + 1
            Me%WorkSize2D%JUB   = Me%Size2D%JUB - 1

            Me%Size3D       = T_Size3D(null_int, null_int, null_int, null_int, null_int, null_int)
            Me%Dim          = Dim2D
            Me%TypeZUV      = TypeZUV
            Me%PropertyID   = PropertyID
            
            Me%Matrix2D     => Matrix2D
            
            where (PointsToFill2D == WaterPoint) Me%Matrix2D = null_real

            if (Me%TypeZUV == TypeU_) then
                Me%Size2D%JUB       = Me%Size2D%JUB + 1
                Me%WorkSize2D%JUB   = Me%WorkSize2D%JUB + 1
            endif

            if (Me%TypeZUV == TypeV_) then
                Me%Size2D%IUB       = Me%Size2D%IUB + 1
                Me%WorkSize2D%IUB   = Me%WorkSize2D%IUB + 1
            endif

            call ReadOptions (ExtractType, PointsToFill2D = PointsToFill2D)
            
            nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
            if (nUsers == 0) stop 'ConstructFillMatrix2D - ModuleFillMatrix - ERR01' 

            
            if(Me%TimeEvolution == None)then
                PropertyID%SolutionFromFile  = .false.
            else 
                PropertyID%SolutionFromFile  = .true.
            end if

            !Returns ID
            PropertyID%ObjFillMatrix = Me%InstanceID

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
                                     PointsToFill3D, Matrix3D, TypeZUV, FillMatrix, STAT)

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
        real   , optional, intent(IN )                  :: FillMatrix
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        real                                            :: FillMatrix_
        integer                                         :: STAT_, nUsers

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mFillMatrix_)) then
            nullify (FirstObjFillMatrix)
            call RegisterModule (mFillMatrix_) 
        endif

        call Ready(PropertyID%ObjFillMatrix, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%ObjEnterData      = AssociateInstance (mENTERDATA_,      EnterDataID     )
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )

            Me%Size3D%ILB       = lbound(PointsToFill3D, dim = 1)
            Me%Size3D%IUB       = ubound(PointsToFill3D, dim = 1)
            Me%Size3D%JLB       = lbound(PointsToFill3D, dim = 2)
            Me%Size3D%JUB       = ubound(PointsToFill3D, dim = 2)
            Me%Size3D%KLB       = lbound(PointsToFill3D, dim = 3)
            Me%Size3D%KUB       = ubound(PointsToFill3D, dim = 3)
            Me%WorkSize3D%ILB   = Me%Size3D%ILB + 1
            Me%WorkSize3D%IUB   = Me%Size3D%IUB - 1
            Me%WorkSize3D%JLB   = Me%Size3D%JLB + 1
            Me%WorkSize3D%JUB   = Me%Size3D%JUB - 1
            Me%WorkSize3D%KLB   = Me%Size3D%KLB + 1
            Me%WorkSize3D%KUB   = Me%Size3D%KUB - 1


            Me%Size2D       = T_Size2D(null_int, null_int, null_int, null_int)
            Me%Dim          = Dim3D
            Me%TypeZUV      = TypeZUV
            Me%PropertyID   = PropertyID

            Me%Matrix3D     => Matrix3D


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


            call ReadOptions (ExtractType, PointsToFill3D = PointsToFill3D)
            
            nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
            if (nUsers == 0) stop 'ConstructFillMatrix3D - ModuleFillMatrix - ERR01' 

            if(Me%TimeEvolution == None)then
                PropertyID%SolutionFromFile  = .false.
            else 
                PropertyID%SolutionFromFile  = .true.
            end if

            !Returns ID
            PropertyID%ObjFillMatrix = Me%InstanceID

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

    subroutine ReadOptions (ExtractType, PointsToFill2D, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        integer                                         :: ExtractType
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D

        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: iflag
        character(len=StringLength)                     :: AuxString

        !Reads Time Evolution
        call GetData(AuxString, Me%ObjEnterData,  iflag,                            &
                     SearchType     = ExtractType,                                  &
                     keyword        = 'FILE_IN_TIME',                               &
                     default        = "None",                                       &
                     ClientModule   = 'ModuleFillMatrix',                           &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR01'

        select case (trim(adjustl(AuxString)))
            case ("None",       "NONE", "none")
                Me%TimeEvolution    = None
            case ("Hdf",        "HDF",          "hdf")
                Me%TimeEvolution    = ReadHDF
            case ("Timeserie",  "TIMESERIE",    "timeserie",    "TimeSerie")
                Me%TimeEvolution    = ReadTimeSerie
            case ("Profile_Timeserie",  "PROFILE_TIMESERIE",    "profile_timeserie",    "Profile_TimeSerie")
                Me%TimeEvolution    = ProfileTimeSerie
            case default
                write(*,*)'Invalid option for keyword FILE_IN_TIME'
                stop 'ReadOptions - ModuleFillMatrix - ERR02'
        end select


        if(Me%TimeEvolution == None)then

            call GetData(AuxString, Me%ObjEnterData,  iflag,                            &
                         SearchType     = ExtractType,                                  &
                         keyword        = 'INITIALIZATION_METHOD',                      &
                         default        = "Constant",                                   &
                         ClientModule   = 'ModuleFillMatrix',                           &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR03'


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
                case default
                    write(*,*)'Invalid option for keyword INITIALIZATION_METHOD'
                    stop 'ReadOptions - ModuleFillMatrix - ERR04'
            end select


            call GetData(Me%RemainsConstant, Me%ObjEnterData,  iflag,                   &
                         SearchType     = ExtractType,                                  &
                         keyword        = 'REMAIN_CONSTANT',                            &
                         default        = .false.,                                      &
                         ClientModule   = 'ModuleFillMatrix',                           &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR05'

            call GetData(AuxString, Me%ObjEnterData,  iflag,                            &
                         SearchType     = ExtractType,                                  &
                         keyword        = 'INITIALIZATION_DEFAULT',                     &
                         default        = "Constant",                                   &
                         ClientModule   = 'ModuleFillMatrix',                           &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR03'

            select case (trim(adjustl(AuxString)))
                case ("Constant",   "CONSTANT",   "constant")
                    Me%InitializationDefault = Constant
                case ("Profile_Timeserie",  "PROFILE_TIMESERIE",    "profile_timeserie",    "Profile_TimeSerie")
                    Me%InitializationDefault = ProfileTimeSerie
            end select
        
        end if

        !Gets the default value
        call GetData(Me%DefaultValue, Me%ObjEnterData,  iflag,                          &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'DEFAULTVALUE',                                   &
                     ClientModule   = 'ModuleFillMatrix',                               &
                     default        = FillValueReal,                                    &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR06'
        if (iflag == 0) then
            write(*,*)'Please define default value for property'//trim(Me%PropertyID%Name)
            stop 'ReadOptions - ModuleFillMatrix - ERR07'
        end if

        !Accumulitve Property (e.g. Rain from gauges)
        call GetData(Me%AccumulatedValue, Me%ObjEnterData,  iflag,                      &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'NO_INTERPOLATION',                               &
                     default        = .false.,                                          &
                     ClientModule   = 'ModuleFillMatrix',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR06b'

        !When to shut down DT?
        call GetData(Me%MinForDTDecrease, Me%ObjEnterData,  iflag,                      &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'MIN_FOR_DT_DECREASE',                            &
                     default        = AllMostZero,                                      &
                     ClientModule   = 'ModuleFillMatrix',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFillMatrix - ERR06c'


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
                            stop 'ReadOptions - ModuleFillMatrix - ERR07'
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
                            stop 'ReadOptions - ModuleFillMatrix - ERR08'
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

                    case(Profile)

                        if (Me%Dim == Dim2D) then
                            write(*,*)'Cannot Initialise 2D Matrix by Profile'
                            stop 'ReadOptions - ModuleFillMatrix - ERR09'
                        else
                            call ConstructSpaceProfile(ExtractType, PointsToFill3D)
                        endif


                    case(AnalyticProfile)

                        if (Me%Dim == Dim2D) then
                            write(*,*)'Cannot Initialise 2D Matrix by Profile'
                            stop 'ReadOptions - ModuleFillMatrix - ERR19'
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

                        call ConstructHDFInput (ExtractType)

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

                call ConstructHDFInput (ExtractType)

            case(ProfileTimeSerie)

                if (Me%Dim == Dim2D) then
                    write(*,*)'Cannot read 2D matrix by profile time serie'
                    stop 'ReadOptions - ModuleFillMatrix - ERR10'
                else
                    call ConstructProfileTimeSerie (PointsToFill3D, ExtractType)
                endif

        end select


    end subroutine ReadOptions

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

        !Gets boxes Values
        allocate (Me%Boxes%Values(BoxesNumber))

        call GetData(Me%Boxes%Values,                                            &
                     Me%ObjEnterData , iflag,                                    &
                     SearchType   = ExtractType,                                 &
                     keyword      = 'BOXES_VALUES',                              &
                     ClientModule = 'ModuleFillMatrix',                          &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceBox - ModuleFillMatrix - ERR10'

        if       (STAT_CALL .EQ. SIZE_ERR_)  then
            write(*,*) 'Incorrect number of boxes'
            stop 'ConstructSpaceBox - ModuleFillMatrix - ERR11'
        else if ((STAT_CALL .NE. SIZE_ERR_) .AND.  (STAT_CALL .NE. SUCCESS_)) then
            stop 'ConstructSpaceBox - ModuleFillMatrix - ERR12'
        end if           
        if (iflag==0) then
            write(*,*) 'Boxes Values not given'           
            stop       'ConstructSpaceBox - ModuleFillMatrix - ERR13'
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
            write(*,*)'ASCII File Name not given'
            stop 'ConstructSpaceASCIIFile - ModuleFillMatrix - ERR02'
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

    subroutine ConstructSpaceTimeSerie (ExtractType)

        !Arguments-------------------------------------------------------------
        integer                                         :: ExtractType

        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: iflag

        call GetData(Me%TimeSerie%FileName,                                      &
                     Me%ObjEnterData , iflag,                                    &
                     SearchType   = ExtractType,                                 &
                     keyword      = 'FILENAME',                                  &
                     ClientModule = 'ModuleFillMatrix',                          &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR01'

        if (iflag==0)then
            write(*,*)'Time Serie File Name not given'
            stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR02'
        endif

        call GetData(Me%TimeSerie%Column,                                        &
                     Me%ObjEnterData , iflag,                                    &
                     SearchType   = ExtractType,                                 &
                     keyword      = 'DATA_COLUMN',                               &
                     ClientModule = 'ModuleFillMatrix',                          &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR03'

        if (iflag==0)then
            write(*,*)'Data Column not given'
            stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR04'
        endif

        !Starts Time Serie
        call StartTimeSerieInput(Me%TimeSerie%ObjTimeSerie,                      &
                                 Me%TimeSerie%FileName,                          &
                                 Me%ObjTime,                                     &
                                 STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSpaceTimeSerie - ModuleFillMatrix - ERR05'

    end subroutine ConstructSpaceTimeSerie

    !--------------------------------------------------------------------------


    subroutine ConstructHDFInput (ExtractType, PointsToFill2D, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        integer                                         :: ExtractType
        integer, dimension(:, :),    pointer, optional  :: PointsToFill2D
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: iflag, HDF5_READ
        logical                                         :: exist
        type(T_Time)                                    :: Now

        !Local-----------------------------------------------------------------
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB
        type(T_Time)                                    :: LastInstantTime, EndTime
        logical                                         :: FoundSecondInstant, LastGroupEqualField
        real                                            :: Year, Month, Day, Hour, Minute, Second

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

        call GetData(Me%HDF%FileName,                                                   &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'FILENAME',                                         &
                     ClientModule = 'ModuleFillMatrix',                                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR90'

        if (iflag==0)then
            write(*,*)'HDF filename not given'
            stop 'ConstructHDFInput - ModuleFillMatrix - ERR100'
        endif


        inquire (file=trim(Me%HDF%FileName), exist = exist)
        if (.not. exist) then
            write(*,*)'Could not find file '//trim(Me%HDF%FileName)
            stop 'ConstructHDFInput - ModuleFillMatrix - ERR110'
        endif

        
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


        call ConstructHDF5 (Me%HDF%ObjHDF5, trim(Me%HDF%FileName), HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR120'


        call GetHDF5GroupNumberOfItems(Me%HDF%ObjHDF5, "/Time", &
                                       Me%HDF%NumberOfInstants, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR130'

        !if only one instant is found then values remain constant
        if(Me%HDF%NumberOfInstants == 1) Me%RemainsConstant = .true.

        call GetComputeCurrentTime(Me%ObjTime, Now, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR140'

        call GetComputeTimeLimits(Me%ObjTime, EndTime = EndTime, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR150'


i1:     if (.not. Me%HDF%Generic4D%ON) then
            
i2:         if(Me%HDF%NumberOfInstants > 1)then

                Me%HDF%PreviousInstant  = 1
                Me%HDF%NextInstant      = Me%HDF%PreviousInstant
                Me%HDF%PreviousTime     = HDF5TimeInstant(Me%HDF%PreviousInstant) 
                
                call CheckCyclicMonths(Me%HDF%PreviousTime, RefTime = Now, CyclicTimeON = Me%HDF%CyclicTimeON)

i3:             if (Me%HDF%CyclicTimeON) then
            
                    if (Me%HDF%NumberOfInstants /= 12) stop 'ConstructHDFInput - ModuleFillMatrix - ERR160'

                else i3

                    if(Me%HDF%PreviousTime .gt. Now)then
                        write(*,*)
                        write(*,*)'Could not read solution from HDF5 file'
                        write(*,*)'First file instant greater than current time'
                        write(*,*)'Matrix name: '//trim(Me%HDF%FieldName)
                        stop      'ConstructHDFInput - ModuleFillMatrix - ERR170'
                    end if
            
                    LastInstantTime = HDF5TimeInstant(Me%HDF%NumberOfInstants)
                    
                    call CheckCyclicMonths(LastInstantTime, RefTime = EndTime)


                    if(Me%TimeEvolution .ne. None)then
                        if(LastInstantTime .lt. EndTime)then
                            write(*,*)
                            write(*,*)'Could not read solution from HDF5 file'
                            write(*,*)'Last instant in file lower than simulation ending time'
                            write(*,*)'Matrix name: '//trim(Me%HDF%FieldName)
                            stop      'ConstructHDFInput - ModuleFillMatrix - ERR180'
                        end if
                    end if

                endif i3

                FoundSecondInstant = .false.
            
                !if number of instants greater than 1 then 
                !find first and second instants
d2:             do while(.not. FoundSecondInstant)
                
                    Me%HDF%PreviousInstant  = Me%HDF%NextInstant
                    Me%HDF%NextInstant      = Me%HDF%NextInstant + 1

                    if (Me%HDF%CyclicTimeON .and. Me%HDF%NextInstant .gt. Me%HDF%NumberOfInstants) then
                        Me%HDF%NextInstant  = 1
                    end if


                    Me%HDF%NextTime         = HDF5TimeInstant(Me%HDF%NextInstant)

                    call CheckCyclicMonths(Me%HDF%NextTime, RefTime = Now,              &
                                           CyclicTimeON = Me%HDF%CyclicTimeON)

                    if(Me%HDF%PreviousTime .le. Now .and. Me%HDF%NextTime .ge. Now) then
                        FoundSecondInstant  = .true.
                        exit
                    end if

                    Me%HDF%PreviousTime            = Me%HDF%NextTime

                    if(Me%HDF%NextInstant .gt. Me%HDF%NumberOfInstants .and. .not. Me%HDF%CyclicTimeON) then
                        write(*,*)
                        write(*,*)'Could not read solution from HDF5 file'
                        write(*,*)'Could not find second instant in file'
                        write(*,*)'Matrix name: '//trim(Me%HDF%FieldName)
                        stop      'ConstructHDFInput - ModuleFillMatrix - ERR190'
                    end if

                end do d2

            elseif(Me%HDF%NumberOfInstants == 1)then i2

                Me%HDF%PreviousInstant  = 1
                Me%HDF%NextInstant      = Me%HDF%PreviousInstant

                Me%HDF%PreviousTime     = HDF5TimeInstant(Me%HDF%PreviousInstant)

                call CheckCyclicMonths(Me%HDF%PreviousTime, RefTime = Now)


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
                stop 'ConstructHDFInput - ModuleFillMatrix - ERR200'

            end if i2


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

                if (Me%HDF%PreviousInstant /= Me%HDF%NextInstant) then

                    call InterpolateMatrix3DInTime(ActualTime       = Now,                         &
                                                   Size             = Me%WorkSize3D,               &
                                                   Time1            = Me%HDF%PreviousTime,         &
                                                   Matrix1          = Me%HDF%PreviousField3D,      &
                                                   Time2            = Me%HDF%NextTime,             &
                                                   Matrix2          = Me%HDF%NextField3D,          &
                                                   MatrixOut        = Me%Matrix3D,                 &
                                                   PointsToFill3D   = PointsToFill3D)

                else

                    !Prev and next are equal (last instant?)
                    Me%Matrix3D(:,:,:)  = Me%HDF%NextField3D(:,:,:)

                endif

            end if i4

        endif i1


    end subroutine ConstructHDFInput

    !--------------------------------------------------------------------------
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
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleFillMatrix - ERR10'


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
        
        call HDF5SetLimits  (Me%HDF%ObjHDF5, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = Me%HDF%ObjHDF5,                           &
                             GroupName      = "/Time",                                  &
                             Name           = "Time",                                   &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleFillMatrix - ERR01'

        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))

                                     
        deallocate(TimeVector)

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
    subroutine ReadHDF5Values2D (Instant, Field)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        real, dimension(:,:), pointer           :: Field
        
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (Me%HDF%ObjHDF5, Me%WorkSize2D%ILB, Me%WorkSize2D%IUB,      &
                             Me%WorkSize2D%JLB, Me%WorkSize2D%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR10'
        

        call HDF5ReadData(Me%HDF%ObjHDF5, trim(Me%HDF%VGroupPath),                      &
                          trim(Me%HDF%FieldName),                                       &
                          Array2D = Field, OutputNumber = Instant, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values2D - ModuleFillMatrix - ERR40'

        if(Me%HDF%HasMultiplyingFactor)then
            Field = Field * Me%HDF%MultiplyingFactor
        end if

        if(Me%HDF%HasAddingFactor)then
            Field = Field + Me%HDF%AddingFactor
        end if

    end subroutine ReadHDF5Values2D
    
    
    !--------------------------------------------------------------------------

    
    subroutine ReadHDF5Values3D (Instant, Field)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        real, dimension(:,:,:), pointer         :: Field

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL, i, j, k, ILB, IUB, JLB, JUB, KLB, KUB

        !Begin-----------------------------------------------------------------


         ILB = Me%WorkSize3D%ILB
         IUB = Me%WorkSize3D%IUB
         JLB = Me%WorkSize3D%JLB
         JUB = Me%WorkSize3D%JUB

        if (Me%HDF%From2Dto3D) then
            KLB = 1
            KUB = 1
        else
            KLB = Me%WorkSize3D%KLB
            KUB = Me%WorkSize3D%KUB
                     
            Me%HDF%ReadField3D => Field
        endif


        call HDF5SetLimits  (Me%HDF%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR10'
             
        call HDF5ReadData(Me%HDF%ObjHDF5, trim(Me%HDF%VGroupPath),                      &
                          trim(Me%HDF%FieldName),                                       &
                          Array3D = Me%HDF%ReadField3D, OutputNumber = Instant, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleFillMatrix - ERR02'

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
            Field = Field * Me%HDF%MultiplyingFactor
        end if

        if(Me%HDF%HasAddingFactor)then
            Field = Field + Me%HDF%AddingFactor
        end if

    end subroutine ReadHDF5Values3D


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

    subroutine GetFillMatrixDTPrediction (FillMatrixID, PredictedDT, DTForNextEvent, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: FillMatrixID
        real, intent(OUT)                               :: PredictedDT, DTForNextEvent
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FillMatrixID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            PredictedDT     = Me%PredictedDT
            DTForNextEvent  = Me%DTForNextEvent

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetFillMatrixDTPrediction


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

                end select

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
        real                                            :: Generic_4D_Value_

        !Local----------------------------------------------------------------

        !Begin----------------------------------------------------------------

        if (Me%HDF%Generic4D%ON) then

            call ModifyHDFInput3DGeneric4D(PointsToFill3D, Generic_4D_Value_)

        else
            
            call ModifyHDFInput3DTime(PointsToFill3D)

        endif

    end subroutine ModifyHDFInput3D

    !----------------------------------------------------------------------------

    subroutine ModifyHDFInput3DTime(PointsToFill3D)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL, n
        type (T_Time)                                   :: Now

        !Begin----------------------------------------------------------------


        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyHDFInput3D - ModuleFillMatrix - ERR10'

        if(Now .ge. Me%HDF%NextTime)then

            n = 0

            do while (Now .ge. Me%HDF%NextTime) 
        
                Me%HDF%PreviousInstant  = Me%HDF%NextInstant
        
                if(Me%HDF%NextInstant .lt. Me%HDF%NumberOfInstants) then
                    Me%HDF%NextInstant  = Me%HDF%NextInstant + 1
                else
                    if (Me%HDF%CyclicTimeON) then
                        Me%HDF%NextInstant  = 1
                    else
                        exit
                    endif
                end if
        
                Me%HDF%PreviousTime     = Me%HDF%NextTime
                Me%HDF%NextTime         = HDF5TimeInstant(Me%HDF%NextInstant)

                n = n + 1

            enddo

            if(Now .gt. Me%HDF%NextTime)then
                write(*,*)
                write(*,*)'Could not read solution from HDF5 file'
                write(*,*)'Time instants inconsistency.'
                stop      'ModifyHDFInput3D - ModuleFillMatrix - ERR20'
            end if

            if (n==1) then
                call SetMatrixValue(Me%HDF%PreviousField3D, Me%WorkSize3D, Me%HDF%NextField3D)
            else
                call ReadHDF5Values3D(Me%HDF%PreviousInstant, Me%HDF%PreviousField3D)
            endif

            call ReadHDF5Values3D(Me%HDF%NextInstant, Me%HDF%NextField3D)

        end if

        if (Me%HDF%PreviousInstant /= Me%HDF%NextInstant) then

            call InterpolateMatrix3DInTime(ActualTime       = Now,                         &
                                           Size             = Me%WorkSize3D,               &
                                           Time1            = Me%HDF%PreviousTime,         &
                                           Matrix1          = Me%HDF%PreviousField3D,      &
                                           Time2            = Me%HDF%NextTime,             &
                                           Matrix2          = Me%HDF%NextField3D,          &
                                           MatrixOut        = Me%Matrix3D,                 &
                                           PointsToFill3D   = PointsToFill3D)

        else
            
            !Prev and next are equal (last instant?)
            call SetMatrixValue(Me%Matrix3D, Me%WorkSize3D, Me%HDF%NextField3D)

        endif


    end subroutine ModifyHDFInput3DTime


    !--------------------------------------------------------------------------


    subroutine ModifyHDFInput3DGeneric4D(PointsToFill3D, Generic_4D_Value_)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        real                                            :: Generic_4D_Value_

        !Local----------------------------------------------------------------
        integer                                         :: PrevI, NextI
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
        real                                            :: Generic_4D_Value_

        !Local----------------------------------------------------------------

        !Begin----------------------------------------------------------------

        if (Me%HDF%Generic4D%ON) then

            call ModifyHDFInput2DGeneric4D(PointsToFill2D, Generic_4D_Value_)
            
        else

            call ModifyHDFInput2DTime(PointsToFill2D)

        endif

    end subroutine ModifyHDFInput2D

    !----------------------------------------------------------------------------



    subroutine ModifyHDFInput2DTime(PointsToFill2D)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL, n
        type (T_Time)                                   :: Now

        !Begin----------------------------------------------------------------

        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyHDFInput2DTime - ModuleFillMatrix - ERR10'

        if (Now .ge. Me%HDF%NextTime)then

            n = 0

            do while (Now .ge. Me%HDF%NextTime) 

                Me%HDF%PreviousInstant  = Me%HDF%NextInstant

                if(Me%HDF%NextInstant .lt. Me%HDF%NumberOfInstants)then
                    Me%HDF%NextInstant  = Me%HDF%NextInstant + 1
                else
                    if (Me%HDF%CyclicTimeON) then 
                        Me%HDF%NextInstant  = 1
                    else
                        exit
                    endif
                end if

                Me%HDF%PreviousTime     = Me%HDF%NextTime
                Me%HDF%NextTime         = HDF5TimeInstant(Me%HDF%NextInstant)

                n = n + 1

            enddo

            if(Now .gt. Me%HDF%NextTime)then
                write(*,*)
                write(*,*)'Could not read solution from HDF5 file'
                write(*,*)'Time instants inconsistency.'
                stop      'ModifyHDFInput2DTime - ModuleFillMatrix - ERR20'
            end if

            if (n==1) then 
                call SetMatrixValue(Me%HDF%PreviousField2D, Me%WorkSize2D, Me%HDF%NextField2D, PointsToFill2D)
            else
                call ReadHDF5Values2D(Me%HDF%PreviousInstant, Me%HDF%PreviousField2D)
            endif

            call ReadHDF5Values2D(Me%HDF%NextInstant, Me%HDF%NextField2D)

        end if

        if (Me%HDF%PreviousInstant /= Me%HDF%NextInstant) then

            if (Me%AccumulatedValue) then       !For Rain
                Me%Matrix2D = Me%HDF%NextField2D / (Me%HDF%NextTime - Me%HDF%PreviousTime)
    
                call PredictDTForHDF (PointsToFill2D, Me%HDF%PreviousTime, Me%HDF%NextTime, Now)
                
            else

                !Interpolates the two matrixes in time
                call InterpolateMatrix2DInTime(ActualTime       = Now,                         &
                                               Size             = Me%WorkSize2D,               &
                                               Time1            = Me%HDF%PreviousTime,         &
                                               Matrix1          = Me%HDF%PreviousField2D,      &
                                               Time2            = Me%HDF%NextTime,             &
                                               Matrix2          = Me%HDF%NextField2D,          &
                                               MatrixOut        = Me%Matrix2D,                 &
                                               PointsToFill2D   = PointsToFill2D)
            endif
            
        else

            !Prev and next are equal (last instant?)
            if (.not. Me%AccumulatedValue) then

                call SetMatrixValue(Me%Matrix2D, Me%WorkSize2D, Me%HDF%PreviousField2D, PointsToFill2D)

            else
                !do nothing. Values will not be used 
            endif

        endif

    end subroutine ModifyHDFInput2DTime

    !--------------------------------------------------------------------------

    subroutine ModifyHDFInput2DGeneric4D(PointsToFill2D, Generic_4D_Value_)
        
        !Arguments------------------------------------------------------------
        integer, dimension(:, :), pointer               :: PointsToFill2D
        real                                            :: Generic_4D_Value_

        !Local----------------------------------------------------------------
        integer                                         :: PrevI, NextI
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

    subroutine ModifyProfileTimeSerie(PointsToFill3D)    
    
        !Arguments------------------------------------------------------------
        integer, dimension(:, :, :), pointer            :: PointsToFill3D

        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        type (T_Time)                                   :: Now

        !Begin----------------------------------------------------------------

        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyHDFInput3D - ModuleFillMatrix - ERR01'

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
                stop      'ModifyProfileTimeSerie - ModuleFillMatrix - ERR01'
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
        if (STAT_CALL /= SUCCESS_) stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR01'


        !Gets Value for current Time
        call GetTimeSerieValue (Me%TimeSerie%ObjTimeSerie, Now,                  &
                                Me%TimeSerie%Column,                             &
                                Time1, Value1, Time2, Value2, TimeCycle,         &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR02'

        
        if (Me%PropertyID%IDNumber /= WindDirection_) then

            if (TimeCycle) then
                NewValue = Value1

            else

                if (Me%AccumulatedValue) then       !For Rain
                    NewValue = Value2 / (Time2 - Time1)
    
                    call GetTimeSerieDTForNextEvent (Me%TimeSerie%ObjTimeSerie,          &
                                                     NewValue, Me%TimeSerie%Column, Now, &
                                                     Me%PredictedDT, Me%DTForNextEvent,  &
                                                     STAT  = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR03'
                
                else
                    !Interpolates Value for current instant
                    call InterpolateValueInTime(Now, Time1, Value1, Time2, Value2, NewValue)
                endif
            endif

        else

            !Gets Grid Angle
            call GetGridAngle(Me%ObjHorizontalGrid, Angle, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifySpaceTimeSerie - ModuleFillMatrix - ERR03'

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

        if (Me%Dim == Dim2D) then
            call SetMatrixValue(Me%Matrix2D, Me%WorkSize2D, NewValue, PointsToFill2D)
        else
            call SetMatrixValue(Me%Matrix3D, Me%WorkSize3D, NewValue, PointsToFill3D)
        endif

    end subroutine ModifySpaceTimeSerie

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

                if (Me%HDF%ObjHDF5 /= 0) then

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
                    
                    call KillHDF5(Me%HDF%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillFillMatrix - ModuleFillMatrix - ERR30'

  
                    if (Me%HDF%Generic4D%ReadFromTimeSerie) then
                        call KillTimeSerie(Me%HDF%Generic4D%ObjTimeSerie, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillFillMatrix - ModuleFillMatrix - ERR40'
                    endif
  
                endif

                if (Me%ObjGeometry /= 0) then
                    nUsers = DeassociateInstance (mGEOMETRY_, Me%ObjGeometry)
                    if (nUsers == 0) stop 'KillFillMatrix - ModuleFillMatrix - ERR50'
                endif

                nUsers = DeassociateInstance (mTIME_,           Me%ObjTime           )
                if (nUsers == 0) stop 'KillFillMatrix - ModuleFillMatrix - ERR60'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid )
                if (nUsers == 0) stop 'KillFillMatrix - ModuleFillMatrix - ERR70'
                


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
