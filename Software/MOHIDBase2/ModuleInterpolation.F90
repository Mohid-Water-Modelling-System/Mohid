!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : Interpolation
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : «Date»
! REVISION      : Ricardo Lemos & Guillaume Riflet - v4.0
! DESCRIPTION   : «Description»
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

!DataFile
!   TYPE_ZUV           : integer              -        !Where pointa are defined in the cell (Z - center; U - Face U; V : Face V) 
!   METHODOLOGY        : integer          [4]          !The methodology use in the interpolation process
                                                            !ConvolutionConserv_    =  1
                                                            !ConvolutionNonConserv_ =  2
                                                            !Triangulation_         =  3
                                                            !Bilinear_              =  4
                                                            !Spline2D_              =  5
                                                            !InverseWeight_         =  6
!   POLI_DEGREE_VERT   : integer          [1]          !The order of the polinomial use to interpolate in the vertical 
!   N_DIM              : integer          [3]          !The number of dimensions
!   MAX_DISTANCE       : real             [*]          !Max distance for points to be consider in the inverse weight interpolation
!   IWD_N              : real             [2.0]        !Coefficent use in the inverse weight interpolation
!   KERNEL_TYPE        : char             [CharGaussian_] ! Type of kernel use in the convolution interpolations 
!   EXTRAPOLATE_PROFILE: logical          [true]       !Chek if the user wants to extrapolate in the vertical
!   EXTRAPOLATE_2D     : integer          [0] NotActive_   !??? Paulo?
!   NON_CONSERVATIVE   : logical          [true]        !Checks if the user wants to use the NonConservative convolution process
!   NC_TYPE            : integer          [2]           !Cheks what class of NonConservative convolution process to use
!   PHI                : real             [.9]          !Smoothing parameter. Gives the degree of smoothing in the interpolated 
!                                                       !field. Its range is ]0,1].
!   N_GROUPS           : integer          [1]           !Number of groups generated for each dimension in the data-oriented 
!                                                       !convolution.
!   SAMPLE_SIZE        : integer          [50] SampleSize_ !Number of observations needed for the logistic regression in the 
!                                                          !data-oriented convolution.
!   MAX_ITERATIONS     : integer          [20] MaxIterations_ !Maximum number of iterations allowed in the logistic regression 
!                                                             !in the data-oriented convolution.

module ModuleInterpolation
        
    use ModuleGlobalData
    use ModuleTime
    use ModuleFunctions,        only : RodaXY, InterpolateProfileR8, PolIntProfile
    use ModuleTriangulation
    use ModuleDrawing
    use ModuleEnterData,        only : GetData
    use ModuleFunctions,        only : ConstructPropertyID, IsOdd
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, UnGetHorizontalGrid, ConstructFatherGridLocation, &
                                       InterpolRegularGrid, GetGridOrigin, GetGridAngle, GetHorizontalGridSize, &
                                       GetGridLatitudeLongitude

    use ModuleHorizontalMap,    only : GetWaterPoints2D, UngetHorizontalMap
    use ModuleMap,              only : ConstructMap, GetWaterPoints3D, GetOpenPoints3D, UngetMap
    use ModuleGeometry,         only : ConstructGeometry, ComputeInitialGeometry, GetGeometryDistances, &
                                       GetGeometrySize,   GetGeometryInputFile,   UnGetGeometry


    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: Construct_Interpolation3D

    private ::      AllocateInstance
    private ::      ReadKeywordOptions

    private ::      Construct_NonConvoStruct
    private ::          ConstructCLocalFatherGrid

    private ::      Construct_NonConvoNonStruct
    !private ::         Paulo: aqui deves por as tuas rotinas de construcao das antigas rotinas de interpolacao

    private ::      Construct_ConvoNCStructured

    private ::      Construct_ConvoNCNonStructured
    private ::          AllocateExternal_ConvoNonStruct
    private ::          AllocateTXYZP_points
    private ::          ConstructPoints2Grid_3D
    private ::              sorter

    private ::      Construct_ConvoConservative
    private ::          PopulateArrays

    private ::      AllocateExternal_ConvoStruct
    private ::      AllocateExternal
    private ::      AllocateGrid
    private ::      AllocateField
    private ::      ConstructGrid
    private ::          ReadLockExternalVar
    private ::              ReadLockExternalZ
    private ::              ReadLockExternalU
    private ::              ReadLockExternalV
    private ::          PopulateParameters
    private ::      ConstructNC_User3D
    private ::      ConstructNC_Smoothing3D
    private ::      ConstructNC_Data3D
    private ::      ConstructLocalFather
    private ::          AllocateLocalFatherGrid
    private ::          AllocateLocalFatherField

!    private ::          Relate_FatherSon
!    private ::              Gold
!    private ::              Corners
!    private ::          InFatherCell

    !Selector
    public  ::  GetInterpolationPointer
    public  ::  GetInterpolationInteger
    public  ::  UnGetInterpolation3D_I
    public  ::  UnGetInterpolation3D_R8

    !Modifier
    public  ::  ModifyInterpolator
    !private ::      NonConservativeModifier
    private ::          NCInterpolation_3D
    private ::              NC_User_3D
    private ::              NC_Smoothing_3D
    private ::              NC_Data_3D
    private ::          GetSample_3D
    private ::          LogRegr_3D
    private ::          MasterG_3D
    private ::          FinalConvolution_3D
!    private ::      ConservativeModifier
    private ::          CellSelector
    private ::          FatherToSonGrid
    private ::          Integrator
    private ::      ClassicInterpolaStruct3D
    private ::      ClassicInterpolaStruct2D
    private ::          StructureBilinear
    private ::          StructureSpline2D
    private ::          StructureTriangulation
    private ::          StructureInverseWeight 
!    private ::      ClassicInterpolaUnStruct3D
    private ::      ClassicInterpolaUnStruct2D
!    private ::          UnStructureBilinear
!    private ::          UnStructureSpline2D
    private ::          UnStructureTriangulation
    private ::          UnStructureInverseWeight 

    !Destructor (LIFO)
    public  ::  Kill_Interpolation3D

    private ::      Kill_NonConvoStruct
!    private ::      KillFatherGridLocation
!    private ::          KillCLocalFatherGrid

    private ::      Kill_NonConvoNonStruct

    private ::      Kill_ConvoNCStructured
    private ::          KillGrid
    private ::              ReadUnlockExternalVar
    private ::              ReadUnlockExternalZUV
    private ::          Kill_NC
    private ::          DeallocateExternal_ConvoStruct
    private ::              DeallocateField
    private ::              DeallocateGrid
    private ::          DeallocateExternal

    private ::       Kill_ConvoNCNonStructured
    private ::          KillPoints2Grid_3D
    private ::              KillLocalFather
    private ::                  DeallocateLocalFatherField
    private ::                  DeallocateLocalFatherGrid
    private ::          DeallocateTXYZP_points
    private ::          DeallocExtern_ConvoNonStruct

    private ::      Kill_ConvoConservative
    private ::      DeallocateArrays

    private ::      UnreadKeywordOptions
    private ::      DeallocateInstance

    !Management
    private ::  Ready
    private ::      LocateObjInterpolation
    
    !Interfaces----------------------------------------------------------------

    interface       ModifyInterpolator
        module procedure ModifyInterpolator2D
        module procedure ModifyInterpolator3D
    end interface   ModifyInterpolator


    !interface       ConservativeModifier
        !module procedure ConservativeModifier_2D
        !module procedure ConservativeModifier_3D
    !end interface   ConservativeModifier

    !interface       NonConservativeModifier
    !    module procedure NonConservativeModifier_2D
    !   module procedure NonConservativeModifier_3D
    !end interface   NonConservativeModifier


    interface       CellSelector
!        module procedure CellSelector_2D
        module procedure CellSelector_3D
    end interface   CellSelector

    interface       FatherToSonGrid
!        module procedure FatherToSonGrid_2D
        module procedure FatherToSonGrid_3D
    end interface   FatherToSonGrid

    interface       Integrator
!        module procedure Integrator_2D
        module procedure Integrator_3D
    end interface   Integrator

    interface       Relate_FatherSon
        module procedure Relate_FatherSon_3D
!        module procedure PrivateInterface2
    end interface   Relate_FatherSon

    interface           Gold
        module procedure Gold_3D
!        module procedure PrivateInterface2
    end interface       Gold

    interface           Corners
        module procedure Corners_3D
!        module procedure PrivateInterface2
    end interface       Corners

    interface       InFatherCell
        module procedure InFatherCell_3D
!        module procedure PrivateInterface2
    end interface   InFatherCell

    !Parameter-----------------------------------------------------------------
    integer, parameter     :: SampleSize_           =  50
    integer, parameter     :: MaxIterations_        =  20

    !Interpolation
    integer,    parameter  :: NotActive_                  =  0                 

    !Interpolation 3D
    integer,    parameter  :: ConvolutionConserv_         =  1
    integer,    parameter  :: ConvolutionNonConserv_      =  2

    !Interpolation 2D
    integer,    parameter  :: Triangulation_              =  3
    integer,    parameter  :: Bilinear_                   =  4
    integer,    parameter  :: Spline2D_                   =  5
    integer,    parameter  :: InverseWeight_              =  6

    !Interpolation 1D
    integer,    parameter  :: LinearProfile_              =  1
    integer,    parameter  :: PolinomialProfile_          =  2


    !Extrapolation 2D
    integer,    parameter  :: MediumTriang                =  1
    integer,    parameter  :: HighTriang                  =  2
    integer,    parameter  :: NearestNeighbor             =  3

    !NonConservative
    integer, parameter :: NC_User_                   = 1
    integer, parameter :: NC_Smoothing_              = 2
    integer, parameter :: NC_Data_                   = 3
    real,    parameter      :: Default_Phi_               =  .9d0




    character(len=StringLength), parameter :: CharGaussian_    = "Gaussian"
    character(len=StringLength), parameter :: CharExponential_ = "Exponential"

    !Types---------------------------------------------------------------------
    private :: T_NC
    type       T_NC
        real,   pointer, dimension(:)           :: TauM
        real,   pointer, dimension(:)           :: GroupCenterX
        real,   pointer, dimension(:)           :: GroupCenterY
        real,   pointer, dimension(:)           :: GroupCenterZ
        real,   pointer, dimension(:,:)         :: Tau, TauR
        real,   pointer, dimension(:,:,:)       :: Response, Predictor
    end type   T_NC

    private :: T_Grid
    type       T_Grid

        real,    pointer, dimension(:,:)        :: DX, DY, CX, CY
        real,    pointer, dimension(:,:,:)      :: DZ, CZ, SZZ
        real,    pointer, dimension(:, :)       :: XX_IE, YY_IE
        real,    pointer, dimension(:,:)        :: CenterX, CenterY
        real,    pointer, dimension(:,:,:)      :: CenterZ
        integer, pointer, dimension(:,:,:)      :: OpenPoints3D
        integer, pointer, dimension(:,:,:)      :: WaterPoints3D
        integer, pointer, dimension(:,:)        :: WaterPoints2D
        real,    pointer, dimension(:,:,:)      :: p

        Type (T_Size3D), pointer                :: Size3D
        Type (T_Size3D), pointer                :: WorkSize3D
        Type (T_Size2D), pointer                :: Size2D
        Type (T_Size2D), pointer                :: WorkSize2D
 
        integer                                 :: ObjTime                  = 0

        integer                                 :: ObjHorizontalGrid        = 0

        integer                                 :: ObjHorizontalMap         = 0

        integer                                 :: ObjBathymetry            = 0

        integer                                 :: ObjGeometry              = 0

        integer                                 :: ObjMap                   = 0

    end type T_Grid

    type       T_Field
        character(len=StringLength)             :: Name
        character(len=StringLength)             :: Units
        integer                                 :: IDNumber
!        type(T_Time)                            :: Time
        real, dimension(:,:  ),     pointer     :: Values2D
        real, dimension(:,:,:),     pointer     :: Values3D
        type(T_Field),              pointer     :: Next
    end type  T_Field

    private :: T_External
    type       T_External
        Type (T_Grid),      pointer             :: FatherGrid, SonGrid
        Type (T_Field),     pointer             :: FatherField, SonField
        Type (T_XYZPoints), pointer             :: TXYZP_Points
    end type T_External

    private :: T_Options
    type       T_Options
        integer                                 :: TypeZUV
        integer                                 :: Methodology, PoliDegreeVert, Extrapolate2D
        logical                                 :: ExtrapolateProfile
        real                                    :: IWDn, MaxDistance
        character(len=StringLength)             :: KernelType
        integer                                 :: NC_type
        integer                                 :: n_groups
        integer                                 :: n_dimensions
        integer                                 :: sample_size
        real                                    :: phi
        integer                                 :: max_iterations
    end type T_Options

    private ::  T_Interpolation
    type        T_Interpolation     

        integer                                 :: InstanceID

        integer                                 :: ObjEnterData
        integer                                 :: ObjTriangulation

        Type (T_External), pointer              :: ExternalVar
        Type (T_NC),       pointer              :: NC

        integer,  pointer,dimension(:,:,:)   :: MatrixUpdate3D
        integer,  pointer,dimension(:,:  )   :: FatherSonX
        integer,  pointer,dimension(:,:  )   :: FatherSonY
        integer,  pointer,dimension(:,:,:)   :: FatherSonZ
        character(3),pointer,dimension(:,:,:)   :: InsideFatherCell

        logical                                 :: ConvolutionApproach
        logical                                 :: NonConservative
        logical                                 :: OkZ, OkU, OkV
        logical                                 :: StructuredData
        Type (T_Field),     pointer             :: LocalFatherField
        Type (T_Grid),      pointer             :: LocalFatherGrid
        Type(T_Options)                         :: ComputeOptions

        Type (T_Interpolation), pointer         :: Next

    end type    T_Interpolation

    !Global Module Variables---------------------------------------------------
    type (T_Interpolation), pointer             :: FirstObjInterpolation
    type (T_Interpolation), pointer             :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Construct_Interpolation3D( InterpolationID,                              &
                                          EnterDataID,                                  &
                                          FromWhere,                                    &
                                          TimeSonID,                                    &
                                          HorizontalGridSonID,                          &
                                          HorizontalMapSonID,                           &    
                                          GeometrySonID,                                &
                                          MapSonID,                                     &
                                          BathymetrySonID,                              &
                                          TXYZP_Points,                                 &
                                          HorizontalGridFatherID,                       &
                                          HorizontalMapFatherID,                        &
                                          GeometryFatherID,                             &
                                          MapFatherID,                                  &
                                          STAT)

        !Arguments---------------------------------------------------------------
        integer           , intent(INOUT)               :: InterpolationID
        integer           , intent(IN )                 :: EnterDataID
        integer           , intent(IN )                 :: FromWhere
        integer           , intent(IN )                 :: TimeSonID
        integer           , intent(IN )                 :: HorizontalGridSonID
        integer           , intent(IN )                 :: GeometrySonID
        integer           , intent(IN )                 :: MapSonID
        integer           , intent(IN )                 :: HorizontalMapSonID
        integer, optional,  intent(IN )                 :: BathymetrySonID
 
        Type (T_XYZPoints), pointer, optional           :: TXYZP_Points

        integer, optional,  intent(IN )                 :: HorizontalGridFatherID
        integer, optional,  intent(IN )                 :: HorizontalMapFatherID
        integer, optional,  intent(IN )                 :: GeometryFatherID
        integer, optional,  intent(IN )                 :: MapFatherID
        integer, optional,  intent(OUT)                 :: STAT


        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: ready_         
        integer                                         :: STAT_

        !------------------------------------------------------------------------
    
        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mInterpolation_)) then
            nullify (FirstObjInterpolation)
            call RegisterModule (mInterpolation_) 
        endif

        call Ready(InterpolationID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            !Allocates a new Instance
            call AllocateInstance

            !Fills the ComputeOptions entries (n_dimensions, phi, n_obs, etc..)
            call ReadKeywordOptions (EnterDataID, FromWhere) 

            if( present(TXYZP_Points)) then 
                Me%StructuredData = .false.
            else
                Me%StructuredData = .true.
            end if

            !Classical approaches-------------------------------------------------------
cd33:       if(.not. Me%ConvolutionApproach) then

cd34:           if (Me%StructuredData) then

                    call Construct_NonConvoStruct(TimeSonID,                                &
                                              HorizontalGridSonID,                          &
                                              BathymetrySonID,                              &
                                              HorizontalMapSonID,                           &    
                                              HorizontalGridFatherID,                       &
                                              STAT_CALL)
                    if(STAT_CALL /= SUCCESS_) stop 'Construct_Interpolation3D - ModuleInterpolation - ERR00'

                else

                    call Construct_NonConvoNonStruct(STAT_CALL)
                    if(STAT_CALL /= SUCCESS_) stop 'Construct_Interpolation3D - ModuleInterpolation - ERR00b'

                end if cd34

           !Convolution processes------------------------------------------------------
            else
                
                !NonConservative convolution processes---------------------
cd55:           if(Me%NonConservative) then

                    allocate(Me%NC)

                    !Structured data-----------------------
cd56:               if(Me%StructuredData) then

                        call Construct_ConvoNCStructured (HorizontalGridSonID,                          &
                                          HorizontalMapSonID,                          &
                                          GeometrySonID,                                &
                                          MapSonID,                                     &
                                          HorizontalGridFatherID,                       &
                                          HorizontalMapFatherID,                       &
                                          GeometryFatherID,                             &
                                          MapFatherID,                                  &
                                          STAT_CALL)
                        if(STAT_CALL /= SUCCESS_) stop 'Construct_Interpolation3D - ModuleInterpolation - ERR01'

                    !NonStructured data-------------------
                    else                       

                        call Construct_ConvoNCNonStructured (HorizontalGridSonID,       &
                                          HorizontalMapSonID,                           &
                                          GeometrySonID,                                &
                                          MapSonID,                                     &
                                          TXYZP_Points,                                 &
                                          STAT_CALL)
                        if(STAT_CALL /= SUCCESS_) stop 'Construct_Interpolation3D - ModuleInterpolation - ERR02'

                    end if cd56

                !Conservative convolution processes-------------------------
                else

                    call Construct_ConvoConservative(HorizontalGridSonID,               &
                                          HorizontalMapSonID,                           &
                                          GeometrySonID,                                &
                                          MapSonID,                                     &
                                          HorizontalGridFatherID,                       &
                                          HorizontalMapFatherID,                        &
                                          GeometryFatherID,                             &
                                          MapFatherID,                                  &
                                          STAT_CALL)
                   if(STAT_CALL /= SUCCESS_) stop 'Construct_Interpolation3D - ModuleInterpolation - ERR03'

                end if cd55

            end if cd33

            !Returns ID
            InterpolationID = Me%InstanceID
            STAT_ = SUCCESS_

        else 
            
            stop 'Construct_Interpolation - ModuleInterpolation - ERR99' 

            STAT_ = ID_ERR_

        end if cd0

        if (present(STAT)) STAT = STAT_
    
    end subroutine Construct_Interpolation3D

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !PAULO: ainda ha trabalho para fazeres aqui ...
    subroutine Construct_NonConvoStruct(  TimeSonID,                                    &
                                          HorizontalGridSonID,                          &
                                          BathymetrySonID,                              &
                                          HorizontalMapSonID,                           &    
                                          HorizontalGridFatherID,                       &
                                          STAT)

        !Arguments-----------------------------------------------------------
        integer                                         ::  TimeSonID,                                    &
                                                            HorizontalGridSonID,                          &
                                                            BathymetrySonID,                              &
                                                            HorizontalMapSonID,                           &    
                                                            HorizontalGridFatherID
        integer, optional, intent(OUT)                  ::  STAT

        !Locals--------------------------------------------------------------
        integer                                         :: STAT_, STAT_CALL

        STAT_ = UNKNOWN_

        call ConstructCLocalFatherGrid (TimeSonID,       HorizontalGridSonID,          &
                                                   BathymetrySonID, HorizontalMapSonID)

        if (Me%ComputeOptions%Methodology == Bilinear_) then
 
          call ConstructFatherGridLocation(HorizontalGridSonID, HorizontalGridFatherID, & 
                                           OkZ = Me%OkZ, OkU = Me%OkU, OkV = Me%OkV, STAT = STAT_CALL)                  
          if (STAT_CALL /= SUCCESS_) stop 'Construct_Interpolation - ModuleInterpolation - ERR89' 

        endif

        STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

    end subroutine Construct_NonConvoStruct

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !PAULO: ainda ha trabalho para fazeres aqui ...
    subroutine Construct_NonConvoNonStruct(STAT)

        !Arguments-----------------------------------------------------------
        integer, optional, intent(OUT)                  ::  STAT

        !Locals--------------------------------------------------------------
        integer                                         :: STAT_

        STAT_ = UNKNOWN_

        !Aqui ha coisas a fazer ...

        !E necessario usar a variavel STAT_ para detectar erros no programa ...
        STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

    end subroutine Construct_NonConvoNonStruct

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine Construct_ConvoConservative (HorizontalGridSonID,                        &
                                            HorizontalMapSonID,                         &
                                            GeometrySonID,                              &
                                            MapSonID,                                   &
                                            HorizontalGridFatherID,                     &
                                            HorizontalMapFatherID,                      &
                                            GeometryFatherID,                           &
                                            MapFatherID,                                &
                                            STAT)

        !Arguments-----------------------------------------------------------
        integer                                         ::  HorizontalGridSonID,                          &
                                                            HorizontalMapSonID,                          &
                                                            GeometrySonID,                                &
                                                            MapSonID,                                     &
                                                            HorizontalGridFatherID,                       &
                                                            HorizontalMapFatherID,                       &
                                                            GeometryFatherID,                             &
                                                            MapFatherID                                  
        integer, optional, intent(OUT)                  ::  STAT

        !Locals--------------------------------------------------------------
        integer                                         :: STAT_
        type(T_Size3D), pointer                         :: SizeSon, SizeFather

        STAT_ = UNKNOWN_

        call AllocateExternal_ConvoStruct(STAT_)

        call ConstructGrid(Me%ExternalVar%FatherGrid,                           &
                                                    HorizontalGridFatherID,                         &
                                                    HorizontalMapFatherID,                         &
                                                    GeometryFatherID,                               &
                                                    MapFatherID, STAT_)

        call ConstructGrid(Me%ExternalVar%SonGrid,                           &
                                                    HorizontalGridSonID,                         &
                                                    HorizontalMapSonID,                         &
                                                    GeometrySonID,                               &
                                                    MapSonID, STAT_)

        !Shorten variable names
        SizeFather => Me%ExternalVar%FatherGrid%WorkSize3D
        SizeSon => Me%ExternalVar%SonGrid%WorkSize3D

        call ConstructLocalFather(  SizeFather%ILB,   &
                                    SizeFather%JLB,   &
                                    SizeFather%KLB,   &
                                    SizeFather%IUB,   &
                                    SizeFather%JUB,   &
                                    SizeFather%KUB,   &
                                    STAT_)

        !Allocate Me%Matrixupdate etc...
        allocate(Me%MatrixUpdate3D( SizeFather%ILB:SizeFather%IUB, &
                                    SizeFather%JLB:SizeFather%JUB, &
                                    SizeFather%KLB:SizeFather%KUB))    

        allocate(Me%FatherSonX(     SizeSon%ILB:SizeSon%IUB, &
                                    SizeSon%JLB:SizeSon%JUB))    

        allocate(Me%FatherSonY(     SizeSon%ILB:SizeSon%IUB, &
                                    SizeSon%JLB:SizeSon%JUB))    

        allocate(Me%FatherSonZ(     SizeSon%ILB:SizeSon%IUB, &
                                    SizeSon%JLB:SizeSon%JUB, &
                                    SizeSon%KLB:SizeSon%KUB))    

        allocate(Me%InsideFatherCell(   SizeSon%ILB:SizeSon%IUB, &
                                        SizeSon%JLB:SizeSon%JUB, &
                                        SizeSon%KLB:SizeSon%KUB))
                                        
        call PopulateArrays    

        STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

    end subroutine Construct_ConvoConservative

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine Construct_ConvoNCStructured (HorizontalGridSonID,                          &
                                          HorizontalMapSonID,                          &                                          
                                          GeometrySonID,                                &
                                          MapSonID,                                     &
                                          HorizontalGridFatherID,                       &
                                          HorizontalMapFatherID,                       &
                                          GeometryFatherID,                             &
                                          MapFatherID,                                  &
                                          STAT)

        !Arguments-----------------------------------------------------------
        integer                                         ::  HorizontalGridSonID,                          &
                                                            HorizontalMapSonID,                          &
                                                            GeometrySonID,                                &
                                                            MapSonID,                                     &
                                                            HorizontalGridFatherID,                       &
                                                            HorizontalMapFatherID,                       &
                                                            GeometryFatherID,                             &
                                                            MapFatherID                                  
        integer, optional, intent(OUT)                  ::  STAT

        !Locals--------------------------------------------------------------
        integer                                         ::  STAT_

        STAT_ = UNKNOWN_

        select case (Me%ComputeOptions%NC_Type)

            case (NC_User_)

                call AllocateExternal_ConvoStruct(STAT_)

                call ConstructGrid(Me%ExternalVar%FatherGrid,                           &
                                                    HorizontalGridFatherID,                         &
                                                    HorizontalMapFatherID,                         &
                                                    GeometryFatherID,                               &
                                                    MapFatherID, STAT_)

                call ConstructGrid(Me%ExternalVar%SonGrid,                           &
                                                    HorizontalGridSonID,                         &
                                                    HorizontalMapSonID,                         &
                                                    GeometrySonID,                               &
                                                    MapSonID, STAT_)

                call ConstructNC_User3D(Me%NC, Me%ExternalVar%SonGrid%WorkSize3D, STAT_)

            case (NC_Smoothing_)

                call AllocateExternal_ConvoStruct(STAT_)

                call ConstructGrid(Me%ExternalVar%FatherGrid,                           &
                                                    HorizontalGridFatherID,                         &
                                                    HorizontalMapFatherID,                         &
                                                    GeometryFatherID,                               &
                                                    MapFatherID, STAT_)

                call ConstructGrid(Me%ExternalVar%SonGrid,                           &
                                                    HorizontalGridSonID,                         &
                                                    HorizontalMapSonID,                         &
                                                    GeometrySonID,                               &
                                                    MapSonID, STAT_)

                call ConstructNC_Smoothing3D(Me%NC, Me%ExternalVar%SonGrid%WorkSize3D, STAT_)

            case (NC_Data_)

                call AllocateExternal_ConvoStruct(STAT_)

                call ConstructGrid(Me%ExternalVar%FatherGrid,                           &
                                                    HorizontalGridFatherID,                         &
                                                    HorizontalMapFatherID,                         &
                                                    GeometryFatherID,                               &
                                                    MapFatherID, STAT_)

                call ConstructGrid(Me%ExternalVar%SonGrid,                           &
                                                    HorizontalGridSonID,                         &
                                                    HorizontalMapSonID,                         &
                                                    GeometrySonID,                               &
                                                    MapSonID, STAT_)

                call ConstructNC_Data3D(Me%NC, Me%ExternalVar%SonGrid%WorkSize3D, STAT_)

            end select

        if (present(STAT)) STAT = STAT_

    end subroutine Construct_ConvoNCStructured

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine Construct_ConvoNCNonStructured (HorizontalGridSonID,                     &
                                               HorizontalMapSonID,                      &
                                               GeometrySonID,                           &
                                               MapSonID,                                &
                                               TXYZP_Points,                            &
                                               STAT)

        !Arguments-----------------------------------------------------------
        integer                                         ::  HorizontalGridSonID,        &
                                                            HorizontalMapSonID,         &
                                                            GeometrySonID,              &
                                                            MapSonID
        type(T_XYZpoints), pointer                      ::  TXYZP_Points
        integer, optional, intent(OUT)                  ::  STAT

        !Locals--------------------------------------------------------------
        integer                                         ::  STAT_

        STAT_ = UNKNOWN_

        select case (Me%ComputeOptions%NC_Type)

            case (NC_User_)

                call AllocateExternal_ConvoNonStruct(STAT_)

                call AllocateTXYZP_points(Me%ExternalVar%TXYZP_points,TXYZP_Points,STAT_)

                call ConstructGrid(Me%ExternalVar%SonGrid,                           &
                                                    HorizontalGridSonID,                         &
                                                    HorizontalMapSonID,                         &
                                                    GeometrySonID,                               &
                                                    MapSonID, STAT_)

                call ConstructNC_User3D(Me%NC, Me%ExternalVar%SonGrid%WorkSize3D, STAT_)

            case (NC_Smoothing_)

                call AllocateExternal_ConvoNonStruct(STAT_)

                call AllocateTXYZP_points(Me%ExternalVar%TXYZP_points,TXYZP_Points,STAT_)

                call ConstructGrid(Me%ExternalVar%SonGrid,                           &
                                                    HorizontalGridSonID,                         &
                                                    HorizontalMapSonID,                         &
                                                    GeometrySonID,                               &
                                                    MapSonID, STAT_)

                call ConstructNC_Smoothing3D(Me%NC, Me%ExternalVar%SonGrid%WorkSize3D, STAT_)

            case (NC_Data_)

                call AllocateExternal_ConvoNonStruct(STAT_)
                
                call AllocateTXYZP_points(Me%ExternalVar%TXYZP_points,TXYZP_Points,STAT_)

                call ConstructPoints2Grid_3D(STAT_)

                call ConstructGrid(Me%ExternalVar%SonGrid,                           &
                                                    HorizontalGridSonID,                         &
                                                    HorizontalMapSonID,                         &
                                                    GeometrySonID,                               &
                                                    MapSonID, STAT_)

                call ConstructNC_Data3D(Me%NC, Me%ExternalVar%SonGrid%WorkSize3D, STAT_)

            end select

        if (present(STAT)) STAT = STAT_

    end subroutine Construct_ConvoNCNonStructured

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine AllocateExternal_ConvoStruct(STAT)

        !Argument---------------------------------------------------------------
        integer, optional, intent(OUT)                         :: STAT

        !Local---------------------------------------------------------------
        integer                                                 :: STAT_

        !Begin Shorten variables names-------------------------------------------

        STAT_ = UNKNOWN_
        
        !Allocates the external variables
        call AllocateExternal(STAT_)
        call AllocateGrid(Me%ExternalVar%FatherGrid, STAT_) 
        call AllocateGrid(Me%ExternalVar%SonGrid, STAT_)
        call AllocateField(Me%ExternalVar%FatherField, STAT_)
        call AllocateField(Me%ExternalVar%SonField, STAT_)

        if (present(STAT)) STAT = STAT_

    end subroutine AllocateExternal_ConvoStruct

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine AllocateExternal_ConvoNonStruct(STAT)

        !Argument---------------------------------------------------------------
        integer, optional, intent(OUT)                         :: STAT

        !Local---------------------------------------------------------------
        integer                                                 :: STAT_

        !Begin Shorten variables names-------------------------------------------

        STAT_ = UNKNOWN_
        
        !Allocates the external variables
        call AllocateExternal(STAT_)
        call AllocateGrid(Me%ExternalVar%SonGrid, STAT_)
        call AllocateField(Me%ExternalVar%SonField, STAT_)

        if (present(STAT)) STAT = STAT_

    end subroutine AllocateExternal_ConvoNonStruct

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine AllocateTXYZP_points(TXYZP,TXYZP_Points, STAT)

        !Argument---------------------------------------------------------------
        Type (T_XYZpoints), pointer                     :: TXYZP, TXYZP_points
        integer, optional, intent(OUT)                  :: STAT

        !Local---------------------------------------------------------------
        integer                                         :: STAT_
 
        STAT_ = UNKNOWN_

        if(associated(TXYZP_Points)) then
        
            TXYZP => TXYZP_Points
            STAT_ = SUCCESS_
        
        end if

        if(present(STAT)) STAT = STAT_

    end subroutine AllocateTXYZP_points

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine AllocateNC(NC, STAT)

        !Argument---------------------------------------------------------------
        Type (T_NC), pointer                            :: NC
        integer, optional, intent(OUT)                  :: STAT

        !Local---------------------------------------------------------------
        Type (T_NC), pointer                            :: NewNC
        integer                                         :: STAT_
 
        STAT_ = UNKNOWN_

        allocate(NewNC)

        NC => NewNC

        if(associated(NC)) STAT_ = SUCCESS_

        if(present(STAT)) STAT = STAT_

    end subroutine AllocateNC

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ConstructNC_User3D(NC, Size, STAT)

        !Argument---------------------------------------------------------------
        Type (T_NC), pointer                            :: NC
        Type (T_Size3D), pointer                        :: Size
        integer, optional, intent(OUT)                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_

        STAT_ = UNKNOWN_

        Allocate(NC%Tau(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        if (associated(NC%Tau)) STAT_ = SUCCESS_

        if(present(STAT)) STAT = STAT_

    end subroutine ConstructNC_User3D

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ConstructNC_Smoothing3D(NC, Size, STAT)

        !Argument---------------------------------------------------------------
        Type (T_NC), pointer                            :: NC
        Type (T_Size3D), pointer                        :: Size
        integer, optional, intent(OUT)                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: aux
        integer                                         :: n_dim
        integer                                         :: STAT_

        STAT_ = UNKNOWN_

        !Begin shorten variable names-------------------------------------------
        n_dim = Me%ComputeOptions%n_dimensions
        
        aux = (Size%IUB - Size%ILB + 1)*(Size%JUB - Size%JLB + 1)*(Size%KUB - Size%KLB + 1)

        Allocate(NC%Tau(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        Allocate(NC%TauR(aux*n_dim, n_dim))
        Allocate(NC%TauM(n_dim))

        STAT_ = SUCCESS_

        if(present(STAT)) STAT = STAT_

    end subroutine ConstructNC_Smoothing3D

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ConstructNC_Data3D(NC, Size, STAT)

        !Argument---------------------------------------------------------------
        Type (T_NC), pointer                            :: NC
        Type (T_Size3D), pointer                        :: Size
        integer, optional, intent(OUT)                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: aux
        integer                                         :: n_dim
        integer                                         :: n_groups, SampleSize
        integer                                         :: STAT_

        STAT_ = UNKNOWN_

        !Begin shorten variable names-------------------------------------------
        n_dim = Me%ComputeOptions%n_dimensions
        n_groups = Me%ComputeOptions%n_groups
        SampleSize = Me%ComputeOptions%sample_size
        
        aux = (Size%IUB - Size%ILB + 1)*(Size%JUB - Size%JLB + 1)*(Size%KUB - Size%KLB + 1)

        Allocate(NC%GroupCenterX(n_groups**n_dim))
        Allocate(NC%GroupCenterY(n_groups**n_dim))
        Allocate(NC%GroupCenterZ(n_groups**n_dim))
        Allocate(NC%Tau(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        Allocate(NC%TauR(n_groups**n_dim, n_dim))
        Allocate(NC%TauM(n_dim))
        Allocate(NC%Response(n_groups**n_dim, n_dim, SampleSize))
        Allocate(NC%Predictor(n_groups**n_dim, n_dim, SampleSize))

        STAT_ = SUCCESS_

        if(present(STAT)) STAT = STAT_

    end subroutine ConstructNC_Data3D

    !--------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    
    subroutine PopulateArrays

        !Local --------------------------------------------------------------
        integer                             :: i, j, k, h, iub, jub, kub, n_odd
        logical                             :: odd(5)
        integer, pointer, dimension(:,:)    :: WaterPoints2D
        integer, pointer, dimension(:,:,:)  :: MatrixUpdate3D, OpenPoints3D

        !Begin shorten variables---------------------------------------------
        iub             = Me%ExternalVar%FatherGrid%WorkSize3D%IUB
        jub             = Me%ExternalVar%FatherGrid%WorkSize3D%JUB
        kub             = Me%ExternalVar%FatherGrid%WorkSize3D%KUB
        WaterPoints2D   => Me%ExternalVar%FatherGrid%WaterPoints2D
        OpenPoints3D    => Me%ExternalVar%FatherGrid%OpenPoints3D
        MatrixUpdate3D  => Me%MatrixUpdate3D
    
        !populating MatrixUpdate3D
        odd(1)=.true.
        do i=1,iub
            odd(2)=IsOdd(i)
            do j=1,jub
                if(WaterPoints2D(i,j).eq.0) then
                    MatrixUpdate3D(i,j,:)=0
                else
                    odd(3)=IsOdd(j)
                    do k=1,kub
                        if(OpenPoints3D(i,j,k).eq.0) then
                            MatrixUpdate3D(i,j,k)=0
                        else
                            odd(4)=IsOdd(k)
                            n_odd=0
                            do h=2,4
                                if(odd(h)) n_odd=n_odd+1
                            enddo
                            odd(5)=IsOdd(n_odd)
                            if (odd(1).eqv.odd(5)) then
                                MatrixUpdate3D(i,j,k)=1
                            else
                                MatrixUpdate3D(i,j,k)=0
                            endif
                        endif
                    enddo
                endif
            enddo
        enddo

    end subroutine PopulateArrays

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadKeywordOptions(EnterDataID, FromWhere) 

        !Arguments-------------------------------------------------------------
        integer                                     :: EnterDataID, FromWhere 

        !Local-----------------------------------------------------------------
        integer                                     :: iflag, STAT_CALL, nUsers

        !Begin-----------------------------------------------------------------

        
        Me%ObjEnterData = AssociateInstance(mENTERDATA_, EnterDataID)

        !The number of dimensions
        call GetData(Me%ComputeOptions%n_dimensions,                                    &
                Me%ObjEnterData, iflag,                                                 &
                SearchType   = FromWhere,                                               &
                keyword      = 'N_DIM',                                                 &
                default      = 3,                                                       &
                ClientModule = 'ModuleInterpolation',                                   &
                STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR10'

        !The methodology use in the interpolation process
        call GetData(Me%ComputeOptions%Methodology,                                     &
                Me%ObjEnterData, iflag,                                                 &
                SearchType   = FromWhere,                                               &
                keyword      = 'METHODOLOGY',                                           &
                default      = ConvolutionConserv_,                                               &
!                default      = Bilinear_,                                               &
                ClientModule = 'ModuleInterpolation',                                   &
                STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR20'

        if (Me%ComputeOptions%Methodology == ConvolutionNonConserv_ .or.                &
            Me%ComputeOptions%Methodology == ConvolutionConserv_) then

            Me%ConvolutionApproach = .true.     !Algoritmos do Ricardo

        else

            Me%ConvolutionApproach = .false.    !Algoritmos classicos

        endif

        if (.not. Me%ConvolutionApproach) then

            !The degree of the polinomial use for vertical 1D interpolations
            call GetData(Me%ComputeOptions%PoliDegreeVert,                              &
                    Me%ObjEnterData, iflag,                                             &
                    SearchType   = FromWhere,                                           &
                    keyword      = 'POLI_DEGREE_VERT',                                  &
                    default      = LinearProfile_,                                      &
                    ClientModule = 'ModuleInterpolation',                               &
                    STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR30'

            if (Me%ComputeOptions%PoliDegreeVert > 10)                                  &
                stop 'ReadKeywordOptions - ModuleInterpolation - ERR31'

            if (Me%ComputeOptions%PoliDegreeVert < 1)                                   &
                stop 'ReadKeywordOptions - ModuleInterpolation - ERR32'

            !Chek if the user wants to extrapolate in the vertical 
            call GetData(Me%ComputeOptions%ExtrapolateProfile,                          &
                    Me%ObjEnterData, iflag,                                             &
                    SearchType   = FromWhere,                                           &
                    keyword      = 'EXTRAPOLATE_PROFILE',                               &
                    default      = .true.,                                              &
                    ClientModule = 'ModuleInterpolation',                               &
                    STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR33'


            !Chek if the user wants to extrapolate in the vertical 
            call GetData(Me%ComputeOptions%Extrapolate2D,                               &
                    Me%ObjEnterData, iflag,                                             &
                    SearchType   = FromWhere,                                           &
                    keyword      = 'EXTRAPOLATE_2D',                                    &
                    default      = NotActive_,                                          &
                    ClientModule = 'ModuleInterpolation',                               &
                    STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR34'

        endif


        !Where points are defined in the cell (Z - center; U - Face U; V : Face V) 
        call GetData(Me%ComputeOptions%TypeZUV,                                         &
                Me%ObjEnterData, iflag,                                                 &
                SearchType   = FromWhere,                                               &
                keyword      = 'TYPE_ZUV',                                              &
                default      = TypeZ_,                                                  &
                ClientModule = 'ModuleInterpolation',                                   &
                STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR40'

        Me%OkZ = .false.; Me%OkU  = .false.;  Me%OkV = .false.;

        if (Me%ComputeOptions%TypeZUV == TypeZ_) Me%OkZ = .true.
        if (Me%ComputeOptions%TypeZUV == TypeU_) Me%OkU = .true.
        if (Me%ComputeOptions%TypeZUV == TypeV_) Me%OkV = .true.

inv:    if (Me%ComputeOptions%Methodology == InverseWeight_) then

            call GetData(Me%ComputeOptions%MaxDistance,                                 &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromWhere,                                      &
                         keyword      = 'MAX_DISTANCE',                                 &
                         ClientModule = 'ModuleInterpolation',                          &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR50'

            if (iflag == 0) then
                write(*,*)'Max Distance not given'
                write(*,*)'Use Keyword MAX_DISTANCE'
                stop 'ReadKeywordOptions - ModuleInterpolation - ERR60'
            endif
            
            call GetData(Me%ComputeOptions%IWDn,                                        &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromWhere,                                      &
                         keyword      = 'IWD_N',                                        &
                         default      = 2.0,                                            &
                         ClientModule = 'ModuleInterpolation',                          &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR70'

        endif inv

        if (Me%ConvolutionApproach) then

            call GetData(Me%ComputeOptions%KernelType,                                  &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromWhere,                                      &
                         keyword      = 'KERNEL_TYPE',                                  &
                         default      = CharGaussian_,                                  &
                         ClientModule = 'ModuleInterpolation',                          &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR80'

            if(Me%ComputeOptions%Methodology == 2) then
                Me%NonConservative = .true.
            else
                Me%NonConservative = .false.
            end if

            if (Me%NonConservative) then

                call GetData(Me%ComputeOptions%NC_Type,                                  &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromWhere,                                      &
                             keyword      = 'NC_TYPE',                                  &
                             default      = NC_Smoothing_,                                  &
                             ClientModule = 'ModuleInterpolation',                          &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR80'

                if ( Me%ComputeOptions%NC_Type == NC_Smoothing_) then

                    call GetData(Me%ComputeOptions%phi,                                  &
                                 Me%ObjEnterData, iflag,                                        &
                                 SearchType   = FromWhere,                                      &
                                 keyword      = 'PHI',                                  &
                                 default      = Default_Phi_,                                  &
                                 ClientModule = 'ModuleInterpolation',                          &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR80'

                end if

                if ( Me%ComputeOptions%NC_Type == NC_Data_) then

                    call GetData(Me%ComputeOptions%n_groups,                                  &
                                 Me%ObjEnterData, iflag,                                        &
                                 SearchType   = FromWhere,                                      &
                                 keyword      = 'N_GROUPS',                                  &
                                 default      = 1,                                  &
                                 ClientModule = 'ModuleInterpolation',                          &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR80'

                    call GetData(Me%ComputeOptions%sample_size,                                  &
                                 Me%ObjEnterData, iflag,                                        &
                                 SearchType   = FromWhere,                                      &
                                 keyword      = 'SAMPLE_SIZE',                                  &
                                 default      = SampleSize_,                                  &
                                 ClientModule = 'ModuleInterpolation',                          &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR80'

                    call GetData(Me%ComputeOptions%max_iterations,                                  &
                                 Me%ObjEnterData, iflag,                                        &
                                 SearchType   = FromWhere,                                      &
                                 keyword      = 'MAX_ITERATIONS',                                  &
                                 default      = MaxIterations_,                                  &
                                 ClientModule = 'ModuleInterpolation',                          &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadKeywordOptions - ModuleInterpolation - ERR80'

                end if

            end if

        endif


        nUsers = DeAssociateInstance (mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'ReadKeywordOptions - ModuleInterpolation - ERR300'

    end subroutine ReadKeywordOptions

    !----------------------------------------------------------------------------

    !--------------------------------------------------------------------------


    subroutine Relate_FatherSon_3D

        !Arguments---------------------------------------------------------

        !External----------------------------------------------------------

        !Internal----------------------------------------------------------
        integer                              :: i,j,k
        integer                              :: lb(3),ub(3),lbF(3),ubF(3)
        integer                              :: PossibleF(3), res(3)
        real, pointer, dimension (:,:)          :: XX_IE, YY_IE
        real, pointer, dimension (:,:)          :: CenterX, CenterY
        real, pointer, dimension (:,:,:)        :: SZZ
        real, pointer, dimension (:,:,:)        :: CenterZ
        integer, pointer, dimension (:,:)    :: FatherSonX, FatherSonY
        integer, pointer, dimension (:,:,:)  :: FatherSonZ

        ! Begin shorten variable names-------------------------------------

        lb(1)           = Me%ExternalVar%SonGrid%WorkSize3D%ILB
        lb(2)           = Me%ExternalVar%SonGrid%WorkSize3D%JLB
        lb(3)           = Me%ExternalVar%SonGrid%WorkSize3D%KLB
        ub(1)           = Me%ExternalVar%SonGrid%WorkSize3D%IUB
        ub(2)           = Me%ExternalVar%SonGrid%WorkSize3D%JUB
        ub(3)           = Me%ExternalVar%SonGrid%WorkSize3D%KUB
        CenterX         => Me%ExternalVar%SonGrid%CenterX
        CenterY         => Me%ExternalVar%SonGrid%CenterY
        CenterZ         => Me%ExternalVar%SonGrid%CenterZ
        FatherSonX      => Me%FatherSonX
        FatherSonY      => Me%FatherSonY
        FatherSonZ      => Me%FatherSonZ
        XX_IE            => Me%ExternalVar%FatherGrid%XX_IE
        YY_IE            => Me%ExternalVar%FatherGrid%YY_IE
        SZZ             => Me%ExternalVar%FatherGrid%SZZ

        ! Begin -----------------------------------------------------------    
    
        call Corners_3D(lbF,ubF)
    
        PossibleF(:)=lbF(:)

        do i=lb(1),ub(1)
        do j=lb(2),ub(2)
        do k=lb(3),ub(3)

            if (YY_IE(PossibleF(1),j).gt.CenterY(i,j).or.      &
                XX_IE(i,PossibleF(2)).gt.CenterX(i,j).or.      &
                SZZ(i,j,PossibleF(3)).gt.CenterZ(i,j,k).or.   &
                YY_IE(PossibleF(1)+1,j).lt.CenterY(i,j).or.    &
                XX_IE(i,PossibleF(2)+1).lt.CenterX(i,j).or.    &
                SZZ(i,j,PossibleF(3)+1).lt.CenterZ(i,j,k)) then

                call Gold_3D(i,j,k,lbF,ubF,res, clue = PossibleF)
                PossibleF(:)=res(:)

            endif
            
            FatherSonY(i,j)=PossibleF(1)
            FatherSonX(i,j)=PossibleF(2)
            FatherSonZ(i,j,k)=PossibleF(3)

        enddo
        enddo
        enddo

    end subroutine Relate_FatherSon_3D
    
        !------------------------------------------------------------------

        !------------------------------------------------------------------

    subroutine Corners_3D(lbF,ubF)
    
        !Arguments---------------------------------------------------------

        !External----------------------------------------------------------
        integer                          :: lbF(3),ubF(3)

        !Internal----------------------------------------------------------
        integer                          :: a,b,c
        integer                          :: h,i,j,k
        integer                          :: corner(3)
        integer                          :: lb(3),ub(3)

        ! Begin shorten variable names-------------------------------------

        lb(1)           = Me%ExternalVar%FatherGrid%WorkSize3D%ILB
        lb(2)           = Me%ExternalVar%FatherGrid%WorkSize3D%JLB
        lb(3)           = Me%ExternalVar%FatherGrid%WorkSize3D%KLB
        ub(1)           = Me%ExternalVar%FatherGrid%WorkSize3D%IUB
        ub(2)           = Me%ExternalVar%FatherGrid%WorkSize3D%JUB
        ub(3)           = Me%ExternalVar%FatherGrid%WorkSize3D%KUB

        ! Begin -----------------------------------------------------------    
        
        lbF(:)=ub(:)
        ubF(:)=lb(:)

        do a=0,1
            i=lb(1)+a*(ub(1)-lb(1))
            do b=0,1
                j=lb(2)+b*(ub(2)-lb(2))
                do c=0,1
                    k=lb(3)+c*(ub(3)-lb(3))
                    call Gold_3D(i,j,k,lb,ub,corner)
                    do h=1,3
                        if (corner(h).lt.lbF(h)) lbF(h)=corner(h)
                        if (corner(h).gt.ubF(h)) lbF(h)=corner(h)
                    enddo
                enddo
            enddo
        enddo

    end subroutine Corners_3D

   !------------------------------------------------------------------

   !------------------------------------------------------------------

    subroutine Gold_3D(i,j,k,lb,ub,search,clue)

        !Arguments---------------------------------------------------------
        integer                          :: i,j,k
        integer                          :: lb(3),ub(3)
        integer, optional                :: clue(3)

        !External----------------------------------------------------------
        integer                          :: search(3)

        !Internal----------------------------------------------------------
        real(8)                             :: GoldPhi
        real(8)                             :: alphaL(3), alphaH(3)
        logical                             :: NotReached(3)
        real, pointer, dimension (:,:)      :: XX_IE, YY_IE
        real, pointer, dimension (:,:)      :: CenterX, CenterY
        real, pointer, dimension (:,:,:)    :: SZZ
        real, pointer, dimension (:,:,:)    :: CenterZ

        ! Begin shorten variable names-------------------------------------

        CenterX         => Me%ExternalVar%SonGrid%CenterX
        CenterY         => Me%ExternalVar%SonGrid%CenterY
        CenterZ         => Me%ExternalVar%SonGrid%CenterZ
        XX_IE            => Me%ExternalVar%SonGrid%XX_IE
        YY_IE            => Me%ExternalVar%SonGrid%YY_IE
        SZZ             => Me%ExternalVar%SonGrid%SZZ
        
        ! Begin -----------------------------------------------------------

        GoldPhi=(Sqrt(5d0)-1d0)/2d0
        alphaL(:)=real(lb(:))
        alphaH(:)=real(ub(:))

        if(present(clue)) then
            search(:)=clue(:)
            else
            search(:)=int(GoldPhi*alphaL(:)+(1d0-GoldPhi)*alphaH(:))
        endif

        do while(NotReached(1))
            if(YY_IE(search(1),j).lt.CenterY(i,j)) then
                if (YY_IE(search(1)+1,j).lt.CenterY(i,j)) then
                    alphaL(1)=real(search(1)+1)
                    search(1)=int(GoldPhi*alphaL(1)+(1d0-GoldPhi)*alphaH(1))
                    else
                    NotReached(1)=.false.
                endif
            else
                alphaH(1)=real(search(1))
                search(1)=int((1d0-GoldPhi)*alphaL(1)+GoldPhi*alphaH(1))
            endif
        enddo

        do while(NotReached(2))
            if(XX_IE(i,search(2)).lt.CenterX(i,j)) then
                if (XX_IE(i,search(2)+1).lt.CenterX(i,j)) then
                    alphaL(2)=real(search(2)+1)
                    search(2)=int(GoldPhi*alphaL(2)+(1d0-GoldPhi)*alphaH(2))
                    else
                    NotReached(2)=.false.
                endif
            else
                alphaH(2)=real(search(2))
                search(2)=int((1d0-GoldPhi)*alphaL(2)+GoldPhi*alphaH(2))
            endif
        enddo

        do while(NotReached(2))        
            if(SZZ(i,j,search(3)).lt.CenterZ(i,j,k)) then
                if (SZZ(i,j,search(3)+1).lt.CenterZ(i,j,k)) then
                    alphaL(3)=real(search(3)+1)
                    search(3)=int(GoldPhi*alphaL(3)+(1d0-GoldPhi)*alphaH(3))
                    else
                    NotReached(3)=.false.
                endif
            else
                alphaH(3)=real(search(3))
                search(3)=int((1d0-GoldPhi)*alphaL(3)+GoldPhi*alphaH(3))
            endif

        enddo

    end subroutine Gold_3D

   !------------------------------------------------------------------


   !------------------------------------------------------------------

    subroutine InFatherCell_3D(i,j,k)
 
        !Arguments---------------------------------------------------------------
        integer(8)                              :: i,j,k

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                              :: n_dimensions
        real(8)                                 :: angle(3)
        real(8)                                 :: X, Y, Z, XC, YC, ZC
        real(8)                                 :: X_CRITICAl, Y_CRITICAL
        real, pointer, dimension(:,:)           :: DX, DY, XX_IE, YY_IE
        real, pointer, dimension(:,:,:)         :: DZ, SZZ
        character(3), pointer, dimension(:,:,:) :: InsideFatherCell


        !Begin shorten variables names-------------------------------------------
        
        n_dimensions    = Me%ComputeOptions%n_dimensions
        Y               = Me%ExternalVar%SonGrid%CenterY(i,j)
        X               = Me%ExternalVar%SonGrid%CenterX(i,j)
        Z               = Me%ExternalVar%SonGrid%CenterZ(i,j,k)
        YC              = Me%ExternalVar%FatherGrid%CenterY(i,j)
        XC              = Me%ExternalVar%FatherGrid%CenterX(i,j)
        ZC              = Me%ExternalVar%FatherGrid%CenterZ(i,j,k)

        XX_IE                => Me%ExternalVar%FatherGrid%XX_IE
        YY_IE                => Me%ExternalVar%FatherGrid%YY_IE
        DX                  => Me%ExternalVar%FatherGrid%DX
        DY                  => Me%ExternalVar%FatherGrid%DY
        DZ                  => Me%ExternalVar%FatherGrid%DZ
        SZZ                 => Me%ExternalVar%FatherGrid%SZZ
        InsideFatherCell    => Me%InsideFatherCell

        !Begin-------------------------------------------------------------------

        select case(n_dimensions)

        case(1)
            if(Y.lt.YC) then
                InsideFatherCell(i,j,k)="Dl1"
                else
                InsideFatherCell(i,j,k)="Du1"
            endif
        case(2)
            angle(1)=asin((Y-YC)/(((X-XC)**2d0+(Y-YC)**2d0)**0.5d0))
            angle(2)=asin((YY_IE(i+1,j+1)-YC)/(((XX_IE(i+1,j+1)-XC)**2d0        &
                +(YY_IE(i+1,j+1)-YC)**2d0)**0.5d0))
            
            if (X.gt.XC) then
                if(abs(angle(1)).lt.angle(2)) then
                    InsideFatherCell(i,j,k)="Au2"
                else
                    if(Y.gt.YC) then
                        InsideFatherCell(i,j,k)="Au1"
                    else
                        InsideFatherCell(i,j,k)="Al2"
                    endif
                endif
            else
                if(abs(angle(1)).lt.angle(2)) then
                    InsideFatherCell(i,j,k)="Al2"
                else
                    if(Y.gt.YC) then
                        InsideFatherCell(i,j,k)="Au1"
                    else
                        InsideFatherCell(i,j,k)="Al1"
                    endif
                endif
            endif

        case(3)
                
            ZC=(SZZ(i,j,k-1)+SZZ(i,j,k))/2d0

            Y_CRITICAL=YC+abs(Z-ZC)*DY(i,j)/DZ(i,j,k)
            X_CRITICAL=XC+abs(Z-ZC)*DX(i,j)/DZ(i,j,k)
        
            if (abs(X-XC).lt.X_CRITICAL.and.abs(Y-YC).lt.Y_CRITICAL) then
                if(Z.lt.ZC) then
                    InsideFatherCell(i,j,k)="Vl3"
                else
                    InsideFatherCell(i,j,k)="Vu3"
                endif
            else
                angle(1)=asin((Y-YC)/(((X-XC)**2d0+(Y-YC)**2d0)**0.5d0))
                angle(2)=asin((YY_IE(i+1,j+1)-YC)/(((XX_IE(2,2)-XC)**2d0    &
                    +(YY_IE(2,2)-YC)**2d0)**0.5d0))
                if (X.gt.XC) then
                    if(abs(angle(1)).lt.angle(2)) then
                        InsideFatherCell(i,j,k)="Vu2"
                    else
                        if(Y.gt.YC) then
                            InsideFatherCell(i,j,k)="Vu1"
                        else !Al1
                            InsideFatherCell(i,j,k)="Vl1"
                        endif
                    endif
                endif
            endif
        
        case default
            stop 'InFatherCell_3D - ModuleInterpolation - ERR01' 

        end select
                
    end subroutine InFatherCell_3D

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ConstructGrid(GridPointer,HorizontalGridID,  HorizontalMapID,           &
                             GeometryID, MapID, STAT)

        !Arguments---------------------------------------------------------------
        type (T_Grid), pointer                          :: GridPointer
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: GeometryID
        integer                                         :: MapID
        integer, optional, intent(OUT)                  :: STAT

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        STAT_ = UNKNOWN_
    
        GridPointer%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
        GridPointer%ObjGeometry       = AssociateInstance (mGEOMETRY_, GeometryID)
        GridPointer%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_, HorizontalMapID)
        GridPointer%ObjMap            = AssociateInstance (mMAP_, MapID)

        call ReadLockExternalVar(GridPointer, STAT_)

        ! Construct the variables common to all modules  
        call PopulateParameters(GridPointer, STAT_)

        if (present(STAT)) STAT = STAT_

    end subroutine ConstructGrid

    !--------------------------------------------------------------------------
    subroutine ConstructCLocalFatherGrid (TimeSonID,       HorizontalGridSonID,          &
                                         BathymetrySonID, HorizontalMapSonID)

        !Arguments---------------------------------------------------------------
        integer                                 :: TimeSonID
        integer                                 :: HorizontalGridSonID
        integer                                 :: BathymetrySonID
        integer                                 :: HorizontalMapSonID

        !Local-------------------------------------------------------------
        real,       dimension(:,:  ), pointer   :: SurfaceElevation
        type(T_Time), pointer                   :: Time
        integer                                 :: STAT_CALL
        character(Len=StringLength)             :: FatherGeometryFile

        !Begin-----------------------------------------------------------------


        write(*,*)'Constructing local FatherGrid grid...'
        
        Me%LocalFatherGrid%ObjTime           = AssociateInstance (mTIME_,           TimeSonID)
        Me%LocalFatherGrid%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridSonID)
        Me%LocalFatherGrid%ObjBathymetry     = AssociateInstance (mGRIDDATA_,       BathymetrySonID) 
        Me%LocalFatherGrid%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapSonID )    
        

        call GetGeometryInputFile(Me%ExternalVar%FatherGrid%ObjGeometry, FatherGeometryFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructCLocalFatherGrid - ModuleInterpolation - ERR010'

        call GetComputeCurrentTime(Me%LocalFatherGrid%ObjTime, Time, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructCLocalFatherGrid - ModuleInterpolation - ERR011'

        call ConstructGeometry      (GeometryID       = Me%LocalFatherGrid%ObjGeometry,             &
                                     GridDataID       = Me%LocalFatherGrid%ObjBathymetry,           &
                                     HorizontalGridID = Me%LocalFatherGrid%ObjHorizontalGrid,       &
                                     HorizontalMapID  = Me%LocalFatherGrid%ObjHorizontalMap,        &
                                     NewDomain        = FatherGeometryFile,                     &
                                     ActualTime       = Time,                     &
                                     STAT             = STAT_CALL)

        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructCLocalFatherGrid -  ModuleInterpolateGrids - ERR20'

        call GetGeometrySize(GeometryID     = Me%LocalFatherGrid%ObjGeometry,                       &
                             Size           = Me%LocalFatherGrid%Size3D,                            &
                             WorkSize       = Me%LocalFatherGrid%WorkSize3D,                        &
                             STAT           = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructCLocalFatherGrid -  ModuleInterpolateGrids - ERR30'

        call ConstructMap ( Map_ID          = Me%LocalFatherGrid%ObjMap,                            &
                            GeometryID      = Me%LocalFatherGrid%ObjGeometry,                       &
                            HorizontalMapID = Me%LocalFatherGrid%ObjHorizontalMap,                  &
                            TimeID          = Me%LocalFatherGrid%ObjTime,                           &
                            STAT            = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructCLocalFatherGrid -  ModuleInterpolateGrids - ERR40'

        allocate(SurfaceElevation (Me%LocalFatherGrid%WorkSize3D%ILB:Me%LocalFatherGrid%WorkSize3D%IUB,         &
                                   Me%LocalFatherGrid%WorkSize3D%JLB:Me%LocalFatherGrid%WorkSize3D%JUB))
        SurfaceElevation(:,:) = 0.

        call GetWaterPoints3D(Me%LocalFatherGrid%ObjMap, Me%LocalFatherGrid%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ConstructCLocalFatherGrid - ModuleInterpolateGrids - ERR50'

        call ComputeInitialGeometry(GeometryID       = Me%LocalFatherGrid%ObjGeometry,              &
                                    WaterPoints3D    = Me%LocalFatherGrid%WaterPoints3D,            &
                                    SurfaceElevation = SurfaceElevation,            &
                                    ActualTime       = Time,                        &
                                    STAT             = STAT_CALL )
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructCLocalFatherGrid -  ModuleInterpolateGrids - ERR60'

        call UnGetMap(Me%LocalFatherGrid%ObjMap, Me%LocalFatherGrid%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ConstructCLocalFatherGrid - ModuleInterpolateGrids - ERR70'

        deallocate(SurfaceElevation)


    end subroutine ConstructCLocalFatherGrid

    !----------------------------------------------------------------------------

    subroutine ConstructPoints2Grid_3D(STAT)

        !Arguments-------------------------------------------------------------------
        integer, optional, intent(OUT)           :: STAT

        !Local-----------------------------------------------------------------------
        integer                                  :: n_dimensions, n_obs
        integer                                  :: dim_X,dim_Y,dim_Z
        integer                                  :: new_dimX,new_dimY,new_dimZ
        integer                                  :: work_X,work_Y,work_Z
        integer                                  :: a,i,j,k
        real,pointer,dimension(:)                   :: x,y,z,w
        real,allocatable,dimension(:)               :: xx,yy,zz
        logical,allocatable,dimension(:,:,:)        :: CellHasPoint,new_CHP
        real(8),allocatable,dimension(:,:,:)        :: CellValue,new_CV
        real(8),allocatable,dimension(:)            :: Cell_XBounds,Cell_YBounds,Cell_ZBounds
        real(8),allocatable,dimension(:)            :: new_XBounds,new_YBounds,new_ZBounds
        integer,allocatable,dimension(:)         :: originX,originY,originZ
        integer,allocatable,dimension(:)         :: destinationX,destinationY,destinationZ

        logical                                     :: ok_X_merge, ok_Y_merge, ok_Z_merge
        logical                                     :: not_finished
        logical                                     :: warning_X, warning_Y,warning_Z
        logical, allocatable, dimension(:)          :: merge_X, merge_Y, merge_Z
        integer                                     :: STAT_CALL, STAT_
        
        STAT_ = UNKNOWN_

        !provavelmente a propriedade não se encontra em PropertyValue - corrigir isto
        x               => Me%ExternalVar%TXYZP_Points%X
        y               => Me%ExternalVar%TXYZP_Points%Y
        z               => Me%ExternalVar%TXYZP_Points%Z
        !CUIDADO!!!!
        w               => Me%ExternalVar%TXYZP_Points%Z
        n_obs           =  Me%ExternalVar%TXYZP_Points%Count
        n_dimensions    =  Me%ComputeOptions%n_dimensions

        select case (n_dimensions)

        case (1)
            allocate(yy(n_obs))
            allocate(originY(n_obs))
            allocate(destinationY(n_obs))

            call sorter(y,yy,n_obs,originY,destinationY)

            dim_Y=n_obs
            dim_X=1
            dim_Z=1

            allocate (CellHasPoint(dim_Y,1,1),CellValue(dim_Y,1,1))
            allocate (Cell_YBounds(dim_Y+1))
            allocate (merge_Y(dim_Y))
        
            Cell_YBounds(1)=yy(1)-0.5d0*(yy(2)-yy(1))
            Cell_YBounds(n_obs+1)=yy(n_obs)-0.5d0*(yy(n_obs)-yy(n_obs-1))

            do a=1,n_obs
                CellHasPoint(a,1,1)=.true.
                CellValue(a,1,1)=w(originY(a))
            enddo

            do a=2,n_obs
                Cell_YBounds(a)=0.5d0*(yy(a-1)+yy(a))
            enddo

            ok_Y_merge=.true.
            warning_Y=.false.
            ok_X_merge=.false.
            ok_Z_merge=.false.

        case (2)
            allocate(xx(n_obs),yy(n_obs))
            allocate(originX(n_obs),originY(n_obs))
            allocate(destinationX(n_obs),destinationY(n_obs))
        
            call sorter(x,xx,n_obs,originX,destinationX)
            call sorter(y,yy,n_obs,originY,destinationY)

            dim_X=n_obs
            dim_Y=n_obs
            dim_Z=1

            allocate (CellHasPoint(dim_Y,dim_X,1),CellValue(dim_Y,dim_X,1))
            allocate (Cell_XBounds(dim_X+1),Cell_YBounds(dim_Y+1))
            allocate (merge_X(dim_X), merge_Y(dim_Y))
        
            Cell_XBounds(1)=xx(1)-0.5d0*(xx(2)-xx(1))
            Cell_YBounds(1)=yy(1)-0.5d0*(yy(2)-yy(1))

            Cell_XBounds(n_obs+1)=xx(n_obs)-0.5d0*(xx(n_obs)-xx(n_obs-1))
            Cell_YBounds(n_obs+1)=yy(n_obs)-0.5d0*(yy(n_obs)-yy(n_obs-1))

            CellHasPoint(:,:,:)=.false.
            CellValue(:,:,:)=0d0
            do a=1,n_obs
                CellHasPoint(a,destinationX(originY(a)),1)=.true.
                CellValue(a,destinationX(originY(a)),1)=w(originY(a))
            enddo

            do a=2,n_obs
                Cell_XBounds(a)=0.5d0*(xx(a-1)+xx(a))
                Cell_YBounds(a)=0.5d0*(yy(a-1)+yy(a))
            enddo

            ok_Y_merge=.true.
            warning_Y=.false.
            ok_X_merge=.true.
            warning_X=.false.
            ok_Z_merge=.false.

        case (3)
            allocate(xx(n_obs),yy(n_obs),zz(n_obs))
            allocate(originX(n_obs),originY(n_obs),originZ(n_obs))
            allocate(destinationX(n_obs),destinationY(n_obs),destinationZ(n_obs))
        
            call sorter(x,xx,n_obs,originX,destinationX)
            call sorter(y,yy,n_obs,originY,destinationY)
            call sorter(z,zz,n_obs,originZ,destinationZ)

            dim_X=n_obs
            dim_Y=n_obs
            dim_Z=n_obs

            allocate (CellHasPoint(dim_Y,dim_X,dim_Z))
            allocate (Cell_XBounds(dim_X+1),Cell_YBounds(dim_Y+1),Cell_ZBounds(dim_Z+1))
            allocate (merge_X(dim_X), merge_Y(dim_Y), merge_Z(dim_Z))
        
            Cell_XBounds(1)=xx(1)-0.5d0*(xx(2)-xx(1))
            Cell_YBounds(1)=yy(1)-0.5d0*(yy(2)-yy(1))
            Cell_ZBounds(1)=zz(1)-0.5d0*(zz(2)-zz(1))

            Cell_XBounds(n_obs+1)=xx(n_obs)-0.5d0*(xx(n_obs)-xx(n_obs-1))
            Cell_YBounds(n_obs+1)=yy(n_obs)-0.5d0*(yy(n_obs)-yy(n_obs-1))
            Cell_ZBounds(n_obs+1)=zz(n_obs)-0.5d0*(zz(n_obs)-zz(n_obs-1))

            CellHasPoint(:,:,:)=.false.
            CellValue(:,:,:)=0d0
            do a=1,n_obs
                CellHasPoint(a,destinationX(originY(a)),destinationZ(originY(a)))=.true.
                CellValue(a,destinationX(originY(a)),destinationZ(originY(a)))=w(originY(a))
            enddo

            do a=2,n_obs
                Cell_XBounds(a)=0.5d0*(xx(a-1)+xx(a))
                Cell_YBounds(a)=0.5d0*(yy(a-1)+yy(a))
                Cell_ZBounds(a)=0.5d0*(zz(a-1)+zz(a))
            enddo
        
            ok_Y_merge=.true.
            warning_Y=.false.
            ok_X_merge=.true.
            warning_X=.false.
            ok_Z_merge=.true.
            warning_Z=.false.

        end select
        
        not_finished=.true.

        do while (not_finished)
            do a=0,1
                
                if(ok_Y_merge) then
                    merge_Y(:)=.false.
                    i=1+a
                    new_dimY=dim_Y
                    do while(i.lt.dim_Y)
                        merge_Y(i)=.true.
                        j=1
                        do while(j.le.dim_X.and.merge_Y(i))
                            k=1
                            do while(k.le.dim_Z.and.merge_Y(i))
                                if(CellHasPoint(i,j,k).and.CellHasPoint(i+1,j,k)) merge_Y(i)=.false.
                                k=k+1
                            enddo
                            j=j+1
                        enddo
                        if(merge_Y(i)) new_dimY=new_dimY-1
                        i=i+2
                    enddo
                    if(new_dimY.lt.dim_Y) then
                        warning_Y=.false.
                        allocate (new_CHP(new_dimY,dim_X,dim_Z),new_CV(new_dimY,dim_X,dim_Z),new_YBounds(new_dimY+1))
                        work_Y=1
                        i=1
                        do while(i.le.dim_Y)
                            if(merge_Y(i)) then
                                do j=1,dim_X
                                do k=1,dim_Z
                                    new_CHP(work_Y,j,k)=CellHasPoint(i,j,k).or.CellHasPoint(i+1,j,k)
                                    new_CV(work_Y,j,k)=CellValue(i,j,k)+CellValue(i+1,j,k)
                                enddo
                                enddo
                                new_YBounds(work_Y)=Cell_YBounds(i)
                                new_YBounds(work_Y+1)=Cell_YBounds(i+2)
                                i=i+2
                                else
                                do j=1,dim_X
                                do k=1,dim_Z
                                    new_CHP(work_Y,j,k)=CellHasPoint(i,j,k)
                                    new_CV(work_Y,j,k)=CellValue(i,j,k)
                                enddo
                                enddo
                                new_YBounds(work_Y)=Cell_YBounds(i)
                                new_YBounds(work_Y+1)=Cell_YBounds(i+1)
                                i=i+1
                            endif
                            work_Y=work_Y+1
                        enddo
                        deallocate (CellHasPoint,CellValue,Cell_YBounds,merge_Y)
                        dim_Y=new_dimY
                        allocate (CellHasPoint(dim_Y,dim_X,dim_Z),CellValue(dim_Y,dim_X,dim_Z),Cell_YBounds(dim_Y+1),merge_Y(dim_Y))
                        CellHasPoint(:,:,:)=new_CHP(:,:,:)
                        CellValue(:,:,:)=new_CV(:,:,:)
                        Cell_YBounds(:)=new_YBounds(:)
                        deallocate(new_CHP,new_CV,new_YBounds)
                        else
                        if(warning_Y) then 
                            ok_Y_merge=.false.
                            else
                            warning_Y=.true.
                        endif
                    endif
                endif

                if(ok_X_merge) then
                    merge_X(:)=.false.
                    j=1+a
                    new_dimX=dim_X
                    do while(j.lt.dim_X)
                        merge_X(j)=.true.
                        i=1
                        do while(i.le.dim_Y.and.merge_X(j))
                            k=1
                            do while(k.le.dim_Z.and.merge_X(j))
                                if(CellHasPoint(i,j,k).and.CellHasPoint(i,j+1,k)) merge_X(j)=.false.
                                k=k+1
                            enddo
                            i=i+1
                        enddo
                        if(merge_X(j)) new_dimX=new_dimX-1
                        j=j+2
                    enddo
                    if(new_dimX.lt.dim_X) then
                        warning_X=.false.
                        allocate (new_CHP(dim_Y,new_dimX,dim_Z),new_CV(dim_Y,new_dimX,dim_Z),new_XBounds(new_dimX+1))
                        work_X=1
                        j=1
                        do while(j.le.dim_X)
                            if(merge_X(j)) then
                                do i=1,dim_Y
                                do k=1,dim_Z
                                    new_CHP(i,work_X,k)=CellHasPoint(i,j,k).or.CellHasPoint(i,j+1,k)
                                    new_CV(i,work_X,k)=CellValue(i,j,k)+CellValue(i,j+1,k)
                                enddo
                                enddo
                                new_XBounds(work_X)=Cell_XBounds(j)
                                new_XBounds(work_X+1)=Cell_XBounds(j+2)
                                j=j+2
                                else
                                do i=1,dim_Y
                                do k=1,dim_Z
                                    new_CHP(i,work_X,k)=CellHasPoint(i,j,k)
                                    new_CV(i,work_X,k)=CellValue(i,j,k)
                                enddo
                                enddo
                                new_XBounds(work_X)=Cell_XBounds(j)
                                new_XBounds(work_X+1)=Cell_XBounds(j+1)
                                j=j+1
                            endif
                            work_X=work_X+1
                        enddo
                        deallocate (CellHasPoint,CellValue,Cell_XBounds,merge_X)
                        dim_X=new_dimX
                        allocate (CellHasPoint(dim_Y,dim_X,dim_Z),CellValue(dim_Y,dim_X,dim_Z),Cell_XBounds(dim_X+1),merge_X(dim_X))
                        CellHasPoint(:,:,:)=new_CHP(:,:,:)
                        CellValue(:,:,:)=new_CV(:,:,:)
                        Cell_XBounds(:)=new_XBounds(:)
                        deallocate(new_CHP,new_CV,new_XBounds)
                        else
                        if(warning_X) then 
                            ok_X_merge=.false.
                            else
                            warning_X=.true.
                        endif
                    endif
                endif

                if(ok_Z_merge) then
                    merge_Z(:)=.false.
                    k=1+a
                    new_dimZ=dim_Z
                    do while(k.lt.dim_Z)
                        merge_Z(k)=.true.
                        i=1
                        do while(i.le.dim_Y.and.merge_Z(k))
                            j=1
                            do while(j.le.dim_X.and.merge_Z(k))
                                if(CellHasPoint(i,j,k).and.CellHasPoint(i,j+1,k)) merge_Z(k)=.false.
                                j=j+1
                            enddo
                            i=i+1
                        enddo
                        if(merge_Z(k)) new_dimZ=new_dimZ-1
                        k=k+2
                    enddo
                    if(new_dimZ.lt.dim_Z) then
                        warning_Z=.false.
                        allocate (new_CHP(dim_Y,dim_X,new_dimZ),new_CV(dim_Y,dim_X,new_dimZ),new_ZBounds(new_dimZ+1))
                        work_Z=1
                        k=1
                        do while(k.le.dim_Z)
                            if(merge_Z(k)) then
                                do i=1,dim_Y
                                do j=1,dim_X
                                    new_CHP(i,j,work_Z)=CellHasPoint(i,j,k).or.CellHasPoint(i,j,k+1)
                                    new_CV(i,j,work_Z)=CellValue(i,j,k)+CellValue(i,j,k+1)
                                enddo
                                enddo
                                new_ZBounds(work_Z)=Cell_ZBounds(k)
                                new_ZBounds(work_Z+1)=Cell_ZBounds(k+2)
                                k=k+2
                                else
                                do i=1,dim_Y
                                do j=1,dim_X
                                    new_CHP(i,j,work_Z)=CellHasPoint(i,j,k)
                                    new_CV(i,j,work_Z)=CellValue(i,j,k)
                                enddo
                                enddo
                                new_ZBounds(work_Z)=Cell_ZBounds(k)
                                new_ZBounds(work_Z+1)=Cell_ZBounds(k+1)
                                k=k+1
                            endif
                            work_Z=work_Z+1
                        enddo
                        deallocate (CellHasPoint,CellValue,Cell_ZBounds,merge_Z)
                        dim_Z=new_dimZ
                        allocate (CellHasPoint(dim_Y,dim_X,dim_Z),CellValue(dim_Y,dim_X,dim_Z),Cell_ZBounds(dim_Z+1),merge_Z(dim_Z))
                        CellHasPoint(:,:,:)=new_CHP(:,:,:)
                        CellValue(:,:,:)=new_CV(:,:,:)
                        Cell_ZBounds(:)=new_ZBounds(:)
                        deallocate(new_CHP,new_CV,new_ZBounds)
                        else
                        if(warning_Z) then 
                            ok_Z_merge=.false.
                            else
                            warning_Z=.true.
                        endif
                    endif
                endif

            enddo

            not_finished=ok_X_merge.or.ok_Y_merge.or.ok_Z_merge
        enddo

        !Call ConstructLocalFather
        call ConstructLocalFather(1,1,1, dim_X, dim_Y, dim_Z, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructPoints2Grid - ModuleInterpolation - ERR01'
            
        select case (n_dimensions)
        case (1)
            do i=1,dim_Y
                if(CellHasPoint(i,j,1)) then
                    Me%LocalFatherGrid%OpenPoints3D(i,1,1) = 1
                    Me%LocalFatherGrid%CenterY(i,1) = 0.5d0*(Cell_YBounds(i)+Cell_YBounds(i+1))
                    Me%LocalFatherField%Values3D(i,1,1) = CellValue(i,1,1)
                endif
            enddo
            deallocate(yy,originY,destinationY,Cell_YBounds,merge_Y,CellHasPoint,CellValue)
        case (2)
            do i=1,dim_Y
            do j=1,dim_X
                if(CellHasPoint(i,j,1)) then
                    Me%LocalFatherGrid%OpenPoints3D(i,j,1) = 1
                    Me%LocalFatherGrid%CenterY(i,j) = 0.5d0*(Cell_YBounds(i)+Cell_YBounds(i+1))
                    Me%LocalFatherGrid%CenterX(i,j) = 0.5d0*(Cell_XBounds(i)+Cell_XBounds(i+1))
                    Me%LocalFatherField%Values3D(i,j,1) = CellValue(i,j,1)
                endif
            enddo
            enddo
            deallocate(yy,originY,destinationY,Cell_YBounds,merge_Y,CellHasPoint,CellValue)
            deallocate(xx,originX,destinationX,Cell_XBounds,merge_X)
        case (3)
            do i=1,dim_Y
            do j=1,dim_X
            do k=1,dim_Z
                if(CellHasPoint(i,j,k)) then
                    Me%LocalFatherGrid%OpenPoints3D(i,j,k) = 1
                    Me%LocalFatherGrid%CenterY(i,j) = 0.5d0*(Cell_YBounds(i)+Cell_YBounds(i+1))
                    Me%LocalFatherGrid%CenterX(i,j) = 0.5d0*(Cell_XBounds(i)+Cell_XBounds(i+1))
                    Me%LocalFatherGrid%CenterZ(i,j,k) = 0.5d0*(Cell_ZBounds(i)+Cell_ZBounds(i+1))
                    Me%LocalFatherField%Values3D(i,j,k) = CellValue(i,j,k)
                endif
            enddo
            enddo
            enddo
            deallocate(yy,originY,destinationY,Cell_YBounds,merge_Y,CellHasPoint,CellValue)
            deallocate(xx,originX,destinationX,Cell_XBounds,merge_X)
            deallocate(zz,originZ,destinationZ,Cell_ZBounds,merge_Z)
        end select

        call NC_Data_3D (.true.)

        nullify(x,y,z,w)

        STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

    end subroutine ConstructPoints2Grid_3D

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    subroutine ConstructLocalFather(ILB,IUB,JLB,JUB,KLB,KUB,STAT)

        !Arugments---------------------------------------------------------------
        integer                                 :: ILB,IUB,JLB,JUB,KLB,KUB
        integer, optional, intent(OUT)          :: STAT

        !Locals------------------------------------------------------------------
        integer                                 :: STAT_, STAT_CALL
        
        STAT_ = UNKNOWN_
        
        call AllocateGrid(Me%LocalFatherGrid)
        call AllocateField(Me%LocalFatherField)
        
        Me%LocalFatherGrid%WorkSize3D%ILB = ILB
        Me%LocalFatherGrid%WorkSize3D%JLB = JLB
        Me%LocalFatherGrid%WorkSize3D%KLB = KLB
        Me%LocalFatherGrid%WorkSize3D%IUB = IUB
        Me%LocalFatherGrid%WorkSize3D%JUB = JUB
        Me%LocalFatherGrid%WorkSize3D%KUB = KUB

        call AllocateLocalFatherGrid(Me%LocalFatherGrid, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructLocalFather - ModuleInterpolation - ERR01'

        call AllocateLocalFatherField(Me%LocalFatherField, Me%LocalFatherGrid%WorkSize3D, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructLocalFather - ModuleInterpolation - ERR01'

        Me%LocalFatherGrid%OpenPoints3D(:,:,:) = 0
        Me%LocalFatherGrid%CenterX(:,:) = null_real
        Me%LocalFatherGrid%CenterY(:,:) = null_real
        Me%LocalFatherGrid%CenterZ(:,:,:) = null_real
        Me%ExternalVar%FatherField%Values3D(:,:,:) = null_real

        STAT_ = SUCCESS_

        if(present(STAT)) STAT = STAT_

    end subroutine ConstructLocalFather
    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    
    subroutine AllocateLocalFatherGrid(Gridpointer, STAT)

        !Arguments--------------------------------------------------------------
        Type(T_Grid), pointer                           :: Gridpointer
        integer, optional, intent(OUT)                  :: STAT

        !Locals-----------------------------------------------------------------
        integer                                         :: STAT_
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB

        !Begin shorten varibles names-------------------------------------------
        ILB = Gridpointer%WorkSize3D%ILB
        IUB = Gridpointer%WorkSize3D%IUB
        JLB = Gridpointer%WorkSize3D%JLB
        JUB = Gridpointer%WorkSize3D%JUB
        KLB = Gridpointer%WorkSize3D%KLB
        KUB = Gridpointer%WorkSize3D%KUB

        STAT_ = UNKNOWN_

        allocate(Gridpointer%OpenPoints3D(ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Gridpointer%WaterPoints3D(ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Gridpointer%WaterPoints2D(ILB:IUB, JLB:JUB))
        allocate(Gridpointer%CenterX(ILB:IUB, JLB:JUB))
        allocate(Gridpointer%CenterY(ILB:IUB, JLB:JUB))
        allocate(Gridpointer%CenterZ(ILB:IUB, JLB:JUB, KLB:KUB))

        STAT_ = SUCCESS_

        if(present(STAT)) STAT = STAT_

    end subroutine AllocateLocalFatherGrid

    !----------------------------------------------------------------------------
    
    !----------------------------------------------------------------------------
    
    subroutine AllocateLocalFatherField(Fieldpointer, WorkSize3D, STAT)

        !Arguments--------------------------------------------------------------
        Type(T_Field), pointer                          :: Fieldpointer
        Type(T_Size3D), pointer                         :: WorkSize3D
        integer, optional, intent(OUT)                  :: STAT

        !Locals-----------------------------------------------------------------
        integer                                         :: STAT_
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB

        !Begin shorten varibles names-------------------------------------------
        ILB = WorkSize3D%ILB
        IUB = WorkSize3D%IUB
        JLB = WorkSize3D%JLB
        JUB = WorkSize3D%JUB
        KLB = WorkSize3D%KLB
        KUB = WorkSize3D%KUB

        STAT_ = UNKNOWN_

        allocate(Fieldpointer%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))
        if (associated(Fieldpointer%Values3D)) STAT_ = SUCCESS_

        if(present(STAT)) STAT = STAT_
 
    end subroutine AllocateLocalFatherField

    !----------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Local-----------------------------------------------------------------
        type (T_Interpolation), pointer           :: NewInterpolation
        type (T_Interpolation), pointer           :: PreviousInterpolation

        !Allocates new instance
        allocate (NewInterpolation)
        nullify  (NewInterpolation%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjInterpolation)) then
            FirstObjInterpolation      => NewInterpolation
            Me                           => NewInterpolation
        else
            PreviousInterpolation     => FirstObjInterpolation
            Me                           => FirstObjInterpolation%Next
            do while (associated(Me))
                PreviousInterpolation  => Me
                Me                       => Me%Next
            enddo
            Me                           => NewInterpolation
            PreviousInterpolation%Next => NewInterpolation
        endif

        Me%InstanceID = RegisterNewInstance (mINTERPOLATION_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine AllocateExternal(STAT)

        !Arguments-------------------------------------------------------------
        integer, optional, intent(OUT)          :: STAT

        !Local-----------------------------------------------------------------
        type (T_External), pointer              :: NewExternal
        integer                                 :: STAT_
 
        STAT_ = UNKNOWN_

        !Allocates new instance
        allocate (NewExternal)
 
        Me%ExternalVar => NewExternal

        if (associated(Me%ExternalVar)) STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

    end subroutine AllocateExternal

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine AllocateGrid(GridPointer, STAT)

        !Argument--------------------------------------------------------------
        type (T_Grid), pointer           :: GridPointer
        integer, optional, intent(OUT)   :: STAT

        !Local-----------------------------------------------------------------
        type (T_Grid), pointer           :: NewGrid
        integer                          :: STAT_

        STAT_ = UNKNOWN_

        !Allocates new instance
        allocate (NewGrid)
        allocate (NewGrid%Size3D)
        allocate (NewGrid%WorkSize3D)
        allocate (NewGrid%Size2D)
        allocate (NewGrid%WorkSize2D)

        GridPointer => NewGrid

        if(associated(GridPointer)) STAT_ = SUCCESS_

        if(present(STAT)) STAT = STAT_

    end subroutine AllocateGrid

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine AllocateField(FieldPointer, STAT)

        !Argument--------------------------------------------------------------
        type (T_Field), pointer           :: FieldPointer
        integer, optional, intent(OUT)   :: STAT

        !Local-----------------------------------------------------------------
        type (T_Field), pointer          :: NewField
        integer                          :: STAT_

        STAT_ = UNKNOWN_

        !Allocates new instance
        allocate (NewField)
        nullify(NewField%Values2D)
        nullify(NewField%Values3D)

        FieldPointer => NewField

        if(present(STAT)) STAT = STAT_

    end subroutine AllocateField

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine PopulateParameters(GridPointer, STAT)

        !Arguments------------------------------------------------
        type (T_Grid), pointer              :: GridPointer
        integer, optional, intent(OUT)      :: STAT

        !Local----------------------------------------------------------------------
        integer                             :: STAT_
        
        STAT_ = UNKNOWN_

        !Nullify fields
        nullify(Me%MatrixUpdate3D)

        select case(Me%ComputeOptions%TypeZUV)

            case(TypeZ_)                

            case(TypeU_)
              GridPointer%Size3D%JUB        = GridPointer%Size3D%JUB + 1
              GridPointer%WorkSize3D%JUB        = GridPointer%WorkSize3D%JUB + 1
                
            case(TypeV_)
              GridPointer%Size3D%IUB        = GridPointer%Size3D%IUB + 1
              GridPointer%WorkSize3D%IUB        = GridPointer%WorkSize3D%IUB + 1
           
            case default
                
                write(*,*)'Invalid type ZUV in grid data '
                stop 'PopulateParameters - ModuleInterpolation - ERR90'

        end select

        if(present(STAT)) STAT = STAT_

    end subroutine PopulateParameters
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadLockExternalVar(GridPointer, STAT)
        
        !Argument--------------------------------------------------------------
        type (T_Grid), pointer                  :: GridPointer
        integer, optional, intent(OUT)          :: STAT

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL, STAT_
                
        STAT_ = UNKNOWN_

            !WaterPoints2D
            call GetWaterPoints2D(GridPointer%ObjHorizontalMap, GridPointer%WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR10'

            !OpenPoints3D
            call GetOpenPoints3D(GridPointer%ObjMap, GridPointer%OpenPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR20'

            STAT_ = SUCCESS_
      
            !Size2D and WorkSize2D
            call GetHorizontalGridSize(GridPointer%ObjHorizontalGrid,                            &
                                       WorkSize = GridPointer%WorkSize2D,                        &
                                       Size     = GridPointer%Size2D,                            &
                                       STAT     = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ReadLockExternalVar -  ModuleInterpolateGrids - ERR30'

            STAT_ = SUCCESS_

        !FatherGrid
        select case (Me%ComputeOptions%TypeZUV)

            case(TypeZ_)
                call ReadLockExternalZ(GridPointer, STAT_)

            case(TypeU_)
                call ReadLockExternalU(GridPointer, STAT_)

            case(TypeV_)
                call ReadLockExternalV(GridPointer, STAT_)
            
            case default
                
                write(*,*)'Invalid type ZUV in grid data '
                stop 'ReadLockExternalVar - ModuleInterpolation - ERR40'

        end select

        if (present(STAT)) STAT = STAT_

    end subroutine ReadLockExternalVar

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    Subroutine ReadLockExternalZ(GridPointer, STAT)
    
        !Arguments------------------------------------------------
        type (T_Grid), pointer              :: GridPointer
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: STAT_

        STAT_ = UNKNOWN_

            !XX_IE and YY_IE
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid,                                   &
                                    XX_IE = GridPointer%XX_IE,                           &
                                    YY_IE = GridPointer%YY_IE,                           &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR04'

            STAT_ = SUCCESS_

            !DX
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid, DUX = GridPointer%DX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR08'
        
            !DY
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid, DVY = GridPointer%DY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR09'

            !Size3D and WorkSize3D
            call GetGeometrySize(GridPointer%ObjGeometry,                     &
                                 Size    = GridPointer%Size3D,                &
                                 WorkSize    = GridPointer%WorkSize3D,        &
                                 STAT    = STAT_CALL)             
            if (STAT_CALL /= SUCCESS_)                              &
            stop 'ReadLockExternalVar - ModuleInterpolation - ERR01'

            !DZ
            call GetGeometryDistances (GridPointer%ObjGeometry, DWZ = GridPointer%DZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR07b'

            !CX
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid, DZX = GridPointer%CX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR09'

            !CY
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid, DZY = GridPointer%CY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR09'

            !CZ
            call GetGeometryDistances (GridPointer%ObjGeometry, DZZ = GridPointer%CZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR09'

            !SZZ
            call GetGeometryDistances (GridPointer%ObjGeometry, SZZ = GridPointer%SZZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR10'

            STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

    end subroutine ReadLockExternalZ

    !--------------------------------------------------------------------------

    Subroutine ReadLockExternalU(GridPointer, STAT)
    
        !Arguments------------------------------------------------
        type (T_Grid), pointer              :: GridPointer
        integer, optional, intent(OUT)      :: STAT

        !Local--------------------------------------------------------------
        integer                                 :: STAT_CALL, STAT_
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                             :: i,j

        !These variables point to object Geometry. They are used to define DX, XX_IE and YY_IE.
        !Formulas for U cell are:
        !       DX = DZX(j-1)
        !       XX_IE = 0.5 * [ XX_IE + XX_IE(j-1) ]
        !       YY_IE = 0.5 * [ YY_IE + YY_IE(j-1) ]
        !See "MOHID Interpolation Manual" for a figure and more details
        real, dimension(:,:), pointer           :: XX_IE, YY_IE, DZX

        STAT_ = UNKNOWN_

        !Begin shorten variables names---------------------------------------
        ILB   = GridPointer%WorkSize3D%ILB
        IUB   = GridPointer%WorkSize3D%IUB

        JLB   = GridPointer%WorkSize3D%JLB 
        JUB   = GridPointer%WorkSize3D%JUB

        KLB   = GridPointer%WorkSize3D%KLB
        KUB   = GridPointer%WorkSize3D%KUB

        !XX_IE and YY_IE--------------------------------------------------------------------

            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid,                                   &
                                XX_IE = XX_IE,                           &
                                    YY_IE = YY_IE,                           &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR04'
        
            STAT_ = SUCCESS_


        !XX_IE
        do j = JLB, JUB
        do i = ILB, IUB
            if (j == JLB) then
                GridPointer%XX_IE(i,JLB) = 1.5 * XX_IE(i,JLB) - 0.5 * XX_IE(i,JLB + 1)
            else
                GridPointer%XX_IE(i,j) = 0.5 * ( XX_IE(i,j) + XX_IE(i,j-1) )
            endif
        enddo
        enddo

        !YY_IE
        do j = JLB, JUB
        do i = ILB, IUB
            if (j == JLB) then
                GridPointer%YY_IE(i,JLB) = 1.5 * YY_IE(i,JLB) - 0.5 * YY_IE(i,JLB + 1)
            else
                GridPointer%YY_IE(i,j) = 0.5 * ( YY_IE(i,j) + YY_IE(i,j-1) )
            endif
        enddo
        enddo

        !DX-------------------------------------------------------------------------------
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid, DZX = DZX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR08'

            do j = JLB, JUB
            do i = ILB, IUB
                if(j == JLB) then
                    GridPointer%DX(i,JLB) = DZX(i,JLB)
                else
                    GridPointer%DX(i,j) = DZX(i,j-1)
                endif
            enddo
            enddo
        
            !DY-------------------------------------------------------------------------------
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid, DYY = GridPointer%DY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR09'

            !DZ--------------------------------------------------------------------------------
            call GetGeometryDistances (GridPointer%ObjGeometry, DUZ = GridPointer%DZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR07b'

            !CX--------------------------------------------------------------------------------
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid, DUX = GridPointer%CX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR09'

            !CY--------------------------------------------------------------------------------
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid, DUY = GridPointer%CY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR09'

            !CZ--------------------------------------------------------------------------------
            call GetGeometryDistances (GridPointer%ObjGeometry, DZE = GridPointer%CZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR09'

            !SZZ-------------------------------------------------------------------------------
            call GetGeometryDistances (GridPointer%ObjGeometry, SZZ = GridPointer%SZZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR10'

            STAT_ = SUCCESS_

        if(present(STAT)) STAT = STAT_

    end subroutine ReadLockExternalU

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    Subroutine ReadLockExternalV(GridPointer, STAT)
    
        !Arguments------------------------------------------------
        type (T_Grid), pointer              :: GridPointer
        integer, optional, intent(OUT)      :: STAT

        !Local--------------------------------------------------------------
        integer                             :: STAT_CALL
        integer                             :: STAT_
        integer                             :: i,j
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB

        !These variables point to object Geometry. They are used to define DX, XX_IE and YY_IE.
        !Formulas for U cell are:
        !       DY = DXX(j-1)
        !       XX_IE = 0.5 * [ XX_IE + XX_IE(i-1) ]
        !       YY_IE = 0.5 * [ YY_IE + YY_IE(i-1) ]
        !See "MOHID Interpolation Manual" for a figure and more details
        real, dimension(:,:), pointer           :: XX_IE, YY_IE, DZY

        STAT_ = UNKNOWN_

        !Begin shorten variables names---------------------------------------
        ILB   = GridPointer%WorkSize3D%ILB
        IUB   = GridPointer%WorkSize3D%IUB

        JLB   = GridPointer%WorkSize3D%JLB 
        JUB   = GridPointer%WorkSize3D%JUB

        KLB   = GridPointer%WorkSize3D%KLB
        KUB   = GridPointer%WorkSize3D%KUB

        !XX_IE and YY_IE--------------------------------------------------------------------
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid,                                   &
                                    XX_IE = XX_IE,                           &
                                    YY_IE = YY_IE,                           &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR04'

            STAT_ = SUCCESS_
        !XX_IE
        do j = JLB, JUB
        do i = ILB, IUB
            if (i == ILB) then
                GridPointer%XX_IE(ILB,j) = 1.5 * XX_IE(ILB,j) - 0.5 * XX_IE(ILB + 1,j)
            else
                GridPointer%XX_IE(i,j) = 0.5 * ( XX_IE(i,j) + XX_IE(i-1,j) )
            endif
        enddo
        enddo

        !YY_IE
        do j = JLB, JUB
        do i = ILB, IUB
            if (i == ILB) then
                GridPointer%YY_IE(ILB,j) = 1.5 * YY_IE(ILB,j) - 0.5 * YY_IE(ILB + 1,j)
            else
                GridPointer%YY_IE(i,j) = 0.5 * ( YY_IE(i,j) + YY_IE(i-1,j) )
            endif
        enddo
        enddo

            !DX-------------------------------------------------------------------------------
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid, DXX = GridPointer%DX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR08'
        
            !DY-------------------------------------------------------------------------------
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid, DZY = DZY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR09'

            do j = JLB, JUB
            do i = ILB, IUB
                if(i == ILB) then
                    GridPointer%DY(ILB,j) = DZY(ILB,j)
                else
                    GridPointer%DY(i,j) = DZY(i-1,j)
                endif
            enddo
            enddo

            !DZ--------------------------------------------------------------------------------
            call GetGeometryDistances (GridPointer%ObjGeometry, DVZ = GridPointer%DZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR07b'

            !CX--------------------------------------------------------------------------------
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid, DVX = GridPointer%CX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR09'

            !CY--------------------------------------------------------------------------------
            call GetHorizontalGrid (GridPointer%ObjHorizontalGrid, DVY = GridPointer%CY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR09'

            !CZ--------------------------------------------------------------------------------
            call GetGeometryDistances (GridPointer%ObjGeometry, DZI = GridPointer%CZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR09'

            !SZZ-------------------------------------------------------------------------------
            call GetGeometryDistances (GridPointer%ObjGeometry, SZZ = GridPointer%SZZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterpolation - ERR09'
            
            STAT_ = SUCCESS_

        if(present(STAT)) STAT = STAT_

    end subroutine ReadLockExternalV

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    !--------------------------------------------------------------------------
 
    subroutine GetInterpolationPointer (ObjInterpolationID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjInterpolationID
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjInterpolationID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mInterpolation_, Me%InstanceID)

!            Matrix => Me%Matrix

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetInterpolationPointer
    
    !--------------------------------------------------------------------------
    
    subroutine GetInterpolationInteger (ObjInterpolationID, Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjInterpolationID
        real                                            :: Int
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjInterpolationID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Int = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetInterpolationInteger

    !--------------------------------------------------------------------------

    subroutine UnGetInterpolation3D_I(ObjInterpolationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjInterpolationID
        integer, dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjInterpolationID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mInterpolation_, Me%InstanceID, "UnGetInterpolation3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetInterpolation3D_I

    !--------------------------------------------------------------------------

    subroutine UnGetInterpolation3D_R8(ObjInterpolationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjInterpolationID
        real(8), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjInterpolationID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mInterpolation_, Me%InstanceID,  "UnGetInterpolation3D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetInterpolation3D_R8

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !----------------------------------------------------------------------------

    subroutine ModifyInterpolator3D(ObjInterpolationID, SonField3D, FatherField3D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjInterpolationID
        real,    dimension(:,:,:), pointer          :: SonField3D
        real,    dimension(:,:,:), pointer, optional:: FatherField3D
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        logical                                     :: ErrorON

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjInterpolationID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            ErrorON = .false.

            if (present(FatherField3D)) Me%ExternalVar%FatherField%Values3D => FatherField3D
            Me%ExternalVar%SonField%Values3D => SonField3D

c1:         Select Case (Me%ComputeOptions%Methodology)

                case (ConvolutionConserv_)

                    !call ConservativeModifier   (Me%ExternalVar%FatherField%Values3D) !here we go to either 2D or 3D subroutines

                case (ConvolutionNonConserv_)

                    !call NonConservativeModifier(Me%ExternalVar%FatherField%Values3D) !here we go to either 2D or 3D subroutines

                case default

                    if (Me%StructuredData) then
                        call ClassicInterpolaStruct3D
                    else
                        !call ClassicInterpolaUnStruct3D
                    endif

            end select c1

            STAT_ = SUCCESS_

        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        !call interface


    end subroutine ModifyInterpolator3D


    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine ModifyInterpolator2D(ObjInterpolationID, GridData2D, Data1D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjInterpolationID
        real,    dimension(:,:  ), pointer, optional:: GridData2D
        real,    dimension(:    ), pointer, optional:: Data1D
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        logical                                     :: ErrorON

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjInterpolationID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            ErrorON = .false.

            if (present(GridData2D)) then
                Me%ExternalVar%FatherField%Values2D => GridData2D
            endif
            if (present(Data1D))     then
                !CUIDADO!!!!
                Me%ExternalVar%TXYZP_Points%Z       => Data1D
            endif
            
c1:         Select Case (Me%ComputeOptions%Methodology)

                case (ConvolutionConserv_)

                    !call ConservativeModifier   (Me%ExternalVar%FatherField%Values2D) !here we go to either 2D or 3D subroutines

                case (ConvolutionNonConserv_)

                    !call NonConservativeModifier(Me%ExternalVar%FatherField%Values2D) !here we go to either 2D or 3D subroutines

                case default

                    if (Me%StructuredData) then
                        call ClassicInterpolaStruct2D
                    else
                        call ClassicInterpolaUnStruct2D
                    endif

            end select c1

        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        !call interface


    end subroutine ModifyInterpolator2D


    !----------------------------------------------------------------------------



    !----------------------------------------------------------------------------

    subroutine ClassicInterpolaUnStruct2D

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        !From FatherGrid grid to SonGrid grid
        
        if(Me%ComputeOptions%Methodology == Triangulation_)then

            call UnStructureTriangulation (Me%ExternalVar%TXYZP_Points, Me%ExternalVar%SonField, Me%ExternalVar%SonGrid)

        elseif(Me%ComputeOptions%Methodology == InverseWeight_)then

            call UnStructureInverseWeight (Me%ExternalVar%TXYZP_Points, Me%ExternalVar%SonField, Me%ExternalVar%SonGrid)

        else

            write(*,*) 'Unknown type of interpolation'
            stop       'StartInterpolateGrids - ModuleInterpolateGrids - ERR250' 

        end if

    end subroutine ClassicInterpolaUnStruct2D

    !--------------------------------------------------------------------------

    subroutine ClassicInterpolaStruct3D

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,    dimension(:,:), pointer            :: InValues2D, OutValues2D
        integer                                     :: k


        !----------------------------------------------------------------------

        !From FatherGrid grid to local FatherGrid grid
        
        do k=Me%ExternalVar%FatherGrid%WorkSize3D%KLB,Me%ExternalVar%FatherGrid%WorkSize3D%KUB

            Me%LocalFatherField%Values2D            (:,:) = FillValueReal

            Me%ExternalVar%FatherField%Values2D     (:,:) = Me%ExternalVar%FatherField%Values3D (:,:,k)
            Me%ExternalVar%FatherGrid%WaterPoints2D     (:,:) = Me%ExternalVar%FatherGrid%WaterPoints3D (:,:,k)
            Me%LocalFatherGrid%WaterPoints2D            (:,:) = Me%LocalFatherGrid%WaterPoints3D        (:,:,k)

            if(Me%ComputeOptions%Methodology == Bilinear_)then

                call StructureBilinear(Me%ExternalVar%FatherField, Me%LocalFatherField, Me%LocalFatherGrid, k)

            elseif(Me%ComputeOptions%Methodology == Spline2D_)then

                call StructureSpline2D(Me%ExternalVar%FatherField, Me%LocalFatherField, Me%LocalFatherGrid)

            elseif(Me%ComputeOptions%Methodology == Triangulation_)then

                call StructureTriangulation (Me%ExternalVar%FatherField, Me%LocalFatherField, Me%LocalFatherGrid)

            elseif(Me%ComputeOptions%Methodology == InverseWeight_)then

                call StructureInverseWeight (Me%ExternalVar%FatherField, Me%LocalFatherField, Me%LocalFatherGrid)

            else

                write(*,*) 'Unknown type of interpolation'
                stop       'StartInterpolateGrids - ModuleInterpolateGrids - ERR250' 

            end if

            Me%LocalFatherField%Values3D       (:,:,k) = Me%LocalFatherField%Values2D      (:,:)

        enddo

        !Horizontal Extrapolation procedure
        if (Me%ComputeOptions%Extrapolate2D /= NotActive_) then

            allocate(InValues2D (Me%ExternalVar%FatherGrid%Size2D%ILB:Me%ExternalVar%FatherGrid%Size2D%IUB,&
                                 Me%ExternalVar%FatherGrid%Size2D%JLB:Me%ExternalVar%FatherGrid%Size2D%JUB))
            allocate(OutValues2D(Me%LocalFatherGrid%Size2D%ILB       :Me%LocalFatherGrid%Size2D%IUB,       &
                                 Me%LocalFatherGrid%Size2D%JLB       :Me%LocalFatherGrid%Size2D%JUB))

            do k=Me%ExternalVar%FatherGrid%WorkSize3D%KLB,Me%ExternalVar%FatherGrid%WorkSize3D%KUB

                InValues2D (:,:) = Me%ExternalVar%FatherField%Values3D(:,:,k)
                OutValues2D(:,:) = Me%LocalFatherField%Values3D       (:,:,k)
        
                if      (Me%ComputeOptions%Extrapolate2D == MediumTriang) then
        
                    call ExtraPol2DFieldsTriang(Me%ExternalVar%FatherField, Me%LocalFatherGrid, .false., &
                                                InValues2D, OutValues2D, k)
       
                else if (Me%ComputeOptions%Extrapolate2D == HighTriang) then

                    call ExtraPol2DFieldsTriang(Me%ExternalVar%FatherField, Me%LocalFatherGrid, .true.,&
                                                InValues2D, OutValues2D, k)

                else if (Me%ComputeOptions%Extrapolate2D == NearestNeighbor) then

                    call ExtraPol2DFieldsNearest(Me%ExternalVar%FatherField, Me%LocalFatherGrid, &
                                                 InValues2D, OutValues2D, k)


                endif

                Me%LocalFatherField%Values3D       (:,:,k) = OutValues2D(:,:)

            enddo
            
            deallocate(InValues2D )
            deallocate(OutValues2D)


        endif


        !From local FatherGrid grid to SonGrid grid
        call VerticalInterpolation(Me%LocalFatherField, Me%ExternalVar%SonField )


    end subroutine ClassicInterpolaStruct3D

    !----------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ClassicInterpolaStruct2D

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------


        if(Me%ComputeOptions%Methodology == Bilinear_)then

            call StructureBilinear      (Me%ExternalVar%FatherField, Me%ExternalVar%SonField, Me%ExternalVar%SonGrid)

        elseif(Me%ComputeOptions%Methodology == Spline2D_)then

            call StructureSpline2D      (Me%ExternalVar%FatherField, Me%ExternalVar%SonField, Me%ExternalVar%SonGrid)

        elseif(Me%ComputeOptions%Methodology == Triangulation_)then

            call StructureTriangulation (Me%ExternalVar%FatherField, Me%ExternalVar%SonField, Me%ExternalVar%SonGrid)

        elseif(Me%ComputeOptions%Methodology == InverseWeight_)then

            call StructureInverseWeight (Me%ExternalVar%FatherField, Me%ExternalVar%SonField, Me%ExternalVar%SonGrid)

        else

            write(*,*) 'Unknown type of interpolation'
            stop       'ClassicInterpolaStruct2D - ModuleInterpolation - ERR250' 

        endif

        if (Me%ComputeOptions%Extrapolate2D /= NotActive_) then
        
            if      (Me%ComputeOptions%Extrapolate2D == MediumTriang) then
        
                call ExtraPol2DFieldsTriang(Me%ExternalVar%FatherField, Me%ExternalVar%SonGrid, .false., &
                                            Me%ExternalVar%FatherField%Values2D, Me%ExternalVar%SonField%Values2D, &
                                            Me%ExternalVar%FatherGrid%WorkSize3D%KUB)
       
            else if (Me%ComputeOptions%Extrapolate2D == HighTriang) then

                call ExtraPol2DFieldsTriang(Me%ExternalVar%FatherField, Me%ExternalVar%SonGrid, .true.,&
                                            Me%ExternalVar%FatherField%Values2D, Me%ExternalVar%SonField%Values2D, &
                                            Me%ExternalVar%FatherGrid%WorkSize3D%KUB)

            else if (Me%ComputeOptions%Extrapolate2D == NearestNeighbor) then

                call ExtraPol2DFieldsNearest(Me%ExternalVar%FatherField, Me%ExternalVar%SonGrid, &
                                             Me%ExternalVar%FatherField%Values2D, Me%ExternalVar%SonField%Values2D, &
                                             Me%ExternalVar%FatherGrid%WorkSize3D%KUB)


            endif

        endif




    end subroutine ClassicInterpolaStruct2D

    !----------------------------------------------------------------------------

    !--------------------------------------------------------------------------


    subroutine VerticalInterpolation(LocalFatherField, SonField )

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                  :: SonField, LocalFatherField

        !Local-------------------------------------------------------------
        !real,    dimension(:,:,:), pointer      :: Values3D
        real,    dimension(:,:,:), pointer      :: ZCellCenter
        real,    dimension(:,:,:), pointer      :: ZCellCenterLocalFather
        real(8), dimension(:    ), pointer      :: Depth, Values
        real(8)                                 :: INDepth, dz, Error
        integer                                 :: i,j,k, PoliDegree, NDEPTHS, Aux, STAT_CALL
        logical                                 :: PoliIsEven, FoundBottom, FoundSurface

        !Begin-----------------------------------------------------------------

        allocate(Depth (Me%ExternalVar%FatherGrid%WorkSize3D%KLB: Me%ExternalVar%FatherGrid%WorkSize3D%KUB))
        allocate(Values(Me%ExternalVar%FatherGrid%WorkSize3D%KLB: Me%ExternalVar%FatherGrid%WorkSize3D%KUB))

        call GetGeometryDistances(Me%ExternalVar%SonGrid%ObjGeometry, ZCellCenter = ZCellCenter, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalInterpolation - ModuleInterpolateGrids - ERR10'

        call GetGeometryDistances(Me%LocalFatherGrid%ObjGeometry, ZCellCenter = ZCellCenterLocalFather, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalInterpolation - ModuleInterpolateGrids - ERR20'

        do i = Me%ExternalVar%SonGrid%WorkSize3D%ILB, Me%ExternalVar%SonGrid%WorkSize3D%IUB
        do j = Me%ExternalVar%SonGrid%WorkSize3D%JLB, Me%ExternalVar%SonGrid%WorkSize3D%JUB

            if (Me%ExternalVar%SonGrid%WaterPoints3D(i, j, Me%ExternalVar%SonGrid%WorkSize3D%KUB   ) == 1) then
            
            if (Me%LocalFatherGrid%WaterPoints3D(i, j, Me%ExternalVar%FatherGrid%WorkSize3D%KUB) == 1) then

                Aux=Me%ExternalVar%FatherGrid%WorkSize3D%KLB

                do k=Me%ExternalVar%FatherGrid%WorkSize3D%KUB, Me%ExternalVar%FatherGrid%WorkSize3D%KLB, -1
                    dz = abs(ZCellCenterLocalFather(i,j, k) - ZCellCenterLocalFather(i,j, k-1))
                    if(Me%LocalFatherGrid%WaterPoints3D(i, j, k) == 1                               &
                    !In the Hycom model instead of having a bottom mapping to the bottom cells are 
                    !given the same depth
                    .and.  (dz < 0.0001 .or. dz > -FillvalueReal/10.)) then
                        Aux = k
                        exit
                    endif
                enddo

                Depth (Aux:Me%ExternalVar%FatherGrid%WorkSize3D%KUB) = &
                -  ZCellCenterLocalFather(i,j, Aux:Me%ExternalVar%FatherGrid%WorkSize3D%KUB)
                Values(Aux:Me%ExternalVar%FatherGrid%WorkSize3D%KUB) = &
                LocalFatherField%Values3D(i,j, Aux:Me%ExternalVar%FatherGrid%WorkSize3D%KUB)
        
                NDEPTHS = Me%ExternalVar%FatherGrid%WorkSize3D%KUB - Aux  + 1

                do k=Me%ExternalVar%SonGrid%WorkSize3D%KLB,Me%ExternalVar%SonGrid%WorkSize3D%KUB

                    if (Me%ExternalVar%SonGrid%WaterPoints3D(i, j, k) == 1) then

                        INDepth = - ZCellCenter (i,j,k)

                        if (Me%ComputeOptions%PoliDegreeVert == LinearProfile_) then

                            SonField%Values3D (i,j,k) = InterpolateProfileR8 (INDepth, NDEPTHS, &
                                                        Depth (Aux:Me%ExternalVar%FatherGrid%WorkSize3D%KUB),   &
                                                        Values(Aux:Me%ExternalVar%FatherGrid%WorkSize3D%KUB),   &
                                                        FoundBottom, FoundSurface)

                        else if (Me%ComputeOptions%PoliDegreeVert >= 2) then 


                            if ( NDEPTHS == 1) then
                                !Uniform profile is assumed when there only one layer
                                SonField%Values3D (i,j,k) = LocalFatherField%Values3D (i,j, k)
                            else
                                !Interpolation n degree
                                PoliDegree = min(Me%ComputeOptions%PoliDegreeVert, NDEPTHS-1)

                                if(IsOdd(PoliDegree))then
                                    PoliIsEven = .false.
                                else
                                    PoliIsEven = .true.
                                endif


                                SonField%Values3D (i,j,k) = PolIntProfile (INDepth, NDEPTHS,                     &
                                                                           Depth (Aux:Me%ExternalVar%FatherGrid%WorkSize3D%KUB), &
                                                                           Values(Aux:Me%ExternalVar%FatherGrid%WorkSize3D%KUB), &
                                                                           PoliDegree, PoliIsEven, Error)

                            endif

                        endif

                        if (.not. Me%ComputeOptions%ExtrapolateProfile) then

                            if (INDepth< Depth (Me%ExternalVar%FatherGrid%WorkSize3D%KUB)) SonField%Values3D (i,j,k) = FillValueReal
                            if (INDepth> Depth (Aux                     )) SonField%Values3D (i,j,k) = FillValueReal

                        endif

                    else

                        SonField%Values3D (i,j,k) = FillValueReal

                    endif

                enddo

            else

                stop 'VerticalInterpolation - ModuleInterpolateGrids - ERR30'                

            endif
            endif

        enddo
        enddo

        deallocate(Depth )
        deallocate(Values)

        call UnGetGeometry(Me%ExternalVar%SonGrid%ObjGeometry, ZCellCenter, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalInterpolation - ModuleInterpolateGrids - ERR40'

        call UnGetGeometry(Me%LocalFatherGrid%ObjGeometry, ZCellCenter, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalInterpolation - ModuleInterpolateGrids - ERR50'


    end subroutine VerticalInterpolation


    !------------------------------------------------------------------------

    subroutine StructureBilinear(FatherField, SonField, SonGrid, k)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                  :: FatherField, SonField
        type (T_Grid)                           :: SonGrid
        integer, optional                       :: k

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: k_

        !----------------------------------------------------------------------

        write(*,*)'Interpolating : '//trim(FatherField%Name), FatherField%IDNumber

        if (present(k)) then

            k_ = k

        else

            k_ = 1

        endif

        if(Me%ComputeOptions%n_dimensions == 3)then

            call InterpolRegularGrid(SonGrid%ObjHorizontalGrid,                         &
                                     Me%ExternalVar%FatherGrid%ObjHorizontalGrid,           &
                                     FatherField%Values2D,                              &
                                     SonField%Values2D,                                 &
                                     Compute = Me%ComputeOptions%TypeZUV,               &
                                     KUBFather = k_,                                    &
                                     ComputeFather = Me%ExternalVar%FatherGrid%WaterPoints3D,&
                                     STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)stop 'StructureBilinear - ModuleInterpolation - ERR10'

        else if(Me%ComputeOptions%n_dimensions == 2)then

            call InterpolRegularGrid(SonGrid%ObjHorizontalGrid,                         &
                                     Me%ExternalVar%FatherGrid%ObjHorizontalGrid,           &
                                     FatherField%Values2D,                              &
                                     SonField%Values2D,                                 &
                                     Compute = Me%ComputeOptions%TypeZUV,               &
                                     KUBFather = k_,                                    &
                                     STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)stop 'StructureBilinear - ModuleInterpolation - ERR20'

        else

            stop 'StructureBilinear - ModuleInterpolation - ERR30'

        end if


    end subroutine StructureBilinear

    
    !------------------------------------------------------------------------


    subroutine StructureSpline2D(FatherField, SonField, SonGrid)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                  :: FatherField, SonField
        type (T_Grid)                           :: SonGrid

        !Local-----------------------------------------------------------------
        real, dimension(:  ),       pointer             :: Father_XX_Z, Father_YY_Z
        real, dimension(:  ),       pointer             :: New_XX_Z, New_YY_Z
        real                                            :: Xorig, Yorig, GridAngle
        integer                                         :: STAT_CALL, i, j
        real, dimension(:,:),       pointer             :: function_int
        real, dimension(:  ),       pointer             :: Father_X, Father_Y
        real, dimension(:  ),       pointer             :: New_X, New_Y

        !----------------------------------------------------------------------
        
        call GetHorizontalGrid(Me%ExternalVar%FatherGrid%ObjHorizontalGrid, &
                               XX_Z = Father_XX_Z,          &
                               YY_Z = Father_YY_Z,          &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'StructureSpline2D - ModuleInterpolateGrids - ERR10'

        call GetGridOrigin(Me%ExternalVar%FatherGrid%ObjHorizontalGrid, Xorig, Yorig, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'StructureSpline2D - ModuleInterpolateGrids - ERR20'
        
        call GetGridAngle(Me%ExternalVar%FatherGrid%ObjHorizontalGrid, GridAngle, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'StructureSpline2D - ModuleInterpolateGrids - ERR30'

        allocate(Father_X(Me%ExternalVar%FatherGrid%Size2D%JLB:Me%ExternalVar%FatherGrid%Size2D%JUB))
        allocate(Father_Y(Me%ExternalVar%FatherGrid%Size2D%ILB:Me%ExternalVar%FatherGrid%Size2D%IUB))

        do i = Me%ExternalVar%FatherGrid%Size2D%ILB, Me%ExternalVar%FatherGrid%Size2D%IUB
        do j = Me%ExternalVar%FatherGrid%Size2D%JLB, Me%ExternalVar%FatherGrid%Size2D%JUB

            Father_X(j) = Father_XX_Z(j)
            Father_Y(i) = Father_YY_Z(i)

            call RodaXY(Xorig, Yorig, GridAngle, Father_X(j), Father_Y(i))

        end do
        end do

        allocate(New_X   (SonGrid%Size2D%JLB:SonGrid%Size2D%JUB))
        allocate(New_Y   (SonGrid%Size2D%ILB:SonGrid%Size2D%IUB))
        
        call GetGridOrigin(SonGrid%ObjHorizontalGrid, Xorig, Yorig, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'StructureSpline2D - ModuleInterpolateGrids - ERR40'
        
        call GetGridAngle(SonGrid%ObjHorizontalGrid, GridAngle, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'StructureSpline2D - ModuleInterpolateGrids - ERR50'

        call GetHorizontalGrid(SonGrid%ObjHorizontalGrid,                               &
                               XX_Z = New_XX_Z,                                         &
                               YY_Z = New_YY_Z,                                         &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'StructureSpline2D - ModuleInterpolateGrids - ERR60'

        do i = SonGrid%Size2D%ILB, SonGrid%Size2D%IUB
        do j = SonGrid%Size2D%JLB, SonGrid%Size2D%JUB

            New_X(j) = New_XX_Z(j)
            New_Y(i) = New_YY_Z(i)

            call RodaXY(Xorig, Yorig, GridAngle, New_X(j), New_Y(i))

        end do
        end do

        allocate(function_int(Me%ExternalVar%FatherGrid%Size2D%ILB:Me%ExternalVar%FatherGrid%Size2D%IUB,                &
                              Me%ExternalVar%FatherGrid%Size2D%JLB:Me%ExternalVar%FatherGrid%Size2D%JUB))


        write(*,*)'Interpolating : '//trim(FatherField%Name), FatherField%IDNumber

        call splie2(npmax   = max(Me%ExternalVar%FatherGrid%Size2D%IUB,Me%ExternalVar%FatherGrid%Size2D%JUB),   &
                    x1a     = Father_Y,                                                 &
                    x2a     = Father_X,                                                 &
                    ya      = FatherField%Values2D,                                     &
                    m       = Me%ExternalVar%FatherGrid%Size2D%IUB,                                 &
                    n       = Me%ExternalVar%FatherGrid%Size2D%JUB,                                 &
                    y2a     = function_int)
               
        do i = SonGrid%Size2D%ILB,  SonGrid%Size2D%IUB
        do j = SonGrid%Size2D%JLB , SonGrid%Size2D%JUB


            call splin2(npmax   = max(Me%ExternalVar%FatherGrid%Size2D%IUB,Me%ExternalVar%FatherGrid%Size2D%JUB),&
                        x1a     = Father_Y,                                             &
                        x2a     = Father_X,                                             &
                        ya      = FatherField%Values2D,                                 &
                        y2a     = function_int,                                         &
                        m       = Me%ExternalVar%FatherGrid%Size2D%IUB,                             &
                        n       = Me%ExternalVar%FatherGrid%Size2D%JUB,                             &
                        x1      = New_Y(i),                                             &
                        x2      = New_X(j),                                             &
                        yc      = SonField%Values2D(i,j))
           
        enddo
        enddo
        
        deallocate(Father_X, Father_Y, New_X, New_Y, function_int)
        nullify   (Father_X, Father_Y, New_X, New_Y, function_int)

        call UnGetHorizontalGrid(Me%ExternalVar%FatherGrid%ObjHorizontalGrid, Father_XX_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'StructureSpline2D - ModuleInterpolateGrids - ERR80'

        call UnGetHorizontalGrid(Me%ExternalVar%FatherGrid%ObjHorizontalGrid, Father_YY_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'StructureSpline2D - ModuleInterpolateGrids - ERR90'

        call UnGetHorizontalGrid(SonGrid%ObjHorizontalGrid, Father_XX_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'StructureSpline2D - ModuleInterpolateGrids - ERR100'

        call UnGetHorizontalGrid(SonGrid%ObjHorizontalGrid, Father_YY_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'StructureSpline2D - ModuleInterpolateGrids - ERR110'


    end subroutine StructureSpline2D
    
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------


    subroutine splie2(npmax,x1a,x2a,ya,m,n,y2a)
        
        !Arguments-------------------------------------------------------------
        integer                         :: npmax, m,n
        real, dimension(:,:), pointer   :: y2a
        real, dimension(:,:), pointer   :: ya
        real, dimension(:  ), pointer   :: x2a
        real, dimension(:  ), pointer   :: x1a

        !Local-----------------------------------------------------------------
        real                            :: y2tmp(npmax),ytmp(npmax)
        integer                         :: j, k
        
        !Begin-----------------------------------------------------------------
        
        !Just to avoid compiler warning
        x1a = x1a

        do j = 1, m
            
            do k = 1, n
                ytmp(k)=ya(j,k)
            end do

            call spline(npmax, x2a, ytmp, n, 1.e30, 1.e30, y2tmp)
            
            do k=1, n
                y2a(j,k) = y2tmp(k)
            enddo

        enddo
      
    end subroutine splie2

    !------------------------------------------------------------------------
    
    subroutine spline(npmax,xb,yb,n,yp1,ypn,y2)
        
        !Arguments-------------------------------------------------------------
        integer                     :: npmax, n
        real                        :: yp1,ypn,yb(npmax),y2(npmax)

        
        !Local-----------------------------------------------------------------
        real, dimension(:), pointer :: xb
        real                        :: p,qn,sig,un,uu(npmax)
        integer                     :: i, k
        !Begin-----------------------------------------------------------------

        if (yp1.gt..99e30) then
            y2(1) = 0.
            uu(1) = 0.
        else
            y2(1) = -0.5
            uu(1) = (3./(xb(2)-xb(1)))*((yb(2)-yb(1))/(xb(2)-xb(1))-yp1)
        endif
      
        do i=2,n-1
            sig     = (xb(i)-xb(i-1))/(xb(i+1)-xb(i-1))
            p       = sig*y2(i-1)+2.
            y2(i)   = (sig-1.)/p
            uu(i)   = (6.*((yb(i+1)-yb(i))/(xb(i+1)-xb(i))-(yb(i)-yb(i-1))/ &
                      (xb(i)-xb(i-1)))/(xb(i+1)-xb(i-1))-sig* uu(i-1))/p
        enddo


        if (ypn.gt..99e30) then
            qn = 0.
            un = 0.
        else
            qn = 0.5
            un = (3./(xb(n)-xb(n-1)))*(ypn-(yb(n)-yb(n-1))/(xb(n)-xb(n-1)))
        endif

        y2(n) =(un-qn*uu(n-1))/(qn*y2(n-1)+1.)
        
        do k=n-1,1,-1
            y2(k) = y2(k)*y2(k+1)+uu(k)
        enddo

    end subroutine spline

    !------------------------------------------------------------------------

    subroutine splint_local(npmax,xa,ya,y2a,n,xc,yc)
        
        !Arguments-------------------------------------------------------------
        integer                     :: npmax, n
        real                        :: xc,yc
        real, dimension(:), pointer :: xa
        real                        :: y2a(npmax),ya(npmax)
        
        !Local-----------------------------------------------------------------
        real                        :: a,b,hh
        integer                     :: klo, khi,k

        !Begin-----------------------------------------------------------------


        klo = 1
        khi = n

        do while(khi-klo.gt.1)

            k = (khi + klo)/2
            
            if(xa(k).gt.xc)then
                khi = k
            else
                klo = k
            endif

        end do

        hh = xa(khi)-xa(klo)

        if (hh.eq.0.) stop 'bad xa input in splint_local - ModuleInterpolateGrids'
        
        a   = (xa(khi)-xc)/hh
        b   = (xc-xa(klo))/hh
        
        yc  = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(hh**2)/6.
        

    end subroutine splint_local

    !------------------------------------------------------------------------

    subroutine splin2(npmax,x1a,x2a,ya,y2a,m,n,x1,x2,yc)

        !Arguments-------------------------------------------------------------
        integer                         :: npmax, m, n
        real                            :: x1,x2,yc
        real, dimension(:,:), pointer   :: ya,y2a
        real, dimension(:  ), pointer   :: x2a
        real, dimension(:  ), pointer   :: x1a
        
        !Local-----------------------------------------------------------------
        real                            :: y2tmp(npmax),ytmp(npmax),yytmp(npmax)
        integer                         :: j, k

        !Begin-----------------------------------------------------------------

        do j=1,m
            do k=1,n
                ytmp(k) = ya(j,k)
                y2tmp(k)= y2a(j,k)
            enddo

            call splint_local(npmax,x2a,ytmp,y2tmp,n,x2,yytmp(j))
        enddo

        call spline(npmax,x1a,yytmp,m,1.e30,1.e30,y2tmp)
        call splint_local(npmax,x1a,yytmp,y2tmp,m,x1,yc)
           

    end subroutine splin2

    !----------------------------------------------------------------------------

    subroutine StructureTriangulation (FatherField, SonField, SonGrid)

        !Arguments-------------------------------------------------------------
        type(T_Field),           pointer            :: FatherField, SonField
        type(T_Grid )                               :: SonGrid

        !Local-----------------------------------------------------------------
        real,        dimension(:,:),   pointer      :: SonCenterX, SonCenterY
        real,        dimension(:),     pointer      :: NodeX, NodeY, NodeZ
        integer                                     :: STAT_CALL
        integer                                     :: NumberOfNodes, Count, i, j
        logical                                     :: FillOutsidePoints   = .false.
        !integer,    dimension(:,:  ), pointer       :: WaterPoints2D
        integer                                     :: ILB, IUB, JLB, JUB
        
        !Begin-----------------------------------------------------------------
        ILB = Me%ExternalVar%FatherGrid%Size2D%ILB
        IUB = Me%ExternalVar%FatherGrid%Size2D%IUB
        JLB = Me%ExternalVar%FatherGrid%Size2D%JLB
        JUB = Me%ExternalVar%FatherGrid%Size2D%JUB

        NumberOfNodes =  Sum(Me%ExternalVar%FatherGrid%WaterPoints2D(ILB:IUB, JLB:JUB))

iN:     if (NumberOfNodes >= 3) then

        allocate(NodeX(NumberOfNodes))
        allocate(NodeY(NumberOfNodes))
        allocate(NodeZ(NumberOfNodes))

        Count = 0

        do j = Me%ExternalVar%FatherGrid%Size2D%JLB, Me%ExternalVar%FatherGrid%Size2D%JUB
        do i = Me%ExternalVar%FatherGrid%Size2D%ILB, Me%ExternalVar%FatherGrid%Size2D%IUB

            if (Me%ExternalVar%FatherGrid%WaterPoints2D(i, j) == WaterPoint) then

                Count           = Count + 1

                NodeX(Count) = ((Me%ExternalVar%FatherGrid%XX_IE(i, j  ) + Me%ExternalVar%FatherGrid%XX_IE(i+1, j  ))/2. + &
                                (Me%ExternalVar%FatherGrid%XX_IE(i, j+1) + Me%ExternalVar%FatherGrid%XX_IE(i+1, j+1))/2.)/2.
        
                NodeY(Count) = ((Me%ExternalVar%FatherGrid%YY_IE(i, j  ) + Me%ExternalVar%FatherGrid%YY_IE(i+1, j  ))/2. + &
                                (Me%ExternalVar%FatherGrid%YY_IE(i, j+1) + Me%ExternalVar%FatherGrid%YY_IE(i+1, j+1))/2.)/2.

                NodeZ(Count) = FatherField%Values2D(i, j)
                

            endif

        enddo
        enddo

            !Constructs Triangulation
        call ConstructTriangulation (Me%ObjTriangulation,   &
                                     NumberOfNodes,         &
                                     NodeX,              &
                                     NodeY,              &
                                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Triangulator - ModuleInterpolateGrids - ERR10'

        call SetHeightValues(Me%ObjTriangulation, NodeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Triangulator - ModuleInterpolateGrids - ERR20'


        allocate(SonCenterX(SonGrid%Size2D%ILB:SonGrid%Size2D%IUB, &
                               SonGrid%Size2D%JLB:SonGrid%Size2D%JUB))
    
        allocate(SonCenterY(SonGrid%Size2D%ILB:SonGrid%Size2D%IUB, &
                               SonGrid%Size2D%JLB:SonGrid%Size2D%JUB))

        do j = SonGrid%Size2D%JLB, SonGrid%Size2D%JUB
        do i = SonGrid%Size2D%ILB, SonGrid%Size2D%IUB
                
            !Find SonGrid cell center
            SonCenterX(i,j) = ((SonGrid%XX_IE(i, j  ) + SonGrid%XX_IE(i+1, j  ))/2. + &
                               (SonGrid%XX_IE(i, j+1) + SonGrid%XX_IE(i+1, j+1))/2.)/2.
    
            SonCenterY(i,j) = ((SonGrid%YY_IE(i, j  ) + SonGrid%YY_IE(i+1, j  ))/2. + &
                               (SonGrid%YY_IE(i, j+1) + SonGrid%YY_IE(i+1, j+1))/2.)/2.

        enddo
        enddo


        do j = SonGrid%Size2D%JLB, SonGrid%Size2D%JUB
        do i = SonGrid%Size2D%ILB, SonGrid%Size2D%IUB

        if(SonGrid%WaterPoints2D(i, j) == WaterPoint) then

                SonField%Values2D(i, j) = InterPolation(Me%ObjTriangulation,            &
                                                        SonCenterX(i,j),             &
                                                        SonCenterY(i,j),             &
                                                        FillOutsidePoints,              &
                                                        Default = null_real,            &
                                                        STAT    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Triangulator - ModuleInterpolateGrids - ERR30'

        end if

        enddo
        enddo

        call KillTriangulation (Me%ObjTriangulation, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Triangulator - ModuleInterpolateGrids - ERR40'

        deallocate(SonCenterX)
        deallocate(SonCenterY)

        deallocate(NodeX)
        deallocate(NodeY)
        deallocate(NodeZ)

        endif iN


    end subroutine StructureTriangulation

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine UnStructureTriangulation (Points_XYP, SonField, SonGrid)

        !Arguments-------------------------------------------------------------
        type(T_XYZPoints),       pointer            :: Points_XYP
        type(T_Field),           pointer            :: SonField
        type(T_Grid )                               :: SonGrid

        !Local-----------------------------------------------------------------
        real,        dimension(:,:),   pointer      :: SonCenterX, SonCenterY
        real,        dimension(:),     pointer      :: NodeX, NodeY, NodeZ
        integer                                     :: STAT_CALL
        integer                                     :: NumberOfNodes, i, j
        logical                                     :: FillOutsidePoints   = .false.
        !integer,    dimension(:,:  ), pointer       :: WaterPoints2D
        
        !Begin-----------------------------------------------------------------


        NumberOfNodes =  Points_XYP%Count

iN:     if (NumberOfNodes >= 3) then

        NodeX => Points_XYP%X 
        NodeY => Points_XYP%Y
        !CUIDADO!!!!!!! 
        NodeZ => Points_XYP%Z


            !Constructs Triangulation
        call ConstructTriangulation (Me%ObjTriangulation,   &
                                     NumberOfNodes,         &
                                     NodeX,              &
                                     NodeY,              &
                                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Triangulator - ModuleInterpolateGrids - ERR10'

        call SetHeightValues(Me%ObjTriangulation, NodeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Triangulator - ModuleInterpolateGrids - ERR20'


        allocate(SonCenterX(SonGrid%Size2D%ILB:SonGrid%Size2D%IUB, &
                               SonGrid%Size2D%JLB:SonGrid%Size2D%JUB))
    
        allocate(SonCenterY(SonGrid%Size2D%ILB:SonGrid%Size2D%IUB, &
                               SonGrid%Size2D%JLB:SonGrid%Size2D%JUB))

        do j = SonGrid%Size2D%JLB, SonGrid%Size2D%JUB
        do i = SonGrid%Size2D%ILB, SonGrid%Size2D%IUB
                
            !Find SonGrid cell center
            SonCenterX(i,j) = ((SonGrid%XX_IE(i, j  ) + SonGrid%XX_IE(i+1, j  ))/2. + &
                                  (SonGrid%XX_IE(i, j+1) + SonGrid%XX_IE(i+1, j+1))/2.)/2.
    
            SonCenterY(i,j) = ((SonGrid%YY_IE(i, j  ) + SonGrid%YY_IE(i+1, j  ))/2. + &
                                  (SonGrid%YY_IE(i, j+1) + SonGrid%YY_IE(i+1, j+1))/2.)/2.

        enddo
        enddo

        do j = SonGrid%Size2D%JLB, SonGrid%Size2D%JUB
        do i = SonGrid%Size2D%ILB, SonGrid%Size2D%IUB

        if(SonGrid%WaterPoints2D(i, j) == WaterPoint) then

            SonField%Values2D(i, j) = InterPolation(Me%ObjTriangulation,                &
                                                    SonCenterX(i,j),                    &
                                                    SonCenterY(i,j),                    &
                                                    FillOutsidePoints,                  &
                                                    Default = null_real,                &
                                                    STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Triangulator - ModuleInterpolateGrids - ERR30'

        end if

        enddo
        enddo

        call KillTriangulation (Me%ObjTriangulation, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Triangulator - ModuleInterpolateGrids - ERR40'

        deallocate(SonCenterX)
        deallocate(SonCenterY)

        nullify(NodeX, NodeY, NodeZ)

        endif iN


    end subroutine UnStructureTriangulation

    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------

    subroutine StructureInverseWeight (FatherField, SonField, SonGrid)

        !Arguments-------------------------------------------------------------
        type(T_Field),           pointer            :: FatherField, SonField
        type(T_Grid )                               :: SonGrid

        !Local-----------------------------------------------------------------
        real,        dimension(:,:),   pointer      :: SonCenterX, SonCenterY
        real,        dimension(:),     pointer      :: NodeX, NodeY, NodeZ

        integer                                     :: NumberOfNodes, Count, i, j, ip
        integer                                     :: Denominator, Nominator, dist
            !integer,    dimension(:,:  ), pointer       :: WaterPoints2D
        
        !Begin-----------------------------------------------------------------


        NumberOfNodes =  Sum(Me%ExternalVar%FatherGrid%WaterPoints2D                        &
                             (Me%ExternalVar%FatherGrid%Size2D%ILB:Me%ExternalVar%FatherGrid%Size2D%IUB, &
                              Me%ExternalVar%FatherGrid%Size2D%JLB:Me%ExternalVar%FatherGrid%Size2D%JUB))

iN:     if (NumberOfNodes >= 3) then

        allocate(NodeX(NumberOfNodes))
        allocate(NodeY(NumberOfNodes))
        allocate(NodeZ(NumberOfNodes))

        Count = 0

        do j = Me%ExternalVar%FatherGrid%Size2D%JLB, Me%ExternalVar%FatherGrid%Size2D%JUB
        do i = Me%ExternalVar%FatherGrid%Size2D%ILB, Me%ExternalVar%FatherGrid%Size2D%IUB

            if (Me%ExternalVar%FatherGrid%WaterPoints2D(i, j) == WaterPoint) then

                Count           = Count + 1

                NodeX(Count) = ((Me%ExternalVar%FatherGrid%XX_IE(i, j  ) + Me%ExternalVar%FatherGrid%XX_IE(i+1, j  ))/2. + &
                                (Me%ExternalVar%FatherGrid%XX_IE(i, j+1) + Me%ExternalVar%FatherGrid%XX_IE(i+1, j+1))/2.)/2.
        
                NodeY(Count) = ((Me%ExternalVar%FatherGrid%YY_IE(i, j  ) + Me%ExternalVar%FatherGrid%YY_IE(i+1, j  ))/2. + &
                                (Me%ExternalVar%FatherGrid%YY_IE(i, j+1) + Me%ExternalVar%FatherGrid%YY_IE(i+1, j+1))/2.)/2.

                NodeZ(Count) = FatherField%Values2D(i, j)
                

            endif

        enddo
        enddo


        allocate(SonCenterX(SonGrid%Size2D%ILB:SonGrid%Size2D%IUB, &
                            SonGrid%Size2D%JLB:SonGrid%Size2D%JUB))
    
        allocate(SonCenterY(SonGrid%Size2D%ILB:SonGrid%Size2D%IUB, &
                            SonGrid%Size2D%JLB:SonGrid%Size2D%JUB))

        do j = SonGrid%Size2D%JLB, SonGrid%Size2D%JUB
        do i = SonGrid%Size2D%ILB, SonGrid%Size2D%IUB
                
            !Find SonGrid cell center
            SonCenterX(i,j) = ((SonGrid%XX_IE(i, j  ) + SonGrid%XX_IE(i+1, j  ))/2. + &
                               (SonGrid%XX_IE(i, j+1) + SonGrid%XX_IE(i+1, j+1))/2.)/2.
    
            SonCenterY(i,j) = ((SonGrid%YY_IE(i, j  ) + SonGrid%YY_IE(i+1, j  ))/2. + &
                               (SonGrid%YY_IE(i, j+1) + SonGrid%YY_IE(i+1, j+1))/2.)/2.

        enddo
        enddo


        do j = SonGrid%Size2D%JLB, SonGrid%Size2D%JUB
        do i = SonGrid%Size2D%ILB, SonGrid%Size2D%IUB

            if(SonGrid%WaterPoints2D(i, j) == WaterPoint) then

                Nominator   = 0.0
                Denominator = 0.0
DoPoints:       do ip = 1, NumberOfNodes
                    dist = sqrt((SonCenterX(i,j) - NodeX(ip))**2.0 + &
                                (SonCenterY(i,j) - NodeY(ip))**2.0)

                    if (dist > 0.0) then    
                        if (dist < Me%ComputeOptions%MaxDistance) then 
                            Nominator   = Nominator   + NodeZ(ip) / (dist ** Me%ComputeOptions%IWDn)
                            Denominator = Denominator + 1.0 / (dist ** Me%ComputeOptions%IWDn)
                        endif
                    else
                        Nominator   = NodeZ(ip)
                        Denominator = 1.0
                        exit DoPoints
                    endif
                enddo DoPoints

                if (Denominator .lt. AllMostZero) then
                    write (*,*)'Insufficient data avaliable'
                    write (*,*)'Increase MAX_DISTANCE or get more data ;-)'
                    write (*,*)'Point [i, j]', i, j
                    stop 'StructureInverseWeight - ModuleInterpolation - ERR10'
                endif

                SonField%Values2D(i, j) =  Nominator / Denominator

            end if

        enddo
        enddo


        deallocate(SonCenterX)
        deallocate(SonCenterY)

        deallocate(NodeX)
        deallocate(NodeY)
        deallocate(NodeZ)

        endif iN


    end subroutine StructureInverseWeight

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine UnStructureInverseWeight (Points_XYP, SonField, SonGrid)

        !Arguments-------------------------------------------------------------
        type(T_XYZPoints),       pointer            :: Points_XYP
        type(T_Field),           pointer            :: SonField
        type(T_Grid )                               :: SonGrid

        !Local-----------------------------------------------------------------
        real,        dimension(:,:),   pointer      :: SonCenterX, SonCenterY
        real,        dimension(:),     pointer      :: NodeX, NodeY, NodeZ
        integer                                     :: NumberOfNodes, i, j, ip
        integer                                     :: Denominator, Nominator, dist
            !integer,    dimension(:,:  ), pointer       :: WaterPoints2D
        
        !Begin-----------------------------------------------------------------



        NumberOfNodes =  Points_XYP%Count

iN:     if (NumberOfNodes >= 3) then

        NodeX => Points_XYP%X 
        NodeY => Points_XYP%Y
        !CUIDADO!!!!!!! 
        NodeZ => Points_XYP%Z

        allocate(SonCenterX(SonGrid%Size2D%ILB:SonGrid%Size2D%IUB, &
                            SonGrid%Size2D%JLB:SonGrid%Size2D%JUB))
    
        allocate(SonCenterY(SonGrid%Size2D%ILB:SonGrid%Size2D%IUB, &
                            SonGrid%Size2D%JLB:SonGrid%Size2D%JUB))

        do j = SonGrid%Size2D%JLB, SonGrid%Size2D%JUB
        do i = SonGrid%Size2D%ILB, SonGrid%Size2D%IUB
                
            !Find SonGrid cell center
            SonCenterX(i,j) = ((SonGrid%XX_IE(i, j  ) + SonGrid%XX_IE(i+1, j  ))/2. + &
                               (SonGrid%XX_IE(i, j+1) + SonGrid%XX_IE(i+1, j+1))/2.)/2.
    
            SonCenterY(i,j) = ((SonGrid%YY_IE(i, j  ) + SonGrid%YY_IE(i+1, j  ))/2. + &
                               (SonGrid%YY_IE(i, j+1) + SonGrid%YY_IE(i+1, j+1))/2.)/2.

        enddo
        enddo


        do j = SonGrid%Size2D%JLB, SonGrid%Size2D%JUB
        do i = SonGrid%Size2D%ILB, SonGrid%Size2D%IUB

            if(SonGrid%WaterPoints2D(i, j) == WaterPoint) then

                Nominator   = 0.0
                Denominator = 0.0
DoPoints:       do ip = 1, NumberOfNodes
                    dist = sqrt((SonCenterX(i,j) - NodeX(ip))**2.0 + &
                                (SonCenterY(i,j) - NodeY(ip))**2.0)

                    if (dist > 0.0) then    
                        if (dist < Me%ComputeOptions%MaxDistance) then 
                            Nominator   = Nominator   + NodeZ(ip) / (dist ** Me%ComputeOptions%IWDn)
                            Denominator = Denominator + 1.0 / (dist ** Me%ComputeOptions%IWDn)
                        endif
                    else
                        Nominator   = NodeZ(ip)
                        Denominator = 1.0
                        exit DoPoints
                    endif
                enddo DoPoints

                if (Denominator .lt. AllMostZero) then
                    write (*,*)'Insufficient data avaliable'
                    write (*,*)'Increase MAX_DISTANCE or get more data ;-)'
                    write (*,*)'Point [i, j]', i, j
                    stop 'StructureInverseWeight - ModuleInterpolation - ERR10'
                endif

                SonField%Values2D(i, j) =  Nominator / Denominator

            end if

        enddo
        enddo


        deallocate(SonCenterX)
        deallocate(SonCenterY)

        nullify(NodeX)
        nullify(NodeY)
        nullify(NodeZ)

        endif iN


    end subroutine UnStructureInverseWeight

    !----------------------------------------------------------------------------

    subroutine ExtraPol2DFieldsTriang (FatherField, NewGrid, FillOutsidePoints, &
                                       InValues2D, OutValues2D, k)

        !Arguments-------------------------------------------------------------
        type(T_Field),           pointer            :: FatherField
        type(T_Grid )                               :: NewGrid
        logical                                     :: FillOutsidePoints
        real,         dimension(:,:), pointer       :: InValues2D, OutValues2D
        integer                                     :: k

        !Local-----------------------------------------------------------------
        real,        dimension(:,:),   pointer      :: SonCenterX, SonCenterY
        real,        dimension(:),     pointer      :: NodeX, NodeY, NodeZ
        integer                                     :: STAT_CALL
        integer                                     :: NumberOfNodes, Count, i, j
        real,       dimension(:,:),   pointer       :: LatitudeZ, LongitudeZ
        integer                                     :: ILB, IUB, JLB, JUB
        !Begin-----------------------------------------------------------------

        
        !Begin-----------------------------------------------------------------
        ILB = Me%ExternalVar%FatherGrid%Size2D%ILB
        IUB = Me%ExternalVar%FatherGrid%Size2D%IUB
        JLB = Me%ExternalVar%FatherGrid%Size2D%JLB
        JUB = Me%ExternalVar%FatherGrid%Size2D%JUB


        NumberOfNodes =  Sum(Me%ExternalVar%FatherGrid%WaterPoints3D(ILB:IUB,JLB:JUB,k))


iN:     if (NumberOfNodes > 3) then

        allocate(NodeX(NumberOfNodes))
        allocate(NodeY(NumberOfNodes))
        allocate(NodeZ(NumberOfNodes))

        Count = 0

        do j = Me%ExternalVar%FatherGrid%WorkSize3D%JLB, Me%ExternalVar%FatherGrid%WorkSize3D%JUB
        do i = Me%ExternalVar%FatherGrid%WorkSize3D%ILB, Me%ExternalVar%FatherGrid%WorkSize3D%IUB

            if (Me%ExternalVar%FatherGrid%WaterPoints3D(i, j, k) == WaterPoint) then

                Count           = Count + 1

                NodeX(Count) = ((Me%ExternalVar%FatherGrid%XX_IE(i, j  ) + Me%ExternalVar%FatherGrid%XX_IE(i+1, j  ))/2. + &
                                   (Me%ExternalVar%FatherGrid%XX_IE(i, j+1) + Me%ExternalVar%FatherGrid%XX_IE(i+1, j+1))/2.)/2.
        
                NodeY(Count) = ((Me%ExternalVar%FatherGrid%YY_IE(i, j  ) + Me%ExternalVar%FatherGrid%YY_IE(i+1, j  ))/2. + &
                                   (Me%ExternalVar%FatherGrid%YY_IE(i, j+1) + Me%ExternalVar%FatherGrid%YY_IE(i+1, j+1))/2.)/2.

                NodeZ(Count) = InValues2D(i, j)
                

            endif

        enddo
        enddo

            !Constructs Triangulation
        call ConstructTriangulation (Me%ObjTriangulation,   &
                                     NumberOfNodes,         &
                                     NodeX,                 &
                                     NodeY,                 &
                                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR10'

        call SetHeightValues(Me%ObjTriangulation, NodeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR20'


        allocate(SonCenterX(NewGrid%WorkSize3D%ILB:NewGrid%WorkSize3D%IUB, &
                               NewGrid%WorkSize3D%JLB:NewGrid%WorkSize3D%JUB))
    
        allocate(SonCenterY(NewGrid%WorkSize3D%ILB:NewGrid%WorkSize3D%IUB, &
                               NewGrid%WorkSize3D%JLB:NewGrid%WorkSize3D%JUB))

        call GetGridLatitudeLongitude(NewGrid%ObjHorizontalGrid,                        &
                                      GridLongitude = LongitudeZ,                       &
                                      GridLatitude  = LatitudeZ,                        &
                                      STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR30'

        do j = NewGrid%WorkSize3D%JLB, NewGrid%WorkSize3D%JUB
        do i = NewGrid%WorkSize3D%ILB, NewGrid%WorkSize3D%IUB

            !Find SonGrid cell center
            SonCenterX(i,j)   = LongitudeZ(i,j)
            SonCenterY(i,j)   = LatitudeZ(i,j)

        enddo
        enddo

        call UngetHorizontalGrid(NewGrid%ObjHorizontalGrid, LongitudeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR40'
            
        call UngetHorizontalGrid(NewGrid%ObjHorizontalGrid, LatitudeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR50'

        write(*,*)'Extrapolating : '//trim(FatherField%Name), FatherField%IDNumber

        do j = NewGrid%WorkSize3D%JLB, NewGrid%WorkSize3D%JUB
        do i = NewGrid%WorkSize3D%ILB, NewGrid%WorkSize3D%IUB

            if( NewGrid%WaterPoints3D(i, j, k) /= WaterPoint .or. &
               (NewGrid%WaterPoints3D(i, j, k) == WaterPoint .and. OutValues2D(i, j) == FillValueReal)) then

                OutValues2D(i, j) = InterPolation(Me%ObjTriangulation,                  &
                                                        SonCenterX(i,j),             &
                                                        SonCenterY(i,j),             &
                                                        FillOutsidePoints,              &
                                                        Default = null_real,            &
                                                        STAT    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR60'

            end if

        enddo
        enddo

        call KillTriangulation (Me%ObjTriangulation, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR70'

        deallocate(SonCenterX)
        deallocate(SonCenterY)

        deallocate(NodeX)
        deallocate(NodeY)
        deallocate(NodeZ)

        endif iN


    end subroutine ExtraPol2DFieldsTriang



    subroutine ExtraPol2DFieldsNearest (FatherField, NewGrid, InValues2D, OutValues2D, k)

        !Arguments-------------------------------------------------------------
        type(T_Field),           pointer            :: FatherField
        type(T_Grid )                               :: NewGrid
        real,         dimension(:,:), pointer       :: InValues2D, OutValues2D
        integer                                     :: k
        !logical                                     :: FillOutsidePoints

        !Local-----------------------------------------------------------------
        real,        dimension(:,:),   pointer      :: SonCenterX, SonCenterY
        real,        dimension(:),     pointer      :: NodeX, NodeY, NodeZ
        integer                                     :: STAT_CALL
        integer                                     :: NumberOfNodes, Count, i, j
        real                                        :: Distance, Aux
        real,       dimension(:,:),   pointer       :: LatitudeZ, LongitudeZ
        integer                                     :: ILB, IUB, JLB, JUB
        
        !Begin-----------------------------------------------------------------
        
        ILB = Me%ExternalVar%FatherGrid%Size2D%ILB
        IUB = Me%ExternalVar%FatherGrid%Size2D%IUB
        JLB = Me%ExternalVar%FatherGrid%Size2D%JLB
        JUB = Me%ExternalVar%FatherGrid%Size2D%JUB

        NumberOfNodes =  Sum(Me%ExternalVar%FatherGrid%WaterPoints3D(ILB:IUB,JLB:JUB,k))


iN:     if (NumberOfNodes > 0) then

        allocate(NodeX(NumberOfNodes))
        allocate(NodeY(NumberOfNodes))
        allocate(NodeZ(NumberOfNodes))

        Count = 0

        do j = Me%ExternalVar%FatherGrid%WorkSize3D%JLB, Me%ExternalVar%FatherGrid%WorkSize3D%JUB
        do i = Me%ExternalVar%FatherGrid%WorkSize3D%ILB, Me%ExternalVar%FatherGrid%WorkSize3D%IUB

            if (Me%ExternalVar%FatherGrid%WaterPoints3D(i, j, k) == WaterPoint) then

                Count           = Count + 1

                NodeX(Count) = ((Me%ExternalVar%FatherGrid%XX_IE(i, j  ) + Me%ExternalVar%FatherGrid%XX_IE(i+1, j  ))/2. + &
                                   (Me%ExternalVar%FatherGrid%XX_IE(i, j+1) + Me%ExternalVar%FatherGrid%XX_IE(i+1, j+1))/2.)/2.
        
                NodeY(Count) = ((Me%ExternalVar%FatherGrid%YY_IE(i, j  ) + Me%ExternalVar%FatherGrid%YY_IE(i+1, j  ))/2. + &
                                   (Me%ExternalVar%FatherGrid%YY_IE(i, j+1) + Me%ExternalVar%FatherGrid%YY_IE(i+1, j+1))/2.)/2.

                NodeZ(Count) = InValues2D(i, j)
                

            endif

        enddo
        enddo


        allocate(SonCenterX(NewGrid%WorkSize3D%ILB:NewGrid%WorkSize3D%IUB, &
                               NewGrid%WorkSize3D%JLB:NewGrid%WorkSize3D%JUB))
    
        allocate(SonCenterY(NewGrid%WorkSize3D%ILB:NewGrid%WorkSize3D%IUB, &
                               NewGrid%WorkSize3D%JLB:NewGrid%WorkSize3D%JUB))

        call GetGridLatitudeLongitude(NewGrid%ObjHorizontalGrid,                        &
                                      GridLongitude = LongitudeZ,                       &
                                      GridLatitude  = LatitudeZ,                        &
                                      STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsNearest - ModuleInterpolateGrids - ERR30'

        do j = NewGrid%WorkSize3D%JLB, NewGrid%WorkSize3D%JUB
        do i = NewGrid%WorkSize3D%ILB, NewGrid%WorkSize3D%IUB

                !Find SonGrid cell center
                SonCenterX(i,j)   = LongitudeZ(i,j)
                SonCenterY(i,j)   = LatitudeZ(i,j)

        enddo
        enddo

        call UngetHorizontalGrid(NewGrid%ObjHorizontalGrid, LongitudeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsNearest - ModuleInterpolateGrids - ERR40'
            
        call UngetHorizontalGrid(NewGrid%ObjHorizontalGrid, LatitudeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsNearest - ModuleInterpolateGrids - ERR50'


        write(*,*)'Extrapolating : '//trim(FatherField%Name), FatherField%IDNumber


        do j = NewGrid%WorkSize3D%JLB, NewGrid%WorkSize3D%JUB
        do i = NewGrid%WorkSize3D%ILB, NewGrid%WorkSize3D%IUB

        if ( NewGrid%WaterPoints3D(i, j, k) /= WaterPoint  .or. &
            (NewGrid%WaterPoints3D(i, j, k) == WaterPoint .and. OutValues2D(i, j) == FillValueReal)) then

            Aux = - FillValueReal

            do Count = 1, NumberOfNodes

                Distance = sqrt((NodeX(Count) - SonCenterX(i,j))**2 +                &
                                (NodeY(Count) - SonCenterY(i,j))**2)

                !Distance = abs(NodeX(Count) - SonCenterX(i,j))       +            &
                !           abs(NodeY(Count) - SonCenterY(i,j))


                if (Distance < Aux .and. NodeZ(Count) /= FillValueReal) then

                    OutValues2D(i, j) = NodeZ(Count)

                    Aux = Distance

                endif

            enddo

        end if

        enddo
        enddo

        deallocate(SonCenterX)
        deallocate(SonCenterY)

        deallocate(NodeX)
        deallocate(NodeY)
        deallocate(NodeZ)

        endif iN


    end subroutine ExtraPol2DFieldsNearest

    !---------------------------------------------------------------------------

    subroutine sorter (vec,res,dim,origin,destination)
        
        integer              :: dim
        real                    :: vec(dim),res(dim)
        integer              :: origin(dim), destination(dim)

        integer              :: aux0,aux2,a,b
        real(8)                 :: aux1

        do a=1,dim
            origin(a)=a
        enddo

        res(:)=vec(:)
        do a=1,dim-1
            aux0=a
            aux1=res(a)
            do b=a+1,dim
                if(res(b).lt.aux1) then
                    aux0=b
                    aux1=res(b)
                endif
            enddo

            aux2=origin(a)
            origin(a)=origin(aux0)
            origin(aux0)=aux2
            destination(origin(a))=a

            aux1=res(a)
            res(a)=res(aux0)
            res(aux0)=aux1
        enddo
        aux2=origin(dim)
        origin(dim)=origin(dim)
        origin(dim)=aux2
        destination(origin(dim))=dim

    end subroutine sorter

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine Points2GridSmooth_3D

        integer                                  :: n_dimensions, n_obs
        integer                                  :: a,b,d,ri
        integer                                  :: i,j,k
        integer                                  :: ilb,jlb,klb
        integer                                  :: iub,jub,kub
        real(8)                                     :: DistSum,EE,phi
        real(8)                                     :: coord1(3),coord2(3)
        real,pointer,dimension(:)                   :: X,Y,Z,w
        real,pointer,dimension(:)                   :: tauM
        real,pointer,dimension(:,:)                 :: CenterX,CenterY
        real,pointer,dimension(:,:)                 :: tauR,tau
        real,pointer,dimension(:,:,:)               :: CenterZ,v
        real(8),allocatable,dimension(:,:)          :: G,G2,MinDist
        real(8),allocatable,dimension(:)            :: distance

        n_dimensions    = Me%ComputeOptions%n_dimensions

        ilb             = Me%ExternalVar%SonGrid%WorkSize3D%ILB
        jlb             = Me%ExternalVar%SonGrid%WorkSize3D%JLB
        klb             = Me%ExternalVar%SonGrid%WorkSize3D%KLB
        iub             = Me%ExternalVar%SonGrid%WorkSize3D%IUB
        jub             = Me%ExternalVar%SonGrid%WorkSize3D%JUB
        kub             = Me%ExternalVar%SonGrid%WorkSize3D%KUB

        CenterX         => Me%ExternalVar%SonGrid%CenterX
        CenterY         => Me%ExternalVar%SonGrid%CenterY
        CenterZ         => Me%ExternalVar%SonGrid%CenterZ

        X               => Me%ExternalVar%TXYZP_Points%X
        Y               => Me%ExternalVar%TXYZP_Points%Y
        Z               => Me%ExternalVar%TXYZP_Points%Z
        !CUIDADO!!!!!!!!
        w               => Me%ExternalVar%TXYZP_Points%Z
        n_obs           =  Me%ExternalVar%TXYZP_Points%Count

        tau             => Me%NC%tau
        tauM            => Me%NC%tauM
        tauR            => Me%NC%tauR
        phi             =  Me%ComputeOptions%phi

        allocate(MinDist(n_obs,n_dimensions),G((iub-ilb+1)*(jub-jlb+1)*(kub-klb+1),n_obs))
        allocate(distance(n_obs),G((iub-ilb+1)*(jub-jlb+1)*(kub-klb+1),n_obs))

        select case (Me%ComputeOptions%KernelType)
        case (CharGaussian_)
            EE=2d0
        case (CharExponential_)
            EE=1d0
        case default
            stop 'Points2GridSmooth_3D - ModuleInterpolation - ERR01'
        end select

        do a=1,n_obs
            coord1(1)=Y(a)
            coord1(2)=X(a)
            coord1(3)=Z(a)
            do b=1,n_obs
                if(b.ne.a) then
                    coord2(1)=Y(b)
                    coord2(2)=X(b)
                    coord2(3)=Z(b)
                    do d=1,n_dimensions
                        if (Abs(coord1(d)-coord2(d)).lt.MinDist(a,d)) MinDist(a,d)=Abs(coord1(d)-coord2(d))
                    enddo
                endif
            enddo
            do d=1,n_dimensions
                tauR(a,d)=log(-1d0+2d0/phi)/(MinDist(a,d)**EE)
            enddo
        enddo
        
        do d=1,n_dimensions
            tauM(d)=sum(tauR(:,d))/real(n_obs)
        enddo
    
        ri=1
        do i=ilb,iub
        do j=jlb,jub
        do k=klb,kub
            DistSum=0d0
            do a=1,n_obs
                distance(a)=exp(-tauM(1)*(abs(CenterY(i,j)-Y(a))**EE) &
                                -tauM(2)*(abs(CenterX(i,j)-X(a))**EE) &
                                -tauM(3)*(abs(CenterZ(i,j,k)-Z(a))**EE))
                DistSum=DistSum+distance(a)
            enddo
            G2(ri,:)=0d0
            do a=1,n_obs
                G2(ri,a)=G2(ri,a)+distance(a)/DistSum
            enddo
            tau(ri,:)=0d0
            do a=1,n_obs
                tau(a,:)=tau(a,:)+G2(ri,a)*tauR(a,:)
            enddo
            ri=ri+1
        enddo
        enddo
        enddo

        ri=1
        do i=ilb,iub
        do j=jlb,jub
        do k=klb,kub
            DistSum=0d0
            do a=1,n_obs
 !falta aqui o teste ao openpoints3d do pai
                distance(a)=exp(-tau(ri,1)*(abs(CenterY(i,j)-Y(a))**EE)    &
                                -tau(ri,2)*(abs(CenterX(i,j)-X(a))**EE)    &
                                -tau(ri,3)*(abs(CenterZ(i,j,k)-Z(a))**EE)  &
                                )
                DistSum=DistSum+distance(a)
            enddo
            G(ri,:)=0d0
            do a=1,n_obs
                G(ri,a)=G(ri,a)+distance(a)/DistSum
            enddo
            v(i,j,k)=0d0
            do a=1,n_obs
                v(i,j,k)=v(i,j,k)+G(ri,a)*w(a)
            enddo
            ri=ri+1
        enddo
        enddo
        enddo

        nullify(x,y,z,w)


    end subroutine Points2GridSmooth_3D

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine NCInterpolation_3D
    
        select case (Me%ComputeOptions%NC_type)
        case(NC_User_)
            call NC_User_3D
        case(NC_Smoothing_)
            call NC_Smoothing_3D
        case(NC_Data_)
            call NC_Data_3D(.false.)
        case default
            stop 'NCInterpolation_3D - ModuleInterpolation - ERR01'
        end select

    end subroutine NCInterpolation_3D

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine NC_User_3D
        call FinalConvolution_3D (.false.)
    end subroutine NC_User_3D

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine NC_Smoothing_3D
        
        !Local-------------------------------------------------------------------
        integer                              :: n_dimensions
        integer                              :: i,j,k,fi,fj,fk
        integer                              :: n,nj,ri
        integer                              :: ilb, jlb, klb
        integer                              :: iub, jub, kub
        integer                              :: ilbF, jlbF, klbF
        integer                              :: iubF, jubF, kubF
        real(8)                                 :: EE,delta(3),DistSum,phi
        real,pointer,dimension(:)               :: tauM
        real,pointer,dimension(:,:)             :: tau,tauR
        real,pointer,dimension(:,:)             :: CenterXF
        real,pointer,dimension(:,:)             :: CenterYF
        real,pointer,dimension(:,:,:)           :: CenterZF
        real,pointer,dimension(:,:)             :: CenterX
        real,pointer,dimension(:,:)             :: CenterY
        real,pointer,dimension(:,:,:)           :: CenterZ
        real(8),allocatable,dimension(:)        :: distance
        real(8),allocatable,dimension(:,:)      :: G2

        !Begin shorten variables names-------------------------------------------
        n_dimensions    = Me%ComputeOptions%n_dimensions

        ilb             = Me%ExternalVar%SonGrid%WorkSize3D%ILB
        jlb             = Me%ExternalVar%SonGrid%WorkSize3D%JLB
        klb             = Me%ExternalVar%SonGrid%WorkSize3D%KLB
        iub             = Me%ExternalVar%SonGrid%WorkSize3D%IUB
        jub             = Me%ExternalVar%SonGrid%WorkSize3D%JUB
        kub             = Me%ExternalVar%SonGrid%WorkSize3D%KUB

        ilbF            = Me%ExternalVar%FatherGrid%WorkSize3D%ILB
        jlbF            = Me%ExternalVar%FatherGrid%WorkSize3D%JLB
        klbF            = Me%ExternalVar%FatherGrid%WorkSize3D%KLB
        iubF            = Me%ExternalVar%FatherGrid%WorkSize3D%IUB
        jubF            = Me%ExternalVar%FatherGrid%WorkSize3D%JUB
        kubF            = Me%ExternalVar%FatherGrid%WorkSize3D%KUB

        CenterXF        => Me%ExternalVar%FatherGrid%CenterX
        CenterYF        => Me%ExternalVar%FatherGrid%CenterY
        CenterZF        => Me%ExternalVar%FatherGrid%CenterZ
        CenterX         => Me%ExternalVar%SonGrid%CenterX
        CenterY         => Me%ExternalVar%SonGrid%CenterY
        CenterZ         => Me%ExternalVar%SonGrid%CenterZ

        tau             => Me%NC%tau
        tauM            => Me%NC%tauM
        tauR            => Me%NC%tauR
        phi             =  Me%ComputeOptions%phi
                
        !Begin ------------------------------------------------------------------

        if(phi.eq.1d0) then

            tau(:,:)=0d0

            else

            n=(iubF-ilbF+1)*(jubF-jlbF+1)*(kub-klb+1)

            select case (Me%ComputeOptions%KernelType)
            case (CharGaussian_)
                EE=2d0
            case (CharExponential_)
                EE=1d0
            case default
                stop 'GetSample_3D - ModuleInterpolation - ERR01'
            end select

            select case(n_dimensions)
            case(1)
                do fi=ilbF,iubF
                    if (fi .eq. ilbF ) then
                        delta(1)=abs(CenterYF(fi,1)-CenterYF(fi+1,1))
                    else if (fi .eq. iubF) then
                        delta(1)=abs(CenterYF(fi,1)-CenterYF(fi-1,1))
                    else
                        delta(1)=min(abs(CenterYF(fi,1)-CenterYF(fi+1,1)),  &
                                     abs(CenterYF(fi,1)-CenterYF(fi-1,1)))
                    endif
                    tauR(fi,1)=log(-1d0+2d0/phi)/(delta(1)**EE)
                enddo

                tauM(1)=sum(tauR(:,1))/real(iubF-ilbF+1)

                do i=ilb,iub
                    DistSum=0d0
                    do fi=iubF,ilbF
                        distance(fi)=exp(-tauM(1)*(abs(CenterY(i,1)-CenterYF(fi,1))**EE))
                        DistSum=DistSum+distance(fi)
                    enddo
                    G2(i,:)=0d0
                    do fi=iubF,ilbF
                        G2(i,fi)=G2(i,fi)+distance(fi)/DistSum
                    enddo
                    tau(i,1)=0d0
                    do fi=iubF,ilbF
                        tau(i,1)=tau(i,1)+G2(i,fi)*tauR(fi,1)
                    enddo
                enddo
                call FinalConvolution_3D (.false.)

            case(2)
            
                nj=1
                do fi=ilbF,iub
                do fj=jlbF,jub
                    if(fi.eq.ilbF) then
                        delta(1)=abs(CenterYF(fi,fj)-CenterYF(fi+1,fj))
                    elseif(fi .eq. iubF) then
                        delta(1)=abs(CenterYF(fi,fj)-CenterYF(fi-1,fj))
                    else
                        delta(1)=min(abs(CenterYF(fi,fj)-CenterYF(fi+1,fj)),  &
                                     abs(CenterYF(fi,fj)-CenterYF(fi-1,fj)))
                    endif
                        
                    if(fj .eq. jlbF) then
                        delta(2)=abs(CenterXF(fi,fj)-CenterXF(fi,fj+1))
                    elseif(fj .eq. jubF) then
                        delta(2)=abs(CenterXF(fi,fj)-CenterXF(fi,fj-1))
                    else
                        delta(2)=min(abs(CenterXF(fi,fj)-CenterXF(fi,fj+1)),  &
                                     abs(CenterXF(fi,fj)-CenterXF(fi,fj-1)))
                    endif

                    tauR(nj,1)=log(-1d0+2d0/phi)/delta(1)
                    tauR(nj,2)=log(-1d0+2d0/phi)/delta(2)

                    nj=nj+1
                enddo
                enddo

                tauM(1)=sum(tauR(:,1))/real(n)
                tauM(2)=sum(tauR(:,2))/real(n)
            
                ri=1
                do i=ilb,iub
                do j=jlb,jub
                    nj=1
                    DistSum=0d0
                    do fi=iubF,ilbF
                    do fj=jubF,jlbF
                        distance(nj)=exp(-tauM(1)*(abs(CenterY(i,j)-CenterYF(fi,fj))**EE) &
                                         -tauM(2)*(abs(CenterX(i,j)-CenterXF(fi,fj))**EE))
                        DistSum=DistSum+distance(nj)
                        nj=nj+1
                    enddo
                    enddo
                    G2(ri,:)=0d0
                    do nj=1,n
                        G2(ri,nj)=G2(ri,nj)+distance(nj)/DistSum
                    enddo
                    tau(ri,:)=0d0
                    do nj=1,n
                        tau(nj,:)=tau(nj,:)+G2(ri,nj)*tauR(nj,:)
                    enddo
                    ri=ri+1
                enddo
                enddo
                call FinalConvolution_3D (.false.)

            case(3)
            
                nj=1
                do fi=ilbF,iub
                do fj=jlbF,jub
                    if(fi.eq.ilbF) then
                        delta(1)=abs(CenterYF(fi,fj)-CenterYF(fi+1,fj))
                    elseif(fi.eq.iubF) then
                        delta(1)=abs(CenterYF(fi,fj)-CenterYF(fi-1,fj))
                    else
                        delta(1)=min(abs(CenterYF(fi,fj)-CenterYF(fi+1,fj)),  &
                                     abs(CenterYF(fi,fj)-CenterYF(fi-1,fj)))
                    endif
                        
                    if(fj.eq.jlbF) then
                        delta(2)=abs(CenterXF(fi,fj)-CenterXF(fi,fj+1))
                    elseif(fj.eq.jubF) then
                        delta(2)=abs(CenterXF(fi,fj)-CenterXF(fi,fj-1))
                    else
                        delta(2)=min(abs(CenterXF(fi,fj)-CenterXF(fi,fj+1)),  &
                                     abs(CenterXF(fi,fj)-CenterXF(fi,fj-1)))
                    endif

                    tauR(nj,1)=log(-1d0+2d0/phi)/delta(1)
                    tauR(nj,2)=log(-1d0+2d0/phi)/delta(2)

                    do fk=klbF,kubF
                        if(fk.eq.klbF) then
                            delta(3)=abs(CenterZF(fi,fj,fk)-CenterZF(fi,fj,fk+1))
                        elseif(fk.eq.kubF) then
                            delta(3)=abs(CenterZF(fi,fj,fk)-CenterZF(fi,fj,fk-1))
                        else
                            delta(3)=min(abs(CenterZF(fi,fj,fk)-CenterZF(fi,fj,fk+1)), &
                                         abs(CenterZF(fi,fj,fk)-CenterZF(fi,fj,fk-1)))
                        endif
                        tauR(nj,3)=log(-1d0+2d0/phi)/delta(3)
                        nj=nj+1
                    enddo
                enddo
                enddo

                tauM(1)=sum(tauR(:,1))/real(n)
                tauM(2)=sum(tauR(:,2))/real(n)
                tauM(3)=sum(tauR(:,3))/real(n)
            
                ri=1
                do i=ilb,iub
                do j=jlb,jub
                do k=klb,kub
                    nj=1
                    DistSum=0d0
                    do fi=iubF,ilbF
                    do fj=jubF,jlbF
                    do fk=kubF,klbF
                        distance(nj)=exp(-tauM(1)*(abs(CenterY(i,j)-CenterYF(fi,fj))**EE) &
                                         -tauM(2)*(abs(CenterX(i,j)-CenterXF(fi,fj))**EE) &
                                         -tauM(3)*(abs(CenterZ(i,j,k)-CenterZF(fi,fj,fk))**EE))
                        DistSum=DistSum+distance(nj)
                        nj=nj+1
                    enddo
                    enddo
                    enddo
                    G2(ri,:)=0d0
                    do nj=1,n
                        G2(ri,nj)=G2(ri,nj)+distance(nj)/DistSum
                    enddo
                    tau(ri,:)=0d0
                    do nj=1,n
                        tau(nj,:)=tau(nj,:)+G2(ri,nj)*tauR(nj,:)
                    enddo
                    ri=ri+1
                enddo
                enddo
                enddo
                call FinalConvolution_3D (.false.)
            end select
        endif

    end subroutine NC_Smoothing_3D


    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine NC_Data_3D (FatherGrid)
        
        !Arguments---------------------------------------------------------------
        logical                                 :: FatherGrid

        !Local-------------------------------------------------------------------
        integer                              :: n_dimensions, n_groups
        integer                              :: h,i,j,k
        integer                              :: ilb, jlb, klb
        integer                              :: iub, jub, kub

        !Begin shorten variables names-------------------------------------------
        n_dimensions    = Me%ComputeOptions%n_dimensions
        n_groups        = Me%ComputeOptions%n_groups

        ilb             = Me%ExternalVar%SonGrid%WorkSize3D%ILB
        jlb             = Me%ExternalVar%SonGrid%WorkSize3D%JLB
        klb             = Me%ExternalVar%SonGrid%WorkSize3D%KLB
        iub             = Me%ExternalVar%SonGrid%WorkSize3D%IUB
        jub             = Me%ExternalVar%SonGrid%WorkSize3D%JUB
        kub             = Me%ExternalVar%SonGrid%WorkSize3D%KUB
        
        !Begin ------------------------------------------------------------------

        do h=1, n_groups**n_dimensions
            k=klb+int(real(h-1)/real(n_groups**2))
            j=jlb+int(real(h-(k-klb)*(n_groups**2)-1)/real(n_groups))
            i=ilb+h-(k-klb)*(n_groups**2)-(j-jlb)*n_groups-1
            call GetSample_3D(h,i,j,k)
            call LogRegr_3D(h)
        enddo
        call MasterG_3D
        call FinalConvolution_3D (FatherGrid)

    end subroutine NC_Data_3D

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine GetSample_3D(h,i,j,k)

        !Arguments---------------------------------------------------------------
        integer                              :: h,i,j,k
        
        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                              :: n_dimensions
        integer                              :: n_iterations, n_groups
        integer                              :: d
        integer                              :: sample_size
        integer                              :: ilb, jlb, klb
        integer                              :: iub, jub, kub
        integer                              :: lb(3),ub(3)
        integer                              :: point1(3), point2(3), point3(3)
        integer,allocatable,dimension(:)     :: sample
        real(8)                                 :: EE, rnd(3)
        real(8)                                 :: w1,w2,w3
        real(8)                                 :: d21,d23
        logical                                 :: point3_notfound
        real,pointer,dimension(:)               :: GroupCenterX, GroupCenterY
        real,pointer,dimension(:)               :: GroupCenterZ
        integer,pointer,dimension(:,:,:)     :: OpenPoints3D
        real,pointer,dimension(:,:)             :: CenterX, CenterY
        real,pointer,dimension(:,:,:)           :: CenterZ
        real,pointer,dimension(:,:,:)           :: predictor,response,w
        type (T_Grid),  pointer                 :: FatherSource
        type (T_Field), pointer                 :: FatherField
        
        !Begin shorten variables names-------------------------------------------
        n_dimensions    = Me%ComputeOptions%n_dimensions
        n_groups        = Me%ComputeOptions%n_groups
        sample_size     = Me%ComputeOptions%sample_size
        allocate(sample(n_dimensions))

        FatherSource => Me%ExternalVar%FatherGrid
        FatherField  => Me%ExternalVar%FatherField

        ilb             = FatherSource%WorkSize3D%ILB
        jlb             = FatherSource%WorkSize3D%JLB
        klb             = FatherSource%WorkSize3D%KLB
        iub             = FatherSource%WorkSize3D%IUB
        jub             = FatherSource%WorkSize3D%JUB
        kub             = FatherSource%WorkSize3D%KUB
        CenterX         => FatherSource%CenterX
        CenterY         => FatherSource%CenterY
        CenterZ         => FatherSource%CenterZ
        OpenPoints3D    => FatherSource%OpenPoints3D
        w               => FatherField%Values3D

        GroupCenterX    => Me%NC%GroupCenterX
        GroupCenterY    => Me%NC%GroupCenterY
        GroupCenterZ    => Me%NC%GroupCenterZ
        response        => Me%NC%response
        predictor       => Me%NC%predictor

        call random_seed
    
        select case (Me%ComputeOptions%KernelType)
        case (CharGaussian_)
            EE=2d0
        case (CharExponential_)
            EE=1d0
        case default
            stop 'GetSample_3D - ModuleInterpolation - ERR01'
        end select

        lb(1)=ilb+(i-1)*int(real(iub-ilb+1)/real(n_groups))
        lb(2)=jlb+(j-1)*int(real(jub-jlb+1)/real(n_groups))
        lb(3)=klb+(k-1)*int(real(kub-klb+1)/real(n_groups))
        ub(1)=iub-(n_groups-i)*int(real(iub-ilb+1)/real(n_groups))
        ub(2)=jub-(n_groups-j)*int(real(jub-jlb+1)/real(n_groups))
        ub(3)=kub-(n_groups-k)*int(real(kub-klb+1)/real(n_groups))

        sample(:)=0
        GroupCenterX(h)=0d0
        GroupCenterY(h)=0d0
        GroupCenterZ(h)=0d0
    
        do d=1, n_dimensions
            do while(sample(d).lt.sample_size)
                call random_number (rnd(:))
                point1(:)=lb(:)+int(rnd(:)*real(ub(:)-lb(:)+1))
                point2(:)=point1(:)
                point3(:)=point1(:)
        
                point3_notfound=.true.
                n_iterations=0
                do while(point3_notfound.and.n_iterations.lt.20)
                    call random_number (rnd(1))
                    point3(d)=lb(d)+int(rnd(1)*real(ub(d)-lb(d)+1))
                    if(abs(point3(d)-point1(d)).ge.2) point3_notfound=.false.
                    n_iterations=n_iterations+1
                enddo
                if(n_iterations.eq.20) stop '*********** - ModuleInterpolation - ERR01'

                call random_number (rnd(1))
                point2(d)=1+min(point1(d),point3(d))+int(rnd(1)*real(abs(point3(d)-point1(d)-1)))
    
                w1=w(point1(1),point1(2),point1(3))
                w2=w(point2(1),point2(2),point2(3))
                w3=w(point3(1),point3(2),point3(3))

                select case(d)
                case(1)
                    d21=abs(CenterY(point2(1),point2(2))-CenterY(point1(1),point1(2)))**EE
                    d23=abs(CenterY(point2(1),point2(2))-CenterY(point3(1),point3(2)))**EE
                case(2)
                    d21=abs(CenterX(point2(1),point2(2))-CenterX(point1(1),point1(2)))**EE
                    d23=abs(CenterX(point2(1),point2(2))-CenterX(point3(1),point3(2)))**EE
                case(3)
                    d21=abs(CenterZ(point2(1),point2(2),point2(3))-CenterZ(point1(1),point1(2),point1(3)))**EE
                    d23=abs(CenterZ(point2(1),point2(2),point2(3))-CenterZ(point3(1),point3(2),point3(3)))**EE
                end select

                if ((w1.le.w2.and.w2.le.w3).or.(w1.ge.w2.and.w2.ge.w3)) then
                    if(w1.lt.w3) then
                        response(h,d,sample(d))=(w2-w1)/(w3-w1)
                        predictor(h,d,sample(d))=d23-d21
                        GroupCenterY(h)=GroupCenterY(h)+CenterY(point1(1),point1(2))/real(sample_size)
                        GroupCenterX(h)=GroupCenterX(h)+CenterX(point1(1),point1(2))/real(sample_size)
                        GroupCenterZ(h)=GroupCenterZ(h)+CenterZ(point1(1),point1(2),point1(3))/real(sample_size)
                        else
                        response(h,d,sample(d))=(w2-w3)/(w1-w3)
                        predictor(h,d,sample(d))=d21-d23
                        GroupCenterY(h)=GroupCenterY(h)+CenterY(point3(1),point3(2))/real(sample_size)
                        GroupCenterX(h)=GroupCenterX(h)+CenterX(point3(1),point3(2))/real(sample_size)
                        GroupCenterZ(h)=GroupCenterZ(h)+CenterZ(point3(1),point3(2),point1(3))/real(sample_size)
                    endif
                    sample(d)=sample(d)+1
                endif
            enddo
        enddo
        deallocate(sample)
    end subroutine GetSample_3D

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine LogRegr_3D(h)

        !Arguments---------------------------------------------------------------
        integer                              :: h
        
        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                 :: n_dimensions, n_groups, n_iterations
        integer                                 :: a,d, sample_size
        integer                                 :: max_iterations
        real(8)                                 :: aux,eps,pi,summation(2)
        real,pointer,dimension(:,:)             :: tauR
        real,pointer,dimension(:,:,:)           :: predictor,response

        !Begin shorten variables names-------------------------------------------
        
        n_dimensions    = Me%ComputeOptions%n_dimensions
        n_groups        = Me%ComputeOptions%n_groups
        max_iterations  = Me%ComputeOptions%max_iterations
        sample_size     = Me%ComputeOptions%sample_size

        predictor       => Me%NC%predictor
        response        => Me%NC%response
        tauR            => Me%NC%tauR

        !Begin ------------------------------------------------------------------

        tauR(:,:)=0d0
    do d=1,n_dimensions
            n_iterations=0
            do while (eps.gt.0.000001d0.and.n_iterations.lt.max_iterations)
                summation(1) = 0d0
                summation(2) = 0d0
                do a = 1,sample_size
                    pi = exp(tauR(h,d) * predictor(h,d,a)) / (1 + exp(tauR(h,d) * predictor(h,d,a)))
                    summation(1) = summation(1) + predictor(h,d,a) * (response(h,d,a) - pi)
                    summation(2) = summation(2) + (predictor(h,d,a)**2d0) * pi * (1d0 - pi)
                enddo
                aux = tauR(h,d) + summation(1) / summation(2)
                eps = abs(1d0 - tauR(h,d) / aux)
                tauR(h,d) = aux
                n_iterations = n_iterations + 1
            enddo
            if(n_iterations.eq.max_iterations) stop '*********** - ModuleInterpolation - ERR01'
        enddo

    end subroutine LogRegr_3D

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine MasterG_3D
 
        !Arguments---------------------------------------------------------------
        
        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                              :: d,h,i,j,k,r,ri
        integer                              :: n_dimensions, n_groups
        integer                              :: ilb, jlb, klb
        integer                              :: iub, jub, kub
        real(8)                                 :: DistSum, EE
        real(8),allocatable,dimension(:)        :: distance
        real,pointer,dimension(:)               :: GroupCenterX
        real,pointer,dimension(:)               :: GroupCenterY
        real,pointer,dimension(:)               :: GroupCenterZ
        real,allocatable,dimension(:,:)         :: G2
        real,pointer,dimension(:)               :: tauM
        real,pointer,dimension(:,:)             :: tau,tauR
        real,pointer,dimension(:,:)             :: CenterX
        real,pointer,dimension(:,:)             :: CenterY
        real,pointer,dimension(:,:,:)           :: CenterZ

        !Begin shorten variables names-------------------------------------------
        
        n_dimensions    = Me%ComputeOptions%n_dimensions
        n_groups        = Me%ComputeOptions%n_groups

        ilb             = Me%ExternalVar%SonGrid%WorkSize3D%ILB
        jlb             = Me%ExternalVar%SonGrid%WorkSize3D%JLB
        klb             = Me%ExternalVar%SonGrid%WorkSize3D%KLB
        iub             = Me%ExternalVar%SonGrid%WorkSize3D%IUB
        jub             = Me%ExternalVar%SonGrid%WorkSize3D%JUB
        kub             = Me%ExternalVar%SonGrid%WorkSize3D%KUB

        tau             => Me%NC%tau
        tauM            => Me%NC%tauM
        tauR            => Me%NC%tauR

        CenterX         => Me%ExternalVar%SonGrid%CenterX
        CenterY         => Me%ExternalVar%SonGrid%CenterY
        CenterZ         => Me%ExternalVar%SonGrid%CenterZ

        GroupCenterX    => Me%NC%GroupCenterX
        GroupCenterY    => Me%NC%GroupCenterY
        GroupCenterZ    => Me%NC%GroupCenterZ

        !Begin ------------------------------------------------------------------

        r=(iub-ilb+1)*(jub-jlb+1)*(kub-klb+1)

        allocate(G2(r,n_dimensions))
        allocate(distance(n_groups**n_dimensions))

        select case (Me%ComputeOptions%KernelType)
        case (CharGaussian_)
            EE=2d0
        case (CharExponential_)
            EE=1d0
        case default
            stop 'GetSample_3D - ModuleInterpolation - ERR01'
        end select

        !master convolution kernel's shape parameters
        do d=1,n_dimensions
            tauM(d)=sum(tauR(:,d))/real(n_groups**n_dimensions)
        enddo

        ri=1
        do i=ilb,iub
        do j=jlb,jub
        do k=klb,kub
            DistSum=0d0
            do h=1,n_groups**n_dimensions
                distance(h)=exp(-tauM(1)*(abs(CenterY(i,j)-GroupCenterY(h))**EE)     &
                                -tauM(2)*(abs(CenterX(i,j)-GroupCenterX(h))**EE)     &
                                -tauM(3)*(abs(CenterZ(i,j,k)-GroupCenterZ(h))**EE)   &
                                )
                DistSum=DistSum+distance(h)
            enddo
            G2(ri,:)=0d0
            do h=1,n_groups**n_dimensions
                G2(ri,h)=G2(ri,h)+distance(h)/DistSum
            enddo
            tau(ri,:)=0d0
            do h=1,n_groups**n_dimensions
                tau(ri,:)=tau(ri,:)+G2(ri,h)*tauR(h,:)
            enddo
            ri=ri+1
        enddo
        enddo
        enddo

        deallocate(distance,G2)

    end subroutine MasterG_3D

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine FinalConvolution_3D (FatherGrid)

        !Arguments---------------------------------------------------------------
        logical                                 :: FatherGrid

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                              :: a,b,c,i,j,k,n,nj,r,ri
        integer                              :: ilb, jlb, klb
        integer                              :: iub, jub, kub
        integer                              :: ilbF, jlbF, klbF
        integer                              :: iubF, jubF, kubF
        real(8)                                 :: DistSum,EE
        real(8),allocatable,dimension(:)        :: distance
        real(8),allocatable,dimension(:,:)      :: G
        real,pointer,dimension(:,:)             :: tau
        real,pointer,dimension(:,:)             :: CenterXF
        real,pointer,dimension(:,:)             :: CenterYF
        real,pointer,dimension(:,:,:)           :: CenterZF
        real,pointer,dimension(:,:)             :: CenterX
        real,pointer,dimension(:,:)             :: CenterY
        real,pointer,dimension(:,:,:)           :: CenterZ
        real,pointer,dimension(:,:,:)           :: v,w
        type (T_Grid),  pointer                 :: FatherSource
        type (T_Field), pointer                 :: FatherField

        !Begin shorten variables names-------------------------------------------
        
        if(FatherGrid) then
            FatherSource => Me%ExternalVar%FatherGrid
            FatherField  => Me%ExternalVar%FatherField
        else
            FatherSource => Me%ExternalVar%FatherGrid
            FatherField  => Me%ExternalVar%FatherField
        endif

        ilbF            = FatherSource%WorkSize3D%ILB
        jlbF            = FatherSource%WorkSize3D%JLB
        klbF            = FatherSource%WorkSize3D%KLB
        iubF            = FatherSource%WorkSize3D%IUB
        jubF            = FatherSource%WorkSize3D%JUB
        kubF            = FatherSource%WorkSize3D%KUB

        w               => FatherField%Values3D

        CenterXF        => FatherSource%CenterX
        CenterYF        => FatherSource%CenterY
        CenterZF        => FatherSource%CenterZ
        
        tau             => Me%NC%tau
        
        ilb             = Me%ExternalVar%SonGrid%WorkSize3D%ILB
        jlb             = Me%ExternalVar%SonGrid%WorkSize3D%JLB
        klb             = Me%ExternalVar%SonGrid%WorkSize3D%KLB
        iub             = Me%ExternalVar%SonGrid%WorkSize3D%IUB
        jub             = Me%ExternalVar%SonGrid%WorkSize3D%JUB
        kub             = Me%ExternalVar%SonGrid%WorkSize3D%KUB

        v               => Me%ExternalVar%SonField%Values3D

        CenterX         => Me%ExternalVar%SonGrid%CenterX
        CenterY         => Me%ExternalVar%SonGrid%CenterY
        CenterZ         => Me%ExternalVar%SonGrid%CenterZ

        !Begin-------------------------------------------------------------------
        
        n=(iubF-ilbF+1)*(jubF-jlbF+1)*(kubF-klbF+1)
        r=(iub-ilb+1)*(jub-jlb+1)*(kub-klb+1)
        select case (Me%ComputeOptions%KernelType)
        case (CharGaussian_)
            EE=2d0
        case (CharExponential_)
            EE=1d0
        case default
            stop 'GetSample_3D - ModuleInterpolation - ERR01'
        end select

        allocate(distance(n),G(r,n))

        ri=1
        do i=ilb,iub
        do j=jlb,jub
        do k=klb,kub
            DistSum=0d0
            nj=1
            do a=ilbF,iubF
            do b=jlbF,jubF
            do c=klbF,kubF
 !falta aqui o teste ao openpoints3d do pai
                distance(nj)=exp(-tau(ri,1)*(abs(CenterY(i,j)-CenterYF(a,b))**EE)     &
                                 -tau(ri,2)*(abs(CenterX(i,j)-CenterXF(a,b))**EE)     &
                                 -tau(ri,3)*(abs(CenterZ(i,j,k)-CenterZF(a,b,c))**EE) &
                                 )
                DistSum=DistSum+distance(nj)
                nj=nj+1
            enddo
            enddo
            enddo
            G(ri,:)=0d0
            do nj=1,n
                G(ri,nj)=G(ri,nj)+distance(nj)/DistSum
            enddo
            v(i,j,k)=0d0
            nj=1
            do a=ilbF,iubF
            do b=jlbF,jubF
            do c=klbF,kubF
                v(i,j,k)=v(i,j,k)+G(ri,nj)*w(a,b,c)
                nj=nj+1
            enddo
            enddo
            enddo
            ri=ri+1
        enddo
        enddo
        enddo

        deallocate(distance, G)

    end subroutine FinalConvolution_3D

    !----------------------------------------------------------------------------

        

    !----------------------------------------------------------------------------
        
!    subroutine ConservativeModifier
!        call CellSelector
!        call FatherToSonGrid
!        call Integrator
!    end subroutine

    !----------------------------------------------------------------------------
    !subroutine ConservativeModifier_2D(p_2D)
        
        !Arguments----------------------------------------------------------
    !    real, pointer, dimension(:,:)               :: p_2D

!        call CellSelector_2D
!        call FatherToSonGrid_2D
!        call Integrator_2D

    !end subroutine

    !----------------------------------------------------------------------------
        
    subroutine ConservativeModifier_3D()

        !Arguments----------------------------------------------------------
        !real, pointer, dimension(:,:,:)               :: p_3D

        call CellSelector_3D
        call FatherToSonGrid_3D
        call Integrator_3D

    end subroutine

    !----------------------------------------------------------------------------

    !subroutine NonConservativeModifier_3D(p_3D)


        !Arguments----------------------------------------------------------
    !    real, pointer, dimension(:,:,:)               :: p_3D


        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------


        !call interface


    !end subroutine NonConservativeModifier_3D


    !----------------------------------------------------------------------------

    !subroutine NonConservativeModifier_2D(p_2D)


        !Arguments----------------------------------------------------------
    !    real, pointer, dimension(:,:)               :: p_2D


        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------


        !call interface


    !end subroutine NonConservativeModifier_2D


    !----------------------------------------------------------------------------

    
    subroutine CellSelector_3D
 
        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                              :: n_runs
        integer                              :: iub
        integer                              :: jub
        integer                              :: kub
        integer                              :: h,i,j,k
        integer, pointer, dimension(:,:)     :: WaterPoints2D
        integer, pointer, dimension(:,:,:)   :: OpenPoints3D
        integer, pointer, dimension(:,:,:)   :: MatrixUpdate3D

        !Begin shorten variables names-------------------------------------------
        iub             = Me%ExternalVar%FatherGrid%WorkSize3D%IUB
        jub             = Me%ExternalVar%FatherGrid%WorkSize3D%JUB
        kub             = Me%ExternalVar%FatherGrid%WorkSize3D%KUB
        WaterPoints2D   => Me%ExternalVar%FatherGrid%WaterPoints2D
        OpenPoints3D    => Me%ExternalVar%FatherGrid%OpenPoints3D
        MatrixUpdate3D  => Me%MatrixUpdate3D       

        !Begin-------------------------------------------------------------------
        
        ! for a given iteration, selects a cell and sees if its P should be updated; 
        
        n_runs = 4
        do h=1,n_runs
            do i=1,iub
                do j=1,jub
                    if(WaterPoints2D(i,j).eq.1) then
                        do k=1,kub
                            if(OpenPoints3D(i,j,k).eq.1) then
                                if (MatrixUpdate3D(i,j,k).eq.mod(h,2)) then
                                    call UpdateP_3D(i,j,k)
                                endif
                            endif
                        enddo
                    endif
                enddo
            enddo
        enddo
            

    end subroutine CellSelector_3D

    !----------------------------------------------------------------------------

    subroutine UpdateP_3D(i,j,k)
 
        !Arguments---------------------------------------------------------------
        integer                          :: i,j,k

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                          :: n_dimensions
        real, pointer, dimension(:,:)       :: CX,CY,DX,DY
        real, pointer, dimension(:,:,:)     :: CZ,DZ
        real, pointer, dimension(:,:,:)     :: OpenPoints3D, p, w

        !Begin shorten variables names-------------------------------------------

        n_dimensions = Me%ComputeOptions%n_dimensions
        
        w           => Me%ExternalVar%FatherField%Values3D
        p           => Me%ExternalVar%FatherGrid%p
        CX          => Me%ExternalVar%FatherGrid%CX
        CY          => Me%ExternalVar%FatherGrid%CY
        CZ          => Me%ExternalVar%FatherGrid%CZ
        DX          => Me%ExternalVar%FatherGrid%DX
        DY          => Me%ExternalVar%FatherGrid%DY
        DZ          => Me%ExternalVar%FatherGrid%DZ


        !Begin-------------------------------------------------------------------
        
        select case (n_dimensions)

        case (1)
            if(OpenPoints3D(i-1,j,k).eq.0) p(i-1,j,k)=p(i,j,k)
            if(OpenPoints3D(i+1,j,k).eq.0) p(i+1,j,k)=p(i,j,k)

            p(i,j,k)=(CY(i-1,j)*DY(i,j)*p(i+1,j,k)+CY(i,j)*(DY(i,j)         &
                *p(i-1,j,k)-8d0*CY(i-1,j)*w(i,j,k)))/(CY(i,j)*DY(i,j)       &
                +CY(i-1,j)*(-8d0*CY(i,j)+DY(i,j))) !Quinta, 23 Junho
            
        case(2)
            if(OpenPoints3D(i-1,j,k).eq.0) p(i-1,j,k)=p(i,j,k)
            if(OpenPoints3D(i,j-1,k).eq.0) p(i,j-1,k)=p(i,j,k)
            if(OpenPoints3D(i+1,j,k).eq.0) p(i+1,j,k)=p(i,j,k)
            if(OpenPoints3D(i,j+1,k).eq.0) p(i,j+1,k)=p(i,j,k)

            p(i,j,k)=(CX(i,j)*CY(i-1,j)*CY(i,j)*DX(i,j)*p(i,j-1,k)          &
                +CX(i,j-1)*(CY(i-1,j)*CY(i,j)*DX(i,j)*p(i,1+j,k)+CX(i,j)    &
                *(CY(i-1,j)*DY(i,j)*p(1+i,j,k)+CY(i,j)*(DY(i,j)*p(i-1,j,k)  &
                -12d0*CY(i-1,j)*w(i,j,k)))))/(CX(i,j)*CY(i-1,j)*CY(i,j)     &
                *DX(i,j)+CX(i,j-1)*(CY(i-1,j)*CY(i,j)*DX(i,j)+CX(i,j)       &
                *(CY(i,j)*DY(i,j)+CY(i-1,j)*(-12d0*CY(i,j)+DY(i,j)))))
            !Quarta 22/06

        case(3)
            if(OpenPoints3D(i-1,j,k).eq.0) p(i-1,j,k)=p(i,j,k)
            if(OpenPoints3D(i,j-1,k).eq.0) p(i,j-1,k)=p(i,j,k)
            if(OpenPoints3D(i,j,k-1).eq.0) p(i,j,k-1)=p(i,j,k)
            if(OpenPoints3D(i+1,j,k).eq.0) p(i+1,j,k)=p(i,j,k)
            if(OpenPoints3D(i,j+1,k).eq.0) p(i,j+1,k)=p(i,j,k)
            if(OpenPoints3D(i,j,k+1).eq.0) p(i,j,k+1)=p(i,j,k)

            p(i,j,k)=(CX(i,j)*CY(i-1,j)*CY(i,j)*CZ(i,j,k-1)*CZ(i,j,k)           &
                *DX(i,j)*p(i,j-1,k)+CX(i,j-1)*(CY(i-1,j)*CY(i,j)                &
                *CZ(i,j,k-1)*CZ(i,j,k)*DX(i,j)*p(i,j+1,k)+CX(i,j)               &
                *(CY(i-1,j)*CZ(i,j,k-1)*CZ(i,j,k)*DY(i,j)*p(i+1,j,k)            &
                +CY(i,j)*(CY(i-1,j)*CZ(i,j,k)*DZ(i,j,k)*p(i,j,k-1)              &
                +CZ(i,j,k-1)*(CY(i-1,j)*DZ(i,j,k)*p(i,j,k+1)+CZ(i,j,k)          &
                *(DY(i,j)*p(i-1,j,k)-16d0*CY(i-1,j)*w(i,j,k)))))))              &
                /(CX(i,j)*CY(i-1,j)*CY(i,j)*CZ(i,j,k-1)*CZ(i,j,k)*DX(i,j)       &
                +CX(i,j-1)*(CY(i-1,j)*CY(i,j)*CZ(i,j,k-1)*CZ(i,j,k)*DX(i,j)     &
                +CX(i,j)*(CY(i,j)*CZ(i,j,k-1)*CZ(i,j,k)*DY(i,j)+CY(i-1,j)       &
                *(CZ(i,j,k-1)*CZ(i,j,k)*DY(i,j)+CY(i,j)*(CZ(i,j,k)              &
                *DZ(i,j,k)+CZ(i,j,k-1)*(-16d0*CZ(i,j,k)+DZ(i,j,k))))))) 
                !Quarta 22/06

        case default
            stop 'UpdateP_3D - ModuleInterpolation - ERR01' 
        end select

    end subroutine UpdateP_3D
    
    !----------------------------------------------------------------------------

    subroutine FatherToSonGrid_3D
 
        !Arguments---------------------------------------------------------------
        
        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                              :: n_dimensions

        integer                              :: a,b,c
        integer                              :: i,j,k
        integer                              :: ilb, jlb, klb
        integer                              :: iub, jub, kub

        integer, pointer, dimension(:,:)     :: WaterPoints2D
        integer, pointer, dimension(:,:,:)   :: OpenPoints3D
        real, pointer, dimension(:,:,:)         :: pSon

        integer, pointer, dimension(:,:)     :: FatherSonX, FatherSonY
        integer, pointer, dimension(:,:,:)   :: FatherSonZ
        
        character(3), pointer, dimension(:,:,:) :: InsideFatherCell

        real, pointer, dimension(:,:)           :: CenterX, CenterY
        real, pointer, dimension(:,:,:)         :: CenterZ
        real, pointer, dimension(:,:)           :: X, Y
        real, pointer, dimension(:,:,:)         :: Z
        real, pointer, dimension(:,:)           :: CX, CY
        real, pointer, dimension(:,:,:)         :: CZ
        real, pointer, dimension(:,:,:)         :: p

        !Begin shorten variables names-------------------------------------------
        
        n_dimensions    = Me%ComputeOptions%n_dimensions

        ilb             = Me%ExternalVar%SonGrid%WorkSize3D%ILB
        jlb             = Me%ExternalVar%SonGrid%WorkSize3D%JLB
        klb             = Me%ExternalVar%SonGrid%WorkSize3D%KLB
        iub             = Me%ExternalVar%SonGrid%WorkSize3D%IUB
        jub             = Me%ExternalVar%SonGrid%WorkSize3D%JUB
        kub             = Me%ExternalVar%SonGrid%WorkSize3D%KUB

        X               => Me%ExternalVar%SonGrid%CenterX
        Y               => Me%ExternalVar%SonGrid%CenterY
        Z               => Me%ExternalVar%SonGrid%CenterZ

        CenterX         => Me%ExternalVar%FatherGrid%CenterX
        CenterY         => Me%ExternalVar%FatherGrid%CenterY
        CenterZ         => Me%ExternalVar%FatherGrid%CenterZ
        p               => Me%ExternalVar%FatherGrid%p
        pSon            => Me%ExternalVar%SonGrid%p

        FatherSonX       => Me%FatherSonX
        FatherSonY       => Me%FatherSonY
        FatherSonZ       => Me%FatherSonZ
        InsideFatherCell => Me%InsideFatherCell
        
        CX              => Me%ExternalVar%FatherGrid%CX
        CY              => Me%ExternalVar%FatherGrid%CY
        CZ              => Me%ExternalVar%FatherGrid%CZ

        WaterPoints2D   => Me%ExternalVar%SonGrid%WaterPoints2D
        OpenPoints3D    => Me%ExternalVar%SonGrid%OpenPoints3D
        

        !Begin-------------------------------------------------------------------

        do a=ilb,iub
        do b=jlb,jub
        if(WaterPoints2D(a,b).eq.1) then
            i=FatherSonY(a,b)
            j=FatherSonX(a,b)
            do c=klb,kub
                if(OpenPoints3D(a,b,c).eq.1) then
                    k=FatherSonZ(a,b,c)
                    select case(n_dimensions)
                    case(1)
                        select case(InsideFatherCell(a,b,c))
                        case("Dl1")
! falta pôr as protecções contra celulas adjacentes vazias
                            pSon(a,b,c)=p(i,j,k)+(CenterY(i,j)-Y(i,j))        &
                                *(p(i-1,j,k)-p(i,j,k))/CY(i-1,j)
                        case("Du1")
                            pSon(a,b,c)=p(i,j,k)+(Y(i,j)-CenterY(i,j))        &
                                *(p(i+1,j,k)-p(i,j,k))/CY(i,j)
                        case default
                            stop 'FatherToSonGrid_3D - ModuleInterpolation - ERR01' 
                        end select
                    case(2)
                        select case(InsideFatherCell(a,b,c))
                        case("Al1")
                            pSon(a,b,c)=p(i,j,k)+(CenterY(i,j)-Y(i,j))        &
                                *(p(i-1,j,k)-p(i,j,k))/CY(i-1,j)
                        case("Al2")
                            pSon(a,b,c)=p(i,j,k)+(CenterX(i,j)-X(i,j))        &
                                *(p(i,j-1,k)-p(i,j,k))/CX(i,j-1)
                        case("Au1")
                            pSon(a,b,c)=p(i,j,k)+(Y(i,j)-CenterY(i,j))        &
                                *(p(i+1,j,k)-p(i,j,k))/CY(i,j)
                        case("Au2")
                            pSon(a,b,c)=p(i,j,k)+(X(i,j)-CenterX(i,j))        &
                                *(p(i,j+1,k)-p(i,j,k))/CX(i,j)
                        case default
                            stop 'FatherToSonGrid_3D - ModuleInterpolation - ERR02'
                        end select
                    case(3)
                        select case(InsideFatherCell(a,b,c))
                        case("Vl1")
                            pSon(a,b,c)=p(i,j,k)+(CenterY(i,j)-Y(i,j))        &
                                *(p(i-1,j,k)-p(i,j,k))/CY(i-1,j)
                        case("Vl2")
                            pSon(a,b,c)=p(i,j,k)+(CenterX(i,j)-X(i,j))        &
                                *(p(i,j-1,k)-p(i,j,k))/CX(i,j-1)
                        case("Vl3")
                            pSon(a,b,c)=p(i,j,k)+(CenterZ(i,j,k)-Z(i,j,k))    &
                                *(p(i,j,k-1)-p(i,j,k))/CZ(i,j,k-1)
                        case("Vu1")
                            pSon(a,b,c)=p(i,j,k)+(Y(i,j)-CenterY(i,j))        &
                                *(p(i+1,j,k)-p(i,j,k))/CY(i,j)
                        case("Vu2")
                            pSon(a,b,c)=p(i,j,k)+(X(i,j)-CenterX(i,j))        &
                                *(p(i,j+1,k)-p(i,j,k))/CX(i,j)
                        case("Vu3")
                            pSon(a,b,c)=p(i,j,k)+(Z(i,j,k)-CenterZ(i,j,k))    &
                                *(p(i,j,k+1)-p(i,j,k))/CZ(i,j,k)
                        case default
                            stop 'FatherToSonGrid_3D - ModuleInterpolation - ERR03'
                        end select
                    case default
                        stop 'FatherToSonGrid_3D - ModuleInterpolation - ERR04'
                    end select
                endif
            enddo
        endif
        enddo
        enddo

    end subroutine FatherToSonGrid_3D

    !----------------------------------------------------------------------------

    subroutine Integrator_3D

        !Arguments---------------------------------------------------------------
        
        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                              :: n_dimensions
        integer                              :: i,j,k
        integer                              :: ilb, jlb, klb
        integer                              :: iub, jub, kub
        integer, pointer, dimension(:,:)     :: WaterPoints2D
        integer, pointer, dimension(:,:,:)   :: OpenPoints3D
        real,       pointer, dimension(:,:,:)   :: p, v
        real,       pointer, dimension(:,:)     :: CX, CY, DX, DY
        real,       pointer, dimension(:,:,:)   :: CZ, DZ

        !Begin shorten variables names-------------------------------------------
        
        n_dimensions    = Me%ComputeOptions%n_dimensions

        ilb             = Me%ExternalVar%SonGrid%WorkSize3D%ILB
        jlb             = Me%ExternalVar%SonGrid%WorkSize3D%JLB
        klb             = Me%ExternalVar%SonGrid%WorkSize3D%KLB
        iub             = Me%ExternalVar%SonGrid%WorkSize3D%IUB
        jub             = Me%ExternalVar%SonGrid%WorkSize3D%JUB
        kub             = Me%ExternalVar%SonGrid%WorkSize3D%KUB

        v               => Me%ExternalVar%SonField%Values3D
        p               => Me%ExternalVar%SonGrid%p

        CX              => Me%ExternalVar%SonGrid%CX
        CY              => Me%ExternalVar%SonGrid%CY
        CZ              => Me%ExternalVar%SonGrid%CZ
        DX              => Me%ExternalVar%SonGrid%DX
        DY              => Me%ExternalVar%SonGrid%DY
        DZ              => Me%ExternalVar%SonGrid%DZ

        WaterPoints2D   => Me%ExternalVar%SonGrid%WaterPoints2D
        OpenPoints3D    => Me%ExternalVar%SonGrid%OpenPoints3D

        !Begin-------------------------------------------------------------------

        do i=ilb,iub
        do j=jlb,jub
        if(WaterPoints2D(i,j).eq.1) then
            do k=1,kub
                if(OpenPoints3D(i,j,k).eq.1) then
                    select case (n_dimensions)
                    case(1)
                        v(i,j,k)=(CY(i,j)*(DY(i,j)*(p(i-1,j,k)-p(i,j,k))        &
                            +8d0*CY(i-1,j)*p(i,j,k))+CY(i-1,j)*DY(i,j)          &
                            *(-p(i,j,k)+p(i+1,j,k)))/(8d0*CY(i-1,j)*CY(i,j))
                    case(2)
                        v(i,j,k)=(CX(i,j)*CY(i-1,j)*CY(i,j)*DX(i,j)*(p(i,j-1,k) &
                            -p(i,j,k))+CX(i,j-1)*(CY(i-1,j)*CY(i,j)*DX(i,j)     &
                            *(-p(i,j,k)+p(i,j+1,k))+CX(i,j)*(CY(i,j)*(DY(i,j)   &
                            *(p(i-1,j,k)-p(i,j,k))+12d0*CY(i-1,j)*p(i,j,k))     &
                            +CY(i-1,j)*DY(i,j)*(-p(i,j,k)+p(i+1,j,k)))))/(12d0  &
                            *CX(i,j-1)*CX(i,j)*CY(i-1,j)*CY(i,j))
                    case(3)
                        v(i,j,k)=((DZ(i,j,k)*p(i,j,k-1))/CZ(i,j,k-1)+(DY(i,j)   &
                            *(p(i-1,j,k)-p(i,j,k)))/CY(i-1,j)+(DX(i,j)          &
                            *(p(i,j-1,k)-p(i,j,k)))/CX(i,j-1)+16d0*p(i,j,k)     &
                            -(DX(i,j)*p(i,j,k))/CX(i,j)-(DY(i,j)*p(i,j,k))      &
                            /CY(i,j)-(DZ(i,j,k)*p(i,j,k))/CZ(i,j,k-1)           &
                            -(DZ(i,j,k)*p(i,j,k))/CZ(i,j,k)+(DZ(i,j,k)          &
                            *p(i,j,k+1))/CZ(i,j,k)+(DX(i,j)*p(i,j+1,k))         &
                            /CX(i,j)+(DY(i,j)*p(i+1,j,k))/CY(i,j))/16d0
                    case default
                        stop 'Integrator_3D - ModuleInterpolation - ERR01'
                    end select
                endif
            enddo
        endif
        enddo
        enddo

    end subroutine Integrator_3D

    !----------------------------------------------------------------------------



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !----------------------------------------------------------------------------

    subroutine Kill_Interpolation3D(ObjInterpolationID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjInterpolationID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjInterpolationID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mInterpolation_,  Me%InstanceID)

            if (nUsers == 0) then

                !Classical methods-------------------------------------------------------------------
cd22 :          if(.not. Me%ConvolutionApproach) then

                    !Classical and Structured--------------------------------------------------------
cd23 :              if(Me%StructuredData) then

                        call Kill_NonConvoStruct

                    !Classical and Nonstructured-----------------------------------------------------
                    else

                        call Kill_NonConvoNonStruct

                    end if cd23

                !Convolution methods-----------------------------------------------------------------
                else

                    !NonConservative Convolution-----------------------------------------------------
cd24 :              if(Me%NonConservative) then

                        !NonConservative and Structured--------------------------------------------------------
cd25 :                  if(Me%StructuredData) then

                            call Kill_ConvoNCStructured

                        !NonConservative and Nonstructured-----------------------------------------------------
                        else
        
                            call Kill_ConvoNCNonStructured

                        end if cd25

                    !Conservative Convolution-----------------------------------------------------
                    else

                        call Kill_ConvoConservative

                    end if cd24

                end if cd22
            
                !Unassociate keywords
                call UnreadKeywordOptions

                !Deallocates Instance
                call DeallocateInstance ()

                ObjInterpolationID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           
    end subroutine Kill_Interpolation3D

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine Kill_ConvoConservative

        call DeallocateArrays
        call KillLocalFather
        call KillGrid(Me%ExternalVar%FatherGrid)
        call KillExternalField(Me%ExternalVar%FatherField)
        call KillGrid(Me%ExternalVar%SonGrid)
        call KillExternalField(Me%ExternalVar%SonField)
        call DeallocateExternal_ConvoStruct

    end subroutine Kill_ConvoConservative

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine DeallocateExternal_ConvoStruct

        call DeallocateField(Me%ExternalVar%FatherField)
        call DeallocateGrid(Me%ExternalVar%FatherGrid)
        call DeallocateField(Me%ExternalVar%SonField)
        call DeallocateGrid(Me%ExternalVar%SonGrid)
        call DeallocateExternal

    end subroutine DeallocateExternal_ConvoStruct

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine KillExternalField (FieldPointer)

        !Arguments---------------------------------------------------------------
        type (T_Field), pointer                          :: FieldPointer
      
        if (associated(Fieldpointer%Values2D)) then
            nullify(Fieldpointer%Values2D)
        end if

        if (associated(Fieldpointer%Values3D)) then
            nullify(Fieldpointer%Values3D)
        end if

    end subroutine KillExternalField

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine KillGrid(Gridpointer)

        !Arguments---------------------------------------------------------------
        type (T_Grid), pointer                          :: GridPointer

        call DeassociateGridInstances(Gridpointer)
        call ReadUnlockExternalVar(Gridpointer)
        
        if (associated(Gridpointer%Size3D)) then
            deallocate(Gridpointer%Size3D)
        end if

        if (associated(Gridpointer%WorkSize3D)) then
            deallocate(Gridpointer%WorkSize3D)
        end if

        if (associated(Gridpointer%Size2D)) then
            deallocate(Gridpointer%Size2D)
        end if

        if (associated(Gridpointer%WorkSize2D)) then
            deallocate(Gridpointer%WorkSize2D)
        end if

    end subroutine KillGrid

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine KillLocalFather

        call DeallocateLocalFatherGrid(Me%LocalFatherGrid)
        call DeallocateLocalFatherField(Me%LocalFatherField)

        deallocate(Me%LocalFatherGrid)
        deallocate(Me%LocalFatherField)

        nullify(Me%LocalFatherGrid)
        nullify(Me%LocalFatherField)

    end subroutine KillLocalFather

    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------

    subroutine DeallocateLocalFatherField(Fieldpointer)

    !Locals----------------------------------------------------
    type(T_Field)                                :: Fieldpointer
    integer                                      :: STAT_CALL

        if(associated(Fieldpointer%Values3D)) then
            deallocate(Fieldpointer%Values3D, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateNC - Interpolation - ERR01'
        end if

        if(associated(Fieldpointer%Values2D)) then
            deallocate(Fieldpointer%Values2D, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateNC - Interpolation - ERR01'
        end if

    end subroutine DeallocateLocalFatherField

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine DeallocateLocalFatherGrid(Gridpointer)

    !Locals----------------------------------------------------
    type(T_Grid),pointer                         :: Gridpointer
    integer                                      :: STAT_CALL

        if(associated(Gridpointer%OpenPoints3D)) then
            deallocate(Gridpointer%OpenPoints3D, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateLocalFatherGrid - Interpolation - ERR01'
        end if

        if(associated(Gridpointer%WaterPoints3D)) then
            deallocate(Gridpointer%WaterPoints3D, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateLocalFatherGrid - Interpolation - ERR01'
        end if

        if(associated(Gridpointer%WaterPoints2D)) then
            deallocate(Gridpointer%WaterPoints2D, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateLocalFatherGrid - Interpolation - ERR01'
        end if

        if(associated(Gridpointer%CenterX)) then
            deallocate(Gridpointer%CenterX, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateLocalFatherGrid - Interpolation - ERR01'
        end if

        if(associated(Gridpointer%CenterY)) then
            deallocate(Gridpointer%CenterY, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateLocalFatherGrid - Interpolation - ERR01'
        end if

        if(associated(Gridpointer%CenterZ)) then
            deallocate(Gridpointer%CenterZ, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateLocalFatherGrid - Interpolation - ERR01'
        end if
        
        if(associated(Gridpointer%Size3D)) then
            deallocate(Gridpointer%Size3D, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateLocalFatherGrid - Interpolation - ERR01'
        end if

        if(associated(Gridpointer%WorkSize3D)) then
            deallocate(Gridpointer%WorkSize3D, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateLocalFatherGrid - Interpolation - ERR01'
        end if

        if(associated(Gridpointer%Size2D)) then
            deallocate(Gridpointer%Size2D, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateLocalFatherGrid - Interpolation - ERR01'
        end if

        if(associated(Gridpointer%WorkSize2D)) then
            deallocate(Gridpointer%WorkSize2D, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateLocalFatherGrid - Interpolation - ERR01'
        end if

    end subroutine DeallocateLocalFatherGrid

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine Kill_ConvoNCStructured

        call Kill_NC(Me%NC)
        nullify(Me%NC)
        call KillGrid(Me%ExternalVar%FatherGrid)   
        call KillExternalField(Me%ExternalVar%FatherField)
        call KillGrid(Me%ExternalVar%SonGrid)   
        call KillExternalField(Me%ExternalVar%SonField)
        call DeallocateExternal_ConvoStruct

    end subroutine Kill_ConvoNCStructured

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine Kill_NC(NC)
        
    !Locals----------------------------------------------------
    type(T_NC), pointer                          :: NC
    integer                                      :: STAT_CALL

        if(associated(NC%GroupCenterX)) then
            Deallocate(NC%GroupCenterX, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillNC - Interpolation - ERR01'
        end if

        if(associated(NC%GroupCenterY)) then
            Deallocate(NC%GroupCenterY, stat = STAT_CALL)
        end if

        if(associated(NC%GroupCenterZ)) then
            Deallocate(NC%GroupCenterZ, stat = STAT_CALL)
        end if

        if(associated(NC%Tau)) then
            Deallocate(NC%Tau, stat = STAT_CALL)
        end if

        if(associated(NC%TauR)) then
            Deallocate(NC%TauR, stat = STAT_CALL)
        end if

        if(associated(NC%TauM)) then
            Deallocate(NC%TauM, stat = STAT_CALL)
        end if

        if(associated(NC%Response)) then
            Deallocate(NC%Response, stat = STAT_CALL)
        end if

        if(associated(NC%Predictor)) then
            Deallocate(NC%Predictor, stat = STAT_CALL)
        end if

    end subroutine Kill_NC

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine Kill_ConvoNCNonStructured

        call Kill_NC(Me%NC)
        call KillGrid(Me%ExternalVar%SonGrid)
        call KillExternalField(Me%ExternalVar%SonField)
        call KillPoints2Grid_3D
        call DeallocateTXYZP_Points
        call DeallocExtern_ConvoNonStruct

    end subroutine Kill_ConvoNCNonStructured

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine DeallocExtern_ConvoNonStruct

        call DeallocateField(Me%ExternalVar%SonField)
        call DeallocateGrid(Me%ExternalVar%SonGrid)
        call DeallocateExternal

    end subroutine DeallocExtern_ConvoNonStruct

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine DeallocateTXYZP_Points

        if (associated(Me%ExternalVar%TXYZP_Points)) nullify(Me%ExternalVar%TXYZP_Points)

    end subroutine DeallocateTXYZP_Points

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine KillPoints2Grid_3D

        call KillLocalFather

    end subroutine KillPoints2Grid_3D

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine Kill_NonConvoNonStruct

    !TOWRITE

    end subroutine Kill_NonConvoNonStruct

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine Kill_NonConvoStruct

    !TOWRITE

    end subroutine Kill_NonConvoStruct

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine UnreadKeywordOptions

    !TOWRITE

    end subroutine UnreadKeywordOptions

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    subroutine DeassociateGridInstances (GridPointer)
        !Arguments---------------------------------------------------------------
        type (T_Grid), pointer                          :: GridPointer

        !Local-------------------------------------------------------------------
        integer                                         :: nUsers

        !------------------------------------------------------------------------
    
        nUsers = DeAssociateInstance (mHORIZONTALGRID_, GridPointer%ObjHorizontalGrid)
        if (nUsers == 0) stop 'DeassociateGridInstances - ModuleInterpolation - ERR10'

        nUsers = DeAssociateInstance (mGEOMETRY_, GridPointer%ObjGeometry)
        if (nUsers == 0) stop 'DeassociateGridInstances - ModuleInterpolation - ERR20'

        nUsers = DeAssociateInstance (mHORIZONTALMAP_, GridPointer%ObjHorizontalMap)
        if (nUsers == 0) stop 'DeassociateGridInstances - ModuleInterpolation - ERR30'

        nUsers = DeAssociateInstance (mMAP_, GridPointer%ObjMap)
        if (nUsers == 0) stop 'DeassociateGridInstances - ModuleInterpolation - ERR30'
 
    end subroutine DeassociateGridInstances
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    subroutine DeassociateLocalFatherGrid 

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: nUsers
        !Begin-----------------------------------------------------------------

       
        nUsers = DeAssociateInstance (mTIME_,           Me%LocalFatherGrid%ObjTime)
        if (nUsers == 0) stop 'DeassociateLocalFatherGrid - ModuleInterpolation - ERR10'

        nUsers = DeAssociateInstance (mHORIZONTALGRID_, Me%LocalFatherGrid%ObjHorizontalGrid)
        if (nUsers == 0) stop 'DeassociateLocalFatherGrid - ModuleInterpolation - ERR20'

        nUsers = DeAssociateInstance (mGRIDDATA_,       Me%LocalFatherGrid%ObjBathymetry)
        if (nUsers == 0) stop 'DeassociateLocalFatherGrid - ModuleInterpolation - ERR30'
         
        nUsers = DeAssociateInstance (mHORIZONTALMAP_,  Me%LocalFatherGrid%ObjHorizontalMap)    
        if (nUsers == 0) stop 'DeassociateLocalFatherGrid - ModuleInterpolation - ERR40'

    end subroutine DeassociateLocalFatherGrid 
        
   !------------------------------------------------------------------------
        
   !------------------------------------------------------------------------

    subroutine DeallocateNC(NC)

        !Arguments-------------------------------------------------------------
        type (T_NC), pointer                            :: NC

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Deallocates 
        if (associated(NC%GroupCenterX)) then
            deallocate (NC%GroupCenterX, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateNC - Interpolation - ERR01'
        endif

        if (associated(NC%GroupCenterY)) then
            deallocate (NC%GroupCenterY, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateNC - Interpolation - ERR02'
        endif

        if (associated(NC%GroupCenterZ)) then
            deallocate (NC%GroupCenterZ, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateNC - Interpolation - ERR03'
        endif

        if (associated(NC%TauR)) then
            deallocate (NC%TauR, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateNC - Interpolation - ERR04'
        endif

        if (associated(NC%TauM)) then
            deallocate (NC%TauM, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateNC - Interpolation - ERR05'
        endif

        if (associated(NC%Tau)) then
            deallocate (NC%Tau, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateNC - Interpolation - ERR06'
        endif

        if (associated(NC%Response)) then
            deallocate (NC%Response, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateNC - Interpolation - ERR07'
        endif

        if (associated(NC%Predictor)) then
            deallocate (NC%Predictor, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateNC - Interpolation - ERR08'
        endif

        if (associated(NC)) then
            deallocate (NC, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateNC - Interpolation - ERR09'
        endif

    end subroutine DeallocateNC

    !--------------------------------------------------------------------------

   !------------------------------------------------------------------------

    subroutine DeallocateArrays

        !Parameter-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                         :: STAT_CALL

        !Deallocates 
        if (associated(Me%MatrixUpdate3D)) then
            deallocate (Me%MatrixUpdate3D, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateArrays - Interpolation - ERR01'
        endif

        if (associated(Me%FatherSonX)) then
            deallocate (Me%FatherSonX, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateArrays - Interpolation - ERR01'
        endif

        if (associated(Me%FatherSonY)) then
            deallocate (Me%FatherSonY, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateArrays - Interpolation - ERR02'
        endif

        if (associated(Me%FatherSonZ)) then
            deallocate (Me%FatherSonZ, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateArrays - Interpolation - ERR03'
        endif

        if (associated(Me%InsideFatherCell)) then
            deallocate (Me%InsideFatherCell, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateArrays - Interpolation - ERR04'
        endif

    end subroutine DeallocateArrays

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine DeallocateGridArrays(GridPointer)

        !Arguments-------------------------------------------------------------
        type (T_Grid), pointer                          :: GridPointer

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Deallocates 
        if (associated(GridPointer%XX_IE)) then
            deallocate (GridPointer%XX_IE, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateGridArrays - Interpolation - ERR01'
        endif

        !Deallocates 
        if (associated(GridPointer%YY_IE)) then
            deallocate (GridPointer%YY_IE, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateGridArrays - Interpolation - ERR01'
        endif

        !Deallocates 
        if (associated(GridPointer%CenterX)) then
            deallocate (GridPointer%CenterX, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateGridArrays - Interpolation - ERR01'
        endif

        !Deallocates 
        if (associated(GridPointer%CenterY)) then
            deallocate (GridPointer%CenterY, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateGridArrays - Interpolation - ERR01'
        endif

        !Deallocates 
        if (associated(GridPointer%CenterZ)) then
            deallocate (GridPointer%CenterZ, stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateGridArrays - Interpolation - ERR01'
        endif

        !Deallocates 
        if (associated(GridPointer%DX)) then
             deallocate (GridPointer%DX, stat = STAT_CALL)
             if (STAT_CALL /= SUCCESS_) stop 'DeallocateGridArrays - Interpolation - ERR01'
        endif
                
        !Deallocates 
        if (associated(GridPointer%DY)) then
           deallocate (GridPointer%DY, stat = STAT_CALL)
           if (STAT_CALL /= SUCCESS_) stop 'DeallocateGridArrays - Interpolation - ERR01'
        endif
           
        !Deallocates 
        if (associated(GridPointer%CX)) then
           deallocate (GridPointer%CX, stat = STAT_CALL)
           if (STAT_CALL /= SUCCESS_) stop 'DeallocateGridArrays - Interpolation - ERR01'
        endif

        !Deallocates 
        if (associated(GridPointer%CY)) then
           deallocate (GridPointer%CY, stat = STAT_CALL)
           if (STAT_CALL /= SUCCESS_) stop 'DeallocateGridArrays - Interpolation - ERR01'
        endif

        !Deallocates 
        if (associated(GridPointer%DZ)) then
           deallocate (GridPointer%DZ, stat = STAT_CALL)
           if (STAT_CALL /= SUCCESS_) stop 'DeallocateGridArrays - Interpolation - ERR01'
        endif

        !Deallocates 
        if (associated(GridPointer%CZ)) then
           deallocate (GridPointer%CZ, stat = STAT_CALL)
           if (STAT_CALL /= SUCCESS_) stop 'DeallocateGridArrays - Interpolation - ERR01'
        endif

        !Deallocates 
        if (associated(GridPointer%SZZ)) then
           deallocate (GridPointer%SZZ, stat = STAT_CALL)
           if (STAT_CALL /= SUCCESS_) stop 'DeallocateGridArrays - Interpolation - ERR01'
        endif

    end subroutine DeallocateGridArrays

    !--------------------------------------------------------------------------

    !------------------------------------------------------------------------
       
    subroutine DeallocateGrid (GridPointer)

        !Arguments-------------------------------------------------------------
        type (T_Grid), pointer                          :: GridPointer

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL


        if (associated(GridPointer)) then
            deallocate (GridPointer, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateGrid - Interpolation - ERR01'
            nullify (GridPointer)
        end if
            
    end subroutine DeallocateGrid

    !--------------------------------------------------------------------------

    !------------------------------------------------------------------------
       
    subroutine DeallocateField (FieldPointer)

        !Arguments-------------------------------------------------------------
        type (T_Field), pointer                          :: FieldPointer

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL


        if (associated(FieldPointer)) then
            deallocate (FieldPointer, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateField - Interpolation - ERR01'
            nullify (FieldPointer)
        end if
            
    end subroutine DeallocateField

    !--------------------------------------------------------------------------

    !------------------------------------------------------------------------
       
    subroutine DeallocateExternal

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        if (associated(Me%ExternalVar)) then
            deallocate (Me%ExternalVar, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeallocateExternal - Interpolation - ERR01'
            nullify (Me%ExternalVar)
        end if
            
    end subroutine DeallocateExternal

    !--------------------------------------------------------------------------

    !------------------------------------------------------------------------
       
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Interpolation), pointer          :: AuxObjInterpolation
        type (T_Interpolation), pointer          :: PreviousObjInterpolation

        !Updates pointers
        if (Me%InstanceID == FirstObjInterpolation%InstanceID) then
            FirstObjInterpolation => FirstObjInterpolation%Next
        else
            PreviousObjInterpolation => FirstObjInterpolation
            AuxObjInterpolation      => FirstObjInterpolation%Next
            do while (AuxObjInterpolation%InstanceID /= Me%InstanceID)
                PreviousObjInterpolation => AuxObjInterpolation
                AuxObjInterpolation      => AuxObjInterpolation%Next
            enddo

            !Now update linked list
            PreviousObjInterpolation%Next => AuxObjInterpolation%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadUnlockExternalVar(GridPointer)
        
        !External--------------------------------------------------------------
        type (T_Grid), pointer                  :: GridPointer
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        !WaterPoints2D
        call UngetHorizontalMap(GridPointer%ObjHorizontalMap, GridPointer%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleInterpolation - ERR02'

        !OpenPoints3D
        call UngetMap(GridPointer%ObjMap, GridPointer%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleInterpolation - ERR03'
        
        !FatherGrid
        select case (Me%ComputeOptions%TypeZUV)

            case(TypeZ_)
                call ReadUnlockExternalZUV(GridPointer)

            case(TypeU_)
                call ReadUnlockExternalZUV(GridPointer)

            case(TypeV_)
                call ReadUnlockExternalZUV(GridPointer)
            
            case default
                
                write(*,*)'Invalid type ZUV in grid data '
                stop 'ReadUnlockExternalVar - ModuleInterpolation - ERR90'

        end select

    end subroutine ReadUnlockExternalVar

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    Subroutine ReadUnlockExternalZUV(GridPointer)
    
        !Arguments------------------------------------------------
        type (T_Grid), pointer              :: GridPointer


        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !XX_IE and YY_IE
        if(associated(GridPointer%XX_IE)) then
            call UnGetHorizontalGrid (GridPointer%ObjHorizontalGrid,            &
                                    GridPointer%XX_IE,                           &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleInterpolation - ERR04'
        end if
      
        if(associated(GridPointer%YY_IE)) then
            call UnGetHorizontalGrid (GridPointer%ObjHorizontalGrid,            &
                                    GridPointer%YY_IE,                           &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleInterpolation - ERR04'
        end if

        !DX
        if(associated(GridPointer%DX)) then
            call UnGetHorizontalGrid (GridPointer%ObjHorizontalGrid, GridPointer%DX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleInterpolation - ERR08'
        end if
        
        !DY
        if(associated(GridPointer%DY)) then
            call UnGetHorizontalGrid (GridPointer%ObjHorizontalGrid, GridPointer%DY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleInterpolation - ERR09'
        end if

        !DZ
        if(associated(GridPointer%DZ)) then
            call UnGetGeometry (GridPointer%ObjGeometry, GridPointer%DZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleInterpolation - ERR07b'
        end if

        !CX
        if(associated(GridPointer%CX)) then
            call UnGetHorizontalGrid (GridPointer%ObjHorizontalGrid, GridPointer%CX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleInterpolation - ERR09'
        end if

        !CY
        if(associated(GridPointer%CY)) then
            call UnGetHorizontalGrid (GridPointer%ObjHorizontalGrid, GridPointer%CY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleInterpolation - ERR09'
        end if

        !CZ
        if(associated(GridPointer%CZ)) then
            call UnGetGeometry (GridPointer%ObjGeometry, GridPointer%CZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleInterpolation - ERR09'
        end if

        !SZZ
        if(associated(GridPointer%SZZ)) then
            call UnGetGeometry (GridPointer%ObjGeometry, GridPointer%SZZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleInterpolation - ERR09'
        end if

    end subroutine ReadUnlockExternalZUV
  
    !--------------------------------------------------------------------------
  
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjInterpolation_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjInterpolation_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjInterpolation_ID > 0) then
            call LocateObjInterpolation (ObjInterpolation_ID)
            ready_ = VerifyReadLock (mInterpolation_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjInterpolation (ObjInterpolationID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjInterpolationID

        !Local-----------------------------------------------------------------

        Me => FirstObjInterpolation
        do while (associated (Me))
            if (Me%InstanceID == ObjInterpolationID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleInterpolation - LocateObjInterpolation - ERR02'

    end subroutine LocateObjInterpolation


end module ModuleInterpolation

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
