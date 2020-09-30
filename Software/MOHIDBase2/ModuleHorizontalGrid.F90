!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : HorizontalGrid
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to read and calculate a horizontal Grid
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

Module ModuleHorizontalGrid

    use ModuleGlobalData

    use ModuleFunctions, only   : Rodaxy, FromCartesianToGrid, FromGridToCartesian,     &
                                  UTMToLatLon, LatLonToUTM, ComputeGridZone,            &
                                  LatLonToLambertSP2, RelativePosition4VertPolygon,     &
                                  CHUNK_J, WGS84toGoogleMaps, AngleFromFieldToGrid,     &
                                  AngleFromGridToField, THOMAS_2D, THOMAS_3D,           &
                                  SphericalToCart
#ifdef _USE_PROJ4
    use ModuleFunctions, only   : GeographicToCartesian, CartesianToGeographic
#endif _USE_PROJ4

#ifdef _USE_MPI
    use ModuleFunctions, only   : MPIKind
    use mpi
#endif _USE_MPI

    use ModuleStopWatch, only   : StartWatch, StopWatch
    use ModuleEnterData
    use ModuleDrawing
    use ModuleHDF5


    implicit none
    private

    !Subroutine----------------------------------------------------------------

    !Constructor
    public  :: ConstructHorizontalGrid
    private ::      AllocateInstance
    private ::      ConstructGlobalVariables
    private ::          AllocateVariables
    private ::          Mercator
    private ::              ConvertCoordenates     !Fromer CONVCO_BAT
    private ::          GridPointsLocation1D
    private ::          GridPointsLocation2D
    private ::          ComputeDistances
    private ::          ComputeRotation
    private ::          ComputeGridCellArea
    private ::          ComputeCoriolis
    public  :: ConstructFatherGridLocation
    private ::      ConstructNewFatherGrid1D
    private ::      ConstructNewFatherGrid2D
    private ::      DetermineMaxRatio
    private ::      Add_FatherGrid
    private :: CheckGridBorder
    private :: DefineBorderPolygons
    private ::      DefinesBorderPoly
    public  :: ConstructP2C_IWD
    public  ::      ConstructP2C_Avrg
    private ::      ConstructIWDVel
    
    !Modifier
    public  :: WriteHorizontalGrid
    public  :: WriteHorizontalGrid_UV
    public  :: LocateCell
    public  :: LocateCell1D
    public  :: LocateCellPolygons
    public  :: RecenterHorizontalGrid
    public  :: Add_MPI_ID_2_Filename
    public  :: No_BoxMap_HaloArea

    !Selector
    public  :: GetHorizontalGridSize
    public  :: GetHorizontalGrid
    public  :: GetGridOrigin
    public  :: GetGridAngle
    public  :: GetGridZone
    public  :: GetLatitudeLongitude
    public  :: GetGridCoordType
    public  :: GetCoordTypeList
    public  :: GetGridMeanLatLong
    public  :: GetCoriolisFrequency
    public  :: GetGridLatitudeLongitude
    public  :: GetComputeZUV
    public  :: GetGridFileName
    public  :: GetCheckDistortion
    public  :: GetGridRotation
    public  :: GetCellRotation
    public  :: GetDefineCellsMap
    public  :: GetNotDefinedCells
    public  :: GetZCoordinates
    public  :: GetCornersCoordinates
    public  :: GetFatherGridID
    public  :: GetGridCellArea
    private ::      Search_FatherGrid
    public  :: GetCellIDfromIJ
    public  :: GetCellIJfromID
    public  :: GetCellIDfromIJArray
    public  :: GetXYCellZ
    public  :: GetXYCellZ_ThreadSafe
    public  :: GetXYArrayIJ
    public  :: GetCellZ_XY
    public  :: GetCellZInterceptByLine
    public  :: GetCellZInterceptByPolygon
    public  :: GetCellZInterceptByXYZPoints

    public  :: GetGridBorderPolygon
    public  :: GetGridOutBorderPolygon
    public  :: GetGridBorderCartPolygon
    public  :: GetGridOutBorderCartPolygon
    public  :: GetGridBorderLimits
    public  :: GetGridOutBorderCartLimits

    public  :: GetXYInsideDomain
    private ::      InsideDomainPolygon

    public  :: GetGridBorderType
    public  :: GetDDecompOpenBorders
    public  :: GetDDecompParameters
    public  :: GetDDecompSlaves
    public  :: GetDDecompSlavesSize
    public  :: GetDDecompWorkSize2D
    public  :: GetDDecompMapping2D
    public  :: GetDDecompMPI_ID
    public  :: GetDDecompON

    public  :: GetSonWindow

    public  :: GetConnections
    public  :: UnGetConnections


    public  :: UnGetHorizontalGrid

    public  :: InterpolXYPoint

    public  :: ReturnsIntersectionCorners
    private ::      ComputesIntersectionCorners

    public  :: WindowIntersectDomain
    private ::      WindowCellsIntersection

#ifdef _USE_PROJ4
    public  :: FromGeo2SpherMercator1D
    public  :: FromGeo2SpherMercatorScalar
    public  :: FromGeographic2SphericMercator
#endif

    !Destructor
    public  :: KillHorizontalGrid
    private ::      DeallocateInstance
    private ::      KillFatherGridList


    !Management
    private ::      Ready

    !Auxilary
    private :: USCONVCO
    private :: UTGP83
    private :: GPUT83
    private :: IUTPC
    private :: IUTGP
    private :: DRGPUT
    private :: DRUTGP
    private :: TODMS
    private :: DATUMM
    private :: TMGRID
    private :: TCONST
    private :: TCONPC
    private :: TMGEOD

    private :: wgs84_to_rd
    private :: wgs842bessel_
    private :: bessel2rd_
    private :: rd_to_wgs84
    private :: rd2bessel_
    private :: bessel2wgs84_

    !Interfaces----------------------------------------------------------------

    interface  ConstructHorizontalGrid
        module procedure ConstructHorizontalGridV1
        module procedure ConstructHorizontalGridV2
        module procedure ConstructHorizontalGridV3
    end interface  ConstructHorizontalGrid

    interface  LocateCell1D
        module procedure LocateCell1DR4
        module procedure LocateCell1DR8
    end interface LocateCell1D


    private :: UnGetHorizontalGrid1d
    private :: UnGetHorizontalGrid2d
    private :: UngetHorizontalGrid1DInt
    private :: UngetHorizontalGrid2DInt
    private :: UnGetHorizontalGridPolygon
    interface  UnGetHorizontalGrid
        module procedure UnGetHorizontalGrid1D
        module procedure UnGetHorizontalGrid2D
        module procedure UnGetHorizontalGrid1DInt
        module procedure UnGetHorizontalGrid2DInt
        module procedure UnGetHorizontalGridPolygon
    end interface  UnGetHorizontalGrid


    public :: InterpolRegularGrid
    interface  InterpolRegularGrid
        module procedure InterpolRegularGrid2D
        module procedure InterpolRegularGrid3D
        module procedure InterpolRegularGrid3D8
    end interface  InterpolRegularGrid

    private:: RotateVectorFieldToGrid2D
    private:: RotateVectorFieldToGrid3D
    public :: RotateVectorFieldToGrid
    interface  RotateVectorFieldToGrid
        module procedure RotateVectorFieldToGrid2D
        module procedure RotateVectorFieldToGrid3D
    end interface  RotateVectorFieldToGrid

    private:: RotateVectorGridToField2DR4
    private:: RotateVectorGridToField2DR8
    private:: RotateVectorGridToField3D
    public :: RotateVectorGridToField
    interface  RotateVectorGridToField
        module procedure RotateVectorGridToField2DR4
        module procedure RotateVectorGridToField2DR8
        module procedure RotateVectorGridToField3D
    end interface  RotateVectorGridToField

    private:: RotateAngleFieldToGrid2D
    private:: RotateAngleFieldToGrid3D
    public :: RotateAngleFieldToGrid
    interface  RotateAngleFieldToGrid
        module procedure RotateAngleFieldToGrid2D
        module procedure RotateAngleFieldToGrid3D
    end interface  RotateAngleFieldToGrid

    private:: RotateAngleGridToField2D
    private:: RotateAngleGridToField3D
    public :: RotateAngleGridToField
    interface  RotateAngleGridToField
        module procedure RotateAngleGridToField2D
        module procedure RotateAngleGridToField3D
    end interface  RotateAngleGridToField

    private:: ComputeAngleFromGridComponents2D
    private:: ComputeAngleFromGridComponents3D
    public :: ComputeAngleFromGridComponents
    interface  ComputeAngleFromGridComponents
        module procedure ComputeAngleFromGridComponents2D
        module procedure ComputeAngleFromGridComponents3D
    end interface  ComputeAngleFromGridComponents

#ifdef _USE_MPI
    public :: THOMAS_DDecompHorizGrid
    interface THOMAS_DDecompHorizGrid
        module procedure THOMAS_2D_DDecompHorizGrid
        module procedure THOMAS_3D_DDecompHorizGrid
    end interface

    public :: JoinGridData
    interface JoinGridData
        module procedure JoinGridData_1D
        module procedure JoinGridData_2D_R4
        module procedure JoinGridData_2D_R8
        module procedure JoinGridData_2Dint
        module procedure JoinGridData_3D_R4
        module procedure JoinGridData_3D_R8
        module procedure JoinGridData_3Dint
    end interface

    public :: JoinGridData_In
    interface JoinGridData_In
        module procedure JoinGridData_1D_In
        module procedure JoinGridData_2D_R4_In
        module procedure JoinGridData_2D_R8_In
        module procedure JoinGridData_2Dint_In
        module procedure JoinGridData_3D_R4_In
        module procedure JoinGridData_3D_R8_In
        module procedure JoinGridData_3Dint_In
    end interface

    public :: BroadcastGridData
    interface BroadcastGridData
        module procedure BroadcastGridData2D_R4
        module procedure BroadcastGridData2D_R8
        module procedure BroadcastGridData3D_R4
        module procedure BroadcastGridData3D_R8
    end interface

    public :: BroadcastGridData_In
    interface BroadcastGridData_In
        module procedure BroadcastGridData2D_R4_In
        module procedure BroadcastGridData2D_R8_In
        module procedure BroadcastGridData3D_R4_In
        module procedure BroadcastGridData3D_R8_In
    end interface

    public :: ReceiveSendProperitiesMPI
    interface ReceiveSendProperitiesMPI
        module procedure ReceiveSendProperities3DMPIr4
        module procedure ReceiveSendProperities2DMPIr4
        module procedure ReceiveSendProperities3DMPIr8
        module procedure ReceiveSendProperities2DMPIr8
    endinterface

    public :: ReceiveSendLogicalMPI
    interface ReceiveSendLogicalMPI
        module procedure AtLeastOneDomainIsTrue
    endinterface
    
    public :: GetKfloorZminMPI

#endif _USE_MPI



    !Parameter-----------------------------------------------------------------


    !Block information
    character(LEN = StringLength), parameter :: BeginXX         = '<BeginXX>'
    character(LEN = StringLength), parameter :: EndXX           = '<EndXX>'
    character(LEN = StringLength), parameter :: BeginYY         = '<BeginYY>'
    character(LEN = StringLength), parameter :: EndYY           = '<EndYY>'
    character(LEN = StringLength), parameter :: BeginCornersXY  = '<CornersXY>'
    character(LEN = StringLength), parameter :: EndCornersXY    = '<'//backslash//'CornersXY>'
    character(LEN = StringLength), parameter :: BeginCartCornersXY  = '<CartCornersXY>'
    character(LEN = StringLength), parameter :: EndCartCornersXY    = '<'//backslash//'CartCornersXY>'


    !Calculation Points
    integer, parameter :: ComputeZ_     = 1
    integer, parameter :: ComputeU_     = 2
    integer, parameter :: ComputeV_     = 3
    integer, parameter :: ComputeCross_ = 4
    integer, parameter :: ComputeZU_    = 5
    integer, parameter :: ComputeZV_    = 6


    integer, dimension(1:4), parameter :: diFace = (/0,1 ,0,-1/), djFace = (/-1, 0, 1,  0/)

    integer, parameter :: LargeScaleModel_ = 1
    integer, parameter :: Assimila_        = 2
    integer, parameter :: SurfaceMM5_      = 3


    !Input / Output
    integer, parameter :: FileOpen = 1, FileClose = 0


    !Type----------------------------------------------------------------------
    type T_BorderLimits
        real,    dimension(4)            :: Values = FillValueReal
        logical                          :: ON     = .false.
    end type T_BorderLimits    
    
    type T_Compute
        real,    dimension(:),   pointer :: XX_Z => null()
        real,    dimension(:),   pointer :: YY_Z => null()
        real,    dimension(:),   pointer :: XX_U => null()
        real,    dimension(:),   pointer :: YY_U => null()
        real,    dimension(:),   pointer :: XX_V => null()
        real,    dimension(:),   pointer :: YY_V => null()
        real,    dimension(:),   pointer :: XX_Cross => null()
        real,    dimension(:),   pointer :: YY_Cross => null()
        real,    dimension(:,:), pointer :: XX2D_Z => null()
        real,    dimension(:,:), pointer :: YY2D_Z => null()
        real,    dimension(:,:), pointer :: XX2D_U => null()
        real,    dimension(:,:), pointer :: YY2D_U => null()
        real,    dimension(:,:), pointer :: XX2D_V => null()
        real,    dimension(:,:), pointer :: YY2D_V => null()
    end type T_Compute

    type T_Window
        logical                          :: ON  = .false. !initialization: jauch
        integer                          :: ILB = null_int !initialization: jauch
        integer                          :: IUB = null_int !initialization: jauch
        integer                          :: JLB = null_int !initialization: jauch
        integer                          :: JUB = null_int !initialization: jauch
    end type T_Window

    type T_FatherGrid
        integer                          :: GridID = null_int !initialization: jauch - or should be set to 0 (zero)?
        real,    dimension(:,:), pointer :: XX_Z     => null()
        real,    dimension(:,:), pointer :: YY_Z     => null()
        real,    dimension(:,:), pointer :: XX_U     => null()
        real,    dimension(:,:), pointer :: YY_U     => null()
        real,    dimension(:,:), pointer :: XX_V     => null()
        real,    dimension(:,:), pointer :: YY_V     => null()
        real,    dimension(:,:), pointer :: XX_Cross => null()
        real,    dimension(:,:), pointer :: YY_Cross => null()
        logical                          :: OkZ            = .false. !initialization: jauch
        logical                          :: OkU            = .false. !initialization: jauch
        logical                          :: OkV            = .false. !initialization: jauch
        logical                          :: OkCross        = .false. !initialization: jauch
        logical                          :: CornersXYInput = .false.
        type (T_Window)                  :: Window
        real                             :: Fhc = null_real !initialization: jauch
        integer                          :: JX  = null_int !initialization: jauch
        integer                          :: IY  = null_int !initialization: jauch
        integer, dimension(:,:), pointer :: IZ     => null()
        integer, dimension(:,:), pointer :: JZ     => null()
        integer, dimension(:,:), pointer :: IU     => null()
        integer, dimension(:,:), pointer :: JU     => null()
        integer, dimension(:,:), pointer :: IV     => null()
        integer, dimension(:,:), pointer :: JV     => null()
        integer, dimension(:,:), pointer :: ICross => null()
        integer, dimension(:,:), pointer :: JCross => null()
        
        integer, dimension(:,:), pointer :: ILinkZ => null()
        integer, dimension(:,:), pointer :: JLinkZ => null()
        integer, dimension(:,:), pointer :: ILinkU => null()
        integer, dimension(:,:), pointer :: JLinkU => null()
        integer, dimension(:,:), pointer :: ILinkV => null()
        integer, dimension(:,:), pointer :: JLinkV => null()
        
        type (T_Size2D)                  :: MPI_Window
        type (T_FatherGrid),     pointer :: Next => null()
        type (T_FatherGrid),     pointer :: Prev => null()
    end type T_FatherGrid

    type T_Border
        !Grid boundary
        type(T_Polygon),          pointer       :: Polygon_ => null()
        integer                                 :: Type_    = null_int
    end type T_Border

    private :: T_Coef2D
    type       T_Coef2D
        logical                                 :: AllocateON       = .false.
        real(8), pointer, dimension(:)          :: VECG             => null()
        real(8), pointer, dimension(:)          :: VECW             => null()
        real,    pointer, dimension(:,:)        :: Results2D        => null()
        real,    pointer, dimension(:,:)        :: D                => null()
        real(8), pointer, dimension(:,:)        :: E                => null()
        real   , pointer, dimension(:,:)        :: F                => null()
        real   , pointer, dimension(:,:)        :: Ti               => null()
    end type T_Coef2D

    private :: T_Coef3D
    type       T_Coef3D
        logical                                 :: AllocateON       = .false.
        real(8), pointer, dimension(:)          :: VECG             => null()
        real(8), pointer, dimension(:)          :: VECW             => null()
        real,    pointer, dimension(:,:,:)      :: Results3D        => null()
        real,    pointer, dimension(:,:,:)      :: D                => null()
        real(8), pointer, dimension(:,:,:)      :: E                => null()
        real,    pointer, dimension(:,:,:)      :: F                => null()
        real,    pointer, dimension(:,:,:)      :: Ti               => null()
        real,    pointer, dimension(:,:,:)      :: C                => null()
        real,    pointer, dimension(:,:,:)      :: G                => null()
        integer                                 :: KLB              = FillValueInt
        integer                                 :: KUB              = FillValueInt
    end type T_Coef3D

    private :: T_DDecomp
    type       T_DDecomp
        logical                                 :: ON               = .false.
        logical                                 :: Master           = .false.
        logical                                 :: MasterOrSlave    = .false.
        logical                                 :: Auto             = .false.
        integer                                 :: Master_MPI_ID    = null_int
        integer                                 :: Nslaves          = 0
        integer, dimension(:), pointer          :: Slaves_MPI_ID    => null()
        type (T_Size2D), dimension(:), pointer  :: Slaves_Size      => null()
        type (T_Size2D), dimension(:), pointer  :: Slaves_Inner     => null()
        type (T_Size2D), dimension(:), pointer  :: Slaves_Mapping   => null()
        type (T_Size2D), dimension(:), pointer  :: Slaves_HaloMap   => null()
        integer                                 :: MPI_ID           = null_int
        type (T_Size2D)                         :: Global
        type (T_Size2D)                         :: Mapping
        type (T_Size2D)                         :: Inner
        type (T_Size2D)                         :: HaloMap
        integer                                 :: NeighbourSouth = null_int
        integer                                 :: NeighbourWest  = null_int
        integer                                 :: NeighbourEast  = null_int
        integer                                 :: NeighbourNorth = null_int
        integer                                 :: Halo_Points    = null_int
        integer, pointer, dimension(:,:)        :: Interfaces     => null()
        integer                                 :: NInterfaces    = null_int
        logical, dimension(1:4)                 :: OpenBordersON  = .true.
        character(PathLength)                   :: FilesListName  = "DecomposedFiles.dat"
        logical                                 :: FilesListOpen  = .false.
        integer                                 :: FilesListID    = null_int
        character(PathLength)                   :: ModelPath      = null_str
        type (T_Coef2D)                         :: Coef2D
        type (T_Coef3D)                         :: Coef3D
        logical                                 :: AutomaticLines = .false.
    end type T_DDecomp


    type T_HorizontalGrid
        integer                                 :: InstanceID = null_int !initialization: jauch - or should be set to 0 (zero)?

        !Former Bathymetry
        real, pointer, dimension(:  )           :: XX         => null()
        real, pointer, dimension(:  )           :: YY         => null()
        real                                    :: Xorig      = null_real
        real                                    :: Yorig      = null_real
        real                                    :: Latitude   = null_real
        real                                    :: Longitude  = null_real
        real                                    :: GRID_ANGLE = null_real

        integer                                 :: ProjType   = null_int
        real                                    :: SP1        = null_real
        real                                    :: SP2        = null_real
        real                                    :: Easting    = null_real !initialization: jauch
        real                                    :: Northing   = null_real !initialization: jauch
        integer                                 :: Datum      = null_int  !ellipsoid
        logical                                 :: UseLambert = .FALSE.
        integer                                 :: CoordType  = null_int
        integer                                 :: ZoneLong   = null_int
        integer                                 :: ZoneLat    = null_int
        integer, dimension(2)                   :: Grid_Zone

        integer, dimension(:, :), allocatable   :: Connections_U
        integer, dimension(:, :), allocatable   :: Connections_V
        integer, dimension(:, :), allocatable   :: Connections_Z
        real, pointer, dimension(:)             :: IWD_Distances_U   => null()
        real, pointer, dimension(:)             :: IWD_Distances_V   => null()
        real, dimension   (:),    allocatable   :: IWD_Distances_Z
        logical                                 :: UsedIWD_2Way      = .false.
        integer                                 :: IWD_Nodes_Z       = null_int
        integer                                 :: IWD_Nodes_U       = null_int
        integer                                 :: IWD_Nodes_V       = null_int

        type(T_Compute)                         :: Compute

        type (T_FatherGrid),     pointer        :: FirstFatherGrid => null()
        type (T_FatherGrid),     pointer        :: LastFatherGrid  => null()

        !Grid boundary
        type(T_Border),          pointer        :: GridBorderCart      => null()
        type(T_Border),          pointer        :: GridBorderCoord     => null()
        type(T_Border),          pointer        :: GridOutBorderCoord  => null()
        type(T_Border),          pointer        :: GridBorderAlongGrid => null()
        type(T_Border),          pointer        :: GridOutBorderCart   => null()


        !Distances (DXX... DVY)
        real, dimension(:, :), pointer          :: DXX           => null()
        real, dimension(:, :), pointer          :: DYY           => null()
        real, dimension(:, :), pointer          :: DZX           => null()
        real, dimension(:, :), pointer          :: DZY           => null()
        real, dimension(:, :), pointer          :: DUX           => null()
        real, dimension(:, :), pointer          :: DUY           => null()
        real, dimension(:, :), pointer          :: DVX           => null()
        real, dimension(:, :), pointer          :: DVY           => null()
        real, dimension(:, :), pointer          :: XX_IE         => null()
        real, dimension(:, :), pointer          :: YY_IE         => null()
        real, dimension(:, :), pointer          :: XX_AlongGrid  => null()
        real, dimension(:, :), pointer          :: YY_AlongGrid  => null()
        real, dimension(:, :), pointer          :: LatitudeConn  => null()
        real, dimension(:, :), pointer          :: LongitudeConn => null()

        logical                                 :: ReadCartCorners = .false.

        real, dimension(:, :), pointer          :: RotationX => null()
        real, dimension(:, :), pointer          :: RotationY => null()

        logical                                 :: CornersXYInput  = .false.
        logical                                 :: Distortion      = .false.
        logical                                 :: RegularRotation = .false.

        integer, dimension(:,:), pointer        :: DefineCellsMap  => null()
        integer, dimension(:,:), pointer        :: DefineFacesUMap => null()
        integer, dimension(:,:), pointer        :: DefineFacesVMap => null()
        integer, dimension(:,:), pointer        :: DefineCrossMap  => null()

        logical                                 :: NotDefinedCells = .false.

        !Latitude, Longitude
        real, dimension(:, :), pointer          :: LatitudeZ  => null()
        real, dimension(:, :), pointer          :: LongitudeZ => null()

        !Coriolis Factor
        real, dimension(:, :), pointer          :: F => null()

        !Area
        real, dimension(:, :), pointer          :: GridCellArea => null()

        !1D aux
        real, dimension(:   ), pointer          :: XX1D_Aux => null()
        real, dimension(:   ), pointer          :: YY1D_Aux => null()

        !Other
        type (T_Size2D)                         :: Size
        type (T_Size2D)                         :: WorkSize
        type (T_Size2D)                         :: GlobalWorkSize

        character(PathLength)                   :: FileName = null_str

        type(T_DDecomp)                         :: DDecomp
        type (T_BorderLimits)                   :: BorderLimits

        !Instances
        integer                                 :: ObjHDF5       = 0
        integer                                 :: ObjEnterData  = 0
        integer                                 :: ObjEnterData2 = 0

        type(T_Polygon),          pointer       :: AuxPolygon

        type (T_HorizontalGrid), pointer        :: GlobalGrid

        type (T_HorizontalGrid), pointer        :: Next => null()


    end type T_HorizontalGrid

    !Global Module Variables
    type (T_HorizontalGrid), pointer            :: FirstHorizontalGrid  => null()
    type (T_HorizontalGrid), pointer            :: Me                   => null()


    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructHorizontalGridV1(HorizontalGridID, DataFile,                    &
                                         MPI_ID, MasterID, LastSlaveID, ModelPath, STAT)

        !Arguments-------------------------------------------------------------
        integer                                 :: HorizontalGridID
        character(len=*), optional              :: DataFile
        integer, optional,  intent(IN)          :: MPI_ID
        integer, optional,  intent(IN)          :: MasterID
        integer, optional,  intent(IN)          :: LastSlaveID
        character(len=*), optional,  intent(IN) :: ModelPath
        integer, optional,  intent(OUT)         :: STAT

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, ready_
        integer                                 :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mHorizontalGrid_)) then
            nullify (FirstHorizontalGrid)
            call RegisterModule (mHorizontalGrid_)
        endif

        call Ready(HorizontalGridID, ready_)

cd2 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            nullify (Me%FirstFatherGrid)
            nullify (Me%LastFatherGrid )

            !Reads data file
            if (.not. present(DataFile)) then
                call ReadFileName('IN_BATIM', Me%FileName, Message = "Grid File", STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHorizontalGridV1 - HorizontalGrid - ERR01'
            else
                Me%FileName = DataFile
            endif

            if (present(MPI_ID)) then
                Me%DDecomp%MPI_ID = MPI_ID
            endif


            if (present(MasterID)) then
                Me%DDecomp%Master_MPI_ID = MasterID
                if (MasterID == MPI_ID) then
                    Me%DDecomp%Master = .true.
                endif
            endif

            if (present(LastSlaveID)) then
                if (MPI_ID > null_int .and. MasterID > null_int) then
                    if (LastSlaveID > null_int) then
                        Me%DDecomp%ON            = .true.
                        Me%DDecomp%MasterOrSlave = .true.
                        Me%DDecomp%Nslaves       = LastSlaveID - MasterID
                        if (present(ModelPath)) then
                            Me%DDecomp%ModelPath  = ModelPath
                        else
                            Me%DDecomp%ModelPath  = null_str
                        endif
                    endif
                endif
            endif



            !Construct the variable common to all module
            call ConstructGlobalVariables

            call GenerateGrid

            !Returns ID
            HorizontalGridID    = Me%InstanceID

            STAT_ = SUCCESS_
        else

            stop 'HorizontalGrid - ConstructHorizontalGridV1 - ERR99'

        end if cd2


        if (present(STAT)) STAT = STAT_

    end subroutine ConstructHorizontalGridV1

    !--------------------------------------------------------------------------
    subroutine ConstructHorizontalGridV2(HorizontalGridID, LatitudeConn, LongitudeConn, &
                                         XX, YY, Xorig, Yorig, Latitude, Longitude,     &
                                         ILB, IUB, JLB, JUB, STAT)


        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real, dimension(:, :), pointer, optional    :: LatitudeConn, LongitudeConn
        real, dimension(:   ), pointer              :: XX, YY
        real                          , optional    :: Xorig, Yorig
        real                                        :: Latitude, Longitude
        integer                                     :: ILB, IUB, JLB, JUB
        integer, optional,  intent(OUT)             :: STAT

        !Local-----------------------------------------------------------------
        real, dimension(:, :), pointer              :: LatitudeConn_, LongitudeConn_
        real                                        :: Xorig_, Yorig_
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mHorizontalGrid_)) then
            nullify (FirstHorizontalGrid)
            call RegisterModule (mHorizontalGrid_)
        endif

        call Ready(HorizontalGridID, ready_)

cd2 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            nullify (Me%FirstFatherGrid)
            nullify (Me%LastFatherGrid )

            if (present(LongitudeConn)) then
                LongitudeConn_ => LongitudeConn
            else
                nullify(LongitudeConn_)
            endif

            if (present(LatitudeConn)) then
                LatitudeConn_ => LatitudeConn
            else
                nullify(LatitudeConn_)
            endif

            if (present(Xorig)) then
                Xorig_ = Xorig
            else
                Xorig_ = 0.
            endif

            if (present(Yorig)) then
                Yorig_ = Yorig
            else
                Yorig_ = 0.
            endif

            call ConstructGlobalVariablesV1(LatitudeConn_, LongitudeConn_, XX, YY,      &
                                            Xorig_, Yorig_, Latitude, Longitude,        &
                                            ILB, IUB, JLB, JUB)

            call GenerateGrid

            !Returns ID
            HorizontalGridID    = Me%InstanceID

            STAT_ = SUCCESS_
        else

            stop 'HorizontalGrid - ConstructHorizontalGridV2 - ERR99'

        end if cd2


        if (present(STAT)) STAT = STAT_

    end subroutine ConstructHorizontalGridV2

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine ConstructHorizontalGridV3(HorizontalGridID, HDF5ID, STAT)


        !Arguments-------------------------------------------------------------
        integer                                 :: HorizontalGridID
        integer                                 :: HDF5ID
        integer, optional,  intent(OUT)         :: STAT

        !Local-----------------------------------------------------------------
        integer                             :: STAT_, ready_, nUsers

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mHorizontalGrid_)) then
            nullify (FirstHorizontalGrid)
            call RegisterModule (mHorizontalGrid_)
        endif

        call Ready(HorizontalGridID, ready_)

cd2 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            !Associates Horizontal Grid
            Me%ObjHDF5 = AssociateInstance (mHDF5_, HDF5ID)

            nullify (Me%FirstFatherGrid)
            nullify (Me%LastFatherGrid )

            call ConstructGlobalVariablesV2()

            !Associates Horizontal Grid
            nUsers     = DeassociateInstance (mHDF5_, HDF5ID)
            if (nUsers == 0) stop 'HorizontalGrid - ConstructHorizontalGridV3 - ERR10'


            call GenerateGrid

            !Returns ID
            HorizontalGridID    = Me%InstanceID

            STAT_ = SUCCESS_
        else

            stop 'HorizontalGrid - ConstructHorizontalGridV3 - ERR20'

        end if cd2


        if (present(STAT)) STAT = STAT_

    end subroutine ConstructHorizontalGridV3

    !--------------------------------------------------------------------------


    subroutine GenerateGrid

           !Constructs XX_IE, YY_IE, LatitudeZ and LongitudeZ
            call Mercator

            if (.not. Me%CornersXYInput) then

                !Constructs XX_Z, YY_Z, XX_U, YY_U, XX_V, YY_V, XX_Cross, YY_Cross
                call GridPointsLocation1D

            endif

            !Constructs XX2D_Z, YY2D_Z, XX2D_U, YY2D_U, XX2D_V, YY2D_V, XX2D_Cross, YY2D_Cross
            call GridPointsLocation2D


            !Check grid border type
            call CheckGridBorder

            !Computes DXX, DYY, DZX, DZY, DUX, DUY, DVX, DVY, XX_AlongGrid, YY_AlongGrid
            call ComputeDistances

            !Computes RotationX, RotationY
            call ComputeRotation

            !Computes Area
            call ComputeGridCellArea

            !Computes Coriolis
            call ComputeCoriolis

            !Defines the grid border polygon
            call DefineBorderPolygons

            !Intialization of domain decomposition procedure
            !call ConstructDDecomp

    end subroutine GenerateGrid


    subroutine ConstructDDecomp

#ifdef _USE_MPI

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len = PathLength)             :: DDFile
        character(len = StringLength)           :: Message
        integer                                 :: STAT_CALL, i, iflag
        integer                                 :: AuxHalo
        !type (T_Size2D)                         :: Mapping
        logical                                 :: Exist
        integer, pointer, dimension(:)          :: Aux1D
        integer                                 :: Source, Destination
        integer                                 :: iSize
        integer, save                           :: Precision
        integer                                 :: status(MPI_STATUS_SIZE)
        !----------------------------------------------------------------------

        Me%DDecomp%Global%JLB = Me%GlobalWorkSize%JLB
        Me%DDecomp%Global%JUB = Me%GlobalWorkSize%JUB
        Me%DDecomp%Global%ILB = Me%GlobalWorkSize%ILB
        Me%DDecomp%Global%IUB = Me%GlobalWorkSize%IUB

        allocate(Me%DDecomp%Slaves_MPI_ID(Me%DDecomp%Nslaves))

        do i=1, Me%DDecomp%Nslaves
            Me%DDecomp%Slaves_MPI_ID(i) = Me%DDecomp%Master_MPI_ID + i
        enddo

        ! ---> ASCII file used to define domain decomposition
        Message   ='ASCII file used to define domain decomposition.'
        Message   = trim(Message)

        !call ReadFileName('D_DECOMP', DDFile, Message = Message, STAT = STAT_CALL)

        call GetData(DDFile,Me%ObjEnterData, iflag,                                     &
                     keyword      = 'D_DECOMP',                                         &
                     default      = null_str,                                           &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
           call SetError(FATAL_, INTERNAL_, "ConstructDDecomp - Hydrodynamic - ERR10")

ii:     if (iflag == 0) then
            Me%DDecomp%Auto = .true.
        else ii
            inquire(FILE = DDFile, EXIST = Exist)
iE:         if  (Exist) then
                !open()
                call ConstructEnterData(Me%ObjEnterData2, DDFile, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                     &
                   call SetError(FATAL_, INTERNAL_, "ConstructDDecomp - Hydrodynamic - ERR20")

                write(*,*) "Read from file Domain Decomposition mapping"
                write(*,*) "Present MPI ID =", Me%DDecomp%MPI_ID
                write(*,*) "Master MPI_ID = ", Me%DDecomp%Master_MPI_ID
                write(*,*) "Number of domain Slaves = ",Me%DDecomp%Nslaves
                do i=1, Me%DDecomp%Nslaves
                    write(*,*) "ID of Slave number ",i, "is =", Me%DDecomp%Slaves_MPI_ID(i)
                enddo

                call OptionsDDecomp

                call KillEnterData(Me%ObjEnterData2, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                     &
                   call SetError(FATAL_, INTERNAL_, "ConstructDDecomp - Hydrodynamic - ERR30")

            else iE
                write(*,*) 'keyword D_DECOMP do not exit file =',trim(DDFile)
                write(*,*) 'Automatic decomposition will be assumed'
                Me%DDecomp%Auto = .true.
            endif iE

        endif ii

        if (Me%DDecomp%MasterOrSlave) then

            if (Me%DDecomp%Auto) then
                call AutomaticDDecomp
            endif

            if (Me%DDecomp%NeighbourSouth /= null_int) then
                AuxHalo = Me%DDecomp%Halo_Points
                !1 - South; 2 - North; 3 - West; 4 - East
                Me%DDecomp%OpenBordersON(1)  = .false.
            else
                AuxHalo = 0.
            endif

            Me%DDecomp%HaloMap%ILB  = Me%DDecomp%Mapping%ILB - AuxHalo
            Me%WorkSize%ILB         = 1
            Me%DDecomp%Inner%ILB    = Me%WorkSize%ILB        + AuxHalo

            if (Me%DDecomp%NeighbourNorth /= null_int) then
                AuxHalo = Me%DDecomp%Halo_Points
                !1 - South; 2 - North; 3 - West; 4 - East
                Me%DDecomp%OpenBordersON(2)  = .false.
            else
                AuxHalo = 0.
            endif

            Me%DDecomp%HaloMap%IUB  = Me%DDecomp%Mapping%IUB + AuxHalo
            Me%WorkSize%IUB         = Me%DDecomp%HaloMap%IUB - Me%DDecomp%HaloMap%ILB + 1
            Me%DDecomp%Inner%IUB    = Me%WorkSize%IUB        - AuxHalo


            if (Me%DDecomp%HaloMap%IUB - Me%DDecomp%HaloMap%ILB /=  &
                Me%WorkSize%IUB        - Me%WorkSize%ILB) then

                write(*,*) "Decomposition domains is inconsistent with the grid data input - Lines "
                write(*,*) 'WorkSize ILB,IUB, =',Me%WorkSize%ILB,Me%WorkSize%IUB
                write(*,*) 'HaloMap  ILB,IUB, =',Me%DDecomp%HaloMap%ILB,Me%DDecomp%HaloMap%IUB
                stop "ConstructDDecomp - ModuleHorizontalGrid - ERR40"

            endif

            if (Me%DDecomp%NeighbourWest /= null_int) then
                AuxHalo = Me%DDecomp%Halo_Points
                !1 - South; 2 - North; 3 - West; 4 - East
                Me%DDecomp%OpenBordersON(3)  = .false.
            else
                AuxHalo = 0.
            endif

            Me%DDecomp%HaloMap%JLB  = Me%DDecomp%Mapping%JLB - AuxHalo
            Me%WorkSize%JLB         = 1
            Me%DDecomp%Inner%JLB    = Me%WorkSize%JLB        + AuxHalo


            if (Me%DDecomp%NeighbourEast /= null_int) then
                AuxHalo = Me%DDecomp%Halo_Points
                !1 - South; 2 - North; 3 - West; 4 - East
                Me%DDecomp%OpenBordersON(4)  = .false.
            else
                AuxHalo = 0.
            endif

            Me%DDecomp%HaloMap%JUB  = Me%DDecomp%Mapping%JUB + AuxHalo
            Me%WorkSize%JUB         = Me%DDecomp%HaloMap%JUB - Me%DDecomp%HaloMap%JLB + 1
            Me%DDecomp%Inner%JUB    = Me%WorkSize%JUB        - AuxHalo

            if (Me%DDecomp%HaloMap%JUB - Me%DDecomp%HaloMap%JLB /=  &
                Me%WorkSize%JUB        - Me%WorkSize%JLB) then

                write(*,*) "Decomposition domains is inconsistent with the grid data input - Columns "
                write(*,*) 'WorkSize JLB,JUB, =',Me%WorkSize%JLB,Me%WorkSize%JUB
                write(*,*) 'HaloMap  JLB,JUB, =',Me%DDecomp%HaloMap%JLB,Me%DDecomp%HaloMap%JUB
                stop "ConstructDDecomp - ModuleHorizontalGrid - ERR50"

            endif


            if (Me%DDecomp%Master) then

                !if (Me%DDecomp%Nslaves /= 7) write(*,*) 'e44'

                !allocate(Me%DDecomp%Slaves_MPI_ID (Me%DDecomp%Nslaves))
                allocate(Me%DDecomp%Slaves_Inner  (Me%DDecomp%Nslaves))
                allocate(Me%DDecomp%Slaves_Size   (Me%DDecomp%Nslaves))
                allocate(Me%DDecomp%Slaves_Mapping(Me%DDecomp%Nslaves))
                allocate(Me%DDecomp%Slaves_HaloMap(Me%DDecomp%Nslaves))

            endif

            !PCL - Master slave mapping
            iSize = 16
            allocate(Aux1D(iSize))

            if (Me%DDecomp%Master) then

                do i=1, Me%DDecomp%Nslaves

                    Precision = MPI_INTEGER
                    Source    = Me%DDecomp%Slaves_MPI_ID(i)

                    write(*,*) 'mpi_receive from to ', i, Me%DDecomp%MPI_ID

                    call MPI_Recv (Aux1D(1:iSize), iSize, Precision,   Source, 30001, MPI_COMM_WORLD, status, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop "ConstructDDecomp - ModuleHorizontalGrid - ERR60"

                    Me%DDecomp%Slaves_Inner  (i)%ILB = Aux1D(1)
                    Me%DDecomp%Slaves_Inner  (i)%IUB = Aux1D(2)
                    Me%DDecomp%Slaves_Inner  (i)%JLB = Aux1D(3)
                    Me%DDecomp%Slaves_Inner  (i)%JUB = Aux1D(4)

                    Me%DDecomp%Slaves_HaloMap(i)%ILB = Aux1D(5)
                    Me%DDecomp%Slaves_HaloMap(i)%IUB = Aux1D(6)
                    Me%DDecomp%Slaves_HaloMap(i)%JLB = Aux1D(7)
                    Me%DDecomp%Slaves_HaloMap(i)%JUB = Aux1D(8)

                    Me%DDecomp%Slaves_Size   (i)%ILB = Aux1D(9)
                    Me%DDecomp%Slaves_Size   (i)%IUB = Aux1D(10)
                    Me%DDecomp%Slaves_Size   (i)%JLB = Aux1D(11)
                    Me%DDecomp%Slaves_Size   (i)%JUB = Aux1D(12)

                    Me%DDecomp%Slaves_Mapping(i)%ILB = Aux1D(13)
                    Me%DDecomp%Slaves_Mapping(i)%IUB = Aux1D(14)
                    Me%DDecomp%Slaves_Mapping(i)%JLB = Aux1D(15)
                    Me%DDecomp%Slaves_Mapping(i)%JUB = Aux1D(16)

                enddo

            else

                Aux1D(1) = Me%DDecomp%Inner%ILB
                Aux1D(2) = Me%DDecomp%Inner%IUB
                Aux1D(3) = Me%DDecomp%Inner%JLB
                Aux1D(4) = Me%DDecomp%Inner%JUB

                Aux1D(5) = Me%DDecomp%HaloMap%ILB
                Aux1D(6) = Me%DDecomp%HaloMap%IUB
                Aux1D(7) = Me%DDecomp%HaloMap%JLB
                Aux1D(8) = Me%DDecomp%HaloMap%JUB

                Aux1D(9) = Me%WorkSize%ILB
                Aux1D(10)= Me%WorkSize%IUB
                Aux1D(11)= Me%WorkSize%JLB
                Aux1D(12)= Me%WorkSize%JUB

                Aux1D(13) = Me%DDecomp%Mapping%ILB
                Aux1D(14) = Me%DDecomp%Mapping%IUB
                Aux1D(15) = Me%DDecomp%Mapping%JLB
                Aux1D(16) = Me%DDecomp%Mapping%JUB

                Precision   = MPI_INTEGER
                Destination = Me%DDecomp%Master_MPI_ID

                write(*,*) 'mpi_send', Me%DDecomp%MPI_ID

                call MPI_Send (Aux1D(1:iSize), iSize, Precision, Destination, 30001, MPI_COMM_WORLD, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop "ConstructDDecomp - ModuleHorizontalGrid - ERR80"

            endif

            deallocate(Aux1d)


        endif

#endif _USE_MPI

    end subroutine ConstructDDecomp

    !End----------------------------------------------------------------

#if _USE_MPI

    subroutine OptionsDDecomp()

        !Arguments------------------------------------------------------------

        !Local----------------------------------------------------------------
        character(len = StringLength), parameter    :: BeginBlock1 ="<BeginSubDD>"
        character(len = StringLength), parameter    :: EndBlock1   ="<EndSubDD>"
        character(len = StringLength), parameter    :: BeginBlock2 ="<BeginInterfaceSN>"
        character(len = StringLength), parameter    :: EndBlock2   ="<EndInterfaceSN>"
        character(len = StringLength), parameter    :: BeginBlock3 ="<BeginInterfaceWE>"
        character(len = StringLength), parameter    :: EndBlock3   ="<EndInterfaceWE>"


        integer, dimension(:), allocatable          :: Aux1D
        logical                                     :: BlockFound
        integer                                     :: ClientNumber, STAT_CALL, iflag
        integer                                     :: in, line, FirstLine, LastLine, i, ii, jj
        integer                                     :: SN_N_Interfaces, WE_N_Interfaces, MPI_ID
        logical                                     :: MissMatchID

        !Begin----------------------------------------------------------------



        call RewindBuffer(Me%ObjEnterData2, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR120'


        Me%DDecomp%Auto = .false.

iSl:    do i =1, Me%DDecomp%Nslaves + 1

            !Searches sub-domains blocks
            call ExtractBlockFromBuffer (Me%ObjEnterData2, ClientNumber,                &
                                         BeginBlock1, EndBlock1,                        &
                                         BlockFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR130'
            if (.not. BlockFound     ) then
                Me%DDecomp%Auto = .true.
                write(*,*) 'OptionsDDecomp  - ModuleHorizontalGrid - WRN140'
                write(*,*) 'Domain Decomposition in automatic mode'
                exit
            endif

            call GetData(Value          = MPI_ID,                                       &
                         EnterDataID    = Me%ObjEnterData2,                             &
                         flag           = iflag,                                        &
                         keyword        = 'MPI_ID',                                     &
                         SearchType     = FromBlock,                                    &
                         ClientModule   = 'ModuleHorizontalGrid',                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR150'
            if (iflag     == 0       ) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR160'

            MissMatchID = .true.

            do ii =1, Me%DDecomp%Nslaves
                if (Me%DDecomp%Slaves_MPI_ID(ii) == MPI_ID) then
                    MissMatchID = .false.
                endif
            enddo

            if (Me%DDecomp%Master_MPI_ID == MPI_ID) then
                MissMatchID = .false.
            endif

            if (MissMatchID) then
                write(*,*) 'MPI ID of present Domain', Me%DDecomp%MPI_ID
                write(*,*) 'Domain -', MPI_ID, ' is not one of decomposition domains'
                write(*,*) "All MPI_ID - Slaves"
                do ii =1, Me%DDecomp%Nslaves
                    write(*,*) "MPI_ID =", Me%DDecomp%Slaves_MPI_ID(ii)
                enddo
                write(*,*) "MPI_ID Master=",Me%DDecomp%Master_MPI_ID
                stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR165'
            endif

            if (MPI_ID == Me%DDecomp%MPI_ID) then

                Me%DDecomp%MasterOrSlave = .true.

                allocate(Aux1D(1:2))

                call GetData(Vector         = Aux1D,                                        &
                             EnterDataID    = Me%ObjEnterData2,                             &
                             flag           = iflag,                                        &
                             keyword        = 'ILB_IUB',                                    &
                             SearchType     = FromBlock,                                    &
                             ClientModule   = 'ModuleHorizontalGrid',                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= KEYWORD_NOT_FOUND_ERR_) then
                    stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR170'
                endif
                if (iflag     /= 2) then
                    if (iflag == 0) then
                        Aux1D(1) = Me%DDecomp%Global%ILB
                        Aux1D(2) = Me%DDecomp%Global%IUB
                    else
                        stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR180'
                    endif
                endif

                Me%DDecomp%Mapping%ILB = Aux1D(1)
                Me%DDecomp%Mapping%IUB = Aux1D(2)

                call GetData(Vector         = Aux1D,                                        &
                             EnterDataID    = Me%ObjEnterData2,                             &
                             flag           = iflag,                                        &
                             keyword        = 'JLB_JUB',                                    &
                             SearchType     = FromBlock,                                    &
                             ClientModule   = 'ModuleHorizontalGrid',                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= KEYWORD_NOT_FOUND_ERR_) then
                    stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR190'
                endif
                if (iflag     /= 2) then
                    if (iflag == 0) then
                        Aux1D(1) = Me%DDecomp%Global%JLB
                        Aux1D(2) = Me%DDecomp%Global%JUB
                    else
                        stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR200'
                    endif
                endif

                Me%DDecomp%Mapping%JLB = Aux1D(1)
                Me%DDecomp%Mapping%JUB = Aux1D(2)


                deallocate(Aux1D)

                exit

            endif

        enddo iSl

        call Block_Unlock(Me%ObjEnterData2, ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR210'

iAuto:  if (Me%DDecomp%Auto) then

            call GetData(Value          = Me%DDecomp%AutomaticLines,                    &
                         EnterDataID    = Me%ObjEnterData2,                             &
                         flag           = iflag,                                        &
                         keyword        = 'AUTOMATIC_LINES',                            &
                         SearchType     = FromFile,                                     &
                         default        = .false.,                                      &
                         ClientModule   = 'ModuleHorizontalGrid',                       &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR230'

        else

            call GetData(Value          = Me%DDecomp%NInterfaces,                       &
                         EnterDataID    = Me%ObjEnterData2,                                 &
                         flag           = iflag,                                            &
                         keyword        = 'INTERFACES_NUMBER',                              &
                         SearchType     = FromFile,                                         &
                         ClientModule   = 'ModuleHorizontalGrid',                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR230'

            allocate(Me%DDecomp%Interfaces(Me%DDecomp%NInterfaces,3))

            allocate(Aux1D(1:2))

            call RewindBuffer(Me%ObjEnterData2, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR220'

            !Searches sub-domains blocks
            call ExtractBlockFromBuffer (Me%ObjEnterData2, ClientNumber,                    &
                                         BeginBlock2, EndBlock2, BlockFound,                &
                                         FirstLine = FirstLine, LastLine = LastLine,        &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR240'

            if (.not. BlockFound) then
                SN_N_Interfaces = 0
            endif

            in = 0

            if (BlockFound) then

                SN_N_Interfaces = LastLine - FirstLine - 1

                if (LastLine == FirstLine) then
                    SN_N_Interfaces = 0
                endif

                if (SN_N_Interfaces > Me%DDecomp%NInterfaces) then
                    stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR250'
                endif

                do line = FirstLine + 1, LastLine - 1

                    call GetData(Vector         = Aux1D,                                            &
                                 EnterDataID    = Me%ObjEnterData2,                                 &
                                 flag           = iflag,                                            &
                                 SearchType     = FromBlock,                                        &
                                 ClientModule   = 'ModuleHorizontalGrid',                             &
                                 Buffer_Line    = line,                                             &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR260'
                    if (iflag     /= 2       ) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR270'

                    in = in + 1

                    do jj = 1, 2
                        MissMatchID = .true.

                        do ii =1, Me%DDecomp%Nslaves
                            if (Me%DDecomp%Slaves_MPI_ID(ii) == Aux1D(jj)) then
                                MissMatchID = .false.
                            endif
                        enddo

                        if (Me%DDecomp%Master_MPI_ID == Aux1D(jj)) then
                            MissMatchID = .false.
                        endif

                        if (MissMatchID) then
                            write(*,*) 'Domain -', Aux1D(jj), ' is not one of decomposition domains'
                            stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR275'
                        endif
                    enddo


                    Me%DDecomp%Interfaces(in,1) = Aux1D(1)
                    Me%DDecomp%Interfaces(in,2) = Aux1D(2)
                    Me%DDecomp%Interfaces(in,3) = SouthNorth_

                    if (Me%DDecomp%MPI_ID == Aux1D(1)) then
                        Me%DDecomp%NeighbourNorth = Aux1D(2)
                    endif

                    if (Me%DDecomp%MPI_ID == Aux1D(2)) then
                        Me%DDecomp%NeighbourSouth = Aux1D(1)
                    endif

                enddo

             endif

            call RewindBuffer(Me%ObjEnterData2, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR275'

            !Searches sub-domains blocks
            call ExtractBlockFromBuffer (Me%ObjEnterData2, ClientNumber,                    &
                                         BeginBlock3, EndBlock3, BlockFound,                &
                                         FirstLine = FirstLine, LastLine = LastLine,        &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR280'

            if (.not. BlockFound) then
                WE_N_Interfaces = 0
            endif

            if (BlockFound) then

                WE_N_Interfaces = LastLine - FirstLine - 1

                if (LastLine == FirstLine) then
                    WE_N_Interfaces = 0
                endif

                if (SN_N_Interfaces + WE_N_Interfaces /= Me%DDecomp%NInterfaces) then
                    stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR290'
                endif


                do line = FirstLine + 1, LastLine - 1

                    call GetData(Vector         = Aux1D,                                            &
                                 EnterDataID    = Me%ObjEnterData2,                                 &
                                 flag           = iflag,                                            &
                                 SearchType     = FromBlock,                                        &
                                 ClientModule   = 'ModuleHorizontalGrid',                             &
                                 Buffer_Line    = line,                                             &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR300'
                    if (iflag     /= 2       ) stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR310'

                    in = in + 1

                    do jj = 1, 2
                        MissMatchID = .true.

                        do ii =1, Me%DDecomp%Nslaves
                            if (Me%DDecomp%Slaves_MPI_ID(ii) == Aux1D(jj)) then
                                MissMatchID = .false.
                            endif
                        enddo

                        if (Me%DDecomp%Master_MPI_ID == Aux1D(jj)) then
                            MissMatchID = .false.
                        endif

                        if (MissMatchID) then
                            write(*,*) 'Domain -', Aux1D(jj), ' is not one of decomposition domains'
                            stop 'OptionsDDecomp  - ModuleHorizontalGrid - ERR315'
                        endif
                    enddo

                    Me%DDecomp%Interfaces(in,1) = Aux1D(1)
                    Me%DDecomp%Interfaces(in,2) = Aux1D(2)
                    Me%DDecomp%Interfaces(in,3) = WestEast_

                    if (Me%DDecomp%MPI_ID == Aux1D(1)) then
                        Me%DDecomp%NeighbourEast = Aux1D(2)
                    endif

                    if (Me%DDecomp%MPI_ID == Aux1D(2)) then
                        Me%DDecomp%NeighbourWest = Aux1D(1)
                    endif

                enddo

            endif

            deallocate(Aux1D)

        endif iAuto

    end subroutine OptionsDDecomp

    !End------------------------------------------------------------------


    subroutine AutomaticDDecomp()

        !Arguments------------------------------------------------------------

        !Local----------------------------------------------------------------
        !Begin----------------------------------------------------------------


        write(*,*) 'halo_points', Me%DDecomp%Halo_Points

        if (Me%DDecomp%Global%IUB  > Me%DDecomp%Global%JUB .or. Me%DDecomp%AutomaticLines) then
            call AutomaticDDecompLines  ()
        else
            call AutomaticDDecompColumns()
        endif



    end subroutine AutomaticDDecomp

    !End------------------------------------------------------------------

    subroutine AutomaticDDecompLines()

        !Arguments------------------------------------------------------------

        !Local----------------------------------------------------------------


        integer                                     :: i, iSouth, iNorth
        integer                                     :: iub_map, ilb_map, ILB, IUB, ND, is
        real                                        :: IDD

        !Begin----------------------------------------------------------------


        !In automatic mode MOHID slices the domains along the matrix lines (ILB-IUB)

        Me%DDecomp%Mapping%JLB = Me%DDecomp%Global%JLB
        Me%DDecomp%Mapping%JUB = Me%DDecomp%Global%JUB

        ILB  = Me%DDecomp%Global%ILB
        IUB  = Me%DDecomp%Global%IUB
        ND   = Me%DDecomp%Nslaves + 1
        IDD  = real (IUB - ILB + 1) / real(ND)
        is   = Me%DDecomp%Master_MPI_ID

        do i=1, Me%DDecomp%Nslaves + 1
            if (is + i - 1 == Me%DDecomp%MPI_ID) then
                ilb_map = ILB     + int(real(i-1)*IDD)
                iub_map = ILB - 1 + int(real(i  )*IDD)
                exit
            endif
        enddo

        Me%DDecomp%Mapping%ILB = ilb_map
        Me%DDecomp%Mapping%IUB = iub_map
        write(*,*) 'Limits ', Me%DDecomp%MPI_ID, iub_map, ilb_map

        Me%DDecomp%NInterfaces = Me%DDecomp%Nslaves

        allocate(Me%DDecomp%Interfaces(Me%DDecomp%NInterfaces,3))

        do i = 1, Me%DDecomp%NInterfaces

            iSouth = Me%DDecomp%Master_MPI_ID + i - 1
            iNorth = Me%DDecomp%Master_MPI_ID + i

            Me%DDecomp%Interfaces(i,1) = iSouth
            Me%DDecomp%Interfaces(i,2) = iNorth
            Me%DDecomp%Interfaces(i,3) = SouthNorth_
            write(*,*) 'Interface ', i, Me%DDecomp%Interfaces(i,1), Me%DDecomp%Interfaces(i,2)

             if (Me%DDecomp%MPI_ID == iSouth) then
                Me%DDecomp%NeighbourNorth = iNorth
            endif

            if (Me%DDecomp%MPI_ID == iNorth) then
                Me%DDecomp%NeighbourSouth = iSouth
            endif

        enddo


    end subroutine AutomaticDDecompLines

    !End------------------------------------------------------------------


    subroutine AutomaticDDecompColumns()

        !Arguments------------------------------------------------------------

        !Local----------------------------------------------------------------


        integer                                     :: i, jWest, jEast
        integer                                     :: jub_map, jlb_map, JLB, JUB, ND, is
        real                                        :: IDD

        !Begin----------------------------------------------------------------


        !In automatic mode MOHID slices the domains along the matrix columns (JLB-JUB)

        Me%DDecomp%Mapping%ILB = Me%DDecomp%Global%ILB
        Me%DDecomp%Mapping%IUB = Me%DDecomp%Global%IUB

        JLB  = Me%DDecomp%Global%JLB
        JUB  = Me%DDecomp%Global%JUB
        ND   = Me%DDecomp%Nslaves + 1
        IDD  = real (JUB - JLB + 1) / real(ND)
        is   = Me%DDecomp%Master_MPI_ID

        do i=1, Me%DDecomp%Nslaves + 1
            if (is + i - 1 == Me%DDecomp%MPI_ID) then
                jlb_map = JLB     + int(real(i-1)*IDD)
                jub_map = JLB - 1 + int(real(i  )*IDD)
                exit
            endif
        enddo

        Me%DDecomp%Mapping%JLB = jlb_map
        Me%DDecomp%Mapping%JUB = jub_map
        write(*,*) 'Limits ', Me%DDecomp%MPI_ID, jub_map, jlb_map

        Me%DDecomp%NInterfaces = Me%DDecomp%Nslaves

        allocate(Me%DDecomp%Interfaces(Me%DDecomp%NInterfaces,3))

        do i = 1, Me%DDecomp%NInterfaces

            jWest = Me%DDecomp%Master_MPI_ID + i - 1
            jEast = Me%DDecomp%Master_MPI_ID + i

            Me%DDecomp%Interfaces(i,1) = jWest
            Me%DDecomp%Interfaces(i,2) = jEast
            Me%DDecomp%Interfaces(i,3) = WestEast_
            write(*,*) 'Interface ', i, Me%DDecomp%Interfaces(i,1), Me%DDecomp%Interfaces(i,2)

             if (Me%DDecomp%MPI_ID == jWest) then
                Me%DDecomp%NeighbourEast = jEast
            endif

            if (Me%DDecomp%MPI_ID == jEast) then
                Me%DDecomp%NeighbourWest = jWest
            endif

        enddo


    end subroutine AutomaticDDecompColumns

    !End------------------------------------

#endif _USE_MPI


    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_HorizontalGrid), pointer            :: NewObjHorizontalGrid
        type (T_HorizontalGrid), pointer            :: PreviousObjHorizontalGrid


        !Allocates new instance
        allocate (NewObjHorizontalGrid)
        nullify  (NewObjHorizontalGrid%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstHorizontalGrid)) then
            FirstHorizontalGrid     => NewObjHorizontalGrid
            Me                      => NewObjHorizontalGrid
        else
            PreviousObjHorizontalGrid   => FirstHorizontalGrid
            Me                          => FirstHorizontalGrid%Next
            do while (associated(Me))
                PreviousObjHorizontalGrid  => Me
                Me                         => Me%Next
            enddo
            Me                             => NewObjHorizontalGrid
            PreviousObjHorizontalGrid%Next => NewObjHorizontalGrid
        endif

        Me%InstanceID = RegisterNewInstance (mHORIZONTALGRID_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ConstructFatherGridLocation(HorizontalGridID, HorizontalGridFatherID, &
                                           GridID, OkCross, OkZ, OkU, OkV, Window, STAT)

        !Arguments-------------------------------------------------------------
        integer                                   :: HorizontalGridID
        integer                                   :: HorizontalGridFatherID
        integer, optional,           intent (IN)  :: GridID
        logical, optional,           intent (IN)  :: OkCross, OkZ, OkU, OkV
        type (T_Size2D), optional                 :: Window
        integer, optional,           intent (OUT) :: STAT

        !Local-----------------------------------------------------------------
        type (T_HorizontalGrid), pointer          :: ObjHorizontalGridFather
        type (T_FatherGrid), pointer              :: NewFatherGrid
        integer                                   :: STAT_, ready_, GridID_
        logical                                   :: OkZ_, OkU_, OkV_, OkCross_
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            if (present(GridID)) then
                GridID_ = GridID
            else
                GridID_ = LargeScaleModel_
            endif


            if (present(OkCross)) then
                OkCross_ = OkCross
            else
                OkCross_ = .false.
            endif

            if (present(OkZ)) then
                OkZ_ = OkZ
            else
                OkZ_ = .true.
            endif

            if (present(OkU)) then
                OkU_ = OkU
            else
                OkU_ = .true.
            endif


            if (present(OkV)) then
                OkV_ = OkV
            else
                OkV_ = .true.
            endif

            allocate(NewFatherGrid)

            nullify(ObjHorizontalGridFather)
            call LocateObjFather        (ObjHorizontalGridFather, HorizontalGridFatherID)

            if (present(Window)) then
                NewFatherGrid%MPI_Window = Window
            else
                NewFatherGrid%MPI_Window = ObjHorizontalGridFather%WorkSize
            endif

            if (ObjHorizontalGridFather%CornersXYInput) then

                call ConstructNewFatherGrid2D (ObjHorizontalGridFather, NewFatherGrid, GridID_, &
                                             OkZ_, OkU_, OkV_, OkCross_)
            else
                call ConstructNewFatherGrid1D (ObjHorizontalGridFather, NewFatherGrid, GridID_, &
                                             OkZ_, OkU_, OkV_, OkCross_)
            endif

            if (Me%CoordType == SIMPLE_GEOG_ .and. .not. Me%ReadCartCorners .and. Me%ProjType == PAULO_PROJECTION_) then
                if (ObjHorizontalGridFather%Longitude /= Me%Longitude) then
                    write(*,*) 'LONGITUDE (in the bathymetry file) must be the same in Father and Son Domains'
                    stop 'ConstructFatherGridLocation - ModuleHorizontalGrid - ERR10'
                endif
                if (ObjHorizontalGridFather%Latitude /= Me%Latitude) then
                    write(*,*) 'LATITUDE (in the bathymetry file) must be the same in Father and Son Domains'
                    stop 'ConstructFatherGridLocation - ModuleHorizontalGrid - ERR20'
                endif

            endif


            call Add_FatherGrid         (NewFatherGrid)

            STAT_ = SUCCESS_

        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_


    end subroutine ConstructFatherGridLocation




    !--------------------------------------------------------------------------

    subroutine GetSonWindow(HorizontalGridID, HorizontalGridFatherID, Window, STAT)

        !Arguments-------------------------------------------------------------
        integer                                   :: HorizontalGridID
        integer                                   :: HorizontalGridFatherID
        type (T_Size2D),             intent (OUT) :: Window
        integer, optional,           intent (OUT) :: STAT

        !Local-----------------------------------------------------------------
        type (T_HorizontalGrid), pointer          :: ObjHorizontalGridFather
        integer                                   :: STAT_, ready_
        real,   dimension(:,:), pointer           :: XX2D, YY2D
        real,      dimension(:,:), pointer        :: WindowInXY
        integer,   dimension(:,:), pointer        :: WindowOutJI
        logical                                   :: WindowWithData
        real                                      :: West, East, South, North
        integer                                   :: ILB, IUB, JLB, JUB

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            nullify(ObjHorizontalGridFather)
            call LocateObjFather        (ObjHorizontalGridFather, HorizontalGridFatherID)

            allocate (WindowInXY (1:2,1:2))
            allocate (WindowOutJI(1:2,1:2))

            WindowOutJI(:,:) = null_int
            WindowInXY (:,:) = null_real

            West    = Me%GridOutBorderCoord%Polygon_%Limits%Left
            East    = Me%GridOutBorderCoord%Polygon_%Limits%Right
            South   = Me%GridOutBorderCoord%Polygon_%Limits%Bottom
            North   = Me%GridOutBorderCoord%Polygon_%Limits%Top

            WindowInXY(2,1) = South
            WindowInXY(2,2) = North
            WindowInXY(1,1) = West
            WindowInXY(1,2) = East

            ILB     = ObjHorizontalGridFather%WorkSize%ILB
            IUB     = ObjHorizontalGridFather%WorkSize%IUB
            JLB     = ObjHorizontalGridFather%WorkSize%JLB
            JUB     = ObjHorizontalGridFather%WorkSize%JUB

            if (Me%CoordType == SIMPLE_GEOG_ .or. Me%CoordType == GEOG_) then

                XX2D        => ObjHorizontalGridFather%LongitudeConn
                YY2D        => ObjHorizontalGridFather%LatitudeConn

            else

                XX2D        => ObjHorizontalGridFather%XX_IE
                YY2D        => ObjHorizontalGridFather%YY_IE

            endif

            call ArrayPolygonWindow(XX              = XX2D,                             &
                                    YY              = YY2D,                             &
                                    WIn             = WindowInXY,                       &
                                    ILB             = ILB,                              &
                                    IUB             = IUB+1,                            &
                                    JLB             = JLB,                              &
                                    JUB             = JUB+1,                            &
                                    WOut            = WindowOutJI,                      &
                                    WindowWithData  = WindowWithData)

            if (.not.WindowWithData) then
                write(*,*) 'Father grid do not intersect Nested model'
                stop       'GetSonWindow - ModuleHoriuzontalGrid - ERR20'
            endif

            Window%ILB = max(WindowOutJI(1,1) - 1,ILB)
            Window%IUB = min(WindowOutJI(1,2) + 1,IUB)
            Window%JLB = max(WindowOutJI(2,1) - 1,JLB)
            Window%JUB = min(WindowOutJI(2,2) + 1,JUB)

            deallocate (WindowInXY )
            deallocate (WindowOutJI)

            STAT_ = SUCCESS_

        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_


    end subroutine GetSonWindow

    !--------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Gets connections matrixes for the different twoway methods  
    !>@param[in] HorizontalGridID, Connections_U, Connections_V, Connections_Z, IWD_Distances_U, IWD_Distances_V, &
    !> IWD_Distances_Z, IWD_Nodes_Z, IWD_Nodes_U, IWD_Nodes_V, STAT 
    subroutine GetConnections (HorizontalGridID, Connections_U, Connections_V, Connections_Z, &
                             IWD_Distances_U, IWD_Distances_V, IWD_Distances_Z, IWD_Nodes_Z, IWD_Nodes_U, &
                             IWD_Nodes_V, STAT)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                       :: HorizontalGridID
        integer, optional,                          intent(OUT)   :: IWD_Nodes_Z, IWD_Nodes_U, IWD_Nodes_V
        integer, optional,                          intent(INOUT) :: STAT
        integer, dimension(:,:), pointer, optional, intent(OUT)   :: Connections_U, Connections_V, Connections_Z
        real,    dimension(:  ), pointer, optional, intent(OUT)   :: IWD_Distances_U, IWD_Distances_V, IWD_Distances_Z
        !Local-----------------------------------------------------------------
        integer                                                 :: STAT_, ready_
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            if (present(Connections_U))then
                Connections_U => Me%Connections_U
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            if (present(Connections_V))then
                Connections_V => Me%Connections_V
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            if (present(Connections_Z))then
                Connections_Z => Me%Connections_Z
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            if (present(IWD_Distances_U))then
                IWD_Distances_U => Me%IWD_Distances_U
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            if (present(IWD_Distances_V))then
                IWD_Distances_V => Me%IWD_Distances_V
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            if (present(IWD_Distances_Z))then
                IWD_Distances_Z => Me%IWD_Distances_Z
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            if (present(IWD_Nodes_U))then
                IWD_Nodes_U = Me%IWD_Nodes_U
            endif
            if (present(IWD_Nodes_V))then
                IWD_Nodes_V = Me%IWD_Nodes_V
            endif
            if (present(IWD_Nodes_Z))then
                IWD_Nodes_Z = Me%IWD_Nodes_Z
            endif

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT)) &
            STAT = STAT_

                             end subroutine GetConnections
    !--------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Gets connections matrixes for the different twoway methods  
    !>@param[in] HorizontalGridID, ConnectionsU, ConnectionsV, ConnectionsZ, IWD_Distances_U, IWD_Distances_V, &
    !> IWD_Distances_Z, STAT 
    subroutine UnGetConnections(HorizontalGridID, ConnectionsU, ConnectionsV, ConnectionsZ, &
                             IWD_Distances_U, IWD_Distances_V, IWD_Distances_Z, STAT)
        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer, optional,        intent (OUT)      :: STAT
        integer,  dimension(:,:), pointer, optional :: ConnectionsU, ConnectionsV, ConnectionsZ
        real,     dimension(:  ), pointer           :: IWD_Distances_U, IWD_Distances_V, IWD_Distances_Z
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            nullify(ConnectionsU)
            call Read_UnLock(mHORIZONTALGRID_, Me%InstanceID, "UnGetConnections - ERR010")
            nullify(ConnectionsV)
            call Read_UnLock(mHORIZONTALGRID_, Me%InstanceID, "UnGetConnections - ERR020")
            nullify(ConnectionsZ)
            call Read_UnLock(mHORIZONTALGRID_, Me%InstanceID, "UnGetConnections - ERR030")
            nullify(IWD_Distances_U)
            call Read_UnLock(mHORIZONTALGRID_, Me%InstanceID, "UnGetConnections - ERR040")
            nullify(IWD_Distances_V)
            call Read_UnLock(mHORIZONTALGRID_, Me%InstanceID, "UnGetConnections - ERR050")
            nullify(IWD_Distances_Z)
            call Read_UnLock(mHORIZONTALGRID_, Me%InstanceID, "UnGetConnections - ERR060")

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT)) &
            STAT = STAT_

    end subroutine UnGetConnections

    !--------------------------------------------------------------------------

    subroutine ConstructNewFatherGrid1D(ObjHorizontalGridFather, NewFatherGrid, GridID,    &
                                      OkZ, OkU, OkV, OkCross, CheckRotation)

        !Arguments-------------------------------------------------------------
        type (T_HorizontalGrid), pointer          :: ObjHorizontalGridFather
        type (T_FatherGrid)    , pointer          :: NewFatherGrid
        integer,                intent (IN)       :: GridID
        logical,                intent (IN)       :: OkZ, OkU, OkV, OkCross
        logical,     optional,  intent (IN)       :: CheckRotation

        !Local-----------------------------------------------------------------
        real,    dimension(:,:), pointer          :: XX_Z, YY_Z, XX_U, YY_U, XX_V, YY_V, XX_Cross, YY_Cross
        integer, dimension(:,:), pointer          :: IZ,  JZ,  IU,  JU,  IV,  JV, ICross, JCross
        real                                      :: Xorig   , Yorig   , Rotation
        integer                                   :: ILB, IUB, JLB, JUB, i, j
        integer                                   :: ILBwork, IUBwork, JLBwork, JUBwork
        logical                                   :: CheckRotation_
        integer                                   :: ILBFAther, IUBFather, JLBFather, JUBFather, U, V

        !----------------------------------------------------------------------

        !if (Me%CornersXYInput) stop 'ConstructNewFatherGrid1D - ModuleHoriuzontalGrid - ERR01'

        if (present(CheckRotation)) then
            CheckRotation_ = CheckRotation
        else
            CheckRotation_ = .true.
        endif

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        ILBwork = Me%WorkSize%ILB
        IUBwork = Me%WorkSize%IUB

        JLBwork = Me%WorkSize%JLB
        JUBwork = Me%WorkSize%JUB

        NewFatherGrid%GridID            = GridID
        NewFatherGrid%OkZ               = OkZ
        NewFatherGrid%OkU               = OkU
        NewFatherGrid%OkV               = OkV
        NewFatherGrid%OkCross           = OkCross

        NewFatherGrid%CornersXYInput    = .false.

        nullify (NewFatherGrid%XX_Z )
        nullify (NewFatherGrid%YY_Z )
        nullify (NewFatherGrid%XX_U )
        nullify (NewFatherGrid%YY_U )
        nullify (NewFatherGrid%XX_V )
        nullify (NewFatherGrid%YY_V )

        nullify (NewFatherGrid%Next )
        nullify (NewFatherGrid%Prev )


        !Gets the bounds from the Father

        ILBFAther = NewFatherGrid%MPI_Window%ILB
        IUBFather = NewFatherGrid%MPI_Window%IUB
        JLBFather = NewFatherGrid%MPI_Window%JLB
        JUBFather = NewFatherGrid%MPI_Window%JUB

        NewFatherGrid%MPI_Window%ILB    = -null_int
        NewFatherGrid%MPI_Window%IUB    =  null_int
        NewFatherGrid%MPI_Window%JLB    = -null_int
        NewFatherGrid%MPI_Window%JUB    =  null_int


        !Compute points location ib the father grid
        if (NewFatherGrid%OkZ) then

            !Allocates Variables
            allocate (NewFatherGrid%XX_Z (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%YY_Z (ILB:IUB, JLB:JUB))

            !Initialize Values
            NewFatherGrid%XX_Z (:,:)  = FillValueReal
            NewFatherGrid%YY_Z (:,:)  = FillValueReal


            allocate (NewFatherGrid%IZ   (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%JZ   (ILB:IUB, JLB:JUB))
            
            allocate (NewFatherGrid%ILinkZ (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%JLinkZ (ILB:IUB, JLB:JUB))

            NewFatherGrid%IZ   (:,:)  = FillValueInt
            NewFatherGrid%JZ   (:,:)  = FillValueInt
            
            NewFatherGrid%ILinkZ (:,:) = FillValueInt
            NewFatherGrid%JLinkZ (:,:) = FillValueInt

            XX_Z  => NewFatherGrid%XX_Z
            YY_Z  => NewFatherGrid%YY_Z
            IZ    => NewFatherGrid%IZ
            JZ    => NewFatherGrid%JZ

        endif

        if (NewFatherGrid%OkU) then

            allocate (NewFatherGrid%XX_U (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%YY_U (ILB:IUB, JLB:JUB))

            NewFatherGrid%XX_U (:,:)  = FillValueReal
            NewFatherGrid%YY_U (:,:)  = FillValueReal

            allocate (NewFatherGrid%IU   (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%JU   (ILB:IUB, JLB:JUB))
            
            allocate (NewFatherGrid%ILinkU (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%JLinkU (ILB:IUB, JLB:JUB))

            NewFatherGrid%IU   (:,:)  = FillValueInt
            NewFatherGrid%JU   (:,:)  = FillValueInt
            
            NewFatherGrid%ILinkU (:,:) = FillValueInt
            NewFatherGrid%JLinkU (:,:) = FillValueInt

            XX_U  => NewFatherGrid%XX_U
            YY_U  => NewFatherGrid%YY_U
            IU    => NewFatherGrid%IU
            JU    => NewFatherGrid%JU

        endif

        if (NewFatherGrid%OkV) then

            allocate (NewFatherGrid%XX_V (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%YY_V (ILB:IUB, JLB:JUB))

            NewFatherGrid%XX_V (:,:)  = FillValueReal
            NewFatherGrid%YY_V (:,:)  = FillValueReal

            allocate (NewFatherGrid%IV   (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%JV   (ILB:IUB, JLB:JUB))
            
            allocate (NewFatherGrid%ILinkV (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%JLinkV (ILB:IUB, JLB:JUB))

            NewFatherGrid%IV   (:,:)  = FillValueInt
            NewFatherGrid%JV   (:,:)  = FillValueInt
            
            NewFatherGrid%ILinkV (:,:) = FillValueInt
            NewFatherGrid%JLinkV (:,:) = FillValueInt

            YY_V  => NewFatherGrid%YY_V
            XX_V  => NewFatherGrid%XX_V
            IV    => NewFatherGrid%IV
            JV    => NewFatherGrid%JV

        endif


        if (NewFatherGrid%OkCross) then

            allocate (NewFatherGrid%XX_Cross (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%YY_Cross (ILB:IUB, JLB:JUB))

            NewFatherGrid%XX_Cross (:,:)  = FillValueReal
            NewFatherGrid%YY_Cross (:,:)  = FillValueReal

            allocate (NewFatherGrid%ICross   (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%JCross   (ILB:IUB, JLB:JUB))

            NewFatherGrid%ICross   (:,:)  = FillValueInt
            NewFatherGrid%JCross   (:,:)  = FillValueInt

            YY_Cross  => NewFatherGrid%YY_Cross
            XX_Cross  => NewFatherGrid%XX_Cross
            ICross    => NewFatherGrid%ICross
            JCross    => NewFatherGrid%JCross

        endif

        if (Me%CornersXYInput) then

            Xorig    = - ObjHorizontalGridFather%Xorig
            Yorig    = - ObjHorizontalGridFather%Yorig

            Rotation = - ObjHorizontalGridFather%Grid_Angle

        else

            Xorig    = Me%Xorig      - ObjHorizontalGridFather%Xorig
            Yorig    = Me%Yorig      - ObjHorizontalGridFather%Yorig

            Rotation = Me%Grid_Angle - ObjHorizontalGridFather%Grid_Angle

            if (Rotation /= 0..and. CheckRotation_) then
                write(*,*) 'The submodels do not work with referentials with different rotations'
                write(*,*) 'ConstructNewFatherGrid1D - HorizontalGrid - WARN10'
            endif

            call RODAXY(0., 0., -Me%Grid_Angle, Xorig, Yorig)


            if (Xorig < 0 .and. Yorig < 0) then
                write(*,*) 'The submodels must be locate inside the father domain'
                write(*,*) 'ConstructNewFatherGrid1D - HorizontalGrid - WARN20'
            endif
        endif

        !Compute points location ib the father grid
        if (NewFatherGrid%OkZ) then

            !Rotation of the points
do3 :       do j = JLBwork, JUBwork
do4 :       do i = ILBwork, IUBwork

                if (Me%CornersXYInput) then

                    if (Me%DefineCellsMap(i,j) == 0) then

                        IZ(i, j) = FillValueInt
                        JZ(i, j) = FillValueInt
                        cycle

                    endif

                    XX_Z(i, j) = Me%Compute%XX2D_Z(i, j)
                    YY_Z(i, j) = Me%Compute%YY2D_Z(i, j)

                else

                    XX_Z(i, j) = Me%Compute%XX_Z(j)
                    YY_Z(i, j) = Me%Compute%YY_Z(i)

                endif

                call RODAXY(Xorig, Yorig, Rotation, XX_Z(i, j), YY_Z(i, j))

                call LocateCell (ObjHorizontalGridFather%Compute%XX_Z,               &
                                 ObjHorizontalGridFather%Compute%YY_Z,               &
                                 XX_Z(i, j), YY_Z(i, j),                             &
                                 ILBFAther, IUBFather, JLBFather, JUBFather,         &
                                 IZ(i, j), JZ(i, j))
                
                call LocateCell_2(ObjHorizontalGridFather%Compute%XX_U,                &
                                 ObjHorizontalGridFather%Compute%YY_V,                 &
                                 XX_Z(i, j), YY_Z(i, j),                               &
                                 ILBFAther, IUBFather, JLBFather, JUBFather,           &
                                 NewFatherGrid%ILinkZ(i, j), NewFatherGrid%JLinkZ(i, j))

                !Window
                NewFatherGrid%MPI_Window%ILB = min(NewFatherGrid%MPI_Window%ILB, IZ(i, j))
                NewFatherGrid%MPI_Window%IUB = max(NewFatherGrid%MPI_Window%IUB, IZ(i, j))
                NewFatherGrid%MPI_Window%JLB = min(NewFatherGrid%MPI_Window%JLB, JZ(i, j))
                NewFatherGrid%MPI_Window%JUB = max(NewFatherGrid%MPI_Window%JUB, JZ(i, j))

            end do do4
            end do do3

            nullify(XX_Z , YY_Z )
            nullify(IZ   , JZ   )

        endif

        if (NewFatherGrid%OkU) then

            !Rotation of the points
do1 :       do j = JLBwork, JUBwork+1
do2 :       do i = ILBwork, IUBwork

                if (Me%CornersXYInput) then

                    if (Me%DefineFacesUMap(i,j) == 0) then

                        IU(i, j) = FillValueInt
                        JU(i, j) = FillValueInt
                        cycle

                    endif

                    XX_U(i, j) = Me%Compute%XX2D_U(i, j)
                    YY_U(i, j) = Me%Compute%YY2D_U(i, j)

                else

                    XX_U(i, j) = Me%Compute%XX_U(j)
                    YY_U(i, j) = Me%Compute%YY_U(i)

                endif

                call RODAXY(Xorig, Yorig, Rotation, XX_U(i, j), YY_U(i, j))


                call LocateCell (ObjHorizontalGridFather%Compute%XX_U,               &
                                 ObjHorizontalGridFather%Compute%YY_U,               &
                                 XX_U(i, j), YY_U(i, j),                             &
                                 ILBFAther, IUBFather, JLBFather, JUBFather + 1,     &
                                 IU(i, j), JU(i, j))
                
                call LocateCell_2(ObjHorizontalGridFather%Compute%XX_Z,                &
                                 ObjHorizontalGridFather%Compute%YY_V,                 &
                                 XX_U(i, j), YY_U(i, j),                               &
                                 ILBFAther, IUBFather, JLBFather, JUBFather,           &
                                 NewFatherGrid%ILinkU(i, j), NewFatherGrid%JLinkU(i, j), U = U)

                !Window
                NewFatherGrid%MPI_Window%ILB = min(NewFatherGrid%MPI_Window%ILB, IU(i, j))
                NewFatherGrid%MPI_Window%IUB = max(NewFatherGrid%MPI_Window%IUB, IU(i, j))
                NewFatherGrid%MPI_Window%JLB = min(NewFatherGrid%MPI_Window%JLB, JU(i, j))
                NewFatherGrid%MPI_Window%JUB = max(NewFatherGrid%MPI_Window%JUB, JU(i, j))

            end do do2
            end do do1

            nullify(XX_U , YY_U )
            nullify(IU   , JU   )

        endif

        if (NewFatherGrid%OkV) then

            !Rotation of the points
do5:       do j = JLBwork, JUBwork
do6:       do i = ILBwork, IUBwork+1

                if (Me%CornersXYInput) then

                    if (Me%DefineFacesVMap(i,j) == 0) then

                        IV(i, j) = FillValueInt
                        JV(i, j) = FillValueInt
                        cycle

                    endif

                    XX_V(i, j) = Me%Compute%XX2D_V(i, j)
                    YY_V(i, j) = Me%Compute%YY2D_V(i, j)

                else

                    XX_V(i, j) = Me%Compute%XX_V(j)
                    YY_V(i, j) = Me%Compute%YY_V(i)

                endif

                call RODAXY(Xorig, Yorig, Rotation, XX_V(i, j), YY_V(i, j))

                call LocateCell (ObjHorizontalGridFather%Compute%XX_V,               &
                                 ObjHorizontalGridFather%Compute%YY_V,               &
                                 XX_V(i, j), YY_V(i, j),                             &
                                 ILBFAther, IUBFather + 1, JLBFather, JUBFather,     &
                                 IV(i, j), JV(i, j))
                
                call LocateCell_2(ObjHorizontalGridFather%Compute%XX_U,                                &
                                 ObjHorizontalGridFather%Compute%YY_Z,                                 &
                                 XX_V(i, j), YY_V(i, j),                                               &
                                 ILBFAther, IUBFather, JLBFather, JUBFather,                           &
                                 NewFatherGrid%ILinkV(i, j), NewFatherGrid%JLinkV(i, j), V = V)

                !MPI_Window
                NewFatherGrid%MPI_Window%ILB = min(NewFatherGrid%MPI_Window%ILB, IV(i, j))
                NewFatherGrid%MPI_Window%IUB = max(NewFatherGrid%MPI_Window%IUB, IV(i, j))
                NewFatherGrid%MPI_Window%JLB = min(NewFatherGrid%MPI_Window%JLB, JV(i, j))
                NewFatherGrid%MPI_Window%JUB = max(NewFatherGrid%MPI_Window%JUB, JV(i, j))

            end do do6
            end do do5

            nullify(XX_V , YY_V )
            nullify(IV   , JV   )

        endif


        if (NewFatherGrid%OkCross) then

            !Rotation of the points
do7:       do j = JLBwork, JUBwork
do8:       do i = ILBwork, IUBwork

                if (Me%CornersXYInput) then

                    if (Me%DefineCrossMap(i,j) == 0) then

                        ICross(i, j) = FillValueInt
                        JCross(i, j) = FillValueInt
                        cycle

                    endif

                    XX_Cross(i, j) = Me%XX_IE(i, j)
                    YY_Cross(i, j) = Me%YY_IE(i, j)

                else

                    XX_Cross(i, j) = Me%Compute%XX_Z(j)
                    YY_Cross(i, j) = Me%Compute%YY_Z(i)

                endif

                call RODAXY(Xorig, Yorig, Rotation, XX_Cross(i, j), YY_Cross(i, j))


                call LocateCell (ObjHorizontalGridFather%Compute%XX_Cross,           &
                                 ObjHorizontalGridFather%Compute%YY_Cross,           &
                                 XX_Cross(i, j), YY_Cross(i, j),                     &
                                 ILBFAther + 1, IUBFather + 1, JLBFather, JUBFather, &
                                 ICross(i, j), JCross(i, j))

                !MPI_Window
                NewFatherGrid%MPI_Window%ILB = min(NewFatherGrid%MPI_Window%ILB, ICross(i, j))
                NewFatherGrid%MPI_Window%IUB = max(NewFatherGrid%MPI_Window%IUB, ICross(i, j))
                NewFatherGrid%MPI_Window%JLB = min(NewFatherGrid%MPI_Window%JLB, JCross(i, j))
                NewFatherGrid%MPI_Window%JUB = max(NewFatherGrid%MPI_Window%JUB, JCross(i, j))

            end do do8
            end do do7

            nullify(XX_Cross, YY_Cross)
            nullify(ICross  , JCross  )

        endif

        !Increase MPI Window Size by one (Fluxes at the upper boundary)
        NewFatherGrid%MPI_Window%IUB = NewFatherGrid%MPI_Window%IUB + 1
        NewFatherGrid%MPI_Window%JUB = NewFatherGrid%MPI_Window%JUB + 1

    end subroutine ConstructNewFatherGrid1D


    !--------------------------------------------------------------------------

    subroutine ConstructNewFatherGrid2D(ObjHorizontalGridFather, NewFatherGrid, GridID,    &
                                      OkZ, OkU, OkV, OkCross)

        !Arguments-------------------------------------------------------------
        type (T_HorizontalGrid), pointer          :: ObjHorizontalGridFather
        type (T_FatherGrid)    , pointer          :: NewFatherGrid
        integer,                intent (IN)       :: GridID
        logical,                intent (IN)       :: OkZ, OkU, OkV, OkCross

        !Local-----------------------------------------------------------------
        real                                      :: XPoint, YPoint
        integer, dimension(:,:), pointer          :: IZ,  JZ,  IU,  JU,  IV,  JV,       &
                                                     ICross, JCross, DefinedPoint
        integer                                   :: ILB, IUB, JLB, JUB, i, j
        integer                                   :: ILBwork, IUBwork, JLBwork, JUBwork
        integer                                   :: ILBFAther, IUBFather, JLBFather, JUBFather

        !----------------------------------------------------------------------

        if (.not. Me%CornersXYInput) stop 'ConstructNewFatherGrid2D - ModuleHorizontalGrid - ERR01'

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        ILBwork = Me%WorkSize%ILB
        IUBwork = Me%WorkSize%IUB

        JLBwork = Me%WorkSize%JLB
        JUBwork = Me%WorkSize%JUB

        NewFatherGrid%GridID            = GridID
        NewFatherGrid%OkZ               = OkZ
        NewFatherGrid%OkU               = OkU
        NewFatherGrid%OkV               = OkV
        NewFatherGrid%OkCross           = OkCross

        !Gets the bounds from the Father

        ILBFAther = NewFatherGrid%MPI_Window%ILB
        IUBFather = NewFatherGrid%MPI_Window%IUB
        JLBFather = NewFatherGrid%MPI_Window%JLB
        JUBFather = NewFatherGrid%MPI_Window%JUB

        NewFatherGrid%MPI_Window%ILB    = -null_int
        NewFatherGrid%MPI_Window%IUB    =  null_int
        NewFatherGrid%MPI_Window%JLB    = -null_int
        NewFatherGrid%MPI_Window%JUB    =  null_int

        NewFatherGrid%CornersXYInput    = .true.

        nullify (NewFatherGrid%Next )
        nullify (NewFatherGrid%Prev )

        write(*,*) 'The submodels in curvilinear grids must have identical distortion'
        write(*,*) 'of the coarser grid'
        write(*,*) 'ConstructNewFatherGrid2D - HorizontalGrid - WARN01'


        !Compute points location ib the father grid
        if (NewFatherGrid%OkZ) then

            !Allocates Variables
            allocate (NewFatherGrid%IZ   (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%JZ   (ILB:IUB, JLB:JUB))

            !Initialize Values

            NewFatherGrid%IZ   (:,:)  = FillValueInt
            NewFatherGrid%JZ   (:,:)  = FillValueInt

            IZ    => NewFatherGrid%IZ
            JZ    => NewFatherGrid%JZ

            NewFatherGrid%XX_Z => Me%Compute%XX2D_Z
            NewFatherGrid%YY_Z => Me%Compute%YY2D_Z

            DefinedPoint    => ObjHorizontalGridFather%DefineCellsMap

        endif

        if (NewFatherGrid%OkU) then

            allocate (NewFatherGrid%IU   (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%JU   (ILB:IUB, JLB:JUB))

            NewFatherGrid%IU   (:,:)  = FillValueInt
            NewFatherGrid%JU   (:,:)  = FillValueInt

            IU    => NewFatherGrid%IU
            JU    => NewFatherGrid%JU

            NewFatherGrid%XX_U => Me%Compute%XX2D_U
            NewFatherGrid%YY_U => Me%Compute%YY2D_U


            DefinedPoint    => ObjHorizontalGridFather%DefineFacesUMap

        endif

        if (NewFatherGrid%OkV) then

            allocate (NewFatherGrid%IV   (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%JV   (ILB:IUB, JLB:JUB))

            NewFatherGrid%IV   (:,:)  = FillValueInt
            NewFatherGrid%JV   (:,:)  = FillValueInt

            IV    => NewFatherGrid%IV
            JV    => NewFatherGrid%JV

            NewFatherGrid%XX_V => Me%Compute%XX2D_V
            NewFatherGrid%YY_V => Me%Compute%YY2D_V

            DefinedPoint    => ObjHorizontalGridFather%DefineFacesVMap

        endif


        if (NewFatherGrid%OkCross) then

            allocate (NewFatherGrid%ICross   (ILB:IUB, JLB:JUB))
            allocate (NewFatherGrid%JCross   (ILB:IUB, JLB:JUB))

            NewFatherGrid%ICross   (:,:)  = FillValueInt
            NewFatherGrid%JCross   (:,:)  = FillValueInt

            ICross    => NewFatherGrid%ICross
            JCross    => NewFatherGrid%JCross

            NewFatherGrid%XX_Cross => Me%XX_IE
            NewFatherGrid%YY_Cross => Me%YY_IE

            DefinedPoint    => ObjHorizontalGridFather%DefineCrossMap

        endif


        !Compute points location ib the father grid
        if (NewFatherGrid%OkZ) then

            !Rotation of the points
do3 :       do j = JLBwork, JUBwork
do4 :       do i = ILBwork, IUBwork

                if (Me%CornersXYInput) then

                    XPoint = Me%Compute%XX2D_Z(i, j)
                    YPoint = Me%Compute%YY2D_Z(i, j)

                else

                    XPoint = Me%Compute%XX_Z(j)
                    YPoint = Me%Compute%YY_Z(i)

                endif


                call LocateCellPolygons(ObjHorizontalGridFather%Compute%XX2D_Z,         &
                                        ObjHorizontalGridFather%Compute%YY2D_Z,         &
                                        XPoint, YPoint, DefinedPoint,                   &
                                        ILBFAther, IUBFather, JLBFather, JUBFather,     &
                                        IZ(i, j), JZ(i, j))


                !Window
                NewFatherGrid%MPI_Window%ILB = min(NewFatherGrid%MPI_Window%ILB, IZ(i, j))
                NewFatherGrid%MPI_Window%IUB = max(NewFatherGrid%MPI_Window%IUB, IZ(i, j))
                NewFatherGrid%MPI_Window%JLB = min(NewFatherGrid%MPI_Window%JLB, JZ(i, j))
                NewFatherGrid%MPI_Window%JUB = max(NewFatherGrid%MPI_Window%JUB, JZ(i, j))

            end do do4
            end do do3

            nullify(IZ   , JZ   )

        endif

        if (NewFatherGrid%OkU) then

            !Rotation of the points
do1 :       do j = JLBwork, JUBwork+1
do2 :       do i = ILBwork, IUBwork

                call LocateCellPolygons(ObjHorizontalGridFather%Compute%XX2D_U,         &
                                        ObjHorizontalGridFather%Compute%YY2D_U,         &
                                        Me%Compute%XX2D_U(i, j),                        &
                                        Me%Compute%YY2D_U(i, j),                        &
                                        DefinedPoint,                                   &
                                        ILBFAther, IUBFather, JLBFather, JUBFather + 1, &
                                        IU(i, j), JU(i, j))

                !Window
                NewFatherGrid%MPI_Window%ILB = min(NewFatherGrid%MPI_Window%ILB, IU(i, j))
                NewFatherGrid%MPI_Window%IUB = max(NewFatherGrid%MPI_Window%IUB, IU(i, j))
                NewFatherGrid%MPI_Window%JLB = min(NewFatherGrid%MPI_Window%JLB, JU(i, j))
                NewFatherGrid%MPI_Window%JUB = max(NewFatherGrid%MPI_Window%JUB, JU(i, j))

            end do do2
            end do do1

            nullify(IU   , JU   )

        endif

        if (NewFatherGrid%OkV) then

            !Rotation of the points
do5:       do j = JLBwork, JUBwork
do6:       do i = ILBwork, IUBwork+1

                call LocateCellPolygons(ObjHorizontalGridFather%Compute%XX2D_V,         &
                                        ObjHorizontalGridFather%Compute%YY2D_V,         &
                                        Me%Compute%XX2D_V(i, j),                        &
                                        Me%Compute%YY2D_V(i, j),                        &
                                        DefinedPoint,                                   &
                                        ILBFAther, IUBFather + 1, JLBFather, JUBFather, &
                                        IV(i, j), JV(i, j))

                !MPI_Window
                NewFatherGrid%MPI_Window%ILB = min(NewFatherGrid%MPI_Window%ILB, IV(i, j))
                NewFatherGrid%MPI_Window%IUB = max(NewFatherGrid%MPI_Window%IUB, IV(i, j))
                NewFatherGrid%MPI_Window%JLB = min(NewFatherGrid%MPI_Window%JLB, JV(i, j))
                NewFatherGrid%MPI_Window%JUB = max(NewFatherGrid%MPI_Window%JUB, JV(i, j))

            end do do6
            end do do5

            nullify(IV   , JV   )

        endif


        if (NewFatherGrid%OkCross) then

            !Rotation of the points
do7:       do j = JLBwork, JUBwork
do8:       do i = ILBwork, IUBwork

                call LocateCellPolygons(ObjHorizontalGridFather%XX_IE,                  &
                                        ObjHorizontalGridFather%YY_IE,                  &
                                        Me%XX_IE(i, j),                                 &
                                        Me%YY_IE(i, j),                                 &
                                        DefinedPoint,                                   &
                                        ILBFAther + 1, IUBFather + 1, JLBFather, JUBFather,&
                                        ICross(i, j), JCross(i, j))

                !MPI_Window
                NewFatherGrid%MPI_Window%ILB = min(NewFatherGrid%MPI_Window%ILB, ICross(i, j))
                NewFatherGrid%MPI_Window%IUB = max(NewFatherGrid%MPI_Window%IUB, ICross(i, j))
                NewFatherGrid%MPI_Window%JLB = min(NewFatherGrid%MPI_Window%JLB, JCross(i, j))
                NewFatherGrid%MPI_Window%JUB = max(NewFatherGrid%MPI_Window%JUB, JCross(i, j))

            end do do8
            end do do7

            nullify(ICross  , JCross  )

        endif

        !Increase MPI Window Size by one (Fluxes at the upper boundary)
        NewFatherGrid%MPI_Window%IUB = NewFatherGrid%MPI_Window%IUB + 1
        NewFatherGrid%MPI_Window%JUB = NewFatherGrid%MPI_Window%JUB + 1

                                      end subroutine ConstructNewFatherGrid2D

    !-------------------------------------------------------------------------------------

    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Builds connection matrix between parent and child grid on a IWD interpolation
    !>@param[in] FatherID, SonID
    subroutine ConstructP2C_IWD(FatherID, SonID)

        !Arguments--------------------------------------------------------------
        type (T_HorizontalGrid), pointer    :: ObjHorizontalFather
        integer                             :: FatherID, SonID
        !Local------------------------------------------------------------------
        real                                :: FatherCenterX, FatherCenterY, DistanceToFather, SonCenterX, SonCenterY
        real                                :: SearchRadious, MaxRatio
        integer                             :: index, Nbr_Connections, minJ, minI, maxJ, maxI
        integer                             :: i, j, i2, j2, ready_
        integer, dimension (:, :), pointer  :: FatherLinkI, FatherLinkJ
        !-------------------------------------------------------------------------
        index = 1
        call Ready (SonID, ready_)

        call LocateObjFather        (ObjHorizontalFather, FatherID)

        call DetermineMaxRatio(ObjHorizontalFather, MaxRatio)

        FatherLinkI      => Me%LastFatherGrid%ILinkZ
        FatherLinkJ      => Me%LastFatherGrid%JLinkZ
            
        minJ = min(FatherLinkJ(1,1), FatherLinkJ(1, Me%Size%JUB - 1), &
                    FatherLinkJ(Me%Size%IUB - 1, 1), FatherLinkJ(Me%Size%IUB - 1, Me%Size%JUB - 1))
        minI = min(FatherLinkI(1,1), FatherLinkI(1, Me%Size%JUB - 1), &
                    FatherLinkI(Me%Size%IUB - 1, 1), FatherLinkI(Me%Size%IUB - 1, Me%Size%JUB - 1))
        maxJ = max(FatherLinkJ(1,1), FatherLinkJ(1, Me%Size%JUB - 1), &
                    FatherLinkJ(Me%Size%IUB - 1, 1), FatherLinkJ(Me%Size%IUB - 1, Me%Size%JUB - 1))
        maxI = max(FatherLinkI(1,1), FatherLinkI(1, Me%Size%JUB - 1), &
                    FatherLinkI(Me%Size%IUB - 1, 1), FatherLinkI(Me%Size%IUB - 1, Me%Size%JUB - 1))
        
        nullify(FatherLinkI)
        nullify(FatherLinkJ)
        
        !Uses the maxRatio to avoid allocating too few indexes. 2nd term is to account for search radious
        Nbr_Connections   = (maxJ - minJ + 2) * (maxI - minI + 2) * (MaxRatio + 4 * (int(sqrt(MaxRatio))) + 4)
        !Vectorials and scalars
        allocate (Me%Connections_Z (Nbr_Connections, 4))
        allocate (Me%IWD_Distances_Z   (Nbr_Connections))
        Me%Connections_Z(:, :) = FillValueInt
        Me%IWD_Distances_Z(:)      = FillValueReal
        allocate (Me%Connections_U(Nbr_Connections, 4))
        allocate (Me%IWD_Distances_U  (Nbr_Connections))
        Me%Connections_U(:, :) = FillValueInt
        Me%IWD_Distances_U(:)      = FillValueReal
        allocate (Me%Connections_V (Nbr_Connections, 4))
        allocate (Me%IWD_Distances_V   (Nbr_Connections))
        Me%Connections_V(:, :) = FillValueInt
        Me%IWD_Distances_V(:)      = FillValueReal
        
        !i and j   -> father cell
        !i2 and j2 -> son cell
        do j = minJ, maxJ
        do i = minI, maxI

            !Find Father's cell center coordinates
            FatherCenterX = (( ObjHorizontalFather%XX_IE(i, j  ) +  ObjHorizontalFather%XX_IE(i+1, j  ))/2. + &
                                ( ObjHorizontalFather%XX_IE(i, j+1) +  ObjHorizontalFather%XX_IE(i+1, j+1))/2.)/2.
            FatherCenterY = (( ObjHorizontalFather%YY_IE(i, j  ) +  ObjHorizontalFather%YY_IE(i+1, j  ))/2. + &
                                ( ObjHorizontalFather%YY_IE(i, j+1) +  ObjHorizontalFather%YY_IE(i+1, j+1))/2.)/2.
            
            SearchRadious = (1.01+(1/(Sqrt(MaxRatio)))) * Sqrt((FatherCenterX - ObjHorizontalFather%XX(j))**2 + &
                                                                (FatherCenterY - ObjHorizontalFather%YY(i))**2)

                !Find and build matrix of correspondent son cells
            do j2 = 1, Me%Size%JUB - 1
            do i2 = 1, Me%Size%IUB - 1
                SonCenterX = (( Me%XX_IE(i2, j2  ) +  Me%XX_IE(i2+1, j2  ))/2. + &
                                ( Me%XX_IE(i2, j2+1) +  Me%XX_IE(i2+1, j2+1))/2.)/2.
                SonCenterY = (( Me%YY_IE(i2, j2  ) +  Me%YY_IE(i2+1, j2  ))/2. + &
                                ( Me%YY_IE(i2, j2+1) +  Me%YY_IE(i2+1, j2+1))/2.)/2.

                DistanceToFather = Sqrt((SonCenterX - FatherCenterX)**2.0 + &
                                        (SonCenterY - FatherCenterY)**2.0)
                if (DistanceToFather <= SearchRadious) then
                    Me%Connections_Z(index, 1) = i
                    Me%Connections_Z(index, 2) = j
                    Me%Connections_Z(index, 3) = i2
                    Me%Connections_Z(index, 4) = j2

                    if (DistanceToFather == 0)then
                        Me%IWD_Distances_Z(index) = 1.e-5
                    else
                        Me%IWD_Distances_Z(index) = DistanceToFather
                    endif
                    index = index + 1
                endif
            enddo
            enddo

        enddo
        enddo
        Me%IWD_Nodes_Z = index - 1
        
        call  ConstructIWDVel(ObjHorizontalFather, minI, minJ, maxI, maxJ, MaxRatio)

    end subroutine ConstructP2C_IWD

    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Builds connection matrix for velocities cells between son and father grid on a IWD interpolation
    !>@param[in] ObjHorizontalFather, IWD_Distances, IWD_connections, minI, minJ, maxI, maxJ, MaxRatio
    subroutine ConstructIWDVel(ObjHorizontalFather, minI, minJ, maxI, maxJ, MaxRatio)
        !Arguments--------------------------------------------------------------
        type (T_HorizontalGrid), pointer    :: ObjHorizontalFather
        integer                             :: minJ, minI, maxJ, maxI
        real                                :: MaxRatio
        !Local------------------------------------------------------------------
        real                                :: FatherCenterX_U, FatherCenterX_V, FatherCenterY_U, FatherCenterY_V, &
                                               FatherCenterX_Z, FatherCenterY_Z, SonCenterX_U, SonCenterX_V,       &
                                               SonCenterY_U, SonCenterY_V, DistanceToFather_U, DistanceToFather_V, &
                                               SearchRadious_U, SearchRadious_V
        integer                             :: index_U, index_V, i, j, i2, j2
        !-------------------------------------------------------------------------
        index_U = 1
        index_V = 1
        do j = minJ, maxJ
        do i = minI, maxI

            !Find Father cell center for U cell and V cell

            FatherCenterX_U = (ObjHorizontalFather%XX_IE(i, j  ))/2. +  (ObjHorizontalFather%XX_IE(i+1, j  ))/2.
            FatherCenterY_U = (ObjHorizontalFather%YY_IE(i, j  ))/2. +  (ObjHorizontalFather%YY_IE(i+1, j  ))/2.

            FatherCenterX_V = (ObjHorizontalFather%XX_IE(i, j  ))/2. +  (ObjHorizontalFather%XX_IE(i  , j+1))/2.
            FatherCenterY_V = (ObjHorizontalFather%YY_IE(i, j  ))/2. +  (ObjHorizontalFather%YY_IE(i  , j+1))/2.

            FatherCenterX_Z = (( ObjHorizontalFather%XX_IE(i, j  ) +  ObjHorizontalFather%XX_IE(i+1, j  ))/2. + &
                               ( ObjHorizontalFather%XX_IE(i, j+1) +  ObjHorizontalFather%XX_IE(i+1, j+1))/2.)/2.
            FatherCenterY_Z = (( ObjHorizontalFather%YY_IE(i, j  ) +  ObjHorizontalFather%YY_IE(i+1, j  ))/2. + &
                               ( ObjHorizontalFather%YY_IE(i, j+1) +  ObjHorizontalFather%YY_IE(i+1, j+1))/2.)/2.

            SearchRadious_U = (1.01+(1/(Sqrt(MaxRatio)))) * &
                              Sqrt((FatherCenterX_U - FatherCenterX_Z)**2 + &
                                   (FatherCenterY_U - ObjHorizontalFather%YY_IE(i, j))**2)
            SearchRadious_V = (1.01+(1/(Sqrt(MaxRatio)))) * &
                              Sqrt((FatherCenterX_V - ObjHorizontalFather%XX_IE(i, j))**2 + &
                                   (FatherCenterY_V - FatherCenterY_Z)**2)

            ! Find and build matrix of correspondent son cells
            do j2 = 1, Me%Size%JUB - 1
            do i2 = 1, Me%Size%IUB - 1
                SonCenterX_U = (Me%XX_IE(i2, j2  ))/2. +  (Me%XX_IE(i2+1, j2  ))/2.
                SonCenterY_U = (Me%YY_IE(i2, j2  ))/2. +  (Me%YY_IE(i2+1, j2  ))/2.

                SonCenterX_V = (Me%XX_IE(i2, j2  ))/2. +  (Me%XX_IE(i2  , j2+1))/2.
                SonCenterY_V = (Me%YY_IE(i2, j2  ))/2. +  (Me%YY_IE(i2  , j2+1))/2.

                DistanceToFather_U = Sqrt((SonCenterX_U - FatherCenterX_U)**2.0 + &
                                          (SonCenterY_U - FatherCenterY_U)**2.0)
                DistanceToFather_V = Sqrt((SonCenterX_V - FatherCenterX_V)**2.0 + &
                                          (SonCenterY_V - FatherCenterY_V)**2.0)

                if (DistanceToFather_U <= SearchRadious_U) then
                    Me%Connections_U(index_U, 1) = i
                    Me%Connections_U(index_U, 2) = j
                    Me%Connections_U(index_U, 3) = i2
                    Me%Connections_U(index_U, 4) = j2

                    if (DistanceToFather_U == 0)then
                        !The 0.001 is the reference distance. The if also avoids /0 in module functions
                        Me%IWD_Distances_U(index_U) = 1.e-5
                    else
                        Me%IWD_Distances_U(index_U) = DistanceToFather_U
                    endif

                    index_U = index_U + 1
                endif

                if (DistanceToFather_V <= SearchRadious_V) then
                    Me%Connections_V(index_V, 1) = i
                    Me%Connections_V(index_V, 2) = j
                    Me%Connections_V(index_V, 3) = i2
                    Me%Connections_V(index_V, 4) = j2

                    if (DistanceToFather_V == 0)then
                        !The 0.001 is the reference distance. The if also avoids /0 in module functions
                        Me%IWD_Distances_V(index_V) = 1.e-5
                    else
                        Me%IWD_Distances_V(index_V) = DistanceToFather_V
                    endif

                    index_V = index_V + 1
                endif
            enddo
            enddo

        enddo
        enddo
        Me%IWD_Nodes_U = index_U - 1
        Me%IWD_Nodes_V = index_V - 1
    end subroutine ConstructIWDVel
    !---------------------------------------------------------------------------
    
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>builds Connection matrix for Z cells, which will be used to detect TwoWay momentum discharges
    !>@param[in] FatherID, SonID     
    subroutine ConstructP2C_Avrg(FatherID, SonID)
        !Arguments--------------------------------------------------------------
        type (T_HorizontalGrid), pointer    :: ObjHorizontalFather
        !Local------------------------------------------------------------------
        integer                             :: FatherID, SonID
        integer                             :: index, Nbr_Connections, minJ, minI, maxJ, maxI
        integer                             :: i, j, i2, j2, ready_
        integer, dimension (:, :), pointer  :: FatherLinkI, FatherLinkJ
        !-------------------------------------------------------------------------
        index = 1
        call Ready (SonID, ready_)

        call LocateObjFather        (ObjHorizontalFather, FatherID)
        
        FatherLinkI      => Me%LastFatherGrid%ILinkZ
        FatherLinkJ      => Me%LastFatherGrid%JLinkZ
            
        minJ = min(FatherLinkJ(1,1), FatherLinkJ(1, Me%Size%JUB - 1), &
                    FatherLinkJ(Me%Size%IUB - 1, 1), FatherLinkJ(Me%Size%IUB - 1, Me%Size%JUB - 1))
        minI = min(FatherLinkI(1,1), FatherLinkI(1, Me%Size%JUB - 1), &
                    FatherLinkI(Me%Size%IUB - 1, 1), FatherLinkI(Me%Size%IUB - 1, Me%Size%JUB - 1))
        maxJ = max(FatherLinkJ(1,1), FatherLinkJ(1, Me%Size%JUB - 1), &
                    FatherLinkJ(Me%Size%IUB - 1, 1), FatherLinkJ(Me%Size%IUB - 1, Me%Size%JUB - 1))
        maxI = max(FatherLinkI(1,1), FatherLinkI(1, Me%Size%JUB - 1), &
                    FatherLinkI(Me%Size%IUB - 1, 1), FatherLinkI(Me%Size%IUB - 1, Me%Size%JUB - 1))

        
        Nbr_Connections = Me%Size%JUB * Me%Size%IUB + 1
        
        allocate (Me%Connections_Z (Nbr_Connections, 4))
        !i and j   -> father cell
        !i2 and j2 -> son cell
        do j = minJ, maxJ
        do i = minI, maxI

                !Find and build matrix of correspondent son cells
            do j2 = 1, Me%Size%JUB - 1
            do i2 = 1, Me%Size%IUB - 1

                if (FatherLinkI(i2, j2) == i) then
                    if (FatherLinkJ(i2, j2) == j) then
                        
                        Me%Connections_Z(index, 1) = i
                        Me%Connections_Z(index, 2) = j
                        Me%Connections_Z(index, 3) = i2
                        Me%Connections_Z(index, 4) = j2

                        index = index + 1
                    endif
                endif
                
            enddo
            enddo

        enddo
        enddo
        
        nullify(FatherLinkI)
        nullify(FatherLinkJ)
        
    end subroutine ConstructP2C_Avrg
    
    !---------------------------------------------------------------------------

    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Computes Maximum grid ratio between father and son domains
    !>@param[in] ObjHorizontalGridFather, Ratio
    subroutine DetermineMaxRatio(ObjHorizontalGridFather, MaxRatio)
        !Argumments---------------------------------------------------------------------
        type (T_HorizontalGrid), pointer    :: ObjHorizontalGridFather
        real, intent(OUT)                   :: MaxRatio
        !local--------------------------------------------------------------------------
        integer                             ::  i, j
        real                                :: aux

        !-------------------------------------------------------------------------------
        MaxRatio = 1.0
        aux      = 1.0

        do j = 1, Me%Size%JUB - 1
        do i = 1, Me%Size%IUB - 1
            aux = ObjHorizontalGridFather%GridCellArea(Me%LastFatherGrid%IZ(i, j), Me%LastFatherGrid%JZ(i, j)) / &
                  Me%GridCellArea(i, j)
            if (aux > MaxRatio) then
                MaxRatio = aux
            endif
        enddo
        enddo

    end subroutine DetermineMaxRatio
    !--------------------------------------------------------------------------
    ! This subroutine adds a new grid to the father grid List

    subroutine Add_FatherGrid(NewFatherGrid)

        !Arguments-------------------------------------------------------------
        type(T_FatherGrid),         pointer     :: NewFatherGrid

        !----------------------------------------------------------------------

        ! Add to the WaterFatherGrid List a new father grid
        if (.not. associated(Me%FirstFatherGrid)) then

            Me%FirstFatherGrid    => NewFatherGrid
            Me%LastFatherGrid     => NewFatherGrid
        else
            NewFatherGrid%Prev     => Me%LastFatherGrid
            Me%LastFatherGrid%Next => NewFatherGrid
            Me%LastFatherGrid      => NewFatherGrid

        end if

        !----------------------------------------------------------------------

    end subroutine Add_FatherGrid

    !--------------------------------------------------------------------------

    subroutine Search_FatherGrid(NewFatherGrid, GridID)

        !Arguments-------------------------------------------------------------
        type(T_FatherGrid),         pointer     :: NewFatherGrid
        integer,                    intent (IN) :: GridID


        !Local-----------------------------------------------------------------

        NewFatherGrid => Me%FirstFatherGrid

do1 :   do while (associated(NewFatherGrid))

cd1 :       if (NewFatherGrid%GridID == GridID) then
                exit do1
            else
                NewFatherGrid => NewFatherGrid%Next
            end if cd1

        end do do1

        if (.not.Associated(NewFatherGrid))                                              &
            call SetError(FATAL_, INTERNAL_, 'Search_FatherGrid - HorizontalGrid - ERR01')

        !----------------------------------------------------------------------

    end subroutine Search_FatherGrid

    !--------------------------------------------------------------------------

    subroutine ConstructGlobalVariables

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, ZoneLong
        integer, dimension(2)               :: AuxInt
        real,    dimension(2)               :: AuxReal
        real                                :: DX, DY, Aux, XY_Aux
        integer                             :: flag, flag1, flag2
        integer                             :: ClientNumber
        logical                             :: BlockFound, ConstantSpacingX, ConstantSpacingY
        integer                             :: FirstLine, LastLine, line, i, j, ii, jj, iflag

        !----------------------------------------------------------------------

        !Nullify T_Distances
        nullify (Me%XX )
        nullify (Me%YY )
        nullify (Me%DXX)
        nullify (Me%DYY)
        nullify (Me%DZX)
        nullify (Me%DZY)
        nullify (Me%DUX)
        nullify (Me%DUY)
        nullify (Me%DVX)
        nullify (Me%DVY)
        nullify (Me%XX_IE)
        nullify (Me%YY_IE)
        nullify (Me%Compute%XX_Z )
        nullify (Me%Compute%YY_Z )
        nullify (Me%Compute%XX_U )
        nullify (Me%Compute%YY_U )
        nullify (Me%Compute%XX_V )
        nullify (Me%Compute%YY_V )
        nullify (Me%Compute%XX_Cross)
        nullify (Me%Compute%YY_Cross)
        nullify (Me%XX_AlongGrid)
        nullify (Me%YY_AlongGrid)

        !Nullify Other
        nullify (Me%LatitudeZ    )
        nullify (Me%LongitudeZ   )
        nullify (Me%F            )
        nullify (Me%GridCellArea )
        nullify (Me%LatitudeConn )
        nullify (Me%LongitudeConn)


        !Opens Data File
        call ConstructEnterData(Me%ObjEnterData, Me%FileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR10'


        !Reads ILB_IUB
        call GetData(AuxInt,Me%ObjEnterData, flag,                                      &
                     keyword      = 'ILB_IUB',                                          &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR20'
        if (flag      /= 2       ) stop 'ConstructGlobalVariables - HorizontalGrid - ERR30'

        Me%GlobalWorkSize%ILB = AuxInt(1)
        Me%GlobalWorkSize%IUB = AuxInt(2)

        !Reads JLB_JUB
        call GetData(AuxInt,Me%ObjEnterData, flag,                                      &
                     keyword      = 'JLB_JUB',                                          &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR40'
        if (flag      /= 2       ) stop 'ConstructGlobalVariables - HorizontalGrid - ERR50'

        Me%GlobalWorkSize%JLB = AuxInt(1)
        Me%GlobalWorkSize%JUB = AuxInt(2)

        call GetData(Value          = Me%DDecomp%Halo_Points,                           &
                     EnterDataID    = Me%ObjEnterData,                                  &
                     flag           = iflag,                                            &
                     keyword        = 'HALOPOINTS',                                     &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'HorizontalGrid',                                 &
                     default        = 3,                                                &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR440'

        if (Me%DDecomp%Halo_Points < 2) then
            write(*,*) 'HALOPOINTS must equal to 2 or bigger'
            stop 'ConstructGlobalVariables - HorizontalGrid - ERR445'
        endif

        !Intialization of domain decomposition procedure
        if (Me%DDecomp%ON) then
            call ConstructDDecomp
        endif


        if (Me%DDecomp%MasterOrSlave)  then

            Me%Size%ILB     = Me%WorkSize%ILB-1
            Me%Size%IUB     = Me%WorkSize%IUB+1

            if (Me%GlobalWorkSize%ILB /= Me%DDecomp%Global%ILB) then
                write(*,*) "Me%GlobalWorkSize%ILB",Me%GlobalWorkSize%ILB
                write(*,*) "Me%DDecomp%Global%ILB",Me%DDecomp%Global%ILB
                stop 'ConstructGlobalVariables - HorizontalGrid - ERR32'
            endif

            if (Me%GlobalWorkSize%IUB /= Me%DDecomp%Global%IUB) then
                stop 'ConstructGlobalVariables - HorizontalGrid - ERR34'
            endif

        else

            Me%Size%ILB     = Me%GlobalWorkSize%ILB-1
            Me%Size%IUB     = Me%GlobalWorkSize%IUB+1

            Me%WorkSize%ILB = Me%GlobalWorkSize%ILB
            Me%WorkSize%IUB = Me%GlobalWorkSize%IUB

        endif


        if (Me%DDecomp%MasterOrSlave)  then

            Me%Size%JLB     = Me%WorkSize%JLB-1
            Me%Size%JUB     = Me%WorkSize%JUB+1

            if (Me%GlobalWorkSize%JLB /= Me%DDecomp%Global%JLB) then
                stop 'ConstructGlobalVariables - HorizontalGrid - ERR52'
            endif

            if (Me%GlobalWorkSize%JUB /= Me%DDecomp%Global%JUB) then
                stop 'ConstructGlobalVariables - HorizontalGrid - ERR54'
            endif

        else

            Me%Size%JLB     = Me%GlobalWorkSize%JLB-1
            Me%Size%JUB     = Me%GlobalWorkSize%JUB+1

            Me%WorkSize%JLB = Me%GlobalWorkSize%JLB
            Me%WorkSize%JUB = Me%GlobalWorkSize%JUB

        endif


        !Reads COORD_TIP
        call GetData(Me%CoordType, Me%ObjEnterData, flag,                               &
                     keyword      = 'COORD_TIP',                                        &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR60'

        if (.NOT. ((Me%CoordType == GEOG_       ) .OR.                                  &
                   (Me%CoordType == UTM_        ) .OR.                                  &
                   (Me%CoordType == MIL_PORT_   ) .OR.                                  &
                   (Me%CoordType == SIMPLE_GEOG_) .OR.                                  &
                   (Me%CoordType == GRID_COORD_ ) .OR.                                  &
                   (Me%CoordType == CIRCULAR_   ) .OR.                                  &
                   (Me%CoordType == NLRD_       ))) then
            call SetError (FATAL_, INTERNAL_, "Invalid Coordinates")
        endif


        !Reads Latitude (Reference latitude. If externally specificed is the center of domain)
        call GetData(Me%Latitude, Me%ObjEnterData, flag,                                &
                     keyword      = 'LATITUDE',                                         &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR70'

        !Reads Longitude (Reference longitude. If externally specificed is the center of domain)
        call GetData(Me%Longitude, Me%ObjEnterData, flag,                                  &
                     keyword      = 'LONGITUDE',                                        &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR80'

        call GetData(Me%Datum, Me%ObjEnterData, flag,                                      &
                     keyword      = 'DATUM',                                            &
                     default      = WGS_84_DATUM,                                       &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR90'

        !Valor metrico a acrescentar as coordenadas na direccao x
        call GetData(Me%Easting, Me%ObjEnterData, flag,                                    &
                     keyword      = 'EASTING',                                          &
                     default      = 0.0,                                                &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR100'

        !Valor metrico a acrescentar as coordenadas na direccao x
        call GetData(Me%Northing, Me%ObjEnterData, flag,                                   &
                     keyword      = 'NORTHING',                                         &
                     default      = 0.0,                                                &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR110'

        !Projection Type to convert to cartesian coordinates
        call GetData(Me%ProjType, Me%ObjEnterData, flag,                                   &
                     keyword      = 'PROJ_TYPE',                                        &
                     default      = PAULO_PROJECTION_,                                  &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR120'

        !Lambert Conformal Conic Projection
        call GetData(Me%SP1, Me%ObjEnterData, flag1,                                       &
                     keyword      = 'SP1',                                              &
                     default      = 30.,                                                &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR120'

        !Lambert Conformal Conic Projection
        call GetData(Me%SP2, Me%ObjEnterData, flag2,                                       &
                     keyword      = 'SP2',                                              &
                     default      = 60.,                                                &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR130'

#ifdef _USE_PROJ4
#else
        if (flag1 + flag2 == 2) then  !estao la
            Me%UseLambert = .TRUE.
        else
            Me%UseLambert = .FALSE.
        endif
#endif

         !Reads ZoneLong
        if (Me%CoordType == UTM_ ) then
            call ComputeGridZone (dble(Me%Longitude), dble(Me%Latitude), Me%Grid_Zone)

            Me%ZoneLong = Me%Grid_Zone(1)
            Me%ZoneLat  = Me%Grid_Zone(2)

            call GetData(ZoneLong, Me%ObjEnterData, flag,                                  &
                         keyword      = 'ZONE',                                         &
                         ClientModule = 'HorizontalGrid',                               &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR140'
            if (flag == 1) then
                if (ZoneLong /= Me%ZoneLong) stop 'ConstructGlobalVariables - HorizontalGrid - ERR150'
            endif

        else if (Me%CoordType .EQ. MIL_PORT_) then
            Me%ZoneLong = PORTUGUESE_UTM_ZONE_
        else
            Me%ZoneLong = null_int
        endif


        !Reads ORIGIN
        call GetData(AuxReal,Me%ObjEnterData, flag,                                        &
                     keyword      = 'ORIGIN',                                           &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR160'
        if (flag      /= 2       ) stop 'ConstructGlobalVariables - HorizontalGrid - ERR170'

        Me%Xorig = AuxReal(1)
        Me%Yorig = AuxReal(2)

        !Reads GridAngle
        call GetData(Me%Grid_Angle, Me%ObjEnterData, flag,                                 &
                     keyword      = 'GRID_ANGLE',                                       &
                     default      = 0.0,                                                &
                     ClientModule = 'HorizontalGrid',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR180'


        !Allocates XX and YY
        allocate(Me%XX(Me%Size%JLB+1 : Me%Size%JUB+1))
        allocate(Me%YY(Me%Size%ILB+1 : Me%Size%IUB+1))

        allocate(Me%XX_IE(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))
        allocate(Me%YY_IE(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))

        allocate(Me%LatitudeConn(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))
        allocate(Me%LongitudeConn(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))

        allocate(Me%DefineCellsMap(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))

        allocate(Me%XX_AlongGrid(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))
        allocate(Me%YY_AlongGrid(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))


        Me%XX    = null_real
        Me%YY    = null_real
        Me%XX_IE = null_real
        Me%YY_IE = null_real

        Me%LatitudeConn  = null_real
        Me%LongitudeConn = null_real


        Me%DefineCellsMap(:,:)  = 1

        !Reads XX_YE, YY_IE
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR190'

        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                 &
                                    BeginCornersXY, EndCornersXY, BlockFound,   &
                                    FirstLine = FirstLine, LastLine = LastLine, &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR200'

BF:     if (BlockFound) then

            Me%CornersXYInput = .true.

            if (Me%CoordType == CIRCULAR_ .and. Me%CoordType == GEOG_) then
                stop 'ConstructGlobalVariables - HorizontalGrid - ERR210'
            endif


            line = FirstLine

            do i = Me%GlobalWorkSize%ILB, Me%GlobalWorkSize%IUB+1
            do j = Me%GlobalWorkSize%JLB, Me%GlobalWorkSize%JUB+1

                line = line+1

                !Last line found before end?
                if (line == LastLine) then
                    write(*,*)
                    write(*,*) 'Error reading Corners X Y'
                    stop 'ConstructGlobalVariables - HorizontalGrid - ERR220'
                end if

                !Reads XX_IE , YY_IE
                call GetData(AuxReal, Me%ObjEnterData, flag,                                &
                             Buffer_Line  = Line,                                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR230'

                if (Me%DDecomp%MasterOrSlave) then
                    if (i>= Me%DDecomp%HaloMap%ILB .and. &
                        i<= Me%DDecomp%HaloMap%IUB+1) then
                        ii = i - Me%DDecomp%HaloMap%ILB + 1
                    else
                        cycle
                    endif
                else
                    ii = i
                endif

                if (Me%DDecomp%MasterOrSlave) then
                    if (j>= Me%DDecomp%HaloMap%JLB .and. &
                        j<= Me%DDecomp%HaloMap%JUB+1) then
                        jj = j - Me%DDecomp%HaloMap%JLB + 1
                    else
                        cycle
                    endif
                else
                    jj = j
                endif

                Me%XX_IE(ii, jj) = AuxReal(1)
                Me%YY_IE(ii, jj) = AuxReal(2)

            end do
            end do


            if (Me%CoordType == CIRCULAR_ .or. Me%CornersXYInput) then

                Me%Distortion      = .true.
                Me%RegularRotation = .false.

            else

                Me%Distortion = .false.

                if (abs(Me%Grid_Angle)>0.) then
                    Me%RegularRotation = .true.
                else
                    Me%RegularRotation = .false.
                endif

            endif


            allocate(Me%DefineFacesUMap(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))
            allocate(Me%DefineFacesVMap(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))
            allocate(Me%DefineCrossMap (Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))


            Me%DefineFacesUMap(:,:) = 1
            Me%DefineFacesVMap(:,:) = 1
            Me%DefineCrossMap (:,:) = 1

            Me%NotDefinedCells       = .false.

            !If one corner has FillValue Real, the point is not to consider
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB

                if (Me%XX_IE(i  , j  ) <= FillValueReal / 2. .or.       &
                    Me%XX_IE(i  , j+1) <= FillValueReal / 2. .or.       &
                    Me%XX_IE(i+1, j+1) <= FillValueReal / 2. .or.       &
                    Me%XX_IE(i+1, j  ) <= FillValueReal / 2.) then
                    Me%DefineCellsMap(i, j)  = 0
                    Me%NotDefinedCells       = .true.
                endif

            end do
            end do


            do i = Me%WorkSize%ILB  , Me%WorkSize%IUB

                if (Me%DefineCellsMap(i, Me%WorkSize%JLB) == 0) Me%DefineFacesUMap(i, Me%WorkSize%JLB  )  = 0
                if (Me%DefineCellsMap(i, Me%WorkSize%JUB) == 0) Me%DefineFacesUMap(i, Me%WorkSize%JUB+1)  = 0

                do j = Me%WorkSize%JLB+1, Me%WorkSize%JUB


                    !if (Me%XX_IE(i  , j  ) <= FillValueReal / 2.  .or.      &
                    !    Me%XX_IE(i+1, j  ) <= FillValueReal / 2. ) then
                    !    Me%DefineFacesUMap(i,j)  = 0
                    !endif
                    if (Me%DefineCellsMap(i, j-1) == 0 .and. Me%DefineCellsMap(i, j) == 0) Me%DefineFacesUMap(i,j)  = 0

                end do
            end do


            do j = Me%WorkSize%JLB,   Me%WorkSize%JUB

                if (Me%DefineCellsMap(Me%WorkSize%ILB, j) == 0) Me%DefineFacesVMap(Me%WorkSize%ILB,j  )  = 0
                if (Me%DefineCellsMap(Me%WorkSize%IUB, j) == 0) Me%DefineFacesVMap(Me%WorkSize%IUB+1,j)  = 0

                do i = Me%WorkSize%ILB+1, Me%WorkSize%IUB

                    !if (Me%XX_IE(i  , j+1) <= FillValueReal / 2.  .or.      &
                    !    Me%XX_IE(i  , j  ) <= FillValueReal / 2. ) then
                    !    Me%DefineFacesVMap(i,j)  = 0
                    !endif
                    if (Me%DefineCellsMap(i-1, j) == 0 .and. Me%DefineCellsMap(i, j) == 0) Me%DefineFacesVMap(i,j)  = 0

                end do
            end do

            do i = Me%WorkSize%ILB+1, Me%WorkSize%IUB
            do j = Me%WorkSize%JLB+1, Me%WorkSize%JUB


                !if (Me%XX_IE(i  , j  ) <= FillValueReal / 2. ) then
                !    Me%DefineCrossMap(i, j)  = 0
                !endif
                if (Me%DefineCellsMap(i-1, j  ) == 0 .and. Me%DefineCellsMap(i, j  ) == 0 .or. &
                    Me%DefineCellsMap(i-1, j-1) == 0 .and. Me%DefineCellsMap(i, j-1) == 0) then
                    Me%DefineCrossMap(i, j)  = 0
                endif

            end do
            end do

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR240'

        else BF

            if (abs(Me%Grid_Angle)>0.) then
                Me%RegularRotation = .true.
            else
                Me%RegularRotation = .false.
            endif

           !Check if the spacing in Y is constant
            call GetData(ConstantSpacingY,                                              &
                         Me%ObjEnterData, flag,                                         &
                         SearchType   = FromFile,                                       &
                         keyword      ='CONSTANT_SPACING_Y',                            &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR250'

iconY:      if (ConstantSpacingY) Then
                !Get grid origin
                call GetData(DY,                                                        &
                             Me%ObjEnterData, flag,                                     &
                             SearchType   = FromFile,                                   &
                             keyword      ='DY',                                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR260'

                XY_Aux   = -DY
                ii       = 0
                do i = Me%GlobalWorkSize%ILB, Me%GlobalWorkSize%IUB + 1

                    XY_Aux   = XY_Aux + DY

                    if (Me%DDecomp%MasterOrSlave) then
                        if (i>= Me%DDecomp%HaloMap%ILB .and.                &
                            i<= Me%DDecomp%HaloMap%IUB+1) then
                            ii = ii + 1
                        else
                            cycle
                        endif
                    else
                        ii = i
                    endif

                    Me%YY(ii) = XY_Aux
                end do

            else iconY

                !Reads YY
                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR270'

                call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                 &
                                            BeginYY, EndYY, BlockFound,                 &
                                            FirstLine = FirstLine, LastLine = LastLine, &
                                            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR280'

                if (BlockFound) then

                    line = FirstLine
                    ii   = 0
                    do i = Me%GlobalWorkSize%ILB, Me%GlobalWorkSize%IUB+1
                        line = line+1

                        !Last line found before end?
                        if (line == LastLine) then
                            write(*,*)
                            write(*,*) 'Error reading YY'
                            stop 'ConstructGlobalVariables - HorizontalGrid - ERR290'
                        end if

                        call GetData(Aux, Me%ObjEnterData, flag,                           &
                                     Buffer_Line  = Line,                               &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR300'

                        if (Me%DDecomp%MasterOrSlave) then
                            if (i>= Me%DDecomp%HaloMap%ILB .and. &
                                i<= Me%DDecomp%HaloMap%IUB+1) then
                                ii = ii + 1
                            else
                                cycle
                            endif
                        else
                            ii = i
                        endif

                        Me%YY(ii) = Aux

                    end do
                end if

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR310'

            end if iconY

           !Check if the spacing in X is constant
            call GetData(ConstantSpacingX,                                              &
                         Me%ObjEnterData, flag,                                         &
                         SearchType   = FromFile,                                       &
                         keyword      ='CONSTANT_SPACING_X',                            &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR320'

iconX:      if (ConstantSpacingX) Then
                !Get grid spacing dx
                call GetData(DX,                                                        &
                             Me%ObjEnterData, flag,                                     &
                             SearchType   = FromFile,                                   &
                             keyword      ='DX',                                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR330'

                XY_Aux   = -DX
                jj       = 0
                do j = Me%GlobalWorkSize%JLB, Me%GlobalWorkSize%JUB + 1

                    XY_Aux = XY_Aux + DX

                    if (Me%DDecomp%MasterOrSlave) then
                        if (j>= Me%DDecomp%HaloMap%JLB .and. &
                            j<= Me%DDecomp%HaloMap%JUB+1) then
                            jj = jj + 1
                        else
                            cycle
                        endif
                    else
                        jj = j
                    endif

                    Me%XX(jj) = XY_Aux

                end do

            else iconX


                !Reads XX
                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR340'

                call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                 &
                                            BeginXX, EndXX, BlockFound,                 &
                                            FirstLine = FirstLine, LastLine = LastLine, &
                                            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR350'

                if (BlockFound) then

                    line = FirstLine

                    jj = 0
                    do j = Me%GlobalWorkSize%JLB, Me%GlobalWorkSize%JUB+1

                        line = line+1

                        !Last line found before end?
                        if (line == LastLine) then
                            write(*,*)
                            write(*,*) 'Error reading XX'
                            stop 'ConstructGlobalVariables - HorizontalGrid - ERR360'
                        end if

                        call GetData(Aux, Me%ObjEnterData, flag,                           &
                                     Buffer_Line  = Line,                               &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR370'


                        if (Me%DDecomp%MasterOrSlave) then
                            if (j>= Me%DDecomp%HaloMap%JLB .and. &
                                j<= Me%DDecomp%HaloMap%JUB+1) then
                                jj = jj + 1
                            else
                                cycle
                            endif
                        else
                            jj = j
                        endif

                        Me%XX(jj) = Aux

                    end do
                end if

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR380'


            endif iconX


        end if BF


        !Reads XX_YE, YY_IE
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR390'

        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                         &
                                    BeginCartCornersXY, EndCartCornersXY,               &
                                    Me%ReadCartCorners,                                 &
                                    FirstLine = FirstLine, LastLine = LastLine,         &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR400'

BF1:    if (Me%ReadCartCorners) then

            if (Me%CoordType /= SIMPLE_GEOG_) then
                stop 'ConstructGlobalVariables - HorizontalGrid - ERR410'
            endif


            line = FirstLine

            Me%LatitudeConn  (:, :) = Me%YY_IE(:, :)
            Me%LongitudeConn (:, :) = Me%XX_IE(:, :)



            ii = 0
            jj = 0
            do i = Me%GlobalWorkSize%ILB, Me%GlobalWorkSize%IUB+1
            do j = Me%GlobalWorkSize%JLB, Me%GlobalWorkSize%JUB+1

                line = line+1

                !Last line found before end?
                if (line == LastLine) then
                    write(*,*)
                    write(*,*) 'Error reading Cartesian Corners X Y'
                    stop 'ConstructGlobalVariables - HorizontalGrid - ERR420'
                end if

                !Reads XX_IE , YY_IE
                call GetData(AuxReal, Me%ObjEnterData, flag,                                &
                             Buffer_Line  = Line,                                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR430'

                if (Me%DDecomp%MasterOrSlave) then
                    if (i>= Me%DDecomp%HaloMap%ILB .and. &
                        i<= Me%DDecomp%HaloMap%IUB+1) then
                        ii = ii + 1
                    else
                        cycle
                    endif
                else
                    ii = i
                endif

                if (Me%DDecomp%MasterOrSlave) then
                    if (j>= Me%DDecomp%HaloMap%JLB .and. &
                        j<= Me%DDecomp%HaloMap%JUB+1) then
                        jj = jj + 1
                    else
                        cycle
                    endif
                else
                    jj = j
                endif

                Me%XX_IE(ii, jj) = AuxReal(1)
                Me%YY_IE(ii, jj) = AuxReal(2)

            end do
            end do

        else  BF1



        endif BF1
        

        !Impose boder ,limits do not compute automatically from grid
        !West, East, South, North
        call GetData(Me%BorderLimits%Values,                                            &
                     Me%ObjEnterData ,  flag,                                           &
                     keyword      = 'BORDER_LIMITS',                                    &
                     ClientModule = 'ModuleHorizontalGrid',                             &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR440'        
        
        if (flag == 4) then
            Me%BorderLimits%ON = .true.
        else
            Me%BorderLimits%ON = .false.
        endif            


        !Closes Data File
        call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - HorizontalGrid - ERR450'


        !Allocates variables common to the module
        call AllocateVariables

    end subroutine ConstructGlobalVariables

    !--------------------------------------------------------------------------


    subroutine ConstructGlobalVariablesV1(LatitudeConn, LongitudeConn, XX, YY,          &
                                          Xorig, Yorig, Latitude, Longitude,            &
                                          ILB, IUB, JLB, JUB)

        !Arguments-------------------------------------------------------------
        real, dimension(:, :), pointer          :: LatitudeConn, LongitudeConn
        real, dimension(:   ), pointer          :: XX, YY
        real                                    :: Xorig, Yorig
        real                                    :: Latitude, Longitude
        integer                                 :: ILB, IUB, JLB, JUB


        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        !Nullify T_Distances
        nullify (Me%XX )
        nullify (Me%YY )
        nullify (Me%DXX)
        nullify (Me%DYY)
        nullify (Me%DZX)
        nullify (Me%DZY)
        nullify (Me%DUX)
        nullify (Me%DUY)
        nullify (Me%DVX)
        nullify (Me%DVY)
        nullify (Me%XX_IE)
        nullify (Me%YY_IE)
        nullify (Me%Compute%XX_Z )
        nullify (Me%Compute%YY_Z )
        nullify (Me%Compute%XX_U )
        nullify (Me%Compute%YY_U )
        nullify (Me%Compute%XX_V )
        nullify (Me%Compute%YY_V )
        nullify (Me%Compute%XX_Cross)
        nullify (Me%Compute%YY_Cross)
        nullify (Me%XX_AlongGrid)
        nullify (Me%YY_AlongGrid)

        !Nullify Other
        nullify (Me%LatitudeZ    )
        nullify (Me%LongitudeZ   )
        nullify (Me%F            )
        nullify (Me%GridCellArea )
        nullify (Me%LatitudeConn )
        nullify (Me%LongitudeConn)


        Me%Size%ILB     = ILB-1
        Me%Size%IUB     = IUB+1

        Me%WorkSize%ILB = ILB
        Me%WorkSize%IUB = IUB

        Me%Size%JLB     = JLB-1
        Me%Size%JUB     = JUB+1

        Me%WorkSize%JLB = JLB
        Me%WorkSize%JUB = JUB


        !Reads COORD_TIP
        Me%CoordType = SIMPLE_GEOG_


        Me%Latitude   = Latitude
        Me%Longitude  = Longitude

        Me%Datum = WGS_84_DATUM


        Me%ProjType = PAULO_PROJECTION_

        Me%Xorig = Xorig
        Me%Yorig = Yorig

        Me%Grid_Angle = 0

        Me%Distortion      = .false.
        Me%RegularRotation = .false.
        Me%CornersXYInput  = .false.

        if (associated(LongitudeConn))                                                  &
            Me%CornersXYInput = .true.



        !Allocates XX and YY
        allocate(Me%XX(Me%Size%JLB+1 : Me%Size%JUB+1))
        allocate(Me%YY(Me%Size%ILB+1 : Me%Size%IUB+1))

        allocate(Me%XX_IE(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))
        allocate(Me%YY_IE(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))

        allocate(Me%LatitudeConn(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))
        allocate(Me%LongitudeConn(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))

        allocate(Me%DefineCellsMap(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))

        allocate(Me%XX_AlongGrid(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))
        allocate(Me%YY_AlongGrid(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))


        Me%XX    = null_real
        Me%YY    = null_real
        Me%XX_IE = null_real
        Me%YY_IE = null_real

        Me%LatitudeConn  = null_real
        Me%LongitudeConn = null_real


        Me%DefineCellsMap(:,:)  = 1

        if (Me%CornersXYInput) then

            Me%YY_IE(:,:) = LatitudeConn (:,:)
            Me%XX_IE(:,:) = LongitudeConn(:,:)

            Me%Distortion      = .true.

        else
            Me%XX(:) = XX(:)
            Me%YY(:) = YY(:)

        endif


        !Allocates variables common to the module
        call AllocateVariables

    end subroutine ConstructGlobalVariablesV1

    !--------------------------------------------------------------------------

    subroutine ConstructGlobalVariablesV2()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                 :: Imax, Jmax, STAT_CALL
        logical                 :: Exist
        !----------------------------------------------------------------------

        !Nullify T_Distances
        nullify (Me%XX )
        nullify (Me%YY )
        nullify (Me%DXX)
        nullify (Me%DYY)
        nullify (Me%DZX)
        nullify (Me%DZY)
        nullify (Me%DUX)
        nullify (Me%DUY)
        nullify (Me%DVX)
        nullify (Me%DVY)
        nullify (Me%XX_IE)
        nullify (Me%YY_IE)
        nullify (Me%Compute%XX_Z )
        nullify (Me%Compute%YY_Z )
        nullify (Me%Compute%XX_U )
        nullify (Me%Compute%YY_U )
        nullify (Me%Compute%XX_V )
        nullify (Me%Compute%YY_V )
        nullify (Me%Compute%XX_Cross)
        nullify (Me%Compute%YY_Cross)
        nullify (Me%XX_AlongGrid)
        nullify (Me%YY_AlongGrid)

        !Nullify Other
        nullify (Me%LatitudeZ    )
        nullify (Me%LongitudeZ   )
        nullify (Me%F            )
        nullify (Me%GridCellArea )
        nullify (Me%LatitudeConn )
        nullify (Me%LongitudeConn)

        call GetHDF5DataSetExist (HDF5ID        = Me%ObjHDF5,                           &
                                  DataSetName   = "/Grid/Latitude",                     &
                                  Exist         = Exist,                                &
                                  STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariablesV2 - HorizontalGrid - ERR10'

        if (Exist) then

            call GetHDF5ArrayDimensions (HDF5ID = Me%ObjHDF5, GroupName = "/Grid",      &
                                        ItemName = "Latitude", Imax = Imax, Jmax = Jmax,&
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariablesV2 - HorizontalGrid - ERR20'

            !Reads COORD_TIP
            Me%CoordType = SIMPLE_GEOG_

            Me%Datum     = WGS_84_DATUM

            Me%ProjType = PAULO_PROJECTION_
        else

            call GetHDF5DataSetExist (HDF5ID        = Me%ObjHDF5,                       &
                                      DataSetName   = "/Grid/ConnectionX",              &
                                      Exist         = Exist,                            &
                                      STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariablesV2 - HorizontalGrid - ERR30'

                if (Exist) then
                call GetHDF5ArrayDimensions (HDF5ID = Me%ObjHDF5, GroupName = "/Grid",      &
                                            ItemName = "ConnectionX", Imax = Imax, Jmax = Jmax,&
                                            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariablesV2 - HorizontalGrid - ERR40'

                !Reads COORD_TIP
                Me%CoordType = GRID_COORD_

            else
                stop 'ConstructGlobalVariablesV2 - HorizontalGrid - ERR50'
            endif


        endif

        Me%WorkSize%ILB = 1
        Me%WorkSize%IUB = Imax - 1

        Me%WorkSize%JLB = 1
        Me%WorkSize%JUB = Jmax - 1


        Me%Size%ILB     = Me%WorkSize%ILB-1
        Me%Size%IUB     = Me%WorkSize%IUB+1

        Me%Size%JLB     = Me%WorkSize%JLB-1
        Me%Size%JUB     = Me%WorkSize%JUB+1


        Me%Xorig = 0
        Me%Yorig = 0

        Me%Grid_Angle = 0

        Me%Distortion      = .false.
        Me%RegularRotation = .false.
        Me%CornersXYInput  = .false.

        Me%CornersXYInput = .true.

        !Allocates XX and YY
        allocate(Me%XX(Me%Size%JLB+1 : Me%Size%JUB+1))
        allocate(Me%YY(Me%Size%ILB+1 : Me%Size%IUB+1))

        allocate(Me%XX_IE(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))
        allocate(Me%YY_IE(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))

        allocate(Me%LatitudeConn(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))
        allocate(Me%LongitudeConn(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))

        allocate(Me%DefineCellsMap(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))

        allocate(Me%XX_AlongGrid(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))
        allocate(Me%YY_AlongGrid(Me%Size%ILB : Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))


        Me%XX    = null_real
        Me%YY    = null_real
        Me%XX_IE = null_real
        Me%YY_IE = null_real

        Me%LatitudeConn  = null_real
        Me%LongitudeConn = null_real

        Me%DefineCellsMap(:,:)  = 1

        call ReadHDF5HorizontalGrid()

        if (Me%CoordType == SIMPLE_GEOG_) then

            Me%YY_IE = Me%LatitudeConn
            Me%XX_IE = Me%LongitudeConn

            Me%Latitude   = Me%LatitudeConn (1,1)
            Me%Longitude  = Me%LongitudeConn(1,1)

        endif

        !Allocates variables common to the module
        call AllocateVariables

    end subroutine ConstructGlobalVariablesV2

    !--------------------------------------------------------------------------


    subroutine AllocateVariables()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                             :: STATUS
        integer                             :: ILB, IUB, JLB, JUB

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        !LatitudeZ
        allocate (Me%LatitudeZ(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR10'

        !LongitudeZ
        allocate (Me%LongitudeZ(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR20'

        !Allocate Distances
        allocate (Me%DXX(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR30'

        allocate (Me%DYY(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR40'

        allocate (Me%DZX(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR50'

        allocate (Me%DZY(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR60'

        allocate (Me%DUX(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR70'

        allocate (Me%DUY(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR80'

        allocate (Me%DVX(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR90'

        allocate (Me%DVY(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR100'

        !Coriolis
        allocate (Me%F(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR110'

        allocate (Me%GridCellArea(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR120'


        !Compute points location

        if (.not. Me%CornersXYInput) then

            !XX_Z
            allocate (Me%Compute%XX_Z(JLB:JUB), stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR170'

            !YY_Z
            allocate (Me%Compute%YY_Z(ILB:IUB), stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR180'

            !XX_U
            allocate (Me%Compute%XX_U(JLB:JUB), stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR190'

            !YY_V
            allocate (Me%Compute%YY_V(ILB:IUB), stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR200'

            !XX_V
            allocate (Me%Compute%XX_V(JLB:JUB), stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR210'

            !YY_U
            allocate (Me%Compute%YY_U(ILB:IUB), stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR220'

            !XX_Cross
            allocate (Me%Compute%XX_Cross(JLB:JUB), stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR230'

            !YY_Cross
            allocate (Me%Compute%YY_Cross(ILB:IUB), stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR240'

            Me%Compute%XX_Z     (:) = FillValueReal
            Me%Compute%YY_Z     (:) = FillValueReal
            Me%Compute%XX_U     (:) = FillValueReal
            Me%Compute%YY_V     (:) = FillValueReal
            Me%Compute%XX_V     (:) = FillValueReal
            Me%Compute%YY_U     (:) = FillValueReal
            Me%Compute%XX_Cross (:) = FillValueReal
            Me%Compute%YY_Cross (:) = FillValueReal

        endif

        !XX2D_Z
        allocate (Me%Compute%XX2D_Z(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR250'

        !YY2D_Z
        allocate (Me%Compute%YY2D_Z(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR260'

        !XX2D_U
        allocate (Me%Compute%XX2D_U(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR270'

        !YY2D_V
        allocate (Me%Compute%YY2D_V(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR280'

        !XX2D_V
        allocate (Me%Compute%XX2D_V(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR290'

        !YY2D_U
        allocate (Me%Compute%YY2D_U(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR300'

        Me%Compute%XX2D_Z     (:,:) = FillValueReal
        Me%Compute%YY2D_Z     (:,:) = FillValueReal
        Me%Compute%XX2D_U     (:,:) = FillValueReal
        Me%Compute%YY2D_V     (:,:) = FillValueReal
        Me%Compute%XX2D_V     (:,:) = FillValueReal
        Me%Compute%YY2D_U     (:,:) = FillValueReal


        !Initialize Values
        Me%LatitudeZ    (:,:)    = FillValueReal
        Me%LongitudeZ   (:,:)    = FillValueReal
        Me%DXX          (:,:)    = FillValueReal
        Me%DYY          (:,:)    = FillValueReal
        Me%DZX          (:,:)    = FillValueReal
        Me%DZY          (:,:)    = FillValueReal
        Me%DUX          (:,:)    = FillValueReal
        Me%DUY          (:,:)    = FillValueReal
        Me%DVX          (:,:)    = FillValueReal
        Me%DVY          (:,:)    = FillValueReal
        Me%F            (:,:)    = FillValueReal
        Me%GridCellArea (:,:)    = FillValueReal

        if (Me%CoordType == CIRCULAR_ .or. Me%CornersXYInput) then

            Me%Distortion      = .true.

            allocate (Me%RotationX(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR310'

            Me%RotationX(:,:) = FillValueReal

            allocate (Me%RotationY(ILB:IUB, JLB:JUB), stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR320'

            Me%RotationY(:,:) = FillValueReal

        endif

        allocate (Me%XX1D_Aux(JLB:JUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR330'


        allocate (Me%YY1D_Aux(ILB:IUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR330'

        allocate(Me%AuxPolygon)
        Me%AuxPolygon%Count = 5
        allocate(Me%AuxPolygon%VerticesF(1:5))

    end subroutine AllocateVariables


    !--------------------------------------------------------------------------
    !This subroutine needs to be reviewed for the case that the type of the coordinates
    !is equal to one (Geographic coordinates) - Frank Abr 99.

    subroutine Mercator()
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
#ifdef _USE_PROJ4
        real(8), dimension(:,:), pointer    :: XX_Pointer, YY_Pointer
#endif
        real(8), dimension(:,:), allocatable, target :: XX_aux, YY_aux
        real, dimension(:,:), pointer       :: XX_IE, YY_IE
        real, dimension(:  ), pointer       :: XX, YY
        integer                             :: AuxCoordTip
        real(8)                             :: DB_LAT, DB_LONG
        real                                :: X_PONTO, Y_PONTO
        real(8), dimension(:), pointer      :: DLZONE
        integer, dimension(:), pointer      :: IZONE1
        integer                             :: ILB, IUB, JLB, JUB
        integer                             :: i, j, STATUS, k
        integer, dimension(2)               :: GridUTMmax, GridUTMmin, GridUTM

        real(8)                             :: radians, Radius, Angle
        !real(8)                             :: EarthRadius, Rad_Lat, CosenLat

        real                                :: MaxLon, MinLon, MaxLat, MinLat, MinX, MinY


        !Worksize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        XX_IE => Me%XX_IE
        YY_IE => Me%YY_IE

        XX    => Me%XX
        YY    => Me%YY

        !By default
do5 :   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do6 :   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
           Me%LatitudeZ (i, j) = Me%Latitude
           Me%LongitudeZ(i, j) = Me%Longitude
        end do do6
        end do do5


        !Projected coordenates
        !Calculates the coordinates in the lower left corner (IE)

cd1:    if (Me%CoordType == UTM_            .or.    &
            Me%CoordType == MIL_PORT_       .or.    &
            Me%CoordType == GRID_COORD_     .or.    &
            Me%CoordType == NLRD_           .or.    &
            Me%CoordType == SIMPLE_GEOG_) then

Inp:        if (.not. Me%CornersXYInput) then

do1 :           do j = JLB, JUB + 1
do2 :           do i = ILB, IUB + 1
                    XX_IE(I, J) = XX(J)
                    YY_IE(I, J) = YY(I)
                end do do2
                end do do1


                !Rotation of the points
do3 :           do j = JLB, JUB + 1
do4 :           do i = ILB, IUB + 1
                    call RODAXY(Me%Xorig, Me%Yorig,            &
                                Me%Grid_Angle, XX_IE(i, j), YY_IE(i, j))
                end do do4
                end do do3

            endif Inp

        endif cd1

            !Computes the Geographic Coordinates in the center od the cells (to use later,
            !when computing Coriolis and the elevation of the sun
cd2:    if ((Me%CoordType == UTM_) .OR. (Me%CoordType == MIL_PORT_)) then

do7 :       do j = JLB, JUB
do8 :       do i = ILB, IUB

                !X_PONTO = (XX_IE(i  , j+1) + XX_IE(i, j)) / 2.
                !Y_PONTO = (YY_IE(i+1, j  ) + YY_IE(i, j)) / 2.

                if (Me%DefineCellsMap(i, j) == 1) then

                    Y_PONTO  = (YY_IE(i  , j+1) + YY_IE(i, j) + YY_IE(i+1, j+1) + YY_IE(i+1, j)) / 4.
                    X_PONTO  = (XX_IE(i  , j+1) + XX_IE(i, j) + XX_IE(i+1, j+1) + XX_IE(i+1, j)) / 4.
                else
                    Y_PONTO  = 0
                    X_PONTO  = 0
                endif


                if      (Me%CoordType == MIL_PORT_) then

                    !from Portuguese coord. to geographic coord.
                    AuxCoordTip = 5

                    call USCONVCO (AuxCoordTip, Me%ZoneLong, DB_LAT, DB_LONG,        &
                                   dble(X_PONTO), dble(Y_PONTO))

                else if (Me%CoordType == UTM_) then

                    !DATUM - WGS84
                    call UTMToLatLon (dble(X_PONTO), dble(Y_PONTO),                     &
                                      DB_LONG, DB_LAT, Me%Grid_Zone, WGS_84_DATUM)
                endif

                Me%LatitudeZ(i, j)  = real(DB_LAT)
                Me%LongitudeZ(i, j) = real(DB_LONG)

            end do do8
            end do do7

        endif cd2



cd21:   if(Me%CoordType == NLRD_) then

            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%DefineCellsMap(i, j) == 1) then
                    Y_PONTO  = (YY_IE(i  , j+1) + YY_IE(i, j) + YY_IE(i+1, j+1) + YY_IE(i+1, j)) / 4.
                    X_PONTO  = (XX_IE(i  , j+1) + XX_IE(i, j) + XX_IE(i+1, j+1) + XX_IE(i+1, j)) / 4.
                else
                    Y_PONTO  = 0.
                    X_PONTO  = 0.
                endif

                call rd_to_wgs84(dble(X_PONTO), dble(Y_PONTO), DB_LAT, DB_LONG)


                Me%LatitudeZ(i, j)  = real(DB_LAT)
                Me%LongitudeZ(i, j) = real(DB_LONG)

            end do
            end do

        endif cd21

        !Geographic coordinates
cd22:   if (Me%CoordType == GEOG_               .or.  &
            Me%CoordType == SIMPLE_GEOG_) then


            !Computes the Geographic Coordinates in the center od the cells (to use later,
            !when computing Coriolis and the elevation of the sun

            if (.not. Me%CornersXYInput .and. Me%CoordType == GEOG_) then

do9  :          do j = JLB, JUB
do10 :          do i = ILB, IUB
                   Me%LatitudeZ (i, j) = Me%Yorig + (YY(i+1) + YY(i)) / 2.0
                   Me%LongitudeZ(i, j) = Me%Xorig + (XX(j+1) + XX(j)) / 2.0
                end do do10
                end do do9

            else

do19 :          do j = JLB, JUB
do20 :          do i = ILB, IUB
                   Me%LatitudeZ (i, j) = (YY_IE(i  , j+1) + YY_IE(i, j) + YY_IE(i+1, j+1) + YY_IE(i+1, j)) / 4.
                   Me%LongitudeZ(i, j) = (XX_IE(i  , j+1) + XX_IE(i, j) + XX_IE(i+1, j+1) + XX_IE(i+1, j)) / 4.
                end do do20
                end do do19

            endif

cd36:       if (Me%CoordType == GEOG_) then

Inp3:           if (.not. Me%CornersXYInput) then

cd33 :              if (Me%Grid_Angle .NE. 0.0) then

                        write(*,*) 'The Geographic coordinates can not have the bathymetry rotated'
                        write(*,*)'Please select SIMPLE GEOGRAPHIC coordinates.'
                        stop  'Subroutine Mercator; ModuleHorizontalGrid. ERR35'

                    endif cd33

                    !Convertes geodesic coordenates into UTM
                    nullify (DLZONE, IZONE1)
                    allocate(DLZONE(Me%Size%ILB:Me%Size%IUB))
                    allocate(IZONE1(Me%Size%JLB:Me%Size%JUB))

                    call ConvertCoordenates(XX, YY, Me%Xorig, Me%Yorig, IZONE1, DLZONE)

                    radians = 180.0 / PI
                    do i = ILB, IUB + 1
                        XX_IE(i, 1) = 0.
                        YY_IE(i, 1) = (Me%Yorig + YY(i)) * 6378000. / radians
                        do j = JLB+1, JUB+1
                            XX_IE(i,j) = abs((XX(j) - XX(1))* DLZONE(i)) / 6.
                            YY_IE(i,j) = YY_IE(I,1)
                        enddo
                    enddo

                    deallocate(DLZONE, stat = STATUS)
                        if (STATUS /= SUCCESS_) stop 'Mercator - HorizontalGrid - ERR13'
                    deallocate(IZONE1, stat = STATUS)
                        if (STATUS /= SUCCESS_) stop 'Mercator - HorizontalGrid - ERR14'
                    nullify(DLZONE,IZONE1)

                else

                    write(*,*)'Cannot define curvilinear grid in GEOGRAPHIC coordinates.'
                    write(*,*)'Please select SIMPLE GEOGRAPHIC coordinates.'
                    stop 'Mercator - HorizontalGrid - ERR14a'

                end if Inp3


            else if (Me%CoordType == SIMPLE_GEOG_ .and. .not. Me%ReadCartCorners) then cd36

                allocate (XX_aux(ILB-1:IUB+1,JLB-1:JUB+1))
                allocate (YY_aux(ILB-1:IUB+1,JLB-1:JUB+1))

                XX_aux(:,:) = FillValueReal
                YY_aux(:,:) = FillValueReal

                Me%LatitudeConn  (:, :) = YY_IE(:, :)
                Me%LongitudeConn (:, :) = XX_IE(:, :)

ipp:            if (Me%ProjType == PAULO_PROJECTION_) then



                    !radians      = Pi / 180.0
                    !EarthRadius  = 6378000.
                    !
                    do j = JLB, JUB + 1
                    do i = ILB, IUB + 1
                    !    Rad_Lat     = Me%LatitudeConn(i,j)* radians
                    !    CosenLat    = cos(Rad_Lat)
                    !    XX_IE(i, j) = CosenLat * EarthRadius * (Me%LongitudeConn(i, j) - Me%Longitude) * radians
                    !    YY_IE(i, j) =            EarthRadius * (Me%LatitudeConn (i, j) - Me%Latitude ) * radians
                        call FromSphericalToCart(Lat = Me%LatitudeConn (i, j),              &
                                                 Lon = Me%LongitudeConn(i, j),              &
                                                 X   = XX_IE(i, j),                         &
                                                 Y   = YY_IE(i, j))    
                    enddo
                    enddo

                else ipp

#ifdef _USE_PROJ4
ifU:                if (Me%ProjType < 0) then
#else
ifU:                if (.not. Me%UseLambert) then
#endif

                        GridUTMmax(1) = -999
                        GridUTMmax(2) = -999
                        GridUTMmin(1) = +999
                        GridUTMmin(2) = +999

                        do j = JLB, JUB + 1
                        do i = ILB, IUB + 1

                            if (XX_IE(i,j) > FillValueReal/3.) then

                                !DATUM - WGS84
                                call LatLonToUTM (dble(XX_IE(i,j)), dble(YY_IE(i,j)),   &
                                                  DB_LONG, DB_LAT, GridUTM, WGS_84_DATUM)

                                YY_aux(i,j) = real(DB_LAT)
                                XX_aux(i,j) = real(DB_LONG)

                                do k=1,2
                                    if (GridUTM(k) > GridUTMmax(k)) GridUTMmax(k) = GridUTM(k)
                                    if (GridUTM(k) < GridUTMmin(k)) GridUTMmin(k) = GridUTM(k)
                                enddo

                            endif
                        enddo
                        enddo


                        do k=1,2
                            if (GridUTMmax(k) /= GridUTMmin(k)) then

                                write(*,*)'The grid is not inside one UTM cell'
!                               write(*,*)'Cannot define curvilinear grid in SIMPLE GEOGRAPHIC coordinates.'
!                               write(*,*)'The grid must be fully inside one UTM cell'
!                               stop 'Mercator - HorizontalGrid - ERR140'

                                !find central reference coordinates

                                MaxLon = null_real
                                MaxLat = null_real

                                MinLon = - null_real
                                MinLat = - null_real

                                do j = JLB, JUB + 1
                                do i = ILB, IUB + 1

                                    if (XX_IE(i,j) .LT. MinLon) MinLon = XX_IE(i,j)
                                    if (XX_IE(i,j) .GT. MaxLon) MaxLon = XX_IE(i,j)

                                    if (YY_IE(i,j) .LT. MinLat) MinLat = YY_IE(i,j)
                                    if (YY_IE(i,j) .GT. MaxLat) MaxLat = YY_IE(i,j)

                                enddo
                                enddo

                                write(*,*)'Max, Min Lon', MaxLon, MinLon
                                write(*,*)'Max, Min Lat', MaxLat, MinLat

                                Me%Longitude = (MaxLon + MinLon) / 2.
                                Me%Latitude  = (MaxLat + MinLat) / 2.

                                Me%Easting  = null_real
                                Me%Northing = null_real


#ifdef _USE_PROJ4
                                Me%ProjType = SPHERE_MERCATOR_
                                Me%SP1      = 0.0
#else
                                Me%UseLambert = .true.

#endif

                            endif
                        enddo

                    endif ifU
#ifdef _USE_PROJ4
ifP:                if (Me%ProjType == LAMB_CONF_CONIC_) then
#else
ifP:                if (Me%UseLambert) then
#endif

                        write(*,*)'Using Lambert Conformal Conic Projection'
                        write(*,*)'with parameters :'
                        write(*,*)'SP1                  = ', Me%SP1
                        write(*,*)'SP2                  = ', Me%SP2
                        write(*,*)'Reference Longitude  = ', Me%Longitude
                        write(*,*)'Reference Latitude   = ', Me%Latitude
                        write(*,*)

                        do j = JLB, JUB + 1
                        do i = ILB, IUB + 1
                            if (XX_IE(i,j) > FillValueReal/3.) then

                                call LatLonToLambertSP2 (dble(YY_IE(i,j)),dble(XX_IE(i,j)),      &
                                                         dble(Me%Latitude), dble(Me%Longitude),  &
                                                         dble(Me%SP1), dble(Me%SP2),             &
                                                         Me%Datum, DB_LONG, DB_LAT)

                                YY_aux(i,j) = real(DB_LAT)
                                XX_aux(i,j) = real(DB_LONG)

                            endif
                        enddo
                        enddo
#ifdef _USE_PROJ4
                    elseif (Me%ProjType == SPHERE_MERCATOR_) then  !ifP

                        XX_Pointer => XX_aux
                        YY_Pointer => YY_aux

                        call FromGeographic2SphericMercator(XX_IE, YY_IE, Me%WorkSize, XX_Pointer, YY_Pointer)

                        nullify (XX_Pointer, YY_Pointer)
#endif
                    endif ifP


#ifdef _USE_PROJ4
                if (Me%ProjType > 0  .and. Me%Easting .LT. 0.0) then
#else
                if (Me%UseLambert .and. Me%Easting .LT. 0.0) then
#endif

                    MinX = - null_real
                    MinY = - null_real

                    do j = JLB, JUB + 1
                    do i = ILB, IUB + 1
                        if (XX_aux(i,j) .LT. MinX) MinX = XX_aux(i,j)
                        if (XX_aux(i,j) .LT. MinY) MinY = XX_aux(i,j)
                    enddo
                    enddo

                    Me%Easting  = - MinX
                    Me%Northing = - MinY
                end if

                XX_aux = XX_aux + Me%Easting
                YY_aux = YY_aux + Me%Northing

                XX_IE(:, :) = XX_aux(:, :)
                YY_IE(:, :) = YY_aux(:, :)

                deallocate (XX_aux)
                deallocate (YY_aux)

           endif ipp
           endif cd36

        endif cd22

        !Circular coordinates
cd23:   if (Me%CoordType == CIRCULAR_) then

            radians      = PI / 180.0


            do  i = ILB, IUB + 1

                Angle  = Me%Yorig * radians + YY(i)* radians

                do  j = JLB, JUB + 1

                    Radius      = Me%Xorig + XX(j)

                    XX_IE(i, j) = Radius * cos (Angle)
                    YY_IE(i, j) = Radius * sin (Angle)

                enddo
            enddo

        endif cd23



        !In the case that the model uses Coordinates different then Grid Coordinates
        !or Circular coordinates write them to file
        if (Me%CoordType == GEOG_) then

            if (.not. Me%CornersXYInput) then

                do j = JLB, JUB + 1
                do i = ILB, IUB + 1
                   Me%LatitudeConn  (i, j) = Me%YORIG + Me%YY(i)
                   Me%LongitudeConn (i, j) = Me%XORIG + Me%XX(j)
                end do
                end do

            end if

        elseif (Me%CoordType == UTM_ .or. Me%CoordType == MIL_PORT_) then

            do j = JLB, JUB + 1
            do i = ILB, IUB + 1

                X_PONTO = Me%XX_IE(i, j)
                Y_PONTO = Me%YY_IE(i, j)

                if (X_PONTO < FillValueReal/2. .or. Y_PONTO < FillValueReal/2.) then
                    X_PONTO = 0.
                    Y_PONTO = 0.
                endif

                if (Me%CoordType == MIL_PORT_) then

                    AuxCoordTip = 5 !from Portuguese coord. to geographic coord.

                    call USCONVCO (AuxCoordTip, Me%ZoneLong, DB_LAT, DB_LONG, &
                                   dble(X_PONTO), dble(Y_PONTO))

                else if (Me%CoordType == UTM_) then

                    !DATUM - WGS84
                    call UTMToLatLon (dble(X_PONTO), dble(Y_PONTO),                     &
                                   DB_LONG, DB_LAT, Me%Grid_Zone,  WGS_84_DATUM)
                endif

                Me%LatitudeConn (i, j) = real(DB_LAT)
                Me%LongitudeConn(i, j) = real(DB_LONG)


            end do
            end do

        elseif(Me%CoordType == NLRD_)then

            do j = JLB, JUB + 1
            do i = ILB, IUB + 1


                X_PONTO = Me%XX_IE(i, j)
                Y_PONTO = Me%YY_IE(i, j)

                if (X_PONTO == FillValueReal .or. Y_PONTO == FillValueReal) then
                    X_PONTO = 0.
                    Y_PONTO = 0.
                endif

                call rd_to_wgs84(dble(X_PONTO), dble(Y_PONTO), DB_LAT, DB_LONG)

                Me%LatitudeConn (i, j) = real(DB_LAT)
                Me%LongitudeConn(i, j) = real(DB_LONG)

            enddo
            enddo

        endif

        nullify(XX, YY)
        nullify(XX_IE, YY_IE)


    end subroutine Mercator

    !--------------------------------------------------------------------------
    !This subroutine convert spherical coordinates in to cartesina coordinates

    subroutine FromSphericalToCart(Lat, Lon, X, Y)
    
        !Arguments-------------------------------------------------------------
        real(8), intent(IN)             :: Lat, Lon    
        real(8), intent(Out)            :: X, Y

        !Local-----------------------------------------------------------------
        !real(8)                         :: radians, EarthRadius, Rad_Lat, CosenLat

        !Begin-----------------------------------------------------------------
                
        !radians      = Pi / 180.0
        !EarthRadius  = 6378000.
        !
        !Rad_Lat     = Lat * radians
        !CosenLat    = cos(Rad_Lat)
        !X           = CosenLat * EarthRadius * (Lon - Me%Longitude) * radians
        !Y           =            EarthRadius * (Lat - Me%Latitude ) * radians             

        call SphericalToCart(Lat, Lon, X, Y, Me%Longitude, Me%Latitude)        

    end subroutine FromSphericalToCart
    
    
                
    subroutine FromCartToSpherical(X, Y, Lat, Lon)
    
        !Arguments-------------------------------------------------------------
        real(8), intent(IN)             :: X, Y   
        real(8), intent(Out)            :: Lat, Lon

        !Local-----------------------------------------------------------------
        real(8)                         :: radians, EarthRadius, Rad_Lat, CosenLat

        !Begin-----------------------------------------------------------------
                
        radians      = Pi / 180.0
        EarthRadius  = 6378000.
        
        Lat          = Y / (EarthRadius * radians) + Me%Latitude
        
        Rad_Lat      = Lat * radians
        CosenLat     = cos(Rad_Lat) 

        if (CosenLat == 0.) then
            stop 'FromCartToSpherical - ModuleHorizontalGrid - ERR10'
        endif    
        
        Lon          = X / (CosenLat * EarthRadius * radians) + Me%Longitude
        
        
    end subroutine FromCartToSpherical    

#ifdef _USE_PROJ4

    subroutine FromGeographic2SphericMercator(XX_IE, YY_IE, WorkSize, XX_aux, YY_aux)

        use proj4

        !Arguments-------------------------------------------------------------
        real,    dimension(:,:), pointer         :: XX_IE, YY_IE
        real(8), dimension(:,:), pointer         :: XX_aux, YY_aux
        type (T_Size2D)                          :: WorkSize

        !Local-----------------------------------------------------------------
        real(8)                                  :: DB_LAT, DB_LONG
        integer                                  :: ILB, IUB, JLB, JUB
        integer                                  :: status, i, j

        character(len=20), dimension(:), pointer :: params
        type(prj90_projection)                   :: proj


        !Begin-----------------------------------------------------------------

        allocate(params(10))
        params(1) = 'proj=merc'
        params(2) = 'lat_ts=0.0'
        params(3) = 'lon_0=0.0'
        params(4) = 'k=1.0'
        params(5) = 'x_0=0.0'
        params(6) = 'y_0=0.0'
        params(7) = 'a=6371000'
        params(8) = 'b=6371000'
        params(9) = 'datum=WGS84'
        params(10)= 'units=m'

        ILB = WorkSize%ILB
        IUB = WorkSize%IUB
        JLB = WorkSize%JLB
        JUB = WorkSize%JUB


        status=prj90_init(proj,params)
        if (status.ne.PRJ90_NOERR) then
            write(*,*) prj90_strerrno(status)
            stop 'FromGeographic2SphericMercator - ModuleHorizontalGrid - ERR10'
        endif

        do j = JLB, JUB + 1
        do i = ILB, IUB + 1
            if (XX_IE(i,j) > FillValueReal/3.) then

                status = prj90_fwd(proj,dble(XX_IE(i,j)),dble(YY_IE(i,j)),DB_LONG, DB_LAT)
                if (status.ne.PRJ90_NOERR) then
                    write(*,*) prj90_strerrno(status)
                    stop 'FromGeographic2SphericMercator - ModuleHorizontalGrid - ERR20'
                end if


                YY_aux(i,j) = real(DB_LAT)
                XX_aux(i,j) = real(DB_LONG)

            endif
        enddo
        enddo

        status = prj90_free(proj)
        if (status.ne.PRJ90_NOERR) then
            write(*,*) prj90_strerrno(status)
            stop 'FromGeographic2SphericMercator - ModuleHorizontalGrid - ERR30'
        end if

        deallocate(params)


    end subroutine FromGeographic2SphericMercator

    subroutine FromGeo2SpherMercator1D(X1D_Geo, Y1D_Geo, ILB, IUB, X1D_Out, Y1D_Out)

        use proj4

        !Arguments-------------------------------------------------------------
        real(8), dimension(:  ), pointer         :: X1D_Geo, Y1D_Geo
        real(8), dimension(:  ), pointer         :: X1D_Out, Y1D_Out
        integer                                  :: ILB, IUB

        !Local-----------------------------------------------------------------
        real(8)                                  :: DB_LAT, DB_LONG
        integer                                  :: status, i

        character(len=20), dimension(:), pointer :: params
        type(prj90_projection)                   :: proj


        !Begin-----------------------------------------------------------------

        allocate(params(8))
        params(1) = 'proj=merc'
        params(2) = 'lat_ts=0.0'
        params(3) = 'lon_0=0.0'
        params(4) = 'k=1.0'
        params(5) = 'x_0=0.0'
        params(6) = 'y_0=0.0'
        params(7) = 'a=6371000'
        params(8) = 'b=6371000'


        status=prj90_init(proj,params)
        if (status.ne.PRJ90_NOERR) then
            write(*,*) prj90_strerrno(status)
            stop 'FromGeo2SpherMercator1D - ModuleHorizontalGrid - ERR10'
        endif

        do i = ILB, IUB
            status = prj90_fwd(proj,dble(X1D_Geo(i)),dble(Y1D_Geo(i)),DB_LONG, DB_LAT)
            if (status.ne.PRJ90_NOERR) then
                write(*,*) prj90_strerrno(status)
                stop 'FromGeo2SpherMercator1D - ModuleHorizontalGrid - ERR20'
            end if


            Y1D_Out(i) = real(DB_LAT)
            X1D_Out(i) = real(DB_LONG)
        enddo

        status = prj90_free(proj)
        if (status.ne.PRJ90_NOERR) then
            write(*,*) prj90_strerrno(status)
            stop 'FromGeo2SpherMercator1D - ModuleHorizontalGrid - ERR30'
        end if

        deallocate(params)


    end subroutine FromGeo2SpherMercator1D

    subroutine FromGeo2SpherMercatorScalar(X_Geo, Y_Geo, X_Out, Y_Out)

        use proj4

        !Arguments-------------------------------------------------------------
        real(8)                                  :: X_Geo, Y_Geo, X_Out, Y_Out

        !Local-----------------------------------------------------------------
        real(8)                                  :: DB_LAT, DB_LONG
        integer                                  :: status

        character(len=20), dimension(:), pointer :: params
        type(prj90_projection)                   :: proj


        !Begin-----------------------------------------------------------------

        allocate(params(8))
        params(1) = 'proj=merc'
        params(2) = 'lat_ts=0.0'
        params(3) = 'lon_0=0.0'
        params(4) = 'k=1.0'
        params(5) = 'x_0=0.0'
        params(6) = 'y_0=0.0'
        params(7) = 'a=6371000'
        params(8) = 'b=6371000'


        status=prj90_init(proj,params)
        if (status.ne.PRJ90_NOERR) then
            write(*,*) prj90_strerrno(status)
            stop 'FromGeo2SpherMercatorScalar - ModuleHorizontalGrid - ERR10'
        endif

        status = prj90_fwd(proj,dble(X_Geo),dble(Y_Geo),DB_LONG, DB_LAT)
        if (status.ne.PRJ90_NOERR) then
            write(*,*) prj90_strerrno(status)
            stop 'FromGeo2SpherMercatorScalar - ModuleHorizontalGrid - ERR20'
        end if


        Y_Out = real(DB_LAT)
        X_Out = real(DB_LONG)


        status = prj90_free(proj)
        if (status.ne.PRJ90_NOERR) then
            write(*,*) prj90_strerrno(status)
            stop 'FromGeo2SpherMercatorScalar - ModuleHorizontalGrid - ERR30'
        end if

        deallocate(params)


    end subroutine FromGeo2SpherMercatorScalar


#endif

    subroutine CheckGridBorder()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        nullify (Me%GridBorderCart     )
        nullify (Me%GridBorderCoord    )
        nullify (Me%GridOutBorderCoord )
        nullify (Me%GridBorderAlongGrid)
        nullify (Me%GridOutBorderCart  )


        allocate(Me%GridBorderCart     )
        allocate(Me%GridBorderCoord    )
        allocate(Me%GridOutBorderCoord )
        allocate(Me%GridBorderAlongGrid)
        allocate(Me%GridOutBorderCart  )


        Me%GridBorderCart%Type_         = Rectang_
        Me%GridBorderCoord%Type_        = Rectang_
        Me%GridOutBorderCoord%Type_     = Rectang_
        Me%GridBorderAlongGrid%Type_    = Rectang_
        Me%GridOutBorderCart%Type_      = Rectang_

Inp:    if (Me%CornersXYInput) then

            Me%GridBorderCart%Type_      = ComplexPolygon_
            Me%GridBorderCoord%Type_     = ComplexPolygon_
            Me%GridOutBorderCoord%Type_  = ComplexPolygon_
            Me%GridBorderAlongGrid%Type_ = ComplexPolygon_
            Me%GridOutBorderCart%Type_   = ComplexPolygon_

        else  Inp

            if (abs(Me%Grid_Angle) > 0.) then

                Me%GridBorderCart%Type_      = RotatedRectang_
                Me%GridBorderCoord%Type_     = RotatedRectang_
                Me%GridOutBorderCoord%Type_  = RotatedRectang_
                Me%GridBorderAlongGrid%Type_ = RotatedRectang_
                Me%GridOutBorderCart%Type_   = RotatedRectang_


            endif


        endif Inp

        if (Me%CoordType == SIMPLE_GEOG_ .or. Me%CoordType == GEOG_) then

            Me%GridBorderCart%Type_      = ComplexPolygon_
            Me%GridBorderAlongGrid%Type_ = ComplexPolygon_
            Me%GridOutBorderCart%Type_   = ComplexPolygon_
        endif


    end subroutine CheckGridBorder


    !--------------------------------------------------------------------------
    !This subroutine defines the border polygon

    subroutine DefineBorderPolygons()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,   dimension(:,:), pointer             :: XX2D, YY2D

        !Begin-----------------------------------------------------------------


        if (Me%CoordType == SIMPLE_GEOG_ .or. Me%CoordType == GEOG_) then

            XX2D        => Me%LongitudeConn
            YY2D        => Me%LatitudeConn

        else

            XX2D        => Me%XX_IE
            YY2D        => Me%YY_IE

        endif

        call DefinesBorderPoly(XX2D, YY2D, Me%GridBorderCoord)

        call DefinesBorderPoly(XX2D, YY2D, Me%GridOutBorderCoord, Outer = .true.)

        nullify(XX2D, YY2D)

        call DefinesBorderPoly(Me%XX_IE, Me%YY_IE, Me%GridBorderCart)

        call DefinesBorderPoly(Me%XX_IE, Me%YY_IE, Me%GridOutBorderCart, Outer = .true.)

        call DefinesBorderPoly(Me%XX_AlongGrid, Me%YY_AlongGrid, Me%GridBorderAlongGrid)

    end subroutine DefineBorderPolygons
    !--------------------------------------------------------------------------
    !This subroutine defines the border polygon


    subroutine DefinesBorderPoly(XX2D, YY2D, GridBorder, Outer)

        !Arguments-------------------------------------------------------------
        type (T_Border),  pointer                   :: GridBorder
        real,   dimension(:,:), pointer             :: XX2D, YY2D
        logical, optional, intent(IN)               :: Outer

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: Nvert, i, j, di
        logical                                     :: Outer_

        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        if (present(Outer)) then
            Outer_ = Outer
        else
            Outer_ = .false.
        endif


        if      (GridBorder%Type_ == ComplexPolygon_) then

            Nvert = 0

            Nvert = 2*(IUB-ILB+1)

            if (Outer_) then
                Nvert = Nvert + 4
            endif

            Nvert = Nvert + 2*(JUB-JLB+1)

            if (Outer_) then
                Nvert = Nvert + 4
            endif

        else

            Nvert = 4

        endif

        !The last vertix equal to the first
        Nvert = Nvert + 1



        allocate(GridBorder%Polygon_)
        allocate(GridBorder%Polygon_%VerticesF(1:Nvert))

        if (Outer_) then
            di = 0
        else
            di = 1
        endif

        if      (GridBorder%Type_ == ComplexPolygon_) then

            GridBorder%Polygon_%Count = 0

            !West boundary
            do i = ILB + di, IUB-di
                GridBorder%Polygon_%Count = GridBorder%Polygon_%Count + 1
                GridBorder%Polygon_%VerticesF(GridBorder%Polygon_%Count)%X = XX2D(i, JLB+di)
                GridBorder%Polygon_%VerticesF(GridBorder%Polygon_%Count)%Y = YY2D(i, JLB+di)
            enddo

            !North boundary
            do j = JLB+di, JUB-di + 1
                GridBorder%Polygon_%Count = GridBorder%Polygon_%Count + 1
                GridBorder%Polygon_%VerticesF(GridBorder%Polygon_%Count)%X = XX2D(IUB+1-di, j)
                GridBorder%Polygon_%VerticesF(GridBorder%Polygon_%Count)%Y = YY2D(IUB+1-di, j)
            enddo

            !East boundary
            do i =  IUB-di, ILB + di, -1
                GridBorder%Polygon_%Count = GridBorder%Polygon_%Count + 1
                GridBorder%Polygon_%VerticesF(GridBorder%Polygon_%Count)%X = XX2D(i, JUB+1-di)
                GridBorder%Polygon_%VerticesF(GridBorder%Polygon_%Count)%Y = YY2D(i, JUB+1-di)
            enddo


            !South boundary
            do j =  JUB-di, JLB+di, -1
                GridBorder%Polygon_%Count = GridBorder%Polygon_%Count + 1
                GridBorder%Polygon_%VerticesF(GridBorder%Polygon_%Count)%X = XX2D(ILB+di, j)
                GridBorder%Polygon_%VerticesF(GridBorder%Polygon_%Count)%Y = YY2D(ILB+di, j)
            enddo

!        else if (GridBorder%Type_ == RotatedRectang_) then

        else

            GridBorder%Polygon_%Count = Nvert

            GridBorder%Polygon_%VerticesF(1)%X = XX2D(ILB+di,     JLB+di)
            GridBorder%Polygon_%VerticesF(1)%Y = YY2D(ILB+di,     JLB+di)


            GridBorder%Polygon_%VerticesF(2)%X = XX2D(ILB+di,     JUB+1-di)
            GridBorder%Polygon_%VerticesF(2)%Y = YY2D(ILB+di,     JUB+1-di)

            GridBorder%Polygon_%VerticesF(3)%X = XX2D(IUB+1-di,   JUB+1-di)
            GridBorder%Polygon_%VerticesF(3)%Y = YY2D(IUB+1-di,   JUB+1-di)

            GridBorder%Polygon_%VerticesF(4)%X = XX2D(IUB+1-di,   JLB+di)
            GridBorder%Polygon_%VerticesF(4)%Y = YY2D(IUB+1-di,   JLB+di)

            !Last vertex equal to first
            GridBorder%Polygon_%VerticesF(5)   = GridBorder%Polygon_%VerticesF(1)

        endif

        call SetLimits(GridBorder%Polygon_)



    end subroutine DefinesBorderPoly

    !--------------------------------------------------------------------------

    !***********************************************************************
    !                                                                      *
    !  Esta subroutina calcula 2 vectores e um escalar para a conversão    *
    !  de coordenadas geodésicas em UTM:                                   *
    !                                                                      *
    !        comprimento do fuso:                    DLZONE(I)             *
    !        nº do fuso:                             IZONE1(J)             *
    !        nº do fuso de origem:                   IZONE_ORIG            *
    !                                                                      *
    !  O meridiano de Greenwich está no fuso 31. Considera-se como origem  *
    !  o ínicio do 1º fuso da malha. Cada fuso tem 6 graus e há 60 fusos.  *
    !                                                                      *
    !                                                                      *
    !    AIRES: 30/6/1977                                                  *
    !    Pedro Montero 20/11/97   Se pasan las variables por argumentos    *
    !***********************************************************************

    subroutine ConvertCoordenates(XX, YY, XORIG, YORIG, IZONE1, DLZONE)

        !Arguments-------------------------------------------------------------
        real, dimension(:), pointer         :: XX, YY
        real                                :: XORIG, YORIG
        integer, dimension(:), pointer      :: IZONE1
        real(8), dimension(:), pointer      :: DLZONE

        !Local-----------------------------------------------------------------

        integer :: i, j
        integer :: ILB, IUB, JLB, JUB
        integer :: IZONE_ORIG

        real    :: radians, DB_LAT, DB_LONG, DB_LONG_ORIG

        !----------------------------------------------------------------------

        !Worksize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        radians = 180.0 / PI

do1:    do i = ILB, IUB+1
            DB_LAT       = YORIG + YY(i)
            DB_LONG_ORIG = XORIG
            IZONE_ORIG   = int((DB_LONG_ORIG + 360.) / 6) + 31 - 60
            IZONE1(1)    = IZONE_ORIG
            DLZONE(i)    = 6378000. * COS(DB_LAT / radians) * 6. / radians

do2:        do j = JLB+1, JUB+1
                DB_LONG   = XORIG + XX(j)
                IZONE1(j) = int((DB_LONG + 360.) / 6) + 31 - 60
            enddo do2
        enddo do1

    !--------------------------------------------------------------------------

    end subroutine ConvertCoordenates

    !--------------------------------------------------------------------------

    subroutine GridPointsLocation1D

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        real, dimension(:  ), pointer       :: XX, YY
        integer                             :: ILB, IUB, JLB, JUB
        integer                             :: i, j
        real, dimension(:  ), pointer       :: XX_Z, YY_Z, XX_U , YY_V , XX_V, YY_U, XX_Cross, YY_Cross


        !----------------------------------------------------------------------

        !Worksize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        XX        => Me%XX
        YY        => Me%YY
        XX_Z      => Me%Compute%XX_Z
        YY_Z      => Me%Compute%YY_Z
        XX_U      => Me%Compute%XX_U
        YY_V      => Me%Compute%YY_V
        XX_V      => Me%Compute%XX_V
        YY_U      => Me%Compute%YY_U
        XX_Cross  => Me%Compute%XX_Cross
        YY_Cross  => Me%Compute%YY_Cross


do1 :   do j = JLB, JUB

            XX_Z(J) = (XX(J) + XX(J+1)) / 2.

            XX_V(J) = XX_Z(J)

        end do do1


do2 :   do i = ILB, IUB

            YY_Z(I) = (YY(I) + YY(I+1)) / 2.

            YY_U(I) = YY_Z(I)

        end do do2

do3 :   do j = JLB, JUB + 1

            XX_U(J)     = XX(J)

            XX_Cross(J) = XX(J)

        end do do3

do4 :   do i = ILB, IUB + 1

            YY_V    (I) = YY(I)

            YY_Cross(I) = YY(I)

        end do do4


        nullify(XX   )
        nullify(YY   )
        nullify(XX_Z )
        nullify(YY_Z )
        nullify(XX_U )
        nullify(YY_V )
        nullify(XX_V )
        nullify(YY_U )

        nullify(XX_Cross)
        nullify(YY_Cross)


    end subroutine GridPointsLocation1D

    !--------------------------------------------------------------------------
    subroutine GridPointsLocation2D

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        real, dimension(:,:), pointer       :: XX, YY
        integer                             :: ILB, IUB, JLB, JUB
        integer                             :: i, j
        real, dimension(:,:), pointer       :: XX_Z, YY_Z, XX_U , YY_V , XX_V, YY_U


        !----------------------------------------------------------------------

        !Worksize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB



        if (Me%CoordType == SIMPLE_GEOG_ .or. Me%CoordType == GEOG_) then

            XX        => Me%LongitudeConn
            YY        => Me%LatitudeConn

        else

            XX        => Me%XX_IE
            YY        => Me%YY_IE

        endif

        XX_Z      => Me%Compute%XX2D_Z
        YY_Z      => Me%Compute%YY2D_Z
        XX_U      => Me%Compute%XX2D_U
        YY_V      => Me%Compute%YY2D_V
        XX_V      => Me%Compute%XX2D_V
        YY_U      => Me%Compute%YY2D_U

do1 :   do j = JLB, JUB
do2 :   do i = ILB, IUB

            XX_Z(I,J) = (XX(I,J) + XX(I,J+1) + XX(I+1,J) + XX(I+1,J+1)) / 4.
            YY_Z(I,J) = (YY(I,J) + YY(I,J+1) + YY(I+1,J) + YY(I+1,J+1)) / 4.

        end do do2
        end do do1


do5:   do j = JLB, JUB + 1
do6:   do i = ILB, IUB

            XX_U(I,J) = (XX(I,J) + XX(I+1,J)) / 2.
            YY_U(I,J) = (YY(I,J) + YY(I+1,J)) / 2.

        end do do6
        end do do5

do7:   do j = JLB, JUB
do8:   do i = ILB, IUB + 1

            XX_V(I,J) = (XX(I,J) + XX(I,J+1)) / 2.
            YY_V(I,J) = (YY(I,J) + YY(I,J+1)) / 2.

        end do do8
        end do do7



        nullify(XX   )
        nullify(YY   )
        nullify(XX_Z )
        nullify(YY_Z )
        nullify(XX_U )
        nullify(YY_V )
        nullify(XX_V )
        nullify(YY_U )



    end subroutine GridPointsLocation2D

    !--------------------------------------------------------------------------

    subroutine ComputeDistances

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB
        integer                             :: i, j
        real(8)                             :: AuxDble

        !----------------------------------------------------------------------

        !Worksize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB


cd1:    if (Me%Grid_Angle /= 0. .or. Me%CoordType == CIRCULAR_ .or. Me%CornersXYInput) then

            !DXX
            do j = JLB, JUB
            do i = ILB, IUB + 1
                AuxDble =                                                &
                    sqrt( (dble(Me%XX_IE(i, j+1))    -    &
                           dble(Me%XX_IE(i, j)))**2. +    &
                          (dble(Me%YY_IE(i, j+1))    -    &
                           dble(Me%YY_IE(i, j)))**2. )

                Me%DXX(i, j) = real (AuxDble)
            enddo
            enddo

            !DYY
            do j = JLB, JUB + 1
            do i = ILB, IUB
                AuxDble =                                                &
                    sqrt( (dble(Me%XX_IE(i+1, j))    -    &
                           dble(Me%XX_IE(i, j)))**2. +    &
                          (dble(Me%YY_IE(i+1, j))    -    &
                           dble(Me%YY_IE(i, j)))**2. )

                Me%DYY(i, j) = real (AuxDble)
            enddo
            enddo

        else  cd1 !Grid_Angle = 0. No grid rotation. Or coordinate type not circular

            !DXX
            do j = JLB, JUB
            do i = ILB, IUB + 1
                Me%DXX(i, j) = Me%XX_IE(i  , j+1) -    &
                                              Me%XX_IE(i  , j  )
            enddo
            enddo

            !DYY
            do j = JLB, JUB + 1
            do i = ILB, IUB
                Me%DYY(i, j) = Me%YY_IE(i+1, j  ) -    &
                                              Me%YY_IE(i  , j  )
            enddo
            enddo


        endif cd1

        !DUX, DVY
        do j = JLB, JUB
        do i = ILB, IUB
            Me%DUX(i, j) = (Me%DXX(i, j) +       &
                                           Me%DXX(i+1, j)) / 2.
            Me%DVY(i, j) = (Me%DYY(i, j) +       &
                                           Me%DYY(i, j+1)) / 2.
        enddo
        enddo

        !DVX
        do j = JLB, JUB - 1
        do i = ILB, IUB + 1
           Me%DVX(i, j) = (Me%DXX(i, j) +        &
                                          Me%DXX(i, j+1)) / 2.
        enddo
        enddo

        !DZX
        do j = JLB, JUB - 1
        do i = ILB, IUB
            Me%DZX(i, j) = (Me%DUX(i, j) +       &
                                           Me%DUX(i, j+1)) / 2.
        enddo
        enddo

        !DUY
        do j = JLB, JUB + 1
        do i = ILB, IUB - 1
            Me%DUY(i, j) = (Me%DYY(i, j) +       &
                                           Me%DYY(i+1, j)) / 2.

        enddo
        enddo

        !DZY
        do j = JLB, JUB
        do i = ILB, IUB - 1
            Me%DZY(i, j) = (Me%DVY(i, j) +       &
                                           Me%DVY(i+1, j)) / 2.
        enddo
        enddo


        if (Me%GridBorderAlongGrid%Type_ == Rectang_ .or. Me%CoordType == SIMPLE_GEOG_) then
            Me%XX_AlongGrid(:, :)  = Me%XX_IE(:,:)
        else

            Me%XX_AlongGrid(ILB:IUB+1, JLB)  = 0.

            do j = JLB, JUB
            do i = ILB, IUB + 1
                Me%XX_AlongGrid(i, j+1) = Me%XX_AlongGrid(i, j) + Me%DXX(i, j)
            enddo
            enddo

        endif

        if (Me%GridBorderAlongGrid%Type_ == Rectang_ .or. Me%CoordType == SIMPLE_GEOG_) then
            Me%YY_AlongGrid(:,:)  = Me%YY_IE(:,:)
        else

            Me%YY_AlongGrid(ILB, JLB:JUB+1)  = 0.

            do j = JLB, JUB + 1
            do i = ILB, IUB
                Me%YY_AlongGrid(i+1, j) = Me%YY_AlongGrid(i, j) + Me%DYY(i, j)
            enddo
            enddo
        endif

    end subroutine ComputeDistances


    !--------------------------------------------------------------------------

    subroutine ComputeRotation

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB
        integer                             :: i, j
        real(8)                             :: x1, x2, y1, y2, dx, dy

        !----------------------------------------------------------------------

        !Worksize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB


cd1:    if (Me%CornersXYInput .or. Me%CoordType == CIRCULAR_) then


            do j = JLB, JUB
            do i = ILB, IUB

                !Rotation X

                x1 = (Me%XX_IE(i, j)   + Me%XX_IE(i+1, j)) / 2.
                x2 = (Me%XX_IE(i, j+1) + Me%XX_IE(i+1, j+1)) / 2.

                y1 = (Me%YY_IE(i, j)   + Me%YY_IE(i+1, j)) / 2.
                y2 = (Me%YY_IE(i, j+1) + Me%YY_IE(i+1, j+1)) / 2.

                dx = x2 - x1
                dy = y2 - y1

                Me%RotationX(i, j) = atan2(dy, dx)

                !Rotation Y

                x1 = (Me%XX_IE(i,   j) + Me%XX_IE(i,   j+1)) / 2.
                x2 = (Me%XX_IE(i+1, j) + Me%XX_IE(i+1, j+1)) / 2.

                y1 = (Me%YY_IE(i,   j) + Me%YY_IE(i  , j+1)) / 2.
                y2 = (Me%YY_IE(i+1, j) + Me%YY_IE(i+1, j+1)) / 2.

                dx = x2 - x1
                dy = y2 - y1

                Me%RotationY(i, j) = atan2(dy, dx)

            enddo
            enddo

         endif cd1

    end subroutine ComputeRotation

    !--------------------------------------------------------------------------

    subroutine ComputeGridCellArea

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB
        integer                             :: i, j

        !Worksize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        do j = JLB, JUB
        do i = ILB, IUB

            Me%GridCellArea(i, j) = Me%DUX(i, j) * Me%DVY(i, j)

        enddo
        enddo


    end subroutine ComputeGridCellArea

    !--------------------------------------------------------------------------

    subroutine ComputeCoriolis

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB
        integer                             :: i, j
        real                                :: OMEGA

        !----------------------------------------------------------------------

        !Worksize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !OMEGA
        !OMEGA = 2. * PI / (24.016 * 3600.)
        !rad/s
        OMEGA = 7.2921e-5

        do j = JLB, JUB
        do i = ILB, IUB
            Me%F(i, j) = 2. * OMEGA * sin(Me%LatitudeZ(i, j)*PI / 180.)
        end do
        end do

    end subroutine ComputeCoriolis

    !--------------------------------------------------------------------------

    subroutine ReCenterHorizontalGrid(HorizontalGridID, Xcenter, Ycenter, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real,              intent(IN)               :: Xcenter, Ycenter
        integer, optional, intent(OUT)              :: STAT


        !Local-----------------------------------------------------------------
        real                                        :: DX, DY, TranslationX, TranslationY
        integer                                     :: STAT_, ready_
        integer                                     :: ILB, IUB, JLB, JUB, ic, jc
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            if (Me%CornersXYInput) then

                ic = int((IUB - ILB + 2)/2.)
                jc = int((JUB - JLB + 2)/2.)

                TranslationX = Xcenter - Me%XX_IE(ic,jc)
                TranslationY = Ycenter - Me%YY_IE(ic,jc)

                Me%XX_IE(:,:) = Me%XX_IE(:,:) + TranslationX
                Me%YY_IE(:,:) = Me%YY_IE(:,:) + TranslationY

            else
                DX = Me%XX(JUB+1) - Me%XX(JLB)
                DY = Me%YY(IUB+1) - Me%YY(ILB)

                TranslationX = Xcenter - Me%Xorig - 0.5*DX
                TranslationY = Ycenter - Me%Yorig - 0.5*DY

                Me%Xorig = Me%Xorig + TranslationX
                Me%Yorig = Me%Yorig + TranslationY

            endif

            call GenerateGrid

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ReCenterHorizontalGrid

    !--------------------------------------------------------------------------

#ifdef _USE_MPI

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! When using a domain decomposition approach this subroutine aggregates in one global  !
    ! domain 1D vectors associated with a specific sub-domain                              !
    !                                                                                      !
    ! Input : specific domain 1D vector                                                    !
    ! OutPut: global domain 1D vector                                                      !
    ! Author: Paulo Chambel (2015/11)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine JoinGridData_1D(HorizontalGridID, In1D, Out1D, Type1DIJ, dj, di, STAT)


        !Arguments------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real,   dimension(:), pointer               :: In1D
        real,   dimension(:), pointer               :: Out1D
        integer                                     :: Type1DIJ
        integer             , optional              :: di, dj
        integer             , optional              :: STAT

        !Local---------------------------------------------------------------
        integer                                     :: di_, dj_
        integer                                     :: STAT_, ready_

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(di     )) then
                di_ = di
            else
                di_ = 0
            endif

            if (present(dj     )) then
                dj_ = dj
            else
                dj_ = 0
            endif

            call JoinGridData_In(In1D, Out1D, Type1DIJ, dj_, di_)

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

    end subroutine JoinGridData_1D

    !-----------------------------------------------------------------------------------

    subroutine JoinGridData_1D_In(In1D, Out1D, Type1DIJ, dj_, di_)


        !Arguments------------------------------------------------------------
        real,   dimension(:), pointer               :: In1D
        real,   dimension(:), pointer               :: Out1D
        integer                                     :: Type1DIJ
        integer                                     :: di_, dj_

        !Local---------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB, imin, imax, i

        integer                                     :: Source, Destination
        integer                                     :: iSize
        integer, save                               :: Precision
        integer                                     :: status(MPI_STATUS_SIZE)

        type(T_Size2D)                              :: Inner, Mapping

        !Begin---------------------------------------------------------------


        IUB     = Me%DDecomp%Global%IUB
        ILB     = Me%DDecomp%Global%ILB
        JUB     = Me%DDecomp%Global%JUB
        JLB     = Me%DDecomp%Global%JLB

        if (Type1DIJ /= Type1DI .and. Type1DIJ /= Type1DJ) then
            stop 'JoinGridData_1D - ModuleHorizontalGrid - ERR10'
        endif

        if (Me%DDecomp%Master) then

            Inner   = Me%DDecomp%Inner
            Mapping = Me%DDecomp%Mapping

            !Copy form Master (input) to Global (output)
            if (Type1DIJ == Type1DI) then
                if (.not. associated(Out1D)) allocate(Out1D(ILB:IUB+di_))
                Out1D(Mapping%ILB:Mapping%IUB+di_) = In1D(Inner%ILB:Inner%IUB+di_)
            endif
            if (Type1DIJ == Type1DJ) then
                if (.not. associated(Out1D)) allocate(Out1D(JLB:JUB+dj_))
                Out1D(Mapping%JLB:Mapping%JUB+dj_) = In1D(Inner%JLB:Inner%JUB+dj_)
            endif

            !Receive from slaves (input) to Global (output)
            do i=1, Me%DDecomp%Nslaves

                Mapping = Me%DDecomp%Slaves_Mapping(i)

                if     (Type1DIJ == Type1DI) then
                    imin      = Mapping%ILB
                    imax      = Mapping%IUB + di_
                elseif (Type1DIJ == Type1DJ) then
                    imin      = Mapping%JLB
                    imax      = Mapping%JUB + dj_
                endif

                iSize     = (imax-imin+1)

                Precision = MPIKind(In1D)

                Source    = Me%DDecomp%Slaves_MPI_ID(i)

                call MPI_Recv (Out1D(imin:imax), iSize, Precision, Source, 9001, MPI_COMM_WORLD, status, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_1D_In - ModuleHorizontalGrid - ERR40'

            enddo

        else

            Inner       = Me%DDecomp%Inner

            if     (Type1DIJ == Type1DI) then
                imin      = Inner%ILB
                imax      = Inner%IUB + di_
            elseif (Type1DIJ == Type1DJ) then
                imin      = Inner%JLB
                imax      = Inner%JUB + dj_
            endif

            iSize     = (imax-imin+1)

            Precision   = MPIKind(In1D)

            Destination = Me%DDecomp%Master_MPI_ID

            call MPI_Send (In1D(imin:imax), iSize, Precision, Destination, 9001, MPI_COMM_WORLD, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_1D_In - ModuleHorizontalGrid - ERR50'

        endif


    end subroutine JoinGridData_1D_In

    !-----------------------------------------------------------------------------------


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! When using a domain decomposition approach this subroutine aggregates in one global  !
    ! domain 2D single precision matrixes associated with a specific sub-domain            !
    !                                                                                      !
    ! Input : specific domain 2D matrixes                                                  !
    ! OutPut: global domain 2D matrixes                                                    !
    ! Author: Paulo Chambel (2015/11)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine JoinGridData_2D_R4(HorizontalGridID, In2D, Out2D, dj, di, STAT)


        !Arguments------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real(4),dimension(:,:), pointer             :: In2D
        real(4),dimension(:,:), pointer             :: Out2D
        integer               , optional            :: di, dj
        integer               , optional            :: STAT

        !Local---------------------------------------------------------------
        integer                                     :: di_, dj_
        integer                                     :: STAT_, ready_

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(di     )) then
                di_ = di
            else
                di_ = 0
            endif

            if (present(dj     )) then
                dj_ = dj
            else
                dj_ = 0
            endif

            call JoinGridData_In(In2D, Out2D, dj_, di_)


            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

    end subroutine JoinGridData_2D_R4

    !-----------------------------------------------------------------------------------

    subroutine JoinGridData_2D_R4_In(In2D, Out2D, dj_, di_)


        !Arguments------------------------------------------------------------
        real(4),dimension(:,:), pointer             :: In2D
        real(4),dimension(:,:), pointer             :: Out2D
        integer                                     :: di_, dj_
        !Local---------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: imin, imax, jmin, jmax

        integer                                     :: Source, Destination
        integer                                     :: iSize, i
        integer, save                               :: Precision
        integer                                     :: status(MPI_STATUS_SIZE)

        type(T_Size2D)                              :: Inner, Mapping

        !Begin---------------------------------------------------------------

        IUB     = Me%DDecomp%Global%IUB
        ILB     = Me%DDecomp%Global%ILB
        JUB     = Me%DDecomp%Global%JUB
        JLB     = Me%DDecomp%Global%JLB

        if (Me%DDecomp%Master) then

            !Copy form Master (input) to Global (output)
            Inner   = Me%DDecomp%Inner
            Mapping = Me%DDecomp%Mapping

            if (.not. associated(Out2D)) allocate(Out2D(ILB:IUB+di_,JLB:JUB+dj_))

            Out2D(Mapping%ILB:Mapping%IUB+di_,Mapping%JLB:Mapping%JUB+dj_) =        &
            In2D (  Inner%ILB:  Inner%IUB+di_,  Inner%JLB:  Inner%JUB+dj_)


            !Receive from slaves (input) to Global (output)
            do i=1, Me%DDecomp%Nslaves

                Mapping   = Me%DDecomp%Slaves_Mapping(i)

                imax      = Mapping%IUB+di_
                imin      = Mapping%ILB
                jmax      = Mapping%JUB+dj_
                jmin      = Mapping%JLB

                iSize     = (imax-imin+1)*(jmax-jmin+1)

                Precision = MPI_REAL

                Source    =  Me%DDecomp%Slaves_MPI_ID(i)

                !Receive from sub-domanins to global output
                call MPI_Recv (Out2D(imin:imax,jmin:jmax),                          &
                               iSize, Precision, Source, 9002, MPI_COMM_WORLD, status, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_2D_R4_In - ModuleHorizontalGrid - ERR40'

            enddo

        else

            Inner       = Me%DDecomp%Inner

            imax        = Inner%IUB+di_
            imin        = Inner%ILB
            jmax        = Inner%JUB+dj_
            jmin        = Inner%JLB

            iSize       = (imax-imin+1)*(jmax-jmin+1)

            Precision   = MPI_REAL

            Destination = Me%DDecomp%Master_MPI_ID

            call MPI_Send (In2D(imin:imax,jmin:jmax), iSize, Precision,             &
                           Destination, 9002, MPI_COMM_WORLD, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_2D - ModuleHorizontalGrid - ERR50'

        endif


    end subroutine JoinGridData_2D_R4_In

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! When using a domain decomposition approach this subroutine aggregates in one global  !
    ! domain 2D single precision matrixes associated with a specific sub-domain            !
    !                                                                                      !
    ! Input : specific domain 2D matrixes                                                  !
    ! OutPut: global domain 2D matrixes                                                    !
    ! Author: Paulo Chambel (2015/11)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine JoinGridData_2D_R8(HorizontalGridID, In2D, Out2D, dj, di, STAT)


        !Arguments------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real(8),dimension(:,:), pointer             :: In2D
        real(8),dimension(:,:), pointer             :: Out2D
        integer               , optional            :: di, dj
        integer               , optional            :: STAT

        !Local---------------------------------------------------------------
        integer                                     :: di_, dj_
        integer                                     :: STAT_, ready_

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(di     )) then
                di_ = di
            else
                di_ = 0
            endif

            if (present(dj     )) then
                dj_ = dj
            else
                dj_ = 0
            endif

            call JoinGridData_In(In2D, Out2D, dj_, di_)


            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

    end subroutine JoinGridData_2D_R8

    !-----------------------------------------------------------------------------------

    subroutine JoinGridData_2D_R8_In(In2D, Out2D, dj_, di_)


        !Arguments------------------------------------------------------------
        real(8),dimension(:,:), pointer             :: In2D
        real(8),dimension(:,:), pointer             :: Out2D
        integer                                     :: di_, dj_
        !Local---------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: imin, imax, jmin, jmax

        integer                                     :: Source, Destination
        integer                                     :: iSize, i
        integer, save                               :: Precision
        integer                                     :: status(MPI_STATUS_SIZE)

        type(T_Size2D)                              :: Inner, Mapping

        !Begin---------------------------------------------------------------

        IUB     = Me%DDecomp%Global%IUB
        ILB     = Me%DDecomp%Global%ILB
        JUB     = Me%DDecomp%Global%JUB
        JLB     = Me%DDecomp%Global%JLB

        if (Me%DDecomp%Master) then

            !Copy form Master (input) to Global (output)
            Inner   = Me%DDecomp%Inner
            Mapping = Me%DDecomp%Mapping

            if (.not. associated(Out2D)) allocate(Out2D(ILB:IUB+di_,JLB:JUB+dj_))

            Out2D(Mapping%ILB:Mapping%IUB+di_,Mapping%JLB:Mapping%JUB+dj_) =        &
            In2D (  Inner%ILB:  Inner%IUB+di_,  Inner%JLB:  Inner%JUB+dj_)


            !Receive from slaves (input) to Global (output)
            do i=1, Me%DDecomp%Nslaves

                Mapping   = Me%DDecomp%Slaves_Mapping(i)

                imax      = Mapping%IUB+di_
                imin      = Mapping%ILB
                jmax      = Mapping%JUB+dj_
                jmin      = Mapping%JLB

                iSize     = (imax-imin+1)*(jmax-jmin+1)

                Precision = MPI_DOUBLE_PRECISION

                Source    =  Me%DDecomp%Slaves_MPI_ID(i)

                !Receive from sub-domanins to global output
                call MPI_Recv (Out2D(imin:imax,jmin:jmax),                          &
                               iSize, Precision, Source, 9003, MPI_COMM_WORLD, status, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_2D_R8_In - ModuleHorizontalGrid - ERR40'

            enddo

        else

            Inner       = Me%DDecomp%Inner

            imax        = Inner%IUB+di_
            imin        = Inner%ILB
            jmax        = Inner%JUB+dj_
            jmin        = Inner%JLB

            iSize       = (imax-imin+1)*(jmax-jmin+1)

            Precision   = MPI_DOUBLE_PRECISION

            Destination = Me%DDecomp%Master_MPI_ID

            call MPI_Send (In2D(imin:imax,jmin:jmax), iSize, Precision,             &
                           Destination, 9003, MPI_COMM_WORLD, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_2D_R8_In - ModuleHorizontalGrid - ERR50'

        endif


    end subroutine JoinGridData_2D_R8_In


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! When using a domain decomposition approach this subroutine aggregates in one global  !
    ! domain 2D matrixes associated with a specific sub-domain                             !
    !                                                                                      !
    ! Input : specific domain 2D matrixes (integer)                                        !
    ! OutPut: global domain 2D matrixes (integer)                                          !
    ! Author: Paulo Chambel (2015/11)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine JoinGridData_2Dint(HorizontalGridID, In2D, Out2D, dj, di, STAT)


        !Arguments------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer,dimension(:,:), pointer             :: In2D
        integer,dimension(:,:), pointer             :: Out2D
        integer               , optional            :: di, dj
        integer               , optional            :: STAT

        !Local---------------------------------------------------------------
        integer                                 :: di_, dj_
        integer                                 :: STAT_, ready_

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(di     )) then
                di_ = di
            else
                di_ = 0
            endif

            if (present(dj     )) then
                dj_ = dj
            else
                dj_ = 0
            endif

            call JoinGridData_In(In2D, Out2D, dj_, di_)

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

    end subroutine JoinGridData_2Dint

    !-----------------------------------------------------------------------------------

    subroutine JoinGridData_2Dint_In(In2D, Out2D, dj_, di_)


        !Arguments------------------------------------------------------------
        integer,dimension(:,:), pointer         :: In2D
        integer,dimension(:,:), pointer         :: Out2D
        integer                                 :: di_, dj_

        !Local---------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: imin, imax, jmin, jmax, i

        integer                                 :: Source, Destination
        integer                                 :: iSize
        integer, save                           :: Precision
        integer                                 :: status(MPI_STATUS_SIZE)

        type(T_Size2D)                          :: Inner, Mapping

        !Begin---------------------------------------------------------------

        IUB     = Me%DDecomp%Global%IUB
        ILB     = Me%DDecomp%Global%ILB
        JUB     = Me%DDecomp%Global%JUB
        JLB     = Me%DDecomp%Global%JLB

        if (Me%DDecomp%Master) then

            Inner   = Me%DDecomp%Inner
            Mapping = Me%DDecomp%Mapping

            !Copy form Master (input) to Global (output)
            if (.not. associated(Out2D)) allocate(Out2D(ILB:IUB+di_,JLB:JUB+dj_))

            Out2D(Mapping%ILB:Mapping%IUB+di_,Mapping%JLB:Mapping%JUB+dj_) =        &
            In2D (  Inner%ILB:  Inner%IUB+di_,  Inner%JLB:  Inner%JUB+dj_)

            !Receive from slaves (input) to Global (output)
            do i=1, Me%DDecomp%Nslaves

                Mapping   = Me%DDecomp%Slaves_Mapping(i)

                imax      = Mapping%IUB+di_
                imin      = Mapping%ILB
                jmax      = Mapping%JUB+dj_
                jmin      = Mapping%JLB

                iSize     = (imax-imin+1)*(jmax-jmin+1)

                Precision = MPI_INTEGER

                Source    =  Me%DDecomp%Slaves_MPI_ID(i)

                !Receive from sub-domanins to global output
                call MPI_Recv (Out2D(imin:imax,jmin:jmax),                          &
                               iSize, Precision, Source, 9004, MPI_COMM_WORLD, status, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_2Dint_In - ModuleHorizontalGrid - ERR40'

            enddo

        else

            Inner       = Me%DDecomp%Inner

            imax        = Inner%IUB+di_
            imin        = Inner%ILB
            jmax        = Inner%JUB+dj_
            jmin        = Inner%JLB

            iSize       = (imax-imin+1)*(jmax-jmin+1)

            Precision   = MPI_INTEGER

            Destination = Me%DDecomp%Master_MPI_ID

            call MPI_Send (In2D(imin:imax,jmin:jmax), iSize, Precision,             &
                           Destination, 9004, MPI_COMM_WORLD, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_2Dint_In - ModuleHorizontalGrid - ERR50'

        endif


    end subroutine JoinGridData_2Dint_In


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! When using a domain decomposition approach this subroutine aggregates in one global  !
    ! domain 3D matrixes associated with a specific sub-domain                             !
    !                                                                                      !
    ! Input : specific domain 3D matrixes                                                  !
    ! OutPut: global domain 3D matrixes                                                    !
    ! Author: Paulo Chambel (2015/11)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine JoinGridData_3D_R4(HorizontalGridID, In3D, Out3D, KLB, KUB, dj, di, STAT)


        !Arguments------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real(4),dimension(:,:,:), pointer           :: In3D
        real(4),dimension(:,:,:), pointer           :: Out3D
        integer                                     :: KLB, KUB
        integer                 , optional          :: di, dj
        integer                 , optional          :: STAT


        !Local---------------------------------------------------------------
        integer                                     :: di_, dj_
        integer                                     :: STAT_, ready_

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(di     )) then
                di_ = di
            else
                di_ = 0
            endif

            if (present(dj     )) then
                dj_ = dj
            else
                dj_ = 0
            endif

            call JoinGridData_In(In3D, Out3D, KLB, KUB, dj_, di_)

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_


    end subroutine JoinGridData_3D_R4

    !-----------------------------------------------------------------------------------


    subroutine JoinGridData_3D_R4_In(In3D, Out3D, KLB, KUB, dj_, di_)


        !Arguments------------------------------------------------------------
        real(4),dimension(:,:,:), pointer           :: In3D
        real(4),dimension(:,:,:), pointer           :: Out3D
        integer                                     :: KLB, KUB
        integer                                     :: di_, dj_


        !Local---------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: imin, imax, jmin, jmax, i
        integer                                     :: Source, Destination
        integer                                     :: iSize
        integer, save                               :: Precision
        integer                                     :: status(MPI_STATUS_SIZE)

        type(T_Size2D)                              :: Inner, Mapping

        !Begin---------------------------------------------------------------

        IUB     = Me%DDecomp%Global%IUB
        ILB     = Me%DDecomp%Global%ILB
        JUB     = Me%DDecomp%Global%JUB
        JLB     = Me%DDecomp%Global%JLB

        if (Me%DDecomp%Master) then

            Inner   = Me%DDecomp%Inner
            Mapping = Me%DDecomp%Mapping

            !Copy form Master (input) to Global (output)
            if (.not. associated(Out3D)) allocate(Out3D(ILB:IUB+di_,JLB:JUB+dj_,KLB:KUB))

            Out3D(Mapping%ILB:Mapping%IUB+di_,Mapping%JLB:Mapping%JUB+dj_,KLB:KUB) =        &
            In3D (  Inner%ILB:  Inner%IUB+di_,  Inner%JLB:  Inner%JUB+dj_,KLB:KUB)

            !Receive from slaves (input) to Global (output)
            do i=1, Me%DDecomp%Nslaves

                Mapping   = Me%DDecomp%Slaves_Mapping(i)

                imax      = Mapping%IUB+di_
                imin      = Mapping%ILB
                jmax      = Mapping%JUB+dj_
                jmin      = Mapping%JLB

                iSize     = (imax-imin+1)*(jmax-jmin+1)*(KUB-KLB+1)

                Precision = MPI_REAL

                Source    =  Me%DDecomp%Slaves_MPI_ID(i)

                !Receive from sub-domanins to global output
                call MPI_Recv (Out3D(imin:imax,jmin:jmax,KLB:KUB),                  &
                               iSize, Precision, Source, 9005, MPI_COMM_WORLD, status, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_3D_R4_In - ModuleHorizontalGrid - ERR40'

            enddo

        else

            Inner       = Me%DDecomp%Inner

            imax        = Inner%IUB+di_
            imin        = Inner%ILB
            jmax        = Inner%JUB+dj_
            jmin        = Inner%JLB

            iSize       = (imax-imin+1)*(jmax-jmin+1)*(KUB-KLB+1)

            Precision   = MPI_REAL

            Destination = Me%DDecomp%Master_MPI_ID

            call MPI_Send (In3D(imin:imax,jmin:jmax,KLB:KUB), iSize, Precision,             &
                           Destination, 9005, MPI_COMM_WORLD, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_3D_R4_In - ModuleHorizontalGrid - ERR50'

        endif


    end subroutine JoinGridData_3D_R4_In

    !-----------------------------------------------------------------------------------

    subroutine JoinGridData_3D_R8(HorizontalGridID, In3D, Out3D, KLB, KUB, dj, di, STAT)


        !Arguments------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real(8),dimension(:,:,:), pointer           :: In3D
        real(8),dimension(:,:,:), pointer           :: Out3D
        integer                                     :: KLB, KUB
        integer                 , optional          :: di, dj
        integer                 , optional          :: STAT


        !Local---------------------------------------------------------------
        integer                                     :: di_, dj_
        integer                                     :: STAT_, ready_

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(di     )) then
                di_ = di
            else
                di_ = 0
            endif

            if (present(dj     )) then
                dj_ = dj
            else
                dj_ = 0
            endif

            call JoinGridData_In(In3D, Out3D, KLB, KUB, dj_, di_)

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_


    end subroutine JoinGridData_3D_R8

    !-----------------------------------------------------------------------------------


    subroutine JoinGridData_3D_R8_In(In3D, Out3D, KLB, KUB, dj_, di_)


        !Arguments------------------------------------------------------------
        real(8),dimension(:,:,:), pointer           :: In3D
        real(8),dimension(:,:,:), pointer           :: Out3D
        integer                                     :: KLB, KUB
        integer                                     :: di_, dj_


        !Local---------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: imin, imax, jmin, jmax, i
        integer                                     :: Source, Destination
        integer                                     :: iSize
        integer, save                               :: Precision
        integer                                     :: status(MPI_STATUS_SIZE)

        type(T_Size2D)                              :: Inner, Mapping

        !Begin---------------------------------------------------------------


        IUB     = Me%DDecomp%Global%IUB
        ILB     = Me%DDecomp%Global%ILB
        JUB     = Me%DDecomp%Global%JUB
        JLB     = Me%DDecomp%Global%JLB

        if (Me%DDecomp%Master) then

            Inner   = Me%DDecomp%Inner
            Mapping = Me%DDecomp%Mapping

            !Copy form Master (input) to Global (output)
            if (.not. associated(Out3D)) allocate(Out3D(ILB:IUB+di_,JLB:JUB+dj_,KLB:KUB))

            Out3D(Mapping%ILB:Mapping%IUB+di_,Mapping%JLB:Mapping%JUB+dj_,KLB:KUB) =        &
            In3D (  Inner%ILB:  Inner%IUB+di_,  Inner%JLB:  Inner%JUB+dj_,KLB:KUB)

            !Receive from slaves (input) to Global (output)
            do i=1, Me%DDecomp%Nslaves

                Mapping   = Me%DDecomp%Slaves_Mapping(i)

                imax      = Mapping%IUB+di_
                imin      = Mapping%ILB
                jmax      = Mapping%JUB+dj_
                jmin      = Mapping%JLB

                iSize     = (imax-imin+1)*(jmax-jmin+1)*(KUB-KLB+1)

                Precision = MPI_DOUBLE_PRECISION

                Source    =  Me%DDecomp%Slaves_MPI_ID(i)

                !Receive from sub-domanins to global output
                call MPI_Recv (Out3D(imin:imax,jmin:jmax,KLB:KUB),                  &
                               iSize, Precision, Source, 9006, MPI_COMM_WORLD, status, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_3D_R8_In - ModuleHorizontalGrid - ERR40'

            enddo

        else

            Inner       = Me%DDecomp%Inner

            imax        = Inner%IUB+di_
            imin        = Inner%ILB
            jmax        = Inner%JUB+dj_
            jmin        = Inner%JLB

            iSize       = (imax-imin+1)*(jmax-jmin+1)*(KUB-KLB+1)

            Precision   = MPI_DOUBLE_PRECISION

            Destination = Me%DDecomp%Master_MPI_ID

            call MPI_Send (In3D(imin:imax,jmin:jmax,KLB:KUB), iSize, Precision,             &
                           Destination, 9006, MPI_COMM_WORLD, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_3D_R8_In - ModuleHorizontalGrid - ERR50'

        endif


    end subroutine JoinGridData_3D_R8_In

    !-----------------------------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! When using a domain decomposition approach this subroutine aggregates in one global  !
    ! domain 3D integer matrixes associated with a specific sub-domain                     !
    !                                                                                      !
    ! Input : specific domain 3D matrixes (integer)                                        !
    ! OutPut: global domain 3D matrixes   (integer)                                        !
    ! Author: Paulo Chambel (2015/11)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine JoinGridData_3Dint(HorizontalGridID, In3D, Out3D, KLB, KUB, dj, di, STAT)


        !Arguments------------------------------------------------------------
        integer                                         :: HorizontalGridID
        integer, dimension(:,:,:), pointer              :: In3D
        integer, dimension(:,:,:), pointer              :: Out3D
        integer                                         :: KLB, KUB
        integer                  , optional             :: di, dj
        integer                  , optional             :: STAT

        !Local---------------------------------------------------------------

        integer                                         :: di_, dj_
        integer                                         :: STAT_, ready_

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(di     )) then
                di_ = di
            else
                di_ = 0
            endif

            if (present(dj     )) then
                dj_ = dj
            else
                dj_ = 0
            endif

            call JoinGridData_In(In3D, Out3D, KLB, KUB, dj_, di_)

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

    end subroutine JoinGridData_3Dint

    !------------------------------------------------------------------------------

    subroutine JoinGridData_3Dint_In(In3D, Out3D, KLB, KUB, dj_, di_)


        !Arguments------------------------------------------------------------
        integer, dimension(:,:,:), pointer              :: In3D
        integer, dimension(:,:,:), pointer              :: Out3D
        integer                                         :: KLB, KUB
        integer                                         :: di_, dj_

        !Local---------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: ILB, IUB, JLB, JUB

        integer                                         :: imin, imax, jmin, jmax, i

        integer                                         :: Source, Destination
        integer                                         :: iSize
        integer, save                                   :: Precision
        integer                                         :: status(MPI_STATUS_SIZE)

        type(T_Size2D)                                  :: Inner, Mapping

        !Begin---------------------------------------------------------------


        IUB     = Me%DDecomp%Global%IUB
        ILB     = Me%DDecomp%Global%ILB
        JUB     = Me%DDecomp%Global%JUB
        JLB     = Me%DDecomp%Global%JLB

        if (Me%DDecomp%Master) then

            Inner   = Me%DDecomp%Inner
            Mapping = Me%DDecomp%Mapping

            !Copy form Master (input) to Global (output)
            if (.not. associated(Out3D)) allocate(Out3D(ILB:IUB+di_,JLB:JUB+dj_,KLB:KUB))

            Out3D(Mapping%ILB:Mapping%IUB+di_,Mapping%JLB:Mapping%JUB+dj_,KLB:KUB) =        &
            In3D (  Inner%ILB:  Inner%IUB+di_,  Inner%JLB:  Inner%JUB+dj_,KLB:KUB)

            !Receive from slaves (input) to Global (output)
            do i=1, Me%DDecomp%Nslaves

                Mapping   = Me%DDecomp%Slaves_Mapping(i)

                imax      = Mapping%IUB+di_
                imin      = Mapping%ILB
                jmax      = Mapping%JUB+dj_
                jmin      = Mapping%JLB

                iSize     = (imax-imin+1)*(jmax-jmin+1)*(KUB-KLB+1)

                Precision = MPI_INTEGER

                Source    =  Me%DDecomp%Slaves_MPI_ID(i)

                !Receive from sub-domanins to global output
                call MPI_Recv (Out3D(imin:imax,jmin:jmax,KLB:KUB),                  &
                               iSize, Precision, Source, 9007, MPI_COMM_WORLD, status, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_3Dint_In - ModuleHorizontalGrid - ERR40'

            enddo

        else

            Inner       = Me%DDecomp%Inner

            imax        = Inner%IUB+di_
            imin        = Inner%ILB
            jmax        = Inner%JUB+dj_
            jmin        = Inner%JLB

            iSize       = (imax-imin+1)*(jmax-jmin+1)*(KUB-KLB+1)

            Precision   = MPI_INTEGER

            Destination = Me%DDecomp%Master_MPI_ID

            call MPI_Send (In3D(imin:imax,jmin:jmax,KLB:KUB), iSize, Precision,             &
                           Destination, 9007, MPI_COMM_WORLD, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'JoinGridData_3Dint_In - ModuleHorizontalGrid - ERR50'

        endif

    end subroutine JoinGridData_3Dint_In

    !------------------------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! This subroutine broadcast for each domain a 2D matrix                                !
    !                                                                                      !
    ! Input : global domain 2D matrixes                                                    !
    ! OutPut: specific domain 2D matrixes                                                  !
    ! Author: Paulo Chambel (2015/12)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine BroadcastGridData2D_R4(HorizontalGridID, In2D, Out2D, STAT)


        !Arguments------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real(4),   dimension(:,:), pointer          :: In2D
        real(4),   dimension(:,:), pointer          :: Out2D
        integer               , optional            :: STAT

        !Local---------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call BroadcastGridData_In(In2D, Out2D)

            STAT_ = SUCCESS_

        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

    end subroutine BroadcastGridData2D_R4

    !------------------------------------------------------------------------------

    subroutine BroadcastGridData2D_R4_In(In2D, Out2D)


        !Arguments------------------------------------------------------------
        real(4),   dimension(:,:), pointer          :: In2D
        real(4),   dimension(:,:), pointer          :: Out2D

        !Local---------------------------------------------------------------
        real(8), dimension(:,:), pointer            :: Aux2D
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB

        integer                                     :: Source, Destination
        integer                                     :: iSize, i
        integer, save                               :: Precision
        integer                                     :: status(MPI_STATUS_SIZE)

        type(T_Size2D)                              :: WorkSize, HaloMap

        !Begin---------------------------------------------------------------

        IUB     = Me%DDecomp%Global%IUB
        ILB     = Me%DDecomp%Global%ILB
        JUB     = Me%DDecomp%Global%JUB
        JLB     = Me%DDecomp%Global%JLB

        if (Me%DDecomp%Master) then

            WorkSize = Me%WorkSize
            HaloMap  = Me%DDecomp%HaloMap

            Out2D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB) =           &
                In2D(HaloMap%ILB:HaloMap%IUB,HaloMap%JLB:HaloMap%JUB)

            do i=1, Me%DDecomp%Nslaves

                WorkSize =  Me%DDecomp%Slaves_Size   (i)
                HaloMap  =  Me%DDecomp%Slaves_HaloMap(i)

                allocate(Aux2D (WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB))

                Aux2D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB) =       &
                    In2D(HaloMap%ILB:HaloMap%IUB, HaloMap%JLB:HaloMap%JUB)

                iSize       = (WorkSize%IUB-WorkSize%ILB+1) * (WorkSize%JUB-WorkSize%JLB+1)

                Precision   = MPI_DOUBLE_PRECISION

                Destination =  Me%DDecomp%Slaves_MPI_ID(i)

                call MPI_Send (Aux2D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB), iSize, Precision,   &
                               Destination, 90016, MPI_COMM_WORLD, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'BroadcastGridData2D_R4_In - ModuleHorizontalGrid - ERR10'

                deallocate(Aux2D)
                nullify   (Aux2D)

            enddo

        else
            WorkSize    = Me%WorkSize

            allocate(Aux2D (WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB))


            iSize       = (WorkSize%IUB-WorkSize%ILB+1) * (WorkSize%JUB-WorkSize%JLB+1)
            Precision   = MPI_DOUBLE_PRECISION
            Source      = Me%DDecomp%Master_MPI_ID

            call MPI_Recv (Aux2D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB), iSize, Precision,   &
                           Source, 90016, MPI_COMM_WORLD, status, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'BroadcastGridData2D_R4_In - ModuleHorizontalGrid - ERR20'

            Out2D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB) = &
                Aux2D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB)

            deallocate(Aux2D)
            nullify   (Aux2D)

        endif

    end subroutine BroadcastGridData2D_R4_In

    !------------------------------------------------------------------------------

    subroutine BroadcastGridData2D_R8(HorizontalGridID, In2D, Out2D, STAT)


        !Arguments------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real(8),   dimension(:,:), pointer          :: In2D
        real(8),   dimension(:,:), pointer          :: Out2D
        integer               , optional            :: STAT

        !Local---------------------------------------------------------------

        integer                                     :: STAT_, ready_

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call BroadcastGridData_In(In2D, Out2D)

            STAT_ = SUCCESS_

        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

    end subroutine BroadcastGridData2D_R8

    !------------------------------------------------------------------------------

    subroutine BroadcastGridData2D_R8_In(In2D, Out2D)


        !Arguments------------------------------------------------------------
        real(8),   dimension(:,:), pointer          :: In2D
        real(8),   dimension(:,:), pointer          :: Out2D

        !Local---------------------------------------------------------------
        real(8), dimension(:,:), pointer            :: Aux2D
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB

        integer                                     :: Source, Destination
        integer                                     :: iSize, i
        integer, save                               :: Precision
        integer                                     :: status(MPI_STATUS_SIZE)

        type(T_Size2D)                              :: WorkSize, HaloMap

        !Begin---------------------------------------------------------------

        IUB     = Me%DDecomp%Global%IUB
        ILB     = Me%DDecomp%Global%ILB
        JUB     = Me%DDecomp%Global%JUB
        JLB     = Me%DDecomp%Global%JLB

        if (Me%DDecomp%Master) then

            WorkSize = Me%WorkSize
            HaloMap  = Me%DDecomp%HaloMap

            Out2D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB) =           &
                In2D(HaloMap%ILB:HaloMap%IUB,HaloMap%JLB:HaloMap%JUB)

            do i=1, Me%DDecomp%Nslaves

                WorkSize =  Me%DDecomp%Slaves_Size   (i)
                HaloMap  =  Me%DDecomp%Slaves_HaloMap(i)

                allocate(Aux2D (WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB))

                Aux2D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB) =       &
                    In2D(HaloMap%ILB:HaloMap%IUB, HaloMap%JLB:HaloMap%JUB)

                iSize       = (WorkSize%IUB-WorkSize%ILB+1) * (WorkSize%JUB-WorkSize%JLB+1)

                Precision   = MPI_DOUBLE_PRECISION

                Destination =  Me%DDecomp%Slaves_MPI_ID(i)

                call MPI_Send (Aux2D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB), iSize, Precision,   &
                               Destination, 90026, MPI_COMM_WORLD, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'BroadcastGridData2D_R8_In - ModuleHorizontalGrid - ERR10'

                deallocate(Aux2D)
                nullify   (Aux2D)

            enddo

        else
            WorkSize    = Me%WorkSize

            allocate(Aux2D (WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB))


            iSize       = (WorkSize%IUB-WorkSize%ILB+1) * (WorkSize%JUB-WorkSize%JLB+1)
            Precision   = MPI_DOUBLE_PRECISION
            Source      = Me%DDecomp%Master_MPI_ID

            call MPI_Recv (Aux2D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB), iSize, Precision,   &
                           Source, 90026, MPI_COMM_WORLD, status, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'BroadcastGridData2D_R8_In - ModuleHorizontalGrid - ERR20'

            Out2D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB) = &
                Aux2D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB)

            deallocate(Aux2D)
            nullify   (Aux2D)

        endif

    end subroutine BroadcastGridData2D_R8_In
    !------------------------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! This subroutine broadcast for each domain a 3D matrix                                !
    !                                                                                      !
    ! Input : global domain 3D matrixes                                                    !
    ! OutPut: specific domain 3D matrixes                                                  !
    ! Author: Paulo Chambel (2015/12)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine BroadcastGridData3D_R4(HorizontalGridID, In3D, Out3D, KLB, KUB, STAT)


        !Arguments------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real(4),   dimension(:,:,:), pointer        :: In3D
        real(4),   dimension(:,:,:), pointer        :: Out3D
        integer                                     :: KLB, KUB
        integer               , optional            :: STAT

        !Local---------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call BroadcastGridData_In(In3D, Out3D, KLB, KUB)

            STAT_ = SUCCESS_

        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

    end subroutine BroadcastGridData3D_R4

    !-----------------------------------------------------------------------------------


    subroutine BroadcastGridData3D_R4_In(In3D, Out3D, KLB, KUB)


        !Arguments------------------------------------------------------------
        real(4),   dimension(:,:,:), pointer        :: In3D
        real(4),   dimension(:,:,:), pointer        :: Out3D
        integer                                     :: KLB, KUB
        !Local---------------------------------------------------------------
        real(8), dimension(:,:,:), pointer          :: Aux3D
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB

        integer                                     :: Source, Destination
        integer                                     :: iSize, i
        integer, save                               :: Precision
        integer                                     :: status(MPI_STATUS_SIZE)

        type(T_Size2D)                              :: WorkSize, HaloMap

        !Begin---------------------------------------------------------------

        IUB     = Me%DDecomp%Global%IUB
        ILB     = Me%DDecomp%Global%ILB
        JUB     = Me%DDecomp%Global%JUB
        JLB     = Me%DDecomp%Global%JLB

        if (Me%DDecomp%Master) then

            WorkSize = Me%WorkSize
            HaloMap  = Me%DDecomp%HaloMap

            Out3D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB) =  &
                In3D(HaloMap%ILB:HaloMap%IUB,HaloMap%JLB:HaloMap%JUB, KLB:KUB)

            do i=1, Me%DDecomp%Nslaves

                WorkSize =  Me%DDecomp%Slaves_Size   (i)
                HaloMap  =  Me%DDecomp%Slaves_HaloMap(i)

                allocate(Aux3D (WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB))

                Aux3D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB) = &
                    In3D(HaloMap%ILB:HaloMap%IUB, HaloMap%JLB:HaloMap%JUB, KLB:KUB)

                iSize       = (WorkSize%IUB-WorkSize%ILB+1) * (WorkSize%JUB-WorkSize%JLB+1)*(KUB-KLB+1)

                Precision   = MPI_DOUBLE_PRECISION

                Destination =  Me%DDecomp%Slaves_MPI_ID(i)

                call MPI_Send (Aux3D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB), iSize, Precision,   &
                               Destination, 90018, MPI_COMM_WORLD, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'BroadcastGridData3D_R4_In - ModuleHorizontalGrid - ERR10'

                deallocate(Aux3D)
                nullify   (Aux3D)

            enddo

        else
            WorkSize    = Me%WorkSize

            allocate(Aux3D (WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB))


            iSize       = (WorkSize%IUB-WorkSize%ILB+1) * (WorkSize%JUB-WorkSize%JLB+1) * (KUB-KLB+1)
            Precision   = MPI_DOUBLE_PRECISION
            Source      = Me%DDecomp%Master_MPI_ID

            call MPI_Recv (Aux3D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB), iSize, Precision,   &
                           Source, 90018, MPI_COMM_WORLD, status, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'BroadcastGridData3D_R4_In - ModuleHorizontalGrid - ERR20'

            Out3D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB) = &
                Aux3D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB)

            deallocate(Aux3D)
            nullify   (Aux3D)

        endif

    end subroutine BroadcastGridData3D_R4_In

    !-----------------------------------------------------------------------------------

    subroutine BroadcastGridData3D_R8(HorizontalGridID, In3D, Out3D, KLB, KUB, STAT)


        !Arguments------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real(8),   dimension(:,:,:), pointer        :: In3D
        real(8),   dimension(:,:,:), pointer        :: Out3D
        integer                                     :: KLB, KUB
        integer               , optional            :: STAT

        !Local---------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call BroadcastGridData_In(In3D, Out3D, KLB, KUB)

            STAT_ = SUCCESS_

        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

    end subroutine BroadcastGridData3D_R8

    !-----------------------------------------------------------------------------------


    subroutine BroadcastGridData3D_R8_In(In3D, Out3D, KLB, KUB)


        !Arguments------------------------------------------------------------
        real(8),   dimension(:,:,:), pointer        :: In3D
        real(8),   dimension(:,:,:), pointer        :: Out3D
        integer                                     :: KLB, KUB
        !Local---------------------------------------------------------------
        real(8), dimension(:,:,:), pointer          :: Aux3D
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB

        integer                                     :: Source, Destination
        integer                                     :: iSize, i
        integer, save                               :: Precision
        integer                                     :: status(MPI_STATUS_SIZE)

        type(T_Size2D)                              :: WorkSize, HaloMap

        !Begin---------------------------------------------------------------

        IUB     = Me%DDecomp%Global%IUB
        ILB     = Me%DDecomp%Global%ILB
        JUB     = Me%DDecomp%Global%JUB
        JLB     = Me%DDecomp%Global%JLB

        if (Me%DDecomp%Master) then

            WorkSize = Me%WorkSize
            HaloMap  = Me%DDecomp%HaloMap

            Out3D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB) =  &
                In3D(HaloMap%ILB:HaloMap%IUB,HaloMap%JLB:HaloMap%JUB, KLB:KUB)

            do i=1, Me%DDecomp%Nslaves

                WorkSize =  Me%DDecomp%Slaves_Size   (i)
                HaloMap  =  Me%DDecomp%Slaves_HaloMap(i)

                allocate(Aux3D (WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB))

                Aux3D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB) = &
                    In3D(HaloMap%ILB:HaloMap%IUB, HaloMap%JLB:HaloMap%JUB, KLB:KUB)

                iSize       = (WorkSize%IUB-WorkSize%ILB+1) * (WorkSize%JUB-WorkSize%JLB+1)*(KUB-KLB+1)

                Precision   = MPI_DOUBLE_PRECISION

                Destination =  Me%DDecomp%Slaves_MPI_ID(i)

                call MPI_Send (Aux3D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB), iSize, Precision,   &
                               Destination, 90019, MPI_COMM_WORLD, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'BroadcastGridData3D_R8_In - ModuleHorizontalGrid - ERR10'

                deallocate(Aux3D)
                nullify   (Aux3D)

            enddo

        else
            WorkSize    = Me%WorkSize

            allocate(Aux3D (WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB))


            iSize       = (WorkSize%IUB-WorkSize%ILB+1) * (WorkSize%JUB-WorkSize%JLB+1) * (KUB-KLB+1)
            Precision   = MPI_DOUBLE_PRECISION
            Source      = Me%DDecomp%Master_MPI_ID

            call MPI_Recv (Aux3D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB), iSize, Precision,   &
                           Source, 90019, MPI_COMM_WORLD, status, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'BroadcastGridData3D_R8_In - ModuleHorizontalGrid - ERR20'

            Out3D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB) = &
                Aux3D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, KLB:KUB)

            deallocate(Aux3D)
            nullify   (Aux3D)

        endif

    end subroutine BroadcastGridData3D_R8_In

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! This subroutine computes the thomaz 2D in parallel following a                       !
    ! domain decomposition approach                                                        !
    !                                                                                      !
    ! Input : Coefficients of the  linear system equation                                  !
    ! OutPut: Water level                                                                  !
    ! Author: Paulo Chambel (2013/2)                                                       !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine THOMAS_2D_DDecompHorizGrid(HorizontalGridID, DCoef_2D, FCoef_2D,         &
                                          TiCoef_2D, ECoef_2D, Results_2D, di, dj, STAT)

        !Arguments------------------------------------------------------------
        integer                          :: HorizontalGridID
        real,    dimension(:,:), pointer :: DCoef_2D, FCoef_2D, TiCoef_2D
        real(8), dimension(:,:), pointer :: ECoef_2D
        real,    dimension(:,:), pointer :: Results_2D
        integer                          :: di,    dj
        integer, optional                :: STAT

        !Local---------------------------------------------------------------
        integer                          :: STAT_, ready_
        real(8), pointer, dimension(:)   :: VECG
        real(8), pointer, dimension(:)   :: VECW
        real,    pointer, dimension(:,:) :: R2D
        real,    pointer, dimension(:,:) :: D
        real(8), pointer, dimension(:,:) :: E
        real,    pointer, dimension(:,:) :: F
        real,    pointer, dimension(:,:) :: Ti
        integer                          :: IUB, ILB, JUB, JLB
        integer                          :: IJmin, IJmax
        integer                          :: JImin, JImax

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call AllocateMasterCoef2D

            !Joins all the sub-domain coefficients in a global domain
            call JoinGridData_In(In2D = DCoef_2D,  Out2D = Me%DDecomp%Coef2D%D , dj_= 0, di_ = 0)
            call JoinGridData_In(In2D = ECoef_2D,  Out2D = Me%DDecomp%Coef2D%E , dj_= 0, di_ = 0)
            call JoinGridData_In(In2D = FCoef_2D,  Out2D = Me%DDecomp%Coef2D%F , dj_= 0, di_ = 0)
            call JoinGridData_In(In2D = TiCoef_2D, Out2D = Me%DDecomp%Coef2D%Ti, dj_= 0, di_ = 0)

            nullify(R2D)

if2:        if (Me%DDecomp%Master) then

                IUB = Me%DDecomp%Global%IUB
                ILB = Me%DDecomp%Global%ILB
                JUB = Me%DDecomp%Global%JUB
                JLB = Me%DDecomp%Global%JLB

                D   => Me%DDecomp%Coef2D%D
                E   => Me%DDecomp%Coef2D%E
                F   => Me%DDecomp%Coef2D%F
                Ti  => Me%DDecomp%Coef2D%Ti
                R2D => Me%DDecomp%Coef2D%Results2D

                VECG    => Me%DDecomp%Coef2D%VECG
                VECW    => Me%DDecomp%Coef2D%VECW

                IJmin = ILB * dj + JLB * di
                IJmax = IUB * dj + JUB * di

                JImin = ILB * di + JLB * dj
                JImax = IUB * di + JUB * dj

                !call THOMAS_2D(IJmin, IJmax, JImin, JImax, di, dj, D, E, F, Ti, R2D, VECG, VECW)

                call THOMAS_2D(IJmin       = IJmin,                                    &
                               IJmax       = IJmax,                                    &
                               JImin       = JImin,                                    &
                               JImax       = JImax,                                    &
                               di          = di,                                       &
                               dj          = dj,                                       &
                               DCoef_2D    = D,                                        &
                               ECoef_2D    = E,                                        &
                               FCoef_2D    = F,                                        &
                               TiCoef_2D   = Ti,                                       &
                               ANSWER      = R2D,                                      &
                               VECG        = VECG,                                     &
                               VECW        = VECW)

                nullify(D, E, F, Ti, VECG, VECW)

            endif if2

            !Send for each the water level result
            call BroadcastGridData_In(In2D = R2D, Out2D = Results_2D)

            nullify(R2D)

            STAT_ = SUCCESS_

        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_


    end subroutine THOMAS_2D_DDecompHorizGrid

    !------------------------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! This subroutine computes the thomaz 3D in parallel following a                       !
    ! domain decomposition approach                                                        !
    !                                                                                      !
    ! Input : Coefficients of the  linear system equation                                  !
    ! OutPut: Water level                                                                  !
    ! Author: Paulo Chambel (2013/2)                                                       !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine THOMAS_3D_DDecompHorizGrid(HorizontalGridID, DCoef_3D, FCoef_3D,         &
                                          TiCoef_3D, ECoef_3D, Results_3D, di, dj,      &
                                          KLB, KUB, STAT)

        !Arguments------------------------------------------------------------
        integer                             :: HorizontalGridID
        real,    dimension(:,:,:), pointer  :: DCoef_3D, FCoef_3D, TiCoef_3D
        real(8), dimension(:,:,:), pointer  :: ECoef_3D
        real,    dimension(:,:,:), pointer  :: Results_3D
        integer                             :: di,    dj
        integer                             :: KLB, KUB
        integer, optional                   :: STAT

        !Local---------------------------------------------------------------
        integer                             :: STAT_, ready_
        real(8), pointer, dimension(:)      :: VECG
        real(8), pointer, dimension(:)      :: VECW
        real,    pointer, dimension(:,:,:)  :: R3D
        real,    pointer, dimension(:,:,:)  :: D
        real(8), pointer, dimension(:,:,:)  :: E
        real,    pointer, dimension(:,:,:)  :: F
        real,    pointer, dimension(:,:,:)  :: Ti
        integer                             :: IUB, ILB, JUB, JLB
        integer                             :: IJmin, IJmax
        integer                             :: JImin, JImax

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Me%DDecomp%Coef3D%KLB = KLB
            Me%DDecomp%Coef3D%KUB = KUB

            call AllocateMasterCoef3D

            !Joins all the sub-domain coefficients in a global domain
            call JoinGridData_In(In3D = DCoef_3D,  Out3D = Me%DDecomp%Coef3D%D , KLB = KLB, KUB = KUB, dj_= 0, di_ = 0)
            call JoinGridData_In(In3D = ECoef_3D,  Out3D = Me%DDecomp%Coef3D%E , KLB = KLB, KUB = KUB, dj_= 0, di_ = 0)
            call JoinGridData_In(In3D = FCoef_3D,  Out3D = Me%DDecomp%Coef3D%F , KLB = KLB, KUB = KUB, dj_= 0, di_ = 0)
            call JoinGridData_In(In3D = TiCoef_3D, Out3D = Me%DDecomp%Coef3D%Ti, KLB = KLB, KUB = KUB, dj_= 0, di_ = 0)

            nullify(R3D)

if2:        if (Me%DDecomp%Master) then

                IUB = Me%DDecomp%Global%IUB
                ILB = Me%DDecomp%Global%ILB
                JUB = Me%DDecomp%Global%JUB
                JLB = Me%DDecomp%Global%JLB

                D   => Me%DDecomp%Coef3D%D
                E   => Me%DDecomp%Coef3D%E
                F   => Me%DDecomp%Coef3D%F
                Ti  => Me%DDecomp%Coef3D%Ti
                R3D => Me%DDecomp%Coef3D%Results3D

                VECG    => Me%DDecomp%Coef3D%VECG
                VECW    => Me%DDecomp%Coef3D%VECW

                IJmin = ILB * dj + JLB * di
                IJmax = IUB * dj + JUB * di

                JImin = ILB * di + JLB * dj
                JImax = IUB * di + JUB * dj

                call THOMAS_3D(IJmin       = IJmin,                                    &
                               IJmax       = IJmax,                                    &
                               JImin       = JImin,                                    &
                               JImax       = JImax,                                    &
                               Kmin        = KLB,                                      &
                               Kmax        = KUB,                                      &
                               di          = di,                                       &
                               dj          = dj,                                       &
                               DCoef_3D    = D,                                        &
                               ECoef_3D    = E,                                        &
                               FCoef_3D    = F,                                        &
                               TiCoef_3D   = Ti,                                       &
                               ANSWER      = R3D,                                      &
                               VECG        = VECG,                                     &
                               VECW        = VECW)

                nullify(D, E, F, Ti, VECG, VECW)

            endif if2

            !Send for each the water level result
            call BroadcastGridData_In(In3D = R3D, Out3D = Results_3D, KLB = KLB, KUB = KUB)

            nullify(R3D)

            STAT_ = SUCCESS_

        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_


    end subroutine THOMAS_3D_DDecompHorizGrid

    !------------------------------------------------------------------------------


    subroutine AllocateMasterCoef2D

        !Arguments------------------------------------------------------------

        !Local---------------------------------------------------------------
        integer                          :: IUB, ILB, JUB, JLB
        integer                          :: IUS, ILS, JUS, JLS
        integer                          :: JImin, JImax

        !Begin---------------------------------------------------------------

        if  (.not. Me%DDecomp%Coef2D%AllocateON .and. Me%DDecomp%Master) then

            IUB = Me%DDecomp%Global%IUB
            ILB = Me%DDecomp%Global%ILB
            JUB = Me%DDecomp%Global%JUB
            JLB = Me%DDecomp%Global%JLB

            ILS = ILB-1
            IUS = IUB+1
            JLS = JLB-1
            JUS = JUB+1

            write(*,*) 'ILS, IUS, JLS, JUS', ILS, IUS, JLS, JUS

            allocate(Me%DDecomp%Coef2D%D        (ILS:IUS, JLS:JUS))
            allocate(Me%DDecomp%Coef2D%E        (ILS:IUS, JLS:JUS))
            allocate(Me%DDecomp%Coef2D%F        (ILS:IUS, JLS:JUS))
            allocate(Me%DDecomp%Coef2D%Ti       (ILS:IUS, JLS:JUS))
            allocate(Me%DDecomp%Coef2D%Results2D(ILS:IUS, JLS:JUS))

            Me%DDecomp%Coef2D%D        (:,:) = 0.
            Me%DDecomp%Coef2D%E        (:,:) = 1.
            Me%DDecomp%Coef2D%F        (:,:) = 0.
            Me%DDecomp%Coef2D%Ti       (:,:) = 0.
            Me%DDecomp%Coef2D%Results2D(:,:) = 0.


            JImin = MIN(ILS,JLS)
            JImax = MAX(IUS,JUS)

            allocate(Me%DDecomp%Coef2D%VECG(JImin:JImax))
            allocate(Me%DDecomp%Coef2D%VECW(JImin:JImax))

            Me%DDecomp%Coef2D%VECG(:) = 0
            Me%DDecomp%Coef2D%VECW(:) = 0

            Me%DDecomp%Coef2D%AllocateON = .true.

        endif

    end subroutine AllocateMasterCoef2D



    subroutine AllocateMasterCoef3D

        !Arguments------------------------------------------------------------

        !Local---------------------------------------------------------------
        integer                          :: IUB, ILB, JUB, JLB
        integer                          :: IUS, ILS, JUS, JLS
        integer                          :: JImin, JImax
        integer                          :: KLB, KUB, KLS, KUS

        !Begin---------------------------------------------------------------

        if  (.not. Me%DDecomp%Coef3D%AllocateON .and. Me%DDecomp%Master) then

            IUB = Me%DDecomp%Global%IUB
            ILB = Me%DDecomp%Global%ILB
            JUB = Me%DDecomp%Global%JUB
            JLB = Me%DDecomp%Global%JLB

            KLB = Me%DDecomp%Coef3D%KLB
            KUB = Me%DDecomp%Coef3D%KUB

            ILS = ILB-1
            IUS = IUB+1
            JLS = JLB-1
            JUS = JUB+1
            KLS = KLB-1
            KUS = KUB+1

            allocate(Me%DDecomp%Coef3D%D        (ILS:IUS, JLS:JUS, KLS:KUS))
            allocate(Me%DDecomp%Coef3D%E        (ILS:IUS, JLS:JUS, KLS:KUS))
            allocate(Me%DDecomp%Coef3D%F        (ILS:IUS, JLS:JUS, KLS:KUS))
            allocate(Me%DDecomp%Coef3D%Ti       (ILS:IUS, JLS:JUS, KLS:KUS))
            allocate(Me%DDecomp%Coef3D%Results3D(ILS:IUS, JLS:JUS, KLS:KUS))

            Me%DDecomp%Coef3D%D        (:,:,:) = 0.
            Me%DDecomp%Coef3D%E        (:,:,:) = 1.
            Me%DDecomp%Coef3D%F        (:,:,:) = 0.
            Me%DDecomp%Coef3D%Ti       (:,:,:) = 0.
            Me%DDecomp%Coef3D%Results3D(:,:,:) = 0.


            JImin = MIN(ILS,JLS)
            JImax = MAX(IUS,JUS)

            allocate(Me%DDecomp%Coef3D%VECG(JImin:JImax))
            allocate(Me%DDecomp%Coef3D%VECW(JImin:JImax))

            Me%DDecomp%Coef3D%VECG(:) = 0
            Me%DDecomp%Coef3D%VECW(:) = 0

            Me%DDecomp%Coef3D%AllocateON = .true.

        endif

    end subroutine AllocateMasterCoef3D


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! This subroutine exchange 3D properties between decompose domains                     !
    ! coefficients of all domains                                                          !
    !                                                                                      !
    ! Input : Property3D halo region                                                       !
    ! OutPut: Property3D halo region - boundary domains                                    !
    ! Author: Paulo Chambel (2015/11)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ReceiveSendProperities3DMPIr8(HorizontalGridID, Property3D, KLB, KUB, di, dj, STAT)

        !Arguments------------------------------------------------------------
        integer                            :: HorizontalGridID
        real(8), dimension(:,:,:), pointer :: Property3D
        integer                            :: KLB, KUB
        integer                 , optional :: di, dj
        integer                 , optional :: STAT
        !Local---------------------------------------------------------------
        integer                            :: STAT_CALL, IUB, ILB, JUB, JLB
        integer                            :: IUB_R, ILB_R, JUB_R, JLB_R
        integer                            :: IUB_S, ILB_S, JUB_S, JLB_S
        integer                            :: Bandwidth
        integer                            :: DomainA, DomainB, ifd, Direction

        integer                            :: Source, Destination
        integer                            :: iSize
        integer, save                      :: Precision
        integer                            :: status(MPI_STATUS_SIZE)
        integer                            :: di_, dj_
        integer                            :: STAT_, ready_

        !Begin---------------------------------------------------------------

       STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
    
    idd:    if (Me%DDecomp%MasterOrSlave) then

                IUB = Me%WorkSize%IUB
                ILB = Me%WorkSize%ILB
                JUB = Me%WorkSize%JUB
                JLB = Me%WorkSize%JLB

                if (present(di)) then
                    di_ = di
                else
                    di_ = 0
                endif

                if (present(dj)) then
                    dj_ = dj
                else
                    dj_ = 0
                endif

                !Precision   = MPIKind(Property3D)
                Precision = MPI_DOUBLE_PRECISION

        difd:   do ifd = 1, Me%DDecomp%NInterfaces

                    DomainA   = Me%DDecomp%Interfaces(ifd,1)
                    DomainB   = Me%DDecomp%Interfaces(ifd,2)
                    Direction = Me%DDecomp%Interfaces(ifd,3)

        iSN:        if (Direction == SouthNorth_) then

                        !Then North border communication
        iN:             if (Me%DDecomp%MPI_ID == DomainA) then

                            Bandwidth   = Me%DDecomp%Halo_Points + di_
                            iSize       = Bandwidth * (JUB-JLB+1) * (KUB-KLB+1)
                            Source      = DomainB
                            Destination = Source

                            IUB_R = IUB                     + di_
                            ILB_R = IUB_R - Bandwidth   + 1
                            IUB_S = ILB_R               - 1 + di_
                            ILB_S = IUB_S - Bandwidth   + 1

                            !Receive
                            call MPI_Recv (Property3D(ILB_R:IUB_R, JLB:JUB, KLB:KUB), iSize, Precision,      &
                                           Source, 180001, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR20'

                            !Send
                            call MPI_Send (Property3D(ILB_S:IUB_S, JLB:JUB, KLB:KUB), iSize, Precision,      &
                                           Destination, 180002, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR30'

                        endif iN

                        !Then South border communication
        iS:             if (Me%DDecomp%MPI_ID == DomainB) then

                            Bandwidth   = Me%DDecomp%Halo_Points + di_
                            iSize       = Bandwidth * (JUB-JLB+1) * (KUB-KLB+1)
                            Source      = DomainA
                            Destination = Source

                            ILB_R = ILB
                            IUB_R = ILB_R + Bandwidth - 1
                            ILB_S = IUB_R             + 1 - di_
                            IUB_S = ILB_S + Bandwidth - 1

                            !Send
                            call MPI_Send (Property3D(ILB_S:IUB_S, JLB:JUB, KLB:KUB), iSize, Precision,      &
                                           Destination, 180001, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR40'

                            !Receive
                            call MPI_Recv (Property3D(ILB_R:IUB_R, JLB:JUB, KLB:KUB), iSize, Precision,      &
                                           Source, 180002, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR50'

                        endif iS

                    endif iSN

        iWE:        if (Direction == WestEast_) then

                        !Then East border communication
        iE:             if (Me%DDecomp%MPI_ID == DomainA) then

                            Bandwidth   = Me%DDecomp%Halo_Points + dj_
                            iSize       = Bandwidth * (IUB-ILB+1) * (KUB-KLB+1)
                            Source      = DomainB
                            Destination = Source

                            !Receive
                            JUB_R = JUB                     + dj_
                            JLB_R = JUB_R - Bandwidth   + 1
                            JUB_S = JLB_R               - 1 + dj_
                            JLB_S = JUB_S - Bandwidth   + 1

                            call MPI_Recv (Property3D(ILB:IUB, JLB_R:JUB_R, KLB:KUB), iSize, Precision,      &
                                           Source, 180005, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR60'

                            !Send
                            call MPI_Send (Property3D(ILB:IUB, JLB_S:JUB_S, KLB:KUB), iSize, Precision,      &
                                           Destination, 180006, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR70'

                        endif iE

                        !Then West border communication
        iW:             if (Me%DDecomp%MPI_ID == DomainB) then

                            Bandwidth   = Me%DDecomp%Halo_Points + dj_
                            iSize       = Bandwidth * (IUB-ILB+1) * (KUB-KLB+1)
                            Source      = DomainA
                            Destination = Source

                            !Receive
                            JLB_R = JLB
                            JUB_R = JLB_R + Bandwidth - 1
                            JLB_S = JUB_R             + 1 - dj_
                            JUB_S = JLB_S + Bandwidth - 1

                            !Send
                            call MPI_Send (Property3D(ILB:IUB, JLB_S:JUB_S, KLB:KUB), iSize, Precision,      &
                                           Destination, 180005, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR80'

                            call MPI_Recv (Property3D(ILB:IUB, JLB_R:JUB_R, KLB:KUB), iSize, Precision,      &
                                           Source, 180006, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR90'

                        endif iW

                    endif iWE

                enddo difd

            endif idd
            
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


    end subroutine ReceiveSendProperities3DMPIr8

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! This subroutine exchange 3D properties between decompose domains                     !
    ! coefficients of all domains                                                          !
    !                                                                                      !
    ! Input : Property3D halo region                                                       !
    ! OutPut: Property3D halo region - boundary domains                                    !
    ! Author: Paulo Chambel (2015/11)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ReceiveSendProperities3DMPIr4(HorizontalGridID, Property3D, KLB, KUB, di, dj, STAT)

        !Arguments------------------------------------------------------------
        integer                            :: HorizontalGridID
        real(4), dimension(:,:,:), pointer :: Property3D
        integer                            :: KLB, KUB
        integer                 , optional :: di, dj
        integer                 , optional :: STAT
        !Local---------------------------------------------------------------
        integer                            :: STAT_CALL, IUB, ILB, JUB, JLB
        integer                            :: IUB_R, ILB_R, JUB_R, JLB_R
        integer                            :: IUB_S, ILB_S, JUB_S, JLB_S
        integer                            :: Bandwidth
        integer                            :: DomainA, DomainB, ifd, Direction

        integer                            :: Source, Destination
        integer                            :: iSize
        integer, save                      :: Precision
        integer                            :: status(MPI_STATUS_SIZE)
        integer                            :: di_, dj_
        integer                            :: STAT_, ready_

        !Begin---------------------------------------------------------------

       STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then


    idd:    if (Me%DDecomp%MasterOrSlave) then

                IUB = Me%WorkSize%IUB
                ILB = Me%WorkSize%ILB
                JUB = Me%WorkSize%JUB
                JLB = Me%WorkSize%JLB

                if (present(di)) then
                    di_ = di
                else
                    di_ = 0
                endif

                if (present(dj)) then
                    dj_ = dj
                else
                    dj_ = 0
                endif

                !Precision   = MPIKind(Property3D)
                Precision = MPI_REAL

        difd:   do ifd = 1, Me%DDecomp%NInterfaces

                    DomainA   = Me%DDecomp%Interfaces(ifd,1)
                    DomainB   = Me%DDecomp%Interfaces(ifd,2)
                    Direction = Me%DDecomp%Interfaces(ifd,3)

        iSN:        if (Direction == SouthNorth_) then

                        !Then North border communication
        iN:             if (Me%DDecomp%MPI_ID == DomainA) then

                            Bandwidth   = Me%DDecomp%Halo_Points + di_
                            iSize       = Bandwidth * (JUB-JLB+1) * (KUB-KLB+1)
                            Source      = DomainB
                            Destination = Source

                            IUB_R = IUB                     + di_
                            ILB_R = IUB_R - Bandwidth   + 1
                            IUB_S = ILB_R               - 1 + di_
                            ILB_S = IUB_S - Bandwidth   + 1

                            !Receive
                            call MPI_Recv (Property3D(ILB_R:IUB_R, JLB:JUB, KLB:KUB), iSize, Precision,      &
                                           Source, 180001, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR20'

                            !Send
                            call MPI_Send (Property3D(ILB_S:IUB_S, JLB:JUB, KLB:KUB), iSize, Precision,      &
                                           Destination, 180002, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR30'

                        endif iN

                        !Then South border communication
        iS:             if (Me%DDecomp%MPI_ID == DomainB) then

                            Bandwidth   = Me%DDecomp%Halo_Points + di_
                            iSize       = Bandwidth * (JUB-JLB+1) * (KUB-KLB+1)
                            Source      = DomainA
                            Destination = Source

                            ILB_R = ILB
                            IUB_R = ILB_R + Bandwidth - 1
                            ILB_S = IUB_R             + 1 - di_
                            IUB_S = ILB_S + Bandwidth - 1

                            !Send
                            call MPI_Send (Property3D(ILB_S:IUB_S, JLB:JUB, KLB:KUB), iSize, Precision,      &
                                           Destination, 180001, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR40'

                            !Receive
                            call MPI_Recv (Property3D(ILB_R:IUB_R, JLB:JUB, KLB:KUB), iSize, Precision,      &
                                           Source, 180002, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR50'

                        endif iS

                    endif iSN

        iWE:        if (Direction == WestEast_) then

                        !Then East border communication
        iE:             if (Me%DDecomp%MPI_ID == DomainA) then

                            Bandwidth   = Me%DDecomp%Halo_Points + dj_
                            iSize       = Bandwidth * (IUB-ILB+1) * (KUB-KLB+1)
                            Source      = DomainB
                            Destination = Source

                            !Receive
                            JUB_R = JUB                     + dj_
                            JLB_R = JUB_R - Bandwidth   + 1
                            JUB_S = JLB_R               - 1 + dj_
                            JLB_S = JUB_S - Bandwidth   + 1

                            call MPI_Recv (Property3D(ILB:IUB, JLB_R:JUB_R, KLB:KUB), iSize, Precision,      &
                                           Source, 180005, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR60'

                            !Send
                            call MPI_Send (Property3D(ILB:IUB, JLB_S:JUB_S, KLB:KUB), iSize, Precision,      &
                                           Destination, 180006, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR70'

                        endif iE

                        !Then West border communication
        iW:             if (Me%DDecomp%MPI_ID == DomainB) then

                            Bandwidth   = Me%DDecomp%Halo_Points + dj_
                            iSize       = Bandwidth * (IUB-ILB+1) * (KUB-KLB+1)
                            Source      = DomainA
                            Destination = Source

                            !Receive
                            JLB_R = JLB
                            JUB_R = JLB_R + Bandwidth - 1
                            JLB_S = JUB_R             + 1 - dj_
                            JUB_S = JLB_S + Bandwidth - 1

                            !Send
                            call MPI_Send (Property3D(ILB:IUB, JLB_S:JUB_S, KLB:KUB), iSize, Precision,      &
                                           Destination, 180005, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR80'

                            call MPI_Recv (Property3D(ILB:IUB, JLB_R:JUB_R, KLB:KUB), iSize, Precision,      &
                                           Source, 180006, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities3DMPI - ModuleHorizontalGrid - ERR90'

                        endif iW

                    endif iWE

                enddo difd

            endif idd

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


    end subroutine ReceiveSendProperities3DMPIr4



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! This subroutine exchange 2D properties between decompose domains                     !
    ! coefficients of all domains                                                          !
    !                                                                                      !
    ! Input : Property3D halo region                                                       !
    ! OutPut: Property3D halo region - boundary domains                                    !
    ! Author: Paulo Chambel (2015/11)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ReceiveSendProperities2DMPIr8(HorizontalGridID, Property2D, STAT)

        !Arguments------------------------------------------------------------
        integer                            :: HorizontalGridID
        real(8), dimension(:,:)  , pointer :: Property2D
        integer                 , optional :: STAT
        !Local---------------------------------------------------------------
        integer                            :: STAT_CALL, IUB, ILB, JUB, JLB
        integer                            :: IUB_R, ILB_R, JUB_R, JLB_R
        integer                            :: IUB_S, ILB_S, JUB_S, JLB_S
        integer                            :: Bandwidth
        integer                            :: DomainA, DomainB, ifd, Direction

        integer                            :: Source, Destination
        integer                            :: iSize
        integer, save                      :: Precision
        integer                            :: status(MPI_STATUS_SIZE)
        integer                            :: STAT_, ready_

        !Begin---------------------------------------------------------------

       STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then


    idd:    if (Me%DDecomp%MasterOrSlave) then

                IUB = Me%WorkSize%IUB
                ILB = Me%WorkSize%ILB
                JUB = Me%WorkSize%JUB
                JLB = Me%WorkSize%JLB

                !Precision   = MPIKind(Property2D)
                Precision = MPI_DOUBLE_PRECISION


        difd:   do ifd = 1, Me%DDecomp%NInterfaces

                    DomainA   = Me%DDecomp%Interfaces(ifd,1)
                    DomainB   = Me%DDecomp%Interfaces(ifd,2)
                    Direction = Me%DDecomp%Interfaces(ifd,3)

        iSN:        if (Direction == SouthNorth_) then

                        !Then North border communication
        iN:             if (Me%DDecomp%MPI_ID == DomainA) then

                            Bandwidth   = Me%DDecomp%Halo_Points
                            iSize       = Bandwidth * (JUB-JLB+1)
                            Source      = DomainB
                            Destination = Source

                            IUB_R = IUB
                            ILB_R = IUB_R - Bandwidth   + 1
                            IUB_S = ILB_R               - 1
                            ILB_S = IUB_S - Bandwidth   + 1

                            !Receive
                            call MPI_Recv (Property2D(ILB_R:IUB_R, JLB:JUB), iSize, Precision,      &
                                           Source, 281001, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR20'

                            !Send
                            call MPI_Send (Property2D(ILB_S:IUB_S, JLB:JUB), iSize, Precision,      &
                                           Destination, 281002, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR30'

                        endif iN

                        !Then South border communication
        iS:             if (Me%DDecomp%MPI_ID == DomainB) then

                            Bandwidth   = Me%DDecomp%Halo_Points
                            iSize       = Bandwidth * (JUB-JLB+1)
                            Source      = DomainA
                            Destination = Source

                            ILB_R = ILB
                            IUB_R = ILB_R + Bandwidth - 1
                            ILB_S = IUB_R             + 1
                            IUB_S = ILB_S + Bandwidth - 1

                            !Send
                            call MPI_Send (Property2D(ILB_S:IUB_S, JLB:JUB), iSize, Precision,      &
                                           Destination, 281001, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR40'

                            !Receive
                            call MPI_Recv (Property2D(ILB_R:IUB_R, JLB:JUB), iSize, Precision,      &
                                           Source, 281002, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR50'

                        endif iS

                    endif iSN

        iWE:        if (Direction == WestEast_) then

                        !Then East border communication
        iE:             if (Me%DDecomp%MPI_ID == DomainA) then

                            Bandwidth   = Me%DDecomp%Halo_Points
                            iSize       = Bandwidth * (IUB-ILB+1)
                            Source      = DomainB
                            Destination = Source

                            !Receive
                            JUB_R = JUB
                            JLB_R = JUB_R - Bandwidth   + 1
                            JUB_S = JLB_R               - 1
                            JLB_S = JUB_S - Bandwidth   + 1

                            call MPI_Recv (Property2D(ILB:IUB, JLB_R:JUB_R), iSize, Precision, &
                                           Source, 281005, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR60'

                            !Send
                            call MPI_Send (Property2D(ILB:IUB, JLB_S:JUB_S), iSize, Precision,      &
                                           Destination, 281006, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR70'

                        endif iE

                        !Then West border communication
        iW:             if (Me%DDecomp%MPI_ID == DomainB) then

                            Bandwidth   = Me%DDecomp%Halo_Points
                            iSize       = Bandwidth * (IUB-ILB+1)
                            Source      = DomainA
                            Destination = Source

                            !Receive
                            JLB_R = JLB
                            JUB_R = JLB_R + Bandwidth - 1
                            JLB_S = JUB_R             + 1
                            JUB_S = JLB_S + Bandwidth - 1

                            !Send
                            call MPI_Send (Property2D(ILB:IUB, JLB_S:JUB_S), iSize, Precision,      &
                                           Destination, 281005, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR80'

                            call MPI_Recv (Property2D(ILB:IUB, JLB_R:JUB_R), iSize, Precision,      &
                                           Source, 281006, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR90'

                        endif iW

                    endif iWE

                enddo difd

            endif idd
            
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


    end subroutine ReceiveSendProperities2DMPIr8

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! This subroutine exchange 2D properties between decompose domains                     !
    ! coefficients of all domains                                                          !
    !                                                                                      !
    ! Input : Property3D halo region                                                       !
    ! OutPut: Property3D halo region - boundary domains                                    !
    ! Author: Paulo Chambel (2015/11)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ReceiveSendProperities2DMPIr4(HorizontalGridID, Property2D, STAT)

        !Arguments------------------------------------------------------------
        integer                            :: HorizontalGridID
        real(4), dimension(:,:)  , pointer :: Property2D
        integer                 , optional :: STAT
        !Local---------------------------------------------------------------
        integer                            :: STAT_CALL, IUB, ILB, JUB, JLB
        integer                            :: IUB_R, ILB_R, JUB_R, JLB_R
        integer                            :: IUB_S, ILB_S, JUB_S, JLB_S
        integer                            :: Bandwidth
        integer                            :: DomainA, DomainB, ifd, Direction

        integer                            :: Source, Destination
        integer                            :: iSize
        integer, save                      :: Precision
        integer                            :: status(MPI_STATUS_SIZE)
        integer                            :: STAT_, ready_

        !Begin---------------------------------------------------------------

       STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then


    idd:    if (Me%DDecomp%MasterOrSlave) then

                IUB = Me%WorkSize%IUB
                ILB = Me%WorkSize%ILB
                JUB = Me%WorkSize%JUB
                JLB = Me%WorkSize%JLB

                !Precision   = MPIKind(Property2D)
                Precision = MPI_REAL

        difd:   do ifd = 1, Me%DDecomp%NInterfaces

                    DomainA   = Me%DDecomp%Interfaces(ifd,1)
                    DomainB   = Me%DDecomp%Interfaces(ifd,2)
                    Direction = Me%DDecomp%Interfaces(ifd,3)

        iSN:        if (Direction == SouthNorth_) then

                        !Then North border communication
        iN:             if (Me%DDecomp%MPI_ID == DomainA) then

                            Bandwidth   = Me%DDecomp%Halo_Points
                            iSize       = Bandwidth * (JUB-JLB+1)
                            Source      = DomainB
                            Destination = Source

                            IUB_R = IUB
                            ILB_R = IUB_R - Bandwidth   + 1
                            IUB_S = ILB_R               - 1
                            ILB_S = IUB_S - Bandwidth   + 1

                            !Receive
                            call MPI_Recv (Property2D(ILB_R:IUB_R, JLB:JUB), iSize, Precision,      &
                                           Source, 281001, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR20'

                            !Send
                            call MPI_Send (Property2D(ILB_S:IUB_S, JLB:JUB), iSize, Precision,      &
                                           Destination, 281002, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR30'

                        endif iN

                        !Then South border communication
        iS:             if (Me%DDecomp%MPI_ID == DomainB) then

                            Bandwidth   = Me%DDecomp%Halo_Points
                            iSize       = Bandwidth * (JUB-JLB+1)
                            Source      = DomainA
                            Destination = Source

                            ILB_R = ILB
                            IUB_R = ILB_R + Bandwidth - 1
                            ILB_S = IUB_R             + 1
                            IUB_S = ILB_S + Bandwidth - 1

                            !Send
                            call MPI_Send (Property2D(ILB_S:IUB_S, JLB:JUB), iSize, Precision,      &
                                           Destination, 281001, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR40'

                            !Receive
                            call MPI_Recv (Property2D(ILB_R:IUB_R, JLB:JUB), iSize, Precision,      &
                                           Source, 281002, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR50'

                        endif iS

                    endif iSN

        iWE:        if (Direction == WestEast_) then

                        !Then East border communication
        iE:             if (Me%DDecomp%MPI_ID == DomainA) then

                            Bandwidth   = Me%DDecomp%Halo_Points
                            iSize       = Bandwidth * (IUB-ILB+1)
                            Source      = DomainB
                            Destination = Source

                            !Receive
                            JUB_R = JUB
                            JLB_R = JUB_R - Bandwidth   + 1
                            JUB_S = JLB_R               - 1
                            JLB_S = JUB_S - Bandwidth   + 1

                            call MPI_Recv (Property2D(ILB:IUB, JLB_R:JUB_R), iSize, Precision, &
                                           Source, 281005, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR60'

                            !Send
                            call MPI_Send (Property2D(ILB:IUB, JLB_S:JUB_S), iSize, Precision,      &
                                           Destination, 281006, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR70'

                        endif iE

                        !Then West border communication
        iW:             if (Me%DDecomp%MPI_ID == DomainB) then

                            Bandwidth   = Me%DDecomp%Halo_Points
                            iSize       = Bandwidth * (IUB-ILB+1)
                            Source      = DomainA
                            Destination = Source

                            !Receive
                            JLB_R = JLB
                            JUB_R = JLB_R + Bandwidth - 1
                            JLB_S = JUB_R             + 1
                            JUB_S = JLB_S + Bandwidth - 1

                            !Send
                            call MPI_Send (Property2D(ILB:IUB, JLB_S:JUB_S), iSize, Precision,      &
                                           Destination, 281005, MPI_COMM_WORLD, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR80'

                            call MPI_Recv (Property2D(ILB:IUB, JLB_R:JUB_R), iSize, Precision,      &
                                           Source, 281006, MPI_COMM_WORLD, status, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ReceiveSendProperities2DMPI - ModuleHorizontalGrid - ERR90'

                        endif iW

                    endif iWE

                enddo difd

            endif idd

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


    end subroutine ReceiveSendProperities2DMPIr4

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! This subroutine send a logical variable from the slaves domains to the master one    !
    ! GridData File                                                                        !
    !                                                                                      !
    ! Input : Property3D halo region                                                       !
    ! OutPut: Property3D halo region - boundary domains                                    !
    ! Author: Paulo Chambel (2015/11)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine AtLeastOneDomainIsTrue(HorizontalGridID, LogicalIn, LogicalOut, STAT)

        !Arguments------------------------------------------------------------
        integer            , intent(IN)    :: HorizontalGridID
        logical            , intent(IN)    :: LogicalIn
        logical            , intent(OUT)   :: LogicalOut
        integer            , optional      :: STAT
        !Local---------------------------------------------------------------
        integer                            :: STAT_CALL, i
        logical                            :: LogicalAux
        integer                            :: Source, Destination
        integer                            :: iSize
        integer, save                      :: Precision
        integer                            :: status(MPI_STATUS_SIZE)
        integer                            :: STAT_, ready_

        !Begin---------------------------------------------------------------

       STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            LogicalOut = .false.

            if (Me%DDecomp%Master) then

                if (LogicalIn) LogicalOut = .true.

                iSize     = 1

                Precision = MPI_LOGICAL


                !Receive from slaves
                do i=1, Me%DDecomp%Nslaves

                    Source    =  Me%DDecomp%Slaves_MPI_ID(i)

                    !Receive logical from slaves
                    call MPI_Recv (LogicalAux, iSize, Precision, Source, 999003, MPI_COMM_WORLD, status, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'AtLeastOneDomainIsTrue - ModuleHorizontalGrid - ERR10'

                    if (LogicalAux) LogicalOut = .true.

                enddo

                !send to slaves
                do i=1, Me%DDecomp%Nslaves

                    Destination  =  Me%DDecomp%Slaves_MPI_ID(i)

                    !Send logical to slaves
                    call MPI_Send (LogicalOut, iSize, Precision, Destination, 999004, MPI_COMM_WORLD, status, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'AtLeastOneDomainIsTrue - ModuleHorizontalGrid - ERR20'

                enddo

            else

                iSize       = 1

                Precision   = MPI_INTEGER

                Destination = Me%DDecomp%Master_MPI_ID

                !Send logical to master
                call MPI_Send (LogicalIn, iSize, Precision, Destination, 999003, MPI_COMM_WORLD, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'AtLeastOneDomainIsTrue - ModuleHorizontalGrid - ERR30'

                Source      = Me%DDecomp%Master_MPI_ID

                !Receive logical from master
                call MPI_Recv (LogicalOut, iSize, Precision, Source, 999004, MPI_COMM_WORLD, status, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'AtLeastOneDomainIsTrue - ModuleHorizontalGrid - ERR40'

            endif

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


    end subroutine AtLeastOneDomainIsTrue
    
    !------------------------------------------------------------------------------

    subroutine GetKfloorZminMPI(HorizontalGridID, KminZin, KminZout, STAT)
    
        !Arguments------------------------------------------------------------
        integer            , intent(IN)    :: HorizontalGridID
        integer            , intent(IN)    :: KminZin        
        integer            , intent(OUT)   :: KminZout
        integer            , optional      :: STAT
        !Local---------------------------------------------------------------
        integer                            :: STAT_CALL, i
        integer                            :: KminZAux
        integer                            :: Source, Destination
        integer                            :: iSize
        integer, save                      :: Precision
        integer                            :: status(MPI_STATUS_SIZE)
        integer                            :: STAT_, ready_

        !Begin---------------------------------------------------------------    

       STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
    
            KminZout = KminZin
            
            if (Me%DDecomp%ON) then

                if (Me%DDecomp%Master) then

                    iSize     = 1

                    Precision = MPI_INTEGER
                
                    !Receive from slaves
                    do i=1, Me%DDecomp%Nslaves

                        Source    =  Me%DDecomp%Slaves_MPI_ID(i)

                        !Receive logical from slaves
                        call MPI_Recv (KminZAux, iSize, Precision, Source, 999903, MPI_COMM_WORLD, status, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'GetKfloorZminMPI - ModuleHorizontalGrid - ERR10'

                        if (KminZAux < KminZout) KminZout = KminZAux

                    enddo

                    !send to slaves
                    do i=1, Me%DDecomp%Nslaves

                        Destination  =  Me%DDecomp%Slaves_MPI_ID(i)

                        !Send logical to slaves
                        call MPI_Send (KminZout, iSize, Precision, Destination, 999904, MPI_COMM_WORLD, status, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'GetKfloorZminMPI - ModuleHorizontalGrid - ERR20'

                    enddo

                else

                    iSize       = 1

                    Precision   = MPI_INTEGER

                    Destination = Me%DDecomp%Master_MPI_ID

                    !Send integer to master
                    call MPI_Send (KminZin, iSize, Precision, Destination, 999903, MPI_COMM_WORLD, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetKfloorZminMPI - ModuleHorizontalGrid - ERR30'

                    Source      = Me%DDecomp%Master_MPI_ID

                    !Receive logical from master
                    call MPI_Recv (KminZout, iSize, Precision, Source, 999904, MPI_COMM_WORLD, status, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetKfloorZminMPI - ModuleHorizontalGrid - ERR40'

                endif
                
            endif                

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetKfloorZminMPI
    !------------------------------------------------------------------------------
    

#endif _USE_MPI

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! This subroutine add the MPI ID to a Filename                                         !
    !                                                                                      !
    ! Input : Filename                                                                     !
    ! OutPut: MPI_X_Filename                                                               !
    ! Author: Paulo Chambel (2015/12)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Add_MPI_ID_2_Filename(HorizontalGridID, FileName, STAT)

        !Arguments------------------------------------------------------------
        integer            , intent(IN)    :: HorizontalGridID
        character(len=*)   , intent(INOUT) :: FileName
        integer            , optional      :: STAT
        !Local---------------------------------------------------------------
        integer                            :: I, ipath, iFN
        character(LEN = StringLength)      :: AuxChar
        integer                            :: STAT_, ready_

        !Begin---------------------------------------------------------------

       STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (Me%DDecomp%ON) then
                if (STAT_ == SUCCESS_) then
                    if (Me%DDecomp%MPI_ID > null_int) then
                        write(AuxChar,fmt='(i5)') Me%DDecomp%MPI_ID
                        Auxchar = "MPI_"//trim(adjustl(Auxchar))//"_"
                        iFN = len_trim(Filename)
                        ipath = 0
                        do i = iFN, 1, -1
                            if (Filename(i:i) == '/' .or. Filename(i:i) == backslash) then
                                ipath = i
                                exit
                            endif
                        enddo
                        if (ipath > 0) then
                            Filename = Filename(1:ipath)//trim(Auxchar)//Filename(ipath+1:iFN)
                        endif
                    endif
                endif
            endif

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


    end subroutine Add_MPI_ID_2_Filename

   !--------------------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! This subroutine assumes null box mapping in the halo area                            !
    !                                                                                      !
    ! Input : Box2D, Box3D                                                                 !
    ! OutPut: Box2D, Box3D - no map in halo area                                           !
    ! Author: Paulo Chambel (2015/12)                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine No_BoxMap_HaloArea(HorizontalGridID, Boxes2D, Boxes3D, KLB, KUB, STAT)

        !Arguments------------------------------------------------------------
        integer                                     , intent(IN)    :: HorizontalGridID
        integer, dimension(:,:  ), pointer, optional, intent(INOUT) :: Boxes2D
        integer, dimension(:,:,:), pointer, optional, intent(INOUT) :: Boxes3D
        integer                           , optional, intent(IN)    :: KLB, KUB
        integer                           , optional, intent(OUT)   :: STAT
        !Local---------------------------------------------------------------
        integer                                                     :: STAT_, ready_
        integer                                                     :: NeighbourSouth, NeighbourWest
        integer                                                     :: NeighbourEast, NeighbourNorth
        integer                                                     :: Halo_Points, WILB, WIUB, WJLB, WJUB

        !Begin---------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (Me%DDecomp%ON) then

                WILB = Me%WorkSize%ILB
                WIUB = Me%WorkSize%IUB
                WJLB = Me%WorkSize%JLB
                WJUB = Me%WorkSize%JUB

                NeighbourSouth = Me%DDecomp%NeighbourSouth
                NeighbourWest  = Me%DDecomp%NeighbourWest
                NeighbourEast  = Me%DDecomp%NeighbourEast
                NeighbourNorth = Me%DDecomp%NeighbourNorth
                Halo_Points    = Me%DDecomp%Halo_Points

                if (present(Boxes2D)) then

                    if (NeighbourSouth /= null_int) then
                        Boxes2D(WILB:WILB+Halo_Points-1,WJLB:WJUB) = -99
                    endif

                    if (NeighbourNorth /= null_int) then
                        Boxes2D(WIUB-Halo_Points+2:WIUB,WJLB:WJUB) = -99
                    endif

                    if (NeighbourWest /= null_int) then
                        Boxes2D(WILB:WIUB, WJLB:WJLB+Halo_Points-1) = -99
                    endif

                    if (NeighbourEast /= null_int) then
                        Boxes2D(WILB:WIUB, WJUB-Halo_Points+2:WJUB) = -99
                    endif

                endif

                if (present(Boxes3D)) then

                    if (present(KLB).and.present(KUB)) then
                        if (NeighbourSouth /= null_int) then
                            Boxes3D(WILB:WILB+Halo_Points-1,WJLB:WJUB,KLB:KUB) = -99
                        endif

                        if (NeighbourNorth /= null_int) then
                            Boxes3D(WIUB-Halo_Points+2:WIUB,WJLB:WJUB,KLB:KUB) = -99
                        endif

                        if (NeighbourWest /= null_int) then
                            Boxes3D(WILB:WIUB, WJLB:WJLB+Halo_Points-1,KLB:KUB) = -99
                        endif

                        if (NeighbourEast /= null_int) then
                            Boxes3D(WILB:WIUB, WJUB-Halo_Points+2:WJUB,KLB:KUB) = -99
                        endif
                    else
                        stop "HorizontalGrid - No_BoxMap_HaloArea - ERR10"
                    endif

                endif

            endif

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


    end subroutine No_BoxMap_HaloArea


    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine GetHorizontalGridSize(HorizontalGridID, Size, WorkSize, GlobalWorkSize, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        type (T_Size2D), optional, intent(OUT)      :: Size
        type (T_Size2D), optional, intent(OUT)      :: WorkSize
        type (T_Size2D), optional, intent(OUT)      :: GlobalWorkSize
        integer,         optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !Size
            if (present(Size            )) Size           = Me%Size
            if (present(WorkSize        )) WorkSize       = Me%WorkSize
            if (present(GlobalWorkSize  )) GlobalWorkSize = Me%GlobalWorkSize

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetHorizontalGridSize

    !--------------------------------------------------------------------------

    subroutine GetHorizontalGrid(HorizontalGridID, XX_IE, YY_IE, XX_Z, YY_Z,            &
                                 XX_U, YY_U, XX_V, YY_V, XX_Cross, YY_Cross,            &
                                 DXX, DYY, DZX, DZY, DUX, DUY, DVX, DVY, XX, YY,        &
                                 XX2D_Z, YY2D_Z, XX2D_U, YY2D_U, XX2D_V, YY2D_V,        &
                                 IV, JV, IU, JU, IZ, JZ,                                &
                                 ILinkV, JLinkV, ILinkU, JLinkU, ILinkZ, JLinkZ,        &
                                 STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real, dimension(:, :), pointer, optional    :: XX_IE, YY_IE, XX2D_Z, YY2D_Z
        real, dimension(:, :), pointer, optional    :: XX2D_U, YY2D_U, XX2D_V, YY2D_V
        real, dimension(:   ), pointer, optional    :: XX_Z, YY_Z, XX_U, YY_U, XX_V, YY_V, XX_Cross, YY_Cross
        real, dimension(:, :), pointer, optional    :: DXX, DYY, DZX, DZY
        real, dimension(:, :), pointer, optional    :: DUX, DUY, DVX, DVY
        integer, dimension(:, :), pointer, optional :: ILinkV, JLinkV, ILinkU, JLinkU, ILinkZ, JLinkZ
        integer, dimension(:, :), pointer, optional :: IV, JV, IU, JU, IZ, JZ
        real, dimension(:   ), pointer, optional    :: XX, YY
        integer, optional,  intent(OUT)             :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !XX_IE, YY_IE
            if (present(XX_IE)) then
                XX_IE => Me%XX_IE
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(YY_IE)) then
                YY_IE => Me%YY_IE
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            !XX2D_Z, YY2D_Z
            if (present(XX2D_Z)) then
                XX2D_Z => Me%Compute%XX2D_Z
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(YY2D_Z)) then
                YY2D_Z => Me%Compute%YY2D_Z
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif


            !XX_Z, YY_Z
            if (present(XX_Z)) then
                XX_Z => Me%Compute%XX_Z
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(YY_Z)) then
                YY_Z => Me%Compute%YY_Z
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif


            !XX_U, YY_U
            if (present(XX_U)) then
                XX_U => Me%Compute%XX_U
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(YY_U)) then
                YY_U => Me%Compute%YY_U
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif


            !XX_V, YY_V
            if (present(XX_V)) then
                XX_V => Me%Compute%XX_V
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(YY_V)) then
                YY_V => Me%Compute%YY_V
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            !XX_Cross, YY_Cross
            if (present(XX_Cross)) then
                XX_Cross => Me%Compute%XX_Cross
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(YY_Cross)) then
                YY_Cross => Me%Compute%YY_Cross
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif


            !DXX, DYY
            if (present(DXX)) then
                DXX => Me%DXX
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(DYY)) then
                DYY => Me%DYY
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            !DZX, DZY
            if (present(DZX)) then
                DZX => Me%DZX
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(DZY)) then
                DZY => Me%DZY
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            !DUX, DUY, DVX, DVY
            if (present(DUX)) then
                DUX => Me%DUX
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(DUY)) then
                DUY => Me%DUY
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(DVX)) then
                DVX => Me%DVX
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(DVY)) then
                DVY => Me%DVY
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            !XX, YY
            if (present(XX)) then
                XX => Me%XX
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(YY)) then
                YY => Me%YY
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            !XX2D_U, YY2D_U
            if (present(XX2D_U)) then
                XX2D_U => Me%Compute%XX2D_U
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(YY2D_U)) then
                YY2D_U => Me%Compute%YY2D_U
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            !XX2D_V, YY2D_V
            if (present(XX2D_V)) then
                XX2D_V => Me%Compute%XX2D_V
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            if (present(YY2D_V)) then
                YY2D_V => Me%Compute%YY2D_V
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            !IV. Father cell inside which is each son cell (row).
            if (present(ILinkV)) then
                ILinkV => Me%LastFatherGrid%ILinkV
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            !JV. Father cell inside which is each son cell(column).
            if (present(JLinkV)) then
                JLinkV => Me%LastFatherGrid%JLinkV
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            !IV. Father cell inside which is each son cell (row).
            if (present(ILinkU)) then
                ILinkU => Me%LastFatherGrid%ILinkU
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            !JV. Father cell inside which is each son cell(column).
            if (present(JLinkU)) then
                JLinkU => Me%LastFatherGrid%JLinkU
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            !IV. Father cell inside which is each son cell (row).
            if (present(JLinkZ)) then
                JLinkZ => Me%LastFatherGrid%JLinkZ
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            !JV. Father cell inside which is each son cell(column).
            if (present(ILinkZ)) then
                ILinkZ => Me%LastFatherGrid%ILinkZ
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            !JV. Father cell inside which is each son cell(column).
            if (present(IZ)) then
                IZ => Me%LastFatherGrid%IZ
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            !JV. Father cell inside which is each son cell(column).
            if (present(JZ)) then
                JZ => Me%LastFatherGrid%JZ
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            !JV. Father cell inside which is each son cell(column).
            if (present(JU)) then
                JU => Me%LastFatherGrid%JU
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            !JV. Father cell inside which is each son cell(column).
            if (present(JV)) then
                JV => Me%LastFatherGrid%JV
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            !JV. Father cell inside which is each son cell(column).
            if (present(IU)) then
                IU => Me%LastFatherGrid%IU
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif
            !JV. Father cell inside which is each son cell(column).
            if (present(IV)) then
                IV => Me%LastFatherGrid%IV
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            endif

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) &
            STAT = STAT_

    end subroutine GetHorizontalGrid

    !--------------------------------------------------------------------------

    subroutine GetGridCellArea (HorizontalGridID, GridCellArea, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real, dimension(:, :), pointer, optional    :: GridCellArea
        integer, optional,  intent(OUT)             :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            GridCellArea => Me%GridCellArea
            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)


            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

    end subroutine GetGridCellArea

    !--------------------------------------------------------------------------

    subroutine GetGridOrigin(HorizontalGridID, Xorig, Yorig, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real,              intent(OUT)              :: Xorig, Yorig
        integer, optional, intent(OUT)              :: STAT


        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Xorig = Me%Xorig
            Yorig = Me%Yorig

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetGridOrigin

    !--------------------------------------------------------------------------

    subroutine GetGridAngle(HorizontalGridID, Angle, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real,              intent(OUT)              :: Angle
        integer, optional, intent(OUT)              :: STAT


        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Angle = Me%Grid_Angle

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetGridAngle

        !--------------------------------------------------------------------------

    subroutine GetLatitudeLongitude(HorizontalGridID, Latitude, Longitude, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real,    optional, intent(OUT)              :: Latitude
        real,    optional, intent(OUT)              :: Longitude
        integer, optional, intent(OUT)              :: STAT


        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if(present(Latitude)) Latitude  = Me%Latitude
            if(present(Longitude))Longitude = Me%Longitude

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetLatitudeLongitude

    !--------------------------------------------------------------------------

    subroutine GetGridZone(HorizontalGridID, Zone, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer,           intent(OUT)              :: Zone
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Zone = Me%ZoneLong

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetGridZone

    !--------------------------------------------------------------------------

    subroutine GetGridCoordType(HorizontalGridID, CoordType, ReadCartCorners, ProjType, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer,           intent(OUT)              :: CoordType
        logical, optional, intent(OUT)              :: ReadCartCorners
        integer, optional, intent(OUT)              :: ProjType
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            CoordType = Me%CoordType

            if (present(ReadCartCorners)) then
                ReadCartCorners = Me%ReadCartCorners
            endif

            if (present(ProjType)) then
                ProjType = Me%ProjType
            endif

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetGridCoordType

    !--------------------------------------------------------------------------

    subroutine GetCoordTypeList(GEOG, UTM, MIL_PORT, SIMPLE_GEOG, GRID_COORD, CIRCULAR, NLRD, &
                                LAMB_CONF_CONIC, PORTUGUESE_UTM_ZONE)

        !Arguments-------------------------------------------------------------

        integer, optional, intent(OUT) :: GEOG                  !Coordenadas Geograficas
        integer, optional, intent(OUT) :: SIMPLE_GEOG           !Coordenadas Geograficas Simplificadas
        integer, optional, intent(OUT) :: UTM                   !Coordenadas (UTM)
        integer, optional, intent(OUT) :: MIL_PORT              !Coordenadas Militares Portuguesas
        integer, optional, intent(OUT) :: GRID_COORD            !Coordenadas da malha
        integer, optional, intent(OUT) :: CIRCULAR              !Coordenadas circulares (XX - raio, YY - anglo em graus)
        integer, optional, intent(OUT) :: NLRD                  !Coordenadas Netherlands RD
        integer, optional, intent(OUT) :: LAMB_CONF_CONIC
        integer, optional, intent(OUT) :: PORTUGUESE_UTM_ZONE   !

        !----------------------------------------------------------------------

        if (present(GEOG                )) GEOG                 = GEOG_
        if (present(SIMPLE_GEOG         )) SIMPLE_GEOG          = SIMPLE_GEOG_
        if (present(UTM                 )) UTM                  = UTM_
        if (present(MIL_PORT            )) MIL_PORT             = MIL_PORT_
        if (present(GRID_COORD          )) GRID_COORD           = GRID_COORD_
        if (present(CIRCULAR            )) CIRCULAR             = CIRCULAR_
        if (present(NLRD                )) NLRD                 = NLRD_
        if (present(LAMB_CONF_CONIC     )) LAMB_CONF_CONIC      = LAMB_CONF_CONIC_

        if (present(PORTUGUESE_UTM_ZONE )) PORTUGUESE_UTM_ZONE  = PORTUGUESE_UTM_ZONE_

        !----------------------------------------------------------------------

    end subroutine GetCoordTypeList

    !--------------------------------------------------------------------------


    subroutine GetGridMeanLatLong(HorizontalGridID, Latitude, Longitude, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real,              intent(OUT)              :: Latitude, Longitude
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Latitude  = Me%Latitude
            Longitude = Me%Longitude

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetGridMeanLatLong

    !--------------------------------------------------------------------------

    subroutine GetCoriolisFrequency(HorizontalGridID, F, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real, dimension(:, :), pointer              :: F
        integer, optional,  intent(OUT)             :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            F => Me%F
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if

        if (present(STAT)) STAT = SUCCESS_

    end subroutine GetCoriolisFrequency

    !----------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetDefineCellsMap(HorizontalGridID, DefineCellsMap, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer, dimension(:, :), pointer           :: DefineCellsMap
        integer, optional,  intent(OUT)             :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            DefineCellsMap => Me%DefineCellsMap
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if

        if (present(STAT)) STAT = SUCCESS_

    end subroutine GetDefineCellsMap

    !----------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetNotDefinedCells(HorizontalGridID, NotDefinedCells, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        logical                                     :: NotDefinedCells
        integer, optional,  intent(OUT)             :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            NotDefinedCells = Me%NotDefinedCells
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if

        if (present(STAT)) STAT = SUCCESS_

    end subroutine GetNotDefinedCells

    !----------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine GetCheckDistortion(HorizontalGridID, Distortion, CornersXYInput, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        logical          ,  intent(OUT)             :: Distortion
        logical, optional,  intent(OUT)             :: CornersXYInput
        integer, optional,  intent(OUT)             :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(CornersXYInput)) CornersXYInput = Me%CornersXYInput

            Distortion = Me%Distortion

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if

        if (present(STAT)) STAT = SUCCESS_

    end subroutine GetCheckDistortion

    !----------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetGridRotation(HorizontalGridID, RotationX, RotationY, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real, dimension(:, :), pointer              :: RotationX, RotationY
        integer, optional,  intent(OUT)             :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            RotationX  => Me%RotationX

            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            RotationY  => Me%RotationY

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if

        if (present(STAT)) STAT = SUCCESS_

    end subroutine GetGridRotation

    !----------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    !Get cell rotation in radians
    subroutine GetCellRotation(HorizontalGridID, i, j, CellRotationX, CellRotationY, STAT)

        !Arguments-------------------------------------------------------------
        integer                                 :: HorizontalGridID
        integer                                 :: i, j
        real                                    :: CellRotationX
        real,    optional                       :: CellRotationY
        integer, optional                       :: STAT
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, ready_

        !Begin--------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then

            
            CellRotationX = 0.
            
            if      (Me%Distortion) then

                CellRotationX = Me%RotationX(i, j)
                
            else if (Me%RegularRotation) then

                CellRotationX = Me%Grid_Angle * Pi / 180.
                
            endif
            
            if (present(CellRotationY)) then
                
                CellRotationY = Pi/2.
            
                if      (Me%Distortion) then

                    CellRotationY = Me%RotationY(i, j)                

                else if (Me%RegularRotation) then

                    CellRotationY = CellRotationX + Pi/2.

                endif            
            endif

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1

        if (present(STAT))  STAT = STAT_

    end subroutine GetCellRotation

    !--------------------------------------------------------------------------

    subroutine GetGridLatitudeLongitude(HorizontalGridID,                       &
                                        GridLatitude,                           &
                                        GridLongitude,                          &
                                        GridLatitudeConn,                       &
                                        GridLongitudeConn,                      &
                                        STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real, dimension(:,:), optional, pointer     :: GridLongitude, GridLatitude
        real, dimension(:,:), optional, pointer     :: GridLongitudeConn, GridLatitudeConn
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ == IDLE_ERR_     ) .OR.                                     &
            (ready_ == READ_LOCK_ERR_)) then

cd1 :       if (present(GridLongitude)) then
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
                GridLongitude => Me%LongitudeZ
            end if cd1


cd2 :       if (present(GridLatitude)) then
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
                GridLatitude => Me%LatitudeZ
            end if cd2

cd3 :       if (present(GridLatitudeConn)) then
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
                GridLatitudeConn => Me%LatitudeConn
            end if cd3


cd4 :       if (present(GridLongitudeConn)) then
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
                GridLongitudeConn => Me%LongitudeConn
            end if cd4


            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT))                                                      &
            STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine GetGridLatitudeLongitude

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------


    subroutine GetZCoordinates(HorizontalGridID, CoordX, CoordY, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real, dimension(:,:),  pointer              :: CoordX, CoordY
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ == IDLE_ERR_     ) .OR.                                     &
            (ready_ == READ_LOCK_ERR_)) then

            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)

            if (Me%CoordType == GEOG_ .or. Me%CoordType == SIMPLE_GEOG_) then
                CoordX => Me%LongitudeZ
                CoordY => Me%LatitudeZ
            else
                CoordX => Me%Compute%XX2D_Z
                CoordY => Me%Compute%YY2D_Z
            endif

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT))                                                      &
            STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine GetZCoordinates

    subroutine GetCornersCoordinates(HorizontalGridID, CoordX, CoordY, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real, dimension(:,:), optional, pointer     :: CoordX, CoordY
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ == IDLE_ERR_     ) .OR.                                     &
            (ready_ == READ_LOCK_ERR_)) then

            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)

            if (Me%CoordType == GEOG_ .or. Me%CoordType == SIMPLE_GEOG_) then
                CoordX => Me%LongitudeConn
                CoordY => Me%LatitudeConn
            else
                CoordX => Me%XX_IE
                CoordY => Me%YY_IE
            endif

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT))                                                      &
            STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine GetCornersCoordinates

    !----------------------------------------------------------------------------

    subroutine GetComputeZUV(HorizontalGridID, ComputeZ, ComputeU, ComputeV,             &
                             ComputeCross,  ComputeZU, ComputeZV, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer, optional,  intent(OUT)             :: ComputeZ, ComputeU, ComputeV
        integer, optional,  intent(OUT)             :: ComputeCross, ComputeZU, ComputeZV
        integer, optional,  intent(OUT)  :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ == IDLE_ERR_     ) .OR.                                              &
            (ready_ == READ_LOCK_ERR_)) then

cd1 :       if (present(ComputeZ)) then
                ComputeZ = ComputeZ_
            end if cd1

cd2 :       if (present(ComputeU)) then
                ComputeU = ComputeU_
            end if cd2

cd3 :       if (present(ComputeV)) then
                ComputeV = ComputeV_
            end if cd3

cd4 :       if (present(ComputeCross)) then
                ComputeCross = ComputeCross_
            end if cd4

cd5 :       if (present(ComputeZU)) then
                ComputeZU = ComputeZU_
            end if cd5

cd6 :       if (present(ComputeZV)) then
                ComputeZV = ComputeZV_
            end if cd6

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT))                                                      &
            STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine GetComputeZUV

    !----------------------------------------------------------------------------

    subroutine GetFatherGridID(HorizontalGridID, LargeScaleModel, Assimila, SurfaceMM5, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer, optional,  intent(OUT)             :: LargeScaleModel, Assimila, SurfaceMM5
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ == IDLE_ERR_     ) .OR.                                              &
            (ready_ == READ_LOCK_ERR_)) then

cd1 :       if (present(LargeScaleModel)) then
                LargeScaleModel = LargeScaleModel_
            end if cd1

cd2 :       if (present(Assimila)) then
                Assimila = Assimila_
            end if cd2

cd3 :       if (present(SurfaceMM5)) then
                SurfaceMM5 = SurfaceMM5_
            end if cd3

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT))                                                      &
            STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine GetFatherGridID

    !--------------------------------------------------------------------------

    subroutine GetGridFileName (HorizontalGridID, FileName, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: HorizontalGridID
        character(len=*)                            :: FileName
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ == IDLE_ERR_     ) .OR.                                              &
            (ready_ == READ_LOCK_ERR_)) then

            FileName = Me%FileName
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetGridFileName
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Maps the IJ coordinates of a cell to a unique ID
    !---------------------------------------------------------------------------
    subroutine GetCellIDfromIJ(HorizontalGridID, I, J, ID)
        !Arguments---------------------------------------------------------------
        integer,            intent(IN)              :: HorizontalGridID
        integer,            intent(in)              :: I, J
        integer,            intent(OUT)             :: ID
        !Local-------------------------------------------------------------------
        integer                                     :: ready_
        !------------------------------------------------------------------------
        call Ready(HorizontalGridID, ready_)
        if ((ready_ == IDLE_ERR_) .or. (ready_ == READ_LOCK_ERR_)) then
            ID = (I-1)*(Me%Size%JUB - Me%Size%JLB) + J
        end if

    end subroutine GetCellIDfromIJ
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the IJ coordinates of a cell given a unique ID
    !---------------------------------------------------------------------------
    subroutine GetCellIJfromID(HorizontalGridID, I, J, ID)
        !Arguments---------------------------------------------------------------
        integer,            intent(IN)              :: HorizontalGridID
        integer,            intent(OUT)             :: I, J
        integer,            intent(in)              :: ID
        !Local-------------------------------------------------------------------
        integer                                     :: ready_
        !------------------------------------------------------------------------
        call Ready(HorizontalGridID, ready_)
        if ((ready_ == IDLE_ERR_) .or. (ready_ == READ_LOCK_ERR_)) then
            J = mod(ID, (Me%Size%JUB - Me%Size%JLB))
            I = (ID-J)/(Me%Size%JUB - Me%Size%JLB) + 1            
        end if

    end subroutine GetCellIJfromID
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> Gets the ID coordinates array of an IJ array
    !---------------------------------------------------------------------------
    subroutine GetCellIDfromIJArray(HorizontalGridID, ArrayIJ, ArrayID)
        !Arguments---------------------------------------------------------------
        integer,            intent(IN)              :: HorizontalGridID
        integer, dimension(:,:), intent(in)         :: ArrayIJ
        integer, dimension(:), intent(inout)        :: ArrayID
        !Local-------------------------------------------------------------------
        integer                                     :: ready_, i
        !------------------------------------------------------------------------
        call Ready(HorizontalGridID, ready_)
        if ((ready_ == IDLE_ERR_) .or. (ready_ == READ_LOCK_ERR_)) then
            do i=1, size(ArrayID)
                call GetCellIDfromIJ(HorizontalGridID, ArrayIJ(i,1), ArrayIJ(i,2), ArrayID(i))
            end do
        end if

    end subroutine GetCellIDfromIJArray

    !--------------------------------------------------------------------------

    subroutine GetXYCellZ(HorizontalGridID, XPoint, YPoint, I, J, PercI, PercJ,         &
        Referential, Iold, Jold, STAT)
    !Arguments---------------------------------------------------------------
    integer,            intent(IN)              :: HorizontalGridID
    real,               intent(IN)              :: XPoint, YPoint
    integer,            intent(OUT)             :: I, J
    real,    optional,  intent(OUT)             :: PercI, PercJ
    integer, optional,  intent(IN)              :: Referential
    integer, optional,  intent(IN)              :: Iold, Jold
    integer, optional,  intent(OUT)             :: STAT
    !Local-------------------------------------------------------------------
    real,   dimension(:,:), pointer             :: XX2D, YY2D
    real,   dimension(:  ), pointer             :: XX1D, YY1D
    real,   dimension(:  ), pointer             :: XX1D_Aux, YY1D_Aux
    integer                                     :: STAT_
    integer                                     :: ready_
    real                                        :: XPoint2, YPoint2, Xorig2, Yorig2
    integer                                     :: Referential_, GetGridBorderType
    integer                                     :: ILB, IUB, JLB, JUB
    logical                                     :: CellLocated
    integer                                     :: Iold_, Jold_
    !------------------------------------------------------------------------

    STAT_ = UNKNOWN_
    call Ready(HorizontalGridID, ready_)

    if ((ready_ == IDLE_ERR_) .or. (ready_ == READ_LOCK_ERR_)) then
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        XX2D        => Me%XX_IE
        YY2D        => Me%YY_IE
        XX1D        => Me%Compute%XX_Cross
        YY1D        => Me%Compute%YY_Cross
        allocate(XX1D_Aux(JLB:JUB))
        allocate(YY1D_Aux(ILB:IUB))

        if (present(Referential)) then
            Referential_ = Referential
        else
            Referential_ = GridCoord_ !default  behaviour
        end if
        if (Referential_ == Cartesian_) GetGridBorderType = Me%GridBorderCart%Type_
        if (Referential_ == GridCoord_) then
            GetGridBorderType = Me%GridBorderCoord%Type_
            if (Me%CoordType == SIMPLE_GEOG_ .or. Me%CoordType == GEOG_) then
                XX2D        => Me%LongitudeConn
                YY2D        => Me%LatitudeConn
            end if
        end if
        if (Referential_ == AlongGrid_) then
            GetGridBorderType = Me%GridBorderAlongGrid%Type_
            XX2D        => Me%XX_AlongGrid
            YY2D        => Me%YY_AlongGrid
        end if

        if (GetGridBorderType == ComplexPolygon_) then
            if (present(Iold) .and. present(Jold)) then
                Iold_ = Iold
                Jold_ = Jold
            else
                Iold_ = FillValueInt
                Jold_ = FillValueInt
            end if

            call LocateCellPolygons(XX2D,                                           &
                YY2D,                                           &
                XPoint, YPoint, Me%DefineCellsMap,              &
                !                                        Me%WorkSize%ILB, Me%WorkSize%IUB + 1,           &
                !                                        Me%WorkSize%JLB, Me%WorkSize%JUB + 1,           &
                Me%WorkSize%ILB, Me%WorkSize%IUB,           &
                Me%WorkSize%JLB, Me%WorkSize%JUB,           &
                I, J, CellLocated, Iold_, Jold_)

            if (I < 0 .or. J < 0  .or. .not. CellLocated) then
                STAT_ = OUT_OF_BOUNDS_ERR_
                !stop 'GetXYCellZ - ModuleHorizontalGrid - ERR10'
            else
                if (present(PercI) .and. present(PercJ)) then
                    call RelativePosition4VertPolygon(Xa = XX2D(I+1, J  ), Ya = YY2D(I+1, J  ), &
                        Xb = XX2D(I+1, J+1), Yb = YY2D(I+1, J+1), &
                        Xc = XX2D(I  , J  ), Yc = YY2D(I  , J  ), &
                        Xd = XX2D(I  , J+1), Yd = YY2D(I  , J+1), &
                        Xe = XPoint,         Ye = YPoint,         &
                        Xex= PercJ,          Yey= PercI)
                end if
            end if
        else
            XPoint2 = XPoint
            YPoint2 = YPoint

            if (GetGridBorderType == Rectang_ .or. Referential_ == AlongGrid_) then
                XX1D_Aux(JLB:JUB) = XX2D(ILB+1  ,JLB:JUB)
                YY1D_Aux(ILB:IUB) = YY2D(ILB:IUB,JLB+1  )
            else
                Xorig2  = Me%Xorig
                Yorig2  = Me%Yorig
                XPoint2 =  XPoint - Xorig2
                YPoint2 =  YPoint - Yorig2
                call RODAXY(0., 0., -Me%Grid_Angle, XPoint2, YPoint2)
                XX1D_Aux(JLB+1:JUB) = Me%XX(JLB+1:JUB)
                YY1D_Aux(ILB+1:IUB) = Me%YY(ILB+1:IUB)
            end if
            
            call LocateCell (XX1D_Aux,                                           &
                YY1D_Aux,                                           &
                XPoint2, YPoint2,                                      &
                Me%WorkSize%ILB, Me%WorkSize%IUB + 1,                  &
                Me%WorkSize%JLB, Me%WorkSize%JUB + 1,                  &
                I, J)

            if (present(PercI)) then
                if (I < 0) then
                    STAT_ = OUT_OF_BOUNDS_ERR_
                else
                    PercI  = (YPoint2 - YY1D_Aux(I)) / (YY1D_Aux(I+1) - YY1D_Aux(I))
                end if
            end if

            if (present(PercJ)) then
                if (J < 0) then
                    STAT_ = OUT_OF_BOUNDS_ERR_
                else
                    PercJ  = (XPoint2 - XX1D_Aux(J)) / (XX1D_Aux(J+1) - XX1D_Aux(J))
                end if
            end if
        end if

        deallocate(XX1D_Aux)
        deallocate(YY1D_Aux)
        nullify(XX1D_Aux)
        nullify(YY1D_Aux)
        nullify(XX2D)
        nullify(YY2D)

        if (STAT_ == UNKNOWN_) STAT_ = SUCCESS_

    else
        STAT_ = ready_
    end if

    if (present(STAT)) STAT = STAT_

    end subroutine GetXYCellZ

    !--------------------------------------------------------------------------

    subroutine GetXYCellZ_ThreadSafe(HorizontalGridID, XPoint, YPoint, I, J, PercI, PercJ, Referential, STAT)

        !Arguments---------------------------------------------------------------
        integer,            intent(IN)              :: HorizontalGridID
        real,               intent(IN)              :: XPoint, YPoint
        integer,            intent(OUT)             :: I, J
        real,    optional,  intent(OUT)             :: PercI, PercJ
        integer, optional,  intent(IN)              :: Referential
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        real,   dimension(:,:), pointer             :: XX2D, YY2D
        real,   dimension(:  ), pointer             :: XX1D, YY1D
        real,   dimension(:  ), pointer             :: XX1D_Aux, YY1D_Aux
        integer                                     :: STAT_
        integer                                     :: ready_
        real                                        :: XPoint2, YPoint2, Xorig2, Yorig2
        integer                                     :: Referential_, GetGridBorderType
        integer                                     :: ILB, IUB, JLB, JUB
        type (T_HorizontalGrid), pointer            :: LocalMe
        logical                                     :: CellLocated
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        LocalMe => Ready_ThreadSafe(HorizontalGridID, ready_)

i1:     if ((ready_ == IDLE_ERR_     ) .OR.                                             &
            (ready_ == READ_LOCK_ERR_)) then

            ILB = LocalMe%Size%ILB
            IUB = LocalMe%Size%IUB

            JLB = LocalMe%Size%JLB
            JUB = LocalMe%Size%JUB

            if (present(Referential)) then

                Referential_ = Referential

            else

                Referential_ = GridCoord_

            endif


            XX2D        => LocalMe%XX_IE
            YY2D        => LocalMe%YY_IE

            XX1D        => LocalMe%Compute%XX_Cross
            YY1D        => LocalMe%Compute%YY_Cross

            allocate(XX1D_Aux(JLB:JUB))
            allocate(YY1D_Aux(ILB:IUB))

            if (Referential_ == Cartesian_) then

                GetGridBorderType = LocalMe%GridBorderCart%Type_

            endif

            if (Referential_ == GridCoord_) then

                GetGridBorderType = LocalMe%GridBorderCoord%Type_

                if (LocalMe%CoordType == SIMPLE_GEOG_ .or. LocalMe%CoordType == GEOG_) then

                    XX2D        => LocalMe%LongitudeConn
                    YY2D        => LocalMe%LatitudeConn

                endif

            endif

            if (Referential_ == AlongGrid_) then

                GetGridBorderType = LocalMe%GridBorderAlongGrid%Type_

                XX2D        => LocalMe%XX_AlongGrid
                YY2D        => LocalMe%YY_AlongGrid

            endif


i2:         if (GetGridBorderType == ComplexPolygon_) then

                call LocateCellPolygons(XX2D,                                           &
                                        YY2D,                                           &
                                        XPoint, YPoint, LocalMe%DefineCellsMap,              &
!                                        LocalMe%WorkSize%ILB, LocalMe%WorkSize%IUB + 1,           &
!                                        LocalMe%WorkSize%JLB, LocalMe%WorkSize%JUB + 1,           &
                                        LocalMe%WorkSize%ILB, LocalMe%WorkSize%IUB,           &
                                        LocalMe%WorkSize%JLB, LocalMe%WorkSize%JUB,           &
                                        I, J, CellLocated)

                if (I < 0 .or. J < 0 .or. .not. CellLocated) then
                    STAT_ = OUT_OF_BOUNDS_ERR_
                    !stop 'GetXYCellZ - ModuleHorizontalGrid - ERR10'

                else

                    if (present(PercI) .and. present(PercJ)) then
                        !
                        call RelativePosition4VertPolygon(Xa = XX2D(I+1, J  ), Ya = YY2D(I+1, J  ), &
                                                          Xb = XX2D(I+1, J+1), Yb = YY2D(I+1, J+1), &
                                                          Xc = XX2D(I  , J  ), Yc = YY2D(I  , J  ), &
                                                          Xd = XX2D(I  , J+1), Yd = YY2D(I  , J+1), &
                                                          Xe = XPoint,         Ye = YPoint,         &
                                                          Xex= PercJ,          Yey= PercI)
                    endif

                endif

            else i2

                XPoint2 = XPoint
                YPoint2 = YPoint

                if (GetGridBorderType == Rectang_ .or. Referential_ == AlongGrid_) then

                    XX1D_Aux(JLB:JUB) = XX2D(ILB+1  ,JLB:JUB)
                    YY1D_Aux(ILB:IUB) = YY2D(ILB:IUB,JLB+1  )

                else

                    Xorig2  = LocalMe%Xorig
                    Yorig2  = LocalMe%Yorig

                    XPoint2 =  XPoint - Xorig2
                    YPoint2 =  YPoint - Yorig2

                    call RODAXY(0., 0., -LocalMe%Grid_Angle, XPoint2, YPoint2)

                    XX1D_Aux(JLB+1:JUB) = LocalMe%XX(JLB+1:JUB)
                    YY1D_Aux(ILB+1:IUB) = LocalMe%YY(ILB+1:IUB)

                endif

                call LocateCell (XX1D_Aux,                                           &
                                 YY1D_Aux,                                           &
                                 XPoint2, YPoint2,                                      &
                                 LocalMe%WorkSize%ILB, LocalMe%WorkSize%IUB + 1,                  &
                                 LocalMe%WorkSize%JLB, LocalMe%WorkSize%JUB + 1,                  &
                                 I, J)

                if (present(PercI)) then
                    if (I < 0) then
                        STAT_ = OUT_OF_BOUNDS_ERR_
                    else
                        PercI  = (YPoint2 - YY1D_Aux(I)) / (YY1D_Aux(I+1) - YY1D_Aux(I))
                    endif
                endif

                if (present(PercJ)) then
                    if (J < 0) then
                        STAT_ = OUT_OF_BOUNDS_ERR_
                    else
                        PercJ  = (XPoint2 - XX1D_Aux(J)) / (XX1D_Aux(J+1) - XX1D_Aux(J))
                    endif
                endif

            endif i2

            deallocate(XX1D_Aux)
            deallocate(YY1D_Aux)

            nullify(XX1D_Aux)
            nullify(YY1D_Aux)

            nullify(XX2D)
            nullify(YY2D)

            if (STAT_ == UNKNOWN_) STAT_ = SUCCESS_

        else    i1

            STAT_ = ready_

        end if  i1

        if (present(STAT)) STAT = STAT_

    end subroutine GetXYCellZ_ThreadSafe
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - Bentley Systems
    !> @brief
    !> mapps and array of XY point corrdinates to an array of IJ cell coordinates
    !---------------------------------------------------------------------------
    subroutine GetXYArrayIJ(HorizontalGridID, mapArrayXY, mapArrayIJ)
    integer, intent(inout) :: HorizontalGridID
    real, dimension(:,:), intent(inout) :: mapArrayXY
    integer, dimension(:,:), intent(inout) :: mapArrayIJ
    integer :: i
    
    mapArrayIJ = null_int
    do i=1, size(mapArrayXY, 1)
        if (GetXYInsideDomain(HorizontalGridID, mapArrayXY(i,1), mapArrayXY(i,2))) then
            call GetXYCellZ(HorizontalGridID, mapArrayXY(i,1), mapArrayXY(i,2),         &
                                              mapArrayIJ(i,1), mapArrayIJ(i,2))
        endif
        !print*, 'id=',i, 'x=',mapArrayXY(i,1), 'y=',mapArrayXY(i,2)
        !print*, 'I=',mapArrayIJ(i,1), 'J=',mapArrayIJ(i,2)
    end do
    
    end subroutine GetXYArrayIJ

    !--------------------------------------------------------------------------


    subroutine GetCellZ_XY(HorizontalGridID, I, J, PercI, PercJ, XPoint, YPoint, &
                           Xin, Yin, Referential, STAT)

        !Arguments---------------------------------------------------------------
        integer,            intent(OUT)             :: HorizontalGridID
        integer,            intent(IN )             :: I, J
        real,               intent(IN )             :: PercI, PercJ
        real,               intent(OUT)             :: XPoint, YPoint
        integer, optional,  intent(IN )             :: Referential
        real,    optional,  intent(IN )             :: Xin, Yin
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        real,   dimension(:,:), pointer             :: XX2D, YY2D
        integer                                     :: STAT_
        integer                                     :: ready_
        real                                        :: xac, yac, xbd, ybd
        integer                                     :: Referential_
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

i1:     if ((ready_ == IDLE_ERR_     ) .OR.                                             &
            (ready_ == READ_LOCK_ERR_)) then

            if (present(Referential)) then

                Referential_ = Referential

            else

                Referential_ = GridCoord_

            endif

            if (Me%CoordType == SIMPLE_GEOG_ .and. present(Xin) .and. present(Yin)) then 
                
                if (Referential_ == GridCoord_) then
                    
                    call FromCartToSpherical(Lat = YPoint,                              &
                                             Lon = XPoint,                              &
                                             X   = Xin,                                 &
                                             Y   = Yin)    
                else
            
                    call FromSphericalToCart(Lat = Yin,                                 &
                                             Lon = Xin,                                 &
                                             X   = XPoint,                              &
                                             Y   = YPoint)    
                endif
                
            else    
            XX2D        => Me%XX_IE
            YY2D        => Me%YY_IE

                if (Referential_ == GridCoord_.and. (Me%CoordType == GEOG_ .or. Me%CoordType == SIMPLE_GEOG_)) then

                XX2D        => Me%LongitudeConn
                YY2D        => Me%LatitudeConn

            endif

            if (Referential_ == AlongGrid_) then

                XX2D        => Me%XX_AlongGrid
                YY2D        => Me%YY_AlongGrid

            endif

            xac    = XX2D(I+1, J) * PercI + XX2D(I, J) * (1. - PercI)
            yac    = YY2D(I+1, J) * PercI + YY2D(I, J) * (1. - PercI)

            xbd    = XX2D(I+1, J+1) * PercI + XX2D(I, J+1) * (1. - PercI)
            ybd    = YY2D(I+1, J+1) * PercI + YY2D(I, J+1) * (1. - PercI)


            XPoint = xbd * PercJ + xac * (1. - PercJ)
            YPoint = ybd * PercJ + yac * (1. - PercJ)

            nullify(XX2D)
            nullify(YY2D)
            endif                       

            STAT_ = SUCCESS_
        else    i1
            STAT_ = ready_
        end if  i1

        if (present(STAT)) STAT = STAT_

    end subroutine GetCellZ_XY

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    logical function GetXYInsideDomain(HorizontalGridID, XPoint, YPoint, Referential, STAT)

        !Arguments---------------------------------------------------------------
        integer,            intent(IN)              :: HorizontalGridID
        real,               intent(IN)              :: XPoint, YPoint
        integer, optional,  intent(IN)              :: Referential
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        real,   dimension(:,:), pointer             :: XX2D, YY2D
        integer                                     :: Referential_
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: STAT_
        integer                                     :: ready_
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

i1:     if ((ready_ == IDLE_ERR_     ) .OR.                                             &
            (ready_ == READ_LOCK_ERR_)) then

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB
            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            if (present(Referential)) then

                Referential_ = Referential

            else

                Referential_ = GridCoord_

            endif

iR:         if (Referential_ == GridCoord_) then

                if (Me%CoordType == SIMPLE_GEOG_ .or. Me%CoordType == GEOG_) then

                    XX2D        => Me%LongitudeConn
                    YY2D        => Me%LatitudeConn

                else

                    XX2D        => Me%XX_IE
                    YY2D        => Me%YY_IE

                endif

                if (Me%GridBorderCoord%Type_ == Rectang_) then
                    if (XPoint >= XX2D(ILB+1, JLB+1) .and. XPoint <= XX2D(IUB, JUB) .and. &
                        YPoint >= YY2D(ILB+1, JLB+1) .and. YPoint <= YY2D(IUB, JUB)) then
                        GetXYInsideDomain = .true.
                    else
                        GetXYInsideDomain = .false.
                    endif

                else

                    GetXYInsideDomain = InsideDomainPolygon(Me%GridBorderCoord%Polygon_, XPoint, YPoint)
                endif


            else if (Referential_ == Cartesian_) then  iR

                if (Me%GridBorderCart%Type_ == Rectang_) then
                    if (XPoint >= Me%XX_IE(ILB+1, JLB+1) .and. XPoint <= Me%XX_IE(IUB, JUB) .and. &
                        YPoint >= Me%YY_IE(ILB+1, JLB+1) .and. YPoint <= Me%YY_IE(IUB, JUB)) then
                        GetXYInsideDomain = .true.
                    else
                        GetXYInsideDomain = .false.
                    endif

                else

                    GetXYInsideDomain = InsideDomainPolygon(Me%GridBorderCart%Polygon_, XPoint, YPoint)
                endif

            else if (Referential_ == AlongGrid_) then  iR

                if (Me%GridBorderAlongGrid%Type_ == Rectang_) then
                    if (XPoint >= Me%XX_AlongGrid(ILB+1, JLB+1) .and. XPoint <= Me%XX_AlongGrid(IUB, JUB) .and. &
                        YPoint >= Me%YY_AlongGrid(ILB+1, JLB+1) .and. YPoint <= Me%YY_AlongGrid(IUB, JUB)) then
                        GetXYInsideDomain = .true.
                    else
                        GetXYInsideDomain = .false.
                    endif

                else

                    GetXYInsideDomain = InsideDomainPolygon(Me%GridBorderAlongGrid%Polygon_, XPoint, YPoint)
                endif

            endif iR

            STAT_ = SUCCESS_
        else    i1
            STAT_ = ready_
        end if  i1

        if (present(STAT)) STAT = STAT_

    end function GetXYInsideDomain

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    integer function GetGridBorderType(HorizontalGridID, Referential, STAT)

        !Arguments---------------------------------------------------------------
        integer,            intent(IN)              :: HorizontalGridID
        integer,            intent(IN)              :: Referential
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

i1:     if ((ready_ == IDLE_ERR_     ) .OR.                                             &
            (ready_ == READ_LOCK_ERR_)) then

            if      (Referential == GridCoord_) then

                GetGridBorderType = Me%GridBorderCoord%Type_

            else if (Referential == Cartesian_) then

                GetGridBorderType = Me%GridBorderCart%Type_

            else if (Referential == AlongGrid_) then

                GetGridBorderType = Me%GridBorderAlongGrid%Type_

            endif

            STAT_ = SUCCESS_
        else    i1
            STAT_ = ready_
        end if  i1

        if (present(STAT)) STAT = STAT_

    end function GetGridBorderType

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetGridBorderPolygon(HorizontalGridID, Polygon, STAT)

        !Arguments---------------------------------------------------------------
        integer,                  intent(IN)        :: HorizontalGridID
        type(T_Polygon),  pointer                   :: Polygon
        integer, optional,        intent(OUT)       :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

i1:     if ((ready_ == IDLE_ERR_     ) .OR.                                             &
            (ready_ == READ_LOCK_ERR_)) then

            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)

            Polygon => Me%GridBorderCoord%Polygon_

            STAT_ = SUCCESS_
        else    i1
            STAT_ = ready_
        end if  i1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGridBorderPolygon

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetGridOutBorderPolygon(HorizontalGridID, Polygon, STAT)

        !Arguments---------------------------------------------------------------
        integer,                  intent(IN)        :: HorizontalGridID
        type(T_Polygon),  pointer                   :: Polygon
        integer, optional,        intent(OUT)       :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

i1:     if ((ready_ == IDLE_ERR_     ) .OR.                                             &
            (ready_ == READ_LOCK_ERR_)) then

            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)

            Polygon => Me%GridOutBorderCoord%Polygon_

            STAT_ = SUCCESS_
        else    i1
            STAT_ = ready_
        end if  i1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGridOutBorderPolygon

    !-------------------------------------------------------------------------

    subroutine GetGridBorderCartPolygon(HorizontalGridID, Polygon, STAT)

        !Arguments---------------------------------------------------------------
        integer,                  intent(IN)        :: HorizontalGridID
        type(T_Polygon),  pointer                   :: Polygon
        integer, optional,        intent(OUT)       :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

i1:     if ((ready_ == IDLE_ERR_     ) .OR.                                             &
            (ready_ == READ_LOCK_ERR_)) then

            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)

            Polygon => Me%GridBorderCart%Polygon_

            STAT_ = SUCCESS_
        else    i1
            STAT_ = ready_
        end if  i1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGridBorderCartPolygon

    !-------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine GetGridOutBorderCartPolygon(HorizontalGridID, Polygon, STAT)

        !Arguments---------------------------------------------------------------
        integer,                  intent(IN)        :: HorizontalGridID
        type(T_Polygon),  pointer                   :: Polygon
        integer, optional,        intent(OUT)       :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

i1:     if ((ready_ == IDLE_ERR_     ) .OR.                                             &
            (ready_ == READ_LOCK_ERR_)) then

            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)

            Polygon => Me%GridOutBorderCart%Polygon_

            STAT_ = SUCCESS_
        else    i1
            STAT_ = ready_
        end if  i1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGridOutBorderCartPolygon

    !-------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine GetGridBorderLimits(HorizontalGridID, West, East, South, North, STAT)

        !Arguments---------------------------------------------------------------
        integer,                  intent(IN)        :: HorizontalGridID
        real,    optional,        intent(OUT)       :: West, East, South, North
        integer, optional,        intent(OUT)       :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

i1:     if ((ready_ == IDLE_ERR_     ) .OR.                                             &
            (ready_ == READ_LOCK_ERR_)) then
    
            if (Me%BorderLimits%ON) then
                
                West    = Me%BorderLimits%Values(1)
                East    = Me%BorderLimits%Values(2)
                South   = Me%BorderLimits%Values(3)
                North   = Me%BorderLimits%Values(4)
                
            else                

                West    = Me%GridBorderCoord%Polygon_%Limits%Left
                East    = Me%GridBorderCoord%Polygon_%Limits%Right
                South   = Me%GridBorderCoord%Polygon_%Limits%Bottom
                North   = Me%GridBorderCoord%Polygon_%Limits%Top
                
            endif                

            STAT_ = SUCCESS_
        else    i1
            STAT_ = ready_
        end if  i1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGridBorderLimits

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetGridOutBorderCartLimits(HorizontalGridID, West, East, South, North, STAT)

        !Arguments---------------------------------------------------------------
        integer,                  intent(IN)        :: HorizontalGridID
        real,    optional,        intent(OUT)       :: West, East, South, North
        integer, optional,        intent(OUT)       :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

i1:     if ((ready_ == IDLE_ERR_     ) .OR.                                             &
            (ready_ == READ_LOCK_ERR_)) then

            West    = Me%GridOutBorderCart%Polygon_%Limits%Left
            East    = Me%GridOutBorderCart%Polygon_%Limits%Right
            South   = Me%GridOutBorderCart%Polygon_%Limits%Bottom
            North   = Me%GridOutBorderCart%Polygon_%Limits%Top

            STAT_ = SUCCESS_
        else    i1
            STAT_ = ready_
        end if  i1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGridOutBorderCartLimits

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine GetCellZInterceptByXYZPoints(HorizontalGridID, XYZPoints,                &
                                            VectorI, VectorJ, VectorK, nCell, STAT)

        !Arguments---------------------------------------------------------------
        integer,            intent(IN)              :: HorizontalGridID
        type (T_XYZPoints),      pointer            :: XYZPoints
        integer, dimension(:),   pointer            :: VectorI, VectorJ, VectorK
        integer                                     :: nCell
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        type (T_XYZPoints), pointer                 :: CurrentXYZPoints
        integer, dimension(:), pointer              :: AuxI, AuxJ
        integer                                     :: nCell_, l
        integer                                     :: STAT_, STAT_CALL
        integer                                     :: ready_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

i1:     if ((ready_ == IDLE_ERR_     ) .OR.                                             &
            (ready_ == READ_LOCK_ERR_)) then

            nCell_ = (Me%WorkSize%JUB-Me%WorkSize%JLB) * (Me%WorkSize%IUB-Me%WorkSize%ILB)

            allocate(AuxI(nCell_)); allocate(AuxJ(nCell_));

            nCell_ = 1

            CurrentXYZPoints => XYZPoints

dw:         do while (associated(CurrentXYZPoints))

                do l = 1, CurrentXYZPoints%Count

                    call GetXYCellZ(Me%InstanceID, CurrentXYZPoints%X(l),               &
                                    CurrentXYZPoints%Y(l), AuxI(nCell_), AuxJ(nCell_),  &
                                    STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= OUT_OF_BOUNDS_ERR_) then
                        stop 'GetCellZInterceptByXYZPoints - ModuleHorizontalGrid - ERR10'
                    endif

                    if (STAT_CALL == OUT_OF_BOUNDS_ERR_) then
                        cycle
                    else
                        nCell_ = nCell_ + 1
                    endif

                enddo

                CurrentXYZPoints => CurrentXYZPoints%Next

            enddo dw

            nullify   (CurrentXYZPoints)

            nCell = nCell_ - 1

            allocate(VectorI(nCell), VectorJ(nCell), VectorK(nCell))

            VectorI(:) = AuxI(1:nCell)
            VectorJ(:) = AuxJ(1:nCell)
            VectorK(:) = FillValueInt

            deallocate(AuxI); deallocate(AuxJ);

            STAT_ = SUCCESS_
        else    i1
            STAT_ = ready_
        end if  i1

        if (present(STAT)) STAT = STAT_

    end subroutine GetCellZInterceptByXYZPoints

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetCellZInterceptByLine(HorizontalGridID, Line, WaterPoints2D, VectorI, VectorJ, VectorK, nCell, STAT)

        !Arguments---------------------------------------------------------------
        integer,            intent(IN)              :: HorizontalGridID
        type (T_Lines), pointer                     :: Line
        integer, dimension(:,:), pointer            :: WaterPoints2D
        integer, dimension(:),   pointer            :: VectorI, VectorJ, VectorK
        integer                                     :: nCell
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        real,    dimension(:,:), pointer            :: XX2D, YY2D
        type (T_Lines), pointer                     :: CurrentLine
        type (T_Segment), pointer                   :: Segment
        type (T_Polygon), pointer                   :: Polygon
        integer, dimension(:), pointer              :: AuxI, AuxJ
        integer                                     :: nCell_, i, j, l, k
        integer                                     :: STAT_
        integer                                     :: ready_
        logical                                     :: DoesCount

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

i1:     if ((ready_ == IDLE_ERR_     ) .OR.                                             &
            (ready_ == READ_LOCK_ERR_)) then

            if (Me%CoordType == SIMPLE_GEOG_ .or. Me%CoordType == GEOG_) then

                XX2D        => Me%LongitudeConn
                YY2D        => Me%LatitudeConn

            else

                XX2D        => Me%XX_IE
                YY2D        => Me%YY_IE

            endif

            allocate(Polygon)
            allocate(Segment)
            allocate(Segment%StartAt)
            allocate(Segment%EndAt)

            Polygon%Count = 5

            allocate(Polygon%VerticesF(Polygon%Count))

            nCell_ = (Me%WorkSize%JUB-Me%WorkSize%JLB) * (Me%WorkSize%IUB-Me%WorkSize%ILB)

            allocate(AuxI(nCell_)); allocate(AuxJ(nCell_));

            AuxI(:) = FillValueInt
            AuxJ(:) = FillValueInt

            nCell_ = 0

            CurrentLine => Line

dw:         do while (associated(CurrentLine))

d1:             do l = 1, CurrentLine%nNodes - 1

                    Segment%StartAt%X = CurrentLine%X(l)
                    Segment%StartAt%Y = CurrentLine%Y(l)

                    Segment%EndAt%X = CurrentLine%X(l+1)
                    Segment%EndAt%Y = CurrentLine%Y(l+1)


d2:                 do j = Me%WorkSize%JLB, Me%WorkSize%JUB
d3:                 do i = Me%WorkSize%ILB, Me%WorkSize%IUB

i2:                     if (Me%DefineCellsMap(i, j) == 1 .and. WaterPoints2D(i,j) == WaterPoint) then

                            Polygon%VerticesF(1)%X = XX2D(i, j)
                            Polygon%VerticesF(1)%Y = YY2D(i, j)

                            Polygon%VerticesF(2)%X = XX2D(i, j+1)
                            Polygon%VerticesF(2)%Y = YY2D(i, j+1)

                            Polygon%VerticesF(3)%X = XX2D(i+1, j+1)
                            Polygon%VerticesF(3)%Y = YY2D(i+1, j+1)

                            Polygon%VerticesF(4)%X = XX2D(i+1, j)
                            Polygon%VerticesF(4)%Y = YY2D(i+1, j)

                            Polygon%VerticesF(5)%X = Polygon%VerticesF(1)%X
                            Polygon%VerticesF(5)%Y = Polygon%VerticesF(1)%Y

                            if (Intersect2D_SegPoly(Segment, Polygon)) then

                                DoesCount =.true.
                                do k=1, nCell_
                                    if (i == AuxI(k) .and. j==AuxJ(k)) then
                                        DoesCount =.false.
                                        exit
                                    endif
                                enddo

                                if (DoesCount) then
                                    nCell_ = nCell_ + 1
                                    AuxI(nCell_) = i
                                    AuxJ(nCell_) = j
                                endif

                            endif


                        endif i2

                    enddo d3
                    enddo d2

                enddo d1

                CurrentLine => CurrentLine%Next

            enddo dw

            deallocate(Polygon%VerticesF)
            nullify   (Polygon%VerticesF)

            deallocate(Polygon)
            nullify   (Polygon)

            deallocate(Segment%StartAt)
            nullify   (Segment%StartAt)

            deallocate(Segment%EndAt)
            nullify   (Segment%EndAt)

            deallocate(Segment)
            nullify   (Segment)

            nullify   (CurrentLine)

            nCell = nCell_

            allocate(VectorI(nCell), VectorJ(nCell), VectorK(nCell))

            VectorI(:) = AuxI(1:nCell)
            VectorJ(:) = AuxJ(1:nCell)
            VectorK(:) = FillValueInt

            deallocate(AuxI); deallocate(AuxJ);

            STAT_ = SUCCESS_
        else    i1
            STAT_ = ready_
        end if  i1

        if (present(STAT)) STAT = STAT_

    end subroutine GetCellZInterceptByLine

    !--------------------------------------------------------------------------

   !--------------------------------------------------------------------------

    subroutine GetCellZInterceptByPolygon(HorizontalGridID, Polygon, WaterPoints2D,     &
                                          VectorI, VectorJ, VectorK, nCell, STAT)

        !Arguments---------------------------------------------------------------
        integer,            intent(IN)              :: HorizontalGridID
        type (T_Polygon), pointer                   :: Polygon
        integer, dimension(:,:), pointer            :: WaterPoints2D
        integer, dimension(:),   pointer            :: VectorI, VectorJ, VectorK
        integer                                     :: nCell
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        real,    dimension(:,:), pointer            :: XX2D, YY2D
        type (T_Polygon), pointer                   :: CurrentPolygon
        type (T_PointF),   pointer                  :: Point
        integer, dimension(:), pointer              :: AuxI, AuxJ
        integer                                     :: nCell_, i, j, k
        integer                                     :: STAT_
        integer                                     :: ready_
        logical                                     :: DoesCount

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

i1:     if ((ready_ == IDLE_ERR_     ) .OR.                                             &
            (ready_ == READ_LOCK_ERR_)) then

            allocate(Point)

            nCell_ = (Me%WorkSize%JUB-Me%WorkSize%JLB) * (Me%WorkSize%IUB-Me%WorkSize%ILB)

            allocate(AuxI(nCell_)); allocate(AuxJ(nCell_));

            nCell_ = 0

            if (Me%CoordType == SIMPLE_GEOG_ .or. Me%CoordType == GEOG_) then

                XX2D        => Me%LongitudeConn
                YY2D        => Me%LatitudeConn

            else

                XX2D        => Me%XX_IE
                YY2D        => Me%YY_IE

            endif

            CurrentPolygon => Polygon

dw:         do while (associated(CurrentPolygon))

d2:             do j = Me%WorkSize%JLB, Me%WorkSize%JUB
d3:             do i = Me%WorkSize%ILB, Me%WorkSize%IUB

i2:                 if (Me%DefineCellsMap(i, j) == 1 .and. WaterPoints2D(i,j) == WaterPoint) then

                        Point%X  = (XX2D(i  , j+1) + XX2D(i, j) + XX2D(i+1, j+1) + XX2D(i+1, j)) / 4.
                        Point%Y  = (YY2D(i  , j+1) + YY2D(i, j) + YY2D(i+1, j+1) + YY2D(i+1, j)) / 4.

                        if (IsPointInsidePolygon(Point, CurrentPolygon)) then

                            DoesCount =.true.
                            do k=1, nCell_
                                if (i == AuxI(k) .and. j==AuxJ(k)) then
                                    DoesCount =.false.
                                    exit
                                endif
                            enddo

                            if (DoesCount) then
                                nCell_ = nCell_ + 1
                                AuxI(nCell_) = i
                                AuxJ(nCell_) = j
                            endif

                        endif


                    endif i2

                enddo d3
                enddo d2

                CurrentPolygon => CurrentPolygon%Next

            enddo dw


            nCell = nCell_

            allocate(VectorI(nCell), VectorJ(nCell), VectorK(nCell))

            VectorI(:) = AuxI(1:nCell)
            VectorJ(:) = AuxJ(1:nCell)
            VectorK(:) = FillValueInt

            deallocate(AuxI); deallocate(AuxJ);
            deallocate(Point)

            nullify   (AuxI); nullify   (AuxJ);
            nullify   (Point)

            nullify(CurrentPolygon)

            STAT_ = SUCCESS_
        else    i1
            STAT_ = ready_
        end if  i1

        if (present(STAT)) STAT = STAT_

    end subroutine GetCellZInterceptByPolygon

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    logical function WindowIntersectDomain(HorizontalGridID, WorkSize2D, STAT)


        !Arguments-------------------------------------------------------------
        integer,                     intent(IN )    :: HorizontalGridID
        type (T_Size2D),             intent(IN )    :: WorkSize2D
        integer, optional,           intent(OUT)    :: STAT


        !External--------------------------------------------------------------

        integer :: ready_

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (Me%DDecomp%MasterOrSlave) then

                WindowIntersectDomain = WindowCellsIntersection(Me%DDecomp%HaloMap, WorkSize2D)

            else

                WindowIntersectDomain = WindowCellsIntersection(Me%WorkSize, WorkSize2D)

            endif

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end function WindowIntersectDomain

    !----------------------------------------------------------------------

    logical function WindowCellsIntersection(WorkSizeA, WorkSizeB)


        !Arguments-------------------------------------------------------------
        type (T_Size2D), intent(IN )    :: WorkSizeA, WorkSizeB

        !Local-----------------------------------------------------------------
        logical                         :: LineInterSection, ColumnInterSection

        !Begin-----------------------------------------------------------------

        LineInterSection        = .true.
        ColumnInterSection      = .true.

        if (WorkSizeA%ILB < WorkSizeB%ILB .and. WorkSizeA%IUB < WorkSizeB%ILB) then
            LineInterSection   = .false.
        endif

        if (WorkSizeA%ILB > WorkSizeB%IUB .and. WorkSizeA%IUB > WorkSizeB%IUB) then
            LineInterSection   = .false.
        endif

        if (WorkSizeA%JLB < WorkSizeB%JLB .and. WorkSizeA%JUB < WorkSizeB%JLB) then
            ColumnInterSection = .false.
        endif

        if (WorkSizeA%JLB > WorkSizeB%JUB .and. WorkSizeA%JUB > WorkSizeB%JUB) then
            ColumnInterSection = .false.
        endif

        if (LineInterSection .and. ColumnInterSection) then
            WindowCellsIntersection = .true.
        else
            WindowCellsIntersection = .false.
        endif
        !----------------------------------------------------------------------

    end function WindowCellsIntersection

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    type (T_Size2D) function ReturnsIntersectionCorners(HorizontalGridID, WorkSize2D, STAT)


        !Arguments-------------------------------------------------------------
        integer,                     intent(IN )    :: HorizontalGridID
        type (T_Size2D),             intent(IN )    :: WorkSize2D
        integer, optional,           intent(OUT)    :: STAT


        !External--------------------------------------------------------------

        integer :: ready_

        !Local-----------------------------------------------------------------
        type (T_Size2D)                             :: WorkSize2DAux
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (Me%DDecomp%MasterOrSlave) then

                WorkSize2DAux                  = ComputesIntersectionCorners(Me%DDecomp%HaloMap, WorkSize2D)

                ReturnsIntersectionCorners%ILB = WorkSize2DAux%ILB - Me%DDecomp%HaloMap%ILB + 1
                ReturnsIntersectionCorners%IUB = WorkSize2DAux%IUB - Me%DDecomp%HaloMap%ILB + 1

                ReturnsIntersectionCorners%JLB = WorkSize2DAux%JLB - Me%DDecomp%HaloMap%JLB + 1
                ReturnsIntersectionCorners%JUB = WorkSize2DAux%JUB - Me%DDecomp%HaloMap%JLB + 1

            else

                ReturnsIntersectionCorners = ComputesIntersectionCorners(Me%WorkSize, WorkSize2D)

            endif

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end function ReturnsIntersectionCorners

    !----------------------------------------------------------------------

    type (T_Size2D) function ComputesIntersectionCorners(WorkSizeA, WorkSizeB)


        !Arguments-------------------------------------------------------------
        type (T_Size2D), intent(IN )    :: WorkSizeA, WorkSizeB


        !Begin-----------------------------------------------------------------

        ComputesIntersectionCorners%ILB = max(WorkSizeA%ILB,WorkSizeB%ILB)
        ComputesIntersectionCorners%IUB = min(WorkSizeA%IUB,WorkSizeB%IUB)
        ComputesIntersectionCorners%JLB = max(WorkSizeA%JLB,WorkSizeB%JLB)
        ComputesIntersectionCorners%JUB = min(WorkSizeA%JUB,WorkSizeB%JUB)

        !----------------------------------------------------------------------

    end function ComputesIntersectionCorners

    !--------------------------------------------------------------------------

    subroutine GetDDecompParameters(HorizontalGridID, ON, Master,           &
                                                Master_MPI_ID, MasterOrSlave,           &
                                                NInterfaces, Interfaces, Halo_Points,   &
                                                MPI_ID, Global, Mapping, Inner, HaloMap,&
                                                STAT)


        !Arguments-------------------------------------------------------------
        integer,                     intent(IN )    :: HorizontalGridID
        logical, optional,           intent(OUT)    :: ON
        logical, optional,           intent(OUT)    :: Master
        integer, optional,           intent(OUT)    :: Master_MPI_ID
        logical, optional,           intent(OUT)    :: MasterOrSlave
        integer, optional,           intent(OUT)    :: NInterfaces
        integer, optional, dimension(:,:), pointer  :: Interfaces
        integer, optional,           intent(OUT)    :: Halo_Points
        integer, optional,           intent(OUT)    :: MPI_ID
        type (T_Size2D), optional,   intent(OUT)    :: Global
        type (T_Size2D), optional,   intent(OUT)    :: Mapping
        type (T_Size2D), optional,   intent(OUT)    :: Inner
        type (T_Size2D), optional,   intent(OUT)    :: HaloMap
        integer, optional,           intent(OUT)    :: STAT


        !External--------------------------------------------------------------

        integer :: ready_

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(ON)) then
                ON            = Me%DDecomp%ON
            endif

            if (present(Master)) then
              Master          = Me%DDecomp%Master
            endif

            if (present(Master_MPI_ID)) then
              Master_MPI_ID   = Me%DDecomp%Master_MPI_ID
            endif

            if (present(MasterOrSlave)) then
                MasterOrSlave = Me%DDecomp%MasterOrSlave
            endif

            if (present(NInterfaces)) then
                NInterfaces   = Me%DDecomp%NInterfaces
            endif

            if (present(Interfaces)) then
                call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
                Interfaces   => Me%DDecomp%Interfaces
            endif

            if (present(Halo_Points)) then
                Halo_Points   = Me%DDecomp%Halo_Points
            endif

            if (present(MPI_ID)) then
                MPI_ID        = Me%DDecomp%MPI_ID
            endif

            if (present(Global)) then
                Global       = Me%DDecomp%Global
            endif

            if (present(Mapping)) then
                Mapping      = Me%DDecomp%Mapping
            endif

            if (present(Inner)) then
                Inner        = Me%DDecomp%Inner
            endif

            if (present(HaloMap)) then
                HaloMap       = Me%DDecomp%HaloMap
            endif

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDDecompParameters

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetDDecompSlaves(HorizontalGridID, Nslaves, Slaves_MPI_ID, STAT)


        !Arguments-------------------------------------------------------------
        integer,                     intent(IN )    :: HorizontalGridID
        integer,                     intent(OUT)    :: Nslaves
        integer, dimension(:),       pointer        :: Slaves_MPI_ID
        integer, optional,           intent(OUT)    :: STAT


        !External--------------------------------------------------------------

        integer :: ready_

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Nslaves          = Me%DDecomp%Nslaves

            call Read_Lock(mHORIZONTALGRID_, Me%InstanceID)
            Slaves_MPI_ID   => Me%DDecomp%Slaves_MPI_ID

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDDecompSlaves

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetDDecompSlavesSize(HorizontalGridID, iSlave,               &
                                                Slaves_Inner, Slaves_Size,              &
                                                Slaves_Mapping, Slaves_HaloMap, STAT)


        !Arguments-------------------------------------------------------------
        integer,                     intent(IN )    :: HorizontalGridID
        integer,                     intent(IN )    :: iSlave
        type (T_Size2D), optional,   intent(OUT)    :: Slaves_Inner
        type (T_Size2D), optional,   intent(OUT)    :: Slaves_Size
        type (T_Size2D), optional,   intent(OUT)    :: Slaves_Mapping
        type (T_Size2D), optional,   intent(OUT)    :: Slaves_HaloMap
        integer, optional,           intent(OUT)    :: STAT


        !External--------------------------------------------------------------

        integer :: ready_

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(Slaves_Inner)) then
                Slaves_Inner   = Me%DDecomp%Slaves_Inner   (iSlave)
            endif

            if (present(Slaves_Size)) then
                Slaves_Size    = Me%DDecomp%Slaves_Size    (iSlave)
            endif

           if (present(Slaves_Mapping)) then
                Slaves_Mapping = Me%DDecomp%Slaves_Mapping (iSlave)
            endif

            if (present(Slaves_HaloMap)) then
                Slaves_HaloMap = Me%DDecomp%Slaves_HaloMap (iSlave)
            endif

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDDecompSlavesSize

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    integer function GetDDecompMPI_ID(HorizontalGridID, STAT)

        !Arguments-------------------------------------------------------------
        integer,                     intent(IN )    :: HorizontalGridID
        integer,   optional,         intent(OUT)    :: STAT

        !External--------------------------------------------------------------

        integer :: ready_

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        GetDDecompMPI_ID = null_int

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            GetDDecompMPI_ID  = Me%DDecomp%MPI_ID

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end function GetDDecompMPI_ID

    !--------------------------------------------------------------------------
     !--------------------------------------------------------------------------

    logical function GetDDecompON(HorizontalGridID, STAT)

        !Arguments-------------------------------------------------------------
        integer,                     intent(IN )    :: HorizontalGridID
        integer,   optional,         intent(OUT)    :: STAT

        !External--------------------------------------------------------------

        integer :: ready_

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        GetDDecompON = .false.

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            GetDDecompON  = Me%DDecomp%MasterOrSlave

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end function GetDDecompON

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    function GetDDecompOpenBorders(HorizontalGridID, STAT)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )              :: HorizontalGridID
        integer, optional, intent(OUT)              :: STAT


        !External--------------------------------------------------------------
        logical,  dimension(1:4)                    :: GetDDecompOpenBorders
        integer :: ready_

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            GetDDecompOpenBorders(1:4) = Me%DDecomp%OpenBordersON(1:4)

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end function GetDDecompOpenBorders

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine GetDDecompWorkSize2D(HorizontalGridID, WorkSize, STAT)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )              :: HorizontalGridID
        type(T_Size2D)   , intent(OUT)              :: WorkSize
        integer, optional, intent(OUT)              :: STAT


        !External--------------------------------------------------------------
        integer :: ready_

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            WorkSize%ILB = Me%DDecomp%HaloMap%ILB
            WorkSize%IUB = Me%DDecomp%HaloMap%IUB
            WorkSize%JLB = Me%DDecomp%HaloMap%JLB
            WorkSize%JUB = Me%DDecomp%HaloMap%JUB

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDDecompWorkSize2D

    !--------------------------------------------------------------------------

    subroutine GetDDecompMapping2D(HorizontalGridID, Mapping2D, STAT)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )              :: HorizontalGridID
        type(T_Size2D)   , intent(OUT)              :: Mapping2D
        integer, optional, intent(OUT)              :: STAT


        !External--------------------------------------------------------------
        integer :: ready_

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Mapping2D%ILB = Me%DDecomp%Mapping%ILB
            Mapping2D%IUB = Me%DDecomp%Mapping%IUB
            Mapping2D%JLB = Me%DDecomp%Mapping%JLB
            Mapping2D%JUB = Me%DDecomp%Mapping%JUB

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDDecompMapping2D

    !--------------------------------------------------------------------------

    subroutine UngetHorizontalGrid1D(HorizontalGridID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real, pointer, dimension(:)                 :: Array
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ == IDLE_ERR_     ) .OR.                                              &
            (ready_ == READ_LOCK_ERR_)) then

            nullify(Array)
            call Read_UnLock(mHORIZONTALGRID_, Me%InstanceID, "UngetHorizontalGrid1D")

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetHorizontalGrid1D

    !--------------------------------------------------------------------------

    subroutine UngetHorizontalGrid2D(HorizontalGridID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        real, pointer, dimension(:,:)               :: Array
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ == IDLE_ERR_     ) .OR.                                              &
            (ready_ == READ_LOCK_ERR_)) then

            nullify(Array)

            call Read_UnLock(mHORIZONTALGRID_, Me%InstanceID, "UngetHorizontalGrid2D")

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetHorizontalGrid2D

    !--------------------------------------------------------------------------

    subroutine UngetHorizontalGrid2DInt(HorizontalGridID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer, pointer, dimension(:,:)            :: Array
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ == IDLE_ERR_     ) .OR.                                              &
            (ready_ == READ_LOCK_ERR_)) then

            nullify(Array)

            call Read_UnLock(mHORIZONTALGRID_, Me%InstanceID, "UngetHorizontalGrid2DInt")

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetHorizontalGrid2DInt

    !--------------------------------------------------------------------------

    subroutine UngetHorizontalGridPolygon(HorizontalGridID, polygon, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        type (T_Polygon), pointer                   :: polygon
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ == IDLE_ERR_     ) .OR.                                              &
            (ready_ == READ_LOCK_ERR_)) then

            nullify(polygon)

            call Read_UnLock(mHORIZONTALGRID_, Me%InstanceID, "UngetHorizontalGridPolygon")

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetHorizontalGridPolygon

    subroutine UngetHorizontalGrid1DInt(HorizontalGridID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer, pointer, dimension(:)              :: Array
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_
        integer                                     :: ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ == IDLE_ERR_     ) .OR.                                              &
            (ready_ == READ_LOCK_ERR_)) then

            nullify(Array)
            call Read_UnLock(mHORIZONTALGRID_, Me%InstanceID, "UngetHorizontalGrid1DInt")

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetHorizontalGrid1DInt

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    real  function InterpolXYPoint(HorizontalGridID, Field2DFather,                     &
                                   ComputeFather, XInput, YInput,                       &
                                   Compute, STAT)

        !Arguments-------------------------------------------------------------
        integer                                             :: HorizontalGridID
        real,              dimension(:,:  ), pointer        :: Field2DFather
        integer, optional, dimension(:,:  ), pointer        :: ComputeFather
        real,                                intent(IN)     :: XInput, YInput
        integer, optional,                   intent(IN)     :: Compute
        integer, optional,                   intent(OUT)    :: STAT

        !Local-----------------------------------------------------------------
        real   , pointer, dimension(:,:)                    :: DXFather, DYFather, XXFather2D, YYFather2D
        real   , pointer, dimension(:  )                    :: XXFather, YYFather
        integer, pointer, dimension(:,:)                    :: DefinedPoint
        logical                                             :: InsideDomain

        real                                                :: YYUpper, YYLower, XXUpper, XXLower,  &
                                                               PropLowLeft,  PropUpLeft, PropLowRight, PropUpRight

        integer                                             :: ONLowLeft, ONUpLeft, ONLowRight, ONUpRight
        integer                                             :: ready_

        integer                                             :: Jlower, Jupper, Ilower, Iupper
        integer                                             :: STAT_

        integer                                             :: I, J
        real                                                :: PercI, PercJ
        integer                                             :: Compute_
        logical                                             :: InterPolOK
        real                                                :: XPosition, YPosition, InterpolatedValue
        integer                                             :: Dij, Dji, JUBFather, IUBFather

        !Begin------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then


            if (.NOT. associated(Field2DFather))                                        &
                call SetError(FATAL_, INTERNAL_,                                        &
                             "InterpolXYPoint; HorizontalGrid. ERR10")
            if (present(Compute)) then
                Compute_ = Compute
            else
                Compute_ = ComputeZ_
            endif

            InsideDomain = InsideDomainPolygon(Me%GridBorderCoord%Polygon_, XInput, YInput)

id:         if (InsideDomain) then

                call  GetXYCellZ(Me%InstanceID, XInput, YInput, I, J, PercI = PercI, PercJ = PercJ)

                if      (Compute_ == ComputeU_)  then

                    JUBFather = JUBFather + 1

                    XXFather => Me%Compute%XX_U
                    YYFather => Me%Compute%YY_U

                    XXFather2D => Me%Compute%XX2D_U
                    YYFather2D => Me%Compute%YY2D_U

                    DXFather => Me%DXX
                    DYFather => Me%DZY

                    Dij = 0
                    Dji = 1

                    Jlower = J
                    Jupper = J + 1

                    if (PercI < 0.5) then
                        Ilower = I-1
                    else
                        Ilower = I
                    endif

                    Iupper = Ilower + 1

                    if (Me%CornersXYInput) DefinedPoint => Me%DefineFacesUMap


                else if (Compute_ == ComputeV_)  then

                    IUBFather = IUBFather + 1

                    XXFather => Me%Compute%XX_V
                    YYFather => Me%Compute%YY_V

                    XXFather2D => Me%Compute%XX2D_V
                    YYFather2D => Me%Compute%YY2D_V

                    DXFather => Me%DZX
                    DYFather => Me%DYY

                    Dij = 1
                    Dji = 0

                    if (Me%CornersXYInput) DefinedPoint => Me%DefineFacesVMap

                    Ilower = I
                    Iupper = I + 1

                    if (PercJ < 0.5) then
                        Jlower = J-1
                    else
                        Jlower = J
                    endif

                    Jupper = Jlower + 1


                else if (Compute_ == ComputeZ_)  then

                    XXFather => Me%Compute%XX_Z
                    YYFather => Me%Compute%YY_Z

                    XXFather2D => Me%Compute%XX2D_Z
                    YYFather2D => Me%Compute%YY2D_Z

                    DXFather => Me%DVX
                    DYFather => Me%DUY

                    Dji = 1
                    Dij = 1

                    if (Me%CornersXYInput) DefinedPoint => Me%DefineCellsMap

                    if (PercI < 0.5) then
                        Ilower = I-1
                    else
                        Ilower = I
                    endif

                    Iupper = Ilower + 1

                    if (PercJ < 0.5) then
                        Jlower = J-1
                    else
                        Jlower = J
                    endif

                    Jupper = Jlower + 1


                else if (Compute_ == ComputeCross_)  then

                    IUBFather = IUBFather + 1
                    JUBFather = JUBFather + 1

                    XXFather => Me%Compute%XX_Cross
                    YYFather => Me%Compute%YY_Cross

                    XXFather2D => Me%XX_IE
                    YYFather2D => Me%YY_IE

                    DXFather => Me%DUX
                    DYFather => Me%DVY

                    Dji = 0
                    Dij = 0

                    if (Me%CornersXYInput) DefinedPoint => Me%DefineCrossMap

                    Ilower = I
                    Iupper = I + 1
                    Jlower = J
                    Jupper = J + 1


                endif

                if (.not. Me%CornersXYInput) then

                    XXLower  = XXFather(Jlower)
                    XXUpper  = XXFather(Jupper)

                    YYLower  = YYFather(Ilower)
                    YYUpper  = YYFather(Iupper)

                    XPosition = XInput
                    YPosition = YInput

                    call RODAXY(-Me%Xorig, -Me%Yorig, -Me%Grid_Angle, XPosition, YPosition)

                else

                    call CellReferential(DXFather, DYFather, XXFather2D, YYFather2D,    &
                                         Me%RotationX,                &
                                         Me%RotationY,                &
                                         XInput, YInput,                                &
                                         Ilower, Iupper, Jlower, Jupper, dij, dji,      &
                                         XPosition, YPosition,                          &
                                         YYUpper, YYLower, XXUpper, XXLower)

                endif

                PropLowLeft = Field2DFather(Ilower, Jlower)
                PropUpLeft  = Field2DFather(Iupper, Jlower)
                PropLowRight= Field2DFather(Ilower, Jupper)
                PropUpRight = Field2DFather(Iupper, Jupper)

                if (present(ComputeFather)) then

                    ONLowLeft   = ComputeFather(Ilower, Jlower)
                    ONUpLeft    = ComputeFather(Iupper, Jlower)
                    ONLowRight  = ComputeFather(Ilower, Jupper)
                    ONUpRight   = ComputeFather(Iupper, Jupper)

                else

                    ONLowLeft   = 1
                    ONUpLeft    = 1
                    ONLowRight  = 1
                    ONUpRight   = 1

                endif

                if (ONLowLeft + ONUpLeft + ONLowRight + ONUpRight == 4) then

                    call InterpolPoint(XPosition, YPosition,                             &
                                       YYUpper, YYLower, XXUpper, XXLower,               &
                                       PropLowLeft, PropUpLeft, PropLowRight, PropUpRight,&
                                       ONLowLeft,   ONUpLeft,   ONLowRight,   ONUpRight, &
                                       InterpolatedValue, InterPolOK)

                    if (InterPolOK)   then

                        InterpolXYPoint = InterpolatedValue

                    else

                        InterpolXYPoint = FillValueReal

                    endif

                else

                    InterpolXYPoint = FillValueReal

                endif

                nullify(XXFather,   YYFather  )
                nullify(XXFather2D, YYFather2D)
                nullify(DXFather,   DYFather  )

            else  id

                InterpolXYPoint = FillValueReal

            endif id

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end function InterpolXYPoint

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine InterpolRegularGrid2D(HorizontalGridSonID, HorizontalGridFatherID,        &
                                     Field2DFather, Field2DSon,                          &
                                     ComputeFather, Compute,                             &
                                     KUBFather, GridID, STAT)



        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridSonID
        integer                                     :: HorizontalGridFatherID
        real,    dimension(:,:  ), pointer          :: Field2DFather, Field2DSon
        integer, dimension(:,:,:), pointer, optional:: ComputeFather
        integer, optional, intent(IN)               :: Compute, KUBFather
        integer, optional, intent(IN)               :: GridID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        type(T_HorizontalGrid),    pointer          :: ObjHorizontalGridFather
        type (T_FatherGrid),       pointer          :: FatherGrid
        real   , pointer, dimension(:,:)            :: XXSon, YYSon, DXFather, DYFather, XXFather2D, YYFather2D
        real   , pointer, dimension(:  )            :: XXFather, YYFather
        integer, pointer, dimension(:,:)            :: JSon, ISon, DefinedPoint


        real    :: YYUpper, YYLower, XXUpper, XXLower,                                   &
                   PropLowLeft,  PropUpLeft, PropLowRight, PropUpRight

        integer :: ONLowLeft, ONUpLeft, ONLowRight, ONUpRight
        integer :: ready_, GridID_

        integer :: JLBSon, JUBSon, ILBSon, IUBSon
        integer :: JLBFather, JUBFather, ILBFather, IUBFather
        integer :: Jlower, Jupper, Ilower, Iupper
        integer :: STAT_

        integer :: i, j
        integer :: Compute_, KUBFather_
        logical :: InterPolOK, MapFatherGrid

        real    :: XPosition, YPosition
        integer :: Dij, Dji

        !$ integer :: CHUNK

        !Begin------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleHorizontalGrid", "InterpolRegularGrid2D")

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridSonID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then


            if (.NOT. associated(Field2DFather))                                         &
                call SetError(FATAL_, INTERNAL_,                                           &
                             "InterpolRegularGrid2D; HorizontalGrid. ERR01")

            if (.NOT. associated(Field2DSon))                                            &
                call SetError(FATAL_, INTERNAL_,                                           &
                             "InterpolRegularGrid2D; HorizontalGrid. ERR02")

            nullify(ObjHorizontalGridFather)
            call LocateObjFather (ObjHorizontalGridFather, HorizontalGridFatherID)

            JLBSon = Me%WorkSize%JLB
            JUBSon = Me%WorkSize%JUB
            ILBSon = Me%WorkSize%ILB
            IUBSon = Me%WorkSize%IUB

            if (present(GridID)) then
                GridID_ = GridID
            else
                GridID_ = LargeScaleModel_
            endif

            if (present(ComputeFather)) then
                MapFatherGrid = .true.
            else
                MapFatherGrid = .false.
            endif

            if (present(Compute)) then
                Compute_ = Compute
            else
                Compute_ = ComputeZ_
            endif

            if (present(KUBFather)) then
                KUBFather_ = KUBFather
            else
                KUBFather_ = 1
            endif

            call Search_FatherGrid (FatherGrid, GridID_)


            !Gets the bounds from the Father
            ILBFAther = ObjHorizontalGridFather%WorkSize%ILB
            IUBFather = ObjHorizontalGridFather%WorkSize%IUB
            JLBFather = ObjHorizontalGridFather%WorkSize%JLB
            JUBFather = ObjHorizontalGridFather%WorkSize%JUB

            if      (Compute_ == ComputeU_)  then

                JUBFather = JUBFather + 1
                JUBSon    = JUBSon    + 1

                XXSon    => FatherGrid%XX_U
                YYSon    => FatherGrid%YY_U

                JSon     => FatherGrid%JU
                ISon     => FatherGrid%IU

                XXFather => ObjHorizontalGridFather%Compute%XX_U
                YYFather => ObjHorizontalGridFather%Compute%YY_U

                XXFather2D => ObjHorizontalGridFather%Compute%XX2D_U
                YYFather2D => ObjHorizontalGridFather%Compute%YY2D_U

                DXFather => ObjHorizontalGridFather%DXX
                DYFather => ObjHorizontalGridFather%DZY

                Dij = 0
                Dji = 1

                if (Me%CornersXYInput) DefinedPoint => Me%DefineFacesUMap


            else if (Compute_ == ComputeV_)  then

                IUBFather = IUBFather + 1
                IUBSon    = IUBSon    + 1

                XXSon => FatherGrid%XX_V
                YYSon => FatherGrid%YY_V

                JSon     => FatherGrid%JV
                ISon     => FatherGrid%IV

                XXFather => ObjHorizontalGridFather%Compute%XX_V
                YYFather => ObjHorizontalGridFather%Compute%YY_V

                XXFather2D => ObjHorizontalGridFather%Compute%XX2D_V
                YYFather2D => ObjHorizontalGridFather%Compute%YY2D_V

                DXFather => ObjHorizontalGridFather%DZX
                DYFather => ObjHorizontalGridFather%DYY

                Dij = 1
                Dji = 0

                if (Me%CornersXYInput) DefinedPoint => Me%DefineFacesVMap

            else if (Compute_ == ComputeZ_)  then

                XXSon => FatherGrid%XX_Z
                YYSon => FatherGrid%YY_Z

                XXFather => ObjHorizontalGridFather%Compute%XX_Z
                YYFather => ObjHorizontalGridFather%Compute%YY_Z

                JSon     => FatherGrid%JZ
                ISon     => FatherGrid%IZ

                XXFather2D => ObjHorizontalGridFather%Compute%XX2D_Z
                YYFather2D => ObjHorizontalGridFather%Compute%YY2D_Z

                DXFather => ObjHorizontalGridFather%DVX
                DYFather => ObjHorizontalGridFather%DUY

                Dji = 1
                Dij = 1

                if (Me%CornersXYInput) DefinedPoint => Me%DefineCellsMap

            else if (Compute_ == ComputeCross_)  then

                IUBFather = IUBFather + 1
                JUBFather = JUBFather + 1

                IUBSon    = IUBSon

                XXSon => FatherGrid%XX_Z
                YYSon => FatherGrid%YY_Z

                XXFather => ObjHorizontalGridFather%Compute%XX_Cross
                YYFather => ObjHorizontalGridFather%Compute%YY_Cross

                JSon     => FatherGrid%JCross
                ISon     => FatherGrid%ICross

                XXFather2D => ObjHorizontalGridFather%XX_IE
                YYFather2D => ObjHorizontalGridFather%YY_IE

                DXFather => ObjHorizontalGridFather%DUX
                DYFather => ObjHorizontalGridFather%DVY

                Dji = 0
                Dij = 0

                if (Me%CornersXYInput) DefinedPoint => Me%DefineCrossMap

            else
                write(*,*) 'Compute ==', Compute
                stop "InterpolRegularGrid2D; HorizontalGrid. ERR10"

            endif

            if (.not.MapFatherGrid) then

                ONLowLeft   = 1
                ONUpLeft    = 1
                ONLowRight  = 1
                ONUpRight   = 1

            endif


            !$ CHUNK = CHUNK_J(JLBSon,JUBSon)

            !griflet
            !! $OMP PARALLEL PRIVATE( j,i,                      &
            !! $OMP                   Jlower, Jupper,             &
            !! $OMP                   Ilower,Iupper,              &
            !! $OMP                   XXLower,XXUpper,            &
            !! $OMP                   YYLower,YYUpper,            &
            !! $OMP                   XPosition,YPosition,        &
            !! $OMP                   PropLowLeft,PropUpLeft,     &
            !! $OMP                   PropLowRight,PropUpRight,   &
            !! $OMP                   ONLowLeft,ONUpLeft,         &
            !! $OMP                   ONLowRight,ONUpRight,       &
            !! $OMP                   InterPolOK)
            !! $OMP DO SCHEDULE(DYNAMIC,CHUNK)
doj:        do j = JLBSon, JUBSon
doi:        do i = ILBSon, IUBSon

                if (Me%CornersXYInput) then
                    if(DefinedPoint(i,j) == 0) cycle
                endif

                Jlower = JSon(i, j)

                !Father domain smaller than son domain
                if (Jlower < -100) then
                    Field2DSon(i, j) = FillValueReal
                    cycle
                endif

                Jupper = JSon(i, j) + 1

                Ilower = ISon(i, j)

                !Father domain smaller than son domain
                if (Ilower < -100) then
                    Field2DSon(i, j) = FillValueReal
                    cycle
                endif

                Iupper = ISon(i, j) + 1

                if (Jupper > JUBFather) Jupper = JUBFather
                if (Iupper > IUBFather) Iupper = IUBFather
                if (Jlower < JLBFather) Jlower = JLBFather
                if (Ilower < ILBFather) Ilower = ILBFather

                if (Jupper < JLBFather) Jupper = JLBFather
                if (Iupper < ILBFather) Iupper = ILBFather

                if (.not. ObjHorizontalGridFather%CornersXYInput) then

                    XXLower  = XXFather(Jlower)
                    XXUpper  = XXFather(Jupper)

                    YYLower  = YYFather(Ilower)
                    YYUpper  = YYFather(Iupper)

                    XPosition = XXSon(i, j)
                    YPosition = YYSon(i, j)

                else

                    !griflet: All write argument variables are scalars and are all declared
                    !griflet: PRIVATE. Clear! Ok!
                    call CellReferential(DXFather, DYFather, XXFather2D, YYFather2D,    &
                                         ObjHorizontalGridFather%RotationX,             &
                                         ObjHorizontalGridFather%RotationY,             &
                                         XXSon(i,j), YYSon(i,j),                        &
                                         Ilower, Iupper, Jlower, Jupper, dij, dji,      &
                                         XPosition, YPosition,                          &
                                         YYUpper, YYLower, XXUpper, XXLower)

                endif

                PropLowLeft = Field2DFather(Ilower, Jlower)
                PropUpLeft  = Field2DFather(Iupper, Jlower)
                PropLowRight= Field2DFather(Ilower, Jupper)
                PropUpRight = Field2DFather(Iupper, Jupper)

                if (MapFatherGrid) then

                    ONLowLeft   = ComputeFather(Ilower, Jlower, KUBFather_)
                    ONUpLeft    = ComputeFather(Iupper, Jlower, KUBFather_)
                    ONLowRight  = ComputeFather(Ilower, Jupper, KUBFather_)
                    ONUpRight   = ComputeFather(Iupper, Jupper, KUBFather_)

                endif

                !griflet: All argument variables were declared PRIVATE: clear! Ok.
                call InterpolPoint(XPosition, YPosition,                             &
                                   YYUpper, YYLower, XXUpper, XXLower,               &
                                   PropLowLeft, PropUpLeft, PropLowRight, PropUpRight,&
                                   ONLowLeft,   ONUpLeft,   ONLowRight,   ONUpRight, &
                                   Field2DSon(i, j), InterPolOK)

ifc:            if (.not. InterPolOK)   then

                    Field2DSon(i, j) = FillValueReal

                endif ifc


            enddo doi
            enddo doj
            !! $OMP END DO NOWAIT
            !! $OMP END PARALLEL

            nullify(XXSon,    YYSon   )
            nullify(JSon,     ISon    )
            nullify(XXFather, YYFather)


            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        if (MonitorPerformance) call StopWatch ("ModuleHorizontalGrid", "InterpolRegularGrid2D")

        !----------------------------------------------------------------------

    end subroutine InterpolRegularGrid2D

    !--------------------------------------------------------------------------


    subroutine InterpolRegularGrid3D(HorizontalGridSonID, HorizontalGridFatherID,         &
                                      Field3DFather, Field3DSon,                          &
                                      ComputeFather, Compute,                             &
                                      KLBFather, KUBFather, KUBSon,                       &
                                      FluxType, GridID, STAT)



        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridSonID
        integer                                     :: HorizontalGridFatherID
        real(4), dimension(:,:,:), pointer          :: Field3DFather, Field3DSon
        integer, dimension(:,:,:), pointer          :: ComputeFather
        integer,           intent(IN)               :: Compute, KLBFather, KUBFather, KUBSon
        logical, optional, intent(IN)               :: FluxType
        integer, optional, intent(IN)               :: GridID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        type(T_HorizontalGrid),    pointer          :: ObjHorizontalGridFather
        type (T_FatherGrid),       pointer          :: FatherGrid
        real,    dimension(:,:  ), pointer          :: DXFather, DYFather, XXFather2D, YYFather2D
        real   , pointer, dimension(:,:)            :: XXSon, YYSon, DYXFather
        real   , pointer, dimension(:  )            :: XXFather, YYFather
        integer, pointer, dimension(:,:)            :: JSon, ISon, DefinedPoint

        real    :: YYUpper, YYLower, XXUpper, XXLower,                                   &
                   PropLowLeft,  PropUpLeft, PropLowRight, PropUpRight, PropSon

        integer :: ONLowLeft, ONUpLeft, ONLowRight, ONUpRight
        integer :: ready_, GridID_

        integer :: JLBSon, JUBSon, ILBSon, IUBSon
        integer :: JLBFather, JUBFather, ILBFather, IUBFather
        integer :: Jlower, Jupper, Ilower, Iupper
        integer :: STAT_

        integer :: i, j, k
        logical :: InterPolOK, FluxType_

        real    :: XPosition, YPosition
        integer :: Dij, Dji

        !$ integer  :: CHUNK

        !Begin------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleHorizontalGrid", "InterpolRegularGrid3D")

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridSonID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then

            if (.NOT. associated(Field3DFather))                                         &
                call SetError(FATAL_, INTERNAL_,                                           &
                             "InterpolRegularGrid3D; HorizontalGrid. ERR01")

            if (.NOT. associated(Field3DSon))                                            &
                call SetError(FATAL_, INTERNAL_,                                           &
                             "InterpolRegularGrid3D; HorizontalGrid. ERR02")

            nullify(ObjHorizontalGridFather)
            call LocateObjFather (ObjHorizontalGridFather, HorizontalGridFatherID)

            JLBSon = Me%WorkSize%JLB
            JUBSon = Me%WorkSize%JUB
            ILBSon = Me%WorkSize%ILB
            IUBSon = Me%WorkSize%IUB

            if (present(GridID)) then

                GridID_ = GridID

            else

                GridID_ = LargeScaleModel_

            endif

            call Search_FatherGrid (FatherGrid, GridID_)


            !Gets the bounds from the Father
            ILBFAther = ObjHorizontalGridFather%WorkSize%ILB
            IUBFather = ObjHorizontalGridFather%WorkSize%IUB
            JLBFather = ObjHorizontalGridFather%WorkSize%JLB
            JUBFather = ObjHorizontalGridFather%WorkSize%JUB


            if (present(FluxType)) then
                FluxType_ = FluxType
            else
                FluxType_ = .false.
            endif

            if      (Compute == ComputeU_)  then

                JUBFather = JUBFather + 1
                JUBSon    = JUBSon    + 1

                XXSon    => FatherGrid%XX_U
                YYSon    => FatherGrid%YY_U

                JSon     => FatherGrid%JU
                ISon     => FatherGrid%IU


                XXFather => ObjHorizontalGridFather%Compute%XX_U
                YYFather => ObjHorizontalGridFather%Compute%YY_U

                DYXFather=> ObjHorizontalGridFather%DYY

                XXFather2D => ObjHorizontalGridFather%Compute%XX2D_U
                YYFather2D => ObjHorizontalGridFather%Compute%YY2D_U

                DXFather => ObjHorizontalGridFather%DXX
                DYFather => ObjHorizontalGridFather%DZY

                Dij = 0
                Dji = 1

                if (Me%CornersXYInput) DefinedPoint => Me%DefineFacesUMap

            else if (Compute == ComputeV_)  then

                IUBFather = IUBFather + 1
                IUBSon    = IUBSon    + 1

                XXSon => FatherGrid%XX_V
                YYSon => FatherGrid%YY_V

                JSon     => FatherGrid%JV
                ISon     => FatherGrid%IV

                XXFather => ObjHorizontalGridFather%Compute%XX_V
                YYFather => ObjHorizontalGridFather%Compute%YY_V

                DYXFather=> ObjHorizontalGridFather%DXX

                XXFather2D => ObjHorizontalGridFather%Compute%XX2D_V
                YYFather2D => ObjHorizontalGridFather%Compute%YY2D_V

                DXFather => ObjHorizontalGridFather%DZX
                DYFather => ObjHorizontalGridFather%DYY

                Dij = 1
                Dji = 0

                if (Me%CornersXYInput) DefinedPoint => Me%DefineFacesVMap


            else if (Compute == ComputeZ_)  then

                if (FluxType_) call SetError(FATAL_, INTERNAL_,                            &
                             "InterpolRegularGrid3D; HorizontalGrid. ERR04")

                XXSon => FatherGrid%XX_Z
                YYSon => FatherGrid%YY_Z

                XXFather => ObjHorizontalGridFather%Compute%XX_Z
                YYFather => ObjHorizontalGridFather%Compute%YY_Z

                JSon     => FatherGrid%JZ
                ISon     => FatherGrid%IZ

                XXFather2D => ObjHorizontalGridFather%Compute%XX2D_Z
                YYFather2D => ObjHorizontalGridFather%Compute%YY2D_Z

                DXFather => ObjHorizontalGridFather%DVX
                DYFather => ObjHorizontalGridFather%DUY

                Dji = 1
                Dij = 1

                if (Me%CornersXYInput) DefinedPoint => Me%DefineCellsMap

            else

                write(*,*) 'Compute ==', Compute
                stop "InterpolRegularGrid3D; HorizontalGrid. ERR10"

            endif

            !$ CHUNK = CHUNK_J(JLBSon,JUBSon)

            !griflet
            !! $OMP PARALLEL PRIVATE( k,j,i,                     &
            !! $OMP                   InterPolOK, Jlower, Jupper, &
            !! $OMP                   Ilower,Iupper,              &
            !! $OMP                   XXLower,XXUpper,            &
            !! $OMP                   YYLower,YYUpper,            &
            !! $OMP                   XPosition,YPosition,        &
            !! $OMP                   PropLowLeft,PropUpLeft,     &
            !! $OMP                   PropLowRight,PropUpRight,   &
            !! $OMP                   ONLowLeft,ONUpLeft,         &
            !! $OMP                   ONLowRight,ONUpRight,PropSon)
            !! $OMP           SHARED ( DXFather,DYFather,         &
            !! $OMP                   XXFather2D, YYFather2D,     &
            !! $OMP                   XXSon, YYSon)
dok:        do k = KLBFather   , KUBFather
            !! $OMP DO SCHEDULE(DYNAMIC,CHUNK)
doj:        do j = JLBSon, JUBSon
doi:        do i = ILBSon, IUBSon

                if (Me%CornersXYInput) then
                    if(DefinedPoint(i,j) == 0) cycle
                endif

                InterPolOK    = .false.

                Jlower = JSon(i, j)

                !Father domain smaller than son domain
                if (Jlower < -100) then
                    Field3DSon(i, j, k) = FillValueReal
                    cycle
                endif


                Jupper = JSon(i, j) + 1

                Ilower = ISon(i, j)

                !Father domain smaller than son domain
                if (Ilower < -100) then
                    Field3DSon(i, j, k) = FillValueReal
                    cycle
                endif


                Iupper = ISon(i, j) + 1

                if (.not. ObjHorizontalGridFather%CornersXYInput) then

                    XXLower  = XXFather(Jlower)
                    XXUpper  = XXFather(Jupper)

                    YYLower  = YYFather(Ilower)
                    YYUpper  = YYFather(Iupper)

                    XPosition = XXSon(i,j)
                    YPosition = YYSon(i,j)


                else

                    !griflet: no problem: this subroutine only alters private scalar variables
                    !griflet: All the arrays are for reading-only.
                    !griflet: Written: XPosition, YPosition, XX/YYUpper, XX/YYLower
                    call CellReferential(DXFather, DYFather, XXFather2D, YYFather2D,    &
                                         ObjHorizontalGridFather%RotationX,             &
                                         ObjHorizontalGridFather%RotationY,             &
                                         XXSon(i,j), YYSon(i,j),                        &
                                         Ilower, Iupper, Jlower, Jupper, dij, dji,      &
                                         XPosition, YPosition,                          &
                                         YYUpper, YYLower, XXUpper, XXLower)

                endif


                if (FluxType_) then

                    PropLowLeft = Field3DFather(Ilower, Jlower, k) / DYXFather(Ilower, Jlower)
                    PropUpLeft  = Field3DFather(Iupper, Jlower, k) / DYXFather(Iupper, Jlower)
                    PropLowRight= Field3DFather(Ilower, Jupper, k) / DYXFather(Ilower, Jupper)
                    PropUpRight = Field3DFather(Iupper, Jupper, k) / DYXFather(Iupper, Jupper)

                else

                    PropLowLeft = Field3DFather(Ilower, Jlower, k)
                    PropUpLeft  = Field3DFather(Iupper, Jlower, k)
                    PropLowRight= Field3DFather(Ilower, Jupper, k)
                    PropUpRight = Field3DFather(Iupper, Jupper, k)

                endif

                ONLowLeft   = ComputeFather(Ilower, Jlower, k)
                ONUpLeft    = ComputeFather(Iupper, Jlower, k)
                ONLowRight  = ComputeFather(Ilower, Jupper, k)
                ONUpRight   = ComputeFather(Iupper, Jupper, k)

                !griflet: No problem: all variables are private. Written: PropSon & InterpolOk only.
                call InterpolPoint(XPosition, YPosition,                             &
                                   YYUpper, YYLower, XXUpper, XXLower,               &
                                   PropLowLeft, PropUpLeft, PropLowRight, PropUpRight,&
                                   ONLowLeft,   ONUpLeft,   ONLowRight,   ONUpRight, &
                                   PropSon, InterPolOK)

                if (.not. InterPolOK) then

                    PropSon      = FillValueReal

                endif

                if (KUBFather > 1) then

                    Field3DSon(i, j, k) = PropSon


                else

                    Field3DSon(i, j, KUBSon) = PropSon

                endif



            enddo doi
            enddo doj
            !! $OMP END DO NOWAIT
            enddo dok
            !! $OMP END PARALLEL

            nullify(XXSon,    YYSon   )
            nullify(JSon,     ISon    )
            nullify(XXFather, YYFather)
            nullify(DYXFather         )


            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        if (MonitorPerformance) call StopWatch ("ModuleHorizontalGrid", "InterpolRegularGrid3D")

        !----------------------------------------------------------------------

    end subroutine InterpolRegularGrid3D


    !--------------------------------------------------------------------------

    subroutine InterpolRegularGrid3D8(HorizontalGridSonID, HorizontalGridFatherID,       &
                                      Field3DFather, Field3DSon,                         &
                                      ComputeFather, Compute,                            &
                                      KLBFather, KUBFather, KUBSon,                      &
                                      FluxType, GridID, STAT)



        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridSonID
        integer                                     :: HorizontalGridFatherID
        real(8), dimension(:,:,:), pointer          :: Field3DFather, Field3DSon
        integer, dimension(:,:,:), pointer          :: ComputeFather
        integer,           intent(IN)               :: Compute, KLBFather, KUBFather, KUBSon
        logical, optional, intent(IN)               :: FluxType
        integer, optional, intent(IN)               :: GridID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        type(T_HorizontalGrid),    pointer          :: ObjHorizontalGridFather
        type (T_FatherGrid),       pointer          :: FatherGrid
        real,    dimension(:,:  ), pointer          :: DXFather, DYFather, XXFather2D, YYFather2D
        real   , pointer, dimension(:,:)            :: XXSon, YYSon, DYXFather
        real   , pointer, dimension(:  )            :: XXFather, YYFather
        integer, pointer, dimension(:,:)            :: JSon, ISon, DefinedPoint


        real    :: YYUpper, YYLower, XXUpper, XXLower,                                   &
                   PropLowLeft,  PropUpLeft, PropLowRight, PropUpRight, PropSon

        integer :: ONLowLeft, ONUpLeft, ONLowRight, ONUpRight
        integer :: ready_, GridID_

        integer :: JLBSon, JUBSon, ILBSon, IUBSon
        integer :: JLBFather, JUBFather, ILBFather, IUBFather
        integer :: Jlower, Jupper, Ilower, Iupper
        integer :: STAT_

        integer :: i, j, k
        logical :: InterPolOK, FluxType_

        real    :: XPosition, YPosition
        integer :: Dij, Dji
        !$ integer  :: CHUNK

        !Begin------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleHorizontalGrid", "InterpolRegularGrid3D8")

        STAT_ = UNKNOWN_


        call Ready(HorizontalGridSonID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then


            if (.NOT. associated(Field3DFather))                                         &
                call SetError(FATAL_, INTERNAL_,                                           &
                             "InterpolRegularGrid3D8; HorizontalGrid. ERR01")

            if (.NOT. associated(Field3DSon))                                            &
                call SetError(FATAL_, INTERNAL_,                                           &
                             "InterpolRegularGrid3D8; HorizontalGrid. ERR02")

            nullify(ObjHorizontalGridFather)
            call LocateObjFather (ObjHorizontalGridFather, HorizontalGridFatherID)


            if (present(FluxType)) then

                FluxType_ = FluxType

            else

                FluxType_ = .false.

            endif

            JLBSon = Me%WorkSize%JLB
            JUBSon = Me%WorkSize%JUB
            ILBSon = Me%WorkSize%ILB
            IUBSon = Me%WorkSize%IUB

            if (present(GridID)) then
                GridID_ = GridID
            else
                GridID_ = LargeScaleModel_
            endif


            call Search_FatherGrid (FatherGrid, GridID_)


            !Gets the bounds from the Father
            ILBFAther = ObjHorizontalGridFather%WorkSize%ILB
            IUBFather = ObjHorizontalGridFather%WorkSize%IUB
            JLBFather = ObjHorizontalGridFather%WorkSize%JLB
            JUBFather = ObjHorizontalGridFather%WorkSize%JUB


            if      (Compute == ComputeU_)  then

                JUBFather = JUBFather + 1
                JUBSon    = JUBSon    + 1


                XXSon    => FatherGrid%XX_U
                YYSon    => FatherGrid%YY_U

                JSon     => FatherGrid%JU
                ISon     => FatherGrid%IU


                XXFather => ObjHorizontalGridFather%Compute%XX_U
                YYFather => ObjHorizontalGridFather%Compute%YY_U

                DYXFather=> ObjHorizontalGridFather%DYY

                XXFather2D => ObjHorizontalGridFather%Compute%XX2D_U
                YYFather2D => ObjHorizontalGridFather%Compute%YY2D_U

                DXFather => ObjHorizontalGridFather%DXX
                DYFather => ObjHorizontalGridFather%DZY

                Dij = 0
                Dji = 1

                if (Me%CornersXYInput) DefinedPoint => Me%DefineFacesUMap

            else if (Compute == ComputeV_)  then

                IUBFather = IUBFather + 1
                IUBSon    = IUBSon    + 1

                XXSon => FatherGrid%XX_V
                YYSon => FatherGrid%YY_V

                JSon     => FatherGrid%JV
                ISon     => FatherGrid%IV

                XXFather => ObjHorizontalGridFather%Compute%XX_V
                YYFather => ObjHorizontalGridFather%Compute%YY_V

                DYXFather=> ObjHorizontalGridFather%DXX

                XXFather2D => ObjHorizontalGridFather%Compute%XX2D_V
                YYFather2D => ObjHorizontalGridFather%Compute%YY2D_V

                DXFather => ObjHorizontalGridFather%DZX
                DYFather => ObjHorizontalGridFather%DYY

                Dij = 1
                Dji = 0

                if (Me%CornersXYInput) DefinedPoint => Me%DefineFacesVMap


            else if (Compute == ComputeZ_)  then

                if (FluxType_) call SetError(FATAL_, INTERNAL_,                            &
                             "InterpolRegularGrid3D8; HorizontalGrid. ERR04")

                XXSon => FatherGrid%XX_Z
                YYSon => FatherGrid%YY_Z

                XXFather => ObjHorizontalGridFather%Compute%XX_Z
                YYFather => ObjHorizontalGridFather%Compute%YY_Z

                JSon     => FatherGrid%JZ
                ISon     => FatherGrid%IZ

                XXFather2D => ObjHorizontalGridFather%Compute%XX2D_Z
                YYFather2D => ObjHorizontalGridFather%Compute%YY2D_Z

                DXFather => ObjHorizontalGridFather%DVX
                DYFather => ObjHorizontalGridFather%DUY

                Dji = 1
                Dij = 1

                if (Me%CornersXYInput) DefinedPoint => Me%DefineCellsMap

            else

                write(*,*) 'Compute ==', Compute
                stop "InterpolRegularGrid3D8; HorizontalGrid. ERR10"

            endif


            !$ CHUNK = CHUNK_J(JLBSon,JUBSon)

            !griflet
            !! $OMP PARALLEL PRIVATE( k,j,i,                     &
            !! $OMP                   InterPolOK, Jlower, Jupper, &
            !! $OMP                   Ilower,Iupper,              &
            !! $OMP                   XXLower,XXUpper,            &
            !! $OMP                   YYLower,YYUpper,            &
            !! $OMP                   XPosition,YPosition,        &
            !! $OMP                   PropLowLeft,PropUpLeft,     &
            !! $OMP                   PropLowRight,PropUpRight,   &
            !! $OMP                   ONLowLeft,ONUpLeft,         &
            !! $OMP                   ONLowRight,ONUpRight,PropSon)
            !! $OMP           SHARED ( DXFather,DYFather,         &
            !! $OMP                   XXFather2D, YYFather2D,     &
            !! $OMP                   XXSon, YYSon)
dok:        do k = KLBFather, KUBFather
            !! $OMP DO SCHEDULE(DYNAMIC,CHUNK)
doj:        do j = JLBSon, JUBSon
doi:        do i = ILBSon, IUBSon

                if (Me%CornersXYInput) then
                    if (DefinedPoint(i,j) == 0) cycle
                endif

                InterPolOK    = .false.

                Jlower = JSon(i, j)

                !Father domain smaller than son domain
                if (Jlower < -100) then
                    Field3DSon(i, j, k) = FillValueReal
                    cycle
                endif

                Jupper = JSon(i, j) + 1

                Ilower = ISon(i, j)

                !Father domain smaller than son domain
                if (Jlower < -100) then
                    Field3DSon(i, j, k) = FillValueReal
                    cycle
                endif


                Iupper = ISon(i, j) + 1

                if (Jupper > JUBFather) Jupper = JUBFather
                if (Iupper > IUBFather) Iupper = IUBFather
                if (Jlower < JLBFather) Jlower = JLBFather
                if (Ilower < ILBFather) Ilower = ILBFather

                if (Jupper < JLBFather) Jupper = JLBFather
                if (Iupper < ILBFather) Iupper = ILBFather

                if (.not. ObjHorizontalGridFather%CornersXYInput) then


                    XXLower  = XXFather(Jlower)
                    XXUpper  = XXFather(Jupper)

                    YYLower  = YYFather(Ilower)
                    YYUpper  = YYFather(Iupper)

                    XPosition = XXSon(i,j)
                    YPosition = YYSon(i,j)

                else

                    call CellReferential(DXFather, DYFather, XXFather2D, YYFather2D,    &
                                         ObjHorizontalGridFather%RotationX,             &
                                         ObjHorizontalGridFather%RotationY,             &
                                         XXSon(i,j), YYSon(i,j),                        &
                                         Ilower, Iupper, Jlower, Jupper, dij, dji,      &
                                         XPosition, YPosition,                          &
                                         YYUpper, YYLower, XXUpper, XXLower)

                endif


                if (FluxType_) then

                    PropLowLeft = Field3DFather(Ilower, Jlower, k) / DYXFather(Ilower, Jlower)
                    PropUpLeft  = Field3DFather(Iupper, Jlower, k) / DYXFather(Iupper, Jlower)
                    PropLowRight= Field3DFather(Ilower, Jupper, k) / DYXFather(Ilower, Jupper)
                    PropUpRight = Field3DFather(Iupper, Jupper, k) / DYXFather(Iupper, Jupper)

                else

                    PropLowLeft = Field3DFather(Ilower, Jlower, k)
                    PropUpLeft  = Field3DFather(Iupper, Jlower, k)
                    PropLowRight= Field3DFather(Ilower, Jupper, k)
                    PropUpRight = Field3DFather(Iupper, Jupper, k)

                endif

                ONLowLeft   = ComputeFather(Ilower, Jlower, k)
                ONUpLeft    = ComputeFather(Iupper, Jlower, k)
                ONLowRight  = ComputeFather(Ilower, Jupper, k)
                ONUpRight   = ComputeFather(Iupper, Jupper, k)


                call InterpolPoint(XPosition, YPosition,                             &
                                   YYUpper, YYLower, XXUpper, XXLower,               &
                                   PropLowLeft, PropUpLeft, PropLowRight, PropUpRight,&
                                   ONLowLeft,   ONUpLeft,   ONLowRight,   ONUpRight, &
                                   PropSon, InterPolOK)


                if (.not. InterPolOK) then
                    PropSon     = FillValueReal
                endif


                if (KUBFather > 1) then
                    Field3DSon(i, j, k) =  dble(PropSon)
                else
                    Field3DSon(i, j, KUBSon) =  dble(PropSon)
                endif


            enddo doi
            enddo doj
            !! $OMP END DO NOWAIT
            enddo dok
            !! $OMP END PARALLEL


            nullify(XXSon,    YYSon   )
            nullify(JSon,     ISon    )
            nullify(XXFather, YYFather)
            nullify(DYXFather         )

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        if (MonitorPerformance) call StopWatch ("ModuleHorizontalGrid", "InterpolRegularGrid3D8")

        !----------------------------------------------------------------------

    end subroutine InterpolRegularGrid3D8



    !--------------------------------------------------------------------------
    subroutine CellReferential(DXFather, DYFather, XXFather2D, YYFather2D,    &
                               RotationX, RotationY,                          &
                               XXSon, YYSon,                                  &
                               Ilower, Iupper, Jlower, Jupper, dij, dji,      &
                               XPosition, YPosition,                          &
                               YYUpper, YYLower, XXUpper, XXLower)

        !Arguments-------------------------------------------------------------
        real   , pointer, dimension(:,:), intent(IN)  :: DXFather, DYFather,            &
                                                         XXFather2D, YYFather2D,        &
                                                         RotationX, RotationY
        real,                             intent(IN)  :: XXSon, YYSon
        integer,                          intent(IN)  :: Jlower, Jupper, Ilower, Iupper, dij, dji
        real,                             intent(OUT) :: YYUpper, YYLower, XXUpper, XXLower, &
                                                         XPosition, YPosition
        !Local-----------------------------------------------------------------
        real                              :: XlocalOrigin, YlocalOrigin, dx, dy, dteta


        !Begin------------------------------------------------------------------

            XXLower  = 0.
            XXUpper  = DXFather(Ilower + dji, Jlower)

            YYLower  = 0.
            YYUpper  = DYFather(Ilower, Jlower + dij)

            !West-East
            XlocalOrigin = (XXFather2D(Ilower, Jlower) + XXFather2D(Iupper, Jlower)) / 2.
            YlocalOrigin = (YYFather2D(Ilower, Jlower) + YYFather2D(Iupper, Jlower)) / 2.
            dx           = XXSon - XlocalOrigin
            dy           = YYSon - YlocalOrigin
            dteta        = atan2(dy, dx)
            dteta        = dteta - (RotationX(Ilower, Jlower) + RotationX(Ilower + dji, Jlower + dij)) / 2.
            XPosition    = sqrt(dx**2 + dy**2)
            XPosition    = XPosition * cos(dteta)
            XPosition    = abs(XPosition)

            !South-North
            XlocalOrigin = (XXFather2D(Ilower, Jlower) + XXFather2D(Ilower, Jupper)) / 2.
            YlocalOrigin = (YYFather2D(Ilower, Jlower) + YYFather2D(Ilower, Jupper)) / 2.
            dx           = XXSon - XlocalOrigin
            dy           = YYSon - YlocalOrigin
            dteta        = atan2(dy, dx)
            dteta        = dteta - (RotationY(Ilower, Jlower) + RotationY(Ilower + dji, Jlower + dij)) / 2.
            YPosition    = sqrt(dx**2 + dy**2)
            YPosition    = YPosition * cos(dteta)
            YPosition    = abs(YPosition)

    end subroutine CellReferential

    !----------------------------------------------------------------------------


    subroutine WriteHorizontalGrid (HorizontalGridID, ObjHDF5, OutputNumber, WorkSize,  &
                                    WindowGrid, GlobalWorkSizeWindow, STAT)


        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer                                     :: ObjHDF5
        integer,        optional                    :: OutputNumber
        type(T_Size2D), optional                    :: WorkSize
        logical       , optional                    :: WindowGrid
        type(T_Size2D), optional                    :: GlobalWorkSizeWindow
        integer,        optional                    :: STAT

        !Local-----------------------------------------------------------------
#if _GOOGLEMAPS
        real(8),     dimension(:,:), pointer        :: XX_aux, YY_aux
#endif
        integer,     dimension(:  ), pointer        :: AuxInt4
        type(T_Size2D)                              :: WorkSize_
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: STAT_, ready_, STAT_CALL, ilen, iFile, i
        character(len=PathLength)                   :: FileName, AuxFile
        character(len=StringLength)                 :: AuxChar
        logical                                     :: WindowGrid_
        type(T_Size2D)                              :: GlobalWorkSizeWindow_

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then

            !WorkSize
            if (present(WorkSize)) then
                WorkILB = WorkSize%ILB
                WorkIUB = WorkSize%IUB

                WorkJLB = WorkSize%JLB
                WorkJUB = WorkSize%JUB

                WorkSize_ = WorkSize
            else
                WorkILB = Me%WorkSize%ILB
                WorkIUB = Me%WorkSize%IUB

                WorkJLB = Me%WorkSize%JLB
                WorkJUB = Me%WorkSize%JUB

                WorkSize_ = Me%WorkSize
            endif

            if (present(WindowGrid)) then
                WindowGrid_ = WindowGrid
            else
                WindowGrid_ = .false.
            endif

            if (present(GlobalWorkSizeWindow)) then
                GlobalWorkSizeWindow_ = GlobalWorkSizeWindow
            endif

            !Sets limits for next write operations
            call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB+1, WorkJLB, WorkJUB+1,      &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR02'

            if (present(OutputNumber)) then

                call HDF5WriteData   (ObjHDF5, "/Grid/ConnectionX", "ConnectionX", "m", &
                                      Array2D = Me%XX_IE,                               &
                                      OutputNumber = OutputNumber,                      &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR03'

                call HDF5WriteData   (ObjHDF5, "/Grid/ConnectionY", "ConnectionY", "m", &
                                      Array2D = Me%YY_IE,                               &
                                      OutputNumber = OutputNumber,                      &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR04'

                call HDF5WriteData   (ObjHDF5, "/Grid/Longitude", "Longitude", "º",     &
                                      Array2D = Me%LongitudeConn,                       &
                                      OutputNumber = OutputNumber,                      &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR05'

                call HDF5WriteData   (ObjHDF5, "/Grid/Latitude", "Latitude", "º",       &
                                      Array2D = Me%LatitudeConn,                        &
                                      OutputNumber = OutputNumber,                      &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR06'

                if (Me%CornersXYInput) then

                    call HDF5WriteData   (ObjHDF5, "/Grid/Define_Cells", "Define Cells", "-",&
                                          Array2D = Me%DefineCellsMap,                  &
                                          OutputNumber = OutputNumber,                  &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR09'


                endif

            else

                call HDF5WriteData   (ObjHDF5, "/Grid", "ConnectionX", "m",             &
                                      Array2D = Me%XX_IE,                               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR03'

                call HDF5WriteData   (ObjHDF5, "/Grid", "ConnectionY", "m",             &
                                      Array2D = Me%YY_IE,                               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR04'

                call HDF5WriteData   (ObjHDF5, "/Grid", "Longitude", "º",               &
                                      Array2D = Me%LongitudeConn,                       &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR05'

                call HDF5WriteData   (ObjHDF5, "/Grid", "Latitude", "º",                &
                                      Array2D = Me%LatitudeConn,                        &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR06'

                if (Me%CornersXYInput) then

                    call HDF5WriteData   (ObjHDF5, "/Grid", "Define Cells", "-",        &
                                          Array2D = Me%DefineCellsMap,                  &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR09'


                endif

#if _GOOGLEMAPS

                !Geographic coordinates
                if (Me%CoordType == GEOG_               .or.  &
                    Me%CoordType == SIMPLE_GEOG_) then

                        allocate(XX_aux(WorkILB-1:WorkIUB+1, WorkJLB-1:WorkJUB+1))
                        allocate(YY_Aux(WorkILB-1:WorkIUB+1, WorkJLB-1:WorkJUB+1))

                        call WGS84toGoogleMaps(Me%LongitudeConn, Me%LatitudeConn,       &
                                               WorkSize_%ILB, WorkSize_%IUB,            &
                                               WorkSize_%JLB, WorkSize_%JUB,            &
                                               XX_aux, YY_aux)

                        !Sets limits for next write operations
                        call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB+1, WorkJLB, WorkJUB+1,      &
                                              STAT = STAT_CALL)

                        call HDF5WriteData   (ObjHDF5, "/Grid", "googlemaps_x", "-",    &
                                              Array2D = XX_aux,                         &
                                              STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR20'

                        call HDF5WriteData   (ObjHDF5, "/Grid", "googlemaps_y", "-",&
                                              Array2D = YY_aux,                         &
                                              STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR30'

                        deallocate(XX_aux,YY_Aux)

                endif

#endif

                if (Me%DDecomp%MasterOrSlave) then

                    allocate(AuxInt4(4))

                    !Sets limits for next write operations
                    call HDF5SetLimits   (ObjHDF5, 1, 4, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR40'


                    if (WindowGrid_) then
                    !For the case of the window hdf5 output
                        AuxInt4(1) = GlobalWorkSizeWindow_%ILB
                        AuxInt4(2) = GlobalWorkSizeWindow_%IUB
                        AuxInt4(3) = GlobalWorkSizeWindow_%JLB
                        AuxInt4(4) = GlobalWorkSizeWindow_%JUB
                    else
                        AuxInt4(1) = Me%DDecomp%Global%ILB
                        AuxInt4(2) = Me%DDecomp%Global%IUB
                        AuxInt4(3) = Me%DDecomp%Global%JLB
                        AuxInt4(4) = Me%DDecomp%Global%JUB
                    endif

                    call HDF5WriteData   (ObjHDF5, "/Grid/Decomposition/Global", "ILB_IUB_JLB_JUB", &
                                          "-", Array1D = AuxInt4, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR50'

                    AuxInt4(1) = Me%DDecomp%HaloMap%ILB + WorkILB - 1
                    AuxInt4(2) = Me%DDecomp%HaloMap%ILB + WorkIUB - 1
                    AuxInt4(3) = Me%DDecomp%HaloMap%JLB + WorkJLB - 1
                    AuxInt4(4) = Me%DDecomp%HaloMap%JLB + WorkJUB - 1

                    call HDF5WriteData   (ObjHDF5, "/Grid/Decomposition/Mapping", "ILB_IUB_JLB_JUB",&
                                          "-", Array1D = AuxInt4, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR60'

                    if (AuxInt4(1) == Me%DDecomp%HaloMap%ILB) then
                        AuxInt4(1) =  Me%DDecomp%Mapping%ILB
                    endif

                    if (AuxInt4(2) == Me%DDecomp%HaloMap%IUB) then
                        AuxInt4(2) =  Me%DDecomp%Mapping%IUB
                    endif

                    if (AuxInt4(3) == Me%DDecomp%HaloMap%JLB) then
                        AuxInt4(3) =  Me%DDecomp%Mapping%JLB
                    endif

                    if (AuxInt4(4) == Me%DDecomp%HaloMap%JUB) then
                        AuxInt4(4) =  Me%DDecomp%Mapping%JUB
                    endif

                    call HDF5WriteData   (ObjHDF5, "/Grid/Decomposition/InnerMapping", "ILB_IUB_JLB_JUB",&
                                          "-", Array1D = AuxInt4, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR70'

                    deallocate(AuxInt4)

                    if (.not. Me%DDecomp%FilesListOpen) then

                        call UnitsManager(Me%DDecomp%FilesListID, FileOpen, STAT = STAT_CALL)

                        if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR80'

                        write(AuxChar,fmt='(i5)') Me%DDecomp%MPI_ID

                        Me%DDecomp%FilesListName = trim(adjustl(Me%DDecomp%FilesListName))
                        Me%DDecomp%ModelPath     = trim(adjustl(Me%DDecomp%ModelPath    ))

                        Me%DDecomp%FilesListName = "MPI_"//trim(adjustl(Auxchar))//"_"//trim(Me%DDecomp%FilesListName)

                        ilen = len_trim(Me%DDecomp%ModelPath)

                        Me%DDecomp%ModelPath(ilen-2: ilen) = "res"

                        !windows path
                        AuxFile = trim(adjustl(Me%DDecomp%ModelPath))//"/"//trim(adjustl(Me%DDecomp%FilesListName))

                        open(file   = AuxFile,                                          &
                             unit   = Me%DDecomp%FilesListID,                           &
                             status = "unknown",                                        &
                             form   = "formatted",                                      &
                             IOSTAT = STAT_CALL)

                        if (STAT_CALL /= SUCCESS_) then
                            !linux path
                            AuxFile = trim(adjustl(Me%DDecomp%ModelPath))//backslash//trim(adjustl(Me%DDecomp%FilesListName))

                            open(file   = AuxFile,                                      &
                                 unit   = Me%DDecomp%FilesListID,                       &
                                 status = "unknown",                                    &
                                 form   = "formatted",                                  &
                                 IOSTAT = STAT_CALL)

                            if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR90'
                        endif

                        Me%DDecomp%FilesListOpen = .true.

                    endif

                    if (Me%DDecomp%FilesListOpen) then

                        call GetHDF5FileName (ObjHDF5, FileName, STAT= STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid - HorizontalGrid - ERR100'
                        iFile = 1
                        ilen  = len_trim(FileName)
                        do i = ilen,1,-1
                            if (FileName(i:i) == '/' .or. FileName(i:i) == backslash) then
                                iFile = i+1
                                exit
                            endif
                        enddo

                        write(Me%DDecomp%FilesListID,'(A)') trim(FileName(iFile:ilen))

                    endif

                endif

            endif

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT))  STAT = STAT_


    end subroutine WriteHorizontalGrid

    !--------------------------------------------------------------------------


    subroutine WriteHorizontalGrid_UV (HorizontalGridID, ObjHDF5, WorkSize, STAT)


        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer                                     :: ObjHDF5
        type(T_Size2D)                              :: WorkSize
        integer,        optional                    :: STAT

        !Local-----------------------------------------------------------------
        real(8),     dimension(:,:), pointer        :: Aux2D
        integer,     dimension(:  ), pointer        :: AuxInt4
        type(T_Size2D)                              :: WorkSize_
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: STAT_, ready_, STAT_CALL, ilen, iFile, i
        character(len=PathLength)                   :: FileName, AuxFile
        character(len=StringLength)                 :: AuxChar

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then

            !WorkSize
            WorkILB = WorkSize%ILB
            WorkIUB = WorkSize%IUB

            WorkJLB = WorkSize%JLB
            WorkJUB = WorkSize%JUB

            WorkSize_ = WorkSize

            !Sets limits for next write operations
            call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB+2, WorkJLB, WorkJUB+2,      &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR10'

            allocate(Aux2D(WorkILB:WorkIUB+2, WorkJLB:WorkJUB+2))

            Aux2D(WorkILB:WorkIUB+1, WorkJLB:WorkJUB+1) = Me%XX_IE(WorkILB:WorkIUB+1, WorkJLB:WorkJUB+1)
            Aux2D(WorkIUB+2        , WorkJLB:WorkJUB+1) = Me%XX_IE(WorkIUB+1        , WorkJLB:WorkJUB+1)
            Aux2D(WorkILB:WorkIUB+2,         WorkJUB+2) = Aux2D   (WorkILB:WorkIUB+2,         WorkJUB+1)

            call HDF5WriteData   (ObjHDF5, "/Grid", "ConnectionX", "m",             &
                                  Array2D = Aux2D,                                  &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR20'

            Aux2D(WorkILB:WorkIUB+1, WorkJLB:WorkJUB+1) = Me%YY_IE(WorkILB:WorkIUB+1, WorkJLB:WorkJUB+1)
            Aux2D(WorkIUB+2        , WorkJLB:WorkJUB+1) = Me%YY_IE(WorkIUB+1        , WorkJLB:WorkJUB+1)
            Aux2D(WorkILB:WorkIUB+2,         WorkJUB+2) = Aux2D   (WorkILB:WorkIUB+2,         WorkJUB+1)

            call HDF5WriteData   (ObjHDF5, "/Grid", "ConnectionY", "m",             &
                                  Array2D = Aux2D,                                  &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR30'

            Aux2D(WorkILB:WorkIUB+1, WorkJLB:WorkJUB+1) = Me%LongitudeConn(WorkILB:WorkIUB+1, WorkJLB:WorkJUB+1)
            Aux2D(WorkIUB+2        , WorkJLB:WorkJUB+1) = Me%LongitudeConn(WorkIUB+1        , WorkJLB:WorkJUB+1)
            Aux2D(WorkILB:WorkIUB+2,         WorkJUB+2) = Aux2D           (WorkILB:WorkIUB+2,         WorkJUB+1)

            call HDF5WriteData   (ObjHDF5, "/Grid", "Longitude", "º",               &
                                  Array2D = Aux2D,                                  &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR40'

            Aux2D(WorkILB:WorkIUB+1, WorkJLB:WorkJUB+1) = Me%LatitudeConn(WorkILB:WorkIUB+1, WorkJLB:WorkJUB+1)
            Aux2D(WorkIUB+2        , WorkJLB:WorkJUB+1) = Me%LatitudeConn(WorkIUB+1        , WorkJLB:WorkJUB+1)
            Aux2D(WorkILB:WorkIUB+2,         WorkJUB+2) = Aux2D          (WorkILB:WorkIUB+2,         WorkJUB+1)

            call HDF5WriteData   (ObjHDF5, "/Grid", "Latitude", "º",                &
                                  Array2D = Aux2D,                                  &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR50'

            deallocate(Aux2D)

            if (Me%CornersXYInput) then

                call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB+1, WorkJLB, WorkJUB+1,      &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR60'

                call HDF5WriteData   (ObjHDF5, "/Grid", "Define Cells", "-",        &
                                      Array2D = Me%DefineCellsMap,                  &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR70'


            endif

            if (Me%DDecomp%MasterOrSlave) then

                allocate(AuxInt4(4))

                !Sets limits for next write operations
                call HDF5SetLimits   (ObjHDF5, 1, 4, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR80'

                AuxInt4(1) = Me%DDecomp%Global%ILB
                AuxInt4(2) = Me%DDecomp%Global%IUB + 1
                AuxInt4(3) = Me%DDecomp%Global%JLB
                AuxInt4(4) = Me%DDecomp%Global%JUB + 1

                call HDF5WriteData   (ObjHDF5, "/Grid/Decomposition/Global", "ILB_IUB_JLB_JUB", &
                                      "-", Array1D = AuxInt4, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR90'

                AuxInt4(1) = Me%DDecomp%HaloMap%ILB + WorkILB - 1
                AuxInt4(2) = Me%DDecomp%HaloMap%ILB + WorkIUB - 1 + 1
                AuxInt4(3) = Me%DDecomp%HaloMap%JLB + WorkJLB - 1
                AuxInt4(4) = Me%DDecomp%HaloMap%JLB + WorkJUB - 1 + 1

                call HDF5WriteData   (ObjHDF5, "/Grid/Decomposition/Mapping", "ILB_IUB_JLB_JUB",&
                                      "-", Array1D = AuxInt4, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR100'

                if (AuxInt4(1) == Me%DDecomp%HaloMap%ILB) then
                    AuxInt4(1) =  Me%DDecomp%Mapping%ILB
                endif

                if (AuxInt4(2) == Me%DDecomp%HaloMap%IUB) then
                    AuxInt4(2) =  Me%DDecomp%Mapping%IUB + 1
                endif

                if (AuxInt4(3) == Me%DDecomp%HaloMap%JLB) then
                    AuxInt4(3) =  Me%DDecomp%Mapping%JLB
                endif

                if (AuxInt4(4) == Me%DDecomp%HaloMap%JUB) then
                    AuxInt4(4) =  Me%DDecomp%Mapping%JUB + 1
                endif

                call HDF5WriteData   (ObjHDF5, "/Grid/Decomposition/InnerMapping", "ILB_IUB_JLB_JUB",&
                                      "-", Array1D = AuxInt4, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR110'

                deallocate(AuxInt4)

                if (.not. Me%DDecomp%FilesListOpen) then

                    call UnitsManager(Me%DDecomp%FilesListID, FileOpen, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR120'

                    write(AuxChar,fmt='(i5)') Me%DDecomp%MPI_ID

                    Me%DDecomp%FilesListName = trim(adjustl(Me%DDecomp%FilesListName))
                    Me%DDecomp%ModelPath     = trim(adjustl(Me%DDecomp%ModelPath    ))

                    Me%DDecomp%FilesListName = "MPI_"//trim(adjustl(Auxchar))//"_"//trim(Me%DDecomp%FilesListName)

                    ilen = len_trim(Me%DDecomp%ModelPath)

                    Me%DDecomp%ModelPath(ilen-2: ilen) = "res"

                    !windows path
                    AuxFile = trim(adjustl(Me%DDecomp%ModelPath))//"/"//trim(adjustl(Me%DDecomp%FilesListName))

                    open(file   = AuxFile,                                          &
                         unit   = Me%DDecomp%FilesListID,               &
                         status = "unknown",                                        &
                         form   = "formatted",                                      &
                         IOSTAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) then
                        !linux path
                        AuxFile = trim(adjustl(Me%DDecomp%ModelPath))//backslash//trim(adjustl(Me%DDecomp%FilesListName))

                        open(file   = AuxFile,                                      &
                             unit   = Me%DDecomp%FilesListID,           &
                             status = "unknown",                                    &
                             form   = "formatted",                                  &
                             IOSTAT = STAT_CALL)

                        if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR130'
                    endif

                    Me%DDecomp%FilesListOpen = .true.

                endif

                if (Me%DDecomp%FilesListOpen) then

                    call GetHDF5FileName (ObjHDF5, FileName, STAT= STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteHorizontalGrid_UV - HorizontalGrid - ERR140'
                    iFile = 1
                    ilen  = len_trim(FileName)
                    do i = ilen,1,-1
                        if (FileName(i:i) == '/' .or. FileName(i:i) == backslash) then
                            iFile = i+1
                            exit
                        endif
                    enddo

                    write(Me%DDecomp%FilesListID,'(A)') trim(FileName(iFile:ilen))

                endif

            endif

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT))  STAT = STAT_


    end subroutine WriteHorizontalGrid_UV

    !--------------------------------------------------------------------------

    subroutine ReadHDF5HorizontalGrid ()


        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: STAT_CALL
        logical                                     :: Exist
        !Begin-----------------------------------------------------------------


        !WorkSize
        WorkILB = Me%WorkSize%ILB
        WorkIUB = Me%WorkSize%IUB
        WorkJLB = Me%WorkSize%JLB
        WorkJUB = Me%WorkSize%JUB


        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB+1, WorkJLB, WorkJUB+1,      &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5HorizontalGrid - HorizontalGrid - ERR10'

        if (Me%CoordType == SIMPLE_GEOG_) then

            call HDF5ReadData   (Me%ObjHDF5, "/Grid", "Longitude",                          &
                                  Array2D = Me%LongitudeConn,                               &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5HorizontalGrid - HorizontalGrid - ERR20'

            call HDF5ReadData   (Me%ObjHDF5, "/Grid", "Latitude",                           &
                                  Array2D = Me%LatitudeConn,                                &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5HorizontalGrid - HorizontalGrid - ERR30'

        else

            call HDF5ReadData   (Me%ObjHDF5, "/Grid", "ConnectionX",                        &
                                  Array2D = Me%XX_IE,                                       &
                                  STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5HorizontalGrid - HorizontalGrid - ERR40'

            call HDF5ReadData   (Me%ObjHDF5, "/Grid", "ConnectionY",                        &
                                  Array2D = Me%YY_IE,                                       &
                                  STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5HorizontalGrid - HorizontalGrid - ERR50'

        endif


        if (Me%CornersXYInput) then

            call GetHDF5DataSetExist (Me%ObjHDF5, DataSetName ="/Grid/Define Cells",    &
                                      Exist = Exist, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5HorizontalGrid - HorizontalGrid - ERR60'

            if (Exist) then

                call HDF5ReadData   (Me%ObjHDF5, "/Grid", "Define Cells",               &
                                      Array2D = Me%DefineCellsMap,                      &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5HorizontalGrid - HorizontalGrid - ERR70'

            endif

        endif


    end subroutine ReadHDF5HorizontalGrid

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------


    subroutine RotateVectorFieldToGrid2D(HorizontalGridID, VectorInX, VectorInY,        &
                                         VectorOutX, VectorOutY, WaterPoints2D,         &
                                         RotateX, RotateY, STAT)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, :), pointer       :: VectorInX, VectorInY, VectorOutX, VectorOutY
        integer, dimension(:, :), pointer       :: WaterPoints2D
        integer                                 :: HorizontalGridID
        logical                                 :: RotateX, RotateY
        integer, optional                       :: STAT
        !Local-----------------------------------------------------------------
        real                                    :: XGrid, YGrid, AngleX, AngleY
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j
        integer                                 :: STAT_, ready_
        real                                    :: GridRotationRadians

        !Begin--------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then



if1:        if  (Me%Distortion .or. Me%RegularRotation) then

                !Bounds
                ILB = Me%WorkSize%ILB
                IUB = Me%WorkSize%IUB

                JLB = Me%WorkSize%JLB
                JUB = Me%WorkSize%JUB

                GridRotationRadians = Me%Grid_Angle * Pi/180.

                !Rotate
                do j = JLB, JUB
                do i = ILB, IUB

if2:                if (WaterPoints2D(i, j) == WaterPoint) then

                        !VectorX = X_Grid * cos(AngleX) + Y_Grid * cos(AngleY)
                        !VectorY = X_Grid * sin(AngleX) + Y_Grid * sin(AngleY)

                        AngleX = 0.
                        AngleY = Pi / 2.

                        if      (Me%Distortion) then

                            AngleX = Me%RotationX(i, j)
                            AngleY = Me%RotationY(i, j)

                        else if (Me%RegularRotation) then

                            AngleX = GridRotationRadians
                            AngleY = GridRotationRadians + Pi / 2.

                        endif

                        call FromCartesianToGrid (VectorInX(i, j), VectorInY(i, j),           &
                                                  AngleX, AngleY, Xgrid, Ygrid)

                        if (RotateX) then

                            VectorOutX(i, j) = Xgrid

                        endif

                        if (RotateY) then

                            VectorOutY(i, j) = Ygrid

                        endif


                    endif if2

                enddo
                enddo

            else  if1

                if (RotateX) then

                    VectorOutX(:, :) = VectorInX(:, :)

                endif

                if (RotateY) then

                    VectorOutY(:, :) = VectorInY(:, :)

                endif

            endif if1

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT))  STAT = STAT_

    end subroutine RotateVectorFieldToGrid2D
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------


    subroutine RotateVectorFieldToGrid3D(HorizontalGridID, VectorInX, VectorInY,        &
                                         VectorOutX, VectorOutY, WaterPoints3D,         &
                                         RotateX, RotateY, KLB, KUB, STAT)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, :, :), pointer    :: VectorInX, VectorInY, VectorOutX, VectorOutY
        integer, dimension(:, :, :), pointer    :: WaterPoints3D
        logical                                 :: RotateX, RotateY
        integer                                 :: HorizontalGridID, KLB, KUB
        integer, optional                       :: STAT
        !Local-----------------------------------------------------------------
        real                                    :: XGrid, YGrid, AngleX, AngleY
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j, k
        integer                                 :: STAT_, ready_
        real                                    :: GridRotationRadians

        !Begin--------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then



if1:        if  (Me%Distortion .or. Me%RegularRotation) then

                !Bounds
                ILB = Me%WorkSize%ILB
                IUB = Me%WorkSize%IUB

                JLB = Me%WorkSize%JLB
                JUB = Me%WorkSize%JUB

                GridRotationRadians = Me%Grid_Angle * Pi/180.

                !Rotate
                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB

if2:                if (WaterPoints3D(i, j, k) == WaterPoint) then

                        !VectorX = X_Grid * cos(AngleX) + Y_Grid * cos(AngleY)
                        !VectorY = X_Grid * sin(AngleX) + Y_Grid * sin(AngleY)

                        AngleX = 0.
                        AngleY = Pi / 2.

                        if      (Me%Distortion) then

                            AngleX = Me%RotationX(i, j)
                            AngleY = Me%RotationY(i, j)

                        else if (Me%RegularRotation) then

                            AngleX = GridRotationRadians
                            AngleY = GridRotationRadians + Pi / 2.

                        endif

                        call FromCartesianToGrid (VectorInX(i, j, k), VectorInY(i, j, k),           &
                                                  AngleX, AngleY, Xgrid, Ygrid)

                        if (RotateX) then

                            VectorOutX(i, j, k) = Xgrid

                        endif

                        if (RotateY) then

                            VectorOutY(i, j, k) = Ygrid

                        endif


                    endif if2

                enddo
                enddo
                enddo

            else  if1


                if (RotateX) then

                    VectorOutX(:, :, :) = VectorInX(:, :, :)

                endif

                if (RotateY) then

                    VectorOutY(:, :, :) = VectorInY(:, :, :)

                endif

            endif if1

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT))  STAT = STAT_

    end subroutine RotateVectorFieldToGrid3D
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine RotateAngleFieldToGrid2D(HorizontalGridID, AngleIn, InReferential,       &
                                         AngleOut, WaterPoints2D,                      &
                                         Rotate, STAT)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, :), pointer       :: AngleIn, AngleOut
        integer, dimension(:, :), pointer       :: WaterPoints2D
        integer                                 :: HorizontalGridID
        logical                                 :: Rotate
        integer, optional                       :: STAT
        integer                                 :: InReferential
        !Local-----------------------------------------------------------------
        real                                    :: AngleOutCell, GridAngle
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j
        integer                                 :: STAT_, ready_

        !Begin--------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then



!if1:        if  (Me%Distortion .or. Me%RegularRotation) then

                !Bounds
                ILB = Me%WorkSize%ILB
                IUB = Me%WorkSize%IUB

                JLB = Me%WorkSize%JLB
                JUB = Me%WorkSize%JUB

                !Rotate
                do j = JLB, JUB
                do i = ILB, IUB

if2:                if (WaterPoints2D(i, j) == WaterPoint) then

                        !VectorX = X_Grid * cos(AngleX) + Y_Grid * cos(AngleY)
                        !VectorY = X_Grid * sin(AngleX) + Y_Grid * sin(AngleY)

                        GridAngle = 0.

                        if      (Me%Distortion) then

                            !with distortion, the cell trigonometric circle origin is coincident with cell "x plane"
                            !and rotation Y is not accounted
                            GridAngle = Me%RotationX(i, j) / 180. * Pi

                        else if (Me%RegularRotation) then

                            GridAngle = Me%Grid_Angle

                        endif

                        call AngleFromFieldToGrid (AngleIn(i, j), InReferential,           &
                                                         GridAngle, AngleOutCell)

                        if (Rotate) then

                            AngleOut(i, j) = AngleOutCell

                        endif


                    endif if2

                enddo
                enddo

            !else  if1
            !
            !    if (Rotate) then
            !
            !        AngleOut(:, :) = AngleIn(:, :)
            !
            !    endif
            !
            !
            !endif if1

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT))  STAT = STAT_

    end subroutine RotateAngleFieldToGrid2D

    !--------------------------------------------------------------------------

    subroutine RotateAngleFieldToGrid3D(HorizontalGridID, AngleIn, InReferential,       &
                                         AngleOut, WaterPoints3D,         &
                                         Rotate, KLB, KUB, STAT)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, : , :), pointer   :: AngleIn, AngleOut
        integer, dimension(:, : , :), pointer   :: WaterPoints3D
        integer                                 :: HorizontalGridID, KLB, KUB
        logical                                 :: Rotate
        integer, optional                       :: STAT
        integer                                 :: InReferential
        !Local-----------------------------------------------------------------
        real                                    :: AngleOutCell, GridAngle
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j, k
        integer                                 :: STAT_, ready_

        !Begin--------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then



!if1:        if  (Me%Distortion .or. Me%RegularRotation) then

                !Bounds
                ILB = Me%WorkSize%ILB
                IUB = Me%WorkSize%IUB

                JLB = Me%WorkSize%JLB
                JUB = Me%WorkSize%JUB

                !Rotate
                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB

if2:                if (WaterPoints3D(i, j, k) == WaterPoint) then

                        !VectorX = X_Grid * cos(AngleX) + Y_Grid * cos(AngleY)
                        !VectorY = X_Grid * sin(AngleX) + Y_Grid * sin(AngleY)

                        GridAngle = 0.

                        if      (Me%Distortion) then

                            !with distortion, the cell trigonometric circle origin is coincident with cell "x plane"
                            !and rotation Y is not accounted
                            GridAngle = Me%RotationX(i, j) / 180. * Pi

                        else if (Me%RegularRotation) then

                            GridAngle = Me%Grid_Angle

                        endif

                        call AngleFromFieldToGrid (AngleIn(i, j, k), InReferential,           &
                                                         GridAngle, AngleOutCell)

                        if (Rotate) then

                            AngleOut(i, j, k) = AngleOutCell

                        endif


                    endif if2

                enddo
                enddo
                enddo

            !else  if1
            !
            !    if (Rotate) then
            !
            !        AngleOut(:, :, :) = AngleIn(:, : , :)
            !
            !    endif
            !
            !
            !endif if1

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT))  STAT = STAT_

    end subroutine RotateAngleFieldToGrid3D

    !--------------------------------------------------------------------------

    subroutine RotateVectorGridToField2DR4(HorizontalGridID, VectorInX, VectorInY,      &
                                           VectorOutX, VectorOutY, WaterPoints2D,       &
                                           RotateX, RotateY, STAT)

        !Arguments-------------------------------------------------------------
        real(4), dimension(:, :), pointer       :: VectorInX, VectorInY, VectorOutX, VectorOutY
        integer, dimension(:, :), pointer       :: WaterPoints2D
        integer                                 :: HorizontalGridID
        logical                                 :: RotateX, RotateY
        integer, optional                       :: STAT
        !Local-----------------------------------------------------------------
        real(4)                                 :: XGrid, YGrid, AngleX, AngleY
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j
        integer                                 :: STAT_, ready_
        real                                    :: GridRotationRadians

        !Begin--------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then



if1:        if  (Me%Distortion .or. Me%RegularRotation) then

                !Bounds
                ILB = Me%WorkSize%ILB
                IUB = Me%WorkSize%IUB

                JLB = Me%WorkSize%JLB
                JUB = Me%WorkSize%JUB

                GridRotationRadians = Me%Grid_Angle * Pi / 180.

                !Rotate
                do j = JLB, JUB
                do i = ILB, IUB

if2:                if (WaterPoints2D(i, j) == WaterPoint) then

                        !VectorX = X_Grid * cos(AngleX) + Y_Grid * cos(AngleY)
                        !VectorY = X_Grid * sin(AngleX) + Y_Grid * sin(AngleY)

                        AngleX = 0.
                        AngleY = Pi / 2.

                        if      (Me%Distortion) then

                            AngleX = Me%RotationX(i, j)
                            AngleY = Me%RotationY(i, j)

                        else if (Me%RegularRotation) then

                            AngleX = GridRotationRadians
                            AngleY = GridRotationRadians + Pi / 2.

                        endif

                        call FromGridToCartesian (VectorInX(i, j), VectorInY(i, j),           &
                                                  AngleX, AngleY, Xgrid, Ygrid)

                        if (RotateX) then

                            VectorOutX(i, j) = Xgrid

                        endif

                        if (RotateY) then

                            VectorOutY(i, j) = Ygrid

                        endif


                    endif if2

                enddo
                enddo

            else  if1


                if (RotateX) then

                    VectorOutX(:, :) = VectorInX(:, :)

                endif

                if (RotateY) then

                    VectorOutY(:, :) = VectorInY(:, :)

                endif


            endif if1

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT))  STAT = STAT_

    end subroutine RotateVectorGridToField2DR4
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine RotateVectorGridToField2DR8(HorizontalGridID, VectorInX, VectorInY,      &
                                           VectorOutX, VectorOutY, WaterPoints2D,       &
                                           RotateX, RotateY, STAT)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:, :), pointer       :: VectorInX, VectorInY, VectorOutX, VectorOutY
        integer, dimension(:, :), pointer       :: WaterPoints2D
        integer                                 :: HorizontalGridID
        logical                                 :: RotateX, RotateY
        integer, optional                       :: STAT
        !Local-----------------------------------------------------------------
        real(8)                                 :: XGrid, YGrid, AngleX, AngleY
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j
        integer                                 :: STAT_, ready_
        real                                    :: GridRotationRadians

        !Begin--------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then



if1:        if  (Me%Distortion .or. Me%RegularRotation) then

                !Bounds
                ILB = Me%WorkSize%ILB
                IUB = Me%WorkSize%IUB

                JLB = Me%WorkSize%JLB
                JUB = Me%WorkSize%JUB

                GridRotationRadians = Me%Grid_Angle * Pi / 180.

                !Rotate
                do j = JLB, JUB
                do i = ILB, IUB

if2:                if (WaterPoints2D(i, j) == WaterPoint) then

                        !VectorX = X_Grid * cos(AngleX) + Y_Grid * cos(AngleY)
                        !VectorY = X_Grid * sin(AngleX) + Y_Grid * sin(AngleY)

                        AngleX = 0.
                        AngleY = Pi / 2.

                        if      (Me%Distortion) then

                            AngleX = Me%RotationX(i, j)
                            AngleY = Me%RotationY(i, j)

                        else if (Me%RegularRotation) then

                            AngleX = GridRotationRadians
                            AngleY = GridRotationRadians + Pi / 2.

                        endif

                        call FromGridToCartesian (VectorInX(i, j), VectorInY(i, j),           &
                                                  AngleX, AngleY, Xgrid, Ygrid)

                        if (RotateX) then

                            VectorOutX(i, j) = Xgrid

                        endif

                        if (RotateY) then

                            VectorOutY(i, j) = Ygrid

                        endif


                    endif if2

                enddo
                enddo

            else  if1


                if (RotateX) then

                    VectorOutX(:, :) = VectorInX(:, :)

                endif

                if (RotateY) then

                    VectorOutY(:, :) = VectorInY(:, :)

                endif


            endif if1

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT))  STAT = STAT_

    end subroutine RotateVectorGridToField2DR8
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------


    subroutine RotateVectorGridToField3D(HorizontalGridID, VectorInX, VectorInY,        &
                                         VectorOutX, VectorOutY, WaterPoints3D,         &
                                         RotateX, RotateY, KLB, KUB, STAT)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, :, :), pointer    :: VectorInX, VectorInY, VectorOutX, VectorOutY
        integer, dimension(:, :, :), pointer    :: WaterPoints3D
        logical                                 :: RotateX, RotateY
        integer                                 :: HorizontalGridID, KLB, KUB
        integer, optional                       :: STAT
        !Local-----------------------------------------------------------------
        real                                    :: XGrid, YGrid, AngleX, AngleY
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j, k
        integer                                 :: STAT_, ready_
        real                                    :: GridRotationRadians

        !Begin--------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then



if1:        if  (Me%Distortion .or. Me%RegularRotation) then

                !Bounds
                ILB = Me%WorkSize%ILB
                IUB = Me%WorkSize%IUB

                JLB = Me%WorkSize%JLB
                JUB = Me%WorkSize%JUB

                GridRotationRadians = Me%Grid_Angle * Pi / 180.

                !Rotate
                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB

if2:                if (WaterPoints3D(i, j, k) == WaterPoint) then

                        !VectorX = X_Grid * cos(AngleX) + Y_Grid * cos(AngleY)
                        !VectorY = X_Grid * sin(AngleX) + Y_Grid * sin(AngleY)

                        AngleX = 0.
                        AngleY = Pi / 2.

                        if      (Me%Distortion) then

                            AngleX = Me%RotationX(i, j)
                            AngleY = Me%RotationY(i, j)

                        else if (Me%RegularRotation) then

                            AngleX = GridRotationRadians
                            AngleY = GridRotationRadians + Pi / 2.

                        endif

                        call FromGridToCartesian(VectorInX(i, j, k), VectorInY(i, j, k),           &
                                                  AngleX, AngleY, Xgrid, Ygrid)


                        if (RotateX) then

                            VectorOutX(i, j, k) = Xgrid

                        endif

                        if (RotateY) then

                            VectorOutY(i, j, k) = Ygrid

                        endif


                    endif if2

                enddo
                enddo
                enddo

            else  if1

                if (RotateX) then

                    VectorOutX(:, :, :) = VectorInX(:, :, :)

                endif

                if (RotateY) then

                    VectorOutY(:, :, :) = VectorInY(:, :, :)

                endif

            endif if1

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT))  STAT = STAT_

    end subroutine RotateVectorGridToField3D
    !--------------------------------------------------------------------------

    subroutine RotateAngleGridToField2D(HorizontalGridID, AngleIn, OutReferential,       &
                                         AngleOut, WaterPoints2D,         &
                                         Rotate, STAT)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, :), pointer       :: AngleIn, AngleOut
        integer, dimension(:, :), pointer       :: WaterPoints2D
        integer                                 :: HorizontalGridID
        logical                                 :: Rotate
        integer, optional                       :: STAT
        integer                                 :: OutReferential
        !Local-----------------------------------------------------------------
        real                                    :: AngleOutCell, GridAngle
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j
        integer                                 :: STAT_, ready_

        !Begin--------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then



!if1:        !if  (Me%Distortion .or. Me%RegularRotation) then

                !Bounds
                ILB = Me%WorkSize%ILB
                IUB = Me%WorkSize%IUB

                JLB = Me%WorkSize%JLB
                JUB = Me%WorkSize%JUB


                !Rotate
                do j = JLB, JUB
                do i = ILB, IUB

if2:                if (WaterPoints2D(i, j) == WaterPoint) then

                        !VectorX = X_Grid * cos(AngleX) + Y_Grid * cos(AngleY)
                        !VectorY = X_Grid * sin(AngleX) + Y_Grid * sin(AngleY)

                        GridAngle = 0.

                        if      (Me%Distortion) then

                            !with distortion, the cell trigonometric circle origin is coincident with cell "x plane"
                            !and rotation Y is not accounted
                            GridAngle = Me%RotationX(i, j) * 180. / Pi

                        else if (Me%RegularRotation) then

                            GridAngle = Me%Grid_Angle

                        endif

                        call AngleFromGridToField (AngleIn(i, j), OutReferential,           &
                                                         GridAngle, AngleOutCell)

                        if (Rotate) then

                            AngleOut(i, j) = AngleOutCell

                        endif


                    endif if2

                enddo
                enddo

            !else  if1
            !
            !    if (Rotate) then
            !
            !        AngleOut(:, :) = AngleIn(:, :)
            !
            !    endif
            !
            !
            !endif if1

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT))  STAT = STAT_

    end subroutine RotateAngleGridToField2D

    !--------------------------------------------------------------------------

    subroutine RotateAngleGridToField3D(HorizontalGridID, AngleIn, OutReferential,       &
                                         AngleOut, WaterPoints3D,         &
                                         Rotate, KLB, KUB, STAT)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, : , :), pointer   :: AngleIn, AngleOut
        integer, dimension(:, : , :), pointer   :: WaterPoints3D
        integer                                 :: HorizontalGridID, KLB, KUB
        logical                                 :: Rotate
        integer                                 :: OutReferential
        integer, optional                       :: STAT
        !Local-----------------------------------------------------------------
        real                                    :: AngleOutCell, GridAngle
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j, k
        integer                                 :: STAT_, ready_

        !Begin--------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then



!if1:        !if  (Me%Distortion .or. Me%RegularRotation) then

                !Bounds
                ILB = Me%WorkSize%ILB
                IUB = Me%WorkSize%IUB

                JLB = Me%WorkSize%JLB
                JUB = Me%WorkSize%JUB

                !Rotate
                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB

if2:                if (WaterPoints3D(i, j, k) == WaterPoint) then

                        !VectorX = X_Grid * cos(AngleX) + Y_Grid * cos(AngleY)
                        !VectorY = X_Grid * sin(AngleX) + Y_Grid * sin(AngleY)

                        GridAngle = 0.

                        if      (Me%Distortion) then

                            !with distortion, the cell trigonometric circle origin is coincident with cell "x plane"
                            !and rotation Y is not accounted
                            GridAngle = Me%RotationX(i, j) * 180. / Pi

                        else if (Me%RegularRotation) then

                            GridAngle = Me%Grid_Angle

                        endif

                        call AngleFromGridToField (AngleIn(i, j, k), OutReferential,           &
                                                         GridAngle, AngleOutCell)

                        if (Rotate) then

                            AngleOut(i, j, k) = AngleOutCell

                        endif


                    endif if2

                enddo
                enddo
                enddo

            !else  if1
            !
            !    if (Rotate) then
            !
            !        AngleOut(:, :, :) = AngleIn(:, : , :)
            !
            !    endif
            !
            !
            !endif if1

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT))  STAT = STAT_

    end subroutine RotateAngleGridToField3D

    !--------------------------------------------------------------------------

    subroutine ComputeAngleFromGridComponents2D(HorizontalGridID, VectorU, VectorV,       &
                                         AngleOutField, AngleOutGrid, WaterPoints2D,  OutReferential,             &
                                         Rotate, STAT)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, :), pointer       :: VectorU, VectorV, AngleOutField, AngleOutGrid
        integer, dimension(:, :), pointer       :: WaterPoints2D
        integer                                 :: HorizontalGridID, OutReferential
        logical                                 :: Rotate
        integer, optional                       :: STAT
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, ready_, STAT_CALL
        integer                                 :: ILB, IUB, JLB, JUB, i, j

        !Begin--------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then


            !Bounds
            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            do j = JLB, JUB
            do i = ILB, IUB

                if (WaterPoints2D(i, j) == WaterPoint) then
                    AngleOutGrid(i,j) = atan2(VectorV(i,j), VectorU(i,j)) * 180. / Pi
                endif

            enddo
            enddo

            !rotate the angles to input referential
            call RotateAngleGridToField2D(HorizontalGridID,                         &
                                            AngleIn = AngleOutGrid,                 &
                                            OutReferential = OutReferential,        &
                                            AngleOut = AngleOutField,               &
                                            WaterPoints2D = WaterPoints2D,          &
                                            Rotate = Rotate,                        &
                                            STAT = STAT_CALL)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT))  STAT = STAT_

    end subroutine ComputeAngleFromGridComponents2D

    !--------------------------------------------------------------------------

    subroutine ComputeAngleFromGridComponents3D(HorizontalGridID, VectorU, VectorV,       &
                                         AngleOutField, AngleOutGrid, WaterPoints3D, OutReferential,         &
                                         Rotate, KLB, KUB, STAT)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, :, :), pointer    :: VectorU, VectorV, AngleOutField, AngleOutGrid
        integer, dimension(:, :, :), pointer    :: WaterPoints3D
        integer                                 :: HorizontalGridID, OutReferential
        integer                                 :: KLB, KUB
        logical                                 :: Rotate
        integer, optional                       :: STAT
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, ready_, STAT_CALL
        integer                                 :: ILB, IUB, JLB, JUB, i, j, k
        !Begin--------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ == IDLE_ERR_ .or. ready_ == READ_LOCK_ERR_) then


            !Bounds
            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB


            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB

                if (WaterPoints3D(i, j, k) == WaterPoint) then
                    AngleOutGrid(i,j, k) = atan2(VectorV(i,j,k), VectorU(i,j,k)) * 180. / Pi
                endif

            enddo
            enddo
            enddo

            !rotate the angles to input referential
            call RotateAngleGridToField3D(HorizontalGridID,                         &
                                            AngleIn = AngleOutGrid,                 &
                                            OutReferential = OutReferential,        &
                                            AngleOut = AngleOutField,               &
                                            WaterPoints3D = WaterPoints3D,          &
                                            Rotate = Rotate,                        &
                                            KLB    = KLB,                           &
                                            KUB    = KUB,                           &
                                            STAT = STAT_CALL)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT))  STAT = STAT_

    end subroutine ComputeAngleFromGridComponents3D

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillHorizontalGrid(HorizontalGridID, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: HorizontalGridID
        integer, optional                   :: STAT

        !Local-----------------------------------------------------------------
        integer                             :: STAT_, ready_, STATUS
        integer                             :: nUsers

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(HorizontalGridID, ready_)

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mHORIZONTALGRID_,  Me%InstanceID)

            if (nUsers == 0) then

                !XX
                deallocate (Me%XX, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR10'

                !YY
                deallocate (Me%YY, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR20'

                !XX_IE
                deallocate (Me%XX_IE, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR30'

                !YY_IE
                deallocate (Me%YY_IE, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR40'

                !LatitudeConn
                deallocate (Me%LatitudeConn, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR50'

                !LongitudeConn
                deallocate (Me%LongitudeConn, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR60'

                !LatitudeZ
                deallocate (Me%LatitudeZ, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR70'

                !LongitudeZ
                deallocate (Me%LongitudeZ, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR80'

                ! X,Y along grid
                deallocate (Me%XX_AlongGrid, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR320'

                deallocate (Me%YY_AlongGrid, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'AllocateVariables - HorizontalGrid - ERR330'


                !Distances
                deallocate (Me%DXX, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR90'

                deallocate (Me%DYY, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR100'

                deallocate (Me%DZX, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR110'

                deallocate (Me%DZY, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR120'

                deallocate (Me%DUX, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR130'

                deallocate (Me%DUY, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR140'

                deallocate (Me%DVX, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR150'

                deallocate (Me%DVY, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR160'

                !Coriolis
                deallocate (Me%F, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR170'

                !Compute points location

                if (.not. Me%CornersXYInput) then
                    !XX_Z
                    deallocate (Me%Compute%XX_Z, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR180'

                    !YY_Z
                    deallocate (Me%Compute%YY_Z, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR190'

                    !XX_U
                    deallocate (Me%Compute%XX_U, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR200'

                    !YY_V
                    deallocate (Me%Compute%YY_V, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR210'


                    !XX_V
                    deallocate (Me%Compute%XX_V, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR220'

                    !YY_U
                    deallocate (Me%Compute%YY_U, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR230'


                    !XX_Cross
                    deallocate (Me%Compute%XX_Cross, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR240'

                    !YY_Cross
                    deallocate (Me%Compute%YY_Cross, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR250'

                endif

                !XX2D_Z
                deallocate (Me%Compute%XX2D_Z, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR260'

                !YY2D_Z
                deallocate (Me%Compute%YY2D_Z, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR270'

                !XX2D_U
                deallocate (Me%Compute%XX2D_U, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR280'

                !YY2D_V
                deallocate (Me%Compute%YY2D_V, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR290'


                !XX2D_V
                deallocate (Me%Compute%XX2D_V, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR300'

                !YY2D_U
                deallocate (Me%Compute%YY2D_U, stat = STATUS)
                if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR310'


                !DefineCellsMap
                deallocate(Me%DefineCellsMap)

                if (Me%CoordType == CIRCULAR_ .or. Me%CornersXYInput) then

                    deallocate (Me%RotationX, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR340'

                    deallocate (Me%RotationY, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR350'

                endif


                if (Me%CornersXYInput) then

                    if (associated(Me%DefineFacesUMap)) deallocate(Me%DefineFacesUMap, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR360'

                    if (associated(Me%DefineFacesVMap)) deallocate(Me%DefineFacesVMap, stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR370'

                    if (associated(Me%DefineCrossMap )) deallocate(Me%DefineCrossMap,  stat = STATUS)
                    if (STATUS /= SUCCESS_) stop 'KillHorizontalGrid - HorizontalGrid - ERR380'

                endif

                !Border polygons

!                if (Me%GridBorderCoord%Type_ /= Rectang_) then
                    deallocate(Me%GridBorderCoord%Polygon_%VerticesF)
                    nullify   (Me%GridBorderCoord%Polygon_%VerticesF)

                    deallocate(Me%GridBorderCoord%Polygon_)
                    nullify   (Me%GridBorderCoord%Polygon_)

!                endif

                deallocate(Me%GridOutBorderCoord%Polygon_%VerticesF)
                nullify   (Me%GridOutBorderCoord%Polygon_%VerticesF)

                deallocate(Me%GridOutBorderCoord%Polygon_)
                nullify   (Me%GridOutBorderCoord%Polygon_)

                deallocate(Me%GridBorderCart%Polygon_%VerticesF)
                nullify   (Me%GridBorderCart%Polygon_%VerticesF)

                deallocate(Me%GridBorderCart%Polygon_)
                nullify   (Me%GridBorderCart%Polygon_)

                deallocate(Me%GridOutBorderCart%Polygon_%VerticesF)
                nullify   (Me%GridOutBorderCart%Polygon_%VerticesF)

                deallocate(Me%GridOutBorderCart%Polygon_)
                nullify   (Me%GridOutBorderCart%Polygon_)

                deallocate(Me%GridBorderAlongGrid%Polygon_%VerticesF)
                nullify   (Me%GridBorderAlongGrid%Polygon_%VerticesF)

                deallocate(Me%GridBorderAlongGrid%Polygon_)
                nullify   (Me%GridBorderAlongGrid%Polygon_)


                deallocate(Me%GridBorderCart     )
                deallocate(Me%GridOutBorderCart  )
                deallocate(Me%GridBorderCoord    )
                deallocate(Me%GridOutBorderCoord )
                deallocate(Me%GridBorderAlongGrid)


                nullify   (Me%GridBorderCart     )
                nullify   (Me%GridOutBorderCart  )
                nullify   (Me%GridBorderCoord    )
                nullify   (Me%GridOutBorderCoord )
                nullify   (Me%GridBorderAlongGrid)


                deallocate (Me%XX1D_Aux)
                deallocate (Me%YY1D_Aux)

                nullify    (Me%XX1D_Aux)
                nullify    (Me%YY1D_Aux)

                deallocate(Me%AuxPolygon%VerticesF)
                deallocate(Me%AuxPolygon)

                if (allocated(Me%Connections_U)) then
                    deallocate(Me%Connections_U)
                    deallocate(Me%IWD_Distances_U)
                    nullify   (Me%IWD_Distances_U)
                endif
                if (allocated(Me%Connections_V)) then
                    deallocate(Me%Connections_V)
                    deallocate(Me%IWD_Distances_V)
                    nullify   (Me%IWD_Distances_V)
                endif
                if (allocated(Me%Connections_Z)) then
                    deallocate(Me%Connections_Z)
                endif
                
                if (allocated(Me%IWD_Distances_Z)) then
                    deallocate(Me%IWD_Distances_Z)
                endif

                call KillFatherGridList

                call KillDDecomp


                !ObjHorizontalGrid
                call DeallocateInstance

                HorizontalGridID = 0
                STAT_            = SUCCESS_

            end if

        else

            STAT_ = ready_

        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine KillHorizontalGrid

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_HorizontalGrid), pointer          :: AuxObjHorizontalGrid
        type (T_HorizontalGrid), pointer          :: PreviousObjHorizontalGrid

        !Updates pointers
        if (Me%InstanceID == FirstHorizontalGrid%InstanceID) then
            FirstHorizontalGrid => FirstHorizontalGrid%Next
        else
            PreviousObjHorizontalGrid => FirstHorizontalGrid
            AuxObjHorizontalGrid      => FirstHorizontalGrid%Next
            do while (AuxObjHorizontalGrid%InstanceID /= Me%InstanceID)
                PreviousObjHorizontalGrid => AuxObjHorizontalGrid
                AuxObjHorizontalGrid      => AuxObjHorizontalGrid%Next
            enddo

            !Now update linked list
            PreviousObjHorizontalGrid%Next => AuxObjHorizontalGrid%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me)

    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------

    subroutine KillFatherGridList

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_FatherGrid)    , pointer    :: FatherGrid
        integer                             :: status

        ! Deallocates all the water properties

        FatherGrid => Me%FirstFatherGrid

do1 :   do while(associated(FatherGrid))

            if (FatherGrid%OkZ) then


                if (.not. FatherGrid%CornersXYInput) then

                    !XX_Z
                    deallocate(FatherGrid%XX_Z, STAT = status)

                    if (status /= SUCCESS_)                                             &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR10")

                    !YY_Z
                    deallocate(FatherGrid%YY_Z, STAT = status)

                    if (status /= SUCCESS_)                                             &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR20")

                endif

                nullify   (FatherGrid%XX_Z)
                nullify   (FatherGrid%YY_Z)

                !IZ
                deallocate(FatherGrid%IZ, STAT = status)
                if (status /= SUCCESS_)                                                 &
                    call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR30")

                nullify   (FatherGrid%IZ)

                !JZ
                deallocate(FatherGrid%JZ, STAT = status)
                if (status /= SUCCESS_)                                                 &
                    call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR40")

                nullify   (FatherGrid%JZ)
                
                if (associated(FatherGrid%ILinkZ)) then
                    !ILinkZ
                    deallocate(FatherGrid%ILinkZ, STAT = status)
                    if (status /= SUCCESS_)                                                 &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR41") 

                    nullify   (FatherGrid%ILinkZ)
                endif
                
                if (associated(FatherGrid%JLinkZ)) then
                    !JLinkZ
                    deallocate(FatherGrid%JLinkZ, STAT = status)
                    if (status /= SUCCESS_)                                                 &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR42") 

                    nullify   (FatherGrid%JLinkZ) 
                    
                endif                    

            endif

            if (FatherGrid%OkU) then

                if (.not. FatherGrid%CornersXYInput) then

                    !XX_U
                    deallocate(FatherGrid%XX_U, STAT = status)

                    if (status /= SUCCESS_)                                             &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR50")

                    !YY_U
                    deallocate(FatherGrid%YY_U, STAT = status)

                    if (status /= SUCCESS_)                                             &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR60")

                endif

                nullify   (FatherGrid%XX_U)
                nullify   (FatherGrid%YY_U)


                !IU
                deallocate(FatherGrid%IU, STAT = status)

                if (status /= SUCCESS_)                                                 &
                    call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR70")

                nullify   (FatherGrid%IU)

                !JU
                deallocate(FatherGrid%JU, STAT = status)

                if (status /= SUCCESS_)                                                 &
                    call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR80")

                nullify   (FatherGrid%JU)
                
                if (associated(FatherGrid%ILinkU)) then
                    !ILinkU
                    deallocate(FatherGrid%ILinkU, STAT = status)

                    if (status /= SUCCESS_)                                                 &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR81") 

                    nullify   (FatherGrid%ILinkU)
                    
                endif                                        
                
                if (associated(FatherGrid%JLinkU)) then                

                    deallocate(FatherGrid%JLinkU, STAT = status)

                    if (status /= SUCCESS_)                                                 &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR82") 

                    nullify   (FatherGrid%JLinkU)
                endif
            endif

            if (FatherGrid%OkV) then

                if (.not. FatherGrid%CornersXYInput) then

                    !YY_V
                    deallocate(FatherGrid%YY_V, STAT = status)

                    if (status /= SUCCESS_)                                             &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR90")

                    !XX_V
                    deallocate(FatherGrid%XX_V, STAT = status)

                    if (status /= SUCCESS_)                                             &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR100")
                endif

                nullify   (FatherGrid%XX_V)
                nullify   (FatherGrid%YY_V)

                !IV
                deallocate(FatherGrid%IV, STAT = status)

                if (status /= SUCCESS_)                                                 &
                    call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR110")

                nullify   (FatherGrid%IV)

                !JV
                deallocate(FatherGrid%JV, STAT = status)

                if (status /= SUCCESS_)                                                 &
                    call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR120")

                nullify   (FatherGrid%JV)
                
                if (associated(FatherGrid%ILinkV)) then             
                    !ILinkV
                    deallocate(FatherGrid%ILinkV, STAT = status)

                    if (status /= SUCCESS_)                                                 &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR111") 

                    nullify   (FatherGrid%ILinkV)
                endif
                
                if (associated(FatherGrid%JLinkV)) then
                    !JLinkV
                    deallocate(FatherGrid%JLinkV, STAT = status)

                    if (status /= SUCCESS_)                                                 &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR112") 

                    nullify   (FatherGrid%JLinkV)
                endif
                
            endif

            if (FatherGrid%OkCross) then

                if (.not. FatherGrid%CornersXYInput) then

                    !XX_Cross
                    deallocate(FatherGrid%XX_Cross, STAT = status)

                    if (status /= SUCCESS_)                                             &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR130")

                    !YY_Cross
                    deallocate(FatherGrid%YY_Cross, STAT = status)

                    if (status /= SUCCESS_)                                             &
                        call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR140")

                endif

                nullify   (FatherGrid%XX_Cross)
                nullify   (FatherGrid%YY_Cross)



                !ICross
                deallocate(FatherGrid%ICross, STAT = status)

                if (status /= SUCCESS_)                                                 &
                    call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR150")

                nullify   (FatherGrid%ICross)

                !JCross
                deallocate(FatherGrid%JCross, STAT = status)

                if (status /= SUCCESS_)                                                 &
                    call SetError(FATAL_, INTERNAL_, "KillFatherGridList; HorizontalGrid. ERR160")

                nullify   (FatherGrid%JCross)

            endif

            FatherGrid => FatherGrid%Next

        end do do1

        nullify   (Me%FirstFatherGrid,Me%LastFatherGrid)

    end subroutine KillFatherGridList

    subroutine KillDDecomp

#ifdef _USE_MPI

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL


        !Begin-----------------------------------------------------------------

        if (associated(Me%DDecomp%Slaves_MPI_ID))    then
            deallocate(Me%DDecomp%Slaves_MPI_ID)
            nullify   (Me%DDecomp%Slaves_MPI_ID)
        endif

        if (associated(Me%DDecomp%Slaves_Size))      then
            deallocate(Me%DDecomp%Slaves_Size)
            nullify   (Me%DDecomp%Slaves_Size)
        endif

        if (associated(Me%DDecomp%Slaves_Inner))     then
            deallocate(Me%DDecomp%Slaves_Inner)
            nullify   (Me%DDecomp%Slaves_Inner)
        endif

        if (associated(Me%DDecomp%Slaves_Mapping))   then
            deallocate(Me%DDecomp%Slaves_Mapping)
            nullify   (Me%DDecomp%Slaves_Mapping)
        endif

        if (associated(Me%DDecomp%Slaves_HaloMap))   then
            deallocate(Me%DDecomp%Slaves_HaloMap)
            nullify   (Me%DDecomp%Slaves_HaloMap)
        endif

        if (associated(Me%DDecomp%Interfaces))       then
            deallocate(Me%DDecomp%Interfaces)
            nullify   (Me%DDecomp%Interfaces)
        endif

        if (Me%DDecomp%FilesListOpen) then

            call UnitsManager(Me%DDecomp%FilesListID, FileClose, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'KillDDecomp - HorizontalGrid - ERR10'

        endif

        call DeAllocateMasterCoef2D

        call DeAllocateMasterCoef3D

        !----------------------------------------------------------------------

#endif _USE_MPI

    end subroutine KillDDecomp


    subroutine DeAllocateMasterCoef2D

        !Arguments------------------------------------------------------------

        !Local---------------------------------------------------------------

        !Begin---------------------------------------------------------------

        if  (Me%DDecomp%Coef2D%AllocateON .and. Me%DDecomp%Master) then

            deallocate(Me%DDecomp%Coef2D%D        )
            deallocate(Me%DDecomp%Coef2D%E        )
            deallocate(Me%DDecomp%Coef2D%F        )
            deallocate(Me%DDecomp%Coef2D%Ti       )
            deallocate(Me%DDecomp%Coef2D%Results2D)
            deallocate(Me%DDecomp%Coef2D%VECG     )
            deallocate(Me%DDecomp%Coef2D%VECW     )

        endif

    end subroutine DeAllocateMasterCoef2D

    !End----------------------------------------------------------------


    subroutine DeAllocateMasterCoef3D

        !Arguments------------------------------------------------------------

        !Local---------------------------------------------------------------

        !Begin---------------------------------------------------------------

        if  (Me%DDecomp%Coef3D%AllocateON .and. Me%DDecomp%Master) then

            deallocate(Me%DDecomp%Coef3D%D        )
            deallocate(Me%DDecomp%Coef3D%E        )
            deallocate(Me%DDecomp%Coef3D%F        )
            deallocate(Me%DDecomp%Coef3D%Ti       )
            deallocate(Me%DDecomp%Coef3D%Results3D)
            deallocate(Me%DDecomp%Coef3D%VECG     )
            deallocate(Me%DDecomp%Coef3D%VECW     )

        endif

    end subroutine DeAllocateMasterCoef3D

    !End----------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Ready (ObjHorizontalGrid_ID, ready_)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjHorizontalGrid_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjHorizontalGrid_ID > 0) then
            call LocateObjHorizontalGrid (ObjHorizontalGrid_ID)
            ready_ = VerifyReadLock (mHORIZONTALGRID_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    function Ready_ThreadSafe (ObjHorizontalGrid_ID, ready_) result(LocalMe)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )              :: ObjHorizontalGrid_ID
        integer,           intent(OUT)              :: ready_
        type (T_HorizontalGrid), pointer            :: LocalMe

        !----------------------------------------------------------------------

        nullify (LocalMe)

cd1:    if (ObjHorizontalGrid_ID > 0) then
            LocalMe => LocateObjHorizontalGrid_ThreadSafe(ObjHorizontalGrid_ID)
            ready_ = VerifyReadLock (mHORIZONTALGRID_, LocalMe%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end function Ready_ThreadSafe

    !--------------------------------------------------------------------------

    subroutine LocateObjHorizontalGrid (ObjHorizontalGridID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjHorizontalGridID

        !Local-----------------------------------------------------------------

        Me => FirstHorizontalGrid
        do while (associated (Me))
            if (Me%InstanceID == ObjHorizontalGridID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'HorizontalGrid - LocateObjHorizontalGrid - ERR01'

    end subroutine LocateObjHorizontalGrid

    !--------------------------------------------------------------------------

    function LocateObjHorizontalGrid_ThreadSafe (ObjHorizontalGridID) result(LocatedMe)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )              :: ObjHorizontalGridID
        type (T_HorizontalGrid), pointer            :: LocatedMe

        !Local-----------------------------------------------------------------

        LocatedMe => FirstHorizontalGrid
        do while (associated (LocatedMe))
            if (LocatedMe%InstanceID == ObjHorizontalGridID) exit
            LocatedMe => LocatedMe%Next
        enddo

        if (.not. associated(LocatedMe)) stop 'HorizontalGrid - LocateObjHorizontalGrid_ThreadSafe - ERR01'

    end function LocateObjHorizontalGrid_ThreadSafe

    !--------------------------------------------------------------------------
    subroutine LocateObjFather (ObjHorizontalGrid, ObjHorizontalGridID)

        !Arguments-------------------------------------------------------------
        type (T_HorizontalGrid), pointer            :: ObjHorizontalGrid
        integer                                     :: ObjHorizontalGridID

        !Local-----------------------------------------------------------------

        ObjHorizontalGrid => FirstHorizontalGrid
        do while (associated (ObjHorizontalGrid))
            if (ObjHorizontalGrid%InstanceID == ObjHorizontalGridID) exit
            ObjHorizontalGrid => ObjHorizontalGrid%Next
        enddo

        if (.not. associated(ObjHorizontalGrid)) then
            stop 'HorizontalGrid - LocateObjFather - ERR01'
        endif

    end subroutine LocateObjFather


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !AUXILARY AUXILARY AUXILARY AUXILARY AUXILARY AUXILARY AUXILARY AUXILARY AU

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine wgs84_to_rd(lat_wgs84, lon_wgs84, x_rd, y_rd)

        !Source code converted from C code by Ejo Schrama <e.j.o.schrama@lr.tudelft.nl>
        !The source is provided as is and has only been tested within the range of
        !the RD coordinate system and the accuracy of the conversion is no better than 50 cm


        !Arguments-------------------------------------------------------------
        real(8), intent(in)                :: lat_wgs84, lon_wgs84
        real(8), intent(out)               :: x_rd, y_rd

        !Local-----------------------------------------------------------------
        real(8)                             :: lambes, phibes

        !Begin-----------------------------------------------------------------


        call wgs842bessel_(lat_wgs84, lon_wgs84, phibes, lambes)

        call bessel2rd_(phibes, lambes, x_rd, y_rd)

    end subroutine


    subroutine wgs842bessel_(phiwgs, lamwgs, phibes, lambes)

        !Source code converted from C code by Ejo Schrama <e.j.o.schrama@lr.tudelft.nl>
        !The source is provided as is and has only been tested within the range of
        !the RD coordinate system and the accuracy of the conversion is no better than 50 cm

        !Arguments-------------------------------------------------------------
        real(8), intent(in)                 :: phiwgs, lamwgs
        real(8), intent(out)                :: phibes, lambes

        !Local-----------------------------------------------------------------
        real(8)                             :: dlam, dphi, lamcor, phicor

        !Begin-----------------------------------------------------------------

        dphi = phiwgs - 52.
        dlam = lamwgs - 5.
        phicor = (-96.862 - dphi * 11.714 - dlam * .125) * 1e-5
        lamcor = (dphi * 0.329 - 37.902 - dlam * 14.667) * 1e-5
        phibes = phiwgs - phicor
        lambes = lamwgs - lamcor

    end subroutine wgs842bessel_


    subroutine rd_to_wgs84(x_rd, y_rd, lat_wgs84, lon_wgs84)

        !Source code converted from C code by Ejo Schrama <e.j.o.schrama@lr.tudelft.nl>
        !The source is provided as is and has only been tested within the range of
        !the RD coordinate system and the accuracy of the conversion is no better than 50 cm

        !Arguments-------------------------------------------------------------
        real(8), intent(in)                :: x_rd, y_rd
        real(8), intent(out)               :: lat_wgs84, lon_wgs84

        !Local-----------------------------------------------------------------
        real(8)                             :: phibes, lambes

        !Begin-----------------------------------------------------------------

        call rd2bessel_(x_rd, y_rd, phibes, lambes)

        call bessel2wgs84_(phibes, lambes, lat_wgs84, lon_wgs84)

    end subroutine



    subroutine rd2bessel_(x, y, phi, lambda)

        !Source code converted from C code by Ejo Schrama <e.j.o.schrama@lr.tudelft.nl>
        !The source is provided as is and has only been tested within the range of
        !the RD coordinate system and the accuracy of the conversion is no better than 50 cm


        !Arguments-------------------------------------------------------------
        real(8)                             :: x, y, phi, lambda

        !Local-----------------------------------------------------------------
        real(8)                             :: d_1, d_2

        real(8)                             :: bigr, cpsi, spsi
        real(8)                             :: phiprime, a, b, e
        real(8)                             :: k, m, n, q, r, w
        real(8)                             :: b0, l0, x0, y0, ca, cb, dl, sa, sb
        real(8)                             :: dq, pi, lambda0, sdl, psi, phi0
        integer                             :: i

        !Begin-----------------------------------------------------------------

        x0 = 1.55e5
        y0 = 4.63e5
        k = .9999079
        bigr = 6382644.571
        m = .003773953832
        n = 1.00047585668
        pi = atan(1.) * 4.
        lambda0 = pi * .029931327161111111
        phi0 = pi * .28975644753333335
        l0 = pi * .029931327161111111
        b0 = pi * .28956165138333334
        e = .08169683122
        a = 6377397.155




        d_1 = x - x0
        d_2 = y - y0
        r = sqrt(d_1 * d_1 + d_2 * d_2)


        if(r .ne. 0.)then
            sa = (x - x0) / r
            ca = (y - y0) / r
        else
            sa = 0.
            ca = 0.
        endif


        psi = atan2(r, k * 2. * bigr) * 2.
        cpsi = cos(psi)
        spsi = sin(psi)
        sb = ca * cos(b0) * spsi + sin(b0) * cpsi
        d_1 = sb
        cb = sqrt(1. - d_1 * d_1)
        b = acos(cb)
        sdl = sa * spsi / cb
        dl = asin(sdl)
        lambda = dl / n + lambda0
        w = log(tan(b / 2. + pi / 4.))
        q = (w - m) / n
        phiprime = atan(exp(q)) * 2. - pi / 2.


        do i = 1, 4

            dq = e / 2. * log((e * sin(phiprime) + 1.) / (1. - e * sin(phiprime)))
            phi = atan(exp(q + dq)) * 2. - pi / 2.
            phiprime = phi

        enddo

        lambda = lambda / pi * 180.
        phi    = phi / pi * 180.


    end subroutine



    subroutine bessel2wgs84_(phibes, lambes, phiwgs, lamwgs)

        !Source code converted from C code by Ejo Schrama <e.j.o.schrama@lr.tudelft.nl>
        !The source is provided as is and has only been tested within the range of
        !the RD coordinate system and the accuracy of the conversion is no better than 50 cm


        !Arguments-------------------------------------------------------------
        real(8)                             :: phibes, lambes, phiwgs, lamwgs

        !Local-----------------------------------------------------------------
        real(8)                             :: dlam, dphi, lamcor, phicor


        !Begin-----------------------------------------------------------------

        dphi = phibes - 52.
        dlam = lambes - 5.
        phicor = (-96.862 - dphi * 11.714 - dlam * .125) * 1e-5
        lamcor = (dphi * 0.329 - 37.902 - dlam * 14.667) * 1e-5
        phiwgs = phibes + phicor
        lamwgs = lambes + lamcor

    end subroutine bessel2wgs84_


    subroutine bessel2rd_(argphi, arglam, x, y)

        !Source code converted from C code by Ejo Schrama <e.j.o.schrama@lr.tudelft.nl>
        !The source is provided as is and has only been tested within the range of
        !the RD coordinate system and the accuracy of the conversion is no better than 50 cm


        !Arguments-------------------------------------------------------------
        real(8)                             :: argphi, arglam, x, y

        !Local-----------------------------------------------------------------

        real(8)                             :: bigr, cpsi, cpsihalf, spsi
        real(8)                             :: spsihalf, tpsihalf, a, b, e
        real(8)                             :: s2psihalf, k, m, n, q, r, w
        real(8)                             :: b0, l0, x0, y0, ca, lambda, dl, sa
        real(8)                             :: dq, pi, qprime, lambda0, phi, phi0
        real(8)                             :: d_1, d_2

        !Begin-----------------------------------------------------------------

        x0 = 1.55e5
        y0 = 4.63e5
        k = .9999079
        bigr = 6382644.571
        m = .003773953832
        n = 1.00047585668
        pi = atan(1.) * 4.
        lambda0 = pi * .029931327161111111
        phi0 = pi * .28975644753333335
        l0 = pi * .029931327161111111
        b0 = pi * .28956165138333334
        e = .08169683122
        a = 6377397.155





        phi = argphi / 180. * pi
        lambda = arglam / 180. * pi
        qprime = log(tan(phi / 2. + pi / 4.))
        dq = e / 2. * log((e * sin(phi) + 1.) / (1. - e * sin(phi)))
        q = qprime - dq
        w = n * q + m
        b = atan(exp(w)) * 2. - pi / 2.
        dl = n * (lambda - lambda0)

        !/* Computing 2nd power */
        d_1 = sin((b - b0) / 2.)
        !/* Computing 2nd power */
        d_2 = sin(dl / 2.)
        s2psihalf = d_1 * d_1 + d_2 * d_2 * cos(b) * cos(b0)
        cpsihalf = sqrt(1. - s2psihalf)
        spsihalf = sqrt(s2psihalf)
        tpsihalf = spsihalf / cpsihalf
        spsi = spsihalf * 2. * cpsihalf
        cpsi = 1. - s2psihalf * 2.
        sa = sin(dl) * cos(b) / spsi
        ca = (sin(b) - sin(b0) * cpsi) / (cos(b0) * spsi)
        r = k * 2. * bigr * tpsihalf
        x = r * sa + x0
        y = r * ca + y0



    end subroutine bessel2rd_



!**********************************************************************
!                                                                     *
!     THIS PROGRAM IS SHAREWARE AND HAS BEEN DOWNLOADED FROM INTERNET *
!                                                                     *
!     THIS PROGRAM CONVERTS GPS TO UNIVERSIAL TRANSVERSE MERACTOR     *
!     COORDINATES AND VICE VERSA FOR THE NAD27 AND NAD83 DATUM.       *
!     THIS PROGRAM WAS WRITTEN BY E. CARLSON                          *
!     SUBROUTINES TMGRID, TCONST, TMGEOD, TCONPC,                     *
!     WERE WRITTEN BY T. VINCENTY, NGS, IN JULY 1984 .                *
!     THE ORGINAL PROGRAM WAS WRITTEN IN SEPTEMBER OF 1988.           *
!                                                                     *
!     THIS PROGRAM WAS UPDATED ON FEBUARY 16, 1989.  THE UPDATE WAS   *
!     HAVING THE OPTION OF SAVING AND *81* RECORD FILE.               *
!                                                                     *
!                                                                     *
!     THIS PROGRAM WAS UPDATED ON APRIL 3, 1990.  THE FOLLOWING UPDATE*
!     WERE MADE:                                                      *
!      1. CHANGE FROM JUST A CHOICE OF NAD27 OF NAD83 REFERENCE       *
!         ELLIPSOIDS TO; CLARKE 1866, GRS80/WGS84, INTERNATIONAL, AND *
!         ALLOW FOR USER DEFINED OTHER.                               *
!      2. ALLOW USE OF LATITUDES IN SOUTHERN HEMISPHERE AND LONGITUDES*
!         UP TO 360 DEGREES WEST.                                     *
!                                                                     *
!**********************************************************************
!                  DISCLAIMER                                         *
!                                                                     *
!   THIS PROGRAM AND SUPPORTING INFORMATION IS FURNISHED BY THE       *
! GOVERNMENT OF THE UNITED STATES OF AMERICA, AND IS ACCEPTED AND     *
! USED BY THE RECIPIENT WITH THE UNDERSTANDING THAT THE UNITED STATES *
! GOVERNMENT MAKES NO WARRANTIES, EXPRESS OR IMPLIED, CONCERNING THE  *
! ACCURACY, COMPLETENESS, RELIABILITY, OR SUITABILITY OF THIS         *
! PROGRAM, OF ITS CONSTITUENT PARTS, OR OF ANY SUPPORTING DATA.       *
!                                                                     *
!   THE GOVERNMENT OF THE UNITED STATES OF AMERICA SHALL BE UNDER NO  *
! LIABILITY WHATSOEVER RESULTING FROM ANY USE OF THIS PROGRAM.  THIS  *
! PROGRAM SHOULD NOT BE RELIED UPON AS THE SOLE BASIS FOR SOLVING A   *
! PROBLEM WHOSE INCORRECT SOLUTION COULD RESULT IN INJURY TO PERSON   *
! OR PROPERTY.                                                        *
!                                                                     *
!   THIS PROGRAM IS PROPERTY OF THE GOVERNMENT OF THE UNITED STATES   *
! OF AMERICA.  THEREFORE, THE RECIPIENT FURTHER AGREES NOT TO ASSERT  *
! PROPRIETARY RIGHTS THEREIN AND NOT TO REPRESENT THIS PROGRAM TO     *
! ANYONE AS BEING OTHER THAN A GOVERNMENT PROGRAM.                    *
!                                                                     *
!**********************************************************************
! TRANSFORMADO EM SUBROTINA DE CHAMADA (FLAVIO)
!
!     SIGNIFICADO DAS VARIAVEIS:
!
!     (in) icode: zona UTM para inum=2
!     (in)valores de inum:
!       ' 1  GEODETIC POSITIONS TO UTM COORDINATES '/,
!       ' 2  UTM COORDINATES TO GEODETIC POSITIONS'/,
!       ' 4  GEODETIC POSITIONS TO PORTUGUESE COORDINATES (HIDROMOD)'/,
!       ' 5  PORTUGUESE COORDINATES TO GEODETIC POSITIONS (HIDROMOD)'/,
!
!     (in) or (out)  REAIS DE DUPLA PRECISAO
!            ext_lat LATITUDE EM GRAUS (POSITIVO PARA NORTE)
!            ext_long LONGITUDE EM GRAUS (POSITIVO PARA ESTE)
!            ext_x COORDENADA X EM METROS (POSITIVO PARA ESTE)
!            ext_y COORDENADA Y EM METROS (POSITIVO PARA NORTE)
!***********************************************************************
      SUBROUTINE USCONVCO (inum, icode, ext_lat, ext_long, ext_x, ext_y)
!
!      IMPLICIT DOUBLEPRECISION (A - H, O - Z)
        !Arguments-------------------------------------------------------------
        real(8)                                     :: ext_lat, ext_long
        real(8)                                     :: ext_x, ext_y
        integer                                     :: inum, icode, LD, LM, LOD, LOM
        character(1)                                :: DATNUM, NORS, EORW
        logical                                     :: MILP
        !Local-----------------------------------------------------------------
        real(8)                                     :: PI, RAD, SLAT, SLON
        real(8)                                     :: north, east, lat, lon
        real(8)                                     :: ER, RF, F, ESQ

        COMMON / MIL / MILP
        COMMON / CONST / RAD, ER, RF, ESQ, PI
        COMMON / DATUM / DATNUM
        COMMON / XY / NORTH, EAST



      PI = 4.D0 * DATAN (1.D0)
      RAD = 180.D0 / PI
      MILP = .FALSE.

! passagem para as variaveis internas
      IF (INUM.EQ.1.OR.INUM.EQ.4) THEN
         IF (ext_lat.gt.0.) then
            NORS = 'N'
         ELSE
            NORS = 'S'
         ENDIF
         IF (dabs (ext_long) .lt.180) then
            IF (ext_long.lt.0.) then
               EORW = 'W'
            ELSE
               EORW = 'E'
            ENDIF
         ELSEIF (dabs (ext_long) .lt.360) then
            IF (ext_long.lt.0.) then
               ext_long = 360. - dabs (ext_long)
               EORW = 'E'
            ELSE
               ext_long = 360. - dabs (ext_long)
               EORW = 'W'
            ENDIF
         ELSE
            WRITE ( * , * ) 'erro na longitude em usconvco'
            STOP
         ENDIF
         ext_lat = DABS (ext_lat) / RAD
         ext_long = DABS (ext_long) / RAD
         CALL TODMS (ext_lat, LD, LM, SLAT)
         CALL TODMS (ext_long, LOD, LOM, SLON)
      ELSEIF (INUM.EQ.2.OR.INUM.EQ.5) then
         north = ext_y
         east = ext_x
      ENDIF
!***
                                                     !escolha do elipsoi
      IF (INUM.NE.6) CALL DATUMM (ER, RF, F, ESQ, DATNUM)
!
!**   USE THE NUM TO DO THE CORRECT FUNCTION

      IF (INUM.EQ.1.OR.INUM.EQ.4) THEN
                                   !Geodetic p. -> Portuguese Mil. Coord
         IF (INUM.EQ.4) MILP = .TRUE.
         CALL GPUT83 (LD, LM, SLAT, NORS, LOD, LOM, SLON, EORW)
      ELSEIF (INUM.EQ.2.OR.INUM.EQ.5) THEN
                                   !Portuguese Mil. Coordinates -> Geode
         IF (INUM.EQ.5) MILP = .TRUE.
         CALL UTGP83 (ICODE, lat, lon)
      ENDIF
      IF (INUM.EQ.1.OR.INUM.EQ.4) THEN
         ext_x = east
         ext_y = north
      ELSE
         ext_lat  = lat * rad
         ext_long = lon * rad
      ENDIF

      END SUBROUTINE USCONVCO

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    SUBROUTINE UTGP83 (ICODE, lat, lon)

        !IMPLICIT DOUBLEPRECISION (A - H, O - Z)
        !Arguments-------------------------------------------------------------
        integer                                     :: ICODE
        real(8)                                     :: lat, lon

        !Local-----------------------------------------------------------------
        real(8)                                     :: RAD, ER, RF, ESQ, PI
        common / CONST / RAD, ER, RF, ESQ, PI


        CALL IUTPC (ICODE, lat, lon)

    END SUBROUTINE UTGP83

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    SUBROUTINE GPUT83 (LD, LM, SLAT, NORS, LOD, LOM, SLON, EORW)

        !IMPLICIT DOUBLEPRECISION (A - H, O - Z)
        !Arguments-------------------------------------------------------------
        real(8)                                     :: SLAT, SLON
        integer                                     :: LD, LM, LOD, LOM
        character(1)                                :: NORS, EORW

        !Local-----------------------------------------------------------------
        real(8)                                     :: RAD, ER, RF, ESQ, PI
        COMMON / CONST / RAD, ER, RF, ESQ, PI


        CALL IUTGP (LD, LM, SLAT, NORS, LOD, LOM, SLON, EORW)

    END SUBROUTINE GPUT83

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    SUBROUTINE IUTPC (ICODE, lat, lon)

        !IMPLICIT DOUBLEPRECISION (A - H, O - Z)
        !Arguments-------------------------------------------------------------
        integer                                     :: ICODE
        real(8)                                     :: lat, lon

        !Local-----------------------------------------------------------------
        logical                                     :: FILFLAG, FILPRT
        character(1)                                :: DATNUM
        character(30)                               :: NAME
        character(80)                               :: CARDR
        real(8)                                     :: north, east
        logical                                     :: MILP

        COMMON / MIL / MILP
        COMMON / XY / NORTH, EAST
        COMMON / DATUM / DATNUM


        FILFLAG = .FALSE.
        FILPRT = .FALSE.
        NAME = '                              '
        WRITE (CARDR, 41) NAME
        41 FORMAT(T15,A30)

        ! Only to give it a number!
        IF (MILP)                                                             &
            ICODE = 29

        CALL DRUTGP (ICODE, lat, lon)

    END SUBROUTINE IUTPC

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

      SUBROUTINE IUTGP (LD, LM, SLAT, NORS, LOD, LOM, SLON, EORW)

      !IMPLICIT DOUBLEPRECISION (A - H, O - Z)
        !Arguments-------------------------------------------------------------
        real(8)                                     :: SLAT, SLON
        integer                                     :: LD, LM, LOD, LOM, ISEC, JSEC
        character(1)                                :: EORW
        character(1)                                :: NORS
        character(11)                               :: CLAT
        character(12)                               :: CLON
        character(30)                               :: NAME
        character(80)                               :: CARDR
        logical FILFLAG
        logical FIL81FL

      FILFLAG = .FALSE.
      FIL81FL = .FALSE.
      NAME = '                              '

      ISEC = SLAT * 1.0D5 + 0.5D0
      JSEC = SLON * 1.0D5 + 0.5D0

      WRITE (CLAT, 67) LD, LM, ISEC
   67 FORMAT(I6.2,I6.2,I14.7)

      WRITE (CLON, 68) LOD, LOM, JSEC
   68 FORMAT(I6.3,I6.2,I14.7)

      WRITE (CARDR, 70) NAME, CLAT, NORS, CLON, EORW
   70 FORMAT(T7,'*80*',T15,A30,T45,A11,A1,T57,A12,A1)

      CALL DRGPUT (CARDR)

      RETURN
      END SUBROUTINE IUTGP
!***********************************************************************
      SUBROUTINE DRGPUT (CARDR)
!*********************************************************************
!
!
!      THIS IS THE DRIVER TO COMPUTE UTM NORTHINGS AND EASTINGS
!      FOR EACH PRIMARY ZONE AND THE AJACENT ZONE IF THE LONGITUDE
!      IS WITH 5 MINUTES OF THE ZONE BOUNDARIES
!
!      THE OUTPUT IS FOR THE DATA SHEET PROGRAM
!
!      VARIABLES
!      CARDR = A MODIFIED 80 RECORD CARD WITH A LENGTH OF 211 COLS
!      ER = EQUATORIAL RADIUS OF THE ELLIPSOID (SEMI-MAJOR AXIS)
!      RF = RECIPROCAL OF FLATTING OF THE ELLIPSOD
!      ESQ= E SQUARED
!      RAD = RADIAN CONVERSION FACTOR
!      CM = CENTRAL MERIDIAN ( COMPUTED USING THE LONGITUDE)
!      SF = SCALE FACTOR OF CENTRAL MERIDIAN ( ALWAYS .9996 FOR UTM)
!      OR = SOUTHERNMOST PARALLEL OF LATITUDE ( ALWAYS ZERO FOR UTM)
!      R, A, B, C, U, V, W = ELLIPSOID CONSTANTS USED FOR COMPUTING
!                            MERIDIONAL DISTANCE FROM LATITUDE
!      SO = MERIDIONAL DISTANCE (MULTIPLIED BY SCALE FACTOR )
!           FROM THE EQUATOR TO THE SOUTHERNMOST PARALLEL OF LATITUDE
!           ( ALWAYS ZERO FOR UTM)
!

      !IMPLICIT DOUBLEPRECISION (A - H, O - Z)
        !Arguments-------------------------------------------------------------
        character(80)                               :: CARDR

        !Local-----------------------------------------------------------------
        character(1)                                :: EORW
        character(1)                                :: NORS
        character(4)                                :: ZONE
        integer                                     :: LD, LM, LOD, LOM, FOUND, IZ, ICM
        real(8)                                     :: SLAT, SLON, FI, LAM, LCM, UCM
        real(8)                                     :: KP
        real(8)                                     :: TD, TM, ND, NM
        real(8)                                     :: north, east
        real(8)                                     :: LOD1, R, RAD, CM, TOL, FN, FE, SF, OR
        real(8)                                     :: ER, RF, ESQ, EPS, A, B, C, U, V, W, SO
        real(8)                                     :: CONV, PI
        logical                                     :: MILP
        COMMON / MIL / MILP
        COMMON / CONST / RAD, ER, RF, ESQ, PI
        COMMON / XY / NORTH, EAST

      READ (CARDR, 50) LD, LM, SLAT, NORS, LOD, LOM, SLON, EORW
   50 FORMAT (T45,I2,I2,F10.5,A1,I3,I2,F10.5,A1)


!
!      CONVERT THE LATITUDE AND LONGITUDE TO PI AND LAM
!
      TD = DBLE (FLOAT (LD) )
      TM = DBLE (FLOAT (LM) )
      FI = (TD+ (TM + SLAT / 60.D0) / 60.D0) / RAD

      ND = DBLE (FLOAT (LOD) )
      NM = DBLE (FLOAT (LOM) )

      IF ( (EORW.EQ.'E') .OR. (EORW.EQ.'e') ) THEN
         LAM = (360.D0 - (ND+ (NM + SLON / 60.D0) / 60.D0) ) / RAD
         LOD1 = (360.D0 - (ND+ (NM + SLON / 60.D0) / 60.D0) )
         LOD = DINT (LOD1)
      ENDIF

      IF ( (EORW.EQ.'W') .OR. (EORW.EQ.'w') ) THEN
         LAM = (ND+ (NM + SLON / 60.D0) / 60.D0) / RAD
         LOD = LOD
      ENDIF



!
!      FIND THE ZONE FOR LONGITUDE LESS THAN 180 DEGREES
!
      IF (LOD.LT.180) THEN
         IZ = LOD / 6
         IZ = 30 - IZ

         ICM = (183 - (6 * IZ) )
         CM = DBLE (FLOAT (ICM) ) / RAD
         UCM = (ICM + 3) / RAD
         LCM = (ICM - 3) / RAD
      ENDIF
!
!       FIND THE ZONE FOR LONGITUDE GREATER THAN 180 DEGREES
!
      IF (LOD.GE.180) THEN
         IZ = (LOD) / 6
         IZ = 90 - IZ
         ICM = (543 - (6 * IZ) )
         CM = DBLE (FLOAT (ICM) ) / RAD
         UCM = (ICM + 3) / RAD
         LCM = (ICM - 3) / RAD
      ENDIF



      TOL = (5.0D0 / 60.0D0) / RAD

      FN = 0.D0

      IF ( (NORS.EQ.'S') .OR. (NORS.EQ.'s') ) THEN
         FN = 10000000.D0
      ENDIF

      IF ( (NORS.EQ.'N') .OR. (NORS.EQ.'n') ) THEN
         FN = 0.D0
      ENDIF


      FE = 500000.0D0
      SF = 0.9996D0
      OR = 0.0D0

      FOUND = 0

      CALL TCONST (ER, RF, SF, OR, ESQ, EPS, R, A, B, C, U, V, W, SO)

!
!      COMPUTE THE NORTH AND EASTINGS
!
                                          ! Long. of origin of Portugues
      IF (MILP) CM = DBLE (8.13190611) / RAD
                                          ! coordinates (8d 7' 54.862")
!
  200 CALL TMGRID (FI, LAM, NORTH, EAST, CONV, KP, ER, ESQ, EPS, CM, FE,&
      FN, SF, SO, R, A, B, C)
!
! --> Computes the transformations from UTM's to Portuguese Mil. coordin
!
      IF (MILP) THEN
         NORTH = NORTH / SF - 4092593.43D0
         EAST = EAST / SF - 300200.D0
      ENDIF
!
!      WRITE THE ZONE NUMBER
!
      IF (IZ.GT.9) THEN
         WRITE (ZONE, 600) IZ
  600 FORMAT   (1X,I2)
      ELSE
         WRITE (ZONE, 605) IZ
  605 FORMAT   (1X,I6.2)
      ENDIF

!
!      DO THE TEST TO SEE IF THE LONGITUDE IS WITHIN 5 MINUTES
!      OF THE BOUNDARIES FOR THE ZONE AND IF SO COMPUTE THE
!      NORHT AND EASTING FOR THE ADJACENT ZONE
!
      IF (FOUND.NE.0) THEN
         RETURN
      ENDIF

      IF (DABS (UCM - LAM) .LE.TOL) THEN
         CM = DBLE (FLOAT (ICM + 6) ) / RAD
         IZ = IZ - 1
         IF (IZ.EQ.0) IZ = 60
         FOUND = FOUND+1
         GOTO 200
      ENDIF

      IF (DABS (LCM - LAM) .LE.TOL) THEN
         CM = DBLE (FLOAT (ICM - 6) ) / RAD
         IZ = IZ + 1
         IF (IZ.EQ.61) IZ = 1
         FOUND = FOUND+1
         GOTO 200
      ENDIF


      RETURN
      END SUBROUTINE DRGPUT

!******************************************************************
      SUBROUTINE DRUTGP (ICODE, lat, lon)
!
!
!
!      THIS IS THE DRIVER TO COMPUTE LATITUDES AND LONGITUDES FROM
!      THE UTMs FOR EACH ZONE
!
!
!      VARIABLES
!      CARDR = A MODIFIED 80 RECORD CARD WITH A LENGTH OF 211 COLS
!      ER = EQUATORIAL RADIUS OF THE ELLIPSOID (SEMI-MAJOR AXIS)
!      RF = RECIPROCAL OF FLATTING OF THE ELLIPSOD
!      ESQ= E SQUARED
!      RAD = RADIAN CONVERSION FACTOR
!      CM = CENTRAL MERIDIAN ( COMPUTED USEING THE LONGITUDE)
!      SF = SCALE FACTOR OF CENTRAL MERIDIAN ( ALWAYS .9996 FOR UTM)
!      OR = SOUTHERNMOST PARALLEL OF LATITUDE ( ALWAYS ZERO FOR UTM)
!      R, A, B, C, U, V, W = ELLIPSOID CONSTANTS USED FOR COMPUTING
!                            MERIDIONAL DISTANCE FROM LATITUDE
!      SO = MERIDIONAL DISTANCE (MULTIPLIED BY SCALE FACTOR )
!           FROM THE EQUATOR TO THE SOUTHERNMOST PARALLEL OF LATITUDE
!           ( ALWAYS ZERO FOR UTM)
!

      !IMPLICIT DOUBLEPRECISION (A - H, O - Z)
        !Arguments-------------------------------------------------------------
        integer                                     :: ICODE
        real(8)                                     :: lat, lon

        !Local-----------------------------------------------------------------
        integer                                     :: IZ, ICM
        character(1)                                :: NORS
        character(4)                                :: ZONE
        integer                                     :: FOUND
        real(8)                                     :: LCM, UCM
        real(8)                                     :: KP
        real(8)                                     :: NORTH, EAST, R
        character(1)                                :: EORW
        logical                                     :: MILP
        real(8)                                     :: CM, RAD, TOL, FN, FE, SF, OR, CONV
        real(8)                                     :: EPS, SO, V0, V2, V4, V6, ER, ESQ, RF
        real(8)                                     :: PI
        COMMON / MIL / MILP
        COMMON / CONST / RAD, ER, RF, ESQ, PI
        COMMON / XY / NORTH, EAST


!
!      FIND THE CENTRAL MERIDAIN IF THE ZONE NUMBER IS LESS THAN 30
!
      IF (ICODE.LT.30) THEN
         IZ = ICODE
         ICM = (183 - (6 * IZ) )
         CM = DBLE (FLOAT (ICM) ) / RAD
         UCM = (ICM + 3) / RAD
         LCM = (ICM - 3) / RAD
      ENDIF
!
!       FIND THE CENTRAL MERIDAN IF THE ZONE NUMBER IS LARGER THAN 30
!
      IF (ICODE.GE.30) THEN
         IZ = ICODE
         ICM = (543 - (6 * IZ) )
         CM = DBLE (FLOAT (ICM) ) / RAD
         UCM = (ICM + 3) / RAD
         LCM = (ICM - 3) / RAD
      ENDIF

      TOL = (5.0D0 / 60.0D0) / RAD

      IF (NORTH.GT.10000000.0) THEN
         FN = 10000000.0D0
         NORS = 'S'
      ELSE
         FN = 0.D0
         NORS = 'N'
      ENDIF

      FE = 500000.0D0
      SF = 0.9996D0
      OR = 0.0D0

      FOUND = 0

      CALL TCONPC (SF, OR, EPS, R, SO, V0, V2, V4, V6, ER, ESQ, RF)

!
!      COMPUTE THE LATITUDES AND LONGITUDES
!
                                          ! Long. of origin of Portugues
      IF (MILP) CM = DBLE (8.13190611) / RAD
                                          ! coordinates (8d 7' 54.862")
!
!
! --> Computes the transformations from UTM's to Portuguese Mil. coordin
!
      IF (MILP) THEN
         NORTH = (NORTH + 4092593.43D0) * SF
         EAST = (EAST + 300200.D0) * SF
      ENDIF
  200 CALL TMGEOD (NORTH, EAST, LAT, LON, EPS, CM, FE, SF, SO, R, V0,   &
      V2, V4, V6, FN, ER, ESQ, CONV, KP)

!
!      WRITE THE ZONE NUMBER
!
      IF (IZ.GT.9) THEN
         WRITE (ZONE, 600) IZ
  600 FORMAT   (1X,I2)
      ELSE
         WRITE (ZONE, 605) IZ
  605 FORMAT   (1X,I6.2)
      ENDIF

      EORW = 'W'
      IF (nors.eq.'s'.or.nors.eq.'S') lat = - 1.d0 * lat
      IF (eorw.eq.'w'.or.eorw.eq.'W') lon = - 1.d0 * lon
!
!      DO THE TEST TO SEE IF THE LONGITUDE IS WITHIN 5 MINUTES
!      OF THE BOUNDARIES FOR THE ZONE AND IF SO COMPUTE THE
!      NORHT AND EASTING FOR THE ADJACENT ZONE
!
      IF (FOUND.NE.0) THEN
         RETURN
      ENDIF

      IF (DABS (UCM) .LE.TOL) THEN
         CM = DBLE (FLOAT (ICM + 6) ) / RAD
         IZ = IZ - 1
         IF (IZ.EQ.0) IZ = 60
         FOUND = FOUND+1
         GOTO 200
      ENDIF

      IF (DABS (LCM) .LE.TOL) THEN
         CM = DBLE (FLOAT (ICM - 6) ) / RAD
         IZ = IZ + 1
         IF (IZ.EQ.61) IZ = 1
         FOUND = FOUND+1
         GOTO 200
      ENDIF


      RETURN
      END SUBROUTINE DRUTGP

!********************************************************************
      SUBROUTINE TODMS (RAD, IDG, MIN, SEC)
!********************************************************************
!     RADIANS TO DEGREES,MINUTES AND SECONDS
!
      INTEGER IDG, MIN
      DOUBLEPRECISION RAD, SEC, RHOSEC
      DATA RHOSEC / 2.062648062471D05 /
      SEC = RAD * RHOSEC
      IDG = SEC / 3600.D0
      SEC = SEC - DBLE (IDG * 3600)
      MIN = SEC / 60.D0
      SEC = SEC - DBLE (MIN * 60)
      IF ( (60.D0 - DABS (SEC) ) .GT.5.D-6) GOTO 100
      SEC = SEC - DSIGN (60.D0, SEC)
      MIN = MIN + ISIGN (1, MIN)
  100 IF (IABS (MIN) .LT.60) GOTO 101
      MIN = MIN - ISIGN (60, MIN)
      IDG = IDG + ISIGN (1, IDG)
  101 MIN = IABS (MIN)
      SEC = DABS (SEC)
      IF (RAD.GE.0.D0) GOTO 102
      IF (IDG.EQ.0) MIN = - MIN
      IF (IDG.EQ.0.AND.MIN.EQ.0) SEC = - SEC
  102 RETURN
      END SUBROUTINE TODMS
!***********************************************************
      SUBROUTINE DATUMM (ER, RF, F, ESQ, DATNUM)

      CHARACTER(1) ANS, DATNUM
      DOUBLEPRECISION ER
      DOUBLEPRECISION RF
      DOUBLEPRECISION F
      DOUBLEPRECISION ESQ



!   50 PRINT *, '                                          '
!      PRINT *, '   WHICH ELLIPSOID DO YOU WANT ANSWER:    '
!      PRINT *, '   1.  CLARKE 1866                       '
!      PRINT *, '   2.  GRS80/WGS84                       '
!      PRINT *, '   3.  INTERNATIONAL 1910                '
!      PRINT *, '   4.  WGS72                             '
!      PRINT *, '   5.  OTHER ELLIPSOID                   '
!      PRINT *, '   TYPE NUMBER:  '
!      READ(*,FMT='(A1)') ANS
                !faz sempre para internacional
   50 ans = '3'
!*
!*     FIND THE RIGHT SEMI MAJOR AXIS AND FLATTING
!*
      IF (ANS.EQ.'1') THEN
!*      FOR THE NAD 27 DATUM
!*
         ER = 6378206.4D0
         RF = 294.978698D0
         F = 1.D0 / RF
         ESQ = (F + F - F * F)

      ELSEIF (ANS.EQ.'2') THEN
!*     FOR THE NAD83 DATUM
!*
         ER = 6378137.D0
         RF = 298.257222101D0
         F = 1.D0 / RF
         ESQ = (F + F - F * F)

      ELSEIF (ANS.EQ.'3') THEN
!*     FOR THE INT24 DATUM
!*
         ER = 6378388.D0
         RF = 297.0D0
         F = 1.D0 / RF
         ESQ = (F + F - F * F)

      ELSEIF (ANS.EQ.'4') THEN
!*     FOR THE WGS72
!*
         ER = 6378135.D0
         RF = 298.26D0
         F = 1.D0 / RF
         ESQ = (F + F - F * F)

      ELSEIF (ANS.EQ.'5') THEN
!*     FOR THE OTHER DATUM
!*
   10 PRINT * , '                    '
      PRINT * , '  SEMIMAJOR AXIS (meters)  '
      PRINT * , '  TYPE VALUE NOW:  '
         READ ( * , FMT = '(F12.0)') ER

         IF ( (ER.LE.6376400.D0) .OR. (ER.GT.6378500.D0) ) THEN
      PRINT * , '                  '
            PRINT * , ' SEMIMAJOR AXIS IS OUT OF RANGE - DO YOU'
            PRINT * , ' WANT TO TRY AGAIN '
            PRINT * , ' TYPE Y OR N '
            READ ( * , FMT = '(A1)') ANS
            IF ( (ANS.EQ.'Y') .OR. (ANS.EQ.'y') ) THEN
               GOTO 10
            ELSE
               GOTO 50
            ENDIF
         ENDIF

   20 PRINT * , '                  '
      PRINT * , ' FLATTENING            '
      PRINT * , ' TYPE VALUE NOW:  '
         READ ( * , FMT = '(F11.0)') RF

         IF ( (RF.LE.290.D0) .OR. (RF.GT.302.D0) ) THEN
      PRINT * , '                  '
            PRINT * , ' FLATTENING IS OUT OF RANGE - DO YOU '
            PRINT * , ' WANT TO TRY AGAIN '
            PRINT * , ' TYPE Y OR N '
            READ ( * , FMT = '(A1)') ANS
            IF ( (ANS.EQ.'Y') .OR. (ANS.EQ.'y') ) THEN
               GOTO 20
            ELSE
               GOTO 50
            ENDIF
         ENDIF

         F = 1.D0 / RF
         ESQ = (F + F - F * F)

      ELSE
      PRINT * , ' YOU TYPED THE INCORRECT NUMBER   '
         PRINT * , ' SO LET''S TRY AGAIN '
      PRINT * , '                  '
         GOTO 50
      ENDIF

      DATNUM = ANS

      RETURN
      END SUBROUTINE DATUMM

      SUBROUTINE TMGRID (FI, LAM, NORTH, EAST, CONV, KP, ER, ESQ, EPS,  &
                         CM, FE, FN, SF, SO, R, A, B, C)

        !IMPLICIT DOUBLEPRECISION (A - H, K - Z)
        !Arguments-------------------------------------------------------------
        real(8)                                     :: FI, LAM, NORTH, EAST
        real(8)                                     :: CONV, KP, ER, ESQ, EPS
        real(8)                                     :: CM, FE, FN, SF, SO, R
        real(8)                                     :: A, B, C

        !Local-----------------------------------------------------------------
        REAL(8)                                     :: A2, A1, A4, A6, A3, A5, A7
        REAL(8)                                     :: C1, C3, C5, F2, F4
        real(8)                                     :: OM, S, SINFI, COSFI, TN, TS
        real(8)                                     :: ETS, L, LS, RN
!
!*****  TRANSVERSE MERCATOR PROJECTION
!       CONVERSION OF GEODETIC COORDINATES TO GRID COORDINATES
!*****  Programmed by T. Vincenty, NGS, in July 1984.
!*****************  SYMBOLS AND DEFINITIONS *************************
!   Latitude positive north, longitude positive west.  All angles are
!     in radian measure.
!   N, E are northing and easting coordinates respectively.
!   LAT, LON are latitude and longitude respectively.
!   CONV is convergence.
!   KP is point scale factor.
!   ER is equatorial radius of the ellipsoid (= major semiaxis).
!   ESQ is the square of first eccentricity of the ellipsoid.
!   EPS is the square of second eccentricity of the ellipsoid.
!   CM is the central meridian of the projection zone.
!   FE is false easting value at the central meridian.
!   FN is "false northing" at the southernmost latitude, usually zero.
!   SF is scale factor at the central meridian.
!   SO is meridional distance (multiplied by the scale factor) from
!     the equator to the southernmost parallel of latitude for the zone.
!   R is the radius of the rectifying sphere (used for computing
!     meridional distance from latitude and vice versa).
!   A, B, C, U, V, W are other precomputed constants for determination
!     of meridional distance from latitude and vice versa.
!
!   The formula used in this subroutine gives geodetic accuracy within
!   zones of 7 degrees in east-west extent.  Within State transverse
!   Mercator projection zones, several minor terms of the equations
!   may be omitted (see a separate NGS publication).  If programmed
!   in full, the subroutine can be used for computations in surveys
!   extending over two zones.
!
!*********************************************************************
      OM = FI + A * SIN (2. * FI) + B * SIN (4. * FI) + C * SIN (6. *   &
      FI)
      S = R * OM * SF
      SINFI = SIN (FI)
      COSFI = COS (FI)
      TN = SINFI / COSFI
      TS = TN**2
      ETS = EPS * COSFI**2
      L = (LAM - CM) * COSFI
      LS = L * L
      RN = SF * ER / SQRT (1. - ESQ * SINFI**2)
!
      A2 = RN * TN / 2.
      A4 = (5. - TS + ETS * (9. + 4. * ETS) ) / 12.
      A6 = (61. + TS * (TS - 58.) + ETS * (270. - 330. * TS) ) / 360.
      A1 = - RN
      A3 = (1. - TS + ETS) / 6.
      A5 = (5. + TS * (TS - 18.) + ETS * (14. - 58. * TS) ) / 120.
      A7 = (61. - 479. * TS + 179. * TS**2 - TS**3) / 5040.
      NORTH = S - SO + A2 * LS * (1. + LS * (A4 + A6 * LS) ) - FN
      EAST = FE+A1 * L * (1. + LS * (A3 + LS * (A5 + A7 * LS) ) )

      IF (NORTH.LT.0.0) THEN
         NORTH = NORTH * ( - 1.0D0)
      ENDIF
!
!*** CONVERGENCE
      C1 = - TN
      C3 = (1. + 3. * ETS + 2. * ETS**2) / 3.
      C5 = (2. - TS) / 15.
      CONV = C1 * L * (1. + LS * (C3 + C5 * LS) )
!
!*** POINT SCALE FACTOR
      F2 = (1. + ETS) / 2.
      F4 = (5. - 4. * TS + ETS * (9. - 24. * TS) ) / 12.
      KP = SF * (1. + F2 * LS * (1. + F4 * LS) )
!
      RETURN
      END SUBROUTINE TMGRID
!********************************************************************
      SUBROUTINE TCONST (ER, RF, SF, OR, ESQ, EPS, R, A, B, C, U, V, W, &
                         SO)
        !IMPLICIT DOUBLEPRECISION (A - H, O - Z)
        !Arguments-------------------------------------------------------------
        real(8)                                     :: ER, RF, SF, OR, ESQ, EPS
        real(8)                                     :: R, A, B, C, U, V, W, SO

        !Local-----------------------------------------------------------------
        real(8)                                     :: F, PR, EN, OMO
!
!***** TRANSVERSE MERCATOR PROJECTION
!      PRECOMPUTATION OF CONSTANTS
!***** Programmed by T. Vincenty, NGS, in July 1984.
!******************** SYMBOLS AND DEFINITIONS  **********************
!   ER is equatorial radius of the ellipsoid (= major semiaxis).
!   RF is reciprocal of flattening of the ellipsoid.
!   SF is scale factor of the central meridian.
!   OR is southernmost parallel of latitude (in radians) for which
!     the northing coordinate is zero at the central meridian.
!   R, A, B, C, U, V, W are ellipsoid constants used for computing
!     meridional distance from latitude and vice versa.
!   SO is meridional distance (multiplied by the scale factor) from
!     the equator to the southernmost parallel of latitude.
!******************************************************************
!
      F = 1. / RF
      ESQ = (F + F - F**2)
      EPS = ESQ / (1. - ESQ)
      PR = (1. - F) * ER
      EN = (ER - PR) / (ER + PR)
      A = - 1.5D0 * EN + (9. / 16.) * EN**3
      B = 0.9375D0 * EN**2 - (15. / 32.) * EN**4
      C = - (35. / 48.) * EN**3
      U = 1.5D0 * EN - (27. / 32.) * EN**3
      V = 1.3125D0 * EN**2 - (55. / 32.) * EN**4
      W = (151. / 96.) * EN**3
      R = ER * (1. - EN) * (1. - EN**2) * (1. + 2.25D0 * EN**2 +        &
      (225. / 64.) * EN**4)
      OMO = OR + A * SIN (2. * OR) + B * SIN (4. * OR) + C * SIN (6. *  &
      OR)
      SO = SF * R * OMO
!
      RETURN
      END SUBROUTINE TCONST

      SUBROUTINE TCONPC (SF, OR, EPS, R, SO, V0, V2, V4, V6, ER, ESQ,   &
                         RF)

        !IMPLICIT DOUBLEPRECISION (A - H, O - Z)
        !Arguments-------------------------------------------------------------
        real(8)                                     :: SF, OR, EPS, R, SO, RF
        real(8)                                     :: V0, V2, V4, V6, ER, ESQ

        !Local-----------------------------------------------------------------
        real(8)                                     :: F, PR, EN, EN2, EN3, EN4
        real(8)                                     :: C2, C4, C6, C8, U0, U2, U4, U6
        real(8)                                     :: COSOR, OMO

!**          TRANSVERSE MERCATOR PROJECTION               ***
!** CONVERSION OF GRID COORDS TO GEODETIC COORDS
!** REVISED SUBROUTINE OF T. VINCENTY  FEB. 25, 1985
!************* SYMBOLS AND DEFINITIONS ***********************
!** ER IS THE SEMI-MAJOR AXIS FOR GRS-80
!** SF IS THE SCALE FACTOR AT THE CM
!** SO IS THE MERIDIONAL DISTANCE (TIMES THE SF) FROM THE
!**       EQUATOR TO SOUTHERNMOST PARALLEL OF LAT. FOR THE ZONE
!** R IS THE RADIUS OF THE RECTIFYING SPHERE
!** U0,U2,U4,U6,V0,V2,V4,V6 ARE PRECOMPUTED CONSTANTS FOR
!**   DETERMINATION OF MERIDIONAL DIST. FROM LATITUDE
!** OR IS THE SOUTHERNMOST PARALLEL OF LATITUDE FOR WHICH THE
!**       NORTHING COORD IS ZERO AT THE CM
!*************************************************************


      F = 1.D0 / RF
      EPS = ESQ / (1.D0 - ESQ)
      PR = (1.D0 - F) * ER
      EN = (ER - PR) / (ER + PR)
      EN2 = EN * EN
      EN3 = EN * EN * EN
      EN4 = EN2 * EN2

      C2 = - 3.D0 * EN / 2.D0 + 9.D0 * EN3 / 16.D0
      C4 = 15.D0 * EN2 / 16.D0 - 15.D0 * EN4 / 32.D0
      C6 = - 35.D0 * EN3 / 48.D0
      C8 = 315.D0 * EN4 / 512.D0
      U0 = 2.D0 * (C2 - 2.D0 * C4 + 3.D0 * C6 - 4.D0 * C8)
      U2 = 8.D0 * (C4 - 4.D0 * C6 + 10.D0 * C8)
      U4 = 32.D0 * (C6 - 6.D0 * C8)
      U6 = 128.D0 * C8

      C2 = 3.D0 * EN / 2.D0 - 27.D0 * EN3 / 32.D0
      C4 = 21.D0 * EN2 / 16.D0 - 55.D0 * EN4 / 32.D0
      C6 = 151.D0 * EN3 / 96.D0
      C8 = 1097.D0 * EN4 / 512.D0
      V0 = 2.D0 * (C2 - 2.D0 * C4 + 3.D0 * C6 - 4.D0 * C8)
      V2 = 8.D0 * (C4 - 4.D0 * C6 + 10.D0 * C8)
      V4 = 32.D0 * (C6 - 6.D0 * C8)
      V6 = 128.D0 * C8

      R = ER * (1.D0 - EN) * (1.D0 - EN * EN) * (1.D0 + 2.25D0 * EN *   &
      EN + (225.D0 / 64.D0) * EN4)
      COSOR = DCOS (OR)
      OMO = OR + DSIN (OR) * COSOR * (U0 + U2 * COSOR * COSOR + U4 *    &
      COSOR**4 + U6 * COSOR**6)
      SO = SF * R * OMO

      RETURN
      END SUBROUTINE TCONPC

      !------------------------------------------------------------------------

      SUBROUTINE TMGEOD (N, E, LAT, LON, EPS, CM, FE, SF, SO, R, V0, V2,&
                         V4, V6, FN, ER, ESQ, CONV, KP)

        !IMPLICIT DOUBLEPRECISION (A - H, K - Z)
        !Arguments-------------------------------------------------------------
        real(8)                                     :: N, E
        real(8)                                     :: LAT, LON, EPS, CM, FE
        real(8)                                     :: SF, SO, R, V0, V2, V4, V6
        real(8)                                     :: FN, ER, ESQ, CONV, KP

        !Local----------------------------------------------------------------
        real(8)                                     :: B1, B2, C1, TN, TS, ETS, RN, Q
        real(8)                                     :: OM, COSOM, FOOT, SINF, COSF
        real(8)                                     :: QS, B4, B6, B3, B5, B7, L, FI, LAM
        real(8)                                     :: SINFI, COSFI, L1, LS, C3, C5, F2, F4

!**          TRANSVERSE MERCATOR PROJECTION               ***
!** CONVERSION OF GRID COORDS TO GEODETIC COORDS
!** REVISED SUBROUTINE OF T. VINCENTY  FEB. 25, 1985
!************* SYMBOLS AND DEFINITIONS ***********************
!** LATITUDE POSITIVE NORTH, LONGITUDE POSITIVE WEST.  ALL
!**          ANGLES ARE IN RADIAN MEASURE.
!** LAT,LON ARE LAT. AND LONG. RESPECTIVELY
!** N,E ARE NORTHING AND EASTING COORDINATES RESPECTIVELY
!** K IS POINT SCALE FACTOR
!** ER IS THE SEMI-MAJOR AXIS OF THE ELLIPSOID
!** ESQ IS THE SQUARE OF THE 1ST ECCENTRICITY
!** E IS THE 1ST ECCENTRICITY
!** CM IS THE CENTRAL MERIDIAN OF THE PROJECTION ZONE
!** FE IS THE FALSE EASTING VALUE AT THE CM
!** CONV IS CONVERGENCE
!** EPS IS THE SQUARE OF THE 2ND ECCENTRICITY
!** SF IS THE SCALE FACTOR AT THE CM
!** SO IS THE MERIDIONAL DISTANCE (TIMES THE SF) FROM THE
!**       EQUATOR TO SOUTHERNMOST PARALLEL OF LAT. FOR THE ZONE
!** R IS THE RADIUS OF THE RECTIFYING SPHERE
!** U0,U2,U4,U6,V0,V2,V4,V6 ARE PRECOMPUTED CONSTANTS FOR
!**   DETERMINATION OF MERIDIANAL DIST. FROM LATITUDE
!**
!** THE FORMULA USED IN THIS SUBROUTINE GIVES GEODETIC ACCURACY
!** WITHIN ZONES OF 7 DEGREES IN EAST-WEST EXTENT.  WITHIN STATE
!** TRANSVERSE MERCATOR PROJECTION ZONES, SEVERAL MINOR TERMS OF
!** THE EQUATIONS MAY BE OMITTED (SEE A SEPARATE NGS PUBLICATION).
!** IF PROGRAMMED IN FULL, THE SUBROUTINE CAN BE USED FOR
!** COMPUTATIONS IN SURVEYS EXTENDING OVER TWO ZONES.
!**********************************************************************
      OM = (N - FN + SO) / (R * SF)
      COSOM = DCOS (OM)
      FOOT = OM + DSIN (OM) * COSOM * (V0 + V2 * COSOM * COSOM + V4 *   &
      COSOM**4 + V6 * COSOM**6)
      SINF = DSIN (FOOT)
      COSF = DCOS (FOOT)
      TN = SINF / COSF
      TS = TN * TN
      ETS = EPS * COSF * COSF
      RN = ER * SF / DSQRT (1.D0 - ESQ * SINF * SINF)
      Q = (E-FE) / RN
      QS = Q * Q
      B2 = - TN * (1.D0 + ETS) / 2.D0
      B4 = - (5.D0 + 3.D0 * TS + ETS * (1.D0 - 9.D0 * TS) - 4.D0 * ETS *&
      ETS) / 12.D0
      B6 = (61.D0 + 45.D0 * TS * (2.D0 + TS) + ETS * (46.D0 - 252.D0 *  &
      TS - 60.D0 * TS * TS) ) / 360.D0
      B1 = 1.D0
      B3 = - (1.D0 + TS + TS + ETS) / 6.D0
      B5 = (5.D0 + TS * (28.D0 + 24.D0 * TS) + ETS * (6.D0 + 8.D0 * TS) &
      ) / 120.D0
      B7 = - (61.D0 + 662.D0 * TS + 1320.D0 * TS * TS + 720.D0 * TS**3) &
      / 5040.D0
      LAT = FOOT + B2 * QS * (1.D0 + QS * (B4 + B6 * QS) )
      L = B1 * Q * (1.D0 + QS * (B3 + QS * (B5 + B7 * QS) ) )
      LON = - L / COSF + CM

!*********************************************************************
!     COMPUTE CONVERENCE AND SCALE FACTOR
      FI = LAT
      LAM = LON
      SINFI = SIN (FI)
      COSFI = COS (FI)
      L1 = (LAM - CM) * COSFI
      LS = L1 * L1

      TN = SINFI / COSFI
      TS = TN * TN
!
!*** CONVERGENCE
      C1 = - TN
      C3 = (1. + 3. * ETS + 2. * ETS**2) / 3.
      C5 = (2. - TS) / 15.
      CONV = C1 * L1 * (1. + LS * (C3 + C5 * LS) )
!
!*** POINT SCALE FACTOR
      F2 = (1. + ETS) / 2.
      F4 = (5. - 4. * TS + ETS * (9. - 24. * TS) ) / 12.
      KP = SF * (1. + F2 * LS * (1. + F4 * LS) )

      RETURN
      END SUBROUTINE TMGEOD

    !--------------------------------------------------------------------------
!-------------------------------------------------------------------------

!*************************************************************************
!                           SUBROTINA INTERF                             *
!*************************************************************************

    Subroutine Interf2D (FatherGrid, X, Y, PP, ILB, JLB, IUB, JUB, Abs1, Ord1, F)

        !Arguments-------------------------------------------------------------

        type (T_FatherGrid),     pointer  :: FatherGrid
        real, dimension (:,:),   pointer  :: F
        real, dimension (:  ),   pointer  :: Abs1, Ord1
        real                              :: X, Y, PP
        integer                           :: IUB, JUB, ILB, JLB

        !Local------------------------------------------------------------
        integer                           :: JX, IY, J1, J2, I1, I2
        real                              :: Fhc, RJ, PASX, RI, PASY, DF, FJ2, FJ1

        !Begin------------------------------------------------------------

        JX  = FatherGrid%JX
        IY  = FatherGrid%IY
        Fhc = FatherGrid%Fhc


      100 Continue
          IF (X.EQ.ABS1(JUB)) Then
             JX = JUB
             J1 = JUB-1
             J2 = JUB
             Goto 200
          Endif
          IF (X.LT.ABS1(JX)) Goto 110
          JX = JX+1
          IF (JX.GE.JUB) Goto 120
          Goto 100
      110 Continue
          JX = JX-1
          IF (JX.LT.JLB) Goto 330
          IF (X.LT.ABS1(JX)) Goto 110
          J1 = JX
          J2 = JX+1
          Goto 200
      120 Continue
          IF (JX.NE.JUB) Goto 330
          J1 = JX-1
          J2 = JX
      200 Continue
          IF (Y.EQ.ORD1(IUB)) Then
             IY = IUB
             I1 = IUB-1
             I2 = IUB
             Goto 300
          Endif
          IF (Y.LT.ORD1(IY)) Goto 210
          IY = IY+1
          IF (IY.GE.IUB) Goto 220
          Goto 200
      210 Continue
          IY = IY-1
          IF (IY.LT.ILB) Goto 330
          IF (Y.LT.ORD1(IY)) Goto 210
          I1 = IY
          I2 = IY+1
          Goto 300
      220 Continue
          IF (IY.NE.IUB) Goto 330
          I1 = IY-1
          I2 = IY
      300 Continue

          RJ   = X-ABS1(J1)
          PASX = ABS1(J2)-ABS1(J1)
          RI   = Y-ORD1(I1)
          PASY = ORD1(I2)-ORD1(I1)
          DF   = (F(I2,J2)-F(I1,J2))/PASY

          FJ2  = F(I1,J2)+DF*RI
          DF   = (F(I2,J1)-F(I1,J1))/PASY
          FJ1  = F(I1,J1)+DF*RI
          DF   = (FJ2-FJ1)/PASX
          PP   = FJ1+RJ*DF
          IF (JX.GT.JUB) JX=JUB
          IF (IY.GT.IUB) IY=IUB

          RETURN

      330 Continue
          PP = FHC
          IF (JX.LT.JLB) JX=JLB
          IF (JX.GT.JUB) JX=JUB
          IF (IY.LT.ILB) IY=ILB
          IF (IY.GT.IUB) IY=IUB

    end Subroutine Interf2D

!------------------------------------------------------------------------------

    Subroutine Interf3D (FatherGrid, X, Y, PP, ILB, JLB, IUB, JUB, K, Abs1, Ord1, F)

        !Arguments-------------------------------------------------------------

        type (T_FatherGrid),     pointer  :: FatherGrid
        real, dimension (:,:,:), pointer  :: F
        real, dimension (:  ),   pointer  :: Abs1, Ord1
        real                              :: X, Y, PP
        integer                           :: IUB, JUB, ILB, JLB, K

        !Local------------------------------------------------------------
        integer                           :: JX, IY, J1, J2, I1, I2
        real                              :: Fhc, RJ, PASX, RI, PASY, DF, FJ2, FJ1

        !Begin------------------------------------------------------------

        JX  = FatherGrid%JX
        IY  = FatherGrid%IY
        Fhc = FatherGrid%Fhc


      100 Continue
          IF (X.EQ.ABS1(JUB)) Then
             JX = JUB
             J1 = JUB-1
             J2 = JUB
             Goto 200
          Endif
          IF (X.LT.ABS1(JX)) Goto 110
          JX = JX+1
          IF (JX.GE.JUB) Goto 120
          Goto 100
      110 Continue
          JX = JX-1
          IF (JX.LT.JLB) Goto 330
          IF (X.LT.ABS1(JX)) Goto 110
          J1 = JX
          J2 = JX+1
          Goto 200
      120 Continue
          IF (JX.NE.JUB) Goto 330
          J1 = JX-1
          J2 = JX
      200 Continue
          IF (Y.EQ.ORD1(IUB)) Then
             IY = IUB
             I1 = IUB-1
             I2 = IUB
             Goto 300
          Endif
          IF (Y.LT.ORD1(IY)) Goto 210
          IY = IY+1
          IF (IY.GE.IUB) Goto 220
          Goto 200
      210 Continue
          IY = IY-1
          IF (IY.LT.ILB) Goto 330
          IF (Y.LT.ORD1(IY)) Goto 210
          I1 = IY
          I2 = IY+1
          Goto 300
      220 Continue
          IF (IY.NE.IUB) Goto 330
          I1 = IY-1
          I2 = IY
      300 Continue

          RJ   = X-ABS1(J1)
          PASX = ABS1(J2)-ABS1(J1)
          RI   = Y-ORD1(I1)
          PASY = ORD1(I2)-ORD1(I1)
          DF   = (F(I2,J2, K)-F(I1,J2, K))/PASY

          FJ2  = F(I1,J2,K)+DF*RI
          DF   = (F(I2,J1,K)-F(I1,J1,K))/PASY
          FJ1  = F(I1,J1,K)+DF*RI
          DF   = (FJ2-FJ1)/PASX
          PP   = FJ1+RJ*DF
          IF (JX.GT.JUB) JX=JUB
          IF (IY.GT.IUB) IY=IUB

          return

      330 Continue
          PP = FHC
          IF (JX.LT.JLB) JX=JLB
          IF (JX.GT.JUB) JX=JUB
          IF (IY.LT.ILB) IY=ILB
          IF (IY.GT.IUB) IY=IUB



    end Subroutine Interf3D

!------------------------------------------------------------------------------

    Subroutine Interf3D8 (FatherGrid, X, Y, PP, ILB, JLB, IUB, JUB, K, Abs1, Ord1, F)

        !Arguments-------------------------------------------------------------

        type (T_FatherGrid),        pointer  :: FatherGrid
        real(8), dimension (:,:,:), pointer  :: F
        real, dimension (:  ),      pointer  :: Abs1, Ord1
        real                                 :: X, Y
        real(8)                              :: PP
        integer                              :: IUB, JUB, ILB, JLB, K

        !Local------------------------------------------------------------
        integer                              :: JX, IY, J1, J2, I1, I2
        real                                 :: Fhc, RJ, PASX, RI, PASY, DF, FJ2, FJ1

        !Begin------------------------------------------------------------

        JX  = FatherGrid%JX
        IY  = FatherGrid%IY
        Fhc = FatherGrid%Fhc


      100 Continue
          IF (X.EQ.ABS1(JUB)) Then
             JX = JUB
             J1 = JUB-1
             J2 = JUB
             Goto 200
          Endif
          IF (X.LT.ABS1(JX)) Goto 110
          JX = JX+1
          IF (JX.GE.JUB) Goto 120
          Goto 100
      110 Continue
          JX = JX-1
          IF (JX.LT.JLB) Goto 330
          IF (X.LT.ABS1(JX)) Goto 110
          J1 = JX
          J2 = JX+1
          Goto 200
      120 Continue
          IF (JX.NE.JUB) Goto 330
          J1 = JX-1
          J2 = JX
      200 Continue
          IF (Y.EQ.ORD1(IUB)) Then
             IY = IUB
             I1 = IUB-1
             I2 = IUB
             Goto 300
          Endif
          IF (Y.LT.ORD1(IY)) Goto 210
          IY = IY+1
          IF (IY.GE.IUB) Goto 220
          Goto 200
      210 Continue
          IY = IY-1
          IF (IY.LT.ILB) Goto 330
          IF (Y.LT.ORD1(IY)) Goto 210
          I1 = IY
          I2 = IY+1
          Goto 300
      220 Continue
          IF (IY.NE.IUB) Goto 330
          I1 = IY-1
          I2 = IY
      300 Continue

          RJ   = X-ABS1(J1)
          PASX = ABS1(J2)-ABS1(J1)
          RI   = Y-ORD1(I1)
          PASY = ORD1(I2)-ORD1(I1)
          DF   = (F(I2,J2, K)-F(I1,J2, K))/PASY

          FJ2  = F(I1,J2,K)+DF*RI
          DF   = (F(I2,J1,K)-F(I1,J1,K))/PASY
          FJ1  = F(I1,J1,K)+DF*RI
          DF   = (FJ2-FJ1)/PASX
          PP   = FJ1+RJ*DF
          IF (JX.GT.JUB) JX=JUB
          IF (IY.GT.IUB) IY=IUB

          return

      330 Continue
          PP = FHC
          IF (JX.LT.JLB) JX=JLB
          IF (JX.GT.JUB) JX=JUB
          IF (IY.LT.ILB) IY=ILB
          IF (IY.GT.IUB) IY=IUB



    end Subroutine Interf3D8


!------------------------------------------------------------------------------

    Subroutine LocateCell  (XX, YY, XPos, YPos,                                         &
                            ILB, IUB, JLB, JUB, ILower, Jlower)


        !Arguments---------------------------------------------------------
        real,   dimension(:)  , pointer      :: XX, YY
        real                  , intent (IN ) :: XPos, YPos
        integer               , intent (IN ) :: ILB, IUB, JLB, JUB
        integer               , intent (OUT) :: ILower, Jlower


        !Local-------------------------------------------------------------
        integer                              :: IMiddle, JMiddle, Dcd1, Dcd2, IUpper, JUpper

        !Begin-------------------------------------------------------------


        ILower = ILB
        IUpper = IUB

        Dcd1 = 2
        Dcd2 = 2

        if(Ypos < YY(ILB) .or. Ypos > YY(IUB))then

            ILower = null_int

        else

            do while (Dcd1 > 1)
                IMiddle = int((ILower + IUpper)/2)
                if (Ypos > YY(IMiddle)) then
                    ILower = IMiddle
                else
                    IUpper = IMiddle
                endif
                Dcd1 = IUpper - ILower
            enddo

        end if


        Jlower = JLB
        JUpper = JUB

        if(Xpos < XX(JLB) .or. Xpos > XX(JUB))then

            Jlower = null_int

        else

            do while (Dcd2 > 1)
                JMiddle = int((Jlower + JUpper)/2)
                if (Xpos > XX(JMiddle)) then
                    Jlower = JMiddle
                else
                    JUpper = JMiddle
                endif
                Dcd2 = JUpper - Jlower
            enddo

        end if

        if (.not. Dcd1>0 .and. Dcd2>0) then

            call SetError(FATAL_, INTERNAL_, 'LocateCell - HorizontalGrid - ERR01')


        endif

    end subroutine LocateCell


!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

    Subroutine LocateCell1DR4 (XX, XPos, JLB, JUB, Jlower)


        !Arguments---------------------------------------------------------
        real(4), dimension(:)  , pointer      :: XX
        real(4)                , intent (IN ) :: XPos
        integer                , intent (IN ) :: JLB, JUB
        integer                , intent (OUT) :: Jlower


        !Local-------------------------------------------------------------
        integer                              :: JMiddle, Dcd2, JUpper

        !Begin-------------------------------------------------------------

        Dcd2 = 2

        Jlower = JLB
        JUpper = JUB

        if(Xpos <= XX(JLB)) then

            Jlower = JLB

        elseif (Xpos >= XX(JUB))then

            Jlower = JUB

        else

            do while (Dcd2 > 1)
                JMiddle = int((Jlower + JUpper)/2)
                if (Xpos > XX(JMiddle)) then
                    Jlower = JMiddle
                else
                    JUpper = JMiddle
                endif
                Dcd2 = JUpper - Jlower
            enddo

        end if

        if (.not. Dcd2>0) then

            call SetError(FATAL_, INTERNAL_, 'LocateCell1DR4 - HorizontalGrid - ERR01')


        endif

    end subroutine LocateCell1DR4


!------------------------------------------------------------------------------

    Subroutine LocateCell1DR8 (XX, XPos, JLB, JUB, Jlower)


        !Arguments---------------------------------------------------------
        real(8), dimension(:)  , pointer      :: XX
        real(8)                , intent (IN ) :: XPos
        integer                , intent (IN ) :: JLB, JUB
        integer                , intent (OUT) :: Jlower


        !Local-------------------------------------------------------------
        integer                              :: JMiddle, Dcd2, JUpper

        !Begin-------------------------------------------------------------

        Dcd2 = 2

        Jlower = JLB
        JUpper = JUB

        if(Xpos <= XX(JLB)) then

            Jlower = JLB

        elseif (Xpos >= XX(JUB))then

            Jlower = JUB

        else

            do while (Dcd2 > 1)
                JMiddle = int((Jlower + JUpper)/2)
                if (Xpos > XX(JMiddle)) then
                    Jlower = JMiddle
                else
                    JUpper = JMiddle
                endif
                Dcd2 = JUpper - Jlower
            enddo

        end if

        if (.not. Dcd2>0) then

            call SetError(FATAL_, INTERNAL_, 'LocateCell1DR8 - HorizontalGrid - ERR01')


        endif

    end subroutine LocateCell1DR8


!------------------------------------------------------------------------------
    Subroutine LocateCell_2  (XX, YY, XPos, YPos,                                         &
                            ILB, IUB, JLB, JUB, ILink, JLink, U, V)

        !Arguments---------------------------------------------------------
        real,   dimension(:)  , pointer      :: XX, YY
        real                  , intent (IN ) :: XPos, YPos
        integer               , intent (IN ) :: ILB, IUB, JLB, JUB
        integer               , intent (OUT) :: ILink, JLink
        integer, optional                    :: U, V

        !Local-------------------------------------------------------------
        integer                              :: IMiddle, JMiddle, Dcd1, Dcd2, IUpper, JUpper

        !Begin-------------------------------------------------------------
        ILink  = ILB
        IUpper = IUB

        Dcd1 = 2
        Dcd2 = 2

        if(Ypos < YY(ILB) .or. Ypos > YY(IUB))then

            ILink = null_int

        elseif (present(V))then
            do while (Dcd1 > 1)
                IMiddle = int((ILink + IUpper)/2)
                if (Ypos > YY(IMiddle - 1)) then
                    ILink = IMiddle
                else
                    IUpper = IMiddle
                endif
                Dcd1 = IUpper - ILink
            enddo
        else
            do while (Dcd1 > 1)
                IMiddle = int((ILink + IUpper)/2)
                if (Ypos > YY(IMiddle)) then
                    ILink = IMiddle
                else
                    IUpper = IMiddle
                endif
                Dcd1 = IUpper - ILink
            enddo

        end if


        JLink = JLB
        JUpper = JUB

        if(Xpos < XX(JLB) .or. Xpos > XX(JUB))then

            JLink = null_int

        elseif (present(U))then
            do while (Dcd2 > 1)
                JMiddle = int((JLink + JUpper)/2)
                if (Xpos > XX(JMiddle - 1)) then
                    JLink = JMiddle
                else
                    JUpper = JMiddle
                endif
                Dcd2 = JUpper - JLink
            enddo
        else

            do while (Dcd2 > 1)
                JMiddle = int((JLink + JUpper)/2)
                if (Xpos > XX(JMiddle)) then
                    JLink = JMiddle
                else
                    JUpper = JMiddle
                endif
                Dcd2 = JUpper - JLink
            enddo

        end if

        if (.not. Dcd1>0 .and. Dcd2>0) then

            call SetError(FATAL_, INTERNAL_, 'LocateCell_2 - HorizontalGrid - ERR01')


        endif


        end subroutine LocateCell_2
!------------------------------------------------------------------------------

    real  function Bilinear (YYUpper, YYLower, YPos,                                     &
                            XXUpper, XXLower, XPos,                                      &
                            PropLowLeft,  PropUpLeft, PropLowRight, PropUpRight,         &
                            ONLowLeft,    ONUpLeft,   ONLowRight,   ONUpRight)

        !Arguments---------------------------------------------------------
        real    :: YYUpper, YYLower, YPos
        real    :: XXUpper, XXLower, XPos
        real    :: PropLowLeft, PropUpLeft, PropLowRight, PropUpRight
        integer :: ONLowLeft,   ONUpLeft,   ONLowRight,   ONUpRight

        !Local-------------------------------------------------------------
        real    :: Prop1, Prop2, Dx3, Dx2, Dy2, Dy3, D1, D2, D3
        integer :: ON1, ON2

        !Begin-------------------------------------------------------------

        Dy2   = YYUpper - YPos
        Dy3   = YPos    - YYLower

        Dx2   = XXUpper - XPos
        Dx3   = XPos    - XXLower

        D1    =  Dy2 * ONLowLeft  + Dy3 * ONUpLeft
        D2    =  Dy2 * ONLowRight + Dy3 * ONUpRight

cd1:    if (D1 /= 0) then

            Prop1 = (Dy2 * ONLowLeft * PropLowLeft  + Dy3 * ONUpLeft * PropUpLeft) /  D1

        else if (ONLowLeft == 1) then cd1

            Prop1 = PropLowLeft

        else if (ONUpLeft  == 1) then cd1

            Prop1 = PropUpLeft

        else  cd1

            Prop1 = 0.

        endif cd1

cd2:    if (D2 /= 0) then

            Prop2 = (Dy2 * ONLowRight * PropLowRight  + Dy3 * ONUpRight * PropUpRight) /  D2

        else if (ONLowRight == 1) then cd2

            Prop2 = PropLowRight

        else if (ONUpRight  == 1) then cd2

            Prop2 = PropUpRight

        else  cd2

            Prop2 = 0.

        endif cd2


        ON1   = ONUpLeft   + ONLowLeft

        if (ON1 > 1) ON1 = 1

        ON2   = ONLowRight + ONUpRight

        if (ON2 > 1) ON2 = 1

        D3 = Dx2 * ON1 + Dx3 * ON2


cd3:    if (D3 /= 0) then

            Bilinear  = (Dx2 * ON1 * Prop1 + Dx3 * ON2 * Prop2) / D3

        else if (ON1 == 1 ) then cd3

            Bilinear = Prop1

        else if (ON2  == 1) then cd3

            Bilinear = Prop2

        else  cd3

            Bilinear = 0.

        endif cd3



    end function Bilinear

!------------------------------------------------------------------------------

    subroutine InterpolPoint(XXSon, YYSon,                                              &
                             YYUpper, YYLower, XXUpper, XXLower,                        &
                             PropLowLeft,  PropUpLeft, PropLowRight, PropUpRight,       &
                             ONLowLeft,    ONUpLeft,   ONLowRight,   ONUpRight,         &
                             PropSon, OK)

        !Arguments ---------------------------------------------------------
        real   ,                 intent(IN)  :: XXSon, YYSon,                           &
                                                YYUpper, YYLower, XXUpper, XXLower,     &
                                                PropLowLeft,  PropUpLeft, PropLowRight, &
                                                PropUpRight
        integer,                 intent(IN ) :: ONLowLeft, ONUpLeft, ONLowRight,        &
                                                ONUpRight
        real   ,                 intent(OUT) :: PropSon
        logical                , intent(OUT) :: OK

        !Local -------------------------------------------------------------
        integer                            :: SumON

        !Begin -------------------------------------------------------------

        SumON = ONLowLeft + ONUpLeft + ONLowRight + ONUpRight

cd1:    if      (SumON > 0) then

            PropSon = Bilinear(YYUpper, YYLower, YYSon,                                 &
                               XXUpper, XXLower, XXSon,                                 &
                               PropLowLeft,  PropUpLeft, PropLowRight, PropUpRight,     &
                               ONLowLeft,    ONUpLeft,   ONLowRight,   ONUpRight)

            OK = .True.

        else cd1

            OK = .False.

        endif cd1


    end subroutine InterpolPoint

   !------------------------------------------------------------------------

    subroutine LocateCellPolygons(XX2D_Z, YY2D_Z, XPoint, YPoint,                       &
                                  DefinedPoint, ILB, IUB, JLB, JUB, IZ, JZ,             &
                                  CellLocated, Iold, Jold)

        !Arguments ---------------------------------------------------------
        real,    dimension(:,:)  , pointer              :: XX2D_Z, YY2D_Z
        integer, dimension(:,:)  , pointer              :: DefinedPoint
        real   ,                 intent(IN )            :: XPoint, YPoint
        integer,                 intent(IN )            :: ILB, IUB, JLB, JUB
        integer,                 intent(OUT)            :: IZ, JZ
        logical,                 intent(OUT), optional  :: CellLocated
        integer,                 intent(IN ), optional  :: Iold, Jold

        !Local -------------------------------------------------------------
        integer                                         :: ICenter, JCenter, Iupper, ILower,&
                                                           Jupper, JLower
        type(T_PointF ),          pointer               :: Point
        type(T_Polygon),          pointer               :: Polygon
        integer                                         :: i, j, pi
        integer                                         :: I1, I2, I3, I4
        integer                                         :: J1, J2, J3, J4
        logical                                         :: IsPointInside, SearchCell, CellFound
        logical                                         :: UnknownInitialIJ

        !Begin -------------------------------------------------------------

        CellFound        = .false.
        SearchCell       = .true.
        UnknownInitialIJ = .true.

        allocate(Point)
        Point%X = XPoint
        Point%Y = YPoint

        !check if
f1:     if (present(Iold) .and. present(Jold)) then

f12:        if (Iold>=ILB .and. Iold<=IUB .and. Jold>=JLB .and. Jold<=JUB) then

                UnknownInitialIJ = .false.

                Me%AuxPolygon%VerticesF(1)%X  = XX2D_Z(Iold  ,Jold  )
                Me%AuxPolygon%VerticesF(1)%Y  = YY2D_Z(Iold  ,Jold  )

                Me%AuxPolygon%VerticesF(2)%X  = XX2D_Z(Iold+1,Jold  )
                Me%AuxPolygon%VerticesF(2)%Y  = YY2D_Z(Iold+1,Jold  )

                Me%AuxPolygon%VerticesF(3)%X  = XX2D_Z(Iold+1,Jold+1)
                Me%AuxPolygon%VerticesF(3)%Y  = YY2D_Z(Iold+1,Jold+1)

                Me%AuxPolygon%VerticesF(4)%X  = XX2D_Z(Iold  ,Jold+1)
                Me%AuxPolygon%VerticesF(4)%Y  = YY2D_Z(Iold  ,Jold+1)

                !close polygon
                Me%AuxPolygon%VerticesF(5)%X  = Me%AuxPolygon%VerticesF(1)%X
                Me%AuxPolygon%VerticesF(5)%Y  = Me%AuxPolygon%VerticesF(1)%Y

                call SetLimits(Me%AuxPolygon)

                if (IsPointInsidePolygon(Point, Me%AuxPolygon)) then

                    SearchCell = .false.
                    CellFound  = .true.

                    IZ         = Iold
                    JZ         = Jold

                endif

            else f12

                UnknownInitialIJ = .true.

            endif f12

        endif f1

f2:     if (SearchCell) then

            !search adjacent cells
f3:         if (present(Iold) .and. present(Jold) .and. (.not. UnknownInitialIJ)) then

                ILower = Iold - 1
                IUpper = Iold + 1
                JLower = Jold - 1
                JUpper = Jold + 1

                if(JLower .lt. JLB)then
                    JLower = JLB
                end if

                if(JUpper .gt. JUB)then
                    JUpper = JUB
                end if

                if(ILower .lt. ILB)then
                    ILower = ILB
                end if

                if(IUpper .gt. IUB)then
                    IUpper = IUB
                end if

                do j = JLower, JUpper
                do i = ILower, IUpper

                    if(i == iold .and. j == jold) cycle

                    if(DefinedPoint(i,j) == 1)then

                        Me%AuxPolygon%VerticesF(1)%X  = XX2D_Z(i  ,j  )
                        Me%AuxPolygon%VerticesF(1)%Y  = YY2D_Z(i  ,j  )

                        Me%AuxPolygon%VerticesF(2)%X  = XX2D_Z(i+1,j  )
                        Me%AuxPolygon%VerticesF(2)%Y  = YY2D_Z(i+1,j  )

                        Me%AuxPolygon%VerticesF(3)%X  = XX2D_Z(i+1,j+1)
                        Me%AuxPolygon%VerticesF(3)%Y  = YY2D_Z(i+1,j+1)

                        Me%AuxPolygon%VerticesF(4)%X  = XX2D_Z(i  ,j+1)
                        Me%AuxPolygon%VerticesF(4)%Y  = YY2D_Z(i  ,j+1)

                        !close polygon
                        Me%AuxPolygon%VerticesF(5)%X  = Me%AuxPolygon%VerticesF(1)%X
                        Me%AuxPolygon%VerticesF(5)%Y  = Me%AuxPolygon%VerticesF(1)%Y

                        call SetLimits(Me%AuxPolygon)

                        if (IsPointInsidePolygon(Point, Me%AuxPolygon)) then

                            SearchCell = .false.
                            CellFound  = .true.

                            IZ         = i
                            JZ         = j

                            exit

                        endif

                    endif

                enddo
                if(CellFound) exit  !exit do J loop
                enddo

            end if f3

            !if still has not found the grid cell
f4:         if(SearchCell)then

                Iupper = IUB
                ILower = ILB

                Jupper = JUB
                JLower = JLB

                allocate(Polygon)
                allocate(Polygon%VerticesF(1: 2*(IUB-ILB+1+JUB-JLB+1)))


                do while (SearchCell)

                    ICenter = int(real((Iupper - ILower)/ 2.)) + ILower
                    JCenter = int(real((Jupper - JLower)/ 2.)) + JLower

                    IsPointInside = .false.

                    !Construct 4 polygons NW, NE, SW, SE
                    do pi = 1, 4

                        !All polygons are defined anti-Clockwise
                        if      (pi==1) then

                            if (ILower == ICenter) cycle
                            if (JLower == JCenter) cycle

                            I1 = ILower;  J1 = JLower
                            I2 = ILower;  J2 = JCenter
                            I3 = ICenter; J3 = JCenter
                            I4 = ICenter; J4 = JLower

                        else if (pi==2) then

                            if (ILower == ICenter) cycle

                            I1 = ILower;  J1 = JCenter;
                            I2 = ILower;  J2 = Jupper;
                            I3 = ICenter; J3 = Jupper;
                            I4 = ICenter; J4 = JCenter;

                        else if (pi==3) then

                            if (JLower == JCenter) cycle

                            I1 = ICenter; J1 = JLower;
                            I2 = ICenter; J2 = JCenter;
                            I3 = IUpper ; J3 = JCenter;
                            I4 = IUpper;  J4 = JLower;

                        else if (pi==4) then


                            I1 = ICenter; J1 = JCenter;
                            I2 = ICenter; J2 = JUpper;
                            I3 = IUpper ; J3 = JUpper;
                            I4 = IUpper;  J4 = JCenter;

                        endif

                        Polygon%Count = 0.

                        do j = J1, J2

                            if (DefinedPoint(I1,j) == 1) then

                                Polygon%Count = Polygon%Count + 1

                                Polygon%VerticesF(Polygon%Count)%X  = XX2D_Z(I1,j)
                                Polygon%VerticesF(Polygon%Count)%Y  = YY2D_Z(I1,j)

                            endif

                        end do

                        do i = I2, I3

                            if (DefinedPoint(i,J2) == 1) then

                                Polygon%Count = Polygon%Count + 1

                                Polygon%VerticesF(Polygon%Count)%X  = XX2D_Z(i,J2)
                                Polygon%VerticesF(Polygon%Count)%Y  = YY2D_Z(i,J2)

                            endif

                        end do

                        do j = J3, J4, -1

                            if (DefinedPoint(I3,j) == 1) then

                                Polygon%Count = Polygon%Count + 1

                                Polygon%VerticesF(Polygon%Count)%X  = XX2D_Z(I3,j)
                                Polygon%VerticesF(Polygon%Count)%Y  = YY2D_Z(I3,j)

                            endif

                        end do

                        do i = I4, I1, -1

                            if (DefinedPoint(i,J4) == 1) then

                                Polygon%Count = Polygon%Count + 1

                                Polygon%VerticesF(Polygon%Count)%X  = XX2D_Z(i,J4)
                                Polygon%VerticesF(Polygon%Count)%Y  = YY2D_Z(i,J4)

                            endif

                        end do

                        Polygon%Count = Polygon%Count + 1

                        !Close polygon
                        Polygon%VerticesF(Polygon%Count)%X  = Polygon%VerticesF(1)%X
                        Polygon%VerticesF(Polygon%Count)%Y  = Polygon%VerticesF(1)%Y

                        call SetLimits(Polygon)

                        IsPointInside = IsPointInsidePolygon(Point, Polygon)

                        if (IsPointInside) exit

                    enddo

                    Iupper  = I3
                    ILower  = I1

                    Jupper  = J2
                    JLower  = J1

                    if ((Iupper - ILower) == 1 .and. (Jupper - JLower) == 1) then

                        SearchCell = .false.

                        if (IsPointInside) then
                            !Test if the cell was found
                             CellFound = .true.
                        endif
                    endif

                end do

                IZ = ILower
                JZ = JLower

                deallocate(Polygon%VerticesF)
                nullify   (Polygon%VerticesF)

                deallocate(Polygon)
                nullify   (Polygon)

            endif f4

        endif f2

        if (present(CellLocated)) CellLocated = CellFound

        deallocate(Point)
        nullify(Point)


    end subroutine LocateCellPolygons


   !------------------------------------------------------------------------

    logical function InsideDomainPolygon(DomainPolygon, XPoint, YPoint)

        !Arguments ---------------------------------------------------------
        type(T_Polygon), pointer             :: DomainPolygon
        real   ,                 intent(IN ) :: XPoint, YPoint

        !Local -------------------------------------------------------------
        type(T_PointF ),          pointer    :: Point
        !Begin -------------------------------------------------------------



        allocate(Point)



        Point%X = XPoint
        Point%Y = YPoint

        if (IsPointInsidePolygon(Point, DomainPolygon)) then

            InsideDomainPolygon = .true.

        else

            InsideDomainPolygon = .false.

        endif


        deallocate(Point)
        nullify   (Point)


    end function InsideDomainPolygon


#ifdef _OPENMI_


    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetIUB
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETIUB"::GetIUB
    !DEC$ ENDIF
    integer function GetIUB(HorizontalGridID)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            GetIUB = Me%WorkSize%IUB

        else
            GetIUB = - 99
        end if

        return

    end function GetIUB

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetJUB
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETJUB"::GetJUB
    !DEC$ ENDIF
    integer function GetJUB(HorizontalGridID)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            GetJUB = Me%WorkSize%JUB

        else
            GetJUB = - 99
        end if

        return

    end function GetJUB



    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetCenterXCoordinate
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETCENTERXCOORDINATE"::GetCenterXCoordinate
    !DEC$ ENDIF
    real(8) function GetCenterXCoordinate(HorizontalGridID, i, j)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer                                     :: i, j

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            GetCenterXCoordinate = Me%Compute%XX2D_Z(i, j)
        else
            GetCenterXCoordinate = - 99.0
        end if

        return

    end function GetCenterXCoordinate

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetCenterYCoordinate
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETCENTERYCOORDINATE"::GetCenterYCoordinate
    !DEC$ ENDIF
    real(8) function GetCenterYCoordinate(HorizontalGridID, i, j)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer                                     :: i, j

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            GetCenterYCoordinate = Me%Compute%YY2D_Z(i, j)
        else
            GetCenterYCoordinate = - 99.0
        end if

        return

    end function GetCenterYCoordinate


    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetGridCellCoordinates
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETGRIDCELLCOORDINATES"::GetGridCellCoordinates
    !DEC$ ENDIF
    logical function GetGridCellCoordinates(HorizontalGridID, i, j, xCoords, yCoords)

        !Arguments-------------------------------------------------------------
        integer                                     :: HorizontalGridID
        integer                                     :: i, j
        real(8), dimension(5)                       :: xCoords
        real(8), dimension(5)                       :: yCoords

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_

        call Ready(HorizontalGridID, ready_)

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then


            if (Me%CoordType == SIMPLE_GEOG_ .or. Me%CoordType == GEOG_) then

                !Anticlockwise, closed
                xCoords(1) = Me%LongitudeConn(i, j)
                xCoords(2) = Me%LongitudeConn(i, j+1)
                xCoords(3) = Me%LongitudeConn(i+1, j+1)
                xCoords(4) = Me%LongitudeConn(i+1, j)
                xCoords(5) = Me%LongitudeConn(i, j)

                !Anticlockwise, closed
                yCoords(1) = Me%LatitudeConn(i, j)
                yCoords(2) = Me%LatitudeConn(i, j+1)
                yCoords(3) = Me%LatitudeConn(i+1, j+1)
                yCoords(4) = Me%LatitudeConn(i+1, j)
                yCoords(5) = Me%LatitudeConn(i, j)


            else

                !Anticlockwise, closed
                xCoords(1) = Me%XX_IE(i, j)
                xCoords(2) = Me%XX_IE(i, j+1)
                xCoords(3) = Me%XX_IE(i+1, j+1)
                xCoords(4) = Me%XX_IE(i+1, j)
                xCoords(5) = Me%XX_IE(i, j)

                !Anticlockwise, closed
                yCoords(1) = Me%YY_IE(i, j)
                yCoords(2) = Me%YY_IE(i, j+1)
                yCoords(3) = Me%YY_IE(i+1, j+1)
                yCoords(4) = Me%YY_IE(i+1, j)
                yCoords(5) = Me%YY_IE(i, j)

            endif

            GetGridCellCoordinates = .true.
        else
            call PlaceErrorMessageOnStack("Horizontal Grid not ready")
            GetGridCellCoordinates = .false.
        end if





    end function GetGridCellCoordinates

#endif


end module ModuleHorizontalGrid

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon.
!----------------------------------------------------------------------------------------------------------
