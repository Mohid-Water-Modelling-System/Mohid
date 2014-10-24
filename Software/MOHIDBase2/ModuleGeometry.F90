!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : Geometry
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to calculate the vertical Geometry 
!
!------------------------------------------------------------------------------
! IMPERMEABILITY                : 0/1               -           !Consider impermeable cell faces
! IMPER_COEF_U                  : real             [1]          !
! IMPER_COEFX_U                 : real             [0]          !
! IMPER_COEF_V                  : real             [1]          !
! IMPER_COEFX_V                 : real             [0]          !

!<begindomain>
!   ID                          : int               -           !Domain ID
!   TYPE                        : char              -           !Type of vertical coordinate of the domain
!                                                               !Multiple options: FIXSPACING, SIGMA,
!                                                               !CARTESIAN, FIXSEDIMENT
!   LAYERS                      : int               -           !Number of layers
!   EQUIDISTANT                 : real             [0]          !Equidistant layers spacing in meters
!   LAYERTHICKNESS              : real vector       -           !If not equidistant specifies layers thickness
!                                                               !starting from bottom layer (e.g. 50. 20. 10. 5.)
!   TOLERANCEDEPTH              : real            [0.05]        !Just for SIGMA,ISOPYCNIC coordinates
!   TOTALTHICKNESS              : real              -           !Total domain thickness 
!                                                               !(Just for FIXSPACING, FIXSEDIMENT, SOIL_TOPLAYER)
!   EMPTY_TOP_LAYERS            : int              [0]          !Number of empty layers counting from top
!   DOMAINDEPTH                 : real
!   LAGRANGIAN                  : 0/1              [0]          !Use lagrangian approach for distorting grometry? 
!                                                               !Layers are displaced with vertical velocity
!   MINEVOLVELAYERTHICKNESS     : real            [0.5]         !Allowed distortion in percentage of initial thickness
!                                                               !(if LAGRANGIAN : 1)
!   DISPLACEMENT_LIMIT          : real           [1000]         !Maximum displacement in meters (if LAGRANGIAN : 1)

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

Module ModuleGeometry

    use ModuleGlobalData
    use ModuleTime                 
    use ModuleEnterData          
    use ModuleGridData,         only: GetGridData, UngetGridData, GetMaximumValue,      &
                                      GetGridDataFileName, WriteGridData,               &
                                      GetGridDataEvolution     
    use ModuleHorizontalMap,    only: GetWaterPoints2D, GetExteriorBoundaryFaces,       &
                                      UnGetHorizontalMap
    use ModuleHorizontalGrid,   only: GetHorizontalGridSize, GetHorizontalGrid,         &
                                      GetGridCoordType, GetGridAngle, GetGridZone,      &
                                      GetLatitudeLongitude, GetGridOrigin,              &
                                      GetGridLatitudeLongitude, GetCoordTypeList,       &
                                      GetCheckDistortion, UnGetHorizontalGrid,          &
                                      GetDDecompWorkSize2D
    use ModuleFunctions,        only: SetMatrixValue, SetMatrixValueAllocatable,        &
                                      Chunk_J, Chunk_K, GetPointer
    use ModuleHDF5
    use ModuleStopWatch,        only : StartWatch, StopWatch         

    implicit none
    private

    !Subroutine----------------------------------------------------------------


    !Constructor
    public  :: ConstructGeometry
    private ::      AllocateInstance
    private ::      ConstructGlobalVariables
    private ::          GetDomainsFromFile
    private ::              Add_Domain
    private ::          ComputeLayers
    private ::          AllocateVariables
    private ::          VerifyBathymetry
    private ::      ConstructKFloor

#ifdef _USE_SEQASSIMILATION
    !Subroutine to point to memory space of geometry variables
    public  :: PointToGeometryState
#endif _USE_SEQASSIMILATION

    !Modifier
    public  :: ComputeInitialGeometry       !To call when the initial SurfaceElevation is known
    public  :: ComputeVerticalGeometry      !To call in the "working cycle"
    public  :: UpdateKfloor
    private ::      ComputeSZZ
    private ::          ComputeFixSpacing
!    private ::          ComputeFixSediment
    private ::          ComputeInitSediment
    private ::          ComputeCartesian
!    private ::          ComputeHarmonic
    private ::          ComputeSigma
    private ::          ComputeLagrangianNew
    private ::      ComputeZCellCenter
    private ::      ComputeDistances
    private ::      ComputeAreas
    private ::      ComputeVolumes
    private ::      StoreVolumeZOld

#ifdef _USE_SEQASSIMILATION
    !Copy subroutines usable in sequential data assimilation to change variables' value
    public  :: CopyGeometryDistances
    public  :: CopyGeometryAreas
    public  :: CopyGeometryVolumes
    public  :: CopyGeometryWaterColumn
#endif _USE_SEQASSIMILATION

    !Input / Output
    public  :: ReadGeometryBin
    public  :: WriteGeometryBin
    public  :: ReadGeometryHDF
    public  :: WriteGeometryHDF

    !Selector
    public  :: GetGeometryDistances                 !SZZ, DZZ, DWZ, ZCellCenter, DUZ, DVZ, DWZ_Xgrad, DWZ_Ygrad, 
    public  :: GetGeometryAreas
    public  :: GetGeometryVolumes
    public  :: GetGeometryKFloor
    public  :: GetGeometryKTop
    public  :: GetGeometryWaterColumn               !HT
    public  :: GetGeometryInputFile
    public  :: GetGeometrySize
    public  :: GetGeometryMinWaterColumn
    public  :: GetLayer4Level
    public  :: UnGetGeometry

#ifdef _USE_SEQASSIMILATION
    !Set subroutines usable to point variables to external variables/memory space
    public  :: SetGeometryDistances
    public  :: SetGeometryAreas
    public  :: SetGeometryVolumes
    public  :: SetGeometryWaterColumn

    !Reset subroutine usable to reestablish variables to internal memory space
    public  :: ReSetGeometry
#endif _USE_SEQASSIMILATION

    !Destructor
    public  :: KillGeometry
    private ::      DeallocateInstance
    private ::      DeallocateVariables

#ifdef _USE_SEQASSIMILATION
    !Subroutine to point to memory space of geometry properties
    public  :: NullifyGeometryStatePointer
#endif _USE_SEQASSIMILATION

    !Management
    private ::      Ready
    private ::          LocateObjGeometry

    !Interfaces----------------------------------------------------------------

    private :: UnGetGeometry2Dreal
    private :: UnGetGeometry3Dreal4
    private :: UnGetGeometry3Dreal8
    private :: UnGetGeometry2Dinteger
    interface  UnGetGeometry
        module procedure UnGetGeometry2Dreal
        module procedure UnGetGeometry3Dreal4
        module procedure UnGetGeometry3Dreal8
        module procedure UnGetGeometry2Dinteger
    end interface  UnGetGeometry

    private :: ConstructGeometryV1
    private :: ConstructGeometryV2
    interface  ConstructGeometry
        module procedure ConstructGeometryV1
        module procedure ConstructGeometryV2
    end interface 

    !Parameter-----------------------------------------------------------------

    !STAT
    real   , parameter :: FillValueDouble = -9.9e15

    !Domain types
    integer, parameter :: FixSpacing            = 1
    integer, parameter :: Sigma                 = 2
    integer, parameter :: Isopycnic             = 3
    !integer, parameter :: Lagrangian            = 4
    integer, parameter :: Cartesian             = 5
    !integer, parameter :: Harmonic              = 6
    integer, parameter :: FixSediment           = 7
    integer, parameter :: SigmaTop              = 8
    integer, parameter :: CartesianTop          = 9


    !Other
    integer, parameter :: INITIALGEOMETRY       = 1
    integer, parameter :: TRANSIENTGEOMETRY     = 2

    !Faces Option
    integer, parameter :: MinTickness           = 3
    integer, parameter :: AverageTickness       = 2


    !Type----------------------------------------------------------------------

    type T_Domain
        integer                                 :: ID                        = FillValueInt
        integer                                 :: DomainType                = FillValueInt
        logical                                 :: IsLagrangian              = .false.      !Lagrangian change w/ vert veloc
        integer                                 :: NumberOfLayers            = FillValueInt
        integer                                 :: EmptyTopLayers            = FillValueInt
        integer                                 :: UpperLayer, LowerLayer    = FillValueInt
        integer                                 :: ActiveUpperLayer          = FillValueInt
        real                                    :: DomainDepth               = FillValueReal
        real                                    :: TotalThickness            = FillValueReal 
        real, dimension(:), pointer             :: LayerThickness            => null()
        real, dimension(:), pointer             :: LayerMinThickness         => null()
        real, dimension(:), pointer             :: LayerMaxThickness         => null()
        real                                    :: ToleranceDepth            = FillValueReal 
        real                                    :: MinInitialLayerThickness  = FillValueReal
        real                                    :: MaxThicknessGrad          = FillValueReal 
        real                                    :: MinEvolveLayerThickness   = FillValueReal
        real                                    :: MinEsp                    = FillValueReal
        real                                    :: BottomLayerThickness      = FillValueReal
        real                                    :: GridMovementDump          = FillValueReal
        real                                    :: DisplacementLimit         = FillValueReal 
        integer                                 :: InitializationMethod      = FillValueInt
        real                                    :: Equidistant               = FillValueReal        
        logical                                 :: RomsDistortion            = .false.
        real                                    :: theta_s                   = null_real, & !initialization: Jauch
                                                   theta_b                   = null_real, & !initialization: Jauch
                                                   Hc                        = null_real    !initialization: Jauch
        type (T_Domain), pointer                :: Next                      => null(), &
                                                   Prev                      => null()
    end type T_Domain

    type T_Distances
        real, dimension(:, :, :), allocatable       :: SZZ 
        real, dimension(:, :, :), allocatable       :: DZZ 
        real, dimension(:, :, :), allocatable       :: DWZ, DUZ, DVZ, DZI, DZE, DWZ_Xgrad, DWZ_Ygrad 
        real, dimension(:, :, :), allocatable       :: InitialSZZ 
        real, dimension(:, :, :), allocatable       :: ZCellCenter  !Distance from the refernce level, 
                                                                ! center of cells, positive upwards
    end type T_Distances

    type T_Areas
        real, dimension(:, :, :), allocatable       :: AreaU, AreaV
        logical                                     :: Impermeability = .false.
        real, dimension(      :), allocatable       :: Coef_U, CoefX_U, Coef_V, CoefX_V  
    end type T_Areas

    type T_Volumes
        real(8), dimension(:, :, :), allocatable    :: VolumeZ, VolumeU, VolumeV, VolumeW, VolumeZOld
        logical                                     :: FirstVolW = .true.
    end type T_Volumes

    type T_KFloor
        integer, dimension(:, :), allocatable       :: Z, U, V, Domain
    end type T_KFloor

    type T_KTop
        integer, dimension(:, :), allocatable       :: Z
    end type T_KTop


    type T_WaterColumn
        real, dimension(:, :), pointer          :: U    => null(), &
                                                   V    => null(), &
                                                   Z    => null()                    !Former HT
        real                                    :: Zmin = FillValueReal !Former Hmin
    end type T_WaterColumn

#ifdef _USE_SEQASSIMILATION
    type T_StatePointer
        real,    dimension(:, :, :), pointer    :: SZZ          => null(), &   
                                                   DWZ          => null(), &
                                                   DUZ          => null(), &
                                                   DVZ          => null(), &
                                                   DZZ          => null(), &
                                                   ZCellCenter  => null(), &
                                                   AreaU        => null(), &
                                                   AreaV        => null()
        real,    dimension(:, :),    pointer    :: WaterColumnU => null(), &        
                                                   WaterColumnV => null(), &
                                                   WaterColumnZ => null()
        real(8), dimension(:, :, :), pointer    :: VolumeZ      => null(), &
                                                   VolumeU      => null(), &
                                                   VolumeV      => null(), &
                                                   VolumeZOld   => null()
    end type T_StatePointer
#endif _USE_SEQASSIMILATION

    type T_External
        logical                                 :: ContinuesCompute = .false.
        logical                                 :: NonHydrostatic   = .false.
        real                                    :: BathymTopoFactor = 1.0
        real, dimension(:,:,:), pointer         :: DecayTime        => null()
    end type T_External

    type T_Geometry
        integer                                 :: InstanceID   = null_int !initialization: Jauch
        integer                                 :: FacesOption  = null_int !initialization: Jauch
        type (T_External)                       :: ExternalVar
        type (T_Distances)                      :: Distances
        type (T_Areas)                          :: Areas
        type (T_Volumes)                        :: Volumes
        type (T_WaterColumn)                    :: WaterColumn      
        type (T_Domain), pointer                :: FirstDomain
        type (T_Domain), pointer                :: LastDomain
        type (T_KFloor)                         :: KFloor
        type (T_KTop)                           :: KTop

        type (T_Time)                           :: ActualTime  
        type (T_Size3D)                         :: Size
        type (T_Size3D)                         :: WorkSize
        
        logical                                 :: IsWindow                 = .false. !initialization: Jauch
        
        logical                                 :: LagrangianLimitsComputed = .false.
        
        logical                                 :: BathymNotCorrect         = .false. 
       
        character(len=Pathlength)               :: InputFile                = null_str !initialization: Jauch

#ifdef _USE_SEQASSIMILATION
        !This variable is used to retain location of original memory space for variables
        !changed in sequential data assimilation (some external memory is used ocasionally)
        type(T_StatePointer  )                  :: AuxPointer
#endif _USE_SEQASSIMILATION    

        !Instance of other modules
        integer                                 :: ObjTopography        = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjHorizontalMap     = 0

        type (T_Geometry), pointer              :: Next                 => null()

    end type T_Geometry

    !Global Module Variables
    type (T_Geometry), pointer                  :: FirstGeometry    => null()
    type (T_Geometry), pointer                  :: Me               => null()
   

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructGeometryV1(GeometryID, GridDataID, HorizontalGridID,            &
                                   HorizontalMapID, ActualTime,                         &
                                   NewDomain, SurfaceElevation, BathymTopoFactor,       &
                                   STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        integer                                     :: GridDataID
        integer                                     :: HorizontalGridID
        integer                                     :: HorizontalMapID
        type (T_Time),                 optional     :: ActualTime
        character (len=*), intent(IN), optional     :: NewDomain
        real, dimension(:,:), pointer, optional     :: SurfaceElevation
        real,              intent(IN), optional     :: BathymTopoFactor
        integer,  intent(OUT),         optional     :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mGeometry_)) then
            nullify (FirstGeometry)
            call RegisterModule (mGeometry_) 
        endif

        call Ready(GeometryID, ready_)    

        if (ready_ .EQ. OFF_ERR_) then

            !Allocates Instance
            call AllocateInstance

            !Associates External Instances
            Me%ObjTopography     = AssociateInstance (mGRIDDATA_,       GridDataID      )
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)

            !Actualize the time
            if (present(ActualTime)) then
                Me%ActualTime = ActualTime
            endif

            !Factor For MOHID Water / MOHID Land
            if (present(BathymTopoFactor)) then
                Me%ExternalVar%BathymTopoFactor = BathymTopoFactor
            else
                Me%ExternalVar%BathymTopoFactor = 1.0
            endif
            
            !Construct the variable common to all module
            if (present(NewDomain)) then
                if (present(SurfaceElevation)) then
                    call ConstructGlobalVariables(NewDomain = NewDomain, SurfaceElevation = SurfaceElevation)
                else
                    call ConstructGlobalVariables(NewDomain = NewDomain)
                endif
            else
                if (present(SurfaceElevation)) then
                    call ConstructGlobalVariables (SurfaceElevation = SurfaceElevation)
                else
                    call ConstructGlobalVariables
                endif
            endif

            !Constructs the Matrixes which contains the Ks - KFloor, etc...
            if (present(SurfaceElevation)) then
                call ConstructKFloor (SurfaceElevation)
            else
                call ConstructKFloor
            endif    

            !Returns ID
            GeometryID    = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            
            stop 'Geometry - ConstructGeometryV1 - ERR99' 

        end if 


        if (present(STAT)) STAT = STAT_


    end subroutine ConstructGeometryV1
 
    !--------------------------------------------------------------------------

    subroutine ConstructGeometryV2(GeometryID, GridDataID, HorizontalGridID,            &
                                   HorizontalMapID, KMAX,                               &
                                   STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        integer                                     :: GridDataID
        integer                                     :: HorizontalGridID
        integer                                     :: HorizontalMapID
        integer                                     :: KMAX        
        integer,  intent(OUT),         optional     :: STAT
        
        !Local-----------------------------------------------------------------
        type (T_Size2D)                             :: WorkSize2D, Size2D
        integer                                     :: STAT_, ready_, STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mGeometry_)) then
            nullify (FirstGeometry)
            call RegisterModule (mGeometry_) 
        endif

        call Ready(GeometryID, ready_)    

        if (ready_ .EQ. OFF_ERR_) then

            !Allocates Instance
            call AllocateInstance

            !Associates External Instances
            Me%ObjTopography     = AssociateInstance (mGRIDDATA_,       GridDataID      )
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
 
 
             !Gets horizontal size from the Bathymetry
            call GetHorizontalGridSize(Me%ObjHorizontalGrid,                                 &
                                       Size     = Size2D,                                    &
                                       WorkSize = WorkSize2D,                                &
                                       STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGeometryV2 - Geometry - ERR10'

            Me%Size%ILB = Size2D%ILB
            Me%Size%IUB = Size2D%IUB
            Me%Size%JLB = Size2D%JLB
            Me%Size%JUB = Size2D%JUB

            Me%WorkSize%ILB = WorkSize2D%ILB
            Me%WorkSize%IUB = WorkSize2D%IUB
            Me%WorkSize%JLB = WorkSize2D%JLB
            Me%WorkSize%JUB = WorkSize2D%JUB


            call AllocateVariables(Kmax)

            !Returns ID
            GeometryID    = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            
            stop 'Geometry - ConstructGeometryV2 - ERR99' 

        end if 


        if (present(STAT)) STAT = STAT_


    end subroutine ConstructGeometryV2
 
    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
        type (T_Geometry), pointer                  :: NewObjGeometry
        type (T_Geometry), pointer                  :: PreviousObjGeometry


        !Allocates new instance
        allocate (NewObjGeometry)
        nullify  (NewObjGeometry%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstGeometry)) then
            FirstGeometry            => NewObjGeometry
            Me                       => NewObjGeometry
        else
            PreviousObjGeometry      => FirstGeometry
            Me                       => FirstGeometry%Next
            do while (associated(Me))
                PreviousObjGeometry  => Me
                Me                   => Me%Next
            enddo
            Me                       => NewObjGeometry
            PreviousObjGeometry%Next => NewObjGeometry
        endif

        Me%InstanceID = RegisterNewInstance (mGEOMETRY_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ConstructGlobalVariables(NewDomain, SurfaceElevation)

        !Arguments-------------------------------------------------------------
        character (len=*), intent(IN), optional                 :: NewDomain
        real, dimension(:,:), pointer, optional                 :: SurfaceElevation
        
        !External--------------------------------------------------------------
        integer                                                 :: STAT_CALL
        character(len = StringLength)                           :: Message
        type (T_Size2D)                                         :: Size2D, WorkSize2D

        !----------------------------------------------------------------------

        !Nullify ObjGeometryVariables
        nullify (Me%FirstDomain)

        !Nullify T_Volumes
!        nullify (Me%Volumes%VolumeZ)
!        nullify (Me%Volumes%VolumeU)
!        nullify (Me%Volumes%VolumeV)
!        nullify (Me%Volumes%VolumeW)
!        nullify (Me%Volumes%VolumeZOld)

        !Nullify T_Areas
!        nullify (Me%Areas%AreaU)
!        nullify (Me%Areas%AreaV)

        !Nullify T_Distances
!        nullify (Me%Distances%SZZ        )
!        nullify (Me%Distances%DZZ        )
!        nullify (Me%Distances%DWZ        )
!        nullify (Me%Distances%DUZ        )
!        nullify (Me%Distances%DZI        )
!        nullify (Me%Distances%DZE        )
!        nullify (Me%Distances%DVZ        )
!        nullify (Me%Distances%InitialSZZ )
!        nullify (Me%Distances%ZCellCenter)
!        nullify (Me%Distances%DWZ_Xgrad  )
!        nullify (Me%Distances%DWZ_Ygrad  )


        !Nullify T_KFloor
!        nullify (Me%KFloor%Z)
!        nullify (Me%KFloor%U)
!        nullify (Me%KFloor%V)
!        nullify (Me%KFloor%Domain)
       
        !Gets horizontal size from the Bathymetry
        call GetHorizontalGridSize(Me%ObjHorizontalGrid,                                 &
                                   Size     = Size2D,                                    &
                                   WorkSize = WorkSize2D,                                &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - Geometry - ERR02'

        Me%Size%ILB = Size2D%ILB
        Me%Size%IUB = Size2D%IUB
        Me%Size%JLB = Size2D%JLB
        Me%Size%JUB = Size2D%JUB

        Me%WorkSize%ILB = WorkSize2D%ILB
        Me%WorkSize%IUB = WorkSize2D%IUB
        Me%WorkSize%JLB = WorkSize2D%JLB
        Me%WorkSize%JUB = WorkSize2D%JUB

        if (present(NewDomain)) then
            Me%InputFile = NewDomain
        else
            !Reads the file name of the domain properties
            Message   ='File of the domain properties.'
            call ReadFileName('DOMAIN', Me%InputFile, Message = Message, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - Geometry - ERR03'
        endif 

        !Reads the Domain Properties
        call GetDomainsFromFile

        !Calculates / Verifies Layers / Computes the bounds in K
        call ComputeLayers

        call ConstructImpermeability

        !Checks if Bathymetry is consistent with the tolerance depth - Fromer REBAIXA
        !and if changes bathymetry if cartasian domain type exists
        
        if (present(SurfaceElevation)) then
            call VerifyBathymetry (SurfaceElevation)
        else
            call VerifyBathymetry
        endif
            
        !Allocates variables
        call AllocateVariables

    end subroutine ConstructGlobalVariables

    !--------------------------------------------------------------------------

    subroutine ConstructImpermeability

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ObjEnterData = 0   
        integer                                     :: STATUS
        integer                                     :: iflag
        integer                                     :: WorkKLB, WorkKUB

        !Begin-----------------------------------------------------------------

        call ConstructEnterData(ObjEnterData, Me%InputFile, STAT = STATUS)           
        if (STATUS /= SUCCESS_) stop "ConstructImpermeability - Geometry - ERR01"

            !Searches for the FacesOption
        call GetData(Me%Areas%ImPermeability,                                           &
                     ObjEnterData, iflag,                                               &
                     SearchType     = FromFile,                                         &
                     keyword        = 'FACES_IMPERMEABILITY',                           &
                     Default        = .false.,                                          &  
                     ClientModule   = 'ModuleGeometry',                                 &
                     STAT           = STATUS)

        if (Me%Areas%ImPermeability) then

            allocate (Me%Areas%Coef_U (Me%Size%KLB:Me%Size%KUB))
            allocate (Me%Areas%CoefX_U(Me%Size%KLB:Me%Size%KUB))

            allocate (Me%Areas%Coef_V (Me%Size%KLB:Me%Size%KUB))
            allocate (Me%Areas%CoefX_V(Me%Size%KLB:Me%Size%KUB))

            WorkKLB = Me%WorkSize%KLB
            WorkKUB = Me%WorkSize%KUB

            call GetData(Me%Areas%Coef_U(WorkKLB:WorkKUB),                              &
                         ObjEnterData, iflag,                                           &
                         SearchType     = FromFile,                                     &
                         keyword        = 'IMPER_COEF_U',                               &
                         Default        = 1.,                                           &  
                         ClientModule   ='ModuleGeometry',                              &
                         STAT           = STATUS)

            call GetData(Me%Areas%CoefX_U(WorkKLB:WorkKUB),                             &
                         ObjEnterData, iflag,                                           &
                         SearchType     = FromFile,                                     &
                         keyword        = 'IMPER_COEFX_U',                              &
                         Default        = 0.,                                           &  
                         ClientModule   ='ModuleGeometry',                              &
                         STAT           = STATUS)

            call GetData(Me%Areas%Coef_V(WorkKLB:WorkKUB),                              &
                         ObjEnterData, iflag,                                           &
                         SearchType     = FromFile,                                     &
                         keyword        = 'IMPER_COEF_V',                               &
                         Default        = 1.,                                           &  
                         ClientModule   ='ModuleGeometry',                              &
                         STAT           = STATUS)

            call GetData(Me%Areas%CoefX_V(WorkKLB:WorkKUB),                             &
                         ObjEnterData, iflag,                                           &
                         SearchType     = FromFile,                                     &
                         keyword        = 'IMPER_COEFX_V',                              &
                         Default        = 0.,                                           &  
                         ClientModule   ='ModuleGeometry',                              &
                         STAT           = STATUS)

        endif

        call KillEnterData(ObjEnterData, STAT = STATUS)

    end subroutine ConstructImpermeability

    !--------------------------------------------------------------------------

    subroutine AllocateVariables(Kmax)

        !Parameter-------------------------------------------------------------
        integer, optional               :: Kmax
        !Local-----------------------------------------------------------------
        integer                         :: STATUS
        integer                         :: ILB, IUB, JLB, JUB, KLB, KUB 


        
        if (present(KMAX)) then

            Me%WorkSize%KLB = 1
            Me%WorkSize%KUB = Kmax

            Me%Size%KLB = Me%WorkSize%KLB - 1
            Me%Size%KUB = Me%WorkSize%KUB + 1

        endif


        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB        

        !Allocates T_Volumes
        allocate (Me%Volumes%VolumeZ(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR10'
        Me%Volumes%VolumeZ = FillValueDouble

        allocate (Me%Volumes%VolumeU(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR20'
        Me%Volumes%VolumeU = FillValueDouble

        allocate (Me%Volumes%VolumeV(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR30'
        Me%Volumes%VolumeV = FillValueDouble

        allocate (Me%Volumes%VolumeZOld(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR40'
        Me%Volumes%VolumeZOld = FillValueDouble

        !Allocate T_Areas
        allocate (Me%Areas%AreaU(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR50'
        Me%Areas%AreaU = FillValueReal

        allocate (Me%Areas%AreaV(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR60'
        Me%Areas%AreaV = FillValueReal

        !Allocate T_Distances
        allocate (Me%Distances%SZZ(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR70'
        Me%Distances%SZZ = FillValueReal

        allocate (Me%Distances%DZZ(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR80'
        Me%Distances%DZZ = FillValueReal

        allocate (Me%Distances%DWZ(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR90'
        Me%Distances%DWZ = FillValueReal

        allocate (Me%Distances%DUZ(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR100'
        Me%Distances%DUZ = FillValueReal

        allocate (Me%Distances%DVZ(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR110'
        Me%Distances%DVZ = FillValueReal

        allocate (Me%Distances%InitialSZZ(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR120'
        Me%Distances%SZZ = FillValueReal

        allocate (Me%Distances%ZCellCenter(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR130'
        Me%Distances%ZCellCenter = FillValueReal

        allocate (Me%Distances%DZI(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR140'
        Me%Distances%DZI = FillValueReal

        allocate (Me%Distances%DZE(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR150'
        Me%Distances%DZE = FillValueReal

        allocate (Me%Distances%DWZ_Xgrad(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR160'
        Me%Distances%DWZ_Xgrad = FillValueReal

        allocate (Me%Distances%DWZ_Ygrad(ILB:IUB, JLB:JUB, KLB:KUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR170'
        Me%Distances%DWZ_Ygrad = FillValueReal

        !Allocate T_KFloor
        allocate (Me%KFloor%Z(ILB:IUB, JLB:JUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR180'
        Me%KFloor%Z = FillValueInt

        allocate (Me%KFloor%U(ILB:IUB, JLB:JUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR190'
        Me%KFloor%U = FillValueInt

        allocate (Me%KFloor%V(ILB:IUB, JLB:JUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR200'
        Me%KFloor%V = FillValueInt

        allocate (Me%KFloor%Domain(ILB:IUB, JLB:JUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR210'
        Me%KFloor%Domain = FillValueInt

        allocate (Me%WaterColumn%Z(ILB:IUB, JLB:JUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR220'
        Me%WaterColumn%Z(:,:) = FillValueReal

        allocate (Me%WaterColumn%U(ILB:IUB, JLB:JUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR230'
        Me%WaterColumn%U(:,:) = FillValueReal

        allocate (Me%WaterColumn%V(ILB:IUB, JLB:JUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR240'
        Me%WaterColumn%V(:,:) = FillValueReal
        
        allocate (Me%KTop%Z(ILB:IUB, JLB:JUB), stat = STATUS)
        if (STATUS /= SUCCESS_) stop 'AllocateVariables - Geometry - ERR250'
        Me%KTop%Z(:,:) = FillValueInt

    end subroutine AllocateVariables

    !--------------------------------------------------------------------------

    subroutine GetDomainsFromFile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STATUS
        integer                                     :: FromBlock, FromFile
        logical                                     :: BlockFound
        logical                                     :: BlockLayersFound
        integer                                     :: iflag, NLayer
        integer                                     :: LBo, UBo
        integer                                     :: ClientNumber

        character(len = StringLength), parameter    :: block_begin = '<begindomain>'
        character(len = StringLength), parameter    :: block_end   = '<enddomain>'

        character(len = StringLength), parameter    :: beginlayers = '<<beginlayers>>'
        character(len = StringLength), parameter    :: endlayers   = '<<endlayers>>'
        
        integer                                     :: LagrangianOld_flag

        character(len = StringLength)               :: DomainType
        real, dimension(:), allocatable             :: AuxVector
        type (T_Domain), pointer                    :: NewDomain
        integer                                     :: ObjEnterData = 0
        integer                                     :: i, ActualID, ID, LastLine, FirstLine, Line

        
        LagrangianOld_flag = 0
        
        !Get Enter data parameter
        call GetExtractType(FromBlock = FromBlock, FromFile = FromFile)

        call ConstructEnterData(ObjEnterData, Me%InputFile, STAT = STATUS)           
        if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR10"

        !Searches for the MinWaterColumn
        call GetData(Me%IsWindow, ObjEnterData, iflag,                          &
                     SearchType     = FromFile,                                         &
                     keyword        = 'WINDOW',                                         &
                     Default        = .false.,                                          &  
                     ClientModule   = 'ModuleGeometry',                                 &
                     STAT           = STATUS)
        if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR20"

        !Searches for the MinWaterColumn
        call GetData(Me%WaterColumn%Zmin, ObjEnterData, iflag,                          &
                     SearchType     = FromFile,                                         &
                     keyword        = 'MINIMUMDEPTH',                                   &
                     Default        =  0.1,                                             &  
                     ClientModule   = 'ModuleGeometry',                                 &
                     STAT           = STATUS)
        if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR20"

        !Searches for the FacesOption
        call GetData(Me%FacesOption, ObjEnterData, iflag,                               &
                     SearchType     = FromFile,                                         &
                     keyword        = 'FACES_OPTION',                                   &
                     Default        = AverageTickness,                                  & 
                     ClientModule   ='ModuleGeometry',                                  &
                     STAT = STATUS)
        if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR30"

        if (Me%FacesOption /= AverageTickness .and. Me%FacesOption /= MinTickness) then
            
            write(*,*) "The option FACES_OPTION : ", int(Me%FacesOption), " is not valid"
            write(*,*) "The only valid options are FACES_OPTION : ", AverageTickness, " advised option"
            write(*,*) "The only valid options are FACES_OPTION : ", MinTickness, " advanced user option"
            stop "GetDomainsFromFile - Geometry - ERR40"

        endif

        !The MinTickness option gives bad results espeially in the large scales.
        ! At the estuary scale bad results have also been identified. 
        !It is necessary in the future to understand the reason of this bad results.

        !Rewinds Buffer
        call RewindBuffer(ObjEnterData, STAT = STATUS)
        if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR50"


        !Read the defined domains
        NLayer   = 0
        ActualID = 1
        BlockFound = .true.
        do while (BlockFound)

            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                     &
                                        block_begin, block_end, BlockFound,             &
                                        STAT = STATUS)
            if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR60"
            
Block:      if (BlockFound) then

                !Searches for the ID of the domain
                call GetData(ID, ObjEnterData, iflag,                                   &
                             keyword        = 'ID',                                     &
                             SearchType     = FromBlock,                                &
                             ClientModule   = 'ModuleGeometry',                         &
                             STAT           = STATUS)
                if (STATUS /= SUCCESS_ .or. iflag == 0)                                 &
                    stop "GetDomainsFromFile - Geometry - ERR70"

CorretID:       if (ID == ActualID) then

                    allocate (NewDomain)
                    nullify  (NewDomain%Next)
                    nullify  (NewDomain%Prev)
                    NewDomain%ID = ActualID

                    !Searches for the type of the domain        
                    call GetData(DomainType, ObjEnterData, iflag,                       &
                                 keyword       = 'TYPE',                                &
                                 SearchType    = FromBlock,                             &
                                 ClientModule  = 'ModuleGeometry',                      &
                                 STAT          = STATUS)
                    if (STATUS /= SUCCESS_ .or. iflag == 0)                             &
                        stop "GetDomainsFromFile - Geometry - ERR80"

                    select case (trim(adjustl(DomainType)))

                    case ("FIXSPACING", "Fixspacing", "fixspacing")
                        NewDomain%DomainType = FixSpacing   
                    case ("SIGMA", "Sigma", "sigma")
                        NewDomain%DomainType = Sigma        
                    case ("LAGRANGIAN", "Lagrangian", "lagrangian")
                        !NewDomain%DomainType = Lagrangian
                        !for backward compatibility mark this domain as having lagrangian method
                        !the domain type itself will be taken from intialization method some lines below
                        LagrangianOld_flag     = 1 
                        NewDomain%IsLagrangian = .true.
                        
                        write(*,*)
                        write(*,*)
                        write(*,*) '-------------------------------------------------------------'
                        write(*,*) 'WARNING: Lagrangian domains are deprecated'
                        write(*,*) 'new Keyword LAGRANGIAN : 1 in sigma or cartesian'
                        write(*,*) 'domain blocks is now used. However, for backward compability'
                        write(*,*) 'your domain will be processed with old keywords'
                        write(*,*) '-------------------------------------------------------------'
                        write(*,*)
                        write(*,*)
                           
                    case ("CARTESIAN", "Cartesian", "cartesian")
                        NewDomain%DomainType = Cartesian
                    !case ("HARMONIC", "Harmonic", "harmonic")
                    !    NewDomain%DomainType = Harmonic
                    case ("FIXSEDIMENT", "Fixsediment", "fixsediment")
                        NewDomain%DomainType = FixSediment
                    case ("SIGMATOP", "Sigmatop", "sigmatop")
                        NewDomain%DomainType = SigmaTop
                    case ("CARTESIANTOP", "Cartesiantop", "cartesiantop")
                        NewDomain%DomainType = CartesianTop
                    case default
                        stop "GetDomainsFromFile - Geometry - ERR90"
                    end select

                    !Searches for the number of Layers
                    call GetData(NewDomain%NumberOfLayers, ObjEnterData, iflag,         &
                                 keyword        = 'LAYERS',                             &
                                 SearchType     = FromBlock,                            &
                                 ClientModule   = 'ModuleGeometry',                     &
                                 STAT           = STATUS)
                    if (STATUS /= SUCCESS_ .or. iflag == 0)                             &
                        stop "GetDomainsFromFile - Geometry - ERR100"

                    NLayer = NLayer + NewDomain%NumberOfLayers

                    !Allocates LayerThickness
                    nullify  (NewDomain%LayerThickness)
                    LBo = NLayer - NewDomain%NumberOfLayers + 1                     !This way is usefull
                    UBo = NLayer                                                    !in the ComputeSZZ routines
                    allocate (NewDomain%LayerThickness(LBo:UBo), STAT = STATUS)
                    if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR110"

                    nullify  (NewDomain%LayerMinThickness)
                    nullify  (NewDomain%LayerMaxThickness)
                    
                    allocate (NewDomain%LayerMinThickness(LBo:UBo), STAT = STATUS)
                    if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR111"

                    allocate (NewDomain%LayerMaxThickness(LBo:UBo), STAT = STATUS)
                    if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR112"
                    
                    ! Allows the definition of equidistant layers, layer thickness = constant
                    call GetData(NewDomain%Equidistant,                                 &
                                     ObjEnterData, iflag,                               &
                                     keyword      = 'EQUIDISTANT',                      &
                                     SearchType   = FromBlock,                          &
                                     ClientModule = 'ModuleGeometry',                   &
                                     DEFAULT      = 0. ,                                &
                                     STAT         = STATUS)
                    if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR120"

cd0:                if(iflag == 0) then

                        call ExtractBlockFromBlock(ObjEnterData, ClientNumber,          &
                                                   beginlayers, endlayers,              &
                                                   BlockLayersFound,                    &
                                                   FirstLine = FirstLine,               &
                                                   LastLine  = LastLine,                &
                                                   STAT      = STATUS)
                            
cd1 :                   if (STATUS .EQ. SUCCESS_      ) then    
cd2 :                       if (BlockLayersFound) then                        

                                if ((UBo - LBo + 1)/= (LastLine - FirstLine - 1))       &
                                    stop "GetDomainsFromFile - Geometry - ERR130"

                                Line = FirstLine + 1

                                do  i = UBo, LBo, -1

                                    call GetData(NewDomain%LayerThickness(i), ObjEnterData, &
                                                 iflag, Buffer_Line  = Line, STAT = STATUS)

                                    if (STATUS /= SUCCESS_ .or. iflag /= 1)             &
                                        stop "GetDomainsFromFile - Geometry - ERR140"

                                    Line = Line + 1

                                enddo


                            else

                                allocate (AuxVector(1:UBo-LBo+1), STAT = STATUS)

                                !Searches for the Layer Thickness
                                call GetData(AuxVector, ObjEnterData, iflag,            &
                                            SearchType   = FromBlock,                   &
                                            keyword      = 'LAYERTHICKNESS',            &
                                            ClientModule = 'ModuleGeometry',            &
                                            STAT = STATUS)
                                if (STATUS /= SUCCESS_ .or. iflag == 0)                 &
                                    stop "GetDomainsFromFile - Geometry - ERR150"

                                do i = LBo, UBo
                                    NewDomain%LayerThickness(i) = AuxVector(i-LBo+1)
                                enddo

                                deallocate (AuxVector)

                            endif cd2

                        else if (STATUS .EQ. BLOCK_END_ERR_) then cd1

                            stop "GetDomainsFromFile - Geometry - ERR160"

                        endif cd1

                    else  cd0
                     
                        do i = LBo, UBo
                            NewDomain%LayerThickness(i) = NewDomain%Equidistant
                        end do

                    endif cd0
                    
                    !New keyword
                    !if geometry is lagrangian than SZZ is changed by vertical velocity
                    !and Lagragian approach can be used in any type of coordinates
                    if (LagrangianOld_flag == 0) then    
                        call GetData(NewDomain%IsLagrangian, ObjEnterData, iflag,       &
                                     SearchType   = FromBlock,                          &  
                                     keyword      = 'LAGRANGIAN',                       &
                                     ClientModule = 'ModuleGeometry',                   &
                                     Default      = .false.,                            &
                                     STAT         = STATUS)
                        if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR161"
                   endif
                   
                    !Searches for the Tolerance Depth (TYPE == SIGMA or ISOPYCNIC)
                    if ((NewDomain%DomainType == Sigma)     .or.                        &
                        (NewDomain%DomainType == IsoPycnic)) then
                        
                        call GetData(NewDomain%ToleranceDepth, ObjEnterData, iflag,     &
                                     SearchType   = FromBlock,                          &  
                                     keyword      = 'TOLERANCEDEPTH',                   &
                                     ClientModule = 'ModuleGeometry',                   &
                                     Default      = 0.05,                               &
                                     STAT         = STATUS)
                        if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR170"
                    endif

                    !Searches for Total Thickness / Bottom Surface
                    if (NewDomain%DomainType == FixSpacing  .or.                        &
                        NewDomain%DomainType == SigmaTop) then

                        call GetData(NewDomain%TotalThickness, ObjEnterData, iflag,     &
                                     SearchType     = FromBlock,                        &
                                     keyword        = 'TOTALTHICKNESS',                 &
                                     ClientModule   = 'ModuleGeometry',                 &
                                     STAT           = STATUS)
                        if (STATUS /= SUCCESS_ .or. iflag == 0)                         &
                            stop "GetDomainsFromFile - Geometry - ERR180"

                    elseif(NewDomain%DomainType == FixSediment)then


                        call GetData(NewDomain%TotalThickness, ObjEnterData, iflag,     &
                                     SearchType     = FromBlock,                        &
                                     keyword        = 'TOTALTHICKNESS',                 &
                                     ClientModule   = 'ModuleGeometry',                 &
                                     STAT           = STATUS)
                        if (STATUS /= SUCCESS_)stop "GetDomainsFromFile - Geometry - ERR181"
                        if (iflag == 1)then
                            write(*,*)"TOTALTHICKNESS keyword is now obsolete when using FIXSEDIMENT."
                            write(*,*)"TOTALTHICKNESS keyword will be ignored."
                            write(*,*)"Sediment thickness is now given by a file. Check your options!"
                            stop "GetDomainsFromFile - Geometry - WRN180"
                        endif

                        call GetData(NewDomain%EmptyTopLayers, ObjEnterData, iflag,     &
                                     SearchType     = FromBlock,                        &
                                     keyword        = 'EMPTY_TOP_LAYERS',               &
                                     ClientModule   = 'ModuleGeometry',                 &
                                     STAT           = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetDomainsFromFile - Geometry - ERR190"

                    else

                        if (NewDomain%DomainType /= CartesianTop) then
                            call GetData(NewDomain%DomainDepth, ObjEnterData, iflag,    &
                                         SearchType     = FromBlock,                    &
                                         keyword        = 'DOMAINDEPTH',                &
                                         ClientModule   = 'ModuleGeometry',             &
                                         STAT           = STATUS)
                            if (STATUS /= SUCCESS_ .or. iflag == 0)                     &
                                stop "GetDomainsFromFile - Geometry - ERR200"
                        endif

                    endif

                    !Searches for the coeficient which indicates how much a Lagrangian layer 
                    !can deform 
                    if (NewDomain%IsLagrangian) then
                        
                        !The percentage of initial layer thickness that a layer can deform
                        !it may colapse to a size as (1 - MINEVOLVELAYERTHICKNESS) * InitialThickness
                        !and expand to (1 + MINEVOLVELAYERTHICKNESS)* InitialThickness
                        call GetData(NewDomain%MinEvolveLayerThickness,                 &
                                     ObjEnterData, iflag,                               &
                                     SearchType     = FromBlock,                        &
                                     keyword        = 'MINEVOLVELAYERTHICKNESS',        &
                                     ClientModule   = 'ModuleGeometry',                 &
                                     Default        = 0.5,                              &
                                     STAT           = STATUS)
                        if (STATUS /= SUCCESS_ .or. iflag == 0)                         &
                            stop "GetDomainsFromFile - Geometry - ERR210"
                        
                        if (NewDomain%MinEvolveLayerThickness < 0. .or. NewDomain%MinEvolveLayerThickness >= 1) then
                            write(*,*)
                            write(*,*)
                            write(*,*) 'MINEVOLVELAYERTHICKNESS keyword in Geometry Doamin'
                            write(*,*) 'can not be < 0 or >= 1. This keyword represents'
                            write(*,*) 'a percentage of initial defined thickness and the' 
                            write(*,*) 'limit for geometry variation'
                            write(*,*) 'Please verify values'
                            write(*,*)
                            stop "GetDomainsFromFile - Geometry - ERR212"
                        endif
                        
!                        call GetData(NewDomain%GridMovementDump,                        &
!                                     ObjEnterData, iflag,                               &
!                                     SearchType     = FromBlock,                        &
!                                     keyword        = 'GRIDMOVEMENTDUMP',               &
!                                     ClientModule   = 'ModuleGeometry',                 &
!                                     Default        = 0.0,                              &
!                                     STAT           = STATUS)
!                        if (STATUS /= SUCCESS_)                                         &
!                            stop "GetDomainsFromFile - Geometry - ERR220"

                        !maximum displacement in meters. If too high will not limit more
                        !than MINEVOLVELAYERTHICKNESS factor
                        call GetData(NewDomain%DisplacementLimit,                       &
                                     ObjEnterData, iflag,                               &
                                     keyword        = 'DISPLACEMENT_LIMIT',             &
                                     Default        = 1000.,                            &
                                     SearchType     = FromBlock,                        &
                                     ClientModule   = 'ModuleGeometry',                 &
                                     STAT           = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetDomainsFromFile - Geometry - ERR230"

                       if (LagrangianOld_flag == 1) then
                            call GetData(DomainType,                                        &
                                         ObjEnterData, iflag,                               &
                                         SearchType     = FromBlock,                        &
                                         keyword        = 'INITIALIZATION_METHOD ',         &
                                         ClientModule   = 'ModuleGeometry',                 &
                                         Default        = "SIGMA",                          &
                                         STAT           = STATUS)
                            if (STATUS /= SUCCESS_)                                         &
                                stop "GetDomainsFromFile - Geometry - ERR240"
                            
                            !Lagrangian is no longer a domain but a process.
                            !As so the domaintype here will be used as the domain type
                            !and the previous keywords are used to compute lagrangian method
                            !in any of the two types
                            
                            if     (DomainType .EQ. "SIGMA"     ) then
                                !NewDomain%InitializationMethod = Sigma
                                NewDomain%DomainType = Sigma
                            elseif (DomainType .EQ. "CARTESIAN" ) then
                                !NewDomain%InitializationMethod = Cartesian
                                NewDomain%DomainType = Cartesian
                            else
                                write (*,*) "Initialization Method invalid"
                                stop "GetDomainsFromFile - Geometry - ERR250"
                            endif                       
                        endif
                    endif


                   !Searches for the coeficient which indicates the minimal thickness in % 
                   !of the bottom cells if MinInitialLayerThickness = 1 - classic cartesian coordinates
                   !if < 1 - cartesian with shave cells
                   !In cartesian coordinates the reference layer thickness is compute for the 
                   !maximum depth
!                    if (NewDomain%DomainType == Cartesian     .or.                      &
!                        NewDomain%DomainType == Harmonic) then
                    if (NewDomain%DomainType == Cartesian) then

                        call GetData(NewDomain%MinInitialLayerThickness,                &
                                     ObjEnterData, iflag,                               &
                                     SearchType     = FromBlock,                        &
                                     keyword        = 'MININITIALLAYERTHICKNESS',       &
                                     ClientModule   = 'ModuleGeometry',                 &
                                     Default        = 0.05,                             &
                                     STAT           = STATUS)
                        if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR270"


                        call GetData(NewDomain%MaxThicknessGrad,                        &
                                     ObjEnterData, iflag,                               &
                                     SearchType     = FromBlock,                        &
                                     keyword        = 'MAX_THICKNESS_GRAD',             &
                                     ClientModule   = 'ModuleGeometry',                 &
                                     Default        = 2./NewDomain%MinInitialLayerThickness, &
                                     STAT           = STATUS)
                        if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR275"

                        if (NewDomain%MaxThicknessGrad < 1.) stop "GetDomainsFromFile - Geometry - ERR277"

                    endif

                    if (NewDomain%DomainType == CartesianTop) then
                        NewDomain%MinInitialLayerThickness = 0.0
                    endif


                    !Seraches for the minimum thickness of colapsing cells of the Harmonic domain
!                    if (NewDomain%DomainType == Harmonic) then
!
!                        call GetData(NewDomain%MinEsp,                                  &
!                                     ObjEnterData, iflag,                               &
!                                     SearchType     = FromBlock,                        &
!                                     keyword        = 'MIN_TOP_THICKNESS',              &
!                                     ClientModule   = 'ModuleGeometry',                 &
!                                     Default        = Me%WaterColumn%ZMin,              &
!                                     STAT           = STATUS)
!                        if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR280"
!
!                    endif

                    
                    !Seraches for the minimum thickness of bottom layer
                    if (NewDomain%DomainType == CartesianTop) then

                        call GetData(NewDomain%BottomLayerThickness,                    &
                                     ObjEnterData, iflag,                               &
                                     SearchType     = FromBlock,                        &
                                     keyword        = 'MIN_BOTTOM_THICKNESS',           &
                                     ClientModule   = 'ModuleGeometry',                 &
                                     Default        = Me%WaterColumn%ZMin,              &
                                     STAT           = STATUS)
                        if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR281"

                    endif


                    if (NewDomain%DomainType == Sigma) then

                        call GetData(NewDomain%RomsDistortion,                          &
                                     ObjEnterData, iflag,                               &
                                     SearchType     = FromBlock,                        &
                                     keyword        = 'ROMS_DISTORTION',                &
                                     ClientModule   = 'ModuleGeometry',                 &
                                     Default        = .false.,                          &
                                     STAT           = STATUS)
                        if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR300"
                        
                        if (NewDomain%RomsDistortion) then

                            !theta_s = terrain following coodinates bottom control parameter = 5
                            call GetData(NewDomain%theta_s,                             &
                                         ObjEnterData, iflag,                           &
                                         SearchType     = FromBlock,                    &
                                         keyword        = 'THETA_S',                    &
                                         ClientModule   = 'ModuleGeometry',             &
                                         Default        = 5.,                           &
                                         STAT           = STATUS)
                            if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR310"

                            !theta_b = terrain following coordinates surface control parameter = 0.4
                            call GetData(NewDomain%theta_b,                             &
                                         ObjEnterData, iflag,                           &
                                         SearchType     = FromBlock,                    &
                                         keyword        = 'THETA_B',                    &
                                         ClientModule   = 'ModuleGeometry',             &
                                         Default        = 0.4,                          &
                                         STAT           = STATUS)
                            if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR320"

                            !hc = s coordinate parameter critical depth = 25
                            call GetData(NewDomain%Hc,                                  &
                                         ObjEnterData, iflag,                           &
                                         SearchType     = FromBlock,                    &
                                         keyword        = 'HC',                         &
                                         ClientModule   = 'ModuleGeometry',             &
                                         Default        = 25.,                          &
                                         STAT           = STATUS)
                            if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR330"                        
                        
                        endif
                    endif


                    !Inserts new domain into the domain list
                    !call AddDomainToGeometry(NewDomain)

                    call Add_Domain(NewDomain)

                    !Next ID to search for
                    ActualID = ActualID + 1

                    !Rewinds Buffer
                    call Block_Unlock(ObjEnterData, ClientNumber)
                    call RewindBuffer(ObjEnterData, STAT = STATUS)
                    if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR390"

                endif CorretID

            endif Block


        enddo 

        call Block_Unlock(ObjEnterData, ClientNumber)


        call KillEnterData(ObjEnterData, STAT = STATUS)
        if (STATUS /= SUCCESS_) stop "GetDomainsFromFile - Geometry - ERR400"

    end subroutine GetDomainsFromFile

    !--------------------------------------------------------------------------

    subroutine ComputeLayers

        !Parameter-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Domain), pointer        :: CurrentDomain
        real                            :: RelativSum, Error,           &
                                           BottomDepth, TopDepth
        real(8)                         :: Sum, DomainDif
        integer                         :: iLayer, LayersBelow
        !Begin-----------------------------------------------------------------

        !Computes limits for K
        Me%WorkSize%KUB = 0
        CurrentDomain => Me%FirstDomain
        do while (associated(CurrentDomain))
            Me%WorkSize%KUB = Me%WorkSize%KUB + CurrentDomain%NumberOfLayers
            CurrentDomain => CurrentDomain%Next
        enddo

        Me%WorkSize%KLB = 1
        Me%Size%KLB     = 0
        Me%Size%KUB     = Me%WorkSize%KUB + 1

        LayersBelow = 0        
        CurrentDomain => Me%FirstDomain
        do while (associated(CurrentDomain))

            !Computes UpperLayer and LowerLayer
            CurrentDomain%LowerLayer = LayersBelow + 1
            CurrentDomain%UpperLayer = LayersBelow + CurrentDomain%NumberOfLayers
            if (CurrentDomain%ID == Me%LastDomain%ID) CurrentDomain%ActiveUpperLayer = CurrentDomain%UpperLayer
            LayersBelow = LayersBelow + CurrentDomain%NumberOfLayers

            !Verifies and corrects relativ layer thickness
            !This is just done in the case of Sigma, Fixspacing, FixSediment 
            !initialization
            if (CurrentDomain%DomainType  == Sigma              .or.                     &
                CurrentDomain%DomainType  == Fixspacing         .or.                     &
                CurrentDomain%DomainType  == FixSediment        .or.                     &
                CurrentDomain%DomainType  == SigmaTop           .or.                     &
                CurrentDomain%DomainType  == IsoPycnic) then

                RelativSum = 0.
                do iLayer = CurrentDomain%LowerLayer, CurrentDomain%UpperLayer
                    RelativSum = RelativSum + CurrentDomain%LayerThickness(iLayer)
                enddo
                Error = abs(1.-RelativSum)
                if ( Error > 0.01) then

                    write(*,*)'Layers incorrectly defined - Domain number ',CurrentDomain%ID
                    stop 'ComputeLayers - Geometry - ERR10.'
!PCL
                else if (Error<=0.01 .and. Error > 1.e-6) then
                    do iLayer = CurrentDomain%LowerLayer, CurrentDomain%UpperLayer
                        CurrentDomain%LayerThickness(iLayer) =                               &
                                                      CurrentDomain%LayerThickness(iLayer) - &
                                                      (1. - RelativSum) /                    &
                                                      real(CurrentDomain%NumberOfLayers)
                    enddo
                endif

            else if (CurrentDomain%DomainType  == Cartesian) then 

                Sum = 0.
                do iLayer = CurrentDomain%LowerLayer, CurrentDomain%UpperLayer
                    Sum = Sum + dble(CurrentDomain%LayerThickness(iLayer))
                enddo

                if (CurrentDomain%ID == Me%LastDomain%ID) then
                    TopDepth = 0.
                else
                    TopDepth = CurrentDomain%DomainDepth
                endif

                
                if (CurrentDomain%ID > 1) then
                    BottomDepth = CurrentDomain%Prev%DomainDepth
                else
                    call GetMaximumValue(Me%ObjTopography, BottomDepth)
                endif
                
                DomainDif = BottomDepth - TopDepth   

                !Error = abs(DomainDif - Sum)

                if ((DomainDif - Sum)>1e-7) then
                    if (.not. Me%IsWindow) then
                        write(*,*)'Layers incorrectly defined - Domain number ',CurrentDomain%ID
                        stop 'ComputeLayers - Geometry - ERR20.'
                    endif
                endif

            endif

            !Verifies if there is a fixspacing domain, other then the lower domain
            if (CurrentDomain%ID /= 1 .and. CurrentDomain%DomainType == FixSpacing) then
                write(*,*)'Just the first domain can be of the type FixSpacing'
                write(*,*)'Verify domain number :',CurrentDomain%ID
                stop 'ComputeLayers - Geometry - ERR30.'
            endif

            !Verifies if there is a fixsediment domain, other then the lower domain
            if (CurrentDomain%ID /= 1 .and. CurrentDomain%DomainType == FixSediment) then
                write(*,*)'Just the first domain can be of the type FixSediment and it cant coexist with other domains'
                write(*,*)'Verify domain number :',CurrentDomain%ID
                stop 'ComputeLayers - Geometry - ERR40.'
            endif

            !Verifies if a cartesian domain is above a fixspacing domain
            if (CurrentDomain%ID > 1 .and. CurrentDomain%DomainType == Cartesian .and.   &
                                 (Me%FirstDomain%DomainType == FixSpacing .or.           &
                                  Me%FirstDomain%DomainType == FixSediment)) then
                write(*,*)'A Cartesian domain type cant overlay a Fixspacing domain'
                write(*,*)'Verify domain number :',CurrentDomain%ID
                stop 'ComputeLayers - Geometry - ERR50.'
            endif

!            !Verifies if a harmonic domain is above a fixspacing domain
!            if (CurrentDomain%ID > 1 .and. CurrentDomain%DomainType == Harmonic .and.    &
!                                 (Me%FirstDomain%DomainType == FixSpacing .or.           &
!                                  Me%FirstDomain%DomainType == FixSediment)) then
!                write(*,*)'A Harmonic domain type cant overlay a Fixspacing or FixSediment domain'
!                write(*,*)'Verify domain number :',CurrentDomain%ID
!                stop 'ComputeLayers - Geometry - ERR60.'
!            endif

            !Verifies if a Lagrangian domain is above a fixspacing domain
            if (CurrentDomain%ID > 1 .and. CurrentDomain%IsLagrangian .and.  &
                                 (Me%FirstDomain%DomainType == FixSpacing .or.           &
                                  Me%FirstDomain%DomainType == FixSediment)) then
                write(*,*)'A Largragian domain type cant overlay a Fixspacing or FixSediment domain'
                write(*,*)'Verify domain number :',CurrentDomain%ID
                stop 'ComputeLayers - Geometry - ERR70.'
            endif

            !Verifies if SigmaTop domain is topmost domain
            if (CurrentDomain%DomainType == SigmaTop .and.                               &
                CurrentDomain%ID         /= Me%LastDomain%ID) then
                write(*,*)'A SigmaTop domain type can only be the topmost domain'
                write(*,*)'Verify domain number :',CurrentDomain%ID
                stop 'ComputeLayers - Geometry - ERR10.'
            endif

            !Verifies if there is a sigma domain below a SigmaTop domain (just one layer?)
            if (CurrentDomain%DomainType      == SigmaTop .and. CurrentDomain%ID == 1) then
                write(*,*)'A SigmaTop domain needs a Sigma domain below'
                write(*,*)'Verify domain number :',CurrentDomain%ID
                stop 'ComputeLayers - Geometry - ERR10.'
            endif

            !Verifies if there is a sigma domain below a SigmaTop domain (sigma below)
            if (CurrentDomain%DomainType      == SigmaTop) then
                if (CurrentDomain%Prev%DomainType /= SIGMA) then
                    write(*,*)'A SigmaTop domain needs a Sigma domain below'
                    write(*,*)'Verify domain number :',CurrentDomain%ID
                    stop 'ComputeLayers - Geometry - ERR11.'
                endif
            endif

            !Verifies if CartesianTop domain is topmost domain
            if (CurrentDomain%DomainType == CartesianTop) then
                if (CurrentDomain%ID /= Me%LastDomain%ID .or. CurrentDomain%ID /= 1) then
                    write(*,*)'A CartesianTop domain type can only be the topmost domain'
                    write(*,*)'A CartesianTop domain type can only exist alone'
                    write(*,*)'Verify domain number :',CurrentDomain%ID
                    stop 'ComputeLayers - Geometry - ERR11a.'
                endif
            endif

            !Verifies if Lagrangian method is used with SIGMA or Cartesian
            if (CurrentDomain%IsLagrangian) then
                if (CurrentDomain%DomainType /= Cartesian .and. CurrentDomain%DomainType /= Sigma) then
                    write(*,*)'Lagrangian method can be used with SIGMA'
                    write(*,*)'or CARTESIAN domains. Verify if other is used'
                    stop 'ComputeLayers - Geometry - ERR20.'
                endif
            endif

!            !Verifies if the Minimal Layer thickness of a lagrangian domain is smalles then
!            !50% of the initial layer thickness
!            if (CurrentDomain%DomainType == Lagrangian) then
!                    Error = 1.e10
!                    do iLayer = CurrentDomain%LowerLayer, CurrentDomain%UpperLayer
!                        if (Error > CurrentDomain%LayerThickness(iLayer)) then
!                            Error = CurrentDomain%LayerThickness(iLayer)
!                        endif
!                    enddo
!                    if (CurrentDomain%MinEvolveLayerThickness > 0.50) then
!                        write(*,*)'The MinimalThickness of the layers should not be greater then 50% of'
!                        write(*,*)'of the initial layer thickness.'
!                        write(*,*)'Verify domain number :',CurrentDomain%ID
!                        write(*,*)'Layer                :',iLayer
!                        stop 'ComputeLayers - Geometry - ERR80.'
!                    else if (CurrentDomain%DisplacementLimit > 0.499 * Error) then
!                        write(*,*)'The DisplacementLimit of the layers is higher than half the minimum LayerThickness'
!                        stop 'ComputeLayers - Geometry - ERR90.'
!                   endif
!            endif

            !Verifies if the inital Layer thickness of the cartesian coordinates is not small then 
            !then 1% 
!            if  (CurrentDomain%DomainType           == Cartesian .or.                 &
!                 CurrentDomain%DomainType           == Harmonic ) then
            if  (CurrentDomain%DomainType           == Cartesian) then

                    if (CurrentDomain%MinInitialLayerThickness < 0.01) then
                        write(*,*)'The MinimalThickness of the layers should not be smaller then 1% of'
                        write(*,*)'of the reference layer thickness'
                        write(*,*)'Verify domain number :',CurrentDomain%ID
                        write(*,*)'Layer                :',iLayer
                        stop 'ComputeLayers - Geometry - ERR100.'
                    endif
            endif


        CurrentDomain => CurrentDomain%Next
        enddo


    end subroutine ComputeLayers

    !--------------------------------------------------------------------------
    !In the case that the domain is Sigma, Isopycnic the Bathymetry is tested in
    !order to avoid layers with zero height.
    !In the case of a cartesian domain, the Bathymetry is rearranged, in order to avoid layers to thin
    subroutine VerifyBathymetry (SurfaceElevation)

        !Arguments-------------------------------------------------------------
        real, dimension(:, :), pointer, optional    :: SurfaceElevation

        !Local-----------------------------------------------------------------
        real, dimension(:, :), pointer              :: Bathymetry, NewBathymetry, XX_IE, YY_IE
        real, dimension(:), pointer                 :: XX, YY
        real                                        :: GRID_ANGLE, Latitude, Longitude, XORIG, YORIG
        integer                                     :: ICOORD_TIP, Zone
        integer                                     :: GEOG, UTM, MIL_PORT, SIMPLE_GEOG, NLRD, GRID_COORD

                                                    
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: i, j, iLayer, STAT_CALL
        real                                        :: DomainDepth, Tolerance
        real                                        :: TopDepth, DomainThickness
        real                                        :: LayerTopDepth, LayerMinBottomDepth
        real                                        :: LayerTop, LayerBottom
        real                                        :: DistToBottom, DistToTop
        real                                        :: MinimalThickness, AllmostZero_ 
        real                                        :: BottomDepth
        character(len=StringLength)                 :: BathymetryFile
        character(len=StringLength)                 :: Comment1, Comment2
        logical                                     :: WriteNewBathymetry = .false., Distortion, ConvertsWaterInLand
        type (T_Domain), pointer                    :: CurrentDomain
        type (T_Size2D)                             :: Size2D
        integer                                     :: LengthWithoutExt
        real                                        :: T1, Taux, Tmax, d1, d2
        character(len=StringLength)                 :: FileVersion
        integer                                     :: FileVersionInt

        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize%ILB;  Size2D%ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB;  Size2D%IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB;  Size2D%JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB;  Size2D%JUB = Me%WorkSize%JUB

        call GetGridData(Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR10'
        
        nullify (NewBathymetry)
        allocate(NewBathymetry(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        
        !Copy bathymetry array
        !NewBathymetry(:,:) = Bathymetry(:,:)
        call SetMatrixValue(NewBathymetry, Size2D, Bathymetry)


        CurrentDomain => Me%FirstDomain
do1:    do while (associated(CurrentDomain))
            
cd1:        if ((CurrentDomain%DomainType == Sigma)      .or.                            &
                (CurrentDomain%DomainType == IsoPycnic)) then

                if (CurrentDomain%ID == Me%LastDomain%ID) exit do1 !Dont correct the upper domain
                DomainDepth = CurrentDomain%DomainDepth
                Tolerance   = CurrentDomain%ToleranceDepth
                do j = JLB, JUB
                do i = ILB, IUB
                    if ((Bathymetry(i, j) >  DomainDepth) .and. &
                        (Bathymetry(i, j) <  DomainDepth + Tolerance)) then
                            write(*,*) NewBathymetry(i, j), DomainDepth, Tolerance
                            NewBathymetry(i, j) = DomainDepth
                            WriteNewBathymetry = .true.
                    endif
                enddo
                enddo

            !Correction in the case of cartesian domain
!            else if (CurrentDomain%DomainType           == Cartesian .or.                &
!                     CurrentDomain%DomainType           == Harmonic) then
            else if (CurrentDomain%DomainType           == Cartesian) then

                !Gets the upper and the lower depth of the domain and calculates its
                !thickness
                if (CurrentDomain%ID == Me%LastDomain%ID) then
                    TopDepth = 0.
                else
                    TopDepth = CurrentDomain%DomainDepth
                endif
                if (CurrentDomain%ID > 1) then
                    BottomDepth = CurrentDomain%Prev%DomainDepth
                else
                    call GetMaximumValue(Me%ObjTopography, BottomDepth)
                endif
                DomainThickness = BottomDepth - TopDepth

                !Verifies if the depth of the Bathymetry is between the upper and the lower
                !part of the layer, if so corrects the Bathymetry
                !write(*,*)'Checking if Bathymetry is between LayerTopDepth and LayerMinBottomDepth'

                LayerTopDepth   = TopDepth

doLayer:        do iLayer = CurrentDomain%UpperLayer, CurrentDomain%LowerLayer, -1
                    
                    !MinimalThickness is now defined as a percentage of the initial thickness of each layer  
                    MinimalThickness = CurrentDomain%MinInitialLayerThickness

                    !The LayerThickness is now given in meters, so dont multiply by DomainTickness
                    !Frank - Jan 2001
                    !LayerMinBottomDepth = LayerTopDepth + CurrentDomain%LayerThickness(iLayer) * &
                    !                      DomainThickness * MinimalThickness

                    LayerMinBottomDepth = LayerTopDepth + CurrentDomain%LayerThickness(iLayer) * &
                                          MinimalThickness
                    
doj:                do j = JLB, JUB
doi:                do i = ILB, IUB
                        ! Avoid round-off errors - Hernani Sept 2004
                        ! if ((Bathymetry(i, j) > LayerTopDepth) .and. &
                        !     (Bathymetry(i, j) < LayerMinBottomDepth)) then
                        ! GRiflet, NVerelst - Usar como erro de truncatura na escrita do ficheiro de batimetria
                        ! metade da menor divisao: i.e. 5e-4 m

                        ConvertsWaterInLand = .false.

                        d1 = LayerTopDepth 
                        d2 = LayerTopDepth + CurrentDomain%LayerThickness(iLayer)

                        AllmostZero_ = 5e-4
                        if (Bathymetry(i, j) >= d2) cycle
                        if (Bathymetry(i, j) <= d1) cycle

                        if ((Bathymetry(i, j) - LayerTopDepth       >  AllmostZero_) .and. &
                            (Bathymetry(i, j) - LayerMinBottomDepth < -AllmostZero_)) then
                            ConvertsWaterInLand = .true.
                        endif

                        if (.not. ConvertsWaterInLand) then
                            Tmax =   FillValueReal

                            T1   = Bathymetry(i, j)  - LayerTopDepth 

                            if (T1 <=AllmostZero_) cycle 

                            Taux = Bathymetry(i-1,j) - LayerTopDepth 

                            if (Taux > CurrentDomain%LayerThickness(iLayer)) then
                                Taux = CurrentDomain%LayerThickness(iLayer)
                            endif
                            if (               Taux > Tmax) Tmax = Taux
                            
                            Taux = Bathymetry(i+1,j) - LayerTopDepth 

                            if (Taux > CurrentDomain%LayerThickness(iLayer)) then
                                Taux = CurrentDomain%LayerThickness(iLayer)
                            endif                             
                            if (               Taux > Tmax) Tmax = Taux


                            Taux = Bathymetry(i,j+1) - LayerTopDepth 

                            if (Taux > CurrentDomain%LayerThickness(iLayer)) then
                                Taux = CurrentDomain%LayerThickness(iLayer)
                            endif                          
                               
                            if (               Taux > Tmax) Tmax = Taux

                            Taux = Bathymetry(i,j-1) - LayerTopDepth 

                            if (Taux > CurrentDomain%LayerThickness(iLayer)) then
                                Taux = CurrentDomain%LayerThickness(iLayer)
                            endif                             
                            if (               Taux > Tmax) Tmax = Taux
                             
                            if (Tmax/T1 > CurrentDomain%MaxThicknessGrad) then
                                ConvertsWaterInLand = .true. 
                            endif
                        endif



                        if (ConvertsWaterInLand) then
                            NewBathymetry(i, j) = LayerTopDepth
                            write(*,*)'i              = ', i
                            write(*,*)'j              = ', j
                            write(*,*)'Bathymetry     = ', Bathymetry(i, j)
                            write(*,*)'New Bathymetry = ',NewBathymetry(i, j)
                            WriteNewBathymetry = .true.
                        endif

                    enddo doi
                    enddo doj
                    
                    !
                    !The Layer thickness is now given in meters, so dont multiply by DomainThickness
                    !Frank Jan 2001
                    LayerTopDepth = LayerTopDepth + CurrentDomain%LayerThickness(iLayer)
                enddo doLayer

            elseif (CurrentDomain%DomainType == CartesianTop) then

                do j = JLB, JUB
                do i = ILB, IUB

                    if (Bathymetry(i, j) > -55.) then

                        TopDepth        = SurfaceElevation(i, j)
                        BottomDepth     = Bathymetry(i, j)

                        iLayer          = CurrentDomain%UpperLayer        
                        LayerTop        = TopDepth
                        LayerBottom     = LayerTop - CurrentDomain%LayerThickness(iLayer)
                        do while (iLayer >= CurrentDomain%LowerLayer)
                        
                            AllmostZero_ = AllmostZeroFraction * CurrentDomain%LayerThickness(iLayer)

!                            if (LayerBottom - AllmostZero_ <= BottomDepth  .and. LayerTop + AllmostZero_ >= BottomDepth) then
                            DistToBottom = BottomDepth - LayerBottom
                            DistToTop    = LayerTop - BottomDepth

                            if ((BottomDepth > LayerBottom) .and. (BottomDepth < LayerTop) .and. &
                                (abs(DistToTop) > AllmostZero_) .and. (abs(DistToBottom) > AllmostZero_)) then

                                if ((iLayer .eq. CurrentDomain%UpperLayer) .or. &
                                    (DistToBottom .le. DistToTop)) then

                                    NewBathymetry(i, j)         = LayerBottom
                                    WriteNewBathymetry          = .true.
                                    write(*,*)'i                = ', i
                                    write(*,*)'j                = ', j
                                    write(*,*)'Bathymetry       = ',Bathymetry(i, j)
                                    write(*,*)'New Bathymetry   = ',NewBathymetry(i, j)
                            
                                else 
                                    NewBathymetry(i, j)         = LayerTop
                                    WriteNewBathymetry          = .true.
                                    write(*,*)'i                = ', i
                                    write(*,*)'j                = ', j
                                    write(*,*)'Bathymetry       = ',Bathymetry(i, j)
                                    write(*,*)'New Bathymetry   = ',NewBathymetry(i, j)
                                endif                                
                                
                                exit 

                            else
                            
                                iLayer       = iLayer - 1
                                
                                if (iLayer > 0) then
                                    LayerTop     = LayerBottom
                                    LayerBottom  = LayerTop - CurrentDomain%LayerThickness(iLayer)
                                else
                                    exit
                                endif
                                
                            endif

                        enddo
                    
                    endif

                enddo
                enddo

            endif cd1

            CurrentDomain => CurrentDomain%Next

        enddo do1

        if (WriteNewBathymetry) then
            
            !Gets XX and YY
            call GetHorizontalGrid(Me%ObjHorizontalGrid, XX = XX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR20'

            call GetHorizontalGrid(Me%ObjHorizontalGrid, YY = YY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR30'


            call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,             &
                                   SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD,      &
                                   NLRD = NLRD)

            !Gets the type of Coordinates
            call GetGridCoordType(Me%ObjHorizontalGrid, ICOORD_TIP, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR40'
 
            if    (ICOORD_TIP == SIMPLE_GEOG)then
        
                call GetGridLatitudeLongitude(Me%ObjHorizontalGrid,                         &
                                                GridLatitudeConn  = YY_IE,                  &
                                                GridLongitudeConn = XX_IE,                  &
                                                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR50'

            elseif(ICOORD_TIP == UTM        .or. ICOORD_TIP == MIL_PORT .or. &
                   ICOORD_TIP == GRID_COORD .or. ICOORD_TIP == NLRD)then

                !Gets XX_IE and YY_IE
                call GetHorizontalGrid(Me%ObjHorizontalGrid, XX_IE = XX_IE, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR60'

                call GetHorizontalGrid(Me%ObjHorizontalGrid, YY_IE = YY_IE, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR70'

            else

                write(*,*)'GEOG coordinate type cannot be used in digital terrain generation'
                stop 'VerifyBathymetry - Geometry - ERR80'

            end if



            !Gets the Grid angle
            call GetGridAngle    (Me%ObjHorizontalGrid, GRID_ANGLE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR90'

            !Gets the Latitude / Longitude
            call GetLatitudeLongitude(Me%ObjHorizontalGrid, Latitude, Longitude, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR100'

            !Gets Zone of the Bathymetry
            call GetGridZone         (Me%ObjHorizontalGrid, Zone, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR110'

            !Gets Origin of the Bathymetry
            call GetGridOrigin      (Me%ObjHorizontalGrid, Xorig, Yorig, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR120'

            !Writes new Bathymetry
            call GetGridDataFileName(Me%ObjTopography, FileName = BathymetryFile, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR130'

            BathymetryFile  = adjustl(BathymetryFile)
            Comment1        = "Automatic Generated Grid Data File"
            Comment2        = "Based On "//trim(BathymetryFile)
            LengthWithoutExt= len_trim(BathymetryFile) - 4
            !BathymetryFile  = BathymetryFile(1:LengthWithoutExt)//"_"//".new"
            
            !check if already has version in name (2 carachters and can go up to 99 versions)
            if (BathymetryFile(LengthWithoutExt-3:LengthWithoutExt-2) == "_v") then
                
                !read existing file version
                FileVersion = trim(adjustl(BathymetryFile(LengthWithoutExt-1:LengthWithoutExt)))
                read(FileVersion, *, IOSTAT = STAT_CALL) FileVersionInt
                if (STAT_CALL /= 0) then
                    !user used a letter after "_v". it may happen in a original filename
                    FileVersionInt = 0
                endif
                
                !update file version
                FileVersionInt = FileVersionInt + 1
                write(FileVersion, "(i2)") FileVersionInt
                if (FileVersionInt < 10) then
                    FileVersion = "0"//FileVersion(2:2)
                endif
            else 
                FileVersionInt = 1
            endif
            
            if (FileVersionInt == 1) then
                BathymetryFile  = BathymetryFile(1:LengthWithoutExt)//"_v"//"01.dat"
            else
                BathymetryFile  = BathymetryFile(1:LengthWithoutExt-2)//FileVersion//".dat"
            endif
            
            
            call GetCheckDistortion (Me%ObjHorizontalGrid, Distortion, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR140'

            if (.not. Distortion) then

                call WriteGridData  (FileName       = BathymetryFile,       &
                                     XX             = XX,                   &
                                     YY             = YY,                   &
                                     COMENT1        = Comment1,             &
                                     COMENT2        = Comment2,             &
                                     WorkSize       = Size2D,               & 
                                     CoordType      = ICOORD_TIP,           &
                                     Xorig          = Xorig,                &
                                     Yorig          = Yorig,                &
                                     Zone           = Zone,                 &
                                     GRID_ANGLE     = GRID_ANGLE,           &
                                     Latitude       = Latitude,             &
                                     Longitude      = Longitude,            &
                                     GridData2D_Real= NewBathymetry,        &
                                     FillValue      = -99.,                 &
                                     Overwrite      = ON,                   &
                                     STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR150'

            else

                call WriteGridData  (FileName       = BathymetryFile,       &
                                     ConnectionX    = XX_IE,                &
                                     ConnectionY    = YY_IE,                &
                                     COMENT1        = "***",                &
                                     COMENT2        = "***",                &
                                     WorkSize       = Size2D,               & 
                                     CoordType      = ICOORD_TIP,           &
                                     Xorig          = Xorig,                &
                                     Yorig          = Yorig,                &
                                     Zone           = Zone,                 &
                                     GRID_ANGLE     = GRID_ANGLE,           &
                                     Latitude       = Latitude,             &
                                     Longitude      = Longitude,            &
                                     GridData2D_Real= NewBathymetry,        &
                                     FillValue      = -99.,                 &
                                     Overwrite      = ON,                   &
                                     STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR160'
                
            endif

            deallocate(NewBathymetry)
            nullify   (NewBathymetry)

            !Displays Message to inform
            write(*,*)'A new Bathymetry has been created, which consists with the geometry'
            write(*,*)'Modify the file Nomfich.dat and Re-run the model'
            write(*,*)'New Bathymetry file : ', trim(BathymetryFile)
            stop
        endif

        !Disposes pointer to the Bathymetry
        call UngetGridData(Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBathymetry - Geometry - ERR170'

        !----------------------------------------------------------------------

    end subroutine VerifyBathymetry


    !--------------------------------------------------------------------------

    subroutine ConstructKFloor (SurfaceElevation)

        !Arguments-------------------------------------------------------------
        real, dimension(:, :), pointer, optional    :: SurfaceElevation

        !Local-----------------------------------------------------------------
        real, dimension(:, :), pointer              :: Bathymetry
        integer, dimension(:, :), pointer           :: WaterPoints2D
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: iLayer, i, j, LayersBelow
        real(8)                                     :: BottomDepth, TopDepth, LayerTopDepth
        real(8)                                     :: LayerBottomDepthMin, LayerBottomDepthMax 
        real(8)                                     :: MinimalThickness, DomainThickness
        real(8)                                     :: LayerThicknessMin, LayerThicknessMax
        logical                                     :: FoundKFloor
        type (T_Domain), pointer                    :: CurrentDomain
        integer, dimension(:,:), pointer            :: ExteriorFacesU, ExteriorFacesV
        integer                                     :: STAT_CALL
        logical                                     :: NeedToStop 
        real(8)                                     :: MaxDomainThickness
        real(8)                                     :: TotalLayerThickness, AuxDepth
        real                                        :: Aux4, AllmostZero_
                                                    
        !WorkSize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        !Gets a pointer to the Bathymetry
        call GetGridData(Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructKFloor - Geometry - ERR01'

        !Gets a pointer to 2D WaterPoints
        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructKFloor - Geometry - ERR02'

        !Exterior Boundaries
        call GetExteriorBoundaryFaces(Me%ObjHorizontalMap,                  &
                                      BoundaryPointsFaceU = ExteriorFacesU, &
                                      BoundaryPointsFaceV = ExteriorFacesV, &
                                      STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructKFloor - ModuleGeometry - ERR03'


        !Checks which is the domain close to the bottom
        if (Me%FirstDomain%DomainType == FixSpacing  .or.                   &
            Me%FirstDomain%DomainType == FixSediment .or.                   &
            Me%FirstDomain%DomainType == CartesianTop) then

            do j = JLB, JUB
            do i = ILB, IUB
    
                if (WaterPoints2D(i, j) == WaterPoint) then
                    Me%KFloor%Domain(i, j) = Me%FirstDomain%ID
                endif

            enddo
            enddo

        else
            
            do j = JLB, JUB
            do i = ILB, IUB
    
                if (WaterPoints2D(i, j) == WaterPoint) then

                    Me%KFloor%Domain(i, j) = Me%LastDomain%ID
                    CurrentDomain => Me%FirstDomain
                    do while (associated(CurrentDomain))
                        if (Bathymetry(i, j) >  CurrentDomain%DomainDepth) then
                            Me%KFloor%Domain(i, j) = CurrentDomain%ID
                            exit
                        endif
                        CurrentDomain => CurrentDomain%Next
                    enddo

                endif

            enddo
            enddo

        endif

         
        !Computes for KFloor%Z
        NeedToStop         = .false.
        MaxDomainThickness =  null_real
doj:    do j = JLB, JUB
doi:    do i = ILB, IUB
                
iw:         if (WaterPoints2D(i, j) == WaterPoint) then
                
                LayersBelow = 0
                CurrentDomain => Me%FirstDomain

                do while (CurrentDomain%ID /= Me%KFloor%Domain(i, j))
                    LayersBelow = LayersBelow + CurrentDomain%NumberOfLayers
                    CurrentDomain => CurrentDomain%Next
                enddo

                
                !Set Current domain so that it points to the domain close to the bottom.
                CurrentDomain => Me%FirstDomain
                do while (CurrentDomain%ID /= Me%KFloor%Domain(i, j))
                    CurrentDomain => CurrentDomain%Next
                enddo

                !In the case of a sigma or fixspacing domain
                !KFloorZ is always equal to LayersBelow + 1
                if ((CurrentDomain%DomainType /= Cartesian   )     .and.                 &
!                    (CurrentDomain%DomainType /= Harmonic    )     .and.                 &
                    (CurrentDomain%DomainType /= CartesianTop)) then
                    Me%KFloor%Z(i, j) = 1 + LayersBelow
                else
                    !Domain Top Depth
                    if (CurrentDomain%ID == Me%LastDomain%ID) then
                        if (CurrentDomain%DomainType /= CartesianTop) then
                            TopDepth = 0.
                        else
                            if (present(SurfaceElevation)) then
                                TopDepth = dble(-1.0 * SurfaceElevation(i, j))
                            else
                                TopDepth = 0.
                            endif
                        endif
                    else
                        TopDepth = dble(CurrentDomain%DomainDepth)
                    endif

                    !Domain Bottom Depth
                    if (CurrentDomain%ID > 1) then
                        BottomDepth = dble(CurrentDomain%Prev%DomainDepth)
                    else
                        if (CurrentDomain%DomainType /= CartesianTop) then 
                            call GetMaximumValue(Me%ObjTopography, Aux4)
                            BottomDepth = dble(Aux4)
                        else

                            !To minimize roundoff errors when BathymTopoFactor = 1
                            if (Me%ExternalVar%BathymTopoFactor/=1.0) then
                                BottomDepth = dble(Me%ExternalVar%BathymTopoFactor) * dble(Bathymetry(i, j))
                            else
                                BottomDepth = dble(Bathymetry(i, j))
                            endif
                        endif
                    endif
                    
                    !Domain Thickness
                    DomainThickness = BottomDepth - TopDepth

                    !Verifies if the depth of the Bathymetry is between the upper and the lower
                    !part of the layer
                    FoundKFloor = .false.
                    LayerTopDepth = TopDepth
                    iLayer = CurrentDomain%UpperLayer
                    do while (iLayer >= CurrentDomain%LowerLayer .and. .not. FoundKFloor)

                        !MinimalThickness is now defined as a percentage of the initial thickness of each layer  
                        MinimalThickness = dble(CurrentDomain%MinInitialLayerThickness)

                        
                        !In the case of Cartesian, Harmonic domain
                        !the layer thickness is given in meters
                        if (CurrentDomain%DomainType            == Cartesian    .or.     &
!                            CurrentDomain%DomainType            == Harmonic     .or.     &
                            CurrentDomain%DomainType            == CartesianTop) then

                            LayerThicknessMax   = dble(CurrentDomain%LayerThickness(iLayer))
                            LayerThicknessMin   = dble(CurrentDomain%LayerThickness(iLayer)) * &
                                                  MinimalThickness
                        else
                            LayerThicknessMax   = dble(CurrentDomain%LayerThickness(iLayer)) * DomainThickness
                            LayerThicknessMin   = dble(CurrentDomain%LayerThickness(iLayer)) * &
                                                  DomainThickness * MinimalThickness
                        endif

                        AllmostZero_ = AllmostZeroFraction * CurrentDomain%LayerThickness(iLayer)

                        LayerBottomDepthMin = LayerTopDepth + LayerThicknessMin - AllmostZero_
                        LayerBottomDepthMax = LayerTopDepth + LayerThicknessMax + AllmostZero_
                        
                        AuxDepth = 0.0

                        !To minimize roundoff errors when BathymTopoFactor = 1
                        if (Me%ExternalVar%BathymTopoFactor/=1.0) then

                            AuxDepth = dble(Me%ExternalVar%BathymTopoFactor) * dble(Bathymetry(i, j))
                        
                        else

                            AuxDepth = dble(Bathymetry(i, j))

                        endif

                        if (LayerBottomDepthMin <= AuxDepth .and. LayerBottomDepthMax >= AuxDepth) then

                            Me%KFloor%Z(i, j) = iLayer
                            FoundKFloor = .true.
                            
                        else if (Me%BathymNotCorrect .and. LayerBottomDepthMin > AuxDepth) then

                            Me%KFloor%Z(i, j) = iLayer + 1
                            FoundKFloor = .true.

                        else

                            LayerTopDepth = LayerTopDepth +  LayerThicknessMax

                        endif
                        iLayer = iLayer - 1
                    enddo

                    !Points with Bathymetry > TopDepth
                    if (.not. FoundKFloor) then
                        if (CurrentDomain%DomainType /= CartesianTop) then
                            Me%KFloor%Z(i, j) = CurrentDomain%UpperLayer
                        else
                            MaxDomainThickness  = max(MaxDomainThickness, DomainThickness)
                            TotalLayerThickness = LayerTopDepth - TopDepth
                            NeedToStop = .true.
                        endif
                    endif
                
                endif
            
            else iw

                Me%KFloor%Z(i, j) = FillValueInt
                
            endif iw

        enddo doi
        enddo doj

        if (NeedToStop) then
            write(*,*)'Domain Thickness Exceeds Layer Definition'
            write(*,*)'Maximum DomainThickness : ',MaxDomainThickness
            write(*,*)'Total Layer Thickness   : ',TotalLayerThickness
            stop 'ConstructKFloor - ModuleGeometry - ERR99'
        endif

        !Computes KFloor%U
        do j = JLB, JUB + 1
        do i = ILB, IUB

            !Water in (i,j) and in (i,j-1)? or current face is an exterior face
            if ( ((Me%KFloor%Z(i, j    ) > FillValueInt) .and.           &
                  (Me%KFloor%Z(i, j - 1) > FillValueInt)) .or.           &
                 (ExteriorFacesU(i,j) == 1) ) then
                    Me%KFloor%U(i, j) = max(Me%KFloor%Z(i, j), Me%KFloor%Z(i, j - 1))
            else
                Me%KFloor%U(i, j) = FillValueInt
            endif
        
        enddo
        enddo


        !Computes KFloor%V
        do j = JLB, JUB
        do i = ILB, IUB + 1

            !Water in (i,j) and in (i-1,j)? or current face is an exterior face
            if ( ((Me%KFloor%Z(i,     j) > FillValueInt) .and.           &
                  (Me%KFloor%Z(i - 1, j) > FillValueInt)) .or.           &
                 (ExteriorFacesV(i,j) == 1)) then
                    Me%KFloor%V(i, j) = max(Me%KFloor%Z(i, j), Me%KFloor%Z(i - 1, j))
            else
                Me%KFloor%V(i, j) = FillValueInt
            endif
        
        enddo
        enddo

        call UnGetHorizontalMap(Me%ObjHorizontalMap, ExteriorFacesU, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructKFloor - ModuleGeometry - ERR04'
        
        call UnGetHorizontalMap(Me%ObjHorizontalMap, ExteriorFacesV, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructKFloor - ModuleGeometry - ERR05'


        !Disposes pointer to the Bathymetry
        call UngetGridData(Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructKFloor - ModuleGeometry - ERR06'


        !Disposes pointer to WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructKFloor - ModuleGeometry - ERR07'

        !----------------------------------------------------------------------

    end subroutine ConstructKFloor

    !--------------------------------------------------------------------------

    subroutine Add_Domain(NewDomain)

        !Arguments-------------------------------------------------------------
        type(T_Domain),           pointer     :: NewDomain

        !----------------------------------------------------------------------

        ! Add to the WaterDomain List a new property
        if (.not.associated(Me%FirstDomain)) then
            Me%FirstDomain       => NewDomain
            Me%LastDomain        => NewDomain
        else
            NewDomain%Prev       => Me%LastDomain
            Me%LastDomain%Next   => NewDomain
            Me%LastDomain        => NewDomain
        end if 

        !----------------------------------------------------------------------

    end subroutine Add_Domain 

    !--------------------------------------------------------------------------

#ifdef _USE_SEQASSIMILATION

    subroutine PointToGeometryState(GeometryID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                      :: GeometryID
        integer,            optional, intent(OUT)    :: STAT

        !Local-------------------------------------------------------------------
        integer                                      :: ready_          
        integer                                      :: STAT_    

        !------------------------------------------------------------------------

        !This is a compilation of points (one for each variable) for internal memory spaces

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_) 

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            Me%AuxPointer%SZZ => Me%Distances%SZZ
            
            Me%AuxPointer%DWZ => Me%Distances%DWZ

            Me%AuxPointer%DUZ => Me%Distances%DUZ

            Me%AuxPointer%DVZ => Me%Distances%DVZ

            Me%AuxPointer%DZZ => Me%Distances%DZZ

            Me%AuxPointer%ZCellCenter => Me%Distances%ZCellCenter

            Me%AuxPointer%AreaU => Me%Areas%AreaU

            Me%AuxPointer%AreaV => Me%Areas%AreaV

            Me%AuxPointer%WaterColumnU => Me%WaterColumn%U

            Me%AuxPointer%WaterColumnV => Me%WaterColumn%V

            Me%AuxPointer%WaterColumnZ => Me%WaterColumn%Z

            Me%AuxPointer%VolumeZ => Me%Volumes%VolumeZ

            Me%AuxPointer%VolumeU => Me%Volumes%VolumeU

            Me%AuxPointer%VolumeV => Me%Volumes%VolumeV

            Me%AuxPointer%VolumeZOld => Me%Volumes%VolumeZOld

            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine PointToGeometryState

#endif _USE_SEQASSIMILATION

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ComputeInitialGeometry(GeometryID, WaterPoints3D, SurfaceElevation,       &
                                      ActualTime, ContinuesCompute, NonHydrostatic,      &
                                      SZZ, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        integer, dimension(:, :, :), pointer        :: WaterPoints3D
        real, dimension(:, :), pointer, optional    :: SurfaceElevation
        type(T_Time),        optional               :: ActualTime
        logical, intent(in), optional               :: ContinuesCompute
        logical, intent(in), optional               :: NonHydrostatic
        real, dimension(:, :, :), pointer, optional :: SZZ
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            !Actualize the time
            if (present(ActualTime)) then
                Me%ActualTime = ActualTime
            endif

            if (present(NonHydrostatic)) then
                Me%ExternalVar%NonHydrostatic = NonHydrostatic
            else
                Me%ExternalVar%NonHydrostatic = .false.
            endif

            if (present(ContinuesCompute)) then
                Me%ExternalVar%ContinuesCompute = ContinuesCompute
            else
                Me%ExternalVar%ContinuesCompute = .false.
            endif

cd2 :       if (Me%ExternalVar%ContinuesCompute) then 

                !Computes the Distances
                call ComputeDistances   (WaterPoints3D)

                !Computes the Areas
                call ComputeAreas      

                !Computes the Volumes
                call ComputeVolumes     

                !It is necessary for the Soil model
                call ComputeZCellCenter 

            else cd2

                if (present(SZZ)) then    
                    call SetMatrixValue( GetPointer(Me%Distances%SZZ), Me%Size, SZZ )
                else
                   !Constructs SZZ with the initial surface elevation
                    call ComputeSZZ         (SurfaceElevation, INITIALGEOMETRY, WaterPoints3D = WaterPoints3D)
                endif

                !Computes the Distances
                call ComputeDistances   (WaterPoints3D)

                !Computes the Volumes
                call ComputeVolumes     

                !Computes the Areas
                call ComputeAreas       

                !Stores VolumeZOLD
                !In this case VolumeZOld will be equal to VolumeZ
                call StoreVolumeZOld    

                !It is necessary for the Soil model
                call ComputeZCellCenter 

            end if cd2

            !Computes the WaterColumn
            if (present(SurfaceElevation)) then
                call ComputeWaterColumn (SurfaceElevation)
            endif

            STAT_ = SUCCESS_

        else cd1              

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine ComputeInitialGeometry

    subroutine UpdateKfloor(GeometryID, SurfaceElevation, BathymNotCorrect, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        real, dimension(:, :), pointer, optional    :: SurfaceElevation
        logical,                        optional    :: BathymNotCorrect
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            
            if (present(BathymNotCorrect)) then
                Me%BathymNotCorrect = BathymNotCorrect
            else
                Me%BathymNotCorrect = .false. 
            endif

           !Updates the Matrixes which contains the Ks - KFloor, etc...
            if (present(SurfaceElevation)) then
                call ConstructKFloor (SurfaceElevation)
            else
                call ConstructKFloor
            endif    

            STAT_ = SUCCESS_

        else cd1              

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine UpdateKfloor

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! This subroutine computes the VerticalGeometry
    ! It should be called during the working cycle

    subroutine ComputeVerticalGeometry(GeometryID, WaterPoints3D, SurfaceElevation,      &
                                       ActualTime, VerticalVelocity, DT_Waterlevel,      &
                                       SZZ, DecayTime, KTop, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        integer, dimension(:, :, :), pointer        :: WaterPoints3D
        real, dimension(:, :), pointer, optional    :: SurfaceElevation
        type (T_Time),            optional          :: ActualTime
        real, dimension(:, :, :), pointer, optional :: VerticalVelocity             !Gives the vertical variation
        real, intent(in), optional                  :: DT_Waterlevel                !for the lagragean coordinate
        real, dimension(:, :, :), pointer, optional :: SZZ, DecayTime
        integer, dimension(:, :), pointer, optional :: KTop
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            !Actualize the time
            if (present(ActualTime)) then
                Me%ActualTime = ActualTime
            endif


            if (present(SZZ)) then    
                call SetMatrixValue( GetPointer(Me%Distances%SZZ), Me%Size, SZZ )
            else
               !Computes SZZ
                call ComputeSZZ(SurfaceElevation, TRANSIENTGEOMETRY, VerticalVelocity, DT_Waterlevel, WaterPoints3D)
            endif  

            if(present(KTop))then
                Me%KTop%Z(:,:) = KTop(:,:)
            endif

            if(present(DecayTime)) then
                Me%ExternalVar%DecayTime => DecayTime
            endif

            !Stores VolumeZOld
            call StoreVolumeZOld

            !Computes Distances
            call ComputeDistances(WaterPoints3D)

            !Computes Volumes
            call ComputeVolumes

            !Computes Areas
            call ComputeAreas

            !It is necessary for the Soil model
            call ComputeZCellCenter

            !Computes the WaterColumn
            if (Me%FirstDomain%DomainType /= FixSediment) then
                call ComputeWaterColumn(SurfaceElevation)
            endif

            nullify(Me%Externalvar%DecayTime)


            STAT_ = SUCCESS_

        else cd1              

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine ComputeVerticalGeometry

    !--------------------------------------------------------------------------
    
    subroutine ComputeKTop (SedimentDomain)

        !Arguments-------------------------------------------------------------
        type(T_Domain), pointer             :: SedimentDomain

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer    :: WaterPoints2D
        integer                             :: ILB, IUB, JLB, JUB
        integer                             :: i, j

        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Gets a pointer to 2D WaterPoints
        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeKTop - ModuleGeometry - ERR01'

        do j = JLB, JUB
        do i = ILB, IUB

            if (WaterPoints2D(i, j) == WaterPoint) then
                Me%KTop%Z(i, j) = Me%WorkSize%KUB - SedimentDomain%EmptyTopLayers
            else
                Me%KTop%Z(i, j) = Me%WorkSize%KUB
            endif

        enddo
        enddo

        !Disposes pointer to WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeKTop - ModuleGeometry - ERR02'

    end subroutine ComputeKTop

    !--------------------------------------------------------------------------
    !Computes the WaterColumn (Bathymetry + SurfaceElevation)

    subroutine ComputeWaterColumn(SurfaceElevation)

        !Arguments-------------------------------------------------------------
        real, dimension(:, :), pointer          :: SurfaceElevation

        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), pointer         :: DUZ, DVZ
        real, dimension(:, : ), pointer         :: Bathymetry,   WaterColumnZ,           &
                                                   WaterColumnU, WaterColumnV

        integer, dimension(:, :), pointer       :: WaterPoints2D, KFloorU, KFloorV
        integer                                 :: ILB, IUB, JLB, JUB, KUB
        integer                                 :: i, j, k, STAT_CALL, kbottom
        integer                                 :: CHUNK

        !Worksize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KUB = Me%WorkSize%KUB


        WaterColumnZ => Me%WaterColumn%Z
        WaterColumnU => Me%WaterColumn%U
        WaterColumnV => Me%WaterColumn%V

        KFloorU      => Me%KFloor%U
        KFloorV      => Me%KFloor%V

        DUZ          => Me%Distances%DUZ
        DVZ          => Me%Distances%DVZ



        !Gets Bathymetry
        call GetGridData(Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeWaterColumn - Geometry - ERR01'

        !Gets WaterPoints2D
        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'ComputeWaterColumn - Geometry - ERR02'

        if (MonitorPerformance) call StartWatch ("ModuleGeometry", "ComputeWaterColumn")

        CHUNK = Chunk_J(JLB,JUB)

        !Computes WaterColumn
        !$OMP PARALLEL PRIVATE(i,j,k,kbottom)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (WaterPoints2D(i, j) ==  WaterPoint)                                      &
                WaterColumnZ(i, j) = Bathymetry(i, j) + SurfaceElevation(i, j)
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        

        !Computes WaterColumnU
        !$OMP PARALLEL PRIVATE(i,j,k,kbottom)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB+1, JUB
        do i = ILB  , IUB

            if (WaterPoints2D(i, j - 1) == WaterPoint .and.                              &
                WaterPoints2D(i, j    ) == WaterPoint) then

                WaterColumnU (i, j) = 0.

                kbottom             = KFloorU(i, j)

                do k = kbottom, KUB
                
                    WaterColumnU (i, j) =  WaterColumnU (i, j) + DUZ(i, j, k)

                enddo

            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        !Computes WaterColumnV
        !$OMP PARALLEL PRIVATE(i,j,k,kbottom)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB  , JUB
        do i = ILB+1, IUB

            if (WaterPoints2D(i - 1, j) == WaterPoint .and.                              &
                WaterPoints2D(i    , j) == WaterPoint) then

                WaterColumnV (i, j) = 0.

                kbottom             = KFloorV(i, j)

                do k = kbottom, KUB
                
                    WaterColumnV (i, j) =  WaterColumnV (i, j) + DVZ(i, j, k)

                enddo

            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleGeometry", "ComputeWaterColumn")

        nullify(WaterColumnZ, WaterColumnU, WaterColumnV)
        nullify(KFloorU, KFloorV)
        nullify(DUZ, DVZ)

        !UnGets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'ComputeWaterColumn - Geometry - ERR03'


        !UnGets Bathymetry
        call UngetGridData(Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'ComputeWaterColumn - Geometry - ERR04'

        !----------------------------------------------------------------------

     end subroutine ComputeWaterColumn

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    !Computes the Distances (DWZ + DZZ + DUZ + DVZ + DZI + DZE)

    subroutine ComputeDistances(WaterPoints3D)

        !Arguments-------------------------------------------------------------
        integer, dimension(:, :, :), pointer    :: WaterPoints3D

        !Local-----------------------------------------------------------------
        real   , dimension(:, :, :), pointer    :: DWZ, DZZ, DUZ, DVZ, SZZ, DZI, DZE, DWZ_Xgrad, DWZ_Ygrad 
        real   , dimension(:, :   ), pointer    :: DUX, DVY
        real                                    :: aux
        integer                                 :: FacesOption
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: i, j, k, STAT_CALL
        integer                                 :: CHUNK
        
        !Begin-----------------------------------------------------------------

        !Worksize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        DZZ => Me%Distances%DZZ
        DWZ => Me%Distances%DWZ
        DUZ => Me%Distances%DUZ
        DVZ => Me%Distances%DVZ
        SZZ => Me%Distances%SZZ
        DZI => Me%Distances%DZI
        DZE => Me%Distances%DZE

        DWZ_Xgrad  => Me%Distances%DWZ_Xgrad
        DWZ_Ygrad  => Me%Distances%DWZ_Ygrad

        FacesOption = Me%FacesOption

        if (MonitorPerformance) call StartWatch ("ModuleGeometry", "ComputeDistances")

        CHUNK = Chunk_J(JLB,JUB)

        !Computes DWZ
        !$OMP PARALLEL PRIVATE(i,j,k)
        do k = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB

            if (WaterPoints3D(i, j, k) == WaterPoint) then

                DWZ(i, j, k) = SZZ(i, j, k-1) - SZZ(i, j, k  )

            endif
        enddo
        enddo
        !$OMP END DO NOWAIT
        enddo
        !$OMP END PARALLEL

cd1:    if (FacesOption == MinTickness) then

            !Computes DUZ
            !$OMP PARALLEL PRIVATE(i,j,k)
            do k = KLB  , KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = JLB+1, JUB
            do i = ILB  , IUB

                if (WaterPoints3D(i, j - 1, k) == WaterPoint .and.                            &
                    WaterPoints3D(i, j    , k) == WaterPoint) then
                    DUZ(i, j, k) =  min(DWZ(i, j-1, k), DWZ(i, j, k)) 
                endif

            enddo
            enddo
            !$OMP END DO NOWAIT
            enddo
            !$OMP END PARALLEL

            !Computes DVZ
            !$OMP PARALLEL PRIVATE(i,j,k)
            do k = KLB  , KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = JLB  , JUB
            do i = ILB+1, IUB

                if (WaterPoints3D(i - 1, j, k) == WaterPoint .and.                            &
                    WaterPoints3D(i    , j, k) == WaterPoint) then

                    DVZ(i, j, k) =  min(DWZ(i-1, j, k), DWZ(i, j, k)) 

                endif

            enddo
            enddo
            !$OMP END DO NOWAIT
            enddo
            !$OMP END PARALLEL

        else if (FacesOption == AverageTickness) then cd1

            !! $OMP MASTER
            !Gets DZX, DZY
            call GetHorizontalGrid(Me%ObjHorizontalGrid, DUX = DUX, DVY = DVY, & 
                                   STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ComputeAreas - Geometry - ERR01'
            !! $OMP END MASTER

            !! $OMP BARRIER

            !Computes DUZ
            !$OMP PARALLEL PRIVATE(i,j,k)
            do k = KLB  , KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = JLB+1, JUB
            do i = ILB  , IUB

                if (WaterPoints3D(i, j - 1, k) == WaterPoint .and.                            &
                    WaterPoints3D(i, j    , k) == WaterPoint) then

                    DUZ(i, j, k) =  (DWZ(i, j, k)     * DUX(i, j - 1) +                  &
                                     DWZ(i, j - 1, k) * DUX(i, j))    /                  &
                                    (DUX(i, j - 1)    + DUX(i, j)) 
                
                endif

            enddo
            enddo
            !$OMP END DO NOWAIT
            enddo
            !$OMP END PARALLEL


            !Computes DVZ
            !$OMP PARALLEL PRIVATE(i,j,k)
            do k = KLB  , KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = JLB  , JUB
            do i = ILB+1, IUB

                if (WaterPoints3D(i - 1, j, k) == WaterPoint .and.                            &
                    WaterPoints3D(i    , j, k) == WaterPoint) then

                    DVZ(i, j, k) =  (DWZ(i, j, k)     * DVY(i - 1, j) +                  &
                                     DWZ(i - 1, j, k) * DVY(i, j))    /                  &
                                    (DVY(i - 1, j)    + DVY(i, j))

                endif

            enddo
            enddo
            !$OMP END DO NOWAIT
            enddo
            !$OMP END PARALLEL

            !! $OMP MASTER
            !Nullifies auxilary pointers
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, DUX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ComputeDistances - Geometry - ERR02'


            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, DVY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ComputeDistances - Geometry - ERR03'
            !! $OMP END MASTER

            !! $OMP BARRIER

        endif cd1

        !Computes DZZ
        !$OMP PARALLEL PRIVATE(i,j,k)
        do k = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB

            if (WaterPoints3D(i, j, k) == WaterPoint) then
                DZZ(i, j, k) =  (DWZ(i, j, k+1) + DWZ(i, j, k)) / 2.0 
            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        enddo
        !$OMP END PARALLEL

        !Computes DZE
        !$OMP PARALLEL PRIVATE(i,j,k)
        do k = KLB  , KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB+1, JUB
        do i = ILB  , IUB

            if (WaterPoints3D(i, j - 1, k) == WaterPoint .and.                            &
                WaterPoints3D(i, j    , k) == WaterPoint) then

                DZE(i, j, k) =  (DUZ(i, j, k+1) + DUZ(i, j, k)) / 2.
            
            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        enddo
        !$OMP END PARALLEL

        !Computes DZI
        !$OMP PARALLEL PRIVATE(i,j,k)
        do k = KLB  , KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB  , JUB
        do i = ILB+1, IUB

            if (WaterPoints3D(i - 1, j, k) == WaterPoint .and.                            &
                WaterPoints3D(i    , j, k) == WaterPoint) then

                DZI(i, j, k) =  (DVZ(i, j, k+1) + DVZ(i, j, k)) / 2.
            
            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        enddo
        !$OMP END PARALLEL

!Computes DWZ_Xgrad
        !$OMP PARALLEL PRIVATE(i,j,k,aux)
        do k = KLB  , KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB  , JUB
        do i = ILB  , IUB

            if (WaterPoints3D(i, j - 1, k) == WaterPoint .and.                            &
                WaterPoints3D(i, j    , k) == WaterPoint) then
                
                aux=(DWZ(i, j-1, k) + DWZ(i, j, k))

                if (DUZ(i, j, k) > 0. .and. aux > 0.) then
                    DWZ_Xgrad(i, j, k) = DWZ(i, j, k) / aux
                else
                    DWZ_Xgrad(i, j, k) = 0.5
                endif

            else if (WaterPoints3D(i, j - 1, k) == WaterPoint) then
 
                DWZ_Xgrad(i, j, k) =  0.
 
            else if (WaterPoints3D(i, j    , k) == WaterPoint) then
            
                DWZ_Xgrad(i, j, k) =   1.

            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        enddo
        !$OMP END PARALLEL

!Computes DWZ_Ygrad
        !$OMP PARALLEL PRIVATE(i,j,k,aux)
        do k = KLB  , KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB  , JUB
        do i = ILB  , IUB

            if (WaterPoints3D(i - 1, j, k) == WaterPoint .and.                            &
                WaterPoints3D(i    , j, k) == WaterPoint) then
                
                aux= (DWZ(i, j, k) + DWZ(i-1, j, k))
                
                if (DVZ(i, j, k) > 0. .and. aux > 0.) then
                    DWZ_Ygrad(i, j, k) =  DWZ(i, j, k)/aux
                else
                    DWZ_Ygrad(i, j, k) = 0.5
                endif
                
            else if (WaterPoints3D(i - 1, j, k) == WaterPoint) then
            
                DWZ_Ygrad(i, j, k) =  0.
            
            else if (WaterPoints3D(i    , j, k) == WaterPoint) then

                DWZ_Ygrad(i, j, k) =   1.
            
            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        enddo
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleGeometry", "ComputeDistances")

        nullify(DWZ, DZZ, DUZ, DVZ, DZI, DZE, SZZ, DWZ_Xgrad, DWZ_Ygrad)

    end subroutine ComputeDistances

    !--------------------------------------------------------------------------
    !Computes the Volumes (VolumeZ, VolumeU, VolumeV)

    subroutine ComputeVolumes

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:, :), pointer          :: DUX, DVY, DZX, DZY
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: i, j, k, STAT_CALL
        integer                                 :: CHUNK

        !Worksize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        !Gets DUX, DVY
        call GetHorizontalGrid(Me%ObjHorizontalGrid, DUX = DUX, DVY = DVY,      &
                               DZX = DZX, DZY = DZY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            stop 'ComputeVolumes - Geometry - ERR01'

        if (MonitorPerformance) call StartWatch ("ModuleGeometry", "ComputeVolumes")

        CHUNK = Chunk_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(i,j,k)

        !Computes VolumeZ
        do k = KLB, KUB
        !$OMP DO SCHEDULE(STATIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            !Calculo dos volumes das celulas dos Z's como a area projectada * a altura media
            ! de cada prisma 
            Me%Volumes%VolumeZ(i, j ,k) = dble(Me%Distances%DWZ(i, j, k)) * &
                                          dble(DUX(i, j)) * dble(DVY(i, j))
        enddo
        enddo
        !$OMP END DO
        enddo

        !Computes VolumeU
        do k = KLB  , KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB+1, JUB
        do i = ILB  , IUB
            !Calculo dos volumes das celulas dos U's como a area projectada * a altura media
            Me%Volumes%VolumeU(i, j, k) = (Me%Volumes%VolumeZ(i, j-1, k)  + &
                                           Me%Volumes%VolumeZ(i, j, k  )) / &
                                           2.0 
        enddo
        enddo
        !$OMP END DO
        enddo

        !Computes VolumeV
        do k = KLB  , KUB
        !$OMP DO SCHEDULE(STATIC, CHUNK)
        do j = JLB  , JUB
        do i = ILB+1, IUB
            !Calculo dos volumes das celulas dos V's 
            Me%Volumes%VolumeV(i, j, k) = (Me%Volumes%VolumeZ(i-1, j, k)  + &
                                           Me%Volumes%VolumeZ(i, j, k  )) / &
                                           2.0
        enddo
        enddo
        !$OMP END DO NOWAIT
        enddo

        if (Me%ExternalVar%NonHydrostatic) then

            !$OMP MASTER
            if (Me%Volumes%FirstVolW) then

                allocate (Me%Volumes%VolumeW(Me%Size%ILB:Me%Size%IUB, &
                                             Me%Size%JLB:Me%Size%JUB, &
                                             Me%Size%KLB:Me%Size%KUB),&
                                             stat = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeVolumes - Geometry - ERR42'
                Me%Volumes%VolumeW = FillValueDouble

                Me%Volumes%FirstVolW = .false.
            endif
            !$OMP END MASTER
            
            !$OMP BARRIER

            !Computes VolumeW
            do k = KLB + 1 , KUB
            !$OMP DO SCHEDULE(STATIC, CHUNK)
            do j = JLB     , JUB
            do i = ILB     , IUB
                !Calculo dos volumes das celulas dos W's 
                Me%Volumes%VolumeW(i, j, k) = (Me%Volumes%VolumeZ(i, j, k-1)  + &
                                               Me%Volumes%VolumeZ(i, j, k  )) / &
                                               2.0
            enddo
            enddo
            !$OMP END DO
            enddo

        endif

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleGeometry", "ComputeVolumes")

        !Nullifies auxilary pointers
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, DUX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            stop 'ComputeVolumes - Geometry - ERR02'


        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, DVY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            stop 'ComputeVolumes - Geometry - ERR03'


        !Nullifies auxilary pointers
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, DZX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            stop 'ComputeVolumes - Geometry - ERR04'


        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, DZY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            stop 'ComputeVolumes - Geometry - ERR05'

        !----------------------------------------------------------------------

    end subroutine ComputeVolumes

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !Computes the Areas(AreaU, AreaV)

    subroutine ComputeAreas

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:, :, :), pointer       :: DUZ, DVZ
        real, dimension(:, :),    pointer       :: DXX, DYY
        real                                    :: WaterLevel
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: i, j, k, STAT_CALL
        integer                                 :: CHUNK

        !Worksize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        DUZ => Me%Distances%DUZ
        DVZ => Me%Distances%DVZ

        !Gets DZX, DZY
        call GetHorizontalGrid(Me%ObjHorizontalGrid, DXX = DXX, DYY = DYY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeAreas - Geometry - ERR01'

        if (MonitorPerformance) call StartWatch ("ModuleGeometry", "ComputeAreas")

        CHUNK = Chunk_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(i,j,k,WaterLevel)            

        do k = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            !Calculo das areas nas faces nos pontos de calculo U's (AREA_U(I,J,K))
!PCL
            Me%Areas%AreaU(i, j, k) = DUZ(i, j, k) * DYY (i, j)  

!PCL
            Me%Areas%AreaV(i, j, k) = DVZ(i, j, k) * DXX (i, j)

            if (Me%Areas%Impermeability) then

                WaterLevel = max(-Me%Distances%SZZ(i, j, KUB), -Me%Distances%SZZ(i, j+1, KUB))

                Me%Areas%AreaU(i, j, k) = Me%Areas%AreaU(i, j, k) *             &
                                         (Me%Areas%Coef_U(k) + WaterLevel * Me%Areas%CoefX_U(k)) 
                                                
                WaterLevel = max(-Me%Distances%SZZ(i, j, KUB), -Me%Distances%SZZ(i+1, j, KUB))

                Me%Areas%AreaV(i, j, k) = Me%Areas%AreaV(i, j, k) *             &
                                         (Me%Areas%Coef_V(k) + WaterLevel * Me%Areas%CoefX_V(k))
                                                   


            endif

        enddo
        enddo
        !$OMP END DO
        enddo

        !$OMP END PARALLEL          

        if (MonitorPerformance) call StopWatch ("ModuleGeometry", "ComputeAreas")

        !Nullifies auxilary pointers
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, DXX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeAreas - Geometry - ERR02'
            


        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, DYY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeAreas - Geometry - ERR03'
            
        
        !nullify local pointer
        nullify(DUZ, DVZ)

        !----------------------------------------------------------------------
    
    end subroutine ComputeAreas

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    !Stores the VolumeZOld (VolumeZOLD = VolumeZ)

    subroutine StoreVolumeZOld

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: i, j, k
        integer                                 :: CHUNK

        !Worksize
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        if (MonitorPerformance) call StartWatch ("ModuleGeometry", "StoreVolumeZOld")

        CHUNK = Chunk_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(i,j,k) 

        do k = KLB, KUB
        !$OMP DO SCHEDULE(STATIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            Me%Volumes%VolumeZOld(i, j ,k) = Me%Volumes%VolumeZ(i, j ,k) 
        enddo
        enddo
        !$OMP END DO
        enddo

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleGeometry", "StoreVolumeZOld")

    end subroutine StoreVolumeZOld


    !--------------------------------------------------------------------------
    !For every domain calls the respective computation rotine
    subroutine ComputeSZZ (SurfaceElevation, ComputionType, VerticalVelocity, DT_Waterlevel, WaterPoints3D)

        !Parameter-------------------------------------------------------------
        real, dimension(:, :), pointer                 :: SurfaceElevation
        integer                                        :: ComputionType
        real, dimension(:, :, :), optional, pointer    :: VerticalVelocity
        real, intent(in), optional                     :: DT_Waterlevel
        integer, dimension(:, :, :), optional, pointer :: WaterPoints3D

        !Esternal--------------------------------------------------------------

        integer :: STAT_CALL

        !Local-----------------------------------------------------------------
        real,    dimension(:, :), pointer           :: Bathymetry
        integer, dimension(:, :), pointer           :: WaterPoints2D
        integer                                     :: i, j, STAT, KLandFace
        logical                                     :: BathymEvolution
        type (T_Domain), pointer                    :: CurrentDomain

        !Gets Bathymetry
        call GetGridData(Me%ObjTopography, Bathymetry, STAT)
        if (STAT /= SUCCESS_) stop 'ComputeSZZ - Geometry - ERR01'

        call GetGridDataEvolution(Me%ObjTopography, BathymEvolution, STAT)
        if (STAT /= SUCCESS_) stop 'ComputeSZZ - Geometry - ERR01'

        !Gets a pointer to 2D WaterPoints
        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeSZZ - Geometry - ERR02'



        if (ComputionType == INITIALGEOMETRY .or. BathymEvolution) then

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (WaterPoints2D(i, j) == WaterPoint) then

                    KLandFace = Me%KFloor%Z(i, j) - 1

                    !To minimize roundoff errors when BathymTopoFactor = 1
                    if (Me%ExternalVar%BathymTopoFactor/=1.0) then

                        Me%Distances%SZZ(i, j, KLandFace) = Me%ExternalVar%BathymTopoFactor * Bathymetry(i, j)

                    else

                        Me%Distances%SZZ(i, j, KLandFace) = Bathymetry(i, j)

                    endif

                endif

            enddo
            enddo
        endif

        !When this rotine is called from the constructor, all domains are initialized, so the
        !variable ComputionType should be equal to INITIALGEOMETRY. In the case that this 
        !rotine is called from the modifier, the variable ComputionType should be equal to 
        !TRANSIENTGEOMETRY 

        CurrentDomain => Me%FirstDomain
        do while (associated(CurrentDomain))
                  
            select case (CurrentDomain%DomainType)

                case (FixSpacing)
                    if (ComputionType == INITIALGEOMETRY)                                &
                        call ComputeFixSpacing  (CurrentDomain)

                case (FixSediment)
                    if (ComputionType == INITIALGEOMETRY) then
                        call ComputeKTop        (CurrentDomain)
                        call ComputeInitSediment (CurrentDomain)
                    else
                        !call ComputeFixSediment (CurrentDomain)
                    endif

                case (SigmaTop)

                    if (ComputionType == INITIALGEOMETRY) then
                        call ComputeSigma (SurfaceElevation, CurrentDomain)
                    else
                        !Do Nothing
                    endif
                    
                case (Sigma)
                    
                    if ((ComputionType == INITIALGEOMETRY) .or.                         &
                        ((CurrentDomain%ID == Me%LastDomain%ID) .and.                   &
                         (.not. CurrentDomain%IsLagrangian)))  then        
                        
                        call ComputeSigma(SurfaceElevation, CurrentDomain)  
                    
                    else 
                        call ComputeLagrangianNew(SurfaceElevation, VerticalVelocity,   &
                                                  DT_Waterlevel, CurrentDomain)
                    endif
                    
                case (Isopycnic)
!                    call ConstructIsopycnic(ObjGeometry, CurrentDomain)
                
                !Lagrangian moved to cartesian and sigma. 
                !Lagrangian is not a geometry but a computation method
!                case (Lagrangian)
!                    if (ComputionType == INITIALGEOMETRY) then
!                        if (CurrentDomain%InitializationMethod == Sigma) then
!                            call ComputeSigma(SurfaceElevation, CurrentDomain)
!                        else
!                            call ComputeCartesian(SurfaceElevation, CurrentDomain, ComputionType)
!                        endif
!                    else
!                        call ComputeLagrangianNew(SurfaceElevation, VerticalVelocity,       &
!                                                  DT_Waterlevel, CurrentDomain)
!                    endif

                case (Cartesian)
                    
                    if (.not. CurrentDomain%IsLagrangian) then
                        if ((ComputionType == INITIALGEOMETRY) .or.                            &
                            (CurrentDomain%ID == Me%LastDomain%ID)) then          
                            
                            call ComputeCartesian(SurfaceElevation, CurrentDomain, ComputionType)
                        
                        endif
                    else
                        if (ComputionType == INITIALGEOMETRY) then
                            
                            if (CurrentDomain%ID == Me%LastDomain%ID) then 
                                !cartesian domain has to adapt to surface elevation
                                call ComputeCartesianNew(SurfaceElevation, CurrentDomain)    
                            else
                                call ComputeCartesian(SurfaceElevation, CurrentDomain, ComputionType)
                            endif
                            
                        else        
                            call ComputeLagrangianNew(SurfaceElevation, VerticalVelocity,    &
                                                      DT_Waterlevel, CurrentDomain)
                        endif
                    endif
                    
!                !This geometry will be discontinued
!                case (Harmonic)
!
!                    call ComputeHarmonic(SurfaceElevation, CurrentDomain)

                case (CartesianTop)
                    if (ComputionType == INITIALGEOMETRY) then
                        if (present(WaterPoints3D)) then
                            call ComputeCartesian(SurfaceElevation, CurrentDomain, ComputionType, WaterPoints3D)
                        else
                            call ComputeCartesian(SurfaceElevation, CurrentDomain, ComputionType)
                        endif
                    else
                        !Do Nothing
                    endif

                case default           
                    write(*,*)'UNKNOWN_ Domain type - review data file'
                    stop 'ComputeSZZ - Geometry - ERR01.'

            end select

            CurrentDomain => CurrentDomain%Next

        enddo

        !Stores the initial SZZ
        if (ComputionType == INITIALGEOMETRY) then
            call SetMatrixValueAllocatable(Me%Distances%InitialSZZ, Me%Size, Me%Distances%SZZ)
        endif

        !Disposes pointer to the Bathymetry
        call UngetGridData(Me%ObjTopography, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            stop 'ComputeSZZ - Geometry - ERR03'

       !UnGets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            stop 'ComputeSZZ - Geometry - ERR04'


    end subroutine ComputeSZZ


    !--------------------------------------------------------------------------
    !Computes SZZ for a FixSpacing Domain

    subroutine ComputeFixSpacing(Domain)

        !Parameter-------------------------------------------------------------
        type (T_Domain), pointer                :: Domain

        !Local-----------------------------------------------------------------
        real                                    :: LayerThickness
        integer                                 :: i, j, k, ILB, IUB, JLB, JUB

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        do k = Domain%LowerLayer, Domain%UpperLayer

            LayerThickness = Domain%TotalThickness * Domain%LayerThickness(k)

            do j = JLB, JUB
            do i = ILB, IUB
                Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k-1) - LayerThickness
                                                     
            enddo
            enddo

        enddo

    end subroutine ComputeFixSpacing

    !--------------------------------------------------------------------------
    !Computes SZZ for a FixSpacing Domain

    subroutine ComputeInitSediment(Domain)

        !Parameter-------------------------------------------------------------
        type (T_Domain), pointer                :: Domain

        !Local-----------------------------------------------------------------
        real                                    :: LayerThickness
        integer                                 :: i, j, k, ILB, IUB, JLB, JUB
        real, dimension(:, :), pointer          :: SedThickness
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        call GetGridData(Me%ObjTopography, SedThickness, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeInitSediment - ModuleGeometry - ERR01'

        do k = Domain%LowerLayer, Domain%UpperLayer

            do j = JLB, JUB
            do i = ILB, IUB
                
                LayerThickness            = SedThickness(i,j) * Domain%LayerThickness(k)

                Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k-1) - LayerThickness
                                                     
            enddo
            enddo

        enddo
        
        call UngetGridData(Me%ObjTopography, SedThickness, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeInitSediment - ModuleGeometry - ERR02'
            

    end subroutine ComputeInitSediment
    

    !--------------------------------------------------------------------------
    !Computes SZZ for a FixSediment Domain

    subroutine ComputeFixSediment (Domain, VerticalVelocity, DT)


        !Parameter-------------------------------------------------------------
        type (T_Domain), pointer                :: Domain
        real, dimension(:, :, :), pointer       :: VerticalVelocity
        real                                    :: DT

        !Local-----------------------------------------------------------------
        integer                                 :: i, j, k, ILB, IUB, JLB, JUB
        integer                                 :: STAT_CALL
        integer, dimension(:,:), pointer        :: WaterPoints2D
        real                                    :: Displacement
        real                                    :: IntegratedDisplacement

        !Gets a pointer to 2D WaterPoints
        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeFixSediment - Geometry - ERR01'


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        do j = JLB, JUB
        do i = ILB, IUB

            if (WaterPoints2D(i, j) == WaterPoint) then

                IntegratedDisplacement = 0.
                do k = Domain%LowerLayer, Domain%UpperLayer 
                    Displacement            = -1. * VerticalVelocity(i, j, k+1) * DT
                    IntegratedDisplacement  = IntegratedDisplacement + Displacement
                    Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k) + &
                                                         IntegratedDisplacement
                enddo

            endif

        enddo
        enddo
        

        !UnGets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            stop 'ComputeFixSediment - Geometry - ERR02'
        
    end subroutine ComputeFixSediment

    
    !--------------------------------------------------------------------------
    !Computes SZZ for a Cartesian Domain

    subroutine ComputeCartesian(SurfaceElevation, Domain, ComputionType, WaterPoints3D)

        !Parameter-------------------------------------------------------------
        real, dimension(:, :), pointer                 :: SurfaceElevation
        type (T_Domain), pointer                       :: Domain
        integer                                        :: ComputionType
        integer, dimension(:, :, :), optional, pointer :: WaterPoints3D

        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer        :: WaterPoints2D
        integer                                 :: i, j, k, ILB, IUB, JLB, JUB
        integer                                 :: STAT_CALL, kbottom
        real                                    :: TopDepth
        real                                    :: LayerThickness, AllmostZero_
        integer                                 :: ktop, KFloorZ

        nullify(WaterPoints2D)
        
        AllmostZero_ = 5e-4
        
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Gets a pointer to 2D WaterPoints
        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeCartesian - Geometry - ERR01'


cd1:    if (ComputionType == INITIALGEOMETRY) then

            kbottom = Domain%LowerLayer
     
            do j = JLB, JUB
            do i = ILB, IUB

                if (WaterPoints2D(i, j) == WaterPoint .and.                              &
                    Domain%UpperLayer >= Me%KFloor%Z(i, j)) then

                    !Domain Top Depth
                    if (Domain%ID == Me%LastDomain%ID) then
                        if (Domain%DomainType /= CartesianTop) then
                            TopDepth = 0.
                        else
                            TopDepth = -1. * SurfaceElevation(i, j)
                        endif
                    else
                        TopDepth = Domain%DomainDepth
                    endif

                    !Inits SZZ
                    Me%Distances%SZZ(i, j, Domain%UpperLayer) = TopDepth
                   
                    
                    !griflet: I must pass the structured data variables into scalar variables
                    !otherwise the model returns an access violation in openmp mode
                    !for no apparent reason...
                    ktop = Domain%UpperLayer
                    kbottom = Domain%LowerLayer
                    KFloorZ = Me%KFloor%Z(i,j)
                    if ( kbottom < KFloorZ ) kbottom = KFloorZ
                    
                    do k = ktop - 1, kbottom, -1

                        !The Layerthickness is now given in meters
                        !Frank - Jan 2001
                        LayerThickness = Domain%LayerThickness(k+1)
                        Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k+1) + LayerThickness

                        if ((Domain%DomainType == CartesianTop) .and. (k == kbottom)) then
                            if (abs(Me%Distances%SZZ(i, j, k) - Me%Distances%SZZ(i, j, k - 1)) <= AllmostZero_) then
                                Me%Distances%SZZ(i, j, k - 1) = FillValueDouble
                                Me%KFloor%Z(i, j) = Me%KFloor%Z(i, j) + 1
                                
                                if (present(WaterPoints3D)) then
                                    WaterPoints3D(i, j, k) = 0
                                endif
                            endif
                        endif

                   enddo

                endif

            enddo
            enddo

        endif  cd1

        if (Domain%UpperLayer == Me%WorkSize%KUB) then 

            do j = JLB, JUB
            do i = ILB, IUB

                if (WaterPoints2D(i, j) == WaterPoint)  then
                                
                      !In the cartesian domain the layer close to the surface have variable thickness
                      Me%Distances%SZZ(i, j, Me%WorkSize%KUB) = -1. * SurfaceElevation(i, j)

                endif

            enddo
            enddo

        endif


       !UnGets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            stop 'ComputeCartesian - Geometry - ERR02'



    end subroutine ComputeCartesian
    

    !--------------------------------------------------------------------------
    !Computes SZZ for a Harmonic Domain

    subroutine ComputeHarmonic(SurfaceElevation, Domain)

        !Parameter-------------------------------------------------------------
        real, dimension(:, :), pointer          :: SurfaceElevation
        type (T_Domain), pointer                :: Domain

        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer        :: WaterPoints2D
        integer                                 :: i, j, k, ILB, IUB, JLB, JUB, KUB
        integer                                 :: STAT_CALL, kbottom
        real                                    :: MinEsp, MeanLayerThickness, LayerThickness

        !Size
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        MinEsp = Domain%MinEsp

        !Gets a pointer to 2D WaterPoints
        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeHarmonic - Geometry - ERR01'

        !Cicle
        do j = JLB, JUB
        do i = ILB, IUB

            if (WaterPoints2D(i, j) == WaterPoint) then

                !Inits SZZ
                Me%Distances%SZZ(i, j, Domain%UpperLayer) = 0.0

                !Lowest layer in the domain
                kbottom = max(Domain%LowerLayer, Me%KFloor%Z(i, j))

                !Upper layer in the domain
                KUB     = Domain%UpperLayer

                !Calculates first SZZ approach
                do k = Domain%UpperLayer - 1, kbottom, -1

                    LayerThickness = Domain%LayerThickness(k+1)
                    Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k+1) + LayerThickness

                enddo

                !Sets SZZ(KUB) igual to waterlevel
                Me%Distances%SZZ(i, j, KUB)   = -1. * SurfaceElevation(i, j)

                !Calculates harmonic down
                k              = KUB - 1
                LayerThickness = Me%Distances%SZZ(i, j, k) - Me%Distances%SZZ(i, j, k+1)
                do while (LayerThickness < MinEsp .and. k >= kbottom)
                    Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k+1) + MinEsp
                    k              = k - 1
                    LayerThickness = Me%Distances%SZZ(i, j, k) - Me%Distances%SZZ(i, j, k+1)
                enddo

                !If the thickness of the cell close to the bottom is smaller then MinEsp
                !cells have to be redistributed. Two Options:
                !1. n Layers * MinEsp <= WaterColumn => Redistribute all cells equally
                !2. Increase Bottom Cell until every cell thickness > MinEsp
                if (Me%Distances%SZZ(i, j, kbottom-1) - Me%Distances%SZZ(i, j, kbottom) < MinEsp) then

                    if (Me%Distances%SZZ(i, j, kbottom-1)  - Me%Distances%SZZ(i, j, KUB) <= (KUB - kbottom + 1) * MinEsp) then
                        MeanLayerThickness = (Me%Distances%SZZ(i, j, kbottom-1)  -           &
                                              Me%Distances%SZZ(i, j, KUB)      ) /           &
                                              float(KUB - kbottom + 1)
                        do k = KUB - 1, kbottom, -1
                            Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k+1) + MeanLayerThickness
                        enddo
                    else
                        Me%Distances%SZZ(i, j, kbottom) = Me%Distances%SZZ(i, j, kbottom-1) - MinEsp
                        k = kbottom
                        do while (Me%Distances%SZZ(i, j, k) - Me%Distances%SZZ(i, j, k + 1) <  MinEsp)
                            Me%Distances%SZZ(i, j, k + 1) = Me%Distances%SZZ(i, j, k) - MinEsp 
                            k = k + 1
                        enddo
                    endif
                endif
            endif
        enddo
        enddo

        !UnGets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeHarmonic - Geometry - ERR02'


    end subroutine ComputeHarmonic

    !--------------------------------------------------------------------------
    !Computes SZZ for a Cartesian Domain at initial. Layers have to adapt to surface level

    subroutine ComputeCartesianNew(SurfaceElevation, Domain)

        !Parameter-------------------------------------------------------------
        real, dimension(:, :), pointer          :: SurfaceElevation
        type (T_Domain), pointer                :: Domain

        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer        :: WaterPoints2D
        integer                                 :: i, j, k, ILB, IUB, JLB, JUB, UpperLayer
        integer                                 :: STAT_CALL, LowerLayer
        real                                    :: MeanLayerThickness, LayerThickness
        real                                    :: MinThick, TotalMinThick

        !Size
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Gets a pointer to 2D WaterPoints
        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeCartesianNew - Geometry - ERR01'

        !Cycle
        do j = JLB, JUB
        do i = ILB, IUB

            if (WaterPoints2D(i, j) == WaterPoint.and.                              &
                Domain%UpperLayer >= Me%KFloor%Z(i, j)) then
                
                !Lowest layer in the domain
                LowerLayer = max(Domain%LowerLayer, Me%KFloor%Z(i, j))

                !Upper layer in the domain
                UpperLayer = Domain%UpperLayer
                
                !Inits SZZ
                Me%Distances%SZZ(i, j, UpperLayer) = 0.0

                !Calculates first SZZ approach. It will be rewritten in columns with
                !level negative (below Hydrographyc Zero)
                do k = UpperLayer - 1, LowerLayer, -1

                    LayerThickness = Domain%LayerThickness(k+1)
                    Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k+1) + LayerThickness

                enddo

                !Sets SZZ(KUB) to waterlevel
                Me%Distances%SZZ(i, j, UpperLayer)   = -1. * SurfaceElevation(i, j)
                
                !in negative levels, upper layer will be attached to surface but below need to
                !be adjusted
                if (SurfaceElevation(i, j) .lt. 0.) then
                
                    !Calculates cartesian down (as harmonic)
                    k              = UpperLayer - 1
                    LayerThickness = Me%Distances%SZZ(i, j, k) - Me%Distances%SZZ(i, j, k+1)
                    MinThick       = Domain%LayerThickness(k) * Domain%MinInitialLayerThickness 
                    do while (LayerThickness < MinThick .and. k >= LowerLayer)
                        Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k+1) + MinThick
                        k                         = k - 1
                        LayerThickness            = Me%Distances%SZZ(i, j, k) - Me%Distances%SZZ(i, j, k+1)
                        MinThick                  = Domain%LayerThickness(k) * Domain%MinInitialLayerThickness
                    enddo

                    !If the thickness of the cell close to the bottom is smaller then MinThick
                    !cells have to be redistributed. Two Options:
                    !1. TotalMinThick <= WaterColumn => Redistribute all cells equally
                    !2. else, put MinThick in all from bottom until every cell has thickness > MinThick
                    !In these cases velocity is forgotten but this only occurs in extreme events with
                    !very low water columns
                    LayerThickness   = (Me%Distances%SZZ(i,j,LowerLayer-1) - Me%Distances%SZZ(i,j,LowerLayer)) 
                    MinThick         = Domain%LayerThickness(LowerLayer) * Domain%MinInitialLayerThickness
                    if (LayerThickness < MinThick) then
                        
                        TotalMinThick = 0.
                        do k = LowerLayer, UpperLayer
                            TotalMinThick = TotalMinThick                                              & 
                                            + (Domain%LayerThickness(k) * Domain%MinInitialLayerThickness)
                        enddo

                        if (Me%Distances%SZZ(i, j, LowerLayer-1)  - Me%Distances%SZZ(i, j, UpperLayer) &
                             <= TotalMinThick) then
                            
                            MeanLayerThickness = (Me%Distances%SZZ(i, j, LowerLayer-1)  -              &
                                                  Me%Distances%SZZ(i, j, UpperLayer)      ) /          &
                                                  float(UpperLayer - LowerLayer + 1)
                            
                            do k = UpperLayer - 1, LowerLayer, -1
                                Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k+1) + MeanLayerThickness
                            enddo
                        
                        else
                            
                            k = LowerLayer
                            do while (Me%Distances%SZZ(i, j, k-1) - Me%Distances%SZZ(i, j, k) <  MinThick)
                                Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k-1) - MinThick 
                                k = k + 1
                                MinThick   = Domain%LayerThickness(k) * Domain%MinInitialLayerThickness
                            enddo
                        endif
                    endif
                endif            
            endif
        enddo
        enddo

        !UnGets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeCartesianNew - Geometry - ERR02'


    end subroutine ComputeCartesianNew

    !--------------------------------------------------------------------------
    !Computes SZZ for a Sigma Domain

    subroutine ComputeSigma(SurfaceElevation, Domain)

        !Parameter-------------------------------------------------------------
        real, dimension(:, :), pointer          :: SurfaceElevation
        type (T_Domain), pointer                :: Domain

        !Local-----------------------------------------------------------------
        integer, dimension(:, :), pointer       :: WaterPoints2D
        integer                                 :: i, j, k, ILB, IUB, JLB, JUB
        integer                                 :: LowerLayer, UpperLayer
        integer                                 :: LowerLayerTwoDomain, UpperLayerTwoDomain
        integer                                 :: STAT_CALL
        real                                    :: TopDepth, BottomDepth, DomainThickness
        real                                    :: MinorThickness, LayerThickness
        real                                    :: MeanLayerThicknessBelow, TwoDomainThickness
        real                                    :: A, B, C, s
        integer                                 :: CHUNK

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Gets a pointer to 2D WaterPoints
        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeSigma - Geometry - ERR01'


        !Lower and Upper Layer of the Domain
        LowerLayer = Domain%LowerLayer
        UpperLayer = Domain%UpperLayer

        if (MonitorPerformance) call StartWatch ("ModuleGeometry", "ComputeSigma")

        CHUNK = Chunk_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(i,j, TopDepth, BottomDepth, DomainThickness) &            
        !$OMP PRIVATE(LayerThickness, MeanLayerThicknessBelow, MinorThickness) &
        !$OMP PRIVATE(LowerLayerTwoDomain, UpperLayerTwoDomain, TwoDomainThickness) &
        !$OMP PRIVATE(s, A, B, C) 


        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj:    do j = JLB, JUB
doi:    do i = ILB, IUB

            
cd4:        if (WaterPoints2D(i, j) == WaterPoint .and. &
                LowerLayer >= Me%KFloor%Z(i, j) ) then
                
                !Upper Limit (m)
                if (Domain%ID < Me%LastDomain%ID) then
                    if (associated(Domain%Next)) then
                        if (Domain%Next%DomainType == SigmaTop) then
                            TopDepth = -1.* SurfaceElevation(i, j) + Domain%Next%TotalThickness
                        else
                            TopDepth = Domain%DomainDepth
                        endif
                    else
                        TopDepth = Domain%DomainDepth
                    endif
                else
                    TopDepth = -1.* SurfaceElevation(i, j)
                endif

                !Lower Limit
                BottomDepth = Me%Distances%SZZ(i, j, LowerLayer-1)
        
                !DomainThickness 
                DomainThickness = BottomDepth - TopDepth

                !Frank
                ! The test below is now done in the hydrodynamic module. It can't be done here,
                ! due to the Volumebalance of the cell. 
                ! Paulo Chambel Leitão
                !OLD WAY - is correct because in this way the volume is never zero. 
                !if the volume is zero all the operations that devided the volume will have a overflow
!                if (DomainThickness <= Me%WaterColumn%ZMin) then
!                    write(*,*)'Low domain thickness', i, j
!                    DomainThickness = Me%WaterColumn%ZMin
!                endif

                ! Paulo Chambel Leitão
                !NEW WAY - is wrong            
                ! NOTE: If the Waterlevel in a certain cell is lower then the WaterColumn%ZMin, water
                ! can just enter to the cell. So it will be no "danger" (I hope) to calculate the
                ! volumes of the cell, since the waterlevel is positive. 
                ! If the WaterLevel is negative, the volumes and areas are set to zero, instead of 
                ! setting them to WaterColumn%ZMin (like until now), so the erro is smaller. 
                ! Frank May 99
                !if (DomainThickness < 0.) DomainThickness = 0.

                MinorThickness = 1.e5
                s = -1.
                !Normal computation
                do k = LowerLayer, UpperLayer

                    if (Domain%RomsDistortion) then
                    !ROMS Stretching Function
                        !Song, Y. and D. B. Haidvogel, 1994: A semi-implicit ocean circulation model 
                        !using a generalized topography-following coordinate system, J. Comp. Phys., 115 (1), 228-244. (PDF)
                        !https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
                        s = s + Domain%LayerThickness(k)
                        
                        A = (1.-Domain%theta_b)*sinh(Domain%theta_s * s) / sinh(Domain%theta_s)
                        
                        B = tanh(Domain%theta_s*(s+0.5))/tanh(0.5*Domain%theta_s)
                        
                        C = A  + Domain%theta_b * 0.5*(B-1.)
                        
                        Me%Distances%SZZ(i, j, k) = - (Domain%Hc * s + (DomainThickness - Domain%Hc) * C) + Topdepth           
                    
                        cycle
                    endif                
                    
                
                    LayerThickness = DomainThickness * Domain%LayerThickness(k)
                    Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k-1) - LayerThickness
                    if (LayerThickness < MinorThickness) MinorThickness = LayerThickness
                enddo
                
                 Me%Distances%SZZ(i, j, UpperLayer) = TopDepth
                
                !Redefinition of the current domain and the lower domain in the case that the lower
                !one is of the type FixSpacing and, at the same time (the DomainThickness is less 
                !or equal to the WaterColumn%ZMin) or (the Minor layer thichness is less or equal to
                !the mean layer thickness of the fixspacing domain). 
                !In the case this happens, the thickness of the two layers is added and then
                !evenly distributed
                !
                !The if's are separated, to avoid that the array indexes run out of bound
                !

cd1:            if (Domain%ID > 1) then
cd2:                if (Domain%Prev%DomainType == FixSpacing) then
                        MeanLayerThicknessBelow = Domain%Prev%TotalThickness /           &
                                                  float(Domain%Prev%NumberOfLayers)

                        !The total thickness of the two domains
                        LowerLayerTwoDomain = Domain%Prev%LowerLayer
                        UpperLayerTwoDomain = Domain%UpperLayer
                        TwoDomainThickness  = Me%Distances%SZZ(i, j, LowerLayerTwoDomain-1) &
                                            - Me%Distances%SZZ(i, j, UpperLayerTwoDomain) 
                        
cd3:                    if (TwoDomainThickness <= (Me%WaterColumn%ZMin + &
                            Domain%Prev%TotalThickness)) then

                            !Redefinition of the FixSpacing Domain
                            do k = Domain%Prev%LowerLayer, Domain%Prev%UpperLayer
                                LayerThickness = TwoDomainThickness * &
                                                 Domain%Prev%LayerThickness(k) / 2.
                                Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k-1) - LayerThickness
                            enddo

                            !Redefinition of the Sigma Domain
                            do k = Domain%LowerLayer, Domain%UpperLayer
                                LayerThickness = TwoDomainThickness * &
                                                 Domain%LayerThickness(k) / 2.
                                Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k-1) - LayerThickness
                            enddo

                        endif cd3
                    endif cd2   
                endif cd1

            endif cd4

        enddo doi
        enddo doj
        !$OMP END DO

        !$OMP END PARALLEL
        
        if (MonitorPerformance) call StopWatch ("ModuleGeometry", "ComputeSigma")       

       !UnGets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            stop 'ComputeSigma - Geometry - ERR02'

    end subroutine ComputeSigma


    !--------------------------------------------------------------------------
    !Computes SZZ for a Lagrangian Domain
    subroutine ComputeLagrangian(SurfaceElevation, VerticalVelocity, DT_Waterlevel, Domain)

        !Parameter-------------------------------------------------------------
        real, dimension(:, :), pointer          :: SurfaceElevation
        real, dimension(:, :, :), pointer       :: VerticalVelocity
        real, intent(in)                        :: DT_Waterlevel
        type (T_Domain), pointer                :: Domain

        !Local-----------------------------------------------------------------
        integer                                 :: i, j, k, ILB, IUB, JLB, JUB
        integer                                 :: LowerLayer, UpperLayer
        integer, dimension(:, :), pointer       :: WaterPoints2D
        real                                    :: TopDepth, BottomDepth, DomainThickness
        real                                    :: MinimalThickness, GridMovementDump
        real                                    :: MaximumDisplacementUp, MaximumSzzUp, MaximumUndumpSzzUp
        real                                    :: MaximumDisplacementDown, MaximumSzzDown, MaximumUndumpSzzDown
        real                                    :: FreeSZZ, NewSZZ, MaximumDisplacement, FreeGridVelocity, UpSZZ
        real                                    :: DisplacementLimit
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Shorten variables name 

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Dumping
        GridMovementDump = Domain%GridMovementDump

        !Upper Layer of the Domain
        UpperLayer       = Domain%UpperLayer

        !Gets a pointer to 2D WaterPoints
        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeLagrangian - Geometry - ERR10'

doj:    do j = JLB, JUB
doi:    do i = ILB, IUB

cd0 :       if (WaterPoints2D(i, j) == WaterPoint) then

                !Upper Limit (m)
                if (Domain%ID < Me%LastDomain%ID) then
                    TopDepth = Domain%DomainDepth
                else
                    TopDepth = -1.* SurfaceElevation(i, j)
                endif

                !Lower Layer
                LowerLayer = max(Domain%LowerLayer, Me%KFloor%Z(i, j))

                !Lower Limit
                BottomDepth = Me%Distances%SZZ(i, j, LowerLayer-1)
            
                !DomainThickness 
                DomainThickness = BottomDepth - TopDepth
            
                !See ComputeSigma for more information about the next statement
                !This correction is now done in the hydrodynamic module - Frank
!                if (DomainThickness <= Me%WaterColumn%ZMin) &
!                    DomainThickness = Me%WaterColumn%ZMin

                !Resizes the layers, according to the vertical velocity
                !This process is done from the BOTTOM to the TOP

                !MinimalThickness is now defined as a percentage of the initial thickness of each layer  
                MinimalThickness   = Domain%MinEvolveLayerThickness


                !DisplacementLimit is the maximum displacement that the model allow cell 
                !faces to move vertically in meters
                DisplacementLimit  = Domain%DisplacementLimit

                do k = LowerLayer, UpperLayer-1

                    if (k+1 == Me%WorkSize%KUB) then
                        !Free surface and do not have a rigid reference 
                        UpSZZ = Me%Distances%SZZ(i, j, k+1)
                    else
                        !Need a rigid reference 
                        UpSZZ = Me%Distances%InitialSZZ(i, j, k+1)
                    endif

                    MaximumDisplacementUp   = (Me%Distances%InitialSZZ(i, j, k) -   &
                                               UpSZZ)* &
                                               (1.-MinimalThickness)/2.



                    MaximumDisplacementDown = (Me%Distances%InitialSZZ(i, j, k-1) - &
                                               Me%Distances%InitialSZZ(i, j, k)) *  &
                                              (1.-MinimalThickness)/2.

                    MaximumDisplacement     = min(MaximumDisplacementUp,                     &
                                                  MaximumDisplacementDown,                   &
                                                  DisplacementLimit)


                    MaximumSzzUp            = Me%Distances%InitialSZZ(i, j, k) -    &
                                              MaximumDisplacement

                    MaximumUndumpSzzUp      = Me%Distances%InitialSZZ(i, j, k) -    &
                                              MaximumDisplacement   *                        &
                                              (1.0 - GridMovementDump)

                    MaximumSzzDown          = Me%Distances%InitialSZZ(i, j, k) +    &
                                              MaximumDisplacement

                    MaximumUndumpSzzDown    = Me%Distances%InitialSZZ(i, j, k) +    &
                                              MaximumDisplacement     *                      &
                                              (1.0 - GridMovementDump)

                    !For the Layer k the Vertical Velocity to consider is the velocity in 
                    !k+1(Upper Face) once the vertical velocity has a different index than 
                    !the SZZ
                    FreeGridVelocity = VerticalVelocity(i, j, k+1)

                    
                    FreeSZZ = Me%Distances%SZZ(i, j, k) -                           &
                              FreeGridVelocity * DT_Waterlevel

                    if (FreeSZZ < MaximumUndumpSzzUp) then !Dumping to the maximum level

                        NewSZZ = MaximumUndumpSzzUp - (MaximumUndumpSzzUp - MaximumSzzUp) *  &
                                 (1. - 1./(MaximumUndumpSzzUp - FreeSZZ + 1))

                    elseif (FreeSZZ > MaximumUndumpSzzDown) then !Dumping to the minimum level

                        NewSZZ = MaximumUndumpSzzDown + (MaximumSzzDown -                    &
                                 MaximumUndumpSzzDown) * (1. - 1./(FreeSZZ -                 &
                                 MaximumUndumpSzzDown + 1.))
                    else

                        NewSZZ = FreeSZZ

                    endif

                    Me%Distances%SZZ(i, j, k) = NewSZZ

                    !Nudging
                    if (associated(Me%ExternalVar%DecayTime)) then
                        Me%Distances%SZZ(i, j, k) =  Me%Distances%SZZ(i, j, k) +        &
                                                    (Me%Distances%InitialSZZ(i, j, k) - &
                                                     Me%Distances%SZZ(i, j, k))       * &
                                                     DT_Waterlevel / Me%ExternalVar%DecayTime(i, j, k)
                                                    
                    endif

                enddo
                !Actualizes SZZ at the surface if the domain exists (DomainThickness>=0) 
                !in this I,J cell 
                if (DomainThickness>=0)                                                 &
                    Me%Distances%SZZ(i, j, UpperLayer) = TopDepth

            endif cd0

        enddo doi
        enddo doj

        !UnGets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ComputeLagrangian - Geometry - ERR50'


    end subroutine ComputeLagrangian

    !--------------------------------------------------------------------------
    !Computes SZZ for a Lagrangian Domain - from up
    subroutine ComputeLagrangianNew(SurfaceElevation, VerticalVelocity, DT_Waterlevel, Domain)

        !Parameter-------------------------------------------------------------
        real, dimension(:, :), pointer          :: SurfaceElevation
        real, dimension(:, :, :), pointer       :: VerticalVelocity
        real, intent(in)                        :: DT_Waterlevel
        type (T_Domain), pointer                :: Domain

        !Local-----------------------------------------------------------------
        integer                                 :: i, j, k, ILB, IUB, JLB, JUB
        integer                                 :: LowerLayer, UpperLayer
        integer, dimension(:, :), pointer       :: WaterPoints2D
!        real   , dimension(:, :,:), pointer     :: OldSZZ
        real                                    :: TopDepth, BottomDepth, DomainThickness
        real                                    :: MinimalThickness
        real                                    :: DisplacementLimit
        real                                    :: CompressionLimitSZZ, ExpansionLimitSZZ, DisplacementLimitSZZ
        real                                    :: MaxDisplacementSZZ
        real                                    :: FreeSZZ, NewSZZ, FreeGridVelocity
        integer                                 :: STAT_CALL
        real                                    :: MinThick, TotalMinThick, MeanLayerThickness
        real                                    :: Thickness

        !Begin-----------------------------------------------------------------

        !Shorten variables name 

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Upper Layer of the Domain
        UpperLayer         = Domain%UpperLayer

        !MinimalThickness is now defined as a percentage of the initial thickness of each layer  
        MinimalThickness   = Domain%MinEvolveLayerThickness

        !DisplacementLimit is the maximum displacement that the model allow cell 
        !faces to move vertically in meters in one iteration
        DisplacementLimit  = Domain%DisplacementLimit        
        
        !cartesian domain limits are computed once (not changing)
        if (Domain%DomainType == Cartesian .and. .not. Me%LagrangianLimitsComputed) then
            do k = Domain%LowerLayer, UpperLayer
                Domain%LayerMinThickness(k) = Domain%LayerThickness(k) * (1. - Domain%MinEvolveLayerThickness)
                Domain%LayerMaxThickness(k) = Domain%LayerThickness(k) * (1. + Domain%MinEvolveLayerThickness)
                Me%LagrangianLimitsComputed = .true.
            enddo
        endif

        !Gets a pointer to 2D WaterPoints
        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeLagrangianNew - ModuleGeometry - ERR10'

doj:    do j = JLB, JUB
doi:    do i = ILB, IUB

cd0 :       if (WaterPoints2D(i, j) == WaterPoint .and.                              &
                    Domain%UpperLayer >= Me%KFloor%Z(i, j)) then

                !Upper Limit (m)
                if (Domain%ID < Me%LastDomain%ID) then
                    TopDepth = Domain%DomainDepth
                else
                    TopDepth = -1.* SurfaceElevation(i, j)
                endif

                !Lower Layer
                LowerLayer       = max(Domain%LowerLayer, Me%KFloor%Z(i, j))

                !Lower Limit
                BottomDepth      = Me%Distances%SZZ(i, j, LowerLayer-1)
            
                !DomainThickness 
                DomainThickness  = BottomDepth - TopDepth
                
                !Actualizes SZZ at the surface if the domain exists (DomainThickness>=0) 
                !in this I,J cell 
                if (DomainThickness >= 0)  &
                    Me%Distances%SZZ(i, j, UpperLayer) = TopDepth

                !actualize layer limits in case of sigma domain.
                if (Domain%DomainType == Sigma) then
                    do k = LowerLayer, UpperLayer
                        Domain%LayerMinThickness(k) = Domain%LayerThickness(k) * DomainThickness &
                                                      * (1. - MinimalThickness)
                        Domain%LayerMaxThickness(k) = Domain%LayerThickness(k) * DomainThickness &
                                                      * (1. + MinimalThickness)
                    enddo 
                endif

do1 :           do k = UpperLayer-1, LowerLayer, -1
                    
                    
                    !For the Layer k the Vertical Velocity to consider is the velocity in 
                    !k+1(Upper Face) once the vertical velocity has a different index than 
                    !the SZZ
                    FreeGridVelocity  = VerticalVelocity(i, j, k+1)

                    FreeSZZ           = Me%Distances%SZZ(i, j, k) -           &
                                        FreeGridVelocity * DT_Waterlevel
                    
                    NewSZZ            = FreeSZZ
                                       
                    !verify if the displacement was not over the limits (compressed cell min thick
                    !and expanded cell maximum thickness
                    if (NewSZZ > Me%Distances%SZZ(i, j, k)) then !lowered layer, compresses below layer and expands above layer
                    
                        CompressionLimitSZZ     = Me%Distances%SZZ(i,j,k-1) - &
                                                  Domain%LayerMinThickness(k)
                                              
                        ExpansionLimitSZZ       = Me%Distances%SZZ(i,j,k+1) + &
                                                  Domain%LayerMaxThickness(k+1)
                        
                        DisplacementLimitSZZ    = Me%Distances%SZZ(i, j, k) + DisplacementLimit
                        
                        
                        MaxDisplacementSZZ      = min (CompressionLimitSZZ,   &
                                                       ExpansionLimitSZZ,     &
                                                       DisplacementLimitSZZ)
      
                        
                        if (NewSZZ > MaxDisplacementSZZ) NewSZZ = MaxDisplacementSZZ
 
                        
                    elseif  (NewSZZ < Me%Distances%SZZ(i, j, k)) then !raised layer, compresses above layer and expands below layer
                                               
                        CompressionLimitSZZ     = Me%Distances%SZZ(i,j,k+1) + &
                                                  Domain%LayerMinThickness(k+1)
                                            
                        ExpansionLimitSZZ       = Me%Distances%SZZ(i,j,k-1) - &
                                                  Domain%LayerMaxThickness(k)
                                             
                        DisplacementLimitSZZ    = Me%Distances%SZZ(i, j, k) - DisplacementLimit

                        
                        MaxDisplacementSZZ      = max (CompressionLimitSZZ,   &
                                                       ExpansionLimitSZZ,     &
                                                       DisplacementLimitSZZ)
                        
                        
                        if (NewSZZ < MaxDisplacementSZZ) NewSZZ = MaxDisplacementSZZ
                        
                        
                    endif
                    
                    Me%Distances%SZZ(i, j, k)     = NewSZZ
                    
                    !Layers can not have thickness lower than minimum. push them down 
                    !Avoid very thin or negative thickness. May happen when surface water 
                    !is lowering but water is travelling horizontal in the cell so the
                    !above process does not change enough layers position
                    if (Me%Distances%SZZ(i, j, k) - Me%Distances%SZZ(i, j, k+1) <  &
                        Domain%LayerMinThickness(k+1) ) then
                        
                        Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k+1) +  &
                                                    Domain%LayerMinThickness(k+1)
                    endif
                    
                    
                    !Nudging - layers need to have relaxation to initial layer definition
                    if (associated(Me%ExternalVar%DecayTime)) then
                        Me%Distances%SZZ(i, j, k) =  Me%Distances%SZZ(i, j, k) +        &
                                                    (Me%Distances%InitialSZZ(i, j, k) - &
                                                     Me%Distances%SZZ(i, j, k)) *       &
                                                     DT_Waterlevel / Me%ExternalVar%DecayTime(i, j, k)
                                                    
                    endif

                enddo do1


                !If the thickness of the cell close to the bottom is smaller then MinThick
                !cells have to be redistributed. Two Options:
                !1. TotalMinThick <= WaterColumn => Redistribute all cells equally
                !2. else, put MinThick in all from bottom until every cell has thickness > MinThick
                !In these cases velocity is forgotten but this only occurs in extreme events with
                !very low water columns
                Thickness   = (Me%Distances%SZZ(i,j,LowerLayer-1) - Me%Distances%SZZ(i,j,LowerLayer)) 
                MinThick    = Domain%LayerMinThickness(LowerLayer)
                if (Thickness < MinThick) then
                    
                    TotalMinThick = 0.
                    do k = LowerLayer, UpperLayer
                        TotalMinThick = TotalMinThick + Domain%LayerMinThickness(k)
                    enddo

                    if (Me%Distances%SZZ(i, j, LowerLayer-1)  - Me%Distances%SZZ(i, j, UpperLayer) &
                         <= TotalMinThick) then
                        
                        MeanLayerThickness = (Me%Distances%SZZ(i, j, LowerLayer-1)  -              &
                                              Me%Distances%SZZ(i, j, UpperLayer)      ) /          &
                                              float(UpperLayer - LowerLayer + 1)
                        
                        do k = UpperLayer - 1, LowerLayer, -1
                            Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k+1) + MeanLayerThickness
                        enddo
                    
                    else
                        
                        k = LowerLayer
                        do while (Me%Distances%SZZ(i, j, k-1) - Me%Distances%SZZ(i, j, k) <  MinThick)
                            Me%Distances%SZZ(i, j, k) = Me%Distances%SZZ(i, j, k-1) - MinThick 
                            k = k + 1
                            MinThick   = Domain%LayerMinThickness(k)
                        enddo
                    endif
                endif                    

            endif cd0

        enddo doi
        enddo doj

        !UnGets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ComputeLagrangianNew - ModuleGeometry - ERR200'
        
!        deallocate (OldSZZ)
!        nullify (OldSZZ)
        
    end subroutine ComputeLagrangianNew

    !--------------------------------------------------------------------------
    !Reads SZZ and VolumeZOld

    subroutine ReadGeometryBin(GeometryID, Unit, STAT)

        !Parameter-------------------------------------------------------------
        integer                                     :: GeometryID
        integer, intent(in)                         :: Unit
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_
        integer                                     :: i, j, k, ILB, IUB, JLB, JUB, KLB, KUB

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            !Size
            ILB = Me%Size%ILB
            IUB = Me%Size%IUB

            JLB = Me%Size%JLB
            JUB = Me%Size%JUB

            KLB = Me%Size%KLB
            KUB = Me%Size%KUB

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                read(Unit) Me%Distances%SZZ(i, j, k)
            enddo
            enddo
            enddo

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                read(Unit) Me%Distances%InitialSZZ(i, j, k)
            enddo
            enddo
            enddo

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                read(Unit) Me%Distances%DWZ(i, j, k)
            enddo
            enddo
            enddo

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                read(Unit) Me%Volumes%VolumeZOld(i, j, k)
            enddo
            enddo
            enddo

            STAT_ = SUCCESS_

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT)) &
            STAT = STAT_

    end subroutine ReadGeometryBin


    !--------------------------------------------------------------------------
    !Reads SZZ and VolumeZOld

    subroutine WriteGeometryBin(GeometryID, Unit, STAT)

        !Parameter-------------------------------------------------------------
        integer                                     :: GeometryID
        integer, intent(in)                         :: Unit
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_
        integer                                     :: i, j, k, ILB, IUB, JLB, JUB, KLB, KUB

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            !Size
            ILB = Me%Size%ILB
            IUB = Me%Size%IUB

            JLB = Me%Size%JLB
            JUB = Me%Size%JUB

            KLB = Me%Size%KLB
            KUB = Me%Size%KUB

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                write(Unit) Me%Distances%SZZ(i, j, k)
            enddo
            enddo
            enddo

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                write(Unit) Me%Distances%InitialSZZ(i, j, k)
            enddo
            enddo
            enddo

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                write(Unit) Me%Distances%DWZ(i, j, k)
            enddo
            enddo
            enddo

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                write(Unit)Me%Volumes%VolumeZOld(i, j, k)
            enddo
            enddo
            enddo

            STAT_ = SUCCESS_

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT)) &
            STAT = STAT_

    end subroutine WriteGeometryBin
    
    !----------------------------------------------------------------------------


        !--------------------------------------------------------------------------
    !Reads SZZ and VolumeZOld

    subroutine WriteGeometryHDF(GeometryID, ObjHDF5, STAT)

        !Parameter-------------------------------------------------------------
        integer                                     :: GeometryID
        integer, intent(in)                         :: ObjHDF5
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then
            
            !WorkSize
            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            KLB = Me%WorkSize%KLB
            KUB = Me%WorkSize%KUB
            
        
            !Sets limits for next write operations
            call HDF5SetLimits(ObjHDF5, ILB, IUB, JLB, JUB, KLB-1, KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteGeometryHDF - Geometry - ERR01'

            call HDF5WriteData(ObjHDF5,"/Geometry", "VerticalZ", "m",               &
                               Array3D = GetPointer(Me%Distances%SZZ),              &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteGeometryHDF - Geometry - ERR02'

            call HDF5WriteData(ObjHDF5,"/Geometry", "InitialSZZ", "m",              &
                               Array3D = GetPointer(Me%Distances%InitialSZZ),       &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteGeometryHDF - Geometry - ERR03'

            !Sets limits for next write operations
            call HDF5SetLimits(ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteGeometryHDF - Geometry - ERR01'


            call HDF5WriteData(ObjHDF5,"/Geometry", "DWZ", "m",                     &
                               Array3D = GetPointer(Me%Distances%DWZ),              &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteGeometryHDF - Geometry - ERR04'

            call HDF5WriteData(ObjHDF5,"/Geometry", "VolumeZ", "m3",                &
                               Array3D = GetPointer(Me%Volumes%VolumeZOld),         &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteGeometryHDF - Geometry - ERR05'
            
            call HDF5WriteData(ObjHDF5,"/Geometry", "KTop", "-",                    &
                               Array2D = GetPointer(Me%KTop%Z),                     &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteGeometryHDF - Geometry - ERR06'

            call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteGeometryHDF - Geometry - ERR07'

            STAT_ = SUCCESS_

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT)) &
            STAT = STAT_

    end subroutine WriteGeometryHDF
    
    !----------------------------------------------------------------------------

    subroutine ReadGeometryHDF(GeometryID, HDF5FileName, MasterOrSlave, STAT)

        !Parameter-------------------------------------------------------------
        integer                                 :: GeometryID
        character(len=*)                        :: HDF5FileName
        logical, intent(in ), optional          :: MasterOrSlave
        integer, intent(out), optional          :: STAT

        !Local-----------------------------------------------------------------
        type (T_Size2D)                         :: WindowLimitsJI
        integer                                 :: ready_   
        integer                                 :: STAT_
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: IUW, JUW, ILW, JLW        
        integer                                 :: Imax, Jmax, Kmax 
        integer                                 :: STAT_CALL
        integer                                 :: ObjHDF5 = 0
        integer                                 :: HDF5_READ
        logical                                 :: MasterOrSlave_
        integer, dimension(:,:  ), pointer      :: Aux2DInt
        real,    dimension(:,:,:), pointer      :: Aux3DReal
        real(8), dimension(:,:,:), pointer      :: Aux3DR8        

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            if (present(MasterOrSlave)) then
                MasterOrSlave_ = MasterOrSlave
            else
                MasterOrSlave_ = .false.
            endif
            
            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
            
            !Opens HDF File
            call ConstructHDF5      (ObjHDF5,                                            &
                                     trim(HDF5FileName),                                 &
                                     HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGeometryHDF - Geometry - ERR10'
            
            call GetHDF5ArrayDimensions (HDF5ID = ObjHDF5, GroupName = "/Grid",     &
                                        ItemName = "WaterPoints3D",                    &
                                        Imax = Imax, Jmax = Jmax, Kmax = Kmax, STAT = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGeometryHDF - Geometry - ERR20'
                                        
            ILB = Me%WorkSize%ILB 
            IUB = Me%WorkSize%IUB 
            
            JLB = Me%WorkSize%JLB 
            JUB = Me%WorkSize%JUB 

            KLB = Me%WorkSize%KLB 
            KUB = Me%WorkSize%KUB 
            
            
    ifMS:   if (MasterOrSlave) then
    
                call GetDDecompWorkSize2D(HorizontalGridID = Me%ObjHorizontalGrid, &
                                          WorkSize         = WindowLimitsJI,       &
                                          STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadGeometryHDF - Geometry - ERR30'
                
                ILW = WindowLimitsJI%ILB
                IUW = WindowLimitsJI%IUB

                JLW = WindowLimitsJI%JLB
                JUW = WindowLimitsJI%JUB
                                                      
            else ifMS

                ILW = ILB 
                IUW = IUB

                JLW = JLB 
                JUW = JUB 

            endif ifMS                                            
            
            if (ILB < 1   ) stop 'ReadGeometryHDF - Geometry - ERR40'
            if (IUB > Imax) stop 'ReadGeometryHDF - Geometry - ERR50'
            
            if (JLB < 1   ) stop 'ReadGeometryHDF - Geometry - ERR60'
            if (JUB > Jmax) stop 'ReadGeometryHDF - Geometry - ERR70'

            if (KLB < 1   ) stop 'ReadGeometryHDF - Geometry - ERR80'
            if (KUB > Kmax) stop 'ReadGeometryHDF - Geometry - ERR90'
            
            allocate(Aux3DReal(ILW:IUW,JLW:JUW,KLB-1:KUB))
            
            !Sets limits for next read operations
            call HDF5SetLimits(ObjHDF5, ILW, IUW, JLW, JUW, KLB-1, KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGeometryHDF - Geometry - ERR100'
            
            call HDF5ReadWindow(HDF5ID        = ObjHDF5,                                &                    
                                GroupName     = "/Geometry",                            &
                                Name          = "VerticalZ",                            &
                                Array3D       = Aux3DReal,                              &
                                OffSet3       = 0,                                      &
                                STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGeometryHDF - Geometry - ERR110'
            
            Me%Distances%SZZ(ILB:IUB,JLB:JUB,KLB-1:KUB) = Aux3DReal(ILW:IUW,JLW:JUW,KLB-1:KUB)

            call HDF5ReadWindow(HDF5ID        = ObjHDF5,                                &                    
                                GroupName     = "/Geometry",                            &
                                Name          = "InitialSZZ",                           &
                                Array3D       = Aux3DReal,                              &
                                OffSet3       = 0,                                      &
                                STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGeometryHDF - Geometry - ERR120'
            
            Me%Distances%InitialSZZ(ILB:IUB,JLB:JUB,KLB-1:KUB) = Aux3DReal(ILW:IUW,JLW:JUW,KLB-1:KUB)

            !Sets limits for next read operations
            call HDF5SetLimits(ObjHDF5, ILW, IUW, JLW, JUW, KLB, KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGeometryHDF - Geometry - ERR130'
            
            deallocate(Aux3DReal)            
            allocate  (Aux3DReal(ILW:IUW,JLW:JUW,KLB:KUB))
            allocate  (Aux3DR8  (ILW:IUW,JLW:JUW,KLB:KUB))
            allocate  (Aux2DInt (ILW:IUW,JLW:JUW))

            call HDF5ReadWindow(HDF5ID        = ObjHDF5,                                &
                                GroupName     = "/Geometry",                            &
                                Name          = "DWZ",                                  &
                                Array3D       = Aux3DReal,                              &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGeometryHDF - Geometry - ERR140'
            
            Me%Distances%DWZ(ILB:IUB,JLB:JUB,KLB:KUB) = Aux3DReal(ILW:IUW,JLW:JUW,KLB:KUB)

            call HDF5ReadWindow(HDF5ID        = ObjHDF5,                                &
                                GroupName     = "/Geometry",                            &
                                Name          = "VolumeZ",                              &
                                Array3D       = Aux3DR8,                                &
                                STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGeometryHDF - Geometry - ERR150'
            
            Me%Volumes%VolumeZOld(ILB:IUB,JLB:JUB,KLB:KUB) = Aux3DR8(ILW:IUW,JLW:JUW,KLB:KUB)            

            call HDF5ReadWindow(HDF5ID        = ObjHDF5,                                &
                                GroupName     = "/Geometry",                            &
                                Name          = "KTop",                                 &
                                Array2D       = Aux2DInt,                               &
                                STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGeometryHDF - Geometry - ERR160'
            
            Me%KTop%Z(ILB:IUB,JLB:JUB) = Aux2DInt(ILW:IUW,JLW:JUW)

            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGeometryHDF - Geometry - ERR170'
            
            deallocate  (Aux3DReal)
            deallocate  (Aux3DR8  )
            deallocate  (Aux2DInt )

            STAT_ = SUCCESS_

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT)) &
            STAT = STAT_

    end subroutine ReadGeometryHDF
    
    !----------------------------------------------------------------------------



    !----------------------------------------------------------------------------

    subroutine ComputeZCellCenter

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        real, dimension(:, :, :), pointer       :: SZZ 
        real, dimension(:, :, :), pointer       :: ZCellCenter 
        integer                                 :: I, J, K
        integer                                 :: CHUNK

        !------------------------------------------------------------------------

        SZZ         => Me%Distances%SZZ
        ZCellCenter => Me%Distances%ZCellCenter

        if (MonitorPerformance) call StartWatch ("ModuleGeometry", "ComputeZCellCenter")

        CHUNK = Chunk_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
        !$OMP PARALLEL PRIVATE(I,J,K)            

do1 :   do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(STATIC, CHUNK)
do2 :   do J = Me%WorkSize%JLB, Me%WorkSize%JUB
do3 :   do I = Me%WorkSize%ILB, Me%WorkSize%IUB

            ZCellCenter(I,J,K) = -1.0 * (SZZ(I,J,K-1) + SZZ(I,J,K)) / 2.0

        end do do3
        end do do2
        !$OMP END DO
        end do do1

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleGeometry", "ComputeZCellCenter")

        nullify(SZZ        )
        nullify(ZCellCenter)

        !------------------------------------------------------------------------

    end subroutine ComputeZCellCenter

#ifdef _USE_SEQASSIMILATION

    !--------------------------------------------------------------------------

    subroutine CopyGeometryDistances(GeometryID, SZZ, DWZ, DUZ, DVZ, DZZ,   &
                                     ZCellCenter, STAT)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )              :: GeometryID
        real,    dimension(:, :, :), pointer        :: SZZ, DWZ, DUZ, DVZ
        real,    dimension(:, :, :), pointer        :: DZZ
        real,    dimension(:, :, :), pointer        :: ZCellCenter
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer :: ready_        

        !Local-----------------------------------------------------------------
        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_) 

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call SetMatrixValue(Me%Distances%SZZ, Me%Size, SZZ)

            call SetMatrixValue(Me%Distances%DWZ, Me%Size, DWZ)

            call SetMatrixValue(Me%Distances%DUZ, Me%Size, DUZ)

            call SetMatrixValue(Me%Distances%DVZ, Me%Size, DVZ)

            call SetMatrixValue(Me%Distances%DZZ, Me%Size, DZZ)

            call SetMatrixValue(Me%Distances%ZCellCenter, Me%Size, ZCellCenter)

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine CopyGeometryDistances

    !--------------------------------------------------------------------------

    subroutine CopyGeometryAreas(GeometryID, AreaU, AreaV, STAT)

        !Arguments-------------------------------------------------------------
        integer,              intent(IN )           :: GeometryID
        real, dimension(:, :, :), pointer           :: AreaU, AreaV
        integer, optional,    intent(OUT)           :: STAT

        !External--------------------------------------------------------------
        integer :: ready_        

        !Local-----------------------------------------------------------------
        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call SetMatrixValue(Me%Areas%AreaU, Me%Size, AreaU)

            call SetMatrixValue(Me%Areas%AreaV, Me%Size, AreaV)

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine CopyGeometryAreas

    !--------------------------------------------------------------------------

    subroutine CopyGeometryVolumes(GeometryID, VolumeZ, VolumeZOld,         &
                                   VolumeU, VolumeV, STAT)

        !Arguments-------------------------------------------------------------
        integer,              intent(IN )           :: GeometryID
        real(8), dimension(:, :, :), pointer        :: VolumeZ, VolumeU
        real(8), dimension(:, :, :), pointer        :: VolumeV, VolumeZOld
        integer, optional,    intent(OUT)           :: STAT

        !External--------------------------------------------------------------
        integer :: ready_        

        !Local-----------------------------------------------------------------
        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call SetMatrixValue(Me%Volumes%VolumeZ, Me%Size, VolumeZ)

            call SetMatrixValue(Me%Volumes%VolumeZOld, Me%Size, VolumeZOld)

            call SetMatrixValue(Me%Volumes%VolumeU, Me%Size, VolumeU)

            call SetMatrixValue(Me%Volumes%VolumeV, Me%Size, VolumeV)

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine CopyGeometryVolumes

    !--------------------------------------------------------------------------

    subroutine CopyGeometryWaterColumn(GeometryID, WaterColumnZ, WaterColumnU,  &
                                       WaterColumnV, STAT)

        !Arguments-------------------------------------------------------------
        integer,              intent(IN )           :: GeometryID
        real,    dimension(:, :), pointer           :: WaterColumnZ
        real,    dimension(:, :), pointer           :: WaterColumnU, WaterColumnV
        integer, optional,    intent(OUT)           :: STAT

        !External--------------------------------------------------------------
        integer :: ready_        

        !Local-----------------------------------------------------------------
        integer :: STAT_              !Auxiliar local variable
        type (T_Size2D)                             :: Size2D

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            Size2D%ILB = Me%Size%ILB
            Size2D%IUB = Me%Size%IUB
            Size2D%JLB = Me%Size%JLB
            Size2D%JUB = Me%Size%JUB

            call SetMatrixValue(Me%WaterColumn%Z, Size2D, WaterColumnZ)

            call SetMatrixValue(Me%WaterColumn%U, Size2D, WaterColumnU)

            call SetMatrixValue(Me%WaterColumn%V, Size2D, WaterColumnV)

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine CopyGeometryWaterColumn

    !--------------------------------------------------------------------------

#endif _USE_SEQASSIMILATION

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetGeometrySize(GeometryID , Size, WorkSize, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        type (T_Size3D), optional                   :: Size
        type (T_Size3D), optional                   :: WorkSize
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    

cd1 :   if ((ready_ .EQ. IDLE_ERR_) .OR. &
           (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(Size       )) Size        = Me%Size
            if (present(WorkSize   )) WorkSize    = Me%WorkSize

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetGeometrySize
    
    !--------------------------------------------------------------------------

    subroutine GetGeometryMinWaterColumn(GeometryID, MinWaterColumn, STAT)

        !Parameter-------------------------------------------------------------
        integer                                     :: GeometryID
        real                                        :: MinWaterColumn
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

             MinWaterColumn = Me%WaterColumn%ZMin

             STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGeometryMinWaterColumn

    !--------------------------------------------------------------------------

    subroutine GetGeometryDistances(GeometryID, SZZ, DZZ, DWZ, DUZ, DVZ, DZI, DZE,        &
                                    ZCellCenter, DWZ_Xgrad, DWZ_Ygrad, ActualTime, STAT)

        !Parameter-------------------------------------------------------------
        integer                                     :: GeometryID
        real, dimension(:,:,:), pointer, optional   :: SZZ, DZZ, DWZ, DUZ, DVZ, DWZ_Xgrad, DWZ_Ygrad 
        real, dimension(:,:,:), pointer, optional   :: DZI, DZE, ZCellCenter
        type(T_Time),                    optional   :: ActualTime
        integer, intent(OUT),            optional   :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !Verifies if the distances are up to date
            if (present(ActualTime)) then
                if (ActualTime .ne. Me%ActualTime) &
                    stop 'GetGeometryDistances - Geometry - ERR01'
            endif

            !SZZ
            if (present(SZZ)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                SZZ => Me%Distances%SZZ
            endif

            !DZZ
            if (present(DZZ)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                DZZ => Me%Distances%DZZ
            endif

            !DWZ
            if (present(DWZ)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                DWZ => Me%Distances%DWZ
            endif


            !DUZ
            if (present(DUZ)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                DUZ => Me%Distances%DUZ
            endif

            !DVZ
            if (present(DVZ)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                DVZ => Me%Distances%DVZ
            endif

            !DZI
            if (present(DZI)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                DZI => Me%Distances%DZI
            endif

            !DZE
            if (present(DZE)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                DZE => Me%Distances%DZE
            endif


            !ZCellCenter
            if (present(ZCellCenter)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                ZCellCenter => Me%Distances%ZCellCenter
            endif

            !DWZ_Xgrad
            if (present(DWZ_Xgrad)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                DWZ_Xgrad => Me%Distances%DWZ_Xgrad
            endif

            !DWZ_Ygrad
            if (present(DWZ_Ygrad)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                DWZ_Ygrad => Me%Distances%DWZ_Ygrad
            endif

            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGeometryDistances

    !--------------------------------------------------------------------------

    subroutine GetGeometryAreas(GeometryID, AreaU, AreaV, ActualTime, STAT)

        !Parameter-------------------------------------------------------------
        integer                                     :: GeometryID
        real, dimension(:, :, :), pointer, optional :: AreaU, AreaV
        type (T_Time), optional                     :: ActualTime
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !Verifies if the areas are up to date
            if (present(ActualTime)) then
                if (ActualTime .ne. Me%ActualTime) &
                    stop 'GetGeometryAreas - Geometry - ERR01'
            endif

            !AreaU
            if (present(AreaU)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                AreaU => Me%Areas%AreaU
            endif

            !AreaV
            if (present(AreaV)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                AreaV => Me%Areas%AreaV
            endif

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGeometryAreas

    !--------------------------------------------------------------------------

    subroutine GetGeometryVolumes(GeometryID, VolumeZ, VolumeU, VolumeV, &
                                  VolumeW, VolumeZOld, ActualTime, STAT)

        !Parameter-------------------------------------------------------------
        integer                                         :: GeometryID
        real(8), dimension(:, :, :), pointer, optional  :: VolumeZ, VolumeU, VolumeV, VolumeW, VolumeZOld
        type (T_Time), optional                         :: ActualTime
        integer, intent(out), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: ready_   
        integer                                         :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !Verifies if the volumes are up to date
            if (present(ActualTime)) then
                if (ActualTime .ne. Me%ActualTime) &
                    stop 'GetGeometryVolumes - Geometry - ERR01'
            endif

            !VolumeZ
            if (present(VolumeZ)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                VolumeZ => Me%Volumes%VolumeZ
            endif

            !VolumeU
            if (present(VolumeU)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                VolumeU => Me%Volumes%VolumeU
            endif

            !VolumeV
            if (present(VolumeV)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                VolumeV => Me%Volumes%VolumeV
            endif

            !VolumeW
            if (present(VolumeW)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                VolumeW => Me%Volumes%VolumeW
            endif

            !VolumeZ
            if (present(VolumeZOld)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                VolumeZOld => Me%Volumes%VolumeZOld
            endif

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGeometryVolumes


    !--------------------------------------------------------------------------

    subroutine GetGeometryKTop(GeometryID, KTopZ, STAT)

        !Parameter-------------------------------------------------------------
        integer                                     :: GeometryID
        integer, dimension(:, :), pointer, optional :: KTopZ
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !KTopZ
            if (present(KTopZ)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                KTopZ => Me%KTop%Z
            endif

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGeometryKTop
    
    !--------------------------------------------------------------------------

    subroutine GetGeometryKFloor(GeometryID, Z, U, V, Domain, STAT)

        !Parameter-------------------------------------------------------------
        integer                                     :: GeometryID
        integer, dimension(:, :), pointer, optional :: Z, U, V, Domain
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !Z
            if (present(Z)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                Z => Me%KFloor%Z
            endif

            !U
            if (present(U)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                U => Me%KFloor%U
            endif

            !V
            if (present(V)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                V => Me%KFloor%V
            endif

            !Domain
            if (present(Domain)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                Domain => Me%KFloor%Domain
            endif

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGeometryKFloor

    !--------------------------------------------------------------------------

    subroutine GetGeometryWaterColumn(GeometryID, WaterColumn , WaterColumnU,           &
                                                   WaterColumnV, ActualTime, STAT)

        !Parameter-------------------------------------------------------------
        integer                                     :: GeometryID
        real, dimension(:, :), pointer, optional    :: WaterColumn
        real, dimension(:, :), pointer, optional    :: WaterColumnU, WaterColumnV
        type (T_Time), optional                     :: ActualTime
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !Verifies if the areas are up to date
            if (present(ActualTime)) then
                if (ActualTime .ne. Me%ActualTime) &
                    stop 'GetGeometryWaterColumn - Geometry - ERR01'
            endif

            if (present(WaterColumn)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                WaterColumn => Me%WaterColumn%Z
            endif

            if (present(WaterColumnU)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                WaterColumnU => Me%WaterColumn%U
            endif

            if (present(WaterColumnV)) then
                call Read_Lock(mGEOMETRY_, Me%InstanceID)
                WaterColumnV => Me%WaterColumn%V
            endif

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGeometryWaterColumn

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetGeometryInputFile(GeometryID, InputFile, STAT)

        !Parameter-------------------------------------------------------------
        integer                                     :: GeometryID
        character(len=*)                            :: InputFile
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            InputFile = Me%InputFile

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetGeometryInputFile

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    integer function GetLayer4Level(GeometryID, i, j, DepthLevel, STAT)

        !Parameter-------------------------------------------------------------
        integer, intent (IN)                        :: GeometryID, i, j
        real   , intent (IN)                        :: DepthLevel
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_
        integer                                     :: k, kbottom, KUB

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            KUB     = Me%WorkSize%KUB
            kbottom = Me%KFloor%Z(i, j)
            
i0:         if (kbottom >=1 .and. kbottom <= KUB) then

i1:             if      (DepthLevel < Me%Distances%SZZ(i,j,KUB)) then

                    GetLayer4Level = Me%WorkSize%KUB

                else if (DepthLevel > Me%Distances%SZZ(i,j,kbottom-1)) then i1

                    GetLayer4Level = kbottom

                else i1
                
                    do k = kbottom, KUB
                        if (DepthLevel <= Me%Distances%SZZ(i,j,k-1) .and.                    &
                            DepthLevel >= Me%Distances%SZZ(i,j,k  )) exit
                    enddo

                    GetLayer4Level = k

                endif i1
                
            else i0
            
                GetLayer4Level = KUB
            
            endif i0

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end function GetLayer4Level

        !----------------------------------------------------------------------

    subroutine UnGetGeometry2Dinteger(GeometryID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        integer, pointer, dimension(:,:)            :: Array
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_UnLock(mGEOMETRY_, Me%InstanceID, "UnGetGeometry2Dinteger")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UnGetGeometry2Dinteger

    !--------------------------------------------------------------------------
    
    subroutine UnGetGeometry2Dreal(GeometryID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        real, pointer, dimension(:,:)               :: Array
        integer, optional, intent(OUT)              :: STAT



        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_


        call Ready(GeometryID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_UnLock(mGEOMETRY_, Me%InstanceID, "UnGetGeometry2Dreal")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UnGetGeometry2Dreal

    !--------------------------------------------------------------------------
    
    subroutine UnGetGeometry3Dreal4(GeometryID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        real(4), pointer, dimension(:,:,:)          :: Array
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mGEOMETRY_, Me%InstanceID, "UnGetGeometry3Dreal4")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UnGetGeometry3Dreal4

    !--------------------------------------------------------------------------
    
    subroutine UnGetGeometry3Dreal8(GeometryID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        real(8), pointer, dimension(:,:,:)          :: Array
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_UnLock(mGEOMETRY_, Me%InstanceID, "UnGetGeometry3Dreal8")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UnGetGeometry3Dreal8

    !--------------------------------------------------------------------------

#ifdef _USE_SEQASSIMILATION

    subroutine SetGeometryDistances(GeometryID, SZZ, DWZ, DUZ, DVZ, DZZ,    &
                                    ZCellCenter, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        real,    dimension(:, :, :), pointer        :: SZZ, DWZ, DUZ, DVZ
        real,    dimension(:, :, :), pointer        :: DZZ
        real,    dimension(:, :, :), pointer        :: ZCellCenter
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_) then

            Me%Distances%SZZ => SZZ

            Me%Distances%DWZ => DWZ

            Me%Distances%DUZ => DUZ

            Me%Distances%DVZ => DVZ

            Me%Distances%DZZ => DZZ

            Me%Distances%ZCellCenter => ZCellCenter
            !CAUTION: pointing to an external variable/memory space!

            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if

        if (present(STAT))STAT = STAT_
            
    end subroutine SetGeometryDistances

    !--------------------------------------------------------------------------

    subroutine SetGeometryAreas(GeometryID, AreaU, AreaV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        real, dimension(:, :, :), pointer           :: AreaU, AreaV
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_) then

            Me%Areas%AreaU => AreaU

            Me%Areas%AreaV => AreaV

            !CAUTION: pointing to an external variable/memory space!

            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if

        if (present(STAT))STAT = STAT_
            
    end subroutine SetGeometryAreas

    !--------------------------------------------------------------------------

    subroutine SetGeometryVolumes(GeometryID, VolumeZ, VolumeZOld,          &
                                  VolumeU, VolumeV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        real(8), dimension(:, :, :), pointer        :: VolumeZ, VolumeU
        real(8), dimension(:, :, :), pointer        :: VolumeV, VolumeZOld
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_) then

            Me%Volumes%VolumeZ => VolumeZ

            Me%Volumes%VolumeZOld => VolumeZOld

            Me%Volumes%VolumeU => VolumeU

            Me%Volumes%VolumeV => VolumeV

            !CAUTION: pointing to an external variable/memory space!

            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if

        if (present(STAT))STAT = STAT_
            
    end subroutine SetGeometryVolumes

    !--------------------------------------------------------------------------

    subroutine SetGeometryWaterColumn(GeometryID, WaterColumnZ, WaterColumnU,   &
                                      WaterColumnV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        real, dimension(:, :), pointer              :: WaterColumnZ
        real, dimension(:, :), pointer              :: WaterColumnU, WaterColumnV
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_) then

            Me%WaterColumn%Z => WaterColumnZ

            Me%WaterColumn%U => WaterColumnU

            Me%WaterColumn%V => WaterColumnV

            !CAUTION: pointing to an external variable/memory space!

            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if

        if (present(STAT))STAT = STAT_
            
    end subroutine SetGeometryWaterColumn

    !--------------------------------------------------------------------------

    subroutine ReSetGeometry(GeometryID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                      :: GeometryID
        integer,            optional, intent(OUT)    :: STAT

        !Local-------------------------------------------------------------------
        integer                                      :: ready_          
        integer                                      :: STAT_    

        !------------------------------------------------------------------------

        !This is a compilation of sets (one for each variable) for internal memory spaces

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_) 

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            Me%Distances%SZZ => Me%AuxPointer%SZZ

            Me%Distances%DWZ => Me%AuxPointer%DWZ

            Me%Distances%DUZ => Me%AuxPointer%DUZ

            Me%Distances%DVZ => Me%AuxPointer%DVZ

            Me%Distances%DZZ => Me%AuxPointer%DZZ

            Me%Distances%ZCellCenter => Me%AuxPointer%ZCellCenter

            Me%Areas%AreaU => Me%AuxPointer%AreaU

            Me%Areas%AreaV => Me%AuxPointer%AreaV

            Me%Volumes%VolumeZ => Me%AuxPointer%VolumeZ

            Me%Volumes%VolumeZOld => Me%AuxPointer%VolumeZOld

            Me%Volumes%VolumeU => Me%AuxPointer%VolumeU

            Me%Volumes%VolumeV => Me%AuxPointer%VolumeV

            Me%WaterColumn%Z => Me%AuxPointer%WaterColumnZ

            Me%WaterColumn%U => Me%AuxPointer%WaterColumnU

            Me%WaterColumn%V => Me%AuxPointer%WaterColumnV

            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if cd1


        if (present(STAT))STAT = STAT_

    end subroutine ReSetGeometry

#endif _USE_SEQASSIMILATION

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillGeometry(GeometryID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GeometryID
        integer, intent(out), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_
        type (T_Domain), pointer                    :: CurrentDomain, DomainToDelete
        integer                                     :: nUsers

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mGEOMETRY_,  Me%InstanceID)

            if (nUsers == 0) then

                !Deassociates External Instances
                nUsers = DeassociateInstance (mGRIDDATA_,       Me%ObjTopography)
                if (nUsers == 0) stop 'KillGeometry - Geometry - ERR10'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillGeometry - Geometry - ERR20'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillGeometry - Geometry - ERR30'

                !Deallocates local variables
                call DeallocateVariables

                !Delete all domains
                CurrentDomain => Me%FirstDomain
                do while (associated(CurrentDomain))
                    DomainToDelete => CurrentDomain
                    CurrentDomain  => CurrentDomain%Next
                    call DeleteDomain(DomainToDelete)
                enddo

                if (Me%Areas%Impermeability) then

                    deallocate (Me%Areas%Coef_U )
                    deallocate (Me%Areas%CoefX_U)

                    deallocate (Me%Areas%Coef_V )
                    deallocate (Me%Areas%CoefX_V)

                endif

                call DeallocateInstance

                GeometryID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
              
        !----------------------------------------------------------------------

    end subroutine KillGeometry

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Geometry), pointer          :: AuxObjGeometry
        type (T_Geometry), pointer          :: PreviousObjGeometry

        !Updates pointers
        if (Me%InstanceID == FirstGeometry%InstanceID) then
            FirstGeometry => FirstGeometry%Next
        else
            PreviousObjGeometry => FirstGeometry
            AuxObjGeometry      => FirstGeometry%Next
            do while (AuxObjGeometry%InstanceID /= Me%InstanceID)
                PreviousObjGeometry => AuxObjGeometry
                AuxObjGeometry      => AuxObjGeometry%Next
            enddo

            !Now update linked list
            PreviousObjGeometry%Next => AuxObjGeometry%Next

        endif
            
        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------

    subroutine DeallocateVariables

        !Parameter-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                         :: STATUS

        !Deallocates T_Volumes
        if (allocated(Me%Volumes%VolumeZ)) then
            deallocate (Me%Volumes%VolumeZ, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR10'
        endif

        if (allocated(Me%Volumes%VolumeU)) then
            deallocate (Me%Volumes%VolumeU, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR20'
        endif


        if (allocated(Me%Volumes%VolumeV)) then
            deallocate (Me%Volumes%VolumeV, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR30'
        endif

        if (allocated(Me%Volumes%VolumeW)) then
            deallocate (Me%Volumes%VolumeW, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR40'
        endif


        if (allocated(Me%Volumes%VolumeZOld)) then
            deallocate (Me%Volumes%VolumeZOld, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR50'
        endif

        !Allocate T_Areas
        if (allocated(Me%Areas%AreaU)) then
            deallocate (Me%Areas%AreaU, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR60'
        endif

        if (allocated(Me%Areas%AreaV)) then
            deallocate (Me%Areas%AreaV, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR70'
        endif

        !Allocate T_Distances
        if (allocated(Me%Distances%SZZ)) then
            deallocate (Me%Distances%SZZ, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR80'
        endif

        if (allocated(Me%Distances%DZZ)) then
            deallocate (Me%Distances%DZZ, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR90'
        endif

        if (allocated(Me%Distances%DWZ)) then
            deallocate (Me%Distances%DWZ, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR100'
        endif

        if (allocated(Me%Distances%DUZ)) then
            deallocate (Me%Distances%DUZ, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR110'
        endif

        if (allocated(Me%Distances%DVZ)) then
            deallocate (Me%Distances%DVZ, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR120'
        endif

        if (allocated(Me%Distances%DZI)) then
            deallocate (Me%Distances%DZI, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR130'
        endif

        if (allocated(Me%Distances%DZE)) then
            deallocate (Me%Distances%DZE, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR140'
        endif


        if (allocated(Me%Distances%InitialSZZ)) then
            deallocate (Me%Distances%InitialSZZ, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR150'
        endif


        if (allocated(Me%Distances%ZCellCenter)) then
            deallocate (Me%Distances%ZCellCenter, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR160'
        endif


        if (allocated(Me%Distances%DWZ_Xgrad)) then
            deallocate (Me%Distances%DWZ_Xgrad, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR170'
        endif

        if (allocated(Me%Distances%DWZ_Ygrad)) then
            deallocate (Me%Distances%DWZ_Ygrad, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR180'
        endif

        !deallocate T_KFloor
        if (allocated(Me%KFloor%Z)) then
            deallocate (Me%KFloor%Z, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR190'
        endif

        if (allocated(Me%KFloor%U)) then
            deallocate (Me%KFloor%U, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR200'
        endif

        if (allocated(Me%KFloor%V)) then
            deallocate (Me%KFloor%V, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR210'
        endif

        if (allocated(Me%KFloor%Domain)) then
            deallocate (Me%KFloor%Domain, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR220'
        endif

        !deallocate T_KTop
        if (allocated(Me%KTop%Z)) then
            deallocate (Me%KTop%Z, stat = STATUS)
            if (STATUS /= SUCCESS_) stop 'DeallocateVariables - Geometry - ERR230'
        endif
        

    end subroutine DeallocateVariables

    !--------------------------------------------------------------------------

    subroutine DeleteDomain (DomainToDelete)

        !Arguments-------------------------------------------------------------
        type (T_Domain), pointer                    :: DomainToDelete

        !Local-----------------------------------------------------------------
        type (T_Domain), pointer                    :: CurrentDomain
        type (T_Domain), pointer                    :: Previous
        type (T_Domain), pointer                    :: Next


        CurrentDomain => Me%FirstDomain
        do while (associated(CurrentDomain))
            if (CurrentDomain%ID == DomainToDelete%ID) then

                Previous       => CurrentDomain%Prev
                Next           => CurrentDomain%Next

                !Updates foward pointer
                if (associated(CurrentDomain%Prev)) then
                    Previous%Next  => Next
                endif

                !Updates backward pointer
                if (associated(CurrentDomain%Next)) then
                    Next%Prev      => Previous
                endif

                !Updates first domain
                if (DomainToDelete%ID == Me%FirstDomain%ID) then
                    Me%FirstDomain => Next
                endif

                !Updates last domain
                if (DomainToDelete%ID == Me%LastDomain%ID) then
                    Me%LastDomain => Previous
                endif

                !Deletes domain
                if (associated(DomainToDelete%LayerThickness)) then
                    deallocate(DomainToDelete%LayerThickness)
                    nullify   (DomainToDelete%LayerThickness)
                endif

                nullify   (DomainToDelete%Next)
                nullify   (DomainToDelete%Prev)
                deallocate(DomainToDelete)
                nullify   (DomainToDelete)

                return

            endif

            CurrentDomain => CurrentDomain%Next
        enddo

    end subroutine DeleteDomain

#ifdef _USE_SEQASSIMILATION

    !--------------------------------------------------------------------------

    subroutine NullifyGeometryStatePointer(GeometryID, STAT)

        !Arguments-------------------------------------------------------------

        integer,               intent(IN ) :: GeometryID
        integer, optional,     intent(OUT) :: STAT

        !External--------------------------------------------------------------

        integer :: ready_   

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GeometryID, ready_) 

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            nullify(Me%AuxPointer%SZZ)

            nullify(Me%AuxPointer%DWZ)

            nullify(Me%AuxPointer%DUZ)

            nullify(Me%AuxPointer%DVZ)

            nullify(Me%AuxPointer%DZZ)

            nullify(Me%AuxPointer%ZCellCenter)

            nullify(Me%AuxPointer%AreaU)

            nullify(Me%AuxPointer%AreaV)

            nullify(Me%AuxPointer%VolumeZ)

            nullify(Me%AuxPointer%VolumeZOld)

            nullify(Me%AuxPointer%VolumeU)

            nullify(Me%AuxPointer%VolumeV)

            nullify(Me%AuxPointer%WaterColumnZ)

            nullify(Me%AuxPointer%WaterColumnU)

            nullify(Me%AuxPointer%WaterColumnV)

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine NullifyGeometryStatePointer

#endif _USE_SEQASSIMILATION


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine Ready (ObjGeometry_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjGeometry_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjGeometry_ID > 0) then
            call LocateObjGeometry(ObjGeometry_ID)
            ready_ = VerifyReadLock (mGEOMETRY_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjGeometry (ObjGeometryID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjGeometryID

        !Local-----------------------------------------------------------------

        Me => FirstGeometry
        do while (associated (Me))
            if (Me%InstanceID == ObjGeometryID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'Geometry - LocateObjGeometry - ERR01'

    end subroutine LocateObjGeometry

    !--------------------------------------------------------------------------

end module ModuleGeometry

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

