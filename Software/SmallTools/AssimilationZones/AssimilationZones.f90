!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Batim Filter
! PROJECT       : AssimilationZones
! MODULE        : -
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2004
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Program to create assimilation coefs files
!------------------------------------------------------------------------------

!Data file - default name 'AssimilationZones.dat' (must be placed in working directory)

!   IN_BATIM                    : char              -           !Original bathymetry file name
!   OUT_FILE                    : char              -           !Grid Data file with the assimilation cells
!   SURROUND                    : 0/1              [1]          !Define assimilation zone in grid outer limits
!   NUMBER_CELLS                : int               -           !Number of assimilation cells in XX direction
!   ASSIMILATION_COEFS          : real vector       -           !Assimilation decay coefs (outside-inside)
!   TYPE_Z                      : 0/1              [1]          !Output refers to cells center
!   TYPE_U                      : 0/1              [0]          !Output refers to cells faces in XX direction
!   TYPE_V                      : 0/1              [0]          !Output refers to cells faces in YY direction
!   NO_ASSIMILATION_ZONE        : 0/1              [0]          !Use a polygon file to define where assimilation 
!                                                               !is not to be defined
!   NO_ASSIMILATION_FILE        : char              -           !Path to polygon file
!   DEFAULTVALUE                : real              -           !Default value assumed for assimilation coef.
!   WRITE_AS_GRID_DATA          : 0/1              [1]          !Write GridData regular format
!   WRITE_IN_INDICES_FORMAT     : 0/1              [1]          !Write as [i],[j],[value] format
!   SURFACE_FORCING             : 0/1              [0]          !If the user wants to assimilate values for the surface layer
!   SURFACE_COEF                : real              -           !Assimilation decay coef for the surface layer
!   LAYERS_NUMBER               : integer           -           !Number of layers
!   MINIMUM_DOMINATES           : 0/1              [1]          !By default the lower decay times dominate over the higher values
!   WEST                        : 0/1              [1]          !Decay values in the west boundary
!   EAST                        : 0/1              [1]          !Decay values in the east boundary
!   SOUTH                       : 0/1              [1]          !Decay values in the south boundary
!   NOTTH                       : 0/1              [1]          !Decay values in the north boundary


program AssimilationZones

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData     
    use ModuleHorizontalGrid                                     
    use ModuleGridData  
    use ModuleHorizontalMap
    use ModuleDrawing
             
    implicit none

    !Instances
    integer                                     :: ObjEnterData             = 0
    integer                                     :: ObjBathymetry            = 0
    integer                                     :: ObjHorizontalGrid        = 0
    integer                                     :: ObjHorizontalMap         = 0

    !Dummy variable time
    type (T_Time  )                             :: CurrentTime

    !Size
    type (T_Size2D)                             :: Size
    type (T_Size2D)                             :: WorkSize

    !Files
    character(len=PathLength)                   :: BathymetryFile
    character(len=PathLength)                   :: OutputFile
    character(len=PathLength)                   :: DataFile                 = 'AssimilationZones.dat'
    character(len=PathLength)                   :: PolygonFile

    !Input variables
    logical                                     :: Surround                 = .true.
    logical                                     :: SurfaceForcing           = .false.
    logical                                     :: TypeZ                    = .true.
    logical                                     :: TypeU                    = .false.
    logical                                     :: TypeV                    = .false.
    logical                                     :: DefineNoAssimilationZone = .false.
    logical                                     :: WriteAsGridData          = .true.
    logical                                     :: WriteAsGridIndices       = .true.

    real                                        :: DefaultCoef              = 1e32
    real                                        :: SurfaceCoef              = null_real
    integer                                     :: nCells                   = null_int
    integer                                     :: LayersNumber             = null_int
    real,           dimension(:  ), pointer     :: AssimilationCoefs
    logical                                     :: MinimumDominantes        = .true.
    logical                                     :: West                     = .true.
    logical                                     :: East                     = .true.
    logical                                     :: South                    = .true.
    logical                                     :: North                    = .true.
    logical                                     :: WriteDefaults            = .false.


    !Working variables
    integer,        dimension(:,:), pointer     :: WaterPoints2D
    real,           dimension(:,:), pointer     :: Bathymetry
    real,           dimension(:,:), pointer     :: Zones2D_U
    real,           dimension(:,:), pointer     :: Zones2D_V
    real,           dimension(:,:), pointer     :: Zones2D_Z
    real,           dimension(:,:,:), pointer   :: Zones3D_U
    real,           dimension(:,:,:), pointer   :: Zones3D_V
    real,           dimension(:,:,:), pointer   :: Zones3D_Z
    real,           dimension(:,:), pointer     :: Zone_XX
    real,           dimension(:,:), pointer     :: Zone_YY
    type(T_PointF ),dimension(:,:), pointer     :: GridPoints
    type(T_PointF ),dimension(:,:), pointer     :: GridFacesX
    type(T_PointF ),dimension(:,:), pointer     :: GridFacesY
    type(T_Polygon),                pointer     :: NoAssimilationZone
    integer                                     :: OutputUnit

    !External variables
    integer                                     :: CoordType
    real,    dimension(:,:), pointer            :: XX_IE, YY_IE
    integer                                     :: GEOG, UTM, SIMPLE_GEOG
    integer                                     :: MIL_PORT, GRID_COORD, NLRD

    !Begin---------------------------------------------------------------------

    call StartUpMohid("Mohid Define Assimilation Zones")

    call ReadOptions

    call ConstructDomain

    call AllocateVariables
	
    call DefineZones

    if (SurfaceForcing) then
        call DefineSurfaceForcing
    endif

    call WriteAssimilationZones

    call DeallocateVariables

    call KillAssimilationZones

    contains

    !--------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Reading options...'
        write(*,*)

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR10'

        call GetData(BathymetryFile,                                    &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'IN_BATIM',                         &
                     ClientModule = 'AssimilationZones',                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR20'

        call GetData(OutputFile,                                        &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'OUT_FILE',                         &
                     ClientModule = 'AssimilationZones',                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR30'

        call GetData(Surround,                                          &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'SURROUND',                         &
                     ClientModule = 'AssimilationZones',                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR40'

        if(Surround)then

            call GetData(nCells,                                        &
                         ObjEnterData, iflag   ,                        &
                         SearchType   = FromFile,                       &
                         keyword      = 'NUMBER_CELLS',                 &
                         ClientModule = 'AssimilationZones',            &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR50'


            nullify (AssimilationCoefs)
            allocate(AssimilationCoefs(1:nCells))

            call GetData(AssimilationCoefs,                             &
                         ObjEnterData, iflag   ,                        &
                         SearchType   = FromFile,                       &
                         keyword      = 'ASSIMILATION_COEFS',           &
                         ClientModule = 'AssimilationZones',            &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR60'

        end if

        call GetData(SurfaceForcing,                                    &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'SURFACE_FORCING',                  &
                     ClientModule = 'AssimilationZones',                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR70'

        if(SurfaceForcing)then

            call GetData(SurfaceCoef,                                   &
                         ObjEnterData, iflag   ,                        &
                         SearchType   = FromFile,                       &
                         keyword      = 'SURFACE_COEF',                 &
                         ClientModule = 'AssimilationZones',            &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR80'

            call GetData(LayersNumber,                                  &
                         ObjEnterData, iflag   ,                        &
                         SearchType   = FromFile,                       &
                         keyword      = 'LAYERS_NUMBER',                &
                         ClientModule = 'AssimilationZones',            &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR90'



        end if

        call GetData(TypeZ,                                             &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'TYPE_Z',                           &
                     Default      = OFF,                                &
                     ClientModule = 'AssimilationZones',                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR100'

        call GetData(TypeU,                                             &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'TYPE_U',                           &
                     Default      = OFF,                                &
                     ClientModule = 'AssimilationZones',                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR110'


        call GetData(TypeV,                                             &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'TYPE_V',                           &
                     Default      = OFF,                                &
                     ClientModule = 'AssimilationZones',                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR120'


        call GetData(DefineNoAssimilationZone,                          &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'NO_ASSIMILATION_ZONE',             &
                     Default      = OFF,                                &
                     ClientModule = 'AssimilationZones',                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR130'

        if(DefineNoAssimilationZone)then

            call GetData(PolygonFile,                                       &
                         ObjEnterData, iflag   ,                            &
                         SearchType   = FromFile,                           &
                         keyword      = 'NO_ASSIMILATION_FILE',             &
                         ClientModule = 'AssimilationZones',                &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR140'


            call New(NoAssimilationZone, trim(PolygonFile))

        endif
        
        call GetData(DefaultCoef,                                           &
                     ObjEnterData, iflag   ,                                &
                     SearchType   = FromFile,                               &
                     keyword      = 'DEFAULTVALUE',                         &
                     Default      = 1e32,                                   &
                     ClientModule = 'AssimilationZones',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR150'

        call GetData(WriteAsGridData,                                       &
                     ObjEnterData, iflag   ,                                &
                     SearchType   = FromFile,                               &
                     keyword      = 'WRITE_AS_GRID_DATA',                   &
                     Default      = ON,                                     &
                     ClientModule = 'AssimilationZones',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR160'

        if (.not. TypeZ .and. WriteAsGridData) then
            write(*,*) 'The program is only able to write TypeZ data in GridData format'
            stop 'ReadOptions - AssimilationZones - ERR165'
        endif

        call GetData(WriteAsGridIndices,                                    &
                     ObjEnterData, iflag   ,                                &
                     SearchType   = FromFile,                               &
                     keyword      = 'WRITE_IN_INDICES_FORMAT',              &
                     Default      = ON,                                     &
                     ClientModule = 'AssimilationZones',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR170'


        call GetData(MinimumDominantes,                                     &
                     ObjEnterData, iflag   ,                                &
                     SearchType   = FromFile,                               &
                     keyword      = 'MINIMUM_DOMINATES',                    &
                     Default      = ON,                                     &
                     ClientModule = 'AssimilationZones',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR180'

        call GetData(West,                                                  &
                     ObjEnterData, iflag   ,                                &
                     SearchType   = FromFile,                               &
                     keyword      = 'WEST',                                 &
                     Default      = ON,                                     &
                     ClientModule = 'AssimilationZones',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR190'

        call GetData(East,                                                  &
                     ObjEnterData, iflag   ,                                &
                     SearchType   = FromFile,                               &
                     keyword      = 'EAST',                                 &
                     Default      = ON,                                     &
                     ClientModule = 'AssimilationZones',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR200'


        call GetData(South,                                                 &
                     ObjEnterData, iflag   ,                                &
                     SearchType   = FromFile,                               &
                     keyword      = 'SOUTH',                                &
                     Default      = ON,                                     &
                     ClientModule = 'AssimilationZones',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR210'


        call GetData(North,                                                 &
                     ObjEnterData, iflag   ,                                &
                     SearchType   = FromFile,                               &
                     keyword      = 'NORTH',                                &
                     Default      = ON,                                     &
                     ClientModule = 'AssimilationZones',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR220'


        call GetData(WriteDefaults,                                         &
                     ObjEnterData, iflag   ,                                &
                     SearchType   = FromFile,                               &
                     keyword      = 'WRITE_DEFAULTS',                       &
                     Default      = OFF,                                    &
                     ClientModule = 'AssimilationZones',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR230'

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - AssimilationZones - ERR240'

    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    subroutine ConstructDomain

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
    
        write(*,*)
        write(*,*)"Constructing horizontal grid..."
        write(*,*)

        call SetDate(CurrentTime, 2004, 1, 1, 0, 0, 0)

        !Horizontal Grid
        call ConstructHorizontalGrid(HorizontalGridID = ObjHorizontalGrid,          &
                                     DataFile         = BathymetryFile,             &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - AssimilationZones - ERR01'


        write(*,*)
        write(*,*)"Constructing bathymetry..."
        write(*,*)


        !Horizontal Grid Data - Water Column (Bathymetry)
        call ConstructGridData      (GridDataID       = ObjBathymetry,              &
                                     HorizontalGridID = ObjHorizontalGrid,          &
                                     FileName         = BathymetryFile,             &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - AssimilationZones - ERR02'


        write(*,*)
        write(*,*)"Constructing mapping..."
        write(*,*)

        !Horizontal Map
        call ConstructHorizontalMap (HorizontalMapID  = ObjHorizontalMap,           &
                                     GridDataID       = ObjBathymetry,              &
                                     HorizontalGridID = ObjHorizontalGrid,          &
                                     ActualTime       = CurrentTime,                &
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - AssimilationZones - ERR03'

        write(*,*)
        write(*,*)"Compiling information..."
        write(*,*)

        call GetWaterPoints2D(HorizontalMapID  = ObjHorizontalMap,                  &
                              WaterPoints2D    = WaterPoints2D,                     &
                              STAT             = STAT_CALL)  

        call GetHorizontalGridSize(ObjHorizontalGrid, Size, WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - AssimilationZones - ERR04'

        call GetGridData(ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - AssimilationZones - ERR05'



        call GetCoordTypeList (GEOG         = GEOG,                             &
                               UTM          = UTM,                              &
                               MIL_PORT     = MIL_PORT,                         &
                               SIMPLE_GEOG  = SIMPLE_GEOG,                      &
                               GRID_COORD   = GRID_COORD,                       &
                               NLRD         = NLRD)

        !Gets Coordinates in use
        call GetGridCoordType(ObjHorizontalGrid, CoordType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - AssimilationZones - ERR06'

        if    (CoordType == SIMPLE_GEOG)then

            call GetGridLatitudeLongitude(ObjHorizontalGrid,                    &
                                          GridLatitudeConn  = YY_IE,            &
                                          GridLongitudeConn = XX_IE,            &
                                          STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - AssimilationZones - ERR07'

        elseif(CoordType == UTM        .or. CoordType == MIL_PORT .or. &
               CoordType == GRID_COORD .or. CoordType == NLRD)then

            call GetHorizontalGrid(ObjHorizontalGrid,                           &
                                   XX_IE = XX_IE,                               &
                                   YY_IE = YY_IE,                               &
                                   STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - AssimilationZones - ERR08'

        else

            write(*,*)'GEOG coordinate type cannot be used in digital terrain generation'
            stop 'ConstructDomain - AssimilationZones - ERR09'

        end if


        write(*,*)
        write(*,*)"Domain successfully constructed..."
        write(*,*)


    end subroutine ConstructDomain

    !--------------------------------------------------------------------------


    subroutine AllocateVariables

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        real                                        :: XSW, YSW
        real                                        :: XSE, YSE
        real                                        :: XNE, YNE
        real                                        :: XNW, YNW

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)"Allocating space in memory..."
        write(*,*)

        nullify (Zones2D_Z  ); allocate(Zones2D_Z   (Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        nullify (Zones2D_U  ); allocate(Zones2D_U   (Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        nullify (Zones2D_V  ); allocate(Zones2D_V   (Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        nullify (Zone_XX    ); allocate(Zone_XX     (Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        nullify (Zone_YY    ); allocate(Zone_YY     (Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        nullify (GridPoints ); allocate(GridPoints  (Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        nullify (GridFacesX ); allocate(GridFacesX  (Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        nullify (GridFacesY ); allocate(GridFacesY  (Size%ILB:Size%IUB, Size%JLB:Size%JUB))

        if (SurfaceForcing) then

            nullify (Zones3D_Z  ); allocate(Zones3D_Z   (Size%ILB:Size%IUB, Size%JLB:Size%JUB,1:LayersNumber))
            nullify (Zones3D_U  ); allocate(Zones3D_U   (Size%ILB:Size%IUB, Size%JLB:Size%JUB,1:LayersNumber))
            nullify (Zones3D_V  ); allocate(Zones3D_V   (Size%ILB:Size%IUB, Size%JLB:Size%JUB,1:LayersNumber))

        endif

        if(DefineNoAssimilationZone)then

            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB
                       
                XSW = XX_IE(i,      j    )
                YSW = YY_IE(i,      j    )
                XSE = XX_IE(i,      j + 1)
                YSE = YY_IE(i,      j + 1)
                XNE = XX_IE(i + 1,  j + 1)
                YNE = YY_IE(i + 1,  j + 1)
                XNW = XX_IE(i + 1,  j    )
                YNW = YY_IE(i + 1,  j    )

                GridPoints(i,j)%X = (XSW+XNW+XNE+XSE) / 4.
                GridPoints(i,j)%Y = (YSW+YNW+YNE+YSE) / 4.

            end do
            end do


            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB + 1
                       
                GridFacesX(i,j)%X = XX_IE(i,j)
                GridFacesX(i,j)%Y = YY_IE(i,j)

            end do
            end do

            do i = WorkSize%ILB, WorkSize%IUB + 1
            do j = WorkSize%JLB, WorkSize%JUB
                       
                GridFacesY(i,j)%X = XX_IE(i,j)
                GridFacesY(i,j)%Y = YY_IE(i,j)

            end do
            end do

        end if



    end subroutine AllocateVariables

    !--------------------------------------------------------------------------


    subroutine DefineZones
 
        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)"Defining assimilation zones..."
        write(*,*)

        if(TypeZ) call ComputeZonesZ
        if(TypeU) call ComputeZonesU
        if(TypeV) call ComputeZonesV

    end subroutine DefineZones


    !--------------------------------------------------------------------------

    subroutine DefineSurfaceForcing

        !Local-----------------------------------------------------------------
        integer                         :: i, j
 
        !Begin-----------------------------------------------------------------

        if(TypeZ) then

            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB

                Zones3D_Z(i,j,:) = Zones2D_Z(i,j)
                if (MinimumDominantes) then
                    Zones3D_Z(i,j,LayersNumber) = min(Zones2D_Z(i,j), SurfaceCoef)
                else
                    Zones3D_Z(i,j,LayersNumber) = max(Zones2D_Z(i,j), SurfaceCoef)
                endif
                    
            enddo
            enddo

        endif

        if(TypeU) then 

            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB + 1

                Zones3D_U(i,j,:) = Zones2D_U(i,j)
                if (MinimumDominantes) then
                    Zones3D_U(i,j,LayersNumber) = min(Zones2D_U(i,j), SurfaceCoef)
                else
                    Zones3D_U(i,j,LayersNumber) = max(Zones2D_U(i,j), SurfaceCoef)
                endif
                    
            enddo
            enddo

        endif

        if(TypeV) then 

            do i = WorkSize%ILB, WorkSize%IUB + 1
            do j = WorkSize%JLB, WorkSize%JUB

                Zones3D_V(i,j,:) = Zones2D_V(i,j)
                if (MinimumDominantes) then
                    Zones3D_V(i,j,LayersNumber) = min(Zones2D_V(i,j), SurfaceCoef)
                else
                    Zones3D_V(i,j,LayersNumber) = max(Zones2D_V(i,j), SurfaceCoef)
                endif
                    
            enddo
            enddo

        endif




    end subroutine DefineSurfaceForcing

    !--------------------------------------------------------------------------

    subroutine ComputeZonesZ

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        type(T_PointF), pointer                     :: OneGridPoint

        !Begin-----------------------------------------------------------------

        Zone_XX = DefaultCoef; Zone_YY = DefaultCoef; Zones2D_Z = DefaultCoef


        !Fill in XX direction
        do j = WorkSize%JLB, nCells
            if (West) Zone_XX(:,j)                    = AssimilationCoefs(j)
            if (East) Zone_XX(:,WorkSize%JUB - j + 1) = AssimilationCoefs(j)
        enddo


        !Fill in YY direction
        do i = WorkSize%ILB, nCells
            
            if (South) Zone_YY(i,:                   ) = AssimilationCoefs(i)
            if (North) Zone_YY(WorkSize%IUB - i + 1,:) = AssimilationCoefs(i)

        enddo
        
        do j = WorkSize%JLB, WorkSize%JUB
        do i = WorkSize%ILB, WorkSize%IUB

            if (WaterPoints2D(i,j) == WaterPoint) then

                OneGridPoint => GridPoints(i,j)

                if(IsVisible(NoAssimilationZone, OneGridPoint))then
                    
                    Zones2D_Z(i,j) = DefaultCoef

                else
                    if (MinimumDominantes) then                    
                        Zones2D_Z(i,j) = min (Zone_XX(i,j), Zone_YY(i,j))
                    else
                        Zones2D_Z(i,j) = max (Zone_XX(i,j), Zone_YY(i,j))
                    endif
                endif

            end if

        enddo
        enddo 

        nullify(OneGridPoint)

    end subroutine ComputeZonesZ


    !--------------------------------------------------------------------------


    subroutine ComputeZonesU

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        type(T_PointF), pointer                     :: OneGridFace

        !Begin-----------------------------------------------------------------

        Zone_XX = DefaultCoef; Zone_YY = DefaultCoef; Zones2D_U = DefaultCoef


        !Fill in XX direction
        do j = WorkSize%JLB, nCells
            if (West) Zone_XX(:,j)                        = AssimilationCoefs(j)
            if (East) Zone_XX(:,WorkSize%JUB + 1 - j + 1) = AssimilationCoefs(j)
        enddo


        !Fill in YY direction
        do i = WorkSize%ILB, nCells
            
            if (South) Zone_YY(i,:                   ) = AssimilationCoefs(i)
            if (North) Zone_YY(WorkSize%IUB - i + 1,:) = AssimilationCoefs(i)

        enddo
        
        do j = WorkSize%JLB, WorkSize%JUB + 1
        do i = WorkSize%ILB, WorkSize%IUB

            if (WaterPoints2D(i,j) == WaterPoint .or. WaterPoints2D(i,j-1) == WaterPoint) then

                OneGridFace => GridFacesX(i,j)

                if(IsVisible(NoAssimilationZone, OneGridFace))then
                    
                    Zones2D_U(i,j) = DefaultCoef

                else
                    if (MinimumDominantes) then                     
                        Zones2D_U(i,j) = min (Zone_XX(i,j), Zone_YY(i,j))
                    else    
                        Zones2D_U(i,j) = max (Zone_XX(i,j), Zone_YY(i,j))
                    endif
                endif

            end if

        enddo
        enddo 
        
        nullify   (OneGridFace)

    end subroutine ComputeZonesU


    !--------------------------------------------------------------------------


    subroutine ComputeZonesV

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        type(T_PointF), pointer                     :: OneGridFace

        !Begin-----------------------------------------------------------------

        Zone_XX = DefaultCoef; Zone_YY = DefaultCoef; Zones2D_V = DefaultCoef


        !Fill in XX direction
        do j = WorkSize%JLB, nCells
            if (West) Zone_XX(:,j                   ) = AssimilationCoefs(j)
            if (East) Zone_XX(:,WorkSize%JUB - j + 1) = AssimilationCoefs(j)
        enddo


        !Fill in YY direction
        do i = WorkSize%ILB, nCells
            
            if (South) Zone_YY(i,:                       ) = AssimilationCoefs(i)
            if (North) Zone_YY(WorkSize%IUB + 1 - i + 1,:) = AssimilationCoefs(i)

        enddo
        
        do j = WorkSize%JLB, WorkSize%JUB
        do i = WorkSize%ILB, WorkSize%IUB + 1
            
            if (WaterPoints2D(i,j) == WaterPoint.or. WaterPoints2D(i-1,j) == WaterPoint) then

                OneGridFace => GridFacesY(i,j)

                if(IsVisible(NoAssimilationZone, OneGridFace))then
                    
                    Zones2D_V(i,j) = DefaultCoef

                else
                    if (MinimumDominantes) then                                         
                        Zones2D_V(i,j) = min (Zone_XX(i,j), Zone_YY(i,j))
                    else
                        Zones2D_V(i,j) = max (Zone_XX(i,j), Zone_YY(i,j))
                    endif
                endif

            end if

        enddo
        enddo 

        nullify(OneGridFace)

    end subroutine ComputeZonesV

    !--------------------------------------------------------------------------

    subroutine WriteAssimilationZones

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: i, j
        character(len=PathLength)                   :: NewOutputFileName
        integer                                     :: NameLength, k

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)"Writing zones file..."
        write(*,*)

        NameLength          = len_trim(OutputFile)
        NewOutputFileName   = OutputFile(:NameLength-4)

        if(WriteAsGridData)then

            if (SurfaceForcing) then

                call WriteGridData (FileName            = trim(OutputFile),                 &
                                    COMENT1             = "Bathymetry filtered by Mohid",   &
                                    COMENT2             = "*",                              &
                                    HorizontalGridID    = ObjHorizontalGrid,                &
                                    FillValue           = -99.,                             &
                                    Overwrite           = .true.,                           &
                                    GridData3D_Real     = Zones3D_Z,                        &
                                    KLB                 = 1,                                &
                                    KUB                 = LayersNumber,                     &
                                    STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteAssimilationZones - AssimilationZones - ERR10'

            else

                call WriteGridData (FileName            = trim(OutputFile),                 &
                                    COMENT1             = "Bathymetry filtered by Mohid",   &
                                    COMENT2             = "*",                              &
                                    HorizontalGridID    = ObjHorizontalGrid,                &
                                    FillValue           = -99.,                             &
                                    Overwrite           = .true.,                           &
                                    GridData2D_Real     = Zones2D_Z,                        &
                                    STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteAssimilationZones - AssimilationZones - ERR20'

            endif

        end if


        if(WriteAsGridIndices .and. TypeZ)then


            call UnitsManager(OutputUnit, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteAssimilationZones - AssimilationZones - ERR30'

            open(Unit = OutputUnit, File = trim(NewOutputFileName)//"_Z.dat", Status = "unknown")


            write(OutputUnit,*)"<BeginGridData3D>"

            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB

                if (SurfaceForcing) then

                    do k=1, LayersNumber
                        if (.not. WriteDefaults .and. Zones3D_Z(i,j,k) == DefaultCoef) cycle
                        write(OutputUnit,20)i, j, k, Zones3D_Z(i,j,k)
                    enddo
                
                else
                    if (.not. WriteDefaults .and. Zones2D_Z(i,j) == DefaultCoef) cycle
                    write(OutputUnit,10)i, j,    Zones2D_Z(i,j)
                endif

            enddo
            enddo

            write(OutputUnit,*)"<EndGridData3D>"

            call UnitsManager(OutputUnit, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteAssimilationZones - AssimilationZones - ERR40'

        end if


        if(TypeU)then

            call UnitsManager(OutputUnit, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteAssimilationZones - AssimilationZones - ERR50'

            open(Unit = OutputUnit, File = trim(NewOutputFileName)//"_U.dat", Status = "unknown")

            call WriteDataLine(OutputUnit, "TYPE_ZUV", "U")

            write(OutputUnit,*)"<BeginGridData3D>"

            do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB + 1

                if (SurfaceForcing) then

                    do k=1, LayersNumber
                        if (.not. WriteDefaults .and. Zones3D_U(i,j,k) == DefaultCoef) cycle
                        write(OutputUnit,20)i, j, k, Zones3D_U(i,j,k)
                    enddo

                else

                    if (.not. WriteDefaults .and. Zones2D_U(i,j) == DefaultCoef) cycle
                    write(OutputUnit,10) i, j,   Zones2D_U(i,j)

                endif

            enddo
            enddo

            write(OutputUnit,*)"<EndGridData3D>"

            call UnitsManager(OutputUnit, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteAssimilationZones - AssimilationZones - ERR60'

        end if



        if(TypeV)then

            call UnitsManager(OutputUnit, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteAssimilationZones - AssimilationZones - ERR70'

            open(Unit = OutputUnit, File = trim(NewOutputFileName)//"_V.dat", Status = "unknown")

            call WriteDataLine(OutputUnit, "TYPE_ZUV", "V")

            write(OutputUnit,*)"<BeginGridData3D>"

            do i = WorkSize%ILB, WorkSize%IUB + 1
            do j = WorkSize%JLB, WorkSize%JUB

                if (SurfaceForcing) then

                    do k=1, LayersNumber
                        if (.not. WriteDefaults .and. Zones3D_V(i,j,k) == DefaultCoef) cycle
                        write(OutputUnit,20)i, j, k, Zones3D_V(i,j,k)
                    enddo

                else
                    if (.not. WriteDefaults .and. Zones2D_V(i,j) == DefaultCoef) cycle
                    write(OutputUnit,10)i, j, Zones2D_V(i,j)
                endif

            enddo
            enddo

            write(OutputUnit,*)"<EndGridData3D>"

            call UnitsManager(OutputUnit, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteAssimilationZones - AssimilationZones - ERR80'

        end if

        write(*,*)
        write(*,*)"Done..."
        write(*,*)

        10 format(i5, 1x, i5, 1x, e12.3)
        20 format(i5, 1x, i5, 1x, i5, 1x, e12.3)

    end subroutine WriteAssimilationZones

   
    !--------------------------------------------------------------------------


    subroutine DeallocateVariables

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)"Deallocating variables..."
        write(*,*)

        deallocate(Zones2D_Z        )   ; nullify (Zones2D_Z        )
        deallocate(Zones2D_U        )   ; nullify (Zones2D_U        )
        deallocate(Zones2D_V        )   ; nullify (Zones2D_V        )
        deallocate(Zone_XX          )   ; nullify (Zone_XX          )
        deallocate(Zone_YY          )   ; nullify (Zone_YY          )
        deallocate(AssimilationCoefs)   ; nullify (AssimilationCoefs)
        deallocate(GridPoints       )   ; nullify (GridPoints       )
        deallocate(GridFacesX       )   ; nullify (GridFacesX       )
        deallocate(GridFacesY       )   ; nullify (GridFacesY       )

    end subroutine DeallocateVariables


    !--------------------------------------------------------------------------

    subroutine KillAssimilationZones

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)"Finishing freeing memory..."
        write(*,*)

        call KillHorizontalMap  (ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillAssimilationZones - AssimilationZones - ERR01'

        call KillGridData       (ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillAssimilationZones - AssimilationZones - ERR02'

        call KillHorizontalGrid (ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillAssimilationZones - AssimilationZones - ERR03'

    end subroutine KillAssimilationZones


end program AssimilationZones
