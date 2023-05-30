!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Digital Terrain Creator
! PROGRAM       : Digital Terrain Creator
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2004
! REVISION      : Paulo Chambel / Luis Fernandes / Frank Braunschweig - v4.0
! DESCRIPTION   : Program to create Bathymetries or Basins
!
!------------------------------------------------------------------------------

!DataFile
!
!   <BeginXYZPointsFiles>
!                               : char              -           !Path to the xyz information files
!                                                               !Each line defines one file
!   <EndXYZPointsFiles>
!
!   <BeginLandAreaFiles>
!                               : char              -           !Path to the polygons files defining land
!                                                               !Each line defines one file
!   <EndLandAreaFiles>
!
!   BATIM_FILE                  : char              -           !Path to the bathymetry file to be created
!   GRID_FILE                   : char              -           !Path to the grid file that originates the bathymetry
!   BATIM_INI                   : char              -           !Path to the bathymetry file to be actualized
!   WINDOW_OUTPUT               : logical           -           !Window from Batim_Ini from each you want to extract a window
!   WINDOW_LIMITS               : 4*integer         -           ! ilb, jlb, iub, jub
!   MATCH_IN                    : 0/1              [1]          !It is true if the input grid data has the 
!                                                               !same grid of the new grid data
!   SMOOTH                      : 0/1               -           !Use radius method to fill unknown values
!       RADIUS                  : real             [0]          !Radius within which points are used to compute 
!                                                               !bathymetry value
!   FILL_METHOD                 : 1/2              [1]          !1 - Average value, 2 - Minimum value
!   NO_DATA_POINT               : real            [-99]         !Value to attribute to points with no value
!   LAND_POINT                  : real            [-99]         !Value to attribute to land points
!   INVALID_Z_VALUE             : real        [FillValueReal]   !Z value in points collection that will not be processed

!   DEPTH_REFERENTIAL           : real            [ 1]          ! If the option is 1 the bathymetry referential is use if
                                                                ! is -1 than is used the topographci referential. 

!   UNIT_FACTOR                 : real            [1.]          !This is the factor be each the z values of xyz files is multiply.

!   ADD_FACTOR                  : real            [0.]          !This is the factor be each the z values of xyz files is add.

!   INTERPOLATION               : char       [No Interpolation] !Interpolation method 

!       POINTS_FOR_INTERPOLATION: int              [1]          !Method to select base points for interpolation
!                                                               !1-AllPoints; 2-GridCellCenterPoints; 3-PointsInsideGridLimits
!       EXPAND_GRID_LIMITS      : 0/1              [0]          !Extend grid limits for points selection
!       GRID_LIMITS_PERCENTAGE  : real             [0]          !Percentage of grid limits to consider points
!                                                                   !for triangulation
!       CANONIC_SPACE           : 0/1              [0]          !This option is only valid for curvilinear grids
                                                                !This option makes triangulation in space where 
                                                                !The faces are always normal or paralel to the vectores defining 
                                                                !the space
!       
!       MAINTAIN_RIVER_BED      : 0/1                           !Check if the user wants to maintain the river bed larger depths.
!                                                               !Linear interpolation between cross sections where maximum depths
                                                                !of the river bed is near of opposite margins 
                                                                !tend to decrease maximum 
                                                                !of the inpolated sections and consequently increases 
                                                                !artificially the velocities 
!       RIVER_CHANGE_FACTOR     : 1.                            ! 1   - All section will be changed
                                                                ! 0.5 - Only depths greater than Hmax - 0.5*(Hmax-Hmin)
                                                                !       will be changed

!       DREDGE_DISTRIBUTION     : 1/2              [2]          !Distribuiton admitted for the bathymetry changes
!                                                               !1 - Uniform, 2 - Linear

!       RIVER_FLOW_DIRECTION    : i/j              [j]          !The user specified allong each grid direction (i or j)     

!       SECTIONS_LOCATION       : 1/2              [1]          !Methodology used to defined section location in the grid
                                                                !1 - cell with the maximum depth;
                                                                !2 - This option considers a section starting in the first cell 
                                                                !    with data until the last cell with data 

!   FILTER_DEPTHS               : 0/1              [0]          !Filter bathymetry
!       FILTER_RADIUS_I         : real             [1]          !Filter radius along i
!       FILTER_RADIUS_J         : real             [1]          !Filter radius along j

!   FILTER_SPIKES               : 0/1              [0]          !Remove spikes from the bathymetry
!       SPIKES_FACTOR           : real             [10]         !Parameter where the user controles the intensity 
!                                                               !of spikes to be removed. 


!                                                               !No Interpolation; Triangulation; Spline
!       FILL_OUTSIDE_POINTS     : 0/1              [0]          !Use the closest points for points which lie
!                                                               !outside the polygon
!       TOLERANCE               : real             [0]          !Tolerance distance between nodes (coordinates units)
!       WRITE_TRIANGLES         : 0/1              [0]          !Writes a polygon file with the created triangles
!       TRIANGLES_FILE          : char              -           !Path to the triangles file to create   
!   ASSUME_LAND_DEPTH           : 0/1              [0]          !Assumed a specified depth to land boundary points
!       LAND_DEPTH              : real              -           !Depth to be assumed in land boundary points
!   OVERLAPPING                 : 0/1              [0]          !Overlap information based on grid data and polygons files
!   OVERLAPPING_NUMBER          : int               -           !Number of grid data files used to overlap information

!   CONVERT_2_XBEACH            : 0/1              [0]          !Write also the bathymetry in the XBEACH format

!
!   <BeginGridDataInfo>
!       LEVEL                   : int               -           !Order of importance of the grid data file
!       GRIDATA_FILE            : char              -           !Grid data file to be used in overlapping
!       H_MIN                   : real            [-1e9]        !Minimum value to be used in overlapping
!       H_MAX                   : real            [ 1e9]        !Maximum value to be used in overlapping
!       PERCENTAGE              : real             [1]          !Percentage (0-1) of importance given to this Grid Data
!       <<BeginAreaFiles>>  
!                               : char              -           !Path to the polygons files defining area of the
!                                                               !grid data file to be used. Each line defines one file
!       <<EndAreaFiles>>
!   <EndGridDataInfo>
!
!
!!!Create river section dredging!!!
!
!   SECTION_COMPUTE            : 0/1               [0]           !1 - Compute sections from given xyz points 
!                                                                   (pairs of section begin and end)
!   SECTION_MIDDLE_POINT_METHOD : int              [1]           !1 - use middle height; 2 - use Z in original xyz collection
!   SECTION_HEIGHT             : real              [0.]          !Section middle height (from first point Z)  
!                                                                   - read if SECTION_MIDDLE_POINT_METHOD : 1
!   SECTION_INTERPOL_TYPE      : int               [1]           !1 - linear; 2 - power; 3 - root
!   SECTION_INTERPOL_ORDER     : real              [1]           !power or root order (2 - power of 2, squared root)
!   SECTION_POINTS_TO_ADD      : int               [7]           !points to add to each section (beside section end and start)
!   SECTION_XYZ_FILE           : char               -            !output file with original points (section margins) 
!                                                                  and created section points
!   CREATE_MARGIN_POINTS       : 0/1               [0]           !1 - add more points to margins to make the trinagulation smoother
!                                                                 and only need the user to define main river changes 0 - don't add
!         MARGIN_POINTS_TO_ADD : int               [7]           !how many points do add between defined points (read if 
!                                                                  CREATE_MARGIN_POINTS : 1)
!
!   !areas to select where to triangulate (if active)
!   <BeginAreaFiles>
!   <EndAreaFiles>

program DigitalTerrainCreator
    
    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime
    use ModuleDrawing
    use ModuleFunctions
    use ModuleHorizontalGrid
    use ModuleGridData
    use ModuleTriangulation

    implicit none

    !Parameters----------------------------------------------------------------
    integer, parameter                              :: NoInterpolation              = 0
    integer, parameter                              :: Triangulation                = 1
    integer, parameter                              :: Spline                       = 2
    
    integer, parameter                              :: AllPoints                    = 1
    integer, parameter                              :: GridCellCenterPoints         = 2
    integer, parameter                              :: PointsInsideGridLimits       = 3
    integer, parameter                              :: OnlyXYZPoints                = 4

    !River axes
    integer, parameter                              :: AlongLine                    = 1
    integer, parameter                              :: AlongColumn                  = 2

    integer, parameter                              :: Uniform                      = 1
    integer, parameter                              :: Linear                       = 2


    integer, parameter                              :: MaximumDepth                 = 1
    integer, parameter                              :: DefinedDepths                = 2


    character(LEN=StringLength), parameter          :: Char_GridDataType            = 'grid data'
    character(LEN=StringLength), parameter          :: Char_CellsType               = 'cells'

    integer, parameter                              :: GridDataType                 = 1
    integer, parameter                              :: CellsType                    = 2
    
    integer, parameter                              :: Average_                     = 1
    integer, parameter                              :: Minimum_                     = 2 
    integer, parameter                              :: Maximum_                     = 3
    
    !Z interpolation
    integer, parameter                              :: LinearInterpol_              = 1
    integer, parameter                              :: PowerInterpol_               = 2
    integer, parameter                              :: RootInterpol_                = 3
    
    !MiddlePoint Definition
    integer, parameter                              :: DefinedHeight_               = 1
    integer, parameter                              :: DefinedLevel_                = 2

    !Globals-------------------------------------------------------------------
   
    type     T_ExternalVar
        real, pointer, dimension(:,:)               :: XX_IE, YY_IE 
        real, pointer, dimension(:  )               :: XX, YY
        integer, pointer, dimension(:,:)            :: DefineCellsMap
        type (T_Size2D)                             :: WorkSize
        type (T_Size2D)                             :: Size
        logical                                     :: GridDistortion
        real                                        :: OriginX, OriginY, Rotation
    end type T_ExternalVar
    
    type     T_Triang
        logical                                     :: FillOutsidePoints            = .false.
        real                                        :: Tolerance                    = null_real
        logical                                     :: OutputTriangles              = .false.
        character(len=PathLength)                   :: File
    end type T_Triang


    type     T_Cells
        integer                                     :: i                            = null_int
        integer                                     :: j                            = null_int
        real                                        :: Depth                        = null_real
    end type T_Cells

    type     T_DataInfo
        type(T_Polygon),pointer                     :: Areas
        real                                        :: Hmin,Hmax                    = null_real
        real                                        :: Percentage                   = 1.
        integer                                     :: ObjGridData                  = 0
        type (T_Cells), dimension(:), pointer       :: Cells
        integer                                     :: CellsNumber                  = null_int
        integer                                     :: InfoType                     = null_int
    end type T_DataInfo

        
    type T_Overlapping
        logical                                     :: ON                           = .false.
        integer                                     :: Number                       = 0
        type (T_DataInfo), dimension(:), pointer    :: DataInfo
    end type T_Overlapping

    type T_River
        logical                                     :: CanonicSpace, MaintainBed
        integer                                     :: MainAxe, SectionsNumber
        integer                                     :: Dredge
        real                                        :: FactorLimit
        integer                                     :: LocationMethod
        integer, pointer, dimension(:,:)            :: SectionsLocation
        real   , pointer, dimension(:  )            :: BedDepth
    end type T_River
    
    type T_BorderLimits
        !Cell limits
        integer                                     :: ILB    = null_int
        integer                                     :: IUB    = null_int        
        integer                                     :: JLB    = null_int
        integer                                     :: JUB    = null_int
    end type T_BorderLimits
    
    type T_DD
        logical                                         :: ON   = .false.
        !Number of domains in the J/X direction
        integer                                         :: nJX  = FillValueInt  
        !Number of domains in the I/Y direction
        integer                                         :: nIY  = FillValueInt  
        !Number of domains
        integer                                         :: nD   = FillValueInt
        !Borders of the each domain
        type (T_Border),       dimension(:), pointer    :: Border       => null()
        !Borders Limits of the each domain
        type (T_BorderLimits), dimension(:), pointer    :: BorderLimits => null()        
        
    end type T_DD        
    
    type T_Global
        integer                                     :: ObjGrid                      = 0
        integer                                     :: ObjGridIn                    = 0
        integer                                     :: ObjGridDataIn                = 0
        integer                                     :: ObjTriangulation             = 0
        integer                                     :: ObjEnterData                 = 0
        real,           dimension(:,:),    pointer  :: Depth
        type(T_Polygon),                   pointer  :: LandArea, Rect
        type(T_XYZPoints),                 pointer  :: XYZPoints, XYZPointsOriginal, CurrentXYZPoints
        type(T_PointF), dimension(:,:),    pointer  :: GridPoint 
        type(T_PointF),                    pointer  :: AuxPoint
        real,           dimension(:,:),    pointer  :: PointsWithData, PointsNoData
        integer                                     :: TotalXYZ                     = 0
        integer                                     :: TotalWithData                = 0
        integer                                     :: TotalNoData                  = 0
        character(len=line_length)                  :: BatimFilePathIn
        character(len=line_length)                  :: BatimFilePathOut
        character(len=line_length)                  :: BaseAndComputedInformation   = "BaseAndComputedInformation.xyz"
        character(len=line_length)                  :: PointsWithNoInformation      = "PointsWithNoInformation.xyz"
        integer                                     :: PointsForInterpolation
        logical                                     :: UseSmooth, BatimInON, BatimInMatch
        logical                                     :: ExpandGridLimits             = .false.
        real                                        :: ExpandGridLimitsFraction     = null_real
        real                                        :: Radius                       = null_real
        real                                        :: LandDepth                    = null_real
        real                                        :: UnitFactor, AddFactor
        logical                                     :: AssumeLandDepth
        type (T_Time)                               :: InitialSystemTime, FinalSystemTime
        integer, dimension(8)                       :: F95Time
        integer                                     :: InterpolType
        real                                        :: NoDataPoint, LandPoint, DepthReferential
        real                                        :: InvalidZValue
        integer                                     :: FilterRadiusI, FilterRadiusJ
        logical                                     :: Filter, FilterSpikes
        real                                        :: SpikesFactor

        logical                                     :: WindowOutPut
        integer, dimension(4)                       :: WindowLimits
        
        integer                                     :: FillMethod
        
        logical                                     :: ComputeSections
        integer                                     :: SectionPointsToAdd
        integer                                     :: SectionInterpolType
        integer                                     :: SectionMiddlePointMethod
        real                                        :: SectionInterpolOrder
        real                                        :: SectionHeight
        character(len=line_length)                  :: SectionXYZFile
        type(T_Polygon),pointer                     :: SectionAreas
        logical                                     :: CreateMarginPoints
        integer                                     :: MarginPointsToAdd
        logical                                     :: ExportXYZ
        logical                                     :: Convert2Xbeach
        logical                                     :: OutputSupportFiles
        
        real                                        :: TriangScale = 1000.

        type(T_ExternalVar)                         :: ExtVar
        type(T_Triang     )                         :: Triang
        type(T_Overlapping)                         :: Overlapping
        type(T_River      )                         :: River
        type(T_Limits     )                         :: GridLimits
        type(T_DD         )                         :: DD
        type(T_Border     ), pointer                :: AllBorder
        
    end type T_Global                                 
                                                     
    type (T_Global)                                 :: Me

    !Begin---------------------------------------------------------------------
    call OpenProject
    call RunProject  
    call CloseProject

    contains
    
    !--------------------------------------------------------------------------
    
    subroutine OpenProject
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        character(PathLength), parameter            :: ProjectFilePath  = "CreateBathymetry.dat"

        !Begin-----------------------------------------------------------------


        call StartupMohid ("Digital Terrain Creator")

        !Gets the actual time
        call date_and_time(Values = Me%F95Time)
        call SetDate      (Me%InitialSystemTime, float(Me%F95Time(1)), float(Me%F95Time(2)), &
                                                 float(Me%F95Time(3)), float(Me%F95Time(5)), &
                                                 float(Me%F95Time(6)), float(Me%F95Time(7))+ &
                                                 Me%F95Time(8)/1000.)


        call ConstructEnterData(Me%ObjEnterData, ProjectFilePath, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - DigitalTerrainCreator - ERR01'

        call ReadGridFilesNames
        
        call SetGridLimits

        call ConstructGlobalOptions

        call ConstructLandZones

        call ConstructBaseInformation


        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - DigitalTerrainCreator - ERR190'
        

    end subroutine OpenProject

    
    !--------------------------------------------------------------------------

    
    subroutine ConstructBaseInformation
    
        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL
        character(len=line_length)          :: XYZFilePath
        character(len=line_length)          :: FullBufferLine
        integer                             :: ClientNumber, StartLine, EndLine
        integer                             :: CurrentLineNumber
        real                                :: Factor
        logical                             :: BlockFound
        
        !Begin-----------------------------------------------------------------

        call ExtractBlockFromBuffer(Me%ObjEnterData,                                     &
                                    ClientNumber,                                        &
                                    '<BeginXYZPointsFiles>',                             &
                                    '<EndXYZPointsFiles>',                               &
                                    BlockFound,                                          &
                                    STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructBaseInformation - DigitalTerrainCreator - ERR10'

!        if (Me%BatimInON .and. Me%BatimInMatch .and. BlockFound) stop 'ConstructBaseInformation - DigitalTerrainCreator - ERR15'
        
        call GetBlockSize(Me%ObjEnterData,                                               &
                          ClientNumber,                                                  &
                          StartLine,                                                     &
                          EndLine,                                                       &
                          FromBlock,                                                     &
                          STAT = STAT_CALL)
        
        do CurrentLineNumber = StartLine + 1 , EndLine - 1

            call GetFullBufferLine(Me%ObjEnterData,                                      &
                                   CurrentLineNumber,                                    &
                                   FullBufferLine,                                       &
                                   STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructBaseInformation - DigitalTerrainCreator - ERR20'
            
            !Get LandArea files paths
            call GetData(XYZFilePath,                                                    &
                         Me%ObjEnterData,                                                &
                         flag,                                                           &
                         SearchType   = FromBlock_,                                      &
                         Buffer_Line  = CurrentLineNumber,                               &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructBaseInformation - DigitalTerrainCreator - ERR30'
            
            write(*,*)"Constructing XYZ..."
            
            if (mod(CurrentLineNumber, 10000) == 0) then
                write(*,*)'Reading line :', CurrentLineNumber, ' of ', EndLine
            endif

            Factor = Me%DepthReferential * Me%UnitFactor

            !Construct XYZ points collection
            call New(Me%XYZPoints, trim(XYZFilePath), DepthReferential = Factor, AddFactor = Me%AddFactor)
            
!            !copy points so that Z is saved (will be middle point level)
!            if (Me%ComputeSections .and. Me%SectionMiddlePointMethod == DefinedLevel_) then
!                call New(Me%XYZPointsOriginal, trim(XYZFilePath), DepthReferential = Factor, AddFactor = Me%AddFactor)
!            endif

        end do

        Me%TotalXYZ = 0

        Me%CurrentXYZPoints => Me%XYZPoints

        do while(associated(Me%CurrentXYZPoints))

            Me%TotalXYZ = Me%TotalXYZ + Me%CurrentXYZPoints%Count

            Me%CurrentXYZPoints => Me%CurrentXYZPoints%Next

        enddo 

        nullify(Me%CurrentXYZPoints)

                 
        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructBaseInformation - DigitalTerrainCreator - ERR40'


    
    end subroutine ConstructBaseInformation


    !--------------------------------------------------------------------------

    
    subroutine ConstructLandZones

        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL
        character(len=line_length)          :: LandAreaFilePath
        character(len=line_length)          :: FullBufferLine
        integer                             :: ClientNumber, StartLine, EndLine
        integer                             :: CurrentLineNumber
        logical                             :: BlockFound
        type (T_Polygon), pointer           :: AllPolygons, PolygonAux, PolygonAux2
        !Begin-----------------------------------------------------------------
        
        call ExtractBlockFromBuffer(Me%ObjEnterData,                                    &
                                    ClientNumber,                                       &
                                    '<BeginLandAreaFiles>',                             &
                                    '<EndLandAreaFiles>',                               &
                                    BlockFound,                                         &
                                    STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructLandZones - DigitalTerrainCreator - ERR10'

        if (Me%BatimInON .and. Me%BatimInMatch .and. BlockFound) stop 'ConstructLandZones - DigitalTerrainCreator - ERR15'    
    
        call GetBlockSize(Me%ObjEnterData,                                              &
                          ClientNumber,                                                 &
                          StartLine,                                                    &
                          EndLine,                                                      &
                          FromBlock_,                                                   &
                          STAT = STAT_CALL)

        do CurrentLineNumber = StartLine + 1 , EndLine - 1

            call GetFullBufferLine(Me%ObjEnterData,                                     &
                                   CurrentLineNumber,                                   &
                                   FullBufferLine,                                      &
                                   STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructLandZones - DigitalTerrainCreator - ERR20'
        
            !Get LandArea files paths
            call GetData(LandAreaFilePath,                                              &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType  = FromBlock_,                                      &
                         Buffer_Line = CurrentLineNumber,                               &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructLandZones - DigitalTerrainCreator - ERR30'

            write(*,*)"Constructing Land Area..."

            !Construct LandArea collection
            call New(AllPolygons,  trim(LandAreaFilePath))

        end do
        
        PolygonAux => AllPolygons
        
        nullify(Me%LandArea)
       
        do while (associated(PolygonAux))
            
            allocate(PolygonAux2)
            
            PolygonAux2 = PolygonAux
            
            nullify(PolygonAux2%Next)
            
            call SetLimits(PolygonAux2)            
            
            if (Polygon_Intersect_Grid(PolygonAux2)) then
                call Add(Me%LandArea,  PolygonAux2)
            endif
            
            PolygonAux => PolygonAux%Next
            
        enddo
    
        
        nullify(AllPolygons, PolygonAux, PolygonAux2)
        
        if (Me%OutputSupportFiles) then
            call writeItem (Me%LandArea, 'Final_Land.xy')
        endif
            

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructLandZones - DigitalTerrainCreator - ERR40'

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructLandZones - DigitalTerrainCreator - ERR50'

    end subroutine ConstructLandZones

    !--------------------------------------------------------------------------
    
    logical function Polygon_Intersect_Grid (PolygonX)
    
        !Arguments-------------------------------------------------------------
        type (T_Polygon), pointer           :: PolygonX
    
        !Local-----------------------------------------------------------------
        real                                :: XminA, XmaxA, YminA, YmaxA
        real                                :: XminB, XmaxB, YminB, YmaxB
        !Begin-----------------------------------------------------------------    
        
        
        XminA = Me%GridLimits%Left  
        XmaxA = Me%GridLimits%Right 
        YminA = Me%GridLimits%Bottom
        YmaxA = Me%GridLimits%Top   
        
        call GetPolyLimits(PolygonX, XminB, XmaxB, YminB, YmaxB)
        
        if (.not. XminB > XmaxA .and. .not. XmaxB < XminA .and. &
            .not. YminB > YmaxA .and. .not. YmaxB < YminA) then
            Polygon_Intersect_Grid = .true.
        else
            Polygon_Intersect_Grid = .false.                    
        endif          
    
    end function Polygon_Intersect_Grid
    
    !--------------------------------------------------------------------------    

    subroutine ReadGridFilesNames
        
        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL
        character(len=line_length)          :: GridFilePath

        !Begin-----------------------------------------------------------------


        !Get grid file path
        call GetData(Me%BatimFilePathOut,                                               &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='BATIM_FILE',                                        &
                     ClientModule ='DigitalTerrainCreator',                             &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - DigitalTerrainCreator - ERR10'
        
        
        !Get grid file path
        call GetData(GridFilePath,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='GRID_FILE',                                         &
                     ClientModule ='DigitalTerrainCreator',                             &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - DigitalTerrainCreator - ERR20'
        
        !Construct grid
        call ConstructHorizontalGrid(Me%ObjGrid, trim(GridFilePath), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - DigitalTerrainCreator - ERR30'

        !Initial bathymetry to be updated
        call GetData(Me%BatimFilePathIn,                                                &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     ClientModule ='DigitalTerrainCreator',                             &
                     keyword      ='BATIM_INI',                                         &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - DigitalTerrainCreator - ERR40'

        if (flag /= 0) then
            Me%BatimInON = .true.
        else
            Me%BatimInON = .false.
        endif


    end subroutine ReadGridFilesNames


    !--------------------------------------------------------------------------

    
    subroutine ConstructGlobalOptions
        
        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL
        character(len=line_length)          :: AuxChar

        !Begin-----------------------------------------------------------------


        call GetData(Me%PointsForInterpolation,                                         &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='POINTS_FOR_INTERPOLATION',                          &
                     ClientModule ='DigitalTerrainCreator',                             &
                     Default      = AllPoints,                                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR10'

        call GetData(Me%River%CanonicSpace,                                             &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='CANONIC_SPACE',                                     &
                     ClientModule ='DigitalTerrainCreator',                             &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR20'

        if (Me%River%CanonicSpace .and. Me%PointsForInterpolation /= GridCellCenterPoints) then
            stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR30'
        endif

        
        !Could not understand how to use this with river sections to dredge so created a new way
        !to define section by creating it from section margin (XYZ file) and middle depth
        !David 12-2013
i1:     if (Me%River%CanonicSpace) then

            call GetData(Me%River%MaintainBed,                                          &
                         Me%ObjEnterData, flag,                                         &
                         SearchType   = FromFile_,                                      &
                         keyword      ='MAINTAIN_RIVER_BED',                            &
                         ClientModule ='DigitalTerrainCreator',                         &
                         Default      = .false.,                                        &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR40'

            call GetData(Me%River%FactorLimit,                                          &
                         Me%ObjEnterData, flag,                                         &
                         SearchType   = FromFile_,                                      &
                         keyword      ='RIVER_CHANGE_FACTOR',                           &
                         ClientModule ='DigitalTerrainCreator',                         &
                         Default      = 1.,                                             &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR45'


            call GetData(AuxChar,                                                       &
                         Me%ObjEnterData, flag,                                         &
                         SearchType   = FromFile_,                                      &
                         keyword      ='RIVER_FLOW_DIRECTION',                          &
                         ClientModule ='DigitalTerrainCreator',                         &
                         Default      = 'j',                                            &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR50'

i2:         if      (trim(AuxChar) == 'j') then

                Me%River%MainAxe = AlongLine

            else if (trim(AuxChar) == 'i') then i2
                
                Me%River%MainAxe = AlongColumn

            else  i2

                stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR60'

            endif i2

            call GetData(Me%River%Dredge,                                               &
                         Me%ObjEnterData, flag,                                         &
                         SearchType   = FromFile_,                                      &
                         keyword      ='DREDGE_DISTRIBUTION',                           &
                         ClientModule ='DigitalTerrainCreator',                         &
                         Default      = Linear,                                         &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR65'

            call GetData(Me%River%LocationMethod,                                       &
                         Me%ObjEnterData, flag,                                         &
                         SearchType   = FromFile_,                                      &
                         keyword      ='SECTIONS_LOCATION',                             &
                         ClientModule ='DigitalTerrainCreator',                         &
                         Default      = MaximumDepth,                                   &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66'

        endif i1
        
        !The latter canonic space does not seem usefull in creating real river sections so 
        !created a new method to compute river sections from given points. The given points should be
        !even, created in pairs, being the first of each pair the section begin and the second the section end
        !As these points are constructed in GUI the z is undefined and will be filled with bathymetry in
        !The sections are computed linearly, root or power from section begin given middle point height
        call GetData(Me%ComputeSections,                                                &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='SECTION_COMPUTE',                                   &
                     ClientModule ='DigitalTerrainCreator',                             &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66a'
        
        if (Me%ComputeSections) then
        
            if (.not. Me%BatimInON) then
                write (*,*)
                write (*,*) 'Define the keyword BATIM_INI, defining the path to the bathymetry file'
                write (*,*) 'to be able to compute sections margin elevations'
                stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66b'                
            endif
            
            !interpolation type to create the section points z's (linear, root or power producing
            !a V, U and "hart" section shapes respectively
            call GetData(Me%SectionInterpolType,                                            &
                         Me%ObjEnterData, flag,                                             &
                         SearchType   = FromFile_,                                          &
                         keyword      ='SECTION_INTERPOL_TYPE',                             &
                         ClientModule ='DigitalTerrainCreator',                             &
                         Default      = LinearInterpol_,                                    &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66c'
            
            if (Me%SectionInterpolType .ne. LinearInterpol_         &
                .and. Me%SectionInterpolType .ne. PowerInterpol_  &
                .and. Me%SectionInterpolType .ne. RootInterpol_) then
                write (*,*)
                write (*,*) 'Define the keyword SECTION_INTERPOL_TYPE with value 1'
                write (*,*) '(linear) or 2 (power) or 3 (root)'
                stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66c2'              
            endif 
            
            !for power and root define the order (2 - power to 2 or squared root)
            !by default 1 (same as linear)
            call GetData(Me%SectionInterpolOrder,                                           &
                         Me%ObjEnterData, flag,                                             &
                         SearchType   = FromFile_,                                          &
                         keyword      ='SECTION_INTERPOL_ORDER',                            &
                         ClientModule ='DigitalTerrainCreator',                             &
                         Default      = 1.,                                                 &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66c3'
            
            !if the user will define a section depth for all (method 1) or a middle point level (in first point Z) - method 2
            call GetData(Me%SectionMiddlePointMethod,                                       &
                         Me%ObjEnterData, flag,                                             &
                         SearchType   = FromFile_,                                          &
                         keyword      ='SECTION_MIDDLE_POINT_METHOD',                       &
                         ClientModule ='DigitalTerrainCreator',                             &
                         Default      = DefinedHeight_,                                     &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66d'
            
            if (Me%SectionMiddlePointMethod .ne. DefinedHeight_ .and. Me%SectionMiddlePointMethod .ne. DefinedLevel_) then
                write (*,*) 
                write (*,*) 'SECTION_MIDDLE_POINT_METHOD can only be 1 (defined height) or 2 (defined level)'
                stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66d2'            
            endif
            
            if (Me%SectionMiddlePointMethod == DefinedHeight_) then
                !height for middle point of section and it is measured in relation to first point Z of the pair
                call GetData(Me%SectionHeight,                                                  &
                             Me%ObjEnterData, flag,                                             &
                             SearchType   = FromFile_,                                          &
                             keyword      ='SECTION_HEIGHT',                                    &
                             ClientModule ='DigitalTerrainCreator',                             &
                             Default      = 0.,                                                 &
                             STAT         = STAT_CALL)        
                if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66d'
                
                if (Me%SectionHeight < 0.0) then
                    write (*,*) 
                    write (*,*) 'SECTION_HEIGHT negative. This value need to be positive.'
                    stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66e'
                endif
            endif
            
            !How many points to add in the section beside the margin ones (start and end)
            call GetData(Me%SectionPointsToAdd,                                             &
                         Me%ObjEnterData, flag,                                             &
                         SearchType   = FromFile_,                                          &
                         keyword      ='SECTION_POINTS_TO_ADD',                             &
                         ClientModule ='DigitalTerrainCreator',                             &
                         Default      = 7,                                                  &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66f'
            
            !Output file with xyz from original file (but with z filled) and also the added points
            call GetData(Me%SectionXYZFile,                                                 &
                         Me%ObjEnterData, flag,                                             &
                         SearchType   = FromFile_,                                          &
                         keyword      ='SECTION_XYZ_FILE',                                  &
                         ClientModule ='DigitalTerrainCreator',                             &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66g'    
            
            if (flag == 0) then
                write (*,*)
                write (*,*) 'Define the keyword SECTION_XYZ_FILE, defining the path to the file'
                write (*,*) 'to place the original xyz points and the created by the section method'
                stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66h'
            endif

            !Create more margin points? to be able to give only the main river changes and the
            !program creates the remaining so that the the interpolation is smooth
            call GetData(Me%CreateMarginPoints,                                             &
                         Me%ObjEnterData, flag,                                             &
                         SearchType   = FromFile_,                                          &
                         keyword      ='CREATE_MARGIN_POINTS',                              &
                         ClientModule ='DigitalTerrainCreator',                             &
                         Default      = .false.,                                            &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66i'
            
            if (Me%CreateMarginPoints) then
                !How many points to add in the margin beside the original ones (start and end)
                call GetData(Me%MarginPointsToAdd,                                              &
                             Me%ObjEnterData, flag,                                             &
                             SearchType   = FromFile_,                                          &
                             keyword      ='MARGIN_POINTS_TO_ADD',                              &
                             ClientModule ='DigitalTerrainCreator',                             &
                             Default      = 7,                                                  &
                             STAT         = STAT_CALL)        
                if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR66j'
            endif
            
            call ConstructSelectedArea
                   
        endif

        call GetData(Me%Filter,                                                     &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile_,                                      &
                     keyword      ='FILTER_DEPTHS',                                 &
                     ClientModule ='DigitalTerrainCreator',                         &
                     Default      = .false.,                                        &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR67'

        call GetData(Me%FilterRadiusJ,                                              &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile_,                                      &
                     keyword      ='FILTER_RADIUS_J',                               &
                     ClientModule ='DigitalTerrainCreator',                         &
                     Default      = 1,                                              &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR68'

        call GetData(Me%FilterSpikes,                                               &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile_,                                      &
                     keyword      ='FILTER_SPIKES',                                 &
                     ClientModule ='DigitalTerrainCreator',                         &
                     Default      = .false.,                                        &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR72'

        call GetData(Me%SpikesFactor,                                               &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile_,                                      &
                     keyword      ='SPIKES_FACTOR',                                 &
                     ClientModule ='DigitalTerrainCreator',                         &
                     Default      = 10.,                                            &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR74'



        call GetData(Me%FilterRadiusI,                                              &
                     Me%ObjEnterData, flag,                                         &
                     SearchType   = FromFile_,                                      &
                     keyword      ='FILTER_RADIUS_I',                               &
                     ClientModule ='DigitalTerrainCreator',                         &
                     Default      = 1,                                              &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR76'

        call GetData(Me%ExpandGridLimits,                                               &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     ClientModule ='DigitalTerrainCreator',                             &
                     keyword      ='EXPAND_GRID_LIMITS',                                &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR78'

        if(Me%ExpandGridLimits)then
            
            call GetData(Me%ExpandGridLimitsFraction,                                   &
                         Me%ObjEnterData, flag,                                         &
                         SearchType   = FromFile_,                                      &
                         keyword      ='GRID_LIMITS_PERCENTAGE',                        &
                         ClientModule ='DigitalTerrainCreator',                         &
                         Default      = 0.,                                             &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR80'
            
            if(Me%ExpandGridLimitsFraction < 0.)then
                stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR90'
            end if

        end if

        call GetData(Me%UseSmooth,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='SMOOTH',                                            &
                     Default      = .false.,                                            &
                     ClientModule ='DigitalTerrainCreator',                             &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR100'

        call GetData(Me%NoDataPoint,                                                    &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='NO_DATA_POINT',                                     &
                     ClientModule ='DigitalTerrainCreator',                             &
                     Default      = -999.,                                              &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR110'


        call GetData(Me%InvalidZValue,                                                  &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='INVALID_Z_VALUE',                                   &
                     ClientModule ='DigitalTerrainCreator',                             &
                     Default      = FillValueReal,                                      &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR115'
        
        call GetData(Me%LandPoint,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='LAND_POINT',                                        &
                     ClientModule ='DigitalTerrainCreator',                             &
                     Default      = -99.,                                               &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR120'

        call GetData(Me%DepthReferential,                                               &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='DEPTH_REFERENTIAL',                                 &
                     ClientModule ='DigitalTerrainCreator',                             &
                     Default      = 1.,                                                 &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR130'


        if (Me%DepthReferential /= 1. .and. Me%DepthReferential /= -1.)                 &
            stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR140'


        call GetData(Me%UnitFactor,                                                     &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='UNIT_FACTOR',                                       &
                     ClientModule ='DigitalTerrainCreator',                             &
                     Default      = 1.,                                                 &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR150'

        call GetData(Me%AddFactor,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='ADD_FACTOR',                                        &
                     ClientModule ='DigitalTerrainCreator',                             &
                     Default      = 0.,                                                 &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR160'



        !Get type of interpolation
        call GetData(AuxChar,                                                           &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='INTERPOLATION',                                     &
                     Default      = "No Interpolation",                                 &
                     ClientModule ='DigitalTerrainCreator',                             &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR170'

        select case (trim(AuxChar))
            case ("No Interpolation")
                Me%InterpolType = NoInterpolation
            case ("Triangulation")
                Me%InterpolType = Triangulation
!           case ("Spline")
!               Me%InterpolType = Spline
            case default 
                stop 'Interpolation method not valid'
        end select

        if (Me%InterpolType == Triangulation) call ReadTriangulationOptions


        if(Me%UseSmooth)then
            
            call GetData(Me%Radius,                                                     &
                         Me%ObjEnterData, flag,                                         &
                         SearchType   = FromFile_,                                      &
                         keyword      ='RADIUS',                                        &
                        ClientModule  ='DigitalTerrainCreator',                         &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR180'


        end if


        !Gets if wants to associate depth with land points
        call GetData(Me%AssumeLandDepth,                                                &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='ASSUME_LAND_DEPTH',                                 &
                     ClientModule ='DigitalTerrainCreator',                             &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR190'

        if(Me%AssumeLandDepth)then

            call GetData(Me%LandDepth,                                                  &
                         Me%ObjEnterData, flag,                                         &
                         SearchType   = FromFile_,                                      &
                         keyword      ='LAND_DEPTH',                                    &
                         ClientModule ='DigitalTerrainCreator',                         &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - DigitalTerrainCreator - ERR200'

            if (flag == 0) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR210'

        end if

        !Get grid file path
        call GetData(Me%Overlapping%ON,                                                  &
                     Me%ObjEnterData, flag,                                              &
                     SearchType   = FromFile_,                                           &
                     Default      = .false.,                                             &
                     keyword      ='OVERLAPPING',                                        &
                     ClientModule ='DigitalTerrainCreator',                              &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR220'

        if (Me%Overlapping%ON) call ConstructOverlapping

        !Checks if the input grid data file has a grid that matchs the new grid data file
        call GetData(Me%BatimInMatch,                                                    &
                     Me%ObjEnterData, flag,                                              &
                     SearchType   = FromFile_,                                           &
                     Default      = .true.,                                             &
                     keyword      ='MATCH_IN',                                           &
                     ClientModule ='DigitalTerrainCreator',                              &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR230'

        if (Me%BatimInON) then

            !Check if the user wants to extract a Window
            call GetData(Me%WindowOutput,                                               &
                         Me%ObjEnterData, flag,                                         &
                         SearchType   = FromFile_,                                      &
                         ClientModule ='DigitalTerrainCreator',                         &
                         keyword      ='WINDOW_OUTPUT',                                 &
                         Default      = .false.,                                        &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR240'

            if (Me%WindowOutput .and. .not. Me%BatimInMatch)                            &
                stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR250'

            if (Me%WindowOutput) then

                !Defines the limits of the window
                call GetData(Me%WindowLimits,                                           &
                             Me%ObjEnterData, flag,                                     &
                             SearchType   = FromFile_,                                  &
                             ClientModule ='DigitalTerrainCreator',                     &
                             keyword      ='WINDOW_LIMITS',                             &
                             STAT         = STAT_CALL)        
                if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR260'

            endif

        endif


        call GetData(Me%FillMethod,                                                     &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='FILL_METHOD',                                       &
                     Default      = Average_,                                           &
                     ClientModule ='DigitalTerrainCreator',                             &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR100'
        
        if (Me%FillMethod .ne. Average_ .and.   &
            Me%FillMethod .ne. Minimum_ .and.   &
            Me%FillMethod .ne. Maximum_) then
            
            write (*,*)
            write (*,*) 'FILL_METHOD can only be 1 (average), 2 (minimum) or 3 (maximum)'
            stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR110'
            
        endif
        

        call GetData(Me%ExportXYZ,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='EXPORT_XYZ',                                        &
                     Default      = .false.,                                            &
                     ClientModule ='DigitalTerrainCreator',                             &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR110'

        
        call GetData(Me%Convert2Xbeach,                                                 &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='CONVERT_2_XBEACH',                                  &
                     Default      = .false.,                                            &
                     ClientModule ='DigitalTerrainCreator',                             &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR120'        
        

        call GetData(Me%DD%ON,                                                          &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='DOMAIN_DECOMPOSITION',                              &
                     Default      = .false.,                                            &
                     ClientModule ='DigitalTerrainCreator',                             &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR130'    
        
        
        call GetData(Me%OutputSupportFiles,                                             &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='OUTPUT_SUPPORT_FILES',                              &
                     Default      = .false.,                                            &
                     ClientModule ='DigitalTerrainCreator',                             &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ConstructGlobalOptions - DigitalTerrainCreator - ERR140'    
        
        if (Me%DD%ON) then
            call Construct_DD            
        endif
        
        !call Construct_AllBorder
        

    end subroutine ConstructGlobalOptions
    
    !--------------------------------------------------------------------------

    subroutine Construct_DD

        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL
        integer                             :: i, j, nlines, ncolumns, iAux, jAux, n
        integer                             :: ILB, IUB, JLB, JUB
        type (T_Border),        pointer     :: BorDerX        
        !Begin-----------------------------------------------------------------
        
        call GetData(Me%DD%nJX,                                                         &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='DD_COLUMNS',                                        &
                     ClientModule ='DigitalTerrainCreator',                             &
                     Default      = 1,                                                  &   
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'Construct_DD - DigitalTerrainCreator - ERR010'
        
        if (flag <= 0) then
            stop 'Construct_DD - DigitalTerrainCreator - ERR020'
        endif
        
        if (Me%DD%nJX <= 0) then
            stop 'Construct_DD - DigitalTerrainCreator - ERR030'
        endif        
        
        call GetData(Me%DD%nIY,                                                         &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='DD_LINES',                                          &
                     ClientModule ='DigitalTerrainCreator',                             &
                     Default      = 1,                                                  &   
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'Construct_DD - DigitalTerrainCreator - ERR040'
        
        if (flag <= 0) then
            stop 'Construct_DD - DigitalTerrainCreator - ERR050'
        endif        
        
        if (Me%DD%nIY <= 0) then
            stop 'Construct_DD - DigitalTerrainCreator - ERR060'
        endif            
        
        Me%DD%nD = Me%DD%nJX * Me%DD%nIY
        
        allocate(Me%DD%Border      (1:Me%DD%nD))

        allocate(Me%DD%BorderLimits(1:Me%DD%nD))
        
        call Read_Lock_External_Var
        
        nlines      = Me%ExtVar%WorkSize%IUB - Me%ExtVar%WorkSize%ILB + 1
        ncolumns    = Me%ExtVar%WorkSize%JUB - Me%ExtVar%WorkSize%JLB + 1

        iAux        = nlines   / Me%DD%nIY
        jAux        = ncolumns / Me%DD%nJX
        
        if (iAux < 10) then
            stop 'Construct_DD - DigitalTerrainCreator - ERR050'
        endif
        
        if (jAux < 10) then
            stop 'Construct_DD - DigitalTerrainCreator - ERR060'
        endif        
        
        n = 1
        
        do i = 1, Me%DD%nIY
        do j = 1, Me%DD%nJX
            
            ILB =  Me%ExtVar%WorkSize%ILB     + iAux * (i - 1)
            IUB =  Me%ExtVar%WorkSize%ILB     + iAux * (i    ) - 1
            
            JLB =  Me%ExtVar%WorkSize%JLB     + jAux * (j - 1)
            JUB =  Me%ExtVar%WorkSize%JLB     + jAux * (j    ) - 1
            
            if (i == Me%DD%nIY) then
                IUB =  Me%ExtVar%WorkSize%IUB
            endif
            
            if (j == Me%DD%nJX) then
                JUB =  Me%ExtVar%WorkSize%JUB
            endif            
            
            Me%DD%BorderLimits(n)%ILB = ILB
            Me%DD%BorderLimits(n)%IUB = IUB
            
            Me%DD%BorderLimits(n)%JLB = JLB
            Me%DD%BorderLimits(n)%JUB = JUB     
            
            n = n + 1
        enddo
        enddo
        
        do n = 1, Me%DD%nD
            
            allocate(BorDerX)
            
            call PolygonBoundGridCurvSet(XX2D         = Me%ExtVar%XX_IE,                &
                                         YY2D         = Me%ExtVar%YY_IE,                & 
                                         GridBorder   = BorDerX,                        &
                                         ILB          = Me%DD%BorderLimits(n)%ILB,      &
                                         IUB          = Me%DD%BorderLimits(n)%IUB,      &
                                         JLB          = Me%DD%BorderLimits(n)%JLB,      &
                                         JUB          = Me%DD%BorderLimits(n)%JUB)
            
            Me%DD%Border(n) = BorDerX
            
            deallocate(BorDerX)            
            
        enddo
        
        call Read_UnLock_External_Var
        
    end subroutine Construct_DD
    
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    subroutine Construct_AllBorder

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB
        !Begin-----------------------------------------------------------------
        
        call Read_Lock_External_Var
            
        ILB =  Me%ExtVar%WorkSize%ILB
        IUB =  Me%ExtVar%WorkSize%IUB
            
        JLB =  Me%ExtVar%WorkSize%JLB
        JUB =  Me%ExtVar%WorkSize%JUB
        
        allocate(Me%AllBorDer)
        
        call PolygonBoundGridCurvSet(XX2D         = Me%ExtVar%XX_IE,                    &
                                     YY2D         = Me%ExtVar%YY_IE,                    & 
                                     GridBorder   = Me%AllBorDer,                       &
                                     ILB          = ILB,                                &
                                     IUB          = IUB,                                &
                                     JLB          = JLB,                                &
                                     JUB          = JUB)
            
        
        call Read_UnLock_External_Var
        
    end subroutine Construct_AllBorder
    
    !--------------------------------------------------------------------------
    

    
    !--------------------------------------------------------------------------
    
    subroutine ConstructSelectedArea

        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL
        character(len=line_length)          :: AreaFilePath
        character(len=line_length)          :: FullBufferLine
        integer                             :: ClientNumber, StartLine, EndLine
        integer                             :: CurrentLineNumber
        logical                             :: BlockFound
        !Begin-----------------------------------------------------------------

        call ExtractBlockFromBuffer(Me%ObjEnterData,                         &
                                    ClientNumber,                           &
                                    '<BeginAreaFiles>',                     &
                                    '<EndAreaFiles>',                       &
                                    BlockFound,                             &
                                    STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructSelectedArea - DigitalTerrainCreator - ERR130'

        if (BlockFound) then

            call GetBlockSize(Me%ObjEnterData,                              &
                              ClientNumber,                                 &
                              StartLine,                                    &
                              EndLine,                                      &
                              FromBlock_,                            &
                              STAT = STAT_CALL)

            do CurrentLineNumber = StartLine + 1 , EndLine - 1

                call GetFullBufferLine(Me%ObjEnterData,                     &
                                       CurrentLineNumber,                   &
                                       FullBufferLine,                      &
                                       STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop 'ConstructSelectedArea - DigitalTerrainCreator - ERR140'

                !Get LandArea files paths
                call GetData(AreaFilePath,                                  &
                             Me%ObjEnterData,                               &
                             flag,                                          &
                             SearchType   = FromBlockInBlock_,              &
                             Buffer_Line  = CurrentLineNumber,              &
                             STAT         = STAT_CALL)        
                if(STAT_CALL .ne. SUCCESS_) stop 'ConstructSelectedArea - DigitalTerrainCreator - ERR150'

                !Construct overlapping Area collection
                call New(Me%SectionAreas,  trim(AreaFilePath))

            end do

        endif

        
        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructSelectedArea - DigitalTerrainCreator - ERR160'

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructSelectedArea - DigitalTerrainCreator - ERR170'        
    
    end subroutine ConstructSelectedArea
    
    !--------------------------------------------------------------------------
    
    subroutine ReadTriangulationOptions
        
        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL

        !Begin-----------------------------------------------------------------
        
        !Reads if the user wants to use the closest points for points which lie
        !outside the polygon
        call GetData (Me%Triang%FillOutsidePoints, Me%ObjEnterData, flag,           &
                      SearchType   = FromFile_,                                     &
                      keyword      ='FILL_OUTSIDE_POINTS',                          &
                      default      = .false.,                                       &
                      ClientModule ='DigitalTerrainCreator',                        &
                      STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_ ) stop 'ReadTriangulationOptions - DigitalTerrainCreator - ERR10'


        !Reads the tolerance distance between nodes
        call GetData (Me%Triang%Tolerance, Me%ObjEnterData, flag,                   &
                      SearchType   = FromFile_,                                     &
                      keyword      ='TOLERANCE',                                    &
                      default      = 0.0,                                           &
                      ClientModule ='DigitalTerrainCreator',                        &
                      STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadTriangulationOptions - DigitalTerrainCreator - ERR20'


        !Reads the tolerance distance between nodes
        call GetData (Me%Triang%OutputTriangles, Me%ObjEnterData, flag,             &
                      SearchType   = FromFile_,                                     &
                      keyword      ='WRITE_TRIANGLES',                              &
                      default      = .false.,                                       &
                      ClientModule ='DigitalTerrainCreator',                        &
                      STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadTriangulationOptions - DigitalTerrainCreator - ERR30'


        !Reads the tolerance distance between nodes
        call GetData (Me%Triang%File, Me%ObjEnterData, flag,                        &
                      SearchType   = FromFile_,                                     &
                      keyword      ='TRIANGLES_FILE',                               &
                      ClientModule ='DigitalTerrainCreator',                        &
                      STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadTriangulationOptions - DigitalTerrainCreator - ERR40'

    end subroutine ReadTriangulationOptions
    
    
    !--------------------------------------------------------------------------

    subroutine ConstructOverlapping
        
        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL
        character(len=line_length)          :: AreaFilePath
        character(len=line_length)          :: FullBufferLine, AuxChar
        integer                             :: ClientNumber, StartLine, EndLine
        integer                             :: CurrentLineNumber, CountLevel, AuxLevel, ic
        logical                             :: BlockFound
        logical, dimension(:), pointer      :: LevelIn
        real,    dimension(3)               :: Aux3
        !Begin-----------------------------------------------------------------

        
        !Get grid file path
        call GetData(Me%Overlapping%Number,                                              &
                     Me%ObjEnterData, flag,                                              &
                     SearchType   = FromFile_,                                           &
                     keyword      ='OVERLAPPING_NUMBER',                                 &
                     ClientModule ='DigitalTerrainCreator',                              &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructOverlapping - DigitalTerrainCreator - ERR10'

        if (flag == 0) stop 'ConstructOverlapping - DigitalTerrainCreator - ERR20'

        allocate (Me%Overlapping%DataInfo(Me%Overlapping%Number))
        allocate (LevelIn(Me%Overlapping%Number))

        LevelIn(:) = .false.

        CountLevel = 0

        do 
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                 &
                                        ClientNumber,                                    &
                                        '<BeginGridDataInfo>',                           &
                                        '<EndGridDataInfo>',                             &
                                        BlockFound,                                      &
                                        STAT = STAT_CALL)
        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructOverlapping - DigitalTerrainCreator - ERR30'
    

ibf:        if (BlockFound) then

                CountLevel = CountLevel + 1


                call GetData(AuxLevel,                                                   &
                             Me%ObjEnterData, flag,                                      &
                             SearchType   = FromBlock_,                                  &
                             keyword      ='LEVEL',                                      &
                             ClientModule ='DigitalTerrainCreator',                      &
                             STAT         = STAT_CALL)        
                if(STAT_CALL /= SUCCESS_) stop 'ConstructOverlapping - DigitalTerrainCreator - ERR40'

                if (flag == 0) stop 'ConstructOverlapping - DigitalTerrainCreator - ERR50'

                if (AuxLevel < 1 .or. AuxLevel > Me%Overlapping%Number)                  &
                               stop 'ConstructOverlapping - DigitalTerrainCreator - ERR60'

                if (LevelIn(AuxLevel)) then

                    stop 'ConstructOverlapping - DigitalTerrainCreator - ERR70'

                else

                    LevelIn(AuxLevel) = .true.

                endif

                call GetData(AuxChar,                                                    &
                             Me%ObjEnterData, flag,                                      &
                             SearchType   = FromBlock_,                                  &
                             keyword      ='OVERLAP_TYPE',                               &
                             Default      = trim(Char_GridDataType),                     &
                             ClientModule ='DigitalTerrainCreator',                      &
                             STAT         = STAT_CALL)        

                if      (Auxchar == trim(Char_GridDataType)) then

                    Me%Overlapping%DataInfo(AuxLevel)%InfoType = GridDataType

                else if (Auxchar == trim(Char_CellsType   )) then

                    Me%Overlapping%DataInfo(AuxLevel)%InfoType = CellsType

                else

                    stop 'ConstructOverlapping - DigitalTerrainCreator - ERR80'

                endif


ift:            if (Me%Overlapping%DataInfo(AuxLevel)%InfoType == GridDataType) then


                    call GetData(AuxChar,                                               &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType   = FromBlock_,                             &
                                 keyword      ='GRIDATA_FILE',                          &
                                 ClientModule ='DigitalTerrainCreator',                 &
                                 STAT         = STAT_CALL)        
                    if(STAT_CALL /= SUCCESS_) stop 'ConstructOverlapping - DigitalTerrainCreator - ERR80'

                    if (flag == 0) stop 'ConstructOverlapping - DigitalTerrainCreator - ERR90'


                    call ConstructGridData (Me%Overlapping%DataInfo(AuxLevel)%ObjGridData, Me%ObjGrid,  &
                                            FileName = AuxChar, STAT = STAT_CALL)        
                    if(STAT_CALL /= SUCCESS_) stop 'ConstructOverlapping - DigitalTerrainCreator - ERR100'

                    call GetData(Me%Overlapping%DataInfo(AuxLevel)%Hmin,                                                           &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType   = FromBlock_,                             &
                                 keyword      ='H_MIN',                                 &
                                 ClientModule ='DigitalTerrainCreator',                 &
                                 default      = -1e9,                                   &
                                 STAT         = STAT_CALL)        
                    if(STAT_CALL /= SUCCESS_) stop 'ConstructOverlapping - DigitalTerrainCreator - ERR110'


                    call GetData(Me%Overlapping%DataInfo(AuxLevel)%Hmax,                                                           &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType   = FromBlock_,                             &
                                 keyword      ='H_MAX',                                 &
                                 ClientModule ='DigitalTerrainCreator',                 &
                                 default      =  1e9,                                   &
                                 STAT         = STAT_CALL)        
                    if(STAT_CALL /= SUCCESS_) stop 'ConstructOverlapping - DigitalTerrainCreator - ERR120'


                    call GetData(Me%Overlapping%DataInfo(AuxLevel)%Percentage,          &
                                 Me%ObjEnterData, flag,                                 &
                                 SearchType   = FromBlock_,                             &
                                 keyword      ='PERCENTAGE',                            &
                                 ClientModule ='DigitalTerrainCreator',                 &
                                 default      =  1.,                                    &
                                 STAT         = STAT_CALL)        
                    if(STAT_CALL /= SUCCESS_) stop 'ConstructOverlapping - DigitalTerrainCreator - ERR120'

                    call ExtractBlockFromBlock(Me%ObjEnterData,                         &
                                                ClientNumber,                           &
                                                '<<BeginAreaFiles>>',                   &
                                                '<<EndAreaFiles>>',                     &
                                                BlockFound,                             &
                                                STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - DigitalTerrainCreator - ERR130'

                    if (BlockFound) then
    
    
                        call GetBlockSize(Me%ObjEnterData,                              &
                                          ClientNumber,                                 &
                                          StartLine,                                    &
                                          EndLine,                                      &
                                          FromBlockInBlock_,                            &
                                          STAT = STAT_CALL)

                        do CurrentLineNumber = StartLine + 1 , EndLine - 1

                            call GetFullBufferLine(Me%ObjEnterData,                     &
                                                   CurrentLineNumber,                   &
                                                   FullBufferLine,                      &
                                                   STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - DigitalTerrainCreator - ERR140'
        
                            !Get LandArea files paths
                            call GetData(AreaFilePath,                                  &
                                         Me%ObjEnterData,                               &
                                         flag,                                          &
                                         SearchType   = FromBlockInBlock_,              &
                                         Buffer_Line  = CurrentLineNumber,              &
                                         STAT         = STAT_CALL)        
                            if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - DigitalTerrainCreator - ERR150'
        
                            !Construct overlapping Area collection
                            call New(Me%Overlapping%DataInfo(AuxLevel)%Areas,  trim(AreaFilePath))

                        end do

                    endif

                else if (Me%Overlapping%DataInfo(AuxLevel)%InfoType == CellsType) then ift

                    call ExtractBlockFromBlock(Me%ObjEnterData,                         &
                                                ClientNumber,                           &
                                                '<<BeginCells>>',                       &
                                                '<<EndCells>>',                         &
                                                BlockFound,                             &
                                                STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - DigitalTerrainCreator - ERR160'

                    if (BlockFound) then
    
    
                        call GetBlockSize(Me%ObjEnterData,                              &
                                          ClientNumber,                                 &
                                          StartLine,                                    &
                                          EndLine,                                      &
                                          FromBlockInBlock_,                            &
                                          STAT = STAT_CALL)

                        Me%Overlapping%DataInfo(AuxLevel)%CellsNumber = EndLine - StartLine - 1

                        allocate(Me%Overlapping%DataInfo(AuxLevel)%Cells(Me%Overlapping%DataInfo(AuxLevel)%CellsNumber))

                        ic = 1

                        do CurrentLineNumber = StartLine + 1 , EndLine - 1

                            call GetFullBufferLine(Me%ObjEnterData,                     &
                                                   CurrentLineNumber,                   &
                                                   FullBufferLine,                      &
                                                   STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - DigitalTerrainCreator - ERR170'
        
                            !Get cells depths
                            call GetData(Aux3,                                          &
                                         Me%ObjEnterData,                               &
                                         flag,                                          &
                                         SearchType   = FromBlockInBlock_,              &
                                         Buffer_Line  = CurrentLineNumber,              &
                                         STAT         = STAT_CALL)        
                            if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - DigitalTerrainCreator - ERR180'
        
                            Me%Overlapping%DataInfo(AuxLevel)%Cells(ic)%i     = int(Aux3(1))
                            Me%Overlapping%DataInfo(AuxLevel)%Cells(ic)%j     = int(Aux3(2)) 
                            Me%Overlapping%DataInfo(AuxLevel)%Cells(ic)%Depth =     Aux3(3) 

                            ic = ic + 1
                        end do

                    endif


                endif ift

            else ibf

                if (CountLevel /= Me%Overlapping%Number)                                &
                    stop 'ConstructOverlapping - DigitalTerrainCreator - ERR190'

                exit

            endif ibf

        enddo

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - DigitalTerrainCreator - ERR200'


        deallocate (LevelIn)

    end subroutine ConstructOverlapping
    

    !--------------------------------------------------------------------------
    
    
    subroutine RunProject

        !Begin-----------------------------------------------------------------
       
        write(*,*)"Running..."

        call Read_Lock_External_Var

        call AllocateTempVariables
        
       
        if (Me%BatimInON) then
            call ConstructInitialGridData
        else
            Me%Depth(:,:) = Me%NoDataPoint
        endif

        call DefineGridPoints

        !if true, XYZ data incomplete (only original section margins). 
        !Need to add more points XY to define section and compute Z.
        !Needs to be after BatimInON to fill bathim
        if (Me%ComputeSections) then
            
            !first, check if user wants to increase the number of margin points (interpolating between original)
            if (Me%CreateMarginPoints) then
                call AddMarginPoints
            endif
            
            !copy the original Z's (mid point Z's) because they will be changed in FillingPoints
            !to the DTM values and this info needs to be saved (after adding Margins or not)
            if (Me%SectionMiddlePointMethod == DefinedLevel_) then
                Me%CurrentXYZPoints => Me%XYZPoints
                do while(associated(Me%CurrentXYZPoints))
                    call New(Me%XYZPointsOriginal, Me%CurrentXYZPoints)
                    Me%CurrentXYZPoints => Me%CurrentXYZPoints%Next
                end do
            endif
            
            !second, fill margin points Z with known grid data
            call FillingPoints        
            
            !third, compute section points (position and Z)
            call ComputeSections
            
            !if interpolation is next needs to have no data points
            if (Me%InterpolType == Triangulation) call SelectAreaForInterpolation
        endif

        if (Me%PointsForInterpolation /= OnlyXYZPoints) then
!            if (.not.Me%BatimInOn .or. (Me%BatimInOn .and. .not. Me%BatimInMatch)) then

                if (Me%ExtVar%GridDistortion) then
                    call FillingCells
                else
                    call FillingRegularGrid
                endif

                if(Me%UseSmooth)then
                    call SmoothGridData
                end if

!            endif
        endif

        if(Me%Overlapping%ON)then
            call OverlappingGridData
        end if

        select case (Me%InterpolType)

            case (NoInterpolation)        
                
            case (Triangulation)
                call SelectPoints
                if (Me%DD%ON) then
                    call Triangulator_DD
                else
                    call Triangulator
                endif

                if (Me%River%CanonicSpace) then
                    if (Me%River%MaintainBed) then
                        call DredgeRiverBed
                    endif
                endif

            !case (Spline)
                !SelectPoints
                !call Spline
            case default
                stop 'Unknown interpolation method - RunProject - DigitalTerrainCreator - ERR20'
        end select

        if (Me%FilterSpikes) then
            call DoFilterSpikes
        endif

        if (Me%Filter) then
            call DoFilterDepths
        endif
        
        if (Triangulation) then
            call FillingLand
        endif

        if (Me%WindowOutPut) then
            call WriteWindow
        else
            call WriteDigitalTerrain
        endif

        call Read_UnLock_External_Var

        call DeallocateTempVariables


    end subroutine RunProject

    !--------------------------------------------------------------------------

    subroutine AddMarginPoints
    
        !Local-----------------------------------------------------------------
        integer                                          :: NumberOfPoints
        real,    dimension(:  ), pointer                 :: XVector, YVector, ZVector
        integer                                          :: ipoint, CurrentPoint, NextPoint
        integer                                          :: CurrentMarginPoint
        real                                             :: CurrentX_M1, CurrentY_M1, NextX_M1, NextY_M1
        real                                             :: CurrentX_M2, CurrentY_M2, NextX_M2, NextY_M2
        real                                             :: TotalLength_M1, CurrentLength_M1, StepX_M1, StepY_M1
        real                                             :: TotalLength_M2, CurrentLength_M2, StepX_M2, StepY_M2        
        type(T_XYZPoints),                 pointer       :: XYZPoints
        !Begin-----------------------------------------------------------------
               
        !second, add new points 
        !loop through all xyz points
        Me%CurrentXYZPoints => Me%XYZPoints
        
        write (*,*)
        write (*,*) 'Warning: Using the method to create more margin points' 
        write (*,*) 'between given points. The point collection need to assure'
        write (*,*) 'that the first points of each pair belong to same margin'
        write (*,*) 'and same for the second points of each pair in opposite margin'
        write (*,*) 'WARNING 001 - DigitalTerrainCreator - AddMarginPoints - WARNING 001'
        write (*,*)
        
        !loop trough point collections
        do while(associated(Me%CurrentXYZPoints))
            
            !points to add (number of points defined by the user in the intervals between points)
            NumberOfPoints = ((Me%CurrentXYZPoints%Count / 2.) - 1) * Me%MarginPointsToAdd * 2
            
            allocate(XVector(NumberOfPoints))
            allocate(YVector(NumberOfPoints))                       
            allocate(ZVector(NumberOfPoints))
            
            ipoint = 1
            
            !Need to do the pair on each margin on the same time so that the location in the 
            !new list the points mantained consecutive
            !jump every two points because the original points are always in pairs     
            do CurrentPoint = 1, (Me%CurrentXYZPoints%Count - 3), 2
                
                !Need to duplicate variables for M1 (Margin 1) and M2 (Margin 2)
                !define margin beggin (current point). 
                !and margin end (next point + 1 or currentPoint + 2 to be in the same margin)
                NextPoint   = CurrentPoint + 2
                CurrentX_M1 = Me%CurrentXYZPoints%X(CurrentPoint)
                CurrentX_M2 = Me%CurrentXYZPoints%X(CurrentPoint + 1)
                CurrentY_M1 = Me%CurrentXYZPoints%Y(CurrentPoint)
                CurrentY_M2 = Me%CurrentXYZPoints%Y(CurrentPoint + 1)
                NextX_M1    = Me%CurrentXYZPoints%X(NextPoint)
                NextX_M2    = Me%CurrentXYZPoints%X(NextPoint + 1)
                NextY_M1    = Me%CurrentXYZPoints%Y(NextPoint)
                NextY_M2    = Me%CurrentXYZPoints%Y(NextPoint + 1)
                
                !define length of margin
                TotalLength_M1 = sqrt((NextX_M1 - CurrentX_M1)**2 + (NextY_M1 - CurrentY_M1)**2)
                TotalLength_M2 = sqrt((NextX_M2 - CurrentX_M2)**2 + (NextY_M2 - CurrentY_M2)**2)
                
                !define step to compute new xyz points (same number as points needed - always odd
                StepX_M1 = (NextX_M1 - CurrentX_M1) / (Me%MarginPointsToAdd + 1)
                StepX_M2 = (NextX_M2 - CurrentX_M2) / (Me%MarginPointsToAdd + 1)
                StepY_M1 = (NextY_M1 - CurrentY_M1) / (Me%MarginPointsToAdd + 1)
                StepY_M2 = (NextY_M2 - CurrentY_M2) / (Me%MarginPointsToAdd + 1)
                
                !go trough all new points to define their coordinates
                do CurrentMarginPoint = 1, Me%MarginPointsToAdd
                    
                    !new point coordinates
                    XVector (ipoint)     = CurrentX_M1 + StepX_M1 * CurrentMarginPoint
                    XVector (ipoint + 1) = CurrentX_M2 + StepX_M2 * CurrentMarginPoint
                    YVector (ipoint)     = CurrentY_M1 + StepY_M1 * CurrentMarginPoint
                    YVector (ipoint + 1) = CurrentY_M2 + StepY_M2 * CurrentMarginPoint
                    
                    !Define Z if middle point Z is given by user
                    if (Me%SectionMiddlePointMethod == DefinedLevel_) then
                        
                        CurrentLength_M1 = sqrt((XVector (ipoint) - CurrentX_M1)**2 + (YVector (ipoint) - CurrentY_M1)**2)
                        if (CurrentLength_M1 > TotalLength_M1) stop 'DigitalTerrainCreator - AddMarginPoints - ERR01'
                        
                        CurrentLength_M2 = sqrt((XVector (ipoint + 1) - CurrentX_M2)**2 + (YVector (ipoint + 1) - CurrentY_M2)**2)
                        if (CurrentLength_M2 > TotalLength_M2) stop 'DigitalTerrainCreator - AddMarginPoints - ERR02'
                       
                        ZVector (ipoint)     = LinearInterpolation (0., Me%CurrentXYZPoints%Z(CurrentPoint),          &
                                                                    TotalLength_M1, Me%CurrentXYZPoints%Z(NextPoint), &
                                                                    CurrentLength_M1)
                        
                        ZVector (ipoint + 1) = LinearInterpolation (0., Me%CurrentXYZPoints%Z(CurrentPoint + 1),          &
                                                                    TotalLength_M2, Me%CurrentXYZPoints%Z(NextPoint + 1), &
                                                                    CurrentLength_M2)     
                    else
                        ZVector (ipoint)     = 0.0
                        ZVector (ipoint + 1) = 0.0
                    endif
                    
                    !go for next point (two points were just made)
                    ipoint = ipoint + 2
                    
                enddo
                
            end do
            
            !Construct XYZ points collection. For now only local so that this process is not inifinite
            !if the collection was added to Me%XYZPoints it would be next evaluated again
            call New(XYZPoints, XVector, YVector, ZVector)

            deallocate(XVector)
            deallocate(YVector)                       
            deallocate(ZVector)
            
            Me%CurrentXYZPoints => Me%CurrentXYZPoints%Next
            
        end do
        
        !only after the computation add to the Me%XYZPoints or the previous cycle would be infinite
        !and the margin points added will be in pairs so that the section creation is correct
        Me%CurrentXYZPoints => XYZPoints
        do while(associated(Me%CurrentXYZPoints))
        
            call New(Me%XYZPoints, Me%CurrentXYZPoints)
            
            Me%CurrentXYZPoints => Me%CurrentXYZPoints%Next
            
        enddo
        
        !nullify    (CurrentXYZPoints)
        deallocate (XYZPoints)
        
    
    end subroutine AddMarginPoints
    
    !--------------------------------------------------------------------------
    
    subroutine ComputeSections
    
        !Local-----------------------------------------------------------------
        integer                                          :: NumberOfPoints
        real,    dimension(:  ), pointer                 :: XVector, YVector, ZVector
        integer                                          :: ipoint, CurrentPoint, NextPoint
        integer                                          :: CurrentSectionPoint
        integer                                          :: UnitNumber
        real                                             :: CurrentX, CurrentY, NextX, NextY
        real                                             :: MiddlePointX, MiddlePointY, HalfWidth
        real                                             :: DistanceToMiddle, Factor, StepX, StepY
        integer                                          :: STAT_CALL
        type(T_XYZPoints),                 pointer       :: XYZPoints !, XYZPointsOriginal  !local temp storage 
        real                                             :: Exponent, SectionHeight
        !Begin-----------------------------------------------------------------
        
        
        !open the file to write down the points
        call UnitsManager (UnitNumber, OPEN_FILE, STAT = STAT_CALL)
        open (unit=UnitNumber, status = 'unknown', file = trim(adjustl(Me%SectionXYZFile)))

        !write the new points group
        write(UnitNumber,*)'<begin_xyz>'
        
        !second, add new points and compute Z
        !loop through all xyz points
        Me%CurrentXYZPoints => Me%XYZPoints
        
!        !get original z's (same point collections as Me%XYZPoints but original z's)
!        if (Me%SectionMiddlePointMethod == DefinedLevel_) then
!            XYZPointsOriginal => Me%XYZPointsOriginal
!        endif

         !Only the first collection for now. Because the New fucntion will add a new to the list
         !and then it is infinite adding a new one and analyzing it and adding
        do while(associated(Me%CurrentXYZPoints))
            
            !verify that the xyz points are even so that every one has a pair (section start and end)
            if (MOD(Me%CurrentXYZPoints%Count, 2) /= 0) then
                write (*,*)
                write (*,*) 'In section method need to provide even number of'
                write (*,*) 'xyz points since they define section start and end.'
                stop ' ComputeSections - DigitalTerrainCreator - ERR01'
            endif
           
            
            !points to add (number of points defined bt the user per pair of xyz points)
            NumberOfPoints = (Me%CurrentXYZPoints%Count / 2.) * Me%SectionPointsToAdd
            
            allocate(XVector(NumberOfPoints))
            allocate(YVector(NumberOfPoints))                       
            allocate(ZVector(NumberOfPoints))
            
            ipoint = 0
                       
            !go with step 2 because it is assumed that pairs of points
            !defining the section begin and end were created
            do CurrentPoint = 1, Me%CurrentXYZPoints%Count, 2
                
                !define section beggin (current point)
                !and section end (next point)
                NextPoint = CurrentPoint + 1
                CurrentX = Me%CurrentXYZPoints%X(CurrentPoint)
                CurrentY = Me%CurrentXYZPoints%Y(CurrentPoint)
                NextX    = Me%CurrentXYZPoints%X(NextPoint)
                NextY    = Me%CurrentXYZPoints%Y(NextPoint)
                
                !Begin section point
                write(UnitNumber,*)Me%CurrentXYZPoints%X(CurrentPoint), Me%CurrentXYZPoints%Y(CurrentPoint), &
                                    Me%CurrentXYZPoints%Z(CurrentPoint)
                                                    
                !define middle point (here the depth will be the one defined by the user)
                MiddlePointX = (CurrentX + NextX) / 2.
                MiddlePointY = (CurrentY + NextY) / 2.
                HalfWidth = sqrt((MiddlePointX - CurrentX)**2 + (MiddlePointY - CurrentY)**2)
                
                !define step to compute new xyz points (same sumber as points needed - always odd
                StepX = (NextX - CurrentX) / (Me%SectionPointsToAdd + 1)
                StepY = (NextY - CurrentY) / (Me%SectionPointsToAdd + 1)
                
                
                !go trough all new points to define their coordinates
                do CurrentSectionPoint = 1, Me%SectionPointsToAdd
                    
                    ipoint = ipoint + 1
                    !new point coordinates
                    XVector (ipoint) = CurrentX + StepX * CurrentSectionPoint
                    YVector (ipoint) = CurrentY + StepY * CurrentSectionPoint
                    
                    DistanceToMiddle = sqrt((MiddlePointX - XVector (ipoint))**2 + (MiddlePointY - YVector (ipoint))**2)
                    
                    !shape of the section to compute point Z (linear V, root U, power "heart") 
                    if (Me%SectionInterpolType == LinearInterpol_)then
                        Exponent = 1
                    else if (Me%SectionInterpolType == PowerInterpol_)then
                        Exponent = Me%SectionInterpolOrder
                    else if (Me%SectionInterpolType == RootInterpol_)then
                        Exponent = 1. / Me%SectionInterpolOrder
                    endif
                    
                    !factor of middle height to compute point Z
                    Factor = LinearInterpolation (HalfWidth, 0., 0., 1., DistanceToMiddle)**Exponent
                    
                    !if between first point and middle
                    if (CurrentSectionPoint .le. Me%SectionPointsToAdd / 2.) then
                        
                        if (Me%SectionMiddlePointMethod == DefinedLevel_) then
                        
                            SectionHeight = Me%CurrentXYZPoints%Z(CurrentPoint) -                  &
                                            Me%XYZPointsOriginal%Z(CurrentPoint)
                            if (SectionHeight .lt. 0.0) then
                                write(*,*) 'Negative section height. Make sure that the Z defined for the first point'
                                write(*,*) '(middle point level) is not higher than the section margin levels'
                                write(*,*) ' Occured in Pair: ',CurrentX, CurrentY, NextX, NextY
                                stop 'ComputeSections - DigitalTerrainCreator - ERR005'
                            endif
                        
                        elseif (Me%SectionMiddlePointMethod == DefinedHeight_) then
                            SectionHeight = Me%SectionHeight
                        endif
                        
                        !Point Z is the first point (section begin) minus the height computed
                        ZVector(ipoint) = Me%CurrentXYZPoints%Z(CurrentPoint) - Factor * SectionHeight
                    
                    else !if between middle and second point

                        if (Me%SectionMiddlePointMethod == DefinedLevel_) then
                        
                            SectionHeight = Me%CurrentXYZPoints%Z(NextPoint) -                      &
                                            Me%XYZPointsOriginal%Z(CurrentPoint)
                            if (SectionHeight .lt. 0.0) then
                                write(*,*) 'Negative section height. Make sure that the Z defined for the first point'
                                write(*,*) '(middle point level) is not higher than the section margin levels'
                                write(*,*) ' Occured in Pair: ',CurrentX, CurrentY, NextX, NextY
                                stop 'ComputeSections - DigitalTerrainCreator - ERR006'
                            endif
                        
                        elseif (Me%SectionMiddlePointMethod == DefinedHeight_) then
                        
                            SectionHeight = Me%CurrentXYZPoints%Z(NextPoint)      -   &
                                            (Me%CurrentXYZPoints%Z(CurrentPoint)  -   &
                                             Me%SectionHeight                   )
                                            
                            if (SectionHeight .lt. 0.0) then
                                write(*,*) 'Negative section height. Make sure that the second point of the section'
                                write(*,*) 'is not in an area lower than the middle point (height specified by the user)'
                                write(*,*) ' Occured in Pair: ',CurrentX, CurrentY, NextX, NextY
                                stop 'ComputeSections - DigitalTerrainCreator - ERR010'
                            endif
                        endif

                        !Point Z is the second point (section end) minus the height computed
                        ZVector(ipoint) = Me%CurrentXYZPoints%Z(Nextpoint) - Factor * SectionHeight

                    endif
                    
                    !added point
                    write(UnitNumber,*)XVector (ipoint), YVector (ipoint), ZVector(ipoint)
                    
                enddo
                
                !end section point
                write(UnitNumber,*)Me%CurrentXYZPoints%X(NextPoint), Me%CurrentXYZPoints%Y(NextPoint),       &
                                    Me%CurrentXYZPoints%Z(NextPoint)

                
            end do
            
            !Construct XYZ points collection. For now only local so that this process is not inifinite
            !if the collection was added to Me%XYZPoints it would be next evaluated again
            call New(XYZPoints, XVector, YVector, ZVector)

            deallocate(XVector)
            deallocate(YVector)                       
            deallocate(ZVector)
                        
            Me%CurrentXYZPoints => Me%CurrentXYZPoints%Next
            
            if (Me%SectionMiddlePointMethod == DefinedLevel_) then
                Me%XYZPointsOriginal => Me%XYZPointsOriginal%Next
            endif
            
        end do


        write(UnitNumber,*)'<end_xyz>'
        
        !close the xyz file
        call UnitsManager (UnitNumber, CLOSE_FILE, STAT = STAT_CALL)        
        
        
        !only after the computation add to the Me%XYZPoints or the previous cycle would be infinite
        Me%CurrentXYZPoints => XYZPoints
        do while(associated(Me%CurrentXYZPoints))
        
            call New(Me%XYZPoints, Me%CurrentXYZPoints)
            
            Me%CurrentXYZPoints => Me%CurrentXYZPoints%Next
            
        enddo
        
        !nullify    (CurrentXYZPoints)
        deallocate (XYZPoints)
        
    
    end subroutine ComputeSections
    
    !--------------------------------------------------------------------------

    subroutine FillingPoints

        !Local-----------------------------------------------------------------
        type (T_PointF),   pointer          :: GridPoint
        integer                             :: CurrentPoint
        integer                             :: i, j
        real                                :: XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW
        logical                             :: FoundPoint

        !Begin-----------------------------------------------------------------

        write(*,*)"Filling points with known data..."


        !loop through all xyz points
        Me%CurrentXYZPoints => Me%XYZPoints

        do while(associated(Me%CurrentXYZPoints))

            do CurrentPoint = 1, Me%CurrentXYZPoints%Count

                if(.not. Me%CurrentXYZPoints%Inside(CurrentPoint))then
                    
                    Me%AuxPoint%X = Me%CurrentXYZPoints%X(CurrentPoint)
                    Me%AuxPoint%Y = Me%CurrentXYZPoints%Y(CurrentPoint)

                    FoundPoint = .false.
                    
                    do i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB
                    do j = Me%ExtVar%WorkSize%JLB , Me%ExtVar%WorkSize%JUB
                                   
                        if (Me%ExtVar%DefineCellsMap(i, j)==1 .and. .not. FoundPoint) then

                            GridPoint => Me%GridPoint(i, j)
                            
                            XSW = Me%ExtVar%XX_IE(i, j)
                            YSW = Me%ExtVar%YY_IE(i, j)
                            XSE = Me%ExtVar%XX_IE(i, j + 1)
                            YSE = Me%ExtVar%YY_IE(i, j + 1)
                            XNE = Me%ExtVar%XX_IE(i + 1, j + 1)
                            YNE = Me%ExtVar%YY_IE(i + 1, j + 1)
                            XNW = Me%ExtVar%XX_IE(i + 1, j)
                            YNW = Me%ExtVar%YY_IE(i + 1, j)

                            Me%Rect%VerticesF(1)%X = XSW
                            Me%Rect%VerticesF(1)%Y = YSW
                            Me%Rect%VerticesF(2)%X = XSE
                            Me%Rect%VerticesF(2)%Y = YSE
                            Me%Rect%VerticesF(3)%X = XNE
                            Me%Rect%VerticesF(3)%Y = YNE
                            Me%Rect%VerticesF(4)%X = XNW
                            Me%Rect%VerticesF(4)%Y = YNW

                            Me%Rect%VerticesF(5)%X = Me%Rect%VerticesF(1)%X
                            Me%Rect%VerticesF(5)%Y = Me%Rect%VerticesF(1)%Y
                            
                            call SetLimits(Me%Rect)
                            
                            if(IsPointInsidePolygon(Me%AuxPoint, Me%Rect))then
                                Me%CurrentXYZPoints%Z(CurrentPoint) = Me%Depth(i,j)
                                FoundPoint = .true.
                            end if
                            
                        endif
                    enddo
                    enddo
                    
                    if (.not. FoundPoint) then
                        write (*,*)
                        write (*,*) 'Point x y outside grid',  Me%AuxPoint%X, Me%AuxPoint%Y
                        write (*,*) 'Please remove all points outside grid'
                        stop 'FillingPoints - DigitalTerrainTool - ERR10'
                    endif

                end if

            end do

            Me%CurrentXYZPoints => Me%CurrentXYZPoints%Next
        end do
        
    end subroutine FillingPoints


    !--------------------------------------------------------------------------

    subroutine SelectAreaForInterpolation

        !Local-----------------------------------------------------------------
        type (T_PointF),   pointer          :: GridPoint
        integer                             :: i, j    
        !Begin----------------------------------------------------------------
        
        do i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB
        do j = Me%ExtVar%WorkSize%JLB , Me%ExtVar%WorkSize%JUB

            if (Me%ExtVar%DefineCellsMap(i, j)==1) then

                GridPoint => Me%GridPoint(i, j)

                if (associated(Me%SectionAreas)) then
                
                    if (IsVisible(Me%SectionAreas, GridPoint)) then
                        Me%Depth(i,j) = Me%NoDataPoint
                    endif     
                
                endif
            endif
        enddo
        enddo
    
    end subroutine SelectAreaForInterpolation

    !--------------------------------------------------------------------------

    subroutine SelectPoints

        call SelectPointsWithData

        if (Me%OutputSupportFiles) then
            call WriteItem(Me%PointsWithData, Me%TotalWithData, Me%BaseAndComputedInformation)
        endif

        call SelectPointsWithNoData

        if (Me%OutputSupportFiles) then
            call WriteItem(Me%PointsNoData,   Me%TotalNoData,   Me%PointsWithNoInformation)
        endif

    end subroutine SelectPoints


    !--------------------------------------------------------------------------
    
    subroutine WriteDigitalTerrain
        
        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        character(len=StringLength)         :: Coment1, Coment2 
        
        !Begin-----------------------------------------------------------------

        write(*,*)"Writing bathymetry..."

        Coment1 = 'File generated by'
        Coment2 = 'Mohid Digital Terrain Creator'


        call WriteGridData(FileName         = trim(Me%BatimFilePathOut),        &
                           COMENT1          = Coment1,                          &
                           COMENT2          = Coment2,                          &
                           HorizontalGridID = Me%ObjGrid,                       &
                           FillValue        = Me%LandPoint,                     &
                           Overwrite        = .true.,                           &
                           GridData2D_Real  = Me%Depth,                         &
                           STAT             = STAT_CALL) 

            if (Me%ExportXYZ) then
                
                call Export_To_XYZ(Me%Depth, Me%LandPoint, Me%BatimFilePathOut)
                
            endif 
                        
        if (Me%Convert2Xbeach) then
                
            call Convert_2_Xbeach(Me%Depth, Me%BatimFilePathOut)
                
        endif         
                        

    end subroutine WriteDigitalTerrain
    
    
    !--------------------------------------------------------------------------
    subroutine Export_To_XYZ(Depth, LandPoint, Filename)

        !Arguments---------------------------------------------------------------
        real, dimension(:,:), pointer                   :: Depth
        real                                            :: LandPoint
        character(len=*)                                :: Filename        
        !Local-------------------------------------------------------------------
        real,  dimension(:,:), pointer                  :: CoordX, CoordY
        type (T_Size2D)                                 :: WorkSize
        integer                                         :: i, j, STAT_CALL, Unit, SplitByExtension
      
        
        !------------------------------------------------------------------------

        write(*,*)"Writing MOHID .xyz file..."

        call GetZCoordinates(Me%ObjGrid, CoordX, CoordY)
        
        call GetHorizontalGridSize(Me%ObjGrid, WorkSize = WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Export_To_XYZ - DigitalTerrainTool - ERR10'
              
        SplitByExtension = scan(trim(Filename),'.',Back = .true.)
        
        call UnitsManager (Unit, OPEN_FILE, STAT = STAT_CALL)        
                    
        open(   Unit   = Unit,                                                          &
                File   = trim(Filename(1:SplitByExtension-1)//'.xyz'),                  &
                Form   = 'FORMATTED',                                                   &
                STATUS = 'UNKNOWN',                                                     &
                Action = 'WRITE',                                                       &
                IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Export_To_XYZ - DigitalTerrainTool - ERR20'
                    
        write(Unit,*)'<begin_xyz>'
        
        do i = WorkSize%ILB, WorkSize%IUB
            do j = WorkSize%JLB, WorkSize%JUB
                if (Depth(i,j) == FillValueReal) then
                    write(Unit,*)CoordX(i, j), CoordY(i, j), LandPoint
                else
                    write(Unit,*)CoordX(i, j), CoordY(i, j), Depth(i,j)
                endif
            enddo
        enddo
        
        write(Unit,*)'<end_xyz>'
                        
        call UnitsManager(Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Export_To_XYZ - DigitalTerrainTool - ERR30'
        
    end subroutine Export_To_XYZ
 
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    subroutine Convert_2_Xbeach(Depth, Filename)

        !Arguments---------------------------------------------------------------
        real, dimension(:,:), pointer                   :: Depth
        character(len=*)                                :: Filename        
        !Local-------------------------------------------------------------------
        real,  dimension(:,:), pointer                  :: CoordX, CoordY
        type (T_Size2D)                                 :: WorkSize
        integer                                         :: i, j, STAT_CALL, Ub, Ux, Uy, SplitByExtension
      
        
        !------------------------------------------------------------------------

        write(*,*)"Writing MOHID XBEACH bathymetry files file..."

        call GetZCoordinates(Me%ObjGrid, CoordX, CoordY)
        
        call GetHorizontalGridSize(Me%ObjGrid, WorkSize = WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Convert_2_Xbeach - DigitalTerrainTool - ERR10'
              
        SplitByExtension = scan(trim(Filename),'.',Back = .true.)
        
        call UnitsManager (Ub, OPEN_FILE, STAT = STAT_CALL)
                    
        open(   Unit   = Ub,                                                            &
                File   = trim(Filename(1:SplitByExtension-1)//'.dep'),                  &
                Form   = 'FORMATTED',                                                   &
                STATUS = 'UNKNOWN',                                                     &
                Action = 'WRITE',                                                       &
                IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Convert_2_Xbeach - DigitalTerrainTool - ERR20'
        
                    
        do i = WorkSize%ILB, WorkSize%IUB
            write(Ub,'(10000f8.2)') (Depth (i, j),j=WorkSize%JLB,WorkSize%JUB)
        enddo
        
        call UnitsManager(Ub, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Convert_2_Xbeach - DigitalTerrainTool - ERR30'

        call UnitsManager (Ux, OPEN_FILE, STAT = STAT_CALL)
                    
        open(   Unit   = Ux,                                                            &
                File   = trim(Filename(1:SplitByExtension-1)//'_x.grd'),                &
                Form   = 'FORMATTED',                                                   &
                STATUS = 'UNKNOWN',                                                     &
                Action = 'WRITE',                                                       &
                IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Convert_2_Xbeach - DigitalTerrainTool - ERR40'
        
        do i = WorkSize%ILB, WorkSize%IUB
            write(Ux,'(10000f12.2)')  (CoordX(i, j),j=WorkSize%JLB,WorkSize%JUB)
        enddo
        
        call UnitsManager(Ux, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Convert_2_Xbeach - DigitalTerrainTool - ERR50'
        
        call UnitsManager (Uy, OPEN_FILE, STAT = STAT_CALL)        
        
        open(   Unit   = Uy,                                                            &
                File   = trim(Filename(1:SplitByExtension-1)//'_y.grd'),                &
                Form   = 'FORMATTED',                                                   &
                STATUS = 'UNKNOWN',                                                     &
                Action = 'WRITE',                                                       &
                IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Convert_2_Xbeach - DigitalTerrainTool - ERR60'
                    
        do i = WorkSize%ILB, WorkSize%IUB
            write(Uy,'(10000f12.2)') (CoordY(i, j),j=WorkSize%JLB,WorkSize%JUB)
        enddo
        
        call UnitsManager(Uy, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Convert_2_Xbeach - DigitalTerrainTool - ERR70'        
        
    end subroutine Convert_2_Xbeach
 
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine WriteWindow
        
        !Local-----------------------------------------------------------------

        real,    dimension(:,:), pointer    :: Aux2D, XX_IE, YY_IE, XX_IE_Aux, YY_IE_Aux 
        real,    dimension(:  ), pointer    :: XX, YY, XX_Aux, YY_Aux
        type (T_Size2D)                     :: WorkSize
        logical                             :: DistortionYes, Overwrite
        integer                             :: GEOG, UTM, MIL_PORT, SIMPLE_GEOG, NLRD
        integer                             :: GRID_COORD, CoordType, Zone
        integer                             :: STAT_CALL
        real                                :: Xorig, Yorig, GRID_ANGLE, Latitude, Longitude
        character(len=StringLength)         :: Coment1, Coment2 
        
        !Begin-----------------------------------------------------------------

        write(*,*)"Writing bathymetry..."

        Coment1 = 'File generated by'
        Coment2 = 'Mohid Digital Terrain Creator'


        call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,             &
                               SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD,      &
                               NLRD = NLRD)

        !Gets Coordinates in use
        call GetGridCoordType(Me%ObjGridIN, CoordType, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'WriteWindow - DigitalTerrainTool - ERR10'
    

        call GetCheckDistortion(Me%ObjGridIN, DistortionYes, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'WriteWindow - DigitalTerrainTool - ERR2'

        if    (CoordType == SIMPLE_GEOG)then
        
            call GetGridLatitudeLongitude(Me%ObjGridIn,                                 &
                                            GridLatitudeConn  = YY_IE,                  &
                                            GridLongitudeConn = XX_IE,                  &
                                            STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_)stop 'WriteWindow - DigitalTerrainTool - ERR40'

        elseif(CoordType == UTM .or. CoordType == MIL_PORT .or. &
                CoordType == GRID_COORD .or. CoordType == NLRD)then

            call GetHorizontalGrid(Me%ObjGridIn,                                        &
                                    XX_IE = XX_IE,                                      &
                                    YY_IE = YY_IE,                                      &
                                    STAT  = STAT_CALL)
            if(STAT_CALL /= SUCCESS_)stop 'WriteWindow - DigitalTerrainTool - ERR50'

        else

            write(*,*)'GEOG coordinate type cannot be used in digital terrain generation'
            stop 'WriteWindow - DigitalTerrainTool - ERR60'

        end if

        Xorig = XX_IE(Me%WindowLimits(1), Me%WindowLimits(2))
        Yorig = YY_IE(Me%WindowLimits(1), Me%WindowLimits(2))

        call GetHorizontalGrid(Me%ObjGridIn,                                            &
                                XX = XX,                                                &
                                YY = YY,                                                &
                                STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'WriteWindow - DigitalTerrainTool - ERR70'


        WorkSize%ILB = 1 
        WorkSize%JLB = 1 
        WorkSize%IUB = Me%WindowLimits(3) - Me%WindowLimits(1) + 1
        WorkSize%JUB = Me%WindowLimits(4) - Me%WindowLimits(2) + 1

        call GetGridAngle(Me%ObjGridIn, GRID_ANGLE, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'WriteWindow - DigitalTerrainTool - ERR80'

        call GetGridZone(Me%ObjGridIn, Zone, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'WriteWindow - DigitalTerrainTool - ERR90'

        call GetLatitudeLongitude(Me%ObjGridIn, Latitude, Longitude, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'WriteWindow - DigitalTerrainTool - ERR100'

        Overwrite = .true.

        allocate(Aux2D(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB))

        Aux2D(:,:) = Me%Depth(Me%WindowLimits(1):Me%WindowLimits(3), Me%WindowLimits(2):Me%WindowLimits(4))

        if (DistortionYes) then

            allocate(XX_IE_aux(WorkSize%ILB:WorkSize%IUB+1, WorkSize%JLB:WorkSize%JUB+1))
            allocate(YY_IE_aux(WorkSize%ILB:WorkSize%IUB+1, WorkSize%JLB:WorkSize%JUB+1))

            XX_IE_aux(:,:) = XX_IE(Me%WindowLimits(1):Me%WindowLimits(3) + 1, Me%WindowLimits(2):Me%WindowLimits(4) + 1)
            YY_IE_aux(:,:) = YY_IE(Me%WindowLimits(1):Me%WindowLimits(3) + 1, Me%WindowLimits(2):Me%WindowLimits(4) + 1)

            call WriteGridData(FileName             = trim(Me%BatimFilePathOut),        &
                               ConnectionX          = XX_IE_aux,                        &
                               ConnectionY          = YY_IE_aux,                        &
                               COMENT1              = COMENT1,                          &
                               COMENT2              = COMENT2,                          &
                               WorkSize             = WorkSize,                         &  
                               CoordType            = CoordType,                        &
                               Xorig                = Xorig,                            &
                               Yorig                = Yorig,                            &
                               Zone                 = Zone,                             &
                               GRID_ANGLE           = GRID_ANGLE,                       &
                               Latitude             = Latitude,                         &
                               Longitude            = Longitude,                        &
                               FillValue            = -99.,                             &
                               Overwrite            = Overwrite,                        &
                               GridData2D_Real      = Aux2D,                            &
                               STAT                 = STAT_CALL) 

            deallocate(XX_IE_aux)
            deallocate(YY_IE_aux)


        else

            allocate(XX_aux(WorkSize%JLB:WorkSize%JUB+1))
            allocate(YY_aux(WorkSize%ILB:WorkSize%IUB+1))

            XX_aux(:) = XX(Me%WindowLimits(2):Me%WindowLimits(4) + 1) - XX(Me%WindowLimits(2))
            YY_aux(:) = YY(Me%WindowLimits(1):Me%WindowLimits(3) + 1) - YY(Me%WindowLimits(1))


            call WriteGridData(FileName             = trim(Me%BatimFilePathOut),        &
                               XX                   = XX_aux,                           &
                               YY                   = YY_aux,                           &
                               COMENT1              = COMENT1,                          &
                               COMENT2              = COMENT2,                          &
                               WorkSize             = WorkSize,                         &  
                               CoordType            = CoordType,                        &
                               Xorig                = Xorig,                            &
                               Yorig                = Yorig,                            &
                               Zone                 = Zone,                             &
                               GRID_ANGLE           = GRID_ANGLE,                       &
                               Latitude             = Latitude,                         &
                               Longitude            = Longitude,                        &
                               FillValue            = -99.,                             &
                               Overwrite            = Overwrite,                        &
                               GridData2D_Real      = Aux2D,                            &
                               STAT                 = STAT_CALL) 

            deallocate(XX_aux)
            deallocate(YY_aux)


        endif

        deallocate(Aux2D)


        call UnGetHorizontalGrid(Me%ObjGridIn, XX_IE, STAT  = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'WriteWindow - DigitalTerrainTool - ERR200'

        call UnGetHorizontalGrid(Me%ObjGridIn, YY_IE, STAT =  STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'WriteWindow - DigitalTerrainTool - ERR210'

        if (.not.DistortionYes) then
            call UnGetHorizontalGrid(Me%ObjGridIn, XX, STAT  = STAT_CALL)
            if(STAT_CALL /= SUCCESS_)stop 'WriteWindow - DigitalTerrainTool - ERR220'

            call UnGetHorizontalGrid(Me%ObjGridIn, YY, STAT  = STAT_CALL)
            if(STAT_CALL /= SUCCESS_)stop 'WriteWindow - DigitalTerrainTool - ERR230'
        endif


    end subroutine WriteWindow
    
    
    !--------------------------------------------------------------------------



    subroutine ConstructInitialGridData
        
        !Local-----------------------------------------------------------------
        real,    dimension(:,:), pointer    :: Aux2D, XX_IE, YY_IE
        integer, dimension(:,:), pointer    :: DefineCellsMap
        real,    dimension(:  ), pointer    :: XVector, YVector, ZVector
        type (T_Size2D)                     :: WorkSize
        real                                :: XSW, XSE, XNE, XNW
        real                                :: YSW, YSE, YNE, YNW
        integer                             :: GEOG, UTM, MIL_PORT, SIMPLE_GEOG, NLRD
        integer                             :: GRID_COORD, CoordType
        integer                             :: STAT_CALL, i, j, ic, Counter

        !Begin-----------------------------------------------------------------

        if (Me%BatimInMatch) then

            Me%ObjGridIn = Me%ObjGrid

        else

            call ConstructHorizontalGrid(Me%ObjGridIn, DataFile = Me%BatimFilePathIn, STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR05'

        endif

        call ConstructGridData(Me%ObjGridDataIn, Me%ObjGridIn, FileName = Me%BatimFilePathIn, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR10'


        call GetGridData(Me%ObjGridDataIn, Aux2D, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR20'



        
        if (Me%BatimInMatch) then
            Me%Depth(:,:) = Aux2D(:,:)
        else
            Me%Depth(:,:) = Me%NoDataPoint 
            
            call GetHorizontalGridSize(Me%ObjGridIn, WorkSize = WorkSize, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR30'


            call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,             &
                                   SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD,      &
                                   NLRD = NLRD)

            !Gets Coordinates in use
            call GetGridCoordType(Me%ObjGridIN, CoordType, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR40'
        

            if    (CoordType == SIMPLE_GEOG)then
            
                call GetGridLatitudeLongitude(Me%ObjGridIn,                             &
                                              GridLatitudeConn  = YY_IE,                &
                                              GridLongitudeConn = XX_IE,                &
                                              STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR50'

            elseif(CoordType == UTM .or. CoordType == MIL_PORT .or.                     &
                   CoordType == GRID_COORD .or. CoordType == NLRD)then

                call GetHorizontalGrid(Me%ObjGridIn,                                    &
                                       XX_IE = XX_IE,                                   &
                                       YY_IE = YY_IE,                                   &
                                       STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR60'

            else

                write(*,*)'GEOG coordinate type cannot be used in digital terrain generation'
                stop 'ConstructInitialGridData - DigitalTerrainTool - ERR70'

            end if

            call GetDefineCellsMap(Me%ObjGridIn, DefineCellsMap, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR80'

            Counter = 0

            do i = WorkSize%ILB,  WorkSize%IUB
            do j = WorkSize%JLB , WorkSize%JUB

idef1:      if (DefineCellsMap(i, j)==1 .and. Aux2D(i, j) /= Me%LandPoint) then

                    Counter = Counter + 1

                endif idef1 
            enddo
            enddo

            Me%TotalXYZ = Me%TotalXYZ + Counter

            allocate(XVector(Counter))
            allocate(YVector(Counter))                       
            allocate(ZVector(Counter))


            ic = 1

            write(45,*) "<begin_xyz>"

            do i = WorkSize%ILB,  WorkSize%IUB
            do j = WorkSize%JLB , WorkSize%JUB

idef2:      if (DefineCellsMap(i, j)==1 .and. Aux2D(i, j) /= Me%LandPoint) then
                

                XSW = XX_IE(i, j)
                YSW = YY_IE(i, j)
                XSE = XX_IE(i, j + 1)
                YSE = YY_IE(i, j + 1)
                XNE = XX_IE(i + 1, j + 1)
                YNE = YY_IE(i + 1, j + 1)
                XNW = XX_IE(i + 1, j)
                YNW = YY_IE(i + 1, j)

                XVector(ic) = (XSW + XSE + XNE + XNW) / 4.
                YVector(ic) = (YSW + YSE + YNE + YNW) / 4.
                ZVector(ic) = Aux2D(i, j)

                write(45,*) XVector(ic), YVector(ic), ZVector(ic) 

                ic = ic + 1

            endif idef2

            end do
            end do

            write(45,*) "<end_xyz>"
            close(45)

            !Construct XYZ points collection
            call New(Me%XYZPoints, XVector, YVector, ZVector)

            deallocate(XVector)
            deallocate(YVector)                       
            deallocate(ZVector)

            call UnGetHorizontalGrid(Me%ObjGridIn, XX_IE, STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR90'

            call UnGetHorizontalGrid(Me%ObjGridIn, YY_IE, STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR100'

            call UnGetHorizontalGrid(Me%ObjGridIn, DefineCellsMap, STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR110'

        endif

        call UnGetGridData(Me%ObjGridDataIn, Aux2D, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR30'

        call KillGridData(Me%ObjGridDataIn, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR40'

        if (.not. Me%BatimInMatch) then
            call KillHorizontalGrid(Me%ObjGridIn, STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_)stop 'ConstructInitialGridData - DigitalTerrainTool - ERR50'
        endif

        

    end subroutine ConstructInitialGridData

    !--------------------------------------------------------------------------


    subroutine AllocateTempVariables

        !Begin-----------------------------------------------------------------
        
        allocate(Me%Depth    (Me%ExtVar%Size%ILB:Me%ExtVar%Size%IUB, &
                              Me%ExtVar%Size%JLB:Me%ExtVar%Size%JUB))
        
        allocate(Me%GridPoint(Me%ExtVar%Size%ILB:Me%ExtVar%Size%IUB, &
                              Me%ExtVar%Size%JLB:Me%ExtVar%Size%JUB))

        allocate(Me%AuxPoint)
        allocate(Me%Rect)
        Me%Rect%Count = 5

        allocate(Me%Rect%VerticesF(1:Me%Rect%Count))

i1:     if (Me%River%CanonicSpace) then

i2:         if      (Me%River%MainAxe == AlongLine  ) then

                allocate(Me%River%SectionsLocation(Me%ExtVar%WorkSize%JUB - Me%ExtVar%WorkSize%JLB + 1, 3))
                allocate(Me%River%BedDepth        (Me%ExtVar%WorkSize%JUB - Me%ExtVar%WorkSize%JLB + 1   ))

                Me%River%SectionsLocation(:,:) = FillValueInt
                Me%River%BedDepth        (:  ) = FillValueReal

            else if (Me%River%MainAxe == AlongColumn) then i2

                allocate(Me%River%SectionsLocation(Me%ExtVar%WorkSize%IUB - Me%ExtVar%WorkSize%ILB + 1, 3))
                allocate(Me%River%BedDepth        (Me%ExtVar%WorkSize%IUB - Me%ExtVar%WorkSize%ILB + 1   ))

                Me%River%SectionsLocation(:,:) = FillValueInt
                Me%River%BedDepth        (:  ) = FillValueReal

            else  i2

                stop 'AllocateTempVariables - DigitalTerrainCreator - ERR60'

            endif i2

        endif i1

   
    end subroutine AllocateTempVariables

    !--------------------------------------------------------------------------

    subroutine DeallocateTempVariables
        
        !Begin-----------------------------------------------------------------
        
        deallocate(Me%GridPoint     )
        deallocate(Me%AuxPoint      )
        deallocate(Me%Rect%VerticesF)
        deallocate(Me%Rect          )

            
        if (Me%River%CanonicSpace) then    
            deallocate(Me%River%SectionsLocation)
            deallocate(Me%River%BedDepth        )
        endif


   
    end subroutine DeallocateTempVariables

    !--------------------------------------------------------------------------


    subroutine Read_Lock_External_Var
        
        !External-----------------------------------------------------------------
        integer                             :: STAT_CALL
        integer                             :: GEOG, UTM, MIL_PORT, SIMPLE_GEOG
        integer                             :: GRID_COORD, CoordType, NLRD

        !Begin-----------------------------------------------------------------

        call GetGridOrigin    (Me%ObjGrid, Me%ExtVar%OriginX, Me%ExtVar%OriginY, STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - DigitalTerrainTool - ERR00'

        call GetGridAngle     (Me%ObjGrid, Me%ExtVar%Rotation, STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - DigitalTerrainTool - ERR01'

        call GetCheckDistortion(Me%ObjGrid, Me%ExtVar%GridDistortion, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - DigitalTerrainTool - ERR02'

        call GetHorizontalGridSize(Me%ObjGrid, WorkSize = Me%ExtVar%WorkSize, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - DigitalTerrainTool - ERR10'

        call GetHorizontalGridSize(Me%ObjGrid, Size = Me%ExtVar%Size, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - DigitalTerrainTool - ERR15'


        call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,             &
                               SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD,      &
                               NLRD = NLRD)

        !Gets Coordinates in use
        call GetGridCoordType(Me%ObjGrid, CoordType, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - DigitalTerrainTool - ERR11'
        

        if    (CoordType == SIMPLE_GEOG)then
            
            call GetGridLatitudeLongitude(Me%ObjGrid,                           &
                                          GridLatitudeConn  = Me%ExtVar%YY_IE,  &
                                          GridLongitudeConn = Me%ExtVar%XX_IE,  &
                                          STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - DigitalTerrainTool - ERR12'

        elseif(CoordType == UTM .or. CoordType == MIL_PORT .or.                 &
               CoordType == GRID_COORD .or. CoordType == NLRD)then

            call GetHorizontalGrid(Me%ObjGrid,                                  &
                                   XX_IE = Me%ExtVar%XX_IE,                     &
                                   YY_IE = Me%ExtVar%YY_IE,                     &
                                   STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - DigitalTerrainTool - ERR13'

        else

            write(*,*)'GEOG coordinate type cannot be used in digital terrain generation'
            stop 'Read_Lock_External_Var - DigitalTerrainTool - ERR13'

        end if

        call GetHorizontalGrid(Me%ObjGrid, XX = Me%ExtVar%XX, YY = Me%ExtVar%YY,  STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - DigitalTerrainTool - ERR40'

        call GetDefineCellsMap(Me%ObjGrid, Me%ExtVar%DefineCellsMap, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - DigitalTerrainTool - ERR50'

   
    end subroutine Read_Lock_External_Var

    !--------------------------------------------------------------------------



    subroutine Read_UnLock_External_Var

        !External-----------------------------------------------------------------
        integer                             :: STAT_CALL

        !Begin-----------------------------------------------------------------


        call UngetHorizontalGrid(Me%ObjGrid, Me%ExtVar%XX_IE, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_UnLock_External_Var - DigitalTerrainTool - ERR10'

        call UngetHorizontalGrid(Me%ObjGrid, Me%ExtVar%YY_IE, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_UnLock_External_Var - DigitalTerrainTool - ERR20'

        call UngetHorizontalGrid(Me%ObjGrid, Me%ExtVar%XX, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_UnLock_External_Var - DigitalTerrainTool - ERR30'

        call UngetHorizontalGrid(Me%ObjGrid, Me%ExtVar%YY, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_UnLock_External_Var - DigitalTerrainTool - ERR40'

        call UngetHorizontalGrid(Me%ObjGrid, Me%ExtVar%DefineCellsMap, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_UnLock_External_Var - DigitalTerrainTool - ERR50'

    end subroutine Read_UnLock_External_Var

    
    !--------------------------------------------------------------------------
    
    
    subroutine DefineGridPoints

        !Local-----------------------------------------------------------------
        integer                             :: i, j
        real                                :: XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW

        !Begin-----------------------------------------------------------------

        write(*,*)"Defining grid points..."

        do i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB
        do j = Me%ExtVar%WorkSize%JLB , Me%ExtVar%WorkSize%JUB
                       
            XSW = Me%ExtVar%XX_IE(i, j)
            YSW = Me%ExtVar%YY_IE(i, j)
            XSE = Me%ExtVar%XX_IE(i, j + 1)
            YSE = Me%ExtVar%YY_IE(i, j + 1)
            XNE = Me%ExtVar%XX_IE(i + 1, j + 1)
            YNE = Me%ExtVar%YY_IE(i + 1, j + 1)
            XNW = Me%ExtVar%XX_IE(i + 1, j)
            YNW = Me%ExtVar%YY_IE(i + 1, j)

            Me%GridPoint(i,j)%X = (XSW+XNW+XNE+XSE) / 4.
            Me%GridPoint(i,j)%Y = (YSW+YNW+YNE+YSE) / 4.


        end do
        end do

    end subroutine DefineGridPoints

    
    !--------------------------------------------------------------------------


    subroutine FillingCells

        !Local-----------------------------------------------------------------
        type (T_PointF),   pointer          :: GridPoint
        real                                :: SumOfDepths, MinimumDepth, MaximumDepth
        integer                             :: nPointsInside, CurrentPoint
        integer                             :: i, j
        real                                :: XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW

        !Begin-----------------------------------------------------------------

        write(*,*)"Filling cells with known data..."
        
        do i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB
        do j = Me%ExtVar%WorkSize%JLB , Me%ExtVar%WorkSize%JUB
                       
idef:       if (Me%ExtVar%DefineCellsMap(i, j)==1 .and. Me%Depth(i, j) == Me%NoDataPoint) then

                GridPoint => Me%GridPoint(i, j)
                
                if(IsVisible(Me%LandArea, GridPoint))then

                    Me%Depth(i,j) = Me%LandPoint

                else

                    XSW = Me%ExtVar%XX_IE(i, j)
                    YSW = Me%ExtVar%YY_IE(i, j)
                    XSE = Me%ExtVar%XX_IE(i, j + 1)
                    YSE = Me%ExtVar%YY_IE(i, j + 1)
                    XNE = Me%ExtVar%XX_IE(i + 1, j + 1)
                    YNE = Me%ExtVar%YY_IE(i + 1, j + 1)
                    XNW = Me%ExtVar%XX_IE(i + 1, j)
                    YNW = Me%ExtVar%YY_IE(i + 1, j)

                    Me%Rect%VerticesF(1)%X = XSW
                    Me%Rect%VerticesF(1)%Y = YSW
                    Me%Rect%VerticesF(2)%X = XSE
                    Me%Rect%VerticesF(2)%Y = YSE
                    Me%Rect%VerticesF(3)%X = XNE
                    Me%Rect%VerticesF(3)%Y = YNE
                    Me%Rect%VerticesF(4)%X = XNW
                    Me%Rect%VerticesF(4)%Y = YNW

                    Me%Rect%VerticesF(5)%X = Me%Rect%VerticesF(1)%X
                    Me%Rect%VerticesF(5)%Y = Me%Rect%VerticesF(1)%Y

                    call SetLimits(Me%Rect)

                    SumOfDepths   = 0.
                    nPointsInside = 0
                    MinimumDepth  = - FillValueReal 
                    MaximumDepth  = FillValueReal                   

                    !loop through all xyz points
                    Me%CurrentXYZPoints => Me%XYZPoints
    
                    !write(*,*)'Analysing i, j', i, j

                    do while(associated(Me%CurrentXYZPoints))

                        do CurrentPoint = 1, Me%CurrentXYZPoints%Count

                            if(.not. Me%CurrentXYZPoints%Inside(CurrentPoint))then
                                
                                Me%AuxPoint%X = Me%CurrentXYZPoints%X(CurrentPoint)
                                Me%AuxPoint%Y = Me%CurrentXYZPoints%Y(CurrentPoint)
    
                                if(IsPointInsidePolygon(Me%AuxPoint, Me%Rect))then
                                    !avoid computation with points not defined
                                    if (Me%CurrentXYZPoints%Z(CurrentPoint) /= Me%InvalidZValue) then
                                        SumOfDepths   = SumOfDepths + Me%CurrentXYZPoints%Z(CurrentPoint)
                                        if (Me%CurrentXYZPoints%Z(CurrentPoint) < MinimumDepth) then
                                            MinimumDepth = Me%CurrentXYZPoints%Z(CurrentPoint)
                                        endif
                                        if (Me%CurrentXYZPoints%Z(CurrentPoint) > MaximumDepth) then
                                            MaximumDepth = Me%CurrentXYZPoints%Z(CurrentPoint)
                                        endif                                        
                                        nPointsInside = nPointsInside + 1
                                        Me%CurrentXYZPoints%Inside(CurrentPoint) = .true.
                                    endif
                                end if

                            end if

                        end do

                        Me%CurrentXYZPoints => Me%CurrentXYZPoints%Next
                    end do

                    if(nPointsInside .eq. 0)then
                        !Assume the default vaule Me%NoDataPoint or BATIM_INI
                        !Me%Depth(i, j) = Me%NoDataPoint
                    else
                        if      (Me%FillMethod == Average_) then
                            Me%Depth(i, j) = SumOfDepths / nPointsInside
                        elseif  (Me%FillMethod == Minimum_) then
                            Me%Depth(i, j) = MinimumDepth
                        elseif   (Me%FillMethod == Maximum_) then
                            Me%Depth(i, j) = MaximumDepth
                        endif
                    end if

                end if

            else idef

                Me%Depth(i, j) = Me%LandPoint

            endif idef     

            nullify(GridPoint)
        end do
        end do

    end subroutine FillingCells


    !--------------------------------------------------------------------------
    
    
    subroutine FillingLand
    
        !Local-----------------------------------------------------------------
        type (T_PointF),   pointer          :: GridPoint
        integer                             :: i, j

        !Begin-----------------------------------------------------------------    
    
        write(*,*)"Filling Land ..."
        
        do i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB
        do j = Me%ExtVar%WorkSize%JLB , Me%ExtVar%WorkSize%JUB
                       
idef:       if (Me%ExtVar%DefineCellsMap(i, j)==1) then

                GridPoint => Me%GridPoint(i, j)
                
                if(IsVisible(Me%LandArea, GridPoint))then

                    Me%Depth(i,j) = Me%LandPoint
                
                endif

            else idef

                Me%Depth(i, j) = Me%LandPoint

            endif idef     

            nullify(GridPoint)
        end do
        end do

    end subroutine FillingLand

    !--------------------------------------------------------------------------
    
    subroutine FillingRegularGrid

        !Local-----------------------------------------------------------------
        type (T_PointF),         pointer    :: GridPoint
        real   , dimension(:,:), pointer    :: SumOfDepths  
        real   , dimension(:,:), pointer    :: MinimumDepth
        real   , dimension(:,:), pointer    :: MaximumDepth
        integer, dimension(:,:), pointer    :: nPointsInside
        integer                             :: i, j, CurrentPoint
        integer                             :: ILB, IUB, JLB, JUB

        !Begin-----------------------------------------------------------------


        write(*,*)"Filling cells with known data..."

        ILB = Me%ExtVar%WorkSize%ILB
        IUB = Me%ExtVar%WorkSize%IUB
        JLB = Me%ExtVar%WorkSize%JLB
        JUB = Me%ExtVar%WorkSize%JUB

        allocate(SumOfDepths  (ILB:IUB, JLB:JUB))
        allocate(MinimumDepth (ILB:IUB, JLB:JUB))
        allocate(MaximumDepth (ILB:IUB, JLB:JUB))
        allocate(nPointsInside(ILB:IUB, JLB:JUB))

        SumOfDepths  (:, :) = 0.
        nPointsInside(:, :) = 0
        
        MinimumDepth (:, :) = - FillValueReal
        MaximumDepth (:, :) = FillValueReal

        !loop through all xyz points
        Me%CurrentXYZPoints => Me%XYZPoints

        !write(*,*)'Analysing i, j', i, j

        do while(associated(Me%CurrentXYZPoints))

            do CurrentPoint = 1, Me%CurrentXYZPoints%Count
                
                if(.not. Me%CurrentXYZPoints%Inside(CurrentPoint))then

                    Me%AuxPoint%X = Me%CurrentXYZPoints%X(CurrentPoint)
                    Me%AuxPoint%Y = Me%CurrentXYZPoints%Y(CurrentPoint)

                    !Translate Point
                    call RODAXY(-Me%ExtVar%OriginX, -Me%ExtVar%OriginY, 0.0, Me%AuxPoint%X, Me%AuxPoint%Y)  

                    !Rotate Point
                    call RODAXY(0.0, 0.0, -Me%ExtVar%Rotation, Me%AuxPoint%X, Me%AuxPoint%Y)  

                    call LocateCell (Me%ExtVar%XX, Me%ExtVar%YY,           &
                                     Me%AuxPoint%X, Me%AuxPoint%Y,         &
                                     Me%ExtVar%WorkSize%ILB,               &
                                     Me%ExtVar%WorkSize%IUB+1,             &
                                     Me%ExtVar%WorkSize%JLB,               &
                                     Me%ExtVar%WorkSize%JUB+1,             &
                                     i, j)
                    
                    if (i .ne. null_int .and. j .ne. null_int)then
                        if (i <= IUB .and. i >= ILB .and. j <= JUB .and. j >= JLB) then
                            !avoid computation with points not defined
                            if (Me%CurrentXYZPoints%Z(CurrentPoint) /= Me%InvalidZValue) then
                                SumOfDepths  (i, j) = SumOfDepths  (i, j) + Me%CurrentXYZPoints%Z(CurrentPoint)
                                nPointsInside(i, j) = nPointsInside(i, j) + 1
                                if (Me%CurrentXYZPoints%Z(CurrentPoint) < MinimumDepth (i, j)) then
                                    MinimumDepth (i, j) = Me%CurrentXYZPoints%Z(CurrentPoint)
                                endif
                                if (Me%CurrentXYZPoints%Z(CurrentPoint) > MaximumDepth (i, j)) then
                                    MaximumDepth (i, j) = Me%CurrentXYZPoints%Z(CurrentPoint)
                                endif                                
                                Me%CurrentXYZPoints%Inside(CurrentPoint) = .true.
                            endif
                        endif
                    end if

                endif

            end do

            Me%CurrentXYZPoints => Me%CurrentXYZPoints%Next

        end do
        
        do j = JLB , JUB
        do i = ILB,  IUB

            if(nPointsInside(i, j) .eq. 0)then
                !Assume the default vaule Me%NoDataPoint or BATIM_INI
                !Me%Depth(i, j) = Me%NoDataPoint
            else
                if      (Me%FillMethod == Average_) then
                    Me%Depth(i, j) = SumOfDepths(i,j) / nPointsInside(i, j)
                elseif  (Me%FillMethod == Minimum_) then
                    Me%Depth(i, j) = MinimumDepth(i, j)
                elseif  (Me%FillMethod == Maximum_) then
                    Me%Depth(i, j) = MaximumDepth(i, j)
                endif
            end if
            

        end do
        end do

        deallocate(SumOfDepths  )
        deallocate(nPointsInside)
        deallocate(MinimumDepth )

        do j = JLB, JUB
        do i = ILB, IUB
                       
            if (Me%ExtVar%DefineCellsMap(i, j)==1) then

                GridPoint => Me%GridPoint(i, j)
                
                if(IsVisible(Me%LandArea, GridPoint))then

                    Me%Depth(i,j) = Me%LandPoint

                endif

            endif
        enddo
        enddo

    end subroutine FillingRegularGrid

    !--------------------------------------------------------------------------

    subroutine OverlappingGridData

        !Local-----------------------------------------------------------------
        type (T_PointF),      pointer       :: GridPoint
        real, dimension(:,:), pointer       :: OverlappingMatrix
        integer                             :: i, j, ig, STAT_CALL, ic
        logical                             :: AreaOK, DepthOK

        !Begin-----------------------------------------------------------------

        write(*,*)"Overlapping grid data files..."

        do ig = Me%Overlapping%Number, 1, -1

ift:        if (Me%Overlapping%DataInfo(ig)%InfoType == GridDataType) then

                call GetGridData(Me%Overlapping%DataInfo(ig)%ObjGridData, OverlappingMatrix, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'OverlappingGridData - DigitalTerrainTool - ERR10'


                do i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB
                do j = Me%ExtVar%WorkSize%JLB , Me%ExtVar%WorkSize%JUB

idef:               if (Me%ExtVar%DefineCellsMap(i, j)==1) then

                        GridPoint => Me%GridPoint(i, j)

                        if (associated(Me%Overlapping%DataInfo(ig)%Areas)) then
                        
                            AreaOK = IsVisible(Me%Overlapping%DataInfo(ig)%Areas, GridPoint) 

                        else 
                            !By default is assume the entrie domain
                            AreaOK = .true.

                        endif 

                        if (OverlappingMatrix(i, j) < Me%Overlapping%DataInfo(ig)%Hmax .and. &
                            OverlappingMatrix(i, j) > Me%Overlapping%DataInfo(ig)%Hmin) then

                            DepthOK = .true.

                        else

                            DepthOK = .false.

                        endif
                
                        if(AreaOK .and. DepthOK .and. OverlappingMatrix(i, j) /= Me%NoDataPoint) then

                            Me%Depth(i,j) = Me%Depth(i,j)           * (1. - Me%Overlapping%DataInfo(ig)%Percentage) + &
                                            OverlappingMatrix(i, j) *       Me%Overlapping%DataInfo(ig)%Percentage

                        end if
                        

                        nullify(GridPoint)

                    endif idef

                end do
                end do

                call UnGetGridData(Me%Overlapping%DataInfo(ig)%ObjGridData, OverlappingMatrix, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'OverlappingGridData - DigitalTerrainTool - ERR20'

            else if (Me%Overlapping%DataInfo(ig)%InfoType == CellsType) then ift

                
                do ic = 1, Me%Overlapping%DataInfo(ig)%CellsNumber

                    i             = Me%Overlapping%DataInfo(ig)%Cells(ic)%i
                    j             = Me%Overlapping%DataInfo(ig)%Cells(ic)%j
                    Me%Depth(i,j) = Me%Overlapping%DataInfo(ig)%Cells(ic)%Depth

                enddo

            endif ift

        enddo

    end subroutine OverlappingGridData

    
    !--------------------------------------------------------------------------


    subroutine SmoothGridData

        !Local-----------------------------------------------------------------
        type (T_PointF),   pointer          :: GridPoint
        real                                :: SumOfDepths
        integer                             :: nPointsInside, CurrentPoint
        integer                             :: i, j

        !Begin-----------------------------------------------------------------

        write(*,*)"Filling unknown data using known values inside a circle..."
        
        do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB
        do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB
        
            SumOfDepths   = 0.
            nPointsInside = 0

            if(Me%Depth(i, j) /= Me%LandPoint) then

                Me%CurrentXYZPoints => Me%XYZPoints

                do while(associated(Me%CurrentXYZPoints))

                    do CurrentPoint = 1, Me%CurrentXYZPoints%Count
                        
                        Me%AuxPoint%X = Me%CurrentXYZPoints%X(CurrentPoint)
                        Me%AuxPoint%Y = Me%CurrentXYZPoints%Y(CurrentPoint)

                        GridPoint => Me%GridPoint(i, j)

                        if(IsPointInsideCircle(Me%AuxPoint, GridPoint, Me%Radius))then
                            SumOfDepths   = SumOfDepths + Me%CurrentXYZPoints%Z(CurrentPoint)
                            nPointsInside = nPointsInside + 1
                        end if

                        nullify(GridPoint)

                    end do

                    Me%CurrentXYZPoints => Me%CurrentXYZPoints%Next
                end do

                if(nPointsInside .eq. 0)then
                    !Assume the default vaule Me%NoDataPoint or BATIM_INI
                    !Me%Depth(i, j) = Me%NoDataPoint
                else
                    Me%Depth(i, j) = SumOfDepths / nPointsInside
                end if
            
            end if

        enddo
        enddo

    end subroutine SmoothGridData

    
    !--------------------------------------------------------------------------


    subroutine SelectPointsWithData
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer       :: AuxVector
        integer                             :: Count, CellsNumber

        !Begin-----------------------------------------------------------------

        write(*,*)'Selecting points for interpolation...'

        CellsNumber =   (Me%ExtVar%WorkSize%IUB - Me%ExtVar%WorkSize%ILB + 1) *         &
                        (Me%ExtVar%WorkSize%JUB - Me%ExtVar%WorkSize%JLB + 1)

        allocate (AuxVector(1:(Me%TotalXYZ + CellsNumber),1:3)) 
        Count = 0

        select case(Me%PointsForInterpolation)

            case(AllPoints)

                call SelectAllOutsidePoints    (AuxVector, Count)
                call SelectGridCellCenterPoints(AuxVector, Count)


            case(GridCellCenterPoints)

                call SelectGridCellCenterPoints(AuxVector, Count)

            case(PointsInsideGridLimits)

                call SelectGridCellCenterPoints(AuxVector, Count)
                call SelectInsideGridPoints    (AuxVector, Count)
                
            case(OnlyXYZPoints) 
                call SelectAllXYZPoints        (AuxVector, Count)
            
        end select

        Me%TotalWithData = Count
                
        allocate (Me%PointsWithData(1:Count,1:3))
        Me%PointsWithData(1:Count,1:3) = AuxVector(1:Count,1:3)

        deallocate (AuxVector)

        if (Me%PointsForInterpolation == GridCellCenterPoints .and. Me%River%CanonicSpace) then

            call CanonicVersusCartesian (PointsGroup = Me%PointsWithData, Count = Count, ToCanonic = .true.)

            if (Me%River%MaintainBed) then
                call FindKnownRiverSections
            endif


        endif


    end subroutine SelectPointsWithData

    !--------------------------------------------------------------------------
    
    
    subroutine SelectInsideGridPoints(AuxVector, Count)
        
        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer       :: AuxVector

        !Local-----------------------------------------------------------------
        integer                             :: CurrentPoint
        integer                             :: Count

        !Begin-----------------------------------------------------------------

        Me%CurrentXYZPoints => Me%XYZPoints

        do while(associated(Me%CurrentXYZPoints))

            do CurrentPoint = 1, Me%CurrentXYZPoints%Count

                if(.not. Me%CurrentXYZPoints%Inside(CurrentPoint))then

                    if(Me%CurrentXYZPoints%X(CurrentPoint) .ge. Me%GridLimits%Left      .and. &
                       Me%CurrentXYZPoints%X(CurrentPoint) .le. Me%GridLimits%Right     .and. &
                       Me%CurrentXYZPoints%Y(CurrentPoint) .ge. Me%GridLimits%Bottom    .and. &
                       Me%CurrentXYZPoints%Y(CurrentPoint) .le. Me%GridLimits%Top)then
                        
                        Count = Count + 1

                        AuxVector(Count,1) = Me%CurrentXYZPoints%X(CurrentPoint)
                        AuxVector(Count,2) = Me%CurrentXYZPoints%Y(CurrentPoint)
                        AuxVector(Count,3) = Me%CurrentXYZPoints%Z(CurrentPoint)

                    end if

                end if

            end do

            Me%CurrentXYZPoints => Me%CurrentXYZPoints%Next

        end do
    
    
    end subroutine SelectInsideGridPoints
    

    !--------------------------------------------------------------------------
    
    
    subroutine SelectAllOutsidePoints(AuxVector, Count)
        
        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer       :: AuxVector

        !Local-----------------------------------------------------------------
        integer                             :: CurrentPoint
        integer                             :: Count

        !Begin-----------------------------------------------------------------

        Me%CurrentXYZPoints => Me%XYZPoints

        do while(associated(Me%CurrentXYZPoints))

            do CurrentPoint = 1, Me%CurrentXYZPoints%Count

                if(.not. Me%CurrentXYZPoints%Inside(CurrentPoint))then
                    Count = Count + 1
                    AuxVector(Count,1) = Me%CurrentXYZPoints%X(CurrentPoint)
                    AuxVector(Count,2) = Me%CurrentXYZPoints%Y(CurrentPoint)
                    AuxVector(Count,3) = Me%CurrentXYZPoints%Z(CurrentPoint)

                end if

            end do

            Me%CurrentXYZPoints => Me%CurrentXYZPoints%Next

        end do
    
    
    end subroutine SelectAllOutsidePoints
    

    !--------------------------------------------------------------------------
    
    
    !--------------------------------------------------------------------------
    
    
    subroutine SelectAllXYZPoints(AuxVector, Count)
        
        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer       :: AuxVector

        !Local-----------------------------------------------------------------
        integer                             :: CurrentPoint
        integer                             :: Count

        !Begin-----------------------------------------------------------------

        Me%CurrentXYZPoints => Me%XYZPoints

        do while(associated(Me%CurrentXYZPoints))

            do CurrentPoint = 1, Me%CurrentXYZPoints%Count

                Count = Count + 1
                AuxVector(Count,1) = Me%CurrentXYZPoints%X(CurrentPoint)
                AuxVector(Count,2) = Me%CurrentXYZPoints%Y(CurrentPoint)
                AuxVector(Count,3) = Me%CurrentXYZPoints%Z(CurrentPoint)

            end do

            Me%CurrentXYZPoints => Me%CurrentXYZPoints%Next

        end do
    
    
    end subroutine SelectAllXYZPoints
    

    !--------------------------------------------------------------------------
    subroutine SelectGridCellCenterPoints(AuxVector, Count)
        
        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer       :: AuxVector

        !Local-----------------------------------------------------------------
        integer                             :: i, j
        integer                             :: Count

        !Begin-----------------------------------------------------------------

        do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB
        do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB

            if(Me%Depth(i, j) > Me%LandPoint / 2.) then
           
                Count = Count + 1

                AuxVector(Count,1) =  Me%GridPoint(i,j)%X
                AuxVector(Count,2) =  Me%GridPoint(i,j)%Y
                AuxVector(Count,3) =  Me%Depth(i, j)

            elseif(Me%Depth(i, j) == Me%LandPoint .and. Me%ExtVar%DefineCellsMap(i, j) == 1)then

                if(Me%AssumeLandDepth)then

                    if(IsNearWaterOrPointWithNoValue(i, j))then

                        Count = Count + 1

                        AuxVector(Count,1) =  Me%GridPoint(i,j)%X
                        AuxVector(Count,2) =  Me%GridPoint(i,j)%Y
                        AuxVector(Count,3) =  Me%LandDepth

                    end if

                end if

            end if

        enddo
        enddo

    
    
    end subroutine SelectGridCellCenterPoints
    
    !--------------------------------------------------------------------------

    subroutine SelectPointsWithNoData
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer   :: AuxVector
        integer                             :: i, j, Count, GridPoints

        !Begin-----------------------------------------------------------------

        write(*,*)'Selecting unfilled points...'

        GridPoints = (Me%ExtVar%WorkSize%IUB - Me%ExtVar%WorkSize%ILB + 1) * &
                     (Me%ExtVar%WorkSize%JUB - Me%ExtVar%WorkSize%JLB + 1)
        

        allocate (AuxVector(1:GridPoints,1:3))
        Count = 0

        do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB
        do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB

            if(Me%Depth(i, j) == Me%NoDataPoint) then
                   
                Count = Count + 1

                AuxVector(Count,1) =  Me%GridPoint(i,j)%X
                AuxVector(Count,2) =  Me%GridPoint(i,j)%Y
                AuxVector(Count,3) =  Me%NoDataPoint

            end if

        enddo
        enddo

        Me%TotalNoData = Count
                
        allocate (Me%PointsNoData(1:Count,1:3)) 

        Me%PointsNoData(1:Count,1:3) = AuxVector(1:Count,1:3)

        deallocate (AuxVector)


        if (Me%PointsForInterpolation == GridCellCenterPoints .and. Me%River%CanonicSpace) then

            call CanonicVersusCartesian (PointsGroup = Me%PointsNoData, Count = Count, ToCanonic = .true.)

        endif



    end subroutine SelectPointsWithNoData


    !--------------------------------------------------------------------------


    subroutine SetGridLimits
        
        !Local-----------------------------------------------------------------
        !type (T_PointF), pointer            :: LowerLeft, UpperLeft
        !type (T_PointF), pointer            :: LowerRight, UpperRight
        !integer                             :: i, j
        real                                :: DX, DY
        !real                                :: OriginX, OriginY, Rotation
        integer                             :: STAT_CALL

        !Begin-----------------------------------------------------------------

        write(*,*)"Setting grid limits..."
        
        call GetGridBorderLimits(HorizontalGridID   = Me%ObjGrid,                       &
                                 West               = Me%GridLimits%Left  ,             &
                                 East               = Me%GridLimits%Right ,             & 
                                 South              = Me%GridLimits%Bottom,             &
                                 North              = Me%GridLimits%Top   ,             &
                                 STAT               = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverlappingGridData - DigitalTerrainTool - ERR20'        

        if(Me%ExpandGridLimits)then

            DX = (Me%GridLimits%Right - Me%GridLimits%Left  ) * Me%ExpandGridLimitsFraction
            DY = (Me%GridLimits%Top   - Me%GridLimits%Bottom) * Me%ExpandGridLimitsFraction
            
            Me%GridLimits%Left   = Me%GridLimits%Left   - DX 
            Me%GridLimits%Right  = Me%GridLimits%Right  + DX 
            Me%GridLimits%Bottom = Me%GridLimits%Bottom - DY
            Me%GridLimits%Top    = Me%GridLimits%Top    + DY

        endif
        

    end subroutine SetGridLimits
    
    !---------------------------------------------------------------------------

    subroutine Triangulator

        !Local------------------------------------------------------------------
        real, dimension(:), pointer         :: NodeX, NodeY, NodeZ
        integer                             :: i, j, STAT_CALL, Count, nNodes
        integer                             :: UnitNumber, iT, nTriangles
        real,    dimension(:), pointer      :: XT, YT, ZT
        integer, dimension(:), pointer      :: V1, V2, V3
        real,    dimension(:), allocatable  :: Aux2

        !Begin------------------------------------------------------------------


        write(*,*)"Performing triangulation..."

        NodeX =>  Me%PointsWithData(:,1)
        NodeY =>  Me%PointsWithData(:,2)
        NodeZ =>  Me%PointsWithData(:,3)
        
        NodeX = NodeX * Me%TriangScale 
        NodeY = NodeY * Me%TriangScale

        !Constructs Triangulation
        call ConstructTriangulation (Me%ObjTriangulation, Me%TotalWithData, NodeX, NodeY, NodeZ,   &
                                     Me%Triang%Tolerance, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Triangulator - Digital Terrain Creator - ERR10'

        Count = 1

        allocate(Aux2(1:2))

        do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB
        do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB

            if(Me%Depth(i, j) == Me%NoDataPoint) then

                Aux2(1) = Me%GridPoint(i,j)%X
                Aux2(2) = Me%GridPoint(i,j)%Y

                if (Me%PointsForInterpolation == GridCellCenterPoints .and. Me%River%CanonicSpace) then

                    call CanonicVersusCartesian (Point = Aux2, ToCanonic = .true.)

                endif
                
                Aux2(1) = Aux2(1) * Me%TriangScale
                Aux2(2) = Aux2(2) * Me%TriangScale

                Me%Depth(i, j) = InterPolation(Me%ObjTriangulation,         &
                                               Aux2(1), Aux2(2),            &
                                               Me%Triang%FillOutsidePoints, &
                                               Default  = Me%NoDataPoint,   &
                                               STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Triangulator - Digital Terrain Creator - ERR20'

                Me%PointsNoData(Count,3) = Me%Depth(i, j)

                Count = Count + 1

            end if

        enddo
        enddo

        deallocate(Aux2)



        !Opens resulting triangle file
        if (Me%Triang%OutputTriangles) then

            write(*,*)"Writing triangles file..."

            !Get the number of triangles
            call GetNumberOfTriangles   (Me%ObjTriangulation, nTriangles)

            !Allocates space for the Triangle vertices and gets them
            allocate(V1(nTriangles))
            allocate(V2(nTriangles))
            allocate(V3(nTriangles))

            call GetTriangleList (Me%ObjTriangulation, v1, v2, v3, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Triangulator - Digital Terrain Creator - ERR21'


            !Gets nodes effictive used and the reordered nodes 
            call GetNumberOfNodes (Me%ObjTriangulation, nNodes, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Triangulator - Digital Terrain Creator - ERR30'

            allocate(XT(nNodes))
            allocate(YT(nNodes))
            allocate(ZT(nNodes))

            call GetNodesList   (Me%ObjTriangulation, XT, YT, ZT, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunTriangulation - ERR14b'   
            
            XT = XT / Me%TriangScale
            YT = YT / Me%TriangScale


            call UnitsManager (UnitNumber, OPEN_FILE, STAT = STAT_CALL)
            open (unit=UnitNumber, status = 'unknown', file = trim(adjustl(Me%Triang%File)))
            do iT = 1, nTriangles
                write(UnitNumber,*)'<beginpolygon>'
                write(UnitNumber,*)XT(V1(iT)), YT(V1(iT))
                write(UnitNumber,*)XT(V2(iT)), YT(V2(iT))
                write(UnitNumber,*)XT(V3(iT)), YT(V3(iT))
                write(UnitNumber,*)XT(V1(iT)), YT(V1(iT))
                write(UnitNumber,*)'<endpolygon>'
            enddo
            call UnitsManager (UnitNumber, CLOSE_FILE, STAT = STAT_CALL)

            deallocate(XT, YT, ZT, v1, v2, v3)
            
        endif

        call KillTriangulation (Me%ObjTriangulation, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Triangulator - Digital Terrain Creator - ERR40'

    end subroutine Triangulator

    !--------------------------------------------------------------------------
    
    !---------------------------------------------------------------------------

    subroutine Triangulator_DD

        !Local------------------------------------------------------------------
        real, dimension(:), pointer         :: NodeX, NodeY, NodeZ
        real, dimension(:), pointer         :: Aux_NodeX, Aux_NodeY, Aux_NodeZ        
        integer                             :: i, j, STAT_CALL, icount, nNodes
        integer                             :: n, nD, DomainNodes
        integer                             :: ILB, IUB, JLB, JUB
        real                                :: DX, DY, Left, Right, Bottom, Top, DXmin, DYmin 
        real,    dimension(:), allocatable  :: Aux2
        logical                             :: LandArea, InterpolOk, DomainOk

        !Begin------------------------------------------------------------------

        nNodes = size(Array = Me%PointsWithData, Dim = 1) 
            
        allocate (Aux_NodeX(nNodes), Aux_NodeY(nNodes), Aux_NodeZ(nNodes))
        
        DXmin = (Me%GridLimits%Right - Me%GridLimits%Left  ) * Me%ExpandGridLimitsFraction / Me%DD%nJX 
        DYmin = (Me%GridLimits%Top   - Me%GridLimits%Bottom) * Me%ExpandGridLimitsFraction / Me%DD%nIY       
        
        if (associated(Me%LandArea)) then
            LandArea = .true.
        else
            LandArea = .false.
        endif
        
d1:     do n = 1, Me%DD%nD
    
            InterpolOk = .true.
            DomainOk   = .true.
    
            if (.not. associated(Me%DD%Border(n)%Polygon_)) then
                DomainOk    = .false.
                InterpolOk  = .false.  
            endif
            
            if (DomainOk.and. LandArea) then
    
                if (VertPolygonInsidePolygon(PolygonA = Me%DD%Border(n)%Polygon_, &
                                             PolygonB = Me%LandArea)) then
                    InterpolOk  = .false.    
                endif
                                             
            endif
            
            if (InterpolOk) then
    
                write(*,*)"Performing triangulation of Domain =", n

                NodeX =>  Me%PointsWithData(:,1)
                NodeY =>  Me%PointsWithData(:,2)
                NodeZ =>  Me%PointsWithData(:,3)
            
                Left   = Me%DD%Border(n)%Polygon_%Limits%Left  
                Right  = Me%DD%Border(n)%Polygon_%Limits%Right 
                Bottom = Me%DD%Border(n)%Polygon_%Limits%Bottom
                Top    = Me%DD%Border(n)%Polygon_%Limits%Top               

                DX = (Right - Left  ) * Me%ExpandGridLimitsFraction
                DY = (Top   - Bottom) * Me%ExpandGridLimitsFraction
                
                DX = max (DX, DXmin)
                DY = max (DY, DYmin)

                Left   = Left   - DX 
                Right  = Right  + DX 
                Bottom = Bottom - DY
                Top    = Top    + DY
            
                icount = 0
                do nD = 1, nNodes

                    if (NodeX(nD) > Left    .and. NodeX(nD) < Right  .and.          &
                        NodeY(nD) > Bottom  .and. NodeY(nD) < Top) then
                    
                        icount            = icount + 1
                        Aux_NodeX(icount) = NodeX(nD)
                        Aux_NodeY(icount) = NodeY(nD)
                        Aux_NodeZ(icount) = NodeZ(nD)
                    
                    endif
                
                enddo
            
                nullify(NodeX, NodeY, NodeZ)
                DomainNodes = icount
        
                allocate (NodeX(DomainNodes), NodeY(DomainNodes), NodeZ(DomainNodes))
            
                NodeX(1:DomainNodes) = Aux_NodeX(1:DomainNodes)
                NodeY(1:DomainNodes) = Aux_NodeY(1:DomainNodes)
                NodeZ(1:DomainNodes) = Aux_NodeZ(1:DomainNodes)
            
                NodeX = NodeX * Me%TriangScale 
                NodeY = NodeY * Me%TriangScale
            
                Me%ObjTriangulation = 0
                
                if (DomainNodes >= 10) then
            
                    !Constructs Triangulation
                    call ConstructTriangulation (Me%ObjTriangulation, DomainNodes, NodeX, NodeY, NodeZ,   &
                                                    Me%Triang%Tolerance, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Triangulator_DD - Digital Terrain Creator - ERR10'
        
                    ILB = Me%DD%BorderLimits(n)%ILB
                    IUB = Me%DD%BorderLimits(n)%IUB 
                    JLB = Me%DD%BorderLimits(n)%JLB
                    JUB = Me%DD%BorderLimits(n)%JUB
            
                    allocate(Aux2(1:2))

                    do i = ILB, IUB
                    do j = JLB, JUB

                        if(Me%Depth(i, j) == Me%NoDataPoint) then

                            Aux2(1) = Me%GridPoint(i,j)%X
                            Aux2(2) = Me%GridPoint(i,j)%Y
                
                            Aux2(1) = Aux2(1) * Me%TriangScale
                            Aux2(2) = Aux2(2) * Me%TriangScale

                            Me%Depth(i, j) = InterPolation(Me%ObjTriangulation,         &
                                                           Aux2(1), Aux2(2),            &
                                                           Me%Triang%FillOutsidePoints, &
                                                           Default  = Me%NoDataPoint,   &
                                                           STAT     = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Triangulator_DD - Digital Terrain Creator - ERR20'

                        end if

                    enddo
                    enddo

                    deallocate(Aux2)

                    call KillTriangulation (Me%ObjTriangulation, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Triangulator_DD - Digital Terrain Creator - ERR40'
                    
                endif
            
                deallocate (NodeX, NodeY, NodeZ)
                
            endif
            
        enddo d1

        deallocate (Aux_NodeX, Aux_NodeY, Aux_NodeZ)

    end subroutine Triangulator_DD


    !--------------------------------------------------------------------------
    
    subroutine CanonicVersusCartesian (PointsGroup, Count, Point, ToCanonic)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer, optional :: PointsGroup
        real, dimension(:)           , optional :: Point
        integer                      , optional :: Count
        logical                                 :: ToCanonic

        !Local-----------------------------------------------------------------
        real                                    :: XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW
        integer                                 :: Count_, i, j, CurrentPoint
        logical                                 :: OnePoint

        !Begin-----------------------------------------------------------------

        if (present(Point)) then
            Count_   = 1
            OnePoint = .true.
        else
            Count_   = Count
            OnePoint = .false.
        endif


        do i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB
        do j = Me%ExtVar%WorkSize%JLB , Me%ExtVar%WorkSize%JUB
                       
idef:       if (Me%ExtVar%DefineCellsMap(i, j)==1) then

                if (ToCanonic) then

                    XSW = Me%ExtVar%XX_IE(i, j)
                    YSW = Me%ExtVar%YY_IE(i, j)
                    XSE = Me%ExtVar%XX_IE(i, j + 1)
                    YSE = Me%ExtVar%YY_IE(i, j + 1)
                    XNE = Me%ExtVar%XX_IE(i + 1, j + 1)
                    YNE = Me%ExtVar%YY_IE(i + 1, j + 1)
                    XNW = Me%ExtVar%XX_IE(i + 1, j)
                    YNW = Me%ExtVar%YY_IE(i + 1, j)

                else

                    XSW = real(j)
                    YSW = real(i)
                    XSE = real(j) + 1.
                    YSE = real(i)
                    XNE = real(j) + 1.
                    YNE = real(i) + 1.
                    XNW = real(j)
                    YNW = real(i) + 1.

                endif

                Me%Rect%VerticesF(1)%X = XSW
                Me%Rect%VerticesF(1)%Y = YSW
                Me%Rect%VerticesF(2)%X = XSE
                Me%Rect%VerticesF(2)%Y = YSE
                Me%Rect%VerticesF(3)%X = XNE
                Me%Rect%VerticesF(3)%Y = YNE
                Me%Rect%VerticesF(4)%X = XNW
                Me%Rect%VerticesF(4)%Y = YNW

                Me%Rect%VerticesF(5)%X = Me%Rect%VerticesF(1)%X
                Me%Rect%VerticesF(5)%Y = Me%Rect%VerticesF(1)%Y

                call SetLimits(Me%Rect)


                do CurrentPoint = 1, Count_
                
                    if (OnePoint) then
                        Me%AuxPoint%X = Point(1)
                        Me%AuxPoint%Y = Point(2)
                    else
                        Me%AuxPoint%X = PointsGroup(CurrentPoint, 1)
                        Me%AuxPoint%Y = PointsGroup(CurrentPoint, 2)
                    endif

                    if(IsPointInsidePolygon(Me%AuxPoint, Me%Rect))then

                        if (ToCanonic) then
                            if (OnePoint) then
                                Point(1) = real(J) + 0.5
                                Point(2) = real(I) + 0.5
                            else
                                PointsGroup(CurrentPoint, 1) = real(J) + 0.5
                                PointsGroup(CurrentPoint, 2) = real(I) + 0.5
                            endif

                        else

                            XSW = Me%ExtVar%XX_IE(i, j)
                            YSW = Me%ExtVar%YY_IE(i, j)
                            XSE = Me%ExtVar%XX_IE(i, j + 1)
                            YSE = Me%ExtVar%YY_IE(i, j + 1)
                            XNE = Me%ExtVar%XX_IE(i + 1, j + 1)
                            YNE = Me%ExtVar%YY_IE(i + 1, j + 1)
                            XNW = Me%ExtVar%XX_IE(i + 1, j)
                            YNW = Me%ExtVar%YY_IE(i + 1, j)

                            if (OnePoint) then
                                Point(1) = (XSW + XSE + XNE + XNW) / 4.
                                Point(2) = (YSW + YSE + YNE + YNW) / 4.
                            else
                                PointsGroup(CurrentPoint, 1) = (XSW + XSE + XNE + XNW) / 4.
                                PointsGroup(CurrentPoint, 2) = (YSW + YSE + YNE + YNW) / 4.
                            endif

                        endif
                        
                    end if
                end do

            endif idef     

        end do
        end do

    end subroutine CanonicVersusCartesian

    !--------------------------------------------------------------------------
    subroutine FindKnownRiverSections 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                 :: i, j, Counter
        logical                                 :: FoundSection, EndSection

        !Begin-----------------------------------------------------------------

        Counter                    = 0
        FoundSection               = .false.

i1:     if (Me%River%MainAxe == AlongColumn) then

                       
d1:         do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB

                !Find section begin              
i2:             if (.not. FoundSection) then
                    
d2:                 do  j = Me%ExtVar%WorkSize%JLB,  Me%ExtVar%WorkSize%JUB
                        
i3:                     if (Me%Depth(i,j) /= Me%NoDataPoint.and. Me%Depth(i,j) /= Me%LandPoint) then
                             Counter                              = Counter + 1
                             Me%River%BedDepth(Counter)           = FillValueReal
                             Me%River%SectionsLocation(Counter,1) = i
                             FoundSection                         = .true.
                             exit
                        endif i3

                    enddo d2

                endif i2

i4:             if (FoundSection) then

                    !Find max bed depth
d3:                 do  j = Me%ExtVar%WorkSize%JLB,  Me%ExtVar%WorkSize%JUB
                        if (Me%Depth(i,j) /= Me%NoDataPoint.and. Me%Depth(i,j) /= Me%LandPoint) then
                            !Me%River%BedDepth(Counter) = max(Me%River%BedDepth(Counter), Me%Depth(i,j))
                            if (Me%Depth(i,j)> Me%River%BedDepth(Counter)) then
                                Me%River%BedDepth        (Counter)   = Me%Depth(i,j)
                                Me%River%SectionsLocation(Counter,3) = i
                            endif
                        endif
                    enddo d3

                    EndSection = .true.

                    !Find section end
d4:                 do  j = Me%ExtVar%WorkSize%JLB,  Me%ExtVar%WorkSize%JUB

                        if (Me%Depth(i,j) /= Me%NoDataPoint.and. Me%Depth(i,j) /= Me%LandPoint) then
                             EndSection   = .false.
                             exit
                        endif

                    enddo d4

i5:                 if (EndSection) then
                        FoundSection                         = .false.
                        Me%River%SectionsLocation(Counter,2) = i
                    endif i5

                endif i4

            enddo d1

        else  if (Me%River%MainAxe == AlongLine) then i1 

d5:         do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB

                !Find section begin              
i6:             if (.not. FoundSection) then
                    
d6:                 do  i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB
                        
i7:                     if (Me%Depth(i,j) /= Me%NoDataPoint.and. Me%Depth(i,j) /= Me%LandPoint) then
                             Counter                              = Counter + 1
                             Me%River%BedDepth        (Counter)   = FillValueReal
                             Me%River%SectionsLocation(Counter,1) = j
                             FoundSection                         = .true.
                             exit                             
                        endif i7

                    enddo d6

                endif i6

i8:             if (FoundSection) then


                    !Find max bed depth
d7:                 do  i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB
                        if (Me%Depth(i,j) /= Me%NoDataPoint.and. Me%Depth(i,j) /= Me%LandPoint) then
                            !Me%River%BedDepth(Counter) = max(Me%River%BedDepth(Counter), Me%Depth(i,j))
                            if (Me%Depth(i,j)> Me%River%BedDepth(Counter)) then
                                Me%River%BedDepth        (Counter)   = Me%Depth(i,j)
                                Me%River%SectionsLocation(Counter,3) = j
                            endif
                        endif
                    enddo d7

                    EndSection = .true.

                    !Find section end
d8:                 do  i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB

                        if (Me%Depth(i,j) /= Me%NoDataPoint.and. Me%Depth(i,j) /= Me%LandPoint) then
                             EndSection   = .false.
                             exit
                        endif

                    enddo d8

i9:                 if (EndSection) then
                        FoundSection                         = .false.
                        Me%River%SectionsLocation(Counter,2) = j
                    endif i9

                endif i8

            enddo d5

        endif i1

        Me%River%SectionsNumber = Counter

    end subroutine FindKnownRiverSections

    !--------------------------------------------------------------------------
    subroutine DredgeRiverBed 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                 :: i, i_Previous, i_Next
        integer                                 :: j, j_Previous, j_Next
        integer                                 :: cs_Previous, cs_Next
        logical                                 :: InsideSection
        real                                    :: MaxDepth, MinDepth, DhMax, AverageMax, Hlimit

        !Begin-----------------------------------------------------------------


i1:     if (Me%River%MainAxe == AlongColumn) then

                                   
d1:         do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB
              
                call LocateSection(i, cs_Previous, cs_Next, InsideSection) 

i2:             if (.not. InsideSection) then

                    i_Previous =   Me%River%SectionsLocation(cs_Previous,2)
                    i_Next     =   Me%River%SectionsLocation(cs_Next,    1)

                    MaxDepth =   FillValueReal
                    MinDepth = - FillValueReal

d2:                 do  j = Me%ExtVar%WorkSize%JLB,  Me%ExtVar%WorkSize%JUB

                        if (Me%Depth(i,j) /= Me%NoDataPoint.and. Me%Depth(i,j) /= Me%LandPoint) then
                             MaxDepth = max(Me%Depth(i,j), MaxDepth)
                             MinDepth = min(Me%Depth(i,j), MinDepth)
                        endif

                    enddo d2

                    AverageMax = (Me%River%BedDepth(cs_Next    ) * real(i-i_Previous) +  &
                                  Me%River%BedDepth(cs_Previous) * real(i_Next-i    ))/  &
                                  real(i_Next - i_Previous)

                    DhMax = AverageMax - MaxDepth

                    Hlimit = MaxDepth - Me%River%FactorLimit * (MaxDepth - MinDepth)

i3:                 if (DhMax > 0.) then

d3:                     do  j = Me%ExtVar%WorkSize%JLB,  Me%ExtVar%WorkSize%JUB

i4:                         if (Me%Depth(i,j) /= Me%NoDataPoint.and. Me%Depth(i,j) /= Me%LandPoint) then
                            
i5:                             if (Me%Depth(i,j) > Hlimit) then

i10:                                if      (Me%River%Dredge == Linear) then
                                    
                                        Me%Depth(i,j) = Me%Depth(i,j) + DhMax * (Me%Depth(i,j) - Hlimit) /  (MaxDepth - Hlimit)

                                    else if (Me%River%Dredge == Uniform) then i10

                                        Me%Depth(i,j) = Me%Depth(i,j) + DhMax

                                    endif i10

                                   
                                endif i5

                            endif i4

                        enddo d3

                    endif i3


                endif i2

            enddo  d1

        else  if (Me%River%MainAxe == AlongLine) then i1

d4:         do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB
              
                call LocateSection(j, cs_Previous, cs_Next, InsideSection) 
                
i6:             if (.not. InsideSection) then

                    j_Previous =   Me%River%SectionsLocation(cs_Previous,2)
                    j_Next     =   Me%River%SectionsLocation(cs_Next,    1)

                    MaxDepth   =   FillValueReal
                    MinDepth   = - FillValueReal

d5:                 do  i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB

                        if (Me%Depth(i,j) /= Me%NoDataPoint.and. Me%Depth(i,j) /= Me%LandPoint) then
                             MaxDepth = max(Me%Depth(i,j), MaxDepth)
                             MinDepth = min(Me%Depth(i,j), MinDepth)
                        endif

                    enddo d5

                    AverageMax = (Me%River%BedDepth(cs_Next    ) * real(j-j_Previous) +  &
                                  Me%River%BedDepth(cs_Previous) * real(j_Next-j    ))/  &
                                  real(j_Next - j_Previous)

                    DhMax = AverageMax - MaxDepth

                    Hlimit = MaxDepth - Me%River%FactorLimit * (MaxDepth - MinDepth)

i7:                 if (DhMax > 0.) then

d6:                     do  i = Me%ExtVar%WorkSize%ILB,  Me%ExtVar%WorkSize%IUB

i8:                         if (Me%Depth(i,j) /= Me%NoDataPoint.and. Me%Depth(i,j) /= Me%LandPoint) then
                            
i9:                             if (Me%Depth(i,j) > Hlimit) then

i11:                                if      (Me%River%Dredge == Linear) then
                                    
                                        Me%Depth(i,j) = Me%Depth(i,j) + DhMax * (Me%Depth(i,j) - Hlimit) /  (MaxDepth - Hlimit)

                                    else if (Me%River%Dredge == Uniform) then i11

                                        Me%Depth(i,j) = Me%Depth(i,j) + DhMax

                                    endif i11
                                   
                                endif i9

                            endif i8

                        enddo d6

                    endif i7


                endif i6

            enddo  d4

        endif i1



    end subroutine DredgeRiverBed

    !--------------------------------------------------------------------------
    subroutine DoFilterDepths 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer           :: DepthAux
        integer                                 :: i, ii, iw
        integer                                 :: j, jj, jw
        integer                                 :: Counter, NValues, I50
        real,   dimension(:), pointer           :: AuxBat

        !Begin-----------------------------------------------------------------

        allocate(DepthAux(Me%ExtVar%Size%ILB : Me%ExtVar%Size%IUB,                      &
                          Me%ExtVar%Size%JLB : Me%ExtVar%Size%JUB))
        
        NValues = ((Me%FilterRadiusI + 1)*(Me%FilterRadiusJ + 1))**2 
            
        allocate (AuxBat (1:NValues))         

                                  
d1:     do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB
d2:     do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB

i4:         if (Me%Depth(i,j) /= Me%LandPoint) then
                
                Counter = 0

                AuxBat(:) = - FillValueReal                     
                    
                do jj = j - Me%FilterRadiusJ, j + Me%FilterRadiusJ
                do ii = i - Me%FilterRadiusI, i + Me%FilterRadiusI
                        
                    iw = max(ii, Me%ExtVar%WorkSize%ILB) 
                    jw = max(jj, Me%ExtVar%WorkSize%JLB)
                        
                    iw = min(iw, Me%ExtVar%WorkSize%IUB) 
                    jw = min(jw, Me%ExtVar%WorkSize%JUB)  
                        
                    if (Me%Depth(iw,jw) /= Me%NoDataPoint.and. Me%Depth(iw,jw) /= Me%LandPoint) then
                        Counter = Counter + 1
                        AuxBat(Counter) =  Me%Depth(iw, jw)
                    endif

                enddo
                enddo

                call Insertion_Sort(AuxBat)
                i50 = int(real(Counter)/2.)
                if (i50 > 0) then
                    DepthAux(i, j) = AuxBat(i50)
                endif
            else
                
                DepthAux(i, j) = Me%LandPoint
                
            endif i4

        enddo d2
        enddo d1


        Me%Depth(:,:) = DepthAux(:,:)

        deallocate(DepthAux)


    end subroutine DoFilterDepths

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    subroutine DoFilterSpikes 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer           :: NeighborsAverage
        integer                                 :: i, iw
        integer                                 :: j, jw
        real                                    :: Counter, AuxSum, AverageDif, CellDif 


        !Begin-----------------------------------------------------------------

        allocate(NeighborsAverage(Me%ExtVar%Size%ILB : Me%ExtVar%Size%IUB,              &
                                  Me%ExtVar%Size%JLB : Me%ExtVar%Size%JUB))

        NeighborsAverage(:,:) = FillValueReal

                                  
d1:     do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB
d2:     do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB

i4:         if (Me%Depth(i,j) /= Me%NoDataPoint.and. Me%Depth(i,j) /= Me%LandPoint) then

                Counter = 0
                AuxSum  = 0.
                do iw = i - 1, i + 1
                do jw = j - 1, j + 1
                    if (iw /= i .and. jw /= j              .and.                    &
                        Me%Depth(iw,jw)  /= Me%NoDataPoint .and.                    &
                        Me%Depth(iw,jw)  /= Me%LandPoint) then

                        Counter = Counter + 1
                        AuxSum  = AuxSum  + Me%Depth(iw, jw)

                    endif

                enddo
                enddo

                if (Counter > 0.) NeighborsAverage(i,j) = AuxSum / real(Counter)
                
        
            endif i4

        enddo d2
        enddo d1



d4:     do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB
d5:     do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB

            if (Me%Depth(i,j) /= Me%NoDataPoint.and. Me%Depth(i,j) /= Me%LandPoint) then


                Counter = 0
                AuxSum  = 0.
                do iw = i - 1, i + 1
                do jw = j - 1, j + 1

                    if (iw /= i .and. jw /= j              .and.                    &
                        Me%Depth(iw,jw)  /= Me%NoDataPoint .and.                    &
                        Me%Depth(iw,jw)  /= Me%LandPoint) then

                        Counter = Counter + 1
                        AuxSum  = AuxSum + abs(Me%Depth(iw, jw) - NeighborsAverage(i,j))

                    endif

                enddo
                enddo

                if (Counter > 0.) then
                    AverageDif  = AuxSum / real(Counter)
                    CellDif     = abs(Me%Depth(i, j) - NeighborsAverage(i,j))
                   
                    if (CellDif > AverageDif * Me%SpikesFactor) Me%Depth(i,j) = NeighborsAverage(i,j)
                    
                endif


            endif
        enddo d5
        enddo d4


        deallocate(NeighborsAverage)


    end subroutine DoFilterSpikes

    !--------------------------------------------------------------------------
    
    subroutine LocateSection(ij_Input, ij_Previous, ij_Next, InsideSection) 

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                     :: ij_Input
        integer, intent(OUT)                    :: ij_Previous, ij_Next
        logical, intent(OUT)                    :: InsideSection

        !Local-----------------------------------------------------------------
        integer                                 :: ij

        !Begin-----------------------------------------------------------------

            InsideSection = .false.                      

            ij_Previous   = -99
            ij_Next       = -99

i1:         if (Me%River%LocationMethod == DefinedDepths) then

                do ij = 1, Me%River%SectionsNumber

                    if (ij_Input >= Me%River%SectionsLocation(ij,1) .and. &
                        ij_Input <= Me%River%SectionsLocation(ij,2)) then
                        InsideSection = .true.
                        ij_Next     = ij
                        ij_Previous = ij
                        exit
                    endif


                    if (Me%River%SectionsLocation(ij,2) >  ij_Input) then
                        ij_Next     = ij
                        ij_Previous = ij - 1
                        exit
                    endif

                enddo

            else if (Me%River%LocationMethod == MaximumDepth) then i1

                do ij = 1, Me%River%SectionsNumber

                    if (ij_Input == Me%River%SectionsLocation(ij,3)) then
                        InsideSection = .true.
                        ij_Next     = ij
                        ij_Previous = ij
                        exit
                    endif

                    if (Me%River%SectionsLocation(ij,3) >  ij_Input) then
                        ij_Next     = ij
                        ij_Previous = ij - 1
                        exit
                    endif

                enddo

            endif i1

            if (ij_Next    == -99)  InsideSection = .true.
            if (ij_Previous   < 1)  InsideSection = .true.

    end subroutine LocateSection 



    logical function IsNearWaterOrPointWithNoValue(Grid_i, Grid_j)
        
        !Arguments-------------------------------------------------------------
        integer, intent(in)      :: Grid_i, Grid_j

        !Local-----------------------------------------------------------------
        integer                  :: i, j
        
        !Begin-----------------------------------------------------------------


        IsNearWaterOrPointWithNoValue = .false.

        do i = Grid_i - 1, Grid_i + 1
        do j = Grid_j - 1, Grid_j + 1

            if(i .ge. Me%ExtVar%WorkSize%ILB .and. i .le. Me%ExtVar%WorkSize%IUB .and. &
               j .ge. Me%ExtVar%WorkSize%JLB .and. j .le. Me%ExtVar%WorkSize%JUB) then

                if(Me%Depth(i,j) /= Me%LandPoint)then
                    IsNearWaterOrPointWithNoValue = .true.
                    return
                end if

            end if

        end do
        end do


    end function IsNearWaterOrPointWithNoValue
    
    !--------------------------------------------------------------------------
    
    subroutine CloseProject

        !Local-----------------------------------------------------------------
        real                                        :: ElapsedSeconds
        real                                        :: TotalCPUTime
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if(associated(Me%LandArea ))deallocate(Me%LandArea )
        if(associated(Me%XYZPoints))deallocate(Me%XYZPoints)
        if(associated(Me%XYZPointsOriginal))deallocate(Me%XYZPointsOriginal)
        
        if (Me%Overlapping%ON) call KillOverlapping

        call KillHorizontalGrid(HorizontalGridID= Me%ObjGrid, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'CloseProject - DigitalTerrainCreator - ERR10'

        deallocate(Me%Depth)

        call date_and_time(Values = Me%F95Time)
        call SetDate      (Me%FinalSystemTime,float(Me%F95Time(1)), float(Me%F95Time(2)),      &
                                              float(Me%F95Time(3)), float(Me%F95Time(5)),      &
                                              float(Me%F95Time(6)), float(Me%F95Time(7))+      &
                                              Me%F95Time(8)/1000.)
        call cpu_time(TotalCPUTime)
        ElapsedSeconds = Me%FinalSystemTime - Me%InitialSystemTime

        call ShutdownMohid ("Digital Terrain Creator", ElapsedSeconds, TotalCPUTime)


    end subroutine CloseProject

    !--------------------------------------------------------------------------
    
    subroutine KillOverlapping

        !Local-----------------------------------------------------------------
        integer                             :: ig, STAT_CALL

        !Begin-----------------------------------------------------------------

        do ig = 1, Me%Overlapping%Number

            if      (Me%Overlapping%DataInfo(ig)%InfoType == GridDataType) then

                call KillGridData(Me%Overlapping%DataInfo(ig)%ObjGridData, STAT = STAT_CALL) 
            
                if (STAT_CALL /= SUCCESS_) stop 'KillOverlapping - DigitalTerrainTool - ERR10'

            else if (Me%Overlapping%DataInfo(ig)%InfoType == CellsType   ) then

                deallocate(Me%Overlapping%DataInfo(ig)%Cells)
            
                if (STAT_CALL /= SUCCESS_) stop 'KillOverlapping - DigitalTerrainTool - ERR20'

            endif

        enddo        

        deallocate(Me%Overlapping%DataInfo)


    end subroutine KillOverlapping

    !--------------------------------------------------------------------------

end program DigitalTerrainCreator

