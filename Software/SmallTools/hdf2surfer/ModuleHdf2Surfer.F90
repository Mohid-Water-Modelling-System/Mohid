!------------------------------------------------------------------------------
!        Hidromod & IST , Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Formats conversion
! MODULE        : Hdf2Surfer
! URL           : http://www.mohid.com
! AFFILIATION   : Hidromod
! DATE          : June 2012
! REVISION      : Paulo Leitão
! DESCRIPTION   : Module to serve as Hdf2Surfer to create new modules
!
!------------------------------------------------------------------------------


Module ModuleHdf2Surfer

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModuleDrawing     
    use ModuleHDF5
    use ModuleHorizontalGrid
    use ModuleGridData
    use ModuleHorizontalMap 
    use ModuleGeometry     
    use ModuleMap 
    use ModuleStatistic
    use ModuleField4D
    use ModuleTimeSerie
    USE nrtype; USE nrutil; USE nr


    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructHdf2Surfer

    !Selector
    
    !Modifier
    public  :: ModifyHdf2Surfer

    !Destructor
    public  :: KillHdf2Surfer                                                     
    
    !Parameters----------------------------------------------------------------
    real,    parameter :: nodata   = 0.1701410E+39
    integer, parameter :: TSnumber = 10

    !Types---------------------------------------------------------------------
    private :: T_Unit
    type       T_Unit 
        integer                                     :: bat          = null_int
        integer                                     :: field        = null_int
        integer                                     :: fit          = null_int
        integer                                     :: list         = null_int
        integer                                     :: r            = null_int
        integer                                     :: rr           = null_int
        integer                                     :: bias         = null_int
        integer                                     :: rmse         = null_int
        integer                                     :: evol         = null_int
        integer                                     :: surf         = null_int
    end type   T_Unit

    private :: T_StatParam
    type       T_StatParam
        real                                        :: rcorr        = null_real
        real                                        :: rcorr_quad   = null_real
        real                                        :: bias         = null_real
        real                                        :: rmse         = null_real
        real                                        :: z_fisher     = null_real
        real                                        :: alfa         = null_real
        real                                        :: beta_1       = null_real
        real                                        :: Am           = null_real
        real                                        :: Bm           = null_real
    end type   T_StatParam
    
    private :: T_Stat
    type       T_Stat
        type (T_StatParam), dimension(:), pointer   :: Param                    => null()
        type (T_StatParam)                          :: GlobalParam
    end type   T_Stat
    
    
    private :: T_HDFSolution
    type       T_HDFSolution
        integer                                     :: InstanceID           = 0
        character(Len = StringLength)               :: Name                 = null_str
        character(Len = PathLength  )               :: FileIn               = null_str
        integer                                     :: MaskDim              = null_int
    end type  T_HDFSolution


    
    private :: T_Hdf2Surfer
    type       T_Hdf2Surfer
        integer                                     :: InstanceID           = 0
        type (T_Size2D)                             :: Size2D, WorkSize2D        
        type (T_Size3D)                             :: Size3D, WorkSize3D
        character(Len = StringLength)               :: Regiao               = null_str
        character(Len = StringLength)               :: Run                  = null_str
        character(Len = PathLength  )               :: BatimOut             = null_str
        logical                                     :: BatimOutON           = .false. 
        character(Len = PathLength  )               :: GeoOut               = null_str
        character(Len = PathLength  )               :: FileOutHDF           = null_str
        character(Len = PathLength  )               :: FileOutTS            = null_str        
        real                                        :: LatDefault           = null_real
        real                                        :: LongDefault          = null_real
        type (T_HDFSolution), dimension(:), pointer :: HDFSolution          => null()
        integer                                     :: HDFSolNumber         = null_int
        real,    dimension(:, :, :),  pointer       :: Depth                => null()
        real,    dimension(:, :   ),  pointer       :: Lat                  => null()
        real,    dimension(:, :   ),  pointer       :: Long                 => null()
        real,    dimension(:, :   ),  pointer       :: Bat                  => null()
        real,    dimension(:      ),  pointer       :: X                    => null()
        real,    dimension(:      ),  pointer       :: Y                    => null()
        real,    dimension(:      ),  pointer       :: Z                    => null()
        real,    dimension(:      ),  pointer       :: Xout                 => null()
        real,    dimension(:      ),  pointer       :: Yout                 => null()
        real,    dimension(:      ),  pointer       :: Zout                 => null()
        real,    dimension(:      ),  pointer       :: PropA                => null()
        real,    dimension(:      ),  pointer       :: PropB                => null()
        real,    dimension(:      ),  pointer       :: PropAout             => null()
        real,    dimension(:      ),  pointer       :: PropBout             => null()
        logical, dimension(:      ),  pointer       :: NodataA              => null()
        logical, dimension(:      ),  pointer       :: NodataB              => null()
        integer, dimension(:, :, :),  pointer       :: WaterPoints3D        => null()
        integer, dimension(:, :   ),  pointer       :: WaterPoints2D        => null()
        type (T_Stat)                               :: Statistic            
        integer                                     :: NPoints              = null_int
        integer                                     :: NPointsout           = null_int
        type (T_Time)                               :: BeginTime
        type (T_Time)                               :: EndTime
        type (T_Time)                               :: InitialSystemTime
        type (T_Time)                               :: CurrentTime
        real                                        :: DT                   = null_real

        logical                                     :: BeginTimeIn          = .false. 
        logical                                     :: EndTimeIn            = .false. 

        
        integer                                     :: Nout                 = null_int
        logical                                     :: HDFOut               = .false.
        logical                                     :: SurferOut            = .false.
        logical                                     :: Out3D                = .false. 
        character(Len = StringLength)               :: ScalarProperty       = null_str
        character(Len = StringLength)               :: ScalarPropUnits      = null_str
        character(Len = StringLength)               :: ScalarPropOut        = null_str
        integer                                     :: ScalarPropertyID     = null_int
        real                                        :: ValueMin             = null_real
        real                                        :: ValueMax             = null_real
        type (T_Unit)                               :: Unit
        logical                                     :: Extrapolate          = .false. 
        logical                                     :: AngleProperty        = .false. 
        integer                                     :: ObjEnterData         = 0
        integer                                     :: ObjTime              = 0
        integer                                     :: ObjBathymetry        = 0  
        integer                                     :: ObjHorizontalMap     = 0
        integer                                     :: ObjHorizontalGrid    = 0
        integer                                     :: ObjGeometry          = 0
        integer                                     :: ObjMap               = 0
        integer                                     :: ObjHDF5Out           = 0
        integer                                     :: ObjTimeSerieOut      = 0
        
    end type  T_Hdf2Surfer

    !Global Module Variables
    type (T_Hdf2Surfer), pointer                         :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructHdf2Surfer(ObjHdf2SurferID, InitialSystemTime, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjHdf2SurferID 
        type (T_Time)                                   :: InitialSystemTime
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, STAT_CALL

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_
        
        allocate(Me)

        Me%InstanceID = 1
        
        Me%InitialSystemTime = InitialSystemTime

        call ConstructEnterData(Me%ObjEnterData, "Hdf2Surfer.dat", STAT = STAT_CALL) 
        
        if (STAT_CALL /= SUCCESS_) then
            stop "ConstructHdf2Surfer - ModuleHdf2Surfer - ERR10"
        endif
         
        call ConstrucTimeInfo
        call ConstructOptions
            
        !read input options
        call ConstructInPut

        !read output options
        call ConstructOutPut   
       
       
        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
        
        if (STAT_CALL /= SUCCESS_) then
            stop "ConstructHdf2Surfer - ModuleHdf2Surfer - ERR20"
        endif
        
        call ConstructTS
        
        call Construct1DAnalysis
        

        !Returns ID
        ObjHdf2SurferID          = Me%InstanceID

        STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructHdf2Surfer

    !--------------------------------------------------------------------------

    subroutine ConstructOptions
        !Arguments---------------------------------------------------------------
        !Local-------------------------------------------------------------------
        real                                            :: aux, DT
        integer                                         :: STAT_CALL, flag, FirstLine, LastLine
        

        !Begin-------------------------------------------------------------------
        
        
        
        call GetData(   Me%BatimOut,                                                    &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'BATIM_OUT',                                     &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOptions - ModuleHdf2Surfer - ERR20"
        endif

        if (flag == 0) then            
            Me%BatimOutON = .false.        
        else
            Me%BatimOutON = .true.
        endif

        call GetData(   Me%GeoOut,                                                      &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'GEO_OUT',                                       &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOptions - ModuleHdf2Surfer - ERR40"
        endif

        if (flag == 0) then            
            Me%Out3D = .false.
        else
            Me%Out3D = .true.
        endif

        
        call GetData(   Me%HDFOut,                                                      &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'HDF_OUT_ON',                                    &
                        default      = .false.,                                         &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOptions - ModuleHdf2Surfer - ERR110"
        endif

        call GetData(   Me%SurferOut,                                                   &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'SURFER_OUT',                                    &
                        default      = .false.,                                         &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOptions - ModuleHdf2Surfer - ERR120"
        endif

        if (.not. Me%SurferOut .and. .not. Me%HDFOut) then
            write(*,*) 'You need to choose at least one output option Surfer or HDF or both'
            stop "ConstructOptions - ModuleHdf2Surfer - ERR130"
        endif


        call GetData(   Me%ScalarProperty,                                              &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'SCALAR_PROPERTY',                               &
                        STAT         = STAT_CALL)
        

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOptions - ModuleHdf2Surfer - ERR140"
        endif

        if (flag == 0) then            
            stop "ConstructOptions - ModuleHdf2Surfer - ERR150"
        endif  
        
        Me%ScalarPropertyID = GetPropertyIDNumber (Me%ScalarProperty)


        call GetData(   Me%AngleProperty,                                               &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'ANGLE_PROPERTY',                                &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOptions - ModuleHdf2Surfer - ERR160"
        endif

        if (flag == 0) then            
            Me%AngleProperty = Check_Angle_Property(Property = Me%ScalarPropertyID)
        endif  

    end subroutine ConstructOptions
    
     
    !--------------------------------------------------------------------------

    subroutine ConstructOutPut

        !Arguments---------------------------------------------------------------
        !Local-------------------------------------------------------------------
        real,   dimension(:,:), pointer                 :: SurfaceElevation
        real                                            :: aux, DT
        integer                                         :: STAT_CALL, flag, FirstLine, LastLine
        

        !Begin-------------------------------------------------------------------
        
        if (Me%BatimOutON) then

            call ConstructHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,       &
                                         DataFile         = Me%BatimOut,                &
                                         STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR10'
            endif

            !Horizontal Grid Data - Water Column (Bathymetry)
            call ConstructGridData      (GridDataID       = Me%ObjBathymetry,           &
                                         HorizontalGridID = Me%ObjHorizontalGrid,       &
                                         TimeID           = Me%ObjTime,                 &
                                         FileName         = Me%BatimOut,                &
                                         STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR20'
            endif
            
        else
        
            !it = 1 - model solution                
            call GetField4DGridID (Field4DID = Me%HDFSolution(1)%InstanceID,            &
                                   GridID    = Me%ObjHorizontalGrid,                    &
                                   STAT      = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR30'
            endif

            !it = 1 - model solution                
            call GetField4DBathymID (Field4DID  = Me%HDFSolution(1)%InstanceID,         &
                                     BathymID   = Me%ObjBathymetry,                     &
                                     STAT       = STAT_CALL)
                                     
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR40'
            endif
            
                
        endif
                    
        call GetLatitudeLongitude(HorizontalGridID   = Me%ObjHorizontalGrid,            &
                                  Latitude           = Me%LatDefault,                   &
                                  Longitude          = Me%LongDefault,                  &
                                  STAT               = STAT_CALL)        

        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructOutPut - ModuleHdf2Surfer - ERR50'
        endif

        
        call GetHorizontalGridSize(Me%ObjHorizontalGrid, Size = Me%Size2D,              &
                                   WorkSize = Me%WorkSize2D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructOutPut - ModuleHdf2Surfer - ERR60'
        endif


        !Horizontal Map
        call ConstructHorizontalMap (HorizontalMapID  = Me%ObjHorizontalMap,            &
                                     GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     ActualTime       = Me%BeginTime,                   &
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructOutPut - ModuleHdf2Surfer - ERR70'
        endif

        if (Me%Out3D) then

            !Geometry - Water Column
            call ConstructGeometry      (GeometryID       = Me%ObjGeometry,             &
                                         GridDataID       = Me%ObjBathymetry,           &
                                         HorizontalGridID = Me%ObjHorizontalGrid,       &
                                         HorizontalMapID  = Me%ObjHorizontalMap,        &
                                         ActualTime       = Me%BeginTime,               &
                                         NewDomain        = Me%GeoOut,                  &
                                         STAT             = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR80'
            endif


            call ConstructMap           (Map_ID           = Me%ObjMap,                  &
                                         GeometryID       = Me%ObjGeometry,             &
                                         HorizontalMapID  = Me%ObjHorizontalMap,        &
                                         TimeID           = Me%ObjTime,                 &
                                         GridDataID       = Me%ObjBathymetry,           &
                                         HorizontalGridID = Me%ObjHorizontalGrid,       &
                                         STAT             = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR90'
            endif
            
            call GetWaterPoints3D(Me%ObjMap, Me%WaterPoints3D, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)  then
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR100'
            endif
            
            allocate(SurfaceElevation(Me%WorkSize2D%ILB:Me%WorkSize2D%IUB,Me%WorkSize2D%JLB:Me%WorkSize2D%JUB))
            
            SurfaceElevation(:,:) = 0.
            
            call ComputeInitialGeometry(GeometryID       = Me%ObjGeometry,              &
                                        SurfaceElevation = SurfaceElevation,            &
                                        WaterPoints3D    = Me%WaterPoints3D,            &
                                        STAT             = STAT_CALL)
                                                 
            if (STAT_CALL /= SUCCESS_)  then
                stop "Construct1DAnalysis - ModuleHdf2Surfer - ERR110"
            endif
            
            deallocate(SurfaceElevation)

            call UnGetMap(Me%ObjMap, Me%WaterPoints3D, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)  then
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR120'
            endif

            call GetGeometrySize(Me%ObjGeometry, Size = Me%Size3D,                      &
                                 WorkSize = Me%WorkSize3D,  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR130'
            endif

        
        endif
        
        allocate(Me%Statistic%Param(Me%Nout))
        
        if (Me%SurferOut) then
        
            call GetData(   Me%Regiao,                                                      &
                            Me%ObjEnterData, flag,                                          &
                            SearchType   = FromFile,                                        &
                            ClientModule = 'ModuleHdf2Surfer',                              &
                            keyword      = 'REGIAO',                                        &
                            STAT         = STAT_CALL)
            

            if (STAT_CALL /= SUCCESS_)  then
                stop "ConstructOutPut - ModuleHdf2Surfer - ERR140"
            endif

            if (flag == 0)  then
                stop "ConstructOutPut - ModuleHdf2Surfer - ERR150"
            endif        

            call GetData(   Me%Run,                                                         &
                            Me%ObjEnterData, flag,                                          &
                            SearchType   = FromFile,                                        &
                            ClientModule = 'ModuleHdf2Surfer',                              &
                            keyword      = 'RUN',                                           &
                            STAT         = STAT_CALL)
            

            if (STAT_CALL /= SUCCESS_)  then
                stop "ConstructOutPut - ModuleHdf2Surfer - ERR160"
            endif

            if (flag == 0)  then
                stop "ConstructOutPut - ModuleHdf2Surfer - ERR170"
            endif        
        
        endif


        call GetData(   Me%ScalarPropOut,                                               &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'SCALAR_PROP_OUT',                               &
                        STAT         = STAT_CALL)
        

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR180"
        endif

        if (flag == 0) then            
            Me%ScalarPropOut = Me%ScalarProperty
        endif  

        call GetData(   Me%ScalarPropUnits,                                             &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'SCALAR_PROP_UNITS',                             &
                        STAT         = STAT_CALL)
        

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR190"
        endif

        if (flag == 0) then            
            Me%ScalarPropUnits = ' '
        endif  

        call GetData(   Me%ValueMin,                                                    &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'VALUE_MIN',                                     &
                        default      = FillValueReal,                                   &
                        STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR200"
        endif

        call GetData(   Me%ValueMax,                                                    &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'VALUE_MAX',                                     &
                        default      = -FillValueReal,                                  &
                        STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR210"
        endif

        call GetData(   Me%Extrapolate,                                                 &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'EXTRAPOLATE',                                   &
                        default      = .false.,                                         &
                        STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR220"
        endif
        
        if (Me%HDFOut) then
        
            call GetData(   Me%FileOutHDF,                                              &
                            Me%ObjEnterData, flag,                                      &
                            SearchType   = FromFile,                                    &
                            ClientModule = 'ModuleHdf2Surfer',                          &
                            keyword      = 'HDF_OUT',                                   &
                            STAT         = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)  then
                stop "ConstructOutPut - ModuleHdf2Surfer - ERR230"
            endif

            if (flag == 0) then            
                stop "ConstructOutPut - ModuleHdf2Surfer - ERR240"
            endif          

            call GetData(   Me%FileOutTS,                                               &
                            Me%ObjEnterData, flag,                                      &
                            SearchType   = FromFile,                                    &
                            ClientModule = 'ModuleHdf2Surfer',                          &
                            keyword      = 'TS_OUT',                                    &
                            STAT         = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)  then
                stop "ConstructOutPut - ModuleHdf2Surfer - ERR250"
            endif

            if (flag == 0) then            
                stop "ConstructOutPut - ModuleHdf2Surfer - ERR260"
            endif          

        endif
        


        !----------------------------------------------------------------------

    end subroutine ConstructOutPut
 
    !--------------------------------------------------------------------------

    subroutine ConstructTS  
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        integer                                             :: STAT_CALL

        !----------------------------------------------------------------------   
        
        allocate(PropertyList(1:TSnumber))
        
        PropertyList( 1) = "R"     
        PropertyList( 2) = "R*R"
        PropertyList( 3) = "Bias"
        PropertyList( 4) = "RMSE" 
        PropertyList( 5) = "z_fisher"  
        PropertyList( 6) = "a"
        PropertyList( 7) = "b"
        PropertyList( 8) = "Am"     
        PropertyList( 9) = "Bm"
        PropertyList(10) = "n.total"

        if (Me%Out3D) then

            call StartTimeSerie(TimeSerieID         = Me%ObjTimeSerieOut,               &
                                ObjTime             = Me%ObjTime,                       &
                                TimeSerieDataFile   = "Hdf2Surfer.dat",                 &
                                PropertyList        = PropertyList,                     &
                                Extension           = "",                               &
                                WaterPoints3D       = Me%WaterPoints3D,                 &
                                ResultFileName      = trim(Me%FileOutTS),               &                                
                                STAT                = STAT_CALL)
            if (STAT_CALL /= 0) stop 'ConstructTS - ModuleHdf2Surfer - ERR10'

        else

            call StartTimeSerie(TimeSerieID         = Me%ObjTimeSerieOut,               &
                                ObjTime             = Me%ObjTime,                       &
                                TimeSerieDataFile   = "Hdf2Surfer.dat",                 &
                                PropertyList        = PropertyList,                     &
                                Extension           = "",                               &
                                WaterPoints2D       = Me%WaterPoints2D,                 &
                                ResultFileName      = trim(Me%FileOutTS),               &                                  
                                STAT                = STAT_CALL)
            if (STAT_CALL /= 0) stop 'ConstructTS - ModuleHdf2Surfer - ERR20'        
        
        endif
        
        deallocate(PropertyList)

        !----------------------------------------------------------------------

    end subroutine ConstructTS
    
    !--------------------------------------------------------------------------    
    
   !----------------------------------------------------------------------------


    type(T_Time) function HDF5TimeInstant(Instant, HDF5ID_)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        integer, optional                       :: HDF5ID_
        

        !Local-----------------------------------------------------------------
        real,    dimension(:), pointer          :: TimeVector
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        allocate(TimeVector(6))
        
        call HDF5SetLimits  (HDF5ID_, 1, 6, STAT = STAT_CALL)        

        call HDF5ReadWindow (HDF5ID         = HDF5ID_,                                  &
                             GroupName      = "/Time",                                  &
                             Name           = "Time",                                   &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleField4D - ERR01'

        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))

                                     
        deallocate(TimeVector)

    end function HDF5TimeInstant

    
    !--------------------------------------------------------------------------    
    
    subroutine ConstrucTimeInfo    
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(Len=PathLength  )                 :: FileName
        character(Len=StringLength)                 :: solution_begin, solution_end
        real                                        :: aux
        integer                                     :: flag, ObjHDF5In, iend, STAT_CALL
        integer                                     :: ClientNumber, HDF5_READ     
        logical                                     :: PropertyFound 

        !----------------------------------------------------------------------     
    

        call GetData(   Me%BeginTime,                                                   &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'START',                                         &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstrucTimeInfo - ModuleHdf2Surfer - ERR10"
        endif

        if (flag == 0) then            
            Me%BeginTimeIn = .false.
        else
            Me%BeginTimeIn = .true.                    
        endif        

        call GetData(   Me%EndTime,                                                     &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'END',                                           &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstrucTimeInfo - ModuleHdf2Surfer - ERR20"
        endif

        if (flag == 0) then            
            Me%EndTimeIn = .false.
        else
            Me%EndTimeIn = .true.
        endif 
        
        
        
        if (.not. Me%BeginTimeIn .or. .not. Me%EndTimeIn) then
            
            solution_begin = "<beginobservations>"
            solution_end   = "<endobservations>"
        
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                  &
                                        solution_begin, solution_end,                   &
                                        PropertyFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop "ConstrucTimeInfo - ModuleHdf2Surfer - ERR30"
            
            if (.not. PropertyFound) then

                write(*,*) 'Start Block'
                write(*,*) trim(solution_begin)
                write(*,*) trim(solution_end  )
                write(*,*) 'End Block'
                stop "ConstructInPut - ModuleHdf2Surfer - ERR20" 
            
            endif
            
            call GetData(   Filename,                                                   &
                            Me%ObjEnterData, flag,                                      &
                            SearchType   = FromBlock,                                   &
                            ClientModule = 'ModuleHdf2Surfer',                          &
                            keyword      = 'HDF_IN',                                    &
                            STAT         = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)  then
                stop "ConstrucTimeInfo - ModuleHdf2Surfer - ERR40"
            endif
            

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop "ConstructInPut - ModuleHdf2Surfer - ERR30" 

            !Prepares file for a new block search throughout the entire file
            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop "ConstructInPut - ModuleHdf2Surfer - ERR40" 
                                
            

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
            
            ObjHDF5In = 0

            !Opens HDF5 File
            call ConstructHDF5      (ObjHDF5In, Filename, HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop "ConstrucTimeInfo - ModuleHdf2Surfer - ERR50"
            
           
            if (.not. Me%BeginTimeIn) then
                Me%BeginTime = HDF5TimeInstant(1, ObjHDF5In)
            endif
                    
            if (.not. Me%EndTimeIn) then
                call GetHDF5GroupNumberOfItems(ObjHDF5In, "/Time", iend, STAT = STAT_CALL)        
                Me%EndTime = HDF5TimeInstant(iend, ObjHDF5In)
            endif
            
            !Kill HDF5 File
            call KillHDF5      (ObjHDF5In, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop "ConstrucTimeInfo - ModuleHdf2Surfer - ERR60"
 
        
        endif
        

        call GetData(   Me%DT,                                                          &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'DT',                                            &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstrucTimeInfo - ModuleHdf2Surfer - ERR80"
        endif

        if (flag == 0) then            
            stop "ConstrucTimeInfo - ModuleHdf2Surfer - ERR90"
        endif     
        
        aux  = (Me%EndTime - Me%BeginTime)/Me%DT
        
        if (aux /= real(int(aux))) then
            write(*,*) "DT must be a submultiple of the analysed period"
            stop "ConstrucTimeInfo - ModuleHdf2Surfer - ERR100"
        endif
           
        Me%Nout = int(aux) + 1        


        call StartComputeTime(TimeID            = Me%ObjTime,                           &
                              InitialSystemTime = Me%InitialSystemTime,                 &
                              BeginTime         = Me%BeginTime,                         &
                              EndTime           = Me%EndTime,                           &
                              DT                = Me%DT,                                &
                              VariableDT        = .false.,                              &
                              STAT              = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstrucTimeInfo - ModuleHdf2Surfer - ERR110'
        endif
        
            
        
    end subroutine ConstrucTimeInfo                     

    !--------------------------------------------------------------------------        

    subroutine ConstructInPut

        !Arguments---------------------------------------------------------------
        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL, flag, it, ClientNumber
        character(StringLength)                         :: solution_begin, solution_end
        logical                                         :: PropertyFound
        type (T_Time)                                   :: StartTime, EndTime
        real                                            :: aux

        !Begin-------------------------------------------------------------------
        
        !Number of input solutions = 2
        Me%HDFSolNumber = 2
        
        allocate(Me%HDFSolution(1:Me%HDFSolNumber))
        
        
        do it=1,Me%HDFSolNumber
        
            if     (it == 1) then
                solution_begin = "<beginmodel>"
                solution_end   = "<endmodel>"
            elseif (it == 2) then
                solution_begin = "<beginobservations>"
                solution_end   = "<endobservations>"
            endif
        
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                  &
                                        solution_begin, solution_end,                   &
                                        PropertyFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop "ConstructInPut - ModuleHdf2Surfer - ERR10"
            
i1:         if (PropertyFound) then
                call ConstructOneSolution(it)
            else
                write(*,*) 'Start Block'
                write(*,*) trim(solution_begin)
                write(*,*) trim(solution_end  )
                write(*,*) 'End Block'
                stop "ConstructInPut - ModuleHdf2Surfer - ERR20" 
            endif i1
            
        enddo

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop "ConstructInPut - ModuleHdf2Surfer - ERR30" 

        !Prepares file for a new block search throughout the entire file
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "ConstructInPut - ModuleHdf2Surfer - ERR40" 
                
        
    end subroutine ConstructInPut
    
    !--------------------------------------------------------------------------    
                 
    subroutine ConstructOneSolution(it)     

        !Arguments---------------------------------------------------------------
        integer                                         :: it
        
        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL, flag

        !Begin-------------------------------------------------------------------    
        
        call GetData(   Me%HDFSolution(it)%FileIn,                                       &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromBlock,                                       &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'HDF_IN',                                        &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOneSolution - ModuleHdf2Surfer - ERR10"
        endif

        if (flag == 0) then            
            stop "ConstructOneSolution - ModuleHdf2Surfer - ERR20"
        endif  
        
       
        call GetData(   Me%HDFSolution(it)%MaskDim,                                     &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromBlock,                                       &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'MASK_DIM',                                      &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOneSolution - ModuleHdf2Surfer - ERR50"
        endif

        if (flag == 0) then            
            stop "ConstructOneSolution - ModuleHdf2Surfer - ERR60"
        endif

        call ConstructField4D(Field4DID     = Me%HDFSolution(it)%InstanceID,            &
                              EnterDataID   = Me%ObjEnterData,                          &
                              ExtractType   = FromBlock,                                &
                              FileName      = Me%HDFSolution(it)%FileIn,                &
                              TimeID        = Me%ObjTime,                               &   
                              MaskDim       = Me%HDFSolution(it)%MaskDim,               &
                              LatReference  = Me%LatDefault,                            &
                              LonReference  = Me%LongDefault,                           & 
                              Extrapolate   = Me%Extrapolate,                           &    
                              STAT          = STAT_CALL)
                              
        if (STAT_CALL /= SUCCESS_) then
            stop "ConstructOneSolution - ModuleHdf2Surfer - ERR70"
        endif            
        
            
    end subroutine ConstructOneSolution

    !----------------------------------------------------------------------------
    
    subroutine Construct1DAnalysis    

        !Arguments---------------------------------------------------------------
        integer                                         :: i, j, k, p
        
        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Begin-------------------------------------------------------------------    
    
        if (Me%Out3D) then
            call GetWaterPoints3D(Me%ObjMap, Me%WaterPoints3D, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)  then
                stop "Construct1DAnalysis - ModuleHdf2Surfer - ERR10"
            endif
                
            Me%NPoints = sum(Me%WaterPoints3D)
            
        else
            call GetWaterPoints2D(Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
            
            if (STAT_CALL /= SUCCESS_)  then
                stop "Construct1DAnalysis - ModuleHdf2Surfer - ERR20"
            endif
            
            Me%NPoints = sum(Me%WaterPoints2D)
        endif
        
        allocate(Me%X(1:Me%NPoints), Me%Y(1:Me%NPoints), Me%Z(1:Me%NPoints))
        
        allocate(Me%PropA(1:Me%NPoints), Me%PropB(1:Me%NPoints))
        
        Me%PropA(:) = FillValueReal
        Me%PropB(:) = FillValueReal

        allocate(Me%NodataA(1:Me%NPoints), Me%NodataB(1:Me%NPoints))        
        
        call GetGridLatitudeLongitude(HorizontalGridID  = Me%ObjHorizontalGrid,         &
                                      GridLatitude      = Me%Lat,                       &
                                      GridLongitude     = Me%Long,                      &        
                                      STAT              = STAT_CALL)
            
        if (STAT_CALL /= SUCCESS_)  then
            stop "Construct1DAnalysis - ModuleHdf2Surfer - ERR30"
        endif
        
        if (Me%Out3D) then
        
            call GetGeometryDistances(GeometryID    = Me%ObjGeometry,               &
                                      ZCellCenter   = Me%Depth,                     &
                                      STAT          = STAT_CALL)
            
            if (STAT_CALL /= SUCCESS_)  then
                stop "Construct1DAnalysis - ModuleHdf2Surfer - ERR40"
            endif
        
            p = 1
        
            do k=Me%WorkSize3D%KLB,Me%WorkSize3D%KUB
            do j=Me%WorkSize3D%JLB,Me%WorkSize3D%JUB
            do i=Me%WorkSize3D%ILB,Me%WorkSize3D%IUB                        
                
                if (Me%WaterPoints3D(i, j, k) == WaterPoint)  then
                
                    Me%X(p) =   Me%Long (i, j)
                    Me%Y(p) =   Me%Lat  (i, j)                    
                    Me%Z(p) = - Me%Depth(i, j, k)

                    p = p + 1
                endif
                
            enddo
            enddo
            enddo
        
        else
        
            p = 1
        
            do j=Me%WorkSize2D%JLB,Me%WorkSize2D%JUB
            do i=Me%WorkSize2D%ILB,Me%WorkSize2D%IUB                        
                
                if (Me%WaterPoints2D(i, j) == WaterPoint)  then
                
                    Me%X(p) =   Me%Long (i, j)
                    Me%Y(p) =   Me%Lat  (i, j)   
                    Me%Z(p) =   0.                 

                    p = p + 1
                endif
                
            enddo
            enddo
        
        endif
        
        !initialization
        Me%Statistic%GlobalParam%rcorr         = 0.
        Me%Statistic%GlobalParam%rcorr_quad    = 0.
        Me%Statistic%GlobalParam%bias          = 0.
        Me%Statistic%GlobalParam%rmse          = 0.
        Me%Statistic%GlobalParam%z_fisher      = 0.
        Me%Statistic%GlobalParam%Am            = 0.
        Me%Statistic%GlobalParam%Bm            = 0.
        Me%Statistic%GlobalParam%alfa          = 0.
        Me%Statistic%GlobalParam%BETA_1        = 0.        
        
    
    end  subroutine Construct1DAnalysis   
    
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyHdf2Surfer(ObjHdf2SurferID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjHdf2SurferID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, STAT_CALL, it, io, i

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_
        

        Me%CurrentTime = Me%BeginTime
        
        if (Me%SurferOut) then
            call WriteBathymetry
            call WriteGlobalIntial
        endif
        
        if (Me%HDFOut   ) then        
            call WriteGlobalIntialHDF
        endif            
        
        it = 1
        
        do while (Me%CurrentTime <=Me%EndTime)

            Me%NoDataA(:) = .true.
            Me%NoDataB(:) = .true.
            
    
            call ModifyField4DXYZ(Field4DID             = Me%HDFSolution(1)%InstanceID, &
                                  PropertyIDNumber      = Me%ScalarPropertyID,          &
                                  CurrentTime           = Me%CurrentTime,               &
                                  X                     = Me%X,                         &
                                  Y                     = Me%Y,                         &
                                  Z                     = Me%Z,                         &
                                  Field                 = Me%PropA,                     &
                                  NoData                = Me%NoDataA,                   &
                                  STAT                  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  then
                stop "ModifyHdf2Surfer - ModuleHdf2Surfer - ERR10"
            endif
                                  

            call ModifyField4DXYZ(Field4DID             = Me%HDFSolution(2)%InstanceID, &
                                  PropertyIDNumber      = Me%ScalarPropertyID,          &
                                  CurrentTime           = Me%CurrentTime,               &
                                  X                     = Me%X,                         &
                                  Y                     = Me%Y,                         &
                                  Z                     = Me%Z,                         &
                                  Field                 = Me%PropB,                     &
                                  NoData                = Me%NoDataB,                   &
                                  STAT                  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  then
                stop "ModifyHdf2Surfer - ModuleHdf2Surfer - ERR20"
            endif

                                  
            io = 0                                  
            do i=1, Me%NPoints
                if (.not.Me%NoDataA(i))then
                    if (Me%PropA(i) < Me%ValueMin .or. Me%PropA(i) > Me%ValueMax) then
                        Me%NoDataA(i) = .true.
                    endif
                endif
                if (.not.Me%NoDataB(i))then
                    if (Me%PropB(i) < Me%ValueMin .or. Me%PropB(i) > Me%ValueMax) then
                        Me%NoDataB(i) = .true.
                    endif
                endif
                    
                if (Me%NoDataA(i) .or. Me%NoDataB(i)) then
                    cycle
                else
                    io = io + 1
                endif
            enddo 
            Me%NPointsOut = io
            
            
            allocate(Me%Xout    (1:Me%NPointsOut))
            allocate(Me%Yout    (1:Me%NPointsOut))
            allocate(Me%Zout    (1:Me%NPointsOut))
            allocate(Me%PropAout(1:Me%NPointsOut))
            allocate(Me%PropBout(1:Me%NPointsOut))
            
            io = 0
            
            do i=1, Me%NPoints
                if (Me%NoDataA(i) .or. Me%NoDataB(i)) then
                    cycle
                else
                    io = io + 1
                    Me%Xout    (io) = Me%X    (i)
                    Me%Yout    (io) = Me%Y    (i)
                    Me%Zout    (io) = Me%Z    (i)
                    Me%PropAout(io) = Me%PropA(i)
                    Me%PropBout(io) = Me%PropB(i)
                endif                    
            enddo             

            call ModifyHdfInstant(it)
            
            Me%CurrentTime =  Me%CurrentTime + Me%DT
            
            call ActualizeCurrentTime(TimeID    = Me%ObjTime,                           &
                                      DT_Global = Me%DT,                                &
                                      Current   = Me%CurrentTime,                       &         
                                      STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  then
                stop "ModifyHdf2Surfer - ModuleHdf2Surfer - ERR30"
            endif
            
            it = it + 1

            deallocate(Me%Xout    )
            deallocate(Me%Yout    )
            deallocate(Me%Zout    )
            deallocate(Me%PropAout)
            deallocate(Me%PropBout)            
        enddo
        
        if (Me%SurferOut) then
            call WriteGlobalEnd
        endif
        
        if (Me%HDFOut   ) then

            call KillHDF5 (Me%ObjHDF5Out, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop "ModifyHdf2Surfer - ModuleHdf2Surfer - ERR40"
            endif        
            call KillTimeSerie(Me%ObjTimeSerieOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop "ModifyHdf2Surfer - ModuleHdf2Surfer - ERR50"
            endif                       
            
            call HdfStatisticAnalysis
        endif
        
        

        STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyHdf2Surfer
    
    !---------------------------------------------------------------------------
    
    
    subroutine WriteGlobalIntialHDF    
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,   dimension(:,:), pointer                     :: Bathymetry
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        integer                                             :: STAT_CALL
        integer                                             :: WorkILB, WorkIUB
        integer                                             :: WorkJLB, WorkJUB
        integer                                             :: WorkKLB, WorkKUB
        integer                                             :: HDF5_CREATE

        !----------------------------------------------------------------------   
        if (Me%Out3D) then

            WorkJLB = Me%WorkSize3D%JLB
            WorkJUB = Me%WorkSize3D%JUB
            WorkILB = Me%WorkSize3D%ILB
            WorkIUB = Me%WorkSize3D%IUB           
            WorkKLB = Me%WorkSize3D%KLB
            WorkKUB = Me%WorkSize3D%KUB           
        
        else
        
            WorkJLB = Me%WorkSize2D%JLB
            WorkJUB = Me%WorkSize2D%JUB
            WorkILB = Me%WorkSize2D%ILB
            WorkIUB = Me%WorkSize2D%IUB           
    
        endif

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        Me%ObjHDF5Out = 0

        !Opens HDF5 File
        call ConstructHDF5      (Me%ObjHDF5Out, trim(Me%FileOutHDF),                    &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGlobalIntialHDF - ModuleHdf2Surfer - ERR10'
        
            
        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5Out,                   &
                                 WorkSize = Me%WorkSize2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGlobalIntialHDF - ModuleHdf2Surfer - ERR20'

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5Out, WorkILB, WorkIUB, WorkJLB,                 &
                              WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGlobalIntialHDF - ModuleHdf2Surfer - ERR30'
        
        !Horizontal Grid Data - Water Column (Bathymetry)
        call GetGridData      (GridDataID       = Me%ObjBathymetry,                     &
                               GridData2D       = Bathymetry,                           &
                               STAT             = STAT_CALL)
                               
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalIntialHDF - ModuleHdf2Surfer - ERR40'
        endif            


        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5Out, "/Grid", "Bathymetry", "m",                &
                              Array2D = Bathymetry,                                     &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGlobalIntialHDF - ModuleHdf2Surfer - ERR50'
        
        call UnGetGridData   (GridDataID       = Me%ObjBathymetry,                      &
                              Array            = Bathymetry,                            &
                              STAT             = STAT_CALL)
                               
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalIntialHDF - ModuleHdf2Surfer - ERR60'
        endif            
        
        
        if (Me%Out3D) then
        
            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5Out, WorkILB, WorkIUB, WorkJLB,             &
                                  WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteGlobalIntialHDF - ModuleHdf2Surfer - ERR70'


            call HDF5WriteData   (Me%ObjHDF5Out, "/Grid", "WaterPoints3D", "-",         &
                                  Array3D = Me%WaterPoints3D,                           &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteGlobalIntialHDF - ModuleHdf2Surfer - ERR80'
            
        else


            call HDF5WriteData   (Me%ObjHDF5Out, "/Grid", "WaterPoints2D", "-",         &
                                  Array2D = Me%WaterPoints2D,                           &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteGlobalIntialHDF - ModuleHdf2Surfer - ERR90'        
        
        endif            

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5Out, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGlobalIntialHDF - ModuleHdf2Surfer - ERR100'
        
        !----------------------------------------------------------------------

    end subroutine WriteGlobalIntialHDF
    
  
    !--------------------------------------------------------------------------
    
    subroutine WriteGlobalIntial
    
        !Arguments-------------------------------------------------------------    

        !Local-----------------------------------------------------------------
        character(Len=PathLength)                   :: FileName
        integer                                     :: STAT_CALL


        call UnitsManager(Me%Unit%List, OPEN_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalIntial - ModuleHdf2Surfer - ERR10'
        endif    	    

        !FileName = trim(Me%ScalarProperty)//"_"//trim(DateName(Me%CurrentTime))//".dat"
        FileName = "SST"//"_"//trim(DateName(Me%CurrentTime))//".dat"
        
    
		!open(Me%Unit%List, file='lista_'//trim(Me%ScalarProperty)//trim(adjustl(Me%Run))//'.txt', status='unknown')
		open(Me%Unit%List, file='lista_'//"SST"//trim(adjustl(Me%Run))//'.txt', status='unknown')
		write(Me%Unit%List,*) Me%Nout


        call UnitsManager(Me%Unit%Evol, OPEN_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalIntial - ModuleHdf2Surfer - ERR20'
        endif               

		open (Me%Unit%Evol, file='evolucao_semana'//trim(adjustl(Me%Run))//'_'//trim(Me%Regiao)//'.dat', status='unknown')
		
		write(Me%Unit%Evol,'(A)') 'dia      Time     R     R*R   Bias   RMSE  z_fisher  a       b      Am     Bm    nx  ny  n.total'

        call UnitsManager(Me%Unit%r, OPEN_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalIntial - ModuleHdf2Surfer - ERR30'
        endif               

		open(Me%Unit%r,    file='evolucao_R'   //trim(adjustl(Me%Run))//'.bln', status='unknown')

        call UnitsManager(Me%Unit%rr, OPEN_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalIntial - ModuleHdf2Surfer - ERR40'
        endif               
		
		open(Me%Unit%rr,   file='evolucao_RR'  //trim(adjustl(Me%Run))//'.bln', status='unknown')

        call UnitsManager(Me%Unit%bias, OPEN_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalIntial - ModuleHdf2Surfer - ERR50'
        endif               
		
		open(Me%Unit%bias, file='evolucao_bias'//trim(adjustl(Me%Run))//'.bln', status='unknown')
		
        call UnitsManager(Me%Unit%rmse, OPEN_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalIntial - ModuleHdf2Surfer - ERR60'
        endif               
		
		open(Me%Unit%rmse, file='evolucao_rmse'//trim(adjustl(Me%Run))//'.bln', status='unknown')

		write(Me%Unit%r,   *) Me%Nout, 1
		write(Me%Unit%rr,  *) Me%Nout, 1
		write(Me%Unit%bias,*) Me%Nout, 1
		write(Me%Unit%rmse,*) Me%Nout, 1
    
    end subroutine WriteGlobalIntial
    
    !---------------------------------------------------------------------------    
    
    
    subroutine WriteGlobalEnd
    
        !Arguments-------------------------------------------------------------    

        !Local-----------------------------------------------------------------
        character(Len=StringLength)                 :: produto_name,  unidades_produto_name, &
                                                       variavel_name, unidades_variavel_name,&
                                                       comentario_title, rodape_1, rodape_2, rodape_3               
        integer                                     :: ngrelha, nmercator, nregiao, nmohid, ntipo
        integer                                     :: nnormaliza, LOGARITHM_COLOR_SCALE
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call UnitsManager(Me%Unit%List, CLOSE_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalEnd - ModuleHdf2Surfer - ERR10'
        endif    	 
        
		write(Me%Unit%Evol,'(A)') '_________________________________________________________________________________________________'
	
		write(Me%Unit%Evol,'(A, 5(f12.4,1x),16x, 2(f12.4,1x),4x, i6,A)') '      Average: ',Me%Statistic%GlobalParam%rcorr,&
		                                                                                 Me%Statistic%GlobalParam%rcorr_quad, &
		                                                                                 Me%Statistic%GlobalParam%bias, &
		                                                                                 Me%Statistic%GlobalParam%rmse, &
		                                                                                 Me%Statistic%GlobalParam%z_fisher, &
		                                                                                 Me%Statistic%GlobalParam%Am, &
		                                                                                 Me%Statistic%GlobalParam%Bm, Me%Nout, ' days'
        

        call UnitsManager(Me%Unit%Evol, CLOSE_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalEnd - ModuleHdf2Surfer - ERR20'
        endif    	 
        

        call UnitsManager(Me%Unit%r, CLOSE_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalEnd - ModuleHdf2Surfer - ERR30'
        endif    	 
        

        call UnitsManager(Me%Unit%rr, CLOSE_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalEnd - ModuleHdf2Surfer - ERR40'
        endif    	 
                

        call UnitsManager(Me%Unit%bias, CLOSE_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalEnd - ModuleHdf2Surfer - ERR50'
        endif    	 
        

        call UnitsManager(Me%Unit%rmse, CLOSE_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalEnd - ModuleHdf2Surfer - ERR60'
        endif    	 
            


        call UnitsManager(Me%Unit%surf, OPEN_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalEnd - ModuleHdf2Surfer - ERR10'
        endif    	 

    
		open(Me%Unit%surf,file='run_surfer.dat',status='unknown')
		
		ngrelha   = 0
		nmercator = 0
		nregiao   = 0
		nmohid    = 1
		ntipo     = 1
		
		produto_name           = trim(Me%ScalarPropOut)
		unidades_produto_name  = trim(Me%ScalarPropUnits)
		variavel_name          = trim(Me%ScalarPropOut)
		unidades_variavel_name = trim(Me%ScalarPropUnits)
		comentario_title       = "(Microwave + Infra-red)(*)"
        rodape_1               = "(*): Microwave OI SST data are produced by Remote Sensing Systems and sponsored by National Oceanographic Partnership Program (NOPP),"
        rodape_2               = "the NASA Earth Science Physical  Oceanography Program, and the NASA REASoN DISCOVER Project."
        rodape_3               = "Data are available at www.remss.com."
		
		nnormaliza             = 0
		LOGARITHM_COLOR_SCALE  = 0

		write(Me%Unit%surf,'(A,A,i1,A,i1,A,i1,A,i1)')trim(Me%run),',',ngrelha,',',nmercator,',',nregiao,',',nmohid
		write(Me%Unit%surf,'(A)') adjustl(trim(Me%regiao))
		write(Me%Unit%surf,*)Me%WorkSize2D%JUB
		write(Me%Unit%surf,*)Me%WorkSize2D%IUB
		write(Me%Unit%surf,*) ntipo
		write(Me%Unit%surf,'(A)') adjustl(trim(produto_name))
		write(Me%Unit%surf,'(A)') adjustl(trim(unidades_produto_name))

		write(Me%Unit%surf,'(A)') adjustl(trim(variavel_name)) 
		write(Me%Unit%surf,'(A)') adjustl(trim(unidades_variavel_name))

		write(Me%Unit%surf,'(A)') adjustl(trim(comentario_title))
		write(Me%Unit%surf,'(A)') adjustl(trim(rodape_1))
		write(Me%Unit%surf,'(A)') adjustl(trim(rodape_2))
		write(Me%Unit%surf,'(A)') adjustl(trim(rodape_3))
		write(Me%Unit%surf,*) nnormaliza
		write(Me%Unit%surf,*) LOGARITHM_COLOR_SCALE


        call UnitsManager(Me%Unit%surf, CLOSE_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteGlobalEnd - ModuleHdf2Surfer - ERR10'
        endif    	 

    
    end subroutine WriteGlobalEnd
    
    !---------------------------------------------------------------------------        
    
    subroutine WriteBathymetry    
    
        !Arguments--------------------------------------------------------------
        real,       dimension(:,:), pointer     :: Bat
        
        !Local------------------------------------------------------------------
        integer                                 :: STAT_CALL, i, j
        
        !Begin------------------------------------------------------------------
    
        !Horizontal Grid Data - Water Column (Bathymetry)
        call GetGridData      (GridDataID       = Me%ObjBathymetry,                     &
                               GridData2D       = Bat,                                  &
                               STAT             = STAT_CALL)
                               
        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteBathymetry - ModuleHdf2Surfer - ERR10'
        endif    
    
	    where(bat.lt.0.) bat=0
	    bat=-bat
	    
        call UnitsManager(Me%Unit%bat, OPEN_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteBathymetry - ModuleHdf2Surfer - ERR20'
        endif    	    

	    open(Me%Unit%bat,file='batimetria_'//adjustl(trim(Me%regiao))//'.dat', status='unknown')
	    
	    do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
	    do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
		    write(Me%Unit%bat,*) Me%long(i,j), Me%lat(i,j), bat(i,j)
	    enddo
	    enddo
	    
        call UnitsManager(Me%Unit%bat, CLOSE_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteBathymetry - ModuleHdf2Surfer - ERR20'
        endif    

        call UnGetGridData    (GridDataID       = Me%ObjBathymetry,                     &
                               Array            = Bat,                                  &
                               STAT             = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteBathymetry - ModuleHdf2Surfer - ERR30'
        endif    

    end subroutine WriteBathymetry
    
    !---------------------------------------------------------------------------    
    
    subroutine ModifyHdfInstant(it) 

        !Arguments-------------------------------------------------------------
        integer                             :: it

        !Local-----------------------------------------------------------------

	    real(SP), dimension(:), allocatable :: A, B
	    integer(I4B)                        :: npt	    
	    real(SP)                            :: rcorr, rcorr_quad, bias, rmse, z_fisher   
	    real(SP)                            :: alfa, beta_1, Am, Bm
        real,    dimension(6    ), target   :: AuxTime
        real,    dimension(:    ), pointer  :: TimePtr	 
        integer                             :: STAT_CALL   
	    
		
        !Begin-----------------------------------------------------------------
        
        allocate(A(1:Me%NPointsout))
        allocate(B(1:Me%NPointsout))        
        
        A(1:Me%NPointsout) = Me%PropAout(1:Me%NPointsout)
        B(1:Me%NPointsout) = Me%PropBout(1:Me%NPointsout)        
    
		npt         = Me%NPointsout


		call estatistica(A(1:Me%NPointsout),                                            &
		                 B(1:Me%NPointsout),                                            &
		                 npt,                                                           &
		                 rcorr,                                                         &
		                 rcorr_quad,                                                    &
		                 bias,                                                          &
		                 rmse,                                                          &
		                 z_fisher,                                                      &
		                 alfa,                                                          &
		                 beta_1,                                                        &
		                 Am,                                                            &
		                 Bm)
		
		Me%Statistic%GlobalParam%rcorr         = Me%Statistic%GlobalParam%rcorr      +  &
		                                         rcorr        / real(Me%Nout)
		Me%Statistic%GlobalParam%rcorr_quad    = Me%Statistic%GlobalParam%rcorr_quad +  &
		                                         rcorr_quad   / real(Me%Nout)
		Me%Statistic%GlobalParam%bias          = Me%Statistic%GlobalParam%bias       +  &
		                                         bias         / real(Me%Nout)
		Me%Statistic%GlobalParam%rmse          = Me%Statistic%GlobalParam%rmse       +  &
		                                         rmse         / real(Me%Nout)
		Me%Statistic%GlobalParam%z_fisher      = Me%Statistic%GlobalParam%z_fisher   +  &
		                                         z_fisher     / real(Me%Nout)
		Me%Statistic%GlobalParam%Am            = Me%Statistic%GlobalParam%Am         +  &
		                                         Am           / real(Me%Nout)
		Me%Statistic%GlobalParam%Bm            = Me%Statistic%GlobalParam%Bm         +  &
		                                         Bm           / real(Me%Nout)		
		                                         
		Me%Statistic%Param(it)%rcorr           = rcorr      
		Me%Statistic%Param(it)%rcorr_quad      = rcorr_quad 
		Me%Statistic%Param(it)%bias            = bias       
		Me%Statistic%Param(it)%rmse            = rmse       
		Me%Statistic%Param(it)%z_fisher        = z_fisher   
		Me%Statistic%Param(it)%alfa            = alfa       
		Me%Statistic%Param(it)%beta_1          = beta_1     
		Me%Statistic%Param(it)%Am              = Am         
		Me%Statistic%Param(it)%Bm		       = Bm                                       

        
        if (Me%SurferOut) then
            call WriteInstant(it)
        endif
        
        if (Me%HDFOut   ) then  
            !Writes current time
            call ExtractDate   (Me%CurrentTime, AuxTime(1), AuxTime(2), AuxTime(3),     &
                                AuxTime(4), AuxTime(5), AuxTime(6))
                                
            TimePtr => AuxTime
            
            call HDF5SetLimits  (HDF5ID       = Me%ObjHDF5Out,                          &
                                 ILB          = 1,                                      &
                                 IUB          = 6,                                      &
                                 STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'ModifyHdfInstant - ModuleHdf2Surfer - ERR10'
            endif

            call HDF5WriteData  (HDF5ID       = Me%ObjHDF5Out,                          &
                                 GroupName    = "/Time",                                &
                                 Name         = "Time",                                 &    
                                 Units        = "YYYY/MM/DD HH:MM:SS",                  &
                                 Array1D      = TimePtr,                                &
                                 OutputNumber = it,                                     &
                                 STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'ModifyHdfInstant - ModuleHdf2Surfer - ERR20'
            endif
              
            call WriteInstantHDF(it)
            call WriteInstantTS (it)
        endif            
        
        deallocate(A)
        deallocate(B)

    end subroutine ModifyHdfInstant

    !---------------------------------------------------------------------------
    
    subroutine WriteInstantHDF(it)    
    
        !Arguments--------------------------------------------------------------
        integer                                 :: it
        
        !Local------------------------------------------------------------------
        real,   dimension(:,:,:), pointer       :: Array3D_A, Array3D_B, Array3D_Dif
        real,   dimension(:,:  ), pointer       :: Array2D_A, Array2D_B, Array2D_Dif
        character(Len=StringLength)             :: GroupNameA, GroupNameB, GroupNameDif
        character(Len=StringLength)             :: Name, Units
        integer                                 :: STAT_CALL, i, j, nbln, p, k
        integer                                             :: WorkILB, WorkIUB
        integer                                             :: WorkJLB, WorkJUB
        integer                                             :: WorkKLB, WorkKUB
        integer                                             :: HDF5_CREATE

        
        !Begin------------------------------------------------------------------
        
        if (Me%Out3D) then

            WorkJLB = Me%WorkSize3D%JLB
            WorkJUB = Me%WorkSize3D%JUB
            WorkILB = Me%WorkSize3D%ILB
            WorkIUB = Me%WorkSize3D%IUB           
            WorkKLB = Me%WorkSize3D%KLB
            WorkKUB = Me%WorkSize3D%KUB           
        
        else
        
            WorkJLB = Me%WorkSize2D%JLB
            WorkJUB = Me%WorkSize2D%JUB
            WorkILB = Me%WorkSize2D%ILB
            WorkIUB = Me%WorkSize2D%IUB           
    
        endif


        Name    = Me%ScalarProperty
        Units   = Me%ScalarPropUnits
        
        p = 0
        
        if (Me%Out3D) then
        
            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5Out,                                        &
                                  Me%WorkSize3D%ILB, Me%WorkSize3D%IUB,                 &
                                  Me%WorkSize3D%JLB, Me%WorkSize3D%JUB,                 &
                                  Me%WorkSize3D%KLB, Me%WorkSize3D%KUB,                 &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR10'
            endif
                    
        
            allocate(Array3D_A  (Me%Size3D%ILB:Me%Size3D%IUB,                           &
                                 Me%Size3D%JLB:Me%Size3D%JUB,                           &
                                 Me%Size3D%KLB:Me%Size3D%KUB))
            allocate(Array3D_B  (Me%Size3D%ILB:Me%Size3D%IUB,                           &
                                 Me%Size3D%JLB:Me%Size3D%JUB,                           &
                                 Me%Size3D%KLB:Me%Size3D%KUB))
            allocate(Array3D_Dif(Me%Size3D%ILB:Me%Size3D%IUB,                           &
                                 Me%Size3D%JLB:Me%Size3D%JUB,                           &
                                 Me%Size3D%KLB:Me%Size3D%KUB))
        
	        do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
	        do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
	        do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
	            if (Me%WaterPoints3D(i,j,k) == WaterPoint) then
	                p = p +1 
	                if (.not.(Me%NoDataA(p).or. Me%NoDataB(p))) then

	                    Array3D_A  (i, j, k) = Me%PropA(p)
	                    Array3D_B  (i, j, k) = Me%PropB(p)
	                    Array3D_Dif(i, j, k) = Me%PropA(p)- Me%PropB(p)
	                    
		            else
	                    Array3D_A  (i, j, k) = null_real 
	                    Array3D_B  (i, j, k) = null_real 
	                    Array3D_Dif(i, j, k) = null_real
		            endif
		        else
                    Array3D_A  (i, j, k) = null_real 
                    Array3D_B  (i, j, k) = null_real 
                    Array3D_Dif(i, j, k) = null_real
		        endif
	        enddo
	        enddo
	        enddo
	        
            call HDF5WriteData (HDF5ID       = Me%ObjHDF5Out,                           & 
                                GroupName    = "/Results/"//"Model/"//trim(Name),       &
                                Name         = trim(Name),                              &
                                Units        = trim(Units),                             &
                                Array3D      = Array3D_A,                               &
                                OutputNumber = it,                                      &
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR20'
            endif
    
            
            call HDF5WriteData (HDF5ID       = Me%ObjHDF5Out,                           & 
                                GroupName    = "/Results/"//"Observations/"//trim(Name),&
                                Name         = trim(Name),                              &
                                Units        = trim(Units),                             &
                                Array3D      = Array3D_B,                               &
                                OutputNumber = it,                                      &                                
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR30'
            endif

            call HDF5WriteData (HDF5ID       = Me%ObjHDF5Out,                           & 
                                GroupName    = "/Results/"//"Differences/"//trim(Name), &
                                Name         = trim(Name),                              &
                                Units        = trim(Units),                             &
                                Array3D      = Array3D_Dif,                             &
                                OutputNumber = it,                                      &                                
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR40'
            endif
                      
                
            deallocate(Array3D_A  )
            deallocate(Array3D_B  )
            deallocate(Array3D_Dif)

	                    
        else

            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5Out, WorkILB, WorkIUB, WorkJLB,             &
                                  WorkJUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR50'
            endif
        
        
            allocate(Array2D_A  (Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))
            allocate(Array2D_B  (Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))
            allocate(Array2D_Dif(Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))              
        
	        do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
	        do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
	            if (Me%WaterPoints2D(i,j) == WaterPoint) then
	                p = p +1 
	                if (.not.(Me%NoDataA(p).or. Me%NoDataB(p))) then

	                    Array2D_A  (i, j) = Me%PropA(p)
	                    Array2D_B  (i, j) = Me%PropB(p)
	                    Array2D_Dif(i, j) = Me%PropA(p)- Me%PropB(p)
	                    
		            else
	                    Array2D_A  (i, j) = null_real 
	                    Array2D_B  (i, j) = null_real 
	                    Array2D_Dif(i, j) = null_real
		            endif
		        else
                    Array2D_A  (i, j) = null_real 
                    Array2D_B  (i, j) = null_real 
                    Array2D_Dif(i, j) = null_real
		        endif
	        enddo
	        enddo

            call HDF5WriteData (HDF5ID       = Me%ObjHDF5Out,                           & 
                                GroupName    = "/Results/"//"Model/"//trim(Name),       &
                                Name         = trim(Name),                              &
                                Units        = trim(Units),                             &
                                Array2D      = Array2D_A,                               &
                                OutputNumber = it,                                      &                                                                
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR60'
            endif
    
            
            call HDF5WriteData (HDF5ID       = Me%ObjHDF5Out,                           & 
                                GroupName    = "/Results/"//"Observations/"//trim(Name),&
                                Name         = trim(Name),                              &
                                Units        = trim(Units),                             &
                                Array2D      = Array2D_B,                               &
                                OutputNumber = it,                                      &                                                                
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR70'
            endif                      

            call HDF5WriteData (HDF5ID       = Me%ObjHDF5Out,                           & 
                                GroupName    = "/Results/"//"Differences/"//trim(Name), &
                                Name         = trim(Name),                              &
                                Units        = trim(Units),                             &
                                Array2D      = Array2D_Dif,                             &
                                OutputNumber = it,                                      &                                                                
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR80'
            endif                      
                
	        
            deallocate(Array2D_A  )
            deallocate(Array2D_B  )
            deallocate(Array2D_Dif)
	        
        endif        
       
        
	    
    end subroutine WriteInstantHDF           


    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    
    subroutine HdfStatisticAnalysis
    
        !Arguments--------------------------------------------------------------
        integer                                 :: it
        
        !Local------------------------------------------------------------------
        real,   dimension(:,:,:,:), pointer     :: Array4D_A, Array4D_B
        real,   dimension(:,:,:  ), pointer     :: Array3D_A, Array3D_B, Array3D
        real,   dimension(:,:,:  ), pointer     :: Array3D_R, Array3D_Bias, Array3D_RMSE
        real,   dimension(:,:    ), pointer     :: Array2D_R, Array2D_Bias, Array2D_RMSE        
        real,   dimension(:,:    ), pointer     :: Array2D        
        character(Len=StringLength)             :: GroupNameA, GroupNameB
        character(Len=StringLength)             :: Name, Units
	    real(SP), dimension(:), allocatable     :: A, B
	    integer(I4B)                            :: npt	    
	    real(SP)                                :: rcorr, rcorr_quad, bias, rmse, z_fisher   
	    real(SP)                                :: alfa, beta_1, Am, Bm
        real,    dimension(6    ), target       :: AuxTime
        real,    dimension(:    ), pointer      :: TimePtr	 
        integer                                 :: STAT_CALL, HDF5_READWRITE        
        integer                                 :: i, j, k, n  
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
         
		
        !Begin-----------------------------------------------------------------
        
        if (Me%Out3D) then

            JLB = Me%WorkSize3D%JLB
            JUB = Me%WorkSize3D%JUB
            ILB = Me%WorkSize3D%ILB
            IUB = Me%WorkSize3D%IUB           
            KLB = Me%WorkSize3D%KLB
            KUB = Me%WorkSize3D%KUB           
        
        else
        
            JLB = Me%WorkSize2D%JLB
            JUB = Me%WorkSize2D%JUB
            ILB = Me%WorkSize2D%ILB
            IUB = Me%WorkSize2D%IUB           
    
        endif

        
        
        allocate(A(1:Me%Nout))
        allocate(B(1:Me%Nout))        
        
        Name    = Me%ScalarProperty
        Units   = Me%ScalarPropUnits
        
        call GetHDF5FileAccess  (HDF5_READWRITE = HDF5_READWRITE)
        
        Me%ObjHDF5Out = 0

        !Opens HDF5 File
        call ConstructHDF5      (Me%ObjHDF5Out, trim(Me%FileOutHDF),                    &
                                 HDF5_READWRITE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGlobalIntialHDF - ModuleHdf2Surfer - ERR10'
        
        
        if (Me%Out3D) then
        
            allocate(Array4D_A   (Me%Size3D%ILB:Me%Size3D%IUB,                          &
                                  Me%Size3D%JLB:Me%Size3D%JUB,                          &
                                  Me%Size3D%KLB:Me%Size3D%KUB,                          &
                                  1:Me%Nout))
            allocate(Array4D_B   (Me%Size3D%ILB:Me%Size3D%IUB,                          &
                                  Me%Size3D%JLB:Me%Size3D%JUB,                          &
                                  Me%Size3D%KLB:Me%Size3D%KUB,                          &
                                  1:Me%Nout))

            allocate(Array3D_R   (Me%Size3D%ILB:Me%Size3D%IUB,                          &
                                  Me%Size3D%JLB:Me%Size3D%JUB,                          &
                                  Me%Size3D%KLB:Me%Size3D%KUB))

            allocate(Array3D_Bias(Me%Size3D%ILB:Me%Size3D%IUB,                          &
                                  Me%Size3D%JLB:Me%Size3D%JUB,                          &
                                  Me%Size3D%KLB:Me%Size3D%KUB))

            allocate(Array3D_RMSE(Me%Size3D%ILB:Me%Size3D%IUB,                          &
                                  Me%Size3D%JLB:Me%Size3D%JUB,                          &
                                  Me%Size3D%KLB:Me%Size3D%KUB))
                                  
            allocate(Array3D     (Me%Size3D%ILB:Me%Size3D%IUB,                          &
                                  Me%Size3D%JLB:Me%Size3D%JUB,                          &
                                  Me%Size3D%KLB:Me%Size3D%KUB))
                                  
                                  
            call HDF5SetLimits  (Me%ObjHDF5Out, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues3D - ModuleField4D - ERR30'
                                  
        
	        do n=1,Me%Nout
	        
	            

                call HDF5ReadData  (HDF5ID       = Me%ObjHDF5Out,                           & 
                                    GroupName    = "/Results/"//"Model/"//trim(Name),       &
                                    Name         = trim(Name),                              &
                                    Array3D      = Array3D,                                 &
                                    OutputNumber = n,                                       &
                                    STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR10'
                endif         
                
                Array4D_A(:,:,:,n) = Array3D(:,:,:)
                
                call HDF5ReadData  (HDF5ID       = Me%ObjHDF5Out,                           & 
                                    GroupName    = "/Results/"//"Observations/"//trim(Name),&
                                    Name         = trim(Name),                              &
                                    Array3D      = Array3D,                                 &
                                    OutputNumber = n,                                       &
                                    STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR20'
                endif                      
                
	            Array4D_B(:,:,:,n) = Array3D(:,:,:)
                
            enddo
            
	        do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
	        do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
	        do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
	        
	            if (Me%WaterPoints3D(i,j,k) == WaterPoint) then

                    A(1:Me%Nout) = Array4D_A(i,j,k,1:Me%Nout)
                    B(1:Me%Nout) = Array4D_B(i,j,k,1:Me%Nout)        

	                call estatistica(A(1:Me%Nout),                                      &
	                                 B(1:Me%Nout),                                      &
	                                 Me%Nout,                                           &
	                                 rcorr,                                             &
	                                 rcorr_quad,                                        &
	                                 bias,                                              &
	                                 rmse,                                              &
	                                 z_fisher,                                          &
	                                 alfa,                                              &
	                                 beta_1,                                            &
	                                 Am,                                                &
	                                 Bm)      

                     Array3D_R   (i,j,k) = rcorr
                     Array3D_Bias(i,j,k) = bias
                     Array3D_RMSE(i,j,k) = rmse                                          
	                                  	        
	            endif
	                   
            enddo
	        enddo
	        enddo

            call HDF5SetLimits  (Me%ObjHDF5Out, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadValues3D - ModuleField4D - ERR30'	        

            call HDF5WriteData (HDF5ID       = Me%ObjHDF5Out,                           & 
                                GroupName    = "/Results/"//"R/"//trim(Name),           &
                                Name         = trim(Name),                              &
                                Units        = '-',                                     &
                                Array3D      = Array3D_R,                               &
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR60'
            endif     
	        
            call HDF5WriteData (HDF5ID       = Me%ObjHDF5Out,                           & 
                                GroupName    = "/Results/"//"Bias/"//trim(Name),        &
                                Name         = trim(Name),                              &
                                Units        = trim(Units),                             &
                                Array3D      = Array3D_Bias,                            &
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR60'
            endif     


            call HDF5WriteData (HDF5ID       = Me%ObjHDF5Out,                           & 
                                GroupName    = "/Results/"//"RMSE/"//trim(Name),        &
                                Name         = trim(Name),                              &
                                Units        = trim(Units),                             &
                                Array3D      = Array3D_RMSE,                            &
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR60'
            endif     
            
            deallocate(Array4D_A   )
            deallocate(Array4D_B   )
            
            deallocate(Array3D_R   )
            deallocate(Array3D_Bias)
            deallocate(Array3D_RMSE)
            deallocate(Array3D     )             
	                    
        else

            allocate(Array3D_A   (Me%Size2D%ILB:Me%Size2D%IUB,                          &
                                  Me%Size2D%JLB:Me%Size2D%JUB,                          &
                                  1:Me%Nout))
            allocate(Array3D_B   (Me%Size2D%ILB:Me%Size2D%IUB,                          &
                                  Me%Size2D%JLB:Me%Size2D%JUB,                          &
                                  1:Me%Nout))

            allocate(Array2D_R   (Me%Size2D%ILB:Me%Size2D%IUB,                          &
                                  Me%Size2D%JLB:Me%Size2D%JUB))

            allocate(Array2D_Bias(Me%Size2D%ILB:Me%Size2D%IUB,                          &
                                  Me%Size2D%JLB:Me%Size2D%JUB))

            allocate(Array2D_RMSE(Me%Size2D%ILB:Me%Size2D%IUB,                          &
                                  Me%Size2D%JLB:Me%Size2D%JUB))
                                  
            allocate(Array2D     (Me%Size2D%ILB:Me%Size2D%IUB,                          &
                                  Me%Size2D%JLB:Me%Size2D%JUB))
                                  
            call HDF5SetLimits  (Me%ObjHDF5Out, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues3D - ModuleField4D - ERR30'
                                  
        
	        do n=1,Me%Nout
	        
                call HDF5ReadData  (HDF5ID       = Me%ObjHDF5Out,                           & 
                                    GroupName    = "/Results/"//"Model/"//trim(Name),       &
                                    Name         = trim(Name),                              &
                                    Array2D      = Array2D,                                 &
                                    OutputNumber = n,                                       &
                                    STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR10'
                endif    
                
                Array3D_A(:,:,n) = Array2D(:,:)     
                
                call HDF5ReadData  (HDF5ID       = Me%ObjHDF5Out,                           & 
                                    GroupName    = "/Results/"//"Observations/"//trim(Name),&
                                    Name         = trim(Name),                              &
                                    Array2D      = Array2D,                                 &
                                    OutputNumber = n,                                       &
                                    STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR20'
                endif                      
                
	            Array3D_B(:,:,n) = Array2D(:,:)
                
            enddo
            
	        do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
	        do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
	        
	            if (Me%WaterPoints2D(i,j) == WaterPoint) then

                    A(1:Me%Nout) = Array3D_A(i,j,1:Me%Nout)
                    B(1:Me%Nout) = Array3D_B(i,j,1:Me%Nout)        

	                call estatistica(A(1:Me%Nout),                                      &
	                                 B(1:Me%Nout),                                      &
	                                 Me%Nout,                                           &
	                                 rcorr,                                             &
	                                 rcorr_quad,                                        &
	                                 bias,                                              &
	                                 rmse,                                              &
	                                 z_fisher,                                          &
	                                 alfa,                                              &
	                                 beta_1,                                            &
	                                 Am,                                                &
	                                 Bm)      

                     Array2D_R   (i,j) = rcorr
                     Array2D_Bias(i,j) = bias
                     Array2D_RMSE(i,j) = rmse                                          
	                                  	        
	            endif
	                   
	        enddo
	        enddo

            call HDF5SetLimits  (Me%ObjHDF5Out, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadValues3D - ModuleField4D - ERR30'	        

            call HDF5WriteData (HDF5ID       = Me%ObjHDF5Out,                           & 
                                GroupName    = "/Results/"//"R/"//trim(Name),           &
                                Name         = trim(Name),                              &
                                Units        = '-',                                     &
                                Array2D      = Array2D_R,                               &
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR60'
            endif     
	        
            call HDF5WriteData (HDF5ID       = Me%ObjHDF5Out,                           & 
                                GroupName    = "/Results/"//"Bias/"//trim(Name),        &
                                Name         = trim(Name),                              &
                                Units        = trim(Units),                             &
                                Array2D      = Array2D_Bias,                            &
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR60'
            endif     


            call HDF5WriteData (HDF5ID       = Me%ObjHDF5Out,                           & 
                                GroupName    = "/Results/"//"RMSE/"//trim(Name),        &
                                Name         = trim(Name),                              &
                                Units        = trim(Units),                             &
                                Array2D      = Array2D_RMSE,                            &
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'WriteInstantHDF - ModuleHdf2Surfer - ERR60'
            endif     
            
            deallocate(Array3D_A   )
            deallocate(Array3D_B   )
            
            deallocate(Array2D_R   )
            deallocate(Array2D_Bias)
            deallocate(Array2D_RMSE)
            deallocate(Array2D     )            
	        
        endif        
       
        !Close HDF5 File
        call KillHDF5      (Me%ObjHDF5Out, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGlobalIntialHDF - ModuleHdf2Surfer - ERR10'

        
	    
    end subroutine HdfStatisticAnalysis           


    !---------------------------------------------------------------------------    
    
    subroutine WriteInstantTS(it)    
    
        !Arguments--------------------------------------------------------------
        integer                                 :: it
        
        !Local------------------------------------------------------------------
        real,   dimension(:), pointer           :: DataBuffer
        integer                                 :: STAT_CALL
        
        !Begin------------------------------------------------------------------

        allocate(DataBuffer(1:TSnumber)) 
        
        DataBuffer( 1) = Me%Statistic%Param(it)%rcorr    
        DataBuffer( 2) = Me%Statistic%Param(it)%rcorr_quad
        DataBuffer( 3) = Me%Statistic%Param(it)%bias     
        DataBuffer( 4) = Me%Statistic%Param(it)%rmse     
        DataBuffer( 5) = Me%Statistic%Param(it)%z_fisher
        DataBuffer( 6) = Me%Statistic%Param(it)%alfa
        DataBuffer( 7) = Me%Statistic%Param(it)%beta_1
        DataBuffer( 8) = Me%Statistic%Param(it)%Am
        DataBuffer( 9) = Me%Statistic%Param(it)%Bm
        DataBuffer(10) = Me%Nout

        
        call WriteTimeSerieLine(Me%ObjTimeSerieOut, DataBuffer, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteInstantTS - ModuleHdf2Surfer - ERR10'
        endif                      
                
        deallocate(DataBuffer) 

    end subroutine WriteInstantTS
    
    !---------------------------------------------------------------------------
    
    subroutine WriteInstant(it)    
    
        !Arguments--------------------------------------------------------------
        integer                                 :: it
        real,       dimension(:,:), pointer     :: Bat
        
        !Local------------------------------------------------------------------
        character(Len=StringLength)             :: nome_recta_fit, Aux
        character(Len=PathLength)               :: FileName        
        real                                    :: Pmax, PmaxA, PmaxB
        real                                    :: Pmin, PminA, PminB   
        real                                    :: xi, xdelta, yy
        integer                                 :: STAT_CALL, i, j, nbln, p, k
        
        !Begin------------------------------------------------------------------

        call UnitsManager(Me%Unit%field, OPEN_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteInstant - ModuleHdf2Surfer - ER10'
        endif    	    
        
       

        !FileName = trim(Me%ScalarProperty)//"_"//trim(DateName(Me%CurrentTime))//".dat"
        FileName = "SST"//"_"//trim(DateName(Me%CurrentTime))//".dat"
        
	    open(Me%Unit%field, file=trim(FileName), status='unknown') 
	      
	    write(Me%Unit%List,'(A)') trim(FileName)  ! ficheiro com a listagem dos nomes dos ficheiros
        
        p = 0
        
        if (Me%Out3D) then
        
	        do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
	        do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
	        do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
	            if (Me%WaterPoints3D(i,j,k) == WaterPoint) then
	                p = p +1 
	                if (.not.(Me%NoDataA(p).or. Me%NoDataB(p))) then
		                write(Me%Unit%field,*) Me%long(i,j), Me%lat(i,j), Me%PropA(p), Me%PropB(p), Me%PropB(p)- Me%PropA(p)
		            else
		                write(Me%Unit%field,*) Me%long(i,j), Me%lat(i,j), NoData, NoData, NoData
		            endif
		        else
		            write(Me%Unit%field,*) Me%long(i,j), Me%lat(i,j), NoData, NoData, NoData
		        endif
	        enddo
	        enddo
	        enddo
            
        else
	        do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
	        do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
	            if (Me%WaterPoints2D(i,j) == WaterPoint) then
	                p = p +1 
	                if (.not.(Me%NoDataA(p).or. Me%NoDataB(p))) then
		                write(Me%Unit%field,*) Me%long(i,j), Me%lat(i,j), Me%PropA(p), Me%PropB(p), Me%PropB(p)- Me%PropA(p)
		            else
		                write(Me%Unit%field,*) Me%long(i,j), Me%lat(i,j), NoData, NoData, NoData
		            endif
		        else
		            write(Me%Unit%field,*) Me%long(i,j), Me%lat(i,j), NoData, NoData, NoData
		        endif
	        enddo
	        enddo
        endif        
        
	    
        call UnitsManager(Me%Unit%field, CLOSE_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteInstant - ModuleHdf2Surfer - ERR20'
        endif    


		PmaxA = maxval(Me%PropAout)
		PmaxB = maxval(Me%PropBout)
        Pmax  = max   (PmaxA, PmaxB)
        
		PminA = minval(Me%PropAout)
		PminB = minval(Me%PropBout)
        Pmin  = min   (PminA, PminB)


        call UnitsManager(Me%Unit%fit, OPEN_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteInstant - ModuleHdf2Surfer - ERR20'
        endif    	    

	    open(Me%Unit%fit,file='recta_fit_'//trim(DateName(Me%CurrentTime))//'.bln', status='unknown')

		nbln=30

		xdelta=(Pmax-Pmin)/float(nbln)
		nbln=0
		do xi=Pmin, Pmax, xdelta
			nbln=nbln+1
		enddo

		write(Me%Unit%fit,*) nbln,1
			
		do xi=Pmin, Pmax, xdelta
				
			YY= Me%Statistic%Param(it)%alfa + Me%Statistic%Param(it)%beta_1*xi
			write(Me%Unit%fit,*) xi, yy
		enddo

        call UnitsManager(Me%Unit%fit, CLOSE_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteInstant - ModuleHdf2Surfer - ERR20'
        endif    	    
        
        Aux = DateName(Me%CurrentTime)

        write(Me%Unit%Evol,'(i3,1x,A,1x, 5(f12.4,1x), 2(f12.4,1x), 2(f12.4,1x) ,2(i5,1x),i6)') &
                             it, Aux(1:4)//"-"//Aux(5:6)//"-"//Aux(7:8),                &
                             Me%Statistic%Param(it)%rcorr,                              & 
                             Me%Statistic%Param(it)%rcorr_quad,                         &     
                             Me%Statistic%Param(it)%bias,                               &
                             Me%Statistic%Param(it)%rmse,                               &
                             Me%Statistic%Param(it)%z_fisher,                           &        
                             Me%Statistic%Param(it)%alfa,                               &
                             Me%Statistic%Param(it)%beta_1,                             &
                             Me%Statistic%Param(it)%Am,                                 &
                             Me%Statistic%Param(it)%Bm,                                 &
                             Me%WorkSize2D%JUB,                                         &
                             Me%WorkSize2D%IUB,                                         &
                             Me%Nout        
                             
        write(Me%Unit%r,   *) it, Me%Statistic%Param(it)%rcorr   
        write(Me%Unit%rr,  *) it, Me%Statistic%Param(it)%rcorr_quad
        write(Me%Unit%bias,*) it, Me%Statistic%Param(it)%bias
        write(Me%Unit%rmse,*) it, Me%Statistic%Param(it)%rmse


    end subroutine WriteInstant
    
    !---------------------------------------------------------------------------    


    subroutine estatistica(A, B, npt, rcorr, rcorr_quad, bias, rmse,z_fisher, alfa, beta_1,Am,Bm)

	    !Arguments----------------------------------------------------------------------
	    REAL(SP), DIMENSION(:), INTENT(IN) :: A, B
	    INTEGER(I4B), INTENT(IN)           :: npt	    
	    REAL(SP), INTENT(OUT)              :: rcorr, rcorr_quad, bias, rmse,z_fisher, alfa, beta_1,Am,Bm

	    !Local--------------------------------------------------------------------------	    

	    real(SP) :: cov, var_A, var_B, sdev_A, sdev_B, adev_A, adev_B
	    real(SP) :: skew_A, skew_B, curt_A, curt_B, abdev
	    REAL(SP) :: ave,sdev,var, prob, YY, Pmin, Pmax
	    INTEGER(I4B) :: l_min(1),l_max(1),n1,n2, nbln, i
	    character :: nome_parametros*100, nome_recta_fit*100

    !____________________________________________________________________________________________________
    !
    !	A --> variavel observacaoes (satelite)
    !	B --> variavel modelo (mohid)
    !
    !___________________________________________________________________________________________________


        Am = 0.
        Bm = 0.


	    if (npt > 1) then
	    
            call moment_mohid(npt,A,Am,adev_A,sdev_A,var_A,skew_A,curt_A)
            call moment_mohid(npt,B,Bm,adev_B,sdev_B,var_B,skew_B,curt_B)

	        ! calculo do rms e da bias

	        bias=sum((A(:)-B(:)))
	        bias=bias/float(npt)
	        rmse=sqrt( sum( ( (A(:)-B(:) )**2) )/float(npt))
	        
	        if (Am > null_real/1.e4 .and. Bm > null_real/1.e4) then

	            call pearsn(A,B,rcorr,prob,z_fisher)
	            
	            rcorr_quad=rcorr*rcorr

	            call medfit(B,A,alfa,beta_1,abdev) 
            
            else
                rcorr       = null_real
                prob        = null_real
                z_fisher    = null_real
	            rcorr_quad  = null_real
                alfa        = null_real
                beta_1      = null_real
                abdev	    = null_real
            endif

        else
        
            rcorr       = FillValueReal
            rcorr_quad  = FillValueReal
            bias        = FillValueReal
            rmse        = FillValueReal
            z_fisher    = FillValueReal
            alfa        = FillValueReal
            beta_1      = FillValueReal
            Am          = FillValueReal
            Bm          = FillValueReal        

        
        endif

    end subroutine estatistica
    
    !-----------------------------------------------------------------------------------
    
	subroutine moment_mohid(n,data,ave,adev,sdev,var,skew,curt)
	

        !Arguments----------------------------------------------------------------------
	    integer,                INTENT(IN)  :: n	    
	    REAL(SP), DIMENSION(:), INTENT(IN)  :: data
	    REAL(SP),               INTENT(OUT) :: ave,adev,sdev,var,skew,curt
        
        !Local--------------------------------------------------------------------------       
        REAL(SP), DIMENSION(n)              :: p,s
	    REAL(SP)                            :: ep
        
        !Begin--------------------------------------------------------------------------
        	
	    !if (n <= 1) call nrerror('moment: n must be at least 2')
	    ave=sum(data(:))/n
	    s(:)=data(:)-ave
	    ep=sum(s(:))
	    adev=sum(abs(s(:)))/n
	    p(:)=s(:)*s(:)
	    var=sum(p(:))
	    p(:)=p(:)*s(:)
	    skew=sum(p(:))
	    p(:)=p(:)*s(:)
	    curt=sum(p(:))
	    var=(var-ep**2/n)/(n-1)
	    sdev=sqrt(var)
	    if (var /= 0.0) then
		    skew=skew/(n*sdev**3)
		    curt=curt/(n*var**2)-3.0_sp
	    else
		    !call nrerror('moment: no skew or kurtosis when zero variance')
            skew=0.
		    curt=0.
	    end if
	end subroutine moment_mohid    

    !-----------------------------------------------------------------------------------


    character(len=8) function DateName(DateTime)
    
        !Arguments----------------------------------------------------------------------
        type (T_Time)           :: DateTime
        
        !Local--------------------------------------------------------------------------       
         real                   :: Year, Month, Day
        
        !Begin--------------------------------------------------------------------------
        
        call ExtractDate(Time1 = DateTime, Year = Year, Month = Month, Day = Day)
        
        write(DateName(1:4),'(I4)') int(Year )
        if (Month>9) then
            write(DateName(5:6),'(I2)') int(Month)
        else
            write(DateName(5:5),'(I1)') 0
            write(DateName(6:6),'(I1)') int(Month)
        endif
        if (Day>9) then        
            write(DateName(7:8),'(I2)') int(Day)
        else
            write(DateName(7:7),'(I1)') 0
            write(DateName(8:8),'(I1)') int(Day)
        endif
    
    end function DateName

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillHdf2Surfer(ObjHdf2SurferID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjHdf2SurferID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, STAT_CALL, it

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_
        
        
        
        if (Me%Out3D) then
            call UnGetMap(Me%ObjMap, Me%WaterPoints3D, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)  then
                stop "KillHdf2Surfer - ModuleHdf2Surfer - ERR10"
            endif
                
        else
            call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
            
            if (STAT_CALL /= SUCCESS_)  then
                stop "KillHdf2Surfer - ModuleHdf2Surfer - ERR20"
            endif
            
        endif
        
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%Lat, STAT = STAT_CALL)
            
        if (STAT_CALL /= SUCCESS_)  then
            stop "KillHdf2Surfer - ModuleHdf2Surfer - ERR30"
        endif
        
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%Long, STAT = STAT_CALL)
            
        if (STAT_CALL /= SUCCESS_)  then
            stop "KillHdf2Surfer - ModuleHdf2Surfer - ERR40"
        endif        
        
        if (Me%Out3D) then
        
            call UnGetGeometry(Me%ObjGeometry, Me%Depth, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  then
                stop "KillHdf2Surfer - ModuleHdf2Surfer - ERR50"
            endif        
            
        endif            
        
        deallocate(Me%X, Me%Y, Me%Z)
        deallocate(Me%PropA,   Me%PropB)
        deallocate(Me%NodataA, Me%NodataB)
        
        
        do it=1,Me%HDFSolNumber
            call KillField4D(Field4DID     = Me%HDFSolution(it)%InstanceID,             &
                             STAT          = STAT_CALL)
                                  
            if (STAT_CALL /= SUCCESS_) then
                stop "KillHdf2Surfer - ModuleHdf2Surfer - ERR60"
            endif            
        enddo


        if (Me%Out3D) then
            !Kills Map
            call KillMap            (Me%ObjMap,           STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)  then
                stop "KillHdf2Surfer - ModuleHdf2Surfer - ERR70"
            endif

            !Kills Geometry
            call KillGeometry       (Me%ObjGeometry,      STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)  then
                stop "KillHdf2Surfer - ModuleHdf2Surfer - ERR80"
            endif


        endif            

        !Kills HorizontalMap
        call KillHorizontalMap  (Me%ObjHorizontalMap, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "KillHdf2Surfer - ModuleHdf2Surfer - ERR90"
        endif


        !Kills Bathymetry
        call KillGridData       (Me%ObjBathymetry,    STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "KillHdf2Surfer - ModuleHdf2Surfer - ERR100"
        endif


        !Kills HorizontalGrid
        call KillHorizontalGrid (Me%ObjHorizontalGrid,      STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "KillHdf2Surfer - ModuleHdf2Surfer - ERR110"
        endif


        !Kills Compute Time
        call KillComputeTime    (Me%ObjTime,                STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "KillHdf2Surfer - ModuleHdf2Surfer - ERR120"
        endif


        ObjHdf2SurferID = 0
        
        STAT_      = SUCCESS_

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillHdf2Surfer
        

    !------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------

end module ModuleHdf2Surfer
