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
    USE nrtype; USE nr


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
    real, parameter :: nodata=0.1701410E+39

    !Types---------------------------------------------------------------------
    private :: T_Unit
    type       T_Unit 
        integer     :: bat, field, fit, list, r, rr, bias, rmse, evol, surf
    end type   T_Unit

    private :: T_StatParam
    type       T_StatParam
        real        :: rcorr, rcorr_quad, bias, rmse,z_fisher,alfa,beta_1,Am,Bm
    end type   T_StatParam
    
    private :: T_Stat
    type       T_Stat
        type (T_StatParam), dimension(:), pointer :: Param
        type (T_StatParam)                        :: GlobalParam
    end type   T_Stat
    
    
    private :: T_HDFSolution
    type       T_HDFSolution
        integer                                     :: InstanceID
        character(Len = StringLength)               :: Name
        character(Len = PathLength  )               :: FileIn, FileOut
        integer                                     :: MaskDim
    end type  T_HDFSolution


    
    private :: T_Hdf2Surfer
    type       T_Hdf2Surfer
        integer                                     :: InstanceID
        type (T_Size2D)                             :: Size2D, WorkSize2D        
        type (T_Size3D)                             :: Size3D, WorkSize3D
        character(Len = StringLength)               :: Regiao, Run        
        character(Len = PathLength)                 :: BatimOut, GeoOut
        character(Len = PathLength)                 :: SurferOut1, SurferOut2, SurferOut3
        real                                        :: LatDefault, LongDefault
        type (T_HDFSolution), dimension(:), pointer :: HDFSolution
        integer                                     :: HDFSolNumber
        real,    dimension(:, :, :),  pointer       :: Depth
        real,    dimension(:, :   ),  pointer       :: Lat, Long, Bat
        
        real,    dimension(:      ),  pointer       :: X, Y, Z, Xout, Yout, Zout
        real,    dimension(:      ),  pointer       :: PropA, PropB, PropAout, PropBout
        logical, dimension(:      ),  pointer       :: NodataA, NodataB
        integer, dimension(:, :, :),  pointer       :: WaterPoints3D
        integer, dimension(:, :   ),  pointer       :: WaterPoints2D  
        type (T_Stat)                               :: Statistic      
        integer                                     :: NPoints, NPointsout
        type (T_Time)                               :: BeginTime, EndTime
        type (T_Time)                               :: InitialSystemTime, CurrentTime
        real                                        :: DT
        integer                                     :: Nout
        logical                                     :: HDFOut, SurferOut, Out3D
        character(Len = StringLength)               :: ScalarProperty, ScalarPropUnits, ScalarPropOut
        integer                                     :: ScalarPropertyID
        real                                        :: ValueMin, ValueMax        
        type (T_Unit)                               :: Unit
        logical                                     :: Extrapolate
        
        integer                                     :: ObjEnterData
        integer                                     :: ObjTime          
        integer                                     :: ObjBathymetry      
        integer                                     :: ObjHorizontalMap 
        integer                                     :: ObjHorizontalGrid
        integer                                     :: ObjGeometry      
        integer                                     :: ObjMap           
        
        
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
         
        !read output options
        call ConstructOutPut   
            
        !read input options
        call ConstructInPut
       
       
        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
        
        if (STAT_CALL /= SUCCESS_) then
            stop "ConstructHdf2Surfer - ModuleHdf2Surfer - ERR20"
        endif
        
        call Construct1DAnalysis
        

        !Returns ID
        ObjHdf2SurferID          = Me%InstanceID

        STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructHdf2Surfer
 
    !--------------------------------------------------------------------------

    subroutine ConstructOutPut

        !Arguments---------------------------------------------------------------
        !Local-------------------------------------------------------------------
        real,   dimension(:,:), pointer                 :: SurfaceElevation
        real                                            :: aux
        integer                                         :: STAT_CALL, flag, FirstLine, LastLine

        !Begin-------------------------------------------------------------------
        
        
        
        call GetData(   Me%BatimOut,                                                    &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'BATIM_OUT',                                     &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR20"
        endif

        if (flag == 0) then            
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR30"
        endif

        call GetData(   Me%GeoOut,                                                      &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'GEO_OUT',                                       &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR40"
        endif

        if (flag == 0) then            
            Me%Out3D = .false.
        else
            Me%Out3D = .true.
        endif

        call GetData(   Me%BeginTime,                                                   &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'START',                                         &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR50"
        endif

        if (flag == 0) then            
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR60"
        endif        

        call GetData(   Me%EndTime,                                                     &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'END',                                           &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR70"
        endif

        if (flag == 0) then            
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR80"
        endif        

        call GetData(   Me%DT,                                                          &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'DT',                                            &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR90"
        endif

        if (flag == 0) then            
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR100"
        endif     
        
        aux  = (Me%EndTime - Me%BeginTime)/Me%DT
        
        if (aux /= real(int(aux))) then
            write(*,*) "DT must be a submultiple of the analysed period"
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR105"
        endif
           
        Me%Nout = int(aux) + 1
        
        call GetData(   Me%HDFOut,                                                      &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'HDF_OUT_ON',                                    &
                        default      = .false.,                                         &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR110"
        endif

        call GetData(   Me%SurferOut,                                                   &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'SURFER_OUT',                                    &
                        default      = .false.,                                         &
                        STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR120"
        endif

        if (.not. Me%SurferOut .and. .not. Me%HDFOut) then
            write(*,*) 'You need to choose at least one output option Surfer or HDF or both'
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR130"
        endif


        call GetData(   Me%ScalarProperty,                                              &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'SCALAR_PROPERTY',                               &
                        STAT         = STAT_CALL)
        

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR140"
        endif

        if (flag == 0) then            
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR150"
        endif  
        
        Me%ScalarPropertyID = GetPropertyIDNumber (Me%ScalarProperty)
        

        call StartComputeTime(Me%ObjTime, Me%InitialSystemTime,                         &
                              Me%BeginTime, Me%EndTime, Me%DT, VariableDT = .false.,    &
                              STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructOutPut - ModuleHdf2Surfer - ERR160'
        endif
        

        call ConstructHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     DataFile         = Me%BatimOut,                    &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructOutPut - ModuleHdf2Surfer - ERR170'
        endif
        
        call GetLatitudeLongitude(HorizontalGridID   = Me%ObjHorizontalGrid,            &
                                  Latitude           = Me%LatDefault,                   &
                                  Longitude          = Me%LongDefault,                  &
                                  STAT               = STAT_CALL)        

        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructOutPut - ModuleHdf2Surfer - ERR180'
        endif

        
        call GetHorizontalGridSize(Me%ObjHorizontalGrid, Size = Me%Size2D,              &
                                   WorkSize = Me%WorkSize2D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructOutPut - ModuleHdf2Surfer - ERR190'
        endif
        

        !Horizontal Grid Data - Water Column (Bathymetry)
        call ConstructGridData      (GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     TimeID           = Me%ObjTime,                     &
                                     FileName         = Me%BatimOut,                    &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructOutPut - ModuleHdf2Surfer - ERR200'
        endif


        !Horizontal Map
        call ConstructHorizontalMap (HorizontalMapID  = Me%ObjHorizontalMap,            &
                                     GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     ActualTime       = Me%BeginTime,                   &
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructOutPut - ModuleHdf2Surfer - ERR210'
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
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR220'
            endif


            call ConstructMap           (Map_ID           = Me%ObjMap,                  &
                                         GeometryID       = Me%ObjGeometry,             &
                                         HorizontalMapID  = Me%ObjHorizontalMap,        &
                                         TimeID           = Me%ObjTime,                 &
                                         GridDataID       = Me%ObjBathymetry,           &
                                         HorizontalGridID = Me%ObjHorizontalGrid,       &
                                         STAT             = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR230'
            endif
            
            call GetWaterPoints3D(Me%ObjMap, Me%WaterPoints3D, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)  then
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR240'
            endif
            
            allocate(SurfaceElevation(Me%WorkSize2D%ILB:Me%WorkSize2D%IUB,Me%WorkSize2D%JLB:Me%WorkSize2D%JUB))
            
            SurfaceElevation(:,:) = 0.
            
            call ComputeInitialGeometry(GeometryID       = Me%ObjGeometry,              &
                                        SurfaceElevation = SurfaceElevation,            &
                                        WaterPoints3D    = Me%WaterPoints3D,            &
                                        STAT             = STAT_CALL)
                                                 
            if (STAT_CALL /= SUCCESS_)  then
                stop "Construct1DAnalysis - ModuleHdf2Surfer - ERR250"
            endif
            
            deallocate(SurfaceElevation)

            call UnGetMap(Me%ObjMap, Me%WaterPoints3D, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)  then
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR260'
            endif

            call GetGeometrySize(Me%ObjGeometry, Size = Me%Size3D,                      &
                                 WorkSize = Me%WorkSize3D,  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructOutPut - ModuleHdf2Surfer - ERR270'
            endif

        
        endif
        
        allocate(Me%Statistic%Param(Me%Nout))
        
        call GetData(   Me%Regiao,                                                      &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'REGIAO',                                        &
                        STAT         = STAT_CALL)
        

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR280"
        endif

        if (flag == 0)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR290"
        endif        

        call GetData(   Me%Run,                                                         &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'RUN',                                           &
                        STAT         = STAT_CALL)
        

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR300"
        endif

        if (flag == 0)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR310"
        endif        

        call GetData(   Me%ScalarPropOut,                                               &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'SCALAR_PROP_OUT',                               &
                        STAT         = STAT_CALL)
        

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR320"
        endif

        if (flag == 0) then            
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR330"
        endif  

        call GetData(   Me%ScalarPropUnits,                                             &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'SCALAR_PROP_UNITS',                             &
                        STAT         = STAT_CALL)
        

        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR340"
        endif

        if (flag == 0) then            
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR350"
        endif  

        call GetData(   Me%ValueMin,                                                    &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'VALUE_MIN',                                     &
                        default      = FillValueReal,                                   &
                        STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR360"
        endif

        call GetData(   Me%ValueMax,                                                    &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'VALUE_MAX',                                     &
                        default      = -FillValueReal,                                  &
                        STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR370"
        endif

        call GetData(   Me%Extrapolate,                                                 &
                        Me%ObjEnterData, flag,                                          &
                        SearchType   = FromFile,                                        &
                        ClientModule = 'ModuleHdf2Surfer',                              &
                        keyword      = 'EXTRAPOLATE',                                   &
                        default      = .false.,                                         &
                        STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_)  then
            stop "ConstructOutPut - ModuleHdf2Surfer - ERR380"
        endif


        !----------------------------------------------------------------------

    end subroutine ConstructOutPut
 
    !--------------------------------------------------------------------------

    subroutine ConstructInPut

        !Arguments---------------------------------------------------------------
        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL, flag, it, ClientNumber
        character(StringLength)                         :: solution_begin, solution_end
        logical                                         :: PropertyFound

        !Begin-------------------------------------------------------------------
        
        !Number of input solutions = 2
        Me%HDFSolNumber = 2
        
        allocate(Me%HDFSolution(1:Me%HDFSolNumber))
        
        solution_begin = "<beginsolution>"
        solution_end   = "<endsolution>"
        
        do it=1,Me%HDFSolNumber
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                  &
                                        solution_begin, solution_end,                   &
                                        PropertyFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop "ConstructInPut - ModuleHdf2Surfer - ERR10"
            
i1:         if (PropertyFound) then
                call ConstructOneSolution(it)
            else
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
        
        if (Me%HDFOut) then
        
            
            call GetData(   Me%HDFSolution(it)%FileOut,                                 &
                            Me%ObjEnterData, flag,                                      &
                            SearchType   = FromBlock,                                   &
                            ClientModule = 'ModuleHdf2Surfer',                          &
                            keyword      = 'HDF_OUT',                                   &
                            STAT         = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)  then
                stop "ConstructOneSolution - ModuleHdf2Surfer - ERR30"
            endif

            if (flag == 0) then            
                stop "ConstructOneSolution - ModuleHdf2Surfer - ERR40"
            endif          

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
        
            do j=Me%WorkSize3D%JLB,Me%WorkSize3D%JUB
            do i=Me%WorkSize3D%ILB,Me%WorkSize3D%IUB                        
                
                if (Me%WaterPoints2D(i, j) == WaterPoint)  then
                
                    Me%X(p) =   Me%Long (i, j)
                    Me%Y(p) =   Me%Lat  (i, j)   
                    Me%Z(p) =   0.                 

                    p = p + 1
                endif
                
            enddo
            enddo
        
        endif
    
    end  subroutine Construct1DAnalysis   
    
    !----------------------------------------------------------------------------
    
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
        
        call WriteBathymetry
        
        call WriteGlobalIntial
        
        Me%CurrentTime = Me%BeginTime
        
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

            call ModifyField4DXYZ(Field4DID             = Me%HDFSolution(2)%InstanceID, &
                                  PropertyIDNumber      = Me%ScalarPropertyID,          &
                                  CurrentTime           = Me%CurrentTime,               &
                                  X                     = Me%X,                         &
                                  Y                     = Me%Y,                         &
                                  Z                     = Me%Z,                         &
                                  Field                 = Me%PropB,                     &
                                  NoData                = Me%NoDataB,                   &
                                  STAT                  = STAT_CALL)

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
            
            it = it + 1

            deallocate(Me%Xout    )
            deallocate(Me%Yout    )
            deallocate(Me%Zout    )
            deallocate(Me%PropAout)
            deallocate(Me%PropBout)            
        enddo
        
        call WriteGlobalEnd

        STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyHdf2Surfer
    
    !---------------------------------------------------------------------------
    
    subroutine WriteGlobalIntial
    
        !Arguments-------------------------------------------------------------    

        !Local-----------------------------------------------------------------
        character(Len=StringLength)                 :: FileName
        integer                                     :: STAT_CALL


        call UnitsManager(Me%Unit%field, OPEN_FILE, STAT = STAT_CALL)

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
	
		write(Me%Unit%Evol,'(A, 5(f6.3,1x),16x, 2(f6.3,1x),4x, i6,A)') '      Average: ',Me%Statistic%GlobalParam%rcorr,&
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
        integer                                     :: it

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
    
		call estatistica(Me%PropAout(1:Me%NPointsout), Me%PropBout(1:Me%NPointsout),    &
		                 Me%NPointsout,                                                 &
		                 Me%Statistic%Param(it)%rcorr,                                  &
		                 Me%Statistic%Param(it)%rcorr_quad,                             &
		                 Me%Statistic%Param(it)%bias,                                   &
		                 Me%Statistic%Param(it)%rmse,                                   &
		                 Me%Statistic%Param(it)%z_fisher,                               &
		                 Me%Statistic%Param(it)%alfa,                                   &
		                 Me%Statistic%Param(it)%beta_1,                                 &
		                 Me%Statistic%Param(it)%Am,                                     &
		                 Me%Statistic%Param(it)%Bm)
		
		
		Me%Statistic%GlobalParam%rcorr     =Me%Statistic%GlobalParam%rcorr     + Me%Statistic%Param(it)%rcorr     /real(Me%Nout)
		Me%Statistic%GlobalParam%rcorr_quad=Me%Statistic%GlobalParam%rcorr_quad+ Me%Statistic%Param(it)%rcorr_quad/real(Me%Nout)
		Me%Statistic%GlobalParam%bias      =Me%Statistic%GlobalParam%bias      + Me%Statistic%Param(it)%bias      /real(Me%Nout)
		Me%Statistic%GlobalParam%rmse      =Me%Statistic%GlobalParam%rmse      + Me%Statistic%Param(it)%rmse      /real(Me%Nout)
		Me%Statistic%GlobalParam%z_fisher  =Me%Statistic%GlobalParam%z_fisher  + Me%Statistic%Param(it)%z_fisher  /real(Me%Nout)
		Me%Statistic%GlobalParam%Am        =Me%Statistic%GlobalParam%Am        + Me%Statistic%Param(it)%Am        /real(Me%Nout)
		Me%Statistic%GlobalParam%Am        =Me%Statistic%GlobalParam%Am        + Me%Statistic%Param(it)%Am        /real(Me%Nout)		

        call WriteInstant(it)

    end subroutine ModifyHdfInstant

    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    
    subroutine WriteInstant(it)    
    
        !Arguments--------------------------------------------------------------
        integer                                 :: it
        real,       dimension(:,:), pointer     :: Bat
        
        !Local------------------------------------------------------------------
        character(Len=StringLength)             :: nome_recta_fit, FileName, Aux
        real                                    :: Pmax, PmaxA, PmaxB
        real                                    :: Pmin, PminA, PminB   
        real                                    :: xi, xdelta, yy
        integer                                 :: STAT_CALL, i, j, nbln, p, k
        
        !Begin------------------------------------------------------------------

        call UnitsManager(Me%Unit%field, OPEN_FILE, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) then
            stop 'WriteInstant - ModuleHdf2Surfer - ERR20'
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

        write(Me%Unit%Evol,'(i3,1x,A,1x, 5(f6.3,1x), 2(f7.3,1x), 2(f6.3,1x) ,2(i3,1x),i6)') &
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

        USE nrtype; USE nrutil; USE nr

	    REAL(SP), DIMENSION(:), INTENT(IN) :: A, B
	    INTEGER(I4B), INTENT(IN) :: npt	    
	    REAL(SP), INTENT(OUT) :: rcorr_quad, bias, rmse

	    real(SP) :: cov, var_A, var_B, sdev_A, sdev_B, beta_1,alfa, Am,Bm,Adev,skew, curt,abdev
	    REAL(SP) :: ave,sdev,var, Rcorr, prob, z_fisher, YY, Pmin, Pmax
	    INTEGER(I4B) :: l_min(1),l_max(1),n1,n2, nbln, i
	    character :: nome_parametros*100, nome_recta_fit*100

    !____________________________________________________________________________________________________
    !
    !	A --> variavel observacaoes (satelite)
    !	B --> variavel modelo (mohid)
    !
    !___________________________________________________________________________________________________



	    call moment(A,Am,adev,sdev_A,var_A,skew,curt)
	    call moment(B,Bm,adev,sdev_B,var_B,skew,curt)


	    if(.not.(var_B.eq.0.or.var_A.eq.0)) then

	        call pearsn(A,B,Rcorr,prob,z_fisher)

	        call medfit(B,A,alfa,beta_1,abdev) 
	        rcorr_quad=rcorr*rcorr

	        ! calculo do rms e da bias

	        bias=sum((A(:)-B(:)))
	        bias=bias/float(npt)
	        rmse=sqrt( sum( ( (A(:)-B(:) )**2) )/float(npt))
        
        endif

    end subroutine estatistica


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
