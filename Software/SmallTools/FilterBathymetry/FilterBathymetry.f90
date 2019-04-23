!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Batim Filter
! PROJECT       : MohidBatimFilter
! MODULE        : -
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2004
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Program to filter Mohid bathymetries
!------------------------------------------------------------------------------

!Data file - default name 'Filter.dat' (must be placed in working directory)

!   IN_BATIM                    : char              -           !Original bathymetry file name
!   OUT_BATIM                   : char              -           !Filtered bathymetry file name
!   FILTER_RADIUS               : int               -           !Number of influence cells to filter bathymetry
!   FACTOR                      : real              -           !Filtering factor
!   SIGMA_SMOOTH                : logical           -           !Check if the user wants to do an automaitc  smoothing 
!                                                               !in a way that the slope of each cell is below a slope parameter
!   SLOPE_PARAMETER             : real             0.2          ! slope parameter in the bibliography this value for sigma models should be lower than 0.2
!       H_MIN                   : real            [-1e9]        !Minimum value to be filter
!       H_MAX                   : real            [ 1e9]        !Maximum value to be filter


program MohidBatimFilter

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData     
    use ModuleHorizontalGrid                                     
    use ModuleGridData  
    use ModuleHorizontalMap
    use ModuleFunctions
             
    implicit none
    
    !Parameters
    !Filter methods
    integer                                     :: SigmaSmooth_     = 1
    integer                                     :: Average_         = 2    
    integer                                     :: Percentile_      = 3
    

    !Instances
    integer                                     :: ObjEnterData         = 0
    integer                                     :: ObjBathymetry        = 0
    integer                                     :: ObjHorizontalGrid    = 0
    integer                                     :: ObjHorizontalMap     = 0

    !Dummy variable time
    type (T_Time  )                             :: CurrentTime

    !Size
    type (T_Size2D)                             :: Size
    type (T_Size2D)                             :: WorkSize

    !Files
    character(len=PathLength)                   :: BathymetryFile
    character(len=PathLength)                   :: NewBathymetryFile
    character(len=PathLength)                   :: DataFile             = 'Filter.dat'

    !Input variables
    integer                                     :: FilterRadius
    real                                        :: Factor
    logical                                     :: SigmaSmooth
    real                                        :: SlopeLimit, Hmin, Hmax
    integer                                     :: FilterMethod
    logical                                     :: PersistentFilter
    real                                        :: PercentileValue

    !Working variables
    integer, dimension(:,:), pointer            :: WaterPoints2D
    real,    dimension(:,:), pointer            :: Bathymetry
    real,    dimension(:,:), pointer            :: NewBathymetry
    real,    dimension(:,:), pointer            :: SlopeAux

    
    !Begin---------------------------------------------------------------------

    call StartUpMohid("Mohid Bathymetry Filter")

    call ReadOptions

    call ConstructDomain

    call AllocateVariables
    
    if (PersistentFilter) then
        call ApplyPersistentFilter
    else
        call ApplyFilter
    endif

    call WriteNewBathymetry

    call DeallocateVariables

    call KillBathymetryFilter

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
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidBatimFilter - ERR10'

        call GetData(BathymetryFile,                                    &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'IN_BATIM',                         &
                     ClientModule = 'MohidBatimFilter',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidBatimFilter - ERR20'

        call GetData(NewBathymetryFile,                                 &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'OUT_BATIM',                        &
                     ClientModule = 'MohidBatimFilter',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidBatimFilter - ERR30'

        call GetData(FilterRadius,                                      &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'FILTER_RADIUS',                    &
                     ClientModule = 'MohidBatimFilter',                 &
                     Default      = 2,                                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidBatimFilter - ERR40'

        call GetData(Factor,                                            &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'FACTOR',                           &
                     Default      = 0.5,                                &
                     ClientModule = 'MohidBatimFilter',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidBatimFilter - ERR50'

        if(Factor < 0. .or. Factor > 1.)then
            write(*,*)'Factor must be defined between 0 and 1'
            stop 'ReadOptions - MohidBatimFilter - ERR60'
        end if


        call GetData(SigmaSmooth,                                       &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'SIGMA_SMOOTH',                     &
                     Default      = .false.,                            &
                     ClientModule = 'MohidBatimFilter',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidBatimFilter - ERR70'
        
        
        
        call GetData(FilterMethod,                                                      &
                     ObjEnterData, iflag   ,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'FILTER_METHOD',                                    &
                     Default      = SigmaSmooth_,                                       &
                     ClientModule = 'MohidBatimFilter',                                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidBatimFilter - ERR71'
        
        if (FilterMethod /= SigmaSmooth_ .and.                                          &
            FilterMethod /= Average_     .and.                                          &
            FilterMethod /= Percentile_) then
            stop 'ReadOptions - MohidBatimFilter - ERR72'
        endif
            
        if (FilterMethod == Percentile_) then    
            call GetData(PercentileValue,                                               &
                         ObjEnterData, iflag   ,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'PERCENTILE_VALUE',                             &
                         Default      = 0.5,                                            &
                         ClientModule = 'MohidBatimFilter',                             &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidBatimFilter - ERR73'            
        endif
        
        if (PercentileValue < 0 .or. PercentileValue > 1) then
            stop 'ReadOptions - MohidBatimFilter - ERR74'            
        endif
    

        call GetData(PersistentFilter,                                                  &
                     ObjEnterData, iflag   ,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PERSISTENT_FILTER',                                &
                     Default      = .true.,                                             &
                     ClientModule = 'MohidBatimFilter',                                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidBatimFilter - ERR75'    
        
        if (SigmaSmooth) then
            FilterMethod     = SigmaSmooth_
            PersistentFilter = .true.
        endif                

        
        !This parameter is equal to dh/(2h) 
        !The ROMS model community says that for sigma vertical discretizations this values should 
        !not greater than 0.2 
        call GetData(SlopeLimit,                                        &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'SLOPE_PARAMETER',                  &
                     Default      = 0.2,                                &
                     ClientModule = 'MohidBatimFilter',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidBatimFilter - ERR80'

        call GetData(Hmin,                                                              &
                     ObjEnterData, iflag,                                               &
                     SearchType   = FromFile,                                           &
                     keyword      ='H_MIN',                                             &
                     ClientModule ='MohidBatimFilter',                                  &
                     default      = -1e9,                                               &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidBatimFilter - ERR90'


        call GetData(Hmax,                                                              &
                     ObjEnterData, iflag,                                               &
                     SearchType   = FromFile,                                           &
                     keyword      ='H_MAX',                                             &
                     ClientModule ='MohidBatimFilter',                                  &
                     default      =  1e9,                                               &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidBatimFilter - ERR100'



        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidBatimFilter - ERR110'

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
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidBatimFilter - ERR01'


        write(*,*)
        write(*,*)"Constructing bathymetry..."
        write(*,*)


        !Horizontal Grid Data - Water Column (Bathymetry)
        call ConstructGridData      (GridDataID       = ObjBathymetry,              &
                                     HorizontalGridID = ObjHorizontalGrid,          &
                                     FileName         = BathymetryFile,             &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidBatimFilter - ERR02'


        write(*,*)
        write(*,*)"Constructing mapping..."
        write(*,*)

        !Horizontal Map
        call ConstructHorizontalMap (HorizontalMapID  = ObjHorizontalMap,           &
                                     GridDataID       = ObjBathymetry,              &
                                     HorizontalGridID = ObjHorizontalGrid,          &
                                     ActualTime       = CurrentTime,                &
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDsotomain - MohidBatimFilter - ERR03'

        write(*,*)
        write(*,*)"Compiling information..."
        write(*,*)

        call GetWaterPoints2D(HorizontalMapID  = ObjHorizontalMap,                  &
                              WaterPoints2D    = WaterPoints2D,                     &
                              STAT             = STAT_CALL)  

        call GetHorizontalGridSize(ObjHorizontalGrid, Size, WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidBatimFilter - ERR04'

        call GetGridData(GridDataID = ObjBathymetry, GridData2D = Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidBatimFilter - ERR05'

        write(*,*)
        write(*,*)"Domain successfully constructed..."
        write(*,*)


    end subroutine ConstructDomain

    !--------------------------------------------------------------------------


    subroutine AllocateVariables

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)"Allocating new bathymetry in memory..."
        write(*,*)

        nullify (NewBathymetry); allocate(NewBathymetry(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        nullify (SlopeAux     ); allocate(SlopeAux     (Size%ILB:Size%IUB, Size%JLB:Size%JUB))        

        NewBathymetry = -99

    end subroutine AllocateVariables

    !--------------------------------------------------------------------------


    subroutine ApplyFilter 
    
        !Arguments-------------------------------------------------------------    

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, NValues
        real                                        :: AuxSum
        real,   dimension(:), pointer               :: AuxBat
        integer                                     :: iw, jw, Counter, iPer, ii, jj

        !Begin-----------------------------------------------------------------
        
        SlopeAux(:,:) = 0.        
    
        if (FilterMethod == SigmaSmooth_) then

            call SLPMIN(Bathymetry,WorkSize%IUB,WorkSize%JUB,SlopeLimit, WaterPoints2D, SlopeAux)        
            
        else

            NewBathymetry(:,:) = Bathymetry(:,:)
        
            NValues = (2* FilterRadius + 1)**2 
            
            allocate (AuxBat (1:NValues))                                
        

            do j = WorkSize%JLB, WorkSize%JUB
            do i = WorkSize%ILB, WorkSize%IUB

                if (WaterPoints2D(i,j) == WaterPoint .and.                                  &
                    Bathymetry(i, j) < Hmax          .and.                                  &
                    Bathymetry(i, j) > Hmin) then

                    Counter = 0
                    AuxSum  = 0.
                    
                    AuxBat(:) = - FillValueReal 

                    do jj = j - FilterRadius, j + FilterRadius
                    do ii = i - FilterRadius, i + FilterRadius
                        
                        iw = max(ii, WorkSize%ILB) 
                        jw = max(jj, WorkSize%JLB)
                        
                        iw = min(iw, WorkSize%IUB) 
                        jw = min(jw, WorkSize%JUB)                        

                        if (WaterPoints2D(iw, jw) == WaterPoint) then
                        
                            if (Bathymetry(iw, jw) < Hmax                   .and.           &
                                Bathymetry(iw, jw) > Hmin) then                                 

                                Counter = Counter + 1
                                if     (FilterMethod == Average_) then                                
                                    AuxSum  = AuxSum  + Bathymetry(iw, jw)
                                elseif (FilterMethod == Percentile_) then         
                                    AuxBat(Counter) =  Bathymetry(iw, jw)
                                endif                                    
                            endif

                        endif

                    enddo
                    enddo
                    if      (FilterMethod == Average_) then                                
                        NewBathymetry(i, j) = Factor * Bathymetry(i, j) + (1. - Factor) * AuxSum / real(Counter)
                    else if (FilterMethod == Percentile_) then         
                        call Insertion_Sort(AuxBat)
                        iPer = int(real(Counter)*PercentileValue)
                        if (iPer > 0) then
                            NewBathymetry(i, j) = AuxBat(iPer)
                        endif                            
                    endif  
                    
                    SlopeAux(i, j) = FillValueReal
                    
                    if (WaterPoints2D(i  ,j+1) == WaterPoint) then
                        SlopeAux(i, j) = max(SlopeAux(i, j), abs(NewBathymetry(i, j) - NewBathymetry(i  , j+1))/abs(NewBathymetry(i, j)))
                    endif                        
                    if (WaterPoints2D(i  ,j-1) == WaterPoint) then
                        SlopeAux(i, j) = max(SlopeAux(i, j), abs(NewBathymetry(i, j) - NewBathymetry(i  , j-1))/abs(NewBathymetry(i, j)))
                    endif 

                    if (WaterPoints2D(i+1,j  ) == WaterPoint) then
                        SlopeAux(i, j) = max(SlopeAux(i, j), abs(NewBathymetry(i, j) - NewBathymetry(i+1, j  ))/abs(NewBathymetry(i, j)))
                    endif                        
                    if (WaterPoints2D(i-1,j  ) == WaterPoint) then
                        SlopeAux(i, j) = max(SlopeAux(i, j), abs(NewBathymetry(i, j) - NewBathymetry(i-1, j  ))/abs(NewBathymetry(i, j)))
                    endif 
                    
                    
                end if

            enddo
            enddo
            
            deallocate (AuxBat)                                  
            
            Bathymetry(:,:) = NewBathymetry(:,:)
            
        endif        

    end subroutine ApplyFilter

    !--------------------------------------------------------------------------

    subroutine ApplyPersistentFilter

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, jj, jj_old
        logical                                     :: GoalAchived

        !Begin-----------------------------------------------------------------
    
        write(*,*)
        write(*,*)"Filtering original bathymetry..."
        write(*,*)

        GoalAchived = .false.
        
        jj_old = FillValueInt
        
        do while (.not. GoalAchived) 

            GoalAchived = .true.

            call ApplyFilter 

            jj = 0

            do j = WorkSize%JLB+1, WorkSize%JUB-1
            do i = WorkSize%ILB+1, WorkSize%IUB-1

                if (WaterPoints2D(i,j) == WaterPoint) then

                    if (SlopeAux(i,j) > SlopeLimit + .001)  then
                        GoalAchived = .false.
                        jj = jj + 1
                    endif

                endif

            enddo
            enddo

            write(*,*) 'The number of cells where the slope is greater than the limit is ',jj
            
            if (jj_old ==jj) exit
            
            jj_old = jj

        enddo

        NewBathymetry(:,:) = Bathymetry(:,:)

    end subroutine ApplyPersistentFilter

    !--------------------------------------------------------------------------

    subroutine WriteNewBathymetry

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)"Writing new bathymetry..."
        write(*,*)

        call WriteGridData (FileName            = trim(NewBathymetryFile),          &
                            COMENT1             = "Bathymetry filtered by Mohid",   &
                            COMENT2             = "*",                              &
                            HorizontalGridID    = ObjHorizontalGrid,                &
                            FillValue           = -99.,                             &
                            Overwrite           = .true.,                           &
                            GridData2D_Real     = NewBathymetry,                    &
                            STAT                = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteNewBathymetry - MohidBatimFilter - ERR01'

        write(*,*)
        write(*,*)"Done..."
        write(*,*)




    end subroutine WriteNewBathymetry

    !--------------------------------------------------------------------------

    subroutine DeallocateVariables

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)"Deallocating variables..."
        write(*,*)

        deallocate(NewBathymetry); nullify (NewBathymetry)
        deallocate(SlopeAux     ); nullify (SlopeAux     )


    end subroutine DeallocateVariables


    !--------------------------------------------------------------------------

    subroutine KillBathymetryFilter

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)"Finishing freeing memory..."
        write(*,*)

        call KillHorizontalMap  (ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillBathymetryFilter - MohidBatimFilter - ERR01'

        call KillGridData       (ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillBathymetryFilter - MohidBatimFilter - ERR02'

        call KillHorizontalGrid (ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillBathymetryFilter - MohidBatimFilter - ERR03'

    end subroutine KillBathymetryFilter


end program MohidBatimFilter
