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

    !Working variables
    integer, dimension(:,:), pointer            :: WaterPoints2D
    real,    dimension(:,:), pointer            :: Bathymetry
    real,    dimension(:,:), pointer            :: NewBathymetry

    
    !Begin---------------------------------------------------------------------

    call StartUpMohid("Mohid Bathymetry Filter")

    call ReadOptions

    call ConstructDomain

    call AllocateVariables
	
    if (SigmaSmooth) then
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
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidBatimFilter - ERR03'

        write(*,*)
        write(*,*)"Compiling information..."
        write(*,*)

        call GetWaterPoints2D(HorizontalMapID  = ObjHorizontalMap,                  &
                              WaterPoints2D    = WaterPoints2D,                     &
                              STAT             = STAT_CALL)  

        call GetHorizontalGridSize(ObjHorizontalGrid, Size, WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidBatimFilter - ERR04'

        call GetGridData(ObjBathymetry, Bathymetry, STAT = STAT_CALL)
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

        NewBathymetry = -99

    end subroutine AllocateVariables

    !--------------------------------------------------------------------------


    subroutine ApplyFilter

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        real                                        :: AuxSum
        integer                                     :: iw, jw, Counter

        !Begin-----------------------------------------------------------------
    
        write(*,*)
        write(*,*)"Filtering original bathymetry..."
        write(*,*)

        NewBathymetry(:,:) = Bathymetry(:,:)

        do j = WorkSize%JLB, WorkSize%JUB
        do i = WorkSize%ILB, WorkSize%IUB

            if (WaterPoints2D(i,j) == WaterPoint .and.                                  &
                Bathymetry(i, j) < Hmax          .and.                                  &
                Bathymetry(i, j) > Hmin) then

                Counter = 0
                AuxSum  = 0.

                do jw = j - FilterRadius, j + FilterRadius
                do iw = i - FilterRadius, i + FilterRadius

                    if (jw >= WorkSize%JLB .and. jw <= WorkSize%JUB .and.               &
                        iw >= WorkSize%ILB .and. iw <= WorkSize%IUB .and.               &
                        WaterPoints2D(iw, jw) == WaterPoint) then
                        
                        if (Bathymetry(iw, jw) < Hmax                   .and.           &
                            Bathymetry(iw, jw) > Hmin) then                                 

                            Counter = Counter + 1
                            AuxSum  = AuxSum  + Bathymetry(iw, jw)

                        endif

                    endif

                enddo
                enddo

                NewBathymetry(i, j) = Factor * Bathymetry(i, j) + (1. - Factor) * AuxSum / real(Counter)

            end if

        enddo
        enddo

    end subroutine ApplyFilter

    !--------------------------------------------------------------------------

    subroutine ApplyPersistentFilter

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, jj
        logical                                     :: GoalAchived
        real                                        :: Slope,r1,r2,r3,r4
        real,   dimension(:,:), pointer             :: SlopeAux

        !Begin-----------------------------------------------------------------
    
        write(*,*)
        write(*,*)"Filtering original bathymetry..."
        write(*,*)

        allocate(SlopeAux(Size%ILB:Size%IUB,Size%JLB:Size%JUB))

        GoalAchived = .false.

        do while (.not. GoalAchived) 

            GoalAchived = .true.

            SlopeAux(:,:) = 0.

            call SLPMIN(Bathymetry,WorkSize%IUB,WorkSize%JUB,SlopeLimit, WaterPoints2D, SlopeAux)

            jj = 0

            do j = WorkSize%JLB+1, WorkSize%JUB-1
            do i = WorkSize%ILB+1, WorkSize%IUB-1

                if (WaterPoints2D(i,j) == WaterPoint) then

                    r1 = abs((Bathymetry(i,j)-Bathymetry(i,j-1))/(Bathymetry(i,j)+Bathymetry(i,j-1)))
                    r2 = abs((Bathymetry(i,j)-Bathymetry(i,j+1))/(Bathymetry(i,j)+Bathymetry(i,j+1)))
                    r3 = abs((Bathymetry(i,j)-Bathymetry(i-1,j))/(Bathymetry(i,j)+Bathymetry(i-1,j)))
                    r4 = abs((Bathymetry(i,j)-Bathymetry(i+1,j))/(Bathymetry(i,j)+Bathymetry(i+1,j)))

                    Slope = max(r1,r2,r3,r4)

                    if (SlopeAux(i,j) > SlopeLimit + .001)  then
                        GoalAchived = .false.
                        jj = jj + 1
                    endif

                endif

            enddo
            enddo

            write(*,*) 'The number of cells where the slope is greater than the limit is ',jj

        enddo

        NewBathymetry(:,:) = Bathymetry(:,:)

        deallocate(SlopeAux)

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
