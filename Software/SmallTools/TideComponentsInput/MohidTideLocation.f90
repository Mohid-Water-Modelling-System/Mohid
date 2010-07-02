program MohidTideLocation

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData     
    use ModuleHorizontalGrid                                     
    use ModuleGridData  
    use ModuleHorizontalMap
    
             
    implicit none


    integer                                     :: ObjEnterData     = 0
    integer                                     :: ObjBathymetry    = 0
    integer                                     :: ObjHorizontalGrid= 0
    integer                                     :: ObjHorizontalMap = 0
    type (T_Time  )                             :: CurrentTime
    type (T_Size2D)                             :: Size
    type (T_Size2D)                             :: WorkSize
    character(len=StringLength)                 :: BathymetryFile
    real                                        :: MeanLevel
    integer                                     :: CoordType
    integer, dimension(:, :) , pointer          :: BoundaryPoints2D
    real, dimension(:,:)     , pointer          :: GridLongitude, GridLatitude
    real, dimension(:,:)     , pointer          :: GridLongitudeConn, GridLatitudeConn



    call StartMohidTide

    contains

    subroutine StartMohidTide

        call StartUpMohid("Mohid Tide Points Generator")

        call ReadOptions
		
        call PointGenerator

    end subroutine StartMohidTide

    
    !--------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        character(PathLength)                       :: DataFile     = 'MohidTideAction.dat'
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Running MohidTide...'
        write(*,*)

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidTide - ERR01'

        call GetData(BathymetryFile,                                    &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'IN_BATIM',                         &
                     ClientModule = 'MohidTide',                        &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidTide - ERR02'

        call GetData(MeanLevel,                                         &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'MEAN_LEVEL',                       &
                     ClientModule = 'MohidTide',                        &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidTide - ERR03'

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidTide - ERR04'

    end subroutine ReadOptions

    subroutine PointGenerator

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: i, j

        !Begin-----------------------------------------------------------------
    
        write(*,*)"ObjHorizontalGrid"

        call SetDate(CurrentTime, 2004, 1, 1, 0, 0, 0)

        !Horizontal Grid
        call ConstructHorizontalGrid(HorizontalGridID = ObjHorizontalGrid,          &
                                     DataFile         = BathymetryFile,             &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PointGenerator - MohidTide - ERR01'


        write(*,*)"ObjBathymetry"

        !Horizontal Grid Data - Water Column (Bathymetry)
        call ConstructGridData      (GridDataID       = ObjBathymetry,              &
                                     HorizontalGridID = ObjHorizontalGrid,          &
                                     FileName         = BathymetryFile,             &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PointGenerator - MohidTide - ERR02'


        write(*,*)"ObjHorizontalMap"

        !Horizontal Map
        call ConstructHorizontalMap (HorizontalMapID  = ObjHorizontalMap,           &
                                     GridDataID       = ObjBathymetry,              &
                                     HorizontalGridID = ObjHorizontalGrid,          &
                                     ActualTime       = CurrentTime,                &
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'PointGenerator - MohidTide - ERR03'


        write(*,*)"Getting variables"

        call GetBoundaries(HorizontalMapID  = ObjHorizontalMap,                     &
                           BoundaryPoints2D = BoundaryPoints2D,                     &
                           STAT             = STAT_CALL)  

        call GetGridLatitudeLongitude(HorizontalGridID        = ObjHorizontalGrid,  &
                                            GridLatitude      = GridLatitude ,      &
                                            GridLongitude     = GridLongitude,      &
                                            GridLatitudeConn  = GridLatitudeConn,   &
                                            GridLongitudeConn = GridLongitudeConn,  &
                                            STAT              = STAT_CALL )

        call GetHorizontalGridSize(ObjHorizontalGrid, Size, WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PointGenerator - MohidTide - ERR04'

        call GetGridCoordType(ObjHorizontalGrid, CoordType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PointGenerator - MohidTide - ERR05'

        open(10,file="input.dat")
        do i = WorkSize%ILB, WorkSize%IUB
            do j= WorkSize%JLB, WorkSize%JUB
                if (BoundaryPoints2D(i,j) == 1) then

                    select case (CoordType)

                        case (4)

                            write(10,'(4f10.4,f6.2,I3)') GridLongitude(i,j), GridLatitude(i,j),   &
                                                         GridLongitude(i,j), GridLatitude(i,j),   &
                                                         MeanLevel, 0

                        case default

                      end select

                end if
            enddo
        enddo
        close(10)

    end subroutine PointGenerator

end program MohidTideLocation
