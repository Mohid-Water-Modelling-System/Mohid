program TestVerticalCoordinates

    use ModuleGlobalData
    use ModuleEnterData     
    use ModuleTime
    use ModuleHorizontalGrid                                     
    use ModuleGridData  
    use ModuleHorizontalMap
    use ModuleGeometry 
    use ModuleMap
    use ModuleHDF5
    
             
    implicit none


    integer                             :: ObjBathymetry    = 0
    integer                             :: ObjHorizontalGrid= 0
    integer                             :: ObjHorizontalMap = 0
    integer                             :: ObjGeometry      = 0
    integer                             :: ObjMap           = 0
    integer                             :: ObjTime          = 0
    integer                             :: ObjHDF5          = 0
    type(T_Size3D)                      :: Size, WorkSize 
    integer                             :: STAT_CALL, i,j,k
    type(T_Time)                        :: CurrentTime
    character(len=StringLength)         :: BathymetryFile, GeometryFile, OutputFileName
    real, dimension(:,:,:), pointer     :: SZZ, SlopeX, SlopeY
    real, dimension(:,:  ), pointer     :: Bathymetry, SurfaceElevation, DUX, DVY
    integer, dimension(:,:,:), pointer  :: WaterPoints3D
    integer                             :: HDF5_CREATE
    real                                :: MeanOfTheDepths


    call StartUpMohid("Test Vertical Coordinates")


    !Gets the file name of the Bathymetry
    call ReadFileName('IN_BATIM', BathymetryFile, "Bathymetry File", STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR70'

    !Gets the file name of Geometry
    call ReadFileName('IN_GEOMETRY', GeometryFile, "Geometry File", STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR70'
    
    !Gets the output file name
    call ReadFileName('OUT_FILE', OutputFileName, "Output HDF5 File Name", STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR70'

    write(*,*)"ObjHorizontalGrid"

    call SetDate(CurrentTime, 2004, 1, 1, 0, 0, 0)

    call StartComputeTime( ObjTime,  CurrentTime,  CurrentTime,  60.,     &
                           .false., STAT = STAT_CALL)   
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR60'

    !Horizontal Grid
    call ConstructHorizontalGrid(HorizontalGridID = ObjHorizontalGrid,          &
                                 DataFile         = BathymetryFile,             &
                                 STAT             = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR90'


    write(*,*)"ObjBathymetry"

    !Horizontal Grid Data - Water Column (Bathymetry)
    call ConstructGridData      (GridDataID       = ObjBathymetry,              &
                                 HorizontalGridID = ObjHorizontalGrid,          &
                                 FileName         = BathymetryFile,             &
                                 STAT             = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR100'


    write(*,*)"ObjHorizontalMap"

    !Horizontal Map
    call ConstructHorizontalMap (HorizontalMapID  = ObjHorizontalMap,           &
                                 GridDataID       = ObjBathymetry,              &
                                 HorizontalGridID = ObjHorizontalGrid,          &
                                 ActualTime       = CurrentTime,                &
                                 STAT             = STAT_CALL)  
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR110'


    write(*,*)"ObjGeometry"

    !Geometry - Water Column
    call ConstructGeometry      (GeometryID       = ObjGeometry,                &
                                 GridDataID       = ObjBathymetry,              &
                                 HorizontalGridID = ObjHorizontalGrid,          &
                                 HorizontalMapID  = ObjHorizontalMap,           &
                                 ActualTime       = CurrentTime,                &
                                 NewDomain        = GeometryFile,               &
                                 STAT             = STAT_CALL)  
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR120'


    write(*,*)"ObjMap"

    !Map - Water Column            
    call ConstructMap           (Map_ID           = ObjMap,                     &
                                 GeometryID       = ObjGeometry,                &
                                 HorizontalMapID  = ObjHorizontalMap,           &
                                 TimeID           = ObjTime,                    &
                                 STAT             = STAT_CALL)  
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR130'

    write(*,*)"Getting variables"

    call GetGeometrySize(ObjGeometry, Size, WorkSize, STAT_CALL)  
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR130'

    call GetWaterPoints3D(ObjMap, WaterPoints3D = WaterPoints3D, STAT = STAT_CALL)  
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR130'


    allocate(SurfaceElevation(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB))
    SurfaceElevation = 0.

    !Compute new volume 
    call ComputeInitialGeometry(ObjGeometry, WaterPoints3D,                         &
                                SurfaceElevation, .false.,                          &
                                .false.,                                            & 
                                CurrentTime,                                        &
                                STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR130'


    call GetGeometryDistances(ObjGeometry, SZZ = SZZ, STAT = STAT_CALL)  
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR130'

    call GetHorizontalGrid(ObjHorizontalGrid,                                       &
                           DUX = DUX,                                               &
                           DVY = DVY,                                               &
                           STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleWaterProperties - ERR08'



    call GetGridData(ObjBathymetry, Bathymetry, STAT = STAT_CALL)  
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR130'

    write(*,*)"Writing..."

    !Gets File Access Code
    call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
    
    !Opens HDF5 File
    call ConstructHDF5(ObjHDF5, OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR130'


    call HDF5SetLimits  (ObjHDF5,   WorkSize%ILB,   WorkSize%IUB,       &
                           WorkSize%JLB,   WorkSize%JUB, STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR20'

    call HDF5WriteData  (ObjHDF5, "/Grid", "Bathymetry", "m",           &
                          Array2D = Bathymetry, STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR30'            

    call WriteHorizontalGrid (ObjHorizontalGrid, ObjHDF5, STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR40'



    call HDF5SetLimits  ( ObjHDF5,   WorkSize%ILB,   WorkSize%IUB,      &
                           WorkSize%JLB,   WorkSize%JUB, WorkSize%KLB-1,   WorkSize%KUB, STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR50'


    !Writes SZZ
    call HDF5WriteData  ( ObjHDF5, "/Grid/VerticalZ", "Vertical",    &
                         "m", Array3D =  SZZ,            &
                         OutputNumber = 1, STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaterProperties - ERR05'


    call HDF5SetLimits  ( ObjHDF5,   WorkSize%ILB,   WorkSize%IUB,      &
                           WorkSize%JLB,   WorkSize%JUB, WorkSize%KLB,   WorkSize%KUB, STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR50'

    call HDF5WriteData   ( ObjHDF5, "/Grid", "WaterPoints", "-",        &
                          Array3D =   WaterPoints3D,  STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR60'
        

     
    allocate(SlopeX(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, WorkSize%KLB:WorkSize%KUB))
    SlopeX = 0.
    allocate(SlopeY(WorkSize%ILB:WorkSize%IUB, WorkSize%JLB:WorkSize%JUB, WorkSize%KLB:WorkSize%KUB))
    SlopeY = 0.


    do i = WorkSize%ILB,   WorkSize%IUB
    do j = WorkSize%JLB,   WorkSize%JUB-1
    do k = WorkSize%KLB,   WorkSize%KUB

        if(WaterPoints3D(i,j,k) == WaterPoint)then

            MeanOfTheDepths = abs((SZZ(i, j, k-1) + SZZ(i, j + 1, k-1)) /  2.)

            SlopeX(i,j,k)   = abs((SZZ(i, j, k-1) - SZZ(i, j + 1, k-1)) / (2. * MeanOfTheDepths))

        end if


    enddo
    enddo
    enddo


    do i = WorkSize%ILB,   WorkSize%IUB-1
    do j = WorkSize%JLB,   WorkSize%JUB
    do k = WorkSize%KLB,   WorkSize%KUB

        if(WaterPoints3D(i,j,k) == WaterPoint)then

            MeanOfTheDepths = abs((SZZ(i, j, k-1) + SZZ(i+1, j , k-1)) /  2.)

            SlopeY(i,j,k)   = abs((SZZ(i, j, k-1) - SZZ(i+1, j , k-1)) / (2. * MeanOfTheDepths))

        end if

    enddo
    enddo
    enddo

    call HDF5WriteData  (ObjHDF5, "/Results", "SlopeX", "-",           &
                          Array3D = SlopeX, STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR30'            

    call HDF5WriteData  (ObjHDF5, "/Results", "SlopeY", "-",           &
                          Array3D = SlopeY, STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR30'   
    
    !Writes everything to disk
    call HDF5FlushMemory ( ObjHDF5, STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR70'


    !Writes everything to disk
    call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR70'


end program TestVerticalCoordinates
