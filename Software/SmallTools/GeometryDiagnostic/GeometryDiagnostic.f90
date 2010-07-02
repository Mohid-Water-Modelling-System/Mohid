!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Geometry Bathymetry Diagnostic
! PROGRAM       : GeometryDiagnostic
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : February 2008
! REVISION      : Luis Fernandes, Guillaume Riflet 
! DESCRIPTION   : Program that analyzes the bathymetry 
!                 and returns the bottom layer thickness, 
!                 and the X(U) and Y(V) thickness ratios
!
!------------------------------------------------------------------------------

program GeometryDiagnostic

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleHDF5
    use ModuleHorizontalGrid
    use ModuleGridData
    use ModuleHorizontalMap
    use ModuleMap
    use ModuleGeometry

    implicit none

    integer                             :: ObjBathymetry    = 0
    integer                             :: ObjHorizontalGrid= 0
    integer                             :: ObjHorizontalMap = 0
    integer                             :: ObjGeometry      = 0
    integer                             :: ObjMap           = 0
    integer                             :: ObjHDF5          = 0
    type(T_Size3D)                      :: Size, WorkSize 
    type(T_Size2D)                      :: Size2D, WorkSize2D
    integer                             :: STAT_CALL, i,j,k, HDF5_CREATE
    integer                             :: kbottom
    character(len=StringLength)         :: BathymetryFile, GeometryFile, OutputFileName
    character(len=StringLength)         :: OutputGriDataName
    real, dimension(:,:,:), pointer     :: SZZ, DWZ
    real, dimension(:, :), pointer      :: DZX, DZY
    real, dimension(:,:  ), pointer     :: Bathymetry, SurfaceElevation, LowerLayerThickness
    real, dimension(:,:  ), pointer     :: LowerLayerThicknessGradU, LowerLayerThicknessGradV
    real, dimension(:,:  ), pointer     :: TopoStiffnessU, TopoStiffnessV
    integer, dimension(:,:,:), pointer  :: WaterPoints3D
    integer, dimension(:,:  ), pointer  :: WaterPoints2D
    integer, dimension(:,:  ), pointer  :: KFloor_Z, KFloor_U, KFloor_V
    real                                :: aux

    call ConstructGeometryDiagnostic
    call ModifyGeometryDiagnostic
    call KillGeometryDiagnostic

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructGeometryDiagnostic
        
        call StartUpMohid("GeometryDiagnostic")

        call ReadKeywords

        call GeometryConstructors

        call Allocatables

        call GetVariables

    end subroutine ConstructGeometryDiagnostic
    
    !--------------------------------------------------------------------------

    subroutine ModifyGeometryDiagnostic
        
        call MakeOutVariables

        call WriteOutVariables    
    
    end subroutine ModifyGeometryDiagnostic
    
    !--------------------------------------------------------------------------

    subroutine KillGeometryDiagnostic

        call Deallocatables

    end subroutine KillGeometryDiagnostic

    !--------------------------------------------------------------------------
    
    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        character(PathLength)                       :: DataFile
        integer                                     :: STAT_CALL
        integer                                     :: ObjEnterData = 0
        integer                                     :: FromFile
        integer                                     :: flag

        call ReadFileName('IN_MODEL', DataFile, "GeometryDiagnostic", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - GeometryDiagnostic - ERR01'

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - GeometryDiagnostic - ERR02'

        call GetExtractType     (FromFile = FromFile)

        !Get Bathymetry
        call GetData(BathymetryFile,                                                   &
                     ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='BATIM',                                         &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadKeywords - GeometryDiagnostic - ERR20'

        !Get Geometry
        call GetData(GeometryFile,                                                   &
                     ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='GEOMETRY',                                         &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadKeywords - GeometryDiagnostic - ERR20'

        !Get Output filename
        call GetData(OutputFileName,                                                   &
                     ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='OUT_FILE',                                         &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadKeywords - GeometryDiagnostic - ERR20'        

        !Get Output GridData filename
        call GetData(OutputGriDataName,                                              &
                     ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                       &
                     keyword      ='OUT_GRIDDATA3D',                                   &
                     STAT         = STAT_CALL)    
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadKeywords - GeometryDiagnostic - ERR20'        

    end subroutine ReadKeywords

    !--------------------------------------------------------------------------
    
    subroutine GeometryConstructors

        write(*,*) "Constructing horizontal grid"

        !Horizontal Grid
        call ConstructHorizontalGrid(HorizontalGridID = ObjHorizontalGrid,          &
                                     DataFile         = BathymetryFile,             &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GeometryConstructors - GeometryDiagnostic - ERR90'

        call GetHorizontalGridSize(ObjHorizontalGrid, Size2D, WorkSize2D, STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'GeometryConstructors - GeometryDiagnostic - ERR130'

        write(*,*)"Constructing Bathymetry"

        !Horizontal Grid Data - Water Column (Bathymetry)
        call ConstructGridData      (GridDataID       = ObjBathymetry,              &
                                     HorizontalGridID = ObjHorizontalGrid,          &
                                     FileName         = BathymetryFile,             &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GeometryConstructors - GeometryDiagnostic - ERR100'

        write(*,*)"Constructing horizontal map"

        !Horizontal Map
        call ConstructHorizontalMap (HorizontalMapID  = ObjHorizontalMap,           &
                                    GridDataID       = ObjBathymetry,              &
                                    HorizontalGridID = ObjHorizontalGrid,          &
                                    STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'GeometryConstructors - GeometryDiagnostic - ERR110'

        write(*,*)"Constructing geometry"

        !Geometry - Water Column
        call ConstructGeometry      (GeometryID       = ObjGeometry,                &
                                     GridDataID       = ObjBathymetry,              &
                                     HorizontalGridID = ObjHorizontalGrid,          &
                                     HorizontalMapID  = ObjHorizontalMap,           &
                                     NewDomain        = GeometryFile,               &
                                    STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'GeometryConstructors - GeometryDiagnostic - ERR120'


        write(*,*)"Constructing 3D map"

        !Map - Water Column            
        call ConstructMap           (Map_ID           = ObjMap,                     &
                                     GeometryID       = ObjGeometry,                &
                                     HorizontalMapID  = ObjHorizontalMap,           &
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'GeometryConstructors - GeometryDiagnostic - ERR130'

    end subroutine GeometryConstructors

    !--------------------------------------------------------------------------

    subroutine Allocatables
    
        allocate(SurfaceElevation(WorkSize2D%ILB:WorkSize2D%IUB, WorkSize2D%JLB:WorkSize2D%JUB))
        SurfaceElevation = 0.
    
        allocate(LowerLayerThickness(WorkSize2D%ILB:WorkSize2D%IUB, WorkSize2D%JLB:WorkSize2D%JUB))
        LowerLayerThickness = -99.

        !Going to calculate the horizontal gradient of the LowerLayerThickness.
        allocate(LowerLayerThicknessGradU(WorkSize2D%ILB:WorkSize2D%IUB, WorkSize2D%JLB:WorkSize2D%JUB+1))
        LowerLayerThicknessGradU = -99.

        allocate(LowerLayerThicknessGradV(WorkSize2D%ILB:WorkSize2D%IUB+1, WorkSize2D%JLB:WorkSize2D%JUB))
        LowerLayerThicknessGradV = -99.

        !Going to calculate the topographic-stiffness ratio 
        !from Beckmann1993 and Shchepetkin2003.
        !
        ! r = ( h_i - h_(i+1) ) / ( h_i + h_(i+1) )
        !
        allocate(TopoStiffnessU(WorkSize2D%ILB:WorkSize2D%IUB, WorkSize2D%JLB:WorkSize2D%JUB))
        TopoStiffnessU = -99.

        allocate(TopoStiffnessV(WorkSize2D%ILB:WorkSize2D%IUB, WorkSize2D%JLB:WorkSize2D%JUB))
        TopoStiffnessV = -99.

    end subroutine Allocatables

    !--------------------------------------------------------------------------

    subroutine GetVariables

        write(*,*)"Getting variables"
        
        !This is 2D---------------------------------------------
        call GetHorizontalGrid(ObjHorizontalGrid, DZX = DZX, DZY = DZY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetVariables - GeometryDiagnostic - ERR137'
        
        call GetGridData(ObjBathymetry, Bathymetry, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'GetVariables - GeometryDiagnostic - ERR138'

        call GetWaterPoints2D(ObjHorizontalMap, WaterPoints2D = WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetVariables - GeometryDiagnostic - ERR331'

        !This is 3D---------------------------------------------
        call GetGeometrySize(ObjGeometry, Size, WorkSize, STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'GetVariables - GeometryDiagnostic - ERR330'

        call GetWaterPoints3D(ObjMap, WaterPoints3D = WaterPoints3D, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'GetVariables - GeometryDiagnostic - ERR332'

        !Compute new volume 
        call ComputeInitialGeometry(ObjGeometry, WaterPoints3D,                         &
                                    SurfaceElevation, ContinuesCompute = .false.,       &
                                    NonHydrostatic = .false.,                           & 
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetVariables - GeometryDiagnostic - ERR333'

        call GetGeometryKFloor(ObjGeometry, Z = KFloor_Z, U = KFloor_U, V = KFloor_V, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'GetVariables - GeometryDiagnostic - ERR03'

        call GetGeometryDistances(ObjGeometry, SZZ = SZZ, DWZ = DWZ, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'GetVariables - GeometryDiagnostic - ERR136'

    end subroutine GetVariables

    !--------------------------------------------------------------------------

    subroutine MakeOutVariables

        ! Create the LowerLayerThickness
        do i = WorkSize%ILB,   WorkSize%IUB
        do j = WorkSize%JLB,   WorkSize%JUB

            if(WaterPoints3D(i,j,WorkSize%KUB) == 1)then
                kbottom = kFloor_Z(i, j)
                LowerLayerThickness(i,j) = DWZ(i,j, kbottom)
            end if

        enddo
        enddo

        ! Create the LowerLayerThicknessGradU
        do i = WorkSize%ILB,   WorkSize%IUB
        do j = WorkSize%JLB+1,   WorkSize%JUB

            kbottom = kFloor_U(i, j)
            if( (WaterPoints2D(i,j) == 1) .and. (WaterPoints2D(i,j-1) == 1) ) then

                aux = DWZ(i,j, kbottom)/DWZ(i,j-1, kbottom)
                LowerLayerThicknessGradU(i,j) = min( aux, 1/aux)

            end if

        enddo
        enddo

        ! Create the LowerLayerThicknessGradV
        do i = WorkSize%ILB+1,   WorkSize%IUB
        do j = WorkSize%JLB,   WorkSize%JUB

            kbottom = kFloor_V(i, j)
            if( (WaterPoints2D(i,j) == 1) .and. (WaterPoints2D(i-1,j) == 1) ) then

                aux = DWZ(i,j, kbottom)/DWZ(i-1,j, kbottom)
                LowerLayerThicknessGradV(i,j) = min( aux, 1/aux)

            end if

        enddo
        enddo

        ! Create the TopoStiffnessU
        do i = WorkSize%ILB,   WorkSize%IUB
        do j = WorkSize%JLB+1,   WorkSize%JUB

            kbottom = kFloor_V(i, j)
            if( (WaterPoints2D(i,j) == 1) .and. (WaterPoints2D(i-1,j) == 1) ) then

                aux = DWZ(i,j, kbottom)/DWZ(i-1,j, kbottom)
                LowerLayerThicknessGradV(i,j) = min( aux, 1/aux)

            end if

        enddo
        enddo

    end subroutine MakeOutVariables

    !--------------------------------------------------------------------------

    subroutine WriteOutVariables

        call WriteHDFVariables

        call WriteGridDataVariables

    end subroutine WriteOutVariables

    !--------------------------------------------------------------------------

    subroutine WriteHDFVariables

        write(*,*)"Writing output file..."

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
    
        !Opens HDF5 File
        call ConstructHDF5(ObjHDF5, OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDFGridVariables - GeometryDiagnostic - ERR335'

        call WriteHDFGridVariables

        call WriteHDFResultsVariables

        !Writes everything to disk
        call HDF5FlushMemory ( ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR70'

        !Writes everything to disk
        call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - FillMatrix - ERR71'

    end subroutine WriteHDFVariables

    !--------------------------------------------------------------------------

    subroutine WriteHDFGridVariables

        !2D TypeZ mesh
        call HDF5SetLimits  (ObjHDF5,   WorkSize%ILB,   WorkSize%IUB,       &
                               WorkSize%JLB,   WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDFGridVariables - GeometryDiagnostic - ERR20'

        call HDF5WriteData  (ObjHDF5, "/Grid", "Bathymetry", "m",           &
                              Array2D = Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDFGridVariables - GeometryDiagnostic - ERR30'            

        call WriteHorizontalGrid (ObjHorizontalGrid, ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDFGridVariables - GeometryDiagnostic - ERR40'

        !3D TypeZ / Vertical Face Centered mesh (SZZ)
        call HDF5SetLimits  ( ObjHDF5,   WorkSize%ILB,   WorkSize%IUB,      &
                               WorkSize%JLB,   WorkSize%JUB, WorkSize%KLB-1,   WorkSize%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDFGridVariables - GeometryDiagnostic - ERR50'

        call HDF5WriteData  ( ObjHDF5, "/Grid/VerticalZ", "Vertical",    &
                             "m", Array3D =  SZZ,            &
                             OutputNumber = 1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDFGridVariables - GeometryDiagnostic - ERR05'

        !3D TypeZ / Vertical Centered mesh (Waterpoints)
        call HDF5SetLimits  ( ObjHDF5,   WorkSize%ILB,   WorkSize%IUB,      &
                            WorkSize%JLB,   WorkSize%JUB, WorkSize%KLB,   WorkSize%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDFGridVariables - GeometryDiagnostic - ERR50'

        call HDF5WriteData   ( ObjHDF5, "/Grid", "WaterPoints", "-",        &
                              Array3D =   WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDFGridVariables - GeometryDiagnostic - ERR60'

    end subroutine WriteHDFGridVariables

    !--------------------------------------------------------------------------

    subroutine WriteHDFResultsVariables
        
        ! Write the LowerLayerThickness
        call HDF5WriteData  (ObjHDF5, "/Results", "LowerLayerThickness", "-",           &
                            Array2D = LowerLayerThickness, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDFResultsVariables - GeometryDiagnostic - ERR30'            

        ! Write the LowerLayerThicknessGradU
        call HDF5WriteData  (ObjHDF5, "/Results", "LowerLayerThicknessGradU", "-",           &
                              Array2D = LowerLayerThicknessGradU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDFResultsVariables - GeometryDiagnostic - ERR31'            

        ! Write the LowerLayerThicknessGradV
        call HDF5WriteData  (ObjHDF5, "/Results", "LowerLayerThicknessGradV", "-",           &
                              Array2D = LowerLayerThicknessGradV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDFResultsVariables - GeometryDiagnostic - ERR32'            

        ! Topographic stiffness ratio (Shchepetkin2003, Beckmann1993)
        call HDF5WriteData  (ObjHDF5, "/Results", "TopoStiffnessU", "-",           &
                              Array2D = LowerLayerThicknessGradU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDFResultsVariables - GeometryDiagnostic - ERR31'            

        call HDF5WriteData  (ObjHDF5, "/Results", "TopoStiffnessV", "-",           &
                              Array2D = LowerLayerThicknessGradV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteHDFResultsVariables - GeometryDiagnostic - ERR32'            

    end subroutine WriteHDFResultsVariables

    !--------------------------------------------------------------------------

    subroutine WriteGridDataVariables

        !Local------------------------------------
        integer                     :: OutputUnit

        !Writes LoweLayerThickness to ascii file
        call UnitsManager(OutputUnit, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteDataFile - GeometryDiagnostic - ERR300'

        open(Unit = OutputUnit, File = trim(OutputGriDataName), Status = "unknown")

        call WriteGridDataFile( OutputUnit, trim("Z"), LowerLayerThickness,      0, 0 )

        call WriteGridDataFile( OutputUnit, trim("U"), LowerLayerThicknessGradU, 0, 1 )

        call WriteGridDataFile( OutputUnit, trim("V"), LowerLayerThicknessGradV, 1, 0 )

        call UnitsManager(OutputUnit, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteDataFile - GeometryDiagnostic - ERR301'

    end subroutine WriteGridDataVariables

    !--------------------------------------------------------------------------

    subroutine WriteGridDataFile (OutputUnit, TYPE_ZUV, Variable, IC, JC)

        !Arguments-----------------------------------------
        integer                                 :: OutputUnit
        character (Len=1)                       :: TYPE_ZUV
        real, dimension(:,:), pointer           :: Variable
        integer                                 :: IC, JC

        call WriteDataLine(OutputUnit, "TYPE_ZUV", trim(TYPE_ZUV))

        write(OutputUnit,*)"<BeginGridData3D>"

        ! Create the LowerLayerThickness(Grad(U,V))
        do i = WorkSize%ILB+IC,   WorkSize%IUB
        do j = WorkSize%JLB+JC,   WorkSize%JUB

            if(WaterPoints2D(i,j) == 1 .and. WaterPoints2D(i-IC,j-JC) == 1)then
                write(OutputUnit,*)i, j, Variable(i,j)
            end if

        enddo
        enddo

        write(OutputUnit,*)"<EndGridData3D>"
              
        write(OutputUnit,*)" "

    end subroutine WriteGridDataFile

    !--------------------------------------------------------------------------

    subroutine Deallocatables
    
        deallocate(SurfaceElevation)
    
        deallocate(LowerLayerThickness)

        deallocate(LowerLayerThicknessGradU)

        deallocate(LowerLayerThicknessGradV)

        deallocate(TopoStiffnessU)

        deallocate(TopoStiffnessV)

    end subroutine Deallocatables

    !--------------------------------------------------------------------------

end program GeometryDiagnostic
