!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : ConvertGridDataToHDF5
! PROGRAM       : ConvertGridDataToHDF5
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2010
! REVISION      : Eduardo Jauch
! DESCRIPTION   : ConvertGridDataToHDF5 to convert Grid Data files to MOHID HDF5 format
!
!------------------------------------------------------------------------------

program ConvertGridDataToHDF5

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleHorizontalGrid
    use ModuleGridData
    use ModuleTime
    use ModuleEnterData
	
    implicit none

    integer                              :: ObjHDF5           = 0
    integer                              :: ObjHorizontalGrid = 0
    integer                              :: ObjGridData       = 0
    
    integer                              :: HDF5_READWRITE
    integer                              :: HDF5_READ
    integer                              :: HDF5_CREATE
    
    type (T_Size2D)                      :: Size
    type (T_Size2D)                      :: WorkSize
    real, dimension(:,:), pointer        :: Grid

    integer                              :: STAT_CALL
    
    type (T_Time)                        :: InitialSystemTime, FinalSystemTime
    integer, dimension(8)                :: F95Time    
    real, dimension(6), target           :: AuxTime
    real, dimension(:), pointer          :: TimePointer
    type (T_Time)                        :: Time
    character(1024)                      :: InputFile, OutputFile
    real                                 :: ElapsedSeconds, TotalCPUTime

    !----------------------------------------------------------------------------------------------------------

    call date_and_time(Values = F95Time)
    call StartUpMohid("ConvertGridDataToHDF5")
    call SetDate(InitialSystemTime, float(F95Time(1)), float(F95Time(2)), &
                                    float(F95Time(3)), float(F95Time(5)), &
                                    float(F95Time(6)), float(F95Time(7))+ &
                                    F95Time(8)/1000.)
    
    call ReadDataFile    
    call CreateHDF5File

    call date_and_time(Values = F95Time)
    call SetDate (FinalSystemTime, float(F95Time(1)), float(F95Time(2)), &
                                   float(F95Time(3)), float(F95Time(5)), &
                                   float(F95Time(6)), float(F95Time(7))+ &
                                   F95Time(8)/1000.)
    call cpu_time(TotalCPUTime)
    ElapsedSeconds = FinalSystemTime - InitialSystemTime
    call ShutdownMohid ("ConvertGridDataToHDF5", ElapsedSeconds, TotalCPUTime)

    contains
    
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
        
    subroutine ReadDataFile
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer                                     :: ObjEnterData = 0     

        !Begin-----------------------------------------------------------------
        call ConstructEnterData (ObjEnterData, "ConvertGridDataToHDF5.dat", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ConvertGridDataToHDF5 - ERR010'

        call GetData(Time,                                   &
                     ObjEnterData, iflag,                    &
                     SearchType   = FromFile,                &
                     keyword      = 'TIME',                  &
                     ClientModule = 'ConvertGridDataToHDF5', &
                     STAT         = STAT_CALL)
        if ((STAT_CALL .NE. SUCCESS_) .OR. (iflag .EQ. 0)) stop 'ReadDataFile - ConvertGridDataToHDF5 - ERR020'
        
        call GetData(InputFile,                              &
                     ObjEnterData, iflag,                    &
                     SearchType   = FromFile,                &
                     keyword      = 'INPUT_FILE',            &
                     default      = 'GridData.dat',          &                     
                     ClientModule = 'ConvertGridDataToHDF5', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - ConvertGridDataToHDF5 - ERR030'
        
        call GetData(OutputFile,                             &
                     ObjEnterData, iflag,                    &
                     SearchType   = FromFile,                &
                     keyword      = 'OUTPUT_FILE',           &
                     default      = 'HDF.hdf5',              &                     
                     ClientModule = 'ConvertGridDataToHDF5', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - ConvertGridDataToHDF5 - ERR040'
    
        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ConvertGridDataToHDF5 - ERR050'
        !End-------------------------------------------------------------------
    
    end subroutine ReadDataFile

    !--------------------------------------------------------------------------

    subroutine CreateHDF5File

        !Local-----------------------------------------------------------------
        integer,   dimension(:,:), pointer:: MapingPoints
        
        !Begin-----------------------------------------------------------------
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE,      &
                                 HDF5_READ = HDF5_READ,          &
                                 HDF5_READWRITE = HDF5_READWRITE)
        
        !Constructs Horizontal Grid
        call ConstructHorizontalGrid(ObjHorizontalGrid, &
                                     trim(InputFile),   &
                                     STAT = STAT_CALL)           
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR010'

        !Constructs GridData
        call ConstructGridData(ObjGridData,                &
                               ObjHorizontalGrid,          &
                               FileName = trim(InputFile), &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR020'

        call GetHorizontalGridSize(ObjHorizontalGrid, &
                                   Size,              &
                                   WorkSize,          &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR030'
                                       
        !Gets a pointer to Grid
        call GetGridData (ObjGridData, Grid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR040'

        call ConstructHDF5(ObjHDF5,          &
                           trim(OutputFile), &
                           HDF5_CREATE,      &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR050'                       

        !Write the Horizontal Grid
        call WriteHorizontalGrid(ObjHorizontalGrid, &
                                 ObjHDF5,           &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR060'


        call ExtractDate(Time, AuxTime(1), AuxTime(2), &
                               AuxTime(3), AuxTime(4), &
                               AuxTime(5), AuxTime(6))
        TimePointer => AuxTime

        call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR070'

        call HDF5WriteData  (ObjHDF5, "/Time", "Time",    &
                             "YYYY/MM/DD HH:MM:SS",       &
                             Array1D      = TimePointer,  &
                             OutputNumber = 1,            &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR080'

        call HDF5SetLimits(ObjHDF5,             &
                           WorkSize%ILB,        &
                           WorkSize%IUB,        &
                           WorkSize%JLB,        &
                           WorkSize%JUB,        &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR090'

        call HDF5WriteData(ObjHDF5,                     &
                           "//Results/leaf area index", &
                           "leaf area index",           &
                           "m",                         &                           
                           Array2D = Grid,              &
                           OutputNumber = 1,            &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR100'

        allocate(MapingPoints(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        MapingPoints = 1
        call HDF5WriteData   (ObjHDF5, "/Grid", "MapingPoints", "-",      &
                              Array2D = MapingPoints, STAT = STAT_CALL)                              
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR110'

        call HDF5SetLimits(ObjHDF5,             &
                           WorkSize%ILB,        &
                           WorkSize%IUB,        &
                           WorkSize%JLB,        &
                           WorkSize%JUB,        &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR120'

        call HDF5WriteData(ObjHDF5, "/Grid", "Bathymetry", "m",           &
                           Array2D = MapingPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR130'
        deallocate(MapingPoints)


        call KillHDF5(ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR140'

        call KillGridData(ObjGridData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR150'
        
        call KillHorizontalGrid(ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR160'

        !End-------------------------------------------------------------------
    
    end subroutine CreateHDF5File

end program ConvertGridDataToHDF5


