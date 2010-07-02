!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : ReservoirVolumes
! PROGRAM       : MainReservoirVolumes
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig  v4.0
! DESCRIPTION   : To calculate the volume balance of a reservoir
!
!------------------------------------------------------------------------------

program ReservoirVolumes

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModuleDischarges
    use ModuleHorizontalGrid
    use ModuleGridData

    implicit none

    type(T_Time)                        :: BeginTime, EndTime, CurrentTime
    real                                :: DT
    logical                             :: VariableDT
    type (T_Time)                       :: InitialSystemTime, FinalSystemTime
    real                                :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)               :: F95Time
                                        
    integer                             :: ObjTime              = 0
    integer                             :: ObjHorizontalGrid    = 0
    integer                             :: ObjGridData          = 0
    integer                             :: ObjDischarges        = 0
                                        
    real                                :: RefLevel
    logical                             :: AdjustVolumesOn
    real(8)                             :: CurrentArea, CurrentVolume
    real                                :: CurrentLevel
    character(PathLength)               :: BathymetryFile
    character(PathLength)               :: DischargesFile
    character(PathLength)               :: ObservedLevelsFile
    character(len=*), parameter         :: ProjectFilePath  = "ReservoirVolumes.dat"
    integer                             :: UnitFlow, UnitLevel


    type T_AdjustVolumes
        type (T_Time), dimension(:), pointer :: Times
        real, dimension(:), pointer          :: Level
    end type T_AdjustVolumes

    type (T_AdjustVolumes)              :: AdjustVolumes

    !ExternalVar
    type (T_Size2D)                     :: Size, WorkSize
    real, dimension(:,:), pointer       :: GridData
    real, dimension(:,:), pointer       :: GridCellArea
    integer, dimension(:,:), pointer    :: DefineCellsMap
    real                                :: MaximumValue, MinimumValue, FillValue
    real                                :: Year, Month, Day, Hour, Minute, Second
    type (T_Time)                       :: TestDate

    call ConstructReservoirVolumes
    call ModifyReservoirVolumes
    call KillReservoirVolumes

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructReservoirVolumes

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ObjAdjustVolumes
        integer                                     :: ClientNumber    
        logical                                     :: BlockFound
        integer                                     :: StartLine, EndLine
        real, dimension(:), pointer                 :: DataLine
        integer                                     :: count, flag, iLine

        call StartUpMohid("Mohid Reservoir Volumes")

        call StartCPUTime

        call ReadKeywords

        call StartComputeTime(ObjTime, BeginTime, EndTime, DT, VariableDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirVolumes - ReservoirVolumes - ERR01'

        call ConstructHorizontalGrid(ObjHorizontalGrid, BathymetryFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirVolumes - ReservoirVolumes - ERR02'

        call ConstructGridData      (ObjGridData, ObjHorizontalGrid, ObjTime, BathymetryFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirVolumes - ReservoirVolumes - ERR03'

        call Construct_Discharges   (ObjDischarges, ObjTime, DischargesFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirVolumes - ReservoirVolumes - ERR03a'

        !Gets ExternalVar
        call GetGridData            (ObjGridData, GridData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirVolumes - ReservoirVolumes - ERR04'

        call GetMaximumValue        (ObjGridData, MaximumValue, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirVolumes - ReservoirVolumes - ERR05'

        call GetMinimumValue        (ObjGridData, MinimumValue, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirVolumes - ReservoirVolumes - ERR05'

        call GetFillValue           (ObjGridData, FillValue, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirVolumes - ReservoirVolumes - ERR06'

        call GetGridCellArea        (ObjHorizontalGrid, GridCellArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirVolumes - ReservoirVolumes - ERR07'

        call GetDefineCellsMap      (ObjHorizontalGrid, DefineCellsMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirVolumes - ReservoirVolumes - ERR08'

        call GetHorizontalGridSize  (ObjHorizontalGrid, Size, WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirVolumes - ReservoirVolumes - ERR09'

        !Calculates initial Volume 
        call AreaAndVolumeFromLevel(CurrentLevel, CurrentArea, CurrentVolume)

        !Reads Volumes to adjust
        if (AdjustVolumesOn) then

            call ConstructEnterData (ObjAdjustVolumes, ObservedLevelsFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirVolumes - ReservoirVolumes - ERR20'

            call ExtractBlockFromBuffer(ObjAdjustVolumes, ClientNumber,                  &
                                        "<begin_level>", "<end_level>",                  &
                                        BlockFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirVolumes - ReservoirVolumes - ERR21'

            call GetBlockSize(ObjAdjustVolumes, ClientNumber, StartLine, EndLine, FromBlock_, STAT = STAT_CALL)

            allocate (AdjustVolumes%Times       (EndLine-StartLine+1-2))
            allocate (AdjustVolumes%Level       (EndLine-StartLine+1-2))

            allocate (DataLine(7))
            count = 1
            do iLine = StartLine + 1, EndLine - 1
                
                !Get a new line
                call GetData(DataLine,                                                   &
                             ObjAdjustVolumes,                                           &
                             flag,                                                       &
                             Buffer_Line = iLine,                                        &
                             STAT = STAT_CALL) 

                call SetDate (AdjustVolumes%Times(count), DataLine(1), DataLine(2), DataLine(3), DataLine(4), DataLine(5), DataLine(6))
                AdjustVolumes%Level(count) = DataLine(7)

                count = count+1
            enddo

            call KillEnterData      (ObjAdjustVolumes)
            
            call UnitsManager(UnitFlow, OPEN_FILE, STAT = STAT_CALL)
            open(unit=UnitFlow, status='unknown', file ='AdjustedFlows.dat')

        endif

        call UnitsManager(UnitLevel, OPEN_FILE, STAT = STAT_CALL)
        open(unit=UnitLevel, status='unknown', file ='PredictedLevel.dat')


    end subroutine ConstructReservoirVolumes
    
    !--------------------------------------------------------------------------

    subroutine ModifyReservoirVolumes
        
        !Local-----------------------------------------------------------------
        real                                        :: CPUTime, LastCPUTime = 0.
        logical                                     :: Running
        integer                                     :: STAT_CALL
        integer                                     :: nDischarges, iDis
        real                                        :: Flow
        real(8)                                     :: NewVolume
        real                                        :: Year, Month, Day, Hour, Minute, Second
        integer                                     :: iTime
        real                                        :: ReservoirLevel
        real(8)                                     :: ReservoirVolume, ReservoirArea
        real                                        :: AddFlowIn, AddFlowOut, Temp
        logical :: a

        Running      = .true.
        CurrentTime  = BeginTime

        do while (Running)
            
            CurrentTime = CurrentTime + DT
            call ExtractDate (CurrentTime, Year, Month, Day, Hour, Minute, Second)

            call ActualizeCurrentTime(ObjTime, DT, STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_) stop 'ModifyReservoirVolumes - ERR01'


            !Gets number of discharges
            call GetDischargesNumber (ObjDischarges, nDischarges, STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_) stop 'ModifyReservoirVolumes - ERR02'

            !Calculates new volume
            NewVolume = CurrentVolume
            do iDis = 1, nDischarges
                call GetDischargeWaterFlow (ObjDischarges, CurrentTime, iDis, CurrentLevel, Flow, STAT = STAT_CALL)
                if (STAT_CALL .ne. SUCCESS_) stop 'ModifyReservoirVolumes - ERR03' 
                NewVolume = NewVolume + (Flow * DT)
            enddo

            call IterateLevel          (NewVolume, CurrentVolume, CurrentLevel)
            call AreaAndVolumeFromLevel(CurrentLevel, CurrentArea, CurrentVolume)

            !If Program runs to adjust volumes, see which volume should be in the reservoir
            !end calculate the differnence
            if (AdjustVolumesOn) then
            
                iTime = 1
                do while (AdjustVolumes%Times(iTime) .le. CurrentTime)
                    iTime = iTime + 1
                enddo

                call InterpolateValueInTime(CurrentTime, AdjustVolumes%Times(iTime-1), AdjustVolumes%Level(iTime-1), &
                                            AdjustVolumes%Times(iTime), AdjustVolumes%Level(iTime), ReservoirLevel)
            
        
                call AreaAndVolumeFromLevel(ReservoirLevel, ReservoirArea, ReservoirVolume)

                if (ReservoirVolume > CurrentVolume) then
                    AddFlowIn  = 0
                    AddFlowOut = (ReservoirVolume - CurrentVolume) / DT
                else
                    AddFlowIn  = (ReservoirVolume - CurrentVolume) / DT
                    AddFlowOut = 0
                endif

                NewVolume = CurrentVolume + (AddFlowIn + AddFlowOut) * DT

                call IterateLevel          (NewVolume, CurrentVolume, CurrentLevel)
                call AreaAndVolumeFromLevel(CurrentLevel, CurrentArea, CurrentVolume)

                !To Improve
                select case (int(Month))
                    case (1)
                        Temp = 8.267473118
                    case (2)
                        Temp = 8.006696429
                    case (3)
                        Temp = 12.07490909
                    case (4)
                        Temp = 12.37148818
                    case (5)
                        Temp = 18.41962366
                    case (6)
                        Temp = 21.22236111
                    case (7)
                        Temp = 21.72419355
                    case (8)
                        Temp = 25.62150538
                    case (9)
                        Temp = 21.5375
                    case (10)
                        Temp = 15.8995935
                    case (11)
                        Temp = 12
                    case (12)
                        Temp = 10
                end select

                write(unit=UnitFlow, fmt=100)int(Year), int(Month), int(Day), int(Hour), int(Minute), int(Second), &
                                         AddFlowIn, AddFlowOut, Temp     
 100            format (1X, i4, 1x, i2, 1x, i2, 1x, i2, 1x, i2, 1x, i2, 1x, f12.4, 1x, f12.4, 1x, f8.2)

            endif

            write(unit=UnitLevel, fmt=110)int(Year), int(Month), int(Day), int(Hour), int(Minute), int(Second), &
                                    CurrentLevel, CurrentArea, CurrentVolume
 110        format (1X, i4, 1x, i2, 1x, i2, 1x, i2, 1x, i2, 1x, i2, 1x, f8.2, 1x, f14.2, 1x, f14.2)


            if (abs(CurrentTime - EndTime) > DT / 10.) then
                Running = .true.
            else
                Running = .false.
            endif

            call CPU_TIME(CPUTime)
            if (CPUTime - LastCPUTime > 10.) then
                LastCPUTime = CPUTime
                call PrintProgress(ObjTime, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyReservoirVolumes - ERR99'
            endif


        enddo
    
    end subroutine ModifyReservoirVolumes
    
    !--------------------------------------------------------------------------

    subroutine IterateLevel(NewVolume, CurrentVolume, NewLevel)

        !Arguments-------------------------------------------------------------
        real(8), intent(IN)                         :: NewVolume, CurrentVolume
        real, intent(OUT)                           :: NewLevel

        !Local-----------------------------------------------------------------
        real                                        :: Level1, Level2
        real(8)                                     :: Volume1, Volume2
        real(8)                                     :: Area1, Area2

        !Actualizes Current Level
        if (NewVolume > CurrentVolume) then
            Level1 = CurrentLevel
            Level2 = RefLevel - MinimumValue
        else
            Level1 = RefLevel - MaximumValue
            Level2 = CurrentLevel
        endif

        call AreaAndVolumeFromLevel(Level1, Area1, Volume1)
        call AreaAndVolumeFromLevel(Level2, Area2, Volume2)

        do while (abs(Level1-Level2) > 1e-6)    !This value must be quiet small
            
            if (NewVolume - Volume1 > Volume2 - NewVolume) then
                Level1 = ((Level1 + Level2) / 2.0 + Level1) / 2.0
                call AreaAndVolumeFromLevel(Level1, Area1, Volume1)
            else
                Level2 = ((Level1 + Level2) / 2.0 + Level2) / 2.0
                call AreaAndVolumeFromLevel(Level2, Area2, Volume2)
            endif

        enddo

        NewLevel = (Level1 + Level2) / 2.0

    end subroutine IterateLevel

    !--------------------------------------------------------------------------

    subroutine KillReservoirVolumes

        !Closes files
        if (AdjustVolumesOn) call UnitsManager(UnitFlow, CLOSE_FILE)
        call UnitsManager(UnitLevel, CLOSE_FILE)

        call StopCPUTime

        call ShutdownMohid ("ReservoirVolumes", ElapsedSeconds, TotalCPUTime)

    end subroutine KillReservoirVolumes
    
    !--------------------------------------------------------------------------

    subroutine StartCPUTime

        call date_and_time(Values = F95Time)
        
        call SetDate      (InitialSystemTime, float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)

    end subroutine StartCPUTime
    
    !--------------------------------------------------------------------------

    subroutine StopCPUTime

        call date_and_time(Values = F95Time)
        
        call SetDate      (FinalSystemTime,   float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)
        
        call cpu_time(TotalCPUTime)

        ElapsedSeconds = FinalSystemTime - InitialSystemTime

    end subroutine StopCPUTime
    
    !--------------------------------------------------------------------------
    
    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ObjEnterData = 0
        integer                                     :: flag

        call ConstructEnterData(ObjEnterData, ProjectFilePath, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ReservoirVolumes - ERR00'


        !File Name of the Reservoir Grid Data File
        call GetData    (BathymetryFile, ObjEnterData, flag,                     &
                         keyword      = 'IN_GRID_DATA',                          &
                         ClientModule = 'ModuleGridData',                        &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ReservoirVolumes - ERR01'

        !File Name of the Discharges Data File
        call GetData    (DischargesFile, ObjEnterData, flag,                     &
                         keyword      = 'DISCHARG',                              &
                         ClientModule = 'ModuleGridData',                        &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ReservoirVolumes - ERR02'

        call ReadTimeKeyWords   (ObjEnterData, FromFile, BeginTime, EndTime, DT, &
                                 VariableDT, "ReservoirVolumes")
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ReservoirVolumes - ERR03'

        !Gets Initial Level
        call GetData    (CurrentLevel, ObjEnterData, flag,                       &
                         keyword      = 'INITIAL_LEVEL',                         &
                         ClientModule = 'ModuleGridData',                        &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ReservoirVolumes - ERR04'

        !Gets RefLevel
        call GetData    (RefLevel, ObjEnterData, flag,                           &
                         keyword      = 'REFERENCE_LEVEL',                       &
                         ClientModule = 'ModuleGridData',                        &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ReservoirVolumes - ERR05'

        !Gets the option to adjust exiting data
        call GetData    (AdjustVolumesOn, ObjEnterData, flag,                    &
                         keyword      = 'ADJUST_VOLUMES',                        &
                         ClientModule = 'ModuleGridData',                        &
                         default      = .false.,                                 &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ReservoirVolumes - ERR06'

        if (AdjustVolumesOn) then

            !File Name Data File with the observed Levels
            call GetData (ObservedLevelsFile, ObjEnterData, flag,                &
                          keyword      = 'IN_OBSERVED_LEVELS',                   &
                          ClientModule = 'ModuleGridData',                       &
                          STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ReservoirVolumes - ERR07'
            
        endif


        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ReservoirVolumes - ERR03'

    end subroutine ReadKeywords

    !--------------------------------------------------------------------------

    subroutine AreaAndVolumeFromLevel (Level, Area, Volume)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: Level
        real(8), intent(OUT)                        :: Area, Volume

        !Local-----------------------------------------------------------------
        integer                                     :: i, j

        Volume = 0.
        Area   = 0.

        do j = WorkSize%JLB, WorkSize%JUB
        do i = WorkSize%ILB, WorkSize%IUB

            if (GridData(i, j) > FillValue / 2.0 .and. DefineCellsMap(i, j) == 1) then
                if ((RefLevel - GridData(i, j)) < Level) then
                    Area   = Area   + GridCellArea(i, j)
                    Volume = Volume + (Level - (RefLevel - GridData(i, j))) * GridCellArea(i, j)
                endif
            endif

        enddo
        enddo

    end subroutine AreaAndVolumeFromLevel

end program ReservoirVolumes
