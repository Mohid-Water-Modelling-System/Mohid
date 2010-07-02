!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Tide Preview
! PROJECT       : TidePreview
! PROGRAM       : MainTidePreview
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : August 2004
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Program to output tidal elevations preview based on gauge files
!
!------------------------------------------------------------------------------

!Data file - default name 'TidePrevInput.dat' (must be placed in working directory)

!   START                       : YYYY MM DD HH MM SS -         !Start time to compute water level
!   END                         : YYYY MM DD HH MM SS -         !End time to compute water level
!   DT                          : real              -           !Time step to compute water level
!   EXPORT_TO_XYZ               : 0/1               0           !Create a XYZ file with gauge locations
!   XYZ_FILE                    : char              -           !Name of XYZ file to be created

!<begintideprev>
!   IN_TIDES                    : char              -           !Path to gauge file
!   OUT_FILE                    : char              -           !Path to output water level time serie file 
!<endtideprev>

program MohidTidePreview

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModuleGauge

    implicit none

    !Time variables
    type(T_Time)                        :: BeginTime, EndTime, CurrentTime
    type(T_Time)                        :: InitialSystemTime, FinalSystemTime
    logical                             :: VariableDT
    real                                :: DT, TotalCPUTime, CPUTime, LastCPUTime, ElapsedSeconds
    integer, dimension(8)               :: F95Time

    !Types
    type T_TidePrev
        integer                         :: ObjGauge             = 0
        character(len=PathLength)       :: GaugeFile
        character(len=PathLength)       :: OutPutFileName
        character(len=StringLength)     :: Name
        integer                         :: OutPutFileUnit
        real                            :: Longitude, Latitude
        real, dimension(:), pointer     :: ReferenceLevel
        type(T_TidePrev), pointer       :: Next
    end type T_TidePrev
    

    integer                             :: ObjTime              = 0
    type(T_TidePrev), pointer           :: FirstTidePrev
    character(len=PathLength)           :: DataFile             = 'TidePrevInput.dat'
    logical                             :: ExportToXYZ          = OFF
    character(len=PathLength)           :: XYZFile


    call ConstructMohidTidePreview
    call ModifyMohidTidePreview
    call KillMohidTidePreview

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidTidePreview

        !Local-----------------------------------------------------------------
        integer                                     :: NGauges, TotalTidePrevs
        integer                                     :: ObjEnterData         = 0
        integer                                     :: iflag, STAT_CALL, ClientNumber
        logical                                     :: BlockFound
        type(T_TidePrev),   pointer                 :: NewTidePrev
        real, dimension(:), pointer                 :: XLocation, YLocation

        !Begin-----------------------------------------------------------------

        nullify(NewTidePrev, FirstTidePrev)
        
        call StartUpMohid("MohidTidePreview")

        call StartCPUTime
        
        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR01'


        call ReadTimeKeyWords   (ObjEnterData, FromFile, BeginTime, EndTime, DT,        &
                                 VariableDT, "MohidTidePreview")
        
        call StartComputeTime(ObjTime, BeginTime, EndTime, DT = DT,                     &
                              VariableDT = VariableDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR02'

        call GetData(ExportToXYZ,                                                       &
                     ObjEnterData,  iflag,                                              &
                     SearchType   = FromFile,                                           &
                     keyword      = 'EXPORT_TO_XYZ',                                    &
                     Default      = OFF,                                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR03'


        call GetData(XYZFile,                                                           &
                     ObjEnterData,  iflag,                                              &
                     SearchType   = FromFile,                                           &
                     keyword      = 'XYZ_FILE',                                         &
                     Default      = 'GaugeLocation.xyz',                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR04'


        TotalTidePrevs = 0

        write(*,*)
        write(*,*)'Reading gauge(s) information'
        
        do 
            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                         &
                                        '<begintideprev>', '<endtideprev>', BlockFound,     &
                                        STAT = STAT_CALL)
            if (STAT_CALL .EQ. SUCCESS_) then
                
                if (BlockFound) then

                    call AddTidePrev (FirstTidePrev, NewTidePrev)
        

                    call GetData(NewTidePrev%Name,                                          &
                                 ObjEnterData,  iflag,                                      &
                                 SearchType   = FromBlock,                                  &
                                 keyword      = 'NAME',                                     &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR05'

                    call GetData(NewTidePrev%OutPutFileName,                                &
                                 ObjEnterData,  iflag,                                      &
                                 SearchType   = FromBlock,                                  &
                                 keyword      = 'OUT_FILE',                                 &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR06'


                    call GetData(NewTidePrev%GaugeFile,                                     &
                                 ObjEnterData,  iflag,                                      &
                                 SearchType   = FromBlock,                                  &
                                 keyword      = 'IN_TIDES',                                 &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR07'


                    call ConstructGauges(GaugeID    = NewTidePrev%ObjGauge,                 &
                                         TimeID     = ObjTime,                              &
                                         GaugeFile  = NewTidePrev%GaugeFile,                &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR08'


                    !Get the number of gauges in use
                    call GetNGauges(NewTidePrev%ObjGauge, NGauges, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR09'
        
                    if(NGauges .ne. 1)then
                        write(*,*)'Can only compute one gauge per tide preview at the time.'
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR10'
                    end if

                    allocate(NewTidePrev%ReferenceLevel(NGauges))

                    !Get the current elevation at the gauges
                    call GetReferenceLevel(NewTidePrev%ObjGauge,                            &
                                           NewTidePrev%ReferenceLevel,                      &
                                           STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR11'

                    allocate(XLocation(NGauges), YLocation(NGauges))

                    call GetGaugeLocation(NewTidePrev%ObjGauge, XLocation, YLocation, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR12'

                    NewTidePrev%Longitude = XLocation(1)
                    NewTidePrev%Latitude  = YLocation(1)
                    
                    TotalTidePrevs = TotalTidePrevs + 1
                else
                    exit     !No more blocks
                end if 
            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then 
                stop 'ConstructMohidTidePreview - MohidTidePreview - ERR13'
            end if
        end do


        if(ExportToXYZ) call ExportGaugeLocations

        call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR14'

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR15'

        write(*,*)
        write(*,*)'Total number of gauges to preview :', TotalTidePrevs

    end subroutine ConstructMohidTidePreview
    
    !--------------------------------------------------------------------------

    subroutine ExportGaugeLocations
        
        !Local-----------------------------------------------------------------
        type(T_TidePrev), pointer                   :: TidePrev
        integer                                     :: STAT_CALL, XYZUnit
        integer                                     :: ID

        !Begin-----------------------------------------------------------------

        ID = 1

        call UnitsManager (XYZUnit, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExportGaugeLocations - MohidTidePreview - ERR01'

        open(XYZUnit, file = trim(XYZFile), Status = 'unknown')

        write(XYZUnit, *)"<begin_xyz>"

        TidePrev => FirstTidePrev
        
        do while(associated(TidePrev))

            write(XYZUnit, *) TidePrev%Longitude, TidePrev%Latitude, ID, trim(TidePrev%Name)

            ID = ID + 1

            TidePrev => TidePrev%Next
        end do

        write(XYZUnit, *)"<end_xyz>"

        call UnitsManager (XYZUnit, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExportGaugeLocations - MohidTidePreview - ERR02'


    end subroutine ExportGaugeLocations
   
    !--------------------------------------------------------------------------

    subroutine ModifyMohidTidePreview
        
        !Local-----------------------------------------------------------------
        logical                                     :: Running
        integer                                     :: STAT_CALL
        real                                        :: TotalTime
        real,       dimension(:), pointer           :: OpenPoints, WaterLevel
        type(T_TidePrev), pointer                   :: TidePrev

        !Begin-----------------------------------------------------------------


        Running      = .true.
        CurrentTime  = BeginTime

        allocate(OpenPoints(1), WaterLevel(1))
        OpenPoints = 1

        
        TidePrev => FirstTidePrev
        
        do while(associated(TidePrev))

            call UnitsManager (TidePrev%OutPutFileUnit, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidTidePreview - MohidTidePreview - ERR01'

            open(TidePrev%OutPutFileUnit, file = TidePrev%OutPutFileName, Status = 'unknown')

            call WriteDataLine(TidePrev%OutPutFileUnit, "SERIE_INITIAL_DATA", CurrentTime)
            write(TidePrev%OutPutFileUnit, *) "TIME_UNITS              : SECONDS"
            write(TidePrev%OutPutFileUnit, *) "LONGITUDE               : ", TidePrev%Longitude
            write(TidePrev%OutPutFileUnit, *) "LATITUDE                : ", TidePrev%Latitude

            write(TidePrev%OutPutFileUnit, *)
            write(TidePrev%OutPutFileUnit, *)
            write(TidePrev%OutPutFileUnit, *) 'Seconds   Elevation'
            write(TidePrev%OutPutFileUnit, *) '<BeginTimeSerie>'


            TidePrev => TidePrev%Next
        end do


        TotalTime = 0.

        do while (Running)

            TidePrev => FirstTidePrev
        
            do while(associated(TidePrev))

                call GaugeLevel(TidePrev%ObjGauge,                          &
                                WaterLevel,                                 &
                                OpenPoints,                                 &
                                CurrentTime,                                &
                                ReferenceLevel = TidePrev%ReferenceLevel,   &
                                STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidTidePreview - MohidTidePreview - ERR02'

                write(TidePrev%OutPutFileUnit, *) TotalTime, WaterLevel(1) + TidePrev%ReferenceLevel

                TidePrev => TidePrev%Next

            end do

            call cpu_time(CPUTime)
            if (CPUTime - LastCPUTime > 10.) then
                LastCPUTime = CPUTime
                call PrintProgress(ObjTime, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidTidePreview - MohidTidePreview - ERR03'
            endif

            
            CurrentTime = CurrentTime + DT
            TotalTime   = TotalTime   + DT

            call ActualizeCurrentTime(ObjTime, DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidTidePreview - MohidTidePreview - ERR04'


            if (abs(CurrentTime - EndTime) > DT / 10.) then
                Running = .true.
            else
                Running = .false.
            endif

        enddo
    



        TidePrev => FirstTidePrev

        do while(associated(TidePrev))

            write(TidePrev%OutPutFileUnit, *) '<EndTimeSerie>'

            call UnitsManager (TidePrev%OutPutFileUnit, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidTidePreview - MohidTidePreview - ERR05'

            TidePrev => TidePrev%Next

        end do

     
    end subroutine ModifyMohidTidePreview
    
    !--------------------------------------------------------------------------

    subroutine KillMohidTidePreview

        call StopCPUTime

        call ShutdownMohid ("MohidTidePreview", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMohidTidePreview
    
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

    subroutine AddTidePrev (FirstTidePrev, ObjTidePrev)

        !Arguments-------------------------------------------------------------
        type (T_TidePrev), pointer                   :: FirstTidePrev
        type (T_TidePrev), pointer                   :: ObjTidePrev

        !Local-----------------------------------------------------------------
        type (T_TidePrev), pointer                   :: NewTidePrev
        type (T_TidePrev), pointer                   :: PreviousTidePrev
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewTidePrev)
        nullify  (NewTidePrev%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstTidePrev)) then
            FirstTidePrev         => NewTidePrev
            ObjTidePrev           => NewTidePrev
        else
            PreviousTidePrev      => FirstTidePrev
            ObjTidePrev           => FirstTidePrev%Next
            do while (associated(ObjTidePrev))
                PreviousTidePrev  => ObjTidePrev
                ObjTidePrev       => ObjTidePrev%Next
            enddo
            ObjTidePrev           => NewTidePrev
            PreviousTidePrev%Next => NewTidePrev
        endif


    end subroutine AddTidePrev

end program MohidTidePreview
