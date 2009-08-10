program RiverNetwork

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime
    use ModuleDrainageNetwork
    use ModuleFunctions
    use ModuleStopWatch,      only : CreateWatchGroup, KillWatchGroup
    use ModuleTimeSerie

    implicit none


    !Variables-------------------------------------------------------------------
    type(T_Time)                                :: BeginTime, EndTime, CurrentTime
    real                                        :: DT, MaxDT
    logical                                     :: VariableDT
    type (T_Time)                               :: InitialSystemTime, FinalSystemTime
    real                                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)                       :: F95Time

    integer                                     :: ObjDrainageNetwork   = 0
    integer                                     :: ObjTime              = 0    
    integer                                     :: ObjEnterData         = 0
    character(PathLength)                       :: ModelName            = 'RiverNetwork'

    !Atmosphere stuff   
    integer                                     :: ObjTSTopRadiation    = 0
    integer                                     :: ObjTSAirTemperature  = 0
    integer                                     :: ObjTSRelativeHumidity= 0
    integer                                     :: ObjTSCloudCover      = 0
    integer                                     :: ObjTSWindVelocity    = 0

    integer                                     :: ColumnTopRadiation    = 0
    integer                                     :: ColumnAirTemperature  = 0
    integer                                     :: ColumnRelativeHumidity= 0
    integer                                     :: ColumnCloudCover      = 0
    integer                                     :: ColumnWindVelocity    = 0
  
    

    !Begin-----------------------------------------------------------------------  
    
    call ConstructRiverNetwork
    call ModifyRiverNetwork   
    call KillRiverNetwork
    

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructRiverNetwork

         
        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL, flag
        logical                                     :: CheckMass, NeedsRadiation
        logical                                     :: NeedAtmosphere
        character(PathLength)                       :: FileName

        call StartupMohid(trim(ModelName))

        call StartCPUTime

        call ReadKeywords

        !Constructs Time 
        call StartComputeTime (ObjTime, BeginTime, EndTime, DT, VariableDT, MaxDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR01'
        CheckMass = .true. 
        								       
        call ConstructDrainageNetwork (ObjDrainageNetwork, ObjTime, CheckMass = CheckMass, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR02'

        !Check if Radiation is needed
        call GetNeedsRadiation (ObjDrainageNetwork, NeedsRadiation, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR03'

        !Check if Radiation is needed
        call GetNeedsAtmosphere(ObjDrainageNetwork, NeedAtmosphere, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR04'
        
        if (NeedAtmosphere .or. NeedsRadiation) then

            !Solar radiation
            call GetData(FileName,                                                  &   
                         ObjEnterData, flag,                                        &  
                         keyword      = 'SOLARRADIATION_FILE',                      &
                         ClientModule = 'RiverNetwork',                             &
                         SearchType   = FromFile,                                   &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR05'

            if (flag == 0) then
                write (*,*) 'Missing SOLARRADIATION_FILE in Model input file'
                stop 'RiverNetwork - ConstructRiverNetwork - ERR06'
            endif


            call GetData(ColumnTopRadiation,                                        &   
                         ObjEnterData, flag,                                        &  
                         keyword      = 'SOLARRADIATION_COLUMN',                    &
                         ClientModule = 'RiverNetwork',                             &
                         SearchType   = FromFile,                                   &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR07'

            call StartTimeSerieInput(ObjTSTopRadiation, FileName, ObjTime,          &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR08'

        endif


        if (NeedAtmosphere) then

            !AirTemperature
            call GetData(FileName,                                                  &   
                         ObjEnterData, flag,                                        &  
                         keyword      = 'AIR_TEMPERATURE_FILE',                     &
                         ClientModule = 'RiverNetwork',                             &
                         SearchType   = FromFile,                                   &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR09'

            if (flag == 0) then
                write (*,*) 'Missing AIR_TEMPERATURE_FILE in Model input file'
                stop 'RiverNetwork - ConstructRiverNetwork - ERR10'
            endif

            call GetData(ColumnAirTemperature,                                      &   
                         ObjEnterData, flag,                                        &  
                         keyword      = 'AIR_TEMPERATURE_COLUMN',                   &
                         ClientModule = 'RiverNetwork',                             &
                         SearchType   = FromFile,                                   &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR11'

            call StartTimeSerieInput(ObjTSAirTemperature, FileName, ObjTime,        &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR12'

            !CloudCover
            call GetData(FileName,                                                  &   
                         ObjEnterData, flag,                                        &  
                         keyword      = 'CLOUD_COVER_FILE',                         &
                         ClientModule = 'RiverNetwork',                             &
                         SearchType   = FromFile,                                   &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR13'

            if (flag == 0) then
                write (*,*) 'Missing CLOUD_COVER_FILE in Model input file'
                stop 'RiverNetwork - ConstructRiverNetwork - ERR14'
            endif

            call GetData(ColumnCloudCover,                                          &   
                         ObjEnterData, flag,                                        &  
                         keyword      = 'CLOUD_COVER_COLUMN',                       &
                         ClientModule = 'RiverNetwork',                             &
                         SearchType   = FromFile,                                   &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR15'

            call StartTimeSerieInput(ObjTSCloudCover, FileName, ObjTime,        &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR16'

            !RelativeHumidity
            call GetData(FileName,                                                  &   
                         ObjEnterData, flag,                                        &  
                         keyword      = 'RELATIVE_HUMIDITY_FILE',                   &
                         ClientModule = 'RiverNetwork',                             &
                         SearchType   = FromFile,                                   &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR17'

            if (flag == 0) then
                write (*,*) 'Missing RELATIVE_HUMIDITY_FILE in Model input file'
                stop 'RiverNetwork - ConstructRiverNetwork - ERR18'
            endif

            call GetData(ColumnRelativeHumidity,                                    &   
                         ObjEnterData, flag,                                        &  
                         keyword      = 'RELATIVE_HUMIDITY_COLUMN',                 &
                         ClientModule = 'RiverNetwork',                             &
                         SearchType   = FromFile,                                   &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR19'

            call StartTimeSerieInput(ObjTSRelativeHumidity, FileName, ObjTime,      &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR20'

            !WindSpeed
            call GetData(FileName,                                                  &   
                         ObjEnterData, flag,                                        &  
                         keyword      = 'WIND_SPEED_FILE',                          &
                         ClientModule = 'RiverNetwork',                             &
                         SearchType   = FromFile,                                   &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR21'

            if (flag == 0) then
                write (*,*) 'Missing WIND_SPEED_FILE in Model input file'
                stop 'RiverNetwork - ConstructRiverNetwork - ERR22'
            endif

            call GetData(ColumnWindVelocity,                                        &   
                         ObjEnterData, flag,                                        &  
                         keyword      = 'WIND_SPEED_COLUMN',                        &
                         ClientModule = 'RiverNetwork',                             &
                         SearchType   = FromFile,                                   &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR23'


            call StartTimeSerieInput(ObjTSWindVelocity, FileName, ObjTime,      &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR24'

        endif


        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - RiverNetwork - ERR25'


    end subroutine ConstructRiverNetwork

    !--------------------------------------------------------------------------
    
    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        character(PathLength)                       :: DataFile
        integer                                     :: STAT_CALL
        character(PathLength)                       :: WatchFile, DTLogFile


        !Monitor Performance of the model execution?
        call ReadFileName('OUTWATCH', WatchFile, Message = 'Start Watch File', STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            STAT_CALL = CreateWatchGroup (WatchFile)
            MonitorPerformance  = .true.
        else
            MonitorPerformance  = .false.
        endif

        call ReadFileName('DT_LOG', DTLogFile, Message = 'Start DTLog File', STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then            
            MonitorDT  = .true.
        else
            MonitorDT  = .false.
        endif

        if (MonitorDT) then
            call UnitsManager (UnitDT, OPEN_FILE)      
            open(UNIT   = UnitDT, FILE   = DTLogFile, STATUS  = "UNKNOWN", IOSTAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - RiverNetwork - ERR01'
        end if


        call ReadFileName('IN_MODEL', DataFile, ModelName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - RiverNetwork - ERR02'

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - RiverNetwork - ERR03'

        call ReadTimeKeyWords   (ObjEnterData, FromFile, BeginTime, EndTime, DT,         &
                                 VariableDT, ModelName, MaxDT)

        
    end subroutine ReadKeywords
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine ModifyRiverNetwork
        !Local-----------------------------------------------------------------
        logical                                     :: Running
        integer                                     :: STAT_CALL
        real                                        :: CPUTime, LastCPUTime = 0.
        real                                        :: NewDT, AuxDT
        logical                                     :: NeedsRadiation
        logical                                     :: NeedAtmosphere
        real                                        :: SolarRadiation
        real                                        :: AirTemperature
        real                                        :: CloudCover
        real                                        :: RelativeHumidity
        real                                        :: WindVelocity

        write(*, *)"-------------------------- MOHID -------------------------"
        write(*, *)
        write(*, *)"Running MOHID River Network, please wait..."
        write(*, *)                    


        Running      = .true.
        CurrentTime  = BeginTime

        do while (Running)
            
            !Actualize the CurrentTime with Model time interval DT
            call ActualizeCurrentTime (TimeID    = ObjTime,                     &
                                       DT_Global = DT,                          &
                                       STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ModifyRiverNetwork - ERR01'
           
            !Gives the actualized Current time
            call GetComputeCurrentTime(ObjTime, CurrentTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ModifyRiverNetwork - ERR02'
            
            if (CurrentTime.LT.EndTime) then
                Running = .true.
            else
                Running = .false.
            endif

            !Check if Radiation is needed
            call GetNeedsRadiation (ObjDrainageNetwork, NeedsRadiation, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR03'

            !Check if Radiation is needed
            call GetNeedsAtmosphere(ObjDrainageNetwork, NeedAtmosphere, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ConstructRiverNetwork - ERR04'
        
            if (NeedAtmosphere .or. NeedsRadiation) then

                call GetOneTimeSeriesValue (ObjTSTopRadiation, ColumnTopRadiation, SolarRadiation)            

                call SetAtmosphereRiverNet (ObjDrainageNetwork, TopRadiation = SolarRadiation, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)  stop 'RiverNetwork - ConstructRiverNetwork - ERR06'

            endif


            if (NeedAtmosphere) then

                call GetOneTimeSeriesValue (ObjTSAirTemperature,    ColumnAirTemperature,   AirTemperature)            
                call GetOneTimeSeriesValue (ObjTSCloudCover,        ColumnCloudCover,       CloudCover)            
                call GetOneTimeSeriesValue (ObjTSRelativeHumidity,  ColumnRelativeHumidity, RelativeHumidity)            
                call GetOneTimeSeriesValue (ObjTSWindVelocity,      ColumnWindVelocity,     WindVelocity)            

                call SetAtmosphereRiverNet  (ObjDrainageNetwork,                                &
                                             AirTemperature      = AirTemperature,              &
                                             CloudCover          = CloudCover,                  &
                                             RelativeHumidity    = RelativeHumidity,            &
                                             WindSpeed           = WindVelocity,                &
                                             STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CalcEvapoTranspiration - ModuleBasin - ERR016'   

            endif



            call ModifyDrainageNetwork (ObjDrainageNetwork, STAT = STAT_CALL)                
            if (STAT_CALL .ne. SUCCESS_) stop 'RiverNetwork - ModifyRiverNetwork -  ERR03'
                
            if (VariableDT) then

                NewDT = min(DT * 1.1, MaxDT)
                call GetNextDrainageNetDT (ObjDrainageNetwork, AuxDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyBasin - ModuleBasin - ERR03'
                DT = min(NewDT, AuxDT)

                call ActualizeDT(TimeID = ObjTime, DT = DT, STAT = STAT_CALL)     
                if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ModifyRiverNetwork - ERR04'

            endif 

            if (MonitorDT) call WriteDTLog ('RiverNetwork', -99, NewDT)
       
            call CPU_TIME(CPUTime)
            if (CPUTime - LastCPUTime > 30.) then
                LastCPUTime = CPUTime
                call PrintProgress(ObjTime, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - ModifyRiverNetwork - ERR04'
            endif

        enddo
   
    end subroutine ModifyRiverNetwork

    !--------------------------------------------------------------------------


    subroutine GetOneTimeSeriesValue (ObjTimeSerie, DataColumn, ReturnValue)            

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTimeSerie
        integer                                     :: DataColumn
        real                                        :: ReturnValue

        !Local-----------------------------------------------------------------
        type (T_Time)                               :: Time1, Time2
        real                                        :: Value1, Value2
        logical                                     :: TimeCycle
        integer                                     :: STAT_CALL

        !Gets Value for current Time
        call GetTimeSerieValue (ObjTimeSerie, CurrentTime,                              &
                                DataColumn,                                             &
                                Time1, Value1, Time2, Value2, TimeCycle,                &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'RiverNetwork - ConstructRiverNetwork - ERR05'
        
        if (TimeCycle) then
            ReturnValue = Value1
        else
            !Interpolates Value for current instant
            call InterpolateValueInTime(CurrentTime,                                    &
                                        Time1, Value1,                                  &
                                        Time2, Value2, ReturnValue)
        endif      
        
    end subroutine  GetOneTimeSeriesValue             
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine KillRiverNetwork

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        call Write_Errors_Messages

        !Last Progress message
        call PrintProgress(ObjTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillRiverNetwork - ModuleRiverNetwork - ERR01'

        call KillDrainageNetwork (ObjDrainageNetwork)
        
        if (MonitorPerformance) then
            call KillWatchGroup (STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillRiverNetwork - MohidRiverNetwork - ERR02'
        endif

        if (MonitorDT) call UnitsManager (UnitDT, CLOSE_FILE)


        call StopCPUTime

        call ShutdownMohid ("Mohid RiverNetwork", ElapsedSeconds, TotalCPUTime)


    end subroutine KillRiverNetwork

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

    subroutine Write_Errors_Messages

        !Local-----------------------------------------------------------------
        real(8)                                     :: TotalInputVolume
        real(8)                                     :: TotalOutputVolume, TotalStoredVolume        
        character (Len = StringLength)              :: str_mass, string_to_be_written
        integer                                     :: STAT_CALL


        call GetVolumes (ObjDrainageNetwork,                                            &
                         TotalInputVolume  = TotalInputVolume,                          &
                         TotalOutputVolume = TotalOutputVolume,                         &
                         TotalStoredVolume = TotalStoredVolume,                         &
                         STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RiverNetwork - Write_Errors_Messages - ERR01'

       
        str_mass = ''
        
        write(str_mass, '(f20.8)') TotalInputVolume

        string_to_be_written = 'TotalInputVolume = '//trim(adjustl(adjustr(str_mass))) 

        call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)

        str_mass = ''
        
        write(str_mass, '(f20.8)') TotalOutputVolume

        string_to_be_written = 'TotalOutputVolume = '//trim(adjustl(adjustr(str_mass))) 

        call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)

        str_mass = ''
        
        write(str_mass, '(f20.8)') TotalStoredVolume

        string_to_be_written = 'TotalStoredVolume = '//trim(adjustl(adjustr(str_mass))) 

        call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)

    end subroutine Write_Errors_Messages

end program

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
