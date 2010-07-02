!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : CEQUALW2 preprocessor
! PROGRAM       : MainReservoirOptimization
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Ricardo Miranda - v1.0
! DESCRIPTION   : 
!
!------------------------------------------------------------------------------

program MainReservoirOptimization

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModuleReservoirOptimization

    implicit none
    
    type(T_Time)                :: BeginTime, EndTime, CurrentTime
    real                        :: DT
    logical                     :: VariableDT
    type (T_Time)               :: InitialSystemTime, FinalSystemTime
    real                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)       :: F95Time

    integer                     :: ObjTime                  = 0
    integer                     :: ObjReservoirOptimization = 0
    integer                     :: STAT_CALL


    call ConstructMainReservoirOptimization
    call ModifyMainReservoirOptimization
    call KillMainReservoirOptimization

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMainReservoirOptimization
        
        call StartUpMohid("MainReservoirOptimization")

        call StartCPUTime

        call ReadKeywords

        call StartComputeTime(ObjTime, BeginTime, EndTime, DT = DT, VariableDT = VariableDT, STAT = STAT_CALL)
        if (STAT_CALL .ne. SUCCESS_) stop 'ConstructMainReservoirOptimization - ERR01'

        call ConstructReservoirOptimization(ObjReservoirOptimization, STAT = STAT_CALL)
        if (STAT_CALL .ne. SUCCESS_) stop 'ConstructMainReservoirOptimization - ERR02'
        
    end subroutine ConstructMainReservoirOptimization
    
    !--------------------------------------------------------------------------

    subroutine ModifyMainReservoirOptimization
        
        !Local-----------------------------------------------------------------
        logical                                     :: Running


        Running      = .true.
        CurrentTime  = BeginTime

        do while (Running)
            
            CurrentTime = CurrentTime + DT

            call ActualizeCurrentTime(ObjTime, DT, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'ModifyMainReservoirOptimization - ERR01'


            if (abs(CurrentTime - EndTime) > DT / 10.) then
                Running = .true.
            else
                Running = .false.
            endif

            call PrintProgress(ObjTime, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'ModifyMainReservoirOptimization - ERR02'
        enddo
    
        call ModifyReservoirOptimization(ObjReservoirOptimization, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'ModifyMainReservoirOptimization - ERR03'
            
        call WriteW2Bathymetry(ObjReservoirOptimization, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'ModifyMainReservoirOptimization - ERR04'
        
    
    end subroutine     
    
    !--------------------------------------------------------------------------

    subroutine KillMainReservoirOptimization

        call StopCPUTime

        call ShutdownMohid ("MainReservoirOptimization", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMainReservoirOptimization
    
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
        character(PathLength)                       :: DataFile
        integer                                     :: STAT_CALL
        integer                                     :: ObjEnterData = 0
        integer                                     :: FromFile

        call ReadFileName('IN_MODEL', DataFile, "MainReservoirOptimization", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainReservoirOptimization - ERR01'

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainReservoirOptimization - ERR02'

        call GetExtractType     (FromFile = FromFile)

        call ReadTimeKeyWords   (ObjEnterData, FromFile, BeginTime, EndTime, DT,         &
                                 VariableDT, "MainReservoirOptimization")

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainReservoirOptimization - ERR03'

    end subroutine ReadKeywords

end program MainReservoirOptimization

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

