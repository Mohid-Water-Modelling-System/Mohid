!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Shell
! PROGRAM       : MainShell
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig /Luis Fernandes - v4.0
! DESCRIPTION   : Shell to create main program to use MOHID modules
!
!------------------------------------------------------------------------------

program MohidShell

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions

    implicit none

    type(T_Time)                :: BeginTime, EndTime, CurrentTime
    real                        :: DT
    logical                     :: VariableDT
    type (T_Time)               :: InitialSystemTime, FinalSystemTime
    real                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)       :: F95Time

    integer                     :: ObjTime = 0
    integer                     :: STAT_CALL


    call ConstructMohidShell
    call ModifyMohidShell
    call KillMohidShell

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidShell
        
        call StartUpMohid("MohidShell")

        call StartCPUTime

        call ReadKeywords

        call StartComputeTime(ObjTime, BeginTime, EndTime, DT = DT, VariableDT = VariableDT, STAT = STAT_CALL)

    end subroutine ConstructMohidShell
    
    !--------------------------------------------------------------------------

    subroutine ModifyMohidShell
        
        !Local-----------------------------------------------------------------
        logical                                     :: Running


        Running      = .true.
        CurrentTime  = BeginTime

        do while (Running)
            
            CurrentTime = CurrentTime + DT

            call ActualizeCurrentTime(ObjTime, DT, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'ModifyMohidShell - ERR01'


            if (abs(CurrentTime - EndTime) > DT / 10.) then
                Running = .true.
            else
                Running = .false.
            endif

            call PrintProgress(ObjTime, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'ModifyMohidShell - ERR02'

        enddo
    
    
    
    end subroutine ModifyMohidShell
    
    !--------------------------------------------------------------------------

    subroutine KillMohidShell

        call StopCPUTime

        call ShutdownMohid ("MohidShell", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMohidShell
    
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

        call ReadFileName('IN_MODEL', DataFile, "MohidShell", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidShell - ERR01'

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidShell - ERR02'

        call GetExtractType     (FromFile = FromFile)

        call ReadTimeKeyWords   (ObjEnterData, FromFile, BeginTime, EndTime, DT,         &
                                 VariableDT, "MohidShell")

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidShell - ERR03'

    end subroutine ReadKeywords

end program MohidShell
