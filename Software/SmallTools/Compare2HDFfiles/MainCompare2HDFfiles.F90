!------------------------------------------------------------------------------
!        Hidromod & IST, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : HDF 2 Surfer
! PROJECT       : Compare2HDFfiles
! PROGRAM       : MainCompare2HDFfiles
! URL           : http://www.mohid.com
! AFFILIATION   : Hidromod
! DATE          : June 2012
! REVISION      : Paulo Leitão
! DESCRIPTION   : Compare2HDFfiles to create main program to use MOHID modules
!
!------------------------------------------------------------------------------

program Compare2HDFfiles

    use ModuleGlobalData
    use ModuleTime
    use ModuleStopWatch
    use ModuleEnterData
    use ModuleFunctions
    use ModuleDrawing     
    use ModuleHDF5
    use ModuleHorizontalGrid
    use ModuleGridData
    use ModuleHorizontalMap 
    use ModuleGeometry     
    use ModuleMap 
    use ModuleStatistic
    use ModuleTimeSerie
    use ModuleCompare2HDFfiles

    implicit none

    type (T_Time)               :: InitialSystemTime, FinalSystemTime
    real                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)       :: F95Time

    integer                     :: STAT_CALL, ObjCompare2HDFfilesID


    call ConstructMain
    call ModifyMain
    call KillMain

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMain
        
        call StartUpMohid("Compare2HDFfiles")

        call StartCPUTime

        ObjCompare2HDFfilesID = 0
        
        call ConstructCompare2HDFfiles(ObjCompare2HDFfilesID, InitialSystemTime, STAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) then
            stop "ConstructMain - MainCompare2HDFfiles - ERR10"
        endif

    end subroutine ConstructMain
    
    !--------------------------------------------------------------------------

    subroutine ModifyMain
        
        !Local-----------------------------------------------------------------

        call ModifyCompare2HDFfiles(STAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) then
            stop "ModifyMain - MainCompare2HDFfiles - ERR10"
        endif
    
    end subroutine ModifyMain
    
    !--------------------------------------------------------------------------

    subroutine KillMain
    
        call KillCompare2HDFfiles(STAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) then
            stop "KillMain - MainCompare2HDFfiles - ERR10"
        endif
    
        call StopCPUTime

        call ShutdownMohid ("Compare2HDFfiles", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMain
    
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
 

end program Compare2HDFfiles
