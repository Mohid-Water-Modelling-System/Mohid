!------------------------------------------------------------------------------
!        Hidromod & IST, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : HDF 2 Surfer
! PROJECT       : Hdf2Surfer
! PROGRAM       : MainHdf2Surfer
! URL           : http://www.mohid.com
! AFFILIATION   : Hidromod
! DATE          : June 2012
! REVISION      : Paulo Leitão
! DESCRIPTION   : Hdf2Surfer to create main program to use MOHID modules
!
!------------------------------------------------------------------------------

program MohidHdf2Surfer

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
    use ModuleHdf2Surfer

    implicit none

    type (T_Time)               :: InitialSystemTime, FinalSystemTime
    real                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)       :: F95Time

    integer                     :: STAT_CALL, ObjHdf2SurferID


    call ConstructMain
    call ModifyMain
    call KillMain

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMain
        
        call StartUpMohid("Hdf2Surfer")

        call StartCPUTime

        ObjHdf2SurferID = 0
        
        call ConstructHdf2Surfer(ObjHdf2SurferID, InitialSystemTime, STAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) then
            stop "ConstructMain - MainHdf2Surfer - ERR10"
        endif

    end subroutine ConstructMain
    
    !--------------------------------------------------------------------------

    subroutine ModifyMain
        
        !Local-----------------------------------------------------------------

        call ModifyHdf2Surfer(ObjHdf2SurferID, STAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) then
            stop "ModifyMain - MainHdf2Surfer - ERR10"
        endif
    
    end subroutine ModifyMain
    
    !--------------------------------------------------------------------------

    subroutine KillMain
    
        call KillHdf2Surfer(ObjHdf2SurferID, STAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) then
            stop "KillMain - MainHdf2Surfer - ERR10"
        endif
    
        call StopCPUTime

        call ShutdownMohid ("Hdf2Surfer", ElapsedSeconds, TotalCPUTime)

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
 

end program MohidHdf2Surfer
