!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : PercentileComputation
! PROGRAM       : MainPercentileComputation
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig /Luis Fernandes - v4.0
! DESCRIPTION   : PercentileComputation to create main program to use MOHID modules
!
!------------------------------------------------------------------------------

program MohidPercentileComputation

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModulePercentileComputation

    implicit none

    type (T_Time)               :: InitialSystemTime, FinalSystemTime
    real                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)       :: F95Time

    integer                     :: STAT_CALL
    
    integer                     :: ObjPercentileComputationID  = 0


    call ConstructMohidPercentileComputation
    call ModifyMohidPercentileComputation
    call KillMohidPercentileComputation

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidPercentileComputation
        
        call StartUpMohid("MohidPercentileComputation")

        call StartCPUTime
        
        call ConstructPercentileComputation(ObjPercentileComputationID, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructMohidPercentileComputation - MohidPercentileComputation - ERR10'
        endif

    end subroutine ConstructMohidPercentileComputation
    
    !--------------------------------------------------------------------------

    subroutine ModifyMohidPercentileComputation
        
        !Local-----------------------------------------------------------------
        logical                                     :: Running


        call ModifyPercentileComputation(ObjPercentileComputationID, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructMohidPercentileComputation - MohidPercentileComputation - ERR10'
        endif
    
    end subroutine ModifyMohidPercentileComputation
    
    !--------------------------------------------------------------------------

    subroutine KillMohidPercentileComputation
    
        call KillPercentileComputation(ObjPercentileComputationID, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ConstructMohidPercentileComputation - MohidPercentileComputation - ERR10'
        endif    

        call StopCPUTime

        call ShutdownMohid ("MohidPercentileComputation", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMohidPercentileComputation
    
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
    
end program MohidPercentileComputation
