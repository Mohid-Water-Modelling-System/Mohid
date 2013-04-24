!  NetworkStatistics.f90 
!
!  FUNCTIONS:
!  NetworkStatistics - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: NetworkStatistics
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program NetworkStatistics

    use ModuleGlobalData            
    use ModuleTime                  
    use ModuleNetworkStatistics    

    implicit none
    
    !Variables 
    real                            :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)           :: F95Time
    type (T_Time)                   :: FinalSystemTime, InitialSystemTime

    integer                         :: ObjNetworkStatistics = 0
  
    !Begin---------------------------------------------------------

    call StartUpMohid("Network Statistics")
    
    call StartCPUTime
        
    call ConstructNetworkStatistics(ObjNetworkStatistics)
    
    call ModifyNetworkStatistics   (ObjNetworkStatistics)
    
    call KillNetworkStatistics     (ObjNetworkStatistics)
    
    call StopCPUTime
    
    call ShutdownMohid ("Network Statistics", ElapsedSeconds, TotalCPUTime)    
    
    contains
    
    !--------------------------------------------------------------------------

    subroutine StartCPUTime

        call date_and_time(Values = F95Time)
        
        call SetDate      (InitialSystemTime, float(F95Time(1)), float(F95Time(2)),     &
                                              float(F95Time(3)), float(F95Time(5)),     &
                                              float(F95Time(6)), float(F95Time(7))+     &
                                              F95Time(8)/1000.)

    end subroutine StartCPUTime
    
    !--------------------------------------------------------------------------

    subroutine StopCPUTime

        call date_and_time(Values = F95Time)
        
        call SetDate      (FinalSystemTime,   float(F95Time(1)), float(F95Time(2)),     &
                                              float(F95Time(3)), float(F95Time(5)),     &
                                              float(F95Time(6)), float(F95Time(7))+     &
                                              F95Time(8)/1000.)
        
        call cpu_time(TotalCPUTime)

        ElapsedSeconds = FinalSystemTime - InitialSystemTime

    end subroutine StopCPUTime    
    
end program NetworkStatistics