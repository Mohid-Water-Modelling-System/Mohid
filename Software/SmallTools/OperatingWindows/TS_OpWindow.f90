!  TS_OpWindow.f90 
!
!  FUNCTIONS:
!  TS_OpWindow - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: TS_OpWindow
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program TS_OpWindow

    use ModuleGlobalData            
    use ModuleTime    
    !use Module_TS_Synch
    use Module_TS_Operator
    
    !use ModuleTS_OpWindow    

    implicit none
    
    !Variables 
    real                            :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)           :: F95Time
    type (T_Time)                   :: FinalSystemTime, InitialSystemTime

    !integer                         :: ObjTS_OpWindow = 0
    !integer                         :: ObjTS_Synch    = 0    
    integer                         :: Obj_TS_Operator = 0
    
  
    !Begin---------------------------------------------------------

    call StartUpMohid("Time Series Operating Windows")
    
    call StartCPUTime
        
    !call ConstructTS_OpWindow(ObjTS_OpWindow)
    !
    !call ModifyTS_OpWindow   (ObjTS_OpWindow)
    !
    !call KillTS_OpWindow     (ObjTS_OpWindow)
    
    !call Construct_TS_Synch   (ObjTS_Synch   )  
    !call Modify_TS_Synch      (ObjTS_Synch   )
    !call Kill_TS_Synch        (ObjTS_Synch   )    
    
    call Construct_TS_Operator(Obj_TS_Operator)    
    call Modify_TS_Operator   (Obj_TS_Operator) 
    call Kill_TS_Operator     (Obj_TS_Operator) 
    
    call StopCPUTime
    
    call ShutdownMohid ("Time Series Operating Windows", ElapsedSeconds, TotalCPUTime)    
    
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
    
end program TS_OpWindow