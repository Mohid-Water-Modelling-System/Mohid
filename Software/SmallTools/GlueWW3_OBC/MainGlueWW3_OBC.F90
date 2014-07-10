!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : GlueWW3_OBC
! PROGRAM       : MainGlueWW3_OBC
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2014
! REVISION      : Paulo Leitão - v4.0
! DESCRIPTION   : GlueWW3_OBC to glue boundary WW3 files
!
!------------------------------------------------------------------------------

program MainGlueWW3_OBC

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModuleGlueWW3_OBC

    implicit none

    type(T_Time)                :: BeginTime, EndTime, CurrentTime
    real                        :: DT
    logical                     :: VariableDT
    type (T_Time)               :: InitialSystemTime, FinalSystemTime
    real                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)       :: F95Time

    integer                     :: ObjTime          = 0
    integer                     :: ObjGlueWW3_OBCID = 0
    integer                     :: STAT_CALL


    call ConstructMainGlueWW3_OBC
    call ModifyMainGlueWW3_OBC
    call KillMainGlueWW3_OBC

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMainGlueWW3_OBC
        
        call StartUpMohid("MainGlueWW3_OBC")

        call StartCPUTime

        call ReadKeywords

        call StartComputeTime(ObjTime, InitialSystemTime, BeginTime, EndTime, DT = DT,  &
                              VariableDT = VariableDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMainGlueWW3_OBC - MainGlueWW3_OBC - ERR10'
        
        call ConstructGlueWW3_OBC(ObjGlueWW3_OBCID, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMainGlueWW3_OBC - MainGlueWW3_OBC - ERR20'

    end subroutine ConstructMainGlueWW3_OBC
    
    !--------------------------------------------------------------------------

    subroutine ModifyMainGlueWW3_OBC
        
        !Local-----------------------------------------------------------------
        logical                                     :: Running


        call ModifyGlueWW3_OBC(ObjGlueWW3_OBCID, BeginTime, EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyMainGlueWW3_OBC - MainGlueWW3_OBC - ERR10'


    end subroutine ModifyMainGlueWW3_OBC
    
    !--------------------------------------------------------------------------

    subroutine KillMainGlueWW3_OBC
    
        call KillGlueWW3_OBC(ObjGlueWW3_OBCID, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillMainGlueWW3_OBC - MainGlueWW3_OBC - ERR10'

        call StopCPUTime

        call ShutdownMohid ("MainGlueWW3_OBC", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMainGlueWW3_OBC
    
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
        integer                                     :: FromFile, iflag

        !Begin-----------------------------------------------------------------

        call ConstructEnterData (ObjEnterData, FileName = 'GlueWW3_OBC.dat', STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainGlueWW3_OBC - ERR10'

        call GetExtractType     (FromFile = FromFile)
        
        call GetData(BeginTime,                                                         &
                     ObjEnterData,iflag,                                                &
                     SearchType     = FromFile,                                         &
                     keyword        = 'START',                                          &
                     ClientModule   = 'MainGlueWW3_OBC',                                &
                     STAT           = STAT_CALL)              
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainGlueWW3_OBC - ERR20'
        if (iflag     == 0       ) stop 'ReadKeywords - MainGlueWW3_OBC - ERR30' 
                
        call GetData(EndTime,                                                           &
                     ObjEnterData,iflag,                                                &
                     SearchType     = FromFile,                                         &
                     keyword        = 'END',                                            &
                     ClientModule   = 'MainGlueWW3_OBC',                                &
                     STAT           = STAT_CALL)              
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainGlueWW3_OBC - ERR40'
        if (iflag     == 0       ) stop 'ReadKeywords - MainGlueWW3_OBC - ERR50'


        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainGlueWW3_OBC - ERR60'

    end subroutine ReadKeywords

end program MainGlueWW3_OBC
