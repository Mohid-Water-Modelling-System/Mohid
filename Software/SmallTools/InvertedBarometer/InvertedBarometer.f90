!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : InvertedBarometer
! PROGRAM       : InvertedBarometer
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : October 2005
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Creates mean water level time serie for gauge locations
!                 based on the inverted barometer hypothesis
!
!------------------------------------------------------------------------------

!DataFile
!
!   IN_MODEL                : char                  [-]         !Name of input file with
!                                                               !user's instructions                                                                
!   ROOT_SRT                : char                  [-]         !Path of folder where the
!                                                               !input files are and where 
!                                                               !output files will appear
!  (file's name must be 'nomfich.dat')


program InvertedBarometer

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData,            only : ReadFileName
    use ModuleInvertedBarometer,    only : StartInvertedBarometer,          & 
                                           ModifyInvertedBarometerLevel,    &
                                           KillInvertedBarometerLevel

    implicit none

    type (T_Time)                   :: InitialSystemTime, FinalSystemTime
    real                            :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)           :: F95Time

    integer                         :: ObjInvertedBarometerID = 0


    call ConstructInvertedBarometer
    call ModifyInvertedBarometer
    call KillInvertedBarometer

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructInvertedBarometer
        
        call StartUpMohid("InvertedBarometer")

        call StartCPUTime

        call ReadKeywords

    end subroutine ConstructInvertedBarometer
    
    !--------------------------------------------------------------------------

    subroutine ModifyInvertedBarometer
        
        !Local-----------------------------------------------------------------

        call ModifyInvertedBarometerLevel
    
    end subroutine ModifyInvertedBarometer
    
    !--------------------------------------------------------------------------

    subroutine KillInvertedBarometer

        call KillInvertedBarometerLevel

        call StopCPUTime

        write(*,*)

        call ShutdownMohid ("InvertedBarometer", ElapsedSeconds, TotalCPUTime)

    end subroutine KillInvertedBarometer
    
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

        call ReadFileName('IN_MODEL', DataFile, "InvertedBarometer", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - InvertedBarometer - ERR10'

        call StartInvertedBarometer(ObjInvertedBarometerID, DataFile)      
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - InvertedBarometer - ERR20'

    end subroutine ReadKeywords

end program InvertedBarometer
