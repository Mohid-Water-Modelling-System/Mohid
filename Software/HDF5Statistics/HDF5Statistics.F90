!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : HDF5Statistics
! PROGRAM       : HDF5Statistics
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2005
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Creates Statistics HDF5 from HDF5 files
!
!------------------------------------------------------------------------------

!DataFile
!
!   IN_MODEL                : char                  [-]         !Name of input file with
!                                                               !user's instructions                                                                
!   ROOT_SRT                : char                  [-]         !Path of folder where the
!                                                               !input file and HDF5 files
!                                                               !are and where output files
!                                                               !will appear                   
!  (file's name must be 'nomfich.dat')

program HDF5Statistics

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModuleStatistic
    use ModuleHDF5Statistics

    implicit none

    type (T_Time)               :: InitialSystemTime, FinalSystemTime
    real                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)       :: F95Time

    integer                     :: ObjHDF5StatisticsID = 0

    call ConstructHDF5Stats
    call ModifyHDF5Stats
    call KillHDF5Stats

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Stats
        
        call StartUpMohid("HDF5Statistics")

        call StartCPUTime

        call ReadKeywords

    end subroutine ConstructHDF5Stats
    
    !--------------------------------------------------------------------------

    subroutine ModifyHDF5Stats
        
        !Local-----------------------------------------------------------------

  
        call ModifyHDF5Statistics      
     
    end subroutine ModifyHDF5Stats
    
    !--------------------------------------------------------------------------

    subroutine KillHDF5Stats

        call KillHDF5Statistics

        call StopCPUTime

        call ShutdownMohid ("HDF5Statistics", ElapsedSeconds, TotalCPUTime)

    end subroutine KillHDF5Stats
    
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

        call ReadFileName('IN_MODEL', DataFile, "HDF5Statistics", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HDF5Statistics - ERR01'

        call StartHDF5StatsCreator(ObjHDF5StatisticsID, DataFile)      
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HDF5Statistics - ERR02'

    end subroutine ReadKeywords

end program HDF5Statistics
