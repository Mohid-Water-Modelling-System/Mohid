!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : CreateVerticalResolution
! PROGRAM       : CreateVerticalResolution
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : January 2006
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Creates vertical resolution (geometry) based on a property 
!                 vertical profile
!
!------------------------------------------------------------------------------

!DataFile
!
!   IN_MODEL                : char                  [-]         !Name of input file with
!                                                               !user's instructions            
!   ROOT_SRT                : char                  [-]         !Path of folder where the
!                                                               !input files are and where 
!                                                               !output file will appear
!  (file's name must be 'nomfich.dat')

program CreateVerticalResolution

    use ModuleGlobalData,           only: PathLength, SUCCESS_, StartUpMohid, ShutdownMohid 
    use ModuleTime,                 only: T_Time, SetDate, operator(-), operator(+) 
    use ModuleEnterData,            only: ReadFileName
    use ModuleVerticalResolution

    implicit none

    type (T_Time)               :: InitialSystemTime, FinalSystemTime
    real                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)       :: F95Time

    integer                     :: ObjVerticalResolutionID = 0

    call ConstructCreateVertResolution
    call ModifyCreateVertResolution
    call KillCreateVertResolution

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructCreateVertResolution
        
        call StartUpMohid("CreateVerticalResolution")

        call StartCPUTime

        call ReadKeywords

    end subroutine ConstructCreateVertResolution
    
    !--------------------------------------------------------------------------

    subroutine ModifyCreateVertResolution
        
        !Local-----------------------------------------------------------------

        !------------------------------------------------------------------------

        call ModifyVerticalResolution   
    
    end subroutine ModifyCreateVertResolution
    
    !--------------------------------------------------------------------------

    subroutine KillCreateVertResolution

        call KillVerticalResolution

        call StopCPUTime

        write(*,*)

        call ShutdownMohid ("CreateVerticalResolution", ElapsedSeconds, TotalCPUTime)

    end subroutine KillCreateVertResolution
    
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

        !------------------------------------------------------------------------

        call ReadFileName('IN_MODEL', DataFile, "CreateVerticalResolution", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - CreateVerticalResolution - ERR10'

        call StartVerticalResolution(ObjVerticalResolutionID, DataFile)      
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - CreateVerticalResolution - ERR20'

    end subroutine ReadKeywords

end program CreateVerticalResolution
