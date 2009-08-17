!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : HDF5Extractor
! PROGRAM       : HDF5Extractor
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : October 2005
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Creates a HDF5 file from a time window of a supplied 
!                 HDF5 file.
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

program HDF5Extractor

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData,         only : ReadFileName
    use ModuleHDF5Extractor,     only : StartHDF5Extractor, ModifyExtractHDF5,  &
                                        KillExtractHDF5

    implicit none

    type (T_Time)               :: InitialSystemTime, FinalSystemTime
    real                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)       :: F95Time

    integer                     :: ObjHDF5ExtractorID = 0


    call ConstructHDF5Extractor
    call ModifyHDF5Extractor
    call KillHDF5Extractor

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Extractor
        
        call StartUpMohid("HDF5Extractor")

        call StartCPUTime

        call ReadKeywords

    end subroutine ConstructHDF5Extractor
    
    !--------------------------------------------------------------------------

    subroutine ModifyHDF5Extractor
        
        !Local-----------------------------------------------------------------

        call ModifyExtractHDF5     
    
    end subroutine ModifyHDF5Extractor
    
    !--------------------------------------------------------------------------

    subroutine KillHDF5Extractor

        call KillExtractHDF5

        call StopCPUTime

        write(*,*)

        call ShutdownMohid ("HDF5Extractor", ElapsedSeconds, TotalCPUTime)

    end subroutine KillHDF5Extractor
    
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

        call ReadFileName('IN_MODEL', DataFile, "HDF5Extractor", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HDF5Extractor - ERR01'

        call StartHDF5Extractor(ObjHDF5ExtractorID, DataFile)      
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HDF5Extractor - ERR02'

    end subroutine ReadKeywords

end program HDF5Extractor
