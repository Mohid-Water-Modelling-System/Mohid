!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : HDF5Exporter files
! PROGRAM       : ExportToTimeSerie
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : April 2004
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Creates Time Series from HDF5 files
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


program HDF5Exporter

    use ModuleGlobalData              
    use ModuleTime                 
    use ModuleEnterData,                only : ReadFileName
    use ModuleExportHDF5ToTimeSerie,    only : StartExportHDF5ToTimeSerie,  &
                                               ModifyExportHDF5ToTimeSerie, &
                                               KillExportHDF5ToTimeSerie

    implicit none

    type (T_Time)               :: InitialSystemTime, FinalSystemTime
    real                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)       :: F95Time

    integer                     :: ObjExportHDF5ToTimeSerieID = 0

    call ConstructHDF5Exporter
    call ModifyHDF5Exporter
    call KillHDF5Exporter

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Exporter
        
        !Local-----------------------------------------------------------------
        character(PathLength)       :: DataFile
        integer                     :: STAT_CALL

        call StartUpMohid("HDF5Exporter")

        call ReadFileName('IN_MODEL', DataFile, "HDF5Exporter", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HDF5Exporter - ERR10'

        call StartCPUTime

        call StartExportHDF5ToTimeSerie(ObjExportHDF5ToTimeSerieID, DataFile)      
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HDF5Exporter - ERR20'

    end subroutine ConstructHDF5Exporter

    !--------------------------------------------------------------------------

    subroutine ModifyHDF5Exporter

        !Local-----------------------------------------------------------------

        call ModifyExportHDF5ToTimeSerie      
           
    end subroutine ModifyHDF5Exporter
    
    !--------------------------------------------------------------------------

    subroutine KillHDF5Exporter

        !Local-----------------------------------------------------------------

        call KillExportHDF5ToTimeSerie

        call StopCPUTime

        call ShutdownMohid ("HDF5Exporter", ElapsedSeconds, TotalCPUTime)

    end subroutine KillHDF5Exporter
    
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
    
end program HDF5Exporter
