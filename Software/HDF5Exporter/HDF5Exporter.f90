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

    character(PathLength)       :: DataFile  = 'HDF5Exporter.dat'
    logical                     :: ConfigByArgument = .false.
      
    call ReadArguments
    call ConstructHDF5Exporter
    call ModifyHDF5Exporter
    call KillHDF5Exporter

    contains
    
    !--------------------------------------------------------------------------    
      
    subroutine ReadArguments
    
        integer         :: i
        integer         :: n_args
        character(1024) :: arg
        integer         :: last_arg = 0
        
        n_args = command_argument_count()
        
        do i = 1, n_args
           
            call get_command_argument(i, arg)

            select case (arg)
            
            case ('-v', '--version')
            
                call print_version()
                stop
                
            case ('-h', '--help')
               
                call print_help()
                stop
                                               
            case ('-c', '--config')
                if (last_arg > 0) then             
                
                     print *, 'Invalid parameter.'
                
                     if (i > n_args) then
                        call print_help()
                        stop                   
                     endif
                endif
                
                last_arg = 2                
            
            case default
                select case (last_arg)
                    
                case (2)
                   
                    call get_command_argument(i, DataFile)
                    ConfigByArgument = .true.
                    
                case default
                   
                   print *, 'Invalid parameter: ', arg
                   stop
                   
                end select
                
                last_arg = 0
                
            end select
        end do    
    
    end subroutine ReadArguments
    
    subroutine print_help()
        print '(a)', ''
        print '(a)', 'HDF5Exporter usage: HDF5Exporter [OPTIONS]'
        print '(a)', ''
        print '(a)', 'ConvertToHDF5 options:'
        print '(a)', ''
        print '(a)', '  [-v, --version]     : print version information and exit'
        print '(a)', '  [-h, --help]        : Print usage information and exit'
        print '(a)', '  [-c, --config] file : Uses "file" as input configuration'
        print '(a)', ''        
        print '(a)', 'Without -c option, HDF5Exporter uses "nomfich.dat" to know the name of input file.'
        print '(a)', 'The name of the input file is provided using IN_MODEL keyword, inside "nomfich.dat"'
    end subroutine print_help    
    
    subroutine print_version ()
        print '(a)', ''
#if defined (CODEPLEXVERSION)
        print '(a, i0)', 'HDF5Exporter version (codeplex): ', CODEPLEXVERSION
#else
        print '(a)', 'HDF5Exporter version : PERSONAL'
#endif
        print '(a)', ''
    end subroutine print_version
    
    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Exporter
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL

        call StartUpMohid("HDF5Exporter")

        if (.not. ConfigByArgument) then
            call ReadFileName('IN_MODEL', DataFile, "HDF5Exporter", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HDF5Exporter - ERR10'        
        endif
        
        FilesName = DataFile
        
        call StartCPUTime

        call StartExportHDF5ToTimeSerie(ObjExportHDF5ToTimeSerieID, DataFile)      
        !if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HDF5Exporter - ERR20'

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
