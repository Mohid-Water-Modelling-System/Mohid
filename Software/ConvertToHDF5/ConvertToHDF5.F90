!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : ConvertToHDF5
! PROGRAM       : ConvertToHDF5
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : July 2003
! REVISION      : Luis Fernandes
! DESCRIPTION   : ConvertToHDF5 to convert files in various formats to MOHID HDF5 format
!
!------------------------------------------------------------------------------

program ConvertToHDF5

    use ModuleGlobalData 
    use ModuleTime
    use ModuleStopWatch,      only : CreateWatchGroup, KillWatchGroup
    use ModuleEnterData
    use ModuleFunctions
    use ModuleHDF5
    use ModuleMM5Format
    use ModuleEUCenterFormat
#ifndef _NO_NETCDF
    use ModuleERA40Format
    use ModuleARPSFormat
    use ModuleMercatorFormat
    use ModuleARPSToWW3
    use ModuleHYCOMFormat
    use ModuleNetCDFCF_2_HDF5MOHID
#endif
    use ModuleLevitusFormat
    use ModuleHellermanRosensteinAscii
    use ModuleTecnoceanAscii
#ifndef _NO_NETCDF
    use ModuleCowamaAsciiWind
#endif
    use ModuleSWAN
    use ModuleReadSWANNonStationary
    use ModuleHDF5ToASCIIandBIN
    use ModuleGFSasciiWind
#ifndef _NO_NETCDF    
#ifdef _USE_PROJ4 
    use ModuleWRFFormat
    use ModuleCALMETFormat
#endif
#endif
#ifdef _USE_MODIS
!    use ModuleConvertModisL3Mapped
!    use ModuleConvertModisL3_V2
    use ModuleConvertModisL2
    use ModuleConvertOceanColorL2
#endif
#ifndef _NO_NETCDF
    use ModuleWOAFormat
#endif
    use ModuleInterpolateGrids
    use ModuleInterpolateTime
    use ModuleGlueHDF5Files
    use ModulePatchHDF5Files
    use ModuleUpscaleHDF5
#ifndef _NO_NETCDF
    use ModuleCFFormat
    use ModuleCFPolcomFormat
    use ModuleAladinFormat
    use ModuleMOG2DFormat
    use ModuleIHRadarFormat
#endif
    use ModuleDelft3D_2_Mohid

    implicit none

    type (T_Time)                           :: InitialSystemTime, FinalSystemTime
    real                                    :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)                   :: F95Time
    
    character(len = StringLength), parameter:: ConvertEUCenterFormatToHDF5  = 'CONVERT EU CENTER FORMAT'
    character(len = StringLength), parameter:: ConvertMM5FormatToHDF5       = 'CONVERT MM5 FORMAT'
    character(len = StringLength), parameter:: ConvertToARPSFormat          = 'CONVERT ARPS FORMAT'
    character(len = StringLength), parameter:: ConvertARPSToWW3Format       = 'CONVERT ARPS_TO_WW3 FORMAT'
    character(len = StringLength), parameter:: ConvertToCFFormat            = 'CONVERT CF FORMAT'
    character(len = StringLength), parameter:: ConvertToCFPolcomFormat      = 'CONVERT CF POLCOM FORMAT'
    character(len = StringLength), parameter:: ConvertERA40FormatToHDF5     = 'CONVERT ERA40 FORMAT'
    character(len = StringLength), parameter:: ConvertToMercatorFormat      = 'CONVERT MERCATOR FORMAT'    
    character(len = StringLength), parameter:: ConvertToNetCDFCF            = 'CONVERT NETCDF CF TO HDF5 MOHID'        
    character(len = StringLength), parameter:: ConvertToHYCOMFormat         = 'CONVERT HYCOM FORMAT'
    character(len = StringLength), parameter:: ConvertToLevitusFormat       = 'CONVERT LEVITUS FORMAT'
    character(len = StringLength), parameter:: ConvertToHellermanRosenstein = 'CONVERT HELLERMAN ROSENSTEIN ASCII'
    character(len = StringLength), parameter:: ConvertToTecnocean           = 'CONVERT TECNOCEAN ASCII'
    character(len = StringLength), parameter:: ConvertToCowama              = 'CONVERT COWAMA ASCII WIND'
    character(len = StringLength), parameter:: ConvertToAndFromSWAN         = 'CONVERT TO AND FROM SWAN'
    character(len = StringLength), parameter:: ConvertFromSWANNonStationary = 'CONVERT FROM SWAN NONSTATIONARY'
    character(len = StringLength), parameter:: ConvertHDF5ToSWANorMOHID     = 'CONVERT FROM HDF5 TO SWAN OR MOHID'
    character(len = StringLength), parameter:: ConvertToAladinFormat        = 'CONVERT ALADIN FORMAT'    
    character(len = StringLength), parameter:: ConvertMOG2DFormatToHDF5     = 'CONVERT MOG2D FORMAT'    
    character(len = StringLength), parameter:: ConvertGFS_from_ASCII2HDF    = 'CONVERT GFS FORMAT'  
    character(len = StringLength), parameter:: ConvertWRFFormatToHDF5       = 'CONVERT WRF FORMAT'
    character(len = StringLength), parameter:: ConvertCALMETFormatToHDF5    = 'CONVERT CALMET FORMAT'

#ifdef _USE_MODIS
!    character(len = StringLength), parameter:: ConvertModisL3Mapped         = 'CONVERT MODIS L3MAPPED'
    character(len = StringLength), parameter:: ConvertModisL2               = 'CONVERT MODIS L2'
    character(len = StringLength), parameter:: ConvertOCL2                  = 'CONVERT OCEAN COLOR L2'
!    character(len = StringLength), parameter:: ConvertModisL3               = 'CONVERT MODIS L3'
#endif
    character(len = StringLength), parameter:: ConvertToWOAFormat           = 'CONVERT WOA FORMAT'
    character(len = StringLength), parameter:: InterpolateGrids             = 'INTERPOLATE GRIDS'
    character(len = StringLength), parameter:: InterpolateTime              = 'INTERPOLATE TIME'    
    character(len = StringLength), parameter:: GluesHD5Files                = 'GLUES HDF5 FILES'
    character(len = StringLength), parameter:: PatchHD5Files                = 'PATCH HDF5 FILES'
    character(len = StringLength), parameter:: UpscaleHDF5                  = 'UPSCALE HDF5 FILES'
    character(len = StringLength), parameter:: ConvertIHRadarFormatToHDF5   = 'CONVERT IH RADAR FORMAT'
    
    character(len = StringLength), parameter:: ConvertDelft3DFormatToHDF5   = 'CONVERT DELFT3D FORMAT'    

    logical               :: WatchPassedAsArgument = .false.
    logical               :: Watch     = .false.
    character(PathLength) :: WatchFile
    character(PathLength) :: DataFile  = null_str
    
    call ReadArguments     
    
    call StartConvertToHDF5; call KillConvertToHDF5


    contains
    
    !--------------------------------------------------------------------------

    subroutine StartConvertToHDF5

        call StartUpMohid("ConvertToHDF5")

        call StartCPUTime

        call ReadOptions

    end subroutine StartConvertToHDF5
    
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
                
            case ('-h', '--help')
               
                call print_help()
                stop
                
            case ('-v', '--version')
            
                call print_version()
                stop                
                
            case ('-w', '--watch')
               
                if (last_arg > 0) then             
                
                     print *, 'Invalid parameter.'
                
                     if (i > n_args) then
                        call print_help()
                        stop                   
                     endif
                endif
                
                last_arg = 1
                
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
                   
                case (1)
                   
                    WatchPassedAsArgument = .true.
                    call get_command_argument(i, WatchFile)                
                    Watch = .true.  
                    
                case (2)
                   
                    call get_command_argument(i, DataFile)
                    
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
        print '(a)', 'ConvertToHDF5 usage: ConvertToHDF5 [OPTIONS]'
        print '(a)', ''
        print '(a)', 'Without further options, ConvertToHDF5 uses "ConvertToHDF5Action.dat" as input file.'
        print '(a)', ''
        print '(a)', 'ConvertToHDF5 options:'
        print '(a)', ''
        print '(a)', '  [-v, --version]     print version information and exit'
        print '(a)', '  [-h, --help]        print usage information and exit'
        print '(a)', '  [-w, --watch] file  Monitor performance and save data to "file"'
        print '(a)', '  [-c, --config] file Uses "file" as input configuration'
        print '(a)', ''
    end subroutine print_help    
    
    subroutine print_version ()
        print '(a)', ''
#if defined (CODEPLEXVERSION)
        print '(a, i0)', 'ConvertToHDF5 version (codeplex): ', CODEPLEXVERSION
#else
        print '(a)', 'ConvertToHDF5 version : PERSONAL'
#endif
        print '(a)', ''
    end subroutine print_version    
    
    !--------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        !character(PathLength)                       :: DataFile     = 'ConvertToHDF5Action.dat'
        integer                                     :: STAT_CALL
        integer                                     :: ObjEnterData = 0
        integer                                     :: ClientNumber, iflag
        logical                                     :: BlockFound
        character(len = StringLength)               :: Action
        character(PathLength)                       :: WatchFile

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Running ConvertToHDF5...'
        write(*,*)
        
        
        if (DataFile == null_str) then
        
            !Read input file name from nomfich file
            call ReadFileName('IN_MODEL', DataFile, "ConvertToHDF5", STAT = STAT_CALL)
            
            if     (STAT_CALL == FILE_NOT_FOUND_ERR_) then
                DataFile = 'ConvertToHDF5Action.dat'    
            elseif (STAT_CALL /= SUCCESS_              ) then
                stop 'ReadOptions - ConvertToHDF5 - ERR10'
            endif                        
                           
        endif


        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR20'

        !call GetData(WatchFile,                                             &
        !             ObjEnterData, iflag,                                   &
        !             SearchType   = FromFile,                               &
        !             keyword      = 'OUTWATCH',                             &
        !             ClientModule = 'ConvertToHDF5',                        &
        !             STAT         = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR02'
        !if (iflag == 0)then
        !    MonitorPerformance  = .false.
        !else
        !    STAT_CALL = CreateWatchGroup (WatchFile)
        !    MonitorPerformance  = .true.
        !endif
        
        if (.not. WatchPassedAsArgument) then
            call GetData(WatchFile,                                             &
                         ObjEnterData, iflag,                                   &
                         SearchType   = FromFile,                               &
                         keyword      = 'OUTWATCH',                             &
                         ClientModule = 'ConvertToHDF5',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR30'
            if (iflag == 0)then
                MonitorPerformance  = .false.
            else
                STAT_CALL = CreateWatchGroup (WatchFile)
                MonitorPerformance  = .true.
            endif
        else
            if (Watch) then
                STAT_CALL = CreateWatchGroup (WatchFile)
                MonitorPerformance  = .true.           
            endif
        endif      

do1 :   do
            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                 &
                                        '<begin_file>', '<end_file>', BlockFound,   &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR40'

if1 :       if(STAT_CALL .EQ. SUCCESS_) then    
if2 :           if (BlockFound) then


                    call GetData(Action,                                            &
                                 ObjEnterData, iflag,                               &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'ACTION',                           &
                                 ClientModule = 'ConvertToHDF5',                    &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR50'
                    if (iflag == 0)then
                        write(*,*)'Must specify type of file to convert'
                        stop 'ReadOptions - ConvertToHDF5 - ERR60'
                    end if

                    write(*,*) Action


                    select case(Action)


#ifndef _NO_NETCDF
                        case(ConvertEUCenterFormatToHDF5)

                            call ConvertEUCenterFormat(ObjEnterData, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR50'
#endif
                        
                        case(ConvertMM5FormatToHDF5)
                            
                            call ConvertMM5Format(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR60'

#ifndef _NO_NETCDF
                        case(ConvertToARPSFormat)

                            call ConvertARPSFormat(ObjEnterData, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR70'
#endif
#ifndef _NO_NETCDF
                        case(ConvertARPSToWW3Format)

                            call ConvertARPSWW3Format(ObjEnterData, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR80'
#endif
#ifndef _NO_NETCDF
                        case(ConvertERA40FormatToHDF5)
                            
                            call ConvertERA40Format(ObjEnterData, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR90'
#endif
#ifndef _NO_NETCDF
                        case(ConvertToMercatorFormat)
                            
                            call ConvertMercatorFormat(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR100'
                            
                        case(ConvertToNetCDFCF)
                            
                            call ConvertNetCDFCF_2_HDF5MOHID(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR105'
                            
#endif
#ifndef _NO_NETCDF
                        case(ConvertToHYCOMFormat)
                            
                            call ConvertHYCOMFormat(ObjEnterData, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR110'
#endif

                        case(ConvertToLevitusFormat)
                            
                            call ConvertLevitusFormat(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR120'


                        case(ConvertToHellermanRosenstein)
                            
                            call ConvertHellermanRosensteinAscii(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR130'


                        case(ConvertToTecnocean)
                            
                            call ConvertTecnoceanAscii(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR135'
#ifndef _NO_NETCDF
                        case(ConvertToCowama)
                            
                            call ConvertCowamaAsciiWind(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR138'
#endif

                        case(ConvertToAndFromSWAN)
                            
                            call ConvertSWAN(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR139'

                        case(InterpolateGrids)

                            call StartInterpolateGrids(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR140'

                        case(ConvertHDF5ToSWANorMOHID)

                            call ConvertHDF5ToASCIIandBIN(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR145'

                        case(GluesHD5Files)

                            call StartGlueHDF5Files(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR150'
                            
                        case(ConvertFromSWANNonStationary)
                            
                            call ConvertSWANNonStationary(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR160'


#ifdef _USE_MODIS
!                        case (ConvertModisL3Mapped)
                           
!                            call StartConvertModisL3Mapped(ObjEnterData, ClientNumber,  STAT = STAT_CALL)
!                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR160'                           


                        case (ConvertModisL2)
                           
                            call StartConvertModisL2(ObjEnterData, ClientNumber,  STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR170'
                                                  

                        case (ConvertOCL2)

                            call StartConvertOceanColorL2(ObjEnterData, ClientNumber,  STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR180'


!                        case (ConvertModisL3)

!                            call StartConvertModisL3(ObjEnterData, ClientNumber,  STAT = STAT_CALL)
!                           if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR190'
#endif
#ifndef _NO_NETCDF
                        case (ConvertToWOAFormat)

                          call ConvertWOAFormat(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                          if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR195'
#endif

                        case(PatchHD5Files)

                            call StartPatchHDF5Files(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR200'
#ifndef _NO_NETCDF
                        case(ConvertToCFFormat)

                            call ConvertCFFormat(ObjEnterData, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR205'
#endif
#ifndef _NO_NETCDF
                        case(ConvertToCFPolcomFormat)

                            call ConvertCFPolcomFormat(ObjEnterData, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR206'
#endif
#ifndef _NO_NETCDF
                        case(ConvertToAladinFormat)
                            
                            call ConvertAladinFormat(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR207'
#endif
#ifndef _NO_NETCDF
                        case(ConvertMOG2DFormatToHDF5)
                            
                            call ConvertMOG2DFormat(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR208'
#endif
                        case(ConvertGFS_from_ASCII2HDF)
                            
                            call ConvertGFSasciiWind(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR209'
#ifndef _NO_NETCDF
#ifdef _USE_PROJ4
                        case(ConvertWRFFormatToHDF5)
                            
                            call ConvertWRFFormat(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR210'

                        case(ConvertCALMETFormatToHDF5)
                            
                            call ConvertCALMETFormat(ObjEnterData, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR211'
#endif
#endif
                        case(InterpolateTime)

                            call StartInterpolateTime(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR260'

#ifndef _NO_NETCDF
                        case(ConvertIHRadarFormatToHDF5)
                        
                            call ConvertIHRadarFormat(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR270'
#endif                  

                        case (ConvertDelft3DFormatToHDF5)
                        
                            call ConvertDelft3D_2_Mohid(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            
                        case (UpscaleHDF5)
                        
                            call StartUpscaleHDF5(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR270'
                            
                        case default
                            
                            stop 'Option not known - ReadOptions - ConvertToHDF5 - ERR299'

                    end select


                else
                    call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR20'
                        
                    exit do1    !No more blocks

                end if if2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                if(STAT_CALL .ne. SUCCESS_)stop 'ReadOptions - ConvertToHDF5 - ERR230'
                    
            end if if1
        end do do1


        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR240'

    end subroutine ReadOptions
    
    !--------------------------------------------------------------------------

    subroutine KillConvertToHDF5

        integer                             :: STAT_CALL

        if (MonitorPerformance) then
            call KillWatchGroup (STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillConvertToHDF5 - ConvertToHDF5 - ERR01'
        endif

        call StopCPUTime

        call ShutdownMohid ("ConvertToHDF5", ElapsedSeconds, TotalCPUTime)

    end subroutine KillConvertToHDF5
    
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
    
!Lana
end program ConvertToHDF5
