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
    use ModuleHDF5ToASCIIandBIN
    use ModuleGFSasciiWind
#ifndef _NO_NETCDF    
#ifdef _USE_PROJ4 
    use ModuleWRFFormat
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
#ifndef _NO_NETCDF
    use ModuleCFFormat
    use ModuleCFPolcomFormat
    use ModuleAladinFormat
    use ModuleMOG2DFormat
    use ModuleIHRadarFormat
#endif

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
    character(len = StringLength), parameter:: ConvertHDF5ToSWANorMOHID     = 'CONVERT FROM HDF5 TO SWAN OR MOHID'
    character(len = StringLength), parameter:: ConvertToAladinFormat        = 'CONVERT ALADIN FORMAT'    
    character(len = StringLength), parameter:: ConvertMOG2DFormatToHDF5     = 'CONVERT MOG2D FORMAT'    
    character(len = StringLength), parameter:: ConvertGFS_from_ASCII2HDF    = 'CONVERT GFS FORMAT'  
    character(len = StringLength), parameter:: ConvertWRFFormatToHDF5       = 'CONVERT WRF FORMAT'
   

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
    character(len = StringLength), parameter:: ConvertIHRadarFormatToHDF5   = 'CONVERT IH RADAR FORMAT'


    call StartConvertToHDF5; call KillConvertToHDF5


    contains
    
    !--------------------------------------------------------------------------

    subroutine StartConvertToHDF5

        call StartUpMohid("ConvertToHDF5")

        call StartCPUTime

        call ReadOptions

    end subroutine StartConvertToHDF5
    
    !--------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        character(PathLength)                       :: DataFile     = 'ConvertToHDF5Action.dat'
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

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR01'

        call GetData(WatchFile,                                             &
                     ObjEnterData, iflag,                                   &
                     SearchType   = FromFile,                               &
                     keyword      = 'OUTWATCH',                             &
                     ClientModule = 'ConvertToHDF5',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR02'
        if (iflag == 0)then
            MonitorPerformance  = .false.
        else
            STAT_CALL = CreateWatchGroup (WatchFile)
            MonitorPerformance  = .true.
        endif

do1 :   do
            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                 &
                                        '<begin_file>', '<end_file>', BlockFound,   &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR20'

if1 :       if(STAT_CALL .EQ. SUCCESS_) then    
if2 :           if (BlockFound) then


                    call GetData(Action,                                            &
                                 ObjEnterData, iflag,                               &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'ACTION',                           &
                                 ClientModule = 'ConvertToHDF5',                    &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR30'
                    if (iflag == 0)then
                        write(*,*)'Must specify type of file to convert'
                        stop 'ReadOptions - ConvertToHDF5 - ERR40'
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
#endif
#endif
                        case(InterpolateTime)

                            call StartInterpolateTime(ObjEnterData, ClientNumber, STAT = STAT_CALL)
                            if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - ConvertToHDF5 - ERR260'

                        case(ConvertIHRadarFormatToHDF5)
                        
                            call ConvertIHRadarFormat(ObjEnterData, ClientNumber, STAT = STAT_CALL)
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
