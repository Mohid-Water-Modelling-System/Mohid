!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Land
! PROGRAM       : Mohid Land
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module Basin is the top level of Basin and Aquifer 
!
!------------------------------------------------------------------------------
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License 
!version 2, as published by the Free Software Foundation.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!
!------------------------------------------------------------------------------

#ifdef _OPENMI_
module MohidLand
#else
program MohidLand
#endif

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime
    use ModuleFunctions,      only : ReadTimeKeyWords
    use ModuleBasin
    use ModuleStopWatch,      only : CreateWatchGroup, KillWatchGroup


    implicit none

    !Parameters

    !Instance IDs
    integer                             :: ObjComputeTime       = 0
    integer                             :: ObjBasin             = 0

    !Time Variables
    type (T_Time)                       :: BeginTime, EndTime, CurrentTime
    real                                :: DT                   = null_real
    real                                :: MaxDT                = null_real
    logical                             :: VariableDT           = .false.
    real                                :: GmtReference         = null_real
    real                                :: DTPredictionInterval = null_real

    !Model Name
    character(len=StringLength)         :: ModelName            = null_str
    character(PathLength)               :: OutputFile           = null_str
    logical                             :: SaveOutput           = .false.
    character(PathLength)               :: DataFile             = 'nomfich.dat'
    logical                             :: ConfigByArgument     = .false.    

    !Other Stuff
    type (T_Time)                       :: InitialSystemTime, FinalSystemTime
    integer, dimension(8)               :: F95Time              = null_int
    
    !OpenMI flag
    logical                             :: ModelConstructed     = .false.

#ifndef _OPENMI_

    call ReadArguments
    call ConstructMohidLand
    call ModifyMohidLand
    call KillMohidLand
#endif

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
        print '(a)', 'MohidLand usage: MohidLand [OPTIONS]'
        print '(a)', ''
        print '(a)', 'MohidLand options:'
        print '(a)', ''
        print '(a)', '  [-v, --version]     : print version information and exit'
        print '(a)', '  [-h, --help]        : Print usage information and exit'
        print '(a)', '  [-c, --config] file : Uses "file" as input configuration (nomfich.dat by default)'
        print '(a)', ''        
    end subroutine print_help    
    
    subroutine print_version ()
        print '(a)', ''
#if defined (CODEPLEXVERSION)
        print '(a, i0)', 'MohidLand version (codeplex): ', CODEPLEXVERSION
#else
        print '(a)', 'MohidLand version : PERSONAL'
#endif
        print '(a)', ''
    end subroutine print_version
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidLand

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        if (ConfigByArgument) then
            call SetInputFullPath (DataFile)
        endif

        call StartupMohid ("Mohid Land")

        !Gets the actual time
        call date_and_time(Values = F95Time)
        call SetDate      (InitialSystemTime, float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)
        
        !Reads the Keywords
        call ReadKeywords

        !Constructs Time 
        call StartComputeTime (ObjComputeTime, InitialSystemTime, BeginTime, EndTime, DT, VariableDT, MaxDT, STAT = STAT_CALL)

        !Update Current Time
        CurrentTime  = BeginTime

        !Constructs Basin
        call ConstructBasin   (ObjBasinID = ObjBasin, ObjTime = ObjComputeTime, ModelName = ModelName, STAT = STAT_CALL)

        ModelConstructed = .true.

    end subroutine ConstructMohidLand

    !--------------------------------------------------------------------------

    subroutine ReadKeywords

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(PathLength)                       :: DataFile
        integer                                     :: STAT_CALL
        integer                                     :: ObjEnterData = 0
        integer                                     :: iflag
        character(PathLength)                       :: WatchFile, DTLogFile

        !Monitor Performance of the model execution?
        call ReadFileName('OUTWATCH', WatchFile, Message = 'Start Watch File', STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            STAT_CALL = CreateWatchGroup (WatchFile)
            MonitorPerformance  = .true.
        else
            MonitorPerformance  = .false.
        endif

        !Monitor Performance of the model execution?
        call ReadFileName('OUTPUT_RESULT', OutputFile, Message = 'Start Result File', STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            SaveOutput  = .true.
        else
            SaveOutput  = .false.
        endif

        call ReadFileName('DT_LOG', DTLogFile, Message = 'Start DTLog File', STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then            
            MonitorDT  = .true.
        else
            MonitorDT  = .false.
        endif

        if (MonitorDT) then
            call UnitsManager (UnitDT, OPEN_FILE)      
            open(UNIT   = UnitDT, FILE   = DTLogFile, STATUS  = "UNKNOWN", IOSTAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLand - ERR01'
            write (UnitDT, '(A25, A10,A12,A12,A13,A13,A13,A13,A26)') "ModuleName", "iter", "DT", "DNet", "RunOff", "PorousMedia", "Atmosphere", "DTNextEv", "NextTime"
        end if

        call ReadFileName('IN_MODEL', DataFile, "Mohid Land Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLand - ERR02'

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLand - ERR03'

        call ReadTimeKeyWords   (ObjEnterData, FromFile, BeginTime, EndTime, DT,         &
                                 VariableDT, "Mohid Land", MaxDT, GmtReference,          &
                                 DTPredictionInterval)
                                 
        !Model Name
        call GetData(ModelName,                                                          &
                     ObjEnterData, iflag,                                                &
                     SearchType   = FromFile,                                            &
                     keyword      = 'MODEL_NAME',                                        &
                     default      = 'MOHID Land Model',                                  & 
                     ClientModule = 'MOHIDLand',                                         &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLand - ERR04'


        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLand - ERR05'

    end subroutine ReadKeywords

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyMohidLand 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        !integer                                     :: STAT_CALL
        !real                                        :: CPUTime, LastCPUTime = 0.
        logical                                      :: one_more = .true.

#ifndef _OUTPUT_OFF_
        write(*, *)"-------------------------- MOHID -------------------------"
        write(*, *)
        write(*, *)"Running MOHID Land, please wait..."
        write(*, *)                    
#endif

        do while (one_more)

            one_more = DoOneTimeStep()
            
        enddo

    end subroutine ModifyMohidLand

    !--------------------------------------------------------------------------

    logical function  DoOneTimeStep ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                                        :: NewDT
        integer                                     :: STAT_CALL
        real                                        :: CPUTime, LastCPUTime = 0.

        !Actualize the CurrentTime with Model time interval DT
        call ActualizeCurrentTime (TimeID    = ObjComputeTime,      &
                                   DT_Global = DT,                  &
                                   STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DoOneTimeStep - MohidLand - ERR01'
       
        !Gives the actualized Current time
        call GetComputeCurrentTime(ObjComputeTime, CurrentTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DoOneTimeStep - MohidLand - ERR02'
        
        if (VariableDT) then
                            
            NewDT = min(DT * 1.50, MaxDT)
        
        else
            
            NewDT = DT

        end if
                    
        call ModifyBasin(ObjBasin, NewDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DoOneTimeStep - MohidLand - ERR03'
    
        if (VariableDT) then
                              
            DT = min(NewDT, MaxDT)

!            !Rounds new DT to full decimal second
!            if (DT * 10.0 > AINT(DT*10.0)) then
!               DT = AINT(DT*10.0) + 1.0
!            else
!               DT = max(AINT(DT*10.0), 1.0)
!            endif
!
!            DT = DT / 10.0

            !Fit last Iteration
            if (CurrentTime >= EndTime) then
                    DoOneTimeStep = .false.
                    return            
            elseif (CurrentTime + DT > EndTime) then
                DT = EndTime - CurrentTime
!                if (abs(DT) < 1e-5) then
!                    DoOneTimeStep = .false.
!                    return
!                endif 
            else
                if ((EndTime - (CurrentTime + DT)) < 1e-5) then
                    DT = EndTime - CurrentTime 
                endif
            endif
                
            call ActualizeDT(TimeID = ObjComputeTime, DT = DT, STAT = STAT_CALL)     
            if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidLand - MohidLand - ERR04'
            
        else

            !Fit last Iteration
            if (CurrentTime .GE. EndTime) then
                DoOneTimeStep = .false.
                return
            endif

        endif 

        call CPU_TIME(CPUTime)
        if (CPUTime - LastCPUTime > DTPredictionInterval) then
            LastCPUTime = CPUTime
            call PrintProgress(ObjComputeTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidLand - MohidLand - ERR05'
        endif        

        DoOneTimeStep = .true.
    
    end function


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillMohidLand

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        real                                        :: ElapsedSeconds, TotalCPUTime


        call KillBasin(ObjBasinID = ObjBasin, STAT = STAT_CALL)

        if (MonitorPerformance) then
            call KillWatchGroup (STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillMohidLand - MohidLand - ERR01'
        endif

        if (MonitorDT) call UnitsManager (UnitDT, CLOSE_FILE)

        call date_and_time(Values = F95Time)
        call SetDate      (FinalSystemTime,   float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)
        call cpu_time(TotalCPUTime)
        ElapsedSeconds = FinalSystemTime - InitialSystemTime

        if (SaveOutput) call SaveRunInfo ("Mohid Land", ElapsedSeconds, TotalCPUTime, OutputFile)
        
        call ShutdownMohid ("Mohid Land", ElapsedSeconds, TotalCPUTime)        
        
        ModelConstructed = .false.

    end subroutine KillMohidLand

#ifdef _OPENMI_

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::Initialize
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_INITIALIZE"::Initialize
    !DEC$ ENDIF
    logical function Initialize(workingDirectory)
                     
        !Arguments-------------------------------------------------------------
        character(*)                                :: workingDirectory
        
        !Local-----------------------------------------------------------------
        
!        call SetError(WARNING_, INTERNAL_, "Test Error Message")
!        Initialize = .false.
!        return
        
        FilesName = workingDirectory
        
        call ConstructMohidLand()

        Initialize = .true.

        return
    
    end function Initialize
    
    !--------------------------------------------------------------------------

    !Perform a single time step
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::PerformTimeStep
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_PERFORMTIMESTEP"::PerformTimeStep
    !DEC$ ENDIF
    logical function PerformTimeStep()

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        logical                                     :: dummy
        
        dummy = DoOneTimeStep()
        PerformTimeStep = .true.

    end function PerformTimeStep
    
    !--------------------------------------------------------------------------

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::Finish
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_FINISH"::Finish
    !DEC$ ENDIF
    logical function Finish()

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------

        if (ModelConstructed) then
            call KillMohidLand()
        endif
        
        Finish = .true.

    end function Finish

    !--------------------------------------------------------------------------

!!!    !DEC$ IFDEFINED (VF66)
!!!    !dec$ attributes dllexport::Dispose
!!!    !DEC$ ELSE
!!!    !dec$ attributes dllexport,alias:"_DISPOSE"::Dispose
!!!    !DEC$ ENDIF
!!!    !The dispose function does not do anything. All Clean up is done by de Finish function
!!!    logical function Dispose()
!!!
!!!        !Arguments-------------------------------------------------------------
!!!        
!!!        !Local-----------------------------------------------------------------
!!!
!!!        Dispose = .true.
!!!
!!!    end function Dispose
    
    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetModelID
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETMODELID"::GetModelID
    !DEC$ ENDIF
    logical function GetModelID(id)
    
        !Arguments-------------------------------------------------------------
        character(*)                                :: id       
    
        id = ModelName
        GetModelID = .true.
        return
    
    end function GetModelID

    !--------------------------------------------------------------------------

    !Test Function - Runs the whole model
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::RunSimulation
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_RUNSIMULATION"::RunSimulation
    !DEC$ ENDIF
    !Test method to run the whole simulation once
    logical function RunSimulation()

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------

        call ModifyMohidLand
        
        RunSimulation = .true.
    
    end function RunSimulation

    !--------------------------------------------------------------------------

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetNumberOfMessages
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETNUMBEROFMESSAGES"::GetNumberOfMessages
    !DEC$ ENDIF
    !Return the number of Error Messages
    integer function GetNumberOfMessages()
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------

        GetNumberOfMessages = NumberOfErrorMessages
        
        return
    
    end function GetNumberOfMessages

    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetMessage
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETMESSAGE"::GetMessage
    !DEC$ ENDIF
    logical function GetMessage(Number, Message)

        !Arguments-------------------------------------------------------------
        integer                                     :: Number
        character(len=*)                            :: Message        
        !Local-----------------------------------------------------------------


        if(Number .ge. 1 .and. Number .le. MaxErrorMessages)then
            Message=ErrorMessagesStack(Number)
            GetMessage=.true.
        else
            Message=' '
            GetMessage=.false.
        endif

      end function GetMessage
      
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetStartInstant
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETSTARTINSTANT"::GetStartInstant
    !DEC$ ENDIF
    logical function GetStartInstant(Instant)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: Instant        

        !Local-----------------------------------------------------------------

        Instant = ConvertTimeToString(BeginTime)
        
        GetStartInstant = .true.

      end function GetStartInstant      

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetStopInstant
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETSTOPINSTANT"::GetStopInstant
    !DEC$ ENDIF
    logical function GetStopInstant(Instant)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: Instant        

        !Local-----------------------------------------------------------------

        Instant = ConvertTimeToString(EndTime)
        
        GetStopInstant = .true.

      end function GetStopInstant      

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetCurrentInstant
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETCURRENTINSTANT"::GetCurrentInstant
    !DEC$ ENDIF
    logical function GetCurrentInstant(Instant)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: Instant        

        !Local-----------------------------------------------------------------

        Instant = ConvertTimeToString(CurrentTime)
        
        GetCurrentInstant = .true.

      end function GetCurrentInstant      


    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetCurrentTimeStep
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETCURRENTTIMESTEP"::GetCurrentTimeStep
    !DEC$ ENDIF
    real(8) function GetCurrentTimeStep()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        GetCurrentTimeStep = dble(DT)
        

      end function GetCurrentTimeStep      


    
#endif


#ifdef _OPENMI_
end module MohidLand
#else
end program MohidLand
#endif

