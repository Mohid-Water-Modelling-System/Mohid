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

program MohidLand

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
    real                                :: DT, MaxDT
    logical                             :: VariableDT

    !Other Stuff
    type (T_Time)                       :: InitialSystemTime, FinalSystemTime
    integer, dimension(8)               :: F95Time

    call ConstructMohidLand
    call ModifyMohidLand
    call KillMohidLand

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructMohidLand

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

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
        call StartComputeTime (ObjComputeTime, BeginTime, EndTime, DT, VariableDT, MaxDT, STAT = STAT_CALL)

        call ConstructBasin   (ObjBasinID = ObjBasin, ObjTime = ObjComputeTime, STAT = STAT_CALL)

    end subroutine ConstructMohidLand

    !--------------------------------------------------------------------------

    subroutine ReadKeywords

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(PathLength)                       :: DataFile
        integer                                     :: STAT_CALL
        integer                                     :: ObjEnterData = 0
        character(PathLength)                       :: WatchFile, DTLogFile

        !Monitor Performance of the model execution?
        call ReadFileName('OUTWATCH', WatchFile, Message = 'Start Watch File', STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            STAT_CALL = CreateWatchGroup (WatchFile)
            MonitorPerformance  = .true.
        else
            MonitorPerformance  = .false.
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
        end if

        call ReadFileName('IN_MODEL', DataFile, "Mohid Land Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLand - ERR02'

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLand - ERR03'

        call ReadTimeKeyWords   (ObjEnterData, FromFile, BeginTime, EndTime, DT,         &
                                 VariableDT, "Mohid Land", MaxDT)

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLand - ERR04'

    end subroutine ReadKeywords

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyMohidLand 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        real                                        :: CPUTime, LastCPUTime = 0.
        real                                        :: NewDT


        write(*, *)"-------------------------- MOHID -------------------------"
        write(*, *)
        write(*, *)"Running MOHID Land, please wait..."
        write(*, *)                    

        CurrentTime  = BeginTime

        do

            !Actualize the CurrentTime with Model time interval DT
            call ActualizeCurrentTime (TimeID    = ObjComputeTime,      &
                                       DT_Global = DT,                  &
                                       STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidLand - MohidLand - ERR01'
           
            !Gives the actualized Current time
            call GetComputeCurrentTime(ObjComputeTime, CurrentTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidLand - MohidLand - ERR02'
            
            
            if (VariableDT) then
                                
                NewDT = min(DT * 1.50, MaxDT)
            
            else
                
                NewDT = DT

            end if
                        
            call ModifyBasin(ObjBasin, NewDT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidLand - MohidLand - ERR03'
            
            if (VariableDT) then
                              
                DT = min(NewDT, MaxDT)

                !Fit last Iteration
                if (CurrentTime + DT > EndTime) then
                    DT = EndTime - CurrentTime
                    if (abs(DT) < 1e-5) exit 
                endif

                call ActualizeDT(TimeID = ObjComputeTime, DT = DT, STAT = STAT_CALL)     
                if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidLand - MohidLand - ERR04'
            
            else

                !Fit last Iteration
                if (CurrentTime .GE. EndTime) exit                
    

            endif 

            call CPU_TIME(CPUTime)
            if (CPUTime - LastCPUTime > 60.) then
                LastCPUTime = CPUTime
                call PrintProgress(ObjComputeTime, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidLand - MohidLand - ERR05'
            endif

        enddo

    end subroutine ModifyMohidLand


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

        call ShutdownMohid ("Mohid Land", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMohidLand


end program MohidLand
