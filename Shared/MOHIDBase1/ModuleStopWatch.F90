!-------------------------------------------------------------------------
!        IST/MARETEC, Marine Modelling Group, Mohid2000 modelling system
!-------------------------------------------------------------------------
!BOI
! !TITLE: Mohid2000 Lagrangian model 
! !AUTHORS: Paulo Chambel, Frank Braunschweig
! !AFFILIATION: IST/MARETEC, Marine Modelling Group
! !DATE: Jul2001
! !INTRODUCTION: Time Measure Module
!
!EOI
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

Module ModuleStopWatch

!BOP
!
! !MODULE: ModuleStopWatch

!    !DESCRIPTION: 
!     This model is responsible for measuring the cpu time spend in different modules
 
! !REVISION HISTORY: 
!    Fev2002   Frank Braunschweig  New Implementation
!
! !FILES USED:   
!    USE_WATCH   : data file
!              
!
! !SEE ALSO:    
!  http://194.65.82.103/manuais/default.htm     
!
    use ModuleGlobalData
    use ModuleTime

    implicit none 

    private

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: CreateWatchGroup

    !Modifier
    public  :: StartWatch
    public  :: StopWatch

    !Destructor
    public  :: KillWatchGroup

    !Parameter
    integer, parameter                              :: STOPPED = 1
    integer, parameter                              :: RUNNING = 2

    !Defines one watch
    type T_Watch
        integer                                     :: State
        character(StringLength)                     :: ModuleName
        character(StringLength)                     :: RoutineName
        real                                        :: CpuElapsed
        real                                        :: WallElapsed
        integer, dimension(8)                       :: LastWallStart        
        integer, dimension(8)                       :: LastWallStop
        real                                        :: LastCPUStart
        real                                        :: LastCPUStop
        type (T_Watch), pointer                     :: Next
    end type T_Watch

    !Defines the list of watchgroups
    type T_WatchGroup
        type (T_Watch), pointer                     :: FirstWatch
        character (StringLength)                    :: OutputFileName
    end type T_WatchGroup

    type (T_WatchGroup), pointer                     :: WatchGroup

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    integer function CreateWatchGroup (OutputFile)
        
        !Arguments-------------------------------------------------------------
        character (len=*)                           :: OutputFile

        !Local-----------------------------------------------------------------
        
        nullify  (WatchGroup)
        allocate (WatchGroup)
        nullify  (WatchGroup%FirstWatch)

        WatchGroup%OutputFileName = OutputFile

        CreateWatchGroup = SUCCESS_

    end function CreateWatchGroup


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartWatch (ModuleName, RoutineName)

        !Arguments-------------------------------------------------------------
        character (len=*)                           :: ModuleName
        character (len=*)                           :: RoutineName

        !Local-----------------------------------------------------------------
        type (T_Watch), pointer                     :: NewWatch
        type (T_Watch), pointer                     :: CurrentWatch, PreviousWatch

        !Searches Watch
        CurrentWatch => WatchGroup%FirstWatch
        do while (associated (CurrentWatch))
            
            if (trim(adjustl(CurrentWatch%ModuleName )) == trim(adjustl(ModuleName )) .and.                &
                trim(adjustl(CurrentWatch%RoutineName)) == trim(adjustl(RoutineName)) ) exit
            CurrentWatch  => CurrentWatch%Next
        enddo

        !If the Watch doesnt exists, create a new one
        if (.not. associated (CurrentWatch)) then

            allocate (NewWatch)

            !Inits variables
            NewWatch%CpuElapsed  = 0.
            NewWatch%WallElapsed = 0.
            NewWatch%State       = STOPPED
            NewWatch%ModuleName  = ModuleName
            NewWatch%RoutineName = RoutineName
            nullify (NewWatch%Next)

            !Searches for right position in list
            nullify (PreviousWatch)
            CurrentWatch => WatchGroup%FirstWatch
            do while (associated(CurrentWatch))
                PreviousWatch => CurrentWatch
                CurrentWatch  => CurrentWatch%Next
            enddo

            !Inserts to list
            if (associated (PreviousWatch)) then
                PreviousWatch%Next => NewWatch
            else
                WatchGroup%FirstWatch => NewWatch
            endif

            CurrentWatch => NewWatch

        endif

        !If the current watch is running, something is wrong
        if (CurrentWatch%State == RUNNING) then
            write (*,*)'Watch Module  :', CurrentWatch%ModuleName
            write (*,*)'Watch Routine :', CurrentWatch%RoutineName
            write (*,*)'ALREADY RUNNING'
            stop 'ModuleStopWatch - StartWatch - ERR02'
        endif

        call CPU_TIME       (CurrentWatch%LastCPUStart)
        call DATE_AND_TIME  (values = CurrentWatch%LastWallStart)
        CurrentWatch%State = RUNNING

    end subroutine StartWatch

    !--------------------------------------------------------------------------

    subroutine StopWatch (ModuleName, RoutineName)

        !Arguments-------------------------------------------------------------
        character (len=*)                           :: ModuleName
        character (len=*)                           :: RoutineName

        !Local-----------------------------------------------------------------
        type (T_Watch), pointer                     :: CurrentWatch
        type (T_Time)                               :: StartDate    
        type (T_Time)                               :: EndDate

        !Searches Watch
        CurrentWatch => WatchGroup%FirstWatch
        do while (associated (CurrentWatch))
            if (trim(adjustl(CurrentWatch%ModuleName )) == trim(adjustl(ModuleName )) .and.                &
                trim(adjustl(CurrentWatch%RoutineName)) == trim(adjustl(RoutineName)) ) exit
            CurrentWatch  => CurrentWatch%Next
        enddo

        if (.not. associated (CurrentWatch)) then
            write (*,*)'Watch Module  :', trim(adjustl(ModuleName ))
            write (*,*)'Watch Routine :', trim(adjustl(RoutineName))
            write (*,*)'NOT ASSOCIATED'
            stop 'ModuleStopWatch - StopWatch - ERR01'
        else
            if (CurrentWatch%State == STOPPED) then
                write (*,*)'Watch Module  :', trim(adjustl(ModuleName ))
                write (*,*)'Watch Routine :', trim(adjustl(RoutineName))
                write (*,*)'ALREADY STOPPED'
                stop 'ModuleStopWatch - StopWatch - ERR02'
            endif

            call CPU_TIME       (CurrentWatch%LastCPUStop)
            call DATE_AND_TIME  (values = CurrentWatch%LastWallStop)

            CurrentWatch%CpuElapsed = CurrentWatch%CpuElapsed +                          &
                                     (CurrentWatch%LastCPUStop - CurrentWatch%LastCPUStart)

            call SetDate (StartDate, float(CurrentWatch%LastWallStart(1)),               &
                                     float(CurrentWatch%LastWallStart(2)),               &
                                     float(CurrentWatch%LastWallStart(3)),               &
                                     float(CurrentWatch%LastWallStart(5)),               &
                                     float(CurrentWatch%LastWallStart(6)),               &
                                     float(CurrentWatch%LastWallStart(7)) +              &
                                     float(CurrentWatch%LastWallStart(8))/1000.)

            call SetDate (EndDate,   float(CurrentWatch%LastWallStop(1)),                &
                                     float(CurrentWatch%LastWallStop(2)),                &
                                     float(CurrentWatch%LastWallStop(3)),                &
                                     float(CurrentWatch%LastWallStop(5)),                &
                                     float(CurrentWatch%LastWallStop(6)),                &
                                     float(CurrentWatch%LastWallStop(7)) +               &
                                     float(CurrentWatch%LastWallStop(8))/1000.)

            CurrentWatch%WallElapsed = CurrentWatch%WallElapsed +                        &
                                       (EndDate - StartDate)

            CurrentWatch%State = STOPPED

        endif

    end subroutine StopWatch

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine KillWatchGroup (STAT_CALL)

        !Arguments-------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: Unit
        type (T_Watch), pointer                     :: CurrentWatch, KillWatch
        
        call UnitsManager (Unit, OPEN_FILE, STAT = STAT_CALL)
        open (unit = unit, file = WatchGroup%OutputFileName, status = "unknown")

        write (unit = unit, fmt = 10)"Module Name", "Routine Name", "CPU Elapsed", "Wall Elapsed"


        CurrentWatch => WatchGroup%FirstWatch
        do while (associated (CurrentWatch))

            write (unit = unit, fmt = 20)trim(adjustl(CurrentWatch%ModuleName )),        &
                                         trim(adjustl(CurrentWatch%RoutineName)),        &
                                         CurrentWatch%CpuElapsed,                        &
                                         CurrentWatch%WallElapsed
            CurrentWatch  => CurrentWatch%Next
        enddo

        call UnitsManager (Unit, CLOSE_FILE, STAT = STAT_CALL)

        CurrentWatch => WatchGroup%FirstWatch
        do while (associated (CurrentWatch))
            KillWatch     => CurrentWatch
            CurrentWatch  => CurrentWatch%Next
            deallocate (KillWatch)
            nullify    (KillWatch)
        enddo

        deallocate (WatchGroup)

        STAT_CALL = SUCCESS_

 10 format (1x, a30, 1x, a30, 1x, a12,   1x, a12)
 20 format (1x, a30, 1x, a30, 1x, f12.3, 1x, f12.3)

    end subroutine KillWatchGroup


end Module ModuleStopWatch

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
