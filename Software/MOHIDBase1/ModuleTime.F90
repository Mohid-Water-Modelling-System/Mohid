!------------------------------------------------------------------------------
!        IST/MARETEC, Marine Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE: Mohid Lagrangian model 
! AUTHORS: Paulo Chambel, Frank Braunschweig
! AFFILIATION: IST/MARETEC, Marine Modelling Group
! DATE: May2001
! DESCRIPTION: Global Data / Constants for the whole MOHID Framework
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

Module ModuleTime

    use ModuleGlobalData

    implicit none

    private

    !Constructor
    public  :: StartComputeTime
    private ::      AllocateInstance


    !Selector
    public  :: GetComputeTimeStep
    public  :: GetMaxComputeTimeStep
    public  :: GetComputeCurrentTime
    public  :: GetComputeTimeLimits
    public  :: GetVariableDT
    public  :: GetGmtReference
    public  :: GetBackTracking
    
    !Modifier
    public  :: ActualizeCurrentTime
    public  :: ActualizeDT
    
    !Destructor
    public  :: KillComputeTime
    private ::      DeallocateInstance

    !Management
    private ::      Ready

    !Subroutines & Functions---------------------------------------------------

    public  :: GregorianDayToDate   !Day 1 -> 1/1/0000
    public  :: DateToGregorianDay   !Inverts Gregorian Day
    public  :: NumberOfDaysInMonth  !Return number of days in Month, including for Leap Years
    public  :: JulianDay            !Julday = 1 -> 1st of January
    public  :: JulianDayToMonthDay  !Convert JulianDay into Month and Day of month
    private :: Calendario
    public  :: ExtractDate
    public  :: SetDate
    public  :: TimeHours            !Function -> computes number of hours in a day
    public  :: null_time            !Turns type(time) to FillValueInt
    public  :: PrintProgress        !Writes a message to the screen         !Frank 3-8-99
    public  :: ConvertTimeToString  !Converts T_Time to a String like 2000:01:01:23:59:59


    !Operator------------------------------------------------------------------

    public  :: operator(.LT.) 
    public  :: operator(.LE.)  
    public  :: operator(.GT.)    
    public  :: operator(.GE.)   
    public  :: operator(.EQ.)      
    public  :: operator(.NE.)    

    public  :: operator(+)
    public  :: operator(-) 

    public  :: assignment(=)

    !Parameter-----------------------------------------------------------------

    integer, parameter :: Year_   = 1
    integer, parameter :: Month_  = 2
    integer, parameter :: Day_    = 3
    integer, parameter :: Hour_   = 4
    integer, parameter :: Minute_ = 5
    integer, parameter :: Second_ = 6


    !Type----------------------------------------------------------------------
    public :: T_Time
    type      T_Time
        private
        real, dimension(6) :: Time_ = FillValueReal
    end type T_Time            

    public :: T_OutPutTime
    type      T_OutPutTime
        type (T_Time), dimension(:), pointer     :: OutTime
        integer                                  :: NextOutPut, Number
        logical                                  :: ON
        integer                                  :: KLB, KUB, JLB, JUB, ILB, IUB
    end type T_OutPutTime 

    type      T_ComputeTime
        integer                         :: InstanceID
        type(T_Time   )                 :: InitialSystemTime
        type(T_Time   )                 :: Begin
        type(T_Time   )                 :: Finish
        type(T_Time   )                 :: Current
        real                            :: DT                   !Time step
        real                            :: MaxDT                !Maximum Allowed Time Step (when variable DT)
        real                            :: MinDTSoFar           !Min Time step During Simulation - Output only
        real                            :: MaxDTSoFar           !Max Time step During Simulation - Output only
        real                            :: nIter                !Number of iteration (to calc. average DT)
        real                            :: IterSinceLastPrint   !Number of iteration between prints
        real                            :: LastCpuTime          !CPU Time at last print
        logical                         :: VariableDT
        real                            :: GmtReference         !Time zone from GMT. Ex: Lisbon :: GMT +0; 
                                                                !Madrid :: GMT +1; NYC :: GMT -5
        logical                         :: Backtracking                                                                
        type(T_ComputeTime), pointer    :: Next
    end type T_ComputeTime

    !Interface-----------------------------------------------------------------

    private :: SetDateReal
    private :: SetDateInteger
    interface  SetDate
        module procedure SetDateReal
        module procedure SetDateInteger
    end interface  SetDate

    private :: TimeEarlierThan
    interface operator (.LT.)
        module procedure TimeEarlierThan
    end interface


    private :: TimeEarlierEqualThan
    interface  operator (.LE.)
        module procedure TimeEarlierEqualThan
    end interface

    private :: TimeLaterThan
    interface  operator (.GT.)
        module procedure TimeLaterThan
    end interface


    private :: TimeLaterEqualThan
    interface  operator (.GE.)
        module procedure TimeLaterEqualThan
    end interface


    private :: TimeEqualLogical
    interface  operator (.EQ.)
        module procedure TimeEqualLogical
    end interface


    private :: TimeNotEqualLogical
    interface  operator (.NE.)
        module procedure TimeNotEqualLogical
    end interface

    interface  assignment (=)
        module procedure TimeEqualAssignment
    end interface


    private :: TimePlusDTReal
    private :: TimePlusDTReal4
    private :: TimePlusDTReal8
    private :: TimePlusDTinteger
    interface  operator (+)
        module procedure TimePlusDTReal4
        module procedure TimePlusDTReal8
        module procedure TimePlusDTinteger
    end interface


    private :: TimeLess
    interface  operator (-)
        module procedure TimeLess
    end interface

    !Global Variables
    type (T_ComputeTime), pointer                   :: Me, FirstComputeTime

    !----------------------------------------------------------------------------

    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartComputeTime(TimeID, InitialSystemTime, BeginTime, EndTime,                              &
                                DT, VariableDT, MaxDT, GmtReference, BackTracking, STAT)
                      
        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID
        type(T_Time),  intent(IN)                   :: InitialSystemTime
        type(T_Time),  intent(IN)                   :: BeginTime, EndTime
        real,          intent(IN)                   :: DT
        logical,       intent(IN)                   :: VariableDT
        real,    optional, intent(IN)               :: MaxDT
        real,    optional, intent(IN)               :: GmtReference  
        logical, optional, intent(IN)               :: BackTracking
        integer, optional, intent(OUT)              :: STAT  

        !Local-----------------------------------------------------------------
        integer                                     :: ready_
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mTime_)) then
            nullify (FirstComputeTime)
            call RegisterModule (mTime_) 
        endif

        call Ready(TimeID, ready_)

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance       ()

            Me%InitialSystemTime   = InitialSystemTime
            Me%Finish              = EndTime
            Me%Begin               = BeginTime
            Me%Current             = BeginTime
            Me%DT                  = DT
            Me%VariableDT          = VariableDT

            if (present(MaxDT)) then
                Me%MaxDT           = MaxDT
            else
                Me%MaxDT           = null_real
            endif

            if (present(BackTracking)) then
                Me%BackTracking    = BackTracking
            else
                Me%BackTracking    = .false.
            endif

            if (present(GmtReference)) then
                Me%GmtReference           = GmtReference
            else
                Me%GmtReference           = 0. !Default is the Greenwhich Time.
            endif

            Me%MinDTSoFar          = DT
            Me%MaxDTSoFar          = DT
            Me%nIter               = 0
            Me%IterSinceLastPrint  = 0

!            call GetSystemTime(Me%LastSystemTime)
            call cpu_time(Me%LastCpuTime)

            !Returns ID
            TimeID     = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0

            stop 'ModuleTime - StartComputeTime - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartComputeTime

    !--------------------------------------------------------------------------

    subroutine AllocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_ComputeTime), pointer               :: NewTime
        type (T_ComputeTime), pointer               :: PreviousTime


        !Allocates new instance
        allocate (NewTime)
        nullify  (NewTime%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstComputeTime)) then
            FirstComputeTime    => NewTime
            Me                  => NewTime
        else
            PreviousTime        => FirstComputeTime
            Me                  => FirstComputeTime%Next
            do while (associated(Me))
                PreviousTime    => Me
                Me              => Me%Next
            enddo
            Me                  => NewTime
            PreviousTime%Next   => NewTime
        endif

        Me%InstanceID = RegisterNewInstance (mTIME_)

    end subroutine AllocateInstance

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetComputeTimeStep(TimeID, DT, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID
        real,              intent(OUT)              :: DT
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            DT = Me%DT

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetComputeTimeStep

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetMaxComputeTimeStep(TimeID, MaxDT, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID
        real,              intent(OUT)              :: MaxDT
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            MaxDT = Me%MaxDT

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetMaxComputeTimeStep

    !--------------------------------------------------------------------------

    subroutine GetComputeCurrentTime(TimeID, CurrentTime, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID
        type(T_Time),      intent(OUT)              :: CurrentTime
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            CurrentTime = Me%Current

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetComputeCurrentTime

    !--------------------------------------------------------------------------

    subroutine GetComputeTimeLimits(TimeID, EndTime, BeginTime, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID
        type(T_Time), optional                      :: EndTime
        type(T_Time), optional                      :: BeginTime
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(BeginTime)) BeginTime = Me%Begin
            if (present(EndTime  )) EndTime   = Me%Finish

            if (.not.present(BeginTime).and..not.present(EndTime))            &
                stop 'Sub. GetComputeTimeLimits - ModuleTime -  Error01'

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetComputeTimeLimits

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine GetVariableDT(TimeID, VariableDT, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID
        logical, optional, intent(OUT)              :: VariableDT
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(VariableDT)) VariableDT = Me%VariableDT

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetVariableDT

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetGmtReference(TimeID, GmtReference, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID
        real, optional, intent(OUT)                 :: GmtReference
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(GmtReference)) GmtReference = Me%GmtReference

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetGmtReference

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetBacktracking(TimeID, Backtracking, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID
        logical,        intent(OUT)                 :: Backtracking
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Backtracking = Me%Backtracking

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetBacktracking

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine ActualizeCurrentTime(TimeID, DT_Global, Current, STAT)
        
        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID
        real,                   intent(IN )         :: DT_Global
        type(T_Time), optional, intent(IN )         :: Current
        integer,      optional, intent(OUT)         :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeID, ready_) 


cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            if (present(Current)) then
            
                !Actualize the instant time
                Me%Current   = Current
            
            else

                Me%Current   = Me%Current + DT_Global
                
            endif

            !Increase number of iterations
            Me%nIter                = Me%nIter              + 1
            Me%IterSinceLastPrint   = Me%IterSinceLastPrint + 1
          
            STAT_ = SUCCESS_
        else               
            STAT_ = UNKNOWN_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ActualizeCurrentTime

    !--------------------------------------------------------------------------

    subroutine ActualizeDT(TimeID, DT, STAT)
        
        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID
        real,              intent (IN)              :: DT
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeID, ready_) 

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            !Actualize the time for the next iteration
            Me%DT                   = DT

            !Verifies current Min/Max
            if (Me%DT < Me%MinDTSoFar) Me%MinDTSoFar = Me%DT
            if (Me%DT > Me%MaxDTSoFar) Me%MaxDTSoFar = Me%DT


            STAT_ = SUCCESS_
        else               
            STAT_ = UNKNOWN_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ActualizeDT

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillComputeTime(TimeID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID
        integer, optional, intent(OUT)              :: STAT    

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        integer                                     :: nUsers

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeID, ready_) 
        
cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mTIME_,  Me%InstanceID)

            if (nUsers == 0) then
 
                call null_time(Me%Begin)
                call null_time(Me%Finish)
                call null_time(Me%Current)

                Me%DT =FillValueReal

                !Deallocates Instance
                call DeallocateInstance ()
                
                TimeID = 0

                STAT_ = SUCCESS_

            end if

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine KillComputeTime


    !--------------------------------------------------------------------------

    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_ComputeTime), pointer               :: AuxTime
        type (T_ComputeTime), pointer               :: PreviousTime

        !Updates pointers
        if (Me%InstanceID == FirstComputeTime%InstanceID) then
            FirstComputeTime => FirstComputeTime%Next
        else
            PreviousTime => FirstComputeTime
            AuxTime      => FirstComputeTime%Next
            do while (AuxTime%InstanceID /= Me%InstanceID)
                PreviousTime => AuxTime
                AuxTime      => AuxTime%Next
            enddo

            !Now update linked list
            PreviousTime%Next => AuxTime%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine Ready (TimeID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (TimeID > 0) then
            call LocateMe (TimeID)
            ready_ = VerifyReadLock (mTIME_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateMe (TimeID)

        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID

        !Local-----------------------------------------------------------------

        Me => FirstComputeTime
        do while (associated (Me))
            if (Me%InstanceID == TimeID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleTime - LocateMe - ERR01'

    end subroutine LocateMe

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !OPERATOR OPERATOR OPERATOR OPERATOR OPERATOR OPERATOR OPERATOR OPERATOR OP 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    logical function TimeEarlierThan(Time1, Time2)

        !Arguments-------------------------------------------------------------
        type(T_Time), intent(IN)                    :: Time1, Time2

        !Local-----------------------------------------------------------------
        integer                                     :: I

        !----------------------------------------------------------------------

        TimeEarlierThan = .FALSE.

do1 :   do I = 1, 6
cd1 :       if      (Time1%Time_(I) < Time2%Time_(I)) then
                TimeEarlierThan = .TRUE.
                return
            else if (Time1%Time_(I) > Time2%Time_(I)) then
                return
            end if cd1
        end do do1

        !----------------------------------------------------------------------

    end function  TimeEarlierThan

    !--------------------------------------------------------------------------

    logical function TimeLaterThan(Time1, Time2)

        !Arguments-------------------------------------------------------------
        type(T_Time), intent(IN)                    :: Time1, Time2

        !Local-----------------------------------------------------------------
        integer                                     :: I

        !----------------------------------------------------------------------

        TimeLaterThan = .false.

do1 :   do i = 1, 6
cd1 :       if      (Time1%Time_(I) > Time2%Time_(I)) then
                TimeLaterThan = .TRUE.
                return
            else if (Time1%Time_(I) < Time2%Time_(I)) then
                return
            end if cd1
        end do do1

        !----------------------------------------------------------------------

    end function  TimeLaterThan

    !--------------------------------------------------------------------------

    logical function TimeLaterEqualThan(Time1, Time2)

        !Arguments-------------------------------------------------------------
        type(T_Time), intent(IN)                    :: Time1, Time2

        !Local-----------------------------------------------------------------

cd1 :   if      (Time1 .GT. Time2) then
            TimeLaterEqualThan = .TRUE.
        else if (Time1 .EQ. Time2) then
            TimeLaterEqualThan = .TRUE.
        else
            TimeLaterEqualThan = .FALSE.
        end if cd1
             
        !----------------------------------------------------------------------

    end function  TimeLaterEqualThan

    !--------------------------------------------------------------------------

    logical function TimeEarlierEqualThan(Time1, Time2)

        !Arguments-------------------------------------------------------------
        type(T_Time), intent(IN)                    :: Time1, Time2

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

cd1 :   if      (Time1 .LT. Time2) then
            TimeEarlierEqualThan = .TRUE.
        else if (Time1 .EQ. Time2) then
            TimeEarlierEqualThan = .TRUE.
        else
            TimeEarlierEqualThan = .FALSE.
        end if cd1

        !----------------------------------------------------------------------

    end function  TimeEarlierEqualThan

    !--------------------------------------------------------------------------

    subroutine TimeEqualAssignment(Time1, Time2)

        !Arguments-------------------------------------------------------------
        type(T_Time), intent(OUT)                   :: Time1
        type(T_Time), intent(IN)                    :: Time2

        !Local-----------------------------------------------------------------
        integer                                     :: I

        !----------------------------------------------------------------------

do1 :   do I = 1, 6
            Time1%Time_(I) = Time2%Time_(I)
        end do do1

        !----------------------------------------------------------------------

    end subroutine TimeEqualAssignment

    !--------------------------------------------------------------------------

    logical function TimeEqualLogical(Time1, Time2)

        !Arguments-------------------------------------------------------------
        type(T_Time), intent(IN)                    :: Time1, Time2

        !Local-----------------------------------------------------------------
        integer                                     :: I

        !----------------------------------------------------------------------

        TimeEqualLogical = .TRUE.

do1 :   do I = 1, 6
            if (Time1%Time_(I) .NE. Time2%Time_(I)) &
                TimeEqualLogical = .FALSE.
        end do do1

        !----------------------------------------------------------------------

    end function TimeEqualLogical

    !--------------------------------------------------------------------------

    logical function  TimeNotEqualLogical(Time1, Time2)

        !Arguments-------------------------------------------------------------
        type(T_Time), intent(IN)                    :: Time1, Time2

        !Local-----------------------------------------------------------------
        integer                                     :: I

        !----------------------------------------------------------------------

        TimeNotEqualLogical = .FALSE.

do1 :   do i = 1, 6
            if (Time1%Time_(I) .NE. Time2%Time_(I))                           &
                TimeNotEqualLogical = .TRUE.
        end do do1

        !----------------------------------------------------------------------

    end function  TimeNotEqualLogical

    !--------------------------------------------------------------------------

    function TimePlusDTreal(Time1, DT)

        type(T_Time) :: TimePlusDTreal

        !Arguments-------------------------------------------------------------
        type(T_Time), intent(IN)                    :: Time1
        real,         intent(IN)                    :: DT

        !Local-----------------------------------------------------------------
        integer                                     :: GregDay2, GregDay1, Days
        real                                        :: DT_seconds
        type(T_Time)                                :: TimeAux
        real                                        :: Hour, Minutes, Seconds

        !----------------------------------------------------------------------

        TimeAux = Time1

        if (DT >= 86400.0) then

            Days = int(DT/86400.0)

            call DateToGregorianDay(TimeAux,GregDay1)

            Hour    = TimeAux%Time_(4)
            Minutes = TimeAux%Time_(5)
            Seconds = TimeAux%Time_(6)


            GregDay2 = GregDay1 + Days 
            call GregorianDayToDate(GregDay2, TimeAux)

            TimeAux%Time_(4)  = Hour
            TimeAux%Time_(5)  = Minutes
            TimeAux%Time_(6)  = Seconds

            DT_seconds = DT - real(Days) * 86400.0
        else 
            DT_seconds = DT
        end if 
     
        TimeAux%Time_(6) = TimeAux%Time_(6) + DT_seconds 

        call Calendario(TimeAux)   
        TimePlusDTreal = TimeAux

        !----------------------------------------------------------------------

    end function TimePlusDTreal

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    function TimePlusDTReal8(Time1, DT)

        type(T_Time) :: TimePlusDTReal8

        !Arguments-------------------------------------------------------------
        type(T_Time), intent(IN)                    :: Time1
        real(8),      intent(IN)                    :: DT

        !Local-----------------------------------------------------------------
        integer                                     :: GregDay2, GregDay1, Days
        real(8)                                     :: DT_seconds
        type(T_Time)                                :: TimeAux
        real                                        :: Hour, Minutes, Seconds

        !----------------------------------------------------------------------

        TimeAux = Time1

        if (DT >= 86400.0) then

            Days = int(DT/86400.0)

            call DateToGregorianDay(TimeAux,GregDay1)

            Hour    = TimeAux%Time_(4)
            Minutes = TimeAux%Time_(5)
            Seconds = TimeAux%Time_(6)


            GregDay2 = GregDay1 + Days 
            call GregorianDayToDate(GregDay2, TimeAux)

            TimeAux%Time_(4)  = Hour
            TimeAux%Time_(5)  = Minutes
            TimeAux%Time_(6)  = Seconds

            DT_seconds = DT - real(Days) * 86400.0
        else 
            DT_seconds = DT
        end if 
     
        TimeAux%Time_(6) = TimeAux%Time_(6) + DT_seconds 

        call Calendario(TimeAux)   
        TimePlusDTReal8 = TimeAux

        !----------------------------------------------------------------------

    end function TimePlusDTReal8

    !--------------------------------------------------------------------------

    function TimePlusDTinteger(Time1, DT)

        type(T_Time) :: TimePlusDTinteger

        !Arguments-------------------------------------------------------------
        type(T_Time), intent(IN)                    :: Time1
        integer,      intent(IN)                    :: DT

        !Local-----------------------------------------------------------------
        real                                        :: DT_seconds

        !----------------------------------------------------------------------

        DT_seconds          = DT
        TimePlusDTinteger   = TimePlusDTreal(Time1, DT_seconds)

        !----------------------------------------------------------------------

    end function TimePlusDTinteger

    !--------------------------------------------------------------------------

    function TimePlusDTReal4(Time1, DT)

        type(T_Time) :: TimePlusDTReal4

        !Arguments-------------------------------------------------------------
        type(T_Time), intent(IN)                    :: Time1
        real(4),      intent(IN)                    :: DT

        !Local-----------------------------------------------------------------
        real                                        :: DT_seconds
        !----------------------------------------------------------------------

        DT_seconds      = DT
        TimePlusDTReal4 = TimePlusDTreal(Time1, DT_seconds)

        !----------------------------------------------------------------------

    end function TimePlusDTReal4

    !--------------------------------------------------------------------------

    real(8) function TimeLess(Time1, Time2)

        !Arguments-------------------------------------------------------------
        type(T_Time), intent(IN)                    :: Time1, Time2

        !Local-----------------------------------------------------------------
        real(8)                                     :: Aux
        integer                                     :: GregDay1, GregDay2
        
        !----------------------------------------------------------------------

        call DateToGregorianDay(Time1,GregDay1)
        call DateToGregorianDay(Time2,GregDay2)

        Aux =       (dble(GregDay1      )-dble(GregDay2      ))*86400.
        Aux = Aux + (dble(Time1%Time_(4))-dble(Time2%Time_(4)))*3600.
        Aux = Aux + (dble(Time1%Time_(5))-dble(Time2%Time_(5)))*60.
        Aux = Aux + (dble(Time1%Time_(6))-dble(Time2%Time_(6)))
        
        TimeLess = Aux
       
        !----------------------------------------------------------------------

    end function TimeLess

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SUBROUTINE SUBROUTINE SUBROUTINE SUBROUTINE SUBROUTINE SUBROUTINE SUBROUTI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine DateToGregorianDay(TimeIn, GregDay)
    
        !Arguments-------------------------------------------------------------
        type (T_Time), intent(IN)                   :: TimeIn
        integer, intent(OUT)                        :: GregDay

        !Local-----------------------------------------------------------------
        real                                        :: Year, Month, Day, hour, minute, second
        integer                                     :: n400, n100, n4, nx
        
        !Gets Date in format YY MM DD hh mm ss
        call ExtractDate(TimeIn, Year, Month, Day, hour, minute, second)

        !Number of 400 yrs period
        n400 = int(Year/400)
        
        !Number of 100 yrs period
        n100 = int(Year/100) - n400
        
        if (mod(nint(Year),100) == 0 .and. mod(nint(Year), 400) /= 0) then
            n100 = n100 - 1
        endif
        
        !Number of 4 yrs period
        n4   = int(Year/4) - n400 - n100
        
        if (mod(nint(Year),4) == 0) then
            n4 = n4 - 1
        endif
        
        !Number of other years
        nx   = Year - n400 - n100 - n4
        
        !Total number of days until last day of last century
        !Sum one at the end because year 0 was a leap year
        GregDay= (n400 + n4) * 366 + (n100 + nx) * 365 + 1
        
        !Add days until beginning end of last month
        GregDay = GregDay + DaysInYearBeforeMonth(nint(Year), nint(Month))
        
        !Add days of current months
        GregDay = GregDay + Day

    end subroutine DateToGregorianDay
    
    !--------------------------------------------------------------------------
    
    subroutine GregorianDayToDate (GreorgianDay, TimeOut)

        !Arguments-------------------------------------------------------------
        type (T_Time), intent(OUT)                  :: TimeOut
        integer, intent(IN)                         :: GreorgianDay

        !Local-----------------------------------------------------------------
        integer                                     :: GregDay, i

        real                                        :: Year, Month, Day
        integer                                     :: n400, n100, n4, nx
        integer                                     :: DaysInMonth

        !Local Var - intent in attribute
        GregDay = GreorgianDay

        !-1 one since year zero was a leap year... consistent with previous version
        GregDay = GregDay - 1
        
        !Number of 400 year periods
        !A 400 year period has 366 + 3*365 + (100-4)*366 + (400-100)* 365 days = 146097 days
        n400 = int(GregDay/146097)
        
        !400 years period 
        GregDay = GregDay - n400 * 146097

        if (GregDay == 0) then
            Year  = n400 * 400
            Month = 1
            Day   = 1
            call SetDate(TimeOut, Year, Month, Day, 0., 0., 0.)
            return
        endif
        
       
        
        !Number of 100 year periods in reminder
        !A 100 year period has 365 + (25-1)*366 + (100-25)* 365 days = 36524 days
        n100 = int(GregDay/36524)

        GregDay = GregDay - n100 * 36524
        
        if (GregDay == 0) then
            Year  = n400 * 400 + n100 * 100 - 1
            Month = 12
            Day   = 31
            call SetDate(TimeOut, Year, Month, Day, 0., 0., 0.)
            return
        endif
        
        !Number of 4 year periods in reminder
        !A 4 year period has 3*365 + 366 = 1461 days
        n4   = int(GregDay/1461)
        GregDay = GregDay - n4 * 1461
        
       
        !Number of years in reminder
        nx   = int(GregDay/365)
        
        if (mod(GregDay,365) == 0 .and. nx > 0) then
            nx = nx - 1
        endif
        
        GregDay = GregDay - nx * 365
        
        !Calculates Year
        Year = n400 * 400 + n100 * 100 + n4 * 4 + nx
        
        if (IsLeapYear(nint(Year))) then
            GregDay = GregDay +1
        endif

       
        !Calculates Month
do1:    do i=1, 12
            DaysInMonth = NumberOfDaysInMonth(nint(Year), i)
            if (GregDay <= DaysInMonth) then
                exit do1
            else
                GregDay = GregDay - DaysInMonth
            endif
        enddo do1
        
        Month   = i
        Day     = GregDay
        
        call SetDate(TimeOut, Year, Month, Day, 0., 0., 0.)
       
    end subroutine GregorianDayToDate
    
    !--------------------------------------------------------------------------

    logical function IsLeapYear (Year)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: Year

        
        if ( (mod(Year,4) == 0 .and. mod(Year,100) /= 0) .or. (mod(Year,400) == 0))  then
            IsLeapYear =.true.
        else
            IsLeapYear = .false.
        endif

    
    
    end function IsLeapYear

    !--------------------------------------------------------------------------

    integer function NumberOfDaysInMonth(Year, Month)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: Year, Month
        
        !Local-----------------------------------------------------------------
        integer                                     :: NDM(12)

        DATA NDM/31,28,31,30,31,30,31,31,30,31,30,31/
        
        if (Month /= 2) then
            NumberOfDaysInMonth = NDM(Month)
        else
            if (IsLeapYear(Year)) then
                NumberOfDaysInMonth = 29
            else
                NumberOfDaysInMonth = 28
            endif
        endif
            
    end function NumberOfDaysInMonth
    
    !--------------------------------------------------------------------------

    integer function DaysInYearBeforeMonth(Year, Month)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: Year, Month
        
        !Local-----------------------------------------------------------------
        integer                                     :: NDP(12)

        DATA NDP/0,31,59,90,120,151,181,212,243,273,304,334/
        
        DaysInYearBeforeMonth = NDP(Month)
        
        if (IsLeapYear(Year) .and. Month > 2) then
            DaysInYearBeforeMonth = DaysInYearBeforeMonth + 1
        endif
    
    end function DaysInYearBeforeMonth      

    !--------------------------------------------------------------------------

    subroutine JulianDay(Time1, JulDay)

        !Arguments-------------------------------------------------------------

        type(T_Time), intent(IN) :: Time1

        integer :: JulDay

        !External--------------------------------------------------------------

        integer            :: GregDay1, GregDay2
        real               :: Year
        real,    parameter :: Month  = 1
        real,    parameter :: Day    = 1
        real,    parameter :: Hour   = 0
        real,    parameter :: Minute = 0
        real,    parameter :: Second = 0

        !Local-----------------------------------------------------------------

        type(T_Time) :: Time2

        !----------------------------------------------------------------------


        call ExtractDate(Time1, Year = Year)
        call SetDate    (Time2, Year, Month, Day, Hour, Minute, Second)

        call DateToGregorianDay(Time1, GregDay1)
        call DateToGregorianDay(Time2, GregDay2)

        JulDay = GregDay1 - GregDay2 + 1

        !----------------------------------------------------------------------

    end subroutine JulianDay

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine JulianDayToMonthDay(year, JulDay, TimeMonthDay)

        !Arguments-------------------------------------------------------------

        type(T_Time), intent(OUT) :: TimeMonthDay

        integer, intent(IN)       :: JulDay
        integer, intent(IN)       :: year

        !External--------------------------------------------------------------

        integer            :: GregDay

        !Local-----------------------------------------------------------------

        !Get a reference GregorianDay 
        !for the first day of the first month
        call SetDate    (TimeMonthDay, year, 1, 1, 0, 0, 0)
        call DateToGregorianDay(TimeMonthDay, GregDay)

        !Convert the JulDay to a valid GregDay
        GregDay = GregDay + JulDay - 1

        !Get the corresponding Month and Day of the Julian Day
        call GregorianDayToDate(GregDay, TimeMonthDay)

        !----------------------------------------------------------------------

    end subroutine JulianDayToMonthDay

    !--------------------------------------------------------------------------

    !Esta subroutina actualiza a data apos um incremento.   

    subroutine Calendario(Time1)

        !Arguments-------------------------------------------------------------
        type (T_Time), intent(InOut)                :: Time1

        !Local-----------------------------------------------------------------
        integer, Dimension(12)                      :: NDay
        integer                                     :: J
        integer                                     :: nTimes

        DATA NDay/31,28,31,30,31,30,31,31,30,31,30,31/

    
        !Leap year ?
        if (IsLeapYear(nint(Time1%Time_(1)))) then
            NDay(2)=29
        else
            NDay(2)=28
        endif


        !Adjusts Seconds
        if     (Time1%Time_(6).GE.60.) then
            nTimes = int(Time1%Time_(6)/ 60.)
            Time1%Time_(6) = Time1%Time_(6) - nTimes * 60.
            Time1%Time_(5) = Time1%Time_(5) + nTimes
        elseif (Time1%Time_(6).LT. 0.) then
            nTimes = -int(Time1%Time_(6)/ 60.) + 1
            Time1%Time_(6) = Time1%Time_(6) + nTimes * 60.
            Time1%Time_(5) = Time1%Time_(5) - nTimes
            if (Time1%Time_(6) == 60.0) then
                Time1%Time_(6) = 0.0
                Time1%Time_(5) = Time1%Time_(5) + 1
            endif
        endif

        !Adjusts Minutes
        if     (Time1%Time_(5).GE.60.) then
            nTimes = int(Time1%Time_(5)/ 60.)
            Time1%Time_(5) = Time1%Time_(5) - nTimes * 60.
            Time1%Time_(4) = Time1%Time_(4) + nTimes
        elseif (Time1%Time_(5).LT. 0.) then
            nTimes = -int(Time1%Time_(5)/ 60.) + 1
            Time1%Time_(5) = Time1%Time_(5) + nTimes * 60.
            Time1%Time_(4) = Time1%Time_(4) - nTimes
            if (Time1%Time_(5) == 60.0) then
                Time1%Time_(5) = 0.0
                Time1%Time_(4) = Time1%Time_(4) + 1
            endif
        endif


        !Adjusts Hours
        if     (Time1%Time_(4).GE.24.) then
            nTimes = int(Time1%Time_(4)/ 24.)
            Time1%Time_(4) = Time1%Time_(4) - nTimes * 24.
            Time1%Time_(3) = Time1%Time_(3) + nTimes
        elseif (Time1%Time_(4).LT. 0.) then
            nTimes = -int(Time1%Time_(4)/ 24.) + 1
            Time1%Time_(4) = Time1%Time_(4) + nTimes * 24.
            Time1%Time_(3) = Time1%Time_(3) - nTimes
            if (Time1%Time_(4) == 24.0) then
                Time1%Time_(4) = 0.0
                Time1%Time_(3) = Time1%Time_(3) + 1
            endif
        endif

        !
        J=int(Time1%Time_(2))
cd1 :   if ((J .LT. 1) .OR. (J .GT. 12)) then
            write(*,*) 
            write(*,*) 'Date out of bounds. Month = ', J,'.'
            stop       'Subroutine Calendario; Module Time. ERR01.'
        end if cd1
        IF (int(Time1%Time_(3)).GT.NDay(J)) Then
           Time1%Time_(3)=Time1%Time_(3)-real(NDay(J))
           Time1%Time_(2)=Time1%Time_(2)+1.
        EndIf
        IF (Time1%Time_(3).LT.1.) Then
           J = int(Time1%Time_(2))-1
           IF (J.EQ.0) J = 12
           Time1%Time_(3)=Time1%Time_(3)+real(NDay(J))
           Time1%Time_(2)=Time1%Time_(2)-1.
        EndIf

        IF (Time1%Time_(2).GT.12.) Then
           Time1%Time_(2)=Time1%Time_(2)-12
           Time1%Time_(1)=Time1%Time_(1)+1
        EndIf

        IF (Time1%Time_(2).LT.1.) THEN
           Time1%Time_(2)=Time1%Time_(2)+12
           Time1%Time_(1)=Time1%Time_(1)-1
        EndIf

        !----------------------------------------------------------------------

    end subroutine Calendario

    !--------------------------------------------------------------------------







    !--------------------------------------------------------------------------

    subroutine ExtractDate(Time1, Year, Month, Day, Hour, Minute, Second)

        !Arguments-------------------------------------------------------------

        real, optional, intent(OUT) :: Year, Month, Day, Hour, Minute, Second

        type(T_Time),      intent(IN ) :: Time1

        !----------------------------------------------------------------------

        if (present(Year  )) Year   = Time1%Time_(Year_  )
        if (present(Month )) Month  = Time1%Time_(Month_ )
        if (present(Day   )) Day    = Time1%Time_(Day_   )
        if (present(Hour  )) Hour   = Time1%Time_(Hour_  )
        if (present(Minute)) Minute = Time1%Time_(Minute_)
        if (present(Second)) Second = Time1%Time_(Second_)

        !----------------------------------------------------------------------

    end subroutine ExtractDate

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine SetDateInteger(Time1, Year, Month, Day, Hour, Minute, Second)

        !Arguments-------------------------------------------------------------

        integer,      intent(IN ) :: Year, Month, Day, Hour, Minute, Second

        type(T_Time), intent(OUT) :: Time1

        !----------------------------------------------------------------------

        Time1%Time_(Year_  ) = Year  
        Time1%Time_(Month_ ) = Month 
        Time1%Time_(Day_   ) = Day   
        Time1%Time_(Hour_  ) = Hour  
        Time1%Time_(Minute_) = Minute
        Time1%Time_(Second_) = Second

        !----------------------------------------------------------------------

    end subroutine SetDateInteger

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine SetDateReal(Time1, Year, Month, Day, Hour, Minute, Second)

        !Arguments-------------------------------------------------------------

        Real,         intent(IN ) :: Year, Month, Day, Hour, Minute, Second

        type(T_Time), intent(OUT) :: Time1

        !----------------------------------------------------------------------

        Time1%Time_(Year_  ) = Year  
        Time1%Time_(Month_ ) = Month 
        Time1%Time_(Day_   ) = Day   
        Time1%Time_(Hour_  ) = Hour  
        Time1%Time_(Minute_) = Minute
        Time1%Time_(Second_) = Second

        !----------------------------------------------------------------------

    end subroutine SetDateReal

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine null_time(Time1)

        !Arguments-------------------------------------------------------------

        type(T_Time), intent(OUT) :: Time1

        !----------------------------------------------------------------------

        Time1%Time_(Year_  ) = FillValueInt  
        Time1%Time_(Month_ ) = FillValueInt 
        Time1%Time_(Day_   ) = FillValueInt   
        Time1%Time_(Hour_  ) = FillValueInt  
        Time1%Time_(Minute_) = FillValueInt
        Time1%Time_(Second_) = FillValueInt

        !----------------------------------------------------------------------

    end subroutine null_time

    !--------------------------------------------------------------------------




    !--------------------------------------------------------------------------

    function TimeHours(Time1)

        real :: TimeHours

        !Arguments-------------------------------------------------------------

        type(T_Time), intent(INOUT) :: Time1

        !Local-----------------------------------------------------------------

        real :: TimeHours_ 

        !----------------------------------------------------------------------

        TimeHours_ = Time1%Time_(4) + Time1%Time_(5) / 60.0 + Time1%Time_(6) / 3600.0

cd1 :   if (TimeHours_ .GT. 23.99999999) then
            TimeHours_ = 0.0
        end if cd1

        TimeHours = TimeHours_


        !----------------------------------------------------------------------

    end function TimeHours

    !--------------------------------------------------------------------------

    subroutine PrintProgress(TimeID, STAT)
                      
        !Arguments-------------------------------------------------------------
        integer                                     :: TimeID
        integer, optional, intent(OUT)              :: STAT    

        !Local-----------------------------------------------------------------
        integer                                     :: ready_
        integer                                     :: STAT_             
        real(4)                                     :: StillToRun, ExecutionTime
        real(4)                                     :: TotalExecutionTime
        real                                        :: TimeSimulated, TotalSimulationTime
        real                                        :: Coeficient
        integer, dimension(8)                       :: F95Time
        type(T_Time)                                :: CurrentSystemTime, FinishSystemTime

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TimeID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            !Time in seconds of the simulation
            TimeSimulated       = Me%Current - Me%Begin
            TotalSimulationTime = Me%Finish  - Me%Begin                        
            
            !Gets the current system time
            call DATE_AND_TIME(Values = F95Time)
            call SetDate(CurrentSystemTime, F95Time(1), F95Time(2), F95Time(3), F95Time(5), F95Time(6), F95Time(7))

            !!Gets the CPU time used until now
            !call CPU_Time(ExecutionTime)
            !CPU execution time gets biased when multi-threading.
            !It's safer to estimate based on the initialSystemTime
            ExecutionTime = CurrentSystemTime - Me%InitialSystemTime

            if (TimeSimulated > 0.) then

                !Calculates a preview total ExecutionTime
                TotalExecutionTime = TotalSimulationTime / TimeSimulated * ExecutionTime

                !Calculates the time which is still missing to end the simulation
                StillToRun = TotalExecutionTime - ExecutionTime

                !Calculates the reason CPU / Model
                Coeficient = TotalExecutionTime / TotalSimulationTime

                !System Time when the model will stop
                FinishSystemTime = CurrentSystemTime + StillToRun

            endif

            !Prints Messages to the Screen
#ifndef _OUTPUT_OFF_
            write(*, *)"-----Current Simulation Instant---------------------------"
            write(*,100)int(Me%Current%Time_(1)), int(Me%Current%Time_(2)), &
                        int(Me%Current%Time_(3)), int(Me%Current%Time_(4)), &
                        int(Me%Current%Time_(5)), int(Me%Current%Time_(6))

            if (Me%VariableDT) then
                write(*,105)Me%DT
                write(*,106)Me%MinDTSoFar
                write(*,107)Me%MaxDTSoFar
                write(*,108)TimeSimulated / Me%nIter
            else
                write(*,*)
            endif

            if (TimeSimulated > 0.) then
                write(*, *)"-----CPU Time---------------------------------------------"
                write(*,110) int(ExecutionTime )
                write(*,120) int(StillToRun )
                write(*,125) (ExecutionTime / (ExecutionTime + StillToRun))*100.0
                write(*,130) Coeficient
                if (Me%IterSinceLastPrint /= 0) then
                    write(*,131)(ExecutionTime - Me%LastCpuTime) /Me%IterSinceLastPrint
                    Me%IterSinceLastPrint = 0
                    Me%LastCpuTime        = ExecutionTime
                endif
            endif

            write(*, *)"-----System Time------------------------------------------"
            write(*,140)int(CurrentSystemTime%Time_(1)), int(CurrentSystemTime%Time_(2)), &
                        int(CurrentSystemTime%Time_(3)), int(CurrentSystemTime%Time_(4)), &
                        int(CurrentSystemTime%Time_(5)), int(CurrentSystemTime%Time_(6))

            if (TimeSimulated > 0.) then
                write(*,150)int(FinishSystemTime%Time_(1)), int(FinishSystemTime%Time_(2)), &
                            int(FinishSystemTime%Time_(3)), int(FinishSystemTime%Time_(4)), &
                            int(FinishSystemTime%Time_(5)), int(FinishSystemTime%Time_(6))
            endif
#endif

            100 format(1X, "Time Instant           : ",(i4,":"),4(i2, ":"), i2)
            105 format(1X, "Time Step              : ",f12.2,"s")
            106 format(1X, "Min Time Step so far   : ",f12.2,"s")
            107 format(1X, "Max Time Step so far   : ",f12.2,"s")
            108 format(1X, "Average Time Step      : ",f12.2,"s",/)
            110 format(1X, "Elapsed                : ",i12,"s")
            120 format(1X, "Remaining (aprox.)     : ",i12,"s")
            125 format(1X, "Completed (%)          : ",f14.4)
            130 format(1x, "Coeficient CPU / Model : ",f14.4)
            131 format(1X, "Seconds per Iteration  : ",f14.4,"s")
            140 format(1x, "System time            : ",(i4,":"),4(i2, ":"), i2)
            150 format(1x, "End of the run         : ",(i4,":"),4(i2, ":"), i2,6/)
            160 format(1x, "Number of threads      : ", i12,/)

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine PrintProgress
    
    !--------------------------------------------------------------------------

    character(len=30) function ConvertTimeToString(Time1)
    
        !Arguments-------------------------------------------------------------
        type(T_Time),      intent(IN )              :: Time1

        !Local-----------------------------------------------------------------
        character(len=30)                           :: auxStr
        
        !----------------------------------------------------------------------

        write(auxStr, 10)int(Time1%Time_(1)), int(Time1%Time_(2)), int(Time1%Time_(3)),  &
                         int(Time1%Time_(4)), int(Time1%Time_(5)), Time1%Time_(6)


!        write(*,*)'TIME CONVERSIONS'
!        write(*,*)Time1%Time_(1), Time1%Time_(2), Time1%Time_(3), Time1%Time_(4), Time1%Time_(5), Time1%Time_(6)
!        write(*,*)auxStr
        
        ConvertTimeToString = auxStr
 

     10 format((i4,":"),4(i2, ":"), f12.8)


        !2000:12:12:23:59:59.123345678

    end function ConvertTimeToString         

    

end module ModuleTime

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------



