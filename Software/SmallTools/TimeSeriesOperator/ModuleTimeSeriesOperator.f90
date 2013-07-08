!------------------------------------------------------------------------------
!        HIDROMOD : Modelação em Engenharia
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Time Serie
! URL           : http://www.mohid.com
! AFFILIATION   : HIDROMOD
! DATE          : Nov2012
! REVISION      : Paulo Leitão - v1.0
! DESCRIPTION   : Module to do operations between time series (e.g. mass balance)
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

Module ModuleTimeSeriesOperator

    use ModuleGlobalData
    use ModuleFunctions
    use ModuleEnterData
    use ModuleTime
    use ModuleTimeSerie
    use nr; use nrtype; use nrutil

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructTimeSeriesOperator
    private ::      AllocateInstance

    !Selector
    public  :: GetTimeSeriesOperatorPointer
    public  :: GetTimeSeriesOperatorInteger
                     
    
    !Modifier
    public  :: ModifyTimeSeriesOperator

    !Destructor
    public  :: KillTimeSeriesOperator                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjTimeSeriesOperator 
    
    !Interfaces----------------------------------------------------------------


    !Parameter-----------------------------------------------------------------
    
    real(8), parameter  :: Pi_ = 3.1415926535897932384626433832795
    !Input / Output
    integer, parameter  :: FileOpen = 1, FileClose = 0
    

    integer, parameter  :: Day_ = 1, Week_ = 2, Month_ = 3    
    
    integer, parameter  :: WeekEnd_ = 1, WeekWay_ = 2
    


    
    !Types---------------------------------------------------------------------
    type T_TimeSeriesOperator
    
        integer                                                 :: InstanceID
    
        type(T_Time)                                            :: BeginTime, EndTime
        real                                                    :: StartNightHour, EndNightHour
        
        real                                                    :: DT
        logical                                                 :: VariableDT
        real                                                    :: LevelVolX, LevelVol

	    integer                                                 :: ObjEnterData          = 0
	    integer                                                 :: ObjTime               = 0
	    integer, dimension(:), pointer                          :: ObjTimeSerieInFlux    
	    integer, dimension(:), pointer                          :: ObjTimeSerieOutFlux   
	    integer, dimension(:), pointer                          :: ObjTimeSerieMass      
	    integer                                                 :: ObjTimeSerieOut       = 0	    
	    integer                                                 :: ObjTimeSerieOutNight  = 0      	    
	    
	    integer                                                 :: TimeSerieColumn       = 2
	    
	    integer                                                 :: NFluxIn, NFluxOut, NMass
	    logical                                                 :: FluxInON, FluxOutON, MassON
	    real                                                    :: GapLimit
	    logical                                                 :: FoundOneGap
	    real                                                    :: NotValidMask
	    real                                                    :: MinNightValuesRacio
	    
	    character(len=PathLength)                               :: OptionsFile = "TimeSeriesOperator.dat", OutPutFile
	    
        type(T_TimeSeriesOperator), pointer                     :: Next

    end type T_TimeSeriesOperator    

    !Global Variables
    type (T_TimeSeriesOperator), pointer                        :: FirstObjTimeSeriesOperator
    type (T_TimeSeriesOperator), pointer                        :: Me    


    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructTimeSeriesOperator(ObjTimeSeriesOperatorID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                             :: ObjTimeSeriesOperatorID 
        integer, optional, intent(OUT)                      :: STAT     

        !External----------------------------------------------------------------
        integer                                             :: ready_         

        !Local-------------------------------------------------------------------
        integer                                             :: STAT_, STAT_CALL
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        character(len=StringLength)                         :: Extension

        

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mTimeSeriesOperator_)) then
            nullify (FirstObjTimeSeriesOperator)
            call RegisterModule (mTimeSeriesOperator_) 
        endif

        call Ready(ObjTimeSeriesOperatorID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            call ReadKeywords
            
            call ConstructInputInFlux
           
            call ConstructInputOutFlux
            
            call ConstructInputMass            

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - ConstructTimeSeriesOperator - ERR10'   
            
            allocate(PropertyList(3)) 
            
            PropertyList(1) = "Flux"
            PropertyList(2) = "Flux_percentage"
            PropertyList(3) = "Flux_Filter_P50_Window"
            
            Extension       = " "

            call StartTimeSerie(TimeSerieID             = Me%ObjTimeSerieOut,           &
                                ObjTime                 = Me%ObjTime,                   &
                                TimeSerieDataFile       = Me%OptionsFile,               &
                                PropertyList            = PropertyList,                 &
                                Extension               = Extension,                    &
                                ResultFileName          = Me%OutPutFile,                &
                                HavePath                = .true.,                       &
                                STAT                    = STAT_CALL)       
            if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - ConstructTimeSeriesOperator - ERR20'
            
            deallocate(PropertyList) 
            
            allocate(PropertyList(2)) 
            
            PropertyList(1) = "Night_Average_Flux"
            
            PropertyList(2) = "Night_P50_Flux"            
            
            Extension       = " "

            call StartTimeSerie(TimeSerieID             = Me%ObjTimeSerieOutNight,      &
                                ObjTime                 = Me%ObjTime,                   &
                                TimeSerieDataFile       = Me%OptionsFile,               &
                                PropertyList            = PropertyList,                 &
                                Extension               = Extension,                    &
                                ResultFileName          = "Night_"//trim(Me%OutPutFile),&
                                HavePath                = .true.,                       &
                                STAT                    = STAT_CALL)       
            if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - ConstructTimeSeriesOperator - ERR30'

            deallocate(PropertyList) 
                     
            !Returns ID
            ObjTimeSeriesOperatorID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleTimeSeriesOperator - ConstructTimeSeriesOperator - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructTimeSeriesOperator
 
    !--------------------------------------------------------------------------

    subroutine ReadKeywords

        !Local--------------------------------------------------------------
        integer                                             :: status, flag, STAT_CALL
        
        !Begin--------------------------------------------------------------

    
        Me%ObjEnterData = 0
        
        call ConstructEnterData(Me%ObjEnterData, Me%OptionsFile, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - ReadKeywords - ERR10'
        
        
        
        call ReadTimeKeyWords(ObjEnterData = Me%ObjEnterData,                           &     
                              ExtractTime  = FromFile,                                  &
                              BeginTime    = Me%BeginTime,                              &
                              EndTime      = Me%EndTime,                                &
                              DT           = Me%DT,                                     &
                              VariableDT   = Me%VariableDT,                             &
                              ClientModule = "ModuleTimeSeriesOperator")

        call StartComputeTime(Me%ObjTime, Me%BeginTime, Me%BeginTime, Me%EndTime, Me%DT, &
                              Me%VariableDT, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - ReadKeywords - ERR20'

        
        call GetData(Me%OutputFile,                                                     &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='OUTPUT_FILE',                                       &
                     ClientModule ='ModuleTimeSeriesOperator',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - ReadKeywords - ERR30'
        if (flag == 0) then
            write(*,*) 'Needs the output file'
            stop 'ModuleTimeSeriesOperator - ReadKeywords - ERR40'
        endif

        call GetData(Me%GapLimit,                                                       &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='GAP_LIMIT',                                         &
                     ClientModule ='ModuleTimeSeriesOperator',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - ReadKeywords - ERR50'
        if (flag == 0) then
            write(*,*) 'Needs the time gap limit'
            stop 'ModuleTimeSeriesOperator - ReadKeywords - ERR60'
        endif

        call GetData(Me%NotValidMask,                                                   &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='NOT_VALID_MASK',                                    &
                     ClientModule ='ModuleTimeSeriesOperator',                          &
                     default      = FillValueReal,                                      &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - ReadKeywords - ERR70'


        call GetData(Me%StartNightHour,                                                 &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='START_NIGHT_HOUR',                                  &
                     ClientModule ='ModuleTimeSeriesOperator',                          &
                     default      = 0.,                                                 &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - ReadKeywords - ERR80'

        call GetData(Me%EndNightHour,                                                   &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='END_NIGHT_HOUR',                                    &
                     ClientModule ='ModuleTimeSeriesOperator',                          &
                     default      = 4.,                                                 &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - ReadKeywords - ERR90'

        
        if (Me%StartNightHour >= Me%EndNightHour) then
            Me%EndNightHour = Me%EndNightHour + 24.
        endif
        
        call GetData(Me%MinNightValuesRacio,                                            &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='MIN_NIGHT_VAlUES_RACIO',                            &
                     ClientModule ='ModuleTimeSeriesOperator',                          &
                     default      = 0.75,                                               &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - ReadKeywords - ERR110'
        

    end subroutine ReadKeywords
    
    !-------------------------------------------------------------------------
 
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_TimeSeriesOperator), pointer                         :: NewObjTimeSeriesOperator
        type (T_TimeSeriesOperator), pointer                         :: PreviousObjTimeSeriesOperator


        !Allocates new instance
        allocate (NewObjTimeSeriesOperator)
        nullify  (NewObjTimeSeriesOperator%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjTimeSeriesOperator)) then
            FirstObjTimeSeriesOperator         => NewObjTimeSeriesOperator
            Me                    => NewObjTimeSeriesOperator
        else
            PreviousObjTimeSeriesOperator      => FirstObjTimeSeriesOperator
            Me                    => FirstObjTimeSeriesOperator%Next
            do while (associated(Me))
                PreviousObjTimeSeriesOperator  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjTimeSeriesOperator
            PreviousObjTimeSeriesOperator%Next => NewObjTimeSeriesOperator
        endif

        Me%InstanceID = RegisterNewInstance (mTimeSeriesOperator_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ConstructInputInFlux

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                     :: GroupExist 
        integer                                     :: i, ClientNumber, line, iflag
        integer                                     :: FirstLine, LastLine, STAT_CALL
        logical                                     :: BlockFound
        character(len=PathLength)                   :: InputFile
      
        !Begin-----------------------------------------------------------------
        

        call ExtractBlockFromBuffer(Me%ObjEnterData,                                    &
                                    ClientNumber    = ClientNumber,                     &
                                    block_begin     = '<BeginInFlux>',                  &
                                    block_end       = '<EndInFlux>',                    &
                                    BlockFound      = BlockFound,                       &
                                    FirstLine       = FirstLine,                        &
                                    LastLine        = LastLine,                         &
                                    STAT            = STAT_CALL)
IS:     if (STAT_CALL == SUCCESS_) then

BF:         if (BlockFound) then

                Me%FluxInON = .true.
                 
                Me%NFluxIn = LastLine - FirstLine - 1
                
                allocate(Me%ObjTimeSerieInFlux(Me%NFluxIn))
                Me%ObjTimeSerieInFlux(1:Me%NFluxIn) = 0
                
                i=0
                do line=FirstLine +1, LastLine-1
                    i = i + 1
                    call GetData(InputFile, EnterDataID = Me%ObjEnterData, flag = iflag,     &
                                 Buffer_Line = line, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructInputInFlux - ModuleTimeSeriesOperator - ERR10'
                    if (iflag == 0) stop 'ConstructInputInFlux - ModuleTimeSeriesOperator - ERR20'
                    
                    
                    
                    call StartTimeSerieInput(TimeSerieID       = Me%ObjTimeSerieInFlux(i), &
                                             TimeSerieDataFile = InputFile,                &
                                             ObjTime           = Me%ObjTime,               &
                                             CheckDates        = .false.,                  &
                                             STAT              = STAT_CALL)
                enddo
            
            else BF
                Me%FluxInON = .false.
            endif BF
             
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)                  
        
        else IS
            
            stop 'ConstructInputInFlux - ModuleTimeSeriesOperator - ERR30'
        
        endif IS           
        
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructInputInFlux - ModuleTimeSeriesOperator - ERR40'
        
        
    end subroutine ConstructInputInFlux

    !--------------------------------------------------------------------------

    subroutine ConstructInputOutFlux

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                     :: GroupExist 
        integer                                     :: i, ClientNumber, line, iflag
        integer                                     :: FirstLine, LastLine, STAT_CALL
        logical                                     :: BlockFound
        character(len=PathLength)                   :: InputFile
      
        !Begin-----------------------------------------------------------------
        

        call ExtractBlockFromBuffer(Me%ObjEnterData,                                    &
                                    ClientNumber    = ClientNumber,                     &
                                    block_begin     = '<BeginOutFlux>',                 &
                                    block_end       = '<EndOutFlux>',                   &
                                    BlockFound      = BlockFound,                       &
                                    FirstLine       = FirstLine,                        &
                                    LastLine        = LastLine,                         &
                                    STAT            = STAT_CALL)
IS:     if (STAT_CALL == SUCCESS_) then

BF:         if (BlockFound) then

                Me%FluxOutON = .true.
                 
                Me%NFluxOut = LastLine - FirstLine - 1
                
                allocate(Me%ObjTimeSerieOutFlux(Me%NFluxOut))
                Me%ObjTimeSerieOutFlux       (1:Me%NFluxOut) = 0
                
                i=0
                do line=FirstLine +1, LastLine-1
                    i = i + 1
                    call GetData(InputFile, EnterDataID = Me%ObjEnterData, flag = iflag,   &
                                 Buffer_Line = line, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructInputOutFlux - ModuleTimeSeriesOperator - ERR10'
                    if (iflag == 0) stop 'ConstructInputOutFlux - ModuleTimeSeriesOperator - ERR20'
                    
                    
                    
                    call StartTimeSerieInput(TimeSerieID       = Me%ObjTimeSerieOutFlux(i),&
                                             TimeSerieDataFile = InputFile,                &
                                             ObjTime           = Me%ObjTime,               &
                                             CheckDates        = .false.,                  &
                                             STAT              = STAT_CALL)
                enddo
            
            else BF
                Me%FluxOutON = .false.
            endif BF
             
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)                  
        
        else IS
            
            stop 'ConstructInputOutFlux - ModuleTimeSeriesOperator - ERR30'
        
        endif IS           
        
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructInputOutFlux - ModuleTimeSeriesOperator - ERR40'
        
        
    end subroutine ConstructInputOutFlux


    !--------------------------------------------------------------------------

    subroutine ConstructInputMass

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                     :: GroupExist 
        integer                                     :: i, ClientNumber, line, iflag
        integer                                     :: FirstLine, LastLine, STAT_CALL
        logical                                     :: BlockFound
        character(len=PathLength)                   :: InputFile
      
        !Begin-----------------------------------------------------------------
        

        call ExtractBlockFromBuffer(Me%ObjEnterData,                                    &
                                    ClientNumber    = ClientNumber,                     &
                                    block_begin     = '<BeginMass>',                    &
                                    block_end       = '<EndMass>',                      &
                                    BlockFound      = BlockFound,                       &
                                    FirstLine       = FirstLine,                        &
                                    LastLine        = LastLine,                         &
                                    STAT            = STAT_CALL)
IS:     if (STAT_CALL == SUCCESS_) then

BF:         if (BlockFound) then

                Me%MassON = .true.
                 
                Me%NMass = LastLine - FirstLine - 1
                
                allocate(Me%ObjTimeSerieMass(Me%NMass))
                Me%ObjTimeSerieMass(1:Me%NMass) = 0
                
                i=0
                do line=FirstLine +1, LastLine-1
                    i = i + 1
                    call GetData(InputFile, EnterDataID = Me%ObjEnterData, flag = iflag,   &
                                 Buffer_Line = line, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructInputMass - ModuleTimeSeriesOperator - ERR10'
                    if (iflag == 0) stop 'ConstructInputMass - ModuleTimeSeriesOperator - ERR20'
                    
                    
                    
                    call StartTimeSerieInput(TimeSerieID       = Me%ObjTimeSerieMass(i),   &
                                             TimeSerieDataFile = InputFile,                &
                                             ObjTime           = Me%ObjTime,               &
                                             CheckDates        = .false.,                  &
                                             STAT              = STAT_CALL)
                enddo
            
            else BF
                Me%MassON = .false.
            endif BF
             
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)                  
        
        else IS
            
            stop 'ConstructInputMass - ModuleTimeSeriesOperator - ERR30'
        
        endif IS           
        
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructInputMass - ModuleTimeSeriesOperator - ERR40'
        
        if (Me%MassON) then

            call GetData(Me%LevelVolX,                                                  &
                         Me%ObjEnterData,                                               &
                         iflag,                                                         &
                         SearchType   = FromFile,                                       &
                         keyword      ='LEVEL_VOL_X',                                   &
                         default      = 1.,                                             &
                         ClientModule ='ModuleTimeSeriesOperator',                      &
                         STAT         = STAT_CALL)        
            if(STAT_CALL /= SUCCESS_) stop 'ConstructInputMass - ModuleTimeSeriesOperator - ERR50'
        
        
            call GetData(Me%LevelVol,                                                   &
                         Me%ObjEnterData,                                               &
                         iflag,                                                         &
                         SearchType   = FromFile,                                       &
                         keyword      ='LEVEL_VOL',                                     &
                         default      = 0.,                                             &
                         ClientModule ='ModuleTimeSeriesOperator',                      &
                         STAT         = STAT_CALL)        
            if(STAT_CALL /= SUCCESS_) stop 'ConstructInputMass - ModuleTimeSeriesOperator - ERR60'        
        
        endif
        
    end subroutine ConstructInputMass

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    
    !--------------------------------------------------------------------------
    subroutine GetTimeSeriesOperatorPointer (ObjTimeSeriesOperatorID, Matrix, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjTimeSeriesOperatorID
        real(8), dimension(:, :, :),  pointer           :: Matrix
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTimeSeriesOperatorID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mTimeSeriesOperator_, Me%InstanceID)

            !Matrix => Me%Matrix

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetTimeSeriesOperatorPointer
    
    !--------------------------------------------------------------------------
    
    subroutine GetTimeSeriesOperatorInteger (ObjTimeSeriesOperatorID, Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjTimeSeriesOperatorID
        real                                            :: Int
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTimeSeriesOperatorID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Int = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetTimeSeriesOperatorInteger

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyTimeSeriesOperator(ObjTimeSeriesOperatorID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTimeSeriesOperatorID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTimeSeriesOperatorID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            call ComputeTimeSeriesOperator

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyTimeSeriesOperator
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine ComputeTimeSeriesOperator


        !Arguments-------------------------------------------------------------

        
        !Local-----------------------------------------------------------------
        type (T_Time)                  :: FluxTime, MassTime, CurrentTime, StartNightTime, EndNightTime, OutNightTime
        real                           :: OutPutFlux, RelativeOutFlux, NewMass, OldMass, InFlux, Year, Month, Day, NightPeriod, NightPeriodValid
        real                           :: NightAvFlux, NightP50Flux
        real, dimension(:,:), pointer  :: TotalArray
        real, dimension(:  ), pointer  :: WriteAux, WriteAuxNight, SortArray, AuxSort
        integer                        :: STAT_CALL, iSort, iN, i50, iP50, jmin, jmax, j, k, iTotal, i, imax
        logical                        :: NightPeriodON, NightPeriodOut
        
        
        !Begin-----------------------------------------------------------------
        
        allocate(WriteAux     (3))
        allocate(WriteAuxNight(2))
        
        if (Me%MassON) then
            OldMass = MassValue(Me%BeginTime)
        endif

        FluxTime        = Me%BeginTime + Me%DT/2.
        MassTime        = Me%BeginTime + Me%DT
        
        NightPeriod      = (Me%EndNightHour - Me%StartNightHour) * 3600.
        NightPeriodON    = .false.
       
        NightPeriodOut   = .false.
        NightAvFlux      = 0.
        NightPeriodValid = 0.    
        
        iN               = int(NightPeriod / Me%DT) + 1
        iSort            = 0
        
        iTotal           = int((Me%EndTime - Me%BeginTime)/Me%DT) + 1
        iP50             = int(3600. /Me%DT / 2.) + 1
        
        allocate(TotalArray(iTotal,3))
        allocate(SortArray (iN      ))

        i = 0
        
        call ExtractDate(FluxTime, Year, Month, Day)           

        do while(MassTime <= Me%EndTime)
        
            i = i + 1
        
            call SetDate(StartNightTime, Year, Month, Day, Me%StartNightHour, 0., 0.)
            call SetDate(EndNightTime  , Year, Month, Day, Me%EndNightHour  , 0., 0.)
            call SetDate(OutNightTime  , Year, Month, Day, (Me%StartNightHour + Me%EndNightHour) / 2. , 0., 0.)            
        
            OutPutFlux      = 0.
            
            Me%FoundOneGap     = .false.
            
            if (Me%FluxInON) then
                InFlux          = FluxInValue(FluxTime)
                OutPutFlux      = InFlux
                
            endif
            
            if (Me%FluxOutON) then
                OutPutFlux      = OutPutFlux - FluxOutValue(FluxTime)
            endif
        
            if (Me%MassON) then
                NewMass         = MassValue(MassTime)
                !Center in time
                OutPutFlux      = OutPutFlux + (NewMass - OldMass) / Me%DT
                OldMass         = NewMass
            endif
            
            if (Me%FoundOneGap) then
            
                TotalArray(i,1) = Me%NotValidMask
                TotalArray(i,2) = Me%NotValidMask
                TotalArray(i,3) = Me%NotValidMask
            
            else

                TotalArray(i,1) = OutPutFlux
                
                if (Me%FluxInON .and. abs(InFlux)>0.) then
                    TotalArray(i,2)      = OutPutFlux / InFlux * 100.
                endif
            
            endif
            
            if (FluxTime > StartNightTime .and. FluxTime < EndNightTime) Then
                if (.not. Me%FoundOneGap) then
                    NightPeriodValid = NightPeriodValid + Me%DT
                    NightAvFlux      = NightAvFlux + OutPutFlux * Me%DT
                    NightPeriodON    = .true.
                    iSort            = iSort + 1
                    SortArray(iSort) = OutPutFlux
                    
                endif
            else
                if (NightPeriodON) then
                    if (NightPeriodValid/NightPeriod > Me%MinNightValuesRacio) then
                        NightPeriodOut = .true.
                        if (NightPeriodValid > 0.) then
                            NightAvFlux  = NightAvFlux / NightPeriodValid
                        endif
                        WriteAuxNight(1) = NightAvFlux
                        allocate(AuxSort(1:iSort))
                        AuxSort(1:iSort) = SortArray(1:iSort)
                        call sort(AuxSort)
                        i50              = int(iSort/2) + 1 
                        NightP50Flux     = AuxSort(i50)            
                        WriteAuxNight(2) = NightP50Flux
                        deallocate(AuxSort)
                    endif
                NightPeriodON    = .false.                                            
                endif
                iSort            = 0       
                call ExtractDate(FluxTime, Year, Month, Day)                   
                               
            endif
            
            if (NightPeriodOut) then

                call WriteTimeSerieLine(TimeSerieID         = Me%ObjTimeSerieOutNight,  &
                                        DataLine            = WriteAuxNight,            &
                                        ExternalCurrentTime = OutNightTime,             &
                                        STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeTimeSeriesOperator - ModuleHydrodynamic - ERR10'
            
                NightPeriodOut   = .false.
                NightAvFlux      = 0.
                NightPeriodValid = 0.   

            endif
            
            !Center in time
            FluxTime = FluxTime + Me%DT
            !Forward in time
            MassTime = MassTime + Me%DT
        
        enddo
        
        imax          = i
        MassTime      = Me%BeginTime + Me%DT
        FluxTime      = Me%BeginTime + Me%DT/2.        
        i             = 0

        allocate(AuxSort(iP50*2+2))

        do while(MassTime <= Me%EndTime)

            i = i + 1

            WriteAux(1:2) = TotalArray(i,1:2)
            
            jmin = max(1,   i-iP50)
            jmax = min(imax,i+iP50)
            k    = 0
            do j = jmin, jmax
            
                if (TotalArray(j,1) /= Me%NotValidMask) then
                
                    k          = k + 1
                    AuxSort(k) = TotalArray(j,1)
                
                endif
            
            enddo
            
            if (k > iP50 .and. jmax-jmin == iP50*2) then

                call sort(AuxSort(1:k))
                i50         = int(k/2) + 1 
                WriteAux(3) = AuxSort(i50)
            else
                
                WriteAux(3) = Me%NotValidMask
                
            endif
                            
            
            call WriteTimeSerieLine(TimeSerieID         = Me%ObjTimeSerieOut,           &
                                    DataLine            = WriteAux,                     &
                                    ExternalCurrentTime = FluxTime,                     &
                                    STAT                = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeTimeSeriesOperator - ModuleHydrodynamic - ERR10'
            

            !Forward in time
            MassTime = MassTime + Me%DT
            FluxTime = FluxTime + Me%DT

        enddo        
        
        deallocate(AuxSort      )
        
        deallocate(WriteAux     )
        deallocate(WriteAuxNight)       
        deallocate(SortArray    ) 
        deallocate(TotalArray   )
                
    end subroutine ComputeTimeSeriesOperator
    
    
   
    !--------------------------------------------------------------------------
    
    

    !--------------------------------------------------------------------------

    real function FluxInValue(Now)    
        !Arguments--------------------------------------------------------------
        type (T_Time)                                   :: Now
        !Local------------------------------------------------------------------
        real                                            :: FluxAux
        integer                                         :: i

        !Begin------------------------------------------------------------------

        FluxAux = 0.
        
        do i=1, Me%NFluxIn
            FluxAux = FluxAux + TimeSerieValue(Me%ObjTimeSerieInFlux(i), Now) 
        enddo

        FluxInValue = FluxAux
        

    end function FluxInValue

    !--------------------------------------------------------------------------
    
    
    !--------------------------------------------------------------------------

    real function FluxOutValue(Now)    
        !Arguments--------------------------------------------------------------
        type (T_Time)                                   :: Now
        !Local------------------------------------------------------------------
        real                                            :: FluxAux
        integer                                         :: i
        !Begin------------------------------------------------------------------

        FluxAux = 0.
        
        do i=1, Me%NFluxOut
            FluxAux = FluxAux + TimeSerieValue(Me%ObjTimeSerieOutFlux(i), Now) 
        enddo
        
        FluxOutValue = FluxAux
        

    end function FluxOutValue

    !--------------------------------------------------------------------------
        
    !--------------------------------------------------------------------------

    real function MassValue(Now)    
        !Arguments--------------------------------------------------------------
        type (T_Time)                                   :: Now
        !Local------------------------------------------------------------------
        real                                            :: MassAux
        integer                                         :: i
        !Begin------------------------------------------------------------------

        MassAux = 0.
        
        do i=1, Me%NMass
            MassAux = MassAux + TimeSerieValue(Me%ObjTimeSerieMass(i), Now) * Me%LevelVolX + Me%LevelVol
        enddo
        
        MassValue = MassAux
        

    end function MassValue

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    real function TimeSerieValue(ObjTimeSerie, Now)
        !Arguments--------------------------------------------------------------
        integer                                         :: ObjTimeSerie
        type (T_Time)                                   :: Now
        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL
        type (T_Time)                                   :: Time1, Time2
        real                                            :: Value1, Value2
        logical                                         :: TimeCycle, ValidInstant

        !Begin------------------------------------------------------------------

        call GetTimeSerieCycle(ObjTimeSerie, TimeCycle, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'TimeSerieValue - ModuleTimeSeriesOperator - ERR10'
     
        if (TimeCycle) then

            !Gets Value for current Time
            call GetTimeSerieValue (ObjTimeSerie, Now, Me%TimeSerieColumn,              &
                                    Time1, Value1, Time2, Value2, TimeCycle,            &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TimeSerieValue - ModuleTimeSeriesOperator - ERR20'

            TimeSerieValue = Value1

        else

            ValidInstant = GetTimeSerieCheckDate (ObjTimeSerie, Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TimeSerieValue - ModuleTimeSeriesOperator - ERR30'
            
            if (ValidInstant) then

                !Gets Value for current Time
                call GetTimeSerieValue (ObjTimeSerie, Now, Me%TimeSerieColumn,          &
                                        Time1, Value1, Time2, Value2, TimeCycle,        &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'TimeSerieValue - ModuleTimeSeriesOperator - ERR40'


                !Interpolates Value for current instant
                call InterpolateValueInTime(Now, Time1, Value1, Time2, Value2, TimeSerieValue)
                
                if ((Time2-Time1) > Me%GapLimit) then
                    Me%FoundOneGap = .true.
                endif
                
                if (Value1 == Me%NotValidMask .or. Value2 == Me%NotValidMask) then
                    Me%FoundOneGap = .true.
                    TimeSerieValue = Me%NotValidMask
                endif
                
            else
                Me%FoundOneGap = .true.
                TimeSerieValue = Me%NotValidMask
                
            endif                

        endif

    end function TimeSerieValue

    !--------------------------------------------------------------------------
    

    function FreqAnalysis(SortArray,SizeArray, Percentil)    

        !Arguments-----------------------------
        real :: Percentil, FreqAnalysis
        integer :: SizeArray
        real, dimension(:) :: SortArray
        !Local---------------------------------
        real       :: Aux, Raux
        integer    :: Iaux
        
        !Begin---------------------------------
        if (SizeArray==1) then
            FreqAnalysis=SortArray(1)
        else
            Aux = real(SizeArray-1)*Percentil
            Iaux = int(Aux)
            Raux = Aux-real(Iaux)
            if (Percentil == 0 .or. Iaux == 0) then
                FreqAnalysis = SortArray(1)
            else
                if (Raux>0.) then
                    FreqAnalysis = SortArray(Iaux+1)*Raux+SortArray(Iaux)*(1.-Raux)
                else
                    FreqAnalysis = SortArray(Iaux)
                endif
            endif
        endif
    
    
    end function FreqAnalysis
    
    
    logical function WeekWayDay (Date)
    
        !Arguments-----------------------------    
        type (T_Time) :: Date
           
        !Local---------------------------------
        type (T_Time) :: AuxDate
        real(8)       :: AuxSeconds, Weeks, AuxWeek 
        !Begin---------------------------------

        call SetDate (AuxDate, 2013.,2.,11.,0.,0.,0.)
        
        AuxSeconds = Date-AuxDate
        
        if (AuxSeconds >= 0.) then
            Weeks    = AuxSeconds / (168.*3600.)
            AuxWeek  = Weeks - int(Weeks)
            if (AuxWeek <5./7.) then
                WeekWayDay = .true.
            else
                WeekWayDay = .false.
            endif    
        else
            AuxSeconds = - AuxSeconds
            Weeks    = AuxSeconds / (168.*3600.)
            AuxWeek  = Weeks - int(Weeks)
            if (AuxWeek >2./7.) then
                WeekWayDay = .true.
            else
                WeekWayDay = .false.
            endif             
        endif 
    
    end function WeekWayDay 
    
    !----------------------------------------------------------------

    logical function WeekEndDay (Date)
    
        !Arguments-----------------------------    
        type (T_Time) :: Date
           
        !Local---------------------------------
        !Begin---------------------------------

        if (WeekWayDay(Date)) then
            
            WeekEndDay = .false.
            
        else
        
            WeekEndDay = .true.        
        
        endif
    
    end function WeekEndDay 
    
    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillTimeSeriesOperator(ObjTimeSeriesOperatorID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjTimeSeriesOperatorID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTimeSeriesOperatorID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mTimeSeriesOperator_,  Me%InstanceID)

            if (nUsers == 0) then
            
                call KillVariablesAndFiles
            
                !Deallocates Instance
                call DeallocateInstance ()

                ObjTimeSeriesOperatorID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillTimeSeriesOperator
        

    !------------------------------------------------------------------------
    



    subroutine KillVariablesAndFiles
    
    
        !Local--------------------------------------------------------------------------
        integer         :: STAT_CALL, i
        
        !Begin--------------------------------------------------------------------------
        
        !Kill Time Serie File Out
        call KillTimeSerie(Me%ObjTimeSerieOut, STAT = STAT_CALL)       
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - KillVariablesAndFiles - ERR10'
        
        !Kill Time Serie File Out Night
        call KillTimeSerie(Me%ObjTimeSerieOutNight, STAT = STAT_CALL)       
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - KillVariablesAndFiles - ERR15'        
        
        if (Me%FluxInON) then        
            do i=1, Me%NFluxIn
                !Kill Time Serie File In flux
                call KillTimeSerie(Me%ObjTimeSerieInFlux(i), STAT = STAT_CALL)       
                if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - KillVariablesAndFiles - ERR20'
            enddo
            deallocate(Me%ObjTimeSerieInFlux)
        endif

        if (Me%FluxOutON) then        
            do i=1, Me%NFluxOut
                !Kill Time Serie File Out Flux
                call KillTimeSerie(Me%ObjTimeSerieOutFlux(i), STAT = STAT_CALL)       
                if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - KillVariablesAndFiles - ERR30'
            enddo
            deallocate(Me%ObjTimeSerieOutFlux)
        endif
        
        if (Me%MassON) then
            do i=1, Me%NMass
                !Kill Time Serie File Mass
                call KillTimeSerie(Me%ObjTimeSerieMass   (i), STAT = STAT_CALL)       
                if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - KillVariablesAndFiles - ERR40'
            enddo
            deallocate(Me%ObjTimeSerieMass)
        endif
        
        call KillComputeTime(Me%ObjTime, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesOperator - KillVariablesAndFiles - ERR50'

    end subroutine KillVariablesAndFiles    
    
    !--------------------------------------------------------------------------    
    
   !--------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_TimeSeriesOperator), pointer          :: AuxObjTimeSeriesOperator
        type (T_TimeSeriesOperator), pointer          :: PreviousObjTimeSeriesOperator

        !Updates pointers
        if (Me%InstanceID == FirstObjTimeSeriesOperator%InstanceID) then
            FirstObjTimeSeriesOperator => FirstObjTimeSeriesOperator%Next
        else
            PreviousObjTimeSeriesOperator => FirstObjTimeSeriesOperator
            AuxObjTimeSeriesOperator      => FirstObjTimeSeriesOperator%Next
            do while (AuxObjTimeSeriesOperator%InstanceID /= Me%InstanceID)
                PreviousObjTimeSeriesOperator => AuxObjTimeSeriesOperator
                AuxObjTimeSeriesOperator      => AuxObjTimeSeriesOperator%Next
            enddo

            !Now update linked list
            PreviousObjTimeSeriesOperator%Next => AuxObjTimeSeriesOperator%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjTimeSeriesOperator_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTimeSeriesOperator_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjTimeSeriesOperator_ID > 0) then
            call LocateObjTimeSeriesOperator (ObjTimeSeriesOperator_ID)
            ready_ = VerifyReadLock (mTimeSeriesOperator_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjTimeSeriesOperator (ObjTimeSeriesOperatorID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTimeSeriesOperatorID

        !Local-----------------------------------------------------------------

        Me => FirstObjTimeSeriesOperator
        do while (associated (Me))
            if (Me%InstanceID == ObjTimeSeriesOperatorID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleTimeSeriesOperator - LocateObjTimeSeriesOperator - ERR01'

    end subroutine LocateObjTimeSeriesOperator

    !--------------------------------------------------------------------------

    end module ModuleTimeSeriesOperator