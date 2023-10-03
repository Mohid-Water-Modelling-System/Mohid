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
! DESCRIPTION   : Module to analyse (filter, interpolate, identify patterns, compare) Time Series
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

Module Module_TS_Synch

    use ModuleGlobalData
    use ModuleFunctions
    use ModuleEnterData
    use ModuleTime
    use ModuleTimeSerie
    
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: Construct_TS_Synch
    private ::      AllocateInstance

    !Selector
    public  :: Get_TS_Synch_Values
    public  :: Get_TS_Synch_ValueX
    public  :: Get_TS_SynchID
    public  :: Get_TS_SynchDT
    public  :: Get_TS_SyncName
    public  :: Get_TS_Synch_nValues
    public  :: UnGet_TS_Synch
                     
    
    !Modifier
    public  :: Modify_TS_Synch

    !Destructor
    public  :: Kill_TS_Synch                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObj_TS_Synch 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGet_TS_Synch1D_R8
    interface  UnGet_TS_Synch
        module procedure UnGet_TS_Synch1D_R8
    end interface  UnGet_TS_Synch


    !Parameter-----------------------------------------------------------------
    
    real(8), parameter  :: Pi_ = 3.1415926535897932384626433832795
    !Input / Output
    integer, parameter  :: FileOpen = 1, FileClose = 0
    
    !Interpolation methods
    integer, parameter  :: LinearTS_ = 1, BackwardTS_ = 2, CumulativeTS_ = 3
    
    !Types---------------------------------------------------------------------

    type T__TS_Synch
    
        integer                                                 :: InstanceID       =  0
                                                                
        real(8),                    dimension(:),   pointer     :: TimeSerie        => null()
	    real,                       dimension(:,:), pointer     :: DataMatrix       => null()
	    type (T_Time),              dimension(:)  , pointer     :: TimeTS           => null()
        type (T_Time),              dimension(:)  , pointer     :: TimeTSOutPut     => null()   
        character(len=Stringlength)                             :: Name             =  null_str     
        
	    character(Len=PathLength  )                             :: InterpolFile     =  null_str    
	    character(len=PathLength  )                             :: InputFileName    =  null_str    
	    character(len=PathLength  )                             :: OutputFileName   =  null_str         
        character(len=Stringlength)	                            :: TimeSerieName    = "Time Serie"        
        
	    type (T_Time)                                           :: InitialData
        type (T_Time)                                           :: BeginTime
        type (T_Time)                                           :: EndTime
        type (T_Time)                                           :: BeginTimeDefault
        type (T_Time)                                           :: EndTimeDefault        
        
	    real                                                    :: CoordX           = null_real
	    real                                                    :: CoordY           = null_real
	    real                                                    :: DT_Synch         = null_real
	    real                                                    :: DT_Synch_Default = null_real        
        real                                                    :: GapLimit         = null_real
        real                                                    :: FillValue        = null_real
        real                                                    :: FillValueDefault = null_real
	    

	    integer                                                 :: ObjEnterData     = 0
	    integer                                                 :: ObjTimeSerie     = 0
        
	    integer                                                 :: DataValues       = null_int    
        integer                                                 :: DataColumns      = null_int     
	    integer                                                 :: NValues          = null_int    
        integer                                                 :: DataColumn       = null_int

        integer                                                 :: AngleUnits       = null_int          
        integer                                                 :: AngleReferential = null_int    
        integer                                                 :: iInterpol        = null_int    
        integer                                                 :: InterpolInTime   = null_int
        
        logical                                                 :: AngleProp        = .false. 

        type(T__TS_Synch), pointer                     :: Next

    end type T__TS_Synch    

    !Global Variables
    type (T__TS_Synch), pointer                        :: FirstObj_TS_Synch
    type (T__TS_Synch), pointer                        :: Me    


    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Construct_TS_Synch(Obj_TS_SynchID, EnterDataID, ExtractType,             &
                                  BeginTime, EndTime, DT_Synch, FillValue, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: Obj_TS_SynchID 
        integer      ,          optional, intent(IN )   :: EnterDataID
        integer      ,          optional, intent(IN )   :: ExtractType      
        type (T_Time),          optional, intent(IN )   :: BeginTime
        type (T_Time),          optional, intent(IN )   :: EndTime
	    real         ,          optional, intent(IN )   :: DT_Synch
        real         ,          optional, intent(IN )   :: FillValue
        integer      ,          optional, intent(OUT)   :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, ExtractType_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(m_TS_Synch_)) then
            nullify (FirstObj_TS_Synch)
            call RegisterModule (m_TS_Synch_) 
        endif

        call Ready(Obj_TS_SynchID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            if (present(EnterDataID)) then
                Me%ObjEnterData      = AssociateInstance (mENTERDATA_,      EnterDataID     )            
            endif                    
                
            if (present(ExtractType))  then
                ExtractType_ =  ExtractType
            else 
                ExtractType_ = FromFile
            endif
                    
            if (present(BeginTime))  then
                Me%BeginTimeDefault =  BeginTime
            else 
               call null_time(Me%BeginTimeDefault)
            endif
            
            if (present(EndTime))  then
                Me%EndTimeDefault =  EndTime
            else 
               call null_time(Me%EndTimeDefault)
            endif

            if (present(DT_Synch))  then
                Me%DT_Synch_Default =  DT_Synch
            endif            
            
            if (present(FillValue))  then
                Me%FillValueDefault =  FillValue
            endif            
            
                
            call ReadKeywords(ExtractType_)    
                
            call ConstructRawTimeSerie
            
            call ConstructInterpolTS
            
            !Returns ID
            Obj_TS_SynchID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'Module_TS_Synch - Construct_TS_Synch - ERR40' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine Construct_TS_Synch
 
    !--------------------------------------------------------------------------

    subroutine ReadKeywords(ExtractType)

        !Arguments----------------------------------------------------------
        integer                 :: ExtractType
        !Local--------------------------------------------------------------
        integer                 :: flag, STAT_CALL
        
        !Begin--------------------------------------------------------------
        
        if (ExtractType == FromFile) then
        
            Me%ObjEnterData = 0
        
            call ConstructEnterData(Me%ObjEnterData, "TS_Synch.dat", STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR10'                
            
        endif

        call GetData(Me%DT_Synch,                                                       &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='DT_SYNCHRONISATION',                                &
                     ClientModule ='Module_TS_Synch',                                   &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR20'
        
        if (flag == 0) then
            Me%DT_Synch = Me%DT_Synch_Default
        endif   
        
        if (Me%DT_Synch <= 0) then
            stop 'Module_TS_Synch - ReadKeywords - ERR25'
        endif

        call GetData(Me%InputFileName,                                                  &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='INPUT_FILE',                                        &
                     ClientModule ='Module_TS_Synch',                                   &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR30'
        if (flag == 0)            stop 'Module_TS_Synch - ReadKeywords - ERR40'
        
        call GetData(Me%Name,                                                           &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='NAME',                                              &
                     ClientModule ='Module_TS_Synch',                                   &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR50'
        if (flag == 0)            stop 'Module_TS_Synch - ReadKeywords - ERR60'        
        
        call GetData(Me%DataColumn,                                                     &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='DATA_COLUMN',                                       &
                     ClientModule ='Module_TS_Synch',                                   &
                     default      = 2,                                                  &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR70'
        
        call GetData(Me%GapLimit,                                                       &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='TIME_GAP_LIMIT',                                    &
                     Default      = -FillValueReal,                                     &
                     ClientModule ='Module_TS_Synch',                                   &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR80'
        

        !Reads Begin Time
        call GetData(Me%BeginTime,                                                      &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='START',                                             &
                     ClientModule ='Module_TS_Synch',                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR90'
        
        if (flag == 0) then
            Me%BeginTime = Me%BeginTimeDefault
        endif
        
        !Reads End Time
        call GetData(Me%EndTime,                                                        &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='END',                                               &
                     ClientModule ='Module_TS_Synch',                                   &
                     STAT         = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR100'

        if (flag == 0) then
            Me%EndTime = Me%EndTimeDefault
        endif
        
        call GetData(Me%FillValue,                                                      &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='FILL_VALUE',                                        &
                     ClientModule ='Module_TS_Synch',                                   &
                     STAT         = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR110'
        
        if (flag == 0) then
            Me%FillValue = Me%FillValueDefault
        endif        
        
        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='OUTPUT_FILE',                                       &
                     ClientModule ='Module_TS_Synch',                                   &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR120'
        if (flag == 0)            stop 'Module_TS_Synch - ReadKeywords - ERR130'        
        
        call GetData(Me%CoordX,                                                         &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='COORD_X',                                           &
                     default      = FillValueReal,                                      &            
                     ClientModule ='Module_TS_Synch',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR140'
        
        call GetData(Me%CoordY,                                                         &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='COORD_Y',                                           &
                     default      = FillValueReal,                                      &            
                     ClientModule ='Module_TS_Synch',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR150'

        call GetData(Me%AngleProp,                                                      &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='ANGLE_PROPERTY',                                    &
                     default      = .false.,                                            &
                     ClientModule ='Module_TS_Synch',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR160'
        
        if (Me%AngleProp) then 
            call GetData(Me%AngleReferential,                                           &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = ExtractType,                                    &
                         keyword      ='ANGLE_REFERENTIAL',                             &
                         ClientModule ='Module_TS_Synch',                               &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR170'
            
            if (Me%AngleReferential /= NauticalWind_      .and.                         &
                Me%AngleReferential /= NauticalCurrent_   .and.                         &
                Me%AngleReferential /= CartesianDir_    ) then
                stop 'Module_TS_Synch - ReadKeywords - ERR180'
            endif            

            call GetData(Me%AngleUnits,                                                 &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = ExtractType,                                    &
                         keyword      ='ANGLE_UNITS',                                   &
                         ClientModule ='Module_TS_Synch',                               &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR190'            
            
            if (Me%AngleUnits /= Degree_      .and.                                     &
                Me%AngleUnits /= Radian_          ) then
                stop 'Module_TS_Synch - ReadKeywords - ERR200'
            endif                
        endif
        
        call GetData(Me%InterpolInTime,                                                 &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='INTERPOLATION_IN_TIME',                             &
                     Default      = LinearTS_,                                          &
                     ClientModule ='Module_TS_Synch',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR210'

        if (Me%InterpolInTime /= LinearTS_      .and.                                   &
            Me%InterpolInTime /= BackwardTS_    .and.                                   &
            Me%InterpolInTime /= CumulativeTS_) then
            stop 'Module_TS_Synch - ReadKeywords - ERR220'
        endif
            
        if (Me%InterpolInTime /= LinearTS_ .and. Me%AngleProp) then
            stop 'Module_TS_Synch - ReadKeywords - ERR230'
        endif
            
        if (ExtractType == FromFile) then
        
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ReadKeywords - ERR240'        
            
        endif
                
    
    end subroutine ReadKeywords
    
    !-------------------------------------------------------------------------
    !
    !character(len=PathLength) function AddString2FileName(Filename, AddString)
    !
    !    !Arguments------------------------------------------------------------    
    !    character(len=*)                :: Filename, AddString
    !    
    !    !Local----------------------------------------------------------------
    !    integer                         :: n, i, k
    !
    !    !Begin----------------------------------------------------------------    
    !    
    !    n = len_trim(Filename)
    !    
    !    k = FillValueInt
    !    
    !    do i=n,1,-1
    !        if (Filename(i:i) == "/" .or. Filename(i:i) == "\") then
    !            k = i
    !            exit
    !        endif                
    !    enddo
    !    
    !    if (k > FillValueInt) then
    !        AddString2FileName = Filename(1:k)//trim(AddString)//Filename(k+1:n)
    !    else
    !        AddString2FileName = trim(AddString)//trim(Filename)
    !    endif            
    !    
    !
    !end function AddString2FileName

!-------------------------------------------------------------------------    

    subroutine ConstructInterpolTS
    
        !Local----------------------------------------------------------------
        integer                         :: STAT_CALL

        !Begin----------------------------------------------------------------
        
        !Me%InterpolFile = "InterpolOut_"//Me%InputFileName 
        !Me%InterpolFile = AddString2FileName(Me%InputFileName, "InterpolOut_")
        
        !Open Output files
        call UnitsManager(Me%iInterpol, FileOpen, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Synch - ConstructInterpolTS - ERR10'
        
        open(unit = Me%iInterpol, file =trim(Me%OutputFileName), form = 'FORMATTED', status = 'UNKNOWN')

        
    end subroutine ConstructInterpolTS        
    
    !-------------------------------------------------------------------------        
    
    subroutine ConstructRawTimeSerie
                 
        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL, i
        type (T_Time)                                   :: auxBeginTime, auxEndTime

        !Begin----------------------------------------------------------------

        Me%ObjTimeSerie = 0

        call StartTimeSerieInput     (Me%ObjTimeSerie, TimeSerieDataFile = Me%InputFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "Module_TS_Synch - ConstructRawTimeSerie - ERR10"

        call GetTimeSerieInitialData (Me%ObjTimeSerie, Me%InitialData, STAT= STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "Module_TS_Synch - ConstructRawTimeSerie - ERR20"

        call GetTimeSerieDataValues  (Me%ObjTimeSerie, Me%DataValues,  STAT= STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "Module_TS_Synch - ConstructRawTimeSerie - ERR30"

        call GetTimeSerieDataColumns (Me%ObjTimeSerie, Me%DataColumns, STAT= STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop "Module_TS_Synch - ConstructRawTimeSerie - ERR40"
                       
        
        !Allocate matrixes
        allocate(Me%DataMatrix(Me%DataValues,Me%DataColumns))
        allocate(Me%TimeTS    (Me%DataValues))   
        
        call GetTimeSerieDataMatrix (Me%ObjTimeSerie, Me%DataMatrix, STAT= STAT_CALL) 
        
        auxBeginTime = Me%InitialData + Me%DataMatrix(1,1) 
        auxEndTime   = Me%InitialData + Me%DataMatrix(Me%DataValues,1) 
        
        if (Me%BeginTime < auxBeginTime) then
            write(*,*) "Start  synchronisation date can not be older than the sart date of the raw input time serie"
            stop "Module_TS_Synch - ConstructRawTimeSerie - ERR50"
        endif        
         
        if (Me%EndTime > auxEndTime) then
            write(*,*) "End synchronisation date can not be newer than the end date of the raw input time serie"
            stop "Module_TS_Synch - ConstructRawTimeSerie - ERR60"
        endif

        
        do i=1,Me%DataValues
            Me%TimeTS(i) = Me%InitialData + Me%DataMatrix(i,1)
            if (i>1) then
                if ( Me%TimeTS(i-1) >  Me%TimeTS(i)) then
                    write(*,*) 'Instant number', i, 'older than instant', i-1
                    write(*,*) 'Time series time must be order in ascending way'
                    stop "Module_TS_Synch - ConstructRawTimeSerie - ERR70"
                endif
            endif
        enddo
        
    end subroutine ConstructRawTimeSerie        
        
    !-------------------------------------------------------------------------    


    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T__TS_Synch), pointer                         :: NewObj_TS_Synch
        type (T__TS_Synch), pointer                         :: PreviousObj_TS_Synch


        !Allocates new instance
        allocate (NewObj_TS_Synch)
        nullify  (NewObj_TS_Synch%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObj_TS_Synch)) then
            FirstObj_TS_Synch         => NewObj_TS_Synch
            Me                    => NewObj_TS_Synch
        else
            PreviousObj_TS_Synch      => FirstObj_TS_Synch
            Me                    => FirstObj_TS_Synch%Next
            do while (associated(Me))
                PreviousObj_TS_Synch  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObj_TS_Synch
            PreviousObj_TS_Synch%Next => NewObj_TS_Synch
        endif

        Me%InstanceID = RegisterNewInstance (m_TS_Synch_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    
    !--------------------------------------------------------------------------
    subroutine Get_TS_Synch_Values (Obj_TS_SynchID, ValuesTS, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Obj_TS_SynchID
        real(8), dimension(:),  pointer                 :: ValuesTS
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_SynchID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(m_TS_Synch_, Me%InstanceID)

            ValuesTS => Me%TimeSerie

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine Get_TS_Synch_Values
    
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    subroutine Get_TS_Synch_ValueX (Obj_TS_SynchID, n, ValueX, FillValue, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Obj_TS_SynchID
        integer                                         :: n
        real(8)                                         :: ValueX
        logical                                         :: FillValue
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_SynchID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(m_TS_Synch_, Me%InstanceID)

            ValueX = Me%TimeSerie(n)
            
            if (ValueX == Me%FillValue) then
                FillValue = .true.
            else
                FillValue = .false.
            endif

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine Get_TS_Synch_ValueX
    
    !--------------------------------------------------------------------------
    
    
    subroutine Get_TS_SynchID (Obj_TS_SynchID, ID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Obj_TS_SynchID
        real                                            :: ID
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_SynchID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            ID = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine Get_TS_SynchID

    !--------------------------------------------------------------------------

    subroutine Get_TS_SynchDT (Obj_TS_SynchID, DT_Synch, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Obj_TS_SynchID
        real                                            :: DT_Synch
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_SynchID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            DT_Synch = Me%DT_Synch

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine Get_TS_SynchDT

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine Get_TS_SyncName (Obj_TS_SynchID, Name, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Obj_TS_SynchID
        character(len=*)                                :: Name
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_SynchID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Name = Me%Name

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine Get_TS_SyncName

    !--------------------------------------------------------------------------    
    
    subroutine Get_TS_Synch_nValues (Obj_TS_SynchID, nValues, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Obj_TS_SynchID
        integer                                         :: nValues
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_SynchID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            nValues = Me%nValues

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine Get_TS_Synch_nValues

    !--------------------------------------------------------------------------

    subroutine UnGet_TS_Synch1D_R8(Obj_TS_SynchID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Obj_TS_SynchID
        real(8), dimension(:), pointer                  :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_SynchID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(m_TS_Synch_, Me%InstanceID,  "UnGet_TS_Synch3D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGet_TS_Synch1D_R8

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Modify_TS_Synch(Obj_TS_SynchID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: Obj_TS_SynchID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_SynchID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            
            call ModifyInterpolTimeSeries
            
            
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine Modify_TS_Synch
    
    !--------------------------------------------------------------------------
    
    subroutine ModifyInterpolTimeSeries

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, STAT_CALL
        type (T_Time)                               :: Time1, Time2
        real                                        :: Value1, Value2, NewValue
        real                                        :: Year, Month, Day, hour, minute, second
        real                                        :: DT_Gap
        logical                                     :: TimeCycle
        
        !----------------------------------------------------------------------
    
        !Sample the input time serie with a regular time step DT_Synch
        
        !most be mutiple of 2
        Me%nValues = (Me%EndTime - Me%BeginTime) / Me%DT_Synch + 1
        
        allocate(Me%TimeSerie    (1:Me%nValues))

        
        allocate(Me%TimeTSOutPut (1:Me%nValues))

        Me%TimeTSOutPut(1) = Me%BeginTime
        j                  = 1
        
        do i=1, Me%nValues 
        
            if (i<Me%nValues) then
                Me%TimeTSOutPut(i+1) = Me%TimeTSOutPut(i) + Me%DT_Synch             
            endif

        enddo

        do i=1, Me%nValues                             

            call GetTimeSerieValue(Me%ObjTimeSerie, Me%TimeTSOutPut(i), Me%DataColumn, Time1, Value1,   &
                                    Time2, Value2, TimeCycle, STAT= STAT_CALL) 
            
            if (STAT_CALL /= SUCCESS_) stop "Module_TS_Synch_TS_Synch - ModifyInterpolTimeSeries - ERR10"
            
            if (TimeCycle) then
                Me%TimeSerie(i) =  Value1            
                cycle
            endif
            
            DT_Gap =  Time2 - Time1            
            
            if (DT_Gap > Me%GapLimit) then
                
                NewValue = Me%FillValue
                
            elseif (Value1 ==  Me%FillValue .or. Value2 == Me%FillValue) then
                
                write(*,*) ' Raw time series can not have FillValues or abnormal values'
                write(*,*) ' please filter the raw time serie with the NAME =', trim(Me%Name)
                stop "Module_TS_Synch_TS_Synch - ModifyInterpolTimeSeries - ERR20"
                
            else
            
                !Interpolates Value for current instant
                if      (Me%InterpolInTime == LinearTS_) then
                    call InterpolateValueInTimeAngleProof(Me%TimeTSOutPut(i), Time1, Value1, Time2,  &
                                                            Value2, NewValue)
                else if (Me%InterpolInTime == BackwardTS_) then
                    NewValue = Value1
                else if (Me%InterpolInTime == CumulativeTS_) then
                    
                    if (i==1) then
                        
                        NewValue = Me%FillValue 
                    else
                    
                        NewValue = GetTimeSerieCumulativeValue(TimeSerieID    = Me%ObjTimeSerie,      &                
                                                                StartTime      = Me%TimeTSOutPut(i-1), &
                                                                EndTime        = Me%TimeTSOutPut(i  ), &
                                                                DataColumn     = Me%DataColumn,        &
                                                                STAT           = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop "Module_TS_Synch_TS_Synch - ModifyInterpolTimeSeries - ERR30"                    
                        
                    endif
                    
                endif                                                
                
            endif
        
            Me%TimeSerie(i) =  NewValue

        enddo            
        
       
        write(Me%iInterpol,*) "NAME                    : ", trim(Me%TimeSerieName)
        write(Me%iInterpol,*) "LOCALIZATION_I          : -999999"
        write(Me%iInterpol,*) "LOCALIZATION_J          : -999999"
        write(Me%iInterpol,*) "LOCALIZATION_K          : -999999"
        
        call ExtractDate(Me%BeginTime, Year, Month, Day, hour, minute, second)
        
        write(Me%iInterpol,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
        write(Me%iInterpol,*) "TIME_UNITS              : SECONDS"
        write(Me%iInterpol,*) "COORD_X                 : ", Me%CoordX
        write(Me%iInterpol,*) "COORD_Y                 : ", Me%CoordY
        write(Me%iInterpol,*) "FILL_VALUE              : ", Me%FillValue
        
        write(Me%iInterpol,'(A21)') "Time InterpolatedData"
        
        write(Me%iInterpol,*) "<BeginTimeSerie>"         
        
        do i=1, Me%nValues        
            write(Me%iInterpol,*) Me%TimeTSOutPut(i)-Me%BeginTime, Me%TimeSerie(i)
        enddo
        
        write(Me%iInterpol,*) "<EndTimeSerie>"      
        

        !Open Output files
        call UnitsManager(Me%iInterpol, FileClose, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "Module_TS_Synch_TS_Synch - ModifyInterpolTimeSeries - ERR40"


        
        
    end subroutine ModifyInterpolTimeSeries
    
    !--------------------------------------------------------------------------    
    
    subroutine InterpolateValueInTimeAngleProof(ActualTime, Time1, Value1, Time2,  &
                                                Value2, NewValue)
    
    
        !Arguments-------------------------------------------------------------
        type(T_Time),      intent(IN)               :: ActualTime
        type(T_Time),      intent(IN)               :: Time1
        real,              intent(IN)               :: Value1
        type(T_Time),      intent(IN)               :: Time2
        real,              intent(IN)               :: Value2
        
        real,              intent(OUT)              :: NewValue

        !Local-----------------------------------------------------------------
        real                                        :: Aux1   , Aux2
        real                                        :: XValue1, XValue2
        real                                        :: YValue1, YValue2        
        real                                        :: XNewValue, YNewValue
        real                                        :: aux
    
        !Begin-----------------------------------------------------------------    
        
   
        if (Me%AngleProp) then
   
            if (Me%AngleUnits == Degree_) then

                call AngleFromFieldToGrid (AngleInReferential = Value1,                 &
                                           Referential        = Me%AngleReferential,    &
                                           GridAngle          = 0.,                     &
                                           AngleOutGrid       = Aux1)

                call AngleFromFieldToGrid (AngleInReferential = Value2,                 &
                                           Referential        = Me%AngleReferential,    &
                                           GridAngle          = 0.,                     &
                                           AngleOutGrid       = Aux2)            
            
                Aux1 = Aux1 * Pi_ / 180.
                Aux2 = Aux2 * Pi_ / 180.
                
            else
                Aux1 = Value1
                Aux2 = Value2
            endif                    
            
            !X component
            XValue1 = cos(Aux1)
            XValue2 = cos(Aux2)

            !Y component
            YValue1 = sin(Aux1)
            YValue2 = sin(Aux2)                
        
            !X component
            call InterpolateValueInTime(ActualTime, Time1, XValue1, Time2,  &
                                        XValue2, XNewValue)
            !Y component
            call InterpolateValueInTime(ActualTime, Time1, YValue1, Time2,  &
                                        YValue2, YNewValue)
                                        
            NewValue = atan2(YNewValue, XNewValue)
            
            if (Me%AngleUnits == Degree_) then

                aux =  NewValue * 180. / Pi_
            
                call AngleFromGridToField (AngleInGrid        = aux,               &
                                           Referential        = Me%AngleReferential,&
                                           GridAngle          = 0.,                 &
                                           AngleOutReferential= NewValue)             
            endif
            
        else                                        

        
            call InterpolateValueInTime(ActualTime, Time1, Value1, Time2,  &
                                        Value2, NewValue)                  
        endif
                                            
    end subroutine InterpolateValueInTimeAngleProof
    
    !--------------------------------------------------------------------------
    
    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine Kill_TS_Synch(Obj_TS_SynchID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: Obj_TS_SynchID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_SynchID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(m_TS_Synch_,  Me%InstanceID)

            if (nUsers == 0) then
            
                call KillVariablesAndFiles
            
                !Deallocates Instance
                call DeallocateInstance ()

                Obj_TS_SynchID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine Kill_TS_Synch
        

    !------------------------------------------------------------------------
    


    subroutine KillVariablesAndFiles
    
    
        !Local--------------------------------------------------------------------------
        !integer         :: STAT_CALL
        
        !Begin--------------------------------------------------------------------------
    
        deallocate(Me%DataMatrix    )
        deallocate(Me%TimeTS        )        
        deallocate(Me%TimeTSOutPut  )
        deallocate(Me%TimeSerie     )
    
    end subroutine KillVariablesAndFiles    
    
    !--------------------------------------------------------------------------    
    
   !--------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T__TS_Synch), pointer          :: AuxObj_TS_Synch
        type (T__TS_Synch), pointer          :: PreviousObj_TS_Synch

        !Updates pointers
        if (Me%InstanceID == FirstObj_TS_Synch%InstanceID) then
            FirstObj_TS_Synch => FirstObj_TS_Synch%Next
        else
            PreviousObj_TS_Synch => FirstObj_TS_Synch
            AuxObj_TS_Synch      => FirstObj_TS_Synch%Next
            do while (AuxObj_TS_Synch%InstanceID /= Me%InstanceID)
                PreviousObj_TS_Synch => AuxObj_TS_Synch
                AuxObj_TS_Synch      => AuxObj_TS_Synch%Next
            enddo

            !Now update linked list
            PreviousObj_TS_Synch%Next => AuxObj_TS_Synch%Next

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

    subroutine Ready (Obj_TS_Synch_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: Obj_TS_Synch_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (Obj_TS_Synch_ID > 0) then
            call LocateObj_TS_Synch (Obj_TS_Synch_ID)
            ready_ = VerifyReadLock (m_TS_Synch_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObj_TS_Synch (Obj_TS_SynchID)

        !Arguments-------------------------------------------------------------
        integer                                     :: Obj_TS_SynchID

        !Local-----------------------------------------------------------------

        Me => FirstObj_TS_Synch
        do while (associated (Me))
            if (Me%InstanceID == Obj_TS_SynchID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'Module_TS_Synch - LocateObj_TS_Synch - ERR01'

    end subroutine LocateObj_TS_Synch

    !--------------------------------------------------------------------------

    end module Module_TS_Synch