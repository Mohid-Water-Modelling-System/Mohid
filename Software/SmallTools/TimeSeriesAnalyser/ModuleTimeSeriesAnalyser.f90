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

Module ModuleTimeSeriesAnalyser

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
    public  :: ConstructTimeSeriesAnalyser
    private ::      AllocateInstance

    !Selector
    public  :: GetTimeSeriesAnalyserPointer
    public  :: GetTimeSeriesAnalyserInteger
    public  :: UnGetTimeSeriesAnalyser
                     
    
    !Modifier
    public  :: ModifyTimeSeriesAnalyser

    !Destructor
    public  :: KillTimeSeriesAnalyser                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjTimeSeriesAnalyser 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetTimeSeriesAnalyser3D_I
    private :: UnGetTimeSeriesAnalyser3D_R8
    interface  UnGetTimeSeriesAnalyser
        module procedure UnGetTimeSeriesAnalyser3D_I
        module procedure UnGetTimeSeriesAnalyser3D_R8
    end interface  UnGetTimeSeriesAnalyser


    !Parameter-----------------------------------------------------------------
    
    real(8), parameter  :: Pi_ = 3.1415926535897932384626433832795
    !Input / Output
    integer, parameter  :: FileOpen = 1, FileClose = 0
    
    !Interpolation methods
    integer, parameter  :: LinearTS_ = 1, BackwardTS_ = 2

    
    !Types---------------------------------------------------------------------
    type T_TimeSeriesAnalyser
    
        integer                                                 :: InstanceID
    
        complex (SPC),              dimension(:),   allocatable :: FFT
        real     (SP),              dimension(:),   allocatable :: TimeSerie, AuxTimeSerie, amplitude, phase, frequency, FlagTimeSerie 
        real     (SP),              dimension(:),   allocatable :: Tempo, AuxFreq
        
        
	    real,                       dimension(:,:), pointer     :: Hourly_Per
	    real,                       dimension(:,:), pointer     :: DataMatrix, OutData, DataCollection, Aux2D
	    integer,                    dimension(:)  , allocatable :: CollectionSize
	    character(len=Stringlength),dimension(:)  , pointer     :: PropertyList
	    character(len=PathLength  )                             :: TimeSerieDataFile, TimeSerieFilterIn
	    real,                       dimension(:)  , allocatable :: TSGap, TSFilter
	    type (T_Time),              dimension(:)  , allocatable :: TimeTS, TimeTSOutPut
	    integer      ,              dimension(:)  , allocatable :: FlagFilter
	    type (T_Time)                                           :: InitialData, Time1, Time2, BeginTime, EndTime, CurrentTime
        character(len=Stringlength)	                            :: TimeSerieName = "Time Serie"
	    real                                                    :: CoordX = 0.
	    real                                                    :: CoordY = 0.
	    real                                                    :: Value1, Value2, NewValue, DT_Analysis, Aux,ave,adev,sdev,var,skew
	    real                                                    :: curt, maxvalues, minvalues, t, hour
        real                                                    :: stdv, RMS, Error, Average    
        real    (SP)                                            :: ofac,hifac, prob

        integer (I4B)                                           :: isign, jmax

	    integer                                                 :: ObjEnterData          = 0

	    integer                                                 :: ObjTimeSerie          = 0
        integer                                                 :: ObjTimeSerieFilterIn  = 0
	    integer                                                 :: ObjTimeSerieInterpol  = 0

	    integer                                                 :: ObjTimeSerie2         = 0
	    integer                                                 :: ObjTimeSerieFilter2   = 0
	    integer                                                 :: ObjTimeSerieInterpol2 = 0
    	
	    integer                                                 :: ObjTimeSerieDiference = 0	
    	
	    integer                                                 :: STAT_CALL, DataValues, DataColumns 
	    integer                                                 :: NValues, DataColumn, nValuesAux
	    integer                                                 :: ObjTimeSerieOut, ObjTime, k, FlagColumn, nFTT
	    integer                                                 :: JulDayStart, JulDayEnd, JulDayWindow, JulDayCenter, JulDayTS
        integer                                                 :: FrequencyWindow
        integer                                                 :: FilterFlagColumn
        real                                                    :: FilterFlagLimit
        logical                                                 :: FilterFlagLimitAbove
        !Files Units
        integer                                                 :: iS, iP, iPE, iFilter, iGap, iInterpol
        
	    logical                                                 :: TimeCycle, OutWindow, SpectralAnalysis, PercentileAnalysis, PercentileEvolution
	    logical                                                 :: FilterTimeSerie
	    logical                                                 :: TimeSerieFilterInON
	    real                                                    :: FilterMinValue, FilterMaxValue, FilterMaxRateValue
        
        integer                                                 :: InterpolInTime = FillValueInt
        
        real                                                    :: GapLimit

        type(T_TimeSeriesAnalyser), pointer                     :: Next

    end type T_TimeSeriesAnalyser    

    !Global Variables
    type (T_TimeSeriesAnalyser), pointer                        :: FirstObjTimeSeriesAnalyser
    type (T_TimeSeriesAnalyser), pointer                        :: Me    


    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructTimeSeriesAnalyser(ObjTimeSeriesAnalyserID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjTimeSeriesAnalyserID 
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mTimeSeriesAnalyser_)) then
            nullify (FirstObjTimeSeriesAnalyser)
            call RegisterModule (mTimeSeriesAnalyser_) 
        endif

        call Ready(ObjTimeSeriesAnalyserID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            call ReadKeywords
            
            call ConstructRawTimeSerie
            
            call ConstructFilterTS
            
            call ConstructGapsTS
            
            call ConstructInterpolTS
            
            call ConstructPatternsTS
            

            !Returns ID
            ObjTimeSeriesAnalyserID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleTimeSeriesAnalyser - ConstructTimeSeriesAnalyser - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructTimeSeriesAnalyser
 
    !--------------------------------------------------------------------------

    subroutine ReadKeywords

        !Local--------------------------------------------------------------
        integer                 :: status, flag, STAT_CALL
        
        !Begin--------------------------------------------------------------

    
        Me%ObjEnterData = 0
        
        call ConstructEnterData(Me%ObjEnterData, "TimeSeriesAnalyser.dat", STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR10'
        
        call GetData(Me%DT_Analysis,                                                    &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='DT_ANALYSIS',                                       &
                     Default      = 60.,                                                &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR20'
        

        call GetData(Me%TimeSerieDataFile,                                              &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='INPUT_FILE',                                        &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR30'
        if (flag == 0)            stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR40'
        
        call GetData(Me%DataColumn,                                                     &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='DATA_COLUMN',                                       &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR50'
        if (flag == 0)            stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR60'
        
        call GetData(Me%FlagColumn,                                                     &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='FLAG_COLUMN',                                       &
                     default      = FillValueInt,                                       &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR70'
        !if (flag == 0)            stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR80'
        
       
        call GetData(Me%SpectralAnalysis,                                               &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='SPECTRAL_ANALYSIS',                                 &
                     Default      = .true.,                                             &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR90'

       
        call GetData(Me%PercentileAnalysis,                                             &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='PERCENTILE_ANALYSIS',                               &
                     Default      = .true.,                                             &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR100'
        
        call GetData(Me%JulDayWindow,                                                   &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='JULDAY_WINDOW',                                     &
                     Default      = 365,                                                &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR110'   
        
        if (Me%JulDayWindow < 1  ) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR120'   
        if (Me%JulDayWindow > 365) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR130'   
        
        call GetData(Me%FrequencyWindow,                                                &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='FREQUENCY_WINDOW',                                  &
                     !Number of hour of a week 
                     Default      = 168,                                                &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR140'
        
        if (Me%FrequencyWindow /= 24 .and. Me%FrequencyWindow /= 168) then
            stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR150'
        endif
        
        call GetData(Me%PercentileEvolution,                                            &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='PERCENTILE_EVOLUTION',                              &
                     Default      = .true.,                                             &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR160'

        call GetData(Me%InterpolInTime,                                                 &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='INTERPOLATION_IN_TIME',                             &
                     Default      = LinearTS_,                                          &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR170'

        if (Me%InterpolInTime /= LinearTS_ .and. Me%InterpolInTime /= BackwardTS_) then
            stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR180'
        endif


        call GetData(Me%GapLimit,                                                       &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='TIME_GAP_LIMIT',                                    &
                     Default      = -FillValueReal,                                     &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR190'
        
        call GetData(Me%FilterTimeSerie,                                                &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='FILTER_TIME_SERIE',                                 &
                     Default      = .true.,                                             &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR200'
        
        if (Me%FilterTimeSerie) then
        
            call GetData(Me%FilterMinValue,                                             &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='MIN_VALUE',                                     &
                         Default      = FillValueReal,                                  &
                         ClientModule ='ModuleTimeSeriesAnalyser',                      &
                         STAT         = STAT_CALL)        
            if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR210'

            call GetData(Me%FilterMaxValue,                                             &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='MAX_VALUE',                                     &
                         Default      = - FillValueReal,                                &
                         ClientModule ='ModuleTimeSeriesAnalyser',                      &
                         STAT         = STAT_CALL)        
            if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR220'


            call GetData(Me%FilterMaxRateValue,                                         &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='MAX_RATE',                                      &
                         Default      = - FillValueReal,                                &
                         ClientModule ='ModuleTimeSeriesAnalyser',                      &
                         STAT         = STAT_CALL)        
            if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR230'
            
            Me%TimeSerieFilterInON = .true.

            call GetData(Me%TimeSerieFilterIn,                                          &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='INPUT_FILTER_FILE',                             &
                         ClientModule ='ModuleTimeSeriesAnalyser',                      &
                         STAT         = STAT_CALL)        
            if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR240'

            if (flag == 0) Me%TimeSerieFilterInON = .false.
            
            if (Me%TimeSerieFilterInON) then
            
                call GetData(Me%FilterFlagColumn,                                       &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromFile,                                   &
                             keyword      ='FILTER_FLAG_COLUMN',                        &
                             default      = 2,                                          &
                             ClientModule ='ModuleTimeSeriesAnalyser',                  &
                             STAT         = STAT_CALL)        
                if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR250'

                call GetData(Me%FilterFlagLimit,                                        &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromFile,                                   &
                             keyword      ='FILTER_FLAG_LIMIT',                         &
                             default      = 0.5,                                        &
                             ClientModule ='ModuleTimeSeriesAnalyser',                  &
                             STAT         = STAT_CALL)        
                if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR260'

                call GetData(Me%FilterFlagLimitAbove,                                   &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromFile,                                   &
                             keyword      ='FILTER_FLAG_LIMIT_ABOVE',                   &
                             default      = .true.,                                     &
                             ClientModule ='ModuleTimeSeriesAnalyser',                  &
                             STAT         = STAT_CALL)        
                if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR270'

            endif                                
                    
        endif
        
        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR300'

    
    end subroutine ReadKeywords
    
    !-------------------------------------------------------------------------
 
    subroutine ConstructFilterTS
    
        !Local----------------------------------------------------------------
        character (Len=PathLength)      :: FilterFile
        integer                         :: STAT_CALL

        !Begin----------------------------------------------------------------
    
        if (Me%FilterTimeSerie) then
        
            FilterFile = "FilterOut_"//Me%TimeSerieDataFile 
            !Open Output files
            call UnitsManager(Me%iFilter, FileOpen, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ConstructFilterTS - ERR10'
            
            open(unit   = Me%iFilter, file =trim(FilterFile), form = 'FORMATTED', status = 'UNKNOWN')
            
            if (Me%TimeSerieFilterInON) then
            
                Me%ObjTimeSerieFilterIn = 0

                call StartTimeSerieInput(Me%ObjTimeSerieFilterIn,                       &
                                         TimeSerieDataFile = Me%TimeSerieFilterIn, STAT = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ConstructFilterTS - ERR20'
                
            endif
        endif

        
    end subroutine ConstructFilterTS        
    
    !-------------------------------------------------------------------------    

    subroutine ConstructGapsTS
    
        !Local----------------------------------------------------------------
        character (Len=PathLength)      :: GapsFile
        integer                         :: STAT_CALL

        !Begin----------------------------------------------------------------
    
        GapsFile = "GapsOut_"//Me%TimeSerieDataFile 
        !Open Output files
        call UnitsManager(Me%iGap, FileOpen, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ConstructGapsTS - ERR10'
        
        open(unit = Me%iGap, file =trim(GapsFile), form = 'FORMATTED', status = 'UNKNOWN')
        
        
    end subroutine ConstructGapsTS        
    
    !-------------------------------------------------------------------------        

!-------------------------------------------------------------------------    

    subroutine ConstructInterpolTS
    
        !Local----------------------------------------------------------------
        character (Len=PathLength)      :: InterpolFile
        integer                         :: STAT_CALL

        !Begin----------------------------------------------------------------
        
        InterpolFile = "InterpolOut_"//Me%TimeSerieDataFile 
        !Open Output files
        call UnitsManager(Me%iInterpol, FileOpen, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ConstructInterpolTS - ERR10'
        
        open(unit = Me%iInterpol, file =trim(InterpolFile), form = 'FORMATTED', status = 'UNKNOWN')

        
    end subroutine ConstructInterpolTS        
    
    !-------------------------------------------------------------------------        

    subroutine ConstructPatternsTS
    
        !Local----------------------------------------------------------------
        character (Len=PathLength)      :: SpectralAnalysisFile, PercentileFile, PerEvolutionFile
        integer                         :: STAT_CALL

        !Begin----------------------------------------------------------------
    
        if (Me%SpectralAnalysis) then
        
            SpectralAnalysisFile = "SpectralOut_"//Me%TimeSerieDataFile 
            !Open Output files
            call UnitsManager(Me%iS, FileOpen, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - Main - ERR10'
            
            open(unit   = Me%iS      , file =trim(SpectralAnalysisFile), form = 'FORMATTED',   &
                 status = 'UNKNOWN')
            
        endif

        if (Me%PercentileAnalysis) then
        
            PercentileFile       = "PercentileOut_"//Me%TimeSerieDataFile 
            
            !Open Output files
            call UnitsManager(Me%iP, FileOpen, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - Main - ERR20' 
            
            open(unit   = Me%iP      , file =trim(PercentileFile), form = 'FORMATTED',     &
                 status = 'UNKNOWN')
            
        endif

        if (Me%PercentileEvolution) then
        
            PerEvolutionFile       = "PerEvolutionOut_"//Me%TimeSerieDataFile 
            
            !Open Output files
            call UnitsManager(Me%iPE, FileOpen, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - Main - ERR20' 
            
            open(unit   = Me%iPE      , file =trim(PerEvolutionFile), form = 'FORMATTED',  &
                 status = 'UNKNOWN')
            
        endif
        
    end subroutine ConstructPatternsTS        
    
    !-------------------------------------------------------------------------    
    
    subroutine ConstructRawTimeSerie
                 
        !Local----------------------------------------------------------------
        integer      :: STAT_CALL, i

        !Begin----------------------------------------------------------------

        Me%ObjTimeSerie = 0

        call StartTimeSerieInput     (Me%ObjTimeSerie, TimeSerieDataFile = Me%TimeSerieDataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "ModuleTimeSeriesAnalyser - ConstructRawTimeSerie - ERR10"

        call GetTimeSerieInitialData (Me%ObjTimeSerie, Me%InitialData, STAT= STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "ModuleTimeSeriesAnalyser - ConstructRawTimeSerie - ERR20"

        call GetTimeSerieDataValues  (Me%ObjTimeSerie, Me%DataValues,  STAT= STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "ModuleTimeSeriesAnalyser - ConstructRawTimeSerie - ERR30"

        call GetTimeSerieDataColumns (Me%ObjTimeSerie, Me%DataColumns, STAT= STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop "ModuleTimeSeriesAnalyser - ConstructRawTimeSerie - ERR40"
                       
        
        !Allocate matrixes
        allocate(Me%DataMatrix(Me%DataValues,Me%DataColumns), Me%OutData(Me%DataValues,14))
        allocate(Me%TimeTS    (Me%DataValues))   
        allocate(Me%FlagFilter(Me%DataValues))   
        Me%FlagFilter(:) = 1
        allocate(Me%PropertyList(3)) 
        
        call GetTimeSerieDataMatrix (Me%ObjTimeSerie, Me%DataMatrix, STAT= STAT_CALL) 
        
        Me%BeginTime = Me%InitialData + Me%DataMatrix(1,1) 
        Me%EndTime   = Me%InitialData + Me%DataMatrix(Me%DataValues,1) 
        
        do i=1,Me%DataValues
            Me%TimeTS(i) = Me%InitialData + Me%DataMatrix(i,1)
        enddo
        
        
        
    end subroutine ConstructRawTimeSerie        
        
    !-----------------------------------------------------------------------------------        
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_TimeSeriesAnalyser), pointer                         :: NewObjTimeSeriesAnalyser
        type (T_TimeSeriesAnalyser), pointer                         :: PreviousObjTimeSeriesAnalyser


        !Allocates new instance
        allocate (NewObjTimeSeriesAnalyser)
        nullify  (NewObjTimeSeriesAnalyser%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjTimeSeriesAnalyser)) then
            FirstObjTimeSeriesAnalyser         => NewObjTimeSeriesAnalyser
            Me                    => NewObjTimeSeriesAnalyser
        else
            PreviousObjTimeSeriesAnalyser      => FirstObjTimeSeriesAnalyser
            Me                    => FirstObjTimeSeriesAnalyser%Next
            do while (associated(Me))
                PreviousObjTimeSeriesAnalyser  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjTimeSeriesAnalyser
            PreviousObjTimeSeriesAnalyser%Next => NewObjTimeSeriesAnalyser
        endif

        Me%InstanceID = RegisterNewInstance (mTimeSeriesAnalyser_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    
    !--------------------------------------------------------------------------
    subroutine GetTimeSeriesAnalyserPointer (ObjTimeSeriesAnalyserID, Matrix, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjTimeSeriesAnalyserID
        real(8), dimension(:, :, :),  pointer           :: Matrix
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTimeSeriesAnalyserID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mTimeSeriesAnalyser_, Me%InstanceID)

            !Matrix => Me%Matrix

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetTimeSeriesAnalyserPointer
    
    !--------------------------------------------------------------------------
    
    subroutine GetTimeSeriesAnalyserInteger (ObjTimeSeriesAnalyserID, Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjTimeSeriesAnalyserID
        real                                            :: Int
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTimeSeriesAnalyserID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Int = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetTimeSeriesAnalyserInteger

    !--------------------------------------------------------------------------

    subroutine UnGetTimeSeriesAnalyser3D_I(ObjTimeSeriesAnalyserID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjTimeSeriesAnalyserID
        integer, dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTimeSeriesAnalyserID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mTimeSeriesAnalyser_, Me%InstanceID, "UnGetTimeSeriesAnalyser3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetTimeSeriesAnalyser3D_I

    !--------------------------------------------------------------------------

    subroutine UnGetTimeSeriesAnalyser3D_R8(ObjTimeSeriesAnalyserID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjTimeSeriesAnalyserID
        real(8), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTimeSeriesAnalyserID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mTimeSeriesAnalyser_, Me%InstanceID,  "UnGetTimeSeriesAnalyser3D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetTimeSeriesAnalyser3D_R8

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyTimeSeriesAnalyser(ObjTimeSeriesAnalyserID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTimeSeriesAnalyserID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTimeSeriesAnalyserID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            if (Me%FilterTimeSerie) then
                call ModifyFilterTimeSeries
            endif

            call ModifyInterpolTimeSeries
            
            call ModifyTimeSeriesPatterns
            

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyTimeSeriesAnalyser
    
    !--------------------------------------------------------------------------
    

    
    !--------------------------------------------------------------------------
    
    subroutine ModifyFilterTimeSeries

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,         pointer, dimension(:,:)       :: FilterDataMatrix
        type (T_Time)                               :: FilterInitialData, Time1, Time2
        integer                                     :: FilterDataValues, FilterDataColumns
        integer                                     :: i, j, STAT_CALL
        real                                        :: Rate, Value1, Value2
        logical                                     :: search
        
        !----------------------------------------------------------------------
    
        !Sample the input time serie with a regular time step DT_Analysis
d1:     do i=1, Me%DataValues 
            ! Initial value equal to the flag quality of the raw data (if exist)
            if (Me%FlagColumn == FillValueInt) then
        
                Me%FlagFilter(i) = 1            
            
            else                
                Me%FlagFilter(i) = Me%DataMatrix(i,Me%FlagColumn)
            endif
            
            !Filter minimum values
            if (Me%DataMatrix(i,Me%DataColumn) < Me%FilterMinValue) then
                Me%FlagFilter(i) = 0
            endif
            
            !Filter maximum values
            if (Me%DataMatrix(i,Me%DataColumn) > Me%FilterMaxValue) then
                Me%FlagFilter(i) = 0
            endif

            !Filter maximum rate changes
            if (i < Me%DataValues) then
                Rate = abs(Me%DataMatrix(i+1,Me%DataColumn)-Me%DataMatrix(i,Me%DataColumn))/(Me%TimeTS(i+1)-Me%TimeTS(i))
                if (Rate > Me%FilterMaxRateValue) then
                    Me%FlagFilter(i  ) = 0
                    Me%FlagFilter(i+1) = 0
                endif
            endif

        enddo d1
        
i1:     if (Me%TimeSerieFilterInON) then
        
                    
            call GetTimeSerieInitialData (Me%ObjTimeSerieFilterIn, FilterInitialData, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ModifyFilterTimeSeries - ERR30'
            
            call GetTimeSerieDataValues  (Me%ObjTimeSerieFilterIn, FilterDataValues,  STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ModifyFilterTimeSeries - ERR40'

            call GetTimeSerieDataColumns (Me%ObjTimeSerieFilterIn, FilterDataColumns, STAT= STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ModifyFilterTimeSeries - ERR50'
            
            if (Me%FilterFlagColumn > FilterDataColumns) then
                stop 'ModuleTimeSeriesAnalyser - ModifyFilterTimeSeries - ERR60'
            endif

            !Allocate matrixes
            allocate(FilterDataMatrix(FilterDataValues,FilterDataColumns))
    
            call GetTimeSerieDataMatrix (Me%ObjTimeSerieFilterIn, FilterDataMatrix, STAT= STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ModifyFilterTimeSeries - ERR70'
            
            j = 1

d2:         do i=1, Me%DataValues             

                search = .true.
                do while (search)
                
                    Time1 = FilterInitialData + FilterDataMatrix(j,1) 
                    if (j+1 > FilterDataValues) then
                        stop 'ModuleTimeSeriesAnalyser - ModifyFilterTimeSeries - ERR80'
                    endif
                    Time2 = FilterInitialData + FilterDataMatrix(j+1,1)                 

                    if (Time1 <= Me%TimeTS(i) .and. Time2 >= Me%TimeTS(i)) then
                        search = .false. 
                    else
                        j = j + 1
                    endif                            
                enddo
                
                Value1 = FilterDataMatrix(j,   Me%FilterFlagColumn)
                Value2 = FilterDataMatrix(j+1, Me%FilterFlagColumn)       
                
                if (Me%FilterFlagLimitAbove) then
                
                    if (Value1 > Me%FilterFlagLimit .or. Value2 > Me%FilterFlagLimit) then
                        Me%FlagFilter(i  ) = 0                            
                    endif
                    
                else
                
                    if (Value1 < Me%FilterFlagLimit .or. Value2 < Me%FilterFlagLimit) then
                        Me%FlagFilter(i  ) = 0                            
                    endif
                
                endif            
                
            enddo d2
            
        endif i1

    end subroutine ModifyFilterTimeSeries
    
    !--------------------------------------------------------------------------
 
    subroutine ModifyInterpolTimeSeries

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, STAT_CALL
        type (T_Time)                               :: Time1, Time2
        real                                        :: Value1, Value2, NewValue
        real                                        :: Year, Month, Day, hour, minute, second
        logical                                     :: search
        
        !----------------------------------------------------------------------
    
        !Sample the input time serie with a regular time step DT_Analysis
        
        !most be mutiple of 2
        Me%nValues = (Me%EndTime - Me%BeginTime) / Me%DT_Analysis + 1
        
        allocate(Me%TimeSerie    (1:Me%nValues))
        allocate(Me%FlagTimeSerie(1:Me%nValues))
        allocate(Me%amplitude    (1:Me%nValues/2), Me%phase(1:Me%nValues/2), Me%frequency(1:Me%nValues/2), Me%AuxFreq(1:Me%nValues/2))
        allocate(Me%TimeTSOutPut (1:Me%nValues))
        allocate(Me%TSGap        (1:Me%nValues))

        Me%TimeTSOutPut(1) = Me%BeginTime
        j                  = 1
        
        do i=1, Me%nValues 
        
            !0 - bad value, 1 - good value
            if (Me%FilterTimeSerie) then
                search = .true.
                do while (search)
                    if (Me%TimeTS(j)   <= Me%TimeTSOutPut(i) .and.                     &
                        Me%TimeTS(j+1) >= Me%TimeTSOutPut(i)) then
                        search = .false. 
                    else
                        j = j + 1
                    endif                            
                enddo
                
                Value1 = Me%FlagFilter(j  )
                Value2 = Me%FlagFilter(j+1)                    
            else
                if (Me%FlagColumn == FillValueInt) then
            
                    Value1 = 1            
                    Value2 = 1
                
                else
                    call GetTimeSerieValue(Me%ObjTimeSerie, Me%TimeTSOutPut(i), Me%FlagColumn, &
                                           Time1, Value1, Time2, Value2, Me%TimeCycle, STAT= STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ModifyInterpolTimeSeries - ERR10'
                endif                
            endif                                                            
            
            if (Value1 < 1 .or. Value2 < 1) then
                Me%FlagTimeSerie(i) = 0
            else
                Me%FlagTimeSerie(i) = 1
            endif
                
            if (i<Me%nValues) then
                Me%TimeTSOutPut(i+1) = Me%TimeTSOutPut(i) + Me%DT_Analysis             
            endif

        enddo

        do i=1, Me%nValues                                

            if (Me%FlagTimeSerie(i) == 1)  then
                call GetTimeSerieValue(Me%ObjTimeSerie, Me%TimeTSOutPut(i), Me%DataColumn, Time1, Value1,   &
                                       Time2, Value2, Me%TimeCycle, STAT= STAT_CALL) 
            
                                   
                !Interpolates Value for current instant
                if      (Me%InterpolInTime == LinearTS_) then
                    call InterpolateValueInTime(Me%TimeTSOutPut(i), Time1, Value1, Time2,  &
                                                Value2, NewValue)
                else if (Me%InterpolInTime == BackwardTS_) then
                    NewValue = Value1
                endif                                                
        
                Me%TimeSerie(i) =  NewValue
                Me%TSGap    (i) =  Time2 - Time1
            else
                Me%TimeSerie(i) =  FillValueReal
            endif
            
        enddo            
        
        
        do i=1, Me%nValues        
                
            if (Me%FlagTimeSerie(i) == 0)  then
                !lower limit
                j = i-1
                Time1  =  Me%TimeTSOutPut(1)
                Value1 =  FillValueReal
                
                do while (j>=1) 
                    if (Me%FlagTimeSerie(j) == 1) then
                        
                        Time1  =  Me%TimeTSOutPut(j)
                        Value1 =  Me%TimeSerie   (j)
                        exit
                    else
                        j = j - 1
                    endif
                enddo
                
                !upper limit
                j = i+1
                Time2  =  Me%TimeTSOutPut(Me%nValues)
                Value2 =  FillValueReal
                
                do while (j<=Me%nValues) 
                    if (Me%FlagTimeSerie(j) == 1) then
                        
                        Time2  =  Me%TimeTSOutPut(j)
                        Value2 =  Me%TimeSerie   (j)
                        exit
                    else
                        j = j + 1
                    endif
                enddo
                
                if (Value1 > FillValueReal .and. Value2 > FillValueReal) then
                
                    !Interpolates Value for current instant
                    if      (Me%InterpolInTime == LinearTS_) then
                        call InterpolateValueInTime(Me%TimeTSOutPut(i), Time1, Value1,     &
                                                    Time2, Value2, NewValue)     
                    else if (Me%InterpolInTime == BackwardTS_) then
                        NewValue = Value1
                    endif                                                                

                    Me%TSGap    (i) =  Time2 - Time1                                                
                                                           
                else if (Value1 > FillValueReal) then
                    NewValue = Value1
                    Me%TSGap(i) = Me%TimeTSOutPut(i) - Time1
                else if (Value2 > FillValueReal) then
                    NewValue = Value2
                    Me%TSGap(i) = Time2 - Me%TimeTSOutPut(i)
                else
                    stop "TimeSeriesAnalyser - ConstructTimeSeriesAnalyser - ERR30"
                endif                                      
                
                Me%TimeSerie(i) =  NewValue          

            endif            

        enddo
        
        write(Me%iInterpol,*) "NAME                    : ", Me%TimeSerieName
        write(Me%iInterpol,*) "LOCALIZATION_I          : -999999"
        write(Me%iInterpol,*) "LOCALIZATION_J          : -999999"
        write(Me%iInterpol,*) "LOCALIZATION_K          : -999999"
        
        call ExtractDate(Me%BeginTime, Year, Month, Day, hour, minute, second)
        
        write(Me%iInterpol,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
        write(Me%iInterpol,*) "TIME_UNITS              : SECONDS"
        write(Me%iInterpol,*) "MODEL_DOMAIN            : Sintra"
        write(Me%iInterpol,*) "COORD_X    : ", Me%CoordX
        write(Me%iInterpol,*) "COORD_Y    : ", Me%CoordY
        
        write(Me%iInterpol,'(A82)') "Time InterpolatedData"
        
        write(Me%iInterpol,*) "<BeginTimeSerie>"         
        
        do i=1, Me%nValues        
            write(Me%iInterpol,*) Me%TimeTSOutPut(i)-Me%BeginTime, Me%TimeSerie(i)
        enddo
        
        write(Me%iInterpol,*) "<EndTimeSerie>"          
        

        do i=1, Me%nValues        
            if (Me%TSGap(i) > Me%GapLimit) then
                Me%FlagTimeSerie(i) = 0
            endif
        enddo
        
        
        if (Me%FilterTimeSerie) then
        
            write(Me%iFilter,*) "NAME                    : ", Me%TimeSerieName
            write(Me%iFilter,*) "LOCALIZATION_I          : -999999"
            write(Me%iFilter,*) "LOCALIZATION_J          : -999999"
            write(Me%iFilter,*) "LOCALIZATION_K          : -999999"
            call ExtractDate(Me%BeginTime, Year, Month, Day, hour, minute, second)
            write(Me%iFilter,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
            write(Me%iFilter,*) "TIME_UNITS              : SECONDS"
            write(Me%iFilter,*) "MODEL_DOMAIN            : Sintra"
            write(Me%iFilter,*) "COORD_X    : ", Me%CoordX
            write(Me%iFilter,*) "COORD_Y    : ", Me%CoordY
            
            write(Me%iFilter,'(A82)') "Time FilterData"
            
            write(Me%iFilter,*) "<BeginTimeSerie>" 
               
            do i=1, Me%nValues        
                if (Me%FlagTimeSerie(i) == 1) then
                    write(Me%iFilter,*) Me%TimeTSOutPut(i)-Me%BeginTime, Me%TimeSerie(i)
                endif
            enddo
        
            write(Me%iFilter,*) "<EndTimeSerie>"  

        endif
        
        
    end subroutine ModifyInterpolTimeSeries
    
    !--------------------------------------------------------------------------
       
    subroutine ModifyTimeSeriesPatterns
    
        !Local--------------------------------------------------------------
        
        !Begin--------------------------------------------------------------
    
    
        call moment(Me%TimeSerie,Me%ave,Me%adev,Me%sdev,Me%var,Me%skew,Me%curt)    
        
        !Determinação dos 
        if  (Me%SpectralAnalysis) then
            call ComputeSpectralAnalysis
        endif

        if  (Me%PercentileAnalysis) then
            call ComputePercentileAnalysis
        endif        

       !Testing padron based in frequency analsys
        if  (Me%PercentileEvolution) then
            call ComputePercentileEvolution
        endif        
        
    end subroutine ModifyTimeSeriesPatterns    

    !--------------------------------------------------------------------------    
       
    subroutine  ComputeSpectralAnalysis        

        !Local--------------------------------------------------------------
        real                    :: Aux
        integer                 :: j
        !Begin--------------------------------------------------------------
    
        !http://en.wikipedia.org/wiki/Logarithm
        !Aux = log2(nValuesAux)
        Aux = log(real(Me%nValues))/log(2.)
        
        Me%nFTT = 2**int(Aux)
        
        allocate(Me%AuxTimeSerie(1:Me%nFTT), Me%FFT(1:Me%nFTT/2) )    
        
        Me%AuxTimeSerie(1:Me%nFTT) = Me%TimeSerie(1:Me%nFTT)
        
        Me%AuxTimeSerie(1:Me%nFTT) = Me%TimeSerie(1:Me%nFTT)

        call realft_sp(data=Me%AuxTimeSerie,isign=1,zdata=Me%FFT)

        write(Me%iS,'(A95)') "frequency[s-1], Period[days], amplitude[m], amplitude**2, amplitude**2/frequency phase*180./Pi_"
        
        do j=1,Me%nFTT/2
            Me%amplitude(j) = cabs(Me%FFT(j)) / real(Me%nFTT)
            Me%phase(j)     = atan2(imag(Me%FFT(j)),real(Me%FFT(j)))
            Me%frequency(j) = real(j) / real(Me%nValues) / Me%DT_Analysis
            
            write(Me%iS,'(6e18.5)') Me%frequency(j), 1/Me%frequency(j)/3600/24, Me%amplitude(j), Me%amplitude(j)**2, &
                                 Me%amplitude(j)**2/Me%frequency(j), Me%phase(j)*180./Pi_
        enddo   
        
    end subroutine ComputeSpectralAnalysis     
    
    !--------------------------------------------------------------------------    
    
    subroutine ComputePercentileAnalysis             

        !Local--------------------------------------------------------------
        real                        :: Year, Month, Day, hour, minute, second
        integer                     :: j, NPerHour
        !Begin--------------------------------------------------------------

        
        NPerHour = 2*int(Me%nValues/Me%FrequencyWindow) 
        
        allocate(Me%Hourly_Per    (0:Me%FrequencyWindow-1,1:11)    )
        allocate(Me%DataCollection(0:Me%FrequencyWindow-1,NPerHour))
        allocate(Me%CollectionSize(0:Me%FrequencyWindow-1)         )        
        
        call ComputePercentile(Me%EndTime, Me%BeginTime, Me%Hourly_Per) 
        
        write(Me%iP,*) "NAME                    :",  Me%TimeSerieName
        write(Me%iP,*) "LOCALIZATION_I          : -999999"
        write(Me%iP,*) "LOCALIZATION_J          : -999999"
        write(Me%iP,*) "LOCALIZATION_K          : -999999"
        call ExtractDate(Me%BeginTime, Year, Month, Day, hour, minute, second)
        write(Me%iP,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
        write(Me%iP,*) "TIME_UNITS              : HOURS"
        write(Me%iP,*) "COORD_X    : ", Me%CoordX
        write(Me%iP,*) "COORD_Y    : ", Me%CoordY
        
        write(Me%iP,'(A82)') "Time Per_0 Per_10 Per_20 Per_30 Per_40 Per_50 Per_60 Per_70 Per_80 Per_90  Per_100"
        
        write(Me%iP,*) "<BeginTimeSerie>" 
        

        do j=0,Me%FrequencyWindow-1
            write(Me%iP,'(I4,11(f14.6," "))') j, Me%Hourly_Per(j,1:11)
        enddo
        
        write(Me%iP,*) "<EndTimeSerie>"  
        
        deallocate(Me%Hourly_Per    )
        deallocate(Me%DataCollection)
        deallocate(Me%CollectionSize)        

    end subroutine ComputePercentileAnalysis      

    !--------------------------------------------------------------------------    
    
    subroutine ComputePercentileEvolution               

        !Local--------------------------------------------------------------
        real                    :: Year, Month, Day, hour, minute, second, aux, dti
        integer                 :: i, j, NPerHour, k, istart, ifinal, j1, j2
        !Begin--------------------------------------------------------------

        
        NPerHour = 2*int(Me%nValues/Me%FrequencyWindow) 
        
        allocate(Me%DataCollection(0:Me%FrequencyWindow-1,NPerHour))
        allocate(Me%CollectionSize(0:Me%FrequencyWindow-1)         )
        allocate(Me%Hourly_Per    (0:Me%FrequencyWindow-1,1:11)    )
        
        Me%CollectionSize(:) = 0
        
        Me%CurrentTime = Me%BeginTime + 3600.*Me%FrequencyWindow
        
        i = 1
        
        call ComputePercentile(Me%EndTime, Me%BeginTime, Me%Hourly_Per) 

        do while (Me%TimeTS(i)<=Me%EndTime) 
            
            !call ExtractDate(TimeTS(i), hour = hour)
            
            !j = int(hour)

            aux = (Me%TimeTS(i) - Me%BeginTime)/real(Me%FrequencyWindow)/3600.
            aux = (aux-int(aux)) * Me%FrequencyWindow
            j   = int(aux)
            dti = aux-real(j)
            if (dti > 0.5) then
                j1 = j
                j2 = j+1
                if (j1==Me%FrequencyWindow-1) j2=0
                dti = dti-0.5
            else
                j1 = j
                j2 = j-1
                if (j1==0) j2=Me%FrequencyWindow-1
                dti = 0.5 - dti
            endif            
            
            
            Me%OutData(i,1) = Me%DataMatrix(i,1)
            do k = 1, 11
                Me%OutData(i,k+1) = Me%Hourly_Per(j2,k)*dti +  Me%Hourly_Per(j1,k)*(1.-dti)
            enddo
            Me%OutData(i,13) = Me%DataMatrix(i,Me%DataColumn)
            Me%OutData(i,14) = Me%FlagFilter(i)
        
            i = i + 1
            if (i > Me%DataValues) exit
        enddo                
        
        
        
        !do while (Me%TimeTS(j)<= Me%BeginTime + 3600.*Me%FrequencyWindow)     
        
        !    j = j  + 1
        
        !enddo
        
        istart = 1
        ifinal = i - 1               
        
        
        write(Me%iPE,*) "NAME                    : Q1m3h"
        write(Me%iPE,*) "LOCALIZATION_I          : -999999"
        write(Me%iPE,*) "LOCALIZATION_J          : -999999"
        write(Me%iPE,*) "LOCALIZATION_K          : -999999"
        call ExtractDate(Me%InitialData, Year, Month, Day, hour, minute, second)        
        write(Me%iPE,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
        write(Me%iPE,*) "TIME_UNITS              : SECONDS"
        write(Me%iPE,*) "MODEL_DOMAIN            : Sintra"
        write(Me%iPE,*) "COORD_X    :  -4.140000    "
        write(Me%iPE,*) "COORD_Y    :   36.71000     "
        
        write(Me%iPE,'(A256)') "Time Per_0 Per_10 Per_20 Per_30 Per_40 Per_50 Per_60 Per_70 Per_80 Per_90  Per_100 Raw_Data stdv RMS Error[%]"
        
        write(Me%iPE,*) "<BeginTimeSerie>" 
        

        do j=istart,ifinal
        
            if (Me%FlagFilter(j) == 0) cycle
            !Compute Per_50 error
            Me%RMS  = 0
            Me%Average=sum(Me%OutData(1:j,13))/real(j)
            Me%stdv = -Me%Average**2        
            do k=1, j
                if (j==1) then
                    Me%stdv = (Me%OutData(1,13))**2
                else
                    Me%stdv = Me%stdv + Me%OutData(k,13)**2/real(j)
                endif
            enddo    
            do k=istart, j                
                Me%RMS = Me%RMS + (Me%OutData(k,13)-Me%OutData(k,7))**2.
            enddo
            
            Me%RMS  = sqrt(Me%RMS / real(j-istart+1))
            Me%stdv = sqrt(Me%stdv)
            !Error in percentage normalized using the standard deviation multiply by 4. In a normal distribution 
            !corresponds to ~95% of the sample 
            if (Me%stdv >0) then
                Me%Error = Me%RMS / (4*Me%stdv) * 100.
            else
                !standard deviation equal zero. Error is assume equal to 100%
                Me%Error = 100
            endif
            write(Me%iPE,'(f12.2,16(e14.6," "))') Me%OutData(j,1:13), Me%stdv, Me%RMS, Me%Error
        enddo
        
        write(Me%iPE,*) "<EndTimeSerie>"  
        
        deallocate(Me%Hourly_Per    )
        deallocate(Me%DataCollection)
        deallocate(Me%CollectionSize)        

    end subroutine ComputePercentileEvolution  
    

    !------------------------------------------------------------------

    subroutine ComputePercentile(NowTime, ReferenceTime, Hourly_Per_Aux)               

        !Arguments----------------------------------------------------------
        type(T_Time)                            :: NowTime, ReferenceTime
        real,       dimension(:,:), pointer     :: Hourly_Per_Aux
        
        !Local--------------------------------------------------------------
        real   (SP),dimension(:),   allocatable :: SortArray
        real                                    :: aux
        logical                                 :: InvertOrder, OutWindow
        integer                                 :: ii, i, j, JulDayCenter, JulDayStart, JulDayTS, JulDayEnd
        !Begin--------------------------------------------------------------
    
        !Sample data for the percentile analysis 

        call JulianDay(NowTime, JulDayCenter)
        
        JulDayStart = JulDayCenter - int(Me%JulDayWindow/2.)
        if (JulDayStart < 1) JulDayStart = 365 + JulDayStart 
        JulDayEnd   = JulDayCenter + int(Me%JulDayWindow/2.)
        if (JulDayEnd > 365) JulDayEnd = JulDayEnd - 365
        
        if (JulDayEnd < JulDayStart) then
            OutWindow = .true.
        else
            OutWindow = .false.
        endif
        
        Me%CollectionSize(:) = 0
        
        i=1
        
        do while (Me%TimeTSOutPut(i)<=NowTime) 
        
            if (Me%FlagTimeSerie(i) == 1) then
            
                call JulianDay(Me%TimeTSOutPut(i), JulDayTS)
                
                if (OutWindow) then
                    if (JulDayTS > JulDayEnd .and. JulDayTS < JulDayStart) then
                        i = i + 1
                        cycle
                    endif
                else
                    if (JulDayTS > JulDayEnd .or.  JulDayTS < JulDayStart) then
                        i = i + 1
                        cycle
                    endif
                endif
                
                !call ExtractDate(TimeTS(i), hour = hour)
                
                !j = int(hour)
                aux = (Me%TimeTSOutPut(i)-ReferenceTime)/real(Me%FrequencyWindow) / 3600.
                
                if (aux < 0) then
                    aux = abs(aux)
                    
                    InvertOrder = .true.
                else
                    InvertOrder = .false.                                
                endif
                
                j = int((aux-int(aux)) * Me%FrequencyWindow)
                
                Me%CollectionSize(j) = Me%CollectionSize(j) + 1
                
                Me%DataCollection(j,Me%CollectionSize(j)) = Me%TimeSerie(i)
            endif
            
            i = i + 1
            
            if (i > Me%nValues) exit
            
        enddo      
        
        !Compute percentile matrix
        
        Hourly_Per_Aux(:,:) = FillValueReal
        
        do i=0, Me%FrequencyWindow-1
        
            if (Me%CollectionSize(i) > 0) then

                allocate(SortArray(Me%CollectionSize(i)))
                
                SortArray(:) = Me%DataCollection(i,1:Me%CollectionSize(i))
                
                if (InvertOrder) then
                    ii = i
                else
                    ii = Me%FrequencyWindow + i - 1                    
                endif
                    
                call sort(SortArray)     
                do j=1,11
                    Hourly_Per_Aux(i,j) = FreqAnalysis(SortArray,Me%CollectionSize(i),real(j-1)/10)    
                enddo
               
                deallocate(SortArray)
            endif
        enddo
            
    end subroutine ComputePercentile
    
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
    
   


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillTimeSeriesAnalyser(ObjTimeSeriesAnalyserID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjTimeSeriesAnalyserID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTimeSeriesAnalyserID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mTimeSeriesAnalyser_,  Me%InstanceID)

            if (nUsers == 0) then
            
                call KillVariablesAndFiles
            
                !Deallocates Instance
                call DeallocateInstance ()

                ObjTimeSeriesAnalyserID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillTimeSeriesAnalyser
        

    !------------------------------------------------------------------------
    


    subroutine KillVariablesAndFiles
    
    
        !Local--------------------------------------------------------------------------
        integer         :: STAT_CALL
        
        !Begin--------------------------------------------------------------------------
    
        deallocate(Me%DataMatrix, Me%OutData)
        deallocate(Me%TimeTS             )   
        deallocate(Me%PropertyList       ) 
        
        deallocate(Me%TimeTSOutPut       )
        
        if (Me%FilterTimeSerie) then
        
            !closes Output files
            call UnitsManager(Me%iFilter, FileClose, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - KillVariablesAndFiles - ERR10'
        
        
            if (Me%TimeSerieFilterInON) then

                call KillTimeSerie(Me%ObjTimeSerieFilterIn, STAT = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - KillVariablesAndFiles - ERR20'
            
            endif
            
            deallocate(Me%FlagFilter)
                    
        endif

        call UnitsManager(Me%iGap, FileClose, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - KillVariablesAndFiles - ERR30'

        !Open Output files
        call UnitsManager(Me%iInterpol, FileClose, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - KillVariablesAndFiles - ERR40'
        
        
        if (Me%SpectralAnalysis) then
            deallocate(Me%AuxTimeSerie, Me%FFT)
            
            !close Output file
            call UnitsManager(Me%iS, FileClose, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - KillVariablesAndFiles - ERR50'             
        endif
        
        if (Me%PercentileEvolution) then
           
           !close Output file
            call UnitsManager(Me%iPE, FileClose, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - KillVariablesAndFiles - ERR60'
        endif

        if (Me%PercentileAnalysis) then
            
           !close Output file
            call UnitsManager(Me%iP, FileClose, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - KillVariablesAndFiles - ERR70'
        endif
                
    
    end subroutine KillVariablesAndFiles    
    
    !--------------------------------------------------------------------------    
    
   !--------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_TimeSeriesAnalyser), pointer          :: AuxObjTimeSeriesAnalyser
        type (T_TimeSeriesAnalyser), pointer          :: PreviousObjTimeSeriesAnalyser

        !Updates pointers
        if (Me%InstanceID == FirstObjTimeSeriesAnalyser%InstanceID) then
            FirstObjTimeSeriesAnalyser => FirstObjTimeSeriesAnalyser%Next
        else
            PreviousObjTimeSeriesAnalyser => FirstObjTimeSeriesAnalyser
            AuxObjTimeSeriesAnalyser      => FirstObjTimeSeriesAnalyser%Next
            do while (AuxObjTimeSeriesAnalyser%InstanceID /= Me%InstanceID)
                PreviousObjTimeSeriesAnalyser => AuxObjTimeSeriesAnalyser
                AuxObjTimeSeriesAnalyser      => AuxObjTimeSeriesAnalyser%Next
            enddo

            !Now update linked list
            PreviousObjTimeSeriesAnalyser%Next => AuxObjTimeSeriesAnalyser%Next

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

    subroutine Ready (ObjTimeSeriesAnalyser_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTimeSeriesAnalyser_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjTimeSeriesAnalyser_ID > 0) then
            call LocateObjTimeSeriesAnalyser (ObjTimeSeriesAnalyser_ID)
            ready_ = VerifyReadLock (mTimeSeriesAnalyser_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjTimeSeriesAnalyser (ObjTimeSeriesAnalyserID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTimeSeriesAnalyserID

        !Local-----------------------------------------------------------------

        Me => FirstObjTimeSeriesAnalyser
        do while (associated (Me))
            if (Me%InstanceID == ObjTimeSeriesAnalyserID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleTimeSeriesAnalyser - LocateObjTimeSeriesAnalyser - ERR01'

    end subroutine LocateObjTimeSeriesAnalyser

    !--------------------------------------------------------------------------

    end module ModuleTimeSeriesAnalyser