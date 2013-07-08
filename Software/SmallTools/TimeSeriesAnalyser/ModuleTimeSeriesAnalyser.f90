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
    
    
    integer, parameter  :: WeekEnd_ = 1, WeekWay_ = 2
    
    integer, parameter  :: stdv4_   = 1, Average_ = 2
    
    !Types---------------------------------------------------------------------
    type T_MovAveBack
        logical     :: ON
        real        :: P50_Weekly, Raw_Weekly, P50_Daily, Raw_Daily 
    end type T_MovAveBack
    
    type T_TimeSeriesAnalyser
    
        integer                                                 :: InstanceID
    
        complex (SPC),              dimension(:),   allocatable :: FFT
        real     (SP),              dimension(:),   allocatable :: TimeSerie, AuxTimeSerie, amplitude, phase, frequency, FlagTimeSerie 
        real     (SP),              dimension(:),   allocatable :: Tempo, AuxFreq
        
        
	    real,                       dimension(:,:), pointer     :: Hourly_Per
	    real,                       dimension(:,:), pointer     :: DataMatrix, OutData, DataCollection, Aux2D
	    integer,                    dimension(:)  , allocatable :: CollectionSize
	    character(len=Stringlength),dimension(:)  , pointer     :: PropertyList
	    character (Len=PathLength)                              :: InterpolFile, FilterFile
	    character(len=PathLength  )                             :: TimeSerieDataFile, TimeSerieFilterIn, TimeSerieCompareFile
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
	    integer                                                 :: ObjTimeSerieFilterOut = 0	    
	    integer                                                 :: ObjTimeSerieOutNight  = 0

	    integer                                                 :: ObjTimeSerieCompare   = 0
    	
	    integer                                                 :: ObjTimeSerieDiference = 0	
    	
	    integer                                                 :: STAT_CALL, DataValues, DataColumns 
	    integer                                                 :: NValues, DataColumn, nValuesAux
	    integer                                                 :: ObjTimeSerieOut, ObjTime, k, FlagColumn, nFTT
	    integer                                                 :: JulDayStart, JulDayEnd, JulDayWindow, JulDayCenter, JulDayTS
        integer                                                 :: FrequencyWindow
        integer                                                 :: FilterFlagColumn
        real                                                    :: FilterFlagLimit
        logical                                                 :: FilterFlagLimitAbove
        logical                                                 :: CompareTimeSerieOn
        integer                                                 :: CompareColumn
        logical                                                 :: CompareObservations
        !Files Units
        integer                                                 :: iS, iP, iPE, iFilter, iGap, iInterpol, iCompare, iPWeekEnd, iPWeekWay, iMovAveBackWeek, iMovAveBackDay
        
	    logical                                                 :: TimeCycle, OutWindow, SpectralAnalysis, PercentileAnalysis, PercentileEvolution
	    logical                                                 :: FilterTimeSerie
	    logical                                                 :: TimeSerieFilterInON, WeekEndOut
	    real                                                    :: FilterMinValue, FilterMaxValue, FilterMaxRateValue, FilterValue, FilterMinRateValue

        logical                                                 :: RemoveNoisyPeriods, NoValidPoints = .false. 
        real                                                    :: NoisyPeriodAnalysis, NoisyRacioStdevAverage
        integer                                                 :: NoisyPeriodMinSample
                
        integer                                                 :: InterpolInTime = FillValueInt
        
        integer                                                 :: ErrorNormalization 
        
        real                                                    :: GapLimit
        
        real                                                    :: StartNightHour, EndNightHour, MinNightValuesRacio
        
        type (T_MovAveBack)                                     :: MovAverageBackward

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
            
            call ConstructCompareTS
            
            !call ConstructMovingAverageTS
            

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
        integer                 :: status, flag, STAT_CALL, flag_min, flag_max
        
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
        
        if (Me%PercentileAnalysis .and. Me%FrequencyWindow == 168) then
            Me%WeekEndOut = .true.
        else
            Me%WeekEndOut = .false.
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
                         flag_min,                                                      &
                         SearchType   = FromFile,                                       &
                         keyword      ='MIN_VALUE',                                     &
                         Default      = FillValueReal,                                  &
                         ClientModule ='ModuleTimeSeriesAnalyser',                      &
                         STAT         = STAT_CALL)        
            if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR210'

            call GetData(Me%FilterMaxValue,                                             &
                         Me%ObjEnterData,                                               &
                         flag_max,                                                      &
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
            
            if (flag == 1) then
                if (flag_min == 1 .and. flag_max == 1) then
                    Me%FilterMaxRateValue = min (Me%FilterMaxRateValue, (Me%FilterMaxValue - Me%FilterMinValue) / 3600.) 
                endif
            endif
            
            call GetData(Me%FilterMinRateValue,                                         &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='MIN_RATE',                                      &
                         Default      = FillValueReal,                                  &
                         ClientModule ='ModuleTimeSeriesAnalyser',                      &
                         STAT         = STAT_CALL)        
            if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR235'

            if (flag == 1) then
                if (flag_min == 1 .and. flag_max == 1) then
                    Me%FilterMinRateValue = max (Me%FilterMinRateValue, (Me%FilterMaxValue - Me%FilterMinValue) / 3600. * 1e-4) 
                endif
            endif
            
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
            
            call GetData(Me%RemoveNoisyPeriods,                                         &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='REMOVE_NOISY_PERIODS',                          &
                         Default      = .false.,                                        &
                         ClientModule ='ModuleTimeSeriesAnalyser',                      &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR280'
            
            if (Me%RemoveNoisyPeriods) then

                call GetData(Me%NoisyPeriodAnalysis,                                    &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromFile,                                   &
                             keyword      ='NOISY_PERIOD_ANALYSYS',                     &
                             Default      = 86400.,                                     &
                             ClientModule ='ModuleTimeSeriesAnalyser',                  &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR280'
            
                call GetData(Me%NoisyPeriodMinSample,                                   &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromFile,                                   &
                             keyword      ='NOISY_PERIOD_MIN_SAMPLE',                   &
                             Default      = 10,                                         &
                             ClientModule ='ModuleTimeSeriesAnalyser',                  &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR280'
                
                call GetData(Me%NoisyRacioStdevAverage,                                 &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromFile,                                   &
                             keyword      ='NOISY_RACIO_STDEV_AVERAGE',                 &
                             Default      = 0.5,                                        &
                             ClientModule ='ModuleTimeSeriesAnalyser',                  &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR280'
                
            
            endif


        endif
        
        call GetData(Me%FilterValue,                                                    &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='FILTER_VALUE',                                      &
                     Default      = FillValueReal,                                      &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR280'
        
        
        call GetData(Me%CompareTimeSerieOn,                                             &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='COMPARE_TIME_SERIE',                                &
                     default      = .false.,                                            &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR290'
        
        if (Me%CompareTimeSerieOn) then
        
            call GetData(Me%TimeSerieCompareFile,                                       &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='COMPARE_FILE',                                  &
                         ClientModule ='ModuleTimeSeriesAnalyser',                      &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR300'
            if (flag == 0)             stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR310'

            call GetData(Me%CompareColumn,                                              &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='COMPARE_COLUMN',                                &
                         default      = 2,                                              &
                         ClientModule ='ModuleTimeSeriesAnalyser',                      &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR310'
            

            call GetData(Me%CompareObservations,                                        &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='COMPARE_OBSERVATIONS',                          &
                         default      = .true.,                                         &
                         ClientModule ='ModuleTimeSeriesAnalyser',                      &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR320'
            
           
        endif
        
        !options for the error normalization : average value (Average_), 4*standard deviation = 95.4% of the sample in a normal distribution (stdv4_)
        call GetData(Me%ErrorNormalization,                                             &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='ERROR_NORMALIZATION',                               &
                     Default      = stdv4_,                                             &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR330'

        if (Me%ErrorNormalization /= stdv4_ .and. Me%ErrorNormalization /= Average_) then
            stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR340'
        endif

        call GetData(Me%MovAverageBackward%ON,                                          &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='MOVING_AVERAGE_BACKWARD',                           &
                     Default      = .false.,                                            &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR350'
        

        call GetData(Me%StartNightHour,                                                 &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='START_NIGHT_HOUR',                                  &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     default      = 0.,                                                 &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR360'

        call GetData(Me%EndNightHour,                                                   &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='END_NIGHT_HOUR',                                    &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     default      = 4.,                                                 &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR370'

        
        if (Me%StartNightHour >= Me%EndNightHour) then
            Me%EndNightHour = Me%EndNightHour + 24.
        endif

        call GetData(Me%MinNightValuesRacio,                                            &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='MIN_NIGHT_VAlUES_RACIO',                            &
                     ClientModule ='ModuleTimeSeriesAnalyser',                          &
                     default      = 0.75,                                               &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR380'
        
        
        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ReadKeywords - ERR390'

    
    end subroutine ReadKeywords
    
    !-------------------------------------------------------------------------
 
    subroutine ConstructFilterTS
    
        !Local----------------------------------------------------------------
        integer                         :: STAT_CALL

        !Begin----------------------------------------------------------------
    
        if (Me%FilterTimeSerie) then
        
            Me%FilterFile = "FilterOut_"//Me%TimeSerieDataFile 
            !Open Output files
            call UnitsManager(Me%iFilter, FileOpen, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ConstructFilterTS - ERR10'
            
            open(unit   = Me%iFilter, file =trim(Me%FilterFile), form = 'FORMATTED', status = 'UNKNOWN')
            
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
        integer                         :: STAT_CALL

        !Begin----------------------------------------------------------------
        
        Me%InterpolFile = "InterpolOut_"//Me%TimeSerieDataFile 
        !Open Output files
        call UnitsManager(Me%iInterpol, FileOpen, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ConstructInterpolTS - ERR10'
        
        open(unit = Me%iInterpol, file =trim(Me%InterpolFile), form = 'FORMATTED', status = 'UNKNOWN')

        
    end subroutine ConstructInterpolTS        
    
    !-------------------------------------------------------------------------        

    subroutine ConstructPatternsTS
    
        !Local----------------------------------------------------------------
        character (Len=PathLength)      :: SpectralAnalysisFile, PercentileFile,        &
                                           PerEvolutionFile, PercentileWeekEndFile,     &
                                           PercentileWeekWayFile
        
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
                 
            if (Me%WeekEndOut) then

                PercentileWeekEndFile       = "PercentileWeekEndOut_"//Me%TimeSerieDataFile 
                
                !Open Output files
                call UnitsManager(Me%iPWeekEnd, FileOpen, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - Main - ERR30' 
                
                open(unit   = Me%iPWeekEnd, file =trim(PercentileWeekEndFile), form = 'FORMATTED',     &
                     status = 'UNKNOWN')

                PercentileWeekWayFile       = "PercentileWeekWayOut_"//Me%TimeSerieDataFile 
                
                !Open Output files
                call UnitsManager(Me%iPWeekWay, FileOpen, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - Main - ERR40' 
                
                open(unit   = Me%iPWeekWay, file =trim(PercentileWeekWayFile), form = 'FORMATTED',     &
                     status = 'UNKNOWN')
            
            endif                 
            
        endif

        if (Me%PercentileEvolution) then
        
            PerEvolutionFile       = "PerEvolutionOut_"//Me%TimeSerieDataFile 
            
            !Open Output files
            call UnitsManager(Me%iPE, FileOpen, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - Main - ERR50' 
            
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
            if (i>1) then
                if ( Me%TimeTS(i-1) >  Me%TimeTS(i)) then
                    write(*,*) 'Instant number', i, 'older than instant', i-1
                    write(*,*) 'Time series time must be order in ascending way'
                    stop "ModuleTimeSeriesAnalyser - ConstructRawTimeSerie - ERR40"
                endif
            endif
        enddo
        
        
        
    end subroutine ConstructRawTimeSerie        
        
    !-----------------------------------------------------------------------------------        
    
    !-------------------------------------------------------------------------
 
    subroutine ConstructCompareTS
    
        !Local----------------------------------------------------------------
        character (len = PathLength)    :: CompareOutFile
        integer                         :: STAT_CALL

        !Begin----------------------------------------------------------------
    
        if (Me%CompareTimeSerieOn) then
            
            Me%ObjTimeSerieCompare = 0

            call StartTimeSerieInput(Me%ObjTimeSerieCompare,                            &
                                     TimeSerieDataFile = Me%TimeSerieCompareFile,       &
                                     STAT = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ConstructCompareTS - ERR10'
            
            CompareOutFile = "CompareOut_"//Me%TimeSerieDataFile 
            !Open Output files
            call UnitsManager(Me%iCompare, FileOpen, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ConstructCompareTS - ERR10'
            
            open(unit = Me%iCompare, file =trim(CompareOutFile), form = 'FORMATTED', status = 'UNKNOWN')
            
                
        endif

        
    end subroutine ConstructCompareTS        
    
    !-------------------------------------------------------------------------    

    !-------------------------------------------------------------------------
 
    subroutine ConstructMovAverageBackwardTS
    
        !Local----------------------------------------------------------------
        character (len = PathLength)    :: MovAveBackOutFile
        integer                         :: STAT_CALL

        !Begin----------------------------------------------------------------
    
        if (Me%MovAverageBackward%ON) then
            
            MovAveBackOutFile = "MovAveBackDay_"//Me%TimeSerieDataFile 
            !Open Output files
            call UnitsManager(Me%iMovAveBackDay, FileOpen, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ConstructMovAverageBackwardTS - ERR10'
            
            open(unit = Me%iMovAveBackDay, file =trim(MovAveBackOutFile), form = 'FORMATTED', status = 'UNKNOWN')
            
            MovAveBackOutFile = "MovAveBackWeek_"//Me%TimeSerieDataFile 
            !Open Output files
            call UnitsManager(Me%iMovAveBackWeek, FileOpen, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ConstructMovAverageBackwardTS - ERR20'
            
            open(unit = Me%iMovAveBackWeek, file =trim(MovAveBackOutFile), form = 'FORMATTED', status = 'UNKNOWN')
                
        endif

        
    end subroutine ConstructMovAverageBackwardTS        
    
    !-------------------------------------------------------------------------    

    
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
            
            call ModifyTimeSeriesCompare          
            
            call ModifyTimeSeriesMovAveBack
            
            call ModifyTimeSeriesNightMedian
            
            call ModifyTimeSeriesFillAll

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyTimeSeriesAnalyser
    
    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------
    
    subroutine ModifyTimeSeriesNightMedian


        !Arguments-------------------------------------------------------------

        
        !Local-----------------------------------------------------------------
        type (T_Time)                  :: CurrentTime, StartNightTime, EndNightTime, OutNightTime
        real                           :: TS_Value, RelativeOutFlux, NewMass, OldMass, InFlux, Year, Month, Day, NightPeriod, NightPeriodValid
        real                           :: NightAvFlux, NightP50Flux, hour, minute, second
        real, dimension(:  ), pointer  :: WriteAuxNight, SortArray, AuxSort
        integer                        :: STAT_CALL, iN, i50, iP50, jmin, jmax, j, k, iTotal, i, iNight, FilterValues, iSort
        logical                        :: NightPeriodON, NightPeriodOut, TS_Gap
        
        
        !Begin-----------------------------------------------------------------

        call UnitsManager(iNight, FileOpen, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ConstructFilterTS - ERR10'
        
        open(unit   = iNight, file ="Night_"//trim(Me%TimeSerieDataFile), form = 'FORMATTED', status = 'UNKNOWN')
        
        
        write(iNight,*) "NAME                    : ", trim(Me%TimeSerieName)
        write(iNight,*) "LOCALIZATION_I          : -999999"
        write(iNight,*) "LOCALIZATION_J          : -999999"
        write(iNight,*) "LOCALIZATION_K          : -999999"
        call ExtractDate(Me%BeginTime, Year, Month, Day, hour, minute, second)
        write(iNight,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
        write(iNight,*) "TIME_UNITS              : SECONDS"
        write(iNight,*) "COORD_X    : ", Me%CoordX
        write(iNight,*) "COORD_Y    : ", Me%CoordY
        
        FilterValues = Me%DataValues - sum(Me%FlagFilter(1:Me%DataValues))
        
        write(iNight,*) "TOTAL_VALUES  : ", Me%DataValues 
        write(iNight,*) "FILTER_VALUES : ", FilterValues
        write(iNight,'(A25,f6.1)') "GOOD_VALUES_PERCENTAGE : ", real(Me%DataValues - FilterValues) / real(Me%DataValues) * 100.
        
        write(iNight,'(A82)') "Time Night_Average Night_P50"
        
        write(iNight,*) "<BeginTimeSerie>" 
           
        allocate(WriteAuxNight(2))
        
        NightPeriod      = (Me%EndNightHour - Me%StartNightHour) * 3600.
        NightPeriodON    = .false.
       
        NightPeriodOut   = .false.
        NightAvFlux      = 0.
        NightPeriodValid = 0.    
        
        iN               = int(NightPeriod / Me%DT_Analysis) + 1
        iSort            = 0
        
        allocate(SortArray (iN      ))
        
        call ExtractDate(CurrentTime, Year, Month, Day)     
        
        CurrentTime = Me%TimeTSOutPut(1)

        do i=1, Me%nValues 
        
            CurrentTime = Me%TimeTSOutPut(i)
        
            call SetDate(StartNightTime, Year, Month, Day, Me%StartNightHour, 0., 0.)
            call SetDate(EndNightTime  , Year, Month, Day, Me%EndNightHour  , 0., 0.)
            call SetDate(OutNightTime  , Year, Month, Day, (Me%StartNightHour + Me%EndNightHour) / 2. , 0., 0.)            
        
            TS_Value    = Me%TimeSerie(i)
            
            if      (Me%FlagTimeSerie(i) == 0) then
                TS_Gap     = .true.
            else if (Me%FlagTimeSerie(i) == 1) then
                TS_Gap     = .false.
            endif

            if (CurrentTime > StartNightTime .and. CurrentTime < EndNightTime) Then
                if (.not. TS_Gap) then
                    NightPeriodValid = NightPeriodValid + Me%DT_Analysis
                    NightAvFlux      = NightAvFlux + TS_Value * Me%DT_Analysis
                    NightPeriodON    = .true.
                    iSort            = iSort + 1
                    SortArray(iSort) = TS_Value
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
                call ExtractDate(CurrentTime, Year, Month, Day)                   
                               
            endif
            
            if (NightPeriodOut) then
                
                write(iNight,*) CurrentTime-Me%BeginTime, WriteAuxNight(1:2)

                NightPeriodOut   = .false.
                NightAvFlux      = 0.
                NightPeriodValid = 0.   

            endif
 
        enddo
        
    
        write(iNight,*) "<EndTimeSerie>"        
        
        !closes Output files
        call UnitsManager(Me%iFilter, FileClose, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ModifyTimeSeriesNightMedian - ERR20'
        
        deallocate(WriteAuxNight)       
        deallocate(SortArray    ) 
                
    end subroutine ModifyTimeSeriesNightMedian
    
    
   
    !--------------------------------------------------------------------------
        

    !--------------------------------------------------------------------------
    
    subroutine ModifyTimeSeriesFillAll


        !Arguments-------------------------------------------------------------

        
        !Local-----------------------------------------------------------------
        type (T_Time)                  :: CurrentTime, StartNightTime, EndNightTime, OutNightTime
        real                           :: TS_Value, RelativeOutFlux, NewMass, OldMass, InFlux, Year, Month, Day, NightPeriod, NightPeriodValid
        real                           :: NightAvFlux, NightP50Flux, hour, minute, second
        real, dimension(:  ), pointer  :: WriteAux, SortArray, AuxSort
        integer                        :: STAT_CALL, iSort, iN, i50, iP50, jmin, jmax, j, k, iTotal, i, imax, iAll, FilterValues
        logical                        :: NightPeriodON, NightPeriodOut, TS_Gap
        
        
        !Begin-----------------------------------------------------------------

        call UnitsManager(iAll, FileOpen, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ConstructFilterTS - ERR10'
        
        open(unit   = iAll, file ="FillAll_"//trim(Me%TimeSerieDataFile), form = 'FORMATTED', status = 'UNKNOWN')
        
        
        write(iAll,*) "NAME                    : ", trim(Me%TimeSerieName)
        write(iAll,*) "LOCALIZATION_I          : -999999"
        write(iAll,*) "LOCALIZATION_J          : -999999"
        write(iAll,*) "LOCALIZATION_K          : -999999"
        call ExtractDate(Me%BeginTime, Year, Month, Day, hour, minute, second)
        write(iAll,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
        write(iAll,*) "TIME_UNITS              : SECONDS"
        write(iAll,*) "COORD_X    : ", Me%CoordX
        write(iAll,*) "COORD_Y    : ", Me%CoordY
        
        FilterValues = Me%DataValues - sum(Me%FlagFilter(1:Me%DataValues))
        
        write(iAll,*) "TOTAL_VALUES  : ", Me%DataValues 
        write(iAll,*) "FILTER_VALUES : ", FilterValues
        write(iAll,'(A25,f6.1)') "GOOD_VALUES_PERCENTAGE : ", real(Me%DataValues - FilterValues) / real(Me%DataValues) * 100.
        
        write(iAll,'(A82)') "Time FillAll FillAll_P50"
        
        write(iAll,*) "<BeginTimeSerie>" 
           
        allocate(WriteAux(2))
        
        NightPeriod      = (Me%EndNightHour - Me%StartNightHour) * 3600.
        NightPeriodON    = .false.
       
        NightPeriodOut   = .false.
        NightAvFlux      = 0.
        NightPeriodValid = 0.    
        
        iN               = int(NightPeriod / Me%DT_Analysis) + 1
        iSort            = 0
        
        iP50             = int(3600. /Me%DT_Analysis / 2.) + 1
        
        allocate(SortArray (iN      ))
        
        call ExtractDate(CurrentTime, Year, Month, Day)     
        
        CurrentTime = Me%TimeTSOutPut(1)

        allocate(AuxSort(iP50*2+2))

        do i=1, Me%nValues 
        
            CurrentTime = Me%TimeTSOutPut(i)
            TS_Value    = Me%TimeSerie(i)
            
            if      (Me%FlagTimeSerie(i) == 0) then
                TS_Gap     = .true.
            else if (Me%FlagTimeSerie(i) == 1) then
                TS_Gap     = .false.
            endif

            if (TS_Gap) then
                WriteAux(1) = Me%FilterValue
            else
                WriteAux(1) = TS_Value
            endif                
            
            jmin = max(1         ,i-iP50)
            jmax = min(Me%nValues,i+iP50)
            k    = 0
            do j = jmin, jmax
            
                if (Me%FlagTimeSerie(j) == 1) then
                
                    k          = k + 1
                    AuxSort(k) = Me%TimeSerie(j)
                
                endif
            
            enddo
            
            if (k > iP50 .and. jmax-jmin == iP50*2) then

                call sort(AuxSort(1:k))
                i50         = int(k/2) + 1 
                WriteAux(2) = AuxSort(i50)
                
            else
                
                WriteAux(2) = Me%FilterValue
                
            endif
                            
            write(iAll,*) CurrentTime-Me%BeginTime, WriteAux(1:2)

 
        enddo
        
    
        write(iAll,*) "<EndTimeSerie>"        
        
        !closes Output files
        call UnitsManager(Me%iFilter, FileClose, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - ModifyTimeSeriesNightMedian - ERR20'
        
        deallocate(WriteAux )       
        deallocate(SortArray) 
        deallocate(AuxSort  ) 
                
    end subroutine ModifyTimeSeriesFillAll
    
    
   
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    
    subroutine ModifyFilterTimeSeries

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,         pointer, dimension(:,:)       :: FilterDataMatrix
        real,         pointer, dimension(:)         :: AuxValues
        integer,      pointer, dimension(:)         :: AuxI       
        type (T_Time)                               :: FilterInitialData, Time1, Time2, TFinal
        integer                                     :: FilterDataValues, FilterDataColumns
        integer                                     :: i, j, STAT_CALL, k
        real                                        :: Rate, Value1, Value2, Ave, Stv, DT
        logical                                     :: search, StartPeriod, FilterPeriod                
        
        !----------------------------------------------------------------------

        ! Initial value equal to the flag quality of the raw data (if exist)
d21:    do i=1, Me%DataValues 
            if (Me%FlagColumn < 0) then
        
                Me%FlagFilter(i) = 1            
            
            else                
                Me%FlagFilter(i) = Me%DataMatrix(i,Me%FlagColumn)
            endif
        
        enddo d21            
        
        !Filter by value 
d0:     do i=1, Me%DataValues 
            if (Me%FlagFilter(i) == 1) then
                if (Me%DataMatrix(i,Me%DataColumn) == Me%FilterValue) then
                    Me%FlagFilter(i) = 0
                endif
            endif                                
        enddo d0
        

        !Filter by Extreme Values 
d100:   do i=1, Me%DataValues 
            if (Me%FlagFilter(i) == 1) then
                if (abs(Me%DataMatrix(i,Me%DataColumn)) > abs(Me%FilterValue)) then
                    Me%DataMatrix(i,Me%DataColumn) = Me%FilterValue
                    Me%FlagFilter(i) = 0
                endif
            endif                                
        enddo d100

        
i23:    if (Me%RemoveNoisyPeriods) then
        
            StartPeriod = .true.
            
            allocate(AuxValues(1:Me%DataValues), AuxI(1:Me%DataValues))

            !Filter noisy periods
d1:         do i=1, Me%DataValues 

                if (Me%FlagFilter(i) == 1) then

                    if (StartPeriod) then
                        TFinal = Me%TimeTS(i) + Me%NoisyPeriodAnalysis
                        j = 0
                    endif
                    
                    StartPeriod = .false.
                    
                    j = j + 1
                    AuxValues(j) = Me%DataMatrix(i,Me%DataColumn)
                    AuxI     (j) = i

                    if (i < Me%DataValues) then
                        if (Me%TimeTS(i+1) >= TFinal) then                
                            StartPeriod = .true.
                        endif
                    else
                        StartPeriod = .true.                                        
                    endif
                    
                    if (StartPeriod) then
                    
                        FilterPeriod = .false.                
                        
                        if (j < Me%NoisyPeriodMinSample) then
                            FilterPeriod = .true.
                        else
                            Ave = sum(AuxValues(1:j)) / real(j)
                            stv = 0
                            do k=1,j
                                stv = stv + (AuxValues(k) - Ave)**2
                            enddo
                            stv = sqrt(stv / real(j))
                            if (Ave > 0) then
                                if (stv/Ave > Me%NoisyRacioStdevAverage) then
                                    FilterPeriod = .true.
                                endif
                            endif
                        endif
                        
                        if (FilterPeriod) then
                            do k=1,j
                                Me%FlagFilter(AuxI(k)) = 0
                            enddo
                        endif
                                                
                    endif
                                                                            
                endif
                
            enddo d1
            
            deallocate(AuxValues, AuxI)
        
        endif i23        
    
        !Filter points min., máx., spicks, filter period where the sensor sticks to a value (3 consecutive values exactly equal) 
d2:     do i=1, Me%DataValues 
           
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
                DT = Me%TimeTS(i+1)-Me%TimeTS(i)
                if (DT>0.) then
                    Rate = abs(Me%DataMatrix(i+1,Me%DataColumn)-Me%DataMatrix(i,Me%DataColumn))/DT
                    if (Rate > Me%FilterMaxRateValue) then
                        Me%FlagFilter(i  ) = 0
                        Me%FlagFilter(i+1) = 0
                    endif
                    if (Rate < Me%FilterMinRateValue) then
                        Me%FlagFilter(i  ) = 0
                        Me%FlagFilter(i+1) = 0
                    endif
                else
                    Me%FlagFilter(i  ) = 0
                    Me%FlagFilter(i+1) = 0
                endif
            endif

            !Filter persistent constant values
            if (i < Me%DataValues-3) then
                if (Me%DataMatrix(i,Me%DataColumn) == Me%DataMatrix(i+1,Me%DataColumn) .and. &
                    Me%DataMatrix(i,Me%DataColumn) == Me%DataMatrix(i+2,Me%DataColumn)) then 
                    Me%FlagFilter(i  ) = 0
                    Me%FlagFilter(i+1) = 0
                    Me%FlagFilter(i+2) = 0
                endif
            endif

        enddo d2
        
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

d3:         do i=1, Me%DataValues             

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
                
            enddo d3
            
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
        real     (SP),  dimension(:),   allocatable :: TimeSerieFilter
        integer                                     :: FilterValues
        
        !----------------------------------------------------------------------
    
        !Sample the input time serie with a regular time step DT_Analysis
        
        !most be mutiple of 2
        Me%nValues = (Me%EndTime - Me%BeginTime) / Me%DT_Analysis + 1
        
        allocate(Me%TimeSerie    (1:Me%nValues))
        allocate(Me%FlagTimeSerie(1:Me%nValues))
        
        Me%FlagTimeSerie (:) = 0
        
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
                if (Me%FlagColumn < 0) then
            
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
                Me%TimeSerie(i) =  Me%FilterValue
            endif
            
        enddo            
        
        
        do i=1, Me%nValues        
                
            if (Me%FlagTimeSerie(i) == 0)  then
                !lower limit
                j = i-1
                Time1  =  Me%TimeTSOutPut(1)
                Value1 =  Me%FilterValue
                
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
                Value2 =  Me%FilterValue
                
                do while (j<=Me%nValues) 
                    if (Me%FlagTimeSerie(j) == 1) then
                        
                        Time2  =  Me%TimeTSOutPut(j)
                        Value2 =  Me%TimeSerie   (j)
                        exit
                    else
                        j = j + 1
                    endif
                enddo
                
                if (Value1 > Me%FilterValue .and. Value2 > Me%FilterValue) then
                
                    !Interpolates Value for current instant
                    if      (Me%InterpolInTime == LinearTS_) then
                        call InterpolateValueInTime(Me%TimeTSOutPut(i), Time1, Value1,     &
                                                    Time2, Value2, NewValue)     
                    else if (Me%InterpolInTime == BackwardTS_) then
                        NewValue = Value1
                    endif                                                                

                    Me%TSGap    (i) =  Time2 - Time1                                                
                                                           
                else if (Value1 > Me%FilterValue) then
                    NewValue = Value1
                    Me%TSGap(i) = Me%TimeTSOutPut(i) - Time1
                else if (Value2 > Me%FilterValue) then
                    NewValue    = Value2
                    Me%TSGap(i) = Time2 - Me%TimeTSOutPut(i)
                else
                    NewValue    = Value2
                    Me%TimeSerie(i) = 0
                    Me%TSGap(i) = Time2 - Time1
                    Me%nValues  = 0
                    Me%NoValidPoints = .true.
                    exit
                endif                                      
                
                Me%TimeSerie(i) =  NewValue          

            endif            

        enddo
        
        if (.not. Me%NoValidPoints) then
        
            !compute the global statiscal paramters of the time serie
            call moment(Me%TimeSerie,Me%ave,Me%adev,Me%sdev,Me%var,Me%skew,Me%curt)            
        
        
        endif
        
        write(Me%iInterpol,*) "NAME                    : ", trim(Me%TimeSerieName)
        write(Me%iInterpol,*) "LOCALIZATION_I          : -999999"
        write(Me%iInterpol,*) "LOCALIZATION_J          : -999999"
        write(Me%iInterpol,*) "LOCALIZATION_K          : -999999"
        
        call ExtractDate(Me%BeginTime, Year, Month, Day, hour, minute, second)
        
        write(Me%iInterpol,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
        write(Me%iInterpol,*) "TIME_UNITS              : SECONDS"
        write(Me%iInterpol,*) "COORD_X                 : ", Me%CoordX
        write(Me%iInterpol,*) "COORD_Y                 : ", Me%CoordY
        
        write(Me%iInterpol,*) "AVERAGE                 : ",Me%ave
        write(Me%iInterpol,*) "AVERAGE_DEVIATION       : ",Me%adev
        write(Me%iInterpol,*) "STANDARD_DEVIATION      : ",Me%sdev
        write(Me%iInterpol,*) "VARIANCE                : ",Me%var
        write(Me%iInterpol,*) "SKEWNESS                : ",Me%skew
        write(Me%iInterpol,*) "KURTOSIS                : ",Me%curt   
        
        write(Me%iInterpol,'(A21)') "Time InterpolatedData"
        
        write(Me%iInterpol,*) "<BeginTimeSerie>"         
        
        do i=1, Me%nValues        
            write(Me%iInterpol,*) Me%TimeTSOutPut(i)-Me%BeginTime, Me%TimeSerie(i)
        enddo
        
        write(Me%iInterpol,*) "<EndTimeSerie>"      
        

        !Open Output files
        call UnitsManager(Me%iInterpol, FileClose, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "TimeSeriesAnalyser - ConstructTimeSeriesAnalyser - ERR40"

!        call StartTimeSerieInput(Me%ObjTimeSerieInterpol,                       &
!                                 TimeSerieDataFile = Me%InterpolFile, STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) stop "TimeSeriesAnalyser - ConstructTimeSeriesAnalyser - ERR50"        
                    

        do i=1, Me%nValues        
            if (Me%TSGap(i) > Me%GapLimit) then
                Me%FlagTimeSerie(i) = 0
            endif
        enddo
        
        if (Me%FilterTimeSerie) then
        
            write(Me%iFilter,*) "NAME                    : ", trim(Me%TimeSerieName)
            write(Me%iFilter,*) "LOCALIZATION_I          : -999999"
            write(Me%iFilter,*) "LOCALIZATION_J          : -999999"
            write(Me%iFilter,*) "LOCALIZATION_K          : -999999"
            call ExtractDate(Me%BeginTime, Year, Month, Day, hour, minute, second)
            write(Me%iFilter,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
            write(Me%iFilter,*) "TIME_UNITS              : SECONDS"
            write(Me%iFilter,*) "COORD_X    : ", Me%CoordX
            write(Me%iFilter,*) "COORD_Y    : ", Me%CoordY
            
            FilterValues = Me%DataValues - sum(Me%FlagFilter(1:Me%DataValues))
            
            write(Me%iFilter,*) "TOTAL_VALUES  : ", Me%DataValues 
            write(Me%iFilter,*) "FILTER_VALUES : ", FilterValues
            write(Me%iFilter,'(A25,f6.1)') "GOOD_VALUES_PERCENTAGE : ", real(Me%DataValues - FilterValues) / real(Me%DataValues) * 100.
            
            write(Me%iFilter,'(A82)') "Time FilterData"
            
            write(Me%iFilter,*) "<BeginTimeSerie>" 
               
            do i=1, Me%nValues        
                if (Me%FlagTimeSerie(i) == 1) then
                    write(Me%iFilter,*) Me%TimeTSOutPut(i)-Me%BeginTime, Me%TimeSerie(i)
                endif
            enddo
        
            write(Me%iFilter,*) "<EndTimeSerie>"  
            
            !closes Output files
            call UnitsManager(Me%iFilter, FileClose, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop "TimeSeriesAnalyser - ConstructTimeSeriesAnalyser - ERR60"
        
            !call StartTimeSerieInput(Me%ObjTimeSerieFilterOut,                          &
            !                         TimeSerieDataFile = Me%FilterFile, STAT = STAT_CALL)        
            !if (STAT_CALL /= SUCCESS_) stop "TimeSeriesAnalyser - ConstructTimeSeriesAnalyser - ERR70"     
            

        endif
        
        
    end subroutine ModifyInterpolTimeSeries
    
    !--------------------------------------------------------------------------
    
    subroutine ModifyTimeSeriesCompare
    
        real    (SP), dimension(:), allocatable :: TimeCompareTS, CompareTS, DifTS, TS
        real                                    :: ave,adev,sdev,var,skew,curt, RMSE
        type (T_Time)                           :: Time1, Time2
        real                                    :: Value1, Value2, NewValue, AverageCompare, AverageObs
        real                                    :: Year, Month, Day, hour, minute, second
        real                                    :: a, b, d, t1, t2, NashS, Skill
        logical                                 :: TimeCycle
        integer                                 :: i, c, STAT_CALL, Ncompare
        
        !Begin-------------------------------------------------------------------------        
             
i1:     if (Me%CompareTimeSerieOn) then

            allocate(TimeCompareTS(1: Me%nValues), TS(1:Me%nValues),CompareTS (1:Me%nValues), DifTS(1:Me%nValues))   

            c = 0
            AverageCompare = 0.

            do i=1, Me%nValues        

                if (Me%FilterTimeSerie) then
                    if (Me%FlagTimeSerie(i) == 0) cycle
                endif
                
                c = c + 1

                call GetTimeSerieValue(Me%ObjTimeSerieCompare, Me%TimeTSOutPut(i), Me%CompareColumn, Time1, Value1,   &
                                       Time2, Value2, TimeCycle, STAT= STAT_CALL) 
                
                !Interpolates Value for current instant
                if      (Me%InterpolInTime == LinearTS_) then
                    call InterpolateValueInTime(Me%TimeTSOutPut(i), Time1, Value1, Time2,  &
                                                Value2, NewValue)
                else if (Me%InterpolInTime == BackwardTS_) then
                    NewValue = Value1
                endif   
                
                TimeCompareTS(c) =  Me%TimeTSOutPut(i) -Me%BeginTime
                CompareTS    (c) =  NewValue
                AverageCompare   =  AverageCompare + CompareTS    (c)
                TS           (c) =  Me%TimeSerie(i)
                DifTS        (c) =  TS(c) - CompareTS(c)
                
            enddo                        

            Ncompare = c        
            
            if (Ncompare > 0)AverageCompare = AverageCompare / real(Ncompare)
                
            call moment (DifTS(1:Ncompare),ave,adev,sdev,var,skew,curt)
            
            RMSE = 0.
            a    = 0
            b    = 0
            d    = 0
            if (Me%CompareObservations) then
                AverageObs = AverageCompare
            else
                AverageObs = Me%ave
            endif
            
            do c=1,Ncompare
            
                a = a + DifTS(c)**2
                
                t1= CompareTS(c)-AverageObs
                t2= TS(c)       -AverageObs
                
                if (Me%CompareObservations) then
                    b = b + t1**2
                else
                    b = b + t2**2
                endif
                
                d = d + (abs(t1)+abs(t2))**2
                
            enddo
            
            RMSE = sqrt(a/NCompare)
            
            if (b>0.) then
                NashS  = 1 - a/b
            else
                NashS  = 0.
            endif            
            
            
            if (d>0.) then
                Skill = 1 - a/d
            else
                Skill = 0.
            endif

            write(Me%iCompare,*) "NAME                    : ", trim(Me%TimeSerieName)
            write(Me%iCompare,*) "LOCALIZATION_I          : -999999"
            write(Me%iCompare,*) "LOCALIZATION_J          : -999999"
            write(Me%iCompare,*) "LOCALIZATION_K          : -999999"
            
            call ExtractDate(Me%BeginTime, Year, Month, Day, hour, minute, second)
            
            write(Me%iCompare,'(A26,5F6.0,1f8.2)')                                          &
                                 "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
            write(Me%iCompare,*) "TIME_UNITS              : SECONDS"
            write(Me%iCompare,*) "COORD_X                 : ", Me%CoordX
            write(Me%iCompare,*) "COORD_Y                 : ", Me%CoordY
            
            write(Me%iCompare,*) "BIAS                    : ",ave
            write(Me%iCompare,*) "RMSE                    : ",RMSE
            write(Me%iCompare,*) "AVERAGE_DEVIATION       : ",adev
            write(Me%iCompare,*) "STANDARD_DEVIATION      : ",sdev
            write(Me%iCompare,*) "VARIANCE                : ",var
            write(Me%iCompare,*) "SKEWNESS                : ",skew
            write(Me%iCompare,*) "KURTOSIS                : ",curt   
            write(Me%iCompare,*) "NASH–SUTCLIFFE          : ",NashS
            write(Me%iCompare,*) "SKILL                   : ",Skill         
            
            write(Me%iCompare,'(A35)') "Time TimeSerie CompareTimeSerie Dif"
            
            write(Me%iCompare,*) "<BeginTimeSerie>"         
            
            do c=1, NCompare        
                write(Me%iCompare,*) TimeCompareTS(c), TS(c), CompareTS(c), DifTS(c)
            enddo
            
            write(Me%iCompare,*) "<EndTimeSerie>"      
            
            deallocate(TimeCompareTS, CompareTS, DifTS, TS)   
            
            !Close Output files
            call UnitsManager(Me%iCompare, FileClose, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop "TimeSeriesAnalyser - ModifyTimeSeriesCompare - ERR10"        

        endif i1
    
    end subroutine ModifyTimeSeriesCompare
    
    !--------------------------------------------------------------------------    
       
    subroutine ModifyTimeSeriesPatterns
    
        !Local--------------------------------------------------------------
        
        !Begin--------------------------------------------------------------
        
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
        type (T_Time), dimension(:),   pointer :: TimeInstant 
        real,          dimension(:,:), pointer :: Hourly_StandardHyd
        real,          dimension(:),   pointer :: ValuesOut, ValuesIn
        real                                   :: Year, Month, Day, hour, minute, second
        integer                                :: i, j, NPerHour
        
        !Begin--------------------------------------------------------------

        NPerHour = 2*int(Me%nValues/Me%FrequencyWindow) 
        
        allocate(Me%Hourly_Per    (0:Me%FrequencyWindow-1,1:11)    )
        allocate(Me%DataCollection(0:Me%FrequencyWindow-1,NPerHour))
        allocate(Me%CollectionSize(0:Me%FrequencyWindow-1)         )        
        
        call ComputePercentile(Me%EndTime, Me%BeginTime, Me%Hourly_Per) 
        
        write(Me%iP,*) "NAME                    :",  trim(Me%TimeSerieName)
        write(Me%iP,*) "LOCALIZATION_I          : -999999"
        write(Me%iP,*) "LOCALIZATION_J          : -999999"
        write(Me%iP,*) "LOCALIZATION_K          : -999999"
        call ExtractDate(Me%BeginTime, Year, Month, Day, hour, minute, second)
        write(Me%iP,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
        write(Me%iP,*) "TIME_UNITS              : HOURS"
        write(Me%iP,*) "COORD_X    : ", Me%CoordX
        write(Me%iP,*) "COORD_Y    : ", Me%CoordY
        write(Me%iP,*) "TIME_CYCLE : 1"
        
        write(Me%iP,'(A82)') "Time Per_0 Per_10 Per_20 Per_30 Per_40 Per_50 Per_60 Per_70 Per_80 Per_90  Per_100"
        
        write(Me%iP,*) "<BeginTimeSerie>" 
        

        do j=0,Me%FrequencyWindow-1
            write(Me%iP,'(I4,11(f14.6," "))') j, Me%Hourly_Per(j,1:11)
        enddo
        
        write(Me%iP,*) "<EndTimeSerie>"  
        
        
        if (Me%WeekEndOut) then

            allocate(Hourly_StandardHyd    (0:23,1:11))

            allocate(TimeInstant           (0:Me%FrequencyWindow-1))
            allocate(ValuesOut             (0:23)                  )
            allocate(ValuesIn              (0:Me%FrequencyWindow-1))
            
            do j=0,Me%FrequencyWindow-1
                TimeInstant(j) = Me%BeginTime + real(j) * 3600.
            enddo
            
            do i=1,11
                ValuesIn(:) = Me%Hourly_Per(:,i)
                call DailyStandardHydrographs(Me%FrequencyWindow, TimeInstant, ValuesIn, ValuesOut, Option= WeekEnd_)
                Hourly_StandardHyd(:,i) = ValuesOut(:)
            enddo
            
            write(Me%iPWeekEnd,*) "NAME                    :",  trim(Me%TimeSerieName)
            write(Me%iPWeekEnd,*) "LOCALIZATION_I          : -999999"
            write(Me%iPWeekEnd,*) "LOCALIZATION_J          : -999999"
            write(Me%iPWeekEnd,*) "LOCALIZATION_K          : -999999"
            call ExtractDate(Me%BeginTime, Year, Month, Day, hour, minute, second)
            write(Me%iPWeekEnd,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, 0., 0., 0.
            write(Me%iPWeekEnd,*) "TIME_UNITS              : HOURS"
            write(Me%iPWeekEnd,*) "COORD_X    : ", Me%CoordX
            write(Me%iPWeekEnd,*) "COORD_Y    : ", Me%CoordY
            
            write(Me%iPWeekEnd,'(A82)') "Time Per_0 Per_10 Per_20 Per_30 Per_40 Per_50 Per_60 Per_70 Per_80 Per_90  Per_100"
            
            write(Me%iPWeekEnd,*) "<BeginTimeSerie>" 

            
            do j=0,23
                write(Me%iPWeekEnd,'(I4,11(f14.6," "))') j, Hourly_StandardHyd(j,1:11)
            enddo
            
            write(Me%iPWeekEnd,*) "<EndTimeSerie>"  
            

            do i=1,11
                ValuesIn(:) = Me%Hourly_Per(:,i)
                call DailyStandardHydrographs(Me%FrequencyWindow, TimeInstant, ValuesIn, ValuesOut, Option= WeekWay_)
                Hourly_StandardHyd(:,i) = ValuesOut(:)
            enddo
            
            write(Me%iPWeekWay,*) "NAME                    :",  trim(Me%TimeSerieName)
            write(Me%iPWeekWay,*) "LOCALIZATION_I          : -999999"
            write(Me%iPWeekWay,*) "LOCALIZATION_J          : -999999"
            write(Me%iPWeekWay,*) "LOCALIZATION_K          : -999999"
            call ExtractDate(Me%BeginTime, Year, Month, Day, hour, minute, second)
            write(Me%iPWeekWay,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, 0., 0., 0.
            write(Me%iPWeekWay,*) "TIME_UNITS              : HOURS"
            write(Me%iPWeekWay,*) "COORD_X    : ", Me%CoordX
            write(Me%iPWeekWay,*) "COORD_Y    : ", Me%CoordY
            
            write(Me%iPWeekWay,'(A82)') "Time Per_0 Per_10 Per_20 Per_30 Per_40 Per_50 Per_60 Per_70 Per_80 Per_90  Per_100"
            
            write(Me%iPWeekWay,*) "<BeginTimeSerie>" 

            
            do j=0,23
                write(Me%iPWeekWay,'(I4,11(f14.6," "))') j, Hourly_StandardHyd(j,1:11)
            enddo
            
            write(Me%iPWeekWay,*) "<EndTimeSerie>"  


            deallocate(Hourly_StandardHyd)
            deallocate(TimeInstant       )
            deallocate(ValuesOut         )
            deallocate(ValuesIn          )
        
        endif

        deallocate(Me%Hourly_Per    )
        deallocate(Me%DataCollection)
        deallocate(Me%CollectionSize)

    end subroutine ComputePercentileAnalysis      

    !--------------------------------------------------------------------------    
    
        

    !--------------------------------------------------------------------------    
    
    subroutine DailyStandardHydrographs(NValues, TimeInstant, ValuesIn, ValuesOut, Option)             

        !Arguments----------------------------------------------------------
        type (T_Time), dimension(:),   pointer :: TimeInstant 
        real,          dimension(:),   pointer :: ValuesOut, ValuesIn
        integer                                :: NValues, Option

        !Local--------------------------------------------------------------
        real                                   :: Hour
        integer, dimension(0:23)               :: Counter
        integer                                :: i, j

        
        !Begin--------------------------------------------------------------
        
        Counter  (:) = 0 
        ValuesOut(:) = 0
        
        do i=0, NValues-1
            
            call ExtractDate(TimeInstant(i), Hour = Hour)
            if      (Option == WeekWay_) then
                if (WeekEndDay(TimeInstant(i))) cycle
            else if (Option == WeekEnd_) then
                if (WeekWayDay(TimeInstant(i))) cycle
            else
                write(*,*) "Not valid option"
                stop "DailyStandardHydrographs - TimeSeriesAnalyser - ERR10"
            endif
            
            j            = int(Hour)
            
            ValuesOut(j) = ValuesOut(j) + ValuesIn(i)
            Counter  (j) = Counter(j) + 1

        enddo
        
        do j=0,23
            if (Counter(j) > 0) then
                ValuesOut(j) = ValuesOut(j) / real(Counter(j))
            endif    
        enddo
        
        
    end subroutine DailyStandardHydrographs


    !--------------------------------------------------------------------------    
    
    subroutine ComputePercentileEvolution               

        !Local--------------------------------------------------------------
        real                    :: Year, Month, Day, hour, minute, second, aux, dti
        integer                 :: i, j, NPerHour, k, istart, ifinal, j1, j2, iValues
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
        write(Me%iPE,*) "COORD_X    :  -4.140000    "
        write(Me%iPE,*) "COORD_Y    :   36.71000     "
        
        write(Me%iPE,'(A256)') "Time Per_0 Per_10 Per_20 Per_30 Per_40 Per_50 Per_60 Per_70 Per_80 Per_90  Per_100 Raw_Data stdv RMS Error[%]"
        
        write(Me%iPE,*) "<BeginTimeSerie>" 

        iValues = 0
        
        do j=istart,ifinal
        
            if (Me%FlagFilter(j) == 0) then
                cycle
            else
                iValues = iValues + 1
            endif
            
            !Compute Per_50 error
            Me%RMS  = 0
            Me%stdv = 0
            
            Me%Average=sum(Me%OutData(1:j,13)*Me%FlagFilter(1:j))/real(iValues)
            
            do k=1, j
                if (Me%FlagFilter(k) == 0) cycle

                Me%stdv = Me%stdv + (Me%OutData(k,13)-Me%Average     )**2
                Me%RMS  = Me%RMS  + (Me%OutData(k,13)-Me%OutData(k,7))**2
                
            enddo    
            
            Me%RMS  = sqrt(Me%RMS / real(iValues))
            Me%stdv = sqrt(Me%stdv/ real(iValues))
            !Error in percentage normalized using the standard deviation multiply by 4. In a normal distribution 
            !corresponds to ~95% of the sample 
            if (Me%ErrorNormalization == stdv4_) then
                if (Me%stdv >0) then
                    Me%Error = Me%RMS / (4*Me%stdv) * 100.
                else
                    !standard deviation equal zero. Error is assume equal to 100%
                    Me%Error = 100
                endif
            elseif (Me%ErrorNormalization == Average_) then
                if (abs(Me%Average) >0) then
                    Me%Error = Me%RMS / abs(Me%Average)* 100. 
                else
                    Me%Error = 100
                endif            
            endif    
                
            write(Me%iPE,'(f12.2,16(e14.6," "))') Me%OutData(j,1:13), Me%stdv, Me%RMS, Me%Error
        enddo
        
        write(Me%iPE,*) "<EndTimeSerie>"  
        
        deallocate(Me%Hourly_Per    )
        deallocate(Me%DataCollection)
        deallocate(Me%CollectionSize)        

    end subroutine ComputePercentileEvolution  
    

    !------------------------------------------------------------------
    

    !--------------------------------------------------------------------------    
    
    subroutine ModifyTimeSeriesMovAveBack               

    end subroutine ModifyTimeSeriesMovAveBack
    
    !-----------------------------------------------------------------------  
        

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
        
        Hourly_Per_Aux(:,:) = Me%FilterValue
        
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
        
!        call KillTimeSerie(Me%ObjTimeSerieInterpol, STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) stop "TimeSeriesAnalyser - KillVariablesAndFiles - ERR10"          
        
        if (Me%FilterTimeSerie) then

            !call KillTimeSerie(Me%ObjTimeSerieFilterOut, STAT = STAT_CALL)        
            !if (STAT_CALL /= SUCCESS_) stop "TimeSeriesAnalyser - KillVariablesAndFiles - ERR20"
        
            if (Me%TimeSerieFilterInON) then

                call KillTimeSerie(Me%ObjTimeSerieFilterIn, STAT = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - KillVariablesAndFiles - ERR30'
                
            endif
            
            deallocate(Me%FlagFilter)
                    
        endif

        call UnitsManager(Me%iGap, FileClose, STAT = STAT_CALL)
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
            
            if (Me%WeekEndOut) then
               !close Output file
                call UnitsManager(Me%iPWeekEnd, FileClose, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - KillVariablesAndFiles - ERR80'
                
               !close Output file
                call UnitsManager(Me%iPWeekWay, FileClose, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - KillVariablesAndFiles - ERR90'                
            endif
        endif
        
        if (Me%CompareTimeSerieOn) then
            
            call KillTimeSerie(Me%ObjTimeSerieCompare, STAT = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTimeSeriesAnalyser - KillVariablesAndFiles - ERR100'
                
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