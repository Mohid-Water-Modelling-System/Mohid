!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : Gauge
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to calculate waterlevel at gauges
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

Module ModuleGauge

    use ModuleGlobalData
    use ModuleTime                  
    use ModuleEnterData
    use ModuleToga               
    use ModuleHorizontalGrid,   only : GetGridAngle, GetGridOrigin, GetHorizontalGrid,   &
                                       GetXYCellZ, UnGetHorizontalGrid
    use ModuleTimeSerie,        only : StartTimeSerieInput, GetTimeSerieValue, KillTimeSerie
    use ModuleFunctions,        only : RodaXY
    use ModuleTask2000,         only : Task2000Level          

    implicit none

    private 

    !Subroutine----------------------------------------------------------------

    !Constructor
    public  :: ConstructGauges
    private ::      ReadGaugeData
    private ::          Calc_Decimal_geo_coord
    private ::          ReadHarmonics
    private ::          NewWave
    private ::      ReadOldForm
    private ::      VerifyGaugeLocation


    !Selector
    public  :: GetNGauges
    public  :: GetNComponents
    public  :: GetMetricGaugeLocation
    public  :: GetGaugeLocation
    private ::      GetNumberWaves
    public  :: GetReferenceLevel
    public  :: GetLevelEvolution
    public  :: GetIJWaterLevel
    public  :: GetIJReferenceLevel
    public  :: GetTriangGaugesON
    public  :: GetVelEvolution

    !Modifier
    public  :: GaugeLevel


    !Destructor
    public  :: KillGauge
    private ::      KillGaugeList
    private ::          KillTidalWaveList


    !Management
    private ::      Ready
    private ::          LocateObjGauge

    !Parameter
    integer, parameter :: WaveNameLength = 5
    integer, parameter :: NComponents   = 146

    !Evolution Options
    character(LEN = StringLength), parameter :: Char_Harmonics   = 'Harmonics',          &
                                                Char_TimeSerie   = 'Time Serie',         &
                                                Char_NoEvolution = 'No',                 &
                                                Char_Constant    = 'Constant'
    integer, parameter                       :: Constant = 0, Harmonics = 1,             &
                                                TimeSerie = 2, NoEvolution = 3

    !prevision method 
    integer, parameter                       :: Task2000_ = 1, Toga_ = 2
    
    !Type
    type T_TidalWave
        character(LEN = WaveNameLength) :: Name         = null_str
        real                            :: Amplitude
        real                            :: Phase
        type(T_TidalWave), pointer      :: Prev
        type(T_TidalWave), pointer      :: Next
    end type T_TidalWave


    type T_TideGauge 
        character (len=StringLength)        :: Name             = null_str
        real, dimension(3)                  :: Longitude        = FillValueReal
        real, dimension(3)                  :: Latitude         = FillValueReal
        real                                :: DecimalLatitude  = FillValueReal
        real                                :: DecimalLongitude = FillValueReal

        real                                :: Metric_X         = FillValueReal
        real                                :: Metric_Y         = FillValueReal

        logical                             :: DefinedMetric    = .false.

        integer                             :: Grid_I           = FillValueInt
        integer                             :: Grid_J           = FillValueInt

        real                                :: ReferenceLevel   = FillValueReal
        real                                :: WaterLevel       = FillValueReal
        real                                :: TimeReference    = FillValueReal
        integer                             :: LevelEvolution   = FillValueInt
        integer                             :: RefLevelEvolution= FillValueInt
        real                                :: RefLevelAmpFactor= FillValueReal        
        real                                :: RefLevelPhase    = FillValueReal                
        integer                             :: LevelColumn      = FillValueInt
        integer                             :: RefLevelColumn   = FillValueInt

        logical                             :: NoCoveredColumn

        integer                             :: VelEvolution     = FillValueInt
        integer                             :: VelUColumn       = FillValueInt
        integer                             :: VelVColumn       = FillValueInt


        integer                             :: CoveredColumn    = FillValueInt
        character (len=StringLength)        :: TimeSerieDataFile= null_str
        
        type(T_TidalWave), pointer          :: FirstWave
        type(T_TidalWave), pointer          :: LastWave

        type(T_TideGauge), pointer          :: Prev
        type(T_TideGauge), pointer          :: Next 

        !Instance of ModuleToga
        integer                             :: ObjToga      = 0
        integer                             :: ObjTimeSerie = 0
    end type T_TideGauge


    type       T_Gauge
        integer                             :: InstanceID
        integer                             :: UnitGauge    = FillValueInt
        character(LEN = StringLength)       :: FileName     = null_str

        !Instance of ModuleTime
        integer                             :: ObjTime = 0

        !Instance of Module_EnterData
        integer                             :: ObjEnterData = 0
        
        logical                             :: Triangulation
        
        integer                             :: TidePrevision
       
        !Link list
        type (T_TideGauge ), pointer        :: FirstGauge 
        type (T_TideGauge ), pointer        :: LastGauge

        type (T_Gauge), pointer             :: Next

    end type T_Gauge


    !Global Module Variables
    type (T_Gauge), pointer                 :: FirstGauge
    type (T_Gauge), pointer                 :: Me

    !--------------------------------------------------------------------------

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructGauges(GaugeID, HorizontalGridID, TimeID, GaugeFile, STAT)

        !Arguments------------------------------------------------------------
        integer                                     :: GaugeID
        integer,          optional,  intent(IN)     :: HorizontalGridID
        integer                                     :: TimeID
        character(len=*), optional,  intent(IN)     :: GaugeFile
        integer,          optional,  intent(OUT)    :: STAT     

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: NWaves, iflag
        logical                                     :: BlockFound
        character(LEN = WaveNameLength), pointer, dimension(:) :: WaveName
        type(T_TideGauge), pointer                  :: PresentGauge
        type(T_TidalWave), pointer                  :: PresentWave
        integer                                     :: I, ClientNumber
        integer                                     :: STAT_, ready_
        character(LEN = StringLength   ), parameter :: block_begin = '<begingauge>'
        character(LEN = StringLength   ), parameter :: block_end   = '<endgauge>'

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_


        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mGauge_)) then
            nullify (FirstGauge)
            call RegisterModule (mGauge_) 
        endif

        call Ready(GaugeID, ready_)    

        if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            !Nullify Applications
            nullify (Me%FirstGauge)
            nullify (Me%LastGauge )

            !Associates other instances
            Me%ObjTime = AssociateInstance (mTIME_, TimeID)

            if(present(GaugeFile))then

                Me%FileName = GaugeFile
                    
            else
                
                !Reads filename of the Tides file
                call ReadFileName('IN_TIDES', Me%FileName, Message = 'Tides file', STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauges - ERR10' 

            end if

            !Opens file with the construct data 
            call ConstructEnterData(Me%ObjEnterData, Me%FileName, STAT = STAT_CALL)           
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauges - ERR20' 
            
            call GetData(Me%Triangulation,                                                  &
                         Me%ObjEnterData, iflag,                                            &
                         SearchType = FromFile,                                             &
                         keyword    = 'TRIANGULATION',                                      &
                         Default    = .true.,                                               &                                           
                         ClientModule ='ModuleGauge',                                       &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauges - ERR30' 
                
            call GetData(Me%TidePrevision,                                                  &
                         Me%ObjEnterData, iflag,                                            &
                         SearchType = FromFile,                                             &
                         keyword    = 'PREVISION',                                          &
                         Default    = Task2000_,                                            &                                           
                         ClientModule ='ModuleGauge',                                       &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauges - ERR35' 
            
            if (Me%TidePrevision /= Task2000_ .and. Me%TidePrevision /= Toga_) then
                write(*,*) 'The tide prevision method is not known'
                stop       'ConstructGauges - ModuleGauges - ERR38'
            endif


            !New Gauger file format
do2:        do 
                call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,  &
                                            block_begin, block_end, BlockFound, &
                                            STAT = STAT_CALL)
if4 :           if (STAT_CALL .EQ. 0) then
if5 :               if (BlockFound) then           
                        call ReadGaugeData(PresentGauge)

if2 :                   if (.NOT. associated(Me%FirstGauge)) then         !If list is not initialized
                            Me%FirstGauge => PresentGauge                 !FirstGauge element is always FirstGauge
                            Me%LastGauge  => PresentGauge                 !FirstGauge element is also the LastGauge one

                            nullify(Me%FirstGauge%Prev)
                            nullify(Me%LastGauge%Next )
                        else 
                            PresentGauge%Prev => Me%LastGauge

                            Me%LastGauge%Next => PresentGauge 
                            Me%LastGauge      => PresentGauge

                            nullify(Me%LastGauge%Next)
                        end if if2
                    else
                        exit do2    !No more blocks
                    end if if5
                else if (STAT_CALL .EQ. BLOCK_END_ERR_) then 
                    stop 'ConstructGauges - ModuleGauges - ERR40' 
                end if if4
            end do do2


            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauges - ERR50' 


            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauges - ERR60' 

            !If no Gaugers were found in the new file format tryes to read former file format.
if3 :       if (.NOT. associated(Me%FirstGauge)) then
do3:            do 
                    call ReadOldForm(BlockFound, PresentGauge)

if6 :               if (BlockFound) then
if7 :                   if      (.NOT. associated(Me%FirstGauge)) then    !If list is not initialized
                            Me%FirstGauge => PresentGauge                 !FirstGauge element is always FirstGauge
                            Me%LastGauge  => PresentGauge                 !FirstGauge element is also the LastGauge one

                            nullify(Me%FirstGauge%Prev)
                            nullify(Me%LastGauge%Next )

                        else if (      associated(Me%FirstGauge)) then
                            stop 'ConstructGauges - ModuleGauges - ERR70' 
                        end if if7
                    else
                        exit do3    !No more blocks
                    end if if6
                end do do3
            end if if3



if8 :       if (.NOT. associated(Me%FirstGauge)) then
                stop 'ConstructGauges - ModuleGauges - ERR80' 
            else
                PresentGauge => Me%FirstGauge
do4 :           do while (associated(PresentGauge))
                    call GetNumberWaves(PresentGauge, NWaves)

                    allocate(WaveName(NWaves), STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauge - ERR90.'


                    I = 0
                    PresentWave => PresentGauge%FirstWave
do5 :               do while (associated(PresentWave))
                        I = I + 1

                        WaveName(I) = PresentWave%Name

                        PresentWave => PresentWave%Next
                    end do do5

                    call ConstructToga(PresentGauge%ObjToga, NWaves, WaveName, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauge - ERR100'

                    deallocate(WaveName)
                    nullify   (WaveName)

    
                    !Verifies the location of teh Gauge
                    if(present(HorizontalGridID)) call VerifyGaugeLocation(PresentGauge, HorizontalGridID)
                    
                    !Points to the next Gauge in the list
                    PresentGauge => PresentGauge%Next
                end do do4
            end if if8

            call LocationConsistence

            !Returns ID
            GaugeID = Me%InstanceID
            STAT_   = SUCCESS_

        else 
            
            stop 'ModuleGauge - ConstructGauges - ERR110' 

        end if 

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructGauges

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
        type (T_Gauge), pointer             :: NewObjGauge
        type (T_Gauge), pointer             :: PreviousObjGauge


        !Allocates new instance
        allocate (NewObjGauge)
        nullify  (NewObjGauge%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstGauge)) then
            FirstGauge            => NewObjGauge
            Me                    => NewObjGauge
        else
            PreviousObjGauge      => FirstGauge
            Me                    => FirstGauge%Next
            do while (associated(Me))
                PreviousObjGauge  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjGauge
            PreviousObjGauge%Next => NewObjGauge
        endif

        Me%InstanceID = RegisterNewInstance (mGAUGE_)

    end subroutine AllocateInstance
    
    !--------------------------------------------------------------------------

    subroutine Calc_Decimal_geo_coord (PresentGauge)

        !Arguments------------------------------------------------------------
        type(T_TideGauge), pointer :: PresentGauge

        !----------------------------------------------------------------------

        if(PresentGauge%latitude (1) .ge. 0)then
            PresentGauge%DecimalLatitude  = PresentGauge%latitude (1) + PresentGauge%latitude (2)/60.0 + &
                                            PresentGauge%latitude (3)/3600.0

        else
            PresentGauge%DecimalLatitude  = PresentGauge%latitude (1) - PresentGauge%latitude (2)/60.0 - &
                                            PresentGauge%latitude (3)/3600.0

        endif


        if(PresentGauge%longitude (1) .ge. 0)then

            PresentGauge%DecimalLongitude = PresentGauge%longitude(1) + PresentGauge%longitude(2)/60.0 + &
                                            PresentGauge%longitude(3)/3600.0

        else
            
            PresentGauge%DecimalLongitude = PresentGauge%longitude(1) - PresentGauge%longitude(2)/60.0 - &
                                            PresentGauge%longitude(3)/3600.0

        end if



        !----------------------------------------------------------------------

    end subroutine Calc_Decimal_geo_coord

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine ReadGaugeData(PresentGauge)

        !Arguments------------------------------------------------------------
        type(T_TideGauge), pointer                  :: PresentGauge

        !Local-----------------------------------------------------------------
        real                                        :: DT, ReadTimeSerieDT
        integer, pointer, dimension(:)              :: ColumnsToRead
        integer                                     :: flag
        integer                                     :: status
        integer                                     :: FromBlock
        integer                                     :: control = 0, NColumns = 0.
        character(LEN = StringLength)               :: text, TimeSerieDataFile, AuxChar
        

        !----------------------------------------------------------------------

        call GetExtractType(FromBlock = FromBlock)

        allocate(PresentGauge)
        nullify (PresentGauge%Prev     )
        nullify (PresentGauge%Next     )
        nullify (PresentGauge%FirstWave)
        nullify (PresentGauge%LastWave )

        call GetData(PresentGauge%Name,                                                  &
                     Me%ObjEnterData, flag,                                              &
                     keyword    = 'NAME',                                                &
                     default    = 'Tidal Gauge X',                                       &
                     SearchType = FromBlock,                                             &
                     text       = 'Tide Gauge name missing',                             &
                     ClientModule ='ModuleGauge',                                        &
                     STAT       = status)
        if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR10'


        call GetData(PresentGauge%Longitude,                                             &
                     Me%ObjEnterData, flag,                                              &
                     keyword    = 'LONGITUDE',                                           &
                     SearchType = FromBlock,                                             &
                     text       = 'Tide Gauge Longitude missing',                        &
                     ClientModule ='ModuleGauge',                                        &
                     STAT       = status)
        if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR20'

        if ((flag > 0).and.(flag < 3)) then
            write (*,*)     '3 values for LONGITUDE were expected!'
            write (*,'(A)') trim(PresentGauge%name)
            stop 'ReadGaugeData - ModuleGauge - ERR30'
        endif
        if (flag == 3) control=control + 1

        call GetData(PresentGauge%Latitude,                                              &
                     Me%ObjEnterData, flag,                                              &
                     keyword    = 'LATITUDE',                                            &
                     SearchType = FromBlock,                                             &
                     text       = 'Tide Gauge Latitude missing',                         &
                     ClientModule ='ModuleGauge',                                        &
                     STAT       = status)
        if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR40'

if5 :   if ((flag > 0).and.(flag < 3)) then
            write (*,*)     '3 values for LATITUDE were expected!'
            write (*,'(A)') trim(PresentGauge%name)
            stop 'ReadGaugeData - ModuleGauge - ERR04'
        end if if5

        call Calc_Decimal_geo_coord(PresentGauge)

        call GetData(PresentGauge%Metric_X,                                              &
                     Me%ObjEnterData, flag,                                              &
                     keyword    = 'METRIC_X',                                            &
                     default    = FillValueReal,                                         &
                     SearchType = FromBlock,                                             &
                     ClientModule ='ModuleGauge',                                        &
                     STAT       = status)
        if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR50'


if1:    if (flag == 1) then 

            call GetData(PresentGauge%Metric_Y,                                          &
                         Me%ObjEnterData, flag,                                      &
                         keyword    = 'METRIC_Y',                                        &
                         default    = FillValueReal,                                     &
                         SearchType = FromBlock,                                         &
                         ClientModule ='ModuleGauge',                                    &
                         STAT       = status)
            if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR60'

            if (flag == 1) then
                PresentGauge%DefinedMetric = .true.
                Control                    = Control + 1
            endif

        endif if1

        call GetData(PresentGauge%GRID_I,                                                &
                     Me%ObjEnterData, flag,                                          &
                     keyword    = 'GRID_I',                                              &
                     default    = FillValueInt,                                          &
                     SearchType = FromBlock,                                             &
                     ClientModule ='ModuleGauge',                                        &
                     STAT       = status)
        if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR70'

if3:    if (flag == 1) then

            call GetData(PresentGauge%GRID_J,                                            &
                         Me%ObjEnterData, flag,                                      &
                         keyword    = 'GRID_J',                                          &
                         default    = FillValueInt,                                      &
                         SearchType = FromBlock,                                         &
                         ClientModule ='ModuleGauge',                                    &
                         STAT       = status)
            if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR80'


            if (flag == 1) then
                Control = Control + 1
            endif

        endif if3

            
if2 :   if (control == 0) then
            write (*,*) 'No Tide Gauge location provided'
            stop 'ReadGaugeData - ModuleGauge - ERR90'
        end if if2


        call GetData(PresentGauge%ReferenceLevel,                                        &
                     Me%ObjEnterData, flag,                                              &
                     keyword    = 'REF_LEVEL',                                           &
                     default    = 0.,                                                    &
                     SearchType = FromBlock,                                             &
                     text       = 'Warning. Tide Gauge Level missing. Zero assumed',     &
                     ClientModule ='ModuleGauge',                                        &
                     STAT       = status)
        if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR100'


        call GetData(PresentGauge%WaterLevel,                                            &
                     Me%ObjEnterData, flag,                                              &
                     keyword    = 'WATER_LEVEL',                                         &
                     default    = 0.,                                                    &
                     SearchType = FromBlock,                                             &
                     text       = text,                                                  &
                     ClientModule ='ModuleGauge',                                        &
                     STAT       = status)
        if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR110'

        call GetData(PresentGauge%TimeReference,                                         &
                     Me%ObjEnterData, flag,                                              &
                     keyword    = 'TIME_REF',                                            &
                     default    = 0.,                                                    &
                     SearchType = FromBlock,                                             &
                     text       = 'Warning. Tide Gauge TIME_REF missing',                &
                     ClientModule ='ModuleGauge',                                        &
                     STAT       = status)
        if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR120'

        if (flag == 1) PresentGauge%TimeReference = - PresentGauge%TimeReference !rcm


        call GetData(AuxChar ,                                                           &
                  Me%ObjEnterData, flag,                                                 &
                  keyword    = 'EVOLUTION',                                              &
                  default    = Char_Harmonics,                                           &
                  SearchType = FromBlock,                                                &
                  text       = 'Warning. Evolution missing. Harmonics assumed',          &
                  ClientModule ='ModuleGauge',                                           &
                  STAT       = status)
        if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR130'


        if      (AuxChar == Char_Harmonics) then
        
            PresentGauge%LevelEvolution = Harmonics

        else if (AuxChar == Char_TimeSerie) then

            PresentGauge%LevelEvolution = TimeSerie        

        else if (AuxChar == Char_Constant) then

            PresentGauge%LevelEvolution = Constant
        else 

            call SetError(FATAL_, INTERNAL_, "ReadGaugeData - ModuleGauge - ERR140")

        endif



        call GetData(AuxChar ,                                                           &
                  Me%ObjEnterData, flag,                                                 &
                  keyword    = 'EVOLUTION_VEL',                                          &
                  default    = Char_NoEvolution,                                         &
                  SearchType = FromBlock,                                                &
                  ClientModule ='ModuleGauge',                                           &
                  STAT       = status)
        if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR150'



        if      (AuxChar == Char_Harmonics) then
        
            PresentGauge%VelEvolution = Harmonics

        else if (AuxChar == Char_TimeSerie) then

            PresentGauge%VelEvolution = TimeSerie        

        else if (AuxChar == Char_NoEvolution) then

            PresentGauge%VelEvolution = NoEvolution        

        else 

            call SetError(FATAL_, INTERNAL_, "ReadGaugeData - ModuleGauge - ERR160")

        endif

        call GetData(AuxChar ,                                                           &
                  Me%ObjEnterData, flag,                                                 &
                  keyword    = 'REF_EVOLUTION',                                          &
                  default    = Char_Constant,                                            &
                  SearchType = FromBlock,                                                &
                  ClientModule ='ModuleGauge',                                           &
                  STAT       = status)
        if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR170'



        if      (AuxChar == Char_Constant) then
        
            PresentGauge%RefLevelEvolution = Constant

        else if (AuxChar == Char_TimeSerie) then

            PresentGauge%RefLevelEvolution = TimeSerie        

        else 

            call SetError(FATAL_, INTERNAL_, "ReadGaugeData - ModuleGauge - ERR180")

        endif


        call GetData(PresentGauge%RefLevelPhase ,                                       &
                Me%ObjEnterData, flag,                                                  &
                keyword    = 'REF_PHASE',                                               &
                default    = 0.,                                                        &
                SearchType = FromBlock,                                                 &
                ClientModule ='ModuleGauge',                                            &
                STAT       = status)
        if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR190'
        
        call GetData(PresentGauge%RefLevelAmpFactor,                                    &
                Me%ObjEnterData, flag,                                                  &
                keyword    = 'REF_AMP_FACTOR',                                          &
                default    = 1.,                                                        &
                SearchType = FromBlock,                                                 &
                ClientModule ='ModuleGauge',                                            &
                STAT       = status)
        if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR200'            

        if (PresentGauge%LevelEvolution    == TimeSerie .or.                            &
            PresentGauge%VelEvolution      == TimeSerie .or.                            &
            PresentGauge%RefLevelEvolution == TimeSerie) then

            call GetData(TimeSerieDataFile,                                             &
                      Me%ObjEnterData, flag,                                            &
                      keyword    = 'TIME_SERIE_FILE',                                   &
                      default    = Char_Harmonics,                                      &
                      SearchType = FromBlock,                                           &
                      ClientModule ='ModuleGauge',                                      &
                      STAT       = status)
            if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR210'

            PresentGauge%TimeSerieDataFile = TimeSerieDataFile

            call GetData(PresentGauge%CoveredColumn,                                    &
                      Me%ObjEnterData,  flag,                                           &
                      keyword    = 'COVERED_COLUMN',                                    &
                      default    = FillValueInt,                                        &
                      SearchType = FromBlock,                                           &
                      ClientModule ='ModuleGauge',                                      &
                      STAT       = status)
            if (status /= SUCCESS_) stop 'ReadGaugeData - ModuleGauge - ERR220'

            !flag = 1 The covered column is defined 
            if (flag == 1) then 
                NColumns = 1
                PresentGauge%NoCoveredColumn = .false.
            else
                NColumns = 0
                PresentGauge%NoCoveredColumn = .true.
            endif

            if (PresentGauge%LevelEvolution == TimeSerie) then

                call GetData(PresentGauge%LevelColumn,                                  &
                          Me%ObjEnterData, flag,                                        &
                          keyword    = 'LEVEL_COLUMN',                                  &
                          default    = FillValueInt,                                    &
                          SearchType = FromBlock,                                       &
                          text       = text,                                            &
                          ClientModule ='ModuleGauge',                                  &
                          STAT       = status)

                if (status /= SUCCESS_)                                                 &
                   call SetError(FATAL_, INTERNAL_, "ReadGaugeData - ModuleGauge - ERR230")

                NColumns = NColumns + 1

            endif

            if (PresentGauge%VelEvolution == TimeSerie) then

                call GetData(PresentGauge%VelUColumn,                                   &
                          Me%ObjEnterData, flag,                                        &
                          keyword    = 'VELU_COLUMN',                                   &
                          default    = FillValueInt,                                    &
                          SearchType = FromBlock,                                       &
                          text       = text,                                            &
                          ClientModule ='ModuleGauge',                                  &
                          STAT       = status)

                if (status /= SUCCESS_)                                                 &
                   call SetError(FATAL_, INTERNAL_, "ReadGaugeData - ModuleGauge - ERR240")


                call GetData(PresentGauge%VelVColumn,                                   &
                          Me%ObjEnterData, flag,                                        &
                          keyword    = 'VELV_COLUMN',                                   &
                          default    = FillValueInt,                                    &
                          SearchType = FromBlock,                                       &
                          text       = text,                                            &
                          ClientModule ='ModuleGauge',                                  &
                          STAT       = status)

                if (status /= SUCCESS_)                                                 &
                   call SetError(FATAL_, INTERNAL_, "ReadGaugeData - ModuleGauge - ERR250")


                NColumns = NColumns + 2

            endif

            if (PresentGauge%RefLevelEvolution == TimeSerie) then

                call GetData(PresentGauge%RefLevelColumn,                               &
                          Me%ObjEnterData, flag,                                        &
                          keyword    = 'REFLEVEL_COLUMN',                               &
                          default    = FillValueInt,                                    &
                          SearchType = FromBlock,                                       &
                          text       = text,                                            &
                          ClientModule ='ModuleGauge',                                  &
                          STAT       = status)

                if (status /= SUCCESS_)                                                 &
                   call SetError(FATAL_, INTERNAL_, "ReadGaugeData - ModuleGauge - ERR260")


                NColumns = NColumns + 1

            endif


            allocate(ColumnsToRead(NColumns))

            if (PresentGauge%NoCoveredColumn) then
                
                NColumns = 0

            else

                NColumns = 1

                ColumnsToRead(NColumns) = PresentGauge%CoveredColumn

            endif
            
            if (PresentGauge%LevelEvolution == TimeSerie) then
            
                NColumns = NColumns + 1

                ColumnsToRead(NColumns) = PresentGauge%LevelColumn

            endif 


            if (PresentGauge%VelEvolution == TimeSerie) then

                NColumns = NColumns + 1

                ColumnsToRead(NColumns) = PresentGauge%VelUColumn

                NColumns = NColumns + 1

                ColumnsToRead(NColumns) = PresentGauge%VelVColumn


            endif

            if (PresentGauge%RefLevelEvolution == TimeSerie) then

                NColumns = NColumns + 1

                ColumnsToRead(NColumns) = PresentGauge%RefLevelColumn

            endif


            call GetComputeTimeStep(Me%ObjTime, DT, STAT = status)
        
            if (status /= SUCCESS_)                                                     &
               call SetError(FATAL_, INTERNAL_, "ReadGaugeData - ModuleGauge - ERR280")


            call GetData(ReadTimeSerieDT,                                               &
                      Me%ObjEnterData,  flag,                                           &
                      keyword    = 'DT_SERIE',                                          &
                      default    = DT,                                                  &
                      SearchType = FromBlock,                                           &
                      text       = text,                                                &
                      ClientModule ='ModuleGauge',                                      &
                      STAT       = status)

            if (status /= SUCCESS_)                                                     &
               call SetError(FATAL_, INTERNAL_, "ReadGaugeData - ModuleGauge - ERR290")
          
        endif



        if (PresentGauge%LevelEvolution == Harmonics) then

            call ReadHarmonics(PresentGauge)

        endif 

        if (PresentGauge%LevelEvolution    == TimeSerie .or.                            &
            PresentGauge%RefLevelEvolution == TimeSerie .or.                            &
            PresentGauge%VelEvolution      == TimeSerie) then

            call StartTimeSerieInput(PresentGauge%ObjTimeSerie,                         &
                                     PresentGauge%TimeSerieDataFile,                    &
                                     Me%ObjTime,                                        &
                                     ColumnsToRead = ColumnsToRead,                     &
                                     DT_TimeSerie  = ReadTimeSerieDT,                   &
                                     CheckDates    = .true.,                            &
                                     ReadAll       = .true.,                            &
                                     STAT = status)

            if (status /= SUCCESS_)                                                     &
               call SetError(FATAL_, INTERNAL_, "ReadGaugeData - ModuleGauge - ERR300")


            deallocate(ColumnsToRead)

        endif


        !----------------------------------------------------------------------

    end subroutine ReadGaugeData

    !--------------------------------------------------------------------------

    subroutine ReadHarmonics(PresentGauge)

        !Arguments-------------------------------------------------------------
        type(T_TideGauge), pointer                  :: PresentGauge

        !Local-----------------------------------------------------------------
        integer                                     :: I, J

        character(LEN = WaveNameLength), dimension(NComponents), save :: WaveNameList
        logical, save                               :: start = .TRUE.

        !External--------------------------------------------------------------
        character(LEN = StringLength   ), parameter :: block_begin = '<beginharmonics>'
        character(LEN = StringLength   ), parameter :: block_end   = '<endharmonics>'
        character(LEN = WaveNameLength) :: WaveName
        

        real,dimension(2) :: AuxWave

        integer :: FromBlock
        integer :: STAT_CALL
        integer :: flag

        real :: Amplitude
        real :: Phase

        !Wave names------------------------------------------------------------

if1 :   if (start) then
            WaveNameList(  1) = 'Z0'
            WaveNameList(  2) = 'SA'
            WaveNameList(  3) = 'SSA'
            WaveNameList(  4) = 'MSM'
            WaveNameList(  5) = 'MM'
            WaveNameList(  6) = 'MSF'
            WaveNameList(  7) = 'MF'
            WaveNameList(  8) = 'ALP1'
            WaveNameList(  9) = '2Q1'
            WaveNameList( 10) = 'SIG1'
            WaveNameList( 11) = 'Q1'
            WaveNameList( 12) = 'RHO1'
            WaveNameList( 13) = 'O1'
            WaveNameList( 14) = 'TAU1'
            WaveNameList( 15) = 'BET1'
            WaveNameList( 16) = 'NO1'
            WaveNameList( 17) = 'CHI1'
            WaveNameList( 18) = 'PI1'
            WaveNameList( 19) = 'P1'
            WaveNameList( 20) = 'S1'
            WaveNameList( 21) = 'K1'
            WaveNameList( 22) = 'PSI1'
            WaveNameList( 23) = 'PHI1'
            WaveNameList( 24) = 'THE1'
            WaveNameList( 25) = 'J1'
            WaveNameList( 26) = 'OO1'
            WaveNameList( 27) = 'UPS1'
            WaveNameList( 28) = 'OQ2'
            WaveNameList( 29) = 'EPS2'
            WaveNameList( 30) = '2N2'
            WaveNameList( 31) = 'MU2'
            WaveNameList( 32) = 'N2'
            WaveNameList( 33) = 'NU2'
            WaveNameList( 34) = 'GAM2'
            WaveNameList( 35) = 'H1'
            WaveNameList( 36) = 'M2'
            WaveNameList( 37) = 'H2'
            WaveNameList( 38) = 'LDA2'
            WaveNameList( 39) = 'L2'
            WaveNameList( 40) = 'T2'
            WaveNameList( 41) = 'S2'
            WaveNameList( 42) = 'R2'
            WaveNameList( 43) = 'K2'
            WaveNameList( 44) = 'ETA2'
            WaveNameList( 45) = 'M3'
            WaveNameList( 46) = '2PO1'
            WaveNameList( 47) = 'SO1'
            WaveNameList( 48) = 'ST36'
            WaveNameList( 49) = '2NS2'
            WaveNameList( 50) = 'ST37'
            WaveNameList( 51) = 'ST1'
            WaveNameList( 52) = 'ST2'
            WaveNameList( 53) = 'ST3'
            WaveNameList( 54) = 'O2'
            WaveNameList( 55) = 'ST4'
            WaveNameList( 56) = 'SNK2'
            WaveNameList( 57) = 'OP2'
            WaveNameList( 58) = 'MKS2'
            WaveNameList( 59) = 'ST5'
            WaveNameList( 60) = 'ST6'
            WaveNameList( 61) = '2SK2'
            WaveNameList( 62) = 'MSN2'
            WaveNameList( 63) = 'ST7'
            WaveNameList( 64) = '2SM2'
            WaveNameList( 65) = 'ST38'
            WaveNameList( 66) = 'SKM2'
            WaveNameList( 67) = '2SN2'
            WaveNameList( 68) = 'NO3'
            WaveNameList( 69) = 'MO3'
            WaveNameList( 70) = 'NK3'
            WaveNameList( 71) = 'SO3'
            WaveNameList( 72) = 'MK3'
            WaveNameList( 73) = 'SP3'
            WaveNameList( 74) = 'SK3'
            WaveNameList( 75) = 'ST8'
            WaveNameList( 76) = 'N4'
            WaveNameList( 77) = '3MS4'
            WaveNameList( 78) = 'ST39'
            WaveNameList( 79) = 'MN4'
            WaveNameList( 80) = 'ST40'
            WaveNameList( 81) = 'ST9'
            WaveNameList( 82) = 'M4'
            WaveNameList( 83) = 'ST10'
            WaveNameList( 84) = 'SN4'
            WaveNameList( 85) = 'KN4'
            WaveNameList( 86) = 'MS4'
            WaveNameList( 87) = 'MK4'
            WaveNameList( 88) = 'SL4'
            WaveNameList( 89) = 'S4'
            WaveNameList( 90) = 'SK4'
            WaveNameList( 91) = 'MNO5'
            WaveNameList( 92) = '2MO5'
            WaveNameList( 93) = '3MP5'
            WaveNameList( 94) = 'MNK5'
            WaveNameList( 95) = '2MP5'
            WaveNameList( 96) = '2MK5'
            WaveNameList( 97) = 'MSK5'
            WaveNameList( 98) = '3KM5'
            WaveNameList( 99) = '2SK5'
            WaveNameList(100) = 'ST11'
            WaveNameList(101) = '2NM6'
            WaveNameList(102) = 'ST12'
            WaveNameList(103) = 'ST41'
            WaveNameList(104) = '2MN6'
            WaveNameList(105) = 'ST13'
            WaveNameList(106) = 'M6'
            WaveNameList(107) = 'MSN6'
            WaveNameList(108) = 'MKN6'
            WaveNameList(109) = '2MS6'
            WaveNameList(110) = '2MK6'
            WaveNameList(111) = 'NSK6'
            WaveNameList(112) = '2SM6'
            WaveNameList(113) = 'MSK6'
            WaveNameList(114) = 'ST42'
            WaveNameList(115) = 'S6'
            WaveNameList(116) = 'ST14'
            WaveNameList(117) = 'ST15'
            WaveNameList(118) = 'M7'
            WaveNameList(119) = 'ST16'
            WaveNameList(120) = '3MK7'
            WaveNameList(121) = 'ST17'
            WaveNameList(122) = 'ST18'
            WaveNameList(123) = '3MN8'
            WaveNameList(124) = 'ST19'
            WaveNameList(125) = 'M8'
            WaveNameList(126) = 'ST20'
            WaveNameList(127) = 'ST21'
            WaveNameList(128) = '3MS8'
            WaveNameList(129) = '3MK8'
            WaveNameList(130) = 'ST22'
            WaveNameList(131) = 'ST23'
            WaveNameList(132) = 'ST24'
            WaveNameList(133) = 'ST25'
            WaveNameList(134) = 'ST26'
            WaveNameList(135) = '4MK9'
            WaveNameList(136) = 'ST27'
            WaveNameList(137) = 'ST28'
            WaveNameList(138) = 'M10'
            WaveNameList(139) = 'ST29'
            WaveNameList(140) = 'ST30'
            WaveNameList(141) = 'ST31'
            WaveNameList(142) = 'ST32'
            WaveNameList(143) = 'ST33'
            WaveNameList(144) = 'M12'
            WaveNameList(145) = 'ST34'
            WaveNameList(146) = 'ST35'
        end if if1

        !----------------------------------------------------------------------

        call GetExtractType(FromBlock = FromBlock)

do1 :   do I = 1, NComponents
do2 :       do J = 1, WaveNameLength
                WaveName(J:J) = space
            end do do2

            WaveName = WaveNameList(I)
            WaveName = adjustl(WaveName)
            WaveName = trim   (WaveName)

            call GetData( AuxWave,                                                           &
                          Me%ObjEnterData,                                               &
                          flag,                                                              &
                          SearchType = FromBlock,                                            &
                          keyword=WaveName,                                                  &
                          ClientModule ='ModuleGauge',                                       &
                          STAT = STAT_CALL)
     
            if (STAT_CALL /= 0) &
                stop 'Error calling GetData, SUBROUTINE ReadHarmonics; Module ModuleGauge. ERR01.'

if2 :       if (flag .EQ. 2) then
                Amplitude = AuxWave(1)
                Phase     = AuxWave(2) / 360.0

                call NewWave(PresentGauge, WaveName, Amplitude, Phase)
            end if if2
        end do do1

if3 :   if (.NOT. associated(PresentGauge%FirstWave)) then
            write(*,*)
            write(*,*) 'Gauge Name: '
            write(*,'(A)') PresentGauge%Name
    
            stop 'Subroutine ReadHarmonics; Module ModuleGauge. ERR02'
        end if if3

        !----------------------------------------------------------------------

    end subroutine ReadHarmonics

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine NewWave(PresentGauge, WaveName, Amplitude, Phase)

        !Arguments------------------------------------------------------------

        type(T_TideGauge), pointer :: PresentGauge

        !External--------------------------------------------------------------

        integer                        :: STAT_CALL

        character(LEN = *), intent(IN) :: WaveName

        real,               intent(IN) :: Amplitude
        real,               intent(IN) :: Phase

        !Local-----------------------------------------------------------------

        type(T_TidalWave), pointer :: PresentWave

        !----------------------------------------------------------------------

        allocate(PresentWave     , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                           &
            stop 'SUBROUTINE NewWave; Module ModuleGauge. ERR10' 

        nullify (PresentWave%Prev)
        nullify (PresentWave%Next)


        PresentWave%Name = WaveName
        PresentWave%Name = adjustl(PresentWave%Name)
        PresentWave%Name = trim   (PresentWave%Name)


        PresentWave%Amplitude = Amplitude
        PresentWave%Phase     = Phase


if2 :   if (.NOT. associated(PresentGauge%FirstWave)) then      !If list is not initialized
            PresentGauge%FirstWave => PresentWave               !FirstGauge element is always FirstGauge
            PresentGauge%LastWave  => PresentWave               !FirstGauge element is also the LastGauge one

            nullify(PresentGauge%FirstWave%Prev)
            nullify(PresentGauge%LastWave%Next )
        else 
            PresentWave%Prev           => PresentGauge%LastWave

            PresentGauge%LastWave%Next => PresentWave 
            PresentGauge%LastWave      => PresentWave

            nullify(PresentGauge%LastWave%Next)
        end if if2

        !----------------------------------------------------------------------

    end subroutine NewWave

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine ReadOldForm(BlockFound, PresentGauge)

        !Arguments------------------------------------------------------------
        type(T_TideGauge), pointer                  :: PresentGauge
        Logical, intent(OUT)                        :: BlockFound

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        real                                        :: Amplitude
        real                                        :: Phase
        character(LEN = WaveNameLength)             :: WaveName
        integer                                     :: NMaregMar, NMaregRio, NWaves
        integer, dimension(3)                       :: Longitude
        integer, dimension(3)                       :: Latitude
        integer                                     :: Grid_I
        integer                                     :: Grid_J
        real                                        :: ReferenceLevel
        real                                        :: TimeReference
        character(LEN = StringLength   )            :: Coment1, Coment2 
        character(LEN = StringLength   )            :: String

        !----------------------------------------------------------------------

        BlockFound = .False. 

        call UnitsManager(Me%UnitGauge, OPEN_FILE, STAT = STAT_CALL) 

        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'SUBROUTINE ReadOldForm; Module ModuleGauge. ERR01' 

        open(Unit = Me%UnitGauge, File = Me%FileName, form = 'FORMATTED', status = 'OLD',         &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'SUBROUTINE ReadOldForm; Module ModuleGauge. ERR02' 

       

        rewind(Me%UnitGauge)
        read  (Me%UnitGauge,'(A)',end=100) Coment1
        read  (Me%UnitGauge,'(A)',end=100) Coment2
        read  (Me%UnitGauge,*,    end=100) NMaregMar,NMaregRio
        read  (Me%UnitGauge,*,    end=100) Latitude
        read  (Me%UnitGauge,*,    end=100) Longitude
        read  (Me%UnitGauge,*,    end=100) Grid_I, Grid_J
        read  (Me%UnitGauge,*,    end=100) TimeReference
        read  (Me%UnitGauge,*,    end=100) ReferenceLevel
        read  (Me%UnitGauge,*,    end=100) NWaves

        call Calc_Decimal_geo_coord(PresentGauge)


        allocate(PresentGauge          , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'SUBROUTINE ReadOldForm; Module ModuleGauge. ERR03' 

        nullify (PresentGauge%Prev     )
        nullify (PresentGauge%Next     )
        nullify (PresentGauge%FirstWave)
        nullify (PresentGauge%LastWave )

        PresentGauge%Latitude       = Latitude
        PresentGauge%Longitude      = Longitude
        PresentGauge%Grid_I         = Grid_I
        PresentGauge%Grid_J         = Grid_J
        PresentGauge%TimeReference  = TimeReference
        PresentGauge%ReferenceLevel = ReferenceLevel


        
do1 :   do 
            read(Me%UnitGauge,*, end=10) WaveName, String

            read(String,*) Amplitude, Phase
            Phase = Phase / 360.0

            call NewWave(PresentGauge, WaveName, Amplitude, Phase)
        end do do1

10      continue


        BlockFound = .True. 
        Write (*,60)

        return

100     continue

        call UnitsManager(Me%UnitGauge, CLOSE_FILE, STAT = STAT_CALL) 

        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'SUBROUTINE ReadOldForm; Module ModuleGauge. ERR04' 



        !Format----------------------------------------------------------------

60      format(///T10,'A T T E N T I O N',//T10,                       & 
              'The TIDE components file is written in an OLD Format.', &
              //T10,'It is advisable that you change it to the NEW Format.',////)

        !----------------------------------------------------------------------

    end subroutine ReadOldForm

    !
    ! This routine verifies if the position GridI and GridJ are consistent
    ! with the MetricX and MetrixY coordinates (Latitude and Longitude should be
    ! added later)
    ! This routine should be modified in a way that the user can especify directly
    ! the MetricX and MetricY coordinates, rather than especifing the grid 
    ! coordinates, or the latitude and longitude.
    ! Also, its necessary to verify if the position of the gauge is really in the 
    ! inferior left corner, or if it will be located in the center of the cell.
    !
    subroutine VerifyGaugeLocation(Gauge, HorizontalGridID)

        !Arguments-------------------------------------------------------------
        type (T_TideGauge), pointer         :: Gauge
        integer                             :: HorizontalGridID

        !Local-----------------------------------------------------------------
        real, dimension(:   ), pointer      :: XX, YY
        real                                :: XORIG, YORIG, GridRotation
        integer                             :: status

        !Begin-----------------------------------------------------------------
        if (.not. Gauge%DefinedMetric) then

            !Gets the Grid angle
            call GetGridAngle(HorizontalGridID, GridRotation, STAT = status)
            if (status /= SUCCESS_) stop 'VerifyGaugeLocation - ModuleGauge - ERR01'

            !Gets Origin of the Bathymetry
            call GetGridOrigin(HorizontalGridID, Xorig, Yorig, STAT = status)
            if (status /= SUCCESS_) stop 'VerifyGaugeLocation - ModuleGauge - ERR02'

            !Gets XX and YY
            call GetHorizontalGrid(HorizontalGridID, XX = XX, YY = YY, STAT = status)
            if (status /= SUCCESS_) stop 'VerifyGaugeLocation - ModuleGauge - ERR03'

            !Verifies Position in X
            Gauge%Metric_X    = (XX(Gauge%Grid_J) + XX(Gauge%Grid_J + 1)) / 2.

            !Verifies Position in Y
            Gauge%Metric_Y    = (YY(Gauge%Grid_I) + YY(Gauge%Grid_I + 1)) / 2.

            !Localized the cell center using the bathymetry coordinates 
            call RODAXY (XORIG, YORIG, GridRotation, Gauge%Metric_X, Gauge%Metric_Y)

            !Ungets XX and YY
            call UnGetHorizontalGrid(HorizontalGridID, XX, stat = status)
            if (status /= SUCCESS_) stop 'VerifyGaugeLocation - ModuleGauge - ERR04'
            call UnGetHorizontalGrid(HorizontalGridID, YY, stat = status)
            if (status /= SUCCESS_) stop 'VerifyGaugeLocation - ModuleGauge - ERR05'

        endif

    end subroutine VerifyGaugeLocation


    !--------------------------------------------------------------------------

    subroutine LocationConsistence

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_TideGauge), pointer :: PresentGauge
        logical                    :: DefinedMetric

        !----------------------------------------------------------------------

        !----------------------------------------------------------------------

        PresentGauge => Me%FirstGauge

        DefinedMetric = PresentGauge%DefinedMetric

do2 :   do while (associated(PresentGauge)) 

            if (.not. (DefinedMetric .EQV. PresentGauge%DefinedMetric)) then

                stop 'LocationConsistence - ModuleGauge - ERR01'

            endif
                
            PresentGauge => PresentGauge%Next

        end do do2

        !----------------------------------------------------------------------

    end subroutine LocationConsistence

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetNumberWaves(PresentGauge, NWaves)

        !Arguments-------------------------------------------------------------
        type(T_TideGauge), pointer                  :: PresentGauge 
        integer, intent(OUT)                        :: NWaves

        !Local-----------------------------------------------------------------
        type(T_TidalWave), pointer                  :: PresentWave

        !----------------------------------------------------------------------

        NWaves = 0
        PresentWave => PresentGauge%FirstWave

do2 :   do while (associated(PresentWave)) 
            NWaves = NWaves + 1
            PresentWave => PresentWave%Next
        end do do2

        !----------------------------------------------------------------------

    end subroutine GetNumberWaves

    !--------------------------------------------------------------------------

    subroutine GetNGauges(GaugeID, NGauges, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GaugeID
        integer,           intent(OUT)              :: NGauges
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        type(T_TideGauge), pointer                  :: PresentGauge

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            NGauges = 0

            PresentGauge => Me%FirstGauge

do2 :       do while (associated(PresentGauge)) 
                NGauges = NGauges + 1
                    
                PresentGauge => PresentGauge%Next
            end do do2

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetNGauges

    !--------------------------------------------------------------------------

    subroutine GetNComponents(GaugeID, NumberComp, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GaugeID
        integer,           intent(OUT)              :: NumberComp
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
 
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            NumberComp = NComponents

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetNComponents

    !--------------------------------------------------------------------------

    subroutine GetMetricGaugeLocation(GaugeID, XLocation, YLocation, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: GaugeID
        real, dimension(:), pointer     :: XLocation, YLocation
        integer, optional, intent(OUT)  :: STAT

        !Local-----------------------------------------------------------------
        integer                         :: STAT_, ready_
        type(T_TideGauge), pointer      :: PresentGauge
        integer                         :: i

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            i = 0
            PresentGauge => Me%FirstGauge
            do while (associated(PresentGauge)) 
                i = i + 1
            
                XLocation(i) = PresentGauge%Metric_X
                YLocation(i) = PresentGauge%Metric_Y

                PresentGauge => PresentGauge%Next
             enddo

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetMetricGaugeLocation


    !--------------------------------------------------------------------------

    subroutine GetGaugeLocation(GaugeID, XLocation, YLocation, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: GaugeID
        real, dimension(:), pointer     :: XLocation, YLocation
        integer, optional, intent(OUT)  :: STAT

        !Local-----------------------------------------------------------------
        integer                         :: STAT_, ready_
        type(T_TideGauge), pointer      :: PresentGauge
        integer                         :: i

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            i = 0
            PresentGauge => Me%FirstGauge
            do while (associated(PresentGauge)) 
                i = i + 1
            
                XLocation(i) = PresentGauge%DecimalLongitude
                YLocation(i) = PresentGauge%DecimalLatitude

                PresentGauge => PresentGauge%Next
             enddo

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetGaugeLocation
    
    !--------------------------------------------------------------------------

    subroutine GetIJReferenceLevel(GaugeID, HorizontalGridID, I, J, ReferenceLevel, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: GaugeID
        integer                         :: HorizontalGridID
        integer,           intent(IN )  :: I, J
        real,              intent(OUT)  :: ReferenceLevel
        integer, optional, intent(OUT)  :: STAT

        !Local-----------------------------------------------------------------
        real                            :: PX, PY             
        integer                         :: STAT_, ready_
        type(T_TideGauge), pointer      :: PresentGauge
        logical                         :: Found
        real                            :: PercI, PercJ 
        integer                         :: Iaux, Jaux        
        

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then



            PresentGauge => Me%FirstGauge

            Found = .false.

            do while (associated(PresentGauge) .and. .not. Found) 

                if (PresentGauge%Grid_I == I .and. PresentGauge%Grid_J == J) then

                    ReferenceLevel = PresentGauge%ReferenceLevel                    

                    Found = .true.

                endif
                
                PX = PresentGauge%Metric_X
                PY = PresentGauge%Metric_Y
                
                call GetXYCellZ(HorizontalGridID, PX, PY, Iaux, Jaux, STAT = STAT_)

  
                if (I == Iaux .and. J == Jaux) then
                
                    call GetXYCellZ(HorizontalGridID, PX, PY, Iaux, Jaux, PercI = PercI, PercJ = PercJ, STAT = STAT_)

                    if (abs(PercI-0.5) < 0.02 .and. abs(PercJ-0.5) < 0.02) then

                        ReferenceLevel = PresentGauge%ReferenceLevel                     

                        Found = .true.
                        
                    endif
                
                endif

                PresentGauge => PresentGauge%Next

             enddo


            if (Found) then

                STAT_ = SUCCESS_

            else

                STAT_ = NOT_FOUND_ERR_

            endif
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetIJReferenceLevel

    !----------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine GetIJWaterLevel(GaugeID, HorizontalGridID, I, J, WaterLevel, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: GaugeID
        integer                         :: HorizontalGridID
        integer,           intent(IN )  :: I, J
        real,              intent(OUT)  :: WaterLevel
        integer, optional, intent(OUT)  :: STAT

        !Local-----------------------------------------------------------------
        real                            :: PX, PY       
        integer                         :: STAT_, ready_
        type(T_TideGauge), pointer      :: PresentGauge
        logical                         :: Found
        real                            :: PercI, PercJ 
        integer                         :: Iaux, Jaux 

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then



            PresentGauge => Me%FirstGauge

            Found = .false.

            do while (associated(PresentGauge) .and. .not. Found) 

                if (PresentGauge%Grid_I == I .and. PresentGauge%Grid_J == J) then

                    WaterLevel = PresentGauge%WaterLevel                    

                    Found      = .true.

                endif
                
                PX = PresentGauge%Metric_X
                PY = PresentGauge%Metric_Y                
                
                call GetXYCellZ(HorizontalGridID, PX, PY, Iaux, Jaux, STAT = STAT_)

  
                if (I == Iaux .and. J == Jaux) then

                    call GetXYCellZ(HorizontalGridID, PX, PY, Iaux, Jaux, PercI = PercI, PercJ = PercJ, STAT = STAT_)

                    if (abs(PercI-0.5) < 0.02 .and. abs(PercJ-0.5) < 0.02) then

                        WaterLevel = PresentGauge%WaterLevel                    

                        Found = .true.
                        
                    endif

                endif

                PresentGauge => PresentGauge%Next

             enddo

            if (Found) then

                STAT_ = SUCCESS_

            else

                STAT_ = NOT_FOUND_ERR_

            endif
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetIJWaterLevel

    !----------------------------------------------------------------------


    subroutine GetWaterLevel(GaugeID, WaterLevel, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: GaugeID
        real, dimension(:), pointer     :: WaterLevel
        integer, optional, intent(OUT)  :: STAT

        !Local-----------------------------------------------------------------
        integer                         :: STAT_, ready_
        type(T_TideGauge), pointer      :: PresentGauge
        integer                         :: i

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            i = 0
            PresentGauge => Me%FirstGauge
            do while (associated(PresentGauge)) 
                i = i + 1
            
                WaterLevel(i) = PresentGauge%WaterLevel

                PresentGauge => PresentGauge%Next
             enddo

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetWaterLevel

    !----------------------------------------------------------------------

    !----------------------------------------------------------------------


    subroutine GetReferenceLevel(GaugeID, ReferenceLevel, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: GaugeID
        real, dimension(:), pointer     :: ReferenceLevel
        integer, optional, intent(OUT)  :: STAT

        !Local-----------------------------------------------------------------
        integer                         :: STAT_, ready_
        type(T_TideGauge), pointer      :: PresentGauge
        integer                         :: i

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            i = 0
            PresentGauge => Me%FirstGauge
            do while (associated(PresentGauge)) 
                i = i + 1
            
                ReferenceLevel(i) = PresentGauge%ReferenceLevel

                PresentGauge => PresentGauge%Next
             enddo

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetReferenceLevel

    !----------------------------------------------------------------------

    subroutine GetLevelEvolution(GaugeID, LevelEvolution, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: GaugeID
        integer, dimension(:), pointer  :: LevelEvolution
        integer, optional, intent(OUT)  :: STAT

        !Local-----------------------------------------------------------------
        integer                         :: STAT_, ready_
        type(T_TideGauge), pointer      :: PresentGauge
        integer                         :: i

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            i = 0
            PresentGauge => Me%FirstGauge
            do while (associated(PresentGauge)) 
                i = i + 1
            
                LevelEvolution(i) =  PresentGauge%LevelEvolution

                PresentGauge      => PresentGauge%Next
             enddo

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetLevelEvolution

    !--------------------------------------------------------------------------

    subroutine GetVelEvolution(GaugeID, VelEvolution, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: GaugeID 
        integer, dimension(:), pointer  :: VelEvolution
        integer, optional, intent(OUT)  :: STAT

        !Local-----------------------------------------------------------------
        integer                         :: STAT_, ready_
        type(T_TideGauge), pointer      :: PresentGauge
        integer                         :: i

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            i = 0
            PresentGauge => Me%FirstGauge
            do while (associated(PresentGauge)) 
                i = i + 1
            
                VelEvolution(i) =  PresentGauge%VelEvolution

                PresentGauge    => PresentGauge%Next
             enddo

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetVelEvolution

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine GetTriangGaugesON(GaugeID, TriangulationON, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GaugeID
        logical,           intent(OUT)              :: TriangulationON
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
 
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            TriangulationON = Me%Triangulation

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetTriangGaugesON

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GaugeLevel(GaugeID, WaterLevel, OpenPoints, Time,                        &
                          ReferenceLevel, VelocityU, VelocityV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: GaugeID 
        type(T_Time), intent(IN)                    :: Time
        real   , pointer, dimension(:)              :: WaterLevel
        real   , pointer, dimension(:)              :: OpenPoints
        real   , pointer, optional, dimension(:)    :: VelocityU, VelocityV, ReferenceLevel
        integer, optional, intent(OUT  )            :: STAT

        !Local-----------------------------------------------------------------
        type(T_Time)                                :: TimeAux        
        integer                                     :: ready_             
        integer                                     :: STAT_
        integer                                     :: I
        type(T_TideGauge), pointer                  :: PresentGauge

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            I = 0
            PresentGauge => Me%FirstGauge
do1 :       do while (associated(PresentGauge))
                I = I + 1

if2 :           if      (PresentGauge%LevelEvolution == Harmonics) then 
                    
                    call HarmonicsWaterLevel(PresentGauge, WaterLevel(I), Time)


                else if (PresentGauge%LevelEvolution == TimeSerie) then if2
                        
                    call TimeSerieValue(PresentGauge, PresentGauge%WaterLevel, PresentGauge%LevelColumn,     Time)
                    
                    WaterLevel(I) = PresentGauge%WaterLevel

                else if (PresentGauge%LevelEvolution == Constant)  then if2
                
                    WaterLevel(I) = PresentGauge%WaterLevel

                endif if2

                if (present(ReferenceLevel)) then

if3:                if (PresentGauge%RefLevelEvolution == TimeSerie) then 

                        TimeAux = Time + PresentGauge%RefLevelPhase
                        
                        call TimeSerieValue(PresentGauge, PresentGauge%ReferenceLevel, PresentGauge%RefLevelColumn, TimeAux)

                        ReferenceLevel (I) = PresentGauge%ReferenceLevel

                    else if (PresentGauge%RefLevelEvolution == Constant) then 
                
                        ReferenceLevel (I) = PresentGauge%ReferenceLevel
                    
                    endif if3
                    
                    ReferenceLevel (I) = ReferenceLevel (I) * PresentGauge%RefLevelAmpFactor

                endif

                    
if4:            if (PresentGauge%LevelEvolution    == TimeSerie   .or.                   &
                    PresentGauge%VelEvolution      == TimeSerie   .or.                   &
                    PresentGauge%RefLevelEvolution == TimeSerie       ) then 

                    if (PresentGauge%NoCoveredColumn) then
                    
                        OpenPoints(I) = 1
                    else

                        call TimeSerieValue(PresentGauge, OpenPoints(I), PresentGauge%CoveredColumn, Time)

                        if (OpenPoints(I) < 1) OpenPoints(I) = 0

                    endif

                else

                    OpenPoints(I) = 1

                endif if4


if5:            if (PresentGauge%VelEvolution == TimeSerie)  then

                    if (present(VelocityU))                                              &
                        call TimeSerieValue(PresentGauge, VelocityU(I), PresentGauge%VelUColumn, Time)

                    if (present(VelocityV))                                              &
                        call TimeSerieValue(PresentGauge, VelocityV(I), PresentGauge%VelVColumn, Time)

                else if (PresentGauge%VelEvolution == NoEvolution)  then

                    if (present(VelocityU)) VelocityU(I) = 0.
                
                    if (present(VelocityV)) VelocityV(I) = 0.

                endif if5

                PresentGauge => PresentGauge%Next

            end do do1
            

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GaugeLevel

    !--------------------------------------------------------------------------

    subroutine HarmonicsWaterLevel(PresentGauge, WaterLevel, TimeIn)

        !Arguments-------------------------------------------------------------
        type(T_TideGauge), pointer          :: PresentGauge
        type(T_Time), intent(IN)            :: TimeIn
        real                                :: WaterLevel

        !Local-----------------------------------------------------------------

        integer :: status
        integer :: NWaves

        real, pointer, dimension(:) :: WaveAmplitude
        real, pointer, dimension(:) :: WavePhase
        

        integer :: J

        type(T_TidalWave), pointer :: PresentWave

        real                       :: WaterLevelRef
        
        real, dimension(1:115)     :: WaveAmplitude_, WavePhase_
        character(LEN = WaveNameLength), dimension(1:115) :: WaveName_

        
        !----------------------------------------------------------------------                         


        !----------------------------------------------------------------------

        call GetNumberWaves(PresentGauge, NWaves)

        allocate(WaveAmplitude(NWaves), STAT = status)
        if (status /= SUCCESS_)                                                          &
            call SetError(FATAL_, INTERNAL_, "HarmonicsWaterLevel - ModuleGauge - ERR10")

        allocate(WavePhase    (NWaves), STAT = status)

        if (status /= SUCCESS_)                                                          &
            call SetError(FATAL_, INTERNAL_, "HarmonicsWaterLevel - ModuleGauge - ERR20")

        J = 0
        PresentWave => PresentGauge%FirstWave
do5 :   do while (associated(PresentWave))
            J = J + 1

            WaveAmplitude_(J) = PresentWave%Amplitude
            WavePhase_    (J) = PresentWave%Phase
            WaveName_     (J) = PresentWave%Name

            WaveAmplitude (J) = PresentWave%Amplitude
            WavePhase     (J) = PresentWave%Phase

            PresentWave => PresentWave%Next
        end do do5

        ! Reference level is add in the module open boundary
        ! The water level is the Gauge model is compute with relation to 
        ! the hydrographic zero
        WaterLevelRef = 0.
        
        if (Me%TidePrevision == Toga_) then

            call TogaLevel(PresentGauge%ObjToga,                                        &  
                           PresentGauge%WaterLevel,                                     &
                           PresentGauge%DecimalLatitude,                                &
                           PresentGauge%TimeReference,                                  &
                           WaterLevelRef,                                               &
                           NWaves, WaveAmplitude, WavePhase,                            &
                           TimeIn,                                                      &
                           STAT = status)

            if (status /= SUCCESS_)                                                         &
                call SetError(FATAL_, INTERNAL_, "HarmonicsWaterLevel - ModuleGauge - ERR30")
                           
        
        else if (Me%TidePrevision == Task2000_) then
        
            call Task2000Level(PresentGauge%WaterLevel,                                 &
                               PresentGauge%TimeReference,                              &
                               NWaves, WaveAmplitude_,                                  &
                               WavePhase_, WaveName_, TimeIn, STAT = status)

            if (status /= SUCCESS_)                                                         &
                call SetError(FATAL_, INTERNAL_, "HarmonicsWaterLevel - ModuleGauge - ERR40")

        endif
        
        WaterLevel = PresentGauge%WaterLevel

        deallocate(WaveAmplitude, STAT = status)

        if (status /= SUCCESS_)                                                         &
            call SetError(FATAL_, INTERNAL_, "HarmonicsWaterLevel - ModuleGauge - ERR50")

        deallocate(WavePhase    , STAT = status)

        if (status /= SUCCESS_)                                                         &
            call SetError(FATAL_, INTERNAL_, "HarmonicsWaterLevel - ModuleGauge - ERR60")

        nullify   (WaveAmplitude)
        nullify   (WavePhase    )

    end subroutine HarmonicsWaterLevel

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    subroutine TimeSerieValue(PresentGauge, Value, Column, Time)

        !Arguments-------------------------------------------------------------
        type(T_TideGauge), pointer          :: PresentGauge
        type(T_Time)                        :: Time
        real                                :: Value
        integer                             :: Column

        !External--------------------------------------------------------------

        integer :: status

        !Local-----------------------------------------------------------------
        type(T_Time)                        :: Time1, Time2
        real                                :: Value1, Value2
        logical                             :: TimeCycle

        !----------------------------------------------------------------------


        call GetTimeSerieValue(PresentGauge%ObjTimeSerie, Time, Column,                 &
                     Time1, Value1, Time2, Value2, TimeCycle, STAT = status) 
       
        if (Time < Time1)       stop "TimeSerieValue - ModuleGauge - ERR10"
               
        Value = ((Time - Time1) * Value2 + (Time2 - Time) * Value1) / (Time2 - Time1)

        if (status /= SUCCESS_) stop "TimeSerieValue - ModuleGauge - ERR20"

        if (TimeCycle)          stop "TimeSerieValue - ModuleGauge - ERR30"

    end subroutine TimeSerieValue


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillGauge(GaugeID, STAT)  

        !Arguments-------------------------------------------------------------
        integer                                     :: GaugeID 
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: nUsers
        integer                                     :: STAT_, ready_ 

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(GaugeID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mGAUGE_,  Me%InstanceID)

            if (nUsers == 0) then

                !Deassociates External Instances
                nUsers = DeassociateInstance (mTIME_,     Me%ObjTime)
                if (nUsers == 0) stop 'KillGauge - Gauge - ERR01'

                if (associated(Me%FirstGauge)) call KillGaugeList(Me%FirstGauge, Me%LastGauge)
                
                call DeallocateInstance

                GaugeID = 0
                STAT_   = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine KillGauge

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Gauge), pointer          :: AuxObjGauge
        type (T_Gauge), pointer          :: PreviousObjGauge

        !Updates pointers
        if (Me%InstanceID == FirstGauge%InstanceID) then
            FirstGauge => FirstGauge%Next
        else
            PreviousObjGauge => FirstGauge
            AuxObjGauge      => FirstGauge%Next
            do while (AuxObjGauge%InstanceID /= Me%InstanceID)
                PreviousObjGauge => AuxObjGauge
                AuxObjGauge      => AuxObjGauge%Next
            enddo

            !Now update linked list
            PreviousObjGauge%Next => AuxObjGauge%Next

        endif
            
        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------

    subroutine KillTidalWaveList(FirstWave, LastWave)

        !Arguments------------------------------------------------------------

        type(T_TidalWave), pointer :: FirstWave
        type(T_TidalWave), pointer :: LastWave

        !External--------------------------------------------------------------
                       
        !----------------------------------------------------------------------

do1 :   do while (associated(LastWave%Prev))  !Deallocates all except FirstWave.
            LastWave => LastWave%Prev
            deallocate(LastWave%Next)
        end do do1   


        !Deallocates FirstWave.
if1 :   if (associated(FirstWave)) then
            deallocate(FirstWave)
            nullify   (LastWave )
        end if if1

        !----------------------------------------------------------------------

    end subroutine KillTidalWaveList

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine KillGaugeList(FirstGauge, LastGauge)

        !Arguments------------------------------------------------------------

        type(T_TideGauge), pointer :: FirstGauge 
        type(T_TideGauge), pointer :: LastGauge 

        !External--------------------------------------------------------------
                       
        !----------------------------------------------------------------------

do1 :   do while (associated(LastGauge%Prev))  !Deallocates all except FirstWave.
            LastGauge => LastGauge%Prev

            if (associated(LastGauge%Next%FirstWave)) &
                call KillTidalWaveList(LastGauge%Next%FirstWave, LastGauge%Next%LastWave)



            if (LastGauge%Next%ObjToga      /= 0) call KillToga(LastGauge%Next%ObjToga)
            if (LastGauge%Next%ObjToga      /= 0) call KillToga(LastGauge%Next%ObjToga)
            if (LastGauge%Next%ObjTimeSerie /= 0) call KillTimeSerie(LastGauge%Next%ObjTimeSerie)


            deallocate(LastGauge%Next)

        end do do1   


        !Deallocates FirstGauge.
if1 :   if (associated(FirstGauge)) then
            if (associated(FirstGauge%FirstWave)) &
                call KillTidalWaveList(FirstGauge%FirstWave, FirstGauge%LastWave)


            if (FirstGauge%ObjToga      /= 0) call KillToga(FirstGauge%ObjToga)
            if (FirstGauge%ObjTimeSerie /= 0) call KillTimeSerie(FirstGauge%ObjTimeSerie)

            deallocate(FirstGauge)
            nullify   (LastGauge )
        end if if1

        !----------------------------------------------------------------------

    end subroutine KillGaugeList

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Ready (ObjGauge_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjGauge_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjGauge_ID > 0) then
            call LocateObjGauge(ObjGauge_ID)
            ready_ = VerifyReadLock (mGAUGE_,  Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjGauge (ObjGaugeID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjGaugeID

        !Local-----------------------------------------------------------------

        Me => FirstGauge
        do while (associated (Me))
            if (Me%InstanceID == ObjGaugeID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleGauge - LocateObjGauge - ERR01'

    end subroutine LocateObjGauge

    !--------------------------------------------------------------------------

end module ModuleGauge

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
