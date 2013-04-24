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
                                       GetXYCellZ, GetXYCellZ_ThreadSafe, UnGetHorizontalGrid
    use ModuleTimeSerie,        only : StartTimeSerieInput, GetTimeSerieValue, KillTimeSerie
    use ModuleFunctions,        only : RodaXY
    use ModuleTask2000,         only : Task2000Level, NTask2000          

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
    public  :: GetIJWaterLevel_ThreadSafe
    public  :: GetIJReferenceLevel
    public  :: GetIJReferenceLevel_ThreadSafe
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
    integer, parameter :: NComponents    = 146
    integer, parameter :: NAdmit         = 19

    !Evolution Options
    character(LEN = StringLength), parameter :: Char_Harmonics   = 'Harmonics',          &
                                                Char_TimeSerie   = 'Time Serie',         &
                                                Char_NoEvolution = 'No',                 &
                                                Char_Constant    = 'Constant'
    integer, parameter                       :: Constant = 0, Harmonics = 1,             &
                                                TimeSerie = 2, NoEvolution = 3

    !PREDICTION method 
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
        
        integer                             :: TidePREDICTION
        
        logical                             :: ComputeAdmittance
       
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
            
            call GetData(Me%Triangulation,                                              &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType = FromFile,                                         &
                         keyword    = 'TRIANGULATION',                                  &
                         Default    = .true.,                                           &
                         ClientModule ='ModuleGauge',                                   &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauges - ERR30' 
                
            call GetData(Me%TidePREDICTION,                                             &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType = FromFile,                                         &
                         keyword    = 'PREVISION',                                      &
                         Default    = Task2000_,                                        &
                         ClientModule ='ModuleGauge',                                   &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauges - ERR40' 

            if (iflag == 0) then
                call GetData(Me%TidePREDICTION,                                         &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType = FromFile,                                     &
                             keyword    = 'PREDICTION',                                 &
                             Default    = Task2000_,                                    &
                             ClientModule ='ModuleGauge',                               &
                             STAT       = STAT_CALL)            
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauges - ERR50'             
            endif
            
            if (Me%TidePREDICTION /= Task2000_ .and. Me%TidePREDICTION /= Toga_) then
                write(*,*) 'The tide PREDICTION method is not known'
                stop       'ConstructGauges - ModuleGauges - ERR38'
            endif

            call GetData(Me%ComputeAdmittance,                                          &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType = FromFile,                                         &
                         keyword    = 'ADMITTANCE',                                     &
                         Default    = .false.,                                          &
                         ClientModule ='ModuleGauge',                                   &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauges - ERR50' 

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
                    stop 'ConstructGauges - ModuleGauges - ERR60' 
                end if if4
            end do do2


            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauges - ERR70' 


            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauges - ERR80' 

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
                            stop 'ConstructGauges - ModuleGauges - ERR90' 
                        end if if7
                    else
                        exit do3    !No more blocks
                    end if if6
                end do do3
            end if if3



if8 :       if (.NOT. associated(Me%FirstGauge)) then
                stop 'ConstructGauges - ModuleGauges - ERR100' 
            else
                PresentGauge => Me%FirstGauge
do4 :           do while (associated(PresentGauge))
                    call GetNumberWaves(PresentGauge, NWaves)

                    allocate(WaveName(NWaves), STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauge - ERR110.'


                    I = 0
                    PresentWave => PresentGauge%FirstWave
do5 :               do while (associated(PresentWave))
                        I = I + 1

                        WaveName(I) = PresentWave%Name

                        PresentWave => PresentWave%Next
                    end do do5

                    call ConstructToga(PresentGauge%ObjToga, NWaves, WaveName, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructGauges - ModuleGauge - ERR120'

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
            
            stop 'ModuleGauge - ConstructGauges - ERR130' 

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

                if (AuxWave(1) < 0.) then
                    write(*,*) 'Tidal component ', trim(WaveName),' in the tidal gauge named ', trim(PresentGauge%Name)
                    write(*,*) 'has an amplitude below ZERO'
                    stop
                endif
                if (AuxWave(1) > 20.) then
                    write(*,*) 'Tidal component ', trim(WaveName),' in the tidal gauge named ', trim(PresentGauge%Name)
                    write(*,*) 'has an amplitude above 20'
                    stop
                endif

                if (AuxWave(2) < -360.) then
                    write(*,*) 'Tidal component ', trim(WaveName),' in the tidal gauge named ', trim(PresentGauge%Name)
                    write(*,*) 'has a phase below -360'
                    stop
                endif
                if (AuxWave(2) >  360.) then
                    write(*,*) 'Tidal component ', trim(WaveName),' in the tidal gauge named ', trim(PresentGauge%Name)
                    write(*,*) 'has a phase above +360'
                    stop
                endif


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
    !--------------------------------------------------------------------------

    subroutine GetWave(PresentGauge, WaveName, WaveExists, WaveOut)

        !Arguments-------------------------------------------------------------
        type(T_TideGauge), pointer                  :: PresentGauge 
        character(len=*),            intent (IN )   :: WaveName
        logical,                     intent (OUT)   :: WaveExists        
        type(T_TidalWave), pointer                  :: WaveOut

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        WaveOut => PresentGauge%FirstWave
        WaveExists = .false.

do2 :   do while (associated(WaveOut)) 
            if (trim(WaveOut%Name)==trim(WaveName)) then
                WaveExists = .true.
                exit
            endif
            WaveOut => WaveOut%Next
        end do do2

        !----------------------------------------------------------------------

    end subroutine GetWave

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

    subroutine GetIJReferenceLevel_ThreadSafe(GaugeID, HorizontalGridID, I, J, ReferenceLevel, STAT)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )  :: GaugeID
        integer,           intent(IN )  :: HorizontalGridID
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
        type (T_Gauge), pointer         :: LocalMe

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        LocalMe => Ready_ThreadSafe(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            PresentGauge => LocalMe%FirstGauge

            Found = .false.

            do while (associated(PresentGauge) .and. .not. Found) 

                if (PresentGauge%Grid_I == I .and. PresentGauge%Grid_J == J) then

                    ReferenceLevel = PresentGauge%ReferenceLevel                    

                    Found = .true.

                endif
                
                PX = PresentGauge%Metric_X
                PY = PresentGauge%Metric_Y
                
                call GetXYCellZ_ThreadSafe(HorizontalGridID, PX, PY, Iaux, Jaux, STAT = STAT_)

  
                if (I == Iaux .and. J == Jaux) then
                
                    call GetXYCellZ_ThreadSafe(HorizontalGridID, PX, PY, Iaux, Jaux, PercI = PercI, PercJ = PercJ, STAT = STAT_)

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

    end subroutine GetIJReferenceLevel_ThreadSafe

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

    subroutine GetIJWaterLevel_ThreadSafe(GaugeID, HorizontalGridID, I, J, WaterLevel, STAT)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )  :: GaugeID
        integer,           intent(IN )  :: HorizontalGridID
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
        type (T_Gauge), pointer         :: LocalMe

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        LocalMe => Ready_ThreadSafe(GaugeID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            PresentGauge => LocalMe%FirstGauge

            Found = .false.

            do while (associated(PresentGauge) .and. .not. Found) 

                if (PresentGauge%Grid_I == I .and. PresentGauge%Grid_J == J) then

                    WaterLevel = PresentGauge%WaterLevel                    

                    Found      = .true.

                endif
                
                PX = PresentGauge%Metric_X
                PY = PresentGauge%Metric_Y                
                
                call GetXYCellZ_ThreadSafe(HorizontalGridID, PX, PY, Iaux, Jaux, STAT = STAT_)

  
                if (I == Iaux .and. J == Jaux) then

                    call GetXYCellZ_ThreadSafe(HorizontalGridID, PX, PY, Iaux, Jaux, PercI = PercI, PercJ = PercJ, STAT = STAT_)

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

    end subroutine GetIJWaterLevel_ThreadSafe
    
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
        type(T_TideGauge), pointer                              :: PresentGauge
        type(T_Time), intent(IN)                                :: TimeIn
        real                                                    :: WaterLevel

        !Local-----------------------------------------------------------------
        real,   pointer, dimension(:)                           :: WaveAmplitude
        real,   pointer, dimension(:)                           :: WavePhase
        integer                                                 :: status, NWaves, J, AdmitNWaves, Nmax
        type(T_TidalWave), pointer                              :: PresentWave
        real                                                    :: WaterLevelRef
        character(LEN=WaveNameLength), pointer, dimension(:)    :: WaveName


        !Begin------------------------------------------------------------------                         

        call GetNumberWaves(PresentGauge, NWaves)
        
        ! JPW: Define to prevent using random value
        Nmax = NWaves
        if (Me%ComputeAdmittance) then
            Nmax = NWaves + NAdmit
        endif
        
        if (Me%TidePREDICTION == Task2000_) then
            Nmax = max(Nmax, NTask2000)
        endif

        allocate(WaveAmplitude(Nmax), STAT = status)
        if (status /= SUCCESS_)                                                         &
            call SetError(FATAL_, INTERNAL_, "HarmonicsWaterLevel - ModuleGauge - ERR10")

        allocate(WavePhase    (Nmax), STAT = status)
        if (status /= SUCCESS_)                                                         &
            call SetError(FATAL_, INTERNAL_, "HarmonicsWaterLevel - ModuleGauge - ERR20")

        allocate(WaveName     (Nmax), STAT = status)
        if (status /= SUCCESS_)                                                         &
            call SetError(FATAL_, INTERNAL_, "HarmonicsWaterLevel - ModuleGauge - ERR30")

        J = 0
        PresentWave => PresentGauge%FirstWave
do5 :   do while (associated(PresentWave))
            J = J + 1

            WaveAmplitude (J) = PresentWave%Amplitude
            WavePhase     (J) = PresentWave%Phase
            WaveName      (J) = PresentWave%Name

            PresentWave => PresentWave%Next
        end do do5

        ! Reference level is add in the module open boundary
        ! The water level is the Gauge model is compute with relation to 
        ! the hydrographic zero
        WaterLevelRef = 0.
        
        if (Me%ComputeAdmittance) then
        
            call NewComponentsByAdmittance(PresentGauge, NWaves, WaveAmplitude,         &
                                           WavePhase, WaveName, AdmitNWaves)

            NWaves = NWaves + AdmitNWaves
        
        endif
        
        if (Me%TidePREDICTION == Toga_) then

            call TogaLevel(PresentGauge%ObjToga,                                        &  
                           PresentGauge%WaterLevel,                                     &
                           PresentGauge%DecimalLatitude,                                &
                           PresentGauge%TimeReference,                                  &
                           WaterLevelRef,                                               &
                           NWaves, WaveAmplitude, WavePhase,                            &
                           TimeIn,                                                    &
                           STAT = status)

            if (status /= SUCCESS_)                                                     &
                call SetError(FATAL_, INTERNAL_, "HarmonicsWaterLevel - ModuleGauge - ERR30")
                           
        
        else if (Me%TidePREDICTION == Task2000_) then
        
            call Task2000Level(PresentGauge%WaterLevel,                                 &
                               PresentGauge%TimeReference,                              &
                               NWaves, WaveAmplitude,                                   &
                               WavePhase, WaveName, TimeIn, STAT = status)

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
            
        deallocate(WaveName)

        nullify   (WaveAmplitude)
        nullify   (WavePhase    )
        nullify   (WaveName     )

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

    
    !Convert amplitude and phase (degrees / 360º) in complex number (Sreal+ i * Simag)
   
    subroutine AmpPhase_To_Complex (Amplitude, Phase, Sreal, Simag)

        !Arguments-------------------------------------------------------------
        real :: Amplitude, Phase, Sreal, Simag
        !Begin-----------------------------------------------------------------

        Sreal = Amplitude*cos(Phase * 2 * Pi)
        Simag = Amplitude*sin(Phase * 2 * Pi)

    end Subroutine AmpPhase_To_Complex
    
    
    !Convert complex number (Sreal+ i * Simag) in amplitude and phase (degrees / 360º)
    subroutine Complex_to_AmpPhase (Sreal,Simag,Amplitude,Phase)

        !Arguments-------------------------------------------------------------    
        real :: Amplitude, Phase, Sreal, Simag

        !Begin-----------------------------------------------------------------    
        Amplitude = sqrt(Sreal**2.+Simag**2.)
        Phase = Atan (Simag/Sreal)

        if(Sreal < 0 .And. Simag < 0)Then
           Phase = Phase + Pi
        end If

        if(Sreal < 0 .And. Simag > 0)Then
          Phase = Phase - Pi
        end If
        
        Phase = Phase / (2 * Pi)

    end subroutine Complex_to_AmpPhase

!______________________________________________________
!  2Q1
!  computes the real part of 2Q1 from Q1 and O1
!------------------------------------------------------

    real function Real2Q1_Comp (realQ1,realO1)

        !Arguments-------------------------------------------------------------    
        real:: realQ1,realO1 
        !Begin-----------------------------------------------------------------    
        
        Real2Q1_Comp = 0.263*realQ1-0.0252*realO1

    end function Real2Q1_Comp

!______________________________________________________
!  2Q1
!  computes the imaginary part of 2Q1 from Q1 and O1
!------------------------------------------------------

    real function Imag2Q1_Comp (ImagQ1,ImagO1)

        !Arguments-------------------------------------------------------------    
        real:: ImagQ1,ImagO1 
        !Begin-----------------------------------------------------------------            

        Imag2Q1_Comp = 0.263*ImagQ1-0.0252*ImagO1

    end function Imag2Q1_Comp

!______________________________________________________
!  Sigma1
!  computes the real part of Sigma1 from Q1 and O1
!------------------------------------------------------

    real function realSigma1_Comp (realQ1,realO1)

       !Arguments-------------------------------------------------------------    
        real:: realQ1,realO1 
        !Begin-----------------------------------------------------------------            

        realSigma1_Comp = 0.297*realQ1-0.0264*realO1

    end function

!______________________________________________________
!  Sigma1
!computes the imaginary part of Sigma1 from Q1 and O1
!------------------------------------------------------

    real function ImagSigma1_Comp (ImagQ1,ImagO1)

        !Arguments-------------------------------------------------------------    
        real:: ImagQ1,ImagO1 
        !Begin-----------------------------------------------------------------            

        ImagSigma1_Comp = 0.297*ImagQ1-0.0264*ImagO1

    end function ImagSigma1_Comp

!______________________________________________________
!  rho1
!computes the real part of rho1 from Q1 and O1
!------------------------------------------------------

    real function realrho1_Comp (realQ1,realO1)

       !Arguments-------------------------------------------------------------    
        real:: realQ1,realO1 
       !Begin-----------------------------------------------------------------            

        realrho1_Comp = 0.164*realQ1+0.0048*realO1

    end function realrho1_Comp

!______________________________________________________
!  rho1
!computes the imaginary part of rho1 from Q1 and O1
!------------------------------------------------------

    real function Imagrho1_Comp (ImagQ1,ImagO1)

       !Arguments-------------------------------------------------------------    
        real:: ImagQ1,ImagO1 
       !Begin-----------------------------------------------------------------            

        Imagrho1_Comp = 0.164*ImagQ1+0.0048*ImagO1

    end  function Imagrho1_Comp

!______________________________________________________
!  m11
!computes the real part of M11 from O1 and K1
!------------------------------------------------------

    real function realm11_Comp (realO1,realK1)

       !Arguments-------------------------------------------------------------    
        real:: realO1,realK1
       !Begin-----------------------------------------------------------------            

        realm11_Comp = 0.014*realO1+0.0101*realK1

    end  function realm11_Comp

!______________________________________________________
!  m11
!computes the imaginary part of M11 from O1 and K1
!------------------------------------------------------

    real function Imagm11_Comp (ImagO1,ImagK1)

       !Arguments-------------------------------------------------------------    
        Real:: ImagO1,ImagK1 
       !Begin-----------------------------------------------------------------            


        Imagm11_Comp = 0.014*ImagO1+0.0101*ImagK1

    end function Imagm11_Comp

!______________________________________________________
!  m12
!computes the real part of M12 from O1 and K1
!------------------------------------------------------

    real function realm12_Comp (realO1,realK1)

       !Arguments-------------------------------------------------------------    
        real:: realO1,realK1 
       !Begin-----------------------------------------------------------------            
       
        realm12_Comp = 0.0389*realO1+0.0282*realK1

    end function realm12_Comp

!______________________________________________________
!  m12
!computes the imaginary part of M12 from O1 and K1
!------------------------------------------------------


    real function Imagm12_Comp (ImagO1,ImagK1)

       !Arguments-------------------------------------------------------------    
        real:: ImagO1,ImagK1
       !Begin-----------------------------------------------------------------            
       
        Imagm12_Comp = 0.0389*ImagO1+0.0282*ImagK1

    end function Imagm12_Comp

!______________________________________________________
!  chi1
!computes the real part of chi1 from O1 and K1
!------------------------------------------------------

    real function realchi1_Comp (realO1,realK1)

       !Arguments-------------------------------------------------------------    
        real:: realO1,realK1
       !Begin-----------------------------------------------------------------            

        realchi1_Comp = 0.0064*realO1+0.0060*realK1

    end function realchi1_Comp
!______________________________________________________
!  chi1
!computes the real part of chi1 from O1 and K1
!------------------------------------------------------

    real function Imagchi1_Comp (ImagO1,ImagK1)

       !Arguments-------------------------------------------------------------    
        real:: ImagK1,ImagO1 
       !Begin-----------------------------------------------------------------            

        Imagchi1_Comp = 0.0064*ImagO1+0.0060*ImagK1

    end  function Imagchi1_Comp

!______________________________________________________
!  pi1
!computes the real part of pi1 from O1 and K1
!------------------------------------------------------

    real function realpi1_Comp (realO1,realK1)

       !Arguments-------------------------------------------------------------    
        real:: realO1,realK1 
       !Begin-----------------------------------------------------------------            
       
        realpi1_Comp = 0.0030*realO1+0.0171*realK1

    end function realpi1_Comp

!______________________________________________________
!  pi1
!computes the imaginary part of pi1 from O1 and K1
!------------------------------------------------------

    real function Imagpi1_Comp (ImagO1,ImagK1)

       !Arguments-------------------------------------------------------------    
        real:: ImagK1,ImagO1 
       !Begin-----------------------------------------------------------------            

        Imagpi1_Comp = 0.0030*ImagO1+0.0171*ImagK1

    end function Imagpi1_Comp

!______________________________________________________
!  phi1
!computes the real part of phi1 from O1 and K1
!------------------------------------------------------

    real function realphi1_Comp (realO1,realK1)

       !Arguments-------------------------------------------------------------    
        real:: realO1,realK1 
       !Begin-----------------------------------------------------------------            

        realphi1_Comp = -0.0015*realO1+0.0152*realK1

    end function realphi1_Comp 

!______________________________________________________
!  phi1
!computes the imaginary part of phi1 from O1 and K1
!------------------------------------------------------

    real function Imagphi1_Comp (ImagO1,ImagK1)

       !Arguments-------------------------------------------------------------    
        real:: ImagK1,ImagO1 
       !Begin-----------------------------------------------------------------            

        Imagphi1_Comp = -0.0015*ImagO1+0.0152*ImagK1

    end Function Imagphi1_Comp

!______________________________________________________
!  theta1
!computes the real part of theta1 from O1 and K1
!------------------------------------------------------

    real function realtheta1_Comp (realO1,realK1)
       
       !Arguments-------------------------------------------------------------    
        real::  realO1,realK1 
       !Begin-----------------------------------------------------------------            
       
        realtheta1_Comp = -0.0065*realO1+0.0155*realK1

    end function realtheta1_Comp

!______________________________________________________
!  theta1
!computes the imaginary part of theta1 from O1 and K1
!------------------------------------------------------

    real function Imagtheta1_Comp (ImagO1,ImagK1)

       !Arguments-------------------------------------------------------------    
        real::  ImagK1,ImagO1 
       !Begin-----------------------------------------------------------------            

        Imagtheta1_Comp = -0.0065*ImagO1+0.0155*ImagK1

    end function Imagtheta1_Comp

!______________________________________________________
!  J1
!computes the real part of J1 from O1 and K1
!------------------------------------------------------

    real function realJ1_Comp (realO1,realK1)

       !Arguments-------------------------------------------------------------    
        real::  realO1,realK1 
       !Begin-----------------------------------------------------------------            
       
        realJ1_Comp = -0.0389*realO1+0.0836*realK1


    end function realJ1_Comp
    
!______________________________________________________
!  J1
!computes the imaginary part of J1 from O1 and K1
!------------------------------------------------------    

    real function ImagJ1_Comp (ImagO1,ImagK1)

       !Arguments-------------------------------------------------------------    
        real::   ImagK1,ImagO1 
       !Begin-----------------------------------------------------------------            

        ImagJ1_Comp = -0.0389*ImagO1+0.0836*ImagK1

    end function ImagJ1_Comp


!______________________________________________________
!  001
!computes the real part of OO1 from O1 and K1
!------------------------------------------------------

    real function realOO1_Comp (realO1,realK1)

       !Arguments-------------------------------------------------------------    
        real::   realO1,realK1 
       !Begin-----------------------------------------------------------------            

        realOO1_Comp = -0.0431*realO1+0.0613*realK1

    end  function realOO1_Comp

!______________________________________________________
!  001
!computes the imaginary part of OO1 from O1 and K1
!------------------------------------------------------

    real function ImagOO1_Comp (ImagO1,ImagK1)

       !Arguments-------------------------------------------------------------    
        real:: ImagK1,ImagO1 
       !Begin-----------------------------------------------------------------            
       
        ImagOO1_Comp = -0.0431*ImagO1+0.0613*ImagK1

    end function ImagOO1_Comp

!______________________________________________________
!  Eps2
!computes the real part of Eps2 from 2N2 and N2
!------------------------------------------------------

    real function realEpsilon2_Comp (real2N2,realN2)

       !Arguments-------------------------------------------------------------    
        real:: real2N2,realN2 
       !Begin-----------------------------------------------------------------            
       
        realEpsilon2_Comp = 0.53285*real2N2-0.03304*realN2

    end function realEpsilon2_Comp

!______________________________________________________
!  Eps2
!computes the imaginary part of Eps2 from 2N2 and N2
!------------------------------------------------------

    real function ImagEpsilon2_Comp (Imag2N2,ImagN2)

       !Arguments-------------------------------------------------------------    
        real:: Imag2N2,ImagN2 
       !Begin-----------------------------------------------------------------            
       
        ImagEpsilon2_Comp = 0.53285*Imag2N2-0.03304*ImagN2

    end function ImagEpsilon2_Comp 

!______________________________________________________
!  Eta2
!computes the real part of Eta2 from M2 and K2
!------------------------------------------------------

   real function realEta2_Comp (realM2,realK2)

       !Arguments-------------------------------------------------------------    
        real    :: realM2,realK2
       !Begin-----------------------------------------------------------------            

        realEta2_Comp = -0.0034925*realM2+0.0831707*realK2


    end function realEta2_Comp 

!______________________________________________________
!  Eta2
!computes the imaginary part of Eta2 from M2 and K2
!------------------------------------------------------

   real function ImagEta2_Comp (ImagM2,ImagK2)

       !Arguments-------------------------------------------------------------    
        real    :: ImagM2,ImagK2 
       !Begin-----------------------------------------------------------------            

        ImagEta2_Comp =  -0.0034925*ImagM2+0.0831707*ImagK2

    end function ImagEta2_Comp

!______________________________________________________
!  P1
!computes the real part of P1 from Q1, O1 and K1
!------------------------------------------------------

    real function realP1_Comp (realQ1,realO1,realK1)

       !Arguments-------------------------------------------------------------    
        real    :: realQ1,realO1,realK1
       !Local-----------------------------------------------------------------            
        real    :: freqQ1,freqO1,freqK1,freqP1,alp1,alq1,frbar1,deno1,aap1
        real    :: alo1,alk1,bbp1,ccp1

       !Begin-----------------------------------------------------------------            

        freqQ1 = 13.39866087990d0
        freqO1 = 13.94303558000d0
        freqK1 = 15.04106864000d0
        freqP1 = 14.95893136000d0

        alq1  = 0.073017
        alo1  = 0.3771366
        alk1  = 0.5300728
        alp1  = 0.1750754

        frbar1=1/3.*(freqQ1 + freqO1 + freqK1)
        deno1=frbar1**2 - 1/3.*(freqQ1**2 + freqO1**2 + freqK1**2)

        aap1 = alp1/3./alq1*(1.-(freqP1-frbar1)*(freqQ1-frbar1)/deno1)
        bbp1 = alp1/3./alo1*(1.- (freqP1-frbar1)*(freqO1-frbar1)/deno1)
        ccp1 = alp1/3./alk1*(1.- (freqP1-frbar1)*(freqK1-frbar1)/deno1)

        realP1_Comp = aap1*realQ1 + bbp1*realO1+ ccp1*realK1
 
    end function realP1_Comp 
    
!______________________________________________________
!  P1
!computes the imaginary part of P1 from Q1, O1 and K1
!------------------------------------------------------

    real function ImagP1_Comp (ImagQ1,ImagO1,ImagK1)


        !Arguments-------------------------------------------------------------    
        real        :: ImagQ1,ImagO1,ImagK1
        !Local-----------------------------------------------------------------            
        real        :: freqQ1,freqO1,freqK1,freqP1,alp1,alq1,frbar1,deno1,aap1
        real        :: alo1,alk1,bbp1,ccp1

        !Begin-----------------------------------------------------------------            

        freqQ1 = 13.39866087990d0
        freqO1 = 13.94303558000d0
        freqK1 = 15.04106864000d0
        freqP1 = 14.95893136000d0

        alq1  = 0.073017
        alo1  = 0.3771366
        alk1  = 0.5300728
        alp1  = 0.1750754

        frbar1=1/3.*(freqQ1 + freqO1 + freqK1)
        deno1=frbar1**2 - 1/3.*(freqQ1**2 + freqO1**2 + freqK1**2)

        aap1 = alp1/3./alq1*(1.-(freqP1-frbar1)*(freqQ1-frbar1)/deno1)
        bbp1 = alp1/3./alo1*(1.- (freqP1-frbar1)*(freqO1-frbar1)/deno1)
        ccp1 = alp1/3./alk1*(1.- (freqP1-frbar1)*(freqK1-frbar1)/deno1)

        ImagP1_Comp = aap1*ImagQ1 + bbp1*ImagO1+ ccp1*ImagK1

    end function ImagP1_Comp

!______________________________________________________
!  Mu2
!computes the real part of Mu2 from K2, N2 and M2
!------------------------------------------------------

    real function realMu2_Comp (realK2,realN2,realM2)

        !Arguments-------------------------------------------------------------    
        real        :: realM2,realK2,realN2
        !Local-----------------------------------------------------------------            
        real        :: freqM2,freqN2,freqk2,freqmu2,almu2,ck,sk,deno,cn,sn
        real        :: alm2,alk2,aln2,bbmu2,ccmu2,aamu2,cmu2,smu2

        !Begin-----------------------------------------------------------------            

        freqM2 = 28.98410422000d0
        freqN2 = 28.43972952010d0
        freqK2 = 30.08213728000d0
        freqMu2 = 27.96820844000d0

        almu2 = 0.02777
        alk2   = 0.1149327
        aln2   = 0.1758941
        alm2   = 0.9085024

        ck=cos(2*pi*2.*(freqk2-freqm2)/15.) ! CPD
        sk=sin(2*pi*2.*(freqk2-freqm2)/15.) ! CPD

        cn=cos(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        sn=sin(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        deno=sk*(cn-1)-sn*(ck-1)

        cmu2=cos(2*pi*2.*(freqmu2-freqm2)/15.)
        smu2=sin(2*pi*2.*(freqmu2-freqm2)/15.)


        aamu2= (-sn*cmu2 +(cn-1)*smu2 +sn)/deno/alk2*almu2
        bbmu2=(sk*cmu2-(ck-1)*smu2-sk)/deno/aln2*almu2
        ccmu2=(-(sk-sn)*cmu2+(ck-cn)*smu2+sk*cn-sn*ck)/deno/alm2*almu2

        realMu2_Comp = aamu2*realK2 + bbmu2*realN2+ ccmu2*realM2

    end function realMu2_Comp


!______________________________________________________
!  Mu2
!computes the imaginary part of Mu2 from K2, N2 and M2
!------------------------------------------------------

    real function ImagMu2_Comp (ImagK2,ImagN2,ImagM2)

        !Arguments-------------------------------------------------------------    
        real        :: ImagM2,ImagK2,ImagN2
        !Local-----------------------------------------------------------------            
        real        :: freqM2,freqN2,freqK2,freqmu2,almu2,ck,sk,deno,cn,sn
        real        ::  alm2,alk2,aln2,bbmu2,ccmu2,aamu2,cmu2,smu2

        !Begin-----------------------------------------------------------------            

        freqM2 = 28.98410422000d0
        freqN2 = 28.43972952010d0
        freqK2 = 30.08213728000d0
        freqMu2 = 27.96820844000d0

        almu2 = 0.02777
        alk2   = 0.1149327
        aln2   = 0.1758941
        alm2   = 0.9085024

        ck=cos(2*pi*2.*(freqk2-freqm2)/15.) ! CPD
        sk=sin(2*pi*2.*(freqk2-freqm2)/15.) ! CPD

        cn=cos(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        sn=sin(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        deno=sk*(cn-1)-sn*(ck-1)

        cmu2=cos(2*pi*2.*(freqmu2-freqm2)/15.)
        smu2=sin(2*pi*2.*(freqmu2-freqm2)/15.)

        aamu2= (-sn*cmu2 +(cn-1)*smu2 +sn)/deno/alk2*almu2
        bbmu2=(sk*cmu2-(ck-1)*smu2-sk)/deno/aln2*almu2
        ccmu2=(-(sk-sn)*cmu2+(ck-cn)*smu2+sk*cn-sn*ck)/deno/alm2*almu2

        ImagMu2_Comp = aamu2*ImagK2 + bbmu2*ImagN2+ ccmu2*ImagM2

    end function ImagMu2_Comp

!______________________________________________________
!  Nu2
!computes the real part of Nu2 from K2, N2 and M2
!------------------------------------------------------

    real function realNu2_Comp (realK2,realN2,realM2)

        !Arguments-------------------------------------------------------------    
        real        :: realM2,realK2,realN2
        !Local-----------------------------------------------------------------            
        real        :: freqM2,freqN2,freqK2,freqnu2,alnu2,ck,sk,deno,cn,sn
        real        :: alm2,alk2,aln2,bbnu2,ccnu2,aanu2,cnu2,snu2

        !Begin-----------------------------------------------------------------            

        freqM2 = 28.98410422000d0
        freqN2 = 28.43972952010d0
        freqK2 = 30.08213728000d0
        freqNu2 = 28.51258314000d0

        alnu2 = 0.03303
        alk2   = 0.1149327
        aln2   = 0.1758941
        alm2   = 0.9085024

        ck=cos(2*pi*2.*(freqk2-freqm2)/15.) ! CPD
        sk=sin(2*pi*2.*(freqk2-freqm2)/15.) ! CPD

        cn=cos(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        sn=sin(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        deno=sk*(cn-1)-sn*(ck-1)

        cnu2=cos(2*pi*2.*(freqnu2-freqm2)/15.)
        snu2=sin(2*pi*2.*(freqnu2-freqm2)/15.)

        aanu2= (-sn*cnu2 +(cn-1)*snu2 +sn)/deno/alk2*alnu2
        bbnu2=(sk*cnu2-(ck-1)*snu2-sk)/deno/aln2*alnu2
        ccnu2=(-(sk-sn)*cnu2+(ck-cn)*snu2+sk*cn-sn*ck)/deno/alm2*alnu2

        realNu2_Comp = aanu2*realK2 + bbnu2*realN2+ ccnu2*realM2


    end function realNu2_Comp 

!______________________________________________________
!  Nu2
!computes the imaginary part of Nu2 from K2, N2 and M2
!------------------------------------------------------

    real function ImagNu2_Comp (ImagK2,ImagN2,ImagM2)

        !Arguments-------------------------------------------------------------    
        real        :: ImagM2,ImagK2,ImagN2
        !Local-----------------------------------------------------------------            
        real        :: freqM2,freqN2,freqK2,freqnu2,alnu2,ck,sk,deno,cn,sn
        real        :: alm2,alk2,aln2,bbnu2,ccnu2,aanu2,cnu2,snu2

        !Begin-----------------------------------------------------------------            

        freqM2 = 28.98410422000d0
        freqN2 = 28.43972952010d0
        freqK2 = 30.08213728000d0
        freqNu2 = 28.51258314000d0

        alnu2 = 0.03303
        alk2   = 0.1149327
        aln2   = 0.1758941
        alm2   = 0.9085024

        ck=cos(2*pi*2.*(freqk2-freqm2)/15.) ! CPD
        sk=sin(2*pi*2.*(freqk2-freqm2)/15.) ! CPD

        cn=cos(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        sn=sin(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        deno=sk*(cn-1)-sn*(ck-1)

        cnu2=cos(2*pi*2.*(freqnu2-freqm2)/15.)
        snu2=sin(2*pi*2.*(freqnu2-freqm2)/15.)

        aanu2= (-sn*cnu2 +(cn-1)*snu2 +sn)/deno/alk2*alnu2
        bbnu2=(sk*cnu2-(ck-1)*snu2-sk)/deno/aln2*alnu2
        ccnu2=(-(sk-sn)*cnu2+(ck-cn)*snu2+sk*cn-sn*ck)/deno/alm2*alnu2

        ImagNu2_Comp = aanu2*ImagK2 + bbnu2*ImagN2+ ccnu2*ImagM2


    end function ImagNu2_Comp 

!______________________________________________________
!  Lda2
!computes the real part of Lda2 from K2, N2 and M2
!------------------------------------------------------


    real function realLda2_Comp (realK2,realN2,realM2)

        !Arguments-------------------------------------------------------------    
        real        :: realM2,realK2,realN2
        !Local-----------------------------------------------------------------            
        real        :: freqM2,freqN2,freqK2,freqlda2,allda2,ck,sk,deno,cn,sn
        real        :: alm2,alk2,aln2,aalda2,bblda2,cclda2,clda2,slda2

        !Begin-----------------------------------------------------------------            

        freqM2 = 28.98410422000d0
        freqN2 = 28.43972952010d0
        freqK2 = 30.08213728000d0
        freqlda2 = 29.4556253d0

        allda2= 0.0066
        alk2   = 0.1149327
        aln2   = 0.1758941
        alm2   = 0.9085024

        ck=cos(2*pi*2.*(freqk2-freqm2)/15.) ! CPD
        sk=sin(2*pi*2.*(freqk2-freqm2)/15.) ! CPD

        cn=cos(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        sn=sin(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        deno=sk*(cn-1)-sn*(ck-1)

        clda2=cos(2*pi*2.*(freqlda2-freqm2)/15.)
        slda2=sin(2*pi*2.*(freqlda2-freqm2)/15.)

        aalda2=(-sn*clda2+(cn-1)*slda2+sn)/deno/alk2*allda2
        bblda2=(sk*clda2-(ck-1)*slda2-sk)/deno/aln2*allda2
        cclda2=(-(sk-sn)*clda2+(ck-cn)*slda2+sk*cn-sn*ck)/deno/alm2*allda2

        realLda2_Comp = aalda2*realK2 + bblda2*realN2+ cclda2*realM2


    end function realLda2_Comp

!______________________________________________________
!  Lda2
!computes the imaginary part of Lda2 from K2, N2 and M2
!------------------------------------------------------

    real function ImagLda2_Comp (ImagK2,ImagN2,ImagM2)

        !Arguments-------------------------------------------------------------    
        real        :: ImagM2,ImagK2,ImagN2
        !Local-----------------------------------------------------------------            
        real        :: freqM2,freqN2,freqK2,freqlda2,allda2,ck,sk,deno,cn,sn
        real        :: alm2,alk2,aln2,aalda2,bblda2,cclda2,clda2,slda2

        !Begin-----------------------------------------------------------------  
        
        freqM2 = 28.98410422000d0
        freqN2 = 28.43972952010d0
        freqK2 = 30.08213728000d0
        freqlda2 = 29.4556253d0

        allda2= 0.0066
        alk2   = 0.1149327
        aln2   = 0.1758941
        alm2   = 0.9085024

        ck=cos(2*pi*2.*(freqk2-freqm2)/15.) ! CPD
        sk=sin(2*pi*2.*(freqk2-freqm2)/15.) ! CPD

        cn=cos(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        sn=sin(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        deno=sk*(cn-1)-sn*(ck-1)

        clda2=cos(2*pi*2.*(freqlda2-freqm2)/15.)
        slda2=sin(2*pi*2.*(freqlda2-freqm2)/15.)

        aalda2=(-sn*clda2+(cn-1)*slda2+sn)/deno/alk2*allda2
        bblda2=(sk*clda2-(ck-1)*slda2-sk)/deno/aln2*allda2
        cclda2=(-(sk-sn)*clda2+(ck-cn)*slda2+sk*cn-sn*ck)/deno/alm2*allda2


        ImagLda2_Comp = aalda2*ImagK2 + bblda2*ImagN2+ cclda2*ImagM2

    end function ImagLda2_Comp

!______________________________________________________
!  L2
!computes the real part of L2 from K2, N2 and M2
!------------------------------------------------------

    real Function realL2_Comp (realK2,realN2,realM2)
        !Arguments-------------------------------------------------------------    
        real        :: realM2,realK2,realN2
        !Local-----------------------------------------------------------------            
        real        :: freqM2,freqN2,freqK2,freqL2,all2,ck,sk,deno,cn,sn
        real        :: alm2,alk2,aln2,aal2,bbl2,ccl2,cl2,sl2

        !Begin-----------------------------------------------------------------  

        freqM2 = 28.98410422000d0
        freqN2 = 28.43972952010d0
        freqK2 = 30.08213728000d0
        freqL2 = 29.52847892000d0
                 

        all2  = 0.0251
        alk2   = 0.1149327
        aln2   = 0.1758941
        alm2   = 0.9085024

        ck=cos(2*pi*2.*(freqk2-freqm2)/15.) ! CPD
        sk=sin(2*pi*2.*(freqk2-freqm2)/15.) ! CPD

        cn=cos(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        sn=sin(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        deno=sk*(cn-1)-sn*(ck-1)

        cl2=cos(2*pi*2.*(freqL2-freqm2)/15.)
        sl2=sin(2*pi*2.*(freqL2-freqm2)/15.)

        aal2=  (-sn*cl2  +(cn-1)*sl2  +sn)/deno/alk2*all2
        bbl2=(sk*cl2-(ck-1)*sl2-sk)/deno/aln2*all2
        ccl2=(-(sk-sn)*cl2+(ck-cn)*sl2+sk*cn-sn*ck)/deno/alm2*all2

        realL2_Comp = aal2*realK2 + bbl2*realN2+ ccl2*realM2

    end function realL2_Comp 

!______________________________________________________
!  L2
!computes the imaginary part of L2 from K2, N2 and M2
!------------------------------------------------------

    real function ImagL2_Comp (ImagK2,ImagN2,ImagM2)

        !Arguments-------------------------------------------------------------    
        real        :: ImagM2,ImagK2,ImagN2
        !Local-----------------------------------------------------------------            
        real        :: freqM2,freqN2,freqK2,freqL2,all2,ck,sk,deno,cn,sn
        real        :: alm2,alk2,aln2,aal2,bbl2,ccl2,cl2,sl2

        !Begin-----------------------------------------------------------------  
        
        freqM2 = 28.98410422000d0
        freqN2 = 28.43972952010d0
        freqK2 = 30.08213728000d0
        freqL2 = 29.52847892000d0

        all2  = 0.0251 
        alk2   = 0.1149327
        aln2   = 0.1758941
        alm2   = 0.9085024

        ck=cos(2*pi*2.*(freqk2-freqm2)/15.) ! CPD
        sk=sin(2*pi*2.*(freqk2-freqm2)/15.) ! CPD

        cn=cos(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        sn=sin(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        deno=sk*(cn-1)-sn*(ck-1)

        cl2=cos(2*pi*2.*(freqL2-freqm2)/15.)
        sl2=sin(2*pi*2.*(freqL2-freqm2)/15.)

        aal2=  (-sn*cl2  +(cn-1)*sl2  +sn)/deno/alk2*all2
        bbl2=(sk*cl2-(ck-1)*sl2-sk)/deno/aln2*all2
        ccl2=(-(sk-sn)*cl2+(ck-cn)*sl2+sk*cn-sn*ck)/deno/alm2*all2

        ImagL2_Comp = aal2*ImagK2 + bbl2*ImagN2+ ccl2*ImagM2

    end function ImagL2_Comp

!______________________________________________________
!  T2
!computes the real part of T2 from K2, N2 and M2
!------------------------------------------------------
    real function realT2_Comp (realK2,realN2,realM2)

        !Arguments-------------------------------------------------------------    
        real        :: realK2,realN2,realM2
        !Local-----------------------------------------------------------------            
        real        :: freqM2,freqN2,freqK2,freqT2,alT2,ck,sk,deno,cn,sn
        real        :: alm2,alk2,aln2,aaT2,bbT2,ccT2,cT2,sT2

        !Begin-----------------------------------------------------------------  

        freqM2 = 28.98410422000d0
        freqN2 = 28.43972952010d0
        freqK2 = 30.08213728000d0
        freqT2 = 29.95893332010d0
                 
        alt2  = 0.0247766
        alk2   = 0.1149327
        aln2   = 0.1758941
        alm2   = 0.9085024


        ck=cos(2*pi*2.*(freqk2-freqm2)/15.) ! CPD
        sk=sin(2*pi*2.*(freqk2-freqm2)/15.) ! CPD

        cn=cos(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        sn=sin(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        deno=sk*(cn-1)-sn*(ck-1)

        cT2=cos(2*pi*2.*(freqT2-freqm2)/15.)
        sT2=sin(2*pi*2.*(freqT2-freqm2)/15.)

        aat2=  (-sn*ct2  +(cn-1)*st2  +sn)/deno/alk2*alt2
        bbt2=(sk*ct2-(ck-1)*st2-sk)/deno/aln2*alt2
        cct2=(-(sk-sn)*ct2+(ck-cn)*st2+sk*cn-sn*ck)/deno/alm2*alt2

        realT2_Comp = aaT2*realK2 + bbT2*realN2+ ccT2*realM2

    end function realT2_Comp

!______________________________________________________
!  T2
!computes the imaginary part of T2 from K2, N2 and M2
!------------------------------------------------------
    
    real function ImagT2_Comp (ImagK2,ImagN2,ImagM2)

        !Arguments-------------------------------------------------------------    
        real        :: ImagM2,ImagK2,ImagN2
        !Local-----------------------------------------------------------------            
        real        :: freqM2,freqN2,freqK2,freqT2,alT2,ck,sk,deno,cn,sn
        real        :: alm2,alk2,aln2,aaT2,bbT2,ccT2,cT2,sT2

        !Begin-----------------------------------------------------------------  

        freqM2 = 28.98410422000d0
        freqN2 = 28.43972952010d0
        freqK2 = 30.08213728000d0
        freqT2 = 29.95893332010d0

        alt2  = 0.0247766
        alk2   = 0.1149327
        aln2   = 0.1758941
        alm2   = 0.9085024

        ck=cos(2*pi*2.*(freqk2-freqm2)/15.) ! CPD
        sk=sin(2*pi*2.*(freqk2-freqm2)/15.) ! CPD

        cn=cos(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        sn=sin(2*pi*2.*(freqn2-freqm2)/15.) ! CPD
        deno=sk*(cn-1)-sn*(ck-1)

        cT2=cos(2*pi*2.*(freqT2-freqm2)/15.)
        sT2=sin(2*pi*2.*(freqT2-freqm2)/15.)

        aat2=  (-sn*ct2  +(cn-1)*st2  +sn)/deno/alk2*alt2
        bbt2=(sk*ct2-(ck-1)*st2-sk)/deno/aln2*alt2
        cct2=(-(sk-sn)*ct2+(ck-cn)*st2+sk*cn-sn*ck)/deno/alm2*alt2

        ImagT2_Comp = aaT2*ImagK2 + bbT2*ImagN2+ ccT2*ImagM2

    end function

!______________________________________________________
!  End of admittance functions
!------------------------------------------------------


    integer function SearchName(WaveName, n, Name)
    
        !Arguments-------------------------------------------------------------    
        character(StringLength), dimension(:)      :: WaveName 
        integer                                    :: n
        character(*)                               :: Name 
        !Local-----------------------------------------------------------------            
        integer                                    :: i, nl, nl1

        !Begin-----------------------------------------------------------------  
    
        nl = len_Trim(Name)

        do i=1,n
            nl1=len_trim(WaveName(i))

            if (WaveName(i)(1:nl)==Name(1:nl) .and. nl==nl1) exit
        enddo

        SearchName = i

    end function SearchName


    !Add new components to the tidal gauge by atmittance
    subroutine NewComponentsByAdmittance(PresentGauge, NWaves, WaveAmplitude,           &
                                         WavePhase, WaveName, AdmitNWaves)

        !Arguments------------------------------------------------
        type(T_TideGauge), pointer                      :: PresentGauge
        integer,            intent(IN)                  :: NWaves
        real,   pointer, dimension(:)                   :: WaveAmplitude
        real,   pointer, dimension(:)                   :: WavePhase
        character(LEN=*), pointer, dimension(:)         :: WaveName
        integer,            intent(OUT)                 :: AdmitNWaves
       
        !Local----------------------------------------------------        
        real                                            :: RealNew, ImagNew
        integer                                         :: n
        real                                            :: realQ1, realO1, realK1, real2N2, realN2, realM2, realK2
        real                                            :: ImagQ1, ImagO1, ImagK1, Imag2N2, ImagN2, ImagM2, ImagK2
        logical                                         :: Q1, O1, K1, var_2N2, N2, M2, K2, ExistComp
        type(T_TidalWave), pointer                      :: WaveOut
 
        !Begin----------------------------------------------------        

        call GetWave(PresentGauge, 'Q1', Q1, WaveOut)
        if (Q1) then
            call AmpPhase_To_Complex (WaveOut%Amplitude, WaveOut%Phase, realQ1, ImagQ1)
        else
            realQ1 = FillValueReal
            ImagQ1 = FillValueReal
        endif

        call GetWave(PresentGauge, 'O1', O1, WaveOut)
        if (O1) then
            call AmpPhase_To_Complex (WaveOut%Amplitude, WaveOut%Phase, realO1, ImagO1)
        else
            realO1 = FillValueReal
            ImagO1 = FillValueReal
        endif

        call GetWave(PresentGauge, 'K1', K1, WaveOut)
        if (K1) then
            call AmpPhase_To_Complex (WaveOut%Amplitude, WaveOut%Phase, realK1, ImagK1)
        else
            realK1 = FillValueReal
            ImagK1 = FillValueReal
        endif

        call GetWave(PresentGauge, '2N2', var_2N2, WaveOut)
        if (var_2N2) then
            call AmpPhase_To_Complex (WaveOut%Amplitude, WaveOut%Phase, real2N2, Imag2N2)
        else
            real2N2 = FillValueReal
            Imag2N2 = FillValueReal
        endif

        call GetWave(PresentGauge, 'N2', N2, WaveOut)
        if (N2) then
            call AmpPhase_To_Complex (WaveOut%Amplitude, WaveOut%Phase, realN2, ImagN2)
        else
            realN2 = FillValueReal
            ImagN2 = FillValueReal
        endif

        call GetWave(PresentGauge, 'K2', K2, WaveOut)
        if (K2) then
            call AmpPhase_To_Complex (WaveOut%Amplitude, WaveOut%Phase, realK2, ImagK2)
        else
            realK2 = FillValueReal
            ImagK2 = FillValueReal
        endif

        call GetWave(PresentGauge, 'M2', M2, WaveOut)
        if (M2) then
            call AmpPhase_To_Complex (WaveOut%Amplitude, WaveOut%Phase, realM2, ImagM2)
        else
            realM2 = FillValueReal
            ImagM2 = FillValueReal
        endif
        
       
        AdmitNWaves = 0
        
Q1O1:   if (Q1 .and. O1) then
            call GetWave(PresentGauge, '2Q1', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component 2Q1 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = Real2Q1_Comp (realQ1,realO1)
                ImagNew = Imag2Q1_Comp (ImagQ1,ImagO1)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = '2Q1'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif

            call GetWave(PresentGauge, 'SIG1', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component SIG1 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realSigma1_Comp (realQ1,realO1)
                ImagNew = ImagSigma1_Comp (ImagQ1,ImagO1)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'SIG1'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif
            
            call GetWave(PresentGauge, 'RHO1', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component RHO1 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realrho1_Comp (realQ1,realO1)
                ImagNew = Imagrho1_Comp (ImagQ1,ImagO1)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'RHO1'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif
            
        endif Q1O1
        
O1K1:   if (O1 .and. K1) then

            call GetWave(PresentGauge, 'CHI1', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component CHI1 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realchi1_Comp (realO1,realK1)
                ImagNew = Imagchi1_Comp (ImagO1,ImagK1)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'CHI1'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif

            call GetWave(PresentGauge, 'PI1', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component PI1 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realpi1_Comp (realO1,realK1)
                ImagNew = Imagpi1_Comp (ImagO1,ImagK1)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'PI1'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif
            
            call GetWave(PresentGauge, 'PHI1', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component PHI1 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realphi1_Comp (realO1,realK1)
                ImagNew = Imagphi1_Comp (ImagO1,ImagK1)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'PHI1'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif            
            
            call GetWave(PresentGauge, 'THE1', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component THE1 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realtheta1_Comp (realO1,realK1)
                ImagNew = Imagtheta1_Comp (ImagO1,ImagK1)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'THE1'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif            

            call GetWave(PresentGauge, 'J1', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component J1 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realJ1_Comp (realO1,realK1)
                ImagNew = ImagJ1_Comp (ImagO1,ImagK1)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'J1'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif           

            call GetWave(PresentGauge, 'OO1', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component OO1 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realOO1_Comp (realO1,realK1)
                ImagNew = ImagOO1_Comp (ImagO1,ImagK1)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'OO1'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif    

            if (Me%TidePREDICTION == Toga_) then
                !M12 component not recognized by task2000
                call GetWave(PresentGauge, 'M12', ExistComp, WaveOut)
                if (.not. ExistComp) then
                !Add component M12 -----------------------------------
                    AdmitNWaves = AdmitNWaves + 1
                    
                    RealNew = realm12_Comp (realO1,realK1)
                    ImagNew = Imagm12_Comp (ImagO1,ImagK1)
                    
                    n       = NWaves + AdmitNWaves

                    WaveName(n) = 'M12'
                    
                    call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
                endif                
            endif
                                    
        endif O1K1

N22N2:  if (Var_2N2 .and. N2) then

            call GetWave(PresentGauge, 'EPS2', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component EPS2 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realEpsilon2_Comp (real2N2,realN2)
                ImagNew = ImagEpsilon2_Comp (Imag2N2,ImagN2)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'EPS2'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif   

        endif N22N2


M2K2:   if (M2 .and. K2) then

            call GetWave(PresentGauge, 'ETA2', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component ETA2 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realEta2_Comp (realM2,realK2)
                ImagNew = ImagEta2_Comp (ImagM2,ImagK2)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'ETA2'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif   

        endif M2K2


Q1O1K1: if (Q1 .and. O1 .and. K1) then

            call GetWave(PresentGauge, 'P1', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component P1 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realP1_Comp (realQ1,realO1,realK1)
                ImagNew = ImagP1_Comp (ImagQ1,ImagO1,ImagK1)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'P1'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif   

        endif Q1O1K1


K2N2M2: if (K2 .and. N2 .and. M2) then

            call GetWave(PresentGauge, 'MU2', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component MU2 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realMu2_Comp (realK2,realN2,realM2)
                ImagNew = ImagMu2_Comp (ImagK2,ImagN2,ImagM2)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'MU2'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif   

            call GetWave(PresentGauge, 'NU2', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component NU2 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realNu2_Comp (realK2,realN2,realM2)
                ImagNew = ImagNu2_Comp (ImagK2,ImagN2,ImagM2)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'NU2'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif   

            call GetWave(PresentGauge, 'LDA2', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component LDA2 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realLda2_Comp (realK2,realN2,realM2)
                ImagNew = ImagLda2_Comp (ImagK2,ImagN2,ImagM2)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'LDA2'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif               

            call GetWave(PresentGauge, 'L2', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component L2 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realL2_Comp (realK2,realN2,realM2)
                ImagNew = ImagL2_Comp (ImagK2,ImagN2,ImagM2)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'L2'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif 

            call GetWave(PresentGauge, 'T2', ExistComp, WaveOut)
            if (.not. ExistComp) then
            !Add component T2 -----------------------------------
                AdmitNWaves = AdmitNWaves + 1
                
                RealNew = realT2_Comp (realK2,realN2,realM2)
                ImagNew = ImagT2_Comp (ImagK2,ImagN2,ImagM2)
                
                n       = NWaves + AdmitNWaves

                WaveName(n) = 'T2'
                
                call Complex_To_AmpPhase(RealNew, ImagNew, WaveAmplitude(n), WavePhase(n))
            endif 
            
        endif K2N2M2


end subroutine NewComponentsByAdmittance

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
            ready_ = VerifyReadLock (mGAUGE_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------
    
    function Ready_ThreadSafe (ObjGauge_ID, ready_) result(LocalMe)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )              :: ObjGauge_ID
        integer,           intent(OUT)              :: ready_
        type (T_Gauge), pointer                     :: LocalMe
        
        !----------------------------------------------------------------------

        nullify (LocalMe)

cd1:    if (ObjGauge_ID > 0) then
            LocalMe => LocateObjGauge_ThreadSafe(ObjGauge_ID)
            ready_ = VerifyReadLock (mGAUGE_, LocalMe%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end function Ready_ThreadSafe    

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

    function LocateObjGauge_ThreadSafe (ObjGaugeID) result(LocatedMe)

        !Arguments-------------------------------------------------------------
        integer,           intent(IN )              :: ObjGaugeID
        type (T_Gauge), pointer                     :: LocatedMe

        !Local-----------------------------------------------------------------

        LocatedMe => FirstGauge
        do while (associated (LocatedMe))
            if (LocatedMe%InstanceID == ObjGaugeID) exit
            LocatedMe => LocatedMe%Next
        enddo

        if (.not. associated(LocatedMe)) stop 'ModuleGauge - LocateObjGauge - ERR01'
        
    end function LocateObjGauge_ThreadSafe
    
    !--------------------------------------------------------------------------

end module ModuleGauge

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
