!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Litter
! URL           : http://www.mohid.com
! AFFILIATION   : HIDROMOD
! DATE          : April 2018
! REVISION      : Paulo Leitao
! DESCRIPTION   : Module manages the beach and killing of lagrangian litter particles
!
!------------------------------------------------------------------------------


Module ModuleLitter

!<BeginLitter>
!   AGE_TO_BEACH            : seconds                     [ Litter age necessary to beach ] 
!   KILL_BEACH_LITTER       : 1(true)/0(false)            [By default (KILL_BEACH_LITTER : 1) 
!                                                          the litter particles removed from the lagrangian model]
!<<BeginBeachArea>>
!   NAME                    : char                        [Beach Area X]
!   DESCRIPTION                : char                     [Area where litter can beach]
!   VEL_THRESHOLD            : m/s                        [velocity below which litter can beach]
!   WAVE_THRESHOLD            : m                         [significant wave heiht limit below which litter can beach]
!   FILENAME                : char                        [file of polygons delimiting the area where litter can beach]
!   COAST_TYPE              : integer                     [1 - Type A, 2 - Type B, etc.]
!   PROBABILITY             : -                           [0 = 0% beach probability , 1 = 100% - beach probability  
!   AGE_LIMIT               : seconds                     [Age limit after beaching above which the particle 
!                                                          is delete from module litter]
!<<EndBeachArea>>
!<EndLitter>

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions    
    use ModuleDrawing
    use ModuleHDF5
    use ModuleHorizontalGrid

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructLitter
    interface  ConstructLitter
        module procedure ConstructLitterWater
        module procedure ConstructLitterLag
    end interface ConstructLitter   
    
    private ::      AllocateInstance

    !Selector
                     
    
    !Modifier
    public  :: ModifyLitter
    interface  ModifyLitter
        module procedure ModifyLitterWater
        module procedure ModifyLitterLag
    end interface ModifyLitter 
    
    private ::      CheckBeachLitter
    private ::      DeleteOldLitter    
    private ::      OutputNumberGrid

  

    !Destructor
    public  :: KillLitter                                                     

    !Management
    private ::      Ready
    private ::          LocateObjLitter 
    
    !Interfaces----------------------------------------------------------------

    !Types---------------------------------------------------------------------

    type T_IndividualArea
        integer                                         :: ID                   = null_int
        character(PathLength)                           :: FileName             = null_str
        character(PathLength)                           :: Name                 = null_str
        character(PathLength)                           :: Description          = null_str        
        type (T_Polygon), pointer                       :: Polygons             => null()
        integer                                         :: CoastType            = null_int
        real                                            :: VelThreshold         = null_real
        real                                            :: WaveThreshold        = null_real        
        real                                            :: Probability          = null_real     
        real                                            :: AgeLimit             = null_real     
    end type T_IndividualArea

    type T_BeachAreas
        type (T_IndividualArea), dimension(:), pointer  :: Individual           => null()           
        integer                                         :: Number               = null_int
    end type T_BeachAreas    
    
    !IO
    type T_Files
        character(PathLength)                           :: ConstructData        = null_str
        character(PathLength)                           :: ResultsHDF           = null_str
        character(PathLength)                           :: Nomfich              = null_str
        character(PathLength)                           :: LitterIni            = null_str        
        character(PathLength)                           :: LitterFin            = null_str                
    end type T_Files    
    
    type T_ExtVar
        integer                                         :: ObjTime              = 0
        type(T_Time)                                    :: CurrentTime
        type(T_Time)                                    :: NextCompute
        type(T_Time)                                    :: StartTime                      
        type(T_Time)                                    :: EndTime              
        logical                                         :: Backtracking         = .false. 
        type (T_polygon), pointer                       :: ModelDomain          => null()
        integer                                         :: nParticles           = null_int
        real(8), dimension(:), pointer                  :: Longitude            => null()
        real(8), dimension(:), pointer                  :: Latitude             => null()
        real(8), dimension(:), pointer                  :: Age                  => null()
        integer, dimension(:), pointer                  :: Origin               => null()
        integer, dimension(:), pointer                  :: ID                   => null()        
        logical, dimension(:), pointer                  :: Beach                => null()
        logical, dimension(:), pointer                  :: KillPartic           => null()
    end type T_ExtVar        
                                                                                
    type T_Particle
        integer                                         :: ID                   = null_int   
        type(T_Time)                                    :: BeachTime
        real(8)                                         :: Longitude            = null_real
        real(8)                                         :: Latitude             = null_real
        real(8)                                         :: Age                  = null_real
        integer                                         :: Origin               = null_int 
        integer                                         :: BeachAreaID          = null_int 
        integer                                         :: CoastType            = null_int         
        real(8)                                         :: BeachPeriod          = null_real
        type (T_Particle),     pointer                  :: Next                 => null()    
        type (T_Particle),     pointer                  :: Prev                 => null()            
    end type T_Particle        
    
    type T_ParticleList
        integer                                         :: Number               = null_int
        type (T_Particle), pointer                      :: First                => null() 
    end type T_ParticleList        
    
    type T_OutPut
        type (T_Time), dimension(:), pointer            :: OutTime
        integer                                         :: NextOutPut           = null_int
        integer                                         :: Number               = null_int
        real,          dimension(:,:), pointer          :: AuxReal2D            => null()
        character (len = PathLength)                    :: OutputFile           = null_str
        character (len = PathLength)                    :: HotStartFile         = null_str        
        character (len = PathLength)                    :: InputGridFile        = null_str        
        integer                                         :: ObjHDF5              = 0
        integer                                         :: ObjHDF5_2            = 0        
        integer                                         :: ObjHorizontalGrid    = 0 
        type (T_Size2D)                                 :: Size
        type (T_Size2D)                                 :: WorkSize
    end type T_OutPut    
    

    type T_OutputGrids
        type (T_Output), dimension(:), pointer          :: Individual           => null()           
        integer                                         :: Number               = null_int
        type (T_Time),   dimension(:), pointer  :: RestartOutTime       => null()
        integer                                 :: NextRestartOutPut    = null_int
        logical                                 :: WriteRestartFile     = .false. 
        logical                                 :: RestartOverwrite     = .false.
    end type T_OutputGrids       
    
    private :: T_Litter
    type       T_Litter
        integer                                         :: InstanceID           = null_int
        
        type (T_ExtVar)                                 :: ExtVar        
        
        real                                            :: AgeToBeach           = null_real
        logical                                         :: KillBeachLitter      = .false.  
        
        type (T_Time)                                   :: LastAtualization

        type (T_Files)                                  :: Files
        type (T_BeachAreas)                             :: BeachAreas
        type (T_ParticleList)                           :: ParticleList
        type (T_OutputGrids)                            :: OutputGrids

        integer                                         :: ClientNumber         = null_int
        integer                                         :: ObjEnterdata         = 0
        integer                                         :: ObjHDF5              = 0        
        
        type(T_Litter), pointer                         :: Next                 => null()
    end type  T_Litter

    !Global Module Variables
    type (T_Litter), pointer                            :: FirstObjLitter       => null()
    type (T_Litter), pointer                            :: Me                   => null()


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructLitterWater(ObjLitterID, TimeID, Nomfich, ModelDomain,          &
                                    ResultsHDF, STAT)

        !Arguments---------------------------------------------------------------
        integer            ,            intent(OUT)     :: ObjLitterID 
        integer            ,            intent(INOUT)   :: TimeID 
        character(len=*)   ,            intent(IN)      :: Nomfich
        type (T_Polygon), pointer,      intent(IN)      :: ModelDomain
        character(len=*)   ,            intent(IN)      :: ResultsHDF
        integer, optional, intent(OUT)                  :: STAT     

        !Local-------------------------------------------------------------------
        integer                                         :: ready_         
        integer                                         :: STAT_
        integer                                         :: STAT_CALL        

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mLitter_)) then
            nullify (FirstObjLitter)
            call RegisterModule (mLitter_) 
        endif

        call Ready(ObjLitterID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            Me%Files%Nomfich        = Nomfich

            Me%ExtVar%ModelDomain   => ModelDomain
            
            Me%ExtVar%ObjTime       = AssociateInstance (mTIME_, TimeID)
            
            Me%Files%ResultsHDF = ResultsHDF
            
            call GetExternalTime
            
            call ConstructFilesNames
            
            !Construct enter data 
            call ConstructEnterData(EnterDataID     = Me%ObjEnterData,                  &
                                    FileName        = Me%Files%ConstructData,           &
                                    ErrorMessage    = "ConstructLitterWater - ModuleLitter", &
                                    STAT            = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructLitterWater - ModuleLitter - ERR10'
            
            
            
            call ConstructFromLitterBlock
            
           
            !Kill enter data 
            call KillEnterData     (EnterDataID     = Me%ObjEnterData,                  &
                                    STAT            = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructLitterWater - ModuleLitter - ERR20'
            
            call ReadLitterBeach            
            
            !Returns ID
            ObjLitterID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleLitter - ConstructLitterWater - ERR100' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructLitterWater
 
    !--------------------------------------------------------------------------
    
        subroutine ConstructLitterLag(ObjLitterID, StartTime, ModelDomain, STAT)

        !Arguments---------------------------------------------------------------
        integer,                intent(INOUT)           :: ObjLitterID 
        integer, dimension(6),  intent(IN)              :: StartTime
        !1 - Xmin, 2 - Xmax, 3 - Ymin, 4 - Ymax
        real(8), dimension(4),  intent(IN)              :: ModelDomain
        integer, optional,      intent(OUT)             :: STAT     

        !Local-------------------------------------------------------------------
        real,   dimension(:), pointer                   :: VectorX, VectorY
        integer                                         :: ready_         
        integer                                         :: STAT_
        integer                                         :: STAT_CALL        

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mLitter_)) then
            nullify (FirstObjLitter)
            call RegisterModule (mLitter_) 
        endif

        call Ready(ObjLitterID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            Me%Files%Nomfich        = "Nomfich.dat"
            
            allocate(VectorX(5), VectorY(5))
            
            VectorX(1)  = ModelDomain(1)
            VectorX(2)  = ModelDomain(1)
            VectorX(3)  = ModelDomain(2)
            VectorX(4)  = ModelDomain(2)
            VectorX(5)  = ModelDomain(1)
            
            VectorY(1)  = ModelDomain(3)
            VectorY(2)  = ModelDomain(4)
            VectorY(3)  = ModelDomain(4)
            VectorY(4)  = ModelDomain(3)
            VectorY(5)  = ModelDomain(3)

            call New(Polygons   = Me%ExtVar%ModelDomain,                                &
                     VectorX    = VectorX,                                              & 
                     VectorY    = VectorY)     
            
            deallocate(VectorX, VectorY)
            
            Me%ExtVar%ObjTime       = 0
        
            !call GetExternalTime
            !Do not work backtracking
            
            Me%ExtVar%BackTracking = .false.
            call SetDate (Time1 = Me%ExtVar%StartTime,                                  &
                          Year  = StartTime(1),                                         &
                          Month = StartTime(2),                                         &
                          Day   = StartTime(3),                                         &
                          Hour  = StartTime(4),                                         &
                          Minute= StartTime(5),                                         &
                          Second= StartTime(6))
            
            Me%ExtVar%CurrentTime   = Me%ExtVar%StartTime                        
            Me%LastAtualization     = Me%ExtVar%StartTime  
            
            call ConstructFilesNames
            
            !Construct enter data 
            call ConstructEnterData(EnterDataID     = Me%ObjEnterData,                  &
                                    FileName        = Me%Files%ConstructData,           &
                                    ErrorMessage    = "ConstructLitterWater - ModuleLitter", &
                                    STAT            = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructLitterLag - ModuleLitter - ERR10'
            
            
            
            call ConstructFromLitterBlock
            
           
            !Kill enter data 
            call KillEnterData     (EnterDataID     = Me%ObjEnterData,                  &
                                    STAT            = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructLitterLag - ModuleLitter - ERR20'
            
            call ReadLitterBeach            
            
            !Returns ID
            ObjLitterID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleLitter - ConstructLitterLag - ERR100' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructLitterLag
 
    !--------------------------------------------------------------------------
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Litter), pointer                         :: NewObjLitter
        type (T_Litter), pointer                         :: PreviousObjLitter


        !Allocates new instance
        allocate (NewObjLitter)
        nullify  (NewObjLitter%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjLitter)) then
            FirstObjLitter         => NewObjLitter
            Me                    => NewObjLitter
        else
            PreviousObjLitter      => FirstObjLitter
            Me                    => FirstObjLitter%Next
            do while (associated(Me))
                PreviousObjLitter  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjLitter
            PreviousObjLitter%Next => NewObjLitter
        endif

        Me%InstanceID = RegisterNewInstance (mLitter_)


    end subroutine AllocateInstance


    !------------------------------------------------------------------------------
    
    subroutine GetExternalTime
    
        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                 :: STAT_CALL
        !Begin-----------------------------------------------------------------
    
        !Gets Time
        call GetComputeTimeLimits(Me%ExtVar%ObjTime,                                    &
                                  BeginTime = Me%ExtVar%StartTime,                      &
                                  EndTime   = Me%ExtVar%EndTime,                        &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetExternalTime - ModuleLitter - ERR10'


        ! Check if the simulation goes backward in time or forward in time (default mode)
        call GetBackTracking(Me%ExtVar%ObjTime,                                         &
                             Me%ExtVar%BackTracking, STAT = STAT_CALL)                    
        if (STAT_CALL /= SUCCESS_) stop 'GetExternalTime - ModuleLitter - ERR20'

        Me%ExtVar%CurrentTime   = Me%ExtVar%StartTime                        
        Me%LastAtualization     = Me%ExtVar%StartTime    
            
            
    end subroutine GetExternalTime
    
    !------------------------------------------------------------------------------    
    

    subroutine ConstructFilesNames

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL !, i
        character(len = StringLength)               :: Message

        !----------------------------------------------------------------------

        
        !Input data file
        Message   ='ASCII file used to construct litter module'
                                        
        call ReadFileName  (KEYWORD     = 'PARTIC_DATA',                                &
                            FILE_NAME   = Me%Files%ConstructData,                       &
                            Message     = Message,                                      &
                            FilesInput  = Me%Files%Nomfich,                             &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructFilesNames - ModuleLitter - ERR10'



cd1 :   if      (STAT_CALL .EQ. FILE_NOT_FOUND_ERR_  ) then
            write(*,*)  
            write(*,*) 'Initial file not found.'
            
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFilesNames - ModuleLitter - ERR20'

        else if (STAT_CALL .EQ. KEYWORD_NOT_FOUND_ERR_) then
            write(*,*)  
            write(*,*) 'Keyword for the initial file not found in nomfich.dat.'
            write(*,*) 'ConstructFilesNames - ModuleLitter - WRN01'
            write(*,*)  

        else if (STAT_CALL .EQ. SUCCESS_             ) then
            continue
        else
            stop 'ConstructFilesNames - ModuleLitter - ERR30'
        end if cd1  
        
        
        !!Transient HDF File
        !Message   ='Instant fields of litter particles in HDF format.'
        !call ReadFileName  (KEYWORD     = 'PARTIC_HDF',                                 &
        !                    FILE_NAME   = Me%Files%ResultsHDF,                          &
        !                    Message     = Message,                                      &
        !                    TIME_END    = Me%ExtVar%EndTime,                            &
        !                    Extension   = 'hdf5',                                       &
        !                    FilesInput  = Me%Files%Nomfich,                             &
        !                    STAT        = STAT_CALL)                           
        !if (STAT_CALL /= SUCCESS_) stop 'ConstructFilesNames - ModuleLitter - ERR40'
        !
        !i = len_trim(Me%Files%ResultsHDF)
        !
        !if (Me%Files%ResultsHDF(i:i) /="5") then
        !    Me%Files%ResultsHDF = trim(Me%Files%ResultsHDF)//"5"
        !endif

        !Initial Particle Litter List HDF File
        Message   ='Initial list of litter particles in HDF format.'
        call ReadFileName  (KEYWORD     = 'LITTER_INI',                                 &
                            FILE_NAME   = Me%Files%LitterIni,                           &
                            Message     = Message,                                      &
                            TIME_END    = Me%ExtVar%StartTime,                          &
                            Extension   = 'hdf5',                                       &
                            FilesInput  = Me%Files%Nomfich,                             &
                            STAT        = STAT_CALL)                           
        if (STAT_CALL /= SUCCESS_) stop 'ConstructFilesNames - ModuleLitter - ERR50'
        
        !Final Particle Litter List HDF File
        Message   ='Final list of litter particles in HDF format.'
        call ReadFileName  (KEYWORD     = 'LITTER_FIN',                                 &
                            FILE_NAME   = Me%Files%LitterFin,                           &
                            Message     = Message,                                      &
                            TIME_END    = Me%ExtVar%EndTime,                            &
                            Extension   = 'hdf5',                                       &
                            FilesInput  = Me%Files%Nomfich,                             &
                            STAT        = STAT_CALL)                           
        if (STAT_CALL /= SUCCESS_) stop 'ConstructFilesNames - ModuleLitter - ERR60'
        

        
        !----------------------------------------------------------------------

    end subroutine ConstructFilesNames

    !--------------------------------------------------------------------------    
    
    subroutine ConstructFromLitterBlock()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: BlockFound
        !----------------------------------------------------------------------

                                                                                
        call ExtractBlockFromBuffer(EnterDataID     = Me%ObjEnterData,                  &
                                    ClientNumber    = Me%ClientNumber,                  &
                                    block_begin     = '<BeginLitter>',                  &
                                    block_end       = '<EndLitter>',                    &
                                    BlockFound      = BlockFound,                       &
                                    STAT            = STAT_CALL)         
        if (STAT_CALL /= SUCCESS_) stop 'ConstructFromLitterBlock - ModuleLitter - ERR10'

BF:     if (BlockFound) then

            call ConstructGlobalOptions
            
            call ReadBeachAreas
            
            call ReadOutputGrids

        else
        
            stop 'ConstructFromLitterBlock - ModuleLitter - ERR20'
          
        endif BF

            
    end subroutine ConstructFromLitterBlock
    
    !------------------------------------------------------------------------------    

    subroutine ConstructGlobalOptions()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, flag
        !----------------------------------------------------------------------
    
       call GetData(value           = Me%AgeToBeach,                                    &
                    EnterDataID     = Me%ObjEnterData,                                  &
                    flag            = flag,                                             &
                    SearchType      = FromBlock,                                        &
                    keyword         = 'AGE_TO_BEACH',                                   &
                    default         = 0.,                                               &
                    ClientModule    = 'ModuleLitter',                                   &
                    STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalOptions - ModuleLitter - ERR010'


       call GetData(value           = Me%KillBeachLitter,                               &
                    EnterDataID     = Me%ObjEnterData,                                  &
                    flag            = flag,                                             &
                    SearchType      = FromBlock,                                        &
                    keyword         = 'KILL_BEACH_LITTER',                              &
                    default         = .true.,                                           &
                    ClientModule    = 'ModuleLitter',                                   &
                    STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalOptions - ModuleLitter - ERR020'        
        
        Me%OutputGrids%WriteRestartFile = .false. 

        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime = Me%ExtVar%CurrentTime,                         &
                           EndTime     = Me%ExtVar%EndTime,                             &
                           keyword     = 'RESTART_FILE_OUTPUT_TIME',                    &
                           SearchType  = FromBlock,                                     &
                           OutPutsTime = Me%OutputGrids%RestartOutTime,                 &
                           OutPutsOn   = Me%OutputGrids%WriteRestartFile,               &
                           STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalOptions - ModuleLitter - ERR030'

        if(Me%OutputGrids%WriteRestartFile)then

            Me%OutputGrids%NextRestartOutput = 1

        end if 

        !Checks wether to overwrite the Restart File OR not
        call GetData(Me%OutputGrids%RestartOverwrite,                                   &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromBlock,                                          &
                     keyword      ='RESTART_FILE_OVERWRITE',                            &
                     ClientModule ='ModuleLitter',                                      &
                     Default      = .true.,                                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalOptions - ModuleLitter - ERR040'
        
        
    end subroutine ConstructGlobalOptions
    
    !--------------------------------------------------------------------------


    subroutine ReadBeachAreas

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------


        call CountIndividualBeachAreas

        call ReadIndividualBeachAreas 
            

    end subroutine ReadBeachAreas

    !--------------------------------------------------------------------------

    subroutine CountIndividualBeachAreas()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, NAreas
        logical                                         :: AreasFound

        !Begin-----------------------------------------------------------------

        NAreas = 0

DOPROP: do 

            call ExtractBlockFromBlock (EnterDataID         = Me%ObjEnterData,          &
                                        ClientNumber        = Me%ClientNumber,          &
                                        block_begin         = '<<BeginBeachArea>>',     &
                                        block_end           = '<<EndBeachArea>>',       &
                                        BlockInBlockFound   = AreasFound,               &
                                        STAT                = STAT_CALL)         
            if (STAT_CALL /= SUCCESS_) stop 'CountIndividualBeachAreas - ModuleLitter - ERR10'
            
i1:         if (AreasFound) then

                NAreas = NAreas + 1
 
            else i1
            
                call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'CountIndividualBeachAreas - ModuleLitter - ERR20'
                exit
            endif i1

        enddo DOPROP
        
        Me%BeachAreas%Number = NAreas

        allocate(Me%BeachAreas%Individual(Me%BeachAreas%Number))

    end subroutine CountIndividualBeachAreas
    !--------------------------------------------------------------------------

    subroutine ReadIndividualBeachAreas()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                         :: AreasFound
        integer                                         :: STAT_CALL, nAreas, flag
        !Begin-----------------------------------------------------------------
        
DONB:   do nAreas = 1, Me%BeachAreas%Number
 
            call ExtractBlockFromBlock (EnterDataID         = Me%ObjEnterData,          &
                                        ClientNumber        = Me%ClientNumber,          &
                                        block_begin         = '<<BeginBeachArea>>',     &
                                        block_end           = '<<EndBeachArea>>',       &
                                        BlockInBlockFound   = AreasFound,               &
                                        STAT                = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualBeachAreas - ModuleLitter - ERR10'
            
i1:         if (AreasFound) then

                call GetData(Me%BeachAreas%Individual(nAreas)%Name,                     &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='NAME',                                      &
                             default      = 'Areas X',                                  &
                             ClientModule ='ModuleLitter',                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualBeachAreas - ModuleLitter - ERR20'
 
                call GetData(Me%BeachAreas%Individual(nAreas)%Description,              &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='DESCRIPTION',                               &
                             default      = 'Areas to contain a oil spill',             &
                             ClientModule ='ModuleLitter',                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualBeachAreas - ModuleLitter - ERR30'
 
                call GetData(Me%BeachAreas%Individual(nAreas)%VelThreshold,             &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='VEL_THRESHOLD',                             &
                             default      = 0.4,                                        &
                             ClientModule ='ModuleLitter',                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualBeachAreas - ModuleLitter - ERR40'

                call GetData(Me%BeachAreas%Individual(nAreas)%WaveThreshold,            &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='WAVE_THRESHOLD',                            &
                             default      = 0.6,                                        &
                             ClientModule ='ModuleLitter',                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualBeachAreas - ModuleLitter - ERR50'

                call GetData(Me%BeachAreas%Individual(nAreas)%Probability,              &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='PROBABILITY',                               &
                             !100% probability
                             default      = 1.0,                                        &
                             ClientModule ='ModuleLitter',                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualBeachAreas - ModuleLitter - ERR60'                
                
               call GetData(Me%BeachAreas%Individual(nAreas)%FileName,                  &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='FILENAME',                                  &
                             ClientModule ='ModuleLitter',                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualBeachAreas - ModuleLitter - ERR70'
                
                if (flag == 0) then
                    write(*,*) 'Areas named ', trim(Me%BeachAreas%Individual(nAreas)%Name),' needs a ".xy" file'
                    stop 'ReadIndividualBeachAreas - ModuleLitter - ERR80'                
                endif
                
                nullify(Me%BeachAreas%Individual(nAreas)%Polygons)
                
                !Generates beach areas from file -> t_polygons that intersect the model domain
                call New(Polygons           = Me%BeachAreas%Individual(nAreas)%Polygons,&
                         PolygonsFileName   = Me%BeachAreas%Individual(nAreas)%FileName,& 
                         PolygonsRef        = Me%ExtVar%ModelDomain)                                
                

                call GetData(Me%BeachAreas%Individual(nAreas)%CoastType,                &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='COAST_TYPE',                                &
                             default      = 1,                                          &
                             ClientModule ='ModuleLitter',                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualBeachAreas - ModuleLitter - ERR90'
                
                Me%BeachAreas%Individual(nAreas)%ID = nAreas           
                
                
                !   AGE_LIMIT       : seconds                     
                !   [Age limit after beaching above which the particle is delete from module litter]                                           
                call GetData(Me%BeachAreas%Individual(nAreas)%AgeLimit,                 &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='AGE_LIMIT',                                 &
                             default      = - null_real,                                &
                             ClientModule ='ModuleLitter',                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualBeachAreas - ModuleLitter - ERR100'
                
                Me%BeachAreas%Individual(nAreas)%ID = nAreas        
                
               
            else i1
            
                stop 'ReadIndividualBeachAreas - ModuleLitter - ERR110'

            endif i1
            
        enddo DONB
        


    end subroutine ReadIndividualBeachAreas
                
    !--------------------------------------------------------------------------

    subroutine ReadOutputGrids

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------


        call CountIndividualOutputGrids

        call ReadIndividualOutputGrids 
            

    end subroutine ReadOutputGrids

    !--------------------------------------------------------------------------

    subroutine CountIndividualOutputGrids()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, NGrids
        logical                                         :: GridsFound

        !Begin-----------------------------------------------------------------

        NGrids = 0

DOPROP: do 

            call ExtractBlockFromBlock (EnterDataID         = Me%ObjEnterData,          &
                                        ClientNumber        = Me%ClientNumber,          &
                                        block_begin         = '<<BeginOutputGrid>>',    &
                                        block_end           = '<<EndOutputGrid>>',      &
                                        BlockInBlockFound   = GridsFound,               &
                                        STAT                = STAT_CALL)         
            if (STAT_CALL /= SUCCESS_) stop 'CountIndividualOutputGrids - ModuleLitter - ERR10'
            
i1:         if (GridsFound) then

                NGrids = NGrids + 1
 
            else i1
            
                call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'CountIndividualOutputGrids - ModuleLitter - ERR20'
                exit
            endif i1

        enddo DOPROP
        
        Me%OutputGrids%Number = NGrids

        allocate(Me%OutputGrids%Individual(Me%OutputGrids%Number))

    end subroutine CountIndividualOutputGrids
    !--------------------------------------------------------------------------

    subroutine ReadIndividualOutputGrids()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,       dimension(:,:), pointer             :: Bathymetry 
        integer,    dimension(:,:), pointer             :: WaterPoints2D    
        type (T_Size2D)                                 :: Size, WorkSize
        logical                                         :: GridsFound, OutputON
        integer                                         :: STAT_CALL, nGrids, flag, HDF5_CREATE
        !Begin-----------------------------------------------------------------
        
DONB:   do nGrids = 1, Me%OutputGrids%Number
 
            call ExtractBlockFromBlock (EnterDataID         = Me%ObjEnterData,          &
                                        ClientNumber        = Me%ClientNumber,          &
                                        block_begin         = '<<BeginOutputGrid>>',    &
                                        block_end           = '<<EndOutputGrid>>',      &
                                        BlockInBlockFound   = GridsFound,               &
                                        STAT                = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_) then
                stop 'ReadIndividualOutputGrids - ModuleLitter - ERR10'
            endif                
            
i1:         if (GridsFound) then

                call GetOutPutTime(Me%ObjEnterData,                                     &
                                   CurrentTime   = Me%ExtVar%CurrentTime,               &
                                   EndTime       = Me%ExtVar%EndTime,                   &
                                   keyword       = 'OUTPUT_TIME',                       &
                                   SearchType    = FromBlockInBlock,                    &
                                   OutPutsTime   = Me%OutputGrids%Individual(nGrids)%OutTime,&
                                   OutPutsOn     = OutputON,                            &
                                   OutPutsNumber = Me%OutputGrids%Individual(nGrids)%Number, &
                                   STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'ReadIndividualOutputGrids - ModuleLitter - ERR20'
                endif  

                    
                if (.not. OutputON) then 
                    stop 'ReadIndividualOutputGrids - ModuleLitter - ERR30'
                endif
                
                Me%OutputGrids%Individual(nGrids)%NextOutPut = 1

                call GetData(Me%OutputGrids%Individual(nGrids)%OutPutFile,              &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='OUTPUT_FILENAME',                           &
                             ClientModule ='ModuleLitter',                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleLitter - ERR40'
                if (flag      ==        0) stop 'ReadIndividualOutputGrids - ModuleLitter - ERR50'
                
                !Gets File Access Code
                call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
    
                !Opens HDF File
                call ConstructHDF5(HDF5ID   = Me%OutputGrids%Individual(nGrids)%ObjHDF5,&
                                   FileName = trim(Me%OutputGrids%Individual(nGrids)%OutPutFile),&
                                   Access   = HDF5_CREATE,                              &
                                   STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleLitter - ERR60'
 
                call GetData(Me%OutputGrids%Individual(nGrids)%InputGridFile,           &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='INPUT_GRID_FILENAME',                       &
                             ClientModule ='ModuleLitter',                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleLitter - ERR70'
                
                if (flag      ==        0) then
                    stop 'ReadIndividualOutputGrids - ModuleLitter - ERR80'
                endif
                
                call ConstructHorizontalGrid(HorizontalGridID = Me%OutputGrids%Individual(nGrids)%ObjHorizontalGrid,& 
                                             DataFile         = Me%OutputGrids%Individual(nGrids)%InputGridFile,    &
                                             STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleLitter - ERR90'  
                

                call WriteHorizontalGrid (HorizontalGridID    = Me%OutputGrids%Individual(nGrids)%ObjHorizontalGrid, &
                                          ObjHDF5             = Me%OutputGrids%Individual(nGrids)%ObjHDF5, &
                                          STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleLitter - ERR100'
                
                !Gets the grid size 
                call GetHorizontalGridSize(HorizontalGridID = Me%OutputGrids%Individual(nGrids)%ObjHorizontalGrid,  &
                                           Size             = Size,                                                 &
                                           WorkSize         = WorkSize,                                             &
                                           STAT             = STAT_CALL)     
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleLitter - ERR110'
                                           
                Me%OutputGrids%Individual(nGrids)%Size      = Size
                Me%OutputGrids%Individual(nGrids)%WorkSize  = WorkSize                                         
                
                allocate(Bathymetry   (Size%ILB:Size%IUB, Size%JLB:Size%JUB))
                allocate(WaterPoints2D(Size%ILB:Size%IUB, Size%JLB:Size%JUB))                
                

                call HDF5SetLimits(HDF5ID   = Me%OutputGrids%Individual(nGrids)%ObjHDF5,&
                                   ILB      = WorkSize%ILB,                             &
                                   IUB      = WorkSize%IUB,                             &
                                   JLB      = WorkSize%JLB,                             &
                                   JUB      = WorkSize%JUB,                             &
                                   STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleLitter - ERR120'
                
                Bathymetry(:,:) = 0.

                call HDF5WriteData(HDF5ID   = Me%OutputGrids%Individual(nGrids)%ObjHDF5,&
                                   GroupName= "/Grid",                                  &
                                   Name     = "Bathymetry",                             &
                                   Units    = "m",                                      &                           
                                   Array2D  = Bathymetry,                               &
                                   STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR130'                

                WaterPoints2D(:,:) = 1

                call HDF5WriteData(HDF5ID   = Me%OutputGrids%Individual(nGrids)%ObjHDF5,&
                                   GroupName= "/Grid",                                  &
                                   Name     = "WaterPoints",                            &
                                   Units    = "-",                                      &                           
                                   Array2D  = WaterPoints2D,                            &
                                   STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR130'                
                
                deallocate(Bathymetry   )
                deallocate(WaterPoints2D)                

                allocate(Me%OutputGrids%Individual(nGrids)%AuxReal2D(Size%ILB:Size%IUB, Size%JLB:Size%JUB))     
                
            else i1
            
                stop 'ReadIndividualOutputGrids - ModuleLitter - ERR200'

            endif i1
            
        enddo DONB
        


    end subroutine ReadIndividualOutputGrids
               
        
    !--------------------------------------------------------------------------
        
    type(T_Time) function HDF5TimeInstant(Instant, HDF5ID)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        integer                                 :: HDF5ID
        

        !Local-----------------------------------------------------------------
!        type(T_Time)                            :: TimeInstant
        real,    dimension(:), pointer          :: TimeVector
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        allocate(TimeVector(6))

        call HDF5SetLimits  (HDF5ID, 1, 6, STAT = STAT_CALL)        

        call HDF5ReadWindow (HDF5ID         = HDF5ID,                                   &
                             GroupName      = "/Time",                                  &
                             Name           = "Time",                                   &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleField4D - ERR01'

        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))

                                     
        deallocate(TimeVector)

    end function HDF5TimeInstant

!--------------------------------------------------------------------------    
        

    subroutine AllocateNewParticle (NewPartic, nP, nArea)

        !Arguments-------------------------------------------------------------
        type (T_Particle), pointer                  :: NewPartic
        integer                                     :: nP, nArea
        
        
        !Local-----------------------------------------------------------------
        
        !Begin-----------------------------------------------------------------        

        nullify  (NewPartic)
        allocate (NewPartic)
        nullify  (NewPartic%Next)
        nullify  (NewPartic%Prev)
        
        NewPartic%BeachTime     = Me%ExtVar%CurrentTime
        
        NewPartic%ID            = Me%ExtVar%ID         (nP)
        NewPartic%Longitude     = Me%ExtVar%Longitude  (nP)
        NewPartic%Latitude      = Me%ExtVar%Latitude   (nP)
        NewPartic%Age           = Me%ExtVar%Age        (nP)
        NewPartic%Origin        = Me%ExtVar%Origin     (nP) 

        NewPartic%BeachAreaID   = Me%BeachAreas%Individual(nArea)%ID
        NewPartic%CoastType     = Me%BeachAreas%Individual(nArea)%CoastType       
        NewPartic%BeachPeriod   = 0

    end subroutine AllocateNewParticle

    !--------------------------------------------------------------------------

    subroutine InsertParticleToBeachList (NewParticle)

        !Arguments-------------------------------------------------------------
        type (T_Particle), pointer                      :: NewParticle

        !Local-----------------------------------------------------------------
        type (T_Particle), pointer                      :: CurrentParticle  => null()
        type (T_Particle), pointer                      :: PreviousParticle => null()

        !Begin-----------------------------------------------------------------

        !Inserts a new property to the list of properties
        if (.not. associated(Me%ParticleList%First)) then
            Me%ParticleList%Number =  0           
            Me%ParticleList%First  => NewParticle
        else
            PreviousParticle => Me%ParticleList%First
            CurrentParticle  => PreviousParticle%Next
            do while (associated(CurrentParticle))
                PreviousParticle => CurrentParticle
                CurrentParticle  => PreviousParticle%Next
            enddo
            PreviousParticle%Next => NewParticle
            NewParticle%Prev      => PreviousParticle
        endif
        

        Me%ParticleList%Number = Me%ParticleList%Number + 1


    end subroutine InsertParticleToBeachList

    !--------------------------------------------------------------------------

    subroutine DeleteParticle (ParticleToDelete, NextParticleOut)

        !Arguments-------------------------------------------------------------
        type (T_Particle), pointer                    :: ParticleToDelete
        type (T_Particle), pointer                    :: NextParticleOut

        !Local-----------------------------------------------------------------
        type (T_Particle), pointer                    :: CurrentPartic => null()
        type (T_Particle), pointer                    :: NextParticle  => null()
        type (T_Particle), pointer                    :: PrevParticle  => null()

        logical                                       :: ParticleDeleted
        integer                                       :: nP
        
        !Begin-----------------------------------------------------------------

        nP = ParticleID(ParticleToDelete)         
        
        if (nP > 0) then
            Me%ExtVar%KillPartic(nP) = .true.
        endif

        ParticleDeleted = .false. 

        CurrentPartic => Me%ParticleList%First

d1:     do while (associated(CurrentPartic)) 
                     
i1:         if (CurrentPartic%ID == ParticleToDelete%ID) then


                PrevParticle => CurrentPartic%Prev
                NextParticle => CurrentPartic%Next

 
                !Updates foward pointer
                if (associated(CurrentPartic%Prev)) then
                    if (associated(NextParticle)) then
                        PrevParticle%Next => NextParticle
                    else
                        nullify(PrevParticle%Next)
                    endif
                endif

                !Updates backward pointer
                if (associated(CurrentPartic%Next)) then
                    if (associated(PrevParticle)) then
                        NextParticle%Prev => PrevParticle
                    else
                        nullify(NextParticle%Prev)
                    endif
                endif

                if (ParticleToDelete%ID == Me%ParticleList%First%ID) then
                    Me%ParticleList%First => NextParticle
                endif

                NextParticleOut => NextParticle

                !Deallocate Particle
                nullify       (ParticleToDelete%Next)
                nullify       (ParticleToDelete%Prev)
                deallocate    (ParticleToDelete)
                nullify       (ParticleToDelete, CurrentPartic)

                !Decreases number of particle
                Me%ParticleList%Number = Me%ParticleList%Number - 1

                if (Me%ParticleList%Number == 0) then
                    nullify (Me%ParticleList%First) 
                endif

                ParticleDeleted = .true. 

                exit

            endif i1

            CurrentPartic => CurrentPartic%Next

        enddo d1

        if (.not. ParticleDeleted) then
            stop 'DeleteParticle - ModuleLitter - ERR10'
        endif
        
    end subroutine DeleteParticle    
    
   
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyLitterWater(ObjLitterID,                                           &
                            nParticles,                                                 &
                            CurrentTime,                                                &
                                 NextCompute,                                           &
                            Longitude,                                                  &
                            Latitude,                                                   &                
                            Age,                                                        &
                            Origin,                                                     &
                            ID,                                                         &
                            Beach,                                                      &
                            KillPartic,                                                 &    
                            STAT)        
        !Arguments-------------------------------------------------------------
        integer                       , intent(IN)      :: ObjLitterID
        integer                       , intent(IN)      :: nParticles
        type (T_Time)                 , intent(IN)      :: CurrentTime  
        type (T_Time)                 , intent(IN)      :: NextCompute
        real(8), dimension(:), pointer, intent(IN)      :: Longitude
        real(8), dimension(:), pointer, intent(IN)      :: Latitude
        real(8), dimension(:), pointer, intent(IN)      :: Age     
        integer, dimension(:), pointer, intent(IN)      :: Origin             
        integer, dimension(:), pointer, intent(IN)      :: ID        
        logical, dimension(:), pointer, intent(INOUT)   :: Beach
        logical, dimension(:), pointer, intent(INOUT)   :: KillPartic                
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjLitterID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
                                    
            Me%ExtVar%nParticles    =  nParticles
            Me%ExtVar%CurrentTime   =  CurrentTime
            Me%ExtVar%NextCompute   =  NextCompute
            Me%ExtVar%Longitude     => Longitude
            Me%ExtVar%Latitude      => Latitude
            Me%ExtVar%Age           => Age     
            Me%ExtVar%Origin        => Origin     
            Me%ExtVar%ID            => ID
            Me%ExtVar%Beach         => Beach
            Me%ExtVar%KillPartic    => KillPartic      
            
            call CheckBeachLitter
            
            call DeleteOldLitter
            
            call OutputNumberGrid
            
            call OutputRestartFile
            
            Me%LastAtualization = Me%ExtVar%CurrentTime

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyLitterWater


    !--------------------------------------------------------------------------
    
    subroutine ModifyLitterLag(ObjLitterID,                                           &
                                 CurrentTime,                                           &
                                 Longitude,                                             &
                                 Latitude,                                              &                
                                 Age,                                                   &
                                 Origin,                                                &
                                 ID,                                                    &
                                 Beach,                                                 &
                                 KillPartic,                                            &    
                                 STAT)        
        !Arguments-------------------------------------------------------------
        integer                       , intent(IN)      :: ObjLitterID
        integer, dimension(6)         , intent(IN)      :: CurrentTime  
        real(8), dimension(:), pointer, intent(IN)      :: Longitude
        real(8), dimension(:), pointer, intent(IN)      :: Latitude
        real(8), dimension(:), pointer, intent(IN)      :: Age     
        integer, dimension(:), pointer, intent(IN)      :: Origin             
        integer, dimension(:), pointer, intent(IN)      :: ID        
        logical, dimension(:), pointer, intent(INOUT)   :: Beach
        logical, dimension(:), pointer, intent(INOUT)   :: KillPartic                
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjLitterID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
                                    
            Me%ExtVar%nParticles    =  size(Longitude)
            
            call SetDate (Time1 = Me%ExtVar%CurrentTime,                                &
                          Year  = CurrentTime(1),                                       &
                          Month = CurrentTime(2),                                       &
                          Day   = CurrentTime(3),                                       &
                          Hour  = CurrentTime(4),                                       &
                          Minute= CurrentTime(5),                                       &
                          Second= CurrentTime(6))
            
            Me%ExtVar%Longitude     => Longitude
            Me%ExtVar%Latitude      => Latitude
            Me%ExtVar%Age           => Age     
            Me%ExtVar%Origin        => Origin     
            Me%ExtVar%ID            => ID
            Me%ExtVar%Beach         => Beach
            Me%ExtVar%KillPartic    => KillPartic      
            
            call CheckBeachLitter
            
            call DeleteOldLitter
            
            call OutputNumberGrid
            
            call OutputRestartFile
            
            Me%LastAtualization = Me%ExtVar%CurrentTime

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyLitterLag


    !--------------------------------------------------------------------------

    subroutine CheckBeachLitter()        
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Particle), pointer                  :: NewPartic
        type (T_PointF),    pointer                 :: Point   
        real                                        :: Rand1
        integer                                     :: nArea, nP, nAreaTotal, nPtotal
        
        !----------------------------------------------------------------------
        
        nPtotal     = Me%ExtVar%nParticles
        nAreaTotal  = Me%BeachAreas%Number
        
        allocate(Point)
        
d1:     do nP    = 1, nPtotal
            
i0:         if (.not. Me%ExtVar%KillPartic(nP)) then
                if (Me%KillBeachLitter) then
                Me%ExtVar%Beach     (nP) = .false.
                endif

i1:             if (Me%ExtVar%Age (nP) > Me%AgeToBeach .and. .not. Me%ExtVar%Beach(nP)) then

d2:                 do nArea = 1, nAreaTotal

                        Point%X = Me%ExtVar%Longitude(nP)
                        Point%Y = Me%ExtVar%Latitude (nP)

i2:                    if (IsVisible(Me%BeachAreas%Individual(nArea)%Polygons, Point)) then

                            call RANDOM_NUMBER(Rand1)
                                
i3:                         if (Me%BeachAreas%Individual(nArea)%Probability >= Rand1) then  
                    
                                Me%ExtVar%Beach(nP) = .true.
                    
                                call AllocateNewParticle       (NewPartic, nP, nArea)
                                call InsertParticleToBeachList (NewPartic)                    
                    
i4:                             if (Me%KillBeachLitter) then
                                    Me%ExtVar%KillPartic(nP) = .true.
                                endif i4
                            endif i3
                                                                
                        endif i2
                    enddo   d2
                endif i1
            endif i0
        enddo   d1
        
        deallocate(Point)

    end subroutine CheckBeachLitter


    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    integer function ParticleID(CurrentPartic)        
        !Arguments-------------------------------------------------------------
        type (T_Particle), pointer                  :: CurrentPartic

        !Local-----------------------------------------------------------------
        integer                                     :: nP, nPtotal
        
        !----------------------------------------------------------------------
        
        ParticleID  = FillValueInt
        
        nPtotal     = Me%ExtVar%nParticles
        
d1:     do nP    = 1, nPtotal
    
i1:         if (Me%ExtVar%Origin(nP) == CurrentPartic%Origin .and.                      &
                Me%ExtVar%ID    (nP) == CurrentPartic%ID    ) then
                ParticleID = nP
                exit
            endif i1

        enddo   d1
        

    end function ParticleID


    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------

    subroutine DeleteOldLitter()        
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Particle), pointer                  :: DeletePartic
        type (T_Particle), pointer                  :: CurrentPartic        
        real                                        :: DT
        integer                                     :: nArea
        
        !----------------------------------------------------------------------
        
        DT = Me%ExtVar%CurrentTime - Me%LastAtualization
        
        CurrentPartic => Me%ParticleList%First

d1:     do while (associated(CurrentPartic)) 

            
            nArea = CurrentPartic%BeachAreaID

            !Check if the particle need to be killed                                    
i2:         if (CurrentPartic%BeachPeriod >= Me%BeachAreas%Individual(nArea)%AgeLimit) then

                DeletePartic   => CurrentPartic
                
                call DeleteParticle(ParticleToDelete = DeletePartic, NextParticleOut = CurrentPartic)
                            
            else i2
            
                !Atualize the particles beach age if not created in this instant (BeachPeriod = 0)            
                CurrentPartic%BeachPeriod = CurrentPartic%BeachPeriod + DT

                CurrentPartic => CurrentPartic%Next

            endif i2                
                    
        enddo   d1
        


    end subroutine DeleteOldLitter


    !--------------------------------------------------------------------------

    
    !--------------------------------------------------------------------------

    subroutine OutputNumberGrid()        
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Particle), pointer                  :: CurrentPartic
        real                                        :: PointX, PointY
        integer                                     :: nGrid, nTotalGrids
        integer                                     :: iOut, STAT_CALL, I, J, Nout
        logical                                     :: HaveDomain
        type (T_Time)                               :: Aux, OutTime
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr      
        real                                        :: TotalTime, AuxPeriod
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !----------------------------------------------------------------------
        
        nTotalGrids = Me%OutputGrids%Number
        
d1:     do nGrid = 1, nTotalGrids        

            iOut = Me%OutputGrids%Individual(nGrid)%NextOutput
            
            OutTime = Me%OutputGrids%Individual(nGrid)%OutTime(iOut)
            Nout    = size(Me%OutputGrids%Individual(nGrid)%OutTime)
            
i1:         if (Me%ExtVar%CurrentTime >=  OutTime .and. iOut <= Nout) then
    
                if (Me%ExtVar%Backtracking) then
                    iOut = Me%OutputGrids%Individual(nGrid)%Number - iOut + 1 
                endif          

                Me%OutputGrids%Individual(nGrid)%AuxReal2D(:,:) = 0
        
                CurrentPartic => Me%ParticleList%First

d2:             do while (associated(CurrentPartic)) 

                    PointX = CurrentPartic%Longitude
                    PointY = CurrentPartic%Latitude
                
                    HaveDomain = GetXYInsideDomain(Me%OutputGrids%Individual(nGrid)%ObjHorizontalGrid,  &
                                                   PointX,                              &
                                                   PointY,                              &
                                                   Referential= GridCoord_,             &
                                                   STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutputNumberGrid - ModuleLitter - ERR10'
                
                    if (HaveDomain) then

                        call GetXYCellZ(HorizontalGridID = Me%OutputGrids%Individual(nGrid)%ObjHorizontalGrid,&
                                        XPoint           = PointX,                      &
                                        YPoint           = PointY,                      &
                                        I                = I,                           &
                                        J                = J,                           &
                                        Referential      = GridCoord_,                  &
                                        STAT             = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutputNumberGrid - ModuleLitter - ERR20'

                        Me%OutputGrids%Individual(nGrid)%AuxReal2D(I,J) = Me%OutputGrids%Individual(nGrid)%AuxReal2D(I,J) + 1

                    endif                

                    CurrentPartic => CurrentPartic%Next
                    
                enddo   d2

                if (Me%ExtVar%Backtracking) then
                    
                    TotalTime = Me%ExtVar%EndTime - Me%ExtVar%StartTime                  
                    AuxPeriod = OutTime           - Me%ExtVar%StartTime
                    AuxPeriod = TotalTime         - AuxPeriod
                    
                    Aux = Me%ExtVar%StartTime   + AuxPeriod
                    
                else
                    
                    Aux = OutTime
                    
                endif 
                
                !Aux = Me%OutputGrids%Individual(nGrid)%OutTime(iOut)
                
                !Writes the Instant - HDF 5
                call ExtractDate   (Aux, AuxTime(1), AuxTime(2), AuxTime(3),          &
                                    AuxTime(4), AuxTime(5), AuxTime(6))
                TimePtr => AuxTime

                call HDF5SetLimits(HDF5ID   = Me%OutputGrids%Individual(nGrid)%ObjHDF5, &
                                   ILB      = 1,                                        &
                                   IUB      = 6,                                        &
                                   STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutputNumberGrid - ModuleLitter - ERR30'

                call HDF5WriteData(HDF5ID       = Me%OutputGrids%Individual(nGrid)%ObjHDF5, &
                                   GroupName    = "/Time",                              &
                                   Name         = "Time",                               &
                                   Units        =  "YYYY/MM/DD HH:MM:SS",               &                           
                                   Array1D      = TimePtr,                              &
                                   OutPutNumber = iOut,                                 &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutputNumberGrid - ModuleLitter - ERR40'                
                
                call HDF5SetLimits(HDF5ID   = Me%OutputGrids%Individual(nGrid)%ObjHDF5,     &
                                   ILB      = Me%OutputGrids%Individual(nGrid)%WorkSize%ILB,&
                                   IUB      = Me%OutputGrids%Individual(nGrid)%WorkSize%IUB,&
                                   JLB      = Me%OutputGrids%Individual(nGrid)%WorkSize%JLB,&
                                   JUB      = Me%OutputGrids%Individual(nGrid)%WorkSize%JUB,&
                                   STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutputNumberGrid - ModuleLitter - ERR50'

                
                call HDF5WriteData(HDF5ID       = Me%OutputGrids%Individual(nGrid)%ObjHDF5, &
                                   GroupName    = "/Results/Number",                        &
                                   Name         = "Number",                                 &
                                   Units        = "-",                                      &                           
                                   Array2D      = Me%OutputGrids%Individual(nGrid)%AuxReal2D,&
                                   OutPutNumber = iOut,                                     &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutputNumberGrid - ModuleLitter - ERR60'
                
                Me%OutputGrids%Individual(nGrid)%NextOutput = Me%OutputGrids%Individual(nGrid)%NextOutput + 1                
                
                
            endif i1                
        enddo   d1
        
    end subroutine OutputNumberGrid


    !--------------------------------------------------------------------------
    
    subroutine OutputRestartFile
        
        !Local-----------------------------------------------------------------
        real                        :: Year, Month, Day, Hour, Minute, Second
        integer                     :: i
    
        !Begin-----------------------------------------------------------------
        
        if (Me%OutputGrids%WriteRestartFile) then
            
            i = Me%OutputGrids%NextRestartOutput
            
            if  (Me%ExtVar%CurrentTime >= Me%OutputGrids%RestartOutTime(i)) then
    
                call ExtractDate(Me%ExtVar%CurrentTime,                                 &
                                 Year = Year, Month  = Month,  Day    = Day,            &
                                 Hour = Hour, Minute = Minute, Second = Second)
    
                call WriteLitterBeach(Final =.false.)
    
                Me%OutputGrids%NextRestartOutput = Me%OutputGrids%NextRestartOutput + 1
    
                call SetError(WARNING_, INTERNAL_, "Litter restart file saved      : ", &
                              Year, Month, Day, Hour, Minute, Second)
            endif                
        endif
    
    end subroutine OutputRestartFile
    
    !--------------------------------------------------------------------------    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillLitter(ObjLitterID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjLitterID              
        integer, optional, intent(OUT)      :: STAT

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           
        integer                             :: ready_              

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjLitterID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            call WriteLitterBeach()

            call DeallocateVariables

            nUsers = DeassociateInstance(mLitter_,  Me%InstanceID)

            if (nUsers == 0) then
                
                !Deallocates Instance
                call DeallocateInstance ()

                ObjLitterID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillLitter
        

    !------------------------------------------------------------------------
    
    subroutine WriteLitterBeach(Final)

        !Arguments---------------------------------------------------------------
        logical , optional                      :: Final
        !Local-------------------------------------------------------------------
        character(len=PathLength)               :: Filename
        character(len=StringLength)             :: AuxChar
        logical                                 :: Final_
        integer                             :: STAT_CALL, HDF5_CREATE, HDF5_READWRITE

        !------------------------------------------------------------------------

        !Write the particle list in a specific file to make more clear the hot start process
        if  (present(Final)) then
            Final_ = Final
        else
            Final_ = .true.
        endif     
        

        !Writes final results in the lagrangian transient results. 
        !This way the user can do quick figures using the lagrangian model grid
        
        Me%ObjHDF5 = 0
        
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READWRITE = HDF5_READWRITE)
    
        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%ResultsHDF), HDF5_READWRITE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteLitterBeach - ModuleLitter - ERR10'


        if (Me%ParticleList%Number > 0) then
            
            if (Final_) then
                
            call WriteAllParticlesBeach(WriteTime = .false.)

            else
                
                AuxChar = trim(TimeToString(Me%ExtVar%CurrentTime))
                call WriteAllParticlesBeach(WriteTime = .false., CharDate = AuxChar)

            endif             
            
        endif            
        
        !Closes HDF File
        call KillHDF5      (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteLitterBeach - ModuleLitter - ERR20'
        
        !Checks if it's at the end of the run 
        !or !if it's supposed to overwrite the final HDF file
        if (Me%OutputGrids%RestartOverwrite .or. Final_) then

            Filename = trim(Me%Files%LitterFin)

        else

            Filename =  ChangeSuffix(Me%Files%LitterFin,                         &
                                "_"//Trim(AuxChar)//".fin")

        endif     
        
        Me%ObjHDF5 = 0

        if (Me%ParticleList%Number > 0) then
            
            
            
            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)            
            
            !Opens HDF File
            call ConstructHDF5      (Me%ObjHDF5, trim(Filename), HDF5_CREATE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteLitterBeach - ModuleLitter - ERR30'            
            
            call WriteAllParticlesBeach(WriteTime = .True.)
            
            !Closes HDF File
            call KillHDF5      (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteLitterBeach - ModuleLitter - ERR40'     
        
        endif            
        
   

    end subroutine WriteLitterBeach

    !------------------------------------------------------------------------
    
    subroutine ReadLitterBeach()

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                             :: STAT_CALL, HDF5_READ
        logical                             :: FileExists
        !------------------------------------------------------------------------
        
        inquire (FILE=trim(Me%Files%LitterIni), EXIST = FileExists)        
        
        if (FileExists) then
            
            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
            
            Me%ObjHDF5 = 0

            !Opens HDF File
            call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%LitterIni), HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLitterBeach - ModuleLitter - ERR30'

            call ReadAllParticlesBeach
        
            !Closes HDF File
            call KillHDF5      (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLitterBeach - ModuleLitter - ERR40'        
            
        endif            

    end subroutine ReadLitterBeach
    
    !------------------------------------------------------------------------    
    
    subroutine WriteAllParticlesBeach(WriteTime, CharDate)

        !Arguments---------------------------------------------------------------
        logical                                 :: WriteTime
        character(len=*), optional              :: CharDate

        !Local-------------------------------------------------------------------
        integer     , dimension(:),   pointer   :: ID          
        real        , dimension(:,:), pointer   :: BeachTime
        real(8)     , dimension(:),   pointer   :: Longitude   
        real(8)     , dimension(:),   pointer   :: Latitude    
        real(8)     , dimension(:),   pointer   :: Age         
        integer     , dimension(:),   pointer   :: Origin      
        integer     , dimension(:),   pointer   :: BeachAreaID 
        integer     , dimension(:),   pointer   :: CoastType   
        real(8)     , dimension(:),   pointer   :: BeachPeriod
        type (T_Particle), pointer              :: CurrentPartic
        integer                                 :: STAT_CALL, nP, nPtotal
        character (len = StringLength)          :: GroupName, PropertyName, UnitsName
        real        , dimension(:),   pointer   :: Time6

        !------------------------------------------------------------------------

        nPtotal = Me%ParticleList%Number
        
        allocate(ID         (1:nPtotal))  
        allocate(BeachTime  (1:nPtotal,1:6))
        allocate(Longitude  (1:nPtotal))
        allocate(Latitude   (1:nPtotal))
        allocate(Age        (1:nPtotal))
        allocate(Origin     (1:nPtotal))
        allocate(BeachAreaID(1:nPtotal))
        allocate(CoastType  (1:nPtotal))         
        allocate(BeachPeriod(1:nPtotal))                 
         
        !Output matrixes
        CurrentPartic   => Me%ParticleList%First
        nP = 0
        do while (associated(CurrentPartic))
            nP = nP + 1
            
            Longitude   (nP)= CurrentPartic%Longitude
            Latitude    (nP)= CurrentPartic%Latitude
            Age         (nP)= CurrentPartic%Age
            ID          (nP)= CurrentPartic%ID
            Origin      (nP)= CurrentPartic%Origin
            BeachAreaID (nP)= CurrentPartic%BeachAreaID
            CoastType   (nP)= CurrentPartic%CoastType
            BeachPeriod (nP)= CurrentPartic%BeachPeriod
            
            call ExtractDate   (Time1   = CurrentPartic%BeachTime,                      &
                                Year    = BeachTime(nP, 1),                             &
                                Month   = BeachTime(nP, 2),                             &
                                Day     = BeachTime(nP, 3),                             &
                                Hour    = BeachTime(nP, 4),                             & 
                                Minute  = BeachTime(nP, 5),                             & 
                                Second  = BeachTime(nP, 6)) 
            
            CurrentPartic => CurrentPartic%Next
        enddo
        
        if (WriteTime) then
    
            !Read time 
            call HDF5SetLimits (HDF5ID  = Me%ObjHDF5,                                       &
                                ILB     = 1,                                                &
                                IUB     = 6,                                                &
                                STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR10'

            GroupName    = "/"
            PropertyName = "Time"    
            UnitsName    = "YY,MM,DD,HH,MM,SS"
            
            allocate(Time6(1:6))
            
            call ExtractDate   (Time1   = Me%ExtVar%CurrentTime,                            &
                                Year    = Time6(1),                                         &
                                Month   = Time6(2),                                         &
                                Day     = Time6(3),                                         &
                                Hour    = Time6(4),                                         & 
                                Minute  = Time6(5),                                         & 
                                Second  = Time6(6))                 
        
            call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                                GroupName   = trim(GroupName)//trim(PropertyName),          &
                                Name        = trim(PropertyName),                           & 
                                Units       = trim(UnitsName),                              &             
                                Array1D     = Time6,                                        &
                                STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR20'
            
            deallocate(Time6)
        
        endif        
        
        
        if (present(CharDate)) then
            GroupName = "/Results/BeachLitter/Restart/"//trim(CharDate)//"/"
        else
        GroupName = "/Results/BeachLitter/"
        endif
        
        !Output matrixes 1D - HDF5       
        call HDF5SetLimits (HDF5ID  = Me%ObjHDF5,                                       &
                            ILB     = 1,                                                &
                            IUB     = Me%ParticleList%Number,                           &
                            STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR30'

        PropertyName = "Longitude"    
        UnitsName    = "o"
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = Longitude,                                    &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR40'

        PropertyName = "Latitude"    
        UnitsName    = "o"
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = Latitude,                                     &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR50'

        PropertyName = "Particle_ID"    
        UnitsName    = "-"        
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = ID,                                           &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR60'    
        

        PropertyName = "Age"    
        UnitsName    = "seconds"        
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = Age,                                          &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR70'          
        

        PropertyName = "Origin_ID"    
        UnitsName    = "-"        
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = Origin,                                       &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR80'          

        PropertyName = "BeachArea_ID"    
        UnitsName    = "-"        
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = BeachAreaID,                                  &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR90'                  

        PropertyName = "Coast_Type"    
        UnitsName    = "-"        
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = CoastType,                                    &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR100'

        PropertyName = "Beach_Period"    
        UnitsName    = "seconds"           
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = BeachPeriod,                                  &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR110'        
        
        !Output matrixes 1D - HDF5       
        call HDF5SetLimits (HDF5ID  = Me%ObjHDF5,                                       &
                            ILB     = 1,                                                &
                            IUB     = Me%ParticleList%Number,                           &
                            JLB     = 1,                                                &
                            JUB     = 6,                                                &
                            STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR120'

        PropertyName = "Beach_Time"    
        UnitsName    = "o"
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array2D     = BeachTime,                                    &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR130'        
        
        deallocate(ID         )  
        deallocate(BeachTime  )
        deallocate(Longitude  )
        deallocate(Latitude   )
        deallocate(Age        )
        deallocate(Origin     )
        deallocate(BeachAreaID)
        deallocate(CoastType  )   
        deallocate(BeachPeriod)           
        

    end subroutine WriteAllParticlesBeach
    
    
    !------------------------------------------------------------------------    
    
    subroutine ReadAllParticlesBeach()

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer     , dimension(:),   pointer   :: ID          
        real        , dimension(:,:), pointer   :: BeachTime
        real(8)     , dimension(:),   pointer   :: Longitude   
        real(8)     , dimension(:),   pointer   :: Latitude    
        real(8)     , dimension(:),   pointer   :: Age         
        integer     , dimension(:),   pointer   :: Origin      
        integer     , dimension(:),   pointer   :: BeachAreaID 
        integer     , dimension(:),   pointer   :: CoastType   
        real(8)     , dimension(:),   pointer   :: BeachPeriod
        type (T_Particle), pointer              :: CurrentPartic
        integer                                 :: STAT_CALL, nP, nPtotal
        character (len = StringLength)          :: GroupName, PropertyName, UnitsName
        logical                                 :: GroupExists
        real        , dimension(:),   pointer   :: Time6
        type (T_Time)                           :: TimeAux

        !------------------------------------------------------------------------
        
        GroupName = "/Results/BeachLitter/"
        
        
        call GetHDF5GroupExist (HDF5ID      = Me%ObjHDF5,                               &
                                GroupName   = trim(GroupName),                          &
                                Exist       = GroupExists,                              & 
                                STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR10'        
        
ex:     if (GroupExists) then
    
            !Read time 
            call HDF5SetLimits (HDF5ID  = Me%ObjHDF5,                                   &
                                ILB     = 1,                                            &
                                IUB     = 6,                                            &
                                STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR20'

            GroupName    = "/"
            PropertyName = "Time"                
            UnitsName    = "o"
            
            allocate(Time6(1:6))
        
            call HDF5ReadData  (HDF5ID      = Me%ObjHDF5,                               &
                                GroupName   = trim(GroupName)//trim(PropertyName),      &
                                Name        = trim(PropertyName),                       & 
                                Array1D     = Time6,                                    &
                                STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR30'
            
            call SetDate   (Time1   = TimeAux,                                          &
                            Year    = Time6(1),                                         &
                            Month   = Time6(2),                                         &
                            Day     = Time6(3),                                         &
                            Hour    = Time6(4),                                         & 
                            Minute  = Time6(5),                                         & 
                            Second  = Time6(6))          
            
            if (TimeAux /= Me%ExtVar%StartTime) then
                write(*,*) 'Litter fin file date different from the start time of the run'
                write(*,*) 'Litter fin date ', ConvertTimeToString(TimeAux, ":")
                write(*,*) 'Run start time  ', ConvertTimeToString(Me%ExtVar%StartTime, ":")                
                stop 'ReadAllParticlesBeach - ModuleLitter - ERR40'
            endif                
            
            deallocate(Time6)

            
            GroupName = "/Results/BeachLitter/"        
            PropertyName = "Longitude"    
            UnitsName    = "o"
        
            call GetHDF5ArrayDimensions(HDF5ID    = Me%ObjHDF5,                         &
                                        GroupName = trim(GroupName)//trim(PropertyName),&
                                        ItemName  = trim(PropertyName),                 &
                                        Imax      = Me%ParticleList%Number,             & 
                                        STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR20'        

            nPtotal = Me%ParticleList%Number
        
            allocate(ID         (1:nPtotal))  
            allocate(BeachTime  (1:nPtotal,1:6))
            allocate(Longitude  (1:nPtotal))
            allocate(Latitude   (1:nPtotal))
            allocate(Age        (1:nPtotal))
            allocate(Origin     (1:nPtotal))
            allocate(BeachAreaID(1:nPtotal))
            allocate(CoastType  (1:nPtotal))         
            allocate(BeachPeriod(1:nPtotal))     
            

            !Output matrixes 1D - HDF5       
            call HDF5SetLimits (HDF5ID  = Me%ObjHDF5,                                       &
                                ILB     = 1,                                                &
                                IUB     = Me%ParticleList%Number,                           &
                                STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR30'

            PropertyName = "Longitude"    
            UnitsName    = "o"
        
            call HDF5ReadData  (HDF5ID      = Me%ObjHDF5,                                   &
                                GroupName   = trim(GroupName)//trim(PropertyName),          &
                                Name        = trim(PropertyName),                           & 
                                Array1D     = Longitude,                                    &
                                STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR40'

            PropertyName = "Latitude"    
            UnitsName    = "o"
        
            call HDF5ReadData  (HDF5ID      = Me%ObjHDF5,                                   &
                                GroupName   = trim(GroupName)//trim(PropertyName),          &
                                Name        = trim(PropertyName),                           & 
                                Array1D     = Latitude,                                     &
                                STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR50'

            PropertyName = "Particle_ID"    
            UnitsName    = "-"        
        
            call HDF5ReadData  (HDF5ID      = Me%ObjHDF5,                                   &
                                GroupName   = trim(GroupName)//trim(PropertyName),          &
                                Name        = trim(PropertyName),                           & 
                                Array1D     = ID,                                           &
                                STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR60'    
        

            PropertyName = "Age"    
            UnitsName    = "seconds"        
        
            call HDF5ReadData  (HDF5ID      = Me%ObjHDF5,                                   &
                                GroupName   = trim(GroupName)//trim(PropertyName),          &
                                Name        = trim(PropertyName),                           & 
                                Array1D     = Age,                                          &
                                STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR70'          
        

            PropertyName = "Origin_ID"    
            UnitsName    = "-"        
        
            call HDF5ReadData  (HDF5ID      = Me%ObjHDF5,                                   &
                                GroupName   = trim(GroupName)//trim(PropertyName),          &
                                Name        = trim(PropertyName),                           & 
                                Array1D     = Origin,                                       &
                                STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR80'          

            PropertyName = "BeachArea_ID"    
            UnitsName    = "-"        
        
            call HDF5ReadData  (HDF5ID      = Me%ObjHDF5,                                   &
                                GroupName   = trim(GroupName)//trim(PropertyName),          &
                                Name        = trim(PropertyName),                           & 
                                Array1D     = BeachAreaID,                                  &
                                STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR90'                  

            PropertyName = "Coast_Type"    
            UnitsName    = "-"        
        
            call HDF5ReadData  (HDF5ID      = Me%ObjHDF5,                                   &
                                GroupName   = trim(GroupName)//trim(PropertyName),          &
                                Name        = trim(PropertyName),                           & 
                                Array1D     = CoastType,                                    &
                                STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR100'

            PropertyName = "Beach_Period"    
            UnitsName    = "seconds"           
        
            call HDF5ReadData  (HDF5ID      = Me%ObjHDF5,                                   &
                                GroupName   = trim(GroupName)//trim(PropertyName),          &
                                Name        = trim(PropertyName),                           & 
                                Array1D     = BeachPeriod,                                  &
                                STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR110'        
        
            !Output matrixes 1D - HDF5       
            call HDF5SetLimits (HDF5ID  = Me%ObjHDF5,                                       &
                                ILB     = 1,                                                &
                                IUB     = Me%ParticleList%Number,                           &
                                JLB     = 1,                                                &
                                JUB     = 6,                                                &
                                STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR120'

            PropertyName = "Beach_Time"    
            UnitsName    = "o"
        
            call HDF5ReadData (HDF5ID      = Me%ObjHDF5,                                    &
                                GroupName   = trim(GroupName)//trim(PropertyName),          &
                                Name        = trim(PropertyName),                           & 
                                Array2D     = BeachTime,                                    &
                                STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadAllParticlesBeach - ModuleLitter - ERR130'        
            
         
            do nP = 1, nPtotal
            
                nullify  (CurrentPartic)
                allocate (CurrentPartic)
                nullify  (CurrentPartic%Next)
                nullify  (CurrentPartic%Prev)            
            
                CurrentPartic%Longitude      = Longitude   (nP) 
                CurrentPartic%Latitude       = Latitude    (nP) 
                CurrentPartic%Age            = Age         (nP) 
                CurrentPartic%ID             = ID          (nP) 
                CurrentPartic%Origin         = Origin      (nP) 
                CurrentPartic%BeachAreaID    = BeachAreaID (nP) 
                CurrentPartic%CoastType      = CoastType   (nP) 
                CurrentPartic%BeachPeriod    = BeachPeriod (nP) 
            
                call SetDate   (Time1   = CurrentPartic%BeachTime,                          &
                                Year    = BeachTime(nP, 1),                                 &
                                Month   = BeachTime(nP, 2),                                 &
                                Day     = BeachTime(nP, 3),                                 &
                                Hour    = BeachTime(nP, 4),                                 & 
                                Minute  = BeachTime(nP, 5),                                 & 
                                Second  = BeachTime(nP, 6)) 

                call InsertParticleToBeachList (CurrentPartic)
        
            enddo        
            
            deallocate(ID         )  
            deallocate(BeachTime  )
            deallocate(Longitude  )
            deallocate(Latitude   )
            deallocate(Age        )
            deallocate(Origin     )
            deallocate(BeachAreaID)
            deallocate(CoastType  )   
            deallocate(BeachPeriod)           
        
        endif ex
            
    end subroutine ReadAllParticlesBeach
    
    
    !------------------------------------------------------------------------    
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Litter), pointer          :: AuxObjLitter
        type (T_Litter), pointer          :: PreviousObjLitter

        !Updates pointers
        if (Me%InstanceID == FirstObjLitter%InstanceID) then
            FirstObjLitter => FirstObjLitter%Next
        else
            PreviousObjLitter => FirstObjLitter
            AuxObjLitter      => FirstObjLitter%Next
            do while (AuxObjLitter%InstanceID /= Me%InstanceID)
                PreviousObjLitter => AuxObjLitter
                AuxObjLitter      => AuxObjLitter%Next
            enddo

            !Now update linked list
            PreviousObjLitter%Next => AuxObjLitter%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance
    
    !--------------------------------------------------------------------------    

    
    
    subroutine DeallocateVariables ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                         :: STAT_CALL, iG, nUsers
        !Begin-----------------------------------------------------------------        

        !Deallocates variables
        if (associated(Me%BeachAreas%Individual)) then
            deallocate(Me%BeachAreas%Individual)            
            nullify   (Me%BeachAreas%Individual)                    
        endif            
        
        do iG = 1, Me%OutputGrids%Number

            call KillHorizontalGrid (HorizontalGridID = Me%OutputGrids%Individual(iG)%ObjHorizontalGrid, &
                                        STAT             = STAT_CALL)      

            if (STAT_CALL /= SUCCESS_) then
                stop 'DeallocateVariables - ModuleLitter - ERR10'
            endif                                               
            
            call KillHDF5 (HDF5ID =  Me%OutputGrids%Individual(iG)%ObjHDF5,             &
                            STAT   = STAT_CALL)
                                         
            if (STAT_CALL /= SUCCESS_) then
                stop 'DeallocateVariables - ModuleLitter - ERR20'
            endif              

            deallocate(Me%OutputGrids%Individual(iG)%AuxReal2D)
            
        enddo            
        
        if (associated(Me%OutputGrids%Individual)) then            
            deallocate(Me%OutputGrids%Individual)            
            nullify   (Me%OutputGrids%Individual)                    
        endif           
        
        if (Me%ExtVar%ObjTime > 0) then
        
        nUsers = DeassociateInstance (mTIME_,   Me%ExtVar%ObjTime)
        if (nUsers == 0) stop 'DeallocateVariables - ModuleLitter - ERR30'        
        
        endif
            
    end subroutine DeallocateVariables
    

    
    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjLitter_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjLitter_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjLitter_ID > 0) then
            call LocateObjLitter (ObjLitter_ID)
            ready_ = VerifyReadLock (mLitter_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjLitter (ObjLitterID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjLitterID

        !Local-----------------------------------------------------------------

        Me => FirstObjLitter
        do while (associated (Me))
            if (Me%InstanceID == ObjLitterID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleLitter - LocateObjLitter - ERR01'

    end subroutine LocateObjLitter

    !--------------------------------------------------------------------------

end module ModuleLitter

