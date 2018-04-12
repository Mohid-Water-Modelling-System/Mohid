!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Litter
! URL           : http://www.mohid.com
! AFFILIATION   : HIDROMOD & ??????
! DATE          : April 2018
! REVISION      : Paulo Leitao & ?????? 
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
!<<EndBeachArea>>
!<EndLitter>

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions    
    use ModuleDrawing
    use ModuleHDF5

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructLitter
    private ::      AllocateInstance

    !Selector
                     
    
    !Modifier
    public  :: ModifyLitter

    !Destructor
    public  :: KillLitter                                                     
    private ::      DeAllocateInstance

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
    end type T_Files    
    
    type T_ExtVar
        type(T_Time)                                    :: CurrentTime
        type(T_Time)                                    :: EndTime              
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

        type (T_Particle),     pointer                  :: Next                 => null()    
        type (T_Particle),     pointer                  :: Prev                 => null()            
    end type T_Particle        
    
    type T_ParticleList
        integer                                         :: Number               = null_int
        type (T_Particle), pointer                      :: First                => null() 
    end type T_ParticleList        
    
    private :: T_Litter
    type       T_Litter
        integer                                         :: InstanceID           = null_int
        
        type (T_ExtVar)                                 :: ExtVar        
        
        real                                            :: AgeToBeach           = null_real
        logical                                         :: KillBeachLitter      = .false.     

        type (T_Files)                                  :: Files
        type (T_BeachAreas)                             :: BeachAreas
        type (T_ParticleList)                           :: ParticleList

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

    subroutine ConstructLitter(ObjLitterID, Nomfich, EndTime, ModelDomain, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjLitterID 
        character(len=*)                                :: Nomfich
        type (T_Time)                                   :: EndTime
        type (T_Polygon), pointer                       :: ModelDomain
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
            Me%ExtVar%EndTime       = EndTime
            Me%ExtVar%ModelDomain   => ModelDomain
            
            call ConstructFilesNames
            
            !Construct enter data 
            call ConstructEnterData(EnterDataID     = Me%ObjEnterData,                  &
                                    FileName        = Me%Files%ConstructData,           &
                                    ErrorMessage    = "ConstructLitter - ModuleLitter", &
                                    STAT            = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructLitter - ModuleLitter - ERR10'
            
            
            call ConstructFromLitterBlock
            
           
            !Kill enter data 
            call KillEnterData     (EnterDataID     = Me%ObjEnterData,                  &
                                    STAT            = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructLitter - ModuleLitter - ERR20'
            
            !Returns ID
            ObjLitterID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleLitter - ConstructLitter - ERR100' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructLitter
 
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

    subroutine ConstructFilesNames

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, i
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
        
        
        !Transient HDF File
        Message   ='Instant fields of litter particles in HDF format.'
        call ReadFileName  (KEYWORD     = 'PARTIC_HDF',                                 &
                            FILE_NAME   = Me%Files%ResultsHDF,                          &
                            Message     = Message,                                      &
                            TIME_END    = Me%ExtVar%EndTime,                            &
                            Extension   = 'hdf5',                                       &
                            FilesInput  = Me%Files%Nomfich,                             &
                            STAT        = STAT_CALL)                           
        if (STAT_CALL /= SUCCESS_) stop 'ConstructFilesNames - ModuleLitter - ERR40'

        i = len_trim(Me%Files%ResultsHDF)
        
        if (Me%Files%ResultsHDF(i:i) /="5") then
            Me%Files%ResultsHDF = trim(Me%Files%ResultsHDF)//"5"
        endif

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
                call New(Me%BeachAreas%Individual(nAreas)%Polygons, Me%BeachAreas%Individual(nAreas)%FileName, Me%ExtVar%ModelDomain)                                
                

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
                
               
            else i1
            
                stop 'ReadIndividualBeachAreas - ModuleLitter - ERR100'

            endif i1
            
        enddo DONB
        


    end subroutine ReadIndividualBeachAreas
                
    !--------------------------------------------------------------------------

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

    subroutine DeleteParticle (ParticleToDelete)

        !Arguments-------------------------------------------------------------
        type (T_Particle), pointer                    :: ParticleToDelete

        !Local-----------------------------------------------------------------
        type (T_Particle), pointer                    :: CurrentPartic => null()
        type (T_Particle), pointer                    :: NextParticle  => null()
        type (T_Particle), pointer                    :: PrevParticle  => null()

        logical                                       :: ParticleDeleted
        !Begin-----------------------------------------------------------------

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


                !Deallocate Particle
!                nullify       (ParticleToDelete%Next)
!                nullify       (ParticleToDelete%Prev)
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

    subroutine ModifyLitter(ObjLitterID,                                                &
                            nParticles,                                                 &
                            CurrentTime,                                                &
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
            Me%ExtVar%Longitude     => Longitude
            Me%ExtVar%Latitude      => Latitude
            Me%ExtVar%Age           => Age     
            Me%ExtVar%Origin        => Origin     
            Me%ExtVar%ID            => ID
            Me%ExtVar%Beach         => Beach
            Me%ExtVar%KillPartic    => KillPartic      
            
            call CheckBeachLitter

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyLitter


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

            Me%ExtVar%Beach     (nP) = .false.
            Me%ExtVar%KillPartic(nP) = .false.

i1:         if (Me%ExtVar%Age (nP) > Me%AgeToBeach) then

d2:             do nArea = 1, nAreaTotal

                    Point%X = Me%ExtVar%Longitude(nP)
                    Point%Y = Me%ExtVar%Latitude (nP)

i2:                 if (IsVisible(Me%BeachAreas%Individual(nArea)%Polygons, Point)) then

                        call RANDOM_NUMBER(Rand1)
                                
i3:                     if (Me%BeachAreas%Individual(nArea)%Probability >= Rand1) then  
                    
                            Me%ExtVar%Beach(nP) = .true.
                    
                            call AllocateNewParticle       (NewPartic, nP, nArea)
                            call InsertParticleToBeachList (NewPartic)                    
                    
i4:                         if (Me%KillBeachLitter) then
                                Me%ExtVar%KillPartic(nP) = .true.
                            endif i4
                        endif i3
                                                                
                    endif i2
                enddo   d2
            endif i1                                
        enddo   d1
        
        deallocate(Point)

    end subroutine CheckBeachLitter


    !--------------------------------------------------------------------------
    
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

            call WriteLitterBeach

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
    
    subroutine WriteLitterBeach()

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                             :: STAT_CALL, HDF5_READWRITE

        !------------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READWRITE = HDF5_READWRITE)
    
        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%ResultsHDF), HDF5_READWRITE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteLitterBeach - ModuleLitter - ERR10'


        if (Me%ParticleList%Number > 0) then
            call WriteAllParticlesBeach
        endif            
        
        !Closes HDF File
        call KillHDF5      (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteLitterBeach - ModuleLitter - ERR20'
        

    end subroutine WriteLitterBeach

    !------------------------------------------------------------------------
    
    subroutine WriteAllParticlesBeach()

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
        type (T_Particle), pointer              :: CurrentPartic
        integer                                 :: STAT_CALL, nP, nPtotal
        character (len = StringLength)          :: GroupName, PropertyName, UnitsName

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
            
            call ExtractDate   (Time1   = CurrentPartic%BeachTime,                      &
                                Year    = BeachTime(nP, 1),                             &
                                Month   = BeachTime(nP, 2),                             &
                                Day     = BeachTime(nP, 3),                             &
                                Hour    = BeachTime(nP, 4),                             & 
                                Minute  = BeachTime(nP, 5),                             & 
                                Second  = BeachTime(nP, 6)) 
            
            CurrentPartic => CurrentPartic%Next
        enddo
        
        GroupName = "/Results/BeachLitter/"
        
        !Output matrixes 1D - HDF5       
        call HDF5SetLimits (HDF5ID  = Me%ObjHDF5,                                       &
                            ILB     = 1,                                                &
                            IUB     = Me%ParticleList%Number,                           &
                            STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR10'

        PropertyName = "Longitude"    
        UnitsName    = "o"
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = Longitude,                                    &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR20'

        PropertyName = "Latitude"    
        UnitsName    = "o"
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = Latitude,                                     &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR30'

        PropertyName = "Particle_ID"    
        UnitsName    = "-"        
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = ID,                                           &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR40'    
        

        PropertyName = "Age"    
        UnitsName    = "seconds"        
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = Age,                                          &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR50'          
        

        PropertyName = "Origin_ID"    
        UnitsName    = "-"        
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = Origin,                                       &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR60'          

        PropertyName = "BeachArea_ID"    
        UnitsName    = "-"        
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = BeachAreaID,                                  &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR70'                  

        PropertyName = "Coast_Type"    
        UnitsName    = "-"        
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array1D     = CoastType,                                    &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR80'
        
        !Output matrixes 1D - HDF5       
        call HDF5SetLimits (HDF5ID  = Me%ObjHDF5,                                       &
                            ILB     = 1,                                                &
                            IUB     = Me%ParticleList%Number,                           &
                            JLB     = 1,                                                &
                            JUB     = 6,                                                &
                            STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR90'

        PropertyName = "Beach_Time"    
        UnitsName    = "o"
        
        call HDF5WriteData (HDF5ID      = Me%ObjHDF5,                                   &
                            GroupName   = trim(GroupName)//trim(PropertyName),          &
                            Name        = trim(PropertyName),                           & 
                            Units       = trim(UnitsName),                              & 
                            Array2D     = BeachTime,                                    &
                            STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteAllParticlesBeach - ModuleLitter - ERR100'        
        
        deallocate(ID         )  
        deallocate(BeachTime  )
        deallocate(Longitude  )
        deallocate(Latitude   )
        deallocate(Age        )
        deallocate(Origin     )
        deallocate(BeachAreaID)
        deallocate(CoastType  )   
        

    end subroutine WriteAllParticlesBeach
    
    
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
        
        !Begin-----------------------------------------------------------------        

        !Deallocates variables
        deallocate(Me%BeachAreas%Individual)            
        nullify   (Me%BeachAreas%Individual)                    

            
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

