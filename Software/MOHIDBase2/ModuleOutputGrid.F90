!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : OutputGrid
! URL           : http://www.mohid.com
! AFFILIATION   : HIDROMOD
! DATE          : August 2020
! REVISION      : Paulo Leitao
! DESCRIPTION   : Module manages the output of grid results 
!
!------------------------------------------------------------------------------


Module ModuleOutputGrid

!<BeginOutputGrid>
!   AGE_TO_BEACH            : seconds                     [ OutputGrid age necessary to beach ] 
!   KILL_BEACH_OutputGrid       : 1(true)/0(false)            [By default (KILL_BEACH_OutputGrid : 1) 
!                                                          the OutputGrid particles removed from the lagrangian model]
!<EndOutputGrid>

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
    public  :: ConstructOutputGrid
    
    private ::      AllocateInstance

    !Selector
                     
    
    !Modifier
    public  :: ModifyOutputGrid

    
    private ::      OutputNumberGrid

  

    !Destructor
    public  :: KillOutputGrid

    !Management
    private ::      Ready
    private ::          LocateObjOutputGrid 
    
    !Interfaces----------------------------------------------------------------

    !Types---------------------------------------------------------------------
 
    
    !IO
    type T_Files
        character(PathLength)                           :: ConstructData        = null_str
        character(PathLength)                           :: block_begin          = null_str
        character(PathLength)                           :: block_end            = null_str
    end type T_Files    
    
    type T_ExtVar
        integer                                         :: ObjTime              = 0
        type(T_Time)                                    :: CurrentTime
        type(T_Time)                                    :: NextCompute
        type(T_Time)                                    :: StartTime                      
        type(T_Time)                                    :: EndTime              
        integer                                         :: nParticles           = null_int
        logical                                         :: Backtracking         = .false. 
        real(8),   dimension(:), pointer                :: Latitude             => null()
        real(8),   dimension(:), pointer                :: Longitude            => null()
    end type T_ExtVar        
                                                                                
    type T_OutPut
        type (T_Time), dimension(:), pointer            :: OutTime
        integer                                         :: NextOutPut           = null_int
        integer                                         :: Number               = null_int
        real,          dimension(:,:), pointer          :: AuxReal2D            => null()
        character (len = PathLength)                    :: OutputFile           = null_str
        character (len = PathLength)                    :: InputGridFile        = null_str        
        integer                                         :: ObjHDF5              = 0
        integer                                         :: ObjHDF5_2            = 0        
        integer                                         :: ObjHorizontalGrid    = 0 
        type (T_Size2D)                                 :: Size
        type (T_Size2D)                                 :: WorkSize
    end type T_OutPut    
    

    type T_OutputGrid
        type (T_Output), dimension(:), pointer          :: Individual           => null()           
        integer                                         :: Number               = null_int
        integer                                         :: InstanceID           = null_int
        
        integer                                         :: FromWathBlock        = null_int 
        
        logical                                         :: OutputGridON         = .false. 
        
        type (T_ExtVar)                                 :: ExtVar        
        
        type (T_Time)                                   :: LastAtualization

        type (T_Files)                                  :: Files

        integer                                         :: ClientNumber         = null_int
        integer                                         :: ObjEnterdata         = 0
        integer                                         :: ObjHDF5              = 0        
        
        type(T_OutputGrid), pointer                     :: Next                 => null()
        
    end type T_OutputGrid
    

    !Global Module Variables
    type (T_OutputGrid), pointer                            :: FirstObjOutputGrid       => null()
    type (T_OutputGrid), pointer                            :: Me                   => null()


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructOutputGrid(ObjOutputGridID, TimeID, ConstructData,          &
                                   block_begin, block_end, FromWathBlock,           &
                                   OutputGridON, STAT)

        !Arguments---------------------------------------------------------------
        integer            ,            intent(OUT)     :: ObjOutputGridID 
        integer            ,            intent(INOUT)   :: TimeID 
        character(len=*)   ,            intent(IN)      :: ConstructData
        character(len=*)   ,            intent(IN)      :: block_begin
        character(len=*)   ,            intent(IN)      :: block_end  
        integer            ,            intent(IN)      :: FromWathBlock
        logical            ,            intent(OUT)     :: OutputGridON
        integer, optional,              intent(OUT)     :: STAT     

        !Local-------------------------------------------------------------------
        integer                                         :: ready_         
        integer                                         :: STAT_
        integer                                         :: STAT_CALL        
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mOutputGrid_)) then
            nullify (FirstObjOutputGrid)
            call RegisterModule (mOutputGrid_) 
        endif

        call Ready(ObjOutputGridID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            Me%Files%ConstructData  = ConstructData 
            Me%Files%block_begin    = block_begin
            Me%Files%block_end      = block_end
            
            Me%FromWathBlock        = FromWathBlock
           
            Me%ExtVar%ObjTime       = AssociateInstance (mTIME_, TimeID)
            
            call GetExternalTime
            
            !Construct enter data 
            call ConstructEnterData(EnterDataID     = Me%ObjEnterData,                  &
                                    FileName        = Me%Files%ConstructData,           &
                                    ErrorMessage    = "ConstructOutputGrid - ModuleOutputGrid", &
                                    STAT            = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOutputGrid - ModuleOutputGrid - ERR10'
            
            
            
            call ConstructFromOutputGridBlock
            
           
            !Kill enter data 
            call KillEnterData     (EnterDataID     = Me%ObjEnterData,                  &
                                    STAT            = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOutputGrid - ModuleOutputGrid - ERR20'
            
            OutputGridON = Me%OutputGridON
                        
            !Returns ID
            ObjOutputGridID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleOutputGrid - ConstructOutputGrid - ERR100' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructOutputGrid
 
    !--------------------------------------------------------------------------
                                    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_OutputGrid), pointer                         :: NewObjOutputGrid
        type (T_OutputGrid), pointer                         :: PreviousObjOutputGrid


        !Allocates new instance
        allocate (NewObjOutputGrid)
        nullify  (NewObjOutputGrid%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjOutputGrid)) then
            FirstObjOutputGrid         => NewObjOutputGrid
            Me                    => NewObjOutputGrid
        else
            PreviousObjOutputGrid      => FirstObjOutputGrid
            Me                    => FirstObjOutputGrid%Next
            do while (associated(Me))
                PreviousObjOutputGrid  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjOutputGrid
            PreviousObjOutputGrid%Next => NewObjOutputGrid
        endif

        Me%InstanceID = RegisterNewInstance (mOutputGrid_)


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
        if (STAT_CALL /= SUCCESS_) stop 'GetExternalTime - ModuleOutputGrid - ERR10'

        ! Check if the simulation goes backward in time or forward in time (default mode)
        call GetBackTracking(Me%ExtVar%ObjTime,                                         &
                             Me%ExtVar%BackTracking, STAT = STAT_CALL)                    
        if (STAT_CALL /= SUCCESS_) stop 'GetExternalTime - ModuleOutputGrid - ERR20'

        Me%ExtVar%CurrentTime   = Me%ExtVar%StartTime                        
        Me%LastAtualization     = Me%ExtVar%StartTime    
            
    end subroutine GetExternalTime
    

    !--------------------------------------------------------------------------    
    
    subroutine ConstructFromOutputGridBlock()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

            
        call ReadOutputGrids
            
    end subroutine ConstructFromOutputGridBlock
    
    !------------------------------------------------------------------------------    

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
        logical                                         :: BlockFound

        !Begin-----------------------------------------------------------------

        NGrids = 0

DOPROP: do 

            call ExtractBlockFromBuffer(EnterDataID   = Me%ObjEnterData,                &
                                        ClientNumber  = Me%ClientNumber,                &
                                        block_begin   = trim(Me%Files%block_begin),     &
                                        block_end     = trim(Me%Files%block_end),       &
                                        BlockFound    = BlockFound,                     &
                                        STAT          = STAT_CALL)         
            if (STAT_CALL /= SUCCESS_) stop 'CountIndividualOutputGrids - ModuleOutputGrid - ERR10'
            
i1:         if (BlockFound) then

                NGrids = NGrids + 1
                
                Me%OutputGridON = .true. 
 
            else i1
            
                call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'CountIndividualOutputGrids - ModuleOutputGrid - ERR20'
                exit
            endif i1

        enddo DOPROP
        
        Me%Number = NGrids

        allocate(Me%Individual(Me%Number))

    end subroutine CountIndividualOutputGrids
    !--------------------------------------------------------------------------

    subroutine ReadIndividualOutputGrids()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,       dimension(:,:), pointer             :: Bathymetry 
        integer,    dimension(:,:), pointer             :: WaterPoints2D    
        type (T_Size2D)                                 :: Size, WorkSize
        logical                                         :: BlockFound, OutputON
        integer                                         :: STAT_CALL, nGrids, flag, HDF5_CREATE
        !Begin-----------------------------------------------------------------
        
DONB:   do nGrids = 1, Me%Number
 
            call ExtractBlockFromBuffer(EnterDataID         = Me%ObjEnterData,          &
                                        ClientNumber        = Me%ClientNumber,          &
                                        block_begin         = trim(Me%Files%block_begin),&
                                        block_end           = trim(Me%Files%block_end), &
                                        BlockFound          = BlockFound,               &
                                        STAT                = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) then
                stop 'ReadIndividualOutputGrids - ModuleOutputGrid - ERR10'
            endif                
            
i1:         if (BlockFound) then

                call GetOutPutTime(Me%ObjEnterData,                                     &
                                   CurrentTime   = Me%ExtVar%CurrentTime,               &
                                   EndTime       = Me%ExtVar%EndTime,                   &
                                   keyword       = 'OUTPUT_TIME',                       &
                                   SearchType    = Me%FromWathBlock,                    &
                                   OutPutsTime   = Me%Individual(nGrids)%OutTime,&
                                   OutPutsOn     = OutputON,                            &
                                   OutPutsNumber = Me%Individual(nGrids)%Number, &
                                   STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'ReadIndividualOutputGrids - ModuleOutputGrid - ERR20'
                endif  

                    
                if (.not. OutputON) then 
                    stop 'ReadIndividualOutputGrids - ModuleOutputGrid - ERR30'
                endif
                
                Me%Individual(nGrids)%NextOutPut = 1

                call GetData(Me%Individual(nGrids)%OutPutFile,              &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = Me%FromWathBlock,                           &
                             keyword      ='OUTPUT_FILENAME',                           &
                             ClientModule ='ModuleOutputGrid',                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleOutputGrid - ERR40'
                if (flag      ==        0) stop 'ReadIndividualOutputGrids - ModuleOutputGrid - ERR50'
                
                !Gets File Access Code
                call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
    
                !Opens HDF File
                call ConstructHDF5(HDF5ID   = Me%Individual(nGrids)%ObjHDF5,            &
                                   FileName = trim(Me%Individual(nGrids)%OutPutFile),   &
                                   Access   = HDF5_CREATE,                              &
                                   STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleOutputGrid - ERR60'
 
                call GetData(Me%Individual(nGrids)%InputGridFile,                       &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = Me%FromWathBlock,                           &
                             keyword      ='INPUT_GRID_FILENAME',                       &
                             ClientModule ='ModuleOutputGrid',                          &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleOutputGrid - ERR70'
                
                if (flag      ==        0) then
                    stop 'ReadIndividualOutputGrids - ModuleOutputGrid - ERR80'
                endif
                
                call ConstructHorizontalGrid(HorizontalGridID = Me%Individual(nGrids)%ObjHorizontalGrid,& 
                                             DataFile         = Me%Individual(nGrids)%InputGridFile,    &
                                             STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleOutputGrid - ERR90'  
                

                call WriteHorizontalGrid (HorizontalGridID    = Me%Individual(nGrids)%ObjHorizontalGrid, &
                                          ObjHDF5             = Me%Individual(nGrids)%ObjHDF5, &
                                          STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleOutputGrid - ERR100'
                
                !Gets the grid size 
                call GetHorizontalGridSize(HorizontalGridID = Me%Individual(nGrids)%ObjHorizontalGrid,  &
                                           Size             = Size,                                                 &
                                           WorkSize         = WorkSize,                                             &
                                           STAT             = STAT_CALL)     
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleOutputGrid - ERR110'
                                           
                Me%Individual(nGrids)%Size      = Size
                Me%Individual(nGrids)%WorkSize  = WorkSize                                         
                
                allocate(Bathymetry   (Size%ILB:Size%IUB, Size%JLB:Size%JUB))
                allocate(WaterPoints2D(Size%ILB:Size%IUB, Size%JLB:Size%JUB))                
                

                call HDF5SetLimits(HDF5ID   = Me%Individual(nGrids)%ObjHDF5,&
                                   ILB      = WorkSize%ILB,                             &
                                   IUB      = WorkSize%IUB,                             &
                                   JLB      = WorkSize%JLB,                             &
                                   JUB      = WorkSize%JUB,                             &
                                   STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualOutputGrids - ModuleOutputGrid - ERR120'
                
                Bathymetry(:,:) = 0.

                call HDF5WriteData(HDF5ID   = Me%Individual(nGrids)%ObjHDF5,            &
                                   GroupName= "/Grid",                                  &
                                   Name     = "Bathymetry",                             &
                                   Units    = "m",                                      &                           
                                   Array2D  = Bathymetry,                               &
                                   STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR130'                

                WaterPoints2D(:,:) = 1

                call HDF5WriteData(HDF5ID   = Me%Individual(nGrids)%ObjHDF5,&
                                   GroupName= "/Grid",                                  &
                                   Name     = "WaterPoints",                            &
                                   Units    = "-",                                      &                           
                                   Array2D  = WaterPoints2D,                            &
                                   STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConvertGridDataToHDF5 - ERR130'                
                
                deallocate(Bathymetry   )
                deallocate(WaterPoints2D)                

                allocate(Me%Individual(nGrids)%AuxReal2D(Size%ILB:Size%IUB, Size%JLB:Size%JUB))     
                
            else i1
            
                stop 'ReadIndividualOutputGrids - ModuleOutputGrid - ERR200'

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

    subroutine ModifyOutputGrid(ObjOutputGridID,                                        &
                                nParticles,                                             &
                                CurrentTime,                                            &
                                Longitude,                                              &
                                Latitude,                                               &   
                                STAT)        
        !Arguments-------------------------------------------------------------
        integer                       , intent(IN)      :: ObjOutputGridID
        integer                       , intent(IN)      :: nParticles
        type (T_Time)                 , intent(IN)      :: CurrentTime
        real(8), dimension(:), pointer, intent(IN)      :: Longitude
        real(8), dimension(:), pointer, intent(IN)      :: Latitude
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjOutputGridID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
            
            Me%ExtVar%Longitude     => Longitude
            Me%ExtVar%Latitude      => Latitude
            Me%ExtVar%nParticles    =  nParticles
            Me%ExtVar%CurrentTime   =  CurrentTime
                                    
            call OutputNumberGrid
            
            Me%LastAtualization = Me%ExtVar%CurrentTime
            
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyOutputGrid

    !--------------------------------------------------------------------------

    subroutine OutputNumberGrid()        
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                                        :: PointX, PointY
        integer                                     :: nGrid, nTotalGrids
        integer                                     :: iOut, STAT_CALL, I, J, Nout, nP
        logical                                     :: HaveDomain
        type (T_Time)                               :: Aux, OutTime
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr      
        real                                        :: TotalTime, AuxPeriod
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !----------------------------------------------------------------------
        
        nTotalGrids = Me%Number
        
d1:     do nGrid = 1, nTotalGrids        

            iOut = Me%Individual(nGrid)%NextOutput
            
            OutTime = Me%Individual(nGrid)%OutTime(iOut)
            Nout    = size(Me%Individual(nGrid)%OutTime)
            
i1:         if (Me%ExtVar%CurrentTime >=  OutTime .and. iOut <= Nout) then
    
                if (Me%ExtVar%Backtracking) then
                    iOut = Me%Individual(nGrid)%Number - iOut + 1 
                endif          

                Me%Individual(nGrid)%AuxReal2D(:,:) = 0
        
d2:             do nP = 1, Me%ExtVar%nParticles

                    PointX = Me%ExtVar%Longitude(nP)
                    PointY = Me%ExtVar%Latitude(nP)
                
                    HaveDomain = GetXYInsideDomain(Me%Individual(nGrid)%ObjHorizontalGrid,  &
                                                   PointX,                              &
                                                   PointY,                              &
                                                   Referential= GridCoord_,             &
                                                   STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutputNumberGrid - ModuleOutputGrid - ERR10'
                
                    if (HaveDomain) then

                        call GetXYCellZ(HorizontalGridID = Me%Individual(nGrid)%ObjHorizontalGrid,&
                                        XPoint           = PointX,                      &
                                        YPoint           = PointY,                      &
                                        I                = I,                           &
                                        J                = J,                           &
                                        Referential      = GridCoord_,                  &
                                        STAT             = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutputNumberGrid - ModuleOutputGrid - ERR20'

                        Me%Individual(nGrid)%AuxReal2D(I,J) = Me%Individual(nGrid)%AuxReal2D(I,J) + 1

                    endif                
                    
                enddo   d2

                if (Me%ExtVar%Backtracking) then
                    
                    TotalTime = Me%ExtVar%EndTime - Me%ExtVar%StartTime                  
                    AuxPeriod = OutTime           - Me%ExtVar%StartTime
                    AuxPeriod = TotalTime         - AuxPeriod
                    
                    Aux = Me%ExtVar%StartTime   + AuxPeriod
                    
                else
                    
                    Aux = OutTime
                    
                endif 
                
                !Aux = Me%Individual(nGrid)%OutTime(iOut)
                
                !Writes the Instant - HDF 5
                call ExtractDate   (Aux, AuxTime(1), AuxTime(2), AuxTime(3),          &
                                    AuxTime(4), AuxTime(5), AuxTime(6))
                TimePtr => AuxTime

                call HDF5SetLimits(HDF5ID   = Me%Individual(nGrid)%ObjHDF5, &
                                   ILB      = 1,                                        &
                                   IUB      = 6,                                        &
                                   STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutputNumberGrid - ModuleOutputGrid - ERR30'

                call HDF5WriteData(HDF5ID       = Me%Individual(nGrid)%ObjHDF5, &
                                   GroupName    = "/Time",                              &
                                   Name         = "Time",                               &
                                   Units        =  "YYYY/MM/DD HH:MM:SS",               &                           
                                   Array1D      = TimePtr,                              &
                                   OutPutNumber = iOut,                                 &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutputNumberGrid - ModuleOutputGrid - ERR40'                
                
                call HDF5SetLimits(HDF5ID   = Me%Individual(nGrid)%ObjHDF5,     &
                                   ILB      = Me%Individual(nGrid)%WorkSize%ILB,&
                                   IUB      = Me%Individual(nGrid)%WorkSize%IUB,&
                                   JLB      = Me%Individual(nGrid)%WorkSize%JLB,&
                                   JUB      = Me%Individual(nGrid)%WorkSize%JUB,&
                                   STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutputNumberGrid - ModuleOutputGrid - ERR50'

                
                call HDF5WriteData(HDF5ID       = Me%Individual(nGrid)%ObjHDF5, &
                                   GroupName    = "/Results/Number",                        &
                                   Name         = "Number",                                 &
                                   Units        = "-",                                      &                           
                                   Array2D      = Me%Individual(nGrid)%AuxReal2D,&
                                   OutPutNumber = iOut,                                     &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutputNumberGrid - ModuleOutputGrid - ERR60'
                
                Me%Individual(nGrid)%NextOutput = Me%Individual(nGrid)%NextOutput + 1                
                
                
            endif i1                
        enddo   d1
        
    end subroutine OutputNumberGrid
    
    !--------------------------------------------------------------------------



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillOutputGrid(ObjOutputGridID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjOutputGridID              
        integer, optional, intent(OUT)      :: STAT

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           
        integer                             :: ready_              

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjOutputGridID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            call DeallocateVariables

            nUsers = DeassociateInstance(mOutputGrid_,  Me%InstanceID)

            if (nUsers == 0) then
                
                !Deallocates Instance
                call DeallocateInstance ()

                ObjOutputGridID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillOutputGrid
        
    !------------------------------------------------------------------------    
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_OutputGrid), pointer          :: AuxObjOutputGrid
        type (T_OutputGrid), pointer          :: PreviousObjOutputGrid

        !Updates pointers
        if (Me%InstanceID == FirstObjOutputGrid%InstanceID) then
            FirstObjOutputGrid => FirstObjOutputGrid%Next
        else
            PreviousObjOutputGrid => FirstObjOutputGrid
            AuxObjOutputGrid      => FirstObjOutputGrid%Next
            do while (AuxObjOutputGrid%InstanceID /= Me%InstanceID)
                PreviousObjOutputGrid => AuxObjOutputGrid
                AuxObjOutputGrid      => AuxObjOutputGrid%Next
            enddo

            !Now update linked list
            PreviousObjOutputGrid%Next => AuxObjOutputGrid%Next

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
        do iG = 1, Me%Number

            call KillHorizontalGrid (HorizontalGridID = Me%Individual(iG)%ObjHorizontalGrid, &
                                        STAT             = STAT_CALL)      

            if (STAT_CALL /= SUCCESS_) then
                stop 'DeallocateVariables - ModuleOutputGrid - ERR10'
            endif                                               
            
            call KillHDF5 (HDF5ID =  Me%Individual(iG)%ObjHDF5,             &
                            STAT   = STAT_CALL)
                                         
            if (STAT_CALL /= SUCCESS_) then
                stop 'DeallocateVariables - ModuleOutputGrid - ERR20'
            endif              

            deallocate(Me%Individual(iG)%AuxReal2D)
            
        enddo            
        
        if (associated(Me%Individual)) then            
            deallocate(Me%Individual)            
            nullify   (Me%Individual)                    
        endif 
        
        if (Me%ExtVar%ObjTime > 0) then
        
            nUsers = DeassociateInstance (mTIME_,   Me%ExtVar%ObjTime)
            if (nUsers == 0) stop 'DeallocateVariables - ModuleOutputGrid - ERR30'        
        
        endif
        
    end subroutine DeallocateVariables
    

    
    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjOutputGrid_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjOutputGrid_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjOutputGrid_ID > 0) then
            call LocateObjOutputGrid (ObjOutputGrid_ID)
            ready_ = VerifyReadLock (mOutputGrid_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjOutputGrid (ObjOutputGridID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjOutputGridID

        !Local-----------------------------------------------------------------

        Me => FirstObjOutputGrid
        do while (associated (Me))
            if (Me%InstanceID == ObjOutputGridID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleOutputGrid - LocateObjOutputGrid - ERR01'

    end subroutine LocateObjOutputGrid

    !--------------------------------------------------------------------------

end module ModuleOutputGrid

