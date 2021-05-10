!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : MeshGlue
! URL           : http://www.mohid.com
! AFFILIATION   : Hidromod
! DATE          : December 2020
! REVISION      : Paulo Leitao - v1.0
! DESCRIPTION   : Module to serve as MeshGlue to create new modules
!
!------------------------------------------------------------------------------


Module ModuleMeshGlue

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleHorizontalGrid

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructMeshGlue
    private ::      AllocateInstance

    !Selector
                     
    
    !Modifier
    public  :: ModifyMeshGlue

    !Destructor
    public  :: KillMeshGlue                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjMeshGlue 
    
    !Interfaces----------------------------------------------------------------

    !Types---------------------------------------------------------------------
    private ::  T_Grid
    type        T_Grid
        ! Filename of the input grid data file
        character(len=PathLength)   :: FileName         = null_str
        !Instance of ModuleHorizontalGrid
        integer                     :: ObjHorizontalGrid= 0
        integer                     :: IO               = null_int
        integer                     :: JO               = null_int
        logical                     :: Flip_IJ          = .false.
        logical                     :: RowFlip          = .false.
        logical                     :: ColumnFlip       = .false.
        integer                     :: ILB              = null_int
        integer                     :: IUB              = null_int
        integer                     :: JLB              = null_int
        integer                     :: JUB              = null_int
    end type T_Grid

   
    private :: T_MeshGlue
    type       T_MeshGlue
        type (T_Grid), dimension(:), pointer            :: GridIn
        integer                                         :: GridIn_Number            = null_int        
        ! Filename of the input grid data file
        character(len=PathLength)                       :: FileNameGrid_Out_Start   = null_str
        character(len=PathLength)                       :: FileNameGrid_Out         = null_str        
        !Instance of ModuleHorizontalGrid
        integer                                         :: ObjHorizontalGridOut1    = 0
        integer                                         :: ObjHorizontalGridOut2    = 0
        
        logical                                         :: ResetGrid                = .false.  
        
        real(8), pointer, dimension(:,:)                :: XX_IE                    => null()
        real(8), pointer, dimension(:,:)                :: YY_IE                    => null()  
        
        type(T_Size2D)                                  :: Size, WorkSize
        
        integer                                         :: InstanceID
        !Instance of ObjEnterData
        integer                                         :: ObjEnterData   = 0
        
        type(T_MeshGlue), pointer                       :: Next
        
        
    end type  T_MeshGlue

    !Global Module Variables
    type (T_MeshGlue), pointer                          :: FirstObjMeshGlue
    type (T_MeshGlue), pointer                          :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructMeshGlue(ObjMeshGlueID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjMeshGlueID 
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mMeshGlue_)) then
            nullify (FirstObjMeshGlue)
            call RegisterModule (mMeshGlue_) 
        endif

        call Ready(ObjMeshGlueID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            !Returns ID
            ObjMeshGlueID          = Me%InstanceID
            
            call ReadGridIn
            
            call ReadGridOut            

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleMeshGlue - ConstructMeshGlue - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructMeshGlue
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_MeshGlue), pointer                         :: NewObjMeshGlue
        type (T_MeshGlue), pointer                         :: PreviousObjMeshGlue


        !Allocates new instance
        allocate (NewObjMeshGlue)
        nullify  (NewObjMeshGlue%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjMeshGlue)) then
            FirstObjMeshGlue         => NewObjMeshGlue
            Me                    => NewObjMeshGlue
        else
            PreviousObjMeshGlue      => FirstObjMeshGlue
            Me                    => FirstObjMeshGlue%Next
            do while (associated(Me))
                PreviousObjMeshGlue  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjMeshGlue
            PreviousObjMeshGlue%Next => NewObjMeshGlue
        endif

        Me%InstanceID = RegisterNewInstance (mMeshGlue_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------

    subroutine ReadGridIn

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL, iflag
        integer                                 :: ifl, ClientNumber
        logical                                 :: BlockFound

        !----------------------------------------------------------------------
        
        !Construct enter data
        call ConstructEnterData(Me%ObjEnterData, "MeshGlue.dat", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridIn - ModuleMeshGlue - ERR10'   
            
        call GetData(Me%GridIn_Number,                                                  &
                        Me%ObjEnterData, iflag,                                         &
                        Keyword        = 'GRID_IN_NUMBER',                              &
                        Default        = 1,                                             &
                        SearchType     = FromFile,                                      &
                        ClientModule   = 'ModuleMeshGlue',                              &
                        STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridIn - ModuleMeshGlue - ERR20'   
            
        allocate(Me%GridIn(Me%GridIn_Number))
            
do1 :   do ifl =1, Me%GridIn_Number
    
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                        ClientNumber    = ClientNumber,                 &
                                        block_begin     = "<begin_grid_in>",            &
                                        block_end       = "<end_grid_in>",              &
                                        BlockFound      = BlockFound,                   &
                                        STAT            = STAT_CALL)
cd1 :       if      (STAT_CALL == SUCCESS_     ) then
cd2 :           if (BlockFound) then
    
                    call ReadGridKeywords(FromBlock, Me%GridIn(ifl))
                        
                else cd2
    
                    stop 'ReadGridIn  - ModuleMeshGlue - ERR030'
                    
                endif cd2
                    
            else cd1 
                stop 'ReadGridIn  - ModuleMeshGlue - ERR040'
            endif cd1
                
        enddo do1

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)       
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridIn  - ModuleMeshGlue - ERR050'

    
    end subroutine ReadGridIn
    
    !--------------------------------------------------------------------------   
    
    subroutine ReadGridKeywords(FromWhere, Grid)
    
        !Arguments-------------------------------------------------------------
        integer                                 :: FromWhere
        type (T_Grid)                           :: Grid    

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL, iflag
        type(T_Size2D)                          :: WorkSize

        !----------------------------------------------------------------------
    
        !Read grid filemane of the input grid data file
        call GetData(Grid%FileName,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'FILENAME',                                       &
                     SearchType     = FromWhere,                                        &
                     ClientModule   = 'ModuleMeshGlue',                                 &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR10'

        if (iflag == 0) then
            stop 'ReadGridKeywords - ModuleMeshGlue - ERR20'
        endif
        
        call ConstructHorizontalGrid(HorizontalGridID = Grid%ObjHorizontalGrid,         &
                                     DataFile         = Grid%FileName,                  &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR30'


        call GetHorizontalGridSize(HorizontalGridID  = Grid%ObjHorizontalGrid,          &
                                    WorkSize         = WorkSize,                        &
                                    STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR40'            

        call GetData(Grid%IO,                                                           &
                        Me%ObjEnterData, iflag,                                         &
                        Keyword        = 'IO',                                          &
                        SearchType     = FromWhere,                                     &
                        ClientModule   = 'ModuleMeshGlue',                              &
                        STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR50'

        if (iflag /= 1) then
            stop 'ReadGridKeywords - ModuleMeshGlue - ERR60'
        endif

        call GetData(Grid%JO,                                                           &
                        Me%ObjEnterData, iflag,                                         &
                        Keyword        = 'JO',                                          &
                        SearchType     = FromWhere,                                     &
                        ClientModule   = 'ModuleMeshGlue',                              &
                        STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR70'

        if (iflag /= 1) then
            stop 'ReadGridKeywords - ModuleMeshGlue - ERR80'
        endif
        
        call GetData(Grid%Flip_IJ,                                                      &    
                        Me%ObjEnterData, iflag,                                         &
                        Keyword        = 'FLIP_IJ',                                     &
                        SearchType     = FromWhere,                                     &
                        default        = .false.,                                       &                    
                        ClientModule   = 'ModuleMeshGlue',                              &
                        STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR90'
    
        call GetData(Grid%RowFlip,                                                      &
                        Me%ObjEnterData, iflag,                                         &
                        Keyword        = 'ROW_FLIP',                                    &
                        SearchType     = FromWhere,                                     &
                        default        = .false.,                                       &
                        ClientModule   = 'ModuleMeshGlue',                              &
                        STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR100'
        
        call GetData(Grid%ColumnFlip,                                                   &
                        Me%ObjEnterData, iflag,                                         &
                        Keyword        = 'COLUMN_FLIP',                                 &
                        SearchType     = FromWhere,                                     &
                        default        = .false.,                                       &
                        ClientModule   = 'ModuleMeshGlue',                              &
                        STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR110'        
        
        call GetData(Grid%ILB,                                                          &
                        Me%ObjEnterData, iflag,                                         &
                        Keyword        = 'ILB',                                         &
                        SearchType     = FromWhere,                                     &
                        default        = WorkSize%ILB,                                  &
                        ClientModule   = 'ModuleMeshGlue',                              &
                        STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR120'
    
        call GetData(Grid%IUB,                                                          &
                        Me%ObjEnterData, iflag,                                         &
                        Keyword        = 'IUB',                                         &
                        SearchType     = FromWhere,                                     &
                        default        = WorkSize%IUB,                                  &
                        ClientModule   = 'ModuleMeshGlue',                              &
                        STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR130'            
        
        call GetData(Grid%JLB,                                                          &
                        Me%ObjEnterData, iflag,                                         &
                        Keyword        = 'JLB',                                         &
                        SearchType     = FromWhere,                                     &
                        default        = WorkSize%JLB,                                  &
                        ClientModule   = 'ModuleMeshGlue',                              &
                        STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR140'
        
        call GetData(Grid%JUB,                                                          &
                        Me%ObjEnterData, iflag,                                         &
                        Keyword        = 'JUB',                                         &
                        SearchType     = FromWhere,                                     &
                        default        = WorkSize%JUB,                                  &
                        ClientModule   = 'ModuleMeshGlue',                              &
                        STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR150'        
        
    
    end subroutine ReadGridKeywords
    
    !--------------------------------------------------------------------------    
            
    subroutine ReadGridOut
    
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL, iflag
        integer                                 :: ifl, ClientNumber
        logical                                 :: BlockFound

        !----------------------------------------------------------------------
  
    
        !Read grid filemane of the first iteration of the output grid
        call GetData(Me%FileNameGrid_Out_Start,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'FILENAME_OUT_START',                             &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleMeshGlue',                                 &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR10'

        if (iflag == 0) then
            stop 'ReadGridKeywords - ModuleMeshGlue - ERR20'
        endif
        
        call ConstructHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGridOut1,       &
                                     DataFile         = Me%FileNameGrid_Out_Start,      &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR30'        

        !Read grid filemane of the first iteration of the output grid
        call GetData(Me%FileNameGrid_Out,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'FILENAME_OUT',                                   &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleMeshGlue',                                 &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR40'

        if (iflag == 0) then
            stop 'ReadGridKeywords - ModuleMeshGlue - ERR50'
        endif

        call GetData(Me%ResetGrid,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'RESET_GRID',                                     &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleMeshGlue',                                 &
                     default        = .false.,                                          &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR50'
        
        
         Me%ObjHorizontalGridOut2 = 0
        
         
        !Kill enter data
        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridKeywords - ModuleMeshGlue - ERR60'
    
    end subroutine ReadGridOut
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

    subroutine ModifyMeshGlue(ObjMeshGlueID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjMeshGlueID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjMeshGlueID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call WriteNewGrid

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyMeshGlue

    !--------------------------------------------------------------------------    
            
    subroutine WriteNewGrid
    
        !Local-----------------------------------------------------------------
        real,   pointer, dimension(:,:)         :: XX_IE, YY_IE
        real,   pointer, dimension(:)           :: Dummy
        real                                    :: Latitude, Longitude
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
  
        call GetHorizontalGridSize(HorizontalGridID = Me%ObjHorizontalGridOut1,         &
                                   Size             = Me%Size,                          &
                                   WorkSize         = Me%WorkSize,                      &
                                   STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteNewGrid - ModuleMeshGlue - ERR10'
        
        call GetLatitudeLongitude(HorizontalGridID = Me%ObjHorizontalGridOut1,          &
                                  Latitude         = Latitude,                          &
                                  Longitude        = Longitude,                         &
                                  STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteNewGrid - ModuleMeshGlue - ERR20'
        
        allocate(Me%XX_IE(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%YY_IE(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))        
        
        call GetGridLatitudeLongitude(HorizontalGridID  = Me%ObjHorizontalGridOut1,     &
                                      GridLatitudeConn  = YY_IE,                        &
                                      GridLongitudeConn = XX_IE,                        &
                                      STAT              = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteNewGrid - ModuleMeshGlue - ERR30'
        
        Me%XX_IE(:,:) = XX_IE(:,:)
        Me%YY_IE(:,:) = YY_IE(:,:)
        
        if (Me%ResetGrid) then
            Me%XX_IE(:,:) = FillValueReal
            Me%YY_IE(:,:) = FillValueReal
        endif
        
        call UpdateGrid

        call ConstructHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGridOut2,       &
                                     LatitudeConn     = Me%YY_IE,                       &
                                     LongitudeConn    = Me%XX_IE,                       &
                                     XX               = Dummy,                          &
                                     YY               = Dummy,                          &
                                     Latitude         = Latitude,                       &
                                     Longitude        = Longitude,                      &
                                     ILB              = Me%WorkSize%ILB,                &
                                     IUB              = Me%WorkSize%IUB,                &
                                     JLB              = Me%WorkSize%JLB,                &
                                     JUB              = Me%WorkSize%JUB,                &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteNewGrid - ModuleMeshGlue - ERR40'  
        
        call WriteGridData_ASCII (FileName          = Me%FileNameGrid_Out,              &
                                  COMENT1           = "Mesh Glue",                      &
                                  COMENT2           = "Line 2",                         &
                                  HorizontalGridID  = Me%ObjHorizontalGridOut2,         &
                                  FillValue         = FillValueReal,                    &
                                  Overwrite         = .true.,                           &
                                  DistortionYes     = .true.,                           &    
                                  STAT              = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteNewGrid - ModuleMeshGlue - ERR50'
        
    
    end subroutine WriteNewGrid
    !--------------------------------------------------------------------------      
    
    subroutine UpdateGrid()
    
        !Local-----------------------------------------------------------------
        real,   pointer, dimension(:,:)         :: XX_IE, YY_IE
        real,   pointer, dimension(:,:)         :: XX_Aux, YY_Aux, Aux2D
        real,   pointer, dimension(:)           :: Dummy
        real                                    :: Latitude, Longitude
        integer                                 :: STAT_CALL, ifl, io, jo, di, dj
        integer                                 :: ILB, IUB, JLB, JUB, i, j

        !Begin-----------------------------------------------------------------
        
        
        do ifl =1, Me%GridIn_Number
            
            
            call GetGridLatitudeLongitude(HorizontalGridID  = Me%GridIn(ifl)%ObjHorizontalGrid, &
                                          GridLatitudeConn  = YY_IE,                            &
                                          GridLongitudeConn = XX_IE,                            &
                                          STAT              = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteNewGrid - ModuleMeshGlue - ERR30'            
            
            
            io = Me%GridIn(ifl)%IO
            jo = Me%GridIn(ifl)%JO            
            
            if (Me%GridIn(ifl)%Flip_IJ) then

                ILB = Me%GridIn(ifl)%JLB
                IUB = Me%GridIn(ifl)%JUB
                JLB = Me%GridIn(ifl)%ILB
                JUB = Me%GridIn(ifl)%IUB
                
                allocate(XX_Aux(ILB:IUB+1,JLB:JUB+1))
                allocate(YY_Aux(ILB:IUB+1,JLB:JUB+1))                
                
                do i = ILB, IUB + 1
                do j = JLB, JUB + 1
                    XX_Aux(i, j) = XX_IE(j, i)
                    YY_Aux(i, j) = YY_IE(j, i)
                enddo
                enddo

            else
                
                ILB = Me%GridIn(ifl)%ILB
                IUB = Me%GridIn(ifl)%IUB
                JLB = Me%GridIn(ifl)%JLB
                JUB = Me%GridIn(ifl)%JUB
                
                allocate(XX_Aux(ILB:IUB+1,JLB:JUB+1))
                allocate(YY_Aux(ILB:IUB+1,JLB:JUB+1))                
                
                XX_Aux(ILB:IUB+1,JLB:JUB+1) = XX_IE(ILB:IUB+1,JLB:JUB+1)
                YY_Aux(ILB:IUB+1,JLB:JUB+1) = YY_IE(ILB:IUB+1,JLB:JUB+1)
                
            endif   
            
            if (Me%GridIn(ifl)%RowFlip) then

                allocate(Aux2D(ILB:IUB+1,JLB:JUB+1))
                
                Aux2D(:,:) = XX_Aux(:,:)

                do i = ILB, IUB + 1
                do j = JLB, JUB + 1
                    XX_Aux(i, j) = Aux2D(IUB + 2 - i, j)
                enddo
                enddo                
                
                Aux2D(:,:) = YY_Aux(:,:)

                do i = ILB, IUB + 1
                do j = JLB, JUB + 1
                    YY_Aux(i, j) = Aux2D(IUB + 2 - i, j)
                enddo
                enddo                                

            endif   

            if (Me%GridIn(ifl)%ColumnFlip) then

                allocate(Aux2D(ILB:IUB+1,JLB:JUB+1))
                
                Aux2D(:,:) = XX_Aux(:,:)

                do i = ILB, IUB + 1
                do j = JLB, JUB + 1
                    XX_Aux(i, j) = Aux2D(i, JUB + 2 - j)
                enddo
                enddo                
                
                Aux2D(:,:) = YY_Aux(:,:)

                do i = ILB, IUB + 1
                do j = JLB, JUB + 1
                    YY_Aux(i, j) = Aux2D(i, JUB + 2 - j)
                enddo
                enddo                 

            endif
            
            dj = JUB + 1 - JLB 
            di = IUB + 1 - ILB 
            
            Me%XX_IE(io:io+di,jo:jo+dj) = XX_Aux(ILB:IUB+1, JLB:JUB+1)
            
            Me%YY_IE(io:io+di,jo:jo+dj) = YY_Aux(ILB:IUB+1, JLB:JUB+1)
            
            deallocate(XX_Aux, YY_Aux)
            
            if (associated(Aux2D)) deallocate(Aux2D)
           
        enddo
        

        
    end subroutine UpdateGrid
    
    !--------------------------------------------------------------------------         

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillMeshGlue(ObjMeshGlueID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjMeshGlueID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjMeshGlueID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mMeshGlue_,  Me%InstanceID)

            if (nUsers == 0) then

                !Deallocates Instance
                call DeallocateInstance ()

                ObjMeshGlueID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillMeshGlue
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_MeshGlue), pointer          :: AuxObjMeshGlue
        type (T_MeshGlue), pointer          :: PreviousObjMeshGlue

        !Updates pointers
        if (Me%InstanceID == FirstObjMeshGlue%InstanceID) then
            FirstObjMeshGlue => FirstObjMeshGlue%Next
        else
            PreviousObjMeshGlue => FirstObjMeshGlue
            AuxObjMeshGlue      => FirstObjMeshGlue%Next
            do while (AuxObjMeshGlue%InstanceID /= Me%InstanceID)
                PreviousObjMeshGlue => AuxObjMeshGlue
                AuxObjMeshGlue      => AuxObjMeshGlue%Next
            enddo

            !Now update linked list
            PreviousObjMeshGlue%Next => AuxObjMeshGlue%Next

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

    subroutine Ready (ObjMeshGlue_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjMeshGlue_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjMeshGlue_ID > 0) then
            call LocateObjMeshGlue (ObjMeshGlue_ID)
            ready_ = VerifyReadLock (mMeshGlue_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjMeshGlue (ObjMeshGlueID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjMeshGlueID

        !Local-----------------------------------------------------------------

        Me => FirstObjMeshGlue
        do while (associated (Me))
            if (Me%InstanceID == ObjMeshGlueID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleMeshGlue - LocateObjMeshGlue - ERR01'

    end subroutine LocateObjMeshGlue

    !--------------------------------------------------------------------------

end module ModuleMeshGlue









