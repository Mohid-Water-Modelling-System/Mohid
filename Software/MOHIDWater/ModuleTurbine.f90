!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Shell
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to serve as shell to create new modules
!
!------------------------------------------------------------------------------


Module ModuleTurbine

    use ModuleGlobalData
    
    use ModuleEnterData,        only : ReadFileName, ConstructEnterData, GetData,       &
                                        ExtractBlockFromBuffer, ExtractBlockFromBlock,  &
                                        Block_Unlock, GetOutPutTime, RewindBuffer,      &
                                        GetKeywordFromLine, GetFullBufferLine,          &
                                        ReplaceFullBufferLine, GetBlockSize,            &
                                        KillEnterData
    
    use ModuleHorizontalGrid,   only : GetXYCellZ
    
                                        
    implicit none

    private !What does that mean? all declarations are private if not specified??

    !Subroutines---------------------------------------------------------------

!    public :: ForceTurbine    
    
    !Constructor
    public  :: ConstructTurbine2
    private ::      AllocateInstance
    private ::      AllocateTurbine
    private ::      AllocateDataTurbine
    !Selector
    public  :: GetShellPointer
    public  :: GetShellInteger
    public  :: UnGetShell
                     
    
    !Modifier
    public  :: ModifyShell
    !public  :: CalculateTurbineForces

    !Destructor
    public  :: KillTurbine                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjShell 
    
    !Parameter
    !Parameter
    character(LEN = StringLength), parameter    :: list_begin          = '<beginturbinelist>'
    character(LEN = StringLength), parameter    :: list_end            = '<endturbinelist>'
    character(LEN = StringLength), parameter    :: turbine_begin       = '<<beginturbine>>'
    character(LEN = StringLength), parameter    :: turbine_end         = '<<endturbine>>'
    
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetShell3D_I
    private :: UnGetShell3D_R8
    interface  UnGetShell
        module procedure UnGetShell3D_I
        module procedure UnGetShell3D_R8
    end interface  UnGetShell

    !Types---------------------------------------------------------------------
    private :: T_Turbine_Param
    type       T_Turbine_Param
       real(8)         :: Diameter
       real(8)         :: H
       real(8)         :: Cto
       real(8)         :: LowerVelocity
       real(8)         :: UpperVelocity
       real(8)         :: CD
       real(8)         :: Width
       real(8)         :: Lat_coord
       real(8)         :: Long_coord
       integer         :: TimeSeriesOut             !1:We have output data 0: we don't have outputdata
       integer, dimension(2)            :: Cell_IJ  !Vector 2D where we store the cell where the turbine is placed
       integer         :: ID
       type(T_Turbine_Param), pointer   ::Next
    end type T_Turbine_Param
                

    private :: T_Turbine
    type       T_Turbine
        integer                                     :: InstanceID
        type (T_Size3D)                             :: Size, WorkSize
        real(8), dimension(:, :, :),  pointer       :: Matrix
        
        integer                                     :: nTurbines = 0
        type(T_Turbine_Param), pointer              :: Turbine !Debo verificar si es puntero o no
        type(T_Turbine), pointer                      :: Next
        character(PathLength)                       :: DataFile
        integer                                     :: ObjEnterData     !ID del objeto del modulo de lectura de datos
        integer                                     :: HorizontalGridID !ID de la grid con la que estamos trabajando
    end type  T_Turbine

    
    


    !Global Module Variables
    
    type (T_Turbine), pointer                     :: FirstObjShell
    type (T_Turbine), pointer                     :: Me
    
    integer                                         :: mSHELL_ = 0 !just to compile Se puede borrar esto¿?

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    subroutine ConstructTurbine2(TurbineID, HorizontalGridID, STAT) !ModelName
        
        
        !Arguments---------------------------------------------------------------
        integer, intent (INOUT)         :: TurbineID
        integer, optional, intent(OUT)  :: STAT
        integer, intent(IN)             :: HorizontalGridID
        
        !Local-------------------------------------------------------------------
        
        type (T_Turbine_param), pointer :: NewTurbine
        integer                         :: ready_
        integer                         :: STAT_, STAT_CALL
        integer                         :: ClientNumber
        integer                         :: FirstLine, LastLine
        !integer                         :: flag En principio no hace falta
        logical                         :: BlockFound, BlockInBlockFound
              
        
        STAT_ = UNKNOWN_
        
    
        if (.not. ModuleIsRegistered(mTurbine_)) then 
            nullify (FirstObjShell)
            call RegisterModule (mTurbine_)
        endif
        
        call Ready(TurbineID, ready_)
        
if0 :   if (ready_.EQ. OFF_ERR_) then
            
            call AllocateInstance 
            
            call ReadFileName('TURBINE', Me%DataFile, STAT=STAT_CALL)
            !Ahora tenemos en la variable DataFile la ruta del fichero turbina introducida en el nomfich.dat
            if (STAT_CALL /= SUCCESS_) stop 'ConsrtuctTurbine - TurbineModule - ERR10'
            
            !Construct enter data
            call ConstructEnterData(Me%ObjEnterData, Me%DataFile, STAT=STAT_CALL)
     
            !Quizas aqui se podria poneer otro loop en caso de que haya mas de una lista de turbinas
            call ExtractBlockFromBuffer (Me%ObjEnterData, ClientNumber, list_begin, list_end, BlockFound, &
                                         FirstLine, LastLine, STAT=STAT_CALL)
            
            if (BlockFound) then
                
DOPROP:         do   
                    
                    call ExtractBlockFromBlock (Me%ObjEnterData, ClientNumber, turbine_begin, turbine_end, BlockInBlockFound, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR1140' !Desconozco si el codigo de error es correcto
                    
if1:                if (BlockInBlockFound) then
                        call AllocateDataTurbine(HorizontalGridID,STAT_CALL)    
                    else if1
                        exit DOPROP
                    endif if1
                enddo DOPROP
            endif
        
        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR1600'

        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR14'
        
        !Returns ID
        TurbineID   = Me%InstanceID
        STAT_       = SUCCESS_
        
        else if0
            stop 'ConstructTurbine - ModuleTurbine - ERR99'
        endif if0
        
        if (present(STAT)) STAT = STAT_
        
    end subroutine ConstructTurbine2

   
    subroutine AllocateDataTurbine (HorizontalGridID, STAT)   !Revisar si hace falta inclur el flag
    
        !Arguments---------------------------------------------------------------
        integer, optional, intent(OUT)  :: STAT
        integer, intent (IN)            :: HorizontalGridID  
        !Local-------------------------------------------------------------------
        type (T_Turbine_param), pointer :: NewTurbine
        integer                         :: STAT_CALL
        integer                         :: flag
        !Allocates new instance
        allocate (NewTurbine)
        nullify  (NewTurbine%Next)                        
                    
        call GetData(NewTurbine%Diameter,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &!Debo tratar el flag ¿?
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='DIAMETER',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR1060'
                       
        call GetData(NewTurbine%H,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='HEIGHT',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR1060'
                       
        call GetData(NewTurbine%Cto,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='CTO',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR1060'
                        
        call GetData(NewTurbine%LowerVelocity,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='LOWER_VEL',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR1060'
                       
        call GetData(NewTurbine%UpperVelocity,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='UPPER_VEL',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR1060'
                       
        call GetData(NewTurbine%CD,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='CD',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR1060'
                        
        call GetData(NewTurbine%Lat_coord,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='POS_LAT',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR1060'
                       
        call GetData(NewTurbine%Long_coord,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='POS_LONG',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR1060'
        call GetData(NewTurbine%TimeSeriesOut,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='TIMESERIES',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR1060'
        
        call GetData(NewTurbine%Width,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='WIDTH_STRUCT',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR1060'
        
        
        !Saves the position of the turbine in grid coordinates
        call GetXYCellZ(HorizontalGridID, XPoint = NewTurbine%Lat_coord,Ypoint = NewTurbine%Long_coord,   &
                        I = NewTurbine%Cell_IJ(1), J=NewTurbine%Cell_IJ(2))
        
        call AllocateTurbine(NewTurbine) 
        
        
        if (present(STAT)) STAT = STAT_CALL
        nullify(NewTurbine)
                    
   end subroutine AllocateDataTurbine
                
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance !Guarda espacio de memoria para 
        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Turbine), pointer                         :: NewObjShell
        type (T_turbine), pointer                         :: PreviousObjShell


        !Allocates new instance
        allocate (NewObjShell)              !Se guarda un espacio de memoria para el nuevo objeto NewObjShell que es de tipo TurbineType
        nullify  (NewObjShell%Next)         !Se asegura que el siguente objeto tipo TurbineType no apunta a nada.        

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjShell)) then       !FirstObjShell es una variable global, en esta condicion se le asocia una dirección a este puntero en caso de no tenerla
            FirstObjShell         => NewObjShell
            Me                    => NewObjShell        !Se apuntan ambas variables/punteros al mismo puntero (este es una puntero local de la subrutina)
        else
            PreviousObjShell      => FirstObjShell      !En el caso de que la variable local FirstObjShell ya tenga una dirección asociada, la variable local PreviosuObjShell
            Me                    => FirstObjShell%Next !apunta al primer objeto mientras que la variable Me apunta al siguiente
            do while (associated(Me))
                PreviousObjShell  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjShell
            PreviousObjShell%Next => NewObjShell
        endif

        Me%InstanceID = RegisterNewInstance (mTURBINE_)


    end subroutine AllocateInstance
    
    
    subroutine AllocateTurbine(NewTurbine) !Guarda espacio de memoria para 
        !Arguments-------------------------------------------------------------
        type (T_Turbine_param), pointer :: NewTurbine
    
        !Local-----------------------------------------------------------------
        type (T_Turbine_param), pointer :: CurrentTurbine => null()
        type (T_Turbine_param), pointer :: PreviousTurbine => null()


        !Insert New Instance into list and makes Current point to it
        if (.not. associated(Me%Turbine)) then       !FirstObjShell es una variable global, en esta condicion se le asocia una dirección a este puntero en caso de no tenerla
            !FirstObjShell%Turbine         => NewObjShell
            Me%Turbine                    => NewTurbine        !Se apuntan ambas variables/punteros al mismo puntero (este es una puntero local de la subrutina)
        else
            PreviousTurbine      => Me%Turbine      !En el caso de que la variable local FirstObjShell ya tenga una dirección asociada, la variable local PreviosuObjShell
            CurrentTurbine       => PreviousTurbine%Next !apunta al primer objeto mientras que la variable Me apunta al siguiente
            do while (associated(CurrentTurbine))
                PreviousTurbine  => CurrentTurbine
                CurrentTurbine   => PreviousTurbine%Next
            enddo
            PreviousTurbine%Next => NewTurbine
        endif
        
        Me%nTurbines    = Me%nTurbines + 1  !Se utiliza el numero de turbinas para darles un identificador sequencial 
        NewTurbine%ID   =   Me%nTurbines
        !Me%InstanceID = RegisterNewInstance (mTURBINE_)
        

    end subroutine AllocateTurbine


    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    
    !--------------------------------------------------------------------------
    subroutine GetShellPointer (ObjShellID, Matrix, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjShellID
        real(8), dimension(:, :, :),  pointer           :: Matrix
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjShellID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mSHELL_, Me%InstanceID)

            Matrix => Me%Matrix

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetShellPointer
    
    !--------------------------------------------------------------------------
    
    subroutine GetShellInteger (ObjShellID, Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjShellID
        real                                            :: Int
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjShellID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Int = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetShellInteger

    !--------------------------------------------------------------------------

    subroutine UnGetShell3D_I(ObjShellID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjShellID
        integer, dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjShellID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSHELL_, Me%InstanceID, "UnGetShell3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetShell3D_I

    !--------------------------------------------------------------------------

    subroutine UnGetShell3D_R8(ObjShellID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjShellID
        real(8), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjShellID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSHELL_, Me%InstanceID,  "UnGetShell3D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetShell3D_R8

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !subroutine GetTurbineForces2(Turbines, Velocity_U, Velocity_V, STAT) !Trabajar con matrices
    !!Arguments-------------------------------------------------------------------
    !Type(T_Turbine) pointer             :: Turbines
    !real, dimension(:,:) pointer        :: Velocity_U, Velocity_V, STAT)
    !integer                             :: STAT
    !!Local-----------------------------------------------------------------------
    !real(8)                             :: FU, FV     
    !real(8)                             :: U
    !
    !
    !
    !end subroutine GetTurbineForces2
    !
    !subroutine CalculateTurbineForces(Turbine, Velocity_U, Velocity_V, STAT) !Buscar en moduleHoritzontalData funcion devuelve gridCell
    !!Faltan comprobaciones modulo
    !!Arguments------------------------------------------------------------------
    !Type(T_Turbine_param), pointer       :: Turbine
    !real, pointer                        :: Velocity_U, Velocity_V
    !integer                             :: STAT
    !
    !!Local----------------------------------------------------------------------
    !real(8)                             :: FU, FV     
    !real(8)                             :: VelocityModul
    !
    !U=(SQRT(Velocity_U%OLD*Velocity_V%OLD) + SQRT(Velocity_U%NEW*Velocity_V%NEW))/2.0
    !if (U .GE. Turbine%U_inf .AND. U.LE.Turbine%U_sup) then 
    !    
    !    FU = (0.5*DENSITY*PI*0.25*(Turbine%D_Turbine**2)*Turbine_Cto + 0.5*DENSITY*Turbine%H*Turbine%W*Turbine%CD)*(Velocity_U%OLD*Velocity_U%NEW)
    !    FV = (0.5*DENSITY*PI*0.25*(Turbine%D_Turbine**2)*Turbine_Cto + 0.5*DENSITY*Turbine%H*Turbine%W*Turbine%CD)*(Velocity_V%OLD*Velocity_V%NEW)
    !else
    !    FU = (0.5*DENSITY*Turbine%H*Turbine%W*Turbine%CD)*(Velocity_U%OLD*Velocity_U%NEW)
    !    FV = (0.5*DENSITY*Turbine%H*Turbine%W*Turbine%CD)*(Velocity_V%OLD*Velocity_V%NEW)
    !end if    
    !
    ! 
    !
    !
    !
    !end subroutine CalculateTurbineForces
    
  
    !subroutine ActualiseHydroninamicCoeficient(STAT)
    !
    !end subroutine ActualiseHydronamicCoeficient(STAT)
    
    subroutine ModifyShell(ObjShellID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjShellID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjShellID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then




            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyShell


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillTurbine(TurbineID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: TurbineID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbineID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mTURBINE_,  Me%InstanceID)

            if (nUsers == 0) then

                !Deallocates Instance
                call DeallocateInstance ()

                TurbineID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillTurbine
        

    !------------------------------------------------------------------------
    
    subroutine DeallocateTurbine ()
    
    end subroutine DeallocateTurbine
    
    subroutine DeallocateInstance ()
    
        !Arguments-------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
        type (T_Turbine), pointer          :: AuxObjShell
        type (T_Turbine), pointer          :: PreviousObjShell
    
        !Updates pointers
        if (Me%InstanceID == FirstObjShell%InstanceID) then
            FirstObjShell => FirstObjShell%Next
        else
            PreviousObjShell => FirstObjShell
            AuxObjShell      => FirstObjShell%Next
            do while (AuxObjShell%InstanceID /= Me%InstanceID)
                PreviousObjShell => AuxObjShell
                AuxObjShell      => AuxObjShell%Next
            enddo
    
            !Now update linked list
            PreviousObjShell%Next => AuxObjShell%Next
    
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

    subroutine Ready (ObjTurbine_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTurbine_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjTurbine_ID > 0) then
            call LocateObjShell (ObjTurbine_ID)
            ready_ = VerifyReadLock (mTURBINE_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjShell (ObjShellID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjShellID

        !Local-----------------------------------------------------------------

        Me => FirstObjShell
        do while (associated (Me))
            if (Me%InstanceID == ObjShellID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleShell - LocateObjShell - ERR01'

    end subroutine LocateObjShell

    !--------------------------------------------------------------------------

end module ModuleTurbine