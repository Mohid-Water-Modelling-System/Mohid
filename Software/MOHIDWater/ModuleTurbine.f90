!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Turbine
! PROJECT       : Mohid Water
! MODULE        : ModuleTurbine
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Sept 2017
! REVISION      : Oscar Balsells - v4.0
! DESCRIPTION   : Module for calculating the force exerted by turbines, the power output and the energy extraction.
!
!------------------------------------------------------------------------------


Module ModuleTurbine

    use ModuleGlobalData
    
    use ModuleEnterData,        only : ReadFileName, ConstructEnterData, GetData,       &
                                       ExtractBlockFromBuffer, ExtractBlockFromBlock,  &
                                       Block_Unlock, KillEnterData
    
    use ModuleHorizontalGrid,    only : GetXYCellZ

    use ModuleGeometry,          only : GetGeometrySize, GetGeometryKFloor,               &
                                      GetGeometryDistances, UnGetGeometry
    
    use ModuleTimeSerie,         only : StartTimeSerieTurbine, KillTimeSerie, WriteTimeSerie

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    
    !Constructor
    public  :: ConstructTurbine
    private ::      AllocateInstance
    private ::      AllocateDataTurbine
    private ::          AllocateTurbine
    private ::      AllocateTimeSerieTurbine
    private ::      ConstructTimeSerieTurbine
    
    !Selector
    public  :: GetTurbineAcceleration
    public  :: UnGetTurbineAcceleration

                     
    
    !Modifier
    public  :: ModifyTurbine
    private ::      TurbineVerticalDiscretisation
    public  :: ComputeTurbineEnergy
    public  :: OutPut_turbine

    !Destructor
    public  :: KillTurbine
    private ::      Deassociate_External_Modules
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjShell 
    
    
    !Parameter
    character(LEN = StringLength), parameter    :: list_begin          = '<beginturbinelist>'
    character(LEN = StringLength), parameter    :: list_end            = '<endturbinelist>'
    character(LEN = StringLength), parameter    :: turbine_begin       = '<<beginturbine>>'
    character(LEN = StringLength), parameter    :: turbine_end         = '<<endturbine>>'
    


    !Types---------------------------------------------------------------------
    public :: T_Turbine_Param
    type       T_Turbine_Param
       real(8)              :: Diameter
       real(8)              :: H
       real(8)              :: CP,CT
       real(8)              :: LowerVelocity
       real(8)              :: UpperVelocity
       !real(8)              :: CD !Drag coef for the structure
       !real(8)              :: Width !Width of the structure
       real(8)              :: Lat_coord !Or y coord
       real(8)              :: Long_coord!Or x coord
       integer              :: TimeSerieOut     !1:We have output data 0: we don't have outputdata
       integer              :: I,J  !Vector 2D where we store the cell where the turbine is placed
       real, dimension(:), pointer   :: T_Area
       !real                 :: ForceU,ForceV
       integer              :: ID
       type(T_Turbine_Param), pointer   ::Next
    end type T_Turbine_Param
                

    private :: T_Turbine
    type       T_Turbine
        integer                                     :: InstanceID
        type (T_Size3D)                             :: Size, WorkSize
        real, dimension(:, :, :),  pointer          :: TurbineAcceleration
        !real, dimension(:, :, :), pointer           :: Coef
        real, dimension(:), pointer                 ::  Power, Energy, TurbineVelocity 
        real                                        :: TotalPower
        integer, dimension(:), pointer                 :: TurbineTimeSerieList
        
        integer                                     :: nTurbines = 0, nTurbinesTS = 0
        type(T_Turbine_Param), pointer              :: Turbine 
        type(T_Turbine), pointer                      :: Next
        character(PathLength)                       :: DataFile
        integer                                     :: ObjEnterData = 0    !ID del objeto del modulo de lectura de datos
        integer                                     :: ObjHorizontalGrid = 0!ID de la grid con la que estamos trabajando
        integer                                     :: ObjGeometry  = 0
        integer                                     :: ObjTimeSerie = 0
        integer                                     :: ObjTime      = 0
        logical                                     :: TimeSerieON = .false.
        real, dimension(:,:,:), pointer               :: DUZ      => null()
        integer, dimension(:,:), pointer            :: KFloor_Z   => null()
    end type  T_Turbine

    
    


    !Global Module Variables
    
    type (T_Turbine), pointer                     :: FirstObjShell
    type (T_Turbine), pointer                     :: Me
    
    integer                                         :: mSHELL_ = 0 !

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    
    subroutine ConstructTurbine(TurbineID, HorizontalGridID, GeometryID, TimeID, OutPut, STAT) !ModelName
            
        !Arguments---------------------------------------------------------------
        integer, intent (INOUT)         :: TurbineID
        integer, optional, intent(OUT)  :: STAT
        integer, intent(IN)             :: HorizontalGridID
        integer, intent(IN)             :: GeometryID
        integer, intent(IN)             :: TimeID
        logical, intent(OUT)            :: OutPut
        !Local-------------------------------------------------------------------
        
        integer                         :: ready_, iflag
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
            
            !Obtaining geometry to create the matrix for the forces
            call GetGeometrySize(GeometryID, Me%Size, Me%WorkSize, STAT)
            !Seteamos la matriz con 0
            allocate (Me%TurbineAcceleration (Me%Size%ILB: Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
            
            !Associates External Instances
            Me%ObjHorizontalGrid     = AssociateInstance (mHORIZONTALGRID_,  HorizontalGridID)
            Me%ObjGeometry           = AssociateInstance (mGEOMETRY_,        GeometryID)
            Me%ObjTime               = AssociateInstance (mTIME_,            TimeID)
            
            !Reading input data file
            call ReadFileName('TURBINE', Me%DataFile, STAT=STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConsrtuctTurbine - TurbineModule - ERR01'
            
            !Construct enter data
            call ConstructEnterData(Me%ObjEnterData, Me%DataFile, STAT=STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConsrtuctTurbine - TurbineModule - ERR02'
            
            call GetData(Me%TimeSerieON,                                         &
                     Me%ObjEnterData, iflag,                                            & 
                     Keyword    = 'TIMESERIE',                            &
                     Default    = .false.,                                              &
                     SearchType = FromFile,                                             &
                     ClientModule ='ModuleHydrodynamic',                                &
                     STAT       = STAT_CALL)     
            if (STAT_CALL /= SUCCESS_) stop 'ConsrtuctTurbine - TurbineModule - ERR03'
            
            OutPut = Me%TimeSerieON
            
            call ExtractBlockFromBuffer (Me%ObjEnterData, ClientNumber, list_begin, list_end, BlockFound, &
                                         FirstLine, LastLine, STAT=STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR04'
            
            if (BlockFound) then
                
DOPROP:         do   
                    
                    call ExtractBlockFromBlock (Me%ObjEnterData, ClientNumber,          &
                        turbine_begin, turbine_end, BlockInBlockFound, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR05' 
                    
if1:                if (BlockInBlockFound) then
                        call AllocateDataTurbine    
                    else if1
                        exit DOPROP
                    endif if1
                    
                enddo DOPROP
            endif
        
            if (Me%TimeSerieON) then
                call AllocateTimeSerieTurbine
                call ConstructTimeSerieTurbine
            endif
            
            allocate(Me%Power(Me%nTurbines + 1))    
            allocate(Me%Energy(Me%nTurbines + 1))
            allocate (Me%TurbineVelocity(Me%nTurbines + 1))
            Me%Power(:) = 0.0 
            Me%Energy(:) = 0.0 
            Me%TurbineVelocity(:) = 0.0

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR06'

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTurbine - ModuleTurbine - ERR07'
        
            !Returns ID
            TurbineID   = Me%InstanceID
            STAT_       = SUCCESS_
        
        else if0
            STAT_ = UNKNOWN_
            !stop 'ConstructTurbine - ModuleTurbine - ERR99'
        endif if0
        
            if (present(STAT)) STAT = STAT_
        
    end subroutine ConstructTurbine
    
    subroutine AllocateTimeSerieTurbine
        !Local-----------------------------------------------------------------
        type (T_Turbine_param), pointer     :: Turbine, PreviousTurbine
        integer                             :: i = 1
        !Begin-----------------------------------------------------------------
        Turbine         => Me%Turbine
        PreviousTurbine     => null()
        
        allocate(Me%TurbineTimeSerieList(Me%nTurbinesTS + 1)) !+1 is for the general turbine output
        
        do while (associated(Turbine))
            if (Turbine%TimeSerieOut == 1) then 
                Me%TurbineTimeSerieList(i) = Turbine%ID
                i = i + 1
            endif 
           PreviousTurbine  => Turbine
           Turbine          => PreviousTurbine%Next
        enddo
        Me%TurbineTimeSerieList(i) = Me%nTurbines + 1
        
        nullify(Turbine)            
        nullify(PreviousTurbine)
    
    end subroutine AllocateTimeSerieTurbine
    
    subroutine ConstructTimeSerieTurbine
    
        !Local-----------------------------------------------------------------
         character(len=StringLength), dimension(:), pointer  :: PropertyList
         integer                                             :: STAT_CALL
        !Begin-----------------------------------------------------------------
         allocate(PropertyList(3)) !3 Prop, Power [W], Energy [KJ], Velocity [m/s]
         PropertyList(1) = 'Power [kW]'
         PropertyList(2) = 'Energy [MJ]'
         PropertyList(3) = 'Velocity [m/s]'
         
         call StartTimeSerieTurbine(Me%ObjTimeSerie, Me%ObjTime, Me%ObjEnterData, Me%TurbineTimeSerieList, PropertyList, &
                                    "srh", STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_)                     &
            call SetError (FATAL_, OUT_OF_MEM_, "ConstructTimeSeriesTurbine - ModuleTurbine - ERR08")                      
    
    end subroutine ConstructTimeSerieTurbine

   
    subroutine AllocateDataTurbine   
    
        !Arguments---------------------------------------------------------------
        !Local-------------------------------------------------------------------
        type (T_Turbine_param), pointer :: NewTurbine
        integer                         :: STAT_CALL
        integer                         :: flag
        !Allocates new instance
        allocate (NewTurbine)
        nullify  (NewTurbine%Next)                        
                    
        call GetData(NewTurbine%Diameter,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='DIAMETER',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'AllocateDataTurbine - ModuleTurbine - ERR09'
                       
        call GetData(NewTurbine%H,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='HEIGHT',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'AllocateDataTurbine - ModuleTurbine - ERR10'
                       
        call GetData(NewTurbine%CP,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='CP',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'AllocateDataTurbine - ModuleTurbine - ERR11'
        
        call GetData(NewTurbine%CT,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='CT',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'AllocateDataTurbine - ModuleTurbine - ERR12'        
        
        call GetData(NewTurbine%LowerVelocity,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='LOWER_VEL',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'AllocateDataTurbine - ModuleTurbine - ERR13'
                       
        call GetData(NewTurbine%UpperVelocity,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='UPPER_VEL',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'AllocateDataTurbine - ModuleTurbine - ERR14'
                       
        !call GetData(NewTurbine%CD,                           &
        !            Me%ObjEnterData,                                        &
        !            flag,                                                   &
        !            SearchType   = FromBlockInBlock,                               &
        !            keyword      ='CD',                              &
        !            ClientModule ='Moduleturbine',                          &
        !            STAT         = STAT_CALL)             
        !if (STAT_CALL /= SUCCESS_) stop 'AllocateDataTurbine - ModuleTurbine - ERR15'
                        
        call GetData(NewTurbine%Lat_coord,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='POS_LAT',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'AllocateDataTurbine - ModuleTurbine - ERR16'
                       
        call GetData(NewTurbine%Long_coord,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='POS_LONG',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'AllocateDataTurbine - ModuleTurbine - ERR17'
        call GetData(NewTurbine%TimeSerieOut,                           &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromBlockInBlock,                               &
                    keyword      ='TIMESERIE',                              &
                    ClientModule ='Moduleturbine',                          &
                    STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'AllocateDataTurbine - ModuleTurbine - ERR18'
        
        if (NewTurbine%TimeSerieOut /=0) then
            Me%nTurbinesTS = Me%nTurbinesTS + 1
        endif
        
        !call GetData(NewTurbine%Width,                           &
        !            Me%ObjEnterData,                                        &
        !            flag,                                                   &
        !            SearchType   = FromBlockInBlock,                               &
        !            keyword      ='WIDTH_STRUCT',                              &
        !            ClientModule ='Moduleturbine',                          &
        !            STAT         = STAT_CALL)             
        !if (STAT_CALL /= SUCCESS_) stop 'AllocateDataTurbine - ModuleTurbine - ERR19'
        
        
        !Saves the position of the turbine in grid coordinates
        call GetXYCellZ(Me%ObjHorizontalGrid, XPoint = NewTurbine%Long_coord,Ypoint = NewTurbine%Lat_coord,   &
                        I = NewTurbine%I, J=NewTurbine%J, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AllocateDataTurbine - ModuleTurbine - ERR20'

        call AllocateTurbine(NewTurbine) 
        
        
        !if (present(STAT)) STAT = STAT_CALL
        nullify(NewTurbine)
                    
   end subroutine AllocateDataTurbine
                
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance !Guarda espacio de memoria para 
        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Turbine), pointer                         :: NewObjShell
        type (T_turbine), pointer                         :: PreviousObjShell


        !Allocates new instance
        allocate (NewObjShell)              
        nullify  (NewObjShell%Next)         
        
        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjShell)) then      
            FirstObjShell         => NewObjShell
            Me                    => NewObjShell         
        else
            PreviousObjShell      => FirstObjShell      
            Me                    => FirstObjShell%Next 
            do while (associated(Me))
                PreviousObjShell  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjShell
            PreviousObjShell%Next => NewObjShell
        endif

        Me%InstanceID = RegisterNewInstance (mTURBINE_)


    end subroutine AllocateInstance
    
    
    subroutine AllocateTurbine(NewTurbine) 
        !Arguments-------------------------------------------------------------
        type (T_Turbine_param), pointer :: NewTurbine
    
        !Local-----------------------------------------------------------------
        type (T_Turbine_param), pointer :: CurrentTurbine => null()
        type (T_Turbine_param), pointer :: PreviousTurbine => null()


        !Insert New Instance into list and makes Current point to it
        if (.not. associated(Me%Turbine)) then       
            Me%Turbine                    => NewTurbine        
        else
            PreviousTurbine      => Me%Turbine      
            CurrentTurbine       => PreviousTurbine%Next 
            do while (associated(CurrentTurbine))
                PreviousTurbine  => CurrentTurbine
                CurrentTurbine   => PreviousTurbine%Next
            enddo
            PreviousTurbine%Next => NewTurbine
        endif
        
        Me%nTurbines    = Me%nTurbines + 1  
        NewTurbine%ID   =   Me%nTurbines

    end subroutine AllocateTurbine


    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetTurbineAcceleration (TurbineID, TurbineAcceleration, STAT)
        !Arguments-------------------------------------------------------------
        integer                         :: TurbineID
        real, dimension (:,:,:), pointer :: TurbineAcceleration
        integer, optional               :: STAT
        !Local-----------------------------------------------------------------
        integer                         :: STAT_, ready_
    
        !----------------------------------------------------------------------
    
        STAT_ = UNKNOWN_
    
        call Ready(TurbineID, ready_)
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
        
            call Read_Lock(mTURBINE_, Me%InstanceID)
    
            TurbineAcceleration => Me%TurbineAcceleration
    
            STAT_ = SUCCESS_
    
        else 
            STAT_ = ready_
        end if
    
        if (present(STAT)) STAT = STAT_
    
    end subroutine GetTurbineAcceleration
    
    
    subroutine UnGetTurbineAcceleration(TurbineID, TurbineAcceleration, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: TurbineID
        real, dimension(:, :, :), pointer            :: TurbineAcceleration
        integer, intent(OUT), optional                  :: STAT
    
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
    
        !----------------------------------------------------------------------
    
        STAT_ = UNKNOWN_
    
        call Ready(TurbineID, ready_)
    
        if (ready_ .EQ. READ_LOCK_ERR_) then
    
            nullify(TurbineAcceleration)
            call Read_Unlock(mTURBINE_, Me%InstanceID, "UnGetTurbineAcceleration")
    
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if
    
        if (present(STAT)) STAT = STAT_  
    end subroutine UnGetTurbineAcceleration
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    subroutine ModifyTurbine(TurbineID,VelocityU, VelocityV, VelocityUV, VelocityVU, VolumeUV,  &
                              KFloor_UV, DXX_YY, DUX_VY, di, dj, Density, DT, STAT)
    
        !Arguments------------------------------------------------------------------
        integer                             :: TurbineID
        real, dimension(:,:,:), pointer     :: VelocityUV, VelocityVU, VolumeUV, Density, VelocityU, VelocityV
        integer, dimension(:,:), pointer    :: KFloor_UV
        real, dimension(:,:), pointer       :: DXX_YY, DUX_VY
        integer, optional, intent (OUT)     :: STAT
        integer                             :: di, dj
        real                                :: DT, CT, CP, DensityAv
        !Local-----------------------------------------------------------------
        logical                             :: START = .false.
        real, dimension(:), pointer         :: T_Area
        integer                             :: STAT_CALL   
        real                                :: VelocityModul, VelocityModulTurbine
        type (T_Turbine_param), pointer     :: PreviousTurbine, Turbine
        integer                             :: I,J,K
        integer                             :: KUB
        integer                             :: KBottom
        !integer                             :: iSouth, I_North, J_East, jWest
        
        real                            ::  TotalAreaTurbine, aux
        integer :: ready_, STAT_
        !Inicialización de variables
         STAT_ = UNKNOWN_
        
        call Ready(TurbineID, ready_)
        
cd1 :   if (ready_ .EQ. IDLE_ERR_) then
        
            Turbine                             => Me%Turbine
            PreviousTurbine                     => null()
            T_Area                              => Turbine%T_Area
            KUB = Me%WorkSize%KUB

            !----------Calculate Area per cell in k domain-----------------------
            !Si se escoge esta rutina se deben de sacar estas dos funciones
                                                                        !,así solo la llamamos una ve 
            call GetGeometryKFloor( GeometryID = Me%ObjGeometry, Z = Me%KFloor_Z, STAT = STAT_CALL)
            call GetGeometryDistances(GeometryID = Me%ObjGeometry, DUZ = Me%DUZ, STAT = STAT_CALL)
        
        
            !----------Calculate the velocity modul-----------------------
            Me%TurbineAcceleration(:,:,:) = 0.0 !Reset to 0 because for every direction needs to be calculated again
           
            Me%Power(:)                   = 0.0 !Reset to 0 because the total power (Me%Power(Me%nTurbines+1)) 
                                                !is calculated again in each direction, the result for dj or di
                                                !is almost the same because the modulus of the velocity is calculated
                                                !in the center of the cell. There are other ways to compute the power
           
            Me%TurbineVelocity(:)        = 0.0  !The same reason as above, because the average Velocity modulus 
                                                !is calculated on each iteration, and for the 
                                                !total averga velocity (Me%TurbineVelocity(Me%nTurbines+1)) 
                                                !of all the turbines we need to reset this value
            
            
            do while (associated(Turbine))
                !Initialisation variables
                
                aux = 0.0
                !aux2= 0.0
                TotalAreaTurbine = 0.0
                DensityAv = 0.0

               I        = Turbine%I 
               J        = Turbine%J
               Kbottom  = KFloor_UV(I,J)
           
               allocate (T_Area(0:KUB))
               T_Area(:) =0. 
               !iSouth     = I - di
               !jWest      = J - dj 

               ! This values (i_North and j_east) can only be used to compute the velocity modulus in a face
               !I_North     = I + dj
               !J_East      = J + di
           
               !Vertical Discretisation 3D
               call TurbineVerticalDiscretisation(Turbine, T_Area)           
               
               do K=Kbottom, KUB
                   
                    !This codealculates the velocity modulus in the face, not in the center
                    !VelocityModul2 = Face_Velocity_Modulus(                                &
                    !                      VelocityVU(I_North, jWest, K),     &
                    !                      VelocityVU(I_North, J_East, K),    &
                    !                      VelocityVU(iSouth, jWest, K),      &
                    !                      VelocityVU(iSouth, J_East, K),     &
                    !                      DXX_YY(I_North, jWest),               &
                    !                      DXX_YY(I_North, J_East),              &
                    !                      DXX_YY(iSouth, jWest),                &
                    !                      DXX_YY(iSouth, J_East),               &
                    !                      VelocityUV(I,J,K))             
               
                   
                   VelocityModul = sqrt(((VelocityU(I,J,K)+VelocityU(I,J+1,K))/2.0)**2 +            &
                                    ((VelocityV(I,J,K)+VelocityV(I+1,J,K))/2.0)**2)
                   
                   aux = aux + VelocityModul*T_Area(K)
                   TotalAreaTurbine = TotalAreaTurbine + T_Area(K)
    
                
               end do
               !Calculates the average velocity on the turbine surface as restriction for the parametrisation of CP and CT.
               VelocityModulTurbine             = aux/TotalAreaTurbine
               Me%TurbineVelocity(Turbine%ID)   = VelocityModulTurbine
                  
               !Parametrisation of CP and CT coefficients
               if (VelocityModulTurbine .GE. (0.75*Turbine%LowerVelocity) .and. START &
                   .and. VelocityModulTurbine .LE. Turbine%UpperVelocity) then
                   CT = Turbine%CT
                   CP = Turbine%CP
               
               else if (VelocityModulTurbine .GE. Turbine%LowerVelocity) then
                   START = .true.
                   if (VelocityModulTurbine .LE. Turbine%UpperVelocity) then
                        CT = Turbine%CT
                        CP = Turbine%CP
                   else
                        CT = ((Turbine%CT)*Turbine%UpperVelocity**3)/(VelocityModulTurbine**3) 
                        CP = ((Turbine%CP)*Turbine%UpperVelocity**3)/(VelocityModulTurbine**3)
                   endif
               else
                   CT = 0.0
                   CP = 0.0
                   START = .false.
               end if
               
               !Compute the forces of the turbine
               do K=Kbottom, KUB
                   Me%TurbineAcceleration(I,J,K) = Me%TurbineAcceleration(I,J,K) -      &
                    (0.5*CT*T_Area(K)*(VelocityModul*VelocityUV(I,J,K))/VolumeUV(I,J,K))
                   DensityAv = DensityAv + Density(I,J,K)*T_Area(K)/TotalAreaTurbine
               enddo
                              
               Me%Power(Turbine%ID) = (0.5*CP*TotalAreaTurbine*DensityAv*VelocityModulTurbine**3)/1000
               
               Me%Power(Me%nTurbines+1) = Me%Power(Me%nTurbines+1) + Me%Power(Turbine%ID)
               Me%TurbineVelocity(Me%nTurbines+1) = Me%TurbineVelocity(Me%nTurbines+1) + Me%TurbineVelocity(Turbine%ID)/Me%nTurbines
                              PreviousTurbine  => Turbine
               Turbine          => PreviousTurbine%Next
    
            end do
        
            call UnGetGeometry(Me%ObjGeometry, Me%KFloor_Z)
            call UnGetGeometry(Me%ObjGeometry, Me%DUZ)

            nullify(Turbine)
            nullify(PreviousTurbine)
           
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1
        if (present(STAT)) STAT = STAT_
    
    end subroutine ModifyTurbine
        
    subroutine ComputeTurbineEnergy(TurbineID,DT, STAT)
        !Arguments-------------------------------------------------------------
        integer                         :: TurbineID
        integer, optional, intent (OUT) :: STAT
        real                            :: DT
        !Local-----------------------------------------------------------------
        integer                         :: STAT_, ready_    
        type (T_Turbine_param), pointer :: PreviousTurbine, Turbine
        
        !----------------------------------------------------------------------
    
        STAT_ = UNKNOWN_
    
        call Ready(TurbineID, ready_)
        
        if (ready_ .EQ. IDLE_ERR_) then
            Turbine             => Me%Turbine
            PreviousTurbine     => null()               
            
            do while (associated(Turbine))

                Me%Energy(Turbine%ID) = Me%Energy(Turbine%ID) + Me%Power(Turbine%ID)*DT/1000 ![MJ]
                Previousturbine     => Turbine
                Turbine             => PreviousTurbine%Next
            enddo
            Me%Energy(Me%nTurbines+1) = Me%Energy(Me%nTurbines+1) + Me%Power(Me%nTurbines+1)*DT/1000
            
        else 
            STAT_ = ready_
        end if
    
        if (present(STAT)) STAT = STAT_               
    
    end subroutine ComputeTurbineEnergy

    
    subroutine TurbineVerticalDiscretisation (Turbine,T_Area)
        !Arguments------------------------------------------------------------------
        real, dimension(:), pointer                 :: T_Area
        type(T_Turbine_param)                       :: Turbine
        !integer, optional                           :: STAT
        !Local----------------------------------------------------------------------
        real, parameter                             :: pi = 3.14159
        real                                        :: Theta, Beta, H, dist, Radius, PreviousArea, aux

        integer                                     :: q, I, J, K, Kbottom, KUB
        real, dimension(:,:,:), pointer               :: DUZ
        !ShortenNames---------------------------------------------------------------
        I       = Turbine%I
        J       = Turbine%J
        Radius  = (Turbine%Diameter)/2.0
        H       = Turbine%H
        KUB     = Me%WorkSize%KUB
        Kbottom = Me%KFloor_Z(I,J)
        DUZ     => Me%DUZ
        
        !Begin----------------------------------------------------------------------
        aux             = 0.0
        PreviousArea    = 0.0
        
        if (Kbottom == KUB) then 
            T_Area(Kbottom) = pi*Radius**2
        else                
            do K = Kbottom, KUB
                aux= aux + DUZ(I,J,K)
                if (aux .gt. (H - Radius) .and. aux .lt. (H + Radius)) then
                    dist = abs(H-aux)
                    Theta = 2*acos(dist/Radius)
                    if (aux .gt. H) then
                        q =  - 1
                        Beta = 2*pi - Theta    
                    else
                        q =  1
                        Beta = Theta
                    end if
                    T_Area(K) = ((Radius**2)/2.0)*(Beta) - q*Radius*sin(Beta/2)*dist - PreviousArea
                    PreviousArea = PreviousArea + T_Area(K)
                else if (aux .gt. (H + Radius)) then
                    T_Area(K) = pi*Radius**2 - PreviousArea   
                    exit
                end if
            end do
          end if  
        
    end subroutine TurbineVerticalDiscretisation
    
   
    
    Subroutine OutPut_Turbine(TurbineID, STAT)
        !Arguments-------------------------------------------------------------
        integer                         :: TurbineID
        integer, optional, intent (OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                         :: STAT_, STAT_CALL, ready_
        !Begin-----------------------------------------------------------------
        STAT_ = UNKNOWN_
    
        call Ready(TurbineID, ready_)
        
        if (ready_ .EQ. IDLE_ERR_) then
            call WriteTimeSerie(Me%ObjTimeSerie,                        &
                                Data1D = Me%Power,                      &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Turbine - ModuleTurbine - ERR21'
        
            call WriteTimeSerie(Me%ObjTimeSerie,                        &
                              Data1D = Me%Energy,                      &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Turbine - ModuleTurbine - ERR22'
        
            call WriteTimeSerie(Me%ObjTimeSerie,                        &
                                Data1D = Me%TurbineVelocity,                      &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Turbine - ModuleTurbine - ERR23'        
         else 
            STAT_ = ready_
        end if
    
        if (present(STAT)) STAT = STAT_
       
    end subroutine OutPut_Turbine
    
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
                if (Me%TimeSerieON) then
                    call KillTimeSerie(Me%ObjTimeSerie)
                end if
                call Deassociate_External_Modules

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
    
    subroutine Deassociate_External_Modules ()

        !Arguments-------------------------------------------------------------

        !Local---------------------------------------------------------
        integer                       :: nUsers

        nUsers = DeassociateInstance (mHORIZONTALGRID_,           Me%ObjHorizontalGrid)
        if (nUsers == 0) stop 'KillHydrodynamic - ModuleHydrodynamic - ERR24'

        nUsers = DeassociateInstance (mGEOMETRY_,       Me%ObjGeometry)
        if (nUsers == 0) stop 'KillHydrodynamic - ModuleHydrodynamic - ERR25'
        
        nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
        if (nUsers == 0) stop 'KillHydrodynamic - ModuleHydrodynamic - ERR26'
    
    
    end subroutine Deassociate_External_Modules
    
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

        if (.not. associated(Me)) stop 'ModuleShell - LocateObjShell - ERR26'

    end subroutine LocateObjShell

end module ModuleTurbine