!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : FreeVerticalMovement
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to compute particulate properties transport due to settling velocity
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

!DataFile
!
!<beginproperty>
!   NAME                        : char              -           !Property name
!   UNITS                       : char              -           !Property units
!   DESCRIPTION                 : char              -           !Small description of the property
!   CHS                         : real            [4 g/l]       !Hindered settling threshold concentration  
!   WS_VALUE                    : real          [0.0001 m/s]    !Constant settling velocity
!   WS_TYPE                     : int               [1]         !Settling velocity compute method 
!                                                               !Constant = 1, SPMFunction = 2, WSSecondaryClarifier = 3
!   SALTINT                     : 0/1               [0]         !Use salinity effect on settling velocity
!   SALTINTVALUE                : real            [3 psu]       !Salinity threshold concentration for affecting
!                                                               !settling velocity
!   KL                          : real            [0.006]       !Parameter to compute settling velocity 
!   KL1                         : real             [0.1]        !Parameter to compute settling velocity
!   ML                          : real            [4.65]        !Parameter to compute settling velocity
!   M                           : real              [1]         !Parameter to compute settling velocity
!   FREEVERT_IMPEXP_ADV         : real              [0]         !Time discretization method to compute vertical
!                                                               !transport: 0 - Implicit ; 1 - Explicit
!   DEPOSITION                  : 0/1               [1]         !Property can be deposited in the bottom
!   WITH_COMPRESSION            : 0/1               [1]         !Settling velocity with the compression effect
!<endproperty>

Module ModuleFreeVerticalMovement

    use ModuleGlobalData
    use ModuleFunctions,        only: THOMASZ, ConstructPropertyID, SettlingVelocity,          &
                                      SetMatrixValue, CHUNK_J, CHUNK_K, T_VECGW,                &
                                      T_THOMAS, T_D_E_F, Pad, SettlingVelSecondaryClarifier
    use ModuleTime
    use ModuleHorizontalGrid,   only: GetGridCellArea, UngetHorizontalGrid
    use ModuleGeometry,         only: GetGeometrySize, GetGeometryVolumes, UnGetGeometry,      &
                                      GetGeometryKFloor
    use ModuleMap,              only: GetOpenPoints3D, GetLandPoints3D, UngetMap
    use ModuleEnterData,        only: ReadFileName, ConstructEnterData, GetData, Block_Unlock, &
                                      ExtractBlockFromBuffer, KillEnterData
    use ModuleStopWatch,            only: StartWatch, StopWatch

#ifdef _ENABLE_CUDA
    use ModuleCuda
#endif _ENABLE_CUDA

    !griflet: openmp
    !$ use omp_lib

    implicit none 

    private

    !Subroutines & Functions---------------------------------------------------

    !Constructor
    public  ::  Construct_FreeVerticalMovement
    private ::      AllocateInstance
    private ::      AllocateVariables
    private ::      ReadFreeVertMovFilesName
    private ::      Construct_PropertyList
    private ::          Construct_Property  
    private ::              Construct_PropertyParameters
    private ::          Add_Property


    !Modifier
    public  :: Modify_FreeVerticalMovement
    private ::     FreeVerticalMovementIteration
    private ::         VerticalFreeConvection
    private ::         CalcVerticalFreeConvFlux
    private ::         Vertical_Velocity

    !Selector 
    public  :: Get_FreeVelocity
    public  :: Get_FreeConvFlux
    public  :: GetFreeVertMovOptions
    public  :: FreeVertPropertyExists                   !Function
    public  :: FreeVertPropertyHasDeposition            !Function
    public  :: UngetFreeVerticalMovement

    public  :: SetDepositionProbability

    !Destructor
    public  :: Kill_FreeVerticalMovement
    private ::     DeallocateInstance

    !Management
    private ::     Ready
    private ::         LocateObjFreeVerticalMovement
    
    !Parameters
    character(LEN = StringLength), parameter    :: prop_block_begin    = '<beginproperty>'
    character(LEN = StringLength), parameter    :: prop_block_end      = '<endproperty>'


    !Types---------------------------------------------------------------------
    type       T_Property
        type (T_PropertyID)                     :: ID
        real                                    :: Ws_Value         = FillValueReal
        real                                    :: CHS              = FillValueReal
        real                                    :: KL               = FillValueReal
        real                                    :: KL1              = FillValueReal     
        real                                    :: M                = FillValueReal
        real                                    :: ML               = FillValueReal
        real                                    :: ImpExp_AdvV      = FillValueReal
        logical                                 :: SalinityEffect   = .true.
        real                                    :: SalinityLimit    = FillValueReal
        integer                                 :: Ws_Type          = SPMFunction
        logical                                 :: Deposition       = .true.
        logical                                 :: WithCompression  = .true.
        real                                    :: SVI              = FillValueReal
        real                                    :: Clarification    = FillValueReal
        real, pointer, dimension(:,:,:)         :: FreeConvFlux     => null()
        real, pointer, dimension(:,:,:)         :: Velocity         => null()
        type(T_Property), pointer               :: Next             => null(), &
                                                   Prev             => null()
    end type T_Property

    type       T_DEF
        real,    pointer, dimension(: , : , :)  :: D    => null()
        real(8), pointer, dimension(: , : , :)  :: E    => null()
        real,    pointer, dimension(: , : , :)  :: F    => null()
#ifdef _USE_PAGELOCKED
        type(C_PTR)                             :: DPtr
        type(C_PTR)                             :: EPtr
        type(C_PTR)                             :: FPtr
#endif _USE_PAGELOCKED
        real,    pointer, dimension(: , : , :)  :: D_flux   => null()    !Coeficient to calculate ConvFlux and DifFlux
        real,    pointer, dimension(: , : , :)  :: E_flux   => null()    !Coeficient to calculate ConvFlux and DifFlux
    end type T_DEF

    type       T_External
        integer, pointer, dimension(: , : , :)  :: LandPoints       => null()
        integer, pointer, dimension(: , : , :)  :: OpenPoints3D     => null()
        real, pointer, dimension   (: , : , :)  :: Concentration    => null()
#ifdef _USE_PAGELOCKED
        type(C_PTR)                             :: ConcentrationPtr
#endif
        real, pointer, dimension   (: , : , :)  :: SPM                      => null()
        real, pointer, dimension   (: , : , :)  :: SalinityField            => null()
        real,    pointer, dimension(: , :    )  :: DepositionProbability    => null()
        integer, pointer, dimension(: , :    )  :: KFloor_Z                 => null()
        type(T_Time)                            :: Now
        real                                    :: DTProp           = FillValueReal
        real                                    :: IS_Coef          = FillValueReal
        real                                    :: SPMISCoef        = FillValueReal
        logical                                 :: NoFlux           = .false. !initialization: Jauch
        integer, pointer, dimension(: , : , :)  :: NoFluxW          => null()
    end type T_External

    type       T_Options
        logical                                 :: Salinity = .false. !initialization: Jauch 
        logical                                 :: SPM      = .false. !initialization: Jauch
    end type   T_Options

    type      T_FreeVerticalMovement
        private
        integer                                 :: InstanceID   = null_int !initialization: Jauch
        type(T_Size3D  )                        :: Size
        type(T_Size3D  )                        :: WorkSize
        type(T_External)                        :: ExternalVar
        type(T_DEF   )                          :: COEF3
        type(T_Options )                        :: Needs
        real, pointer, dimension(:,:,:)         :: TICOEF3  => null()
#ifdef _USE_PAGELOCKED
        type(C_PTR)                             :: TICOEF3Ptr
#endif _USE_PAGELOCKED
                                                
        !Auxiliar thomas arrays                 
        real(8), pointer, dimension(:)          :: VECG => null()
        real(8), pointer, dimension(:)          :: VECW => null()
        
        !griflet, openmp
        type(T_THOMAS), pointer                 :: THOMAS       => null()
        integer                                 :: MaxThreads   = null_int !initialization: Jauch

        type(T_Property), pointer               :: FirstProperty    => null()
        type(T_Property), pointer               :: LastProperty     => null()
        integer                                 :: PropertiesNumber = null_int !initialization: Jauch

        character(LEN = PathLength)             :: FileName         = null_str !initialization: Jauch
        
        !Instance of ModuleTime
        integer                                 :: ObjEnterData         = 0

        !Instance of ModuleTime
        integer                                 :: ObjTime              = 0
                                                                        
        !Instance of ModuleHorizontalGrid                               
        integer                                 :: ObjHorizontalGrid    = 0
                                                                        
        !Instance of ModuleGeometry                                     
        integer                                 :: ObjGeometry          = 0
                                                                        
        !Instance of ModuleMap                                          
        integer                                 :: ObjMap               = 0

#ifdef _ENABLE_CUDA        
        !Instance of ModuleCuda
        integer                                 :: ObjCuda              = 0
#endif _ENABLE_CUDA

        !Collection of instances
        type(T_FreeVerticalMovement), pointer   :: Next => null()

    end type T_FreeVerticalMovement

    !Global Module Variables
    type (T_FreeVerticalMovement), pointer      :: FirstObjFreeVerticalMovement => null()
    type (T_FreeVerticalMovement), pointer      :: Me                           => null()


    !--------------------------------------------------------------------------


    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




    Subroutine Construct_FreeVerticalMovement(FreeVerticalMovementID,         &
                                              TimeID, HorizontalGridID,       &
                                              MapID, GeometryID,              &
#ifdef _ENABLE_CUDA
                                              ObjCudaID,                      &
#endif    
                                              STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: FreeVerticalMovementID
        integer                                     :: TimeID
        integer                                     :: HorizontalGridID
        integer                                     :: MapID
        integer                                     :: GeometryID
#ifdef _ENABLE_CUDA
        integer                                     :: ObjCudaID
#endif
        integer, optional, intent(OUT)              :: STAT                
        
        !External--------------------------------------------------------------
        integer                                     :: ready_, STAT_CALL                        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_
 
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_


        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mFreeVerticalMovement_)) then
            nullify (FirstObjFreeVerticalMovement)
            call RegisterModule (mFreeVerticalMovement_) 
        endif

        call Ready(FreeVerticalMovementID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then
                               
            !Allocates a new Instance
            call AllocateInstance

            !Associates External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjMap            = AssociateInstance (mMAP_,            MapID           )
#ifdef _ENABLE_CUDA
            Me%ObjCuda           = AssociateInstance (mCUDA_,           ObjCudaID       )
#endif _ENABLE_CUDA

            call GetGeometrySize(Me%ObjGeometry,                                        &
                                 Size       = Me%Size,                                  &
                                 WorkSize   = Me%WorkSize,                              &
                                 STAT       = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)                                                 &
                stop 'Construct_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR01'

            call ReadFreeVertMovFilesName

            call ConstructEnterData(Me%ObjEnterData, Me%FileName, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)                                                 &
                stop 'Construct_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR02'

            call AllocateVariables

            call Construct_PropertyList

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)                                                 &
                stop 'Construct_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR03'
            
            STAT_ = SUCCESS_

            FreeVerticalMovementID = Me%InstanceID
        
        else cd0
            
            stop 'Construct_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR99' 

        end if cd0

        if (present(STAT)) STAT = STAT_

    end subroutine Construct_FreeVerticalMovement

    !--------------------------------------------------------------------------


    subroutine AllocateInstance

        !Local-----------------------------------------------------------------
        type (T_FreeVerticalMovement), pointer           :: NewFreeVerticalMovement
        type (T_FreeVerticalMovement), pointer           :: PreviousFreeVerticalMovement


        !Allocates new instance
        allocate (NewFreeVerticalMovement)
        nullify  (NewFreeVerticalMovement%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjFreeVerticalMovement)) then
            FirstObjFreeVerticalMovement        => NewFreeVerticalMovement
            Me                                  => NewFreeVerticalMovement
        else
            PreviousFreeVerticalMovement        => FirstObjFreeVerticalMovement
            Me                                  => FirstObjFreeVerticalMovement%Next
            do while (associated(Me))
                PreviousFreeVerticalMovement    => Me
                Me                              => Me%Next
            enddo
            Me                                  => NewFreeVerticalMovement
            PreviousFreeVerticalMovement%Next   => NewFreeVerticalMovement
        endif

        Me%InstanceID   = RegisterNewInstance (mFREEVERTICALMOVEMENT_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine AllocateVariables
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                 :: ILB, IUB 
        integer                                 :: JLB, JUB 
        integer                                 :: KLB, KUB 
        
        !griflet: openmp
        integer                                 :: m
        type(T_VECGW), pointer                  :: VECGW

        !----------------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB   
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB   
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB
        
        nullify(Me%COEF3%D) 
        nullify(Me%COEF3%E)
        nullify(Me%COEF3%F)
        nullify(Me%COEF3%D_flux)    !Coeficient to calculate ConvFlux and DifFlux
        nullify(Me%COEF3%E_flux)    !Coeficient to calculate ConvFlux and DifFlux
        
        nullify(Me%ExternalVar%LandPoints)
        nullify(Me%ExternalVar%Concentration)
        nullify(Me%ExternalVar%SPM)
        nullify(Me%ExternalVar%SalinityField)

        nullify(Me%TICOEF3)       

#ifdef _USE_PAGELOCKED
        ! Allocate pagelocked memory to optimize CUDA transfers
        call Alloc3DPageLocked(Me%ObjCuda, Me%COEF3%DPtr, Me%COEF3%D, IUB + 1, JUB + 1, KUB + 1)        
        call Alloc3DPageLocked(Me%ObjCuda, Me%COEF3%EPtr, Me%COEF3%E, IUB + 1, JUB + 1, KUB + 1)        
        call Alloc3DPageLocked(Me%ObjCuda, Me%COEF3%FPtr, Me%COEF3%F, IUB + 1, JUB + 1, KUB + 1)        
#else
        allocate(Me%COEF3%D      (ILB:Pad(ILB, IUB), JLB:JUB, KLB:KUB), STAT = STAT_CALL)         
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleFreeVerticalMovement - ERR01'
        Me%COEF3%D = FillValueReal


        allocate(Me%COEF3%E      (ILB:Pad(ILB, IUB), JLB:JUB, KLB:KUB), STAT = STAT_CALL)         
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleFreeVerticalMovement - ERR02'
         Me%COEF3%E = FillValueReal


        allocate(Me%COEF3%F      (ILB:Pad(ILB, IUB), JLB:JUB, KLB:KUB), STAT = STAT_CALL)         
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleFreeVerticalMovement - ERR03'
        Me%COEF3%F = FillValueReal
#endif _USE_PAGELOCKED

        allocate(Me%COEF3%D_flux (ILB:Pad(ILB, IUB), JLB:JUB, KLB:KUB), STAT = STAT_CALL)         
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleFreeVerticalMovement - ERR04'
         Me%COEF3%D_flux = FillValueReal


        allocate(Me%COEF3%E_flux (ILB:Pad(ILB, IUB), JLB:JUB, KLB:KUB), STAT = STAT_CALL)         
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleFreeVerticalMovement - ERR05'
         Me%COEF3%E_flux = FillValueReal

#ifdef _USE_PAGELOCKED
        ! Allocate pagelocked memory to optimize CUDA transfers
        call Alloc3DPageLocked(Me%ObjCuda, Me%TICOEF3Ptr, Me%TICOEF3, IUB + 1, JUB + 1, KUB + 1)        
#else
        allocate(Me%TICOEF3      (ILB:Pad(ILB, IUB), JLB:JUB, KLB:KUB), STAT = STAT_CALL)         
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleFreeVerticalMovement - ERR06'
         Me%TICOEF3        = FillValueReal
#endif _USE_PAGELOCKED

        allocate(Me%VECG        (                  KLB:KUB), STAT = STAT_CALL)         
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleFreeVerticalMovement - ERR07'
        Me%VECG = FillValueReal


        allocate(Me%VECW        (                  KLB:KUB), STAT = STAT_CALL)         
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleFreeVerticalMovement - ERR08'
        Me%VECW = FillValueReal

        !griflet: BEGIN this is the alternate version that allows parallel openmp

        Me%MaxThreads = 1
        !$ Me%MaxThreads = omp_get_max_threads()

        allocate(Me%THOMAS)
        allocate(Me%THOMAS%COEF3)
        allocate(Me%THOMAS%VEC(1:Me%MaxThreads))

        do m = 1, Me%MaxThreads

            VECGW => Me%THOMAS%VEC(m)

            allocate(VECGW%G(KLB:KUB))
            allocate(VECGW%W(KLB:KUB))

        enddo

        Me%THOMAS%COEF3%D => Me%COEF3%D
        Me%THOMAS%COEF3%E => Me%COEF3%E
        Me%THOMAS%COEF3%F => Me%COEF3%F
        Me%THOMAS%TI => Me%TICOEF3

        !griflet: END

        !----------------------------------------------------------------------

    end subroutine AllocateVariables   
    
    !----------------------------------------------------------------------
   
    subroutine Construct_PropertyList

        !External----------------------------------------------------------------
        integer                         :: ClientNumber
        integer                         :: STAT_CALL
        logical                         :: BlockFound

        !Local-------------------------------------------------------------------
        type (T_Property), pointer      :: NewProperty

        !------------------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = prop_block_begin,     &
                                        block_end       = prop_block_end,       &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :           if (BlockFound) then                                                  
                    ! Construct a New Property 
                    Call Construct_Property(NewProperty)

                    ! Add new Property to the WaterProperties List 
                    Call Add_Property(NewProperty)
                else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Construct_PropertyList - ModuleFreeVerticalMovement - ERR01'

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'Construct_PropertyList - ModuleFreeVerticalMovement - ERR02'
            else cd1
                stop 'Construct_PropertyList - ModuleFreeVerticalMovement - ERR03'
            end if cd1
        end do do1

        !------------------------------------------------------------------------

    end subroutine Construct_PropertyList

    !----------------------------------------------------------------------

    subroutine Construct_Property(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer       :: NewProperty

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL

        !Local------------------------------------------------------------------
        integer                                 :: ILB, IUB
        integer                                 :: JLB, JUB
        integer                                 :: KLB, KUB
        
        !----------------------------------------------------------------------
        ILB = Me%Size%ILB   
        IUB = Me%Size%IUB   
        JLB = Me%Size%JLB   
        JUB = Me%Size%JUB   
        KLB = Me%Size%KLB   
        KUB = Me%Size%KUB   
     
        allocate (NewProperty, STAT = STAT_CALL)            
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'Construct_Property - ModuleFreeVerticalMovement - ERR01' 


        nullify(NewProperty%Velocity             )
        nullify(NewProperty%FreeConvFlux         )
        nullify(NewProperty%Prev,NewProperty%Next)


        allocate(NewProperty%Velocity     (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)         
        if (STAT_CALL .NE. SUCCESS_)&
            stop 'Construct_Property - ModuleFreeVerticalMovement - ERR02'
        NewProperty%Velocity       = 0.

        allocate(NewProperty%FreeConvFlux (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)         
        if (STAT_CALL .NE. SUCCESS_)&
            stop 'Construct_Property - ModuleFreeVerticalMovement - ERR03'
        NewProperty%FreeConvFlux   = 0.  

        call ConstructPropertyID            (NewProperty%ID, Me%ObjEnterData, FromBlock)  

        !Construct property values
        call Construct_PropertyParameters   (NewProperty)

    end subroutine Construct_Property
    
    !----------------------------------------------------------------------
    
    subroutine Construct_PropertyParameters(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer       :: NewProperty

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                         :: iflag

        !<BeginKeyword>
            !Keyword          : CHS
            !<BeginDescription>       
               ! 
               ! Hindered settling  - CHS
               !
            !<EndDescription>
            !Type             : Real 
            !Default          : 4.0
            !File keyword     : FREE_DAT
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%CHS,                                           &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'CHS',                                      &
                     default      = 4.0,                                        &
                     ClientModule = 'ModuleFreeVerticalMovement',               &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR10'        

        !<BeginKeyword>
            !Keyword          : WS_VALUE
            !<BeginDescription>       
               ! 
               ! Constant settling velocity
               !
            !<EndDescription>
            !Type             : Real 
            !Default          : 0.0001    !m/s
            !File keyword     : FREE_DAT
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%Ws_Value,                                      &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'WS_VALUE',                                 &
                     default      = 0.0001,                                     &
                     ClientModule = 'ModuleFreeVerticalMovement',               &
                     STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR20'        

        if (NewProperty%Ws_Value > 0.) NewProperty%Ws_Value  = - NewProperty%Ws_Value

        !<BeginKeyword>
            !Keyword          : WS_TYPE
            !<BeginDescription>       
               ! 
               ! Settling type
               !
            !<EndDescription>
            !Type             : integer 
            !Default          : WSConstant
            !File keyword     : FREE_DAT
            !Multiple Options : WSConstant = 1, SPMFunction = 2, WSSecondaryClarifier = 3
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%Ws_Type,                                       &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'WS_TYPE',                                  &
                     default      = WSConstant,                                 &
                     ClientModule = 'ModuleFreeVerticalMovement',               &
                     STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR40' 
           
           
        if(NewProperty%Ws_Type == SPMFunction          ) Me%Needs%SPM = .true.       
        
        if(NewProperty%Ws_Type == WSSecondaryClarifier ) Me%Needs%SPM = .true.              
        
        if(NewProperty%Ws_Type == WSConstant)then
            call SetMatrixValue(NewProperty%Velocity, Me%Size, NewProperty%Ws_Value)
        endif
        
        !<BeginKeyword>
            !Keyword          : SALTINT
            !<BeginDescription>       
               ! 
               ! Free vertical movement is function of the salinity ? 
               !
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : FREE_DAT
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%SalinityEffect,                                &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'SALTINT',                                  &
                     Default      = .false.,                                    &
                     ClientModule = 'ModuleFreeVerticalMovement',               &
                     STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR50'

        if(NewProperty%SalinityEffect) Me%Needs%Salinity = .true.
        
        
        If (NewProperty%SalinityEffect) then

            !<BeginKeyword>
                !Keyword          : SALTINTVALUE
                !<BeginDescription>       
                   ! Salinity limit. For salinity values smaller the settling velocity is zero
                   ! For salinity values greater then this limit the settling velocity is compute using 
                   ! the Krone formulation
                !<EndDescription>
                !Type             : Real 
                !Default          : 3.0
                !File keyword     : FREE_DAT
                !Multiple Options : Do not have
                !Search Type      : FromBlock
                !Begin Block      : <beginproperty>
                !End Block        : <endproperty>
            !<EndKeyword>
            call GetData(NewProperty%SalinityLimit,                             &
                         Me%ObjEnterData, iflag,                                &
                         SearchType   = FromBlock,                              &
                         keyword      = 'SALTINTVALUE',                         &
                         Default      = 3.0,                                    &
                         ClientModule = 'ModuleFreeVerticalMovement',           &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR60'        
        
        endif

        ! KL

        !<BeginKeyword>
            !Keyword          : KL
            !<BeginDescription>       
               ! 
               ! Do not have
               ! 
            !<EndDescription>
            !Type             : Real 
            !Default          : 0.006
            !File keyword     : FREE_DAT
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%KL,                                            &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'KL',                                       &
                     Default      = 0.006,                                      &
                     ClientModule = 'ModuleFreeVerticalMovement',               &
                     STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR70'        

        !<BeginKeyword>
            !Keyword          : KL1
            !<BeginDescription>       
               ! 
               ! Do not have
               ! 
            !<EndDescription>
            !Type             : Real 
            !Default          : 0.1
            !File keyword     : FREE_DAT
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%KL1,                                           &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'KL1',                                      &  
                     Default      = 0.1,                                        &
                     ClientModule = 'ModuleFreeVerticalMovement',               &
                     STAT       = STAT_CALL)  
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR80'        

          
        !<BeginKeyword>
            !Keyword          : ML
            !<BeginDescription>       
               ! 
               ! Do not have
               ! 
            !<EndDescription>
            !Type             : Real 
            !Default          : 4.65
            !File keyword     : FREE_DAT
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%ML,                                            &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'ML',                                       &
                     Default      = 4.65,                                       & 
                     ClientModule = 'ModuleFreeVerticalMovement',               &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR90'        


        !<BeginKeyword>
            !Keyword          : M
            !<BeginDescription>       
               ! 
               ! Do not have
               ! 
            !<EndDescription>
            !Type             : Real 
            !Default          : 1.0 
            !File keyword     : FREE_DAT
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%M,                                             &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'M',                                        &
                     Default      = 1.0,                                        &
                     ClientModule = 'ModuleFreeVerticalMovement',               &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR100'        


        !<BeginKeyword>
            !Keyword          : FREEVERT_IMPEXP_ADV
            !<BeginDescription>       
               ! 
               ! Do not have
               ! 
            !<EndDescription>
            !Type             : Real 
            !Default          : 0.0
            !File keyword     : FREE_DAT
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%ImpExp_AdvV,                                   &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'FREEVERT_IMPEXP_ADV',                      &
                     Default      = 0.0,                                        &
                     ClientModule = 'ModuleFreeVerticalMovement',               &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR110'  
            
        !<BeginKeyword>
            !Keyword          : DEPOSITION
            !<BeginDescription>       
               ! 
               ! Specifies if the bottom boundary is considered to be closed
               ! 
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .true.
            !File keyword     : FREE_DAT
            !Multiple Options : 0 - false, 1 - true
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%Deposition,                                    &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'DEPOSITION',                               &
                     Default      = .true.,                                     &
                     ClientModule = 'ModuleFreeVerticalMovement',               &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR120' 
            
        call GetData(NewProperty%WithCompression,                               &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'WITH_COMPRESSION',                         &
                     Default      = .true.,                                     &
                     ClientModule = 'ModuleFreeVerticalMovement',               &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR130'

        call GetData(NewProperty%SVI,                                           &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'SVI',                                      &
                     Default      = 120.,                                       &
                     ClientModule = 'ModuleFreeVerticalMovement',               &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR140'
            
        if (NewProperty%SVI < 10 .or. NewProperty%SVI > 500) then
            write(*,*) 'SVI = ', NewProperty%SVI
            write(*,*) 'out of bounds, min. 10 and máx. 500'
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR140'
        endif

        call GetData(NewProperty%Clarification,                                 &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'CLARIFICATION',                            &
                     Default      = 1.,                                         &
                     ClientModule = 'ModuleFreeVerticalMovement',               &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR150'
            
        if (NewProperty%Clarification < 0.1 .or. NewProperty%Clarification > 1) then
            write(*,*) 'CLARIFICATION = ', NewProperty%Clarification
            write(*,*) 'out of bounds, min. 0.1 and máx. 1'
            stop 'Read_FreeVert_Parameters - ModuleFreeVerticalMovement - ERR160'
        endif
            

            
    end subroutine Construct_PropertyParameters

    !----------------------------------------------------------------------

    subroutine Add_Property(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),              pointer     :: NewProperty

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstProperty)) then
            Me%PropertiesNumber     = 1
            Me%FirstProperty        => NewProperty
            Me%LastProperty         => NewProperty
        else
            NewProperty%Prev        => Me%LastProperty
            Me%LastProperty%Next    => NewProperty
            Me%LastProperty         => NewProperty
            Me%PropertiesNumber     = Me%PropertiesNumber + 1
        end if 

    end subroutine Add_Property 


    !--------------------------------------------------------------------------

    
    subroutine ReadFreeVertMovFilesName

        !External--------------------------------------------------------------
        integer                     :: STAT_CALL

        call ReadFileName('FREE_DAT', Me%FileName, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) &
            stop 'ReadFreeVertMovFilesName - ModuleFreeVerticalMovement - ERR01'

                                                                             
    end subroutine ReadFreeVertMovFilesName

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     subroutine Modify_FreeVerticalMovement(FreeVerticalMovementID, Concentration, SPM,     &
                                            SPMISCoef, SalinityField, PropertyID, IS_Coef,  &
                                            DTProp, CurrentTime, NoFlux, NoFluxW, STAT)
       
        !Arguments-------------------------------------------------------------
        integer                                 :: FreeVerticalMovementID
        real, dimension(:,:,:), pointer         :: Concentration
        real, dimension(:,:,:), pointer         :: SPM
        real                                    :: SPMISCoef
        real, dimension(:,:,:), pointer         :: SalinityField  
        integer                                 :: PropertyID
        real                                    :: IS_Coef, DTProp
        type(T_Time)                            :: CurrentTime
        logical, intent(IN), optional           :: NoFlux
        integer, dimension(:,:,:), pointer, optional :: NoFluxW           
                
        integer, optional                       :: STAT
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: ready_
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_    
        type(T_Property), pointer               :: PropertyX
        
        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(FreeVerticalMovementID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then


            call Search_Property(PropertyX      = PropertyX,                            &
                                 PropertyXID    = PropertyID,                           &
                                 STAT           = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                                &
                stop 'Modify_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR02'
           
            Me%ExternalVar%Now           =  CurrentTime
            Me%ExternalVar%IS_Coef       =  IS_Coef
            Me%ExternalVar%DTProp        =  DTProp 
            Me%ExternalVar%Concentration => Concentration

            if(associated(SalinityField))then
                Me%ExternalVar%SalinityField => SalinityField
            end if

            if(PropertyX%Ws_Type == SPMFunction .or. PropertyX%Ws_Type == WSSecondaryClarifier )then
                Me%ExternalVar%SPMISCoef =  SPMISCoef
                Me%ExternalVar%SPM       => SPM
            end if
            
            if (present(NoFlux)) then
                Me%ExternalVar%NoFlux  =  NoFlux
                Me%ExternalVar%NoFluxW => NoFluxW
            else
                Me%ExternalVar%NoFlux  =  .false. 
            endif
                        
            !Compute free vertical movement
            call FreeVerticalMovementIteration(PropertyX)

            nullify(Me%ExternalVar%Concentration)
            nullify(Me%ExternalVar%SPM)
            nullify(Me%ExternalVar%SalinityField)
            nullify(Me%ExternalVar%NoFluxW)

            call null_time   (Me%ExternalVar%Now)
 
           
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
           

     end subroutine Modify_FreeVerticalMovement

    !--------------------------------------------------------------------------

    subroutine FreeVerticalMovementIteration(PropertyX)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer               :: PropertyX

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                 :: i, j, k, CHUNK 
        
        !----------------------------------------------------------------------

        !OpenPoints3D
        call GetOpenPoints3D(Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL .ne. SUCCESS_) stop 'FreeVerticalMovementIteration - ModuleFreeVerticalMovement - ERR01'

        !LandPoints
        call GetLandPoints3D(Me%ObjMap,Me%ExternalVar%LandPoints, STAT = STAT_CALL)
        if (STAT_CALL .ne. SUCCESS_) stop 'FreeVerticalMovementIteration - ModuleFreeVerticalMovement - ERR02'

        call SetMatrixValue(Me%COEF3%D, Me%Size,                          0.0 )
        call SetMatrixValue(Me%COEF3%E, Me%Size,                     dble(1.0))
        call SetMatrixValue(Me%COEF3%F, Me%Size,                          0.0 )
        call SetMatrixValue(Me%TICOEF3, Me%Size, Me%ExternalVar%Concentration )

        call SetMatrixValue(PropertyX%FreeConvFlux, Me%Size,              0.0 )        
              
        call Vertical_Velocity          (PropertyX)

        call BottomBoundary             (PropertyX)

        call VerticalFreeConvection     (PropertyX)

        !Explicit Fluxes among cells
        call CalcVerticalFreeConvFlux   (PropertyX, PropertyX%ImpExp_AdvV)


        ! At this stage the variable TICOEF3 have null values in the land points
        ! We want that all proprieties in land points have the value of FillValueReal.
        
        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)

        if (MonitorPerformance) call StartWatch ("ModuleFreeVerticalMovement", "FreeVerticalMovementIteration")

        !$OMP PARALLEL SHARED(CHUNK) PRIVATE(I,J,K)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k=Me%WorkSize%KLB, Me%WorkSize%KUB 
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
           Me%TICOEF3(i,j,k) = Me%TICOEF3(i,j,k)                                *    &
                               (1. - real(Me%ExternalVar%LandPoints(i,j,k)))    +    &
                                real(Me%ExternalVar%LandPoints(i,j,k))          *    &
                                FillValueReal
        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleFreeVerticalMovement", "FreeVerticalMovementIteration")

        if (PropertyX%ImpExp_AdvV/=1.) then  
            !Inversao do sistema de equacoes pelo algoritmo de thomas                             
            !griflet: old call
            !CALL THOMASZ(Me%WorkSize%ILB, Me%WorkSize%IUB,  &
            !              Me%WorkSize%JLB, Me%WorkSize%JUB,  &
            !              Me%WorkSize%KLB, Me%WorkSize%KUB,  &
            !              Me%COEF3%D,                        &
            !              Me%COEF3%E,                        &
            !              Me%COEF3%F,                        &
            !              Me%TICOEF3,                        &
            !              Me%ExternalVar%Concentration,      &
            !              Me%VECG,                           &
            !              Me%VECW)
            !griflet: new call
            CALL THOMASZ(Me%WorkSize%ILB, Me%WorkSize%IUB,  &
                         Me%WorkSize%JLB, Me%WorkSize%JUB,  &
                         Me%WorkSize%KLB, Me%WorkSize%KUB,  &
                         Me%THOMAS,                         &
                         Me%ExternalVar%Concentration       &
#ifdef _ENABLE_CUDA
            ! Use CUDA to solve Thomas
                         , Me%ObjCuda,                      &
                         .FALSE.                            &
#endif _ENABLE_CUDA
                         )

        else ! Explicit
        
            call SetMatrixValue(Me%ExternalVar%Concentration, Me%Size, Me%TICOEF3)
        
        endif

        !Implicit Fluxes among cells 
        call CalcVerticalFreeConvFlux(PropertyX, 1.0 - PropertyX%ImpExp_AdvV)


        !OpenPoints3D           
        call UnGetMap(Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FreeVerticalMovementIteration - ModuleFreeVerticalMovement - ERR04'
        
        !LandPoints
        call UngetMap(Me%ObjMap, Me%ExternalVar%LandPoints, STAT = STAT_CALL)
        if (STAT_CALL .ne. SUCCESS_) stop 'FreeVerticalMovementIteration - ModuleFreeVerticalMovement - ERR05'

    end subroutine FreeVerticalMovementIteration

    !----------------------------------------------------------------------

    subroutine VerticalFreeConvection(PropertyX)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer               :: PropertyX

        !Local-----------------------------------------------------------------               
        real,    dimension(:, :   ), pointer    :: GridCellArea
        real(8), dimension(:, :, :), pointer    :: VolumeZ
        real(8)                                 :: W1, AW1, W2, AW2                   
        real(8)                                 :: DT_V 
        real(8)                                 :: D_flux, coef_D
        real(8)                                 :: E_flux, coef_E
        real(8)                                 :: coef_F
        real                                    :: UPDC_Conv
        integer                                 :: I, J, K, STAT_CALL

        !----------------------------------------------------------------------

        !Gets GetGridCellArea
        nullify(GridCellArea)
        call GetGridCellArea(Me%ObjHorizontalGrid, GridCellArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalFreeConvection - ModuleFreeVerticalMovement - ERR01'

        !Gets VolumeZ
        nullify(VolumeZ)
        call GetGeometryVolumes(Me%ObjGeometry,                         &
                                VolumeZ     = VolumeZ,                  &
                                STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalFreeConvection - ModuleFreeVerticalMovement - ERR02'

        !The vertical free fluxes are computed always using upwind discretization
        UPDC_Conv = 1
   
        !If the property sediments the discretization should be implicit
        if (.not.(PropertyX%ImpExp_AdvV == 0.0 .or. PropertyX%ImpExp_AdvV == 1)) then 
            write(*,*)'The module ModuleFreeVerticalMovement is not prepared to run this option'   
            stop 'VerticalFreeConvection - ModuleFreeVerticalMovement - ERR04'
        endif
            
do3 :   do k=Me%WorkSize%KLB, Me%WorkSize%KUB 
do2 :   do j=Me%WorkSize%JLB, Me%WorkSize%JUB
do1 :   do i=Me%WorkSize%ILB, Me%WorkSize%IUB
                   
cd4 :       if (Me%ExternalVar%OpenPoints3D (i, j, Me%WorkSize%KUB) == OpenPoint) then

                ! Variaveis auxiliares
                DT_V = dble(Me%ExternalVar%DTProp) / VolumeZ(i, j, k) 

                W1   = PropertyX%Velocity(i, j, k) * GridCellArea (i,j)    

                if (K < Me%WorkSize%KUB .and. K > 0) then                                             

                    W2   = PropertyX%Velocity(i, j, k+1) * GridCellArea(i,j)

                else  if (k == Me%WorkSize%KUB) then

                    W2 = 0.0

                else

                    stop 'VerticalFreeConvection - ModuleFreeVerticalMovement - ERR05'

                endif


                AW1  = ABS(W1)
                AW2  = ABS(W2)

                D_flux =-        UPDC_Conv  * (W1 + AW1) / 2.0                                             !& 
                       !- ( 1.0 - UPDC_Conv) *  W1        * DWZ(I, J, K) / (DWZ(I, J, K-1) + DWZ(I, J, K))

                coef_D = D_flux * DT_V

                
                E_flux =-       UPDC_Conv  * (W1 - AW1) / 2.0                                              !&
                       !- (1.0 - UPDC_Conv) *  W1        * DWZ(I, J, K-1) / (DWZ(I, J, K) + DWZ(I, J, K-1))

                
                coef_E = (      UPDC_Conv * (W2 + AW2) / 2.0 + E_flux) * DT_V                               !&
                       !+ (1.0 - UPDC_Conv) *  W2        * DWZ(I, J, K+1) / (DWZ(I, J, K) + DWZ(I, J, K+1)) &


                coef_F = (       UPDC_Conv  * (W2 - AW2) / 2.0 )* DT_V                                            !&
                       !+  (1.0 - UPDC_Conv) *  W2        * DWZ(I, J, K) / (DWZ(I, J, K) + DWZ(I, J, K+1))) &


                Me%COEF3%D_flux(I,J,K) = D_flux
                Me%COEF3%E_flux(I,J,K) = E_flux


                !Coefficients 
cd2 :           if (PropertyX%ImpExp_AdvV == 0.0) then     !0 = Implicit                                             

                    Me%COEF3%D(i,j,k) = Me%COEF3%D(i,j,k) + coef_D
                         
                    Me%COEF3%E(i,j,k) = Me%COEF3%E(i,j,k) + coef_E
                         
                    Me%COEF3%F(i,j,k) = Me%COEF3%F(i,j,k) + coef_F

                endif cd2

                !Independent Term
cd3 :           if (PropertyX%ImpExp_AdvV == 1.0) then     !1 = Explicit    
                    Me%TICOEF3(i,j,k) = Me%TICOEF3(i,j,k)  &
                           - (  coef_D * Me%ExternalVar%Concentration(i,j,k-1)           &
                              + coef_E * Me%ExternalVar%Concentration(i,j,k)             &
                              + coef_F * Me%ExternalVar%Concentration(i,j,k+1)) 
        
                endif cd3
            end if cd4
            
        end do do1
        end do do2
        end do do3
       
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, GridCellArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalFreeConvection - ModuleFreeVerticalMovement - ERR06'

        call UnGetGeometry(Me%ObjGeometry, VolumeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalFreeConvection - ModuleFreeVerticalMovement - ERR08'

    end subroutine VerticalFreeConvection

    !--------------------------------------------------------------------------

    subroutine CalcVerticalFreeConvFlux(PropertyX, Weigth)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer           :: PropertyX
        real, intent(IN)                    :: Weigth !Refers to the weigth of Implicit-Explicit calculations

        !Local-----------------------------------------------------------------
        integer                             :: i, j, k  
        !----------------------------------------------------------------------
       
do3 :   do k=Me%WorkSize%KLB, Me%WorkSize%KUB 
do2 :   do j=Me%WorkSize%JLB, Me%WorkSize%JUB
do1 :   do i=Me%WorkSize%ILB, Me%WorkSize%IUB
           
            if (Me%ExternalVar%OpenPoints3D(i,j,k)== OpenPoint) then

                !PCL
                ! [g/s] = c* [m^3/s] * [g/m^3]
                PropertyX%FreeConvFlux(i, j, k) = PropertyX%FreeConvFlux(i, j, k) - Weigth                      &
                                        * (Me%COEF3%D_flux(i, j, k) * Me%ExternalVar%Concentration(i, j, k-1)   &
                                        +  Me%COEF3%E_flux(i, j, k) * Me%ExternalVar%Concentration(i, j, k  ))
           endif
        end do do1
        end do do2
        end do do3

    end subroutine CalcVerticalFreeConvFlux
    
    !--------------------------------------------------------------------------
    
    subroutine Vertical_Velocity(PropertyX)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer               :: PropertyX

        !Local-----------------------------------------------------------------
        real,    dimension(:,:,:), pointer      :: SPM
        real                                    :: CHS,KL,KL1,M,ML
        integer                                 :: I, J, K  
        real                                    :: SPMISCoef
        integer                                 :: CHUNK
        !----------------------------------------------------------------------

        SPM             => Me%ExternalVar%SPM
        SPMISCoef       =  Me%ExternalVar%SPMISCoef

        CHS             =  PropertyX%CHS
        KL              =  PropertyX%KL
        KL1             =  PropertyX%KL1
        M               =  PropertyX%M
        ML              =  PropertyX%ML
                
        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        if (MonitorPerformance) call StartWatch ("ModuleFreeVerticalMovement", "Vertical_Velocity")
        
        if(PropertyX%Ws_Type == SPMFunction)then
        
            if(PropertyX%SalinityEffect)then
                
                !ACanas: Function call (SettlingVelocity) in cycle iteration leads to 
                !ACanas: parallelization error.
                
                do k=Me%WorkSize%KLB, Me%WorkSize%KUB 
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%ExternalVar%OpenPoints3D(i, j, k)== OpenPoint) then
                        if (Me%ExternalVar%SalinityField(i, j, k) .gt. PropertyX%SalinityLimit) then
                            PropertyX%Velocity(i, j, k)=-SettlingVelocity (SPM(i, j, k)*SPMISCoef,  &
                                                                           CHS,KL,KL1,M,ML,I,J,K)
                        else
                            PropertyX%Velocity(i, j, k)=0.
                        endif

                    end if
                enddo
                enddo
                enddo

            else
                do k=Me%WorkSize%KLB, Me%WorkSize%KUB 
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%ExternalVar%OpenPoints3D(i, j, k)== OpenPoint) then
                        PropertyX%Velocity(i, j, k)=-SettlingVelocity (SPM(i, j, k)*SPMISCoef,  &
                                                                       CHS,KL,KL1,M,ML,I,J,K)
                    end if
                enddo
                enddo
                enddo

            end if
            
        elseif(PropertyX%Ws_Type ==  WSSecondaryClarifier) then            
            
            do k=Me%WorkSize%KLB, Me%WorkSize%KUB 
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExternalVar%OpenPoints3D(i, j, k)== OpenPoint) then

                    PropertyX%Velocity(i, j, k)=-SettlingVelSecondaryClarifier (                &
                                                   Cx              = SPM(i, j, k)*SPMISCoef,    &
                                                   WithCompression = PropertyX%WithCompression, &
                                                   SVI             = PropertyX%SVI,             & 
                                                   Clarification   = PropertyX%Clarification)
                else
                    PropertyX%Velocity(i, j, k)= 0.                                            
                end if
            enddo
            enddo
            enddo
        
        elseif(PropertyX%Ws_Type ==  WSConstant) then
            
            if(PropertyX%SalinityEffect)then
                
                !$OMP PARALLEL SHARED(CHUNK, PropertyX, SPM, SPMISCoef, CHS,KL,KL1,M,ML) PRIVATE(I,J,K)
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do k=Me%WorkSize%KLB, Me%WorkSize%KUB 
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%ExternalVar%OpenPoints3D(i, j, k)== OpenPoint) then
                        if (Me%ExternalVar%SalinityField(i, j, k) .gt. PropertyX%SalinityLimit) then
                            PropertyX%Velocity(i, j, k) = PropertyX%WS_Value
                        else
                            PropertyX%Velocity(i, j, k) = 0.
                        endif
                    end if
                enddo
                enddo
                enddo
                !$OMP END DO NOWAIT
                !$OMP END PARALLEL
            else
                call SetMatrixValue(PropertyX%Velocity, Me%Size, PropertyX%Ws_Value)
            endif
        end if
        
        if (Me%ExternalVar%NoFlux) then
           !$OMP PARALLEL SHARED(CHUNK, PropertyX) PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k=Me%WorkSize%KLB, Me%WorkSize%KUB 
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExternalVar%OpenPoints3D(i, j, k) == OpenPoint .and.             &
                    Me%ExternalVar%NoFluxW     (i, j, k) == 1) then
                        PropertyX%Velocity(i, j, k) = 0.
                end if
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL        
        endif        
        
        if (MonitorPerformance) call StopWatch ("ModuleFreeVerticalMovement", "Vertical_Velocity")
        
        nullify(SPM)

    end subroutine Vertical_Velocity

    !--------------------------------------------------------------------------
    
    subroutine BottomBoundary(PropertyX)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer               :: PropertyX

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                 :: i, j, Kbottom 
        integer                                 :: ILB, IUB
        integer                                 :: JLB, JUB
        integer                                 :: KUB
        integer                                 :: CHUNK
        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KUB = Me%WorkSize%KUB

        !Gets KFloorZ
        call GetGeometryKFloor(GeometryID   = Me%ObjGeometry,                   &
                               Z            = Me%ExternalVar%KFloor_Z,          &
                               STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'BottomBoundary - ModuleFreeVerticalMovement - ERR01'

        CHUNK = CHUNK_J(JLB, JUB)

        if (MonitorPerformance) call StartWatch ("ModuleFreeVerticalMovement", "BottomBoundary")
        
        if(PropertyX%Deposition)then

            !$OMP PARALLEL SHARED(CHUNK, PropertyX) PRIVATE(I,J,Kbottom)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j=JLB, JUB
            do i=ILB, IUB
                
                if (Me%ExternalVar%OpenPoints3D (i,j, KUB) == OpenPoint) then

                    Kbottom = Me%ExternalVar%KFloor_Z(i, j)

                    PropertyX%Velocity(i,j,Kbottom) = PropertyX%Velocity(i,j,Kbottom) *     &
                                                      Me%ExternalVar%DepositionProbability(i,j)

                endif
                
            end do
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL

        else

            !$OMP PARALLEL SHARED(CHUNK, PropertyX) PRIVATE(I,J,Kbottom)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j=JLB, JUB
            do i=ILB, IUB
                
                if (Me%ExternalVar%OpenPoints3D (i,j, KUB) == OpenPoint) then

                    Kbottom = Me%ExternalVar%KFloor_Z(i, j)

                    PropertyX%Velocity(i,j,Kbottom) = 0.

                endif
                
            end do
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        
        endif    

        if (MonitorPerformance) call StopWatch ("ModuleFreeVerticalMovement", "BottomBoundary")

        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%KFloor_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'BottomBoundary - ModuleFreeVerticalMovement - ERR02'


    end subroutine BottomBoundary

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine Get_FreeVelocity(FreeVerticalMovementID, PropertyID, Free_Velocity, STAT)

        !Arguments--------------------------------------------------------------
        integer                                 :: FreeVerticalMovementID
        integer                                 :: PropertyID
        real, pointer, dimension(:,:,:)         :: Free_Velocity
        integer, optional, Intent(OUT)          :: STAT
     
        !External--------------------------------------------------------------
        integer                                 :: ready_          

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, STAT_CALL
        type(T_Property), pointer               :: PropertyX

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(FreeVerticalMovementID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mFREEVERTICALMOVEMENT_, Me%InstanceID)

            call Search_Property(PropertyX      = PropertyX,                &
                                 PropertyXID    = PropertyID,               &
                                 STAT           = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                    &
                stop 'Get_FreeVelocity - ModuleFreeVerticalMovement - ERR01'


            Free_Velocity => PropertyX%Velocity

            nullify(PropertyX)


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
            

       !----------------------------------------------------------------------

    end subroutine Get_FreeVelocity


    
    subroutine Get_FreeConvFlux(FreeVerticalMovementID, PropertyID, FreeConvFlux, STAT)

        !Arguments--------------------------------------------------------------
        integer                                 :: FreeVerticalMovementID
        integer                                 :: PropertyID
        real, pointer, dimension(:,:,:)         :: FreeConvFlux
        integer, optional, Intent(OUT)          :: STAT
     
        !External--------------------------------------------------------------
        integer                                 :: ready_          

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, STAT_CALL
        type(T_Property), pointer               :: PropertyX

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(FreeVerticalMovementID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mFREEVERTICALMOVEMENT_, Me%InstanceID)

            call Search_Property(PropertyX      = PropertyX,                &
                                 PropertyXID    = PropertyID,               &
                                 STAT           = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                    &
                stop 'Get_FreeConvFlux - ModuleFreeVerticalMovement - ERR01'


            FreeConvFlux => PropertyX%FreeConvFlux

            nullify(PropertyX)


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
            

       !----------------------------------------------------------------------

    end subroutine Get_FreeConvFlux

    
    !----------------------------------------------------------------------
    
    
    logical function FreeVertPropertyExists(FreeVerticalMovementID, PropertyID)

        !Arguments--------------------------------------------------------------
        integer                                 :: FreeVerticalMovementID
        integer                                 :: PropertyID
     
        !External--------------------------------------------------------------
        integer                                 :: ready_          

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, STAT_CALL
        type(T_Property), pointer               :: PropertyX

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(FreeVerticalMovementID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
                        
            FreeVertPropertyExists = .false.

            call Search_Property(PropertyX      = PropertyX,                &
                                 PropertyXID    = PropertyID,               &
                                 STAT           = STAT_CALL)
            if (STAT_CALL == SUCCESS_)FreeVertPropertyExists = .true.

            nullify(PropertyX)

        else 
            stop 'FreeVertPropertyExists - ModuleFreeVerticalMovement - ERR01'
        end if cd1


    end function FreeVertPropertyExists


    !----------------------------------------------------------------------
    
    
    logical function FreeVertPropertyHasDeposition(FreeVerticalMovementID, PropertyID)

        !Arguments--------------------------------------------------------------
        integer                                 :: FreeVerticalMovementID
        integer                                 :: PropertyID
     
        !External--------------------------------------------------------------
        integer                                 :: ready_          

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, STAT_CALL
        type(T_Property), pointer               :: PropertyX

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(FreeVerticalMovementID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
                        
            FreeVertPropertyHasDeposition = .false.

            call Search_Property(PropertyX      = PropertyX,                &
                                 PropertyXID    = PropertyID,               &
                                 STAT           = STAT_CALL)
            if (STAT_CALL == SUCCESS_)then
                FreeVertPropertyHasDeposition = PropertyX%Deposition
            end if

            nullify(PropertyX)

        else 
            stop 'FreeVertPropertyHasDeposition - ModuleFreeVerticalMovement - ERR01'
        end if cd1


    end function FreeVertPropertyHasDeposition
    
    !--------------------------------------------------------------------------

    subroutine GetFreeVertMovOptions(FreeVerticalMovementID, Salinity, SPM, STAT)

        !Arguments--------------------------------------------------------------
        integer                                 :: FreeVerticalMovementID
        logical, intent(OUT)                    :: Salinity
        logical, intent(OUT)                    :: SPM
        integer, optional, Intent(OUT)          :: STAT
     
        !External--------------------------------------------------------------
        integer                                 :: ready_          

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FreeVerticalMovementID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                               &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            Salinity    = Me%Needs%Salinity
            SPM         = Me%Needs%SPM

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
            

    end subroutine GetFreeVertMovOptions
    
    !--------------------------------------------------------------------------

    
    subroutine UngetFreeVerticalMovement(FreeVerticalMovementID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                 :: FreeVerticalMovementID
        real, dimension(:,:,:), pointer         :: Array
        integer, optional, intent (OUT)         :: STAT

        !External--------------------------------------------------------------
        integer :: ready_   

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FreeVerticalMovementID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            
            nullify(Array)

            call Read_UnLock(mFREEVERTICALMOVEMENT_, Me%InstanceID, "UngetFreeVerticalMovement")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine UngetFreeVerticalMovement

    !--------------------------------------------------------------------------
    
    subroutine Search_Property(PropertyX, PropertyXID, STAT)

        !Arguments-------------------------------------------------------------
        type(T_Property),           pointer             :: PropertyX
        integer         ,           intent (IN)         :: PropertyXID
        integer         , optional, intent (OUT)        :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_ 
        
        !----------------------------------------------------------------------

        STAT_  = UNKNOWN_

        PropertyX => Me%FirstProperty

        do while (associated(PropertyX)) 
            if (PropertyX%ID%IDNumber==PropertyXID) then
                exit        
            else
                PropertyX => PropertyX%Next                 
            end if    
        end do    

       if (associated(PropertyX)) then

            STAT_ = SUCCESS_  

        else
            STAT_  = NOT_FOUND_ERR_  
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine Search_Property

    !----------------------------------------------------------------------
    
    subroutine SetDepositionProbability(FreeVerticalMovementID, DepositionProbability, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: FreeVerticalMovementID
        real,    pointer, dimension(:,:)    :: DepositionProbability
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              
        
        !Local-----------------------------------------------------------------
        integer                             :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(FreeVerticalMovementID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then

            Me%ExternalVar%DepositionProbability => DepositionProbability
            
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SetDepositionProbability


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    Subroutine Kill_FreeVerticalMovement(FreeVerticalMovementID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                 :: FreeVerticalMovementID
        integer, optional, intent(OUT)          :: STAT

        !External--------------------------------------------------------------
        integer                                 :: ready_      
        integer                                 :: STAT_CALL, nUsers       

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_
        type(T_Property), pointer               :: PropertyX

        !griflet
        integer                 :: p
        type(T_VECGW), pointer  :: VECGW

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(FreeVerticalMovementID, ready_)    

cd1:    if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mFREEVERTICALMOVEMENT_,  Me%InstanceID)

            if (nUsers == 0) then

                nUsers = DeassociateInstance(mTIME_,            Me%ObjTime)
                if (nUsers == 0) stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR01'

                nUsers = DeassociateInstance(mHORIZONTALGRID_,  Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR02'

                nUsers = DeassociateInstance(mGEOMETRY_,        Me%ObjGeometry)
                if (nUsers == 0) stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR03'

                nUsers = DeassociateInstance(mMAP_,             Me%ObjMap)
                if (nUsers == 0) stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR04'

#ifdef _USE_PAGELOCKED
                ! FreePageLocked will also nullify the pointers and arrays
                call FreePageLocked(Me%ObjCuda, Me%COEF3%DPtr, Me%COEF3%D)
                call FreePageLocked(Me%ObjCuda, Me%COEF3%EPtr, Me%COEF3%E)
                call FreePageLocked(Me%ObjCuda, Me%COEF3%FPtr, Me%COEF3%F)
#else
                deallocate(Me%COEF3%D,     STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR05'
                nullify(Me%COEF3%D     ) 


                deallocate(Me%COEF3%E,     STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR06'
                nullify(Me%COEF3%E     ) 


                deallocate(Me%COEF3%F,     STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR07'
                nullify(Me%COEF3%F     ) 
#endif _USE_PAGELOCKED

                deallocate(Me%COEF3%D_flux, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR08'
                nullify(Me%COEF3%D_flux) 


                deallocate(Me%COEF3%E_flux, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR09'
                nullify(Me%COEF3%E_flux) 


#ifdef _USE_PAGELOCKED
                ! FreePageLocked will also nullify the pointers and arrays
                call FreePageLocked(Me%ObjCuda, Me%TICOEF3Ptr, Me%TICOEF3)
#else
                deallocate(Me%TICOEF3,      STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR10'
                nullify(Me%TICOEF3     )
#endif _USE_PAGELOCKED

                !Jauch: Missing initialization of PropertyX
                PropertyX => Me%FirstProperty
                do while(associated(PropertyX))

                    deallocate(PropertyX%Velocity,     STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR11'
                    nullify(PropertyX%Velocity)


                    deallocate(PropertyX%FreeConvFlux, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR13'
                    nullify(PropertyX%FreeConvFlux)
                    PropertyX => PropertyX%Next
                end do

                deallocate(Me%VECG, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR14'
                nullify(Me%VECG)


                deallocate(Me%VECW, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'Kill_FreeVerticalMovement - ModuleFreeVerticalMovement - ERR15'
                nullify(Me%VECW)
                
                !griflet
                do p = 1, Me%MaxThreads
                    VECGW => Me%THOMAS%VEC(p)
                    deallocate(VECGW%G)
                    deallocate(VECGW%W)
                enddo 
                deallocate(Me%THOMAS%VEC)
                deallocate(Me%THOMAS%COEF3)
                deallocate(Me%THOMAS)

#ifdef _ENABLE_CUDA                
                !Kills ModuleCuda
                call KillCuda (Me%ObjCuda, STAT = STAT_CALL)
                ! No need to give error yet, Module still has users
#endif _ENABLE_CUDA

                !Deallocates Instance
                call DeallocateInstance 
                
                FreeVerticalMovementID = 0

                STAT_ = SUCCESS_

            end if

        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine Kill_FreeVerticalMovement


    !--------------------------------------------------------------------------

    subroutine DeallocateInstance

        !Local-----------------------------------------------------------------
        type (T_FreeVerticalMovement), pointer           :: AuxFreeVerticalMovement
        type (T_FreeVerticalMovement), pointer           :: PreviousFreeVerticalMovement

        !Updates pointers
        if (Me%InstanceID == FirstObjFreeVerticalMovement%InstanceID) then
            FirstObjFreeVerticalMovement       => FirstObjFreeVerticalMovement%Next
        else
            PreviousFreeVerticalMovement    => FirstObjFreeVerticalMovement
            AuxFreeVerticalMovement         => FirstObjFreeVerticalMovement%Next
            do while (AuxFreeVerticalMovement%InstanceID /= Me%InstanceID)
                PreviousFreeVerticalMovement => AuxFreeVerticalMovement
                AuxFreeVerticalMovement      => AuxFreeVerticalMovement%Next
            enddo

            !Now update linked list
            PreviousFreeVerticalMovement%Next => AuxFreeVerticalMovement%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 
            
    end subroutine DeallocateInstance

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !--------------------------------------------------------------------------

    subroutine Ready (FreeVerticalMovementID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: FreeVerticalMovementID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (FreeVerticalMovementID > 0) then
            call LocateObjFreeVerticalMovement (FreeVerticalMovementID)
            ready_ = VerifyReadLock (mFREEVERTICALMOVEMENT_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjFreeVerticalMovement (FreeVerticalMovementID)

        !Arguments-------------------------------------------------------------
        integer                                     :: FreeVerticalMovementID

        !Local-----------------------------------------------------------------

        Me => FirstObjFreeVerticalMovement
        do while (associated (Me))
            if (Me%InstanceID == FreeVerticalMovementID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                                        &
            stop 'ModuleFreeVerticalMovement - LocateObjFreeVerticalMovement - ERR01'

    end subroutine LocateObjFreeVerticalMovement

end Module ModuleFreeVerticalMovement

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
