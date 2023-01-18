!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid TwoWay nesting
! PROJECT       : Mohid Base 2
! MODULE        : TwoWay
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jul 2018
! REVISION      : Joao Sobrinho - v2.0 - Dec 2018
!>  DESCRIPTION : Module that computes twoway and upscaling operations
!
!> @author
!> Joao Sobrinho
!------------------------------------------------------------------------------

Module ModuleTwoWay

    use ModuleGlobalData
    use ModuleGeometry,         only : GetGeometryVolumes, UnGetGeometry, GetGeometrySize, GetGeometryAreas, &
                                       GetGeometryKFloor, GetGeometryKLink
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, UngetHorizontalGrid, GetHorizontalGridSize, GetConnections, &
                                       UnGetConnections, ConstructP2C_IWD, ConstructP2C_Avrg, GetGridCellArea

    use ModuleHorizontalMap,    only : GetBoundaries, UnGetHorizontalMap
    use ModuleFunctions

    use ModuleMap,              only : GetComputeFaces3D, GetOpenPoints3D, GetWaterPoints3D, UnGetMap
    use ModuleStopWatch,        only : StartWatch, StopWatch

    use ModuleDischarges,        only : GetDischargesNumber, GetDischargeFlowDistribuiton, IsUpscaling, UnGetDischarges
    use ModuleTime
    
    implicit none

    private

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructTwoWay
    private ::  AllocateInstance

    public  :: ConstructTwoWayHydrodynamic
    private ::  Compute_MatrixFilterOB
    public  :: AllocateTwoWayAux
    public  :: ConstructUpscalingDischarges
    public  :: Fill_Upscaling_DecayTime

    !Selector
    public  :: GetUpscalingDischarge
    public  :: GetSonVolInFather
    public  :: UnGetSonVolInFather
    public  :: GetUpscalingIDs
    public  :: GetUpscalingDischargeVelocity
    public  :: ActiveUpscalingMomentumDischarge
    public  :: UnGetUpscalingDischargeVelocity
    !Modifier
    public  :: ModifyTwoWay
    private ::  ComputeAuxMatrixes
    private ::      ComputeSonVolInFather
    private ::          ComputeSonVolInFather_offline
    private ::  Nudging_average
    private ::  Nudging_IWD
    public  :: PrepTwoWay
    private :: Nudging_average_offline
    private :: Flux_Velocity_Offline
    private :: Flux_Velocity_Online
    private :: DischageConc_Online
    public  :: UpscalingVolumeVariation
    public  :: UngetTwoWayExternal_Vars
    public  :: UpscaleDischarge
    public  :: UpscaleDischarge_WP
    public  :: Offline_Upscaling_Discharge
    public  :: Offline_Upscaling_Discharge_WP
    public  :: Offline_Upscaling_Discharge_WP_V2
    public  :: InterpolUpscaling_Velocity

    !Destructor
    public  :: KillTwoWay
    private ::      DeAllocateInstance
    private ::      DeallocateVariables

    !Management
    private ::      Ready
    private ::      ReadyFather
    private ::          LocateObjTwoWay

    !Interfaces----------------------------------------------------------------
    !Parameters----------------------------------------------------------------
    integer, parameter :: DirectionX_      = 1
    integer, parameter :: DirectionY_      = 2
    !Types---------------------------------------------------------------------

    private :: T_Hydro
    type       T_Hydro
        real                                        :: TimeDecay           = null_real
        integer                                     :: IWDn                = null_int
        integer                                     :: InterpolationMethod = null_int
        real                                        :: VelDT               = null_real
        real                                        :: DT                  = null_real
    end type T_Hydro

    private :: T_External
    type       T_External
        integer, dimension(:, :   ), pointer        :: IV               => null()
        integer, dimension(:, :   ), pointer        :: JV               => null()
        integer, dimension(:, :   ), pointer        :: IU               => null()
        integer, dimension(:, :   ), pointer        :: JU               => null()
        integer, dimension(:, :   ), pointer        :: IZ               => null()
        integer, dimension(:, :   ), pointer        :: JZ               => null()
        integer, dimension(:, :, :), pointer        :: KZ               => null()
        real(8),    dimension(:, :, :), pointer     :: VolumeU          => null()
        real(8),    dimension(:, :, :), pointer     :: VolumeV          => null()
        real(8),    dimension(:, :, :), pointer     :: VolumeZ          => null()
        real(8),    dimension(:, :   ), pointer     :: VolumeZ_2D       => null()
        real,    dimension(:, :, :), pointer        :: AreaU            => null()
        real,    dimension(:, :, :), pointer        :: AreaV            => null()
        real,    dimension(:, :   ), pointer        :: AreaZ            => null()
        integer, dimension(:, :, :), pointer        :: Open3D           => null()
        integer, dimension(:, :, :), pointer        :: WaterPoints3D    => null()
        integer, dimension(:, :   ), pointer        :: WaterPoints2D    => null()
        integer, dimension(:, :, :), pointer        :: ComputeFaces3D_U => null()
        integer, dimension(:, :, :), pointer        :: ComputeFaces3D_V => null()
        integer, dimension(:, :   ), pointer        :: BoundaryPoints2D => null()
        integer, dimension(:, :   ), pointer        :: KFloor_U         => null()
        integer, dimension(:, :   ), pointer        :: KFloor_V         => null()
        integer, dimension(:, :   ), pointer        :: KFloor_Z         => null()
        integer, dimension(:, :   ), pointer        :: IWD_Connections_U => null()
        integer, dimension(:, :   ), pointer        :: IWD_Connections_V => null()
        integer, dimension(:, :   ), pointer        :: IWD_Connections_Z => null()
        real, dimension(:      ), pointer           :: IWD_Distances_U   => null()
        real, dimension(:      ), pointer           :: IWD_Distances_V   => null()
        real, dimension(:      ), pointer           :: IWD_Distances_Z   => null()
        integer                                     :: IWD_Nodes_Z       = null_int
        integer                                     :: IWD_Nodes_U       = null_int
        integer                                     :: IWD_Nodes_V       = null_int
        integer, dimension(:, :   ), pointer        :: Connections_Z     => null()
    end type T_External

    private :: T_Discharges
    type       T_Discharges
        integer, dimension(:, :), allocatable       :: U, V, Z
        integer                                     :: n_U = 0, n_V = 0, n_Z = 0
        integer                                     :: Current_U = 0, Current_V = 0, Current_Z = 0
        real(8), dimension(:), allocatable          :: Flow
        real(8), dimension(:), pointer              :: AuxFlow => null()
        integer, dimension(:, :), pointer           :: AuxConnections => null()
        real, dimension(:, :, :), allocatable       :: FlowMatrix
        real, dimension(:, :, :), allocatable       :: VelocityU
        real, dimension(:, :, :), allocatable       :: VelocityV
        real, dimension(:, :, :), allocatable       :: VelocityU_Prev
        real, dimension(:, :, :), allocatable       :: VelocityU_Next
        real, dimension(:, :, :), allocatable       :: VelocityV_Prev
        real, dimension(:, :, :), allocatable       :: VelocityV_Next
        logical                                     :: FirstTimeU = .true.
        logical                                     :: FirstTimeV = .true. 
        logical                                     :: Init = .true.
        logical                                     :: init_U = .true.
        logical                                     :: init_V = .true.
    end type T_Discharges

    private :: T_FatherDomain
    type       T_FatherDomain
        type (T_Size3D)                             :: Size, WorkSize
        type (T_Size2D)                             :: Size2D, WorkSize2D
        type (T_External)                           :: External_Var
        type (T_Discharges)                         :: DischargeCells
        type (T_Discharges)                         :: Discharge
        integer                                     :: InstanceID
        integer                                     :: ObjHorizontalGrid
        integer                                     :: ObjGeometry
        integer                                     :: ObjHorizontalMap
        integer                                     :: ObjMap
        real, dimension (:, :, :), allocatable      :: TotSonIn
        real, dimension (:, :   ), allocatable      :: TotSonIn_2D
        real, dimension (:, :, :), allocatable      :: AuxMatrix
        real, dimension (:, :   ), allocatable      :: AuxMatrix2D
        real, dimension (:, :, :), allocatable      :: IWDNom
        real, dimension (:, :, :), allocatable      :: IWDDenom
    end type T_FatherDomain

    private :: T_TwoWay
    type       T_TwoWay
        integer                                     :: InstanceID
        integer                                     :: DoCycleMethod
        character(PathLength)                       :: ModelName
        real(8), dimension(:, :, :),  pointer       :: Matrix
        integer, dimension(:, :   ),  allocatable   :: IgnoreOBCells

        type (T_External)                           :: External_Var
        type (T_Hydro)                              :: Hydro
        type (T_Discharges)                         :: DischargeCells
        type (T_Discharges)                         :: Discharge
        type (T_Size3D)                             :: Size, WorkSize
        type (T_Size2D)                             :: Size2D, WorkSize2D
        type (T_FatherDomain)                       :: Father
        real, dimension (:, :, :), allocatable      :: TotSonIn
        real, dimension (:, :   ), allocatable      :: TotSonIn_2D
        !real, dimension (:, :, :), pointer          :: TotSonIn => null()
        !real, dimension (:, :   ), pointer          :: TotSonIn_2D => null()
        type (T_TwoWay), pointer                    :: Next

        !Instance of ModuleHorizontalGrid
        integer                                     :: ObjHorizontalGrid = 0
        !Instance of ModuleGeometry
        integer                                     :: ObjGeometry       = 0
        !Instance of ModuleMap
        integer                                     :: ObjMap            = 0
        !Instance of HorizontalMap
        integer                                     :: ObjHorizontalMap  = 0

    end type  T_TwoWay

    !Global Module Variables
    type (T_TwoWay), pointer                         :: FirstObjTwoWay  => null()
    type (T_TwoWay), pointer                         :: Me              => null()

    !--------------------------------------------------------------------------

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Contructs TwoWay global pointer structure.
    !>@param[in] ModelName, TwoWayID, HorizontalGridID, GeometryID, HorizontalMapID, MapID, IntMethod, STAT
    subroutine ConstructTwoWay(ModelName, TwoWayID, HorizontalGridID, GeometryID, HorizontalMapID, MapID, IntMethod, &
    STAT)
        !Arguments---------------------------------------------------------------
        character(Len=*), optional                      :: ModelName
        integer         , intent(OUT)                   :: TwoWayID
        integer         , intent(IN)                    :: HorizontalGridID, GeometryID, HorizontalMapID, MapID
        integer, optional, intent(OUT)                  :: STAT
        integer, optional, intent(IN)                   :: IntMethod
        !External----------------------------------------------------------------
        integer                                         :: ready_
        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, STAT_CALL
        !------------------------------------------------------------------------
        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mTwoWay_)) then
            nullify (FirstObjTwoWay)
            call RegisterModule (mTwoWay_)
        endif

        call Ready(TwoWayID, ready_)

        if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            if (present(ModelName)) Me%ModelName = ModelName
            
            if (Me%ObjHorizontalGrid == 0) &
                Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            if (Me%ObjGeometry == 0) &
                Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            if (Me%ObjMap == 0) &
                Me%ObjMap            = AssociateInstance (mMAP_,            MapID           )
            if (Me%ObjHorizontalMap == 0) &
                Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            !Returns ID
            TwoWayID          = Me%InstanceID

            if (present(IntMethod)) Me%Hydro%InterpolationMethod = IntMethod

            call GetGeometrySize(GeometryID = Me%ObjGeometry, Size = Me%Size, WorkSize = Me%WorkSize, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - ConstructTwoWay - ERR01'

           call GetHorizontalGridSize (HorizontalGridID=Me%ObjHorizontalGrid, Size=Me%Size2D, WorkSize=Me%WorkSize2D, &
                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - ConstructTwoWay - ERR02'

            STAT_ = SUCCESS_
        else
            stop 'ModuleTwoWay - ConstructTwoWay - ERR01'
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ConstructTwoWay

    !--------------------------------------------------------------------------

    subroutine AllocateInstance
        !Local-----------------------------------------------------------------
        type (T_TwoWay), pointer                         :: NewObjTwoWay
        type (T_TwoWay), pointer                         :: PreviousObjTwoWay
        !----------------------------------------------------------------------
        !Allocates new instance
        allocate (NewObjTwoWay)
        nullify  (NewObjTwoWay%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjTwoWay)) then
            FirstObjTwoWay         => NewObjTwoWay
            Me                    => NewObjTwoWay
        else
            PreviousObjTwoWay      => FirstObjTwoWay
            Me                    => FirstObjTwoWay%Next
            do while (associated(Me))
                PreviousObjTwoWay  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjTwoWay
            PreviousObjTwoWay%Next => NewObjTwoWay
        endif

        Me%InstanceID = RegisterNewInstance (mTwoWay_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Contructs hydrodynamic variables of MOHIDWater with TwoWay nesting
    !>@param[in] TwoWayID, TimeDecay, IntMethod, VelDT, DT, NumCellsToIgnore, IWDn, DoCycleMethod, STAT
    subroutine ConstructTwoWayHydrodynamic (TwoWayID, TimeDecay, IntMethod, VelDT, DT, NumCellsToIgnore, IWDn, DoCycleMethod, STAT)
        !Arguments------------------------------------------------------------
        integer                                     :: TwoWayID, & ! ID
                                                       IntMethod
        integer                                     :: NumCellsToIgnore ! number of cells ignored counting from an open boundary
        real                                        :: TimeDecay ! Decay factor in seconds, in the nudging equation
        real                                        :: VelDT, DT
        integer, optional, intent(OUT)              :: STAT
        integer, optional                           :: IWDn
        integer,           intent(IN)               :: DoCycleMethod
        !Local----------------------------------------------------------------
        integer                                     :: STAT_CALL, ready_, STAT_
        !---------------------------------------------------------------------
        call Ready(TwoWayID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Me%Hydro%InterpolationMethod = IntMethod
            !grid interpolation method from child to father grid 1-Volume Averaged; 2-Inverse Weighted Distance
            Me%Hydro%TimeDecay           = TimeDecay
            Me%Hydro%VelDT               = VelDT
            Me%Hydro%DT                  = DT
            Me%DoCycleMethod             = DoCycleMethod

            if (IntMethod == 2) Me%Hydro%IWDn = IWDn

            call GetBoundaries(HorizontalMapID = Me%ObjHorizontalMap,                               &
                               BoundaryPoints2D = Me%External_Var%BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - ConstructTwoWayHydrodynamic - ERR01'

            if (.not. allocated(Me%IgnoreOBCells)) &
                call Compute_MatrixFilterOB (NumCellsToIgnore, Me%External_Var%BoundaryPoints2D)
            
            call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%External_Var%BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - ConstructTwoWayHydrodynamic - ERR02'

            STAT_ = SUCCESS_
        else
            stop 'ModuleTwoWay - ConstructTwoWayHydrodynamic - ERR02'
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ConstructTwoWayHydrodynamic

    !-------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Creates matrix that filters cells close to the openboundary
    !>@param[in] IgnoreOBNumCells, BoundaryPoint
    subroutine Compute_MatrixFilterOB (IgnoreOBNumCells, BoundaryPoint)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                          :: IgnoreOBNumCells
        integer, dimension(:,:), pointer, intent(IN) :: BoundaryPoint
        !Locals ---------------------------------------------------------------
        integer                                      :: ILB, IUB, JLB, JUB, i, j, AuxIgnoreOBNumCells
        !Begin ----------------------------------------------------------------
        ILB = Me%WorkSize2D%ILB; JLB = Me%WorkSize2D%JLB
        IUB = Me%WorkSize2D%IUB; JUB = Me%WorkSize2D%JUB

        !Me%IgnoreOBCells is 0 for when the cell must be ignored.
        allocate (Me%IgnoreOBCells(Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB))
        Me%IgnoreOBCells(:,:) = 1
        Me%IgnoreOBCells(Me%Size2D%ILB,:) = 0
        Me%IgnoreOBCells(Me%Size2D%IUB,:) = 0
        Me%IgnoreOBCells(:,Me%Size2D%JLB) = 0
        Me%IgnoreOBCells(:,Me%Size2D%JUB) = 0

        !This if is here to solve the issue of a east/west channel simulation with only 3 cells
        if ((IUB - IgnoreOBNumCells) <= 3)then
            AuxIgnoreOBNumCells = 1
        else
            AuxIgnoreOBNumCells = IgnoreOBNumCells
        endif
        !compute south border
        do j = JLB, JUB
        do i = ILB, AuxIgnoreOBNumCells
            Me%IgnoreOBCells(i, j) = 1 - BoundaryPoint(ILB, j)
        enddo
        enddo
        AuxIgnoreOBNumCells = IgnoreOBNumCells

        !compute North border
        if ((IUB - IgnoreOBNumCells) <= 3)then
            AuxIgnoreOBNumCells = 0
        else
            AuxIgnoreOBNumCells = IgnoreOBNumCells
        endif

        do j = JLB, JUB
        do i = IUB - AuxIgnoreOBNumCells, IUB
            Me%IgnoreOBCells(i, j) = 1 - BoundaryPoint(IUB, j)
        enddo
        enddo

        AuxIgnoreOBNumCells = IgnoreOBNumCells
        !compute west border
        do j = JLB, AuxIgnoreOBNumCells
        do i = ILB, IUB
            Me%IgnoreOBCells(i, j) = 1 - BoundaryPoint(i, JLB)
        enddo
        enddo
        !compute east border
        if (IgnoreOBNumCells == JUB)then
            AuxIgnoreOBNumCells = AuxIgnoreOBNumCells - 1
        endif
        do j = JUB - AuxIgnoreOBNumCells, JUB
        do i = ILB, IUB
            Me%IgnoreOBCells(i, j) = 1 - BoundaryPoint(i, JUB)
        enddo
        enddo

    end subroutine Compute_MatrixFilterOB
    !-------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Allocates auxiliar matrixes
    !>@param[in] FatherTwoWayID, TwoWayID, CallerID
    subroutine AllocateTwoWayAux(FatherTwoWayID, TwoWayID, CallerID)

        !Arguments-------------------------------------------------------------
        integer         , intent(IN)        :: FatherTwoWayID, TwoWayID, CallerID
        !Local-----------------------------------------------------------------
        integer                             :: ready_, ready_father, ILB, IUB, JLB, JUB, KLB, KUB, STAT_CALL
        logical                             :: isIWD
        type (T_TwoWay), pointer            :: ObjFather
        !----------------------------------------------------------------------
        call Ready (TwoWayID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            isIWD = .false.
            Me%Father%InstanceID = FatherTwoWayID

            call ReadyFather(FatherTwoWayID, ObjFather, ready_father)

            Me%Father%ObjHorizontalGrid = ObjFather%ObjHorizontalGrid
            Me%Father%ObjGeometry       = ObjFather%ObjGeometry
            Me%Father%ObjHorizontalMap  = ObjFather%ObjHorizontalMap
            Me%Father%ObjMap            = ObjFather%ObjMap

            call GetGeometrySize (Me%Father%ObjGeometry, Size = Me%Father%Size, WorkSize = Me%Father%WorkSize, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - AllocateTwoWayAux - ERR10'

            call GetHorizontalGridSize (HorizontalGridID = Me%Father%ObjHorizontalGrid,  Size = Me%Father%Size2D, &
                                       WorkSize = Me%Father%WorkSize2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - AllocateTwoWayAux - ERR20'

            ILB = Me%Father%WorkSize%ILB; JLB = Me%Father%WorkSize%JLB; KUB = Me%Father%WorkSize%KUB
            IUB = Me%Father%WorkSize%IUB; JUB = Me%Father%WorkSize%JUB
            !adjust to number of layers in son domain
            KLB = Me%Father%WorkSize%KLB + (Me%Father%WorkSize%KUB - Me%WorkSize%KUB) !Sobrinho

            if (Me%Hydro%InterpolationMethod == 1)then
                allocate(Me%Father%TotSonIn   (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%Father%AuxMatrix  (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%Father%TotSonIn_2D(ILB:IUB, JLB:JUB))
                allocate(Me%Father%AuxMatrix2D(ILB:IUB, JLB:JUB))

                Me%Father%TotSonIn(:,:,:)  = 0.0
                Me%Father%TotSonIn_2D(:,:) = 0.0
                Me%Father%AuxMatrix(:,:,:) = 0.0
                Me%Father%AuxMatrix2D(:,:) = 0.0

                call ConstructP2C_Avrg(Me%Father%ObjHorizontalGrid, Me%ObjHorizontalGrid)
            elseif (Me%Hydro%InterpolationMethod == 2) then
                call ConstructP2C_IWD(Me%Father%ObjHorizontalGrid, Me%ObjHorizontalGrid)

                allocate(Me%Father%IWDNom  (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%Father%IWDDenom(ILB:IUB, JLB:JUB, KLB:KUB))
                Me%Father%IWDDenom(:,:,:) = 0.0
                Me%Father%IWDNom(:,:,:) = 0.0
            elseif (Me%Hydro%InterpolationMethod == 3) then
                if (CallerID == mFILLMATRIX_) then
                    !Only upscaling Discharge - no nudging
                    allocate(Me%Father%Discharge%VelocityU   (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate(Me%Father%Discharge%VelocityV   (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate(Me%Father%Discharge%VelocityU_Prev  (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate(Me%Father%Discharge%VelocityU_Next  (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate(Me%Father%Discharge%VelocityV_Prev  (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate(Me%Father%Discharge%VelocityV_Next  (ILB:IUB, JLB:JUB, KLB:KUB))
                    Me%Father%Discharge%VelocityU(:,:,:) = 0.0
                    Me%Father%Discharge%VelocityV(:,:,:) = 0.0
                    Me%Father%Discharge%VelocityU_Prev(:,:,:) = 0.0
                    Me%Father%Discharge%VelocityU_Next(:,:,:) = 0.0
                    Me%Father%Discharge%VelocityV_Prev(:,:,:) = 0.0
                    Me%Father%Discharge%VelocityV_Next(:,:,:) = 0.0
                else
                    if (.not. allocated(ObjFather%Discharge%FlowMatrix)) then
                      allocate(ObjFather%Discharge%FlowMatrix   (ILB:IUB, JLB:JUB, KLB:KUB))
                      call SetMatrixValueAllocatable(ObjFather%Discharge%FlowMatrix, ObjFather%WorkSize, 0.0)
                    endif
                endif
                
                call ConstructP2C_Avrg(Me%Father%ObjHorizontalGrid, Me%ObjHorizontalGrid)
            endif
        else
            stop 'ModuleTwoWay - AllocateTwoWayAux - ERR30'
        endif

    end subroutine AllocateTwoWayAux

    ! ------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Searches for next son domain IDs
    !>@param[in] FatherTwoWayID, TwoWayID, HorizontalGridID, HorizontalMapID, LastID, FoundDomain, STAT
    subroutine GetUpscalingIDs (FatherTwoWayID, TwoWayID, HorizontalGridID, HorizontalMapID, LastID, FoundDomain, STAT)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                 :: FatherTwoWayID
        integer, intent(INOUT)              :: LastID
        integer, intent(OUT)                :: TwoWayID, HorizontalGridID, HorizontalMapID, STAT
        logical, intent(OUT)                :: FoundDomain
        !Local-----------------------------------------------------------------
        integer                             :: ready_
        type (T_TwoWay), pointer            :: SonObj
        !----------------------------------------------------------------------
        STAT = UNKNOWN_
        
        call Ready(FatherTwoWayID, ready_)
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            SonObj => FirstObjTwoWay
            
            do while (associated(SonObj))
                
                if (SonObj%InstanceID <= LastID) then
                    !Skip
                elseif (SonObj%Father%InstanceID == FatherTwoWayID) then
                    HorizontalGridID    = SonObj%ObjHorizontalGrid
                    HorizontalMapID     = SonObj%ObjHorizontalMap
                    TwoWayID            = SonObj%InstanceID
                    
                    FoundDomain         = .true.
                    LastID              = SonObj%InstanceID
                    exit
                endif
                
                SonObj => SonObj%Next
            enddo
            
            if (.not. associated(SonObj)) FoundDomain = .false.
            
            STAT = SUCCESS_
        else
            STAT = ready_
        end if
        
    end subroutine GetUpscalingIDs
    ! ------------------------------------------------------------------------
    
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Searches and sets discharge faces of father cell, and finds son cells to use in upscaling
    !>@param[in] SonID, dI, dJ, Connections, SonWaterPoints, FatherWaterPoints, IZ, JZ, Task
    subroutine ConstructUpscalingDischarges(SonID, dI, dJ, Connections, SonWaterPoints, FatherWaterPoints, IZ, JZ, &
    Task, Flag)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                              :: SonID, dI, dJ !dI & dJ = Cell discharge location
        integer, intent(IN)                              :: Task
        integer,  dimension(:,:), pointer, intent(IN)    :: Connections, SonWaterPoints, FatherWaterPoints
        integer, dimension(:,:), pointer , intent(IN)    :: IZ, JZ !Connection between a son Zcell and its father Zcell
        logical, optional                , intent(INOUT) :: Flag
        !Local-----------------------------------------------------------------
        integer                                       :: ready_, STAT_CALL
        !----------------------------------------------------------------------
        call Ready (SonID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call GetGeometryKFloor(GeometryID = Me%Father%ObjGeometry, Z = Me%Father%External_Var%KFloor_Z, &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructUpscalingDischarges - Failed to get Father KfloorZ'

            if (Task == 1) then !find number of lines to allocate
                if (DischargeIsAssociated (Connections, dI, dJ)) then
                    call SearchDischargeFace(Connections, SonWaterPoints, FatherWaterPoints, dI, dJ,  IZ, JZ, &
                                             Me%Father%DischargeCells%n_U, Me%Father%DischargeCells%n_V, &
                                             Me%Father%DischargeCells%n_Z, Me%Father%External_Var%KFloor_Z, &
                                             Me%Father%WorkSize%KUB)
                    !If this call did not generate an error, then a discharge has been found
                    Flag = .true.
                endif

                !Me%Father%DischargeCells%n_U/V - son cells to be included in the calculation
            elseif (Task == 2) then ! allocate

                if (Me%Father%DischargeCells%n_U > 0) then
                    allocate (Me%Father%DischargeCells%Z(Me%Father%DischargeCells%n_Z, 3))
                    allocate (Me%Father%DischargeCells%Flow(Me%Father%DischargeCells%n_Z))
                elseif (Me%Father%DischargeCells%n_V > 0) then
                    allocate (Me%Father%DischargeCells%Z(Me%Father%DischargeCells%n_Z, 3))
                    allocate (Me%Father%DischargeCells%Flow(Me%Father%DischargeCells%n_Z))
                else
                    write(*,*)'No upscaling discharge faces were found between SonID: ', SonID, 'and its father'
                endif

            else ! Fill upscaling matrixes which have the son cells need to be included in the calculation

                !This update is done for each new discharge cell, in order to fill the matrix of discharge cells
                !Save discharge cell IDs into a matrix
                if (Me%Father%DischargeCells%n_Z > 0) then
                    if (DischargeIsAssociated (Connections, dI, dJ)) then
                        call UpdateDischargeConnections(Me%Father%DischargeCells%Current_Z, &
                                                        Me%Father%DischargeCells%Z, Me%Father%External_Var%KFloor_Z, &
                                                        Me%Father%WorkSize%KUB, dI, dJ)
                    endif
                endif

            endif
                call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%KFloor_Z,   STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructUpscalingDischarges : failed to unget KFloor_Z.'
        else
            stop 'Construct_Upscaling_Discharges - ModuleTwoWay -  Failed ready function'
        endif

    end subroutine ConstructUpscalingDischarges


    !-------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Gets upscaling discharge flow
    !>@param[in] TwoWayID, UpscaleDischarge, STAT
    !--------------------------------------------------------------------------
    subroutine GetUpscalingDischarge(TwoWayID, UpscaleDischarge, Done, STAT)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                 :: TwoWayID
        real(8),  dimension(:,:, :), pointer                :: UpscaleDischarge
        integer, optional                                   :: STAT
        logical, optional                   , intent(OUT)   :: Done
        !Local-----------------------------------------------------------------
        type (T_TwoWay), pointer                            :: SonObj
        integer                                             :: ready_, STAT_, line, MaxSize, i, j, k, MaxSize2
        integer                                             :: N_Domains_byDischarge, N_Domains_Nudging
        !----------------------------------------------------------------------
        call Ready(TwoWayID, ready_)
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            N_Domains_byDischarge = 0
            N_Domains_Nudging = 0
            SonObj => FirstObjTwoWay
            do while (associated (SonObj))
                if (SonObj%Father%InstanceID == TwoWayID) then
                    if (SonObj%Hydro%InterpolationMethod == 3) then
                        N_Domains_byDischarge = N_Domains_byDischarge + 1
                    else
                        N_Domains_Nudging = N_Domains_Nudging
                    endif
                endif
                SonObj => SonObj%Next
            enddo
            if (N_Domains_byDischarge > 0 .and. N_Domains_Nudging > 0) then
                write (*,*) 'Model not ready for two different upscaling methods in the same simulation'
            elseif (N_Domains_byDischarge > 0) then
                !Me%Discharge%Flow has the flow from every upscaling domain.
                call SetMatrixValue (UpscaleDischarge, Me%WorkSize, Me%Discharge%FlowMatrix, MaskValue = 0.0)
                if (present(Done)) Done = .true.
            else
                MaxSize = Size(Me%DischargeCells%AuxFlow)
                
                do line = 1, MaxSize
                    i = Me%DischargeCells%AuxConnections(line, 1)
                    j = Me%DischargeCells%AuxConnections(line, 2)
                    k = Me%DischargeCells%AuxConnections(line, 3)
                    UpscaleDischarge(i, j, k) = Me%DischargeCells%AuxFlow(line)
                enddo
            endif
            nullify (SonObj)
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if
        if (present(STAT)) STAT = STAT_

    end subroutine GetUpscalingDischarge

    !---------------------------------------------------------

    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Gets the 2D or 3D integrated son volume inside each father cell
    !>@param[in] TwoWayID, Matrix3D, Matrix2D, STAT
    !--------------------------------------------------------------------------
    subroutine GetSonVolInFather(TwoWayID, Matrix3D, Matrix2D, STAT)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                         :: TwoWayID
        real,  dimension(:,:, :), pointer, optional,  intent(OUT)   :: Matrix3D
        real,  dimension(:,:   ), pointer, optional,  intent(OUT)   :: Matrix2D
        integer, optional                                           :: STAT
        !Local-----------------------------------------------------------------
        integer                                             :: ready_, STAT_
        !----------------------------------------------------------------------

        call Ready(TwoWayID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mTwoWay_, Me%InstanceID)

            if (present(Matrix2D)) Matrix2D => Me%TotSonIn_2D

            if (present(Matrix3D)) Matrix3D => Me%TotSonIn

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetSonVolInFather

    !--------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> UnGets the 2D or 3D integrated son volume inside each father cell
    !>@param[in] TwoWayID, Matrix3D, Matrix2D, STAT
    !--------------------------------------------------------------------------
    subroutine UnGetSonVolInFather(TwoWayID, Matrix3D, Matrix2D, STAT)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                         :: TwoWayID
        real,  dimension(:,:, :), pointer, optional,  intent(OUT)   :: Matrix3D
        real,  dimension(:,:   ), pointer, optional,  intent(OUT)   :: Matrix2D
        integer, optional                                           :: STAT
        !Local-----------------------------------------------------------------
        integer                                             :: ready_, STAT_
        !----------------------------------------------------------------------

        call Ready(TwoWayID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(Matrix2D)) nullify(Matrix2D)
            if (present(Matrix3D)) nullify(Matrix3D)
            if (present(Matrix2D)) Me%TotSonIn_2D(:,:) = 0.0
            if (present(Matrix3D)) Me%TotSonIn(:,:,:)  = 0.0
            call Read_UnLock(mTwoWay_, Me%InstanceID, "UnGetSonVolInFather")

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
    end subroutine UnGetSonVolInFather
     !--------------------------------------------------------------------------
    
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Gets a pointer to the Velocity U or V matrix, to be used in a momemtum discharge
    !>@param[in] TwoWayID, VelocityU, VelocityV
    !--------------------------------------------------------------------------
    subroutine GetUpscalingDischargeVelocity(TwoWayID, VelocityU, VelocityV, STAT)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                         :: TwoWayID
        real,  dimension(:,:, :), pointer, optional,  intent(OUT)   :: VelocityU, VelocityV
        integer, optional                                           :: STAT
        !Local-----------------------------------------------------------------
        integer                                         :: ready_, STAT_, STAT_CALL
        integer                                         :: FatherI_max, FatherJ_max, FatherI_min, FatherJ_min, i, j, k
        integer, dimension(:, :   ), pointer            :: IZ, JZ
        type (T_TwoWay), pointer                        :: SonObj
        !----------------------------------------------------------------------

        call Ready(TwoWayID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mTwoWay_, Me%InstanceID)
            
            if(present(VelocityU)) call SetMatrixValueAllocatable(Me%Discharge%VelocityU, Me%WorkSize, 0.0)
            if(present(VelocityV)) call SetMatrixValueAllocatable(Me%Discharge%VelocityV, Me%WorkSize, 0.0)
            
            SonObj => FirstObjTwoWay
            do while (associated (SonObj))
                !if current object is a son domain then update discharge flow
                if (SonObj%Father%InstanceID == TwoWayID) then
                    call GetHorizontalGrid(HorizontalGridID = SonObj%ObjHorizontalGrid, &
                                           ILinkZ = IZ, JLinkZ = JZ, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetUpscalingDischargeVelocity - Failed to get Son-Father link matrixes'
                    
                    FatherI_max = maxval(IZ)
                    FatherJ_max = maxval(JZ)
                    FatherI_min = minval(IZ, MASK=IZ .GT. 0.0)
                    FatherJ_min = minval(JZ, MASK=JZ .GT. 0.0)
                    if (present(VelocityU)) then
                        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                        do j = FatherJ_min, FatherJ_max
                        do i = FatherI_min, FatherI_max
                            if (SonObj%Father%Discharge%VelocityU(i, j, k) /= 0) then
                                Me%Discharge%VelocityU(i, j, k) = SonObj%Father%Discharge%VelocityU(i, j, k)
                            endif
                        enddo
                        enddo
                        enddo
                    elseif (present(VelocityV)) then
                        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                        do j = FatherJ_min, FatherJ_max
                        do i = FatherI_min, FatherI_max
                            if (SonObj%Father%Discharge%VelocityV(i, j, k) /= 0) then
                                Me%Discharge%VelocityV(i, j, k) = SonObj%Father%Discharge%VelocityV(i, j, k)
                            endif
                        enddo
                        enddo
                        enddo
                    endif
                    call UngetHorizontalGrid(SonObj%ObjHorizontalGrid, IZ,    STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetUpscalingDischargeVelocity-TwoWay-ERR01'
                    call UngetHorizontalGrid(SonObj%ObjHorizontalGrid, JZ,    STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetUpscalingDischargeVelocity-TwoWay-ERR02'
                endif
                
                SonObj => SonObj%Next
                
            enddo
            
            if (present(VelocityU)) VelocityU => Me%Discharge%VelocityU
            if (present(VelocityV)) VelocityV => Me%Discharge%VelocityV
            nullify (SonObj)

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
    end subroutine GetUpscalingDischargeVelocity
    !--------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> checks if upscaling dischage velocity matrix is allocated
    !>@param[in] ID  
    logical function ActiveUpscalingMomentumDischarge (ID, Direction)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: ID, Direction
        !Locals----------------------------------------------------------------
        integer                                         :: ready_
        !Begin-----------------------------------------------------------------
        call Ready(ID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            ActiveUpscalingMomentumDischarge = .false.
            if (direction == DirectionX_) then
                if (allocated(Me%Discharge%VelocityU)) ActiveUpscalingMomentumDischarge = .true.
            else
                if (allocated(Me%Discharge%VelocityV)) ActiveUpscalingMomentumDischarge = .true.
            endif
        else
            write (*,*) 'Cannot acess domain ID', ID
            stop 'ActiveUpscalingMomentumDischarge : objtwoway pointer not ready'
        endif
                
    end function ActiveUpscalingMomentumDischarge
    
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Nullifies pointer to the Velocity U or V matrix
    !>@param[in] TwoWayID, Matrix, STAT
    !--------------------------------------------------------------------------
    subroutine UnGetUpscalingDischargeVelocity(TwoWayID, Matrix, STAT)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: TwoWayID
        real,  dimension(:,:, :), pointer               :: Matrix
        integer, optional                               :: STAT
        !Local-----------------------------------------------------------------
        integer                                         :: ready_, STAT_
        !----------------------------------------------------------------------

        call Ready(TwoWayID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            nullify(Matrix)
            call Read_UnLock(mTwoWay_, Me%InstanceID, "UnGetUpscalingDischargeVelocity")

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
    end subroutine UnGetUpscalingDischargeVelocity
    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Updates father grid domain with son's results
    !>@param[in] SonID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, DischargeFlux, CallerID, VelocityID, TD, STAT
    subroutine ModifyTwoWay(SonID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, CallerID, &
    VelocityID, TD, STAT)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                         :: SonID, CallerID
        integer, optional, intent(IN)                               :: VelocityID
        real,    optional, intent(IN)                               :: TD !TimeDecay for twoway
        integer, optional, intent(OUT)                              :: STAT
        real, dimension(:, :, :), pointer, optional, intent(INOUT)  :: FatherMatrix
        real, dimension(:, :, :), pointer, optional, intent(IN)     :: SonMatrix
        real, dimension(:, :),    pointer, optional, intent(INOUT)  :: FatherMatrix2D
        real, dimension(:, :),    pointer, optional, intent(IN)     :: SonMatrix2D
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, InterpolMethod
        real                                        :: TimeDecay
        logical                                     :: Offline, Online
        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(SonID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Offline = .false.

            if (MonitorPerformance) call StartWatch ("ModuleTwoWay", "ModifyTwoWay")

            if (CallerID == mHydrodynamic_) then
                TimeDecay = Me%Hydro%TimeDecay
                InterpolMethod = Me%Hydro%InterpolationMethod
                Online = .true.
            elseif (CallerID == mWATERPROPERTIES_) then
                ! for now interpolation method is set by the hydrodynamic module. The if is here for when
                ! each property can set its own interpolation method
                InterpolMethod = Me%Hydro%InterpolationMethod
                TimeDecay = TD
                Online = .true.
            elseif (CallerID == mField4D_) then
                InterpolMethod = Me%Hydro%InterpolationMethod !Was Defined in constructTwoWay, via fillmatrix.
                Offline =  .true.
            endif

            !if it is a 3D matrix
            if (present(FatherMatrix) .and. Me%Hydro%InterpolationMethod /= 3) then
                if (present(VelocityID) .and. Online)then
                    if (VelocityID == VelocityU_) then
                        !if 3D matrixes were sent. Even 2D domains allocate a 3D matrix (only one vertical layer)
                        !Type_U
                        call ComputeAuxMatrixes (Volume_3D = Me%External_Var%VolumeU, InterpolMethod = InterpolMethod,&
                                                 Ilink = Me%External_Var%IU, Jlink = Me%External_Var%JU, &
                                                 Klink = Me%External_Var%KZ, VelocityID = VelocityID,    &
                                                 Offline = Offline)
                    else
                        !Type_V
                        call ComputeAuxMatrixes (Volume_3D = Me%External_Var%VolumeV, InterpolMethod = InterpolMethod,&
                                                 Ilink = Me%External_Var%IV, Jlink = Me%External_Var%JV,              &
                                                 Klink = Me%External_Var%KZ, VelocityID = VelocityID,                 &
                                                 Offline = Offline)
                    endif
                else
                    !Type Z - specific for offline upscaling
                    call ComputeAuxMatrixes     (Volume_3D = Me%External_Var%VolumeZ, InterpolMethod = InterpolMethod,&
                                                 Ilink = Me%External_Var%IZ, Jlink = Me%External_Var%JZ, &
                                                 Klink = Me%External_Var%KZ, Offline = Offline)
                endif
            elseif (Me%Hydro%InterpolationMethod /= 3) then
                !if a 2D matrix was sent (specific for waterLevel - at least for MohidWater).
                call ComputeAuxMatrixes (Volume_2D = Me%External_Var%VolumeZ_2D, InterpolMethod = InterpolMethod, &
                                         Ilink = Me%External_Var%IZ, Jlink = Me%External_Var%JZ, Offline = Offline)
            endif

            if (CallerID == mField4D_) then
                if (InterpolMethod == 1) then
                    call Nudging_average_offline (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D)
                elseif (InterpolMethod == 3) then
                    if (present(VelocityID)) then
                        call Flux_Velocity_Offline (VelocityID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D)
                    else
                        call DischageConc_Offline (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D)
                    endif
                endif
            elseif (CallerID == mHydrodynamic_) then
                if (InterpolMethod == 1) then
                    call Nudging_average (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, TimeDecay)
                elseif (InterpolMethod == 2) then
                    call Nudging_IWD (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, TimeDecay)
                elseif (InterpolMethod == 3) then
                    call Flux_Velocity_Online (VelocityID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D)
                endif
            else
                if (InterpolMethod == 1) then
                    call Nudging_average (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, TimeDecay)
                elseif (InterpolMethod == 2) then
                    call Nudging_IWD (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, TimeDecay)
                elseif (InterpolMethod == 3) then
                    call DischageConc_Online (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D)
                endif
            endif
            

            if (MonitorPerformance) call StopWatch ("ModuleTwoWay", "ModifyTwoWay")
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyTwoWay

    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Gets external variables
    !>@param[in] SonID, CallerID, STAT
    subroutine PrepTwoWay (SonID, CallerID, STAT)
        !Arguments--------------------------------------------------------------
        integer,           intent(IN)               :: SonID 
        integer, optional, intent(IN)               :: CallerID
        integer, optional, intent(OUT)              :: STAT
        !Locals-----------------------------------------------------------------
        integer                                     :: STAT_CALL, ready_, STAT_, CallerID_
        !Begin------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(SonID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            CallerID_ = 0

            if (present(callerID)) CallerID_ = callerID
            if (CallerID_ == 1) then
                !For future developments (when other modules call for twoway)
            endif

            call GetHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid, &
                                   ILinkV = Me%External_Var%IV, JLinkV = Me%External_Var%JV, &
                                   ILinkU = Me%External_Var%IU, JLinkU = Me%External_Var%JU, &
                                   ILinkZ = Me%External_Var%IZ, JLinkZ = Me%External_Var%JZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Failed to get Son-Father link matrixes'
            
            call GetGeometryKLink(GeometryID = Me%ObjGeometry, KLinkZ = Me%External_Var%KZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Failed to get Son-Father Klink matrix'

            Call GetGeometryAreas(GeometryID = Me%ObjGeometry, AreaU = Me%External_Var%AreaU, &
                                  AreaV = Me%External_Var%AreaV, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Could not get Son AreaU or AreaV matrix'

            call GetGeometryVolumes(GeometryID = Me%ObjGeometry, VolumeU = Me%External_Var%VolumeU, &
                                    VolumeV = Me%External_Var%VolumeV, VolumeZ = Me%External_Var%VolumeZ, &
                                    VolumeZ_2D = Me%External_Var%VolumeZ_2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Failed to get son U/V volume matrixes'

            call GetOpenPoints3D   (Map_ID = Me%ObjMap, OpenPoints3D = Me%External_Var%Open3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Failed to get Son OpenPoints3D'

            call GetComputeFaces3D (Map_ID = Me%ObjMap, ComputeFacesU3D = Me%External_Var%ComputeFaces3D_U, &
                                    ComputeFacesV3D = Me%External_Var%ComputeFaces3D_V, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Failed to get Son U/V ComputeFaces3D'

            call GetBoundaries(HorizontalMapID = Me%ObjHorizontalMap,                               &
                               BoundaryPoints2D = Me%External_Var%BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Failed to get Son BoundaryPoints2D'

            call GetWaterPoints3D(Map_ID = Me%ObjMap, WaterPoints3D = Me%External_Var%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'PrepTwoWay - Failed to get Son WaterPoints3D'

            if (Me%Hydro%InterpolationMethod == 2) then
                call GetConnections (HorizontalGridID = Me%ObjHorizontalGrid, &
                                     Connections_U = Me%External_Var%IWD_Connections_U, &
                                     IWD_Distances_U = Me%External_Var%IWD_Distances_U, &
                                     Connections_V = Me%External_Var%IWD_Connections_V, &
                                     IWD_Distances_V = Me%External_Var%IWD_Distances_V, &
                                     Connections_Z = Me%External_Var%IWD_Connections_Z, &
                                     IWD_Distances_Z = Me%External_Var%IWD_Distances_Z, &
                                     IWD_Nodes_Z  = Me%External_Var%IWD_Nodes_Z,        &
                                     IWD_Nodes_U = Me%External_Var%IWD_Nodes_U,         &
                                     IWD_Nodes_V = Me%External_Var%IWD_Nodes_V, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Failed to get IWD connections'
            endif

            call GetGeometryAreas(GeometryID = Me%Father%ObjGeometry, AreaU = Me%Father%External_Var%AreaU, &
                                  AreaV = Me%Father%External_Var%AreaV, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Failed to get father AreaU or AreaV matrix'

            call GetGridCellArea(HorizontalGridID = Me%Father%ObjHorizontalGrid, &
                                 GridCellArea = Me%Father%External_Var%AreaZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Failed to get Father AreaZ matrix'

            call GetOpenPoints3D (Map_ID = Me%Father%ObjMap, OpenPoints3D = Me%Father%External_Var%Open3D, &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Failed to get Father OpenPoints3D'

            call GetGeometryVolumes(GeometryID = Me%Father%ObjGeometry, VolumeU = Me%Father%External_Var%VolumeU, &
                                    VolumeV = Me%Father%External_Var%VolumeV, VolumeZ = Me%Father%External_Var%VolumeZ,&
                                    VolumeZ_2D = Me%Father%External_Var%VolumeZ_2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Failed to get Father Volume U/V/2D matrixes'

            call GetComputeFaces3D(Map_ID = Me%Father%ObjMap, ComputeFacesU3D = Me%Father%External_Var%ComputeFaces3D_U, &
                                   ComputeFacesV3D = Me%Father%External_Var%ComputeFaces3D_V, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Failed to get Father ComputeFaces3D U/V'

            call GetGeometryKFloor(GeometryID = Me%Father%ObjGeometry, U = Me%Father%External_Var%KFloor_U, &
                                   V = Me%Father%External_Var%KFloor_V,                                &
                                   Z = Me%Father%External_Var%KFloor_Z, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PrepTwoWay - Failed to get Father KfloorU/V'

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        endif

        if (present(STAT)) STAT = STAT_

    end subroutine PrepTwoWay

    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Computes water volume variation from velocities at the Z cell faces (applied to the father domain).
    !>@param[in] MatrixNew, MatrixOld, DT, VelocityID
    subroutine UpscalingVolumeVariation (MatrixNew, MatrixOld, DT, VelocityID)
    !Arguments-------------------------------------------------------------
    real, dimension(:, :, :), allocatable, intent(INOUT) :: MatrixNew
    real, dimension(:, :, :), pointer,     intent(IN)    :: MatrixOld
    real,                                  intent(IN)    :: DT
    integer,                               intent(IN)    :: VelocityID
    !Locals----------------------------------------------------------------
    integer                                          :: i, j, k, ILB, IUB, JLB, JUB, KLB, KUB, Kbottom, CHUNK
    real                                             :: AuxWest, AuxEast, AuxSouth, AuxNorth, AuxBottom, AuxUp
    !Begin-----------------------------------------------------------------
    ILB = Me%Father%WorkSize%ILB ; JLB = Me%Father%WorkSize%JLB ; KLB = Me%Father%WorkSize%KLB
    IUB = Me%Father%WorkSize%IUB ; JUB = Me%Father%WorkSize%JUB ; KUB = Me%Father%WorkSize%KUB
    CHUNK = CHUNK_K(KLB, KUB)
    if (VelocityID == VelocityU_) then
        !$OMP PARALLEL PRIVATE(i,j,k,AuxWest, AuxEast)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%Father%External_Var%ComputeFaces3D_U(i,j,k) == 1)then
                AuxWest = (MatrixNew(i, j  , k) - MatrixOld(i, j  , k)) * Me%Father%External_Var%AreaU(i, j  , k)
                AuxEast = (MatrixNew(i, j+1, k) - MatrixOld(i, j+1, k)) * Me%Father%External_Var%AreaU(i, j+1, k)
                !m3                = m3/s                * s
                MatrixNew(i, j, k) = (AuxWest + AuxEast) * DT
            endif
        enddo
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    elseif (VelocityID == VelocityV_) then
        !$OMP PARALLEL PRIVATE(i,j,k,AuxSouth, AuxNorth)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%Father%External_Var%ComputeFaces3D_V(i,j,k) == 1)then
                AuxSouth = (MatrixNew(i  , j, k) - MatrixOld(i  , j, k)) * Me%Father%External_Var%AreaV(i  , j, k)
                AuxNorth = (MatrixNew(i+1, j, k) - MatrixOld(i+1, j, k)) * Me%Father%External_Var%AreaV(i+1, j, k)
                !m3                = m3                 +         m3/s          * s
                MatrixNew(i, j, k) = MatrixNew(i, j, k) + (AuxSouth + AuxNorth) * DT
            endif
        enddo
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    elseif (VelocityID == VelocityW_) then
        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j,k,AuxUp, AuxBottom, Kbottom)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%Father%External_Var%Open3D(i,j,KUB) == 1)then
                Kbottom = Me%Father%External_Var%KFloor_Z(i, j)
                AuxBottom = (MatrixNew(i, j, Kbottom) - MatrixOld(i, j, Kbottom)) * Me%Father%External_Var%AreaZ(i, j)
                do k = Kbottom+1, KUB
                    AuxUp     = (MatrixNew(i, j, k) - MatrixOld(i, j, k)) * Me%Father%External_Var%AreaZ(i, j)
                    !m3                = m3                 +        m3/s         * s
                    MatrixNew(i, j, k) = MatrixNew(i, j, k) + (AuxBottom + AuxUp) * DT
                    AuxBottom = AuxUp
                enddo
            endif
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    endif

    end subroutine UpscalingVolumeVariation

    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Computes auxiliar matrixes for the feedback
    !>@param[in] Volume_3D, Volume_2D, VelocityID, InterpolMethod, Ilink, Jlink, Offline
    subroutine ComputeAuxMatrixes(Volume_3D, Volume_2D, VelocityID, InterpolMethod, Ilink, Jlink, Klink, Offline)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                        :: interpolMethod
        real(8), dimension(:, :, :), pointer, optional, intent(IN) :: Volume_3D
        real(8), dimension(:, :),    pointer, optional, intent(IN) :: Volume_2D
        integer, dimension(:, :),    pointer, intent(IN)           :: Ilink, Jlink
        integer, dimension(:, :, :), pointer, optional, intent(IN) :: Klink
        integer, optional,                              intent(IN) :: VelocityID
        logical,                                        intent(IN) :: Offline
        !Local-----------------------------------------------------------------
        !Begin -----------------------------------------------------------------
        if (present(Volume_3D)) then
            !Goes for 3D
            if  (interpolMethod == 1)then
                Me%Father%TotSonIn (:,:,:) = 0.0
                Me%Father%AuxMatrix(:,:,:) = 0.0
                ! Volume Weighted average
                if (present(VelocityID))then
                    if (VelocityID == VelocityU_)then
                        call ComputeSonVolInFather(Volume_3D = Volume_3D, Ilink = Ilink, Jlink = Jlink, Klink = Klink,&
                                                SonComputeFaces = Me%External_Var%ComputeFaces3D_U, Offline=Offline, &
                                                VelocityID = VelocityID)
                    else
                        call ComputeSonVolInFather(Volume_3D = Volume_3D, Ilink = Ilink, Jlink = Jlink, Klink = Klink,&
                                                SonComputeFaces = Me%External_Var%ComputeFaces3D_V, Offline=Offline, &
                                                VelocityID = VelocityID)
                    endif
                else
                    call ComputeSonVolInFather   (Volume_3D = Volume_3D, Ilink = Ilink, Jlink = Jlink, Klink = Klink, &
                                                  Offline = Offline)
                endif
            else
                Me%Father%IWDNom = 0.0
                Me%Father%IWDDenom = 0.0
            endif
        else
        !Goes for 2D
            if  (interpolMethod == 1)then

                Me%Father%AuxMatrix2D(:,:) = 0.0
                Me%Father%TotSonIn_2D(:,:) = 0.0
                ! Volume Weighted average
                call ComputeSonVolInFather   (Volume_2D = Volume_2D, Ilink = Ilink, Jlink = Jlink, Offline = Offline)
            else
                Me%Father%IWDNom = 0.0
                Me%Father%IWDDenom = 0.0
            endif
        endif

    end subroutine ComputeAuxMatrixes

    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Computes Son cells volume inside each father cell
    !>@param[in] Volume_3D, Volume_2D, Ilink, Jlink, SonComputeFaces
    subroutine ComputeSonVolInFather (Volume_3D, Volume_2D, Ilink, Jlink, Klink, SonComputeFaces, Offline, VelocityID)
        !Arguments--------------------------------------------------------------------------------
        real(8), dimension(:, :, :), pointer, optional :: Volume_3D
        real(8), dimension(:, :),    pointer, optional :: Volume_2D
        integer, dimension(:, :), pointer              :: Ilink, Jlink
        integer, dimension(:, :, :), pointer, optional :: SonComputeFaces, Klink
        logical,                            intent(IN) :: Offline
        integer, optional                              :: VelocityID
        !Local variables--------------------------------------------------------------------------
        integer                                 :: i, j, k, ifather, jfather, kfather, CHUNK, Flag, Flag2
        !Begin------------------------------------------------------------------------------------

        if (Offline) then
            if (present(Klink)) then
                call ComputeSonVolInFather_offline(Volume_3D, Volume_2D, Ilink, Jlink, Klink = Klink)
            else
                stop 'ComputeSonVolInFather - ModuleTwoWay - Klink matrix not present'
            endif
        endif

        if (present(Volume_3D) .and. .not. offline) then
            
            CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
            !!$OMP PARALLEL PRIVATE(i,j,k,Flag, ifather, jfather, kfather)
            if (present(SonComputeFaces))then
                if (VelocityID == VelocityU_) then
                !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                    do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        Flag = Me%External_Var%Open3D(i, j, k) + SonComputeFaces(i, j, k) + Me%IgnoreOBCells(i, j)
                        if ((Flag == 3) .and. Klink(i, j, k) /= FillValueInt) then
                            !if ((Me%IgnoreOBCells(i, j) == 0 .and. Me%IgnoreOBCells(i, j + 1) == 1) &
                            !      .or. (Me%IgnoreOBCells(i, j) == 1 .and. Me%IgnoreOBCells(i, j + 1) == 1)) then
                            ifather = ILink(i, j); jfather = JLink(i, j); kfather = Klink(i, j, k)
                            Me%Father%TotSonIn(ifather, jfather, kfather) = Me%Father%TotSonIn(ifather, jfather, kfather) &
                                                                            + Volume_3D(i, j, k)
                            !endif
                        endif
                    enddo
                    enddo
                    enddo
                
                else
                    do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        Flag = Me%External_Var%Open3D(i, j, k) + Me%IgnoreOBCells(i, j) + SonComputeFaces(i, j, k)
                        if ((Flag == 3) .and. Klink(i, j, k) /= FillValueInt) then
                            ifather = ILink(i, j); jfather = JLink(i, j); kfather = Klink(i, j, k)
                            Me%Father%TotSonIn(ifather, jfather, kfather) = Me%Father%TotSonIn(ifather, jfather, kfather) &
                                                                            + Volume_3D(i, j, k)
                        endif
                    enddo
                    enddo
                    enddo
                endif
                
                
                !!$OMP END DO
            else
                !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    Flag = Me%External_Var%Open3D(i, j, k) + Me%IgnoreOBCells(i, j)
                    if (Flag == 2 .and. Klink(i, j, k) /= FillValueInt) then
                        ifather = ILink(i, j); jfather = JLink(i, j); kfather = Klink(i, j, k)
                        Me%Father%TotSonIn(ifather, jfather, kfather) = Me%Father%TotSonIn(ifather, jfather, kfather) &
                                                                        + Volume_3D(i, j, k)
                    endif
                enddo
                enddo
                enddo
                !!$OMP END DO
            endif
            !!$OMP END PARALLEL
        elseif (present(Volume_2D) .and. .not. offline) then
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                Flag = Me%External_Var%Open3D(i, j, Me%WorkSize%KUB) + Me%IgnoreOBCells(i, j)
                if (Flag == 2) then
                    Me%Father%TotSonIn_2D(ILink(i, j), JLink(i, j)) = Me%Father%TotSonIn_2D(ILink(i, j), JLink(i, j)) &
                                                                        + Volume_2D(i, j)
                endif
            enddo
            enddo
        endif

    end subroutine ComputeSonVolInFather

    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Computes Son cells volume inside each father cell for offline simulations
    !>@param[in] Volume_3D, Volume_2D, Ilink, Jlink, SonComputeFaces
    subroutine ComputeSonVolInFather_offline (Volume_3D, Volume_2D, Ilink, Jlink, Klink)
        !Arguments--------------------------------------------------------------------------------
        real(8), dimension(:, :, :), pointer, optional :: Volume_3D
        real(8), dimension(:, :),    pointer, optional :: Volume_2D
        integer, dimension(:, :),    pointer           :: Ilink, Jlink
        integer, dimension(:, :, :), pointer, optional :: Klink
        !Local variables--------------------------------------------------------------------------
        integer                                        :: i, j, k, ifather, jfather, kfather, KLB, KUB
        integer                                        :: CHUNK
        !Begin------------------------------------------------------------------------------------
        KLB = Me%WorkSize%KLB; KUB = Me%WorkSize%KUB
        if (present(Volume_3D)) then
            CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)

            !!$OMP PARALLEL PRIVATE(i,j,k, ifather, jfather, kfather)
            !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Klink(i, j, k) /= FillValueInt) then
                    ifather = ILink(i, j); jfather = JLink(i, j); kfather = Klink(i, j, k)
                    Me%Father%TotSonIn(ifather,jfather,kfather) = Me%Father%TotSonIn(ifather, jfather, kfather) &
                                                                + Volume_3D(i, j, k)                            &
                                                                * Me%External_Var%Open3D(i, j, k)
                endif
            enddo
            enddo
            enddo
            !!$OMP END DO
            !!$OMP END PARALLEL

        elseif (present(Volume_2D)) then
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                ifather = ILink(i, j); jfather = JLink(i, j)
                Me%Father%TotSonIn_2D(ifather, jfather) = Me%Father%TotSonIn_2D(ifather, jfather) &
                                                        + Volume_2D(i, j) * Me%External_Var%Open3D(i, j, KUB)
            enddo
            enddo
        endif

    end subroutine ComputeSonVolInFather_offline

    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Calls Feedback routines present in mondule functions, based on the type of input matrixes
    !>@param[in] FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, LocalTimeDecay
    subroutine Nudging_average (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, TimeDecay)
        !Arguments-------------------------------------------------------------
        integer, optional                                          :: VelocityID
        real, dimension(:, :, :), pointer, optional, intent(INOUT) :: FatherMatrix
        real, dimension(:, :, :), pointer, optional, intent(IN)    :: SonMatrix
        real, dimension(:, :),    pointer, optional, intent(INOUT) :: FatherMatrix2D
        real, dimension(:, :),    pointer, optional, intent(IN)    :: SonMatrix2D
        real                                                       :: TimeDecay
        !Local-----------------------------------------------------------------
        if (present(FatherMatrix)) then

            if (present(VelocityID))then
                !If velocity U
                if (VelocityID == VelocityU_) then
                    call FeedBack_Avrg_UV (FatherMatrix         = FatherMatrix,                           &
                                           SonMatrix            = SonMatrix,                              &
                                           Open3DFather         = Me%Father%External_Var%Open3D,           &
                                           Open3DSon            = Me%External_Var%Open3D,                  &
                                           FatherComputeFaces3D = Me%Father%External_Var%ComputeFaces3D_U, &
                                           SonComputeFaces3D    = Me%External_Var%ComputeFaces3D_U,        &
                                           SizeFather           = Me%Father%WorkSize,                     &
                                           SizeSon              = Me%WorkSize,                            &
                                           Ilink                = Me%External_Var%IU,                      &
                                           Jlink                = Me%External_Var%JU,                      &
                                           Klink                = Me%External_Var%KZ,                       &
                                           DecayTime            = TimeDecay,                         &
                                           DT                   = Me%Hydro%VelDT,                         &
                                           SonVolInFather       = GetPointer(Me%Father%TotSonIn),                  &
                                           AuxMatrix            = GetPointer(Me%Father%AuxMatrix),                    &
                                           VolumeSon            = Me%External_Var%VolumeU,                 &
                                           VolumeFather         = Me%Father%External_Var%VolumeU,         &
                                           IgnoreOBCells        = GetPointer(Me%IgnoreOBCells))
                !If velocity V
                else
                    call FeedBack_Avrg_UV (FatherMatrix         = FatherMatrix,                           &
                                           SonMatrix            = SonMatrix,                              &
                                           Open3DFather         = Me%Father%External_Var%Open3D,           &
                                           Open3DSon            = Me%External_Var%Open3D,                  &
                                           FatherComputeFaces3D = Me%Father%External_Var%ComputeFaces3D_V, &
                                           SonComputeFaces3D    = Me%External_Var%ComputeFaces3D_V,        &
                                           SizeFather           = Me%Father%WorkSize,                     &
                                           SizeSon              = Me%WorkSize,                            &
                                           Ilink                = Me%External_Var%IV,                      &
                                           Jlink                = Me%External_Var%JV,                      &
                                           Klink                = Me%External_Var%KZ,                       &
                                           DecayTime            = TimeDecay,                         &
                                           DT                   = Me%Hydro%VelDT,                         &
                                           SonVolInFather       = GetPointer(Me%Father%TotSonIn),                  &
                                           AuxMatrix            = GetPointer(Me%Father%AuxMatrix),                    &
                                           VolumeSon            = Me%External_Var%VolumeV,                 &
                                           VolumeFather         = Me%Father%External_Var%VolumeV,          &
                                           IgnoreOBCells        = GetPointer(Me%IgnoreOBCells))
                endif
            else
                !compute nudging Z type cell
                call FeedBack_Avrg   (FatherMatrix     = FatherMatrix,                    &
                                      SonMatrix        = SonMatrix,                       &
                                      Open3DFather     = Me%Father%External_Var%Open3D,    &
                                      Open3DSon        = Me%External_Var%Open3D,           &
                                      SizeFather       = Me%Father%WorkSize,              &
                                      SizeSon          = Me%WorkSize,                     &
                                      Ilink            = Me%External_Var%IZ,               &
                                      Jlink            = Me%External_Var%JZ,               &
                                      Klink            = Me%External_Var%KZ,                       &
                                      DecayTime        = TimeDecay,                  &
                                      DT               = Me%Hydro%DT,                     &
                                      SonVolInFather   = GetPointer(Me%Father%TotSonIn),           &
                                      AuxMatrix        = GetPointer(Me%Father%AuxMatrix),             &
                                      VolumeSon        = Me%External_Var%VolumeZ,          &
                                      VolumeFather     = Me%Father%External_Var%VolumeZ, &
                                      IgnoreOBCells    = GetPointer(Me%IgnoreOBCells))
            endif

        else

            call FeedBack_Avrg_WL (FatherMatrix2D   = FatherMatrix2D,                  &
                                  SonMatrix2D      = SonMatrix2D,                     &
                                  Open3DFather     = Me%Father%External_Var%Open3D,    &
                                  Open3DSon        = Me%External_Var%Open3D,           &
                                  SizeFather       = Me%Father%WorkSize,              &
                                  SizeSon          = Me%WorkSize,                     &
                                  Ilink            = Me%External_Var%IZ,               &
                                  Jlink            = Me%External_Var%JZ,               &
                                  DecayTime        = TimeDecay,                  &
                                  DT               = Me%Hydro%DT,                     &
                                  SonVolInFather2D = GetPointer(Me%Father%TotSonIn_2D),        &
                                  AuxMatrix2D      = GetPointer(Me%Father%AuxMatrix2D),           &
                                  VolumeSon2D      = Me%External_Var%VolumeZ_2D,       &
                                  VolumeFather2D   = Me%Father%External_Var%VolumeZ_2D, &
                                  IgnoreOBCells    = GetPointer(Me%IgnoreOBCells))

        endif

    end subroutine Nudging_average

    !------------------------------------------------------------------------------

    !>@author Joao Sobrinho +Atlantic
    !>@Brief
    !>Calls offline Feedback routines present in mondule functions, based on the type of input matrixes
    !>Matrix output is the product of concentration and volume.
    !>@param[in] FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, LocalTimeDecay
    subroutine Nudging_average_offline(FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D)
        !Arguments-------------------------------------------------------------
        real, dimension(:, :, :), pointer, optional, intent(INOUT)  :: FatherMatrix
        real, dimension(:, :, :), pointer, optional, intent(IN)     :: SonMatrix
        real, dimension(:, :),    pointer, optional, intent(INOUT)  :: FatherMatrix2D
        real, dimension(:, :),    pointer, optional, intent(IN)     :: SonMatrix2D
        !Local-----------------------------------------------------------------
        type (T_TwoWay), pointer                                    :: ObjFather
        integer                                                     :: ready_father
        !Begin ----------------------------------------------------------------
        if (present(FatherMatrix)) then
            !3D
            !compute nudging Z type cell
            call Upscaling_Avrg (   FatherMatrix     = FatherMatrix,                        &
                                    SonMatrix        = SonMatrix,                           &
                                    FatherMask       = Me%Father%External_Var%Open3D,       &
                                    SonMask          = Me%External_Var%Open3D,              &
                                    SizeFather       = Me%Father%WorkSize,                  &
                                    SizeSon          = Me%WorkSize,                         &
                                    Ilink            = Me%External_Var%IZ,                  &
                                    Jlink            = Me%External_Var%JZ,                  &
                                    Klink            = Me%External_Var%KZ,                  &
                                    VolumeSon        = Me%External_Var%VolumeZ,             &
                                    TotSonIn         = GetPointer(Me%Father%TotSonIn))

            call ReadyFather(Me%Father%InstanceID, ObjFather, ready_father)
            !if (.not. associated(ObjFather%TotSonIn)) then
            !    ObjFather%TotSonIn => Me%Father%TotSonIn
            !else
            !    where (Me%Father%TotSonIn(:,:,:) /= 0) ObjFather%TotSonIn(:,:,:) = Me%Father%TotSonIn
            !endif
            if (.not. allocated(ObjFather%TotSonIn)) then
                allocate(ObjFather%TotSonIn(ObjFather%WorkSize%ILB:ObjFather%WorkSize%IUB, &
                                            ObjFather%WorkSize%JLB:ObjFather%WorkSize%JUB, &
                                            ObjFather%WorkSize%KLB:ObjFather%WorkSize%KUB))
            endif
            where (Me%Father%TotSonIn(:,:,:) /= 0) ObjFather%TotSonIn(:,:,:) = Me%Father%TotSonIn(:,:,:)

        else
            call Upscaling_Avrg_WL (    FatherMatrix2D   = FatherMatrix2D,                      &
                                        SonMatrix2D      = SonMatrix2D,                         &
                                        FatherMask       = Me%Father%External_Var%Open3D,       &
                                        SonMask          = Me%External_Var%Open3D,              &
                                        SizeFather       = Me%Father%WorkSize,                  &
                                        SizeSon          = Me%WorkSize,                         &
                                        Ilink            = Me%External_Var%IZ,                  &
                                        Jlink            = Me%External_Var%JZ,                  &
                                        VolumeSon2D      = Me%External_Var%VolumeZ_2D,          &
                                        TotSonIn2D       = GetPointer(Me%Father%TotSonIn_2D))

            call ReadyFather(Me%Father%InstanceID, ObjFather, ready_father)
            !if (.not. associated(ObjFather%TotSonIn_2D)) then
            !    ObjFather%TotSonIn_2D => Me%Father%TotSonIn_2D
            !else
            !     where (ObjFather%TotSonIn_2D(:,:) /= 0) ObjFather%TotSonIn_2D(:,:) = Me%Father%TotSonIn_2D
            !endif
            if (.not. allocated(ObjFather%TotSonIn_2D)) then
                allocate(ObjFather%TotSonIn_2D( ObjFather%WorkSize%ILB:ObjFather%WorkSize%IUB, &
                                                ObjFather%WorkSize%JLB:ObjFather%WorkSize%JUB))
            endif
            where (ObjFather%TotSonIn_2D(:,:) /= 0) ObjFather%TotSonIn_2D(:,:) = Me%Father%TotSonIn_2D
            
        endif

    end subroutine Nudging_average_offline
    
    !>@author Joao Sobrinho +Atlantic
    !>@Brief
    !>Computes flux and velocity from an upscaling domain at the PD land boundary crossed by the CD
    !>@param[in] VelocityID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D
    subroutine Flux_Velocity_Offline(VelocityID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D)
        !Arguments-------------------------------------------------------------
        real, dimension(:, :, :), pointer, optional, intent(INOUT)  :: FatherMatrix
        real, dimension(:, :, :), pointer, optional, intent(IN)     :: SonMatrix
        real, dimension(:, :),    pointer, optional, intent(INOUT)  :: FatherMatrix2D
        real, dimension(:, :),    pointer, optional, intent(IN)     :: SonMatrix2D
        integer,                                     intent(IN)     :: VelocityID
        !Local-----------------------------------------------------------------
        type (T_TwoWay), pointer                                    :: ObjFather
        integer                                                     :: ready_father, i, j, k, ison,json,kson
        integer                                                     :: ILBSon, IUBSon, JLBSon, JUBSon, KUBSon, KLBSon
        integer                                                     :: FatherI_max, FatherJ_max, FatherK_max
        integer                                                     :: FatherI_min, FatherJ_min, FatherK_min
        integer                                                     :: Number_Cells, ifather, jfather
        integer                                                     :: CellFaceOpen, belongs
        integer                                                     :: jfather_Water, ifather_Water
        !integer                                                     :: Enter_Flag
        real, dimension(:, :, :), pointer                           :: SonArea_U, SonArea_V, Velocity
        logical                                                     :: FoundFirstColumn, FoundFirstLine
        integer, dimension(:, :, :), pointer                        :: SonMask, FatherMask, KLink
        integer, dimension(:, :), pointer                           :: ILink, JLink
        type (T_Size3D)                                             :: SizeAux
        !Begin ----------------------------------------------------------------
        if (present(FatherMatrix)) then
            !3D
            ILBSon = Me%WorkSize%ILB; JLBSon = Me%WorkSize%JLB; KUBSon = Me%WorkSize%KUB
            IUBSon = Me%WorkSize%IUB; JUBSon = Me%WorkSize%JUB; KLBSon = Me%WorkSize%KLB
            
            SonArea_U  => Me%External_Var%AreaU
            SonArea_V  => Me%External_Var%AreaV
            SonMask    => Me%External_Var%Open3D
            FatherMask => Me%Father%External_Var%Open3D
            ILink      => Me%External_Var%IZ
            JLink      => Me%External_Var%JZ
            KLink      => Me%External_Var%KZ
        
            FatherI_max = maxval(ILink)             ;   FatherI_min = minval(ILink, MASK=ILink .GT. 0.0)
            FatherJ_max = maxval(JLink)             ;   FatherJ_min = minval(JLink, MASK=JLink .GT. 0.0)
            FatherK_max = maxval(KLink)             ;   FatherK_min = minval(KLink, MASK=KLink .GT. 0.0)
            
            SizeAux%IUB = FatherI_max               ;   SizeAux%ILB = FatherI_min
            SizeAux%JUB = FatherJ_max               ;   SizeAux%JLB = FatherJ_min
            SizeAux%KUB = Me%Father%WorkSize%KUB    ;   SizeAux%KLB = Me%Father%WorkSize%KLB
            
            if (Me%Discharge%init .and. (.not. Me%Discharge%init_U) .and. (.not. Me%Discharge%init_V)) then
                Me%Discharge%init = .false.
            endif
            
            if (VelocityID == VelocityU_) then
                if (Me%Discharge%init) then
                    if (Me%Discharge%FirstTimeU) then
                        call SetMatrixValueAllocatable(Me%Father%Discharge%VelocityU_Prev, SizeAux, 0.0)
                        Velocity => Me%Father%Discharge%VelocityU_Prev
                        Me%Discharge%FirstTimeU = .false.
                    else
                        call SetMatrixValueAllocatable(Me%Father%Discharge%VelocityU_Next, SizeAux, 0.0)
                        Velocity => Me%Father%Discharge%VelocityU_Next
                        Me%Discharge%init_U = .false.
                    endif
                else
                    do k = Me%Father%WorkSize%KLB, Me%Father%WorkSize%KUB
                    do j = FatherJ_min, FatherJ_max
                    do i = FatherI_min, FatherI_max
                        Me%Father%Discharge%VelocityU_Prev(i,j,k) = Me%Father%Discharge%VelocityU_Next(i,j,k)
                        Me%Father%Discharge%VelocityU_Next(i,j,k) = 0.0
                        FatherMatrix(i, j, k) = 0.0
                    enddo
                    enddo
                    enddo
                    Velocity => Me%Father%Discharge%VelocityU_Next
                endif
            elseif (VelocityID == VelocityV_) then
                if (Me%Discharge%init) then
                    if (Me%Discharge%FirstTimeV) then
                        call SetMatrixValueAllocatable(Me%Father%Discharge%VelocityV_Prev, SizeAux, 0.0)
                        Velocity => Me%Father%Discharge%VelocityV_Prev
                        Me%Discharge%FirstTimeV = .false.
                    else
                        call SetMatrixValueAllocatable(Me%Father%Discharge%VelocityV_Next, SizeAux, 0.0)
                        Velocity => Me%Father%Discharge%VelocityV_Next
                        Me%Discharge%init_V = .false.
                    endif
                else
                    do k = Me%Father%WorkSize%KLB, Me%Father%WorkSize%KUB
                    do j = FatherJ_min, FatherJ_max
                    do i = FatherI_min, FatherI_max
                        Me%Father%Discharge%VelocityV_Prev(i,j,k) = Me%Father%Discharge%VelocityV_Next(i,j,k)
                        Me%Father%Discharge%VelocityV_Next(i,j,k) = 0.0
                        FatherMatrix(i, j, k) = 0.0
                    enddo
                    enddo
                    enddo
                    Velocity => Me%Father%Discharge%VelocityV_Next
                endif
            endif
            
            if (VelocityID == VelocityU_) then
                do k = FatherK_min, Me%Father%WorkSize%KUB
                do j = FatherJ_min, FatherJ_max
                do i = FatherI_min, FatherI_max
                    Number_Cells = 0
                    if (FatherMask(i, j , FatherK_max) == 1 .and. FatherMask(i, j + 1, FatherK_max) == 0) then
                        FoundFirstColumn = .false.
                        !land on the right side
                do1:    do json = JLBSon, JUBSon
                            if (FoundFirstColumn) then
                                !Only doing this computation for the first son column of cells inside PD land Cell
                                exit do1
                            endif
                        do ison = ILBSon, IUBSon
                            
                            ifather = ILink(ison, json)
                            jfather = JLink(ison, json)
                            !Check if son cell is inside PD land cell
                            belongs = 1 + abs(ifather - i) + abs(jfather - (j + 1))
                            CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison, json - 1, KUBSon)
                            if (belongs == 1 .and. CellFaceOpen == 2) then
                                FoundFirstColumn = .true.
                                do kson = KLBSon, KUBSon
                                    !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                    !because Klink of a land cell is FillValueReal
                                    if (KLink(ison,json-1,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                        Number_Cells = Number_Cells + 1
                                        !m3/s
                                        FatherMatrix(i, j, k) = FatherMatrix(i, j, k) &
                                                                - SonMatrix(ison,json,kson) * SonArea_U(ison,json,kson)
                                        Velocity(i, j, k) = Velocity(i, j, k) + SonMatrix(ison,json,kson)
                                    endif
                                enddo
                                !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                !because Klink of a land cell is FillValueReal
                                !This code below considers that both PD and CD use the same geometry implementation.
                                !Enter_Flag=SonMask(ison, json-1, k)+SonMask(ison, json, k)+FatherMask(i, j, k)
                                !if (Enter_Flag == 3) then
                                !    Number_Cells = Number_Cells + 1
                                !    !m3/s
                                !    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) &
                                !                            - SonMatrix(ison,json+1,k) * SonArea_U(ison,json+1,k)
                                !    Velocity(i, j, k) = Velocity(i, j, k) + SonMatrix(ison,json+1,k)
                                !endif
                            endif
                        enddo
                        enddo do1
                    end if    
                    if (FatherMask(i, j , FatherK_max) == 1 .and. FatherMask(i, j - 1, FatherK_max) == 0) then
                        FoundFirstColumn = .false.
                        !land on the left side
                do2:    do json = JLBSon, JUBSon
                            if (FoundFirstColumn) then
                                !Only doing this computation for the first son column of cells inside PD land Cell
                                exit do2
                            endif
                        do ison = ILBSon, IUBSon
                            
                            ifather = ILink(ison, json)
                            jfather = JLink(ison, json)
                            jfather_Water = jfather + 1
                            !Check if son cell is inside PD land cell. Must select CD open cell nearest to PD land cell
                            belongs = 1 + abs(ifather - i) + abs(jfather - (j - 1)) &
                                        + abs(jfather_Water - JLink(ison, json + 1))
                            CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison, json + 1, KUBSon)
                            if (belongs == 1 .and. CellFaceOpen == 2) then
                                FoundFirstColumn = .true.
                                do kson = KLBSon, KUBSon
                                    !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                    !because Klink of a land cell is FillValueReal
                                    if (KLink(ison,json-1,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                        Number_Cells = Number_Cells + 1
                                        !m3/s
                                        FatherMatrix(i, j, k) = FatherMatrix(i, j, k) &
                                                                + SonMatrix(ison,json,kson) * SonArea_U(ison,json,kson)
                                        Velocity(i, j, k) = Velocity(i, j, k) + SonMatrix(ison,json,kson)
                                    endif
                                enddo
                                !Enter_Flag=SonMask(ison, json-1, k)+SonMask(ison, json+1, k)+FatherMask(i, j, k)
                                !if (Enter_Flag == 3) then
                                !    Number_Cells = Number_Cells + 1
                                !    !m3/s
                                !    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) &
                                !                            + SonMatrix(ison,json-1,k) * SonArea_U(ison,json-1,k)
                                !    Velocity(i, j, k) = Velocity(i, j, k) + SonMatrix(ison,json-1,k)
                                !endif
                            endif
                        enddo
                        enddo do2
                    endif
                    if (Number_Cells > 0)Velocity(i, j, k) = Velocity(i, j, k) / Number_Cells
                enddo
                enddo
                enddo
            elseif (VelocityID == VelocityV_) then
                do k = FatherK_min, Me%Father%WorkSize%KUB
                do j = FatherJ_min, FatherJ_max
                do i = FatherI_min, FatherI_max
                    Number_Cells = 0
                    if (FatherMask(i, j, FatherK_max) == 1 .and. FatherMask(i + 1, j, FatherK_max) == 0) then
                        FoundFirstLine = .false.
                        !land to the North
                do3:    do ison = ILBSon, IUBSon
                            if (FoundFirstLine) then
                                !Only doing this computation for the first son column of cells inside PD land Cell
                                exit do3
                            endif
                        do json = JLBSon, JUBSon
                            
                            ifather = ILink(ison, json)
                            jfather = JLink(ison, json)
                            !Check if son cell is inside PD land cell
                            belongs = 1 + abs(ifather - (i + 1)) + abs(jfather - j)
                            CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison - 1, json, KUBSon)
                            if (belongs == 1 .and. CellFaceOpen == 2) then
                                FoundFirstLine = .true.
                                do kson = KLBSon, KUBSon
                                    !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                    !because Klink of a land cell is FillValueReal
                                    if (KLink(ison - 1,json,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                        Number_Cells = Number_Cells + 1
                                        !m3/s
                                        FatherMatrix(i, j, k) = FatherMatrix(i, j, k) &
                                                                - SonMatrix(ison,json,kson) * SonArea_V(ison,json,kson)
                                        Velocity(i, j, k) = Velocity(i, j, k) + SonMatrix(ison,json,kson)
                                    endif
                                enddo
                                !Enter_Flag=SonMask(ison-1, json, k)+SonMask(ison+1, json, k)+FatherMask(i, j, k)
                                !if (Enter_Flag == 3) then
                                !    Number_Cells = Number_Cells + 1
                                !    !m3/s
                                !    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) &
                                !                        - SonMatrix(ison+1,json,k) * SonArea_V(ison+1,json,k)
                                !    Velocity(i, j, k) = Velocity(i, j, k) + SonMatrix(ison+1,json,k)
                                !endif
                            endif
                        enddo
                        enddo do3
                    end if    
                    if (FatherMask(i, j, FatherK_max) == 1 .and. FatherMask(i - 1, j, FatherK_max) == 0) then
                        FoundFirstLine = .false.
                        !land to the South
                do4:    do ison = ILBSon, IUBSon
                            if (FoundFirstLine) then
                                !Only doing this computation for the first son column of cells inside PD land Cell
                                exit do4
                            endif
                        do json = JLBSon, JUBSon
                            
                            ifather = ILink(ison, json)
                            jfather = JLink(ison, json)
                            ifather_Water = ifather + 1
                            !Check if son cell is inside PD land cell. Must select CD open cell nearest to PD land cell
                            belongs = 1 + abs(ifather - (i - 1)) + abs(jfather - j) &
                                        + abs(ifather_Water - ILink(ison + 1, json))
                            CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison + 1, json, KUBSon)
                            if (belongs == 1 .and. CellFaceOpen == 2) then
                                FoundFirstLine = .true.
                                do kson = KLBSon, KUBSon
                                    !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                    !because Klink of a land cell is FillValueReal
                                    if (KLink(ison + 1,json,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                        Number_Cells = Number_Cells + 1
                                        !m3/s
                                        FatherMatrix(i, j, k) = FatherMatrix(i, j, k) &
                                                                + SonMatrix(ison,json,kson) * SonArea_V(ison,json,kson)
                                        Velocity(i, j, k) = Velocity(i, j, k) + SonMatrix(ison,json,kson)
                                    endif
                                enddo
                                !Enter_Flag=SonMask(ison-1, json, k)+SonMask(ison+1, json, k)+FatherMask(i, j, k)
                                !if (Enter_Flag == 3) then
                                !    Number_Cells = Number_Cells + 1
                                !    !m3/s
                                !    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) &
                                !                            + SonMatrix(ison-1,json,k) * SonArea_V(ison-1,json,k)
                                !    Velocity(i, j, k) = Velocity(i, j, k) + SonMatrix(ison-1,json,k)
                                !endif
                            endif
                        enddo
                        enddo do4
                    endif
                    if (Number_Cells > 0) Velocity(i, j, k) = Velocity(i, j, k) / Number_Cells
                enddo
                enddo
                enddo
            endif
            
            call ReadyFather(Me%Father%InstanceID, ObjFather, ready_father)
            
            if (VelocityID == VelocityU_) then
                if (.not. allocated(ObjFather%Discharge%VelocityU)) then
                    allocate(ObjFather%Discharge%VelocityU( ObjFather%WorkSize%ILB:ObjFather%WorkSize%IUB, &
                                                            ObjFather%WorkSize%JLB:ObjFather%WorkSize%JUB, &
                                                            ObjFather%WorkSize%KLB:ObjFather%WorkSize%KUB))
                endif
                
            elseif (VelocityID == VelocityV_) then
                if (.not. allocated(ObjFather%Discharge%VelocityV)) then
                    allocate(ObjFather%Discharge%VelocityV( ObjFather%WorkSize%ILB:ObjFather%WorkSize%IUB, &
                                                            ObjFather%WorkSize%JLB:ObjFather%WorkSize%JUB, &
                                                            ObjFather%WorkSize%KLB:ObjFather%WorkSize%KUB))
                endif
            endif
            
            nullify(SonArea_U, SonArea_V, SonMask, ILink, JLink, KLink, Velocity, FatherMask)
        else
            write (*,*) 'Upscaling discharge not yet ready for 2D'
            stop
        endif

    end subroutine Flux_Velocity_Offline
    
    !>@author Joao Sobrinho +Atlantic
    !>@Brief
    !>Recieves temporal information in order to perform interpolation on Upscaling discharge velocities
    !>@param[in] ActualTime, Size, Time1, Time2, PointsToFill3D
    subroutine InterpolUpscaling_Velocity (TwoWayID, PropertyNumber, ActualTime, Time1, Time2, PointsToFill3D, STAT)
        !Arguments-------------------------------------------------------------
        integer ,          intent(IN)                   :: TwoWayID, PropertyNumber
        type(T_Time),      intent(IN)                   :: ActualTime, Time1, Time2
        integer, dimension(:, :, :), pointer            :: PointsToFill3D
        integer,           optional, intent(OUT)        :: STAT
        !Local-----------------------------------------------------------------
        integer, dimension(:, :   ), pointer            :: IZ, JZ
        integer                                         :: STAT_CALL, ready_, STAT_
        type (T_Size3D)                                 :: SizeAux
        !Begin-----------------------------------------------------------------
        call Ready(TwoWayID, ready_)

        if ((ready_ .EQ. IDLE_ERR_  ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call GetHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid, &
                                    ILinkZ = IZ, JLinkZ = JZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InterpolUpscaling_Velocity - Failed to get Son-Father link matrixes'
            
            SizeAux%IUB = maxval(IZ)                ;   SizeAux%ILB = minval(IZ, MASK=IZ .GT. 0.0)
            SizeAux%JUB = maxval(JZ)                ;   SizeAux%JLB = minval(JZ, MASK=JZ .GT. 0.0)
            SizeAux%KUB = Me%Father%WorkSize%KUB    ;   SizeAux%KLB = Me%Father%WorkSize%KLB
            
            if (PropertyNumber == VelocityU_) then
                call SetMatrixValueAllocatable(Me%Father%Discharge%VelocityU, SizeAux, 0.0)
                
                call InterpolateMatrix3DInTime_Alloc(   ActualTime      = ActualTime,                          &
                                                        Size            = Me%Father%WorkSize,                  &
                                                        Time1           = Time1,                               &
                                                        Matrix1         = Me%Father%Discharge%VelocityU_Prev,  &
                                                        Time2           = Time2,                               &
                                                        Matrix2         = Me%Father%Discharge%VelocityU_Next,  &
                                                        MatrixOut       = Me%Father%Discharge%VelocityU,       &
                                                        PointsToFill3D  = PointsToFill3D)
            elseif (PropertyNumber == VelocityV_) then
                call SetMatrixValueAllocatable(Me%Father%Discharge%VelocityV, SizeAux, 0.0)
                
                call InterpolateMatrix3DInTime_Alloc(   ActualTime      = ActualTime,                          &
                                                        Size            = Me%Father%WorkSize,                  &
                                                        Time1           = Time1,                               &
                                                        Matrix1         = Me%Father%Discharge%VelocityV_Prev,  &
                                                        Time2           = Time2,                               &
                                                        Matrix2         = Me%Father%Discharge%VelocityV_Next,  &
                                                        MatrixOut       = Me%Father%Discharge%VelocityV,       &
                                                        PointsToFill3D  = PointsToFill3D)
            endif
            
            call UngetHorizontalGrid(Me%ObjHorizontalGrid, IZ,    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InterpolUpscaling_Velocity-TwoWay-ERR01'
            call UngetHorizontalGrid(Me%ObjHorizontalGrid, JZ,    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InterpolUpscaling_Velocity-TwoWay-ERR01'
            
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if
        if (present(STAT)) STAT = STAT_
        
    end subroutine InterpolUpscaling_Velocity
    
    !>@author Joao Sobrinho +Atlantic
    !>@Brief
    !>Computes concentration from an upscaling domain at the PD land boundary crossed by the CD
    !>@param[in] FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D
    subroutine DischageConc_Offline(FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D)
        !Arguments-------------------------------------------------------------
        real, dimension(:, :, :), pointer, optional, intent(INOUT)  :: FatherMatrix
        real, dimension(:, :, :), pointer, optional, intent(IN)     :: SonMatrix
        real, dimension(:, :),    pointer, optional, intent(INOUT)  :: FatherMatrix2D
        real, dimension(:, :),    pointer, optional, intent(IN)     :: SonMatrix2D
        !Local-----------------------------------------------------------------
        integer                                                     :: i, j, k, ison,json,kson
        integer                                                     :: ILBSon, IUBSon, JLBSon, JUBSon, KUBSon, KLBSon
        integer                                                     :: FatherI_max, FatherJ_max, FatherK_max
        integer                                                     :: FatherI_min, FatherJ_min, FatherK_min
        integer                                                     :: Number_Cells, ifather, jfather
        integer                                                     :: CellFaceOpen, belongs
        integer                                                     :: jfather_Water, ifather_Water
        logical                                                     :: FoundFirstColumn, FoundFirstLine
        integer, dimension(:, :, :), pointer                        :: SonMask, FatherMask, KLink
        integer, dimension(:, :), pointer                           :: ILink, JLink
        type (T_Size3D)                                             :: SizeAux
        !Begin ----------------------------------------------------------------
        if (present(FatherMatrix)) then
            !3D
            ILBSon = Me%WorkSize%ILB; JLBSon = Me%WorkSize%JLB; KLBSon = Me%WorkSize%KLB
            IUBSon = Me%WorkSize%IUB; JUBSon = Me%WorkSize%JUB; KUBSon = Me%WorkSize%KUB
            
            SonMask    => Me%External_Var%Open3D
            FatherMask => Me%Father%External_Var%Open3D
            ILink      => Me%External_Var%IZ
            JLink      => Me%External_Var%JZ
            KLink      => Me%External_Var%KZ
        
            FatherI_max = maxval(ILink)             ;  FatherI_min = minval(ILink, MASK=ILink .GT. 0.0)
            FatherJ_max = maxval(JLink)             ;  FatherJ_min = minval(JLink, MASK=JLink .GT. 0.0)
            FatherK_max = maxval(KLink)             ;  FatherK_min = minval(KLink, MASK=KLink .GT. 0.0)
            
            SizeAux%IUB = FatherI_max               ;  SizeAux%ILB = FatherI_min
            SizeAux%JUB = FatherJ_max               ;  SizeAux%JLB = FatherJ_min
            SizeAux%KUB = Me%Father%WorkSize%KUB    ;  SizeAux%KLB = Me%Father%WorkSize%KLB
            
            call SetMatrixValue(FatherMatrix, SizeAux, 0.0)
            
            do k = FatherK_min, Me%Father%WorkSize%KUB
            do j = FatherJ_min, FatherJ_max
            do i = FatherI_min, FatherI_max
                Number_Cells = 0
                if (FatherMask(i, j , FatherK_max) == 1 .and. FatherMask(i, j + 1, FatherK_max) == 0) then
                    FoundFirstColumn = .false.
                    FoundFirstLine = .false.
                    !land on the right side
            do1:    do json = JLBSon, JUBSon
                        if (FoundFirstColumn) then
                            !Only doing this computation for the first son column of cells inside PD land Cell
                            exit do1
                        endif
                    do ison = ILBSon, IUBSon
                            
                        ifather = ILink(ison, json)
                        jfather = JLink(ison, json)
                        !Check if son cell is inside PD land cell
                        belongs = 1 + abs(ifather - i) + abs(jfather - (j + 1))
                        CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison, json - 1, KUBSon)
                        if (belongs == 1 .and. CellFaceOpen == 2) then
                            FoundFirstColumn = .true.
                            do kson = KLBSon, KUBSon
                                !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                !because Klink of a land cell is FillValueReal
                                if (KLink(ison,json-1,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                    Number_Cells = Number_Cells + 1
                                    !m3/s
                                    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,kson)
                                endif
                            enddo
                            !This code below is only for when both PD and CD have the same geometry implementation
                            !if (SonMask(ison, json, k) == 1 .and. SonMask(ison, json-1, k) == 1 .and. FatherMask(i, j, k) == 1) then
                            !    Number_Cells = Number_Cells + 1
                            !    !m3/s
                            !    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,k)
                            !endif
                        endif
                    enddo
                    enddo do1
                end if    
                if (FatherMask(i, j , FatherK_max) == 1 .and. FatherMask(i, j - 1, FatherK_max) == 0) then
                    FoundFirstColumn = .false.
                    FoundFirstLine = .false.
                    !land on the left side
            do2:    do json = JLBSon, JUBSon
                        if (FoundFirstColumn) then
                            !Only doing this computation for the first son column of cells inside PD land Cell
                            exit do2
                        endif
                    do ison = ILBSon, IUBSon
                            
                        ifather = ILink(ison, json)
                        jfather = JLink(ison, json)
                        jfather_Water = jfather + 1
                        !Check if son cell is inside PD land cell. Must select CD open cell nearest to PD land cell
                        belongs = 1 + abs(ifather - i) + abs(jfather - (j - 1)) &
                                    + abs(jfather_Water - JLink(ison, json + 1))
                        CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison, json + 1, KUBSon)
                        if (belongs == 1 .and. CellFaceOpen == 2) then
                            FoundFirstColumn = .true.
                            do kson = KLBSon, KUBSon
                                !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                !because Klink of a land cell is FillValueReal
                                if (KLink(ison,json+1,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                    Number_Cells = Number_Cells + 1
                                    !m3/s
                                    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,kson)
                                endif
                            enddo
                            !if (SonMask(ison, json, k) == 1 .and. SonMask(ison, json+1, k) == 1 .and. FatherMask(i, j, k) == 1) then
                            !    Number_Cells = Number_Cells + 1
                            !    !m3/s
                            !    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,k)
                            !endif
                        endif
                    enddo
                enddo do2
                end if
                if (FatherMask(i, j, FatherK_max) == 1 .and. FatherMask(i + 1, j, FatherK_max) == 0) then
                    FoundFirstLine = .false.
                    !land to the North
            do3:    do ison = ILBSon, IUBSon
                        if (FoundFirstLine) then
                            !Only doing this computation for the first son line of cells inside PD land Cell
                            exit do3
                        endif
                    do json = JLBSon, JUBSon
                            
                        ifather = ILink(ison, json)
                        jfather = JLink(ison, json)
                        !Check if son cell is inside PD land cell
                        belongs = 1 + abs(ifather - (i + 1)) + abs(jfather - j)
                        CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison - 1, json, KUBSon)
                        if (belongs == 1 .and. CellFaceOpen == 2) then
                            FoundFirstLine = .true.
                            do kson = KLBSon, KUBSon
                                !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                !because Klink of a land cell is FillValueReal
                                if (KLink(ison - 1,json,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                    Number_Cells = Number_Cells + 1
                                    !m3/s
                                    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,kson)
                                endif
                            enddo
                            !if (SonMask(ison, json, k) == 1 .and. SonMask(ison-1, json, k) == 1 .and. FatherMask(i, j, k) == 1) then
                            !    Number_Cells = Number_Cells + 1
                            !    !m3/s
                            !    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,k)
                            !endif
                        endif
                    enddo
                    enddo do3
                end if
                if (FatherMask(i, j, FatherK_max) == 1 .and. FatherMask(i - 1, j, FatherK_max) == 0) then
                    FoundFirstLine = .false.
                    !land to the South
            do4:    do ison = ILBSon, IUBSon
                        if (FoundFirstLine) then
                            !Only doing this computation for the first son line of cells inside PD land Cell
                            exit do4
                        endif
                    do json = JLBSon, JUBSon
                            
                        ifather = ILink(ison, json)
                        jfather = JLink(ison, json)
                        ifather_Water = ifather + 1
                        !Check if son cell is inside PD land cell. Must select CD open cell nearest to PD land cell
                        belongs = 1 + abs(ifather - (i - 1)) + abs(jfather - j) &
                                    + abs(ifather_Water - ILink(ison + 1, json))
                        CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison + 1, json, KUBSon)
                        if (belongs == 1 .and. CellFaceOpen == 2) then
                            FoundFirstLine = .true.
                            do kson = KLBSon, KUBSon
                                !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                !because Klink of a land cell is FillValueReal
                                if (KLink(ison + 1,json,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                    Number_Cells = Number_Cells + 1
                                    !m3/s
                                    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,kson)
                                endif
                            enddo
                            !if (SonMask(ison, json, k) == 1 .and. SonMask(ison+1, json, k) == 1 .and. FatherMask(i, j, k) == 1) then
                            !    Number_Cells = Number_Cells + 1
                            !    !m3/s
                            !    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,k)
                            !endif
                        endif
                    enddo
                    enddo do4
                endif
                if (Number_Cells > 0) then
                    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) / Number_Cells
                endif
            enddo
            enddo
            enddo
            nullify(SonMask, FatherMask, ILink, JLink, KLink)
        else
            write (*,*) 'Upscaling discharge not yet ready for 2D'
            stop
        endif

    end subroutine DischageConc_Offline
    
    !>@author Joao Sobrinho +Atlantic
    !>@Brief
    !>Computes concentration from an Online upscaling domain at the PD land boundary crossed by the CD
    !>@param[in] FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D
    subroutine DischageConc_Online(FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D)
        !Arguments-------------------------------------------------------------
        real, dimension(:, :, :), pointer, optional, intent(INOUT)  :: FatherMatrix
        real, dimension(:, :, :), pointer, optional, intent(IN)     :: SonMatrix
        real, dimension(:, :),    pointer, optional, intent(INOUT)  :: FatherMatrix2D
        real, dimension(:, :),    pointer, optional, intent(IN)     :: SonMatrix2D
        !Local-----------------------------------------------------------------
        integer                                                     :: i, j, k, ison,json,kson
        integer                                                     :: ILBSon, IUBSon, JLBSon, JUBSon, KUBSon, KLBSon
        integer                                                     :: FatherI_max, FatherJ_max, FatherK_max
        integer                                                     :: FatherI_min, FatherJ_min, FatherK_min
        integer                                                     :: Number_Cells, ifather, jfather
        integer                                                     :: CellFaceOpen, belongs
        integer                                                     :: jfather_Water, ifather_Water
        logical                                                     :: FoundFirstColumn, FoundFirstLine
        integer, dimension(:, :, :), pointer                        :: SonMask, FatherMask, KLink
        integer, dimension(:, :), pointer                           :: ILink, JLink
        type (T_Size3D)                                             :: SizeAux
        !Begin ----------------------------------------------------------------
        if (present(FatherMatrix)) then
            !3D
            ILBSon = Me%WorkSize%ILB; JLBSon = Me%WorkSize%JLB; KLBSon = Me%WorkSize%KLB
            IUBSon = Me%WorkSize%IUB; JUBSon = Me%WorkSize%JUB; KUBSon = Me%WorkSize%KUB
            
            SonMask    => Me%External_Var%Open3D
            FatherMask => Me%Father%External_Var%Open3D
            ILink      => Me%External_Var%IZ
            JLink      => Me%External_Var%JZ
            KLink      => Me%External_Var%KZ
        
            FatherI_max = maxval(ILink)             ;  FatherI_min = minval(ILink, MASK=ILink .GT. 0.0)
            FatherJ_max = maxval(JLink)             ;  FatherJ_min = minval(JLink, MASK=JLink .GT. 0.0)
            FatherK_max = maxval(KLink)             ;  FatherK_min = minval(KLink, MASK=KLink .GT. 0.0)
            
            SizeAux%IUB = FatherI_max               ;  SizeAux%ILB = FatherI_min
            SizeAux%JUB = FatherJ_max               ;  SizeAux%JLB = FatherJ_min
            SizeAux%KUB = Me%Father%WorkSize%KUB    ;  SizeAux%KLB = Me%Father%WorkSize%KLB
            
            call SetMatrixValue(FatherMatrix, SizeAux, 0.0)
            
            do k = FatherK_min, FatherK_max
            do j = FatherJ_min, FatherJ_max
            do i = FatherI_min, FatherI_max
                Number_Cells = 0
                if (FatherMask(i, j , FatherK_max) == 1 .and. FatherMask(i, j + 1, FatherK_max) == 0) then
                    FoundFirstColumn = .false.
                    !land on the right side
            do1:    do json = JLBSon, JUBSon
                        if (FoundFirstColumn) then
                            !Only doing this computation for the first son column of cells inside PD land Cell
                            exit do1
                        endif
                    do ison = ILBSon, IUBSon
                            
                        ifather = ILink(ison, json)
                        jfather = JLink(ison, json)
                        !Check if son cell is inside PD land cell
                        belongs = 1 + abs(ifather - i) + abs(jfather - (j + 1))
                        CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison, json - 1, KUBSon)
                        if (belongs == 1 .and. CellFaceOpen == 2) then
                            FoundFirstColumn = .true.
                            do kson = KLBSon, KUBSon
                                !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                !because Klink of a land cell is FillValueReal
                                if (KLink(ison,json-1,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                    Number_Cells = Number_Cells + 1
                                    !m3/s
                                    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,kson)
                                endif
                            enddo
                        endif
                    enddo
                    enddo do1
                end if        
                if (FatherMask(i, j , FatherK_max) == 1 .and. FatherMask(i, j - 1, FatherK_max) == 0) then
                    FoundFirstColumn = .false.
                    !land on the left side
            do2:    do json = JLBSon, JUBSon
                        if (FoundFirstColumn) then
                            !Only doing this computation for the first son column of cells inside PD land Cell
                            exit do2
                        endif
                    do ison = ILBSon, IUBSon
                            
                        ifather = ILink(ison, json)
                        jfather = JLink(ison, json)
                        jfather_Water = jfather + 1
                        !Check if son cell is inside PD land cell. Must select CD open cell nearest to PD land cell
                        belongs = 1 + abs(ifather - i) + abs(jfather - (j - 1)) &
                                    + abs(jfather_Water - JLink(ison, json + 1))
                        CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison, json + 1, KUBSon)
                        if (belongs == 1 .and. CellFaceOpen == 2) then
                            FoundFirstColumn = .true.
                            do kson = KLBSon, KUBSon
                                !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                !because Klink of a land cell is FillValueReal
                                if (KLink(ison,json+1,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                    Number_Cells = Number_Cells + 1
                                    !m3/s
                                    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,kson)
                                endif
                            enddo
                        endif
                    enddo
                    enddo do2
                end if
                if (FatherMask(i, j, FatherK_max) == 1 .and. FatherMask(i + 1, j, FatherK_max) == 0) then
                    FoundFirstLine = .false.
                    !land to the North
            do3:    do ison = ILBSon, IUBSon
                        if (FoundFirstLine) then
                            !Only doing this computation for the first son line of cells inside PD land Cell
                            exit do3
                        endif
                    do json = JLBSon, JUBSon
                            
                        ifather = ILink(ison, json)
                        jfather = JLink(ison, json)
                        !Check if son cell is inside PD land cell
                        belongs = 1 + abs(ifather - (i + 1)) + abs(jfather - j)
                        CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison - 1, json, KUBSon)
                        if (belongs == 1 .and. CellFaceOpen == 2) then
                            FoundFirstLine = .true.
                            do kson = KLBSon, KUBSon
                                !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                !because Klink of a land cell is FillValueReal
                                if (KLink(ison - 1,json,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                    Number_Cells = Number_Cells + 1
                                    !m3/s
                                    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,kson)
                                endif
                            enddo
                        endif
                    enddo
                    enddo do3
                end if 
                if (FatherMask(i, j, FatherK_max) == 1 .and. FatherMask(i - 1, j, FatherK_max) == 0) then
                    FoundFirstLine = .false.
                    !land to the South
            do4:    do ison = ILBSon, IUBSon
                        if (FoundFirstLine) then
                            !Only doing this computation for the first son line of cells inside PD land Cell
                            exit do4
                        endif
                    do json = JLBSon, JUBSon
                            
                        ifather = ILink(ison, json)
                        jfather = JLink(ison, json)
                        ifather_Water = ifather + 1
                        !Check if son cell is inside PD land cell. Must select CD open cell nearest to PD land cell
                        belongs = 1 + abs(ifather - (i - 1)) + abs(jfather - j) &
                                    + abs(ifather_Water - ILink(ison + 1, json))
                        CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison + 1, json, KUBSon)
                        if (belongs == 1 .and. CellFaceOpen == 2) then
                            FoundFirstLine = .true.
                            do kson = KLBSon, KUBSon
                                !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                !because Klink of a land cell is FillValueReal
                                if (KLink(ison + 1,json,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                    Number_Cells = Number_Cells + 1
                                    !m3/s
                                    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,kson)
                                endif
                            enddo
                        endif
                    enddo
                    enddo do4
                endif
                if (Number_Cells > 0) then
                    FatherMatrix(i, j, k) = FatherMatrix(i, j, k) / Number_Cells
                endif
            enddo
            enddo
            enddo
            nullify(SonMask, FatherMask, ILink, JLink, KLink)
        else
            write (*,*) 'Upscaling discharge not yet ready for 2D'
            stop
        endif

    end subroutine DischageConc_Online
    
    !>@author Joao Sobrinho +Atlantic
    !>@Brief
    !>Computes flux and velocity from an online coupling upscaling domain at the PD land boundary crossed by the CD
    !>@param[in] VelocityID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D
    subroutine Flux_Velocity_Online(VelocityID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D)
        !Arguments-------------------------------------------------------------
        real, dimension(:, :, :), pointer, optional, intent(INOUT)  :: FatherMatrix
        real, dimension(:, :, :), pointer, optional, intent(IN)     :: SonMatrix
        real, dimension(:, :),    pointer, optional, intent(INOUT)  :: FatherMatrix2D
        real, dimension(:, :),    pointer, optional, intent(IN)     :: SonMatrix2D
        integer,                                     intent(IN)     :: VelocityID
        !Local-----------------------------------------------------------------
        type (T_TwoWay), pointer                                    :: ObjFather
        integer                                                     :: ready_father, i, j, k, ison,json,kson
        integer                                                     :: ILBSon, IUBSon, JLBSon, JUBSon, KUBSon, KLBSon
        integer                                                     :: FatherI_max, FatherJ_max, FatherK_max
        integer                                                     :: FatherI_min, FatherJ_min, FatherK_min
        integer                                                     :: Number_Cells, ifather, jfather
        integer                                                     :: CellFaceOpen, belongs
        integer                                                     :: jfather_Water, ifather_Water
        real, dimension(:, :, :), pointer                           :: SonArea_U, SonArea_V
        logical                                                     :: FoundFirstColumn, FoundFirstLine
        integer, dimension(:, :, :), pointer                        :: SonMask, FatherMask, KLink
        integer, dimension(:, :), pointer                           :: ILink, JLink
        type (T_Size3D)                                             :: SizeAux
        !Begin ----------------------------------------------------------------
        if (present(FatherMatrix)) then
            !3D
            ILBSon = Me%WorkSize%ILB; JLBSon = Me%WorkSize%JLB; KUBSon = Me%WorkSize%KUB
            IUBSon = Me%WorkSize%IUB; JUBSon = Me%WorkSize%JUB; KLBSon = Me%WorkSize%KLB
            
            SonArea_U  => Me%External_Var%AreaU
            SonArea_V  => Me%External_Var%AreaV
            SonMask    => Me%External_Var%Open3D
            FatherMask => Me%Father%External_Var%Open3D
            ILink      => Me%External_Var%IZ
            JLink      => Me%External_Var%JZ
            KLink      => Me%External_Var%KZ
        
            FatherI_max = maxval(ILink)             ;   FatherI_min = minval(ILink, MASK=ILink .GT. 0.0)
            FatherJ_max = maxval(JLink)             ;   FatherJ_min = minval(JLink, MASK=JLink .GT. 0.0)
            FatherK_max = maxval(KLink)             ;   FatherK_min = minval(KLink, MASK=KLink .GT. 0.0)
            
            SizeAux%IUB = FatherI_max               ;   SizeAux%ILB = FatherI_min
            SizeAux%JUB = FatherJ_max               ;   SizeAux%JLB = FatherJ_min
            SizeAux%KUB = Me%Father%WorkSize%KUB    ;   SizeAux%KLB = Me%Father%WorkSize%KLB
            
            call ReadyFather(Me%Father%InstanceID, ObjFather, ready_father)
            
            if (VelocityID == VelocityU_) then
                call SetMatrixValue(FatherMatrix, SizeAux, 0.0)
                call SetMatrixValueAllocatable(ObjFather%Discharge%FlowMatrix, SizeAux, 0.0)
            elseif (VelocityID == VelocityV_) then
                call SetMatrixValue(FatherMatrix, SizeAux, 0.0)
            endif
            
            if (VelocityID == VelocityU_) then
                do k = FatherK_min, FatherK_max
                do j = FatherJ_min, FatherJ_max
                do i = FatherI_min, FatherI_max
                    Number_Cells = 0
                    if (FatherMask(i, j , FatherK_max) == 1 .and. FatherMask(i, j + 1, FatherK_max) == 0) then
                        FoundFirstColumn = .false.
                        !land on the right side
                do1:    do json = JLBSon, JUBSon
                            if (FoundFirstColumn) then
                                !Only doing this computation for the first son column of cells inside PD land Cell
                                exit do1
                            endif
                        do ison = ILBSon, IUBSon
                            
                            ifather = ILink(ison, json)
                            jfather = JLink(ison, json)
                            !Check if son cell is inside PD land cell
                            belongs = 1 + abs(ifather - i) + abs(jfather - (j + 1))
                            CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison, json - 1, KUBSon)
                            if (belongs == 1 .and. CellFaceOpen == 2) then
                                FoundFirstColumn = .true.
                               do kson = KLBSon, KUBSon
                                   !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                   !because Klink of a land cell is FillValueReal
                                   if (KLink(ison,json-1,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                        Number_Cells = Number_Cells + 1
                                        !m3/s
                                        ObjFather%Discharge%FlowMatrix(i, j, k) = &
                                        ObjFather%Discharge%FlowMatrix(i, j, k) - SonMatrix(ison,json,kson) &
                                                                                * SonArea_U(ison,json,kson)
                                        FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,kson)
                                   endif
                               enddo
                            endif
                        enddo
                        enddo do1
                    end if    
                    if (FatherMask(i, j , FatherK_max) == 1 .and. FatherMask(i, j - 1, FatherK_max) == 0) then
                        FoundFirstColumn = .false.
                        !land on the left side
                do2:    do json = JLBSon, JUBSon
                            if (FoundFirstColumn) then
                                !Only doing this computation for the first son column of cells inside PD land Cell
                                exit do2
                            endif
                        do ison = ILBSon, IUBSon
                            
                            ifather = ILink(ison, json)
                            jfather = JLink(ison, json)
                            jfather_Water = jfather + 1
                            !Check if son cell is inside PD land cell. Must select CD open cell nearest to PD land cell
                            belongs = 1 + abs(ifather - i) + abs(jfather - (j - 1)) &
                                        + abs(jfather_Water - JLink(ison, json + 1))
                            CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison, json + 1, KUBSon)
                            if (belongs == 1 .and. CellFaceOpen == 2) then
                                FoundFirstColumn = .true.
                               do kson = KLBSon, KUBSon
                                   !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                   !because Klink of a land cell is FillValueReal
                                   if (KLink(ison,json+1,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                        Number_Cells = Number_Cells + 1
                                        !m3/s
                                        ObjFather%Discharge%FlowMatrix(i, j, k) = &
                                        ObjFather%Discharge%FlowMatrix(i, j, k) + SonMatrix(ison,json,kson) &
                                                                                * SonArea_U(ison,json,kson)
                                        FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,kson)
                                   endif
                               enddo
                            endif
                        enddo
                        enddo do2
                    endif
                    if (Number_Cells > 0) FatherMatrix(i, j, k) = FatherMatrix(i, j, k) / Number_Cells
                enddo
                enddo
                enddo
            elseif (VelocityID == VelocityV_) then
                
                do k = FatherK_min, FatherK_max
                do j = FatherJ_min, FatherJ_max
                do i = FatherI_min, FatherI_max
                    Number_Cells = 0
                    if (FatherMask(i, j, FatherK_max) == 1 .and. FatherMask(i + 1, j, FatherK_max) == 0) then
                        FoundFirstLine = .false.
                        !land to the North
                do3:    do ison = ILBSon, IUBSon
                            if (FoundFirstLine) then
                                !Only doing this computation for the first son column of cells inside PD land Cell
                                exit do3
                            endif
                        do json = JLBSon, JUBSon
                            
                            ifather = ILink(ison, json)
                            jfather = JLink(ison, json)
                            !Check if son cell is inside PD land cell
                            belongs = 1 + abs(ifather - (i + 1)) + abs(jfather - j)
                            CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison - 1, json, KUBSon)
                            if (belongs == 1 .and. CellFaceOpen == 2) then
                                FoundFirstLine = .true.
                               do kson = KLBSon, KUBSon
                                   !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                   !because Klink of a land cell is FillValueReal
                                   if (KLink(ison - 1,json,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                        Number_Cells = Number_Cells + 1
                                        !m3/s
                                        ObjFather%Discharge%FlowMatrix(i, j, k) = &
                                        ObjFather%Discharge%FlowMatrix(i, j, k) - SonMatrix(ison,json,kson) &
                                                                                * SonArea_V(ison,json,kson)
                                        FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,kson)
                                   endif
                               enddo
                            endif
                        enddo
                        enddo do3
                    end if    
                    if (FatherMask(i, j, FatherK_max) == 1 .and. FatherMask(i - 1, j, FatherK_max) == 0) then
                        FoundFirstLine = .false.
                        !land to the South
                do4:    do ison = ILBSon, IUBSon
                            if (FoundFirstLine) then
                                !Only doing this computation for the first son column of cells inside PD land Cell
                                exit do4
                            endif
                        do json = JLBSon, JUBSon
                            
                            ifather = ILink(ison, json)
                            jfather = JLink(ison, json)
                            ifather_Water = ifather + 1
                            !Check if son cell is inside PD land cell. Must select CD open cell nearest to PD land cell
                            belongs = 1 + abs(ifather - (i - 1)) + abs(jfather - j) &
                                        + abs(ifather_Water - ILink(ison + 1, json))
                            CellFaceOpen = SonMask(ison, json, KUBSon) + SonMask(ison + 1, json, KUBSon)
                            if (belongs == 1 .and. CellFaceOpen == 2) then
                                FoundFirstLine = .true.
                               do kson = KLBSon, KUBSon
                                   !Check if CD cell in vertical belongs to PD cell (using left PD cell
                                   !because Klink of a land cell is FillValueReal
                                   if (KLink(ison + 1,json,kson) == k .and. SonMask(ison, json, kson) == 1) then
                                        Number_Cells = Number_Cells + 1
                                        !m3/s
                                        ObjFather%Discharge%FlowMatrix(i, j, k) = &
                                        ObjFather%Discharge%FlowMatrix(i, j, k) + SonMatrix(ison,json,kson) &
                                                                                * SonArea_V(ison,json,kson)
                                        FatherMatrix(i, j, k) = FatherMatrix(i, j, k) + SonMatrix(ison,json,kson)
                                   endif
                               enddo
                            endif
                        enddo
                        enddo do4
                    endif
                    if (Number_Cells > 0) FatherMatrix(i, j, k) = FatherMatrix(i, j, k) / Number_Cells
                enddo
                enddo
                enddo
            endif
            
            nullify(SonArea_U, SonArea_V, SonMask, ILink, JLink, KLink, FatherMask)
        else
            write (*,*) 'Upscaling discharge not yet ready for 2D'
            stop
        endif

    end subroutine Flux_Velocity_Online
    
    !------------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Calls Feedback routines present in mondule functions, based on the type of input matrixes - IWD Method
    !>@param[in] FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, LocalTimeDecay
    subroutine Nudging_IWD (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, LocalTimeDecay)
        !Arguments-------------------------------------------------------------
        integer, optional                                        :: VelocityID
        real, dimension(:, :, :), pointer, optional, intent(OUT) :: FatherMatrix
        real, dimension(:, :, :), pointer, optional, intent(IN)  :: SonMatrix
        real, dimension(:, :),    pointer, optional, intent(OUT) :: FatherMatrix2D
        real, dimension(:, :),    pointer, optional, intent(IN)  :: SonMatrix2D
        real                                                     :: LocalTimeDecay
        !Local-----------------------------------------------------------------

        if (present(FatherMatrix)) then

            if (present(VelocityID))then
                !If velocity U
                if (VelocityID == VelocityU_) then
                    call FeedBack_IWD_UV (FatherMatrix         = FatherMatrix,                           &
                                           SonMatrix            = SonMatrix,                              &
                                           Open3DFather         = Me%Father%External_Var%Open3D,           &
                                           Open3DSon            = Me%External_Var%Open3D,                  &
                                           FatherComputeFaces3D = Me%Father%External_Var%ComputeFaces3D_U, &
                                           SonComputeFaces3D    = Me%External_Var%ComputeFaces3D_U,        &
                                           SizeSon              = Me%WorkSize,                            &
                                           Ilink                = Me%External_Var%IU,                      &
                                           Jlink                = Me%External_Var%JU,                      &
                                           Connections          = Me%External_Var%IWD_Connections_U,         &
                                           Dist                 = Me%External_Var%IWD_Distances_U,          &
                                           DecayTime            = LocalTimeDecay,                         &
                                           DT                   = Me%Hydro%VelDT,                         &
                                           IgnoreOBCells        = GetPointer(Me%IgnoreOBCells),           &
                                           Nodes                = Me%External_Var%IWD_Nodes_U,      &
                                           IWDn                 = Me%Hydro%IWDn,                    &
                                           Nom                  = GetPointer(Me%Father%IWDNom),                 &
                                           Denom                = GetPointer(Me%Father%IWDDenom))
                !If velocity V
                else
                    call FeedBack_IWD_UV (FatherMatrix         = FatherMatrix,                           &
                                           SonMatrix            = SonMatrix,                              &
                                           Open3DFather         = Me%Father%External_Var%Open3D,           &
                                           Open3DSon            = Me%External_Var%Open3D,                  &
                                           FatherComputeFaces3D = Me%Father%External_Var%ComputeFaces3D_V, &
                                           SonComputeFaces3D    = Me%External_Var%ComputeFaces3D_V,        &
                                           SizeSon              = Me%WorkSize,                            &
                                           Ilink                = Me%External_Var%IV,                      &
                                           Jlink                = Me%External_Var%JV,                      &
                                           Connections          = Me%External_Var%IWD_Connections_U,         &
                                           Dist             = Me%External_Var%IWD_Distances_U,          &
                                           DecayTime            = LocalTimeDecay,                         &
                                           DT                   = Me%Hydro%VelDT,                         &
                                           IgnoreOBCells        = GetPointer(Me%IgnoreOBCells),           &
                                           Nodes                = Me%External_Var%IWD_Nodes_V,      &
                                           IWDn                 = Me%Hydro%IWDn,                    &
                                           Nom                  = GetPointer(Me%Father%IWDNom),                 &
                                           Denom                = GetPointer(Me%Father%IWDDenom))
                endif
            else
                !compute nudging Z type cell
                call FeedBack_IWD   (FatherMatrix     = FatherMatrix,                    &
                                      SonMatrix        = SonMatrix,                       &
                                      Open3DFather     = Me%Father%External_Var%Open3D,    &
                                      Open3DSon        = Me%External_Var%Open3D,           &
                                      SizeSon          = Me%WorkSize,                     &
                                      Ilink            = Me%External_Var%IZ,               &
                                      Jlink            = Me%External_Var%JZ,               &
                                      Connections      = Me%External_Var%IWD_Connections_Z,         &
                                      Dist             = Me%External_Var%IWD_Distances_Z,          &
                                      DecayTime        = LocalTimeDecay,                  &
                                      DT               = Me%Hydro%DT,                     &
                                      IgnoreOBCells    = GetPointer(Me%IgnoreOBCells),           &
                                      Nodes            = Me%External_Var%IWD_Nodes_Z,      &
                                      IWDn             = Me%Hydro%IWDn,                    &
                                      Nom              = GetPointer(Me%Father%IWDNom),                 &
                                      Denom            = GetPointer(Me%Father%IWDDenom))
            endif

        else

            call FeedBack_IWD_WL (FatherMatrix2D   = FatherMatrix2D,                  &
                                  SonMatrix2D      = SonMatrix2D,                     &
                                  Open3DFather     = Me%Father%External_Var%Open3D,    &
                                  Open3DSon        = Me%External_Var%Open3D,           &
                                  SizeSon          = Me%WorkSize,                     &
                                  Ilink            = Me%External_Var%IZ,               &
                                  Jlink            = Me%External_Var%JZ,               &
                                  Connections      = Me%External_Var%IWD_Connections_Z,         &
                                  Dist             = Me%External_Var%IWD_Distances_Z,          &
                                  DecayTime        = LocalTimeDecay,                  &
                                  DT               = Me%Hydro%DT,                     &
                                  IgnoreOBCells    = GetPointer(Me%IgnoreOBCells),           &
                                  Nodes            = Me%External_Var%IWD_Nodes_Z,      &
                                  IWDn             = Me%Hydro%IWDn,                    &
                                  Nom              = GetPointer(Me%Father%IWDNom),                 &
                                  Denom            = GetPointer(Me%Father%IWDDenom))

        endif

        !----------------------------------------------------------------------
    end subroutine Nudging_IWD

    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Routine responsible for the computation of discharge volume by upscaling
    !>@param[in] SonID, Father_old, Father, VelocityID, STAT
    subroutine UpscaleDischarge(SonID, Father_old, Father, VelocityID, STAT)
        !Arguments-------------------------------------------------------------
        integer                          , intent(IN)     :: SonID, VelocityID
        real, dimension(:, :, :), allocatable, intent(IN) :: Father_old
        real, dimension(:, :, :), pointer, intent(IN)     :: Father
        integer                          , intent(OUT)    :: STAT
        type (T_TwoWay), pointer                          :: ObjFather
        !Local-----------------------------------------------------------------
        integer                                           :: ready_, ready_father
        !----------------------------------------------------------------------
        STAT = UNKNOWN_
        call Ready(SonID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_))then
            if (VelocityID == VelocityU_) then
                call DischargeFluxU(FatherU_old = Father_old, FatherU = Father, AreaU = Me%Father%External_Var%AreaU, &
                                    Flow = Me%Father%DischargeCells%Flow, &
                                    DischargeConnection = Me%Father%DischargeCells%Z)
            elseif (VelocityID == VelocityV_) then
                call DischargeFluxV(FatherV_old = Father_old, FatherV = Father, AreaV = Me%Father%External_Var%AreaV, &
                                    Flow = Me%Father%DischargeCells%Flow, &
                                    DischargeConnection = Me%Father%DischargeCells%Z)
            else
                Write(*,*) 'Variable ID does not exist', VelocityID
                stop 'UpscaleDischarge - ModuleTwoWay'
            endif

            call ReadyFather (Me%Father%InstanceID, ObjFather, ready_father)

            ObjFather%DischargeCells%AuxFlow => Me%Father%DischargeCells%Flow
            ObjFather%DischargeCells%AuxConnections => Me%Father%DischargeCells%Z
            STAT = SUCCESS_
        else
            STAT = ready_
        endif

    end subroutine UpscaleDischarge
!---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Routine responsible for filling the concentration for each upscaling discharge cell
    !>@param[in] FatherID, Prop, PropVector, Flow, FlowVector, dI, dJ, dK, Kmin, Kmax, FirstTime
    subroutine UpscaleDischarge_WP (FatherID, Connections_Z, Prop, PropVector, Flow, FlowVector, dI, dJ, dK, Kmin, &
    Kmax, FirstTime)
        !Arguments-------------------------------------------------------------
        integer                          , intent(IN)     :: FatherID
        integer, dimension(:,:    ), pointer, intent(IN)  :: Connections_Z
        real(8), dimension(:, :, :), pointer, intent(IN)  :: Flow
        real,    dimension(:, :, :), pointer, intent(IN)  :: Prop
        real, dimension(:)      , pointer, intent(IN)     :: FlowVector, PropVector
        integer, dimension(:), pointer, intent(INOUT)     :: dI, dJ, dK, Kmin, Kmax
        logical                       , intent(IN)        :: FirstTime
        !Local-----------------------------------------------------------------
        integer                                           :: STAT_CALL, DischargeNumber, AuxCell, AuxKmin, AuxKmax, dis, &
                                                             nCells, i
        integer, dimension(:), pointer                    :: VectorI, VectorJ, VectorK
        !----------------------------------------------------------------------

        call GetDischargesNumber(FatherID, DischargeNumber, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'UpscaleDischarge_WP - failed to get discharges number of a father'

        AuxCell = 0
        do dis = 1 , DischargeNumber
            !write(*,*) 'Discharge number : ', dis
            call GetDischargeFlowDistribuiton(FatherID, DischargeIDNumber = dis, nCells = nCells, &
                                                VectorI = VectorI, VectorJ = VectorJ, VectorK = VectorK, &
                                                kmin = AuxKmin, kmax = AuxKmax, STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'UpscaleDischarge_WP - failed to get dischargeflowdistribution'

            !write(*,*) 'Discharge is associated : ', dis, VectorI(1), VectorJ(1)
            if (IsUpscaling(FatherID, dis)) then
                !using VectorI/J(1) because the value is the same for the entire vector (only the K value changes)
                if (DischargeIsAssociated (Connections_Z, VectorI(1), VectorJ(1))) then
                    if (FirstTime)then
                        do i = 1, nCells
                            AuxCell = AuxCell + 1
                            dI        (AuxCell) = VectorI(i)
                            dJ        (AuxCell) = VectorJ(i)
                            dK        (AuxCell) = VectorK(i)
                            Kmin      (AuxCell) = AuxKmin
                            Kmax      (AuxCell) = AuxKmax
                            FlowVector(AuxCell) = Flow(VectorI(i), VectorJ(i), VectorK(i))
                            PropVector(AuxCell) = Prop(VectorI(i), VectorJ(i), VectorK(i))
                        enddo
                    else
                        do i = 1, nCells
                            AuxCell = AuxCell + 1
                            PropVector(AuxCell) = Prop(VectorI(i), VectorJ(i), VectorK(i))
                        enddo
                    endif
                endif
            else
                !skip the discharge - will be taken caro of in the main body of the code, outside the two-way
                AuxCell = AuxCell + nCells
            end if
            
            call UnGetDischarges(FatherID, VectorI, STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'UpscaleDischarge_WP - failed to UnGetDischarges Vector I'

            call UnGetDischarges(FatherID, VectorJ, STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'UpscaleDischarge_WP - failed to UnGetDischarges Vector J'

            call UnGetDischarges(FatherID, VectorK, STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'UpscaleDischarge_WP - failed to UnGetDischarges Vector K'
            !if (DischargeIsAssociated (Connections_Z, VectorI(1), VectorJ(1))) then
            !    !write(*,*) 'Discharge is associated : ', dis, VectorI(1), VectorJ(1)
            !    if (IsUpscaling(FatherID, dis)) then
            !        if (FirstTime)then
            !            do i = 1, nCells
            !                AuxCell = AuxCell + 1
            !                dI        (AuxCell) = VectorI(i)
            !                dJ        (AuxCell) = VectorJ(i)
            !                dK        (AuxCell) = VectorK(i)
            !                Kmin      (AuxCell) = AuxKmin
            !                Kmax      (AuxCell) = AuxKmax
            !                FlowVector(AuxCell) = Flow(VectorI(i), VectorJ(i), VectorK(i))
            !                PropVector(AuxCell) = Prop(VectorI(i), VectorJ(i), VectorK(i))
            !            enddo
            !        else
            !            do i = 1, nCells
            !                AuxCell = AuxCell + 1
            !                PropVector(AuxCell) = Prop(VectorI(i), VectorJ(i), VectorK(i))
            !            enddo
            !        endif
            !    else
            !        AuxCell = AuxCell + nCells
            !    endif
            !endif
        enddo

    end subroutine UpscaleDischarge_WP

    !--------------------------------------------------------------------------
    !>@author Joao Sobrinho +Atlantic
    !>@Brief
    !>Fills discharge flow matrix for an offline upscaling discharge
    !>@param[in] TwoWayID, Flow, VelFather, VelSon, DecayTime, CoefCold, VelID, VelDT, STAT
    subroutine Offline_Upscaling_Discharge (TwoWayID, Flow, VelFather, VelSon, DecayTime, CoefCold, VelID, VelDT, &
    SonVolInFather, FatherVolume, STAT)
        !Arguments--------------------------------------------------------------
        integer,                            intent(IN   )     :: TwoWayID
        real,    dimension(:,:,:), pointer, intent(INOUT)     :: Flow
        real,    dimension(:,:,:), pointer, intent(IN   )     :: VelFather, SonVolInFather, FatherVolume
        real,    dimension(:,:,:), allocatable, intent(IN   ) :: VelSon
        real,    dimension(:,:  ), allocatable, intent(IN   ) :: DecayTime
        real                              , intent(IN   ) :: CoefCold, VelDT
        integer                           , intent(IN   ) :: VelID
        integer, optional                 , intent(OUT  ) :: STAT
        !Locals-----------------------------------------------------------------
        integer                                           :: line, i, j, k, STAT_CALL, STAT_, status, ready_
        real,    dimension(:,:,:), pointer                :: AreaU, AreaV
        type (T_TwoWay), pointer                          :: SonObj
        integer                                           :: MaxSize
        !Begin------------------------------------------------------------------
        
        call Ready(TwoWayID, ready_)

        if ((ready_ .EQ. IDLE_ERR_ ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            SonObj => FirstObjTwoWay
            do while (associated (SonObj))
                !if current object is a son domain then update discharge flow
                if (SonObj%Father%InstanceID == TwoWayID) then
                    
                    Call GetGeometryAreas(  GeometryID  = Me%ObjGeometry, AreaU = AreaU, AreaV = AreaV, &
                                            STAT        = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Offline_Upscaling_Discharge - failed get Son AreaU/V matrix'
                    
                    if (VelID == VelocityU_) then
                        
                        call Offline_DischargeFluxU(Flow                = SonObj%Father%DischargeCells%Flow,    &
                                                    DischargeConnection = SonObj%Father%DischargeCells%Z,       &
                                                    VelFather           = VelFather,                            &
                                                    VelSon              = VelSon,                               &
                                                    AreaU               = AreaU,                                &
                                                    DecayTime           = DecayTime,                            &
                                                    VelDT               = VelDT,                                &
                                                    SonVolInFather      = SonVolInFather,                       &
                                                    FatherVolume        = FatherVolume,                         &
                                                    CoefCold            = CoefCold)
                    else
                        
                        call Offline_DischargeFluxV(Flow                = SonObj%Father%DischargeCells%Flow,    &
                                                    DischargeConnection = SonObj%Father%DischargeCells%Z,       &
                                                    VelFather           = VelFather,                            &
                                                    VelSon              = VelSon,                               &
                                                    AreaV               = AreaV,                                &
                                                    DecayTime           = DecayTime,                            &
                                                    VelDT               = VelDT,                                &
                                                    SonVolInFather      = SonVolInFather,                       &
                                                    FatherVolume        = FatherVolume,                         &
                                                    CoefCold            = CoefCold)
                    endif
                    
                    !---------------------- Update Father discharge matrix --------------------------------------------
                    MaxSize = Size(SonObj%Father%DischargeCells%Flow)
                    do line = 1, MaxSize
                        i = SonObj%Father%DischargeCells%Z(line, 1)
                        j = SonObj%Father%DischargeCells%Z(line, 2)
                        k = SonObj%Father%DischargeCells%Z(line, 3)

                        Flow(i, j, k) = Flow(i, j, k) + SonObj%Father%DischargeCells%Flow(line)
                    enddo
                    !---------------------- Finish Update Father discharge matrix -------------------------------------
                    
                    call UnGetGeometry(Me%ObjGeometry, AreaU, STAT = status)
                    if (status /= SUCCESS_) stop 'Offline_Upscaling_Discharge-TwoWay-ERR01.'
                    call UnGetGeometry(Me%ObjGeometry, AreaV, STAT = status)
                    if (status /= SUCCESS_) stop 'Offline_Upscaling_Discharge-TwoWay-ERR02.'
                    
                endif
                
                SonObj => SonObj%Next
                
            enddo
            
            nullify (SonObj)
            
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        endif

        if (present(STAT)) STAT = STAT_
    end subroutine Offline_Upscaling_Discharge
    
    !--------------------------------------------------------------------------
    !>@author Joao Sobrinho +Atlantic
    !>@Brief
    !>Fills discharge flow matrix for an offline upscaling discharge
    !>@param[in] FatherID, PropAssimilation, Prop, PropVector, Flow, FlowVector, &
    !> DecayTime, CoefCold, DTProp, FatherVolume, dI, dJ, dK, Kmin, Kmax, AuxKmin, AuxKmax, CellID, nCells, VectorI, VectorJ, &
    !> VectorK, FoundDomain, STAT
    subroutine Offline_Upscaling_Discharge_WP (FatherID, PropAssimilation, Prop, PropVector, Flow, FlowVector, &
    DecayTime, CoefCold, DTProp, FatherVolume, dI, dJ, dK, Kmin, Kmax, AuxKmin, AuxKmax, CellID, nCells, VectorI, VectorJ, &
    VectorK, FoundDomain, STAT)
        !Arguments--------------------------------------------------------------
        integer,                            intent(IN )     :: FatherID, CellID, nCells
        real,    dimension(:,:,:), pointer, intent(IN )     :: Flow, Prop, PropAssimilation, FatherVolume
        real,    dimension(:    ), pointer, intent(INOUT)   :: FlowVector, PropVector
        real,    dimension(:,:  ), pointer, intent(IN)      :: DecayTime
        integer, dimension(:    ), pointer, intent(INOUT)   :: Kmin, Kmax
        integer, dimension(:    ), pointer, intent(IN )     :: VectorI, VectorJ, VectorK
        integer,                            intent(IN )     :: AuxKmin, AuxKmax
        real,                               intent(IN)      :: CoefCold, DTProp
        integer, dimension(:    ), pointer, intent(INOUT)   :: dI, dJ, dK
        logical                           , intent(OUT )    :: FoundDomain
        integer, optional                 , intent(OUT)     :: STAT
        !Locals-----------------------------------------------------------------
        integer                                             :: i, STAT_, ready_
        type (T_TwoWay), pointer                            :: SonObj
        integer                                             :: Aux, iSon, ifather, jfather, kfather
        real                                                :: Vol_Rat, TimeCoef
        integer, dimension(:,: ), pointer                   :: AuxConnections
        !Begin------------------------------------------------------------------
        STAT_ = UNKNOWN_
        
        call Ready(FatherID, ready_)

        if ((ready_ .EQ. IDLE_ERR_ ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            SonObj => FirstObjTwoWay
            do while (associated (SonObj))
                !if current object is a son domain then update discharge flow
                if (SonObj%Father%InstanceID == FatherID) then
                    AuxConnections => SonObj%Father%DischargeCells%Z
                    !call GetConnections(SonObj%ObjHorizontalGrid, Connections_Z = Connections_Z, STAT = STAT_CALL)
                    !if (STAT_CALL /= SUCCESS_) stop 'Offline_Upscaling_Discharge_WP - Failed to get Connections matrix'
                    !Cheks if current discharge is inside current upscaling domain
                    if (DischargeIsAssociated(AuxConnections, VectorI(1), VectorJ(1))) then
                        
                        Aux = CellID + 1
                        if (CellID == 0) Aux = 1

                        do i = Aux, nCells + CellID
                            iSon          = i-CellID
                            dI        (i) = VectorI(iSon)
                            dJ        (i) = VectorJ(iSon)
                            dK        (i) = VectorK(iSon)
                            Kmin      (i) = AuxKmin
                            Kmax      (i) = AuxKmax
                            FlowVector(i) = Flow(VectorI(iSon), VectorJ(iSon), VectorK(iSon))
                        
                            if (FlowVector(i) >= 0) then
                                ifather = VectorI(iSon) ; jfather = VectorJ(iSon) ; kfather = VectorK(iSon)
                                Vol_Rat = Me%TotSonIn(ifather,jfather,kfather) / FatherVolume(ifather ,jfather,kfather)
                                TimeCoef = (DTProp / DecayTime(ifather, jfather)) * Vol_Rat * CoefCold
                                PropVector(i) = (Prop            (ifather, jfather, kfather) + &
                                                 PropAssimilation(ifather, jfather, kfather) * TimeCoef) / (1 + TimeCoef)
                            else
                                PropVector(i) = Prop(VectorI(iSon), VectorJ(iSon), VectorK(iSon))
                            endif
                        enddo
                        FoundDomain = .true.
                    endif
                endif
                SonObj => SonObj%Next
            enddo
            nullify (SonObj)
            
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        endif

        if (present(STAT)) STAT = STAT_
        
    end subroutine Offline_Upscaling_Discharge_WP
    
    !--------------------------------------------------------------------------
    !>@author Joao Sobrinho +Atlantic
    !>@Brief
    !>Fills discharge flow matrix for an offline upscaling discharge
    !>@param[in] FatherID, PropAssimilation, Prop, PropVector, Flow, FlowVector, &
    !> FatherID, PropAssimilation, Prop, PropVector, Flow, FlowVector, &
    !> dI, dJ, dK, Kmin, Kmax, AuxKmin, AuxKmax, CellID, nCells, VectorI, VectorJ, VectorK, FoundDomain, STAT
    subroutine Offline_Upscaling_Discharge_WP_V2 (FatherID, PropAssimilation, Prop, PropVector, Flow, FlowVector, &
    dI, dJ, dK, Kmin, Kmax, AuxKmin, AuxKmax, CellID, nCells, VectorI, VectorJ, VectorK, FoundDomain, STAT)
        !Arguments--------------------------------------------------------------
        integer,                            intent(IN )     :: FatherID, CellID, nCells
        real,    dimension(:,:,:), pointer, intent(IN )     :: Flow, Prop, PropAssimilation
        real,    dimension(:    ), pointer, intent(INOUT)   :: FlowVector, PropVector
        integer, dimension(:    ), pointer, intent(INOUT)   :: Kmin, Kmax
        integer, dimension(:    ), pointer, intent(IN )     :: VectorI, VectorJ, VectorK
        integer,                            intent(IN )     :: AuxKmin, AuxKmax
        integer, dimension(:    ), pointer, intent(INOUT)   :: dI, dJ, dK
        logical                           , intent(OUT )    :: FoundDomain
        integer, optional                 , intent(OUT)     :: STAT
        !Locals-----------------------------------------------------------------
        integer                                             :: i, STAT_, ready_
        type (T_TwoWay), pointer                            :: SonObj
        integer                                             :: Aux, iSon
        integer, dimension(:,: ), pointer                   :: AuxConnections
        !Begin------------------------------------------------------------------
        STAT_ = UNKNOWN_
        
        call Ready(FatherID, ready_)

        if ((ready_ .EQ. IDLE_ERR_ ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            SonObj => FirstObjTwoWay
            do while (associated (SonObj))
                !if current object is a son domain then update discharge flow
                if (SonObj%Father%InstanceID == FatherID) then
                    AuxConnections => SonObj%Father%DischargeCells%Z
                    !Cheks if current discharge is inside current upscaling domain
                    if (DischargeIsAssociated(AuxConnections, VectorI(1), VectorJ(1))) then
                        Aux = CellID + 1
                        if (CellID == 0) Aux = 1

                        do i = Aux, nCells + CellID
                            iSon          = i-CellID
                            dI        (i) = VectorI(iSon)
                            dJ        (i) = VectorJ(iSon)
                            dK        (i) = VectorK(iSon)
                            Kmin      (i) = AuxKmin
                            Kmax      (i) = AuxKmax
                            FlowVector(i) = Flow(VectorI(iSon), VectorJ(iSon), VectorK(iSon))
                        
                            if (FlowVector(i) >= 0) then
                                PropVector(i) = PropAssimilation(VectorI(iSon), VectorJ(iSon), VectorK(iSon))
                            else
                                PropVector(i) = Prop(VectorI(iSon), VectorJ(iSon), VectorK(iSon))
                            endif
                        enddo
                        FoundDomain = .true.
                    endif
                endif
                SonObj => SonObj%Next
            enddo
            nullify (SonObj)
            
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        endif

        if (present(STAT)) STAT = STAT_
        
    end subroutine Offline_Upscaling_Discharge_WP_V2
    
    !--------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Fills domain matrix with user-provided upscaling decay time only in the overlapped area, excluding n cells
    !>from the border
    !>@param[in] SonID, Matrix2D, DecayValue, NumCellsToIgnore, STAT
    subroutine Fill_Upscaling_DecayTime(SonID, Matrix2D, DecayValue, NumCellsToIgnore, STAT)
        !Arguments--------------------------------------------------------------
        integer,                                     intent(IN)     :: SonID, NumCellsToIgnore
        real, dimension(:, :),    pointer,           intent(INOUT)  :: Matrix2D
        integer,                           optional, intent(OUT)    :: STAT
        real,                                        intent(IN)     :: DecayValue
        !Locals-----------------------------------------------------------------
        integer                                                     :: ready_, STAT_, STAT_CALL
        integer                                                     :: i, j, ILB, JLB, IUB, JUB, Flag, IFather, JFather
        !Begin------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(SonID, ready_)

        if ((ready_ .EQ. IDLE_ERR_ ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            !Son matrix bounds
            ILB = Me%WorkSize%ILB ; JLB = Me%WorkSize%JLB
            IUB = Me%WorkSize%IUB ; JUB = Me%WorkSize%JUB

            call PrepTwoWay(SonID, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Failed get Upscaling ExternalVars - Fill_Upscaling_DecayTime - TwoWay'

            if (.not. allocated(Me%IgnoreOBCells)) &
                call Compute_MatrixFilterOB (NumCellsToIgnore, Me%External_Var%BoundaryPoints2D)

            !A lot of repetitions because there are many son cells per father cell.. but its only for the construct
            !and the matrix is 2D so there is not much performance loss.
            do j = JLB, JUB
            do i = ILB, IUB
                Flag = Me%External_Var%Open3D(i, j, Me%WorkSize%KUB) + Me%IgnoreOBCells(i, j)
                if (Flag == 2) then
                    IFather = Me%External_Var%IZ(i, j)  ;  JFather = Me%External_Var%JZ(i, j)

                    Matrix2D(IFather, JFather) = DecayValue
                endif
            enddo
            enddo

            call UngetTwoWayExternal_Vars (SonID, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Failed unget Upscaling ExternalVars - Fill_Upscaling_DecayTime - TwoWay'

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        endif

        if (present(STAT)) STAT = STAT_

    end subroutine Fill_Upscaling_DecayTime

    !------------------------------------------------------------------------------

    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Ungets all external vars
    !>@param[in] SonID, CallerID, STAT
    subroutine UngetTwoWayExternal_Vars(SonID, CallerID, STAT)
        !Arguments--------------------------------------------------------------
        integer, intent(IN)                         :: SonID
        integer, optional, intent(IN)               :: CallerID
        integer, optional, intent(OUT)              :: STAT
        !Locals-----------------------------------------------------------------
        integer                                     :: ready_, status, STAT_, callerID_
        !Begin------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(SonID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            callerID_ = 0
            if (present(callerID)) callerID_ = callerID_

            if (callerID == 1) then
                !For future developments (when other modules call for twoway)
            endif
            !Unget son
            call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IV,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR01'
            call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%JV,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR10'
            call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IU,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR20'
            call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%JU,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR30'
            call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IZ,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR40'
            call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%JZ,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR50'
            
            call UngetGeometry(Me%ObjGeometry, Me%External_Var%KZ,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR60'
            call UnGetGeometry(Me%ObjGeometry, Me%External_Var%VolumeU,     STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR70'
            call UnGetGeometry(Me%ObjGeometry, Me%External_Var%VolumeV,     STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR80'
            call UnGetGeometry(Me%ObjGeometry, Me%External_Var%VolumeZ,     STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR90'
            call UnGetGeometry(Me%ObjGeometry, Me%External_Var%VolumeZ_2D,  STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR100'
            call UnGetGeometry(Me%ObjGeometry, Me%External_Var%AreaU, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR110.'
            call UnGetGeometry(Me%ObjGeometry, Me%External_Var%AreaV, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR120.'
            
            call UnGetMap(Me%ObjMap, Me%External_Var%Open3D,           STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR130'
            call UnGetMap(Me%ObjMap, Me%External_Var%ComputeFaces3D_U, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR140'
            call UnGetMap(Me%ObjMap, Me%External_Var%ComputeFaces3D_V, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR150'
            call UnGetMap(Me%ObjMap, Me%External_Var%WaterPoints3D, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR160'            

            call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%External_Var%BoundaryPoints2D, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR170'

            if (Me%Hydro%InterpolationMethod == 2) then
                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IWD_Connections_U, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR180'
                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IWD_Connections_V, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR190'
                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IWD_Connections_Z, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR200'
                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IWD_Distances_U, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR210'
                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IWD_Distances_V, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR220'
                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IWD_Distances_Z, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR230'
            endif

            !Unget father
            call UnGetMap(Me%Father%ObjMap, Me%Father%External_Var%Open3D,                STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR240'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%VolumeZ,     STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR250'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%VolumeU,     STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR260'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%VolumeV,     STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR270'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%VolumeZ_2D,  STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR280'
            call UnGetMap(Me%Father%ObjMap, Me%Father%External_Var%ComputeFaces3D_U,      STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR290'
            call UnGetMap(Me%Father%ObjMap, Me%Father%External_Var%ComputeFaces3D_V,      STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR300'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%AreaU,       STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR310.'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%AreaV,       STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR320.'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%KFloor_U,    STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR330.'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%KFloor_V,    STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR340.'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%KFloor_Z,    STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR350.'
            call UnGetHorizontalGrid(Me%Father%ObjHorizontalGrid, Me%Father%External_Var%AreaZ, STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR360.'

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        endif

        if (present(STAT)) STAT = STAT_

    end subroutine UngetTwoWayExternal_Vars
    !---------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillTwoWay(ObjTwoWayID, STAT)
        !Arguments---------------------------------------------------------------
        integer                             :: ObjTwoWayID
        integer, optional, intent(OUT)      :: STAT
        !External----------------------------------------------------------------
        integer                             :: ready_
        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers
        !------------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(ObjTwoWayID, ready_)

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mTwoWay_,  Me%InstanceID)

            if (nUsers == 0) then
                call DeallocateVariables

                nUsers = DeassociateInstance (mHORIZONTALGRID_,     Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillTwoWay - ModuleTwoWay - ERR01'
                nUsers = DeassociateInstance (mGEOMETRY_,           Me%ObjGeometry)
                if (nUsers == 0) stop 'KillTwoWay - ModuleTwoWay - ERR02'
                nUsers = DeassociateInstance (mMAP_,                Me%ObjMap)
                if (nUsers == 0) stop 'KillTwoWay - ModuleTwoWay - ERR03'
                nUsers = DeassociateInstance (mHORIZONTALMAP_,      Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillTwoWay - ModuleTwoWay - ERR04'

                !Deallocates Instance
                call DeallocateInstance ()
                ObjTwoWayID = 0
                STAT_      = SUCCESS_
            end if
        else
            STAT_ = ready_
        end if cd1
        if (present(STAT)) STAT = STAT_

    end subroutine KillTwoWay

    !------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Deallocates TwoWay module variables
    subroutine DeallocateVariables

        if (allocated(Me%Father%AuxMatrix)) deallocate(Me%Father%AuxMatrix)
        if (allocated(Me%Father%AuxMatrix2D)) deallocate(Me%Father%AuxMatrix2D)
        if (allocated(Me%Father%TotSonIn)) deallocate(Me%Father%TotSonIn)
        if (allocated(Me%Father%TotSonIn_2D)) deallocate(Me%Father%TotSonIn_2D)
        if (allocated(Me%IgnoreOBCells)) deallocate(Me%IgnoreOBCells)
        if (allocated(Me%Father%IWDNom)) deallocate(Me%Father%IWDNom)
        if (allocated(Me%Father%IWDDenom)) deallocate(Me%Father%IWDDenom)

    end subroutine DeallocateVariables

    !-------------------------------------------------------------------------

    subroutine DeallocateInstance ()
        type (T_TwoWay), pointer          :: AuxObjTwoWay
        type (T_TwoWay), pointer          :: PreviousObjTwoWay

        !Updates pointers
        if (Me%InstanceID == FirstObjTwoWay%InstanceID) then
            FirstObjTwoWay => FirstObjTwoWay%Next
        else
            PreviousObjTwoWay => FirstObjTwoWay
            AuxObjTwoWay      => FirstObjTwoWay%Next
            do while (AuxObjTwoWay%InstanceID /= Me%InstanceID)
                PreviousObjTwoWay => AuxObjTwoWay
                AuxObjTwoWay      => AuxObjTwoWay%Next
            enddo

            !Now update linked list
            PreviousObjTwoWay%Next => AuxObjTwoWay%Next
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

    subroutine Ready (ObjTwoWay_ID, ready_)
        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTwoWay_ID
        integer                                     :: ready_
        !----------------------------------------------------------------------
        nullify (Me)

cd1:    if (ObjTwoWay_ID > 0) then
            call LocateObjTwoWay (ObjTwoWay_ID)
            ready_ = VerifyReadLock (mTwoWay_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

    end subroutine Ready

    subroutine ReadyFather(TwoWayID, TwoWayFather, ready_)
        !Arguments-------------------------------------------------------------
        integer                                     :: TwoWayID
        type (T_TwoWay), pointer                    :: TwoWayFather
        integer                                     :: ready_
        !----------------------------------------------------------------------
        nullify (TwoWayFather)
cd1:    if (TwoWayID > 0) then
            TwoWayFather => FirstObjTwoWay
            do while (associated (TwoWayFather))
                if (TwoWayFather%InstanceID == TwoWayID) exit
                TwoWayFather => TwoWayFather%Next
            enddo

            if (.not. associated(TwoWayFather))                                             &
                stop 'ModuleTwoWay - ReadyFather - father (Me) object not associated'
            ready_ = VerifyReadLock (mTwoWay_, TwoWayFather%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1
    end subroutine ReadyFather

    !--------------------------------------------------------------------------

    subroutine LocateObjTwoWay (ObjTwoWayID)
        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTwoWayID
        !Local-----------------------------------------------------------------
        Me => FirstObjTwoWay
        do while (associated (Me))
            if (Me%InstanceID == ObjTwoWayID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleTwoWay - LocateObjTwoWay - ERR01'

    end subroutine LocateObjTwoWay
    !--------------------------------------------------------------------------

end module ModuleTwoWay
