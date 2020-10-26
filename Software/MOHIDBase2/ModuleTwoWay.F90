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
                                       GetGeometryKFloor
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, UngetHorizontalGrid, GetHorizontalGridSize, GetConnections, &
                                       UnGetConnections, ConstructP2C_IWD, ConstructP2C_Avrg, GetGridCellArea

    use ModuleHorizontalMap,    only : GetBoundaries, UnGetHorizontalMap
    use ModuleFunctions

    use ModuleMap,              only : GetComputeFaces3D, GetOpenPoints3D, GetWaterPoints3D, UnGetMap
    use ModuleStopWatch,        only : StartWatch, StopWatch

    use ModuleDischarges,        only : GetDischargesNumber, GetDischargeFlowDistribuiton, IsUpscaling

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
    !Modifier
    public  :: ModifyTwoWay
    private ::  ComputeAuxMatrixes
    private ::      ComputeSonVolInFather
    private ::          ComputeSonVolInFather_offline
    private ::  Nudging_average
    private ::  Nudging_IWD
    public  :: PrepTwoWay
    public  :: UpscalingVolumeVariation
    public  :: UngetTwoWayExternal_Vars
    public  :: UpscaleDischarge
    public  :: UpscaleDischarge_WP
    public  :: Offline_Upscaling_Discharge
    public  :: Offline_Upscaling_Discharge_WP

    !Destructor
    public  :: KillTwoWay
    private ::      DeAllocateInstance
    private ::      DeallocateVariables

    !Management
    private ::      Ready
    private ::      ReadyFather
    private ::          LocateObjTwoWay

    !Interfaces----------------------------------------------------------------

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
    end type T_Discharges

    private :: T_FatherDomain
    type       T_FatherDomain
        type (T_Size3D)                             :: Size, WorkSize
        type (T_Size2D)                             :: Size2D, WorkSize2D
        type (T_External)                           :: External_Var
        type (T_Discharges)                         :: DischargeCells
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
        type (T_Size3D)                             :: Size, WorkSize
        type (T_Size2D)                             :: Size2D, WorkSize2D
        type (T_FatherDomain)                       :: Father
        real, dimension (:, :, :), pointer          :: TotSonIn => null()
        real, dimension (:, :   ), pointer          :: TotSonIn_2D => null()
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
    !>@param[in] ModelName, ObjTwoWayID
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
        allocate (Me%IgnoreOBCells(ILB:IUB, JLB:JUB))
        Me%IgnoreOBCells(:,:) = 1

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
    !>@param[in] FatherTwoWayID, TwoWayID
    subroutine AllocateTwoWayAux(FatherTwoWayID, TwoWayID)

        !Arguments-------------------------------------------------------------
        integer                            :: FatherTwoWayID, TwoWayID
        !Local-----------------------------------------------------------------
        integer                            :: ready_, ready_father, ILB, IUB, JLB, JUB, KLB, KUB, STAT_CALL
        logical                            :: isIWD
        type (T_TwoWay), pointer           :: ObjFather
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
            else
                call ConstructP2C_IWD(Me%Father%ObjHorizontalGrid, Me%ObjHorizontalGrid)

                allocate(Me%Father%IWDNom  (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%Father%IWDDenom(ILB:IUB, JLB:JUB, KLB:KUB))
                Me%Father%IWDDenom(:,:,:) = 0.0
                Me%Father%IWDNom(:,:,:) = 0.0
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
    subroutine GetUpscalingDischarge(TwoWayID, UpscaleDischarge, STAT)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                 :: TwoWayID
        real(8),  dimension(:,:, :), pointer, intent(INOUT) :: UpscaleDischarge
        integer, optional                                   :: STAT
        !Local-----------------------------------------------------------------
        integer                                             :: ready_, STAT_, line, MaxSize, i, j, k
        !----------------------------------------------------------------------
        call Ready(TwoWayID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            MaxSize = Size(Me%DischargeCells%AuxFlow)
            do line = 1, MaxSize
                i = Me%DischargeCells%AuxConnections(line, 1)
                j = Me%DischargeCells%AuxConnections(line, 2)
                k = Me%DischargeCells%AuxConnections(line, 3)

                UpscaleDischarge(i, j, k) = Me%DischargeCells%AuxFlow(line)
            enddo

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

            call Read_UnLock(mTwoWay_, Me%InstanceID, "UnGetSonVolInFather")

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
    end subroutine UnGetSonVolInFather

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Updates father grid domain with son's results
    !>@param[in] SonID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, CallerID, VelocityID, TD, STAT
    subroutine ModifyTwoWay(SonID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, CallerID, VelocityID, TD, STAT)
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
        logical                                     :: Offline
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
            elseif (CallerID == mWATERPROPERTIES_) then
                ! for now interpolation method is set by the hydrodynamic module. The if is here for when
                ! each property can set its own interpolation method
                InterpolMethod = Me%Hydro%InterpolationMethod
                TimeDecay = TD
            elseif (CallerID == mField4D_) then
                InterpolMethod = 1
                Offline =  .true.
            endif

            !if it is a 3D matrix
            if (present(FatherMatrix)) then
                if (present(VelocityID))then
                    if (VelocityID == VelocityU_) then
                        !if 3D matrixes were sent. Even 2D domains allocate a 3D matrix (only one vertical layer)
                        !Type_U
                        call ComputeAuxMatrixes (Volume_3D = Me%External_Var%VolumeU, InterpolMethod = InterpolMethod,&
                                                 Ilink = Me%External_Var%IU, Jlink = Me%External_Var%JU, &
                                                 VelocityID = VelocityID, Offline = Offline)
                    else
                        !Type_V
                        call ComputeAuxMatrixes (Volume_3D = Me%External_Var%VolumeV, InterpolMethod = InterpolMethod,&
                                                 Ilink = Me%External_Var%IV, Jlink = Me%External_Var%JV, &
                                                 VelocityID = VelocityID, Offline = Offline)
                    endif
                else
                    !Type Z
                    call ComputeAuxMatrixes     (Volume_3D = Me%External_Var%VolumeZ, InterpolMethod = InterpolMethod,&
                                                 Ilink = Me%External_Var%IZ, Jlink = Me%External_Var%JZ, &
                                                 Offline = Offline)
                endif
            else
                !if a 2D matrix was sent (specific for waterLevel - at least for MohidWater).
                call ComputeAuxMatrixes (Volume_2D = Me%External_Var%VolumeZ_2D, InterpolMethod = InterpolMethod, &
                                         Ilink = Me%External_Var%IZ, Jlink = Me%External_Var%JZ, Offline = Offline)
            endif

            if (InterpolMethod == 1) then
                if (CallerID == mField4D_) then
                    call Nudging_average_offline (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D)
                else
                    call Nudging_average (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, TimeDecay)
                endif

            elseif (InterpolMethod == 2) then
                call Nudging_IWD (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, TimeDecay)
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
    !>@param[in] Volume_3D, Volume_2D, VelocityID, InterpolMethod, Ilink, Jlink
    subroutine ComputeAuxMatrixes(Volume_3D, Volume_2D, VelocityID, InterpolMethod, Ilink, Jlink, Offline)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                        :: interpolMethod
        real(8), dimension(:, :, :), pointer, optional, intent(IN) :: Volume_3D
        real(8), dimension(:, :),    pointer, optional, intent(IN) :: Volume_2D
        integer, dimension(:, :), pointer, intent(IN)              :: Ilink, Jlink
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
                        call ComputeSonVolInFather(Volume_3D = Volume_3D, Ilink = Ilink, Jlink = Jlink, &
                                                   SonComputeFaces = Me%External_Var%ComputeFaces3D_U, Offline=Offline)
                    else
                        call ComputeSonVolInFather(Volume_3D = Volume_3D, Ilink = Ilink, Jlink = Jlink, &
                                                   SonComputeFaces = Me%External_Var%ComputeFaces3D_V, Offline=Offline)
                    endif
                else
                    call ComputeSonVolInFather   (Volume_3D = Volume_3D, Ilink = Ilink, Jlink = Jlink, Offline=Offline)
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
    subroutine ComputeSonVolInFather (Volume_3D, Volume_2D, Ilink, Jlink, SonComputeFaces, Offline)
        !Arguments--------------------------------------------------------------------------------
        real(8), dimension(:, :, :), pointer, optional :: Volume_3D
        real(8), dimension(:, :),    pointer, optional :: Volume_2D
        integer, dimension(:, :), pointer              :: Ilink, Jlink
        integer, dimension(:, :, :), pointer, optional :: SonComputeFaces
        logical,                            intent(IN) :: Offline
        !Local variables--------------------------------------------------------------------------
        integer                                 :: i, j, k, ifather, jfather, kfather, k_difference, CHUNK, Flag
        !Begin------------------------------------------------------------------------------------

        if (Offline) call ComputeSonVolInFather_offline(Volume_3D, Volume_2D, Ilink, Jlink)

        if (present(Volume_3D) .and. .not. offline) then
            !account for different amount of vertical layers between domains.
            k_difference = Me%Father%WorkSize%KUB-Me%WorkSize%KUB
            CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
            !$OMP PARALLEL PRIVATE(i,j,k,Flag, ifather, jfather, kfather)
            if (present(SonComputeFaces))then
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    Flag = Me%External_Var%Open3D(i, j, k) + Me%IgnoreOBCells(i, j) + SonComputeFaces(i, j, k)
                    if (Flag == 3) then
                        ifather = ILink(i, j); jfather = JLink(i, j); kfather = k + k_difference
                        Me%Father%TotSonIn(ifather, jfather, kfather) = Me%Father%TotSonIn(ifather, jfather, kfather) &
                                                                        + Volume_3D(i, j, k)
                    endif
                enddo
                enddo
                enddo
                !$OMP END DO
            else
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    Flag = Me%External_Var%Open3D(i, j, k) + Me%IgnoreOBCells(i, j)
                    if (Flag == 2) then
                        ifather = ILink(i, j); jfather = JLink(i, j); kfather = k + k_difference
                        Me%Father%TotSonIn(ifather, jfather, kfather) = Me%Father%TotSonIn(ifather, jfather, kfather) &
                                                                        + Volume_3D(i, j, k)
                    endif
                enddo
                enddo
                enddo
                !$OMP END DO
            endif
            !$OMP END PARALLEL
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
    subroutine ComputeSonVolInFather_offline (Volume_3D, Volume_2D, Ilink, Jlink)
        !Arguments--------------------------------------------------------------------------------
        real(8), dimension(:, :, :), pointer, optional :: Volume_3D
        real(8), dimension(:, :),    pointer, optional :: Volume_2D
        integer, dimension(:, :), pointer              :: Ilink, Jlink
        !Local variables--------------------------------------------------------------------------
        integer                                        :: i, j, k, ifather, jfather, kfather, k_difference, KLB, KUB
        integer                                        :: CHUNK
        !Begin------------------------------------------------------------------------------------
        KLB = Me%WorkSize%KLB; KUB = Me%WorkSize%KUB

        if (present(Volume_3D)) then
            !account for different amount of vertical layers between domains.
            k_difference = Me%Father%WorkSize%KUB-Me%WorkSize%KUB
            CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)

            !$OMP PARALLEL PRIVATE(i,j,k, ifather, jfather, kfather)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                ifather = ILink(i, j); jfather = JLink(i, j); kfather = k + k_difference
                Me%Father%TotSonIn(ifather,jfather,kfather) = Me%Father%TotSonIn(ifather, jfather, kfather) &
                                                            + Volume_3D(i, j, k) * Me%External_Var%WaterPoints3D(i,j,k)
            enddo
            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

        elseif (present(Volume_2D)) then
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                ifather = ILink(i, j); jfather = JLink(i, j)
                Me%Father%TotSonIn_2D(ifather, jfather) = Me%Father%TotSonIn_2D(ifather, jfather) &
                                                        + Volume_2D(i, j) * Me%External_Var%WaterPoints3D(i, j, KUB)
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
                                    SonMask          = Me%External_Var%WaterPoints3D,       &
                                    SizeFather       = Me%Father%WorkSize,                  &
                                    SizeSon          = Me%WorkSize,                         &
                                    Ilink            = Me%External_Var%IZ,                  &
                                    Jlink            = Me%External_Var%JZ,                  &
                                    VolumeSon        = Me%External_Var%VolumeZ,             &
                                    TotSonIn         = GetPointer(Me%Father%TotSonIn))

            call ReadyFather(Me%Father%InstanceID, ObjFather, ready_father)
            if (.not. associated(ObjFather%TotSonIn)) ObjFather%TotSonIn => Me%Father%TotSonIn

        else
            call Upscaling_Avrg_WL (    FatherMatrix2D   = FatherMatrix2D,                      &
                                        SonMatrix2D      = SonMatrix2D,                         &
                                        SonMask          = Me%External_Var%WaterPoints3D,       &
                                        SizeSon          = Me%WorkSize,                         &
                                        Ilink            = Me%External_Var%IZ,                  &
                                        Jlink            = Me%External_Var%JZ,                  &
                                        VolumeSon2D      = Me%External_Var%VolumeZ_2D,          &
                                        TotSonIn2D       = GetPointer(Me%Father%TotSonIn_2D))

            call ReadyFather(Me%Father%InstanceID, ObjFather, ready_father)
            if (.not. associated(ObjFather%TotSonIn_2D)) ObjFather%TotSonIn_2D => Me%Father%TotSonIn_2D
        endif

    end subroutine Nudging_average_offline

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

            call GetDischargeFlowDistribuiton(FatherID, DischargeIDNumber = dis, nCells = nCells, &
                                                VectorI = VectorI, VectorJ = VectorJ, VectorK = VectorK, &
                                                kmin = AuxKmin, kmax = AuxKmax, STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'UpscaleDischarge_WP - failed to get dischargeflowdistribution'

            !using VectorI/J(1) because the value is the same for the entire vector (only the K value changes)
            if (DischargeIsAssociated (Connections_Z, VectorI(1), VectorJ(1))) then
                if (IsUpscaling(FatherID, dis)) then
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
                else
                    AuxCell = AuxCell + nCells
                endif
            endif
        enddo

    end subroutine UpscaleDischarge_WP

    !--------------------------------------------------------------------------
    !>@author Joao Sobrinho +Atlantic
    !>@Brief
    !>Fills discharge flow matrix for an offline upscaling discharge
    !>@param[in] TwoWayID, Flow, VelFather, VelSon, DecayTime, CoefCold, VelID, VelDT, STAT
    subroutine Offline_Upscaling_Discharge (TwoWayID, Flow, VelFather, VelSon, DecayTime, CoefCold, VelID, VelDT, STAT)
        !Arguments--------------------------------------------------------------
        integer,                            intent(IN   ) :: TwoWayID
        real,    dimension(:,:,:), pointer, intent(INOUT) :: Flow
        real,    dimension(:,:,:), pointer, intent(IN   ) :: VelFather, VelSon
        real,    dimension(:,:  ), pointer, intent(IN   ) :: DecayTime
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
                                                    CoefCold            = CoefCold)
                    else
                        
                        call Offline_DischargeFluxV(Flow                = SonObj%Father%DischargeCells%Flow,    &
                                                    DischargeConnection = SonObj%Father%DischargeCells%Z,       &
                                                    VelFather           = VelFather,                            &
                                                    VelSon              = VelSon,                               &
                                                    AreaV               = AreaV,                                &
                                                    DecayTime           = DecayTime,                            &
                                                    VelDT               = VelDT,                                &
                                                    CoefCold            = CoefCold)
                    endif
                    
                    !---------------------- Update Father discharge matrix --------------------------------------------
                    MaxSize = Size(SonObj%Father%DischargeCells%Flow)
                    do line = 1, MaxSize
                        i = SonObj%Father%DischargeCells%Z(line, 1)
                        j = SonObj%Father%DischargeCells%Z(line, 2)
                        k = SonObj%Father%DischargeCells%Z(line, 3)

                        Flow(i, j, k) = SonObj%Father%DischargeCells%Flow(line)
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
    !>@param[in] FatherID, PropAssimilation, Prop, PropVector, Flow, FlowVector, dI, dJ, dK, Kmin, Kmax, AuxKmin, AuxKmax, STAT
    subroutine Offline_Upscaling_Discharge_WP (FatherID, PropAssimilation, Prop, PropVector, Flow, FlowVector, dI, dJ, &
    dK, Kmin, Kmax, AuxKmin, AuxKmax, CellID, nCells, VectorI, VectorJ, VectorK, STAT)
        !Arguments--------------------------------------------------------------
        integer,                            intent(IN )     :: FatherID, CellID, nCells
        real,    dimension(:,:,:), pointer, intent(IN )     :: Flow, Prop, PropAssimilation
        real,    dimension(:    ), pointer, intent(INOUT)     :: FlowVector, PropVector
        integer, dimension(:    ), pointer, intent(INOUT)     :: Kmin, Kmax
        integer, dimension(:    ), pointer, intent(IN )     :: VectorI, VectorJ, VectorK
        integer,                            intent(IN )     :: AuxKmin, AuxKmax
        integer, dimension(:    ), pointer, intent(INOUT)     :: dI, dJ, dK
        integer, optional                 , intent(OUT)     :: STAT
        !Locals-----------------------------------------------------------------
        integer                                             :: i, STAT_CALL, STAT_, ready_
        type (T_TwoWay), pointer                            :: SonObj
        integer                                             :: Aux
        integer, dimension(:,: ), pointer                   :: Connections_Z
        !Begin------------------------------------------------------------------
        STAT_ = UNKNOWN_
        
        call Ready(FatherID, ready_)

        if ((ready_ .EQ. IDLE_ERR_ ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            SonObj => FirstObjTwoWay
            do while (associated (SonObj))
                !if current object is a son domain then update discharge flow
                if (SonObj%Father%InstanceID == FatherID) then
                    
                    call GetConnections(SonObj%ObjHorizontalGrid, Connections_Z = Connections_Z, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Modify_Upscaling_Discharges - Failed to get Connections matrix'
                    !Cheks if current discharge is inside current upscaling domain
                    if (DischargeIsAssociated(Connections_Z, VectorI(1), VectorJ(1))) then
                        
                        Aux = CellID + 1
                        if (CellID == 0) Aux = 1

                        do i = Aux, nCells + CellID
                            dI        (i) = VectorI(i)
                            dJ        (i) = VectorJ(i)
                            dK        (i) = VectorK(i)
                            Kmin      (i) = AuxKmin
                            Kmax      (i) = AuxKmax
                            FlowVector(i) = Flow(VectorI(i), VectorJ(i), VectorK(i))
                        
                            if (FlowVector(i) >= 0) then
                                PropVector(i) = PropAssimilation(VectorI(i), VectorJ(i), VectorK(i))
                            else
                                PropVector(i) = Prop(VectorI(i), VectorJ(i), VectorK(i))
                            endif
                        enddo
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
                Flag = Me%External_Var%WaterPoints3D(i, j, Me%WorkSize%KUB) + Me%IgnoreOBCells(i, j)
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
            call UnGetGeometry(Me%ObjGeometry, Me%External_Var%VolumeU,     STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR60'
            call UnGetGeometry(Me%ObjGeometry, Me%External_Var%VolumeV,     STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR70'
            call UnGetGeometry(Me%ObjGeometry, Me%External_Var%VolumeZ,     STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR80'
            call UnGetGeometry(Me%ObjGeometry, Me%External_Var%VolumeZ_2D,  STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR90'
            call UnGetMap(Me%ObjMap, Me%External_Var%Open3D,           STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR100'
            call UnGetMap(Me%ObjMap, Me%External_Var%ComputeFaces3D_U, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR110'
            call UnGetMap(Me%ObjMap, Me%External_Var%ComputeFaces3D_V, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR120'
            call UnGetGeometry(Me%ObjGeometry, Me%External_Var%AreaU, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR121.'
            call UnGetGeometry(Me%ObjGeometry, Me%External_Var%AreaV, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR122.'
            call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%External_Var%BoundaryPoints2D, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR123'
            call UnGetMap(Me%ObjMap, Me%External_Var%WaterPoints3D, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR124'

            if (Me%Hydro%InterpolationMethod == 2) then
                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IWD_Connections_U, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR130'
                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IWD_Connections_V, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR140'
                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IWD_Connections_Z, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR150'
                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IWD_Distances_U, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR160'
                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IWD_Distances_V, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR170'
                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%External_Var%IWD_Distances_Z, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR180'
            endif

            !Unget father
            call UnGetMap(Me%Father%ObjMap, Me%Father%External_Var%Open3D,                STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR190'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%VolumeZ,     STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR200'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%VolumeU,     STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR210'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%VolumeV,     STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR220'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%VolumeZ_2D,  STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR230'
            call UnGetMap(Me%Father%ObjMap, Me%Father%External_Var%ComputeFaces3D_U,      STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR240'
            call UnGetMap(Me%Father%ObjMap, Me%Father%External_Var%ComputeFaces3D_V,      STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR250'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%AreaU,       STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR260.'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%AreaV,       STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR270.'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%KFloor_U,    STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR275.'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%KFloor_V,    STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR280.'
            call UnGetGeometry(Me%Father%ObjGeometry, Me%Father%External_Var%KFloor_Z,    STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR285.'
            call UnGetHorizontalGrid(Me%Father%ObjHorizontalGrid, Me%Father%External_Var%AreaZ, STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-TwoWay-ERR290.'

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
