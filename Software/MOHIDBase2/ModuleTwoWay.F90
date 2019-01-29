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
!>  DESCRIPTION : Module that computes twoway operations for MohidWater and MohidLand
!
!> @author
!> Joao Sobrinho
!------------------------------------------------------------------------------


Module ModuleTwoWay

    use ModuleGlobalData
    use ModuleGeometry,         only : GetGeometryVolumes, UnGetGeometry, GetGeometrySize
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, UngetHorizontalGrid, GetHorizontalGridSize, GetConnections, &
                                       UnGetConnections, ConstructP2C_IWD, ConstructP2C_Avrg
    
    use ModuleHorizontalMap,    only : GetBoundaries, UnGetHorizontalMap
    use ModuleFunctions
    use ModuleMap,              only : GetComputeFaces3D, GetOpenPoints3D, UnGetMap
    use ModuleStopWatch,        only : StartWatch, StopWatch
    !use ModuleTwoWayDischarges
    
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructTwoWay
    private ::  AllocateInstance
    
    
    public  :: ConstructTwoWayHydrodynamic
    private ::  Compute_MatrixFilterOB
    public  :: AllocateTwoWayAux
    public  :: Construct_Upscaling_Discharges

    !Selector                
    
    !Modifier
    public  :: ModifyTwoWay
    private ::  ComputeAuxMatrixes
    private ::    ComputeSonVolInFather
    private ::  Nudging_average
    private ::  Nudging_IWD
    public  :: PrepTwoWay
    public  :: UngetTwoWayExternal_Vars

    !Destructor
    public  :: KillTwoWay                                                     
    private ::      DeAllocateInstance
    private ::      DeallocateVariables

    !Management
    private ::      Ready
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
        integer, dimension(:, :, :), pointer        :: Open3D           => null()
        integer, dimension(:, :, :), pointer        :: WaterPoints3D    => null()
        integer, dimension(:, :, :), pointer        :: ComputeFaces3D_U => null()
        integer, dimension(:, :, :), pointer        :: ComputeFaces3D_V => null()
        integer, dimension(:, :   ), pointer        :: BoundaryPoints2D => null()
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
    
    private :: T_FatherDomain
    type       T_FatherDomain
        type (T_Size3D)                             :: Size, WorkSize
        type (T_Size2D)                             :: Size2D, WorkSize2D
        type (T_External)                           :: External_Var
        integer                                     :: InstanceID
        real, dimension (:, :, :), allocatable      :: TotSonVolIn
        real, dimension (:, :   ), allocatable      :: TotSonVolIn_2D
        real, dimension (:, :, :), allocatable      :: AuxMatrix
        real, dimension (:, :   ), allocatable      :: AuxMatrix2D
        real, dimension (:, :, :), allocatable      :: IWDNom
        real, dimension (:, :, :), allocatable      :: IWDDenom
    end type T_FatherDomain

    private :: T_TwoWay
    type       T_TwoWay
        integer                                     :: InstanceID
        character(PathLength)                       :: ModelName
        real(8), dimension(:, :, :),  pointer       :: Matrix
        integer, dimension(:, :   ),  allocatable   :: IgnoreOBCells

        type(T_External)                            :: External_Var
        type (T_Hydro)                              :: Hydro
        type (T_Size3D)                             :: Size, WorkSize
        type (T_Size2D)                             :: Size2D, WorkSize2D
        type (T_FatherDomain)                       :: Father
        type(T_TwoWay), pointer                     :: Next
        
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
    subroutine ConstructTwoWay(ModelName,               &
                               TwoWayID,                &
                               HorizontalGridID,        & 
                               GeometryID,              &
                               HorizontalMapID,         &
                               MapID,                   &
                               STAT)

        !Arguments---------------------------------------------------------------
        character(Len=*)                                :: ModelName
        integer                                         :: TwoWayID, HorizontalGridID
        integer                                         :: GeometryID, HorizontalMapID, MapID
        integer, optional                               :: STAT
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

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            Me%ModelName = ModelName
            
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjMap            = AssociateInstance (mMAP_,            MapID           )
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            !Returns ID
            TwoWayID          = Me%InstanceID
            
            call GetGeometrySize (GeometryID = Me%ObjGeometry, &
                                  Size       = Me%Size,        &
                                  WorkSize   = Me%WorkSize,    &
                                  STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - ConstructTwoWay - ERR01'
            
           call GetHorizontalGridSize (HorizontalGridID = Me%ObjHorizontalGrid, &
                                       Size             = Me%Size2D,            &
                                       WorkSize         = Me%WorkSize2D,        &
                                       STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - ConstructTwoWay - ERR02'

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleTwoWay - ConstructTwoWay - ERR01' 

        end if cd0

        if (present(STAT)) STAT = STAT_

    end subroutine ConstructTwoWay
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------                                            
        !Local-----------------------------------------------------------------
        type (T_TwoWay), pointer                         :: NewObjTwoWay
        type (T_TwoWay), pointer                         :: PreviousObjTwoWay


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
    !>@param[in] TwoWayID, TimeDecay, IntMethod, VelDT, DT, IgnoreOBNumCells
    subroutine ConstructTwoWayHydrodynamic (TwoWayID, TimeDecay, IntMethod, VelDT, DT, IgnoreOBNumCells, IWDn, STAT)
    
        !Arguments------------------------------------------------------------
        integer                                     :: TwoWayID, & ! ID
                                                       IntMethod !Method for grid interpolation from child to father grid
                                                              ! 1 - Volume Averaged; 2 - Inverse Weighted Distance
        integer                                     :: IgnoreOBNumCells ! number of Lines and columns ignored
                                                                      ! counting from an open boundary
        real                                        :: TimeDecay ! Decay factor in seconds, in the nudging equation
        real                                        :: VelDT, DT
        integer, optional, intent(OUT)              :: STAT 
        integer, optional                           :: IWDn
        !Local----------------------------------------------------------------
        integer                                     :: STAT_CALL, ready_, STAT_
        !---------------------------------------------------------------------

        call Ready(TwoWayID, ready_)    

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            Me%Hydro%InterpolationMethod = IntMethod
            Me%Hydro%TimeDecay           = TimeDecay
            Me%Hydro%VelDT               = VelDT
            Me%Hydro%DT                  = DT
            
            if (IntMethod == 2)then
                Me%Hydro%IWDn            = IWDn
            endif
            
            call GetBoundaries(HorizontalMapID     = TwoWayID,                              &
                               BoundaryPoints2D    = Me%External_Var%BoundaryPoints2D,      &
                               STAT                = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop "ModuleTwoWay - ConstructTwoWayHydrodynamic - ERR01"
            
            call Compute_MatrixFilterOB (IgnoreOBNumCells)
            
            call UnGetHorizontalMap(TwoWayID, Me%External_Var%BoundaryPoints2D, STAT = STAT_CALL)
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
    !>@param[in] IgnoreOBNumCells       
    subroutine Compute_MatrixFilterOB (IgnoreOBNumCells)
        !Arguments-------------------------------------------------------------
        integer                                     :: IgnoreOBNumCells, AuxIgnoreOBNumCells
        !Locals ---------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB, i, j
        !Begin ----------------------------------------------------------------
        ILB = Me%WorkSize2D%ILB 
        IUB = Me%WorkSize2D%IUB 
        JLB = Me%WorkSize2D%JLB 
        JUB = Me%WorkSize2D%JUB
        
        allocate (Me%IgnoreOBCells(ILB:IUB, JLB:JUB))
        call SetMatrixValue (GetPointer(Me%IgnoreOBCells), Me%WorkSize2D, 1)
        AuxIgnoreOBNumCells = IgnoreOBNumCells
        !compute south border
        do j = JLB, JUB
        do i = ILB, AuxIgnoreOBNumCells
            Me%IgnoreOBCells(i, j) = 0
        enddo
        enddo
        
        !compute North border
        if (IgnoreOBNumCells == IUB)then
            AuxIgnoreOBNumCells = AuxIgnoreOBNumCells - 1
        endif
        
        do j = JLB, JUB
        do i = IUB - AuxIgnoreOBNumCells, IUB
            Me%IgnoreOBCells(i, j) = 0
        enddo
        enddo
        AuxIgnoreOBNumCells = IgnoreOBNumCells
        !compute west border
        do j = JLB, AuxIgnoreOBNumCells
        do i = ILB, IUB
            Me%IgnoreOBCells(i, j) = 0
        enddo
        enddo
        !compute east border
        if (IgnoreOBNumCells == JUB)then
            AuxIgnoreOBNumCells = AuxIgnoreOBNumCells - 1
        endif        
        do j = JUB - AuxIgnoreOBNumCells, JUB
        do i = ILB, IUB
            Me%IgnoreOBCells(i, j) = 0
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
        integer                            :: ready_, ILB, IUB, JLB, JUB, KLB, KUB, STAT_CALL
        logical                            :: isIWD
        !----------------------------------------------------------------------
            
        call Ready (TwoWayID, ready_)
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            isIWD = .false.
            
            Me%Father%InstanceID = FatherTwoWayID
            
            call GetGeometrySize (FatherTwoWayID,                &
                                  Size     = Me%Father%Size,     &
                                  WorkSize = Me%Father%WorkSize, &
                                  STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                  &
                stop 'ModuleTwoWay - AllocateTwoWayAux - ERR10'
            
           call GetHorizontalGridSize (HorizontalGridID = FatherTwoWayID, &
                                       Size             = Me%Father%Size2D,            &
                                       WorkSize         = Me%Father%WorkSize2D,        &
                                       STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - AllocateTwoWayAux - ERR20'            
            
            ILB = Me%Father%WorkSize%ILB 
            IUB = Me%Father%WorkSize%IUB 
            JLB = Me%Father%WorkSize%JLB 
            JUB = Me%Father%WorkSize%JUB 
            KLB = Me%Father%WorkSize%KLB 
            KUB = Me%Father%WorkSize%KUB
            
            if (Me%Hydro%InterpolationMethod == 1)then
                allocate(Me%Father%TotSonVolIn(ILB:IUB, JLB:JUB, KLB:KUB))
                Me%Father%TotSonVolIn(:,:,:) = 0.0
            
                allocate(Me%Father%TotSonVolIn_2D(ILB:IUB, JLB:JUB))
                Me%Father%TotSonVolIn_2D(:,:) = 0.0
                
                allocate(Me%Father%AuxMatrix(ILB:IUB, JLB:JUB, KLB:KUB))
                Me%Father%AuxMatrix(:,:,:) = 0.0
                
                allocate(Me%Father%AuxMatrix2D(ILB:IUB, JLB:JUB))
                Me%Father%AuxMatrix2D(:,:) = 0.0
                
                !construct connection matrix of father-son for use in two-way discharges
                call ConstructP2C_Avrg(FatherTwoWayID, TwoWayID)
            else
                call ConstructP2C_IWD(FatherTwoWayID, TwoWayID)
                
                allocate(Me%Father%IWDNom(ILB:IUB, JLB:JUB, KLB:KUB))
                Me%Father%IWDNom(:,:,:) = 0.0
                allocate(Me%Father%IWDDenom(ILB:IUB, JLB:JUB, KLB:KUB))
                Me%Father%IWDDenom(:,:,:) = 0.0
            endif
        else
            stop 'ModuleTwoWay - AllocateTwoWayAux - ERR30'               
        endif

    end subroutine AllocateTwoWayAux
    
    ! ------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Detects and builds matrix with two-way discharges  
    !>@param[in] FatherTwoWayID, TwoWayID    
    subroutine Construct_Upscaling_Discharges(FatherTwoWayID, TwoWayID)
        !Arguments-------------------------------------------------------------
        integer                            :: FatherTwoWayID, TwoWayID
        !Local-----------------------------------------------------------------
        integer                            :: ready_, ILB, IUB, JLB, JUB, KLB, KUB, STAT_CALL
        !----------------------------------------------------------------------
        !call Ready (TwoWayID, ready_)
        !
        !if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
        !    (ready_ .EQ. READ_LOCK_ERR_)) then
            
            !Present = .false.
            !call GetConnections(TwoWayID, Me%External_Var%Connections_Z, STAT             = STAT_CALL)
            !if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Upscaling_Discharges - ModuleTwoWay - ERR10'
            !
            !call GetWaterPoints3D(TwoWayID,       Me%External_Var%WaterPoints3D,        STAT = STAT_CALL)
            !if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Upscaling_Discharges - ModuleTwoWay - ERR20'
            !
            !call GetWaterPoints3D(FatherTwoWayID, Me%Father%External_Var%WaterPoints3D, STAT = STAT_CALL)
            !if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Upscaling_Discharges - ModuleTwoWay - ERR30'          
            
            
            
            !call SearchDischargeFace(Me%External_Var%Connections_Z, Me%External_Var%WaterPoints3D, &
            !                         Me%Father%External_Var%WaterPoints3D, Me%Size2D, Me%Father%Size2D, &
            !                         SonLandPoints2D, FatherLandPoints2D, Present)
            !
            !
            !call UnGetHorizontalGrid(TwoWayID, Me%External_Var%Connections_Z, STAT_CALL)
            !if (STAT_CALL /= SUCCESS_) stop 'Construct_Upscaling_Discharges - ModuleTwoWay - ERR60'
            !
            !call UnGetMap(TwoWayID, Me%External_Var%WaterPoints3D, STAT = STAT_CALL)
            !if (status /= SUCCESS_) stop 'Construct_Upscaling_Discharges - ModuleTwoWay - ERR70'
            !call UnGetMap(TwoWayID, Me%Father%External_Var%WaterPoints3D, STAT = STAT_CALL)
            !if (status /= SUCCESS_) stop 'Construct_Upscaling_Discharges - ModuleTwoWay - ERR80'  
            
        !    
        !    call deallocateConnections_Z ! Adicionar isto no horizontal grid( a matrix já nao deve ser necessaria
        !                                 ! depois da verificaçao das descargas de quantidade de movimento
        !    
        !else  
        !    stop 'Construct_Upscaling_Discharges -ModuleTwoWay -  ERR90'               
        !endif
            
        
    end subroutine Construct_Upscaling_Discharges
    
    !-------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Updates father grid domain with son's results
    !>@param[in] SonID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, CallerID, TD, STAT
    subroutine ModifyTwoWay(SonID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, CallerID, VelocityID, TD, &
                            STAT)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: SonID, CallerID
        integer, optional                           :: VelocityID
        real,    optional                           :: TD !TimeDecay for twoway
        integer, optional, intent(OUT)              :: STAT
        real, dimension(:, :, :), pointer, optional :: FatherMatrix, SonMatrix
        real, dimension(:, :),    pointer, optional :: FatherMatrix2D, SonMatrix2D
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, InterpolMethod
        real                                        :: LocalTimeDecay

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SonID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (MonitorPerformance) call StartWatch ("ModuleTwoWay", "ModifyTwoWay")
            
            if (CallerID == mHydrodynamic_) then
                LocalTimeDecay = Me%Hydro%TimeDecay
                InterpolMethod = Me%Hydro%InterpolationMethod
            elseif (CallerID == mWATERPROPERTIES_) then
                ! for now interpolation method is set by the hydrodynamic module. The if is here for when
                ! each property can set its own interpolation method
                InterpolMethod = Me%Hydro%InterpolationMethod
                LocalTimeDecay = TD
            endif
            
            !if it is a 3D matrix
            if (present(FatherMatrix)) then
                if (present(VelocityID))then
                    if (VelocityID == VelocityU_) then
                        !if 3D matrixes were sent. Even 2D domains allocate a 3D matrix (only one vertical layer)
                        !Type_U
                        call ComputeAuxMatrixes (Volume_3D        = Me%External_Var%VolumeU,         &
                                                 InterpolMethod   = InterpolMethod,    &
                                                 Ilink            = Me%External_Var%IU, &
                                                 Jlink            = Me%External_Var%JU, &
                                                 VelocityID       = VelocityID)
                        
                    else
                        !Type_V
                        call ComputeAuxMatrixes (Volume_3D        = Me%External_Var%VolumeV,         &
                                                 InterpolMethod   = InterpolMethod,    &
                                                 Ilink            = Me%External_Var%IV, &
                                                 Jlink            = Me%External_Var%JV, &
                                                 VelocityID       = VelocityID)
                    endif
                else
                    !Type Z
                    call ComputeAuxMatrixes (Volume_3D        = Me%External_Var%VolumeZ,         &
                                             InterpolMethod   = InterpolMethod,    &
                                             Ilink            = Me%External_Var%IZ, &
                                             Jlink            = Me%External_Var%JZ)
                endif
                
            else
                !if a 2D matrix was sent (specific for waterLevel - at least for MohidWater).
                call ComputeAuxMatrixes (Volume_2D      = Me%External_Var%VolumeZ_2D,       &
                                         InterpolMethod   = InterpolMethod,    &
                                         Ilink            = Me%External_Var%IZ, &
                                         Jlink            = Me%External_Var%JZ)
            endif

            if (InterpolMethod == 1) then
                call Nudging_average (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, LocalTimeDecay)
            elseif (InterpolMethod == 2) then
                call Nudging_IWD (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, LocalTimeDecay)
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
    !>@param[in] SonID, FatherID, CallerID, STAT   
    subroutine PrepTwoWay (SonID, FatherID, CallerID, STAT)
    
        !Arguments--------------------------------------------------------------
        integer, intent(IN)                         :: SonID, FatherID, CallerID
        integer, optional                           :: STAT
        !Locals-----------------------------------------------------------------
        integer                                     :: STAT_CALL, ready_, STAT_
        !Begin------------------------------------------------------------------
        STAT_ = UNKNOWN_      
        
        call Ready(SonID, ready_)
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            if (callerID == 1) then
                !For future developments (when other modules call for twoway)
            endif
            
            call GetHorizontalGrid (HorizontalGridID = SonID,                         &
                                    ILinkV            = Me%External_Var%IV,                &
                                    JLinkV            = Me%External_Var%JV,                &
                                    ILinkU            = Me%External_Var%IU,                &
                                    JLinkU            = Me%External_Var%JU,                &
                                    ILinkZ            = Me%External_Var%IZ,                &
                                    JLinkZ            = Me%External_Var%JZ,                &
                                    STAT              = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - PrepTwoWay - ERR01'
            
            call GetGeometryVolumes(GeometryID     = SonID,                            &
                                    VolumeU        = Me%External_Var%VolumeU,           &
                                    VolumeV        = Me%External_Var%VolumeV,           &
                                    VolumeZ        = Me%External_Var%VolumeZ,           &
                                    VolumeZ_2D     = Me%External_Var%VolumeZ_2D,        &
                                    STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - PrepTwoWay - ERR02'
            
            call GetOpenPoints3D   (Map_ID         = SonID,                            &
                                    OpenPoints3D   = Me%External_Var%Open3D,            &
                                    STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - PrepTwoWay - ERR03'
            
            call GetComputeFaces3D(Map_ID          = SonID,                           &
                                   ComputeFacesU3D = Me%External_Var%ComputeFaces3D_U, &
                                   ComputeFacesV3D = Me%External_Var%ComputeFaces3D_V, &
                                   STAT            = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - PrepTwoWay - ERR04'
            
            if (Me%Hydro%InterpolationMethod == 2) then
                
                call GetConnections (HorizontalGridID       = SonID,   &
                                     Connections_U          = Me%External_Var%IWD_Connections_U, &
                                     IWD_Distances_U        = Me%External_Var%IWD_Distances_U,   &
                                     Connections_V          = Me%External_Var%IWD_Connections_V, &
                                     IWD_Distances_V        = Me%External_Var%IWD_Distances_V,   &
                                     Connections_Z          = Me%External_Var%IWD_Connections_Z, &
                                     IWD_Distances_Z        = Me%External_Var%IWD_Distances_Z,   &
                                     IWD_Nodes_Z            = Me%External_Var%IWD_Nodes_Z,       &
                                     IWD_Nodes_U            = Me%External_Var%IWD_Nodes_U,       &
                                     IWD_Nodes_V            = Me%External_Var%IWD_Nodes_V,       &
                                     STAT                   = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - PrepTwoWay - ERR05'

            endif                           

            
            call GetOpenPoints3D   (Map_ID         = FatherID,                         &
                                    OpenPoints3D   = Me%Father%External_Var%Open3D,     &
                                    STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - PrepTwoWay - ERR06'
        
            call GetGeometryVolumes(GeometryID     = FatherID,                         &
                                    VolumeU        = Me%Father%External_Var%VolumeU,    &
                                    VolumeV        = Me%Father%External_Var%VolumeV,    &
                                    VolumeZ        = Me%Father%External_Var%VolumeZ,    &
                                    VolumeZ_2D     = Me%Father%External_Var%VolumeZ_2D, &
                                    STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - PrepTwoWay - ERR07'
            
            call GetComputeFaces3D(Map_ID          = FatherID,                               &
                                   ComputeFacesU3D = Me%Father%External_Var%ComputeFaces3D_U, &
                                   ComputeFacesV3D = Me%Father%External_Var%ComputeFaces3D_V, &
                                   STAT            = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - PrepTwoWay - ERR08'
            
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        endif
        
        if (present(STAT)) STAT = STAT_
        
    end subroutine PrepTwoWay
    
    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Computes auxiliar matrixes for the feedback
    !>@param[in] Volume_3D, Volume_2D, VelocityID, InterpolMethod   
    subroutine ComputeAuxMatrixes(Volume_3D, Volume_2D, VelocityID, InterpolMethod, Ilink, Jlink)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                            :: interpolMethod
        real(8), dimension(:, :, :), pointer, optional :: Volume_3D
        real(8), dimension(:, :),    pointer, optional :: Volume_2D
        integer, dimension(:, :), pointer              :: Ilink, Jlink
        integer, optional                              :: VelocityID
        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------    
        if (present(Volume_3D)) then
            !Goes for 3D
            if     (interpolMethod == 1)then
                ! set the matrix to 0.001 so it can be divided without giving an error
                call SetMatrixValue (GetPointer(Me%Father%TotSonVolIn), Me%Father%WorkSize, 0.001)
                call SetMatrixValue (GetPointer(Me%Father%AuxMatrix), Me%Father%WorkSize, 0.0)                  
                ! Volume Weighted average
                if (present(VelocityID))then
                    if (VelocityID == VelocityU_)then
                        call ComputeSonVolInFather(Volume_3D      = Volume_3D,                     &
                                                   Ilink          = Ilink,                         &
                                                   Jlink          = Jlink,                         &
                                                   SonComputeFaces= Me%External_Var%ComputeFaces3D_U)
                    else
                        call ComputeSonVolInFather(Volume_3D      = Volume_3D,                     &
                                                   Ilink          = Ilink,                         &
                                                   Jlink          = Jlink,                         &
                                                   SonComputeFaces= Me%External_Var%ComputeFaces3D_V)
                    endif
                    
                else
                    
                    call ComputeSonVolInFather   (Volume_3D      = Volume_3D,      &
                                                  Ilink          = Ilink,          &
                                                  Jlink          = Jlink)
                endif
                
            else
                call SetMatrixValue (GetPointer(Me%Father%IWDNom), Me%Father%WorkSize, 0.0)
                call SetMatrixValue (GetPointer(Me%Father%IWDDenom), Me%Father%WorkSize, 0.0)
            endif
        else
        !Goes for 2D             
            if     (interpolMethod == 1)then
                call SetMatrixValue (GetPointer(Me%Father%AuxMatrix2D), Me%Father%WorkSize2D, 0.0)
                ! set the matrix to 0.001 so it can be divided without giving an error
                call SetMatrixValue (GetPointer(Me%Father%TotSonVolIn_2D), Me%Father%WorkSize2D, 0.001)
                
                ! Volume Weighted average
                call ComputeSonVolInFather   (Volume_2D    = Volume_2D,     &
                                              Ilink          = Ilink,           &
                                              Jlink          = Jlink)
            else
                call SetMatrixValue (GetPointer(Me%Father%IWDNom), Me%Father%WorkSize, 0.0)
                call SetMatrixValue (GetPointer(Me%Father%IWDDenom), Me%Father%WorkSize, 0.0)
            endif                
        endif
          
    end subroutine ComputeAuxMatrixes
    
    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Computes Son cells volume inside each father cell
    !>@param[in] SonMatrix, SonMatrix2D, Ilink, Jlink, SonComputeFaces
    subroutine ComputeSonVolInFather (Volume_3D, Volume_2D, Ilink, Jlink, SonComputeFaces)
    
        !Arguments--------------------------------------------------------------------------------
        real(8), dimension(:, :, :), pointer, optional :: Volume_3D
        real(8), dimension(:, :),    pointer, optional :: Volume_2D
        integer, dimension(:, :), pointer              :: Ilink, Jlink
        integer, dimension(:, :, :), pointer, optional :: SonComputeFaces
        !Local variables--------------------------------------------------------------------------
        integer                                 :: i, j, k, NThreads, OMPmethod, CHUNK
        !Begin------------------------------------------------------------------------------------
        OMPmethod = 2
        NThreads = openmp_num_threads
        if (NThreads > 1) then
            if (NThreads - Me%WorkSize%KUB > 1) then
                OMPmethod = 1
            endif
        endif
        
        if (present(Volume_3D)) then
            if (present(SonComputeFaces))then
                if (OMPmethod == 2) then
                    CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB, NThreads)
                    !$OMP PARALLEL PRIVATE(i,j,k)
                    !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                    do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                            Me%Father%TotSonVolIn(ILink(i, j), JLink(i, j), k) =                             &
                            Me%Father%TotSonVolIn(ILink(i, j), JLink(i, j), k) + Volume_3D(i, j, k) *        &
                            Me%External_Var%Open3D(i, j, k) * Me%IgnoreOBCells(i, j) * SonComputeFaces(i, j, k)
                    enddo        
                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL
                else
                    CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%KUB, NThreads)
                    !$OMP PARALLEL PRIVATE(i,j,k)
                    do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                    !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                            Me%Father%TotSonVolIn(ILink(i, j), JLink(i, j), k) =                             &
                            Me%Father%TotSonVolIn(ILink(i, j), JLink(i, j), k) + Volume_3D(i, j, k) *        &
                            Me%External_Var%Open3D(i, j, k) * Me%IgnoreOBCells(i, j) * SonComputeFaces(i, j, k)
                    enddo        
                    enddo
                    !$OMP END DO
                    enddo
                    !$OMP END PARALLEL                    
                endif
                
            else
                if (OMPmethod == 2) then
                    CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB, NThreads)
                    !$OMP PARALLEL PRIVATE(i,j,k)
                    !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                    do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                            Me%Father%TotSonVolIn(ILink(i, j), JLink(i, j), k) =                      &
                            Me%Father%TotSonVolIn(ILink(i, j), JLink(i, j), k) + Volume_3D(i, j, k) * &
                            Me%External_Var%Open3D(i, j, k) * Me%IgnoreOBCells(i, j)
                    enddo        
                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL
                else
                    CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB, NThreads)
                    !$OMP PARALLEL PRIVATE(i,j,k)
                    do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                    !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                            Me%Father%TotSonVolIn(ILink(i, j), JLink(i, j), k) =                      &
                            Me%Father%TotSonVolIn(ILink(i, j), JLink(i, j), k) + Volume_3D(i, j, k) * &
                            Me%External_Var%Open3D(i, j, k) * Me%IgnoreOBCells(i, j)
                    enddo        
                    enddo
                    !$OMP END DO
                    enddo
                    !$OMP END PARALLEL
                endif
                
            endif
            
        else
            CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB, NThreads)
            !$OMP PARALLEL PRIVATE(i,j,k)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                Me%Father%TotSonVolIn_2D(ILink(i, j), JLink(i, j)) = &
                Me%Father%TotSonVolIn_2D(ILink(i, j), JLink(i, j)) + Volume_2D(i, j) * &
                Me%External_Var%Open3D(i, j, Me%WorkSize%KUB) * Me%IgnoreOBCells(i, j)
            enddo        
            enddo            
            !$OMP END DO
            !$OMP END PARALLEL
        endif
          
    end subroutine ComputeSonVolInFather
    
    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Calls Feedback routines present in mondule functions, based on the type of input matrixes
    !>@param[in] FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, LocalTimeDecay
    subroutine Nudging_average (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, LocalTimeDecay)
        !Arguments-------------------------------------------------------------
        integer, optional                           :: VelocityID
        real, dimension(:, :, :), pointer, optional :: FatherMatrix, SonMatrix
        real, dimension(:, :),    pointer, optional :: FatherMatrix2D, SonMatrix2D
        real                                        :: LocalTimeDecay
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
                                           DecayTime            = LocalTimeDecay,                         &
                                           DT                   = Me%Hydro%VelDT,                         &
                                           SonVolInFather       = GetPointer(Me%Father%TotSonVolIn),                  &
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
                                           DecayTime            = LocalTimeDecay,                         &
                                           DT                   = Me%Hydro%VelDT,                         &
                                           SonVolInFather       = GetPointer(Me%Father%TotSonVolIn),                  &
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
                                      DecayTime        = LocalTimeDecay,                  &
                                      DT               = Me%Hydro%DT,                     &
                                      SonVolInFather   = GetPointer(Me%Father%TotSonVolIn),           &
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
                                  DecayTime        = LocalTimeDecay,                  &
                                  DT               = Me%Hydro%DT,                     &
                                  SonVolInFather2D = GetPointer(Me%Father%TotSonVolIn_2D),        &
                                  AuxMatrix2D      = GetPointer(Me%Father%AuxMatrix2D),           &
                                  VolumeSon2D      = Me%External_Var%VolumeZ_2D,       &
                                  VolumeFather2D   = Me%Father%External_Var%VolumeZ_2D, &
                                  IgnoreOBCells    = GetPointer(Me%IgnoreOBCells))               
                
        endif
  
    end subroutine Nudging_average
    
    !------------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Calls Feedback routines present in mondule functions, based on the type of input matrixes - IWD Method
    !>@param[in] FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, LocalTimeDecay
    subroutine Nudging_IWD (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, LocalTimeDecay)
        !Arguments-------------------------------------------------------------
        integer, optional                           :: VelocityID
        real, dimension(:, :, :), pointer, optional :: FatherMatrix, SonMatrix
        real, dimension(:, :),    pointer, optional :: FatherMatrix2D, SonMatrix2D
        real                                        :: LocalTimeDecay
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

    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Computes momentum discharges provided by a nested domain
    !>@param[in] FatherID, UpscalingMomentum, ComputeFaces3D, Velocity, KFloor, DischargesVelUV, di, dj
    subroutine ModifyUpscalingDischarge(FatherID, UpscalingMomentum, ComputeFaces3D, Velocity, KFloor, &
                                        DischargesVelUV, di, dj, I, J, STAT)
        !Arguments-------------------------------------------------------------
        integer                                        :: FatherID, di, dj, STAT, nCells, I, J
        real, dimension(:, :, :), pointer,             :: Velocity, FatherAreaU, FatherAreaV
        real, dimension(:, :, :), pointer, intent(OUT) :: UpscalingMomentum, DischargesVelUV
        integer, dimension(:, :, :), pointer,          :: ComputeFaces3D
        integer, dimension(:, :),    pointer,          :: KFloor
        !Local-----------------------------------------------------------------
        integer                                        :: ready_
        
        STAT = UNKNOWN_
        
        call Ready(FatherID, ready_)
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call GetDischargeFlowDistribuiton(Me%ObjDischarges, DischargeNumber, nCells, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - ModifyUpscalingDischarge - ERR01'
            
            if (nCells == 1)then
                
                !Buscar a velocidade proveniente da média ponderada do filho : ChildVelocity
                
                !DischargeVelocity => Me%Son%DischargeVelocity

                !Call GetGeometryAreas(FatherID, AreaU = FatherAreaU, AreaU = FatherAreaV, STAT_CALL)
                !if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - ModifyUpscalingDischarge - ERR10'
                
                !Child velocity => Me%Son%DischargeVelocityU
                Call ComputeUpscalingVelocity(ComputeFaces3D, 
                                              KFloor, 
                                              Velocity, 
                                              Me%Size, 
                                              AreaV, 
                                              DischargeVelocity, 
                                              di, 
                                              dj, 
                                              I, 
                                              J)
                    

            else
                write(*,*)'Model is not yet ready to accept any other discharge type than a point discharge'
                stop 'ModuleTwoWay - ModifyUpscalingDischarge - ERR20'
            endif
            
            iNorth = i+di
            jEast =  j+dj
                    
            
            
            
            STAT = SUCCESS_
        else
            STAT = ready_
            
        endif
        
    
    
    
    
    end subroutine ModifyUpscalingDischarge
    
    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Ungets all external vars
    !>@param[in] SonID, FatherID, CallerID, STAT    
    subroutine UngetTwoWayExternal_Vars(SonID, FatherID, CallerID, STAT)
    
        !Arguments--------------------------------------------------------------
        integer, intent(IN)                         :: SonID, FatherID, CallerID
        integer, optional, intent(OUT)              :: STAT
        !Locals-----------------------------------------------------------------
        integer                                     :: ready_, status, STAT_
        !Begin------------------------------------------------------------------
        
        STAT_ = UNKNOWN_      
        
        call Ready(SonID, ready_)        
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            if (callerID == 1) then
                !For future developments (when other modules call for twoway)
            endif            
            !Unget son
            call UngetHorizontalGrid(SonID, Me%External_Var%IV,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR01'
            call UngetHorizontalGrid(SonID, Me%External_Var%JV,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR10'
            call UngetHorizontalGrid(SonID, Me%External_Var%IU,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR20'
            call UngetHorizontalGrid(SonID, Me%External_Var%JU,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR30'
            call UngetHorizontalGrid(SonID, Me%External_Var%IZ,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR40'
            call UngetHorizontalGrid(SonID, Me%External_Var%JZ,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR50'
            call UnGetGeometry(SonID, Me%External_Var%VolumeU,     STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR60'
            call UnGetGeometry(SonID, Me%External_Var%VolumeV,     STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR70'
            call UnGetGeometry(SonID, Me%External_Var%VolumeZ,     STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR80'
            call UnGetGeometry(SonID, Me%External_Var%VolumeZ_2D,  STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR90'
            call UnGetMap(SonID, Me%External_Var%Open3D,           STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR100'
            call UnGetMap(SonID, Me%External_Var%ComputeFaces3D_U, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR110'
            call UnGetMap(SonID, Me%External_Var%ComputeFaces3D_V, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR120'
            
            if (Me%Hydro%InterpolationMethod == 2) then
                
                !call UnGetConnections (HorizontalGridID = SonID,   &
                !                     Connections_U     = Me%External_Var%IWD_Connections_U, &
                !                     IWD_Distances_U   = Me%External_Var%IWD_Distances_U,   &
                !                     Connections_V     = Me%External_Var%IWD_Connections_V, &
                !                     IWD_Distances_V   = Me%External_Var%IWD_Distances_V,   &
                !                     Connections_Z     = Me%External_Var%IWD_Connections_Z, &
                !                     IWD_Distances_Z   = Me%External_Var%IWD_Distances_Z,   &
                !                     STAT              = status)
                !if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-Hydrodynamic-ERR014'
                
                call UnGetHorizontalGrid(SonID, Me%External_Var%IWD_Connections_U, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR130'
                call UnGetHorizontalGrid(SonID, Me%External_Var%IWD_Connections_V, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR140'
                call UnGetHorizontalGrid(SonID, Me%External_Var%IWD_Connections_Z, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR150'
                call UnGetHorizontalGrid(SonID, Me%External_Var%IWD_Distances_U, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR160'
                call UnGetHorizontalGrid(SonID, Me%External_Var%IWD_Distances_V, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR170'
                call UnGetHorizontalGrid(SonID, Me%External_Var%IWD_Distances_Z, status)
                if (status /= SUCCESS_) stop 'UngetTwoWayExternal_Vars-TwoWay-ERR180'
                
            endif             
            
            !Unget father
            call UnGetMap(FatherID, Me%Father%External_Var%Open3D,           STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-Hydrodynamic-ERR190'
            call UnGetGeometry(FatherID, Me%Father%External_Var%VolumeZ,     STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-Hydrodynamic-ERR200'
            call UnGetGeometry(FatherID, Me%Father%External_Var%VolumeU,     STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-Hydrodynamic-ERR210'
            call UnGetGeometry(FatherID, Me%Father%External_Var%VolumeV,     STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-Hydrodynamic-ERR220'
            call UnGetGeometry(FatherID, Me%Father%External_Var%VolumeZ_2D,  STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-Hydrodynamic-ERR230'
            call UnGetMap(FatherID, Me%Father%External_Var%ComputeFaces3D_U, STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-Hydrodynamic-ERR240'
            call UnGetMap(FatherID, Me%Father%External_Var%ComputeFaces3D_V, STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-Hydrodynamic-ERR250'
            
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
    
        if (allocated(Me%Father%AuxMatrix)) then
            deallocate(Me%Father%AuxMatrix)
        endif
        if (allocated(Me%Father%AuxMatrix2D)) then
            deallocate(Me%Father%AuxMatrix2D)
        endif
        if (allocated(Me%Father%TotSonVolIn)) then
            deallocate(Me%Father%TotSonVolIn)
        endif
        if (allocated(Me%Father%TotSonVolIn_2D)) then
            deallocate(Me%Father%TotSonVolIn_2D)
        endif
        if (allocated(Me%IgnoreOBCells)) then
            deallocate(Me%IgnoreOBCells)
        endif
        if (allocated(Me%Father%IWDNom)) then
            deallocate(Me%Father%IWDNom)
        endif
        if (allocated(Me%Father%IWDDenom)) then
            deallocate(Me%Father%IWDDenom)
        endif

    end subroutine DeallocateVariables
    
    !-------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
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









