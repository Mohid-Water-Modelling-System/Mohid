!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid TwoWay nesting
! PROJECT       : Mohid Base 1
! MODULE        : TwoWay
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jul 2018
! REVISION      : Joao Sobrinho - v1.0
!>  DESCRIPTION : Module that computes twoway operations for MohidWater and MohidLand
!
!> @author
!> Joao Sobrinho
!------------------------------------------------------------------------------


Module ModuleTwoWay

    use ModuleGlobalData
    use ModuleGeometry
    use ModuleHorizontalGrid
    use ModuleFuntions
    use ModuleStopWatch,        only : StartWatch, StopWatch
    
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructTwoWay
    private ::      AllocateInstance
    private ::      Add_Property
    
    
    public  :: ConstructTwoWayHydrodynamic
    public  :: Allocate2WayAuxiliars_Hydrodynamic    
    public  :: ConstructTwoWayWP

    

    !Selector
    public  :: GetTwoWayPointer
    public  :: GetTwoWayInteger
    public  :: UnGetTwoWay
                     
    
    !Modifier
    public  :: ModifyTwoWay
    private ::  ComputeAuxMatrixes
    private ::    ComputeSonVolInFather
    private ::    ComputeAuxMatrixes_RWAvg
    private ::  Nudging_average
    private ::  Nudging_IWD
    public  :: PrepTwoWay
    public  :: UngetTwoWayExternalVars
 

    !Destructor
    public  :: KillTwoWay                                                     
    private ::      DeAllocateInstance
    private ::      DeallocateVariables

    !Management
    private ::      Ready
    private ::          LocateObjTwoWay 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetTwoWay3D_I
    private :: UnGetTwoWay3D_R8
    interface  UnGetTwoWay
        module procedure UnGetTwoWay3D_I
        module procedure UnGetTwoWay3D_R8
    end interface  UnGetTwoWay

    !Types---------------------------------------------------------------------
    
    private :: T_Hydro
    type       T_Hydro
        real                                        :: TimeDecay           = 3600.
        integer                                     :: InterpolationMethod = 1
        real                                        :: VelDT               = null_real
        real                                        :: DT                  = null_real
    end type T_Hydro
    
    private :: T_WP
    type       T_WP
        
        type(T_Property)                            :: Property
        integer                                     :: NumberOfProperties
        
    end type T_WP
    
    private :: T_Property
    type       T_Property
        integer                                     :: ID                  = null_int
        character(len=PathLength)                   :: Name                = null_str
        logical                                     :: TwoWay              = .false.
        integer                                     :: InterpolationMethod = 1
        real                                        :: TimeDecay           = 3600.
        type (T_Property), pointer                  :: Next                => null()
    end type T_Property
    
    private :: T_FatherDomain
    type       T_FatherDomain
        type (T_Size3D)                             :: Size, WorkSize
        type (T_External)                           :: ExternalVar
        integer                                     :: InstanceID
        real, dimension (:, :, :), pointer          :: TotSonVolIn          => null()
        real, dimension (:, :   ), pointer          :: TotSonVolIn2D        => null()
        real, dimension (:, :, :), pointer          :: AuxMatrix            => null()
        real, dimension (:, :   ), pointer          :: AuxMatrix2D          => null()
    end type T_FatherDomain
    
    private :: T_External
    type       T_External
        integer, dimension(:, :   ), pointer        :: IV
        integer, dimension(:, :   ), pointer        :: JV
        integer, dimension(:, :   ), pointer        :: IU
        integer, dimension(:, :   ), pointer        :: JU
        integer, dimension(:, :   ), pointer        :: IZ
        integer, dimension(:, :   ), pointer        :: JZ
        real,    dimension(:, :, :), pointer        :: VolumeU
        real,    dimension(:, :, :), pointer        :: VolumeV
        real,    dimension(:, :, :), pointer        :: VolumeZ
        real,    dimension(:, :, :), pointer        :: VolumeZ_2D
        real,    dimension(:, :, :), pointer        :: Open3D
    end type T_External
    
    
    private :: T_TwoWay
    type       T_TwoWay
        integer                                     :: InstanceID
        character(PathLength)                       :: ModelName
        real(8), dimension(:, :, :),  pointer       :: Matrix

        type(T_ExternalVar)                         :: ExternalVar
        type(T_Property), pointer                   :: FirstProperty
        type(T_Property), pointer                   :: LastProperty
        type (T_Hydro)                              :: Hydro
        type (T_WP)                                 :: WP
        type (T_Size3D)                             :: Size, WorkSize
        type (T_FatherDomain)                       :: Father
        type(T_TwoWay), pointer                     :: Next
        
        !Instance of ModuleHorizontalGrid
        integer                                     :: ObjHorizontalGrid = 0
        !Instance of ModuleGeometry
        integer                                     :: ObjGeometry       = 0 
        
    end type  T_TwoWay
    


    !Global Module Variables
    type (T_TwoWay), pointer                         :: FirstObjTwoWay
    type (T_TwoWay), pointer                         :: Me

    integer                                         :: mTwoWay_ = 0 !just to compile

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
    subroutine ConstructTwoWay(ModelName, ObjTwoWayID, STAT)

        !Arguments---------------------------------------------------------------
        character(Len=*)                                :: ModelName
        integer                                         :: ObjTwoWayID 
        integer                                         :: STAT     
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

        call Ready(ObjTwoWayID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            Me%ModelName = ModelName
            
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            
            !Returns ID
            ObjTwoWayID          = Me%InstanceID
            
            call GetGeometrySize(Me%ObjGeometry,         &
                                 Size     = Me%Size,     &
                                 WorkSize = Me%WorkSize, &
                                 STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                  &
                stop 'ModuleTwoWay - ConstructTwoWay - ERR01'            

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
    !>@param[in] TwoWayID, TimeDecay, IntMethod, VelDT, DT
    subroutine ConstructTwoWayHydrodynamic (TwoWayID, TimeDecay, IntMethod, VelDT, DT, STAT)
    
        !Arguments------------------------------------------------------------
        integer                                     :: TwoWayID, & ! ID
                                                       IntMethod !Method for grid interpolation from child to father grid
                                                              ! 1 - Volume Averaged; 2 - Inverse Weighted Distance
        real                                        :: TimeDecay ! Decay factor in seconds, in the nudging equation
        real                                        :: VelDT, DT
        
        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL, TwoWayID, ready_, STAT_
        !---------------------------------------------------------------------

        call Ready(TwoWayID, ready_)    

        if (ready_ .EQ. OFF_ERR_) then
            
            Me%Hydro%InterpolationMethod = IntMethod
            Me%Hydro%TimeDecay           = TimeDecay
            Me%Hydro%VelDT               = VelDT
            Me%Hydro%DT                  = DT
            
            STAT_ = SUCCESS_

        else
            stop 'ModuleTwoWay - ConstructTwoWayHydrodynamic - ERR02' 

        end if

        if (present(STAT)) STAT = STAT_        
        
    end subroutine ConstructTwoWayHydrodynamic
    
    !-------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Allocates auxiliar matrixes  
    !>@param[in] FatherTwoWayID, TwoWayID   
    subroutine Allocate2WayAuxiliars_Hydrodynamic(FatherTwoWayID, TwoWayID, IntMethod)
    
        !Arguments-------------------------------------------------------------
        integer                            :: FatherTwoWayID, TwoWayID, IntMethod
        !Local-----------------------------------------------------------------

        integer                            :: ready_, VarType

        !----------------------------------------------------------------------
            
        call Ready (TwoWayID, ready_)
        
        if (ready_ .EQ. OFF_ERR_)then
            Me%Father%InstanceID = FatherID
            
        call GetGeometrySize(FatherTwoWayID,                   &
                                Size     = Me%Father%Size,     &
                                WorkSize = Me%Father%WorkSize, &
                                STAT     = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                  &
            stop 'ModuleTwoWay - Allocate2WayAuxiliars_Hydrodynamic - ERR01'  
            
            ILB = Me%WorkSize%ILB 
            IUB = Me%WorkSize%IUB 
            JLB = Me%WorkSize%JLB 
            JUB = Me%WorkSize%JUB 
            KLB = Me%WorkSize%KLB 
            KUB = Me%WorkSize%KUB            
            
            allocate(Me%Father%AuxMatrix(ILB:IUB, JLB:JUB, KLB:KUB))
            Me%Father%AuxMatrix(:,:,:) = 0.0
                
            allocate(Me%Father%AuxMatrix2D(ILB:IUB, JLB:JUB))
            Me%Father%AuxMatrix2D(:,:) = 0.0                
            
            if (Me%Hydro%InterpolationMethod == 1)then
            
                allocate(Me%Father%TotSonVolIn(ILB:IUB, JLB:JUB, KLB:KUB))
                Me%Father%TotSonVolIn(:,:,:) = 0.0
            
                allocate(Me%Father%TotSonVolIn2D(ILB:IUB, JLB:JUB))
                Me%Father%TotSonVolIn2D(:,:) = 0.0
            else
                VarType = 2 !One day will a user can define this
                call ConstructIWDTwoWay (FatherID, TwoWayID, VarType)         
            endif

        else
            stop 'ModuleTwoWay - Allocate2WayAuxiliars_Hydrodynamic - ERR01'               
        endif

    end subroutine Allocate2WayAuxiliars_Hydrodynamic
    
    !-------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Contructs WaterProperties variables of MOHIDWater with TwoWay nesting  
    !>@param[in] TwoWayID, WaitPeriod, TimeDecay, Method, Continuous, FatherSize3D
                                            
    subroutine ConstructTwoWayWP
    
        !Arguments------------------------------------------------------------
        
        !Locals---------------------------------------------------------------
        type (T_Property), pointer                      :: NewProperty
        
        !------------------------------------------------------------------------
        
        !Will only get the necessary information : ID, Name, TwoWay, WaitPeriod
        call GetNumberOfProperties    (WaterPropertiesID      = ObjTwoWayID,               &
                                       PropertiesNumber       = Me%WP%NumberOfProperties,  &
                                       STAT                   = STAT_CALL)
        if (STAT_CALL /= SUCCESS) then
            write(*,*) 'Could not get number of properties in WaterProperties module:'
            stop 'ModuleTwoWay - ConstructMohidWaterVariables - ERR02'
        endif
        
        !Colects properties information relevant for TwoWay and saves it in this module
        nullify (Me%FirstProperty)
        nullify (Me%LastProperty)
        do i = 1,  Me%WP%NumberOfProperties
            
            allocate (NewProperty)
            nullify  (NewProperty%PropertyID,   &
                      NewProperty%PropertyName, &
                      NewProperty%TwoWay,       &
                      NewProperty%WaitPeriod,   &
                      NewProperty%TimeDecay)
            
            call GetWaterPropertiesOptions(WaterPropertiesID  = ObjTwoWayID,                &
                                           PropertyID         = NewProperty%PropertyID,     &
                                           PropertyName       = NewProperty%PropertyName,   &
                                           TwoWay             = NewProperty%TwoWay,         &
                                           WaitPeriod         = NewProperty%WaitPeriod,     &
                                           TimeDecay          = NewProperty%TimeDecay,      &
                                           STAT               = STAT_CALL)
            if (STAT_CALL /= SUCCESS) then
                write(*,*) 'Could not get information from one property in WaterProperties module:'
                write(*,*) 'ID, Name, TwoWay, WaitPeriod, TimeDecay, TwoWay_Method, continuous'
                stop 'ModuleTwoWay - ConstructMohidWaterVariables - ERR03'
            endif
            
            call Add_Property(NewProperty)

        enddo
         
        
    end subroutine ConstructTwoWayWP
    
    !-------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Adds a new property into the properties list
    !>@param[in] NewProperty
    subroutine Add_Property(NewProperty)
    
        !Arguments-------------------------------------------------------------
        type(T_Property),           pointer     :: NewProperty
        !----------------------------------------------------------------------

        if (.not.associated(Me%FirstProperty)) then
            Me%FirstProperty    => NewProperty
            Me%LastProperty     => NewProperty
        else
            NewProperty%Prev     => Me%LastProperty
            Me%LastProperty%Next => NewProperty
            Me%LastProperty      => NewProperty
        end if     
    
    end subroutine Add_Property
    
    !-------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !--------------------------------------------------------------------------
    subroutine GetTwoWayPointer (ObjTwoWayID, Matrix, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjTwoWayID
        real(8), dimension(:, :, :),  pointer           :: Matrix
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTwoWayID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mTwoWay_, Me%InstanceID)

            Matrix => Me%Matrix

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetTwoWayPointer
    
    !--------------------------------------------------------------------------
    
    subroutine GetTwoWayInteger (ObjTwoWayID, Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjTwoWayID
        real                                            :: Int
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTwoWayID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Int = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetTwoWayInteger

    !--------------------------------------------------------------------------

    subroutine UnGetTwoWay3D_I(ObjTwoWayID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjTwoWayID
        integer, dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTwoWayID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mTwoWay_, Me%InstanceID, "UnGetTwoWay3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetTwoWay3D_I

    !--------------------------------------------------------------------------

    subroutine UnGetTwoWay3D_R8(ObjTwoWayID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjTwoWayID
        real(8), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjTwoWayID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mTwoWay_, Me%InstanceID,  "UnGetTwoWay3D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetTwoWay3D_R8

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Updates father grid domain with son's results
    !>@param[in] SonID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, CallerID, STAT
    subroutine ModifyTwoWay(SonID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, CallerID, VelocityID, STAT)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: SonID, CallerID
        integer, optional                           :: VelocityID
        integer, intent(OUT)                        :: STAT
        real, dimension(:, :, :), pointer, optional :: FatherMatrix, SonMatrix
        real, dimension(:, :),    pointer, optional :: FatherMatrix2D, SonMatrix2D
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, InterpolMethod
        real                                        :: LocalTimeDecay

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SonID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            if (MonitorPerformance) call StartWatch ("ModuleHydrodynamic", "Modify_Hydrodynamic")
            
            if (CallerID == mHydrodynamic_) then
                LocalTimeDecay = Me%Hydro%TimeDecay
                InterpolMethod = Me%Hydro%InterpolationMethod
            endif
            
            !if it is a 3D matrix
            if (present(FatherMatrix)) then
                if (present(VelocityID))then
                    if (VelocityID == VelocityU_) then
                        !if 3D matrixes were sent. Even2D domains allocate a 3D matrix (only one vertical layer)
                        call ComputeAuxMatrixes (SonID            = SonID,             &
                                                 FatherMatrix     = FatherMatrix,      &
                                                 SonMatrix        = SonMatrix,         &
                                                 InterpolMethod   = InterpolMethod,    &
                                                 Ilink            = Me%ExternalVar%IU, &
                                                 Jlink            = Me%ExternalVar%JU)
                    else
                        
                        call ComputeAuxMatrixes (SonID            = SonID,             &
                                                 FatherMatrix     = FatherMatrix,      &
                                                 SonMatrix        = SonMatrix,         &
                                                 InterpolMethod   = InterpolMethod,    &
                                                 Ilink            = Me%ExternalVar%IV, &
                                                 Jlink            = Me%ExternalVar%JV)
                    endif
                else
                    
                    call ComputeAuxMatrixes (SonID            = SonID,             &
                                             FatherMatrix     = FatherMatrix,      &
                                             SonMatrix        = SonMatrix,         &
                                             InterpolMethod   = InterpolMethod,    &
                                             Ilink            = Me%ExternalVar%IZ, &
                                             Jlink            = Me%ExternalVar%JZ)
                endif
                
            else
                !if a 2D matrix was sent (specific for waterLevel - at least for MohidWater).
                call ComputeAuxMatrixes (SonID            = SonID,             &
                                         FatherMatrix2D   = FatherMatrix2D,    &
                                         SonMatrix2D      = SonMatrix2D,       &
                                         InterpolMethod   = InterpolMethod,    &
                                         Ilink            = Me%ExternalVar%IZ, &
                                         Jlink            = Me%ExternalVar%JZ)   
            endif
            !Mandar esta parte para outra routina
            if (InterpolMethod == 1) then
                call Nudging_average (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, LocalTimeDecay)
            elseif (InterpolMethod == 2) then
                call Nudging_IWD (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, LocalTimeDecay)
            endif

           
            if (MonitorPerformance) call StopWatch ("ModuleHydrodynamic", "Modify_Hydrodynamic")  

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
        integer                                     :: STAT, STAT_
        !Locals-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        !Begin------------------------------------------------------------------
        
        call Ready(SonID, ready_)
        
        if (ready_ .EQ. IDLE_ERR_) then
            if (callerID == 1) then
                !For future developments (when other modules call for twoway)
            endif
            
            call GetHorizontalGrid (HorizontalGridID = SonID,                         &
                                    IV            = Me%ExternalVar%IV,                &
                                    JV            = Me%ExternalVar%JV,                &
                                    IU            = Me%ExternalVar%IU,                &
                                    JU            = Me%ExternalVar%JU,                &
                                    IZ            = Me%ExternalVar%IZ,                &
                                    JZ            = Me%ExternalVar%JZ,                &
                                    STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - PrepTwoWay - ERR01'
            
            call GetGeometryVolumes(GeometryID    = SonID,                            &
                                    VolumeU       = Me%ExternalVar%VolumeU,           &
                                    VolumeV       = Me%ExternalVar%VolumeV,           &
                                    VolumeZ       = Me%ExternalVar%VolumeZ,           &
                                    VolumeZ_2D    = Me%ExternalVar%VolumeZ_2D,        &
                                    STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - PrepTwoWay - ERR02'
            
            call GetOpenPoints3D   (Map_ID        = SonID,                            &
                                    OpenPoint3D   = Me%ExternalVar%Open3D,            &
                                    STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - PrepTwoWay - ERR03'
                           

            
            call GetOpenPoints3D   (Map_ID        = FatherID,                         &
                                    OpenPoint3D   = Me%Father%ExternalVar%Open3D,     &
                                    STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - PrepTwoWay - ERR04'
        
            call GetGeometryVolumes(GeometryID    = SonID, &
                                    VolumeU       = Me%Father%ExternalVar%VolumeU,    &
                                    VolumeV       = Me%Father%ExternalVar%VolumeV,    &
                                    VolumeZ       = Me%Father%ExternalVar%VolumeZ,    &
                                    VolumeZ_2D    = Me%Father%ExternalVar%VolumeZ_2D, &
                                    STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleTwoWay - PrepTwoWay - ERR05'           
        else
            STAT_ = ready_
        endif
        
        if present(STAT) STAT = STAT_
        
    end subroutine PrepTwoWay
    
    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Computes auxiliar matrixes for the feedback
    !>@param[in] SonID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, InterpolMethod    
    subroutine ComputeAuxMatrixes(SonID, FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, InterpolMethod, &
                                  Ilink, Jlink)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: SonID, interpolMethod
        integer, intent(OUT)                        :: STAT
        real, dimension(:, :, :), pointer, optional :: FatherMatrix, SonMatrix
        real, dimension(:, :),    pointer, optional :: FatherMatrix2D, SonMatrix2D
        integer, dimension(:, :), pointer           :: Ilink, Jlink
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        real                                        :: LocalTimeDecay

        !---------------------------------------------------------------------- 
        
        if (present(FatherMatrix)) then
            call SetMatrixValue (Me%Father%AuxMatrix, Me%WorkSize, 0.0)           
            
            !Goes for 3D
            if     (interpolMethod == 1)then
                
                call SetMatrixValue (Me%Father%TotSonVolIn, Me%WorkSize, 0.0)
                
                ! Volume Weighted average
                call ComputeSonVolInFather   (FatherMatrix   = FatherMatrix,   &
                                              SonMatrix      = SonMatrix,      &
                                              Ilink          = Ilink,          &
                                              Jlink          = Jlink) 
                
            elseif (interpolMethod == 2) then
                
                ! Inverse weighted distance method (includes a search radiuous)
                call GetAuxMatrixes_IWD  (    SonID          = SonID,           &
                                              FatherMatrix   = FatherMatrix,    &
                                              SonMatrix      = SonMatrix) 
                
            endif      
        !Goes for 2D 
        else
            call SetMatrixValue (Me%Father%AuxMatrix2D, Me%WorkSize, 0.0)
            
            if     (interpolMethod == 1)then
                
                call SetMatrixValue (Me%Father%TotSonVolIn2D, Me%WorkSize, 0.0)
                
                ! Volume Weighted average
                call ComputeSonVolInFather   (FatherMatrix2D = FatherMatrix2D,  &
                                              SonMatrix2D    = SonMatrix2D,     &
                                              Ilink          = Ilink,           &
                                              Jlink          = Jlink)

            elseif (interpolMethod == 2) then
                
                ! Inverse weighted distance method (includes a search radiuous)
                call GetAuxMatrixes_IWD      (SonID          = SonID,            &
                                              FatherMatrix2D = FatherMatrix2D,   &
                                              SonMatrix2D    = SonMatrix2D) 
                
            endif                
        endif
        
        
    end subroutine ComputeAuxMatrixes
    
    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Computes Son cells volume inside each father cell
    !>@param[in] FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D      
    subroutine ComputeSonVolInFather (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, Ilink, Jlink)
    
        !Arguments--------------------------------------------------------------------------------
        real, dimension(:, :, :), pointer, optional :: FatherMatrix, SonMatrix
        real, dimension(:, :),    pointer, optional :: FatherMatrix2D, SonMatrix2D
        integer, dimension(:, :), pointer           :: Ilink, Jlink
        !Local variables--------------------------------------------------------------------------
        integer                                 :: i, j, k
        !Begin------------------------------------------------------------------------------------
        if (present(FatherMatrix)) then
            
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    Me%Father%TotSonVolIn(ILink(i, j)+1, JLink(i, j)+1, k) =                   &
                    Me%Father%TotSonVolIn(ILink(i, j)+1, JLink(i, j)+1, k) + Matrix(i, j, k) * &
                    Me%ExternalVar%Open3D(i, j, k)
            enddo        
            enddo
            enddo
            
        else
            
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                Me%Father%TotSonVolIn2D(ILink(i, j)+1, JLink(i, j)+1) = &
                Me%Father%TotSonVolIn2D(ILink(i, j)+1, JLink(i, j)+1) + Matrix(i, j) * &
                Me%ExternalVar%Open3D(i, j, WorkSize%KUB)
            enddo        
            enddo            
            
        endif
          
    end subroutine ComputeSonVolInFather
    
    !---------------------------------------------------------------------------
    subroutine Nudging_average (FatherMatrix, SonMatrix, FatherMatrix2D, SonMatrix2D, VelocityID, LocalTimeDecay)
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: SonID, CallerID
        integer, optional                           :: VelocityID
        integer, intent(OUT)                        :: STAT
        real, dimension(:, :, :), pointer, optional :: FatherMatrix, SonMatrix
        real, dimension(:, :),    pointer, optional :: FatherMatrix2D, SonMatrix2D
        real                                        :: LocalTimeDecay
        !Local-----------------------------------------------------------------
        if (present(FatherMatrix)) then
                
            if (present(VelocityID))then
                !If velocity U
                if (VelocityID == VelocityU_)then
                    call FeedBack_Avrg (FatherMatrix   = FatherMatrix,                 &
                                        SonMatrix      = SonMatrix,                    &
                                        Open3DFather   = Me%Father%ExternalVar%Open3D, &
                                        Open3DSon      = Me%ExternalVar%Open3D,        &
                                        SizeFather     = Me%Father%WorkSize,           &
                                        SonWorkSize    = Me%WorkSize,                  &
                                        Ilink          = Me%ExternalVar%IU,            &
                                        Jlink          = Me%ExternalVar%JU,            &
                                        DecayTime      = LocalTimeDecay,               &
                                        DT             = Me%Hydro%VelDT,               &
                                        SonVolInFather = Me%Father%TotSonVolIn,        &
                                        AuxMatrix      = Me%Father%AuxMatrix,          &
                                        VolumeSon      = Me%ExternalVar%VolumeU,       &
                                        VolumeFather   = Me%Father%ExternalVar%VolumeU)
                !If velocity V    
                else
                    call FeedBack_Avrg (FatherMatrix   = FatherMatrix,                 &
                                        SonMatrix      = SonMatrix,                    &
                                        Open3DFather   = Me%Father%ExternalVar%Open3D, &
                                        Open3DSon      = Me%ExternalVar%Open3D,        &
                                        SizeFather     = Me%Father%WorkSize,           &
                                        SonWorkSize    = Me%WorkSize,                  &
                                        Ilink          = Me%ExternalVar%IV,            &
                                        Jlink          = Me%ExternalVar%JV,            &
                                        DecayTime      = LocalTimeDecay,               &
                                        DT             = Me%Hydro%VelDT,               &
                                        SonVolInFather = Me%Father%TotSonVolIn,        &
                                        AuxMatrix      = Me%Father%AuxMatrix,          &
                                        VolumeSon      = Me%ExternalVar%VolumeV,       &
                                        VolumeFather   = Me%Father%ExternalVar%VolumeV)                        
                endif
            endif
            !compute nudging Z type cell
            call FeedBack_Avrg   (FatherMatrix     = FatherMatrix,                    &
                                  SonMatrix        = SonMatrix,                       &
                                  Open3DFather     = Me%Father%ExternalVar%Open3D,    &
                                  Open3DSon        = Me%ExternalVar%Open3D,           &
                                  SizeFather       = Me%Father%WorkSize,              &
                                  SonWorkSize      = Me%WorkSize,                     &
                                  Ilink            = Me%ExternalVar%IZ,               &
                                  Jlink            = Me%ExternalVar%JZ,               &
                                  DecayTime        = LocalTimeDecay,                  &
                                  DT               = Me%Hydro%DT,                     &
                                  SonVolInFather   = Me%Father%TotSonVolIn,           &
                                  AuxMatrix        = Me%Father%AuxMatrix,             &
                                  VolumeSon        = Me%ExternalVar%VolumeZ,          &
                                  VolumeFather     = Me%Father%ExternalVar%VolumeZ)
                
        else
                
            call FeedBack_AvrgWL (FatherMatrix2D   = FatherMatrix2D,                  &
                                  SonMatrix2D      = SonMatrix2D,                     &
                                  Open3DFather     = Me%Father%ExternalVar%Open3D,    &
                                  Open3DSon        = Me%ExternalVar%Open3D,           &
                                  SizeFather       = Me%Father%WorkSize,              &
                                  SonWorkSize      = Me%WorkSize,                     &
                                  Ilink            = Me%ExternalVar%IZ,               &
                                  Jlink            = Me%ExternalVar%JZ,               &
                                  DecayTime        = LocalTimeDecay,                  &
                                  DT               = Me%Hydro%DT,                     &
                                  SonVolInFather   = Me%Father%TotSonVolIn_2D,        &
                                  AuxMatrix2D      = Me%Father%AuxMatrix2D,           &
                                  VolumeSon2D      = Me%ExternalVar%VolumeZ_2D,       &
                                  VolumeFather2D   = Me%Father%ExternalVar%VolumeZ_2D)               
                
        endif

        !----------------------------------------------------------------------    
    end subroutine Nudging_average
    
    !---------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !>Ungets all external vars
    !>@param[in] SonID, FatherID, CallerID, STAT    
    subroutine UngetTwoWayExternalVars(SonID, FatherID, CallerID, STAT)
    
        !Arguments--------------------------------------------------------------
        integer, intent(IN)                         :: SonID, FatherID, CallerID
        integer, intent(OUT)                        :: STAT
        !Locals-----------------------------------------------------------------
        integer                                     :: ready_, status, STAT_
        !Begin------------------------------------------------------------------
        
        if (callerID == 1) then
            !For future developments (when other modules call for twoway)
        endif
        call Ready(SonID, ready_)
        if (ready_ .EQ. IDLE_ERR_) then 
            
            !Unget son
            call UngetHorizontalGrid(SonID, Me%ExternalVar%IV,   STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternalVars-TwoWay-ERR01'
            call UngetHorizontalGrid(SonID, Me%ExternalVar%JV,   STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternalVars-TwoWay-ERR02'
            call UngetHorizontalGrid(SonID, Me%ExternalVar%IU,   STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternalVars-TwoWay-ERR03'
            call UngetHorizontalGrid(SonID, Me%ExternalVar%JU,   STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternalVars-TwoWay-ERR04'
            call UngetHorizontalGrid(SonID, Me%ExternalVar%IZ,   STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternalVars-TwoWay-ERR05'
            call UngetHorizontalGrid(SonID, Me%ExternalVar%JZ,   STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternalVars-TwoWay-ERR06'
            call UnGetGeometry(SonID, Me%ExternalVar%VolumeU,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternalVars-TwoWay-ERR07'
            call UnGetGeometry(SonID, Me%ExternalVar%VolumeV,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternalVars-TwoWay-ERR08'
            call UnGetGeometry(SonID, Me%ExternalVar%VolumeZ,    STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternalVars-TwoWay-ERR09'
            call UnGetGeometry(SonID, Me%ExternalVar%VolumeZ_2D, STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternalVars-TwoWay-ERR10'
            call UnGetMap(SonID, Me%ExternalVar%Open3D,          STAT = status)
            if (status /= SUCCESS_) stop 'UngetTwoWayExternalVars-TwoWay-ERR11'
            
            !Unget father
            call UnGetMap(FatherID, Me%Father%ExternalVar%Open3D,          STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-Hydrodynamic-ERR12'
            call UnGetGeometry(FatherID, Me%Father%ExternalVar%VolumeZ,    STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-Hydrodynamic-ERR13'
            call UnGetGeometry(FatherID, Me%Father%ExternalVar%VolumeU,    STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-Hydrodynamic-ERR14'
            call UnGetGeometry(FatherID, Me%Father%ExternalVar%VolumeV,    STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-Hydrodynamic-ERR15'
            call UnGetGeometry(FatherID, Me%Father%ExternalVar%VolumeZ_2D, STAT = status)
            if (status /= SUCCESS_) stop 'UnGetExternal2WayAuxVariables-Hydrodynamic-ERR16'
        else
            STAT_ = ready_
        endif
        
        if (present(STAT)) STAT = STAT_
        
    end subroutine UnGetExternal2WayAuxVariables    
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
        if (allocated(Me%Father%TotSonVolIn2D)) then
            deallocate(Me%Father%TotSonVolIn2D)
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

        !----------------------------------------------------------------------

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









