!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Cuda Module
! PROJECT       : Mohid Base 1
! MODULE        : ModuleCude
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : September 2011
! REVISION      : Jonathan van der Wielen - v0.1
! DESCRIPTION   : Wrapper module to have access to C/C++ methods that use CUDA (project ModuleCuda in Libs)
! IMPLEMENTED   : Currently only ThomasZ, Thomas2D and Thomas3D are implemented in CUDA
!
!------------------------------------------------------------------------------

Module ModuleCuda

    use ModuleGlobalData

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructCuda
    private ::      AllocateInstance

    !Destructor
    public  :: KillCuda                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjCuda 
    
    public  :: Alloc3DPageLocked
    public  :: FreePageLocked
    
    !Thomas
    public  :: InitializeThomas
    !public  :: SolveThomas2D
    !public  :: SolveThomas3D
    public  :: SolveThomas  
    public  :: SaveThomas
    
    interface ConstructCudaBinding_C
        subroutine ConstructCudaBinding_C bind(c, name="ConstructCudaBinding_C")
            use, intrinsic  :: Iso_C_Binding
            implicit none
        end subroutine ConstructCudaBinding_C
    end interface ConstructCudaBinding_C
    
    interface KillCudaBinding_C
        subroutine KillCudaBinding_C bind(c, name="KillCudaBinding_C")
            use, intrinsic  :: Iso_C_Binding
            implicit none
        end subroutine KillCudaBinding_C
    end interface KillCudaBinding_C
    
    interface Alloc3DPageLocked_C
        subroutine Alloc3DPageLocked_C(Ptr, XDim, YDim, ZDim) bind(c, name="Alloc3DPageLocked_C")
            use, intrinsic          :: Iso_C_Binding
            implicit none
            
            integer(C_INT)          :: XDim, YDim, ZDim
            type(C_PTR)             :: Ptr
        end subroutine Alloc3DPageLocked_C
    end interface Alloc3DPageLocked_C
    
    interface FreePageLocked_C
        subroutine FreePageLocked_C(Ptr) bind(c, name="FreePageLocked_C")
            use, intrinsic          :: Iso_C_Binding
            implicit none
            
            type(C_PTR)             :: Ptr
        end subroutine FreePageLocked_C
    end interface FreePageLocked_C
    
    interface InitializeThomas_C
        subroutine InitializeThomas_C(ObjCudaID, Size) bind(c, name="InitializeThomas_C")
            use, intrinsic                  :: Iso_C_Binding
            ! import is needed to know the T_Size3D type within the interface
            import                          :: T_Size3D
            implicit none
            
            integer(C_INT)                  :: ObjCudaID
            type(T_Size3D)                  :: Size
        end subroutine InitializeThomas_C
    end interface InitializeThomas_C
    
    interface KillThomas_C
        subroutine KillThomas_C(ObjCudaID) bind(c, name="KillThomas_C")
            use, intrinsic  :: Iso_C_Binding
            implicit none

            integer(C_INT)      :: ObjCudaID
        end subroutine KillThomas_C
    end interface KillThomas_C
    
    interface SolveThomas_C
        subroutine SolveThomas_C(ObjCudaID, ILB, IUB, JLB, JUB, KLB, KUB, D, E, F, TI, Res, Dim) bind(c, name="SolveThomas_C")
            use, intrinsic  :: Iso_C_Binding
            implicit none
            
            integer(C_INT)                                  :: ObjCudaID, ILB, IUB, JLB, JUB, KLB, KUB, Dim
            real(C_DOUBLE), dimension(0:10000, 0:10000, *)  :: D, E, F, TI, Res
        end subroutine SolveThomas_C
    end interface SolveThomas_C
    
    interface SaveThomas_C
        subroutine SaveThomas_C(ObjCudaID, Res, Dim) bind(c, name="SaveThomas_C")
            use, intrinsic  :: Iso_C_Binding
            implicit none
            
            integer(C_INT)                                  :: ObjCudaID, Dim
            real(C_DOUBLE), dimension(0:10000, 0:10000, *)  :: Res
        end subroutine SaveThomas_C
    end interface
    
    !Types---------------------------------------------------------------------
    
    private :: T_Cuda
    type       T_Cuda
        integer                                     :: InstanceID
        type (T_Size3D)                             :: Size, WorkSize
        real(8), dimension(:, :, :),  pointer       :: Matrix
        type(T_Cuda), pointer                       :: Next
    end type  T_Cuda

    !Global Module Variables
    type (T_Cuda), pointer                          :: FirstObjCuda
    type (T_Cuda), pointer                          :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructCuda(ObjCudaID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjCudaID 
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mCUDA_)) then
            nullify (FirstObjCuda)
            call RegisterModule (mCUDA_) 
        endif

        call Ready(ObjCudaID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            !Returns ID
            ObjCudaID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleCuda - ConstructCuda - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructCuda
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Cuda), pointer                         :: NewObjCuda
        type (T_Cuda), pointer                         :: PreviousObjCuda


        !Allocates new instance
        allocate (NewObjCuda)
        nullify  (NewObjCuda%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjCuda)) then
            FirstObjCuda         => NewObjCuda
            Me                    => NewObjCuda
            
            ! If this is the first cuda module construct the binding to C/C++/Cuda. This initializes the CUDA device.
            call ConstructCudaBinding_C()
        else
            PreviousObjCuda      => FirstObjCuda
            Me                    => FirstObjCuda%Next
            do while (associated(Me))
                PreviousObjCuda  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjCuda
            PreviousObjCuda%Next => NewObjCuda
        endif

        Me%InstanceID = RegisterNewInstance (mCUDA_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillCuda(ObjCudaID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjCudaID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjCudaID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mCUDA_,  Me%InstanceID)

            if (nUsers == 0) then

                ! Kill the C++ Thomas instance that is associated with this ID.
                ! C++ code will check if an instance with this ID actually exists
                call KillThomas_C(ObjCudaID)
                
                !Deallocates Instance. Also kills CUDA if this is the last CUDA module
                call DeallocateInstance ()
                
                ObjCudaID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillCuda

    !------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Cuda), pointer          :: AuxObjCuda
        type (T_Cuda), pointer          :: PreviousObjCuda

        !Updates pointers
        if (Me%InstanceID == FirstObjCuda%InstanceID) then
            FirstObjCuda => FirstObjCuda%Next
            ! If this is the last cuda object, kill the cuda binding and reset the device
            if(.not. associated(FirstObjCuda)) then
                call KillCudaBinding_C()
            endif
        else
            PreviousObjCuda => FirstObjCuda
            AuxObjCuda      => FirstObjCuda%Next
            do while (AuxObjCuda%InstanceID /= Me%InstanceID)
                PreviousObjCuda => AuxObjCuda
                AuxObjCuda      => AuxObjCuda%Next
            enddo

            !Now update linked list
            PreviousObjCuda%Next => AuxObjCuda%Next

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

    subroutine Ready (ObjCuda_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjCuda_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjCuda_ID > 0) then
            call LocateObjCuda (ObjCuda_ID)
            ready_ = VerifyReadLock (mCUDA_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjCuda (ObjCudaID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjCudaID

        !Local-----------------------------------------------------------------

        Me => FirstObjCuda
        do while (associated (Me))
            if (Me%InstanceID == ObjCudaID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleCuda - LocateObjCuda - ERR01'

    end subroutine LocateObjCuda

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !ALLOCATE ALLOCATE ALLOCATE ALLOCATE ALLOCATE ALLOCATE ALLOCATE ALLOCATE AL

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------
    !JPW: TODO: use lowerbound / upperbound instead of _Dim
    subroutine Alloc3DPageLocked (ObjCudaID, Ptr, Arr, XDim, YDim, ZDim)
        
        !Arguments-------------------------------------------------------------
        
        integer(C_INT), intent(IN)                  :: ObjCudaID, XDim, YDim, ZDim
        type(C_PTR)                                 :: Ptr
        real(C_DOUBLE), dimension(:,:,:), pointer   :: Arr
        
        !External----------------------------------------------------------------
        integer                             :: ready_              
        
        !Local-----------------------------------------------------------------
        real(C_DOUBLE), dimension(:,:,:), pointer   :: TmpArr
        !----------------------------------------------------------------------
       
        call Ready(ObjCudaID, ready_)
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mCUDA_, ObjCudaID)
            ! Allocate memory to the pointer. XDim will be altered to have the correct pitch!
            call Alloc3DPageLocked_C(Ptr, XDim, YDim, ZDim)
            ! Associate the array with the pointer (default lower bound = 1)
            call c_f_pointer(Ptr, TmpArr, [XDim, YDim, ZDim])
            ! Remap from 1 to 0
            ! TODO: theoretically speaking this should be ILB, JLB, KLB instead of 0
#ifdef _USE_PAGELOCKED
            !Fortran 2003 onwards only. Requires ifort version >= 12.0
            Arr(0:, 0:, 0:) => TmpArr
#else
            !griflet - use this just to make it compilable with ifort version <= 11.1
            Arr => TmpArr
#endif _USE_PAGELOCKED
        end if

    end subroutine
    
    subroutine FreePageLocked (ObjCudaID, Ptr, Arr)
        
        !Arguments-------------------------------------------------------------
        
        integer                                     :: ObjCudaID
        type(C_PTR)                                 :: Ptr  
        real(C_DOUBLE), dimension(:,:,:), pointer   :: Arr
        
        !External----------------------------------------------------------------
        integer                             :: ready_              
        
        !----------------------------------------------------------------------
       
        call Ready(ObjCudaID, ready_)
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mCUDA_, ObjCudaID)
            
            call FreePageLocked_C(Ptr)
            
            nullify(Arr)
        end if

    end subroutine
    
    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !THOMAS THOMAS THOMAS THOMAS THOMAS THOMAS THOMAS THOMAS THOMAS THOMAS THOM

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------
    
    subroutine InitializeThomas (ObjCudaID, Size, STAT)
        
        !Arguments-------------------------------------------------------------
        
        integer                             :: ObjCudaID
        type(T_Size3D)                      :: Size
        integer, optional, intent(OUT)      :: STAT
        
        !External----------------------------------------------------------------
        integer                             :: ready_              
        
        !Local-------------------------------------------------------------------
        integer                             :: STAT_
        
        !----------------------------------------------------------------------
       
        call Ready(ObjCudaID, ready_)
        
        STAT_ = UNKNOWN_
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mCUDA_, ObjCudaID)
            
            ! Initialize Thomas in CUDA/C, use geometry size including boundary.
            ! NOTE: For simplicity we assume xLB is always 1 when sending the arrays.
            ! Dim: 0 = X, 1 = Y, 2 = Z
            call InitializeThomas_C(ObjCudaID, Size)
            
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) &
            STAT = STAT_

    end subroutine
    
    !--------------------------------------------------------------------------
    
    subroutine SolveThomas(ObjCudaID,               &
                            ILB, IUB,               &
                            JLB, JUB,               &
                            KLB, KUB,               &
                            D, E, F, TI, Res, Dim)
        !Arguments-------------------------------------------------------------
        integer(C_INT)                          :: ObjCudaID, Dim  
        integer(C_INT)                          :: ILB, IUB
        integer(C_INT)                          :: JLB, JUB
        integer(C_INT)                          :: KLB, KUB
        real(C_DOUBLE), dimension(:, :, :)      :: D, E, F, TI, Res
        
        !External--------------------------------------------------------------
        integer                         :: ready_              
        
        ! Select the correct module
        call Ready(ObjCudaID, ready_)
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mCUDA_, ObjCudaID)

            ! Solve the Thomas in CUDA
            call SolveThomas_C(ObjCudaID, ILB, IUB, JLB, JUB, KLB, KUB, D, E, F, TI, Res, Dim)
        
        end if

    end subroutine
    
    subroutine SaveThomas(ObjCudaID, Res, Dim)
        
        !Arguments-------------------------------------------------------------
        integer(C_INT)                          :: ObjCudaID, Dim 
        real(C_DOUBLE), dimension(:, :, :)      :: Res
        
        !External--------------------------------------------------------------
        integer                         :: ready_              
        
        ! Select the correct module
        call Ready(ObjCudaID, ready_)
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mCUDA_, ObjCudaID)

            ! Save the results to a file
            call SaveThomas_C(ObjCudaID, Res, Dim)
        
        end if
        
    end subroutine
    
end module ModuleCuda