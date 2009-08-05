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


Module ModuleShell

    use ModuleGlobalData

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructShell
    private ::      AllocateInstance

    !Selector
    public  :: GetShellPointer
    public  :: GetShellInteger
    public  :: UnGetShell
                     
    
    !Modifier
    public  :: ModifyShell

    !Destructor
    public  :: KillShell                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjShell 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetShell3D_I
    private :: UnGetShell3D_R8
    interface  UnGetShell
        module procedure UnGetShell3D_I
        module procedure UnGetShell3D_R8
    end interface  UnGetShell

    !Types---------------------------------------------------------------------
    
    private :: T_Shell
    type       T_Shell
        integer                                     :: InstanceID
        type (T_Size3D)                             :: Size, WorkSize
        real(8), dimension(:, :, :),  pointer       :: Matrix
        type(T_Shell), pointer                      :: Next
    end type  T_Shell

    !Global Module Variables
    type (T_Shell), pointer                         :: FirstObjShell
    type (T_Shell), pointer                         :: Me

    integer                                         :: mSHELL_ = 0 !just to compile

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructShell(ObjShellID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjShellID 
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mSHELL_)) then
            nullify (FirstObjShell)
            call RegisterModule (mSHELL_) 
        endif

        call Ready(ObjShellID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            !Returns ID
            ObjShellID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleShell - ConstructShell - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructShell
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Shell), pointer                         :: NewObjShell
        type (T_Shell), pointer                         :: PreviousObjShell


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

        Me%InstanceID = RegisterNewInstance (mSHELL_)


    end subroutine AllocateInstance


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



    subroutine KillShell(ObjShellID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjShellID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjShellID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mSHELL_,  Me%InstanceID)

            if (nUsers == 0) then

                !Deallocates Instance
                call DeallocateInstance ()

                ObjShellID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillShell
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Shell), pointer          :: AuxObjShell
        type (T_Shell), pointer          :: PreviousObjShell

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

    subroutine Ready (ObjShell_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjShell_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjShell_ID > 0) then
            call LocateObjShell (ObjShell_ID)
            ready_ = VerifyReadLock (mSHELL_, Me%InstanceID)
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

end module ModuleShell









