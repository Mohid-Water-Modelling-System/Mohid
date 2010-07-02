!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : ClusterTimeSeriesByProp
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to serve as ClusterTimeSeriesByProp to create new modules
!
!------------------------------------------------------------------------------


Module ModuleClusterTimeSeriesByProp

    use ModuleGlobalData
    use ModuleTime                 
    use ModuleEnterData            
    use ModuleTimeSerie


    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructClusterTimeSeriesByProp
    private ::      AllocateInstance

    !Selector
    public  :: GetClusterTimeSeriesByPropInteger
                
    
    !Modifier
    public  :: ModifyClusterTimeSeriesByProp

    !Destructor
    public  :: KillClusterTimeSeriesByProp                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjClusterTimeSeriesByProp 
    
    !Interfaces----------------------------------------------------------------

    !Types---------------------------------------------------------------------
    
    private :: T_ClusterTimeSeriesByProp
    type       T_ClusterTimeSeriesByProp
        integer                                     :: InstanceID
        integer                                     :: ObjEnterData
        integer                                     :: ObjTimeSerieIn
        integer                                     :: ObjTimeSerieOut
        type(T_ClusterTimeSeriesByProp), pointer    :: Next

    end type  T_ClusterTimeSeriesByProp

    !Global Module Variables
    type (T_ClusterTimeSeriesByProp), pointer                         :: FirstObjClusterTimeSeriesByProp
    type (T_ClusterTimeSeriesByProp), pointer                         :: Me

    integer                                         :: mClusterTimeSeriesByProp_ = 0 !just to compile

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructClusterTimeSeriesByProp

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL

        !------------------------------------------------------------------------


            !Constructs EnterData
            call ConstructEnterData(Me%ObjEnterData, 'ClusterTimeSeriesByProp.dat', STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartTimeSerie - ModuleTimeSerie - ERR01'


            call ConstructTimeSeriesIn

            call ConstructTimeSeriesOut


        !----------------------------------------------------------------------

    end subroutine ConstructClusterTimeSeriesByProp
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_ClusterTimeSeriesByProp), pointer                         :: NewObjClusterTimeSeriesByProp
        type (T_ClusterTimeSeriesByProp), pointer                         :: PreviousObjClusterTimeSeriesByProp


        !Allocates new instance
        allocate (NewObjClusterTimeSeriesByProp)
        nullify  (NewObjClusterTimeSeriesByProp%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjClusterTimeSeriesByProp)) then
            FirstObjClusterTimeSeriesByProp         => NewObjClusterTimeSeriesByProp
            Me                    => NewObjClusterTimeSeriesByProp
        else
            PreviousObjClusterTimeSeriesByProp      => FirstObjClusterTimeSeriesByProp
            Me                    => FirstObjClusterTimeSeriesByProp%Next
            do while (associated(Me))
                PreviousObjClusterTimeSeriesByProp  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjClusterTimeSeriesByProp
            PreviousObjClusterTimeSeriesByProp%Next => NewObjClusterTimeSeriesByProp
        endif

        Me%InstanceID = RegisterNewInstance (mClusterTimeSeriesByProp_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    !--------------------------------------------------------------------------
    
    subroutine GetClusterTimeSeriesByPropInteger (ObjClusterTimeSeriesByPropID, Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjClusterTimeSeriesByPropID
        real                                            :: Int
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjClusterTimeSeriesByPropID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Int = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetClusterTimeSeriesByPropInteger

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyClusterTimeSeriesByProp(ObjClusterTimeSeriesByPropID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjClusterTimeSeriesByPropID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjClusterTimeSeriesByPropID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyClusterTimeSeriesByProp


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillClusterTimeSeriesByProp(ObjClusterTimeSeriesByPropID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjClusterTimeSeriesByPropID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjClusterTimeSeriesByPropID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mClusterTimeSeriesByProp_,  Me%InstanceID)

            if (nUsers == 0) then

                !Deallocates Instance
                call DeallocateInstance ()

                ObjClusterTimeSeriesByPropID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillClusterTimeSeriesByProp
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_ClusterTimeSeriesByProp), pointer          :: AuxObjClusterTimeSeriesByProp
        type (T_ClusterTimeSeriesByProp), pointer          :: PreviousObjClusterTimeSeriesByProp

        !Updates pointers
        if (Me%InstanceID == FirstObjClusterTimeSeriesByProp%InstanceID) then
            FirstObjClusterTimeSeriesByProp => FirstObjClusterTimeSeriesByProp%Next
        else
            PreviousObjClusterTimeSeriesByProp => FirstObjClusterTimeSeriesByProp
            AuxObjClusterTimeSeriesByProp      => FirstObjClusterTimeSeriesByProp%Next
            do while (AuxObjClusterTimeSeriesByProp%InstanceID /= Me%InstanceID)
                PreviousObjClusterTimeSeriesByProp => AuxObjClusterTimeSeriesByProp
                AuxObjClusterTimeSeriesByProp      => AuxObjClusterTimeSeriesByProp%Next
            enddo

            !Now update linked list
            PreviousObjClusterTimeSeriesByProp%Next => AuxObjClusterTimeSeriesByProp%Next

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

    subroutine Ready (ObjClusterTimeSeriesByProp_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjClusterTimeSeriesByProp_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjClusterTimeSeriesByProp_ID > 0) then
            call LocateObjClusterTimeSeriesByProp (ObjClusterTimeSeriesByProp_ID)
            ready_ = VerifyReadLock (mClusterTimeSeriesByProp_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjClusterTimeSeriesByProp (ObjClusterTimeSeriesByPropID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjClusterTimeSeriesByPropID

        !Local-----------------------------------------------------------------

        Me => FirstObjClusterTimeSeriesByProp
        do while (associated (Me))
            if (Me%InstanceID == ObjClusterTimeSeriesByPropID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleClusterTimeSeriesByProp - LocateObjClusterTimeSeriesByProp - ERR01'

    end subroutine LocateObjClusterTimeSeriesByProp

    !--------------------------------------------------------------------------

end module ModuleClusterTimeSeriesByProp









