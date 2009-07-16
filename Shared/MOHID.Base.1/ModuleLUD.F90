!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : LUD
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : solve linear equations system by Gauss elimination method
!                 
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

module ModuleLUD

    use ModuleGlobalData

    implicit none

    private

    !Subroutine----------------------------------------------------------------

    !Contructor
    public  :: StartLUD
    private ::      AllocateInstance

    !Modifier
    public  :: LUD
    private ::      Pivot
    private ::      solve
    private ::      decomp
    private ::      order


    !Destructor
    public  :: KillLUD
    private ::      DeallocateInstance
    

    !Management
    private ::      Ready
    private ::          LocateObjLUD

    !Type----------------------------------------------------------------------
    
    private :: T_External
    type       T_External
        double precision, pointer, dimension(:,:) :: Coefficients
        real,             pointer, dimension(:  ) :: IndTerm
    end type T_External


    private :: T_LUD
    type       T_LUD  
        integer                                 :: InstanceID
        type(T_Size1D  )                        :: Size
        type(T_Size1D  )                        :: WorkSize
        type(T_External)                        :: ExternalVar
        integer, pointer, dimension(:)          :: o             !Order vector
        real,    pointer, dimension(:)          :: x             !Unknowns
        real,    pointer, dimension(:)          :: s             !Scale vector
        type(T_LUD     ), pointer :: Next
    end type   T_LUD


    !Global Module Variables
    type (T_LUD), pointer                       :: FirstObjLUD
    type (T_LUD), pointer                       :: Me


    !--------------------------------------------------------------------------

    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartLUD(LUD_ID, SizeLB, SizeUB, WorkSizeLB, WorkSizeUB, STAT)
                        

        !Arguments-------------------------------------------------------------
        integer                         :: LUD_ID
        integer,           intent(IN )  :: SizeLB,     SizeUB
        integer,           intent(IN )  :: WorkSizeLB, WorkSizeUB 
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                         :: ready_         

        !Local-----------------------------------------------------------------
        integer                         :: STAT_              !Auxiliar local variable
 
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mLUD_)) then
            nullify (FirstObjLUD)
            call RegisterModule (mLUD_) 
        endif
        
        call Ready(LUD_ID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
                               
            Me%Size%ILB     = SizeLB
            Me%Size%IUB     = SizeUB
            Me%WorkSize%ILB = WorkSizeLB
            Me%WorkSize%IUB = WorkSizeUB

            allocate(Me%o(Me%Size%ILB:Me%Size%IUB))
            allocate(Me%s(Me%Size%ILB:Me%Size%IUB))
            allocate(Me%x(Me%Size%ILB:Me%Size%IUB))

            Me%o = null_int
            Me%s = null_real
            Me%x = null_real

            !Returns ID
            LUD_ID   = Me%InstanceID

            STAT_    = SUCCESS_

        else cd0
            
            stop 'ModuleLUD - StartLUD - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartLUD
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine AllocateInstance


        !Local-----------------------------------------------------------------
        type (T_LUD), pointer                   :: NewObjLUD
        type (T_LUD), pointer                   :: PreviousObjLUD


        !Allocates new instance
        allocate (NewObjLUD)
        nullify  (NewObjLUD%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjLUD)) then
            FirstObjLUD         => NewObjLUD
            Me                  => NewObjLUD
        else
            PreviousObjLUD      => FirstObjLUD
            Me                  => FirstObjLUD%Next
            do while (associated(Me))
                PreviousObjLUD  => Me
                Me              => Me%Next
            enddo
            Me                  => NewObjLUD
            PreviousObjLUD%Next => NewObjLUD
        endif

        Me%InstanceID = RegisterNewInstance (mLUD_)


    end subroutine AllocateInstance

 
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    !Rotina de resolucao de sistemas de equacoes lineares por eliminacao de Gauss
    !Baseado em 'Numerical Methods for Engineers', Steven C. Chapra, Raymond P. Canale, 
    !second edition, pag. 282.
    !
    !--------------------------------------------------------------------------
    !
    ! Just as for Gauss slimination, LU decomposition algorithms must employ partial pivoting to 
    ! avoid division by zero and to minimize round-off error. The pivoting implemented imediately
    ! after computing each column of [L].
    ! In contrast to Gauss elimination, the process is complicated by the fact that the 
    ! right-hand-side vector {C} is not operated simultaniously with [A]. This means that the 
    ! computer algorithm must keep track of any row swiches that occur during the decomposition 
    ! step. Note that a vector, {o}, is used to keep track of the row interchanges.
    !
    !--------------------------------------------------------------------------
    !
    !   SizeLB,     SizeUB                    -> Allocated Lower Bound, Upper Bound
    !   WorkSizeLB, WorkSizeUB                -> Work Lower Bound, Upper Bound
    !   o      (0:nPropMax)                   -> Order vector
    !   Coefficients (0:nPropMax, 0:nPropMax) -> Matrix of coefficients
    !   IndTerm(0:nPropMax)                   -> Right-hand-side vector
    !   x      (0:nPropMax)                   -> Unknowns
    !   s      (0:nPropMax)                   -> Scale vector
    !
    !
    !   WARNING - ATTENTION TO KNOWN BUG - The coefficients matrix (Coefficients) 
    !                                      values are changed inside this subroutine.
    !                                      This means that every time the subroutine 
    !                                      is called the coefficients matrix must be
    !                                      re-filled with the values.
    !

    subroutine LUD(LUD_ID, Coefficients, IndTerm, x, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: LUD_ID
        double precision, pointer, dimension(: , :) :: Coefficients !Matrix of coefficients
        real,             pointer, dimension(:    ) :: IndTerm      !Right-hand-side vector
        real,             pointer, dimension(:    ) :: x            !Unknowns
        integer, optional, intent(OUT)              :: STAT

        !External---------------------------------------------------------------
        integer                                     :: ready_

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LUD_ID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            Me%s = 0.0
            Me%o = 0.0
            Me%x = 0.0

            Me%ExternalVar%Coefficients  => Coefficients
            if (.NOT. associated(Me%ExternalVar%Coefficients))                      &
                stop 'Subroutine LUD; Module ModuleLUD. ERR01' 

            
            Me%ExternalVar%IndTerm => IndTerm
            if (.NOT. associated(Me%ExternalVar%IndTerm))                           &
                stop 'Subroutine LUD; Module ModuleLUD. ERR02' 
           
            call order 
            call decomp
            call solve 

            x => Me%x       
            if (.NOT. associated(x))                                                &
                stop 'Subroutine LUD; Module ModuleLUD. ERR03' 

            nullify(Me%ExternalVar%Coefficients)
            nullify(Me%ExternalVar%IndTerm     )

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT)) &
            STAT = STAT_


        !----------------------------------------------------------------------

    end subroutine LUD

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine order

        !Local-----------------------------------------------------------------

        integer :: i, j

        !----------------------------------------------------------------------

do1 :   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            Me%o(i) = i
            Me%s(i) = abs(Me%ExternalVar%Coefficients(i, 1))
do2 :       do j = Me%WorkSize%ILB+1, Me%WorkSize%IUB
cd1 :           if (abs(Me%ExternalVar%Coefficients(i, j)) .GT. Me%s(i)) then
                    Me%s(i) = abs(Me%ExternalVar%Coefficients(i, j))
                end if cd1
            end do do2
        end do do1

        !----------------------------------------------------------------------

    end subroutine order

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine decomp

        !Local-----------------------------------------------------------------

        real :: sum

        integer :: i, j, k

        !----------------------------------------------------------------------

        j = 1
        call Pivot(j)
do1 :   do j = Me%WorkSize%ILB+1, Me%WorkSize%IUB
            Me%ExternalVar%Coefficients(Me%o(1), j) = Me%ExternalVar%Coefficients(Me%o(1), j) &
                                                    / Me%ExternalVar%Coefficients(Me%o(1), 1)
        end do do1

do2 :   do j = Me%WorkSize%ILB+1, Me%WorkSize%IUB-1
do6 :       do i = j, Me%WorkSize%IUB
                sum = 0.
do3 :           do k = 1, (j - 1)
                    sum = sum + Me%ExternalVar%Coefficients(Me%o(i), k)                       &
                              * Me%ExternalVar%Coefficients(Me%o(k), j)
                end do do3
                Me%ExternalVar%Coefficients(Me%o(i), j) = Me%ExternalVar%Coefficients(Me%o(i), j) - sum
            end do do6

            call Pivot(j)
do4 :       do k = j+1, Me%WorkSize%IUB
                sum = 0.
do5 :           do i = 1, (j - 1)
                    sum = sum + Me%ExternalVar%Coefficients(Me%o(j), i)                       &
                              * Me%ExternalVar%Coefficients(Me%o(i), k)
                end do do5
                Me%ExternalVar%Coefficients(Me%o(j), k) = (Me%ExternalVar%Coefficients(Me%o(j), k) - sum) &
                                                       /  Me%ExternalVar%Coefficients(Me%o(j), j)
            end do do4
        end do do2

        sum = 0.0
do7 :   do k = Me%WorkSize%ILB, Me%WorkSize%IUB-1
            sum = sum + Me%ExternalVar%Coefficients(Me%o(Me%WorkSize%IUB), k)                 &
                      * Me%ExternalVar%Coefficients(Me%o(k), Me%WorkSize%IUB)
        end do do7

        Me%ExternalVar%Coefficients(Me%o(Me%WorkSize%IUB), Me%WorkSize%IUB) =                 &
                    Me%ExternalVar%Coefficients(Me%o(Me%WorkSize%IUB), Me%WorkSize%IUB) - sum

        !----------------------------------------------------------------------

    end subroutine decomp

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine solve

        !Local-----------------------------------------------------------------

        real :: sum

        integer :: i, j

        !----------------------------------------------------------------------

        Me%x(1) = Me%ExternalVar%IndTerm (Me%o(1))                                  &
                / Me%ExternalVar%Coefficients(Me%o(1), 1)

do1 :   do i = Me%WorkSize%ILB+1, Me%WorkSize%IUB
            sum = 0.
do2 :       do j =  1, (i - 1)
                sum = sum + Me%ExternalVar%Coefficients(Me%o(i), j) * Me%x(j)
            end do do2

            Me%x(i) = (Me%ExternalVar%IndTerm (Me%o(i)) - sum)                      &
                    /  Me%ExternalVar%Coefficients(Me%o(i), i)
        end do do1

do3 :   do i = i-1, 1, -1
            sum = 0.
do4 :       do j = (i + 1), Me%WorkSize%IUB 
                sum = sum + Me%ExternalVar%Coefficients(Me%o(i), j) * Me%x(j)
            end do do4

            Me%x(i) = Me%x(i) - sum
        end do do3

        !----------------------------------------------------------------------
    
    end subroutine solve


    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------
    subroutine Pivot(j)

        !Arguments-------------------------------------------------------------

        integer, intent(IN) :: j

        !Local-----------------------------------------------------------------

        real :: big, dummy

        integer :: i, PivIt, iDum

        !----------------------------------------------------------------------

        PivIt = j
        big   = abs(Me%ExternalVar%Coefficients(Me%o(j), j)) / Me%s(Me%o(j))
do1 :   do i = j+1, Me%WorkSize%IUB
            dummy = abs(Me%ExternalVar%Coefficients(Me%o(i), j) / Me%s(Me%o(i)))
cd1 :       if (dummy .GT. big) then
                big   = dummy
                PivIt = i
            end if cd1
        end do do1

        iDum        = Me%o(PivIt)
        Me%o(PivIt) = Me%o(j)
        Me%o(j)     = iDum      

        !----------------------------------------------------------------------

    end subroutine Pivot



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillLUD(LUD_ID, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: LUD_ID
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                             :: STAT_, nUsers

        !----------------------------------------------------------------------                         


        STAT_ = UNKNOWN_


        call Ready(LUD_ID, ready_)
        
cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mLUD_,  Me%InstanceID)

            if (nUsers == 0) then

                deallocate(Me%o, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                  &
                    stop 'Subroutine KillLUD; module ModuleLUD. ERR03.'



                deallocate(Me%s, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                  &
                    stop 'Subroutine KillLUD; module ModuleLUD. ERR04.'



                deallocate(Me%x, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                  &
                    stop 'Subroutine KillLUD; module ModuleLUD. ERR01.'

                !Deallocates Instance
                call DeallocateInstance

                STAT_ = SUCCESS_

            end if

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine KillLUD

        !------------------------------------------------------------------------

    subroutine DeallocateInstance 


        !Local-----------------------------------------------------------------
        type (T_LUD), pointer                       :: AuxObjLUD
        type (T_LUD), pointer                       :: PreviousObjLUD

        !Updates pointers
        if (Me%InstanceID == FirstObjLUD%InstanceID) then
            FirstObjLUD => FirstObjLUD%Next
        else
            PreviousObjLUD => FirstObjLUD
            AuxObjLUD      => FirstObjLUD%Next
            do while (AuxObjLUD%InstanceID /= Me%InstanceID)
                PreviousObjLUD => AuxObjLUD
                AuxObjLUD      => AuxObjLUD%Next
            enddo

            !Now update linked list
            PreviousObjLUD%Next => AuxObjLUD%Next

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

    subroutine Ready (LUD_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: LUD_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (LUD_ID > 0) then
            call LocateObjLUD (LUD_ID)
            ready_ = VerifyReadLock (mLUD_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjLUD (LUD_ID)

        !Arguments-------------------------------------------------------------
        integer                                     :: LUD_ID

        !Local-----------------------------------------------------------------

        Me => FirstObjLUD
        do while (associated (Me))
            if (Me%InstanceID == LUD_ID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                                        &
            stop 'ModuleLUD - LocateObjLUD - ERR01'

    end subroutine LocateObjLUD

    !--------------------------------------------------------------------------

end module ModuleLUD

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------





