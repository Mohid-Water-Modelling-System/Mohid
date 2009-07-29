!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : HydroIntegration
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to integrate fluxes in time
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
!
!   Time convention (ex. 4 Intervals)
!
!
!       1   2   3   4   5   2a  3a      5x
!       |---|---|---|---|---|---| ... --|
!
!       +---+               ->  Integration%DT_ComputeStep
!
!       +---------------+   ->  Integration%DT
!
!
!       1 - Just called from the constructor
!               Sets Integration%InitialTime to 2 (CurrentTime + DT_ComputeStep)
!
!       2 - Reinitializes Integration (Stores VolumeZOld, which refers to 1)
!               Sets Integration%FinalTime to 5 (CurrentTime + Integration%DT - Integration%DT_ComputeStep)
!         - Integrates one Time Step
!
!       3 - Integrates one Time Step
!
!       4 - Integrates one Time Step
!
!       5 - Integrates one Time Step
!         - Finalizes Integration
!               Sets Integration%InitialTime to 2a (CurrentTime + DT_ComputeStep)
!


Module ModuleHydroIntegration

    use ModuleGlobalData
    use ModuleTime

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartHydroIntegration
    private :: AllocateInstance
    public  :: StartHydroIntegrationList
    private ::      NewIntegrationStep
    private ::          SearchForCurrentIntegration

    !Selector
    public  :: GetHydroIntegrationWaterFluxes
    public  :: GetHydroIntegrationComputeFaces
    public  :: GetHydroIntegrationDischarges
    public  :: GetHydroIntegrationVolumeZOld
    public  :: UnGetHydroIntegration
                     
    
    !Modifier
    public  :: ModifyHydroIntegration
    private ::      ReInitalizeIntegration
    private ::      OneIntegrationStep
    private ::      EndIntegrationStep


    !Destructor
    public  ::  KillHydroIntegration                                                     
    

    !Management
    private ::      Ready
    private ::          LocateObjHydroIntegration 
    !Interfaces----------------------------------------------------------------

    private :: UnGetHydroIntegration3D_I
    private :: UnGetHydroIntegration3D_R8
    interface  UnGetHydroIntegration
        module procedure UnGetHydroIntegration3D_I
        module procedure UnGetHydroIntegration3D_R8
    end interface  UnGetHydroIntegration

    !Types---------------------------------------------------------------------

    type T_OneIntegration
        real                                        :: DT
        real                                        :: DT_ComputeStep
        type (T_Time)                               :: InitialTime
        type (T_Time)                               :: FinalTime
        type (T_Time)                               :: LastActualization
        integer                                     :: n             !number of iterations in each integration time-step
        integer                                     :: Users
        type (T_Size3D)                             :: Size
        integer, dimension(:, :), pointer           :: BoundaryPoints2D

        real(8), dimension(:, :, :), pointer        :: InitialVolume

        real(8), dimension(:, :, :), pointer        :: WfluxX
        real(8), dimension(:, :, :), pointer        :: WfluxY
        real(8), dimension(:, :, :), pointer        :: WfluxZ

        real(8), dimension(:, :, :), pointer        :: Discharges

        integer, dimension(:, :, :), pointer        :: ComputeFacesU
        integer, dimension(:, :, :), pointer        :: ComputeFacesV
        integer, dimension(:, :, :), pointer        :: ComputeFacesW
        integer, dimension(:, :, :), pointer        :: OpenPoints3D

        type (T_OneIntegration), pointer            :: Next
    end type T_OneIntegration


    private :: T_HydroIntegration
    type       T_HydroIntegration
        integer                                     :: InstanceID
        type(T_OneIntegration  ), pointer           :: FirstIntegration

        type(T_HydroIntegration), pointer           :: Next
    end type  T_HydroIntegration

    !Global Module Variables
    type (T_HydroIntegration), pointer              :: FirstObjHydroIntegration
    type (T_HydroIntegration), pointer              :: Me
    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine StartHydroIntegration(ObjHydroIntegrationID, STAT)
                                     

        !Arguments---------------------------------------------------------------
        integer                                 :: ObjHydroIntegrationID
        integer, optional, intent(OUT)          :: STAT     

        !External----------------------------------------------------------------
        integer                                 :: ready_         

        !Local-------------------------------------------------------------------
        integer                                 :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mHydroIntegration_)) then
            nullify (FirstObjHydroIntegration)
            call RegisterModule (mHydroIntegration_) 
        endif

        call Ready(ObjHydroIntegrationID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            nullify(Me%FirstIntegration)

            !Returns ID
            ObjHydroIntegrationID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleHydroIntegration - StartHydroIntegration - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartHydroIntegration

    !--------------------------------------------------------------------------


    subroutine StartHydroIntegrationList(ObjHydroIntegrationID, Size,     &
                                         DT, CurrentTime, DT_ComputeStep, &
                                         BoundaryPoints2D, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: ObjHydroIntegrationID
        type (T_Size3D), intent(in)                 :: Size
        real, intent(IN)                            :: DT
        type(T_Time)                                :: CurrentTime
        real, intent(IN)                            :: DT_ComputeStep
        integer, dimension(:, :), pointer           :: BoundaryPoints2D
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        type (T_OneIntegration), pointer            :: CurrentIntegration
        integer                                     :: STAT_, ready_
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjHydroIntegrationID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !Searches for a integration with the current caracteristics
            nullify(CurrentIntegration)
            call SearchForCurrentIntegration(CurrentIntegration, DT, Size)
                                             
            !If there is no integration with the current caracteristics, create a new one
            if (.not. associated(CurrentIntegration)) then
                call NewIntegrationStep(CurrentIntegration, DT, Size, CurrentTime, &
                                        DT_ComputeStep, BoundaryPoints2D)                                        
            endif

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine StartHydroIntegrationList

    !--------------------------------------------------------------------------

    
    !--------------------------------------------------------------------------
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_HydroIntegration), pointer      :: NewObjHydroIntegration
        type (T_HydroIntegration), pointer      :: PreviousObjHydroIntegration


        !Allocates new instance
        allocate (NewObjHydroIntegration)
        nullify  (NewObjHydroIntegration%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjHydroIntegration)) then
            FirstObjHydroIntegration        => NewObjHydroIntegration
            Me                              => NewObjHydroIntegration
        else
            PreviousObjHydroIntegration     => FirstObjHydroIntegration
            Me                              => FirstObjHydroIntegration%Next
            do while (associated(Me))
                PreviousObjHydroIntegration => Me
                Me                          => Me%Next
            enddo
            Me                              => NewObjHydroIntegration
            PreviousObjHydroIntegration%Next=> NewObjHydroIntegration
        endif

        Me%InstanceID = RegisterNewInstance (mHYDROINTEGRATION_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------


    subroutine NewIntegrationStep(CurrentIntegration, DT, AuxSize,              &
                                  CurrentTime, DT_ComputeStep, BoundaryPoints2D)

        !Arguments-------------------------------------------------------------
        type(T_OneIntegration),   pointer           :: CurrentIntegration
        real                                        :: DT
        type(T_Size3D)                              :: AuxSize
        type(T_Time)                                :: CurrentTime
        real                                        :: DT_ComputeStep
        integer, dimension(:, :), pointer           :: BoundaryPoints2D

        !Local-----------------------------------------------------------------
        type(T_OneIntegration), pointer             :: PointerIntegration
        type(T_OneIntegration), pointer             :: PreviousIntegration
        integer                                     :: ILB, IUB
        integer                                     :: JLB, JUB
        integer                                     :: KLB, KUB
        integer                                     :: i, j   

        !Shorten Variables Names
        ILB = AuxSize%ILB - 1
        IUB = AuxSize%IUB + 1

        JLB = AuxSize%JLB - 1
        JUB = AuxSize%JUB + 1
        
        KLB = AuxSize%KLB - 1
        KUB = AuxSize%KUB + 1

        !Allocates new integration
        allocate(CurrentIntegration)
        nullify (CurrentIntegration%Next)

        !Allocates Variables of this integration
        allocate(CurrentIntegration%BoundaryPoints2D (ILB:IUB, JLB:JUB         ))
        allocate(CurrentIntegration%InitialVolume    (ILB:IUB, JLB:JUB, KLB:KUB))

        allocate(CurrentIntegration%WfluxX           (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(CurrentIntegration%WfluxY           (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(CurrentIntegration%WfluxZ           (ILB:IUB, JLB:JUB, KLB:KUB))

        allocate(CurrentIntegration%Discharges       (ILB:IUB, JLB:JUB, KLB:KUB))

        allocate(CurrentIntegration%ComputeFacesU    (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(CurrentIntegration%ComputeFacesV    (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(CurrentIntegration%ComputeFacesW    (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(CurrentIntegration%OpenPoints3D     (ILB:IUB, JLB:JUB, KLB:KUB))


        !Initializes data
        CurrentIntegration%DT            = DT
        CurrentIntegration%DT_ComputeStep= DT_ComputeStep
        CurrentIntegration%Size          = AuxSize
        CurrentIntegration%Users         = 1
        CurrentIntegration%n             = 0
        CurrentIntegration%InitialTime   = CurrentTime + DT_ComputeStep
        CurrentIntegration%FinalTime     = CurrentTime                              

        CurrentIntegration%InitialVolume = null_real

        CurrentIntegration%WfluxX        = 0.
        CurrentIntegration%WfluxY        = 0.
        CurrentIntegration%WfluxZ        = 0.

        CurrentIntegration%Discharges    = 0.

        CurrentIntegration%ComputeFacesU = 0
        CurrentIntegration%ComputeFacesV = 0
        CurrentIntegration%ComputeFacesW = 0
        CurrentIntegration%OpenPoints3D  = 0

        do j = JLB, JUB
        do i = ILB, IUB
            CurrentIntegration%BoundaryPoints2D(i, j) = BoundaryPoints2D(i, j)
        enddo
        enddo
            


        !Insert new Integration into the list
        if (.not. associated(Me%FirstIntegration)) then
            Me%FirstIntegration => CurrentIntegration
        else
            PointerIntegration => Me%FirstIntegration
            do while (associated(PointerIntegration))
                PreviousIntegration => PointerIntegration
                PointerIntegration  => PointerIntegration%Next
            enddo
            PreviousIntegration%Next => CurrentIntegration
        endif

    end subroutine NewIntegrationStep        


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine GetHydroIntegrationWaterFluxes(ObjHydroIntegrationID, DT, Size,  &
                                              WaterFluxX, WaterFluxY,           &
                                              WaterFluxZ, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjHydroIntegrationID
        real                                            :: DT
        type (T_Size3D)                                 :: Size
        real(8), dimension(:, :, :), optional, pointer  :: WaterFluxX
        real(8), dimension(:, :, :), optional, pointer  :: WaterFluxY
        real(8), dimension(:, :, :), optional, pointer  :: WaterFluxZ
        integer, intent(OUT), optional                  :: STAT


        !Local-----------------------------------------------------------------
        type (T_OneIntegration), pointer                :: CurrentIntegration
        integer                                         :: STAT_, ready_
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjHydroIntegrationID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !Searches for a integration with the current caracteristics
            nullify(CurrentIntegration)
            call SearchForCurrentIntegration(CurrentIntegration, DT, Size)
                                             

            !If there is no integration with the current caracteristics, create a new one
            if (.not. associated(CurrentIntegration))                                    &
                call SetError(FATAL_, INTERNAL_, "GetHydroIntegrationWaterFluxes - ModuleHydroIntegration - ERR01")


            if (present(WaterFluxX)) then
                call Read_Lock(mHYDROINTEGRATION_, Me%InstanceID)
                WaterFluxX => CurrentIntegration%WfluxX
            endif

            if (present(WaterFluxY)) then
                call Read_Lock(mHYDROINTEGRATION_, Me%InstanceID)
                WaterFluxY => CurrentIntegration%WfluxY
            endif

            if (present(WaterFluxZ)) then
                call Read_Lock(mHYDROINTEGRATION_, Me%InstanceID)
                WaterFluxZ => CurrentIntegration%WfluxZ
            endif

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetHydroIntegrationWaterFluxes

    !--------------------------------------------------------------------------

    subroutine GetHydroIntegrationComputeFaces(ObjHydroIntegrationID, DT, Size,            &
                                               ComputeFacesU,ComputeFacesV, ComputeFacesW, &
                                               OpenPoints3D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjHydroIntegrationID
        real                                            :: DT
        type (T_Size3D)                                 :: Size
        integer, dimension(:, :, :), optional, pointer  :: ComputeFacesU
        integer, dimension(:, :, :), optional, pointer  :: ComputeFacesV
        integer, dimension(:, :, :), optional, pointer  :: ComputeFacesW
        integer, dimension(:, :, :), optional, pointer  :: OpenPoints3D
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        type (T_OneIntegration), pointer                :: CurrentIntegration
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_


        call Ready(ObjHydroIntegrationID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !Searches for a integration with the current caracteristics
            nullify(CurrentIntegration)
            call SearchForCurrentIntegration(CurrentIntegration, DT, Size)
                                             

            !If there is no integration with the current caracteristics, create a new one
            if (.not. associated(CurrentIntegration))                                    &
                call SetError(FATAL_, INTERNAL_, "GetHydroIntegrationWaterFluxes - ModuleHydroIntegration - ERR01")


            if (present(ComputeFacesU)) then
                call Read_Lock(mHYDROINTEGRATION_, Me%InstanceID)
                ComputeFacesU => CurrentIntegration%ComputeFacesU
            endif

            if (present(ComputeFacesV)) then
                call Read_Lock(mHYDROINTEGRATION_, Me%InstanceID)
                ComputeFacesV => CurrentIntegration%ComputeFacesV
            endif

            if (present(ComputeFacesW)) then
                call Read_Lock(mHYDROINTEGRATION_, Me%InstanceID)
                ComputeFacesW => CurrentIntegration%ComputeFacesW
            endif

            if (present(OpenPoints3D)) then
                call Read_Lock(mHYDROINTEGRATION_, Me%InstanceID)
                OpenPoints3D => CurrentIntegration%OpenPoints3D
            endif

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetHydroIntegrationComputeFaces

    !--------------------------------------------------------------------------

    subroutine GetHydroIntegrationVolumeZOld (ObjHydroIntegrationID, DT, Size, &
                                              VolumeZOld, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjHydroIntegrationID
        real                                            :: DT
        type (T_Size3D)                                 :: Size
        real(8), dimension(:, :, :), optional, pointer  :: VolumeZOld
        integer, intent(OUT), optional                  :: STAT


        !Local-----------------------------------------------------------------
        type (T_OneIntegration), pointer                :: CurrentIntegration
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjHydroIntegrationID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            !Searches for a integration with the current caracteristics
            nullify(CurrentIntegration)
            call SearchForCurrentIntegration(CurrentIntegration, DT, Size)
                                             

            !If there is no integration with the current caracteristics, create a new one
            if (.not. associated(CurrentIntegration))                                    &
                call SetError(FATAL_, INTERNAL_, "GetHydroIntegrationWaterFluxes - ModuleHydroIntegration - ERR01")


            if (present(VolumezOld)) then
                call Read_Lock(mHYDROINTEGRATION_, Me%InstanceID)
                VolumeZOld => CurrentIntegration%InitialVolume
            endif

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetHydroIntegrationVolumeZOld

    !--------------------------------------------------------------------------

    subroutine GetHydroIntegrationDischarges (ObjHydroIntegrationID, DT, Size, Discharges, STAT)
                                             

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjHydroIntegrationID
        real                                            :: DT
        type (T_Size3D)                                 :: Size
        real(8), dimension(:, :, :), optional, pointer  :: Discharges
        integer, intent(OUT), optional                  :: STAT


        !Local-----------------------------------------------------------------
        type (T_OneIntegration), pointer                :: CurrentIntegration
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjHydroIntegrationID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !Searches for a integration with the current caracteristics
            nullify(CurrentIntegration)
            call SearchForCurrentIntegration(CurrentIntegration, DT, Size)

            !If there is no integration with the current caracteristics, create a new one
            if (.not. associated(CurrentIntegration))                                    &
                call SetError(FATAL_, INTERNAL_, "GetHydroIntegrationDischarges - ModuleHydroIntegration - ERR01")


            if (present(Discharges)) then
                call Read_Lock(mHYDROINTEGRATION_, Me%InstanceID)
                Discharges => CurrentIntegration%Discharges
            endif

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetHydroIntegrationDischarges

    !--------------------------------------------------------------------------

    subroutine UnGetHydroIntegration3D_I(ObjHydroIntegrationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjHydroIntegrationID
        integer, dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjHydroIntegrationID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mHYDROINTEGRATION_, Me%InstanceID, "UnGetHydroIntegration3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetHydroIntegration3D_I

    !--------------------------------------------------------------------------

    subroutine UnGetHydroIntegration3D_R8(ObjHydroIntegrationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjHydroIntegrationID
        real(8), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjHydroIntegrationID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mHYDROINTEGRATION_, Me%InstanceID, "UnGetHydroIntegration3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetHydroIntegration3D_R8

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyHydroIntegration(ObjHydroIntegrationID, WaterFluxX, WaterFluxY,     &
                                      Discharges, ComputeFacesU, ComputeFacesV,          &
                                      WaterPoints3D, VolumeZ, VolumeZOld, CurrentTime, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjHydroIntegrationID
        real(8), dimension(:, :, :), pointer        :: WaterFluxX, WaterFluxY
        real(8), dimension(:, :, :), pointer        :: Discharges
        integer, dimension(:, :, :), pointer        :: ComputeFacesU, ComputeFacesV
        integer, dimension(:, :, :), pointer        :: WaterPoints3D
        real(8), dimension(:, :, :), pointer        :: VolumeZ, VolumeZOld
        type(T_Time)                                :: CurrentTime
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        type (T_OneIntegration), pointer            :: CurrentIntegration
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjHydroIntegrationID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !Steps through the integration list
            CurrentIntegration => Me%FirstIntegration
            do while (associated(CurrentIntegration))


                if (CurrentTime == CurrentIntegration%InitialTime) then
                    call ReInitalizeIntegration(CurrentIntegration, CurrentTime, VolumeZOld)
                endif

                call OneIntegrationStep(CurrentIntegration, WaterFluxX, WaterFluxY, Discharges, ComputeFacesU,               &
                                        ComputeFacesV)

                if (CurrentTime == CurrentIntegration%FinalTime) then
                    call EndIntegrationStep(CurrentIntegration, VolumeZ, WaterPoints3D)
                endif
            
                CurrentIntegration => CurrentIntegration%Next
            enddo

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyHydroIntegration

    !--------------------------------------------------------------------------

    subroutine ReInitalizeIntegration(CurrentIntegration, CurrentTime, VolumeZOld)
            
        !Arguments-------------------------------------------------------------
        type(T_OneIntegration), pointer             :: CurrentIntegration
        type(T_Time)                                :: CurrentTime
        real(8), dimension(:, :, :), pointer        :: VolumeZOld

        !Local-----------------------------------------------------------------

        CurrentIntegration%Users         = 1
        CurrentIntegration%n             = 0

        CurrentIntegration%InitialVolume = VolumeZOld

        CurrentIntegration%WfluxX        = 0.
        CurrentIntegration%WfluxY        = 0.
        CurrentIntegration%WfluxZ        = 0.

        CurrentIntegration%Discharges    = 0.

        CurrentIntegration%ComputeFacesU = 0
        CurrentIntegration%ComputeFacesV = 0
        CurrentIntegration%ComputeFacesW = 0
        CurrentIntegration%OpenPoints3D  = 0

        CurrentIntegration%FinalTime = CurrentTime +                                     &
                                      (CurrentIntegration%DT -                           &
                                       CurrentIntegration%DT_ComputeStep)

    end subroutine ReInitalizeIntegration


    !--------------------------------------------------------------------------

    subroutine SearchForCurrentIntegration(CurrentIntegration, DT, AuxSize)
                                           

        !Arguments-------------------------------------------------------------
        type(T_OneIntegration),   pointer           :: CurrentIntegration
        real                                        :: DT
        type(T_Size3D)                              :: AuxSize

        !Local-----------------------------------------------------------------
        type(T_OneIntegration), pointer             :: PointerIntegration

        
        PointerIntegration => Me%FirstIntegration
        do while (associated(PointerIntegration))
            
            if (PointerIntegration%DT       == DT            .and.                       &
                PointerIntegration%Size%ILB == AuxSize%ILB   .and.                       &
                PointerIntegration%Size%IUB == AuxSize%IUB   .and.                       &
                PointerIntegration%Size%JLB == AuxSize%JLB   .and.                       &
                PointerIntegration%Size%JUB == AuxSize%JUB   .and.                       &
                PointerIntegration%Size%KLB == AuxSize%KLB   .and.                       &
                PointerIntegration%Size%KUB == AuxSize%KUB) then

                CurrentIntegration => PointerIntegration
                CurrentIntegration%Users = CurrentIntegration%Users + 1
                exit

            else

                PointerIntegration => PointerIntegration%Next

            endif
        
        enddo

    end subroutine SearchForCurrentIntegration

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine OneIntegrationStep(CurrentIntegration, WaterFluxX, WaterFluxY,            &
                                  Discharges, ComputeFacesU, ComputeFacesV)

        !Arguments-------------------------------------------------------------
        type(T_OneIntegration),   pointer           :: CurrentIntegration
        real(8), dimension(:, :, :), pointer        :: WaterFluxX
        real(8), dimension(:, :, :), pointer        :: WaterFluxY
        real(8), dimension(:, :, :), pointer        :: Discharges
        integer, dimension(:, :, :), pointer        :: ComputeFacesU
        integer, dimension(:, :, :), pointer        :: ComputeFacesV

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB
        integer                                     :: JLB, JUB
        integer                                     :: KLB, KUB
        integer                                     :: i, j, k
        integer                                     :: n

        !Shorten Variables Names
        ILB = CurrentIntegration%Size%ILB
        IUB = CurrentIntegration%Size%IUB

        JLB = CurrentIntegration%Size%JLB
        JUB = CurrentIntegration%Size%JUB
        
        KLB = CurrentIntegration%Size%KLB
        KUB = CurrentIntegration%Size%KUB


        CurrentIntegration%n = CurrentIntegration%n + 1
        n                    = CurrentIntegration%n

        do j = JLB, JUB + 1
        do i = ILB, IUB + 1
        do k = KLB, KUB
            CurrentIntegration%WfluxX(i, j, k) = (CurrentIntegration%WfluxX(i, j, k)  &
                                                 * (n-1) + WaterFluxX(i, j, k)) / n       !/ AverageFactor
            CurrentIntegration%WfluxY(i, j, k) = (CurrentIntegration%WfluxY(i, j, k)  &
                                                 * (n-1) + WaterFluxY(i, j, k)) / n       !/ AverageFactor
        end do
        end do
        end do

        do j = JLB, JUB
        do i = ILB, IUB
        do k = KLB, KUB
            CurrentIntegration%Discharges(i, j, k) = (CurrentIntegration%Discharges(i,j,k) &
                                                   * (n-1) + Discharges(i,j,k)) / n     !/ AverageFactor
        end do
        end do
        end do

        !Integration of Mapping
        do j = JLB, JUB + 1
        do i = ILB, IUB + 1
        do k = KLB, KUB
            if (ComputeFacesU(i, j, k) > 0) CurrentIntegration%ComputeFacesU(i, j, k) = 1
            if (ComputeFacesV(i, j, k) > 0) CurrentIntegration%ComputeFacesV(i, j, k) = 1
        end do
        end do
        end do


    end subroutine OneIntegrationStep

    !--------------------------------------------------------------------------

    subroutine EndIntegrationStep(CurrentIntegration, VolumeZ, WaterPoints3D)

        !Arguments-------------------------------------------------------------
        type(T_OneIntegration),   pointer           :: CurrentIntegration
        real(8), dimension(:, :, :), pointer        :: VolumeZ
        integer, dimension(:, :, :), pointer        :: WaterPoints3D

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB
        integer                                     :: JLB, JUB
        integer                                     :: KLB, KUB
        integer                                     :: i, j, k
        real(8)                                     :: dVdt

        !Shorten Variables Names
        ILB = CurrentIntegration%Size%ILB
        IUB = CurrentIntegration%Size%IUB

        JLB = CurrentIntegration%Size%JLB
        JUB = CurrentIntegration%Size%JUB
        
        KLB = CurrentIntegration%Size%KLB
        KUB = CurrentIntegration%Size%KUB

            
        !Fluxes divergence
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            
            dVdt = (VolumeZ(i, j, k) - CurrentIntegration%InitialVolume(i, j, k)) /  &
                    dble(CurrentIntegration%DT)

            CurrentIntegration%WfluxZ(i, j, k+1) =                               &
                        ( CurrentIntegration%WfluxZ(i,  j,  k)                   &  !Bottom face
                        + CurrentIntegration%WfluxX(i,  j,  k)                   &  !West   face
                        - CurrentIntegration%WfluxX(i,  j+1,k)                   &  !East   face 
                        + CurrentIntegration%WfluxY(i,  j,  k)                   &  !South  face
                        - CurrentIntegration%WfluxY(i+1,j,  k)                   &  !North  face 
                        - dVdt + CurrentIntegration%Discharges(i, j, k))         &
                        * (1.0 - CurrentIntegration%BoundaryPoints2D(i, j))

        enddo
        enddo
        enddo


        !ComputeFacesW
        do j = JLB, JUB
        do i = ILB, IUB
            if (CurrentIntegration%ComputeFacesU(i, j,   KUB) +                        &
                CurrentIntegration%ComputeFacesU(i, j+1, KUB) +                        &
                CurrentIntegration%ComputeFacesV(i, j,   KUB) +                        &
                CurrentIntegration%ComputeFacesV(i+1, j, KUB) > 0) then
                                    
                do k = KLB+1, KUB
                    if (WaterPoints3D(i, j, k-1) == WaterPoint)                         &
                        CurrentIntegration%ComputeFacesW(i, j, k) = 1
                enddo

            endif

        enddo
        enddo


        !OpenPoints
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            if (CurrentIntegration%ComputeFacesU(i, j, k)   +                        &
                CurrentIntegration%ComputeFacesU(i, j+1, k) +                        & 
                CurrentIntegration%ComputeFacesV(i, j, k)   +                        &
                CurrentIntegration%ComputeFacesV(i+1, j, k) +                        & 
                CurrentIntegration%ComputeFacesW(i, j, k)   +                        &
                CurrentIntegration%ComputeFacesW(i, j, k+1) > 0 )                    &
                CurrentIntegration%OpenPoints3D(i, j, k) = 1
        enddo
        enddo
        enddo

        CurrentIntegration%InitialTime = CurrentIntegration%FinalTime +                  &
                                         CurrentIntegration%DT_ComputeStep

    end subroutine EndIntegrationStep

    !----------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillHydroIntegration(ObjHydroIntegrationID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjHydroIntegrationID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_, nUsers              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_          

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjHydroIntegrationID, ready_)    


cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mHYDROINTEGRATION_,  Me%InstanceID)

            if (nUsers == 0) then

                !Deallocates Instance
                call DeallocateInstance

                ObjHydroIntegrationID = 0

                STAT_ = SUCCESS_

            end if

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                      &
            STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine KillHydroIntegration
    
    
    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_HydroIntegration), pointer          :: AuxObjHydroIntegration
        type (T_HydroIntegration), pointer          :: PreviousObjHydroIntegration

        !Updates pointers
        if (Me%InstanceID == FirstObjHydroIntegration%InstanceID) then
            FirstObjHydroIntegration => FirstObjHydroIntegration%Next
        else
            PreviousObjHydroIntegration => FirstObjHydroIntegration
            AuxObjHydroIntegration      => FirstObjHydroIntegration%Next
            do while (AuxObjHydroIntegration%InstanceID /= Me%InstanceID)
                PreviousObjHydroIntegration => AuxObjHydroIntegration
                AuxObjHydroIntegration      => AuxObjHydroIntegration%Next
            enddo

            !Now update linked list
            PreviousObjHydroIntegration%Next => AuxObjHydroIntegration%Next

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

    subroutine Ready (ObjHydroIntegration_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjHydroIntegration_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

        
cd1:    if (ObjHydroIntegration_ID > 0) then
            call LocateObjHydroIntegration (ObjHydroIntegration_ID)
            ready_ = VerifyReadLock (mHYDROINTEGRATION_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjHydroIntegration (ObjHydroIntegrationID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjHydroIntegrationID

        !Local-----------------------------------------------------------------

        Me => FirstObjHydroIntegration
        do while (associated (Me))
            if (Me%InstanceID == ObjHydroIntegrationID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                          &
            stop 'ModuleHydroIntegration - LocateObjHydroIntegration - ERR01'

    end subroutine LocateObjHydroIntegration

    !--------------------------------------------------------------------------

end module ModuleHydroIntegration

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------






