!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Profile
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : April 2005
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module output vertical profile results in HDF5
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
!DataFile
!   PROFILE_OUTPUT_TIME         : real             [-]          !Time step to perform profile output
!<begin_profile>
    !NAME                       : char              -           !Name of the profile output time
    !LOCALIZATION_I             : int               -           !Profile localization I    
    !LOCALIZATION_J             : int               -           !Profile localization J
!<end_profile>


Module ModuleProfile

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleHDF5

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartProfile
    private ::      AllocateInstance
    private ::      ConstructProfileList
    private ::          AddProfile
    private ::          ConstructProfile
    private ::      VerifyProfileLocation

    !Selector
    public  :: GetProfileNextOutputTime                 
    
    !Modifier
    public  :: WriteProfile

    !Destructor
    public  :: KillProfile
    private :: WriteFinalFile
    private ::      DeallocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjProfile 
    
    !Interfaces----------------------------------------------------------------
    private :: WriteProfile_R4
    private :: WriteProfile_R8
    interface  WriteProfile
        module procedure WriteProfile_R4
        module procedure WriteProfile_R8
    end interface WriteProfile

    !Types---------------------------------------------------------------------
    
    type      T_ExternalVar
        real, dimension(:,:,:), pointer                     :: SZZ
        type(T_Time)                                        :: Now
    end type  T_ExternalVar

    type      T_ID
        integer                                             :: ID               = null_int
        character(len=PathLength)                           :: Name
    end type  T_ID

    type       T_OneProfile
        type(T_ID)                                          :: ID
        integer                                             :: ObjHDF5          = 0
        integer                                             :: LocalizationI    = null_int
        integer                                             :: LocalizationJ    = null_int
        real(8), dimension(:,:,:), pointer                  :: Values
        real(8), dimension(:,:  ),  pointer                 :: VerticalZ
        type(T_OneProfile), pointer                         :: Next
    end type  T_OneProfile

    type       T_Profile
        integer                                             :: InstanceID
        type (T_Size3D)                                     :: Size, WorkSize
        integer                                             :: ObjTime          = 0
        integer                                             :: ObjEnterData     = 0
        integer                                             :: KUB
        integer                                             :: PropertyCount
        integer                                             :: NumberOfProperties
        integer                                             :: OutputNumber
        integer                                             :: nInstants
        type (T_Time), dimension(:  ),  pointer             :: OutTime
        real(8),       dimension(:,:),  pointer             :: Time2D
        real(8),       dimension(:,:),  pointer             :: Depth
        character(len=StringLength)                         :: ClientName
        character(len=StringLength), dimension(:,:), pointer:: PropertyList
        character(len=PathLength)                           :: RootPath
        type(T_ExternalVar)                                 :: ExternalVar
        type(T_OneProfile),             pointer             :: FirstProfile
        type(T_Profile),                pointer             :: Next
    end type  T_Profile

    !Global Module Variables
    type (T_Profile), pointer                               :: FirstObjProfile
    type (T_Profile), pointer                               :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartProfile(ProfileID, ObjTime, ProfileDataFile, WaterPoints2D, &
                            nProperties, KUB, PropertyList, ClientName,         &
                            VerifyLocation, OutTime, STAT)

        !Arguments---------------------------------------------------------------
        integer                                             :: ProfileID
        integer                                             :: ObjTime
        character(len=*)                                    :: ProfileDataFile
        integer,                    dimension(:,:), pointer :: WaterPoints2D
        integer                                             :: KUB, nProperties
        character(len=StringLength),dimension(:,:), pointer :: PropertyList
        character(len=*)                                    :: ClientName
        logical,        optional                            :: VerifyLocation
        type (T_Time),  optional,   dimension(:  ), pointer :: OutTime
        integer,        optional,   intent(OUT)             :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL        

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_
        logical                                         :: VerifyLocation_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mProfile_)) then
            nullify (FirstObjProfile)
            call RegisterModule (mProfile_) 
        endif

        call Ready(ProfileID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            !Associates module Time
            Me%ObjTime                   =  AssociateInstance (mTIME_, ObjTime)
            Me%KUB                       =  KUB
            Me%NumberOfProperties        =  nProperties
            Me%ClientName                =  ClientName

            if(present(VerifyLocation))then
                VerifyLocation_ = VerifyLocation
            else
                VerifyLocation_ = .true.
            end if

            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)   
            if (STAT_CALL .NE. SUCCESS_) stop 'StartProfile - ModuleProfile - ERR00'

            !Constructs EnterData
            call ConstructEnterData(Me%ObjEnterData, ProfileDataFile, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartProfile - ModuleProfile - ERR01'

            call ReadGlobalOptions(OutTime)

            call ConstructPropertyList(PropertyList)

            call ConstructProfileList

            if(VerifyLocation_)then
                call VerifyProfileLocation(WaterPoints2D)
            endif

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartProfile - ModuleProfile - ERR02'

            !Returns ID
            ProfileID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleProfile - StartProfile - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartProfile
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Profile), pointer                         :: NewObjProfile
        type (T_Profile), pointer                         :: PreviousObjProfile


        !Allocates new instance
        allocate (NewObjProfile)
        nullify  (NewObjProfile%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjProfile)) then
            FirstObjProfile         => NewObjProfile
            Me                    => NewObjProfile
        else
            PreviousObjProfile    => FirstObjProfile
            Me                    => FirstObjProfile%Next
            do while (associated(Me))
                PreviousObjProfile  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjProfile
            PreviousObjProfile%Next => NewObjProfile
        endif

        Me%InstanceID = RegisterNewInstance (mProfile_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ReadGlobalOptions(OutTime)
        
        !Arguments-------------------------------------------------------------
        type (T_Time),  optional,   dimension(:  ), pointer :: OutTime

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        logical                                             :: FoundDT
        type(T_Time)                                        :: BeginTime, EndTime

        !Begin-----------------------------------------------------------------

        Me%OutputNumber = 1

        call GetComputeTimeLimits(Me%ObjTime,                                   &
                                  BeginTime = BeginTime,                        &
                                  EndTime   = EndTime,                          &
                                  STAT      = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             &
            stop 'ReadGlobalOptions - ModuleProfile - ERR00'

        call GetOutPutTime(Me%ObjEnterData,                                     &
                           CurrentTime = Me%ExternalVar%Now,                    &
                           EndTime     = EndTime,                               &
                           keyword     = 'PROFILE_OUTPUT_TIME',                 &
                           SearchType  = FromFile,                              &
                           OutPutsTime = Me%OutTime,                            &
                           OutPutsOn   = FoundDT,                               &
                           STAT        = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             &
            stop 'ReadGlobalOptions - ModuleProfile - ERR01'

        if (.not. FoundDT) then

            write(*,*)'Trying to construct vertical profiles output'
            write(*,*)'but keyword PROFILE_OUTPUT_TIME is missing.'
            write(*,*)'Client module/program    : ', trim(Me%ClientName)

            if(present(OutTime))then

                write(*,*)'Assuming output times set by client module'
                
                Me%nInstants = size(OutTime)

                allocate(Me%OutTime(1:Me%nInstants))
                
                Me%OutTime = OutTime

            else
                stop 'ReadGlobalOptions - ModuleProfile - ERR02'
            end if 

        else

            Me%nInstants = size(Me%OutTime)

        endif

        
        allocate(Me%Depth    (1:Me%KUB, 1:Me%nInstants))
        allocate(Me%Time2D   (1:Me%KUB, 1:Me%nInstants))

        !Gets the root path from the file nomfich.dat
        call ReadFileName("ROOT_SRT", Me%RootPath, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            call ReadFileName("ROOT", Me%RootPath, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                call ReadFileName("RAIZ", Me%RootPath, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleProfile - ERR04'
            endif
        endif

    end subroutine ReadGlobalOptions

    !--------------------------------------------------------------------------

    subroutine ConstructPropertyList(PropertyList)
       
        !Arguments-------------------------------------------------------------
        character(len=StringLength), dimension(:,:), pointer:: PropertyList
        
        !Begin-----------------------------------------------------------------
        integer                                             :: u,n
        !Begin-----------------------------------------------------------------
        
        allocate(Me%PropertyList(1:Me%NumberOfProperties, 1:2))

        do n = 1, Me%NumberOfProperties
        do u = 1, 2
            Me%PropertyList(n,u) =  PropertyList(n,u)
        enddo
        enddo


    end subroutine ConstructPropertyList
    
    !--------------------------------------------------------------------------

    subroutine ConstructProfileList
        
        !Local-----------------------------------------------------------------
        type (T_OneProfile),    pointer                     :: NewProfile
        integer                                             :: ClientNumber, STAT_CALL
        logical                                             :: BlockFound

        !Begin-----------------------------------------------------------------


        do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<begin_profile>',    &
                                        block_end       = '<end_profile>',      &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
            if(STAT_CALL .EQ. SUCCESS_)then
            
                if (BlockFound) then                                                  
                
                    call AddProfile          (NewProfile)

                    call ConstructProfile    (NewProfile)

                    nullify(NewProfile)

                else 
                    
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop       'ConstructProfileList - ModuleProfile - ERR01'

                    exit 

                end if 

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                    stop       'ConstructProfileList - ModuleProfile - ERR02'
            else 
                    stop       'ConstructProfileList - ModuleProfile - ERR03'
            end if 
        end do 

    end subroutine ConstructProfileList

    !--------------------------------------------------------------------------

    subroutine AddProfile (ObjProfile)

        !Arguments-------------------------------------------------------------
        type (T_OneProfile),      pointer           :: ObjProfile
        
        !Local-----------------------------------------------------------------
        type (T_OneProfile),      pointer           :: PreviousProfile
        type (T_OneProfile),      pointer           :: NewProfile
        integer, save                               :: NextProfileID = 1
        
        !Begin-----------------------------------------------------------------

        !Allocates new Profile
        allocate (NewProfile)
        nullify  (NewProfile%Next)

        !Insert new Profile into list and makes current algae point to it
        if (.not. associated(Me%FirstProfile)) then
            Me%FirstProfile            => NewProfile
            ObjProfile                 => NewProfile
        else
            PreviousProfile            => Me%FirstProfile
            ObjProfile                 => Me%FirstProfile%Next

            do while (associated(ObjProfile))
                PreviousProfile        => ObjProfile
                ObjProfile             => ObjProfile%Next
            enddo
            ObjProfile                 => NewProfile
            PreviousProfile%Next       => NewProfile
        endif

        !Attributes ID
        ObjProfile%ID%ID               = NextProfileID

        NextProfileID                  = NextProfileID + 1


    end subroutine AddProfile

    !--------------------------------------------------------------------------

    subroutine ConstructProfile(NewProfile)

        !Arguments-------------------------------------------------------------
        type (T_OneProfile),      pointer                   :: NewProfile
 
        !Local-----------------------------------------------------------------
        integer                                             :: iflag, STAT_CALL
        character(len=3)                                    :: AuxI, AuxJ
        character(len=StringLength)                         :: AuxName

        !Begin-----------------------------------------------------------------


        call GetData(NewProfile%LocalizationI,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     KeyWord      = 'LOCALIZATION_I',                                   &
                     SearchType   = FromBlock,                                          &
                     ClientModule = 'ModuleProfile',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfile - ModuleProfile - ERR01' 
        
        
        call GetData(NewProfile%LocalizationJ,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     KeyWord      = 'LOCALIZATION_J',                                   &
                     SearchType   = FromBlock,                                          &
                     ClientModule = 'ModuleProfile',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfile - ModuleProfile - ERR02' 

        !Constructs the default name of the file
        AuxI = "   "
        AuxJ = "   "
        write(AuxI,'(i3)')NewProfile%LocalizationI
        write(AuxJ,'(i3)')NewProfile%LocalizationJ


        call GetData(AuxName,                                                           &
                     Me%ObjEnterData, iflag,                                            &
                     KeyWord      = 'NAME',                                             &
                     SearchType   = FromBlock,                                          &
                     ClientModule = 'ModuleProfile',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfile - ModuleProfile - ERR03'


        if(iflag == 1) then
            NewProfile%ID%Name = trim(adjustl(Me%RootPath))  //                         &
                                 trim(adjustl(Me%ClientName))//"_"//                    &
                                 trim(adjustl(AuxName))//".hdf5"
        else
            NewProfile%ID%Name = trim(adjustl(Me%RootPath))     //                      &
                                 trim(adjustl(Me%ClientName))   //"_"//                 &
                                 trim(adjustl(AuxI))//"_"       //                      &
                                 trim(adjustl(AuxJ))//".hdf5"

        end if

        allocate(NewProfile%Values   (1:Me%KUB, 1:Me%nInstants, 1:Me%NumberOfProperties))
        allocate(NewProfile%VerticalZ(0:Me%KUB, 1:Me%nInstants))
    
    end subroutine ConstructProfile

    !--------------------------------------------------------------------------

    subroutine VerifyProfileLocation(WaterPoints2D)

        !Arguments-------------------------------------------------------------
        integer, dimension(:,:), pointer        :: WaterPoints2D

        !Local-----------------------------------------------------------------
        integer                                 :: i, j
        logical                                 :: AllTestsPassed
        type(T_OneProfile), pointer             :: Profile
        !----------------------------------------------------------------------

        AllTestsPassed = .true.

        Profile => Me%FirstProfile

        do while(associated(Profile))

            i = Profile%LocalizationI
            j = Profile%LocalizationJ

            !Checks the bounds of the matrix
            if (i < lbound(WaterPoints2D, dim = 1) .or. &
                i > ubound(WaterPoints2D, dim = 1) .or. &
                j < lbound(WaterPoints2D, dim = 2) .or. &
                j > ubound(WaterPoints2D, dim = 2)) then
                write(*,*)"A profile output was defined in the following point,"
                write(*,*)"meanwhile it is outside the domain"
                write(*,*)"I = ", i
                write(*,*)"J = ", j
                write(*,*)"Please review the data files"
                AllTestsPassed = .false.
            endif

            if (.not. AllTestsPassed) stop 'VerifyProfileLocation - ModuleProfile - ERR01'

            !Checks the Waterpoints
            if (WaterPoints2D(i, j) /= 1) then
                write(*,*)"A profile output was defined in the following point,"
                write(*,*)"meanwhile it is not a water point"
                write(*,*)"I = ", i
                write(*,*)"J = ", j
                write(*,*)"Please review the data files"
                AllTestsPassed = .false.
            endif

            Profile => Profile%Next

        enddo

        if (.not. AllTestsPassed) stop 'VerifyProfileLocation - ModuleProfile - ERR02'
    
    end subroutine VerifyProfileLocation


    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine GetProfileNextOutputTime(ProfileID, NextOutputTime, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ProfileID
        type (T_Time), optional                     :: NextOutputTime
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_   
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ProfileID, ready_)    

cd1 :   if ((ready_ .EQ. IDLE_ERR_) .OR. &
           (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(NextOutputTime)) then

                NextOutputTime = Me%OutTime(Me%OutputNumber)

            endif

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetProfileNextOutputTime


    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine WriteProfile_R4(ObjProfileID, Data3D, SZZ, LocI, LocJ, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjProfileID
        real(4), dimension(:,:,:), pointer          :: Data3D
        real   , dimension(:,:,:), pointer          :: SZZ
        integer, intent(IN ), optional              :: LocI, LocJ
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, STAT_CALL
        type(T_OneProfile),        pointer          :: Profile

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjProfileID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_) stop 'WriteProfile_R4 - ModuleProfile - ERR01'
            
            if (Me%ExternalVar%Now .ge. Me%OutTime(Me%OutputNumber)) then

                Me%PropertyCount = Me%PropertyCount + 1

                Profile => Me%FirstProfile

                do while(associated(Profile))

                    if(present(LocI))then
                        Profile%LocalizationI = LocI
                    end if

                    if(present(LocJ))then
                        Profile%LocalizationJ = LocJ
                    end if

                    Profile%Values(1:Me%KUB, Me%OutputNumber, Me%PropertyCount) = &
                    Data3D(Profile%LocalizationI, Profile%LocalizationJ, 1:Me%KUB)


                    Profile => Profile%Next

                end do

                if(Me%PropertyCount == Me%NumberOfProperties)then

                    Profile => Me%FirstProfile

                    do while(associated(Profile))
                
                        Profile%VerticalZ(:,Me%OutputNumber) = &
                                        SZZ(Profile%LocalizationI, Profile%LocalizationJ, 0:Me%KUB)

                        Profile => Profile%Next

                    end do

                    Me%PropertyCount = 0

                    Me%OutputNumber  = Me%OutputNumber + 1

                end if
                
            end if


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine WriteProfile_R4

    !------------------------------------------------------------------------

    subroutine WriteProfile_R8(ObjProfileID, Data3D, SZZ, LocI, LocJ, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjProfileID
        real(8), dimension(:,:,:), pointer          :: Data3D
        real   , dimension(:,:,:), pointer          :: SZZ
        integer, intent(IN ), optional              :: LocI, LocJ
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, STAT_CALL
        type(T_OneProfile),        pointer          :: Profile

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjProfileID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_) stop 'WriteProfile_R4 - ModuleProfile - ERR01'
            
            if (Me%ExternalVar%Now .ge. Me%OutTime(Me%OutputNumber)) then

                Me%PropertyCount = Me%PropertyCount + 1

                Profile => Me%FirstProfile

                do while(associated(Profile))

                    if(present(LocI))then
                        Profile%LocalizationI = LocI
                    end if

                    if(present(LocJ))then
                        Profile%LocalizationJ = LocJ
                    end if

                    Profile%Values(1:Me%KUB, Me%OutputNumber, Me%PropertyCount) = &
                    Data3D(Profile%LocalizationI, Profile%LocalizationJ, 1:Me%KUB)


                    Profile => Profile%Next

                end do

                if(Me%PropertyCount == Me%NumberOfProperties)then


                    Profile => Me%FirstProfile

                    do while(associated(Profile))
                
                        Profile%VerticalZ(:,Me%OutputNumber) = &
                                    SZZ(Profile%LocalizationI, Profile%LocalizationJ, 0:Me%KUB)

                        Profile => Profile%Next

                    end do

                    Me%PropertyCount = 0

                    Me%OutputNumber  = Me%OutputNumber + 1

                end if
                
            end if


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine WriteProfile_R8

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillProfile(ObjProfileID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjProfileID              
        integer, optional, intent(OUT)      :: STAT

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, ready_, nUsers, STAT_CALL
        type(T_OneProfile), pointer         :: Profile

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjProfileID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mProfile_,  Me%InstanceID)
            
            if (nUsers == 0) then

                call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)   
                if (STAT_CALL /= SUCCESS_) stop 'WriteProfile_R4 - ModuleProfile - ERR01'
                
                call WriteFinalFile

                
                deallocate(Me%Depth         )
                deallocate(Me%Time2D        )
                deallocate(Me%PropertyList  )

                Profile => Me%FirstProfile

                do while(associated(Profile))

                    deallocate(Profile%Values   )
                    deallocate(Profile%VerticalZ)

                    Profile => Profile%Next   
                end do


                nUsers = DeassociateInstance(mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillProfile - ModuleProfile - ERR00'

                !Deallocates Instance
                call DeallocateInstance ()


                ObjProfileID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillProfile
    
    !------------------------------------------------------------------------
    
    
    subroutine WriteFinalFile 

        !Local-------------------------------------------------------------------
        integer                                             :: OutputNumber, STAT_CALL
        type(T_OneProfile), pointer                         :: Profile
        character(len=StringLength)                         :: Name, Units
        integer                                             :: HDF5_CREATE, n
        real(8), dimension(:,:), pointer                    :: Data2D
        real,    dimension(:  ), pointer                    :: TimePtr
        real,    dimension(6  ), target                     :: AuxTime
        integer                                             :: iJulDay
        real(8)                                             :: JDay
        real                                                :: Hours, Minutes, Seconds

        !------------------------------------------------------------------------

        Profile => Me%FirstProfile

        do while(associated(Profile))

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
            !Opens HDF5 File
            call ConstructHDF5(Profile%ObjHDF5, trim(Profile%ID%Name), HDF5_CREATE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleProfile - ERR00'


            do n = 1, Me%NumberOfProperties

                call HDF5SetLimits  (Profile%ObjHDF5, 1, Me%KUB, 1, Me%nInstants, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleProfile - ERR10'


                Data2D      => Profile%Values(1:Me%KUB, 1:Me%nInstants, n)
                Name        =  trim(Me%PropertyList(n,1))
                Units       =  trim(Me%PropertyList(n,2))

                call HDF5WriteData(Profile%ObjHDF5,                                     &
                                   "/Results/"//trim(Name),                             &
                                   trim(Name),                                          &
                                   trim(Units),                                         &
                                   Array2D      = Data2D,                               &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleProfile - ERR50'

            enddo
            
            do n = 1, Me%KUB

                Me%Depth(n,:) =  0.5*(Profile%VerticalZ(n-1,:) + Profile%VerticalZ(n,:))

            enddo  
            
            call HDF5WriteData(Profile%ObjHDF5,                                         &
                               "Grid/Depth",                                            &
                               "Depth",                                                 &
                               "m",                                                     &
                               Array2D      = Me%Depth,                                 &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleProfile - ERR50'
          

            call HDF5SetLimits  (Profile%ObjHDF5, 0, Me%KUB, 1, Me%nInstants, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleProfile - ERR10'

            call HDF5WriteData(Profile%ObjHDF5,                                         &
                               "Grid/VerticalZ",                                        &
                               "VerticalZ",                                             &
                               "m",                                                     &
                               Array2D      = Profile%VerticalZ,                        &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleProfile - ERR50'
            
            do OutputNumber = 1, Me%nInstants
                        
                call HDF5SetLimits  (Profile%ObjHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleProfile - ERR10'

                
                call ExtractDate   (Me%OutTime(OutputNumber), AuxTime(1), AuxTime(2),   &
                                    AuxTime(3), AuxTime(4), AuxTime(5), AuxTime(6))
                
                TimePtr => AuxTime

                call HDF5WriteData  (Profile%ObjHDF5,                                   &
                                     "/Time",                                           &
                                     "Time",                                            &
                                     "YYYY/MM/DD HH:MM:SS",                             &
                                     Array1D        = TimePtr,                          &
                                     OutputNumber   = OutputNumber,                     &
                                     STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleProfile - ERR20'

                call JulianDay(Me%OutTime(OutputNumber), iJulDay)
                
                JDay = dble(iJulDay) - 1.
                
                Hours   = AuxTime(4)
                Minutes = AuxTime(5)
                Seconds = AuxTime(6)

                JDay = JDay + Hours/24. + Minutes/(1440.) + Seconds/86400.

                Me%Time2D(:, OutputNumber) = JDay

            end do
            
            call HDF5SetLimits  (Profile%ObjHDF5, 1, Me%KUB, 1, Me%nInstants, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleProfile - ERR10'

            call HDF5WriteData  (Profile%ObjHDF5,                                   &
                                 "/JulianDay",                                      &
                                 "JulianDay",                                       &
                                 "days",                                            &
                                 Array2D        = Me%Time2D,                        &
                                 STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleProfile - ERR20'


            call KillHDF5(Profile%ObjHDF5, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'WriteFinalFile - ModuleProfile - ERR01'

            Profile => Profile%Next

        enddo

    end subroutine WriteFinalFile 

    !------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Profile), pointer          :: AuxObjProfile
        type (T_Profile), pointer          :: PreviousObjProfile

        !Updates pointers
        if (Me%InstanceID == FirstObjProfile%InstanceID) then
            FirstObjProfile => FirstObjProfile%Next
        else
            PreviousObjProfile => FirstObjProfile
            AuxObjProfile      => FirstObjProfile%Next
            do while (AuxObjProfile%InstanceID /= Me%InstanceID)
                PreviousObjProfile => AuxObjProfile
                AuxObjProfile      => AuxObjProfile%Next
            enddo

            !Now update linked list
            PreviousObjProfile%Next => AuxObjProfile%Next

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

    subroutine Ready (ObjProfile_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjProfile_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjProfile_ID > 0) then
            call LocateObjProfile (ObjProfile_ID)
            ready_ = VerifyReadLock (mProfile_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjProfile (ObjProfileID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjProfileID

        !Local-----------------------------------------------------------------

        Me => FirstObjProfile
        do while (associated (Me))
            if (Me%InstanceID == ObjProfileID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleProfile - LocateObjProfile - ERR01'

    end subroutine LocateObjProfile

    !--------------------------------------------------------------------------

end module ModuleProfile

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------









