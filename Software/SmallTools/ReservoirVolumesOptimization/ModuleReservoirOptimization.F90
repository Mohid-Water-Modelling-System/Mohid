!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : CEQUALW2 preporocessor
! MODULE        : ModuleReservoirOptimization
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Ricardo Miranda - v1.0
! DESCRIPTION   : 
!
!------------------------------------------------------------------------------


module ModuleReservoirOptimization

    use ModuleGlobalData
    use ModuleEnterData

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructReservoirOptimization
    private ::      AllocateInstance
    private ::      Read_W2_Files_Name
    private ::      ReadInputData
    private ::          ReadInputDataFixedFormat
    private ::          ReadInputDataFreeFormat
    private ::              AllocateArrays
    private ::      ReadW2Bathymetry
    private ::          ReadW2BathymetryFixedFormat
    private ::          ReadW2BathymetryFreeFormat
    private ::              CalcSegments
    private ::      PopulateVariables
    private ::      ReadFillShape
    private ::          FillShapeMatrix
    private ::              CalcVol
    private ::                  InterpolateVol
    
    !Selector
    public  :: GetReservoirOptimizationPointer
    public  :: GetReservoirOptimizationInteger
    public  :: UnGetReservoirOptimization
    public  :: WriteW2Bathymetry
                     
    
    !Modifier
    public  :: ModifyReservoirOptimization
    private ::      CalcCurrentVolume
    private ::      CalcDiffVolume
    private ::      RecalculateVolume
    private ::      Verify              !Confirms layer width increases to the surface

    !Destructor
    public  :: KillReservoirOptimization                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjReservoirOptimization 
    
    !Interfaces----------------------------------------------------------------

    private :: UnGetReservoirOptimization3D_I
    private :: UnGetReservoirOptimization3D_R8
    interface  UnGetReservoirOptimization
        module procedure UnGetReservoirOptimization3D_I
        module procedure UnGetReservoirOptimization3D_R8
    end interface  UnGetReservoirOptimization
    
    !Constants-----------------------------------------------------------------

    character(Len = StringLength),  parameter   :: freeFormat   ='free'
    character(Len = StringLength),  parameter   :: fixedFormat  ='fixed'

    !Types---------------------------------------------------------------------
    
    private ::  T_CEQUALW2
    type        T_CEQUALW2
        character(Len = StringLength)                :: fileFormat            !Either "free"  or "fixed"
        character(Len = StringLength)                :: dataFormat            !Either "free"  or "fixed"
        CHARACTER(256), DIMENSION(:),    pointer     :: CommentHeader
        CHARACTER(256), DIMENSION(:),    pointer     :: CommentDLX
        CHARACTER(256), DIMENSION(:),    pointer     :: CommentELWS
        CHARACTER(256), DIMENSION(:),    pointer     :: CommentPHI0
        CHARACTER(256), DIMENSION(:),    pointer     :: CommentFRIC
        CHARACTER(256), DIMENSION(:),    pointer     :: CommentH
        CHARACTER(256), DIMENSION(:,:),  pointer     :: CommentB
        INTEGER                                      :: NWB
        INTEGER,        DIMENSION(:),    pointer     :: BTH
        INTEGER                                      :: CON
        INTEGER                                      :: IMX,    KMX
        INTEGER                                      :: NBR
        INTEGER,        DIMENSION(:),    pointer     :: KTWB,   KB,     US,     DS,     BS,     BE
        INTEGER,        DIMENSION(:),    pointer     :: KTI
        INTEGER,        DIMENSION(:),    pointer     :: UHS,    DHS,    UQB,    DQB,    JBDN,   NPOINT
        INTEGER,        DIMENSION(:,:),  pointer     :: NCCGR,  NCCBR
        CHARACTER(72),  DIMENSION(:),    pointer     :: BTHFN
        REAL(4),        DIMENSION(:),    pointer     :: Z,      DLX,    ELWS,   SHADE,  PHI0,   FRIC
        REAL(4),        DIMENSION(:,:),  pointer     :: H,      EL
        REAL(4),        DIMENSION(:,:),  pointer     :: B
        REAL(4),        DIMENSION(:),    pointer     :: ELBOT,  SLOPE
        REAL(8),        DIMENSION(:),    pointer     :: VOLB,   VOLG
        REAL(4),        DIMENSION(:,:),  pointer     :: SAGR,   SABR,   CVBR,   CVGR
        REAL(4),        DIMENSION(:),    pointer     :: ALPHA
        LOGICAL,        DIMENSION(:),    pointer     :: ZERO_SLOPE     
        LOGICAL,        DIMENSION(:),    pointer     :: UH_EXTERNAL,    DH_EXTERNAL,    UQ_EXTERNAL,    DQ_EXTERNAL,    DAM_FLOW
        LOGICAL,        DIMENSION(:),    pointer     :: UH_INTERNAL,    DH_INTERNAL,    UQ_INTERNAL,    DQ_INTERNAL
    end type    T_CEQUALW2
    
    private ::  T_FillData
    type        T_FillData
        real(4),        dimension(:),    pointer     :: Vol             !Dam Volume correspondig to T_CEQUALW2%H
        real(4),        dimension(:),    pointer     :: VolDiff
        real(4),        dimension(:),    pointer     :: VolError
        real(4)                                      :: DamVolError     !Allowed volume Error, between 0.0 and 1.0
    end type    T_FillData
    
    private ::  T_Files
    type        T_Files
        character(Len = StringLength)           :: FillShape            !File containig Vol/heigth relationship
        character(Len = StringLength)           :: W2Bathymetry         !CEQUAL-W2 bathymetry file
        character(Len = StringLength)           :: InputData            !IN_MODEL
        character(Len = StringLength)           :: ModifiedW2Bathymetry !Modified CEQUAL-W2 bathymetry file
        character(Len = StringLength)           :: W2ControlFile        !CEQUAL-W2 control file
    end type    T_Files

    private ::  T_ReservoirOptimization
    type        T_ReservoirOptimization
        integer                                     :: InstanceID
        type (T_Size3D)                             :: Size, WorkSize
        real(8), dimension(:, :, :),  pointer       :: Matrix
        type (T_Files   )                           :: Files
        type (T_CEQUALW2)                           :: CEQUALW2
        type (T_FillData)                           :: FillData
        !Instance of Module_EnterData           
        integer                                     :: ObjEnterData01         = 0   !InputData
        integer                                     :: ObjEnterData02         = 0   !Reservoir Free Format  File
        integer                                     :: ObjEnterData03         = 0   !Fill Reservoir Data
        type(T_ReservoirOptimization), pointer      :: Next
    end type    T_ReservoirOptimization

    !Global Module Variables
    type (T_ReservoirOptimization), pointer                         :: FirstObjReservoirOptimization
    type (T_ReservoirOptimization), pointer                         :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructReservoirOptimization(ObjReservoirOptimizationID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjReservoirOptimizationID 
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mRESERVOIROPTIMIZATION_)) then
            nullify (FirstObjReservoirOptimization)
            call RegisterModule (mRESERVOIROPTIMIZATION_) 
        endif

        call Ready(ObjReservoirOptimizationID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance   ()
            call Read_W2_Files_Name ()
            call ReadInputData      ()
            call ReadW2Bathymetry   ()
            call PopulateVariables  ()
            call ReadFillShape      ()

            !Returns ID
            ObjReservoirOptimizationID          = Me%InstanceID

            STAT_ = SUCCESS_

        else    cd0
            
            stop 'ConstructReservoirOptimization - ModuleReservoirOptimization - ERR01' 

        endif   cd0

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructReservoirOptimization
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance()

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_ReservoirOptimization), pointer                         :: NewObjReservoirOptimization
        type (T_ReservoirOptimization), pointer                         :: PreviousObjReservoirOptimization


        !Allocates new instance
        allocate (NewObjReservoirOptimization)
        nullify  (NewObjReservoirOptimization%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjReservoirOptimization)) then
            FirstObjReservoirOptimization           => NewObjReservoirOptimization
            Me                                      => NewObjReservoirOptimization
        else
            PreviousObjReservoirOptimization        => FirstObjReservoirOptimization
            Me                                      => FirstObjReservoirOptimization%Next
            do while (associated(Me))
                PreviousObjReservoirOptimization    => Me
                Me                                  => Me%Next
            enddo
            Me                                      => NewObjReservoirOptimization
            PreviousObjReservoirOptimization%Next   => NewObjReservoirOptimization
        endif

        Me%InstanceID = RegisterNewInstance (mRESERVOIROPTIMIZATION_)


    end subroutine AllocateInstance
    
    !--------------------------------------------------------------------------
    
    subroutine Read_W2_Files_Name()

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------

        integer                       :: STAT_CALL 
        character(len = StringLength) :: Message

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        ! ---> ASCII file used to construct new properties
        Message   ='ASCII CEQUAL-W2 bathymetry file.'
        Message   = trim(Message)

        call ReadFileName('W2_BATHYMETRY', Me%Files%W2Bathymetry,                       &
                           Message = Message, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Read_W2_Files_Name - ModuleReservoirOptimization - ERR01' 
            
            

        ! ---> ASCII file used to read reservoir volume/heoght ratio
        Message   ='ASCII Dam volume/height.'
        Message   = trim(Message)

        call ReadFileName('RESERVOIR_DATA', Me%Files%FillShape,                         &
                           Message = Message, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Read_W2_Files_Name - ModuleReservoirOptimization - ERR02' 



        ! ---> ASCII file with reservoir conversion data
        Message   ='Model data.'
        Message   = trim(Message)

        call ReadFileName('IN_MODEL', Me%Files%InputData,                               &
                           Message = Message, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Read_W2_Files_Name - ModuleReservoirOptimization - ERR03' 



        ! ---> ASCII file with modified bathymetry data
        Message   ='Optimized bathymetry file in CEQUAL-W2 format.'
        Message   = trim(Message)

        call ReadFileName('W2_BATHYMETRY_MODIFIED', Me%Files%ModifiedW2Bathymetry,      &
                           Message = Message, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Read_W2_Files_Name - ModuleReservoirOptimization - ERR04' 



        ! ---> ASCII file with modified bathymetry data
        Message   ='CEQUAL-W2 control file.'
        Message   = trim(Message)

        call ReadFileName('W2_CONTROL_FILE', Me%Files%W2ControlFile,                    &
                           Message = Message, STAT = STAT_CALL)
        if ((STAT_CALL .NE. SUCCESS_) .AND. (STAT_CALL .NE. KEYWORD_NOT_FOUND_ERR_))    &
            stop 'Read_W2_Files_Name - ModuleReservoirOptimization - ERR10' 

        !----------------------------------------------------------------------
        
    end subroutine Read_W2_Files_Name
    
    !--------------------------------------------------------------------------
    
    subroutine ReadInputData()

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------

        integer :: STAT_CALL 
        integer :: flag = 0
        integer :: FromFile

        !Local-----------------------------------------------------------------

        integer :: ObjEnterData01 = 0

        !----------------------------------------------------------------------

        call GetExtractType     (FromFile = FromFile)

        call ConstructEnterData(ObjEnterData01, Me%Files%InputData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInputData - ModuleReservoirOptimization - ERR10'




        call GetData           (Me%CEQUALW2%fileFormat,                 &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "W2_BATH_FORMAT")

cd11:   if      (flag .EQ. 0) then
            write(*,*) "Keyword W2_BATH_FORMAT not found in "//Me%Files%InputData
            stop       "ReadInputData - ModuleReservoirOptimization - ERR11"
        endif cd11





        call GetData           (Me%CEQUALW2%dataFormat,                 &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "W2_DATA_FORMAT")

cd12:   if      (flag .EQ. 0) then
            write(*,*) "Keyword W2_DATA_FORMAT not found in "//Me%Files%InputData
            stop       "ReadInputData - ModuleReservoirOptimization - ERR12"
        endif cd12




        call GetData           (Me%FillData%DamVolError,                &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "RESERVOIR_ERROR")

cd30 :  if (flag .EQ. 0) then
            write(*,*) "Keyword RESERVOIR_ERROR not found in "//Me%Files%InputData
            stop       "ReadInputData - ModuleReservoirOptimization - ERR40"
        endif cd30




cf13 :  if      (Me%CEQUALW2%dataFormat .EQ. freeFormat ) then
            call ReadInputDataFreeFormat (ObjEnterData01,               &
                                          FromFile)

        elseif  (Me%CEQUALW2%dataFormat .EQ. fixedFormat) then  cf13
            call ReadInputDataFixedFormat(              )
        endif                                                   cf13



        call KillEnterData     (ObjEnterData01, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInputData - ModuleReservoirOptimization - ERR90'

        !----------------------------------------------------------------------
        
    end subroutine ReadInputData
    
    !--------------------------------------------------------------------------
    
    subroutine ReadInputDataFreeFormat(ObjEnterData01,                  &
                                       FromFile)

        !Arguments-------------------------------------------------------------

        integer,        intent(IN ) :: ObjEnterData01
        integer,        intent(IN ) :: FromFile

        !External--------------------------------------------------------------

        integer :: STAT_CALL 
        integer :: flag = 0

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        call GetData           (Me%CEQUALW2%NWB,                        &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "NWB")

cd5 :   if (flag .EQ. 0) then
            write(*,*) "Keyword NWB not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR14"
        endif cd5




        call GetData           (Me%CEQUALW2%NBR,                        &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "NBR")

cd4 :   if (flag .EQ. 0) then
            write(*,*) "Keyword NBR not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR15"
        endif cd4




        call GetData           (Me%CEQUALW2%IMX,                        &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "IMX")

cd6 :   if (flag .EQ. 0) then
            write(*,*) "Keyword IMX not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR20"
        endif cd6




        call GetData           (Me%CEQUALW2%KMX,                        &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "KMX")

cd7 :   if (flag .EQ. 0) then
            write(*,*) "Keyword KMX not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR25"
        endif cd7





        call AllocateArrays()




        call GetData           (Me%CEQUALW2%BE,                         &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "BE")

cd8 :   if (flag .EQ. 0) then
            write(*,*) "Keyword BE not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR26"
        endif cd8





        call GetData           (Me%CEQUALW2%SLOPE,                      &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "SLOPE")

cd34 :  if (flag .EQ. 0) then
            write(*,*) "Keyword SLOPE not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR34"
        endif cd34




        call GetData           (Me%CEQUALW2%ELBOT,                      &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "ELBOT")

cd35 :  if (flag .EQ. 0) then
            write(*,*) "Keyword ELBOT not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR35"
        endif cd35




        call GetData           (Me%CEQUALW2%JBDN,                       &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "JBDN")

cd350 :  if (flag .EQ. 0) then
            write(*,*) "Keyword JBDN not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR350"
        endif cd350




        call GetData           (Me%CEQUALW2%DHS,                        &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "DHS")

cd36 :  if (flag .EQ. 0) then
            write(*,*) "Keyword DHS not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR35"
        endif cd36




        call GetData           (Me%CEQUALW2%UHS,                        &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "UHS")

cd37 :  if (flag .EQ. 0) then
            write(*,*) "Keyword UHS not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR37"
        endif cd37




        call GetData           (Me%CEQUALW2%US,                         &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "US")

cd40 :  if (flag .EQ. 0) then
            write(*,*) "Keyword US not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR40"
        endif cd40




        call GetData           (Me%CEQUALW2%DS,                         &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "DS")

cd50 :  if (flag .EQ. 0) then
            write(*,*) "Keyword DS not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR50"
        endif cd50




        call GetData           (Me%CEQUALW2%BS,                         &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "BS")

cd60 :  if (flag .EQ. 0) then
            write(*,*) "Keyword BS not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR60"
        endif cd60




        call GetData           (Me%CEQUALW2%UQB,                        &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "UQB")

cd62 :  if (flag .EQ. 0) then
            write(*,*) "Keyword UQB not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR62"
        endif cd62




        call GetData           (Me%CEQUALW2%DQB,                        &
                                ObjEnterData01,                         &
                                flag,                                   &
                                SearchType = FromFile,                  &
                                keyword    = "DQB")

cd64 :  if (flag .EQ. 0) then
            write(*,*) "Keyword DQB not found in "//Me%Files%InputData
            stop       "ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR64"
        endif cd64

        !----------------------------------------------------------------------
        
    end subroutine ReadInputDataFreeFormat
    
    !--------------------------------------------------------------------------
    
    subroutine ReadInputDataFixedFormat()

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
    
        integer                         :: STAT_CALL 
        integer                         :: IERR
        
        !Local-----------------------------------------------------------------
    
        integer                         :: I,  J
        integer                         :: JW, JB 
        integer                         :: auxint

        real(4)                         :: auxreal
        
        CHARACTER(8)                                        :: AID
        CHARACTER(72),     ALLOCATABLE, DIMENSION(:)        :: TITLE
        INTEGER                                             :: NTR,      NST,    NIW,    NWD,    NGT,    NSP,    NPI,    NPU
        INTEGER                                             :: NGC,      NSS,    NAL,    NEP,    NBOD
        INTEGER                                             :: NOD
        INTEGER                                             :: NDLT
        REAL(4),            ALLOCATABLE, DIMENSION(:)       :: auxrealNOD
        INTEGER,            ALLOCATABLE, DIMENSION(:)       :: auxintNBR
        REAL(4),            ALLOCATABLE, DIMENSION(:)       :: auxrealNWB
        CHARACTER(8),       ALLOCATABLE, DIMENSION(:)       :: auxcharNWB

        !----------------------------------------------------------------------

        call UnitsManager(Me%CEQUALW2%CON, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR10'
                
        OPEN (Me%CEQUALW2%CON,FILE=Me%Files%W2ControlFile,STATUS='OLD',IOSTAT=IERR)
cf1 :   IF (IERR == 0) THEN

            ! Title cards

            ALLOCATE (TITLE(11))
            READ (Me%CEQUALW2%CON,'(//A8/(8X,A72))',ERR=400) AID, (TITLE(J),J=1,10)
            IF (AID /= 'TITLE C ')               GO TO 400

    ! Array dimensions

            READ (Me%CEQUALW2%CON,'(/A8/(8X,10I8))')         AID, Me%CEQUALW2%NWB, Me%CEQUALW2%NBR, Me%CEQUALW2%IMX, Me%CEQUALW2%KMX
            IF (AID /= 'GRID    ')               GO TO 401
            READ (Me%CEQUALW2%CON,'(/A8/(8X,10I8))')         AID,   NTR,    NST,    NIW,    NWD,    NGT,    NSP,    NPI,    NPU
            IF (AID /= 'IN/OUTFL')               GO TO 402
            READ (Me%CEQUALW2%CON,'(/A8/(8X,10I8))')         AID,   NGC,    NSS,    NAL,    NEP,    NBOD
            IF (AID /= 'CONSTITU')               GO TO 403
            READ (Me%CEQUALW2%CON,'(/A8/(8X,10I8))')         AID,   NOD
            IF (AID /= 'MISCELL ')               GO TO 404





            call AllocateArrays()
            ALLOCATE (auxrealNOD(NOD            ))
            ALLOCATE (auxcharNWB(Me%CEQUALW2%NWB))
            ALLOCATE (auxrealNWB(Me%CEQUALW2%NWB))
            ALLOCATE (auxintNBR (Me%CEQUALW2%NBR))



            ! Constituent numbers
  
!            NTDS  = 1
!            NGCS  = 2
!            NGCE  = NGCS+NGC-1
!            NSSS  = NGCE+1
!            NSSE  = NSSS+NSS-1
!            NPO4  = NSSE+1
!            NNH4  = NPO4+1
!            NNO3  = NNH4+1
!            NDSI  = NNO3+1
!            NPSI  = NDSI+1
!            NFE   = NPSI+1
!            NLDOM = NFE+1
!            NRDOM = NLDOM+1
!            NLPOM = NRDOM+1
!            NRPOM = NLPOM+1
!            NBODS = NRPOM+1
!            NBODE = NBODS+NBOD-1
!            NAS   = NBODE+1
!            NAE   = NAS+NAL-1
!            NDO   = NAE+1
!            NTIC  = NDO+1
!            NALK  = NTIC+1

            ! Constituent, tributary, and widthdrawal totals

!            NCT  = NALK
!            NTRT = NTR+NGT+NSP+NPI+NPU
!            NWDT = NWD+NGT+NSP+NPI+NPU

            ! Time control cards

!           READ (Me%CEQUALW2%CON,'(/A8/8X,2F8.0,I8)',       ERR=405)   AID,    TMSTRT,     TMEND,      YEAR
            READ (Me%CEQUALW2%CON,'(/A8/8X,2F8.0,I8)',       ERR=405)   AID,    auxreal,    auxreal,    auxint
            IF (AID /= 'TIME CON')                         GO TO 405
!           READ (Me%CEQUALW2%CON,'(/A8/8X,I8,F8.0)',        ERR=406)   AID,    NDLT,   DLTMIN
            READ (Me%CEQUALW2%CON,'(/A8/8X,I8,F8.0)',        ERR=406)   AID,    NDLT,   auxreal
            IF (AID /= 'DLT CON ')                         GO TO 406
!           READ (Me%CEQUALW2%CON,'(/A8/(:8X,9F8.0))',       ERR=407)   AID, (DLTD      (J),J=1,NDLT)
            READ (Me%CEQUALW2%CON,'(/A8/(:8X,9F8.0))',       ERR=407)   AID, (auxrealNOD(J),J=1,NDLT)
            IF (AID /= 'DLT DATE')                         GO TO 407
!           READ (Me%CEQUALW2%CON,'(/A8/(:8X,9F8.0))',       ERR=408)   AID, (DLTMAX    (J),J=1,NDLT)
            READ (Me%CEQUALW2%CON,'(/A8/(:8X,9F8.0))',       ERR=408)   AID, (auxrealNOD(J),J=1,NDLT)
            IF (AID /= 'DLT MAX ')                         GO TO 408
!           READ (Me%CEQUALW2%CON,'(/A8/(:8X,9F8.0))',       ERR=409)   AID, (DLTF      (J),J=1,NDLT)
            READ (Me%CEQUALW2%CON,'(/A8/(:8X,9F8.0))',       ERR=409)   AID, (auxrealNOD(J),J=1,NDLT)
            IF (AID /= 'DLT FRN ')                         GO TO 409
!           READ (Me%CEQUALW2%CON,'(/A8/(8X,2A8))',          ERR=410)  AID, (VISC       (JW),   CELC        (JW), JW=1,Me%CEQUALW2%NWB)
            READ (Me%CEQUALW2%CON,'(/A8/(8X,2A8))',          ERR=410)  AID, (auxcharNWB (JW),   auxcharNWB  (JW), JW=1,Me%CEQUALW2%NWB)
            IF (AID /= 'DLT LIMI')                         GO TO 410


            ! Grid definition cards

!           READ (Me%CEQUALW2%CON,'(/A8/(8X,7I8,F8.3))',    ERR=411)  AID, (Me%CEQUALW2%US   (JB),  &
!                                                                           Me%CEQUALW2%DS   (JB),  &
!                                                                           Me%CEQUALW2%UHS  (JB),  &
!                                                                           Me%CEQUALW2%DHS  (JB),  &
!                                                                           Me%CEQUALW2%UQB  (JB),  &
!                                                                           Me%CEQUALW2%DQB  (JB),  &
!                                                                           NL               (JB),  &
!                                                                           Me%CEQUALW2%SLOPE(JB),  JB=1,Me%CEQUALW2%NBR)
            READ (Me%CEQUALW2%CON,'(/A8/(8X,7I8,F8.3))',    ERR=411)  AID, (Me%CEQUALW2%US   (JB),  &
                                                                            Me%CEQUALW2%DS   (JB),  &
                                                                            Me%CEQUALW2%UHS  (JB),  &
                                                                            Me%CEQUALW2%DHS  (JB),  &
                                                                            Me%CEQUALW2%UQB  (JB),  &
                                                                            Me%CEQUALW2%DQB  (JB),  &
                                                                            auxintNBR        (JB),  &
                                                                            Me%CEQUALW2%SLOPE(JB),  JB=1,Me%CEQUALW2%NBR)
            IF (AID /= 'BRANCH G')                        GO TO 411
!           READ (Me%CEQUALW2%CON,'(/A8/(8X,3F8.0,3I8))',   ERR=412)  AID, (LAT              (JW),  &
!                                                                           LONG             (JW),  &
!                                                                           ELBOT            (JW),  &
!                                                                           BS               (JW),  &
!                                                                           BE               (JW),  &
!                                                                           JBDN             (JW),  JW=1,Me%CEQUALW2%NWB)
            READ (Me%CEQUALW2%CON,'(/A8/(8X,3F8.0,3I8))',   ERR=412)  AID, (auxrealNWB       (JW),  &
                                                                            auxrealNWB       (JW),  &
                                                                            Me%CEQUALW2%ELBOT(JW),  &
                                                                            Me%CEQUALW2%BS   (JW),  &
                                                                            Me%CEQUALW2%BE   (JW),  &
                                                                            Me%CEQUALW2%JBDN (JW),  JW=1,Me%CEQUALW2%NWB)
            IF (AID /= 'LOCATION')                        GO TO 412

                
            call UnitsManager          (Me%CEQUALW2%CON, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadW2BathymetryFixedFormat - ModuleReservoirOptimization - ERR85' 
        ELSE                                                                cf1
            write(*,*)  'Could not open '//Me%Files%W2ControlFile
            stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR90'
        ENDIF                                                               cf1
        
        !----------------------------------------------------------------------
        
        return
        
400     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR400'
        
401     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR401'
        
402     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR402'
        
403     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR403'
        
404     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR404'
        
405     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR405'
        
406     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR406'
        
407     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR407'
        
408     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR408'
        
409     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR409'
        
410     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR410'
        
411     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR411'
        
412     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadInputDataFixedFormat - ModuleReservoirOptimization - ERR412'

        !----------------------------------------------------------------------
        
    end subroutine ReadInputDataFixedFormat

    !--------------------------------------------------------------------------
    
    subroutine AllocateArrays()
        
        !Arguments-------------------------------------------------------------
    
        !External--------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
    
        !----------------------------------------------------------------------
    
        ! Variable initializations
        nullify (Me%CEQUALW2%CommentHeader  )
        nullify (Me%CEQUALW2%CommentDLX     )
        nullify (Me%CEQUALW2%CommentELWS    )
        nullify (Me%CEQUALW2%CommentPHI0    )
        nullify (Me%CEQUALW2%CommentFRIC    )
        nullify (Me%CEQUALW2%CommentH       )
        nullify (Me%CEQUALW2%CommentB       )
        nullify (Me%CEQUALW2%BTH            )
        nullify (Me%CEQUALW2%BTHFN          )
        nullify (Me%CEQUALW2%BS             )
        nullify (Me%CEQUALW2%BE             )
        nullify (Me%CEQUALW2%US             )
        nullify (Me%CEQUALW2%DS             )
        nullify (Me%CEQUALW2%DLX            )
        nullify (Me%CEQUALW2%ELWS           )
        nullify (Me%CEQUALW2%PHI0           )
        nullify (Me%CEQUALW2%FRIC           )
        nullify (Me%CEQUALW2%H              )
        nullify (Me%CEQUALW2%EL             )
        nullify (Me%CEQUALW2%B              )
        nullify (Me%CEQUALW2%KTI            )
        nullify (Me%CEQUALW2%KTWB           )
        nullify (Me%CEQUALW2%KB             )
        nullify (Me%CEQUALW2%ELBOT          )
        nullify (Me%CEQUALW2%SLOPE          )
        nullify (Me%CEQUALW2%ZERO_SLOPE     )
        nullify (Me%CEQUALW2%VOLB           )
        nullify (Me%CEQUALW2%VOLG           )
        nullify (Me%CEQUALW2%UHS            )
        nullify (Me%CEQUALW2%DHS            )
        nullify (Me%CEQUALW2%UQB            )
        nullify (Me%CEQUALW2%DQB            )
        nullify (Me%CEQUALW2%JBDN           )
        nullify (Me%CEQUALW2%NPOINT         )
        nullify (Me%CEQUALW2%DAM_FLOW       )
        nullify (Me%CEQUALW2%UH_EXTERNAL    )
        nullify (Me%CEQUALW2%DH_EXTERNAL    )
        nullify (Me%CEQUALW2%UQ_EXTERNAL    )
        nullify (Me%CEQUALW2%DQ_EXTERNAL    )
        nullify (Me%CEQUALW2%UH_INTERNAL    )
        nullify (Me%CEQUALW2%DH_INTERNAL    )
        nullify (Me%CEQUALW2%UQ_INTERNAL    )
        nullify (Me%CEQUALW2%DQ_INTERNAL    )
        nullify (Me%CEQUALW2%ALPHA          )
        nullify (Me%CEQUALW2%Z              )
        nullify (Me%CEQUALW2%NCCGR          )
        nullify (Me%CEQUALW2%NCCBR          )
        nullify (Me%CEQUALW2%SAGR           )
        nullify (Me%CEQUALW2%SABR           )
        nullify (Me%CEQUALW2%CVBR           )
        nullify (Me%CEQUALW2%CVGR           )
       
        allocate(Me%CEQUALW2%CommentHeader  (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%CommentDLX     (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%CommentELWS    (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%CommentPHI0    (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%CommentFRIC    (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%CommentH       (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%CommentB       (Me%CEQUALW2%NWB, Me%CEQUALW2%KMX   ))
        allocate(Me%CEQUALW2%BTH            (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%BTHFN          (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%BS             (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%BE             (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%US             (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%DS             (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%DLX            (Me%CEQUALW2%IMX                    ))
        allocate(Me%CEQUALW2%ELWS           (Me%CEQUALW2%IMX                    ))
        allocate(Me%CEQUALW2%PHI0           (Me%CEQUALW2%IMX                    ))
        allocate(Me%CEQUALW2%FRIC           (Me%CEQUALW2%IMX                    ))
        allocate(Me%CEQUALW2%H              (Me%CEQUALW2%KMX, Me%CEQUALW2%NWB   ))
        allocate(Me%CEQUALW2%EL             (Me%CEQUALW2%KMX, Me%CEQUALW2%IMX   ))
        allocate(Me%CEQUALW2%B              (Me%CEQUALW2%KMX, Me%CEQUALW2%IMX   ))
        allocate(Me%CEQUALW2%KTI            (Me%CEQUALW2%IMX                    ))
        allocate(Me%CEQUALW2%KTWB           (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%KB             (Me%CEQUALW2%IMX                    ))
        allocate(Me%CEQUALW2%ELBOT          (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%SLOPE          (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%ZERO_SLOPE     (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%VOLB           (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%VOLG           (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%UHS            (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%DHS            (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%UQB            (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%DQB            (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%JBDN           (Me%CEQUALW2%NWB                    ))
        allocate(Me%CEQUALW2%NPOINT         (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%DAM_FLOW       (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%UH_EXTERNAL    (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%DH_EXTERNAL    (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%UH_INTERNAL    (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%DH_INTERNAL    (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%UQ_EXTERNAL    (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%DQ_EXTERNAL    (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%UQ_INTERNAL    (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%DQ_INTERNAL    (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%ALPHA          (Me%CEQUALW2%NBR                    ))
        allocate(Me%CEQUALW2%Z              (Me%CEQUALW2%IMX                    ))
        allocate(Me%CEQUALW2%NCCGR          (Me%CEQUALW2%KMX, Me%CEQUALW2%NWB   ))
        allocate(Me%CEQUALW2%NCCBR          (Me%CEQUALW2%KMX, Me%CEQUALW2%NBR   ))
        allocate(Me%CEQUALW2%SAGR           (Me%CEQUALW2%KMX, Me%CEQUALW2%NWB   ))
        allocate(Me%CEQUALW2%SABR           (Me%CEQUALW2%KMX, Me%CEQUALW2%NBR   ))
        allocate(Me%CEQUALW2%CVBR           (Me%CEQUALW2%KMX, Me%CEQUALW2%NBR   ))
        allocate(Me%CEQUALW2%CVGR           (Me%CEQUALW2%KMX, Me%CEQUALW2%NWB   ))
        
        Me%CEQUALW2%B (:,:)   = 0.0
        
        !----------------------------------------------------------------------
        
    end subroutine AllocateArrays


    !--------------------------------------------------------------------------
    
    subroutine ReadW2Bathymetry()
        
        !Arguments-------------------------------------------------------------
    
        !External--------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
    
        !----------------------------------------------------------------------
    
        Me%CEQUALW2%BTHFN(Me%CEQUALW2%NWB)=Me%Files%W2Bathymetry
        
        ! Bathymetry definition
        
        !'Bathymetry file'
cf11 :  if      (Me%CEQUALW2%fileFormat .EQ. freeFormat ) then
            call ReadW2BathymetryFreeFormat ()

        elseif  (Me%CEQUALW2%fileFormat .EQ. fixedFormat) then  cf11
            call ReadW2BathymetryFixedFormat()
        endif                                                   cf11

        !----------------------------------------------------------------------
        
    end subroutine ReadW2Bathymetry

    !--------------------------------------------------------------------------
    
    subroutine ReadW2BathymetryFixedFormat()
        
        !Arguments-------------------------------------------------------------
    
        !External--------------------------------------------------------------
    
        integer                         :: STAT_CALL 
        integer                         :: IERR
        
        !Local-----------------------------------------------------------------
    
        integer                         :: I, K
        integer                         :: JW 
        
        !----------------------------------------------------------------------

cl1 :   DO JW=1,Me%CEQUALW2%NWB
            call UnitsManager(Me%CEQUALW2%BTH(JW), OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadW2BathymetryFixedFormat - ModuleReservoirOptimization - ERR10'
                
            OPEN (Me%CEQUALW2%BTH(JW),FILE=Me%CEQUALW2%BTHFN(JW),STATUS='OLD',IOSTAT=IERR)
cf1 :       IF (IERR == 0) THEN
                READ  (Me%CEQUALW2%BTH(JW),'(A256)')                                                &
                    Me%CEQUALW2%CommentHeader(JW)
                !'  segment lengths'
                READ  (Me%CEQUALW2%BTH(JW),'(/A256/(10F8.0))',ERR=400)                              &
                    Me%CEQUALW2%CommentDLX(JW),(Me%CEQUALW2%DLX(I),   I= 1, Me%CEQUALW2%IMX)
                !'  water surface elevations'
                READ  (Me%CEQUALW2%BTH(JW),'(/A256/(10F8.0))',ERR=410)                              &
                    Me%CEQUALW2%CommentELWS(JW),(Me%CEQUALW2%ELWS(I), I= 1, Me%CEQUALW2%IMX)
                !'  segment orientation'
                READ  (Me%CEQUALW2%BTH(JW),'(/A256/(10F8.0))',ERR=420)                              &
                    Me%CEQUALW2%CommentPHI0(JW),(Me%CEQUALW2%PHI0(I), I= 1, Me%CEQUALW2%IMX)
                !'  segment bottom friction'
                READ  (Me%CEQUALW2%BTH(JW),'(/A256/(10F8.0))',ERR=430)                              &
                    Me%CEQUALW2%CommentFRIC(JW),(Me%CEQUALW2%FRIC(I), I= 1, Me%CEQUALW2%IMX)
                !'  layer thickness'
                READ  (Me%CEQUALW2%BTH(JW),'(/A256/(10F8.0))',ERR=440)                              &
                    Me%CEQUALW2%CommentH(JW),(Me%CEQUALW2%H(K,JW), K=1,Me%CEQUALW2%KMX)
                !'  segment widths'
cl2 :           DO I= 1, Me%CEQUALW2%IMX
                    READ (Me%CEQUALW2%BTH(JW),'(/A256/(10F8.0))',ERR=450)                           &
                        Me%CEQUALW2%CommentB(JW,I),(Me%CEQUALW2%B(K,I),  K=1,Me%CEQUALW2%KMX)
                ENDDO cl2
                
                call UnitsManager          (Me%CEQUALW2%BTH(JW), CLOSE_FILE, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadW2BathymetryFixedFormat - ModuleReservoirOptimization - ERR20' 
            ELSE    cf1
                write(*,*)  'Could not open '//Me%CEQUALW2%BTHFN(JW)
                stop        'ReadW2BathymetryFixedFormat - ModuleReservoirOptimization - ERR30'
            ENDIF  cf1
        ENDDO cl1
        
        !----------------------------------------------------------------------
        
        return
        
400     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadW2BathymetryFixedFormat - ModuleReservoirOptimization - ERR400'
                        
410     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadW2BathymetryFixedFormat - ModuleReservoirOptimization - ERR410'
                        
420     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadW2BathymetryFixedFormat - ModuleReservoirOptimization - ERR420'
                        
430     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadW2BathymetryFixedFormat - ModuleReservoirOptimization - ERR430'
                        
440     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadW2BathymetryFixedFormat - ModuleReservoirOptimization - ERR440'
                        
450     CONTINUE
        write(*,*)  'Either illegal value or incorrect card somewhere in the following cards'
        stop        'ReadW2BathymetryFixedFormat - ModuleReservoirOptimization - ERR450'
                        
        !----------------------------------------------------------------------
        
    end subroutine ReadW2BathymetryFixedFormat

    !--------------------------------------------------------------------------
    
    subroutine ReadW2BathymetryFreeFormat()
        
        !Arguments-------------------------------------------------------------
    
        !External--------------------------------------------------------------
    
        integer :: STAT_CALL 
        integer :: flag = 0
        integer :: FromFile, FromBlock
        integer :: ClientNumber
        integer :: FirstLine, LastLine
        logical :: BlockFound
        REAL(4),        DIMENSION(:),   pointer     :: H1D
        CHARACTER(256)                              :: Comment
        
        !Local-----------------------------------------------------------------
    
        integer :: JW, I
        integer :: ObjEnterData02 = 0
        
        !----------------------------------------------------------------------

cl1 :   DO JW=1,Me%CEQUALW2%NWB
            call GetExtractType     (FromFile   = FromFile  )
            call GetExtractType     (FromBlock  = FromBlock )

            call ConstructEnterData(ObjEnterData02, Me%Files%W2Bathymetry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR01'




            !------------------------------------------------------------------



            call GetData           (Me%CEQUALW2%CommentHeader(JW),          &
                                    ObjEnterData02,                         &
                                    flag,                                   &
                                    SearchType  = FromFile,                 &
                                    Default     ="CQUAL-W2 Bathymetry",     &
                                    keyword     ="COMMENT_HEADER",          &
                                    STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR08'




            !------------------------------------------------------------------




            call GetData           (Me%CEQUALW2%CommentDLX(JW),             &
                                    ObjEnterData02,                         &
                                    flag,                                   &
                                    SearchType  = FromFile,                 &
                                    Default     ="CQUAL-W2 DLX",            &
                                    keyword     ="COMMENT_DLX",             &
                                    STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR09'




            call GetData           (Me%CEQUALW2%DLX,                        &
                                    ObjEnterData02,                         &
                                    flag,                                   &
                                    SearchType = FromFile,                  &
                                    keyword    = "DLX")

cf10 :      if (flag .EQ. 0) then
                write(*,*) "Keyword DLX not found in "//Me%Files%W2Bathymetry
                stop       "ReadW2BathymetryFreeFormat - ModuleReservoirOptimization - ERR10"
            endif cf10




            !------------------------------------------------------------------




            call GetData           (Me%CEQUALW2%CommentELWS(JW),            &
                                    ObjEnterData02,                         &
                                    flag,                                   &
                                    SearchType  = FromFile,                 &
                                    Default     ="CQUAL-W2 ELWS",           &
                                    keyword     ="COMMENT_ELWS",            &
                                    STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR14'




            call GetData           (Me%CEQUALW2%ELWS,                       &
                                    ObjEnterData02,                         &
                                    flag,                                   &
                                    SearchType = FromFile,                  &
                                    keyword    = "ELWS")

cf15 :      if (flag .EQ. 0) then
                write(*,*) "Keyword ELWS not found in "//Me%Files%W2Bathymetry
                stop       "ReadW2BathymetryFreeFormat - ModuleReservoirOptimization - ERR15"
            endif cf15




            !------------------------------------------------------------------




            call GetData           (Me%CEQUALW2%CommentPHI0(JW),            &
                                    ObjEnterData02,                         &
                                    flag,                                   &
                                    SearchType  = FromFile,                 &
                                    Default     ="CQUAL-W2 PHI0",           &
                                    keyword     ="COMMENT_PHI0",            &
                                    STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR19'




            call GetData           (Me%CEQUALW2%PHI0,                       &
                                    ObjEnterData02,                         &
                                    flag,                                   &
                                    SearchType = FromFile,                  &
                                    keyword    = "PHI0")

cf20 :      if (flag .EQ. 0) then
                write(*,*) "Keyword PHI0 not found in "//Me%Files%W2Bathymetry
                stop       "ReadW2BathymetryFreeFormat - ModuleReservoirOptimization - ERR20"
            endif cf20




            !------------------------------------------------------------------




            call GetData           (Me%CEQUALW2%CommentFRIC(JW),            &
                                    ObjEnterData02,                         &
                                    flag,                                   &
                                    SearchType  = FromFile,                 &
                                    Default     ="CQUAL-W2 FRIC",           &
                                    keyword     ="COMMENT_FRIC",            &
                                    STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR24'




            call GetData           (Me%CEQUALW2%FRIC,                       &
                                    ObjEnterData02,                         &
                                    flag,                                   &
                                    SearchType = FromFile,                  &
                                    keyword    = "FRIC")

cf25 :      if (flag .EQ. 0) then
                write(*,*) "Keyword FRIC not found in "//Me%Files%W2Bathymetry
                stop       "ReadW2BathymetryFreeFormat - ModuleReservoirOptimization - ERR25"
            endif cf25




            !------------------------------------------------------------------




            call GetData           (Me%CEQUALW2%CommentH(JW),               &
                                    ObjEnterData02,                         &
                                    flag,                                   &
                                    SearchType  = FromFile,                 &
                                    Default     ="CQUAL-W2 H",              &
                                    keyword     ="COMMENT_H",               &
                                    STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR29'




            allocate(H1D(Me%CEQUALW2%KMX))

            call GetData           (H1D,                                    &
                                    ObjEnterData02,                         &
                                    flag,                                   &
                                    SearchType = FromFile,                  &
                                    keyword    = "H")

cf30 :      if (flag .EQ. 0) then
                write(*,*) "Keyword H not found in "//Me%Files%W2Bathymetry
                stop       "ReadW2BathymetryFreeFormat - ModuleReservoirOptimization - ERR30"
            endif cf30

            Me%CEQUALW2%H(:,JW)=H1D(:)

            deallocate(H1D)




            !------------------------------------------------------------------




            call GetData           (Comment,                                    &
                                    ObjEnterData02,                             &
                                    flag,                                       &
                                    SearchType  = FromFile,                     &
                                    Default     ="CQUAL-W2 B",                  &
                                    keyword     ="COMMENT_SEGMENT",             &
                                    STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ReadInputDataFreeFormat - ModuleReservoirOptimization - ERR34'




            call ExtractBlockFromBuffer(ObjEnterData02, ClientNumber,           &
                                                        '<Segments>',           &
                                                        '</Segments>',          &
                                                        BlockFound,             &
                                                        FirstLine = FirstLine,  &
                                                        LastLine  = LastLine,   &
                                                        STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ReadW2BathymetryFreeFormat - ModuleReservoirOptimization - ERR35'

            call CalcSegments(ObjEnterData02,                                   &
                                FirstLine,                                      &
                                LastLine,                                       &
                                JW,                                             &
                                Comment)




            !------------------------------------------------------------------




            call KillEnterData     (ObjEnterData02, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadW2BathymetryFreeFormat - ModuleReservoirOptimization - ERR90'
        ENDDO cl1
                        
        !----------------------------------------------------------------------
        
    end subroutine ReadW2BathymetryFreeFormat

    !--------------------------------------------------------------------------
    
    subroutine CalcSegments(ObjEnterData,                                   &
                            FirstLine,                                      &
                            LastLine,                                       &
                            JW,                                             &
                            Comment)
        
        !Arguments-------------------------------------------------------------
    
        integer,        intent(IN ) :: ObjEnterData
        integer,        intent(IN ) :: JW
        integer,        intent(IN ) :: FirstLine    
        integer,        intent(IN ) :: LastLine    
        CHARACTER(256), intent(IN ) :: Comment

        !External--------------------------------------------------------------
    
        integer                 :: STAT_CALL 
        integer                 :: flag
        real(4), dimension(4)   :: vector

        !Local-----------------------------------------------------------------

        integer                 :: colI     = 1
        integer                 :: colB_K2  = 2
        integer                 :: colB_KN  = 3
        integer                 :: colKN    = 4
        integer                 :: line
        integer                 :: I,   KN
        integer                 :: K
        real(4)                 :: B2,  BN

        !----------------------------------------------------------------------

        if (((LastLine-1) - (FirstLine+1) + 1) .NE. Me%CEQUALW2%IMX)        &
            stop 'Wrong nbr of segments - CalcSegments - ModuleReservoirOptimization - ERR01'


cl10 :  do line = (FirstLine+1), (LastLine-1)
            I       = null_int
            KN      = null_int
            B2      = null_real
            BN      = null_real



            call GetData(vector,                                            &
                            ObjEnterData,                                   &
                            flag        = flag,                             &
                            Buffer_Line = line,                             &
                            STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CalcSegments - ModuleReservoirOptimization - ERR10'
            I       = vector(colI   )
            B2      = vector(colB_K2)
            BN      = vector(colB_KN)
            KN      = vector(colKN  )

cf20 :      if ((B2 .NE. 0.0) .AND. (BN .NE. 0.0) .AND. (KN .NE. 0)) then
cl20 :      do K=2, KN
                Me%CEQUALW2%B(K,I) = B2 + (BN-B2)/(float(KN)-2.0)*(float(K)-2.0)
            enddo   cl20
            endif   cf20

            !Comments
            Me%CEQUALW2%CommentB(JW,I) = trim(Comment)//" "
        enddo       cl10

        !----------------------------------------------------------------------
        
    end subroutine CalcSegments

    !--------------------------------------------------------------------------

    subroutine PopulateVariables()
        
        !Arguments-------------------------------------------------------------
    
        !External--------------------------------------------------------------
    
        !Local-----------------------------------------------------------------

        integer :: KTMAX = null_int
        integer :: IFLAG
        integer :: JW, JB, I, K, JJB
        integer :: NUP, NNBR, NCBP, NINTERNAL, JCHECK, NNBROLD

        real    :: ZMIN

        !----------------------------------------------------------------------

        ! Initialize logical control variables
cl40 :  DO JB=1,Me%CEQUALW2%NBR
            Me%CEQUALW2%DAM_FLOW   (JB) = Me%CEQUALW2%UHS(JB) <  -1
            Me%CEQUALW2%UH_EXTERNAL(JB) = Me%CEQUALW2%UHS(JB) == -1
            Me%CEQUALW2%DH_EXTERNAL(JB) = Me%CEQUALW2%DHS(JB) == -1

            Me%CEQUALW2%UH_INTERNAL(JB) = Me%CEQUALW2%UHS(JB) >   0
            Me%CEQUALW2%DH_INTERNAL(JB) = Me%CEQUALW2%DHS(JB) >   0
            Me%CEQUALW2%UQ_EXTERNAL(JB) = Me%CEQUALW2%UHS(JB) ==  0

            Me%CEQUALW2%DQ_EXTERNAL(JB) = Me%CEQUALW2%DHS(JB) ==  0
            Me%CEQUALW2%DQ_INTERNAL(JB) = Me%CEQUALW2%DQB(JB) >   0
            Me%CEQUALW2%UQ_INTERNAL(JB) = Me%CEQUALW2%UQB(JB) >   0 .AND. .NOT. Me%CEQUALW2%DAM_FLOW(JB)

            Me%CEQUALW2%UHS(JB)         = ABS(Me%CEQUALW2%UHS(JB))
            Me%CEQUALW2%DHS(JB)         = ABS(Me%CEQUALW2%DHS(JB))
        ENDDO                               cl40


        ! Grid linkage
cl30 :  DO JB=1,Me%CEQUALW2%NBR
            JCHECK = 0
cf30 :      IF (Me%CEQUALW2%UHS(JB) /= 0 .AND. .NOT. Me%CEQUALW2%UH_EXTERNAL(JB) .AND. .NOT. Me%CEQUALW2%DH_EXTERNAL(JB)) THEN
cl31 :          DO JJB=1,Me%CEQUALW2%NBR
                    IF (ABS(Me%CEQUALW2%UHS(JB)) >= Me%CEQUALW2%US(JJB) .AND. ABS(Me%CEQUALW2%UHS(JB)) <= Me%CEQUALW2%DS(JJB)) EXIT
                ENDDO                       cl31
cf31 :          IF (ABS(Me%CEQUALW2%UHS(JB)) == Me%CEQUALW2%DS(JJB)) THEN
                    IF (Me%CEQUALW2%DHS(JJB) == Me%CEQUALW2%US(JB)) JCHECK = 1
cf32 :              IF (Me%CEQUALW2%UHS(JB) < 0) THEN
cf33 :              IF (JCHECK == 1) THEN
                        WRITE (*,*) 'Error in grid linkage - check [US], [DS], [UHS], [DHS], and/or [JBDN] for waterbody ',JW
                        STOP
                    ENDIF                   cf33
                    ENDIF                   cf32
                ENDIF                       cf31
            ENDIF                           cf30
        ENDDO                               cl30
        Me%CEQUALW2%ZERO_SLOPE(:) = .TRUE.
cl50 :  DO JW=1,Me%CEQUALW2%NWB
cl51 :  DO JB=Me%CEQUALW2%BS(JW),Me%CEQUALW2%BE(JW)
            IF (Me%CEQUALW2%SLOPE(JB) /= 0.0) Me%CEQUALW2%ZERO_SLOPE(JW) = .FALSE.
        ENDDO                               cl51
        ENDDO                               cl50


        ! Layer elevations
        Me%CEQUALW2%EL(:,:) = 0.0

cl01 :  DO JW=1,Me%CEQUALW2%NWB
cf01 :      IF (Me%CEQUALW2%ZERO_SLOPE(JW)) THEN
cl02 :          DO I=Me%CEQUALW2%US(Me%CEQUALW2%BS(JW))-1,Me%CEQUALW2%DS(Me%CEQUALW2%BE(JW))+1
                    Me%CEQUALW2%EL(Me%CEQUALW2%KMX,I) = Me%CEQUALW2%ELBOT(JW)
cl03 :              DO K=Me%CEQUALW2%KMX-1,1,-1
                        Me%CEQUALW2%EL(K,I) = Me%CEQUALW2%EL(K+1,I)+Me%CEQUALW2%H(K,JW)
                    ENDDO                   cl03
                ENDDO                       cl02
            ELSE                            cf01
                Me%CEQUALW2%NPOINT          = 0                                                                                                       !TC 02/11/01
                NUP                         = 0
                NNBR                        = 1
                NCBP                        = 0
                NINTERNAL                   = 0
                JB                          = Me%CEQUALW2%JBDN(JW)
                Me%CEQUALW2%NPOINT(JB)       = 1
                Me%CEQUALW2%EL(Me%CEQUALW2%KMX,Me%CEQUALW2%DS(JB)+1) = Me%CEQUALW2%ELBOT(JW)
cl90 :          DO WHILE (NNBR <= (Me%CEQUALW2%BE(JW)-Me%CEQUALW2%BS(JW)+1))
                    NCBP = NCBP+1
cf90 :              IF (NCBP > Me%CEQUALW2%NBR) THEN
                        WRITE (*,*) 'Model linkage not logical for waterbody ',JW
                        EXIT
                    ENDIF                   cf90
cf91 :              IF (NINTERNAL == 0) THEN
cf92 :                  IF (NUP == 0) THEN
cl91 :                      DO I=Me%CEQUALW2%DS(JB),Me%CEQUALW2%US(JB)-1,-1
cf93 :                          IF (I /= Me%CEQUALW2%DS(JB)) THEN
                                    Me%CEQUALW2%EL(Me%CEQUALW2%KMX,I) = Me%CEQUALW2%EL(Me%CEQUALW2%KMX,I+1)+Me%CEQUALW2%SLOPE(JB)*(Me%CEQUALW2%DLX(I)+Me%CEQUALW2%DLX(I+1))*0.5
                                ELSE        cf93
                                    Me%CEQUALW2%EL(Me%CEQUALW2%KMX,I) = Me%CEQUALW2%EL(Me%CEQUALW2%KMX,I+1)                                                                                !bug?
                                ENDIF       cf93
cl92 :                          DO K=Me%CEQUALW2%KMX-1,1,-1
                                    Me%CEQUALW2%EL(K,I) = Me%CEQUALW2%EL(K+1,I)+Me%CEQUALW2%H(K,JW)*COS(ATAN(Me%CEQUALW2%SLOPE(JB)))
                                ENDDO       cl92
                            ENDDO           cl91
                        ELSE                cf92
cl93 :                      DO I=Me%CEQUALW2%US(JB),Me%CEQUALW2%DS(JB)+1
cf94 :                          IF (I /= Me%CEQUALW2%US(JB)) THEN
                                    Me%CEQUALW2%EL(Me%CEQUALW2%KMX,I) = Me%CEQUALW2%EL(Me%CEQUALW2%KMX,I-1)-Me%CEQUALW2%SLOPE(JB)*(Me%CEQUALW2%DLX(I)+Me%CEQUALW2%DLX(I+1))*0.5
                                ELSE        cf94
                                    Me%CEQUALW2%EL(Me%CEQUALW2%KMX,I) = Me%CEQUALW2%EL(Me%CEQUALW2%KMX,I-1)
                                ENDIF       cf94
cl94 :                          DO K=Me%CEQUALW2%KMX-1,1,-1
                                    Me%CEQUALW2%EL(K,I) = Me%CEQUALW2%EL(K+1,I)+Me%CEQUALW2%H(K,JW)*COS(ATAN(Me%CEQUALW2%SLOPE(JB)))
                                ENDDO       cl94
                            ENDDO           cl93
                            NUP = 0
                        ENDIF               cf92
cl95 :                  DO K=Me%CEQUALW2%KMX,1,-1
cf95 :                      IF (.NOT.Me%CEQUALW2%DAM_FLOW(JB)) THEN
cf96 :                          IF (Me%CEQUALW2%UHS(JB) /= 0) THEN
                                    Me%CEQUALW2%EL(K,Me%CEQUALW2%US(JB)-1) = Me%CEQUALW2%EL(K,Me%CEQUALW2%US(JB))+Me%CEQUALW2%SLOPE(JB)*Me%CEQUALW2%DLX(Me%CEQUALW2%US(JB))
                                ELSE        cf96
                                    Me%CEQUALW2%EL(K,Me%CEQUALW2%US(JB)-1) = Me%CEQUALW2%EL(K,Me%CEQUALW2%US(JB))
                                ENDIF       cf96
                            ELSE            cf95
                                Me%CEQUALW2%EL(K,Me%CEQUALW2%US(JB)-1) = Me%CEQUALW2%EL(K,Me%CEQUALW2%US(JB))
                            ENDIF           cf95
cf97 :                      IF (Me%CEQUALW2%DHS(JB) /= 0) THEN
                                Me%CEQUALW2%EL(K,Me%CEQUALW2%DS(JB)+1) = Me%CEQUALW2%EL(K,Me%CEQUALW2%DS(JB))-Me%CEQUALW2%SLOPE(JB)*Me%CEQUALW2%DLX(Me%CEQUALW2%DS(JB))
                            ELSE            cf97
                                Me%CEQUALW2%EL(K,Me%CEQUALW2%DS(JB)+1) = Me%CEQUALW2%EL(K,Me%CEQUALW2%DS(JB))
                            ENDIF           cf97
                        ENDDO               cl95
                    ELSE                    cf91
cl96 :                  DO K=Me%CEQUALW2%KMX-1,1,-1
                            Me%CEQUALW2%EL(K,Me%CEQUALW2%UHS(JJB)) = Me%CEQUALW2%EL(K+1,Me%CEQUALW2%UHS(JJB))+Me%CEQUALW2%H(K,JW)*COS(ATAN(Me%CEQUALW2%SLOPE(JB)))
                        ENDDO               cl96    
cl97 :                  DO I = Me%CEQUALW2%UHS(JJB)+1, Me%CEQUALW2%DS(JB)
                            Me%CEQUALW2%EL(Me%CEQUALW2%KMX,I) = Me%CEQUALW2%EL(Me%CEQUALW2%KMX,I-1)-Me%CEQUALW2%SLOPE(JB)*(Me%CEQUALW2%DLX(I)+Me%CEQUALW2%DLX(I-1))*0.5
cl98 :                      DO K=Me%CEQUALW2%KMX-1,1,-1
                                Me%CEQUALW2%EL(K,I) = Me%CEQUALW2%EL(K+1,I)+Me%CEQUALW2%H(K,JW)*COS(ATAN(Me%CEQUALW2%SLOPE(JB)))
                            ENDDO           cl98
                        ENDDO               cl97
cl99 :                  DO I=Me%CEQUALW2%UHS(JJB)-1,Me%CEQUALW2%US(JB),-1
                            Me%CEQUALW2%EL(Me%CEQUALW2%KMX,I) = Me%CEQUALW2%EL(Me%CEQUALW2%KMX,I+1)+Me%CEQUALW2%SLOPE(JB)*(Me%CEQUALW2%DLX(I)+Me%CEQUALW2%DLX(I+1))*0.5
cl100 :                     DO K=Me%CEQUALW2%KMX-1,1,-1
                                Me%CEQUALW2%EL(K,I) = Me%CEQUALW2%EL(K+1,I)+Me%CEQUALW2%H(K,JW)*COS(ATAN(Me%CEQUALW2%SLOPE(JB)))
                            ENDDO           cl100
                        ENDDO               cl99
                        NINTERNAL = 0
                    ENDIF                   cf91
                    IF (NNBR == (Me%CEQUALW2%BE(JW)-Me%CEQUALW2%BS(JW)+1)) EXIT
                    NNBROLD = NNBR
cl101 :             DO JB=Me%CEQUALW2%BS(JW),Me%CEQUALW2%BE(JW)
cf98  :                 IF (Me%CEQUALW2%NPOINT(JB) /= 1) THEN
cl102 :                     DO JJB = Me%CEQUALW2%BS(JW), Me%CEQUALW2%BE(JW)
cf99  :                         IF (Me%CEQUALW2%DHS(JB) >= Me%CEQUALW2%US(JJB) .AND. Me%CEQUALW2%DHS(JB) <= Me%CEQUALW2%DS(JJB) .AND. Me%CEQUALW2%NPOINT(JJB) == 1) THEN
                                    Me%CEQUALW2%NPOINT(JB)       = 1
                                    NNBR             = NNBR+1
                                    Me%CEQUALW2%EL(Me%CEQUALW2%KMX,Me%CEQUALW2%DS(JB)+1) = Me%CEQUALW2%EL(Me%CEQUALW2%KMX,Me%CEQUALW2%DHS(JB))+Me%CEQUALW2%SLOPE(JB)*(Me%CEQUALW2%DLX(Me%CEQUALW2%DS(JB))+Me%CEQUALW2%DLX(Me%CEQUALW2%DHS(JB)))*0.5
                                    EXIT
                                ENDIF       cf99
cf100 :                         IF (Me%CEQUALW2%UHS(JJB) == Me%CEQUALW2%DS(JB) .AND. Me%CEQUALW2%NPOINT(JJB) == 1) THEN
                                    Me%CEQUALW2%NPOINT(JB)       = 1
                                    NNBR             = NNBR+1
                                    Me%CEQUALW2%EL(Me%CEQUALW2%KMX,Me%CEQUALW2%DS(JB)+1) = Me%CEQUALW2%EL(Me%CEQUALW2%KMX,Me%CEQUALW2%US(JJB))+Me%CEQUALW2%SLOPE(JB)*(Me%CEQUALW2%DLX(Me%CEQUALW2%US(JJB))+Me%CEQUALW2%DLX(Me%CEQUALW2%DS(JB)))*0.5
                                    EXIT
                                ENDIF       cf100
cf101 :                         IF (Me%CEQUALW2%UHS(JJB) <= Me%CEQUALW2%DS(JB) .AND. Me%CEQUALW2%UHS(JJB) >= Me%CEQUALW2%US(JB) .AND. Me%CEQUALW2%NPOINT(JJB) == 1) THEN
                                    NNBR             = NNBR+1
                                    NINTERNAL        = 1
                                    Me%CEQUALW2%NPOINT(JB)       = 1
                                    Me%CEQUALW2%EL(Me%CEQUALW2%KMX,Me%CEQUALW2%UHS(JJB)) = Me%CEQUALW2%EL(Me%CEQUALW2%KMX,Me%CEQUALW2%US(JJB))+Me%CEQUALW2%SLOPE(JJB)*Me%CEQUALW2%DLX(Me%CEQUALW2%US(JJB))*0.5
                                    EXIT
                                ENDIF      cf101
cf102 :                         IF (Me%CEQUALW2%UHS(JB) <= Me%CEQUALW2%DS(JJB) .AND. Me%CEQUALW2%UHS(JB) >= Me%CEQUALW2%US(JJB) .AND. Me%CEQUALW2%NPOINT(JJB) == 1) THEN
                                    NNBR             = NNBR+1
                                    Me%CEQUALW2%NPOINT(JB)       = 1
                                    Me%CEQUALW2%EL(Me%CEQUALW2%KMX,Me%CEQUALW2%US(JB)-1) = Me%CEQUALW2%EL(Me%CEQUALW2%KMX,Me%CEQUALW2%UHS(JB))-Me%CEQUALW2%SLOPE(JB)*Me%CEQUALW2%DLX(Me%CEQUALW2%US(JB))*0.5
cf103 :                             IF (Me%CEQUALW2%SLOPE(JB) /= 0.0) THEN
                                        WRITE (*,*) 'Branch ',JB,' has non-zero slope and is only'
                                        WRITE (*,*) 'connected through a UHS BC. This may not Me%CEQUALW2%BE'
                                        WRITE (*,*) 'a problem:if it is,change branch to a waterbody'                                   !TC 02/09/01
                                    ENDIF   cf103
                                    NUP = 1
                                    EXIT
                                ENDIF       cf102
                            ENDDO           cl102
                            IF (Me%CEQUALW2%NPOINT(JB) == 1) EXIT
                        ENDIF                   cf98
                    ENDDO                   cl101
cf104 :             IF (NCBP > Me%CEQUALW2%NBR .OR. NNBROLD == NNBR) THEN
                        WRITE (*,*) 'Error in grid linkage - check [US], [DS], [UHS], [DHS], and or [JBDN] for waterbody ',JW
                    ENDIF                   cf104
                ENDDO                       cl90
            ENDIF                           cf01
        ENDDO                               cl01

        ! Grid geometry
        Me%CEQUALW2%KB  (:) = 0
        Me%CEQUALW2%KTWB(:) = 2


cl10 :  DO JW=1,Me%CEQUALW2%NWB
        ZMIN =-1000.0
cl15 :      DO JB=Me%CEQUALW2%BS(JW),Me%CEQUALW2%BE(JW)
                Me%CEQUALW2%ALPHA(JB) = ATAN(Me%CEQUALW2%SLOPE(JB)) 
                IFLAG     = 0
cl20 :          DO I=Me%CEQUALW2%US(JB)-1,Me%CEQUALW2%DS(JB)+1
                    Me%CEQUALW2%KTI(I) = 2
cl25 :              DO WHILE (Me%CEQUALW2%EL(Me%CEQUALW2%KTI(I),I) > Me%CEQUALW2%ELWS(I) .AND. Me%CEQUALW2%KTI(I) < Me%CEQUALW2%KMX)
                        Me%CEQUALW2%KTI(I) = Me%CEQUALW2%KTI(I)+1
cf10 :                  IF (Me%CEQUALW2%KTI(I) > Me%CEQUALW2%KMX-1) THEN
                            WRITE (*,*) 'Initial water surface elevation below grid at segment ',I
                            STOP        'PopulateVariables - ModuleReservoirOptimization - ERR10'
                        ENDIF                           cf10
                    ENDDO                               cl25
                    Me%CEQUALW2%Z(I)        = (Me%CEQUALW2%EL(Me%CEQUALW2%KTI(I),I)-Me%CEQUALW2%ELWS(I))/COS(Me%CEQUALW2%ALPHA(JB))
                    ZMIN                    =  MAX(ZMIN,Me%CEQUALW2%Z(I))
                    KTMAX                   =  MAX(2,       Me%CEQUALW2%KTI (I ))
                    Me%CEQUALW2%KTWB(JW)    =  MAX(KTMAX,   Me%CEQUALW2%KTWB(JW))
                    Me%CEQUALW2%KTI (I )    =  MAX(Me%CEQUALW2%KTI(I)-1,2                   )
cf15 :              IF (Me%CEQUALW2%KTWB(JW) >= Me%CEQUALW2%KMX .AND. IFLAG == 0) THEN
                        IFLAG =  1
                        WRITE (*,*) 'Water surface elevation too low in branch ',JB,' at segment ',I
                    ENDIF                               cf15
cl26 :              DO K=2,Me%CEQUALW2%KMX-1
cf26 :              IF (Me%CEQUALW2%B(K+1,I) == 0.0 .OR. K == Me%CEQUALW2%KMX-1) THEN
                        Me%CEQUALW2%KB(I) = K
                        EXIT cl26
                    ENDIF                               cf26
                    ENDDO                               cl26
                ENDDO                                   cl20
                Me%CEQUALW2%KB(Me%CEQUALW2%US(JB)-1) = Me%CEQUALW2%KB(Me%CEQUALW2%US(JB))
                Me%CEQUALW2%KB(Me%CEQUALW2%DS(JB)+1) = Me%CEQUALW2%KB(Me%CEQUALW2%DS(JB))
            ENDDO                                       cl15
        ENDDO                                       cl10


        ! Boundary cells
cl300 : DO JB=1,Me%CEQUALW2%NBR
cl35 :      DO K=1,Me%CEQUALW2%KMX
                IF (Me%CEQUALW2%B(K,Me%CEQUALW2%US(JB)-1) /= 0.0)                                           &
                    WRITE (*,*) 'Upstream boundary segment width [Me%CEQUALW2%B(',                          &
                                    K,',',Me%CEQUALW2%US(JB)-1,')=',Me%CEQUALW2%B(K,Me%CEQUALW2%US(JB)-1),  &
                                    '] /= 0 for branch ',JB

                IF (Me%CEQUALW2%B(K,Me%CEQUALW2%DS(JB)+1) /= 0.0)                                           &
                    WRITE (*,*) 'Downstream boundary segment width [Me%CEQUALW2%B(',                        &
                                    K,',',Me%CEQUALW2%DS(JB)+1,')=',Me%CEQUALW2%B(K,Me%CEQUALW2%DS(JB)+1),  &
                                    '] /= 0 for branch ',JB
            ENDDO                                                               cl35
cl400 :     DO I=Me%CEQUALW2%US(JB),Me%CEQUALW2%DS(JB)
cl45 :          DO K=2,Me%CEQUALW2%KB(I)
                    IF (Me%CEQUALW2%B(K+1,I) > Me%CEQUALW2%B(K,I) .AND. Me%CEQUALW2%B(K+1,I) /= 0.0)    &
                        WRITE (*,*) 'Cell width [Me%CEQUALW2%B(',K+1,',',I,')=',Me%CEQUALW2%B(K+1,I),   &
                                        '] > [Me%CEQUALW2%B(',K,',',I,')=',Me%CEQUALW2%B(K,I)

                    IF (Me%CEQUALW2%B(K,I) < 5.0 .AND. Me%CEQUALW2%B(K,I) > 0.0)            &
                        WRITE (*,*) 'Cell width [Me%CEQUALW2%B(',                           &
                                        K,',',I,')=',Me%CEQUALW2%B(K,I),                    &
                                        '] < 5m which can cause stability problems'
                ENDDO                                                           cl45
cf45 :          IF      (Me%CEQUALW2%B(1,              I) /= 0.0)   THEN
                    WRITE (*,*) 'Surface boundary layer cell width [Me%CEQUALW2%B(1,',I,')=',Me%CEQUALW2%B(1,I),'] /= 0'

                ELSEIF  (Me%CEQUALW2%B(Me%CEQUALW2%KMX,I) /= 0.0)   THEN        cf45
                    WRITE (*,*) 'Bottom boundary layer cell width [Me%CEQUALW2%B(',         &
                                    Me%CEQUALW2%KMX,',',I,')=',Me%CEQUALW2%B(Me%CEQUALW2%KMX,I),'] /= 0'
                ENDIF                                                           cf45
            ENDDO                                                               cl400
        ENDDO                                                                   cl300



        ! Boundary widths
cl60 :  DO JB=1,Me%CEQUALW2%NBR
cl61 :  DO K =2,Me%CEQUALW2%KMX
            Me%CEQUALW2%B(K,Me%CEQUALW2%US(JB)-1) = Me%CEQUALW2%B(K,Me%CEQUALW2%US(JB))
            Me%CEQUALW2%B(K,Me%CEQUALW2%DS(JB)+1) = Me%CEQUALW2%B(K,Me%CEQUALW2%DS(JB)) 
        ENDDO                                           cl61
        ENDDO                                           cl60
  


        ! Width change
cl70 :  DO JB=1,Me%CEQUALW2%NBR
cl71 :  DO I=Me%CEQUALW2%US(JB),Me%CEQUALW2%DS(JB)-1
cl72 :  DO K=2,Me%CEQUALW2%KB(I)
cf71 :  IF (Me%CEQUALW2%B(K,I+1) /= 0.0 .AND. Me%CEQUALW2%B(K,I) /= 0.0) THEN
            IF (Me%CEQUALW2%B(K,I)/Me%CEQUALW2%B(K,I+1) >= 7.0 .OR. Me%CEQUALW2%B(K,I+1)/Me%CEQUALW2%B(K,I) >= 7.0)     &
                WRITE (*,*) 'Cell width [Me%CEQUALW2%B(',K,',',I,')=',Me%CEQUALW2%B(K,I),                               &
                                '] < or > 7x width [Me%CEQUALW2%B(',K,',',I+1,')=',Me%CEQUALW2%B(K,I+1),']'
        ENDIF                                           cf71
        ENDDO                                           cl72
        ENDDO                                           cl71
        ENDDO                                           cl70



    ! Layer thickness
cl80 :  DO JW=1,Me%CEQUALW2%NWB
cl81 :  DO K =1,Me%CEQUALW2%KMX
cf81 :      IF      (Me%CEQUALW2%H(K,JW) <= 0 ) THEN
                WRITE (*,*) 'Layer thickness [Me%CEQUALW2%H =',Me%CEQUALW2%H(K,JW),         &
                            '] <= 0.0 m for layer ',K,' in waterbody ',JW
            ELSEIF  (Me%CEQUALW2%H(K,JW) < 0.1) THEN
                WRITE (*,*) 'Layer thickness [Me%CEQUALW2%H =',Me%CEQUALW2%H(K,JW),         &
                            '] <  0.1 m for layer ',K,' in waterbody ',JW
            ELSEIF  (Me%CEQUALW2%H(K,JW) > 3.0) THEN
                WRITE (*,*) 'Layer thickness [Me%CEQUALW2%H =',Me%CEQUALW2%H(K,JW),         &
                            '] >  3.0 m for layer ',K,' in waterbody ',JW
            ENDIF                                   cf81
        ENDDO                                       cl81
        ENDDO                                       cl80

       
        !----------------------------------------------------------------------
        
    end subroutine PopulateVariables

    !--------------------------------------------------------------------------
    
    subroutine ReadFillShape()
        
        !Arguments-------------------------------------------------------------
    
        !External--------------------------------------------------------------
    
        integer :: STAT_CALL 
        integer :: ClientNumber
        integer :: FirstLine, LastLine
        logical :: BlockFound

        !Local-----------------------------------------------------------------

        integer :: ObjEnterData03 = 0

        !----------------------------------------------------------------------

        call ConstructEnterData(ObjEnterData03, Me%Files%FillShape, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFillShape - ModuleReservoirOptimization - ERR10'



        call ExtractBlockFromBuffer(ObjEnterData03, ClientNumber,           &
                                                    '<Fill>', '</Fill>',    &
                                                    BlockFound,             &
                                                    FirstLine = FirstLine,  &
                                                    LastLine  = LastLine,   &
                                                    STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFillShape - ModuleReservoirOptimization - ERR20'



        call FillShapeMatrix(ObjEnterData03, FirstLine, LastLine)



        call KillEnterData     (ObjEnterData03, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFillShape - ModuleReservoirOptimization - ERR90'
                        
        !----------------------------------------------------------------------
        
    end subroutine ReadFillShape

    !--------------------------------------------------------------------------
    
    subroutine FillShapeMatrix(ObjEnterData, FirstLine, LastLine)
        
        !Arguments-------------------------------------------------------------
    
        integer, intent(IN) :: ObjEnterData
        integer, intent(IN) :: FirstLine    
        integer, intent(IN) :: LastLine    

        !External--------------------------------------------------------------
    
        integer             :: STAT_CALL 

        !Local-----------------------------------------------------------------

        integer             :: K, JW
        real(4)             :: Height

        !----------------------------------------------------------------------

        nullify (Me%FillData%Vol                         )
        nullify (Me%FillData%VolDiff                     )
        nullify (Me%FillData%VolError                    )
       
        allocate(Me%FillData%Vol     (Me%CEQUALW2%KMX   ))
        allocate(Me%FillData%VolDiff (Me%CEQUALW2%KMX   ))
        allocate(Me%FillData%VolError(Me%CEQUALW2%KMX   ))

        Me%FillData%Vol(:)  = null_real

cl09 :  do JW= 1, Me%CEQUALW2%NWB
            Height = 0

cl10 :      do K=Me%CEQUALW2%KMX-1,2,-1
                Height = Height + Me%CEQUALW2%H(K,JW)

                call CalcVol(ObjEnterData,                                  &
                                Me%FillData%Vol (K  ),                      &
                                Height,                                     &
                                FirstLine,                                  &
                                LastLine)                                   
            enddo cl10
        enddo cl09
        !----------------------------------------------------------------------
        
    end subroutine FillShapeMatrix

    !--------------------------------------------------------------------------
    
    subroutine CalcVol(ObjEnterData,                                    &
                        Volume,                                         &
                        H,                                              &
                        FirstLine,                                      &
                        LastLine)
        
        !Arguments-------------------------------------------------------------
    
        integer, intent(IN )    :: ObjEnterData
        integer, intent(IN )    :: FirstLine    
        integer, intent(IN )    :: LastLine    
        real(4), intent(IN )    :: H
        real(4), intent(OUT)    :: Volume

        !External--------------------------------------------------------------
    
        integer                 :: STAT_CALL 
        integer                 :: flag
        real(4), dimension(2)   :: vector

        !Local-----------------------------------------------------------------

        integer                 :: colH     = 1
        integer                 :: colVol   = 2
        integer                 :: line
        real(4)                 :: H1,      H2
        real(4)                 :: Volume1, Volume2
        logical                 :: Found

        !----------------------------------------------------------------------

        H1      = null_real
        H2      = null_real
        Volume  = null_real
        Volume1 = null_real
        Volume2 = null_real
        Found   = .false.

        if (((LastLine-1) - (FirstLine+1)) < 1)                         &
            stop 'Dam fill data needs at least 2 values - CalcVol - ModuleReservoirOptimization - ERR20'

cl10 :  do line = (FirstLine+1), (LastLine-2)
            call GetData(vector,                                        &
                            ObjEnterData,                               &
                            flag        = flag,                         &
                            Buffer_Line = line,                         &
                            STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CalcVol - ModuleReservoirOptimization - ERR10'
            H1      = vector(colH   )
            Volume1 = vector(colVol )

            call GetData(vector,                                        &
                            ObjEnterData,                               &
                            flag        = flag,                         &
                            Buffer_Line = line+1,                       &
                            STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CalcVol - ModuleReservoirOptimization - ERR15'
            H2      = vector(colH   )
            Volume2 = vector(colVol )



            !It is not necessary to interpolate/extrapolate
cs10 :      if          (H == H1) then
                Volume  = Volume1
                Found   = .true.
                exit cl10

            elseif      (H == H2) then  cs10
                Volume  = Volume2
                Found   = .true.
                exit cl10
            endif                       cs10




cs15 :      if      ((H1 < H) .AND. (H < H2)) then         !H1<H<H2
                call InterpolateVol(Volume,     H,                                  &
                                    H1,         H2,                                 &
                                    Volume1,    Volume2)
                Found   = .true.
                exit cl10

            elseif  ( (H < H1) .AND. (H < H2)                               ) then  cs15   !H < file's data
                call InterpolateVol(Volume,     H,                                  &
                                    H1,         H2,                                 &
                                    Volume1,    Volume2)
                Found   = .true.
                exit cl10

            elseif  (((H1 < H) .AND. (H2 < H)) .AND. (line == (LastLine-2)) ) then  cs15   !file's data < H
                call InterpolateVol(Volume,     H,                                  &
                                    H1,         H2,                                 &
                                    Volume1,    Volume2)
                Found   = .true.
                exit cl10
            endif                                   cs15
        enddo cl10
        
        if (.NOT. Found)                                                        &
            stop 'Dam Volume not found - CalcVol - ModuleReservoirOptimization - ERR90'
        
        !----------------------------------------------------------------------
        
    end subroutine CalcVol

    !--------------------------------------------------------------------------
    
    subroutine InterpolateVol(Volume,   H,                                      &
                                H1,     H2,                                     &
                                Vol1,   Vol2)
        
        !Arguments-------------------------------------------------------------
    
        real(4), intent(OUT)    :: Volume
        real(4), intent(IN )    :: H
        real(4), intent(IN )    :: H1,      H2
        real(4), intent(IN )    :: Vol1,    Vol2

        !External--------------------------------------------------------------
    
        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        Volume = Vol1+(Vol2-Vol1)/(H2-H1)*(H-H1)

        !----------------------------------------------------------------------
        
    end subroutine InterpolateVol

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    !--------------------------------------------------------------------------
    subroutine GetReservoirOptimizationPointer (ObjReservoirOptimizationID, Matrix, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjReservoirOptimizationID
        real(8), dimension(:, :, :),  pointer           :: Matrix
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirOptimizationID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRESERVOIROPTIMIZATION_, Me%InstanceID)

            Matrix => Me%Matrix

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        endif

        if (present(STAT)) STAT = STAT_

    end subroutine GetReservoirOptimizationPointer
    
    !--------------------------------------------------------------------------
    
    subroutine GetReservoirOptimizationInteger (ObjReservoirOptimizationID, Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjReservoirOptimizationID
        real                                            :: Int
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirOptimizationID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Int = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        endif

        if (present(STAT)) STAT = STAT_

    end subroutine GetReservoirOptimizationInteger

    !--------------------------------------------------------------------------

    subroutine UnGetReservoirOptimization3D_I(ObjReservoirOptimizationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjReservoirOptimizationID
        integer, dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirOptimizationID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mRESERVOIROPTIMIZATION_, Me%InstanceID, "UnGetReservoirOptimization3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        endif

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetReservoirOptimization3D_I

    !--------------------------------------------------------------------------

    subroutine UnGetReservoirOptimization3D_R8(ObjReservoirOptimizationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjReservoirOptimizationID
        real(8), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirOptimizationID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mRESERVOIROPTIMIZATION_, Me%InstanceID,  "UnGetReservoirOptimization3D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        endif

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetReservoirOptimization3D_R8

    !--------------------------------------------------------------------------

    subroutine WriteW2Bathymetry(ObjReservoirOptimizationID, STAT)
 
        !Arguments-------------------------------------------------------------

        integer                                     :: ObjReservoirOptimizationID
        integer, intent(OUT), optional              :: STAT

        !External------------------------------------------------------------

        integer                         :: STAT_CALL 
        integer                         :: IERR
        integer                         :: nunit

        !Local-----------------------------------------------------------------

        integer                         :: I, K
        integer                         :: JW, JB
        character(len = StringLength)   :: Message
        integer                         :: STAT_, ready_

        !----------------------------------------------------------------------
        
        STAT_ = UNKNOWN_

        call Ready(ObjReservoirOptimizationID, ready_)

        ! Boundary widths
cl60 :  DO JB=1,Me%CEQUALW2%NBR
cl61 :  DO K =2,Me%CEQUALW2%KMX
            Me%CEQUALW2%B(K,Me%CEQUALW2%US(JB)-1) = 0.0
            Me%CEQUALW2%B(K,Me%CEQUALW2%DS(JB)+1) = 0.0
        ENDDO                                           cl61
        ENDDO                                           cl60


cf3 :   if (ready_ .EQ. IDLE_ERR_) then
cl1 :       DO JW=1,Me%CEQUALW2%NWB
                call UnitsManager(nunit, OPEN_FILE, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'WriteW2Bathymetry - ModuleReservoirOptimization - ERR10'

                OPEN (nunit,FILE=Me%Files%ModifiedW2Bathymetry,STATUS='UNKNOWN',IOSTAT=IERR)
cf1 :           IF (IERR == 0) THEN
                    WRITE (nunit,'(A256)')                                                              &
                        Me%CEQUALW2%CommentHeader(JW)//" - Modified"
                    !'  segment lengths'
                    WRITE (nunit,'(/A256/(10F8.0))')                                                    &
                        Me%CEQUALW2%CommentDLX(JW),(Me%CEQUALW2%DLX(I),                                 &
                        I=Me%CEQUALW2%US(Me%CEQUALW2%BS(JW))-1,Me%CEQUALW2%DS(Me%CEQUALW2%BE(JW))+1)
                    !'  water surface elevations'
                    WRITE (nunit,'(/A256/(10F8.0))')                                                    &
                        Me%CEQUALW2%CommentELWS(JW),(Me%CEQUALW2%ELWS(I),                               &
                        I=Me%CEQUALW2%US(Me%CEQUALW2%BS(JW))-1,Me%CEQUALW2%DS(Me%CEQUALW2%BE(JW))+1)
                    !'  segment orientation'
                    WRITE (nunit,'(/A256/(10F8.0))')                                                    &
                        Me%CEQUALW2%CommentPHI0(JW),(Me%CEQUALW2%PHI0(I),                               &
                        I=Me%CEQUALW2%US(Me%CEQUALW2%BS(JW))-1,Me%CEQUALW2%DS(Me%CEQUALW2%BE(JW))+1)
                    !'  segment bottom friction'
                    WRITE (nunit,'(/A256/(10F8.0))')                                                    &
                        Me%CEQUALW2%CommentFRIC(JW),(Me%CEQUALW2%FRIC(I),                               &
                        I=Me%CEQUALW2%US(Me%CEQUALW2%BS(JW))-1,Me%CEQUALW2%DS(Me%CEQUALW2%BE(JW))+1)
                    !'  layer thickness'
                    WRITE (nunit,'(/A256/(10F8.0))')                                                    &
                        Me%CEQUALW2%CommentH(JW),(Me%CEQUALW2%H(K,JW), K=1,Me%CEQUALW2%KMX)
                    !'  segment widths'
cl2 :               DO I=Me%CEQUALW2%US(Me%CEQUALW2%BS(JW))-1,Me%CEQUALW2%DS(Me%CEQUALW2%BE(JW))+1
                        WRITE(nunit,'(/A256/(10F8.0))')                                                 &
                            Me%CEQUALW2%CommentB(JW,I),(Me%CEQUALW2%B(K,I),  K=1,Me%CEQUALW2%KMX)
                    ENDDO cl2
                
                    call UnitsManager          (nunit, CLOSE_FILE, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'WriteW2Bathymetry - ModuleReservoirOptimization - ERR20' 
                ELSE    cf1
                    write(*,*)  'Could not open '//Me%Files%ModifiedW2Bathymetry
                    stop        'WriteW2Bathymetry - ModuleReservoirOptimization - ERR30'
                ENDIF  cf1
            ENDDO cl1

            STAT_ = SUCCESS_
        else                            cf3
            STAT_ = ready_
        endif                           cf3

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------
        
    end subroutine WriteW2Bathymetry

    !------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyReservoirOptimization(ObjReservoirOptimizationID, STAT)

        !Arguments-------------------------------------------------------------

        integer                                     :: ObjReservoirOptimizationID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------

        integer                                     :: STAT_, ready_
        logical                                     :: CorrectVol

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirOptimizationID, ready_)

cf10 :  if (ready_ .EQ. IDLE_ERR_) then
            call CalcCurrentVolume  (           )
            call CalcDiffVolume     (           )
            call RecalculateVolume  (CorrectVol )
            call Verify             (           )

cf11 :      if (.NOT. CorrectVol) then
                STAT_ = UNKNOWN_
            else            cf11
                STAT_ = SUCCESS_
            endif           cf11
        else                cf10
            STAT_ = ready_
        endif               cf10

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyReservoirOptimization

    !--------------------------------------------------------------------------
    
    subroutine CalcCurrentVolume()
        
        !Arguments-------------------------------------------------------------
    
        !External--------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
        
        integer :: JW, JB, IU, ID, IUT, I, K
        real(4) :: DIFF

        !----------------------------------------------------------------------

        Me%CEQUALW2%VOLB (:)    = 0.0
        Me%CEQUALW2%VOLG (:)    = 0.0
        Me%CEQUALW2%CVBR (:,:)  = 0.0
        Me%CEQUALW2%CVGR (:,:)  = 0.0
        Me%CEQUALW2%SAGR (:,:)  = 0.0
        Me%CEQUALW2%SABR (:,:)  = 0.0
        Me%CEQUALW2%NCCGR(:,:)  = 0
        Me%CEQUALW2%NCCBR(:,:)  = 0

cl10 :  DO JW=1,Me%CEQUALW2%NWB
cl11 :  DO JB=Me%CEQUALW2%BS(JW),Me%CEQUALW2%BE(JW)
            IU  = Me%CEQUALW2%US(JB)
            ID  = Me%CEQUALW2%DS(JB)
cf11 :      IF (.NOT.Me%CEQUALW2%DAM_FLOW(JB)) THEN
cf12 :      IF (Me%CEQUALW2%UHS(JB) /= 0 .AND. Me%CEQUALW2%UHS(JB) > 1) THEN
                DIFF = ABS(Me%CEQUALW2%ELWS(IU)-Me%CEQUALW2%ELWS(Me%CEQUALW2%UHS(JB)))
                IF (DIFF > 1.5)                                                         &
                    WRITE (*,*) 'Water surface elevation difference of',DIFF,           &  
                                ' m between segment ',IU,'and segment ',Me%CEQUALW2%UHS(JB)
            ENDIF                                   cf12
            ENDIF                                   cf11
cf13 :      IF (Me%CEQUALW2%DHS(JB) /= 0 .AND. Me%CEQUALW2%DHS(JB) > 1) THEN
                DIFF = ABS(Me%CEQUALW2%ELWS(ID)-Me%CEQUALW2%ELWS(Me%CEQUALW2%DHS(JB)))
                IF (DIFF > 1.5)                                                         &
                    WRITE (*,*) 'Water surface elevation difference of ',DIFF,          &
                                ' m between segment ',ID,'and segment ', Me%CEQUALW2%DHS(JB)
            ENDIF                                   cf13



            !**** Water surface and bottom layers
cl20 :      DO I=IU-1,ID+1
cf20 :          IF      (Me%CEQUALW2%ELWS(I) <= Me%CEQUALW2%EL(Me%CEQUALW2%KB(I)+1,I)) THEN
                    WRITE (*,*) 'Water surface elevation is below bottom elevation at segment ',I
                ELSEIF  (Me%CEQUALW2%ELWS(I) <=  Me%CEQUALW2%EL(Me%CEQUALW2%KB(I)+1,I)+ &
                         0.4*Me%CEQUALW2%H(Me%CEQUALW2%KB(I),JW)) THEN    cf20
                    WRITE (*,*) 'Water surface elevation is close to bottom elevation at segment ',I
                ENDIF                                                                   cf20
cf21 :          IF (I < ID+1) THEN
                    DIFF = ABS(Me%CEQUALW2%ELWS(I)-Me%CEQUALW2%ELWS(I+1))
                    IF (DIFF > 1.5)                             &
                        WRITE (*,*) 'Water surface elevation difference of ',DIFF,'m between segment ',I,'and segment ',I+1
                ENDIF                                       cf21
            ENDDO                                           cl20



            !**** Branch and grid total volume
cl30 :      DO I=IU,ID
                Me%CEQUALW2%VOLG(JW) = Me%CEQUALW2%VOLG(JW)+Me%CEQUALW2%DLX(I)*Me%CEQUALW2%B(Me%CEQUALW2%KTWB(JW),I)*   &
                                        (Me%CEQUALW2%H(Me%CEQUALW2%KTWB(JW),JW)-Me%CEQUALW2%Z(I))
cl31 :          DO K=Me%CEQUALW2%KTWB(JW)+1,Me%CEQUALW2%KB(I)
                    Me%CEQUALW2%VOLB(JB) = Me%CEQUALW2%VOLB(JB)+Me%CEQUALW2%DLX(I)*Me%CEQUALW2%B(K,I)*Me%CEQUALW2%H(K,JW)
                    Me%CEQUALW2%VOLG(JW) = Me%CEQUALW2%VOLG(JW)+Me%CEQUALW2%DLX(I)*Me%CEQUALW2%B(K,I)*Me%CEQUALW2%H(K,JW)
                ENDDO                               cl31
            ENDDO                                   cl30



            !**** Branch and grid area and volume by layer
cl40 :      DO K=Me%CEQUALW2%KMX-1,2,-1
                Me%CEQUALW2%NCCBR(K,JB) = Me%CEQUALW2%NCCBR(K+1,JB)
                Me%CEQUALW2%CVBR (K,JB) = Me%CEQUALW2%CVBR (K+1,JB)
cl41 :          DO I=IU,ID
cf41 :              IF (K <= Me%CEQUALW2%KB(I)) THEN
                        Me%CEQUALW2%SABR (K,JB) = Me%CEQUALW2%SABR(K,JB)+   &
                                                    Me%CEQUALW2%B(K,I)                    *Me%CEQUALW2%DLX(I)
                        Me%CEQUALW2%CVBR (K,JB) = Me%CEQUALW2%CVBR(K,JB)+   &
                                                    Me%CEQUALW2%B(K,I)*Me%CEQUALW2%H(K,JW)*Me%CEQUALW2%DLX(I)
                        Me%CEQUALW2%NCCBR(K,JB) = Me%CEQUALW2%NCCBR(K,JB)+1
                    ENDIF                               cf41
                ENDDO                                   cl41
                Me%CEQUALW2%SAGR (K,JW) = Me%CEQUALW2%SAGR (K,JW)+Me%CEQUALW2%SABR (K,JB)
                Me%CEQUALW2%NCCGR(K,JW) = Me%CEQUALW2%NCCGR(K,JW)+Me%CEQUALW2%NCCBR(K,JB)
                Me%CEQUALW2%CVGR (K,JW) = Me%CEQUALW2%CVGR (K,JW)+Me%CEQUALW2%CVBR (K,JB)
            ENDDO                                       cl40
        ENDDO                                           cl11
        ENDDO                                           cl10

        !----------------------------------------------------------------------
        
    end subroutine CalcCurrentVolume

    !--------------------------------------------------------------------------
    
    subroutine CalcDiffVolume()
        
        !Arguments-------------------------------------------------------------
    
        !External--------------------------------------------------------------
    
        !Local-----------------------------------------------------------------

        integer :: K, JW
real :: a

        !----------------------------------------------------------------------

        Me%FillData%VolDiff (:) = null_real
        Me%FillData%VolError(:) = null_real

cl11:   DO JW   =1,Me%CEQUALW2%NWB
cl10 :  do K    =Me%CEQUALW2%KMX-1,2,-1
            Me%FillData%VolDiff (K) = Me%FillData%Vol(K) - Me%CEQUALW2%CVGR(K,JW)
            Me%FillData%VolError(K) = Me%FillData%Vol(K) / Me%CEQUALW2%CVGR(K,JW)
        enddo                                               cl10
        enddo                                               cl11

        !----------------------------------------------------------------------
        
    end subroutine CalcDiffVolume

    !--------------------------------------------------------------------------
    
    subroutine RecalculateVolume(CorrectVol)
        
        !Arguments-------------------------------------------------------------

        logical, intent(OUT) :: CorrectVol          
    
        !External--------------------------------------------------------------
    
        !Local-----------------------------------------------------------------

        integer :: JW, JB, IU, ID, I, K

        !----------------------------------------------------------------------

        CorrectVol =.false.


cl10 :  DO JW=1,Me%CEQUALW2%NWB
cl11 :  DO JB=Me%CEQUALW2%BS(JW),Me%CEQUALW2%BE(JW)
        IU  = Me%CEQUALW2%US(JB)
        ID  = Me%CEQUALW2%DS(JB)
        !**** Branch and grid area and volume by layer
cl40 :  DO K=Me%CEQUALW2%KMX-1,2,-1
cf40 :  if (ABS(Me%FillData%VolError(K)) > Me%FillData%DamVolError) then
            CorrectVol =.true.

cl41 :      DO I=IU,ID
cf41 :      IF (K <= Me%CEQUALW2%KB(I)) THEN
                Me%CEQUALW2%B(K,I)      = Me%CEQUALW2%B(K,I)*Me%FillData%VolError(K)
            ENDIF                                           cf41
            ENDDO                                           cl41
        endif                                               cf40
        ENDDO                                               cl40
        ENDDO                                               cl11
        ENDDO                                               cl10




        ! Boundary widths
cl60 :  DO JB=1,Me%CEQUALW2%NBR
cl61 :  DO K =2,Me%CEQUALW2%KMX
            Me%CEQUALW2%B(K,Me%CEQUALW2%US(JB)-1) = Me%CEQUALW2%B(K,Me%CEQUALW2%US(JB))
            Me%CEQUALW2%B(K,Me%CEQUALW2%DS(JB)+1) = Me%CEQUALW2%B(K,Me%CEQUALW2%DS(JB)) 
        ENDDO                                           cl61
        ENDDO                                           cl60

        !----------------------------------------------------------------------
        
    end subroutine RecalculateVolume


    !--------------------------------------------------------------------------
    
    subroutine Verify()
        
        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
    
        !Local-----------------------------------------------------------------

        integer :: JW, JB, IU, ID, I, K

        !----------------------------------------------------------------------

cl10 :  DO JW=1,Me%CEQUALW2%NWB
cl11 :  DO JB=Me%CEQUALW2%BS(JW),Me%CEQUALW2%BE(JW)
        IU  = Me%CEQUALW2%US(JB)
        ID  = Me%CEQUALW2%DS(JB)
        !**** Branch and grid area and volume by layer
cl40 :  DO K=Me%CEQUALW2%KMX-1,3,-1
cl41 :  DO I=IU,ID
cf41 :  IF (K <= Me%CEQUALW2%KB(I)) THEN
cf42 :  if (Me%CEQUALW2%B(K,I) > Me%CEQUALW2%B(K-1,I))  then
            write(*,*)'Warning : Layer width decreases to the top at K=',K,' and I=',   &
                        I,' ModuleReservoirOptimization - Verify - ERR01'
        endif                                               cf42
        ENDIF                                               cf41
        ENDDO                                               cl41
        ENDDO                                               cl40
        ENDDO                                               cl11
        ENDDO                                               cl10

        !----------------------------------------------------------------------
        
    end subroutine Verify

    !------------------------------------------------------------------------

    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillReservoirOptimization(ObjReservoirOptimizationID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjReservoirOptimizationID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirOptimizationID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mRESERVOIROPTIMIZATION_,  Me%InstanceID)

            if (nUsers == 0) then
                deallocate(Me%CEQUALW2%CommentHeader)
                deallocate(Me%CEQUALW2%CommentDLX   )
                deallocate(Me%CEQUALW2%CommentELWS  )
                deallocate(Me%CEQUALW2%CommentPHI0  )
                deallocate(Me%CEQUALW2%CommentFRIC  )
                deallocate(Me%CEQUALW2%CommentH     )
                deallocate(Me%CEQUALW2%CommentB     )
                deallocate(Me%CEQUALW2%BTH          )
                deallocate(Me%CEQUALW2%BTHFN        )
                deallocate(Me%CEQUALW2%BS           )
                deallocate(Me%CEQUALW2%BE           )
                deallocate(Me%CEQUALW2%US           )
                deallocate(Me%CEQUALW2%DS           )
                deallocate(Me%CEQUALW2%DLX          )
                deallocate(Me%CEQUALW2%ELWS         )
                deallocate(Me%CEQUALW2%PHI0         )
                deallocate(Me%CEQUALW2%FRIC         )
                deallocate(Me%CEQUALW2%H            )
                deallocate(Me%CEQUALW2%EL           )
                deallocate(Me%CEQUALW2%B            )
                deallocate(Me%CEQUALW2%KTI          )
                deallocate(Me%CEQUALW2%KTWB         )
                deallocate(Me%CEQUALW2%KB           )
                deallocate(Me%CEQUALW2%ELBOT        )
                deallocate(Me%CEQUALW2%SLOPE        )
                deallocate(Me%CEQUALW2%ZERO_SLOPE   )
                deallocate(Me%CEQUALW2%VOLB         )
                deallocate(Me%CEQUALW2%VOLG         )
                deallocate(Me%CEQUALW2%DHS          )
                deallocate(Me%CEQUALW2%UHS          )
                deallocate(Me%CEQUALW2%DQB          )
                deallocate(Me%CEQUALW2%UQB          )
                deallocate(Me%CEQUALW2%JBDN         )
                deallocate(Me%CEQUALW2%NPOINT       )
                deallocate(Me%CEQUALW2%DAM_FLOW     )
                deallocate(Me%CEQUALW2%UH_EXTERNAL  )
                deallocate(Me%CEQUALW2%DH_EXTERNAL  )
                deallocate(Me%CEQUALW2%UH_INTERNAL  )
                deallocate(Me%CEQUALW2%DH_INTERNAL  )
                deallocate(Me%CEQUALW2%UQ_EXTERNAL  )
                deallocate(Me%CEQUALW2%DQ_EXTERNAL  )
                deallocate(Me%CEQUALW2%UQ_INTERNAL  )
                deallocate(Me%CEQUALW2%DQ_INTERNAL  )
                deallocate(Me%CEQUALW2%ALPHA        )
                deallocate(Me%CEQUALW2%Z            )
                deallocate(Me%CEQUALW2%NCCGR        )
                deallocate(Me%CEQUALW2%NCCBR        )
                deallocate(Me%CEQUALW2%SAGR         )
                deallocate(Me%CEQUALW2%SABR         )
                deallocate(Me%CEQUALW2%CVBR         )
                deallocate(Me%CEQUALW2%CVGR         )
                
                nullify (Me%CEQUALW2%CommentHeader  )
                nullify (Me%CEQUALW2%CommentDLX     )
                nullify (Me%CEQUALW2%CommentELWS    )
                nullify (Me%CEQUALW2%CommentPHI0    )
                nullify (Me%CEQUALW2%CommentFRIC    )
                nullify (Me%CEQUALW2%CommentH       )
                nullify (Me%CEQUALW2%CommentB       )
                nullify (Me%CEQUALW2%BTH            )
                nullify (Me%CEQUALW2%BTHFN          )
                nullify (Me%CEQUALW2%BS             )
                nullify (Me%CEQUALW2%BE             )
                nullify (Me%CEQUALW2%US             )
                nullify (Me%CEQUALW2%DS             )
                nullify (Me%CEQUALW2%DLX            )
                nullify (Me%CEQUALW2%ELWS           )
                nullify (Me%CEQUALW2%PHI0           )
                nullify (Me%CEQUALW2%FRIC           )
                nullify (Me%CEQUALW2%H              )
                nullify (Me%CEQUALW2%EL             )
                nullify (Me%CEQUALW2%B              )
                nullify (Me%CEQUALW2%KTI            )
                nullify (Me%CEQUALW2%KTWB           )
                nullify (Me%CEQUALW2%KB             )
                nullify (Me%CEQUALW2%ELBOT          )
                nullify (Me%CEQUALW2%SLOPE          )
                nullify (Me%CEQUALW2%ZERO_SLOPE     )
                nullify (Me%CEQUALW2%VOLB           )
                nullify (Me%CEQUALW2%VOLG           )
                nullify (Me%CEQUALW2%DHS            )
                nullify (Me%CEQUALW2%UHS            )
                nullify (Me%CEQUALW2%DQB            )
                nullify (Me%CEQUALW2%UQB            )
                nullify (Me%CEQUALW2%JBDN           )
                nullify (Me%CEQUALW2%NPOINT         )
                nullify (Me%CEQUALW2%DAM_FLOW       )
                nullify (Me%CEQUALW2%UH_EXTERNAL    )
                nullify (Me%CEQUALW2%DH_EXTERNAL    )
                nullify (Me%CEQUALW2%UH_INTERNAL    )
                nullify (Me%CEQUALW2%DH_INTERNAL    )
                nullify (Me%CEQUALW2%UQ_EXTERNAL    )
                nullify (Me%CEQUALW2%DQ_EXTERNAL    )
                nullify (Me%CEQUALW2%UQ_INTERNAL    )
                nullify (Me%CEQUALW2%DQ_INTERNAL    )
                nullify (Me%CEQUALW2%ALPHA          )
                nullify (Me%CEQUALW2%Z              )
                nullify (Me%CEQUALW2%NCCGR          )
                nullify (Me%CEQUALW2%NCCBR          )
                nullify (Me%CEQUALW2%SAGR           )
                nullify (Me%CEQUALW2%SABR           )
                nullify (Me%CEQUALW2%CVBR           )
                nullify (Me%CEQUALW2%CVGR           )

        
               
                deallocate(Me%FillData%Vol      )
                deallocate(Me%FillData%VolDiff  )
                deallocate(Me%FillData%VolError )

                nullify (Me%FillData%Vol        )
                nullify (Me%FillData%VolDiff    )
                nullify (Me%FillData%VolError   )
        
                !Deallocates Instance
                call DeallocateInstance ()

                ObjReservoirOptimizationID = 0
                STAT_      = SUCCESS_

            endif
        else 
            STAT_ = ready_
        endif cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillReservoirOptimization
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_ReservoirOptimization), pointer          :: AuxObjReservoirOptimization
        type (T_ReservoirOptimization), pointer          :: PreviousObjReservoirOptimization
       
        !Updates pointers
        if (Me%InstanceID == FirstObjReservoirOptimization%InstanceID) then
            FirstObjReservoirOptimization           => FirstObjReservoirOptimization%Next
        else
            PreviousObjReservoirOptimization        => FirstObjReservoirOptimization
            AuxObjReservoirOptimization             => FirstObjReservoirOptimization%Next
            do while (AuxObjReservoirOptimization%InstanceID /= Me%InstanceID)
                PreviousObjReservoirOptimization    => AuxObjReservoirOptimization
                AuxObjReservoirOptimization         => AuxObjReservoirOptimization%Next
            enddo

            !Now update linked list
            PreviousObjReservoirOptimization%Next   => AuxObjReservoirOptimization%Next

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

    subroutine Ready (ObjReservoirOptimization_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjReservoirOptimization_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjReservoirOptimization_ID > 0) then
            call LocateObjReservoirOptimization (ObjReservoirOptimization_ID)
            ready_ = VerifyReadLock (mRESERVOIROPTIMIZATION_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        endif cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjReservoirOptimization (ObjReservoirOptimizationID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjReservoirOptimizationID

        !Local-----------------------------------------------------------------

        Me => FirstObjReservoirOptimization
        do while (associated (Me))
            if (Me%InstanceID == ObjReservoirOptimizationID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleReservoirOptimization - LocateObjReservoirOptimization - ERR01'

    end subroutine LocateObjReservoirOptimization

    !--------------------------------------------------------------------------

end module ModuleReservoirOptimization

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------








