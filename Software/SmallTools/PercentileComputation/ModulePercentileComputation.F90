!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : PercentileComputation
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to serve as PercentileComputation to create new modules
!
!------------------------------------------------------------------------------


Module ModulePercentileComputation

    use ModuleGlobalData
    use ModuleTime               
    use ModuleEnterData,         only : ConstructEnterData, KillEnterData,              &
                                        GetData, ExtractBlockFromBuffer, Block_Unlock
    use ModuleFunctions         
    use ModuleHDF5
    
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructPercentileComputation
    private ::      AllocateInstance
    
    !Modifier
    public  :: ModifyPercentileComputation

    !Destructor
    public  :: KillPercentileComputation                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjPercentileComputation 

    !Types---------------------------------------------------------------------
    
    private :: T_Conditions
    type       T_Conditions
        character(len=StringLength)                 :: PropName
        logical                                     :: Lower
        real                                        :: Limit
        real                                        :: Percentile
    end type  T_Conditions    
    
    ! Definition of type T_Parameter
    type       T_Parameter
        character(len=StringLength)                 :: DataSetName
        character(len=PathLength  )                 :: GroupName
        character(len=PathLength  )                 :: File        
        integer                                     :: ObjHDF5
        type(T_Parameter), pointer                  :: Next     => null()
    end type  T_Parameter
    
    
    private :: T_PercentileComputation
    type       T_PercentileComputation
        integer                                     :: InstanceID
        integer                                     :: ObjEnterData     = 0
        integer                                     :: IUB, JUB, KUB
        real                                        :: DX              
        real                                        :: Xorig, Yorig
        real                                        :: FillValueIn
        real                                        :: FillValueOut = -99
        character(Len=PathLength)                   :: OutputESRI
        real, dimension(:, :),     pointer          :: OutMatrix2D
        type (T_Parameter), pointer                 :: FirstParameter     
        integer                                     :: ParameterNumber        
        type (T_Conditions),  dimension(:), pointer :: Conditions 
        integer                                     :: NumberCond        
        type(T_PercentileComputation), pointer      :: Next
    end type  T_PercentileComputation

    !Global Module Variables
    type (T_PercentileComputation), pointer         :: FirstObjPercentileComputation
    type (T_PercentileComputation), pointer         :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructPercentileComputation(ObjPercentileComputationID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjPercentileComputationID 
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mGeneric_)) then
            nullify (FirstObjPercentileComputation)
            call RegisterModule (mGeneric_) 
        endif

        call Ready(ObjPercentileComputationID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            call ReadInputFile
            
            call ReadKeywords

            !Returns ID
            ObjPercentileComputationID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModulePercentileComputation - ConstructPercentileComputation - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructPercentileComputation
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_PercentileComputation), pointer                         :: NewObjPercentileComputation
        type (T_PercentileComputation), pointer                         :: PreviousObjPercentileComputation


        !Allocates new instance
        allocate (NewObjPercentileComputation)
        nullify  (NewObjPercentileComputation%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjPercentileComputation)) then
            FirstObjPercentileComputation         => NewObjPercentileComputation
            Me                    => NewObjPercentileComputation
        else
            PreviousObjPercentileComputation      => FirstObjPercentileComputation
            Me                    => FirstObjPercentileComputation%Next
            do while (associated(Me))
                PreviousObjPercentileComputation  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjPercentileComputation
            PreviousObjPercentileComputation%Next => NewObjPercentileComputation
        endif

        Me%InstanceID = RegisterNewInstance (mGeneric_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------
    subroutine ReadInputFile
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: iFile, l, is, ie, AuxInt
        integer                                         :: STAT_CALL
        character(Len=1000)                             :: AuxString

        !----------------------------------------------------------------------
        call UnitsManager(iFile, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadInputFile - ModulePercentileComputation - ERR10'

        open(Unit   = iFile,                                                            &
             File   = 'InputFile.dat',                                                  &
             Form   = 'FORMATTED',                                                      &
             STATUS = 'UNKNOWN',                                                        &
             Action = 'READ',                                                           &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadInputFile - ModulePercentileComputation - ERR20'
                            
        read(iFile,*) Me%NumberCond
        
        allocate(Me%Conditions(1:Me%NumberCond))
        
        do l=1, Me%NumberCond
        
            read(iFile,'(A1000)') AuxString
        
            is = 1
            ie = index(AuxString(is:1000), ',')-1
            Me%Conditions(l)%PropName = trim(adjustl(AuxString(is:ie)))
            
            is = ie+2
            ie = is + index(AuxString(is:1000), ',')-1
            read(AuxString(is:ie),*) AuxInt
            if      (AuxInt == 0) then
                Me%Conditions(l)%Lower = .false.
            elseif  (AuxInt == 1) then
                Me%Conditions(l)%Lower = .true. 
            else
                stop 'ReadInputFile - ModulePercentileComputation - ERR30' 
            endif
            
            is = ie+2
            ie = is + index(AuxString(is:1000), ',')-1
            read(AuxString(is:ie),*) Me%Conditions(l)%Limit

            is = ie+2
            read(AuxString(is:1000),*) Me%Conditions(l)%Percentile            
            
        enddo            

        call UnitsManager(iFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadInputFile - ModulePercentileComputation - ERR40'
        
    end subroutine ReadInputFile
    
    !--------------------------------------------------------------------------        
    
    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call ConstructEnterData (Me%ObjEnterData, "HDFInputFiles.dat", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadKeywords - ModulePercentileComputation - ERR10'

        call ReadGlobalData        

        call ReadParameters

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadKeywords - ModulePercentileComputation - ERR20'

    end subroutine ReadKeywords

    !--------------------------------------------------------------------------

    subroutine ReadGlobalData

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------

        ! DX
        call GetData(Me%DX, Me%ObjEnterData, iflag,                                     &
                     keyword      = 'DX',                                               &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
            stop 'ReadGlobalData - ModulePercentileComputation - ERR10'   

        ! X origin
        call GetData(Me%Xorig, Me%ObjEnterData, iflag,                                  &
                     keyword      = 'XORIG',                                            &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
            stop 'ReadGlobalData - ModulePercentileComputation - ERR20'   

        ! Y origin
        call GetData(Me%Yorig, Me%ObjEnterData, iflag,                                  &
                     keyword      = 'YORIG',                                            &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
            stop 'ReadGlobalData - ModulePercentileComputation - ERR30'          
        
        ! IUB
        call GetData(Me%IUB, Me%ObjEnterData, iflag,                                    &
                     keyword      = 'IUB',                                              &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
            stop 'ReadGlobalData - ModulePercentileComputation - ERR40'          
        
        ! JUB
        call GetData(Me%JUB, Me%ObjEnterData, iflag,                                    &
                     keyword      = 'JUB',                                              &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
            stop 'ReadGlobalData - ModulePercentileComputation - ERR50'
        
        ! KUB
        call GetData(Me%KUB, Me%ObjEnterData, iflag,                                    &
                     keyword      = 'KUB',                                              &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
            stop 'ReadGlobalData - ModulePercentileComputation - ERR55'        
        
        !call GetData(Me%FillValueOut, Me%ObjEnterData, iflag,                           &
        !             keyword      = 'FILLVALUE_OUT',                                    &
        !             SearchType   = FromFile,                                           &
        !             ClientModule = 'ModulePercentileComputation',                      &
        !             STAT         = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
        !    stop 'ReadGlobalData - ModulePercentileComputation - ERR60'

        call GetData(Me%OutputESRI, Me%ObjEnterData, iflag,                             &
                     keyword      = 'OUTPUT_ESRI',                                      &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
            stop 'ReadGlobalData - ModulePercentileComputation - ERR70'
        
        call GetData(Me%FillValueIn, Me%ObjEnterData, iflag,                            &
                     keyword      = 'FILLVALUE_IN',                                     &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
            stop 'ReadGlobalData - ModulePercentileComputation - ERR80'
                
    end subroutine ReadGlobalData

    !--------------------------------------------------------------------------    

    subroutine ReadParameters

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, ClientNumber, i
        type (T_Parameter),       pointer               :: NewParameter
        logical                                         :: BlockFound
        logical                                         :: AtLeastOneBlock = .false.

        !Begin-----------------------------------------------------------------
        

        ! Obtain Parameters for the statistics' calculation
do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                        ClientNumber    = ClientNumber,                 &
                                        block_begin     = '<BeginParameter>',           &
                                        block_end       = '<EndParameter>',             &
                                        BlockFound      = BlockFound,                   &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  

                    
                    call AddParameter           (NewParameter)

                    call ConstructParameters    (NewParameter)

                    nullify                     (NewParameter)

                else cd2
                    call Block_Unlock(Me%ObjEnterData,                                  & 
                                      ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ReadParameters - ModulePercentileComputation - ERR10'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ReadParameters - ModulePercentileComputation - ERR20'
            else cd1
                stop 'ReadParameters - ModulePercentileComputation - ERR30'
            end if cd1

        end do do1


    end subroutine ReadParameters

    
    !--------------------------------------------------------------------------

    subroutine AddParameter (ObjParameter)

        !Arguments-------------------------------------------------------------
        type (T_Parameter),     pointer           :: ObjParameter

        !Local-----------------------------------------------------------------
        type (T_Parameter),     pointer           :: PreviousParameter
        type (T_Parameter),     pointer           :: NewParameter

        !Begin-----------------------------------------------------------------

        !Allocates new Parameter
        allocate (NewParameter)
        nullify  (NewParameter%Next)

        !Insert new Parameter into list and makes current ?? point to it
        if (.not. associated(Me%FirstParameter)) then
            Me%FirstParameter         => NewParameter
            ObjParameter              => NewParameter
            Me%ParameterNumber = 1
        else
            PreviousParameter         => Me%FirstParameter
            ObjParameter              => Me%FirstParameter%Next
            do while (associated(ObjParameter))
                PreviousParameter     => ObjParameter
                ObjParameter          => ObjParameter%Next
            enddo
            ObjParameter              => NewParameter
            PreviousParameter%Next    => NewParameter
            ! Count number of parameters in list
            Me%ParameterNumber = Me%ParameterNumber + 1
        end if

    end subroutine AddParameter

    !--------------------------------------------------------------------------

    subroutine ConstructParameters (NewParameter)

        !Arguments-------------------------------------------------------------
        type (T_Parameter),      pointer          :: NewParameter

        !Local-----------------------------------------------------------------
        integer                                   :: iflag, STAT_CALL, HDF5_READ
        
        !Begin-----------------------------------------------------------------
        
        ! Obtain parameter name
        call GetData(NewParameter%DataSetName,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'PROPERTY',                                         &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'ConstructParameters - ModulePercentileComputation - ERR10'

        if (.not.CheckPropertyName(NewParameter%DataSetName)) then
            write(*,*)
            write(*,*) 'The property name is not recognised by the model:'
            write(*,*) trim(NewParameter%DataSetName)
            !stop 'ConstructParameters - ModuleHDF5Statistics - ERR02'
            !No stop is made because in some HDF5 files parameters are not
            !registred in Module GlobalData
        end if

        ! Obtain parameter group
        call GetData(NewParameter%GroupName,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'HDF_GROUP',                                        &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_ .and. iflag == 0)                                   &
            stop 'ConstructParameters - ModulePercentileComputation - ERR20'
        
        call GetData(NewParameter%File,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FILENAME',                                         &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_ .and. iflag == 0)                                   &
            stop 'ConstructParameters - ModulePercentileComputation - ERR30'                     
        
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
            
        NewParameter%ObjHDF5 = 0

        !Open HDF5 file
        call ConstructHDF5 (NewParameter%ObjHDF5, trim(NewParameter%File),              &
                            HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'ConstructParameters - ModulePercentileComputation - ERR40'                     

    end subroutine ConstructParameters

    !--------------------------------------------------------------------------

    subroutine Search_Parameter(ParameterX, ParameterName, STAT)

        !Arguments-------------------------------------------------------------
        type (T_Parameter),      pointer                :: ParameterX
        character(Len=*),           intent (IN)         :: ParameterName
        integer         , optional, intent (OUT)        :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_

        !----------------------------------------------------------------------

        STAT_  = UNKNOWN_

        ParameterX => Me%FirstParameter

        do while (associated(ParameterX))
            if (trim(adjustl(ParameterX%DataSetName))==trim(adjustl(ParameterName))) then
                exit
            else
                ParameterX => ParameterX%Next
            end if
        end do

       if (associated(ParameterX)) then

            STAT_ = SUCCESS_

        else
            STAT_  = NOT_FOUND_ERR_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine Search_Parameter

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyPercentileComputation(ObjPercentileComputationID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPercentileComputationID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPercentileComputationID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
            
            !Read properties blocks
            call ComputePercentile

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyPercentileComputation
    
    !--------------------------------------------------------------------------
    
    subroutine ComputePercentile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Parameter), pointer :: ParameterX
        integer                     :: n, STAT_CALL, IUB, JUB    
        
        !Begin-----------------------------------------------------------------
        
        
        IUB = Me%IUB
        JUB = Me%JUB
        
        allocate(Me%OutMatrix2D(1:IUB,1:JUB))
        
        Me%OutMatrix2D(1:IUB,1:JUB) = 0.
    
        do n = 1, Me%NumberCond
            
            nullify(ParameterX)
            
            call Search_Parameter(ParameterX, Me%Conditions(n)%PropName, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) then
                write(*,*) 'No hdf5 file found with the parameter name = ',trim(Me%Conditions(n)%PropName)
                stop 'ComputePercentile - ModulePercentileComputation - ERR10'                     
            endif
            
            call ReadFileCheckCondition(ParameterX, n)
            
        enddo
        
        call WriteESRI_GridData
        
        deallocate(Me%OutMatrix2D)
        
    end subroutine ComputePercentile
    
    !--------------------------------------------------------------------------
    
    subroutine ReadFileCheckCondition(ParameterX, n)

        !Arguments-------------------------------------------------------------
        type (T_Parameter), pointer         :: ParameterX
        integer                             :: n
        
        !Local-----------------------------------------------------------------
        real,   dimension(:,:,:), pointer   :: Classes3D, ReadMatrix3D
        integer,dimension(:,:,:), pointer   :: WaterPoints3D
        real,   dimension(:,:  ), pointer   :: Limits2D        
        character(len=StringLength)         :: FieldName, LimitsName        
        character(len=PathLength  )         :: GroupName
        integer                             :: STAT_CALL, IUB, JUB, KUB, ObjHDF5
        integer                             :: nItems, ni  
        
        !Begin-----------------------------------------------------------------
        
        
        IUB = Me%IUB
        JUB = Me%JUB
        KUB = Me%KUB
        
        ObjHDF5     = ParameterX%ObjHDF5
        GroupName   = ParameterX%GroupName
        FieldName   = "Classes"
        LimitsName  = "Classes_Limits"
        
        call GetHDF5GroupNumberOfItems(HDF5ID    = ObjHDF5,                             &
                                       GroupName = trim(GroupName)//'/'//trim(FieldName),&
                                       nItems    = nItems,                              &
                                       STAT      = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR10' 
                                                
        
        allocate(Classes3D   (0:IUB+1,0:JUB+1,0:nItems+1))
        allocate(ReadMatrix3D(0:IUB+1,0:JUB+1,0:KUB+1   ))  
        allocate(Limits2D    (1:2    ,1:nItems          ))
        
        Classes3D   (1:IUB,1:JUB,1:nItems) = 0.
        ReadMatrix3D(1:IUB,1:JUB,1:KUB   ) = 0.
                            
        call HDF5SetLimits  (HDF5ID = ObjHDF5,                                          &
                             ILB = 1, IUB = 2, JLB = 1,JUB = nItems,                    &
                             STAT   = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR30'

        call HDF5ReadData   (HDF5ID        = ObjHDF5,                                   &
                             GroupName     = trim(GroupName),                           &
                             Name          = trim(LimitsName),                          &
                             Array2D       = Limits2D,                                  &
                             STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR40'
        
    
        do ni = 1, nItems
            
            call HDF5SetLimits  (HDF5ID = ObjHDF5, ILB = 1, IUB = IUB,                  &
                                                   JLB = 1, JUB = JUB,                  &
                                                   KLB = 1, KUB = KUB,                  &  
                                 STAT   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR50'

            call HDF5ReadWindow(HDF5ID        = ObjHDF5,                                &
                                GroupName     = trim(GroupName)//'/'//trim(FieldName),  &
                                Name          = trim(FieldName),                        &
                                Array3D       = ReadMatrix3D,                           &
                                OutputNumber  = ni,                                     &  
                                STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR60'

            Classes3D(:,:,ni) = ReadMatrix3D(:,:,KUB)
            
        enddo
        
        allocate(WaterPoints3D(0:IUB+1,0:JUB+1,0:KUB+1   ))         
        
        call HDF5ReadWindow(HDF5ID        = ObjHDF5,                                    &
                            GroupName     = "/Grid",                                    &
                            Name          = "WaterPoints3D",                            &
                            Array3D       = WaterPoints3D,                              &
                            STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR60'        
        
        call CheckCondition (Classes3D, Limits2D, Me%OutMatrix2D, nItems, n, WaterPoints3D)        

        deallocate(ReadMatrix3D )
        deallocate(Classes3D    ) 
        deallocate(WaterPoints3D)

    end subroutine ReadFileCheckCondition
    
    !--------------------------------------------------------------------------
    
    subroutine CheckCondition(InMatrix3D, Limits2D, OutMatrix2D, nItems, n, WaterPoints3D)

        !Arguments-------------------------------------------------------------
        real,   dimension(:,:,:), pointer   :: InMatrix3D
        integer,dimension(:,:,:), pointer   :: WaterPoints3D
        real,   dimension(:,:  ), pointer   :: OutMatrix2D, Limits2D
        integer                             :: nItems, n
        
        !Local-----------------------------------------------------------------
        real                                :: Sum, dP2, dPx, LimitX
        integer                             :: STAT_CALL, IUB, JUB, KUB, ni, i, j
        
        !Begin-----------------------------------------------------------------
        
        
        IUB = Me%IUB
        JUB = Me%JUB
        KUB = Me%KUB
        
        
        do i = 1, IUB
        do j = 1, JUB
            
            if (WaterPoints3D(i,j,KUB) == 0) then
                OutMatrix2D(i, j) = Me%FillValueOut
                cycle
            endif
            
            Sum = 0.
            
            do ni = 1, nItems
                
                Sum = Sum + InMatrix3D(i, j, ni)
                
                if (Sum >= Me%Conditions(n)%Percentile) then 
                    
                    dPx = Sum - Me%Conditions(n)%Percentile
                    dP2 = InMatrix3D(i, j, ni)                    
                    
                    if (dPx == 0. .or. ni == nItems .or. dP2 == 0.) then
                        LimitX = Limits2D(2,ni)
                    else
                        LimitX = (Limits2D(2,ni) - Limits2D(1,ni)) * dPx/dP2 + Limits2D(1,ni)
                    endif
                
                    if (Me%Conditions(n)%Lower) then
                        if (Me%Conditions(n)%Limit >= LimitX) then 
                            OutMatrix2D(i, j) = OutMatrix2D(i, j) + 1
                        endif
                    endif
                    
                    if (.not. Me%Conditions(n)%Lower) then
                        if (Me%Conditions(n)%Limit <= LimitX) then 
                            OutMatrix2D(i, j) = OutMatrix2D(i, j) + 1
                        endif
                    endif                    
                    
                    exit    
                
                endif
                
            enddo
        enddo
        enddo


    end subroutine CheckCondition
    
    !--------------------------------------------------------------------------  
    
    subroutine WriteESRI_GridData
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len=50000)        :: Line
        integer                     :: Unit, STAT_CALL, i, j    
        integer                     :: a, ba, bt
        logical                     :: Found2Blanks        
        
        !Begin-----------------------------------------------------------------    
    
        call UnitsManager(Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WriteESRI_GridData - ModulePercentileComputation - ERR10'

        open(Unit   = Unit,                                                             &
             File   = trim(Me%OutputESRI),                                              &
             Form   = 'FORMATTED',                                                      &
             STATUS = 'UNKNOWN',                                                        &
             Action = 'WRITE',                                                          &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WriteESRI_GridData - ModulePercentileComputation - ERR20'
                

        write(Unit,'(A14,I6)'   ) 'ncols         ', Me%JUB
        write(Unit,'(A14,I6)'   ) 'nrows         ', Me%IUB
        write(Unit,'(A14,f12.6)') 'xllcorner     ', Me%XOrig
        write(Unit,'(A14,f12.6)') 'yllcorner     ', Me%YOrig
        write(Unit,'(A14,f12.6)') 'cellsize      ', Me%DX
        write(Unit,'(A14,f12.6)') 'nodata_value  ', Me%FillValueOut
                
        do i = Me%IUB, 1, -1
            do j = 1, Me%JUB
                if (Me%OutMatrix2D(i,j) <= Me%FillValueOut) Me%OutMatrix2D(i,j) = Me%FillValueOut
            enddo
            write(Line,'(4000(I4,1x))') int(Me%OutMatrix2D(i,1:Me%JUB))
            Line = adjustl(Line)
            Found2Blanks = .true.

            bt = len_Trim(Line)
            ba = bt
            do a = 1, bt-1
                do while (Line(a:a+1)=='  ') 
                    Line(a:ba)=Line(a+1:ba+1)
                    ba = ba - 1
                    if (ba == a) exit
                enddo
                if (ba == a) exit
            enddo
            write(Unit,'(A)') trim(Line)
        enddo
                
        call UnitsManager(Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WriteESRI_GridData - ModulePercentileComputation - ERR30'   
        
    end subroutine WriteESRI_GridData        
    
    !--------------------------------------------------------------------------        

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillPercentileComputation(ObjPercentileComputationID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjPercentileComputationID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPercentileComputationID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mGeneric_,  Me%InstanceID)

            if (nUsers == 0) then

                !Deallocates Instance
                call DeallocateInstance ()

                ObjPercentileComputationID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillPercentileComputation
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_PercentileComputation), pointer          :: AuxObjPercentileComputation
        type (T_PercentileComputation), pointer          :: PreviousObjPercentileComputation

        !Updates pointers
        if (Me%InstanceID == FirstObjPercentileComputation%InstanceID) then
            FirstObjPercentileComputation => FirstObjPercentileComputation%Next
        else
            PreviousObjPercentileComputation => FirstObjPercentileComputation
            AuxObjPercentileComputation      => FirstObjPercentileComputation%Next
            do while (AuxObjPercentileComputation%InstanceID /= Me%InstanceID)
                PreviousObjPercentileComputation => AuxObjPercentileComputation
                AuxObjPercentileComputation      => AuxObjPercentileComputation%Next
            enddo

            !Now update linked list
            PreviousObjPercentileComputation%Next => AuxObjPercentileComputation%Next

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

    subroutine Ready (ObjPercentileComputation_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPercentileComputation_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjPercentileComputation_ID > 0) then
            call LocateObjPercentileComputation (ObjPercentileComputation_ID)
            ready_ = VerifyReadLock (mGeneric_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjPercentileComputation (ObjPercentileComputationID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPercentileComputationID

        !Local-----------------------------------------------------------------

        Me => FirstObjPercentileComputation
        do while (associated (Me))
            if (Me%InstanceID == ObjPercentileComputationID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModulePercentileComputation - LocateObjPercentileComputation - ERR01'

    end subroutine LocateObjPercentileComputation

    !--------------------------------------------------------------------------

end module ModulePercentileComputation








