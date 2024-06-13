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
    use ModuleEnterData
    use ModuleFunctions         
    use ModuleHDF5
    use ModuleDrawing
    use RasterFunctions !new
    !use ModuleHorizontalGrid
    !use ModuleGridData
    
    
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
    
    private :: T_PolyField
    type       T_PolyField
        character(len=StringLength)                 :: PropName
        integer                                     :: LayerON
    end type  T_PolyField        
    
    ! Definition of type T_Parameter
    type       T_Parameter
        character(len=StringLength)                 :: DataSetName
        character(len=PathLength  )                 :: GroupName
        character(len=PathLength  )                 :: File        
        integer                                     :: ObjHDF5
        character(Len=StringLength)                 :: MappingName  
        logical                                     :: PolyField     = .false. 
        type (T_Polygon), pointer                   :: ExclusionAreas            
        type(T_Parameter), pointer                  :: Next         => null()
    end type  T_Parameter
    
    
    !private :: T_PercentileComputation
    type       T_PercentileComputation
        integer                                     :: InstanceID
        integer                                     :: ObjEnterData     = 0
        integer                                     :: ObjGridData      = 0
        integer                                     :: ObjGrid          = 0
        integer                                     :: IUB, JUB, KUB
        real                                        :: DX, DY           
        real                                        :: Xorig, Yorig
        real                                        :: FillValueIn
        real                                        :: FillValueOut = -99
        character(Len=PathLength)                   :: OutputESRI
        character(Len=PathLength)                   :: OutputGeoTiff
        character(Len=PathLength)                   :: BathymFile        
        real, dimension(:, :, :),  pointer          :: OutMatrix3D
        real, dimension(:, :),     pointer          :: Bathym2D
        real                                        :: BathymMin, BathymMax
        type (T_Parameter), pointer                 :: FirstParameter     
        integer                                     :: ParameterNumber        
        type (T_Conditions),  dimension(:), pointer :: Conditions 
        type (T_PolyField ),  dimension(:), pointer :: PolyField        
        integer                                     :: NumberCond
        integer                                     :: NumberPoly
        logical                                     :: WriteAsGeoTiff = .false.
        logical                                     :: WriteAsAsc     = .false.
        logical                                     :: WriteCorners   = .false.
        integer                                     :: EPSG
        type(T_PercentileComputation), pointer      :: Next
    end type  T_PercentileComputation

    !Global Module Variables
    type (T_PercentileComputation), pointer         :: FirstObjPercentileComputation
    type (T_PercentileComputation), pointer         :: Me


    !--------------------------------------------------------------------------
    !allocate(Me) 
    
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

            call ReadBathymFile
            
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
        integer                                         :: STAT_CALL, ios
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
                            
        read(iFile,*) Me%BathymMin, Me%BathymMax
                            
        read(iFile,*) Me%NumberCond
        
        allocate(Me%Conditions(1:Me%NumberCond))
        
        do l=1, Me%NumberCond
        
            read(iFile,'(A1000)') AuxString
        
            is = 1
            ie = index(AuxString(is:1000), ',')-1
            Me%Conditions(l)%PropName = trim(adjustl(AuxString(is:ie)))
            
            is = ie+2
            ie = is + index(AuxString(is:1000), ',')-2
            read(AuxString(is:ie),*) AuxInt
            if      (AuxInt == 0) then
                Me%Conditions(l)%Lower = .false.
            elseif  (AuxInt == 1) then
                Me%Conditions(l)%Lower = .true. 
            else
                stop 'ReadInputFile - ModulePercentileComputation - ERR30' 
            endif
            
            is = ie+2
            ie = is + index(AuxString(is:1000), ',')-2
            read(AuxString(is:ie),*) Me%Conditions(l)%Limit

            is = ie+2
            read(AuxString(is:1000),*) Me%Conditions(l)%Percentile            
            
        enddo            
        
        read(iFile,*,iostat=ios) Me%NumberPoly 
        if (ios == 0) then
        
            allocate(Me%PolyField(1:Me%NumberPoly))        
        
            do l=1, Me%NumberPoly
        
                read(iFile,'(A1000)') AuxString
        
                is = 1
                ie = index(AuxString(is:1000), ',')-1
                Me%PolyField(l)%PropName = trim(adjustl(AuxString(is:ie)))
            
                is = ie+2
                ie = 1000
                read(AuxString(is:ie),*) Me%PolyField(l)%LayerON
                if ( Me%PolyField(l)%LayerON /= 0 .and.  Me%PolyField(l)%LayerON /= 1) then
                    stop 'ReadInputFile - ModulePercentileComputation - ERR40'
                endif

            enddo
        
        else

            Me%NumberPoly = 0        
            
        endif
        
        call UnitsManager(iFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadInputFile - ModulePercentileComputation - ERR50'
        
        
        
    end subroutine ReadInputFile
    
    !--------------------------------------------------------------------------        
    
    subroutine ReadBathymFile
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: iFile, i, j
        integer                                         :: STAT_CALL
        !character(Len=1000)                             :: AuxString

        !----------------------------------------------------------------------
        call UnitsManager(iFile, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadBathymFile - ModulePercentileComputation - ERR10'
        
        
        allocate(Me%Bathym2D(1:Me%IUB,1:Me%JUB))

        open(Unit   = iFile,                                                            &
             File   = Me%BathymFile,                                                    &
             Form   = 'FORMATTED',                                                      &
             STATUS = 'UNKNOWN',                                                        &
             Action = 'READ',                                                           &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadBathymFile - ModulePercentileComputation - ERR20'
        
        do i=1, Me%IUB
        do j=1, Me%JUB            
            read(iFile,*) Me%Bathym2D(i, j)
        enddo
        enddo
            

        call UnitsManager(iFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadBathymFile - ModulePercentileComputation - ERR40'
        
    end subroutine ReadBathymFile
    
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
        
        ! DY
        call GetData(Me%DY, Me%ObjEnterData, iflag,                                     &
                     keyword      = 'DY',                                               &
                     SearchType   = FromFile,                                           &
                     Default      = Me%DX,                                              &                       
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadGlobalData - ModulePercentileComputation - ERR15'           

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
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalData - ModulePercentileComputation - ERR70'
        if (iflag /= 0) Me%WriteAsAsc = .true.
                
        call GetData(Me%OutputGeoTiff, Me%ObjEnterData, iflag,                          &
                     keyword      = 'OUTPUT_TIFF',                                      &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalData - ModulePercentileComputation - ERR71'
        if (iflag /= 0) Me%WriteAsGeoTiff = .true.
        
        if (Me%WriteAsGeoTiff) then
            call GetData(Me%EPSG,                                                           &
                         Me%ObjEnterData, iflag,                                            &
                         SearchType   = FromFile,                                           &
                         keyword      = 'EPSG',                                             &
                         ClientModule = 'HDF5ToASCIIandBIN',                                &
                         default      = 4326,                                               &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePercentileComputation - ERR71.1'
        end if
        
        if (.not. Me%WriteAsGeoTiff .and. .not. Me%WriteAsAsc) then
            write(*,*) 'No output file provided, ASC or GEOTIFF'
            stop 'ReadGlobalData - ModulePercentileComputation - ERR72'
        endif
        
        if (Me%WriteAsGeoTiff) then
            call GetData(Me%WriteCorners, Me%ObjEnterData, iflag,                       &
                     keyword      = 'WRITE_CORNERS',                                    &
                     Default      = .false.,                                            &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePercentileComputation - ERR73'
        end if
        
        
        call GetData(Me%FillValueIn, Me%ObjEnterData, iflag,                            &
                     keyword      = 'FILLVALUE_IN',                                     &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
            stop 'ReadGlobalData - ModulePercentileComputation - ERR80'
                
        call GetData(Me%BathymFile, Me%ObjEnterData, iflag,                             &
                     keyword      = 'BATHYM_FILE',                                      &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
            stop 'ReadGlobalData - ModulePercentileComputation - ERR90'
        
        
        !call ConstructHorizontalGrid(HorizontalGridID = Me%ObjGrid,                     &
        !                             DataFile         = Me%BathymFile,                  &
        !                             STAT             = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)                                                      &
        !    stop 'ReadGlobalData - ModulePercentileComputation - ERR100'        
        !
        !call ConstructGridData(GridDataID       = Me%ObjGridData,                       &
        !                       HorizontalGridID = Me%ObjGrid,                           &
        !                       FileName         = Me%BathymFile,                        &
        !                       STAT             = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)                                                      &
        !    stop 'ReadGlobalData - ModulePercentileComputation - ERR110'
        !
        !call GetGridData(GridDataID     = Me%ObjGridData,                               &
        !                 GridData2D     = Me%Bathym2D,                                  &
        !                 STAT           = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)                                                      &
        !    stop 'ReadGlobalData - ModulePercentileComputation - ERR120'
        !
        
                
    end subroutine ReadGlobalData

    !--------------------------------------------------------------------------    

    subroutine ReadParameters

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, ClientNumber!, i
        type (T_Parameter),       pointer               :: NewParameter
        logical                                         :: BlockFound
        !logical                                         :: AtLeastOneBlock = .false.

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
        
        call GetData(NewParameter%File,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FILENAME',                                         &
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_ .and. iflag == 0)                                   &
            stop 'ConstructParameters - ModulePercentileComputation - ERR30'         
        
        ! Obtain parameter type
        call GetData(NewParameter%PolyField,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'POLY_FILE',                                        &
                     default      = .false.,                                            &    
                     ClientModule = 'ModulePercentileComputation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) then
            stop 'ConstructParameters - ModulePercentileComputation - ERR15'        
        endif
        
        if (.not.NewParameter%PolyField) then
            ! Obtain parameter group
            call GetData(NewParameter%GroupName,                                            &
                         Me%ObjEnterData, iflag,                                            &
                         SearchType   = FromBlock,                                          &
                         keyword      = 'HDF_GROUP',                                        &
                         ClientModule = 'ModulePercentileComputation',                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ .and. iflag == 0)                                   &
                stop 'ConstructParameters - ModulePercentileComputation - ERR20'
        
        
            call GetData(NewParameter%MappingName, Me%ObjEnterData, iflag,                  &
                         keyword      = 'MAPPING',                                          &
                         SearchType   = FromBlock,                                          &
                         ClientModule = 'ModulePercentileComputation',                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
                stop 'ConstructParameters - ModulePercentileComputation - ERR40'   

        
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
            
            NewParameter%ObjHDF5 = 0

            !Open HDF5 file
            call ConstructHDF5 (NewParameter%ObjHDF5, trim(NewParameter%File),              &
                                HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                    &
                stop 'ConstructParameters - ModulePercentileComputation - ERR50'                     
        endif

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
        
        allocate(Me%OutMatrix3D(1:IUB,1:JUB,0:Me%NumberCond))
        
        Me%OutMatrix3D(1:IUB,1:JUB,:) = 0
        if (Me%WriteAsGeoTiff) then
            call ReadHD5Grid        
        endif
        
    
        do n = 1, Me%NumberCond
            
            nullify(ParameterX)
            
            call Search_Parameter(ParameterX, Me%Conditions(n)%PropName, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) then
                write(*,*) 'No hdf5 file found with the parameter name = ',trim(Me%Conditions(n)%PropName)
                stop 'ComputePercentile - ModulePercentileComputation - ERR10'                     
            endif
            
            if (.not. ParameterX%PolyField) then
                call ReadFileCheckCondition(ParameterX, n)
            endif
            
        enddo
        
        call BathymFilter
        
        do n = 1, Me%NumberPoly
            
            if ( Me%PolyField(n)%LayerON == 1) then
            
                nullify(ParameterX)
            
                call Search_Parameter(ParameterX, Me%PolyField(n)%PropName, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) then
                    write(*,*) 'No hdf5 file found with the parameter name = ',trim(Me%PolyField(n)%PropName)
                    stop 'ComputePercentile - ModulePercentileComputation - ERR20'                     
                endif
            
                if (ParameterX%PolyField) then
                    call CheckExclusionAreas(ParameterX)
                endif
                
            endif
            
        enddo        
        
        if (Me%WriteAsAsc) then
            call WriteESRI_GridData
        end if
        
        if (Me%WriteAsGeoTiff) then
            call WriteESRI_GeoTiff
        end if
        
        deallocate(Me%OutMatrix3D)
        
    end subroutine ComputePercentile
    
    !--------------------------------------------------------------------------
    
    subroutine ReadHD5Grid
    
        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                             :: i, j, STAT_CALL
        integer                             :: ObjHDF5
        real, dimension(:,:  ), pointer     :: CoordsX, CoordsY
        real, dimension(:,:,:), allocatable :: Aux3D
        type (T_Parameter)    , pointer     :: ParameterX
        real(kind=c_double)                 :: geotrans(6) = 0.0
        character(len=1000)                 :: projref
        integer                             :: Unit ! to write 4 corners
        !Body------------------------------------------------------------------
        
        ! Get one HDF (the first) to check for conditions
        nullify(ParameterX)
        call Search_Parameter(ParameterX, Me%Conditions(1)%PropName, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ComputePercentile - ModulePercentileComputation - ERR10'
        
        ObjHDF5     = ParameterX%ObjHDF5
        
        ! Is and Js from HDF instead of 
        
        ! Read Coordinates from HDF5
        allocate(CoordsX(0:Me%IUB+2,0:Me%JUB+2)); CoordsX = 0.0
        allocate(CoordsY(0:Me%IUB+2,0:Me%JUB+2)); CoordsY = 0.0
        
        ! CoordX and CoordY have one extra size because its the corners of the HDF cells
        call HDF5SetLimits  (HDF5ID = ObjHDF5, ILB = 1, IUB = Me%IUB + 1, JLB = 1, JUB = Me%JUB + 1, STAT= STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR20'
        
        call HDF5ReadWindow(HDF5ID        = ObjHDF5,                            &
                            GroupName     = 'Grid',                             &
                            Name          = 'Longitude',                        &
                            Array2D       = CoordsX,                            &
                            STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteESRI_GeoTiff - ModulePercentileComputation - ERR30'
        call HDF5ReadWindow(HDF5ID        = ObjHDF5,                            &
                            GroupName     = 'Grid',                             &
                            Name          = 'Latitude',                         &
                            Array2D       = CoordsY,                            &
                            STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteESRI_GeoTiff - ModulePercentileComputation - ERR31'
        
        
        Me%XOrig    = CoordsX  (1,1) ! origin x (bottom left)
        Me%DX       = subtract_with_precision(CoordsX(1,Me%JUB + 1), CoordsX(1,1)) / Me%JUB ! dx
        Me%YOrig    = CoordsY  (1,1) ! origin y (bottom left)
        Me%DY       = subtract_with_precision(CoordsY(Me%IUB + 1,1), CoordsY(1,1)) / Me%IUB ! dy

        deallocate(CoordsX)
        deallocate(CoordsY)
        
    
    end subroutine ReadHD5Grid
    
    !--------------------------------------------------------------------------    
    
    subroutine ReadFileCheckCondition(ParameterX, n)

        !Arguments-------------------------------------------------------------
        type (T_Parameter), pointer         :: ParameterX
        integer                             :: n
        
        !Local-----------------------------------------------------------------
        real,   dimension(:,:,:), pointer   :: Classes3D, ReadMatrix3D
        integer,dimension(:,:,:), pointer   :: WaterPoints3D
        integer,dimension(:,:  ), pointer   :: WaterPoints2D
        real,   dimension(:,:  ), pointer   :: Limits2D, ReadMatrix2D        
        character(len=StringLength)         :: FieldName, LimitsName        
        character(len=PathLength  )         :: GroupName
        integer                             :: STAT_CALL, IUB, JUB, KUB, ObjHDF5
        integer                             :: nItems, ni, NDim
        
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
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR10' 
                                                
        
        call GetHDF5ArrayDimensions(HDF5ID          = ObjHDF5,                          &
                                    GroupName       = trim(GroupName)//'/'//trim(FieldName),& 
                                    ItemName        = trim(FieldName),                  &
                                    OutputNumber    = 1,                                &
                                    NDim            = NDim)        
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR20'
                                                
        
        allocate(Classes3D   (0:IUB+1,0:JUB+1,0:nItems+1))
        allocate(Limits2D    (1:nItems,1:2              ))        
        
        Classes3D   (1:IUB,1:JUB,1:nItems) = 0.
                            
        call HDF5SetLimits  (HDF5ID = ObjHDF5, ILB = 1, IUB = nItems,                   &
                                               JLB = 1, JUB = 2,                        &
                             STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR30'

        call HDF5ReadData   (HDF5ID        = ObjHDF5,                                   &
                             GroupName     = trim(GroupName),                           &
                             Name          = trim(LimitsName),                          &
                             Array2D       = Limits2D,                                  &
                             STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR40'
        
    
        
        if      (NDim == 3) then
            allocate(ReadMatrix3D(0:IUB+1,0:JUB+1,0:KUB+1   ))  
            ReadMatrix3D(1:IUB,1:JUB,1:KUB   ) = 0.
        elseif  (NDim == 2) then
            allocate(ReadMatrix2D(0:IUB+1,0:JUB+1           )) 
            ReadMatrix2D(1:IUB,1:JUB         ) = 0.
        endif
    
        do ni = 1, nItems
            
            if      (NDim == 3) then            
            
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
            
            else if  (NDim == 2) then
                                                    
                call HDF5SetLimits  (HDF5ID = ObjHDF5, ILB = 1, IUB = IUB,              &
                                                       JLB = 1, JUB = JUB,              &
                                     STAT   = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR70'

                call HDF5ReadWindow(HDF5ID        = ObjHDF5,                            &
                                    GroupName     = trim(GroupName)//'/'//trim(FieldName),&
                                    Name          = trim(FieldName),                    &
                                    Array2D       = ReadMatrix2D,                       &
                                    OutputNumber  = ni,                                 &  
                                    STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR80'

                Classes3D(:,:,ni) = ReadMatrix2D(:,:)
                
            endif
            
        enddo
        
        if      (NDim == 3) then
            deallocate(ReadMatrix3D)  
        elseif  (NDim == 2) then
            deallocate(ReadMatrix2D) 
        endif            
        
        allocate(WaterPoints2D(0:IUB+1,0:JUB+1))           
        
        if      (trim(ParameterX%MappingName)=="WaterPoints3D") then 
            
            call HDF5SetLimits  (HDF5ID = ObjHDF5, ILB = 1, IUB = IUB,                  &
                                                   JLB = 1, JUB = JUB,                  &
                                                   KLB = 1, KUB = KUB,                  &  
                                 STAT   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR90'               
        
        allocate(WaterPoints3D(0:IUB+1,0:JUB+1,0:KUB+1   ))         
        
        call HDF5ReadWindow(HDF5ID        = ObjHDF5,                                    &
                            GroupName     = "/Grid",                                    &
                            Name          = "WaterPoints3D",                            &
                            Array3D       = WaterPoints3D,                              &
                            STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR100'  
            
            WaterPoints2D(0:IUB+1,0:JUB+1)  = WaterPoints3D(0:IUB+1,0:JUB+1,KUB) 
            
            deallocate(WaterPoints3D)
            
        else if (trim(ParameterX%MappingName)=="WaterPoints2D") then 

            call HDF5SetLimits  (HDF5ID = ObjHDF5, ILB = 1, IUB = IUB,                  &
                                                   JLB = 1, JUB = JUB,                  &
                                 STAT   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR110'               

            call HDF5ReadWindow(HDF5ID        = ObjHDF5,                                &
                                GroupName     = "/Grid",                                &
                                Name          = "WaterPoints2D",                        &
                                Array2D       = WaterPoints2D,                          &
                                STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadFileCheckCondition - ModulePercentileComputation - ERR120' 
        
        endif    
        
        call CheckCondition (Classes3D, Limits2D, Me%OutMatrix3D, nItems, n, WaterPoints2D)        
        

        deallocate(Classes3D    ) 
        
        deallocate(WaterPoints2D)

    end subroutine ReadFileCheckCondition
    
    !--------------------------------------------------------------------------
    
    subroutine CheckCondition(InMatrix3D, Limits2D, OutMatrix3D, nItems, n, WaterPoints2D)

        !Arguments-------------------------------------------------------------
        real,   dimension(:,:,:), pointer   :: InMatrix3D
        integer,dimension(:,:  ), pointer   :: WaterPoints2D
        real,   dimension(:,:,:), pointer   :: OutMatrix3D
        real,   dimension(:,:  ), pointer   :: Limits2D
        integer                             :: nItems, n
        
        !Local-----------------------------------------------------------------
        real                                :: Sum, dP2, dPx, LimitX, PercentX
        integer                             :: IUB, JUB, ni, i, j !,STAT_CALL, KUB
        real                                :: XPoint, YPoint
        !Begin-----------------------------------------------------------------
        
        
        IUB = Me%IUB
        JUB = Me%JUB
        
        PercentX = 100. - Me%Conditions(n)%Percentile 
        
        
        do i = 1, IUB
        do j = 1, JUB
                                     
iW:         if (WaterPoints2D(i,j) == 1) then
            
                Sum = 0.
            
                do ni = 1, nItems
                
                    Sum = Sum + InMatrix3D(i, j, ni)
                
                    if (Me%Conditions(n)%Lower) then
                
                    if (Sum >= Me%Conditions(n)%Percentile) then 
                    
                        dPx = Sum - Me%Conditions(n)%Percentile
                        dP2 = InMatrix3D(i, j, ni)                    
                    
                        if (dPx == 0. .or. ni == nItems .or. dP2 == 0.) then
                                LimitX = Limits2D(ni,2)
                        else
                                LimitX = (Limits2D(ni,2) - Limits2D(ni,1)) * dPx/dP2 + Limits2D(ni,1)
                        endif
                
                            if (Me%Conditions(n)%Limit >= LimitX) then 
                                OutMatrix3D(i, j, 0) = OutMatrix3D(i, j,0) + 1
                                OutMatrix3D(i, j, n) = 1
                            endif
                        
                            exit
                        
                        endif
                    
                    endif
                
                    
                        if (.not. Me%Conditions(n)%Lower) then
                    
                        if (Sum >= PercentX) then 
                    
                            dPx = Sum - PercentX
                            dP2 = InMatrix3D(i, j, ni)                    
                    
                            if (dPx == 0. .or. ni == nItems .or. dP2 == 0.) then
                                LimitX = Limits2D(ni,2)
                            else
                                LimitX = (Limits2D(ni,2) - Limits2D(ni,2)) * dPx/dP2 + Limits2D(ni,2)
                            endif
                    
                            if (Me%Conditions(n)%Limit <= LimitX) then 
                                OutMatrix3D(i, j, 0) = OutMatrix3D(i, j, 0) + 1
                                OutMatrix3D(i, j, n) = 1
                            endif
                    
                        exit    
                
                    endif
                
                    endif
                
                enddo
                
            endif iW
        enddo
        enddo


    end subroutine CheckCondition
    
    !--------------------------------------------------------------------------  
    
    subroutine CheckExclusionAreas(ParameterX)
    
        !Arguments-------------------------------------------------------------
        type (T_Parameter), pointer         :: ParameterX

        !Local-----------------------------------------------------------------
        integer                             :: IUB, JUB, i, j
        real                                :: XPoint, YPoint
        !Begin-----------------------------------------------------------------
        
        call New(ParameterX%ExclusionAreas,  trim(ParameterX%File))        
        
        do i = 1, Me%IUB
        do j = 1, Me%JUB
            
            if (Me%OutMatrix3D(i, j, 0) > Me%FillValueOut) then
            
                !if(Me%Bathym2D(i,j) < Me%BathymMin .or. Me%Bathym2D(i,j) > Me%BathymMax) then
                !    Me%OutMatrix3D(i, j) = Me%FillValueOut
                !    cycle
                !endif            
            
                if (associated(ParameterX%ExclusionAreas)) then
                    XPoint  = Me%Xorig + Me%DX * (real(j-1) + 0.5)
                    YPoint  = Me%Yorig + Me%DY * (real(i-1) + 0.5)
                    if (PointXYInsidePolySet(ParameterX%ExclusionAreas, XPoint, YPoint)) then
                        Me%OutMatrix3D(i, j, :) = Me%FillValueOut
                    endif
                endif
                
            endif
            
        enddo
        enddo


    end subroutine CheckExclusionAreas
    
    !--------------------------------------------------------------------------  
    
    
    !--------------------------------------------------------------------------  
    
    subroutine BathymFilter()
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                             :: IUB, JUB, i, j
        real                                :: XPoint, YPoint
        !Begin-----------------------------------------------------------------
        
        
        do i = 1, Me%IUB
        do j = 1, Me%JUB
            
            if (Me%OutMatrix3D(i, j, 0) > Me%FillValueOut) then
            
                if(Me%Bathym2D(i,j) < Me%BathymMin .or. Me%Bathym2D(i,j) > Me%BathymMax) then
                    Me%OutMatrix3D(i, j, :) = Me%FillValueOut
                    cycle
                endif            
            
            endif
            
        enddo
        enddo


    end subroutine BathymFilter
    
    !--------------------------------------------------------------------------      
    function subtract_with_precision(x, y) result(res)
        real(kind=c_double), intent(in) :: x, y
        real(kind=c_double) :: res
        integer :: max_decimals

        ! Determine max number of decimals
        max_decimals = max(num_decimals(x), num_decimals(y))

        ! Perform subtraction with appropriate scaling
        res = (x * 10.0_c_double**max_decimals - y * 10.0_c_double**max_decimals) &
              / 10.0_c_double**max_decimals
    end function subtract_with_precision
    
    !function num_decimals(val) result(ndec)
    !    real(kind=c_double), intent(in) :: val
    !    integer :: ndec
    !    real(kind=c_double) :: temp
    !    integer, parameter :: max_iterations = 15
    !
    !    ndec = 0
    !    temp = val
    !
    !    do while (ndec < max_iterations)
    !        if (abs(temp - nint(temp)) < 1.0_c_double * 1.0e-14) exit
    !        temp = temp * 10.0_c_double
    !        ndec = ndec + 1
    !    end do
    !end function num_decimals
    
    function num_decimals(val) result(ndec)
        real(kind=c_double), intent(in) :: val
        integer :: ndec, pos
        character(len=100) :: tempStr
    
        ! Convert the number to a string
        write(tempStr, '(f0.15)') val
    
        ! Search for the decimal point
        pos = index(tempStr, '.')
    
        ! If no decimal point is found, number of decimals is zero
        if (pos == 0) then
            ndec = 0
            return
        end if
    
        ! Determine the length of the decimal part
        ! And make sure to exclude trailing zeros after the decimal point
        ndec = len_trim(tempStr) - pos
        do while(tempStr(pos + ndec:len_trim(tempStr)) == '0' .and. ndec > 0)
            ndec = ndec - 1
        end do
    
    end function num_decimals


    
    subroutine WriteESRI_GeoTiff
        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                             :: i, j, STAT_CALL
        integer                             :: ObjHDF5
        real, dimension(:,:  ), pointer     :: CoordsX, CoordsY
        real, dimension(:,:,:), allocatable :: Aux3D
        type (T_Parameter)    , pointer     :: ParameterX
        real(kind=c_double)                 :: geotrans(6) = 0.0
        character(len=1000)                 :: projref
        integer                             :: Unit ! to write 4 corners
        !Body------------------------------------------------------------------
                
        ! Get one HDF (the first) to check for conditions
        nullify(ParameterX)
        call Search_Parameter(ParameterX, Me%Conditions(1)%PropName, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ComputePercentile - ModulePercentileComputation - ERR10'
        
        ObjHDF5     = ParameterX%ObjHDF5
        
        ! Is and Js from HDF instead of 
        
        ! Read Coordinates from HDF5
        allocate(CoordsX(0:Me%IUB+2,0:Me%JUB+2)); CoordsX = 0.0
        allocate(CoordsY(0:Me%IUB+2,0:Me%JUB+2)); CoordsY = 0.0
        
        ! CoordX and CoordY have one extra size because its the corners of the HDF cells
        call HDF5SetLimits  (HDF5ID = ObjHDF5, ILB = 1, IUB = Me%IUB + 1, JLB = 1, JUB = Me%JUB + 1, STAT= STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteESRI_GeoTiff - ModulePercentileComputation - ERR20'
        
        call HDF5ReadWindow(HDF5ID        = ObjHDF5,                            &
                            GroupName     = 'Grid',                             &
                            Name          = 'Longitude',                        &
                            Array2D       = CoordsX,                            &
                            STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteESRI_GeoTiff - ModulePercentileComputation - ERR30'
        call HDF5ReadWindow(HDF5ID        = ObjHDF5,                            &
                            GroupName     = 'Grid',                             &
                            Name          = 'Latitude',                         &
                            Array2D       = CoordsY,                            &
                            STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteESRI_GeoTiff - ModulePercentileComputation - ERR31'
        
        ! Get GeoTransform from HDF5 (Convert the floating-point resolution to integers by scaling)
        
        geotrans(1) = CoordsX  (1,1) ! origin x (bottom left)
        geotrans(2) = subtract_with_precision(CoordsX(1,Me%JUB + 1), CoordsX(1,1)) / Me%JUB ! dx
        geotrans(3) = 0.0 ! offsetx
        geotrans(4) = CoordsY  (1,1) ! origin y (bottom left)
        geotrans(5) = 0.0 ! offsety
        geotrans(6) = subtract_with_precision(CoordsY(Me%IUB + 1,1), CoordsY(1,1)) / Me%IUB ! dy
        
        ! Get projection - ideally from HDF5 - does it only work in 4326?
        select case (Me%EPSG)
            case (4326)
                projref = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]'
            case default
                write (*,*) 'Unknown projection for GeoTiff'
                stop 'WriteESRI_GeoTiff - ModulePercentileComputation - ERR40'
            end select
            
        ! Apply mask when value below fillval
        do i = Me%IUB, 1, -1
            do j = 1, Me%JUB
                ! Assure fill values
                if (Me%OutMatrix3D(i,j, 0) <= Me%FillValueOut) then
                    Me%OutMatrix3D(i,j,:) = Me%FillValueOut
                endif
            enddo
        enddo
        

        
        ! Set to 3D Matrix
        allocate(Aux3D(1:Me%JUB, 1:Me%IUB, 1:Me%NumberCond+1)); Aux3D = 0
        
        ! Transpose Matrix
        do i=1, Me%IUB
        do j=1, Me%JUB
            Aux3D (j, i,1:Me%NumberCond+1) = Me%OutMatrix3D (i, j, 0:Me%NumberCond)
        enddo
        enddo
        
        ! Create GeoTiff
        call CreateRaster   (FileName       = trim(Me%OutputGeoTiff)  , &
                             DriverName     = "GTiff"                 , &
                             RasterWidth    = Me%JUB                  , &
                             RasterHeight   = Me%IUB                  , &
                             DataType       = GDT_Float64             , &
                             NumberBands    = Me%NumberCond+1         , &
                             Projection     = projref                 , &
                             GeoTrans       = geotrans                , &
                             DataMatrix3D   = Aux3D )
        deallocate(Aux3D)
        
        ! Write corners of file in RasterCorners.txt
        if (Me%WriteCorners) then
            call UnitsManager(Unit, OPEN_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'WriteESRI_GeoTiff - ModulePercentileComputation - ERR100'
            
            open(Unit   = Unit,                                                             &
                 File   = 'RasterCorners.txt',                                              &
                 Form   = 'FORMATTED',                                                      &
                 STATUS = 'UNKNOWN',                                                        &
                 Action = 'WRITE',                                                          &
                 IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'WriteESRI_GeoTiff - ModulePercentileComputation - ERR110'
            
            write(Unit,'(A14,f12.6,A,f12.6)') 'topright:     ', CoordsX(1,Me%JUB + 1),' ', CoordsY(Me%IUB + 1,1)
            write(Unit,'(A14,f12.6,A,f12.6)') 'bottomright:  ', CoordsX(1,Me%JUB + 1),' ', CoordsY(1         ,1)
            write(Unit,'(A14,f12.6,A,f12.6)') 'topleft:      ', CoordsX(1,1         ),' ', CoordsY(Me%IUB + 1,1)
            write(Unit,'(A14,f12.6,A,f12.6)') 'bottomleft:   ', CoordsX(1,1         ),' ', CoordsY(1,1         )
            
            call UnitsManager(Unit, CLOSE_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'WriteESRI_GeoTiff - ModulePercentileComputation - ERR120'
            
        end if 
        
        deallocate(CoordsX)
        deallocate(CoordsY)
        
    
    end subroutine WriteESRI_GeoTiff
    
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
                if (Me%OutMatrix3D(i,j, 0) <= Me%FillValueOut) Me%OutMatrix3D(i,j, 0) = Me%FillValueOut
            enddo
            write(Line,'(4000(I4,1x))') int(Me%OutMatrix3D(i,1:Me%JUB,0))
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
        integer                             :: STAT_, nUsers!, STAT_CALL        

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPercentileComputationID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mGeneric_,  Me%InstanceID)

            if (nUsers == 0) then
                
                !call KillHorizontalGrid(HorizontalGridID = Me%ObjGrid,                          &
                !                        STAT             = STAT_CALL)
                !if (STAT_CALL /= SUCCESS_)                                                      &
                !    stop 'KillPercentileComputation - ModulePercentileComputation - ERR10'
                !
                !call KillGridData(GridDataID       = Me%ObjGridData,                            &
                !                  STAT             = STAT_CALL)
                !if (STAT_CALL /= SUCCESS_)                                                      &
                !    stop 'KillPercentileComputation - ModulePercentileComputation - ERR20'
                
                deallocate(Me%Bathym2D)
                
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








