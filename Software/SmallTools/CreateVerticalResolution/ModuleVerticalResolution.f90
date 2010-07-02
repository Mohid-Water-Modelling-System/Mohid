!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : CreateVerticalResolution
! MODULE        : VerticalResolution
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : January 2005
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Module to create vertical resolution (geometry) based on a property 
!                 vertical profile
!
!------------------------------------------------------------------------------

!DataFile
!
!   PROFILEFILENAME         : char                  [-]         !Path to file with property profile
!
!   OUTPUTFILENAME          : char                  [-]         !Path to output geometry file
!
!   MAXIMUMDEPTH            : real                  [-]         !Maximum bathymetry depth
!
!   LAYERS                  : integer               [-]         !Required number of layers
!                                                               !(specified instead of PROPERTY_STEP)
!
!   PROPERTY_STEP           : real                  [-]         !Required property change per layer
!                                                               !(specified instead of LAYERS)
!
!   DOMAIN_TYPE             : char                  [CARTESIAN] !Vertical domain type intended
!
!   MAXDEPTH_SURFACE        : real                  [100.0]     !Maximum depth for surface layer
!                                                               !(depth until reference level) 
!                                                               !(zero hidrográfico)
!   
!   MAXDEPTH_BOTTOM         : real                  [500.0]     !Maximum depth for bottom layer

Module ModuleVerticalResolution

    use ModuleGlobalData,           only: StringLength, PathLength, SUCCESS_, FromFile,     &
                                          null_str, FromBlock, OPEN_FILE, CLOSE_FILE
    use ModuleEnterData,            only: ConstructEnterData, GetData,                      &
                                          ExtractBlockFromBuffer, Block_Unlock,             &
                                          RewindBuffer, UnitsManager, WriteDataLine

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartVerticalResolution
    private ::      ConstructVerticalResolution
    private ::          ReadKeywords
    private ::      ConstructProfile
    private ::      OpenOutputFile

    !Selector                    
    
    !Modifier
    public  :: ModifyVerticalResolution
    private ::      CheckMaxLayerDepth
    private ::      CorrectForBathymetry
    private ::      CloseOutputFile

    !Destructor
    public  :: KillVerticalResolution                                                     

    !Management
    
    !Interfaces----------------------------------------------------------------

    !Parameter-----------------------------------------------------------------
    
    !Domain types
    integer, parameter                              :: Cartesian = 1
    !integer, parameter                              :: Sigma = 2

    !Block delimitators
    character(LEN = StringLength), parameter :: BeginDepth        = '<BeginDepth>'
    character(LEN = StringLength), parameter :: EndDepth          = '<EndDepth>'
    character(LEN = StringLength), parameter :: BeginProfileValues= '<BeginProfileValues>'
    character(LEN = StringLength), parameter :: EndProfileValues  = '<EndProfileValues>'

    !Types---------------------------------------------------------------------

    type T_Profile
        integer                                     :: ObjProfile = 0
        character(PathLength)                       :: FileName = null_str
        real,           dimension(:), pointer       :: Values, Depths
        integer                                     :: nDepths
    end type T_Profile

    private :: T_VerticalResolution
    type       T_VerticalResolution
        integer                                     :: ObjEnterData = 0
        integer                                     :: ObjOutput = 0
        character(PathLength)                       :: OutputFile, DataFile = null_str
        type(T_Profile)                             :: Profile
        integer                                     :: DomainType
        integer                                     :: NumberLayers
        real                                        :: MaximumDepth
        real                                        :: MaxDepthSurface
        real                                        :: MaxDepthBottom
        real                                        :: PropertyStep
        logical                                     :: DefinitionByNumber
        type(T_VerticalResolution), pointer                      :: Next
    end type  T_VerticalResolution

    !Global Module Variables
    type (T_VerticalResolution), pointer                         :: Me

    integer                                         :: mVerticalResolution_ = 0 !just to compile

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine StartVerticalResolution(ObjVerticalResolutionID, DataFile)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjVerticalResolutionID 
        character(PathLength), intent(IN)               :: DataFile 

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------

        !------------------------------------------------------------------------

        !Assures nullification of the global variable
        allocate(Me)

        !Returns ID
        ObjVerticalResolutionID          = 1

        !Atribute the name of data file            
        Me%DataFile = DataFile

        call ConstructVerticalResolution

        !----------------------------------------------------------------------

    end subroutine StartVerticalResolution

    !--------------------------------------------------------------------------

    subroutine ConstructVerticalResolution

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        ! Read keyword data file
        call ReadKeywords

        ! Read profile values
        call ConstructProfile

        ! Open and write header to Output File
        call OpenOutputFile

        !----------------------------------------------------------------------

    end subroutine ConstructVerticalResolution

    !--------------------------------------------------------------------------

    subroutine ReadKeywords

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: DomainType = null_str
        integer                                     :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------

        call ConstructEnterData (Me%ObjEnterData, Me%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'ReadKeywords - ModuleVerticalResolution - ERR10'

        !Obtain profile file name
        call GetData(Me%Profile%FileName, Me%ObjEnterData, iflag,               &
                     keyword      = 'PROFILEFILENAME',                          &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'VerticalResolution',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                              &
        stop 'ReadKeywords - ModuleVerticalResolution - ERR20'   

        !Obtain output file name
        call GetData(Me%OutputFile, Me%ObjEnterData, iflag,                     &
                     keyword      = 'OUTPUTFILENAME',                           &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'VerticalResolution',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                              &
        stop 'ReadKeywords - ModuleVerticalResolution - ERR30'   

        !Obtain number of layers
        call GetData(Me%NumberLayers, Me%ObjEnterData, iflag,                   &
                     keyword      = 'LAYERS',                                   &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'VerticalResolution',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) then
            !Obtain property step
            call GetData(Me%PropertyStep, Me%ObjEnterData, iflag,               &
                     keyword      = 'PROPERTY_STEP',                            &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'VerticalResolution',                       &
                     STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                          &
            stop 'ReadKeywords - ModuleVerticalResolution - ERR40'

            !Option chosen: calculate layer depth with given max property step            
            Me%DefinitionByNumber = .false.
        
        else

            !Option chosen: calculate layer depth with given layer number            
            Me%DefinitionByNumber = .true.

        endif

        !Obtain domain maximum depth
        call GetData(Me%MaximumDepth, Me%ObjEnterData, iflag,                   &
                     keyword      = 'MAXIMUMDEPTH',                             &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'VerticalResolution',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) then
            write(*,*)'The maximum bathymetry depth is not provided for geometry.'
            stop 'ReadKeywords - ModuleVerticalResolution - ERR50'   
        endif

        !Obtain domain layer type
        call GetData(DomainType, Me%ObjEnterData, iflag,                        &
                     keyword      = 'DOMAIN_TYPE',                              &
                     SearchType   = FromFile,                                   &
                     Default      = trim("CARTESIAN"),                          &
                     ClientModule = 'VerticalResolution',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'ReadKeywords - ModuleVerticalResolution - ERR60'   

        select case (trim(adjustl(DomainType)))
            case ("Cartesian", "CARTESIAN", "cartesian")
                Me%DomainType = Cartesian
            !case ("Sigma", "SIGMA", "sigma")
            !    Me%DomainType    = Sigma
            case default
                write(*,*)'Invalid option for keyword DOMAIN_TYPE'
                stop 'ReadKeywords - ModuleVerticalResolution - ERR70'
        end select

        !Obtain maximum bottom layer depth (meters) 
        call GetData(Me%MaxDepthBottom, Me%ObjEnterData, iflag,                 &
                     keyword      = 'MAXDEPTH_BOTTOM',                          &
                     SearchType   = FromFile,                                   &
                     Default      = 500.0,                                      &
                     ClientModule = 'VerticalResolution',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'ReadKeywords - ModuleVerticalResolution - ERR80'   

        !Obtain maximum bottom layer depth (meters) 
        call GetData(Me%MaxDepthSurface, Me%ObjEnterData, iflag,                &
                     keyword      = 'MAXDEPTH_SURFACE',                         &
                     SearchType   = FromFile,                                   &
                     Default      = 100.0,                                      &
                     ClientModule = 'VerticalResolution',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'ReadKeywords - ModuleVerticalResolution - ERR90'   

        !----------------------------------------------------------------------

    end subroutine ReadKeywords

    !--------------------------------------------------------------------------

    subroutine ConstructProfile

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                     :: ClientNumber
        integer                                     :: FirstLine, LastLine
        logical                                     :: BlockFound
        integer                                     :: line, l, AuxNValues
        integer                                     :: STAT_CALL, iflag
        real                                        :: Aux

        !Begin-----------------------------------------------------------------

        !Open the profile file
        call ConstructEnterData(Me%Profile%ObjProfile, Me%Profile%FileName,     & 
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructProfile - ModuleVerticalResolution - ERR10'

        !Obtain depth values
        call ExtractBlockFromBuffer(Me%Profile%ObjProfile, ClientNumber,        &
                                    BeginDepth, EndDepth, BlockFound,           &
                                    FirstLine = FirstLine, LastLine = LastLine, &
                                    STAT = STAT_CALL)
    
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructProfile - ModuleVerticalResolution - ERR20'

        if (BlockFound) then

            Me%Profile%nDepths =  LastLine - FirstLine - 1

            call GetData(Aux, EnterDataID = Me%Profile%ObjProfile,              & 
                         flag = iflag, SearchType = FromBlock,                  &
                         Buffer_Line = FirstLine+1, STAT = STAT_CALL)

            if(iflag .ne. 1)then
                write(*,*) 'More than one profile is provided in the depth block.'
                write(*,*) 'Only one must be provided.'
                write(*,*) 'FileName = ', trim(Me%Profile%FileName)
                stop 'ConstructProfile - ModuleVerticalResolution - ERR30'

            endif

            !Allocates variable
            allocate (Me%Profile%Depths (1:Me%Profile%nDepths))

            l = 1
            do line = FirstLine + 1, LastLine - 1

                call GetData(Me%Profile%Depths(l),                              &
                             EnterDataID = Me%Profile%ObjProfile, flag = iflag, &
                             SearchType = FromBlock, Buffer_Line = line,        &
                             STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'ConstructProfile - ModuleVerticalResolution - ERR40'

                l = l + 1

            enddo

            call Block_Unlock(Me%Profile%ObjProfile, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructProfile - ModuleVerticalResolution - ERR50'

            call RewindBuffer(Me%Profile%ObjProfile, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructProfile - ModuleVerticalResolution - ERR60'

        else 

            write(*,*) 'Block <BeginDepth>, <EndDepth> not found'
            write(*,*) 'FileName = ', trim(Me%Profile%FileName)
            stop 'ConstructProfile - ModuleVerticalResolution - ERR70'

        endif

        !Obtain profile values
        call ExtractBlockFromBuffer(Me%Profile%ObjProfile, ClientNumber,        &
                                    BeginProfileValues, EndProfileValues,       &
                                    BlockFound, FirstLine = FirstLine,          &
                                    LastLine = LastLine,                        &
                                    STAT = STAT_CALL)
    
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructProfile - ModuleVerticalResolution - ERR80'

        if (BlockFound) then

            call GetData(Aux, EnterDataID = Me%Profile%ObjProfile,              &
                         flag = iflag, SearchType = FromBlock,                  &
                         Buffer_Line = FirstLine+1, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructProfile - ModuleVerticalResolution - ERR90'

            if(iflag .ne. 1)then
                write(*,*) 'More than one profile is provided in the values block.'
                write(*,*) 'Only one must be provided.'
                write(*,*) 'FileName = ', trim(Me%Profile%FileName)
                stop 'ConstructProfile - ModuleVerticalResolution - ERR100'
            endif

            AuxNValues =  LastLine - FirstLine - 1

            if (AuxNValues .ne. Me%Profile%nDepths) then
                write(*,*) 'The number of profile values is not equal to number of depths.'
                write(*,*) 'FileName = ', trim(Me%Profile%FileName)
                stop 'ConstructProfile - ModuleVerticalResolution - ERR110'
            endif

            !Allocates variable
            allocate (Me%Profile%Values (1:Me%Profile%nDepths))

            l = 1
            do line = FirstLine + 1, LastLine - 1

                call GetData(Me%Profile%Values(l),                              &
                             EnterDataID = Me%Profile%ObjProfile, flag = iflag, &
                             SearchType = FromBlock, Buffer_Line = line,        & 
                             STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'ConstructProfile - ModuleVerticalResolution - ERR120'

                l = l + 1

            enddo

            call Block_Unlock(Me%Profile%ObjProfile, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructProfile - ModuleVerticalResolution - ERR130'

        else 

            write(*,*) 'Block <BeginProfileValues>, <EndProfileValues> not found'
            write(*,*) 'FileName = ', trim(Me%Profile%FileName)
            stop 'ConstructProfile - ModuleVerticalResolution - ERR140'

        endif

        !----------------------------------------------------------------------

    end subroutine ConstructProfile

    !--------------------------------------------------------------------------

    subroutine OpenOutputFile

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        character(len=StringLength)                 :: AuxDomainType = null_str

        !Begin-----------------------------------------------------------------

        !Create output file
        call UnitsManager (Me%ObjOutput, OPEN_FILE, STAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OpenOutputFile - ModuleVerticalResolution - ERR10'

        open(Me%ObjOutput, file = Me%OutputFile, Status = 'unknown')

        !Write header
        call WriteDataLine(Me%ObjOutput, "<begindomain>")
        call WriteDataLine(Me%ObjOutput, "ID", 1)

        select case (Me%DomainType)
            case (1)
                AuxDomainType = "CARTESIAN"
            !case (2)
            !    AuxDomainType = "SIGMA"
        end select

        call WriteDataLine(Me%ObjOutput, "TYPE", trim(adjustl(AuxDomainType)))
        call WriteDataLine(Me%ObjOutput, "DOMAINDEPTH", -99.00)
        call WriteDataLine(Me%ObjOutput,  "<<beginlayers>>")

        !----------------------------------------------------------------------

    end subroutine OpenOutputFile

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyVerticalResolution

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: l, lastdepth
        real                                        :: PropStep, LastAuxDepth
        real                                        :: LayerDepth, LastLayerDepth
        real                                        :: AuxValue, AuxDepth
        real                                        :: SumDepths, ProfileStep
        real                                        :: DepthAvailable, addfactor
        integer                                     :: Auxd, d, li
        integer                                     :: lastlayer
        real                                        :: Increment
        real                                        :: MaxLayerDepth

        !Begin-----------------------------------------------------------------

        !Calculate the actual limit value of the profile 
        !(when profile is deeper than maximum bathymetry)
        if (Me%MaximumDepth < Me%Profile%Depths(Me%Profile%nDepths)) then

            write(*,*)
            write(*,*) 'Property profile is deeper than bathymetry.'
            write(*,*) 'New profile end value is calculated interpolating'
            write(*,*) 'for maximum bathymetry depth.'

do1:        do d = 1, Me%Profile%nDepths-1
            
                Auxd = Me%Profile%nDepths-d

                if (Me%MaximumDepth >= Me%Profile%Depths(Auxd)) then

                    !(pseudo profile number of depths and profile values are created)
                    Me%Profile%Values(Auxd+1) = Me%Profile%Values(Auxd) +       &
                        (Me%MaximumDepth - Me%Profile%Depths(Auxd))*            & 
                        (Me%Profile%Values(Auxd) - Me%Profile%Values(Auxd+1))/  &
                        (Me%Profile%Depths(Auxd+1) - Me%Profile%Depths(Auxd))

                    Me%Profile%Depths(Auxd+1) = Me%MaximumDepth

                    Me%Profile%nDepths = Auxd + 1
                    
                    exit do1

                endif

            enddo do1

        endif 


        !Initializate value for layer limit
        AuxValue = Me%Profile%Values(1)

        !Initializate sum of layer depths
        SumDepths = 0.0

        lastdepth = 2

        lastlayer = 0

        !Initializate depth available from previous layer definition 
        AuxDepth = 0.0


        if (Me%DefinitionByNumber) then

            !Calculate property step according to layer number
            PropStep = (Me%Profile%Values(Me%Profile%nDepths)                   &
                        - Me%Profile%Values(1))/Me%NumberLayers
            !(inform the user of property step)
            write(*,*)
            write(*,fmt=1000) PropStep 
1000        format('Layer depths defined for property change of: ', f12.6)
            write(*,*) '(average between profile lower depth and higher depth values)' 
 
        else
        
            PropStep = Me%PropertyStep
            
        endif        


        !Calculate depth available in profile for layer depths
        DepthAvailable = Me%Profile%Depths(Me%Profile%nDepths) -                &
                         Me%Profile%Depths(1)

1001    format(f8.3)

        l = 0

        !Calculation of layer depths given layer property step
do2:    do while (SumDepths < DepthAvailable)

            l = l + 1

            !Find where layer limit is located in profile
do3:        do d = lastdepth, Me%Profile%nDepths

                if (.not. Me%DefinitionByNumber) then

                    ProfileStep = Me%Profile%Values(d) -                        &
                        Me%Profile%Values(lastdepth - 1)

                    !Calculate the value for this layer limit
                    if (ProfileStep > 0) then

                        addfactor = 1.0

                    else

                        addfactor = -1.0
            
                    endif

                else 

                    addfactor = 1.0

                endif
                
                AuxValue = AuxValue + addfactor*PropStep

                if (((AuxValue .le. Me%Profile%Values(d)) .and.                 &
                    (AuxValue .ge. Me%Profile%Values(d-1))) .or.                &
                    ((AuxValue .ge. Me%Profile%Values(d)) .and.                 &
                    (AuxValue .le. Me%Profile%Values(d-1)))) then
                    !Layer limit is located in profile interval before depth d

                    if (l > 1) then

                        if (d /= lastdepth) then
                            !profile interval changed from previous layer definition
                            AuxDepth = AuxDepth + (Me%Profile%Depths(lastdepth) & 
                            - SumDepths - Me%Profile%Depths(1))
                            !(the remaining of previous profile interval is accounted)
                    
                        else
                            !profile interval is the same used for previous layer
                            AuxDepth = - (LastLayerDepth - LastAuxDepth)
                            !(depth occupied by previous layer is discounted)
                        
                        endif  
                        
                    endif

                    !Calculate layer depth
                    Increment = (Me%Profile%Values(d-1) - AuxValue)*            &
                        (Me%Profile%Depths(d) - Me%Profile%Depths(d-1))         &
                        /(Me%Profile%Values(d-1) - Me%Profile%Values(d))
                    !(for comprehension)
                    
                    LayerDepth = AuxDepth + Increment 

                    SumDepths = SumDepths + LayerDepth

                    LastLayerDepth = LayerDepth

                    if ((l == 1) .and.  (.not. Me%DefinitionByNumber)) then
                        !only checked for PROPERTY_STEP option

                        call CheckMaxLayerDepth(l, LayerDepth, Me%MaxDepthSurface)

                    endif

                    if ((SumDepths == DepthAvailable) .or.                      &
                        ((Me%DefinitionByNumber) .and. (l == Me%NumberLayers))) then
                        !this is the last layer

                        if (Me%MaximumDepth >                                   &
                            Me%Profile%Depths(Me%Profile%nDepths)) then
                            !Add bathymetry excedent to last layer depth
                            call CorrectForBathymetry(LayerDepth, SumDepths)

                        endif

                        if (.not. Me%DefinitionByNumber) then
                            !only checked for PROPERTY_STEP option
                            call CheckMaxLayerDepth(l, LayerDepth, Me%MaxDepthBottom)
                        endif

                    endif

                    do li = lastlayer + 1, l

                        !Write depths to output file
                        write(Me%ObjOutput, fmt=1001) LayerDepth

                    enddo

                    LastAuxDepth = AuxDepth
                    !(depth used to discount this layer depth when next d == lastdepth)

                    AuxDepth = 0.0

                    lastdepth = d

                    lastlayer = l

                    exit do3

                else 

                    if (d == Me%Profile%nDepths) then             
                    !(profile end is reached: nowhere else to go)
                                        
                        LayerDepth = Me%Profile%Depths(Me%Profile%nDepths)      &
                                     - SumDepths - Me%Profile%Depths(1)  

                        SumDepths = SumDepths + LayerDepth                         

                        if (Me%MaximumDepth >                                   &
                            Me%Profile%Depths(Me%Profile%nDepths)) then
                            !Add bathymetry excedent to last layer depth
                            call CorrectForBathymetry(LayerDepth, SumDepths)

                        endif

                        if (l == 1) then
                            !Only one layer!
                            !Find suitable maximum depth
                            !(minimum of MaxDepthSurface and MaxDepthBottom) 
                            MaxLayerDepth = min(Me%MaxDepthBottom,              &
                                                Me%MaxDepthSurface)

                        else

                            MaxLayerDepth = Me%MaxDepthBottom

                        endif

                        if (.not. Me%DefinitionByNumber) then 
                            !only checked for PROPERTY_STEP option
                            call CheckMaxLayerDepth(l, LayerDepth, MaxLayerDepth)
                        endif
                          
                        do li = lastlayer + 1, l

                            !Write depths to output file
                            write(Me%ObjOutput, fmt=1001) LayerDepth

                        enddo

                        exit do3
                    
                    elseif ((l == 1) .or. (d > lastdepth)) then
                 
                        AuxDepth = AuxDepth + (Me%Profile%Depths(d) -           &
                           Me%Profile%Depths(d-1)) 

                    endif

                    AuxValue = AuxValue - addfactor*PropStep

                endif

            enddo do3                    

        enddo do2

        Me%NumberLayers = l

        !Terminate output file
        call CloseOutputFile
        
        !----------------------------------------------------------------------

    end subroutine ModifyVerticalResolution

    !--------------------------------------------------------------------------

    subroutine CheckMaxLayerDepth(layer, LayerDepth, MaxLayerDepth)

        !Arguments-------------------------------------------------------------
        integer                                     :: layer
        real                                        :: LayerDepth, MaxLayerDepth
                                                    
        !Local-----------------------------------------------------------------
        integer                                     :: f, auxlayer
        real                                        :: AuxLayerDepth, factor

        !Begin-----------------------------------------------------------------

        f = 1

        auxlayer = layer
        
        if (layer == 1) then
            !(first layer depth checked until the hydrographic zero)  
            factor = 1.0
        else
            factor = 0.0

        endif

        AuxLayerDepth = LayerDepth 

        do while ((AuxLayerDepth + factor*Me%Profile%Depths(1))                 &
                  > MaxLayerDepth)

            f = f + 1
        
            AuxLayerDepth = LayerDepth/f                          

            layer = layer + 1

        enddo          

        if (Me%DefinitionByNumber) then

            Me%NumberLayers = Me%NumberLayers + f - 1

        endif
        
        if (f > 1) then
            write(*,*)
            write(*,fmt=1002) auxlayer, LayerDepth, MaxLayerDepth        
1002        format('Layer ', i3, ' depth (', f8.3,') bigger than allowed (', f8.3, ').')
            write(*,fmt=1003) f, AuxLayerDepth
1003        format('Substituted by ', i3,' layers with depth ',f8.3) 
        
            LayerDepth = AuxLayerDepth

        endif

        !----------------------------------------------------------------------

    end subroutine CheckMaxLayerDepth

    !--------------------------------------------------------------------------

    subroutine CorrectForBathymetry(LayerDepth, SumDepths)

        !Arguments-------------------------------------------------------------
        real                                        :: LayerDepth, SumDepths

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        !The excedent bathymetry is added to last layer depth
        LayerDepth = LayerDepth + (Me%MaximumDepth - SumDepths -                &
                     Me%Profile%Depths(1))

        SumDepths = SumDepths + LayerDepth

        !----------------------------------------------------------------------

    end subroutine CorrectForBathymetry

    !--------------------------------------------------------------------------

    subroutine CloseOutputFile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call WriteDataLine(Me%ObjOutput, "<<endlayers>>")

        call WriteDataLine(Me%ObjOutput, "LAYERS", Me%NumberLayers)

        call WriteDataLine(Me%ObjOutput, "<enddomain>")

        call UnitsManager(Me%ObjOutput, CLOSE_FILE,                             &
                          STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)
            write(*,*)'Error closing created geometry file',                    &
                       trim(Me%OutputFile)
            write(*,*)'CloseOutputFile - ModuleVerticalResolution - WRN10'
        endif

    end subroutine CloseOutputFile

    !--------------------------------------------------------------------------



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillVerticalResolution

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL

        !------------------------------------------------------------------------

        !Kill profile values
        deallocate(Me%Profile%Values, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'KillVerticalResolution - ModuleVerticalResolution - ERR10'
        nullify(Me%Profile%Values)

        !Kill profile depths
        deallocate(Me%Profile%Depths, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'KillVerticalResolution - ModuleVerticalResolution - ERR20'
        nullify(Me%Profile%Depths)

        !Kill global variable
        deallocate(Me)
        nullify(Me)

        !------------------------------------------------------------------------

    end subroutine KillVerticalResolution
        



    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

end module ModuleVerticalResolution









