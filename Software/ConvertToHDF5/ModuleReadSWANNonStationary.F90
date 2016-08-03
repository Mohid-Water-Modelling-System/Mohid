!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : SWAN
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group & MeteoGalicia
! DATE          : January 2016
! REVISION      : Joao Ribeiro - v2.0
! DESCRIPTION   : Module to convert SWAN NonStationary run table format
!                 files into HDF5 format. For reading into Mohid module
!                  HydrodynamicFile
!
!------------------------------------------------------------------------------


Module ModuleReadSWANNonStationary

    use ModuleTime
    use ModuleGlobalData
    use ModuleFunctions
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleHorizontalGrid
    use ModuleGridData
    use ModuleHorizontalMap

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertSWANNonStationary
    private ::      ReadGlobalOptions
    private ::      ReadProperties
    private ::      ConstructBathymetry
    private ::      Open_HDF5_OutPut_File
    private ::      OutputFields
    private ::      KillSWAN

    !Parameters----------------------------------------------------------------

    character(LEN = StringLength), parameter    :: prop_block_begin   = '<<beginproperty>>'
    character(LEN = StringLength), parameter    :: prop_block_end     = '<<endproperty>>'

    character(LEN = StringLength), parameter    :: units_block_begin  = '<<beginunits>>'
    character(LEN = StringLength), parameter    :: units_block_end    = '<<endunits>>'

    integer,                       parameter    :: NoConversion_          = 0
    
    !Types---------------------------------------------------------------------
    
    type       T_OutPut                                 
         type (T_Time), pointer, dimension(:)   :: OutTime
         integer                                :: Number
         integer                                :: NextOutPut
         logical                                :: ON
         logical                                :: Yes                  = .false.
         integer                                :: TotalOutputs
    end type T_OutPut                                   
    
    private :: T_SWAN
    type       T_SWAN
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ClientNumber
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjBathymetry        = 0
        integer                                 :: ObjHorizontalMap     = 0
        integer                                 :: ObjTime              = 0
        integer                                 :: Unit
        integer                                 :: Nautical
        integer                                 :: WavePower
        integer                                 :: ScaleToHS
        character(len=PathLength)               :: FileName
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: SwanFileName
        character(len=PathLength)               :: OutputFileName
        character(len=StringLength)             :: Generic4DProperty
        character(len=StringLength)             :: Generic4DUnits
        
        integer                                 :: NumberDates          = FillValueInt
        integer                                 :: NumberFields         = FillValueInt
        integer                                 :: NumberProps          = FillValueInt
        integer                                 :: NumberUnits          = FillValueInt
        integer                                 :: NumberCovertions     = FillValueInt

        integer                                 :: DirectionReferential = FillValueInt
        
        character(len=StringLength), dimension(:), pointer :: PropsName
        character(len=StringLength), dimension(:), pointer :: PropsUnits
!        integer,                     dimension(:), pointer :: ConvertProp

        type(T_Time)                            :: BeginTime, EndTime
        real,   dimension(:,:,:),         pointer :: Fields
        real,   dimension(:,:,:),         pointer :: HS 
        real,   dimension(:,:,:),         pointer :: WD 
        real,   dimension(:,:,:),         pointer :: TEX
        real,   dimension(:,:,:),         pointer :: TEY
        real,   dimension(:),           pointer :: PropVector

        integer                                 :: Clientumber

        logical                                 :: WriteVelModulus = .false., WriteWindModulus = .false.
        logical                                 :: ReplaceNonComputedValues
        integer                                 :: ReplacementMethod
        real                                    :: FillValue
        real                                    :: Generic4DValue
        integer                                 :: Generic4DPropertyIndeX
        integer                                 :: Generic4D
        integer                                 :: ReadType
        integer, dimension(:,:  ),  pointer     :: WaterPoints2D
        integer, dimension(:,:  ),  pointer     :: BoundaryFacesU
        integer, dimension(:,:  ),  pointer     :: BoundaryFacesV
        type(T_OutPut)                          :: OutPut
        type(T_Size2D)                          :: WorkSize, Size
    end type  T_SWAN

    type(T_SWAN), pointer              :: Me

    !--------------------------------------------------------------------------
    
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertSWANNonStationary(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !Local-------------------------------------------------------------------
!        integer                                         :: l
        integer                                         :: p
        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        Me%ClientNumber = ClientNumber

        call ReadGlobalOptions

        call ConstructBathymetry

        call ReadProperties

        call Open_HDF5_OutPut_File
        
        do p = 1, Me%NumberProps
            
            call ReadFieldFromFile(p)
            
            if (Me%ReplaceNonComputedValues) then
                call ReplaceNonComputedValues(p)
            endif

            call OutputFields (p)
            
        enddo
        
        call OutputWaveVector
        
        if (Me%WavePower) then
            call CalculateWavePower
        endif
        
        call KillSWAN

        STAT = SUCCESS_

    end subroutine ConvertSWANNonStationary

    !------------------------------------------------------------------------

    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag, iflag1

        !Begin-----------------------------------------------------------------

        call GetData(Me%GridFileName,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'INPUT_GRID_FILENAME',                              &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR10'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR20'

        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'SWAN',                                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR30'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR40'

        call GetData(Me%SwanFileName,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'INPUT_SWAN_FILENAME',                              &
                     ClientModule = 'SWAN',                                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR50'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR60'
       
        call GetData(Me%FillValue,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FILL_VALUE',                                       &
                     default      = FillValueReal,                                      &
                     ClientModule = 'SWAN',                                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR70'
        
        call GetData(Me%Generic4D,                                            &
                     Me%ObjEnterData, iflag,                                  &
                     SearchType   = FromBlock,                                &
                     keyword      = 'Generic4D',                              &
                     default      = 0,                                        &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR80'
        
        call GetData(Me%Generic4DProperty,                                    &
                     Me%ObjEnterData, iflag,                                  &
                     SearchType   = FromBlock,                                &
                     keyword      = 'Generic4D_Property',                     &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR90'
        if (iflag     == 0 .AND. Me%Generic4D == 1 )                          &
                stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR100'

        call GetData(Me%BeginTime,                                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'START',                            &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR110'

        call GetData(Me%EndTime,                                        &
                     Me%ObjEnterData, iflag1,                           &
                     SearchType   = FromBlock,                          &
                     keyword      = 'END',                              &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR120'

        if (iflag==0 .and. iflag1==0)                                   &
                stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR130'
                
        call GetOutPutTime(Me%ObjEnterData,                                   &
                     CurrentTime      = Me%BeginTime,                         &
                     EndTime          = Me%EndTime,                           &
                     keyword          = 'OUTPUT_TIME',                        &
                     SearchType       = FromFile,                             &
                     OutPutsTime      = Me%OutPut%OutTime,                    &
                     OutPutsOn        = Me%OutPut%Yes,                        &
                     OutPutsNumber    = Me%OutPut%TotalOutputs,               &
                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                            &
                stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR140' 

        
        !!also use for stationary where BeginTime == EndTime
        !!writed a dummy OUTPUT_TIME : 0 XXX (where XXX exists so that Output exists
        !!Only one output (the first, the second is dummy with same date)
        if (Me%BeginTime .eq. Me%EndTime) then
            Me%OutPut%TotalOutputs = 1;
        endif
        
        call GetData(Me%Nautical,                                             &
                     Me%ObjEnterData, iflag,                                  &
                     SearchType   = FromBlock,                                &
                     keyword      = 'NAUTICAL',                               &
                     default      = 0,                                        &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)                                            &
         stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR150' 

        call GetData(Me%WavePower,                                            &
                     Me%ObjEnterData, iflag,                                  &
                     SearchType   = FromBlock,                                &
                     keyword      = 'WAVEPOWER',                              &
                     default      = 0,                                        &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)                                            &
         stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR160' 

        call GetData(Me%ScaleToHS,                                            &
                     Me%ObjEnterData, iflag,                                  &
                     SearchType   = FromBlock,                                &
                     keyword      = 'SACLE_TO_HS',                            &
                     default      = 0,                                        &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)                                            &
         stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR170' 
        
        call GetData(Me%ReplaceNonComputedValues,                           &
                     Me%ObjEnterData, iflag,                                  &
                     SearchType   = FromBlock,                                &
                     keyword      = 'REPLACE_NON_COMPUTED_VALUES',            &
                     default      = .TRUE.,                                   &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)                                            &
         stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR180' 
        
        ! 1-Null Value 2-Null Gradient
        call GetData(Me%ReplacementMethod,                                    &
                     Me%ObjEnterData, iflag,                                  &
                     SearchType   = FromBlock,                                &
                     keyword      = 'REPLACEMENT_METHOD',                     &
                     default      = 2,                                        &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)                                            &
         stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR190' 
    end subroutine ReadGlobalOptions

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadProperties

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: l, iflag, STAT_CALL, Line, FirstLine, LastLine
        character(len=StringLength)                 :: AuxChar
        logical                                     :: BlockFound
        !Begin-----------------------------------------------------------------


        call ExtractBlockFromBlock (Me%ObjEnterData,                                    &
                                    ClientNumber      = Me%ClientNumber,                &
                                    block_begin       = prop_block_begin,               &
                                    block_end         = prop_block_end,                 &
                                    BlockInBlockFound = BlockFound,                     &
                                    FirstLine         = FirstLine,                      &
                                    LastLine          = LastLine,                       &
                                    STAT              = STAT_CALL)

cd1 :   if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :       if (.not. BlockFound) then                                                  
                stop 'ReadProperties - ModuleReadSWANNonStationary - ERR10'
            end if cd2
        else cd1
            stop 'ReadProperties - ModuleReadSWANNonStationary - ERR20'
        end if cd1


        Me%NumberProps = LastLine - FirstLine - 1

        allocate(Me%PropsName (Me%NumberProps))  
  
        allocate(Me%PropVector(Me%NumberProps))  


d1:     do l= 1, Me%NumberProps

            line = FirstLine + l

            call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleReadSWANNonStationary - ERR30'

            Me%PropsName(l) = AuxChar

            if (.not. CheckPropertyName (AuxChar)) then
                write(*,*) 'The name ',trim(AuxChar),' is not valid name for the MOHID system'
            endif

        enddo d1

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleReadSWANNonStationary - ERR40'

        call ExtractBlockFromBlock (Me%ObjEnterData,                                    &
                                    ClientNumber      = Me%ClientNumber,                &
                                    block_begin       = units_block_begin,              &
                                    block_end         = units_block_end,                &
                                    BlockInBlockFound = BlockFound,                     &
                                    FirstLine         = FirstLine,                      &
                                    LastLine          = LastLine,                       &
                                    STAT              = STAT_CALL)

cd3 :   if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd4 :       if (.not. BlockFound) then                                                  
                stop 'ReadProperties - ModuleReadSWANNonStationary - ERR50'
            end if cd4
        else cd3
            stop 'ReadProperties - ModuleReadSWANNonStationary - ERR60'
        end if cd3

        Me%NumberUnits = LastLine - FirstLine - 1

        allocate(Me%PropsUnits(Me%NumberUnits))

d2:     do l= 1, Me%NumberUnits

            line = FirstLine + l

            call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleReadSWANNonStationary - ERR70'

            Me%PropsUnits(l) = AuxChar

        enddo d2

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleReadSWANNonStationary - ERR80'

    end subroutine ReadProperties 

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadFieldFromFile(p)

    !Arguments-------------------------------------------------------------
    integer                                     :: p
    !Local----------------------------------------------------------------
    integer                                     :: i, j, n, STAT_CALL,k
    !Begin-----------------------------------------------------------------
    
    if (p == 1) then
        allocate(Me%Fields(Me%OutPut%TotalOutputs, Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate(Me%HS(Me%OutPut%TotalOutputs,Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate(Me%WD(Me%OutPut%TotalOutputs,Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        
        if (Me%WavePower) then
            allocate(Me%TEX(Me%OutPut%TotalOutputs,Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%TEY(Me%OutPut%TotalOutputs,Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        endif
        
        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldFromFile - ModuleReadSWANNonStationary - ERR10'

        open(Unit   = Me%Unit,                                                          &
             File   = trim(Me%SwanFileName),                                                  &
             Form   = 'FORMATTED',                                                      &
             STATUS = 'OLD',                                                            &
             Action = 'READ',                                                           &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldFromFile - ModuleReadSWANNonStationary - ERR20'
    endif
    
    Me%Fields(:,:,:) = 0.
    
    rewind (Me%Unit)
    
    do n=1,Me%OutPut%TotalOutputs
    
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB

            read(Me%Unit,*) Me%PropVector

            if(trim(Me%Generic4DProperty) == trim(Me%PropsName(p))) Then 
            
                Me%Generic4DPropertyIndeX=p
            
            end if
                  
            Me%Fields(n, i, j) = Me%PropVector(p)
            
            if (Me%Fields(n, i, j).eq.-99.0 &
                    .OR. Me%Fields(n, i, j).eq.-999.0) Then
                        
                Me%Fields(n, i, j) = FillValueReal
            endif
            
            if(Me%Generic4D==1 .AND. p==Me%Generic4DPropertyIndeX) Then
                Me%Generic4DValue=Me%PropVector(p)
                Me%Generic4DUnits=trim(Me%PropsUnits(p))
            end if

        enddo
        enddo
             
i1:     if (trim(Me%PropsName(p)) == GetPropertyName(MeanWaveDirection_)) then
            Me%WD(n,:,:) = Me%Fields(n, :, :)
        endif i1
i2:     if (trim(Me%PropsName(p)) == GetPropertyName(SignificantWaveHeight_) .and. Me%ScaleToHS) then
            Me%HS(n,:,:) = Me%Fields(n, :, :)
        else
            Me%HS(n,:,:) = 0. 
        endif i2
           
i3:     if (Me%WavePower) then
i4:         if (trim(Me%PropsName(p)) == GetPropertyName(TransportEnergyX_)) then
                Me%TEX(n,:,:) = Me%Fields(n, :, :)
            endif i4
i5:         if (trim(Me%PropsName(p)) == GetPropertyName(TransportEnergyY_)) then
                Me%TEY(n,:,:) = Me%Fields(n, :, :)
            endif i5
        endif i3
    enddo
    
    if (p == Me%NumberProps) then
        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldFromFile - ModuleReadSWANNonStationary - ERR30'
    endif
    
    end subroutine ReadFieldFromFile
    
    !----------------------------------------------------------------------------
   
    !------------------------------------------------------------------------
    
    subroutine ConstructBathymetry
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !Me%GridFileName="NewGrid.dat_.new"
        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleReadSWANNonStationary - ERR10'

        call GetHorizontalGridSize(Me%ObjHorizontalGrid, Me%Size, Me%WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleReadSWANNonStationary - ERR20'

        call ConstructGridData(Me%ObjBathymetry, Me%ObjHorizontalGrid, FileName = Me%GridFileName,&
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleReadSWANNonStationary - ERR30'


        call ConstructHorizontalMap(Me%ObjHorizontalMap, Me%ObjBathymetry, Me%ObjHorizontalGrid, &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleReadSWANNonStationary - ERR40'

        call GetWaterPoints2D   (Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructBathymetry - ModuleReadSWANNonStationary - ERR50'
        
        call GetBoundaryFaces(Me%ObjHorizontalMap,                         &
                              Me%BoundaryFacesU,                           &
                              Me%BoundaryFacesV,                           &
                              STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_)stop 'ConstructBathymetry - ModuleReadSWANNonStationary - ERR60'

    end subroutine ConstructBathymetry
    
    !------------------------------------------------------------------------
    
    !------------------------------------------------------------------------
    
    subroutine OutputFields(p)

        !Arguments-------------------------------------------------------------
        integer                                         :: p
        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer                   :: Aux2D
        real,    dimension(:), pointer                  :: RealArray1D
        REAL, PARAMETER                                 :: RADE= 180./pi
        integer                                         :: STAT_CALL, i, n
        !Begin-----------------------------------------------------------------


        allocate(Aux2D (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        !allocate(Aux2DX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        !allocate(Aux2DY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        
        if (p == 1) then
        
            if (Me%Nautical) then
                write(*,*)'Swan in Nautical Convention...'
                write(*,*)            
            else
                write(*,*)'Swan in Cartesian Convention...'
                write(*,*)
            endif
            
i2:         if (Me%Generic4D==1) then
    
                allocate(RealArray1D(1:1))
                
                RealArray1D(1) = Me%Generic4DValue
                
                call HDF5SetLimits  (Me%ObjHDF5, 1, 1, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleReadSWANNonStationary - ERR30'
                
                call HDF5WriteData  (Me%ObjHDF5, "/Generic4D",                           &
                                     "Generic4D", Me%Generic4DUnits,                     &
                                     Array1D = RealArray1D,                              &
                                     OutputNumber = n, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleReadSWANNonStationary - ERR40'
                
            endif i2
        endif

        write(*,*)"Converting: "//trim(Me%PropsName(p))
        write(*,*)
        
        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,            &
                           Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleReadSWANNonStationary - ERR60'
            
d1:     do n=1,Me%OutPut%TotalOutputs

            Aux2D(:,:) = Me%Fields(n,:,:)             
             
            call HDF5WriteData(Me%ObjHDF5,                                         &
                        "/Results/"//trim(Me%PropsName(p)),                      &
                        trim(Me%PropsName(p)),                                   &
                        trim(Me%PropsUnits(p)),                                  &
                        Array2D      = Aux2D,                                &
                        OutputNumber = n,                                       &
                        STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleReadSWANNonStationary - ERR70'
        
        enddo d1

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleReadSWANNonStationary - ERR90'

        if (p == 1 .and. Me%Generic4D==1) then
            deallocate(RealArray1D)
            nullify(RealArray1D)
        endif
        
        deallocate(Aux2D)
        nullify(Aux2D)

    end subroutine OutputFields
    
    !------------------------------------------------------------------------
    
    !------------------------------------------------------------------------
    
    subroutine OutputWaveVector

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        !real, dimension(:,:), pointer                   :: Aux2DWD, Aux2DHS
        real, dimension(:,:), pointer                   :: Aux2DX, Aux2DY
        REAL, PARAMETER                                 :: RADE= 180./pi
        real                                            :: Angle,Hsig
        integer                                         :: STAT_CALL, ii, jj
        integer                                         :: n
        character(len=StringLength)                     :: FieldNameX, FieldNameY
        !Begin-----------------------------------------------------------------
        
        !allocate(Aux2DWD(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        !allocate(Aux2DHS(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate(Aux2DX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate(Aux2DY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            
d1:     do n=1,Me%OutPut%TotalOutputs

             !Aux2DWD(:,:) = Me%Fields(kk,WD,:,:)
             !if(Me%ScaleToHS)then
             !   Aux2DHS(:,:) = Me%Fields(kk,HS,:,:)
             !endif             

             Aux2DX(:,:) = 0.
             Aux2DY(:,:) = 0.

d2:         do jj= Me%WorkSize%JLB,Me%WorkSize%JUB  
d3:             do ii= Me%WorkSize%ILB,Me%WorkSize%IUB  

                    if (Me%WaterPoints2D(ii,jj) == WaterPoint) then

                        Angle = Me%WD(n,ii,jj)
                        if(Me%ScaleToHS)then
                            Hsig = Me%HS(n,ii,jj)
                        else
                            Hsig = 1.
                        endif
!                        Angle = MOD(630.-(Aux2D(ii,jj)*(pi/180))*RADE,360.)
                        if (Me%Nautical) then
                            Aux2DX(ii,jj) = Hsig*cos((270-Angle) * Pi / 180.)
                            Aux2DY(ii,jj) = Hsig*sin((270-Angle) * Pi / 180.)
                        else
                            Aux2DX(ii,jj) = Hsig*cos(Angle * Pi / 180.) 
                            Aux2DY(ii,jj) = Hsig*sin(Angle * Pi / 180.) 
                        endif

                    endif

                enddo d3
            enddo d2
            
            FieldNameX = GetPropertyName(WaveX_)
            FieldNameY = GetPropertyName(WaveY_)
        
            call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,            &
                               Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputWaveVector - ModuleReadSWANNonStationary - ERR10'

            call HDF5WriteData(Me%ObjHDF5,                                      &
                               "/Results/"//FieldNameX,                         &
                               FieldNameX,                                       &
                               '-',                                             &
                               Array2D      = Aux2DX,                           &
                               OutputNumber = n,                               &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputWaveVector - ModuleReadSWANNonStationary - ERR20'

            call HDF5WriteData(Me%ObjHDF5,                                      &
                               "/Results/"//FieldNameY,                         &
                               FieldNameY,                                       &
                               '-',                                             &
                               Array2D      = Aux2DY,                           &
                               OutputNumber = n,                               &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputWaveVector - ModuleReadSWANNonStationary - ERR30'

        enddo d1

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputWaveVector - ModuleReadSWANNonStationary - ERR40'

        !deallocate(Aux2DWD)
        !deallocate(Aux2DHS)
        deallocate(Aux2DX)
        deallocate(Aux2DY)
        nullify(Aux2DX, Aux2DY)

    end subroutine OutputWaveVector
    
    !------------------------------------------------------------------------
    
    !------------------------------------------------------------------------
    
    subroutine CalculateWavePower

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer                   :: Aux2D
        !real, dimension(:,:), pointer                   :: Aux2DX, Aux2DY
        REAL, PARAMETER                                 :: RADE= 180./pi
        integer                                         :: STAT_CALL, ii, jj
        integer                                         :: n
        character(len=StringLength)                     :: FieldName
        !Begin-----------------------------------------------------------------
        
        allocate(Aux2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

d1:     do n=1,Me%OutPut%TotalOutputs
d2:         do jj= Me%WorkSize%JLB,Me%WorkSize%JUB  
d3:         do ii= Me%WorkSize%ILB,Me%WorkSize%IUB  

                if (Me%WaterPoints2D(ii,jj) == WaterPoint) then
                    Aux2D(ii,jj) = SQRT((Me%TEX(n,ii,jj)*Me%TEX(n,ii,jj)+Me%TEY(n,ii,jj)*Me%TEY(n,ii,jj)))
                endif

            enddo d3
            enddo d2
            
            FieldName = GetPropertyName(WavePower_)
            
            call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,            &
                                Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputWaveVector - ModuleReadSWANNonStationary - ERR10'

            call HDF5WriteData(Me%ObjHDF5,                                      &
                                "/Results/"//FieldName,                          &
                                FieldName,                                       &
                                '-',                                             &
                                Array2D      = Aux2D,                            &
                                OutputNumber = n,                               &
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'CalculateWavePower - ModuleReadSWANNonStationary - ERR20'
        enddo d1

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'CalculateWavePower - ModuleReadSWANNonStationary - ERR30'

        deallocate(Aux2D)
        nullify(Aux2D)

    end subroutine CalculateWavePower

   !----------------------------------------------------------------------
    subroutine ReplaceNonComputedValues(p)

    !Arguments-------------------------------------------------------------
    integer                                     :: p
    !Local----------------------------------------------------------------
    integer                                     :: i, j, n, STAT_CALL
    real                                        :: a1, a2, a3, a4, atotal
    !Begin-----------------------------------------------------------------
    
     do n=1,Me%OutPut%TotalOutputs   
    
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            
            if (Me%Fields(n, i, j).eq.FillValueReal) then
                        
                if (Me%WaterPoints2D(i, j) == 1) then
                                                
                    if (Me%ReplacementMethod == 1) then !NullValue
                        Me%Fields(n, i, j) = 0.
                            
                    elseif (Me%ReplacementMethod == 2) then !NullGradient
                            
                        a1 = 0; a2 = 0; a3 = 0; a4 = 0
                            
                        if (Me%BoundaryFacesU(i, j) == Not_Boundary .and. Me%WaterPoints2D(i, j-1) == 1 .and. &                                
                            Me%Fields(n, i, j-1).ne.FillValueReal) then
                                
                            a1 = 1
                                
                        elseif(Me%BoundaryFacesU(i, j+1) == Not_Boundary .and. Me%WaterPoints2D(i, j+1) == 1 .and. &
                                Me%Fields(n, i, j+1).ne.FillValueReal) then
                                
                            a2 = 1
                                
                        endif
                            
                        if (Me%BoundaryFacesV(i, j) == Not_Boundary .and. Me%WaterPoints2D(i-1, j) == 1 .and. &
                            Me%Fields(n, i-1, j).ne.FillValueReal) then
                                
                            a3 = 1
                                
                        elseif(Me%BoundaryFacesV(i+1, j) == Not_Boundary .and. Me%WaterPoints2D(i+1, j) == 1 .and. &
                                Me%Fields(n, i+1, j).ne.FillValueReal) then
                                
                            a4 = 1
                                
                        endif
                                  
                        atotal = (a1 + a2 + a3 + a4)
                            
                        if (atotal > 0) then
                            Me%Fields(n, i, j) = (a1*Me%Fields(n, i, j-1) + a2*Me%Fields(n, i, j+1) + &
                                                        a3*Me%Fields(n, i-1, j) + a4*Me%Fields(n, i+1, j))/atotal
                        else
                            Me%Fields(n, i, j) = 0.
                        endif
                    endif
                endif                   
            endif
        enddo
    enddo
    enddo         

    end subroutine ReplaceNonComputedValues
    
    !----------------------------------------------------------------------------
   !------------------------------------------------------------------------
    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        real,    dimension(6), target               :: AuxTime
        real,    dimension(:), pointer              :: TimePtr
        real,    dimension(:,:), pointer            :: Bathymetry
        integer                                     :: STAT_CALL,i
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        call GetGridData        (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR10'

      
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR30'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR40'

            
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",                   &
                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR50'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR60'            
   
       
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR70'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",                &
                              Array2D = Me%WaterPoints2D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR80'


        call UnGetGridData      (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR90'


        do i=1,Me%Output%TotalOutputs

             call ExtractDate   (Me%Output%OutTime(i),                               &
                                AuxTime(1), AuxTime(2), AuxTime(3),                 &
                                AuxTime(4), AuxTime(5), AuxTime(6))

             TimePtr => AuxTime

             call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
             if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR100'


             call HDF5WriteData  (Me%ObjHDF5, "/Time",                               &
                                 "Time", "YYYY/MM/DD HH:MM:SS",                     &
                                 Array1D = TimePtr,                                 &
                                 OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR110'
        
        enddo
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR120'

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

   
    !--------------------------------------------------------------------------

    
    subroutine KillSWAN
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call UnGetHorizontalMap (Me%ObjHorizontalMap, Me%BoundaryFacesU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleReadSWANNonStationary - ERR01'
        
        call UnGetHorizontalMap (Me%ObjHorizontalMap, Me%BoundaryFacesV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleReadSWANNonStationary - ERR02'
        
        call UnGetHorizontalMap (Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleReadSWANNonStationary - ERR10'

        call KillHorizontalMap(Me%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleReadSWANNonStationary - ERR20'

        call KillGridData(Me%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleReadSWANNonStationary - ERR30'

        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleReadSWANNonStationary - ERR40'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillSWAN - ModuleReadSWANNonStationary - ERR60'

        deallocate(Me%OutPut%OutTime )
        nullify   (Me%OutPut%OutTime )

        deallocate(Me%PropsName )
        nullify   (Me%PropsName )

        deallocate(Me%PropsUnits)
        nullify   (Me%PropsUnits)

!        deallocate(Me%ConvertProp)
!        nullify   (Me%ConvertProp)

        deallocate(Me%Fields    )
        nullify   (Me%Fields    )
        
        deallocate(Me%HS    )
        nullify   (Me%Hs    )
        deallocate(Me%WD    )
        nullify   (Me%WD    )
        
        if (Me%WavePower) then
            deallocate(Me%TEX    )
            nullify   (Me%TEX    )
            deallocate(Me%TEY    )
            nullify   (Me%TEY    )
        endif

        deallocate(Me%PropVector)
        nullify   (Me%PropVector)

        deallocate(Me)
        nullify   (Me)
    
    end subroutine KillSWAN

    !--------------------------------------------------------------------------
 
end module ModuleReadSWANNonStationary
