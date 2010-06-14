!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : SWAN
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group & MeteoGalicia
! DATE          : September 2003
! REVISION      : Pablo Carracedo - v4.0
! DESCRIPTION   : Module to convert SWAN files into HDF5 format.
!                 For reading into Mohid module HydrodynamicFile
!
!------------------------------------------------------------------------------


Module ModuleSWAN

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
    public  :: ConvertSWAN
    private ::      ReadGlobalOptions
    private ::      ReadDates
    private ::      ReadProperties
    private ::      ReadFields
    private ::      ConstructGlobalOutput
    private ::      ConstructBathymetry
    private ::      Open_HDF5_OutPut_File
    private ::          OutputFields
    private ::      KillSWAN

    !Parameters----------------------------------------------------------------

    character(LEN = StringLength), parameter    :: field_block_begin  = '<<beginfield>>'
    character(LEN = StringLength), parameter    :: field_block_end    = '<<endfield>>'

    character(LEN = StringLength), parameter    :: prop_block_begin   = '<<beginproperty>>'
    character(LEN = StringLength), parameter    :: prop_block_end     = '<<endproperty>>'

    character(LEN = StringLength), parameter    :: units_block_begin  = '<<beginunits>>'
    character(LEN = StringLength), parameter    :: units_block_end    = '<<endunits>>'

    character(LEN = StringLength), parameter    :: date_block_begin   = '<<begindate>>'
    character(LEN = StringLength), parameter    :: date_block_end     = '<<enddate>>'

    character(LEN = StringLength), parameter    :: convert_block_begin= '<<beginconvert>>'
    character(LEN = StringLength), parameter    :: convert_block_end  = '<<endconvert>>'


    integer,                       parameter    :: NoConversion_          = 0
    integer,                       parameter    :: CartToNauticalDegrees_ = 1
    !Types---------------------------------------------------------------------
    
    type       T_OutPut                                 
         type (T_Time), pointer, dimension(:)               :: OutTime
         integer                                            :: Number
         integer                                            :: NextOutPut
         logical                                            :: ON
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
        character(len=PathLength)               :: FileName
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: OutputFileName

        integer                                 :: NumberDates          = FillValueInt
        integer                                 :: NumberFields         = FillValueInt
        integer                                 :: NumberProps          = FillValueInt
        integer                                 :: NumberUnits          = FillValueInt

        integer                                 :: DirectionReferential = FillValueInt

        character(len=StringLength), dimension(:), pointer :: PropsName
        character(len=StringLength), dimension(:), pointer :: PropsUnits
        character(len=PathLength),   dimension(:), pointer :: FilesName
        integer,                     dimension(:), pointer :: ConvertProp

        type(T_Time),     dimension(:), pointer :: InputDates
        real,     dimension(:,:,:),     pointer :: Fields 
        real,     dimension(:),         pointer :: PropVector

        integer                                 :: Clientumber

        logical                                 :: WriteVelModulus = .false., WriteWindModulus = .false.
        real                                    :: FillValue
        integer                                 :: ReadType
        integer, dimension(:,:  ),  pointer     :: WaterPoints2D
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

    subroutine ConvertSWAN(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !Local-------------------------------------------------------------------
        integer                                         :: l, p
        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        Me%ClientNumber = ClientNumber

        call ReadGlobalOptions

        call ConstructBathymetry

        call ReadProperties

        call ReadConversion

        call ReadDates

        call ReadFields

        call ConstructGlobalOutput

        call Open_HDF5_OutPut_File


d1:     do l = 1, Me%NumberDates

            call ReadFieldFromFile(l)

d2:         do p=1, Me%NumberProps

                call OutputFields     (p)

                !if (Me%WriteVelModulus) then
                    !call WriteVelocityModulus(VelocityU_, VelocityV_, VelocityModulus_)
                !endif

                !if (Me%WriteWindModulus) then
                    !call WriteVelocityModulus(WindVelocityX_, WindVelocityY_, WindModulos_)
                !endif

            enddo d2

        enddo d1

  
        call KillSWAN


        STAT = SUCCESS_


    end subroutine ConvertSWAN

    !------------------------------------------------------------------------

    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------

        call GetData(Me%GridFileName,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'INPUT_GRID_FILENAME',                              &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleSWAN - ERR10'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleSWAN - ERR20'

        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleSWAN - ERR30'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleSWAN - ERR40'

       
        call GetData(Me%FillValue,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FILL_VALUE',                                       &
                     default      = -99.999900,                                         &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleSWAN - ERR70'

 
    end subroutine ReadGlobalOptions

    !--------------------------------------------------------------------------

    subroutine ReadDates

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_Time)                                :: AuxTime
        real,       dimension(6)                    :: Aux6
        integer                                     :: l, iflag, STAT_CALL, Line, FirstLine, LastLine
        logical                                     :: BlockFound

        !Begin-----------------------------------------------------------------


        call ExtractBlockFromBlock (Me%ObjEnterData,                                &
                                    ClientNumber      = Me%ClientNumber,            &
                                    block_begin       = date_block_begin,           &
                                    block_end         = date_block_end,             &
                                    BlockInBlockFound = BlockFound,                 &
                                    FirstLine         = FirstLine,                  &
                                    LastLine          = LastLine,                   &
                                    STAT              = STAT_CALL)

cd1 :   if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :       if (.not. BlockFound) then                                                  
                stop 'ReadDates - ModuleSWAN - ERR10'
            end if cd2
        else cd1
            stop 'ReadDates - ModuleSWAN - ERR20'
        end if cd1


        Me%NumberDates = LastLine - FirstLine - 1

        allocate(Me%InputDates(Me%NumberDates))  


d3:     do l= 1, Me%NumberDates

            Line = FirstLine + l

            call GetData(Aux6, EnterDataID = Me%ObjEnterData, flag = iflag,             &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadDates - ModuleSWAN - ERR30'

            call SetDate(AuxTime          ,                                             &
                        Year    = Aux6(1) ,                                             &
                        Month   = Aux6(2) ,                                             & 
                        Day     = Aux6(3) ,                                             &
                        Hour    = Aux6(4) ,                                             &
                        Minute  = Aux6(5) ,                                             &
                        Second  = Aux6(6) )  

            Me%InputDates(l) = AuxTime

        enddo d3

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDates - ModuleSWAN - ERR40'


    end subroutine ReadDates 

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
                stop 'ReadProperties - ModuleSWAN - ERR10'
            end if cd2
        else cd1
            stop 'ReadProperties - ModuleSWAN - ERR20'
        end if cd1


        Me%NumberProps = LastLine - FirstLine - 1

        allocate(Me%PropsName (Me%NumberProps))  
  
        allocate(Me%PropVector(Me%NumberProps))  


d1:     do l= 1, Me%NumberProps

            line = FirstLine + l

            call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleSWAN - ERR30'


            Me%PropsName(l) = AuxChar

            if (.not. CheckPropertyName (AuxChar)) then
                write(*,*) 'The name ',trim(AuxChar),' is not valid name for the MOHID system'
            endif

        enddo d1

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleSWAN - ERR40'


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
                stop 'ReadProperties - ModuleSWAN - ERR50'
            end if cd4
        else cd3
            stop 'ReadProperties - ModuleSWAN - ERR60'
        end if cd3


        Me%NumberUnits = LastLine - FirstLine - 1

        allocate(Me%PropsUnits(Me%NumberUnits))

d2:     do l= 1, Me%NumberUnits

            line = FirstLine + l

            call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleSWAN - ERR70'

            Me%PropsUnits(l) = AuxChar

        enddo d2

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleSWAN - ERR80'


    end subroutine ReadProperties 

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadConversion

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: l, iflag, STAT_CALL, ConversionType, j
        character(len=StringLength)                 :: PropName
        logical                                     :: BlockFound, PropNameWrong
        !Begin-----------------------------------------------------------------


        allocate(Me%ConvertProp(Me%NumberProps))  
  
        Me%ConvertProp(:) = NoConversion_


d1:     do l= 1, Me%NumberProps


            call ExtractBlockFromBlock (Me%ObjEnterData,                                 &
                                        ClientNumber      = Me%ClientNumber,             &
                                        block_begin       = convert_block_begin,         &
                                        block_end         = convert_block_end,           &
                                        BlockInBlockFound = BlockFound,                  &
                                        STAT              = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :           if (.not. BlockFound) then                                                  
                    exit
                end if cd2
            else cd1
                stop 'ReadConversion - ModuleSWAN - ERR10'
            end if cd1

            call GetData(PropName,                                                      &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NAME',                                         &
                         ClientModule = 'ConvertToHDF5',                                &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadConversion - ModuleSWAN - ERR20'

            PropNameWrong = .true.
            do j = 1, Me%NumberProps
                if (trim(Me%PropsName(j))==trim(PropName)) then
                    PropNameWrong = .false. 
                    exit
                endif
            enddo
            
            if (PropNameWrong) stop 'ReadConversion - ModuleSWAN - ERR30'

            call GetData(ConversionType,                                                &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'CONVERSION',                                   &
                         ClientModule = 'ConvertToHDF5',                                &
                         default      = NoConversion_,                                  &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadConversion - ModuleSWAN - ERR40'

            if (ConversionType /= NoConversion_ .and. ConversionType /= CartToNauticalDegrees_) then
                stop 'ReadConversion - ModuleSWAN - ERR50'
            endif

            Me%ConvertProp(j) = ConversionType

        enddo d1

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadConversion - ModuleSWAN - ERR60'


    end subroutine ReadConversion 

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadFields

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: l, iflag, STAT_CALL, Line, FirstLine, LastLine
        character(len=StringLength)                 :: AuxChar
        logical                                     :: BlockFound

        !Begin-----------------------------------------------------------------


        call ExtractBlockFromBlock (Me%ObjEnterData,                                &
                                    ClientNumber      = Me%ClientNumber,            &
                                    block_begin       = field_block_begin,          &
                                    block_end         = field_block_end,            &
                                    BlockInBlockFound = BlockFound,                 &
                                    FirstLine         = FirstLine,                  &
                                    LastLine          = LastLine,                   &
                                    STAT              = STAT_CALL)

cd1 :   if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :       if (.not. BlockFound) then                                                  
                stop 'ReadFields - ModuleSWAN - ERR10'
            end if cd2
        else cd1
            stop 'ReadFields - ModuleSWAN - ERR20'
        end if cd1


        Me%NumberFields = LastLine - FirstLine - 1

        if (Me%NumberFields/=Me%NumberDates) then
            write(*,*) "The number of dates can not be different from the number of files"
            stop 'ReadFields - ModuleSWAN - ERR30'
        endif

        allocate(Me%FilesName(Me%NumberFields))  

        allocate(Me%Fields(Me%NumberProps, Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))


d3:     do l= 1, Me%NumberFields

            Line = FirstLine + l

            call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadFields - ModuleSWAN - ERR40'

            Me%FilesName(l) = AuxChar

        enddo d3

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFields - ModuleSWAN - ERR50'

    end subroutine ReadFields 

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine ConstructGlobalOutput 

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        

        nullify(Me%OutPut%OutTime)

        Me%OutPut%OutTime => Me%InputDates
        Me%OutPut%Number  =  Me%NumberDates

        Me%OutPut%NextOutPut = 1
                                                                               
    end subroutine ConstructGlobalOutput

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadFieldFromFile(l) 

        !Arguments-------------------------------------------------------------
        integer                                     :: l
        !Local----------------------------------------------------------------
        real                                        :: x, y
        integer                                     :: i, j, p, STAT_CALL
        !Begin-----------------------------------------------------------------
        

        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldFromFile - ModuleTecnoceanAscii - ERR10'

        open(Unit   = Me%Unit,                                                          &
             File   = Me%FilesName(l),                                                  &
             Form   = 'FORMATTED',                                                      &
             STATUS = 'OLD',                                                            &
             Action = 'READ',                                                           &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldFromFile - ModuleTecnoceanAscii - ERR20'

        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
    
            read(Me%Unit,*) x, y, Me%PropVector

            do p = 1, Me%NumberProps
 
                Me%Fields(p, i, j) = Me%PropVector(p)

            enddo

        enddo
        enddo


        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldFromFile - ModuleTecnoceanAscii - ERR30'


    end subroutine ReadFieldFromFile
    

    !----------------------------------------------------------------------------



    
    !------------------------------------------------------------------------

    
    subroutine ConstructBathymetry
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !Me%GridFileName="NewGrid.dat_.new"
        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleSWAN - ERR10'

        call GetHorizontalGridSize(Me%ObjHorizontalGrid, Me%Size, Me%WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleSWAN - ERR20'

        call ConstructGridData(Me%ObjBathymetry, Me%ObjHorizontalGrid, FileName = Me%GridFileName,&
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleSWAN - ERR30'


        call ConstructHorizontalMap(Me%ObjHorizontalMap, Me%ObjBathymetry, Me%ObjHorizontalGrid, &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleSWAN - ERR40'

        call GetWaterPoints2D   (Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructBathymetry - ModuleSWAN - ERR50'


    end subroutine ConstructBathymetry

    
    !------------------------------------------------------------------------

    

    
    !------------------------------------------------------------------------
    subroutine OutputFields(p)

        !Arguments-------------------------------------------------------------
        integer                                         :: p
        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        real, dimension(:,:), pointer                   :: Aux2D, Aux2DX, Aux2DY, Aux2DXY
        real                                            :: Angle
        integer                                         :: STAT_CALL, i, ii, jj, k
        character(len=StringLength)                     :: FieldName
        !Begin-----------------------------------------------------------------

        if (p==1) i = Me%OutPut%NextOutPut

        if (i==1 .and. p==1) then
            allocate(Aux2D (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Aux2DX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Aux2DY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        endif

!i0:     if (i<=Me%Output%Number) then

!dw1:    do while (Field%NextTime >= Me%Output%OutTime(i)) 
       
            !call InterpolateField(Field, i)        
            
i1:         if (p==1) then

        !           Dados para escriver uma soa vez cada date:
                call ExtractDate   (Me%Output%OutTime(i),                               &
                                    AuxTime(1), AuxTime(2), AuxTime(3),                 &
                                    AuxTime(4), AuxTime(5), AuxTime(6))

                TimePtr => AuxTime

                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleSWAN - ERR10'


                call HDF5WriteData  (Me%ObjHDF5, "/Time",                               &
                                     "Time", "YYYY/MM/DD HH:MM:SS",                     &
                                     Array1D = TimePtr,                                 &
                                     OutputNumber = i, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleSWAN - ERR20'

            endif i1

   
            call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,            &
                               Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleSWAN - ERR60'

            Aux2D(:,:) = Me%Fields(p,:,:)

            do jj = Me%WorkSize%JLB, Me%WorkSize%JUB
            do ii = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ConvertProp(p) == CartToNauticalDegrees_) then
                    if (Aux2D(ii,jj) /= -999.)  Aux2D(ii,jj) = 270. - Aux2D(ii,jj)
                endif
            enddo
            enddo

            call HDF5WriteData(Me%ObjHDF5,                                              &
                               "/Results/"//trim(Me%PropsName(p)),                      &
                               trim(Me%PropsName(p)),                                   &
                               trim(Me%PropsUnits(p)),                                  &
                               Array2D      = Aux2D,                                    &
                               OutputNumber = i,                                        &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleSWAN - ERR70'

ip:         if (trim(Me%PropsName(p)) == GetPropertyName(MeanWaveDirection_)) then

                Aux2DX(:,:) = 0.
                Aux2DY(:,:) = 0.

d2:             do k=1,2

                    do jj= Me%WorkSize%JLB,Me%WorkSize%JUB  
                    do ii= Me%WorkSize%ILB,Me%WorkSize%IUB  

                        if (Me%WaterPoints2D(ii,jj) == WaterPoint) then

                            Angle = Aux2D(ii,jj)

                            if (k==1) then
                                Aux2DX(ii,jj) = cos(Angle * Pi / 180.) 
                            else
                                Aux2DY(ii,jj) = sin(Angle * Pi / 180.) 
                            endif

                        endif

                    enddo
                    enddo

           
                    if (k==1) then
                        FieldName = trim(Me%PropsName(p))//'_x'
                        Aux2DXY => Aux2DX
                    else
                        FieldName = trim(Me%PropsName(p))//'_y'
                        Aux2DXY => Aux2DY
                    endif

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//FieldName,                          &
                                       FieldName,                                       &
                                       '-',                                             &
                                       Array2D      = Aux2DXY,                          &
                                       OutputNumber = i,                                &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleSWAN - ERR80'

                enddo d2

            endif ip


            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleSWAN - ERR90'

            if (p == Me%NumberProps) then

                Me%OutPut%NextOutPut = i + 1

            endif

!        enddo dw1

!        endif i0


        if (i==Me%NumberDates .and. p==Me%NumberProps) then
            deallocate(Aux2D)
            deallocate(Aux2DX)
            deallocate(Aux2DY)
            nullify(Aux2D, Aux2DX, Aux2DY, Aux2DXY)
        endif

    end subroutine OutputFields


    !----------------------------------------------------------------------

   !------------------------------------------------------------------------
    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        real,    dimension(:,:), pointer            :: Bathymetry
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        call GetGridData        (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleSWAN - ERR10'

      
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleSWAN - ERR30'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleSWAN - ERR40'

            
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",                   &
                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleSWAN - ERR50'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleSWAN - ERR60'            
   
       
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleSWAN - ERR70'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",                &
                              Array2D = Me%WaterPoints2D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleSWAN - ERR80'


        call UnGetGridData      (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleSWAN - ERR90'

        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleSWAN - ERR110'

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

   
    !--------------------------------------------------------------------------

    
    subroutine KillSWAN
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call UnGetHorizontalMap (Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleSWAN - ERR10'

        call KillHorizontalMap(Me%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleSWAN - ERR20'

        call KillGridData(Me%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleSWAN - ERR30'


        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleSWAN - ERR40'


        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillSWAN - ModuleSWAN - ERR60'

        deallocate(Me%InputDates )
        nullify   (Me%InputDates )
        if (Me%OutPut%ON) then
            deallocate(Me%OutPut%OutTime)
        endif
        nullify   (Me%OutPut%OutTime)



        deallocate(Me%PropsName )
        nullify   (Me%PropsName )

        deallocate(Me%PropsUnits)
        nullify   (Me%PropsUnits)

        deallocate(Me%ConvertProp)
        nullify   (Me%ConvertProp)

        deallocate(Me%FilesName )
        nullify   (Me%FilesName )

        deallocate(Me%Fields    )
        nullify   (Me%Fields    )

        deallocate(Me%PropVector)
        nullify   (Me%PropVector)

        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillSWAN

    !--------------------------------------------------------------------------
 
end module ModuleSWAN
