!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : GFSasciiWind
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group & MeteoGalicia
! DATE          : September 2003
! REVISION      : Pablo Carracedo - v4.0
! DESCRIPTION   : Module to convert GFSasciiWind files into HDF5 format.
!                 For reading into Mohid module HydrodynamicFile
!
!------------------------------------------------------------------------------


Module ModuleGFSasciiWind

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
    public  :: ConvertGFSasciiWind
    private ::      ReadGlobalOptions
    private ::      ReadDates
    private ::      ReadProperties
    private ::      ReadFields
    private ::      ConstructGlobalOutput
    private ::      ConstructGrid
    private ::      Open_HDF5_OutPut_File
    private ::          OutputFields
    private ::      KillGFSasciiWind

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

    
    private :: T_GFSasciiWind
    type       T_GFSasciiWind
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
        integer                                 :: NumberOfFiles        = FillValueInt
        integer                                 :: NumberProps          = FillValueInt
        integer                                 :: NumberUnits          = FillValueInt

        integer                                 :: DirectionReferential = FillValueInt

        character(len=StringLength), dimension(:), pointer :: PropsName
        character(len=StringLength), dimension(:), pointer :: PropsUnits
        character(len=PathLength),   dimension(:), pointer :: FilesName
        integer,                     dimension(:), pointer :: ConvertProp

        type(T_Time),     dimension(:), pointer :: StartFieldsDate
        type(T_Time),     dimension(:), pointer :: FieldsInstant
        real,     dimension(:,:,:,:),   pointer :: Fields 

        type(T_Time)                            :: StartTime, EndTime


        integer                                 :: Clientumber

        integer                                 :: DatesPerFile, DatesNumber
        integer                                 :: ReadType
        integer, dimension(:,:  ),  pointer     :: WaterPoints2D
        type(T_OutPut)                          :: OutPut
        type(T_Size2D)                          :: WorkSize, Size
    end type  T_GFSasciiWind

    type(T_GFSasciiWind), pointer              :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertGFSasciiWind(EnterDataID, ClientNumber, STAT)

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

        call ConstructGrid

        call ReadProperties

        call ReadConversion

        call ReadDates

        call ReadFields

        call ConstructGlobalOutput

        call Open_HDF5_OutPut_File


d1:     do l = 1, Me%NumberOfFiles

            call ReadFieldFromFile(l)

d2:         do p=1, Me%NumberProps

                call OutputFields  

            enddo d2


        enddo d1

  
        call KillGFSasciiWind


        STAT = SUCCESS_


    end subroutine ConvertGFSasciiWind

    !--------------------------------------------------------------------------


    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------

        call GetData(Me%GridFileName,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'INPUT_GRID_FILENAME',                              &
                     ClientModule = 'GFSasciiWind',                                     &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleGFSasciiWind - ERR10'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleGFSasciiWind - ERR20'

        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'GFSasciiWind',                                     &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleGFSasciiWind - ERR30'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleGFSasciiWind - ERR40'

       
        call GetData(Me%DatesPerFile,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'DATES_PER_FILE',                                   &
                     default      = 24,                                                 & 
                     ClientModule = 'GFSasciiWind',                                     &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleGFSasciiWind - ERR50'
 
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
                stop 'ReadDates - ModuleGFSasciiWind - ERR10'
            end if cd2
        else cd1
            stop 'ReadDates - ModuleGFSasciiWind - ERR20'
        end if cd1


        Me%NumberDates = LastLine - FirstLine - 1

        allocate(Me%StartFieldsDate(Me%NumberDates))  


d3:     do l= 1, Me%NumberDates

            Line = FirstLine + l

            call GetData(Aux6, EnterDataID = Me%ObjEnterData, flag = iflag,             &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadDates - ModuleGFSasciiWind - ERR30'

            call SetDate(AuxTime          ,                                             &
                        Year    = Aux6(1) ,                                             &
                        Month   = Aux6(2) ,                                             & 
                        Day     = Aux6(3) ,                                             &
                        Hour    = Aux6(4) ,                                             &
                        Minute  = Aux6(5) ,                                             &
                        Second  = Aux6(6) )  

            Me%StartFieldsDate(l) = AuxTime

        enddo d3

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDates - ModuleGFSasciiWind - ERR40'


    end subroutine ReadDates 

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadProperties

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
!        integer                                     :: FirstLine, LastLine
!        character(len=StringLength)                 :: AuxChar
!        logical                                     :: BlockFound
        !Begin-----------------------------------------------------------------


!        call ExtractBlockFromBlock (Me%ObjEnterData,                                    &
!                                    ClientNumber      = Me%ClientNumber,                &
!                                    block_begin       = prop_block_begin,               &
!                                    block_end         = prop_block_end,                 &
!                                    BlockInBlockFound = BlockFound,                     &
!                                    FirstLine         = FirstLine,                      &
!                                    LastLine          = LastLine,                       &
!                                    STAT              = STAT_CALL)

!cd1 :   if      (STAT_CALL .EQ. SUCCESS_     ) then    
!cd2 :       if (.not. BlockFound) then                                                  
!                stop 'ReadProperties - ModuleGFSasciiWind - ERR10'
!            end if cd2
!        else cd1
!            stop 'ReadProperties - ModuleGFSasciiWind - ERR20'
!        end if cd1


        !Me%NumberProps = LastLine - FirstLine - 1
        Me%NumberProps = 3

        allocate(Me%PropsName (Me%NumberProps))  
        allocate(Me%PropsUnits(Me%NumberProps))

        Me%PropsName(1) = GetPropertyName(WindVelocityX_)
        Me%PropsName(2) = GetPropertyName(WindVelocityY_)
        Me%PropsName(3) = GetPropertyName(WindModulus_  )

        Me%PropsUnits(1:3) = 'm/s' 


!d1:     do l= 1, Me%NumberProps

!            line = FirstLine + l

 !           call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
!                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

!            if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleGFSasciiWind - ERR30'


!            Me%PropsName(l) = AuxChar

!            if (.not. CheckPropertyName (AuxChar)) then
!                write(*,*) 'The name ',trim(AuxChar),' is not valid name for the MOHID system'
!            endif

!        enddo d1

!        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleGFSasciiWind - ERR40'


!        call ExtractBlockFromBlock (Me%ObjEnterData,                                    &
!                                    ClientNumber      = Me%ClientNumber,                &
!                                    block_begin       = units_block_begin,              &
!                                    block_end         = units_block_end,                &
!                                    BlockInBlockFound = BlockFound,                     &
!                                    FirstLine         = FirstLine,                      &
!                                    LastLine          = LastLine,                       &
!                                    STAT              = STAT_CALL)

!cd3 :   if      (STAT_CALL .EQ. SUCCESS_     ) then    
!cd4 :       if (.not. BlockFound) then                                                  
!                stop 'ReadProperties - ModuleGFSasciiWind - ERR50'
!            end if cd4
!        else cd3
!            stop 'ReadProperties - ModuleGFSasciiWind - ERR60'
!        end if cd3


!        Me%NumberUnits = LastLine - FirstLine - 1

!        allocate(Me%PropsUnits(Me%NumberUnits))

!d2:     do l= 1, Me%NumberUnits

!            line = FirstLine + l

!            call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
!                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

!            if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleGFSasciiWind - ERR70'

!            Me%PropsUnits(l) = AuxChar

!        enddo d2

!        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleGFSasciiWind - ERR80'


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
                stop 'ReadConversion - ModuleGFSasciiWind - ERR10'
            end if cd1

            call GetData(PropName,                                                      &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NAME',                                         &
                         ClientModule = 'ConvertToHDF5',                                &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadConversion - ModuleGFSasciiWind - ERR20'

            PropNameWrong = .true.
            do j = 1, Me%NumberProps
                if (trim(Me%PropsName(j))==trim(PropName)) then
                    PropNameWrong = .false. 
                    exit
                endif
            enddo
            
            if (PropNameWrong) stop 'ReadConversion - ModuleGFSasciiWind - ERR30'

            call GetData(ConversionType,                                                &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'CONVERSION',                                   &
                         ClientModule = 'ConvertToHDF5',                                &
                         default      = NoConversion_,                                  &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadConversion - ModuleGFSasciiWind - ERR40'

            if (ConversionType /= NoConversion_ .and. ConversionType /= CartToNauticalDegrees_) then
                stop 'ReadConversion - ModuleGFSasciiWind - ERR50'
            endif

            Me%ConvertProp(j) = ConversionType

        enddo d1

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadConversion - ModuleGFSasciiWind - ERR60'


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
                stop 'ReadFields - ModuleGFSasciiWind - ERR10'
            end if cd2
        else cd1
            stop 'ReadFields - ModuleGFSasciiWind - ERR20'
        end if cd1


        Me%NumberOfFiles = LastLine - FirstLine - 1

        if (Me%NumberOfFiles/=Me%NumberDates) then
            write(*,*) "The number of dates can not be different from the number of files"
            stop 'ReadFields - ModuleGFSasciiWind - ERR30'
        endif

        allocate(Me%FilesName(Me%NumberOfFiles))  

        allocate(Me%Fields       (Me%NumberProps, Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB,Me%DatesPerFile))
        allocate(Me%FieldsInstant(Me%DatesPerFile))

d3:     do l= 1, Me%NumberOfFiles

            Line = FirstLine + l

            call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadFields - ModuleGFSasciiWind - ERR40'

            Me%FilesName(l) = AuxChar

        enddo d3

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFields - ModuleGFSasciiWind - ERR50'

    end subroutine ReadFields 

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine ConstructGlobalOutput 

        !Local-----------------------------------------------------------------
        integer                           :: STAT_CALL
        !Begin-----------------------------------------------------------------
        

        nullify(Me%OutPut%OutTime)

        Me%StartTime = Me%StartFieldsDate(1)

        Me%EndTime   = Me%StartFieldsDate(Me%NumberDates) + 3600 * 24 * 15

        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime = Me%StartTime,                                  &
                           EndTime     = Me%EndTime,                                    &
                           keyword     = 'OUTPUT_TIME',                                 &
                           SearchType  = FromBlock,                                     &
                           OutPutsTime = Me%OutPut%OutTime,                             &
                           OutPutsOn   = Me%OutPut%ON,                                  &
                           STAT        = STAT_CALL)                                     
                                                                                        
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructGlobalOutput - ModuleGFSasciiWind - ERR10'              
                                                                            
        if (Me%OutPut%ON) then
            Me%OutPut%NextOutPut = 1
        else
            write(*,*)'Neeed to define Keyword OUTPUT_TIME'
            stop 'ConstructGlobalOutput - ModuleGFSasciiWind - ERR20'
        endif 

                                                                               
    end subroutine ConstructGlobalOutput

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadFieldFromFile(l) 

        !Arguments-------------------------------------------------------------
        integer                                     :: l
        !Local----------------------------------------------------------------
        type (T_Time)                               :: AuxTime 
        real    , dimension(6)                      :: Aux6
!        real                                        :: x, y
        integer                                     :: i, j, STAT_CALL, k, it
        character(len=17)                           :: Aux17
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

        k = 0

        do it=1, Me%DatesPerFile

            read(Me%Unit,'(A17)',END=10) Aux17

            !Year
            read(Aux17( 1: 4),'(f4.0)') Aux6(1)
            !Month
            read(Aux17( 5: 6),'(f2.0)') Aux6(2)
            !Day
            read(Aux17( 7: 8),'(f2.0)') Aux6(3)
            !Hours
            read(Aux17(12:13),'(f2.0)') Aux6(4)
            !Hours
            read(Aux17(14:15),'(f2.0)') Aux6(5)
            !Seconds
            read(Aux17(16:17),'(f2.0)') Aux6(6)


            call SetDate(AuxTime          ,                                             &
                        Year    = Aux6(1) ,                                             &
                        Month   = Aux6(2) ,                                             & 
                        Day     = Aux6(3) ,                                             &
                        Hour    = Aux6(4) ,                                             &
                        Minute  = Aux6(5) ,                                             &
                        Second  = Aux6(6) )  

            Me%FieldsInstant(it) = AuxTime

            do i=Me%WorkSize%ILB, Me%WorkSize%IUB
                read(Me%Unit,'(2000f7.2)') (Me%Fields(1, i, j,it),j=Me%WorkSize%JLB,Me%WorkSize%JUB)
            enddo


            do i=Me%WorkSize%ILB, Me%WorkSize%IUB
                read(Me%Unit,'(2000f7.2)') (Me%Fields(2, i, j,it),j=Me%WorkSize%JLB,Me%WorkSize%JUB)
            enddo

            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                Me%Fields(3, i, j,it) = sqrt(Me%Fields(1, i, j,it)**2 + Me%Fields(2, i, j,it)**2)

            enddo
            enddo

            k = k + 1

        enddo

10      continue

        Me%DatesNumber = k

        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldFromFile - ModuleTecnoceanAscii - ERR30'


    end subroutine ReadFieldFromFile
    

    !----------------------------------------------------------------------------



    
    !------------------------------------------------------------------------

    
    subroutine ConstructGrid
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !Me%GridFileName="NewGrid.dat_.new"
        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleGFSasciiWind - ERR10'

        call GetHorizontalGridSize(Me%ObjHorizontalGrid, Me%Size, Me%WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleGFSasciiWind - ERR20'

        !call ConstructGridData(Me%ObjBathymetry, Me%ObjHorizontalGrid, FileName = Me%GridFileName,&
        !                       STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleGFSasciiWind - ERR30'


        !call ConstructHorizontalMap(Me%ObjHorizontalMap, Me%ObjBathymetry, Me%ObjHorizontalGrid, &
        !                       STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleGFSasciiWind - ERR40'

        !call GetWaterPoints2D   (Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'ConstructGrid - ModuleGFSasciiWind - ERR50'


    end subroutine ConstructGrid

    
    !------------------------------------------------------------------------

    

    
    !------------------------------------------------------------------------
    subroutine OutputFields

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        real, dimension(:,:), pointer                   :: Aux2D
!        real                                            :: Angle
        integer                                         :: STAT_CALL, i, it, p
        !Begin-----------------------------------------------------------------

        i = Me%OutPut%NextOutPut

        allocate(Aux2D (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        do it =1, Me%DatesNumber 
dw1:    do while (Me%FieldsInstant(it) >= Me%Output%OutTime(i)) 
       

    !           Dados para escriver uma soa vez cada date:
            call ExtractDate   (Me%Output%OutTime(i),                               &
                                AuxTime(1), AuxTime(2), AuxTime(3),                 &
                                AuxTime(4), AuxTime(5), AuxTime(6))

            TimePtr => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleGFSasciiWind - ERR10'


            call HDF5WriteData  (Me%ObjHDF5, "/Time",                               &
                                 "Time", "YYYY/MM/DD HH:MM:SS",                     &
                                 Array1D = TimePtr,                                 &
                                 OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleGFSasciiWind - ERR20'

   
            call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,            &
                               Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleGFSasciiWind - ERR60'

            p = 1

            Aux2D(:,:) = Me%Fields(p,:,:,it)

            call HDF5WriteData(Me%ObjHDF5,                                              &
                               "/Results/"//trim(Me%PropsName(p)),                      &
                               trim(Me%PropsName(p)),                                   &
                               trim(Me%PropsUnits(p)),                                  &
                               Array2D      = Aux2D,                                    &
                               OutputNumber = i,                                        &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleGFSasciiWind - ERR70'

            p = 2

            Aux2D(:,:) = Me%Fields(p,:,:,it)

            call HDF5WriteData(Me%ObjHDF5,                                              &
                               "/Results/"//trim(Me%PropsName(p)),                      &
                               trim(Me%PropsName(p)),                                   &
                               trim(Me%PropsUnits(p)),                                  &
                               Array2D      = Aux2D,                                    &
                               OutputNumber = i,                                        &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleGFSasciiWind - ERR70'

            p = 3

            Aux2D(:,:) = Me%Fields(p,:,:,it) 

            call HDF5WriteData(Me%ObjHDF5,                                              &
                               "/Results/"//trim(Me%PropsName(p)),                      &
                               trim(Me%PropsName(p)),                                   &
                               trim(Me%PropsUnits(p)),                                  &
                               Array2D      = Aux2D,                                    &
                               OutputNumber = i,                                        &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleGFSasciiWind - ERR70'


            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleGFSasciiWind - ERR90'

            i                    = i + 1
            Me%OutPut%NextOutPut = i


        enddo dw1
        enddo 

!        endif i0


        deallocate(Aux2D)

    end subroutine OutputFields


    !----------------------------------------------------------------------

   !------------------------------------------------------------------------
    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        !real,    dimension(:,:), pointer            :: Bathymetry
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !call GetGridData        (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleGFSasciiWind - ERR10'

      
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleGFSasciiWind - ERR30'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleGFSasciiWind - ERR40'

            
        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleGFSasciiWind - ERR60'            
   
       
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleGFSasciiWind - ERR70'            

        allocate(Me%WaterPoints2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%WorkSize%JUB))

        Me%WaterPoints2D(:,:) = 1

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",                &
                              Array2D = Me%WaterPoints2D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleGFSasciiWind - ERR80'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",                   &
                              Array2D =  Me%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleGFSasciiWind - ERR50'            

        deallocate(Me%WaterPoints2D)


        !call UnGetGridData      (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleGFSasciiWind - ERR90'

        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleGFSasciiWind - ERR110'

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

   
    !--------------------------------------------------------------------------

    
    subroutine KillGFSasciiWind
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        !call UnGetHorizontalMap (Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'KillGFSasciiWind - ModuleGFSasciiWind - ERR10'

        !call KillHorizontalMap(Me%ObjHorizontalMap, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'KillGFSasciiWind - ModuleGFSasciiWind - ERR20'

        !call KillGridData(Me%ObjBathymetry, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'KillGFSasciiWind - ModuleGFSasciiWind - ERR30'


        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillGFSasciiWind - ModuleGFSasciiWind - ERR40'


        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillGFSasciiWind - ModuleGFSasciiWind - ERR60'

        deallocate(Me%StartFieldsDate )
        nullify   (Me%StartFieldsDate )
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

        deallocate(Me%FieldsInstant)
        nullify   (Me%FieldsInstant)

        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillGFSasciiWind

    !--------------------------------------------------------------------------
 
end module ModuleGFSasciiWind
