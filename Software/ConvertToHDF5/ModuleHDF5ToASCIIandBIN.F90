!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : HDF5ToASCIIandBIN
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group & MeteoGalicia
! DATE          : September 2003
! REVISION      : Pablo Carracedo - v4.0
! DESCRIPTION   : Module to convert HDF5ToASCIIandBIN files into HDF5 format.
!                 For reading into Mohid module HydrodynamicFile
!
!------------------------------------------------------------------------------


Module ModuleHDF5ToASCIIandBIN

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
    public  :: ConvertHDF5ToASCIIandBIN
    private ::      ReadGlobalOptions
    private ::      ReadProperties
    private ::      ConstructGrid
    private ::      Open_HDF5_InPut_File
    private ::      OutputFields
    private ::          OutputMohidBin
    private ::          OutputSwanAscii
    private ::      KillHDF5ToASCIIandBIN

    !Parameters----------------------------------------------------------------

    character(LEN = StringLength), parameter    :: prop_block_begin   = '<<beginproperty>>'
    character(LEN = StringLength), parameter    :: prop_block_end     = '<<endproperty>>'

    character(LEN = StringLength), parameter    :: date_block_begin   = '<<begindate>>'
    character(LEN = StringLength), parameter    :: date_block_end     = '<<enddate>>'

    character(LEN = StringLength), parameter    :: convert_block_begin   = '<<beginconvert>>'
    character(LEN = StringLength), parameter    :: convert_block_end     = '<<endconvert>>'

    integer                      , parameter    :: Swan1_ = 1, Swan2_ = 2, Swan3_ = 3,  & 
                                                   Swan4_ = 4, Swan5_ = 5, Swan6_ = 6,  &
                                                   SwanTable_ = 7, Mohid_ = 8, WW3_ = 9

    integer                      , parameter    :: GeoXY_ = 1, CartesianXY_ = 2


    integer                      , parameter    :: NoConversion_ = 0, NauticalToCartDegrees_ = 1
    !Types---------------------------------------------------------------------
    
    
    type       T_OutPut                                 
         type (T_Time), pointer, dimension(:)               :: OutTime
         integer                                            :: Number, Next
         logical                                            :: ON
    end type T_OutPut                                   


    
    private :: T_HDF5ToASCIIandBIN
    type       T_HDF5ToASCIIandBIN
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
        character(len=PathLength)               :: InPutFileName

        integer                                 :: RequestDates         = FillValueInt
        integer                                 :: TotalDates           = FillValueInt
        integer                                 :: NumberFields         = FillValueInt
        integer                                 :: NumberProps          = FillValueInt
        integer                                 :: NumberUnits          = FillValueInt

        integer                                 :: DirectionReferential = FillValueInt

        type (T_OutPut)                         :: Output

        character(len=StringLength), dimension(:), pointer :: PropsName
        integer,                     dimension(:), pointer :: ConvertProp
        
        character(len=PathLength)               :: OutPutPath

        real,     dimension(:),         pointer :: PropVector
        real,     dimension(:,:),       pointer :: X2D, Y2D

        integer                                 :: OutPutOption

        integer                                 :: Clientumber

        logical                                 :: WriteVelModulus = .false., WriteWindModulus = .false.
        real                                    :: FillValue
        integer                                 :: ReadType
        integer, dimension(:,:  ),  pointer     :: WaterPoints2D

        integer                                 :: OutputXY = FillValueInt

        type(T_Size2D)                          :: WorkSize, Size
    end type  T_HDF5ToASCIIandBIN

    type(T_HDF5ToASCIIandBIN), pointer              :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertHDF5ToASCIIandBIN(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !Local-------------------------------------------------------------------
        real, dimension(:  ), pointer                   :: AuxTime
        real, dimension(:,:), pointer                   :: Aux2DNext, Aux2DPrev, Aux2D
        integer                                         :: l, p, STAT_CALL
        type (T_Time)                                   :: NextTime, PrevTime
        logical                                         :: FirstProp
        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        Me%ClientNumber = ClientNumber

        call ReadGlobalOptions

        if (Me%OutPutOption /= WW3_) then
            call ConstructGlobalOutput
        endif

        call Open_HDF5_InPut_File

        call ConstructGrid

        call ReadProperties

        call ReadConversion

        !call ReadDates

        allocate(Aux2D     (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))  
        allocate(Aux2DNext (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate(Aux2DPrev (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        allocate(AuxTime(1:6))


        if (Me%OutPutOption == WW3_) then
            call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'OutputWW3ASCII - ModuleHDF5ToASCIIandBIN - ERR10'
 

            if (trim(Me%OutPutPath) /= '*') then
                Me%OutPutPath =  'WW3_Input.txt'
            else
                Me%OutPutPath = trim(Me%OutPutPath)//'/'//'WW3_Input.txt'
            endif

            open(Unit   = Me%Unit,                                                          &
                 File   = trim(Me%OutPutPath),                                              &
                 Form   = 'FORMATTED',                                                      &
                 STATUS = 'UNKNOWN',                                                        &
                 Action = 'WRITE',                                                          &
                 IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'OutputWW3ASCII - ModuleHDF5ToASCIIandBIN - ERR20'

        endif

i11:    if (Me%OutPutOption == WW3_) then


d11:        do l = 1, Me%TotalDates

                call HDF5SetLimits(Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputWW3ASCII - ModuleHDF5ToASCIIandBIN - ERR30'


                call HDF5ReadData(Me%ObjHDF5,                                               &
                                   "/Time",                                                 &
                                   "Time",                                                  &
                                   Array1D      = AuxTime,                                  &
                                   OutputNumber = l,                                        &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputWW3ASCII - ModuleHDF5ToASCIIandBIN - ERR40'


    d12:         do p=1, Me%NumberProps


                        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,            &
                                           Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                                           STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)stop 'OutputWW3ASCII - ModuleHDF5ToASCIIandBIN - ERR50'


                        call HDF5ReadData(Me%ObjHDF5,                                                   &
                                           "/Results/"//trim(Me%PropsName(p)),                          &
                                           trim(Me%PropsName(p)),                                       &
                                           Array2D      = Aux2D,                                        &
                                           OutputNumber = l,                                            &
                                           STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)stop 'OutputWW3ASCII - ModuleHDF5ToASCIIandBIN - ERR60'

                        FirstProp = .false. 

                        if (p ==1)             FirstProp = .true.

                        call OutputWW3ASCII  (AuxTime, Aux2D, FirstProp) 

                 enddo d12

            enddo d11

        else 


    d1:     do l = 1, Me%TotalDates - 1

                call HDF5SetLimits(Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ConvertHDF5ToASCIIandBIN - ModuleHDF5ToASCIIandBIN - ERR10'

                if (l==1) then

                    call HDF5ReadData(Me%ObjHDF5,                                               &
                                       "/Time",                                                 &
                                       "Time",                                                  &
                                       Array1D      = AuxTime,                                  &
                                       OutputNumber = l,                                        &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'ConvertHDF5ToASCIIandBIN - ModuleHDF5ToASCIIandBIN - ERR20'


                    call SetDate(PrevTime, AuxTime(1), AuxTime(2), AuxTime(3), AuxTime(4), AuxTime(5), AuxTime(6))

                else

                    PrevTime = NextTime

                endif


                call HDF5ReadData(Me%ObjHDF5,                                               &
                                   "/Time",                                                 &
                                   "Time",                                                  &
                                   Array1D      = AuxTime,                                  &
                                   OutputNumber = l+1,                                       &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5ToASCIIandBIN - ModuleHDF5ToASCIIandBIN - ERR30'


                call SetDate(NextTime, AuxTime(1), AuxTime(2), AuxTime(3), AuxTime(4), AuxTime(5), AuxTime(6))



    d2:         do p=1, Me%NumberProps

                        call OutputFields     (Aux2D, Aux2DNext, Aux2DPrev, NextTime, PrevTime, l, p)
                

                enddo d2

                if (Me%Output%Next > Me%Output%Number) exit

            

            enddo d1

        endif i11
 
        deallocate(Aux2D)
        nullify(Aux2D)

        deallocate(Aux2DNext)
        nullify(Aux2DNext)

        deallocate(Aux2DPrev)
        nullify(Aux2DPrev)

        deallocate(AuxTime)
        nullify(AuxTime)

        call KillHDF5ToASCIIandBIN


        STAT = SUCCESS_


    end subroutine ConvertHDF5ToASCIIandBIN

    !------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ConstructGlobalOutput 

        !Local-----------------------------------------------------------------
        type (T_Time)                               :: StartTime, EndTime
        integer                                     :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------
        

        nullify(Me%OutPut%OutTime)

        call GetData(StartTime,                                                         &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'START',                                            &
                     ClientModule = 'ConvertToHDF5',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalOutput - ModuleCowamaAsciiWind - ERR10'

        call GetData(EndTime,                                                           &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'END',                                              &
                     ClientModule = 'ConvertToHDF5',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalOutput - ModuleCowamaAsciiWind - ERR20'


        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime      = StartTime,                                &
                           EndTime          = EndTime,                                  &
                           keyword          = 'OUTPUT_TIME',                            &
                           SearchType       = FromBlock,                                &
                           OutPutsTime      = Me%OutPut%OutTime,                        &
                           OutPutsOn        = Me%OutPut%ON,                             &
                           OutPutsNumber    = Me%OutPut%Number,                         &
                           STAT             = STAT_CALL)                                
                                                                                        
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructGlobalOutput - ModuleCowamaAsciiWind - ERR30'              

        Me%OutPut%Next = 1

    end subroutine ConstructGlobalOutput

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
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleHDF5ToASCIIandBIN - ERR10'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleHDF5ToASCIIandBIN - ERR20'

        call GetData(Me%InPutFileName,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'INPUTFILENAME',                                    &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleHDF5ToASCIIandBIN - ERR30'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleHDF5ToASCIIandBIN - ERR40'

       
        call GetData(Me%FillValue,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FILL_VALUE',                                       &
                     default      = -99.999900,                                         &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleHDF5ToASCIIandBIN - ERR50'

        call GetData(Me%OutPutPath,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'PATH',                                             &
                     default      = '*',                                                &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleHDF5ToASCIIandBIN - ERR60'

        call GetData(Me%OutPutOption,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUT_OPTION',                                    &
                     default      = Swan1_,                                             &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleHDF5ToASCIIandBIN - ERR70'

        if (Me%OutPutOption == SwanTable_) then

            call GetData(Me%OutPutXY,                                                       &
                         Me%ObjEnterData, iflag,                                            &
                         SearchType   = FromBlock,                                          &
                         keyword      = 'OUTPUT_XY',                                        &
                         default      = GeoXY_,                                             &
                         ClientModule = 'HDF5ToASCIIandBIN',                                &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleHDF5ToASCIIandBIN - ERR80'


            if (Me%OutPutXY /= GeoXY_ .and. Me%OutPutXY /= CartesianXY_) then 

                write(*,*) 'The only available options for the OUTPUT_XY keyword are'
                write(*,*) 'Geographic Long Lat - OUTPUT_XY : 1'
                write(*,*) 'Cartesian  X    Y   - OUTPUT_XY : 2'
                stop 'ReadGlobalOptions - ModuleHDF5ToASCIIandBIN - ERR90'

            endif

        endif
 
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
                stop 'ReadProperties - ModuleHDF5ToASCIIandBIN - ERR10'
            end if cd2
        else cd1
            stop 'ReadProperties - ModuleHDF5ToASCIIandBIN - ERR20'
        end if cd1


        Me%NumberProps = LastLine - FirstLine - 1

        allocate(Me%PropsName (Me%NumberProps))  
  
        allocate(Me%PropVector(Me%NumberProps))  


d1:     do l= 1, Me%NumberProps

            line = FirstLine + l

            call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleHDF5ToASCIIandBIN - ERR30'


            Me%PropsName(l) = AuxChar

            if (.not. CheckPropertyName (AuxChar)) then
                write(*,*) 'The name ',trim(AuxChar),' is not valid name for the MOHID system'
            endif

        enddo d1

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleHDF5ToASCIIandBIN - ERR40'


    end subroutine ReadProperties 

    !--------------------------------------------------------------------------


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
                stop 'ReadConversion - ModuleHDF5ToASCIIandBIN - ERR10'
            end if cd1

            call GetData(PropName,                                                      &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NAME',                                         &
                         ClientModule = 'ConvertToHDF5',                                &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadConversion - ModuleCowamaAsciiWind - ERR20'

            PropNameWrong = .true.
            do j = 1, Me%NumberProps
                if (trim(Me%PropsName(j))==trim(PropName)) then
                    PropNameWrong = .false. 
                    exit
                endif
            enddo
            
            if (PropNameWrong) stop 'ReadConversion - ModuleCowamaAsciiWind - ERR30'

            call GetData(ConversionType,                                                &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'CONVERSION',                                   &
                         ClientModule = 'ConvertToHDF5',                                &
                         default      = NoConversion_,                                  &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadConversion - ModuleCowamaAsciiWind - ERR40'

            if (ConversionType /= NoConversion_ .and. ConversionType /= NauticalToCartDegrees_) then
                stop 'ReadConversion - ModuleCowamaAsciiWind - ERR50'
            endif

            Me%ConvertProp(j) = ConversionType

        enddo d1

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadConversion - ModuleHDF5ToASCIIandBIN - ERR60'


    end subroutine ReadConversion 

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine OutputMohidBin(PropName, TimeName, Aux2D, p) 

        !Arguments-------------------------------------------------------------
        character(len = *)                          :: PropName, TimeName
        real,   dimension(:,:), pointer             :: Aux2D
        integer                                     :: p
        !Local----------------------------------------------------------------
        character(len = PathLength)                 :: FileName
        integer                                     :: i, j, STAT_CALL
        !Begin-----------------------------------------------------------------
        
        if (p == 1) then

            call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'OutputMohidBin - ModuleHDF5ToASCIIandBIN - ERR10'
      

            FileName = "Mohid2DWaves"//'_'//trim(TimeName)//'.BIN'

            if (trim(Me%OutPutPath) /= '*') then
                FileName = trim(Me%OutPutPath)//'/'//trim(FileName)
            endif

            open(Unit   = Me%Unit,                                                          &
                 File   = FileName,                                                         &
                 Form   = 'FORMATTED',                                                      &
                 STATUS = 'UNKNOWN',                                                        &
                 Action = 'WRITE',                                                          &
                 IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'OutputMohidBin - ModuleHDF5ToASCIIandBIN - ERR20'

        endif

        write(*,*) 'write in binary property ',PropName


        write(Me%Unit) ((Aux2D(i, j),j=Me%WorkSize%JLB, Me%WorkSize%JUB),i=Me%WorkSize%ILB, Me%WorkSize%IUB)

 
        if (p == Me%NumberProps) then

            call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'OutputMohidBin - ModuleHDF5ToASCIIandBIN - ERR30'

        endif

    end subroutine OutputMohidBin
    

    !----------------------------------------------------------------------------
    


    !--------------------------------------------------------------------------

    subroutine OutputSwanASCII(PropName, TimeName, Aux2D) 

        !Arguments-------------------------------------------------------------
        character(len = *)                          :: PropName, TimeName
        real,   dimension(:,:), pointer             :: Aux2D
        !Local----------------------------------------------------------------
        character(len = PathLength)                 :: FileName
        integer                                     :: i, j, STAT_CALL
        !Begin-----------------------------------------------------------------
        

        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OutputSwanASCII - ModuleHDF5ToASCIIandBIN - ERR10'

      

        FileName = trim(PropName)//'_'//trim(TimeName)//'.txt'

        if (trim(Me%OutPutPath) /= '*') then
            FileName = trim(Me%OutPutPath)//'/'//trim(FileName)
        endif

        open(Unit   = Me%Unit,                                                          &
             File   = FileName,                                                         &
             Form   = 'FORMATTED',                                                      &
             STATUS = 'UNKNOWN',                                                        &
             Action = 'WRITE',                                                          &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OutputSwanASCII - ModuleHDF5ToASCIIandBIN - ERR20'


        if      (Me%OutPutOption == Swan1_) then

            do i=Me%WorkSize%IUB, Me%WorkSize%ILB,-1
                write(Me%Unit,100) (Aux2D(i, j),j=Me%WorkSize%JLB, Me%WorkSize%JUB)
            enddo

        else if (Me%OutPutOption == Swan2_) then

            do i=Me%WorkSize%IUB, Me%WorkSize%ILB,-1
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB

                write(Me%Unit,10) Aux2D(i, j)

            enddo
            enddo

        else if (Me%OutPutOption == Swan3_) then

            do i=Me%WorkSize%ILB, Me%WorkSize%IUB
                write(Me%Unit,100) (Aux2D(i, j),j=Me%WorkSize%JLB, Me%WorkSize%JUB)
            enddo

        else if (Me%OutPutOption == Swan4_) then


            do i=Me%WorkSize%ILB, Me%WorkSize%IUB
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB

                write(Me%Unit,10) Aux2D(i, j)

            enddo
            enddo


        else if (Me%OutPutOption == Swan5_) then

            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                write(Me%Unit,100) (Aux2D(i, j),i=Me%WorkSize%ILB, Me%WorkSize%IUB)
            enddo

        else if (Me%OutPutOption == Swan6_) then

            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                write(Me%Unit,10) Aux2D(i, j)

            enddo
            enddo

        else if (Me%OutPutOption == SwanTable_) then


            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                write(Me%Unit,30) Me%X2D, Me%Y2D, Aux2D(i, j)

            enddo
            enddo


        else

            stop 'OutputSwanASCII - ModuleHDF5ToASCIIandBIN - ERR30'

        endif

        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OutputSwanASCII - ModuleHDF5ToASCIIandBIN - ERR40'


10  format (F12.6)
30  format (3F12.6)
100 format (2000F12.6)

    end subroutine OutputSwanASCII
    

    !----------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine OutputWW3ASCII(AuxTime, Aux2D, FirstProp) 

        !Arguments-------------------------------------------------------------
        real,   dimension(:)  , pointer             :: AuxTime
        real,   dimension(:,:), pointer             :: Aux2D
        logical                                     :: FirstProp
        !Local----------------------------------------------------------------
        character(len = PathLength)                 :: AuxChar
        integer                                     :: i, j
        !Begin-----------------------------------------------------------------
        
        if (FirstProp) then
            write(AuxChar,'(I4,I2,I2,I5,I2,I2)') (int(AuxTime(i)),i=1,6)      
            do i=1,8
                if (AuxChar(i:i)==' ') AuxChar(i:i) = '0'
            enddo
            do i=12,17
                if (AuxChar(i:i)==' ') AuxChar(i:i) = '0'
            enddo
            write(Me%Unit,'(A17)')   AuxChar(1:17)
        endif


        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            if (abs(Aux2D(i, j)) >= 1000.) Aux2D(i, j) = 0.
        enddo
        enddo


        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
            write(Me%Unit,100) (Aux2D(i, j),j=Me%WorkSize%JLB, Me%WorkSize%JUB)
        enddo



100 format (2000F12.6)

    end subroutine OutputWW3ASCII    
    !------------------------------------------------------------------------

    
    subroutine ConstructGrid
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !Me%GridFileName="NewGrid.dat_.new"
        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleHDF5ToASCIIandBIN - ERR10'

        call GetHorizontalGridSize(Me%ObjHorizontalGrid, Me%Size, Me%WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleHDF5ToASCIIandBIN - ERR20'

!        call ConstructGridData(Me%ObjBathymetry, Me%ObjHorizontalGrid, FileName = Me%GridFileName,&
!                               STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleHDF5ToASCIIandBIN - ERR30'


!        call ConstructHorizontalMap(Me%ObjHorizontalMap, Me%ObjBathymetry, Me%ObjHorizontalGrid, &
!                               STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleHDF5ToASCIIandBIN - ERR40'

!        call GetWaterPoints2D   (Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'ConstructGrid - ModuleHDF5ToASCIIandBIN - ERR50'


        
        if (Me%OutPutOption == SwanTable_) then
            
            if      (Me%OutPutXY == GeoXY_      ) then 

                call GetGridLatitudeLongitude(Me%ObjHorizontalGrid,                     &
                                              GridLatitude  = Me%Y2D,                   &
                                              GridLongitude = Me%X2D,                   &
                                              STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutputSwanASCII - ModuleHDF5ToASCIIandBIN - ERR30'

            else if (Me%OutPutXY == CartesianXY_) then 

                call GetHorizontalGrid(Me%ObjHorizontalGrid, XX2D_Z = Me%X2D,           &
                                                             YY2D_Z = Me%Y2D,           &
                                                             STAT   = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutputSwanASCII - ModuleHDF5ToASCIIandBIN - ERR40'


            endif

        endif
            

    end subroutine ConstructGrid

    
    !------------------------------------------------------------------------

    

    
    !------------------------------------------------------------------------
    subroutine OutputFields(Aux2D, Aux2DNext, Aux2DPrev, NextTime, PrevTime, l, p)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer                   :: Aux2D, Aux2DNext, Aux2DPrev
        type(T_Time)                                    :: NextTime, PrevTime 
        integer                                         :: l, p       
        !Local-----------------------------------------------------------------
        character(len=14)                               :: TimeName
        real, dimension(:  ), pointer                   :: AuxTime
        real                                            :: dt1, dt2
        integer                                         :: STAT_CALL, i, j, n
        !Begin-----------------------------------------------------------------

        allocate(AuxTime(1:6))

   
        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,            &
                           Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleHDF5ToASCIIandBIN - ERR60'


        if (l == 1) then

            call HDF5ReadData(Me%ObjHDF5,                                              &
                               "/Results/"//trim(Me%PropsName(p)),                      &
                               trim(Me%PropsName(p)),                                   &
                               Array2D      = Aux2DPrev,                                &
                               OutputNumber = l,                                        &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleHDF5ToASCIIandBIN - ERR70'

            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                if (Aux2DPrev(i,j) == FillValueReal) then
                    Aux2DPrev(i,j) = -999.
                endif
            enddo
            enddo

        else

            Aux2DPrev(:,:) = Aux2DNext(:,:)

        endif

        call HDF5ReadData(Me%ObjHDF5,                                                   &
                           "/Results/"//trim(Me%PropsName(p)),                          &
                           trim(Me%PropsName(p)),                                       &
                           Array2D      = Aux2DNext,                                    &
                           OutputNumber = l+1,                                          &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleHDF5ToASCIIandBIN - ERR70'

        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            if (Aux2DNext(i,j) == FillValueReal) then
                Aux2DNext(i,j) = -999.
            endif
        enddo
        enddo

 
        n = Me%Output%Next

        if (n <= Me%Output%Number) then

    dw1:    do while (NextTime >= Me%Output%OutTime(n) .and. PrevTime <= Me%Output%OutTime(n)) 

                call ExtractDate(Me%Output%OutTime(n), AuxTime(1), AuxTime(2), AuxTime(3), AuxTime(4), AuxTime(5), AuxTime(6)) 

                write(TimeName(1:14),"(I4,5I2)") int(AuxTime(1)), int(AuxTime(2)), int(AuxTime(3)), &
                                                 int(AuxTime(4)), int(AuxTime(5)), int(AuxTime(6))

                do j = 1, 14
                    if (TimeName(j:j) ==' ') TimeName(j:j) = "0"
                enddo
         
                dt1 = NextTime              - Me%Output%OutTime(n) 
                dt2 = Me%Output%OutTime(n)  - PrevTime     

                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    Aux2D(i,j) =  (Aux2DNext(i,j) * dt2 + Aux2DPrev(i,j) * dt1) / (dt1 + dt2)
                enddo
                enddo

                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%ConvertProp(p) == NauticalToCartDegrees_) then
                        Aux2D(i,j) = 270. - Aux2D(i,j)
                    endif
                enddo
                enddo

                if      (Me%OutPutOption < Mohid_) then
                    call OutputSwanASCII (trim(Me%PropsName(p)), TimeName, Aux2D) 
                else if (Me%OutPutOption == Mohid_) then
                    call OutputMohidBin  (trim(Me%PropsName(p)), TimeName, Aux2D, l) 
                endif

                if (p==1) Me%Output%Next = Me%Output%Next + 1 

                n = Me%Output%Next

                if (n > Me%Output%Number) exit

            enddo dw1

        endif


        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleHDF5ToASCIIandBIN - ERR90'


        deallocate(AuxTime)

    end subroutine OutputFields


    !----------------------------------------------------------------------

   !------------------------------------------------------------------------
    subroutine Open_HDF5_InPut_File

        !Local-----------------------------------------------------------------
        !real,    dimension(:,:), pointer            :: Bathymetry
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_READ

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

!        call GetGridData        (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_InPut_File - ModuleHDF5ToASCIIandBIN - ERR10'

      
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%InPutFileName, HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_InPut_File - ModuleHDF5ToASCIIandBIN - ERR30'
        
        
!        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
!                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_InPut_File - ModuleHDF5ToASCIIandBIN - ERR40'

            
!        call HDF5ReadData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",                   &
!                              Array2D =  Bathymetry, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_InPut_File - ModuleHDF5ToASCIIandBIN - ERR50'            


!        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_InPut_File - ModuleHDF5ToASCIIandBIN - ERR60'            
   
       
!        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
!                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_InPut_File - ModuleHDF5ToASCIIandBIN - ERR70'            

!        call HDF5ReadData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",                &
!                              Array2D = Me%WaterPoints2D,  STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_InPut_File - ModuleHDF5ToASCIIandBIN - ERR80'


!        call UnGetGridData      (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_InPut_File - ModuleHDF5ToASCIIandBIN - ERR90'

        
        !Writes everything to disk
!        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_InPut_File - ModuleHDF5ToASCIIandBIN - ERR110'

         call GetHDF5GroupNumberOfItems (Me%ObjHDF5, "/Time", Me%TotalDates, STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_InPut_File - ModuleHDF5ToASCIIandBIN - ERR120'

        if (Me%TotalDates <= 1) then
            write(*,*) 'The number of instants in the hdf5 must be bigger than 1'
            stop 'Open_HDF5_InPut_File - ModuleHDF5ToASCIIandBIN - ERR130'
        endif



    end subroutine Open_HDF5_InPut_File

    !--------------------------------------------------------------------------

   
    !--------------------------------------------------------------------------

    
    subroutine KillHDF5ToASCIIandBIN
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        !call KillHorizontalMap(Me%ObjHorizontalMap, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'KillHDF5ToASCIIandBIN - ModuleHDF5ToASCIIandBIN - ERR20'

        !call KillGridData(Me%ObjBathymetry, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'KillHDF5ToASCIIandBIN - ModuleHDF5ToASCIIandBIN - ERR30'


        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillHDF5ToASCIIandBIN - ModuleHDF5ToASCIIandBIN - ERR40'


        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillHDF5ToASCIIandBIN - ModuleHDF5ToASCIIandBIN - ERR60'


        deallocate(Me%PropsName )
        nullify   (Me%PropsName )

        deallocate(Me%PropVector)
        nullify   (Me%PropVector)

        deallocate(Me%ConvertProp)
        nullify   (Me%ConvertProp)

        if (associated(Me%OutPut%OutTime)) then
            deallocate(Me%OutPut%OutTime)
            nullify   (Me%OutPut%OutTime)
        endif

        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillHDF5ToASCIIandBIN

    !--------------------------------------------------------------------------
 
end module ModuleHDF5ToASCIIandBIN
