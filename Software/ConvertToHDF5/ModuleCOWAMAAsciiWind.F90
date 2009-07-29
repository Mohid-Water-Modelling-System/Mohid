!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : CowamaAsciiWind
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group & MeteoGalicia
! DATE          : September 2003
! REVISION      : Pablo Carracedo - v4.0
! DESCRIPTION   : Module to convert CowamaAsciiWind files into HDF5 format.
!                 For reading into Mohid module HydrodynamicFile
!
!------------------------------------------------------------------------------


Module ModuleCowamaAsciiWind

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
    public  :: ConvertCowamaAsciiWind
    private ::      ReadGlobalOptions
!    private ::      ReadFieldDates
    private ::      ConstructGlobalOutput
    private ::      ConstructBathymetry
    private ::      Open_HDF5_OutPut_File
    private ::      Construct_FieldList
    private ::          ReadNewField
    private ::              AddField
    private ::          ConvertToMohidUnits   
    private ::          OutputFields
    private ::      KillCowamaAsciiWind

    !Parameters----------------------------------------------------------------

    character(LEN = StringLength), parameter    :: field_block_begin   = '<<beginfield>>'
    character(LEN = StringLength), parameter    :: field_block_end     = '<<endfield>>'
    
    integer,                       parameter    :: UpRight_            = 1
    integer,                       parameter    :: DownRight_          = 2

    !Types---------------------------------------------------------------------
    
    type       T_OutPut                                 
         type (T_Time), pointer, dimension(:)               :: OutTime
         integer                                            :: Number
         logical                                            :: ON
    end type T_OutPut                                   


    private :: T_Field
    type       T_Field
        type(T_PropertyID)                      :: ID
        character(len=StringLength)             :: Units
        character(len=PathLength)               :: File
        !integer                                 :: StartChar
        integer                                 :: FileColumn
        real,    dimension(:,:  ),  pointer     :: Aux2DNext, Aux2DPrevious, Aux2D
        type(T_Time)                            :: NextTime, PreviousTime
        logical                                 :: FirstField
        type(T_Field),              pointer     :: Next, Prev
        integer                                 :: NextOutPut
    end type  T_Field

    
    private :: T_CowamaAsciiWind
    type       T_CowamaAsciiWind
        integer                                          :: ObjEnterData         = 0
        integer                                          :: ObjHDF5              = 0
        integer                                          :: ObjHorizontalGrid    = 0
        integer                                          :: ObjBathymetry        = 0
        integer                                          :: ObjHorizontalMap     = 0
        integer                                          :: ObjTime              = 0
        integer                                          :: Unit
        character(len=PathLength)                        :: FileName
        character(len=PathLength)                        :: GridFileName
        character(len=PathLength)                        :: OutputFileName
        character(len=PathLength)                        :: DatesFile
        integer                                          :: NumDates = FillValueInt
        type(T_Time)                                     :: StartTime, EndTime
        real                                             :: ForecastPeriod, ForecastDT
        type(T_Time),     dimension(:), pointer          :: ExistingDates
        character(len=PathLength), dimension(:), pointer :: FileNames
        !type(T_Time),     dimension(:), pointer          :: InputDates
        !character(len=8), dimension(:), pointer          :: DatesName
        integer                                          :: FieldsNumber
        integer                                          :: ClientNumber
        logical                                          :: FirstInputFiles = .false., FirstMapping = .true.
        logical                                          :: WriteVelModulus = .false., WriteWindModulus = .false.
        logical                                          :: FirstFieldON    = .true.
        integer                                          :: ReadMatrixOrder = UpRight_
        real                                             :: FillValue
        integer, dimension(:,:  ),  pointer              :: WaterPoints2D
        type(T_OutPut)                                   :: OutPut
        type(T_Size2D)                                   :: WorkSize, Size
        type(T_Field),              pointer              :: FirstField         
        type(T_Field),              pointer              :: LastField              
    end type  T_CowamaAsciiWind

    type(T_CowamaAsciiWind), pointer              :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertCowamaAsciiWind(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        Me%ClientNumber = ClientNumber

        call ReadGlobalOptions

        !call ReadFieldDates

        call ConstructGlobalOutput

        call ConstructBathymetry

  
        call Open_HDF5_OutPut_File

        call Construct_FieldList


        if (Me%WriteWindModulus) then
            call WriteVelocityModulus(WindVelocityX_, WindVelocityY_, WindModulus_)
        endif


        call KillCowamaAsciiWind


        STAT = SUCCESS_


    end subroutine ConvertCowamaAsciiWind

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
                     ClientModule = 'CowamaAsciiWind',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleCowamaAsciiWind - ERR10'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleCowamaAsciiWind - ERR20'

        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'CowamaAsciiWind',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleCowamaAsciiWind - ERR30'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleCowamaAsciiWind - ERR40'


        !call GetData(Me%DatesFile,                                                      &
        !             Me%ObjEnterData, iflag,                                            &
        !             SearchType   = FromBlock,                                          &
        !             keyword      = 'FILENAME_DATES',                                   &
        !             ClientModule = 'CowamaAsciiWind',                                   &
        !             STAT         = STAT_CALL)        
        !if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleCowamaAsciiWind - ERR50'
        !if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleCowamaAsciiWind - ERR60'
        
       call GetData(Me%FillValue,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FILL_VALUE',                                       &
                     default      = -99.999900,                                         &
                     ClientModule = 'CowamaAsciiWind',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleCowamaAsciiWind - ERR70'


       call GetData(Me%StartTime,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'START_TIME',                                       &
                     ClientModule = 'CowamaAsciiWind',                                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleCowamaAsciiWind - ERR80'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleCowamaAsciiWind - ERR90'


       call GetData(Me%ForecastPeriod,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FORECAST_PERIOD',                                  &
                     default      = 172800.,                                            &
                     ClientModule = 'CowamaAsciiWind',                                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleCowamaAsciiWind - ERR100'

        Me%EndTime = Me%StartTime + Me%ForecastPeriod

       call GetData(Me%ForecastDT,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FORECAST_DT',                                      &
                     default      = 10800.,                                             &
                     ClientModule = 'CowamaAsciiWind',                                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleCowamaAsciiWind - ERR110'

        
        Me%NumDates = int(Me%ForecastPeriod / Me%ForecastDT) + 1


       call GetData(Me%ReadMatrixOrder,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'READ_MATRIX_ORDER',                                &
                     default      = UpRight_,                                           &
                     ClientModule = 'CowamaAsciiWind',                                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleCowamaAsciiWind - ERR120'
        
        if (Me%ReadMatrixOrder /= UpRight_ .and. Me%ReadMatrixOrder /= DownRight_) then
            write (*,*) 'The option ', Me%ReadMatrixOrder, ' do not exist'
            stop 'ReadGlobalOptions - ModuleCowamaAsciiWind - ERR130'
        endif        
        
    end subroutine ReadGlobalOptions

    !--------------------------------------------------------------------------

!    subroutine ReadFieldDates

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
!        type(T_Time)                                :: AuxTime
!        real                                        :: Year, Month, Day, Hour, Minute, Second  
!        integer                                     :: Aux, l, STAT_CALL
!        character(Len=8)                            :: Char8

        !Begin-----------------------------------------------------------------

!        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
!        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldDates - ModuleCowamaAsciiWind - ERR10'

!        open(Unit   = Me%Unit,                                                          &
!             File   = Me%DatesFile,                                                     &
!             Form   = 'FORMATTED',                                                      &
!             STATUS = 'OLD',                                                            &
!             Action = 'READ',                                                           &
!             IOSTAT = STAT_CALL) 
!        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldDates - ModuleCowamaAsciiWind - ERR20'
    
        !read (Me%Unit, *) Me%NumDates

        !allocate(Me%InputDates(Me%NumDates))  

        !allocate(Me%DatesName  (Me%NumDates))  

        !do l=1, Me%NumDates

!            read(Me%Unit,'(8A)') Char8

!            Me%DatesName(l) =  Char8

!            if (Char8(1:1) =='0') then
!                read (Char8(2:2), '(I1)') Aux
!            else
!                read (Char8(1:2), '(I2)') Aux
!            endif

!            Hour = real(Aux)

!            if (Char8(3:3) =='0') then
!                read (Char8(4:4), '(I1)') Aux
!            else
!                read (Char8(3:4), '(I2)') Aux
!            endif

!            Year = real(Aux) + 2000.

!            if (Char8(5:5) =='0') then
!                read (Char8(6:6), '(I1)') Aux
!            else
!                read (Char8(5:6), '(I2)') Aux
!            endif

!            Month = real(Aux)

!            if (Char8(7:7) =='0') then
!                read (Char8(8:8), '(I1)') Aux
!            else
!                read (Char8(7:8), '(I2)') Aux
!            endif

!            Day   = real(Aux)
  
!            Minute = 0.
!            Second = 0.

!            call SetDate(AuxTime         ,                                              &
!                        Year    = Year   ,                                              &
!                        Month   = Month  ,                                              & 
!                        Day     = Day    ,                                              &
!                        Hour    = Hour   ,                                              &
!                        Minute  = Minute ,                                              &
!                        Second  = Second )  

!            Me%InputDates(l) = AuxTime

!        enddo

!        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
!        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleCowamaAsciiWind - ERR30'


 !   end subroutine ReadFieldDates 

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ConstructGlobalOutput 

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        

        nullify(Me%OutPut%OutTime)


        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime      = Me%StartTime,                             &
                           EndTime          = Me%EndTime,                               &
                           keyword          = 'OUTPUT_TIME',                            &
                           SearchType       = FromFile,                                 &
                           OutPutsTime      = Me%OutPut%OutTime,                        &
                           OutPutsOn        = Me%OutPut%ON,                             &
                           OutPutsNumber    = Me%OutPut%Number,                         &
                           STAT             = STAT_CALL)                                
                                                                                        
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructGlobalOutput - ModuleCowamaAsciiWind - ERR10'              

!        if (.not. Me%OutPut%ON) then
!            Me%OutPut%OutTime => Me%InputDates
!            Me%OutPut%Number  =  Me%NumDates
!        endif

                                                                                

    end subroutine ConstructGlobalOutput

    !--------------------------------------------------------------------------

    
    subroutine Construct_FieldList

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Field), pointer         :: NewField
        integer                         :: STAT_CALL
        logical                         :: BlockFound, FirstField



        !------------------------------------------------------------------------
        FirstField = .true.

do1 :   do
            call ExtractBlockFromBlock (Me%ObjEnterData,                                &
                                        ClientNumber      = Me%ClientNumber,            &
                                        block_begin       = field_block_begin,          &
                                        block_end         = field_block_end,            &
                                        BlockInBlockFound = BlockFound,                 &
                                        STAT              = STAT_CALL)
cd1 :       if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :           if (BlockFound) then                                                  
                    ! Construct a New Field 
                    Call Construct_Field(NewField, FirstField)

                    ! Add new Field to the WaterProperties List 
                    Call Add_Field(NewField)

                    FirstField = .false.


                else cd2
!                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

!                    if (STAT_CALL .NE. SUCCESS_)                                &
!                        stop 'Construct_FieldList - ModuleCowamaAsciiWind - ERR01'

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBlock. '
                stop 'Construct_FieldList - ModuleCowamaAsciiWind - ERR02'
            else cd1
                stop 'Construct_FieldList - ModuleCowamaAsciiWind - ERR03'
            end if cd1
        end do do1

        !------------------------------------------------------------------------

    end subroutine Construct_FieldList

    !----------------------------------------------------------------------------

    subroutine Construct_Field(NewField, FirstField)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer       :: NewField
        logical                      :: FirstField

        !----------------------------------------------------------------------
             
        allocate (NewField)

        NewField%FirstField = FirstField


        !Construct Field ID
        call ConstructPropertyID     (NewField%ID, Me%ObjEnterData, FromBlockinBlock)


        !Construct Field values
        call Construct_InputFiles    (NewField)

        Me%FirstFieldON =.false.
       
    end subroutine Construct_Field



    !----------------------------------------------------------------------------

    subroutine Construct_InputFiles(NewField)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer       :: NewField

        !Local-----------------------------------------------------------------
        character(len=PathLength), dimension(:), pointer :: AuxFileNames
        type(T_Time),              dimension(:), pointer :: AuxExistingDates, AuxOutPutTime
        character(len=PathLength)                        :: FileName
        character(len=4)                                 :: AuxC
        integer                                          :: l, STAT_CALL, iflag, AuxI, i, n, istart, iend
        real                                             :: Year, Month, Day, Hour, Minute, Second
        logical                                          :: exist

        !----------------------------------------------------------------------

        call GetData(NewField%File,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlockInBlock,                                   &
                     keyword      = 'FILE',                                             &
                     ClientModule = 'CowamaAsciiWind',                                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleCowamaAsciiWind - ERR10'
        if (iflag     == 0)        stop 'Construct_InputFiles - ModuleCowamaAsciiWind - ERR20'


        call ExtractDate(Me%StartTime    ,                                              &
                    Year    = Year   ,                                                  &
                    Month   = Month  ,                                                  & 
                    Day     = Day    ,                                                  &
                    Hour    = Hour   ,                                                  &
                    Minute  = Minute ,                                                  &
                    Second  = Second )  

        write(AuxC,'(I4)') int(Year)
        
        NewField%File = trim(NewField%File)//'_'//AuxC

        AuxC='0000'

        if(Month<10) then
            write (AuxC(2:2),'(I1)') int(Month)
        else
            write (AuxC(1:2),'(I2)') int(Month)
        endif

        NewField%File = trim(NewField%File)//AuxC(1:2)

        AuxC='0000'

        if(Day<10) then
            write (AuxC(2:2),'(I1)') int(Day)
        else
            write (AuxC(1:2),'(I2)') int(Day)
        endif

        NewField%File = trim(NewField%File)//AuxC(1:2)


        AuxC='0000'

        if(Hour<10) then
            write (AuxC(2:2),'(I1)') int(Hour)
        else
            write (AuxC(1:2),'(I2)') int(Hour)
        endif

        NewField%File = trim(NewField%File)//'_'//AuxC(1:2)

        call GetData(NewField%FileColumn,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlockInBlock,                                   &
                     keyword      = 'FILE_COLUMN',                                      &
                     ClientModule = 'CowamaAsciiWind',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleCowamaAsciiWind - ERR30'
        if (iflag     == 0)        stop 'Construct_InputFiles - ModuleCowamaAsciiWind - ERR40'


!         call GetData(NewField%StartChar,                                                &
!                     Me%ObjEnterData, iflag,                                            &
!                     SearchType   = FromBlockInBlock,                                   &
!                     keyword      = 'START_CHAR',                                       &
!                     default      = 1,                                                  &
!                     ClientModule = 'CowamaAsciiWind',                                   &
!                     STAT         = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleCowamaAsciiWind - ERR60'

        allocate(NewField%Aux2DNext    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(NewField%Aux2DPrevious(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(NewField%Aux2D        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

        call null_time(NewField%PreviousTime) 
        call null_time(NewField%NextTime    )

        NewField%NextTime = Me%StartTime

i2:     if (Me%FirstFieldON) then

            allocate(AuxFileNames(1:Me%NumDates), AuxExistingDates(1:Me%NumDates))

            n = 0

            do l = 1, Me%NumDates

                AuxI = (l-1) * int(Me%ForecastDT/3600.)
            
                AuxC='0000'
            
                if (AuxI<10) then
                    write(AuxC(2:2),'(I1)') AuxI
                else
                    write(AuxC(1:2),'(I2)') AuxI
                endif

                FileName = trim(NewField%File)//'_'//AuxC(1:2)

                inquire(FILE = FileName, EXIST = exist)

                if (exist) then

                    n = n + 1
            
                    AuxExistingDates(n) = Me%StartTime + Me%ForecastDT * real(l-1)

                    AuxFileNames    (n) = FileName

                else
                    write(*,*) 'It should exist the follow file ', trim(FileName)
                endif

            enddo

            Me%NumDates = n

            allocate(Me%FileNames(1:Me%NumDates), Me%ExistingDates(1:Me%NumDates))

            Me%FileNames    (1:Me%NumDates) = AuxFileNames(1:Me%NumDates)
            Me%ExistingDates(1:Me%NumDates) = AuxExistingDates(1:Me%NumDates)

            deallocate(AuxFileNames, AuxExistingDates)

            allocate(AuxOutPutTime(1:Me%Output%Number))

            AuxOutPutTime(1:Me%Output%Number) = Me%Output%OutTime(1:Me%Output%Number)

            l = 0

            do i=1, Me%Output%Number

                l = l + 1
        
                if (AuxOutPutTime(i) >= Me%ExistingDates(1)) exit 

            enddo

            istart = l

            l = 0
        
            do i=1, Me%Output%Number

                l = l + 1
        
                if (AuxOutPutTime(i) > Me%ExistingDates(Me%NumDates)) then
                    l = l - 1
                    exit 
                endif
            
                if (AuxOutPutTime(i) == Me%ExistingDates(Me%NumDates)) exit
            enddo

            iend = l

            Me%Output%Number = iend - istart + 1

            deallocate(Me%Output%OutTime)

            allocate(Me%Output%OutTime(1:Me%Output%Number))
   
            Me%Output%OutTime(1:Me%Output%Number) = AuxOutPutTime(istart:iend)

            deallocate(AuxOutPutTime)

        endif i2

        call ActualizeNextField     (NewField, Me%FileNames(1))

        NewField%NextOutPut = 1

        do n = 2, Me%NumDates

            call ActualizeNextField     (NewField, Me%FileNames(n))

            NewField%PreviousTime = Me%ExistingDates(n-1)

            NewField%NextTime     = Me%ExistingDates(n  )

            call ConvertToMohidUnits    (NewField)

            call OutputFields           (NewField)

        enddo



    end subroutine Construct_InputFiles


    !--------------------------------------------------------------------------

    subroutine ActualizeNextField (NewField, FileName)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                      :: NewField
        character(len=*)                            :: FileName

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeNextField - ModuleCowamaAsciiWind - ERR10'

        open(Unit   = Me%Unit,                                                      &
             File   = FileName,                                                     &
             Form   = 'FORMATTED',                                                  &
             STATUS = 'OLD',                                                        &
             Action = 'READ',                                                       &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeNextField - ModuleCowamaAsciiWind - ERR20'        

        call ReadNewField           (NewField)

        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeNextField - ModuleCowamaAsciiWind - ERR30'

    end subroutine ActualizeNextField

    !--------------------------------------------------------------------------

    subroutine ReadNewField (NewField)

        !Arguemtns-------------------------------------------------------------
        type(T_Field), pointer                      :: NewField

        !Local----------------------------------------------------------------
        real,  dimension(:), allocatable            :: Aux
        character(len=StringLength)                 :: AuxString
        integer                                     :: i, j, l
        !Begin-----------------------------------------------------------------
        
        l = NewField%FileColumn 

        NewField%Aux2DPrevious(:,:) = NewField%Aux2DNext(:,:)
        
       
        if (Me%ReadMatrixOrder == UpRight_) then

            do i=Me%WorkSize%ILB,Me%WorkSize%IUB 
            do j=Me%WorkSize%JLB,Me%WorkSize%JUB
        
                NewField%Aux2DNext(i,j) = ReadValue(NewField%ID%Name,l)

            enddo
            enddo

        else if (Me%ReadMatrixOrder == DownRight_) then

            do i=Me%WorkSize%IUB,Me%WorkSize%ILB, - 1
            do j=Me%WorkSize%JLB,Me%WorkSize%JUB
        
                NewField%Aux2DNext(i,j) = ReadValue(NewField%ID%Name, l)

            enddo
            enddo
        
        
        endif

    end subroutine ReadNewField
    

    !----------------------------------------------------------------------------

    real function ReadValue (FieldName, l)

        !Arguments---------------------------------------------------------------
        character (len=*)                   :: FieldName
        integer                             :: l

        !Local-------------------------------------------------------------------
        character (len=1000)                :: AuxString
        integer                             :: k, v, iostat

        !Begin-------------------------------------------------------------------

        read(Me%Unit,'(1000A)') AuxString

        AuxString = adjustl(AuxString)

        do k=1,l-1
            do v =1,1000
                if (AuxString(v:v)==' ') exit
            enddo
            AuxString = AuxString(v:1000)
            AuxString = adjustl(AuxString)
       
       enddo


        do v =1,1000
            if (AuxString(v:v)==' ') exit
        enddo

        read(AuxString(1:v-1),*, iostat = iostat) ReadValue

        if (iostat /= 0) then
            ReadValue = 0.
        else
            if (FieldName == "mean wave direction") then
                ReadValue = ReadValue - 180.
                if (ReadValue < 0.) ReadValue = ReadValue + 360. 
            endif
        endif

    end function ReadValue
    !------------------------------------------------------------------------

    
    subroutine ConstructBathymetry
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !Me%GridFileName="NewGrid.dat_.new"
        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleCowamaAsciiWind - ERR10'

        call GetHorizontalGridSize(Me%ObjHorizontalGrid, Me%Size, Me%WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleCowamaAsciiWind - ERR20'

        call ConstructGridData(Me%ObjBathymetry, Me%ObjHorizontalGrid, FileName = Me%GridFileName,&
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleCowamaAsciiWind - ERR30'


        call ConstructHorizontalMap(Me%ObjHorizontalMap, Me%ObjBathymetry, Me%ObjHorizontalGrid, &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleCowamaAsciiWind - ERR40'

        call GetWaterPoints2D   (Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructBathymetry - ModuleCowamaAsciiWind - ERR50'


    end subroutine ConstructBathymetry

    
    !------------------------------------------------------------------------

    
    !------------------------------------------------------------------------

    subroutine ConvertToMohidUnits(NewField)

        !Local-----------------------------------------------------------------
        type (T_Field), pointer                   :: NewField
        
        !Begin-----------------------------------------------------------------

        select case(NewField%ID%IDNumber)

            case(WindVelocityX_)

                NewField%ID%Units     = 'm/s'

            case(WindVelocityY_)

                NewField%ID%Units     = 'm/s'

            case default

        end select



    end subroutine ConvertToMohidUnits

    !------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    ! This subroutine adds a new property to the Water Property List  

    subroutine Add_Field(NewField)

        !Arguments-------------------------------------------------------------
        type(T_Field),           pointer     :: NewField

        !----------------------------------------------------------------------

        ! Add to the WaterField List a new property
        if (.not.associated(Me%FirstField)) then
            Me%FieldsNumber = 1
            Me%FirstField    => NewField
            Me%LastField     => NewField
        else
            NewField%Prev     => Me%LastField
            Me%LastField%Next => NewField
            Me%LastField      => NewField
            Me%FieldsNumber   = Me%FieldsNumber + 1
        end if 

        !----------------------------------------------------------------------

    end subroutine Add_Field 

    !--------------------------------------------------------------------------

    
    subroutine AddField (FirstField, ObjField)

        !Arguments-------------------------------------------------------------
        type (T_Field), pointer                   :: FirstField
        type (T_Field), pointer                   :: ObjField

        !Local-----------------------------------------------------------------
        type (T_Field), pointer                   :: NewField
        type (T_Field), pointer                   :: PreviousField
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewField)
        nullify  (NewField%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstField)) then
            FirstField         => NewField
            ObjField           => NewField
        else
            PreviousField      => FirstField
            ObjField           => FirstField%Next
            do while (associated(ObjField))
                PreviousField  => ObjField
                ObjField       => ObjField%Next
            enddo
            ObjField           => NewField
            PreviousField%Next => NewField
        endif


    end subroutine AddField
    
    
    !------------------------------------------------------------------------

    
    
    !------------------------------------------------------------------------
    subroutine OutputFields(Field)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                          :: Field
        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        integer                                         :: STAT_CALL, i
        !Begin-----------------------------------------------------------------



        i = Field%NextOutPut

i0:     if (i<=Me%Output%Number) then

dw1:    do while (Field%NextTime >= Me%Output%OutTime(i)) 
       
            call InterpolateField(Field, i)        
            
i1:         if (Field%FirstField) then

        !           Dados para escriver uma soa vez cada date:
                call ExtractDate   (Me%Output%OutTime(i),                               &
                                    AuxTime(1), AuxTime(2), AuxTime(3),                 &
                                    AuxTime(4), AuxTime(5), AuxTime(6))

                TimePtr => AuxTime

                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCowamaAsciiWind - ERR10'


                call HDF5WriteData  (Me%ObjHDF5, "/Time",                               &
                                     "Time", "YYYY/MM/DD HH:MM:SS",                     &
                                     Array1D = TimePtr,                                 &
                                     OutputNumber = i, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCowamaAsciiWind - ERR20'

            endif i1

   
            call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,            &
                               Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCowamaAsciiWind - ERR60'


            call HDF5WriteData(Me%ObjHDF5,                                              &
                               "/Results/"//Field%ID%Name,                              &
                               Field%ID%Name,                                           &
                               Field%ID%Units,                                          &
                               Array2D      = Field%Aux2D,                              &
                               OutputNumber = i,                                        &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCowamaAsciiWind - ERR70'

         
            if (GetPropertyIDNumber(Field%ID%Name) == WindVelocityX_) Me%WriteWindModulus = .true.

            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCowamaAsciiWind - ERR80'

            if (i==Me%OutPut%Number) exit

            i = i + 1

        enddo dw1

        endif i0

        Field%NextOutPut = i

    end subroutine OutputFields


    !----------------------------------------------------------------------

    !------------------------------------------------------------------------
    subroutine InterpolateField(Field, i) 

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                          :: Field
        integer                                         :: i
        !Local-----------------------------------------------------------------
        real                                            :: dt1, dt2
        integer                                         :: ii, jj
        !Begin-----------------------------------------------------------------




        dt1 = Field%NextTime        - Me%Output%OutTime(i) 
        dt2 = Me%Output%OutTime(i)  - Field%PreviousTime
       
        Field%Aux2D(:,:) = 0.

        do jj= Me%WorkSize%JLB,Me%WorkSize%JUB  
        do ii= Me%WorkSize%ILB,Me%WorkSize%IUB  

            if (Me%WaterPoints2D(ii,jj) == WaterPoint) then

                Field%Aux2D(ii,jj) = (Field%Aux2DPrevious(ii, jj) * dt1   +             &
                                      Field%Aux2DNext    (ii, jj) * dt2 ) / (dt1 + dt2)


            endif

        enddo
        enddo

    !------------------------------------------------------------------------
    end subroutine InterpolateField

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
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleCowamaAsciiWind - ERR10'

      
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleCowamaAsciiWind - ERR30'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleCowamaAsciiWind - ERR40'

            
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",                   &
                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleCowamaAsciiWind - ERR50'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleCowamaAsciiWind - ERR60'            
   
       
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleCowamaAsciiWind - ERR70'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",                &
                              Array2D = Me%WaterPoints2D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleCowamaAsciiWind - ERR80'


        call UnGetGridData      (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleCowamaAsciiWind - ERR90'

        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleCowamaAsciiWind - ERR110'

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

    subroutine WriteVelocityModulus(VelocityU, VelocityV, VelocityModulus)

        !Arguments-------------------------------------------------------------
        integer                                     :: VelocityU, VelocityV, VelocityModulus

        !Local-----------------------------------------------------------------
        real   , dimension(:,:), pointer            :: VelX, VelY, VelM
        character(len = StringLength)               :: GroupNameX, GroupNameY, GroupNameM
        character(len = StringLength)               :: NameX, NameY, NameM
        integer                                     :: STAT_CALL, n, nItems, i, j
        integer                                     :: HDF5_READWRITE

        !----------------------------------------------------------------------
        
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteVelocityModulus - ModuleCowamaAsciiWind - ERR10'

        allocate(VelX(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(VelY(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(VelM(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

       
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READWRITE = HDF5_READWRITE)

        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_READWRITE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteVelocityModulus - ModuleCowamaAsciiWind - ERR20'

        NameX      = GetPropertyName(VelocityU)
        GroupNameX = '/Results/'//trim(NameX)

        NameY      = GetPropertyName(VelocityV)
        GroupNameY = '/Results/'//trim(NameY)

        NameM      = GetPropertyName(VelocityModulus)
        GroupNameM = '/Results/'//trim(NameM)
        
        call GetHDF5GroupNumberOfItems (Me%ObjHDF5, GroupName = trim(GroupNameX), nItems = nItems, STAT= STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteVelocityModulus - ModuleCowamaAsciiWind - ERR30'

        do n=1, nItems 

            call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,            &
                               Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCowamaAsciiWind - ERR40'

            call HDF5ReadData(Me%ObjHDF5, GroupName = trim(GroupNameX),                 &
                              Name = trim(NameX), Array2D = VelX, OutputNumber = n,     &
                              STAT= STAT_CALL)
        

            call HDF5ReadData(Me%ObjHDF5, GroupName = trim(GroupNameY),                 &
                              Name = trim(NameY), Array2D = VelY, OutputNumber = n,     &
                              STAT= STAT_CALL)

            VelM(:, :) = 0.

            do j= Me%WorkSize%JLB,Me%WorkSize%JUB  
            do i= Me%WorkSize%ILB,Me%WorkSize%IUB  

                if (Me%WaterPoints2D(i, j) == WaterPoint) then

                    VelM(i, j) = sqrt(VelX(i,j)*VelX(i,j) + VelY(i,j)*VelY(i,j))

                endif

            enddo
            enddo

            call HDF5WriteData(Me%ObjHDF5,                                              &
                               trim(GroupNameM),                                        &
                               trim(NameM),                                             &
                               'm/s',                                                   &
                               Array2D      = VelM,                                     &
                               OutputNumber = n,                                        &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCowamaAsciiWind - ERR50'


        enddo

        deallocate(VelX, VelY, VelM)

        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteVelocityModulus - ModuleCowamaAsciiWind - ERR60'


    end subroutine WriteVelocityModulus
   
    !--------------------------------------------------------------------------

    
    subroutine KillCowamaAsciiWind
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer      :: Field, FieldNext
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call UnGetHorizontalMap (Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillCowamaAsciiWind - ModuleCowamaAsciiWind - ERR10'

        call KillHorizontalMap(Me%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillCowamaAsciiWind - ModuleCowamaAsciiWind - ERR20'

        call KillGridData(Me%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillCowamaAsciiWind - ModuleCowamaAsciiWind - ERR30'


        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillCowamaAsciiWind - ModuleCowamaAsciiWind - ERR40'


        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillCowamaAsciiWind - ModuleCowamaAsciiWind - ERR60'

!        deallocate(Me%DatesName  )  
!        nullify   (Me%DatesName  )
!        deallocate(Me%InputDates )
!        nullify   (Me%InputDates )
        if (Me%OutPut%ON) then
            deallocate(Me%OutPut%OutTime)
        endif
        nullify   (Me%OutPut%OutTime)


        Field => Me%FirstField 

        do while (associated(Field)) 

            FieldNext => Field%Next

            deallocate(Field%Aux2DNext    )
            deallocate(Field%Aux2DPrevious)
            deallocate(Field%Aux2D        )

            deallocate(Field)
            nullify   (Field)

            Field     => FieldNext
        enddo

        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillCowamaAsciiWind

    !--------------------------------------------------------------------------
 
end module ModuleCowamaAsciiWind
