!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : TecnoceanAscii
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group & MeteoGalicia
! DATE          : September 2003
! REVISION      : Pablo Carracedo - v4.0
! DESCRIPTION   : Module to convert TecnoceanAscii files into HDF5 format.
!                 For reading into Mohid module HydrodynamicFile
!
!------------------------------------------------------------------------------


Module ModuleTecnoceanAscii

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
    public  :: ConvertTecnoceanAscii
    private ::      ReadGlobalOptions
    private ::      ReadFieldDates
    private ::      ConstructGlobalOutput
    private ::      ConstructBathymetry
    private ::      Open_HDF5_OutPut_File
    private ::      Construct_FieldList
    private ::          ReadNewField
    private ::              AddField
    private ::          ConvertToMohidUnits   
    private ::          OutputFields
    private ::      KillTecnoceanAscii

    !Parameters----------------------------------------------------------------

    character(LEN = StringLength), parameter    :: field_block_begin   = '<<beginfield>>'
    character(LEN = StringLength), parameter    :: field_block_end     = '<<endfield>>'

    integer,                       parameter    :: SmallDomainEnric_   = 1
    integer,                       parameter    :: LargeDomainSonia_   = 2


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
        character(len=StringLength)             :: File, FilePath
        integer                                 :: StartChar
        integer                                 :: FileColumn
        real,    dimension(:,:  ),  pointer     :: Aux2DNext, Aux2DPrevious, Aux2D
        type(T_Time)                            :: NextTime, PreviousTime
        logical                                 :: FirstField
        type(T_Field),              pointer     :: Next, Prev
        integer                                 :: NextOutPut
    end type  T_Field

    
    private :: T_TecnoceanAscii
    type       T_TecnoceanAscii
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjBathymetry        = 0
        integer                                 :: ObjHorizontalMap     = 0
        integer                                 :: ObjTime              = 0
        integer                                 :: Unit
        character(len=PathLength)               :: FileName
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: OutputFileName
        character(len=PathLength)               :: DatesFile
        integer                                 :: NumDates = FillValueInt
        type(T_Time),     dimension(:), pointer :: InputDates
        character(len=8), dimension(:), pointer :: DatesName
        integer                                 :: FieldsNumber
        integer                                 :: ClientNumber
        logical                                 :: FirstInputFiles = .false., FirstMapping = .true.
        logical                                 :: WriteVelModulus = .false., WriteWindModulus = .false.
        real                                    :: FillValue
        integer                                 :: ReadType
        integer, dimension(:,:  ),  pointer     :: WaterPoints2D
        type(T_OutPut)                          :: OutPut
        type(T_Size2D)                          :: WorkSize, Size
        type(T_Field),              pointer     :: FirstField         
        type(T_Field),              pointer     :: LastField              
    end type  T_TecnoceanAscii

    type(T_TecnoceanAscii), pointer              :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertTecnoceanAscii(EnterDataID, ClientNumber, STAT)

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

        call ReadFieldDates

        call ConstructGlobalOutput

        call ConstructBathymetry

  
        call Open_HDF5_OutPut_File

        call Construct_FieldList

        if (Me%WriteVelModulus) then
            call WriteVelocityModulus(VelocityU_, VelocityV_, VelocityModulus_)
        endif

        if (Me%WriteWindModulus) then
            call WriteVelocityModulus(WindVelocityX_, WindVelocityY_, WindModulos_)
        endif


        call KillTecnoceanAscii


        STAT = SUCCESS_


    end subroutine ConvertTecnoceanAscii

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
                     ClientModule = 'TecnoceanAscii',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleTecnoceanAscii - ERR10'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleTecnoceanAscii - ERR20'

        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'TecnoceanAscii',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleTecnoceanAscii - ERR30'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleTecnoceanAscii - ERR40'


        call GetData(Me%DatesFile,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FILENAME_DATES',                                   &
                     ClientModule = 'TecnoceanAscii',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleTecnoceanAscii - ERR50'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleTecnoceanAscii - ERR60'
        
       call GetData(Me%FillValue,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FILL_VALUE',                                       &
                     default      = -99.999900,                                         &
                     ClientModule = 'TecnoceanAscii',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleTecnoceanAscii - ERR70'

       call GetData(Me%ReadType,                                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'READ_TYPE',                                        &
                     default      = SmallDomainEnric_,                                  &
                     ClientModule = 'TecnoceanAscii',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleTecnoceanAscii - ERR80'

    end subroutine ReadGlobalOptions

    !--------------------------------------------------------------------------

    subroutine ReadFieldDates

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_Time)                                :: AuxTime
        real                                        :: Year, Month, Day, Hour, Minute, Second  
        integer                                     :: Aux, l, STAT_CALL
        character(Len=8)                            :: Char8

        !Begin-----------------------------------------------------------------

        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldDates - ModuleTecnoceanAscii - ERR10'

        open(Unit   = Me%Unit,                                                          &
             File   = Me%DatesFile,                                                     &
             Form   = 'FORMATTED',                                                      &
             STATUS = 'OLD',                                                            &
             Action = 'READ',                                                           &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldDates - ModuleTecnoceanAscii - ERR20'
    
        read (Me%Unit, *) Me%NumDates

        allocate(Me%InputDates(Me%NumDates))  

        allocate(Me%DatesName  (Me%NumDates))  

        do l=1, Me%NumDates

            read(Me%Unit,'(8A)') Char8

            Me%DatesName(l) =  Char8

            if (Char8(1:1) =='0') then
                read (Char8(2:2), '(I1)') Aux
            else
                read (Char8(1:2), '(I2)') Aux
            endif

            Hour = real(Aux)

            if (Char8(3:3) =='0') then
                read (Char8(4:4), '(I1)') Aux
            else
                read (Char8(3:4), '(I2)') Aux
            endif

            Year = real(Aux) + 2000.

            if (Char8(5:5) =='0') then
                read (Char8(6:6), '(I1)') Aux
            else
                read (Char8(5:6), '(I2)') Aux
            endif

            Month = real(Aux)

            if (Char8(7:7) =='0') then
                read (Char8(8:8), '(I1)') Aux
            else
                read (Char8(7:8), '(I2)') Aux
            endif

            Day   = real(Aux)
  
            Minute = 0.
            Second = 0.

            call SetDate(AuxTime         ,                                              &
                        Year    = Year   ,                                              &
                        Month   = Month  ,                                              & 
                        Day     = Day    ,                                              &
                        Hour    = Hour   ,                                              &
                        Minute  = Minute ,                                              &
                        Second  = Second )  

            Me%InputDates(l) = AuxTime

        enddo

        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleTecnoceanAscii - ERR30'


    end subroutine ReadFieldDates 

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ConstructGlobalOutput 

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        

        nullify(Me%OutPut%OutTime)


        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime      = Me%InputDates(1),                         &
                           EndTime          = Me%InputDates(Me%NumDates),               &
                           keyword          = 'OUTPUT_TIME',                            &
                           SearchType       = FromFile,                                 &
                           OutPutsTime      = Me%OutPut%OutTime,                        &
                           OutPutsOn        = Me%OutPut%ON,                             &
                           OutPutsNumber    = Me%OutPut%Number,                         &
                           STAT             = STAT_CALL)                                
                                                                                        
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructGlobalOutput - ModuleTecnoceanAscii - ERR10'              

        if (.not. Me%OutPut%ON) then
            Me%OutPut%OutTime => Me%InputDates
            Me%OutPut%Number  =  Me%NumDates
        endif

                                                                                

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
!                        stop 'Construct_FieldList - ModuleTecnoceanAscii - ERR01'

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBlock. '
                stop 'Construct_FieldList - ModuleTecnoceanAscii - ERR02'
            else cd1
                stop 'Construct_FieldList - ModuleTecnoceanAscii - ERR03'
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
       
    end subroutine Construct_Field



    !----------------------------------------------------------------------------

    subroutine Construct_InputFiles(NewField)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer       :: NewField

        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: FileName
        integer                                     :: l, STAT_CALL, iflag

        !----------------------------------------------------------------------

        call GetData(NewField%File,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlockInBlock,                                   &
                     keyword      = 'FILE',                                             &
                     ClientModule = 'TecnoceanAscii',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleTecnoceanAscii - ERR10'
        if (iflag     == 0)        stop 'Construct_InputFiles - ModuleTecnoceanAscii - ERR20'

        call GetData(NewField%FileColumn,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlockInBlock,                                   &
                     keyword      = 'FILE_COLUMN',                                      &
                     ClientModule = 'TecnoceanAscii',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleTecnoceanAscii - ERR30'
        if (iflag     == 0)        stop 'Construct_InputFiles - ModuleTecnoceanAscii - ERR40'


        call GetData(NewField%FilePath,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlockInBlock,                                   &
                     keyword      = 'FILE_PATH',                                        &
                     default      = ' ',                                                &
                     ClientModule = 'TecnoceanAscii',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleTecnoceanAscii - ERR50'


        call GetData(NewField%StartChar,                                                &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlockInBlock,                                   &
                     keyword      = 'START_CHAR',                                       &
                     default      = 1,                                                  &
                     ClientModule = 'TecnoceanAscii',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleTecnoceanAscii - ERR60'

        allocate(NewField%Aux2DNext    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(NewField%Aux2DPrevious(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(NewField%Aux2D        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

        call null_time(NewField%PreviousTime) 
        call null_time(NewField%NextTime    )

        NewField%NextOutPut = 1
    
        do l = 1, Me%NumDates

            call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleTecnoceanAscii - ERR70'

            FileName = trim(NewField%FilePath)//Me%DatesName(l)//trim(NewField%File)

            open(Unit   = Me%Unit,                                                      &
                 File   = FileName,                                                     &
                 Form   = 'FORMATTED',                                                  &
                 STATUS = 'OLD',                                                        &
                 Action = 'READ',                                                       &
                 IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleTecnoceanAscii - ERR80'

            call ReadNewField           (NewField)

            NewField%PreviousTime = NewField%NextTime

            NewField%NextTime     = Me%InputDates(l)

            call ConvertToMohidUnits    (NewField)

            if (l > 1) then
                call OutputFields           (NewField)
            endif

            call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleTecnoceanAscii - ERR90'

        enddo



    end subroutine Construct_InputFiles


    !--------------------------------------------------------------------------

    subroutine ReadNewField (NewField)

        !Arguemtns-------------------------------------------------------------
        type(T_Field), pointer                      :: NewField

        !Local----------------------------------------------------------------
        real,  dimension(:), allocatable            :: Aux
        character(len=StringLength)                 :: AuxString
        integer                                     :: i, j, l, s
        !Begin-----------------------------------------------------------------
        
        l = NewField%FileColumn 

        allocate(Aux(1:l))

        NewField%Aux2DPrevious(:,:) = NewField%Aux2DNext(:,:)

        if      (Me%ReadType == SmallDomainEnric_) then

            do j=Me%WorkSize%JUB, Me%WorkSize%JLB, -1 
            do i=Me%WorkSize%IUB, Me%WorkSize%ILB, -1
        
                read(Me%Unit,'(1000A)') AuxString

                s = NewField%StartChar

                read(AuxString(s:StringLength),*) Aux

                NewField%Aux2DNext(i,j) = Aux(l)

            enddo
            enddo

        else if (Me%ReadType == LargeDomainSonia_) then

            do i=Me%WorkSize%IUB, Me%WorkSize%ILB, -1
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        
                read(Me%Unit,'(1000A)') AuxString

                s = NewField%StartChar

                read(AuxString(s:StringLength),*) Aux

                NewField%Aux2DNext(i,j) = Aux(l)

            enddo
            enddo

        endif

         if (GetPropertyIDNumber(NewField%ID%Name) == VelocityU_  .or.                  &
             GetPropertyIDNumber(NewField%ID%Name) == VelocityV_  )                     &
            NewField%Aux2DNext(:,:) = - NewField%Aux2DNext(:,:)

        deallocate(Aux)

    end subroutine ReadNewField
    

    !----------------------------------------------------------------------------



    
    !------------------------------------------------------------------------

    
    subroutine ConstructBathymetry
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !Me%GridFileName="NewGrid.dat_.new"
        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleTecnoceanAscii - ERR10'

        call GetHorizontalGridSize(Me%ObjHorizontalGrid, Me%Size, Me%WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleTecnoceanAscii - ERR20'

        call ConstructGridData(Me%ObjBathymetry, Me%ObjHorizontalGrid, FileName = Me%GridFileName,&
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleTecnoceanAscii - ERR30'


        call ConstructHorizontalMap(Me%ObjHorizontalMap, Me%ObjBathymetry, Me%ObjHorizontalGrid, &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleTecnoceanAscii - ERR40'

        call GetWaterPoints2D   (Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructBathymetry - ModuleTecnoceanAscii - ERR50'


    end subroutine ConstructBathymetry

    
    !------------------------------------------------------------------------

    
    !------------------------------------------------------------------------

    subroutine ConvertToMohidUnits(NewField)

        !Local-----------------------------------------------------------------
        type (T_Field), pointer                   :: NewField
        
        !Begin-----------------------------------------------------------------

        select case(NewField%ID%IDNumber)

            case(VelocityU_)

                NewField%ID%Units     = 'm/s'

            case(VelocityV_)

                NewField%ID%Units     = 'm/s'

            case(SignificantWaveHeight_)

                NewField%ID%Units     = 'm'

            case(MeanWaveDirection_)

                NewField%ID%Units     = 'º'

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
        real, dimension(:,:), pointer                   :: Aux2D
        real                                            :: Angle
        integer                                         :: STAT_CALL, ii, jj, k, i
        character(len=StringLength)                     :: FieldName
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
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleTecnoceanAscii - ERR10'


                call HDF5WriteData  (Me%ObjHDF5, "/Time",                               &
                                     "Time", "YYYY/MM/DD HH:MM:SS",                     &
                                     Array1D = TimePtr,                                 &
                                     OutputNumber = i, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleTecnoceanAscii - ERR20'

            endif i1

   
            call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,            &
                               Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleTecnoceanAscii - ERR60'


            call HDF5WriteData(Me%ObjHDF5,                                              &
                               "/Results/"//Field%ID%Name,                              &
                               Field%ID%Name,                                           &
                               Field%ID%Units,                                          &
                               Array2D      = Field%Aux2D,                              &
                               OutputNumber = i,                                        &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleTecnoceanAscii - ERR70'

            if (GetPropertyIDNumber(Field%ID%Name) == VelocityU_    ) Me%WriteVelModulus  = .true.
            
            if (GetPropertyIDNumber(Field%ID%Name) == WindVelocityX_) Me%WriteWindModulus = .true.

i2:         if (GetPropertyIDNumber(Field%ID%Name) == MeanWaveDirection_) then

                allocate(Aux2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

d2:             do k=1,2

                    Aux2D(:,:) = 0.

                    do jj= Me%WorkSize%JLB,Me%WorkSize%JUB  
                    do ii= Me%WorkSize%ILB,Me%WorkSize%IUB  

                        if (Me%WaterPoints2D(ii,jj) == WaterPoint) then

                            Angle = - 90 - Field%Aux2D(ii,jj)

                            if (k==1) then
                                Aux2D(ii,jj) = cos(Angle * Pi / 180.) 
                            else
                                Aux2D(ii,jj) = sin(Angle * Pi / 180.) 
                            endif

                        endif

                    enddo
                    enddo
                
                    if (k==1) then
                        FieldName = trim(Field%ID%Name)//'_x'
                    else
                        FieldName = trim(Field%ID%Name)//'_y'
                    endif

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//FieldName,                          &
                                       FieldName,                                       &
                                       '-',                                             &
                                       Array2D      = Aux2D,                            &
                                       OutputNumber = i,                                &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleTecnoceanAscii - ERR80'

                enddo d2

                deallocate(Aux2D)

            endif i2

            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleTecnoceanAscii - ERR80'

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
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleTecnoceanAscii - ERR10'

      
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleTecnoceanAscii - ERR30'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleTecnoceanAscii - ERR40'

            
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",                   &
                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleTecnoceanAscii - ERR50'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleTecnoceanAscii - ERR60'            
   
       
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleTecnoceanAscii - ERR70'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",                &
                              Array2D = Me%WaterPoints2D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleTecnoceanAscii - ERR80'


        call UnGetGridData      (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleTecnoceanAscii - ERR90'

        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleTecnoceanAscii - ERR110'

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
        if (STAT_CALL /= SUCCESS_)stop 'WriteVelocityModulus - ModuleTecnoceanAscii - ERR10'

        allocate(VelX(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(VelY(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(VelM(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

       
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READWRITE = HDF5_READWRITE)

        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_READWRITE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteVelocityModulus - ModuleTecnoceanAscii - ERR20'

        NameX      = GetPropertyName(VelocityU)
        GroupNameX = '/Results/'//trim(NameX)

        NameY      = GetPropertyName(VelocityV)
        GroupNameY = '/Results/'//trim(NameY)

        NameM      = GetPropertyName(VelocityModulus)
        GroupNameM = '/Results/'//trim(NameM)
        
        call GetHDF5GroupNumberOfItems (Me%ObjHDF5, GroupName = trim(GroupNameX), nItems = nItems, STAT= STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteVelocityModulus - ModuleTecnoceanAscii - ERR30'

        do n=1, nItems 

            call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,            &
                               Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleTecnoceanAscii - ERR40'

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
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleTecnoceanAscii - ERR50'


        enddo

        deallocate(VelX, VelY, VelM)

        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteVelocityModulus - ModuleTecnoceanAscii - ERR60'


    end subroutine WriteVelocityModulus
   
    !--------------------------------------------------------------------------

    
    subroutine KillTecnoceanAscii
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer      :: Field, FieldNext
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call UnGetHorizontalMap (Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillTecnoceanAscii - ModuleTecnoceanAscii - ERR10'

        call KillHorizontalMap(Me%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillTecnoceanAscii - ModuleTecnoceanAscii - ERR20'

        call KillGridData(Me%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillTecnoceanAscii - ModuleTecnoceanAscii - ERR30'


        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillTecnoceanAscii - ModuleTecnoceanAscii - ERR40'


        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillTecnoceanAscii - ModuleTecnoceanAscii - ERR60'

        deallocate(Me%DatesName  )  
        nullify   (Me%DatesName  )
        deallocate(Me%InputDates )
        nullify   (Me%InputDates )
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

    
    end subroutine KillTecnoceanAscii

    !--------------------------------------------------------------------------
 
end module ModuleTecnoceanAscii
